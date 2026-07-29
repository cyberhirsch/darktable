/*
    This file is part of darktable,
    Copyright (C) 2026 darktable developers.

    darktable is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    darktable is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with darktable.  If not, see <http://www.gnu.org/licenses/>.
*/

#include "common/lens_share.h"

#include "common/curl_tools.h"
#include "common/darktable.h"
#include "common/file_location.h"
#include "control/conf.h"
#include "control/control.h"
#include "gui/gtk.h"

#include <curl/curl.h>
#include <glib/gstdio.h>
#include <json-glib/json-glib.h>
#include <string.h>

gboolean dt_lens_share_enabled(void)
{
  /* There is no separate master switch. The checkbox beside "save profile"
     in the lensfit panel is the only control -- it persists its own state
     across sessions, and reading that state back is what this checks. A
     second toggle in preferences duplicating the same choice would just be
     two switches wired to the same question. */
  return dt_conf_get_bool("plugins/lens_calib/share_profile");
}

/* Its own file rather than a `dt_conf_set_string` into darktablerc.
 *
 * The config file is only written back to disk on a clean shutdown, so an
 * id that lived only there would be lost to any crash or force-kill --
 * confirmed empirically: two different ids appeared across submissions from
 * this one machine in one session's testing. A fresh id on the next launch
 * makes one person's repeat measurements of the same lens look like
 * independent corroboration from different people, which is exactly what
 * `confidence = agreement * independence * quality` is supposed to guard
 * against and what the project's own README warns about.
 *
 * `g_file_set_contents` writes to a temp file and renames over the target,
 * so a crash mid-write leaves either the old file or the new one, never a
 * half-written id. */
static gchar *_install_id_path(void)
{
  gchar *dir = dt_lens_profile_dir();
  if(!dir) return NULL;
  gchar *path = g_build_filename(dir, ".install_id", NULL);
  g_free(dir);
  return path;
}

const char *dt_lens_share_install_id(void)
{
  static gchar *cached = NULL;
  if(cached) return cached;

  gchar *path = _install_id_path();
  gchar *id = NULL;

  if(path && g_file_test(path, G_FILE_TEST_IS_REGULAR))
  {
    gchar *contents = NULL;
    if(g_file_get_contents(path, &contents, NULL, NULL) && contents)
    {
      g_strstrip(contents);
      if(strlen(contents) >= 8) id = contents;
      else g_free(contents);
    }
  }

  if(!id)
  {
    /* Migrate a pre-existing id out of darktablerc rather than orphaning
       it and minting a second one for the same install. */
    gchar *old = dt_conf_get_string("plugins/lensfit/install_id");
    if(old && strlen(old) >= 8) id = old;
    else g_free(old);
  }

  if(!id)
  {
    /* Random, not derived from anything about this machine. The only thing
       it has to do is stay the same across submissions from one install so
       repeats of the same lens can be recognised as repeats rather than
       counted as independent agreement. */
    id = g_uuid_string_random();
  }

  if(path) g_file_set_contents(path, id, -1, NULL);
  g_free(path);

  cached = id;
  return cached;
}

gchar *dt_lens_share_queue_dir(void)
{
  gchar *base = dt_lens_profile_dir();
  if(!base) return NULL;

  gchar *dir = g_build_filename(base, "queue", NULL);
  g_free(base);

  if(!g_file_test(dir, G_FILE_TEST_IS_DIR))
    g_mkdir_with_parents(dir, 0755);

  return dir;
}

int dt_lens_share_queue_count(void)
{
  gchar *dir = dt_lens_share_queue_dir();
  if(!dir) return 0;

  GDir *d = g_dir_open(dir, 0, NULL);
  g_free(dir);
  if(!d) return 0;

  int n = 0;
  const gchar *name;
  while((name = g_dir_read_name(d)))
    if(g_str_has_suffix(name, ".json")) n++;

  g_dir_close(d);
  return n;
}

/* Context worth having and safe to send: what the lens was mounted on and how
 * big the frame was. Both are properties of equipment, not of a person, and
 * both change how the measurement should be interpreted. The file name, its
 * path, the body serial and any capture time are deliberately absent -- none
 * of them would improve an aggregate and all of them point at someone. */
static void _write_context(JsonBuilder *b, JsonObject *profile)
{
  json_builder_set_member_name(b, "context");
  json_builder_begin_object(b);

  json_builder_set_member_name(b, "client");
  json_builder_add_string_value(b, "darktable");
  json_builder_set_member_name(b, "client_version");
  json_builder_add_string_value(b, darktable_package_version);

  /* In v2 the equipment lives under imaging_system and the frame size under
     source_layout, not at the top level. */
  JsonObject *sys = profile && json_object_has_member(profile, "imaging_system")
    ? json_object_get_object_member(profile, "imaging_system") : NULL;

  if(sys)
  {
    static const char *const carry[][2] =
      { { "id", "lens" }, { "manufacturer", "maker" },
        { "model", "camera_model" }, { "camera_type", "camera_type" },
        { NULL, NULL } };

    for(int i = 0; carry[i][0]; i++)
      if(json_object_has_member(sys, carry[i][0]))
      {
        const char *v = json_object_get_string_member(sys, carry[i][0]);
        /* darktable writes "----" when EXIF carries no lens name. Sending
           it would look like a measurement of a lens called "----". */
        if(v && *v && strcmp(v, "----"))
        {
          json_builder_set_member_name(b, carry[i][1]);
          json_builder_add_string_value(b, v);
        }
      }

    if(json_object_has_member(sys, "crop_factor"))
    {
      json_builder_set_member_name(b, "crop_factor");
      json_builder_add_double_value
        (b, json_object_get_double_member(sys, "crop_factor"));
    }
  }

  JsonObject *lay = profile && json_object_has_member(profile, "source_layout")
    ? json_object_get_object_member(profile, "source_layout") : NULL;

  if(lay)
  {
    static const char *const dims[][2] =
      { { "reference_width", "width" },
        { "reference_height", "height" }, { NULL, NULL } };

    for(int i = 0; dims[i][0]; i++)
      if(json_object_has_member(lay, dims[i][0]))
      {
        json_builder_set_member_name(b, dims[i][1]);
        json_builder_add_int_value
          (b, json_object_get_int_member(lay, dims[i][0]));
      }
  }

  json_builder_end_object(b);
}

static void _write_quality(JsonBuilder *b,
                           const dt_lens_share_quality_t *q)
{
  json_builder_set_member_name(b, "quality");
  json_builder_begin_object(b);

  if(q->have_geometry)
  {
    json_builder_set_member_name(b, "straightness_before_px");
    json_builder_add_double_value(b, q->straightness_before_px);
    json_builder_set_member_name(b, "straightness_after_px");
    json_builder_add_double_value(b, q->straightness_after_px);
    json_builder_set_member_name(b, "points_measured");
    json_builder_add_int_value(b, q->points_measured);
    json_builder_set_member_name(b, "points_stray");
    json_builder_add_int_value(b, q->points_stray);
  }

  if(q->have_vignetting)
  {
    json_builder_set_member_name(b, "vig_coverage");
    json_builder_add_double_value(b, q->vig_coverage);
  }

  json_builder_end_object(b);
}

gboolean dt_lens_share_submit(const char *profile_path,
                              const dt_lens_share_quality_t *quality,
                              GError **error)
{
  if(!profile_path || !quality) return FALSE;

  /* The profile is read back from the file just written rather than
     re-serialised from memory. Whatever ends up in the repository is then
     byte for byte what the user has locally, and there is no second writer
     that could drift from the first. */
  JsonParser *parser = json_parser_new();
  if(!json_parser_load_from_file(parser, profile_path, error))
  {
    g_object_unref(parser);
    return FALSE;
  }

  JsonNode *proot = json_parser_get_root(parser);
  JsonObject *pobj = JSON_NODE_HOLDS_OBJECT(proot)
    ? json_node_get_object(proot) : NULL;

  JsonBuilder *b = json_builder_new();
  json_builder_begin_object(b);

  json_builder_set_member_name(b, "schema");
  json_builder_add_string_value(b, "lensfit-submission-1");

  json_builder_set_member_name(b, "install_id");
  json_builder_add_string_value(b, dt_lens_share_install_id());

  /* Date only. A timestamp to the second, across a handful of submissions,
     starts to describe when someone was at their desk. */
  GDateTime *now = g_date_time_new_now_utc();
  gchar *day = g_date_time_format(now, "%Y-%m-%d");
  json_builder_set_member_name(b, "submitted");
  json_builder_add_string_value(b, day);
  g_free(day);
  g_date_time_unref(now);

  json_builder_set_member_name(b, "license");
  json_builder_add_string_value(b, "CC-BY-SA-4.0");

  _write_context(b, pobj);
  _write_quality(b, quality);

  json_builder_set_member_name(b, "profile");
  json_builder_add_value(b, json_node_copy(proot));

  json_builder_end_object(b);

  JsonGenerator *gen = json_generator_new();
  JsonNode *root = json_builder_get_root(b);
  json_generator_set_root(gen, root);
  json_generator_set_pretty(gen, TRUE);
  json_generator_set_indent(gen, 2);

  gchar *dir = dt_lens_share_queue_dir();
  gboolean ok = FALSE;

  if(dir)
  {
    gchar *stem = g_path_get_basename(profile_path);
    if(g_str_has_suffix(stem, ".json")) stem[strlen(stem) - 5] = '\0';

    /* Named after the profile, so re-saving the same lens overwrites its own
       pending submission instead of queueing the same measurement twice. */
    gchar *safe = g_strdup(stem);
    for(char *c = safe; *c; c++)
      if(!g_ascii_isalnum(*c) && *c != '-' && *c != '_') *c = '_';

    gchar *file = g_strdup_printf("%s.json", safe);
    gchar *out = g_build_filename(dir, file, NULL);

    ok = json_generator_to_file(gen, out, error);

    g_free(out);
    g_free(file);
    g_free(safe);
    g_free(stem);
    g_free(dir);
  }

  json_node_free(root);
  g_object_unref(gen);
  g_object_unref(b);
  g_object_unref(parser);

  if(ok) dt_lens_share_flush();

  return ok;
}

static size_t _discard_body(void *ptr, size_t size, size_t nmemb, void *userp)
{
  return size * nmemb;
}

/* Post one queued submission to the no-account collection service. It does,
 * on the contributor's behalf, exactly what the browser path below asks the
 * contributor to do themselves: validate, then open a pull request. Returns
 * TRUE only on an unambiguous 2xx -- anything else, including "could not
 * connect", leaves the file queued rather than guessing that it went
 * through. */
static gboolean _post_to_collector(const char *endpoint,
                                   const char *body,
                                   gsize len)
{
  CURL *curl = curl_easy_init();
  if(!curl) return FALSE;

  struct curl_slist *hdr =
    curl_slist_append(NULL, "Content-Type: application/json");

  dt_curl_init(curl, FALSE);
  curl_easy_setopt(curl, CURLOPT_URL, endpoint);
  curl_easy_setopt(curl, CURLOPT_POSTFIELDS, body);
  curl_easy_setopt(curl, CURLOPT_POSTFIELDSIZE, (long)len);
  curl_easy_setopt(curl, CURLOPT_HTTPHEADER, hdr);
  curl_easy_setopt(curl, CURLOPT_WRITEFUNCTION, _discard_body);
  curl_easy_setopt(curl, CURLOPT_TIMEOUT, 20L);

  const CURLcode res = curl_easy_perform(curl);
  long code = 0;
  curl_easy_getinfo(curl, CURLINFO_RESPONSE_CODE, &code);

  curl_slist_free_all(hdr);
  curl_easy_cleanup(curl);

  return res == CURLE_OK && code >= 200 && code < 300;
}

/* Hand GitHub a prefilled file rather than call an API. This is the primary
 * path -- see the repo README for why: a pull request costs nothing to
 * host, leaves no operator holding anyone's data, and gives every number in
 * the database an audit trail back to the measurement it came from. The
 * contributor reviews it in the browser and proposes the change themselves.
 *
 * Its one gap is that it requires a GitHub account, which the collection
 * service above exists to cover for whoever doesn't have one. */
static gboolean _open_pr_in_browser(const char *repo,
                                    const char *name,
                                    const char *body,
                                    const char *dir)
{
  gchar *stem = g_strdup(name);
  stem[strlen(stem) - 5] = '\0';

  gchar *value = g_uri_escape_string(body, NULL, FALSE);
  gchar *target = g_strdup_printf("submissions/%s/%s.json",
                                  stem, dt_lens_share_install_id());
  gchar *filename = g_uri_escape_string(target, NULL, FALSE);

  gchar *url = g_strdup_printf
    ("https://github.com/%s/new/main?filename=%s&value=%s",
     repo, filename, value);

  gboolean opened = FALSE;

  /* Browsers and servers both give up somewhere past 8k of URL. A profile
     that big cannot go this way, and silently opening a truncated one would
     be worse than saying so. */
  if(strlen(url) < 7500)
  {
    dt_open_url(url);
    opened = TRUE;
    /* Left queued deliberately: only GitHub knows whether the pull request
       was actually opened, and deleting the local copy here would throw
       the measurement away on a closed browser tab. */
    dt_control_log(_("`%s' opened in your browser to propose"), stem);
  }
  else
  {
    dt_control_log
      (_("`%s' is too large to submit through the browser -- "
         "it is queued at %s"), stem, dir);
  }

  g_free(url);
  g_free(filename);
  g_free(target);
  g_free(value);
  g_free(stem);

  return opened;
}

void dt_lens_share_flush(void)
{
  gchar *endpoint = dt_conf_get_string("plugins/lensfit/endpoint");
  gchar *repo = dt_conf_get_string("plugins/lensfit/repo");
  const gboolean have_endpoint = endpoint && *endpoint;
  const gboolean have_repo = repo && *repo;

  if(!have_endpoint && !have_repo)
  {
    g_free(endpoint);
    g_free(repo);
    return;
  }

  gchar *dir = dt_lens_share_queue_dir();
  if(!dir)
  {
    g_free(endpoint);
    g_free(repo);
    return;
  }

  GDir *d = g_dir_open(dir, 0, NULL);
  if(!d)
  {
    g_free(dir);
    g_free(endpoint);
    g_free(repo);
    return;
  }

  const gchar *name;
  while((name = g_dir_read_name(d)))
  {
    if(!g_str_has_suffix(name, ".json")) continue;

    gchar *path = g_build_filename(dir, name, NULL);
    gchar *body = NULL;
    gsize len = 0;

    if(g_file_get_contents(path, &body, &len, NULL))
    {
      if(have_endpoint && _post_to_collector(endpoint, body, len))
      {
        gchar *stem = g_strdup(name);
        stem[strlen(stem) - 5] = '\0';
        dt_control_log(_("`%s' submitted to lensfit"), stem);
        g_free(stem);
        g_unlink(path);
      }
      else if(have_repo)
      {
        /* Falls through here whenever the collection service is not
           configured, or was configured but did not accept it -- offline,
           a restart, a validation rejection worth a human's attention. The
           measurement still gets somewhere rather than being stranded. */
        _open_pr_in_browser(repo, name, body, dir);
      }

      g_free(body);
    }

    g_free(path);
    /* One at a time: a POST retries silently next flush, but opening a
       browser tab per queued file would be an ambush, not a contribution. */
    break;
  }

  g_dir_close(d);
  g_free(dir);
  g_free(endpoint);
  g_free(repo);
}

// clang-format off
// modelines: These editor modelines have been set for all relevant files by tools/update_modelines.py
// vim: shiftwidth=2 expandtab tabstop=2 cindent
// kate: tab-indents: off; indent-width 2; replace-tabs on; indent-mode cstyle; remove-trailing-spaces modified;
// clang-format on
