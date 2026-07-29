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

/* Which lens is being calibrated, and at what settings.
 *
 * First panel in the workflow, because everything below it is a measurement
 * *of* something and the something has to be named first. It also carries the
 * two values that no amount of measuring can recover: focal length and
 * aperture. A manual lens reports neither over the mount -- no contacts, so
 * the camera writes 0 for both -- and a profile with a focal length of zero
 * can never be matched to an image, which is why the database rejects it
 * outright rather than storing it and failing later.
 *
 * The panel also states where the current numbers came from. That is not
 * decoration: a measurement, a Lensfun conversion and a hand-typed set of
 * coefficients are three different kinds of claim, and only the first is
 * evidence.
 */

#include "bauhaus/bauhaus.h"
#include "common/darktable.h"
#include "common/image_cache.h"
#include "common/lens_lensfun.h"
#include "common/lens_warp.h"
#include "control/conf.h"
#include "control/control.h"
#include "gui/accelerators.h"
#include "gui/gtk.h"
#include "libs/lens_calib_common.h"
#include "libs/lib.h"
#include "libs/lib_api.h"
#include "views/view.h"

#include <string.h>

DT_MODULE(1)

typedef struct dt_lib_lens_calib_lens_t
{
  GtkWidget *maker, *model, *mount;
  GtkWidget *crop, *focal, *aperture;
  GtkWidget *focal_min, *focal_max, *aperture_min, *aperture_max;
  GtkWidget *distance_min, *distance_max;
  GtkWidget *source;

  GtkWidget *profile_maker, *profiles, *open, *new_btn;
  /* the shipped database is well over a thousand profiles -- picking a
     vendor first narrows the model list below to something browsable.
     profile_model_ids is parallel to the entries actually in d->profiles
     (no leading "none" placeholder here, unlike iop/lens.cc's picker --
     this list is never empty-selectable, so index i maps straight to
     profile_model_ids[i]). profile_all_* / profile_maker_keys are the
     cached full listing and unique vendor keys behind it, rebuilt only in
     _reload_profiles(); the maker-change handler just re-filters from the
     cache. */
  gchar **profile_all_names, **profile_all_makers, **profile_all_models;
  gchar **profile_maker_keys;
  gchar **profile_model_ids;
  GtkWidget *contents;
  GtkWidget *lf_search, *lf_results, *lf_import;
  GtkEntryCompletion *lf_completion;
  GtkListStore *lf_store;

  GtkWidget *profile_name, *save_profile, *share;

  gchar **lf_names;      // current search results, "Maker | Model"
  gchar **all_names;     // the whole catalogue, for completion
  gboolean updating;
} dt_lib_lens_calib_lens_t;

const char *name(dt_lib_module_t *self)
{
  return _("lens");
}

dt_view_type_flags_t views(dt_lib_module_t *self)
{
  return DT_VIEW_LENS_CALIB;
}

uint32_t container(dt_lib_module_t *self)
{
  return DT_UI_CONTAINER_PANEL_RIGHT_CENTER;
}

int position(const dt_lib_module_t *self)
{
  return 500;
}

/* ------------------------------------------------------------ identity */

static void _entry_changed(GtkWidget *w, dt_lib_module_t *self)
{
  dt_lib_lens_calib_lens_t *d = self->data;
  if(d->updating) return;

  const char *key = g_object_get_data(G_OBJECT(w), "conf-key");
  if(key) dt_conf_set_string(key, gtk_entry_get_text(GTK_ENTRY(w)));
}

static void _spin_changed(GtkWidget *w, dt_lib_module_t *self)
{
  dt_lib_lens_calib_lens_t *d = self->data;
  if(d->updating) return;

  const char *key = g_object_get_data(G_OBJECT(w), "conf-key");
  if(key)
    dt_conf_set_float(key, gtk_spin_button_get_value(GTK_SPIN_BUTTON(w)));
}

static GtkWidget *_entry(const char *label,
                         const char *key,
                         const char *tip,
                         GtkWidget **out,
                         dt_lib_module_t *self)
{
  *out = gtk_entry_new();
  g_object_set_data_full(G_OBJECT(*out), "conf-key", g_strdup(key), g_free);
  if(tip) gtk_widget_set_tooltip_text(*out, tip);
  g_signal_connect(G_OBJECT(*out), "changed",
                   G_CALLBACK(_entry_changed), self);
  return dt_lens_calib_labelled(label, *out);
}

static GtkWidget *_spin(const char *label,
                        const char *key,
                        const double lo,
                        const double hi,
                        const double step,
                        const int digits,
                        const char *tip,
                        GtkWidget **out,
                        dt_lib_module_t *self)
{
  *out = gtk_spin_button_new_with_range(lo, hi, step);
  gtk_spin_button_set_digits(GTK_SPIN_BUTTON(*out), digits);
  g_object_set_data_full(G_OBJECT(*out), "conf-key", g_strdup(key), g_free);
  if(tip) gtk_widget_set_tooltip_text(*out, tip);
  g_signal_connect(G_OBJECT(*out), "value-changed",
                   G_CALLBACK(_spin_changed), self);
  return dt_lens_calib_labelled(label, *out);
}

/* -------------------------------------------------- the lensfit database */

/* Rebuilds d->profiles (the model picker) from the cached full listing,
   filtered to whatever vendor is currently selected in d->profile_maker,
   and selects `preferred_name` if it survives the filter (NULL means just
   take whatever ends up first, matching the old flat list's behaviour of
   defaulting to entry 0 after a clear). Does not touch d->profile_maker
   itself. */
static void _rebuild_profile_models(dt_lib_module_t *self,
                                    const char *preferred_name)
{
  dt_lib_lens_calib_lens_t *d = self->data;

  const int maker_idx = dt_bauhaus_combobox_get(d->profile_maker);
  const char *key = (maker_idx > 0 && d->profile_maker_keys)
    ? d->profile_maker_keys[maker_idx - 1] : NULL; // NULL == "(all)"

  dt_bauhaus_combobox_clear(d->profiles);

  g_strfreev(d->profile_model_ids);
  GPtrArray *ids = g_ptr_array_new();
  int n = 0, restore = -1;

  for(int i = 0; d->profile_all_names && d->profile_all_names[i]; i++)
  {
    // case-insensitive: the source data spells the same vendor
    // differently across entries often enough ("Leica Camera AG" vs.
    // "LEICA CAMERA AG") that an exact-case filter would silently split
    // one vendor into several buckets
    if(key && g_ascii_strcasecmp(d->profile_all_makers[i], key)) continue;

    dt_bauhaus_combobox_add(d->profiles, d->profile_all_models[i]);
    g_ptr_array_add(ids, g_strdup(d->profile_all_names[i]));
    if(preferred_name && !g_strcmp0(d->profile_all_names[i], preferred_name))
      restore = n;
    n++;
  }

  g_ptr_array_add(ids, NULL);
  d->profile_model_ids = (gchar **)g_ptr_array_free(ids, FALSE);

  gtk_widget_set_sensitive(d->profiles, n > 0);
  gtk_widget_set_sensitive(d->open, n > 0);

  if(n == 0)
    dt_bauhaus_combobox_add(d->profiles, _("no saved profiles"));
  else if(restore >= 0)
    dt_bauhaus_combobox_set(d->profiles, restore);
}

static void _profile_maker_changed(GtkWidget *w, dt_lib_module_t *self)
{
  _rebuild_profile_models(self, NULL);
}

static void _reload_profiles(dt_lib_module_t *self)
{
  dt_lib_lens_calib_lens_t *d = self->data;

  /* This runs on every panel refresh -- after opening a profile, after
     adding a measurement, after anything -- not just when the list of
     saved profiles actually changes. Clearing and repopulating a bauhaus
     combobox resets its selection to entry 0, which without this would
     make the model picker visibly jump back to whichever name sorts
     first every time something elsewhere in the view changed, with no
     relation to what was actually open. Captured by actual profile name
     (via profile_model_ids), not by the combobox's displayed text --
     that text is now just the "model" portion, not the full identifier. */
  const int was_idx = dt_bauhaus_combobox_get(d->profiles);
  gchar *keep = (was_idx >= 0 && d->profile_model_ids && d->profile_model_ids[was_idx])
    ? g_strdup(d->profile_model_ids[was_idx]) : NULL;

  g_strfreev(d->profile_all_names);
  g_strfreev(d->profile_all_makers);
  g_strfreev(d->profile_all_models);
  dt_lens_profile_list_full(&d->profile_all_names, &d->profile_all_makers,
                            &d->profile_all_models);

  // unique maker keys -- the full listing is already sorted by (maker,
  // model), so identical makers are always adjacent
  g_strfreev(d->profile_maker_keys);
  GPtrArray *keys = g_ptr_array_new();
  const char *last = NULL;
  for(int i = 0; d->profile_all_names[i]; i++)
  {
    if(!last || g_ascii_strcasecmp(last, d->profile_all_makers[i]))
    {
      last = d->profile_all_makers[i];
      g_ptr_array_add(keys, g_strdup(last));
    }
  }
  g_ptr_array_add(keys, NULL);
  d->profile_maker_keys = (gchar **)g_ptr_array_free(keys, FALSE);

  dt_bauhaus_combobox_clear(d->profile_maker);
  dt_bauhaus_combobox_add(d->profile_maker, _("(all)"));
  for(int i = 0; d->profile_maker_keys[i]; i++)
    dt_bauhaus_combobox_add(d->profile_maker, d->profile_maker_keys[i][0]
                            ? d->profile_maker_keys[i] : _("(other)"));

  // preselect the vendor the kept selection actually belongs to, so a
  // refresh doesn't silently reset the filter out from under the user
  int maker_idx = 0;
  if(keep)
  {
    const char *cur_maker = NULL;
    for(int i = 0; d->profile_all_names[i]; i++)
      if(!g_strcmp0(d->profile_all_names[i], keep))
      {
        cur_maker = d->profile_all_makers[i];
        break;
      }
    if(cur_maker)
      for(int j = 0; d->profile_maker_keys[j]; j++)
        if(!g_ascii_strcasecmp(d->profile_maker_keys[j], cur_maker))
        {
          maker_idx = j + 1;
          break;
        }
  }
  dt_bauhaus_combobox_set(d->profile_maker, maker_idx);

  _rebuild_profile_models(self, keep);
  g_free(keep);
}

static void _open_clicked(GtkWidget *w, dt_lib_module_t *self)
{
  dt_lib_lens_calib_lens_t *d = self->data;
  dt_view_t *view = dt_lens_calib_get_view();
  if(!view || !darktable.view_manager->proxy.lens_calib.open_profile) return;

  const int idx = dt_bauhaus_combobox_get(d->profiles);
  const char *sel = (idx >= 0 && d->profile_model_ids) ? d->profile_model_ids[idx] : NULL;
  if(sel && *sel)
    darktable.view_manager->proxy.lens_calib.open_profile(view, sel);
}

static void _new_clicked(GtkWidget *w, dt_lib_module_t *self)
{
  dt_view_t *view = dt_lens_calib_get_view();
  if(!view || !darktable.view_manager->proxy.lens_calib.new_profile) return;

  darktable.view_manager->proxy.lens_calib.new_profile(view);
}

/* ------------------------------------------------- the Lensfun catalogue */

static gboolean _completion_match(GtkEntryCompletion *completion,
                                  const gchar *key,
                                  GtkTreeIter *iter,
                                  gpointer user_data)
{
  dt_lib_module_t *self = user_data;
  dt_lib_lens_calib_lens_t *d = self->data;

  gchar *row = NULL;
  gtk_tree_model_get(GTK_TREE_MODEL(d->lf_store), iter, 0, &row, -1);
  if(!row) return FALSE;

  /* Substring rather than prefix: nobody remembers whether the catalogue
     files a lens under its maker or its series. */
  gchar *hay = g_utf8_casefold(row, -1);
  gchar *needle = g_utf8_casefold(key, -1);
  const gboolean hit = strstr(hay, needle) != NULL;

  g_free(hay);
  g_free(needle);
  g_free(row);
  return hit;
}

static void _fill_completion(dt_lib_module_t *self)
{
  dt_lib_lens_calib_lens_t *d = self->data;
  if(d->all_names) return;

  d->all_names = dt_lf_all_lens_names();
  gtk_list_store_clear(d->lf_store);

  for(gchar **p = d->all_names; p && *p; p++)
  {
    GtkTreeIter it;
    gtk_list_store_append(d->lf_store, &it);
    gtk_list_store_set(d->lf_store, &it, 0, *p, -1);
  }
}

static gboolean _lf_focus(GtkWidget *w, GdkEvent *e, dt_lib_module_t *self)
{
  _fill_completion(self);
  return FALSE;
}

/* Populate d->lf_results from the search entry, without importing anything.
 *
 * The entry's completion is built from the same "Maker | Model" strings
 * this searches for (dt_lf_all_lens_names), so picking a suggestion leaves
 * exactly that composite string in the entry. Searching with it verbatim
 * asked Lensfun to fuzzy-match "KMZ | Helios-40 85mm f/1.5" against a
 * lens's *model* alone -- the maker name and the " | " separator
 * contaminated the match, so a name taken straight from the completion
 * list still came back "nothing found". Splitting on " | " searches maker
 * and model separately, which is what the composite string actually means;
 * free text with no separator still searches as a model, unchanged.
 */
static void _lf_run_search(dt_lib_module_t *self)
{
  dt_lib_lens_calib_lens_t *d = self->data;

  const gchar *text = gtk_entry_get_text(GTK_ENTRY(d->lf_search));

  gchar *maker = NULL, *model = NULL;
  if(text && *text)
  {
    const char *sep = strstr(text, " | ");
    if(sep)
    {
      maker = g_strndup(text, sep - text);
      model = g_strdup(sep + 3);
    }
    else
    {
      model = g_strdup(text);
    }
  }

  g_strfreev(d->lf_names);
  d->lf_names = dt_lf_search_lenses(maker, model);
  g_free(maker);
  g_free(model);

  dt_bauhaus_combobox_clear(d->lf_results);
  int n = 0;
  for(gchar **p = d->lf_names; p && *p; p++, n++)
    dt_bauhaus_combobox_add(d->lf_results, *p);

  if(n == 0)
  {
    dt_bauhaus_combobox_add(d->lf_results, _("nothing found"));
    dt_control_log(_("no lens in the Lensfun catalogue matches that"));
  }

  gtk_widget_set_sensitive(d->lf_results, n > 0);
}

/* One button rather than "search" then "import" as two separate steps.
 *
 * The completion popup already narrows candidates while typing, so a
 * second, explicit search felt redundant to a name just picked from that
 * list -- and until the fix above, it did not even work for one. This
 * collapses the common case (a name that resolves to exactly one lens,
 * whether typed in full or chosen from completion) into a single click:
 * search, and if there is exactly one result, import it immediately. An
 * ambiguous query instead populates the results list and stops, so the
 * lens can be picked by hand before pressing the button again -- pressing
 * it with a real selection already showing imports that selection without
 * searching again, so re-clicking after disambiguating does not re-run the
 * search and lose the choice.
 */
static void _lf_action(dt_lib_module_t *self)
{
  dt_lib_lens_calib_lens_t *d = self->data;
  dt_view_t *view = dt_lens_calib_get_view();

  const char *sel = dt_bauhaus_combobox_get_text(d->lf_results);
  gboolean have_result = sel && *sel && strcmp(sel, _("nothing found"));

  if(!have_result)
  {
    _lf_run_search(self);
    sel = dt_bauhaus_combobox_get_text(d->lf_results);
    have_result = sel && *sel && strcmp(sel, _("nothing found"));
    // several candidates: stop and let the dropdown be picked by hand
    if(!have_result || dt_bauhaus_combobox_length(d->lf_results) > 1) return;
  }

  if(view && darktable.view_manager->proxy.lens_calib.import_lensfun
     && sel && *sel)
    darktable.view_manager->proxy.lens_calib.import_lensfun(view, sel);
}

static void _lf_import_clicked(GtkWidget *w, dt_lib_module_t *self)
{
  _lf_action(self);
}

static void _lf_entry_activate(GtkWidget *w, dt_lib_module_t *self)
{
  _lf_action(self);
}

/* Editing the search text invalidates whatever is currently in the results
   list, so the single action button re-searches on its next click instead
   of importing a match left over from a previous, different query. */
static void _lf_entry_changed(GtkWidget *w, dt_lib_module_t *self)
{
  dt_lib_lens_calib_lens_t *d = self->data;
  dt_bauhaus_combobox_clear(d->lf_results);
  gtk_widget_set_sensitive(d->lf_results, FALSE);
}

/* --------------------------------------------------------------- saving */

// the name becomes a filename, so it cannot carry path separators
static gchar *_safe_name(const char *name)
{
  gchar *s = g_strdup(name ? name : "");
  g_strdelimit(s, "/\\:*?\"<>|", '_');
  return s;
}

static void _save_clicked(GtkWidget *w, dt_lib_module_t *self)
{
  dt_lib_lens_calib_lens_t *d = self->data;
  dt_view_t *view = dt_lens_calib_get_view();
  if(!view || !darktable.view_manager->proxy.lens_calib.save_profile) return;

  const gchar *name = gtk_entry_get_text(GTK_ENTRY(d->profile_name));
  if(!name || !*name)
  {
    dt_control_log(_("give the profile a name first"));
    return;
  }

  const gboolean share =
    gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(d->share));
  dt_conf_set_bool("plugins/lens_calib/share_profile", share);

  /* Say what will leave the machine, before it does. An opt out only means
     anything if the person had the chance to see what they would be opting
     out of -- so this keeps appearing until they say they have read it, and
     preferences can bring it back if they want checking again. */
  if(share && dt_conf_get_bool("plugins/lensfit/show_share_notice"))
  {
    GtkWidget *win = dt_ui_main_window(darktable.gui->ui);
    GtkWidget *dlg = gtk_message_dialog_new
      (GTK_WINDOW(win), GTK_DIALOG_DESTROY_WITH_PARENT, GTK_MESSAGE_INFO,
       GTK_BUTTONS_NONE,
       _("this profile will be shared with the lensfit project"));

    gtk_message_dialog_format_secondary_markup
      (GTK_MESSAGE_DIALOG(dlg),
       _("<b>sent:</b> the measured correction, the camera maker, the lens "
         "model and mount, the frame size, and how well the fit went.\n\n"
         "<b>never sent:</b> your name, your files or their paths, the "
         "camera's serial number, where the picture was taken, or any "
         "hardware identifier. submissions are told apart by a random "
         "identifier generated on this machine.\n\n"
         "shared profiles are published under CC-BY-SA-4.0. you can turn "
         "this off with the checkbox below the save button."));

    /* Ticked by default: the notice keeps coming back unless it is actively
       dismissed, rather than vanishing after one showing that may well have
       been clicked through without reading. */
    gtk_dialog_add_buttons(GTK_DIALOG(dlg),
                           _("_don't share"), GTK_RESPONSE_REJECT,
                           _("_share"), GTK_RESPONSE_ACCEPT, NULL);
    gtk_dialog_set_default_response(GTK_DIALOG(dlg), GTK_RESPONSE_ACCEPT);

    GtkWidget *again = gtk_check_button_new_with_label(_("show next time"));
    gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(again), TRUE);
    gtk_widget_set_tooltip_text
      (again,
       _("untick to stop showing this before every shared profile.\n"
         "it can be switched back on in preferences under lensfit"));
    gtk_box_pack_start
      (GTK_BOX(gtk_message_dialog_get_message_area(GTK_MESSAGE_DIALOG(dlg))),
       again, FALSE, FALSE, 0);
    gtk_widget_show_all(dlg);

    const gint r = gtk_dialog_run(GTK_DIALOG(dlg));
    const gboolean show_again =
      gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(again));
    gtk_widget_destroy(dlg);

    /* Unticking only takes effect for someone who went ahead. If they
       declined, the notice stays on: it is what would explain the setting
       again the day they turn sharing back on. */
    dt_conf_set_bool("plugins/lensfit/show_share_notice",
                     show_again || r != GTK_RESPONSE_ACCEPT);

    if(r != GTK_RESPONSE_ACCEPT)
    {
      /* Declining here unticks the checkbox itself, since that checkbox is
         now the only place the choice lives. The notice stays on though --
         it is what would explain the setting again the next time it is
         ticked. */
      dt_conf_set_bool("plugins/lens_calib/share_profile", FALSE);
      gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(d->share), FALSE);
    }
  }

  gchar *safe = _safe_name(name);
  dt_conf_set_string("plugins/lens_calib/profile_name", safe);
  darktable.view_manager->proxy.lens_calib.save_profile(view, safe);
  g_free(safe);
}

/* --------------------------------------------------------------- update */

static void _refresh(dt_lib_module_t *self)
{
  dt_lib_lens_calib_lens_t *d = self->data;
  dt_view_t *view = dt_lens_calib_get_view();

  const gboolean solved =
    view && darktable.view_manager->proxy.lens_calib.has_solution
    && darktable.view_manager->proxy.lens_calib.has_solution(view);

  gtk_widget_set_sensitive(d->save_profile, solved);
}

static const char *_source_label(const dt_lens_source_t s)
{
  switch(s)
  {
    case DT_LENS_SOURCE_MEASURED:     return _("measured here");
    case DT_LENS_SOURCE_MANUFACTURER: return _("from the manufacturer");
    case DT_LENS_SOURCE_LENSFUN:      return _("converted from Lensfun");
    case DT_LENS_SOURCE_AGGREGATED:   return _("aggregated from submissions");
    case DT_LENS_SOURCE_REVERSE_ENG:  return _("reverse engineered");
    case DT_LENS_SOURCE_EDITED:       return _("measured, then edited by hand");
  }
  return "";
}

void gui_update(dt_lib_module_t *self)
{
  dt_lib_lens_calib_lens_t *d = self->data;

  d->updating = TRUE;

  static const struct { GtkWidget **w; const char *key; } entries[] = {
    { NULL, NULL }
  };
  (void)entries;

  struct { GtkWidget *w; const char *key; } text[] = {
    { d->maker, "plugins/lens_calib/maker" },
    { d->model, "plugins/lens_calib/model" },
    { d->mount, "plugins/lens_calib/mount" },
    { d->profile_name, "plugins/lens_calib/profile_name" },
  };

  for(size_t i = 0; i < sizeof(text) / sizeof(text[0]); i++)
  {
    gchar *v = dt_conf_get_string(text[i].key);
    gtk_entry_set_text(GTK_ENTRY(text[i].w), v ? v : "");
    g_free(v);
  }

  gtk_spin_button_set_value(GTK_SPIN_BUTTON(d->crop),
                            dt_conf_get_float("plugins/lens_calib/crop_factor"));
  gtk_spin_button_set_value(GTK_SPIN_BUTTON(d->focal),
                            dt_conf_get_float("plugins/lens_calib/focal"));
  gtk_spin_button_set_value(GTK_SPIN_BUTTON(d->aperture),
                            dt_conf_get_float("plugins/lens_calib/aperture"));

  /* The lens's own range: a typed override if there is one, else whatever
     the session's measurements currently span. d->updating guards
     _spin_changed above, so showing the observed value here does not turn
     it into a false override -- the whole point is that this keeps
     tracking new measurements until someone actually types over it. */
  {
    dt_view_t *rv = dt_lens_calib_get_view();
    struct { GtkWidget *w; const char *conf; int axis; gboolean want_max; }
    ranges[] = {
      { d->focal_min,    "plugins/lens_calib/focal_min",    0, FALSE },
      { d->focal_max,    "plugins/lens_calib/focal_max",    0, TRUE  },
      { d->aperture_min, "plugins/lens_calib/aperture_min", 1, FALSE },
      { d->aperture_max, "plugins/lens_calib/aperture_max", 1, TRUE  },
      { d->distance_min, "plugins/lens_calib/distance_min", 2, FALSE },
      { d->distance_max, "plugins/lens_calib/distance_max", 2, TRUE  },
    };
    for(size_t i = 0; i < sizeof(ranges) / sizeof(ranges[0]); i++)
    {
      float v = dt_conf_get_float(ranges[i].conf);
      if(v <= 0.0f && rv
         && darktable.view_manager->proxy.lens_calib.measurement_range)
      {
        float lo = 0.0f, hi = 0.0f;
        if(darktable.view_manager->proxy.lens_calib
             .measurement_range(rv, ranges[i].axis, &lo, &hi))
          v = ranges[i].want_max ? hi : lo;
      }
      gtk_spin_button_set_value(GTK_SPIN_BUTTON(ranges[i].w), v);
    }
  }

  const dt_lens_source_t src =
    (dt_lens_source_t)dt_conf_get_int("plugins/lens_calib/source");

  gchar *parent = dt_conf_get_string("plugins/lens_calib/parent");
  gchar *txt = (parent && *parent)
    ? g_strdup_printf("%s\n%s", _source_label(src), parent)
    : g_strdup(_source_label(src));
  gtk_label_set_text(GTK_LABEL(d->source), txt);
  g_free(txt);
  g_free(parent);

  gtk_toggle_button_set_active
    (GTK_TOGGLE_BUTTON(d->share),
     dt_conf_get_bool("plugins/lens_calib/share_profile"));

  _reload_profiles(self);
  _refresh(self);

  dt_view_t *view = dt_lens_calib_get_view();
  if(view && darktable.view_manager->proxy.lens_calib.contents_summary)
  {
    gchar *summary =
      darktable.view_manager->proxy.lens_calib.contents_summary(view);
    gtk_label_set_text(GTK_LABEL(d->contents), summary);
    g_free(summary);
  }
  else
  {
    gtk_label_set_text(GTK_LABEL(d->contents), "");
  }

  d->updating = FALSE;
}

/* Fill in whatever EXIF does know, without overwriting anything typed.
 *
 * Called when the panel is built rather than on every update, so a value
 * corrected by hand is not put back to the camera's version on the next
 * refresh. */
static void _prefill_from_exif(void)
{
  dt_view_t *view = dt_lens_calib_get_view();
  if(!view) return;

  const dt_lib_module_t *m = NULL;
  (void)m;

  dt_imgid_t imgid = dt_act_on_get_main_image();
  if(!dt_is_valid_imgid(imgid)) return;

  const dt_image_t *img = dt_image_cache_get(imgid, 'r');
  if(!img) return;

  if(dt_conf_get_float("plugins/lens_calib/focal") <= 0.0f
     && img->exif_focal_length > 0.0f)
    dt_conf_set_float("plugins/lens_calib/focal", img->exif_focal_length);

  if(dt_conf_get_float("plugins/lens_calib/aperture") <= 0.0f
     && img->exif_aperture > 0.0f)
    dt_conf_set_float("plugins/lens_calib/aperture", img->exif_aperture);

  gchar *maker = dt_conf_get_string("plugins/lens_calib/maker");
  if((!maker || !*maker) && img->exif_maker[0])
    dt_conf_set_string("plugins/lens_calib/maker", img->exif_maker);
  g_free(maker);

  /* darktable writes "----" when the camera reports no lens, which is
     exactly the case this panel exists for -- so it is not a name worth
     copying in. */
  gchar *model = dt_conf_get_string("plugins/lens_calib/model");
  if((!model || !*model) && img->exif_lens[0]
     && strcmp(img->exif_lens, "----"))
    dt_conf_set_string("plugins/lens_calib/model", img->exif_lens);
  g_free(model);

  dt_image_cache_read_release(img);
}

void gui_init(dt_lib_module_t *self)
{
  dt_lib_lens_calib_lens_t *d = g_malloc0(sizeof(*d));
  self->data = d;

  self->widget = dt_gui_vbox();

  /* --- open an existing profile --- */
  /* the shipped database alone is well over a thousand entries -- picking
     a vendor first narrows the model list below to something browsable */
  d->profile_maker = dt_bauhaus_combobox_new(NULL);
  dt_bauhaus_widget_set_label(d->profile_maker, NULL, N_("lens vendor"));
  dt_bauhaus_combobox_set_selected_text_align
    (d->profile_maker, DT_BAUHAUS_COMBOBOX_ALIGN_LEFT);
  gtk_widget_set_tooltip_text(d->profile_maker,
                              _("narrows the model list below to one vendor"));
  g_signal_connect(G_OBJECT(d->profile_maker), "value-changed",
                   G_CALLBACK(_profile_maker_changed), self);
  dt_gui_box_add(self->widget, d->profile_maker);

  d->profiles = dt_bauhaus_combobox_new(NULL);
  dt_bauhaus_widget_set_label(d->profiles, NULL, N_("model"));
  dt_bauhaus_combobox_set_selected_text_align
    (d->profiles, DT_BAUHAUS_COMBOBOX_ALIGN_LEFT);
  gtk_widget_set_tooltip_text(d->profiles,
                              _("profiles already in your lensfit database"));

  d->open = dt_action_button_new
    (self, N_("open profile"), _open_clicked, self,
     _("load this profile as the current fit, to inspect it or to add\n"
       "another focal length to it"), 0, (GdkModifierType)0);

  d->new_btn = dt_action_button_new
    (self, N_("new"), _new_clicked, self,
     _("forget the profile currently open or being built, so the next\n"
       "save writes a new file instead of adding to it.\n\n"
       "the chart detection, points and pose are not affected -- this\n"
       "only clears which lens is being profiled, not the photograph\n"
       "being measured"), 0, (GdkModifierType)0);

  dt_gui_box_add(self->widget, d->profiles, d->open, d->new_btn);

  /* --- import from Lensfun --- */
  d->lf_search = gtk_entry_new();
  gtk_entry_set_placeholder_text(GTK_ENTRY(d->lf_search),
                                 _("lens name, or empty for all"));
  gtk_widget_set_tooltip_text
    (d->lf_search,
     _("search the Lensfun catalogue by lens name.\n"
       "an imported calibration is somebody else's measurement:\n"
       "it is marked as converted from Lensfun and keeps that licence"));

  d->lf_store = gtk_list_store_new(1, G_TYPE_STRING);
  d->lf_completion = gtk_entry_completion_new();
  gtk_entry_completion_set_model(d->lf_completion,
                                 GTK_TREE_MODEL(d->lf_store));
  gtk_entry_completion_set_text_column(d->lf_completion, 0);
  gtk_entry_completion_set_minimum_key_length(d->lf_completion, 2);
  gtk_entry_completion_set_inline_completion(d->lf_completion, FALSE);
  gtk_entry_completion_set_match_func(d->lf_completion, _completion_match,
                                      self, NULL);
  gtk_entry_set_completion(GTK_ENTRY(d->lf_search), d->lf_completion);
  g_object_unref(d->lf_completion);

  g_signal_connect(G_OBJECT(d->lf_search), "activate",
                   G_CALLBACK(_lf_entry_activate), self);
  g_signal_connect(G_OBJECT(d->lf_search), "changed",
                   G_CALLBACK(_lf_entry_changed), self);
  g_signal_connect(G_OBJECT(d->lf_search), "focus-in-event",
                   G_CALLBACK(_lf_focus), self);

  d->lf_results = dt_bauhaus_combobox_new(NULL);
  dt_bauhaus_widget_set_label(d->lf_results, NULL, N_("Lensfun match"));
  dt_bauhaus_combobox_set_selected_text_align
    (d->lf_results, DT_BAUHAUS_COMBOBOX_ALIGN_LEFT);
  gtk_widget_set_sensitive(d->lf_results, FALSE);

  /* One button rather than a separate search step -- see _lf_action. */
  d->lf_import = dt_action_button_new
    (self, N_("import from Lensfun"), _lf_import_clicked, self,
     _("type or pick a lens above, then click this.\n\n"
       "an unambiguous name imports straight away; several matches are\n"
       "listed below to choose from first, then click again to import\n"
       "the one selected.\n\n"
       "the result is not a measurement and is recorded as such"), 0,
     (GdkModifierType)0);

  GtkWidget *lf_box = dt_gui_vbox(d->lf_search, d->lf_results, d->lf_import);
  GtkWidget *lf_exp = gtk_expander_new(_("import from Lensfun"));
  gtk_container_add(GTK_CONTAINER(lf_exp), lf_box);
  dt_gui_box_add(self->widget, lf_exp);

  /* --- identity --- */
  dt_gui_box_add
    (self->widget,
     _entry(_("maker"), "plugins/lens_calib/maker",
            _("who made the lens"), &d->maker, self),
     _entry(_("model"), "plugins/lens_calib/model",
            _("the lens being calibrated"), &d->model, self),
     _entry(_("mount"), "plugins/lens_calib/mount",
            _("the mount it attaches by, for example E, RF or PL"),
            &d->mount, self));

  dt_gui_box_add
    (self->widget,
     _spin(_("crop factor"), "plugins/lens_calib/crop_factor",
           0.1, 100.0, 0.01, 3,
           _("crop factor of the sensor this was calibrated on"),
           &d->crop, self),
     _spin(_("focal length"), "plugins/lens_calib/focal",
           0.0, 2000.0, 1.0, 1,
           _("the focal length this calibration was shot at, in mm.\n\n"
             "taken from EXIF when the camera reports one. a manual lens\n"
             "reports nothing over the mount, so this has to be typed --\n"
             "and a profile without a focal length can never be matched\n"
             "to an image"),
           &d->focal, self),
     _spin(_("aperture"), "plugins/lens_calib/aperture",
           0.0, 128.0, 0.1, 1,
           _("the aperture this calibration was shot at.\n\n"
             "vignetting depends on aperture more than on anything else,\n"
             "so a falloff measurement without one cannot be selected\n"
             "correctly later"),
           &d->aperture, self));

  /* --- the lens's own range, as opposed to the single measurement above --- */
  GtkWidget *props_box = dt_gui_vbox();
  dt_gui_box_add
    (props_box,
     _spin(_("focal min"), "plugins/lens_calib/focal_min",
           0.0, 2000.0, 1.0, 1,
           _("the lens's own shortest focal length, in mm -- 24 for a\n"
             "24-70mm zoom, not the shortest one actually measured so far.\n\n"
             "type a value to fix it; left at zero, it tracks the shortest\n"
             "focal length measured in this session so far, and keeps\n"
             "tracking it as more are added"),
           &d->focal_min, self),
     _spin(_("focal max"), "plugins/lens_calib/focal_max",
           0.0, 2000.0, 1.0, 1,
           _("the lens's own longest focal length, in mm.\n\n"
             "type a value to fix it; left at zero, it tracks the longest\n"
             "focal length measured in this session so far"),
           &d->focal_max, self));
  dt_gui_box_add
    (props_box,
     _spin(_("aperture min"), "plugins/lens_calib/aperture_min",
           0.0, 128.0, 0.1, 1,
           _("the lens's own widest aperture, as an f-number -- 2.8 for an\n"
             "f/2.8 lens, not the widest one a vignetting measurement\n"
             "happens to have used so far.\n\n"
             "type a value to fix it; left at zero, it tracks what has\n"
             "actually been measured"),
           &d->aperture_min, self),
     _spin(_("aperture max"), "plugins/lens_calib/aperture_max",
           0.0, 128.0, 0.1, 1,
           _("the lens's own smallest aperture, as an f-number.\n\n"
             "type a value to fix it; left at zero, it tracks what has\n"
             "actually been measured"),
           &d->aperture_max, self));
  dt_gui_box_add
    (props_box,
     _spin(_("focus min"), "plugins/lens_calib/distance_min",
           0.0, 1000.0, 0.01, 2,
           _("the lens's own closest focus distance, in metres.\n\n"
             "type a value to fix it; left at zero, it tracks what has\n"
             "actually been measured"),
           &d->distance_min, self),
     _spin(_("focus max"), "plugins/lens_calib/distance_max",
           0.0, 1000.0, 1.0, 1,
           _("the lens's own farthest focus distance, in metres. leave at\n"
             "a large value for a lens that focuses to infinity.\n\n"
             "type a value to fix it; left at zero, it tracks what has\n"
             "actually been measured"),
           &d->distance_max, self));
  GtkWidget *props_exp = gtk_expander_new(_("lens properties"));
  gtk_expander_set_expanded(GTK_EXPANDER(props_exp), TRUE);
  gtk_container_add(GTK_CONTAINER(props_exp), props_box);
  dt_gui_box_add(self->widget, props_exp);

  /* --- what the loaded profile already holds --- */
  d->contents = gtk_label_new("");
  gtk_label_set_xalign(GTK_LABEL(d->contents), 0.0);
  gtk_label_set_line_wrap(GTK_LABEL(d->contents), TRUE);
  gtk_widget_set_tooltip_text
    (d->contents,
     _("what the current session actually contains: which focal lengths\n"
       "have a geometry fit, and which apertures have a vignetting\n"
       "measurement. updates as you open a profile, import one, or add\n"
       "another measurement in the fit panel below"));
  dt_gui_box_add(self->widget, d->contents);

  /* --- where the numbers came from --- */
  d->source = gtk_label_new("");
  gtk_label_set_xalign(GTK_LABEL(d->source), 0.0);
  gtk_label_set_line_wrap(GTK_LABEL(d->source), TRUE);
  gtk_widget_set_tooltip_text
    (d->source,
     _("what this profile claims to be. only an untouched fit counts as a\n"
       "measurement -- an import or a hand edit is recorded differently,\n"
       "because agreement with a copy of someone else's measurement is\n"
       "not corroboration"));
  dt_gui_box_add(self->widget, d->source);

  /* --- saving --- */
  d->profile_name = gtk_entry_new();
  gtk_entry_set_placeholder_text(GTK_ENTRY(d->profile_name),
                                 _("profile name"));
  g_object_set_data_full(G_OBJECT(d->profile_name), "conf-key",
                         g_strdup("plugins/lens_calib/profile_name"), g_free);
  g_signal_connect(G_OBJECT(d->profile_name), "changed",
                   G_CALLBACK(_entry_changed), self);
  gtk_widget_set_tooltip_text
    (d->profile_name,
     _("name to save the lens under in the vault.\n"
       "saving again with the same name adds this focal length to the\n"
       "existing profile rather than replacing it, which is how a zoom\n"
       "gets calibrated across its range"));
  dt_gui_box_add(self->widget, d->profile_name);

  d->save_profile = dt_action_button_new
    (self, N_("save profile"), _save_clicked, self,
     _("store the fit so the lens correction module can use it"), 0,
     (GdkModifierType)0);
  dt_gui_box_add(self->widget, d->save_profile);

  d->share = gtk_check_button_new_with_label(_("share profile with lensfit"));
  gtk_widget_set_tooltip_text
    (d->share,
     _("contribute this measurement to the open lensfit database, so other\n"
       "people get the correction for this lens the way you get theirs.\n\n"
       "the correction, the camera maker, the lens model and mount, the\n"
       "frame size and the quality of the fit are sent. your name, your\n"
       "files, the camera's serial number and where the picture was taken\n"
       "are not, and neither is any hardware identifier.\n\n"
       "shared profiles are published under CC-BY-SA-4.0.\n"
       "this checkbox is the only place sharing is turned on or off,\n"
       "and it remembers its state for the next profile you save"));
  dt_gui_box_add(self->widget, d->share);

  _prefill_from_exif();
  gui_update(self);
}

void gui_cleanup(dt_lib_module_t *self)
{
  dt_lib_lens_calib_lens_t *d = self->data;
  g_strfreev(d->lf_names);
  g_strfreev(d->all_names);
  g_strfreev(d->profile_all_names);
  g_strfreev(d->profile_all_makers);
  g_strfreev(d->profile_all_models);
  g_strfreev(d->profile_maker_keys);
  g_strfreev(d->profile_model_ids);
  g_free(self->data);
  self->data = NULL;
}

// clang-format off
// modelines: These editor modelines have been set for all relevant files by tools/update_modelines.py
// vim: shiftwidth=2 expandtab tabstop=2 cindent
// kate: tab-indents: off; indent-width 2; replace-tabs on; indent-mode cstyle; remove-trailing-spaces modified;
// clang-format on
