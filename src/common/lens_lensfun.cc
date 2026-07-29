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

/* C++ because Lensfun's real API is C++; the header this implements is C so
 * the view and the settings panel can use it without becoming C++. */

#include "common/lens_lensfun.h"
#include "common/darktable.h"
#include "common/file_location.h"
#include "common/lens_warp.h"

#include <lensfun.h>
#include <math.h>
#include <stdarg.h>
#include <string.h>

/* One database for the whole process, opened on first use.
 *
 * Deliberately separate from the copy the lens iop keeps in its global
 * data: that one is private to the module and not reachable from here, and
 * loading the database twice costs a few megabytes rather than any
 * correctness. Guarded because callers arrive from gui callbacks.
 */
static lfDatabase *_db = NULL;
// GMutex may be statically allocated and needs no explicit init
static GMutex _db_lock;

static lfDatabase *_get_db(void)
{
  if(_db) return _db;

  lfDatabase *db = new lfDatabase;

  /* Same fallback ladder the lens iop uses: the system location first, and
     if that fails -- which it does on installs whose Lensfun was built with
     a path that does not resolve at runtime -- the copy shipped beside our
     own data directory. */
  if(db->Load() != LF_NO_ERROR)
  {
    char datadir[PATH_MAX] = { 0 };
    dt_loc_get_datadir(datadir, sizeof(datadir));

    GFile *file = g_file_parse_name(datadir);
    GFile *parent = g_file_get_parent(file);
    gchar *path = parent ? g_file_get_path(parent) : NULL;
    if(parent) g_object_unref(parent);
    g_object_unref(file);

    gboolean ok = FALSE;
    if(path)
    {
#ifdef LF_MAX_DATABASE_VERSION
      gchar *versioned = g_strdup_printf("%s/lensfun/version_%d",
                                         path, LF_MAX_DATABASE_VERSION);
      g_free(db->HomeDataDir);
      db->HomeDataDir = g_strdup(versioned);
      ok = db->Load() == LF_NO_ERROR;
      g_free(versioned);
#endif
      if(!ok)
      {
        g_free(db->HomeDataDir);
        db->HomeDataDir = g_build_filename(path, "lensfun", NULL);
        ok = db->Load() == LF_NO_ERROR;
      }
      g_free(path);
    }

    if(!ok)
      dt_print(DT_DEBUG_ALWAYS,
               "[lens_lensfun] could not load the Lensfun database");
  }

  _db = db;
  return _db;
}

const char *dt_lf_home_data_dir(void)
{
  g_mutex_lock(&_db_lock);
  lfDatabase *db = _get_db();
  const char *dir = db ? db->HomeDataDir : NULL;
  g_mutex_unlock(&_db_lock);
  return dir;
}

void dt_lf_profile_init(dt_lf_profile_t *p)
{
  if(!p) return;
  memset(p, 0, sizeof(*p));
  p->dist = g_array_new(FALSE, FALSE, sizeof(dt_lf_dist_entry_t));
  p->tca = g_array_new(FALSE, FALSE, sizeof(dt_lf_tca_entry_t));
  p->vign = g_array_new(FALSE, FALSE, sizeof(dt_lf_vign_entry_t));
  p->crop_factor = 1.0f;
  p->aspect_ratio = 1.5f;
}

void dt_lf_profile_cleanup(dt_lf_profile_t *p)
{
  if(!p) return;
  if(p->dist) g_array_free(p->dist, TRUE);
  if(p->tca) g_array_free(p->tca, TRUE);
  if(p->vign) g_array_free(p->vign, TRUE);
  p->dist = p->tca = p->vign = NULL;
}

gchar **dt_lf_search_lenses(const char *maker, const char *model)
{
  g_mutex_lock(&_db_lock);
  lfDatabase *db = _get_db();

  const lfLens **lenses = db->FindLenses(NULL, maker, model,
                                         LF_SEARCH_SORT_AND_UNIQUIFY);
  GPtrArray *out = g_ptr_array_new();

  if(lenses)
  {
    for(int i = 0; lenses[i]; i++)
    {
      const char *lm = lf_mlstr_get(lenses[i]->Maker);
      const char *ln = lf_mlstr_get(lenses[i]->Model);
      if(!ln) continue;
      g_ptr_array_add(out, g_strdup_printf("%s | %s", lm ? lm : "", ln));
    }
    lf_free(lenses);
  }

  g_mutex_unlock(&_db_lock);

  g_ptr_array_add(out, NULL);
  return (gchar **)g_ptr_array_free(out, FALSE);
}

gchar **dt_lf_all_lens_names(void)
{
  g_mutex_lock(&_db_lock);
  lfDatabase *db = _get_db();

  GPtrArray *out = g_ptr_array_new();

  /* GetLenses hands back the database's own list, which -- unlike the
     result of FindLenses -- must not be freed. */
  const lfLens *const *lenses = db->GetLenses();
  for(int i = 0; lenses && lenses[i]; i++)
  {
    const char *lm = lf_mlstr_get(lenses[i]->Maker);
    const char *ln = lf_mlstr_get(lenses[i]->Model);
    if(!ln) continue;
    g_ptr_array_add(out, g_strdup_printf("%s | %s", lm ? lm : "", ln));
  }

  g_mutex_unlock(&_db_lock);

  g_ptr_array_add(out, NULL);
  return (gchar **)g_ptr_array_free(out, FALSE);
}

static dt_lf_dist_model_t _map_dist(const int m)
{
  switch(m)
  {
    case LF_DIST_MODEL_POLY3:  return DT_LF_DIST_POLY3;
    case LF_DIST_MODEL_POLY5:  return DT_LF_DIST_POLY5;
    case LF_DIST_MODEL_PTLENS: return DT_LF_DIST_PTLENS;
    default:                   return DT_LF_DIST_NONE;
  }
}

static dt_lf_tca_model_t _map_tca(const int m)
{
  switch(m)
  {
    case LF_TCA_MODEL_LINEAR: return DT_LF_TCA_LINEAR;
    case LF_TCA_MODEL_POLY3:  return DT_LF_TCA_POLY3;
    default:                  return DT_LF_TCA_NONE;
  }
}

gboolean dt_lf_load_profile(const char *model, dt_lf_profile_t *out)
{
  if(!model || !*model || !out) return FALSE;

  g_mutex_lock(&_db_lock);
  lfDatabase *db = _get_db();
  const lfLens **lenses = db->FindLenses(NULL, NULL, model, 0);

  if(!lenses || !lenses[0])
  {
    if(lenses) lf_free(lenses);
    g_mutex_unlock(&_db_lock);
    return FALSE;
  }

  const lfLens *l = lenses[0];

  const char *maker = lf_mlstr_get(l->Maker);
  const char *name = lf_mlstr_get(l->Model);
  if(maker) g_strlcpy(out->maker, maker, sizeof(out->maker));
  if(name) g_strlcpy(out->model, name, sizeof(out->model));
  if(l->Mounts && l->Mounts[0])
    g_strlcpy(out->mount, l->Mounts[0], sizeof(out->mount));

  out->min_focal = l->MinFocal;
  out->max_focal = l->MaxFocal;
  out->min_aperture = l->MinAperture;
  out->max_aperture = l->MaxAperture;
  out->crop_factor = l->CropFactor;
  out->aspect_ratio = l->AspectRatio;
  out->centre_x = l->CenterX;
  out->centre_y = l->CenterY;
  out->type = (int)l->Type;

  /* The calibration arrays are NULL terminated lists of pointers, one
     entry per focal length -- and for vignetting per focal, aperture and
     distance combination. Copy them all: the point of importing is to see
     everything the profile actually knows. */
  if(l->CalibDistortion)
  {
    for(int i = 0; l->CalibDistortion[i]; i++)
    {
      const lfLensCalibDistortion *c = l->CalibDistortion[i];
      dt_lf_dist_entry_t e;
      memset(&e, 0, sizeof(e));
      e.model = _map_dist((int)c->Model);
      e.focal = c->Focal;
      for(int t = 0; t < 3; t++) e.terms[t] = c->Terms[t];
      g_array_append_val(out->dist, e);
    }
  }

  if(l->CalibTCA)
  {
    for(int i = 0; l->CalibTCA[i]; i++)
    {
      const lfLensCalibTCA *c = l->CalibTCA[i];
      dt_lf_tca_entry_t e;
      memset(&e, 0, sizeof(e));
      e.model = _map_tca((int)c->Model);
      e.focal = c->Focal;
      for(int t = 0; t < 6; t++) e.terms[t] = c->Terms[t];
      g_array_append_val(out->tca, e);
    }
  }

  if(l->CalibVignetting)
  {
    for(int i = 0; l->CalibVignetting[i]; i++)
    {
      const lfLensCalibVignetting *c = l->CalibVignetting[i];
      dt_lf_vign_entry_t e;
      memset(&e, 0, sizeof(e));
      e.focal = c->Focal;
      e.aperture = c->Aperture;
      e.distance = c->Distance;
      for(int t = 0; t < 3; t++) e.k[t] = c->Terms[t];
      g_array_append_val(out->vign, e);
    }
  }

  lf_free(lenses);
  g_mutex_unlock(&_db_lock);

  dt_print(DT_DEBUG_ALWAYS,
           "[lens_lensfun] loaded `%s': %u distortion, %u TCA,"
           " %u vignetting entries",
           out->model, out->dist->len, out->tca->len, out->vign->len);

  return TRUE;
}

/* Numbers in the file have to use a dot regardless of the user's locale,
   so every value goes through g_ascii_formatd rather than printf. */
static void _append_num(GString *s, const char *attr, const float v)
{
  char buf[G_ASCII_DTOSTR_BUF_SIZE];
  g_ascii_formatd(buf, sizeof(buf), "%g", (double)v);
  g_string_append_printf(s, " %s=\"%s\"", attr, buf);
}

static void _append_esc(GString *s, const char *tag, const char *value)
{
  if(!value || !*value) return;
  gchar *esc = g_markup_escape_text(value, -1);
  g_string_append_printf(s, "        <%s>%s</%s>\n", tag, esc, tag);
  g_free(esc);
}

gboolean dt_lf_write_xml(const dt_lf_profile_t *p,
                         const char *path,
                         GError **error)
{
  if(!p || !path) return FALSE;

  GString *s = g_string_new(NULL);

#ifdef LF_MAX_DATABASE_VERSION
  const int version = LF_MAX_DATABASE_VERSION;
#else
  const int version = 2;
#endif

  g_string_append(s, "<!-- written by darktable lens calibration -->\n");
  g_string_append_printf(s, "<lensdatabase version=\"%d\">\n", version);
  g_string_append(s, "    <lens>\n");

  _append_esc(s, "maker", p->maker);
  _append_esc(s, "model", p->model);
  _append_esc(s, "mount", p->mount);

  g_string_append(s, "        <focal");
  _append_num(s, "min", p->min_focal);
  _append_num(s, "max", p->max_focal > p->min_focal ? p->max_focal
                                                    : p->min_focal);
  g_string_append(s, " />\n");

  if(p->min_aperture > 0.0f)
  {
    g_string_append(s, "        <aperture");
    _append_num(s, "min", p->min_aperture);
    _append_num(s, "max", p->max_aperture > p->min_aperture ? p->max_aperture
                                                            : p->min_aperture);
    g_string_append(s, " />\n");
  }

  g_string_append(s, "        <cropfactor>");
  {
    char buf[G_ASCII_DTOSTR_BUF_SIZE];
    g_ascii_formatd(buf, sizeof(buf), "%g", (double)p->crop_factor);
    g_string_append(s, buf);
  }
  g_string_append(s, "</cropfactor>\n");

  if(p->aspect_ratio > 0.0f)
  {
    g_string_append(s, "        <aspect-ratio>");
    char buf[G_ASCII_DTOSTR_BUF_SIZE];
    g_ascii_formatd(buf, sizeof(buf), "%g", (double)p->aspect_ratio);
    g_string_append(s, buf);
    g_string_append(s, "</aspect-ratio>\n");
  }

  /* Lensfun treats an absent <center> as 0,0, so only write one that says
     something. */
  if(p->centre_x != 0.0f || p->centre_y != 0.0f)
  {
    g_string_append(s, "        <center");
    _append_num(s, "x", p->centre_x);
    _append_num(s, "y", p->centre_y);
    g_string_append(s, " />\n");
  }

  if(p->type > 0)
    g_string_append_printf(s, "        <type>%s</type>\n",
                           dt_lf_type_name(p->type));

  const guint entries = (p->dist ? p->dist->len : 0)
                      + (p->tca ? p->tca->len : 0)
                      + (p->vign ? p->vign->len : 0);

  if(entries)
  {
    g_string_append(s, "        <calibration>\n");

    for(guint i = 0; p->dist && i < p->dist->len; i++)
    {
      const dt_lf_dist_entry_t *e =
        &g_array_index(p->dist, dt_lf_dist_entry_t, i);
      if(e->model == DT_LF_DIST_NONE) continue;
      /* Lensfun requires a focal length on every entry and has no "any"
         sentinel; writing 0 would not omit the claim, it would make a real
         and wrong one -- lensfun reads focal="0" as an actual 0mm
         calibration and selects it for every wide shot. Skip instead. */
      if(e->focal <= 0.0f) continue;

      const char *const *names = dt_lf_dist_term_names(e->model);
      g_string_append_printf(s, "            <distortion model=\"%s\"",
                             dt_lf_dist_model_name(e->model));
      _append_num(s, "focal", e->focal);
      for(int k = 0; names[k]; k++) _append_num(s, names[k], e->terms[k]);
      g_string_append(s, " />\n");
    }

    for(guint i = 0; p->tca && i < p->tca->len; i++)
    {
      const dt_lf_tca_entry_t *e =
        &g_array_index(p->tca, dt_lf_tca_entry_t, i);
      if(e->model == DT_LF_TCA_NONE) continue;
      if(e->focal <= 0.0f) continue; // see the distortion loop above

      const char *const *names = dt_lf_tca_term_names(e->model);
      g_string_append_printf(s, "            <tca model=\"%s\"",
                             dt_lf_tca_model_name(e->model));
      _append_num(s, "focal", e->focal);
      for(int k = 0; names[k]; k++) _append_num(s, names[k], e->terms[k]);
      g_string_append(s, " />\n");
    }

    for(guint i = 0; p->vign && i < p->vign->len; i++)
    {
      const dt_lf_vign_entry_t *e =
        &g_array_index(p->vign, dt_lf_vign_entry_t, i);
      if(e->focal <= 0.0f) continue; // see the distortion loop above

      g_string_append(s, "            <vignetting model=\"pa\"");
      _append_num(s, "focal", e->focal);
      _append_num(s, "aperture", e->aperture);
      _append_num(s, "distance", e->distance);
      _append_num(s, "k1", e->k[0]);
      _append_num(s, "k2", e->k[1]);
      _append_num(s, "k3", e->k[2]);
      g_string_append(s, " />\n");
    }

    g_string_append(s, "        </calibration>\n");
  }

  g_string_append(s, "    </lens>\n");
  g_string_append(s, "</lensdatabase>\n");

  const gboolean ok = g_file_set_contents(path, s->str, s->len, error);
  g_string_free(s, TRUE);

  if(ok)
    dt_print(DT_DEBUG_ALWAYS, "[lens_lensfun] wrote `%s'", path);

  return ok;
}

static void _note(GPtrArray *out, const char *fmt, ...)
{
  if(!out) return;
  va_list ap;
  va_start(ap, fmt);
  gchar *msg = g_strdup_vprintf(fmt, ap);
  va_end(ap);
  g_ptr_array_add(out, msg);
}

// ptlens's a,b,c and poly3's/poly5's k1,k2 are all special cases of our
// radial_poly, since the forward converter (tools/lensfun_to_lensfit.py)
// builds radial_poly by rearranging exactly these formulas. Recovering them
// is the same rearrangement run backwards: it either matches one of the two
// sparsity patterns exactly, or the entry did not come from (or does not
// reduce to) a Lensfun model and cannot be expressed.
static gboolean _radial_poly_to_lf(const dt_lens_warp_t *w,
                                   dt_lf_dist_entry_t *e,
                                   const char **why)
{
  static const float TOL = 1e-4f;

  if(fabsf(w->p[0] - 1.0f) > TOL)
  {
    *why = "ellipticity is not 1 -- lensfun's distortion models have no "
           "anisotropic term";
    return FALSE;
  }

  const float c1 = w->p[1], c2 = w->p[2], c3 = w->p[3],
              c4 = w->p[4], c5 = w->p[5];

  if(fabsf(c5) < TOL)
  {
    // ptlens: r_d = r ( a r^3 + b r^2 + c r + (1 - a - b - c) ), which
    // forces c1; poly3 is the special case a = 0.
    const float a = c4, b = c3, c = c2;
    if(fabsf(c1 - (1.0f - a - b - c)) > TOL)
    {
      *why = "coefficients do not satisfy ptlens's identity constraint -- "
             "likely hand-edited since import";
      return FALSE;
    }
    e->model = DT_LF_DIST_PTLENS;
    e->terms[0] = a;
    e->terms[1] = b;
    e->terms[2] = c;
    return TRUE;
  }

  if(fabsf(c2) < TOL && fabsf(c4) < TOL && fabsf(c1 - 1.0f) < TOL)
  {
    e->model = DT_LF_DIST_POLY5;
    e->terms[0] = c3;
    e->terms[1] = c5;
    return TRUE;
  }

  *why = "has an r^5 term but does not match poly5's sparsity pattern -- no "
         "lensfun model can hold both a c2/c4 term and a c5 term";
  return FALSE;
}

gboolean dt_lens_profile_export_lensfun(const struct dt_lens_profile_t *p,
                                        const char *path,
                                        GPtrArray *unrepresentable,
                                        GError **error)
{
  if(!p || !path) return FALSE;

  dt_lf_profile_t out;
  dt_lf_profile_init(&out);

  g_strlcpy(out.maker, p->maker, sizeof(out.maker));
  g_strlcpy(out.model, p->model, sizeof(out.model));
  g_strlcpy(out.mount, p->mount, sizeof(out.mount));
  out.crop_factor = p->crop_factor > 0.0f ? p->crop_factor : 1.0f;

  guint written = 0;
  float min_focal = G_MAXFLOAT, max_focal = 0.0f;
  float min_aperture = G_MAXFLOAT, max_aperture = 0.0f;

  for(guint i = 0; p->warps && i < p->warps->len; i++)
  {
    const dt_lens_warp_t *w = &g_array_index(p->warps, dt_lens_warp_t, i);

    if(w->have_geometry)
    {
      if(w->focal <= 0.0f)
      {
        _note(unrepresentable,
              "distortion: no focal length recorded -- lensfun has no "
              "\"any focal length\" entry, so it cannot be included");
      }
      else if(w->kind != DT_LENS_WARP_RADIAL_POLY)
      {
        _note(unrepresentable,
              "distortion at %.1fmm: %s has no lensfun equivalent",
              w->focal, dt_lens_warp_kind_name(w->kind));
      }
      else
      {
        dt_lf_dist_entry_t e = { DT_LF_DIST_NONE, w->focal, { 0 } };
        const char *why = NULL;
        if(_radial_poly_to_lf(w, &e, &why))
        {
          g_array_append_val(out.dist, e);
          written++;
          min_focal = MIN(min_focal, w->focal);
          max_focal = MAX(max_focal, w->focal);
        }
        else
        {
          _note(unrepresentable, "distortion at %.1fmm: %s", w->focal, why);
        }
      }
    }

    const gboolean tca_identity =
      fabsf(w->tca_r[0]) < 1e-6f && fabsf(w->tca_r[1]) < 1e-6f
      && fabsf(w->tca_r[2] - 1.0f) < 1e-6f
      && fabsf(w->tca_b[0]) < 1e-6f && fabsf(w->tca_b[1]) < 1e-6f
      && fabsf(w->tca_b[2] - 1.0f) < 1e-6f;

    if(!tca_identity)
    {
      if(w->focal <= 0.0f)
      {
        _note(unrepresentable,
              "tca: no focal length recorded -- lensfun has no \"any focal "
              "length\" entry, so it cannot be included");
      }
      else
      {
        // r_channel = r (t0 r^2 + t1 r + t2) is exactly lensfun's poly3/acm
        // form, and linear is the special case t0 = t1 = 0 -- always exact.
        dt_lf_tca_entry_t e = { DT_LF_TCA_NONE, w->focal, { 0 } };
        if(fabsf(w->tca_r[0]) < 1e-6f && fabsf(w->tca_r[1]) < 1e-6f
           && fabsf(w->tca_b[0]) < 1e-6f && fabsf(w->tca_b[1]) < 1e-6f)
        {
          e.model = DT_LF_TCA_LINEAR;
          e.terms[0] = w->tca_r[2];
          e.terms[1] = w->tca_b[2];
        }
        else
        {
          e.model = DT_LF_TCA_POLY3;
          e.terms[0] = w->tca_r[0];
          e.terms[1] = w->tca_r[1];
          e.terms[2] = w->tca_r[2];
          e.terms[3] = w->tca_b[0];
          e.terms[4] = w->tca_b[1];
          e.terms[5] = w->tca_b[2];
        }
        g_array_append_val(out.tca, e);
        written++;
        min_focal = MIN(min_focal, w->focal);
        max_focal = MAX(max_focal, w->focal);
      }
    }

    if(w->have_vig)
    {
      if(w->focal <= 0.0f)
      {
        _note(unrepresentable,
              "vignetting: no focal length recorded -- lensfun has no "
              "\"any focal length\" entry, so it cannot be included");
      }
      else if(fabsf(w->vig_ex - 1.0f) > 1e-4f)
      {
        _note(unrepresentable,
              "vignetting at %.1fmm: ellipticity is not 1 -- lensfun's "
              "vignetting model has no anisotropic term", w->focal);
      }
      else if(w->vig_convention != DT_LENS_VIG_TRANSMISSION)
      {
        _note(unrepresentable,
              "vignetting at %.1fmm: stored as correction gain, and 1/gain "
              "is not a polynomial -- cannot be converted to lensfun's "
              "transmission model without refitting", w->focal);
      }
      else
      {
        const float ap = w->vig_aperture > 0.0f ? w->vig_aperture
                                                  : w->aperture;
        dt_lf_vign_entry_t e = { w->focal, ap, w->focus_distance,
                                 { w->vig_k[0], w->vig_k[1], w->vig_k[2] } };
        g_array_append_val(out.vign, e);
        written++;
        min_focal = MIN(min_focal, w->focal);
        max_focal = MAX(max_focal, w->focal);
        if(ap > 0.0f)
        {
          min_aperture = MIN(min_aperture, ap);
          max_aperture = MAX(max_aperture, ap);
        }
      }
    }
  }

  if(!written)
  {
    _note(unrepresentable,
          "nothing in this profile has a lensfun equivalent -- no file "
          "written");
    dt_lf_profile_cleanup(&out);
    return FALSE;
  }

  out.min_focal = min_focal;
  out.max_focal = max_focal;
  if(min_aperture <= max_aperture)
  {
    out.min_aperture = min_aperture;
    out.max_aperture = max_aperture;
  }

  const gboolean ok = dt_lf_write_xml(&out, path, error);
  dt_lf_profile_cleanup(&out);
  return ok;
}

const char *dt_lf_dist_model_name(const dt_lf_dist_model_t m)
{
  switch(m)
  {
    case DT_LF_DIST_POLY3:  return "poly3";
    case DT_LF_DIST_POLY5:  return "poly5";
    case DT_LF_DIST_PTLENS: return "ptlens";
    default:                return "none";
  }
}

const char *dt_lf_tca_model_name(const dt_lf_tca_model_t m)
{
  switch(m)
  {
    case DT_LF_TCA_LINEAR: return "linear";
    case DT_LF_TCA_POLY3:  return "poly3";
    default:               return "none";
  }
}

const char *const *dt_lf_dist_term_names(const dt_lf_dist_model_t m)
{
  static const char *const poly3[] = { "k1", NULL };
  static const char *const poly5[] = { "k1", "k2", NULL };
  static const char *const ptlens[] = { "a", "b", "c", NULL };
  static const char *const none[] = { NULL };

  switch(m)
  {
    case DT_LF_DIST_POLY3:  return poly3;
    case DT_LF_DIST_POLY5:  return poly5;
    case DT_LF_DIST_PTLENS: return ptlens;
    default:                return none;
  }
}

const char *const *dt_lf_tca_term_names(const dt_lf_tca_model_t m)
{
  static const char *const linear[] = { "kr", "kb", NULL };
  static const char *const poly3[] = { "br", "cr", "vr", "bb", "cb", "vb", NULL };
  static const char *const none[] = { NULL };

  switch(m)
  {
    case DT_LF_TCA_LINEAR: return linear;
    case DT_LF_TCA_POLY3:  return poly3;
    default:               return none;
  }
}

int dt_lf_dist_term_count(const dt_lf_dist_model_t m)
{
  const char *const *n = dt_lf_dist_term_names(m);
  int c = 0;
  while(n[c]) c++;
  return c;
}

int dt_lf_tca_term_count(const dt_lf_tca_model_t m)
{
  const char *const *n = dt_lf_tca_term_names(m);
  int c = 0;
  while(n[c]) c++;
  return c;
}

/* Ordered to match the lfLensType values, so a combo box position is the
   stored type without a translation table. */
const char *const *dt_lf_type_names(void)
{
  static const char *const names[] =
  {
    "unknown",        // LF_UNKNOWN
    "rectilinear",    // LF_RECTILINEAR
    "fisheye",        // LF_FISHEYE
    "panoramic",      // LF_PANORAMIC
    "equirectangular",// LF_EQUIRECTANGULAR
    "orthographic",   // LF_FISHEYE_ORTHOGRAPHIC
    "stereographic",  // LF_FISHEYE_STEREOGRAPHIC
    "equisolid",      // LF_FISHEYE_EQUISOLID
    "fisheye_thoby",  // LF_FISHEYE_THOBY
    NULL
  };
  return names;
}

const char *dt_lf_type_name(const int type)
{
  const char *const *names = dt_lf_type_names();
  int n = 0;
  while(names[n]) n++;
  return (type >= 0 && type < n) ? names[type] : names[0];
}


// clang-format off
// modelines: These editor modelines have been set for all relevant files by tools/update_modelines.py
// vim: shiftwidth=2 expandtab tabstop=2 cindent
// kate: tab-indents: off; indent-width 2; replace-tabs on; indent-mode cstyle; remove-trailing-spaces modified;
// clang-format on
