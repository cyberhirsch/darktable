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

#pragma once

#include <glib.h>

G_BEGIN_DECLS

/* Reading, editing and writing Lensfun profiles.
 *
 * Everything Lensfun can store about a lens is represented here, not a
 * convenient subset: a lens carries a list of distortion, TCA and
 * vignetting entries, each measured at a particular focal length -- and
 * for vignetting also at a particular aperture and focus distance. A zoom
 * can easily carry dozens of vignetting rows, so the model is a list, not
 * a fixed set of coefficients.
 *
 * Note the asymmetry that motivates this whole feature: Lensfun's
 * distortion models are all functions of radius alone, so importing gives
 * a radial starting point which an anamorphic calibration then has to
 * depart from.
 *
 * A profile need not come from the database at all: dt_lf_profile_init
 * leaves a valid empty one, entries can be appended to the arrays, and
 * dt_lf_write_xml serialises the result as a Lensfun database fragment. So
 * the same structure covers importing, hand editing and authoring.
 */

#define DT_LF_DIST_MAX_TERMS 3
#define DT_LF_TCA_MAX_TERMS 6

typedef enum dt_lf_dist_model_t
{
  DT_LF_DIST_NONE = 0,
  DT_LF_DIST_POLY3,   // k1
  DT_LF_DIST_POLY5,   // k1, k2
  DT_LF_DIST_PTLENS   // a, b, c
} dt_lf_dist_model_t;

typedef enum dt_lf_tca_model_t
{
  DT_LF_TCA_NONE = 0,
  DT_LF_TCA_LINEAR, // kr, kb
  DT_LF_TCA_POLY3   // br, cr, vr, bb, cb, vb
} dt_lf_tca_model_t;

typedef struct dt_lf_dist_entry_t
{
  dt_lf_dist_model_t model;
  float focal;
  float terms[3]; // meaning depends on model
} dt_lf_dist_entry_t;

typedef struct dt_lf_tca_entry_t
{
  dt_lf_tca_model_t model;
  float focal;
  float terms[6];
} dt_lf_tca_entry_t;

typedef struct dt_lf_vign_entry_t
{
  float focal;
  float aperture;
  float distance;
  float k[3]; // the "pa" model, the only one Lensfun defines
} dt_lf_vign_entry_t;

typedef struct dt_lf_profile_t
{
  char maker[128];
  char model[128];
  char mount[128];

  float min_focal, max_focal;
  float min_aperture, max_aperture;
  float crop_factor;
  float aspect_ratio;
  float centre_x, centre_y; // optical axis offset as Lensfun records it
  int type;                 // lfLensType, kept as int to avoid the include

  GArray *dist;  // dt_lf_dist_entry_t
  GArray *tca;   // dt_lf_tca_entry_t
  GArray *vign;  // dt_lf_vign_entry_t
} dt_lf_profile_t;

/* Search the Lensfun database. `maker` and `model` may be NULL or partial;
   matching is left to Lensfun. Returns a NULL terminated array of newly
   allocated strings in "Maker | Model" form, free with g_strfreev. */
gchar **dt_lf_search_lenses(const char *maker, const char *model);

/* Every lens in the database, in the same "Maker | Model" form. Used to
   feed the search box completion, which needs the whole vocabulary up
   front rather than one query's worth. Free with g_strfreev. */
gchar **dt_lf_all_lens_names(void);

/* Load one lens by its exact Lensfun model name into `out`. Returns FALSE
   if no lens matched. `out` must be cleaned up with dt_lf_profile_cleanup
   either way. */
gboolean dt_lf_load_profile(const char *model, dt_lf_profile_t *out);

/* Lensfun's own per-user override directory (its `HomeDataDir`, typically
   something like `~/.local/share/lensfun` on Linux) -- a file dropped here
   is picked up by Lensfun itself and overrides the system database, which
   is what "usually saved" means for a hand-written or exported Lensfun XML
   fragment. NULL if the database has not been opened yet. Caller does not
   own the returned string. */
const char *dt_lf_home_data_dir(void);

void dt_lf_profile_init(dt_lf_profile_t *p);
void dt_lf_profile_cleanup(dt_lf_profile_t *p);

/* Write `p` as a standalone Lensfun database file. The result is a
   complete <lensdatabase> document holding one <lens>, which is exactly
   what Lensfun expects to find in a file dropped into its data directory,
   and also what an upstream submission looks like. */
gboolean dt_lf_write_xml(const dt_lf_profile_t *p,
                         const char *path,
                         GError **error);

/* Export a lensfit profile as a Lensfun database fragment.
 *
 * Lensfun's models are all a single radial function with no anisotropy and
 * no way to change output aspect ratio, so this is the lossy direction:
 * whatever does not fit -- an anamorphic warp, an anisotropic vignetting
 * fall-off, a gain-convention vignetting curve (1/t is not a polynomial),
 * an entry with no recorded focal length (Lensfun has no "any" sentinel) --
 * is left out rather than approximated, and reported in `unrepresentable`
 * so the caller can say what was lost. Returns FALSE only if nothing at
 * all survived to write; a partial export still returns TRUE with entries
 * in `unrepresentable`.
 *
 * `unrepresentable`, if non-NULL, is appended with newly-allocated
 * (gchar *) strings the caller owns; g_free each and free the array. */
struct dt_lens_profile_t; // defined in common/lens_warp.h
gboolean dt_lens_profile_export_lensfun(const struct dt_lens_profile_t *p,
                                        const char *path,
                                        GPtrArray *unrepresentable,
                                        GError **error);

// human readable names for the enums, for display in the ui
const char *dt_lf_dist_model_name(const dt_lf_dist_model_t m);
const char *dt_lf_tca_model_name(const dt_lf_tca_model_t m);
// names of the coefficients a model uses, NULL terminated
const char *const *dt_lf_dist_term_names(const dt_lf_dist_model_t m);
const char *const *dt_lf_tca_term_names(const dt_lf_tca_model_t m);
// how many of the terms[] slots a model actually uses
int dt_lf_dist_term_count(const dt_lf_dist_model_t m);
int dt_lf_tca_term_count(const dt_lf_tca_model_t m);

/* Projection types, indexed by the lfLensType value stored in
   dt_lf_profile_t::type, so a combo box can map position to value
   directly. NULL terminated. */
const char *const *dt_lf_type_names(void);
const char *dt_lf_type_name(const int type);

G_END_DECLS

// clang-format off
// modelines: These editor modelines have been set for all relevant files by tools/update_modelines.py
// vim: shiftwidth=2 expandtab tabstop=2 cindent
// kate: tab-indents: off; indent-width 2; replace-tabs on; indent-mode cstyle; remove-trailing-spaces modified;
// clang-format on
