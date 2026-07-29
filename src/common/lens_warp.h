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

/* Two dimensional lens warps.
 *
 * Everything here maps *observed* image coordinates to *corrected* ones --
 * the direction a calibration measures. The pixel pipe needs the opposite
 * direction, which dt_lens_warp_invert supplies iteratively.
 *
 * Coordinates are normalized: the frame centre is the origin and one unit
 * is half the frame diagonal, so u and v run to about +-0.8 and +-0.6 on a
 * 3:2 frame. Normalizing by the diagonal rather than per axis keeps the
 * coordinate system isotropic, which matters because an anamorphic lens is
 * anisotropic and we want that to show up in the model rather than be
 * hidden by the coordinate convention.
 *
 * Three backends, in increasing order of freedom:
 *
 *   POLY         a bivariate polynomial per axis. Makes no assumption of
 *                radial symmetry, so it can express the asymmetry of a
 *                real anamorphic lens, but needs many coefficients and
 *                extrapolates badly beyond the measured area.
 *   ANAM_RADIAL  radial distortion measured in an elliptically scaled
 *                space. Few parameters, extrapolates sanely, and is the
 *                right shape for a spherical lens behind an anamorphic
 *                front element.
 *   TPS          a thin plate spline correction layer on top of a POLY
 *                base. Interpolates the measured residual exactly, so it
 *                captures whatever the parametric models missed -- at the
 *                cost of following measurement noise, and of having
 *                nothing to say outside the measured area.
 *
 * The anamorphic squeeze is deliberately *not* part of any backend. It is
 * a separate final stage, because it is identified from different evidence
 * than the distortion: distortion comes from lines failing to be straight,
 * squeeze from cells failing to be square. Keeping them apart means a
 * squeeze estimate that turns out to be wrong can be corrected without
 * re-solving the distortion.
 */

#define DT_LENS_WARP_MAX_ORDER 5
#define DT_LENS_WARP_MAX_PARAMS 64

typedef enum dt_lens_warp_kind_t
{
  DT_LENS_WARP_POLY = 0,
  DT_LENS_WARP_ANAM_RADIAL = 1,
  DT_LENS_WARP_TPS = 2,

  /* A general radial polynomial, in r itself rather than in r^2:
   *
   *     r_observed = c1 r + c2 r^2 + c3 r^3 + c4 r^4 + c5 r^5
   *
   * with r measured in the elliptically scaled space, as ANAM_RADIAL. p[0] is
   * the ellipticity, p[1..5] are c1..c5.
   *
   * Two things about it are deliberate and both matter.
   *
   * Odd powers: every Lensfun distortion model is of this shape, and four
   * fifths of that database is ptlens, which has r^2, r^3 and r^4 terms. A
   * model in r^2 alone -- ours until now -- cannot express them at all, so
   * most of the world's lens calibration was out of reach.
   *
   *     poly3   c1 = 1-k1,      c3 = k1
   *     poly5   c1 = 1,         c3 = k1,  c5 = k2
   *     ptlens  c1 = 1-a-b-c,   c2 = c,   c3 = b,  c4 = a
   *
   * Direction: this maps *corrected to observed*, the opposite of every other
   * kind here, because that is the direction Lensfun's coefficients are
   * defined in. Converting them to our direction would mean inverting a
   * polynomial into something that is not one, so the direction is carried
   * instead of the error. dt_lens_warp_apply therefore solves it numerically
   * for this kind, and dt_lens_warp_invert evaluates it directly -- the two
   * swap roles.
   */
  DT_LENS_WARP_RADIAL_POLY = 3
} dt_lens_warp_kind_t;

/* Which way round a stored vignetting polynomial reads.
 *
 * Not a detail. Lensfun's header says the polynomial multiplies the source to
 * give the corrected pixel, but every coefficient set in the database has a
 * negative first term, which as a gain would darken the corners -- the
 * direction a lens darkens them, not the direction a correction undoes them.
 * The values settle it where the wording does not: the polynomial is
 * transmission and the correction divides. Since 1/t is not a polynomial, the
 * coefficients cannot be converted, only labelled.
 */
typedef enum dt_lens_vig_convention_t
{
  DT_LENS_VIG_GAIN = 0,        // our own: multiply by the polynomial
  DT_LENS_VIG_TRANSMISSION = 1 // Lensfun "pa": divide by the polynomial
} dt_lens_vig_convention_t;

typedef struct dt_lens_warp_t
{
  dt_lens_warp_kind_t kind;
  int order;            // polynomial order, 2..DT_LENS_WARP_MAX_ORDER

  float focal;          // focal length this was measured at, 0 if unknown

  /* Optical centre, as an offset from the frame centre in the normalized
     units above. The distortion is expanded about this point, not about
     the middle of the sensor. */
  float cx, cy;

  /* Horizontal expansion applied after correction. 2.0 desqueezes a 2x
     anamorphic. 1.0 for a spherical lens. */
  float squeeze;

  /* Lateral chromatic aberration: a per channel radial scaling, as a
     polynomial rather than a constant.
   *
   *     r_channel = r * (t[0] r^2 + t[1] r + t[2])
   *
   * Lensfun's linear model is the case t = {0, 0, k}; its poly3 model uses all
   * three, and a scalar cannot hold it. Green is the reference and is never
   * scaled. */
  float tca_r[3], tca_b[3];

  /* Which stop and which focus distance this was measured at. Vignetting
     depends on aperture strongly and on focus distance slightly, and Lensfun
     indexes by both -- 23757 of its measurements collapse onto 1292 lenses
     without these two axes. Zero means not recorded, which is not the same as
     a measurement of zero. */
  float aperture;
  float focus_distance;

  int nparams;
  float p[DT_LENS_WARP_MAX_PARAMS];

  /* Thin plate spline layer, used when kind is TPS. `tps_src` holds
     tps_count control points as consecutive u,v pairs in observed
     normalized coordinates; the weight arrays hold tps_count spline
     weights followed by three affine terms. */
  int tps_count;
  float *tps_src;
  float *tps_wx;
  float *tps_wy;

  /* Vignetting: the fall in brightness away from the optical centre.
   *
   * Stored as the correction *gain*: g(r) = 1 + k0 r^2 + k1 r^4 + k2 r^6,
   * greater than one where the frame is darker than its centre, and applied
   * by multiplying. Equivalently it is the ratio of centre brightness to
   * brightness at r, which is what the fit measures directly.
   *
   * Lensfun's "pa" model uses the opposite convention -- its coefficients
   * describe transmission, so they come out negative and the correction
   * divides. The polynomial has the same shape, but exporting there is a
   * conversion and not a copy.
   *
   * The radius is measured in the same elliptically scaled space as
   * ANAM_RADIAL, using `vig_ex`. Mechanical vignetting is imaged through the
   * same anamorphic front element as everything else, so its iso-brightness
   * contours are ellipses, not circles -- forcing them circular puts the
   * error back as a cross shaped residual.
   *
   * Vignetting depends on aperture far more strongly than on focal length,
   * which is why the f-number it was measured at is recorded with it. There
   * is no interpolation across aperture: a profile says what it measured.
   */
  /* Whether this entry actually describes geometry.
   *
   * A warp can carry vignetting and nothing else, because vignetting is
   * measured on axes distortion does not have -- one focal length can hold a
   * dozen vignetting measurements at different apertures. Those entries must
   * stay out of the focal-length interpolation: blending a real distortion
   * against an identity one would quietly halve the correction, and it would
   * look like a mild lens rather than a bug. */
  gboolean have_geometry;

  gboolean have_vig;
  float vig_k[3];
  float vig_ex;       // ellipticity of the falloff, 1.0 for circular
  float vig_aperture; // f-number it was measured at, 0 if unknown
  dt_lens_vig_convention_t vig_convention;

  /* Lensfun normalises the radius for distortion and for vignetting
     differently -- it carries an explicit AspectRatioCorrection between the
     two -- so the normalisation belongs on the measurement, not on the
     profile. TRUE means r was normalised by half the frame diagonal, as
     everything else here is. */
  gboolean vig_r_half_diagonal;

  /* How much the correction changes the useful frame, as two scale factors
     on the corrected frame.
   *
   * `overscan` >= 1: how much larger the output must be to keep every
   * recorded pixel. Correcting barrel distortion pushes the edges outward, so
   * rendering into the original frame size throws away real image -- which
   * matters to anyone who intends to reframe or stabilise afterwards, and is
   * why a compositor asks for overscan by name.
   *
   * `underscan` <= 1: how much smaller the output must be to have no empty
   * edge at all. Correcting pincushion pulls the edges inward and leaves
   * nothing behind them; this is the largest centred crop that is entirely
   * real image.
   *
   * They are properties of the warp and so could be recomputed on demand,
   * but they are recorded because they are what a downstream tool needs to
   * be *told*: an STmap on its own does not say how much of the plate it
   * expects, and getting it wrong is a silent crop rather than an error.
   */
  float overscan;
  float underscan;
} dt_lens_warp_t;

/* A Lensfit profile: one lens, and every measurement made of it.
 *
 * Holds several warps, one per measured focal length, so a zoom is described
 * by more than its widest setting.
 *
 * The format is our own rather than Lensfun's because Lensfun cannot hold
 * what we measure. Every model it has is a function of a single radius, with
 * no anisotropic term and no way for a profile to change the output aspect
 * ratio -- so an anamorphic squeeze is inexpressible there, not merely
 * awkward. The relationship runs one way: everything Lensfun can say fits in
 * here (its radial polynomials are this model with ellipticity 1), while the
 * reverse loses the squeeze entirely.
 */
/* Where a profile's numbers came from.
 *
 * This is not bookkeeping. A profile that was measured, one converted from
 * Lensfun and one typed in by hand are three different kinds of claim, and
 * the aggregator weights them differently -- it has to, because agreement
 * between a measurement and a copy of somebody else's measurement is not
 * corroboration. Saving every profile as "measured" regardless, which is
 * what this replaced, made all three indistinguishable at exactly the point
 * where the difference matters.
 *
 * It is also a licensing fact: Lensfun's data is CC-BY-SA-3.0, and a
 * converted profile has to carry that forward rather than present itself as
 * original work.
 */
typedef enum dt_lens_source_t
{
  DT_LENS_SOURCE_MEASURED = 0,   // fitted here, from a chart shot
  DT_LENS_SOURCE_MANUFACTURER,   // published by whoever made the lens
  DT_LENS_SOURCE_LENSFUN,        // converted out of the Lensfun database
  DT_LENS_SOURCE_AGGREGATED,     // combined from several submissions
  DT_LENS_SOURCE_REVERSE_ENG,    // recovered from a camera's own metadata
  DT_LENS_SOURCE_EDITED          // a measurement since altered by hand
} dt_lens_source_t;

const char *dt_lens_source_name(const dt_lens_source_t s);
dt_lens_source_t dt_lens_source_from_name(const char *name);

typedef struct dt_lens_profile_t
{
  char name[128];
  char maker[128];
  char model[128];
  char mount[128];

  int width, height;    // pixel size of the frame it was calibrated on
  float crop_factor;

  /* The lens's own advertised range, not the range of what has been
     measured so far -- a 24-70mm zoom is a 24-70mm zoom even with only
     35mm actually calibrated. 0 means unknown in each case, the same
     convention as everywhere else here.
   *
   * Two ways to arrive at a value: typed in directly, or left alone and
   * derived from the focal/aperture/focus_distance actually present across
   * `warps` -- which one just happened is not distinguished here, because
   * by the time a profile is saved there is nothing left to distinguish:
   * both are "the range", and a reader has no reason to prefer one
   * provenance of it over the other. The distinction matters only in the
   * editor, where it decides whether typing a value should stick against
   * future measurements or keep tracking them -- see the lens panel. */
  float focal_min, focal_max;
  float aperture_min, aperture_max;
  float distance_min, distance_max;

  dt_lens_source_t source;
  char license[64];     // set when the source imposes one, else empty
  char parent[160];     // the profile this was derived from, if any

  GArray *warps;        // dt_lens_warp_t, kept sorted by focal
} dt_lens_profile_t;


void dt_lens_warp_init(dt_lens_warp_t *w,
                       const dt_lens_warp_kind_t kind,
                       const int order);
void dt_lens_warp_cleanup(dt_lens_warp_t *w);
// deep copy, including the spline layer; `dst` must not be live
gboolean dt_lens_warp_copy(dt_lens_warp_t *dst, const dt_lens_warp_t *src);

// number of polynomial coefficients per axis at a given order
int dt_lens_warp_poly_terms(const int order);
// how many parameters a kind and order needs
int dt_lens_warp_param_count(const dt_lens_warp_kind_t kind, const int order);
const char *dt_lens_warp_kind_name(const dt_lens_warp_kind_t kind);

/* observed -> corrected. Always succeeds; a warp with no parameters set is
   the identity. */
void dt_lens_warp_apply(const dt_lens_warp_t *w,
                        const float u, const float v,
                        float *ou, float *ov);

/* corrected -> observed, by fixed point iteration. Returns FALSE if it did
   not converge, in which case the output holds the best guess reached --
   usable for a preview, not for a measurement. */
gboolean dt_lens_warp_invert(const dt_lens_warp_t *w,
                             const float u, const float v,
                             float *ou, float *ov);

/* Radial scale for the red or blue channel at radius `r`, undoing lateral
   chromatic aberration. Green is the reference and returns 1. */
float dt_lens_tca_scale(const dt_lens_warp_t *w,
                        const int channel, // 0 red, 1 green, 2 blue
                        const float r);

/* Work out `overscan` and `underscan` for a warp on a frame of this size,
   and store them in it. Cheap -- it walks the frame border rather than the
   whole image, which is enough because the warp is smooth and monotonic
   outward. */
void dt_lens_warp_measure_scan(dt_lens_warp_t *w,
                               const int width,
                               const int height);

/* Brightness gain to multiply a pixel by, undoing the vignetting. Returns
   1.0 for a warp with no vignetting measured. Clamped, since the sixth
   order term can go negative in the far corners of an over-fitted model and
   a negative gain would invert the image rather than brighten it. */
float dt_lens_vignette_gain(const dt_lens_warp_t *w,
                            const float u, const float v);

/* Install a thin plate spline layer fitted to `count` displacements
   measured at `src` (observed normalized coordinates). `dx`/`dy` are the
   corrections to add on top of whatever the base model already produces.
   `smooth` >= 0 relaxes the exact interpolation towards a least squares
   fit, which is what keeps measurement noise from being reproduced
   faithfully. */
gboolean dt_lens_warp_fit_tps(dt_lens_warp_t *w,
                              const float *const src,
                              const float *const dx,
                              const float *const dy,
                              const int count,
                              const float smooth);

/* ---------------------------------------------------------- profiles */

void dt_lens_profile_init(dt_lens_profile_t *p);
void dt_lens_profile_cleanup(dt_lens_profile_t *p);

/* Add a warp, replacing any existing one at the same focal length. Takes
   ownership of nothing; `w` is deep copied. */
gboolean dt_lens_profile_add(dt_lens_profile_t *p, const dt_lens_warp_t *w);

/* The warp to use at `focal`, interpolating between the two measurements
   that bracket it. Returns FALSE if the profile is empty. `out` must not
   be live and becomes the caller's to clean up.
 *
 * Interpolation is linear in log focal length, which fits how distortion
 * actually varies across a zoom range far better than linear in focal:
 * the change from 16 to 24mm is a much bigger optical step than 200 to
 * 208mm. Outside the measured range the nearest measurement is held
 * rather than extrapolated. */
gboolean dt_lens_profile_eval(const dt_lens_profile_t *p,
                              const float focal,
                              dt_lens_warp_t *out);

/* As above, but also states the aperture and focus distance the frame was
   shot at. Geometry is chosen by focal length alone -- distortion does not
   depend on aperture -- while vignetting is chosen by aperture first, because
   it depends on aperture far more than on anything else. Pass 0 for either to
   mean "unknown", which is not the same as a value of zero. */
gboolean dt_lens_profile_eval_at(const dt_lens_profile_t *p,
                                 const float focal,
                                 const float aperture,
                                 const float distance,
                                 dt_lens_warp_t *out);

gboolean dt_lens_profile_save(const dt_lens_profile_t *p,
                              const char *path,
                              GError **error);
gboolean dt_lens_profile_load(dt_lens_profile_t *p, const char *path);

// directory holding the user's own calibration profiles (writable), created
// if needed. This is always the *save* target -- a shipped database is
// read-only, so a saved or edited profile can only ever land here.
gchar *dt_lens_profile_dir(void);
// NULL terminated list of profile names (file basenames without .json),
// merged from the user directory and the shipped datadir database, deduped
// by name.
gchar **dt_lens_profile_list(void);
// same coverage as dt_lens_profile_list(), but also returns the maker and
// model for each entry (parallel NULL-terminated arrays, same length and
// order as *out_names) so a caller can build a vendor-then-model cascading
// picker instead of one flat list -- the shipped database alone is >1000
// entries, too long to browse flat. A profile with no known maker (rare --
// only user profiles saved before that field existed) gets maker="", model
// set to the profile name itself, so it still shows up unfiltered under an
// "(other)" bucket. Frees with g_strfreev on all three.
void dt_lens_profile_list_full(gchar ***out_names,
                               gchar ***out_makers,
                               gchar ***out_models);
// full path for a profile name, in the user directory, for *writing* --
// does not check whether anything exists there yet.
gchar *dt_lens_profile_path(const char *name);
// full path for a profile name, for *reading* -- searches the user
// directory first, then the shipped datadir database, so a user profile of
// the same name always wins. NULL if the name exists in neither.
gchar *dt_lens_profile_find(const char *name);

/* Automatic profile matching from a camera/lens maker+model, e.g. from
 * EXIF. Three independent strategies, meant to be tried in this order by
 * the caller -- see reload_defaults() in iop/lens.cc:
 *
 *  1. dt_lens_profile_match_user() -- search the user's own saved
 *     profiles. Cheap (a handful of files, not 1290) and deliberately
 *     tried first: a profile someone measured themselves must always win
 *     over a converted one for the same lens.
 *  2. dt_lens_profile_manifest_lookup() -- exact lookup by the name a
 *     profile would be given, generated from `maker`/`model` the same way
 *     the converter names its output. Meant for Lensfun's own maker/model
 *     strings once Lensfun itself has already identified the lens: this
 *     reproduces the profile's name exactly rather than guessing at it.
 *  3. dt_lens_profile_manifest_match() -- fuzzy fallback for raw EXIF
 *     text, when there is no Lensfun match to hand a clean name over. Both
 *     maker and model must partially match, deliberately conservative: a
 *     wrong-but-plausible auto-match that looks corrected is worse than no
 *     automatic match at all.
 *
 * All three return FALSE rather than a weak guess when nothing clears
 * their bar -- there is no "closest match" here on purpose.
 */
gboolean dt_lens_profile_match_user(const char *maker, const char *model,
                                    char *out_name, size_t out_size);
gboolean dt_lens_profile_manifest_lookup(const char *maker, const char *model,
                                         char *out_name, size_t out_size);
gboolean dt_lens_profile_manifest_match(const char *maker, const char *model,
                                        char *out_name, size_t out_size);

G_END_DECLS

// clang-format off
// modelines: These editor modelines have been set for all relevant files by tools/update_modelines.py
// vim: shiftwidth=2 expandtab tabstop=2 cindent
// kate: tab-indents: off; indent-width 2; replace-tabs on; indent-mode cstyle; remove-trailing-spaces modified;
// clang-format on
