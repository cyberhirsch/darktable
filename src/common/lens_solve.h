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

#include "common/lens_warp.h"

#include <glib.h>

G_BEGIN_DECLS

/* Fitting a warp to a photographed grid chart.
 *
 * The objective is the plumb-line one: whatever else a lens does, it must
 * not bend straight lines, so the fit minimises how far each traced chart
 * line deviates from the straight line that best fits it *after*
 * correction. The attraction of this is that it needs no knowledge of
 * where the camera was. Any projective transform maps straight lines to
 * straight lines, so camera pose simply cannot affect the objective, and
 * the chart does not have to be square on, measured, or even flat-fronted
 * -- only flat.
 *
 * What plumb lines cannot see is anything that keeps lines straight. That
 * is exactly the anamorphic squeeze, which is a pure scale along one axis,
 * so it is recovered in a second stage from the shape of the cells rather
 * than the straightness of the lines -- see dt_lens_solve_result_t::squeeze
 * and the caveat on `affine_rms_px`.
 *
 * There is one genuine ambiguity worth knowing about: a projective
 * transform preserves straightness, and a second order polynomial can
 * imitate one over a limited area. So a tilted chart lets a little of the
 * pose leak into the distortion model. Regularisation damps it; shooting
 * the chart square on removes it.
 */

typedef struct dt_lens_solve_point_t
{
  float x, y;   // observed position, image pixels
  int col, row; // lattice index
} dt_lens_solve_point_t;

typedef struct dt_lens_solve_input_t
{
  const dt_lens_solve_point_t *points;
  int count;

  int width, height; // frame size the positions refer to

  /* Width over height of one chart cell. 1.0 for the square cells of the
     usual target. Only used for the squeeze stage. */
  float cell_aspect;
} dt_lens_solve_input_t;

typedef struct dt_lens_solve_options_t
{
  dt_lens_warp_kind_t kind;
  int order;

  gboolean solve_centre;  // let the optical centre move off the frame centre

  /* The chart was photographed square on, so its corrected image should
     be an affine image of the chart -- not merely a projective one.
   *
   * This is worth far more than the squeeze it also enables. Plumb lines
   * are blind to any transform that preserves straightness, and that
   * includes the projective ones; a free enough model will happily absorb
   * a keystone that was never in the lens, fit the lines perfectly, and
   * leave the corrected image visibly skewed. Asserting that the chart
   * was square on removes that freedom, at the cost of being wrong if it
   * was not -- so it is the user's statement to make, not ours.
   */
  gboolean chart_frontal;

  /* Tikhonov weight on the model coefficients, applied more strongly to
     the higher orders. Buys stability at the edges, where there is least
     data and the polynomial is keenest to misbehave, at the cost of a
     slightly worse fit in the middle. */
  float regularization;

  /* The squeeze the lens is known to have -- 1.6 for a 1.6x anamorphic,
     1.0 for a spherical lens. Zero means measure it from the chart instead.
   *
   * Declaring it is nearly always better. The squeeze is invisible to plumb
   * lines, so it cannot come out of the straightness fit at all; it has to
   * be inferred from cell shape, which needs the chart to have been square
   * on and confuses any foreshortening with squeeze. A number printed on the
   * lens barrel is better evidence than either. The measurement is still
   * made and reported, as a check on the declaration.
   */
  float known_squeeze;

  // relaxation of the spline layer towards least squares, TPS only
  float tps_smooth;

  int max_iter;
} dt_lens_solve_options_t;

typedef struct dt_lens_solve_result_t
{
  gboolean ok;

  float rms_px;        // straightness residual after the fit
  float rms_before_px; // and before it, so the improvement is visible

  /* Residual of the affine fit used for the squeeze. Large means the chart
     was not square on, which makes `squeeze` unreliable -- the perspective
     foreshortening of a tilted chart looks exactly like a squeeze. */
  float affine_rms_px;

  /* The squeeze the model carries -- the declared one when there is one. */
  float squeeze;

  /* And what the cell shapes actually implied, always measured even when a
     value was declared. Disagreement between the two is worth seeing: it
     means either the wrong ratio was entered or the chart was not square on. */
  float squeeze_measured;

  int iterations;
  int used_points;
  int lines_used; // rows and columns with enough points to constrain a line
} dt_lens_solve_result_t;

void dt_lens_solve_default_options(dt_lens_solve_options_t *opt);

/* Fit `warp` to the observations. `warp` is initialised by the call and
   becomes the caller's to clean up. */
gboolean dt_lens_solve(const dt_lens_solve_input_t *in,
                       const dt_lens_solve_options_t *opt,
                       dt_lens_warp_t *warp,
                       dt_lens_solve_result_t *res);

/* Per point straightness residual for display, in pixels. `out_dx`/`out_dy`
   receive the offset from the point's corrected position to the line it
   should have been on, and `out_mag` its length; each array holds
   `in->count` entries and any may be NULL. */
void dt_lens_solve_residuals(const dt_lens_solve_input_t *in,
                             const dt_lens_warp_t *warp,
                             float *out_dx,
                             float *out_dy,
                             float *out_mag);

G_END_DECLS

// clang-format off
// modelines: These editor modelines have been set for all relevant files by tools/update_modelines.py
// vim: shiftwidth=2 expandtab tabstop=2 cindent
// kate: tab-indents: off; indent-width 2; replace-tabs on; indent-mode cstyle; remove-trailing-spaces modified;
// clang-format on
