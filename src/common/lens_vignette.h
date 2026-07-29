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

/* Measuring vignetting from a photograph of an evenly lit surface.
 *
 * The model is the correction gain polynomial stored in dt_lens_warp_t:
 *
 *     g(r) = v(0) / v(r) = 1 + k0 r^2 + k1 r^4 + k2 r^6
 *
 * with r measured in the elliptically scaled normalized space, so that an
 * anamorphic lens's elliptical falloff is described rather than averaged
 * away.
 *
 * The fit is exactly linear, which is worth arranging. Brightness is
 * v(r) = v0 / att(r) with v0 unknown, and fitting that directly is nonlinear
 * in the product of v0 and the coefficients. Fitting the *reciprocal*
 *
 *     1/v(r) = c0 + c1 r^2 + c2 r^4 + c3 r^6
 *
 * is linear in four unknowns, and the model falls straight out of it as
 * v0 = 1/c0 and k_i = c_{i+1}/c0. No iteration, no starting guess, no
 * convergence to worry about.
 *
 * What this cannot do is tell vignetting apart from an unevenly lit subject.
 * That is a property of the photograph, not of the method: a wall that is
 * brighter in the middle is indistinguishable from a lens that is darker at
 * the edges. The caller has to supply a fair subject; `out_rms` says how
 * well the radial model described it, which catches gross failures such as a
 * subject with a shadow across it.
 */

/* Fit from a linear luminance plane.
 *
 * `lum` must be linear -- proportional to scene light, black level already
 * removed. Gamma encoded data produces a confidently wrong answer, since the
 * transfer curve is itself a brightness dependent distortion.
 *
 * A mosaiced buffer is fine. Over a radius bin the CFA pattern contributes a
 * fixed mixture of channels, so the bin statistic stays proportional to
 * illumination, and the overall scale is divided out anyway.
 *
 * `cx`/`cy` are the optical centre as a normalized offset from the frame
 * centre, matching dt_lens_warp_t. `ex_hint` seeds the ellipticity search;
 * pass the warp's squeeze, or 1.0 for a spherical lens.
 *
 * Returns FALSE if there was not enough usable data to fit.
 */
gboolean dt_lens_vignette_fit(const float *const lum,
                              const int width,
                              const int height,
                              const int stride,
                              const float cx,
                              const float cy,
                              const float ex_hint,
                              float out_k[3],
                              float *out_ex,
                              float *out_rms,
                              float *out_corner_stops);

/* Fit from brightness sampled at known-good places.
 *
 * Much the better method when a chart lattice exists. The radial binning
 * above has to guess which pixels are the surface and which are markings on
 * it; given the lattice we simply know -- the centre of a cell is paper by
 * construction, so there is no ink to reject and no percentile to tune.
 *
 * `uv` holds `count` normalized coordinate pairs, in the same convention as
 * dt_lens_warp_t: origin at the frame centre, one unit is half the frame
 * diagonal. `value` holds the linear brightness measured at each.
 *
 * The catch is coverage, and it is reported rather than hidden. A chart
 * rarely reaches the frame corners, so the samples span a smaller radius
 * than the picture does and the corner figure becomes an extrapolation of a
 * sixth order polynomial -- which is exactly where such a polynomial is
 * least trustworthy. `out_max_r` gives the largest radius actually measured;
 * compare it against the corner radius before believing the corner.
 */
gboolean dt_lens_vignette_fit_points(const float *const uv,
                                     const float *const value,
                                     const int count,
                                     const int width,
                                     const int height,
                                     const float cx,
                                     const float cy,
                                     const float ex_hint,
                                     float out_k[3],
                                     float *out_ex,
                                     float *out_rms,
                                     float *out_max_r,
                                     float *out_corner_stops);

G_END_DECLS

// clang-format off
// modelines: These editor modelines have been set for all relevant files by tools/update_modelines.py
// vim: shiftwidth=2 expandtab tabstop=2 cindent
// kate: tab-indents: off; indent-width 2; replace-tabs on; indent-mode cstyle; remove-trailing-spaces modified;
// clang-format on
