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

/* STmap export.
 *
 * An STmap is an image whose red and green channels hold, for every output
 * pixel, the normalized coordinates of the source pixel to fetch. It is
 * the lingua franca for handing a lens correction to a compositor, and --
 * more to the point here -- the only interchange format that can carry an
 * asymmetric warp at all. Lensfun's models are radial by construction, so
 * an anamorphic calibration written as Lensfun XML is an approximation
 * whereas the same calibration written as an STmap is exact to the
 * sampling grid.
 *
 * 32 bit float is not optional: half float has about three decimal digits
 * of mantissa, which across a 6000 pixel frame quantises the lookup to
 * several pixels and shows up as visible stepping.
 */

/* Write the map that corrects `w` at `width` x `height`.
 *
 * `bottom_up` selects the row order the coordinates are expressed in.
 * Nuke, Fusion and most compositors put the origin at the bottom left, so
 * this defaults to what they expect; leave it off for tools that follow
 * the image processing convention of a top left origin. Getting it wrong
 * flips the correction vertically, which is obvious on any real lens and
 * invisible on a symmetric test pattern.
 */
gboolean dt_lens_stmap_write(const dt_lens_warp_t *w,
                             const int width,
                             const int height,
                             const gboolean bottom_up,
                             const char *path,
                             GError **error);

G_END_DECLS

// clang-format off
// modelines: These editor modelines have been set for all relevant files by tools/update_modelines.py
// vim: shiftwidth=2 expandtab tabstop=2 cindent
// kate: tab-indents: off; indent-width 2; replace-tabs on; indent-mode cstyle; remove-trailing-spaces modified;
// clang-format on
