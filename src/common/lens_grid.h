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
#include <stdint.h>

/* Detection of a photographed grid chart.
 *
 * The lines of the chart are curved by the very distortion we are trying
 * to measure, and most curved near the frame edges where the measurement
 * matters most. A Hough transform looks for straight lines and would be
 * least reliable exactly there, so instead each line is traced as a
 * polyline: scan across the image, find the line crossings on every scan
 * position to sub-pixel accuracy, then link those into curves.
 *
 * Intersecting traced curve i with traced curve j yields a point that is
 * already indexed (i, j), so the usual correspondence problem -- deciding
 * which detected blob is which lattice node -- does not arise.
 */

// one sample along a traced grid line
typedef struct dt_lens_grid_sample_t
{
  float pos;   // position along the scan axis (the scanline coordinate)
  float coord; // sub-pixel position of the line crossing
} dt_lens_grid_sample_t;

// a single traced grid line
typedef struct dt_lens_grid_curve_t
{
  dt_lens_grid_sample_t *samples;
  int count;
  int allocated;
  gboolean horizontal; // TRUE: a line running left-right, scanned by column
  int ordinal;         // lattice index once assigned, -1 while unknown
  int misses;          // consecutive scan positions with no match, while tracing
} dt_lens_grid_curve_t;

// an intersection of a horizontal and a vertical curve
typedef struct dt_lens_grid_point_t
{
  float x, y; // image coordinates, sub-pixel
  int col, row;
} dt_lens_grid_point_t;

typedef struct dt_lens_grid_t
{
  int width, height;

  // expected chart geometry, in cells. a 28x12 chart has a 29x13 lattice.
  int cells_x, cells_y;

  GList *curves_h; // dt_lens_grid_curve_t*, ordered top to bottom
  GList *curves_v; // dt_lens_grid_curve_t*, ordered left to right

  dt_lens_grid_point_t *points;
  int point_count;
} dt_lens_grid_t;

/* Detect the chart in a single channel image.
 *
 * `lum` is row major, `width` * `height` floats, any range -- it gets
 * locally normalized internally, so vignetting and uneven chart lighting
 * are handled. `dark_lines` selects whether the grid is drawn dark on
 * light (the usual case) or the reverse.
 *
 * Returns FALSE if too few curves were found to be useful. Detecting
 * fewer than the expected number of lines is not itself a failure:
 * partial detection is normal at frame edges and in vignetted corners,
 * and the solver works from whatever subset is found.
 */
gboolean dt_lens_grid_detect(const float *const lum,
                             const int width,
                             const int height,
                             const int cells_x,
                             const int cells_y,
                             const gboolean dark_lines,
                             dt_lens_grid_t *grid);

void dt_lens_grid_cleanup(dt_lens_grid_t *grid);

// clang-format off
// modelines: These editor modelines have been set for all relevant files by tools/update_modelines.py
// vim: shiftwidth=2 expandtab tabstop=2 cindent
// kate: tab-indents: off; indent-width 2; replace-tabs on; indent-mode cstyle; remove-trailing-spaces modified;
// clang-format on
