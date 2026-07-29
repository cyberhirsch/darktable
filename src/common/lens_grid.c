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

#include "common/lens_grid.h"
#include "common/darktable.h"
#include "common/math.h"

#include <math.h>
#include <stdlib.h>
#include <string.h>

/* How steeply a chart line may run across the scan direction, as a slope.
 *
 * This is what limits how far a line may move from one scan position to the
 * next, and getting the scale of it wrong is expensive. An earlier version
 * allowed a fraction of the *cell spacing* per step -- sixty odd pixels on a
 * real chart -- while the steps themselves are four pixels apart. A line was
 * therefore free to lurch sideways by fifteen times the step, and the tracer
 * duly wandered off one line and onto its neighbours, producing "horizontal"
 * traces whose y moved a hundred and fifty pixels while x moved twenty
 * eight. Tying the allowance to a slope instead makes it a statement about
 * the picture: at 1.0 a line may run at up to 45 degrees, and beyond that it
 * would belong to the other family anyway.
 */
#define DT_GRID_MAX_SLOPE 1.0f

/* How far a line may deviate from where its own recent direction predicts
 * it will be, per pixel travelled.
 *
 * This is the discriminating test, and the slope bound above is only the
 * bootstrap for a curve too young to have a direction yet. Allowing any
 * slope up to the bound at every step lets a trace wander diagonally
 * across the whole frame and still look plausible, which is how junk
 * traces survived a span test meant to exclude them. Distortion bends
 * lines *gradually* -- that is the assumption the whole tracing approach
 * rests on -- so a line that suddenly changes direction is two different
 * lines, and refusing the match costs one scan position rather than a
 * whole curve.
 */
#define DT_GRID_MAX_CURVE 0.15f

// slack for noise in the sub-pixel position, on top of the above
#define DT_GRID_LINK_BASE 2.5f

// samples used to estimate a curve's current direction
#define DT_GRID_SLOPE_WINDOW 6

// scan positions are taken every Nth pixel; the lines are smooth so
// there is no need to visit every column
#define DT_GRID_SCAN_STEP 4

// a detected line must span at least this fraction of the image to be
// kept, which throws away short spurious ridges. Fragments are stitched
// back together before this test, so it can afford to be strict about
// what it sees and still keep an interrupted line.
#define DT_GRID_MIN_SPAN 0.35f

/* ...and it must actually have been *seen* along that span.
 *
 * Span alone is not evidence of a line. A curve is retired after several
 * consecutive misses, so one that gets matched sporadically -- once every
 * few scan positions, by coincidence -- never accumulates enough
 * consecutive misses to die, and ends up claiming a span across the whole
 * frame on the strength of twenty samples. On a real chart those made up
 * the great majority of surviving curves: junk reached 2-8% of the
 * positions it spanned where real lines reached 31-87%, so density
 * separates them cleanly and span does not separate them at all.
 */
#define DT_GRID_MIN_DENSITY 0.20f

// two curves closer together than this fraction of the expected line
// spacing are the same line found twice
#define DT_GRID_DUPLICATE 0.25f

// how many consecutive scan positions a line may go undetected before we
// treat it as finished rather than merely interrupted
#define DT_GRID_MAX_GAP 6


static dt_lens_grid_curve_t *_curve_new(const gboolean horizontal)
{
  dt_lens_grid_curve_t *c = calloc(1, sizeof(dt_lens_grid_curve_t));
  if(!c) return NULL;
  c->allocated = 64;
  c->samples = malloc(sizeof(dt_lens_grid_sample_t) * c->allocated);
  if(!c->samples)
  {
    free(c);
    return NULL;
  }
  c->horizontal = horizontal;
  c->ordinal = -1;
  return c;
}

static void _curve_free(void *p)
{
  dt_lens_grid_curve_t *c = p;
  if(!c) return;
  free(c->samples);
  free(c);
}

static gboolean _curve_append(dt_lens_grid_curve_t *c,
                              const float pos,
                              const float coord)
{
  if(c->count >= c->allocated)
  {
    const int grown = c->allocated * 2;
    dt_lens_grid_sample_t *s =
      realloc(c->samples, sizeof(dt_lens_grid_sample_t) * grown);
    if(!s) return FALSE;
    c->samples = s;
    c->allocated = grown;
  }
  c->samples[c->count].pos = pos;
  c->samples[c->count].coord = coord;
  c->count++;
  return TRUE;
}

/* The direction a curve is currently heading, over its last few samples.
   Short enough to follow a bend, long enough not to be dominated by the
   sub-pixel noise on any single sample. */
static float _curve_slope(const dt_lens_grid_curve_t *cv)
{
  if(cv->count < 2) return 0.0f;

  const int k = MIN(cv->count, DT_GRID_SLOPE_WINDOW);
  const dt_lens_grid_sample_t *a = &cv->samples[cv->count - k];
  const dt_lens_grid_sample_t *b = &cv->samples[cv->count - 1];

  const float dp = b->pos - a->pos;
  if(dp < 1e-3f) return 0.0f;

  return (b->coord - a->coord) / dp;
}

/* Local contrast normalization.
 *
 * Vignetting darkens the corners, which is exactly where the longest
 * baseline measurements live, so a global threshold would drop them. A
 * separable box blur gives the local background; subtracting it leaves
 * the line structure at comparable amplitude everywhere.
 */
static float *_normalize_local(const float *const in,
                              const int width,
                              const int height,
                              const int radius)
{
  float *tmp = dt_alloc_align_float((size_t)width * height);
  float *out = dt_alloc_align_float((size_t)width * height);
  if(!tmp || !out)
  {
    dt_free_align(tmp);
    dt_free_align(out);
    return NULL;
  }

  /* Both passes use a running sum: the window only gains one sample and
     loses one per step, so the cost is independent of the radius. A
     direct sum would be O(n*radius), and the radius here is a whole cell
     -- hundreds of pixels -- which is far too slow to run interactively. */

  // horizontal pass, into tmp
  DT_OMP_FOR()
  for(int y = 0; y < height; y++)
  {
    const float *r = in + (size_t)y * width;
    float *w = tmp + (size_t)y * width;

    float sum = 0.0f;
    int x1 = MIN(width - 1, radius);
    for(int i = 0; i <= x1; i++) sum += r[i];
    int x0 = 0;

    for(int x = 0; x < width; x++)
    {
      w[x] = sum / (float)(x1 - x0 + 1);

      // advance the window for the next x
      const int nx0 = MAX(0, x + 1 - radius);
      const int nx1 = MIN(width - 1, x + 1 + radius);
      while(x1 < nx1) sum += r[++x1];
      while(x0 < nx0) sum -= r[x0++];
    }
  }

  // vertical pass, subtracting the background as we go
  DT_OMP_FOR()
  for(int x = 0; x < width; x++)
  {
    float sum = 0.0f;
    int y1 = MIN(height - 1, radius);
    for(int j = 0; j <= y1; j++) sum += tmp[(size_t)j * width + x];
    int y0 = 0;

    for(int y = 0; y < height; y++)
    {
      const float bg = sum / (float)(y1 - y0 + 1);
      const size_t idx = (size_t)y * width + x;
      out[idx] = in[idx] - bg;

      const int ny0 = MAX(0, y + 1 - radius);
      const int ny1 = MIN(height - 1, y + 1 + radius);
      while(y1 < ny1) sum += tmp[(size_t)(++y1) * width + x];
      while(y0 < ny0) { sum -= tmp[(size_t)y0 * width + x]; y0++; }
    }
  }

  dt_free_align(tmp);
  return out;
}

/* Find line crossings along one scanline.
 *
 * `stride` walks the scan direction, `n` is its length. Lines appear as
 * local minima of the normalized signal when they are dark on light. The
 * extremum is refined to sub-pixel by fitting a parabola through the
 * three samples around it, which is what makes the eventual calibration
 * worth anything -- integer peak positions would cap accuracy at about
 * half a pixel.
 */
static int _find_crossings(const float *const sig,
                           const int n,
                           const int stride,
                           const float threshold,
                           const int margin,
                           float *out,
                           float *strength,
                           const int out_max)
{
  int found = 0;

  /* One crossing per connected below-threshold run.
   *
   * A printed chart line is several pixels wide and soft at the edges, so
   * after background removal its profile is a broad trough rather than a
   * spike, and noise on the floor of that trough splits it into two or
   * three neighbouring local minima. Taking the deepest minimum of each
   * *run* -- each stretch where the signal stays below the threshold --
   * merges those duplicates without needing to know how far apart two
   * lines ought to be.
   *
   * That last part matters. The obvious suppression rule, "minima closer
   * together than some fraction of a cell are the same line", quietly
   * assumes the lines are evenly spaced. Perspective breaks that: on a
   * chart photographed even slightly off axis the far rows crowd together,
   * and a spacing-based rule starts deleting real lines exactly where the
   * chart is hardest to measure. A run has no such assumption in it.
   */
  const int lo = MAX(1, margin);
  const int hi = MIN(n - 1, n - margin);

  gboolean in_run = FALSE;
  float best_pos = 0.0f, best_depth = 0.0f;
  gboolean have_best = FALSE;

  for(int i = lo; i < hi; i++)
  {
    const float prev = sig[(size_t)(i - 1) * stride];
    const float cur = sig[(size_t)i * stride];
    const float next = sig[(size_t)(i + 1) * stride];

    // a dark line is a negative excursion after background removal
    const gboolean below = cur < -threshold;

    if(below && !in_run)
    {
      in_run = TRUE;
      have_best = FALSE;
      best_depth = 0.0f;
    }

    if(below && cur <= prev && cur <= next)
    {
      /* Refine to sub-pixel by fitting a parabola through the three
         samples around the minimum. Integer peak positions would cap the
         whole calibration at about half a pixel. */
      const float denom = prev - 2.0f * cur + next;
      if(denom > 0.0f)
      {
        const float delta = 0.5f * (prev - next) / denom;
        if(delta >= -1.0f && delta <= 1.0f && (!have_best || -cur > best_depth))
        {
          best_pos = (float)i + delta;
          best_depth = -cur;
          have_best = TRUE;
        }
      }
    }

    if(!below && in_run)
    {
      in_run = FALSE;
      if(have_best && found < out_max)
      {
        out[found] = best_pos;
        strength[found] = best_depth;
        found++;
      }
    }
  }

  // a run that was still open when the scanline ended
  if(in_run && have_best && found < out_max)
  {
    out[found] = best_pos;
    strength[found] = best_depth;
    found++;
  }

  return found;
}

/* How dark a chart line actually is, in this picture.
 *
 * This wants to be a property of the lines, not of the image as a whole. An
 * earlier version used the mean absolute deviation of the whole frame and
 * set the bar at a fraction of it, which is a bad measure twice over: the
 * MAD depends on how much of the frame the lines happen to cover, and it
 * sits far below the lines themselves, so the fraction that survived the
 * faintest real line also admitted the noise floor. On a real chart the
 * lines came out around thirty times the resulting threshold, and detection
 * produced some sixteen hundred noise fragments for thirteen lines.
 *
 * A high percentile of the negative excursion lands inside the line
 * population instead, whatever fraction of the frame they cover, so the
 * threshold can be set as a fraction of a real line's depth.
 */
#define DT_GRID_DEPTH_BINS 512

static float _estimate_line_depth(const float *const sig,
                                  const int width,
                                  const int height,
                                  const int margin)
{
  /* Skip the border. The frame edge and the chart's own heavy outer border
     are the highest contrast features in the picture, and letting them into
     the statistic pushes the estimate up towards them and away from the
     interior lines that have to be measured. */
  const int x0 = MIN(margin, width / 4);
  const int y0 = MIN(margin, height / 4);

  float peak = 0.0f;
  for(int y = y0; y < height - y0; y += 4)
    for(int x = x0; x < width - x0; x += 4)
      peak = MAX(peak, -sig[(size_t)y * width + x]);

  if(peak <= 0.0f) return 0.0f;

  // a histogram rather than a sort: the percentile is all that is needed
  size_t hist[DT_GRID_DEPTH_BINS] = { 0 };
  size_t total = 0;
  const float scale = (DT_GRID_DEPTH_BINS - 1) / peak;

  for(int y = y0; y < height - y0; y += 4)
    for(int x = x0; x < width - x0; x += 4)
    {
      const float depth = MAX(0.0f, -sig[(size_t)y * width + x]);
      const int bin = MIN(DT_GRID_DEPTH_BINS - 1, (int)(depth * scale));
      hist[bin]++;
      total++;
    }

  if(!total) return 0.0f;

  const size_t want = (size_t)(0.95 * total);
  size_t seen = 0;
  for(int b = 0; b < DT_GRID_DEPTH_BINS; b++)
  {
    seen += hist[b];
    if(seen >= want) return (b + 0.5f) / scale;
  }

  return peak;
}

/* Link crossings from consecutive scanlines into curves.
 *
 * Walks the scan axis once. For every crossing it either continues the
 * nearest existing open curve, or starts a new one. `spacing` is the
 * expected distance between adjacent grid lines and sets how far a line
 * may wander between scan positions before we assume it is a different
 * line -- distortion bends the lines, but only gradually.
 */
static GList *_trace_curves(const float *const sig,
                            const int width,
                            const int height,
                            const gboolean horizontal,
                            const float threshold,
                            const float spacing,
                            const int margin,
                            int *total_crossings)
{
  // scanning by column traces the horizontal lines, and vice versa
  const int scan_n = horizontal ? width : height;   // number of scan positions
  const int line_n = horizontal ? height : width;   // samples per scanline
  const int line_stride = horizontal ? width : 1;
  const int scan_stride = horizontal ? 1 : width;

  /* Cap crossings per scanline. A noisy or badly exposed frame can
     otherwise produce thousands of spurious extrema, and the linking
     below is quadratic in the number of open curves. */
  const int expected_lines =
    (int)(line_n / MAX(1.0f, spacing)) + 2;
  const int max_cross = MIN(line_n, MAX(16, expected_lines * 4));
  float *cross = malloc(sizeof(float) * max_cross);
  if(!cross) return NULL;

  float *strength = malloc(sizeof(float) * max_cross);
  if(!strength)
  {
    free(cross);
    return NULL;
  }

  /* Scratch for the assignment below, sized once for the worst case. */
  const int max_open = max_cross + 8;
  dt_lens_grid_curve_t **openv =
    malloc(sizeof(dt_lens_grid_curve_t *) * max_open);
  double *dp = malloc(sizeof(double) * (size_t)(max_open + 1) * (max_cross + 1));
  char *bt = malloc(sizeof(char) * (size_t)(max_open + 1) * (max_cross + 1));
  int *match = malloc(sizeof(int) * max_cross);

  if(!openv || !dp || !bt || !match)
  {
    free(cross);
    free(strength);
    free(openv);
    free(dp);
    free(bt);
    free(match);
    return NULL;
  }

  GList *open = NULL;   // curves still being extended
  GList *done = NULL;   // curves that ended
  /* Never let the allowance reach a neighbouring line, however long the
     gap: past half the spacing a "continuation" is at least as likely to be
     the next line over, and joining those is worse than losing the line. */
  const float tol_cap = MAX(3.0f, spacing * 0.45f);

  // don't scan right up to the frame edge either; the same border
  // structure that inflates the threshold also traces as a "line"
  const int scan_lo = MIN(margin, scan_n / 4);
  const int scan_hi = MAX(scan_lo + 1, scan_n - scan_lo);

  for(int s = scan_lo; s < scan_hi; s += DT_GRID_SCAN_STEP)
  {
    const float *scanline = sig + (size_t)s * scan_stride;
    const int nc = _find_crossings(scanline, line_n, line_stride,
                                   threshold, margin, cross, strength,
                                   max_cross);
    if(total_crossings) *total_crossings += nc;

    /* Collect the open curves in the same order their last samples run
       across the scanline. Both sequences are then sorted, which is what
       lets the assignment below be a simple monotone alignment. */
    int nopen = 0;
    for(GList *it = open; it && nopen < max_open; it = g_list_next(it))
      openv[nopen++] = it->data;

    for(int a = 1; a < nopen; a++)
    {
      dt_lens_grid_curve_t *key = openv[a];
      const float kv = key->samples[key->count - 1].coord;
      int b = a - 1;
      while(b >= 0 && openv[b]->samples[openv[b]->count - 1].coord > kv)
      {
        openv[b + 1] = openv[b];
        b--;
      }
      openv[b + 1] = key;
    }

    /* Order preserving assignment of crossings to open curves.
     *
     * Grid lines of the same family never cross one another, so the order
     * of the lines along a scanline is the same at every scan position.
     * Matching each crossing to its individually nearest curve throws that
     * away, and on a chart with a few spurious crossings it lets a phantom
     * capture a real line's continuation, splitting the line in two. A
     * monotone alignment cannot produce an order violating match at all,
     * so a spurious crossing costs one skipped position instead of
     * derailing its neighbours.
     *
     * Each curve carries its own tolerance, grown by the slope allowance
     * over however many scan positions it has been missing: a line last
     * seen six steps ago may legitimately have drifted six steps' worth,
     * one seen at the previous step may not have. Match costs are divided
     * by that tolerance so they are comparable across curves, which makes
     * the cost of skipping either side simply 1.
     */
    const int stride_dp = nc + 1;
    dp[0] = 0.0;
    bt[0] = 0;
    for(int j = 1; j <= nc; j++)
    {
      dp[j] = dp[j - 1] + 1.0;
      bt[j] = 2; // came from skipping a crossing
    }

    for(int i = 1; i <= nopen; i++)
    {
      const dt_lens_grid_curve_t *cv = openv[i - 1];

      // how far this particular line could have moved since we last saw it
      const float travelled = MAX((float)DT_GRID_SCAN_STEP,
                                  (float)s - cv->samples[cv->count - 1].pos);

      /* Extrapolate along the curve's own direction and measure the
         deviation from that, not from where it last was. A curve too young
         to have a direction falls back to the plain slope bound. */
      const gboolean young = cv->count < 3;
      const float last = cv->samples[cv->count - 1].coord
        + (young ? 0.0f : _curve_slope(cv) * travelled);

      const float rate = young ? DT_GRID_MAX_SLOPE : DT_GRID_MAX_CURVE;
      const float tol = MIN(tol_cap, DT_GRID_LINK_BASE + rate * travelled);

      dp[(size_t)i * stride_dp] = dp[(size_t)(i - 1) * stride_dp] + 1.0;
      bt[(size_t)i * stride_dp] = 1; // came from skipping a curve

      for(int j = 1; j <= nc; j++)
      {
        double best = dp[(size_t)(i - 1) * stride_dp + j] + 1.0;
        char how = 1;

        const double skip_cross = dp[(size_t)i * stride_dp + j - 1] + 1.0;
        if(skip_cross < best)
        {
          best = skip_cross;
          how = 2;
        }

        const float d = fabsf(last - cross[j - 1]);
        if(d <= tol)
        {
          const double pair =
            dp[(size_t)(i - 1) * stride_dp + j - 1] + (double)(d / tol);
          if(pair < best)
          {
            best = pair;
            how = 3;
          }
        }

        dp[(size_t)i * stride_dp + j] = best;
        bt[(size_t)i * stride_dp + j] = how;
      }
    }

    for(int j = 0; j < nc; j++) match[j] = -1;

    /* Count a miss against every open curve up front, then clear it on the
       ones that get matched. Walking the list rather than openv matters
       when there are more open curves than the scratch holds: those are
       skipped by the assignment, and if they never aged they would never
       retire either. */
    for(GList *it = open; it; it = g_list_next(it))
      ((dt_lens_grid_curve_t *)it->data)->misses++;

    {
      int i = nopen, j = nc;
      while(i > 0 || j > 0)
      {
        const char how = bt[(size_t)i * stride_dp + j];
        if(how == 3)
        {
          match[j - 1] = i - 1;
          i--;
          j--;
        }
        else if(how == 1)
          i--;
        else
          j--;
      }
    }

    for(int j = 0; j < nc; j++)
    {
      if(match[j] >= 0)
      {
        dt_lens_grid_curve_t *cv = openv[match[j]];
        _curve_append(cv, (float)s, cross[j]);
        cv->misses = 0;
      }
      else
      {
        dt_lens_grid_curve_t *cv = _curve_new(horizontal);
        if(cv)
        {
          _curve_append(cv, (float)s, cross[j]);
          open = g_list_prepend(open, cv);
        }
      }
    }

    /* A curve is only retired after several consecutive misses. Grid
       lines cross each other, go soft, and dip under the noise floor, so
       a crossing is routinely missed at the odd scan position -- retiring
       on the first miss shreds every line into fragments too short to
       survive the span test. */
    GList *it = open;
    while(it)
    {
      GList *next = g_list_next(it);
      dt_lens_grid_curve_t *cv = it->data;
      if(cv->misses > DT_GRID_MAX_GAP)
      {
        open = g_list_remove(open, cv);
        done = g_list_prepend(done, cv);
      }
      it = next;
    }
  }

  free(strength);
  free(openv);
  free(dp);
  free(bt);
  free(match);

  // whatever is still open ran to the edge of the image
  for(GList *it = open; it; it = g_list_next(it))
    done = g_list_prepend(done, it->data);
  g_list_free(open);
  free(cross);

  return done;
}

/* Rejoin fragments of the same line.
 *
 * A chart line is routinely interrupted -- it passes behind something, dips
 * under the noise floor in a dark corner, or crosses a heavier line that
 * swallows it for a stretch longer than the tracer will wait. Each
 * interruption ends one curve and starts another, and the pieces are then
 * individually too short to survive the span test, so the line is lost
 * entirely rather than merely dented. Stitching the pieces back together
 * before pruning is what turns a partial detection into a usable one.
 *
 * Two fragments are the same line if the second starts after the first
 * ends, the gap is not enormous, and the two ends agree on where the line
 * was -- that last condition is what stops neighbouring lines being
 * spliced into a zigzag.
 */
static gint _curve_compare_start(gconstpointer a, gconstpointer b)
{
  const dt_lens_grid_curve_t *ca = a;
  const dt_lens_grid_curve_t *cb = b;
  const float pa = ca->samples[0].pos;
  const float pb = cb->samples[0].pos;
  return (pa > pb) - (pa < pb);
}

static GList *_merge_fragments(GList *curves,
                               const float tol,
                               const float max_pos_gap)
{
  curves = g_list_sort(curves, _curve_compare_start);

  GList *kept = NULL;

  for(GList *it = curves; it; it = g_list_next(it))
  {
    dt_lens_grid_curve_t *c = it->data;

    dt_lens_grid_curve_t *best = NULL;
    float best_diff = tol;

    for(GList *k = kept; k; k = g_list_next(k))
    {
      dt_lens_grid_curve_t *t = k->data;

      const float t_end = t->samples[t->count - 1].pos;
      const float c_start = c->samples[0].pos;
      const float gap = c_start - t_end;

      if(gap <= 0.0f || gap > max_pos_gap) continue;

      /* Same slope reasoning as the tracer: over a long gap a line may
         legitimately have moved a long way, over a short one it may not.
         A flat tolerance would either refuse real rejoins across a wide
         interruption or splice neighbouring lines across a narrow one. */
      const float allow = MIN(tol, DT_GRID_LINK_BASE
                                     + DT_GRID_MAX_SLOPE * gap);

      const float diff = fabsf(t->samples[t->count - 1].coord
                               - c->samples[0].coord);
      if(diff <= allow && diff < best_diff)
      {
        best_diff = diff;
        best = t;
      }
    }

    if(best)
    {
      for(int i = 0; i < c->count; i++)
        _curve_append(best, c->samples[i].pos, c->samples[i].coord);
      _curve_free(c);
    }
    else
      kept = g_list_prepend(kept, c);
  }

  g_list_free(curves);
  return kept;
}

/* Where a curve sits at the middle of the frame.
 *
 * Used for ordering and for spacing estimates. Comparing curves by their
 * own middle sample instead would compare them at different places, which
 * on a distorted or perspective-skewed chart is enough to shuffle two
 * neighbouring lines into the wrong order. */
static float _curve_coord_at(const dt_lens_grid_curve_t *cv,
                             const float pos)
{
  if(cv->count < 1) return 0.0f;

  const float lo = cv->samples[0].pos;
  const float hi = cv->samples[cv->count - 1].pos;

  if(pos <= lo) return cv->samples[0].coord;
  if(pos >= hi) return cv->samples[cv->count - 1].coord;

  for(int i = 1; i < cv->count; i++)
    if(pos <= cv->samples[i].pos)
    {
      const float p0 = cv->samples[i - 1].pos, p1 = cv->samples[i].pos;
      const float c0 = cv->samples[i - 1].coord, c1 = cv->samples[i].coord;
      const float t = (p1 > p0) ? (pos - p0) / (p1 - p0) : 0.0f;
      return c0 + t * (c1 - c0);
    }

  return cv->samples[cv->count - 1].coord;
}

/* Keep only the stretch of lattice that the chart could actually occupy.
 *
 * The user has said how big the chart is, and once the indices are snapped
 * to a lattice that becomes usable information: a line at index 34 of a
 * chart with 29 of them is something else in the picture -- a door frame, a
 * shadow, the edge of the board. Which end the surplus sits at is not
 * knowable in advance, so slide a window of the right width across the
 * indices and keep the position holding the most measurement. Sample counts
 * rather than curve counts, because a long dense line is better evidence of
 * where the chart is than two short sparse ones.
 */
static GList *_select_window(GList *kept, const int expected)
{
  int omax = 0;
  for(GList *it = kept; it; it = g_list_next(it))
    omax = MAX(omax, ((dt_lens_grid_curve_t *)it->data)->ordinal);

  if(omax <= expected - 1) return kept;

  int best_off = 0;
  long best_score = -1;

  for(int off = 0; off <= omax - (expected - 1); off++)
  {
    long score = 0;
    for(GList *it = kept; it; it = g_list_next(it))
    {
      const dt_lens_grid_curve_t *cv = it->data;
      if(cv->ordinal >= off && cv->ordinal <= off + expected - 1)
        score += cv->count;
    }
    if(score > best_score)
    {
      best_score = score;
      best_off = off;
    }
  }

  GList *it = kept;
  int dropped = 0;
  while(it)
  {
    GList *next = g_list_next(it);
    dt_lens_grid_curve_t *cv = it->data;

    if(cv->ordinal < best_off || cv->ordinal > best_off + expected - 1)
    {
      kept = g_list_delete_link(kept, it);
      _curve_free(cv);
      dropped++;
    }
    else
      cv->ordinal -= best_off;

    it = next;
  }

  if(dropped)
    dt_print(DT_DEBUG_ALWAYS,
             "[lens_grid] lattice indices reached %d for a chart with %d"
             " lines; kept the %d starting at %d and dropped %d line(s)"
             " that cannot be part of it",
             omax, expected, expected, best_off, dropped);

  return kept;
}

/* Drop curves that are too short to be a chart line, and sort what is
   left across the frame so the ordinals run in geometric order. */
static gint _curve_compare(gconstpointer a, gconstpointer b, gpointer data)
{
  const dt_lens_grid_curve_t *ca = a;
  const dt_lens_grid_curve_t *cb = b;
  const float centre = *(const float *)data;

  // compare both curves at the same place, not at their own midpoints
  const float ma = _curve_coord_at(ca, centre);
  const float mb = _curve_coord_at(cb, centre);
  return (ma > mb) - (ma < mb);
}

static gint _curve_compare_len(gconstpointer a, gconstpointer b)
{
  const dt_lens_grid_curve_t *ca = a;
  const dt_lens_grid_curve_t *cb = b;
  return cb->count - ca->count; // longest first
}

static GList *_prune_and_order(GList *curves,
                               const int scan_extent,
                               const int expected,
                               const float spacing)
{
  const float min_span = scan_extent * DT_GRID_MIN_SPAN;
  GList *kept = NULL;

  for(GList *it = curves; it; it = g_list_next(it))
  {
    dt_lens_grid_curve_t *cv = it->data;
    const float span = cv->count > 1
      ? fabsf(cv->samples[cv->count - 1].pos - cv->samples[0].pos)
      : 0.0f;

    // how many samples a line covering this span should have produced
    const float possible = span / (float)DT_GRID_SCAN_STEP;
    const float density = (possible > 0.0f) ? cv->count / possible : 0.0f;

    if(cv->count >= 4 && span >= min_span && density >= DT_GRID_MIN_DENSITY)
      kept = g_list_prepend(kept, cv);
    else
      _curve_free(cv);
  }
  g_list_free(curves);

  /* Keep the intersection step bounded. If detection went wrong we would
     rather work with the most convincing lines than grind through a few
     hundred spurious ones -- the longest traces are the real chart. */
  const guint before_cap = g_list_length(kept);
  const int cap = MAX(4, expected * 3);

  if(before_cap > (guint)cap)
    dt_print(DT_DEBUG_ALWAYS,
             "[lens_grid] %u lines survived pruning for an expected %d;"
             " keeping the %d longest",
             before_cap, expected, cap);

  if(before_cap > (guint)cap)
  {
    kept = g_list_sort(kept, _curve_compare_len);
    while(g_list_length(kept) > (guint)cap)
    {
      GList *last = g_list_last(kept);
      _curve_free(last->data);
      kept = g_list_delete_link(kept, last);
    }
  }

  float centre = 0.5f * scan_extent;
  kept = g_list_sort_with_data(kept, _curve_compare, &centre);

  /* Collapse duplicates.
   *
   * A wide chart line whose trough splits into two runs consistently is
   * traced twice, giving two long dense curves a couple of pixels apart.
   * Both are real detections of the same line, and keeping both would
   * corrupt the lattice spacing that the indices below are derived from.
   * The one with more samples is the better measurement of the two.
   */
  const float too_close = MAX(2.0f, spacing * DT_GRID_DUPLICATE);

  GList *dup = kept;
  while(dup && g_list_next(dup))
  {
    dt_lens_grid_curve_t *a = dup->data;
    dt_lens_grid_curve_t *b = g_list_next(dup)->data;

    if(fabsf(_curve_coord_at(b, centre) - _curve_coord_at(a, centre))
       < too_close)
    {
      // drop the sparser of the pair and look at this position again
      if(b->count > a->count)
      {
        kept = g_list_remove(kept, a);
        _curve_free(a);
        dup = g_list_find(kept, b);
        if(dup && g_list_previous(dup)) dup = g_list_previous(dup);
      }
      else
      {
        kept = g_list_remove(kept, b);
        _curve_free(b);
      }
      continue;
    }

    dup = g_list_next(dup);
  }

  /* Assign lattice indices by *position*, not by sequence.
   *
   * Numbering the kept curves 0,1,2,... assumes every line was found. Miss
   * one and every line beyond it is labelled one too low, which silently
   * corrupts the lattice the solver is handed: the plumb-line part only
   * cares which points share a line and survives it, but the stage that
   * recovers the squeeze compares the measured grid against an ideal one
   * built from these indices, and a compressed ideal grid comes back as a
   * wrong squeeze. Since a gap of two lines is twice the gap of one, the
   * step can be read off the spacing instead of assumed.
   */
  const guint n = g_list_length(kept);

  if(n >= 3)
  {
    float *gaps = malloc(sizeof(float) * (n - 1));
    if(gaps)
    {
      guint g = 0;
      const dt_lens_grid_curve_t *prev = NULL;
      for(GList *it = kept; it; it = g_list_next(it))
      {
        const dt_lens_grid_curve_t *cv = it->data;
        if(prev)
          gaps[g++] = fabsf(_curve_coord_at(cv, centre)
                            - _curve_coord_at(prev, centre));
        prev = cv;
      }

      // the median gap is the spacing of adjacent lines: a missed line
      // shows up as an outlying large gap and cannot shift the median
      float *sorted = malloc(sizeof(float) * g);
      if(sorted)
      {
        memcpy(sorted, gaps, sizeof(float) * g);
        for(guint i = 1; i < g; i++)
        {
          const float key = sorted[i];
          guint j = i;
          while(j > 0 && sorted[j - 1] > key)
          {
            sorted[j] = sorted[j - 1];
            j--;
          }
          sorted[j] = key;
        }
        const float unit = sorted[g / 2];
        free(sorted);

        if(unit > 1.0f)
        {
          int ordinal = 0;
          guint k = 0;
          for(GList *it = kept; it; it = g_list_next(it))
          {
            dt_lens_grid_curve_t *cv = it->data;
            if(k > 0)
            {
              /* At least one step, however small the gap: two curves this
                 close are more likely one line detected twice than a
                 lattice with no room between its lines, and collapsing
                 them onto the same index would be worse. */
              const int step = MAX(1, (int)lroundf(gaps[k - 1] / unit));
              ordinal += step;
            }
            cv->ordinal = ordinal;
            k++;
          }

          free(gaps);
          return _select_window(kept, expected);
        }
      }
      free(gaps);
    }
  }

  // too few lines to estimate a spacing; sequential is all there is
  int ordinal = 0;
  for(GList *it = kept; it; it = g_list_next(it))
  {
    dt_lens_grid_curve_t *cv = it->data;
    cv->ordinal = ordinal++;
  }

  return _select_window(kept, expected);
}

/* Evaluate a traced curve at a given scan position by linear
   interpolation between its samples. Returns FALSE outside its span. */
static gboolean _curve_eval(const dt_lens_grid_curve_t *cv,
                            const float pos,
                            const float margin,
                            float *coord)
{
  if(cv->count < 2) return FALSE;

  const float lo = cv->samples[0].pos;
  const float hi = cv->samples[cv->count - 1].pos;

  /* A short way past either end the curve is extended flat rather than
     refused. Two curves that genuinely cross near the frame edge often
     have one of them stopping a few scan positions short of the crossing,
     and rejecting those loses precisely the corner measurements that
     constrain the distortion most. Extrapolating the slope would be
     worse than holding the value: it amplifies the noise in the last
     sample right where the curve is bending hardest. */
  if(pos < lo - margin || pos > hi + margin) return FALSE;
  if(pos <= lo)
  {
    *coord = cv->samples[0].coord;
    return TRUE;
  }
  if(pos >= hi)
  {
    *coord = cv->samples[cv->count - 1].coord;
    return TRUE;
  }

  for(int i = 1; i < cv->count; i++)
  {
    if(pos <= cv->samples[i].pos)
    {
      const float p0 = cv->samples[i - 1].pos, p1 = cv->samples[i].pos;
      const float c0 = cv->samples[i - 1].coord, c1 = cv->samples[i].coord;
      const float t = (p1 > p0) ? (pos - p0) / (p1 - p0) : 0.0f;
      *coord = c0 + t * (c1 - c0);
      return TRUE;
    }
  }
  return FALSE;
}

/* Intersect every horizontal curve with every vertical one.
 *
 * A horizontal curve gives y as a function of x, a vertical curve gives
 * x as a function of y, so the crossing is the fixed point of the pair.
 * Two iterations of substitution converge well within a pixel because
 * the curves are nearly axis aligned.
 */
static void _intersect_curves(dt_lens_grid_t *grid)
{
  const int nh = g_list_length(grid->curves_h);
  const int nv = g_list_length(grid->curves_v);
  if(!nh || !nv) return;

  grid->points = calloc((size_t)nh * nv, sizeof(dt_lens_grid_point_t));
  if(!grid->points) return;
  grid->point_count = 0;

  const float margin = 0.03f * MAX(grid->width, grid->height);

  for(GList *ih = grid->curves_h; ih; ih = g_list_next(ih))
  {
    const dt_lens_grid_curve_t *ch = ih->data;
    for(GList *iv = grid->curves_v; iv; iv = g_list_next(iv))
    {
      const dt_lens_grid_curve_t *cv = iv->data;

      // seed from the curve midpoints
      float x = cv->samples[cv->count / 2].coord;
      float y = ch->samples[ch->count / 2].coord;
      gboolean ok = TRUE;

      for(int iter = 0; iter < 3 && ok; iter++)
      {
        float ny, nx;
        ok = _curve_eval(ch, x, margin, &ny)
          && _curve_eval(cv, ny, margin, &nx);
        if(ok)
        {
          y = ny;
          x = nx;
        }
      }

      if(!ok) continue;
      if(x < 0 || y < 0 || x >= grid->width || y >= grid->height) continue;

      dt_lens_grid_point_t *p = &grid->points[grid->point_count++];
      p->x = x;
      p->y = y;
      p->col = cv->ordinal;
      p->row = ch->ordinal;
    }
  }
}

gboolean dt_lens_grid_detect(const float *const lum,
                             const int width,
                             const int height,
                             const int cells_x,
                             const int cells_y,
                             const gboolean dark_lines,
                             dt_lens_grid_t *grid)
{
  if(!lum || !grid || width < 16 || height < 16) return FALSE;
  if(cells_x < 1 || cells_y < 1) return FALSE;

  memset(grid, 0, sizeof(dt_lens_grid_t));
  grid->width = width;
  grid->height = height;
  grid->cells_x = cells_x;
  grid->cells_y = cells_y;

  // expected spacing between lines, used for the link tolerance and to
  // size the normalization window
  const float spacing_x = (float)width / (float)(cells_x + 1);
  const float spacing_y = (float)height / (float)(cells_y + 1);

  // the background window has to be wide enough to span a cell,
  // otherwise the lines suppress each other
  const int radius = (int)fmaxf(4.0f, fminf(spacing_x, spacing_y));

  float *norm = _normalize_local(lum, width, height, radius);
  if(!norm) return FALSE;

  // work in "dark line" convention internally
  if(!dark_lines)
  {
    DT_OMP_FOR()
    for(size_t i = 0; i < (size_t)width * height; i++)
      norm[i] = -norm[i];
  }

  /* The mean absolute deviation is pulled up by high contrast features
     such as the chart's own border, so a factor near 1 sets the bar by
     the border and rejects the much fainter interior lines -- the exact
     opposite of what we want. Keep the factor well below that. */
  /* Keep away from the frame edge. The chart's heavy outer border and the
     frame boundary itself are the strongest features present, and they
     trace as perfectly good "lines" that belong to no lattice. */
  const int margin = MAX(4, (int)(0.012f * MAX(width, height)));

  /* A quarter of a real line's depth. Low enough to keep the fainter lines
     in a vignetted corner -- those measured about half the depth of the
     ones in the middle -- and high enough that the noise floor, an order of
     magnitude below either, never gets in. */
  const float depth = _estimate_line_depth(norm, width, height, margin);
  const float threshold = MAX(1e-6f, depth * 0.25f);

  int cross_h = 0, cross_v = 0;
  GList *raw_h = _trace_curves(norm, width, height, TRUE, threshold,
                               spacing_y, margin, &cross_h);
  GList *raw_v = _trace_curves(norm, width, height, FALSE, threshold,
                               spacing_x, margin, &cross_v);

  dt_free_align(norm);

  /* Stitch fragments before pruning, so a line broken into pieces gets a
     chance to be long enough to keep. The gap allowance is generous
     because the coordinate agreement is doing the real filtering. */
  raw_h = _merge_fragments(raw_h, MAX(3.0f, spacing_y * 0.45f),
                           0.25f * width);
  raw_v = _merge_fragments(raw_v, MAX(3.0f, spacing_x * 0.45f),
                           0.25f * height);

  grid->curves_h = _prune_and_order(raw_h, width, cells_y + 1, spacing_y);
  grid->curves_v = _prune_and_order(raw_v, height, cells_x + 1, spacing_x);

  _intersect_curves(grid);

  const int nh = g_list_length(grid->curves_h);
  const int nv = g_list_length(grid->curves_v);

  // the highest lattice index reached says whether the ordinal snapping
  // found gaps, which a bare count cannot show
  int span_h = 0, span_v = 0;
  for(GList *it = grid->curves_h; it; it = g_list_next(it))
    span_h = MAX(span_h, ((dt_lens_grid_curve_t *)it->data)->ordinal);
  for(GList *it = grid->curves_v; it; it = g_list_next(it))
    span_v = MAX(span_v, ((dt_lens_grid_curve_t *)it->data)->ordinal);

  dt_print(DT_DEBUG_ALWAYS,
           "[lens_grid] %dx%d px, radius %d, margin %d, line depth %.5f,"
           " threshold %.5f: %d/%d raw crossings, traced %d horizontal"
           " and %d vertical lines (expected %d/%d), spanning lattice"
           " indices 0..%d and 0..%d, %d intersections",
           width, height, radius, margin, depth, threshold,
           cross_h, cross_v, nh, nv, cells_y + 1, cells_x + 1,
           span_h, span_v, grid->point_count);

  // two lines in each direction is the minimum that describes anything
  return nh >= 2 && nv >= 2;
}

void dt_lens_grid_cleanup(dt_lens_grid_t *grid)
{
  if(!grid) return;

  g_list_free_full(grid->curves_h, _curve_free);
  g_list_free_full(grid->curves_v, _curve_free);
  grid->curves_h = grid->curves_v = NULL;

  free(grid->points);
  grid->points = NULL;
  grid->point_count = 0;
}

// clang-format off
// modelines: These editor modelines have been set for all relevant files by tools/update_modelines.py
// vim: shiftwidth=2 expandtab tabstop=2 cindent
// kate: tab-indents: off; indent-width 2; replace-tabs on; indent-mode cstyle; remove-trailing-spaces modified;
// clang-format on
