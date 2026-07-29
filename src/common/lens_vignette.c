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

#include "common/lens_vignette.h"

#include "common/darktable.h"
#include "iop/gaussian_elimination.h"

#include <math.h>
#include <stdlib.h>
#include <string.h>

#define DT_VIG_BINS 48
#define DT_VIG_MIN_PER_BIN 32

/* Which value in a bin stands for "the surface".
 *
 * Not the mean and not the median. A calibration frame is usually a chart or
 * a wall with something dark on it, and any measure of central tendency mixes
 * the dark marks into the reading -- worse, it mixes in *more* of them where
 * the marks are denser, which is a spatial pattern that looks exactly like
 * vignetting. Taking a high percentile picks the paper and ignores the ink.
 * Not the maximum, which would pick specular flecks and sensor hot pixels.
 */
#define DT_VIG_PERCENTILE 0.80f

/* How many radius bands the cell samples are collapsed into. Enough to
   describe a smooth falloff, few enough that the outermost band still holds
   several cells on a normal chart. */
#define DT_VIG_RBINS 20

static int _cmp_float(const void *a, const void *b)
{
  const float x = *(const float *)a, y = *(const float *)b;
  return (x < y) ? -1 : ((x > y) ? 1 : 0);
}

static int _cmp_double(const void *a, const void *b)
{
  const double x = *(const double *)a, y = *(const double *)b;
  return (x < y) ? -1 : ((x > y) ? 1 : 0);
}

/* Does the gain only ever rise, out to `r2_end`?
 *
 * Vignetting is monotonic in radius -- no lens is brighter at 80% of the
 * frame than at 60% -- so a fitted curve that turns round has stopped
 * describing the lens and started describing the gap between the last
 * measurement and the frame corner. Checking past the data is the point:
 * that region is exactly where an unconstrained sixth order term misbehaves,
 * and where the correction is largest and least forgiving.
 */
/* Is this a correction a lens could actually need?
 *
 * Independent of how the coefficients were arrived at, and applied to every
 * candidate before it is accepted. No lens loses four stops in the corner, so
 * a fit that says so is describing something other than the lens -- and the
 * failure is worth catching here rather than at the point where a 400x gain
 * turns corner noise into white. A gain below one at the corner means the
 * edges came out brighter than the middle, which is a lighting problem, not
 * a lens.
 */
static gboolean _plausible(const float k[3], const double r2_corner)
{
  const double r = r2_corner;
  const double g = 1.0 + k[0] * r + k[1] * r * r + k[2] * r * r * r;
  return g >= 0.95 && g <= 16.0;
}

/* Plain, unweighted residual of a model against the band values, as a
 * fraction of centre brightness.
 *
 * Kept separate from the objective the fit minimises. The objective is
 * weighted towards the frame edge on purpose, which makes it the right thing
 * to choose a model by and the wrong thing to report: a number that moves
 * because the weighting changed, rather than because the fit did, cannot be
 * compared against a fixed "is this believable" threshold. This is what gets
 * shown.
 */
static double _plain_rms(const double *const r2,
                         const double *const val,
                         const int n,
                         const float k[3])
{
  // the scale that maps the measurements onto the model's own normalization
  double num = 0.0;
  int m = 0;
  for(int i = 0; i < n; i++)
  {
    const double rr = r2[i];
    const double g = 1.0 + k[0] * rr + k[1] * rr * rr + k[2] * rr * rr * rr;
    if(g > 1e-6 && val[i] > 1e-12)
    {
      num += g / val[i];
      m++;
    }
  }
  if(!m) return 1e30;

  const double c0 = num / m;

  double acc = 0.0;
  for(int i = 0; i < n; i++)
  {
    const double rr = r2[i];
    const double g = 1.0 + k[0] * rr + k[1] * rr * rr + k[2] * rr * rr * rr;
    const double d = 1.0 / MAX(1e-6, g) - val[i] * c0;
    acc += d * d;
  }

  return sqrt(acc / n);
}

static gboolean _monotonic(const float k[3], const double r2_end)
{
  double prev = 1.0;

  for(int i = 1; i <= 64; i++)
  {
    const double r = r2_end * i / 64.0;
    const double g = 1.0 + k[0] * r + k[1] * r * r + k[2] * r * r * r;

    // a small tolerance, so numerical flatness is not read as a downturn
    if(g < prev - 1e-5) return FALSE;
    prev = g;
  }

  return prev >= 1.0;
}

/* Bin the frame by elliptical radius and reduce each bin to one number.
   Returns the number of usable bins, filling r2/val. */
static int _profile(const float *const lum,
                    const int width, const int height, const int stride,
                    const double cx, const double cy, const double ex,
                    double *r2out, double *valout)
{
  const double hd = 0.5 * hypot((double)width, (double)height);

  /* The largest radius any pixel reaches, so the bins span the frame rather
     than a circle inscribed in it -- the corners are where vignetting lives
     and dropping them would fit the flat part only. */
  double r2max = 0.0;
  {
    const double corner[4][2] = { { 0, 0 }, { 1, 0 }, { 0, 1 }, { 1, 1 } };
    for(int i = 0; i < 4; i++)
    {
      const double u = (corner[i][0] - 0.5) * width / hd - cx;
      const double v = (corner[i][1] - 0.5) * height / hd - cy;
      const double e = u / ex;
      r2max = MAX(r2max, e * e + v * v);
    }
  }
  if(!(r2max > 1e-6)) return 0;

  float **bin = calloc(DT_VIG_BINS, sizeof(float *));
  int *count = calloc(DT_VIG_BINS, sizeof(int));
  int *cap = calloc(DT_VIG_BINS, sizeof(int));
  if(!bin || !count || !cap)
  {
    free(bin);
    free(count);
    free(cap);
    return 0;
  }

  /* Subsample: a 60 megapixel frame is far more data than a 48 bin radial
     profile can use, and reading all of it costs seconds for no gain. */
  const int step = MAX(1, (int)(sqrt((double)width * height / 400000.0)));

  for(int y = 0; y < height; y += step)
  {
    const float *row = lum + (size_t)y * stride;
    for(int x = 0; x < width; x += step)
    {
      const float value = row[x];
      // a non-positive sample carries no ratio; 1/v is what gets fitted
      if(!(value > 1e-8f) || !isfinite(value)) continue;

      const double u = ((double)x + 0.5 - 0.5 * width) / hd - cx;
      const double v = ((double)y + 0.5 - 0.5 * height) / hd - cy;
      const double e = u / ex;
      const double r2 = e * e + v * v;

      int b = (int)(DT_VIG_BINS * r2 / r2max);
      b = CLAMP(b, 0, DT_VIG_BINS - 1);

      if(count[b] >= cap[b])
      {
        const int nc = cap[b] ? cap[b] * 2 : 64;
        float *nb = realloc(bin[b], sizeof(float) * nc);
        if(!nb) continue;
        bin[b] = nb;
        cap[b] = nc;
      }
      bin[b][count[b]++] = value;
    }
  }

  int n = 0;
  for(int b = 0; b < DT_VIG_BINS; b++)
  {
    if(count[b] >= DT_VIG_MIN_PER_BIN)
    {
      qsort(bin[b], count[b], sizeof(float), _cmp_float);
      const int k = CLAMP((int)(DT_VIG_PERCENTILE * count[b]),
                          0, count[b] - 1);

      // the bin's representative radius, at its centre
      r2out[n] = r2max * (b + 0.5) / DT_VIG_BINS;
      valout[n] = bin[b][k];
      n++;
    }
    free(bin[b]);
  }

  free(bin);
  free(count);
  free(cap);
  return n;
}

/* Fit the physical falloff family, then express it in the stored polynomial.
 *
 * This is the answer to extrapolation, and it is a better answer than any
 * choice of which cells to use. No weighting scheme can manufacture data past
 * the outermost cell; what it can do is constrain the *shape* so that the
 * stretch beyond the chart is inferred from optics rather than invented by a
 * polynomial with nothing to hold it down.
 *
 * Natural vignetting falls as cos^4 of the field angle. With tan(theta) = r/k
 * that is a gain of (1 + r^2/k^2)^2 -- monotone, finite and gently curved
 * everywhere by construction, so it cannot turn round outside the measured
 * range the way a free sixth order term can. Real lenses are steeper than
 * cos^4 because of mechanical cut-off, so one extra factor (1 + m r^2) is
 * allowed on top. Two shape parameters against the polynomial's three, and
 * both mean something.
 *
 * k is scanned; given k the rest is linear, because with c = 1 + r^2/k^2
 *
 *     1/v = (1/v0) c^2 + (m/v0) c^2 r^2
 *
 * is a two unknown least squares problem. The fitted curve is then projected
 * onto 1 + k0 r^2 + k1 r^4 + k2 r^6 across the whole frame radius, so the
 * stored format and the Lensfun mapping are untouched -- only the shape being
 * stored has changed.
 */
static gboolean _fit_cos4(const double *const r2,
                          const double *const val,
                          const double *const wt,
                          const int n,
                          const double r2_end,
                          float out_k[3],
                          double *out_rms)
{
  if(n < 5 || !(r2_end > 1e-9)) return FALSE;

  double best_kk = 0.0, best_m = 0.0, best_a = 0.0, best_rms = 1e30;
  gboolean any = FALSE;

  /* k is the focal length in units of the frame's half diagonal, so it is a
     physical quantity with a physical range: about 0.55 for a 12mm on full
     frame, upwards of 10 for a long tele. The scan has to be confined to that
     range. Allowed to go lower, the family stops being the safe extrapolator
     it was chosen for -- k = 0.22 describes a 154 degree half angle and puts
     nearly nine stops in the corner, and since the fit is free to shrink the
     centre brightness to compensate, that shape can match the measured range
     while being absurd outside it. Monotonicity was never the whole of good
     behaviour; magnitude is the other half. */
  for(int s = 0; s < 48; s++)
  {
    const double kk = 0.3 * pow(120.0 / 0.3, s / 47.0);

    double AtA[4] = { 0 }, Atb[2] = { 0 };

    for(int i = 0; i < n; i++)
    {
      if(!(val[i] > 1e-12)) return FALSE;

      const double c = 1.0 + r2[i] / kk;
      const double b[2] = { c * c, c * c * r2[i] };
      const double y = 1.0 / val[i];
      const double w = wt ? wt[i] : 1.0;

      AtA[0] += w * b[0] * b[0];
      AtA[1] += w * b[0] * b[1];
      AtA[2] += w * b[1] * b[0];
      AtA[3] += w * b[1] * b[1];
      Atb[0] += w * b[0] * y;
      Atb[1] += w * b[1] * y;
    }

    double sol[2] = { Atb[0], Atb[1] };
    double M[4] = { AtA[0], AtA[1], AtA[2], AtA[3] };
    if(!gauss_solve(M, sol, 2)) continue;

    const double a = sol[0];
    if(!(a > 1e-12)) continue;

    /* A negative m would mean mechanical vignetting that *brightens* the
       edges, which no aperture does. Pin it at zero and refit the scale
       alone rather than accept a shape that cannot happen. */
    double m = sol[1] / a;
    double aa = a;

    if(m < 0.0)
    {
      m = 0.0;
      double num = 0.0, den = 0.0;
      for(int i = 0; i < n; i++)
      {
        const double c = 1.0 + r2[i] / kk;
        const double b0 = c * c;
        const double w = wt ? wt[i] : 1.0;
        num += w * b0 * (1.0 / val[i]);
        den += w * b0 * b0;
      }
      if(!(den > 1e-30)) continue;
      aa = num / den;
      if(!(aa > 1e-12)) continue;
    }

    double acc = 0.0, wsum = 0.0;
    for(int i = 0; i < n; i++)
    {
      const double c = 1.0 + r2[i] / kk;
      const double g = c * c * (1.0 + m * r2[i]);
      const double model = 1.0 / MAX(1e-6, g);
      const double obs = val[i] * aa;
      const double w = wt ? wt[i] : 1.0;
      acc += w * (model - obs) * (model - obs);
      wsum += w;
    }
    if(!(wsum > 0.0)) continue;

    const double rms = sqrt(acc / wsum);
    if(rms < best_rms)
    {
      best_rms = rms;
      best_kk = kk;
      best_m = m;
      best_a = aa;
      any = TRUE;
    }
  }

  if(!any) return FALSE;
  (void)best_a;

  /* Project onto the stored polynomial over the whole frame radius, not just
     the measured part -- the point of the exercise is that the corner comes
     out right, so the corner is what the projection has to match. The
     constant is fixed at 1 by construction, so only three terms are fitted. */
  {
    double AtA[9] = { 0 }, Atb[3] = { 0 };

    for(int i = 1; i <= 96; i++)
    {
      const double rr = r2_end * i / 96.0;
      const double c = 1.0 + rr / best_kk;
      const double g = c * c * (1.0 + best_m * rr);

      const double b[3] = { rr, rr * rr, rr * rr * rr };
      const double y = g - 1.0;

      for(int a = 0; a < 3; a++)
      {
        for(int c2 = 0; c2 < 3; c2++) AtA[a * 3 + c2] += b[a] * b[c2];
        Atb[a] += b[a] * y;
      }
    }

    if(!gauss_solve(AtA, Atb, 3)) return FALSE;
    for(int i = 0; i < 3; i++) out_k[i] = (float)Atb[i];
  }

  *out_rms = best_rms;
  return TRUE;
}

/* Least squares fit of 1/v = c0 + c1 r2 + c2 r4 + c3 r6, then convert.
   Returns FALSE if the normal equations are singular or the result is not a
   usable attenuation. */
static gboolean _fit_reciprocal(const double *const r2,
                                const double *const val,
                                const int n,
                                const int nterms,
                                float out_k[3],
                                double *out_rms)
{
  const int nc = CLAMP(nterms, 1, 3) + 1; // plus the constant
  if(n < nc + 2) return FALSE;

  double AtA[16] = { 0 };
  double Atb[4] = { 0 };

  for(int i = 0; i < n; i++)
  {
    if(!(val[i] > 1e-12)) return FALSE;

    const double b[4] = { 1.0, r2[i], r2[i] * r2[i], r2[i] * r2[i] * r2[i] };
    const double y = 1.0 / val[i];

    for(int a = 0; a < nc; a++)
    {
      for(int c = 0; c < nc; c++) AtA[a * nc + c] += b[a] * b[c];
      Atb[a] += b[a] * y;
    }
  }

  if(!gauss_solve(AtA, Atb, nc)) return FALSE;

  const double c0 = Atb[0];
  // c0 is 1/v0, the reciprocal of the brightness at the centre
  if(!(c0 > 1e-12)) return FALSE;

  for(int i = 0; i < 3; i++)
    out_k[i] = (i + 1 < nc) ? (float)(Atb[i + 1] / c0) : 0.0f;

  double acc = 0.0;
  for(int i = 0; i < n; i++)
  {
    const double g = 1.0 + out_k[0] * r2[i]
                         + out_k[1] * r2[i] * r2[i]
                         + out_k[2] * r2[i] * r2[i] * r2[i];
    // compare in brightness, relative to the centre, so the number reads as
    // a fraction of full brightness rather than as a reciprocal
    const double model = 1.0 / MAX(1e-6, g);
    const double obs = val[i] * c0;
    acc += (model - obs) * (model - obs);
  }
  *out_rms = sqrt(acc / n);

  return TRUE;
}

/* Radius at the frame corner, in the elliptical space -- the point the corner
   figure is quoted at, and the one a chart usually fails to reach. */
static double _corner_r2(const int width, const int height,
                         const double cx, const double cy, const double ex)
{
  const double hd = 0.5 * hypot((double)width, (double)height);
  const double u = (0.5 * width) / hd - cx;
  const double v = (0.5 * height) / hd - cy;
  const double e = u / ex;
  return e * e + v * v;
}

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
                                     float *out_corner_stops)
{
  if(!uv || !value || count < 12) return FALSE;

  double *r2 = malloc(sizeof(double) * count);
  double *val = malloc(sizeof(double) * count);
  double *wt = malloc(sizeof(double) * count);
  // per band scratch for the medians, one row of `count` per band
  double *scratch = malloc(sizeof(double) * DT_VIG_RBINS * (size_t)count);
  if(!r2 || !val || !wt || !scratch)
  {
    free(r2);
    free(val);
    free(wt);
    free(scratch);
    return FALSE;
  }

  const float seed = (ex_hint > 0.2f && ex_hint < 5.0f) ? ex_hint : 1.0f;

  gboolean any = FALSE;
  float best_k[3] = { 0, 0, 0 };
  double best_rms = 1e30, best_max_r2 = 0.0, best_plain = 0.0;
  float best_ex = 1.0f;

  for(int s = -8; s <= 8; s++)
  {
    const double ex = seed * (1.0 + 0.06 * s);
    if(ex < 0.2 || ex > 5.0) continue;

    /* Balance the data over radius, not over area.
     *
     * The model is a function of radius alone, so every radius ought to carry
     * the same say in it. Fitting the cells directly does the opposite: a
     * rectangular chart has scores of cells at middling radii and only the
     * four corner cells out at the extreme, so an unweighted fit is decided
     * almost entirely by the middle of the frame and merely grazed by the
     * corners -- which are the part anyone cares about. Collapsing each
     * radius band to one value first gives the outermost band the same weight
     * as the innermost.
     *
     * The band value is the median rather than the mean, which also disposes
     * of the odd cell with a smudge or a speck on it: a wrong reading rather
     * than a noisy one, and one that would drag a mean with it. */
    double bsum[DT_VIG_RBINS];
    int bcnt[DT_VIG_RBINS];

    for(int b = 0; b < DT_VIG_RBINS; b++)
    {
      bsum[b] = 0.0;
      bcnt[b] = 0;
    }

    double max_r2 = 0.0;
    for(int i = 0; i < count; i++)
    {
      if(!(value[i] > 1e-8f) || !isfinite(value[i])) continue;
      const double u = (double)uv[2 * i] - cx;
      const double v = (double)uv[2 * i + 1] - cy;
      const double e = u / ex;
      max_r2 = MAX(max_r2, e * e + v * v);
    }
    if(!(max_r2 > 1e-9)) continue;

    for(int i = 0; i < count; i++)
    {
      if(!(value[i] > 1e-8f) || !isfinite(value[i])) continue;

      const double u = (double)uv[2 * i] - cx;
      const double v = (double)uv[2 * i + 1] - cy;
      const double e = u / ex;
      const double rr = e * e + v * v;

      int b = (int)(DT_VIG_RBINS * rr / max_r2);
      b = CLAMP(b, 0, DT_VIG_RBINS - 1);

      bsum[b] += rr;
      scratch[(size_t)b * count + bcnt[b]] = value[i];
      bcnt[b]++;
    }

    int n = 0;
    for(int b = 0; b < DT_VIG_RBINS; b++)
    {
      if(bcnt[b] < 1) continue;

      double *v = scratch + (size_t)b * count;
      qsort(v, bcnt[b], sizeof(double), _cmp_double);

      r2[n] = bsum[b] / bcnt[b];
      val[n] = v[bcnt[b] / 2];

      /* Weight by radius on top of the equal-per-band weighting. Both the
         area a band governs and the size of the correction there grow with
         radius, so an error at the edge costs more than the same error in the
         middle -- and the edge is what anyone looks at. */
      wt[n] = sqrt(r2[n]) + 0.05;
      n++;
    }

    if(n < 6) continue;

    const double rc = _corner_r2(width, height, cx, cy, ex);
    const double check = MAX(rc, max_r2);

    float k[3] = { 0, 0, 0 };
    double rms = 0.0;
    gboolean got = FALSE;

    /* The physical shape first, since it is the one that extrapolates -- but
       checked like everything else. Projecting it onto the polynomial can
       itself introduce a wobble, so the projection is what gets tested, not
       the family it came from. */
    {
      float kt[3];
      double rt;
      if(_fit_cos4(r2, val, wt, n, check, kt, &rt)
         && _monotonic(kt, check) && _plausible(kt, rc))
      {
        memcpy(k, kt, sizeof(k));
        rms = rt;
        got = TRUE;
      }
    }

    /* Free polynomial as a fallback, with as many terms as extrapolate
       sanely: an r^6 term unconstrained by data will happily turn round and
       dive past the last measurement, so a fit whose gain stops rising before
       the frame corner drops a term and tries again. */
    for(int nt = 3; nt >= 1 && !got; nt--)
    {
      float kt[3];
      double rt;
      if(!_fit_reciprocal(r2, val, n, nt, kt, &rt)) continue;
      if(!_monotonic(kt, check) || !_plausible(kt, rc)) continue;

      memcpy(k, kt, sizeof(k));
      rms = rt;
      got = TRUE;
    }

    if(got && rms < best_rms)
    {
      best_rms = rms;
      // reported separately from the objective the choice was made on
      best_plain = _plain_rms(r2, val, n, k);
      memcpy(best_k, k, sizeof(best_k));
      best_ex = (float)ex;
      best_max_r2 = max_r2;
      any = TRUE;
    }
  }

  free(r2);
  free(val);
  free(wt);
  free(scratch);

  if(!any) return FALSE;

  memcpy(out_k, best_k, sizeof(best_k));
  if(out_ex) *out_ex = best_ex;
  if(out_rms) *out_rms = (float)best_plain;
  if(out_max_r) *out_max_r = (float)sqrt(best_max_r2);

  if(out_corner_stops)
  {
    const double rc = _corner_r2(width, height, cx, cy, best_ex);
    const double g = 1.0 + best_k[0] * rc + best_k[1] * rc * rc
                         + best_k[2] * rc * rc * rc;
    *out_corner_stops = (g > 1e-6) ? (float)log2(g) : 0.0f;
  }

  return TRUE;
}

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
                              float *out_corner_stops)
{
  if(!lum || width < 32 || height < 32 || stride < width) return FALSE;

  double *r2 = malloc(sizeof(double) * DT_VIG_BINS);
  double *val = malloc(sizeof(double) * DT_VIG_BINS);
  if(!r2 || !val)
  {
    free(r2);
    free(val);
    return FALSE;
  }

  /* Search the ellipticity rather than fit it.
   *
   * It enters the model through r^2 and so is nonlinear, and it is a single
   * well behaved parameter over a known range -- a scan is both simpler than
   * an optimiser and immune to landing in a local minimum. The candidate set
   * brackets the hint, since an anamorphic lens's falloff is stretched by
   * roughly its own squeeze but not exactly. */
  const float seed = (ex_hint > 0.2f && ex_hint < 5.0f) ? ex_hint : 1.0f;

  gboolean any = FALSE;
  float best_k[3] = { 0, 0, 0 };
  double best_rms = 1e30;
  float best_ex = 1.0f;

  for(int s = -8; s <= 8; s++)
  {
    const float ex = seed * (1.0f + 0.06f * s);
    if(ex < 0.2f || ex > 5.0f) continue;

    const int n = _profile(lum, width, height, stride, cx, cy, ex, r2, val);
    if(n < 5) continue;

    const double check = _corner_r2(width, height, cx, cy, ex);

    float k[3] = { 0, 0, 0 };
    double rms = 0.0;
    gboolean got = FALSE;

    // same order of preference as the cell path: physics, then polynomial
    {
      float kt[3];
      double rt;
      if(_fit_cos4(r2, val, NULL, n, check, kt, &rt)
         && _monotonic(kt, check) && _plausible(kt, check))
      {
        memcpy(k, kt, sizeof(k));
        rms = rt;
        got = TRUE;
      }
    }

    for(int nt = 3; nt >= 1 && !got; nt--)
    {
      float kt[3];
      double rt;
      if(!_fit_reciprocal(r2, val, n, nt, kt, &rt)) continue;
      if(!_monotonic(kt, check) || !_plausible(kt, check)) continue;

      memcpy(k, kt, sizeof(k));
      rms = rt;
      got = TRUE;
    }

    if(got && rms < best_rms)
    {
      best_rms = rms;
      memcpy(best_k, k, sizeof(best_k));
      best_ex = ex;
      any = TRUE;
    }
  }

  free(r2);
  free(val);

  if(!any) return FALSE;

  memcpy(out_k, best_k, sizeof(best_k));
  if(out_ex) *out_ex = best_ex;
  if(out_rms) *out_rms = (float)best_rms;

  /* How much correction this amounts to, in stops at the frame corner. The
     coefficients are unreadable; "1.8 stops" is a number that can be
     sanity checked against the picture. */
  if(out_corner_stops)
  {
    const double hd = 0.5 * hypot((double)width, (double)height);
    const double u = (0.5 * width) / hd - cx;
    const double v = (0.5 * height) / hd - cy;
    const double e = u / best_ex;
    const double rc = e * e + v * v;

    const double g = 1.0 + best_k[0] * rc + best_k[1] * rc * rc
                         + best_k[2] * rc * rc * rc;
    *out_corner_stops = (g > 1e-6) ? (float)(log2(g)) : 0.0f;
  }

  return TRUE;
}

// clang-format off
// modelines: These editor modelines have been set for all relevant files by tools/update_modelines.py
// vim: shiftwidth=2 expandtab tabstop=2 cindent
// kate: tab-indents: off; indent-width 2; replace-tabs on; indent-mode cstyle; remove-trailing-spaces modified;
// clang-format on
