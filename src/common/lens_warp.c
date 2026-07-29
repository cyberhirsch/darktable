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

#include "common/lens_warp.h"
#include "common/darktable.h"
#include "common/file_location.h"
#include "iop/gaussian_elimination.h"

#include <glib/gstdio.h>
#include <json-glib/json-glib.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

int dt_lens_warp_poly_terms(const int order)
{
  /* All monomials of total degree 2..order. Degree 0 would be a
     translation, which the optical centre already covers, and degree 1 an
     affine map, which is pose and squeeze rather than distortion -- so
     both are excluded and the model cannot fight the stages that own
     them. */
  const int o = CLAMP(order, 2, DT_LENS_WARP_MAX_ORDER);
  return (o + 1) * (o + 2) / 2 - 3;
}

int dt_lens_warp_param_count(const dt_lens_warp_kind_t kind, const int order)
{
  switch(kind)
  {
    case DT_LENS_WARP_ANAM_RADIAL:
      return 5; // ellipticity plus four radial terms
    case DT_LENS_WARP_RADIAL_POLY:
      return 6; // ellipticity plus c1..c5
    case DT_LENS_WARP_POLY:
    case DT_LENS_WARP_TPS:
    default:
      return 2 * dt_lens_warp_poly_terms(order);
  }
}

const char *dt_lens_warp_kind_name(const dt_lens_warp_kind_t kind)
{
  switch(kind)
  {
    case DT_LENS_WARP_ANAM_RADIAL:   return "anam_radial";
    case DT_LENS_WARP_TPS:           return "tps";
    case DT_LENS_WARP_RADIAL_POLY:   return "radial_poly";
    default:                         return "poly";
  }
}

static dt_lens_warp_kind_t _kind_from_name(const char *s)
{
  if(!g_strcmp0(s, "anam_radial")) return DT_LENS_WARP_ANAM_RADIAL;
  if(!g_strcmp0(s, "tps")) return DT_LENS_WARP_TPS;
  if(!g_strcmp0(s, "radial_poly")) return DT_LENS_WARP_RADIAL_POLY;
  return DT_LENS_WARP_POLY;
}

void dt_lens_warp_init(dt_lens_warp_t *w,
                       const dt_lens_warp_kind_t kind,
                       const int order)
{
  if(!w) return;
  memset(w, 0, sizeof(*w));
  w->kind = kind;
  w->order = CLAMP(order, 2, DT_LENS_WARP_MAX_ORDER);
  w->squeeze = 1.0f;
  // {b, c, v}: no radial term, unit scale -- green's own behaviour
  w->tca_r[2] = 1.0f;
  w->tca_b[2] = 1.0f;
  w->vig_ex = 1.0f;
  w->vig_r_half_diagonal = TRUE;
  w->overscan = 1.0f;
  w->underscan = 1.0f;
  w->nparams = MIN(DT_LENS_WARP_MAX_PARAMS,
                   dt_lens_warp_param_count(kind, w->order));

  /* Both radial kinds divide by an ellipticity, so seed the one parameter
     that is a ratio rather than a coefficient. RADIAL_POLY also needs its
     linear term at one, or the identity warp collapses everything to the
     optical centre. */
  if(kind == DT_LENS_WARP_ANAM_RADIAL) w->p[0] = 1.0f;
  if(kind == DT_LENS_WARP_RADIAL_POLY)
  {
    w->p[0] = 1.0f;
    w->p[1] = 1.0f;
  }
}

void dt_lens_warp_cleanup(dt_lens_warp_t *w)
{
  if(!w) return;
  free(w->tps_src);
  free(w->tps_wx);
  free(w->tps_wy);
  w->tps_src = w->tps_wx = w->tps_wy = NULL;
  w->tps_count = 0;
}

gboolean dt_lens_warp_copy(dt_lens_warp_t *dst, const dt_lens_warp_t *src)
{
  if(!dst || !src) return FALSE;

  memcpy(dst, src, sizeof(*dst));
  dst->tps_src = dst->tps_wx = dst->tps_wy = NULL;

  if(src->tps_count > 0 && src->tps_src && src->tps_wx && src->tps_wy)
  {
    const size_t n = src->tps_count;
    dst->tps_src = malloc(sizeof(float) * 2 * n);
    dst->tps_wx = malloc(sizeof(float) * (n + 3));
    dst->tps_wy = malloc(sizeof(float) * (n + 3));

    if(!dst->tps_src || !dst->tps_wx || !dst->tps_wy)
    {
      dt_lens_warp_cleanup(dst);
      dst->tps_count = 0;
      return FALSE;
    }

    memcpy(dst->tps_src, src->tps_src, sizeof(float) * 2 * n);
    memcpy(dst->tps_wx, src->tps_wx, sizeof(float) * (n + 3));
    memcpy(dst->tps_wy, src->tps_wy, sizeof(float) * (n + 3));
  }
  else
    dst->tps_count = 0;

  return TRUE;
}

/* Evaluate every monomial of degree 2..order once, so the two axes can
   share the work -- they use the same basis with different coefficients. */
static int _poly_basis(float *mono,
                       const float u,
                       const float v,
                       const int order)
{
  float up[DT_LENS_WARP_MAX_ORDER + 1];
  float vp[DT_LENS_WARP_MAX_ORDER + 1];
  up[0] = vp[0] = 1.0f;
  for(int i = 1; i <= order; i++)
  {
    up[i] = up[i - 1] * u;
    vp[i] = vp[i - 1] * v;
  }

  int n = 0;
  for(int d = 2; d <= order; d++)
    for(int i = d; i >= 0; i--)
      mono[n++] = up[i] * vp[d - i];

  return n;
}

/* The thin plate spline kernel. r^2 log r is the fundamental solution of
   the biharmonic equation, which is what makes the interpolant the one
   that bends least. */
static inline float _tps_phi(const float r2)
{
  if(r2 <= 1e-12f) return 0.0f;
  return 0.5f * r2 * logf(r2); // == r^2 * log(r)
}

static void _tps_eval(const dt_lens_warp_t *w,
                      const float u,
                      const float v,
                      float *dx,
                      float *dy)
{
  const int n = w->tps_count;
  double ax = w->tps_wx[n] + w->tps_wx[n + 1] * u + w->tps_wx[n + 2] * v;
  double ay = w->tps_wy[n] + w->tps_wy[n + 1] * u + w->tps_wy[n + 2] * v;

  for(int i = 0; i < n; i++)
  {
    const float du = u - w->tps_src[2 * i];
    const float dv = v - w->tps_src[2 * i + 1];
    const float phi = _tps_phi(du * du + dv * dv);
    ax += (double)w->tps_wx[i] * phi;
    ay += (double)w->tps_wy[i] * phi;
  }

  *dx = (float)ax;
  *dy = (float)ay;
}

/* Everything except the final squeeze. Kept separate because the inverse
   can undo the squeeze exactly and only has to iterate on the rest. */
/* r_observed as a function of r_corrected -- the direction RADIAL_POLY is
   stored in, and the only direction in which it is a polynomial. */
static inline float _radial_poly_fwd(const dt_lens_warp_t *w, const float r)
{
  return r * (w->p[1]
              + r * (w->p[2]
                     + r * (w->p[3]
                            + r * (w->p[4] + r * w->p[5]))));
}

static inline float _radial_poly_deriv(const dt_lens_warp_t *w, const float r)
{
  return w->p[1] + r * (2.0f * w->p[2]
                        + r * (3.0f * w->p[3]
                               + r * (4.0f * w->p[4] + r * 5.0f * w->p[5])));
}

/* Solve _radial_poly_fwd(r_c) = r_o for r_c.
 *
 * Newton rather than the fixed point iteration used for the two dimensional
 * kinds: this is one equation in one unknown, so Newton converges quadratically
 * and to machine precision in a handful of steps -- which matters, because for
 * this kind the *forward* pipeline direction is the one that has to be solved,
 * so the cost is paid on every pixel rather than only when inverting.
 */
static float _radial_poly_solve(const dt_lens_warp_t *w, const float r_o)
{
  if(r_o <= 1e-8f) return r_o;

  float r = r_o;
  if(w->p[1] > 1e-3f) r = r_o / w->p[1]; // first order start

  for(int i = 0; i < 12; i++)
  {
    const float f = _radial_poly_fwd(w, r) - r_o;
    const float d = _radial_poly_deriv(w, r);
    if(fabsf(d) < 1e-9f) break;

    const float step = f / d;
    r -= step;
    if(r < 0.0f) r = 0.5f * (r + step); // never cross the centre
    if(fabsf(step) < 1e-9f) break;
  }

  return r;
}

static void _apply_core(const dt_lens_warp_t *w,
                        const float u,
                        const float v,
                        float *ou,
                        float *ov)
{
  const float du = u - w->cx;
  const float dv = v - w->cy;

  float x = du, y = dv;

  if(w->kind == DT_LENS_WARP_RADIAL_POLY)
  {
    /* Observed to corrected, so the stored polynomial has to be solved. The
       ellipticity scales x before the radius is measured, exactly as in
       ANAM_RADIAL, and Lensfun profiles simply carry 1 there. */
    const float sx = MAX(0.05f, w->p[0]);
    const float px = du / sx;
    const float r_o = sqrtf(px * px + dv * dv);

    if(r_o > 1e-8f)
    {
      const float r_c = _radial_poly_solve(w, r_o);
      const float f = r_c / r_o;
      x = du * f;
      y = dv * f;
    }
  }
  else if(w->kind == DT_LENS_WARP_ANAM_RADIAL)
  {
    /* Radial distortion measured in a space stretched along x by the
       ellipticity p[0]. Note the factor multiplies the *unstretched*
       offset: the ellipticity decides where the distortion is strong, it
       does not itself scale the image. Letting it do both would make it
       indistinguishable from the squeeze stage. */
    const float sx = MAX(0.05f, w->p[0]);
    const float px = du / sx;
    const float r2 = px * px + dv * dv;
    const float f = 1.0f
      + w->p[1] * r2
      + w->p[2] * r2 * r2
      + w->p[3] * r2 * r2 * r2
      + w->p[4] * r2 * r2 * r2 * r2;

    x = du * f;
    y = dv * f;
  }
  else
  {
    float mono[(DT_LENS_WARP_MAX_ORDER + 1) * (DT_LENS_WARP_MAX_ORDER + 2) / 2];
    const int nt = _poly_basis(mono, du, dv, w->order);
    const int have = MIN(nt, w->nparams / 2);

    double ax = 0.0, ay = 0.0;
    for(int i = 0; i < have; i++)
    {
      ax += (double)w->p[i] * mono[i];
      ay += (double)w->p[nt + i] * mono[i];
    }

    x = du + (float)ax;
    y = dv + (float)ay;

    if(w->kind == DT_LENS_WARP_TPS && w->tps_count > 0
       && w->tps_src && w->tps_wx && w->tps_wy)
    {
      /* The spline layer is evaluated at the observed position, not at
         the polynomial's output: it was fitted against observed
         coordinates and has to be asked the same question it answered. */
      float sx = 0.0f, sy = 0.0f;
      _tps_eval(w, du, dv, &sx, &sy);
      x += sx;
      y += sy;
    }
  }

  *ou = x + w->cx;
  *ov = y + w->cy;
}

void dt_lens_warp_apply(const dt_lens_warp_t *w,
                        const float u, const float v,
                        float *ou, float *ov)
{
  if(!w)
  {
    *ou = u;
    *ov = v;
    return;
  }

  float x, y;
  _apply_core(w, u, v, &x, &y);

  const float sq = (w->squeeze > 0.01f) ? w->squeeze : 1.0f;
  *ou = x * sq;
  *ov = y;
}

/* Overscan and underscan, by walking the frame border.
 *
 * Overscan asks how far the recorded frame lands outside the corrected one:
 * push every border point through the warp and see how far past the corrected
 * frame's half extents it reaches.
 *
 * Underscan is the other question and needs the other direction: a corrected
 * pixel is real image only if its inverse lands inside the recorded frame. So
 * shrink a centred rectangle until its whole border inverts to somewhere
 * inside, by bisection. The border is enough for both -- distortion is
 * monotonic outward from the centre, so nothing in the interior of a
 * rectangle can be uncovered while its border is covered.
 */
void dt_lens_warp_measure_scan(dt_lens_warp_t *w,
                               const int width,
                               const int height)
{
  if(!w) return;

  w->overscan = 1.0f;
  w->underscan = 1.0f;

  if(width < 2 || height < 2) return;

  const double hd = 0.5 * hypot((double)width, (double)height);
  const double sq = (w->squeeze > 0.01f) ? (double)w->squeeze : 1.0;

  // half extents of the recorded frame, and of the corrected one
  const double su = 0.5 * width / hd, sv = 0.5 * height / hd;
  const double cu = su * sq, cv = sv;

  const int N = 64;

  double need = 1.0;
  for(int i = 0; i <= N; i++)
  {
    const double t = -1.0 + 2.0 * i / N;

    // the four edges of the recorded frame
    const double pts[4][2] = {
      { t * su, -sv }, { t * su, sv }, { -su, t * sv }, { su, t * sv }
    };

    for(int e = 0; e < 4; e++)
    {
      float ou, ov;
      dt_lens_warp_apply(w, (float)pts[e][0], (float)pts[e][1], &ou, &ov);
      need = MAX(need, MAX(fabs((double)ou) / cu, fabs((double)ov) / cv));
    }
  }
  w->overscan = (float)MIN(need, 4.0);

  double lo = 0.05, hi = 1.0;
  for(int it = 0; it < 24; it++)
  {
    const double s = 0.5 * (lo + hi);
    gboolean covered = TRUE;

    for(int i = 0; i <= N && covered; i++)
    {
      const double t = -1.0 + 2.0 * i / N;
      const double pts[4][2] = {
        { t * s * cu, -s * cv }, { t * s * cu, s * cv },
        { -s * cu, t * s * cv }, { s * cu, t * s * cv }
      };

      for(int e = 0; e < 4 && covered; e++)
      {
        float iu, iv;
        dt_lens_warp_invert(w, (float)pts[e][0], (float)pts[e][1], &iu, &iv);

        // a hair of tolerance, so the exact edge is not called uncovered
        if(fabs((double)iu) > su * 1.0001 || fabs((double)iv) > sv * 1.0001)
          covered = FALSE;
      }
    }

    if(covered)
      lo = s;
    else
      hi = s;
  }
  w->underscan = (float)lo;
}

float dt_lens_tca_scale(const dt_lens_warp_t *w,
                        const int channel,
                        const float r)
{
  if(!w || channel == 1) return 1.0f;

  const float *t = (channel == 0) ? w->tca_r : w->tca_b;
  return t[0] * r * r + t[1] * r + t[2];
}

float dt_lens_vignette_gain(const dt_lens_warp_t *w,
                            const float u, const float v)
{
  if(!w || !w->have_vig) return 1.0f;

  const float ex = (w->vig_ex > 0.05f) ? w->vig_ex : 1.0f;
  const float du = (u - w->cx) / ex;
  const float dv = v - w->cy;
  const float r2 = du * du + dv * dv;

  const float poly = 1.0f + w->vig_k[0] * r2
                          + w->vig_k[1] * r2 * r2
                          + w->vig_k[2] * r2 * r2 * r2;

  /* Which way round the polynomial reads is stored with it rather than
     assumed.
   *
   * For our own measurements it *is* the gain: fitted as v_centre / v(r), so
   * greater than one wherever the frame is darker than its middle, which is
   * exactly the factor a dark pixel needs multiplying by. Returning its
   * reciprocal applies the vignetting a second time instead of undoing it, and
   * looks convincing while doing so -- the corners go dark, which is the
   * direction a lens darkens them, so the picture never says the sign is
   * wrong.
   *
   * For a Lensfun profile it is transmission and divides. 1/t is not a
   * polynomial, so those coefficients can be labelled but not converted, and
   * guessing that a negation would do is the plausible wrong answer: it agrees
   * to first order and parts company exactly where the correction is biggest.
   */
  float gain = poly;
  if(w->vig_convention == DT_LENS_VIG_TRANSMISSION)
    gain = (poly > 1e-4f) ? 1.0f / poly : 16.0f;

  /* Outside the region it was fitted in, the sixth order term can run away or
     go negative. A negative gain would invert the image and an unbounded one
     would turn corner noise into blown highlights, so it is capped at four
     stops either way -- more correction than any real lens needs. */
  if(!(gain > 0.0625f)) return 0.0625f;
  return MIN(gain, 16.0f);
}

gboolean dt_lens_warp_invert(const dt_lens_warp_t *w,
                             const float u, const float v,
                             float *ou, float *ov)
{
  if(!w)
  {
    *ou = u;
    *ov = v;
    return TRUE;
  }

  // the squeeze is a pure scale, so it inverts exactly
  const float sq = (w->squeeze > 0.01f) ? w->squeeze : 1.0f;
  const float tu = u / sq;
  const float tv = v;

  /* RADIAL_POLY is stored in this direction, so here there is nothing to
     iterate -- the polynomial is simply evaluated. The two directions swap
     roles for this kind, which is the price of holding Lensfun's coefficients
     as they were measured rather than as an approximation of their inverse. */
  if(w->kind == DT_LENS_WARP_RADIAL_POLY)
  {
    const float ex = MAX(0.05f, w->p[0]);
    const float du = tu - w->cx;
    const float dv = tv - w->cy;
    const float px = du / ex;
    const float r_c = sqrtf(px * px + dv * dv);

    if(r_c <= 1e-8f)
    {
      *ou = tu;
      *ov = tv;
      return TRUE;
    }

    const float f = _radial_poly_fwd(w, r_c) / r_c;
    *ou = w->cx + du * f;
    *ov = w->cy + dv * f;
    return TRUE;
  }

  /* Fixed point iteration. The map is close enough to the identity for
     this to converge quickly in the middle of the frame; under relaxation
     keeps it stable in the corners of a strongly distorted lens, where a
     full step can overshoot and oscillate. */
  float x = tu, y = tv;
  const float relax = 0.85f;

  for(int i = 0; i < 24; i++)
  {
    float cu, cv;
    _apply_core(w, x, y, &cu, &cv);

    const float eu = tu - cu;
    const float ev = tv - cv;

    if(fabsf(eu) < 1e-7f && fabsf(ev) < 1e-7f)
    {
      *ou = x;
      *ov = y;
      return TRUE;
    }

    x += relax * eu;
    y += relax * ev;
  }

  *ou = x;
  *ov = y;
  return FALSE;
}

gboolean dt_lens_warp_fit_tps(dt_lens_warp_t *w,
                              const float *const src,
                              const float *const dx,
                              const float *const dy,
                              const int count,
                              const float smooth)
{
  if(!w || !src || !dx || !dy || count < 3) return FALSE;

  /* Cap the control point count. The solve is cubic in it, and a chart
     with more nodes than this is measuring the same lens twice rather
     than measuring more of it. */
  const int n = MIN(count, 600);
  const int dim = n + 3;

  double *A = calloc((size_t)dim * dim, sizeof(double));
  double *Acopy = calloc((size_t)dim * dim, sizeof(double));
  double *bx = calloc(dim, sizeof(double));
  double *by = calloc(dim, sizeof(double));

  if(!A || !Acopy || !bx || !by)
  {
    free(A);
    free(Acopy);
    free(bx);
    free(by);
    return FALSE;
  }

  for(int i = 0; i < n; i++)
  {
    for(int j = 0; j < n; j++)
    {
      const float du = src[2 * i] - src[2 * j];
      const float dv = src[2 * i + 1] - src[2 * j + 1];
      A[(size_t)i * dim + j] = _tps_phi(du * du + dv * dv);
    }
    /* The smoothing term on the diagonal turns exact interpolation into a
       regularized fit. Without it the spline threads every measured point
       including its noise, and the correction ripples between nodes. */
    A[(size_t)i * dim + i] += smooth;

    A[(size_t)i * dim + n] = 1.0;
    A[(size_t)i * dim + n + 1] = src[2 * i];
    A[(size_t)i * dim + n + 2] = src[2 * i + 1];

    A[(size_t)n * dim + i] = 1.0;
    A[(size_t)(n + 1) * dim + i] = src[2 * i];
    A[(size_t)(n + 2) * dim + i] = src[2 * i + 1];

    bx[i] = dx[i];
    by[i] = dy[i];
  }

  memcpy(Acopy, A, sizeof(double) * (size_t)dim * dim);

  const int okx = gauss_solve(A, bx, dim);
  const int oky = gauss_solve(Acopy, by, dim);

  free(A);
  free(Acopy);

  if(!okx || !oky)
  {
    free(bx);
    free(by);
    return FALSE;
  }

  free(w->tps_src);
  free(w->tps_wx);
  free(w->tps_wy);

  w->tps_src = malloc(sizeof(float) * 2 * n);
  w->tps_wx = malloc(sizeof(float) * dim);
  w->tps_wy = malloc(sizeof(float) * dim);

  if(!w->tps_src || !w->tps_wx || !w->tps_wy)
  {
    dt_lens_warp_cleanup(w);
    free(bx);
    free(by);
    return FALSE;
  }

  memcpy(w->tps_src, src, sizeof(float) * 2 * n);
  for(int i = 0; i < dim; i++)
  {
    w->tps_wx[i] = (float)bx[i];
    w->tps_wy[i] = (float)by[i];
  }
  w->tps_count = n;

  free(bx);
  free(by);
  return TRUE;
}


/* ---------------------------------------------------------- profiles */

static const char *const _source_names[] =
{
  "measured", "manufacturer", "lensfun", "aggregated",
  "reverse_engineered",
  /* The schema has no "edited" source, on purpose: from the database's point
     of view a hand-altered measurement is not a measurement, and the honest
     thing to call it is the same as anything else somebody typed. It stays a
     distinct value in here so the panel can say which it is. */
  "manufacturer"
};

const char *dt_lens_source_name(const dt_lens_source_t s)
{
  const int n = (int)(sizeof(_source_names) / sizeof(_source_names[0]));
  return _source_names[CLAMP((int)s, 0, n - 1)];
}

dt_lens_source_t dt_lens_source_from_name(const char *name)
{
  if(!name) return DT_LENS_SOURCE_MEASURED;
  for(int i = 0; i <= DT_LENS_SOURCE_REVERSE_ENG; i++)
    if(!strcmp(name, _source_names[i])) return (dt_lens_source_t)i;
  return DT_LENS_SOURCE_MEASURED;
}

void dt_lens_profile_init(dt_lens_profile_t *p)
{
  if(!p) return;
  memset(p, 0, sizeof(*p));
  p->crop_factor = 1.0f;
  p->source = DT_LENS_SOURCE_MEASURED;
  p->warps = g_array_new(FALSE, FALSE, sizeof(dt_lens_warp_t));
}

void dt_lens_profile_cleanup(dt_lens_profile_t *p)
{
  if(!p || !p->warps) return;

  for(guint i = 0; i < p->warps->len; i++)
    dt_lens_warp_cleanup(&g_array_index(p->warps, dt_lens_warp_t, i));

  g_array_free(p->warps, TRUE);
  p->warps = NULL;
}

static gint _warp_compare(gconstpointer a, gconstpointer b)
{
  const dt_lens_warp_t *wa = a;
  const dt_lens_warp_t *wb = b;
  return (wa->focal > wb->focal) - (wa->focal < wb->focal);
}

gboolean dt_lens_profile_add(dt_lens_profile_t *p, const dt_lens_warp_t *w)
{
  if(!p || !p->warps || !w) return FALSE;

  /* A measurement is identified by focal length, aperture *and* focus
     distance, not by focal length alone. Keying on focal alone was fine while
     we only measured distortion, which does not depend on the other two -- but
     vignetting does, and a Lensfun lens carries a dozen or more vignetting
     measurements at one focal length. Replacing on focal alone would keep the
     last of them and discard the rest. */
  for(guint i = 0; i < p->warps->len; i++)
  {
    dt_lens_warp_t *e = &g_array_index(p->warps, dt_lens_warp_t, i);
    if(fabsf(e->focal - w->focal) < 1e-3f
       && fabsf(e->aperture - w->aperture) < 1e-3f
       && fabsf(e->focus_distance - w->focus_distance) < 1e-3f)
    {
      dt_lens_warp_cleanup(e);
      g_array_remove_index(p->warps, i);
      break;
    }
  }

  dt_lens_warp_t copy;
  if(!dt_lens_warp_copy(&copy, w)) return FALSE;

  g_array_append_val(p->warps, copy);
  g_array_sort(p->warps, _warp_compare);
  return TRUE;
}

gboolean dt_lens_profile_eval(const dt_lens_profile_t *p,
                              const float focal,
                              dt_lens_warp_t *out)
{
  return dt_lens_profile_eval_at(p, focal, 0.0f, 0.0f, out);
}

/* Pick the vignetting measurement to use.
 *
 * Separate from the geometry selection above, and deliberately so: distortion
 * does not depend on aperture at all, vignetting depends on it strongly. One
 * lens can carry a dozen vignetting measurements at a single focal length, and
 * choosing between them by focal length would pick whichever happened to be
 * last. Aperture is therefore the primary key here and focal the secondary --
 * the reverse of the geometry.
 *
 * Nothing is interpolated. Two measurements at f/2.8 and f/8 do not average
 * into f/4; the falloff between stops is not linear, and inventing a curve
 * where the data has a gap would be a claim we cannot support. The nearest
 * measurement is used, which at least is a measurement.
 */
static void _select_vignetting(const dt_lens_profile_t *p,
                               const float focal,
                               const float aperture,
                               const float distance,
                               dt_lens_warp_t *out)
{
  const dt_lens_warp_t *best = NULL;
  double best_cost = 1e30;

  for(guint i = 0; i < p->warps->len; i++)
  {
    const dt_lens_warp_t *w = &g_array_index(p->warps, dt_lens_warp_t, i);
    if(!w->have_vig) continue;

    /* Distances in stops and in log focal, so the two are comparable rather
       than one being measured in f-numbers and the other in millimetres. */
    double cost = 0.0;

    if(aperture > 0.0f && w->vig_aperture > 0.0f)
      cost += 4.0 * fabs(log2((double)aperture / w->vig_aperture));
    else if(aperture > 0.0f || w->vig_aperture > 0.0f)
      cost += 1.0; // one side unknown: usable, but not a match

    if(focal > 0.0f && w->focal > 0.0f)
      cost += fabs(log((double)focal / w->focal));

    if(distance > 0.0f && w->focus_distance > 0.0f)
      cost += 0.25 * fabs(log((double)distance / w->focus_distance));

    if(cost < best_cost)
    {
      best_cost = cost;
      best = w;
    }
  }

  if(!best)
  {
    out->have_vig = FALSE;
    return;
  }

  out->have_vig = TRUE;
  memcpy(out->vig_k, best->vig_k, sizeof(out->vig_k));
  out->vig_ex = best->vig_ex;
  out->vig_aperture = best->vig_aperture;
  out->vig_convention = best->vig_convention;
  out->vig_r_half_diagonal = best->vig_r_half_diagonal;
}

gboolean dt_lens_profile_eval_at(const dt_lens_profile_t *p,
                                 const float focal,
                                 const float aperture,
                                 const float distance,
                                 dt_lens_warp_t *out)
{
  if(!p || !p->warps || !p->warps->len || !out) return FALSE;

  /* Only entries that describe geometry take part in the focal length
     interpolation. The array also holds vignetting-only entries -- one per
     aperture -- and those have no distortion to contribute. */
  guint *geo = malloc(sizeof(guint) * p->warps->len);
  if(!geo) return FALSE;

  guint n = 0;
  for(guint i = 0; i < p->warps->len; i++)
    if(g_array_index(p->warps, dt_lens_warp_t, i).have_geometry)
      geo[n++] = i;

  if(!n)
  {
    /* Vignetting without geometry is a legitimate profile: a lens can be
       measured for falloff and never for distortion. */
    dt_lens_warp_init(out, DT_LENS_WARP_ANAM_RADIAL, 4);
    out->focal = focal;
    _select_vignetting(p, focal, aperture, distance, out);
    free(geo);
    return out->have_vig;
  }

#define GEO(k) (&g_array_index(p->warps, dt_lens_warp_t, geo[k]))

  const dt_lens_warp_t *first = GEO(0);
  const dt_lens_warp_t *last = GEO(n - 1);

  // outside the measured range, hold rather than extrapolate
  if(n == 1 || focal <= first->focal || first->focal <= 0.0f)
  {
    const gboolean ok = dt_lens_warp_copy(out, first);
    free(geo);
    if(!ok) return FALSE;
    _select_vignetting(p, focal, aperture, distance, out);
    return TRUE;
  }
  if(focal >= last->focal)
  {
    const gboolean ok = dt_lens_warp_copy(out, last);
    free(geo);
    if(!ok) return FALSE;
    _select_vignetting(p, focal, aperture, distance, out);
    return TRUE;
  }

  guint hi = 1;
  while(hi < n - 1 && GEO(hi)->focal < focal) hi++;

  const dt_lens_warp_t *a = GEO(hi - 1);
  const dt_lens_warp_t *b = GEO(hi);
  free(geo);
#undef GEO

  /* Mixing two different model shapes would produce a third that is
     neither, so fall back to the nearer measurement unless the two are
     structurally the same. */
  const gboolean compatible =
    a->kind == b->kind && a->order == b->order
    && a->nparams == b->nparams
    && a->tps_count == b->tps_count;

  const float la = logf(MAX(1e-3f, a->focal));
  const float lb = logf(MAX(1e-3f, b->focal));
  const float lf = logf(MAX(1e-3f, focal));
  const float t = (lb > la) ? CLAMP((lf - la) / (lb - la), 0.0f, 1.0f) : 0.0f;

  if(!compatible)
  {
    if(!dt_lens_warp_copy(out, (t < 0.5f) ? a : b)) return FALSE;
    _select_vignetting(p, focal, aperture, distance, out);
    return TRUE;
  }

  if(!dt_lens_warp_copy(out, a)) return FALSE;

  out->focal = focal;
  out->cx = a->cx + t * (b->cx - a->cx);
  out->cy = a->cy + t * (b->cy - a->cy);
  out->squeeze = a->squeeze + t * (b->squeeze - a->squeeze);
  for(int i = 0; i < 3; i++)
  {
    out->tca_r[i] = a->tca_r[i] + t * (b->tca_r[i] - a->tca_r[i]);
    out->tca_b[i] = a->tca_b[i] + t * (b->tca_b[i] - a->tca_b[i]);
  }

  /* Interpolating these rather than recomputing them is slightly wrong -- the
     scan factors of a blended warp are not the blend of the scan factors --
     but the error is tiny between neighbouring focal lengths, and it keeps
     the numbers consistent with what the profile promised. Anything needing
     them exactly can call dt_lens_warp_measure_scan on the result. */
  out->overscan = a->overscan + t * (b->overscan - a->overscan);
  out->underscan = a->underscan + t * (b->underscan - a->underscan);

  for(int i = 0; i < out->nparams; i++)
    out->p[i] = a->p[i] + t * (b->p[i] - a->p[i]);

  /* Spline weights only mean the same thing when they sit on the same
     control points; if they do not, keep the nearer layer whole. */
  if(out->tps_count > 0 && a->tps_src && b->tps_src)
  {
    gboolean same_nodes = TRUE;
    for(int i = 0; i < 2 * out->tps_count && same_nodes; i++)
      same_nodes = fabsf(a->tps_src[i] - b->tps_src[i]) < 1e-5f;

    if(same_nodes)
    {
      for(int i = 0; i < out->tps_count + 3; i++)
      {
        out->tps_wx[i] = a->tps_wx[i] + t * (b->tps_wx[i] - a->tps_wx[i]);
        out->tps_wy[i] = a->tps_wy[i] + t * (b->tps_wy[i] - a->tps_wy[i]);
      }
    }
    else if(t >= 0.5f)
    {
      dt_lens_warp_cleanup(out);
      return dt_lens_warp_copy(out, b);
    }
  }

  return TRUE;
}


/* -------------------------------------------------------------- json */

/* Write the v2 document.
 *
 * A warp is decomposed on the way out rather than written whole, because the
 * measurements inside it do not share axes: distortion varies with focal
 * length, vignetting with focal length *and* aperture *and* focus distance.
 * Writing them as one record per focal length would force every vignetting
 * measurement to invent a distortion it never made -- and Lensfun alone holds
 * 23757 vignetting measurements against 5640 distortions, so that is not a
 * rounding error in the file size, it is most of the file.
 *
 * dt_lens_profile_load merges them back, keyed on all three axes.
 */
static void _write_conventions(JsonBuilder *b)
{
  json_builder_set_member_name(b, "conventions");
  json_builder_begin_object(b);

  json_builder_set_member_name(b, "coordinates");
  json_builder_add_string_value(b, "normalized_half_diagonal");
  json_builder_set_member_name(b, "origin");
  json_builder_add_string_value(b, "frame_centre");
  json_builder_set_member_name(b, "angles");
  json_builder_add_string_value(b, "radians");
  json_builder_set_member_name(b, "warp_direction");
  json_builder_add_string_value(b, "observed_to_corrected");
  json_builder_set_member_name(b, "map_direction");
  json_builder_add_string_value(b, "destination_to_source");
  json_builder_set_member_name(b, "vignetting");
  json_builder_add_string_value(b, "correction_gain");

  json_builder_end_object(b);
}

static void _write_axes(JsonBuilder *b, const dt_lens_warp_t *w)
{
  /* Absent means not recorded, which is not a measurement of zero. That
     applies to the focal length as much as to the other two: a manual lens
     reports nothing over the mount, and writing the resulting 0 claimed the
     lens had no focal length rather than that nobody had said what it was.
     A reader takes an absent axis to mean the entry applies at any value of
     it, which for a prime whose distortion was measured once is exactly
     right. */
  if(w->focal > 0.0f)
  {
    json_builder_set_member_name(b, "focal");
    json_builder_add_double_value(b, w->focal);
  }
  if(w->aperture > 0.0f)
  {
    json_builder_set_member_name(b, "aperture");
    json_builder_add_double_value(b, w->aperture);
  }
  if(w->focus_distance > 0.0f)
  {
    json_builder_set_member_name(b, "focus_distance");
    json_builder_add_double_value(b, w->focus_distance);
  }
}

static void _write_distortion(JsonBuilder *b, const dt_lens_warp_t *w)
{
  json_builder_begin_object(b);

  json_builder_set_member_name(b, "model");
  json_builder_add_string_value(b, dt_lens_warp_kind_name(w->kind));
  json_builder_set_member_name(b, "order");
  json_builder_add_int_value(b, w->order);
  _write_axes(b, w);

  json_builder_set_member_name(b, "params");
  json_builder_begin_array(b);
  for(int k = 0; k < w->nparams; k++)
    json_builder_add_double_value(b, w->p[k]);
  json_builder_end_array(b);

  json_builder_set_member_name(b, "overscan");
  json_builder_add_double_value(b, w->overscan);
  json_builder_set_member_name(b, "underscan");
  json_builder_add_double_value(b, w->underscan);

  if(w->tps_count > 0 && w->tps_src && w->tps_wx && w->tps_wy)
  {
    json_builder_set_member_name(b, "tps_src");
    json_builder_begin_array(b);
    for(int k = 0; k < 2 * w->tps_count; k++)
      json_builder_add_double_value(b, w->tps_src[k]);
    json_builder_end_array(b);

    json_builder_set_member_name(b, "tps_wx");
    json_builder_begin_array(b);
    for(int k = 0; k < w->tps_count + 3; k++)
      json_builder_add_double_value(b, w->tps_wx[k]);
    json_builder_end_array(b);

    json_builder_set_member_name(b, "tps_wy");
    json_builder_begin_array(b);
    for(int k = 0; k < w->tps_count + 3; k++)
      json_builder_add_double_value(b, w->tps_wy[k]);
    json_builder_end_array(b);
  }

  json_builder_end_object(b);
}

static gboolean _tca_is_identity(const dt_lens_warp_t *w)
{
  return fabsf(w->tca_r[0]) < 1e-12f && fabsf(w->tca_r[1]) < 1e-12f
      && fabsf(w->tca_r[2] - 1.0f) < 1e-12f
      && fabsf(w->tca_b[0]) < 1e-12f && fabsf(w->tca_b[1]) < 1e-12f
      && fabsf(w->tca_b[2] - 1.0f) < 1e-12f;
}

static void _write_tca(JsonBuilder *b, const dt_lens_warp_t *w)
{
  json_builder_begin_object(b);
  _write_axes(b, w);

  json_builder_set_member_name(b, "tca_r");
  json_builder_begin_array(b);
  for(int k = 0; k < 3; k++) json_builder_add_double_value(b, w->tca_r[k]);
  json_builder_end_array(b);

  json_builder_set_member_name(b, "tca_b");
  json_builder_begin_array(b);
  for(int k = 0; k < 3; k++) json_builder_add_double_value(b, w->tca_b[k]);
  json_builder_end_array(b);

  json_builder_end_object(b);
}

static void _write_vignetting(JsonBuilder *b, const dt_lens_warp_t *w)
{
  json_builder_begin_object(b);

  json_builder_set_member_name(b, "k");
  json_builder_begin_array(b);
  for(int k = 0; k < 3; k++) json_builder_add_double_value(b, w->vig_k[k]);
  json_builder_end_array(b);

  json_builder_set_member_name(b, "ellipticity");
  json_builder_add_double_value(b, w->vig_ex);

  json_builder_set_member_name(b, "focal");
  json_builder_add_double_value(b, w->focal);
  if(w->vig_aperture > 0.0f)
  {
    json_builder_set_member_name(b, "aperture");
    json_builder_add_double_value(b, w->vig_aperture);
  }
  if(w->focus_distance > 0.0f)
  {
    json_builder_set_member_name(b, "focus_distance");
    json_builder_add_double_value(b, w->focus_distance);
  }

  // the reader must not have to guess which way round these read
  json_builder_set_member_name(b, "convention");
  json_builder_add_string_value
    (b, w->vig_convention == DT_LENS_VIG_TRANSMISSION
          ? "transmission_divide" : "gain_multiply");
  json_builder_set_member_name(b, "radius_normalization");
  json_builder_add_string_value
    (b, w->vig_r_half_diagonal ? "half_diagonal" : "calibration_sensor");

  json_builder_end_object(b);
}

gboolean dt_lens_profile_save(const dt_lens_profile_t *p,
                              const char *path,
                              GError **error)
{
  if(!p || !p->warps || !path) return FALSE;

  JsonBuilder *b = json_builder_new();
  json_builder_begin_object(b);

  json_builder_set_member_name(b, "format");
  json_builder_add_string_value(b, "lensfit");
  json_builder_set_member_name(b, "version");
  json_builder_add_int_value(b, 2);

  _write_conventions(b);

  json_builder_set_member_name(b, "imaging_system");
  json_builder_begin_object(b);
  json_builder_set_member_name(b, "id");
  json_builder_add_string_value(b, p->name);
  json_builder_set_member_name(b, "manufacturer");
  json_builder_add_string_value(b, p->maker);
  json_builder_set_member_name(b, "model");
  json_builder_add_string_value(b, p->model);
  json_builder_set_member_name(b, "camera_type");
  json_builder_add_string_value(b, "monocular");
  json_builder_set_member_name(b, "channel_count");
  json_builder_add_int_value(b, 1);
  json_builder_set_member_name(b, "mounts");
  json_builder_begin_array(b);
  if(p->mount[0]) json_builder_add_string_value(b, p->mount);
  json_builder_end_array(b);
  json_builder_set_member_name(b, "crop_factor");
  json_builder_add_double_value(b, p->crop_factor > 0.0f ? p->crop_factor : 1.0);

  /* Absent means unknown, same convention as every other optional axis
     value in this format -- not written at all rather than written as 0,
     which would claim a lens with a zero-millimetre end of its range. */
  if(p->focal_min > 0.0f)
  {
    json_builder_set_member_name(b, "focal_min");
    json_builder_add_double_value(b, p->focal_min);
  }
  if(p->focal_max > 0.0f)
  {
    json_builder_set_member_name(b, "focal_max");
    json_builder_add_double_value(b, p->focal_max);
  }
  if(p->aperture_min > 0.0f)
  {
    json_builder_set_member_name(b, "aperture_min");
    json_builder_add_double_value(b, p->aperture_min);
  }
  if(p->aperture_max > 0.0f)
  {
    json_builder_set_member_name(b, "aperture_max");
    json_builder_add_double_value(b, p->aperture_max);
  }
  if(p->distance_min > 0.0f)
  {
    json_builder_set_member_name(b, "focus_distance_min");
    json_builder_add_double_value(b, p->distance_min);
  }
  if(p->distance_max > 0.0f)
  {
    json_builder_set_member_name(b, "focus_distance_max");
    json_builder_add_double_value(b, p->distance_max);
  }

  json_builder_end_object(b);

  json_builder_set_member_name(b, "source_layout");
  json_builder_begin_object(b);
  json_builder_set_member_name(b, "packing");
  json_builder_add_string_value(b, "single");
  json_builder_set_member_name(b, "reference_width");
  json_builder_add_int_value(b, p->width);
  json_builder_set_member_name(b, "reference_height");
  json_builder_add_int_value(b, p->height);
  json_builder_set_member_name(b, "coordinates_relative_to");
  json_builder_add_string_value(b, "active_area");
  json_builder_end_object(b);

  /* The optical centre and the anisotropy belong to the channel, not to each
     measurement: they describe where the axis meets the sensor, which does not
     change with focal length. Taken from the first geometry entry. */
  const dt_lens_warp_t *ref = NULL;
  for(guint i = 0; i < p->warps->len && !ref; i++)
  {
    const dt_lens_warp_t *w = &g_array_index(p->warps, dt_lens_warp_t, i);
    if(w->have_geometry) ref = w;
  }

  json_builder_set_member_name(b, "channels");
  json_builder_begin_array(b);
  json_builder_begin_object(b);

  json_builder_set_member_name(b, "channel_id");
  json_builder_add_int_value(b, 0);

  json_builder_set_member_name(b, "lens_intrinsics");
  json_builder_begin_object(b);
  json_builder_set_member_name(b, "projection_model");
  json_builder_add_string_value(b, "rectilinear");
  json_builder_set_member_name(b, "cx");
  json_builder_add_double_value(b, ref ? ref->cx : 0.0);
  json_builder_set_member_name(b, "cy");
  json_builder_add_double_value(b, ref ? ref->cy : 0.0);
  json_builder_set_member_name(b, "scale_x");
  json_builder_add_double_value(b, 1.0);
  json_builder_set_member_name(b, "scale_y");
  json_builder_add_double_value(b, ref ? ref->squeeze : 1.0);
  json_builder_set_member_name(b, "skew");
  json_builder_add_double_value(b, 0.0);
  json_builder_set_member_name(b, "radial_normalization");
  json_builder_add_string_value(b, "half_diagonal");
  json_builder_end_object(b);

  json_builder_set_member_name(b, "distortion");
  json_builder_begin_array(b);
  for(guint i = 0; i < p->warps->len; i++)
  {
    const dt_lens_warp_t *w = &g_array_index(p->warps, dt_lens_warp_t, i);
    if(w->have_geometry) _write_distortion(b, w);
  }
  json_builder_end_array(b);

  json_builder_set_member_name(b, "tca");
  json_builder_begin_array(b);
  for(guint i = 0; i < p->warps->len; i++)
  {
    const dt_lens_warp_t *w = &g_array_index(p->warps, dt_lens_warp_t, i);
    if(!_tca_is_identity(w)) _write_tca(b, w);
  }
  json_builder_end_array(b);

  json_builder_set_member_name(b, "photometric");
  json_builder_begin_object(b);
  json_builder_set_member_name(b, "vignetting");
  json_builder_begin_array(b);
  for(guint i = 0; i < p->warps->len; i++)
  {
    const dt_lens_warp_t *w = &g_array_index(p->warps, dt_lens_warp_t, i);
    if(w->have_vig) _write_vignetting(b, w);
  }
  json_builder_end_array(b);
  json_builder_end_object(b);

  json_builder_end_object(b);
  json_builder_end_array(b);

  json_builder_set_member_name(b, "provenance");
  json_builder_begin_object(b);
  json_builder_set_member_name(b, "source");
  json_builder_add_string_value(b, dt_lens_source_name(p->source));
  /* Only an untouched fit gets to claim it was measured. An imported or
     hand-altered profile may be perfectly good and still not be evidence
     of anything the aggregator can count. */
  json_builder_set_member_name(b, "measured");
  json_builder_add_boolean_value(b, p->source == DT_LENS_SOURCE_MEASURED);
  if(p->license[0])
  {
    json_builder_set_member_name(b, "license");
    json_builder_add_string_value(b, p->license);
  }
  if(p->parent[0])
  {
    json_builder_set_member_name(b, "parent_profile");
    json_builder_add_string_value(b, p->parent);
  }
  json_builder_set_member_name(b, "software_version");
  json_builder_add_string_value(b, "darktable lensfit");
  json_builder_end_object(b);

  json_builder_end_object(b);

  JsonGenerator *gen = json_generator_new();
  JsonNode *root = json_builder_get_root(b);
  json_generator_set_root(gen, root);
  json_generator_set_pretty(gen, TRUE);

  gsize len = 0;
  gchar *text = json_generator_to_data(gen, &len);
  const gboolean ok = g_file_set_contents(path, text, len, error);

  g_free(text);
  json_node_free(root);
  g_object_unref(gen);
  g_object_unref(b);

  return ok;
}

static int _read_float_array(JsonObject *obj,
                             const char *member,
                             float *out,
                             const int max)
{
  if(!json_object_has_member(obj, member)) return 0;

  JsonArray *a = json_object_get_array_member(obj, member);
  if(!a) return 0;

  const int n = MIN((int)json_array_get_length(a), max);
  for(int i = 0; i < n; i++)
    out[i] = (float)json_array_get_double_element(a, i);
  return n;
}

static double _member_d(JsonObject *o, const char *name, const double dflt)
{
  return json_object_has_member(o, name)
    ? json_object_get_double_member(o, name) : dflt;
}

/* Find the entry a measurement belongs to, or make one.
 *
 * Keyed on all three axes. Distortion arrives with a focal length and nothing
 * else, vignetting with a focal length, an aperture and a focus distance, so
 * the two land in different entries on purpose -- one lens can hold a dozen
 * vignetting measurements at a single focal length and they are not variants
 * of each other, they are measurements of different conditions.
 *
 * Returns an index rather than a pointer: appending to the array can move it.
 */
static int _v2_slot(GArray *warps,
                    const float focal,
                    const float aperture,
                    const float distance)
{
  for(guint i = 0; i < warps->len; i++)
  {
    const dt_lens_warp_t *w = &g_array_index(warps, dt_lens_warp_t, i);
    if(fabsf(w->focal - focal) < 1e-3f
       && fabsf(w->aperture - aperture) < 1e-3f
       && fabsf(w->focus_distance - distance) < 1e-3f)
      return (int)i;
  }

  dt_lens_warp_t w;
  dt_lens_warp_init(&w, DT_LENS_WARP_ANAM_RADIAL, 4);
  w.focal = focal;
  w.aperture = aperture;
  w.focus_distance = distance;
  g_array_append_val(warps, w);
  return (int)warps->len - 1;
}

static gboolean _load_v2(dt_lens_profile_t *p, JsonObject *o)
{
  JsonObject *sys = json_object_has_member(o, "imaging_system")
    ? json_object_get_object_member(o, "imaging_system") : NULL;

  if(sys)
  {
    if(json_object_has_member(sys, "id"))
      g_strlcpy(p->name, json_object_get_string_member(sys, "id"),
                sizeof(p->name));
    if(json_object_has_member(sys, "manufacturer"))
      g_strlcpy(p->maker, json_object_get_string_member(sys, "manufacturer"),
                sizeof(p->maker));
    if(json_object_has_member(sys, "model"))
    {
      g_strlcpy(p->model, json_object_get_string_member(sys, "model"),
                sizeof(p->model));
      if(!p->name[0]) g_strlcpy(p->name, p->model, sizeof(p->name));
    }
    p->crop_factor = _member_d(sys, "crop_factor", 1.0);
    p->focal_min = _member_d(sys, "focal_min", 0.0);
    p->focal_max = _member_d(sys, "focal_max", 0.0);
    p->aperture_min = _member_d(sys, "aperture_min", 0.0);
    p->aperture_max = _member_d(sys, "aperture_max", 0.0);
    p->distance_min = _member_d(sys, "focus_distance_min", 0.0);
    p->distance_max = _member_d(sys, "focus_distance_max", 0.0);

    if(json_object_has_member(sys, "mounts"))
    {
      JsonArray *m = json_object_get_array_member(sys, "mounts");
      if(m && json_array_get_length(m))
        g_strlcpy(p->mount, json_array_get_string_element(m, 0),
                  sizeof(p->mount));
    }
  }

  if(json_object_has_member(o, "source_layout"))
  {
    JsonObject *sl = json_object_get_object_member(o, "source_layout");
    if(sl)
    {
      p->width = (int)_member_d(sl, "reference_width", 0);
      p->height = (int)_member_d(sl, "reference_height", 0);
    }
  }

  /* Provenance survives a load. Reading an imported profile back and having
     it call itself measured would launder its origin through a round trip,
     which is the one thing this field exists to prevent. */
  if(json_object_has_member(o, "provenance"))
  {
    JsonObject *pv = json_object_get_object_member(o, "provenance");
    if(pv)
    {
      if(json_object_has_member(pv, "source"))
        p->source =
          dt_lens_source_from_name(json_object_get_string_member(pv, "source"));
      if(json_object_has_member(pv, "license"))
        g_strlcpy(p->license, json_object_get_string_member(pv, "license"),
                  sizeof(p->license));
      if(json_object_has_member(pv, "parent_profile"))
        g_strlcpy(p->parent,
                  json_object_get_string_member(pv, "parent_profile"),
                  sizeof(p->parent));
    }
  }

  JsonArray *chs = json_object_has_member(o, "channels")
    ? json_object_get_array_member(o, "channels") : NULL;
  if(!chs || !json_array_get_length(chs)) return FALSE;

  /* One channel for now. A dual fisheye needs the pipeline to sample two
     source regions into one output, which is a different kind of module than
     the one reading this -- so the extra channels are parsed when there is
     something able to use them, rather than half-read now. */
  JsonObject *ch = json_array_get_object_element(chs, 0);
  if(!ch) return FALSE;

  float cx = 0.0f, cy = 0.0f, squeeze = 1.0f;

  if(json_object_has_member(ch, "lens_intrinsics"))
  {
    JsonObject *li = json_object_get_object_member(ch, "lens_intrinsics");
    if(li)
    {
      cx = _member_d(li, "cx", 0.0);
      cy = _member_d(li, "cy", 0.0);

      /* The squeeze is the ratio of the two axis scales, not a field of its
         own: that factorisation gives skew for free and cannot disagree with
         itself the way a separate scalar could. */
      const double sx = _member_d(li, "scale_x", 1.0);
      const double sy = _member_d(li, "scale_y", 1.0);
      if(sx > 1e-6) squeeze = (float)(sy / sx);
    }
  }

  if(json_object_has_member(ch, "distortion"))
  {
    JsonArray *arr = json_object_get_array_member(ch, "distortion");
    for(guint i = 0; arr && i < json_array_get_length(arr); i++)
    {
      JsonObject *d = json_array_get_object_element(arr, i);
      if(!d) continue;

      const int slot = _v2_slot(p->warps,
                                _member_d(d, "focal", 0.0),
                                _member_d(d, "aperture", 0.0),
                                _member_d(d, "focus_distance", 0.0));

      const char *kname = json_object_has_member(d, "model")
        ? json_object_get_string_member(d, "model") : "poly";
      const int order = (int)_member_d(d, "order", 4);

      dt_lens_warp_t tmp;
      dt_lens_warp_init(&tmp, _kind_from_name(kname), order);

      const int np = _read_float_array(d, "params", tmp.p,
                                       DT_LENS_WARP_MAX_PARAMS);
      if(np > 0) tmp.nparams = np;

      tmp.cx = cx;
      tmp.cy = cy;
      tmp.squeeze = squeeze;
      tmp.overscan = _member_d(d, "overscan", 1.0);
      tmp.underscan = _member_d(d, "underscan", 1.0);
      if(!(tmp.overscan >= 1.0f)) tmp.overscan = 1.0f;
      if(!(tmp.underscan > 0.0f) || tmp.underscan > 1.0f) tmp.underscan = 1.0f;

      if(json_object_has_member(d, "tps_src")
         && json_object_has_member(d, "tps_wx")
         && json_object_has_member(d, "tps_wy"))
      {
        JsonArray *sa = json_object_get_array_member(d, "tps_src");
        const int count = sa ? (int)json_array_get_length(sa) / 2 : 0;
        if(count > 0)
        {
          tmp.tps_src = malloc(sizeof(float) * 2 * count);
          tmp.tps_wx = malloc(sizeof(float) * (count + 3));
          tmp.tps_wy = malloc(sizeof(float) * (count + 3));
          if(tmp.tps_src && tmp.tps_wx && tmp.tps_wy)
          {
            _read_float_array(d, "tps_src", tmp.tps_src, 2 * count);
            _read_float_array(d, "tps_wx", tmp.tps_wx, count + 3);
            _read_float_array(d, "tps_wy", tmp.tps_wy, count + 3);
            tmp.tps_count = count;
          }
          else
            dt_lens_warp_cleanup(&tmp);
        }
      }

      /* Carry the axes and the vignetting the slot may already hold, then
         take the slot over -- the geometry is what this entry is for. */
      dt_lens_warp_t *w = &g_array_index(p->warps, dt_lens_warp_t, slot);
      tmp.focal = w->focal;
      tmp.aperture = w->aperture;
      tmp.focus_distance = w->focus_distance;
      tmp.have_vig = w->have_vig;
      memcpy(tmp.vig_k, w->vig_k, sizeof(tmp.vig_k));
      tmp.vig_ex = w->vig_ex;
      tmp.vig_aperture = w->vig_aperture;
      tmp.vig_convention = w->vig_convention;
      tmp.vig_r_half_diagonal = w->vig_r_half_diagonal;
      memcpy(tmp.tca_r, w->tca_r, sizeof(tmp.tca_r));
      memcpy(tmp.tca_b, w->tca_b, sizeof(tmp.tca_b));
      tmp.have_geometry = TRUE;

      dt_lens_warp_cleanup(w);
      *w = tmp;
    }
  }

  if(json_object_has_member(ch, "tca"))
  {
    JsonArray *arr = json_object_get_array_member(ch, "tca");
    for(guint i = 0; arr && i < json_array_get_length(arr); i++)
    {
      JsonObject *t = json_array_get_object_element(arr, i);
      if(!t) continue;

      const int slot = _v2_slot(p->warps,
                                _member_d(t, "focal", 0.0),
                                _member_d(t, "aperture", 0.0),
                                _member_d(t, "focus_distance", 0.0));
      dt_lens_warp_t *w = &g_array_index(p->warps, dt_lens_warp_t, slot);

      _read_float_array(t, "tca_r", w->tca_r, 3);
      _read_float_array(t, "tca_b", w->tca_b, 3);
    }
  }

  if(json_object_has_member(ch, "photometric"))
  {
    JsonObject *ph = json_object_get_object_member(ch, "photometric");
    JsonArray *arr = (ph && json_object_has_member(ph, "vignetting"))
      ? json_object_get_array_member(ph, "vignetting") : NULL;

    for(guint i = 0; arr && i < json_array_get_length(arr); i++)
    {
      JsonObject *v = json_array_get_object_element(arr, i);
      if(!v || !json_object_has_member(v, "k")) continue;

      const float ap = _member_d(v, "aperture", 0.0);
      const int slot = _v2_slot(p->warps,
                                _member_d(v, "focal", 0.0),
                                ap,
                                _member_d(v, "focus_distance", 0.0));
      dt_lens_warp_t *w = &g_array_index(p->warps, dt_lens_warp_t, slot);

      _read_float_array(v, "k", w->vig_k, 3);
      w->have_vig = TRUE;
      w->vig_ex = _member_d(v, "ellipticity", 1.0);
      if(!(w->vig_ex > 0.05f)) w->vig_ex = 1.0f;
      w->vig_aperture = ap;
      w->cx = cx;
      w->cy = cy;

      /* Convention and normalisation are read, never assumed. A gain and a
         transmission are reciprocals of each other and there is no way to
         tell them apart from the coefficients alone -- the sign that looks
         like a giveaway is the one that misled us once already. */
      if(json_object_has_member(v, "convention"))
        w->vig_convention =
          !g_strcmp0(json_object_get_string_member(v, "convention"),
                     "transmission_divide")
            ? DT_LENS_VIG_TRANSMISSION : DT_LENS_VIG_GAIN;

      if(json_object_has_member(v, "radius_normalization"))
        w->vig_r_half_diagonal =
          g_strcmp0(json_object_get_string_member(v, "radius_normalization"),
                    "calibration_sensor") != 0;
    }
  }

  if(!p->warps->len) return FALSE;

  g_array_sort(p->warps, _warp_compare);
  return TRUE;
}

gboolean dt_lens_profile_load(dt_lens_profile_t *p, const char *path)
{
  if(!p || !path) return FALSE;

  /* Initialise before anything can fail, so that every return path leaves
     `p` in a state the caller can safely clean up. */
  dt_lens_profile_init(p);

  JsonParser *parser = json_parser_new();
  GError *err = NULL;

  if(!json_parser_load_from_file(parser, path, &err))
  {
    dt_print(DT_DEBUG_ALWAYS, "[lens_warp] cannot read `%s': %s",
             path, err ? err->message : "?");
    if(err) g_error_free(err);
    g_object_unref(parser);
    return FALSE;
  }

  JsonNode *root = json_parser_get_root(parser);
  JsonObject *o = root ? json_node_get_object(root) : NULL;
  if(!o)
  {
    g_object_unref(parser);
    return FALSE;
  }

  /* Version 2 nests the measurements under channels and splits them by kind;
     version 1 is a flat list of warps. Both are read, because profiles saved
     before the format changed are somebody's afternoon of work. */
  const int ver = json_object_has_member(o, "version")
    ? (int)json_object_get_int_member(o, "version") : 1;

  if(ver >= 2)
  {
    const gboolean ok = _load_v2(p, o);
    g_object_unref(parser);
    dt_print(DT_DEBUG_ALWAYS, "[lens_warp] loaded `%s' (v2): %u entries",
             path, p->warps->len);
    return ok;
  }

  if(json_object_has_member(o, "name"))
    g_strlcpy(p->name, json_object_get_string_member(o, "name"),
              sizeof(p->name));
  if(json_object_has_member(o, "maker"))
    g_strlcpy(p->maker, json_object_get_string_member(o, "maker"),
              sizeof(p->maker));
  if(json_object_has_member(o, "model"))
    g_strlcpy(p->model, json_object_get_string_member(o, "model"),
              sizeof(p->model));
  if(json_object_has_member(o, "mount"))
    g_strlcpy(p->mount, json_object_get_string_member(o, "mount"),
              sizeof(p->mount));
  if(json_object_has_member(o, "width"))
    p->width = json_object_get_int_member(o, "width");
  if(json_object_has_member(o, "height"))
    p->height = json_object_get_int_member(o, "height");
  if(json_object_has_member(o, "crop_factor"))
    p->crop_factor = json_object_get_double_member(o, "crop_factor");

  JsonArray *warps = json_object_has_member(o, "warps")
    ? json_object_get_array_member(o, "warps") : NULL;

  const guint nw = warps ? json_array_get_length(warps) : 0;
  for(guint i = 0; i < nw; i++)
  {
    JsonObject *wo = json_array_get_object_element(warps, i);
    if(!wo) continue;

    const char *kname = json_object_has_member(wo, "kind")
      ? json_object_get_string_member(wo, "kind") : "poly";
    const int order = json_object_has_member(wo, "order")
      ? json_object_get_int_member(wo, "order") : 4;

    dt_lens_warp_t w;
    dt_lens_warp_init(&w, _kind_from_name(kname), order);

    if(json_object_has_member(wo, "focal"))
      w.focal = json_object_get_double_member(wo, "focal");
    if(json_object_has_member(wo, "cx"))
      w.cx = json_object_get_double_member(wo, "cx");
    if(json_object_has_member(wo, "cy"))
      w.cy = json_object_get_double_member(wo, "cy");
    if(json_object_has_member(wo, "squeeze"))
      w.squeeze = json_object_get_double_member(wo, "squeeze");
    _read_float_array(wo, "tca_r", w.tca_r, 3);
    _read_float_array(wo, "tca_b", w.tca_b, 3);
    if(json_object_has_member(wo, "aperture"))
      w.aperture = json_object_get_double_member(wo, "aperture");
    if(json_object_has_member(wo, "focus_distance"))
      w.focus_distance = json_object_get_double_member(wo, "focus_distance");
    if(json_object_has_member(wo, "overscan"))
      w.overscan = json_object_get_double_member(wo, "overscan");
    if(json_object_has_member(wo, "underscan"))
      w.underscan = json_object_get_double_member(wo, "underscan");

    // an older profile predates these; neutral is better than zero
    if(!(w.overscan >= 1.0f)) w.overscan = 1.0f;
    if(!(w.underscan > 0.0f) || w.underscan > 1.0f) w.underscan = 1.0f;

    const int np = _read_float_array(wo, "params", w.p,
                                     DT_LENS_WARP_MAX_PARAMS);
    if(np > 0) w.nparams = np;

    if(json_object_has_member(wo, "vignette"))
    {
      JsonObject *vo = json_object_get_object_member(wo, "vignette");
      if(vo && json_object_has_member(vo, "k"))
      {
        _read_float_array(vo, "k", w.vig_k, 3);
        w.have_vig = TRUE;

        if(json_object_has_member(vo, "ellipticity"))
          w.vig_ex = json_object_get_double_member(vo, "ellipticity");
        if(json_object_has_member(vo, "aperture"))
          w.vig_aperture = json_object_get_double_member(vo, "aperture");
        if(!(w.vig_ex > 0.05f)) w.vig_ex = 1.0f;

        if(json_object_has_member(vo, "convention"))
          w.vig_convention =
            !g_strcmp0(json_object_get_string_member(vo, "convention"),
                       "transmission_divide")
              ? DT_LENS_VIG_TRANSMISSION : DT_LENS_VIG_GAIN;

        if(json_object_has_member(vo, "radius_normalization"))
          w.vig_r_half_diagonal =
            g_strcmp0(json_object_get_string_member(vo, "radius_normalization"),
                      "calibration_sensor") != 0;
      }
    }

    if(json_object_has_member(wo, "tps_src")
       && json_object_has_member(wo, "tps_wx")
       && json_object_has_member(wo, "tps_wy"))
    {
      JsonArray *sa = json_object_get_array_member(wo, "tps_src");
      const int count = sa ? (int)json_array_get_length(sa) / 2 : 0;

      if(count > 0)
      {
        w.tps_src = malloc(sizeof(float) * 2 * count);
        w.tps_wx = malloc(sizeof(float) * (count + 3));
        w.tps_wy = malloc(sizeof(float) * (count + 3));

        if(w.tps_src && w.tps_wx && w.tps_wy)
        {
          _read_float_array(wo, "tps_src", w.tps_src, 2 * count);
          _read_float_array(wo, "tps_wx", w.tps_wx, count + 3);
          _read_float_array(wo, "tps_wy", w.tps_wy, count + 3);
          w.tps_count = count;
        }
        else
          dt_lens_warp_cleanup(&w);
      }
    }

    /* Every v1 warp carried geometry -- the format had no way to express a
       measurement that did not -- so the flag is set rather than inferred. */
    w.have_geometry = TRUE;

    // add copies, then release our local one
    dt_lens_profile_add(p, &w);
    dt_lens_warp_cleanup(&w);
  }

  g_object_unref(parser);

  dt_print(DT_DEBUG_ALWAYS, "[lens_warp] loaded `%s': %u warp(s)",
           path, p->warps->len);

  return p->warps->len > 0;
}

gchar *dt_lens_profile_dir(void)
{
  char confdir[PATH_MAX] = { 0 };
  dt_loc_get_user_config_dir(confdir, sizeof(confdir));

  gchar *dir = g_build_filename(confdir, "lensfit", NULL);
  g_mkdir_with_parents(dir, 0755);
  return dir;
}

// the shipped, read-only database installed alongside darktable itself.
// Never created here -- if the package did not install one, this simply
// will not exist and every lookup through it falls straight through.
static gchar *_lens_profile_shared_dir(void)
{
  char datadir[PATH_MAX] = { 0 };
  dt_loc_get_datadir(datadir, sizeof(datadir));
  return g_build_filename(datadir, "lensfit", NULL);
}

// g_ptr_array_sort hands the comparator pointers *to* the elements
static gint _name_compare(gconstpointer a, gconstpointer b)
{
  const gchar *const *sa = a;
  const gchar *const *sb = b;
  return g_ascii_strcasecmp(*sa, *sb);
}

static void _collect_profile_names(const char *dir, GHashTable *seen,
                                    GPtrArray *out)
{
  GDir *d = g_dir_open(dir, 0, NULL);
  if(!d) return;

  const gchar *entry;
  while((entry = g_dir_read_name(d)))
  {
    if(!g_str_has_suffix(entry, ".json")) continue;
    gchar *base = g_strndup(entry, strlen(entry) - 5);
    if(g_hash_table_contains(seen, base))
    {
      g_free(base);
      continue;
    }
    g_hash_table_add(seen, g_strdup(base));
    g_ptr_array_add(out, base);
  }
  g_dir_close(d);
}

gchar **dt_lens_profile_list(void)
{
  GPtrArray *out = g_ptr_array_new();
  GHashTable *seen = g_hash_table_new_full(g_str_hash, g_str_equal, g_free,
                                            NULL);

  // user profiles first so a name clash is deduped in the user's favour --
  // _collect_profile_names skips a name already seen
  gchar *user_dir = dt_lens_profile_dir();
  _collect_profile_names(user_dir, seen, out);
  g_free(user_dir);

  gchar *shared_dir = _lens_profile_shared_dir();
  _collect_profile_names(shared_dir, seen, out);
  g_free(shared_dir);

  g_hash_table_destroy(seen);

  g_ptr_array_sort(out, _name_compare);

  g_ptr_array_add(out, NULL);
  return (gchar **)g_ptr_array_free(out, FALSE);
}

gchar *dt_lens_profile_path(const char *name)
{
  if(!name || !*name) return NULL;

  gchar *dir = dt_lens_profile_dir();
  gchar *file = g_strdup_printf("%s.json", name);
  gchar *path = g_build_filename(dir, file, NULL);

  g_free(dir);
  g_free(file);
  return path;
}

gchar *dt_lens_profile_find(const char *name)
{
  if(!name || !*name) return NULL;

  gchar *file = g_strdup_printf("%s.json", name);

  gchar *user_dir = dt_lens_profile_dir();
  gchar *user_path = g_build_filename(user_dir, file, NULL);
  g_free(user_dir);
  if(g_file_test(user_path, G_FILE_TEST_IS_REGULAR))
  {
    g_free(file);
    return user_path;
  }
  g_free(user_path);

  gchar *shared_dir = _lens_profile_shared_dir();
  gchar *shared_path = g_build_filename(shared_dir, file, NULL);
  g_free(shared_dir);
  g_free(file);
  if(g_file_test(shared_path, G_FILE_TEST_IS_REGULAR))
    return shared_path;

  g_free(shared_path);
  return NULL;
}

/* -------------------------------------------- automatic profile matching */

/* Reproduces tools/lensfun_to_lensfit.py's safe_name() exactly, so a name
 * generated here for Lensfun's maker/model reproduces the same profile
 * filename the converter chose -- verified against the shipped database:
 * 1289 of 1290 filenames regenerate exactly this way (the one miss is a
 * case-only collision already carrying a disambiguating suffix, which
 * falls through to the fuzzy matcher below instead). */
static void _manifest_safe_name(const char *maker, const char *model,
                                char *out, size_t out_size)
{
  gchar *cat = g_strdup_printf("%s_%s", maker ? maker : "", model ? model : "");

  gchar *start = cat;
  while(*start == '_') start++;
  gchar *end = start + strlen(start);
  while(end > start && *(end - 1) == '_') end--;
  gchar *stripped = g_strndup(start, end - start);
  g_free(cat);

  for(gchar *p = stripped; *p; p++)
    if(strchr("/\\:*?\"<>|", *p)) *p = '_';

  gchar **parts = g_strsplit_set(stripped, " \t\n\r\f\v", -1);
  GString *joined = g_string_new(NULL);
  for(gchar **p = parts; *p; p++)
  {
    if(!**p) continue;
    if(joined->len) g_string_append_c(joined, ' ');
    g_string_append(joined, *p);
  }
  g_strfreev(parts);
  g_free(stripped);

  for(gchar *p = joined->str; *p; p++)
    if(*p == ' ') *p = '_';

  g_strlcpy(out, joined->len ? joined->str : "lens", out_size);
  g_string_free(joined, TRUE);
}

typedef struct _manifest_entry_t
{
  gchar *maker, *model;
} _manifest_entry_t;

static void _manifest_entry_free(gpointer p)
{
  _manifest_entry_t *e = p;
  g_free(e->maker);
  g_free(e->model);
  g_free(e);
}

// name -> _manifest_entry_t, lazily loaded from <datadir>/lensfit/index.json
// and kept for the process lifetime -- 1290 short entries, not worth
// reloading per lookup, and the shipped database does not change at runtime
static GHashTable *_manifest = NULL;
static gboolean _manifest_load_attempted = FALSE;

static GHashTable *_get_manifest(void)
{
  if(_manifest_load_attempted) return _manifest;
  _manifest_load_attempted = TRUE;

  char datadir[PATH_MAX] = { 0 };
  dt_loc_get_datadir(datadir, sizeof(datadir));
  gchar *path = g_build_filename(datadir, "lensfit", "index.json", NULL);

  JsonParser *parser = json_parser_new();
  if(!json_parser_load_from_file(parser, path, NULL))
  {
    g_object_unref(parser);
    g_free(path);
    return NULL; // no shipped database installed -- not an error
  }
  g_free(path);

  JsonNode *root_node = json_parser_get_root(parser);
  JsonObject *root = root_node ? json_node_get_object(root_node) : NULL;
  JsonObject *profiles = root && json_object_has_member(root, "profiles")
    ? json_object_get_object_member(root, "profiles") : NULL;

  if(profiles)
  {
    GHashTable *m = g_hash_table_new_full(g_str_hash, g_str_equal,
                                          g_free, _manifest_entry_free);
    GList *names = json_object_get_members(profiles);
    for(GList *l = names; l; l = l->next)
    {
      const char *name = l->data;
      JsonObject *e = json_object_get_object_member(profiles, name);
      if(!e) continue;

      _manifest_entry_t *me = g_new0(_manifest_entry_t, 1);
      me->maker = g_strdup(json_object_has_member(e, "maker")
                           ? json_object_get_string_member(e, "maker") : "");
      me->model = g_strdup(json_object_has_member(e, "model")
                           ? json_object_get_string_member(e, "model") : "");
      g_hash_table_insert(m, g_strdup(name), me);
    }
    g_list_free(names);
    _manifest = m;
  }

  g_object_unref(parser);
  return _manifest;
}

typedef struct { gchar *name, *maker, *model; } _profile_row_t;

static int _profile_row_cmp(const void *a, const void *b)
{
  const _profile_row_t *ra = a, *rb = b;
  int c = g_ascii_strcasecmp(ra->maker, rb->maker);
  if(c) return c;
  return g_ascii_strcasecmp(ra->model, rb->model);
}

void dt_lens_profile_list_full(gchar ***out_names,
                               gchar ***out_makers,
                               gchar ***out_models)
{
  GPtrArray *names = g_ptr_array_new();
  GPtrArray *makers = g_ptr_array_new();
  GPtrArray *models = g_ptr_array_new();
  GHashTable *seen = g_hash_table_new_full(g_str_hash, g_str_equal, g_free,
                                            NULL);

  // user profiles first, same precedence as dt_lens_profile_list() -- a
  // name clash is deduped in the user's favour. There are typically only a
  // handful of these, so loading each file for its maker/model is cheap.
  gchar *user_dir = dt_lens_profile_dir();
  GDir *ud = g_dir_open(user_dir, 0, NULL);
  if(ud)
  {
    const gchar *entry;
    while((entry = g_dir_read_name(ud)))
    {
      if(!g_str_has_suffix(entry, ".json")) continue;
      gchar *base = g_strndup(entry, strlen(entry) - 5);
      if(g_hash_table_contains(seen, base)) { g_free(base); continue; }
      g_hash_table_add(seen, g_strdup(base));

      gchar *path = g_build_filename(user_dir, entry, NULL);
      dt_lens_profile_t prof;
      const gboolean ok = dt_lens_profile_load(&prof, path);
      g_free(path);

      g_ptr_array_add(names, base);
      g_ptr_array_add(makers, g_strdup((ok && prof.maker[0]) ? prof.maker : ""));
      g_ptr_array_add(models, g_strdup((ok && prof.model[0]) ? prof.model : base));
      if(ok) dt_lens_profile_cleanup(&prof);
    }
    g_dir_close(ud);
  }
  g_free(user_dir);

  // shared, read-only database -- maker/model come straight out of the
  // already-loaded manifest, so this stays a hash lookup per entry rather
  // than opening 1000+ files.
  GHashTable *manifest = _get_manifest();
  gchar *shared_dir = _lens_profile_shared_dir();
  GDir *sd = g_dir_open(shared_dir, 0, NULL);
  if(sd)
  {
    const gchar *entry;
    while((entry = g_dir_read_name(sd)))
    {
      if(!g_str_has_suffix(entry, ".json")) continue;
      gchar *base = g_strndup(entry, strlen(entry) - 5);
      if(g_hash_table_contains(seen, base)) { g_free(base); continue; }
      g_hash_table_add(seen, g_strdup(base));

      _manifest_entry_t *me = manifest
        ? g_hash_table_lookup(manifest, base) : NULL;

      g_ptr_array_add(names, base);
      g_ptr_array_add(makers, g_strdup(me ? me->maker : ""));
      g_ptr_array_add(models, g_strdup((me && me->model[0]) ? me->model : base));
    }
    g_dir_close(sd);
  }
  g_free(shared_dir);
  g_hash_table_destroy(seen);

  // sort by (maker, model) together -- each row is self-contained so a
  // plain qsort comparator needs no external state
  const guint n = names->len;
  _profile_row_t *rows = g_new(_profile_row_t, n);
  for(guint i = 0; i < n; i++)
  {
    rows[i].name = g_ptr_array_index(names, i);
    rows[i].maker = g_ptr_array_index(makers, i);
    rows[i].model = g_ptr_array_index(models, i);
  }
  g_ptr_array_free(names, FALSE);
  g_ptr_array_free(makers, FALSE);
  g_ptr_array_free(models, FALSE);

  qsort(rows, n, sizeof(_profile_row_t), _profile_row_cmp);

  gchar **on = g_new0(gchar *, n + 1);
  gchar **om = g_new0(gchar *, n + 1);
  gchar **omo = g_new0(gchar *, n + 1);
  for(guint i = 0; i < n; i++)
  {
    on[i] = rows[i].name;
    om[i] = rows[i].maker;
    omo[i] = rows[i].model;
  }
  g_free(rows);

  *out_names = on;
  *out_makers = om;
  *out_models = omo;
}

// keeps only letters and digits, casefolded -- raw EXIF lens names and
// Lensfun/the manifest's own spelling of the same lens routinely differ in
// punctuation only, e.g. a camera's EXIF "F4-5.6" against Lensfun's own
// "f/4-5.6": same lens, different aperture-range notation. Comparing on
// alnum-only content instead of the literal string sidesteps every such
// formatting difference (slashes, hyphens vs. en-dashes, doubled spaces)
// at once rather than chasing them one at a time.
static gchar *_alnum_casefold(const char *s)
{
  gchar *cf = g_utf8_casefold(s, -1);
  GString *out = g_string_sized_new(strlen(cf));
  for(const gchar *p = cf; *p; p = g_utf8_next_char(p))
  {
    const gunichar c = g_utf8_get_char(p);
    if(g_unichar_isalnum(c)) g_string_append_unichar(out, c);
  }
  g_free(cf);
  return g_string_free(out, FALSE);
}

// substring in either direction on alnum-only content -- deliberately loose
// on wording and punctuation but still requires real overlap, not merely
// "some characters in common"
static gboolean _fuzzy_hit(const char *a, const char *b)
{
  if(!a || !b || !*a || !*b) return FALSE;
  gchar *ca = _alnum_casefold(a);
  gchar *cb = _alnum_casefold(b);
  const gboolean hit = *ca && *cb
    && (strstr(ca, cb) != NULL || strstr(cb, ca) != NULL);
  g_free(ca);
  g_free(cb);
  return hit;
}

gboolean dt_lens_profile_match_user(const char *maker, const char *model,
                                    char *out_name, size_t out_size)
{
  if(!model || !*model) return FALSE;

  gchar *dir = dt_lens_profile_dir();
  GDir *d = g_dir_open(dir, 0, NULL);
  gboolean found = FALSE;

  if(d)
  {
    const gchar *entry;
    while(!found && (entry = g_dir_read_name(d)))
    {
      if(!g_str_has_suffix(entry, ".json")) continue;

      gchar *path = g_build_filename(dir, entry, NULL);
      dt_lens_profile_t prof;
      if(dt_lens_profile_load(&prof, path)
         && _fuzzy_hit(model, prof.model)
         && (!maker || !*maker || !prof.maker[0] || _fuzzy_hit(maker, prof.maker)))
      {
        gchar *base = g_strndup(entry, strlen(entry) - 5); // strip ".json"
        g_strlcpy(out_name, base, out_size);
        g_free(base);
        found = TRUE;
      }
      dt_lens_profile_cleanup(&prof);
      g_free(path);
    }
    g_dir_close(d);
  }

  g_free(dir);
  return found;
}

gboolean dt_lens_profile_manifest_lookup(const char *maker, const char *model,
                                         char *out_name, size_t out_size)
{
  GHashTable *m = _get_manifest();
  if(!m || !maker || !model || !*model) return FALSE;

  char name[256];
  _manifest_safe_name(maker, model, name, sizeof(name));
  if(g_hash_table_contains(m, name))
  {
    g_strlcpy(out_name, name, out_size);
    return TRUE;
  }

  /* The exact form is case-sensitive by construction -- it has to
     reproduce the converter's generated name exactly -- but `maker` is
     not guaranteed to be spelled the way Lensfun itself spells it (this
     already happened once: "SONY" from a camera's EXIF Make tag against
     Lensfun's own "Sony"). Falls back to a case-insensitive scan of the
     generated names, still an *exact* match on content, just not on case,
     before giving up and letting the caller try the fuzzy matcher. */
  gchar *casefold_name = g_utf8_casefold(name, -1);
  GHashTableIter it;
  gpointer key;
  gboolean found = FALSE;
  g_hash_table_iter_init(&it, m);
  while(g_hash_table_iter_next(&it, &key, NULL))
  {
    gchar *ck = g_utf8_casefold((const char *)key, -1);
    const gboolean hit = !strcmp(ck, casefold_name);
    g_free(ck);
    if(hit)
    {
      g_strlcpy(out_name, (const char *)key, out_size);
      found = TRUE;
      break;
    }
  }
  g_free(casefold_name);
  return found;
}

gboolean dt_lens_profile_manifest_match(const char *maker, const char *model,
                                        char *out_name, size_t out_size)
{
  GHashTable *m = _get_manifest();
  // both ends required, and deliberately so: model text alone is where a
  // plausible-looking wrong match would come from ("50mm" matches
  // everything), and this path exists exactly to avoid that
  if(!m || !maker || !*maker || !model || !*model) return FALSE;

  GHashTableIter it;
  gpointer key, val;
  g_hash_table_iter_init(&it, m);
  while(g_hash_table_iter_next(&it, &key, &val))
  {
    const _manifest_entry_t *e = val;
    if(_fuzzy_hit(model, e->model) && _fuzzy_hit(maker, e->maker))
    {
      g_strlcpy(out_name, (const char *)key, out_size);
      return TRUE;
    }
  }
  return FALSE;
}

// clang-format off
// modelines: These editor modelines have been set for all relevant files by tools/update_modelines.py
// vim: shiftwidth=2 expandtab tabstop=2 cindent
// kate: tab-indents: off; indent-width 2; replace-tabs on; indent-mode cstyle; remove-trailing-spaces modified;
// clang-format on
