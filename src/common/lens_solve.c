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

#include "common/lens_solve.h"
#include "common/darktable.h"
#include "iop/gaussian_elimination.h"

#include <math.h>
#include <stdlib.h>
#include <string.h>

#define DT_LENS_SOLVE_MAX_PARAMS (DT_LENS_WARP_MAX_PARAMS + 2)

// a line needs three points before it can be anything but straight
#define DT_LENS_SOLVE_MIN_LINE 3


typedef struct _solve_ctx_t
{
  int n;
  float *u, *v;     // observed, normalized
  float *wu, *wv;   // corrected, normalized -- rewritten on every evaluation

  // groups of point indices, one per chart row and per chart column
  int ngroups;
  int *gstart;      // ngroups + 1 offsets into gidx
  int *gidx;

  int nres_line;
  int nres;

  /* The ideal lattice, and scratch for the corrected positions in the
     form _fit_affine wants them. Only used when the chart is declared
     square on. */
  gboolean gauge;
  double *ix, *iy;
  double *dtx, *dty;

  /* Weight on the gauge block. It has to be small.
   *
   * The gauge contributes two residuals per point against one for the line
   * term, and each is far larger -- a line residual is the sub-pixel wobble
   * left in a chart line, while a gauge residual is the whole distortion,
   * which is what we are here to measure and cannot be removed by an affine.
   * Left at equal weight the gauge outweighs straightness by a factor of
   * fifty, and the solver duly sacrifices straightness to it: the fit comes
   * back with the lines *less* straight than it found them. The gauge is
   * meant to choose between answers that fit equally well, so it is scaled
   * to sit an order of magnitude below the objective it is breaking ties in.
   */
  double gauge_w;

  int npar;
  int nwarp;        // how many of the parameters belong to the warp itself

  float regw[DT_LENS_WARP_MAX_PARAMS];
  float regt[DT_LENS_WARP_MAX_PARAMS]; // the value each is pulled towards

  dt_lens_warp_t work;
  gboolean solve_centre;
} _solve_ctx_t;


/* Bucket the points by lattice row and by lattice column. Groups with too
   few members are dropped rather than kept and ignored, so the residual
   count stays fixed for the whole solve -- Levenberg-Marquardt compares
   costs across iterations and cannot do that if the vector changes length. */
static gboolean _build_groups(_solve_ctx_t *c,
                              const dt_lens_solve_point_t *pts,
                              const int n)
{
  int rmin = G_MAXINT, rmax = G_MININT, cmin = G_MAXINT, cmax = G_MININT;
  for(int i = 0; i < n; i++)
  {
    rmin = MIN(rmin, pts[i].row);
    rmax = MAX(rmax, pts[i].row);
    cmin = MIN(cmin, pts[i].col);
    cmax = MAX(cmax, pts[i].col);
  }
  if(rmax < rmin || cmax < cmin) return FALSE;

  const int nrows = rmax - rmin + 1;
  const int ncols = cmax - cmin + 1;
  const int nbuck = nrows + ncols;

  int *count = calloc(nbuck, sizeof(int));
  if(!count) return FALSE;

  for(int i = 0; i < n; i++)
  {
    count[pts[i].row - rmin]++;
    count[nrows + pts[i].col - cmin]++;
  }

  // keep only the buckets that constrain anything
  int *map = calloc(nbuck, sizeof(int));
  if(!map)
  {
    free(count);
    return FALSE;
  }

  int kept = 0, total = 0;
  for(int b = 0; b < nbuck; b++)
  {
    if(count[b] >= DT_LENS_SOLVE_MIN_LINE)
    {
      map[b] = kept++;
      total += count[b];
    }
    else
      map[b] = -1;
  }

  if(!kept)
  {
    free(count);
    free(map);
    return FALSE;
  }

  c->ngroups = kept;
  c->gstart = calloc(kept + 1, sizeof(int));
  c->gidx = calloc(MAX(1, total), sizeof(int));
  if(!c->gstart || !c->gidx)
  {
    free(count);
    free(map);
    return FALSE;
  }

  for(int b = 0; b < nbuck; b++)
    if(map[b] >= 0) c->gstart[map[b] + 1] = count[b];
  for(int g = 0; g < kept; g++) c->gstart[g + 1] += c->gstart[g];

  int *fill = calloc(kept, sizeof(int));
  if(!fill)
  {
    free(count);
    free(map);
    return FALSE;
  }

  for(int i = 0; i < n; i++)
  {
    const int br = pts[i].row - rmin;
    const int bc = nrows + pts[i].col - cmin;
    if(map[br] >= 0) c->gidx[c->gstart[map[br]] + fill[map[br]]++] = i;
    if(map[bc] >= 0) c->gidx[c->gstart[map[bc]] + fill[map[bc]]++] = i;
  }

  c->nres_line = total;

  free(fill);
  free(count);
  free(map);
  return TRUE;
}

static void _ctx_cleanup(_solve_ctx_t *c)
{
  free(c->u);
  free(c->v);
  free(c->wu);
  free(c->wv);
  free(c->gstart);
  free(c->gidx);
  free(c->ix);
  free(c->iy);
  free(c->dtx);
  free(c->dty);
  dt_lens_warp_cleanup(&c->work);
  memset(c, 0, sizeof(*c));
}

/* Degree of the i-th polynomial coefficient, so regularisation can lean on
   the high orders -- they are the ones with almost no data at the frame
   edge and the most enthusiasm for filling the gap. */
static void _fill_regularization(_solve_ctx_t *c,
                                 const dt_lens_warp_kind_t kind,
                                 const int order)
{
  for(int i = 0; i < DT_LENS_WARP_MAX_PARAMS; i++)
  {
    c->regw[i] = 1.0f;
    c->regt[i] = 0.0f;
  }

  if(kind == DT_LENS_WARP_ANAM_RADIAL)
  {
    /* The ellipticity is a ratio and its neutral value is one, not zero;
       pulling it towards zero would be pulling it towards a singularity. */
    c->regw[0] = 1.0f;
    c->regt[0] = 1.0f;
    for(int i = 1; i < 5; i++) c->regw[i] = (float)(1 << (2 * (i - 1)));
    return;
  }

  const int terms = dt_lens_warp_poly_terms(order);
  int idx = 0;
  for(int d = 2; d <= order; d++)
    for(int i = d; i >= 0; i--)
    {
      const float w = (float)(1 << (2 * (d - 2)));
      c->regw[idx] = w;
      if(terms + idx < DT_LENS_WARP_MAX_PARAMS) c->regw[terms + idx] = w;
      idx++;
    }
}

static gboolean _fit_affine(const double *ix,
                            const double *iy,
                            const double *tx,
                            const double *ty,
                            const int n,
                            double out[6],
                            double *rms);

static void _unpack(_solve_ctx_t *c, const double *x)
{
  for(int i = 0; i < c->nwarp; i++) c->work.p[i] = (float)x[i];

  if(c->solve_centre)
  {
    c->work.cx = (float)x[c->nwarp];
    c->work.cy = (float)x[c->nwarp + 1];
  }
}

/* The residual vector.
 *
 * Line residuals are perpendicular distances to each group's own best fit
 * line, computed in corrected space. Fitting the line to the same points
 * it is measured against is deliberate: we are not testing whether the
 * points lie on a *particular* line, only whether they lie on *a* line.
 */
static void _eval(_solve_ctx_t *c, const double *x, double *r)
{
  _unpack(c, x);

  /* The squeeze is left out on purpose: it is a pure scale along one axis,
     which cannot bend a straight line, so it contributes exactly nothing
     here and would only add an unidentifiable parameter. The shallow copy
     shares the spline layer, which is only ever read. */
  dt_lens_warp_t tmp = c->work;
  tmp.squeeze = 1.0f;

  for(int i = 0; i < c->n; i++)
  {
    float ou, ov;
    dt_lens_warp_apply(&tmp, c->u[i], c->v[i], &ou, &ov);
    c->wu[i] = ou;
    c->wv[i] = ov;
  }

  int out = 0;
  for(int g = 0; g < c->ngroups; g++)
  {
    const int s = c->gstart[g], e = c->gstart[g + 1];
    const int m = e - s;

    double mu = 0.0, mv = 0.0;
    for(int k = s; k < e; k++)
    {
      mu += c->wu[c->gidx[k]];
      mv += c->wv[c->gidx[k]];
    }
    mu /= m;
    mv /= m;

    double suu = 0.0, suv = 0.0, svv = 0.0;
    for(int k = s; k < e; k++)
    {
      const double du = c->wu[c->gidx[k]] - mu;
      const double dv = c->wv[c->gidx[k]] - mv;
      suu += du * du;
      suv += du * dv;
      svv += dv * dv;
    }

    // principal axis of the group, and the unit normal to it
    const double theta = 0.5 * atan2(2.0 * suv, suu - svv);
    const double nu = -sin(theta);
    const double nv = cos(theta);

    for(int k = s; k < e; k++)
    {
      const double du = c->wu[c->gidx[k]] - mu;
      const double dv = c->wv[c->gidx[k]] - mv;
      r[out++] = du * nu + dv * nv;
    }
  }

  /* Gauge fixing.
   *
   * Straightness alone leaves a whole family of equally good answers,
   * since composing the warp with any projective transform bends nothing.
   * A polynomial with enough freedom will wander into that family and
   * come back with a perfect straightness score and a corrected image
   * that is visibly keystoned. Requiring the corrected lattice to be an
   * *affine* image of the chart -- which is what a square on shot of it
   * ought to produce -- pins the answer down.
   */
  if(c->gauge)
  {
    for(int i = 0; i < c->n; i++)
    {
      c->dtx[i] = c->wu[i];
      c->dty[i] = c->wv[i];
    }

    double A[6] = { 1, 0, 0, 0, 1, 0 };
    _fit_affine(c->ix, c->iy, c->dtx, c->dty, c->n, A, NULL);

    const double gw = c->gauge_w;
    for(int i = 0; i < c->n; i++)
    {
      r[out++] = gw * (c->wu[i] - (A[0] * c->ix[i] + A[1] * c->iy[i] + A[2]));
      r[out++] = gw * (c->wv[i] - (A[3] * c->ix[i] + A[4] * c->iy[i] + A[5]));
    }
  }

  // regularization pulls each coefficient towards its neutral value
  for(int i = 0; i < c->nwarp; i++)
    r[out++] = c->regw[i] * (x[i] - c->regt[i]);

  // the centre is free; nothing pulls on it
  if(c->solve_centre)
  {
    r[out++] = 0.0;
    r[out++] = 0.0;
  }
}

static double _cost(const double *r, const int n)
{
  double s = 0.0;
  for(int i = 0; i < n; i++) s += r[i] * r[i];
  return s;
}

/* Levenberg-Marquardt with a numeric Jacobian.
 *
 * Numeric rather than analytic because the objective routes through a line
 * fit per group, and hand differentiating that is a lot of arithmetic to
 * get subtly wrong for a solve that already runs in well under a second.
 */
static int _levmar(_solve_ctx_t *c, double *x, const int max_iter)
{
  const int np = c->npar;
  const int nr = c->nres;

  double *r = calloc(nr, sizeof(double));
  double *rt = calloc(nr, sizeof(double));
  double *J = calloc((size_t)nr * np, sizeof(double));
  double *JtJ = calloc((size_t)np * np, sizeof(double));
  double *A = calloc((size_t)np * np, sizeof(double));
  double *g = calloc(np, sizeof(double));
  double *step = calloc(np, sizeof(double));
  double *xt = calloc(np, sizeof(double));

  if(!r || !rt || !J || !JtJ || !A || !g || !step || !xt)
  {
    free(r); free(rt); free(J); free(JtJ);
    free(A); free(g); free(step); free(xt);
    return 0;
  }

  _eval(c, x, r);
  double cost = _cost(r, nr);
  double lambda = 1e-3;
  int iter = 0;
  gboolean converged = FALSE;

  for(; iter < max_iter && !converged; iter++)
  {
    for(int j = 0; j < np; j++)
    {
      const double h = 1e-5 * MAX(1.0, fabs(x[j]));
      const double keep = x[j];
      x[j] = keep + h;
      _eval(c, x, rt);
      x[j] = keep;

      for(int i = 0; i < nr; i++)
        J[(size_t)i * np + j] = (rt[i] - r[i]) / h;
    }
    // restore the state the unperturbed residuals belong to
    _unpack(c, x);

    memset(JtJ, 0, sizeof(double) * (size_t)np * np);
    memset(g, 0, sizeof(double) * np);

    for(int i = 0; i < nr; i++)
    {
      const double *Ji = J + (size_t)i * np;
      for(int a = 0; a < np; a++)
      {
        g[a] += Ji[a] * r[i];
        for(int b = a; b < np; b++) JtJ[(size_t)a * np + b] += Ji[a] * Ji[b];
      }
    }
    for(int a = 0; a < np; a++)
      for(int b = 0; b < a; b++)
        JtJ[(size_t)a * np + b] = JtJ[(size_t)b * np + a];

    gboolean improved = FALSE;

    for(int attempt = 0; attempt < 12 && !improved; attempt++)
    {
      memcpy(A, JtJ, sizeof(double) * (size_t)np * np);
      for(int a = 0; a < np; a++)
        A[(size_t)a * np + a] += lambda * (JtJ[(size_t)a * np + a] + 1e-9);

      for(int a = 0; a < np; a++) step[a] = -g[a];

      if(!gauss_solve(A, step, np))
      {
        lambda *= 8.0;
        continue;
      }

      for(int a = 0; a < np; a++) xt[a] = x[a] + step[a];

      _eval(c, xt, rt);
      const double ct = _cost(rt, nr);

      if(ct < cost)
      {
        memcpy(x, xt, sizeof(double) * np);
        memcpy(r, rt, sizeof(double) * nr);
        const double gain = cost - ct;
        cost = ct;
        lambda = MAX(1e-9, lambda * 0.3);
        improved = TRUE;

        if(gain < 1e-14) converged = TRUE;
      }
      else
        lambda *= 8.0;
    }

    if(!improved) break;
  }

  _unpack(c, x);

  free(r); free(rt); free(J); free(JtJ);
  free(A); free(g); free(step); free(xt);
  return iter;
}

/* Least squares affine map from the ideal lattice to the corrected points.
   Solved as two independent three parameter systems, which is all an
   affine is once the two output axes are separated. */
static gboolean _fit_affine(const double *ix,
                            const double *iy,
                            const double *tx,
                            const double *ty,
                            const int n,
                            double out[6],
                            double *rms)
{
  if(n < 3) return FALSE;

  double M[9] = { 0 };
  double bx[3] = { 0 }, by[3] = { 0 };

  for(int i = 0; i < n; i++)
  {
    const double b[3] = { ix[i], iy[i], 1.0 };
    for(int a = 0; a < 3; a++)
    {
      for(int c = 0; c < 3; c++) M[a * 3 + c] += b[a] * b[c];
      bx[a] += b[a] * tx[i];
      by[a] += b[a] * ty[i];
    }
  }

  double Mx[9], My[9];
  memcpy(Mx, M, sizeof(M));
  memcpy(My, M, sizeof(M));

  if(!gauss_solve(Mx, bx, 3) || !gauss_solve(My, by, 3)) return FALSE;

  out[0] = bx[0]; out[1] = bx[1]; out[2] = bx[2];
  out[3] = by[0]; out[4] = by[1]; out[5] = by[2];

  double acc = 0.0;
  for(int i = 0; i < n; i++)
  {
    const double ex = out[0] * ix[i] + out[1] * iy[i] + out[2] - tx[i];
    const double ey = out[3] * ix[i] + out[4] * iy[i] + out[5] - ty[i];
    acc += ex * ex + ey * ey;
  }
  if(rms) *rms = sqrt(acc / n);

  return TRUE;
}

void dt_lens_solve_default_options(dt_lens_solve_options_t *opt)
{
  if(!opt) return;
  memset(opt, 0, sizeof(*opt));
  opt->kind = DT_LENS_WARP_POLY;
  opt->order = 4;
  opt->solve_centre = TRUE;
  opt->chart_frontal = TRUE;
  opt->regularization = 1e-4f;
  opt->tps_smooth = 1e-5f;
  opt->max_iter = 60;
}

gboolean dt_lens_solve(const dt_lens_solve_input_t *in,
                       const dt_lens_solve_options_t *opt,
                       dt_lens_warp_t *warp,
                       dt_lens_solve_result_t *res)
{
  if(res) memset(res, 0, sizeof(*res));
  if(!in || !opt || !warp || !in->points || in->count < 8) return FALSE;
  if(in->width < 2 || in->height < 2) return FALSE;

  const double hd = 0.5 * hypot((double)in->width, (double)in->height);
  const double halfw = 0.5 * in->width;
  const double halfh = 0.5 * in->height;

  _solve_ctx_t c;
  memset(&c, 0, sizeof(c));
  c.n = in->count;
  c.u = calloc(c.n, sizeof(float));
  c.v = calloc(c.n, sizeof(float));
  c.wu = calloc(c.n, sizeof(float));
  c.wv = calloc(c.n, sizeof(float));

  if(!c.u || !c.v || !c.wu || !c.wv)
  {
    _ctx_cleanup(&c);
    return FALSE;
  }

  for(int i = 0; i < c.n; i++)
  {
    c.u[i] = (float)((in->points[i].x - halfw) / hd);
    c.v[i] = (float)((in->points[i].y - halfh) / hd);
  }

  if(!_build_groups(&c, in->points, c.n))
  {
    _ctx_cleanup(&c);
    return FALSE;
  }

  dt_lens_warp_init(&c.work, opt->kind, opt->order);
  c.nwarp = c.work.nparams;
  c.solve_centre = opt->solve_centre;
  c.npar = c.nwarp + (opt->solve_centre ? 2 : 0);

  const double aspect = (in->cell_aspect > 0.01f) ? in->cell_aspect : 1.0;

  c.gauge = opt->chart_frontal;
  c.gauge_w = 1.0; // recalibrated below, once both costs are known
  if(c.gauge)
  {
    c.ix = calloc(c.n, sizeof(double));
    c.iy = calloc(c.n, sizeof(double));
    c.dtx = calloc(c.n, sizeof(double));
    c.dty = calloc(c.n, sizeof(double));

    if(!c.ix || !c.iy || !c.dtx || !c.dty)
    {
      _ctx_cleanup(&c);
      return FALSE;
    }

    for(int i = 0; i < c.n; i++)
    {
      c.ix[i] = in->points[i].col * aspect;
      c.iy[i] = in->points[i].row;
    }
  }

  c.nres = c.nres_line + c.npar + (c.gauge ? 2 * c.n : 0);

  if(c.npar < 1 || c.npar > DT_LENS_SOLVE_MAX_PARAMS)
  {
    _ctx_cleanup(&c);
    return FALSE;
  }

  _fill_regularization(&c, opt->kind, opt->order);
  const float reg = MAX(0.0f, opt->regularization);
  for(int i = 0; i < DT_LENS_WARP_MAX_PARAMS; i++) c.regw[i] *= reg;

  double *x = calloc(c.npar, sizeof(double));
  double *r0 = calloc(c.nres, sizeof(double));
  if(!x || !r0)
  {
    free(x);
    free(r0);
    _ctx_cleanup(&c);
    return FALSE;
  }

  for(int i = 0; i < c.nwarp; i++) x[i] = c.work.p[i];

  // straightness of the uncorrected observations, for comparison
  _eval(&c, x, r0);
  const double line_cost = _cost(r0, c.nres_line);
  const double before = (c.nres_line > 0)
    ? sqrt(line_cost / c.nres_line) : 0.0;

  /* Set the gauge weight from the two costs as they actually stand, rather
     than guessing a constant. The gauge only has to break a tie, so it is
     scaled to a tenth of the straightness cost; at equal weight it is fifty
     times larger and the solver optimises the wrong thing. */
  if(c.gauge)
  {
    double gauge_cost = 0.0;
    for(int i = c.nres_line; i < c.nres_line + 2 * c.n; i++)
      gauge_cost += r0[i] * r0[i];

    if(gauge_cost > 1e-30 && line_cost > 1e-30)
      c.gauge_w = CLAMP(sqrt(0.1 * line_cost / gauge_cost), 1e-4, 1.0);
  }

  const int iters = _levmar(&c, x, MAX(1, opt->max_iter));

  _eval(&c, x, r0);
  const double after = (c.nres_line > 0)
    ? sqrt(_cost(r0, c.nres_line) / c.nres_line) : 0.0;

  /* The corrected positions in c.wu/c.wv are the ones belonging to the
     final parameters, since _eval was just run with them. */

  float squeeze = 1.0f, measured = 1.0f;
  double affine_rms = 0.0;
  double A[6] = { 1, 0, 0, 0, 1, 0 };
  gboolean have_affine = FALSE;

  {
    double *ix = calloc(c.n, sizeof(double));
    double *iy = calloc(c.n, sizeof(double));
    double *tx = calloc(c.n, sizeof(double));
    double *ty = calloc(c.n, sizeof(double));

    if(ix && iy && tx && ty)
    {
      for(int i = 0; i < c.n; i++)
      {
        ix[i] = in->points[i].col * aspect;
        iy[i] = in->points[i].row;
        tx[i] = c.wu[i];
        ty[i] = c.wv[i];
      }

      have_affine = _fit_affine(ix, iy, tx, ty, c.n, A, &affine_rms);

      if(have_affine)
      {
        /* The images of the chart's two axes. Their lengths differ when
           the recorded frame is anisotropically scaled, which for a
           square on chart means an anamorphic squeeze -- and for a tilted
           one means foreshortening, which is why affine_rms is reported
           alongside and has to be small for this number to mean anything. */
        const double sx = hypot(A[0], A[3]);
        const double sy = hypot(A[1], A[4]);
        if(sx > 1e-9 && sy > 1e-9)
          measured = (float)CLAMP(sy / sx, 0.1, 10.0);
      }
    }

    free(ix);
    free(iy);
    free(tx);
    free(ty);
  }

  /* A declared ratio wins over a measured one. It is not a fitted quantity in
     any case -- the straightness objective is blind to it -- and the cell
     shape route needs the chart to have been square on, which the barrel
     engraving does not. Measuring it anyway costs nothing and gives the user
     something to check the declaration against. */
  if(opt->known_squeeze > 0.01f)
    squeeze = opt->known_squeeze;
  else if(opt->chart_frontal)
    squeeze = measured;

  // hand the fitted model to the caller
  dt_lens_warp_init(warp, opt->kind, opt->order);
  warp->order = c.work.order;
  warp->nparams = c.work.nparams;
  memcpy(warp->p, c.work.p, sizeof(warp->p));
  warp->cx = c.work.cx;
  warp->cy = c.work.cy;
  warp->squeeze = squeeze;

  /* This entry describes geometry, which is what makes it eligible for the
     focal length interpolation. Without the flag a freshly solved lens would
     be written out as carrying no distortion at all. */
  warp->have_geometry = TRUE;

  /* The spline layer models what is left over once the parametric fit has
     done its best, so it can only be fitted after that fit exists. */
  if(opt->kind == DT_LENS_WARP_TPS && have_affine)
  {
    float *src = calloc(2 * c.n, sizeof(float));
    float *dx = calloc(c.n, sizeof(float));
    float *dy = calloc(c.n, sizeof(float));

    if(src && dx && dy)
    {
      for(int i = 0; i < c.n; i++)
      {
        const double gx = in->points[i].col * aspect;
        const double gy = in->points[i].row;
        const double wantx = A[0] * gx + A[1] * gy + A[2];
        const double wanty = A[3] * gx + A[4] * gy + A[5];

        // control points live in the same centred space the model uses
        src[2 * i] = c.u[i] - warp->cx;
        src[2 * i + 1] = c.v[i] - warp->cy;
        dx[i] = (float)(wantx - c.wu[i]);
        dy[i] = (float)(wanty - c.wv[i]);
      }

      if(!dt_lens_warp_fit_tps(warp, src, dx, dy, c.n,
                               MAX(0.0f, opt->tps_smooth)))
        dt_print(DT_DEBUG_ALWAYS,
                 "[lens_solve] the spline layer did not solve;"
                 " keeping the polynomial fit alone");
    }

    free(src);
    free(dx);
    free(dy);
  }

  if(res)
  {
    res->ok = TRUE;
    res->rms_px = (float)(after * hd);
    res->rms_before_px = (float)(before * hd);
    res->affine_rms_px = (float)(affine_rms * hd);
    res->squeeze = squeeze;
    res->squeeze_measured = measured;
    res->iterations = iters;
    res->used_points = c.n;
    res->lines_used = c.ngroups;
  }

  dt_print(DT_DEBUG_ALWAYS,
           "[lens_solve] %s order %d: %d points on %d lines,"
           " %d iterations, straightness %.3f -> %.3f px,"
           " affine residual %.3f px, squeeze %.4f (measured %.4f)",
           dt_lens_warp_kind_name(opt->kind), opt->order, c.n, c.ngroups,
           iters, before * hd, after * hd, affine_rms * hd, squeeze, measured);

  free(x);
  free(r0);
  _ctx_cleanup(&c);
  return TRUE;
}

void dt_lens_solve_residuals(const dt_lens_solve_input_t *in,
                             const dt_lens_warp_t *warp,
                             float *out_dx,
                             float *out_dy,
                             float *out_mag)
{
  if(!in || !in->points || in->count < 1) return;

  for(int i = 0; i < in->count; i++)
  {
    if(out_dx) out_dx[i] = 0.0f;
    if(out_dy) out_dy[i] = 0.0f;
    if(out_mag) out_mag[i] = 0.0f;
  }

  const double hd = 0.5 * hypot((double)in->width, (double)in->height);
  const double halfw = 0.5 * in->width;
  const double halfh = 0.5 * in->height;

  _solve_ctx_t c;
  memset(&c, 0, sizeof(c));
  c.n = in->count;
  c.u = calloc(c.n, sizeof(float));
  c.v = calloc(c.n, sizeof(float));
  c.wu = calloc(c.n, sizeof(float));
  c.wv = calloc(c.n, sizeof(float));

  if(!c.u || !c.v || !c.wu || !c.wv || !_build_groups(&c, in->points, c.n))
  {
    _ctx_cleanup(&c);
    return;
  }

  for(int i = 0; i < c.n; i++)
  {
    c.u[i] = (float)((in->points[i].x - halfw) / hd);
    c.v[i] = (float)((in->points[i].y - halfh) / hd);
  }

  dt_lens_warp_t tmp;
  memset(&tmp, 0, sizeof(tmp));
  if(warp) tmp = *warp; // shallow: the spline layer is only read
  tmp.squeeze = 1.0f;

  for(int i = 0; i < c.n; i++)
  {
    float ou = c.u[i], ov = c.v[i];
    if(warp) dt_lens_warp_apply(&tmp, c.u[i], c.v[i], &ou, &ov);
    c.wu[i] = ou;
    c.wv[i] = ov;
  }

  /* A point sits on two lines, its row and its column, and is generally
     off both. Accumulate the two offsets: the sum points away from where
     the correction still has work to do, which is what the overlay is
     for. */
  for(int g = 0; g < c.ngroups; g++)
  {
    const int s = c.gstart[g], e = c.gstart[g + 1];
    const int m = e - s;

    double mu = 0.0, mv = 0.0;
    for(int k = s; k < e; k++)
    {
      mu += c.wu[c.gidx[k]];
      mv += c.wv[c.gidx[k]];
    }
    mu /= m;
    mv /= m;

    double suu = 0.0, suv = 0.0, svv = 0.0;
    for(int k = s; k < e; k++)
    {
      const double du = c.wu[c.gidx[k]] - mu;
      const double dv = c.wv[c.gidx[k]] - mv;
      suu += du * du;
      suv += du * dv;
      svv += dv * dv;
    }

    const double theta = 0.5 * atan2(2.0 * suv, suu - svv);
    const double nu = -sin(theta), nv = cos(theta);

    for(int k = s; k < e; k++)
    {
      const int i = c.gidx[k];
      const double du = c.wu[i] - mu;
      const double dv = c.wv[i] - mv;
      const double d = du * nu + dv * nv;

      // the offset back onto the line, in pixels
      if(out_dx) out_dx[i] += (float)(-d * nu * hd);
      if(out_dy) out_dy[i] += (float)(-d * nv * hd);
    }
  }

  if(out_mag && out_dx && out_dy)
    for(int i = 0; i < c.n; i++)
      out_mag[i] = hypotf(out_dx[i], out_dy[i]);

  _ctx_cleanup(&c);
}

// clang-format off
// modelines: These editor modelines have been set for all relevant files by tools/update_modelines.py
// vim: shiftwidth=2 expandtab tabstop=2 cindent
// kate: tab-indents: off; indent-width 2; replace-tabs on; indent-mode cstyle; remove-trailing-spaces modified;
// clang-format on
