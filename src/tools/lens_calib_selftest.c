/* Standalone numerical check of the warp model and the plumb-line solver.
 *
 * Not part of the build. It compiles the two implementation files directly
 * and stubs the handful of darktable symbols they reach for, so it can run
 * without starting the application -- which is the point: these are pure
 * numerical routines whose failures are silent in the gui, showing up as a
 * correction that merely looks a bit wrong.
 *
 * Build and run from the repository root:
 *
 *   gcc -O1 -std=c99 -o selftest src/tools/lens_calib_selftest.c \
 *     -Isrc -Ibuild/bin -Isrc/external \
 *     $(pkg-config --cflags glib-2.0 json-glib-1.0 gtk+-3.0 librsvg-2.0 lcms2) \
 *     $(pkg-config --libs glib-2.0 json-glib-1.0) -lm
 *   ./selftest
 *
 * The method is to synthesize a lattice, push it through a known warp to
 * make "observed" points, and check that the solver recovers what was put
 * in. The isotropy checks at the end are the ones that matter most: they
 * ask whether the corrected chart comes out square, which is the property
 * a user would actually notice, rather than whether any particular
 * coefficient matched.
 */
#include <stdio.h>
#include <math.h>

#include "common/lens_warp.c"
#include "common/lens_solve.c"

/* stubs */
void dt_print_ext(const char *msg, ...) { }
void dt_loc_get_user_config_dir(char *d, size_t n) { g_strlcpy(d, ".", n); }

#define W 6000
#define H 4000
#define CX 29
#define CY 13
#define STEP 0.05

static int fails = 0;

/* The property that actually matters for an image: after the full warp,
   squeeze included, is the lattice a *similarity* of the square chart --
   equal scale on both axes and no shear? The squeeze number itself is
   only a means to that end, and different backends are free to split the
   work between the model and the squeeze stage differently. */
static void _check_isotropy(const char *label,
                            const dt_lens_warp_t *w,
                            const dt_lens_solve_point_t *pts,
                            const int n,
                            const int width,
                            const int height,
                            const double step);
static void check(const char *what, double got, double want, double tol)
{
  const int ok = fabs(got - want) <= tol;
  if(!ok) fails++;
  printf("  %-42s %12.6f (want %10.6f +-%.4g) %s\n",
         what, got, want, tol, ok ? "ok" : "FAIL");
}

int main(void)
{
  /* --- 1. apply/invert round trip on a strong barrel + squeeze --- */
  dt_lens_warp_t w;
  dt_lens_warp_init(&w, DT_LENS_WARP_ANAM_RADIAL, 4);
  w.p[0] = 1.8f;    /* ellipticity */
  w.p[1] = -0.12f;  /* k1 */
  w.p[2] = 0.03f;   /* k2 */
  w.cx = 0.01f; w.cy = -0.008f;
  w.squeeze = 2.0f;

  double worst = 0.0;
  for(int i = -8; i <= 8; i++)
    for(int j = -6; j <= 6; j++)
    {
      const float u = i * 0.1f, v = j * 0.1f;
      float au, av, bu, bv;
      dt_lens_warp_apply(&w, u, v, &au, &av);
      if(!dt_lens_warp_invert(&w, au, av, &bu, &bv)) { printf("  no convergence at %g %g\n", u, v); fails++; continue; }
      const double e = hypot(bu - u, bv - v);
      if(e > worst) worst = e;
    }
  printf("anamorphic round trip\n");
  check("worst inverse error (normalized)", worst, 0.0, 1e-5);

  /* --- 2. solver recovers a known barrel distortion --- */
  /* Synthesize a lattice, push it through a known *inverse* warp to make
     "observed" points, then see if the solver undoes it. */
  const double hd = 0.5 * hypot((double)W, (double)H);
  dt_lens_solve_point_t pts[CX * CY];
  int n = 0;

  dt_lens_warp_t truth;
  dt_lens_warp_init(&truth, DT_LENS_WARP_ANAM_RADIAL, 4);
  truth.p[0] = 1.0f;
  truth.p[1] = 0.09f;   /* correcting barrel */
  truth.p[2] = -0.02f;
  truth.cx = 0.012f; truth.cy = -0.006f;
  truth.squeeze = 1.0f;

  for(int r = 0; r < CY; r++)
    for(int c = 0; c < CX; c++)
    {
      /* ideal (corrected) position of this node, filling ~80% of frame */
      /* square cells: one step for both axes, or the solver is being
         told the chart is square when the test data says otherwise */
      const double iu = (c - (CX - 1) * 0.5) * STEP;
      const double iv = (r - (CY - 1) * 0.5) * STEP;
      float ou, ov;
      /* observed = inverse of the correction */
      dt_lens_warp_invert(&truth, (float)iu, (float)iv, &ou, &ov);
      pts[n].x = (float)(ou * hd + 0.5 * W);
      pts[n].y = (float)(ov * hd + 0.5 * H);
      pts[n].col = c;
      pts[n].row = r;
      n++;
    }

  dt_lens_solve_input_t in = { pts, n, W, H, 1.0f };
  dt_lens_solve_options_t opt;
  dt_lens_solve_default_options(&opt);
  opt.kind = DT_LENS_WARP_ANAM_RADIAL;
  opt.regularization = 1e-6f;

  dt_lens_warp_t fit;
  dt_lens_solve_result_t res;
  const gboolean ok = dt_lens_solve(&in, &opt, &fit, &res);

  printf("solver, anamorphic radial on %d synthetic points\n", n);
  if(!ok) { printf("  solve FAILED\n"); fails++; }
  else
  {
    check("straightness before (px)", res.rms_before_px, 8.0, 1e9);
    check("straightness after  (px)", res.rms_px, 0.0, 0.05);
    check("squeeze", res.squeeze, 1.0, 0.01);
  }

  /* --- 3. same but with a real 2x squeeze in the observations --- */
  n = 0;
  truth.squeeze = 2.0f;
  for(int r = 0; r < CY; r++)
    for(int c = 0; c < CX; c++)
    {
      /* square cells: one step for both axes, or the solver is being
         told the chart is square when the test data says otherwise */
      const double iu = (c - (CX - 1) * 0.5) * STEP;
      const double iv = (r - (CY - 1) * 0.5) * STEP;
      float ou, ov;
      dt_lens_warp_invert(&truth, (float)iu, (float)iv, &ou, &ov);
      pts[n].x = (float)(ou * hd + 0.5 * W);
      pts[n].y = (float)(ov * hd + 0.5 * H);
      pts[n].col = c; pts[n].row = r; n++;
    }

  dt_lens_warp_t fit2;
  dt_lens_solve_result_t res2;
  if(dt_lens_solve(&in, &opt, &fit2, &res2))
  {
    printf("solver, with a 2x squeeze present\n");
    check("straightness after  (px)", res2.rms_px, 0.0, 0.05);
    check("recovered squeeze", res2.squeeze, 2.0, 0.02);
  }
  else { printf("  squeeze solve FAILED\n"); fails++; }

  /* --- 4. polynomial backend on the same data --- */
  opt.kind = DT_LENS_WARP_POLY;
  opt.order = 4;
  dt_lens_warp_t fit3;
  dt_lens_solve_result_t res3;
  if(dt_lens_solve(&in, &opt, &fit3, &res3))
  {
    printf("solver, bivariate polynomial order 4\n");
    check("straightness after  (px)", res3.rms_px, 0.0, 0.30);
    check("recovered squeeze", res3.squeeze, 2.0, 0.05);
  }
  else { printf("  poly solve FAILED\n"); fails++; }

  /* The truth carries a k2 r^4 term, so the displacement it produces is
     degree 5 in the coordinates. An order 4 polynomial structurally cannot
     represent that, and the misfit leaks into the affine stage as a squeeze
     error. This run separates a modelling limit from a bug. */
  opt.order = 5;
  dt_lens_warp_t fit4;
  dt_lens_solve_result_t res4;
  if(dt_lens_solve(&in, &opt, &fit4, &res4))
  {
    printf("solver, bivariate polynomial order 5\n");
    check("straightness after  (px)", res4.rms_px, 0.0, 0.05);
    check("recovered squeeze", res4.squeeze, 2.0, 0.02);
  }
  else { printf("  poly5 solve FAILED\n"); fails++; }

  /* --- 5. spline backend --- */
  opt.kind = DT_LENS_WARP_TPS;
  opt.order = 4;
  dt_lens_warp_t fit5;
  dt_lens_solve_result_t res5;
  if(dt_lens_solve(&in, &opt, &fit5, &res5))
  {
    printf("solver, spline layer on a polynomial base\n");
    check("control points", fit5.tps_count, n, 0);
    check("straightness after  (px)", res5.rms_px, 0.0, 0.30);
  }
  else { printf("  tps solve FAILED\n"); fails++; }

  /* --- 6. profile round trip through json --- */
  {
    dt_lens_profile_t prof;
    dt_lens_profile_init(&prof);
    g_strlcpy(prof.name, "test", sizeof(prof.name));
    prof.width = W; prof.height = H;

    dt_lens_warp_t a, b;
    dt_lens_warp_copy(&a, &fit2); a.focal = 24.0f;
    dt_lens_warp_copy(&b, &fit2); b.focal = 96.0f;
    b.p[1] *= 0.5f;                 /* less distortion at the long end */
    b.squeeze = 1.5f;
    dt_lens_profile_add(&prof, &a);
    dt_lens_profile_add(&prof, &b);

    GError *e = NULL;
    const gboolean wrote = dt_lens_profile_save(&prof, "warptest_profile.json", &e);
    if(!wrote) { printf("  profile save FAILED\n"); fails++; }

    dt_lens_profile_t back;
    if(dt_lens_profile_load(&back, "warptest_profile.json"))
    {
      printf("profile round trip\n");
      check("warps read back", back.warps->len, 2, 0);

      const dt_lens_warp_t *r0 = &g_array_index(back.warps, dt_lens_warp_t, 0);
      check("first focal", r0->focal, 24.0, 1e-4);
      check("first k1 preserved", r0->p[1], a.p[1], 1e-6);

      /* Interpolation is linear in log focal, so the midpoint by that
         measure is sqrt(24*96) = 48, not 60. */
      dt_lens_warp_t mid;
      if(dt_lens_profile_eval(&back, 48.0f, &mid))
      {
        check("squeeze halfway (log focal)", mid.squeeze,
              0.5 * (a.squeeze + b.squeeze), 1e-4);
        check("k1 halfway (log focal)", mid.p[1],
              0.5 * (a.p[1] + b.p[1]), 1e-6);
        dt_lens_warp_cleanup(&mid);
      }
      else { printf("  eval FAILED\n"); fails++; }

      /* Outside the measured range it holds rather than extrapolating. */
      dt_lens_warp_t wide;
      if(dt_lens_profile_eval(&back, 10.0f, &wide))
      {
        check("below range holds the widest", wide.p[1], a.p[1], 1e-6);
        dt_lens_warp_cleanup(&wide);
      }
    }
    else { printf("  profile load FAILED\n"); fails++; }

    dt_lens_warp_cleanup(&a);
    dt_lens_warp_cleanup(&b);
    dt_lens_profile_cleanup(&prof);
    dt_lens_profile_cleanup(&back);
    g_unlink("warptest_profile.json");
  }

  _check_isotropy("anam_radial", &fit2, pts, n, W, H, STEP);
  _check_isotropy("poly order 4", &fit3, pts, n, W, H, STEP);
  _check_isotropy("poly order 5", &fit4, pts, n, W, H, STEP);

  printf("\n%s (%d failures)\n", fails ? "FAILURES" : "all checks passed", fails);
  return fails != 0;
}

static void _check_isotropy(const char *label,
                            const dt_lens_warp_t *w,
                            const dt_lens_solve_point_t *pts,
                            const int n,
                            const int width,
                            const int height,
                            const double step)
{
  const double hd = 0.5 * hypot((double)width, (double)height);

  double *ix = calloc(n, sizeof(double)), *iy = calloc(n, sizeof(double));
  double *tx = calloc(n, sizeof(double)), *ty = calloc(n, sizeof(double));

  for(int i = 0; i < n; i++)
  {
    float ou, ov;
    dt_lens_warp_apply(w,
                       (float)((pts[i].x - 0.5 * width) / hd),
                       (float)((pts[i].y - 0.5 * height) / hd), &ou, &ov);
    tx[i] = ou;
    ty[i] = ov;
    ix[i] = pts[i].col * step;
    iy[i] = pts[i].row * step;
  }

  double A[6], rms = 0.0;
  _fit_affine(ix, iy, tx, ty, n, A, &rms);

  const double sx = hypot(A[0], A[3]);
  const double sy = hypot(A[1], A[4]);
  // angle between the images of the two chart axes; 90 degrees means no shear
  const double shear =
    fabs(90.0 - fabs(atan2(A[3], A[0]) - atan2(A[4], A[1])) * 180.0 / G_PI);

  printf("corrected lattice is square? (%s)\n", label);
  check("  axis scale ratio", sy / sx, 1.0, 0.005);
  check("  shear (degrees)", shear, 0.0, 0.2);
  check("  affine residual (px)", rms * hd, 0.0, 0.5);

  free(ix); free(iy); free(tx); free(ty);
}
