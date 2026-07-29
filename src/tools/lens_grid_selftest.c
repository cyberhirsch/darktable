/* Run the real grid detector on a real photograph, outside the gui.
 *
 * Loads an image with GdkPixbuf, converts to luminance, calls
 * dt_lens_grid_detect and prints what came out -- including a coverage map
 * of which lattice nodes were found, which is the thing the on-screen
 * overlay makes obvious and a bare intersection count does not.
 */
#include <stdio.h>
#include <math.h>
#include <gdk-pixbuf/gdk-pixbuf.h>

#include "common/lens_grid.c"

void *dt_alloc_aligned(size_t size)
{
  return g_malloc(size);
}

void dt_free_align(void *p)
{
  g_free(p);
}

void dt_print_ext(const char *msg, ...)
{
  va_list ap;
  va_start(ap, msg);
  vprintf(msg, ap);
  va_end(ap);
  printf("\n");
}

int main(int argc, char **argv)
{
  if(argc < 4)
  {
    fprintf(stderr, "usage: %s image.jpg cells_x cells_y [long_edge]\n",
            argv[0]);
    return 2;
  }

  const int cells_x = atoi(argv[2]);
  const int cells_y = atoi(argv[3]);
  const int target = (argc > 4) ? atoi(argv[4]) : 3600;

  GError *err = NULL;
  GdkPixbuf *pb = gdk_pixbuf_new_from_file_at_size(argv[1], target, target,
                                                   &err);
  if(!pb)
  {
    fprintf(stderr, "cannot load %s: %s\n", argv[1],
            err ? err->message : "?");
    return 1;
  }

  const int w = gdk_pixbuf_get_width(pb);
  const int h = gdk_pixbuf_get_height(pb);
  const int stride = gdk_pixbuf_get_rowstride(pb);
  const int nch = gdk_pixbuf_get_n_channels(pb);
  const guchar *const px = gdk_pixbuf_get_pixels(pb);

  float *lum = malloc(sizeof(float) * (size_t)w * h);
  for(int y = 0; y < h; y++)
    for(int x = 0; x < w; x++)
    {
      const guchar *p = px + (size_t)y * stride + (size_t)x * nch;
      lum[(size_t)y * w + x] =
        (0.299f * p[0] + 0.587f * p[1] + 0.114f * p[2]) / 255.0f;
    }

  printf("image %dx%d, chart %dx%d cells (%dx%d lattice)\n\n",
         w, h, cells_x, cells_y, cells_x + 1, cells_y + 1);

  /* Diagnostic: what does one scanline actually look like? Counting lines
     tells you something is wrong; the profile tells you what. */
  if(getenv("GRIDTEST_PROFILE"))
  {
    const float spacing_x = (float)w / (cells_x + 1);
    const float spacing_y = (float)h / (cells_y + 1);
    const int radius = (int)fmaxf(4.0f, fminf(spacing_x, spacing_y));
    float *norm = _normalize_local(lum, w, h, radius);
    const int margin = MAX(4, (int)(0.012f * MAX(w, h)));
    const float mad = _estimate_line_depth(norm, w, h, margin);
    const float threshold = MAX(1e-6f, mad * 0.25f);

    printf("expected spacing: %.1f px across, %.1f px down\n",
           spacing_x, spacing_y);
    printf("radius %d, mad %.5f, threshold %.5f\n\n", radius, mad, threshold);

    // scan down the middle column, so we meet the horizontal lines
    const int col = w / 2;
    float cross[4096], strength[4096];
    const int nc = _find_crossings(norm + col, h, w, threshold, margin,
                                   cross, strength, 4096);

    printf("middle column: %d crossings\n", nc);
    for(int i = 0; i < nc; i++)
      printf("  %7.2f  depth %.4f%s\n", cross[i], strength[i],
             i ? "" : "");
    printf("\ngaps between them:\n ");
    for(int i = 1; i < nc; i++) printf(" %.1f", cross[i] - cross[i - 1]);
    printf("\n");

    // and the raw signal across one cell, to see the line shape
    if(nc >= 2)
    {
      const int c0 = (int)cross[nc / 2] - 25;
      printf("\nsignal around crossing %d (row %d..%d):\n",
             nc / 2, c0, c0 + 50);
      for(int y = c0; y < c0 + 50 && y < h; y++)
      {
        if(y < 0) continue;
        const float v = norm[(size_t)y * w + col];
        printf("  %5d % .5f %s%s\n", y, v,
               v < -threshold ? "*" : " ",
               fabsf(v) > 0.0f ? "" : "");
      }
    }

    /* Where do the curves come from? Raw traces, then after stitching,
       with each one's position and extent, so a triplet at the same place
       is distinguishable from three real lines. */
    int nx = 0;
    GList *raw = _trace_curves(norm, w, h, TRUE, threshold, spacing_y,
                               margin, &nx);
    printf("raw horizontal traces: %u\n", g_list_length(raw));

    raw = g_list_sort(raw, _curve_compare_start);
    int shown = 0;
    for(GList *it = raw; it && shown < 45; it = g_list_next(it), shown++)
    {
      const dt_lens_grid_curve_t *cv = it->data;
      printf("  start %6.0f end %6.0f  n=%4d  coord %7.2f -> %7.2f\n",
             cv->samples[0].pos, cv->samples[cv->count - 1].pos, cv->count,
             cv->samples[0].coord, cv->samples[cv->count - 1].coord);
    }

    raw = _merge_fragments(raw, MAX(3.0f, spacing_y * 0.45f), 0.25f * w);
    printf("after stitching: %u\n", g_list_length(raw));

    /* Where do the surviving curves actually sit? Clustered means one real
       line is being represented several times; spread means junk. */
    GList *pruned = _prune_and_order(raw, w, cells_y + 1, spacing_y);
    printf("after pruning: %u\n", g_list_length(pruned));
    printf("  ordinal  coord@centre   span   n\n");
    for(GList *it = pruned; it; it = g_list_next(it))
    {
      const dt_lens_grid_curve_t *cv = it->data;
      printf("  %7d  %11.2f  %6.0f  %4d\n", cv->ordinal,
             _curve_coord_at(cv, 0.5f * w),
             cv->samples[cv->count - 1].pos - cv->samples[0].pos, cv->count);
    }

    g_list_free_full(pruned, _curve_free);

    dt_free_align(norm);
    printf("\n----\n\n");
  }

  dt_lens_grid_t grid;
  const gboolean ok =
    dt_lens_grid_detect(lum, w, h, cells_x, cells_y, TRUE, &grid);

  printf("\ndetect returned %s, %d intersections\n",
         ok ? "TRUE" : "FALSE", grid.point_count);

  if(grid.point_count > 0)
  {
    /* Coverage map. Rows are lattice rows, columns lattice columns; a '#'
       is a node that was found. Holes show where detection failed, which
       is what actually needs fixing. */
    int rmax = 0, cmax = 0;
    for(int i = 0; i < grid.point_count; i++)
    {
      if(grid.points[i].row > rmax) rmax = grid.points[i].row;
      if(grid.points[i].col > cmax) cmax = grid.points[i].col;
    }

    const int rows = rmax + 1, cols = cmax + 1;
    char *map = calloc((size_t)rows * cols, 1);
    for(int i = 0; i < rows * cols; i++) map[i] = '.';
    for(int i = 0; i < grid.point_count; i++)
      map[(size_t)grid.points[i].row * cols + grid.points[i].col] = '#';

    printf("\ncoverage, lattice indices 0..%d x 0..%d"
           " (expected 0..%d x 0..%d):\n",
           cmax, rmax, cells_x, cells_y);
    for(int r = 0; r < rows; r++)
    {
      printf("  %2d ", r);
      for(int c = 0; c < cols; c++) putchar(map[(size_t)r * cols + c]);
      putchar('\n');
    }
    free(map);

    const int expected = (cells_x + 1) * (cells_y + 1);
    printf("\nfound %d of %d expected nodes (%.0f%%)\n",
           grid.point_count, expected,
           100.0 * grid.point_count / expected);
  }

  dt_lens_grid_cleanup(&grid);
  free(lum);
  g_object_unref(pb);
  return 0;
}
