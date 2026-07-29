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

/* Lens calibration view.
 *
 * Builds lens correction profiles from a photograph of a grid chart. The
 * chart's cell counts are user supplied, so any regular grid works; the
 * default matches the 28x12 chart shipped by Chaos.
 *
 * This file owns the view and its interaction. The grid detection, the
 * warp models and the profile i/o live in common/ so they can be reused
 * by the lens iop and by tests.
 */

#include "common/act_on.h"
#include "common/collection.h"
#include "common/darktable.h"
#include "common/debug.h"
#include "common/file_location.h"
#include "common/image_cache.h"
#include "common/lens_grid.h"
#include "common/lens_share.h"
#include "common/lens_lensfun.h"
#include "libs/lens_calib_common.h"
#include "common/lens_solve.h"
#include "common/lens_stmap.h"
#include "common/lens_vignette.h"
#include "common/lens_warp.h"
#include "common/mipmap_cache.h"
#include "control/conf.h"
#include "control/control.h"
#include "develop/format.h"
#include "dtgtk/thumbtable.h"
#include "gui/accelerators.h"
#include "gui/gtk.h"
#include "iop/gaussian_elimination.h"
#include "views/view.h"
#include "views/view_api.h"

#include <gdk/gdkkeysyms.h>
#include <glib/gstdio.h>
#include <float.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

DT_MODULE(1)

// the chart shipped by Chaos, and a sane default for anything else
#define DT_LENS_CALIB_DEFAULT_CELLS_X 28
#define DT_LENS_CALIB_DEFAULT_CELLS_Y 12

// long edge used for grid detection, independent of the on-screen size:
// the interior chart lines need enough pixels to survive downscaling
#define DT_LENS_CALIB_DETECT_SIZE 3600

// upper bound on the zoomed display surface, to keep allocations sane
#define DT_LENS_CALIB_MAX_SURFACE 8000

// zoom limits and step per scroll notch
#define DT_LENS_CALIB_ZOOM_MIN 1.0f
#define DT_LENS_CALIB_ZOOM_MAX 16.0f
#define DT_LENS_CALIB_ZOOM_STEP 1.25f

typedef struct dt_lens_calib_t
{
  // the grid photograph being measured
  dt_imgid_t imgid;
  cairo_surface_t *surface;
  int surf_width, surf_height;
  // image dimensions the surface was requested at, so we know when the
  // centre view changed size and the surface has to be rebuilt
  int req_width, req_height;

  /* Where the image ended up on screen. Every overlay and every hit test
     goes through these four numbers, so zoom and pan are implemented
     purely by changing them -- nothing downstream needs to know. */
  double img_x, img_y, img_w, img_h;

  /* zoom == 1 fits the image to the viewport. pan is the normalized image
     point held at the centre of the viewport. */
  float zoom;
  float pan_x, pan_y;
  gboolean panning;
  double pan_ref_x, pan_ref_y;   // screen position where the pan began
  float pan_ref_px, pan_ref_py;  // pan value when it began

  // viewport size as of the last expose; scroll handling needs it and
  // only expose is told what it is
  int vp_w, vp_h;

  // user supplied chart geometry. cells, not intersections: a 28x12 chart
  // has a 29x13 lattice.
  int cells_x, cells_y;
  float cell_aspect;

  // detection results, in full image coordinates
  dt_lens_grid_t grid;
  gboolean have_grid;
  int detect_img_w, detect_img_h;

  /* Manually placed intersections, in normalized image coordinates so
     they stay valid when the display size or the detection resolution
     changes. Kept separate from the detected set: hand-placed points are
     the user's ground truth and must never be discarded by a re-detect. */
  GArray *manual;
  int hover_manual; // index under the pointer, or -1

  gboolean show_points;
  gboolean show_curves;    // the detector's traced polylines
  gboolean show_mesh;      // the nodes joined by lattice adjacency

  /* Hand editing. The detector's output and the hand placed set are two
     separate things right up until this is switched on, at which point
     the detected points are adopted into the editable set -- there is no
     use for a half-corrected detection that the solver reads from two
     places at once. */
  gboolean manual_edit;
  int drag_idx;            // point being dragged, or -1

  /* Show the image as the fit says it should look rather than as shot.
     The corrected surface is cached because rebuilding it means an
     inverse warp per pixel. */
  gboolean flat_view;

  /* Undo the brightness falloff in the preview. Independent of flat_view:
     the two are separate corrections and each is easier to judge on its
     own. */
  gboolean falloff_view;

  cairo_surface_t *flat_surface;
  gboolean flat_valid;

  /* Pinned chart corners. Four clicks are enough to lay out the whole
     lattice, which beats placing several hundred nodes by hand. */
  float corner_x[4], corner_y[4];
  int corner_count;
  gboolean corner_mode;
  int corner_drag;

  /* Where the chart sits in the frame, as a homography from lattice index
     to image position. Without it the reference grid is a plain rectangle
     over the whole frame, which bears no relation to where the chart
     actually is and so cannot be compared to anything. */
  double pose[9];
  gboolean have_pose;
  float pose_rms; // how well the pose alone explains the points, in pixels

  /* The fit, and the per point residuals that go with it. The residuals
     are kept alongside the points they belong to rather than recomputed
     at draw time, since the solve is where they are cheap to get. */
  dt_lens_warp_t warp;
  gboolean have_warp;

  /* The profile being assembled, and which of its entries `warp` is a
     working copy of.
     
     A lens is not one measurement. A zoom needs one per focal length, and
     vignetting needs one per aperture on top of that -- so the thing being
     built is a collection, while the fit in progress is a single member of
     it. Keeping the collection here rather than only on disk is what lets
     the panel show the other entries and switch between them, instead of
     each save silently merging into a file nobody can see. */
  dt_lens_profile_t work;
  int sel_entry;
  dt_lens_solve_result_t result;

  dt_lens_solve_point_t *solve_pts;
  float *res_dx, *res_dy, *res_mag;
  int solve_count;

  // full image size the solve was run against
  int solve_img_w, solve_img_h;

  gboolean show_residuals;

  // how well the radial vignetting model described the frame it was fitted to
  float vig_rms;
  float vig_stops; // the correction at the frame corner, in stops

  /* How far out the samples actually reached, as a fraction of the corner
     radius, and how many cells were read. A chart that stops short of the
     corners makes the corner figure an extrapolation, and that has to be
     visible rather than inferred from the residual. */
  float vig_max_r;
  int vig_samples;

  /* Undo. Snapshots of the whole editable state rather than per-edit deltas:
     the state is a few hundred nodes and a pose, so a copy is some tens of
     kilobytes, while the operations that need undoing most -- relax, mesh
     layout, align -- touch every node at once and have no compact delta.
     Kept local to the view rather than going through dt_undo, which
     coalesces items by timestamp and would fuse a burst of drags into one
     step. Here one edit is always one step. */
  GPtrArray *undo_stack, *redo_stack;
  gboolean restoring; // suppress recording while an undo is being applied

  /* Composite operations record once and then hold, so that laying out the
     mesh -- which aligns the grid on the way -- is one step to undo and not
     two. Without it the user would have to press ctrl-z twice to escape a
     single button. */
  int undo_hold;
} dt_lens_calib_t;

typedef struct dt_lens_calib_manual_point_t
{
  float x, y;   // normalized 0..1
  int col, row; // lattice index, -1 when not yet known

  /* TRUE for a point that was detected or placed by hand, FALSE for one
     interpolated from its neighbours. The distinction has to survive,
     because an interpolated node is a restatement of what we already
     assumed and carries no evidence about the lens. */
  gboolean measured;
} dt_lens_calib_manual_point_t;


/* ------------------------------------------------------------- undo */

/* A snapshot of everything an edit can change.
 *
 * Points, corners and pose travel together because they are not independent:
 * aligning the grid re-indexes the points, prunes duplicates, moves the
 * corners and replaces the pose, all in one action. Restoring only the points
 * would leave a pose fitted to nodes that no longer exist, which is a worse
 * state than either the before or the after.
 *
 * The solve is deliberately *not* snapshotted. A fit belongs to the points it
 * was made from, so after an undo the honest thing is to have no fit rather
 * than a stale one that still reports its old residuals.
 */
#define DT_LENS_CALIB_UNDO_MAX 64

typedef struct dt_lens_calib_snapshot_t
{
  dt_lens_calib_manual_point_t *points;
  guint count;

  float corner_x[4], corner_y[4];
  int corner_count;

  double pose[9];
  gboolean have_pose;
  float pose_rms;

  int cells_x, cells_y;

  // what the action was, so the log can say what it is undoing
  char *label;
} dt_lens_calib_snapshot_t;

// defined further down, but cleanup() and leave() need them
static void _save_points(dt_lens_calib_t *d);
static void _load_points(dt_lens_calib_t *d);
static void _drop_solution(dt_lens_calib_t *d);
static void _drop_flat(dt_lens_calib_t *d);
static void _auto_align(dt_view_t *self);
static void _set_corner_mode(dt_view_t *self, const gboolean on);

static void _snapshot_free(dt_lens_calib_snapshot_t *s)
{
  if(!s) return;
  free(s->points);
  g_free(s->label);
  free(s);
}

static dt_lens_calib_snapshot_t *_snapshot_take(const dt_lens_calib_t *d,
                                               const char *label)
{
  dt_lens_calib_snapshot_t *s = calloc(1, sizeof(*s));
  if(!s) return NULL;

  s->count = d->manual ? d->manual->len : 0;
  if(s->count)
  {
    s->points = malloc(sizeof(dt_lens_calib_manual_point_t) * s->count);
    if(!s->points)
    {
      free(s);
      return NULL;
    }
    memcpy(s->points, d->manual->data,
           sizeof(dt_lens_calib_manual_point_t) * s->count);
  }

  for(int i = 0; i < 4; i++)
  {
    s->corner_x[i] = d->corner_x[i];
    s->corner_y[i] = d->corner_y[i];
  }
  s->corner_count = d->corner_count;

  memcpy(s->pose, d->pose, sizeof(s->pose));
  s->have_pose = d->have_pose;
  s->pose_rms = d->pose_rms;
  s->cells_x = d->cells_x;
  s->cells_y = d->cells_y;
  s->label = g_strdup(label ? label : "");

  return s;
}

static void _snapshot_apply(dt_lens_calib_t *d,
                            const dt_lens_calib_snapshot_t *s)
{
  if(!d->manual) return;

  g_array_set_size(d->manual, 0);
  if(s->count) g_array_append_vals(d->manual, s->points, s->count);

  for(int i = 0; i < 4; i++)
  {
    d->corner_x[i] = s->corner_x[i];
    d->corner_y[i] = s->corner_y[i];
  }
  d->corner_count = s->corner_count;

  memcpy(d->pose, s->pose, sizeof(d->pose));
  d->have_pose = s->have_pose;
  d->pose_rms = s->pose_rms;

  d->cells_x = s->cells_x;
  d->cells_y = s->cells_y;
  dt_conf_set_int("plugins/lens_calib/cells_x", d->cells_x);
  dt_conf_set_int("plugins/lens_calib/cells_y", d->cells_y);

  d->hover_manual = -1;
  d->drag_idx = -1;
  d->corner_drag = -1;

  /* A fit belongs to the points it was made from, so restoring different
     points leaves us with no fit rather than a stale one still reporting
     the residuals of a set that no longer exists. */
  _drop_solution(d);
  _drop_flat(d);
  _save_points(d);
}

static void _undo_stack_clear(GPtrArray *a)
{
  if(!a) return;
  for(guint i = 0; i < a->len; i++) _snapshot_free(g_ptr_array_index(a, i));
  g_ptr_array_set_size(a, 0);
}

static void _undo_clear(dt_lens_calib_t *d)
{
  _undo_stack_clear(d->undo_stack);
  _undo_stack_clear(d->redo_stack);
}

/* Record the state as it is *before* an edit, labelled with the edit that is
   about to happen. Call it at the top of anything that changes points,
   corners or pose. */
static void _undo_push(dt_lens_calib_t *d, const char *label)
{
  if(d->restoring || d->undo_hold > 0 || !d->undo_stack) return;

  dt_lens_calib_snapshot_t *s = _snapshot_take(d, label);
  if(!s) return;

  g_ptr_array_add(d->undo_stack, s);

  // oldest first, so trimming the front is what drops the least useful step
  while(d->undo_stack->len > DT_LENS_CALIB_UNDO_MAX)
  {
    _snapshot_free(g_ptr_array_index(d->undo_stack, 0));
    g_ptr_array_remove_index(d->undo_stack, 0);
  }

  /* A new edit invalidates the redo branch: what was redone forward from
     here no longer follows from this state. */
  _undo_stack_clear(d->redo_stack);
}

/* Discard the last recorded step, for an operation that pushed and then bailed
   out without changing anything. An undo step that does nothing is worse than
   no step at all: it makes ctrl-z look broken. */
static void _undo_drop_last(dt_lens_calib_t *d)
{
  // inside a composite operation the last step belongs to the caller
  if(d->undo_hold > 0 || !d->undo_stack || !d->undo_stack->len) return;
  _snapshot_free(g_ptr_array_index(d->undo_stack, d->undo_stack->len - 1));
  g_ptr_array_remove_index(d->undo_stack, d->undo_stack->len - 1);
}

static void _undo_apply(dt_view_t *self, const gboolean redo)
{
  dt_lens_calib_t *d = self->data;
  GPtrArray *from = redo ? d->redo_stack : d->undo_stack;
  GPtrArray *to = redo ? d->undo_stack : d->redo_stack;

  if(!from || !from->len)
  {
    dt_control_log(redo ? _("nothing to redo") : _("nothing to undo"));
    return;
  }

  dt_lens_calib_snapshot_t *s = g_ptr_array_index(from, from->len - 1);
  g_ptr_array_remove_index(from, from->len - 1);

  // the state we are leaving becomes the way back, under the same name
  dt_lens_calib_snapshot_t *here = _snapshot_take(d, s->label);
  if(here) g_ptr_array_add(to, here);

  d->restoring = TRUE;
  _snapshot_apply(d, s);
  d->restoring = FALSE;

  if(s->label && *s->label)
    dt_control_log(redo ? _("redone: %s") : _("undone: %s"), s->label);

  _snapshot_free(s);
  dt_control_queue_redraw_center();
}

static void _undo_callback(dt_action_t *action)
{
  dt_view_t *self = darktable.view_manager->proxy.lens_calib.view;
  if(self) _undo_apply(self, FALSE);
}

static void _redo_callback(dt_action_t *action)
{
  dt_view_t *self = darktable.view_manager->proxy.lens_calib.view;
  if(self) _undo_apply(self, TRUE);
}

const char *name(const dt_view_t *self)
{
  return C_("view", "lensfit");
}

uint32_t view(const dt_view_t *self)
{
  return DT_VIEW_LENS_CALIB;
}

uint32_t flags()
{
  return VIEW_FLAGS_NONE;
}

void init(dt_view_t *self)
{
  self->data = calloc(1, sizeof(dt_lens_calib_t));
  dt_lens_calib_t *d = self->data;

  d->imgid = NO_IMGID;
  d->cells_x = dt_conf_get_int("plugins/lens_calib/cells_x");
  d->cells_y = dt_conf_get_int("plugins/lens_calib/cells_y");
  if(d->cells_x < 1) d->cells_x = DT_LENS_CALIB_DEFAULT_CELLS_X;
  if(d->cells_y < 1) d->cells_y = DT_LENS_CALIB_DEFAULT_CELLS_Y;
  d->cell_aspect = 1.0f;
  d->show_points = TRUE;
  d->show_curves = TRUE;
  d->show_mesh = TRUE;
  d->drag_idx = -1;
  d->manual = g_array_new(FALSE, FALSE, sizeof(dt_lens_calib_manual_point_t));
  d->hover_manual = -1;

  /* No element free func: snapshots are moved between the two stacks by
     pointer, and a stack that freed on removal would free them mid-move. */
  d->undo_stack = g_ptr_array_new();
  d->redo_stack = g_ptr_array_new();
  d->zoom = 1.0f;
  d->pan_x = d->pan_y = 0.5f;
  d->show_residuals = TRUE;
  d->corner_drag = -1;
  dt_lens_warp_init(&d->warp, DT_LENS_WARP_POLY, 4);
  dt_lens_profile_init(&d->work);
  d->sel_entry = -1;
}

void cleanup(dt_view_t *self)
{
  dt_lens_calib_t *d = self->data;

  _save_points(d);
  dt_lens_profile_cleanup(&d->work);

  if(d->surface)
    cairo_surface_destroy(d->surface);
  if(d->have_grid)
    dt_lens_grid_cleanup(&d->grid);
  if(d->manual)
    g_array_free(d->manual, TRUE);

  _undo_clear(d);
  if(d->undo_stack) g_ptr_array_free(d->undo_stack, TRUE);
  if(d->redo_stack) g_ptr_array_free(d->redo_stack, TRUE);

  _drop_flat(d);
  dt_lens_warp_cleanup(&d->warp);
  free(d->solve_pts);
  free(d->res_dx);
  free(d->res_dy);
  free(d->res_mag);

  free(self->data);
  self->data = NULL;
}

static void _drop_flat(dt_lens_calib_t *d)
{
  if(d->flat_surface)
  {
    cairo_surface_destroy(d->flat_surface);
    d->flat_surface = NULL;
  }
  d->flat_valid = FALSE;
}

static void _drop_surface(dt_lens_calib_t *d)
{
  if(d->surface)
  {
    cairo_surface_destroy(d->surface);
    d->surface = NULL;
  }
  d->surf_width = d->surf_height = 0;
  d->req_width = d->req_height = 0;
  // the corrected copy was made from the surface that just went away
  _drop_flat(d);
}

/* Pick the image to calibrate from. Whatever the rest of the ui treats as
   acted on wins -- the same hover-then-selection rule the other views
   follow, so the image you were looking at carries over. Falls back to
   the first image of the collection so entering the view never lands on
   an empty canvas. */
static dt_imgid_t _pick_image(void)
{
  const dt_imgid_t acted_on = dt_act_on_get_main_image();
  if(dt_is_valid_imgid(acted_on))
    return acted_on;

  dt_imgid_t imgid = NO_IMGID;
  sqlite3_stmt *stmt;
  DT_DEBUG_SQLITE3_PREPARE_V2
    (dt_database_get(darktable.db),
     "SELECT imgid FROM memory.collected_images ORDER BY rowid LIMIT 1",
     -1, &stmt, NULL);
  if(sqlite3_step(stmt) == SQLITE_ROW)
    imgid = sqlite3_column_int(stmt, 0);
  sqlite3_finalize(stmt);

  return imgid;
}


/* Persistence of hand-placed points.
 *
 * Written as a small sidecar per image under the user's config directory,
 * keyed by image id, rather than embedded in the XMP. Putting custom data
 * into darktable's XMP would mean extending its metadata schema and the
 * exif writer, which is a much larger and more invasive change; a sidecar
 * gives the same "pick up where I left off" behaviour and is trivially
 * inspectable. The chart geometry travels with the points, because a
 * lattice index only means something relative to the chart it came from.
 *
 * Coordinates are normalized, so the file stays valid regardless of what
 * resolution the points were placed at.
 */
static char *_points_path(const dt_imgid_t imgid)
{
  char cfg[PATH_MAX] = { 0 };
  dt_loc_get_user_config_dir(cfg, sizeof(cfg));
  char *dir = g_build_filename(cfg, "lensprofiles", "points", NULL);
  g_mkdir_with_parents(dir, 0755);
  char *name = g_strdup_printf("%d.txt", (int)imgid);
  char *path = g_build_filename(dir, name, NULL);
  g_free(dir);
  g_free(name);
  return path;
}

static void _save_points(dt_lens_calib_t *d)
{
  if(!dt_is_valid_imgid(d->imgid) || !d->manual) return;

  char *path = _points_path(d->imgid);

  // nothing to remember: drop the file rather than leave an empty one
  if(!d->manual->len)
  {
    g_unlink(path);
    g_free(path);
    return;
  }

  GString *out = g_string_new("# darktable lens calibration points v2\n");
  g_string_append_printf(out, "cells %d %d\n", d->cells_x, d->cells_y);
  g_string_append_printf(out, "aspect %.6f\n", d->cell_aspect);

  for(int c = 0; c < d->corner_count; c++)
    g_string_append_printf(out, "corner %.8f %.8f\n",
                           d->corner_x[c], d->corner_y[c]);

  /* `n` carries the lattice index and whether the node was measured; the
     older `p` form had neither and is still read below. */
  for(guint i = 0; i < d->manual->len; i++)
  {
    const dt_lens_calib_manual_point_t *pt =
      &g_array_index(d->manual, dt_lens_calib_manual_point_t, i);
    g_string_append_printf(out, "n %d %d %d %.8f %.8f\n",
                           pt->col, pt->row, pt->measured ? 1 : 0,
                           pt->x, pt->y);
  }

  GError *err = NULL;
  if(!g_file_set_contents(path, out->str, out->len, &err))
  {
    dt_print(DT_DEBUG_ALWAYS, "[lens_calib] could not save points to `%s': %s",
             path, err ? err->message : "unknown error");
    g_clear_error(&err);
  }

  g_string_free(out, TRUE);
  g_free(path);
}

static void _load_points(dt_lens_calib_t *d)
{
  if(!d->manual) return;
  g_array_set_size(d->manual, 0);
  d->hover_manual = -1;
  d->corner_count = 0;
  d->corner_drag = -1;

  /* Loading is a change of subject, not an edit. Undoing across it would
     restore one image's points onto another, so the history starts over. */
  _undo_clear(d);

  if(!dt_is_valid_imgid(d->imgid)) return;

  char *path = _points_path(d->imgid);
  char *content = NULL;
  if(!g_file_get_contents(path, &content, NULL, NULL))
  {
    g_free(path);
    return;
  }

  gchar **lines = g_strsplit(content, "\n", -1);
  for(int i = 0; lines[i]; i++)
  {
    const char *l = lines[i];
    if(l[0] == '#' || l[0] == '\0') continue;

    int cx, cy, col, row, meas;
    float fx, fy;

    if(sscanf(l, "n %d %d %d %f %f", &col, &row, &meas, &fx, &fy) == 5)
    {
      if(fx >= 0.0f && fx <= 1.0f && fy >= 0.0f && fy <= 1.0f)
      {
        const dt_lens_calib_manual_point_t pt =
          { fx, fy, col, row, meas != 0 };
        g_array_append_val(d->manual, pt);
      }
      continue;
    }

    if(sscanf(l, "corner %f %f", &fx, &fy) == 2)
    {
      if(d->corner_count < 4
         && fx >= 0.0f && fx <= 1.0f && fy >= 0.0f && fy <= 1.0f)
      {
        d->corner_x[d->corner_count] = fx;
        d->corner_y[d->corner_count] = fy;
        d->corner_count++;
      }
      continue;
    }

    if(sscanf(l, "cells %d %d", &cx, &cy) == 2)
    {
      if(cx > 0 && cy > 0)
      {
        d->cells_x = cx;
        d->cells_y = cy;
        // keep the panel and conf in step with what the file says
        dt_conf_set_int("plugins/lens_calib/cells_x", cx);
        dt_conf_set_int("plugins/lens_calib/cells_y", cy);
      }
    }
    else if(sscanf(l, "aspect %f", &fx) == 1)
    {
      if(fx > 0.0f) d->cell_aspect = fx;
    }
    else if(sscanf(l, "p %f %f", &fx, &fy) == 2)
    {
      /* Version 1 sidecar: positions only, no lattice index. Left at -1,
         which makes the solver fall back to inferring the indices by
         clustering the way it always did for these. */
      if(fx >= 0.0f && fx <= 1.0f && fy >= 0.0f && fy <= 1.0f)
      {
        const dt_lens_calib_manual_point_t pt = { fx, fy, -1, -1, TRUE };
        g_array_append_val(d->manual, pt);
      }
    }
  }

  g_strfreev(lines);
  g_free(content);
  g_free(path);

  if(d->manual->len)
    dt_print(DT_DEBUG_ALWAYS, "[lens_calib] restored %u points for image %d",
             d->manual->len, (int)d->imgid);
}

static void _drop_grid(dt_lens_calib_t *d)
{
  if(d->have_grid)
  {
    dt_lens_grid_cleanup(&d->grid);
    d->have_grid = FALSE;
  }
  d->detect_img_w = d->detect_img_h = 0;
}

/* Run grid detection on the loaded image.
 *
 * Renders its own cairo surface rather than reading a float mipmap: the
 * format is unambiguous and the buffer is certain to be populated. The
 * detection surface is deliberately much larger than the on-screen one,
 * since thin chart lines do not survive being downscaled to display size.
 */
static void _detect_grid(dt_view_t *self)
{
  dt_lens_calib_t *d = self->data;
  if(!dt_is_valid_imgid(d->imgid)) return;

  _drop_grid(d);
  // a new set of points invalidates any fit made from the old ones
  _drop_solution(d);

  /* Ask for a much larger surface than the one on screen. At display size
     the thin interior chart lines are downscaled to well under a pixel and
     wash out, leaving only the heavy chart border detectable -- which is
     exactly the wrong half of the image. */
  /* Force the mipmap into existence first.
   *
   * dt_view_image_get_surface is best effort: when the mipmap it needs is
   * not cached it queues generation in the background and hands back a
   * smaller one. That is right for drawing, where the next redraw picks up
   * the better version, and wrong here -- detection is a one-shot
   * measurement and would silently run at the wrong resolution, or be
   * refused outright. Asking for the mipmap with BLOCKING first means the
   * surface request below finds it ready, so the first click works rather
   * than merely warming the cache for the second.
   */
  const dt_mipmap_size_t mip =
    dt_mipmap_cache_get_matching_size(DT_LENS_CALIB_DETECT_SIZE,
                                      DT_LENS_CALIB_DETECT_SIZE);
  dt_mipmap_buffer_t warm;
  dt_mipmap_cache_get(&warm, d->imgid, mip, DT_MIPMAP_BLOCKING, 'r');
  const gboolean got_mip = warm.buf != NULL;
  dt_mipmap_cache_release(&warm);

  if(!got_mip)
  {
    dt_control_log(_("could not decode the image for grid detection"));
    return;
  }

  cairo_surface_t *det = NULL;
  const dt_view_surface_value_t res =
    dt_view_image_get_surface(d->imgid, DT_LENS_CALIB_DETECT_SIZE,
                              DT_LENS_CALIB_DETECT_SIZE, &det, TRUE);

  /* A smaller surface than asked for is still usable: the detector reads
     its own dimensions and scales everything to them. It only means this
     image has no mipmap level as large as we would have liked, which for
     a small file is simply the truth. */
  if(!det || res == DT_VIEW_SURFACE_KO)
  {
    if(det) cairo_surface_destroy(det);
    dt_control_log(_("could not read the image for grid detection"));
    return;
  }

  const int w = cairo_image_surface_get_width(det);
  const int h = cairo_image_surface_get_height(det);
  cairo_surface_flush(det);
  const uint8_t *const pix = cairo_image_surface_get_data(det);
  const int stride = cairo_image_surface_get_stride(det);
  if(!pix || w < 16 || h < 16 || stride < w * 4)
  {
    cairo_surface_destroy(det);
    dt_control_log(_("could not read the image for grid detection"));
    return;
  }

  float *lum = dt_alloc_align_float((size_t)w * h);
  if(!lum)
  {
    cairo_surface_destroy(det);
    return;
  }

  // cairo ARGB32 is premultiplied BGRA in memory on little endian; for a
  // luminance estimate the channel order does not matter
  DT_OMP_FOR()
  for(int y = 0; y < h; y++)
  {
    const uint8_t *row = pix + (size_t)y * stride;
    for(int x = 0; x < w; x++)
    {
      const uint8_t b = row[4 * x + 0];
      const uint8_t g = row[4 * x + 1];
      const uint8_t r = row[4 * x + 2];
      lum[(size_t)y * w + x] = (0.299f * r + 0.587f * g + 0.114f * b) / 255.0f;
    }
  }

  const int cells_x = dt_conf_get_int("plugins/lens_calib/cells_x");
  const int cells_y = dt_conf_get_int("plugins/lens_calib/cells_y");
  d->cells_x = MAX(1, cells_x);
  d->cells_y = MAX(1, cells_y);

  const gboolean ok = dt_lens_grid_detect(lum, w, h, d->cells_x, d->cells_y,
                                          TRUE, &d->grid);
  dt_free_align(lum);
  cairo_surface_destroy(det);

  if(ok)
  {
    d->have_grid = TRUE;
    d->detect_img_w = w;
    d->detect_img_h = h;
    dt_control_log(ngettext("found %d grid intersection",
                            "found %d grid intersections",
                            d->grid.point_count), d->grid.point_count);
  }
  else
  {
    dt_lens_grid_cleanup(&d->grid);
    dt_control_log(_("no usable grid found -- check the chart settings"));
  }

  dt_control_queue_redraw_center();
}

static void _chart_changed(dt_view_t *self)
{
  dt_lens_calib_t *d = self->data;
  d->cells_x = MAX(1, dt_conf_get_int("plugins/lens_calib/cells_x"));
  d->cells_y = MAX(1, dt_conf_get_int("plugins/lens_calib/cells_y"));
  // the previous detection was made for a different lattice
  _drop_grid(d);
  _drop_solution(d);
  dt_control_queue_redraw_center();
}

static gboolean _has_points(dt_view_t *self)
{
  const dt_lens_calib_t *d = self->data;
  if(d->have_grid && d->grid.point_count > 0) return TRUE;
  // hand placed points count too; the solver does not care which it gets
  return d->manual && d->manual->len >= 8;
}

/* ------------------------------------------------------------- mesh */

/* Laying out the lattice from four pinned corners.
 *
 * A flat chart photographed by a pinhole camera maps to the frame by a
 * homography, so four corners determine where every node should be -- up
 * to the lens distortion, which is the whole residual we are here to
 * measure. That makes this an excellent starting point and a useless
 * measurement: the generated nodes say only what a distortion-free lens
 * would have produced, so they are marked unmeasured and the solver leaves
 * them out. What they are for is saving several hundred clicks: place four
 * corners, generate, then correct the nodes that are visibly off.
 */
static gboolean _fit_homography(const double *const src,
                                const double *const dst,
                                const int n,
                                double H[9])
{
  if(n < 4) return FALSE;

  /* Eight unknowns, two equations per point pair, accumulated into normal
     equations so this covers both the exactly determined four point case
     and a least squares fit over a whole lattice. h22 is fixed at 1, the
     usual normalization, and safe here because no chart node maps to
     infinity in a photograph of it. */
  double AtA[64] = { 0 };
  double Atb[8] = { 0 };

  for(int i = 0; i < n; i++)
  {
    const double x = src[2 * i], y = src[2 * i + 1];
    const double u = dst[2 * i], v = dst[2 * i + 1];

    const double r0[8] = { x, y, 1.0, 0.0, 0.0, 0.0, -u * x, -u * y };
    const double r1[8] = { 0.0, 0.0, 0.0, x, y, 1.0, -v * x, -v * y };

    for(int a = 0; a < 8; a++)
    {
      for(int c = 0; c < 8; c++)
        AtA[a * 8 + c] += r0[a] * r0[c] + r1[a] * r1[c];
      Atb[a] += r0[a] * u + r1[a] * v;
    }
  }

  if(!gauss_solve(AtA, Atb, 8)) return FALSE;

  for(int i = 0; i < 8; i++) H[i] = Atb[i];
  H[8] = 1.0;
  return TRUE;
}

static void _apply_homography(const double H[9],
                              const double x,
                              const double y,
                              float *ox,
                              float *oy)
{
  const double w = H[6] * x + H[7] * y + H[8];
  const double iw = (fabs(w) > 1e-12) ? 1.0 / w : 0.0;
  *ox = (float)((H[0] * x + H[1] * y + H[2]) * iw);
  *oy = (float)((H[3] * x + H[4] * y + H[5]) * iw);
}

static gboolean _invert_homography(const double H[9], double Hi[9])
{
  const double a = H[0], b = H[1], c = H[2];
  const double e = H[3], f = H[4], g = H[5];
  const double m = H[6], n = H[7], o = H[8];

  const double c00 = f * o - g * n;
  const double c01 = c * n - b * o;
  const double c02 = b * g - c * f;
  const double det = a * c00 + e * c01 + m * c02;
  if(fabs(det) < 1e-14) return FALSE;

  const double id = 1.0 / det;
  Hi[0] = c00 * id;
  Hi[1] = c01 * id;
  Hi[2] = c02 * id;
  Hi[3] = (g * m - e * o) * id;
  Hi[4] = (a * o - c * m) * id;
  Hi[5] = (c * e - a * g) * id;
  Hi[6] = (e * n - f * m) * id;
  Hi[7] = (b * m - a * n) * id;
  Hi[8] = (a * f - b * e) * id;
  return TRUE;
}

/* Put the four pinned corners into a consistent cyclic order, so the user
   can click them in any order they like. Sorting by angle about their
   centroid gives the cycle; rotating that so the top left comes first
   fixes which corner is which. */
static void _order_corners(dt_lens_calib_t *d)
{
  if(d->corner_count != 4) return;

  double cx = 0.0, cy = 0.0;
  for(int i = 0; i < 4; i++)
  {
    cx += d->corner_x[i];
    cy += d->corner_y[i];
  }
  cx *= 0.25;
  cy *= 0.25;

  int order[4] = { 0, 1, 2, 3 };
  double angle[4];
  for(int i = 0; i < 4; i++)
    angle[i] = atan2(d->corner_y[i] - cy, d->corner_x[i] - cx);

  for(int i = 1; i < 4; i++)
  {
    const int key = order[i];
    int j = i - 1;
    while(j >= 0 && angle[order[j]] > angle[key])
    {
      order[j + 1] = order[j];
      j--;
    }
    order[j + 1] = key;
  }

  // start the cycle at the corner nearest the top left of the frame
  int first = 0;
  double best = 1e30;
  for(int i = 0; i < 4; i++)
  {
    const int k = order[i];
    const double dist = (double)d->corner_x[k] * d->corner_x[k]
                      + (double)d->corner_y[k] * d->corner_y[k];
    if(dist < best)
    {
      best = dist;
      first = i;
    }
  }

  float nx[4], ny[4];
  for(int i = 0; i < 4; i++)
  {
    const int k = order[(first + i) % 4];
    nx[i] = d->corner_x[k];
    ny[i] = d->corner_y[k];
  }

  /* atan2 increases counter clockwise, but image y runs downwards, so the
     cycle above comes out clockwise on screen: top left, top right, bottom
     right, bottom left. That is the order the ideal corners below assume. */
  for(int i = 0; i < 4; i++)
  {
    d->corner_x[i] = nx[i];
    d->corner_y[i] = ny[i];
  }
}

// index of an existing node at this lattice position, or -1
static int _node_at_index(const dt_lens_calib_t *d, const int col, const int row)
{
  for(guint i = 0; i < d->manual->len; i++)
  {
    const dt_lens_calib_manual_point_t *p =
      &g_array_index(d->manual, dt_lens_calib_manual_point_t, i);
    if(p->col == col && p->row == row) return (int)i;
  }
  return -1;
}

static void _mesh_from_corners(dt_view_t *self)
{
  dt_lens_calib_t *d = self->data;

  _undo_push(d, _("lay out mesh"));
  d->undo_hold++;

  /* Prefer a pose fitted to the points over the pinned corners: it uses
     every node rather than four, so the lattice it lays out is better
     placed. The corners are only the way in when there are no points yet. */
  _auto_align(self);
  d->undo_hold--;

  /* Nothing to go on. Rather than refuse and leave the user to find a
     checkbox, ask for the one thing that would answer it -- four corners --
     and pick this up again once they are placed. */
  if(!d->have_pose && d->corner_count != 4)
  {
    _undo_drop_last(d);
    _set_corner_mode(self, TRUE);
    dt_control_log(_("click the four chart corners to say where the chart is"));
    return;
  }

  if(!d->have_pose) _order_corners(d);

  const int nx = d->cells_x + 1;
  const int ny = d->cells_y + 1;

  const double src[8] = { 0, 0, (double)d->cells_x, 0,
                          (double)d->cells_x, (double)d->cells_y,
                          0, (double)d->cells_y };
  const double dst[8] = { d->corner_x[0], d->corner_y[0],
                          d->corner_x[1], d->corner_y[1],
                          d->corner_x[2], d->corner_y[2],
                          d->corner_x[3], d->corner_y[3] };

  double H[9];
  if(d->have_pose)
    memcpy(H, d->pose, sizeof(H));
  else if(!_fit_homography(src, dst, 4, H))
  {
    _undo_drop_last(d);
    dt_control_log(_("those four corners are degenerate --"
                     " they must not be collinear"));
    return;
  }

  int added = 0;
  for(int row = 0; row < ny; row++)
    for(int col = 0; col < nx; col++)
    {
      /* Never overwrite a measured node. The generated lattice is a guess
         and anything already measured is better than it. */
      const int existing = _node_at_index(d, col, row);
      if(existing >= 0)
      {
        dt_lens_calib_manual_point_t *p =
          &g_array_index(d->manual, dt_lens_calib_manual_point_t, existing);
        if(p->measured) continue;

        _apply_homography(H, col, row, &p->x, &p->y);
        continue;
      }

      dt_lens_calib_manual_point_t np;
      memset(&np, 0, sizeof(np));
      np.col = col;
      np.row = row;
      np.measured = FALSE;
      _apply_homography(H, col, row, &np.x, &np.y);

      if(np.x < -0.05f || np.x > 1.05f || np.y < -0.05f || np.y > 1.05f)
        continue; // outside the frame; there is nothing to measure there

      g_array_append_val(d->manual, np);
      added++;
    }

  _save_points(d);
  dt_control_log(ngettext("laid out %d lattice node from the corners",
                          "laid out %d lattice nodes from the corners", added),
                 added);
  dt_control_queue_redraw_center();
}

/* Align the reference grid to the measured nodes.
 *
 * The reference overlay is meant to answer "where would these lines be if
 * the lens were perfect", and as an axis-aligned rectangle over the frame
 * it answered a different and useless question -- it did not even sit on
 * the chart. Fitting the pose puts it where the chart is, and then the gap
 * between the orange grid and the green crosses *is* the distortion, drawn
 * to scale and at full size rather than inferred from a residual number.
 *
 * A homography is the right model and its failure is informative: a planar
 * chart under a distortion-free lens maps to the frame by exactly one, so
 * whatever the pose cannot explain is lens, not geometry. That leftover is
 * reported as `pose_rms`.
 */
/* The corners follow from the pose; they are not something to click.
 *
 * A corner is just the lattice node at an extreme index, so once the pose
 * is known the four of them are simply where it puts (0,0), (cx,0),
 * (cx,cy) and (0,cy). Deriving them from a fit over every node beats
 * clicking four, and not only in effort: a least squares pose over several
 * hundred points is far better conditioned than one pinned from four, where
 * every click error goes straight into the result. Pinning survives only as
 * a bootstrap for the case where there are no indexed points at all.
 */
static void _derive_corners(dt_lens_calib_t *d)
{
  if(!d->have_pose) return;

  const int cx = d->cells_x, cy = d->cells_y;
  const int ix[4] = { 0, cx, cx, 0 };
  const int iy[4] = { 0, 0, cy, cy };

  for(int i = 0; i < 4; i++)
    _apply_homography(d->pose, ix[i], iy[i],
                      &d->corner_x[i], &d->corner_y[i]);

  d->corner_count = 4;
  d->corner_drag = -1;
}

/* Give every measured node its lattice index, reading it off the pose.
 *
 * A measured point with no index is invisible to everything that reasons
 * about the lattice. It cannot join a line, so the solver never sees it;
 * _node_at_index cannot find it, so laying out the mesh drops a second,
 * interpolated node on top of it. That is where a doubled node count comes
 * from, and hand placed points have always arrived unindexed.
 *
 * The pose supplies what is missing: it maps lattice coordinates into the
 * frame, so its inverse reads a point's lattice coordinates back out, and
 * the nearest integer site is the index.
 *
 * The tolerance is the load bearing part. Distortion pushes a node off its
 * ideal position -- that is the signal we are here to measure -- but it does
 * not push it most of the way to the next site. Anything further off than
 * that is either a spurious detection or evidence the pose is wrong, and
 * inventing an index for it would quietly corrupt a whole line. Those stay
 * unindexed and simply take no part.
 */
#define DT_LENS_CALIB_SNAP_TOL 0.4f

static int _index_from_pose(dt_lens_calib_t *d)
{
  double Hi[9];
  if(!d->have_pose || !d->manual) return 0;
  if(!_invert_homography(d->pose, Hi)) return 0;

  const int nx = d->cells_x + 1, ny = d->cells_y + 1;
  const size_t sites = (size_t)nx * ny;

  int *owner = malloc(sizeof(int) * sites);
  float *offset = malloc(sizeof(float) * sites);
  char *taken = calloc(sites, 1);
  if(!owner || !offset || !taken)
  {
    free(owner);
    free(offset);
    free(taken);
    return 0;
  }
  for(size_t i = 0; i < sites; i++)
  {
    owner[i] = -1;
    offset[i] = 0.0f;
  }

  int indexed = 0;
  for(guint i = 0; i < d->manual->len; i++)
  {
    dt_lens_calib_manual_point_t *p =
      &g_array_index(d->manual, dt_lens_calib_manual_point_t, i);
    if(!p->measured) continue;

    p->col = p->row = -1;

    float cf, rf;
    _apply_homography(Hi, p->x, p->y, &cf, &rf);

    const int col = (int)lroundf(cf), row = (int)lroundf(rf);
    if(col < 0 || col >= nx || row < 0 || row >= ny) continue;
    if(fabsf(cf - col) > DT_LENS_CALIB_SNAP_TOL
       || fabsf(rf - row) > DT_LENS_CALIB_SNAP_TOL)
      continue;

    const int site = row * nx + col;
    const float off = (cf - col) * (cf - col) + (rf - row) * (rf - row);

    if(owner[site] >= 0)
    {
      /* Two measurements claiming one node. They cannot both be it, and the
         nearer one is the better claim; the other is a duplicate detection
         and stays out rather than fighting it inside the fit. */
      if(off >= offset[site]) continue;

      dt_lens_calib_manual_point_t *loser =
        &g_array_index(d->manual, dt_lens_calib_manual_point_t, owner[site]);
      loser->col = loser->row = -1;
      indexed--;
    }

    owner[site] = (int)i;
    offset[site] = off;
    taken[site] = 1;
    p->col = col;
    p->row = row;
    indexed++;
  }

  /* Now that the measured nodes hold the sites they actually occupy, any
     interpolated node sharing a site is a guess about something already
     measured, and worse than useless: it doubles the point count and draws
     a second marker next to a real one. Walking backwards keeps the indices
     of the not yet visited entries valid. */
  for(int i = (int)d->manual->len - 1; i >= 0; i--)
  {
    const dt_lens_calib_manual_point_t *p =
      &g_array_index(d->manual, dt_lens_calib_manual_point_t, i);
    if(p->measured) continue;

    if(p->col < 0 || p->col >= nx || p->row < 0 || p->row >= ny)
    {
      g_array_remove_index(d->manual, i);
      continue;
    }

    const int site = p->row * nx + p->col;
    if(taken[site])
    {
      g_array_remove_index(d->manual, i);
      continue;
    }
    taken[site] = 1; // a second interpolated node here would go the same way
  }

  free(owner);
  free(offset);
  free(taken);
  return indexed;
}

static void _align_grid(dt_view_t *self)
{
  dt_lens_calib_t *d = self->data;

  if(!d->manual || d->manual->len < 4)
  {
    dt_control_log(_("place at least four points to align the grid"));
    return;
  }

  /* Aligning is not a read-only operation: it re-indexes every measured
     point, drops interpolated duplicates and moves the corners. */
  _undo_push(d, _("align grid"));

  d->have_pose = FALSE;
  d->pose_rms = 0.0f;

  int n = 0;
  for(guint i = 0; i < d->manual->len; i++)
  {
    const dt_lens_calib_manual_point_t *p =
      &g_array_index(d->manual, dt_lens_calib_manual_point_t, i);
    if(p->measured && p->col >= 0 && p->row >= 0) n++;
  }

  /* Bootstrap. Four corners determine a homography exactly, which is what
     makes pinning them useful and also what makes stopping there wrong: an
     exactly determined fit has no residual, so it reports nothing, and every
     node between the corners is left unconstrained. It is a way in, not an
     answer -- so take the pose and carry on below. */
  if(n < 4)
  {
    if(d->corner_count != 4)
    {
      dt_control_log(_("need at least four indexed points, or four pinned"
                       " corners, to align the grid"));
      return;
    }

    _order_corners(d);

    const double src[8] = { 0, 0, (double)d->cells_x, 0,
                            (double)d->cells_x, (double)d->cells_y,
                            0, (double)d->cells_y };
    const double dst[8] = { d->corner_x[0], d->corner_y[0],
                            d->corner_x[1], d->corner_y[1],
                            d->corner_x[2], d->corner_y[2],
                            d->corner_x[3], d->corner_y[3] };

    if(!_fit_homography(src, dst, 4, d->pose))
    {
      dt_control_log(_("those four corners are degenerate --"
                       " they must not be collinear"));
      return;
    }
    d->have_pose = TRUE;
  }

  /* Alternate between indexing the points from the pose and re-fitting the
     pose over everything now indexed. Each step makes the other better: a
     rough pose is enough to name the sites, and once named there are
     hundreds of correspondences instead of four. It settles in two or three
     rounds because the indices are integers -- there is nowhere left to go
     once no point changes site. */
  double *src = NULL, *dst = NULL;
  int fitted = 0;

  for(int pass = 0; pass < 4; pass++)
  {
    if(d->have_pose)
    {
      const int idx = _index_from_pose(d);
      if(idx < 4) break;
      n = idx;
    }

    double *ns = realloc(src, 2 * (size_t)n * sizeof(double));
    double *nd = realloc(dst, 2 * (size_t)n * sizeof(double));
    if(ns) src = ns;
    if(nd) dst = nd;
    if(!ns || !nd) break;

    int k = 0;
    for(guint i = 0; i < d->manual->len && k < n; i++)
    {
      const dt_lens_calib_manual_point_t *p =
        &g_array_index(d->manual, dt_lens_calib_manual_point_t, i);
      if(!p->measured || p->col < 0 || p->row < 0) continue;

      src[2 * k] = p->col;
      src[2 * k + 1] = p->row;
      dst[2 * k] = p->x;
      dst[2 * k + 1] = p->y;
      k++;
    }
    if(k < 4) break;

    double H[9];
    if(!_fit_homography(src, dst, k, H)) break;

    memcpy(d->pose, H, sizeof(H));
    d->have_pose = TRUE;

    // report the leftover in pixels of the frame, not in normalized units
    double acc = 0.0;
    for(int i = 0; i < k; i++)
    {
      float ex, ey;
      _apply_homography(d->pose, src[2 * i], src[2 * i + 1], &ex, &ey);
      const double dx = (ex - dst[2 * i]) * MAX(1, d->surf_width);
      const double dy = (ey - dst[2 * i + 1]) * MAX(1, d->surf_height);
      acc += dx * dx + dy * dy;
    }
    d->pose_rms = (float)sqrt(acc / k);
    fitted = k;
  }

  free(src);
  free(dst);

  if(!fitted)
  {
    dt_control_log(_("could not align the grid to those points"));
    d->have_pose = FALSE;
    dt_control_queue_redraw_center();
    return;
  }

  // and the corners fall out of the pose, so there is nothing to pin
  _derive_corners(d);
  _save_points(d);

  dt_control_log(_("grid aligned to %d points, %.1f px unexplained"),
                 fitted, d->pose_rms);
  dt_control_queue_redraw_center();
}

/* Relax the unmeasured nodes towards their neighbours' average.
 *
 * This is the honest use of relaxation: it fills gaps, with every measured
 * node held fixed. Applying it to the measured nodes as well would be
 * disastrous in a way that looks like an improvement -- a node differs from
 * its neighbours' average precisely because the lens distorted it, so
 * smoothing measured nodes subtracts the signal being measured, and the
 * residual would get better as the fit got worse.
 */
static void _mesh_relax(dt_view_t *self)
{
  dt_lens_calib_t *d = self->data;

  const int nx = d->cells_x + 1;
  const int ny = d->cells_y + 1;
  if(nx < 3 || ny < 3 || !d->manual->len) return;

  const size_t n = (size_t)nx * ny;
  float *px = calloc(n, sizeof(float));
  float *py = calloc(n, sizeof(float));
  char *state = calloc(n, 1); // 0 absent, 1 free, 2 fixed

  if(!px || !py || !state)
  {
    free(px);
    free(py);
    free(state);
    return;
  }

  int fixed = 0;
  for(guint i = 0; i < d->manual->len; i++)
  {
    const dt_lens_calib_manual_point_t *p =
      &g_array_index(d->manual, dt_lens_calib_manual_point_t, i);
    if(p->col < 0 || p->row < 0 || p->col >= nx || p->row >= ny) continue;

    const size_t k = (size_t)p->row * nx + p->col;
    px[k] = p->x;
    py[k] = p->y;
    state[k] = p->measured ? 2 : 1;
    if(p->measured) fixed++;
  }

  if(fixed < 4)
  {
    dt_control_log(_("need at least four measured nodes to relax the rest"));
    free(px);
    free(py);
    free(state);
    return;
  }

  _undo_push(d, _("relax mesh"));

  /* Seed the absent nodes from the measured ones with a plain affine fit,
     so relaxation starts somewhere sensible rather than at the origin and
     converges in a few dozen sweeps instead of thousands. */
  {
    double M[9] = { 0 }, bx[3] = { 0 }, by[3] = { 0 };
    for(int row = 0; row < ny; row++)
      for(int col = 0; col < nx; col++)
      {
        const size_t k = (size_t)row * nx + col;
        if(state[k] != 2) continue;

        const double v[3] = { col, row, 1.0 };
        for(int a = 0; a < 3; a++)
        {
          for(int c = 0; c < 3; c++) M[a * 3 + c] += v[a] * v[c];
          bx[a] += v[a] * px[k];
          by[a] += v[a] * py[k];
        }
      }

    double Mx[9], My[9];
    memcpy(Mx, M, sizeof(M));
    memcpy(My, M, sizeof(M));

    if(gauss_solve(Mx, bx, 3) && gauss_solve(My, by, 3))
      for(int row = 0; row < ny; row++)
        for(int col = 0; col < nx; col++)
        {
          const size_t k = (size_t)row * nx + col;
          if(state[k] == 2) continue;
          px[k] = (float)(bx[0] * col + bx[1] * row + bx[2]);
          py[k] = (float)(by[0] * col + by[1] * row + by[2]);
          if(state[k] == 0) state[k] = 1;
        }
  }

  // Gauss-Seidel sweeps; the field is smooth so this settles quickly
  for(int iter = 0; iter < 400; iter++)
  {
    double worst = 0.0;

    for(int row = 0; row < ny; row++)
      for(int col = 0; col < nx; col++)
      {
        const size_t k = (size_t)row * nx + col;
        if(state[k] != 1) continue;

        double sx = 0.0, sy = 0.0;
        int cnt = 0;

        const int dc[4] = { -1, 1, 0, 0 };
        const int dr[4] = { 0, 0, -1, 1 };
        for(int q = 0; q < 4; q++)
        {
          const int c2 = col + dc[q], r2 = row + dr[q];
          if(c2 < 0 || c2 >= nx || r2 < 0 || r2 >= ny) continue;
          const size_t k2 = (size_t)r2 * nx + c2;
          if(!state[k2]) continue;
          sx += px[k2];
          sy += py[k2];
          cnt++;
        }

        if(!cnt) continue;

        const float ox = px[k], oy = py[k];
        px[k] = (float)(sx / cnt);
        py[k] = (float)(sy / cnt);
        worst = MAX(worst, MAX(fabs(px[k] - ox), fabs(py[k] - oy)));
      }

    if(worst < 1e-7) break;
  }

  // write back, adding nodes that did not exist before
  int filled = 0;
  for(int row = 0; row < ny; row++)
    for(int col = 0; col < nx; col++)
    {
      const size_t k = (size_t)row * nx + col;
      if(state[k] != 1) continue;
      if(px[k] < -0.05f || px[k] > 1.05f
         || py[k] < -0.05f || py[k] > 1.05f) continue;

      const int existing = _node_at_index(d, col, row);
      if(existing >= 0)
      {
        dt_lens_calib_manual_point_t *p =
          &g_array_index(d->manual, dt_lens_calib_manual_point_t, existing);
        if(p->measured) continue;
        p->x = px[k];
        p->y = py[k];
      }
      else
      {
        dt_lens_calib_manual_point_t np;
        memset(&np, 0, sizeof(np));
        np.col = col;
        np.row = row;
        np.measured = FALSE;
        np.x = px[k];
        np.y = py[k];
        g_array_append_val(d->manual, np);
        filled++;
      }
    }

  free(px);
  free(py);
  free(state);

  _save_points(d);
  dt_control_log(ngettext("relaxed the mesh, filling %d node",
                          "relaxed the mesh, filling %d nodes", filled),
                 filled);
  dt_control_queue_redraw_center();
}

/* Complete the lattice in one action.
 *
 * Laying out nodes and settling them were two buttons, which was two buttons
 * too many: a laid out node sits where a distortion-free lens would have put
 * it, and leaving it there is never what anyone wants -- settling it against
 * the measured neighbours is strictly better and costs nothing. Nobody has a
 * reason to stop in between, so there is no longer a way to.
 */
static void _mesh_fill(dt_view_t *self)
{
  dt_lens_calib_t *d = self->data;

  _undo_push(d, _("fill in mesh"));
  d->undo_hold++;
  _mesh_from_corners(self);
  _mesh_relax(self);
  d->undo_hold--;
}


/* ---------------------------------------------------------- solving */

static void _drop_solution(dt_lens_calib_t *d)
{
  dt_lens_warp_cleanup(&d->warp);
  dt_lens_warp_init(&d->warp, DT_LENS_WARP_POLY, 4);
  d->have_warp = FALSE;
  memset(&d->result, 0, sizeof(d->result));

  // there is nothing left to flatten the image with
  _drop_flat(d);
  d->flat_view = FALSE;

  /* The pose was fitted to whatever points are being discarded, so it goes
     with them; leaving it would draw a reference grid aligned to a lattice
     that no longer exists. */
  d->have_pose = FALSE;
  d->pose_rms = 0.0f;

  free(d->solve_pts);
  free(d->res_dx);
  free(d->res_dy);
  free(d->res_mag);
  d->solve_pts = NULL;
  d->res_dx = d->res_dy = d->res_mag = NULL;
  d->solve_count = 0;
}

/* Split a sorted list of coordinates wherever the gap between neighbours
   exceeds `gap`, and hand back the index of the cluster each value fell
   into. This is how hand placed points acquire lattice indices.
 *
 * The solver only needs to know which points share a line, not where that
 * line sits in the chart: an overall offset in the indices is absorbed by
 * the pose stage. So relative grouping is enough, and the user does not
 * have to click the intersections in any particular order or start from a
 * particular corner.
 */
static void _cluster_1d(const float *const values,
                        const int n,
                        const float gap,
                        int *out)
{
  if(n < 1) return;

  int *order = calloc(n, sizeof(int));
  if(!order)
  {
    for(int i = 0; i < n; i++) out[i] = 0;
    return;
  }

  for(int i = 0; i < n; i++) order[i] = i;

  // insertion sort: n here is the number of clicked points, never large
  for(int i = 1; i < n; i++)
  {
    const int key = order[i];
    int j = i - 1;
    while(j >= 0 && values[order[j]] > values[key])
    {
      order[j + 1] = order[j];
      j--;
    }
    order[j + 1] = key;
  }

  int cluster = 0;
  out[order[0]] = 0;
  for(int i = 1; i < n; i++)
  {
    if(values[order[i]] - values[order[i - 1]] > gap) cluster++;
    out[order[i]] = cluster;
  }

  free(order);
}

/* Gather the points to fit: the detector's if it ran, otherwise the hand
   placed ones. Detected points already carry lattice indices; clicked ones
   have to be grouped into rows and columns first. */
static int _build_solve_points(dt_lens_calib_t *d,
                               dt_lens_solve_point_t **out,
                               int *img_w,
                               int *img_h)
{
  *out = NULL;
  *img_w = *img_h = 0;

  int full_w = 0, full_h = 0;
  if(dt_is_valid_imgid(d->imgid))
  {
    const dt_image_t *img = dt_image_cache_get(d->imgid, 'r');
    if(img)
    {
      full_w = img->width;
      full_h = img->height;
      dt_image_cache_read_release(img);
    }
  }

  if(d->have_grid && d->grid.point_count >= 8)
  {
    const int n = d->grid.point_count;
    dt_lens_solve_point_t *pts = calloc(n, sizeof(dt_lens_solve_point_t));
    if(!pts) return 0;

    /* Detection ran on its own surface; express everything in that frame
       so the residuals come out in the pixels the points were measured
       in rather than in scaled ones. */
    for(int i = 0; i < n; i++)
    {
      pts[i].x = d->grid.points[i].x;
      pts[i].y = d->grid.points[i].y;
      pts[i].col = d->grid.points[i].col;
      pts[i].row = d->grid.points[i].row;
    }

    *out = pts;
    *img_w = d->detect_img_w;
    *img_h = d->detect_img_h;
    return n;
  }

  if(!d->manual || !d->manual->len) return 0;

  // hand placed points are normalized, so any consistent frame will do;
  // the real one keeps the reported residuals in real pixels
  const int w = full_w > 1 ? full_w : 4000;
  const int h = full_h > 1 ? full_h : 3000;

  /* Interpolated nodes are excluded by default, and it is worth being
     blunt about why: they were computed *from* the assumption that the
     lattice is smooth, so feeding them to a fit whose job is to measure
     departures from smoothness is circular. Including them lowers the
     reported residual while making the calibration worse, which is the
     most dangerous shape a bug can take. The option exists because a
     nearly complete mesh with a couple of filled holes is a defensible
     thing to fit, but it is off unless asked for. */
  const gboolean use_interp =
    dt_conf_get_bool("plugins/lens_calib/use_interpolated");

  int usable = 0;
  gboolean have_indices = FALSE;

  for(guint i = 0; i < d->manual->len; i++)
  {
    const dt_lens_calib_manual_point_t *p =
      &g_array_index(d->manual, dt_lens_calib_manual_point_t, i);
    if(!p->measured && !use_interp) continue;
    if(p->col >= 0 && p->row >= 0) have_indices = TRUE;
    usable++;
  }

  if(usable < 8) return 0;

  dt_lens_solve_point_t *pts = calloc(usable, sizeof(dt_lens_solve_point_t));
  float *xs = malloc(sizeof(float) * usable);
  float *ys = malloc(sizeof(float) * usable);
  int *cols = malloc(sizeof(int) * usable);
  int *rows = malloc(sizeof(int) * usable);

  if(!pts || !xs || !ys || !cols || !rows)
  {
    free(pts);
    free(xs);
    free(ys);
    free(cols);
    free(rows);
    return 0;
  }

  int k = 0;
  for(guint i = 0; i < d->manual->len; i++)
  {
    const dt_lens_calib_manual_point_t *p =
      &g_array_index(d->manual, dt_lens_calib_manual_point_t, i);
    if(!p->measured && !use_interp) continue;

    xs[k] = p->x;
    ys[k] = p->y;
    cols[k] = p->col;
    rows[k] = p->row;
    k++;
  }

  if(have_indices)
  {
    /* The nodes know where they belong, so nothing has to be guessed.
       Any that do not are given an index of their own so they form no
       line and simply sit the fit out, rather than being lumped in with
       whichever row happens to be nearest. */
    int spare = 100000;
    for(int i = 0; i < usable; i++)
      if(cols[i] < 0 || rows[i] < 0)
      {
        cols[i] = spare;
        rows[i] = spare;
        spare++;
      }
  }
  else
  {
    /* Version 1 sidecars carry positions only. Half a cell is the natural
       split: two points on the same chart line are never that far apart
       across it, and two on adjacent lines always are. This does assume a
       roughly square on chart, which is why it is the fallback and not the
       method. */
    const float gap_x = 0.5f / (float)MAX(1, d->cells_x);
    const float gap_y = 0.5f / (float)MAX(1, d->cells_y);

    _cluster_1d(xs, usable, gap_x, cols);
    _cluster_1d(ys, usable, gap_y, rows);
  }

  for(int i = 0; i < usable; i++)
  {
    pts[i].x = xs[i] * w;
    pts[i].y = ys[i] * h;
    pts[i].col = cols[i];
    pts[i].row = rows[i];
  }

  free(xs);
  free(ys);
  free(cols);
  free(rows);

  *out = pts;
  *img_w = w;
  *img_h = h;
  return usable;
}

static void _solve(dt_view_t *self)
{
  dt_lens_calib_t *d = self->data;

  _drop_solution(d);

  dt_lens_solve_point_t *pts = NULL;
  int w = 0, h = 0;
  const int n = _build_solve_points(d, &pts, &w, &h);

  if(n < 8 || w < 2 || h < 2)
  {
    free(pts);
    dt_control_log(_("need at least eight points to fit a lens model"));
    return;
  }

  dt_lens_solve_input_t in = { pts, n, w, h, d->cell_aspect };

  dt_lens_solve_options_t opt;
  dt_lens_solve_default_options(&opt);
  opt.kind = (dt_lens_warp_kind_t)
    dt_conf_get_int("plugins/lens_calib/model");
  opt.order = CLAMP(dt_conf_get_int("plugins/lens_calib/order"), 2,
                    DT_LENS_WARP_MAX_ORDER);
  opt.regularization =
    (float)dt_conf_get_float("plugins/lens_calib/regularization");
  opt.solve_centre = dt_conf_get_bool("plugins/lens_calib/solve_centre");
  opt.chart_frontal = dt_conf_get_bool("plugins/lens_calib/chart_frontal");

  /* Zero tells the solver to measure the squeeze instead of being told it. */
  opt.known_squeeze = dt_conf_get_bool("plugins/lens_calib/measure_squeeze")
    ? 0.0f : (float)dt_conf_get_float("plugins/lens_calib/squeeze");

  if(!dt_lens_solve(&in, &opt, &d->warp, &d->result))
  {
    free(pts);
    dt_control_log(_("the fit did not converge -- check the placed points"));
    return;
  }

  d->have_warp = TRUE;
  // any corrected preview belongs to the fit that just got replaced
  _drop_flat(d);
  d->solve_pts = pts;
  d->solve_count = n;
  d->solve_img_w = w;
  d->solve_img_h = h;

  /* Overscan and underscan follow from the warp, so they are worked out here
     rather than at save time -- the readout should say what the fit implies
     before anyone commits it to a profile. */
  dt_lens_warp_measure_scan(&d->warp, w, h);

  d->res_dx = calloc(n, sizeof(float));
  d->res_dy = calloc(n, sizeof(float));
  d->res_mag = calloc(n, sizeof(float));
  if(d->res_dx && d->res_dy && d->res_mag)
    dt_lens_solve_residuals(&in, &d->warp, d->res_dx, d->res_dy, d->res_mag);

  /* A fit that leaves the lines less straight than it found them has not
     converged on anything, whatever its other numbers say. Say so, rather
     than let a plausible looking readout stand behind a broken model. */
  if(d->result.rms_px > d->result.rms_before_px)
    dt_control_log(_("the fit made the lines worse: %.2f -> %.2f px."
                     " the model cannot describe this lens -- try a"
                     " different one, or more regularization"),
                   d->result.rms_before_px, d->result.rms_px);
  else
    dt_control_log(_("fitted %d points: straightness %.2f -> %.2f px,"
                     " squeeze %.3f"),
                   n, d->result.rms_before_px, d->result.rms_px,
                   d->result.squeeze);

  dt_control_queue_redraw_center();
}

/* Measure vignetting from the loaded frame.
 *
 * Reads DT_MIPMAP_FULL rather than the display surface. That buffer is the
 * image as decoded, before the pixelpipe touches it -- for a raw it is the
 * sensor data itself. Nothing else available here is linear: the display
 * surface has been through the whole pipe and a display transfer curve, and
 * a brightness measurement made on gamma encoded data is confidently wrong
 * rather than merely noisy.
 */
static void _fit_vignette(dt_view_t *self)
{
  dt_lens_calib_t *d = self->data;
  if(!dt_is_valid_imgid(d->imgid))
  {
    dt_control_log(_("no image loaded"));
    return;
  }

  dt_mipmap_buffer_t buf;
  dt_mipmap_cache_get(&buf, d->imgid, DT_MIPMAP_FULL, DT_MIPMAP_BLOCKING, 'r');

  if(!buf.buf || buf.width < 32 || buf.height < 32)
  {
    dt_mipmap_cache_release(&buf);
    dt_control_log(_("could not read the image data for vignetting"));
    return;
  }

  /* The buffer is not a float RGB image, and assuming it is reads off the
     end of it. For a plain raw it is uint16, one channel, and at *raw sensor*
     dimensions -- larger than the image, with the rawprepare crop offset in
     img->crop_x/crop_y. All three have to be honoured: the element size or
     the read overruns, the channel count or the stride is wrong, and the crop
     or the falloff is centred on the sensor rather than on the picture. */
  const int raw_w = buf.width, raw_h = buf.height;

  int channels = 0;
  int datatype = 0;
  int off_x = 0, off_y = 0, iw = 0, ih = 0;
  float black = 0.0f;
  float aperture = 0.0f, focal = 0.0f;

  dt_image_t *img = dt_image_cache_get(d->imgid, 'r');
  if(img)
  {
    channels = (int)img->buf_dsc.channels;
    datatype = (int)img->buf_dsc.datatype;
    off_x = img->crop_x;
    off_y = img->crop_y;
    iw = img->width;
    ih = img->height;
    /* The black level is not an offset that can be ignored. A pedestal
       flattens every ratio towards one, so the measured falloff comes out
       shallower than it is -- and by a factor that depends on exposure. */
    black = (float)img->raw_black_level;
    aperture = img->exif_aperture;
    focal = img->exif_focal_length;
    dt_image_cache_read_release(img);
  }

  const gboolean is_u16 = (datatype == TYPE_UINT16);
  const gboolean is_flt = (datatype == TYPE_FLOAT);

  if((channels != 1 && channels != 4) || (!is_u16 && !is_flt))
  {
    dt_mipmap_cache_release(&buf);
    dt_control_log(_("this image format cannot be measured for vignetting"));
    return;
  }

  // fall back to the whole buffer when the crop is not known
  if(iw < 32 || ih < 32 || off_x < 0 || off_y < 0
     || off_x + iw > raw_w || off_y + ih > raw_h)
  {
    off_x = off_y = 0;
    iw = raw_w;
    ih = raw_h;
  }

  float *lum = malloc(sizeof(float) * (size_t)iw * ih);
  if(!lum)
  {
    dt_mipmap_cache_release(&buf);
    return;
  }

  const uint16_t *const s16 = (const uint16_t *)buf.buf;
  const float *const s32 = (const float *)buf.buf;
  const size_t rstride = (size_t)raw_w * channels;

  DT_OMP_FOR()
  for(int y = 0; y < ih; y++)
  {
    const size_t srow = (size_t)(y + off_y) * rstride;
    for(int x = 0; x < iw; x++)
    {
      const size_t si = srow + (size_t)(x + off_x) * channels;

      float value;
      if(channels == 1)
      {
        // mosaiced sensor data: one sample per photosite, already linear
        value = is_u16 ? (float)s16[si] : s32[si];
        value -= black;
      }
      else
      {
        if(is_u16)
          value = 0.25f * ((float)s16[si] + 2.0f * (float)s16[si + 1]
                           + (float)s16[si + 2]);
        else
          value = 0.25f * (s32[si] + 2.0f * s32[si + 1] + s32[si + 2]);
      }

      lum[(size_t)y * iw + x] = MAX(0.0f, value);
    }
  }

  dt_mipmap_cache_release(&buf);

  const int w = iw, h = ih;

  /* Centre the falloff on the optical centre when one has been fitted. It is
     the same centre the distortion is expanded about -- both are properties
     of where the axis actually crosses the sensor. */
  const float cx = d->have_warp ? d->warp.cx : 0.0f;
  const float cy = d->have_warp ? d->warp.cy : 0.0f;
  const float hint = d->have_warp ? d->warp.squeeze
    : (float)dt_conf_get_float("plugins/lens_calib/squeeze");

  /* Sample the middle of each chart cell, when there is a lattice to say
     where the middles are.
   *
   * Far better than binning the whole frame by radius. That has to guess
   * which pixels are the surface and which are markings on it -- and the
   * guess fails worst at the edges, where the chart runs out and the darker
   * wall behind it gets read as falloff. A cell centre is paper by
   * construction: no ink to reject, nothing outside the chart included, and
   * one clean reading per cell instead of a percentile over a mixture.
   */
  const double hdi = 0.5 * hypot((double)iw, (double)ih);
  float *uv = NULL, *vals = NULL;
  int ns = 0;

  if(d->manual && d->manual->len > 8)
  {
    const int cnx = d->cells_x, cny = d->cells_y;
    uv = malloc(sizeof(float) * 2 * (size_t)cnx * cny);
    vals = malloc(sizeof(float) * (size_t)cnx * cny);

    if(uv && vals)
      for(int r = 0; r < cny; r++)
        for(int c = 0; c < cnx; c++)
        {
          const int ia = _node_at_index(d, c, r);
          const int ib = _node_at_index(d, c + 1, r);
          const int ic = _node_at_index(d, c, r + 1);
          const int id = _node_at_index(d, c + 1, r + 1);
          if(ia < 0 || ib < 0 || ic < 0 || id < 0) continue;

          const dt_lens_calib_manual_point_t *p[4] = {
            &g_array_index(d->manual, dt_lens_calib_manual_point_t, ia),
            &g_array_index(d->manual, dt_lens_calib_manual_point_t, ib),
            &g_array_index(d->manual, dt_lens_calib_manual_point_t, ic),
            &g_array_index(d->manual, dt_lens_calib_manual_point_t, id)
          };

          double mx = 0.0, my = 0.0;
          for(int i = 0; i < 4; i++)
          {
            mx += p[i]->x;
            my += p[i]->y;
          }
          mx *= 0.25;
          my *= 0.25;

          /* Keep well inside the cell. The patch has to miss the lines on
             every side, and the cells nearest the frame edge are the
             narrowest once the lens has finished with them. */
          const double cw = fabs((double)p[1]->x - p[0]->x) * iw;
          const double ch = fabs((double)p[2]->y - p[0]->y) * ih;
          const int rad = (int)(0.30 * MIN(cw, ch));
          if(rad < 1) continue;

          const int px = (int)(mx * iw), py = (int)(my * ih);
          if(px - rad < 0 || py - rad < 0
             || px + rad >= iw || py + rad >= ih) continue;

          double acc = 0.0;
          int m = 0;
          for(int yy = py - rad; yy <= py + rad; yy++)
            for(int xx = px - rad; xx <= px + rad; xx++)
            {
              acc += lum[(size_t)yy * iw + xx];
              m++;
            }
          if(!m) continue;

          uv[2 * ns] = (float)((mx - 0.5) * iw / hdi);
          uv[2 * ns + 1] = (float)((my - 0.5) * ih / hdi);
          vals[ns] = (float)(acc / m);
          ns++;
        }
  }

  float k[3], ex = 1.0f, rms = 0.0f, stops = 0.0f, maxr = 1.0f;
  gboolean ok;

  const gboolean by_cells = (ns >= 24);
  if(by_cells)
    ok = dt_lens_vignette_fit_points(uv, vals, ns, w, h, cx, cy, hint,
                                     k, &ex, &rms, &maxr, &stops);
  else
    ok = dt_lens_vignette_fit(lum, w, h, w, cx, cy, hint,
                              k, &ex, &rms, &stops);

  free(uv);
  free(vals);
  free(lum);

  if(!ok)
  {
    dt_control_log(_("could not measure vignetting from this frame"));
    return;
  }

  memcpy(d->warp.vig_k, k, sizeof(d->warp.vig_k));
  d->warp.vig_ex = ex;
  d->warp.vig_aperture = aperture;
  d->warp.have_vig = TRUE;
  if(d->warp.focal <= 0.0f) d->warp.focal = focal;

  d->vig_rms = rms;
  d->vig_stops = stops;
  d->vig_max_r = maxr;
  d->vig_samples = by_cells ? ns : 0;

  /* Say when the answer is not believable.
   *
   * Real lenses lose one to three stops in the corners; four and beyond is
   * almost always the subject rather than the lens. The usual cause is a
   * frame whose corners are not the evenly lit surface at all -- a chart that
   * does not reach them, so the outermost radius bins measure the darker wall
   * behind it and the fit reads that as falloff. A large residual says the
   * same thing: the brightness is not a function of radius alone.
   *
   * This is reported rather than corrected because there is nothing to
   * correct. The measurement faithfully describes the photograph it was given.
   */
  if(stops > 3.5f || rms > 0.02f)
    dt_control_log(_("vignetting measured as %.2f stops, residual %.3f --"
                     " this is too much to be the lens. the corners of this"
                     " frame are probably not the evenly lit surface;"
                     " use a frame filled edge to edge"),
                   stops, rms);
  else if(by_cells)
    dt_control_log(_("vignetting: %.2f stops at the corner, from %d cells"
                     " out to %.0f%% of the corner radius,"
                     " ellipticity %.2f, residual %.3f"),
                   stops, ns, 100.0 * maxr, ex, rms);
  else
    dt_control_log(_("vignetting: %.2f stops at the corner, ellipticity %.2f,"
                     " residual %.3f"),
                   stops, ex, rms);

  dt_control_queue_redraw_center();
}

/* Take ownership of whatever the detector produced.
 *
 * The two systems are deliberately separate up to this point -- a
 * re-detect must never silently discard hand placed work -- so this is
 * the one crossing, and it goes one way. Detected points become ordinary
 * editable ones and the detection is dropped, because having the solver
 * read half-corrected points from two places would be worse than either.
 */
static void _adopt_detected(dt_lens_calib_t *d)
{
  if(!d->have_grid || d->grid.point_count < 1) return;
  if(d->detect_img_w < 1 || d->detect_img_h < 1) return;

  _undo_push(d, _("adopt detected points"));

  const double inv_w = 1.0 / (double)d->detect_img_w;
  const double inv_h = 1.0 / (double)d->detect_img_h;

  // a point already placed by hand near a detected one wins; it was put
  // there on purpose
  const float near2 = 0.004f * 0.004f;
  int added = 0;

  for(int i = 0; i < d->grid.point_count; i++)
  {
    /* The detector's indices come across with the positions. They are
       lattice snapped and windowed to the declared chart size, so they mean
       something -- and carrying them is what spares the solver from having
       to guess the grouping by clustering positions. */
    const dt_lens_calib_manual_point_t np =
      { (float)(d->grid.points[i].x * inv_w),
        (float)(d->grid.points[i].y * inv_h),
        d->grid.points[i].col,
        d->grid.points[i].row,
        TRUE };

    gboolean duplicate = FALSE;
    for(guint k = 0; k < d->manual->len && !duplicate; k++)
    {
      const dt_lens_calib_manual_point_t *p =
        &g_array_index(d->manual, dt_lens_calib_manual_point_t, k);

      // same lattice node, or close enough on screen to be the same thing
      if(p->col == np.col && p->row == np.row)
      {
        duplicate = TRUE;
        break;
      }
      const float dx = p->x - np.x, dy = p->y - np.y;
      duplicate = (dx * dx + dy * dy) < near2;
    }

    if(!duplicate)
    {
      g_array_append_val(d->manual, np);
      added++;
    }
  }

  _drop_grid(d);
  _save_points(d);

  dt_control_log(ngettext("adopted %d detected point for editing",
                          "adopted %d detected points for editing", added),
                 added);
}

/* Fit the pose whenever the points change, rather than waiting to be asked.
   The reference grid is misleading until it is aligned -- an unaligned one
   sits over the frame rather than the chart -- so leaving it to a button
   press means the default state of the display is the wrong one. */
static void _auto_align(dt_view_t *self)
{
  dt_lens_calib_t *d = self->data;
  if(!d->manual) return;

  int indexed = 0;
  for(guint i = 0; i < d->manual->len && indexed < 4; i++)
  {
    const dt_lens_calib_manual_point_t *p =
      &g_array_index(d->manual, dt_lens_calib_manual_point_t, i);
    if(p->measured && p->col >= 0 && p->row >= 0) indexed++;
  }

  /* Four pinned corners are enough to get started even with nothing indexed:
     _align_grid uses them to bootstrap and then indexes the points itself. */
  if(indexed >= 4 || d->corner_count == 4) _align_grid(self);
}

static void _set_manual_edit(dt_view_t *self, const gboolean on)
{
  dt_lens_calib_t *d = self->data;

  if(on && !d->manual_edit)
  {
    _adopt_detected(d);
    // the adopted points bring their indices, so the pose can be had at once
    d->undo_hold++;
    _auto_align(self);
    d->undo_hold--;
  }

  d->manual_edit = on;
  d->drag_idx = -1;
  d->hover_manual = -1;
  dt_control_queue_redraw_center();
}

static gboolean _get_manual_edit(dt_view_t *self)
{
  const dt_lens_calib_t *d = self->data;
  return d->manual_edit;
}

static void _set_flat(dt_view_t *self, const gboolean on)
{
  dt_lens_calib_t *d = self->data;

  if(on && !d->have_warp)
  {
    dt_control_log(_("fit a lens model first -- there is nothing to"
                     " flatten the image with yet"));
    return;
  }

  d->flat_view = on;
  dt_control_queue_redraw_center();
}

static gboolean _get_flat(dt_view_t *self)
{
  const dt_lens_calib_t *d = self->data;
  return d->flat_view;
}

static void _set_falloff(dt_view_t *self, const gboolean on)
{
  dt_lens_calib_t *d = self->data;

  if(on && !d->warp.have_vig)
  {
    dt_control_log(_("measure vignetting first -- there is nothing to"
                     " even out yet"));
    return;
  }

  d->falloff_view = on;
  // the cached surface was built with the other setting
  _drop_flat(d);
  dt_control_queue_redraw_center();
}

static gboolean _get_falloff(dt_view_t *self)
{
  const dt_lens_calib_t *d = self->data;
  return d->falloff_view;
}

static void _set_show(dt_view_t *self, const int what, const gboolean on)
{
  dt_lens_calib_t *d = self->data;

  switch(what)
  {
    case DT_LENS_CALIB_SHOW_POINTS:    d->show_points = on; break;
    case DT_LENS_CALIB_SHOW_CURVES:    d->show_curves = on; break;
    case DT_LENS_CALIB_SHOW_RESIDUALS: d->show_residuals = on; break;
    case DT_LENS_CALIB_SHOW_MESH:      d->show_mesh = on; break;
    default: return;
  }

  dt_control_queue_redraw_center();
}

static gboolean _get_show(dt_view_t *self, const int what)
{
  const dt_lens_calib_t *d = self->data;

  switch(what)
  {
    case DT_LENS_CALIB_SHOW_POINTS:    return d->show_points;
    case DT_LENS_CALIB_SHOW_CURVES:    return d->show_curves;
    case DT_LENS_CALIB_SHOW_RESIDUALS: return d->show_residuals;
    case DT_LENS_CALIB_SHOW_MESH:      return d->show_mesh;
    default: return FALSE;
  }
}

static int _point_count(dt_view_t *self)
{
  const dt_lens_calib_t *d = self->data;
  return d->manual ? (int)d->manual->len : 0;
}

static int _measured_count(dt_view_t *self)
{
  const dt_lens_calib_t *d = self->data;
  if(!d->manual) return 0;

  int n = 0;
  for(guint i = 0; i < d->manual->len; i++)
    if(g_array_index(d->manual, dt_lens_calib_manual_point_t, i).measured)
      n++;
  return n;
}

static int _interpolated_count(dt_view_t *self)
{
  const dt_lens_calib_t *d = self->data;
  if(!d->manual) return 0;

  int n = 0;
  for(guint i = 0; i < d->manual->len; i++)
    if(!g_array_index(d->manual, dt_lens_calib_manual_point_t, i).measured)
      n++;
  return n;
}

static int _stray_count(dt_view_t *self)
{
  const dt_lens_calib_t *d = self->data;
  if(!d->manual) return 0;

  int n = 0;
  for(guint i = 0; i < d->manual->len; i++)
  {
    const dt_lens_calib_manual_point_t *p =
      &g_array_index(d->manual, dt_lens_calib_manual_point_t, i);
    if(p->measured && (p->col < 0 || p->row < 0)) n++;
  }
  return n;
}

/* Corner mode is not a setting. It is a question the view asks when it has no
   other way to find the chart, and stops asking the moment it does. */
static void _set_corner_mode(dt_view_t *self, const gboolean on)
{
  dt_lens_calib_t *d = self->data;
  d->corner_mode = on;
  d->corner_drag = -1;
  d->hover_manual = -1;
  dt_control_queue_redraw_center();
}

static gboolean _get_corner_mode(dt_view_t *self)
{
  const dt_lens_calib_t *d = self->data;
  return d->corner_mode;
}

static int _corner_count(dt_view_t *self)
{
  const dt_lens_calib_t *d = self->data;
  return d->corner_count;
}

/* Does this overlay have anything to draw?
 *
 * Several of them are empty most of the time -- the traced curves cease to
 * exist the moment the points are adopted, residuals need a fit, interpolated
 * nodes need gaps -- and a ticked checkbox over an empty layer is a lie about
 * what is on screen. The panel greys these out instead.
 */
static gboolean _layer_has_data(dt_view_t *self, const int what)
{
  const dt_lens_calib_t *d = self->data;

  switch(what)
  {
    case DT_LENS_CALIB_SHOW_POINTS:
    case DT_LENS_CALIB_SHOW_MESH:
      return d->manual && d->manual->len > 0;

    case DT_LENS_CALIB_SHOW_CURVES:
      return d->have_grid && d->grid.point_count > 0;

    case DT_LENS_CALIB_SHOW_RESIDUALS:
      return d->have_warp && d->solve_count > 0;

    default:
      return FALSE;
  }
}

static gboolean _has_pose(dt_view_t *self)
{
  const dt_lens_calib_t *d = self->data;
  return d->have_pose;
}

static float _pose_rms(dt_view_t *self)
{
  const dt_lens_calib_t *d = self->data;
  return d->pose_rms;
}

static void _clear_points(dt_view_t *self)
{
  dt_lens_calib_t *d = self->data;
  if(!d->manual || !d->manual->len) return;

  /* The most destructive button in the panel, and the one that most needs a
     way back -- several hundred hand corrections live in here. */
  _undo_push(d, _("clear points"));

  g_array_set_size(d->manual, 0);
  d->hover_manual = -1;
  d->drag_idx = -1;
  _drop_solution(d);
  _save_points(d);
  dt_control_queue_redraw_center();
}

static gboolean _has_solution(dt_view_t *self)
{
  const dt_lens_calib_t *d = self->data;
  return d->have_warp;
}

static char *_status_text(dt_view_t *self)
{
  const dt_lens_calib_t *d = self->data;
  if(!d->have_warp && !d->warp.have_vig) return NULL;

  GString *s = g_string_new(NULL);

  /* The affine residual is reported because it is the health check on the
     squeeze: a chart shot at an angle is foreshortened, which looks
     exactly like an anamorphic squeeze and would be silently folded into
     the number if nobody said so. */
  if(d->have_warp)
    g_string_append_printf
      (s, _("%d points, %d lines\n"
            "straightness %.2f -> %.2f px\n"
            "squeeze %.4f, chart says %.4f\n"
            "(grid fit %.2f px)\n"
            "%d iterations"),
       d->result.used_points, d->result.lines_used,
       d->result.rms_before_px, d->result.rms_px,
       d->result.squeeze, d->result.squeeze_measured,
       d->result.affine_rms_px,
       d->result.iterations);

  if(d->have_warp)
    g_string_append_printf
      (s, _("\noverscan %.3f, underscan %.3f"),
       d->warp.overscan, d->warp.underscan);

  if(d->warp.have_vig)
  {
    if(s->len) g_string_append_c(s, '\n');
    g_string_append_printf(s, _("vignetting %.2f stops at f/%.1f"),
                           d->vig_stops, d->warp.vig_aperture);
    if(d->vig_rms > 0.0f)
      g_string_append_printf(s, _(" (residual %.3f)"), d->vig_rms);

    if(d->vig_samples)
      g_string_append_printf(s, _("\n  from %d cells, out to %.0f%% of the"
                                  " corner"),
                             d->vig_samples, 100.0 * d->vig_max_r);

    // the readout has to carry the doubt, not just the number
    if(d->vig_stops > 3.5f || d->vig_rms > 0.02f)
      g_string_append(s, _("\n  -- implausible, check the frame"));
    else if(d->vig_samples && d->vig_max_r < 0.9f)
      g_string_append(s, _("\n  -- corner figure is extrapolated"));
  }

  return g_string_free(s, FALSE);
}

static gboolean _has_vignette(dt_view_t *self)
{
  const dt_lens_calib_t *d = self->data;
  return d->warp.have_vig;
}


/* Adopt an existing lensfit profile into this session.
 *
 * Opening a profile is not the same as starting from nothing: it brings the
 * warp in as the current fit so it can be inspected, corrected, or extended
 * with another focal length. What it must not do is quietly relabel someone
 * else's numbers as this session's measurement, so the profile's own
 * provenance comes along with it.
 */
static gboolean _open_profile(dt_view_t *self, const char *name)
{
  dt_lens_calib_t *d = self->data;
  if(!name || !*name) return FALSE;

  gchar *path = dt_lens_profile_find(name);
  if(!path) return FALSE;

  dt_lens_profile_t prof;
  const gboolean loaded = dt_lens_profile_load(&prof, path);
  g_free(path);

  if(!loaded || !prof.warps || !prof.warps->len)
  {
    dt_lens_profile_cleanup(&prof);
    dt_control_log(_("could not open `%s'"), name);
    return FALSE;
  }

  /* The entry closest to the focal length being worked at, so opening a zoom
     profile while set up at 50mm does not hand back the 16mm measurement. */
  const float want = dt_conf_get_float("plugins/lens_calib/focal");
  guint best = 0;
  if(want > 0.0f)
  {
    double bestd = G_MAXDOUBLE;
    for(guint i = 0; i < prof.warps->len; i++)
    {
      const dt_lens_warp_t *w = &g_array_index(prof.warps, dt_lens_warp_t, i);
      const double dd = fabs(log(MAX(w->focal, 1e-3) / want));
      if(dd < bestd) { bestd = dd; best = i; }
    }
  }

  const dt_lens_warp_t *pick = &g_array_index(prof.warps, dt_lens_warp_t, best);

  dt_lens_warp_t copy;
  if(!dt_lens_warp_copy(&copy, pick))
  {
    dt_lens_profile_cleanup(&prof);
    return FALSE;
  }

  _undo_push(d, _("open profile"));

  dt_lens_warp_cleanup(&d->warp);
  d->warp = copy;
  d->have_warp = TRUE;

  /* The whole profile comes in, so every focal length and aperture it holds
     stays visible and editable -- opening a zoom must not silently reduce it
     to the one entry being looked at. */
  dt_lens_profile_cleanup(&d->work);
  dt_lens_profile_init(&d->work);
  for(guint i = 0; i < prof.warps->len; i++)
  {
    dt_lens_warp_t c2;
    if(dt_lens_warp_copy(&c2, &g_array_index(prof.warps, dt_lens_warp_t, i)))
      g_array_append_val(d->work.warps, c2);
  }
  d->sel_entry = (int)best;

  dt_conf_set_string("plugins/lens_calib/profile_name", name);
  dt_conf_set_string("plugins/lens_calib/maker", prof.maker);
  dt_conf_set_string("plugins/lens_calib/model", prof.model);
  dt_conf_set_string("plugins/lens_calib/mount", prof.mount);
  dt_conf_set_float("plugins/lens_calib/crop_factor", prof.crop_factor);
  if(pick->focal > 0.0f)
    dt_conf_set_float("plugins/lens_calib/focal", pick->focal);
  if(pick->aperture > 0.0f)
    dt_conf_set_float("plugins/lens_calib/aperture", pick->aperture);
  dt_conf_set_int("plugins/lens_calib/source", (int)prof.source);
  dt_conf_set_string("plugins/lens_calib/license", prof.license);
  dt_conf_set_string("plugins/lens_calib/parent", prof.parent);

  /* The saved range becomes the override from here on, whether it was
     originally typed or derived -- reopening a profile must not make an
     already-recorded 24-70mm range start silently tracking whatever subset
     happens to be measured again in this session. */
  if(prof.focal_min > 0.0f)
    dt_conf_set_float("plugins/lens_calib/focal_min", prof.focal_min);
  if(prof.focal_max > 0.0f)
    dt_conf_set_float("plugins/lens_calib/focal_max", prof.focal_max);
  if(prof.aperture_min > 0.0f)
    dt_conf_set_float("plugins/lens_calib/aperture_min", prof.aperture_min);
  if(prof.aperture_max > 0.0f)
    dt_conf_set_float("plugins/lens_calib/aperture_max", prof.aperture_max);
  if(prof.distance_min > 0.0f)
    dt_conf_set_float("plugins/lens_calib/distance_min", prof.distance_min);
  if(prof.distance_max > 0.0f)
    dt_conf_set_float("plugins/lens_calib/distance_max", prof.distance_max);

  dt_lens_profile_cleanup(&prof);

  dt_control_queue_redraw_center();
  dt_lens_calib_refresh_panels();
  dt_control_log(_("opened `%s'"), name);
  return TRUE;
}

/* Start a fresh profile: forget the one currently open or being built, so
 * the next save writes a new file rather than merging into whatever this
 * session last opened or imported.
 *
 * Deliberately narrower than _drop_solution(): this clears the profile
 * identity and the warp collection, but leaves the chart detection, the
 * points and the pose alone. Those describe the calibration photograph
 * currently loaded, not which lens is being profiled from it -- clearing
 * them here would throw away a working detection for no reason connected
 * to what "new profile" is asking for.
 */
/* The lens's own range along one axis, as observed across every
 * measurement the session currently holds -- not what has been typed
 * anywhere, just what the data says.
 *
 * axis: 0 focal (from any entry with geometry), 1 aperture, 2 focus
 * distance (both from any entry with a vignetting measurement, since
 * that is the only kind that records either). Returns FALSE if nothing
 * relevant has a known value for that axis yet.
 *
 * Shared between the UI proxy (so the lens panel can show it) and the save
 * path (so a value nobody typed still ends up in the saved profile) --
 * one definition of "the observed range" rather than two that could drift.
 */
static gboolean _axis_range(const dt_lens_calib_t *d, const int axis,
                            float *lo, float *hi)
{
  if(!d->work.warps) return FALSE;

  float mn = FLT_MAX, mx = -FLT_MAX;
  gboolean found = FALSE;

  for(guint i = 0; i < d->work.warps->len; i++)
  {
    const dt_lens_warp_t *w = &g_array_index(d->work.warps, dt_lens_warp_t, i);
    float v = 0.0f;

    if(axis == 0)
    {
      if(!w->have_geometry || w->focal <= 0.0f) continue;
      v = w->focal;
    }
    else if(axis == 1)
    {
      if(!w->have_vig) continue;
      v = w->vig_aperture > 0.0f ? w->vig_aperture : w->aperture;
      if(v <= 0.0f) continue;
    }
    else
    {
      if(!w->have_vig || w->focus_distance <= 0.0f) continue;
      v = w->focus_distance;
    }

    found = TRUE;
    mn = MIN(mn, v);
    mx = MAX(mx, v);
  }

  if(!found) return FALSE;
  if(lo) *lo = mn;
  if(hi) *hi = mx;
  return TRUE;
}

static gboolean _measurement_range(dt_view_t *self, const int axis,
                                   float *lo, float *hi)
{
  const dt_lens_calib_t *d = self->data;
  return _axis_range(d, axis, lo, hi);
}

static gboolean _new_profile(dt_view_t *self)
{
  dt_lens_calib_t *d = self->data;

  _undo_push(d, _("new profile"));

  dt_lens_warp_cleanup(&d->warp);
  dt_lens_warp_init(&d->warp, DT_LENS_WARP_POLY, 4);
  d->have_warp = FALSE;

  dt_lens_profile_cleanup(&d->work);
  dt_lens_profile_init(&d->work);
  d->sel_entry = -1;

  dt_conf_set_string("plugins/lens_calib/profile_name", "");
  dt_conf_set_string("plugins/lens_calib/maker", "");
  dt_conf_set_string("plugins/lens_calib/model", "");
  dt_conf_set_string("plugins/lens_calib/mount", "");
  dt_conf_set_float("plugins/lens_calib/crop_factor", 1.0f);
  dt_conf_set_float("plugins/lens_calib/focal", 0.0f);
  dt_conf_set_float("plugins/lens_calib/aperture", 0.0f);
  dt_conf_set_int("plugins/lens_calib/source", (int)DT_LENS_SOURCE_MEASURED);
  dt_conf_set_string("plugins/lens_calib/license", "");
  dt_conf_set_string("plugins/lens_calib/parent", "");
  dt_conf_set_float("plugins/lens_calib/focal_min", 0.0f);
  dt_conf_set_float("plugins/lens_calib/focal_max", 0.0f);
  dt_conf_set_float("plugins/lens_calib/aperture_min", 0.0f);
  dt_conf_set_float("plugins/lens_calib/aperture_max", 0.0f);
  dt_conf_set_float("plugins/lens_calib/distance_min", 0.0f);
  dt_conf_set_float("plugins/lens_calib/distance_max", 0.0f);

  dt_control_queue_redraw_center();
  dt_lens_calib_refresh_panels();
  dt_control_log(_("new profile"));
  return TRUE;
}

/* Convert one Lensfun distortion entry into our radial polynomial.
 *
 * Lensfun's three models are all polynomials in the undistorted radius that
 * return the distorted one, which is the direction DT_LENS_WARP_RADIAL_POLY
 * already runs in -- so this is a rearrangement of coefficients rather than a
 * fit, and it is exact. The same mapping is used by the lensfit repo's
 * converter, where it was checked against Lensfun's own evaluator across the
 * whole database to 6.7e-16.
 */
static gboolean _lf_dist_to_warp(const dt_lf_dist_entry_t *e,
                                 dt_lens_warp_t *w)
{
  dt_lens_warp_init(w, DT_LENS_WARP_RADIAL_POLY, 2);
  if(w->nparams < 6) return FALSE;

  for(int i = 0; i < 6; i++) w->p[i] = 0.0f;
  w->p[0] = 1.0f;   // ellipticity: Lensfun has no anisotropic term

  switch(e->model)
  {
    case DT_LF_DIST_POLY3:
      w->p[1] = 1.0f - e->terms[0];
      w->p[3] = e->terms[0];
      break;
    case DT_LF_DIST_POLY5:
      w->p[1] = 1.0f;
      w->p[3] = e->terms[0];
      w->p[5] = e->terms[1];
      break;
    case DT_LF_DIST_PTLENS:
      w->p[1] = 1.0f - e->terms[0] - e->terms[1] - e->terms[2];
      w->p[2] = e->terms[2];
      w->p[3] = e->terms[1];
      w->p[4] = e->terms[0];
      break;
    default:
      return FALSE;
  }

  w->focal = e->focal;
  w->have_geometry = TRUE;
  return TRUE;
}

/* Bring a lens in from the Lensfun catalogue as a starting point.
 *
 * This is a convenience, not a measurement: the result is somebody else's
 * calibration expressed in our model, and it is marked as such so that
 * saving it cannot claim otherwise and sharing it carries Lensfun's licence
 * forward.
 */
/* Find or make the entry for one operating point.
 *
 * Mirrors _v2_slot in the profile reader on purpose. Lensfun indexes its
 * three tables independently -- distortion by focal length, vignetting by
 * focal length *and* aperture *and* focus distance -- so a lens with
 * distortion at five focal lengths and vignetting at twelve combinations is
 * not five entries, nor twelve, but however many distinct operating points
 * the union of those tables covers. Most end up carrying only one kind of
 * data, which is exactly what the format allows.
 */
static int _work_slot(dt_lens_calib_t *d,
                      const float focal,
                      const float aperture,
                      const float distance)
{
  for(guint i = 0; i < d->work.warps->len; i++)
  {
    const dt_lens_warp_t *w =
      &g_array_index(d->work.warps, dt_lens_warp_t, i);
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
  g_array_append_val(d->work.warps, w);
  return (int)d->work.warps->len - 1;
}

/* Bring a lens in from the Lensfun catalogue as a starting point.
 *
 * The whole lens, not one operating point of it: a zoom calibrated at five
 * focal lengths arrives as five distortion entries, and its vignetting
 * arrives at every focal length and aperture Lensfun measured. Taking only
 * the nearest match -- which this did at first -- silently threw away most
 * of what the catalogue knew and left a profile that looked complete.
 *
 * The result is somebody else's calibration expressed in our model, and it
 * is marked as such so that saving it cannot claim otherwise and sharing it
 * carries Lensfun's licence forward.
 */
static gboolean _import_lensfun(dt_view_t *self, const char *model)
{
  dt_lens_calib_t *d = self->data;
  if(!model || !*model) return FALSE;

  dt_lf_profile_t lf;
  dt_lf_profile_init(&lf);
  if(!dt_lf_load_profile(model, &lf))
  {
    dt_lf_profile_cleanup(&lf);
    dt_control_log(_("could not read that lens from Lensfun"));
    return FALSE;
  }

  const guint n_dist = lf.dist ? lf.dist->len : 0;
  const guint n_tca = lf.tca ? lf.tca->len : 0;
  const guint n_vign = lf.vign ? lf.vign->len : 0;

  if(!n_dist && !n_tca && !n_vign)
  {
    dt_lf_profile_cleanup(&lf);
    dt_control_log(_("that lens carries no calibration data"));
    return FALSE;
  }

  _undo_push(d, _("import from Lensfun"));

  dt_lens_profile_cleanup(&d->work);
  dt_lens_profile_init(&d->work);
  d->sel_entry = -1;

  int unconverted = 0;

  /* Distortion and TCA are indexed by focal length alone -- neither depends
     on aperture -- so they share a slot at aperture 0, meaning "any". */
  for(guint i = 0; i < n_dist; i++)
  {
    const dt_lf_dist_entry_t *e =
      &g_array_index(lf.dist, dt_lf_dist_entry_t, i);

    dt_lens_warp_t conv;
    if(!_lf_dist_to_warp(e, &conv))
    {
      unconverted++;
      continue;
    }

    const int k = _work_slot(d, e->focal, 0.0f, 0.0f);
    dt_lens_warp_t *slot = &g_array_index(d->work.warps, dt_lens_warp_t, k);

    const float ap = slot->aperture, fd = slot->focus_distance;
    dt_lens_warp_cleanup(slot);
    dt_lens_warp_copy(slot, &conv);
    slot->aperture = ap;
    slot->focus_distance = fd;

    dt_lens_warp_cleanup(&conv);
  }

  for(guint i = 0; i < n_tca; i++)
  {
    const dt_lf_tca_entry_t *e =
      &g_array_index(lf.tca, dt_lf_tca_entry_t, i);

    const int k = _work_slot(d, e->focal, 0.0f, 0.0f);
    dt_lens_warp_t *slot = &g_array_index(d->work.warps, dt_lens_warp_t, k);

    /* Our TCA scale is b r^2 + c r + v per channel. Lensfun's poly3 is
       already in that shape; its linear model is the constant term alone. */
    if(e->model == DT_LF_TCA_POLY3)
    {
      slot->tca_r[0] = e->terms[0];
      slot->tca_r[1] = e->terms[1];
      slot->tca_r[2] = e->terms[2];
      slot->tca_b[0] = e->terms[3];
      slot->tca_b[1] = e->terms[4];
      slot->tca_b[2] = e->terms[5];
    }
    else if(e->model == DT_LF_TCA_LINEAR)
    {
      slot->tca_r[0] = slot->tca_r[1] = 0.0f;
      slot->tca_r[2] = e->terms[0];
      slot->tca_b[0] = slot->tca_b[1] = 0.0f;
      slot->tca_b[2] = e->terms[1];
    }
    else
      unconverted++;
  }

  /* Vignetting keeps its own aperture and focus distance, which is why it
     gets slots of its own rather than being folded onto the distortion
     entries: the same focal length at two stops is two measurements. */
  for(guint i = 0; i < n_vign; i++)
  {
    const dt_lf_vign_entry_t *e =
      &g_array_index(lf.vign, dt_lf_vign_entry_t, i);

    const int k = _work_slot(d, e->focal, e->aperture, e->distance);
    dt_lens_warp_t *slot = &g_array_index(d->work.warps, dt_lens_warp_t, k);

    slot->have_vig = TRUE;
    slot->vig_k[0] = e->k[0];
    slot->vig_k[1] = e->k[1];
    slot->vig_k[2] = e->k[2];
    slot->vig_ex = 1.0f;
    /* Lensfun's "pa" polynomial is a transmission that divides, not a gain
       that multiplies. Since 1/t is not a polynomial the coefficients cannot
       be converted, only labelled -- and our evaluator honours the label. */
    slot->vig_convention = DT_LENS_VIG_TRANSMISSION;
    slot->vig_aperture = e->aperture;
  }

  if(!d->work.warps->len)
  {
    dt_lf_profile_cleanup(&lf);
    dt_control_log(_("nothing in that lens could be converted"));
    return FALSE;
  }

  /* Land on the entry nearest whatever focal length is set, so importing a
     zoom while working at 50mm does not open at its widest end. */
  const float want = dt_conf_get_float("plugins/lens_calib/focal");
  int best = 0;
  double bestd = G_MAXDOUBLE;
  for(guint i = 0; i < d->work.warps->len; i++)
  {
    const dt_lens_warp_t *w =
      &g_array_index(d->work.warps, dt_lens_warp_t, i);
    const double dd = (want > 0.0f)
      ? fabs(log(MAX(w->focal, 1e-3) / want)) : (double)i;
    if(dd < bestd) { bestd = dd; best = (int)i; }
  }

  dt_lens_warp_t copy;
  if(dt_lens_warp_copy(&copy, &g_array_index(d->work.warps,
                                             dt_lens_warp_t, best)))
  {
    dt_lens_warp_cleanup(&d->warp);
    d->warp = copy;
    d->have_warp = TRUE;
    d->sel_entry = best;

    if(d->warp.focal > 0.0f)
      dt_conf_set_float("plugins/lens_calib/focal", d->warp.focal);
    if(d->warp.aperture > 0.0f)
      dt_conf_set_float("plugins/lens_calib/aperture", d->warp.aperture);
  }

  dt_conf_set_string("plugins/lens_calib/maker", lf.maker);
  dt_conf_set_string("plugins/lens_calib/model", lf.model);
  dt_conf_set_string("plugins/lens_calib/mount", lf.mount);
  if(lf.crop_factor > 0.0f)
    dt_conf_set_float("plugins/lens_calib/crop_factor", lf.crop_factor);

  dt_conf_set_int("plugins/lens_calib/source", (int)DT_LENS_SOURCE_LENSFUN);
  dt_conf_set_string("plugins/lens_calib/license", "CC-BY-SA-3.0");
  dt_conf_set_string("plugins/lens_calib/parent", model);

  const guint entries = d->work.warps->len;
  dt_lf_profile_cleanup(&lf);

  d->flat_valid = FALSE;
  dt_control_queue_redraw_center();
  dt_lens_calib_refresh_panels();

  if(unconverted)
    dt_control_log(_("imported %u entries, %d could not be converted"),
                   entries, unconverted);
  else
    dt_control_log(_("imported %u entries from Lensfun"), entries);

  return TRUE;
}


/* The fitted numbers, as a flat editable list.
 *
 * The panel needs to show whatever the current model actually has -- a
 * radial polynomial has different coefficients from a spline, and vignetting
 * or TCA may not be present at all -- so rather than the panel knowing the
 * shape of every model, the view enumerates what exists and the panel draws
 * one row per entry.
 */
typedef struct _val_t
{
  const char *name;
  float *ptr;
} _val_t;

#define DT_LENS_CALIB_MAX_VALUES 32

static int _values(dt_lens_calib_t *d, _val_t *out)
{
  int n = 0;
  if(!d->have_warp) return 0;

  /* Name each coefficient by what it means in the model that is actually
     loaded. A column of rows all called "p" is not a readout, and it is
     precisely the case where a hand edit needs to know which number it is
     touching. */
  static const char *const radial_poly[] =
    { "ellipticity", "c1", "c2", "c3", "c4", "c5" };
  static const char *const anam[] =
    { "ellipticity", "k1", "k2", "k3", "k4" };
  static const char *const generic[] =
    { "p1", "p2", "p3", "p4", "p5", "p6", "p7", "p8",
      "p9", "p10", "p11", "p12", "p13", "p14", "p15", "p16" };

  for(int i = 0; i < d->warp.nparams && n < DT_LENS_CALIB_MAX_VALUES; i++)
  {
    if(d->warp.kind == DT_LENS_WARP_RADIAL_POLY && i < 6)
      out[n].name = radial_poly[i];
    else if(d->warp.kind == DT_LENS_WARP_ANAM_RADIAL && i < 5)
      out[n].name = anam[i];
    else
      out[n].name = generic[MIN(i, 15)];
    out[n].ptr = &d->warp.p[i];
    n++;
  }

  if(n + 2 <= DT_LENS_CALIB_MAX_VALUES)
  {
    out[n].name = "centre x"; out[n].ptr = &d->warp.cx; n++;
    out[n].name = "centre y"; out[n].ptr = &d->warp.cy; n++;
  }

  if(d->warp.have_vig && n + 4 <= DT_LENS_CALIB_MAX_VALUES)
  {
    out[n].name = "vignette k1"; out[n].ptr = &d->warp.vig_k[0]; n++;
    out[n].name = "vignette k2"; out[n].ptr = &d->warp.vig_k[1]; n++;
    out[n].name = "vignette k3"; out[n].ptr = &d->warp.vig_k[2]; n++;
    out[n].name = "vignette exponent"; out[n].ptr = &d->warp.vig_ex; n++;
  }

  if(n + 6 <= DT_LENS_CALIB_MAX_VALUES)
  {
    static const char *const tn[] = { "red b", "red c", "red v",
                                      "blue b", "blue c", "blue v" };
    for(int i = 0; i < 3; i++)
    {
      out[n].name = tn[i]; out[n].ptr = &d->warp.tca_r[i]; n++;
    }
    for(int i = 0; i < 3; i++)
    {
      out[n].name = tn[3 + i]; out[n].ptr = &d->warp.tca_b[i]; n++;
    }
  }

  return n;
}

static int _value_count(dt_view_t *self)
{
  _val_t v[DT_LENS_CALIB_MAX_VALUES];
  return _values(self->data, v);
}

static const char *_value_name(dt_view_t *self, const int i)
{
  _val_t v[DT_LENS_CALIB_MAX_VALUES];
  const int n = _values(self->data, v);
  return (i >= 0 && i < n) ? v[i].name : "";
}

static double _value_get(dt_view_t *self, const int i)
{
  _val_t v[DT_LENS_CALIB_MAX_VALUES];
  const int n = _values(self->data, v);
  return (i >= 0 && i < n) ? *v[i].ptr : 0.0;
}

/* Editing a coefficient stops the profile being a measurement.
 *
 * The quality figures a submission carries -- the straightness residual, the
 * point count -- describe the fit that produced these numbers. Change one by
 * hand and those figures now describe something that no longer exists, so
 * the profile has to stop claiming it was measured. It may still be a better
 * profile than the fit was; it is just no longer evidence.
 */
static void _value_set(dt_view_t *self, const int i, const double val)
{
  dt_lens_calib_t *d = self->data;
  _val_t v[DT_LENS_CALIB_MAX_VALUES];
  const int n = _values(d, v);
  if(i < 0 || i >= n) return;

  if(*v[i].ptr == (float)val) return;

  _undo_push(d, _("edit value"));
  *v[i].ptr = (float)val;

  const int src = dt_conf_get_int("plugins/lens_calib/source");
  if(src == DT_LENS_SOURCE_MEASURED)
    dt_conf_set_int("plugins/lens_calib/source", DT_LENS_SOURCE_EDITED);

  d->flat_valid = FALSE;
  dt_control_queue_redraw_center();
  dt_lens_calib_refresh_panels();
}


/* ------------------------------------------------------ profile entries */

/* Fold the fit in progress back into the collection.
 *
 * `warp` is a working copy of one entry, so anything done to it -- fitting,
 * measuring vignetting, editing a coefficient by hand -- has to be written
 * back before the selection moves or the profile is saved. Doing it at those
 * two moments rather than on every change keeps the live fit cheap. */
static void _entry_sync_out(dt_lens_calib_t *d)
{
  if(!d->have_warp || !d->work.warps) return;

  d->warp.focal = dt_conf_get_float("plugins/lens_calib/focal");
  d->warp.aperture = dt_conf_get_float("plugins/lens_calib/aperture");

  if(d->sel_entry >= 0 && d->sel_entry < (int)d->work.warps->len)
  {
    dt_lens_warp_t *slot =
      &g_array_index(d->work.warps, dt_lens_warp_t, d->sel_entry);
    dt_lens_warp_cleanup(slot);
    dt_lens_warp_copy(slot, &d->warp);
  }
  else
  {
    dt_lens_warp_t copy;
    if(dt_lens_warp_copy(&copy, &d->warp))
    {
      g_array_append_val(d->work.warps, copy);
      d->sel_entry = (int)d->work.warps->len - 1;
    }
  }
}

static int _entry_count(dt_view_t *self)
{
  dt_lens_calib_t *d = self->data;
  return d->work.warps ? (int)d->work.warps->len : 0;
}

static int _entry_selected(dt_view_t *self)
{
  dt_lens_calib_t *d = self->data;
  return d->sel_entry;
}

/* A line the panel can put in a list: what the entry was shot at, and what
   it actually carries. "35mm" alone would not distinguish a geometry fit
   from a vignetting measurement made at the same focal length. */
static char *_entry_label(dt_view_t *self, const int i)
{
  dt_lens_calib_t *d = self->data;
  if(!d->work.warps || i < 0 || i >= (int)d->work.warps->len)
    return g_strdup("");

  const dt_lens_warp_t *w = &g_array_index(d->work.warps, dt_lens_warp_t, i);

  GString *s = g_string_new(NULL);
  if(w->focal > 0.0f) g_string_append_printf(s, "%gmm", (double)w->focal);
  else                g_string_append(s, "? mm");

  if(w->aperture > 0.0f)
    g_string_append_printf(s, " f/%g", (double)w->aperture);

  g_string_append(s, "  ");
  if(w->have_geometry) g_string_append(s, "geometry ");
  if(w->have_vig) g_string_append(s, "vignette ");
  if(!w->have_geometry && !w->have_vig) g_string_append(s, "empty");

  return g_string_free(s, FALSE);
}

static void _append_sorted_unique(GArray *out, const float v)
{
  for(guint i = 0; i < out->len; i++)
    if(fabsf(g_array_index(out, float, i) - v) < 1e-3f) return;

  guint pos = out->len;
  for(guint i = 0; i < out->len; i++)
    if(g_array_index(out, float, i) > v) { pos = i; break; }
  g_array_insert_val(out, pos, v);
}

/* One line saying what a profile already holds, so the lens panel can show
   it without opening the fit panel and stepping through the cascading
   selectors -- "2 focal lengths, geometry at 24-70mm; vignetting at f/2.8,
   f/4, f/5.6" rather than making the user find out by clicking through. */
static char *_contents_summary(dt_view_t *self)
{
  dt_lens_calib_t *d = self->data;
  if(!d->work.warps || !d->work.warps->len)
    return g_strdup(_("no measurements yet"));

  GArray *geo_focals = g_array_new(FALSE, FALSE, sizeof(float));
  GArray *vig_apertures = g_array_new(FALSE, FALSE, sizeof(float));

  for(guint i = 0; i < d->work.warps->len; i++)
  {
    const dt_lens_warp_t *w = &g_array_index(d->work.warps, dt_lens_warp_t, i);
    if(w->have_geometry) _append_sorted_unique(geo_focals, w->focal);
    if(w->have_vig)
      _append_sorted_unique(vig_apertures,
                            w->vig_aperture > 0.0f ? w->vig_aperture
                                                    : w->aperture);
  }

  GString *s = g_string_new(NULL);

  if(geo_focals->len == 1)
  {
    const float f = g_array_index(geo_focals, float, 0);
    if(f > 0.0f) g_string_append_printf(s, _("geometry at %gmm"), (double)f);
    else         g_string_append(s, _("geometry, no focal length recorded"));
  }
  else if(geo_focals->len > 1)
  {
    const float lo = g_array_index(geo_focals, float, 0);
    const float hi = g_array_index(geo_focals, float, geo_focals->len - 1);
    if(lo > 0.0f)
      g_string_append_printf(s, _("%u focal lengths (%g\xe2\x80\x93%gmm)"),
                             geo_focals->len, (double)lo, (double)hi);
    else
      g_string_append_printf(s, _("%u focal lengths"), geo_focals->len);
  }

  if(vig_apertures->len)
  {
    if(s->len) g_string_append(s, _("; "));
    g_string_append(s, _("vignetting at "));
    /* A zoom's vignetting can carry a couple dozen apertures (the Nikon
       80-400 imports ~19 of them) -- naming them all would defeat the
       point of a one line summary. */
    const guint show = MIN(vig_apertures->len, 4);
    for(guint i = 0; i < show; i++)
    {
      if(i) g_string_append(s, ", ");
      const float a = g_array_index(vig_apertures, float, i);
      if(a > 0.0f) g_string_append_printf(s, "f/%g", (double)a);
      else         g_string_append(s, _("unrecorded aperture"));
    }
    if(vig_apertures->len > show)
      g_string_append_printf(s, _(" and %u more"), vig_apertures->len - show);
  }

  if(!s->len) g_string_append(s, _("no measurements yet"));

  g_array_free(geo_focals, TRUE);
  g_array_free(vig_apertures, TRUE);
  return g_string_free(s, FALSE);
}

static void _entry_select(dt_view_t *self, const int i)
{
  dt_lens_calib_t *d = self->data;
  if(!d->work.warps || i < 0 || i >= (int)d->work.warps->len) return;
  if(i == d->sel_entry) return;

  _entry_sync_out(d);

  const dt_lens_warp_t *pick =
    &g_array_index(d->work.warps, dt_lens_warp_t, i);

  dt_lens_warp_t copy;
  if(!dt_lens_warp_copy(&copy, pick)) return;

  dt_lens_warp_cleanup(&d->warp);
  d->warp = copy;
  d->have_warp = TRUE;
  d->sel_entry = i;

  if(pick->focal > 0.0f)
    dt_conf_set_float("plugins/lens_calib/focal", pick->focal);
  if(pick->aperture > 0.0f)
    dt_conf_set_float("plugins/lens_calib/aperture", pick->aperture);

  d->flat_valid = FALSE;
  dt_control_queue_redraw_center();
  dt_lens_calib_refresh_panels();
}

/* Start another entry at the focal length and aperture now set.
 *
 * Seeded from the current fit rather than from nothing: calibrating the next
 * stop of the same lens starts much closer to the previous stop than to a
 * blank model, and anything wrong with the seed gets overwritten by the next
 * fit anyway. */
static void _entry_add(dt_view_t *self)
{
  dt_lens_calib_t *d = self->data;
  if(!d->work.warps) return;

  _entry_sync_out(d);

  dt_lens_warp_t copy;
  if(d->have_warp)
  {
    if(!dt_lens_warp_copy(&copy, &d->warp)) return;
  }
  else
    dt_lens_warp_init(&copy, DT_LENS_WARP_POLY, 4);

  copy.focal = dt_conf_get_float("plugins/lens_calib/focal");
  copy.aperture = dt_conf_get_float("plugins/lens_calib/aperture");

  g_array_append_val(d->work.warps, copy);
  d->sel_entry = (int)d->work.warps->len - 1;

  dt_lens_calib_refresh_panels();
  dt_control_log(_("added an entry at %g mm"), (double)copy.focal);
}

static void _entry_remove(dt_view_t *self, const int i)
{
  dt_lens_calib_t *d = self->data;
  if(!d->work.warps || i < 0 || i >= (int)d->work.warps->len) return;

  dt_lens_warp_t *slot = &g_array_index(d->work.warps, dt_lens_warp_t, i);
  dt_lens_warp_cleanup(slot);
  g_array_remove_index(d->work.warps, i);

  if(d->sel_entry >= (int)d->work.warps->len)
    d->sel_entry = (int)d->work.warps->len - 1;

  dt_lens_calib_refresh_panels();
}

static int _cmp_float(const void *a, const void *b)
{
  const float x = *(const float *)a, y = *(const float *)b;
  return (x > y) - (x < y);
}

/* Distinct values along one axis of the (focal, aperture, distance) key.
 *
 * The entries are genuinely three dimensional -- a Lensfun zoom can carry
 * six focal lengths, nineteen apertures at each and four focus distances at
 * each of those, which is four hundred odd entries. Listing that flat is
 * unusable, and grouping by focal length alone still leaves seventy rows a
 * group. So the panel asks one axis at a time, each filtered by the axes
 * above it, and the three short lists reconstruct the key between them.
 *
 * `axis` is 0 focal, 1 aperture, 2 distance. Returns how many were written,
 * sorted ascending.
 */
static int _axis_values(dt_view_t *self,
                        const int axis,
                        const float focal,
                        const float aperture,
                        float *out,
                        const int max)
{
  dt_lens_calib_t *d = self->data;
  int n = 0;

  if(!d->work.warps || !out || max <= 0) return 0;

  for(guint i = 0; i < d->work.warps->len; i++)
  {
    const dt_lens_warp_t *w =
      &g_array_index(d->work.warps, dt_lens_warp_t, i);

    if(axis >= 1 && fabsf(w->focal - focal) > 1e-3f) continue;
    if(axis >= 2 && fabsf(w->aperture - aperture) > 1e-3f) continue;

    const float v = (axis == 0) ? w->focal
                  : (axis == 1) ? w->aperture
                                : w->focus_distance;

    gboolean seen = FALSE;
    for(int k = 0; k < n; k++)
      if(fabsf(out[k] - v) < 1e-3f) { seen = TRUE; break; }

    if(!seen && n < max) out[n++] = v;
  }

  qsort(out, n, sizeof(float), _cmp_float);
  return n;
}

/* The entry at one operating point, or -1 if there is none. */
static int _entry_find(dt_view_t *self,
                       const float focal,
                       const float aperture,
                       const float distance)
{
  dt_lens_calib_t *d = self->data;
  if(!d->work.warps) return -1;

  for(guint i = 0; i < d->work.warps->len; i++)
  {
    const dt_lens_warp_t *w =
      &g_array_index(d->work.warps, dt_lens_warp_t, i);
    if(fabsf(w->focal - focal) < 1e-3f
       && fabsf(w->aperture - aperture) < 1e-3f
       && fabsf(w->focus_distance - distance) < 1e-3f)
      return (int)i;
  }
  return -1;
}

static gboolean _save_profile(dt_view_t *self, const char *name)
{
  dt_lens_calib_t *d = self->data;
  if(!name || !*name) return FALSE;

  /* Vignetting alone is a profile worth saving. It is measured from a
     different frame and by different means than the distortion, and needing
     a geometry fit first would mean no way to record it for a lens whose
     chart shot never worked out. */
  if(!d->have_warp && !d->warp.have_vig) return FALSE;

  gchar *path = dt_lens_profile_path(name);
  if(!path) return FALSE;

  /* Merge into an existing profile of the same name rather than replacing
     it, so calibrating a zoom is a matter of repeating the measurement at
     each focal length. A second measurement at a focal length already
     present replaces just that one. */
  dt_lens_profile_t prof;
  if(!dt_lens_profile_load(&prof, path))
  {
    dt_lens_profile_cleanup(&prof);
    dt_lens_profile_init(&prof);
  }

  g_strlcpy(prof.name, name, sizeof(prof.name));

  int full_w = d->solve_img_w, full_h = d->solve_img_h;
  float focal = 0.0f;

  float aperture = 0.0f;

  if(dt_is_valid_imgid(d->imgid))
  {
    const dt_image_t *img = dt_image_cache_get(d->imgid, 'r');
    if(img)
    {
      if(img->width > 1) full_w = img->width;
      if(img->height > 1) full_h = img->height;
      focal = img->exif_focal_length;
      aperture = img->exif_aperture;
      dt_image_cache_read_release(img);
    }
  }

  /* What the lens panel says wins over EXIF, because for the lenses this
     tool exists to serve EXIF is usually empty: a manual cine lens has no
     electronic contacts, so the camera records focal length 0 and aperture
     0. A profile carrying a focal length of zero cannot be matched to any
     image later, and is rejected outright by the database -- so the typed
     value is the authoritative one whenever there is one. */
  const float ui_focal = dt_conf_get_float("plugins/lens_calib/focal");
  const float ui_aperture = dt_conf_get_float("plugins/lens_calib/aperture");
  if(ui_focal > 0.0f) focal = ui_focal;
  if(ui_aperture > 0.0f) aperture = ui_aperture;

  gchar *ui_maker = dt_conf_get_string("plugins/lens_calib/maker");
  gchar *ui_model = dt_conf_get_string("plugins/lens_calib/model");
  gchar *ui_mount = dt_conf_get_string("plugins/lens_calib/mount");

  if(ui_maker && *ui_maker) g_strlcpy(prof.maker, ui_maker, sizeof(prof.maker));
  if(ui_model && *ui_model) g_strlcpy(prof.model, ui_model, sizeof(prof.model));
  if(ui_mount && *ui_mount) g_strlcpy(prof.mount, ui_mount, sizeof(prof.mount));

  g_free(ui_maker);
  g_free(ui_model);
  g_free(ui_mount);

  const float ui_crop = dt_conf_get_float("plugins/lens_calib/crop_factor");
  if(ui_crop > 0.0f) prof.crop_factor = ui_crop;

  /* Provenance follows whatever the session actually is. It is set when a
     profile is opened, when one is imported, and when a value is hand
     edited -- saving must not overwrite that with an optimistic default. */
  prof.source =
    (dt_lens_source_t)dt_conf_get_int("plugins/lens_calib/source");

  gchar *lic = dt_conf_get_string("plugins/lens_calib/license");
  gchar *par = dt_conf_get_string("plugins/lens_calib/parent");
  if(lic) g_strlcpy(prof.license, lic, sizeof(prof.license));
  if(par) g_strlcpy(prof.parent, par, sizeof(prof.parent));
  g_free(lic);
  g_free(par);

  /* The lens's own range: whatever was typed as an override, or else
     whatever the measurements so far actually span. Either way the saved
     profile carries a real number rather than forcing a reader to
     recompute it from the warp array -- and a value nobody ever measures
     the far end of (a zoom calibrated only at 24mm and 35mm of a
     24-70mm) still gets its typed range saved rather than silently
     narrowed to what happens to exist. */
  {
    struct { float *out; const char *conf; int axis; gboolean want_max; }
    ranges[] = {
      { &prof.focal_min,    "plugins/lens_calib/focal_min",    0, FALSE },
      { &prof.focal_max,    "plugins/lens_calib/focal_max",    0, TRUE  },
      { &prof.aperture_min, "plugins/lens_calib/aperture_min", 1, FALSE },
      { &prof.aperture_max, "plugins/lens_calib/aperture_max", 1, TRUE  },
      { &prof.distance_min, "plugins/lens_calib/distance_min", 2, FALSE },
      { &prof.distance_max, "plugins/lens_calib/distance_max", 2, TRUE  },
    };
    for(size_t i = 0; i < sizeof(ranges) / sizeof(ranges[0]); i++)
    {
      const float manual = dt_conf_get_float(ranges[i].conf);
      if(manual > 0.0f)
      {
        *ranges[i].out = manual;
        continue;
      }
      float lo = 0.0f, hi = 0.0f;
      if(_axis_range(d, ranges[i].axis, &lo, &hi))
        *ranges[i].out = ranges[i].want_max ? hi : lo;
    }
  }

  prof.width = full_w;
  prof.height = full_h;

  /* Write every entry the session holds, not just the live one. A zoom
     calibrated across three focal lengths is three entries, and saving only
     whichever was last selected would quietly drop the other two. */
  _entry_sync_out(d);

  gboolean added = FALSE;

  if(d->work.warps && d->work.warps->len)
  {
    for(guint i = 0; i < d->work.warps->len; i++)
    {
      const dt_lens_warp_t *src =
        &g_array_index(d->work.warps, dt_lens_warp_t, i);
      dt_lens_warp_t w;
      if(!dt_lens_warp_copy(&w, src)) continue;
      if(w.focal <= 0.0f) w.focal = focal;
      if(w.aperture <= 0.0f && aperture > 0.0f) w.aperture = aperture;
      if(dt_lens_profile_add(&prof, &w)) added = TRUE;
      dt_lens_warp_cleanup(&w);
    }
  }
  else
  {
    dt_lens_warp_t w;
    if(!dt_lens_warp_copy(&w, &d->warp))
    {
      dt_lens_profile_cleanup(&prof);
      g_free(path);
      return FALSE;
    }
    w.focal = focal;
    if(aperture > 0.0f) w.aperture = aperture;
    added = dt_lens_profile_add(&prof, &w);
    dt_lens_warp_cleanup(&w);
  }

  GError *err = NULL;
  const gboolean ok = added && dt_lens_profile_save(&prof, path, &err);

  if(ok)
    dt_control_log(_("saved calibration profile `%s'"), name);
  else
    dt_control_log(_("could not save `%s': %s"), name,
                   err ? err->message : "?");

  if(err) g_error_free(err);
  dt_lens_profile_cleanup(&prof);

  /* Contribute the measurement, if the user has left sharing on. This runs
     after the local save has already succeeded and its result is ignored on
     purpose: calibrating a lens is the user's work and uploading it is ours,
     so ours must not be able to fail theirs. */
  if(ok && dt_lens_share_enabled())
  {
    const int measured = _measured_count(self);

    dt_lens_share_quality_t q = { 0 };
    q.have_geometry = d->have_warp;
    q.have_vignetting = d->warp.have_vig;
    q.straightness_before_px = d->result.rms_before_px;
    q.straightness_after_px = d->result.rms_px;
    q.points_measured = measured;
    q.points_stray = _stray_count(self);
    q.vig_coverage = d->vig_samples ? d->vig_max_r : 0.0;

    GError *serr = NULL;
    if(dt_lens_share_submit(path, &q, &serr))
      dt_control_log(_("`%s' queued for the lensfit database"), name);
    else
      dt_print(DT_DEBUG_ALWAYS, "[lensfit] could not queue `%s': %s",
               name, serr ? serr->message : "?");
    if(serr) g_error_free(serr);
  }

  g_free(path);

  return ok;
}

static gboolean _export_stmap(dt_view_t *self,
                              const char *path,
                              const gboolean bottom_up)
{
  dt_lens_calib_t *d = self->data;
  if(!d->have_warp || !path) return FALSE;

  int w = d->solve_img_w, h = d->solve_img_h;
  if(dt_is_valid_imgid(d->imgid))
  {
    const dt_image_t *img = dt_image_cache_get(d->imgid, 'r');
    if(img)
    {
      /* The map is written at the frame's own resolution: it is a lookup
         per output pixel, so anything smaller would have to be resampled
         by the consumer and would lose the sub-pixel accuracy the whole
         calibration exists to provide. */
      if(img->width > 1) w = img->width;
      if(img->height > 1) h = img->height;
      dt_image_cache_read_release(img);
    }
  }

  GError *err = NULL;
  const gboolean ok =
    dt_lens_stmap_write(&d->warp, w, h, bottom_up, path, &err);

  if(ok)
    dt_control_log(_("wrote a %dx%d STmap"), w, h);
  else
    dt_control_log(_("could not write the STmap: %s"),
                   err ? err->message : "?");

  if(err) g_error_free(err);
  return ok;
}

/* Export the session's current profile as a Lensfun database fragment.
 *
 * Unlike saving, this never reads or merges with a profile already on disk
 * -- it writes exactly what the session holds right now, once, to wherever
 * the caller asked. Whatever does not fit Lensfun's model is reported in
 * `unrepresentable` rather than silently dropped or approximated (see
 * dt_lens_profile_export_lensfun for what that can be); the caller owns the
 * strings appended to it.
 */
static gboolean _export_lensfun(dt_view_t *self,
                                const char *path,
                                GPtrArray *unrepresentable)
{
  dt_lens_calib_t *d = self->data;
  if(!path) return FALSE;
  if(!d->have_warp && !d->warp.have_vig
     && (!d->work.warps || !d->work.warps->len))
    return FALSE;

  dt_lens_profile_t prof;
  dt_lens_profile_init(&prof);

  gchar *ui_maker = dt_conf_get_string("plugins/lens_calib/maker");
  gchar *ui_model = dt_conf_get_string("plugins/lens_calib/model");
  gchar *ui_mount = dt_conf_get_string("plugins/lens_calib/mount");
  if(ui_maker && *ui_maker) g_strlcpy(prof.maker, ui_maker, sizeof(prof.maker));
  if(ui_model && *ui_model) g_strlcpy(prof.model, ui_model, sizeof(prof.model));
  if(ui_mount && *ui_mount) g_strlcpy(prof.mount, ui_mount, sizeof(prof.mount));
  g_free(ui_maker);
  g_free(ui_model);
  g_free(ui_mount);

  const float ui_crop = dt_conf_get_float("plugins/lens_calib/crop_factor");
  prof.crop_factor = ui_crop > 0.0f ? ui_crop : 1.0f;

  float focal = 0.0f, aperture = 0.0f;
  if(dt_is_valid_imgid(d->imgid))
  {
    const dt_image_t *img = dt_image_cache_get(d->imgid, 'r');
    if(img)
    {
      focal = img->exif_focal_length;
      aperture = img->exif_aperture;
      dt_image_cache_read_release(img);
    }
  }
  const float ui_focal = dt_conf_get_float("plugins/lens_calib/focal");
  const float ui_aperture = dt_conf_get_float("plugins/lens_calib/aperture");
  if(ui_focal > 0.0f) focal = ui_focal;
  if(ui_aperture > 0.0f) aperture = ui_aperture;

  gboolean added = FALSE;
  if(d->work.warps && d->work.warps->len)
  {
    for(guint i = 0; i < d->work.warps->len; i++)
    {
      const dt_lens_warp_t *src =
        &g_array_index(d->work.warps, dt_lens_warp_t, i);
      dt_lens_warp_t w;
      if(!dt_lens_warp_copy(&w, src)) continue;
      if(w.focal <= 0.0f) w.focal = focal;
      if(w.aperture <= 0.0f && aperture > 0.0f) w.aperture = aperture;
      if(dt_lens_profile_add(&prof, &w)) added = TRUE;
      dt_lens_warp_cleanup(&w);
    }
  }
  else
  {
    dt_lens_warp_t w;
    if(dt_lens_warp_copy(&w, &d->warp))
    {
      w.focal = focal;
      if(aperture > 0.0f) w.aperture = aperture;
      added = dt_lens_profile_add(&prof, &w);
      dt_lens_warp_cleanup(&w);
    }
  }

  if(!added)
  {
    dt_lens_profile_cleanup(&prof);
    return FALSE;
  }

  GError *err = NULL;
  const gboolean ok =
    dt_lens_profile_export_lensfun(&prof, path, unrepresentable, &err);

  if(!ok && err)
    dt_print(DT_DEBUG_ALWAYS, "[lensfit] could not export lensfun xml: %s",
             err->message);
  if(err) g_error_free(err);

  dt_lens_profile_cleanup(&prof);
  return ok;
}

/* Switch calibration image when a filmstrip thumbnail is activated, the
   same way darkroom does. */
static void _filmstrip_activate_callback(gpointer instance,
                                         const dt_imgid_t imgid,
                                         dt_view_t *self)
{
  dt_lens_calib_t *d = self->data;
  if(!dt_is_valid_imgid(imgid) || imgid == d->imgid) return;

  // never lose hand placed work by navigating away from it
  _save_points(d);

  d->imgid = imgid;
  _drop_surface(d);
  _drop_grid(d);
  _drop_solution(d);
  _load_points(d);
  // a different image means a different framing; start fitted again
  d->zoom = 1.0f;
  d->pan_x = d->pan_y = 0.5f;

  if(dt_conf_get_bool("filmstrip/ui/auto_scroll"))
    dt_thumbtable_set_offset_image(dt_ui_thumbtable(darktable.gui->ui),
                                   imgid, TRUE);
  dt_control_queue_redraw();
}

void enter(dt_view_t *self)
{
  dt_lens_calib_t *d = self->data;

  // carry over whatever the rest of the ui was pointing at
  const dt_imgid_t acted_on = _pick_image();
  const gboolean changed =
    dt_is_valid_imgid(acted_on) && acted_on != d->imgid;
  if(changed)
  {
    _save_points(d);
    d->imgid = acted_on;
    _drop_grid(d);
    _drop_solution(d);
  }

  d->cells_x = MAX(1, dt_conf_get_int("plugins/lens_calib/cells_x"));
  d->cells_y = MAX(1, dt_conf_get_int("plugins/lens_calib/cells_y"));

  /* Reload every time we enter, not just on change: the sidecar may have
     been written by a previous session, and _load_points also restores the
     chart geometry the points were placed against. */
  _load_points(d);

  if(dt_is_valid_imgid(d->imgid))
    dt_thumbtable_set_offset_image(dt_ui_thumbtable(darktable.gui->ui),
                                   d->imgid, TRUE);

  DT_CONTROL_SIGNAL_HANDLE(DT_SIGNAL_VIEWMANAGER_THUMBTABLE_ACTIVATE,
                           _filmstrip_activate_callback);

  _drop_surface(d);
  dt_control_queue_redraw_center();
}

void leave(dt_view_t *self)
{
  dt_lens_calib_t *d = self->data;
  _save_points(d);
  // the surface can be large, don't hold onto it while we're not visible
  _drop_surface(d);
}

void configure(dt_view_t *self, int width, int height)
{
  dt_lens_calib_t *d = self->data;
  // force a rebuild at the new size on the next expose
  _drop_surface(d);
}

/* ------------------------------------------------- coordinate mapping */

/* Points are stored normalized over the frame, 0..1. The warp models work
   in centred half-diagonal units. Only the frame's aspect ratio connects
   the two, which is why this is safe to compute from the display surface:
   the ratio is the same at every resolution. */
static void _norm_to_warp(const dt_lens_calib_t *d,
                          const float nx, const float ny,
                          float *u, float *v)
{
  const double ax = d->surf_width > 0 ? d->surf_width : 3.0;
  const double ay = d->surf_height > 0 ? d->surf_height : 2.0;
  const double hd = 0.5 * hypot(ax, ay);

  *u = (float)(((double)nx - 0.5) * ax / hd);
  *v = (float)(((double)ny - 0.5) * ay / hd);
}

static void _warp_to_norm(const dt_lens_calib_t *d,
                          const float u, const float v,
                          float *nx, float *ny)
{
  const double ax = d->surf_width > 0 ? d->surf_width : 3.0;
  const double ay = d->surf_height > 0 ? d->surf_height : 2.0;
  const double hd = 0.5 * hypot(ax, ay);

  *nx = (float)((double)u * hd / ax + 0.5);
  *ny = (float)((double)v * hd / ay + 0.5);
}

/* Geometry preview: the frame remapped so the chart lines come out straight.
   Distinct from the falloff preview below, which changes brightness and
   leaves every pixel where it is -- two different corrections that happen to
   be measured by the same tool, and worth being able to judge separately. */
static gboolean _flat_active(const dt_lens_calib_t *d)
{
  return d->flat_view && d->have_warp;
}

static gboolean _falloff_active(const dt_lens_calib_t *d)
{
  return d->falloff_view && d->warp.have_vig;
}

// either correction means the displayed pixels are not the stored ones
static gboolean _corrected_active(const dt_lens_calib_t *d)
{
  return _flat_active(d) || _falloff_active(d);
}

/* The desqueeze in effect on the displayed frame.
 *
 * A 1.6x anamorphic records a scene 1.6 times wider than the frame shape
 * suggests, so the corrected frame is 1.6 times wider than the recorded one.
 * That is a change of canvas, not of content: drawing the corrected image
 * into the source's rectangle would put the squeeze straight back.
 */
static double _flat_squeeze(const dt_lens_calib_t *d)
{
  if(!_flat_active(d)) return 1.0;
  return (d->warp.squeeze > 0.01f) ? (double)d->warp.squeeze : 1.0;
}

/* Warp coordinates to a fraction of the *corrected* frame.
 *
 * Distinct from _warp_to_norm, which lands on the recorded frame. The two
 * agreed as long as nothing changed the canvas, and a desqueeze does. */
static void _warp_to_flat_norm(const dt_lens_calib_t *d,
                               const float u, const float v,
                               float *nx, float *ny)
{
  const double ax = d->surf_width > 0 ? d->surf_width : 3.0;
  const double ay = d->surf_height > 0 ? d->surf_height : 2.0;
  const double hd = 0.5 * hypot(ax, ay);

  *nx = (float)((double)u * hd / (ax * _flat_squeeze(d)) + 0.5);
  *ny = (float)((double)v * hd / ay + 0.5);
}

static void _flat_norm_to_warp(const dt_lens_calib_t *d,
                               const float nx, const float ny,
                               float *u, float *v)
{
  const double ax = d->surf_width > 0 ? d->surf_width : 3.0;
  const double ay = d->surf_height > 0 ? d->surf_height : 2.0;
  const double hd = 0.5 * hypot(ax, ay);

  *u = (float)(((double)nx - 0.5) * ax * _flat_squeeze(d) / hd);
  *v = (float)(((double)ny - 0.5) * ay / hd);
}

/* Normalized image position to screen. In flat view the overlays have to
   travel with the picture, or the crosses would sit where the points used
   to be and the whole display would contradict itself. */
static void _norm_to_view(const dt_lens_calib_t *d,
                          const float nx, const float ny,
                          double *sx, double *sy)
{
  float px = nx, py = ny;

  if(_flat_active(d))
  {
    float u, v, ou, ov;
    _norm_to_warp(d, nx, ny, &u, &v);
    dt_lens_warp_apply(&d->warp, u, v, &ou, &ov);
    _warp_to_flat_norm(d, ou, ov, &px, &py);
  }

  *sx = d->img_x + px * d->img_w;
  *sy = d->img_y + py * d->img_h;
}

/* Screen back to a normalized position on the *original* frame, which is
   where points live regardless of which view is on show. Returns FALSE
   off the image. */
static gboolean _view_to_norm(const dt_lens_calib_t *d,
                              const double sx,
                              const double sy,
                              float *nx,
                              float *ny)
{
  if(d->img_w <= 0 || d->img_h <= 0) return FALSE;

  double rx = (sx - d->img_x) / d->img_w;
  double ry = (sy - d->img_y) / d->img_h;
  if(rx < 0.0 || rx > 1.0 || ry < 0.0 || ry > 1.0) return FALSE;

  if(_flat_active(d))
  {
    float u, v, su, sv;
    _flat_norm_to_warp(d, (float)rx, (float)ry, &u, &v);
    dt_lens_warp_invert(&d->warp, u, v, &su, &sv);

    float bx, by;
    _warp_to_norm(d, su, sv, &bx, &by);
    rx = bx;
    ry = by;
  }

  *nx = (float)CLAMP(rx, 0.0, 1.0);
  *ny = (float)CLAMP(ry, 0.0, 1.0);
  return TRUE;
}

/* Corrected copy of the displayed image.
 *
 * Capped well below the display surface: this is a check on the fit, not
 * a surface to place markers on, and an inverse warp per pixel over an
 * 8000 pixel zoomed surface would stall the interface for seconds.
 */
#define DT_LENS_CALIB_FLAT_MAX 2600

static void _build_flat(dt_lens_calib_t *d)
{
  _drop_flat(d);

  if(!d->surface || !_corrected_active(d)) return;
  if(d->surf_width < 2 || d->surf_height < 2) return;

  const gboolean geo = _flat_active(d);
  const gboolean vig = _falloff_active(d);

  cairo_surface_flush(d->surface);
  const uint8_t *const src = cairo_image_surface_get_data(d->surface);
  const int sstride = cairo_image_surface_get_stride(d->surface);
  if(!src) return;

  const int sw = d->surf_width, sh = d->surf_height;

  /* The corrected frame, which a desqueeze makes wider than the recorded
     one. Rendering into an sw x sh canvas would crop the sides off and
     leave the content squeezed exactly as it was found. */
  const double sq = geo && (d->warp.squeeze > 0.01f)
    ? (double)d->warp.squeeze : 1.0;
  const double fw = sw * sq, fh = sh;

  const double longest = MAX(fw, fh);
  const double k = (longest > DT_LENS_CALIB_FLAT_MAX)
    ? DT_LENS_CALIB_FLAT_MAX / longest : 1.0;

  const int dw = MAX(2, (int)(fw * k));
  const int dh = MAX(2, (int)(fh * k));

  cairo_surface_t *out =
    cairo_image_surface_create(CAIRO_FORMAT_RGB24, dw, dh);
  if(cairo_surface_status(out) != CAIRO_STATUS_SUCCESS)
  {
    cairo_surface_destroy(out);
    return;
  }

  uint8_t *const dst = cairo_image_surface_get_data(out);
  const int dstride = cairo_image_surface_get_stride(out);
  const double hd = 0.5 * hypot((double)sw, (double)sh);

  DT_OMP_FOR()
  for(int y = 0; y < dh; y++)
  {
    uint8_t *drow = dst + (size_t)y * dstride;
    for(int x = 0; x < dw; x++)
    {
      const double nx = (x + 0.5) / dw;
      const double ny = (y + 0.5) / dh;

      /* hd stays the *source* half-diagonal: the warp is defined on the
         recorded frame, and only the destination canvas has grown. */
      const float u = (float)((nx - 0.5) * fw / hd);
      const float v = (float)((ny - 0.5) * fh / hd);

      /* Only the geometry preview moves pixels. With just the falloff
         correction on, every pixel stays put and only its value changes. */
      float su = u, sv = v;
      if(geo) dt_lens_warp_invert(&d->warp, u, v, &su, &sv);

      const double fx = su * hd + 0.5 * sw - 0.5;
      const double fy = sv * hd + 0.5 * sh - 0.5;

      const int x0 = (int)floor(fx), y0 = (int)floor(fy);

      if(x0 < 0 || y0 < 0 || x0 >= sw - 1 || y0 >= sh - 1)
      {
        // outside the frame: leave it black rather than smearing the edge
        drow[4 * x + 0] = drow[4 * x + 1] = drow[4 * x + 2] = 0;
        drow[4 * x + 3] = 0xff;
        continue;
      }

      const double tx = fx - x0, ty = fy - y0;
      const uint8_t *const a = src + (size_t)y0 * sstride + 4 * x0;
      const uint8_t *const b = src + (size_t)(y0 + 1) * sstride + 4 * x0;

      /* The gain belongs at the source position: vignetting is a property of
         where the light landed on the sensor, so applying it in corrected
         space would smear the pattern relative to what it corrects. */
      const double gain = vig ? dt_lens_vignette_gain(&d->warp, su, sv) : 1.0;

      for(int c = 0; c < 3; c++)
      {
        const double top = a[c] + tx * (a[4 + c] - a[c]);
        const double bot = b[c] + tx * (b[4 + c] - b[c]);
        double value = top + ty * (bot - top);

        /* This surface is display encoded, and a brightness gain is a
           statement about light. Multiplying the encoded value directly would
           under-correct the shadows and over-correct the highlights, so the
           gain is applied in a rough linear space and encoded back. The 2.2
           power is an approximation to the display curve -- good enough for
           judging a correction by eye, which is all this surface is for. */
        if(vig && gain != 1.0)
        {
          const double lin = pow(CLAMP(value, 0.0, 255.0) / 255.0, 2.2);
          value = 255.0 * pow(CLAMP(lin * gain, 0.0, 1.0), 1.0 / 2.2);
        }

        drow[4 * x + c] = (uint8_t)CLAMP(value, 0.0, 255.0);
      }
      drow[4 * x + 3] = 0xff;
    }
  }

  cairo_surface_mark_dirty(out);
  d->flat_surface = out;
  d->flat_valid = TRUE;
}

/* The mesh: every node joined to its lattice neighbours.
 *
 * This is the overlay whose crossings coincide with the points by
 * construction, because it *is* the points -- it draws the lattice as
 * measured rather than as it ought to be. That makes it the diagnostic the
 * reference grid cannot be: a mis-indexed node shows up immediately as a
 * zigzag or a line doubling back, where in a field of loose crosses it is
 * invisible.
 *
 * The reference grid answers a different question and must not be confused
 * with this one. It is drawn straight, so on a distorted frame it *cannot*
 * pass through the points -- and if it did, there would be no distortion to
 * calibrate. The place the two are expected to coincide is flat view, after
 * a good fit: that is the success criterion, and it lives in corrected
 * space where the comparison is meaningful.
 */
static void _draw_mesh(cairo_t *cr, const dt_lens_calib_t *d)
{
  if(!d->show_mesh || !d->manual || !d->manual->len) return;
  if(d->img_w <= 0 || d->img_h <= 0) return;

  const int nx = d->cells_x + 1;
  const int ny = d->cells_y + 1;
  if(nx < 2 || ny < 2) return;

  // lattice position to point index, so adjacency is a lookup not a search
  int *lookup = malloc(sizeof(int) * (size_t)nx * ny);
  if(!lookup) return;
  for(size_t i = 0; i < (size_t)nx * ny; i++) lookup[i] = -1;

  for(guint i = 0; i < d->manual->len; i++)
  {
    const dt_lens_calib_manual_point_t *p =
      &g_array_index(d->manual, dt_lens_calib_manual_point_t, i);
    if(p->col < 0 || p->row < 0 || p->col >= nx || p->row >= ny) continue;
    lookup[(size_t)p->row * nx + p->col] = (int)i;
  }

  cairo_save(cr);
  cairo_rectangle(cr, d->img_x, d->img_y, d->img_w, d->img_h);
  cairo_clip(cr);
  cairo_set_line_width(cr, DT_PIXEL_APPLY_DPI(1.0));

  for(int row = 0; row < ny; row++)
    for(int col = 0; col < nx; col++)
    {
      const int a = lookup[(size_t)row * nx + col];
      if(a < 0) continue;

      const dt_lens_calib_manual_point_t *pa =
        &g_array_index(d->manual, dt_lens_calib_manual_point_t, a);

      // right and down only, so each edge is drawn once
      const int nb[2] = { (col + 1 < nx) ? lookup[(size_t)row * nx + col + 1] : -1,
                          (row + 1 < ny) ? lookup[(size_t)(row + 1) * nx + col] : -1 };

      for(int q = 0; q < 2; q++)
      {
        if(nb[q] < 0) continue;

        const dt_lens_calib_manual_point_t *pb =
          &g_array_index(d->manual, dt_lens_calib_manual_point_t, nb[q]);

        /* An edge between two measured nodes is evidence; one touching an
           interpolated node is partly our own invention, so it is drawn
           fainter. */
        if(pa->measured && pb->measured)
          cairo_set_source_rgba(cr, 0.25, 0.85, 1.0, 0.75);
        else
          cairo_set_source_rgba(cr, 0.4, 0.6, 0.8, 0.30);

        double ax, ay, bx, by;
        _norm_to_view(d, pa->x, pa->y, &ax, &ay);
        _norm_to_view(d, pb->x, pb->y, &bx, &by);

        cairo_move_to(cr, ax, ay);
        cairo_line_to(cr, bx, by);
        cairo_stroke(cr);
      }
    }

  cairo_restore(cr);
  free(lookup);
}

/* Draw what the detector actually found: the traced curves and the
   intersections derived from them. This is the honest readout -- if the
   tracer latched onto noise instead of chart lines, it shows here. */
static void _draw_detected(cairo_t *cr,
                           const dt_lens_calib_t *d)
{
  if(!d->have_grid) return;
  if(d->detect_img_w <= 0 || d->img_w <= 0) return;

  // detection ran on its own surface; normalize out of it before drawing
  const double inv_w = 1.0 / (double)d->detect_img_w;
  const double inv_h = 1.0 / (double)d->detect_img_h;

  cairo_save(cr);

  if(d->show_curves)
  {
    cairo_set_line_width(cr, DT_PIXEL_APPLY_DPI(1.0));
    cairo_set_source_rgba(cr, 0.2, 0.9, 1.0, 0.55);

    for(int pass = 0; pass < 2; pass++)
    {
      for(GList *it = pass ? d->grid.curves_v : d->grid.curves_h;
          it; it = g_list_next(it))
      {
        const dt_lens_grid_curve_t *cv = it->data;
        for(int i = 0; i < cv->count; i++)
        {
          // horizontal curves are sampled along x, vertical along y
          const double px = cv->horizontal ? cv->samples[i].pos
                                           : cv->samples[i].coord;
          const double py = cv->horizontal ? cv->samples[i].coord
                                           : cv->samples[i].pos;
          double vx, vy;
          _norm_to_view(d, (float)(px * inv_w), (float)(py * inv_h),
                        &vx, &vy);
          if(i == 0) cairo_move_to(cr, vx, vy);
          else       cairo_line_to(cr, vx, vy);
        }
        cairo_stroke(cr);
      }
    }
  }

  if(d->show_points)
  {
    const double r = DT_PIXEL_APPLY_DPI(2.5);
    cairo_set_source_rgba(cr, 1.0, 0.9, 0.2, 0.9);
    for(int i = 0; i < d->grid.point_count; i++)
    {
      const dt_lens_grid_point_t *p = &d->grid.points[i];
      double vx, vy;
      _norm_to_view(d, (float)(p->x * inv_w), (float)(p->y * inv_h),
                    &vx, &vy);
      cairo_arc(cr, vx, vy, r, 0, 2.0 * M_PI);
      cairo_fill(cr);
    }
  }

  cairo_restore(cr);
}


// index of the manual point within `tol` screen pixels, or -1
static int _manual_at(const dt_lens_calib_t *d,
                      const double sx,
                      const double sy,
                      const double tol)
{
  int best = -1;
  double best_d2 = tol * tol;

  for(guint i = 0; i < d->manual->len; i++)
  {
    const dt_lens_calib_manual_point_t *p =
      &g_array_index(d->manual, dt_lens_calib_manual_point_t, i);
    double px, py;
    _norm_to_view(d, p->x, p->y, &px, &py);
    const double dx = px - sx, dy = py - sy;
    const double d2 = dx * dx + dy * dy;
    if(d2 < best_d2)
    {
      best_d2 = d2;
      best = (int)i;
    }
  }
  return best;
}


/* Scroll zooms about the pointer, so the feature under the cursor stays
   under the cursor -- placing a marker means zooming into it, and having
   the target drift away on each notch would be maddening. */
void scrolled(dt_view_t *self, double x, double y, int up, int state)
{
  dt_lens_calib_t *d = self->data;

  /* Pan is a property of the displayed rectangle, so this deliberately
     uses the raw position within it rather than _view_to_norm: in flat
     view the latter un-warps to a source coordinate, which is the right
     answer for a point and the wrong one for a scroll position. */
  float nx = 0.0f, ny = 0.0f;
  gboolean over_image = FALSE;

  if(d->img_w > 0 && d->img_h > 0)
  {
    const double rx = (x - d->img_x) / d->img_w;
    const double ry = (y - d->img_y) / d->img_h;
    over_image = rx >= 0.0 && rx <= 1.0 && ry >= 0.0 && ry <= 1.0;
    nx = (float)rx;
    ny = (float)ry;
  }

  const float old = d->zoom;
  d->zoom = up ? d->zoom * DT_LENS_CALIB_ZOOM_STEP
               : d->zoom / DT_LENS_CALIB_ZOOM_STEP;
  d->zoom = CLAMPF(d->zoom, DT_LENS_CALIB_ZOOM_MIN, DT_LENS_CALIB_ZOOM_MAX);
  if(d->zoom == old) return;

  if(over_image && d->zoom > 1.0f)
  {
    /* Shift pan so the image point under the cursor stays at the same
       screen position. The cursor sits at some fraction of the viewport;
       after the zoom that fraction must still map to (nx, ny). */
    const int vw = MAX(1, d->vp_w);
    const int vh = MAX(1, d->vp_h);
    const double fx = (x - 0.5 * vw) / (double)vw;
    const double fy = (y - 0.5 * vh) / (double)vh;
    // one viewport spans 1/zoom of the image once zoomed
    d->pan_x = nx - (float)(fx / d->zoom);
    d->pan_y = ny - (float)(fy / d->zoom);
    d->pan_x = CLAMPF(d->pan_x, 0.0f, 1.0f);
    d->pan_y = CLAMPF(d->pan_y, 0.0f, 1.0f);
  }

  // the surface has to be regenerated at the new scale
  _drop_surface(d);
  dt_control_queue_redraw_center();
}

int button_released(dt_view_t *self,
                    double x,
                    double y,
                    int which,
                    uint32_t state)
{
  dt_lens_calib_t *d = self->data;

  if(which == 2 && d->panning)
  {
    d->panning = FALSE;
    return 1;
  }

  if(which == 1 && (d->drag_idx >= 0 || d->corner_drag >= 0))
  {
    /* Save on release rather than on every motion event: a drag across
       the frame is hundreds of events and each save rewrites the whole
       sidecar. */
    const gboolean was_corner = d->corner_drag >= 0;
    d->drag_idx = -1;
    d->corner_drag = -1;
    _save_points(d);

    /* The fourth corner answers the question corner mode was asked to answer,
       so the mode ends here rather than waiting to be switched off. From now
       on the corners come from the pose. */
    if(was_corner && d->corner_mode && d->corner_count == 4)
    {
      d->undo_hold++;
      _align_grid(self);
      d->undo_hold--;
      _set_corner_mode(self, FALSE);
    }
    return 1;
  }

  return 0;
}

void mouse_moved(dt_view_t *self,
                 double x,
                 double y,
                 double pressure,
                 int which)
{
  dt_lens_calib_t *d = self->data;

  if(d->panning && d->img_w > 0 && d->img_h > 0)
  {
    // drag moves the image with the cursor, so pan travels the other way
    d->pan_x = d->pan_ref_px - (float)((x - d->pan_ref_x) / d->img_w);
    d->pan_y = d->pan_ref_py - (float)((y - d->pan_ref_y) / d->img_h);
    d->pan_x = CLAMPF(d->pan_x, 0.0f, 1.0f);
    d->pan_y = CLAMPF(d->pan_y, 0.0f, 1.0f);
    dt_control_queue_redraw_center();
    return;
  }

  if(d->corner_drag >= 0 && d->corner_drag < d->corner_count)
  {
    float nx, ny;
    if(_view_to_norm(d, x, y, &nx, &ny))
    {
      d->corner_x[d->corner_drag] = nx;
      d->corner_y[d->corner_drag] = ny;
      dt_control_queue_redraw_center();
    }
    return;
  }

  if(d->drag_idx >= 0 && d->drag_idx < (int)d->manual->len)
  {
    float nx, ny;
    if(_view_to_norm(d, x, y, &nx, &ny))
    {
      dt_lens_calib_manual_point_t *p =
        &g_array_index(d->manual, dt_lens_calib_manual_point_t, d->drag_idx);
      p->x = nx;
      p->y = ny;

      /* Moving a node is a measurement of it, so it stops being a guess and
         starts constraining the fit. This is the one place an interpolated
         node is promoted, and it has to be here: a node the user has
         deliberately positioned is evidence in a way one we invented is
         not. */
      p->measured = TRUE;

      dt_control_queue_redraw_center();
    }
    return;
  }

  const int was = d->hover_manual;
  d->hover_manual = d->manual_edit
    ? _manual_at(d, x, y, DT_PIXEL_APPLY_DPI(9)) : -1;
  if(d->hover_manual != was)
    dt_control_queue_redraw_center();
}

// index of the pinned corner within `tol` screen pixels, or -1
static int _corner_at(const dt_lens_calib_t *d,
                      const double sx,
                      const double sy,
                      const double tol)
{
  int best = -1;
  double best_d2 = tol * tol;

  for(int i = 0; i < d->corner_count; i++)
  {
    double px, py;
    _norm_to_view(d, d->corner_x[i], d->corner_y[i], &px, &py);
    const double dx = px - sx, dy = py - sy;
    const double d2 = dx * dx + dy * dy;
    if(d2 < best_d2)
    {
      best_d2 = d2;
      best = i;
    }
  }
  return best;
}

int button_pressed(dt_view_t *self,
                   double x,
                   double y,
                   double pressure,
                   int which,
                   int type,
                   uint32_t state)
{
  dt_lens_calib_t *d = self->data;

  // middle drag pans while zoomed in
  if(which == 2)
  {
    d->panning = TRUE;
    d->pan_ref_x = x;
    d->pan_ref_y = y;
    d->pan_ref_px = d->pan_x;
    d->pan_ref_py = d->pan_y;
    return 1;
  }

  /* Corner pinning comes before everything else, so that placing the four
     corners never has to compete with the hundreds of nodes for clicks. */
  if(d->corner_mode)
  {
    const int hit = _corner_at(d, x, y, DT_PIXEL_APPLY_DPI(12));

    if(which == 3)
    {
      if(hit >= 0)
      {
        _undo_push(d, _("remove corner"));
        for(int i = hit; i < d->corner_count - 1; i++)
        {
          d->corner_x[i] = d->corner_x[i + 1];
          d->corner_y[i] = d->corner_y[i + 1];
        }
        d->corner_count--;
        d->corner_drag = -1;
        _save_points(d);
        dt_control_queue_redraw_center();
      }
      return 1;
    }

    if(which != 1) return 0;

    if(hit >= 0)
    {
      /* Record on press, not on release: this is the last moment the old
         position still exists. */
      _undo_push(d, _("move corner"));
      d->corner_drag = hit;
      return 1;
    }

    float cnx, cny;
    if(!_view_to_norm(d, x, y, &cnx, &cny)) return 0;

    if(d->corner_count >= 4)
    {
      dt_control_log(_("all four corners are placed;"
                       " drag one, or right click to remove it"));
      return 1;
    }

    _undo_push(d, _("place corner"));
    d->corner_x[d->corner_count] = cnx;
    d->corner_y[d->corner_count] = cny;
    d->corner_drag = d->corner_count;
    d->corner_count++;

    dt_control_queue_redraw_center();
    return 1;
  }

  /* Everything below edits points, which only happens once the user has
     said they want to. Otherwise a stray click on a detection result
     would quietly turn it into hand placed data. */
  if(!d->manual_edit) return 0;

  // right click removes the point under the pointer
  if(which == 3)
  {
    const int idx = _manual_at(d, x, y, DT_PIXEL_APPLY_DPI(9));
    if(idx >= 0)
    {
      _undo_push(d, _("remove point"));
      g_array_remove_index(d->manual, idx);
      d->hover_manual = -1;
      d->drag_idx = -1;
      _save_points(d);
      dt_control_queue_redraw_center();
      return 1;
    }
    return 0;
  }

  if(which != 1) return 0;

  /* Grabbing an existing point takes priority over adding one. The point
     is *not* moved to the click position: it is picked up where it is and
     follows the pointer, so a slightly off grab nudges nothing. */
  const int idx = _manual_at(d, x, y, DT_PIXEL_APPLY_DPI(9));
  if(idx >= 0)
  {
    _undo_push(d, _("move point"));
    d->drag_idx = idx;
    d->hover_manual = idx;
    return 1;
  }

  float nx, ny;
  if(!_view_to_norm(d, x, y, &nx, &ny)) return 0;

  /* A freshly placed point has no lattice index unless the pose can supply
     one. With a pose fitted, inverting it says which node the click is
     nearest, which is exactly the index the point needs -- and means adding
     a node to a mesh does not silently create an unindexed point that sits
     the fit out. */
  int col = -1, row = -1;

  _undo_push(d, _("place point"));

  if(d->have_pose)
  {
    /* The pose maps index to position; search the lattice rather than
       inverting the homography, since the lattice is small and this is
       exact rather than nearly so. */
    double best = 1e30;
    for(int r = 0; r <= d->cells_y; r++)
      for(int c = 0; c <= d->cells_x; c++)
      {
        float ex, ey;
        _apply_homography(d->pose, c, r, &ex, &ey);
        const double dx = ex - nx, dy = ey - ny;
        const double dist = dx * dx + dy * dy;
        if(dist < best)
        {
          best = dist;
          col = c;
          row = r;
        }
      }

    // an existing node at that index is moved rather than duplicated
    const int existing = _node_at_index(d, col, row);
    if(existing >= 0)
    {
      dt_lens_calib_manual_point_t *p =
        &g_array_index(d->manual, dt_lens_calib_manual_point_t, existing);
      p->x = nx;
      p->y = ny;
      p->measured = TRUE;
      d->drag_idx = existing;
      d->hover_manual = existing;
      dt_control_queue_redraw_center();
      return 1;
    }
  }

  const dt_lens_calib_manual_point_t np = { nx, ny, col, row, TRUE };
  g_array_append_val(d->manual, np);

  // a new point is dragged straight away, so placing and refining are one
  d->drag_idx = d->manual->len - 1;
  d->hover_manual = d->drag_idx;

  dt_control_queue_redraw_center();
  return 1;
}

/* Residual vectors after the fit.
 *
 * Drawn exaggerated, because a good fit leaves residuals of well under a
 * pixel and an honest 1:1 arrow would be invisible. The exaggeration
 * factor is stated in the panel readout, and the colour saturates with
 * magnitude so a single bad point stands out from a field of small ones --
 * which is the failure this overlay exists to catch: one misplaced marker
 * dragging the whole fit.
 */
/* The longest residual arrow to draw, in screen pixels. The exaggeration
   is derived from this rather than fixed: a good fit leaves well under a
   pixel and needs magnifying to be visible at all, while a broken one
   leaves a hundred, and magnifying *that* forty times drew lines several
   screens long, crossing the whole frame and hiding the picture. Scaling to
   a fixed maximum keeps the overlay readable at both extremes. */
#define DT_LENS_CALIB_RESIDUAL_MAX_PX 36.0

static void _draw_residuals(cairo_t *cr,
                            const dt_lens_calib_t *d)
{
  if(!d->show_residuals || !d->have_warp) return;
  if(!d->solve_pts || !d->res_dx || !d->res_dy) return;
  if(d->solve_img_w <= 0 || d->img_w <= 0) return;

  const double sx = d->img_w / (double)d->solve_img_w;
  const double sy = d->img_h / (double)d->solve_img_h;
  const double inv_w = 1.0 / (double)d->solve_img_w;
  const double inv_h = 1.0 / (double)d->solve_img_h;

  cairo_save(cr);
  cairo_set_line_width(cr, DT_PIXEL_APPLY_DPI(1.2));

  // a residual at this size is worth looking at rather than ignoring
  const float bad = MAX(0.5f, d->result.rms_px * 3.0f);

  /* Scale the exaggeration so the worst arrow is a readable length, rather
     than magnifying by a constant regardless of what the numbers are. */
  float worst = 0.0f;
  if(d->res_mag)
    for(int i = 0; i < d->solve_count; i++)
      worst = MAX(worst, d->res_mag[i]);

  const double gain = (worst > 1e-6f)
    ? MIN(60.0, DT_LENS_CALIB_RESIDUAL_MAX_PX / (worst * MAX(sx, sy)))
    : 1.0;

  for(int i = 0; i < d->solve_count; i++)
  {
    double px, py;
    _norm_to_view(d, (float)(d->solve_pts[i].x * inv_w),
                  (float)(d->solve_pts[i].y * inv_h), &px, &py);
    const double vx = d->res_dx[i] * sx * gain;
    const double vy = d->res_dy[i] * sy * gain;

    const double t = CLAMP(d->res_mag ? d->res_mag[i] / bad : 0.0, 0.0, 1.0);
    cairo_set_source_rgba(cr, 0.3 + 0.7 * t, 1.0 - 0.8 * t, 0.3, 0.9);

    cairo_move_to(cr, px, py);
    cairo_line_to(cr, px + vx, py + vy);
    cairo_stroke(cr);
  }

  cairo_restore(cr);
}

/* Hand-placed points, drawn as crosses so they read differently from the
   detector's filled dots. */
static void _draw_manual(cairo_t *cr,
                         const dt_lens_calib_t *d)
{
  if(!d->manual) return;

  const double r = DT_PIXEL_APPLY_DPI(5.0);
  cairo_set_line_width(cr, DT_PIXEL_APPLY_DPI(1.4));

  /* These are the points now: once the detected set has been adopted there
     is no other, so the points switch has to govern them. It used to gate
     only the detector's own overlay, which after adoption has nothing left
     to draw -- so the switch appeared to do nothing at all. */
  if(!d->show_points) goto corners;

  for(guint i = 0; i < d->manual->len; i++)
  {
    const dt_lens_calib_manual_point_t *p =
      &g_array_index(d->manual, dt_lens_calib_manual_point_t, i);

    // a guess still needs a real point to anchor the mesh/points overlay to
    if(!p->measured) continue;

    double px, py;
    _norm_to_view(d, p->x, p->y, &px, &py);

    if((int)i == d->drag_idx)
      cairo_set_source_rgba(cr, 1.0, 0.9, 0.3, 1.0);
    else if((int)i == d->hover_manual)
      cairo_set_source_rgba(cr, 1.0, 1.0, 1.0, 1.0);
    else if(!p->measured)
      /* Interpolated: drawn faint and small, because it is a guess and
         must not read as evidence. Anything the fit is actually using
         looks solid. */
      cairo_set_source_rgba(cr, 0.45, 0.75, 1.0, 0.5);
    else if(p->col < 0 || p->row < 0)
      /* Measured, but on no lattice site -- so it joins no line and the
         solver never sees it. It has to look different from a point that
         counts, or it reads as evidence that is silently being ignored.
         Drag it onto its crossing, or right click it away. */
      cairo_set_source_rgba(cr, 1.0, 0.45, 0.2, 0.9);
    else if(d->manual_edit)
      cairo_set_source_rgba(cr, 0.3, 1.0, 0.4, 0.9);
    else
      // not editable at the moment, so don't advertise it as green
      cairo_set_source_rgba(cr, 0.55, 0.7, 0.6, 0.7);

    if(!p->measured)
    {
      const double q = DT_PIXEL_APPLY_DPI(1.6);
      cairo_rectangle(cr, px - q, py - q, 2.0 * q, 2.0 * q);
      cairo_fill(cr);
      continue;
    }

    cairo_move_to(cr, px - r, py);
    cairo_line_to(cr, px + r, py);
    cairo_move_to(cr, px, py - r);
    cairo_line_to(cr, px, py + r);
    cairo_stroke(cr);
  }

corners:
  // corner handles only matter while pinning them by hand
  if(!d->corner_mode) return;

  // large enough to grab, and distinct from a node
  const double cr_r = DT_PIXEL_APPLY_DPI(8.0);
  cairo_set_line_width(cr, DT_PIXEL_APPLY_DPI(2.0));

  for(int i = 0; i < d->corner_count; i++)
  {
    double px, py;
    _norm_to_view(d, d->corner_x[i], d->corner_y[i], &px, &py);

    if(i == d->corner_drag)
      cairo_set_source_rgba(cr, 1.0, 1.0, 1.0, 1.0);
    else
      cairo_set_source_rgba(cr, 1.0, 0.55, 0.1, 0.95);

    cairo_arc(cr, px, py, cr_r, 0, 2.0 * M_PI);
    cairo_stroke(cr);
    cairo_arc(cr, px, py, DT_PIXEL_APPLY_DPI(1.5), 0, 2.0 * M_PI);
    cairo_fill(cr);
  }
}

static void _draw_no_image_hint(cairo_t *cr,
                                const int32_t width,
                                const int32_t height)
{
  PangoFontDescription *desc = pango_font_description_from_string("sans");
  pango_font_description_set_absolute_size(desc, DT_PIXEL_APPLY_DPI(14) * PANGO_SCALE);
  PangoLayout *layout = pango_cairo_create_layout(cr);
  pango_layout_set_font_description(layout, desc);
  pango_layout_set_text(layout,
                        _("select a photograph of a grid chart in lighttable"), -1);

  int pw, ph;
  pango_layout_get_pixel_size(layout, &pw, &ph);
  dt_gui_gtk_set_source_rgb(cr, DT_GUI_COLOR_LIGHTTABLE_FONT);
  cairo_move_to(cr, 0.5 * (width - pw), 0.5 * (height - ph));
  pango_cairo_show_layout(cr, layout);

  g_object_unref(layout);
  pango_font_description_free(desc);
}

void expose(dt_view_t *self,
            cairo_t *cri,
            int32_t width,
            int32_t height,
            int32_t pointerx,
            int32_t pointery)
{
  dt_lens_calib_t *d = self->data;

  d->vp_w = width;
  d->vp_h = height;

  dt_gui_gtk_set_source_rgb(cri, DT_GUI_COLOR_DARKROOM_BG);
  cairo_paint(cri);

  if(!dt_is_valid_imgid(d->imgid))
  {
    _draw_no_image_hint(cri, width, height);
    return;
  }

  // leave a small margin so the chart edges stay visible against the
  // background rather than butting up against the panels
  const int avail_w = MAX(1, width - (int)DT_PIXEL_APPLY_DPI(20));
  const int avail_h = MAX(1, height - (int)DT_PIXEL_APPLY_DPI(20));

  /* Render the surface at the zoomed size rather than upscaling a fitted
     one, so zooming in actually reveals detail instead of magnifying
     display-resolution mush -- the whole point of zooming here is placing
     markers precisely. Capped so a big zoom on a big raw can't ask for an
     absurd allocation. */
  int want_w = (int)(avail_w * d->zoom);
  int want_h = (int)(avail_h * d->zoom);
  const int cap = DT_LENS_CALIB_MAX_SURFACE;
  if(want_w > cap || want_h > cap)
  {
    const double k = (double)cap / MAX(want_w, want_h);
    want_w = (int)(want_w * k);
    want_h = (int)(want_h * k);
  }
  want_w = MAX(1, want_w);
  want_h = MAX(1, want_h);

  if(!d->surface || d->req_width != want_w || d->req_height != want_h)
  {
    _drop_surface(d);
    cairo_surface_t *surf = NULL;
    const dt_view_surface_value_t res =
      dt_view_image_get_surface(d->imgid, want_w, want_h, &surf, TRUE);

    if(res == DT_VIEW_SURFACE_OK && surf)
    {
      d->surface = surf;
      d->surf_width = cairo_image_surface_get_width(surf);
      d->surf_height = cairo_image_surface_get_height(surf);
      d->req_width = want_w;
      d->req_height = want_h;
    }
    else
    {
      if(surf) cairo_surface_destroy(surf);
      // mipmap not ready yet: come back once it has been generated
      dt_control_queue_redraw_center();
      return;
    }
  }

  d->img_w = d->surf_width;
  d->img_h = d->surf_height;

  /* A desqueezed frame is wider, so it has to be drawn wider. Fitting it
     back into the width the recorded frame occupied is the same shape as
     widening it, and keeps it inside the viewport. */
  d->img_h /= _flat_squeeze(d);

  if(d->zoom <= 1.0f)
  {
    // fitted: centre it and keep pan neutral so zooming back out recentres
    d->img_x = 0.5 * (width - d->img_w);
    d->img_y = 0.5 * (height - d->img_h);
    d->pan_x = d->pan_y = 0.5f;
  }
  else
  {
    // hold the panned-to image point at the centre of the viewport
    d->img_x = 0.5 * width - d->pan_x * d->img_w;
    d->img_y = 0.5 * height - d->pan_y * d->img_h;

    // don't allow panning past the edges into empty canvas
    if(d->img_w > width)
      d->img_x = CLAMP(d->img_x, width - d->img_w, 0.0);
    else
      d->img_x = 0.5 * (width - d->img_w);
    if(d->img_h > height)
      d->img_y = CLAMP(d->img_y, height - d->img_h, 0.0);
    else
      d->img_y = 0.5 * (height - d->img_h);
  }

  if(_corrected_active(d))
  {
    if(!d->flat_valid) _build_flat(d);

    if(d->flat_surface)
    {
      /* The corrected copy is generally smaller than the display surface,
         so it is scaled up here rather than rendered at full size -- see
         the cap on DT_LENS_CALIB_FLAT_MAX. */
      const double fw = cairo_image_surface_get_width(d->flat_surface);
      const double fh = cairo_image_surface_get_height(d->flat_surface);

      cairo_save(cri);
      cairo_translate(cri, d->img_x, d->img_y);
      cairo_scale(cri, d->img_w / fw, d->img_h / fh);
      cairo_set_source_surface(cri, d->flat_surface, 0, 0);
      cairo_pattern_set_filter(cairo_get_source(cri), CAIRO_FILTER_BILINEAR);
      cairo_paint(cri);
      cairo_restore(cri);
    }
  }
  else
  {
    cairo_set_source_surface(cri, d->surface, d->img_x, d->img_y);
    cairo_paint(cri);
  }

  _draw_detected(cri, d);
  _draw_mesh(cri, d);
  _draw_residuals(cri, d);
  _draw_manual(cri, d);
}

void gui_init(dt_view_t *self)
{
  dt_action_register(DT_ACTION(self), N_("undo"),
                     _undo_callback, GDK_KEY_z, GDK_CONTROL_MASK);
  dt_action_register(DT_ACTION(self), N_("redo"),
                     _redo_callback, GDK_KEY_y, GDK_CONTROL_MASK);

  // let the settings lib reach us without linking the two modules
  darktable.view_manager->proxy.lens_calib.view = self;
  darktable.view_manager->proxy.lens_calib.chart_changed = _chart_changed;
  darktable.view_manager->proxy.lens_calib.detect = _detect_grid;
  darktable.view_manager->proxy.lens_calib.has_points = _has_points;
  darktable.view_manager->proxy.lens_calib.solve = _solve;
  darktable.view_manager->proxy.lens_calib.has_solution = _has_solution;
  darktable.view_manager->proxy.lens_calib.status_text = _status_text;
  darktable.view_manager->proxy.lens_calib.save_profile = _save_profile;
  darktable.view_manager->proxy.lens_calib.contents_summary = _contents_summary;
  darktable.view_manager->proxy.lens_calib.entry_count = _entry_count;
  darktable.view_manager->proxy.lens_calib.entry_label = _entry_label;
  darktable.view_manager->proxy.lens_calib.entry_selected = _entry_selected;
  darktable.view_manager->proxy.lens_calib.entry_select = _entry_select;
  darktable.view_manager->proxy.lens_calib.entry_add = _entry_add;
  darktable.view_manager->proxy.lens_calib.entry_remove = _entry_remove;
  darktable.view_manager->proxy.lens_calib.axis_values = _axis_values;
  darktable.view_manager->proxy.lens_calib.entry_find = _entry_find;
  darktable.view_manager->proxy.lens_calib.value_count = _value_count;
  darktable.view_manager->proxy.lens_calib.value_name = _value_name;
  darktable.view_manager->proxy.lens_calib.value_get = _value_get;
  darktable.view_manager->proxy.lens_calib.value_set = _value_set;
  darktable.view_manager->proxy.lens_calib.open_profile = _open_profile;
  darktable.view_manager->proxy.lens_calib.new_profile = _new_profile;
  darktable.view_manager->proxy.lens_calib.measurement_range = _measurement_range;
  darktable.view_manager->proxy.lens_calib.import_lensfun = _import_lensfun;
  darktable.view_manager->proxy.lens_calib.export_stmap = _export_stmap;
  darktable.view_manager->proxy.lens_calib.export_lensfun = _export_lensfun;
  darktable.view_manager->proxy.lens_calib.set_manual_edit = _set_manual_edit;
  darktable.view_manager->proxy.lens_calib.get_manual_edit = _get_manual_edit;
  darktable.view_manager->proxy.lens_calib.set_flat = _set_flat;
  darktable.view_manager->proxy.lens_calib.get_flat = _get_flat;
  darktable.view_manager->proxy.lens_calib.set_falloff = _set_falloff;
  darktable.view_manager->proxy.lens_calib.get_falloff = _get_falloff;
  darktable.view_manager->proxy.lens_calib.set_show = _set_show;
  darktable.view_manager->proxy.lens_calib.get_show = _get_show;
  darktable.view_manager->proxy.lens_calib.point_count = _point_count;
  darktable.view_manager->proxy.lens_calib.measured_count = _measured_count;
  darktable.view_manager->proxy.lens_calib.clear_points = _clear_points;
  darktable.view_manager->proxy.lens_calib.get_corner_mode = _get_corner_mode;
  darktable.view_manager->proxy.lens_calib.corner_count = _corner_count;
  darktable.view_manager->proxy.lens_calib.align_grid = _align_grid;
  darktable.view_manager->proxy.lens_calib.has_pose = _has_pose;
  darktable.view_manager->proxy.lens_calib.pose_rms = _pose_rms;
  darktable.view_manager->proxy.lens_calib.mesh_fill = _mesh_fill;
  darktable.view_manager->proxy.lens_calib.layer_has_data = _layer_has_data;
  darktable.view_manager->proxy.lens_calib.interpolated_count = _interpolated_count;
  darktable.view_manager->proxy.lens_calib.stray_count = _stray_count;
  darktable.view_manager->proxy.lens_calib.fit_vignette = _fit_vignette;
  darktable.view_manager->proxy.lens_calib.has_vignette = _has_vignette;
}

// clang-format off
// modelines: These editor modelines have been set for all relevant files by tools/update_modelines.py
// vim: shiftwidth=2 expandtab tabstop=2 cindent
// kate: tab-indents: off; indent-width 2; replace-tabs on; indent-mode cstyle; remove-trailing-spaces modified;
// clang-format on
