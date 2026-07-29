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

/* Chart geometry, automatic detection, and hand editing of the points.
 *
 * Detection and hand editing are two separate systems that meet in one
 * place: the "edit points by hand" switch, which adopts whatever the
 * detector found into the editable set. Keeping them apart the rest of
 * the time is what stops a re-detect from throwing away work.
 */

#include "bauhaus/bauhaus.h"
#include "common/darktable.h"
#include "control/conf.h"
#include "control/control.h"
#include "gui/accelerators.h"
#include "gui/gtk.h"
#include "libs/lens_calib_common.h"
#include "libs/lib.h"
#include "libs/lib_api.h"
#include "views/view.h"

DT_MODULE(1)

typedef struct dt_lib_lens_calib_chart_t
{
  GtkWidget *cells_x, *cells_y;
  GtkWidget *detect_button;
  GtkWidget *status;

  GtkWidget *manual_edit;
  GtkWidget *clear_button;

  GtkWidget *align_button;
  GtkWidget *fill_button;
  GtkWidget *use_interpolated;

  gboolean updating;
} dt_lib_lens_calib_chart_t;

const char *name(dt_lib_module_t *self)
{
  return _("chart and points");
}

dt_view_type_flags_t views(dt_lib_module_t *self)
{
  return DT_VIEW_LENS_CALIB;
}

uint32_t container(dt_lib_module_t *self)
{
  return DT_UI_CONTAINER_PANEL_RIGHT_CENTER;
}

int position(const dt_lib_module_t *self)
{
  return 400;
}

static void _refresh(dt_lib_module_t *self)
{
  dt_lib_lens_calib_chart_t *d = self->data;
  dt_view_t *view = dt_lens_calib_get_view();
  if(!view) return;

  const int count = darktable.view_manager->proxy.lens_calib.point_count
    ? darktable.view_manager->proxy.lens_calib.point_count(view) : 0;
  const int measured = darktable.view_manager->proxy.lens_calib.measured_count
    ? darktable.view_manager->proxy.lens_calib.measured_count(view) : 0;
  const int corners = darktable.view_manager->proxy.lens_calib.corner_count
    ? darktable.view_manager->proxy.lens_calib.corner_count(view) : 0;

  const gboolean editing =
    darktable.view_manager->proxy.lens_calib.get_manual_edit
    && darktable.view_manager->proxy.lens_calib.get_manual_edit(view);

  const gboolean pose =
    darktable.view_manager->proxy.lens_calib.has_pose
    && darktable.view_manager->proxy.lens_calib.has_pose(view);

  const int interp = darktable.view_manager->proxy.lens_calib.interpolated_count
    ? darktable.view_manager->proxy.lens_calib.interpolated_count(view) : 0;

  const int stray = darktable.view_manager->proxy.lens_calib.stray_count
    ? darktable.view_manager->proxy.lens_calib.stray_count(view) : 0;

  const gboolean asking =
    darktable.view_manager->proxy.lens_calib.get_corner_mode
    && darktable.view_manager->proxy.lens_calib.get_corner_mode(view);

  // the lattice is complete when every site has a node
  const int cx = dt_conf_get_int("plugins/lens_calib/cells_x");
  const int cy = dt_conf_get_int("plugins/lens_calib/cells_y");
  const int sites = (cx + 1) * (cy + 1);

  gtk_widget_set_sensitive(d->clear_button, count > 0);
  gtk_widget_set_sensitive(d->align_button, corners == 4 || measured >= 4);
  gtk_widget_set_sensitive(d->fill_button, count < sites);
  gtk_widget_set_sensitive(d->use_interpolated, interp > 0);

  GString *text = g_string_new(NULL);

  if(asking)
    g_string_append_printf(text,
                           _("click the chart corners: %d of 4\n"), corners);

  if(count)
  {
    /* What the fit will actually use, first. A measured point off the
       lattice, or a guessed one, is not evidence about the lens. */
    g_string_append_printf(text, _("%d of %d nodes measured"),
                           measured - stray, sites);
    if(interp) g_string_append_printf(text, _(", %d guessed"), interp);
    if(stray)
      g_string_append_printf(text, ngettext(", %d stray", ", %d strays", stray),
                             stray);
    g_string_append_c(text, '\n');
  }
  else
    g_string_append(text, _("no points yet\n"));

  if(pose && darktable.view_manager->proxy.lens_calib.pose_rms)
    g_string_append_printf(text, _("%.1f px of distortion to explain"),
                           darktable.view_manager->proxy.lens_calib.pose_rms(view));
  else
    g_string_append(text, _("chart position not known yet"));

  if(!editing && count)
    g_string_append(text, _("\nswitch on hand editing to change them"));

  gtk_label_set_text(GTK_LABEL(d->status), text->str);
  g_string_free(text, TRUE);
}

static void _cells_changed(GtkWidget *w, dt_lib_module_t *self)
{
  dt_lib_lens_calib_chart_t *d = self->data;
  if(d->updating) return;

  const int cx = gtk_spin_button_get_value_as_int(GTK_SPIN_BUTTON(d->cells_x));
  const int cy = gtk_spin_button_get_value_as_int(GTK_SPIN_BUTTON(d->cells_y));
  dt_conf_set_int("plugins/lens_calib/cells_x", cx);
  dt_conf_set_int("plugins/lens_calib/cells_y", cy);

  dt_view_t *view = dt_lens_calib_get_view();
  if(view && darktable.view_manager->proxy.lens_calib.chart_changed)
    darktable.view_manager->proxy.lens_calib.chart_changed(view);

  _refresh(self);
}

static void _detect_clicked(GtkWidget *w, dt_lib_module_t *self)
{
  dt_view_t *view = dt_lens_calib_get_view();
  if(view && darktable.view_manager->proxy.lens_calib.detect)
    darktable.view_manager->proxy.lens_calib.detect(view);

  // detection gives the traced-lines overlay something to draw
  dt_lens_calib_refresh_panels();
  _refresh(self);
}

static void _manual_edit_toggled(GtkWidget *w, dt_lib_module_t *self)
{
  dt_lib_lens_calib_chart_t *d = self->data;
  if(d->updating) return;

  dt_view_t *view = dt_lens_calib_get_view();
  if(view && darktable.view_manager->proxy.lens_calib.set_manual_edit)
    darktable.view_manager->proxy.lens_calib.set_manual_edit
      (view, gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(w)));

  _refresh(self);
}

static void _clear_clicked(GtkWidget *w, dt_lib_module_t *self)
{
  dt_view_t *view = dt_lens_calib_get_view();
  if(view && darktable.view_manager->proxy.lens_calib.clear_points)
    darktable.view_manager->proxy.lens_calib.clear_points(view);

  _refresh(self);
}

static void _align_clicked(GtkWidget *w, dt_lib_module_t *self)
{
  dt_view_t *view = dt_lens_calib_get_view();
  if(view && darktable.view_manager->proxy.lens_calib.align_grid)
    darktable.view_manager->proxy.lens_calib.align_grid(view);

  _refresh(self);
}

static void _fill_clicked(GtkWidget *w, dt_lib_module_t *self)
{
  dt_view_t *view = dt_lens_calib_get_view();
  if(view && darktable.view_manager->proxy.lens_calib.mesh_fill)
    darktable.view_manager->proxy.lens_calib.mesh_fill(view);

  _refresh(self);
}

static void _use_interpolated_toggled(GtkWidget *w, dt_lib_module_t *self)
{
  dt_lib_lens_calib_chart_t *d = self->data;
  if(d->updating) return;

  dt_conf_set_bool("plugins/lens_calib/use_interpolated",
                   gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(w)));
}

void gui_update(dt_lib_module_t *self)
{
  dt_lib_lens_calib_chart_t *d = self->data;

  /* Seeding one spin button fires the shared handler, which reads both --
     and the other is still at its range minimum at that moment, so it
     would write a 1 back over the real value. */
  d->updating = TRUE;

  gtk_spin_button_set_value(GTK_SPIN_BUTTON(d->cells_x),
                            dt_conf_get_int("plugins/lens_calib/cells_x"));
  gtk_spin_button_set_value(GTK_SPIN_BUTTON(d->cells_y),
                            dt_conf_get_int("plugins/lens_calib/cells_y"));

  dt_view_t *view = dt_lens_calib_get_view();
  if(view && darktable.view_manager->proxy.lens_calib.get_manual_edit)
    gtk_toggle_button_set_active
      (GTK_TOGGLE_BUTTON(d->manual_edit),
       darktable.view_manager->proxy.lens_calib.get_manual_edit(view));

  gtk_toggle_button_set_active
    (GTK_TOGGLE_BUTTON(d->use_interpolated),
     dt_conf_get_bool("plugins/lens_calib/use_interpolated"));

  d->updating = FALSE;

  _refresh(self);
}

void gui_init(dt_lib_module_t *self)
{
  dt_lib_lens_calib_chart_t *d = g_malloc0(sizeof(*d));
  self->data = d;

  self->widget = dt_gui_vbox();

  GtkWidget *row_x =
    dt_lens_calib_labelled_spin(_("cells across"), &d->cells_x, 1, 200);
  GtkWidget *row_y =
    dt_lens_calib_labelled_spin(_("cells down"), &d->cells_y, 1, 200);

  gtk_widget_set_tooltip_text
    (d->cells_x,
     _("number of cells across the chart, counting cells not lines\n"
       "a 28x12 chart has a 29x13 lattice of intersections"));
  gtk_widget_set_tooltip_text(d->cells_y,
                              _("number of cells down the chart"));

  dt_gui_box_add(self->widget, row_x, row_y);

  g_signal_connect(G_OBJECT(d->cells_x), "value-changed",
                   G_CALLBACK(_cells_changed), self);
  g_signal_connect(G_OBJECT(d->cells_y), "value-changed",
                   G_CALLBACK(_cells_changed), self);

  d->detect_button = dt_action_button_new
    (self, N_("detect grid"), _detect_clicked, self,
     _("trace the chart lines and derive their intersections"), 0,
     (GdkModifierType)0);
  dt_gui_box_add(self->widget, d->detect_button);

  d->manual_edit =
    gtk_check_button_new_with_label(_("edit points by hand"));
  gtk_widget_set_tooltip_text
    (d->manual_edit,
     _("take over the detected points for editing, and allow new ones.\n"
       "left drag moves a point, left click on bare image adds one,\n"
       "right click removes one.\n"
       "this is one way: the points become yours and a fresh detection\n"
       "would be needed to get the automatic ones back"));
  g_signal_connect(G_OBJECT(d->manual_edit), "toggled",
                   G_CALLBACK(_manual_edit_toggled), self);
  dt_gui_box_add(self->widget, d->manual_edit);

  d->status = gtk_label_new("");
  gtk_widget_set_halign(d->status, GTK_ALIGN_START);
  gtk_label_set_ellipsize(GTK_LABEL(d->status), PANGO_ELLIPSIZE_END);
  dt_gui_box_add(self->widget, d->status);

  d->clear_button = dt_action_button_new
    (self, N_("clear points"), _clear_clicked, self,
     _("discard every placed point and the fit made from them"), 0,
     (GdkModifierType)0);
  dt_gui_box_add(self->widget, d->clear_button);

  d->align_button = dt_action_button_new
    (self, N_("re-align grid to points"), _align_clicked, self,
     _("work out where the chart sits in the frame, and give every point\n"
       "its lattice index. this runs by itself whenever the points change,\n"
       "so it is only needed after editing a great many nodes.\n"
       "what is left over afterwards is the distortion: a straight grid\n"
       "cannot pass through points that a lens has bent"), 0,
     (GdkModifierType)0);
  dt_gui_box_add(self->widget, d->align_button);

  d->fill_button = dt_action_button_new
    (self, N_("fill in missing nodes"), _fill_clicked, self,
     _("guess the nodes detection did not find, from the ones it did.\n"
       "a fallback for a chart that came out badly -- guessed nodes are\n"
       "not measurements and the fit ignores them, so this buys you\n"
       "something to drag rather than any extra evidence.\n"
       "with nothing to go on it will ask you to click the four corners"), 0,
     (GdkModifierType)0);
  dt_gui_box_add(self->widget, d->fill_button);

  d->use_interpolated =
    gtk_check_button_new_with_label(_("fit interpolated nodes too"));
  gtk_widget_set_tooltip_text
    (d->use_interpolated,
     _("normally the fit uses only measured nodes.\n"
       "an interpolated node was computed from the assumption that the\n"
       "lattice is smooth, so fitting it is circular: it lowers the\n"
       "reported residual while making the calibration worse.\n"
       "only worth enabling for a mesh that is nearly complete"));
  g_signal_connect(G_OBJECT(d->use_interpolated), "toggled",
                   G_CALLBACK(_use_interpolated_toggled), self);
  dt_gui_box_add(self->widget, d->use_interpolated);

  _refresh(self);
}

void gui_cleanup(dt_lib_module_t *self)
{
  g_free(self->data);
  self->data = NULL;
}

// clang-format off
// modelines: These editor modelines have been set for all relevant files by tools/update_modelines.py
// vim: shiftwidth=2 expandtab tabstop=2 cindent
// kate: tab-indents: off; indent-width 2; replace-tabs on; indent-mode cstyle; remove-trailing-spaces modified;
// clang-format on
