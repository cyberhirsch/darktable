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

/* What the centre view shows: which overlays are drawn, and whether the
   frame is shown as shot or undistorted by the current fit. */

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

typedef struct dt_lib_lens_calib_display_t
{
  GtkWidget *undistort;
  GtkWidget *falloff;
  GtkWidget *points, *curves, *residuals, *mesh;
  gboolean updating;
} dt_lib_lens_calib_display_t;

const char *name(dt_lib_module_t *self)
{
  return _("display");
}

dt_view_type_flags_t views(dt_lib_module_t *self)
{
  return DT_VIEW_LENS_CALIB;
}

uint32_t container(dt_lib_module_t *self)
{
  return DT_UI_CONTAINER_PANEL_LEFT_CENTER;
}

int position(const dt_lib_module_t *self)
{
  return 300;
}

static void _undistort_toggled(GtkWidget *w, dt_lib_module_t *self)
{
  dt_lib_lens_calib_display_t *d = self->data;
  if(d->updating) return;

  dt_view_t *view = dt_lens_calib_get_view();
  if(!view || !darktable.view_manager->proxy.lens_calib.set_flat) return;

  const gboolean want_flat = gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(w));
  darktable.view_manager->proxy.lens_calib.set_flat(view, want_flat);

  /* The view refuses to go flat without a fit to flatten with, so read
     back what it actually did rather than assume. */
  const gboolean actual =
    darktable.view_manager->proxy.lens_calib.get_flat
    && darktable.view_manager->proxy.lens_calib.get_flat(view);

  if(actual != want_flat)
  {
    d->updating = TRUE;
    gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(w), actual);
    d->updating = FALSE;
  }
}

static void _falloff_toggled(GtkWidget *w, dt_lib_module_t *self)
{
  dt_lib_lens_calib_display_t *d = self->data;
  if(d->updating) return;

  dt_view_t *view = dt_lens_calib_get_view();
  if(!view || !darktable.view_manager->proxy.lens_calib.set_falloff) return;

  const gboolean want = gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(w));
  darktable.view_manager->proxy.lens_calib.set_falloff(view, want);

  // the view refuses without a measurement, so read back what it really did
  const gboolean actual =
    darktable.view_manager->proxy.lens_calib.get_falloff
    && darktable.view_manager->proxy.lens_calib.get_falloff(view);

  if(actual != want)
  {
    d->updating = TRUE;
    gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(w), actual);
    d->updating = FALSE;
  }
}

static void _overlay_toggled(GtkWidget *w, dt_lib_module_t *self)
{
  dt_lib_lens_calib_display_t *d = self->data;
  if(d->updating) return;

  dt_view_t *view = dt_lens_calib_get_view();
  if(!view || !darktable.view_manager->proxy.lens_calib.set_show) return;

  const int what = GPOINTER_TO_INT(g_object_get_data(G_OBJECT(w), "what"));
  darktable.view_manager->proxy.lens_calib.set_show
    (view, what, gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(w)));
}

static GtkWidget *_overlay_check(dt_lib_module_t *self,
                                 const char *label,
                                 const int what,
                                 const char *tooltip)
{
  GtkWidget *c = gtk_check_button_new_with_label(label);
  g_object_set_data(G_OBJECT(c), "what", GINT_TO_POINTER(what));
  if(tooltip) gtk_widget_set_tooltip_text(c, tooltip);
  g_signal_connect(G_OBJECT(c), "toggled",
                   G_CALLBACK(_overlay_toggled), self);
  return c;
}

void gui_update(dt_lib_module_t *self)
{
  dt_lib_lens_calib_display_t *d = self->data;
  dt_view_t *view = dt_lens_calib_get_view();
  if(!view || !darktable.view_manager->proxy.lens_calib.get_show) return;

  d->updating = TRUE;

  gtk_toggle_button_set_active
    (GTK_TOGGLE_BUTTON(d->points),
     darktable.view_manager->proxy.lens_calib.get_show
       (view, DT_LENS_CALIB_SHOW_POINTS));
  gtk_toggle_button_set_active
    (GTK_TOGGLE_BUTTON(d->curves),
     darktable.view_manager->proxy.lens_calib.get_show
       (view, DT_LENS_CALIB_SHOW_CURVES));
  gtk_toggle_button_set_active
    (GTK_TOGGLE_BUTTON(d->residuals),
     darktable.view_manager->proxy.lens_calib.get_show
       (view, DT_LENS_CALIB_SHOW_RESIDUALS));
  gtk_toggle_button_set_active
    (GTK_TOGGLE_BUTTON(d->mesh),
     darktable.view_manager->proxy.lens_calib.get_show
       (view, DT_LENS_CALIB_SHOW_MESH));

  /* A toggle over an empty layer misrepresents what is on screen. */
  if(darktable.view_manager->proxy.lens_calib.layer_has_data)
  {
    struct { GtkWidget *w; int what; } layer[] = {
      { d->mesh, DT_LENS_CALIB_SHOW_MESH },
      { d->curves, DT_LENS_CALIB_SHOW_CURVES },
      { d->points, DT_LENS_CALIB_SHOW_POINTS },
      { d->residuals, DT_LENS_CALIB_SHOW_RESIDUALS },
    };
    for(size_t i = 0; i < sizeof(layer) / sizeof(layer[0]); i++)
      gtk_widget_set_sensitive
        (layer[i].w,
         darktable.view_manager->proxy.lens_calib.layer_has_data
           (view, layer[i].what));
  }

  if(darktable.view_manager->proxy.lens_calib.get_flat)
    gtk_toggle_button_set_active
      (GTK_TOGGLE_BUTTON(d->undistort),
       darktable.view_manager->proxy.lens_calib.get_flat(view));

  if(darktable.view_manager->proxy.lens_calib.get_falloff)
    gtk_toggle_button_set_active
      (GTK_TOGGLE_BUTTON(d->falloff),
       darktable.view_manager->proxy.lens_calib.get_falloff(view));

  if(darktable.view_manager->proxy.lens_calib.has_vignette)
    gtk_widget_set_sensitive
      (d->falloff,
       darktable.view_manager->proxy.lens_calib.has_vignette(view));

  d->updating = FALSE;
}

void gui_init(dt_lib_module_t *self)
{
  dt_lib_lens_calib_display_t *d = g_malloc0(sizeof(*d));
  self->data = d;

  self->widget = dt_gui_vbox();

  d->undistort = gtk_check_button_new_with_label(_("undistort"));
  gtk_widget_set_tooltip_text
    (d->undistort,
     _("off — the frame as recorded, which is where points are placed.\n"
       "on — the same frame remapped by the current fit, so the chart\n"
       "  lines should come out straight and parallel to the reference\n"
       "  grid, and a squeeze shows as a wider frame.\n"
       "  needs a fit to exist, and is a check on one rather than a\n"
       "  surface to work on: it is rendered at reduced resolution to\n"
       "  stay responsive"));
  g_signal_connect(G_OBJECT(d->undistort), "toggled",
                   G_CALLBACK(_undistort_toggled), self);
  dt_gui_box_add(self->widget, d->undistort);

  d->falloff = gtk_check_button_new_with_label(_("correct falloff"));
  gtk_widget_set_tooltip_text
    (d->falloff,
     _("even out the brightness using the measured vignetting.\n"
       "independent of the geometry above: this one moves no pixels and\n"
       "changes only their brightness, so the two can be judged\n"
       "separately.\n"
       "needs a vignetting measurement. the preview applies the gain in an\n"
       "approximate linear space, so treat it as a look rather than as the\n"
       "final numbers"));
  g_signal_connect(G_OBJECT(d->falloff), "toggled",
                   G_CALLBACK(_falloff_toggled), self);
  dt_gui_box_add(self->widget, d->falloff);

  d->mesh = _overlay_check
    (self, _("grid"), DT_LENS_CALIB_SHOW_MESH,
     _("the points joined to their lattice neighbours.\n"
       "this one does pass through every point, because it is the points.\n"
       "use it to check the lattice itself: a node with the wrong index\n"
       "shows up as a zigzag or a line doubling back, which is invisible\n"
       "in a field of loose crosses"));

  d->curves = _overlay_check
    (self, _("traced lines"), DT_LENS_CALIB_SHOW_CURVES,
     _("the polylines the detector followed. the honest readout of\n"
       "whether it latched onto chart lines or onto noise"));

  d->points = _overlay_check
    (self, _("points"), DT_LENS_CALIB_SHOW_POINTS,
     _("the measured nodes, as green crosses"));

  d->residuals = _overlay_check
    (self, _("residuals"), DT_LENS_CALIB_SHOW_RESIDUALS,
     _("how far each point still is from the line it should sit on,\n"
       "after the fit. drawn exaggerated, since a good fit leaves\n"
       "well under a pixel and a true-to-scale arrow would be invisible"));

  dt_gui_box_add(self->widget, d->mesh, d->curves, d->points, d->residuals);
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
