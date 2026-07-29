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

/* Fitting a lens model to the placed points.
 *
 * Every option here lives in conf rather than in this module's data, so
 * the view can read them at solve time without the two sharing a header.
 */

#include "bauhaus/bauhaus.h"
#include "common/darktable.h"
#include "common/lens_warp.h"
#include "control/conf.h"
#include "control/control.h"
#include "gui/accelerators.h"
#include "gui/gtk.h"
#include "libs/lens_calib_common.h"
#include "libs/lib.h"
#include "libs/lib_api.h"
#include "views/view.h"

DT_MODULE(1)

typedef struct dt_lib_lens_calib_fit_t
{
  GtkWidget *model, *order, *regular;
  GtkWidget *squeeze, *measure_squeeze;
  GtkWidget *solve_centre, *chart_frontal;
  GtkWidget *solve_button;
  GtkWidget *vignette_button;
  GtkWidget *status;
  GtkWidget *ax_focal, *ax_aperture, *ax_distance;
  GtkWidget *entry_add, *entry_rm;
  float sel_focal, sel_aperture, sel_distance;
  GtkWidget *edit_values, *values_box;
  GPtrArray *value_spins;
  gboolean updating;
} dt_lib_lens_calib_fit_t;

const char *name(dt_lib_module_t *self)
{
  return _("fit");
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
  return 200;
}

static void _refresh(dt_lib_module_t *self)
{
  dt_lib_lens_calib_fit_t *d = self->data;
  dt_view_t *view = dt_lens_calib_get_view();

  const gboolean has_points =
    view && darktable.view_manager->proxy.lens_calib.has_points
    && darktable.view_manager->proxy.lens_calib.has_points(view);

  gtk_widget_set_sensitive(d->solve_button, has_points);

  char *text = (view && darktable.view_manager->proxy.lens_calib.status_text)
    ? darktable.view_manager->proxy.lens_calib.status_text(view) : NULL;

  gtk_label_set_text(GTK_LABEL(d->status), text ? text : _("not fitted"));
  g_free(text);

  /* The polynomial order means nothing to the anamorphic model, which has
     a fixed set of radial terms, so it is disabled rather than left there
     looking as though it did something. */
  gtk_widget_set_sensitive
    (d->order, dt_bauhaus_combobox_get(d->model) != DT_LENS_WARP_ANAM_RADIAL);

  // a declared ratio and a measured one are alternatives, not both
  gtk_widget_set_sensitive
    (d->squeeze,
     !gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(d->measure_squeeze)));
}

static void _options_changed(GtkWidget *w, dt_lib_module_t *self)
{
  dt_lib_lens_calib_fit_t *d = self->data;
  if(d->updating) return;

  dt_conf_set_int("plugins/lens_calib/model",
                  dt_bauhaus_combobox_get(d->model));
  dt_conf_set_int("plugins/lens_calib/order",
                  gtk_spin_button_get_value_as_int(GTK_SPIN_BUTTON(d->order)));
  dt_conf_set_float("plugins/lens_calib/regularization",
                    gtk_spin_button_get_value(GTK_SPIN_BUTTON(d->regular)));
  dt_conf_set_bool("plugins/lens_calib/solve_centre",
                   gtk_toggle_button_get_active
                     (GTK_TOGGLE_BUTTON(d->solve_centre)));
  dt_conf_set_bool("plugins/lens_calib/chart_frontal",
                   gtk_toggle_button_get_active
                     (GTK_TOGGLE_BUTTON(d->chart_frontal)));
  dt_conf_set_float("plugins/lens_calib/squeeze",
                    gtk_spin_button_get_value(GTK_SPIN_BUTTON(d->squeeze)));
  dt_conf_set_bool("plugins/lens_calib/measure_squeeze",
                   gtk_toggle_button_get_active
                     (GTK_TOGGLE_BUTTON(d->measure_squeeze)));

  _refresh(self);
}

static void _solve_clicked(GtkWidget *w, dt_lib_module_t *self)
{
  dt_view_t *view = dt_lens_calib_get_view();
  if(view && darktable.view_manager->proxy.lens_calib.solve)
    darktable.view_manager->proxy.lens_calib.solve(view);

  // a fit enables the undistorted preview and the residual overlay elsewhere
  dt_lens_calib_refresh_panels();
  _refresh(self);
}

static void _vignette_clicked(GtkWidget *w, dt_lib_module_t *self)
{
  dt_view_t *view = dt_lens_calib_get_view();
  if(view && darktable.view_manager->proxy.lens_calib.fit_vignette)
    darktable.view_manager->proxy.lens_calib.fit_vignette(view);

  // a measurement enables the falloff preview, which lives in another module
  dt_lens_calib_refresh_panels();
  _refresh(self);
}



#define DT_LENS_CALIB_MAX_AXIS 512

static void _rebuild_axes(dt_lib_module_t *self);

/* Point the session at the entry the three selectors now describe. */
static void _apply_axis_selection(dt_lib_module_t *self)
{
  dt_lib_lens_calib_fit_t *d = self->data;
  dt_view_t *view = dt_lens_calib_get_view();
  if(!view || !darktable.view_manager->proxy.lens_calib.entry_find) return;

  const int i = darktable.view_manager->proxy.lens_calib.entry_find
    (view, d->sel_focal, d->sel_aperture, d->sel_distance);

  if(i >= 0 && darktable.view_manager->proxy.lens_calib.entry_select)
    darktable.view_manager->proxy.lens_calib.entry_select(view, i);
}

static void _axis_changed(GtkWidget *w, dt_lib_module_t *self)
{
  dt_lib_lens_calib_fit_t *d = self->data;
  if(d->updating) return;

  /* Changing a coarser axis can invalidate the finer ones -- a focal length
     may not offer the aperture that was selected -- so everything below the
     one touched is rebuilt before the selection is applied. */
  d->updating = TRUE;
  _rebuild_axes(self);
  d->updating = FALSE;

  _apply_axis_selection(self);
}

/* Fill one selector, and return the value that ended up chosen.
 *
 * Keeps the previous value when the new list still offers it, so stepping
 * through focal lengths at a fixed aperture does not reset the aperture on
 * every step. */
static float _fill_axis(GtkWidget *combo,
                        dt_view_t *view,
                        const int axis,
                        const float focal,
                        const float aperture,
                        const float want,
                        const char *fmt,
                        const char *any_label,
                        int *out_n)
{
  float vals[DT_LENS_CALIB_MAX_AXIS];
  const int n = darktable.view_manager->proxy.lens_calib.axis_values
    (view, axis, focal, aperture, vals, DT_LENS_CALIB_MAX_AXIS);

  dt_bauhaus_combobox_clear(combo);
  *out_n = n;

  if(n == 0)
  {
    dt_bauhaus_combobox_add(combo, "-");
    gtk_widget_set_sensitive(combo, FALSE);
    return 0.0f;
  }

  int pick = 0;
  for(int i = 0; i < n; i++)
  {
    /* A zero on the aperture or distance axis is not a measurement at
       f/0 -- it is the entry that applies at any of them, which is where
       distortion and chromatic aberration live since neither depends on
       these axes. Saying "any" is the literal meaning of the stored
       value, not a convenience. */
    if(vals[i] <= 0.0f && any_label)
      dt_bauhaus_combobox_add(combo, any_label);
    else
    {
      gchar *txt = g_strdup_printf(fmt, (double)vals[i]);
      dt_bauhaus_combobox_add(combo, txt);
      g_free(txt);
    }
    if(fabsf(vals[i] - want) < 1e-3f) pick = i;
  }

  gtk_widget_set_sensitive(combo, n > 1);
  dt_bauhaus_combobox_set(combo, pick);
  return vals[pick];
}

/* Rebuild all three selectors, coarse to fine. */
static void _rebuild_axes(dt_lib_module_t *self)
{
  dt_lib_lens_calib_fit_t *d = self->data;
  dt_view_t *view = dt_lens_calib_get_view();

  if(!view || !darktable.view_manager->proxy.lens_calib.axis_values)
    return;

  int nf = 0, na = 0, nd = 0;

  d->sel_focal = _fill_axis(d->ax_focal, view, 0, 0.0f, 0.0f,
                            d->sel_focal, "%g mm", NULL, &nf);
  d->sel_aperture = _fill_axis(d->ax_aperture, view, 1, d->sel_focal, 0.0f,
                               d->sel_aperture, "f/%g", _("any"), &na);
  d->sel_distance = _fill_axis(d->ax_distance, view, 2, d->sel_focal,
                               d->sel_aperture, d->sel_distance,
                               "%g m", _("any"), &nd);

  /* Focus distance is a Lensfun nicety: a profile measured here has one
     value for it, so the row would be a permanently disabled control
     explaining nothing. It appears only when there is a choice. */
  gtk_widget_set_visible(d->ax_distance, nd > 1);

  gtk_widget_set_sensitive(d->entry_rm, nf > 0);
}

static void _entry_add_clicked(GtkWidget *w, dt_lib_module_t *self)
{
  dt_view_t *view = dt_lens_calib_get_view();
  if(view && darktable.view_manager->proxy.lens_calib.entry_add)
    darktable.view_manager->proxy.lens_calib.entry_add(view);
}

static void _entry_rm_clicked(GtkWidget *w, dt_lib_module_t *self)
{
  dt_lib_lens_calib_fit_t *d = self->data;
  dt_view_t *view = dt_lens_calib_get_view();
  if(!view || !darktable.view_manager->proxy.lens_calib.entry_find) return;

  const int i = darktable.view_manager->proxy.lens_calib.entry_find
    (view, d->sel_focal, d->sel_aperture, d->sel_distance);

  if(i >= 0 && darktable.view_manager->proxy.lens_calib.entry_remove)
    darktable.view_manager->proxy.lens_calib.entry_remove(view, i);
}

static void _value_spin_changed(GtkWidget *w, dt_lib_module_t *self)
{
  dt_lib_lens_calib_fit_t *d = self->data;
  if(d->updating) return;

  dt_view_t *view = dt_lens_calib_get_view();
  if(!view || !darktable.view_manager->proxy.lens_calib.value_set) return;

  const int i = GPOINTER_TO_INT(g_object_get_data(G_OBJECT(w), "value-index"));
  darktable.view_manager->proxy.lens_calib.value_set
    (view, i, gtk_spin_button_get_value(GTK_SPIN_BUTTON(w)));
}

/* Draw one row per value the fit actually has.
 *
 * Rebuilt rather than updated in place because the set itself changes: a
 * different model has different coefficients, and vignetting rows only exist
 * once a falloff has been measured. A row that lingers after the thing it
 * edits has gone is worse than no row.
 */
static void _rebuild_values(dt_lib_module_t *self)
{
  dt_lib_lens_calib_fit_t *d = self->data;
  dt_view_t *view = dt_lens_calib_get_view();

  for(guint i = 0; i < d->value_spins->len; i++)
    gtk_widget_destroy(GTK_WIDGET(g_ptr_array_index(d->value_spins, i)));
  g_ptr_array_set_size(d->value_spins, 0);

  const int n = (view && darktable.view_manager->proxy.lens_calib.value_count)
    ? darktable.view_manager->proxy.lens_calib.value_count(view) : 0;

  const gboolean editable =
    gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(d->edit_values));

  for(int i = 0; i < n; i++)
  {
    const char *nm =
      darktable.view_manager->proxy.lens_calib.value_name(view, i);
    const double v =
      darktable.view_manager->proxy.lens_calib.value_get(view, i);

    /* Wide range and many digits on purpose: these are raw model
       coefficients, not tuned settings, and clamping them to something
       tidy would silently refuse to represent a legitimate fit. */
    GtkWidget *spin = gtk_spin_button_new_with_range(-1e6, 1e6, 1e-5);
    gtk_spin_button_set_digits(GTK_SPIN_BUTTON(spin), 6);
    gtk_spin_button_set_value(GTK_SPIN_BUTTON(spin), v);
    gtk_widget_set_sensitive(spin, editable);
    g_object_set_data(G_OBJECT(spin), "value-index", GINT_TO_POINTER(i));
    g_signal_connect(G_OBJECT(spin), "value-changed",
                     G_CALLBACK(_value_spin_changed), self);

    GtkWidget *row = dt_lens_calib_labelled(_(nm), spin);
    gtk_box_pack_start(GTK_BOX(d->values_box), row, FALSE, FALSE, 0);
    g_ptr_array_add(d->value_spins, row);
  }

  gtk_widget_show_all(d->values_box);
}

static void _edit_values_toggled(GtkWidget *w, dt_lib_module_t *self)
{
  dt_lib_lens_calib_fit_t *d = self->data;
  if(d->updating) return;

  dt_conf_set_bool("plugins/lens_calib/edit_values",
                   gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(w)));
  _rebuild_values(self);
}

void gui_update(dt_lib_module_t *self)
{
  dt_lib_lens_calib_fit_t *d = self->data;

  d->updating = TRUE;

  dt_bauhaus_combobox_set(d->model,
                          dt_conf_get_int("plugins/lens_calib/model"));
  gtk_spin_button_set_value(GTK_SPIN_BUTTON(d->order),
                            dt_conf_get_int("plugins/lens_calib/order"));
  gtk_spin_button_set_value
    (GTK_SPIN_BUTTON(d->regular),
     dt_conf_get_float("plugins/lens_calib/regularization"));
  gtk_toggle_button_set_active
    (GTK_TOGGLE_BUTTON(d->solve_centre),
     dt_conf_get_bool("plugins/lens_calib/solve_centre"));
  gtk_toggle_button_set_active
    (GTK_TOGGLE_BUTTON(d->chart_frontal),
     dt_conf_get_bool("plugins/lens_calib/chart_frontal"));
  gtk_spin_button_set_value(GTK_SPIN_BUTTON(d->squeeze),
                            dt_conf_get_float("plugins/lens_calib/squeeze"));
  gtk_toggle_button_set_active
    (GTK_TOGGLE_BUTTON(d->measure_squeeze),
     dt_conf_get_bool("plugins/lens_calib/measure_squeeze"));
  gtk_toggle_button_set_active
    (GTK_TOGGLE_BUTTON(d->edit_values),
     dt_conf_get_bool("plugins/lens_calib/edit_values"));

  _rebuild_axes(self);
  _rebuild_values(self);

  d->updating = FALSE;

  _refresh(self);
}

void gui_init(dt_lib_module_t *self)
{
  dt_lib_lens_calib_fit_t *d = g_malloc0(sizeof(*d));
  self->data = d;

  self->widget = dt_gui_vbox();

  d->model = dt_bauhaus_combobox_new(NULL);
  dt_bauhaus_widget_set_label(d->model, NULL, N_("model"));
  dt_bauhaus_combobox_add(d->model, _("polynomial"));
  dt_bauhaus_combobox_add(d->model, _("anamorphic radial"));
  dt_bauhaus_combobox_add(d->model, _("spline"));
  gtk_widget_set_tooltip_text
    (d->model,
     _("polynomial — free two dimensional fit, no symmetry assumed\n"
       "anamorphic radial — few parameters, extrapolates sanely,\n"
       "  the right shape for a spherical lens behind an anamorphic front\n"
       "spline — interpolates the leftover residual exactly, which also\n"
       "  means it reproduces measurement noise"));
  g_signal_connect(G_OBJECT(d->model), "value-changed",
                   G_CALLBACK(_options_changed), self);
  dt_gui_box_add(self->widget, d->model);

  GtkWidget *order_row = dt_lens_calib_labelled_spin
    (_("order"), &d->order, 2, DT_LENS_WARP_MAX_ORDER);
  gtk_widget_set_tooltip_text
    (d->order,
     _("highest polynomial degree the fit may use.\n"
       "higher follows the measurements more closely and\n"
       "misbehaves further outside them"));
  g_signal_connect(G_OBJECT(d->order), "value-changed",
                   G_CALLBACK(_options_changed), self);
  dt_gui_box_add(self->widget, order_row);

  d->regular = gtk_spin_button_new_with_range(0.0, 0.1, 0.0001);
  gtk_spin_button_set_digits(GTK_SPIN_BUTTON(d->regular), 5);
  gtk_widget_set_tooltip_text
    (d->regular,
     _("how hard the fit is pulled towards a simple model.\n"
       "raise it if the correction goes wild near the frame edges,\n"
       "where there is least data to constrain it"));
  g_signal_connect(G_OBJECT(d->regular), "value-changed",
                   G_CALLBACK(_options_changed), self);
  dt_gui_box_add(self->widget,
                 dt_lens_calib_labelled(_("regularization"), d->regular));

  /* The squeeze cannot be fitted. A pure scale along one axis keeps straight
     lines straight, so the plumb-line objective is blind to it by
     construction -- it has to be either declared or inferred from cell
     shape. Declaring it is the better of the two: the number is engraved on
     the lens, whereas cell shape needs the chart to have been square on and
     cannot tell foreshortening from squeeze. */
  d->squeeze = gtk_spin_button_new_with_range(0.5, 4.0, 0.01);
  gtk_spin_button_set_digits(GTK_SPIN_BUTTON(d->squeeze), 3);
  gtk_widget_set_tooltip_text
    (d->squeeze,
     _("the lens's squeeze factor: 1.6 for a 1.6x anamorphic,\n"
       "2.0 for a 2x, 1.0 for a spherical lens.\n"
       "this is what the correction desqueezes by, so getting it wrong\n"
       "leaves the output the wrong shape however good the fit is"));
  g_signal_connect(G_OBJECT(d->squeeze), "value-changed",
                   G_CALLBACK(_options_changed), self);
  dt_gui_box_add(self->widget,
                 dt_lens_calib_labelled(_("squeeze"), d->squeeze));

  d->measure_squeeze =
    gtk_check_button_new_with_label(_("measure squeeze from the chart"));
  gtk_widget_set_tooltip_text
    (d->measure_squeeze,
     _("work the ratio out from the shape of the cells instead of being\n"
       "told it. needs the chart to have been shot square on, since the\n"
       "foreshortening of a tilted chart is indistinguishable from a\n"
       "squeeze. the measurement is reported either way, so leave this\n"
       "off and use it to check the value you entered"));
  g_signal_connect(G_OBJECT(d->measure_squeeze), "toggled",
                   G_CALLBACK(_options_changed), self);
  dt_gui_box_add(self->widget, d->measure_squeeze);

  d->solve_centre = gtk_check_button_new_with_label(_("fit optical centre"));
  gtk_widget_set_tooltip_text
    (d->solve_centre,
     _("let the distortion centre move away from the frame centre,\n"
       "which real lenses generally need"));
  g_signal_connect(G_OBJECT(d->solve_centre), "toggled",
                   G_CALLBACK(_options_changed), self);

  d->chart_frontal =
    gtk_check_button_new_with_label(_("chart shot square on"));
  gtk_widget_set_tooltip_text
    (d->chart_frontal,
     _("assert that the camera faced the chart squarely.\n"
       "this does two things: it recovers the anamorphic squeeze from\n"
       "the shape of the cells, and it stops the fit from absorbing a\n"
       "keystone that was never in the lens -- straight lines stay\n"
       "straight under perspective, so nothing else can rule that out.\n"
       "turn it off only for a chart you know was tilted"));
  g_signal_connect(G_OBJECT(d->chart_frontal), "toggled",
                   G_CALLBACK(_options_changed), self);

  dt_gui_box_add(self->widget, d->solve_centre, d->chart_frontal);

  d->solve_button = dt_action_button_new
    (self, N_("fit lens model"), _solve_clicked, self,
     _("fit the model to the placed points by making the chart lines"
       " straight"), 0, (GdkModifierType)0);
  dt_gui_box_add(self->widget, d->solve_button);

  d->vignette_button = dt_action_button_new
    (self, N_("measure vignetting"), _vignette_clicked, self,
     _("measure the brightness falloff towards the corners, from this\n"
       "frame, and store it in the profile.\n"
       "needs no points: it reads the image itself. what it does need is\n"
       "an evenly lit subject, because a wall that is brighter in the\n"
       "middle is indistinguishable from a lens that is darker at the\n"
       "edges -- the chart shot usually serves, since the dark lines are\n"
       "rejected rather than averaged in.\n"
       "vignetting depends on aperture, so the f-number is recorded with\n"
       "it and a profile describes the aperture it was measured at"), 0,
     (GdkModifierType)0);
  dt_gui_box_add(self->widget, d->vignette_button);

  d->status = gtk_label_new(_("not fitted"));
  gtk_widget_set_halign(d->status, GTK_ALIGN_START);
  gtk_label_set_xalign(GTK_LABEL(d->status), 0.0);
  gtk_label_set_selectable(GTK_LABEL(d->status), TRUE);
  dt_gui_box_add(self->widget, d->status);

  /* --- the measured values --- */
  d->value_spins = g_ptr_array_new();

  /* --- which measurement --- */
  d->ax_focal = dt_bauhaus_combobox_new(NULL);
  dt_bauhaus_widget_set_label(d->ax_focal, NULL, N_("focal length"));
  gtk_widget_set_tooltip_text
    (d->ax_focal,
     _("the focal lengths this profile holds a measurement at"));
  g_signal_connect(G_OBJECT(d->ax_focal), "value-changed",
                   G_CALLBACK(_axis_changed), self);

  d->ax_aperture = dt_bauhaus_combobox_new(NULL);
  dt_bauhaus_widget_set_label(d->ax_aperture, NULL, N_("aperture"));
  gtk_widget_set_tooltip_text
    (d->ax_aperture,
     _("apertures measured at this focal length.\n\n"
       "any is where distortion and chromatic aberration live: neither\n"
       "depends on aperture, so they are stored once for all of them"));
  g_signal_connect(G_OBJECT(d->ax_aperture), "value-changed",
                   G_CALLBACK(_axis_changed), self);

  d->ax_distance = dt_bauhaus_combobox_new(NULL);
  dt_bauhaus_widget_set_label(d->ax_distance, NULL, N_("focus distance"));
  gtk_widget_set_tooltip_text
    (d->ax_distance,
     _("focus distances measured at this focal length and aperture.\n"
       "hidden unless the profile actually varies it"));
  g_signal_connect(G_OBJECT(d->ax_distance), "value-changed",
                   G_CALLBACK(_axis_changed), self);

  d->entry_add = dt_action_button_new
    (self, N_("add measurement"), _entry_add_clicked, self,
     _("start another entry at the focal length and aperture now set\n"
       "in the lens panel, seeded from the current fit"), 0,
     (GdkModifierType)0);

  d->entry_rm = dt_action_button_new
    (self, N_("remove measurement"), _entry_rm_clicked, self,
     _("drop the selected entry from the profile"), 0, (GdkModifierType)0);

  dt_gui_box_add(self->widget, d->ax_focal, d->ax_aperture,
                 d->ax_distance,
                 dt_gui_hbox(dt_gui_expand(d->entry_add), d->entry_rm));

  d->edit_values = gtk_check_button_new_with_label(_("edit values"));
  gtk_widget_set_tooltip_text
    (d->edit_values,
     _("unlock the fitted coefficients so they can be typed over,\n"
       "including any imported from Lensfun.\n\n"
       "editing one stops the profile being a measurement: the\n"
       "residual and point count it carries describe the fit that\n"
       "produced these numbers, and once one is changed by hand they\n"
       "no longer do. the file records that, so a hand made profile\n"
       "is not counted as evidence alongside measured ones"));
  g_signal_connect(G_OBJECT(d->edit_values), "toggled",
                   G_CALLBACK(_edit_values_toggled), self);

  d->values_box = dt_gui_vbox();

  GtkWidget *values_exp = gtk_expander_new(_("measured values"));
  gtk_container_add(GTK_CONTAINER(values_exp),
                    dt_gui_vbox(d->edit_values, d->values_box));
  dt_gui_box_add(self->widget, values_exp);

  _refresh(self);
}

void gui_cleanup(dt_lib_module_t *self)
{
  dt_lib_lens_calib_fit_t *d = self->data;
  if(d && d->value_spins) g_ptr_array_free(d->value_spins, TRUE);
  g_free(self->data);
  self->data = NULL;
}

// clang-format off
// modelines: These editor modelines have been set for all relevant files by tools/update_modelines.py
// vim: shiftwidth=2 expandtab tabstop=2 cindent
// kate: tab-indents: off; indent-width 2; replace-tabs on; indent-mode cstyle; remove-trailing-spaces modified;
// clang-format on
