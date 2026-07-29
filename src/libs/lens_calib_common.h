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

#pragma once

/* Shared scraps for the lens calibration panels.
 *
 * The panel is several sibling lib modules rather than one, so that each
 * part gets its own header in the right hand panel and can be collapsed
 * on its own -- a single module would be one long strip that is either
 * entirely open or entirely shut. What they have in common is only this:
 * how to reach the view, and how to put a label next to a widget.
 */

#include "gui/gtk.h"
#include "libs/lib.h"
#include "views/view.h"

static inline dt_view_t *dt_lens_calib_get_view(void)
{
  return darktable.view_manager->proxy.lens_calib.view;
}

/* Nudge the sibling panels to re-read the view.
 *
 * The five modules share one state and each only refreshes itself, so an
 * action in one leaves the others stale: measuring vignetting enables a
 * checkbox in `display`, detecting a grid enables overlays there, and neither
 * module knows the other did anything. A control that has quietly become
 * available but still looks disabled is indistinguishable from a bug.
 */
static inline void dt_lens_calib_refresh_panels(void)
{
  static const char *const mods[] = {
    "lens_calib_lens", "lens_calib_chart", "lens_calib_display",
    "lens_calib_fit", "lens_calib_output"
  };

  for(size_t i = 0; i < sizeof(mods) / sizeof(mods[0]); i++)
  {
    dt_lib_module_t *m = dt_lib_get_module(mods[i]);
    if(m) dt_lib_gui_queue_update(m);
  }
}

static inline GtkWidget *dt_lens_calib_labelled(const char *text,
                                                GtkWidget *w)
{
  GtkWidget *label = gtk_label_new(text);
  gtk_widget_set_halign(label, GTK_ALIGN_START);
  return dt_gui_hbox(dt_gui_expand(label), w);
}

static inline GtkWidget *dt_lens_calib_labelled_spin(const char *text,
                                                     GtkWidget **spin,
                                                     const int min,
                                                     const int max)
{
  *spin = gtk_spin_button_new_with_range(min, max, 1);
  gtk_spin_button_set_digits(GTK_SPIN_BUTTON(*spin), 0);
  return dt_lens_calib_labelled(text, *spin);
}

// clang-format off
// modelines: These editor modelines have been set for all relevant files by tools/update_modelines.py
// vim: shiftwidth=2 expandtab tabstop=2 cindent
// kate: tab-indents: off; indent-width 2; replace-tabs on; indent-mode cstyle; remove-trailing-spaces modified;
// clang-format on
