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

/* Writing the finished calibration out to something other than the lensfit
   vault: an STmap for a compositor, or a Lensfun database fragment.

   Saving the profile itself, and sharing it with the lensfit project, live
   in the lens panel (lens_calib_lens.c) instead -- they are part of naming
   and identifying the lens, not of exporting it, and putting "save" next
   to "maker"/"model" means there is one place that writes the vault rather
   than two. This panel only ever reads the profile name back out of
   plugins/lens_calib/profile_name to suggest a filename; it does not own
   it. */

#include "bauhaus/bauhaus.h"
#include "common/darktable.h"
#include "common/lens_lensfun.h"
#include "common/lens_warp.h"
#include "control/conf.h"
#include "control/control.h"
#include "gui/accelerators.h"
#include "gui/gtk.h"
#include "libs/lens_calib_common.h"
#include "libs/lib.h"
#include "libs/lib_api.h"
#include "views/view.h"

#include <glib/gstdio.h>
#include <string.h>

DT_MODULE(1)

typedef struct dt_lib_lens_calib_output_t
{
  GtkWidget *export_stmap;
  GtkWidget *bottom_up;
  GtkWidget *export_lensfun;
  gboolean updating;
} dt_lib_lens_calib_output_t;

const char *name(dt_lib_module_t *self)
{
  return _("export");
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
  return 100;
}

static void _refresh(dt_lib_module_t *self)
{
  dt_lib_lens_calib_output_t *d = self->data;
  dt_view_t *view = dt_lens_calib_get_view();

  const gboolean solved =
    view && darktable.view_manager->proxy.lens_calib.has_solution
    && darktable.view_manager->proxy.lens_calib.has_solution(view);

  gtk_widget_set_sensitive(d->export_stmap, solved);
  gtk_widget_set_sensitive(d->export_lensfun, solved);
}

// the name becomes a filename, so it cannot carry path separators
static gchar *_safe_name(const char *name)
{
  gchar *s = g_strdup(name ? name : "");
  g_strdelimit(s, "/\\:*?\"<>|", '_');
  return s;
}

// what the lens panel's profile-name entry currently holds, for a filename
// suggestion -- this panel does not own that value, only reads it
static gchar *_current_name_stem(void)
{
  gchar *name = dt_conf_get_string("plugins/lens_calib/profile_name");
  gchar *stem = _safe_name((name && *name) ? name : "lens");
  g_free(name);
  return stem;
}

static void _stmap_clicked(GtkWidget *w, dt_lib_module_t *self)
{
  dt_lib_lens_calib_output_t *d = self->data;
  dt_view_t *view = dt_lens_calib_get_view();
  if(!view || !darktable.view_manager->proxy.lens_calib.export_stmap) return;

  const gboolean bottom_up =
    gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(d->bottom_up));
  dt_conf_set_bool("plugins/lens_calib/stmap_bottom_up", bottom_up);

  GtkWidget *win = dt_ui_main_window(darktable.gui->ui);
  GtkWidget *chooser = gtk_file_chooser_dialog_new
    (_("export STmap"), GTK_WINDOW(win), GTK_FILE_CHOOSER_ACTION_SAVE,
     _("_cancel"), GTK_RESPONSE_CANCEL,
     _("_save"), GTK_RESPONSE_ACCEPT, NULL);

  gtk_file_chooser_set_do_overwrite_confirmation(GTK_FILE_CHOOSER(chooser),
                                                 TRUE);

  gchar *stem = _current_name_stem();
  gchar *base = g_strdup_printf("%s_stmap.exr", stem);
  gtk_file_chooser_set_current_name(GTK_FILE_CHOOSER(chooser), base);
  g_free(base);
  g_free(stem);

  gchar *dir = dt_lens_profile_dir();
  gtk_file_chooser_set_current_folder(GTK_FILE_CHOOSER(chooser), dir);
  g_free(dir);

  if(gtk_dialog_run(GTK_DIALOG(chooser)) == GTK_RESPONSE_ACCEPT)
  {
    gchar *path = gtk_file_chooser_get_filename(GTK_FILE_CHOOSER(chooser));
    darktable.view_manager->proxy.lens_calib.export_stmap(view, path,
                                                          bottom_up);
    g_free(path);
  }

  gtk_widget_destroy(chooser);
}

static void _lensfun_clicked(GtkWidget *w, dt_lib_module_t *self)
{
  dt_view_t *view = dt_lens_calib_get_view();
  if(!view || !darktable.view_manager->proxy.lens_calib.export_lensfun)
    return;

  GtkWidget *win = dt_ui_main_window(darktable.gui->ui);
  GtkWidget *chooser = gtk_file_chooser_dialog_new
    (_("export Lensfun XML"), GTK_WINDOW(win), GTK_FILE_CHOOSER_ACTION_SAVE,
     _("_cancel"), GTK_RESPONSE_CANCEL,
     _("_save"), GTK_RESPONSE_ACCEPT, NULL);

  gtk_file_chooser_set_do_overwrite_confirmation(GTK_FILE_CHOOSER(chooser),
                                                 TRUE);

  gchar *stem = _current_name_stem();
  gchar *base = g_strdup_printf("%s.xml", stem);
  gtk_file_chooser_set_current_name(GTK_FILE_CHOOSER(chooser), base);
  g_free(base);
  g_free(stem);

  /* Default to Lensfun's own per-user override directory rather than ours
     -- a file saved there is picked up by Lensfun itself on its next
     start, which is what "export" implies. Fall back to our own profile
     directory only if Lensfun's database has never been opened. */
  const char *lf_dir = dt_lf_home_data_dir();
  gchar *dir = (lf_dir && *lf_dir) ? g_strdup(lf_dir) : dt_lens_profile_dir();
  g_mkdir_with_parents(dir, 0755);
  gtk_file_chooser_set_current_folder(GTK_FILE_CHOOSER(chooser), dir);
  g_free(dir);

  if(gtk_dialog_run(GTK_DIALOG(chooser)) == GTK_RESPONSE_ACCEPT)
  {
    gchar *path = gtk_file_chooser_get_filename(GTK_FILE_CHOOSER(chooser));
    GPtrArray *bad = g_ptr_array_new_with_free_func(g_free);

    const gboolean ok = darktable.view_manager->proxy.lens_calib
      .export_lensfun(view, path, bad);

    if(ok)
    {
      gchar *written_name = g_path_get_basename(path);
      if(bad->len)
        dt_control_log(_("wrote `%s' -- %u measurement(s) could not be "
                         "converted, see the log"),
                       written_name, bad->len);
      else
        dt_control_log(_("wrote `%s'"), written_name);
      g_free(written_name);
    }
    else if(bad->len)
    {
      /* Say the actual reason rather than a dead end -- an anamorphic
         profile has no Lensfun equivalent at all, which is expected and
         not a bug, but "nothing converted" alone does not tell anyone
         that. */
      dt_control_log(_("nothing exported: %s"),
                     (const char *)g_ptr_array_index(bad, 0));
    }
    else
    {
      dt_control_log(_("nothing in this profile has a lensfun equivalent"));
    }

    for(guint i = 0; i < bad->len; i++)
      dt_print(DT_DEBUG_ALWAYS, "[lensfit] not exported: %s",
               (const char *)g_ptr_array_index(bad, i));

    g_ptr_array_unref(bad);
    g_free(path);
  }

  gtk_widget_destroy(chooser);
}

void gui_update(dt_lib_module_t *self)
{
  dt_lib_lens_calib_output_t *d = self->data;

  d->updating = TRUE;

  gtk_toggle_button_set_active
    (GTK_TOGGLE_BUTTON(d->bottom_up),
     dt_conf_get_bool("plugins/lens_calib/stmap_bottom_up"));

  d->updating = FALSE;

  _refresh(self);
}

void gui_init(dt_lib_module_t *self)
{
  dt_lib_lens_calib_output_t *d = g_malloc0(sizeof(*d));
  self->data = d;

  self->widget = dt_gui_vbox();

  d->bottom_up =
    gtk_check_button_new_with_label(_("STmap origin bottom left"));
  gtk_widget_set_tooltip_text
    (d->bottom_up,
     _("Nuke, Fusion and most compositors put the origin at the bottom\n"
       "left. turn this off for tools that use a top left origin --\n"
       "getting it wrong corrects the image upside down"));
  dt_gui_box_add(self->widget, d->bottom_up);

  d->export_stmap = dt_action_button_new
    (self, N_("export STmap"), _stmap_clicked, self,
     _("write the correction as a 32 bit OpenEXR STmap.\n"
       "this is the only interchange format that can carry an\n"
       "asymmetric warp exactly -- Lensfun's models cannot"), 0,
     (GdkModifierType)0);
  dt_gui_box_add(self->widget, d->export_stmap);

  d->export_lensfun = dt_action_button_new
    (self, N_("export Lensfun XML"), _lensfun_clicked, self,
     _("write the correction as a Lensfun database fragment.\n"
       "lossy: an anamorphic squeeze, an anisotropic vignetting, a gain\n"
       "convention vignetting curve, and any measurement with no focal\n"
       "length recorded have no lensfun equivalent and are left out --\n"
       "the log says what, if anything, could not be converted"), 0,
     (GdkModifierType)0);
  dt_gui_box_add(self->widget, d->export_lensfun);

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
