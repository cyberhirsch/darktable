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

/* Contributing a measured lens profile back to the lensfit project.
 *
 * A lens profile is only worth as much as the number of lenses it covers, and
 * no single person owns enough lenses to matter. Lensfun grew the way it did
 * because measurements flowed back. So sharing defaults to on.
 *
 * What that obliges us to get right:
 *
 *  - Nothing that identifies a person leaves the machine. Not the file name,
 *    not the path, not the camera's body serial, not the shot's location, and
 *    explicitly not a MAC address -- that would identify the photographer
 *    while failing at the only job it could have had, since a modern OS
 *    randomises it per interface. Installations are told apart by a random
 *    UUID generated locally, which is exactly as useful for deduplication and
 *    carries nothing about who generated it.
 *
 *  - The payload is shown before the first one is ever sent. Opt out is only
 *    honest if the person had a chance to see what they are opting out of.
 *
 *  - A submission that cannot be sent is queued, not dropped, and never
 *    blocks the save. Calibrating a lens is the user's work; uploading it is
 *    ours, and ours must not be able to fail theirs.
 */

#pragma once

#include "common/lens_warp.h"

#include <glib.h>

G_BEGIN_DECLS

/* How well the calibration went. The aggregator weights a measurement by
 * this, so a rough fit contributes without drowning out a careful one. */
typedef struct dt_lens_share_quality_t
{
  double straightness_before_px;  // line residual before the fit
  double straightness_after_px;   // and after; after > before means it failed
  int points_measured;            // chart nodes actually measured
  int points_stray;               // measured but on no lattice site
  double vig_coverage;            // fraction of the corner radius reached, 0 if none
  gboolean have_geometry;
  gboolean have_vignetting;
} dt_lens_share_quality_t;

/* Is the next saved profile going to be shared? This is the persisted state
 * of the checkbox beside "save profile" in the lensfit panel -- the only
 * place sharing is turned on or off. Defaults to TRUE -- see the note
 * above. */
gboolean dt_lens_share_enabled(void);

/* This installation's random identifier, generated on first use and kept in
 * the config. Never derived from hardware. */
const char *dt_lens_share_install_id(void);

/* Directory holding submissions waiting to be sent. */
gchar *dt_lens_share_queue_dir(void);

/* Wrap the profile at `profile_path` in a submission and queue it, then try
 * to send whatever is queued. Returns FALSE only if the submission could not
 * be written at all; a failed upload is not a failure, it stays queued.
 *
 * Never call this without the user having asked for it, either through the
 * preference or the checkbox beside the save button.
 */
gboolean dt_lens_share_submit(const char *profile_path,
                              const dt_lens_share_quality_t *quality,
                              GError **error);

/* Try to send everything queued. Safe to call when nothing is queued or when
 * no endpoint is configured, in which case it does nothing. */
void dt_lens_share_flush(void);

/* How many submissions are waiting. */
int dt_lens_share_queue_count(void);

G_END_DECLS

// clang-format off
// modelines: These editor modelines have been set for all relevant files by tools/update_modelines.py
// vim: shiftwidth=2 expandtab tabstop=2 cindent
// kate: tab-indents: off; indent-width 2; replace-tabs on; indent-mode cstyle; remove-trailing-spaces modified;
// clang-format on
