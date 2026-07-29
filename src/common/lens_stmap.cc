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
*/

/* C++ because OpenEXR's writer is C++; the header is C so the calibration
   view can call it without becoming C++. */

#include "common/lens_stmap.h"
#include "common/darktable.h"

#include <OpenEXR/ImfChannelList.h>
#include <OpenEXR/ImfFrameBuffer.h>
#include <OpenEXR/ImfHeader.h>
#include <OpenEXR/ImfOutputFile.h>
#include <OpenEXR/ImfStringAttribute.h>

#include <cmath>
#include <vector>

gboolean dt_lens_stmap_write(const dt_lens_warp_t *w,
                             const int width,
                             const int height,
                             const gboolean bottom_up,
                             const char *path,
                             GError **error)
{
  if(!w || !path || width < 2 || height < 2) return FALSE;

  const size_t npix = (size_t)width * height;

  std::vector<float> r, g, b;
  try
  {
    r.resize(npix);
    g.resize(npix);
    b.resize(npix, 0.0f);
  }
  catch(const std::exception &e)
  {
    g_set_error(error, G_FILE_ERROR, G_FILE_ERROR_NOMEM,
                "could not allocate a %dx%d map: %s", width, height, e.what());
    return FALSE;
  }

  const double hd = 0.5 * std::hypot((double)width, (double)height);
  const double halfw = 0.5 * width;
  const double halfh = 0.5 * height;

  int failed = 0;

#ifdef _OPENMP
#pragma omp parallel for reduction(+ : failed) schedule(static)
#endif
  for(int y = 0; y < height; y++)
  {
    for(int x = 0; x < width; x++)
    {
      /* Sample at pixel centres. Using the corner instead puts a half
         pixel shift into every correction, which is small enough to pass
         a visual check and large enough to matter when the map is used to
         match a plate. */
      const double ox = x + 0.5;
      const double oy = y + 0.5;

      const float u = (float)((ox - halfw) / hd);
      const float v = (float)((oy - halfh) / hd);

      float su, sv;
      // the map answers "where does this corrected pixel come from"
      if(!dt_lens_warp_invert(w, u, v, &su, &sv)) failed++;

      const double sx = su * hd + halfw;
      const double sy = sv * hd + halfh;

      size_t idx = (size_t)y * width + x;
      double tv = sy / height;
      if(bottom_up)
      {
        /* Two separate flips: the row we write to, and the coordinate we
           write. Doing only one of them produces a map that looks
           plausible and corrects upside down. */
        idx = (size_t)(height - 1 - y) * width + x;
        tv = 1.0 - tv;
      }

      r[idx] = (float)(sx / width);
      g[idx] = (float)tv;
    }
  }

  try
  {
    Imf::Header header(width, height);
    header.channels().insert("R", Imf::Channel(Imf::FLOAT));
    header.channels().insert("G", Imf::Channel(Imf::FLOAT));
    header.channels().insert("B", Imf::Channel(Imf::FLOAT));

    header.insert("darktable:generator",
                  Imf::StringAttribute("darktable lens calibration"));
    header.insert("darktable:origin",
                  Imf::StringAttribute(bottom_up ? "bottom-left" : "top-left"));

    Imf::OutputFile file(path, header);

    Imf::FrameBuffer fb;
    const size_t xs = sizeof(float);
    const size_t ys = sizeof(float) * width;
    fb.insert("R", Imf::Slice(Imf::FLOAT, (char *)r.data(), xs, ys));
    fb.insert("G", Imf::Slice(Imf::FLOAT, (char *)g.data(), xs, ys));
    fb.insert("B", Imf::Slice(Imf::FLOAT, (char *)b.data(), xs, ys));

    file.setFrameBuffer(fb);
    file.writePixels(height);
  }
  catch(const std::exception &e)
  {
    g_set_error(error, G_FILE_ERROR, G_FILE_ERROR_FAILED,
                "could not write `%s': %s", path, e.what());
    return FALSE;
  }

  if(failed)
    dt_print(DT_DEBUG_ALWAYS,
             "[lens_stmap] %d of %zu pixels did not converge; the map is"
             " approximate there, which happens when the model is"
             " extrapolated well past the measured area",
             failed, npix);

  dt_print(DT_DEBUG_ALWAYS, "[lens_stmap] wrote %dx%d map to `%s'",
           width, height, path);

  return TRUE;
}

// clang-format off
// modelines: These editor modelines have been set for all relevant files by tools/update_modelines.py
// vim: shiftwidth=2 expandtab tabstop=2 cindent
// kate: tab-indents: off; indent-width 2; replace-tabs on; indent-mode cstyle; remove-trailing-spaces modified;
// clang-format on
