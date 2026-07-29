#!/usr/bin/env python3
"""Convert the Lensfun XML database into Lensfit JSON profiles.

Lensfun's models are all special cases of ours, so this direction is
lossless; the reverse is not, because nothing in Lensfun can express an
anamorphic squeeze or an anisotropic warp.

Three things this deliberately does *not* do:

  * It does not invent measurements. A Lensfun lens with no calibration at a
    given focal length produces no warp at that focal length, rather than a
    zero-filled one -- a coefficient set of all zeros claims a perfect lens,
    which is a stronger statement than saying nothing.

  * It does not silently reconcile conventions. Lensfun's vignetting model
    describes *transmission*, so its coefficients are negative and the
    correction divides; ours is the correction *gain* and multiplies. The
    conversion is explicit and recorded in the output.

  * It does not merge apertures. Lensfun indexes vignetting by focal length,
    aperture and focus distance together; each such combination becomes its
    own measurement, because collapsing them would average away the
    dependency that makes vignetting data worth having.

Usage:
    lensfun_to_lensfit.py <lensfun-xml-dir> <output-dir> [--limit N]
"""

import argparse
import glob
import json
import math
import os
import sys
import xml.etree.ElementTree as ET

# ---------------------------------------------------------------- conventions

# Recorded in every profile, because a convention that lives only in a comment
# eventually gets read the wrong way round. This file exists because that
# already happened once with the vignetting sign.
CONVENTIONS = {
    "coordinates": "normalized_half_diagonal",
    "origin": "frame_centre",
    "angles": "radians",
    "warp_direction": "observed_to_corrected",
    "map_direction": "destination_to_source",
    "vignetting": "correction_gain",
}

PROJECTION = {
    None: "rectilinear",
    "rectilinear": "rectilinear",
    "fisheye": "equidistant",
    "panoramic": "cylindrical",
    "equirectangular": "equirectangular",
    "orthographic": "orthographic",
    "stereographic": "stereographic",
    "equisolid": "equisolid",
    "fisheye_thoby": "thoby",
}


def _f(node, name, default=None):
    v = node.get(name)
    if v is None or v == "":
        return default
    try:
        return float(v)
    except ValueError:
        return default


def _text(lens, tag):
    """Lensfun repeats <name> with xml:lang; the untagged one is canonical."""
    best = None
    for e in lens.findall(tag):
        if e.get("{http://www.w3.org/XML/1998/namespace}lang") is None:
            return (e.text or "").strip()
        if best is None:
            best = (e.text or "").strip()
    return best


# ------------------------------------------------------------------ distortion


def _lensfun_rd(model, t, ru):
    """Lensfun's own formula, verbatim from lensfun.h, for verification."""
    if model == "poly3":
        k1 = t.get("k1", 0.0)
        return ru * (1.0 - k1 + k1 * ru * ru)
    if model == "poly5":
        return ru * (1.0 + t.get("k1", 0.0) * ru**2 + t.get("k2", 0.0) * ru**4)
    if model == "ptlens":
        a, b, c = t.get("a", 0.0), t.get("b", 0.0), t.get("c", 0.0)
        return ru * (a * ru**3 + b * ru**2 + c * ru + 1.0 - a - b - c)
    if model == "acm":
        return ru * (1.0 + t.get("k1", 0.0) * ru**2 + t.get("k2", 0.0) * ru**4
                     + t.get("k3", 0.0) * ru**6)
    return None


def _radial_poly_eval(p, r):
    """Our radial_poly: r_observed = c1 r + c2 r^2 + ... + c5 r^5."""
    return r * (p[1] + r * (p[2] + r * (p[3] + r * (p[4] + r * p[5]))))


def _distortion(node):
    """Lensfun distortion -> our radial_poly warp.

    All of Lensfun's models are r_d = r_u * P(r_u) with P a polynomial of
    degree <= 3, so all of them are exactly r_d = sum c_i r_u^i. Nothing is
    approximated here -- the coefficients are rearranged, not fitted.

    Note the direction: Lensfun defines corrected -> observed, which is why our
    radial_poly kind is defined that way too. Converting to our usual
    observed -> corrected direction would mean inverting a polynomial into a
    function that is not one.
    """
    model = node.get("model")
    focal = _f(node, "focal")
    if focal is None:
        return None

    # p[0] ellipticity (1 = circular, Lensfun has no anisotropy), p[1..5] = c1..c5
    p = [1.0, 0.0, 0.0, 0.0, 0.0, 0.0]
    terms = {}

    if model == "poly3":
        k1 = _f(node, "k1")
        if k1 is None:
            return None
        terms = {"k1": k1}
        p[1] = 1.0 - k1
        p[3] = k1
    elif model == "poly5":
        k1, k2 = _f(node, "k1", 0.0), _f(node, "k2", 0.0)
        terms = {"k1": k1, "k2": k2}
        p[1] = 1.0
        p[3] = k1
        p[5] = k2
    elif model == "ptlens":
        a, b, c = _f(node, "a", 0.0), _f(node, "b", 0.0), _f(node, "c", 0.0)
        terms = {"a": a, "b": b, "c": c}
        p[1] = 1.0 - a - b - c
        p[2] = c
        p[3] = b
        p[4] = a
    elif model == "acm":
        k1, k2, k3 = (_f(node, "k1", 0.0), _f(node, "k2", 0.0),
                      _f(node, "k3", 0.0))
        terms = {"k1": k1, "k2": k2, "k3": k3}
        p[1] = 1.0
        p[3] = k1
        p[5] = k2
        if abs(k3) > 1e-12:
            # r^7 exceeds our five terms; report rather than silently drop
            return {"_unrepresentable": f"acm k3={k3}"}
    else:
        return {"_unrepresentable": f"unknown distortion model {model!r}"}

    if abs(p[1] - 1.0) < 1e-12 and all(abs(x) < 1e-12 for x in p[2:]):
        return None  # identity says nothing

    # verify against Lensfun's formula over the full frame radius and beyond
    err = 0.0
    for i in range(1, 41):
        ru = 1.4 * i / 40.0
        ref = _lensfun_rd(model, terms, ru)
        err = max(err, abs(_radial_poly_eval(p, ru) - ref))

    return {"model": "radial_poly", "focal": focal, "params": p,
            "source_model": model, "_err": err}


# ------------------------------------------------------------------------- tca


def _tca(node):
    model = node.get("model")
    focal = _f(node, "focal")
    if focal is None:
        return None

    # our form is r * (t[0] r^2 + t[1] r + t[2]) per channel
    if model == "linear":
        kr = _f(node, "kr", 1.0)
        kb = _f(node, "kb", 1.0)
        if abs(kr - 1.0) < 1e-12 and abs(kb - 1.0) < 1e-12:
            return None
        return {"focal": focal, "source_model": "linear",
                "tca_r": [0.0, 0.0, kr], "tca_b": [0.0, 0.0, kb], "_err": 0.0}

    if model in ("poly3", "acm"):
        # Lensfun: r_d = r_u * (b r_u^2 + c r_u + v), which is exactly ours
        r = [_f(node, "br", 0.0), _f(node, "cr", 0.0), _f(node, "vr", 1.0)]
        b = [_f(node, "bb", 0.0), _f(node, "cb", 0.0), _f(node, "vb", 1.0)]
        if (abs(r[0]) < 1e-14 and abs(r[1]) < 1e-14 and abs(r[2] - 1) < 1e-14
                and abs(b[0]) < 1e-14 and abs(b[1]) < 1e-14
                and abs(b[2] - 1) < 1e-14):
            return None
        return {"focal": focal, "source_model": model,
                "tca_r": r, "tca_b": b, "_err": 0.0}

    return {"_unrepresentable": f"unknown tca model {model!r}"}


# ----------------------------------------------------------------- vignetting


def _vignetting(node):
    """Lensfun "pa": transmission t(r) = 1 + k1 r^2 + k2 r^4 + k3 r^6, and the
    correction divides by t. Ours is the gain g(r) = 1/t(r), and it multiplies.

    1/t is not a polynomial, so the coefficients cannot simply be negated.
    They are carried across in Lensfun's own convention and marked as such;
    the reader converts by division. Negating them would be the plausible
    wrong answer -- it agrees with 1/t to first order and diverges exactly
    where the correction is largest.
    """
    focal = _f(node, "focal")
    if focal is None:
        return None

    model = node.get("model") or "pa"
    if model not in ("pa", "acm"):
        return {"_unrepresentable": f"unknown vignetting model {model!r}"}

    k = [_f(node, "k1", 0.0), _f(node, "k2", 0.0), _f(node, "k3", 0.0)]
    if all(abs(x) < 1e-12 for x in k):
        return None

    return {
        "k": k,
        "ellipticity": 1.0,
        "focal": focal,
        "aperture": _f(node, "aperture"),
        "focus_distance": _f(node, "distance"),
        "convention": "transmission_divide",
        "radius_normalization": "calibration_sensor",
        "source_model": model,
        "_err": 0.0,
    }


# ---------------------------------------------------------------------- lenses


def convert_lens(lens, source_file):
    maker = _text(lens, "maker") or ""
    model = _text(lens, "model") or ""
    if not model:
        return None

    mounts = [(m.text or "").strip() for m in lens.findall("mount")]
    calib = lens.find("calibration")

    dist, tca, vig, bad = [], [], [], []
    worst = 0.0
    if calib is not None:
        for tag, fn, out in (("distortion", _distortion, dist),
                             ("tca", _tca, tca),
                             ("vignetting", _vignetting, vig)):
            for n in calib.findall(tag):
                e = fn(n)
                if not e:
                    continue
                if "_unrepresentable" in e:
                    bad.append(f"{tag}: {e['_unrepresentable']}")
                    continue
                worst = max(worst, e.pop("_err", 0.0))
                out.append(e)

    if not (dist or tca or vig):
        return None  # an identity entry is not a calibration

    ltype = lens.find("type")
    proj = PROJECTION.get(
        (ltype.text or "").strip() if ltype is not None else None, "rectilinear"
    )

    cx = _f(lens, "center-x", 0.0)
    cy = _f(lens, "center-y", 0.0)
    cxe = lens.find("center")
    if cxe is not None:
        cx = _f(cxe, "x", cx) or 0.0
        cy = _f(cxe, "y", cy) or 0.0

    return {
        "format": "lensfit",
        "version": 2,
        "conventions": CONVENTIONS,
        "imaging_system": {
            "manufacturer": maker,
            "model": model,
            "camera_type": "monocular",
            "channel_count": 1,
            "mounts": mounts,
            "crop_factor": _f(lens, "cropfactor", 1.0),
            "aspect_ratio": _f(lens, "aspect-ratio"),
        },
        "channels": [
            {
                "channel_id": 0,
                "lens_intrinsics": {
                    "projection_model": proj,
                    "cx": cx,
                    "cy": cy,
                    # Lensfun has no anisotropy at all; stating 1/1/0 is a
                    # fact about the source, not a measurement of the lens
                    "scale_x": 1.0,
                    "scale_y": 1.0,
                    "skew": 0.0,
                    "radial_normalization": "half_shorter_side",
                },
                "distortion": dist,
                "tca": tca,
                "photometric": {"vignetting": vig},
            }
        ],
        "provenance": {
            "source": "lensfun",
            "source_file": os.path.basename(source_file),
            "measured": False,
            "estimated_fields": ["scale_x", "scale_y", "skew"],
            "license": "CC-BY-SA-3.0 (Lensfun database)",
            "converter": "lensfun_to_lensfit.py",
            # largest deviation of our stored model from Lensfun's own formula,
            # sampled over the frame radius. Should be zero: the conversion is
            # a rearrangement, not a fit.
            "max_model_error": worst,
            "unrepresentable": bad,
        },
    }


def safe_name(maker, model):
    s = f"{maker}_{model}".strip("_")
    for ch in '/\\:*?"<>|':
        s = s.replace(ch, "_")
    return " ".join(s.split()).replace(" ", "_")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("xmldir")
    ap.add_argument("outdir")
    ap.add_argument("--limit", type=int, default=0)
    a = ap.parse_args()

    files = sorted(glob.glob(os.path.join(a.xmldir, "*.xml")))
    if not files:
        print(f"no xml in {a.xmldir}", file=sys.stderr)
        return 1

    os.makedirs(a.outdir, exist_ok=True)

    stats = dict(files=0, lenses=0, written=0, skipped=0,
                 dist=0, tca=0, vig=0, collide=0, bad=0)
    worst_overall = 0.0
    worst_who = ""
    bad_kinds = {}
    models = {}
    seen = {}

    for f in files:
        stats["files"] += 1
        try:
            root = ET.parse(f).getroot()
        except ET.ParseError as e:
            print(f"  ! {os.path.basename(f)}: {e}", file=sys.stderr)
            continue

        for lens in root.findall("lens"):
            stats["lenses"] += 1
            prof = convert_lens(lens, f)
            if prof is None:
                stats["skipped"] += 1
                continue

            ch = prof["channels"][0]
            stats["dist"] += len(ch["distortion"])
            stats["tca"] += len(ch["tca"])
            stats["vig"] += len(ch["photometric"]["vignetting"])

            for d in ch["distortion"]:
                models[d["source_model"]] = models.get(d["source_model"], 0) + 1
            for t in ch["tca"]:
                key = "tca:" + t["source_model"]
                models[key] = models.get(key, 0) + 1

            pv = prof["provenance"]
            if pv["max_model_error"] > worst_overall:
                worst_overall = pv["max_model_error"]
                worst_who = prof["imaging_system"]["model"]
            for b in pv["unrepresentable"]:
                stats["bad"] += 1
                bad_kinds[b.split("=")[0]] = bad_kinds.get(b.split("=")[0], 0) + 1

            base = safe_name(prof["imaging_system"]["manufacturer"],
                             prof["imaging_system"]["model"])
            name = base
            if base in seen:
                stats["collide"] += 1
                seen[base] += 1
                name = f"{base}__{seen[base]}"
            else:
                seen[base] = 0

            prof["imaging_system"]["id"] = name
            with open(os.path.join(a.outdir, name + ".json"), "w",
                      encoding="utf8") as fh:
                json.dump(prof, fh, indent=1, ensure_ascii=False)
            stats["written"] += 1

            if a.limit and stats["written"] >= a.limit:
                break
        if a.limit and stats["written"] >= a.limit:
            break

    print(f"xml files            {stats['files']}")
    print(f"lens entries         {stats['lenses']}")
    print(f"profiles written     {stats['written']}")
    print(f"skipped (no calib)   {stats['skipped']}")
    print(f"name collisions      {stats['collide']}")
    print(f"distortion measurements {stats['dist']}")
    print(f"tca measurements        {stats['tca']}")
    print(f"vignetting measurements {stats['vig']}")
    print("source models:")
    for k in sorted(models):
        print(f"   {k:16s} {models[k]}")
    print(f"unrepresentable         {stats['bad']}")
    for k in sorted(bad_kinds):
        print(f"   {k}: {bad_kinds[k]}")
    print(f"max model error         {worst_overall:.3e}"
          + (f"   ({worst_who})" if worst_who else ""))
    return 0


if __name__ == "__main__":
    sys.exit(main())
