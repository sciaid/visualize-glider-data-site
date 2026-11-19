#!/usr/bin/env python3
"""
fetch_currents.py - Ocean currents around each glider that is out now, for the real-time map.

    python scripts/fetch_currents.py .erddap-cache currents

Reads the real-time grids written by fetch_erddap.py (<cache>/ERDDAP/*_erddap-realtime-5m.nc)
and writes <out>/currents.json with, for each glider:

  * its track (decimated) and latest position, and its typical dive depth;
  * model forecast (Copernicus Marine global analysis and forecast, 1/12 deg, daily means,
    yesterday to +5 days): currents on every model level from the surface to the glider's dive
    depth, on a grid around the glider and as a column at its latest position, so the page can
    show the mean over any depth from the surface to the dive depth;
  * observed surface current (HF radar, US West Coast 6 km, mean of the last 24 h) where the
    glider is inside the HF radar domain.

Copernicus needs a free account: set COPERNICUSMARINE_SERVICE_USERNAME and
COPERNICUSMARINE_SERVICE_PASSWORD (GitHub repository secrets in the site build). Without them,
or if a download fails, that layer is left out and the rest is still written. The script never
fails the build.

Requirements: pip install copernicusmarine netCDF4 numpy
"""
from __future__ import annotations

import argparse
import datetime as dt
import json
import os
import sys
import tempfile
import traceback
import urllib.request
from pathlib import Path

import numpy as np
import netCDF4

MODEL_DATASET = "cmems_mod_glo_phy-cur_anfc_0.083deg_P1D-m"      # daily mean uo, vo on depth levels
MODEL_NAME = "Copernicus Marine global analysis and forecast (1/12°, daily means)"
MODEL_DOI = "https://doi.org/10.48670/moi-00016"
HFR_ERDDAP = "https://coastwatch.pfeg.noaa.gov/erddap/griddap"
# US West Coast, hourly, 6 km. CoastWatch serves it under two IDs (longitudes -180..180 and 0..360)
# and either can be briefly unavailable, so both are tried, with retries.
HFR_DATASETS = [("ucsdHfrW6", 0.0), ("ucsdHfrW6_Lon0360", 360.0)]
USER_AGENT = "Mozilla/5.0 (X11; Linux x86_64) visualize-glider-data-site (+https://github.com/sciaid/visualize-glider-data-site)"
HFR_BOX = (30.25, 49.99, -130.36, -115.80)                       # lat0, lat1, lon0, lon1
HALF_BOX_DEG = 1.25        # map box: +/- this many degrees of latitude around the glider
LAYER_DEPTH = 100.0        # the 0-100 m mean shown on the map
ACTIVE_DAYS = 14           # gliders whose last data are older than this are skipped
DAYS_AHEAD = 5


def _f(v):
    a = v[:]
    return np.ma.filled(a.astype("f8"), np.nan) if np.ma.isMaskedArray(a) else np.asarray(a, "f8")


# --------------------------------------------------------------------------- #
# Gliders
# --------------------------------------------------------------------------- #
def gliders(cache):
    """Active real-time gliders from fetch_erddap.py grids."""
    out = []
    now = dt.datetime.now(dt.timezone.utc).timestamp()
    for p in sorted(Path(cache, "ERDDAP").glob("*_erddap-realtime-5m.nc")):
        with netCDF4.Dataset(p) as ds:
            t, la, lo = _f(ds["time"]), _f(ds["latitude"]), _f(ds["longitude"])
            ok = np.isfinite(t) & np.isfinite(la) & np.isfinite(lo)
            if not ok.any():
                continue
            t, la, lo = t[ok], la[ok], lo[ok]
            o = np.argsort(t); t, la, lo = t[o], la[o], lo[o]
            if now - t[-1] > ACTIVE_DAYS * 86400:
                print(f"  {p.name}: last data {(now - t[-1]) / 86400:.0f} days ago, skipped")
                continue
            # typical dive depth: 90th percentile of the deepest sample per profile, last 30 profiles
            if "sample_depth" in ds.variables:
                sp, sd = np.asarray(ds["sample_profile"][:], int), _f(ds["sample_depth"])
                prof = [np.nanmax(sd[sp == k]) for k in np.unique(sp)[-30:] if np.isfinite(sd[sp == k]).any()]
            else:
                d = _f(ds["depth"]); T = _f(ds["temperature"]) if "temperature" in ds.variables else None
                prof = [] if T is None else [d[np.isfinite(T[:, k])].max() for k in range(T.shape[1])[-30:]
                                             if np.isfinite(T[:, k]).any()]
            dive = float(np.round(np.percentile(prof, 90) / 10) * 10) if prof else None
            step = max(1, t.size // 1500)
            out.append({"name": str(getattr(ds, "deployment_name", p.stem.split("_")[0])),
                        "lat": round(float(la[-1]), 5), "lon": round(float(lo[-1]), 5), "time": int(t[-1]),
                        "diveDepth": dive,
                        "track": {"t": [int(x) for x in t[::step]], "lat": [round(float(x), 4) for x in la[::step]],
                                  "lon": [round(float(x), 4) for x in lo[::step]]}})
    return out


def box(g):
    dlon = HALF_BOX_DEG / max(0.2, np.cos(np.radians(g["lat"])))
    return g["lat"] - HALF_BOX_DEG, g["lat"] + HALF_BOX_DEG, g["lon"] - dlon, g["lon"] + dlon


# --------------------------------------------------------------------------- #
# Depth averages
# --------------------------------------------------------------------------- #
def layer_weights(depth, bottom):
    """Thickness of each model level inside [0, bottom] (levels are cell centres)."""
    z = np.asarray(depth, float)
    edges = np.concatenate([[0.0], (z[:-1] + z[1:]) / 2, [z[-1] + (z[-1] - z[-2] if z.size > 1 else 1)]])
    top, bot = np.clip(edges[:-1], 0, bottom), np.clip(edges[1:], 0, bottom)
    return bot - top


def depth_mean(u, depth, bottom):
    """Thickness-weighted mean over 0..bottom of u[..., depth, y, x] (NaN below the seafloor is skipped).
    Returns NaN where less than half of 0..bottom has data."""
    w = layer_weights(depth, bottom).reshape((1,) * (u.ndim - 3) + (-1, 1, 1))
    ok = np.isfinite(u)
    num = np.nansum(np.where(ok, u, 0) * w, axis=-3)
    den = np.sum(np.where(ok, w, 0), axis=-3)
    return np.where(den >= 0.5 * min(bottom, w.sum()), num / np.maximum(den, 1e-9), np.nan)


def grid_block(lat, lon, u, v, step=2):
    """Compact grid for the page: every `step`-th cell, values in cm/s rounded."""
    la, lo = lat[::step], lon[::step]
    uu, vv = u[::step, ::step], v[::step, ::step]
    enc = lambda a: [None if not np.isfinite(x) else int(round(x * 100)) for x in a.ravel()]
    return {"lat0": round(float(la[0]), 4), "dlat": round(float(la[1] - la[0]), 5) if la.size > 1 else 0,
            "lon0": round(float(lo[0]), 4), "dlon": round(float(lo[1] - lo[0]), 5) if lo.size > 1 else 0,
            "ny": int(la.size), "nx": int(lo.size), "u": enc(uu), "v": enc(vv)}


def nearest(lat, lon, la, lo):
    return int(np.abs(lat - la).argmin()), int(np.abs(lon - lo).argmin())


# --------------------------------------------------------------------------- #
# Copernicus Marine model
# --------------------------------------------------------------------------- #
def model_currents(g, workdir):
    """Daily-mean model currents on every model level from the surface to the glider's dive depth,
    on a grid around the glider (every 2nd cell) plus the full-resolution column at the glider.
    The page averages levels down to whatever depth its slider is set to."""
    import copernicusmarine
    lat0, lat1, lon0, lon1 = box(g)
    today = dt.datetime.now(dt.timezone.utc).date()
    start, end = today - dt.timedelta(days=1), today + dt.timedelta(days=DAYS_AHEAD)
    bottom = max(LAYER_DEPTH, g["diveDepth"] or LAYER_DEPTH)
    fn = f"{g['name']}_model.nc"
    copernicusmarine.subset(
        dataset_id=MODEL_DATASET, variables=["uo", "vo"],
        minimum_longitude=lon0, maximum_longitude=lon1, minimum_latitude=lat0, maximum_latitude=lat1,
        minimum_depth=0, maximum_depth=bottom * 1.25 + 20,
        start_datetime=f"{start}T00:00:00", end_datetime=f"{end}T23:59:59",
        output_directory=str(workdir), output_filename=fn, disable_progress_bar=True)
    with netCDF4.Dataset(Path(workdir) / fn) as ds:
        tv = ds["time"]
        times = netCDF4.num2date(tv[:], tv.units, only_use_cftime_datetimes=False, only_use_python_datetimes=True)
        lat, lon, depth = _f(ds["latitude"]), _f(ds["longitude"]), _f(ds["depth"])
        u, v = _f(ds["uo"]), _f(ds["vo"])                        # time, depth, lat, lon
    return model_block(g, times, lat, lon, depth, u, v, bottom, str(today))


def model_block(g, times, lat, lon, depth, u, v, bottom, today):
    # keep levels whose layer starts above the deepest depth the slider can reach
    edges = np.concatenate([[0.0], (depth[:-1] + depth[1:]) / 2, [depth[-1] + 1]])
    nl = max(1, int(np.searchsorted(edges, bottom, side="left")))
    depth, edges = depth[:nl], edges[:nl + 1]
    u, v = u[:, :nl], v[:, :nl]
    j, i = nearest(lat, lon, g["lat"], g["lon"])
    step = 2
    meta = grid_block(lat, lon, u[0, 0], v[0, 0], step=step)
    cm = lambda a: [None if not np.isfinite(x) else int(round(x * 100)) for x in np.asarray(a).ravel()]
    days = []
    for n, tt in enumerate(times):
        days.append({"date": tt.strftime("%Y-%m-%d"),
                     # level-major: u[level][cell], cells row-major (lat, lon), every 2nd cell
                     "u": [cm(u[n, k, ::step, ::step]) for k in range(nl)],
                     "v": [cm(v[n, k, ::step, ::step]) for k in range(nl)],
                     # full-resolution column at the glider, cm/s with one decimal
                     "colU": [None if not np.isfinite(x) else round(float(x) * 100, 1) for x in u[n, :, j, i]],
                     "colV": [None if not np.isfinite(x) else round(float(x) * 100, 1) for x in v[n, :, j, i]]})
    return {"source": MODEL_NAME, "dataset": MODEL_DATASET, "doi": MODEL_DOI, "today": today,
            "levels": [round(float(z), 2) for z in depth], "edges": [round(float(z), 2) for z in edges],
            "maxDepth": float(bottom), "defaultDepth": LAYER_DEPTH,
            "grid": {k: meta[k] for k in ("lat0", "dlat", "lon0", "dlon", "ny", "nx")}, "days": days}


# --------------------------------------------------------------------------- #
# HF radar (observed)
# --------------------------------------------------------------------------- #
def in_hfr(g):
    return HFR_BOX[0] <= g["lat"] <= HFR_BOX[1] and HFR_BOX[2] <= g["lon"] <= HFR_BOX[3]


def hfr_currents(g, workdir):
    import time
    lat0, lat1, lon0, lon1 = box(g)
    lat0, lat1 = max(lat0, HFR_BOX[0]), min(lat1, HFR_BOX[1])
    lon0, lon1 = max(lon0, HFR_BOX[2]), min(lon1, HFR_BOX[3])
    path = Path(workdir) / f"{g['name']}_hfr.nc"
    errors = []
    for attempt in range(3):
        for did, shift in HFR_DATASETS:
            sel = f"[last-23:1:last][({lat0:.3f}):1:({lat1:.3f})][({lon0 + shift:.3f}):1:({lon1 + shift:.3f})]"
            url = f"{HFR_ERDDAP}/{did}.nc?water_u{sel},water_v{sel}"
            try:
                req = urllib.request.Request(url, headers={"User-Agent": USER_AGENT})
                with urllib.request.urlopen(req, timeout=300) as r:
                    path.write_bytes(r.read())
                break
            except Exception as exc:
                errors.append(f"{did}: {exc}")
        else:
            time.sleep(20 * (attempt + 1))
            continue
        break
    else:
        raise RuntimeError("; ".join(errors[-2:]))
    with netCDF4.Dataset(path) as ds:
        tv = ds["time"]; t = _f(tv)
        lat, lon = _f(ds["latitude"]), _f(ds["longitude"])
        lon = np.where(lon > 180, lon - 360, lon)
        u, v = _f(ds["water_u"]), _f(ds["water_v"])             # time, lat, lon
    if lat[0] > lat[-1]:
        lat, u, v = lat[::-1], u[:, ::-1], v[:, ::-1]
    n = np.sum(np.isfinite(u), axis=0)
    mu = np.where(n >= 12, np.nanmean(u, axis=0), np.nan)       # at least 12 of 24 hours
    mv = np.where(n >= 12, np.nanmean(v, axis=0), np.nan)
    iso = lambda s: dt.datetime.fromtimestamp(float(s), dt.timezone.utc).strftime("%Y-%m-%d %H:%M UTC")
    return {"source": "HF radar, US West Coast 6 km (CORDC/Scripps via NOAA CoastWatch ERDDAP)",
            "dataset": "ucsdHfrW6", "start": iso(np.nanmin(t)), "end": iso(np.nanmax(t)),
            "grid": grid_block(lat, lon, mu, mv, step=3)}   # ~18 km, like the model arrows


# --------------------------------------------------------------------------- #
def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("cache", help="fetch_erddap.py output folder (holds ERDDAP/*_erddap-realtime-5m.nc)")
    ap.add_argument("out", help="folder for currents.json")
    a = ap.parse_args(argv)
    out = Path(a.out); out.mkdir(parents=True, exist_ok=True)
    made = dt.datetime.now(dt.timezone.utc).strftime("%Y-%m-%d %H:%M UTC")
    result = {"made": made, "gliders": [], "notes": []}
    gs = gliders(a.cache)
    print(f"{len(gs)} active real-time gliders")
    have_cm = bool(os.environ.get("COPERNICUSMARINE_SERVICE_USERNAME") and os.environ.get("COPERNICUSMARINE_SERVICE_PASSWORD"))
    if gs and not have_cm:
        result["notes"].append("Model currents not available: Copernicus Marine credentials are not set.")
        print("  Copernicus Marine credentials not set: model layer skipped")
    with tempfile.TemporaryDirectory() as tmp:
        for g in gs:
            print(f"  {g['name']}: {g['lat']:.3f}, {g['lon']:.3f}, dive depth {g['diveDepth']} m")
            if have_cm:
                try:
                    g["model"] = model_currents(g, tmp)
                    m = g["model"]
                    print(f"    model: {len(m['days'])} days, {len(m['levels'])} levels")
                    for d in m["days"]:                  # at the glider: surface, 0-100 m, 0-dive
                        u, v = np.array(d["colU"], float), np.array(d["colV"], float)
                        e = np.array(m["edges"])
                        def mean(D):
                            w = np.clip(np.minimum(e[1:], D) - np.minimum(e[:-1], D), 0, None) if D > 0 else np.eye(1, u.size)[0]
                            ok = np.isfinite(u) & np.isfinite(v) & (w > 0)
                            return (np.sum(w[ok] * u[ok]) / w[ok].sum(), np.sum(w[ok] * v[ok]) / w[ok].sum()) if ok.any() else (np.nan, np.nan)
                        txt = "  ".join(f"{lab} {np.hypot(*mean(D)):4.1f} cm/s @{(np.degrees(np.arctan2(*mean(D))) % 360):3.0f}"
                                        for lab, D in (("sfc", 0), ("0-100", 100), (f"0-{m['maxDepth']:.0f}", m["maxDepth"])))
                        print(f"      {d['date']}: {txt}")
                except Exception as exc:
                    traceback.print_exc()
                    result["notes"].append(f"{g['name']}: model currents failed ({type(exc).__name__}).")
            if in_hfr(g):
                try:
                    g["hfr"] = hfr_currents(g, tmp)
                    print(f"    HF radar: {g['hfr']['start']} to {g['hfr']['end']}")
                except Exception as exc:
                    traceback.print_exc()
                    result["notes"].append(f"{g['name']}: HF radar failed ({type(exc).__name__}).")
            result["gliders"].append(g)
    (out / "currents.json").write_text(json.dumps(result, separators=(",", ":")), encoding="utf-8")
    print(f"Wrote {out / 'currents.json'} ({(out / 'currents.json').stat().st_size / 1e3:.0f} kB)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
