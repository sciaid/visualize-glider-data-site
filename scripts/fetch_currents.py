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
import multiprocessing as mp
from concurrent.futures import ProcessPoolExecutor
from concurrent.futures.process import BrokenProcessPool
import urllib.request
from pathlib import Path

import numpy as np
import netCDF4

import fetch_satellite as sat

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
DAC_DAYS = 7          # glider-measured currents shown for the last week


def glider_dac(ds, now):
    """Glider-measured depth-averaged currents (cm/s) for the last DAC_DAYS, one per dive segment,
    plus a status: "ok", "missing" (no values) or "suspect" (values that cannot be real: u equal to v
    for nearly every segment, or every value below 1 cm/s). Suspect values are not shown on the map.
    Positions are lat_uv/lon_uv where filled, otherwise the profile position."""
    if "dac_u" not in ds.variables or "dac_v" not in ds.variables:
        return {"status": "missing", "reason": "u and v are not in the real-time dataset"}
    t, la, lo = _f(ds["time"]), _f(ds["latitude"]), _f(ds["longitude"])
    u, v = _f(ds["dac_u"]) * 100, _f(ds["dac_v"]) * 100
    get = lambda n, fb: np.where(np.isfinite(_f(ds[n])), _f(ds[n]), fb) if n in ds.variables else fb
    tt, lat, lon = get("dac_time_uv", t), get("dac_lat_uv", la), get("dac_lon_uv", lo)
    pos_from = "lat_uv/lon_uv" if "dac_lat_uv" in ds.variables and np.isfinite(_f(ds["dac_lat_uv"])).any() \
        else "profile position"
    ok = np.isfinite(u) & np.isfinite(v) & np.isfinite(tt) & np.isfinite(lat) & np.isfinite(lon) \
        & (np.hypot(u, v) < 150)
    if not ok.any():
        return {"status": "missing", "reason": "no depth-averaged current values in the real-time dataset"}
    o = np.argsort(tt[ok], kind="stable")
    tt, lat, lon, u, v = (a[ok][o] for a in (tt, lat, lon, u, v))
    # one value per dive segment: consecutive profiles often repeat the same segment value
    new = np.r_[True, (np.diff(u) != 0) | (np.diff(v) != 0)]
    tt, lat, lon, u, v = (a[new] for a in (tt, lat, lon, u, v))
    same = float(np.mean(u == v))
    speed = np.hypot(u, v)
    out = {"positions": pos_from, "n": int(u.size)}
    if same > 0.9:
        out.update(status="suspect", reason=f"u equals v in {same:.0%} of segments")
    elif np.percentile(speed, 95) < 1:
        out.update(status="suspect", reason="every value is below 1 cm/s")
    else:
        out["status"] = "ok"
    recent = tt >= now - DAC_DAYS * 86400
    if out["status"] == "ok" and recent.any():
        r = lambda a, k: [round(float(x), k) for x in a[recent]]
        out.update(t=[int(x) for x in tt[recent]], lat=r(lat, 4), lon=r(lon, 4), u=r(u, 1), v=r(v, 1))
    elif out["status"] == "ok":
        out.update(status="missing", reason=f"no values in the last {DAC_DAYS} days")
    return out


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
                                  "lon": [round(float(x), 4) for x in lo[::step]]},
                        "dac": glider_dac(ds, now),
                        "surface": sat.glider_surface(ds, now)})
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
# NOAA WCOFS (West Coast Operational Forecast System)
# --------------------------------------------------------------------------- #
WCOFS_TDS = "https://opendap.co-ops.nos.noaa.gov/thredds/dodsC/NOAA/WCOFS/MODELS"
WCOFS_NAME = "NOAA West Coast Operational Forecast System (WCOFS, 4 km, daily means of 3-hourly output)"
WCOFS_URL = "https://tidesandcurrents.noaa.gov/ofs/wcofs/wcofs.html"
WCOFS_BOX = (18.42, 55.78, -145.38, -112.10)
WCOFS_STRIDE = 2           # every 2nd 0.04 deg cell (~8 km); the page shows every 2nd of those


def in_wcofs(g):
    return WCOFS_BOX[0] + 1 <= g["lat"] <= WCOFS_BOX[1] - 1 and WCOFS_BOX[2] + 1 <= g["lon"] <= WCOFS_BOX[3] - 0.2


def _wcofs_url(run, kind, hour):
    d = run.strftime("%Y/%m/%d")
    return f"{WCOFS_TDS}/{d}/wcofs.t03z.{run:%Y%m%d}.regulargrid.{kind}{hour:03d}.nc"


def wcofs_latest_run(now=None):
    """The newest 03Z run whose last forecast file (f072) is on the server: today's or yesterday's."""
    now = now or dt.datetime.now(dt.timezone.utc)
    for back in range(0, 3):
        run = dt.datetime(now.year, now.month, now.day, 3, tzinfo=dt.timezone.utc) - dt.timedelta(days=back)
        if run > now:
            continue
        try:
            with netCDF4.Dataset(_wcofs_url(run, "f", 72)):
                return run
        except OSError:
            continue
    raise RuntimeError("no complete WCOFS run on the server for the last 2 days")


def wcofs_currents(g):
    """Daily-mean WCOFS currents around the glider, yesterday to +3 days, on every z-level down to the
    glider's dive depth. Output in the same layout as the Copernicus block."""
    from concurrent.futures import ThreadPoolExecutor
    run = wcofs_latest_run()
    lat0, lat1, lon0, lon1 = box(g)
    bottom = max(LAYER_DEPTH, g["diveDepth"] or LAYER_DEPTH)
    with netCDF4.Dataset(_wcofs_url(run, "f", 3)) as ds:
        depth = _f(ds["Depth"])
        lat = _f(ds["Latitude"][:, 0]) if ds["Latitude"].ndim == 2 else _f(ds["Latitude"])
        lon = _f(ds["Longitude"][0, :]) if ds["Longitude"].ndim == 2 else _f(ds["Longitude"])
        t_units = ds["time"].units
    nz = int(np.searchsorted(depth, bottom, side="left")) + 1          # levels down to (and including) bottom
    nz = min(nz, depth.size)
    jj = np.where((lat >= lat0) & (lat <= lat1))[0]; ii = np.where((lon >= lon0) & (lon <= lon1))[0]
    j0, j1, i0, i1 = jj.min(), jj.max() + 1, ii.min(), ii.max() + 1
    sl = (0, slice(0, nz), slice(j0, j1, WCOFS_STRIDE), slice(i0, i1, WCOFS_STRIDE))
    steps = [("n", h) for h in range(3, 25, 3)] + [("f", h) for h in range(3, 73, 3)]

    def get(step):
        kind, h = step
        for attempt in range(3):
            try:
                with netCDF4.Dataset(_wcofs_url(run, kind, h)) as ds:
                    t = netCDF4.num2date(ds["time"][:], t_units, only_use_cftime_datetimes=False,
                                         only_use_python_datetimes=True)[0]
                    u = np.ma.filled(ds["u_eastward"][sl].astype("f8"), np.nan)
                    v = np.ma.filled(ds["v_northward"][sl].astype("f8"), np.nan)
                return t.replace(tzinfo=dt.timezone.utc), u, v
            except OSError:
                if attempt == 2:
                    return None
    with ThreadPoolExecutor(max_workers=4) as ex:
        got = [x for x in ex.map(get, steps) if x is not None]
    if not got:
        raise RuntimeError("no WCOFS files could be read")
    by_day = {}
    for t, u, v in got:
        by_day.setdefault(t.date(), []).append((u, v))
    days = sorted(d for d, xs in by_day.items() if len(xs) >= 6)      # at least 6 of 8 three-hourly steps
    times = [dt.datetime(d.year, d.month, d.day) for d in days]
    U = np.stack([np.nanmean([x[0] for x in by_day[d]], axis=0) for d in days])
    V = np.stack([np.nanmean([x[1] for x in by_day[d]], axis=0) for d in days])
    out = model_block(g, times, lat[j0:j1:WCOFS_STRIDE], lon[i0:i1:WCOFS_STRIDE], depth[:nz], U, V, bottom,
                      str(dt.datetime.now(dt.timezone.utc).date()))
    out.update({"source": WCOFS_NAME, "dataset": f"WCOFS {run:%Y-%m-%d} 03Z run", "doi": WCOFS_URL,
                "run": f"{run:%Y-%m-%d %H:%M} UTC", "steps": len(got)})
    return out


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
def log_at_glider(m):
    """Print the model current at the glider each day: surface, 0-100 m and 0-dive depth (for checking)."""
    e = np.array(m["edges"])
    for d in m["days"]:
        u, v = np.array(d["colU"], float), np.array(d["colV"], float)

        def mean(D):
            w = np.clip(np.minimum(e[1:], D) - np.minimum(e[:-1], D), 0, None) if D > 0 else np.eye(1, u.size)[0]
            ok = np.isfinite(u) & np.isfinite(v) & (w > 0)
            return (np.sum(w[ok] * u[ok]) / w[ok].sum(), np.sum(w[ok] * v[ok]) / w[ok].sum()) if ok.any() else (np.nan, np.nan)
        txt = "  ".join(f"{lab} {np.hypot(*mean(D)):4.1f} cm/s @{(np.degrees(np.arctan2(*mean(D))) % 360):3.0f}"
                        for lab, D in (("sfc", 0), ("0-100", 100), (f"0-{m['maxDepth']:.0f}", m["maxDepth"])))
        print(f"      {d['date']}: {txt}")


def _call(fn, *args):
    """Run fn in the child; send back exceptions as plain RuntimeErrors (some, like HTTPError,
    hold open files and cannot be passed back to the parent process)."""
    try:
        return fn(*args)
    except Exception as exc:
        traceback.print_exc()
        raise RuntimeError(f"{type(exc).__name__}: {exc}") from None


def isolated(fn, *args, tries=2):
    """fn(*args) in a separate process. The Copernicus and WCOFS readers occasionally crash the
    Python process (segmentation fault inside the netCDF/HDF5 libraries); in a child process
    that loses only this one source, and it is tried once more."""
    for k in range(tries):
        with ProcessPoolExecutor(max_workers=1, mp_context=mp.get_context("spawn")) as ex:
            try:
                return ex.submit(_call, fn, *args).result()
            except BrokenProcessPool:
                print(f"    {fn.__name__} crashed (attempt {k + 1} of {tries})", flush=True)
    raise RuntimeError(f"{fn.__name__} crashed {tries} times")


def write_result(out, result):
    (out / "currents.json").write_text(json.dumps(result, separators=(",", ":")), encoding="utf-8")


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
    # the gliders, their tracks and their own measured currents go out first, so the map still
    # lists them if every outside source fails
    write_result(out, dict(result, gliders=gs))
    have_cm = bool(os.environ.get("COPERNICUSMARINE_SERVICE_USERNAME") and os.environ.get("COPERNICUSMARINE_SERVICE_PASSWORD"))
    if gs and not have_cm:
        result["notes"].append("Model currents not available: Copernicus Marine credentials are not set.")
        print("  Copernicus Marine credentials not set: model layer skipped")
    with tempfile.TemporaryDirectory() as tmp:
        for g in gs:
            print(f"  {g['name']}: {g['lat']:.3f}, {g['lon']:.3f}, dive depth {g['diveDepth']} m")
            d = g["dac"]
            print(f"    glider-measured currents: {d['status']}"
                  + (f" ({d['reason']})" if d.get("reason") else "")
                  + (f", {len(d['u'])} segments in the last {DAC_DAYS} days, latest "
                     f"{np.hypot(d['u'][-1], d['v'][-1]):.0f} cm/s" if d.get("u") else ""))
            if have_cm:
                try:
                    g["model"] = isolated(model_currents, g, tmp)
                    m = g["model"]
                    print(f"    model: {len(m['days'])} days, {len(m['levels'])} levels")
                    log_at_glider(m)
                except Exception as exc:
                    traceback.print_exc()
                    result["notes"].append(f"{g['name']}: model currents failed ({type(exc).__name__}).")
            if in_wcofs(g):
                try:
                    g["wcofs"] = isolated(wcofs_currents, g)
                    m = g["wcofs"]
                    print(f"    WCOFS: run {m['run']}, {m['steps']} files, {len(m['days'])} days, {len(m['levels'])} levels")
                    log_at_glider(m)
                except Exception as exc:
                    traceback.print_exc()
                    result["notes"].append(f"{g['name']}: WCOFS currents failed ({type(exc).__name__}).")
            if in_hfr(g):
                try:
                    g["hfr"] = isolated(hfr_currents, g, tmp)
                    print(f"    HF radar: {g['hfr']['start']} to {g['hfr']['end']}")
                except Exception as exc:
                    traceback.print_exc()
                    result["notes"].append(f"{g['name']}: HF radar failed ({type(exc).__name__}).")
            for var, fn, need_cm in (("sst", sat.sst_layer, False), ("chl", sat.chl_layer, True)):
                if need_cm and not have_cm:
                    continue
                try:
                    g.setdefault("sat", {})[var] = b = isolated(fn, g, box(g), out, tmp)
                    at = [d for d in b["days"] if "atGlider" in d]
                    last = at[-1] if at else None
                    print(f"    satellite {var}: {len(b['days'])} days, colour range {b['range']}"
                          + (f", at glider {last['date']}: {last['atGlider']['v']}" + (
                              f" (observed)" if last['atGlider'].get('observed') else
                              " (interpolated)" if last['atGlider'].get('observed') is False else "")
                             + f", glider {g['surface'].get(last['date'], {}).get('T' if var == 'sst' else 'chl')}"
                             if last else "") + (f"; {b['note']}" if b.get("note") else ""))
                except Exception as exc:
                    traceback.print_exc()
                    result["notes"].append(f"{g['name']}: satellite {'SST' if var == 'sst' else 'chlorophyll'} failed ({type(exc).__name__}).")
            result["gliders"].append(g)
            write_result(out, dict(result, gliders=result["gliders"] + gs[len(result["gliders"]):]))
    write_result(out, result)
    print(f"Wrote {out / 'currents.json'} ({(out / 'currents.json').stat().st_size / 1e3:.0f} kB)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
