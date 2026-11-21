#!/usr/bin/env python3
"""Satellite SST and chlorophyll around each real-time glider, for the Real-time map.

Called by fetch_currents.py for every active glider (each source in its own process, like the
current models). For the last SAT_DAYS days it writes:

* one PNG per day and variable (``sat/<glider>_<var>_<date>.png``), drawn on a Web Mercator row
  spacing so it lines up with the Leaflet map, and
* the satellite value at the glider's position each day next to the glider's own near-surface
  value, for the glider pop-up.

SST: MUR SST v4.1 (NASA JPL, 0.01 deg, daily, gap-free foundation temperature) from the
CoastWatch West Coast ERDDAP, with its per-pixel analysis error.

Chlorophyll: Copernicus-GlobColour multi-sensor 4 km daily, near real time. The gap-free (L4)
field is drawn on the map; the observed-only (L3) field on the same grid says which pixels a
satellite actually saw that day. On the map observed pixels are opaque and interpolated ones
faint; at the glider the value is taken from observed pixels where there are any.

Glider values: temperature is the mean of the 0-5 m bin over the day's profiles (0-10 m on days
with no data above 5 m, as in decimated real-time files); chlorophyll is
the mean of the 0-10 m bins over the day's night profiles (sun below the horizon), because
daytime fluorescence is lowered by non-photochemical quenching.
"""
import datetime as dt
import time
import warnings
import urllib.request
from pathlib import Path

import netCDF4
import numpy as np

SAT_DAYS = 7
MATCH_KM = 5.0              # satellite pixels within this distance of the glider's position
NIGHT_SUN_ELEV = 0.0        # degrees; profiles with the sun below this count as night
T_DEPTH, CHL_DEPTH = 5.0, 10.0
T_DEPTH_FALLBACK = 10.0     # real-time files often have no samples above ~6 m (decimated data)

MUR_ERDDAP = "https://coastwatch.pfeg.noaa.gov/erddap/griddap/jplMURSST41"
MUR_NAME = "MUR SST v4.1 (NASA JPL, 0.01°, daily, gap-free), via NOAA CoastWatch ERDDAP"
MUR_URL = "https://podaac.jpl.nasa.gov/dataset/MUR-JPL-L4-GLOB-v4.1"
MUR_STRIDE = 2              # every 2nd 0.01 deg pixel (~2 km)
OSTIA = "METOFFICE-GLO-SST-L4-NRT-OBS-SST-V2"       # Copernicus Marine; fallback when ERDDAP refuses
OSTIA_NAME = "OSTIA SST (UK Met Office, 0.05°, daily, gap-free), via Copernicus Marine (MUR was not available)"
OSTIA_DOI = "https://doi.org/10.48670/moi-00165"
CHL_L4 = "cmems_obs-oc_glo_bgc-plankton_nrt_l4-gapfree-multi-4km_P1D"
CHL_L3 = "cmems_obs-oc_glo_bgc-plankton_nrt_l3-multi-4km_P1D"
CHL_NAME = "Copernicus-GlobColour chlorophyll, multi-sensor 4 km daily (gap-free L4 and observed-only L3)"
CHL_DOI = "https://doi.org/10.48670/moi-00279"
USER_AGENT = "Mozilla/5.0 (X11; Linux x86_64) visualize-glider-data-site (+https://github.com/sciaid/visualize-glider-data-site)"

# cmocean colour maps sampled at 11 stops (Thyng et al. 2016, doi:10.5670/oceanog.2016.66),
# the same stops the curtain viewer uses
CMAPS = {
    "thermal": ["#042333", "#0f326a", "#40349f", "#674396", "#8b538d", "#b15f82", "#d66c6c", "#f3824d", "#fca63c", "#f7d045", "#e8fa5b"],
    "algae": ["#d7f9d0", "#b2e3a5", "#8ecf7c", "#66bb58", "#3ba83f", "#17923b", "#0d7a3f", "#11633b", "#134d33", "#11372a", "#0b2219"],
}


def _f(v):
    a = v[:]
    return np.ma.filled(a.astype("f8"), np.nan) if np.ma.isMaskedArray(a) else np.asarray(a, "f8")


# --------------------------------------------------------------------------- #
# Glider near-surface values per UTC day
# --------------------------------------------------------------------------- #
def sun_elevation(t, lat, lon):
    """Solar elevation (degrees) for unix times t (low-precision almanac, ~0.5 deg)."""
    d = np.asarray(t, "f8") / 86400.0 + 2440587.5 - 2451545.0
    g = np.radians((357.529 + 0.98560028 * d) % 360)
    q = (280.459 + 0.98564736 * d) % 360
    L = np.radians(q + 1.915 * np.sin(g) + 0.020 * np.sin(2 * g))
    e = np.radians(23.439 - 0.00000036 * d)
    ra = np.arctan2(np.cos(e) * np.sin(L), np.cos(L))
    dec = np.arcsin(np.sin(e) * np.sin(L))
    gmst = (18.697374558 + 24.06570982441908 * d) % 24
    ha = np.radians(gmst * 15 + np.asarray(lon, "f8")) - ra
    la = np.radians(np.asarray(lat, "f8"))
    return np.degrees(np.arcsin(np.sin(la) * np.sin(dec) + np.cos(la) * np.cos(dec) * np.cos(ha)))


def glider_daily(t, la, lo, z, T, C, keep=None, sun=None):
    """Per UTC day: the glider's median position and its 0-5 m temperature (all profiles; 0-10 m
    on days with no data above 5 m) and 0-10 m chlorophyll (night profiles).
    t, la, lo: per profile; z: bin centres; T, C: [depth, profile] arrays or None."""
    t, la, lo, z = (np.asarray(a, "f8") for a in (t, la, lo, z))
    ok = np.isfinite(t) & np.isfinite(la) & np.isfinite(lo)
    keep = ok if keep is None else keep & ok
    if not keep.any():
        return {}
    sun = sun_elevation(t, la, lo) if sun is None else np.asarray(sun, "f8")
    night = sun < NIGHT_SUN_ELEV

    def col_mean(a, zmax, use):
        if a is None:
            return None
        a = np.asarray(a, "f8")[z < zmax]                       # depth x profile
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", RuntimeWarning)       # profiles with no shallow data
            prof = np.nanmean(a, axis=0) if a.size else np.full(t.size, np.nan)
        return np.where(use, prof, np.nan)

    T5, T10 = col_mean(T, T_DEPTH, keep), col_mean(T, T_DEPTH_FALLBACK, keep)
    Cn = col_mean(C, CHL_DEPTH, keep & night)
    days = np.array([dt.datetime.fromtimestamp(x, dt.timezone.utc).date().isoformat() if k else ""
                     for x, k in zip(np.nan_to_num(t), keep)])
    out = {}
    for d in sorted(set(days) - {""}):
        m = days == d
        rec = {"lat": round(float(np.median(la[m])), 4), "lon": round(float(np.median(lo[m])), 4),
               "nProfiles": int(m.sum()), "nNight": int((m & night).sum())}
        if T5 is not None and np.isfinite(T5[m]).any():
            rec["T"], rec["Tdepth"] = round(float(np.nanmean(T5[m])), 2), T_DEPTH
        elif T10 is not None and np.isfinite(T10[m]).any():           # no data above 5 m that day
            rec["T"], rec["Tdepth"] = round(float(np.nanmean(T10[m])), 2), T_DEPTH_FALLBACK
        mn = m & night
        if Cn is not None and np.isfinite(Cn[mn]).any():
            rec["chl"] = round(float(np.nanmean(Cn[mn])), 3)
        if mn.any():                                               # chlorophyll is matched at night positions
            rec["latNight"], rec["lonNight"] = round(float(np.median(la[mn])), 4), round(float(np.median(lo[mn])), 4)
        out[d] = rec
    return out


def glider_surface(ds, now):
    """glider_daily for the last SAT_DAYS of a real-time grid file."""
    t = _f(ds["time"])
    get = lambda n: _f(ds[n]) if n in ds.variables else None
    keep = np.isfinite(t) & (t >= now - (SAT_DAYS + 1) * 86400)
    return glider_daily(t, _f(ds["latitude"]), _f(ds["longitude"]), _f(ds["depth"]),
                        get("temperature"), get("chlorophyll"), keep)


# --------------------------------------------------------------------------- #
# Images
# --------------------------------------------------------------------------- #
def _hex(c):
    return [int(c[i:i + 2], 16) for i in (1, 3, 5)]


def _merc(lat):
    return np.log(np.tan(np.pi / 4 + np.radians(lat) / 2))


def write_png(path, lat, lon, values, lo, hi, cmap, log=False, alpha=None):
    """values[lat, lon] -> RGBA PNG, north up, rows resampled to even Web Mercator spacing so a
    Leaflet imageOverlay with the lat/lon bounds lines up. alpha: per-pixel opacity 0..1."""
    from PIL import Image
    v = np.log10(np.clip(values, 1e-3, None)) if log else values
    l0, l1 = (np.log10(lo), np.log10(hi)) if log else (lo, hi)
    x = np.clip((v - l0) / (l1 - l0), 0, 1) * (len(CMAPS[cmap]) - 1)
    stops = np.array([_hex(c) for c in CMAPS[cmap]], "f8")
    i0 = np.clip(np.floor(np.nan_to_num(x)).astype(int), 0, len(stops) - 2)
    fr = (np.nan_to_num(x) - i0)[..., None]
    rgb = stops[i0] * (1 - fr) + stops[i0 + 1] * fr
    a = np.where(np.isfinite(v), 1.0 if alpha is None else alpha, 0.0)
    rgba = np.dstack([rgb, a[..., None] * 255]).astype("u1")
    order = np.argsort(lat)[::-1]                                   # north first
    lat, rgba = lat[order], rgba[order]
    my = np.linspace(_merc(lat[0]), _merc(lat[-1]), lat.size)       # even Mercator rows
    rows = np.abs(_merc(lat)[None, :] - my[:, None]).argmin(axis=1)
    Image.fromarray(rgba[rows], "RGBA").save(path, optimize=True)


def color_range(stack, log=False, pad=0.0):
    v = stack[np.isfinite(stack)]
    if log:
        v = v[v > 0]
    if not v.size:
        return None
    lo, hi = np.percentile(v, [2, 98])
    if hi <= lo:
        hi = lo + (0.1 if not log else lo)
    return [round(float(lo - pad), 3), round(float(hi + pad), 3)]


# --------------------------------------------------------------------------- #
# Values at the glider
# --------------------------------------------------------------------------- #
def near(lat, lon, la, lo, km=MATCH_KM):
    """Boolean [lat, lon] mask of grid cells within km of (la, lo)."""
    dy = (lat[:, None] - la) * 111.2
    dx = (lon[None, :] - lo) * 111.2 * np.cos(np.radians(la))
    return np.hypot(dy, dx) <= km


def _date(t):
    return dt.datetime.fromtimestamp(float(t), dt.timezone.utc).date().isoformat()


# --------------------------------------------------------------------------- #
# MUR SST
# --------------------------------------------------------------------------- #
def _get(url, path, tries=3):
    req = urllib.request.Request(url, headers={"User-Agent": USER_AGENT})
    for k in range(tries):
        try:
            with urllib.request.urlopen(req, timeout=180) as r, open(path, "wb") as f:
                f.write(r.read())
            return
        except Exception:
            if k == tries - 1:
                raise
            time.sleep(10 * (k + 1))


def sst_layer(g, box, out, workdir):
    """MUR from the CoastWatch ERDDAP; if that server refuses (it sometimes answers 403 to the
    build machines), OSTIA from Copernicus Marine instead."""
    try:
        return mur_layer(g, box, out, workdir)
    except Exception as exc:
        print(f"    MUR not available ({type(exc).__name__}: {exc}); using OSTIA", flush=True)
        return ostia_layer(g, box, out, workdir)


def ostia_layer(g, box, out, workdir):
    import copernicusmarine
    lat0, lat1, lon0, lon1 = box
    today = dt.datetime.now(dt.timezone.utc).date()
    fn = f"{g['name']}_ostia.nc"
    copernicusmarine.subset(
        dataset_id=OSTIA, variables=["analysed_sst", "analysis_error"],
        minimum_longitude=lon0, maximum_longitude=lon1, minimum_latitude=lat0, maximum_latitude=lat1,
        start_datetime=f"{today - dt.timedelta(days=SAT_DAYS + 2)}T00:00:00", end_datetime=f"{today}T23:59:59",
        output_directory=str(workdir), output_filename=fn, disable_progress_bar=True)
    with netCDF4.Dataset(Path(workdir) / fn) as ds:
        tv = ds["time"]
        tt = netCDF4.num2date(tv[:], tv.units, only_use_cftime_datetimes=False, only_use_python_datetimes=True)
        t = np.array([x.replace(tzinfo=dt.timezone.utc).timestamp() for x in tt])
        sst, err = _f(ds["analysed_sst"]), _f(ds["analysis_error"])
        if getattr(ds["analysed_sst"], "units", "K").lower().startswith("k") or np.nanmedian(sst) > 200:
            sst = sst - 273.15
        lat, lon = _f(ds["latitude"]), _f(ds["longitude"])
    keep = np.argsort(t)[-SAT_DAYS:]
    b = sat_block(g, "sst", t[keep], lat, lon, sst[keep], err[keep], None, out)
    b.update(source=OSTIA_NAME, url=OSTIA_DOI)
    return b


def mur_layer(g, box, out, workdir):
    lat0, lat1, lon0, lon1 = box
    s = MUR_STRIDE
    # time by index (last-6 .. last); "(last-6)" in parentheses would mean 6 seconds before the last time
    sub = f"[last-{SAT_DAYS - 1}:1:last][({lat0:.2f}):{s}:({lat1:.2f})][({lon0:.2f}):{s}:({lon1:.2f})]"
    url = f"{MUR_ERDDAP}.nc?analysed_sst{sub},analysis_error{sub}"
    fn = Path(workdir) / f"{g['name']}_mur.nc"
    _get(url, fn)
    with netCDF4.Dataset(fn) as ds:
        t, lat, lon = _f(ds["time"]), _f(ds["latitude"]), _f(ds["longitude"])
        sst, err = _f(ds["analysed_sst"]), _f(ds["analysis_error"])
    return sat_block(g, "sst", t, lat, lon, sst, err, None, out)


# --------------------------------------------------------------------------- #
# Chlorophyll
# --------------------------------------------------------------------------- #
def _cm_read(dataset, g, box, workdir, tag):
    import copernicusmarine
    lat0, lat1, lon0, lon1 = box
    today = dt.datetime.now(dt.timezone.utc).date()
    fn = f"{g['name']}_{tag}.nc"
    copernicusmarine.subset(
        dataset_id=dataset, variables=["CHL"],
        minimum_longitude=lon0, maximum_longitude=lon1, minimum_latitude=lat0, maximum_latitude=lat1,
        start_datetime=f"{today - dt.timedelta(days=SAT_DAYS + 3)}T00:00:00", end_datetime=f"{today}T23:59:59",
        output_directory=str(workdir), output_filename=fn, disable_progress_bar=True)
    with netCDF4.Dataset(Path(workdir) / fn) as ds:
        tv = ds["time"]
        tt = netCDF4.num2date(tv[:], tv.units, only_use_cftime_datetimes=False, only_use_python_datetimes=True)
        t = np.array([x.replace(tzinfo=dt.timezone.utc).timestamp() for x in tt])
        return t, _f(ds["latitude"]), _f(ds["longitude"]), _f(ds["CHL"])


def chl_layer(g, box, out, workdir):
    t, lat, lon, l4 = _cm_read(CHL_L4, g, box, workdir, "chl_l4")
    obs = np.zeros(l4.shape, bool)
    note = None
    try:
        t3, lat3, lon3, l3 = _cm_read(CHL_L3, g, box, workdir, "chl_l3")
        if l3.shape[1:] == l4.shape[1:] and np.allclose(lat3, lat) and np.allclose(lon3, lon):
            d3 = [_date(x) for x in t3]
            for k, x in enumerate(t):
                if _date(x) in d3:
                    obs[k] = np.isfinite(l3[d3.index(_date(x))])
        else:
            note = "observed-only (L3) grid does not match the gap-free grid"
    except Exception as exc:                       # the map still shows the gap-free field
        note = f"observed-only (L3) field not available ({type(exc).__name__})"
    keep = np.argsort(t)[-SAT_DAYS:]
    b = sat_block(g, "chl", t[keep], lat, lon, l4[keep], None, obs[keep], out)
    if note:
        b["note"] = note
    return b


# --------------------------------------------------------------------------- #
# Common: images + values at the glider
# --------------------------------------------------------------------------- #
def sat_block(g, var, t, lat, lon, field, err, obs, out):
    """field[time, lat, lon]. err: SST analysis error (same shape) or None. obs: chlorophyll
    observed mask (same shape) or None."""
    log = var == "chl"
    rng = color_range(field, log=log)
    if rng is None:
        raise ValueError(f"no {var} values in the box")
    if log:
        rng = [max(rng[0], 0.01), max(rng[1], rng[0] * 2)]
    sat_dir = Path(out) / "sat"
    sat_dir.mkdir(parents=True, exist_ok=True)
    surf = g.get("surface", {})
    days = []
    for k in np.argsort(t):
        d = _date(t[k])
        f = field[k]
        if not np.isfinite(f).any():
            continue
        img = f"sat/{g['name']}_{var}_{d}.png"
        alpha = None if obs is None else np.where(obs[k], 0.9, 0.4)
        write_png(sat_dir / Path(img).name, lat, lon, f, rng[0], rng[1],
                  "algae" if log else "thermal", log=log, alpha=alpha)
        day = {"date": d, "img": img}
        if obs is not None:
            day["observedFrac"] = round(float(obs[k][np.isfinite(f)].mean()), 3)
        # value at the glider: its median (night, for chlorophyll) position that day, or the
        # latest fix when it has no data for that day
        s = surf.get(d, {})
        la = s.get("latNight" if log else "lat", s.get("lat", g["lat"]))
        lo = s.get("lonNight" if log else "lon", s.get("lon", g["lon"]))
        m = near(lat, lon, la, lo) & np.isfinite(f)
        if not m.any():                                          # pixel coarser than the radius
            r, q = np.abs(lat - la).argmin(), np.abs(lon - lo).argmin()
            m = np.zeros(f.shape, bool); m[r, q] = np.isfinite(f[r, q])
        if m.any():
            at = {"lat": la, "lon": lo, "fromGlider": bool(s)}
            if obs is not None and (m & obs[k]).any():
                at.update(v=round(float(np.median(f[m & obs[k]])), 3), observed=True)
            else:
                at.update(v=round(float(np.median(f[m])), 3 if log else 2))
                if obs is not None:
                    at["observed"] = False
            if err is not None and np.isfinite(err[k][m]).any():
                at["err"] = round(float(np.nanmedian(err[k][m])), 2)
            day["atGlider"] = at
        days.append(day)
    hy, hx = abs(float(np.median(np.diff(lat)))) / 2, abs(float(np.median(np.diff(lon)))) / 2   # pixel edges
    return {"days": days, "range": rng, "log": log,
            "bounds": [[float(lat.min()) - hy, float(lon.min()) - hx], [float(lat.max()) + hy, float(lon.max()) + hx]],
            "source": MUR_NAME if var == "sst" else CHL_NAME,
            "url": MUR_URL if var == "sst" else CHL_DOI}
