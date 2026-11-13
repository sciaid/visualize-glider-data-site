#!/usr/bin/env python3
"""
glider_curtain.py - Interactive 3D "curtain" viewer for one gridded glider deployment.

Builds a self-contained HTML page from a PyGlider-style gridded NetCDF file
(dimensions depth x time, 1-D time/latitude/longitude per profile), e.g.
amlr04-20221205_grid-delayed-5m.nc.

    python glider_curtain.py FILE.nc [FILE2.nc ...] [-o OUT_DIR] [options]

One page per file: <deployment>_curtain.html. Open it in a web browser.

The page shows
  * a 3D curtain hanging below the glider track (km east/north x depth),
  * a 2D time-depth section of the same variable,
  * buttons to switch variable (default: temperature, salinity, chlorophyll, oxygen),
  * colour map, colour range, log scale (default on for chlorophyll),
  * vertical exaggeration and a date-window slider (repeat transects overlap
    in 3D, so narrow the window to see one pass at a time),
  * a readout of time / position / depth / all variables under the cursor.

Requirements: pip install netCDF4 numpy
Optional:     pip install plotly   -> the plotting library is embedded in the
              page, so it works offline. Without it the page loads plotly.js
              from the internet.

Notes
  * Gridded QC variables (*_qc) are bin-averaged flags in these files and are
    not used for masking; values are shown as stored.
"""

from __future__ import annotations

import argparse
import base64
import datetime as dt
import json
import sys
from pathlib import Path

import numpy as np

try:
    import netCDF4
except ImportError:  # pragma: no cover
    sys.exit("netCDF4 is required:  pip install netCDF4 numpy")

PLOTLY_CDN = "https://cdn.plot.ly/plotly-2.35.2.min.js"

# Display defaults per variable: (label, colour map, log scale)
VAR_DEFAULTS = {
    "temperature": ("Temperature", "thermal", False),
    "potential_temperature": ("Potential temperature", "thermal", False),
    "salinity": ("Salinity", "haline", False),
    "chlorophyll": ("Chlorophyll", "algae", True),
    "oxygen_concentration": ("Oxygen", "deep", False),
    "oxygen_saturation": ("Oxygen saturation", "deep", False),
    "density": ("Density", "dense", False),
    "potential_density": ("Potential density", "dense", False),
    "backscatter_700": ("Backscatter 700 nm", "tempo", True),
}

# cmocean colour maps sampled at 11 stops (Thyng et al. 2016, doi:10.5670/oceanog.2016.66)
CMAPS = {
    "thermal": ["#042333", "#0f326a", "#40349f", "#674396", "#8b538d", "#b15f82", "#d66c6c", "#f3824d", "#fca63c", "#f7d045", "#e8fa5b"],
    "haline": ["#2a186c", "#2927a2", "#0d4e96", "#18668c", "#2b7b88", "#3a9087", "#4aa584", "#62ba7b", "#89cd6d", "#bfdc6c", "#fdef9a"],
    "algae": ["#d7f9d0", "#b2e3a5", "#8ecf7c", "#66bb58", "#3ba83f", "#17923b", "#0d7a3f", "#11633b", "#134d33", "#11372a", "#0b2219"],
    "oxy": ["#400505", "#6f0c0b", "#930d0d", "#606060", "#7c7c7b", "#999998", "#b8b7b6", "#d8d7d5", "#ebe952", "#dbbf23", "#dcaf00"],
    "deep": ["#fdfecc", "#c3eaa9", "#8cd5a1", "#62bba2", "#4ea0a0", "#43859d", "#3e6a98", "#424f8d", "#3e3b6a", "#2f2a47", "#281a2c"],
    "dense": ["#e6f1f1", "#bbd8e0", "#93bfd8", "#749fd5", "#6f7dcf", "#7059b8", "#6c3a96", "#5f2171", "#4a1550", "#300e31", "#360e24"],
    "tempo": ["#fff6f4", "#d8ddcb", "#b1c6a6", "#87b08d", "#5d9b7f", "#348378", "#1c6b6f", "#185263", "#183b52", "#162540", "#151d44"],
    "viridis": ["#440154", "#482878", "#3e4989", "#31688e", "#26828e", "#1f9e89", "#35b779", "#6ece58", "#b5de2b", "#fde725", "#fde725"],
    # jet-style rainbow, offered because many glider portals default to it; not perceptually uniform
    "rainbow": ["#00007f", "#0000ff", "#0063ff", "#00d4ff", "#4effa9", "#a9ff4e", "#ffe600", "#ff7a00", "#ff1400", "#b20000", "#7f0000"],
}


def _f(var):
    a = var[:]
    return np.ma.filled(a.astype("f8"), np.nan) if np.ma.isMaskedArray(a) else np.asarray(a, "f8")


def _b64(a):
    return base64.b64encode(np.ascontiguousarray(a, dtype="<f4").tobytes()).decode("ascii")


def load(path, var_names, profile_step, bathy=None):
    ds = netCDF4.Dataset(path)
    for need in ("time", "latitude", "longitude", "depth"):
        if need not in ds.variables:
            sys.exit(f"{path}: missing variable '{need}' - expected a gridded (depth x time) glider file")
    t = _f(ds["time"]); lat = _f(ds["latitude"]); lon = _f(ds["longitude"]); depth = _f(ds["depth"])
    units = getattr(ds["time"], "units", "")
    if units and not units.lower().startswith("seconds since 1970-01-01"):
        dates = netCDF4.num2date(t, units, only_use_cftime_datetimes=False, only_use_python_datetimes=True)
        t = netCDF4.date2num(dates, "seconds since 1970-01-01", calendar="standard")

    good = np.isfinite(t) & np.isfinite(lat) & np.isfinite(lon)
    idx = np.where(good)[0][::max(1, profile_step)]
    idx = idx[np.argsort(t[idx], kind="stable")]

    variables = []
    deepest = 0
    for name in var_names:
        if name not in ds.variables:
            print(f"  note: {name} not in file, skipped")
            continue
        v = ds[name]
        a = _f(v)
        if a.shape != (depth.size, t.size):
            if a.shape == (t.size, depth.size):
                a = a.T
            else:
                print(f"  note: {name} has shape {a.shape}, expected (depth, time); skipped")
                continue
        a = a[:, idx]
        rows = np.where(np.isfinite(a).any(axis=1))[0]
        if rows.size == 0:
            print(f"  note: {name} is all missing, skipped")
            continue
        deepest = max(deepest, rows.max())
        label, cmap, log = VAR_DEFAULTS.get(name, (getattr(v, "long_name", name), "viridis", False))
        units = getattr(v, "units", "")
        if name == "salinity" and units.strip() in ("1e-3", "0.001", "1", "PSU", "psu", "PSS-78"):
            units = "PSU"
        variables.append({"name": name, "label": label, "units": units,
                          "cmap": cmap, "log": log, "arr": a})
    if not variables:
        sys.exit(f"{path}: none of {var_names} found")

    nd = deepest + 1
    lat0, lon0 = float(np.nanmean(lat[idx])), float(np.nanmean(lon[idx]))
    kx = 111.32 * np.cos(np.radians(lat0))
    x = (lon[idx] - lon0) * kx
    y = (lat[idx] - lat0) * 110.57

    meta = {a: str(getattr(ds, a)) for a in ("deployment_name", "glider_name", "institution",
                                              "time_coverage_start", "time_coverage_end") if hasattr(ds, a)}
    name = meta.get("deployment_name") or Path(path).stem
    data = {
        "name": name, "file": Path(path).name, "meta": meta,
        "nProf": int(idx.size), "nDepth": int(nd),
        "depth": [round(float(d), 2) for d in depth[:nd]],
        "time": [int(round(float(s))) for s in t[idx]],
        "lat": [round(float(s), 5) for s in lat[idx]],
        "lon": [round(float(s), 5) for s in lon[idx]],
        "x": _b64(x), "y": _b64(y),
        "vars": [],
    }
    data["aux"] = derived_block(ds, idx, depth, nd, t[idx], lat[idx], lon[idx])
    if bathy is not None:
        data["bathy"] = bathy_block(bathy, lat[idx], lon[idx], lat0, lon0, kx)
    for v in variables:
        arr = v["arr"][:nd]
        fin = arr[np.isfinite(arr)]
        lo, hi = np.percentile(fin, [1, 99])
        pos = fin[fin > 0]
        llo, lhi = (np.percentile(pos, [1, 99]) if pos.size else (lo, hi))
        data["vars"].append({
            "name": v["name"], "label": v["label"], "units": v["units"],
            "cmap": v["cmap"], "log": v["log"],
            "range": [float(f"{lo:.4g}"), float(f"{hi:.4g}")],
            "logRange": [float(f"{llo:.4g}"), float(f"{lhi:.4g}")],
            "data": _b64(arr),  # row-major: depth x profile
        })
    ds.close()
    return data


# --------------------------------------------------------------------------- #
# Derived quantities
# --------------------------------------------------------------------------- #
def sigma0_eos80(S, theta):
    """Potential density anomaly (kg m-3) at p=0, UNESCO EOS-80 (Millero & Poisson 1981).
    Used only for background isopycnals on the T-S diagram."""
    T = theta * 1.00024  # ITS-90 -> IPTS-68
    rw = (999.842594 + 6.793952e-2 * T - 9.095290e-3 * T**2 + 1.001685e-4 * T**3
          - 1.120083e-6 * T**4 + 6.536332e-9 * T**5)
    A = 8.24493e-1 - 4.0899e-3 * T + 7.6438e-5 * T**2 - 8.2467e-7 * T**3 + 5.3875e-9 * T**4
    B = -5.72466e-3 + 1.0227e-4 * T - 1.6546e-6 * T**2
    return rw + A * S + B * np.abs(S) ** 1.5 + 4.8314e-4 * S**2 - 1000.0


def sun_elevation(t_epoch, lat, lon):
    """Solar elevation angle (deg, no refraction), NOAA/Astronomical Almanac low-precision
    formulae; accurate to ~0.01 deg, ample for day/night classification."""
    n = t_epoch / 86400.0 + 2440587.5 - 2451545.0
    L = np.radians((280.460 + 0.9856474 * n) % 360)
    g = np.radians((357.528 + 0.9856003 * n) % 360)
    lam = L + np.radians(1.915) * np.sin(g) + np.radians(0.020) * np.sin(2 * g)
    eps = np.radians(23.439 - 4e-7 * n)
    ra = np.arctan2(np.cos(eps) * np.sin(lam), np.cos(lam))
    dec = np.arcsin(np.sin(eps) * np.sin(lam))
    gmst = np.radians((280.46061837 + 360.98564736629 * n) % 360)
    ha = gmst + np.radians(lon) - ra
    la = np.radians(lat)
    return np.degrees(np.arcsin(np.sin(la) * np.sin(dec) + np.cos(la) * np.cos(dec) * np.cos(ha)))


def mixed_layer_depth(sig, depth, ref=10.0, dsig=0.03):
    """Density-threshold MLD (de Boyer Montegut et al. 2004): depth where sigma-theta first
    exceeds its value at `ref` m by `dsig` kg m-3. NaN where the profile lacks near-surface
    data or never crosses the threshold within the sampled depth."""
    out = np.full(sig.shape[1], np.nan)
    for i in range(sig.shape[1]):
        s = sig[:, i]; f = np.isfinite(s)
        if f.sum() < 3:
            continue
        d, v = depth[f], s[f]
        if d[0] > ref + 5:
            continue
        sref = np.interp(ref, d, v) if d[0] <= ref else v[0]
        below = np.where((d > ref) & (v - sref > dsig))[0]
        if below.size == 0:
            continue
        k = below[0]
        if k == 0:
            out[i] = d[0]; continue
        d0, d1, v0, v1 = d[k - 1], d[k], v[k - 1] - sref, v[k] - sref
        out[i] = d0 + (dsig - v0) * (d1 - d0) / (v1 - v0) if v1 != v0 else d1
    return out


def chl_max_depth(chl, depth, max_depth=300.0, top=20.0, min_bins=8):
    """Depth of the (3-bin smoothed) chlorophyll maximum in the upper `max_depth` m.
    Requires data starting within `top` m. Returns (depth, value)."""
    keep = depth <= max_depth
    dz, c = depth[keep], chl[keep]
    dep = np.full(chl.shape[1], np.nan); val = np.full(chl.shape[1], np.nan)
    for i in range(c.shape[1]):
        x = c[:, i]; f = np.isfinite(x)
        if f.sum() < min_bins or dz[f][0] > top:
            continue
        xs = x.copy()
        num = np.convolve(np.where(f, x, 0), np.ones(3), "same"); den = np.convolve(f.astype(float), np.ones(3), "same")
        xs[f] = (num / np.maximum(den, 1))[f]
        k = int(np.nanargmax(np.where(f, xs, -np.inf)))
        dep[i], val[i] = dz[k], x[k] if np.isfinite(x[k]) else xs[k]
    return dep, val


def sigma_contour_lines(sig, dz, t, start, end, step, smooth=5):
    """Pre-compute sigma-theta contour lines on the (profile, depth) grid so the page does not
    have to contour on every redraw. Lightly smoothed along time (display only)."""
    try:
        import contourpy
    except ImportError:
        return None
    f = np.isfinite(sig).astype(float)
    k = np.ones(smooth)
    num = np.apply_along_axis(lambda r: np.convolve(r, k, "same"), 1, np.where(f > 0, sig, 0.0))
    den = np.apply_along_axis(lambda r: np.convolve(r, k, "same"), 1, f)
    sm = np.ma.masked_where(f == 0, num / np.maximum(den, 1))
    xi = np.arange(sig.shape[1], dtype=float)
    gen = contourpy.contour_generator(xi, dz, sm)
    tt, dd, labels = [], [], []
    for lev in np.arange(start, end + step / 2, step):
        best = None
        for ln in gen.lines(float(lev)):
            if len(ln) < 4:
                continue
            ln = ln[:: max(1, len(ln) // 400)] if len(ln) > 800 else ln
            tt.extend(np.round(np.interp(ln[:, 0], xi, t)).astype(int).tolist() + [None])
            dd.extend(np.round(ln[:, 1], 1).tolist() + [None])
            if best is None or len(ln) > len(best):
                best = ln
        if best is not None:
            m = best[len(best) // 2]
            labels.append({"t": int(round(float(np.interp(m[0], xi, t)))), "d": round(float(m[1]), 1),
                           "txt": f"{lev:.2f}".rstrip("0").rstrip(".")})
    return {"t": tt, "d": dd, "labels": labels} if tt else None


def _nice_step(span, n=8):
    for s in (0.01, 0.02, 0.05, 0.1, 0.2, 0.25, 0.5, 1.0, 2.0):
        if span / s <= n:
            return s
    return 5.0


def derived_block(ds, idx, depth, nd, t, lat, lon):
    v = ds.variables
    grab = lambda name: (_f(v[name])[:, idx][:nd] if name in v and v[name].ndim == 2 else None)
    theta = grab("potential_temperature")
    if theta is None:
        theta = grab("temperature")
    sal = grab("salinity")
    sig = grab("potential_density")
    if sig is not None and np.nanmedian(sig) > 500:
        sig = sig - 1000.0
    if sig is None and theta is not None and sal is not None:
        sig = sigma0_eos80(sal, theta)
    chl = grab("chlorophyll")
    out = {"sun": [round(float(x), 2) for x in sun_elevation(t, lat, lon)]}
    dz = depth[:nd]
    if sig is not None:
        out["sigma"] = _b64(sig)
        lo, hi = np.nanpercentile(sig, [2, 98])
        step = _nice_step(hi - lo)
        out["sigStep"] = step
        out["sigStart"] = float(np.floor(np.nanmin(sig) / step) * step)
        out["sigEnd"] = float(np.ceil(np.nanmax(sig) / step) * step)
        lines = sigma_contour_lines(sig, dz, t, out["sigStart"], out["sigEnd"], step)
        if lines:
            out["sigLines"] = lines
        mld = mixed_layer_depth(sig, dz)
        out["mld"] = [None if not np.isfinite(x) else round(float(x), 1) for x in mld]
    if chl is not None:
        cd, cv = chl_max_depth(chl, dz)
        out["chlmax"] = [None if not np.isfinite(x) else round(float(x), 1) for x in cd]
        out["chlmaxVal"] = [None if not np.isfinite(x) else float(f"{x:.3g}") for x in cv]
    if theta is not None and sal is not None:
        out["theta"] = _b64(theta); out["sal"] = _b64(sal)
        f = np.isfinite(theta) & np.isfinite(sal)
        s0, s1 = np.percentile(sal[f], [0.2, 99.8]); t0, t1 = np.percentile(theta[f], [0.2, 99.8])
        ps, pt = 0.08 * (s1 - s0) + 0.01, 0.08 * (t1 - t0) + 0.05
        sg = np.linspace(s0 - ps, s1 + ps, 90); tg = np.linspace(t0 - pt, t1 + pt, 90)
        SS, TT = np.meshgrid(sg, tg)
        z = sigma0_eos80(SS, TT)
        zl, zh = np.nanmin(z), np.nanmax(z); st = _nice_step(zh - zl, 10)
        out["ts"] = {"s": [round(float(x), 4) for x in sg], "t": [round(float(x), 4) for x in tg],
                     "z": [[round(float(x), 4) for x in row] for row in z],
                     "step": st, "start": float(np.ceil(zl / st) * st), "end": float(np.floor(zh / st) * st),
                     "srange": [float(sg[0]), float(sg[-1])], "trange": [float(tg[0]), float(tg[-1])],
                     "tname": "θ" if "potential_temperature" in v else "T"}
    return out


# --------------------------------------------------------------------------- #
# Bathymetry / topography (GEBCO or ETOPO NetCDF grid)
# --------------------------------------------------------------------------- #
class Bathy:
    """Lazy wrapper around a lat/lon elevation grid (GEBCO 'elevation', ETOPO 'z')."""
    def __init__(self, path):
        self.path = Path(path)
        self.ds = netCDF4.Dataset(path)
        v = self.ds.variables
        self.lat = _f(v[next(k for k in ("lat", "latitude", "y") if k in v)])
        self.lon = _f(v[next(k for k in ("lon", "longitude", "x") if k in v)])
        self.zname = next((k for k in ("elevation", "z", "Band1", "altitude", "height") if k in v), None)
        if self.zname is None:
            sys.exit(f"{path}: no elevation variable (expected 'elevation' or 'z')")
        g = self.ds
        self.source = str(getattr(g, "title", self.path.name)).split(" - ")[0].replace("The ", "", 1)
        self.doi = next((str(getattr(g, k)).replace("DOI: ", "doi:") for k in
                         ("identifier_product_doi", "doi", "DOI", "references") if hasattr(g, k)), "")
        self.z = None   # full grid loaded on first use

    def _grid(self):
        if self.z is None:
            self.z = _f(self.ds[self.zname])
        return self.z

    def subset(self, la0, la1, lo0, lo1):
        i = np.where((self.lat >= la0) & (self.lat <= la1))[0]
        j = np.where((self.lon >= lo0) & (self.lon <= lo1))[0]
        if i.size < 2 or j.size < 2:
            return None
        z = self._grid()[i.min():i.max() + 1, j.min():j.max() + 1]
        return self.lat[i.min():i.max() + 1], self.lon[j.min():j.max() + 1], z

    def sample(self, lat, lon):
        """Nearest-cell elevation at points (full resolution)."""
        def nearest(grid, vals):
            k = np.clip(np.searchsorted(grid, vals), 1, grid.size - 1)
            return np.where(np.abs(vals - grid[k - 1]) <= np.abs(vals - grid[k]), k - 1, k)
        asc_lat = self.lat[0] < self.lat[-1]
        glat = self.lat if asc_lat else self.lat[::-1]
        ii = nearest(glat, lat); ii = ii if asc_lat else self.lat.size - 1 - ii
        jj = nearest(self.lon, lon)
        return self._grid()[ii, jj].astype(float)


def find_bathy(paths):
    """Auto-detect a GEBCO/ETOPO grid next to the glider files."""
    for d in {Path(p).resolve().parent for p in paths}:
        for q in sorted(list(d.glob("*.nc")) + list(d.glob("*/*.nc"))):
            if any(k in q.name.lower() for k in ("gebco", "etopo", "bathy")):
                return q
    return None


def bathy_block(bathy, lat, lon, lat0, lon0, kx, max_cells=220):
    la0, la1, lo0, lo1 = np.nanmin(lat), np.nanmax(lat), np.nanmin(lon), np.nanmax(lon)
    mla = max(0.12, 0.12 * (la1 - la0)); mlo = max(0.25, 0.12 * (lo1 - lo0))
    sub = bathy.subset(la0 - mla, la1 + mla, lo0 - mlo, lo1 + mlo)
    if sub is None:
        print("  note: bathymetry grid does not cover this deployment; skipped")
        return None
    blat, blon, z = sub
    if blat[0] > blat[-1]:
        blat, z = blat[::-1], z[::-1]
    coast = _coastline(blat, blon, z, lat0, lon0, kx)
    sy = max(1, int(np.ceil(blat.size / max_cells))); sx = max(1, int(np.ceil(blon.size / max_cells)))
    zb = z[::sy, ::sx]; blat_s = blat[::sy]; blon_s = blon[::sx]
    bottom = -bathy.sample(lat, lon)                                   # positive depth (m)
    sea = zb[zb < 0]; land = zb[zb > 0]
    return {
        "nx": int(blon_s.size), "ny": int(blat_s.size),
        "x": [round(float(v), 3) for v in (blon_s - lon0) * kx],
        "y": [round(float(v), 3) for v in (blat_s - lat0) * 110.57],
        "z": _b64(zb),
        "maxDepth": float(-sea.min()) if sea.size else 0.0,
        "maxLand": float(land.max()) if land.size else 0.0,
        "bottom": _b64(bottom),
        "coast": coast,
        "source": bathy.source, "doi": bathy.doi,
        "res": f"{abs(float(blat[1] - blat[0])) * 3600:.0f} arc-sec grid, shown every {max(sx, sy)} cells",
    }


def _coastline(blat, blon, z, lat0, lon0, kx):
    try:
        import contourpy
    except ImportError:
        return None
    gen = contourpy.contour_generator(blon, blat, np.nan_to_num(z, nan=-1.0))
    xs, ys = [], []
    for line in gen.lines(0.0):
        if len(line) < 3:
            continue
        xs.extend(((line[:, 0] - lon0) * kx).tolist() + [None])
        ys.extend(((line[:, 1] - lat0) * 110.57).tolist() + [None])
    return {"x": [None if v is None else round(v, 3) for v in xs],
            "y": [None if v is None else round(v, 3) for v in ys]}


def plotly_tag(mode):
    if mode != "cdn":
        try:
            import plotly  # noqa
            js = Path(plotly.__file__).parent / "package_data" / "plotly.min.js"
            if js.exists():
                return "<script>" + js.read_text(encoding="utf-8") + "</script>"
        except ImportError:
            pass
        if mode == "inline":
            sys.exit("--plotly inline needs the plotly Python package:  pip install plotly")
    return f'<script src="{PLOTLY_CDN}"></script>'


def build_site(a, files, bathy):
    site = Path(a.site); (site / "data").mkdir(parents=True, exist_ok=True)
    try:
        import plotly  # noqa
        js = Path(plotly.__file__).parent / "package_data" / "plotly.min.js"
        (site / "plotly.min.js").write_bytes(js.read_bytes())
        tag = '<script src="plotly.min.js"></script>'
    except ImportError:
        tag = f'<script src="{PLOTLY_CDN}"></script>'
    (site / "index.html").write_text(HTML.replace("__PLOTLY__", tag).replace("__CMAPS__", json.dumps(CMAPS))
                                     .replace("__TITLE__", "Glider curtain viewer").replace("__DATA__", "null"),
                                     encoding="utf-8")
    index = []
    idx_path = site / "data" / "index.json"
    if idx_path.exists():                       # keep deployments built earlier
        index = [e for e in json.loads(idx_path.read_text()) if (site / "data" / f"{e['name']}.json").exists()]
    for f in files:
        print(f)
        d = load(f, [s.strip() for s in a.vars.split(",") if s.strip()], a.profile_step, bathy)
        (site / "data" / f"{d['name']}.json").write_text(json.dumps(d, separators=(",", ":")), encoding="utf-8")
        start = dt.datetime.fromtimestamp(d["time"][0], dt.timezone.utc)
        end = dt.datetime.fromtimestamp(d["time"][-1], dt.timezone.utc)
        entry = {"name": d["name"], "start": start.strftime("%Y-%m-%d"), "end": end.strftime("%Y-%m-%d"),
                 "profiles": d["nProf"], "maxDepth": d["depth"][-1], "vars": [v["name"] for v in d["vars"]],
                 "glider": d["meta"].get("glider_name", ""),
                 "latRange": [min(d["lat"]), max(d["lat"])], "lonRange": [min(d["lon"]), max(d["lon"])],
                 "bathy": bool(d.get("bathy"))}
        index = [e for e in index if e["name"] != d["name"]] + [entry]
        print(f"  {d['nProf']} profiles -> {site / 'data' / (d['name'] + '.json')}")
    index.sort(key=lambda e: e["start"])
    idx_path.write_text(json.dumps(index, indent=1), encoding="utf-8")
    print(f"Viewer: {site / 'index.html'}  ({len(index)} deployments in data/index.json)")


def main(argv=None):
    ap = argparse.ArgumentParser(description="3D curtain viewer for gridded glider NetCDF files")
    ap.add_argument("files", nargs="+")
    ap.add_argument("-o", "--out-dir", default="glider_output")
    ap.add_argument("--vars", default="temperature,salinity,chlorophyll,oxygen_concentration",
                    help="comma-separated variable names (default: %(default)s)")
    ap.add_argument("--profile-step", type=int, default=1,
                    help="keep every Nth profile (default 1 = all)")
    ap.add_argument("--plotly", choices=["auto", "inline", "cdn"], default="auto",
                    help="auto: embed plotly.js if the plotly package is installed, else CDN")
    ap.add_argument("--site", metavar="DIR",
                    help="write a shared viewer (DIR/index.html + plotly.min.js) and one data file per "
                         "deployment (DIR/data/<name>.json) for a website, instead of standalone pages")
    ap.add_argument("--bathy", default="auto",
                    help="GEBCO/ETOPO NetCDF grid for seafloor + land; 'auto' looks for a file named "
                         "*gebco*/*etopo*/*bathy*.nc next to the glider files; 'none' to skip")
    a = ap.parse_args(argv)

    files = [f for f in a.files if not any(k in Path(f).name.lower() for k in ("gebco", "etopo", "bathy"))]
    out = Path(a.out_dir); out.mkdir(parents=True, exist_ok=True)
    tag = plotly_tag(a.plotly)
    bpath = None if a.bathy == "none" else (find_bathy(files) if a.bathy == "auto" else Path(a.bathy))
    bathy = Bathy(bpath) if bpath else None
    print(f"Bathymetry: {bpath if bpath else 'none'}")
    if a.site:
        return build_site(a, files, bathy)
    for f in files:
        print(f"{f}")
        data = load(f, [s.strip() for s in a.vars.split(",") if s.strip()], a.profile_step, bathy)
        html = (HTML.replace("__PLOTLY__", tag)
                .replace("__CMAPS__", json.dumps(CMAPS))
                .replace("__TITLE__", data["name"])
                .replace("__DATA__", json.dumps(data, separators=(",", ":"))))
        p = out / f"{data['name']}_curtain.html"
        p.write_text(html, encoding="utf-8")
        print(f"  {data['nProf']} profiles x {data['nDepth']} depth bins, "
              f"{', '.join(v['name'] for v in data['vars'])} -> {p}")


HTML = r"""<!DOCTYPE html>
<html lang="en"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>__TITLE__ - curtain</title>
__PLOTLY__
<style>
  :root { --surface:#fcfcfb; --surface-2:#f3f2ef; --border:#e2e0db; --ink:#0b0b0b; --ink-2:#52514e; --muted:#7a7873; --accent:#2a78d6; }
  * { box-sizing:border-box; }
  html,body { margin:0; height:100%; font:13px/1.4 system-ui,-apple-system,"Segoe UI",sans-serif; color:var(--ink); background:var(--surface); }
  #app { display:grid; grid-template-columns:300px 1fr; height:100%; }
  aside { border-right:1px solid var(--border); overflow:auto; padding:14px 16px; }
  h1 { font-size:16px; margin:0 0 2px; }
  .muted { color:var(--muted); font-size:12px; }
  h2 { font-size:11px; text-transform:uppercase; letter-spacing:.05em; color:var(--ink-2); margin:16px 0 6px; font-weight:600; }
  .seg { display:flex; flex-direction:column; gap:4px; }
  .seg button { text-align:left; padding:6px 10px; border:1px solid var(--border); border-radius:6px; background:var(--surface); cursor:pointer; font:inherit; color:var(--ink); }
  .seg button.on { border-color:var(--accent); background:#eaf2fc; font-weight:600; }
  .seg button .u { color:var(--muted); font-weight:400; }
  label.line { display:flex; align-items:center; gap:8px; margin:4px 0; }
  select, input[type=number] { font:inherit; padding:3px 6px; border:1px solid var(--border); border-radius:5px; background:var(--surface); color:var(--ink); }
  input[type=number] { width:84px; }
  input[type=range] { width:100%; }
  .row2 { display:flex; gap:6px; align-items:center; }
  .link { background:none; border:0; color:var(--accent); cursor:pointer; padding:0; font:inherit; font-size:12px; }
  #readout { margin-top:6px; background:var(--surface-2); border-radius:6px; padding:8px 10px; min-height:96px; font-variant-numeric:tabular-nums; }
  #readout table { border-collapse:collapse; width:100%; }
  #readout td { padding:1px 0; } #readout td:first-child { color:var(--ink-2); padding-right:8px; white-space:nowrap; }
  main { display:grid; grid-template-columns:minmax(0,1fr) minmax(0,32%); grid-template-rows:minmax(0,56fr) minmax(0,44fr); min-width:0; }
  #c3d, #c2d, #ts { min-height:0; min-width:0; }
  #c3d { grid-column:1; grid-row:1; }
  #ts { grid-column:2; grid-row:1; border-left:1px solid var(--border); }
  #c2d { grid-column:1 / 3; grid-row:2; border-top:1px solid var(--border); }
</style></head>
<body><div id="app">
<aside>
  <div id="depRow" style="display:none;margin-bottom:8px"><select id="dep" style="width:100%"></select></div>
  <h1 id="title"></h1>
  <div class="muted" id="sub"></div>

  <h2>Variable</h2>
  <div class="seg" id="vars"></div>

  <h2>Colour</h2>
  <label class="line">Map <select id="cmap"></select></label>
  <label class="line"><input type="checkbox" id="log"> Log scale</label>
  <div class="row2"><input type="number" id="cmin" step="any"> to <input type="number" id="cmax" step="any"></div>
  <button class="link" id="auto">Reset to 1st–99th percentile</button>

  <h2>Date window</h2>
  <div class="muted" id="winLabel"></div>
  <input type="range" id="w0" min="0" value="0">
  <input type="range" id="w1" min="0">
  <button class="link" id="winAll">Whole deployment</button>

  <h2>Vertical exaggeration <span id="exLabel" class="muted"></span></h2>
  <input type="range" id="ex" min="0" max="100" value="50">
  <label class="line"><input type="checkbox" id="track" checked> Show surface track</label>
  <div id="bathyCtl">
    <label class="line"><input type="checkbox" id="bathy" checked> Seafloor &amp; land</label>
    <label class="line"><input type="checkbox" id="clip"> Limit depth axis to glider data</label>
    <div class="muted" id="bathyNote"></div>
  </div>

  <h2>Overlays</h2>
  <label class="line">Section depth <select id="secDepth"><option value="0">all</option><option value="500">0–500 m</option><option value="300">0–300 m</option><option value="150">0–150 m</option></select></label>
  <label class="line"><input type="checkbox" id="ovNight"> Shade night on section</label>
  <label class="line"><input type="checkbox" id="ovSig" checked> Density contours (σθ)</label>
  <label class="line"><input type="checkbox" id="ovMld" checked> Mixed-layer depth</label>
  <label class="line"><input type="checkbox" id="ovChl" checked> Chlorophyll-max depth</label>
  <div class="muted">Lines are 7-profile running medians; the readout gives each profile's own value.</div>
  <div class="muted">Sun strip above the section: <span style="color:#c9a227">■</span> day · <span style="color:#8f9aa6">■</span> civil twilight · <span style="color:#4a5868">■</span> nautical · <span style="color:#1d2733">■</span> night</div>

  <h2>T–S diagram</h2>
  <label class="line">Colour by <select id="tsBy"></select></label>
  <div class="muted">Points from the date window. Grey lines: σθ (EOS-80).</div>

  <h2>Under cursor</h2>
  <div id="readout" class="muted">Hover over the curtain or the section.</div>
</aside>
<main><div id="c3d"></div><div id="ts"></div><div id="c2d"></div></main>
</div>
<script id="data" type="application/json">__DATA__</script>
<script>
(async function () {
  let D;
  const emb = document.getElementById('data').textContent.trim();
  if (emb && emb !== 'null') {
    D = JSON.parse(emb);
  } else {                                   // site mode: shared viewer + data/<name>.json
    document.getElementById('title').textContent = 'Loading…';
    const idx = await fetch('data/index.json').then(r => r.json());
    const want = new URLSearchParams(location.search).get('d');
    const name = idx.some(e => e.name === want) ? want : idx[0].name;
    const sel = document.getElementById('dep');
    idx.forEach(e => sel.add(new Option(`${e.name}  (${e.start} → ${e.end})`, e.name, false, e.name === name)));
    sel.onchange = () => { location.search = '?d=' + encodeURIComponent(sel.value); };
    if (window.self === window.top) document.getElementById('depRow').style.display = '';   // hidden when embedded in a page
    D = await fetch(`data/${encodeURIComponent(name)}.json`).then(r => r.json());
  }
  const CMAPS = __CMAPS__;
  const f32 = b64 => { const s = atob(b64), u = new Uint8Array(s.length); for (let i = 0; i < s.length; i++) u[i] = s.charCodeAt(i); return new Float32Array(u.buffer); };
  const NP = D.nProf, ND = D.nDepth, X = f32(D.x), Y = f32(D.y);
  D.vars.forEach(v => { v.arr = f32(v.data); delete v.data; v.cmin = null; v.cmax = null; });
  const fmtT = s => new Date(s * 1000).toISOString().replace('T', ' ').slice(0, 16) + ' UTC';
  const fmtD = s => new Date(s * 1000).toISOString().slice(0, 10);
  const timesISO = D.time.map(s => new Date(s * 1000).toISOString());

  document.getElementById('title').textContent = D.name;
  document.getElementById('sub').textContent = `${fmtD(D.time[0])} → ${fmtD(D.time[NP - 1])} · ${NP.toLocaleString()} profiles · ${D.depth[ND - 1]} m max`;

  const st = { v: 0, i0: 0, i1: NP - 1, ex: 0, track: true, bathy: !!D.bathy, clip: false };
  const B = D.bathy || null;
  let BZ = null, BOTTOM = null, LAND_SCALE = 1;
  if (B) {
    BZ = f32(B.z); BOTTOM = f32(B.bottom);
    const dataMax = D.depth[ND - 1], span = Math.max(dataMax, B.maxDepth);
    LAND_SCALE = B.maxLand > 0 ? Math.min(1, 0.12 * span / B.maxLand) : 1;
    document.getElementById('bathyNote').innerHTML =
      `Seafloor to ${Math.round(B.maxDepth)} m in view` +
      (B.maxLand > 0 ? ` · land heights ×${LAND_SCALE < 1 ? LAND_SCALE.toFixed(2) : 1} (peak ${Math.round(B.maxLand)} m)` : '') +
      `<br>${B.source}${B.doi ? ' · ' + B.doi : ''}<br>${B.res}`;
  } else {
    document.getElementById('bathyCtl').innerHTML = '<div class="muted">No bathymetry grid found. Put a GEBCO/ETOPO .nc next to the glider files and rerun.</div>';
  }
  const zBottom = () => -(B && st.bathy && !st.clip ? Math.max(D.depth[ND - 1], B.maxDepth) : D.depth[ND - 1]) * 1.02;
  const zTop = () => (B && st.bathy && !st.clip && B.maxLand > 0) ? B.maxLand * LAND_SCALE * 1.05 : 0;
  const A = D.aux || {};
  const SIG = A.sigma ? f32(A.sigma) : null, TH = A.theta ? f32(A.theta) : null, SAL = A.sal ? f32(A.sal) : null;
  const SUN = A.sun || null, MLD = A.mld || null, CHLMAX = A.chlmax || null;
  Object.assign(st, { night: false, secDepth: 0, sig: !!SIG, mld: !!MLD, chl: !!CHLMAX, tsBy: 'depth' });
  const runMedian = (vals, w) => {   // vals may contain null
    const h = w >> 1, out = new Array(vals.length).fill(null);
    for (let j = 0; j < vals.length; j++) {
      if (vals[j] == null) continue;
      const win = [];
      for (let q = Math.max(0, j - h); q <= Math.min(vals.length - 1, j + h); q++) if (vals[q] != null) win.push(vals[q]);
      win.sort((a, b) => a - b); out[j] = win[win.length >> 1];
    }
    return out;
  };
  const sunClass = e => e > -0.833 ? 3 : e > -6 ? 2 : e > -12 ? 1 : 0;
  const sunName = ['night', 'nautical twilight', 'civil twilight', 'day'];
  const cur = () => D.vars[st.v];
  const scaleName = v => v.log ? `log₁₀ ${v.label}` : v.label;
  const tr = (v, x) => v.log ? (x > 0 ? Math.log10(x) : NaN) : x;
  const crange = v => {
    const r = v.log ? v.logRange : v.range;
    const lo = v.cmin ?? r[0], hi = v.cmax ?? r[1];
    return [tr(v, lo), tr(v, hi)];
  };
  const colorscale = v => CMAPS[v.cmap].map((c, i, a) => [i / (a.length - 1), c]);

  // ----- controls
  const varsEl = document.getElementById('vars');
  D.vars.forEach((v, i) => {
    const b = document.createElement('button');
    b.innerHTML = `${v.label} <span class="u">${v.units}</span>`;
    b.onclick = () => { st.v = i; syncControls(); draw(); };
    varsEl.appendChild(b);
  });
  const cmapEl = document.getElementById('cmap');
  Object.keys(CMAPS).forEach(k => cmapEl.add(new Option(k, k)));
  cmapEl.onchange = () => { cur().cmap = cmapEl.value; draw(); };
  const logEl = document.getElementById('log');
  logEl.onchange = () => { cur().log = logEl.checked; cur().cmin = cur().cmax = null; syncControls(); draw(); };
  const cminEl = document.getElementById('cmin'), cmaxEl = document.getElementById('cmax');
  cminEl.onchange = () => { cur().cmin = cminEl.value === '' ? null : +cminEl.value; draw(); };
  cmaxEl.onchange = () => { cur().cmax = cmaxEl.value === '' ? null : +cmaxEl.value; draw(); };
  document.getElementById('auto').onclick = () => { cur().cmin = cur().cmax = null; syncControls(); draw(); };

  const w0 = document.getElementById('w0'), w1 = document.getElementById('w1');
  w0.max = w1.max = NP - 1; w1.value = NP - 1;
  const onWin = () => {
    let a = +w0.value, b = +w1.value; if (a > b) [a, b] = [b, a];
    st.i0 = a; st.i1 = Math.max(b, a + 1); labelWin(); draw(true);
  };
  w0.oninput = w1.oninput = labelWinFromSliders;
  w0.onchange = w1.onchange = onWin;
  function labelWinFromSliders() { let a = +w0.value, b = +w1.value; if (a > b) [a, b] = [b, a]; document.getElementById('winLabel').textContent = `${fmtT(D.time[a])} → ${fmtT(D.time[b])}`; }
  function labelWin() { document.getElementById('winLabel').textContent = `${fmtT(D.time[st.i0])} → ${fmtT(D.time[st.i1])} (${st.i1 - st.i0 + 1} profiles)`; }
  document.getElementById('winAll').onclick = () => { w0.value = 0; w1.value = NP - 1; onWin(); };

  const exEl = document.getElementById('ex');
  const exVal = () => Math.round(Math.pow(10, 1 + 2.5 * exEl.value / 100));   // 10x .. ~3000x
  exEl.oninput = () => { document.getElementById('exLabel').textContent = `(${exVal()}×)`; };
  exEl.onchange = () => { exEl.oninput(); relayoutAspect(); };
  document.getElementById('track').onchange = e => { st.track = e.target.checked; draw(true); };
  if (B) {
    document.getElementById('bathy').onchange = e => { st.bathy = e.target.checked; draw(true); };
    document.getElementById('clip').onchange = e => { st.clip = e.target.checked; draw(true); };
  }

  [['ovNight', 'night', !!SUN], ['ovSig', 'sig', !!SIG], ['ovMld', 'mld', !!MLD], ['ovChl', 'chl', !!CHLMAX]].forEach(([id, key, ok]) => {
    const el = document.getElementById(id);
    if (!ok) { el.parentElement.style.display = 'none'; return; }
    el.onchange = () => { st[key] = el.checked; draw(true); };
  });
  document.getElementById('secDepth').onchange = e => { st.secDepth = +e.target.value; draw(true); };
  const tsByEl = document.getElementById('tsBy');
  tsByEl.add(new Option('Depth', 'depth'));
  D.vars.forEach((v, i) => tsByEl.add(new Option(v.label, String(i))));
  tsByEl.onchange = () => { st.tsBy = tsByEl.value; drawTS(); };

  function syncControls() {
    [...varsEl.children].forEach((b, i) => b.classList.toggle('on', i === st.v));
    const v = cur(); cmapEl.value = v.cmap; logEl.checked = v.log;
    const r = v.log ? v.logRange : v.range;
    cminEl.value = v.cmin ?? r[0]; cmaxEl.value = v.cmax ?? r[1];
  }

  // ----- geometry for the window
  function windowGrid() {
    const i0 = st.i0, i1 = st.i1, n = i1 - i0 + 1;
    // Downsample columns if very large, to keep WebGL responsive.
    const step = Math.max(1, Math.ceil(n * ND / 600000));
    const cols = []; for (let i = i0; i <= i1; i += step) cols.push(i);
    return cols;
  }

  // The 3D frame always spans the whole deployment, so a date window is seen in place.
  let FX0 = Infinity, FX1 = -Infinity, FY0 = Infinity, FY1 = -Infinity;
  for (let i = 0; i < NP; i++) { FX0 = Math.min(FX0, X[i]); FX1 = Math.max(FX1, X[i]); FY0 = Math.min(FY0, Y[i]); FY1 = Math.max(FY1, Y[i]); }
  const padX = Math.max(1, 0.03 * (FX1 - FX0)), padY = Math.max(1, 0.03 * (FY1 - FY0));
  FX0 -= padX; FX1 += padX; FY0 -= padY; FY1 += padY;
  function aspect() {
    const sx = FX1 - FX0, sy = FY1 - FY0, m = Math.max(sx, sy);
    const zkm = (zTop() - zBottom()) / 1000;
    if (!st.exUser) {               // first draw: choose exaggeration so depth ~ 40% of width
      const want = 0.4 * m / zkm;
      exEl.value = Math.max(0, Math.min(100, (Math.log10(want) - 1) / 2.5 * 100));
      exEl.oninput(); st.exUser = true;
    }
    return { x: sx / m, y: sy / m, z: Math.min(3, zkm * exVal() / m) };
  }

  function relayoutAspect() {
    Plotly.relayout('c3d', { 'scene.aspectratio': aspect() });
  }

  let lastCols = null;
  function draw(keepCamera) {
    const v = cur(), cols = windowGrid(); lastCols = cols;
    const [cmin, cmax] = crange(v);
    const nc = cols.length;
    const xs = [], ys = [], zs = [], cs = [];
    for (let k = 0; k < ND; k++) {
      const xr = new Array(nc), yr = new Array(nc), zr = new Array(nc), cr = new Array(nc);
      for (let j = 0; j < nc; j++) {
        const i = cols[j], val = v.arr[k * NP + i], c = tr(v, val);
        xr[j] = X[i]; yr[j] = Y[i]; cr[j] = c;
        zr[j] = Number.isFinite(c) ? -D.depth[k] : NaN;   // NaN z = hole in the curtain
      }
      xs.push(xr); ys.push(yr); zs.push(zr); cs.push(cr);
    }
    const cbTitle = `${scaleName(v)}${v.units ? ' (' + v.units + ')' : ''}`;
    const surf = {
      type: 'surface', x: xs, y: ys, z: zs, surfacecolor: cs,
      colorscale: colorscale(v), cmin, cmax, showscale: true,
      colorbar: { title: { text: cbTitle, side: 'right' }, thickness: 14, len: 0.8 },
      hoverinfo: 'none', lighting: { ambient: 1, diffuse: 0, specular: 0, roughness: 1, fresnel: 0 },
      contours: { x: { highlight: false }, y: { highlight: false }, z: { highlight: false } },
    };
    const traces = [surf];
    if (st.track) {
      const all = []; for (let i = 0; i < NP; i += Math.max(1, Math.floor(NP / 4000))) all.push(i);
      traces.push({ type: 'scatter3d', mode: 'lines', x: all.map(i => X[i]), y: all.map(i => Y[i]), z: all.map(() => 0),
                    line: { color: '#b9b7b1', width: 2 }, hoverinfo: 'skip', showlegend: false });
      traces.push({ type: 'scatter3d', mode: 'lines', x: cols.map(i => X[i]), y: cols.map(i => Y[i]), z: cols.map(() => 0),
                    line: { color: '#52514e', width: 3 }, hoverinfo: 'skip', showlegend: false });
      traces.push({ type: 'scatter3d', mode: 'markers', x: [X[cols[0]]], y: [Y[cols[0]]], z: [0],
                    marker: { size: 4, color: '#0b0b0b' }, hovertext: ['start'], hoverinfo: 'text', showlegend: false });
    }
    const line3d = (arr, color, width, name) => {
      const xs3 = [], ys3 = [], zs3 = [];
      const sm = runMedian(cols.map(i => arr[i]), 7);
      cols.forEach((i, j) => { xs3.push(X[i]); ys3.push(Y[i]); zs3.push(sm[j] == null ? null : -sm[j]); });
      traces.push({ type: 'scatter3d', mode: 'lines', x: xs3, y: ys3, z: zs3, connectgaps: false,
                    line: { color, width }, hoverinfo: 'skip', showlegend: false, name });
    };
    if (st.mld && MLD) line3d(MLD, '#ffffff', 5, 'MLD');
    if (st.chl && CHLMAX) line3d(CHLMAX, '#eb6834', 4, 'Chl max');
    if (B && st.bathy) {
      const zb = zBottom(), z = [], c = [];
      for (let r = 0; r < B.ny; r++) {
        const zr = new Array(B.nx), cr = new Array(B.nx);
        for (let q = 0; q < B.nx; q++) {
          const e = BZ[r * B.nx + q];
          cr[q] = e;
          zr[q] = !Number.isFinite(e) ? NaN : e > 0 ? e * LAND_SCALE : (e < zb ? NaN : e);
        }
        z.push(zr); c.push(cr);
      }
      const lo = -Math.max(B.maxDepth, 1), hi = Math.max(B.maxLand, 1), f0 = -lo / (hi - lo);
      traces.push({ type: 'surface', x: B.x, y: B.y, z, surfacecolor: c, cmin: lo, cmax: hi, showscale: false,
        colorscale: [[0, '#3f4a55'], [Math.max(0, f0 - 1e-4), '#c9d1d8'], [f0, '#d8ccad'], [1, '#7d6e4f']],
        opacity: 1, hoverinfo: 'none',
        lighting: { ambient: 0.55, diffuse: 0.6, specular: 0.05, roughness: 0.9 }, lightposition: { x: 1e5, y: 1e5, z: 2e5 },
        contours: { x: { highlight: false }, y: { highlight: false }, z: { highlight: false } } });
      if (B.coast) traces.push({ type: 'scatter3d', mode: 'lines', x: B.coast.x, y: B.coast.y, z: B.coast.x.map(v => v === null ? null : 0),
        line: { color: '#2b2b29', width: 2 }, hoverinfo: 'skip', showlegend: false, connectgaps: false });
    }
    const gd = document.getElementById('c3d');
    const layout = {
      margin: { l: 0, r: 0, t: 0, b: 0 }, paper_bgcolor: '#fcfcfb', showlegend: false,
      uirevision: 'keep',
      scene: {
        xaxis: { title: { text: 'km east' }, range: [FX0, FX1], showspikes: false, gridcolor: '#e2e0db' },
        yaxis: { title: { text: 'km north' }, range: [FY0, FY1], showspikes: false, gridcolor: '#e2e0db' },
        zaxis: { title: { text: 'depth (m)' }, range: [zBottom(), zTop()], showspikes: false, gridcolor: '#e2e0db',
                 tickvals: niceDepths().map(d => -d), ticktext: niceDepths().map(String) },
        aspectmode: 'manual', aspectratio: aspect(),
        camera: keepCamera && gd.layout ? gd.layout.scene.camera : { eye: { x: 1.15, y: -1.3, z: 0.7 } },
      },
    };
    Plotly.react(gd, traces, layout, { responsive: true, displaylogo: false });
    draw2d(v, cmin, cmax, cbTitle);
  }

  function niceDepths() { const m = -zBottom(), s = m > 4000 ? 1000 : m > 1500 ? 500 : m > 600 ? 200 : m > 250 ? 100 : 50, out = []; for (let d = 0; d <= m; d += s) out.push(d); return out; }

  function draw2d(v, cmin, cmax, cbTitle) {
    // full deployment; the 3D window is shaded
    const step = Math.max(1, Math.ceil(NP / 3000));
    const cols = []; for (let i = 0; i < NP; i += step) cols.push(i);
    const z = [];
    for (let k = 0; k < ND; k++) { const r = new Array(cols.length); for (let j = 0; j < cols.length; j++) { const c = tr(v, v.arr[k * NP + cols[j]]); r[j] = Number.isFinite(c) ? c : null; } z.push(r); }
    const t = cols.map(i => timesISO[i]);
    const trace = { type: 'heatmap', x: t, y: D.depth, z, colorscale: colorscale(v), zmin: cmin, zmax: cmax,
                    showscale: false, hoverinfo: 'none', connectgaps: false };
    const shapes = (st.i0 > 0 || st.i1 < NP - 1) ? [
      { type: 'rect', xref: 'x', yref: 'paper', x0: timesISO[0], x1: timesISO[st.i0], y0: 0, y1: 1, fillcolor: '#fcfcfb', opacity: 0.5, line: { width: 0 } },
      { type: 'rect', xref: 'x', yref: 'paper', x0: timesISO[st.i1], x1: timesISO[NP - 1], y0: 0, y1: 1, fillcolor: '#fcfcfb', opacity: 0.5, line: { width: 0 } },
    ] : [];
    const yMax = st.secDepth ? Math.min(st.secDepth, D.depth[ND - 1] * 1.03) : D.depth[ND - 1] * 1.03, traces2 = [trace];
    if (B) {
      const xs = [], ys = [];
      cols.forEach(i => { xs.push(timesISO[i]); ys.push(Number.isFinite(BOTTOM[i]) ? Math.min(BOTTOM[i], yMax) : yMax); });
      traces2.push({ type: 'scatter', mode: 'lines', x: xs.concat([xs[xs.length - 1], xs[0]]), y: ys.concat([yMax, yMax]),
        fill: 'toself', fillcolor: 'rgba(63,74,85,0.85)', line: { color: '#2b2b29', width: 1 }, hoverinfo: 'skip', showlegend: false });
    }
    if (st.sig && A.sigLines) {
      const L = A.sigLines;
      traces2.push({ type: 'scattergl', mode: 'lines', x: L.t.map(v => v == null ? null : new Date(v * 1000).toISOString()), y: L.d,
        line: { color: 'rgba(11,11,11,0.5)', width: 0.8 }, connectgaps: false, hoverinfo: 'skip', name: 'σθ contours', showlegend: true });
      traces2.push({ type: 'scatter', mode: 'text', x: L.labels.map(l => new Date(l.t * 1000).toISOString()), y: L.labels.map(l => l.d),
        text: L.labels.map(l => l.txt), textfont: { size: 9, color: '#0b0b0b' }, hoverinfo: 'skip', showlegend: false, cliponaxis: true });
    } else if (st.sig && SIG) {
      const zs = [];
      for (let k = 0; k < ND; k++) {      // 5-column running mean along time, display only
        const raw = cols.map(i => SIG[k * NP + i]), r = new Array(cols.length);
        for (let j = 0; j < cols.length; j++) {
          if (!Number.isFinite(raw[j])) { r[j] = null; continue; }
          let sum = 0, n = 0; for (let q = Math.max(0, j - 2); q <= Math.min(cols.length - 1, j + 2); q++) if (Number.isFinite(raw[q])) { sum += raw[q]; n++; }
          r[j] = sum / n;
        }
        zs.push(r);
      }
      traces2.push({ type: 'contour', x: t, y: D.depth, z: zs, showscale: false, hoverinfo: 'skip', connectgaps: false,
        contours: { coloring: 'lines', start: A.sigStart, end: A.sigEnd, size: A.sigStep, showlabels: true,
                    labelfont: { size: 9, color: '#0b0b0b' } },
        line: { color: 'rgba(11,11,11,0.55)', width: 0.7, smoothing: 0.8 }, colorscale: [[0, '#0b0b0b'], [1, '#0b0b0b']], name: 'σθ' });
    }
    if (st.mld && MLD) {
      const y = runMedian(cols.map(i => MLD[i]), 7);
      traces2.push({ type: 'scatter', mode: 'lines', x: t, y, line: { color: '#0b0b0b', width: 3 }, connectgaps: false, hoverinfo: 'skip', showlegend: false });
      traces2.push({ type: 'scatter', mode: 'lines', x: t, y, line: { color: '#ffffff', width: 1.4 }, connectgaps: false, hoverinfo: 'skip',
                     name: 'Mixed-layer depth (Δσθ = 0.03)', showlegend: true });
    }
    if (st.chl && CHLMAX) {
      const y = runMedian(cols.map(i => CHLMAX[i]), 7);
      traces2.push({ type: 'scatter', mode: 'lines', x: t, y, line: { color: '#0b0b0b', width: 3 }, connectgaps: false, hoverinfo: 'skip', showlegend: false });
      traces2.push({ type: 'scatter', mode: 'lines', x: t, y, line: { color: '#eb6834', width: 1.6 }, connectgaps: false, hoverinfo: 'skip',
                     name: 'Chlorophyll-max depth', showlegend: true });
    }
    if (SUN) {
      traces2.push({ type: 'heatmap', x: t, y: [0], z: [cols.map(i => sunClass(SUN[i]))], yaxis: 'y2', zmin: 0, zmax: 3, showscale: false,
        colorscale: [[0, '#1d2733'], [0.25, '#1d2733'], [0.25, '#4a5868'], [0.5, '#4a5868'], [0.5, '#8f9aa6'], [0.75, '#8f9aa6'], [0.75, '#f2d77c'], [1, '#f2d77c']],
        hoverinfo: 'none' });
      if (st.night) {   // shade runs with the sun below the horizon
        let run = null;
        cols.forEach((i, j) => {
          const dark = SUN[i] <= -0.833;
          if (dark && run === null) run = j;
          if ((!dark || j === cols.length - 1) && run !== null) {
            shapes.push({ type: 'rect', xref: 'x', yref: 'y', x0: t[run], x1: t[j], y0: 0, y1: yMax, fillcolor: 'rgba(29,39,51,0.22)', line: { width: 0 }, layer: 'above' });
            run = null;
          }
        });
      }
    }
    Plotly.react('c2d', traces2, {
      margin: { l: 56, r: 20, t: 26, b: 36 }, paper_bgcolor: '#fcfcfb', plot_bgcolor: '#f3f2ef',
      xaxis: { type: 'date', showgrid: false }, yaxis: { range: [yMax, 0], title: { text: 'depth (m)' }, showgrid: false, domain: [0, SUN ? 0.93 : 1] },
      yaxis2: { domain: [0.95, 1], showticklabels: true, tickvals: [0], ticktext: ['sun'], fixedrange: true, showgrid: false, zeroline: false },
      legend: { orientation: 'h', x: 1, xanchor: 'right', y: 1.0, yanchor: 'bottom', bgcolor: 'rgba(0,0,0,0)', font: { size: 11 } },
      shapes, uirevision: 'keep2', dragmode: 'zoom',
    }, { responsive: true, displaylogo: false });
    drawTS();
  }

  // ----- T-S diagram
  function drawTS() {
    const el = document.getElementById('ts');
    if (!TH || !SAL || !A.ts) { el.innerHTML = '<div class="muted" style="padding:16px">No temperature/salinity for a T–S diagram.</div>'; return; }
    const T = A.ts, cols = [];
    const n = st.i1 - st.i0 + 1, step = Math.max(1, Math.ceil(n * ND / 60000));
    for (let i = st.i0; i <= st.i1; i += step) cols.push(i);
    const xs = [], ys = [], cs = [], cd = [];
    const byVar = st.tsBy !== 'depth' ? D.vars[+st.tsBy] : null;
    for (const i of cols) for (let k = 0; k < ND; k++) {
      const q = k * NP + i, s0 = SAL[q], t0 = TH[q];
      if (!Number.isFinite(s0) || !Number.isFinite(t0)) continue;
      let c = D.depth[k];
      if (byVar) { c = tr(byVar, byVar.arr[q]); if (!Number.isFinite(c)) continue; }
      xs.push(s0); ys.push(t0); cs.push(c); cd.push([i, k]);
    }
    let cscale, cmin, cmax, ctitle;
    if (byVar) { cscale = colorscale(byVar); [cmin, cmax] = crange(byVar); ctitle = `${scaleName(byVar)}${byVar.units ? ' (' + byVar.units + ')' : ''}`; }
    else { cscale = CMAPS.deep.map((c, i, a) => [i / (a.length - 1), c]); cmin = 0; cmax = D.depth[ND - 1]; ctitle = 'depth (m)'; }
    const traces = [
      { type: 'contour', x: T.s, y: T.t, z: T.z, showscale: false, hoverinfo: 'skip',
        contours: { coloring: 'lines', start: T.start, end: T.end, size: T.step, showlabels: true, labelfont: { size: 9, color: '#7a7873' } },
        line: { color: '#b9b7b1', width: 0.8 }, colorscale: [[0, '#b9b7b1'], [1, '#b9b7b1']] },
      { type: 'scattergl', mode: 'markers', x: xs, y: ys, customdata: cd, hoverinfo: 'none',
        marker: { size: 3, color: cs, colorscale: cscale, cmin, cmax, opacity: 0.75,
                  colorbar: { title: { text: ctitle, side: 'right' }, thickness: 10, len: 0.8 } } },
    ];
    Plotly.react(el, traces, {
      margin: { l: 52, r: 8, t: 10, b: 40 }, paper_bgcolor: '#fcfcfb', plot_bgcolor: '#fcfcfb', showlegend: false,
      xaxis: { title: { text: 'Salinity (PSU)' }, range: T.srange, showgrid: false, zeroline: false },
      yaxis: { title: { text: `${T.tname} (°C)` }, range: T.trange, showgrid: false, zeroline: false },
      uirevision: 'ts',
    }, { responsive: true, displaylogo: false });
  }

  // ----- readout
  const ro = document.getElementById('readout');
  function readout(i, k) {
    if (i == null || k == null || i < 0 || k < 0) return;
    let rows = `<tr><td>Time</td><td>${fmtT(D.time[i])}</td></tr>
      <tr><td>Lat, Lon</td><td>${D.lat[i].toFixed(4)}, ${D.lon[i].toFixed(4)}</td></tr>
      <tr><td>Depth bin</td><td>${D.depth[k]} m</td></tr>` +
      (B ? `<tr><td>Seafloor</td><td>${Number.isFinite(BOTTOM[i]) ? (BOTTOM[i] > 0 ? Math.round(BOTTOM[i]) + ' m' : 'land (' + Math.round(-BOTTOM[i]) + ' m elev.)') : '—'}</td></tr>` : '');
    const q = k * NP + i;
    if (TH) rows += `<tr><td>${A.ts ? A.ts.tname : 'θ'}</td><td>${Number.isFinite(TH[q]) ? TH[q].toFixed(3) + ' °C' : '—'}</td></tr>`;
    if (SIG) rows += `<tr><td>σθ</td><td>${Number.isFinite(SIG[q]) ? SIG[q].toFixed(3) + ' kg m⁻³' : '—'}</td></tr>`;
    if (MLD) rows += `<tr><td>Mixed layer</td><td>${MLD[i] != null ? MLD[i] + ' m' : '—'}</td></tr>`;
    if (CHLMAX) rows += `<tr><td>Chl max</td><td>${CHLMAX[i] != null ? CHLMAX[i] + ' m (' + A.chlmaxVal[i] + ')' : '—'}</td></tr>`;
    if (SUN) rows += `<tr><td>Sun</td><td>${SUN[i].toFixed(1)}° · ${sunName[sunClass(SUN[i])]}</td></tr>`;
    D.vars.forEach((v, j) => {
      const val = v.arr[k * NP + i];
      rows += `<tr><td${j === st.v ? ' style="color:var(--ink);font-weight:600"' : ''}>${v.label}</td><td>${Number.isFinite(val) ? val.toPrecision(4) + ' ' + v.units : '—'}</td></tr>`;
    });
    ro.classList.remove('muted'); ro.innerHTML = `<table>${rows}</table>`;
  }
  const g3 = document.getElementById('c3d'), g2 = document.getElementById('c2d');
  const depthIdx = d => { let b = 0, bd = Infinity; D.depth.forEach((x, k) => { const q = Math.abs(x - d); if (q < bd) { bd = q; b = k; } }); return b; };
  function bind() {
    g3.on('plotly_hover', e => {
      const p = e.points && e.points[0]; if (!p || p.curveNumber !== 0) return;
      // Locate the cell from the hovered coordinates (robust to plotly's index ordering).
      const k = depthIdx(-p.z);
      let i = lastCols[0], bd = Infinity;
      for (const c of lastCols) { const q = (X[c] - p.x) ** 2 + (Y[c] - p.y) ** 2; if (q < bd) { bd = q; i = c; } }
      readout(i, k);
    });
    document.getElementById('ts').on && document.getElementById('ts').on('plotly_hover', e => {
      const p = e.points && e.points.find(pt => pt.customdata); if (p) readout(p.customdata[0], p.customdata[1]);
    });
    g2.on('plotly_hover', e => {
      const p = e.points && e.points.find(pt => pt.curveNumber === 0); if (!p) return;
      const t = new Date(p.x).getTime() / 1000;
      let lo = 0, hi = NP - 1; while (lo < hi) { const m = (lo + hi) >> 1; if (D.time[m] < t) lo = m + 1; else hi = m; }
      readout(lo, depthIdx(p.y));
    });
    g2.on('plotly_click', e => {   // click the section to centre a 2-day window there
      const p = e.points && e.points.find(pt => pt.curveNumber === 0); if (!p) return;
      const t = new Date(p.x).getTime() / 1000;
      let a = 0; while (a < NP - 1 && D.time[a] < t - 86400) a++;
      let b = a; while (b < NP - 1 && D.time[b] < t + 86400) b++;
      w0.value = a; w1.value = b; onWin();
    });
  }

  syncControls(); labelWin(); draw(false); bind();
})();
</script></body></html>
"""

if __name__ == "__main__":
    main()
