#!/usr/bin/env python3
"""
glider_tracks.py - Extract time/lat/lon from IOOS Glider DAC NetCDF files and
build an interactive, filterable map of the deployments.

Usage
-----
    python glider_tracks.py PATH [PATH ...] [options]

PATH can be a folder (searched recursively for *.nc) or individual .nc files.

Outputs (written to --out-dir, default ./glider_output):
    glider_positions.csv    every valid position: deployment, time (epoch s + ISO UTC), lat, lon
    deployments_summary.csv one row per deployment: start, end, points, files, bounds
    glider_map.html         interactive map; filter by year and/or deployment name.
                            Open it in any web browser (needs internet for the
                            Leaflet library and basemap tiles).

Requirements:  pip install netCDF4 numpy

Handles
-------
* IOOS Glider DAC profile files (one profile per file; variables time/lat/lon,
  or profile_time/profile_lat/profile_lon with --positions profile).
* ERDDAP-style aggregated trajectory files (time/latitude/longitude, optionally
  several trajectories in one file as a contiguous ragged array with rowSize).
* Time in any CF "<units> since <date>" encoding; converted to seconds since
  1970-01-01 UTC. Values already in seconds since 1970 are passed through.

Deployment name comes from (in order, with --group-by auto):
  the 'trajectory' variable (cf_role = trajectory_id), the global attributes
  'deployment_name' / 'deployment_id', then the name of the file's parent folder.
"""

from __future__ import annotations

import argparse
import csv
import datetime as dt
import json
import math
import os
import re
import sys
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np

try:
    import netCDF4
except ImportError:  # pragma: no cover
    sys.exit("netCDF4 is required:  pip install netCDF4 numpy")


EPOCH_UNITS = "seconds since 1970-01-01 00:00:00"
EPOCH_RE = re.compile(r"^seconds since 1970-0?1-0?1([ T]0?0:00(:00(\.0+)?)?)?\s*(z|utc|gmt|\+0?0:?00)?$",
                      re.IGNORECASE)

# Preferred variable names, in order. Fallback is a search by standard_name.
NAMES = {
    "measurement": {
        "time": ["time", "m_present_time"],
        "lat": ["lat", "latitude"],
        "lon": ["lon", "longitude"],
    },
    "profile": {
        "time": ["profile_time", "time"],
        "lat": ["profile_lat", "lat", "latitude"],
        "lon": ["profile_lon", "lon", "longitude"],
    },
}
STD_NAMES = {"time": "time", "lat": "latitude", "lon": "longitude"}


# --------------------------------------------------------------------------- #
# Reading
# --------------------------------------------------------------------------- #
def _find_var(ds, role, mode):
    for name in NAMES[mode][role]:
        if name in ds.variables:
            return ds.variables[name]
    # Fallback: first variable with the right standard_name, skipping
    # depth-averaged-current (_uv) and the other mode's profile_* variables.
    for name, var in ds.variables.items():
        if getattr(var, "standard_name", None) != STD_NAMES[role]:
            continue
        if name.endswith("_uv") or name.endswith("_qc"):
            continue
        if mode == "measurement" and name.startswith("profile_"):
            continue
        return var
    return None


def _as_float(var):
    arr = var[:]
    if np.ma.isMaskedArray(arr):
        return np.ma.filled(arr.astype("f8"), np.nan).ravel()
    return np.asarray(arr, dtype="f8").ravel()


def _to_epoch_seconds(var):
    vals = _as_float(var)
    units = (getattr(var, "units", "") or "").strip()
    if EPOCH_RE.match(units):
        return vals
    if "since" not in units.lower():
        # No CF units: assume the caller's convention (seconds since 1970).
        return vals
    out = np.full(vals.shape, np.nan)
    good = np.isfinite(vals)
    if good.any():
        cal = getattr(var, "calendar", "standard")
        dates = netCDF4.num2date(vals[good], units, calendar=cal,
                                 only_use_cftime_datetimes=False,
                                 only_use_python_datetimes=True)
        out[good] = netCDF4.date2num(dates, EPOCH_UNITS, calendar="standard")
    return out


def _string_values(var):
    """Return a list of strings from a char/str variable (scalar or 1-D)."""
    data = var[:]
    if var.dtype == str or data.dtype.kind in ("U", "O"):
        vals = np.atleast_1d(data).tolist()
    else:  # char array -> join along last dim
        if np.ma.isMaskedArray(data):
            data = data.filled(b"")
        vals = np.atleast_1d(netCDF4.chartostring(data)).tolist()
    return [str(v).strip() for v in vals]


def _deployment_names(ds, path, group_by):
    folder_name = Path(path).parent.name or Path(path).stem
    if group_by == "folder":
        return [folder_name], None
    if group_by == "file":
        return [Path(path).stem], None

    # Text-valued trajectory ID only: some files (e.g. PyGlider grids) tag an
    # integer such as mission_number with cf_role=trajectory_id.
    def is_text(v):
        return v.dtype == str or np.dtype(v.dtype).kind in ("S", "U", "O")
    traj = None
    cands = ([ds.variables["trajectory"]] if "trajectory" in ds.variables else []) + \
            [v for v in ds.variables.values() if getattr(v, "cf_role", None) == "trajectory_id"]
    for v in cands:
        if is_text(v):
            traj = v
            break
    if traj is not None:
        names = [n for n in _string_values(traj)]
        if any(names):
            row_size = ds.variables["rowSize"][:] if "rowSize" in ds.variables else None
            return [n or folder_name for n in names], row_size

    for attr in ("deployment_name", "deployment_id"):
        val = getattr(ds, attr, None)
        if val:
            return [str(val).strip()], None
    return [folder_name], None


def read_file(args):
    """Read one file -> (path, [(name, t, lat, lon)], error_or_None)."""
    path, mode, group_by = args
    try:
        with netCDF4.Dataset(path) as ds:
            ds.set_auto_mask(True)
            tv, latv, lonv = (_find_var(ds, r, mode) for r in ("time", "lat", "lon"))
            missing = [r for r, v in zip(("time", "lat", "lon"), (tv, latv, lonv)) if v is None]
            if missing:
                return path, [], f"missing variable(s): {', '.join(missing)}"
            t = _to_epoch_seconds(tv)
            lat = _as_float(latv)
            lon = _as_float(lonv)
            if not (t.size == lat.size == lon.size):
                # e.g. scalar profile_time with per-measurement lat: broadcast scalars
                n = max(t.size, lat.size, lon.size)
                try:
                    t, lat, lon = (np.broadcast_to(a, (n,)) if a.size == 1 else a
                                   for a in (t, lat, lon))
                    assert t.size == lat.size == lon.size
                except Exception:
                    return path, [], (f"time/lat/lon sizes differ "
                                      f"({t.size}, {lat.size}, {lon.size})")
            names, row_size = _deployment_names(ds, path, group_by)
    except Exception as exc:  # unreadable / not netCDF
        return path, [], f"{type(exc).__name__}: {exc}"

    if len(names) > 1 and row_size is not None and len(row_size) == len(names):
        bounds = np.concatenate([[0], np.cumsum(np.asarray(row_size, dtype=int))])
        chunks = [(names[i], slice(bounds[i], bounds[i + 1])) for i in range(len(names))]
    else:
        chunks = [(names[0], slice(None))]

    out = []
    for name, sl in chunks:
        tt, la, lo = t[sl], lat[sl], lon[sl]
        ok = (np.isfinite(tt) & np.isfinite(la) & np.isfinite(lo)
              & (np.abs(la) <= 90) & (np.abs(lo) <= 180)
              & ~((la == 0) & (lo == 0)))          # 0,0 = failed GPS fix
        if ok.any():
            out.append((name, tt[ok], la[ok], lo[ok]))
    return path, out, None


# --------------------------------------------------------------------------- #
# Output
# --------------------------------------------------------------------------- #
def iso(t):
    return dt.datetime.fromtimestamp(float(t), dt.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def decimate(n, max_points):
    if max_points <= 0 or n <= max_points:
        return np.arange(n)
    return np.unique(np.linspace(0, n - 1, max_points).round().astype(int))


def build(paths, out_dir, mode, group_by, max_map_points, workers, title, link_template=None, map_name="glider_map.html"):
    files = []
    for p in paths:
        p = Path(p)
        if p.is_dir():
            files += sorted(q for q in p.rglob("*") if q.suffix.lower() in (".nc", ".nc4"))
        elif p.is_file():
            files.append(p)
        else:
            print(f"warning: {p} not found", file=sys.stderr)
    if not files:
        sys.exit("No .nc files found.")
    print(f"Reading {len(files)} files ...")

    acc = defaultdict(lambda: {"t": [], "lat": [], "lon": [], "files": 0})
    errors = []
    jobs = [(str(f), mode, group_by) for f in files]
    if workers > 1:
        with ProcessPoolExecutor(max_workers=workers) as ex:
            results = ex.map(read_file, jobs, chunksize=32)
            results = list(_progress(results, len(jobs)))
    else:
        results = _progress(map(read_file, jobs), len(jobs))
    for path, chunks, err in results:
        if err:
            errors.append((path, err))
            continue
        for name, t, la, lo in chunks:
            d = acc[name]
            d["t"].append(t); d["lat"].append(la); d["lon"].append(lo); d["files"] += 1

    if not acc:
        for p, e in errors[:20]:
            print(f"  {p}: {e}", file=sys.stderr)
        sys.exit("No valid positions found in any file.")

    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    deployments = []
    pos_path = out_dir / "glider_positions.csv"
    with open(pos_path, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["deployment", "time_epoch_s", "time_utc", "lat", "lon"])
        for name in sorted(acc):
            d = acc[name]
            t = np.concatenate(d["t"]); la = np.concatenate(d["lat"]); lo = np.concatenate(d["lon"])
            order = np.argsort(t, kind="stable")
            t, la, lo = t[order], la[order], lo[order]
            keep = np.ones(t.size, bool)                 # drop exact duplicates
            keep[1:] = (np.diff(t) != 0) | (np.diff(la) != 0) | (np.diff(lo) != 0)
            t, la, lo = t[keep], la[keep], lo[keep]
            for a, b, c in zip(t, la, lo):
                w.writerow([name, f"{a:.3f}", iso(a), f"{b:.6f}", f"{c:.6f}"])
            idx = decimate(t.size, max_map_points)
            deployments.append({
                "name": name,
                "start": float(t[0]), "end": float(t[-1]),
                "n": int(t.size), "files": d["files"],
                "bounds": [float(la.min()), float(lo.min()), float(la.max()), float(lo.max())],
                "t": [int(round(x)) for x in t[idx]],
                "lat": [round(float(x), 5) for x in la[idx]],
                "lon": [round(float(x), 5) for x in lo[idx]],
            })

    sum_path = out_dir / "deployments_summary.csv"
    with open(sum_path, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["deployment", "start_utc", "end_utc", "duration_days", "n_points", "n_files",
                    "lat_min", "lon_min", "lat_max", "lon_max"])
        for d in deployments:
            w.writerow([d["name"], iso(d["start"]), iso(d["end"]),
                        f"{(d['end'] - d['start']) / 86400:.1f}", d["n"], d["files"],
                        *[f"{b:.5f}" for b in d["bounds"]]])

    map_path = out_dir / map_name
    if link_template:
        for d in deployments:
            d["link"] = link_template.replace("{name}", d["name"])
    payload = json.dumps({"title": title, "deployments": deployments}, separators=(",", ":"))
    map_path.write_text(HTML_TEMPLATE.replace("__TITLE__", _esc(title))
                        .replace("__DATA__", payload.replace("</", "<\\/")), encoding="utf-8")

    print(f"\n{len(deployments)} deployments, {sum(d['n'] for d in deployments):,} positions")
    if errors:
        err_path = out_dir / "skipped_files.txt"
        err_path.write_text("\n".join(f"{p}\t{e}" for p, e in errors), encoding="utf-8")
        print(f"{len(errors)} files skipped - reasons in {err_path}")
    print(f"Wrote:\n  {pos_path}\n  {sum_path}\n  {map_path}")


def _progress(it, total):
    step = max(1, total // 20)
    for i, item in enumerate(it, 1):
        if i % step == 0 or i == total:
            print(f"  {i}/{total}", end="\r", flush=True)
        yield item
    print()


def _esc(s):
    return s.replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;")


# --------------------------------------------------------------------------- #
# Map page
# --------------------------------------------------------------------------- #
HTML_TEMPLATE = r"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>__TITLE__</title>
<link rel="stylesheet" href="https://unpkg.com/leaflet@1.9.4/dist/leaflet.css">
<script src="https://unpkg.com/leaflet@1.9.4/dist/leaflet.js"></script>
<style>
  :root {
    --surface-1:#fcfcfb; --surface-2:#f3f2ef; --border:#e2e0db;
    --text-primary:#0b0b0b; --text-secondary:#52514e; --text-muted:#7a7873;
    --accent:#2a78d6; --highlight:#eb6834;
  }
  @media (prefers-color-scheme: dark) {
    :root { --surface-1:#1a1a19; --surface-2:#242422; --border:#383835;
            --text-primary:#f0efec; --text-secondary:#c3c2b7; --text-muted:#9a988f;
            --accent:#3987e5; --highlight:#d95926; }
  }
  * { box-sizing:border-box; }
  html,body { margin:0; height:100%; font:14px/1.4 system-ui,-apple-system,"Segoe UI",sans-serif;
              color:var(--text-primary); background:var(--surface-1); }
  #app { display:flex; height:100%; }
  #panel { width:360px; min-width:260px; display:flex; flex-direction:column; border-right:1px solid var(--border); }
  #panel header { padding:14px 16px 10px; border-bottom:1px solid var(--border); }
  h1 { font-size:16px; margin:0 0 2px; font-weight:600; }
  .muted { color:var(--text-muted); font-size:12px; }
  .section { padding:10px 16px; border-bottom:1px solid var(--border); }
  .section h2 { font-size:12px; text-transform:uppercase; letter-spacing:.04em; color:var(--text-secondary);
                margin:0 0 6px; font-weight:600; display:flex; justify-content:space-between; align-items:center; }
  .links button { background:none; border:0; color:var(--accent); cursor:pointer; font-size:12px; padding:0 0 0 8px; }
  #years { display:flex; flex-wrap:wrap; gap:6px; }
  .chip { display:inline-flex; align-items:center; gap:6px; padding:3px 9px 3px 7px; border:1px solid var(--border);
          border-radius:999px; cursor:pointer; user-select:none; font-size:13px; background:var(--surface-1); }
  .chip.off { opacity:.45; }
  .swatch { width:10px; height:10px; border-radius:2px; flex:none; }
  #search { width:100%; padding:6px 8px; border:1px solid var(--border); border-radius:6px;
            background:var(--surface-1); color:var(--text-primary); font:inherit; }
  #list { flex:1; overflow:auto; padding:4px 8px 12px; }
  .row { display:flex; align-items:flex-start; gap:8px; padding:5px 8px; border-radius:6px; cursor:pointer; }
  .row input { margin-top:3px; } .row .swatch { margin-top:4px; }
  .row .txt { flex:1; min-width:0; }
  .row:hover, .row.hl { background:var(--surface-2); }
  .row input { margin:0; }
  .row .name { display:block; overflow-wrap:anywhere; }
  .row .dates { display:block; color:var(--text-muted); font-size:12px; }
  #map { flex:1; }
  .leaflet-tooltip { font:12px system-ui,sans-serif; }
  .pop b { display:block; margin-bottom:2px; }
  .pop td { padding:0 8px 0 0; vertical-align:top; } .pop td:first-child { color:#666; }
  .legend { background:var(--surface-1); color:var(--text-primary); padding:6px 9px; border-radius:6px;
            border:1px solid var(--border); font-size:12px; line-height:1.6; max-height:40vh; overflow:auto; }
  .legend div { display:flex; align-items:center; gap:6px; }
  .legend .bar { width:18px; height:3px; border-radius:2px; }
</style>
</head>
<body>
<div id="app">
  <aside id="panel">
    <header><h1>__TITLE__</h1><div class="muted" id="count"></div></header>
    <div class="section">
      <h2>Start year <span class="links"><button id="yAll">All</button><button id="yNone">None</button></span></h2>
      <div id="years"></div>
    </div>
    <div class="section">
      <h2>Deployments <span class="links"><button id="dAll">All</button><button id="dNone">None</button><button id="zoom">Zoom to</button></span></h2>
      <input id="search" type="search" placeholder="Filter by name (comma = OR)">
    </div>
    <div id="list"></div>
  </aside>
  <div id="map"></div>
</div>
<script id="data" type="application/json">__DATA__</script>
<script>
(function () {
  const DATA = JSON.parse(document.getElementById('data').textContent);
  const deps = DATA.deployments;
  // Ordinal blue ramp (steps 250..700) - one hue, light -> dark = earlier -> later year.
  const RAMP = ['#86b6ef','#6da7ec','#5598e7','#3987e5','#2a78d6','#256abf','#1c5cab','#184f95','#104281','#0d366b'];
  const HL = getComputedStyle(document.documentElement).getPropertyValue('--highlight').trim() || '#eb6834';

  const fmtDate = s => new Date(s * 1000).toISOString().slice(0, 10);
  const fmtTime = s => new Date(s * 1000).toISOString().replace('T', ' ').slice(0, 19) + ' UTC';
  const yearOf = s => new Date(s * 1000).getUTCFullYear();

  // Each deployment is filed under its START year (matches the year in DAC/PyGlider names;
  // a Dec 2018 - Mar 2019 austral-summer deployment is a 2018 deployment).
  deps.forEach(d => { d.years = [yearOf(d.start)]; });
  const years = [...new Set(deps.flatMap(d => d.years))].sort((a, b) => a - b);
  const yearColor = {};
  years.forEach((y, i) => {
    const k = years.length === 1 ? 4 : Math.round(i * (RAMP.length - 1) / (years.length - 1));
    yearColor[y] = RAMP[k];
  });
  deps.forEach(d => { d.color = yearColor[d.years[0]]; d.checked = true; });

  const selYears = new Set(years);
  let query = [];

  // ---- map
  const map = L.map('map', { preferCanvas: true, worldCopyJump: true });
  const ocean = L.tileLayer('https://server.arcgisonline.com/ArcGIS/rest/services/Ocean/World_Ocean_Base/MapServer/tile/{z}/{y}/{x}',
    { maxZoom: 13, attribution: 'Tiles &copy; Esri, GEBCO, NOAA, NGS' });
  const light = L.tileLayer('https://{s}.basemaps.cartocdn.com/light_all/{z}/{x}/{y}{r}.png',
    { maxZoom: 19, subdomains: 'abcd', attribution: '&copy; OpenStreetMap contributors &copy; CARTO' });
  ocean.addTo(map);
  L.control.layers({ 'Ocean (Esri)': ocean, 'Light (CARTO)': light }, null, { position: 'topright' }).addTo(map);
  L.control.scale({ imperial: false }).addTo(map);
  const renderer = L.canvas({ tolerance: 6 });   // generous hover/click target

  deps.forEach(d => {
    const pts = d.lat.map((la, i) => [la, d.lon[i]]);
    d.line = L.polyline(pts, { color: d.color, weight: 2, opacity: 0.9, renderer });
    d.startMk = L.circleMarker(pts[0], { radius: 4, color: '#fff', weight: 2, fillColor: d.color, fillOpacity: 1, renderer });
    d.line.bindTooltip(d.name, { sticky: true });
    d.startMk.bindTooltip(d.name + ' (start)');
    d.line.on('mouseover', () => highlight(d, true));
    d.line.on('mouseout', () => highlight(d, false));
    d.line.on('click', e => openPopup(d, e.latlng));
    d.startMk.on('click', e => openPopup(d, e.latlng));
  });

  function highlight(d, on) {
    d.line.setStyle(on ? { color: HL, weight: 4, opacity: 1 } : { color: d.color, weight: 2, opacity: 0.9 });
    if (on) d.line.bringToFront();
    if (d.row) d.row.classList.toggle('hl', on);
  }

  function openPopup(d, ll) {
    // nearest decimated point to the click
    let best = 0, bd = Infinity;
    const cosl = Math.cos(ll.lat * Math.PI / 180);
    for (let i = 0; i < d.lat.length; i++) {
      const dy = d.lat[i] - ll.lat, dx = (d.lon[i] - ll.lng) * cosl, q = dx * dx + dy * dy;
      if (q < bd) { bd = q; best = i; }
    }
    const days = ((d.end - d.start) / 86400).toFixed(1);
    L.popup().setLatLng([d.lat[best], d.lon[best]]).setContent(
      `<div class="pop"><b>${esc(d.name)}</b><table>
       <tr><td>Start</td><td>${fmtTime(d.start)}</td></tr>
       <tr><td>End</td><td>${fmtTime(d.end)}</td></tr>
       <tr><td>Duration</td><td>${days} days</td></tr>
       <tr><td>Positions</td><td>${d.n.toLocaleString()} (${d.files.toLocaleString()} files)</td></tr>
       <tr><td colspan="2" style="padding-top:6px;color:inherit"><b>Nearest point</b></td></tr>
       <tr><td>Time</td><td>${fmtTime(d.t[best])}</td></tr>
       <tr><td>Lat, Lon</td><td>${d.lat[best].toFixed(4)}, ${d.lon[best].toFixed(4)}</td></tr>
       </table>${d.link ? `<div style="margin-top:6px"><a href="${esc(d.link)}" target="_top">Open 3D curtain →</a></div>` : ''}</div>`).openOn(map);
  }
  const esc = s => String(s).replace(/[&<>"]/g, c => ({ '&': '&amp;', '<': '&lt;', '>': '&gt;', '"': '&quot;' }[c]));

  // ---- legend (year colors)
  const legend = L.control({ position: 'bottomright' });
  legend.onAdd = () => {
    const div = L.DomUtil.create('div', 'legend');
    div.innerHTML = '<div style="font-weight:600">Start year</div>' +
      years.map(y => `<div><span class="bar" style="background:${yearColor[y]}"></span>${y}</div>`).join('') +
      `<div><span class="bar" style="background:${HL}"></span>hovered</div>`;
    return div;
  };
  legend.addTo(map);

  // ---- controls
  const yearsEl = document.getElementById('years');
  years.forEach(y => {
    const chip = document.createElement('span');
    chip.className = 'chip';
    chip.innerHTML = `<span class="swatch" style="background:${yearColor[y]}"></span>${y}`;
    chip.onclick = () => { selYears.has(y) ? selYears.delete(y) : selYears.add(y); render(); };
    chip.dataset.year = y;
    yearsEl.appendChild(chip);
  });
  document.getElementById('yAll').onclick = () => { years.forEach(y => selYears.add(y)); render(); };
  document.getElementById('yNone').onclick = () => { selYears.clear(); render(); };
  document.getElementById('search').oninput = e => {
    query = e.target.value.toLowerCase().split(',').map(s => s.trim()).filter(Boolean);
    render();
  };
  document.getElementById('dAll').onclick = () => { shown().forEach(d => d.checked = true); render(); };
  document.getElementById('dNone').onclick = () => { shown().forEach(d => d.checked = false); render(); };
  document.getElementById('zoom').onclick = () => fit(visible());

  const inYears = d => d.years.some(y => selYears.has(y));
  const matches = d => !query.length || query.some(q => d.name.toLowerCase().includes(q));
  const shown = () => deps.filter(d => inYears(d) && matches(d));          // in the list
  const visible = () => shown().filter(d => d.checked);                     // on the map

  const listEl = document.getElementById('list');
  function render() {
    [...yearsEl.children].forEach(c => c.classList.toggle('off', !selYears.has(+c.dataset.year)));
    const s = shown(), v = new Set(visible());
    listEl.innerHTML = '';
    s.sort((a, b) => a.start - b.start).forEach(d => {
      const row = document.createElement('label');
      row.className = 'row';
      row.innerHTML = `<input type="checkbox" ${d.checked ? 'checked' : ''}>
        <span class="swatch" style="background:${d.color}"></span>
        <span class="txt"><span class="name">${esc(d.name)}</span>
        <span class="dates">${fmtDate(d.start)} → ${fmtDate(d.end)}</span></span>`;
      row.querySelector('input').onchange = e => { d.checked = e.target.checked; render(); };
      row.onmouseenter = () => { if (v.has(d)) highlight(d, true); };
      row.onmouseleave = () => { if (v.has(d)) highlight(d, false); };
      row.ondblclick = e => { e.preventDefault(); fit([d]); };
      d.row = row;
      listEl.appendChild(row);
    });
    deps.forEach(d => {
      const on = v.has(d);
      if (on && !map.hasLayer(d.line)) { d.line.addTo(map); d.startMk.addTo(map); }
      if (!on && map.hasLayer(d.line)) { map.removeLayer(d.line); map.removeLayer(d.startMk); }
    });
    document.getElementById('count').textContent =
      `${v.size} of ${deps.length} deployments on map · ${s.length} listed · double-click a name to zoom`;
  }

  function fit(ds) {
    if (!ds.length) return;
    let b = L.latLngBounds([]);
    ds.forEach(d => b.extend(d.line.getBounds()));
    map.fitBounds(b.pad(0.05));
  }

  render();
  fit(deps);
})();
</script>
</body>
</html>
"""


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0],
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("paths", nargs="+", help="folders (searched recursively) and/or .nc files")
    ap.add_argument("-o", "--out-dir", default="glider_output")
    ap.add_argument("--positions", choices=["measurement", "profile"], default="measurement",
                    help="'measurement' = time/lat/lon of every sample (default); "
                         "'profile' = profile_time/profile_lat/profile_lon, one per profile (much smaller)")
    ap.add_argument("--group-by", choices=["auto", "trajectory", "folder", "file"], default="auto",
                    help="how to name deployments (default auto: trajectory variable, then "
                         "deployment_name attribute, then parent folder)")
    ap.add_argument("--max-map-points", type=int, default=2000,
                    help="max points per deployment drawn on the map (CSV keeps all). 0 = no limit")
    ap.add_argument("--workers", type=int, default=max(1, min(8, (os.cpu_count() or 2) - 1)),
                    help="parallel file readers (default: CPU count - 1, max 8)")
    ap.add_argument("--title", default="Glider deployments")
    ap.add_argument("--link-template", default=None,
                    help="URL for each deployment's page, with {name} as placeholder "
                         "(adds an 'Open 3D curtain' link to the map popups)")
    ap.add_argument("--map-name", default="glider_map.html", help="file name for the map page")
    a = ap.parse_args(argv)
    if a.group_by == "trajectory":
        a.group_by = "auto"
    build(a.paths, a.out_dir, a.positions, a.group_by, a.max_map_points, a.workers, a.title,
          a.link_template, a.map_name)


if __name__ == "__main__":
    main()
