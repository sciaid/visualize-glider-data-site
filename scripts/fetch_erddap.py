#!/usr/bin/env python3
"""
fetch_erddap.py - Download glider deployments from the IOOS Glider DAC ERDDAP and grid them
into the same 5 m depth x profile NetCDF layout as the local L3 (PyGlider) files, so the site
build treats them exactly like local files.

    python scripts/fetch_erddap.py "C:\\Claude_work\\PROJECTS\\Glider L-3 ncfiles"

Rules
  * Only datasets from INSTITUTIONS below (the older "Antarctic Ecosystem Research Division"
    datasets are deliberately not included).
  * One dataset per deployment: delayed-mode ("...-delayed") if it exists, otherwise real-time.
    When the delayed dataset appears, the real-time grid is deleted and replaced.
  * A deployment that already has a local gridded file anywhere under the folder (matched on its
    deployment_name, e.g. "calanus-20260824") is skipped: local files win.
  * Datasets are re-downloaded only when ERDDAP reports a change (date_modified or maxTime),
    tracked in <folder>/ERDDAP/manifest.json. Raw downloads are deleted after gridding.

Output: <folder>/ERDDAP/<deployment>_erddap-<delayed|realtime>-5m.nc
Gridding check: for calanus-20260824 the gridded ERDDAP data matched the local L3 file in
73,826 of 73,828 cells (the other 2 differ by one sample on a bin edge, < 0.001 C).

Requirements: pip install netCDF4 numpy
"""
from __future__ import annotations

import argparse
import csv
import io
import json
import re
import sys
import urllib.parse
import urllib.request
from pathlib import Path

import numpy as np
import netCDF4

ERDDAP = "https://gliders.ioos.us/erddap"
INSTITUTIONS = ["NOAA National Marine Fisheries Service Southwest Fisheries Science Center"]
# Variables downloaded (only those present in each dataset are requested)
SAMPLE_VARS = ["temperature", "salinity", "chlorophyll", "oxygen_concentration",
               "potential_temperature", "potential_density", "pressure"]
BASE_VARS = ["profile_id", "time", "latitude", "longitude", "depth"]
KEEP_ATTRS = ["project", "institution", "creator_name", "summary", "title", "glider_name",
              "platform_type", "wmo_id", "date_modified", "license", "acknowledgment"]
DZ, MAX_DEPTH = 5.0, 1200.0
BATHY_KEYS = ("gebco", "etopo", "bathy")
ID_RE = re.compile(r"^(?P<base>.+?-\d{8})T\d{4}(?P<delayed>-delayed)?$")


# --------------------------------------------------------------------------- #
# ERDDAP access
# --------------------------------------------------------------------------- #
def _get(url, timeout=600):
    req = urllib.request.Request(url, headers={"User-Agent": "visualize-glider-data-site"})
    with urllib.request.urlopen(req, timeout=timeout) as r:
        return r.read()


def list_datasets(institutions=INSTITUTIONS):
    """Rows of allDatasets for the given institutions: datasetID, minTime, maxTime, lat/lon."""
    cols = "datasetID,institution,minTime,maxTime,minLatitude,maxLatitude,minLongitude,maxLongitude"
    out = []
    for inst in institutions:
        q = cols + "&institution=" + urllib.parse.quote(json.dumps(inst))
        text = _get(f"{ERDDAP}/tabledap/allDatasets.csv?{q}").decode("utf-8")
        rows = list(csv.reader(io.StringIO(text)))
        head = rows[0]
        for r in rows[2:]:                    # row 1 holds units
            if r:
                out.append(dict(zip(head, r)))
    return out


def dataset_info(dataset_id):
    """Variable names and NC_GLOBAL attributes from /info/<id>/index.csv."""
    text = _get(f"{ERDDAP}/info/{dataset_id}/index.csv").decode("utf-8")
    variables, globals_, units = [], {}, {}
    for r in csv.DictReader(io.StringIO(text)):
        if r["Row Type"] == "variable":
            variables.append(r["Variable Name"])
        elif r["Row Type"] == "attribute" and r["Variable Name"] == "NC_GLOBAL":
            globals_[r["Attribute Name"]] = r["Value"]
        elif r["Row Type"] == "attribute" and r["Attribute Name"] == "units":
            units[r["Variable Name"]] = r["Value"]
    return {"variables": variables, "globals": globals_, "units": units}


def download(dataset_id, variables, dest):
    url = f"{ERDDAP}/tabledap/{dataset_id}.nc?" + ",".join(variables)
    dest.write_bytes(_get(url, timeout=1800))
    return url


# --------------------------------------------------------------------------- #
# Selection
# --------------------------------------------------------------------------- #
def parse_id(dataset_id):
    m = ID_RE.match(dataset_id)
    if not m:
        return None, None
    return m.group("base"), ("delayed" if m.group("delayed") else "realtime")


def choose(datasets):
    """{deployment: row}, delayed preferred over real-time."""
    best = {}
    for row in datasets:
        base, mode = parse_id(row["datasetID"])
        if base is None:
            continue
        row = dict(row, deployment=base, mode=mode)
        if base not in best or (mode == "delayed" and best[base]["mode"] != "delayed"):
            best[base] = row
    return best


def local_deployments(root, erddap_dir):
    """deployment_name of every gridded (depth x time) NetCDF under root, outside erddap_dir."""
    names = {}
    import os
    paths = sorted(Path(dp) / n for dp, _, ns in os.walk(root, followlinks=True) for n in ns if n.lower().endswith(".nc"))
    for p in paths:
        if erddap_dir in p.parents or any(k in p.name.lower() for k in BATHY_KEYS):
            continue
        try:
            with netCDF4.Dataset(p) as ds:
                if "depth" in ds.dimensions and "time" in ds.dimensions:
                    names[str(getattr(ds, "deployment_name", p.stem)).strip()] = p
        except OSError:
            continue
    return names


# --------------------------------------------------------------------------- #
# Gridding
# --------------------------------------------------------------------------- #
def _f(var):
    a = var[:]
    return np.ma.filled(a.astype("f8"), np.nan) if np.ma.isMaskedArray(a) else np.asarray(a, "f8")


def grid_dac(path, variables, dz=DZ, max_depth=MAX_DEPTH):
    """Bin each profile's samples into [z, z+dz) depth bins (bin means): depth x profile arrays."""
    ds = netCDF4.Dataset(path)
    pid, depth = _f(ds["profile_id"]), _f(ds["depth"])
    good = np.isfinite(pid) & np.isfinite(depth) & (depth >= 0) & (depth < max_depth)
    pid, depth = pid[good], depth[good]
    profiles, inv = np.unique(pid, return_inverse=True)
    edges = np.arange(0, max_depth + dz, dz)
    centers = edges[:-1] + dz / 2
    k = np.clip(np.searchsorted(edges, depth, side="right") - 1, 0, centers.size - 1)
    ncell, cell = centers.size * profiles.size, k * profiles.size + inv
    out = {"depth": centers, "profile": profiles}
    for name in ("time", "latitude", "longitude"):           # one value per profile
        v = _f(ds[name])[good]
        f = np.isfinite(v)
        s = np.bincount(inv[f], weights=v[f], minlength=profiles.size)
        n = np.bincount(inv[f], minlength=profiles.size)
        out[name] = np.where(n > 0, s / np.maximum(n, 1), np.nan)
    for name in variables:
        if name not in ds.variables:
            continue
        v = _f(ds[name])[good]
        f = np.isfinite(v)
        s = np.bincount(cell[f], weights=v[f], minlength=ncell)
        n = np.bincount(cell[f], minlength=ncell)
        out[name] = np.where(n > 0, s / np.maximum(n, 1), np.nan).reshape(centers.size, profiles.size)
    ds.close()
    # drop depth bins below the deepest data so files stay small
    rows = [k for k in range(centers.size) if any(np.isfinite(out[v][k]).any() for v in variables if v in out)]
    nd = (max(rows) + 1) if rows else 1
    out["depth"] = centers[:nd]
    for v in variables:
        if v in out:
            out[v] = out[v][:nd]
    return out


def fill_density(g):
    """Potential density (kg m-3) from potential temperature (or temperature) and salinity when the
    dataset lacks it. EOS-80 one-atmosphere equation (Millero & Poisson 1981), as on the site."""
    if "potential_density" in g or "salinity" not in g:
        return
    theta = g.get("potential_temperature", g.get("temperature"))
    if theta is None:
        return
    sys.path.insert(0, str(Path(__file__).resolve().parent))
    from glider_curtain import sigma0_eos80
    g["potential_density"] = sigma0_eos80(g["salinity"], theta) + 1000.0


def write_grid(g, info, row, url, path):
    ds = netCDF4.Dataset(path, "w")
    ds.createDimension("depth", g["depth"].size)
    ds.createDimension("time", g["profile"].size)
    gl = info["globals"]
    for a in KEEP_ATTRS:
        if a in gl:
            ds.setncattr(a, gl[a])
    ds.deployment_name = row["deployment"]
    ds.glider_name = gl.get("glider_name") or row["deployment"].rsplit("-", 1)[0]
    ds.data_mode = row["mode"]
    ds.source_dataset = row["datasetID"]
    ds.source_url = f"{ERDDAP}/tabledap/{row['datasetID']}.html"
    ds.history = f"Gridded to {DZ:g} m bins (bin means per profile) by fetch_erddap.py from {url}"
    t = ds.createVariable("time", "f8", ("time",)); t[:] = g["time"]
    t.units = "seconds since 1970-01-01T00:00:00Z"; t.standard_name = "time"
    for name, std in (("latitude", "latitude"), ("longitude", "longitude")):
        v = ds.createVariable(name, "f8", ("time",)); v[:] = g[name]; v.standard_name = std
        v.units = "degrees_north" if name == "latitude" else "degrees_east"
    v = ds.createVariable("profile", "i4", ("time",)); v[:] = g["profile"].astype("i4")
    v = ds.createVariable("depth", "f4", ("depth",)); v[:] = g["depth"]; v.units = "m"; v.positive = "down"
    for name in SAMPLE_VARS:
        if name in g:
            v = ds.createVariable(name, "f4", ("depth", "time"), fill_value=np.float32(np.nan), zlib=True)
            v[:] = g[name].astype("f4")
            v.units = info["units"].get(name, "kg m-3" if name == "potential_density" else "")
    ds.close()


# --------------------------------------------------------------------------- #
def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("nc_dir", help="root folder of the glider files; grids go to <nc_dir>/ERDDAP")
    ap.add_argument("--dry-run", action="store_true", help="list what would be downloaded, change nothing")
    ap.add_argument("--force", action="store_true", help="re-download every selected dataset")
    a = ap.parse_args(argv)

    root = Path(a.nc_dir)
    out_dir = root / "ERDDAP"
    manifest_path = out_dir / "manifest.json"
    manifest = json.loads(manifest_path.read_text()) if manifest_path.exists() else {}

    chosen = choose(list_datasets())
    local = local_deployments(root, out_dir)
    print(f"{len(chosen)} deployments on ERDDAP; {len(local)} local gridded files")

    plan = []
    for dep, row in sorted(chosen.items()):
        if dep in local:
            print(f"  {dep:28s} local file wins ({Path(local[dep]).name})")
            continue
        prev = manifest.get(dep, {})
        changed = (a.force or prev.get("datasetID") != row["datasetID"]
                   or prev.get("maxTime") != row["maxTime"] or not (out_dir / prev.get("file", "-")).exists())
        print(f"  {dep:28s} {row['mode']:8s} {row['datasetID']:40s} {'download' if changed else 'up to date'}")
        if changed:
            plan.append(row)

    # grids no longer wanted: deployment now local, gone from ERDDAP, or superseded (real-time -> delayed)
    keep = {f"{d}_erddap-{r['mode']}-5m.nc" for d, r in chosen.items() if d not in local}
    stale = [p for p in out_dir.glob("*_erddap-*-5m.nc") if p.name not in keep] if out_dir.exists() else []
    for p in stale:
        print(f"  remove {p.name}")
    if a.dry_run:
        print(f"\nDry run: {len(plan)} to download, {len(stale)} to remove.")
        return

    out_dir.mkdir(exist_ok=True)
    for p in stale:
        p.unlink()
        manifest = {d: m for d, m in manifest.items() if m.get("file") != p.name}
    for row in plan:
        dep, did = row["deployment"], row["datasetID"]
        info = dataset_info(did)
        want = BASE_VARS + [v for v in SAMPLE_VARS if v in info["variables"]]
        raw = out_dir / f"_download_{did}.nc"
        print(f"\nDownloading {did} ({len(want)} variables) ...")
        url = download(did, want, raw)
        g = grid_dac(raw, [v for v in SAMPLE_VARS if v in info["variables"]])
        fill_density(g)
        name = f"{dep}_erddap-{row['mode']}-5m.nc"
        write_grid(g, info, row, url, out_dir / name)
        raw.unlink()
        manifest[dep] = {"datasetID": did, "mode": row["mode"], "maxTime": row["maxTime"],
                         "date_modified": info["globals"].get("date_modified", ""), "file": name}
        print(f"  {g['profile'].size} profiles x {g['depth'].size} bins -> {name}")
        manifest_path.write_text(json.dumps(manifest, indent=1), encoding="utf-8")
    manifest_path.write_text(json.dumps(manifest, indent=1), encoding="utf-8")
    print(f"\nDone: {len(plan)} downloaded, {len(stale)} removed. Next: python scripts/build_site.py \"{root}\"")


if __name__ == "__main__":
    main()
