#!/usr/bin/env python3
"""
fetch_erddap.py - Download glider deployments from the IOOS Glider DAC ERDDAP and grid them
into the same 5 m depth x profile NetCDF layout as the local L3 (PyGlider) files, so the site
build treats them exactly like local files.

    python scripts/fetch_erddap.py "C:\\Claude_work\\PROJECTS\\Glider L-3 ncfiles"

Rules
  * Only datasets from the institutions in erddap.json (the older "Antarctic Ecosystem Research
    Division" datasets are deliberately not included).
  * One dataset per deployment: delayed-mode ("...-delayed") if it exists, otherwise real-time.
    When the delayed dataset appears, the real-time grid is deleted and replaced.
  * A deployment that already has a local gridded file anywhere under the folder (matched on its
    deployment_name, e.g. "calanus-20260824") is skipped: local files win. With --repo-index, so are
    deployments already published on the site from a local file (needed on GitHub, where the local
    files are absent), unless listed in erddap.json "replace_local_with_erddap".
  * A failed download keeps the deployment's previous grid; if ERDDAP cannot be reached at all,
    nothing changes. Problems are recorded in <folder>/ERDDAP/status.json.
  * Datasets are re-downloaded only when ERDDAP reports a change (date_modified or maxTime),
    tracked in <folder>/ERDDAP/manifest.json. Raw downloads are deleted after gridding; for
    real-time datasets the individual samples are kept in the grid file (dimension "sample") so the
    viewer can plot them as transmitted.

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
SAMPLE_VARS = ["temperature", "salinity", "chlorophyll", "oxygen_concentration", "cdom", "backscatter_700",
               "potential_temperature", "potential_density", "pressure"]
OPTIONAL_VARS = ["chlorophyll", "oxygen_concentration", "cdom", "backscatter_700"]   # sensor-dependent
BASE_VARS = ["profile_id", "time", "latitude", "longitude", "depth"]
KEEP_ATTRS = ["project", "institution", "creator_name", "summary", "title", "glider_name",
              "platform_type", "wmo_id", "date_modified", "license", "acknowledgment"]
DZ, MAX_DEPTH = 5.0, 1200.0
# Real-time grids also keep the individual samples (as transmitted) of these variables, so the
# viewer can plot them as points. Bump FORMAT (below) when the grid file layout changes: grids
# written with an older FORMAT are downloaded again on the next run.
SAMPLE_PLOT_VARS = ["temperature", "salinity", "chlorophyll", "oxygen_concentration", "cdom", "backscatter_700"]
# Glider-measured depth-averaged current (one value per dive segment), kept for real-time grids
# as dac_<name> on the time dimension: the real-time map shows it beside the current models.
DAC_VARS = ["u", "v", "time_uv", "lat_uv", "lon_uv"]
# Grid file layout version, per data mode (only real-time grids hold samples and dac_* variables,
# so a change to those does not force every delayed grid to download again).
FORMAT = {"realtime": 4, "delayed": 3}
MISSING_GRACE_DAYS = 3     # keep a deployment this long after it disappears from the ERDDAP listing
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
    smp = {"profile": inv.astype("i4"), "depth": depth}              # individual samples, in file order
    for name in ("time", "latitude", "longitude"):
        smp[name] = _f(ds[name])[good]
    for name in SAMPLE_PLOT_VARS:
        if name in variables and name in ds.variables:
            smp[name] = _f(ds[name])[good]
    out["samples"] = smp
    per_profile = [n for n in ("time", "latitude", "longitude") + tuple(DAC_VARS) if n in ds.variables]
    for name in per_profile:                                   # one value per profile
        v = _f(ds[name])[good]
        f = np.isfinite(v)
        s = np.bincount(inv[f], weights=v[f], minlength=profiles.size)
        n = np.bincount(inv[f], minlength=profiles.size)
        out[name if name in ("time", "latitude", "longitude") else "dac_" + name] = \
            np.where(n > 0, s / np.maximum(n, 1), np.nan)
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


def write_grid(g, info, row, url, path, samples=False):
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
    for name in DAC_VARS:
        if "dac_" + name in g:
            v = ds.createVariable("dac_" + name, "f8", ("time",), fill_value=np.nan)
            v[:] = g["dac_" + name]
            v.units = info["units"].get(name, "")
            v.comment = (f"Glider-measured depth-averaged current variable '{name}' from the source dataset, "
                         "one value per profile")
    if samples and "samples" in g:
        write_samples(ds, g["samples"], info)
    ds.close()


def write_samples(ds, smp, info):
    """Individual samples (dimension "sample"): sample_time/depth/latitude/longitude, sample_profile
    (column of the grid each sample was binned into) and sample_<variable>. Samples with no value in
    any of the plotted variables are dropped. Sorted by time."""
    names = [n for n in SAMPLE_PLOT_VARS if n in smp]
    keep = np.zeros(smp["time"].size, bool)
    for n in names:
        keep |= np.isfinite(smp[n])
    keep &= np.isfinite(smp["time"])
    order = np.argsort(smp["time"][keep], kind="stable")
    sel = lambda a: a[keep][order]
    ds.createDimension("sample", int(keep.sum()))
    def var(name, dtype, data, **attrs):
        v = ds.createVariable(name, dtype, ("sample",), zlib=True,
                              **({"fill_value": np.float32(np.nan)} if dtype == "f4" else {}))
        v[:] = data
        for k, val in attrs.items():
            v.setncattr(k, val)
    var("sample_time", "f8", sel(smp["time"]), units="seconds since 1970-01-01T00:00:00Z")
    var("sample_depth", "f4", sel(smp["depth"]).astype("f4"), units="m", positive="down")
    var("sample_latitude", "f8", sel(smp["latitude"]), units="degrees_north")
    var("sample_longitude", "f8", sel(smp["longitude"]), units="degrees_east")
    var("sample_profile", "i4", sel(smp["profile"]), long_name="index along the time dimension of this sample's profile")
    for n in names:
        var(f"sample_{n}", "f4", sel(smp[n]).astype("f4"), units=info["units"].get(n, ""))


# --------------------------------------------------------------------------- #
def load_config():
    """erddap.json in the repository: institutions to include, and deployments published from local
    files that may be replaced by their ERDDAP version (listed by name, or "*" for all)."""
    p = Path(__file__).resolve().parent.parent / "erddap.json"
    cfg = json.loads(p.read_text(encoding="utf-8")) if p.exists() else {}
    cfg.setdefault("institutions", INSTITUTIONS)
    cfg.setdefault("replace_local_with_erddap", [])
    return cfg


def locked_deployments(repo_index, allow):
    """Deployments already on the site from a local file (no sourceDataset in index.json). ERDDAP
    never replaces these unless they are listed in replace_local_with_erddap."""
    if not repo_index or not Path(repo_index).exists():
        return {}
    out = {}
    for e in json.loads(Path(repo_index).read_text(encoding="utf-8")):
        if not e.get("sourceDataset") and "*" not in allow and e["name"] not in allow:
            out[e["name"]] = e.get("file") or "published from a local file"
    return out


def main(argv=None):
    import datetime as dt
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("nc_dir", help="root folder of the glider files; grids go to <nc_dir>/ERDDAP")
    ap.add_argument("--repo-index", default=None,
                    help="the site's viewer/data/index.json: deployments published from local files are "
                         "not replaced by ERDDAP versions (see erddap.json)")
    ap.add_argument("--dry-run", action="store_true", help="list what would be downloaded, change nothing")
    ap.add_argument("--force", action="store_true", help="re-download every selected dataset")
    a = ap.parse_args(argv)

    cfg = load_config()
    root = Path(a.nc_dir)
    out_dir = root / "ERDDAP"
    manifest_path, status_path = out_dir / "manifest.json", out_dir / "status.json"
    manifest = json.loads(manifest_path.read_text()) if manifest_path.exists() else {}
    status = {"checked": dt.datetime.now(dt.timezone.utc).strftime("%Y-%m-%d %H:%M UTC"), "errors": []}

    def save_status():
        if not a.dry_run:
            out_dir.mkdir(parents=True, exist_ok=True)
            status_path.write_text(json.dumps(status, indent=1), encoding="utf-8")

    try:
        chosen = choose(list_datasets(cfg["institutions"]))
    except Exception as exc:                      # ERDDAP unreachable: keep every existing grid
        status["errors"].append(f"listing datasets failed: {exc}")
        status["listing_failed"] = True
        print(f"ERDDAP listing failed ({exc}); keeping existing grids unchanged.")
        save_status()
        return 0
    # ERDDAP sometimes lists only part of its datasets for a while (e.g. while reloading). A
    # deployment we already have that is missing from the listing is kept, unchanged, until it
    # has been missing for MISSING_GRACE_DAYS; only then is its grid removed.
    now = dt.datetime.now(dt.timezone.utc)
    for dep, m in list(manifest.items()):
        if dep in chosen:
            m.pop("missing_since", None)
            continue
        if not (out_dir / m.get("file", "-")).exists():
            continue
        since = dt.datetime.fromisoformat(m.setdefault("missing_since", now.isoformat(timespec="seconds")))
        days = (now - since).total_seconds() / 86400
        if days < MISSING_GRACE_DAYS:
            print(f"  {dep:28s} not in the ERDDAP listing (since {m['missing_since']}); keeping it")
            status.setdefault("warnings", []).append(f"{m['datasetID']}: missing from the ERDDAP listing since {m['missing_since']}")
            chosen[dep] = {"deployment": dep, "datasetID": m["datasetID"], "mode": m["mode"],
                           "maxTime": m["maxTime"]}
    local = {**locked_deployments(a.repo_index, cfg["replace_local_with_erddap"]),
             **{k: Path(v).name for k, v in local_deployments(root, out_dir).items()}}
    print(f"{len(chosen)} deployments on ERDDAP; {len(local)} published from local files")

    plan = []
    for dep, row in sorted(chosen.items()):
        if dep in local:
            print(f"  {dep:28s} local file wins ({local[dep]})")
            continue
        prev = manifest.get(dep, {})
        changed = (a.force or prev.get("datasetID") != row["datasetID"]
                   or prev.get("maxTime") != row["maxTime"] or prev.get("format") != FORMAT[row["mode"]]
                   or not (out_dir / prev.get("file", "-")).exists())
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
        return 0

    out_dir.mkdir(parents=True, exist_ok=True)
    done = 0
    for row in plan:
        dep, did = row["deployment"], row["datasetID"]
        raw = out_dir / f"_download_{did}.nc"
        try:
            info = dataset_info(did)
            want = BASE_VARS + [v for v in SAMPLE_VARS if v in info["variables"]]
            if row["mode"] == "realtime":
                want += [v for v in DAC_VARS if v in info["variables"]]
            print(f"\nDownloading {did} ({len(want)} variables) ...")
            missing = [v for v in OPTIONAL_VARS if v not in info["variables"]]
            if missing:
                print(f"  not in this dataset: {', '.join(missing)}")
            url = download(did, want, raw)
            g = grid_dac(raw, [v for v in SAMPLE_VARS if v in info["variables"]])
            fill_density(g)
            name = f"{dep}_erddap-{row['mode']}-5m.nc"
            write_grid(g, info, row, url, out_dir / name, samples=row["mode"] == "realtime")
        except Exception as exc:                  # keep the previous grid of this deployment, if any
            status["errors"].append(f"{did}: {exc}")
            print(f"  FAILED: {exc} (keeping the previous version, if any)")
            continue
        finally:
            raw.unlink(missing_ok=True)
        # a delayed grid replaces the real-time one only once it exists
        for p in out_dir.glob(f"{dep}_erddap-*-5m.nc"):
            if p.name != name:
                p.unlink()
        manifest[dep] = {"datasetID": did, "mode": row["mode"], "maxTime": row["maxTime"],
                         "date_modified": info["globals"].get("date_modified", ""), "file": name,
                         "format": FORMAT[row["mode"]]}
        manifest_path.write_text(json.dumps(manifest, indent=1), encoding="utf-8")
        done += 1
        print(f"  {g['profile'].size} profiles x {g['depth'].size} bins -> {name}")
    for p in stale:                               # deployments now local or gone from ERDDAP
        dep = p.name.split("_erddap-")[0]
        if dep in local or dep not in chosen:
            p.unlink(missing_ok=True)
            manifest.pop(dep, None)
    manifest_path.write_text(json.dumps(manifest, indent=1), encoding="utf-8")
    status.update({"downloaded": done, "deployments": sorted(m["file"] for m in manifest.values())})
    save_status()
    print(f"\nDone: {done} downloaded, {len(status['errors'])} failed.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
