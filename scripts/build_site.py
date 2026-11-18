#!/usr/bin/env python3
"""
build_site.py - Rebuild the data behind the website from the glider NetCDF files.

Run this on the computer that has the .nc files, then commit and push:

    python scripts/build_site.py "C:\\Claude_work\\PROJECTS\\Glider L-3 ncfiles"

Several folders can be given (for example one per region). Folders are searched
recursively for *.nc files.

Only new or changed NetCDF files are reprocessed: each file's fingerprint (size plus a
hash of its first and last MB) is stored in viewer/data/index.json. That keeps unchanged
data files byte-for-byte identical, so git does not store new copies of them. Use
--force to rebuild everything, e.g. after changing how the data are processed.

It writes (relative to the repository root):
    map/index.html                 overview map of all deployments (always rebuilt; small)
    viewer/index.html              shared 3D curtain viewer (+ plotly.min.js)
    viewer/data/<name>.json        one data file per deployment
    viewer/data/index.json         list of deployments, with region and source fingerprint
    deployments/<name>.qmd         one page per deployment
    deployments.qmd                table of all deployments, grouped by region
    _summary.md                    per-region summary table shown on the home page

Regions come from regions.json (name + lat/lon box, matched on each deployment's median
position). Bathymetry: every GEBCO/ETOPO grid (file name containing gebco, etopo or bathy)
in the folders, or one folder down, is found; each deployment uses the grid that covers it.

Requirements: pip install netCDF4 numpy contourpy plotly
"""
import argparse
import hashlib
import json
import sys
from collections import OrderedDict
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
sys.path.insert(0, str(HERE))

import glider_curtain  # noqa: E402
import glider_tracks   # noqa: E402

BATHY_KEYS = ("gebco", "etopo", "bathy")

DEPLOYMENT_PAGE = """---
title: "{name}"
page-layout: full
toc: false
---

| Region | Project | Glider | Start | End | Profiles | Max depth | Distance (approx.) | Avg. per day |
|---|---|---|---|---|---|---|---|---|
| {region} | {project} | {glider} | {start} | {end} | {profiles:,} | {maxDepth:.0f} m | {distance} | {per_day} |

{realtime_note}**Area:** {area}{project_note}

**Data:** {source}

<iframe src="../viewer/index.html?d={name}" title="3D curtain viewer: {name}"
        style="width:100%; height:86vh; border:1px solid #e2e0db; border-radius:6px;"></iframe>

[Open the viewer full screen](../viewer/index.html?d={name}){{target="_blank"}} ·
[How these plots are made](../methods.qmd)

<small>Distance is summed between consecutive profile positions, not GPS surface fixes, so treat it as
approximate (it depends on how finely the track is sampled); see [Methods](../methods.qmd#track-distance).</small>
"""


def fmt_distance(e):
    """'≈ 2,000 km', rounded to 10 km; blank if not computed."""
    km = e.get("distanceKm")
    return "" if km is None else f"≈ {round(km, -1):,.0f} km"


def fmt_per_day(e):
    """Average over the whole deployment: distance / (last - first profile time)."""
    km, days = e.get("distanceKm"), e.get("days")
    return f"≈ {km / days:.0f} km/day" if km is not None and days else ""




def file_meta(path):
    """Deployment name and project from a file's global attributes (cheap: no data read).
    Project is shortened to the text before any ':' ("FREEBYRD: Estimate ..." -> "FREEBYRD")."""
    import netCDF4
    with netCDF4.Dataset(path) as ds:
        name = str(getattr(ds, "deployment_name", "") or Path(path).stem).strip()
        full = str(getattr(ds, "project", "") or "").strip()
        mode = str(getattr(ds, "data_mode", "") or "").strip()
        src_id = str(getattr(ds, "source_dataset", "") or "").strip()
        src_url = str(getattr(ds, "source_url", "") or "").strip()
    return {"name": name, "project": full.split(":")[0].strip(), "projectFull": full,
            "mode": mode, "sourceDataset": src_id, "sourceUrl": src_url}


def is_gridded(path):
    """True for depth x profile grids (local L3 or fetch_erddap output); False for anything else,
    e.g. a raw per-sample ERDDAP download left in the folder."""
    import netCDF4
    try:
        with netCDF4.Dataset(path) as ds:
            return ("depth" in ds.dimensions and "time" in ds.dimensions
                    and any(v.dimensions == ("depth", "time") for v in ds.variables.values()))
    except OSError:
        return False


REALTIME_NOTE = """::: {{.callout-warning}}
## Real-time data (provisional)
Data through {end}. Real-time data are sent while the glider is deployed: they are thinned,
not quality-controlled and may have gaps. This page is replaced by the delayed-mode data once
it is published.
:::

"""


def fmt_source(e):
    if e.get("sourceUrl"):
        return f"[{e['sourceDataset']}]({e['sourceUrl']}) on the IOOS Glider DAC ERDDAP" + \
               (" (real-time)" if e.get("mode") == "realtime" else "")
    return "SWFSC ESD processed (L3) file"


def fmt_area(e):
    """'63.76–62.21°S, 61.85–57.32°W' (ranges ordered by absolute value)."""
    def rng(a, b, pos, neg):
        h = pos if (a + b) / 2 >= 0 else neg
        lo, hi = sorted((abs(a), abs(b)))
        return f"{lo:.2f}–{hi:.2f}°{h}"
    return f"{rng(*e['latRange'], 'N', 'S')}, {rng(*e['lonRange'], 'E', 'W')}"


def write_summary(index):
    """_summary.md, included on the home page: deployments per region, years, projects, distance."""
    from collections import OrderedDict as OD
    groups = OD()
    for e in index:
        groups.setdefault(e.get("region") or "Other", []).append(e)
    rows = []
    for region, es in groups.items():
        years = sorted({int(e["start"][:4]) for e in es})
        yr = f"{years[0]}–{years[-1]}" if len(years) > 1 else f"{years[0]}"
        projects = ", ".join(sorted({e.get("project") for e in es if e.get("project")})) or "—"
        km = sum(e.get("distanceKm") or 0 for e in es)
        nrt = sum(e.get("mode") == "realtime" for e in es)
        count = f"{len(es)}" + (f" ({nrt} real-time)" if nrt else "")
        rows.append(f"| {region} | {count} | {yr} | {projects} | ≈ {round(km, -2):,.0f} km |")
    total = sum(e.get("distanceKm") or 0 for e in index)
    if len(groups) > 1:
        rows.append(f"| **All** | **{len(index)}** | | | **≈ {round(total, -2):,.0f} km** |")
    (ROOT / "_summary.md").write_text(
        "| Region | Deployments | Years (start) | Projects | Distance (approx.) |\n"
        "|---|---|---|---|---|\n" + "\n".join(rows) + "\n", encoding="utf-8")


def walk_nc(folder):
    """Every *.nc under folder, following linked folders too (sorted)."""
    import os
    out = []
    for dirpath, _, names in os.walk(folder, followlinks=True):
        out += [Path(dirpath) / n for n in names if n.lower().endswith(".nc")]
    return sorted(out)


def fingerprint(path, chunk=1 << 20):
    """Size + SHA-1 of the first and last MB: cheap, and the same on any computer."""
    p = Path(path)
    size = p.stat().st_size
    h = hashlib.sha1(str(size).encode())
    with open(p, "rb") as fh:
        h.update(fh.read(chunk))
        if size > chunk:
            fh.seek(max(chunk, size - chunk))
            h.update(fh.read(chunk))
    return f"{size}-{h.hexdigest()[:16]}"


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("nc_dirs", nargs="+", help="folder(s) with the glider *.nc files (searched recursively)")
    ap.add_argument("--force", action="store_true", help="rebuild every deployment, not just new or changed ones")
    ap.add_argument("--prune", action="store_true",
                    help="drop deployments whose NetCDF file is not in the given folders")
    a = ap.parse_args()

    files, skipped = [], []
    for d in a.nc_dirs:
        for p in walk_nc(d):
            if any(k in p.name.lower() for k in BATHY_KEYS):
                continue
            (files if is_gridded(p) else skipped).append(str(p))
    for p in skipped:
        print(f"  skipping {Path(p).name}: not a gridded depth x profile file")
    # One file per deployment: a local L3 file wins over an ERDDAP grid of the same deployment.
    by_name = {}
    for f in files:
        m = file_meta(f)
        if m["name"] in by_name and not m["sourceDataset"]:
            print(f"  {m['name']}: using local file, ignoring {Path(by_name[m['name']]).name}")
            by_name[m["name"]] = f
        elif m["name"] in by_name:
            print(f"  {m['name']}: using {Path(by_name[m['name']]).name}, ignoring {Path(f).name}")
        else:
            by_name[m["name"]] = f
    files = sorted(by_name.values())
    if not files:
        sys.exit(f"No glider NetCDF files found in {a.nc_dirs}")
    regions = glider_tracks.load_regions(ROOT / "regions.json")
    data_dir = ROOT / "viewer" / "data"
    idx_path = data_dir / "index.json"
    index = json.loads(idx_path.read_text()) if idx_path.exists() else []

    # 1. Which files need (re)processing?
    fps = {f: fingerprint(f) for f in files}
    known = {e.get("file"): e for e in index if e.get("file")}
    todo = [f for f in files
            if a.force or Path(f).name not in known
            or known[Path(f).name].get("fingerprint") != fps[f]
            or not (data_dir / f"{known[Path(f).name]['name']}.json").exists()]
    print(f"{len(files)} glider files; {len(todo)} new or changed"
          + (" (--force)" if a.force else "") + f"; {len(files) - len(todo)} unchanged")

    meta = {f: file_meta(f) for f in files}

    # 2. Overview map (always rebuilt from all files: a single small page)
    print("\n== Overview map ==")
    glider_tracks.build(files, ROOT / "map", "measurement", "auto", 2000, 1, "Glider deployments",
                        link_template="../deployments/{name}.html", map_name="index.html", regions=regions,
                        extra={m["name"]: m for m in meta.values()})
    for extra in ("glider_positions.csv", "deployments_summary.csv"):
        (ROOT / "map" / extra).unlink(missing_ok=True)

    # 3. Curtain data for new/changed deployments (bathymetry grids found next to all files)
    print("\n== Curtain viewer data ==")
    grids = sorted({str(p) for d in a.nc_dirs for p in walk_nc(d)
                    if any(k in p.name.lower() for k in BATHY_KEYS)})
    bathy_arg = ",".join(grids) or "none"
    glider_curtain.main(todo + ["--site", str(ROOT / "viewer"), "--bathy", bathy_arg] if todo else
                        [files[0], "--site", str(ROOT / "viewer"), "--bathy", "none", "--viewer-only"])

    # 4. Record fingerprints and regions; optionally prune
    index = json.loads(idx_path.read_text())
    by_file = {Path(f).name: f for f in files}
    for e in index:
        if e.get("file") in by_file:
            f = by_file[e["file"]]
            e["fingerprint"] = fps[f]
            for k in ("project", "projectFull", "mode", "sourceDataset", "sourceUrl"):
                e[k] = meta[f][k]
        la = (e["latRange"][0] + e["latRange"][1]) / 2
        lo = (e["lonRange"][0] + e["lonRange"][1]) / 2
        e["region"] = glider_tracks.assign_region([la], [lo], regions) if regions else ""
    if a.prune:
        gone = [e for e in index if e.get("file") not in by_file]
        for e in gone:
            (data_dir / f"{e['name']}.json").unlink(missing_ok=True)
            (ROOT / "deployments" / f"{e['name']}.qmd").unlink(missing_ok=True)
            print(f"  pruned {e['name']}")
        index = [e for e in index if e.get("file") in by_file]
    index.sort(key=lambda e: (e.get("region", ""), e["start"]))
    idx_path.write_text(json.dumps(index, indent=1), encoding="utf-8")

    # 5. Pages
    print("\n== Pages ==")
    (ROOT / "deployments").mkdir(exist_ok=True)
    for e in index:
        page = DEPLOYMENT_PAGE.format(
            name=e["name"], region=e.get("region") or "", project=e.get("project") or "", distance=fmt_distance(e), per_day=fmt_per_day(e), glider=e.get("glider") or "",
            start=e["start"], end=e["end"], profiles=e["profiles"], maxDepth=e["maxDepth"],
            area=fmt_area(e), source=fmt_source(e),
            realtime_note=REALTIME_NOTE.format(end=e["end"]) if e.get("mode") == "realtime" else "",
            project_note=(f" · **Project:** {e['projectFull']}"
                          if e.get("projectFull") and e["projectFull"] != e.get("project") else ""))
        (ROOT / "deployments" / f"{e['name']}.qmd").write_text(page, encoding="utf-8")
    groups = OrderedDict()
    for e in index:
        groups.setdefault(e.get("region") or "Other", []).append(e)
    parts = []
    for region, es in groups.items():
        rows = "\n".join(
            f"| [{e['name']}](deployments/{e['name']}.qmd){' · *real-time*' if e.get('mode') == 'realtime' else ''} | {e.get('project', '')} | {e.get('glider', '')} | {e['start']} | {e['end']} "
            f"| {e['profiles']:,} | {e['maxDepth']:.0f} m | {fmt_distance(e)} | {fmt_per_day(e)} |" for e in es)
        heading = f"## {region} ({len(es)})\n\n" if len(groups) > 1 else ""
        parts.append(f"{heading}| Deployment | Project | Glider | Start | End | Profiles | Max depth | Distance (approx.) | Avg. per day |\n"
                     f"|---|---|---|---|---|---|---|---|---|\n{rows}\n")
    (ROOT / "deployments.qmd").write_text(
        "---\ntitle: \"Deployments\"\n---\n\n"
        f"{len(index)} deployments. Click a name for its 3D curtain, depth-time section and T–S diagram.\n\n"
        + "\n".join(parts), encoding="utf-8")
    write_summary(index)
    print(f"Wrote {len(index)} deployment pages. Next: commit and push (or open a pull request).")


if __name__ == "__main__":
    main()
