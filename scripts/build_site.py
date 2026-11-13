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

| Region | Glider | Start | End | Profiles | Max depth | Area |
|---|---|---|---|---|---|---|
| {region} | {glider} | {start} | {end} | {profiles:,} | {maxDepth:.0f} m | {lat0:.2f}° to {lat1:.2f}°, {lon0:.2f}° to {lon1:.2f}° |

<iframe src="../viewer/index.html?d={name}" title="3D curtain viewer: {name}"
        style="width:100%; height:86vh; border:1px solid #e2e0db; border-radius:6px;"></iframe>

[Open the viewer full screen](../viewer/index.html?d={name}){{target="_blank"}} ·
[How these plots are made](../methods.qmd)
"""


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

    files = []
    for d in a.nc_dirs:
        files += [str(p) for p in sorted(Path(d).rglob("*.nc"))
                  if not any(k in p.name.lower() for k in BATHY_KEYS)]
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

    # 2. Overview map (always rebuilt from all files: a single small page)
    print("\n== Overview map ==")
    glider_tracks.build(files, ROOT / "map", "measurement", "auto", 2000, 1, "Glider deployments",
                        link_template="../deployments/{name}.html", map_name="index.html", regions=regions)
    for extra in ("glider_positions.csv", "deployments_summary.csv"):
        (ROOT / "map" / extra).unlink(missing_ok=True)

    # 3. Curtain data for new/changed deployments (bathymetry grids found next to all files)
    print("\n== Curtain viewer data ==")
    bathy_arg = ",".join(str(p) for p in glider_curtain.find_bathy(files)) or "none"
    glider_curtain.main(todo + ["--site", str(ROOT / "viewer"), "--bathy", bathy_arg] if todo else
                        [files[0], "--site", str(ROOT / "viewer"), "--bathy", "none", "--viewer-only"])

    # 4. Record fingerprints and regions; optionally prune
    index = json.loads(idx_path.read_text())
    by_file = {Path(f).name: f for f in files}
    for e in index:
        if e.get("file") in by_file:
            e["fingerprint"] = fps[by_file[e["file"]]]
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
            name=e["name"], region=e.get("region") or "", glider=e.get("glider") or "",
            start=e["start"], end=e["end"], profiles=e["profiles"], maxDepth=e["maxDepth"],
            lat0=e["latRange"][0], lat1=e["latRange"][1], lon0=e["lonRange"][0], lon1=e["lonRange"][1])
        (ROOT / "deployments" / f"{e['name']}.qmd").write_text(page, encoding="utf-8")
    groups = OrderedDict()
    for e in index:
        groups.setdefault(e.get("region") or "Other", []).append(e)
    parts = []
    for region, es in groups.items():
        rows = "\n".join(
            f"| [{e['name']}](deployments/{e['name']}.qmd) | {e.get('glider', '')} | {e['start']} | {e['end']} "
            f"| {e['profiles']:,} | {e['maxDepth']:.0f} m |" for e in es)
        heading = f"## {region} ({len(es)})\n\n" if len(groups) > 1 else ""
        parts.append(f"{heading}| Deployment | Glider | Start | End | Profiles | Max depth |\n"
                     f"|---|---|---|---|---|---|\n{rows}\n")
    (ROOT / "deployments.qmd").write_text(
        "---\ntitle: \"Deployments\"\n---\n\n"
        f"{len(index)} deployments. Click a name for its 3D curtain, depth-time section and T–S diagram.\n\n"
        + "\n".join(parts), encoding="utf-8")
    print(f"Wrote {len(index)} deployment pages. Next: commit and push (or open a pull request).")


if __name__ == "__main__":
    main()
