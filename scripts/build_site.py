#!/usr/bin/env python3
"""
build_site.py - Rebuild the data behind the website from the glider NetCDF files.

Run this on the computer that has the .nc files, then commit and push:

    python scripts/build_site.py "C:\\Claude_work\\PROJECTS\\Glider L-3 ncfiles"

It writes (relative to the repository root):
    map/index.html                 overview map of all deployments
    viewer/index.html              shared 3D curtain viewer (+ plotly.min.js)
    viewer/data/<name>.json        one data file per deployment
    viewer/data/index.json         list of deployments
    deployments/<name>.qmd         one page per deployment
    deployments.qmd                table of all deployments

A GEBCO/ETOPO grid (*gebco*.nc) in the NetCDF folder, or one subfolder down, is
picked up automatically for seafloor and land.

Requirements: pip install netCDF4 numpy contourpy plotly
"""
import argparse
import json
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
sys.path.insert(0, str(HERE))

import glider_curtain  # noqa: E402
import glider_tracks   # noqa: E402

DEPLOYMENT_PAGE = """---
title: "{name}"
page-layout: full
toc: false
---

| Glider | Start | End | Profiles | Max depth | Area |
|---|---|---|---|---|---|
| {glider} | {start} | {end} | {profiles:,} | {maxDepth:.0f} m | {lat0:.2f}° to {lat1:.2f}°, {lon0:.2f}° to {lon1:.2f}° |

<iframe src="../viewer/index.html?d={name}" title="3D curtain viewer: {name}"
        style="width:100%; height:86vh; border:1px solid #e2e0db; border-radius:6px;"></iframe>

[Open the viewer full screen](../viewer/index.html?d={name}){{target="_blank"}} ·
[How these plots are made](../methods.qmd)
"""


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("nc_dir", help="folder with the glider *.nc files")
    ap.add_argument("--pattern", default="*.nc", help="glob for glider files (default *.nc)")
    a = ap.parse_args()

    nc_dir = Path(a.nc_dir)
    files = sorted(str(p) for p in nc_dir.glob(a.pattern)
                   if not any(k in p.name.lower() for k in ("gebco", "etopo", "bathy")))
    if not files:
        sys.exit(f"No NetCDF files matching {a.pattern} in {nc_dir}")
    print(f"{len(files)} glider files")

    # 1. Overview map, with links to the deployment pages
    print("\n== Overview map ==")
    glider_tracks.build(files, ROOT / "map", "measurement", "auto", 2000, 1,
                        "Glider deployments", link_template="../deployments/{name}.html",
                        map_name="index.html")
    for extra in ("glider_positions.csv", "deployments_summary.csv"):   # keep the repo lean
        (ROOT / "map" / extra).unlink(missing_ok=True)

    # 2. Shared curtain viewer + one data file per deployment
    print("\n== Curtain viewer data ==")
    glider_curtain.main(files + ["--site", str(ROOT / "viewer")])

    # 3. One Quarto page per deployment + the deployments table
    print("\n== Pages ==")
    index = json.loads((ROOT / "viewer" / "data" / "index.json").read_text())
    (ROOT / "deployments").mkdir(exist_ok=True)
    for e in index:
        page = DEPLOYMENT_PAGE.format(
            name=e["name"], glider=e.get("glider") or "", start=e["start"], end=e["end"],
            profiles=e["profiles"], maxDepth=e["maxDepth"],
            lat0=e["latRange"][0], lat1=e["latRange"][1], lon0=e["lonRange"][0], lon1=e["lonRange"][1])
        (ROOT / "deployments" / f"{e['name']}.qmd").write_text(page, encoding="utf-8")
    rows = "\n".join(
        f"| [{e['name']}](deployments/{e['name']}.qmd) | {e.get('glider', '')} | {e['start']} | {e['end']} "
        f"| {e['profiles']:,} | {e['maxDepth']:.0f} m |" for e in index)
    (ROOT / "deployments.qmd").write_text(f"""---
title: "Deployments"
---

{len(index)} deployments. Click a name for its 3D curtain, depth-time section and T–S diagram.

| Deployment | Glider | Start | End | Profiles | Max depth |
|---|---|---|---|---|---|
{rows}
""", encoding="utf-8")
    print(f"Wrote {len(index)} deployment pages. Next: git add -A, commit and push.")


if __name__ == "__main__":
    main()
