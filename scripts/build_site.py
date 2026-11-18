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
import html
import json
import sys
from collections import OrderedDict
from pathlib import Path
from urllib.parse import quote

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




def load_project_names():
    """projects.json: [{"starts_with": <start of the file's project attribute>, "short": <label>}]."""
    p = ROOT / "projects.json"
    return json.loads(p.read_text(encoding="utf-8")) if p.exists() else []


PROJECT_NAMES = load_project_names()


def short_project(full):
    """Short label for a project: first matching entry in projects.json (case-insensitive
    'starts with'), otherwise the text before any ':' ("FREEBYRD: Estimate ..." -> "FREEBYRD")."""
    for entry in PROJECT_NAMES:
        if full.lower().startswith(entry["starts_with"].lower()):
            return entry["short"]
    return full.split(":")[0].strip()


def file_meta(path):
    """Deployment name and project from a file's global attributes (cheap: no data read)."""
    import netCDF4
    with netCDF4.Dataset(path) as ds:
        name = str(getattr(ds, "deployment_name", "") or Path(path).stem).strip()
        full = str(getattr(ds, "project", "") or "").strip()
        mode = str(getattr(ds, "data_mode", "") or "").strip()
        src_id = str(getattr(ds, "source_dataset", "") or "").strip()
        src_url = str(getattr(ds, "source_url", "") or "").strip()
        modified = str(getattr(ds, "date_modified", "") or "").strip()[:10] if src_id else ""
    return {"name": name, "project": short_project(full), "projectFull": full,
            "mode": mode, "sourceDataset": src_id, "sourceUrl": src_url, "modified": modified}


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
        return (f"[{e['sourceDataset']}]({e['sourceUrl']}) on the IOOS Glider DAC ERDDAP"
                + (" (real-time)" if e.get("mode") == "realtime" else "")
                + (f", last modified {e['modified']}" if e.get("modified") else ""))
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
        label = f"[{region}](index.html?region={quote(region)})" if len(groups) > 1 else region
        rows.append(f"| {label} | {count} | {yr} | {projects} | ≈ {round(km, -2):,.0f} km |")
    total = sum(e.get("distanceKm") or 0 for e in index)
    if len(groups) > 1:
        rows.append(f"| **All** | **{len(index)}** | | | **≈ {round(total, -2):,.0f} km** |")
    (ROOT / "_summary.md").write_text(
        "| Region | Deployments | Years (start) | Projects | Distance (approx.) |\n"
        "|---|---|---|---|---|\n" + "\n".join(rows) + "\n", encoding="utf-8")


REGIONBAR = """```{=html}
<div id="regionbar" class="regionbar" data-default="__DEFAULT__" __HIDE__>
  <span class="rb-label">Region</span>
__BUTTONS__
</div>
<iframe id="mapframe" title="Map of glider deployments"
        style="width:100%; height:72vh; border:1px solid #e2e0db; border-radius:6px;"></iframe>
<style>
  .regionbar { display:flex; flex-wrap:wrap; gap:8px; align-items:stretch; margin:0 0 10px; }
  .regionbar .rb-label { align-self:center; font-weight:600; margin-right:4px; }
  .regionbar button { text-align:left; border:1px solid #d6d3cc; background:#fff; border-radius:8px;
                      padding:6px 12px; cursor:pointer; line-height:1.25; font:inherit; color:inherit; }
  .regionbar button:hover { border-color:#2a78d6; }
  .regionbar button[aria-pressed="true"] { border-color:#2a78d6; background:#eaf2fc; box-shadow:inset 0 0 0 1px #2a78d6; }
  .regionbar .rb-name { display:block; font-weight:600; }
  .regionbar .rb-meta { display:block; font-size:.82em; color:#5d5b56; }
  .regionbar .rb-rt { color:#b4501a; font-weight:600; }
</style>
<script>
(function () {
  const bar = document.getElementById('regionbar'), frame = document.getElementById('mapframe');
  const buttons = [...bar.querySelectorAll('button')];
  const asked = new URLSearchParams(location.search).get('region');
  const known = buttons.map(b => b.dataset.region);
  let region = known.includes(asked) ? asked : bar.dataset.default;
  function mark(r) {
    region = r;
    buttons.forEach(b => b.setAttribute('aria-pressed', b.dataset.region === r));
    const u = new URL(location.href);
    if (r === bar.dataset.default) u.searchParams.delete('region'); else u.searchParams.set('region', r);
    history.replaceState(null, '', u);
  }
  frame.src = 'map/index.html?regionbar=' + (bar.hidden ? '0' : '1') + '&region=' + encodeURIComponent(region);
  mark(region);
  buttons.forEach(b => b.onclick = () => {
    mark(b.dataset.region);
    try { frame.contentWindow.setRegion(b.dataset.region); }
    catch (e) { frame.src = 'map/index.html?regionbar=1&region=' + encodeURIComponent(b.dataset.region); }
  });
  window.addEventListener('message', e => {           // keep the buttons in step with the map
    if (e.source === frame.contentWindow && e.data && e.data.glidermapRegion) mark(e.data.glidermapRegion);
  });
})();
</script>
```
"""


def write_regionbar(index):
    """_regionbar.md: region buttons above the home-page map (plus the map iframe itself).

    The map opens on the region whose most recent data are newest; ?region=<name> or
    ?region=All in the page URL overrides that. With a single region the bar is hidden.
    """
    groups = {}
    for e in index:
        groups.setdefault(e.get("region") or "Other", []).append(e)
    latest = max(groups, key=lambda r: max(e["end"] for e in groups[r])) if groups else ""
    def button(name, es, label=None):
        years = sorted({int(e["start"][:4]) for e in es})
        yr = f"{years[0]}–{years[-1]}" if len(years) > 1 else (f"{years[0]}" if years else "")
        nrt = sum(e.get("mode") == "realtime" for e in es)
        rt = f' · <span class="rb-rt">{nrt} real-time</span>' if nrt else ""
        n = len(es)
        return (f'  <button type="button" data-region="{html.escape(name)}" aria-pressed="false">'
                f'<span class="rb-name">{html.escape(label or name)}</span>'
                f'<span class="rb-meta">{n} deployment{"s" if n != 1 else ""} · {yr}{rt}</span></button>')
    order = sorted(groups, key=lambda r: max(e["end"] for e in groups[r]), reverse=True)
    btns = [button(r, groups[r]) for r in order]
    if len(groups) > 1:
        btns.append(button("All", index, "All regions"))
    text = (REGIONBAR.replace("__DEFAULT__", html.escape(latest))
            .replace("__HIDE__", "" if len(groups) > 1 else "hidden")
            .replace("__BUTTONS__", "\n".join(btns)))
    (ROOT / "_regionbar.md").write_text(text, encoding="utf-8")


def walk_nc(folder):
    """Every *.nc under folder, following linked folders too (sorted)."""
    import os
    out = []
    for dirpath, _, names in os.walk(folder, followlinks=True):
        out += [Path(dirpath) / n for n in names if n.lower().endswith(".nc")]
    return sorted(out)


def write_realtime(index, data_dir):
    """realtime.qmd: deployments shown from real-time data, with latest position and data age."""
    import datetime as dt
    now = dt.datetime.now(dt.timezone.utc)
    rows = []
    for e in sorted((e for e in index if e.get("mode") == "realtime"), key=lambda e: e["end"], reverse=True):
        try:
            d = json.loads((data_dir / f"{e['name']}.json").read_text(encoding="utf-8"))
            last_t = dt.datetime.fromtimestamp(d["time"][-1], dt.timezone.utc)
            la, lo = d["lat"][-1], d["lon"][-1]
            pos = f"{abs(la):.3f}°{'N' if la >= 0 else 'S'}, {abs(lo):.3f}°{'E' if lo >= 0 else 'W'}"
            age = (now - last_t).total_seconds() / 86400
            through, ago = last_t.strftime("%Y-%m-%d %H:%M UTC"), f"{age:.1f}"
        except (OSError, KeyError, IndexError, ValueError):
            pos, through, ago = "", e["end"], ""
        rows.append(f"| [{e['name']}](deployments/{e['name']}.qmd) | {e.get('region', '')} | {e.get('project', '')} "
                    f"| {e['start']} | {through} | {ago} | {pos} | {e['profiles']:,} |")
    body = ("| Deployment | Region | Project | Start | Data through | Days since last data | Latest position | Profiles |\n"
            "|---|---|---|---|---|---|---|---|\n" + "\n".join(rows)) if rows else \
        "No deployments are currently shown from real-time data."
    (ROOT / "realtime.qmd").write_text(f"""---
title: "Real-time deployments"
page-layout: full
---

Deployments for which only **real-time (provisional)** data are published on the
[IOOS Glider DAC ERDDAP](https://gliders.ioos.us/erddap/). Real-time data are thinned, not
quality-controlled and may have gaps; ERDDAP is updated every few days while a glider is out.
Each deployment moves to delayed-mode data automatically when that is published.

{{{{< include _status.md >}}}}

{body}
""", encoding="utf-8")


def write_status(status, index):
    """_status.md: when ERDDAP was last checked (shown on the home and real-time pages)."""
    n = sum(1 for e in index if e.get("sourceDataset"))
    msg = f"*ERDDAP last checked: {status.get('checked', 'unknown')} · {n} deployments from ERDDAP"
    if status.get("listing_failed"):
        msg += " · ERDDAP could not be reached; showing the previous data"
    elif status.get("errors"):
        msg += f" · {len(status['errors'])} dataset(s) could not be updated"
    (ROOT / "_status.md").write_text(msg + "*\n", encoding="utf-8")


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
    ap.add_argument("--erddap-status", default=None,
                    help="status.json written by fetch_erddap.py: shown as 'ERDDAP last checked' on the site")
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
        print(f"No glider NetCDF files found in {a.nc_dirs}; keeping the deployments already on the site.")
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
    # Keep tracks of deployments already on the site whose files are not here (e.g. local
    # Antarctic files during a GitHub build, or ERDDAP deployments during a local build).
    present = {m["name"] for m in meta.values()}
    on_site = {e["name"] for e in index}
    previous = [d for d in glider_tracks.read_map_tracks(ROOT / "map" / "index.html")
                if d["name"] in on_site and d["name"] not in present]
    if files:
      glider_tracks.build(files, ROOT / "map", "measurement", "auto", 2000, 1, "Glider deployments",
                          link_template="../deployments/{name}.html", map_name="index.html", regions=regions,
                          extra={m["name"]: m for m in meta.values()}, keep_tracks=previous)
    for extra in ("glider_positions.csv", "deployments_summary.csv"):
        (ROOT / "map" / extra).unlink(missing_ok=True)

    # 3. Curtain data for new/changed deployments (bathymetry grids found next to all files)
    print("\n== Curtain viewer data ==")
    grids = sorted({str(p) for d in a.nc_dirs for p in walk_nc(d)
                    if any(k in p.name.lower() for k in BATHY_KEYS)})
    bathy_arg = ",".join(grids) or "none"
    glider_curtain.main(todo + ["--site", str(ROOT / "viewer"), "--bathy", bathy_arg] if todo else
                        ["none", "--site", str(ROOT / "viewer"), "--bathy", "none", "--viewer-only"])

    # 4. Record fingerprints and regions; optionally prune
    index = json.loads(idx_path.read_text())
    by_file = {Path(f).name: f for f in files}
    for e in index:
        if e.get("file") in by_file:
            f = by_file[e["file"]]
            e["fingerprint"] = fps[f]
            for k in ("project", "projectFull", "mode", "sourceDataset", "sourceUrl", "modified"):
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
    write_regionbar(index)
    write_realtime(index, data_dir)
    if a.erddap_status and Path(a.erddap_status).exists():
        write_status(json.loads(Path(a.erddap_status).read_text(encoding="utf-8")), index)
    print(f"Wrote {len(index)} deployment pages. Next: commit and push (or open a pull request).")


if __name__ == "__main__":
    main()
