#!/usr/bin/env python3
"""_status.md: the badge row at the top of the home and Real-time pages.

    python scripts/site_status.py .erddap-cache/ERDDAP currents viewer/data

Run after the ERDDAP, currents and match-up steps. Badges:

* when this build ran,
* the newest glider data on the Glider DAC (and which glider),
* the newest satellite SST and chlorophyll day on the site.

(The last-commit badge is in the page footer, _quarto.yml.) Badges are shields.io static badges written at build time; the line under them gives details
and any problems (DAC unreachable, datasets not updated, satellite source unavailable).
"""
import datetime as dt
import json
import sys
import urllib.parse
from pathlib import Path
from zoneinfo import ZoneInfo

ROOT = Path(__file__).resolve().parent.parent
TZ = ZoneInfo("America/Los_Angeles")


def badge(label, message, color, alt=None, link=None):
    esc = lambda s: urllib.parse.quote(str(s).replace("-", "--").replace("_", "__"), safe="")
    img = f"![{alt or label}](https://img.shields.io/badge/{esc(label)}-{esc(message)}-{color})"
    return f"[{img}]({link})" if link else img


def pacific(t):
    return t.astimezone(TZ).strftime("%b %-d, %Y %H:%M %Z")


def day(d):
    return dt.date.fromisoformat(d).strftime("%b %-d, %Y")


def load(p):
    try:
        return json.loads(Path(p).read_text(encoding="utf-8"))
    except Exception:
        return None


def main(argv=None):
    argv = sys.argv[1:] if argv is None else argv
    erddap_dir, currents_dir, data_dir = (Path(x) for x in argv[:3])
    now = dt.datetime.now(dt.timezone.utc)
    status = load(erddap_dir / "status.json") or {}
    manifest = load(erddap_dir / "manifest.json") or {}
    currents = load(currents_dir / "currents.json") or {}

    badges = [badge("site updated", pacific(now), "blue")]
    notes, problems = [], []

    # newest glider data on the DAC
    newest = None
    for dep, m in manifest.items():
        try:
            t = dt.datetime.fromisoformat(m["maxTime"].replace("Z", "+00:00"))
        except Exception:
            continue
        if newest is None or t > newest[0]:
            newest = (t, dep, m.get("mode", ""))
    if newest:
        age = (now - newest[0]).total_seconds() / 86400
        color = "brightgreen" if age < 1 else "yellow" if age < 3 else "lightgrey"
        badges.append(badge("newest glider data", pacific(newest[0]), color,
                            link="realtime.qmd" if newest[2] == "realtime" else "deployments.qmd"))
    if newest:
        notes.append(f"newest data from {newest[1]} ({'real-time' if newest[2] == 'realtime' else 'delayed-mode'})")
    checked_note = None
    if status.get("checked"):
        try:
            checked = dt.datetime.strptime(status["checked"], "%Y-%m-%d %H:%M UTC").replace(tzinfo=dt.timezone.utc)
            checked_note = f"Glider DAC checked {pacific(checked)}"
        except ValueError:
            checked_note = f"Glider DAC checked {status['checked']}"
    if checked_note:
        notes.insert(0, checked_note)
    notes.append(f"{len(manifest)} deployments from the Glider DAC")
    if status.get("listing_failed"):
        problems.append("the Glider DAC could not be reached; showing the previous data")
    elif status.get("errors"):
        problems.append(f"{len(status['errors'])} Glider DAC dataset(s) could not be updated")
    if status.get("warnings"):
        problems.append(f"{len(status['warnings'])} deployment(s) missing from the Glider DAC listing (kept)")

    # newest satellite days: the real-time map layer and the curtain match-ups
    sst, chl, sst_src = [], [], set()
    for g in currents.get("gliders", []):
        for var, acc in (("sst", sst), ("chl", chl)):
            b = (g.get("sat") or {}).get(var)
            if b:
                acc += [d["date"] for d in b.get("days", [])]
                if var == "sst":
                    sst_src.add("OSTIA" if "OSTIA" in b.get("source", "") else "MUR")
    for p in data_dir.glob("*.sat.json"):
        s = load(p) or {}
        for d in s.get("days", []):
            if d.get("sst") is not None:
                sst.append(d["date"]); sst_src.add("MUR")
            if d.get("satChl") is not None:
                chl.append(d["date"])
    if sst:
        badges.append(badge("satellite SST", "through " + day(max(sst)), "blue"))
    if chl:
        badges.append(badge("satellite chlorophyll", "through " + day(max(chl)), "blue"))
    if sst or chl:
        notes.append("satellite: " + " and ".join(sorted(sst_src) or ["—"]) + " SST, Copernicus-GlobColour chlorophyll")
    if "OSTIA" in sst_src and currents.get("gliders"):
        problems.append("MUR SST was not available for the real-time map; it shows OSTIA")
    if problems:
        badges.append(badge("problems", str(len(problems)), "orange"))

    text = " ".join(badges) + "\n\n*" + " · ".join(notes) + "*"
    if problems:
        text += "\n\n*⚠ " + "; ".join(problems) + ".*"
    (ROOT / "_status.md").write_text(text + "\n", encoding="utf-8")
    print(text)
    return 0


if __name__ == "__main__":
    try:
        sys.exit(main())
    except Exception as exc:                       # keep the status line from build_site.py
        print(f"site status failed ({exc}); keeping the previous _status.md")
        sys.exit(0)
