#!/usr/bin/env python3
"""
echo_ad2.py - Build the echosounder layer for one glider deployment from Nortek Signature 100
compact-echosounder .ad2 files (echo-telemetry records, ID 0x29 / DF=41) sent home by the glider.

    python scripts/echo_ad2.py AD2_FOLDER DEPLOYMENT [--out-dir echo] [--track-csv FILE]
    python scripts/echo_ad2.py --all          # every deployment listed in echo.json

    python scripts/echo_ad2.py "C:\\Claude_work\\PROJECTS\\AMLR26-27_ad2_files\\amlr08" amlr08-20261119

echo.json (repository root) maps deployment names to their .ad2 folders, for --all:
    {"amlr08-20261119": "C:\\Claude_work\\PROJECTS\\AMLR26-27_ad2_files\\amlr08", ...}
A deployment whose folder has no echo data yet is reported and skipped.

Writes echo/<DEPLOYMENT>_echo.nc: one combined file with every descent profile found in the folder
(one copy per descent profile, see dedup(); rerun it on the whole folder each time new files arrive),
containing
  time, depth (5 m bin centres), latitude, longitude, sv_70khz, sv_120khz (calibrated Sv, dB),
  plus the raw instrument Sv for both frequencies.
build_site.py / glider_curtain.py pick the file up and show it as a separate echo layer on the
deployment's curtain.

Calibration: Sv_calibrated = Sv_instrument + offset (Echoview Calibration Assistant offsets, by
instrument serial number; table below). Bins the instrument reports as -100 dB (blanked/no data)
become NaN.

Positions: each echo profile's latitude/longitude is interpolated in time along the glider track
from the deployment's IOOS Glider DAC dataset on ERDDAP (needs internet), or from a CSV with
columns time,latitude,longitude (--track-csv). The site does its own interpolation along the same
track, so the curtain does not depend on these.

Record layout (reverse-engineered; validated bit-for-bit against Nortek MIDAS NetCDF exports):
  data byte 1 = offset to echograms (16); bytes 2-5 uint32 time (s since 1970); bytes 10-11 uint16
  number of 5 m depth bins (60); then N echograms, each float16 frequency (kHz) + nbins float16 Sv.
  Depth bins are labelled 0, 5, ... 295 m by MIDAS (upper edge, [z, z+5)); stored here as centres
  2.5 ... 297.5 m to match the site's 5 m grids.

Requirements: pip install numpy netCDF4
"""
from __future__ import annotations

import argparse
import csv
import datetime as dt
import io
import json
import re
import struct
import sys
import urllib.parse
import urllib.request
from pathlib import Path

import numpy as np

try:
    import netCDF4
except ImportError:  # pragma: no cover
    sys.exit("netCDF4 is required:  pip install netCDF4 numpy")

ROOT = Path(__file__).resolve().parent.parent
ERDDAP = "https://gliders.ioos.us/erddap"

# Calibration table: serial number -> glider, calibration date, {frequency kHz: Sv offset dB}.
# Offsets from Echoview Calibration Assistant. Edit/add units here.
CAL_TABLE = {
    107475: {"glider": "amlr08", "cal_date": "2026-07-27", 70.0: -63.344, 120.0: -64.844},
    105815: {"glider": "amlr06", "cal_date": "2026-07-28", 70.0: -60.405, 120.0: -65.109},
}
FREQS = {"sv_70khz": 70.0, "sv_120khz": 120.0}     # echograms written (narrowband channels)
FREQ_TOL = 0.5                                      # kHz
FILL = -100.0
DZ = 5.0
DUP_WINDOW = 120    # s: records closer than this are copies of one descent profile


# ----------------------------------------------------------------------------- ad2 decoding
def _checksum(b):
    s = 0xB58C
    for i in range(0, len(b) - 1, 2):
        s += b[i] | (b[i + 1] << 8)
    if len(b) % 2:
        s += b[-1] << 8
    return s & 0xFFFF


def _records(path):
    d = Path(path).read_bytes()
    i = 0
    while i + 10 <= len(d):
        if d[i] != 0xA5 or d[i + 1] not in (10, 12):
            i += 1
            continue
        hs, rid = d[i + 1], d[i + 2]
        ds, dck, hck = struct.unpack("<HHH", d[i + 4:i + 10]) if hs == 10 else struct.unpack("<IHH", d[i + 4:i + 12])
        if _checksum(d[i:i + hs - 2]) != hck:
            i += 1
            continue
        body = d[i + hs:i + hs + ds]
        if len(body) == ds and _checksum(body) == dck:
            yield rid, body
        i += hs + ds


def n_bins(sv):
    """Depth bins with data (not the -100 fill), all echograms."""
    return int(np.sum(np.isfinite(sv) & (sv != FILL)))


def dedup(recs, window=DUP_WINDOW):
    """One record per descent profile. The glider can send a profile more than once (for example
    when it re-sends its whole backlog), and a re-sent copy may differ slightly (timestamp a few
    seconds off, a bin more or less). Records within `window` s are one profile; keep the copy
    with the most depth bins with data (ties: the copy received first). recs are in received order.
    Returns (kept, notes); notes describe groups whose copies were not identical."""
    order = sorted(range(len(recs)), key=lambda i: (recs[i]["time"], i))
    groups = []
    for i in order:
        if groups and recs[i]["time"] - recs[groups[-1][-1]]["time"] <= window:
            groups[-1].append(i)
        else:
            groups.append([i])
    kept, notes = [], []
    for g in groups:
        best = max(g, key=lambda i: (n_bins(recs[i]["sv"]), -i))
        k = recs[best]
        kept.append(k)
        if any(recs[i]["time"] != k["time"] or not np.array_equal(recs[i]["sv"], k["sv"]) for i in g):
            others = "; ".join(f'{recs[i]["src"]} {_ts(recs[i]["time"])} ({n_bins(recs[i]["sv"])} bins)'
                               for i in g if i != best)
            notes.append(f'{_ts(k["time"])}: kept {k["src"]} ({n_bins(k["sv"])} bins); not used: {others}')
    return kept, notes


def _ts(t):
    return dt.datetime.fromtimestamp(t, dt.timezone.utc).strftime("%Y-%m-%d %H:%M:%S")


def read_ad2(files):
    """Echo-telemetry profiles and instrument serial numbers from .ad2 files (in received order,
    i.e. sorted file names), one record per descent profile (see dedup)."""
    recs, sns = [], set()
    for f in sorted(files):
        for rid, b in _records(f):
            if rid == 0xA0:
                m = re.search(rb"SN=(\d+)", b)
                if m:
                    sns.add(int(m.group(1)))
            elif rid == 0x29:
                off = b[1]
                t = struct.unpack("<I", b[2:6])[0]
                nb = struct.unpack("<H", b[10:12])[0]
                v = np.frombuffer(b[off:], "<f2").astype(np.float64)
                ne = len(v) // (nb + 1)
                v = v[:ne * (nb + 1)].reshape(ne, nb + 1)
                recs.append({"time": t, "freq": v[:, 0], "sv": v[:, 1:], "src": Path(f).name})
    kept, notes = dedup(recs)
    for n in notes:
        print(f"  differing copies of one profile - {n}")
    return kept, sns


# ----------------------------------------------------------------------------- glider track
def _get(url, timeout=120):
    req = urllib.request.Request(url, headers={"User-Agent": "visualize-glider-data-site echo_ad2.py"})
    with urllib.request.urlopen(req, timeout=timeout) as r:
        return r.read().decode("utf-8")


def erddap_dataset_id(deployment):
    """The DAC dataset for a deployment name like amlr08-20261119 (delayed preferred)."""
    q = urllib.parse.quote(f'"^{re.escape(deployment)}T[0-9]{{4}}(-delayed)?$"', safe="")
    rows = list(csv.reader(io.StringIO(_get(f"{ERDDAP}/tabledap/allDatasets.csv?datasetID&datasetID=~{q}"))))[2:]
    ids = [r[0] for r in rows if r]
    if not ids:
        raise LookupError(f"no ERDDAP dataset found for {deployment}")
    return sorted(ids, key=lambda s: (not s.endswith("-delayed"), s))[0]


def read_track_csv(text):
    rows = list(csv.DictReader(io.StringIO(text)))
    rows = [r for r in rows if r.get("time") and r["time"] != "UTC"]           # ERDDAP units row
    t = np.array([dt.datetime.fromisoformat(r["time"].replace("Z", "+00:00")).timestamp() for r in rows])
    la = np.array([float(r["latitude"] or "nan") for r in rows])
    lo = np.array([float(r["longitude"] or "nan") for r in rows])
    ok = np.isfinite(t) & np.isfinite(la) & np.isfinite(lo)
    o = np.argsort(t[ok])
    return t[ok][o], la[ok][o], lo[ok][o]


def erddap_track(dataset_id):
    url = (f"{ERDDAP}/tabledap/{dataset_id}.csv?time,latitude,longitude"
           f"&distinct()")
    return read_track_csv(_get(url, timeout=300))


def interp_track(t, tt, la, lo, max_gap=6 * 3600):
    """Linear interpolation in time; NaN outside the track by more than max_gap."""
    lat = np.interp(t, tt, la)
    lon = np.interp(t, tt, lo)
    out = (t < tt[0] - max_gap) | (t > tt[-1] + max_gap)
    lat[out] = np.nan
    lon[out] = np.nan
    return lat, lon


# ----------------------------------------------------------------------------- output
def write_echo(path, deployment, profs, sn, cal, lat, lon, track_source):
    nb = profs[0]["sv"].shape[1]
    nt = len(profs)
    freqs = profs[0]["freq"]
    depth = np.arange(nb) * DZ + DZ / 2
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(".tmp.nc")
    ds = netCDF4.Dataset(tmp, "w")
    ds.title = f"{deployment} echosounder telemetry (Nortek Signature 100, calibrated Sv)"
    ds.deployment_name = deployment
    ds.glider_name = deployment.rsplit("-", 1)[0]
    ds.instrument = "Nortek Signature 100 compact echosounder"
    ds.instrument_serial_number = int(sn)
    ds.calibration_date = cal["cal_date"]
    ds.calibration_source = "Echoview Calibration Assistant"
    ds.processing = ("Descent profiles processed on the instrument (GETTMECHO, DF=41), 5 m bins from 0-300 m, "
                     "sent via Iridium in .ad2 files; decoded and calibrated by echo_ad2.py")
    ds.position_source = track_source
    ds.source_files = ", ".join(sorted({p["src"] for p in profs}))
    last = dt.datetime.fromtimestamp(profs[-1]["time"], dt.timezone.utc)
    ds.history = f"echo_ad2.py; {nt} profiles up to {last:%Y-%m-%dT%H:%M:%SZ}"   # data-only: rerun without new data = same file
    ds.createDimension("time", nt)
    ds.createDimension("depth", nb)
    ds.createDimension("nv", 2)
    v = ds.createVariable("time", "f8", ("time",)); v[:] = [p["time"] for p in profs]
    v.units = "seconds since 1970-01-01T00:00:00Z"; v.standard_name = "time"
    v.comment = "instrument timestamp of the descent profile"
    v = ds.createVariable("depth", "f4", ("depth",)); v[:] = depth
    v.units = "m"; v.positive = "down"; v.bounds = "depth_bnds"; v.comment = "bin centre"
    v = ds.createVariable("depth_bnds", "f4", ("depth", "nv")); v[:] = np.c_[depth - DZ / 2, depth + DZ / 2]
    for name, arr, unit in (("latitude", lat, "degrees_north"), ("longitude", lon, "degrees_east")):
        v = ds.createVariable(name, "f8", ("time",), fill_value=np.nan); v[:] = arr
        v.units = unit; v.standard_name = name; v.comment = f"interpolated in time along the glider track ({track_source})"
    for name, f in FREQS.items():
        e = int(np.argmin(np.abs(freqs - f)))
        if abs(freqs[e] - f) > FREQ_TOL:
            continue
        raw = np.stack([p["sv"][e] for p in profs]).T            # depth x time
        raw = np.where(raw == FILL, np.nan, raw)
        off = cal.get(f)
        for vname, data, desc in ((f"{name}_raw", raw, "as reported by the instrument (no calibration)"),
                                  (name, raw + off if off is not None else raw * np.nan,
                                   f"calibrated: instrument Sv + {off} dB")):
            v = ds.createVariable(vname, "f4", ("depth", "time"), fill_value=np.float32(np.nan), zlib=True)
            v[:] = data.astype("f4")
            v.units = "dB re 1 m-1"
            v.long_name = f"Volume backscattering strength {f:g} kHz, {desc}"
            v.frequency_khz = float(freqs[e])
            if vname == name:
                v.calibration_offset_db = float(off) if off is not None else np.nan
                v.calibration_date = cal["cal_date"]
    ds.close()
    tmp.replace(path)


class Skip(Exception):
    """Nothing to write for this deployment (reported, not an error)."""


def build_one(ad2_dir, deployment, out_dir, track_csv=None, dataset_id=None):
    files = sorted(Path(ad2_dir).rglob("*.ad2"))
    if not files:
        raise Skip(f"no .ad2 files in {ad2_dir}")
    profs, sns = read_ad2(files)
    if len(sns) > 1:
        raise SystemExit(f"{ad2_dir} holds files from more than one instrument (SN {sorted(sns)})")
    if not sns:
        raise SystemExit("no instrument serial number in these files (no configuration string)")
    sn = sns.pop()
    if sn not in CAL_TABLE:
        raise SystemExit(f"SN {sn} is not in CAL_TABLE; add its calibration offsets to {Path(__file__).name}")
    cal = CAL_TABLE[sn]
    glider = deployment.rsplit("-", 1)[0]
    if cal["glider"] != glider:
        print(f"  note: SN {sn} is listed for {cal['glider']}, deployment is {glider}")
    print(f"{deployment}: {len(files)} .ad2 files, {len(profs)} echo profiles, SN {sn} "
          f"({cal['glider']}, calibrated {cal['cal_date']})")
    if not profs:
        raise Skip("no echo-telemetry (0x29) records yet; nothing written")

    t = np.array([p["time"] for p in profs], float)
    lat = lon = np.full(t.size, np.nan)
    try:
        if track_csv:
            tt, la, lo = read_track_csv(Path(track_csv).read_text(encoding="utf-8"))
            src = f"track file {Path(track_csv).name}"
        else:
            did = dataset_id or erddap_dataset_id(deployment)
            tt, la, lo = erddap_track(did)
            src = f"{ERDDAP}/tabledap/{did}"
        lat, lon = interp_track(t, tt, la, lo)
        print(f"  positions from {src}: {np.isfinite(lat).sum()} of {t.size} profiles")
    except Exception as e:                                   # offline etc.: the site still places them
        src = "none (track unavailable when the file was built)"
        print(f"  note: no positions ({e}); the site interpolates along its own track")

    out = Path(out_dir) / f"{deployment}_echo.nc"
    before = out.read_bytes() if out.exists() else None
    write_echo(out, deployment, profs, sn, cal, lat, lon, src)
    t0 = dt.datetime.fromtimestamp(t[0], dt.timezone.utc)
    t1 = dt.datetime.fromtimestamp(t[-1], dt.timezone.utc)
    state = "unchanged" if before == out.read_bytes() else ("new" if before is None else "updated")
    print(f"  {len(profs)} profiles {t0:%Y-%m-%d %H:%M} -> {t1:%Y-%m-%d %H:%M} UTC -> {out} ({state})")


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("ad2_dir", nargs="?", help="folder with the deployment's .ad2 files (searched recursively)")
    ap.add_argument("deployment", nargs="?", help="deployment name as on the site, e.g. amlr08-20261119")
    ap.add_argument("--all", action="store_true", help="every deployment in echo.json (repository root)")
    ap.add_argument("--config", default=str(ROOT / "echo.json"), help="default: %(default)s")
    ap.add_argument("--out-dir", default=str(ROOT / "echo"), help="default: <repo>/echo")
    ap.add_argument("--track-csv", help="CSV with time,latitude,longitude instead of ERDDAP")
    ap.add_argument("--dataset-id", help="ERDDAP dataset ID (default: found from the deployment name)")
    a = ap.parse_args(argv)

    if a.all:
        jobs = list(json.loads(Path(a.config).read_text(encoding="utf-8")).items())
        jobs = [(folder, name) for name, folder in jobs]
    elif a.ad2_dir and a.deployment:
        jobs = [(a.ad2_dir, a.deployment)]
    else:
        ap.error("give AD2_FOLDER DEPLOYMENT, or --all")
    for folder, name in jobs:
        try:
            build_one(folder, name, a.out_dir, a.track_csv if not a.all else None,
                      a.dataset_id if not a.all else None)
        except Skip as e:
            print(f"{name}: skipped - {e}")


if __name__ == "__main__":
    main()
