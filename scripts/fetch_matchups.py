#!/usr/bin/env python3
"""Daily satellite vs glider match-ups for every deployment, for the strip above the curtain.

    python scripts/fetch_matchups.py viewer/data .erddap-cache/matchups

For each deployment's viewer/data/<name>.json this writes viewer/data/<name>.sat.json:

* the glider's daily near-surface values (fetch_satellite.glider_daily: 0-5 m temperature, 0-10 m
  night chlorophyll, and the day's median position), and
* the satellite value there that day:
  - SST: MUR v4.1 from the CoastWatch ERDDAP (median of the pixels within 5 km, with MUR's own
    analysis error). MUR covers 2002 to now.
  - Chlorophyll: observed pixels only (Copernicus-GlobColour L3, 4 km), within 5 km of the day's
    night position: the reprocessed multi-sensor series first, then the near-real-time
    multi-sensor series (it only keeps recent weeks), then the near-real-time OLCI series (from
    2023). A day with no clear pixel is "no clear view"; a day no product covers yet is retried on
    later builds (the reprocessed series is extended every few months).

Satellite values are cached in the cache folder, so each day is fetched once (days in the last
10 days are fetched again, because the near-real-time values are replaced). The work per build is
capped (BUDGET_S); a backlog is finished over the next builds. Never fails the build.
"""
import base64
import concurrent.futures as cf
import csv
import datetime as dt
import io
import json
import sys
import time
import traceback
import urllib.error
import urllib.parse
import urllib.request
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import fetch_satellite as sat                                     # noqa: E402

VERSION = 1
BUDGET_S = 15 * 60
RECENT_DAYS = 10
PAD = 0.07                       # degrees around the day's position (> 5 km)
CHL_PRODUCTS = [
    ("cmems_obs-oc_glo_bgc-plankton_my_l3-multi-4km_P1D", "GlobColour multi-sensor L3, reprocessed"),
    ("cmems_obs-oc_glo_bgc-plankton_nrt_l3-multi-4km_P1D", "GlobColour multi-sensor L3, near real time"),
    ("cmems_obs-oc_glo_bgc-plankton_nrt_l3-olci-4km_P1D", "GlobColour OLCI L3, near real time"),
]
T0 = time.time()


def out_of_time():
    return time.time() - T0 > BUDGET_S


def f32(b64):
    return np.frombuffer(base64.b64decode(b64), "<f4").astype("f8")


def today():
    return dt.datetime.now(dt.timezone.utc).date()


def is_recent(d):
    return (today() - dt.date.fromisoformat(d)).days <= RECENT_DAYS


# --------------------------------------------------------------------------- #
# Glider
# --------------------------------------------------------------------------- #
def glider_days(D):
    var = {v["name"]: v for v in D["vars"]}
    nd, npf = D["nDepth"], D["nProf"]
    get = lambda n: f32(var[n]["data"]).reshape(nd, npf) if n in var else None
    sun = (D.get("aux") or {}).get("sun")
    sun = np.array([np.nan if x is None else x for x in sun], "f8") if sun and len(sun) == npf else None
    return sat.glider_daily(D["time"], D["lat"], D["lon"], D["depth"],
                            get("temperature"), get("chlorophyll"), sun=sun)


# --------------------------------------------------------------------------- #
# SST: MUR from ERDDAP, one small request per day
# --------------------------------------------------------------------------- #
class Blocked(Exception):
    pass


def mur_day(d, la, lo):
    dlon = PAD / max(np.cos(np.radians(la)), 0.2)
    box = f"[({d}T09:00:00Z)][({la - PAD:.3f}):1:({la + PAD:.3f})][({lo - dlon:.3f}):1:({lo + dlon:.3f})]"
    url = f"{sat.MUR_ERDDAP}.csv?analysed_sst{box},analysis_error{box}"
    req = urllib.request.Request(url, headers={"User-Agent": sat.USER_AGENT})
    for k in range(3):
        try:
            with urllib.request.urlopen(req, timeout=120) as r:
                rows = list(csv.reader(io.StringIO(r.read().decode())))
            break
        except urllib.error.HTTPError as e:
            if e.code == 403:
                raise Blocked("CoastWatch ERDDAP answered 403") from None
            if e.code == 404:          # no MUR for that day ("Your query produced no matching results")
                return {"status": "none"}
            if k == 2:
                raise
        except Exception:
            if k == 2:
                raise
        time.sleep(5 * (k + 1))
    head, body = rows[0], rows[2:]
    ix = {h: i for i, h in enumerate(head)}
    lat = np.array([float(r[ix["latitude"]]) for r in body])
    lon = np.array([float(r[ix["longitude"]]) for r in body])
    v = np.array([float(r[ix["analysed_sst"]]) if r[ix["analysed_sst"]] not in ("", "NaN") else np.nan for r in body])
    e = np.array([float(r[ix["analysis_error"]]) if r[ix["analysis_error"]] not in ("", "NaN") else np.nan for r in body])
    dist = np.hypot((lat - la) * 111.2, (lon - lo) * 111.2 * np.cos(np.radians(la)))
    m = (dist <= sat.MATCH_KM) & np.isfinite(v)
    if not m.any():
        return {"status": "none"}
    return {"status": "ok", "v": round(float(np.median(v[m])), 2), "err": round(float(np.nanmedian(e[m])), 2),
            "n": int(m.sum()), "t": body[0][ix["time"]][:10]}


# --------------------------------------------------------------------------- #
# Chlorophyll: observed L3 pixels, via the Copernicus Marine toolbox (lazy, no NetCDF files)
# --------------------------------------------------------------------------- #
_open = {}


def chl_dataset(did):
    if did not in _open:
        import copernicusmarine
        try:
            ds = copernicusmarine.open_dataset(dataset_id=did, variables=["CHL"])
            times = ds["time"].values.astype("datetime64[D]")
            _open[did] = (ds, set(str(x) for x in times))
            print(f"  {did}: {str(times.min())} to {str(times.max())}", flush=True)
        except Exception as exc:
            print(f"  {did}: not available ({type(exc).__name__}: {exc})", flush=True)
            _open[did] = (None, set())
    return _open[did]


def chl_days(want):
    """want: {date: (lat, lon)} -> {date: result}. Dates no product covers are left out."""
    out = {}
    left = dict(want)
    for did, label in CHL_PRODUCTS:
        if not left or out_of_time():
            break
        ds, have = chl_dataset(did)
        todo = sorted(d for d in left if d in have)
        if ds is None or not todo:
            continue
        lats = np.array([left[d][0] for d in todo]); lons = np.array([left[d][1] for d in todo])
        dl = PAD / max(np.cos(np.radians(np.abs(lats).max())), 0.2)
        la0, la1, lo0, lo1 = lats.min() - PAD, lats.max() + PAD, lons.min() - dl, lons.max() + dl
        lat_c, lon_c = ds["latitude"].values, ds["longitude"].values
        ys = slice(la1, la0) if lat_c[0] > lat_c[-1] else slice(la0, la1)
        sub = ds["CHL"].sel(latitude=ys, longitude=slice(lo0, lo1))
        lat, lon = sub["latitude"].values, sub["longitude"].values
        for d in todo:
            if out_of_time():
                break
            la, lo = left[d]
            a = np.asarray(sub.sel(time=np.datetime64(d)).values, "f8")
            m = sat.near(lat, lon, la, lo) & np.isfinite(a)
            out[d] = ({"status": "ok", "v": round(float(np.median(a[m])), 3), "n": int(m.sum()), "product": label}
                      if m.any() else {"status": "cloud", "product": label})
            del left[d]
    return out


# --------------------------------------------------------------------------- #
def main(argv=None):
    argv = sys.argv[1:] if argv is None else argv
    data_dir, cache_dir = Path(argv[0]), Path(argv[1])
    cache_dir.mkdir(parents=True, exist_ok=True)
    files = sorted(p for p in data_dir.glob("*.json") if p.name != "index.json" and not p.name.endswith(".sat.json"))
    # most recent deployments first, so real-time ones are always done within the budget
    deps = []
    for p in files:
        try:
            D = json.loads(p.read_text(encoding="utf-8"))
            deps.append((max(D["time"]), p, D))
        except Exception as exc:
            print(f"{p.name}: cannot read ({exc})")
    deps.sort(key=lambda x: -x[0])
    have_cm = True
    try:
        import copernicusmarine  # noqa: F401
    except ImportError:
        have_cm = False
    mur_blocked = False
    totals = {"sst": 0, "chl": 0, "pending": 0}

    for _, p, D in deps:
        name = D["name"]
        cpath = cache_dir / f"{name}.json"
        cache = json.loads(cpath.read_text()) if cpath.exists() else {}
        if cache.get("version") != VERSION:
            cache = {"version": VERSION, "sst": {}, "chl": {}}
        g = glider_days(D)
        pos = lambda r, night: (r.get("latNight", r["lat"]), r.get("lonNight", r["lon"])) if night else (r["lat"], r["lon"])
        key = lambda la, lo: f"{la:.2f},{lo:.2f}"

        def needed(var, night):
            out = {}
            for d, r in g.items():
                la, lo = pos(r, night)
                c = cache[var].get(d)
                if c and c.get("pos") == key(la, lo) and c["status"] != "pending" and not is_recent(d):
                    continue
                if c and c.get("pos") == key(la, lo) and is_recent(d) and c.get("fetched") == str(today()):
                    continue
                out[d] = (la, lo)
            return out

        # SST
        want = needed("sst", False)
        new_sst = 0
        if want and not mur_blocked and not out_of_time():
            with cf.ThreadPoolExecutor(2) as ex:
                futs = {ex.submit(mur_day, d, *want[d]): d for d in sorted(want)}
                for f in cf.as_completed(futs):
                    d = futs[f]
                    try:
                        r = f.result()
                    except Blocked as exc:
                        if not mur_blocked:
                            print(f"  MUR: {exc}; the rest waits for the next build", flush=True)
                        mur_blocked = True
                        for q in futs:
                            q.cancel()
                        continue
                    except Exception as exc:
                        print(f"  {name} {d}: MUR failed ({type(exc).__name__}: {exc})", flush=True)
                        continue
                    cache["sst"][d] = {**r, "pos": key(*want[d]), "fetched": str(today())}
                    new_sst += 1
                    if out_of_time():
                        for q in futs:
                            q.cancel()
        # chlorophyll (night positions)
        want = {d: v for d, v in needed("chl", True).items() if g[d].get("nNight")}
        new_chl = 0
        if want and have_cm and not out_of_time():
            try:
                res = chl_days(want)
            except Exception as exc:
                traceback.print_exc()
                print(f"  {name}: chlorophyll failed ({type(exc).__name__})")
                res = {}
            for d, (la, lo) in want.items():
                r = res.get(d, {"status": "pending"})       # no product covers this day yet
                cache["chl"][d] = {**r, "pos": key(la, lo), "fetched": str(today())}
                new_chl += d in res
        cpath.write_text(json.dumps(cache, separators=(",", ":")), encoding="utf-8")

        # the strip's data
        days = []
        for d, r in g.items():
            s, c = cache["sst"].get(d, {}), cache["chl"].get(d, {})
            row = {"date": d, "T": r.get("T"), "Tdepth": r.get("Tdepth"), "chl": r.get("chl"),
                   "nProfiles": r["nProfiles"], "nNight": r["nNight"]}
            if s.get("status") == "ok":
                row.update(sst=s["v"], sstErr=s.get("err"))
            row["sstStatus"] = s.get("status", "pending")
            if c.get("status") == "ok":
                row.update(satChl=c["v"], chlProduct=c.get("product"))
            row["chlStatus"] = c.get("status", "pending" if r["nNight"] else "no night profiles")
            days.append(row)
        n_sst = sum(1 for x in days if "sst" in x)
        n_chl = sum(1 for x in days if "satChl" in x)
        pend = sum(1 for x in days if x["sstStatus"] == "pending") + sum(1 for x in days if x["chlStatus"] == "pending")
        totals["sst"] += new_sst; totals["chl"] += new_chl; totals["pending"] += pend
        out = {"version": VERSION, "made": dt.datetime.now(dt.timezone.utc).strftime("%Y-%m-%d %H:%M UTC"),
               "sources": {"sst": {"name": sat.MUR_NAME, "url": sat.MUR_URL},
                           "chl": {"name": "Copernicus-GlobColour chlorophyll, observed pixels only (L3, 4 km daily)",
                                   "url": sat.CHL_DOI}},
               "nightSunElev": sat.NIGHT_SUN_ELEV, "matchKm": sat.MATCH_KM, "days": days}
        (data_dir / f"{name}.sat.json").write_text(json.dumps(out, separators=(",", ":")), encoding="utf-8")
        print(f"{name}: {len(days)} days; SST {n_sst} matched (+{new_sst} new), chlorophyll {n_chl} matched "
              f"(+{new_chl} new), {sum(1 for x in days if x['chlStatus'] == 'cloud')} no clear view, {pend} pending",
              flush=True)
    print(f"\nDone in {time.time() - T0:.0f} s: {totals['sst']} new SST and {totals['chl']} new chlorophyll "
          f"values; {totals['pending']} day-values still pending" + ("; MUR was blocked" if mur_blocked else ""))
    return 0


if __name__ == "__main__":
    try:
        sys.exit(main())
    except Exception:
        traceback.print_exc()
        print("match-up step failed; the curtains show no satellite strip")
        sys.exit(0)
