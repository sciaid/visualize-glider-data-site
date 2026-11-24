#!/usr/bin/env python3
"""
echo_plot.py - Echograms (time x depth) of calibrated Sv at 70 and 120 kHz and the 120 - 70 kHz
difference, from an echo file written by echo_ad2.py.

    python scripts/echo_plot.py echo/amlr08-20261119_echo.nc [-o FILE.png] [options]

Panels: Sv 70 kHz, Sv 120 kHz (same colour scale), Sv(120) - Sv(70) in every bin, and the same
difference only where both frequencies exceed --mask-sv (default -90 dB), i.e. away from the
noise floor. Each profile is drawn as a column up to --col-hours wide (default 3 h); longer gaps
between profiles stay blank.

Requirements: pip install numpy netCDF4 matplotlib
"""
from __future__ import annotations

import argparse
import datetime as dt
from pathlib import Path

import numpy as np
import netCDF4
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from matplotlib.colors import LinearSegmentedColormap

# Sequential (one hue, light -> dark) and diverging (blue <-> red, grey midpoint) ramps
SEQ = LinearSegmentedColormap.from_list("seq_blue", ["#cde2fb", "#86b6ef", "#3987e5", "#1c5cab", "#0d366b"])
DIV = LinearSegmentedColormap.from_list("div", ["#104281", "#3987e5", "#9ec5f4", "#f0efec",
                                               "#f4aeac", "#e34948", "#7d1d1d"])
INK, INK2, MUTED, SURF, BLANK = "#0b0b0b", "#52514e", "#7a7873", "#fcfcfb", "#f3f2ef"


def columns(t, max_halfwidth):
    """Left/right edges for each profile column: halfway to neighbours, capped."""
    left, right = np.empty_like(t), np.empty_like(t)
    for i in range(t.size):
        lo = (t[i] - t[i - 1]) / 2 if i else max_halfwidth
        hi = (t[i + 1] - t[i]) / 2 if i < t.size - 1 else max_halfwidth
        left[i] = t[i] - min(lo, max_halfwidth)
        right[i] = t[i] + min(hi, max_halfwidth)
    return left, right


def panel(ax, t, edges_t, z_edges, data, cmap, vmin, vmax):
    left, right = edges_t
    to_num = lambda s: mdates.date2num(dt.datetime.fromtimestamp(s, dt.timezone.utc))
    im = None
    for i in range(t.size):
        x = [to_num(left[i]), to_num(right[i])]
        im = ax.pcolormesh(x, z_edges, data[:, i:i + 1], cmap=cmap, vmin=vmin, vmax=vmax,
                           shading="flat", rasterized=True)
    ax.set_facecolor(BLANK)
    ax.set_ylim(z_edges[-1], 0)
    ax.set_ylabel("Depth (m)", color=INK2)
    return im


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("echo_nc")
    ap.add_argument("-o", "--out", help="PNG path (default: next to the .nc file)")
    ap.add_argument("--vmin", type=float, default=-115.0, help="Sv colour scale minimum (dB)")
    ap.add_argument("--vmax", type=float, default=-80.0, help="Sv colour scale maximum (dB)")
    ap.add_argument("--dmax", type=float, default=15.0, help="difference colour scale +/- (dB)")
    ap.add_argument("--mask-sv", type=float, default=-90.0,
                    help="last panel: show the difference only where both Sv exceed this (dB)")
    ap.add_argument("--col-hours", type=float, default=3.0, help="maximum column width per profile (h)")
    a = ap.parse_args(argv)

    ds = netCDF4.Dataset(a.echo_nc)
    t = ds["time"][:].astype(float)
    depth = ds["depth"][:].astype(float)
    z_edges = ds["depth_bnds"][:, 0].tolist() + [float(ds["depth_bnds"][-1, 1])]
    s70 = np.ma.filled(ds["sv_70khz"][:].astype(float), np.nan)
    s120 = np.ma.filled(ds["sv_120khz"][:].astype(float), np.nan)
    name = getattr(ds, "deployment_name", Path(a.echo_nc).stem)
    cal = getattr(ds, "calibration_date", "?")
    ds.close()

    delta = s120 - s70
    masked = np.where((s70 > a.mask_sv) & (s120 > a.mask_sv), delta, np.nan)
    edges_t = columns(t, a.col_hours * 3600 / 2)

    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 9, "axes.edgecolor": "#c9c7c1",
                         "axes.labelcolor": INK2, "xtick.color": INK2, "ytick.color": INK2})
    fig, axes = plt.subplots(4, 1, figsize=(11, 10.5), sharex=True, facecolor=SURF,
                             gridspec_kw={"hspace": 0.32})
    specs = [
        (s70, SEQ, a.vmin, a.vmax, "Sv 70 kHz", "Sv (dB re 1 m⁻¹)"),
        (s120, SEQ, a.vmin, a.vmax, "Sv 120 kHz", "Sv (dB re 1 m⁻¹)"),
        (delta, DIV, -a.dmax, a.dmax, "Sv(120) − Sv(70), every bin", "ΔSv (dB)"),
        (masked, DIV, -a.dmax, a.dmax, f"Sv(120) − Sv(70), only where both Sv > {a.mask_sv:g} dB", "ΔSv (dB)"),
    ]
    for ax, (data, cmap, lo, hi, title, cbl) in zip(axes, specs):
        im = panel(ax, t, edges_t, z_edges, data, cmap, lo, hi)
        n = int(np.isfinite(data).sum())
        ax.set_title(title, loc="left", fontsize=10, color=INK, fontweight="bold")
        if "only where" in title:
            ax.text(1.0, 1.02, f"{n} of {data.size} bins", transform=ax.transAxes, ha="right",
                    va="bottom", fontsize=8.5, color=MUTED)
        cb = fig.colorbar(im, ax=ax, pad=0.01, fraction=0.025, extend="both")
        cb.set_label(cbl, color=INK2)
        cb.outline.set_visible(False)
        for s in ("top", "right"):
            ax.spines[s].set_visible(False)
    axes[-1].xaxis.set_major_locator(mdates.DayLocator())
    axes[-1].xaxis.set_major_formatter(mdates.DateFormatter("%b %d"))
    axes[-1].xaxis.set_minor_locator(mdates.HourLocator(byhour=[6, 12, 18]))
    axes[-1].set_xlabel("Date (UTC)", color=INK2)
    t0 = dt.datetime.fromtimestamp(t[0], dt.timezone.utc)
    t1 = dt.datetime.fromtimestamp(t[-1], dt.timezone.utc)
    fig.suptitle(f"{name} · Nortek Signature 100 echo telemetry · {t.size} descent profiles, "
                 f"{t0:%Y-%m-%d} to {t1:%Y-%m-%d}", x=0.06, ha="left", y=0.955, fontsize=11.5, color=INK)
    fig.text(0.06, 0.925, f"Calibrated Sv (Echoview offsets, calibration {cal}); 5 m bins; "
             f"blank = no profile or no data (top 10 m blanked by the instrument).",
             fontsize=8.5, color=MUTED)
    out = Path(a.out) if a.out else Path(a.echo_nc).with_suffix(".png")
    fig.savefig(out, dpi=150, facecolor=SURF, bbox_inches="tight")
    print(f"{out}")


if __name__ == "__main__":
    main()
