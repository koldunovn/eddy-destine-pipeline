"""Render per-day frames of SSH with overlaid eddy tracks (trailing line
+ current-position dot), suitable for ffmpeg to stitch into a video.

One frame per date. Anticyclonic in firebrick, cyclonic in royalblue.
Each track is drawn only for obs within the last TRAIL_DAYS up to the
frame date; the current position is a dot at the end of the trail.

Run in the eddy env:
    /home/koldunovn/envs/eddy/bin/python plot_eddy_frames.py \
        --source ifs-fesom --start 1992-01-01 --end 1992-01-01

    # full year, 60-day trail
    /home/koldunovn/envs/eddy/bin/python plot_eddy_frames.py \
        --source ifs-fesom --start 1992-01-01 --end 1992-12-31 --trail-days 60
"""

from __future__ import annotations

import argparse
from pathlib import Path

import cartopy.crs as ccrs
import cartopy.feature as cfeature
import matplotlib.pyplot as plt
import numpy as np
import xarray as xr

import eddy_paths
OUT_ROOT = eddy_paths.default_root()

REGIONS = {
    # name: (lon0, lon1, lat0, lat1) in -180..180 / -90..90, or None for global
    "global":     None,
    "gulfstream": (-82, -40, 25, 50),
    "agulhas":    (10, 40, -45, -25),
    "kuroshio":   (130, 170, 25, 45),
}

# Per-region SSH colour-scale defaults (symmetric ±vmax). Override with --vmax.
REGION_VMAX = {
    "global":     1.2,
    "gulfstream": 0.8,
    "agulhas":    1.2,
    "kuroshio":   1.2,
}

SOURCES = {
    "ifs-fesom": {
        "adt_dir": OUT_ROOT / "out/adt_ifs-fesom_hist",
        "adt_pattern": "adt_DestinE_ifs-fesom_hist_{ymd}.nc",
        "tracks_root": OUT_ROOT / "out/tracks_ifs-fesom_hist_1990_2014",
        "label": "DestinE IFS-FESOM hist",
    },
    "icon": {
        "adt_dir": OUT_ROOT / "out/adt_icon_hist",
        "adt_pattern": "adt_DestinE_icon_hist_{ymd}.nc",
        "tracks_root": OUT_ROOT / "out/tracks_icon_hist_1990_2014",
        "label": "DestinE ICON hist",
    },
}


def load_track_groups(tracks_nc: Path, t_min: np.datetime64, t_max: np.datetime64):
    """Return list of (time_arr, lon_arr[-180..180], lat_arr) per track, keeping
    only tracks that have at least one obs in [t_min, t_max]."""
    ds = xr.open_dataset(tracks_nc)
    tid = ds["track"].values
    time = ds["time"].values
    lon = ((ds["longitude"].values + 180) % 360) - 180
    lat = ds["latitude"].values

    # Drop tracks entirely outside window. First identify tracks with any obs
    # inside [t_min, t_max].
    in_win = (time >= t_min) & (time <= t_max)
    live_ids = np.unique(tid[in_win])
    keep = np.isin(tid, live_ids)
    tid = tid[keep]; time = time[keep]; lon = lon[keep]; lat = lat[keep]

    # Sort by (track, time) so per-track slices are already time-ordered.
    order = np.lexsort((time, tid))
    tid = tid[order]; time = time[order]; lon = lon[order]; lat = lat[order]

    edges = np.r_[0, np.where(np.diff(tid))[0] + 1, len(tid)]
    groups = []
    for a, b in zip(edges[:-1], edges[1:]):
        groups.append((time[a:b], lon[a:b], lat[a:b]))
    return groups


def open_adt(adt_dir: Path, pattern: str, date: np.datetime64):
    ymd = np.datetime_as_string(date, unit="D").replace("-", "")
    path = adt_dir / pattern.format(ymd=ymd)
    ds = xr.open_dataset(path)
    field = ds["adt"].isel(time=0).values
    lat = ds["latitude"].values
    lon = ds["longitude"].values
    # Shift to [-180, 180] for cartopy.
    lon_pm = ((lon + 180) % 360) - 180
    order = np.argsort(lon_pm)
    return field[:, order], lat, lon_pm[order]


def _split_at_dateline(lon, lat, jump=180.0):
    """Yield (lon, lat) sub-segments where consecutive lon jumps > jump are cuts.
    Prevents long horizontal streaks across the map when a track crosses ±180°."""
    if len(lon) == 0:
        return
    cuts = np.where(np.abs(np.diff(lon)) > jump)[0] + 1
    start = 0
    for c in cuts:
        yield lon[start:c], lat[start:c]
        start = c
    yield lon[start:], lat[start:]


def draw_trails(ax, groups, t_now, trail_days, color):
    n_live = 0
    t_cut = t_now - np.timedelta64(trail_days, "D")
    for t, lo, la in groups:
        i0 = np.searchsorted(t, t_cut, side="left")
        i1 = np.searchsorted(t, t_now, side="right")
        if i1 <= i0:
            continue
        # Only draw if the most recent obs is exactly t_now (eddy alive now).
        if t[i1 - 1] != t_now:
            continue
        seg_lon = lo[i0:i1]
        seg_lat = la[i0:i1]
        for sl, sa in _split_at_dateline(seg_lon, seg_lat):
            if len(sl) > 1:
                ax.plot(sl, sa, "-", transform=ccrs.PlateCarree(),
                        color=color, linewidth=0.45, alpha=0.7)
        ax.plot(seg_lon[-1], seg_lat[-1], ".", transform=ccrs.PlateCarree(),
                color=color, markersize=1.8)
        n_live += 1
    return n_live


def render_frame(date, adt_dir, adt_pattern, anti_groups, cyc_groups,
                 trail_days, label, out_path, extent=None, vmax=1.2):
    field, lat, lon = open_adt(adt_dir, adt_pattern, date)

    if extent is None:
        fig = plt.figure(figsize=(16, 8))
        ax = plt.axes(projection=ccrs.Robinson())
        ax.set_global()
    else:
        fig = plt.figure(figsize=(13, 9))
        ax = plt.axes(projection=ccrs.PlateCarree())
        ax.set_extent(extent, crs=ccrs.PlateCarree())
        ax.gridlines(draw_labels=True, linewidth=0.3, alpha=0.5)

    pcm = ax.pcolormesh(lon, lat, field, transform=ccrs.PlateCarree(),
                        cmap="RdBu_r", vmin=-vmax, vmax=vmax, shading="auto",
                        rasterized=True)
    ax.add_feature(cfeature.LAND, facecolor="lightgray", zorder=2)
    ax.coastlines(linewidth=0.4, zorder=3)

    n_a = draw_trails(ax, anti_groups, date, trail_days, "firebrick")
    n_c = draw_trails(ax, cyc_groups, date, trail_days, "royalblue")

    cb = plt.colorbar(pcm, ax=ax, orientation="vertical", pad=0.02,
                      shrink=0.75, label="SSH (m)")
    cb.ax.tick_params(labelsize=8)
    date_str = np.datetime_as_string(date, unit="D")
    ax.set_title(f"{label} — {date_str}  (trail ≤ {trail_days} d)",
                 fontsize=11)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=120, bbox_inches="tight")
    plt.close(fig)
    return n_a, n_c


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--source", choices=list(SOURCES), default="ifs-fesom")
    ap.add_argument("--start", required=True, help="YYYY-MM-DD inclusive")
    ap.add_argument("--end", required=True, help="YYYY-MM-DD inclusive")
    ap.add_argument("--trail-days", type=int, default=60)
    ap.add_argument("--region", choices=list(REGIONS), default="global",
                    help="global (Robinson) or a regional box (PlateCarree)")
    ap.add_argument("--out-dir", default=None,
                    help="defaults to /bigdisk/work/eddy/figs/frames_<source>_<region>_<yyyy>")
    ap.add_argument("--vmax", type=float, default=None,
                    help="override per-region default (symmetric ±vmax)")
    args = ap.parse_args()

    cfg = SOURCES[args.source]
    t_start = np.datetime64(args.start)
    t_end = np.datetime64(args.end)
    dates = np.arange(t_start, t_end + np.timedelta64(1, "D"),
                      np.timedelta64(1, "D"))

    # Tracks window needs to extend back by trail_days so trails are complete.
    t_back = t_start - np.timedelta64(args.trail_days, "D")
    print(f"Loading tracks from {cfg['tracks_root']} for "
          f"{np.datetime_as_string(t_back, unit='D')}..{args.end}")
    anti = load_track_groups(cfg["tracks_root"] / "anti/Anticyclonic.nc",
                             t_back, t_end)
    cyc = load_track_groups(cfg["tracks_root"] / "cyc/Cyclonic.nc",
                            t_back, t_end)
    print(f"  kept {len(anti)} anti, {len(cyc)} cyc tracks with obs in window")

    extent = REGIONS[args.region]
    vmax = args.vmax if args.vmax is not None else REGION_VMAX[args.region]
    year = args.start[:4]
    out_dir = Path(args.out_dir) if args.out_dir else (
        OUT_ROOT / f"figs/frames_{args.source}_{args.region}_{year}")
    print(f"Writing frames to {out_dir}")

    for d in dates:
        ymd = np.datetime_as_string(d, unit="D").replace("-", "")
        out = out_dir / f"frame_{ymd}.png"
        if out.exists():
            print(f"  skip {ymd} (exists)")
            continue
        n_a, n_c = render_frame(d, cfg["adt_dir"], cfg["adt_pattern"],
                                anti, cyc, args.trail_days, cfg["label"],
                                out, extent=extent, vmax=vmax)
        print(f"  {ymd}: anti={n_a}, cyc={n_c}  -> {out.name}")


if __name__ == "__main__":
    main()
