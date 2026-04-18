"""Summary plots from py-eddy-tracker tracking output.

Works on either a single year or a multi-year window — the time range is
auto-detected from the tracking NetCDFs.

  # single year
  plot_year_summary.py --label "DestinE 2014" \\
      --track-dir $EDDY_ROOT/out/tracks_ifs-fesom_hist_2014 \\
      --fig-dir   $EDDY_ROOT/figs/year_2014_ifs-fesom_hist

  # multi-year (kicked off by pipeline_multiyear.py)
  plot_year_summary.py --label "DestinE 1990-2014" \\
      --track-dir $EDDY_ROOT/out/tracks_ifs-fesom_hist_1990_2014 \\
      --fig-dir   $EDDY_ROOT/figs/window_1990_2014_ifs-fesom_hist

Produces:
  tracks_global.png             — all tracks ≥ MIN_LIFETIME_GLOBAL on a Robinson map
  tracks_<region>.png           — same, zoomed to Gulf Stream, Agulhas, Kuroshio, Brazil
  lifetime_histogram.png        — track-lifetime PDF for both polarities
  count_per_day.png             — eddies per day on a real date axis (anti vs cyc)
  scoreboard.txt                — counts, means, ranges, window length
"""

import argparse
import sys
from pathlib import Path

import cartopy.crs as ccrs
import cartopy.feature as cfeature
import matplotlib.dates as mdates
import matplotlib.pyplot as plt
import numpy as np
import xarray as xr

MIN_LIFETIME_GLOBAL = 30   # days; only "long-lived" tracks get a line
MIN_LIFETIME_REGION = 14   # more permissive in zooms

REGIONS = {
    "gulfstream": ("Gulf Stream",          (-82, -40,  25, 50)),
    "agulhas":    ("Agulhas",              ( 10,  40, -45, -25)),
    "kuroshio":   ("Kuroshio Extension",   (130, 180,  25, 45)),
    "brazil":     ("Brazil-Malvinas",      (-60, -30, -50, -25)),
}


def parse_args():
    ap = argparse.ArgumentParser()
    # --label is the title decoration; --year is kept for backwards compat with
    # pipeline_year.py and is converted to a label internally if --label absent.
    ap.add_argument("--label", help="Title decoration, e.g. 'DestinE 1990-2014'")
    ap.add_argument("--year", type=int, help="Single-year shortcut for --label")
    ap.add_argument("--track-dir", type=Path, required=True)
    ap.add_argument("--fig-dir", type=Path, required=True)
    args = ap.parse_args()
    if args.label is None:
        if args.year is None:
            sys.exit("Pass either --label or --year")
        args.label = str(args.year)
    return args


def load_tracks(path):
    """Return (track_id, lon[-180..180], lat, time_seconds, effective_radius_m)."""
    ds = xr.open_dataset(path)
    tid = ds["track"].values
    lon = ((ds["longitude"].values + 180) % 360) - 180
    lat = ds["latitude"].values
    t   = ds["time"].values  # seconds since 1950-01-01 by py-eddy-tracker convention
    r   = ds["effective_radius"].values
    return tid, lon, lat, t, r


def group_by_track(tid, lon, lat, t, r, min_obs):
    order = np.argsort(tid, kind="stable")
    tid = tid[order]; lon = lon[order]; lat = lat[order]
    t = t[order]; r = r[order]
    edges = np.r_[0, np.where(np.diff(tid))[0] + 1, len(tid)]
    out = []
    for a, b in zip(edges[:-1], edges[1:]):
        if b - a >= min_obs:
            out.append({
                "lon": lon[a:b], "lat": lat[a:b], "rad": r[a:b],
                "t": t[a:b], "duration": b - a,
            })
    return out


def _break_at_dateline(lon, lat):
    """Insert NaN where lon jumps by >180° between consecutive samples.
    matplotlib breaks the line at NaNs, killing the straight-across-globe
    segment that would otherwise span 360° at the antimeridian."""
    if len(lon) < 2:
        return lon, lat
    jumps = np.where(np.abs(np.diff(lon)) > 180)[0]
    if not len(jumps):
        return lon, lat
    return (np.insert(lon, jumps + 1, np.nan),
            np.insert(lat, jumps + 1, np.nan))


def plot_global_tracks(anti_tracks, cyc_tracks, out_path, label):
    fig = plt.figure(figsize=(15, 8))
    ax = plt.axes(projection=ccrs.Robinson())
    ax.set_global()
    ax.add_feature(cfeature.LAND, facecolor="lightgray", zorder=1)
    ax.coastlines(linewidth=0.4, zorder=2)

    for t in anti_tracks:
        lon, lat = _break_at_dateline(t["lon"], t["lat"])
        ax.plot(lon, lat, "-", transform=ccrs.PlateCarree(),
                color="firebrick", linewidth=0.25, alpha=0.4)
    for t in cyc_tracks:
        lon, lat = _break_at_dateline(t["lon"], t["lat"])
        ax.plot(lon, lat, "-", transform=ccrs.PlateCarree(),
                color="royalblue", linewidth=0.25, alpha=0.4)

    import matplotlib.patches as mpatches
    ax.legend(handles=[
        mpatches.Patch(color="firebrick", label=f"Anticyclonic ({len(anti_tracks)})"),
        mpatches.Patch(color="royalblue", label=f"Cyclonic ({len(cyc_tracks)})"),
    ], loc="lower left", fontsize=10, framealpha=0.9)
    ax.set_title(f"{label} — eddy tracks lifetime ≥ {MIN_LIFETIME_GLOBAL} d")
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"Wrote {out_path}")


def plot_region_tracks(anti_tracks, cyc_tracks, region_key, out_path, label):
    name, ext = REGIONS[region_key]
    l0, l1, la0, la1 = ext
    fig = plt.figure(figsize=(13, 8))
    ax = plt.axes(projection=ccrs.PlateCarree())
    ax.set_extent(ext, crs=ccrs.PlateCarree())
    ax.add_feature(cfeature.LAND, facecolor="lightgray", zorder=1)
    ax.coastlines(linewidth=0.5, zorder=2)
    ax.gridlines(draw_labels=True, linewidth=0.3, alpha=0.5)

    def in_box(lo, la):
        return (lo >= l0) & (lo <= l1) & (la >= la0) & (la <= la1)

    n_a = n_c = 0
    for t in anti_tracks:
        m = in_box(t["lon"], t["lat"])
        if m.any():
            ax.plot(t["lon"][m], t["lat"][m], "-", transform=ccrs.PlateCarree(),
                    color="firebrick", linewidth=0.5, alpha=0.55)
            ax.plot(t["lon"][m][0], t["lat"][m][0], ".",
                    transform=ccrs.PlateCarree(), color="firebrick", markersize=2)
            n_a += 1
    for t in cyc_tracks:
        m = in_box(t["lon"], t["lat"])
        if m.any():
            ax.plot(t["lon"][m], t["lat"][m], "-", transform=ccrs.PlateCarree(),
                    color="royalblue", linewidth=0.5, alpha=0.55)
            ax.plot(t["lon"][m][0], t["lat"][m][0], ".",
                    transform=ccrs.PlateCarree(), color="royalblue", markersize=2)
            n_c += 1
    import matplotlib.patches as mpatches
    ax.legend(handles=[
        mpatches.Patch(color="firebrick", label=f"Anticyclonic ({n_a})"),
        mpatches.Patch(color="royalblue", label=f"Cyclonic ({n_c})"),
    ], loc="lower left", fontsize=9, framealpha=0.9)
    ax.set_title(f"{name} — {label} — tracks ≥ {MIN_LIFETIME_REGION} d")
    fig.savefig(out_path, dpi=140, bbox_inches="tight")
    plt.close(fig)
    print(f"Wrote {out_path}")


def plot_lifetime_histogram(anti_durations, cyc_durations, out_path, label):
    fig, ax = plt.subplots(figsize=(10, 6))
    longest = max(anti_durations.max(), cyc_durations.max())
    # Coarser bins when the window is many years (otherwise we get an unreadably
    # spiky tail). 7-day bins ≤ 1 year, 14-day bins ≤ 5 years, 30-day beyond.
    bin_w = 7 if longest <= 365 else (14 if longest <= 5 * 365 else 30)
    bins = np.arange(0, longest + bin_w, bin_w)
    ax.hist(anti_durations, bins=bins, alpha=0.55, color="firebrick",
            label=f"Anticyclonic (n={len(anti_durations)}, median={int(np.median(anti_durations))} d)")
    ax.hist(cyc_durations, bins=bins, alpha=0.55, color="royalblue",
            label=f"Cyclonic (n={len(cyc_durations)}, median={int(np.median(cyc_durations))} d)")
    ax.set_xlabel(f"Track lifetime (days) — {bin_w}-day bins")
    ax.set_ylabel("Number of tracks")
    ax.set_yscale("log")
    ax.set_title(f"Eddy track lifetimes — {label}")
    ax.legend()
    ax.grid(alpha=0.3)
    fig.savefig(out_path, dpi=130, bbox_inches="tight")
    plt.close(fig)
    print(f"Wrote {out_path}")


def plot_count_per_day(anti_path, cyc_path, out_path, label):
    """One point per day — total observations of each polarity that day, on
    a real date axis. Auto-spans whatever window the tracking output covers."""
    a = xr.open_dataset(anti_path)
    c = xr.open_dataset(cyc_path)
    t_a = np.array(a["time"].values, dtype="datetime64[D]")
    t_c = np.array(c["time"].values, dtype="datetime64[D]")
    t_min = min(t_a.min(), t_c.min())
    t_max = max(t_a.max(), t_c.max())
    n_days = int((t_max - t_min).astype(int)) + 1
    days_a = (t_a - t_min).astype(int)
    days_c = (t_c - t_min).astype(int)
    counts_a = np.bincount(days_a, minlength=n_days)[:n_days]
    counts_c = np.bincount(days_c, minlength=n_days)[:n_days]
    dates = t_min + np.arange(n_days).astype("timedelta64[D]")

    fig, ax = plt.subplots(figsize=(14, 5))
    ax.plot(dates, counts_a, color="firebrick", linewidth=0.6, label="Anticyclonic")
    ax.plot(dates, counts_c, color="royalblue", linewidth=0.6, label="Cyclonic")
    ax.set_ylabel("Detected eddies per day")
    ax.set_title(f"Daily eddy count — {label}")
    # Sensible x-axis tick spacing depending on window length
    n_years = (t_max.astype("datetime64[Y]").astype(int) -
               t_min.astype("datetime64[Y]").astype(int)) + 1
    if n_years <= 2:
        ax.xaxis.set_major_locator(mdates.MonthLocator(interval=1 if n_years == 1 else 2))
        ax.xaxis.set_major_formatter(mdates.DateFormatter("%Y-%m"))
    elif n_years <= 10:
        ax.xaxis.set_major_locator(mdates.YearLocator())
        ax.xaxis.set_major_formatter(mdates.DateFormatter("%Y"))
    else:
        ax.xaxis.set_major_locator(mdates.YearLocator(base=5))
        ax.xaxis.set_major_formatter(mdates.DateFormatter("%Y"))
    ax.legend()
    ax.grid(alpha=0.3)
    fig.autofmt_xdate()
    fig.savefig(out_path, dpi=130, bbox_inches="tight")
    plt.close(fig)
    print(f"Wrote {out_path}")


def main():
    args = parse_args()
    args.fig_dir.mkdir(parents=True, exist_ok=True)
    anti_path = args.track_dir / "anti" / "Anticyclonic.nc"
    cyc_path  = args.track_dir / "cyc"  / "Cyclonic.nc"
    if not (anti_path.exists() and cyc_path.exists()):
        raise SystemExit(f"Missing {anti_path} or {cyc_path}")

    a_tid, a_lon, a_lat, a_t, a_r = load_tracks(anti_path)
    c_tid, c_lon, c_lat, c_t, c_r = load_tracks(cyc_path)
    print(f"Loaded: anti={len(a_tid)} obs / cyc={len(c_tid)} obs")

    a_long = group_by_track(a_tid, a_lon, a_lat, a_t, a_r, MIN_LIFETIME_GLOBAL)
    c_long = group_by_track(c_tid, c_lon, c_lat, c_t, c_r, MIN_LIFETIME_GLOBAL)
    a_reg  = group_by_track(a_tid, a_lon, a_lat, a_t, a_r, MIN_LIFETIME_REGION)
    c_reg  = group_by_track(c_tid, c_lon, c_lat, c_t, c_r, MIN_LIFETIME_REGION)

    # Global map
    plot_global_tracks(a_long, c_long, args.fig_dir / "tracks_global.png", args.label)

    # Region zooms
    for key in REGIONS:
        plot_region_tracks(a_reg, c_reg, key,
                           args.fig_dir / f"tracks_{key}.png", args.label)

    # Lifetime histogram (use ALL tracks ≥ 4 d, the EddyTracking minimum,
    # so the tail isn't truncated at our display cutoff).
    a_all = group_by_track(a_tid, a_lon, a_lat, a_t, a_r, 4)
    c_all = group_by_track(c_tid, c_lon, c_lat, c_t, c_r, 4)
    plot_lifetime_histogram(
        np.array([t["duration"] for t in a_all]),
        np.array([t["duration"] for t in c_all]),
        args.fig_dir / "lifetime_histogram.png", args.label)

    # Daily count
    plot_count_per_day(anti_path, cyc_path,
                        args.fig_dir / "count_per_day.png", args.label)

    # Window length from the actual data, for the scoreboard.
    a_t_d = np.array(xr.open_dataset(anti_path)["time"].values, dtype="datetime64[D]")
    c_t_d = np.array(xr.open_dataset(cyc_path)["time"].values, dtype="datetime64[D]")
    t_min = min(a_t_d.min(), c_t_d.min())
    t_max = max(a_t_d.max(), c_t_d.max())
    n_days = int((t_max - t_min).astype(int)) + 1

    sb = args.fig_dir / "scoreboard.txt"
    with sb.open("w") as f:
        f.write(f"{args.label}\n")
        f.write(f"Window: {t_min} .. {t_max}  ({n_days} days)\n")
        f.write(f"Anticyclonic: {len(a_tid)} obs, {len(a_all)} tracks ≥4 d, "
                f"{len(a_long)} tracks ≥{MIN_LIFETIME_GLOBAL} d, "
                f"longest = {max(t['duration'] for t in a_all)} d, "
                f"median lifetime = {int(np.median([t['duration'] for t in a_all]))} d\n")
        f.write(f"Cyclonic:     {len(c_tid)} obs, {len(c_all)} tracks ≥4 d, "
                f"{len(c_long)} tracks ≥{MIN_LIFETIME_GLOBAL} d, "
                f"longest = {max(t['duration'] for t in c_all)} d, "
                f"median lifetime = {int(np.median([t['duration'] for t in c_all]))} d\n")
    print(f"Wrote {sb}")
    print(sb.read_text())


if __name__ == "__main__":
    main()
