"""Per-model statistical summary + ICON vs IFS-FESOM comparison.

Produces two sets of figures in --fig-dir:
  a) per-model diagnostics      (one panel per model per figure, or stacked 2x2)
  b) side-by-side comparisons   (A vs B histograms, difference maps)

Inputs are the usual py-eddy-tracker tracking outputs — a directory with
  anti/Anticyclonic.nc  and  cyc/Cyclonic.nc

Run in the eddy env:
  $EDDY_PYTHON plot_model_comparison.py \\
    --label-a ICON       --track-dir-a $EDDY_ROOT/out/tracks_icon_hist_1990_2014 \\
    --label-b IFS-FESOM  --track-dir-b $EDDY_ROOT/out/tracks_ifs-fesom_hist_1990_2014 \\
    --fig-dir $EDDY_ROOT/figs/compare_icon_vs_fesom_1990_2014
"""

from __future__ import annotations

import argparse
from pathlib import Path

import cartopy.crs as ccrs
import cartopy.feature as cfeature
import matplotlib.colors as mcolors
import matplotlib.dates as mdates
import matplotlib.pyplot as plt
import numpy as np
import xarray as xr

# Spatial grid for density / mean maps.
BIN_DEG = 1.0
LON_EDGES = np.arange(-180.0, 180.0 + BIN_DEG, BIN_DEG)
LAT_EDGES = np.arange(-90.0, 90.0 + BIN_DEG, BIN_DEG)
LON_CTR = 0.5 * (LON_EDGES[:-1] + LON_EDGES[1:])
LAT_CTR = 0.5 * (LAT_EDGES[:-1] + LAT_EDGES[1:])

MIN_COUNT_FOR_MEAN = 5  # hide bins with < N obs from mean-value maps

POLARITIES = [
    ("anti", "Anticyclonic", "firebrick"),
    ("cyc",  "Cyclonic",     "royalblue"),
]


# ---------------------------------------------------------------------------
# loaders
# ---------------------------------------------------------------------------

def load_polarity(nc_path: Path):
    """Return a dict of compact arrays for one polarity of one model.

    Float32 for spatial / scalar fields to keep memory modest; uint32 for
    track ids; datetime64[D] for time (day-precision is enough for
    time-series and season-of-year grouping)."""
    ds = xr.open_dataset(nc_path)
    tid = ds["track"].values.astype(np.uint32, copy=False)
    lon = ((ds["longitude"].values + 180.0) % 360.0) - 180.0
    lat = ds["latitude"].values
    t   = ds["time"].values.astype("datetime64[D]")
    rad = ds["effective_radius"].values  # m
    amp = ds["amplitude"].values         # m
    return {
        "tid": tid,
        "lon": lon.astype(np.float32, copy=False),
        "lat": lat.astype(np.float32, copy=False),
        "t":   t,
        "rad_km": (rad.astype(np.float32, copy=False) / 1e3),
        "amp_cm": (amp.astype(np.float32, copy=False) * 1e2),
    }


def load_model(track_dir: Path):
    anti_p = track_dir / "anti" / "Anticyclonic.nc"
    cyc_p  = track_dir / "cyc"  / "Cyclonic.nc"
    if not (anti_p.exists() and cyc_p.exists()):
        raise SystemExit(f"Missing tracks under {track_dir}")
    return {
        "anti": load_polarity(anti_p),
        "cyc":  load_polarity(cyc_p),
    }


# ---------------------------------------------------------------------------
# derived statistics
# ---------------------------------------------------------------------------

def track_lifetimes(tid: np.ndarray):
    """Per-track lifetime in obs (≈ days since detection cadence is daily)."""
    # tid is not globally sorted, but np.unique with return_counts gives duration
    # directly regardless of order — one obs per (track, day).
    _, counts = np.unique(tid, return_counts=True)
    return counts.astype(np.int32)


def count_grid(lon, lat):
    """2D observation count on the BIN_DEG grid."""
    H, _, _ = np.histogram2d(lat, lon, bins=[LAT_EDGES, LON_EDGES])
    return H.astype(np.float32)


def mean_grid(lon, lat, value, min_count=MIN_COUNT_FOR_MEAN):
    """Mean of `value` per (lat, lon) bin. NaN where count < min_count."""
    sum_grid, _, _ = np.histogram2d(lat, lon, bins=[LAT_EDGES, LON_EDGES],
                                    weights=value)
    cnt_grid, _, _ = np.histogram2d(lat, lon, bins=[LAT_EDGES, LON_EDGES])
    with np.errstate(invalid="ignore", divide="ignore"):
        mean = sum_grid / cnt_grid
    mean[cnt_grid < min_count] = np.nan
    return mean.astype(np.float32)


def daily_counts(t: np.ndarray, t_min: np.datetime64, t_max: np.datetime64):
    """Observations per calendar day across [t_min, t_max]."""
    n = int((t_max - t_min).astype(int)) + 1
    d = (t - t_min).astype(int)
    d = d[(d >= 0) & (d < n)]
    return np.bincount(d, minlength=n)[:n]


# ---------------------------------------------------------------------------
# plotting primitives
# ---------------------------------------------------------------------------

def _global_ax(fig, n_rows=1, n_cols=1, index=1):
    ax = fig.add_subplot(n_rows, n_cols, index, projection=ccrs.Robinson())
    ax.set_global()
    ax.add_feature(cfeature.LAND, facecolor="lightgray", zorder=1)
    ax.coastlines(linewidth=0.3, zorder=2)
    return ax


def _pcm_global(ax, grid_ll, cmap, norm=None, vmin=None, vmax=None):
    """pcolormesh on the LAT_CTR/LON_CTR grid (grid_ll shape [nlat, nlon])."""
    if norm is None:
        return ax.pcolormesh(LON_EDGES, LAT_EDGES, grid_ll,
                             transform=ccrs.PlateCarree(),
                             cmap=cmap, vmin=vmin, vmax=vmax,
                             shading="auto", rasterized=True)
    return ax.pcolormesh(LON_EDGES, LAT_EDGES, grid_ll,
                         transform=ccrs.PlateCarree(),
                         cmap=cmap, norm=norm, shading="auto", rasterized=True)


# ---------------------------------------------------------------------------
# per-model figures
# ---------------------------------------------------------------------------

def plot_density_map_one_model(model_name, data, fig_path):
    """2-panel (anti, cyc) log-density map."""
    fig = plt.figure(figsize=(16, 9))
    for i, (key, label, _) in enumerate(POLARITIES, 1):
        d = data[key]
        grid = count_grid(d["lon"], d["lat"])
        ax = _global_ax(fig, 2, 1, i)
        pcm = _pcm_global(ax, grid, cmap="viridis",
                          norm=mcolors.LogNorm(vmin=1, vmax=max(10, grid.max())))
        cb = plt.colorbar(pcm, ax=ax, orientation="vertical", pad=0.02,
                          shrink=0.75, label="obs / 1° bin (log)")
        cb.ax.tick_params(labelsize=8)
        ax.set_title(f"{model_name} — {label} eddy observation density (1990-2014)")
    fig.savefig(fig_path, dpi=130, bbox_inches="tight")
    plt.close(fig)
    print(f"  wrote {fig_path.name}")


def plot_mean_map_one_model(model_name, data, field, unit, fig_path,
                             vmin=None, vmax=None, cmap="plasma"):
    """2-panel mean of `field` (e.g. 'rad_km' or 'amp_cm') on the spatial grid."""
    fig = plt.figure(figsize=(16, 9))
    for i, (key, label, _) in enumerate(POLARITIES, 1):
        d = data[key]
        grid = mean_grid(d["lon"], d["lat"], d[field])
        ax = _global_ax(fig, 2, 1, i)
        pcm = _pcm_global(ax, grid, cmap=cmap, vmin=vmin, vmax=vmax)
        cb = plt.colorbar(pcm, ax=ax, orientation="vertical", pad=0.02,
                          shrink=0.75, label=f"mean {field.split('_')[0]} ({unit})")
        cb.ax.tick_params(labelsize=8)
        ax.set_title(f"{model_name} — {label} mean {field.split('_')[0]} "
                     f"(bins with ≥{MIN_COUNT_FOR_MEAN} obs)")
    fig.savefig(fig_path, dpi=130, bbox_inches="tight")
    plt.close(fig)
    print(f"  wrote {fig_path.name}")


# ---------------------------------------------------------------------------
# comparison figures
# ---------------------------------------------------------------------------

def plot_density_diff(name_a, a, name_b, b, fig_path):
    """Two panels (anti, cyc): A − B per 1° bin. Diverging scale."""
    fig = plt.figure(figsize=(16, 9))
    for i, (key, label, _) in enumerate(POLARITIES, 1):
        ga = count_grid(a[key]["lon"], a[key]["lat"])
        gb = count_grid(b[key]["lon"], b[key]["lat"])
        diff = ga - gb
        vmax = float(np.nanpercentile(np.abs(diff), 99))
        ax = _global_ax(fig, 2, 1, i)
        pcm = _pcm_global(ax, diff, cmap="RdBu_r", vmin=-vmax, vmax=vmax)
        cb = plt.colorbar(pcm, ax=ax, orientation="vertical", pad=0.02,
                          shrink=0.75, label=f"{name_a} − {name_b}  (obs / 1° bin)")
        cb.ax.tick_params(labelsize=8)
        ax.set_title(f"{label} eddy density difference  ({name_a} − {name_b})")
    fig.savefig(fig_path, dpi=130, bbox_inches="tight")
    plt.close(fig)
    print(f"  wrote {fig_path.name}")


def plot_mean_field_diff(name_a, a, name_b, b, field, unit, fig_path):
    """A mean(field) − B mean(field)."""
    fig = plt.figure(figsize=(16, 9))
    for i, (key, label, _) in enumerate(POLARITIES, 1):
        ma = mean_grid(a[key]["lon"], a[key]["lat"], a[key][field])
        mb = mean_grid(b[key]["lon"], b[key]["lat"], b[key][field])
        diff = ma - mb
        vmax = float(np.nanpercentile(np.abs(diff), 99))
        ax = _global_ax(fig, 2, 1, i)
        pcm = _pcm_global(ax, diff, cmap="RdBu_r", vmin=-vmax, vmax=vmax)
        cb = plt.colorbar(pcm, ax=ax, orientation="vertical", pad=0.02,
                          shrink=0.75, label=f"Δ mean {field.split('_')[0]} ({unit})")
        cb.ax.tick_params(labelsize=8)
        ax.set_title(f"{label} mean-{field.split('_')[0]} difference  "
                     f"({name_a} − {name_b})")
    fig.savefig(fig_path, dpi=130, bbox_inches="tight")
    plt.close(fig)
    print(f"  wrote {fig_path.name}")


def plot_hist_compare(name_a, a, name_b, b, field, bins, xlabel, fig_path,
                      log_y=True):
    """Step-line histograms overlaying both models, for each polarity."""
    fig, axes = plt.subplots(1, 2, figsize=(14, 5), sharey=True)
    for ax, (key, label, _color) in zip(axes, POLARITIES):
        for name, data, ls in [(name_a, a, "-"), (name_b, b, "--")]:
            if field == "lifetime":
                vals = track_lifetimes(data[key]["tid"])
            elif field in data[key]:
                vals = data[key][field]
            else:
                continue
            ax.hist(vals, bins=bins, histtype="step", linewidth=1.6,
                    linestyle=ls,
                    label=f"{name} (n={len(vals):,}, "
                          f"median={float(np.median(vals)):.1f})")
        ax.set_xlabel(xlabel)
        ax.set_title(label)
        if log_y:
            ax.set_yscale("log")
        ax.grid(alpha=0.3)
        ax.legend(fontsize=8)
    axes[0].set_ylabel("count")
    fig.suptitle(f"{xlabel} distribution — {name_a} vs {name_b}")
    fig.savefig(fig_path, dpi=130, bbox_inches="tight")
    plt.close(fig)
    print(f"  wrote {fig_path.name}")


def plot_lat_histogram(name_a, a, name_b, b, fig_path):
    """1° latitude-binned observation counts, one panel per polarity."""
    fig, axes = plt.subplots(1, 2, figsize=(14, 5), sharey=True)
    for ax, (key, label, _) in zip(axes, POLARITIES):
        for name, data, ls in [(name_a, a, "-"), (name_b, b, "--")]:
            h, _ = np.histogram(data[key]["lat"], bins=LAT_EDGES)
            ax.plot(LAT_CTR, h, linestyle=ls, linewidth=1.5,
                    label=f"{name} (n={len(data[key]['lat']):,})")
        ax.set_xlabel("latitude (deg)")
        ax.set_title(label)
        ax.grid(alpha=0.3)
        ax.legend(fontsize=8)
    axes[0].set_ylabel("observations per 1° lat")
    fig.suptitle(f"Latitudinal distribution — {name_a} vs {name_b}")
    fig.savefig(fig_path, dpi=130, bbox_inches="tight")
    plt.close(fig)
    print(f"  wrote {fig_path.name}")


def plot_count_per_day_compare(name_a, a, name_b, b, fig_path):
    """Daily observation counts, overlay anti+cyc for both models."""
    t_min = min(a["anti"]["t"].min(), a["cyc"]["t"].min(),
                b["anti"]["t"].min(), b["cyc"]["t"].min())
    t_max = max(a["anti"]["t"].max(), a["cyc"]["t"].max(),
                b["anti"]["t"].max(), b["cyc"]["t"].max())
    dates = t_min + np.arange(int((t_max - t_min).astype(int)) + 1) \
                    .astype("timedelta64[D]")
    fig, axes = plt.subplots(2, 1, figsize=(14, 7), sharex=True)
    for ax, (key, label, _) in zip(axes, POLARITIES):
        for name, data, ls in [(name_a, a, "-"), (name_b, b, "--")]:
            c = daily_counts(data[key]["t"], t_min, t_max)
            ax.plot(dates, c, linestyle=ls, linewidth=0.6,
                    label=f"{name} mean={c.mean():.0f}/day")
        ax.set_ylabel("obs / day")
        ax.set_title(f"{label}")
        ax.grid(alpha=0.3)
        ax.legend(fontsize=9, loc="upper right")
    axes[-1].xaxis.set_major_locator(mdates.YearLocator(base=2))
    axes[-1].xaxis.set_major_formatter(mdates.DateFormatter("%Y"))
    axes[-1].set_xlabel("date")
    fig.suptitle(f"Daily eddy observation count — {name_a} vs {name_b}")
    fig.autofmt_xdate()
    fig.savefig(fig_path, dpi=130, bbox_inches="tight")
    plt.close(fig)
    print(f"  wrote {fig_path.name}")


def plot_season_cycle(name_a, a, name_b, b, fig_path):
    """Mean observations per month-of-year — shows seasonal cycle."""
    months = np.arange(1, 13)
    fig, axes = plt.subplots(1, 2, figsize=(14, 5), sharey=True)
    for ax, (key, label, _) in zip(axes, POLARITIES):
        for name, data, ls in [(name_a, a, "-"), (name_b, b, "--")]:
            mo = data[key]["t"].astype("datetime64[M]").astype(int) % 12 + 1
            h, _ = np.histogram(mo, bins=np.arange(0.5, 13.5, 1))
            years = (data[key]["t"].astype("datetime64[Y]").astype(int).max()
                     - data[key]["t"].astype("datetime64[Y]").astype(int).min() + 1)
            ax.plot(months, h / years, marker="o", linestyle=ls,
                    label=f"{name} ({years} yr)")
        ax.set_xlabel("month of year")
        ax.set_xticks(months)
        ax.set_title(label)
        ax.grid(alpha=0.3)
        ax.legend(fontsize=9)
    axes[0].set_ylabel("mean obs per month")
    fig.suptitle(f"Seasonal cycle — {name_a} vs {name_b}")
    fig.savefig(fig_path, dpi=130, bbox_inches="tight")
    plt.close(fig)
    print(f"  wrote {fig_path.name}")


# ---------------------------------------------------------------------------
# scoreboard
# ---------------------------------------------------------------------------

def write_scoreboard(name_a, a, name_b, b, path):
    def row(name, data):
        lines = [f"=== {name} ==="]
        for key, label, _ in POLARITIES:
            d = data[key]
            life = track_lifetimes(d["tid"])
            lines.append(
                f"{label:13s}  obs={len(d['tid']):>10,}  "
                f"tracks={len(life):>8,}  "
                f"median life={np.median(life):5.1f} d  "
                f"≥30 d={int((life >= 30).sum()):>7,}  "
                f"longest={life.max():>5} d  "
                f"median rad={np.median(d['rad_km']):5.1f} km  "
                f"median amp={np.median(d['amp_cm']):5.1f} cm"
            )
        return "\n".join(lines)
    with path.open("w") as f:
        f.write(row(name_a, a) + "\n\n")
        f.write(row(name_b, b) + "\n")
    print(f"  wrote {path.name}")
    print(path.read_text())


# ---------------------------------------------------------------------------
# main
# ---------------------------------------------------------------------------

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--label-a", required=True)
    ap.add_argument("--track-dir-a", type=Path, required=True)
    ap.add_argument("--label-b", required=True)
    ap.add_argument("--track-dir-b", type=Path, required=True)
    ap.add_argument("--fig-dir", type=Path, required=True)
    args = ap.parse_args()
    args.fig_dir.mkdir(parents=True, exist_ok=True)

    print(f"Loading {args.label_a} from {args.track_dir_a} ...")
    a = load_model(args.track_dir_a)
    print(f"  anti obs={len(a['anti']['tid']):,}  cyc obs={len(a['cyc']['tid']):,}")
    print(f"Loading {args.label_b} from {args.track_dir_b} ...")
    b = load_model(args.track_dir_b)
    print(f"  anti obs={len(b['anti']['tid']):,}  cyc obs={len(b['cyc']['tid']):,}")

    fig_dir = args.fig_dir
    la, lb = args.label_a, args.label_b
    slug_a = la.lower().replace(" ", "").replace("-", "")
    slug_b = lb.lower().replace(" ", "").replace("-", "")

    # --- per-model maps ---
    print("Per-model spatial maps:")
    for name, data, slug in [(la, a, slug_a), (lb, b, slug_b)]:
        plot_density_map_one_model(name, data,
            fig_dir / f"density_{slug}.png")
        plot_mean_map_one_model(name, data, "amp_cm", "cm",
            fig_dir / f"mean_amplitude_{slug}.png", vmin=0, vmax=15)
        plot_mean_map_one_model(name, data, "rad_km", "km",
            fig_dir / f"mean_radius_{slug}.png", vmin=30, vmax=130)

    # --- comparison: difference maps ---
    print("Comparison spatial maps:")
    plot_density_diff(la, a, lb, b, fig_dir / "density_diff.png")
    plot_mean_field_diff(la, a, lb, b, "amp_cm", "cm",
                         fig_dir / "mean_amplitude_diff.png")
    plot_mean_field_diff(la, a, lb, b, "rad_km", "km",
                         fig_dir / "mean_radius_diff.png")

    # --- distributions ---
    print("Distribution comparisons:")
    # Lifetime — track-level. Use a synthetic field name handled in plot_hist_compare.
    plot_hist_compare(la, a, lb, b, "lifetime",
                      bins=np.arange(0, 400, 7),
                      xlabel="track lifetime (days)",
                      fig_path=fig_dir / "hist_lifetime.png")
    # Effective radius — obs-level, km.
    plot_hist_compare(la, a, lb, b, "rad_km",
                      bins=np.arange(0, 301, 5),
                      xlabel="effective radius (km)",
                      fig_path=fig_dir / "hist_radius.png")
    # Amplitude — obs-level, cm.
    plot_hist_compare(la, a, lb, b, "amp_cm",
                      bins=np.arange(0, 80, 1),
                      xlabel="amplitude (cm)",
                      fig_path=fig_dir / "hist_amplitude.png")

    # --- 1D spatial / temporal ---
    print("1D comparisons:")
    plot_lat_histogram(la, a, lb, b, fig_dir / "lat_histogram.png")
    plot_count_per_day_compare(la, a, lb, b, fig_dir / "count_per_day.png")
    plot_season_cycle(la, a, lb, b, fig_dir / "season_cycle.png")

    # --- scoreboard ---
    write_scoreboard(la, a, lb, b, fig_dir / "scoreboard.txt")


if __name__ == "__main__":
    main()
