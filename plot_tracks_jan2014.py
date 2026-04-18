"""Plot the longest-lived eddy trajectories from Jan 2014 tracking output.

Two figures:
  - Global map of all tracks lasting ≥ 14 days, anti vs cyc colour-coded.
  - Gulf Stream zoom of the same.

Run in the eddy env:
    $EDDY_PYTHON plot_tracks_jan2014.py
"""

from pathlib import Path

import cartopy.crs as ccrs
import cartopy.feature as cfeature
import matplotlib.pyplot as plt
import numpy as np
import xarray as xr

import eddy_paths
ROOT = eddy_paths.out_dir() / "tracks_ifs-fesom_hist"
ANTI = ROOT / "anti" / "Anticyclonic.nc"
CYC = ROOT / "cyc" / "Cyclonic.nc"
FIG_DIR = eddy_paths.fig_dir()
FIG_DIR.mkdir(parents=True, exist_ok=True)

MIN_DURATION = 14   # observations (= days) for "long-lived"


def load_tracks(path):
    ds = xr.open_dataset(path)
    track_id = ds["track"].values
    lon = ((ds["longitude"].values + 180) % 360) - 180
    lat = ds["latitude"].values
    rad = ds["effective_radius"].values
    return track_id, lon, lat, rad


def by_track(track_id, lon, lat, rad, min_obs):
    """Group obs by track id; return list of (lon_traj, lat_traj, mean_rad_km)."""
    order = np.argsort(track_id, kind="stable")
    tid = track_id[order]; lon = lon[order]; lat = lat[order]; rad = rad[order]
    out = []
    edges = np.r_[0, np.where(np.diff(tid))[0] + 1, len(tid)]
    for a, b in zip(edges[:-1], edges[1:]):
        if b - a >= min_obs:
            out.append((lon[a:b], lat[a:b], float(rad[a:b].mean()) / 1e3))
    return out


def plot_map(extent, anti, cyc, out_path, title):
    """extent = (lon0, lon1, lat0, lat1)."""
    fig = plt.figure(figsize=(14, 8))
    proj = ccrs.PlateCarree() if extent != "global" else ccrs.Robinson()
    ax = plt.axes(projection=proj)
    if extent == "global":
        ax.set_global()
    else:
        ax.set_extent(extent, crs=ccrs.PlateCarree())
    ax.add_feature(cfeature.LAND, facecolor="lightgray", zorder=1)
    ax.coastlines(linewidth=0.4, zorder=2)
    if extent != "global":
        ax.gridlines(draw_labels=True, linewidth=0.3, alpha=0.5)

    def in_box(lon, lat):
        if extent == "global":
            return np.ones_like(lon, dtype=bool)
        l0, l1, la0, la1 = extent
        return (lon >= l0) & (lon <= l1) & (lat >= la0) & (lat <= la1)

    n_anti_drawn = 0
    for lon_t, lat_t, _ in anti:
        m = in_box(lon_t, lat_t)
        if m.any():
            ax.plot(lon_t[m], lat_t[m], "-", transform=ccrs.PlateCarree(),
                    color="firebrick", linewidth=0.5, alpha=0.55)
            ax.plot(lon_t[m][0], lat_t[m][0], ".", transform=ccrs.PlateCarree(),
                    color="firebrick", markersize=2)
            n_anti_drawn += 1
    n_cyc_drawn = 0
    for lon_t, lat_t, _ in cyc:
        m = in_box(lon_t, lat_t)
        if m.any():
            ax.plot(lon_t[m], lat_t[m], "-", transform=ccrs.PlateCarree(),
                    color="royalblue", linewidth=0.5, alpha=0.55)
            ax.plot(lon_t[m][0], lat_t[m][0], ".", transform=ccrs.PlateCarree(),
                    color="royalblue", markersize=2)
            n_cyc_drawn += 1

    import matplotlib.patches as mpatches
    ax.legend(handles=[
        mpatches.Patch(color="firebrick", label=f"Anticyclonic ({n_anti_drawn})"),
        mpatches.Patch(color="royalblue", label=f"Cyclonic ({n_cyc_drawn})"),
    ], loc="lower left", fontsize=9, framealpha=0.9)
    ax.set_title(title)
    fig.savefig(out_path, dpi=130, bbox_inches="tight")
    plt.close(fig)
    print(f"Wrote {out_path}  (anti={n_anti_drawn}, cyc={n_cyc_drawn})")


def main():
    a = by_track(*load_tracks(ANTI), min_obs=MIN_DURATION)
    c = by_track(*load_tracks(CYC), min_obs=MIN_DURATION)
    print(f"Tracks ≥ {MIN_DURATION} days: anticyclonic={len(a)}, cyclonic={len(c)}")
    print(f"Track length stats: anti mean={np.mean([len(t[0]) for t in a]):.1f}d, "
          f"max={max(len(t[0]) for t in a)}d ; "
          f"cyc mean={np.mean([len(t[0]) for t in c]):.1f}d, "
          f"max={max(len(t[0]) for t in c)}d")

    plot_map("global", a, c, FIG_DIR / "tracks_global_jan2014.png",
             f"DestinE IFS-FESOM hist Jan 2014 — eddy tracks lasting ≥ {MIN_DURATION} d")
    plot_map((-82, -40, 25, 50), a, c,
             FIG_DIR / "tracks_gulfstream_jan2014.png",
             f"Gulf Stream eddy tracks ≥ {MIN_DURATION} d (Jan 2014)")
    plot_map((10, 40, -45, -25), a, c,
             FIG_DIR / "tracks_agulhas_jan2014.png",
             f"Agulhas region eddy tracks ≥ {MIN_DURATION} d (Jan 2014)")


if __name__ == "__main__":
    main()
