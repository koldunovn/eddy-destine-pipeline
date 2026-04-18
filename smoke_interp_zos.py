"""Smoke test: HEALPix → 0.25° NaN/land-aware interpolation of one daily SSH field.

Steps:
  1. Fetch one day of avg_zos (CLTE o2d, IFS-FESOM hist) from DestinE FDB.
  2. Discover the static ocean land mask from that snapshot.
  3. Build wet-renormalized bilinear weights to the CMEMS-aligned 0.25° grid.
     (cell-centred lon 0.125..359.875, lat -89.875..89.875; ascending lat)
  4. Apply the weights to the snapshot, write NetCDF.
  5. Plot a global map and a Gulf Stream zoom for visual sanity check.

Run from the dedl env:
    $EDDY_DEDL_PYTHON smoke_interp_zos.py
"""

import logging
import sys

import eddy_paths
eddy_paths.inject_fdb_xarray()

import healpy as hp
import numpy as np
import xarray as xr

from fdb_xarray import open_climate_dt

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")
log = logging.getLogger("smoke")

CONFIG = str(eddy_paths.fdb_config())
OUT_DIR = eddy_paths.out_dir()
FIG_DIR = eddy_paths.fig_dir()
OUT_DIR.mkdir(parents=True, exist_ok=True)
FIG_DIR.mkdir(parents=True, exist_ok=True)

NSIDE = 1024            # high-res Climate-DT
RES_DEG = 0.25
DATE = "2014-01-01"     # same anchor MOAAP smoke test used
MODEL = "IFS-FESOM"
EXPERIMENT = "hist"
RESOLUTION = "high"
# DestinE ocean uses a sentinel fill on land cells. Empirically it is +9999
# for FESOM zos, but to be robust we also flag any |SSH|>50 m as non-ocean
# (real SSH is ~±2 m globally; tsunamis don't break 50 m).
SSH_PHYS_MAX = 50.0


# ── 1. Fetch one snapshot ────────────────────────────────────────────
def fetch_zos_one_day():
    log.info("Opening avg_zos %s for %s", MODEL, DATE)
    ds = open_climate_dt(
        portfolio="CLTE",
        levtype="o2d",
        model=MODEL,
        experiment=EXPERIMENT,
        start_date=DATE,
        end_date=DATE,
        resolution=RESOLUTION,
        variables=["avg_zos"],
        config_path=CONFIG,
    )
    da = ds["avg_zos"].sel(scenario=EXPERIMENT).isel(time=0)
    log.info("  dims=%s nside=%s", dict(da.sizes), ds.attrs.get("nside"))
    arr = da.values.astype(np.float32)  # (ncell,)
    nonphys = np.abs(arr) > SSH_PHYS_MAX
    log.info("  raw stats: min=%g max=%g mean=%g  nan=%d  |val|>%g: %d",
             np.nanmin(arr), np.nanmax(arr), np.nanmean(arr),
             int(np.isnan(arr).sum()), SSH_PHYS_MAX, int(nonphys.sum()))
    if nonphys.any():
        sample = arr[nonphys][:3]
        log.info("  non-physical sample values: %s", sample.tolist())
    return arr


# ── 2. Land mask ─────────────────────────────────────────────────────
def discover_land_mask(arr):
    """A HEALPix cell is land iff value is NaN or non-physical (sentinel fill)."""
    land = np.isnan(arr) | (np.abs(arr) > SSH_PHYS_MAX)
    log.info("Mask: %d/%d cells land (%.1f%%)", land.sum(), land.size, 100*land.mean())
    return land


# ── 3. Wet-renormalized bilinear weights ────────────────────────────
def build_target_grid(res_deg=RES_DEG):
    """CMEMS-aligned: cell-centred, ascending lat, 0..360 lon."""
    half = res_deg / 2.0
    lat = np.arange(-90.0 + half, 90.0, res_deg)
    lon = np.arange(0.0 + half, 360.0, res_deg)
    return lat.astype(np.float64), lon.astype(np.float64)


def build_wet_weights(nside, lat_1d, lon_1d, land_mask):
    lon2d, lat2d = np.meshgrid(lon_1d, lat_1d)
    theta = np.deg2rad(90.0 - lat2d.ravel())
    phi = np.deg2rad(lon2d.ravel())

    idx, w = hp.get_interp_weights(nside, theta, phi, nest=True)  # (4, N)
    idx = idx.T.astype(np.int64)         # (N, 4)
    w = w.T.astype(np.float32)           # (N, 4)

    # Zero out weights pointing at land cells; renormalize remaining over wet sum.
    wet_neighbour = ~land_mask[idx]      # (N, 4) True where source is ocean
    w_wet = w * wet_neighbour
    wet_sum = w_wet.sum(axis=1, keepdims=True)
    valid = wet_sum.squeeze(-1) > 0      # at least one wet neighbour
    w_norm = np.zeros_like(w_wet)
    w_norm[valid] = w_wet[valid] / wet_sum[valid]

    log.info("Weights: %d/%d targets have ≥1 wet neighbour (%.1f%%)",
             valid.sum(), valid.size, 100*valid.mean())
    return idx, w_norm, valid


def apply_weights(arr, idx, w, valid, nlat, nlon):
    """arr: (ncell,) source. Returns (nlat, nlon) on the target grid."""
    src = np.where(np.isnan(arr) | (np.abs(arr) > SSH_PHYS_MAX), 0.0, arr).astype(np.float32)
    out = (src[idx] * w).sum(axis=1)
    out[~valid] = np.nan
    return out.reshape(nlat, nlon)


# ── 4. Save NetCDF ───────────────────────────────────────────────────
def save_nc(field, lat, lon, path):
    da = xr.DataArray(
        field, dims=("latitude", "longitude"),
        coords={"latitude": lat, "longitude": lon},
        name="adt",
        attrs={
            "long_name": "Sea surface height (proxy for ADT)",
            "units": "m",
            "source": f"DestinE Climate-DT CLTE o2d avg_zos via FDB ({MODEL} {EXPERIMENT})",
            "interpolation": "bilinear, wet-renormalized HEALPix→0.25°",
        },
    )
    ds = da.to_dataset()
    ds.attrs["date"] = DATE
    ds.to_netcdf(path)
    log.info("Wrote %s", path)


# ── 5. Plots ─────────────────────────────────────────────────────────
def plot_global_and_zoom(field, lat, lon, fig_global, fig_zoom):
    import matplotlib.pyplot as plt
    import cartopy.crs as ccrs
    import cartopy.feature as cfeature

    vmax = np.nanpercentile(np.abs(field), 99)

    # Global
    fig = plt.figure(figsize=(14, 7))
    ax = plt.axes(projection=ccrs.Robinson(central_longitude=0))
    pcm = ax.pcolormesh(lon, lat, field, transform=ccrs.PlateCarree(),
                        cmap="RdBu_r", vmin=-vmax, vmax=vmax, shading="auto")
    ax.add_feature(cfeature.LAND, facecolor="lightgray", zorder=1)
    ax.coastlines(linewidth=0.4, zorder=2)
    ax.set_global()
    plt.colorbar(pcm, ax=ax, orientation="horizontal", pad=0.04, shrink=0.7,
                 label="SSH (m), wet-aware bilinear")
    ax.set_title(f"avg_zos {MODEL} {EXPERIMENT} {DATE} — HEALPix→0.25° (CMEMS-aligned)")
    fig.savefig(fig_global, dpi=130, bbox_inches="tight")
    plt.close(fig)
    log.info("Wrote %s", fig_global)

    # Gulf Stream zoom — the canonical eddy-rich region
    fig = plt.figure(figsize=(11, 8))
    ax = plt.axes(projection=ccrs.PlateCarree())
    ax.set_extent([-82, -40, 25, 50], crs=ccrs.PlateCarree())
    # subset for cleaner pcolormesh extents
    lon_box = ((lon + 180) % 360) - 180  # 0-360 -> -180..180 for selection
    msk_lon = (lon_box >= -82) & (lon_box <= -40)
    msk_lat = (lat >= 25) & (lat <= 50)
    sub = field[np.ix_(msk_lat, msk_lon)]
    sub_lon = lon_box[msk_lon]
    sub_lat = lat[msk_lat]
    order = np.argsort(sub_lon)
    sub = sub[:, order]
    sub_lon = sub_lon[order]
    pcm = ax.pcolormesh(sub_lon, sub_lat, sub, transform=ccrs.PlateCarree(),
                        cmap="RdBu_r", vmin=-vmax, vmax=vmax, shading="auto")
    ax.add_feature(cfeature.LAND, facecolor="lightgray", zorder=1)
    ax.coastlines(linewidth=0.6, zorder=2)
    ax.gridlines(draw_labels=True, linewidth=0.3, alpha=0.5)
    plt.colorbar(pcm, ax=ax, orientation="vertical", pad=0.04, shrink=0.8,
                 label="SSH (m)")
    ax.set_title(f"Gulf Stream — coastal NaN handling check  ({DATE})")
    fig.savefig(fig_zoom, dpi=130, bbox_inches="tight")
    plt.close(fig)
    log.info("Wrote %s", fig_zoom)


def main():
    arr = fetch_zos_one_day()
    land = discover_land_mask(arr)
    np.savez(OUT_DIR / f"land_mask_nside{NSIDE}.npz", land=land, nside=np.int64(NSIDE))

    lat, lon = build_target_grid()
    log.info("Target: %d × %d (lat %.4f..%.4f, lon %.4f..%.4f)",
             len(lat), len(lon), lat[0], lat[-1], lon[0], lon[-1])

    idx, w, valid = build_wet_weights(NSIDE, lat, lon, land)
    np.savez(OUT_DIR / f"weights_ocean_nside{NSIDE}_{RES_DEG}deg.npz",
             indices=idx, weights=w, valid=valid,
             lat=lat, lon=lon, nside=np.int64(NSIDE),
             method="bilinear-wet-renorm")
    log.info("Saved weights file")

    field = apply_weights(arr, idx, w, valid, len(lat), len(lon))
    log.info("Interpolated stats: min=%g max=%g mean=%g  nan=%d (%.1f%%)",
             np.nanmin(field), np.nanmax(field), np.nanmean(field),
             int(np.isnan(field).sum()),
             100 * np.isnan(field).sum() / field.size)

    save_nc(field, lat, lon, OUT_DIR / f"zos_{DATE}_0.25deg.nc")

    plot_global_and_zoom(field, lat, lon,
                         FIG_DIR / f"zos_global_{DATE}.png",
                         FIG_DIR / f"zos_gulfstream_{DATE}.png")


if __name__ == "__main__":
    main()
