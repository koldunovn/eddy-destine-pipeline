"""Fetch a window of daily avg_zos from DestinE, interpolate to 0.25° regular,
and write one per-day NetCDF ready for py-eddy-tracker's EddyId.

  python download_interp_zos.py --start-date 2014-01-01 --end-date 2014-01-10
  python download_interp_zos.py --start-date 2014-01-01 --days 31

Output NetCDF naming:
  <out>/adt_DestinE_<model_lower>_<experiment>_<YYYYMMDD>.nc
with variable name "adt" so the EddyId CLI works unchanged:
  EddyId adt_DestinE_ifs-fesom_hist_20140101.nc 20140101 adt None None longitude latitude ./

Reuses the static ocean weights file built by smoke_interp_zos.py (or
auto-builds it once on first run by sampling --start-date as the mask source).
"""

import argparse
import logging
import time
from datetime import datetime, timedelta
from pathlib import Path

import eddy_paths
eddy_paths.inject_fdb_xarray()

import healpy as hp
import numpy as np
import xarray as xr
from fdb_xarray import open_climate_dt

logging.basicConfig(level=logging.INFO,
                    format="%(asctime)s %(levelname)s %(message)s")
log = logging.getLogger("dl-interp")

DEFAULT_CONFIG = str(eddy_paths.fdb_config())
DEFAULT_OUT = eddy_paths.out_dir()
NSIDE = 1024
RES_DEG = 0.25
SSH_PHYS_MAX = 50.0


# ── Target grid (CMEMS-aligned cell-centred 0.25°) ───────────────────
def build_target_grid(res_deg=RES_DEG):
    half = res_deg / 2.0
    lat = np.arange(-90.0 + half, 90.0, res_deg).astype(np.float64)
    lon = np.arange(0.0 + half, 360.0, res_deg).astype(np.float64)
    return lat, lon


# ── Weights file management ──────────────────────────────────────────
# Weights are model-specific because each ocean model uses a slightly
# different land mask on the HEALPix grid (FESOM ~29.0% land, ICON ~29.3%).
# Reusing one model's weights for another silently mis-masks the ocean.
def weights_path(out_dir: Path, model: str) -> Path:
    return out_dir / f"weights_ocean_{model.lower()}_nside{NSIDE}_{RES_DEG}deg.npz"


def land_mask_path(out_dir: Path, model: str) -> Path:
    return out_dir / f"land_mask_{model.lower()}_nside{NSIDE}.npz"


def build_wet_weights(land_mask, lat, lon):
    lon2d, lat2d = np.meshgrid(lon, lat)
    theta = np.deg2rad(90.0 - lat2d.ravel())
    phi = np.deg2rad(lon2d.ravel())
    idx, w = hp.get_interp_weights(NSIDE, theta, phi, nest=True)
    idx = idx.T.astype(np.int64)
    w = w.T.astype(np.float32)
    wet = ~land_mask[idx]
    w_wet = w * wet
    wet_sum = w_wet.sum(axis=1, keepdims=True)
    valid = wet_sum.squeeze(-1) > 0
    w_norm = np.zeros_like(w_wet)
    w_norm[valid] = w_wet[valid] / wet_sum[valid]
    return idx, w_norm, valid


def ensure_weights(out_dir, model, fetch_one_snapshot):
    """Load model-specific weights file; build it from a fresh FDB snapshot
    if missing. Each (out_dir, model) pair has its own weights file."""
    wpath = weights_path(out_dir, model)
    if wpath.exists():
        log.info("Loading existing weights: %s", wpath.name)
        wz = np.load(wpath, allow_pickle=False)
        return (wz["indices"], wz["weights"], wz["valid"],
                wz["lat"], wz["lon"])

    log.info("Weights file missing for model=%s — sampling one snapshot to build land mask", model)
    arr = fetch_one_snapshot()
    land = np.isnan(arr) | (np.abs(arr) > SSH_PHYS_MAX)
    log.info("Land mask: %d/%d cells (%.1f%%)", land.sum(), land.size,
             100 * land.mean())
    np.savez(land_mask_path(out_dir, model), land=land, nside=np.int64(NSIDE))

    lat, lon = build_target_grid()
    idx, w, valid = build_wet_weights(land, lat, lon)
    log.info("Weights: %d/%d targets wet (%.1f%%)",
             valid.sum(), valid.size, 100 * valid.mean())
    np.savez(wpath, indices=idx, weights=w, valid=valid,
             lat=lat, lon=lon, nside=np.int64(NSIDE),
             method="bilinear-wet-renorm", model=model.lower())
    return idx, w, valid, lat, lon


# ── FDB fetch ────────────────────────────────────────────────────────
def fetch_window(model, experiment, resolution, start_iso, end_iso, config):
    log.info("FDB open avg_zos %s..%s", start_iso, end_iso)
    ds = open_climate_dt(
        portfolio="CLTE",
        levtype="o2d",
        model=model,
        experiment=experiment,
        start_date=start_iso,
        end_date=end_iso,
        resolution=resolution,
        variables=["avg_zos"],
        config_path=config,
    )
    da = ds["avg_zos"].sel(scenario=experiment)  # (time, cell)
    return da, ds["time"].values


def fetch_one_day_snapshot(model, experiment, resolution, date_iso, config):
    da, _ = fetch_window(model, experiment, resolution,
                          date_iso, date_iso, config)
    return da.isel(time=0).values.astype(np.float32)


# ── Interp + write ───────────────────────────────────────────────────
def apply_weights(arr, idx, w, valid, nlat, nlon):
    src = np.where(np.isnan(arr) | (np.abs(arr) > SSH_PHYS_MAX),
                   0.0, arr).astype(np.float32)
    out = (src[idx] * w).sum(axis=1)
    out[~valid] = np.nan
    return out.reshape(nlat, nlon)


def write_day_nc(field, lat, lon, time_value, path, model, experiment):
    da = xr.DataArray(
        field[None, :, :], dims=("time", "latitude", "longitude"),
        coords={"time": [time_value], "latitude": lat, "longitude": lon},
        name="adt",
        attrs={
            "long_name": "Sea surface height (treated as ADT for eddy detection)",
            "units": "m",
            "standard_name": "sea_surface_height_above_geoid",
            "source": (f"DestinE Climate-DT CLTE o2d avg_zos via FDB "
                       f"({model} {experiment})"),
            "interpolation": "bilinear, wet-renormalized HEALPix→0.25°",
        },
    )
    ds = da.to_dataset()
    ds["latitude"].attrs.update(units="degrees_north", standard_name="latitude")
    ds["longitude"].attrs.update(units="degrees_east", standard_name="longitude")
    enc = {"adt": {"zlib": True, "complevel": 4,
                    "_FillValue": np.float32(np.nan)}}
    ds.to_netcdf(path, encoding=enc)


def out_name(model, experiment, day_iso):
    yyyymmdd = day_iso.replace("-", "")
    return f"adt_DestinE_{model.lower()}_{experiment}_{yyyymmdd}.nc"


# ── Main loop ────────────────────────────────────────────────────────
def parse_args():
    ap = argparse.ArgumentParser()
    ap.add_argument("--start-date", required=True, help="YYYY-MM-DD")
    g = ap.add_mutually_exclusive_group()
    g.add_argument("--end-date", help="YYYY-MM-DD inclusive")
    g.add_argument("--days", type=int, help="Number of days from start (incl)")
    ap.add_argument("--model", default="IFS-FESOM")
    ap.add_argument("--experiment", default="hist")
    ap.add_argument("--resolution", default="high")
    ap.add_argument("--out-dir", type=Path, default=DEFAULT_OUT)
    ap.add_argument("--config", default=DEFAULT_CONFIG)
    ap.add_argument("--force", "-f", action="store_true")
    return ap.parse_args()


def main():
    args = parse_args()
    args.out_dir.mkdir(parents=True, exist_ok=True)

    start = datetime.strptime(args.start_date, "%Y-%m-%d")
    if args.end_date:
        end = datetime.strptime(args.end_date, "%Y-%m-%d")
    elif args.days:
        end = start + timedelta(days=args.days - 1)
    else:
        end = start
    n_days = (end - start).days + 1
    log.info("Window: %s .. %s (%d days)  model=%s exp=%s res=%s",
             start.date(), end.date(), n_days, args.model,
             args.experiment, args.resolution)

    # Per-day output dir (NetCDFs go in a subdir keyed by model/experiment).
    nc_dir = args.out_dir / f"adt_{args.model.lower()}_{args.experiment}"
    nc_dir.mkdir(parents=True, exist_ok=True)

    # Idempotency: which days are missing?
    missing = []
    for k in range(n_days):
        d = (start + timedelta(days=k)).strftime("%Y-%m-%d")
        path = nc_dir / out_name(args.model, args.experiment, d)
        if args.force or not path.exists():
            missing.append((d, path))
    if not missing:
        log.info("Nothing to do — every day already present in %s", nc_dir)
        return
    log.info("Missing %d/%d days", len(missing), n_days)

    # Weights (build on first run from start_date snapshot).
    def _sample():
        return fetch_one_day_snapshot(args.model, args.experiment,
                                      args.resolution, args.start_date,
                                      args.config)
    idx, w, valid, lat, lon = ensure_weights(args.out_dir, args.model, _sample)
    nlat, nlon = len(lat), len(lon)

    # Single FDB request spanning the whole window — much cheaper than
    # one request per day. Then iterate through the time axis.
    da, times = fetch_window(args.model, args.experiment, args.resolution,
                              args.start_date, end.strftime("%Y-%m-%d"),
                              args.config)
    log.info("Window size on FDB side: %d timesteps × %d cells",
             da.sizes["time"], da.sizes["cell"])

    # Map iso-date → output path so we can skip already-present days even
    # though the FDB request is contiguous.
    target = {d: p for d, p in missing}
    t0 = time.perf_counter()
    written = skipped = 0
    for ti, tval in enumerate(times):
        day_iso = str(np.datetime64(tval, "D"))
        if day_iso not in target:
            skipped += 1
            continue
        snap = da.isel(time=ti).values.astype(np.float32)
        field = apply_weights(snap, idx, w, valid, nlat, nlon)
        write_day_nc(field, lat, lon, tval, target[day_iso],
                     args.model, args.experiment)
        written += 1
        log.info("  [%d/%d] %s  min=%.3f max=%.3f nan=%.1f%%",
                 written, len(target), day_iso,
                 float(np.nanmin(field)), float(np.nanmax(field)),
                 100 * float(np.isnan(field).mean()))
    dt = time.perf_counter() - t0
    log.info("Wrote %d days, skipped %d already-present, in %.1f s (%.2f s/day)",
             written, skipped, dt, dt / max(written, 1))


if __name__ == "__main__":
    main()
