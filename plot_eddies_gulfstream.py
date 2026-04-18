"""Overlay EddyId output (Anticyclonic_*.nc, Cyclonic_*.nc) on the
interpolated SSH for the Gulf Stream box, as a sanity check that the
detected eddies sit on real SSH features.

Run in the eddy env:
    $EDDY_PYTHON plot_eddies_gulfstream.py
"""

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import xarray as xr
import cartopy.crs as ccrs
import cartopy.feature as cfeature

import eddy_paths

DATE = "20140101"
SSH = eddy_paths.out_dir() / f"adt_ifs-fesom_hist/adt_DestinE_ifs-fesom_hist_{DATE}.nc"
EDDY_DIR = eddy_paths.out_dir() / "eddyid_test"
A_NC = EDDY_DIR / f"Anticyclonic_{DATE}T000000.nc"
C_NC = EDDY_DIR / f"Cyclonic_{DATE}T000000.nc"
OUT = eddy_paths.fig_dir() / f"eddies_gulfstream_{DATE}.png"

# Gulf Stream box (-180..180 convention for plotting)
LON0, LON1, LAT0, LAT1 = -82, -40, 25, 50

ssh = xr.open_dataset(SSH)
field = ssh["adt"].isel(time=0).values
lat = ssh["latitude"].values
lon = ssh["longitude"].values
lon_pm = ((lon + 180) % 360) - 180
mlon = (lon_pm >= LON0) & (lon_pm <= LON1)
mlat = (lat >= LAT0) & (lat <= LAT1)
sub = field[np.ix_(mlat, mlon)]
sub_lon = lon_pm[mlon]
sub_lat = lat[mlat]
order = np.argsort(sub_lon)
sub = sub[:, order]
sub_lon = sub_lon[order]


def load_eddies(path):
    """Return arrays of (centre_lon[-180..180], centre_lat, effective_radius_m)."""
    ds = xr.open_dataset(path)
    print(f"{path.name}: {ds.dims}; vars={list(ds.data_vars)[:8]}...")
    clon = ds["longitude"].values
    clat = ds["latitude"].values
    rad = ds["effective_radius"].values  # effective radius (m)
    clon = ((clon + 180) % 360) - 180
    return clon, clat, rad


a_lon, a_lat, a_rad = load_eddies(A_NC)
c_lon, c_lat, c_rad = load_eddies(C_NC)

# Subset to box
def in_box(lon_, lat_):
    return (lon_ >= LON0) & (lon_ <= LON1) & (lat_ >= LAT0) & (lat_ <= LAT1)

am = in_box(a_lon, a_lat)
cm = in_box(c_lon, c_lat)
print(f"In Gulf Stream box: {am.sum()} anticyclonic, {cm.sum()} cyclonic")

vmax = float(np.nanpercentile(np.abs(field), 99))
fig = plt.figure(figsize=(13, 9))
ax = plt.axes(projection=ccrs.PlateCarree())
ax.set_extent([LON0, LON1, LAT0, LAT1], crs=ccrs.PlateCarree())
pcm = ax.pcolormesh(sub_lon, sub_lat, sub, transform=ccrs.PlateCarree(),
                    cmap="RdBu_r", vmin=-vmax, vmax=vmax, shading="auto")
ax.add_feature(cfeature.LAND, facecolor="lightgray", zorder=2)
ax.coastlines(linewidth=0.6, zorder=3)
ax.gridlines(draw_labels=True, linewidth=0.3, alpha=0.5)

# Eddy circles. Convert radius (m) to degrees latitude by /111e3.
for lo, la, r in zip(a_lon[am], a_lat[am], a_rad[am]):
    ax.add_patch(plt.Circle((lo, la), r / 111e3, transform=ccrs.PlateCarree(),
                              edgecolor="black", facecolor="none",
                              linewidth=0.9, zorder=4))
for lo, la, r in zip(c_lon[cm], c_lat[cm], c_rad[cm]):
    ax.add_patch(plt.Circle((lo, la), r / 111e3, transform=ccrs.PlateCarree(),
                              edgecolor="green", facecolor="none",
                              linewidth=0.9, zorder=4))

import matplotlib.patches as mpatches
ax.legend(handles=[
    mpatches.Patch(edgecolor="black", facecolor="none", label=f"Anticyclonic ({am.sum()})"),
    mpatches.Patch(edgecolor="green", facecolor="none", label=f"Cyclonic ({cm.sum()})"),
], loc="upper right")
plt.colorbar(pcm, ax=ax, orientation="vertical", pad=0.04, shrink=0.8,
             label="SSH (m)")
ax.set_title(f"DestinE IFS-FESOM hist {DATE} — EddyId on derived geostrophic u/v")

OUT.parent.mkdir(parents=True, exist_ok=True)
fig.savefig(OUT, dpi=130, bbox_inches="tight")
print(f"Wrote {OUT}")
