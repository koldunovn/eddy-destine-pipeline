# Output product inventory

Paths below are relative to `$EDDY_ROOT` (default `/bigdisk/work/eddy`).

## `out/adt_<model>_<experiment>/`

Daily ADT NetCDFs, one file per day.

```
adt_DestinE_<model>_<experiment>_YYYYMMDD.nc
```

| dimension / variable | shape              | unit | notes |
|---|---|---|---|
| `time`               | (1,)               | —    | `datetime64[ns]`, 00:00 UTC |
| `latitude`           | (720,)             | deg  | cell-centred, ascending, −89.875 … 89.875 |
| `longitude`          | (1440,)            | deg  | cell-centred, 0.125 … 359.875 (0…360 convention) |
| `adt`                | (time, lat, lon)   | m    | SSH on 0.25° CMEMS-aligned grid |

`adt.attrs` contains `source` (DestinE FDB pointer) and `interpolation:
"bilinear, wet-renormalized HEALPix→0.25°"`. NaN = land / no wet neighbour.

Empirical size: 3.6 MB/file; full year ≈ 1.3 GB; 25 yr ≈ 33 GB.

## `out/eddies_<model>_<experiment>/`

Per-day eddy detections from `EddyId`.

```
Anticyclonic_YYYYMMDD.nc
Cyclonic_YYYYMMDD.nc
```

Variables: `longitude`, `latitude`, `amplitude`, `effective_radius`,
`speed_radius`, `speed_average`, contour arrays (`effective_contour_*`,
`speed_contour_*`), `uavg_profile`, `time`, plus CF-standard metadata.

One file per polarity per day. Empirical size: ~0.2–1 MB/file depending
on eddy count; full year ≈ 1.2 GB per polarity; 25 yr ≈ 30 GB total.

Filenames originally written as `_YYYYMMDDT000000.nc` by `EddyId`; our
`detect_eddies_batch.py` strips the `T000000` suffix so `EddyTracking`
can parse them.

## `out/tracks_<model>_<experiment>_<START>_<END>/`

Output of `EddyTracking`. Two sub-directories, one per polarity.

```
anti/
  Anticyclonic.nc                 ← canonical track database
  Anticyclonic_correspondances.nc
  Anticyclonic_track_too_short.nc ← tracks filtered out by min_lifetime (we use 4 d)
  Anticyclonic_untracked.nc       ← obs that didn't attach to a track
cyc/  (same filenames with Cyclonic_ prefix)
track_anti.yaml                    ← generated EddyTracking config
track_cyc.yaml
```

**The canonical track database is `Anticyclonic.nc` / `Cyclonic.nc`.**

Dimensions: `obs` (one row per track-day), `NbSample` = 50 (contour
vertices). Key variables:

| variable                     | unit  | notes |
|---|---|---|
| `track`                      | —     | uint32 track id (non-unique across polarities; see below) |
| `time`                       | —     | datetime64[ns] |
| `longitude` / `latitude`     | deg   | eddy centroid, 0..360 |
| `amplitude`                  | m     | max − min SSH of speed-based contour |
| `effective_radius`           | m     | circle-equivalent radius of effective contour |
| `speed_radius`               | m     | circle-equivalent radius of max-speed contour |
| `speed_average`              | m/s   | mean current speed at max-speed contour |
| `effective_contour_*`        | deg/m | 50-vertex effective contour |
| `speed_contour_*`            | deg/m | 50-vertex max-speed contour |
| `uavg_profile`               | m/s   | mean |u|(radius) profile |
| `observation_flag`           | —     | 0 = clean detection, 1 = gap-filled |
| `cost_association`           | —     | tracker link cost |
| `inner_contour_height`       | m     | innermost closed contour height |
| `effective_contour_height`   | m     | effective contour height |
| `effective_area`, `speed_area` | m²  | contour areas |

Track ids in `Anticyclonic.nc` and `Cyclonic.nc` are **independent
numberings** — an id overlap between files is coincidental, not a match.

Empirical obs counts at 25 yr global 0.25°:
- IFS-FESOM: 21.7 M anti + 23.4 M cyc obs; ≈ 1.1 M anti + 1.2 M cyc
  tracks ≥ 4 d; ≈ 195 k + 206 k tracks ≥ 30 d; longest anti 1813 d,
  cyc 949 d; median lifetime 9 d.
- ICON: 25.0 M + 27.3 M obs; 1.23 M + 1.38 M tracks ≥ 4 d; 222 k + 237 k
  ≥ 30 d; longest 1746 / 1587 d; median lifetime 10 d.

File sizes: `Anticyclonic.nc` ≈ 6 GB per 25-yr model.

## `out/weights_ocean_<model>_nside1024_0.25deg.npz`

Bilinear HEALPix→regular weights with land neighbours zeroed and
renormalized over wet neighbours:

- `indices`   (N_target × 4, int64)  — source HEALPix cells per target
- `weights`   (N_target × 4, float32) — renormalized weights
- `valid`     (N_target, bool)        — True where ≥ 1 wet neighbour
- `lat`, `lon` — target grid definition
- `nside`     — HEALPix nside (1024 for high-res Climate-DT)
- `method`    — `"bilinear-wet-renorm"`

~35 MB. Auto-built on first fetch if missing; reused across subsequent
fetches for the same model.

## `out/land_mask_<model>_nside1024.npz`

Ocean mask in HEALPix space derived from the first successful snapshot:

- `land` (N_cell, bool) — True for land / sentinel-filled cells
- `nside`                — 1024

FESOM ≈ 29.0 % land, ICON ≈ 29.3 %.

## `figs/window_<start>_<end>_<model>_<experiment>/`

Per-model summary from `plot_year_summary.py`:

| file | what |
|---|---|
| `tracks_global.png`   | long-lived tracks on Robinson map |
| `tracks_gulfstream.png`, `_agulhas.png`, `_kuroshio.png`, `_brazil.png` | regional zooms |
| `lifetime_histogram.png` | log-scale lifetime PDF, both polarities |
| `count_per_day.png`      | daily obs count on a real date axis |
| `scoreboard.txt`         | counts / medians / window length |

## `figs/compare_<a>_vs_<b>_<start>_<end>/`

Two-model comparison from `plot_model_comparison.py`:

| file | what |
|---|---|
| `density_{a,b}.png`            | observation density, log scale, anti + cyc |
| `mean_amplitude_{a,b}.png`     | spatial mean amplitude (≥ 5 obs per bin) |
| `mean_radius_{a,b}.png`        | spatial mean effective radius |
| `density_diff.png`             | A − B, diverging scale |
| `mean_amplitude_diff.png`      | Δ mean amplitude |
| `mean_radius_diff.png`         | Δ mean radius |
| `hist_lifetime.png`            | overlaid track-lifetime histograms |
| `hist_radius.png`              | overlaid obs-level radius histograms |
| `hist_amplitude.png`           | overlaid obs-level amplitude histograms |
| `lat_histogram.png`            | obs per 1° latitude band |
| `count_per_day.png`            | daily count overlay |
| `season_cycle.png`             | mean obs per month-of-year |
| `scoreboard.txt`               | side-by-side medians + counts |

## `figs/frames_<source>_<region>_<yyyy>/`

Per-day video frames from `plot_eddy_frames.py`:

```
frame_YYYYMMDD.png    (one per day in the requested window)
<optional ffmpeg-produced .mp4>
```

Global view is Robinson; regional (`gulfstream`, `agulhas`, `kuroshio`) is
PlateCarree. Trail length is `--trail-days` (default 60).
