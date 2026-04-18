# eddy-destine-pipeline

End-to-end ocean-mesoscale-eddy pipeline for the
**[DestinE Climate Digital Twin](https://destine.ecmwf.int/)** ocean output
(IFS-FESOM, ICON, IFS-NEMO — historical, control, SSP3-7.0), wrapping
**[py-eddy-tracker](https://github.com/AntSimi/py-eddy-tracker)** (Mason et
al. 2014) for detection and tracking, with bespoke fetch/interp and
visualisation/diagnostic scripts.

What the pipeline does:

1. **Fetch** daily `avg_zos` (sea-surface height) from the DestinE FDB via
   [`fdb-xarray`](https://github.com/destination-earth-digital-twins/fdb-xarray),
   wet-aware HEALPix→regular lat/lon bilinear regrid (0.25° CMEMS-aligned
   grid), and write one NetCDF per day named for `EddyId`'s CLI.
2. **Detect** — call `EddyId` in parallel over the daily NetCDFs (geostrophic
   u/v are derived inside py-eddy-tracker from SSH).
3. **Track** — single `EddyTracking` pass per polarity over the whole window
   (tracks cross year boundaries correctly).
4. **Plot / diagnose** — summary figures per model + frame-by-frame video
   frames + two-model comparison plots (density, lifetime, radius,
   amplitude, latitudinal and seasonal distributions).

Empirical numbers on a 64-core box: a full **25-year run (1990–2014), one
model, 0.25° global** takes ≈ 8 h wall and ≈ 47 GB disk (17.8 GB ADT +
16.1 GB detections + 13.1 GB tracks).

---

## Repository layout

```
eddy-destine-pipeline/
├── README.md                     this file
├── LICENSE                       MIT
├── eddy_paths.py                 env-based path resolution helper
│
├── smoke_interp_zos.py           one-day smoke fetch + NetCDF + plot
├── download_interp_zos.py        multi-day fetch + wet-aware HEALPix→0.25° regrid
├── detect_eddies_batch.py        parallel EddyId over a dir of daily NetCDFs
│
├── pipeline_year.py              single-year orchestrator (fetch → detect → track → plot)
├── pipeline_multiyear.py         year-range orchestrator (parallel fetch, one tracking pass)
├── chain_icon_after_fesom.sh     polls FESOM PID, then launches ICON multi-year
│
├── plot_eddies_gulfstream.py     one-day SSH + detections sanity overlay
├── plot_tracks_jan2014.py        single-window trajectory plot (kept as example)
├── plot_year_summary.py          window-end summary (tracks, lifetimes, count/day, scoreboard)
├── plot_eddy_frames.py           per-day SSH + trailing-track frames for ffmpeg
├── plot_model_comparison.py      two-model comparison (density, hist, diff maps, scoreboard)
│
└── docs/
    ├── outputs.md                full product inventory
    └── troubleshooting.md        gotchas + detached-launch recipe + IFS-NEMO routing
```

---

## Requirements

### Hardware

| resource | single-year smoke | 25-yr / 0.25° global |
|---|---|---|
| CPU cores  | 8–16 | 32 (leaves headroom for chained model runs) |
| RAM        | ~5 GB | ~25–40 GB during tracking |
| Disk free  | ~5 GB | ~50 GB per model-window (ADT + eddies + tracks) |
| Network    | ~0.4 GB | ~10 GB per model-window |

Runs on smaller boxes — just lower `--workers` and expect longer walltime.

### Software

Two conda environments — kept separate because their pins disagree
(zarr 3.x for `fdb-xarray`, zarr 2.x for py-eddy-tracker; matplotlib < 3.8
for py-eddy-tracker):

1. **`dedl` env** — fetch + interpolation. DestinE Data Lake env, or any env
   with `z3fdb`, `fdb-xarray`, `healpy`, `xarray`, `zarr`, `numpy`, `cartopy`.
   On DestinE-provisioned hosts it's typically pre-installed at
   `/opt/miniconda3/envs/dedl/`.
2. **`eddy` env** — py-eddy-tracker. Built from the project's own
   `environment.yml`; ships `EddyId` and `EddyTracking` CLIs.

### FDB access

DestinE Climate-DT FDB via a `z3fdb` routing config (`~/config.yaml` by
default). **IFS-NEMO ocean may require a polytope-aware config** — smoke
test before queuing a long run, see `docs/troubleshooting.md`.

### py-eddy-tracker source

We build the eddy env from py-eddy-tracker's own environment.yml:

```bash
git clone https://github.com/AntSimi/py-eddy-tracker.git
```

### fdb-xarray

Not on PyPI; clone alongside the repo:

```bash
git clone https://github.com/destination-earth-digital-twins/fdb-xarray.git
```

---

## Installation

```bash
# 1. Clone this repo + py-eddy-tracker + fdb-xarray side-by-side.
git clone https://github.com/<you>/eddy-destine-pipeline.git
git clone https://github.com/AntSimi/py-eddy-tracker.git
git clone https://github.com/destination-earth-digital-twins/fdb-xarray.git

# 2. Build the eddy env (user prefix — /opt/miniconda3/envs is usually
#    root-owned on DestinE hosts).
mkdir -p $HOME/envs
cd py-eddy-tracker
/opt/miniconda3/envs/dedl/bin/mamba env create \
    --prefix $HOME/envs/eddy \
    -f environment.yml --yes

# 3. Verify both envs can import what they need.
/opt/miniconda3/envs/dedl/bin/python -c "
import sys; sys.path.insert(0, '$HOME/fdb-xarray')
import z3fdb, fdb_xarray, healpy, xarray, zarr, cartopy
print('fetch env OK, zarr', zarr.__version__)
"
$HOME/envs/eddy/bin/python -c "
import py_eddy_tracker; print('eddy env OK')
"
ls $HOME/envs/eddy/bin/EddyId $HOME/envs/eddy/bin/EddyTracking
```

---

## Configuration

All machine-specific paths come from environment variables with sensible
defaults (see `eddy_paths.py`). Set them once (e.g. in your shell rc or a
sourced `.env`):

```bash
# Python interpreters
export EDDY_DEDL_PYTHON=/opt/miniconda3/envs/dedl/bin/python
export EDDY_PYTHON=$HOME/envs/eddy/bin/python

# py-eddy-tracker console-script entry points (defaults derive from
# EDDY_PYTHON; override only if EddyId/EddyTracking live elsewhere)
# export EDDY_EDDYID_BIN=$HOME/envs/eddy/bin/EddyId
# export EDDY_TRACKING_BIN=$HOME/envs/eddy/bin/EddyTracking

# DestinE FDB routing config
export EDDY_FDB_CONFIG=$HOME/config.yaml

# Local fdb-xarray checkout (only if not pip-installed into EDDY_DEDL_PYTHON)
export EDDY_FDB_XARRAY=$HOME/fdb-xarray

# Output root — make sure it has enough free disk (see Requirements)
export EDDY_ROOT=/scratch/$USER/eddy
```

Each orchestrator also accepts `--root /alternate/path` per-invocation.
Resolution, model, experiment, parallelism are explicit CLI flags.

---

## Running the pipeline

### 1. Smoke test (~3 min)

Catches every dependency before you commit hours to a real run:

```bash
# A) Fetch one day, interp, write NetCDF and two PNGs
$EDDY_DEDL_PYTHON smoke_interp_zos.py
# Inspect $EDDY_ROOT/figs/zos_global_2014-01-01.png
#         $EDDY_ROOT/figs/zos_gulfstream_2014-01-01.png

# B) Detect eddies in the one day (eddy env)
$EDDY_PYTHON -c "
import eddy_paths, subprocess, sys
subprocess.check_call([
    str(eddy_paths.eddyid_bin()),
    str(eddy_paths.out_dir() / 'zos_2014-01-01_0.25deg.nc'), '20140101',
    'adt', 'None', 'None', 'longitude', 'latitude',
    '/tmp/eddy_smoke', '-v', 'INFO',
])"
ls /tmp/eddy_smoke/Anticyclonic_*.nc /tmp/eddy_smoke/Cyclonic_*.nc
```

Sanity: interpolated SSH range ~ ±1.2 m; ~29 % of cells land-masked;
several thousand anti/cyc detections globally.

### 2. One year (end-to-end)

```bash
$EDDY_DEDL_PYTHON pipeline_year.py \
    --year 2014 --model IFS-FESOM --experiment hist --workers 16
```

~30 min wall on a 64-core box. Produces:

- `$EDDY_ROOT/out/adt_ifs-fesom_hist/adt_DestinE_*_YYYYMMDD.nc`
- `$EDDY_ROOT/out/eddies_ifs-fesom_hist/{Anti,}cyclonic_YYYYMMDD.nc`
- `$EDDY_ROOT/out/tracks_ifs-fesom_hist_2014/{anti,cyc}/{Anti,}cyclonic.nc`
- `$EDDY_ROOT/figs/year_2014_ifs-fesom_hist/`

### 3. Multi-year (production)

Always launch detached so it survives logout:

```bash
cd $EDDY_ROOT && setsid nohup \
    $EDDY_DEDL_PYTHON $REPO/pipeline_multiyear.py \
        --start-year 1990 --end-year 2014 \
        --model IFS-FESOM --experiment hist \
        --workers 32 --fetch-parallel 8 \
    > $EDDY_ROOT/logs/multiyear_1990_2014_main.log 2>&1 < /dev/null & disown

# Capture the real (python) PID — `$!` points at the bash wrapper:
sleep 3 && pgrep -af pipeline_multiyear.py | grep -v grep | head -1 \
    | awk '{print $1}' > $EDDY_ROOT/logs/multiyear.pid
```

~8 h wall for 25 years at 32 workers. The tracking phase runs once over the
whole window, so tracks crossing New Year are continuous.

### 4. Chained runs (FESOM → ICON → IFS-NEMO)

`chain_icon_after_fesom.sh` polls the FESOM PID file and launches ICON when
the first run exits. Copy it, edit `--model` / PID paths for further chains:

```bash
setsid nohup bash $REPO/chain_icon_after_fesom.sh \
    > $EDDY_ROOT/logs/chain_icon.log 2>&1 < /dev/null & disown
```

### 5. Phase skipping

Every phase is idempotent — rerun the same orchestrator and it skips work
whose outputs already exist. Or skip phases explicitly:

| flag | skips |
|---|---|
| `--skip-fetch`   | FDB pull + regrid |
| `--skip-detect`  | parallel EddyId |
| `--skip-track`   | EddyTracking |
| `--no-plots`     | summary plots |
| `--force, -f`    | force rerun even if outputs exist |

### 6. Summary plots for one model

`pipeline_year.py` / `pipeline_multiyear.py` call this at the end; to
regenerate, do it manually:

```bash
$EDDY_PYTHON plot_year_summary.py \
    --label "ICON 1990-2014" \
    --track-dir $EDDY_ROOT/out/tracks_icon_hist_1990_2014 \
    --fig-dir   $EDDY_ROOT/figs/window_1990_2014_icon_hist
```

Produces track maps (global Robinson + 4 regions), lifetime histogram, daily
count time series, scoreboard.txt.

### 7. Two-model comparison (ICON vs IFS-FESOM)

```bash
$EDDY_PYTHON plot_model_comparison.py \
    --label-a ICON       --track-dir-a $EDDY_ROOT/out/tracks_icon_hist_1990_2014 \
    --label-b IFS-FESOM  --track-dir-b $EDDY_ROOT/out/tracks_ifs-fesom_hist_1990_2014 \
    --fig-dir $EDDY_ROOT/figs/compare_icon_vs_fesom_1990_2014
```

Produces 16 figures + scoreboard: per-model density / mean-amplitude /
mean-radius maps, A-B difference maps, overlaid histograms of lifetime /
radius / amplitude, latitudinal distribution, daily count, seasonal cycle.

### 8. Video frames + ffmpeg

Per-day SSH + trailing eddy tracks for animations:

```bash
# Single frame sanity check
$EDDY_PYTHON plot_eddy_frames.py --source ifs-fesom \
    --start 1992-01-01 --end 1992-01-01 --trail-days 60

# Full year, one region (all 4 regions in parallel is typical)
for r in global gulfstream agulhas kuroshio; do
    $EDDY_PYTHON plot_eddy_frames.py --source icon --region $r \
        --start 1992-01-01 --end 1992-12-31 --trail-days 60 &
done; wait

# ffmpeg composition
cd $EDDY_ROOT/figs/frames_icon_global_1992
ffmpeg -framerate 24 -pattern_type glob -i '*.png' \
    -vf "scale=trunc(iw/2)*2:trunc(ih/2)*2" \
    -c:v libx264 -pix_fmt yuv420p icon_global_1992.mp4
```

Regions available: `global` (Robinson), `gulfstream`, `agulhas`, `kuroshio`
(PlateCarree). Default SSH colour range is per-region
(Gulf Stream ±0.8 m, others ±1.2 m); override with `--vmax`.

---

## Output products

Under `$EDDY_ROOT/out/`:

- `adt_<model>_<experiment>/adt_DestinE_<model>_<exp>_YYYYMMDD.nc` —
  daily SSH interpolated to 0.25° CMEMS-aligned grid, variable `adt`,
  unit m. These are the direct inputs to `EddyId`.
- `eddies_<model>_<experiment>/{Anti,}cyclonic_YYYYMMDD.nc` — per-day
  detections (centre lat/lon, contours, radius, amplitude, …).
- `tracks_<model>_<experiment>_<START>_<END>/{anti,cyc}/` — track files
  produced by `EddyTracking`. `Anticyclonic.nc` / `Cyclonic.nc` are the
  **canonical track databases** (obs-level time / longitude / latitude /
  effective_radius / amplitude / track id / …).
- `weights_ocean_<model>_nside1024_0.25deg.npz` — wet-renormalized
  HEALPix→regular bilinear weights (one file per model; built on first
  fetch).
- `land_mask_<model>_nside1024.npz` — HEALPix ocean mask derived from
  the first fetched snapshot.

Under `$EDDY_ROOT/figs/`:

- `window_<start>_<end>_<model>_<experiment>/` — per-model summary.
- `compare_<model_a>_vs_<model_b>_<start>_<end>/` — two-model comparison.
- `frames_<source>_<region>_<yyyy>/frame_YYYYMMDD.png` — video frames.

Full inventory in [`docs/outputs.md`](docs/outputs.md).

---

## Design choices worth knowing

**Variable fetched.** `avg_zos` (CLTE o2d) is stored as `adt` (the unit is
m, unchanged). We only fetch SSH — py-eddy-tracker computes geostrophic
u/v internally via `grid.add_uv` when the u/v args to `EddyId` are literal
`"None"` strings (this is the documented CLI mode).

**Land sentinel.** DestinE ocean uses `+9999` (positive) on land cells, not
`-9999`. Masking uses `|SSH| > 50 m`, which is robust to any future
sentinel change (true SSH is ±2 m globally).

**HEALPix→regular weights are model-specific** because each ocean model has
a slightly different coastline mask (FESOM 29.0 % land, ICON 29.3 %).
Filenames embed the model name; do not share weights between models.

**Wet-renormalized bilinear.** The stock HEALPix→regular weights include
land neighbours. We zero out land neighbours and renormalize over the wet
set, so coastal cells don't average with zeros/sentinels.

**EddyTracking filename regex.** `EddyId` writes
`Anticyclonic_YYYYMMDDTHHMMSS.nc`, but `EddyTracking`'s hardcoded regex
only parses `_YYYYMMDD.nc`. `detect_eddies_batch.py` renames after each
`EddyId` call.

**YAML `FILES_PATTERN` as an explicit list.** `pipeline_multiyear.py`
builds the per-polarity file list explicitly rather than using a glob, so
eddy detections from prior years outside the tracking window don't get
pulled into a new window.

---

## Troubleshooting

See [`docs/troubleshooting.md`](docs/troubleshooting.md) for concrete
errors and fixes — land-sentinel + wrong-weights combo, IFS-NEMO polytope
routing, detached-launch recipe, `EddyTracking` filename regex.

---

## Contributing

PRs welcome. The pipeline is deliberately thin glue around py-eddy-tracker.
Changes that add *science* should usually go upstream to py-eddy-tracker;
changes that broaden supported data sources (more DestinE scenarios, other
reanalysis products) belong here.

## Citation

If you use this pipeline in a publication, please cite py-eddy-tracker:

> Mason, E., A. Pascual, and J. C. McWilliams (2014). A new sea surface
> height-based code for oceanic mesoscale eddy tracking. *J. Atmos. Ocean.
> Technol.*, **31**, 1181–1188.
> [doi:10.1175/JTECH-D-14-00019.1](https://doi.org/10.1175/JTECH-D-14-00019.1)

and the DestinE Climate Digital Twin documentation for the underlying data.

## License

MIT. See [`LICENSE`](LICENSE).
