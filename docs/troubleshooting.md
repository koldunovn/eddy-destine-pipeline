# Troubleshooting

Concrete errors we hit and how we fixed them. If you add a new one, keep
it in this table format.

| symptom | likely cause | fix |
|---|---|---|
| `ModuleNotFoundError: fdb_xarray` in dedl env | local `fdb-xarray/` not on host, or path unset | `export EDDY_FDB_XARRAY=/path/to/fdb-xarray` — the scripts inject it via `eddy_paths.inject_fdb_xarray()` |
| Fetch returns "no data found" for IFS-NEMO | the databridge endpoints typically don't serve IFS-NEMO ocean | swap to a polytope-aware `config.yaml`; see §IFS-NEMO below |
| `EddyTracking` error: "Several timesteps in grid dataset [0 372 …]" | filenames are `_YYYYMMDDT000000.nc` instead of bare `_YYYYMMDD.nc`; the `EddyTracking` date regex only matches the short form | run via `detect_eddies_batch.py` (it auto-renames), or `cd dir && for f in *T000000.nc; do mv "$f" "${f/T000000/}"; done` |
| Output `max=9999` and `nan=0%` after interp | land sentinel (`+9999` for FESOM ocean) not masked | confirm `SSH_PHYS_MAX = 50.0` mask is active in `download_interp_zos.py` / `smoke_interp_zos.py` — the physical-range mask is robust to sentinel-value changes |
| All target cells NaN after interp | wrong weights file (FESOM weights applied to ICON data) | delete `weights_ocean_<wrongmodel>_*.npz`, rerun — the script auto-builds per-model weights on first fetch |
| Pipeline dies on logout | not detached | always `setsid nohup … > log 2>&1 < /dev/null & disown` and verify with `ps -o pid,stat,tt` — you want `STAT=Ss, TT=?` |
| Detached launch: `$!` points at bash, not python | bash wrapper | after launch, `sleep 3 && pgrep -af pipeline_multiyear.py | grep -v grep | head -1 | awk '{print $1}' > $EDDY_ROOT/logs/multiyear.pid` |
| `multiyear_detect_<lo>_<hi>.log` collides between chained model runs | both runs use the same start/end-year | cosmetic — per-day logs in `eddies_<model>_hist/logs/` are unique; rename the main log if you want clean history |
| `EddyTracking` picks up eddy files outside the tracking window | `FILES_PATTERN` uses a glob | `pipeline_multiyear.py` builds an explicit file list per polarity — don't regress to glob |
| `ValueError` during `EddyTracking` tracking phase about missing frame | usually one-off gap in the daily detections (empty polarity file) | verify `ls eddies_<model>_hist/*_YYYYMMDD.nc | wc -l` equals `window_days * 2`; rerun the missing date's `EddyId` |
| Frame generation slow on global view | `pcolormesh` over 720 × 1440 per frame + thousands of line artists | expected — global frame is ~5× slower than regional. Parallelise across regions (each region its own process); don't try multithreading inside a single plot |
| `matplotlib` version conflicts after `mamba env create` | py-eddy-tracker pins matplotlib < 3.8 | don't upgrade matplotlib in the eddy env; keep the pin from py-eddy-tracker's environment.yml |

---

## IFS-NEMO FDB routing

The standard `config.yaml` (FESOM, ICON) routes to the databridge:

```yaml
type: select
fdbs:
  - select: class=d1,dataset=^climate-dt$,generation=2
    type: remote
    host: databridge-prod-catalogue4-ope.ewctest.link
    port: 10000
    engine: remote
    store: remote
  # ...generation=1, extremes-DT, etc.
```

IFS-NEMO ocean typically requires the polytope endpoint:

```
polytope.mn5.apps.dte.destination-earth.eu
```

**Before queuing a multi-decade IFS-NEMO run, smoke-test:**

```bash
$EDDY_DEDL_PYTHON download_interp_zos.py \
    --start-date 2014-01-01 --days 1 \
    --model IFS-NEMO --experiment hist --resolution high \
    --out-dir /tmp/eddy_nemo_smoke
```

- Succeeds in ~a few seconds → databridge serves IFS-NEMO on this host;
  proceed normally.
- "no data found" / hang / error → swap to a polytope-aware config. See
  `fdb-xarray/config.yaml.example` or ask the DestinE FDB admins for the
  exact polytope stanza (keys are site-specific).

The pipeline scripts don't care which endpoint serves the data — they
read `~/config.yaml` (or `$EDDY_FDB_CONFIG`) and let `z3fdb` handle routing.

---

## Detached launch recipe

```bash
# 1. Launch, redirect both streams, detach from shell session.
cd $EDDY_ROOT && setsid nohup \
    $EDDY_DEDL_PYTHON $REPO/pipeline_multiyear.py \
        --start-year 1990 --end-year 2014 \
        --model IFS-FESOM --experiment hist \
        --workers 32 --fetch-parallel 8 \
    > $EDDY_ROOT/logs/multiyear_1990_2014_main.log 2>&1 < /dev/null &
disown

# 2. Capture the real python PID (NOT $!).
sleep 3
pgrep -af pipeline_multiyear.py | grep -v grep | head -1 \
    | awk '{print $1}' > $EDDY_ROOT/logs/multiyear.pid

# 3. Verify it's detached.
ps -o pid,stat,tt,cmd -p "$(cat $EDDY_ROOT/logs/multiyear.pid)"
# Expected: STAT=Ss, TT=?
```

---

## Known gotchas (saved here so they don't bite again)

1. **Land sentinel is `+9999`** (positive), not `-9999`. Our masking uses
   `|SSH| > 50`, robust to sign and future sentinel changes.

2. **Weights files are model-specific.** FESOM/ICON have slightly
   different ocean masks (29.0 % vs 29.3 % land). Filenames embed the
   model name; don't share across models.

3. **`EddyId` → `EddyTracking` filename rename.** `_YYYYMMDDT000000.nc`
   → `_YYYYMMDD.nc`. Handled automatically in `detect_eddies_batch.py`.

4. **YAML `FILES_PATTERN` as a list.** Prevents eddy detections from
   prior runs outside the tracking window from being pulled in.

5. **Detached launch `$!` ≠ python PID.** `setsid nohup bash-wrapper &`
   gives the bash PID. Always grep for the real one after a few seconds.

6. **Do not upgrade matplotlib / zarr in the eddy env.** py-eddy-tracker
   pins these at old versions; upgrading breaks `EddyId`.

7. **The dedl env and the eddy env cannot be merged.** zarr 3.x (fetch)
   and zarr 2.x (py-eddy-tracker) are incompatible; keep the two envs
   separate as designed.
