"""End-to-end ocean-eddy pipeline for a multi-year window of DestinE data.

Differs from pipeline_year.py in three ways:

  1. Fetch is parallelised across (year, month) tasks (default 8 streams).
  2. Detection is one big parallel pool across every day in the window
     (the per-day NetCDFs are pooled in a single input dir, so detection
     doesn't even need to know which year a day belongs to).
  3. Tracking is ONE continuous EddyTracking call per polarity over the
     full multi-year glob — captures eddies that cross New-Year boundaries
     instead of artificially truncating them at Dec 31.

Usage:
    python pipeline_multiyear.py --start-year 1990 --end-year 2014 \\
        --model IFS-FESOM --workers 32 --fetch-parallel 8

Idempotent at every phase: re-runs skip work whose outputs already exist.
"""

import argparse
import calendar
import logging
import multiprocessing as mp
import subprocess
import sys
import time
from pathlib import Path

import eddy_paths

HERE = Path(__file__).resolve().parent
DEDL_PYTHON = str(eddy_paths.dedl_python())
EDDY_PYTHON = str(eddy_paths.eddy_python())
EDDY_TRACKING_BIN = str(eddy_paths.eddytracking_bin())

DEFAULT_ROOT = eddy_paths.default_root()

logging.basicConfig(level=logging.INFO,
                    format="%(asctime)s %(levelname)s %(message)s")
log = logging.getLogger("multiyear")


# ── helpers ─────────────────────────────────────────────────────────
def parse_args():
    ap = argparse.ArgumentParser()
    ap.add_argument("--start-year", type=int, required=True)
    ap.add_argument("--end-year", type=int, required=True,
                    help="Inclusive end year")
    ap.add_argument("--model", default="IFS-FESOM")
    ap.add_argument("--experiment", default="hist")
    ap.add_argument("--resolution", default="high")
    ap.add_argument("--root", type=Path, default=DEFAULT_ROOT)
    ap.add_argument("--workers", type=int, default=32,
                    help="Parallel EddyId workers (default 32)")
    ap.add_argument("--fetch-parallel", type=int, default=8,
                    help="Parallel FDB month-fetches (default 8)")
    ap.add_argument("--skip-fetch", action="store_true")
    ap.add_argument("--skip-detect", action="store_true")
    ap.add_argument("--skip-track", action="store_true")
    ap.add_argument("--no-plots", action="store_true")
    ap.add_argument("--force", "-f", action="store_true")
    return ap.parse_args()


def fmt_dt(secs):
    return f"{secs:.0f}s ({secs/60:.1f} min, {secs/3600:.2f} h)"


def run_logged(cmd, log_path: Path, label: str) -> int:
    log_path.parent.mkdir(parents=True, exist_ok=True)
    log.info("%s -> %s", label, log_path.name)
    t0 = time.perf_counter()
    with log_path.open("w") as f:
        f.write(" ".join(str(c) for c in cmd) + "\n\n")
        f.flush()
        proc = subprocess.Popen(cmd, stdout=subprocess.PIPE,
                                 stderr=subprocess.STDOUT, text=True, bufsize=1)
        for line in proc.stdout:
            f.write(line); f.flush()
        proc.wait()
    dt = time.perf_counter() - t0
    log.info("  %s rc=%d in %s", label, proc.returncode, fmt_dt(dt))
    return proc.returncode


# ── phase 1 — parallel fetch ────────────────────────────────────────
def _fetch_one_month(task):
    """task = (year, month, model, experiment, resolution, out_root, logs, force)"""
    (year, month, model, experiment, resolution, out_root, logs, force) = task
    last = calendar.monthrange(year, month)[1]
    start_iso = f"{year}-{month:02d}-01"
    end_iso   = f"{year}-{month:02d}-{last:02d}"
    cmd = [DEDL_PYTHON, str(HERE / "download_interp_zos.py"),
           "--start-date", start_iso, "--end-date", end_iso,
           "--model", model,
           "--experiment", experiment,
           "--resolution", resolution,
           "--out-dir", str(out_root)]
    if force:
        cmd.append("--force")
    log_path = Path(logs) / f"fetch_{year}_{month:02d}.log"
    t0 = time.perf_counter()
    with log_path.open("w") as f:
        f.write(" ".join(cmd) + "\n\n"); f.flush()
        rc = subprocess.run(cmd, stdout=f, stderr=subprocess.STDOUT).returncode
    return year, month, rc, time.perf_counter() - t0


def fetch_phase(args, paths) -> int:
    tasks = []
    for year in range(args.start_year, args.end_year + 1):
        for month in range(1, 13):
            tasks.append((year, month, args.model, args.experiment,
                          args.resolution, str(paths["out_root"]),
                          str(paths["logs"]), args.force))

    log.info("Fetch: %d (year, month) tasks across %d FDB streams",
             len(tasks), args.fetch_parallel)
    fail = 0
    t0 = time.perf_counter()
    ctx = mp.get_context("spawn")
    with ctx.Pool(args.fetch_parallel) as pool:
        done = 0
        for year, month, rc, dt in pool.imap_unordered(_fetch_one_month, tasks):
            done += 1
            tag = "ok " if rc == 0 else "FAIL"
            log.info("  [%d/%d] %s fetch %d-%02d  rc=%d  %.1fs",
                     done, len(tasks), tag, year, month, rc, dt)
            if rc != 0:
                fail += 1
    log.info("Fetch wall: %s  (%d failed of %d)",
             fmt_dt(time.perf_counter() - t0), fail, len(tasks))
    return fail


# ── phase 2 — parallel detect (single pool, all years) ──────────────
def detect_phase(args, paths) -> int:
    rc = run_logged(
        [EDDY_PYTHON, str(HERE / "detect_eddies_batch.py"),
         "--input-dir",  str(paths["adt_dir"]),
         "--output-dir", str(paths["eddy_dir"]),
         "--workers",    str(args.workers),
         *(["--force"] if args.force else [])],
        paths["logs"] / f"multiyear_detect_{args.start_year}_{args.end_year}.log",
        label=f"detect {args.start_year}..{args.end_year} ({args.workers} workers)",
    )
    return 0 if rc == 0 else 1


# ── phase 3 — single continuous tracking ────────────────────────────
def _build_file_list(eddy_dir: Path, sign_long: str, year_lo: int, year_hi: int):
    """Explicit list of detection NetCDFs whose YYYYMMDD lies in [year_lo, year_hi].

    Building the list in Python (rather than relying on a glob pattern) keeps
    us safe if the eddy_dir contains files from years outside the requested
    window — which happens whenever pipeline_year was used previously."""
    out = []
    for path in sorted(eddy_dir.glob(f"{sign_long}_????????.nc")):
        # Filename: <sign>_YYYYMMDD.nc — extract YYYY from the 8-digit segment.
        try:
            yyyy = int(path.stem.rsplit("_", 1)[1][:4])
        except (ValueError, IndexError):
            continue
        if year_lo <= yyyy <= year_hi:
            out.append(str(path))
    return out


def track_phase(args, paths) -> int:
    import yaml as _yaml
    paths["track_dir"].mkdir(parents=True, exist_ok=True)
    fail = 0
    for sign_short, sign_long in [("anti", "Anticyclonic"), ("cyc", "Cyclonic")]:
        save_dir = paths["track_dir"] / sign_short
        save_dir.mkdir(parents=True, exist_ok=True)
        out_nc = save_dir / f"{sign_long}.nc"
        if out_nc.exists() and not args.force:
            log.info("track %s already done -> %s", sign_short, out_nc.name)
            continue
        files = _build_file_list(paths["eddy_dir"], sign_long,
                                  args.start_year, args.end_year)
        if not files:
            log.error("No %s detection files found for %d..%d in %s",
                      sign_long, args.start_year, args.end_year, paths["eddy_dir"])
            fail += 1
            continue
        log.info("Tracking %s on %d files (%s .. %s)",
                 sign_long, len(files),
                 Path(files[0]).name, Path(files[-1]).name)
        # YAML list under FILES_PATTERN is the supported way to pass an explicit
        # file list (py-eddy-tracker's track() handles list vs string already).
        yaml_path = paths["track_dir"] / f"track_{sign_short}.yaml"
        cfg = {
            "PATHS": {
                "FILES_PATTERN": files,
                "SAVE_DIR": str(save_dir),
            },
            "VIRTUAL_LENGTH_MAX": 0,
            "TRACK_DURATION_MIN": 4,
            "CLASS": {
                "MODULE": "py_eddy_tracker.featured_tracking.area_tracker",
                "CLASS": "AreaTracker",
            },
        }
        yaml_path.write_text(_yaml.safe_dump(cfg, sort_keys=False))
        rc = run_logged(
            [EDDY_TRACKING_BIN, str(yaml_path), "-v", "INFO"],
            paths["logs"] / f"multiyear_track_{sign_short}_"
                            f"{args.start_year}_{args.end_year}.log",
            label=f"track {sign_long} {args.start_year}..{args.end_year}",
        )
        if rc != 0:
            fail += 1
    return fail


# ── phase 4 — plots ─────────────────────────────────────────────────
def plot_phase(args, paths) -> int:
    label = f"{args.model} {args.experiment} {args.start_year}-{args.end_year}"
    rc = run_logged(
        [EDDY_PYTHON, str(HERE / "plot_year_summary.py"),
         "--label",     label,
         "--track-dir", str(paths["track_dir"]),
         "--fig-dir",   str(paths["figs"] /
                            f"window_{args.start_year}_{args.end_year}_"
                            f"{args.model.lower()}_{args.experiment}")],
        paths["logs"] / f"multiyear_plots_{args.start_year}_{args.end_year}.log",
        label="plots",
    )
    return 0 if rc == 0 else 1


# ── main ────────────────────────────────────────────────────────────
def main():
    args = parse_args()
    if args.end_year < args.start_year:
        sys.exit("--end-year must be >= --start-year")

    tag = f"{args.model.lower()}_{args.experiment}"
    paths = {
        "out_root":  args.root / "out",
        "adt_dir":   args.root / "out" / f"adt_{tag}",
        "eddy_dir":  args.root / "out" / f"eddies_{tag}",
        "track_dir": args.root / "out" /
                     f"tracks_{tag}_{args.start_year}_{args.end_year}",
        "logs":      args.root / "logs",
        "figs":      args.root / "figs",
    }
    for p in paths.values():
        p.mkdir(parents=True, exist_ok=True)

    n_years = args.end_year - args.start_year + 1
    n_days = sum(366 if calendar.isleap(y) else 365
                 for y in range(args.start_year, args.end_year + 1))
    log.info("=" * 60)
    log.info("Multi-year %d..%d  (%d years, %d days)  model=%s exp=%s",
             args.start_year, args.end_year, n_years, n_days,
             args.model, args.experiment)
    log.info("Workers: %d detect, %d fetch streams", args.workers,
             args.fetch_parallel)
    log.info("Outputs under %s", args.root)
    log.info("=" * 60)

    wall_t0 = time.perf_counter()
    summary = []

    if not args.skip_fetch:
        t0 = time.perf_counter()
        n_fail = fetch_phase(args, paths)
        summary.append(("fetch+interp", time.perf_counter() - t0, n_fail))
        if n_fail:
            log.error("Fetch had %d failed month(s) — aborting", n_fail)
            sys.exit(2)

    if not args.skip_detect:
        t0 = time.perf_counter()
        rc = detect_phase(args, paths)
        summary.append(("detect", time.perf_counter() - t0, rc))
        if rc:
            log.error("Detect failed — aborting")
            sys.exit(3)

    if not args.skip_track:
        t0 = time.perf_counter()
        rc = track_phase(args, paths)
        summary.append(("track", time.perf_counter() - t0, rc))
        if rc:
            log.error("Track had %d failed polarity/-ies", rc)
            sys.exit(4)

    if not args.no_plots:
        t0 = time.perf_counter()
        rc = plot_phase(args, paths)
        summary.append(("plots", time.perf_counter() - t0, rc))

    wall = time.perf_counter() - wall_t0
    log.info("=" * 60)
    for name, secs, fail in summary:
        log.info("  %-14s  %s   fail=%d", name, fmt_dt(secs), fail)
    log.info("  %-14s  %s", "TOTAL", fmt_dt(wall))
    log.info("Outputs:")
    for k, v in paths.items():
        log.info("  %-10s : %s", k, v)


if __name__ == "__main__":
    main()
