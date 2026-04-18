"""End-to-end ocean-eddy pipeline for one calendar year of DestinE data.

Phases (all idempotent; reruns skip work whose outputs already exist):

  1. Fetch + interp  : month-by-month FDB pulls of avg_zos, HEALPix→0.25°
                       wet-renormalized bilinear interp, one NetCDF per day.
                       Runs in the dedl env.
  2. Detection       : parallel EddyId across every day in the year.
                       Runs in the eddy env.
  3. Tracking        : one EddyTracking call per polarity over the whole year.
                       Runs in the eddy env.
  4. Plots           : optional summary figures (--no-plots to skip).

Usage:
    python pipeline_year.py --year 2014 --model IFS-FESOM
"""

import argparse
import calendar
import logging
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
log = logging.getLogger("pipeline-year")


# ── helpers ─────────────────────────────────────────────────────────
def parse_args():
    ap = argparse.ArgumentParser()
    ap.add_argument("--year", type=int, required=True)
    ap.add_argument("--model", default="IFS-FESOM")
    ap.add_argument("--experiment", default="hist")
    ap.add_argument("--resolution", default="high")
    ap.add_argument("--root", type=Path, default=DEFAULT_ROOT,
                    help="Root for out/, logs/, figs/")
    ap.add_argument("--workers", type=int, default=8,
                    help="Parallel EddyId workers (default 8; box has 64 cores)")
    ap.add_argument("--months", type=int, nargs="+", default=None,
                    help="Restrict to these months 1..12 (default all)")
    ap.add_argument("--skip-fetch", action="store_true")
    ap.add_argument("--skip-detect", action="store_true")
    ap.add_argument("--skip-track", action="store_true")
    ap.add_argument("--no-plots", action="store_true")
    ap.add_argument("--force", "-f", action="store_true",
                    help="Recompute even if outputs exist")
    return ap.parse_args()


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
    log.info("  %s rc=%d in %.1f s (%.1f min)", label, proc.returncode, dt, dt/60)
    return proc.returncode


# ── phase 1 — fetch + interp ────────────────────────────────────────
def fetch_phase(args, paths) -> int:
    months = args.months if args.months else list(range(1, 13))
    fail = 0
    for m in months:
        last = calendar.monthrange(args.year, m)[1]
        start_iso = f"{args.year}-{m:02d}-01"
        end_iso   = f"{args.year}-{m:02d}-{last:02d}"
        rc = run_logged(
            [DEDL_PYTHON, str(HERE / "download_interp_zos.py"),
             "--start-date", start_iso, "--end-date", end_iso,
             "--model", args.model,
             "--experiment", args.experiment,
             "--resolution", args.resolution,
             "--out-dir", str(paths["out_root"]),
             *(["--force"] if args.force else [])],
            paths["logs"] / f"{args.year}_fetch_{m:02d}.log",
            label=f"fetch {args.year}-{m:02d}",
        )
        if rc != 0:
            fail += 1
            log.error("fetch FAILED for %s-%02d", args.year, m)
    return fail


# ── phase 2 — parallel detection ────────────────────────────────────
def detect_phase(args, paths) -> int:
    rc = run_logged(
        [EDDY_PYTHON, str(HERE / "detect_eddies_batch.py"),
         "--input-dir",  str(paths["adt_dir"]),
         "--output-dir", str(paths["eddy_dir"]),
         "--workers",    str(args.workers),
         *(["--force"] if args.force else [])],
        paths["logs"] / f"{args.year}_detect.log",
        label=f"detect {args.year} ({args.workers} workers)",
    )
    return 0 if rc == 0 else 1


# ── phase 3 — tracking ──────────────────────────────────────────────
TRACKING_YAML_TPL = """\
PATHS:
  FILES_PATTERN: {pattern}
  SAVE_DIR: {save_dir}

VIRTUAL_LENGTH_MAX: 0
TRACK_DURATION_MIN: 4

CLASS:
   MODULE: py_eddy_tracker.featured_tracking.area_tracker
   CLASS: AreaTracker
"""


def track_phase(args, paths) -> int:
    paths["track_dir"].mkdir(parents=True, exist_ok=True)
    fail = 0
    for sign_short, sign_long in [("anti", "Anticyclonic"), ("cyc", "Cyclonic")]:
        save_dir = paths["track_dir"] / sign_short
        save_dir.mkdir(parents=True, exist_ok=True)
        out_nc = save_dir / f"{sign_long}.nc"
        if out_nc.exists() and not args.force:
            log.info("track %s already done -> %s", sign_short, out_nc.name)
            continue
        yaml_path = paths["track_dir"] / f"track_{sign_short}.yaml"
        yaml_path.write_text(TRACKING_YAML_TPL.format(
            pattern=str(paths["eddy_dir"] / f"{sign_long}_*.nc"),
            save_dir=str(save_dir),
        ))
        rc = run_logged(
            [EDDY_TRACKING_BIN, str(yaml_path), "-v", "INFO"],
            paths["logs"] / f"{args.year}_track_{sign_short}.log",
            label=f"track {sign_long}",
        )
        if rc != 0:
            fail += 1
    return fail


# ── phase 4 — plots ─────────────────────────────────────────────────
def plot_phase(args, paths) -> int:
    rc = run_logged(
        [EDDY_PYTHON, str(HERE / "plot_year_summary.py"),
         "--label",     f"{args.model} {args.experiment} {args.year}",
         "--track-dir", str(paths["track_dir"]),
         "--fig-dir",   str(paths["figs"] / f"year_{args.year}_{args.model.lower()}_{args.experiment}")],
        paths["logs"] / f"{args.year}_plots.log",
        label="plots",
    )
    return 0 if rc == 0 else 1


# ── main ────────────────────────────────────────────────────────────
def main():
    args = parse_args()
    tag = f"{args.model.lower()}_{args.experiment}"
    paths = {
        "out_root":  args.root / "out",
        "adt_dir":   args.root / "out" / f"adt_{tag}",
        "eddy_dir":  args.root / "out" / f"eddies_{tag}",
        "track_dir": args.root / "out" / f"tracks_{tag}_{args.year}",
        "logs":      args.root / "logs",
        "figs":      args.root / "figs",
    }
    for p in paths.values():
        p.mkdir(parents=True, exist_ok=True)

    log.info("=" * 60)
    log.info("Year %d  model=%s exp=%s res=%s", args.year, args.model,
             args.experiment, args.resolution)
    log.info("Workers=%d  root=%s", args.workers, args.root)
    log.info("=" * 60)

    wall_t0 = time.perf_counter()
    summary = []

    if not args.skip_fetch:
        t0 = time.perf_counter()
        n_fail = fetch_phase(args, paths)
        summary.append(("fetch+interp", time.perf_counter() - t0, n_fail))
        if n_fail:
            log.error("Fetch phase had %d failed month(s) — aborting", n_fail)
            sys.exit(2)

    if not args.skip_detect:
        t0 = time.perf_counter()
        rc = detect_phase(args, paths)
        summary.append(("detect", time.perf_counter() - t0, rc))
        if rc:
            log.error("Detect phase failed — aborting")
            sys.exit(3)

    if not args.skip_track:
        t0 = time.perf_counter()
        rc = track_phase(args, paths)
        summary.append(("track", time.perf_counter() - t0, rc))
        if rc:
            log.error("Track phase had %d failed polarity/-ies", rc)
            sys.exit(4)

    if not args.no_plots:
        t0 = time.perf_counter()
        rc = plot_phase(args, paths)
        summary.append(("plots", time.perf_counter() - t0, rc))

    wall = time.perf_counter() - wall_t0
    log.info("=" * 60)
    for name, secs, fail in summary:
        log.info("  %-14s  %6.1f s (%5.1f min)   fail=%d", name, secs, secs/60, fail)
    log.info("  %-14s  %6.1f s (%5.1f min)", "TOTAL", wall, wall/60)
    log.info("Outputs:")
    log.info("  ADT      : %s", paths["adt_dir"])
    log.info("  Detections: %s", paths["eddy_dir"])
    log.info("  Tracks   : %s", paths["track_dir"])
    log.info("  Figs     : %s", paths["figs"])


if __name__ == "__main__":
    main()
