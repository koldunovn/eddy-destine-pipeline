"""Run EddyId in parallel over every per-day ADT NetCDF in a directory.

Each input "adt_DestinE_<model>_<exp>_<YYYYMMDD>.nc" produces
"Anticyclonic_<YYYYMMDD>T000000.nc" and "Cyclonic_<YYYYMMDD>T000000.nc"
in the chosen output directory. Idempotent: skip days where both
output NetCDFs already exist.

  $EDDY_PYTHON detect_eddies_batch.py \\
      --input-dir  $EDDY_ROOT/out/adt_ifs-fesom_hist \\
      --output-dir $EDDY_ROOT/out/eddies_ifs-fesom_hist \\
      --workers 8

Geostrophic u/v are derived inside py-eddy-tracker from SSH (passing the
literal "None" strings for the u/v args triggers grid.add_uv).
"""

import argparse
import logging
import multiprocessing as mp
import re
import shutil
import subprocess
import sys
import time
from pathlib import Path

logging.basicConfig(level=logging.INFO,
                    format="%(asctime)s %(levelname)s %(message)s")
log = logging.getLogger("detect")

import eddy_paths
EDDYID_BIN = str(eddy_paths.eddyid_bin())
DATE_RE = re.compile(r"_(\d{8})\.nc$")


def find_inputs(input_dir: Path):
    files = sorted(input_dir.glob("adt_*_????????.nc"))
    out = []
    for f in files:
        m = DATE_RE.search(f.name)
        if m:
            out.append((m.group(1), f))
    return out


def already_done(output_dir: Path, yyyymmdd: str) -> bool:
    """Output files after rename (EddyId writes _T000000.nc; we strip it)."""
    a = output_dir / f"Anticyclonic_{yyyymmdd}.nc"
    c = output_dir / f"Cyclonic_{yyyymmdd}.nc"
    return a.exists() and c.exists()


def _rename_eddyid_output(output_dir: Path, yyyymmdd: str):
    """EddyId names files <sign>_YYYYMMDDTHHMMSS.nc, but EddyTracking's
    hardcoded date_regexp can only parse <sign>_YYYYMMDD.nc. Strip the
    'T000000' time suffix so the tracking step works without a custom regex."""
    for sign in ("Anticyclonic", "Cyclonic"):
        src = output_dir / f"{sign}_{yyyymmdd}T000000.nc"
        dst = output_dir / f"{sign}_{yyyymmdd}.nc"
        if src.exists():
            src.rename(dst)


def run_one(task):
    yyyymmdd, infile, output_dir, log_dir = task
    log_file = Path(log_dir) / f"{yyyymmdd}.log"
    cmd = [EDDYID_BIN, str(infile), yyyymmdd,
           "adt", "None", "None", "longitude", "latitude",
           str(output_dir), "-v", "INFO"]
    t0 = time.perf_counter()
    with log_file.open("w") as f:
        f.write(" ".join(cmd) + "\n\n")
        f.flush()
        rc = subprocess.run(cmd, stdout=f, stderr=subprocess.STDOUT,
                             cwd=output_dir).returncode
    if rc == 0:
        _rename_eddyid_output(Path(output_dir), yyyymmdd)
    dt = time.perf_counter() - t0
    return yyyymmdd, rc, dt


def parse_args():
    ap = argparse.ArgumentParser()
    ap.add_argument("--input-dir", type=Path, required=True)
    ap.add_argument("--output-dir", type=Path, required=True)
    ap.add_argument("--workers", type=int, default=4)
    ap.add_argument("--force", "-f", action="store_true")
    return ap.parse_args()


def main():
    args = parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    log_dir = args.output_dir / "logs"
    log_dir.mkdir(exist_ok=True)

    inputs = find_inputs(args.input_dir)
    if not inputs:
        log.error("No inputs found in %s matching adt_*_YYYYMMDD.nc", args.input_dir)
        sys.exit(2)

    pending = [(d, f) for d, f in inputs if args.force or not already_done(args.output_dir, d)]
    log.info("%d inputs, %d need detection (%d already done)",
             len(inputs), len(pending), len(inputs) - len(pending))
    if not pending:
        return

    tasks = [(d, f, args.output_dir, log_dir) for d, f in pending]
    t0 = time.perf_counter()
    fail = 0
    with mp.get_context("spawn").Pool(args.workers) as pool:
        for d, rc, dt in pool.imap_unordered(run_one, tasks):
            ok = "ok " if rc == 0 else "FAIL"
            if rc != 0:
                fail += 1
            log.info("  [%s] %s  %.1fs  log=%s/%s.log",
                     d, ok, dt, log_dir.name, d)
    wall = time.perf_counter() - t0
    log.info("Detection wall: %.1f s for %d days  (avg %.1f s/day, %d failed)",
             wall, len(pending), wall / len(pending), fail)
    if fail:
        sys.exit(3)


if __name__ == "__main__":
    main()
