"""Path resolution helper.

Centralises the machine-specific paths the pipeline needs, so users only
have to configure them in one place. Each path can be set via an environment
variable (preferred) or overridden on each script's CLI. Reasonable fall-
backs are tried if nothing is explicitly set.

Environment variables
---------------------
EDDY_DEDL_PYTHON    Python with z3fdb + fdb_xarray + healpy + xarray +
                    cartopy (typically the DestinE Data Lake env).
EDDY_PYTHON         Python with py-eddy-tracker installed (usually a
                    user-prefix conda env built from py-eddy-tracker's
                    own environment.yml).
EDDY_EDDYID_BIN     Path to the EddyId CLI entry-point (sibling of
                    EDDY_PYTHON).
EDDY_TRACKING_BIN   Path to the EddyTracking CLI entry-point (sibling of
                    EDDY_PYTHON).
EDDY_FDB_CONFIG     z3fdb / fdb routing config. Default ~/config.yaml.
EDDY_FDB_XARRAY     Local checkout of fdb-xarray; added to sys.path at
                    runtime. Unused if fdb_xarray is pip-installed.
EDDY_ROOT           Root for outputs (creates out/, figs/, logs/ under
                    it). Default /bigdisk/work/eddy.
"""
from __future__ import annotations

import os
import shutil
import sys
from pathlib import Path


def _first_that_exists(candidates):
    for c in candidates:
        if c and Path(c).exists():
            return Path(c)
    return None


def dedl_python() -> Path:
    """Python used for the FDB fetch / interpolation subprocesses."""
    p = os.environ.get("EDDY_DEDL_PYTHON")
    if p:
        return Path(p)
    got = _first_that_exists([
        "/opt/miniconda3/envs/dedl/bin/python",
        Path.home() / "miniconda3" / "envs" / "dedl" / "bin" / "python",
        Path.home() / "mambaforge" / "envs" / "dedl" / "bin" / "python",
    ])
    if got:
        return got
    return Path(shutil.which("python") or sys.executable)


def eddy_python() -> Path:
    """Python used for detection, tracking, plotting (py-eddy-tracker env)."""
    p = os.environ.get("EDDY_PYTHON")
    if p:
        return Path(p)
    got = _first_that_exists([
        Path.home() / "envs" / "eddy" / "bin" / "python",
        "/opt/miniconda3/envs/eddy/bin/python",
        Path.home() / "miniconda3" / "envs" / "eddy" / "bin" / "python",
    ])
    if got:
        return got
    return Path(shutil.which("python") or sys.executable)


def eddyid_bin() -> Path:
    """EddyId console script from py-eddy-tracker."""
    p = os.environ.get("EDDY_EDDYID_BIN")
    if p:
        return Path(p)
    return eddy_python().parent / "EddyId"


def eddytracking_bin() -> Path:
    """EddyTracking console script from py-eddy-tracker."""
    p = os.environ.get("EDDY_TRACKING_BIN")
    if p:
        return Path(p)
    return eddy_python().parent / "EddyTracking"


def fdb_config() -> Path:
    p = os.environ.get("EDDY_FDB_CONFIG")
    if p:
        return Path(p)
    return Path.home() / "config.yaml"


def fdb_xarray_path() -> Path | None:
    """Optional: local fdb-xarray checkout to prepend to sys.path."""
    p = os.environ.get("EDDY_FDB_XARRAY")
    if p:
        return Path(p)
    candidate = Path.home() / "fdb-xarray"
    if (candidate / "fdb_xarray").is_dir():
        return candidate
    return None


def default_root() -> Path:
    p = os.environ.get("EDDY_ROOT")
    return Path(p) if p else Path("/bigdisk/work/eddy")


def out_dir() -> Path:
    return default_root() / "out"


def fig_dir() -> Path:
    return default_root() / "figs"


def log_dir() -> Path:
    return default_root() / "logs"


def inject_fdb_xarray():
    """Call at the top of scripts that import fdb_xarray. No-op if the module
    is already importable or if the path is unset."""
    try:
        import fdb_xarray  # noqa: F401
        return
    except ImportError:
        pass
    fx = fdb_xarray_path()
    if fx and str(fx) not in sys.path:
        sys.path.insert(0, str(fx))
