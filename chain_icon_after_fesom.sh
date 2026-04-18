#!/bin/bash
# Wait for the running IFS-FESOM eddy pipeline (PID in $FESOM_PID_FILE) to
# finish, then launch the same multi-year pipeline for ICON. Designed to be
# launched detached:
#
#   setsid nohup bash chain_icon_after_fesom.sh \
#       > $EDDY_ROOT/logs/chain_icon.log 2>&1 < /dev/null &
#   disown
#
# Idempotent end-to-end: pipeline_multiyear.py skips finished phases on rerun.

set -u

# Paths resolve from $EDDY_ROOT (default /bigdisk/work/eddy) and
# $EDDY_DEDL_PYTHON (default /opt/miniconda3/envs/dedl/bin/python).
# $REPO_DIR should point at this clone of eddy-destine-pipeline.
: "${EDDY_ROOT:=/bigdisk/work/eddy}"
: "${EDDY_DEDL_PYTHON:=/opt/miniconda3/envs/dedl/bin/python}"
: "${REPO_DIR:=$(dirname "$(readlink -f "$0")")}"

FESOM_PID_FILE=$EDDY_ROOT/logs/multiyear.pid
ICON_PID_FILE=$EDDY_ROOT/logs/multiyear_icon.pid
ICON_LOG=$EDDY_ROOT/logs/multiyear_1990_2014_icon_main.log

PYTHON=$EDDY_DEDL_PYTHON
PIPELINE=$REPO_DIR/pipeline_multiyear.py

ts() { date '+%Y-%m-%d %H:%M:%S'; }

FESOM_PID=$(cat "$FESOM_PID_FILE" 2>/dev/null || echo "")
if [ -z "$FESOM_PID" ]; then
    echo "$(ts)  No FESOM PID at $FESOM_PID_FILE — launching ICON immediately."
else
    echo "$(ts)  Waiting for FESOM pipeline (PID $FESOM_PID) to exit..."
    while kill -0 "$FESOM_PID" 2>/dev/null; do
        sleep 60
    done
    echo "$(ts)  FESOM pipeline $FESOM_PID exited. Starting ICON."
fi

cd "$EDDY_ROOT"

# Launch ICON pipeline in the foreground (we are already detached). Capture
# its PID for the user to monitor / kill.
"$PYTHON" "$PIPELINE" \
    --start-year 1990 --end-year 2014 \
    --model ICON --experiment hist \
    --workers 32 --fetch-parallel 8 \
    > "$ICON_LOG" 2>&1 &
ICON_PID=$!
echo "$ICON_PID" > "$ICON_PID_FILE"
echo "$(ts)  Launched ICON pipeline as PID $ICON_PID, log $ICON_LOG"

wait "$ICON_PID"
RC=$?
echo "$(ts)  ICON pipeline exited rc=$RC"
exit $RC
