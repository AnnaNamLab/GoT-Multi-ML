#!/bin/bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
IRONTHRONE_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"
RUNNER="$IRONTHRONE_DIR/ironthrone_multi.sh"
OUTDIR="$SCRIPT_DIR/example_run_output"
PYTHON_BIN="${PYTHON_BIN:-$(command -v python || true)}"

if [ -z "$PYTHON_BIN" ] || [ ! -x "$PYTHON_BIN" ]; then
    echo "Error: active python interpreter not found in PATH" >&2
    exit 1
fi

cd "$SCRIPT_DIR"
PYTHON_BIN="$PYTHON_BIN" bash "$RUNNER" \
    --id FFPE6_subset100k_ironthrone_multi_example \
    --params "$SCRIPT_DIR/input/params.txt" \
    --outdir "$OUTDIR"
