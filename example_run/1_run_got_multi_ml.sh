#!/bin/bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
PROJECT_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"
WRAPPER="$PROJECT_ROOT/run_got_multi_ml.sh"
OUTDIR="$SCRIPT_DIR/got_multi_ml_run"
PYTHON_BIN="${PYTHON_BIN:-$(command -v python || true)}"

if [ -z "$PYTHON_BIN" ] || [ ! -x "$PYTHON_BIN" ]; then
    echo "Error: active python interpreter not found in PATH" >&2
    exit 1
fi

cd "$SCRIPT_DIR"
PYTHON_BIN="$PYTHON_BIN" bash "$WRAPPER" \
    --id FFPE6_subset100k_got_multi_ml_example \
    --mut-params "$SCRIPT_DIR/input/params_gene_seq.txt" \
    --barcode-params "$SCRIPT_DIR/input/params_probe_bc.txt" \
    --outdir "$OUTDIR"
