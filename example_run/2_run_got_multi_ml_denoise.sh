#!/bin/bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
PROJECT_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"
WRAPPER="$PROJECT_ROOT/run_got_multi_ml_denoise.sh"
GOT_MULTI_OUTDIR="$SCRIPT_DIR/got_multi_ml_run"
GEX="$SCRIPT_DIR/input/sc_gex.h5ad"
TARGET_INFO="$SCRIPT_DIR/input/target_info.csv"
CELL_GROUP="cell_group"
NEG_CTRL_GROUP="T-cell,Macrophage,Myo Fibro,Vascular EC"
ADDITIONAL_FEATURES="nCount_RNA,nFeature_RNA,percent.mt"
ML_OUTDIR="$GOT_MULTI_OUTDIR/ml_denoise"
PYTHON_BIN="${PYTHON_BIN:-$(command -v python || true)}"

if [ -z "$PYTHON_BIN" ] || [ ! -x "$PYTHON_BIN" ]; then
    echo "Error: active python interpreter not found in PATH" >&2
    exit 1
fi

cd "$SCRIPT_DIR"
PYTHON_BIN="$PYTHON_BIN" bash "$WRAPPER" \
    --got-multi-outdir "$GOT_MULTI_OUTDIR" \
    --gex "$GEX" \
    --target-info "$TARGET_INFO" \
    --cell-group "$CELL_GROUP" \
    --negative-control-cell-group "$NEG_CTRL_GROUP" \
    --additional-features "$ADDITIONAL_FEATURES" \
    --ml-outdir "$ML_OUTDIR"
