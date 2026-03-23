#!/bin/bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
PREPARE_ML_INPUT_SCRIPT_DEFAULT="$SCRIPT_DIR/3_ironthrone_denoise/prepare_ml_denoise_input.py"
ML_DENOISE_SCRIPT_DEFAULT="$SCRIPT_DIR/3_ironthrone_denoise/ml_denoise.py"

GOT_MULTI_OUTDIR=""
GEX=""
TARGET_INFO=""
CELL_GROUP="cell_group"
NEG_CTRL_GROUP=""
ADDITIONAL_FEATURES=""
ML_OUTDIR=""
PREPARE_ML_INPUT_SCRIPT="$PREPARE_ML_INPUT_SCRIPT_DEFAULT"
ML_DENOISE_SCRIPT="$ML_DENOISE_SCRIPT_DEFAULT"
MAX_ALLOWED_FPR=""
ALPHA=""
RANDOM_SEED=""

show_help() {
    cat <<USAGE
Usage:
  $0 --got-multi-outdir <OUTDIR> --gex <GEX_H5AD> --target-info <TARGET_INFO_CSV> \
     [--cell-group <GEX_OBS_COLUMN>] [--negative-control-cell-group <GROUP1,GROUP2,...>] \
     [--additional-features <GEX_OBS_COLUMNS>] [--ml-outdir <ML_OUTDIR>] \
     [--max-allowed-fpr <FLOAT>] [--alpha <FLOAT>] [--random-seed <INT>]

Description:
  Optional second-stage GoT-Multi-ML denoising wrapper.

  This script:
    1. Reads <got-multi-outdir>/ironthrone_out_pp.csv from run_got_multi_ml.sh
    2. Reads hard-coded GEX metadata columns from .obs: BC and sample
    3. Uses the configured --cell-group column from GEX .obs and standardizes it to cell_group
    4. Uses GEX .obs["experiment"] and .obs["is_non_mut"] if present
    5. Includes any requested --additional-features from GEX .obs
    6. Uses target_info to map expected_sample and includes all remaining target_info columns as ML features
    7. Rebuilds expected, Y, and Y_num cleanly even if the input CSV is an older pre-split artifact
    8. Runs ml_denoise.py on the prepared ML input

Notes:
  - target_info must contain 'target' and either 'expected_sample' or 'sample'
  - a target may map to multiple expected samples, either across repeated rows or as a comma-separated list in one row
  - mutant calls in any expected sample are labeled MUT; mutant calls in non-expected samples or negative-control cell groups are labeled FalsePositive
  - if GEX .obs lacks 'is_non_mut', provide --negative-control-cell-group so it can be derived from the configured cell-group column

Outputs:
  <ML_OUTDIR>/ironthrone_out_pp_ml_input.csv
  <ML_OUTDIR>/got_multi_out_refined.csv
USAGE
}

while [[ "$#" -gt 0 ]]; do
    case "$1" in
        --got-multi-outdir)
            GOT_MULTI_OUTDIR="$2"
            shift
            ;;
        --gex)
            GEX="$2"
            shift
            ;;
        --target-info)
            TARGET_INFO="$2"
            shift
            ;;
        --cell-group)
            CELL_GROUP="$2"
            shift
            ;;
        --negative-control-cell-group)
            NEG_CTRL_GROUP="$2"
            shift
            ;;
        --additional-features)
            ADDITIONAL_FEATURES="$2"
            shift
            ;;
        --ml-outdir)
            ML_OUTDIR="$2"
            shift
            ;;
        --prepare-ml-input-script)
            PREPARE_ML_INPUT_SCRIPT="$2"
            shift
            ;;
        --ml-denoise-script)
            ML_DENOISE_SCRIPT="$2"
            shift
            ;;
        --max-allowed-fpr)
            MAX_ALLOWED_FPR="$2"
            shift
            ;;
        --alpha)
            ALPHA="$2"
            shift
            ;;
        --random-seed)
            RANDOM_SEED="$2"
            shift
            ;;
        --help)
            show_help
            exit 0
            ;;
        *)
            echo "Unknown option: $1" >&2
            show_help
            exit 1
            ;;
    esac
    shift
done

if [[ -z "$GOT_MULTI_OUTDIR" || -z "$GEX" || -z "$TARGET_INFO" ]]; then
    echo "Error: --got-multi-outdir, --gex, and --target-info are required." >&2
    show_help
    exit 1
fi

PYTHON_BIN="${PYTHON_BIN:-$(command -v python || true)}"
if [[ -z "$PYTHON_BIN" || ! -x "$PYTHON_BIN" ]]; then
    echo "Error: active python interpreter not found in PATH" >&2
    exit 1
fi

GOT_MULTI_OUTDIR="$(python -c 'import os,sys; print(os.path.abspath(sys.argv[1]))' "$GOT_MULTI_OUTDIR")"
INPUT_PP="$GOT_MULTI_OUTDIR/ironthrone_out_pp.csv"
ML_OUTDIR="${ML_OUTDIR:-$GOT_MULTI_OUTDIR/ml_denoise}"
ML_INPUT_CSV="$ML_OUTDIR/ironthrone_out_pp_ml_input.csv"

for path in "$INPUT_PP" "$GEX" "$TARGET_INFO" "$PREPARE_ML_INPUT_SCRIPT" "$ML_DENOISE_SCRIPT"; do
    if [[ ! -e "$path" ]]; then
        echo "Error: required path not found: $path" >&2
        exit 1
    fi
done

mkdir -p "$ML_OUTDIR"

echo "Preparing ML denoising input from paired GoT-Multi genotype output..."
PREP_CMD=(
    "$PYTHON_BIN" "$PREPARE_ML_INPUT_SCRIPT"
    --ironthrone-pp "$INPUT_PP"
    --gex "$GEX"
    --target-info "$TARGET_INFO"
    --cell-group "$CELL_GROUP"
    --outdir "$ML_OUTDIR"
)
if [[ -n "$NEG_CTRL_GROUP" ]]; then
    PREP_CMD+=(--negative-control-cell-group "$NEG_CTRL_GROUP")
fi
if [[ -n "$ADDITIONAL_FEATURES" ]]; then
    PREP_CMD+=(--additional-features "$ADDITIONAL_FEATURES")
fi
"${PREP_CMD[@]}"

echo "Running optional ML denoising..."
ML_CMD=(
    "$PYTHON_BIN" "$ML_DENOISE_SCRIPT"
    --input-features "$ML_INPUT_CSV"
    --outdir "$ML_OUTDIR"
)
if [[ -n "$MAX_ALLOWED_FPR" ]]; then
    ML_CMD+=(--max-allowed-fpr "$MAX_ALLOWED_FPR")
fi
if [[ -n "$ALPHA" ]]; then
    ML_CMD+=(--alpha "$ALPHA")
fi
if [[ -n "$RANDOM_SEED" ]]; then
    ML_CMD+=(--random-seed "$RANDOM_SEED")
fi
"${ML_CMD[@]}"

echo "ML denoising completed. Final outputs are under: $ML_OUTDIR"
