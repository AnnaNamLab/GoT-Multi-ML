#!/bin/bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
IRONTHRONE_DIR_DEFAULT="$SCRIPT_DIR/1_ironthrone_multi_py"
PREPARE_MUTATION_HELPER_DEFAULT="$SCRIPT_DIR/utils/prepare_got_multi_mutation_run.py"
PREPARE_BARCODE_HELPER_DEFAULT="$SCRIPT_DIR/utils/prepare_got_multi_barcode_run.py"
PROCESS_PP_SCRIPT_DEFAULT="$SCRIPT_DIR/2_ironthrone_pp/process_ironthrone_output.py"

ID=""
MUT_PARAMS=""
BARCODE_PARAMS=""
OUTDIR=""
IRONTHRONE_DIR="$IRONTHRONE_DIR_DEFAULT"
PREPARE_MUTATION_HELPER="$PREPARE_MUTATION_HELPER_DEFAULT"
PREPARE_BARCODE_HELPER="$PREPARE_BARCODE_HELPER_DEFAULT"
PROCESS_PP_SCRIPT="$PROCESS_PP_SCRIPT_DEFAULT"
RUN_PAIRED_PP=true
FORCE=false

show_help() {
    cat <<USAGE
Usage:
  $0 --id <ID> --mut-params <MUT_PARAMS> --barcode-params <BARCODE_PARAMS> --outdir <OUTDIR> \
     [--skip-paired-pp] [--force]

Description:
  GoT-Multi wrapper around standalone IronThrone-Multi.

  This script hides the FRP-specific two-pass logic required for GoT-Multi:
    1. Use the expected-sequence CSV from the MUT params file when provided
       or generate it from the MUT / gene-sequence configs when missing
    2. Run MUT / gene-sequence IronThrone-Multi normally
    3. Copy processed_input into a BARCODE result directory
    4. Rewrite the cached processed_input params CSV to use BARCODE configs
    5. Run BARCODE / probe-barcode IronThrone-Multi with --skip-input-prep
    6. Optionally run paired post-processing to create the complete paired genotype profile

Required arguments:
  --id                           Run identifier used to name the two IronThrone result folders
  --mut-params                   MUT / gene-sequence params file
  --barcode-params               BARCODE / probe-barcode params file
  --outdir                       Root output directory for this GoT-Multi run

Paired post-processing options:
  --skip-paired-pp               Skip paired post-processing and only run the two IronThrone passes

Other options:
  --ironthrone-dir               Path to standalone IronThrone-Multi directory (default: $IRONTHRONE_DIR_DEFAULT)
  --force                        Remove existing managed output subdirectories before running
  --help                         Show this help message and exit

Outputs:
  <OUTDIR>/GeneSeq_<ID>_Results/
  <OUTDIR>/ProbeBC_<ID>_Results/
  <OUTDIR>/ironthrone_out_pp.csv                  (unless --skip-paired-pp)
  <OUTDIR>/initial_genotype_counts_by_target.pdf    (unless --skip-paired-pp)
  <OUTDIR>/wrapper_input/expected_sequence.csv     (only when auto-generated)
  <OUTDIR>/wrapper_input/params_gene_seq.generated.txt
  <OUTDIR>/wrapper_input/params_probe_bc.generated.txt
USAGE
}

while [[ "$#" -gt 0 ]]; do
    case "$1" in
        --id)
            ID="$2"
            shift
            ;;
        --mut-params)
            MUT_PARAMS="$2"
            shift
            ;;
        --barcode-params)
            BARCODE_PARAMS="$2"
            shift
            ;;
        --outdir)
            OUTDIR="$2"
            shift
            ;;
        --skip-paired-pp)
            RUN_PAIRED_PP=false
            ;;
        --ironthrone-dir)
            IRONTHRONE_DIR="$2"
            shift
            ;;
        --force)
            FORCE=true
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

if [[ -z "$ID" || -z "$MUT_PARAMS" || -z "$BARCODE_PARAMS" || -z "$OUTDIR" ]]; then
    echo "Error: --id, --mut-params, --barcode-params, and --outdir are required." >&2
    show_help
    exit 1
fi

PYTHON_BIN="${PYTHON_BIN:-$(command -v python || true)}"
if [[ -z "$PYTHON_BIN" || ! -x "$PYTHON_BIN" ]]; then
    echo "Error: active python interpreter not found in PATH" >&2
    exit 1
fi

for path in "$MUT_PARAMS" "$BARCODE_PARAMS" "$IRONTHRONE_DIR/ironthrone_multi.sh" "$PREPARE_MUTATION_HELPER" "$PREPARE_BARCODE_HELPER"; do
    if [[ ! -e "$path" ]]; then
        echo "Error: required path not found: $path" >&2
        exit 1
    fi
done
if [[ "$RUN_PAIRED_PP" == true && ! -f "$PROCESS_PP_SCRIPT" ]]; then
    echo "Error: paired post-process script not found: $PROCESS_PP_SCRIPT" >&2
    exit 1
fi

OUTDIR="$(python -c 'import os,sys; print(os.path.abspath(sys.argv[1]))' "$OUTDIR")"
GENE_RUN_ID="GeneSeq_${ID}"
BARCODE_RUN_ID="ProbeBC_${ID}"
GENE_OUTDIR="$OUTDIR/${GENE_RUN_ID}_Results"
BARCODE_OUTDIR="$OUTDIR/${BARCODE_RUN_ID}_Results"
WRAPPER_INPUT_DIR="$OUTDIR/wrapper_input"
GENERATED_EXPECTED_SEQUENCE="$WRAPPER_INPUT_DIR/expected_sequence.csv"
GENERATED_MUTATION_PARAMS="$WRAPPER_INPUT_DIR/params_gene_seq.generated.txt"
GENERATED_BARCODE_PARAMS="$WRAPPER_INPUT_DIR/params_probe_bc.generated.txt"
BARCODE_PROCESSED_INPUT_DIR="$BARCODE_OUTDIR/processed_input"
RUN_METADATA="$OUTDIR/run_metadata.txt"
FINAL_PP_CSV="$OUTDIR/ironthrone_out_pp.csv"

MANAGED_PATHS=(
    "$GENE_OUTDIR"
    "$BARCODE_OUTDIR"
    "$WRAPPER_INPUT_DIR"
    "$FINAL_PP_CSV"
    "$RUN_METADATA"
)

if [[ "$FORCE" == true ]]; then
    for path in "${MANAGED_PATHS[@]}"; do
        if [[ -e "$path" ]]; then
            rm -rf "$path"
        fi
    done
fi

mkdir -p "$OUTDIR"

for path in "$GENE_OUTDIR" "$BARCODE_OUTDIR" "$WRAPPER_INPUT_DIR"; do
    if [[ -e "$path" ]]; then
        echo "Error: output path already exists: $path" >&2
        echo "Use --force to remove existing managed outputs first." >&2
        exit 1
    fi
done
if [[ -e "$FINAL_PP_CSV" || -e "$RUN_METADATA" ]]; then
    echo "Error: output file already exists in $OUTDIR. Use --force to overwrite managed outputs." >&2
    exit 1
fi

mkdir -p "$WRAPPER_INPUT_DIR"

printf 'id=%s\n' "$ID" > "$RUN_METADATA"
printf 'mut_params_input=%s\n' "$MUT_PARAMS" >> "$RUN_METADATA"
printf 'barcode_params_input=%s\n' "$BARCODE_PARAMS" >> "$RUN_METADATA"
printf 'gene_outdir=%s\n' "$GENE_OUTDIR" >> "$RUN_METADATA"
printf 'barcode_outdir=%s\n' "$BARCODE_OUTDIR" >> "$RUN_METADATA"
printf 'auto_generated_expected_sequence_output=%s\n' "$GENERATED_EXPECTED_SEQUENCE" >> "$RUN_METADATA"
printf 'generated_mutation_params=%s\n' "$GENERATED_MUTATION_PARAMS" >> "$RUN_METADATA"
printf 'generated_barcode_params=%s\n' "$GENERATED_BARCODE_PARAMS" >> "$RUN_METADATA"
printf 'python=%s\n' "$PYTHON_BIN" >> "$RUN_METADATA"

echo "Preparing MUT params and expected-sequence input (reuse if provided; otherwise generate from MUT configs only)..."
"$PYTHON_BIN" "$PREPARE_MUTATION_HELPER" \
    --mutation-params "$MUT_PARAMS" \
    --out-params "$GENERATED_MUTATION_PARAMS" \
    --out-expected-sequence-csv "$GENERATED_EXPECTED_SEQUENCE"

echo "Running GoT-Multi IronThrone wrapper: MUT pass first..."
PYTHON_BIN="$PYTHON_BIN" bash "$IRONTHRONE_DIR/ironthrone_multi.sh" \
    --id "$GENE_RUN_ID" \
    --params "$GENERATED_MUTATION_PARAMS" \
    --outdir "$GENE_OUTDIR" \
    --iron_throne_dir "$IRONTHRONE_DIR"

if [[ ! -d "$GENE_OUTDIR/processed_input" ]]; then
    echo "Error: MUT run did not create processed_input: $GENE_OUTDIR/processed_input" >&2
    exit 1
fi

echo "Preparing BARCODE pass by reusing processed_input from the MUT run..."
mkdir -p "$BARCODE_OUTDIR"
cp -a "$GENE_OUTDIR/processed_input" "$BARCODE_OUTDIR/"

"$PYTHON_BIN" "$PREPARE_BARCODE_HELPER" \
    --mutation-params "$GENERATED_MUTATION_PARAMS" \
    --barcode-params "$BARCODE_PARAMS" \
    --processed-input-dir "$BARCODE_PROCESSED_INPUT_DIR" \
    --out-params "$GENERATED_BARCODE_PARAMS"

echo "Running GoT-Multi IronThrone wrapper: BARCODE pass using cached processed_input..."
PYTHON_BIN="$PYTHON_BIN" bash "$IRONTHRONE_DIR/ironthrone_multi.sh" \
    --id "$BARCODE_RUN_ID" \
    --params "$GENERATED_BARCODE_PARAMS" \
    --outdir "$BARCODE_OUTDIR" \
    --iron_throne_dir "$IRONTHRONE_DIR" \
    --skip-input-prep

if [[ "$RUN_PAIRED_PP" == true ]]; then
    export NUMBA_DISABLE_JIT="${NUMBA_DISABLE_JIT:-1}"
    echo "Running paired post-processing to create the complete GoT-Multi genotype profile..."
    "$PYTHON_BIN" "$PROCESS_PP_SCRIPT" \
        --mutdir "$GENE_OUTDIR" \
        --bardir "$BARCODE_OUTDIR" \
        --outdir "$OUTDIR"
    echo "Paired genotyping output written to: $FINAL_PP_CSV"
else
    echo "Skipped paired post-processing."
    echo "Use the following directories as inputs for later paired post-processing:"
    echo "  MUT: $GENE_OUTDIR"
    echo "  BARCODE: $BARCODE_OUTDIR"
fi

echo "GoT-Multi IronThrone wrapper completed."
