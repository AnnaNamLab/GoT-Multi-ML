#!/bin/bash

set -e
# === This script runs ironthrone using input parameters csv file (v4.2 - 250116) === #

# Function to display help message
show_help() {
    echo "Usage: $0 --id <ID> --params <PARAMS_FILE> [-o|--outdir <OUTPUT_DIR>] [--iron_throne_dir <IRONTHRONE_DIR>] [--skip-input-prep]"
    echo ""
    echo "Options (Required):"
    echo "  --id                    ID for the run"
    echo "  --params                Path to the parameters file"
    echo ""
    echo "Options (Optional):"
    echo "  -o, --outdir            Output directory (default: <ID>_Results)"
    echo "  --iron_throne_dir       Directory containing IronThrone scripts (default: ./ironthrone)"
    echo "  --skip-input-prep       Skip the input preparation step"
    echo "  --help                  Show this help message and exit"
    echo ""
    echo "This script runs the IronThrone pipeline using the specified parameters file."
    echo "The parameters file should contain sections with headers and values as described below:"
    echo ""
    echo "[R1-fastqs]"
    echo "Path to R1_001.fastq FASTQ files, one per line (one per lane)"
    echo ""
    echo "[R2-fastqs]"
    echo "Path to R2_001.fastq FASTQ files, one per line (one per lane)"
    echo ""
    echo "[config-files]"
    echo "Path to configuration files, one per line (one per target)"
    echo "Filename convention: <PREFIX>_<target_name>_[mutation|barcode].config"
    echo "e.g. FRP9_GENlib_IRF8_Mut1_mutation.config"
    echo ""
    echo "[barcode-whitelist]"
    echo "Path to the barcode whitelist file (specific to protocol)"
    echo ""
    echo "[expected-sequence-csv]"
    echo "Path to the expected sequence CSV file (optional)"
    echo ""
    echo "[ironthrone-args]"
    echo "Comma-separated key-value pairs for IronThrone arguments"
    echo "e.g. prefix=FRP, suffix=.config, max-mismatch=0.2"
}

# Set default values
IRONTHRONE_DIR="$(dirname $(readlink -f $0))"
SKIP_INPUT_PREP=false

while [[ "$#" -gt 0 ]]; do
    case $1 in
        --id) ID="$2"; shift ;;
        --params) PARAMS_FILE="$2"; shift ;;
        -o|--outdir) OUTDIR="$2"; shift ;;
        --iron_throne_dir) IRONTHRONE_DIR="$2"; shift ;;
        --skip-input-prep) SKIP_INPUT_PREP=true ;;
        --help) show_help; exit 0 ;;
        *) echo "Unknown option: $1"; show_help; exit 1 ;;
    esac
    shift
done

# Check if required arguments are provided
if [ -z "$ID" ] || [ -z "$PARAMS_FILE" ]; then
    echo "Error: --id and --params are required."
    show_help
    exit 1
fi

# If OUTDIR is not set, use the $ID_result in current directory
OUTDIR=${OUTDIR:-"${ID}_Results"}
mkdir -p $OUTDIR

declare -A ironthrone_args
declare -a R1_FASTQ_FILES
declare -a R2_FASTQ_FILES
declare -a CONFIG_FILES
EXPECTED_SEQUENCE_CSV=""

while IFS= read -r line || [ -n "$line" ]; do
    if [[ $line =~ ^\[.*\]$ ]]; then
        section=${line:1:-1}
    elif [[ -n $line ]]; then
        if [[ $section == "ironthrone-args" ]]; then
            IFS=',' read -r key value <<< "$line"
            ironthrone_args[$key]=$value
        elif [[ $section == "R1-fastqs" ]]; then
            R1_FASTQ_FILES+=("$line")
        elif [[ $section == "R2-fastqs" ]]; then
            R2_FASTQ_FILES+=("$line")
        elif [[ $section == "config-files" ]]; then
            CONFIG_FILES+=("$line")
        elif [[ $section == "barcode-whitelist" ]]; then
            BARCODE_WHITELIST="$line"
        elif [[ $section == "expected-sequence-csv" ]]; then
            EXPECTED_SEQUENCE_CSV="$line"
        else
            echo "$section: $line"
        fi
    fi
done < "$PARAMS_FILE"

# Preprocess Arguments
PREFIX=${ironthrone_args[prefix]:-""}
SUFFIX=${ironthrone_args[suffix]:-".config"}
MAX_MISMATCH=${ironthrone_args[max-mismatch]:-0.2}

# IronThrone arguments with default values
UMILEN=${ironthrone_args[umilen]:-12}
BCLEN=${ironthrone_args[bclen]:-16}
MMTCH=${ironthrone_args[mmtch]:-0.2}
DUPCUT=${ironthrone_args[dupcut]:-2}
PCT_READ_THRESHOLD=${ironthrone_args[pct_read_threshold]:-0.75}
JOBS=${ironthrone_args[jobs]:-4}
THREADS_IRONTHRONE=${ironthrone_args[threads_ironthrone]:-8}

if [ "$SKIP_INPUT_PREP" = false ]; then
    # 1. Concatenates multiple FASTQ files into a single FASTQ file or copies the single FASTQ file
    mkdir -p $OUTDIR/processed_input
    
    # Check if first R1 file is gzipped to determine extension
    if [[ "${R1_FASTQ_FILES[0]}" == *.gz ]]; then
        R1_EXT=".fastq.gz"
    else
        R1_EXT=".fastq"
    fi
    if [ ${#R1_FASTQ_FILES[@]} -gt 1 ]; then
        echo "Concatenating R1 FASTQ files..."
        cat "${R1_FASTQ_FILES[@]}" > $OUTDIR/processed_input/${ID}_R1${R1_EXT}
    elif [ ${#R1_FASTQ_FILES[@]} -eq 1 ]; then
        echo "Copying single R1 FASTQ file..."
        cp "${R1_FASTQ_FILES[0]}" $OUTDIR/processed_input/${ID}_R1${R1_EXT}
    fi

    # Check if first R2 file is gzipped to determine extension
    if [[ "${R2_FASTQ_FILES[0]}" == *.gz ]]; then
        R2_EXT=".fastq.gz"
    else
        R2_EXT=".fastq"
    fi
    if [ ${#R2_FASTQ_FILES[@]} -gt 1 ]; then
        echo "Concatenating R2 FASTQ files..."
        cat "${R2_FASTQ_FILES[@]}" > $OUTDIR/processed_input/${ID}_R2${R2_EXT}
    elif [ ${#R2_FASTQ_FILES[@]} -eq 1 ]; then
        echo "Copying single R2 FASTQ file..."
        cp "${R2_FASTQ_FILES[0]}" $OUTDIR/processed_input/${ID}_R2${R2_EXT}
    fi

    echo "Processing config files to get parameter values for IronThrone..."
    bash ${IRONTHRONE_DIR}/IronThrone_input_preparation.sh \
        --fastq-r1 $OUTDIR/processed_input/${ID}_R1${R1_EXT} \
        --fastq-r2 $OUTDIR/processed_input/${ID}_R2${R2_EXT} \
        --prefix $PREFIX \
        --suffix $SUFFIX \
        -o $OUTDIR/processed_input \
        --n_subset ${ironthrone_args[n_subset]:-10} \
        -j 4 \
        --ironthrone-dir $IRONTHRONE_DIR \
        --expected-sequence-csv "$EXPECTED_SEQUENCE_CSV" \
        --max-mismatch "$MAX_MISMATCH" \
        "${CONFIG_FILES[@]}"
else
    echo "Skipping input preparation step..."
fi

# 3. Run IronThrone for each target
while IFS=',' read -r TARGET_NAME TARGET_TOTAL_N_LINES TARGET_N_READS TARGET_LINE CONFIG_FILE; do
    R1_FASTQ_FILE="$OUTDIR/processed_input/R1_R2_pairs/${TARGET_NAME}/${TARGET_NAME}_R1_001.fastq"
    R2_FASTQ_FILE="$OUTDIR/processed_input/R1_R2_pairs/${TARGET_NAME}/${TARGET_NAME}_R2_001.fastq"

    # Make sure all required variables are set
    if [[ -z "$TARGET_NAME" || -z "$TARGET_TOTAL_N_LINES" || -z "$TARGET_LINE" || -z "$CONFIG_FILE" ]]; then
        echo "Error: One or more required variables are empty."
        exit 1
    fi

    # Verify CONFIG_FILE has the SUFFIX
    if [[ "$CONFIG_FILE" != *"$SUFFIX" ]]; then
        echo "Error: CONFIG_FILE does not have the required suffix: $SUFFIX"
        exit 1
    fi

    echo ""
    # Remove "_" prefix if it exists
    SUFFIX_NO_PREFIX="${SUFFIX#_}"
    # Remove ".config" suffix if it exists
    SUFFIX_CLEAN="${SUFFIX_NO_PREFIX%.config}"
    echo "Running IronThrone for ${TARGET_NAME} ${SUFFIX_CLEAN}..."
    bash ${IRONTHRONE_DIR}/IronThroneRunner.sh \
        -f1 ${R1_FASTQ_FILE} -f2 ${R2_FASTQ_FILE} \
        -w ${BARCODE_WHITELIST} \
        -c "$CONFIG_FILE" \
        --target_lines "$TARGET_LINE" --total_lines "$TARGET_TOTAL_N_LINES" \
        --sample ${TARGET_NAME} --log ${TARGET_NAME}.log \
        --outdir ${OUTDIR}/${TARGET_NAME} --threads_ironthrone $THREADS_IRONTHRONE --jobs $JOBS \
        --umilen ${UMILEN} --bclen ${BCLEN} --mmtch ${MMTCH} \
        --dupcut ${DUPCUT} --pct_read_threshold ${PCT_READ_THRESHOLD} \
        --ironthrone_dir $IRONTHRONE_DIR --shuffle-fastq
    echo "IroneThrone run finished for $CONFIG_FILE"

    # Remove unnecessary output files
    rm -r ${OUTDIR}/${TARGET_NAME}/ironthrone_output
    rm -r ${OUTDIR}/${TARGET_NAME}/logs
    rm -r ${OUTDIR}/${TARGET_NAME}/preprocessing_fastqs
    rm -r ${OUTDIR}/${TARGET_NAME}/split_fastq

    echo ""

done < $OUTDIR/processed_input/${PREFIX}ironethrone_params.csv

if [ -d ${OUTDIR}/processed_input/target_fastqs_r2 ]; then
    rm -r ${OUTDIR}/processed_input/target_fastqs_r2
fi
