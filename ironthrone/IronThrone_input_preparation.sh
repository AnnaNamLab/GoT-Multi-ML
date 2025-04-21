#!/bin/bash

set -e

# === This script prepares input files for the Ironethrone pipeline (v4.1 - 250109) === #
# Input: Config files, FASTQ file
# Output: CSV file containing target name, number of reads, target line, and concatenated sequence
# Requirements: seqkit

# Function to display help message
show_help() {
    echo "Usage: $0 -f <FASTQ_R2> -r <FASTQ_R1> -o <OUTPUT_DIRECTORY> [-p <PREFIX>] [-j <NUM_JOBS>] [--help] <CONFIG_FILES...>"
    echo ""
    echo "Options:"
    echo "  -f, --fastq-r2          Input FASTQ R2 file"
    echo "  -r, --fastq-r1          Input FASTQ R1 file"
    echo "  -o, --output-directory  Output directory for CSV and FASTQ files (DEFAULT: current directory)"
    echo "  -i, --ironthrone-dir    Directory containing IronThrone scripts"
    echo "  -s, --suffix            Suffix of the config files (DEFAULT: .config)"
    echo "  -p, --prefix            Prefix of the config files (DEFAULT: No prefix)"
    echo "                          Config file name is expected to be in the format [PREFIX]<TARGET_NAME><SUFFIX> (e.g. FRPtest_PrimaryCells_DUSP2_Mut1_mutation.config)"
    echo "                          (e.g. FRPtest_PrimaryCells_DUSP2_Mut1_mutation.config -> FRPtest_PrimaryCells_)"
    echo "  -n, --n_subset          Number to subset the fastq in the ironthrone run (DEFAULT: 10)"
    echo "  -j, --jobs              Number of parallel jobs (DEFAULT: 2)"
    echo "  --expected-sequence-csv Path to the expected sequence CSV file (optional)"
    echo "  --max-mismatch          Maximum mismatch percentage (optional, default: 0.2)"
    echo "  --help                  Show this help message and exit"
    echo ""
    echo "This script processes config files specified as arguments."
    echo "For each config file, it extracts sequences from the first three lines,"
    echo "concatenates them, and searches for the concatenated sequence in the provided FASTQ file."
    echo "The script calculates the number of reads and determines a target line value."
    echo "The results are written to the specified output CSV file (ironethrone_params.csv) in the format:"
    echo "TARGET_NAME,TARGET_TOTAL_N_LINES,TARGET_N_READS,TARGET_LINE,CONCATENATED_SEQUENCE"
    echo ""
    echo "Example input file format:"
    echo "GCATCCACTTCCTGATTATAATCTT  1  25"
    echo "CTTGCCAGCGTGTTTCCAAGGGAT   27 50"
    echo "G                          26 26"
    echo "T                          26 26"
}

# Default values
OUTPUT_DIR="."
NUM_JOBS=2
SUFFIX=".config"
PREFIX=""
N_SUBSET=10
IRONTHRONE_DIR="$(dirname $(readlink -f $0))"
EXPECTED_SEQUENCE_CSV=""
MAX_MISMATCH=0.2

# Parse command line arguments
while [[ "$#" -gt 0 ]]; do
    case $1 in
        -r|--fastq-r1) FASTQ_R1="$2"; shift 2 ;;
        -f|--fastq-r2) FASTQ_R2="$2"; shift 2 ;;
        -o|--output-directory) OUTPUT_DIR="$2"; shift 2 ;;
        -s|--suffix) SUFFIX="$2"; shift 2 ;;
        -p|--prefix) PREFIX="$2"; shift 2 ;;
        -n|--n_subset) N_SUBSET="$2"; shift 2 ;;
        -j|--jobs) NUM_JOBS="$2"; shift 2 ;;
        -i|--ironthrone-dir) IRONTHRONE_DIR="$2"; shift 2 ;;
        --expected-sequence-csv) EXPECTED_SEQUENCE_CSV="$2"; shift 2 ;;
        --max-mismatch) MAX_MISMATCH="$2"; shift 2 ;;
        --help) show_help; exit 0 ;;
        --) shift; CONFIG_FILES+=("$@"); break ;;
        -*) echo "Unknown option: $1"; show_help; exit 1 ;;
        *) CONFIG_FILES+=("$1"); shift ;;
    esac
done

# Check if required arguments are provided and validate files
if [ -z "$FASTQ_R2" ] || [ -z "$FASTQ_R1" ] || [ ${#CONFIG_FILES[@]} -eq 0 ]; then
    echo "Error: FASTQ R2 file, FASTQ R1 file, and at least one config file are required."
    show_help
    exit 1
fi

# Check if config files exist
for config_file in "${CONFIG_FILES[@]}"; do
    if [ ! -f "$config_file" ]; then
        echo "Error: Config file '$config_file' does not exist."
        exit 1
    fi
done

# Check if FASTQ files exist and have matching compression state
if [ ! -f "$FASTQ_R2" ] || [ ! -f "$FASTQ_R1" ]; then
    echo "Error: FASTQ files do not exist. (FASTQ files should have '_001.fastq/.fastq.gz' suffix)"
    exit 1
fi

# Check if both files have the same compression state
is_r1_compressed=$([[ "$FASTQ_R1" =~ \.gz$ ]] && echo "true" || echo "false")
is_r2_compressed=$([[ "$FASTQ_R2" =~ \.gz$ ]] && echo "true" || echo "false")

if [ "$is_r1_compressed" != "$is_r2_compressed" ]; then
    echo "Error: FASTQ files must have the same compression state (both .gz or both uncompressed)"
    exit 1
fi

# Single variable for compression state
is_compressed=$is_r1_compressed
FASTQ_EXT=$($is_compressed && echo ".fastq.gz" || echo ".fastq")

# Export variables needed by parallel processes
export FASTQ_EXT
export is_compressed
export N_SUBSET

OUTPUT_FILE="${PREFIX}ironethrone_params.csv"
mkdir -p "$OUTPUT_DIR/target_fastqs_r2"

# Function to process each file
process_file() {
    local file="$1"
    local fastq_r2="$2"
    local fastq_r1="$3"
    local output_file="$4"
    local output_dir="$5"
    local prefix="$6"
    local suffix="$7"
    local ironthrone_dir="$8"
    local expected_sequence_csv="$9"
    local max_mismatch="${10}"

    echo ""
    echo "Processing $(basename $file)..."
    # Extract the target name from the file name
    TARGET_NAME=$(basename "$file" "$suffix")
    
    # Exclude the PREFIX part if PREFIX is provided
    if [ -n "$prefix" ]; then
        TARGET_NAME=${TARGET_NAME#$prefix}
    fi

    if [ -n "$expected_sequence_csv" ]; then
        PATTERN_SEQUENCE=$(awk -F, -v target="$TARGET_NAME" '$1 == target {print $2}' "$expected_sequence_csv")
        if [ -z "$PATTERN_SEQUENCE" ]; then
            echo "TARGET $TARGET_NAME not found in the expected sequence CSV file. Using LINE1+LINE3+LINE2 from config file."
            LINE1=$(sed -n '1p' "$file" | cut -f1)
            LINE2=$(sed -n '2p' "$file" | cut -f1)
            LINE3=$(sed -n '3p' "$file" | cut -f1)
            PATTERN_SEQUENCE="${LINE1}${LINE3}${LINE2}"
        else
            echo "Using pattern sequence from expected sequence CSV file for TARGET $TARGET_NAME."
        fi
    else
        echo "Expected sequence CSV file not specified. Using LINE1+LINE3+LINE2 from config file."
        LINE1=$(sed -n '1p' "$file" | cut -f1)
        LINE2=$(sed -n '2p' "$file" | cut -f1)
        LINE3=$(sed -n '3p' "$file" | cut -f1)
        PATTERN_SEQUENCE="${LINE1}${LINE3}${LINE2}"
    fi

    MAX_MISMATCH_COUNT=$(echo "scale=0; ${#PATTERN_SEQUENCE} * $max_mismatch / 1" | bc)
    echo "Using MAX_MISMATCH_COUNT=$MAX_MISMATCH_COUNT to search and subset fastq for TARGET=$TARGET_NAME"

    # Remove existing output files if they exist
    if [ -f "${output_dir}/target_fastqs_r2/${TARGET_NAME}_R2_001.fastq.gz" ]; then
        rm "${output_dir}/target_fastqs_r2/${TARGET_NAME}_R2_001.fastq.gz"
    fi
    if [ -d "${output_dir}/R1_R2_pairs/${TARGET_NAME}" ]; then
        rm -r "${output_dir}/R1_R2_pairs/${TARGET_NAME}"
    fi
    
    echo "Subsetting for the target sequence in R2 FASTQ file..."
    seqkit grep --by-seq --max-mismatch "$MAX_MISMATCH_COUNT" \
        --pattern "$PATTERN_SEQUENCE" \
        "$fastq_r2" \
        -o "${output_dir}/target_fastqs_r2/${TARGET_NAME}_R2_001${FASTQ_EXT}" \
        -j 24

    # Update line counting with single compression state
    if $is_compressed; then
        TARGET_SUBSET_N_LINES=$(zcat "${output_dir}/target_fastqs_r2/${TARGET_NAME}_R2_001${FASTQ_EXT}" | wc -l)
    else
        TARGET_SUBSET_N_LINES=$(wc -l < "${output_dir}/target_fastqs_r2/${TARGET_NAME}_R2_001${FASTQ_EXT}")
    fi

    # Calculate the number of reads
    TARGET_N_READS=$((TARGET_SUBSET_N_LINES / 4))

    # Calculate TARGET_LINE as TARGET_SUBSET_N_LINES divided by n_subset and ensure it is a multiple of 4
    ## ? Adding 3 ensures that when we divide by 4 and then multiply by 4,
    ## ? we round up to the next multiple of 4 if TARGET_N_READS is not already a multiple of 4.
    TARGET_LINE=$(( ( (TARGET_SUBSET_N_LINES / $N_SUBSET + 3) / 4 ) * 4 ))

    # Append the result to the output file
    echo "${TARGET_NAME},${TARGET_SUBSET_N_LINES},${TARGET_N_READS},${TARGET_LINE},${file}" >> "${output_dir}/${output_file}"
}

export -f process_file

# Initialize the output file (no header)
> "${OUTPUT_DIR}/${OUTPUT_FILE}"

# Run the parallel processing for subsetting and counting
parallel -j "$NUM_JOBS" process_file ::: "${CONFIG_FILES[@]}" ::: "$FASTQ_R2" ::: "$FASTQ_R1" ::: "$OUTPUT_FILE" ::: "$OUTPUT_DIR" ::: "$PREFIX" ::: "$SUFFIX" ::: "$IRONTHRONE_DIR" ::: "$EXPECTED_SEQUENCE_CSV" ::: "$MAX_MISMATCH"

echo "Output written to ${OUTPUT_DIR}/${OUTPUT_FILE}"

# Process R1-R2 pairs sequentially (due to memory consumption)
echo ""
while IFS=, read -r TARGET_NAME TARGET_SUBSET_N_LINES TARGET_N_READS TARGET_LINE CONFIG_FILE; do
    echo "Pairing R1-R2 files for ${TARGET_NAME}..."
    mkdir -p "$OUTPUT_DIR/R1_R2_pairs/${TARGET_NAME}"
    
    seqkit pair -1 "$FASTQ_R1" -2 "${OUTPUT_DIR}/target_fastqs_r2/${TARGET_NAME}_R2_001${FASTQ_EXT}" \
        --out-dir "$OUTPUT_DIR/R1_R2_pairs/${TARGET_NAME}" \
        -j 16

    mv "${OUTPUT_DIR}/R1_R2_pairs/${TARGET_NAME}/$(basename "$FASTQ_R1")" \
       "${OUTPUT_DIR}/R1_R2_pairs/${TARGET_NAME}/${TARGET_NAME}_R1_001${FASTQ_EXT}"
    
    # Decompress files if they were compressed
    if $is_compressed; then
        echo "Decompressing output files for ${TARGET_NAME}..."
        zcat "${OUTPUT_DIR}/R1_R2_pairs/${TARGET_NAME}/${TARGET_NAME}_R1_001${FASTQ_EXT}" > \
            "${OUTPUT_DIR}/R1_R2_pairs/${TARGET_NAME}/${TARGET_NAME}_R1_001.fastq"
        zcat "${OUTPUT_DIR}/R1_R2_pairs/${TARGET_NAME}/${TARGET_NAME}_R2_001${FASTQ_EXT}" > \
            "${OUTPUT_DIR}/R1_R2_pairs/${TARGET_NAME}/${TARGET_NAME}_R2_001.fastq"
        # Remove compressed files
        rm "${OUTPUT_DIR}/R1_R2_pairs/${TARGET_NAME}/${TARGET_NAME}_R1_001${FASTQ_EXT}"
        rm "${OUTPUT_DIR}/R1_R2_pairs/${TARGET_NAME}/${TARGET_NAME}_R2_001${FASTQ_EXT}"
    fi
done < "${OUTPUT_DIR}/${OUTPUT_FILE}"

# Sort the output file by TARGET_SUBSET_N_LINES (2nd column) numerically
sort -t',' -k2,2n -o "${OUTPUT_DIR}/${OUTPUT_FILE}" "${OUTPUT_DIR}/${OUTPUT_FILE}"

echo "All processing complete"
