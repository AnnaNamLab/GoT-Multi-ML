#!/bin/bash

# IronThrone Runner Script v2.0
# This script processes paired-end fastq files through the IronThrone pipeline
# with shuffling (optional) and automatic splitting for parallel processing.
#
# Output:
#   - Final merged output: <TARGET>.summTable.concat.umi_collapsed.txt in the output directory
#   - All intermediate files are stored in the output directory
#
# Example usage:
#   ./IronThroneRunner.sh -f1 R1.fastq -f2 R2.fastq -w whitelist.txt -c config.config -s sample1 -o output_dir
#
# Requirements:
#   - Rscript
#   - IronThrone-GoT Perl script

set -e

# Function to display usage information
usage() {
    cat << EOF
Usage: $(basename $0) [options]

Required:
    -f1, --fastqR1              Input R1 fastq file
    -f2, --fastqR2              Input R2 fastq file
    -w,  --whitelist            Barcode whitelist file
    -c,  --config               Config file

Optional:
    -tl, --target_lines         Lines per split (default: 500000)
    -n,  --total_lines          Total number of lines in fastq (optional)
    -s,  --sample               Sample name (default: myGoT)
    -o,  --outdir               Output directory (default: myGot_out/)
    --shuffle-fastq             Shuffle fastq files before splitting (flag)
    -u,  --umilen               UMI length (default: 12)
    -b,  --bclen                Barcode length (default: 16)
    -r,  --run                  Run type (default: linear)
    -m,  --mmtch                Match threshold (default: 0.2)
    -p,  --postP                Post-processing threshold (default: 0.99)
    -d,  --dupcut               Duplication cutoff (default: 2)
    -l,  --log                  Log file (default: myGoT.log)
    -k,  --keepouts             Keep intermediate files (default: 0)
    -v,  --verbose              Verbose output (default: 0)
    --pcr_read_threshold        PCR read threshold (default: 0.75)
    --levenshtein_distance      Levenshtein distance (default: 0.2)
    -t,  --threads_ironthrone   Number of threads for one iron throne run (default: 4)
    -j,  --jobs                 Number of parallel iron throne runs (default: 4)
    -i,  --ironthrone_dir      Directory containing IronThrone executable files (default: ./ironthrone)
EOF
    exit 1
}

# Function to validate inputs
validate_inputs() {
    local errors=0
    
    [[ ! -f "$fastqR1" ]] && echo "Error: R1 fastq file not found" && errors=$((errors+1))
    [[ ! -f "$fastqR2" ]] && echo "Error: R2 fastq file not found" && errors=$((errors+1))
    [[ ! -f "$whitelist" ]] && echo "Error: Whitelist file not found" && errors=$((errors+1))
	[[ ! -f "$config" ]] && echo "Error: Config file not found" && errors=$((errors+1))
    
    if ((errors > 0)); then
        exit 1
    fi
}

# Function to shuffle fastq files
shuffle_fastq_files() {
    echo "Joining fastq files..."
    paste "$fastqR1" "$fastqR2" > "$tmp_dir/combined.fastq"
    
    echo "Shuffling fastq files..."
    if (($(grep ";" "$tmp_dir/combined.fastq" | wc -l) == 0)); then
        awk '{printf("%s%s",$0,(NR%4==0)?"\n":";")}' "$tmp_dir/combined.fastq" | \
            shuf | tr ";" "\n" > "$tmp_dir/combined_shuffled.fastq"
    elif (($(grep "|" "$tmp_dir/combined.fastq" | wc -l) == 0)); then
        awk '{printf("%s%s",$0,(NR%4==0)?"\n":"|")}' "$tmp_dir/combined.fastq" | \
            shuf | tr "|" "\n" > "$tmp_dir/combined_shuffled.fastq"
    else
        echo "Error: Cannot find suitable separator for shuffling"
        exit 1
    fi
    
    echo "Separating shuffled files..."
    cut -f1 -d$'\t' "$tmp_dir/combined_shuffled.fastq" > "$tmp_dir/shuffled.R1.fastq"
    cut -f2 -d$'\t' "$tmp_dir/combined_shuffled.fastq" > "$tmp_dir/shuffled.R2.fastq"
    
    # Update input files to use shuffled versions
    fastqR1="$tmp_dir/shuffled.R1.fastq"
    fastqR2="$tmp_dir/shuffled.R2.fastq"
}

# Function to split fastq files
split_fastq_files() {
    # Use provided total_lines if available, otherwise calculate
    local actual_total_lines
    if [ -n "$total_lines" ]; then
        actual_total_lines=$total_lines
    else
        actual_total_lines=$(wc -l < "$fastqR1")
    fi
    
    local total_reads=$((actual_total_lines / 4))
    local reads_per_split=$((target_lines / 4))
    local remainder_reads=$((total_reads % reads_per_split))
    
    # Adjust split size if remainder is too small
    if ((remainder_reads * 10 < reads_per_split * 9 && remainder_reads != 0)); then
        local total_splits=$((total_reads / reads_per_split))
        local adjustment=$(( (remainder_reads / total_splits) + 1))
        reads_per_split=$((reads_per_split + adjustment))
        target_lines=$((reads_per_split * 4))
    fi
    
    # Split files
    split -d -a 4 -l "$target_lines" "$fastqR1" "$split_dir/split.R1"
    split -d -a 4 -l "$target_lines" "$fastqR2" "$split_dir/split.R2"
    
    # Rename split files
    for f in "$split_dir"/split.R*; do
        mv "$f" "${f}.fastq"
    done
    
    echo "Split into $(ls $split_dir/split.R1*.fastq | wc -l) pieces"
}

# Function to run IronThrone in parallel
run_ironthrone_parallel() {
    local commands_file="${commands_dir}/parallel_commands.txt"
    > "$commands_file"
    
    for i in "$split_dir"/split.R1*.fastq; do
        local suffix=$(basename "$i" | sed 's/split.R1\(.*\).fastq/\1/')
        local r1="$split_dir/split.R1${suffix}.fastq"
        local r2="$split_dir/split.R2${suffix}.fastq"
        local out_subdir="$result_dir/${suffix}"
        mkdir -p "$out_subdir"
        
        echo "perl ${ironthrone_dir}/IronThrone-GoT -r ${run} -f1 ${r1} -f2 ${r2} -c ${config} \
        -w ${whitelist} -u ${umilen} -b ${bclen} -o ${out_subdir} \
        -m ${mmtch} -p ${postP} -d 1 -s ${sample} -l ${log} \
        -k ${keepouts} -v ${verbose} -t ${threads}" >> "$commands_file"
    done
    
    echo "Running $jobs instances of IronThrone in parallel..."
    parallel -j "$jobs" :::: "$commands_file"
}

#Set Up Command Line Option Defaults
target_lines=500000
umilen=12
bclen=16
run=linear
mmtch=0.2
postP=0.99
dupcut=2
sample=myGoT
outdir=myGot_out
log=myGoT.log
keepouts=0
verbose=0
skip_shuf=0
pcr_read_threshold=0.75
levenshtein_distance=0.2
threads=4
jobs=4
ironthrone_dir=./ironthrone
total_lines=""

# Add function to setup directory structure
setup_directories() {
    # Remove existing output directory if it exists
    if [ -d "$outdir" ]; then
        echo "Removing existing output directory: $outdir"
        rm -rf "$outdir"
    fi
    
    # Create main output directory
    mkdir -p "$outdir"
    
    # Define subdirectories relative to outdir
    tmp_dir="${outdir}/preprocessing_fastqs"
    split_dir="${outdir}/split_fastq"
    result_dir="${outdir}/ironthrone_output"
    commands_dir="${outdir}/logs"
    
    # Create subdirectories
    mkdir -p "$tmp_dir" "$split_dir" "$result_dir" "$commands_dir"
}

# Parse command line arguments
while [ "$1" != "" ]; do
    case $1 in
        -f1 | --fastqR1 )        shift
                    fastqR1=$1
                    ;;
        -f2 | --fastqR2 )        shift
                    fastqR2=$1
                    ;;
        -tl | --target_lines )    shift
                    target_lines=$1
                    ;;
        -c | --config )        shift
                    config=$1
                    ;;
        -w | --whitelist )    shift
                    whitelist=$1
                    ;;
        -u | --umilen )        shift
                    umilen=$1
                    ;;
        -b | --bclen )        shift
                    bclen=$1
                    ;;
        -r | --run )        shift
                    run=$1
                    ;;
        -m | --mmtch )        shift
                    mmtch=$1
                    ;;
        -p | --postP )        shift
                    postP=$1
                    ;;
        -d | --dupcut )        shift
                    dupcut=$1
                    ;;
        -s | --sample )        shift
                    sample=$1
                    ;;
        -o | --outdir )        shift
                    outdir=$1
                    ;;
        -l | --log )        shift
                    log=$1
                    ;;
        -k | --keepouts )    shift
                    keepouts=$1
                    ;;
        -v | --verbose )    shift
                    verbose=$1
                    ;;
        -z | --skip_shuf )    shift
                    skip_shuf=$1
                    ;;
        -pcr | --pcr_read_threshold )    shift
                    pcr_read_threshold=$1
                    ;;
        -ld | --levenshtein_distance )    shift
                    levenshtein_distance=$1
                    ;;
        -t | --threads_ironthrone )    shift
                    threads=$1
                    ;;
        -j | --jobs )        shift
                    jobs=$1
                    ;;
        -i | --ironthrone_dir )    shift
                    ironthrone_dir=$1
                    ;;
        -n | --total_lines )   shift
                    total_lines=$1
                    ;;
        --shuffle-fastq )    shuffle_fastq=1 ;;
    esac
    shift
done

# Main execution
validate_inputs
setup_directories
if [ "$shuffle_fastq" = 1 ]; then
    shuffle_fastq_files
fi
split_fastq_files
run_ironthrone_parallel

# Run final concatenation
echo "Combining IronThrone output..."
Rscript ${ironthrone_dir}/Combine_IronThrone_Parallel_Output.R "$result_dir" "$sample" \
    "$pcr_read_threshold" "$levenshtein_distance" "$dupcut" "$threads"

echo "All intermediate files are stored in: $outdir"
echo "Processing complete. Results in: $outdir/${sample}.summTable.concat.umi_collapsed.txt"
