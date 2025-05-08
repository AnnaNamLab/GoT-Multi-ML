#!/bin/bash
# Run IronThrone-PP (Step 2)
# Usage: bash run_ironthrone_pp.sh

# Edit the following variables as needed
MUTDIR="<path_to_mutation_results>"  # e.g. 1_ironthrone_multi/FRP_RUN_2025_05_07/ or similar
BARDIR="<path_to_barcode_results>"   # e.g. 1_ironthrone_multi/FRP_RUN_2025_05_07/ or similar
CELL_GROUP="cell_group"
NEG_CTRL_GROUP=""  # comma-separated if multiple
GEX="<path_to_gex.h5ad>"
TARGET_INFO="<path_to_target_info.csv>"
OUTDIR="2_ironthrone_pp/output"
ADDITIONAL_FEATURES=""  # comma-separated if needed

python3 2_ironthrone_pp/process_ironthrone_output.py \
    --mutdir "$MUTDIR" \
    --bardir "$BARDIR" \
    --cell-group "$CELL_GROUP" \
    --negative-control-cell-group "$NEG_CTRL_GROUP" \
    --gex "$GEX" \
    --target-info "$TARGET_INFO" \
    --outdir "$OUTDIR" \
    --additional-features "$ADDITIONAL_FEATURES"
# Output: $OUTDIR/ironthrone_out_pp.csv
