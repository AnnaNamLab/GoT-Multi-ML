#!/bin/bash
# Run IronThrone-Multi (Step 1)
# Output will be in <ID>_Results/

# Edit the following variables as needed
IRONTHRONE_DIR="/path/to/1_ironthrone_multi"

# Run for Gene Sequence
ID="GeneSeq_$(date +"%Y_%m_%d")"
PARAMS_FILE="input/params_gene_seq.txt"
bash ${IRONTHRONE_DIR}/ironthrone_multi.sh \
    --id "$ID" \
    --params "$PARAMS_FILE" \
    --iron_throne_dir "${IRONTHRONE_DIR}" \
    # --skip-input-prep


# Run for Probe Barcode
ID="ProbeBC_$(date +"%Y_%m_%d")"
PARAMS_FILE="input/params_probe_bc.txt"
bash ${IRONTHRONE_DIR}/ironthrone_multi.sh \
    --id "$ID" \
    --params "$PARAMS_FILE" \
    --iron_throne_dir "${IRONTHRONE_DIR}" \
    # --skip-input-prep
