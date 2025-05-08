#!/bin/bash
# Run IronThrone-Denoise (Step 3)
# Usage: bash run_ironthrone_denoise.sh

# Edit the following variables as needed
INPUT_FEATURES="2_ironthrone_pp/output/ironthrone_out_pp.csv"
OUTDIR="3_ironthrone_denoise/output"
MAX_ALLOWED_FPR=0.03
ALPHA=0.8
MODELS_TO_RUN="logistic_regression,random_forest,knn,naive_bayes,xgboost,mlp,gradient_boosting,hist_gradient_boosting,adaboost"

python3 3_ironthrone_denoise/ml_genotyping.py \
    --input-features "$INPUT_FEATURES" \
    --outdir "$OUTDIR" \
    --max-allowed-fpr $MAX_ALLOWED_FPR \
    --alpha $ALPHA \
    --models-to-run "$MODELS_TO_RUN"
# Output: $OUTDIR/aggregated/final_prediction_results_alpha_${ALPHA}.csv
