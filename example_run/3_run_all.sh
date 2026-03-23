#!/bin/bash
set -euo pipefail

bash 1_run_got_multi_ml.sh
bash 2_run_got_multi_ml_denoise.sh
