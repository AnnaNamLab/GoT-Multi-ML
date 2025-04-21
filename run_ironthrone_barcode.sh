#!/bin/bash

#SBATCH --partition=scu-cpu
#SBATCH --job-name=ironethrone_frp9_barcode_241210
#SBATCH --time=7-00:00:00   # D-HH/MM/SS
#SBATCH --mem=512G
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=mip4018@med.cornell.edu


module load perl-5.34.0-gcc-4.8.5-67vae3w
module load r-4.1.3-gcc-8.2.0-wall5zs

# === This script runs the IronThrone pipeline using the specified parameters file === #
##  Result will be saved in the `<id>_Results` directory

# ironthrone_multi \
bash ./ironthrone/ironthrone_multi.sh \
    --id FRP_MUTATION_$(date +"%Y_%m_%d") \
    --params ./input/params_barcode.txt \
    --iron_throne_dir "./ironthrone" \
     --skip-input-prep
