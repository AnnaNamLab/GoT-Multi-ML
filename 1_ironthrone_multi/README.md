# IronThrone-Multi

<aside>

💡 This is a re-implementation of [IronThrone-GoT](https://github.com/dan-landau/IronThrone-GoT) to facilitate execution

IronThrone-Multi can process one or more genotyped targets simultaneously.

For more information, refer to the original repository:  https://github.com/dan-landau/IronThrone-GoT

</aside>

## Requirements
- [seqkit](https://bioinf.shenwei.me/seqkit/) needs to be installed

## Steps Overview

1. Concatenate fastq files if there are multiple lanes (for R1 & R2)
2. Shuffle & Split reads in the fastq files into multiple subsets for parallelization
3. Assign genotype for each read
4. Aggregate the reads for each cell

→ Output: Read counts for each UMI of each cell (Each row is a cell)

## Input

A csv file containing:

- List of R1 fastq files `[R1-fastq]`
- List of R2 fastq files `[R2-fastq]`
- List of config files `[config-files]`
- Barcode whitelist `[barcode-whitelist]`
- Expected sequence CSV `[expected-sequence-csv]`
    - A CSV file with two columns: the first column is the target name (matching the config file), and the second column is the expected sequence for that target. This file will be used to extract the target sequence for subsetting reads.
    - Example format:
      ```csv
      TARGET_NAME,EXPECTED_SEQUENCE
      DUSP2_Mut1,GCATCCACTTCCTGATTATAATCTTCTTGCCAGCGTGTTTCCAAGGGATG
      IKZF3,ACTGCTGATCGTACGATCGATCGTAGCTAGCTAGCTAGCTAGCTAGCTA
      ...
      ```
- Other arguments for IronThrone `[ironthrone-args]`
    - `prefix`: The prefix in the config file in front of the target name
    - `suffix`: The suffix in the config file after the target name (should either be `_mutation.config` or `_barcode.config`)
        - if the run is for probe-sequence, config file has suffix `_mutation.config`
        - if the run is for probe-barcode, config file has suffix `_barcode.config`
        - For regular GoT, can set it to `.config` or just omit it
    - `n_subset`: How many subsets of fastq file to split into (for parallelization). Recommend a value between 5~10
    - `umilen`: UMI length (default: 12)
    - `bclen`: Barcode length (default: 16)
      > These `umilen` and `bclen` values are specific for FRP
    - `jobs`: Maximum number of split parallel IronThrone runs
- Example csv file:
    
    ```bash
    [R1-fastqs]
    /data/frp/fastq/FRP_GENlib_L001_R1_001.fastq.gz
    /data/frp/fastq/FRP_GENlib_L002_R1_001.fastq.gz
    /data/frp/fastq/FRP_GENlib_L003_R1_001.fastq.gz
    /data/frp/fastq/FRP_GENlib_L004_R1_001.fastq.gz
    
    [R2-fastqs]
    /data/frp/fastq/FRP_GENlib_L001_R2_001.fastq.gz
    /data/frp/fastq/FRP_GENlib_L002_R2_001.fastq.gz
    /data/frp/fastq/FRP_GENlib_L003_R2_001.fastq.gz
    /data/frp/fastq/FRP_GENlib_L004_R2_001.fastq.gz
    
    [config-files]
    /data/frp/input/mutation_configs/FRP_PrimaryCells_DUSP2_Mut1_mutation.config
    /data/frp/input/mutation_configs/FRP_PrimaryCells_IKZF3_mutation.config
    /data/frp/input/mutation_configs/FRP_PrimaryCells_PLCG2_Mut3_mutation.config
    /data/frp/input/mutation_configs/FRP_PrimaryCells_TNFAIP3_Mut2_mutation.config
    
    [barcode-whitelist]
    /data/frp/input/737K-fixed-rna-profiling_barcodes.txt

    [expected-sequence-csv]
    /data/frp/input/expected_sequences.csv

    [ironthrone-args]
    prefix,FRP_PrimaryCells_
    suffix,_mutation.config
    n_subset,10
    umilen,12
    bclen,16
    mmtch,0.2
    dupcut,2
    pcr_read_threshold,0.75
    jobs,12
    ```
    

## Running IronThrone-Multi

Make sure the `ironthrone` directory is in the current directory

`run_ironthrone.sh`:

```bash
#!/bin/bash

module load perl-5.34.0-gcc-4.8.5-67vae3w
module load r-4.1.3-gcc-8.2.0-wall5zs

# === This script runs the IronThrone pipeline using the specified parameters file === #
## ! Make sure `ironthrone` directory is in the same directory as this script
##  Result will be saved in the `<id>_Results` directory

bash ./ironthrone/ironthrone_multi.sh \
    --id FRP9_GEN_$(date +"%Y_%m_%d") \
    --params ./input/params.txt

```

## Output

By default, the output directory (`<ID>_Results`) is created in the current directory (You can also specify the name of the output directory with `-o` option).

- Output directory file structure
    
    ```
    FRP9_GEN_Results/
    ├── <TARGET1>/
    │ 	├── <TARGET1>.summTable.concat.txt
    │ 	└── <TARGET1>.summTable.concat.umi_collapsed.txt    # This is the final merged output for TARGET1
    ├── <TARGET2>/
    │ 	├── ...
    ├── ...    # One directory for each TARGET
    │
    └── processed_input/    # Directory containing processed input files for IronThrone
      	├── <PREFIX>ironethrone_params.csv    # csv file containing the parameters to be used in running IronThrone
 		├── <ID>_R1.fastq    # Concatenated Fastq file
 		├── <ID>_R2.fastq    # if only single pair of fastq files was given, this would be just a copy of it
 		└── R1_R2_pairs/    # Fastq files (shuffled,) subsetted for each target (corresponding to a config file)
 		    ├── <TARGET1>/
 		    │   ├── <TARGET1>_R1_001.fastq
 		    │   ├── <TARGET1>_R2_001.fastq
 		    ├── <TARGET2>/
 		    │   ├── <TARGET2>_R1_001.fastq
 		    │   └── <TARGET2>_R2_001.fastq
 		    └── ...
     ```
    

The merged output file (`<TARGET>.summTable.concat.umi.collapsed.txt`)  is a table where each row is a cell, and for cell, it contains:

- UMI sequences (`UMI`)
- Number of reads supporting each UMI (`num.WT.in.dups`, `num.MUT.in.dups`, `num.amb.in.dups`)
- Number of UMIs supporting each genotype (`WT.calls`, `MUT.calls`, `amb.calls`)
- Example Output
    
    ![image.png](example_result.png)
