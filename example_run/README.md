# GoT-Multi-ML example run

This folder is a lightweight example run for the current GoT-Multi-ML pipeline. It is based on the toy FFPE example used during development, but large files are intentionally omitted so the repository stays GitHub-friendly.

## Download required large input files

Download the toy FASTQ files and the GEX object from:

- [Google Drive example input folder](https://drive.google.com/drive/folders/1XQVWG1rzPFKWlSdDfZLor6HLOUceQiwb?usp=share_link)

Place them here before running:

- `input/fastq/CLL8_GENlib_S4_L001_R1_001_subset100k.fastq.gz`
- `input/fastq/CLL8_GENlib_S4_L001_R2_001_subset100k.fastq.gz`
- `input/fastq/CLL8_GENlib_S4_L002_R1_001_subset100k.fastq.gz`
- `input/fastq/CLL8_GENlib_S4_L002_R2_001_subset100k.fastq.gz`
- `input/sc_gex.h5ad`

## Inputs already included in this folder

These small files are already present:

- `input/configs/gene_seq_configs/`
- `input/configs/probe_bc_configs/`
- `input/params_gene_seq.txt`
- `input/params_probe_bc.txt`
- `input/expected_sequence.csv`
- `input/target_info.csv`
- `input/737K-fixed-rna-profiling_barcodes.txt`

## How to run

Step 1: run GoT-Multi initial genotyping

```bash
cd example_run
bash 1_run_got_multi_ml.sh
```

Step 2: run optional ML denoising

```bash
cd example_run
bash 2_run_got_multi_ml_denoise.sh
```

Run both steps:

```bash
cd example_run
bash 3_run_all.sh
```

## Main outputs

After step 1:

- `got_multi_ml_run/ironthrone_out_pp.csv`
- `got_multi_ml_run/initial_genotype_counts_by_target.pdf`
- `got_multi_ml_run/GeneSeq_*_Results/`
- `got_multi_ml_run/ProbeBC_*_Results/`

After step 2:

- `got_multi_ml_run/ml_denoise/ironthrone_out_pp_ml_input.csv`
- `got_multi_ml_run/ml_denoise/got_multi_out_refined.csv`
- `got_multi_ml_run/ml_denoise/refined_genotype_proportions_in_expected_cell_group.pdf`
- `got_multi_ml_run/ml_denoise/initial_genotype_proportions_in_expected_cell_group.pdf`

The final refined genotype column in `got_multi_out_refined.csv` is `genotype_final`.

## Notes

- This example uses the toy subset FASTQs, not the full FFPE FASTQs.
- Large generated outputs such as `got_multi_ml_run/` and intermediate `processed_input/` folders are intentionally not included here.
- The params files use paths relative to `example_run/`, so run the wrapper scripts from this folder as shown above.
