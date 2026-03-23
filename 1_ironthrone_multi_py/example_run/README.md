# IronThrone-Multi example run

This folder is a lightweight example run for standalone IronThrone-Multi. It uses the toy FFPE subset FASTQs that were used during development, but large files are intentionally omitted so the repository stays GitHub-friendly.

## Download required large input files

Download the toy FASTQ files from:

- [Google Drive example input folder](https://drive.google.com/drive/folders/1XQVWG1rzPFKWlSdDfZLor6HLOUceQiwb?usp=share_link)

Place these files here before running:

- `input/fastq/CLL8_GENlib_S4_L001_R1_001_subset100k.fastq.gz`
- `input/fastq/CLL8_GENlib_S4_L001_R2_001_subset100k.fastq.gz`
- `input/fastq/CLL8_GENlib_S4_L002_R1_001_subset100k.fastq.gz`
- `input/fastq/CLL8_GENlib_S4_L002_R2_001_subset100k.fastq.gz`

## Inputs already included in this folder

These small files are already present:

- `input/configs/mutation_configs/`
- `input/params.txt`
- `input/expected_sequence.csv`
- `input/737K-fixed-rna-profiling_barcodes.txt`

## What this example shows

Standalone IronThrone-Multi runs one sequencing library at a time. This example is a single standalone run using the toy mutation-target panel.

## How to run

```bash
cd example_run
bash 1_run_example.sh
```

## Main outputs

After the run, the output directory contains per-target outputs plus a combined standalone result.

- `example_run_output/`
- `example_run_output/ironthrone_out_pp.csv`
- `example_run_output/initial_genotype_counts_by_target.pdf`
- `example_run_output/processed_input/`
- `example_run_output/<TARGET>/<TARGET>.ironthrone_out_pp.csv`

These standalone outputs contain the single-run `genotype` column. They do not create `genotype_merged`, which is added later only in the paired GoT-Multi wrapper workflow.

## Notes

- This example uses the toy subset FASTQs, not the full FFPE FASTQs.
- Large generated outputs such as `example_run_output/` and intermediate `processed_input/` folders are intentionally not included here.
- The params file uses paths relative to `example_run/`, so run the wrapper script from this folder as shown above.
