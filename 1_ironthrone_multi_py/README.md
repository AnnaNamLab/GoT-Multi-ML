# IronThrone-Multi

IronThrone-Multi is a standalone genotyping workflow for GoT-style targeted single-cell sequencing data.
It takes paired FASTQ files, one or more target config files, and a barcode whitelist, then produces per-target and combined genotype tables.

It is built around the original [IronThrone-GoT](https://github.com/dan-landau/IronThrone-GoT) logic, but extends it for multi-target processing, streamlines input preparation, and adds Python-based post-processing for easier downstream use.

If you are using this module as part of the full GoT-Multi-ML pipeline, see the [GoT-Multi-ML README](../README.md).

## What IronThrone-Multi Produces

For each target in one sequencing run, IronThrone-Multi summarizes evidence for:

- wild-type support
- mutant support
- ambiguous support
- per-cell genotype calls

The main output is a processed CSV table for each target plus a combined `ironthrone_out_pp.csv` across all targets in the run.

## Quick Start

```bash
bash ./ironthrone_multi.sh \
  --id MY_RUN \
  --params ./input/params.txt
```

If you want to choose the result directory explicitly:

```bash
bash ./ironthrone_multi.sh \
  --id MY_RUN \
  --params ./input/params.txt \
  --outdir ./results/MY_RUN
```

## What You Need Before Running

- paired R1 and R2 FASTQ files
- one `.config` file per target, or one directory containing the `.config` files for the run
- a barcode whitelist file
- a params file describing the run
- optionally, an `expected_sequence.csv` file

<details>
<summary>Show requirements</summary>

### External tools

- `seqkit`
- GNU `parallel`
- Python 3 in `PATH`, or set through `PYTHON_BIN`

### Python-side workflow

IronThrone-Multi uses the bundled Python scripts in this repository for:

- expected-sequence handling
- output post-processing
- result aggregation

</details>

## Params File

IronThrone-Multi is driven by one params file.
The file is organized into sections.

```text
[R1-fastqs]
/path/to/R1_L001.fastq.gz
/path/to/R1_L002.fastq.gz

[R2-fastqs]
/path/to/R2_L001.fastq.gz
/path/to/R2_L002.fastq.gz

[config-files]
/path/to/configs/

[barcode-whitelist]
/path/to/barcodes.txt

[expected-sequence-csv]
/path/to/expected_sequence.csv

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

<details>
<summary>Show section-by-section details</summary>

### `[R1-fastqs]` and `[R2-fastqs]`

List the paired FASTQ files for the run, one path per line.
Multiple lanes are allowed.

### `[config-files]`

You can provide either:

- one config path per line
- one directory path containing all `.config` files for the run

Only one directory is allowed in this section.
If you use a directory, IronThrone-Multi automatically collects all `*.config` files in that directory.

### `[barcode-whitelist]`

Path to the valid cell-barcode whitelist for your assay.

### `[expected-sequence-csv]`

Optional.
If provided, it should contain two columns with no header:

```csv
TARGET_NAME,EXPECTED_SEQUENCE
```

If it is missing entirely, IronThrone-Multi derives expected sequences from the config files.
If it is provided but missing some targets, IronThrone-Multi creates a new completed CSV in the output directory and leaves the user-supplied file unchanged.

### `[ironthrone-args]`

This section controls the run.
Common keys include:

- `prefix`
- `suffix`
- `n_subset`
- `umilen`
- `bclen`
- `mmtch`
- `dupcut`
- `pcr_read_threshold`
- `jobs`

</details>

## Config File Format

Each target config file describes the local sequence context around one target.
A standard config has four rows:

1. shared left-side sequence
2. shared right-side sequence or shared sequence block
3. WT-specific sequence
4. MUT-specific sequence

Each row stores:

- sequence
- start position
- end position

Example:

```text
AGCTGGACCA        1   10
CTGGTCAGGGGACTC   11  25
A                 26  26
G                 26  26
```

The exact biological interpretation of the four rows can differ between assay designs.
What matters for IronThrone-Multi is that the config captures the shared and allele-specific sequence blocks with valid positions.

## Expected Sequence Logic

IronThrone-Multi uses an expected sequence during the initial target-specific read-subsetting step.
That sequence is derived from the config file using the shared sequence blocks and their positions.

Current rule:

- find the longest contiguous sequence block shared by WT and MUT
- if two shared blocks are adjacent by position, concatenate them
- if a user-provided CSV is missing a target, derive the missing target from the config file and write a completed CSV copy

A helper is provided in `../utils/`:

```bash
python ../utils/build_expected_sequence_csv.py \
  --output ./input/expected_sequence.csv \
  /path/to/configs/TARGET1.config \
  /path/to/configs/TARGET2.config
```

<details>
<summary>Show expected-sequence examples</summary>

### Continuous shared block example

```text
TTAGTCCCCT       26  35
GAGGAATGGCCTCAG  36  50
G                25  25
T                25  25
```

Expected sequence:

```text
TTAGTCCCCTGAGGAATGGCCTCAG
```

### Non-contiguous shared block behavior

If the shared sequence pieces are not contiguous by position, the helper keeps the longest shared block.

</details>

## Outputs

By default, the output directory is `<ID>_Results`.
You can override that with `--outdir`.

Main outputs:

- `<TARGET>/<TARGET>.summTable.concat.txt`
- `<TARGET>/<TARGET>.summTable.concat.umi_collapsed.txt`
- `<TARGET>/<TARGET>.ironthrone_out_pp.csv`
- `ironthrone_out_pp.csv`
- `initial_genotype_counts_by_target.pdf`
- `processed_input/`

The processed CSV tables contain the single-run `genotype` column.
Standalone IronThrone-Multi does not create `genotype_merged`.

<details>
<summary>Show output details</summary>

### Per-target processed CSV

`<TARGET>.ironthrone_out_pp.csv` contains the cleaned single-target output for one run.
Typical columns include:

- `BC`
- `UMI`
- `wt_calls`
- `mut_calls`
- `amb_calls`
- `wt_reads`
- `mut_reads`
- `amb_reads`
- `genotype`

### Combined processed CSV

`ironthrone_out_pp.csv` stacks all target-level processed outputs into one table.
This is convenient for downstream plotting or integration.

### `processed_input/`

This folder stores the prepared FASTQs and cached per-target parameter CSVs used by the run.
It is especially useful if you need to rerun a later stage without repeating the whole input-prep step.

</details>

## Typical Workflow

1. prepare or collect the target config files
2. prepare the params file
3. optionally prepare `expected_sequence.csv`
4. run `ironthrone_multi.sh`
5. inspect the per-target and combined processed CSV outputs

## Related Resources

- Original IronThrone-GoT: [https://github.com/dan-landau/IronThrone-GoT](https://github.com/dan-landau/IronThrone-GoT)
- GoT-Multi-ML wrapper around IronThrone-Multi: [../README.md](../README.md)
