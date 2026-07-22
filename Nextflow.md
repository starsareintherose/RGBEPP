## Nextflow Workflow Manager

RGBEPP can be run directly with the D executable, or through the provided Nextflow workflow. The Nextflow calls the existing RGBEPP stage commands (`clean`, `assembly`, `map`, `postmap`, `varcall`, `consen`, `align`, `ortholog`, `trim`, and `concat`) as managed workflow processes.

### Nextflow requirements

- Nextflow
- A compiled or installed RGBEPP executable available as `RGBEPP`, or an explicit executable path passed with `--rgbepp`
- The external and internal software listed in [README](./README.md#Requirements) must be available in `PATH` or passed as explicit paths to the Nextflow wrapper.

### Nextflow usage

Prepare the same input layout as the command-line pipeline:

```
.
├── 00_raw
│   ├── sample1_R1.fastq.gz
│   ├── sample1_R2.fastq.gz
│   └── ...
├── list
└── reference.aa.fasta
```

Run the workflow:

```
nextflow run starsareintherose/RGBEPP --reference reference.aa.fasta --list list --raw_dir 00_raw --rgbepp RGBEPP
```

If `genes` is omitted, Nextflow generates it from the reference FASTA headers before running RGBEPP:

```
nextflow run starsareintherose/RGBEPP --reference reference.aa.fasta --list list --genes genes --rgbepp /path/to/rgbepp
```

To enable the optional codon extraction step:

```
nextflow run starsareintherose/RGBEPP --reference reference.aa.fasta --list list --codon --rgbepp RGBEPP
```

To run through a specific RGBEPP function and stop there, use `--function`. Valid values are `all`, `clean`, `assembly`, `map`, `postmap`, `varcall`, `consen`, `codon`, `align`, `ortholog`, `trim`, and `concat`.

```
nextflow run starsareintherose/RGBEPP --reference reference.aa.fasta --list list --function trim --rgbepp RGBEPP
```

To pass an existing RGBEPP config file directly:

```
nextflow run starsareintherose/RGBEPP --reference reference.aa.fasta --list list --config config.example --rgbepp rgbepp
```

By default, Nextflow generates a complete RGBEPP config file for each process from its parameters and forwards it to RGBEPP with `-c`. Tool parameters can be command names available in `PATH` or explicit paths. Command names are resolved with `command -v` before writing the RGBEPP config because RGBEPP validates tool paths as files.

All keys from `config.example` can be provided directly through Nextflow:

```
nextflow run starsareintherose/RGBEPP \
  --reference reference.aa.fasta \
  --list list \
  --rgbepp RGBEPP \
  --fastp fastp \
  --spades spades.py \
  --diamond diamond \
  --sortdiamond sortdiamond \
  --bowtie2 bowtie2 \
  --samtools samtools \
  --bcftools bcftools \
  --exonerate exonerate \
  --macse /usr/share/java/macse.jar \
  --delstop delstop \
  --deltaxa deltaxa \
  --trimal trimal \
  --concataln concataln \
  --raw_dir 00_raw \
  --fastp_dir 01_fastp \
  --spades_dir 02_spades \
  --bowtie2_dir 03_bowtie2 \
  --bam_dir 04_bam \
  --vcf_dir 05_vcf \
  --consen_dir 06_consen \
  --macse_dir 07_macse \
  --trimal_dir 08_trimal \
  --spades_paras "--careful --phred-offset 33" \
  --diamond_blastx_paras "--ultra-sensitive" \
  --bcftools_call_paras "-mv" \
  --bcftools_norm_paras "--check-ref s" \
  --trimal_paras "-gt 0.7"
```

The full set of RGBEPP config keys accepted as Nextflow parameters is: `fastp`, `spades`, `diamond`, `sortdiamond`, `bowtie2`, `samtools`, `bcftools`, `exonerate`, `macse`, `delstop`, `deltaxa`, `trimal`, `concataln`, `raw_dir`, `fastp_dir`, `spades_dir`, `bowtie2_dir`, `bam_dir`, `vcf_dir`, `consen_dir`, `macse_dir`, `trimal_dir`, `fastp_paras`, `spades_paras`, `diamond_makedb_paras`, `diamond_blastx_paras`, `bowtie2_build_paras`, `bowtie2_paras`, `samtools_view_paras`, `samtools_fixmate_paras`, `samtools_sort_paras`, `samtools_markdup_paras`, `samtools_index_paras`, `bcftools_mpileup_paras`, `bcftools_call_paras`, `bcftools_norm_paras`, `bcftools_filter_paras`, `bcftools_index_paras`, `bcftools_consensus_paras`, `exonerate_paras`, `macse_paras`, and `trimal_paras`.

The Nextflow wrapper also accepts the RGBEPP input arguments `--list`, `--genes`, `--reference`, `--threads`, `--memory`, `--config`, `--codon`, and `--function`. The RGBEPP `-f/--functions` argument is managed internally by Nextflow: each workflow process calls the appropriate RGBEPP function for that stage. Users can still run the D executable directly when they want a manual single-stage RGBEPP command.

If `--config` is supplied, the external RGBEPP config file is used instead of the generated Nextflow config. In that case, make sure its directory keys match the corresponding Nextflow directory parameters, because Nextflow uses those parameter values to stage and collect process outputs.

Directory parameters must be relative paths and must not contain `..`. This keeps each RGBEPP step isolated inside its Nextflow work directory and avoids accidental writes outside the task sandbox.

Each process checks that its expected outputs exist before Nextflow marks the step successful. This catches failures from external tools even when RGBEPP itself exits with status 0.

Results are published to `results` by default. Nextflow work directories and cached intermediate task states remain under `work`, allowing failed or interrupted runs to be resumed with:

```
nextflow run starsareintherose/RGBEPP -resume --reference reference.aa.fasta --list list --raw_dir 00_raw --rgbepp RGBEPP
```
