# RGB EPP

Reference Genome based Exon Phylogeny Pipeline

License: GPL-2.0-only

Author: Guoyi Zhang

## Requirements

### External software 

- fastp
- spades.py (provided by spades)
- diamond
- bowtie2
- samtools
- bcftools
- exonerate (optional, only for --codon)
- java
- macse (default path: `/usr/share/java/macse.jar`)
- trimal

### Internal software

- sortdiamond
- delstop
- deltaxa
- concataln

RGBEPP resolves command-name tool settings (including the internal utilities above) using `$PATH`. They do not need to be installed in `/usr/bin`. You can also provide an explicit executable path in `config.example` or with the corresponding command-line option.

For MACSE, you should specify which JAR file you prefer. The default is `/usr/share/java/macse.jar`.

## Build with DUB


Before compiling from source, install `dub`, `ldc`, and `rdmd`.

```bash
# Arch Linux / Manjaro
sudo pacman -S dub ldc dtools

# Debian / Ubuntu
sudo apt install dub ldc
sudo snap install dmd --classic  # provides rdmd

# Fedora
sudo dnf install dub ldc

# Rocky Linux / CentOS
sudo dnf install dub ldc
## you should also install snapd before using snap
sudo snap install dmd --classic # provides rdmd

# openSUSE Tumbleweed
sudo zypper install dub dmd
```

The D tools can be built with LDC through `dub.sdl`:

This builds `RGBEPP`, `RGBEPP_refmix`, `sortdiamond`, `splitfasta`, `countTaxa`, `delstop`, `deltaxa`, and `concataln` into `bin/` using `ldc2 -O2`.

To build all executable:

```
rdmd build_all.d
```

To request static linking on Linux:

```
DUB_BUILD_TYPE=release-static rdmd build_all.d
```

Static linking requires static versions of the D runtime, Phobos, libc, and other system libraries on the build machine.

## Run it directly

To run it directly, you can use the following command with following parameters.

```
RGBEPP -f all -l list -t 8 -r reference.fasta  
```

It also supports Nextflow. Users can run the workflow directly from GitHub:

```
nextflow run starsareintherose/RGBEPP --reference reference.aa.fasta --list list --raw_dir 00_raw --rgbepp RGBEPP
```

Please refer [this tutorial](./Nextflow.md) for more details.

### Details

```
    -c	--config	config file for software path (optional)
    -g	--genes		gene file path (optional, if -r is specified)
    -f	--functions	functions type (optional): all clean assembly map postmap 
      	           	varcall consen codon align ortholog trim concat
    -h	--help		show this information
    -l	--list		list file path
    --log		log file path (optional, default rgbepp.log)
    -m	--memory	memory settings (optional, default 16 GB)
    -r	--reference	reference genome path
    -t	--threads	threads setting (optional, default 8 threads)
    --codon		Only use the codon region (optional)
    --fastp		Fastp path (optional)
    --spades		Spades python path (optional)
    --diamond		Diamond python path (optional)
    --sortdiamond	SortDiamond python path (optional)
    --bowtie2		Bowtie2 path (optional)
    --samtools		Samtools path (optional)
    --bcftools		Bcftools path (optional)
    --exonerate		Exonerate path (optional)
    --macse		Macse jarfile path (optional)
    --delstop		Delstop path (optional)
    --deltaxa		Deltaxa path (optional)
    --trimal		Trimal path (optional)
    --concataln		Concataln path (optional)
    for example: ./RGBEPP -f all -l list -t 8 -r reference.fasta 
```

### Directories Design

The tree below shows the default run layout, not a set of output directories that must all be created by hand. Run RGBEPP from the project/run directory; the pipeline creates the output directories as each corresponding stage runs. For a full run (`-f all`), this includes `01_fastp` through `08_trimal` (and intermediate subdirectories). If you run only one stage, its required inputs from earlier stages must already exist.

You must provide the input data yourself: create `00_raw` (or set `raw_dir` in the config) and put the paired FASTQ files there. The directory listed in `raw_dir` is not created automatically. Also provide the sample list, reference FASTA, and—if needed—the gene list as files; their paths are supplied with `-l`, `-r`, and `-g`. The `genes` entry in the illustration represents the gene-list input file, not an output directory. You do not need to create the illustrated output directories in advance.

By default, output folders are created relative to the current working directory. A config file can change the stage directories (for example, `fastp_dir`, `spades_dir`, and `macse_dir`); relative paths in the config are also relative to the directory from which RGBEPP is run.

```
.
├── 00_raw
├── 01_fastp
│   ├── <sample>_R1.fastq.gz
│   ├── <sample>_R2.fastq.gz
│   └── <sample>.html / <sample>.json
├── 02_spades
│   ├── <sample>/                 # SPAdes per-sample work files
│   ├── scaffolds/<sample>.fasta
│   ├── contigs/<sample>.fasta
│   └── fasta/                    # reference and mapped scaffolds
├── 03_bowtie2
│   ├── index/                    # Bowtie2 indexes
│   └── <sample>.bam              # initial alignments
├── 04_bam
│   └── <sample>.bam[.bai]        # processed BAM and index
├── 05_vcf
│   └── <sample>.vcf.gz[.csi]
├── 06_consen
│   ├── taxa/<sample>.fasta
│   └── gene/<gene>.fasta
├── 07_macse
│   ├── AA/<gene>.fasta
│   ├── NT/<gene>.fasta
│   ├── AA_out/ and NT_out/       # ortholog/stop-codon processing intermediates
│   └── taxa/                     # ortholog-search inputs/results
├── 08_trimal
│   ├── NT/<gene>.fasta           # trimmed gene alignments
│   └── concat.fasta.fasta        # concatenated alignment
├── list
├── genes
├── reference.aa.fasta
└── RGBEPP
```

The entries under the numbered directories are examples of important intermediate and final files; exact files depend on the stages run and input samples/genes. For example, the assembly stage writes per-sample SPAdes work directories as well as `02_spades/scaffolds/` and `02_spades/contigs/`; mapping uses `02_spades/fasta/` and creates indexes under `03_bowtie2/index/`. `00_raw` holds the input FASTQ files and must be prepared by the user; the later numbered directories are generated by the pipeline as those stages run.

### Text Files

`list` is the text file containing all samples, if your raw data is following the style ${list_name}\_R1.fastq.gz and  ${list_name}\_R2.fastq.gz, ${list_name} is what you should list in `list` file. The easy way to get it in Linux/Unix system is the following command

```
cd 00_raw
ls *.fastq.gz | sed "s@_R[12].fastq.gz@@g" | uniq > ../list
cd ..
```

`genes` is the text file containing all gene names from the reference fasta file. The easy way to get it in Linux/Unix system is the following command

```
grep '>' Reference.fasta | sed "s@>@@g" > genes
```

`reference.aa.fasta` can be replaced by another other name, but it must contain reference amino acids genome in fasta format

## Process

### RGBEPP functions

 - Function clean: Quality control + trimming (fastp)
 - Function assembly: de novo assembly (spades)
 - Function map: local nucleic acids alignment search against amino acids subject sequence (diamond, sortdiamond), mapping raw reads to its scaffolds sequences (bowtie2) 
 - Function postmap: Sorting and marking the read read alignment (samtools)
 - Function varcall: variant calling and filtering (bcftools) 
 - Function consen: get consensus fasta file from vcf files (bcftools), then sort sequences based on gene name and taxa name (RGBEPP)
 - Function codon (optional): only extract the exon sequence (exonerate)
 - Function align: multiple sequence align based on condon (macse)
 - Function ortholog: remove paralog based on Reciprocal Best Hits (RBH, diamond) 
 - Function trim: trimming based on codon (trimal, delstop)
 - Function concat: Concatenate sequences based on genes (concataln)

### Arguments reuqirements for functions

| Functions | -g/--gene | -l/--list | -r/--reference |
| --------- | --------- | --------- | -------------- | 
| clean | | ✔ | |
| assembly | | ✔ | |
| map | | ✔ | ✔ |
| postmap | | ✔ | |
| varcall | | ✔ | |
| consen | ✔ | ✔ | |
| codon | ✔ | | ✔ |
| align | ✔ | | |
| ortholog | |  ✔ |  ✔ |
| trim | ✔ | | |
| concat | ✔ | | |


### Downstream process

 - concatenate sequences via SeqCombGo or catsequences or sequencematrix
 - coalescent / concatenated phylogeny

## Inner software

### sortdiamond

Usage: `sortdiamond diamond_output.m8 generated.fasta sseq,qstart,qend,bitscore/evalue,qseq(optional, default 1,6,7,11,17, start from 0) bitscore/evalue(optional, default bitscore)`

Default sseq is column 2, qstart is column 8, etc.

Diamond default output format (--outfmt 6) does not contain qseq, you must custom the output format under output format 6. 

### delstop

Usage: `delstop <fasta_aa> <fasta_nt> --delete`

Delete StopCondon generated by Macse. fasta_aa and fasta_nt should be macse output files, `--delete` should be used when downstream software is tirmal

### splitfasta

Usage: `splitfasta sample.fasta`

It always creates directories in the path that you run the splitfasta, and puts split fasta into the directory.

### deltaxa

Usage: `deltaxa <fasta> <taxa1> [taxa2 taxa3 ...]`

It deletes the one or multiple sequences of one fasta file.

### concataln

Usage: `concataln <prefix> file1 file2 ...`

Concatenates multiple sequence alignments (MSAs). It will generate the `<prefix>.fasta` and `<prefix>.partitions.txt`
