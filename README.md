# RGB EPP

Reference Genome based Exon Phylogeny Pipeline

License: GPL-2.0-only
Author: Guoyi Zhang

# Requirements

## External software 

- GNU Bash (provide cd)
- GNU coreutils (provide cp mv mkdir mv)
- GNU findutils (provide find)
- fastp
- spades.py (provided by spades)
- diamond
- java
- macse (default recognized path: /usr/share/java/macse.jar)
- GNU parallel

## Internal software

- splitfasta (default recognized path: /usr/bin/splitfasta)
- sortdiamond (default recognized path: /usr/bin/sortdiamond)

# Arguments

## Details

```
-c	--contigs	contings type: scaffolds or contigs
-g	--genes		gene file path
-f	--functions	functions type (optional): all clean 
	  		assembly fasta map pre split merge align
-h	--help		show this information
-l	--list		list file path
-m	--memory	memory settings (optional, default 16 GB)
-r	--reference	reference genome path
-t	--threads	threads setting (optional, default 8 threads)
	--macse		Macse jarfile path
	--sortdiamond	sortdiamond file path
	--splitfasta	splitfasta file path
for example: bash RGBEPP.sh -c scaffolds -f all -l list -g genes -r reference.aa.fasta 
```

## Directories Design

```
.
├── 00_raw
├── 01_fastp
├── 02_spades
├── 03_assemblied
├── 04_diamond
├── 05_pre
├── 06_split
├── 07_merge
├── 08_macse
├── genes
├── list
├── reference.aa.fasta
└── RGBEPP.sh
```

Each directory corresponds to each function. 

## Text Files

`list` is the text file containing all samples, if your raw data is following the style ${list_name}\_R1.fastq.gz and  ${list_name}\_R2.fastq.gz, ${list_name} is what you should list in `list` file. The easy way to get it in Linux/Unix system is the following command

```
cd 00_raw
ls | sed "s@_R[12].fastq.gz@@g" > ../list
cd ..
```

`genes` is the text file containing all gene names from the reference fasta file. The easy way to get it in Linux/Unix system is the following command

```
grep '>' Reference.fasta | sed "s@>@@g" > genes
```

`reference.aa.fasta` can be replaced by another other name, but it must contain reference amino acids genome in fasta format

# Progress

## RGBEPP.sh functions

 - Function clean: Quality control + trimming (fastp)
 - Function assembly: de novo assembly (spades)
 - Function fasta: gather all fasta files from assembly directories (RGBEPP.sh)
 - Function map: local nucleic acids alignment search against amino acids subject sequence (diamond)
 - Function pre: generate corresponding sequences based on blast-styled output (sortdiamond) 
 - Function split: splitting fasta sequence to directories based on the reference genome (splitfasta)
 - Function merge: merge different taxa in the same reference exon gene to one fasta (RGBEPP.sh)
 - Function align: multiple sequence align based on Condon (macse)

## Downstream process

 - concatenate sequences via SeqCombGo or catsequences or sequencematrix
 - coalescent / concatenated phylogeny


