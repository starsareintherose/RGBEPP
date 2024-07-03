#!/bin/bash

### Environment Setting

pkgver=0.0.1
DirRaw=00_raw
DirQcTrim=01_fastp
DirAssembly=02_spades
DirFasta=03_contig
DirMap=04_map
DirPre=05_pre
DirSplit=06_split
DirMerge=07_merge
DirAlign=08_align

PathSplitfsata=~/Downloads/PhD/wes/splitfasta-cpp
PathMacse=/usr/share/java/macse.jar
PathSortdiamond=/home/guoyi/Downloads/PhD/wes/sortdiamond

HELP=false

### Get some arrays

ARGS=$(getopt -o c:,f:,g:,h,l:,m:,r:,t: --long contigs:,genes:,functions:,help,list:,memory:,reference:,threads: -n 'RGBEPP.sh' -- "$@")
if [ $? != 0 ]; then
    echo "Failed to parse options." >&2
    exit 1
fi
eval set -- "$ARGS"

while true; do
    case "$1" in
        -c|--contigs)
            case "$2" in
                "") ARG_C='scaffolds'; shift 2 ;;
                *) ARG_C=$2; shift 2 ;;
            esac ;;
        -g|--genes)
            case "$2" in
                "") shift 2 ;;
                *) ARG_G=$2; shift 2 ;;
            esac ;;
        -f|--functions)
            case "$2" in
                "") ARG_F='all'; shift 2 ;;
                *) ARG_F=$2; shift 2 ;;
            esac ;;
        -h|--help)
            echo -e "\t\t\t\t\t\t\tRGB EPP\n\t\t\t\t\tReference Genome based Exon Phylogeny Pipeline\n \
		       Version: $pkgver\n \
		       License: GPL-3.0-only\n \
		       Author: Guoyi Zhang\n \
                      -c\t--contigs\tcontings type: scaffolds or contigs\n \
                      -g\t--genes\t\tgene file path\n \
                      -f\t--functions\tfunctions type (optional): all clean \n \
                        \t  \t\tassembly fasta map pre split merge align\n \
                      -h\t--help\t\thelp: show this information\n \
                      -l\t--list\t\tlist file path\n \
                      -m\t--memory\tmemory settings (optional, default 16 GB)\n \
                      -r\t--reference\treference genome path\n \
                      -t\t--threads\tthreads setting (optional, default 8 threads)\n \
                      for example: bash $0 -c scaffolds -f all -l list -g genes -r Reference.exons.aa.fas \n"
            HELP=true
            shift ;;
        -l|--list)
            case "$2" in
                "") shift 2 ;;
                *) ARG_L=$2; shift 2 ;;
            esac ;;
        -m|--memory)
            case "$2" in
                "") ARG_M=16; shift 2 ;;
                *) ARG_M=$2; shift 2 ;;
            esac ;;
        -r|--reference)
            case "$2" in
                "") shift 2 ;;
                *) ARG_R=$2; shift 2 ;;
            esac ;;
        -t|--threads)
            case "$2" in
                "") ARG_T=8; shift 2 ;;
                *) ARG_T=$2; shift 2 ;;
            esac ;;
        --) shift; break ;;
        *) echo "Internal error!"; exit 1 ;;
    esac
done

### Get and check some arguments

if [ "$HELP" = false ]; then
    if [ -z "$ARG_L" ]; then
        echo "List argument can't be empty"
        exit 1
    fi

    readarray -t full_names < "$ARG_L"
    readarray -t genes < "$ARG_G"
    length_fn=${#full_names[@]}
    length_gn=${#genes[@]}
fi

### Quality control && Trimming

if [ "$ARG_F" = "all" ] || [ "$ARG_F" = "clean" ]; then

	## Prepare
	mkdir -p $DirQcTrim

	## Quality control and trimming using fastp
	for (( i=0; i<$length_fn; i++ )); do
		fastp -i $DirRaw/${full_names[$i]}_R1.fastq.gz -I $DirRaw/${full_names[$i]}_R2.fastq.gz -j $DirQcTrim/${full_names[$i]}.json -h $DirQcTrim/${full_names[$i]}.html -o $DirQcTrim/${full_names[$i]}_R1.fastq.gz -O $DirQcTrim/${full_names[$i]}_R2.fastq.gz -w $ARG_T
	done

fi

### De novo assembly

if [ "$ARG_F" = "all" ] || [ "$ARG_F" = "assembly" ]; then

	## Prepare
	mkdir -p $DirAssembly
	
	## De novo assembly using spades
	for (( i=0; i<$length_fn; i++ )); do
		mkdir -p $DirAssembly/${full_names[$i]} 
		spades.py --pe1-1 $DirQcTrim/${full_names[$i]}_R1.fastq.gz --pe1-2 $DirQcTrim/${full_names[$i]}_R2.fastq.gz -t $ARG_T -m $ARG_M --careful --phred-offset 33 -o $DirAssembly/${full_names[$i]}
			# -k 96,107,117,127 \
	done

fi

### Moving scaffords or Contigs out

if [ "$ARG_F" = "all" ] || [ "$ARG_F" = "fasta" ]; then

	## Check if the contigs type is specified
	if [ -z "$ARG_C" ] ; then
		echo "Argument of contig type missing."
		exit 1
	fi

	## Prepare
	mkdir -p $DirFasta

	## Move the assemblied fasta file to the folder
	if [ "$ARG_C" = "scaffolds" ] || [ "$ARG_C" = "contigs" ] ; then
		for (( i=0; i<$length_fn; i++ )); do
			cp $DirAssembly/${full_names[$i]}/$ARG_C.fasta $DirFasta/$ARG_C/${full_names[$i]}.fasta
		done
	fi

fi

### Mapping

if [ "$ARG_F" = "all" ] || [ "$ARG_F" = "map" ]; then

	## Check if the reference or contigs type is specified
	if [ -z "$ARG_R" ] || [ -z "$ARG_C" ] ; then
		echo "Argument of reference or contig type missing."
		exit 1
	fi

	## Prepare
	mkdir -p $DirMap

	## Index reference database
	cd $DirFasta/$ARG_C
	diamond makedb --db Reference --in $ARG_R
	cd -

	## Blastx for mapping DNA sequences to protein reference sequence
	cd $DirFasta/$ARG_C
	for (( i=0; i<$length_fn; i++ )); do		
		diamond blastx -d Reference.dmnd -q ${full_names[$i]}.fasta -o ${full_names[$i]}.m8 \
		--outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen gaps ppos qframe qseq
	# subject: reference; query: align-aimed 
	#1:	qseqid: Query Seq-id
	#2:	sseqid: Subject Seq - id
	#3:	pident: Percentage of identical matches
	#4:	length: Alignment length
	#5:	mismatch: Number of mismatches
	#6:	gapopen: Number of gap openings
	#7:	qstart: Start of alignment in query
	#8:	qend: End of alignment in query
	#9:	sstart: Start of alignment in subject
	#10:	send: End of alignment in subject
	#11:	evalue: Expect value
	#12:	bitscore: Bit score
	#13:	qlen: Query sequence length 比对序列长度
	#14:	slen: Subject sequence length
	#15:	gaps: Total number of gaps
	#16:	ppos: Percentage of positive - scoring matches
	#17:	qframe: Query frame (frames in ECPP.sh)
	#18:	qseq: Aligned part of query sequence

	done
	cd -

	mv $DirFasta/$ARG_C/*.m8 $DirMap

fi

if [ "$ARG_F" = "all" ] || [ "$ARG_F" = "pre" ]; then
	mkdir -p $DirPre
	for (( i=0; i<$length_fn; i++ )); do	
		$PathSortdiamond $DirMap/${full_names[$i]}.m8 $DirPre/${full_names[$i]}.fasta
	done
fi


if [ "$ARG_F" = "all" ] || [ "$ARG_F" = "split" ]; then
	mkdir -p $DirSplit
	cd $DirPre
	for (( i=0; i<$length_fn; i++ )); do
		$PathSplitfsata ${full_names[$i]}.fasta
	done
	find . -mindepth 1 -maxdepth 1 -type d -exec mv {} ../$DirSplit \;
	cd -
fi

if [ "$ARG_F" = "all" ] || [ "$ARG_F" = "merge" ]; then

	## Check if the genes is specified
	if [ -z "$ARG_G" ] ; then
		echo "Argument of genes list missing."
		exit 1
	fi

	mkdir -p $DirMerge
	cd $DirSplit
	for (( i=0; i<$length_gn; i++ )) 
	do
		cd ${genes[$i]}
		cat * > ../${genes[$i]}.fasta
		cd ..
	done
	mv *.fasta ../$DirMerge
	cd -
	
fi

if [ "$ARG_F" = "all" ] || [ "$ARG_F" = "align" ]; then

	## Check if the genes is specified
	if [ -z "$ARG_G" ] ; then
		echo "Argument of genes list missing."
		exit 1
	fi

	mkdir -p $DirAlign
	mkdir -p $DirAlign/AA && mkdir -p $DirAlign/NT
	cd $DirMerge
	for (( i=0; i<$length_gn; i++ ))
	do
		java -jar $PathMacse -prog alignSequences -seq ${genes}.fasta -out_AA ../$DirAlign/AA/$genes.fasta -out_NT ../$DirAlign/NT/$genes.fasta 
	done
	cd -

fi

