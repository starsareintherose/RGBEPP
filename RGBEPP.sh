#!/bin/bash

### Environment Setting

pkgver=0.0.1

DirHome=$(pwd)
DirRaw=$DirHome/00_raw
DirQcTrim=$DirHome/01_fastp
DirAssembly=$DirHome/02_spades
DirFasta=$DirHome/03_assemblied
DirMap=$DirHome/04_diamond
DirPre=$DirHome/05_pre
DirSplit=$DirHome/06_split
DirMerge=$DirHome/07_merge
DirAlign=$DirHome/08_macse

PathSplitfsata=/usr/bin/splitfasta-cpp
PathMacse=/usr/share/java/macse.jar
PathSortdiamond=/usr/bin/sortdiamond

ARG_C='scaffolds'
ARG_M=16
ARG_T=8

### Get some arrays

show_help(){
#            echo -e "\t\t\t\t\t\t\t\033[0;31mR\033[0m\033[0;92mG\033[0m\033[0;94mB\033[0m \033[0;33mE\033[0m\033[0;94mP\033[0m\033[0;33mP\033[0m\n\t\t\t\t\tReference Genome based Exon Phylogeny Pipeline\n \
            echo -e "\t\t\t\t\t\t\t\033[0;47;31mR\033[0m\033[0;47;92mG\033[0m\033[0;47;94mB\033[0m\033[0;47m \033[0m\033[0;47;33mE\033[0m\033[0;47;94mP\033[0m\033[0;47;33mP\033[0m\n\t\t\t\t\tReference Genome based Exon Phylogeny Pipeline\n \
		       Version: $pkgver\n \
		       License: GPL-2.0-only\n \
		       Author: Guoyi Zhang\n \
                      -c\t--contigs\tcontings type: scaffolds or contigs\n \
                      -g\t--genes\t\tgene file path\n \
                      -f\t--functions\tfunctions type (optional): all clean \n \
                        \t  \t\tassembly fasta map pre split merge align\n \
                      -h\t--help\t\tshow this information\n \
                      -l\t--list\t\tlist file path\n \
                      -m\t--memory\tmemory settings (optional, default 16 GB)\n \
                      -r\t--reference\treference genome path\n \
                      -t\t--threads\tthreads setting (optional, default 8 threads)\n \
                        \t--macse\t\tMacse jarfile path\n \
                        \t--sortdiamond\tsortdiamond file path\n \
                        \t--splitfasta\tsplitfasta file path\n \
                      for example: bash $0 -c scaffolds -f all -l list -g genes \ \n \
			\t    -r reference.aa.fas \n"
}

if [ $# -eq 0 ]; then
    show_help 
    exit 1
else

ARGS=$(getopt -o c:,f:,g:,h,l:,m:,r:,t: --long contigs:,genes:,functions:,help,list:,memory:,reference:,threads:,macse:,sortdiamond:,splitfasta: -n 'RGBEPP.sh' -- "$@")
if [ $? != 0 ]; then
    echo "Failed to parse options." >&2
    exit 1
fi
eval set -- "$ARGS"

while true; do
    case "$1" in
        -c|--contigs)
            case "$2" in
                "") shift 2 ;;
                *) ARG_C=$2; shift 2 ;;
            esac ;;
        -g|--genes)
            case "$2" in
                "") shift 2 ;;
                *) ARG_G=$2; shift 2 ;;
            esac ;;
        -f|--functions)
            case "$2" in
                "") shift 2 ;;
                *) ARG_F=$2; shift 2 ;;
            esac ;;
        -h|--help)
	    show_help
            shift ;;
        -l|--list)
            case "$2" in
                "") shift 2 ;;
                *) ARG_L=$2; shift 2 ;;
            esac ;;
        -m|--memory)
            case "$2" in
                "") shift 2 ;;
                *) ARG_M=$2; shift 2 ;;
            esac ;;
        -r|--reference)
            case "$2" in
                "") shift 2 ;;
                *) ARG_R=$2; shift 2 ;;
            esac ;;
        -t|--threads)
            case "$2" in
                "") shift 2 ;;
                *) ARG_T=$2; shift 2 ;;
            esac ;;
	--macse)
	    case "$2" in
	        "") shift 2 ;;
		*) PathMacse=$2; shift 2 ;;
	    esac ;;
	--sortdiamond)
	    case "$2" in
	        "") shift 2 ;;
		*) PathSortdiamond=$2; shift 2 ;;
	    esac ;;
	--splitfasta)
	    case "$2" in
	        "") shift 2 ;;
		*) PathSplitfsata=$2; shift 2 ;;
	    esac ;;
        --)
	   shift; break ;;
        *) echo "Unknown option: $1" 
	   exit 1 
	   ;;
    esac
done

fi

### Get and check some arguments

check_var() {
    local var_name="$1"
    local var_value="${!var_name}"  # get value

    if [ -z "$var_value" ]; then
        echo "Error: $var_name is not set or is empty"
	exit 1
    else
        echo "$var_name is set to: $var_value"

        case "$var_name" in
            "ARG_G")
                readarray -t genes < "$var_value"
                length_gn=${#genes[@]}
                ;;
            "ARG_L")
                readarray -t full_names < "$var_value"
                length_fn=${#full_names[@]}
                ;;
	    "ARG_F")
		check_command "cp"
		check_command "cd"
    		check_command "mv"
		check_command "find"
		check_command "mkdir"
		;;
        esac

    fi
}

check_path(){
    local path_name="$1"
    local path_value="${!path_name}"  # get value

    # expand ~
    path_value=$(eval echo "$path_value")

    if [ -e "$path_value" ]; then
        echo "$path_name exists at: $path_value"
    else
        echo "Error: $path_name does not exist at: $path_value"
	exit 1
    fi
}

check_command() {
    local cmd_name="$1"

    if command -v "$cmd_name" >/dev/null 2>&1; then
        echo "$cmd_name command exists."
    else
        echo "Error: $cmd_name command not found."
	exit 1
    fi
}


### Quality control && Trimming

if [ "$ARG_F" = "all" ] || [ "$ARG_F" = "clean" ]; then

	## Prepare
	mkdir -p $DirQcTrim

	check_var "ARG_L"
	check_command "fastp"

	    readarray -t full_names < "$ARG_L"
    length_fn=${#full_names[@]}

    readarray -t genes < "$ARG_G"
    length_gn=${#genes[@]}

	## Quality control and trimming using fastp
	for (( i=0; i<$length_fn; i++ )); do
		fastp -i $DirRaw/${full_names[$i]}_R1.fastq.gz -I $DirRaw/${full_names[$i]}_R2.fastq.gz -j $DirQcTrim/${full_names[$i]}.json -h $DirQcTrim/${full_names[$i]}.html -o $DirQcTrim/${full_names[$i]}_R1.fastq.gz -O $DirQcTrim/${full_names[$i]}_R2.fastq.gz -w $ARG_T
	done

fi

### De novo assembly

if [ "$ARG_F" = "all" ] || [ "$ARG_F" = "assembly" ]; then

	## Prepare
	mkdir -p $DirAssembly

	check_var "ARG_L"
	check_command "spades.py"

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
	check_var "ARG_C"
	check_var "ARG_L"

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
	check_var "ARG_C"
	check_var "ARG_R"
	check_command "diamond"

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
	#13:	qlen: Query sequence length
	#14:	slen: Subject sequence length
	#15:	gaps: Total number of gaps
	#16:	ppos: Percentage of positive - scoring matches
	#17:	qframe: Query frame (frames in blast?)
	#18:	qseq: Aligned part of query sequence

	done
	cd -

	mv $DirFasta/$ARG_C/*.m8 $DirMap

fi

if [ "$ARG_F" = "all" ] || [ "$ARG_F" = "pre" ]; then

	mkdir -p $DirPre

	check_var "ARG_L"
	check_path "PathSortdiamond"

	# extract fasta file from diamond balst style output
	for (( i=0; i<$length_fn; i++ )); do	
		$PathSortdiamond $DirMap/${full_names[$i]}.m8 $DirPre/${full_names[$i]}.fasta
	done
fi


if [ "$ARG_F" = "all" ] || [ "$ARG_F" = "split" ]; then
	mkdir -p $DirSplit
	cd $DirPre

	check_var "ARG_L"
	check_path "PathSplitfasta"

	# split fasta into different folders based on sequence names
	for (( i=0; i<$length_fn; i++ )); do
		$PathSplitfsata ${full_names[$i]}.fasta
	done

	# mv to destdir
	find . -mindepth 1 -maxdepth 1 -type d -exec mv {} ../$DirSplit \;
	cd -
fi

if [ "$ARG_F" = "all" ] || [ "$ARG_F" = "merge" ]; then

	## Check if the genes is specified
	check_var "ARG_G"

	mkdir -p $DirMerge
	cd $DirSplit

	# merge different taxa sequences in same gene to one fasta 
	for (( i=0; i<$length_gn; i++ )) 
	do
		cd ${genes[$i]}
		cat * > ../${genes[$i]}.fasta
		cd ..
	done

	# mv to destdir
	mv *.fasta ../$DirMerge

	cd -
	
fi

if [ "$ARG_F" = "all" ] || [ "$ARG_F" = "align" ]; then

	## Check if the genes is specified
	check_var "ARG_G"
	check_command "java"
	check_command "parallel"
	check_path "PathMacse"

	# current_thread=0
	mkdir -p $DirAlign
	mkdir -p $DirAlign/AA && mkdir -p $DirAlign/NT
	cd $DirMerge

	# align the sequence based on codon
	parallel -j $ARG_T java -jar $PathMacse -prog alignSequences -seq {}.fasta -out_AA ../$DirAlign/AA/{}.fasta -out_NT ../$DirAlign/NT/{}.fasta ::: "${genes[@]}"

	cd -

fi

