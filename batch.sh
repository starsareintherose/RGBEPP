### Environment Setting

DirRaw=00_raw
DirQcTrim=01_fastp
DirAssembly=02_spades

### Get some arrays

cd $DirRaw

readarray -t full_names < <(ls | awk -F '_' '{print $1 "_" $2 "_" $3 "_" $4}' | uniq)
readarray -t species_names < <(ls | awk -F '_' '{print $2 "_" $3}' | uniq)
readarray -t output_names < <(ls | awk -F '_' '{print $2 "_" $3 "_" $4}' | uniq)

cd ..

length_fn=${#full_names[@]}
length_sn=${#species_names[@]}
length_on=${#output_names[@]}

### Check the arrays

if [ $length_fn -ne $length_sn ] || [ $length_fn -ne $length_on ] || [ $length_sn -ne $length_on ]
then
  echo "Please check the amount number of arrays"
  exit 0
fi

### Quality control && Trimming

mkdir -p $DirQcTrim

for (( i=0; i<$length_fn; i++ )); do
	fastp -i $DirRaw/${full_names[$i]}_R1.fastq.gz -I $DirRaw/${full_names[$i]}_R2.fastq.gz -j ${species_names[$i]}.json -h ${species_names[$i]}.html -o $DirQcTrim/${output_names[$i]}_R1.fastq.gz -O $DirQcTrim/${output_names[$i]}_R2.fastq.gz -w 4
done

### De novo assembly

mkdir -p $DirAssembly

for (( i=0; i<$length_fn; i++ )); do
	mkdir -p $DirAssembly/${species_names[$i]} 
	spades.py --pe1-1 $DirQcTrim/${output_names[$i]}_R1.fastq.gz --pe1-2 $DirQcTrim/${output_names[$i]}_R2.fastq.gz -t 8 -k 97,107,117,127 -m 14 --careful --phred-offset 33 -o $DirAssembly/${species_names[$i]}
done
