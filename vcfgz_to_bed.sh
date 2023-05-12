# submit as vcfgz_to_bed.sh <vcfgzfile> <name>
set -e

file=$(basename $1 ".vcf.gz")
name=$2

cd $(dirname $1)

source /Users/sreenivasan/opt/anaconda3/etc/profile.d/conda.sh
conda activate prepBED

#unzip
gunzip "${file}.vcf.gz"

# convert to bedfile and extract required columns for CADD-SV. Keep only DEL, DUP, and INS
SURVIVOR vcftobed ${file}.vcf 0 -1 tmp_${file}.bed
gzip ${file}.vcf

cut -f1,2,6,11 tmp_${file}.bed | grep -E 'DEL|INS|DUP' | grep -v -E '_random|_alt|chrUn|chrM'  > ${file}.bed && rm tmp_${file}.bed

# Sort bed file
bedtools sort -i ${file}.bed > ${file}.sorted.bed && rm ${file}.bed

#add name as the 5th column
awk -v name=$name 'BEGIN {OFS="\t"} {$5 = name; print}' ${file}.sorted.bed > ${file}.sorted.named.bed && rm ${file}.sorted.bed

# split into multiple, as needed
# split -l 9000 -d ${file}.sorted.bed "${file}.sorted_chunk"
# for tmp_file in ${file}.sorted_chunk*
# do
#    mv "$tmp_file" "${tmp_file}.bed"
# done

mkdir -p bedfiles_for_caddsv
mv *${file}*.bed bedfiles_for_caddsv/