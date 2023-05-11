# submit as vcfgz_to_bed.sh <vcfgzfile>
set -e

file=$(basename $1 ".vcf.gz")
cd $(dirname $1)

source /Users/sreenivasan/opt/anaconda3/etc/profile.d/conda.sh
conda activate prepBED

#unzip
gunzip "${file}.vcf.gz"

# convert to bedfile and extract required columns for CADD-SV. Keep only DEL, DUP, and INS
SURVIVOR vcftobed ${file}.vcf 0 -1 tmp_${file}.bed
cut -f1,2,6,11 tmp_${file}.bed | grep -E 'DEL|INS|DUP' | grep -v -E '_random|_alt|chrUn|chrM'  > ${file}.bed && rm tmp_${file}.bed

# split into multiple, as needed
split -l 9000 -d ${file}.bed "${file}_chunk"
for tmp_file in ${file}_chunk*
do
    mv "$tmp_file" "${tmp_file}.bed"
done

mkdir -p bedfiles
mv *${file}*.bed bedfiles/