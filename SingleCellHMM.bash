# bash SingleCellHMM.bash  05-007B1_gene_exon_tagged.REF_chr20.bam

INPUT_BAM=$1 #05-007B1_gene_exon_tagged.REF_chr22.bam
PREFIX=`echo ${INPUT_BAM} | rev | cut -d . -f 2- |rev`
TMPDIR=`cat /dev/urandom | tr -dc 'a-zA-Z0-9' | fold -w 32 | head -n 1`
mkdir ${TMPDIR}

exec > >(tee SingleCellHMM_Run_${tmp}.log)
exec 2>&1

echo "INPUT_BAM                 $INPUT_BAM"
echo "temp folder               $TMPDIR"
echo ""
echo "Reads  spread over splicing junction will join HMM blocks"
echo "To avoid that, splict reads into small blocks before input to groHMM"
echo "Splicting reads..."
bedtools bamtobed -i ${INPUT_BAM} -split |LC_ALL=C sort -k1,1V -k2,2n --parallel=30 |gzip > ${PREFIX}_split.bed.gz


echo "Start to run groHMM..."
R --vanilla --slave --args $(pwd) ${PREFIX}_split.bed.gz < SingleCellHMM.R 
gzip ${PREFIX}_split_HMM.bed




echo "Mergeing HMM blocks within 500bp..."
f=${PREFIX}_split_HMM
zcat $f.bed.gz | grep + > ${f}_plus
zcat $f.bed.gz | grep - > ${f}_minus
bedtools merge -s -d 500 -i ${f}_plus > ${f}_plus_merge500
bedtools merge -s -d 500 -i ${f}_minus > ${f}_minus_merge500
rm ${f}_plus ${f}_minus

cat ${f}_plus_merge500 | awk 'BEGIN{OFS="\t"} {print $0, ".", ".", "+"}' > ${f}_merge500
cat ${f}_minus_merge500 | awk 'BEGIN{OFS="\t"} {print $0, ".", ".", "-"}' >> ${f}_merge500
rm ${f}_plus_merge500 ${f}_minus_merge500

echo "Calcuating the coverage..." 
LC_ALL=C sort -k1,1V -k2,2n ${f}_merge500 --parallel=30 > ${f}_merge500.sorted.bed
rm ${f}_merge500

bedtools coverage -a ${f}_merge500.sorted.bed -b ${PREFIX}_split.bed.gz -s -counts -split -sorted > ${f}_merge500.sorted.bed_count

echo "Filtering the HMM blocks by coverage..." 
cat ${f}_merge500.sorted.bed_count | awk 'BEGIN{OFS="\t"} ($7 >= 2){print $1, $2, $3, $4, $5, $6}' | gzip > ${f}_merge500_2reads.bed.gz
cat ${f}_merge500.sorted.bed_count | awk 'BEGIN{OFS="\t"} ($7 >= 5){print $1, $2, $3, $4, $5, $6}' | gzip > ${f}_merge500_5reads.bed.gz

mv ${f}_merge500.sorted.bed ${f}_merge500.sorted.bed_count ${TMPDIR}/.
gzip ${TMPDIR}/*
