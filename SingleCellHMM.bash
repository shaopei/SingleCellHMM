# bash SingleCellHMM.bash  pbmc4k_possorted_genome_bam.bam

INPUT_BAM=$1 #pbmc4k_possorted_genome_bam.bam
PREFIX=`echo ${INPUT_BAM} | rev | cut -d . -f 2- |rev`
TMPDIR=`cat /dev/urandom | tr -dc 'a-zA-Z0-9' | fold -w 32 | head -n 1`
mkdir ${TMPDIR}

exec > >(tee SingleCellHMM_Run_${TMPDIR}.log)
exec 2>&1

echo "INPUT_BAM                 $INPUT_BAM"
echo "temp folder               $TMPDIR"
echo ""
echo "Reads spanning over splicing junction will join HMM blocks"
echo "To avoid that, split reads into small blocks before input to groHMM"
echo "Spliting reads..."
bedtools bamtobed -i ${INPUT_BAM} -split |gzip > ${PREFIX}_split.bed.gz

echo ""
echo "Sorting reads in the background..."
zcat ${PREFIX}_split.bed.gz |LC_ALL=C sort -k1,1V -k2,2n --parallel=30 | gzip > ${PREFIX}_split.sorted.bed.gz &
 
echo ""
echo "Start to run groHMM..."
R --vanilla --slave --args $(pwd) ${PREFIX}_split.bed.gz < SingleCellHMM.R 




echo ""
echo "Merging HMM blocks within 500bp..."
f=${PREFIX}_split_HMM
cat $f.bed | grep + > ${f}_plus
cat $f.bed | grep - > ${f}_minus
bedtools merge -s -d 500 -i ${f}_plus > ${f}_plus_merge500
bedtools merge -s -d 500 -i ${f}_minus > ${f}_minus_merge500
rm ${f}_plus ${f}_minus
gzip ${PREFIX}_split_HMM.bed

cat ${f}_plus_merge500 | awk 'BEGIN{OFS="\t"} {print $0, ".", ".", "+"}' > ${f}_merge500
cat ${f}_minus_merge500 | awk 'BEGIN{OFS="\t"} {print $0, ".", ".", "-"}' >> ${f}_merge500
rm ${f}_plus_merge500 ${f}_minus_merge500

echo ""
echo "Calculating the coverage..." 
LC_ALL=C sort -k1,1V -k2,2n ${f}_merge500 --parallel=30 > ${f}_merge500.sorted.bed
rm ${f}_merge500

wait
bedtools coverage -a ${f}_merge500.sorted.bed -b ${PREFIX}_split.sorted.bed.gz -s -counts -split -sorted > ${f}_merge500.sorted.bed_count

echo ""
echo "Filtering the HMM blocks by coverage..." 
cat ${f}_merge500.sorted.bed_count | awk 'BEGIN{OFS="\t"} ($7 >= 2){print $1, $2, $3, $4, $5, $6}' | gzip > ${f}_merge500_2reads.bed.gz
cat ${f}_merge500.sorted.bed_count | awk 'BEGIN{OFS="\t"} ($7 >= 5){print $1, $2, $3, $4, $5, $6}' | gzip > ${f}_merge500_5reads.bed.gz

rm ${PREFIX}_split.bed.gz
mv ${f}_merge500.sorted.bed ${f}_merge500.sorted.bed_count ${TMPDIR}/.
gzip ${TMPDIR}/*
