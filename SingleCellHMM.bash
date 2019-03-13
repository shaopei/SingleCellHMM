# bash SingleCellHMM.bash  Path_to_bam_file/pbmc4k_possorted_genome_bam.bam Path_to_SingleCellHMM.R

INPUT_BAM=$1 #pbmc4k_possorted_genome_bam.bam
PL=$2
${PL:=/workdir/sc2457/SingleCellHMM/}

PREFIX=`echo ${INPUT_BAM} | rev | cut -d / -f 1 |cut -d . -f 2- |rev`
tmp=`cat /dev/urandom | tr -dc 'a-zA-Z0-9' | fold -w 32 | head -n 1`
TMPDIR=${PREFIX}_${tmp}
mkdir ${TMPDIR}

exec > >(tee SingleCellHMM_Run_${TMPDIR}.log)
exec 2>&1
echo "Path to SingleCellHMM.R   $PL" 
echo "INPUT_BAM                 $INPUT_BAM"
echo "temp folder               $TMPDIR"
echo ""
echo "Reads spanning over splicing junction will join HMM blocks"
echo "To avoid that, split reads into small blocks before input to groHMM"
echo "Spliting and sorting reads..."
bedtools bamtobed -i ${INPUT_BAM} -split |LC_ALL=C sort -k1,1V -k2,2n --parallel=30| awk '{print "chr"$0}' | gzip > ${TMPDIR}/${PREFIX}_split.sorted.bed.gz 

cd ${TMPDIR}
zcat ${PREFIX}_split.sorted.bed.gz  |awk '{print $0 >> $1".bed"}' 
find -name "*.bed" -size -1024k -delete
#wc chr*.bed -l > chr_read_count.txt

echo ""
echo "Start to run groHMM in each individual chromosome..."


wait_a_second() {
	joblist=($(jobs -p))
    while (( ${#joblist[*]} >= 50 ))
	    do
	    sleep 1
	    joblist=($(jobs -p))
	done
}


for f in chr*.bed
do 
wait_a_second
R --vanilla --slave --args $(pwd) ${f}  < ${PL}/SingleCellHMM.R  > ${f}.log 2>&1 &
done
#R --vanilla --slave --args $(pwd) ${PREFIX}_split.sorted.bed.gz  < ${PL}/SingleCellHMM.R 
wait


echo ""
echo "Merging HMM blocks within 500bp..."
for f in chr*_HMM.bed
do	
  LC_ALL=C sort -k1,1V -k2,2n --parallel=30 ${f} > ${f}.sorted.bed
  cat ${f}.sorted.bed | grep + > ${f}_plus
  cat ${f}.sorted.bed | grep - > ${f}_minus
  bedtools merge -s -d 500 -i ${f}_plus > ${f}_plus_merge500 &
  bedtools merge -s -d 500 -i ${f}_minus > ${f}_minus_merge500 &
  wait_a_second
done

wait 
mkdir toremove
for f in chr*_HMM.bed
do	
mv ${f}.sorted.bed ${f}_plus ${f}_minus ${f}_merge500 toremove/.
done

#f=${PREFIX}_split.sorted_HMM
#gzip ${f}.bed &

cat chr*_HMM.bed_plus_merge500 | awk 'BEGIN{OFS="\t"} {print $0, ".", ".", "+"}' > ${PREFIX}_merge500
cat chr*_HMM.bed_minus_merge500 | awk 'BEGIN{OFS="\t"} {print $0, ".", ".", "-"}' >> ${PREFIX}_merge500

for f in chr*_HMM.bed
do	
mv ${f}.sorted.bed ${f}_plus ${f}_minus ${f}_merge500 ${f}_plus_merge500 ${f}_minus_merge500 toremove/.
done


echo ""
echo "Calculating the coverage..." 
f=${PREFIX}
LC_ALL=C sort -k1,1V -k2,2n ${f}_merge500 --parallel=30 > ${f}_merge500.sorted.bed
rm ${f}_merge500


bedtools coverage -a ${f}_merge500.sorted.bed -b <(zcat ${PREFIX}_split.sorted.bed.gz |awk '{print "chr"$0}') -s -counts -split -sorted > ${f}_merge500.sorted.bed_count

echo ""
echo "Filtering the HMM blocks by coverage..." 
cat ${f}_merge500.sorted.bed_count | awk 'BEGIN{OFS="\t"} ($7 >= 2){print $1, $2, $3, $4, $5, $6}' | gzip > ${f}_merge500_2reads.bed.gz
cat ${f}_merge500.sorted.bed_count | awk 'BEGIN{OFS="\t"} ($7 >= 5){print $1, $2, $3, $4, $5, $6}' | gzip > ${f}_merge500_5reads.bed.gz

#rm ${PREFIX}_split.bed.gz
#mv ${f}_merge500.sorted.bed ${f}_merge500.sorted.bed_count ${TMPDIR}/.
cd ..
ln -s ${TMPDIR}/${f}_merge500_5reads.bed.gz .

cd ${TMPDIR}
for f in *
do gzip ${f} &
done



