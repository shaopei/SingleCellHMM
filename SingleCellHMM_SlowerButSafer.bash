# bash SingleCellHMM.bash  Path_to_bam_file/PREFIX.bam numberOfThread 5 500 Path_to_SingleCellHMM.R

INPUT_BAM=$1 #pbmc4k_possorted_genome_bam.bam
CORE=$2
MINCOV=$3
MERGEBP=$4
PL=$5

CORE=${CORE:=1}
MINCOV=${MINCOV:=5}
MERGEBP=${MERGEBP:=500}
PL=${PL:=~}

PREFIX=`echo ${INPUT_BAM} | rev | cut -d / -f 1 |cut -d . -f 2- |rev`
tmp="HMM_features"
TMPDIR=${PREFIX}-${tmp}


exec > >(tee SingleCellHMM_Run_${TMPDIR}.log)
exec 2>&1
date 
echo "Path to SingleCellHMM.R   $PL" 
echo "INPUT_BAM                 $INPUT_BAM"
if [[ $PREFIX == *"-"* ]]; then
  echo "Please use NO - in INPUT_BAM"
  exit 0
fi
echo "temp folder               $TMPDIR"
echo "number Of thread          $CORE"
echo "minimum coverage		      $MINCOV"
echo ""
echo "Reads spanning over splicing junction will join HMM blocks"
echo "To avoid that, split reads into small blocks before input to groHMM"
echo "Spliting and sorting reads..."


wait_a_second() {
  joblist=($(jobs -p))
    while (( ${#joblist[*]} >= ${CORE} ))
      do
      sleep 1
      joblist=($(jobs -p))
  done
}

mkdir ${TMPDIR}
bedtools bamtobed -i ${INPUT_BAM} -split |awk -v p=${PREFIX} -v d=${TMPDIR} 'BEGIN {OFS=""}(substr($1,1,3)=="chr"){print $0 >> d"/"p"-"$1".bed"} (substr($1,1,3)!="chr") {print "chr"$0 >> d"/"p"-chr"$1".bed"}' 

cd ${TMPDIR}
echo "Delete chr bed files that are smaller than 1G ..."
find -name "${PREFIX}-chr*.bed" -size -1024k -delete

echo "Sort each chr bed file that is larger than 1G ..."
for f in ${PREFIX}-chr*.bed
do echo $f
h=`echo $f |rev |cut -d . -f 2-|rev`
echo $h
#sort-bed $f | gzip > ${h}.sorted.bed.gz &
LC_ALL=C sort --parallel=${CORE} -k1,2n -k2,3n $f > ${h}.sorted.bed
#gzip ${h}.sorted.bed &
done

wait

#wc chr*.bed -l > chr_read_count.txt
echo ""
echo "Make ${PREFIX}_split.sorted.bed.gz... from individual chromosomes"
rm ${PREFIX}_split.sorted.bed.gz
for c in `ls ${PREFIX}-chr*.sorted.bed |rev|cut -d . -f 3-| rev | cut -d - -f 2-|LC_ALL=C sort -V`
 do echo $c
 cat ${PREFIX}-${c}.sorted.bed |gzip >>  ${PREFIX}_split.sorted.bed.gz
done

echo ""
echo "Start to run groHMM in each individual chromosome..."

for f in ${PREFIX}-chr*.sorted.bed
do 
R --vanilla --slave --args $(pwd) ${f}  < ${PL}/SingleCellHMM.R  > ${f}.log 2>&1 
done
#R --vanilla --slave --args $(pwd) ${PREFIX}_split.sorted.bed.gz  < ${PL}/SingleCellHMM.R 
wait



echo ""
echo "Merging HMM blocks within ${MERGEBP}bp..."
for f in ${PREFIX}-chr*_HMM.bed
do	
  LC_ALL=C sort -k1,1V -k2,2n --parallel=${CORE} ${f} > ${f}.sorted.bed
  cat ${f}.sorted.bed | grep + > ${f}_plus
  cat ${f}.sorted.bed | grep - > ${f}_minus
  bedtools merge -s -d ${MERGEBP} -i ${f}_plus > ${f}_plus_merge${MERGEBP} &
  bedtools merge -s -d ${MERGEBP} -i ${f}_minus > ${f}_minus_merge${MERGEBP} &
  wait_a_second
done

wait 


#f=${PREFIX}_split.sorted_HMM
#gzip ${f}.bed &

cat ${PREFIX}-chr*_HMM.bed_plus_merge${MERGEBP} | awk 'BEGIN{OFS="\t"} {print $0, ".", ".", "+"}' > ${PREFIX}_merge${MERGEBP}
cat ${PREFIX}-chr*_HMM.bed_minus_merge${MERGEBP} | awk 'BEGIN{OFS="\t"} {print $0, ".", ".", "-"}' >> ${PREFIX}_merge${MERGEBP}

mkdir toremove
for f in ${PREFIX}-chr*_HMM.bed
do	
mv ${f}.sorted.bed ${f}_plus ${f}_minus ${f}_plus_merge${MERGEBP} ${f}_minus_merge${MERGEBP} toremove/.
done


echo ""
echo "Calculating the coverage..." 
LC_ALL=C sort -k1,1V -k2,2n ${PREFIX}_merge${MERGEBP} --parallel=${CORE} > ${PREFIX}_merge${MERGEBP}.sorted.bed
rm ${PREFIX}_merge${MERGEBP}


bedtools coverage -a ${PREFIX}_merge${MERGEBP}.sorted.bed -b <(zcat ${PREFIX}_split.sorted.bed.gz) -s -counts -split -sorted > ${PREFIX}_merge${MERGEBP}.sorted.bed_count

echo ""
echo "Filtering the HMM blocks by coverage..." 
#cat ${PREFIX}_merge${MERGEBP}.sorted.bed_count | awk 'BEGIN{OFS="\t"} ($7 >= 2){print $1, $2, $3, $4, $5, $6}' | gzip > ${PREFIX}_merge${MERGEBP}_2reads.bed.gz
cat ${PREFIX}_merge${MERGEBP}.sorted.bed_count | awk 'BEGIN{OFS="\t"} ($7 >= '$MINCOV'){print $1, $2, $3, $4, $5, $6, $7}' | gzip > ${PREFIX}_merge${MERGEBP}_${MINCOV}reads.bed.gz

echo "" 
echo "#### Please examine if major chromosomes are all present in the final PREFIX_merge${MERGEBP}_${MINCOV}reads.bed.gz file ####"
zcat ${PREFIX}_merge${MERGEBP}_${MINCOV}reads.bed.gz |cut -f 1 |uniq

echo "" 
echo "Link the final PREFIX_merge${MERGEBP}_${MINCOV}reads.bed.gz file to the working directory"
cd ..
ln -s ${TMPDIR}/${PREFIX}_merge${MERGEBP}_${MINCOV}reads.bed.gz .


echo ""
echo "Move intermediate files to  ${TMPDIR}/toremove ..." 
echo ""
echo "${TMPDIR}/toremove can be deleted if no error message in SingleCellHMM_Run log file and "
echo "all major chromosomes are present in the final PREFIX_merge${MERGEBP}_${MINCOV}reads.bed.gz file"
echo ""
echo ""

cd ${TMPDIR}
mv ${PREFIX}-chr* toremove/.

for f in *
do gzip ${f} &
done

cd toremove
for f in *
do gzip ${f} &
done

cd ../..

echo ""
date 
echo "done!"

