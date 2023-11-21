#!/bin/bash


### define variables ###
DIR="../../data/subassembly/PEfl"

if [[ -z "$TMPDIR" ]]; then
  TMPDIR=${DIR}/tmp${RANDOM}
  mkdir -p ${TMPDIR}
fi


### extract barcode count files ###
for FILE in ${DIR}/*.count.gz
do
  FILENAME=`basename ${FILE}`
  unpigz -c $FILE > ${TMPDIR}/${FILENAME%.gz}
done


### extract barcodes with at least 5 reads ###
awk -v OFS='\t' '
    BEGIN {print "barcode", "enhancer", "part", "orientation"}
    FNR == 1 {
      sub(/.*\//, "", FILENAME) 
      split(FILENAME, file, /[_.]/)
    }
    $1 >= 5 {
        print $2, file[3], file[4], file[5]
  }' ${TMPDIR}/*.count \
  >> ${TMPDIR}/barcodes_pPSup_PEfl.tsv


### compress final files and copy to main directory ###
pigz -cf ${TMPDIR}/barcodes_pPSup_PEfl.tsv > ${DIR}/barcodes_pPSup_PEfl.tsv.gz

