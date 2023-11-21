#!/bin/bash

# usage:
# extract_barcode.sh <basename of reads>

# requires:
# - Python (version 3.9.13)
# - cutadapt (version 4.4)
# - PANDAseq (version 2.11)
# - bioawk (version fd40150)


### define variables ###
NSLOTS=${NSLOTS:-1}

SAMPLE=${1}

if [[ ${SAMPLE} == *"/"* ]]; then
  SUBDIR=`dirname ${SAMPLE}`
  SUBDIR="/"${SUBDIR}
  SAMPLE=`basename ${SAMPLE}`
else
  SUBDIR=""
fi

READDIR="../../reads/${SUBDIR}"
OUTDIR="../../data/subassembly${SUBDIR}"
mkdir -p ${OUTDIR}

if [[ -z "$TMPDIR" ]]; then
  TMPDIR=${OUTDIR}/tmp${RANDOM}
  mkdir -p ${TMPDIR}
fi


### print temporary directory name ###
echo "read files are at: ${READDIR}"
echo "final files will be deposited at: ${OUTDIR}"
echo "temporary directory is at: ${TMPDIR}"
echo ""


### trim insert reads ###
echo "trimming insert reads ($(date +"%T"))"
cutadapt --cores ${NSLOTS} --report minimal --nextseq-trim 20 --minimum-length 5 -a AGGAGCTGGCAAGACCCTTC -A GAGTGGAATTCGATATCAAG --output ${TMPDIR}/${SAMPLE}_insert_1.fastq --paired-output ${TMPDIR}/${SAMPLE}_insert_2.fastq \
  ${READDIR}/${SAMPLE}_insert_1.fastq.gz \
  ${READDIR}/${SAMPLE}_insert_2.fastq.gz

echo ""


### assemble inserts and barcodes ###
echo "assembling reads ($(date +"%T"))"

mkdir -p ${TMPDIR}/panda_logs

if [[ ${SAMPLE} == *rev ]]; then
  R1=2
  R2=1
else
  R1=1
  R2=2
fi

pandaseq -N -f ${TMPDIR}/${SAMPLE}_insert_${R1}.fastq -r ${TMPDIR}/${SAMPLE}_insert_${R2}.fastq -d bfsrkm -G ${TMPDIR}/panda_logs/log_${SAMPLE}_insert.txt.bz2 -w ${TMPDIR}/${SAMPLE}_insert.fa
pandaseq -N -f ${READDIR}/${SAMPLE}_barcode_1.fastq.gz -r ${READDIR}/${SAMPLE}_barcode_2.fastq.gz -d bfsrkm -G ${TMPDIR}/panda_logs/log_${SAMPLE}_barcode.txt.bz2 -w ${TMPDIR}/${SAMPLE}_barcode.fa


### join barcodes and inserts ###
echo "joining barcodes and inserts ($(date +"%T"))"

# filter for barcodes of correct length
bioawk -c fastx -t '{
    sub(/:[ACGTN0+]*;.*$/, "", $name)
    if (length($seq) == 15) print $name, $seq
  }' ${TMPDIR}/${SAMPLE}_barcode.fa \
  > ${TMPDIR}/${SAMPLE}_barcode.tsv

# join barcodes and inserts; keep only combinations seen at least twice
bioawk -c fastx -t -v tmpdir=${TMPDIR} -v sample=${SAMPLE} '
    BEGIN{
      while (getline < (tmpdir "/" sample "_barcode.tsv")) {
        split($0, line, "\t")
        barcode[line[1]] = line[2]
      }
      print "barcode", "sequence", "assembly_count"
    } {
      sub(/:[ACGTN0+]*;.*$/, "", $name);
      if ($name in barcode) {counts[barcode[$name]$seq]++; delete barcode[$name]}
    } END{
      for (combo in counts) if (counts[combo] > 1) print substr(combo, 1,  15), substr(combo, 16), counts[combo]
    }
  ' ${TMPDIR}/${SAMPLE}_insert.fa \
  > ${TMPDIR}/raw_subassembly_${SAMPLE}.tsv


### compress final files and copy to main directory ###
echo "copying files ($(date +"%T"))"

for FILE in ${TMPDIR}/raw_subassembly_*.tsv
do
  FILENAME=`basename ${FILE}`
  pigz -c ${FILE} > ${OUTDIR}/${FILENAME}.gz
done

echo "all done ($(date +"%T"))"
