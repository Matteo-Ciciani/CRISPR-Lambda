#! /bin/bash
BASE_PATH='NGS'
TRIMMED=${BASE_PATH}'/fastp'
MAPPING=${BASE_PATH}'/bwa'

cd ${BASE_PATH}'/amplicons'

for INPUT_FORWARD in $TRIMMED/SPIKE*_R1_trimmed.fastq.gz; do
    INPUT_REVERSE=${INPUT_FORWARD/_R1_trimmed.fastq.gz/_R2_trimmed.fastq.gz}
    #echo $INPUT_REVERSE
    N=$(echo ${INPUT_FORWARD} | cut -d '/' -f 11 | cut -d '.' -f 1 )
    #echo ${N/_L001_R1_trimmed/}.sam
    bwa mem -t 32 spike.fasta $INPUT_FORWARD $INPUT_REVERSE > $MAPPING/${N/_R1_trimmed/}.sam
done
