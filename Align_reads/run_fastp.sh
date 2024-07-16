#! /bin/bash

BASE_PATH='NGS'
FASTQ=${BASE_PATH}'/fastq'
TRIMMED=${BASE_PATH}'/fastp'

for INPUT_FORWARD in ${FASTQ}/*_R1_*.fastq.gz; do
	#echo forward $INPUT_FORWARD
	INPUT_REVERSE=${INPUT_FORWARD/R1_001.fastq.gz/R2_001.fastq.gz}
	#echo reverse $INPUT_REVERSE
	N=$(echo ${INPUT_FORWARD} | cut -d '/' -f 11)
	#echo $N
	OUTPUT_FORWARD=${N/_001.fastq.gz/_trimmed.fastq.gz}
	OUTPUT_REVERSE=${N/R1_001.fastq.gz/R2_trimmed.fastq.gz}
	fastp -i $INPUT_FORWARD -I $INPUT_REVERSE -o $TRIMMED/$OUTPUT_FORWARD -O $TRIMMED/$OUTPUT_REVERSE --dont_eval_duplication -3 --disable_adapter_trimming -h $TRIMMED/${N/R1_001.fastq.gz/_fastp.html} --umi --umi_loc read1 --umi_len 15 --umi_prefix UMI --thread 32
done