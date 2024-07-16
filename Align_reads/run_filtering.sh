#! /bin/bash

BASE_PATH='NGS'
MAPPING=${BASE_PATH}'/bwa'

## RUN SAMTOOLS
cd $MAPPING
for bam in *sorted.bam; do
	samtools view --threads 32 -q 50 -F 3844 -f 2 -b $bam > filtered/${bam/sorted.bam/filtered.bam}; samtools index filtered/${bam/sorted.bam/filtered.bam}
done