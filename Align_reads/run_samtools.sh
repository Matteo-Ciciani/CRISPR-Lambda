#! /bin/bash

BASE_PATH='NGS'
MAPPING=${BASE_PATH}'/bwa'

## RUN SAMTOOLS
cd $MAPPING
for sam in *sam; do
	samtools view --threads 32 -bS $sam > ${sam/sam/bam}; samtools sort --threads 32 -o ${sam/sam/sorted.bam} ${sam/sam/bam}; samtools index ${sam/sam/sorted.bam}
done