#!/bin/bash

SAMPLE_LIST="sample_names.txt"
CORES=12
BED="CancerPanel_virtuWES_3bpPad.bed"
REF="hg38.bwa.fa"

TO_RUN=$(sed 's%\(\S\+\)%-t bam/\1.bam%g' ${SAMPLE_LIST} | tr "\n" " ")
#echo ${TO_RUN}
#exit
./svaba/bin/svaba run ${TO_RUN} \
			-p ${CORES} \
			-k ${BED} \
			-L 6 \
			-I \
			-a cgp_sv_calls \
			-G ${REF}
