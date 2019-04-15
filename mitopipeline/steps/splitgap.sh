#!/bin/bash
# $1: filename
# $2: start directory
# $3: out directory
# $4: tools directory
# $5: setup directory
module load intel/17  
module load openmpi/2.0.1
module load samtools

BAMS=$2

#old fastq file used by split_bams.py
samtools fastq $BAMS/$1.bam -n > $3/$1_old.fastq

samtools view $BAMS/$1.bam > $3/$1_bam.txt

#1 is bam file id, 2 is new fastq file to write to
python $5/split_bams.py $3 $1 

#remove bam txt + old fastq file
rm $3/$1_bam.txt
rm $3/$1_old.fastq
#mv to /mnt/rds/txl80/LaframboiseLab/tyk3/TARGET_PIPELINE/fastqs
