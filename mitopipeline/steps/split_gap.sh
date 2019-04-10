#!/bin/bash
# $1: bam file id
# $2: split fastq file


# $1: filename
# $2: start directory
# $3: out directory

module load intel/17  
module load openmpi/2.0.1
module load samtools

BAMS=$2

#old fastq file used by split_bams.py
samtools fastq $BAMS/$1_mito.bam -n > $1_old.fastq

samtools view $BAMS/$1_mito.bam > $1_bam.txt

#1 is bam file id, 2 is new fastq file to write to
python ./steps/split_bams.py $1 

#remove bam txt + old fastq file
rm $1_bam.txt
rm $1_old.fastq
mv $1.fastq $3
#mv to /mnt/rds/txl80/LaframboiseLab/tyk3/TARGET_PIPELINE/fastqs
