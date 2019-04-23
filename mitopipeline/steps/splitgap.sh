#!/bin/bash
# $1: filename
# $2: start directory
# $3: out directory
# $4: tools directory
# $5: setup directory
BAMS=$2
#last string following / delimeter will be name of the previous job
#only set filetype if the value was a step in the pipeline
filetype="_"$(awk -F/ '{print $NF}' <<< "$2" | awk '{print tolower($0)}')
if [ "$filetype" != "_extractmito" ] && [ "$filetype" != "_clipping" ];
then
filetype=""
fi

#old fastq file used by split_bams.py
samtools fastq $BAMS/$1${filetype}.bam -n > $3/$1_old.fastq

samtools view $BAMS/$1${filetype}.bam > $3/$1_bam.txt

#1 is bam file id, 2 is new fastq file to write to
python $5/split_bams.py $3 $1 

#remove bam txt + old fastq file
rm $3/$1_bam.txt
rm $3/$1_old.fastq
mv $3/$1_splitgap.fastq $3/$1_1_splitgap.fastq