#!/bin/bash
### $4 is the file size
### $3 is the Cancer Type
### $2 is the Target file ID
### $1 is the download ID

### $1 IS THE FILENAME
### $2 IS DIRECTORY FOR LOCATION OF FILE (START)
### $3 IS WHERE TO STORE THE FILE (OUT)

START=$2
OUT=$3
module load samtools
module load bwa

filename=$1
echo $filename
cancertype=$(awk -F_ '{print $1}' <<< "$filename")
id=$(awk -F_ '{print $2"_"$3}' <<< "$filename")

echo .
echo .
echo .
echo ----bam2fastq
#samtools fastq doesn't handle paired end, update this to allow for paired end fastq format
samtools fastq -n $START/$1_mito.bam > $OUT/$1_mito.fastq
echo ****BAM2FASTQ DONE.
echo .
echo .
echo .
echo .
echo .
echo .
echo ----Clipping Initiated

if [ -e "$OUT/${cancertype}_${id}_mito_1.fastq" ]
then
echo "PAIRED-END"
        echo "--CLIPPED: Removing first and last 2 base pairs from every read"
        ./tools/seqtk-master/seqtk trimfq -b 2 -e 2 $OUT/$1_mito_1.fastq > $OUT/$1_mito_CLIPPED_1.fastq
        ./tools/seqtk-master/seqtk trimfq -b 2 -e 2 $OUT/$1_mito_2.fastq > $OUT/$1_${id}_mito_CLIPPED_2.fastq
else
echo "SINGLE-END"
        echo "--CLIPPED: Removing first and last 2 base pairs from every read"
        ./tools/seqtk-master/seqtk trimfq -b 2 -e 2 $OUT/$1_mito.fastq > $OUT/$1_mito_CLIPPED.fastq
fi
echo ****Clipping DONE.
echo .
echo .
