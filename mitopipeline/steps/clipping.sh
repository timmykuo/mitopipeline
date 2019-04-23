#!/bin/bash
### $1 IS THE FILENAME
### $2 IS DIRECTORY FOR LOCATION OF FILE (START)
### $3 IS WHERE TO STORE THE FILE (OUT)
### $4 is tools directory
### $5 is steps directory
### $6 is refs directory
START=$2
OUT=$3
TOOLS=$4
filename=$1
#last string following / delimeter will be name of the previous job
filetype="_"$(awk -F/ '{print $NF}' <<< "$2" | awk '{print tolower($0)}')
#if this is the first step in the pipeline
if [ "$filetype" != "_extractmito" ];
then
filetype=""
fi
echo .
echo .
echo .
echo ----bam2fastq
java -jar $5/tools/SamToFastq.jar \
I=$START/$1$filetype.bam \
FASTQ=$OUT/$1_bam2fastq.fastq \
SECOND_END_FASTQ=$OUT/$1_bam2fastq_2.fastq

echo ****BAM2FASTQ DONE.
echo .
echo .
echo .
echo .
echo .
echo .
echo ----Clipping Initiated

if [ -e "$OUT/$1_bam2fastq_2.fastq" ]
then
echo "PAIRED-END"
        echo "--CLIPPED: Removing first and last 2 base pairs from every read"
        $TOOLS/seqtk/seqtk trimfq -b 2 -e 2 $OUT/$1_bam2fastq.fastq > $OUT/$1_1_clipping.fastq
        $TOOLS/seqtk/seqtk trimfq -b 2 -e 2 $OUT/$1_bam2fastq_2.fastq > $OUT/$1_2_clipping.fastq
        rm $OUT/$1_bam2fastq_2.fastq
else
echo "SINGLE-END"
        echo "--CLIPPED: Removing first and last 2 base pairs from every read"
        $TOOLS/seqtk/seqtk trimfq -b 2 -e 2 $OUT/$1_bam2fastq.fastq > $OUT/$1_1_clipping.fastq
fi
rm $OUT/$1_bam2fastq.fastq
echo ****Clipping DONE.
echo .
echo .
