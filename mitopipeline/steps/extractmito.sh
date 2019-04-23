#!/bin/bash
### $1 is file name
### $2 is directory to start at
### $3 is the out directory
### $4 is tools directory
### $5 is steps directory
### $6 is refs directory
START=$2
OUT=$3
TOOLS=$4
FILENAME=$1

echo extracting chrM bam from $FILENAME
        samtools index $START/$FILENAME.bam
        samtools view -H $START/$FILENAME.bam > $START/$FILENAME.header.bam
        grep -i SN:MT $START/$FILENAME.header.bam > $START/$FILENAME.MT.bam
        grep -i SN:chrM_rCRS $START/$FILENAME.header.bam > $START/$FILENAME.chrM_rCRS.bam
        grep -i SN:chrM $START/$FILENAME.header.bam > $START/$FILENAME.chrM.bam
        grep -i SN:M $START/$FILENAME.header.bam > $START/$FILENAME.M.bam
        
        rm $START/$FILENAME.header.bam
        if [ -s "$START/$FILENAME.MT.bam" ]
        then
                echo "The mitochondrial genome is annotated by 'MT'."
                samtools view -b $START/$FILENAME.bam MT > $OUT/${FILENAME}_extractmito.bam
        elif [ -s "$START/$FILENAME.chrM_rCRS.bam" ]
        then
                echo "The mitochondrial genome is annotated by 'chrM_rCRS'."
                samtools view -b $START/$FILENAME.bam chrM_rCRS > $OUT/${FILENAME}_extractmito.bam
        elif [ -s "$START/$FILENAME.chrM.bam" ]
        then
                echo "The mitochondrial genome is annotated by 'chrM'."
                samtools view -b $START/$FILENAME.bam chrM > $OUT/${FILENAME}_extractmito.bam
	elif [ -s "$START/$FILENAME.M.bam" ]
        then
                echo "The mitochondrial genome is annotated by 'M'."
                samtools view -b $START/$FILENAME.bam M > $OUT/${FILENAME}_extractmito.bam
        fi

rm $START/$FILENAME.bam.bai
rm $START/$FILENAME.MT.bam
rm $START/$FILENAME.M.bam
rm $START/$FILENAME.chrM.bam
rm $START/$FILENAME.chrM_rCRS.bam
