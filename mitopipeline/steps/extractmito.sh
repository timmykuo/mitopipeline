#!/bin/sh
### $1 is file name
### $2 is directory to start at
### $3 is the out directory
### $4 is tools directory
### $5 is steps directory
### $6 is refs directory

START=$2
OUT=$3
TOOLS=$4

module load samtools

echo extracting chrM bam file
        samtools index $START/$1.bam
        samtools view -H $START/$1.bam > $START/$1.header.bam
        grep -i SN:MT $START/$1.header.bam > $START/$1.MT.bam
        grep -i SN:chrM_rCRS $START/$1.header.bam > $START/$1.chrM_rCRS.bam
        grep -i SN:chrM $START/$1.header.bam > $START/$1.chrM.bam
        grep -i SN:M $START/$1.header.bam > $START/$1.M.bam
        
        rm $START/$1.header.bam
        if [ -s "$START/$1.MT.bam" ]
        then
                echo "The mitochondrial genome is annotated by 'MT'."
                samtools view -b $START/$1.bam MT > $OUT/$1_extractmito.bam
        elif [ -s "$START/$1.chrM_rCRS.bam" ]
        then
                echo "The mitochondrial genome is annotated by 'chrM_rCRS'."
                samtools view -b $START/$1.bam chrM_rCRS > $OUT/$1_extractmito.bam
        elif [ -s "$START/$1.chrM.bam" ]
        then
                echo "The mitochondrial genome is annotated by 'chrM'."
                samtools view -b $START/$1.bam chrM > $OUT/$1_extractmito.bam
	elif [ -s "$START/$1.M.bam" ]
        then
                echo "The mitochondrial genome is annotated by 'M'."
                samtools view -b $START/$1.bam M > $OUT/$1_extractmito.bam
        fi
rm $START/$1.bam.bai
rm $START/$1.MT.bam
rm $START/$1.M.bam
rm $START/$1.chrM.bam
rm $START/$1.chrM_rCRS.bam
