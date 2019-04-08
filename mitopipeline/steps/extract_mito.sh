#!/bin/sh
### $4 is the file size
### $3 is the Cancer Type
### $2 is the Target file ID
### $1 is the download ID

### $1 is file name
### $2 is directory to start at
### $3 is the out directory

START=$2
OUT=$3

#module load base
#module del python
module load intel/17
module load samtools
module load bwa
module load python2/2.7.13

echo extracting chrM bam file
        samtools view -H $START/$1 > $START/$1.header.bam
        grep -i SN:MT $START/$1.header.bam > $START/$1.MT.bam
        grep -i SN:chrM_rCRS $START/$1.header.bam > $START/$1.chrM_rCRS.bam
        grep -i SN:chrM $START/$1.header.bam > $START/$1.chrM.bam
        grep -i SN:M $START/$1.header.bam > $START/$1.M.bam
        
        rm $START/$1.header.bam
        if [ -s "MT.bam" ]
        then
                echo "The mitochondrial genome is annotated by 'MT'."
                samtools view -b $START/$1 MT > $OUT/$1_mito.bam
                rm $START/$1.MT.bam
        elif [ -s "chrM_rCRS.bam" ]
        then
                echo "The mitochondrial genome is annotated by 'chrM_rCRS'."
                samtools view -b $START/$1 chrM_rCRS > $OUT/$1_mito.bam
                rm $START/$1.chrM_rCRS.bam
        elif [ -s "chrM.bam" ]
        then
                echo "The mitochondrial genome is annotated by 'chrM'."
                samtools view -b $START/$1chrM > $OUT/$1_mito.bam
                rm $START/$1
	elif [ -s "M.bam" ]
        then
                echo "The mitochondrial genome is annotated by 'M'."
                samtools view -b $START/$1 M > $OUT/$1_mito.bam
                rm $START/$1.M.bam
        fi
fi