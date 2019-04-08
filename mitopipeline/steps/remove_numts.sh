#!/bin/bash
### $4 is the file size
### $3 is the Cancer Type
### $3 is the Target file ID
### $1 is the download ID

### $1 IS FILENAME
### $2 IS START DIRECTORY
### $3 IS OUT DIRECTORY WITH FILES

REF="./REFs"
COUNTS=$3"/counts"
FASTQS="/mnt/rds/txl80/LaframboiseLab/tyk3/TARGET_PIPELINE/fastqs"
NEW_BAMS=$3
STOR=$3"/numt_removal_stor"
PILEUPS=$3"/pileups"
START=$2

module load intel/17
module load openmpi/2.0.1
module load samtools
module load bwa

function CountReads {
samtools view -c -F 4 $START/$3_$2$1.sorted.bam > $COUNTS/$3_$2$1.count
}

#Align ${id} ${cancer_type} _cl--rCRS $REF/MitoREFs/rCRS-MT.fa
function Align {
#use bwa aln and samse for shorter sequences, look into options for bwa aln if inaccurate sam files
bwa aln $4 $FASTQS/$2_$1.fastq > $STOR/$2_$1$3.sai
bwa samse $4 $STOR/$2_$1$3.sai $FASTQS/$2_$1.fastq > $STOR/$2_$1$3.sam
}

#AlignNUMTS _cl--rCRS ${id} ${cancer_type} pileup $REF/MitoREFs/rCRS-MT.fa ${cancer_type}_${id}_cl    -NM--lowNUMTs.sam lowNUMTs
function AlignNUMTS {
awk 'FNR==NR{a[$1]++;next}!a[$1]' $STOR/$6 $STOR/$3_$2$1.sam > $STOR/$3_$2$1-$7.sam
samtools view -bS $STOR/$3_$2$1-$7.sam > $STOR/$3_$2$1-$7.bam
samtools sort $STOR/$3_$2$1-$7.bam -o $STOR/$3_$2$1-$7.sorted.bam
#
samtools mpileup -B -C 50 -f $5 $STOR/$3_$2$1-$7.sorted.bam > $PILEUPS/$3_$2$1-$7.$4
}

#$1 = _cl--nuclear.sam
#$2 = ${id}
#$3 = ${cancer_type}
#$4 = _cl-NM
#$5 = _cl--NUMTs.sam
function lowNUMTS {
grep "NM:i:1" $3_$2$1 > $3_$2$4
awk '$6==length($10)"M"' $3_$2$4 > $3_$2$4--NUMTs.sam

#awk will print first word of each line in the file
#pipes will match unique lines (and prefix lines by num occurences)  and sort based on the first column in _cl-NM--NUMTS.sam and _cl-NUMTs.sam
awk '{print $1}' $STOR/$3_$2$5 | uniq -c | sort -k1 > $STOR/$3_$2$4i0uniq
awk '{print $1}' $STOR/$3_$2$4--NUMTs.sam | uniq -c | sort -k1 > $STOR/$3_$2$4i1uniq

#if the first word of each line == 2, print the 2nd word of each line
awk '$1==2{print $2}' $STOR/$3_$2$4i0uniq > $STOR/$3_$2$4i0NUMTsReads
awk '$1==2{print $2}' $STOR/$3_$2$4i1uniq > $STOR/$3_$2$4i1NUMTsReads

awk 'FNR==NR{a[$1]++;next}a[$1]' $STOR/$3_$2$4i0NUMTsReads $STOR/$3_$2$1 > $STOR/$3_$2$4i0NUMTs.sam
awk 'FNR==NR{a[$1]++;next}a[$1]' $STOR/$3_$2$4i1NUMTsReads $STOR/$3_$2$1 > $STOR/$3_$2$4i1NUMTs.sam

awk '$1==1{print $2}' $STOR/$3_$2$4i0uniq > $STOR/$3_$2$4i0NUMTsSingleReads
awk '$1==1{print $2}' $STOR/$3_$2$4i1uniq > $STOR/$3_$2$4i1NUMTsSingleReads

#awk 'FNR==NR{a[$1]++;next}a[$1]' $FASTQS/$3_$2$4i0NUMTsSingleReads $FASTQS/$3_$2$4i1NUMTsSingleReads > $FASTQS/$3_$2$4NUMTsDupsReads
#awk 'FNR==NR{a[$1]++;next}a[$1]' $FASTQS/$3_$2$4NUMTsDupsReads $FASTQS/$3_$2$1 >  $FASTQS/$3_$2$4NUMTsDups.sam

#if [ -e "$FASTQS/$3_$2_mito_1.fastq" ]
#then
#        echo PAIRED-END, removing all mate pairs which align with up to 1 mismatch to nuclear genome
#        cat $FASTQS/$3_$2$4i0NUMTs.sam $FASTQS/$3_$2$4i1NUMTs.sam $FASTQS/$3_$2$4NUMTsDups.sam > $FASTQS/$3_$2$4--lowNUMTs.sam
#else
        echo SINGLE-END, removing all singlet reads which align with up to 1 mismatch to nuclear genome
        cat $STOR/$3_$2$4i0NUMTsSingleReads $STOR/$3_$2$4i1NUMTsSingleReads > $STOR/$3_$2$4--lowNUMTs.sam
#fi
}

#remove / up to the last pathname
filename=$1
echo $filename
cancer_type=$(awk -F_ '{print $1}' <<< "$filename")
id=$(awk -F_ '{print $2"_"$3}' <<< "$filename")

echo .
echo .
echo .
echo ----bwa alignment to rCRS
#aligned file is in STOR
Align ${id} ${cancer_type} _cl--rCRS $REF/MitoREFs/rCRS-MT.fa
echo ****bwa alignment to rCRS DONE.
echo .
echo .
echo .
echo .
echo .
echo .
echo ----bwa alignment to hg38-norCRS --NUMTs: going to make ${cancer_type}_${id}_cl--nuclear.sam
cd $STOR
#use bwa aln and samse for short reads
bwa aln $REF/hg38-nochr.fa $FASTQS/${cancer_type}_${id}.fastq > $STOR/${cancer_type}_${id}_cl--nuclear.sai
bwa samse $REF/hg38-nochr.fa $STOR/${cancer_type}_${id}_cl--nuclear.sai $FASTQS/${cancer_type}_${id}.fastq > $STOR/${cancer_type}_${id}_cl--nuclear.sam 
echo ****bwa alignment to hg38-norCRS DONE.
echo .
echo .
echo .
echo .
echo .
echo .
echo ----Extract Perfect Matches to Nuclear Genome.
grep `awk 'END {print "MD:Z:"length($10)}' ${cancer_type}_${id}_cl--nuclear.sam` ${cancer_type}_${id}_cl--nuclear.sam > ${cancer_type}_${id}_cl--NUMTs.sam
lowNUMTS _cl--nuclear.sam ${id} ${cancer_type} _cl-NM _cl--NUMTs.sam
echo ****Extraction of Perfect Matches to Nuclear Genome DONE.
echo .
echo .
echo .
echo .
echo .
echo .
echo ----Generate ${cancer_type}_${id}_cl--rCRS-lowNUMTS.pileup
AlignNUMTS _cl--rCRS ${id} ${cancer_type} pileup $REF/MitoREFs/rCRS-MT.fa ${cancer_type}_${id}_cl-NM--lowNUMTs.sam lowNUMTs
CountReads _cl--rCRS-lowNUMTs ${id} ${cancer_type}
echo ****Generate ${cancer_type}_${id}_cl--rCRS-lowNUMTS.pileup DONE.
echo ****NUMTs DONE.
echo .
echo .
echo .
echo .
echo .
echo .
echo ----BAM2FASTQ ${cancer_type}_${id}_cl--rCRS-lowNUMTs.sam Initiated
samtools view -bS $STOR/${cancer_type}_${id}_cl--rCRS-lowNUMTs.sam > $STOR/${cancer_type}_${id}_cl--rCRS-lowNUMTs.bam
#/home/sxg501/local/bam2fastq-1.1.0/bam2fastq -o $FASTQS/${cancer_type}_${id}_cl--rCRS-lowNUMTs#.fastq $RAW_BAMS/${cancer_type}_${id}_cl--rCRS-lowNUMTs.bam
samtools fastq $STOR/${cancer_type}_${id}_cl--rCRS-lowNUMTs.bam > $STOR/${cancer_type}_${id}_cl--rCRS-lowNUMTs.fastq
echo ****BAM2FASTQ DONE.
echo .
echo .
echo .
echo .
echo .
echo .
echo ****bwa alignment of ${cancer_type}_${id}_cl--rCRS-lowNUMTs fastqs to hg38
bwa aln $REF/MitoREFs/rCRS-MT.fa $STOR/${cancer_type}_${id}_cl--rCRS-lowNUMTs.fastq > $STOR/${cancer_type}_${id}_cl--rCRS-lowNUMTs.sai
bwa samse $REF/MitoREFs/rCRS-MT.fa $STOR/${cancer_type}_${id}_cl--rCRS-lowNUMTs.sai $STOR/${cancer_type}_${id}_cl--rCRS-lowNUMTs.fastq > $STOR/${cancer_type}_${id}.mito_hg38.sam
echo ****bwa DONE.
echo .
echo .
echo .
echo .
echo .
echo .

echo ----samtools Started
samtools view -bS $STOR/${cancer_type}_${id}.mito_hg38.sam > $STOR/${cancer_type}_${id}.mito_hg38.bam
samtools sort $STOR/${cancer_type}_${id}.mito_hg38.bam -o $STOR/${cancer_type}_${id}.mito_hg38.sorted.bam
samtools index $STOR/${cancer_type}_${id}.mito_hg38.sorted.bam
samtools view -b $STOR/${cancer_type}_${id}.mito_hg38.sorted.bam MT > $NEW_BAMS/${cancer_type}_${id}.mito_hg38.sorted.mt.bam
echo ***** samtools DONE
echo .
echo .
echo .
echo .
echo Moving files
mv $STOR/${cancer_type}_${id}_cl--rCRS-lowNUMTs.fastq $FASTQS
