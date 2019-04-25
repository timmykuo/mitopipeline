#!/bin/bash
### $1 IS FILENAME
### $2 IS START DIRECTORY
### $3 IS OUT DIRECTORY WITH FILES
### $4 is tools directory
### $5 is steps directory
### $6 is the refs directory
REF=$6
COUNTS=$3"/counts"
FASTQS=$3"/fastqs"
NEW_BAMS=$3
STOR=$3"/numt_removal_stor"
PILEUPS=$3"/pileups"
START=$2
filename=$1
#last string following / delimeter will be name of the previous job
#only set filetype if the value was a step in the pipeline
filetype="_"$(awk -F/ '{print $NF}' <<< "$2" | awk '{print tolower($0)}')
if [ "$filetype" != "_extractmito" ] && [ "$filetype" != "_clipping" ] && [ "$filetype" != "_splitgap" ];
then
filetype=""
fi

function CountReads {
samtools view -c -F 4 $STOR/$2$1.sorted.bam > $COUNTS/$2$1.count
}

#Align $filename _cl--rCRS $REF/rCRS-MT.fa
function Align {
#use bwa aln and samse for shorter sequences, look into options for bwa aln if inaccurate sam files
bwa aln $3 $START/$1_1$filetype.fastq > $STOR/$1_$2.sai
bwa samse $3 $STOR/$1_$2.sai $START/$1_1$filetype.fastq > $STOR/$1$2.sam
}

#AlignNUMTS _cl--rCRS $filename pileup $REF/rCRS-MT.fa $filename_cl -NM--lowNUMTs.sam lowNUMTs
function AlignNUMTS {
awk 'FNR==NR{a[$1]++;next}!a[$1]' $STOR/$5 $STOR/$2$1.sam > $STOR/$2$1-$6.sam
samtools view -bS $STOR/$2$1-$6.sam > $STOR/$2$1-$6.bam
samtools sort $STOR/$2$1-$6.bam -o $STOR/$2$1-$6.sorted.bam
samtools mpileup -B -C 50 -f $4 $STOR/$2$1-$6.sorted.bam > $PILEUPS/$2$1-$6.$3
}

#$1 = _cl--nuclear.sam
#$2 = $filename
#$3 = _cl-NM
#$4 = _cl--NUMTs.sam
function lowNUMTS {
echo .
grep "NM:i:1" $STOR/$2$1 > $STOR/$2$3
awk '$6==length($10)"M"' $STOR/$2$3 > $STOR/$2$3--NUMTs.sam
echo .
#awk will print first word of each line in the file
#pipes will match unique lines (and prefix lines by num occurences)  and sort based on the first column in _cl-NM--NUMTS.sam and _cl-NUMTs.sam
awk '{print $1}' $STOR/$2$4 | uniq -c | sort -k1 > $STOR/$2$3i0uniq
echo .
awk '{print $1}' $STOR/$2$3--NUMTs.sam | uniq -c | sort -k1 > $STOR/$2$3i1uniq
echo .
#if the first word of each line == 2, print the 2nd word of each line
awk '$1==2{print $2}' $STOR/$2$3i0uniq > $STOR/$2$3i0NUMTsReads
echo .
awk '$1==2{print $2}' $STOR/$2$3i1uniq > $STOR/$2$3i1NUMTsReads
echo .
awk 'FNR==NR{a[$1]++;next}a[$1]' $STOR/$2$3i0NUMTsReads $STOR/$2$1 > $STOR/$2$3i0NUMTs.sam
echo .
awk 'FNR==NR{a[$1]++;next}a[$1]' $STOR/$2$3i1NUMTsReads $STOR/$2$1 > $STOR/$2$3i1NUMTs.sam
echo .
awk '$1==1{print $2}' $STOR/$2$3i0uniq > $STOR/$2$3i0NUMTsSingleReads
echo .
awk '$1==1{print $2}' $STOR/$2$3i1uniq > $STOR/$2$3i1NUMTsSingleReads
echo SINGLE-END, removing all singlet reads which align with up to 1 mismatch to nuclear genome
cat $STOR/$2$3i0NUMTsSingleReads $STOR/$2$3i1NUMTsSingleReads > $STOR/$2$3--lowNUMTs.sam
}
echo .
echo .
echo .
echo ----bwa alignment to rCRS
#aligned file is in STOR
Align ${filename} _cl--rCRS $REF/rCRS-MT.fa
echo ****bwa alignment to rCRS DONE.
echo .
echo .
echo .
echo .
echo .
echo .
echo ----bwa alignment to hg38-norCRS --NUMTs: going to make ${filename}_cl--nuclear.sam
#use bwa aln and samse for short reads
bwa aln $REF/hg38-nochr.fa $START/${filename}_1${filetype}.fastq > $STOR/${filename}_cl--nuclear.sai
bwa samse $REF/hg38-nochr.fa $STOR/${filename}_cl--nuclear.sai $START/${filename}_1${filetype}.fastq > $STOR/${filename}_cl--nuclear.sam 
echo ****bwa alignment to hg38-norCRS DONE.
echo .
echo .
echo .
echo .
echo .
echo .
echo ----Extract Perfect Matches to Nuclear Genome.
grep `awk 'END {print "MD:Z:"length($10)}' $STOR/${filename}_cl--nuclear.sam` $STOR/${filename}_cl--nuclear.sam > $STOR/${filename}_cl--NUMTs.sam
lowNUMTS _cl--nuclear.sam ${filename} _cl-NM _cl--NUMTs.sam
echo ****Extraction of Perfect Matches to Nuclear Genome DONE.
echo .
echo .
echo .
echo .
echo .
echo .
echo ----Generate ${filename}_cl--rCRS-lowNUMTS.pileup
AlignNUMTS _cl--rCRS ${filename} pileup $REF/rCRS-MT.fa ${filename}_cl-NM--lowNUMTs.sam lowNUMTs
CountReads _cl--rCRS-lowNUMTs ${filename}
echo ****Generate ${filename}_cl--rCRS-lowNUMTS.pileup DONE.
echo ****NUMTs DONE.
echo .
echo .
echo .
echo .
echo .
echo .
echo ----BAM2FASTQ ${filename}_cl--rCRS-lowNUMTs.sam Initiated
samtools view -bS $STOR/${filename}_cl--rCRS-lowNUMTs.sam > $STOR/${filename}_cl--rCRS-lowNUMTs.bam
samtools fastq $STOR/${filename}_cl--rCRS-lowNUMTs.bam > $STOR/${filename}_cl--rCRS-lowNUMTs.fastq
echo ****BAM2FASTQ DONE.
echo .
echo .
echo .
echo .
echo .
echo .
echo ****bwa alignment of ${filename}_cl--rCRS-lowNUMTs fastqs to hg38
bwa aln $REF/rCRS-MT.fa $STOR/${filename}_cl--rCRS-lowNUMTs.fastq > $STOR/${filename}_cl--rCRS-lowNUMTs.sai
bwa samse $REF/rCRS-MT.fa $STOR/${filename}_cl--rCRS-lowNUMTs.sai $STOR/${filename}_cl--rCRS-lowNUMTs.fastq > $STOR/${filename}.mito_hg38.sam
echo ****bwa DONE.
echo .
echo .
echo .
echo .
echo .
echo .
echo ----samtools Started
samtools view -bS $STOR/${filename}.mito_hg38.sam > $STOR/${filename}.mito_hg38.bam
samtools sort $STOR/${filename}.mito_hg38.bam -o $STOR/${filename}.mito_hg38.sorted.bam
samtools index $STOR/${filename}.mito_hg38.sorted.bam
samtools view -b $STOR/${filename}.mito_hg38.sorted.bam MT > $NEW_BAMS/${filename}_removenumts.bam
echo ***** samtools DONE
echo .
echo .
echo .
echo .
echo Moving files
mv $STOR/${filename}_cl--rCRS-lowNUMTs.fastq $FASTQS
