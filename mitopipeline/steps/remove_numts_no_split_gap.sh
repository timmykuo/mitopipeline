#!/bin/bash
#$1 is filename
#$2 is start directory
#$3 is out directory
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
TOOLS=$4
filename=$1
currdir=`pwd`
FAILED="False"
#last string following / delimeter will be name of the previous job
#only set filetype if the value was a step in the pipeline
filetype="_"$(awk -F/ '{print $NF}' <<< "$2" | awk '{print tolower($0)}')
if [ "$filetype" != "_extractmito" ] && [ "$filetype" != "_clipping" ] && [ "$filetype" != "_splitgap" ];
then
filetype=""
fi
LENGTH=`python $5/read_depths.py $START/$1_1$filetype.fastq`
function CountReads {
samtools view -c -F 4 $STOR/$2$1.sorted.bam > $COUNTS/$2$1.count
}

#Align $filename _cl--rCRS $REF/rCRS-MT.fa
function Align {
if [ -e "$START/$1_1$filetype.fastq" ]
then
	echo PAIRED_END 
	if (( $LENGTH > 70 ))
	then
	echo "Average read depth is $LENGTH. Using bwa-mem algorithm"
	bwa mem -M -t 12 $3 $START/$1_1$filetype.fastq $START/$1_2$filetype.fastq > $STOR/$1$2.sam
	if [ $? -ne 0 ]; then
		FAILED="True"
	fi
	else
	echo "Average read depth is $LENGTH. Using bwa-backtrack algorithm"
	bwa aln $3 $START/$1_1$filetype.fastq > $STOR/$1$2_1.sai
	if [ $? -ne 0 ]; then
		FAILED="True"
	fi
	bwa aln $3 $START/$1_2$filetype.fastq > $STOR/$1$2_2.sai
	if [ $? -ne 0 ]; then
		FAILED="True"
	fi
	bwa sampe $3 $STOR/$1$2_1.sai $STOR/$1$2_2.sai $START/$1_1$filetype.fastq $START/$1_2$filetype.fastq > $STOR/$1$2.sam
	if [ $? -ne 0 ]; then
		FAILED="True"
	fi
	fi
else
	echo SINGLE_END
	if (( $LENGTH > 70 ))
	then
	echo "Average read depth is $LENGTH. Using bwa-mem algorithm"
	bwa mem -M -t 12 $3 $START/$1$filetype.fastq > $STOR/$1$2.sam
	if [ $? -ne 0 ]; then
		FAILED="True"
	fi
	else
	echo "Average read depth is $LENGTH. Using bwa-backtrack algorithm"
	bwa aln $3 $START/$1$filetype.fastq > $STOR/$1$2.sai
	if [ $? -ne 0 ]; then
		FAILED="True"
	fi
	bwa samse $3 $STOR/$1$2.sai $START/$1$filetype.fastq > $STOR/$1$2.sam
	if [ $? -ne 0 ]; then
		FAILED="True"
	fi
	fi
fi
}

#AlignNUMTS _cl--rCRS $filename pileup $REF/rCRS-MT.fa $filename_cl-NM--lowNUMTs.sam lowNUMTs
function AlignNUMTS {
awk 'FNR==NR{a[$1]++;next}!a[$1]' $STOR/$5 $STOR/$2$1.sam > $STOR/$2$1-$6.sam
samtools view -bS $STOR/$2$1-$6.sam > $STOR/$2$1-$6.bam
samtools sort $STOR/$2$1-$6.bam > $STOR/$2$1-$6.sorted.bam
samtools mpileup -B -C 50 -q 40 -Q 30 -f $4 $STOR/$2$1-$6.sorted.bam > $PILEUPS/$2$1-$6.$3
}

#$1 = _cl--nuclear.sam
#$2 = $filename
#$3 = _cl-NM
#$4 = _cl--NUMTs.sam
function lowNUMTS {
grep "NM:i:1" $STOR/$2$1 > $STOR/$2$3
awk '$6==length($10)"M"' $STOR/$2$3 >  $STOR/$2$3--NUMTs.sam

awk '{print $1}' $STOR/$2$4 | uniq -c | sort -k1 > $STOR/$2$3i0uniq
awk '{print $1}' $STOR/$2$3--NUMTs.sam | uniq -c | sort -k1 > $STOR/$2$3i1uniq

awk '$1==2{print $2}' $STOR/$2$3i0uniq > $STOR/$2$3i0NUMTsReads
awk '$1==2{print $2}' $STOR/$2$3i1uniq > $STOR/$2$3i1NUMTsReads

awk 'FNR==NR{a[$1]++;next}a[$1]' $STOR/$2$3i0NUMTsReads $STOR/$2$1 > $STOR/$2$3i0NUMTs.sam
awk 'FNR==NR{a[$1]++;next}a[$1]' $STOR/$2$3i1NUMTsReads $STOR/$2$1 > $STOR/$2$3i1NUMTs.sam

awk '$1==1{print $2}' $STOR/$2$3i0uniq > $STOR/$2$3i0NUMTsSingleReads
awk '$1==1{print $2}' $STOR/$2$3i1uniq > $STOR/$2$3i1NUMTsSingleReads

awk 'FNR==NR{a[$1]++;next}a[$1]' $STOR/$2$3i0NUMTsSingleReads $STOR/$2$3i1NUMTsSingleReads > $STOR/$2$3NUMTsDupsReads
awk 'FNR==NR{a[$1]++;next}a[$1]' $STOR/$2$3NUMTsDupsReads $STOR/$2$1 >  $STOR/$2$3NUMTsDups.sam

if [ -e "$START/$2_1$filetype.fastq" ]
then
	echo PAIRED-END, removing all mate pairs which align with up to 1 mismatch to nuclear genome
	cat $STOR/$2$3i0NUMTs.sam $STOR/$2$3i1NUMTs.sam $STOR/$2$3NUMTsDups.sam > $STOR/$2$3--lowNUMTs.sam
else
	echo SINGLE-END, removing all singlet reads which align with up to 1 mismatch to nuclear genome
	cat $STOR/$2$3i0NUMTsSingleReads $STOR/$2$3i1NUMTsSingleReads > $STOR/$2$3--lowNUMTs.sam
fi
}

echo ----bwa alignment to rCRS
Align $1 _cl--rCRS $REF/rCRS-MT.fa
echo ****bwa alignment to rCRS DONE.
echo .
echo .
echo .
echo .
echo .
echo .
echo ----bwa alignment to hg38-norCRS --NUMTs: going to make $1_cl--nuclear.sam
if [ -e "$START/$1_1$filetype.fastq" ]
then
	echo PAIRED-END
	if (( $LENGTH > 70 ))
	then
	echo "Average read depth is $LENGTH. Using bwa-mem algorithm"
	bwa mem -M -t 12 $REF/hg38-nochr.fa $START/$1_1$filetype.fastq $START/$1_2$filetype.fastq > $STOR/$1_cl--nuclear.sam
	if [ $? -ne 0 ]; then
		FAILED="True"
	fi
	else
	echo "Average read depth is $LENGTH. Using bwa-backtrack algorithm"
	bwa aln $REF/hg38-nochr.fa $START/$1_1$filetype.fastq > $STOR/$1_1_hg38-nochr.sai
	if [ $? -ne 0 ]; then
		FAILED="True"
	fi
	bwa aln $REF/hg38-nochr.fa $START/$1_2$filetype.fastq > $STOR/$1_2_hg38-nochr.sai
	if [ $? -ne 0 ]; then
		FAILED="True"
	fi
	bwa sampe $REF/hg38-nochr.fa $STOR/$1_1_hg38-nochr.sai $STOR/$1_2_hg38-nochr.sai $START/$1_1$filetype.fastq $START/$1_2$filetype.fastq > $STOR/$1_cl--nuclear.sam
	if [ $? -ne 0 ]; then
		FAILED="True"
	fi
	fi 
else
	echo SINGLE-END
	if (( $LENGTH > 70 ))
	then
	echo "Average read depth is $LENGTH. Using bwa-mem algorithm"
	bwa mem -M -t 12 $REF/hg38-nochr.fa $START/$1$filetype.fastq > $STOR/$1_cl--nuclear.sam
	if [ $? -ne 0 ]; then
		FAILED="True"
	fi
	else
	echo "Average read depth is $LENGTH. Using bwa-backtrack algorithm"
	bwa aln $REF/hg38-nochr.fa $START/$1$filetype.fastq > $STOR/$1_hg38-nochr.sai
	if [ $? -ne 0 ]; then
		FAILED="True"
	fi
	bwa samse $REF/hg38-nochr.fa $STOR/$1_hg38-nochr.sai $START/$1$filetype.fastq > $STOR/$1_cl--nuclear.sam
	if [ $? -ne 0 ]; then
		FAILED="True"
	fi
	fi
fi
echo ****bwa alignment to hg38-norCRS DONE.
echo .
echo .
echo .
echo .
echo .
echo .
echo ----Extract Perfect Matches to Nuclear Genome.
grep `awk 'END {print "MD:Z:"length($10)}' $STOR/$1_cl--nuclear.sam` $STOR/$1_cl--nuclear.sam > $STOR/$1_cl--NUMTs.sam
lowNUMTS _cl--nuclear.sam $1 _cl-NM _cl--NUMTs.sam
echo ****Extraction of Perfect Matches to Nuclear Genome DONE.
echo .
echo .
echo .
echo .
echo .
echo .
echo ----Generate $1_cl--rCRS-lowNUMTS.pileup
AlignNUMTS _cl--rCRS $1 pileup $REF/rCRS-MT.fa $1_cl-NM--lowNUMTs.sam lowNUMTs
CountReads _cl--rCRS-lowNUMTs $1
echo ****Generate $1_cl--rCRS-lowNUMTS.pileup DONE.
echo ****NUMTs DONE.
echo .
echo .
echo .
echo .
echo .
echo .
echo ----BAM2FASTQ $1_cl--rCRS-lowNUMTs.sam Initiated
samtools view -bS $STOR/$1_cl--rCRS-lowNUMTs.sam > $STOR/$1_cl--rCRS-lowNUMTs.bam
bedtools bamtofastq -i $STOR/$1_cl--rCRS-lowNUMTs.bam -fq $STOR/$1_cl--rCRS-lowNUMTs.fastq -fq2 $STOR/$1_cl--rCRS-lowNUMTs_2.fastq 
echo ****BAM2FASTQ DONE.
echo .
echo .
echo .
echo .
echo .
echo .
LENGTH=`python $5/read_depths.py $STOR/$1_cl--rCRS-lowNUMTs.fastq`
echo ****bwa alignment of $1_cl--rCRS-lowNUMTs fastqs to hg38
if [ -e "$STOR/$1_cl--rCRS-lowNUMTs_2.fastq" ]
then
	echo "PAIRED-END"
	if (( $LENGTH < 70 ))
	then
		echo "average read_depth is $LENGTH. Using bwa-backtrack algorithm"
		bwa aln $REF/hg38.fa $STOR/$1_cl--rCRS-lowNUMTs.fastq > $STOR/$1_cl--rCRS-lowNUMTs_1.sai
		if [ $? -ne 0 ]; then
			FAILED="True"
		fi
		bwa aln $REF/hg38.fa $STOR/$1_cl--rCRS-lowNUMTs_2.fastq > $STOR/$1_cl--rCRS-lowNUMTs_2.sai
		if [ $? -ne 0 ]; then
			FAILED="True"
		fi
		bwa sampe $REF/hg38.fa $STOR/$1_cl--rCRS-lowNUMTs_1.sai $STOR/$1_cl--rCRS-lowNUMTs_2.sai $STOR/$1_cl--rCRS-lowNUMTs.fastq $STOR/$1_cl--rCRS-lowNUMTs_2.fastq > $STOR/$1.mito_hg38.sam
		if [ $? -ne 0 ]; then
			FAILED="True"
		fi
	else
		echo "average read_depth is $LENGTH. Using bwa-mem algorithm"
		bwa mem -M -t 12 $REF/hg38.fa $STOR/$1_cl--rCRS-lowNUMTs.fastq $STOR/$1_cl--rCRS-lowNUMTs_2.fastq > $STOR/$1.mito_hg38.sam
		if [ $? -ne 0 ]; then
			FAILED="True"
		fi
	fi
else
	echo "SINGLE-END"
	if (( $LENGTH < 70 ))
	then
		echo "average read_depth is $LENGTH. Using bwa-backtrack algorithm"
		bwa aln $REF/hg38.fa $STOR/$1_cl--rCRS-lowNUMTs.fastq > $STOR/$1_cl--rCRS-lowNUMTs.sai
		if [ $? -ne 0 ]; then
			FAILED="True"
		fi
		bwa samse $REF/hg38.fa $STOR/$1_cl--rCRS-lowNUMTs.sai $STOR/$1_cl--rCRS-lowNUMTs.fastq > $STOR/$1.mito_hg38.sam
		if [ $? -ne 0 ]; then
			FAILED="True"
		fi
	else
		echo "average read_depth is $LENGTH. Using bwa-mem algorithm"
		bwa mem -M -t 12 $REF/hg38.fa $STOR/$1_cl--rCRS-lowNUMTs.fastq > $STOR/$1.mito_hg38.sam
		if [ $? -ne 0 ]; then
			FAILED="True"
		fi
	fi
fi
echo ****bwa DONE.
echo .
echo .
echo .
echo .
echo .
echo .

echo ----samtools Started
if [ "$FAILED" == "True" ]; then
	echo "RemoveNuMTs failed – likely due to memory capacity. Please increase memory capacity or use the -l slurm option"
else
	samtools view -bS ${STOR}/${1}.mito_hg38.sam > ${STOR}/${1}.mito_hg38.bam
	samtools sort ${STOR}/${1}.mito_hg38.bam > ${STOR}/${1}.mito_hg38.sorted.bam
	samtools index ${STOR}/${1}.mito_hg38.sorted.bam
	samtools view -b ${STOR}/${1}.mito_hg38.sorted.bam MT > ${NEW_BAMS}/${1}_removenumts.bam
fi
#mv $STOR/$1.mito_hg38.sorted.bam ${NEW_BAMS}/${1}_removenumts.bam
echo ****samtools DONE.
echo .
echo .
echo .
echo .
echo .
echo .
echo .
echo .
echo .
echo .
echo .
echo .
echo ----Moving / Deleting Files