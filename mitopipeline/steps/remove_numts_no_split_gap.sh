### $4 is the file size
### $3 is the Cancer Type
### $2 is the Target file ID
### $1 is the download ID


#$1 is filename
#$2 is start directory
#$3 is out directory

REF="./REFs"
COUNTS=$3"/counts"
FASTQS=$3"/fastqs"
NEW_BAMS=$3
STOR=$3"/numt_removal_stor"
PILEUPS=$3"/pileups"
START=$2
TOOLS=$4

module load samtools
module load bwa
module load gdc/1.3.0

function CountReads {
samtools view -c -F 4 $START/$2$1.sorted.bam > $COUNTS/$2$1.count
}

#Align $filename _cl--rCRS $REF/rCRS-MT.fa
function Align {
#use bwa aln and samse for shorter sequences, look into options for bwa aln if inaccurate sam files
bwa aln $3 $FASTQS/$1.fastq > $STOR/$1_$2.sai
bwa samse $3 $STOR/$1_$2.sai $FASTQS/$1.fastq > $STOR/$1_$2.sam
}


function Align {
if [ -e "$FASTQS/$1_1.fastq" ]
then
bwa mem -t 12 $4 $FASTQS/$1_1.fastq $FASTQS/$1_2.fastq > $STOR/$1$2.sam
else
bwa mem -t 12 $4 $FASTQS/$1.fastq > $STOR/$1$2.sam
fi
}

#AlignNUMTS _cl--rCRS ${id} ${cancer_type} pileup $REF/rCRS-MT.fa ${cancer_type}_${id}_cl-NM--lowNUMTs.sam lowNUMTs
#AlignNUMTS _cl--rCRS $filename pileup $REF/rCRS-MT.fa $filename_cl-NM--lowNUMTs.sam lowNUMTs
function AlignNUMTS {
awk 'FNR==NR{a[$1]++;next}!a[$1]' $STOR/$5 $STOR/$2$1.sam > $STOR/$2$1-$6.sam
samtools view -bS $STOR/$2$1-$6.sam > $STOR/$2$1-$6.bam
samtools sort $STOR/$2$1-$6.bam $STOR/$2$1-$6.sorted.bam
samtools mpileup -B -C 50 -q 40 -Q 30 -f $4 $STOR/$2$1-$6.sorted.bam > $PILEUPS/$2$1-$6.$3
}


#$1 = _cl--nuclear.sam
#$2 = ${id}
#$3 = ${cancer_type}
#$4 = _cl-NM
#$5 = _cl--NUMTs.sam

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

if [ -e "$FASTQS/$2_1.fastq" ]
then
	echo PAIRED-END, removing all mate pairs which align with up to 1 mismatch to nuclear genome
	cat $STOR/$2$3i0NUMTs.sam $STOR/$2$3i1NUMTs.sam $STOR/$2$3NUMTsDups.sam > $STOR/$2$3--lowNUMTs.sam
else
	echo SINGLE-END, removing all singlet reads which align with up to 1 mismatch to nuclear genome
	cat $STOR/$2$3i0NUMTsSingleReads $STOR/$2$3i1NUMTsSingleReads > $STOR/$2$3--lowNUMTs.sam
fi
}

echo ----bwa alignment to rCRS
Align _mito_CLIPPED $1 _cl--rCRS $REF/rCRS.fa pileup
echo ****bwa alignment to rCRS DONE.
echo .
echo .
echo .
echo .
echo .
echo .
echo ----bwa alignment to hg38-norCRS --NUMTs: going to make $1_cl--nuclear.sam
if [ -e "$STOR/$1_CLIPPED_1.fastq" ]
then
bwa mem -t 12 $REF/hg38-nochr.fa $FASTQS/$1_mito_CLIPPED_1.fastq $FASTQS/$1_mito_CLIPPED_2.fastq > $STOR/$1_cl--nuclear.sam
else
bwa mem -t 12 $REF/hg38-nochr.fa $FASTQS/$1_mito_CLIPPED.fastq > $STOR/$1_cl--nuclear.sam
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
AlignNUMTS _cl--rCRS $1 pileup $REF/rCRS.fa $1_cl-NM--lowNUMTs.sam lowNUMTs
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
$TOOLS/bam2fastq-1.1.0/bam2fastq -o $STOR/$1_cl--rCRS-lowNUMTs#.fastq $STOR/$1_cl--rCRS-lowNUMTs.bam
echo ****BAM2FASTQ DONE.
echo .
echo .
echo .
echo .
echo .
echo .
echo ****bwa alignment of $1_cl--rCRS-lowNUMTs fastqs to hg38
if [ -e "$STOR/$1_cl--rCRS-lowNUMTs_1.fastq" ]
then
echo "PAIRED-END"
	bwa aln $REF/hg38.fa $STOR/$1_cl--rCRS-lowNUMTs_1.fastq > $STOR/$1_cl--rCRS-lowNUMTs_1.sai
	bwa aln $REF/hg38.fa $STOR/$1_cl--rCRS-lowNUMTs_2.fastq > $STOR/$1_cl--rCRS-lowNUMTs_2.sai
	bwa sampe $REF/hg38.fa $STOR/$1_cl--rCRS-lowNUMTs_1.sai $STOR/$1_cl--rCRS-lowNUMTs_2.sai $STOR/$1_cl--rCRS-lowNUMTs_1.fastq $STOR/$1_cl--rCRS-lowNUMTs_2.fastq > $STOR/$1.mito_hg38.sam
else
echo "SINGLE-END"
	bwa aln $REF/hg38.fa $STOR/$1_cl--rCRS-lowNUMTs.fastq > $STOR/$1_cl--rCRS-lowNUMTs.sai
	bwa samse $REF/hg38.fa $STOR/$1_cl--rCRS-lowNUMTs.sai $STOR/$1_cl--rCRS-lowNUMTs.fastq > $STOR/$1.mito_hg38.sam
fi
echo ****bwa DONE.
echo .
echo .
echo .
echo .
echo .
echo .

echo ----samtools Started
samtools view -bS $STOR/$1.mito_hg38.sam > $STOR/$1.mito_hg38.bam
samtools sort $STOR/$1.mito_hg38.bam $STOR/$1.mito_hg38.sorted.bam
samtools index $SCR/$1.mito_hg38.sorted.bam
samtools view -b $STOR/$1.mito_hg38.sorted.bam MT > $NEW_BAMS/$1_removenumts.bam
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
mv $STOR/$1_cl--rCRS-lowNUMTs.fastq $FASTQS
# mv $STOR/$1.mito_hg38.sorted.bam $OUTdir/gdc-mito-bams/
# mv $SCR/$1.mito_hg38.sorted.bam.bai $OUTdir/gdc-mito-bams/
# mv $SCR/$1_cl--rCRS-lowNUMTs.pileup $OUTdir/gdc-mito-pileups/
# rm $SCR/$1*
# rm $STOR/$1*