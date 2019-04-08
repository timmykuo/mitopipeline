#########################################
### Written by Sneha Grandhi, 8.22.15 ###
### This script downloads a slice of  ###
##  mito and one two-copy nuclear .bam ##
#   $1: Genome Build 
#   $2: Patient ID
#   $3: TCGA Cancer Type
#   $4: Analysis ID
#   $5: Mito Annotation
#   $6: Mito Region
#   $7: Line Number
#   $8: Nuc 2-Copy Region

module load samtools
module load bwa

SCR="/scratch/pbsjobs/sneha/BAMSlicer_OUTPUT"
SCR2="/scratch/pbsjobs/sneha/BAMSlicer_tmp"
REF="/home/sxg501/REFs"
SXG="/home/sxg501"
SCR3="/scratch/pbsjobs/sneha/"

echo ----Samtools on $3_$2.input.bam Initiated
samtools sort $SCR/$3_$2.input.bam $SCR/$3_$2.input.sorted
samtools index $SCR/$3_$2.input.sorted.bam
#rm $SCR/$3_$2.input.bam
echo ****Samtools on $3_$2.input.bam DONE.
echo .
echo .
echo .
echo .
echo .
echo ----Nuclear Copy Number 2 Extraction Initiated
samtools view -b $SCR/$3_$2.input.sorted.bam $8 > $SCR2/$3_$2.Nuc2.bam
samtools view -c -F 4 $SCR2/$3_$2.Nuc2.bam > $SCR/$3_$2_COPY2.count
echo ****Nuclear Copy Number Extraction DONE.
echo .
echo .
echo .
echo .
echo .
echo .
echo ----Mito .bam Extraction Initiated
samtools view -b $SCR/$3_$2.input.sorted.bam $5 > $SCR/$3_$2.mito.bam
echo ****Mito .bam Extraction DONE.
echo .
echo .
echo .
echo .
echo .
echo .
echo ----BAM2FASTQ Initiated
~/local/bam2fastq-1.1.0/bam2fastq -f -o $SCR2/$3_$2-mito#.fastq $SCR/$3_$2.mito.bam 
echo ****BAM2FASTQ DONE.
echo .
echo .
echo .
echo .
echo .
echo .
echo ----Clipping Initiated
if [ -e "/scratch/pbsjobs/sneha/BAMSlicer_tmp/$3_$2-mito_1.fastq" ]
then
echo "PAIRED-END"
        echo "--CLIPPED: Removing first and last 2 base pairs from every read"
        ~/local/seqtk-master/seqtk trimfq -b 2 -e 2 $SCR2/$3_$2-mito_1.fastq > $SCR2/$3_$2-mito_CLIPPED_1.fastq
        ~/local/seqtk-master/seqtk trimfq -b 2 -e 2 $SCR2/$3_$2-mito_2.fastq > $SCR2/$3_$2-mito_CLIPPED_2.fastq
else
echo "SINGLE-END"
        echo "--CLIPPED: Removing first and last 2 base pairs from every read"
        ~/local/seqtk-master/seqtk trimfq -b 2 -e 2 $SCR2/$3_$2-mito.fastq > $SCR2/$3_$2-mito_CLIPPED.fastq
fi
echo ****Clipping DONE.
echo .
echo .
echo .
echo .
echo .
echo .



function CountReads {
samtools view -c -F 4 $SCR2/$3_$2$1.sorted.bam > $SCR/$3_$2$1.count
}

function Align {
if [ -e "$SCR2/$3_$2-mito_1.fastq" ]
then
bwa mem -t 12 $5 $SCR2/$3_$2$1_1.fastq $SCR2/$3_$2$1_2.fastq > $SCR2/$3_$2$4.sam
#samtools view -bS $SCR2/$3_$2$4.sam > $SCR2/$3_$2$4.bam
#samtools sort $SCR2/$3_$2$4.bam $SCR2/$3_$2$4.sorted
#samtools mpileup -B -C 50 -q 40 -Q 30 -f $5 $SCR2/$3_$2$4.sorted.bam > $SCR2/$3_$2$4.$6
else
bwa mem -t 12 $5 $SCR2/$3_$2$1.fastq > $SCR2/$3_$2$4.sam
#samtools view -bS $SCR2/$3_$2$4.sam > $SCR2/$3_$2$4.bam
#samtools sort $SCR2/$3_$2$4.bam $SCR2/$3_$2$4.sorted
#samtools mpileup -B -C 50 -q 40 -Q 30 -f $5 $SCR2/$3_$2$4.sorted.bam > $SCR2/$3_$2$4.$6
fi
}

function AlignNUMTS {
awk 'FNR==NR{a[$1]++;next}!a[$1]' $SCR2/$6 $SCR2/$3_$2$1.sam > $SCR2/$3_$2$1-$7.sam
samtools view -bS $SCR2/$3_$2$1-$7.sam > $SCR2/$3_$2$1-$7.bam
samtools sort $SCR2/$3_$2$1-$7.bam $SCR2/$3_$2$1-$7.sorted
samtools mpileup -B -C 50 -q 40 -Q 30 -f $5 $SCR2/$3_$2$1-$7.sorted.bam > $SCR/$3_$2$1-$7.$4
}

function lowNUMTS {
grep "NM:i:1" $SCR2/$3_$2$1 > $SCR2/$3_$2$4
awk '$6==length($10)"M"' $SCR2/$3_$2$4 >  $SCR2/$3_$2$4--NUMTs.sam

awk '{print $1}' $SCR2/$3_$2$5 | uniq -c | sort -k1 > $SCR2/$3_$2$4i0uniq
awk '{print $1}' $SCR2/$3_$2$4--NUMTs.sam | uniq -c | sort -k1 > $SCR2/$3_$2$4i1uniq

awk '$1==2{print $2}' $SCR2/$3_$2$4i0uniq > $SCR2/$3_$2$4i0NUMTsReads
awk '$1==2{print $2}' $SCR2/$3_$2$4i1uniq > $SCR2/$3_$2$4i1NUMTsReads

awk 'FNR==NR{a[$1]++;next}a[$1]' $SCR2/$3_$2$4i0NUMTsReads $SCR2/$3_$2$1 > $SCR2/$3_$2$4i0NUMTs.sam
awk 'FNR==NR{a[$1]++;next}a[$1]' $SCR2/$3_$2$4i1NUMTsReads $SCR2/$3_$2$1 > $SCR2/$3_$2$4i1NUMTs.sam

awk '$1==1{print $2}' $SCR2/$3_$2$4i0uniq > $SCR2/$3_$2$4i0NUMTsSingleReads
awk '$1==1{print $2}' $SCR2/$3_$2$4i1uniq > $SCR2/$3_$2$4i1NUMTsSingleReads

awk 'FNR==NR{a[$1]++;next}a[$1]' $SCR2/$3_$2$4i0NUMTsSingleReads $SCR2/$3_$2$4i1NUMTsSingleReads > $SCR2/$3_$2$4NUMTsDupsReads
awk 'FNR==NR{a[$1]++;next}a[$1]' $SCR2/$3_$2$4NUMTsDupsReads $SCR2/$3_$2$1 >  $SCR2/$3_$2$4NUMTsDups.sam

if [ -e "$SCR2/$3_$2-mito_1.fastq" ]
then
	echo PAIRED-END, removing all mate pairs which align with up to 1 mismatch to nuclear genome
	cat $SCR2/$3_$2$4i0NUMTs.sam $SCR2/$3_$2$4i1NUMTs.sam $SCR2/$3_$2$4NUMTsDups.sam > $SCR2/$3_$2$4--lowNUMTs.sam
else
	echo SINGLE-END, removing all singlet reads which align with up to 1 mismatch to nuclear genome
	cat $SCR2/$3_$2$4i0NUMTsSingleReads $SCR2/$3_$2$4i1NUMTsSingleReads > $SCR2/$3_$2$4--lowNUMTs.sam
fi
}

echo ----bwa alignment to rCRS
Align -mito_CLIPPED $2 $3 _cl--rCRS $REF/MitoREFs/rCRS.fa pileup
echo ****bwa alignment to rCRS DONE.
echo .
echo .
echo .
echo .
echo .
echo .
echo ----bwa alignment to hg38-norCRS --NUMTs: going to make $3_$2_cl--nuclear.sam
if [ -e "/scratch/pbsjobs/sneha/BAMSlicer_tmp/$3_$2-mito_CLIPPED_1.fastq" ]
then
bwa mem -t 12 $REF/hg38-nochr.fa $SCR2/$3_$2-mito_CLIPPED_1.fastq $SCR2/$3_$2-mito_CLIPPED_2.fastq > $SCR2/$3_$2_cl--nuclear.sam
else
bwa mem -t 12 $REF/hg38-nochr.fa $SCR2/$3_$2-mito_CLIPPED.fastq > $SCR2/$3_$2_cl--nuclear.sam
fi
echo ****bwa alignment to hg38-norCRS DONE.
echo .
echo .
echo .
echo .
echo .
echo .
echo ----Extract Perfect Matches to Nuclear Genome.
grep `awk 'END {print "MD:Z:"length($10)}' $SCR2/$3_$2_cl--nuclear.sam` $SCR2/$3_$2_cl--nuclear.sam > $SCR2/$3_$2_cl--NUMTs.sam
lowNUMTS _cl--nuclear.sam $2 $3 _cl-NM _cl--NUMTs.sam
echo ****Extraction of Perfect Matches to Nuclear Genome DONE.
echo .
echo .
echo .
echo .
echo .
echo .
echo ----Generate $3_$2_cl--rCRS-lowNUMTS.pileup
AlignNUMTS _cl--rCRS $2 $3 pileup $REF/MitoREFs/rCRS.fa $3_$2_cl-NM--lowNUMTs.sam lowNUMTs
CountReads _cl--rCRS-lowNUMTs $2 $3
echo ****Generate $3_$2_cl--rCRS-lowNUMTS.pileup DONE.
echo ****NUMTs DONE.
echo .
echo .
echo .
echo .
echo .
echo .
echo ----BAM2FASTQ $3_$2_cl--rCRS-lowNUMTs.sam Initiated
samtools view -bS $SCR2/$3_$2_cl--rCRS-lowNUMTs.sam > $SCR/$3_$2_cl--rCRS-lowNUMTs.bam
~/local/bam2fastq-1.1.0/bam2fastq -o $SCR2/$3_$2_cl--rCRS-lowNUMTs#.fastq $SCR/$3_$2_cl--rCRS-lowNUMTs.bam
echo ****BAM2FASTQ DONE.
echo .
echo .
echo .
echo .
echo .
echo .
echo ****bwa alignment of $3_$2_cl--rCRS-lowNUMTs fastqs to hg38
if [ -e "/scratch/pbsjobs/sneha/BAMSlicer_tmp/$3_$2_cl--rCRS-lowNUMTs_1.fastq" ]
then
echo "PAIRED-END"
	bwa aln $REF/hg38.fa $SCR2/$3_$2_cl--rCRS-lowNUMTs_1.fastq > $SCR2/$3_$2_cl--rCRS-lowNUMTs_1.sai
	bwa aln $REF/hg38.fa $SCR2/$3_$2_cl--rCRS-lowNUMTs_2.fastq > $SCR2/$3_$2_cl--rCRS-lowNUMTs_2.sai
	bwa sampe $REF/hg38.fa $SCR2/$3_$2_cl--rCRS-lowNUMTs_1.sai $SCR2/$3_$2_cl--rCRS-lowNUMTs_2.sai $SCR2/$3_$2_cl--rCRS-lowNUMTs_1.fastq $SCR2/$3_$2_cl--rCRS-lowNUMTs_2.fastq > $SCR2/$3_$2.mito_hg38.sam
else
echo "SINGLE-END"
	bwa aln $REF/hg38.fa $SCR2/$3_$2_cl--rCRS-lowNUMTs.fastq > $SCR2/$3_$2_cl--rCRS-lowNUMTs.sai
	bwa samse $REF/hg38.fa $SCR2/$3_$2_cl--rCRS-lowNUMTs.sai $SCR2/$3_$2_cl--rCRS-lowNUMTs.fastq > $SCR2/$3_$2.mito_hg38.sam
fi
echo ****bwa DONE.
echo .
echo .
echo .
echo .
echo .
echo .

echo ----samtools Started
samtools view -bS $SCR2/$3_$2.mito_hg38.sam > $SCR2/$3_$2.mito_hg38.bam
samtools sort $SCR2/$3_$2.mito_hg38.bam $SCR/$3_$2.mito_hg38.sorted
samtools index $SCR/$3_$2.mito_hg38.sorted.bam $SCR/$3_$2.mito_hg38.sorted.bam.bai
cp $SCR/$3_$2.mito_hg38.sorted.bam $SCR2/$3_$2.mito_hg38.sorted.bam
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
#mv $SCR/$3_$2_cl--rCRS-lowNUMTs.bam $SCR3/Bams/
#mv $SCR/$3_$2.mito_hg38.sorted.bam $SCR3/Bams/
#mv $SCR/$3_$2_*.count $SCR3/Counts
#mv $SCR/$3_$2_cl--rCRS-lowNUMTs.pileup $SCR3/Pileups/
#rm $SCR/$3_$2*
#rm $SCR2/$3_$2*
