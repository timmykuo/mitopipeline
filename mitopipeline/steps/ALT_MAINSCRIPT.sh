### $4 is the file size
### $3 is the Cancer Type
### $2 is the Target file ID
### $1 is the download ID

OUTdir="/mnt/projects/SOM_GEN_TXL80/users/sxg501/TARGET"
STOR="/scratch/pbsjobs/sneha3/Bams"
SCR="/scratch/pbsjobs/sneha/mitoBAMs"
SCR2="/scratch/pbsjobs/sneha/mitoTEMP"
HOMEdir="/home/sxg501/TARGET/ThroughGDC"
REF="/home/sxg501/REFs/"

#module load base
#module del python
module load samtools
module load bwa
module load gdc/1.3.0

echo creating temporary manifest file 
echo cancer: $3
echo file: $2
echo downloadID: $1
echo .
echo .
	grep $1 $HOMEdir/gdc_manifest.2017-08-08T19-27-28.292690.tsv | cat $HOMEdir/manifest_header - > $HOMEdir/$2_$3manifest_tmp


echo downloading from gdc
	cd $STOR
	gdc-client download -n 25 -m $HOMEdir/$2_$3manifest_tmp -t /home/sxg501/gdc_tokens/gdc-user-token.2017.09.txt
 
rm $HOMEdir/$2_$3manifest_tmp
echo .
echo .
echo download finished
	#cd $STOR/$1/
	#samtools index *bam 
echo .
echo .
echo extracting chrM bam file.
	mkdir $STOR/$1/PluckMito/
        samtools view -H $STOR/$1/*.bam > $STOR/$1/PluckMito/header.bam
        cd $STOR/$1/PluckMito/
        grep -i SN:MT $STOR/$1/PluckMito/header.bam > $STOR/$1/PluckMito/MT.bam
        grep -i SN:chrM_rCRS $STOR/$1/PluckMito/header.bam > $STOR/$1/PluckMito/chrM_rCRS.bam
        grep -i SN:chrM $STOR/$1/PluckMito/header.bam > $STOR/$1/PluckMito/chrM.bam
        grep -i SN:M $STOR/$1/PluckMito/header.bam > $STOR/$1/PluckMito/M.bam

        cd $STOR/$1/PluckMito/
        if [ -s "MT.bam" ]
        then
                echo "The mitochondrial genome is annotated by 'MT'."
                samtools view -b $STOR/$1/*.bam MT > $SCR/$3_$2_mito.bam

        elif [ -s "chrM_rCRS.bam" ]
        then
                echo "The mitochondrial genome is annotated by 'chrM_rCRS'."
                samtools view -b $STOR/$1/*.bam chrM_rCRS > $SCR/$3_$2_mito.bam

        elif [ -s "chrM.bam" ]
        then
                echo "The mitochondrial genome is annotated by 'chrM'."
                samtools view -b $STOR/$1/*.bam chrM > $SCR/$3_$2_mito.bam
	elif [ -s "M.bam" ]
        then
                echo "The mitochondrial genome is annotated by 'M'."
                samtools view -b $STOR/$1/*.bam M > $SCR/$3_$2_mito.bam
        fi


rm -r $STOR/$1


echo ----bam2fastq
/home/sxg501/local/bam2fastq-1.1.0/bam2fastq -f -o $SCR2/$3_$2_mito#.fastq $SCR/$3_$2_mito.bam 
echo ****BAM2FASTQ DONE.
echo .
echo .
echo .
echo .
echo .
echo .
echo ----Clipping Initiated
if [ -e "$SCR2/$3_$2_mito_1.fastq" ]
then
echo "PAIRED-END"
        echo "--CLIPPED: Removing first and last 2 base pairs from every read"
        /home/sxg501/local/seqtk-master/seqtk trimfq -b 2 -e 2 $SCR2/$3_$2_mito_1.fastq > $SCR2/$3_$2_mito_CLIPPED_1.fastq
        /home/sxg501/local/seqtk-master/seqtk trimfq -b 2 -e 2 $SCR2/$3_$2_mito_2.fastq > $SCR2/$3_$2_mito_CLIPPED_2.fastq
else
echo "SINGLE-END"
        echo "--CLIPPED: Removing first and last 2 base pairs from every read"
        /home/sxg501/local/seqtk-master/seqtk trimfq -b 2 -e 2 $SCR2/$3_$2_mito.fastq > $SCR2/$3_$2_mito_CLIPPED.fastq
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
if [ -e "$SCR2/$3_$2_mito_1.fastq" ]
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

if [ -e "$SCR2/$3_$2_mito_1.fastq" ]
then
	echo PAIRED-END, removing all mate pairs which align with up to 1 mismatch to nuclear genome
	cat $SCR2/$3_$2$4i0NUMTs.sam $SCR2/$3_$2$4i1NUMTs.sam $SCR2/$3_$2$4NUMTsDups.sam > $SCR2/$3_$2$4--lowNUMTs.sam
else
	echo SINGLE-END, removing all singlet reads which align with up to 1 mismatch to nuclear genome
	cat $SCR2/$3_$2$4i0NUMTsSingleReads $SCR2/$3_$2$4i1NUMTsSingleReads > $SCR2/$3_$2$4--lowNUMTs.sam
fi
}

echo ----bwa alignment to rCRS
Align _mito_CLIPPED $2 $3 _cl--rCRS $REF/MitoREFs/rCRS.fa pileup
echo ****bwa alignment to rCRS DONE.
echo .
echo .
echo .
echo .
echo .
echo .
echo ----bwa alignment to hg38-norCRS --NUMTs: going to make $3_$2_cl--nuclear.sam
if [ -e "$SCR2//$3_$2_mito_CLIPPED_1.fastq" ]
then
bwa mem -t 12 $REF/hg38-nochr.fa $SCR2/$3_$2_mito_CLIPPED_1.fastq $SCR2/$3_$2_mito_CLIPPED_2.fastq > $SCR2/$3_$2_cl--nuclear.sam
else
bwa mem -t 12 $REF/hg38-nochr.fa $SCR2/$3_$2_mito_CLIPPED.fastq > $SCR2/$3_$2_cl--nuclear.sam
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
/home/sxg501/local/bam2fastq-1.1.0/bam2fastq -o $SCR2/$3_$2_cl--rCRS-lowNUMTs#.fastq $SCR/$3_$2_cl--rCRS-lowNUMTs.bam
echo ****BAM2FASTQ DONE.
echo .
echo .
echo .
echo .
echo .
echo .
echo ****bwa alignment of $3_$2_cl--rCRS-lowNUMTs fastqs to hg38
if [ -e "$SCR2/$3_$2_cl--rCRS-lowNUMTs_1.fastq" ]
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
samtools index $SCR/$3_$2.mito_hg38.sorted.bam
cp $SCR/$3_$2.mito_hg38.sorted.bam $SCR3/Bams/$3_$2.mito_hg38.sorted.bam
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
mv $SCR/$3_$2.mito_hg38.sorted.bam $OUTdir/gdc-mito-bams/
mv $SCR/$3_$2.mito_hg38.sorted.bam.bai $OUTdir/gdc-mito-bams/
mv $SCR/$3_$2_cl--rCRS-lowNUMTs.pileup $OUTdir/gdc-mito-pileups/
rm $SCR/$3_$2*
rm $SCR2/$3_$2*
