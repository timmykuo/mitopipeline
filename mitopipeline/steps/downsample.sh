##########written by Sneha Grandhi
#### 8.25.15
#### takes hg38_aligned mito .bam files from TCGA (4155 files)
#### and converts them to the sanger ASCII encoding if they come from a version < Illumina 1.8
#### runs files through GATK pipeline 
#### produces accurate .vcf file
#### during this process it treats single and paired end data differently as well.
######################
#### $1 = original .bam file (with location)
####      /scratch/pbsjobs/sneha/Bams/BLCA_A0C8-01.mito_hg38.sorted.bam
#### $2 = filename unique identifier
####      BLCA_A0C8-01.mito_hg38
#### $3 = directory where you want the fixed ASCII 64 std OUT and tmp realignment files:
####      /scratch/pbsjobs/sneha/ASCII-scaling-GATK
#######################################

FILE=$1
START=$2
OUT=$3
TOOLS=$4
echo "*****Patient:$1"
echo .
echo .
echo .
$TOOLS/bam2fastq-1.1.0/bam2fastq -o $OUT/$1#.fastq $2/$1

if [ -e "$OUT/$1_1.fastq" ]
then
	echo "SAMPLE IS PAIRED-END"
	perl $TOOLS/DetermineFastqQualityEncoding.pl $OUT/$1_1.fastq > $OUT/$1_ENCODING.txt
	grep "offset by 64" $OUT/$1_ENCODING.txt > $OUT/$1_ENCODING_OLD

	if [ -s "$OUT/$1_ENCODING_OLD" ]
	then
		echo "*****This base encoding is in ASCII +64"
		echo "*****Fixing the base encoding with seqtk"
		$TOOLS/seqtk-master/seqtk seq -Q64 -V $OUT/$1_1.fastq > $OUT/$1_sanger_1.fastq
		$TOOLS/seqtk-master/seqtk seq -Q64 -V $OUT/$1_2.fastq > $OUT/$1_sanger_2.fastq
		echo "*****BWA Alignment to rCRS of fixed encoding fastq files"
		$TOOLS/bwa-0.6.0/bwa aln ../REFs/MitoREFs/rCRS.fa $OUT/$1_sanger_1.fastq > $OUT/$1_1.sai
		$TOOLS/bwa-0.6.0/bwa aln ../REFs/MitoREFs/rCRS.fa $OUT/$1_sanger_2.fastq > $OUT/$1_2.sai
		$TOOLS/bwa-0.6.0/bwa sampe ../REFs/MitoREFs/rCRS.fa $OUT/$1_1.sai $OUT/$1_2.sai $OUT/$1_sanger_1.fastq $OUT/$1_sanger_2.fastq > $OUT/$1_realigned.sam
		echo "*****Samtools to set up the sorted and indexed .bam file"
		samtools view -bS $OUT/$1_realigned.sam > $OUT/$1_realigned.bam
		samtools sort $OUT/$1_realigned.bam $OUT/$1_realigned.sorted
		samtools index $OUT/$1_realigned.sorted.bam
		
		java -Xmx40g -jar $TOOLS/picard-tools-1.93/DownsampleSam.jar INPUT=$OUT/$1_realigned.sorted.bam OUTPUT=$OUT/$1_realigned-0.5.bam PROBABILITY=0.5 VALIDATION_STRINGENCY=SILENT
		echo .	
		echo .
		echo .
	else
		echo "*****This base encoding is in ASCII 33 / Sanger"
		echo "*****Continuing with GATK"
		echo "*****BWA Alignment to rCRS of fixed encoding fastq files"
                $TOOLS/bwa-0.6.0/bwa aln ../REFs/MitoREFs/rCRS.fa $OUT/$1_1.fastq > $OUT/$1_1.sai
                $TOOLS/bwa-0.6.0/bwa aln ../REFs/MitoREFs/rCRS.fa $OUT/$1_2.fastq > $OUT/$1_2.sai
                $TOOLS/bwa-0.6.0/bwa sampe ../REFs/MitoREFs/rCRS.fa $OUT/$1_1.sai $OUT/$1_2.sai $OUT/$1_1.fastq $OUT/$1_2.fastq > $OUT/$1_realigned.sam
                echo "*****Samtools to set up the sorted and indexed .bam file"
                samtools view -bS $OUT/$1_realigned.sam > $OUT/$1_realigned.bam
                samtools sort $OUT/$1_realigned.bam $OUT/$1_realigned.sorted
                samtools index $OUT/$1_realigned.sorted.bam
		
		java -Xmx40g -jar $TOOLS/picard-tools-1.93/DownsampleSam.jar INPUT=$OUT/$1_realigned.sorted.bam OUTPUT=$OUT/$1_realigned-0.5.bam PROBABILITY=0.5 VALIDATION_STRINGENCY=SILENT
        fi

else
	echo "SAMPLE IS SINGLE-END"
        perl $TOOLS/DetermineFastqQualityEncoding.pl $OUT/$1.fastq > $OUT/$1_ENCODING.txt
        grep "offset by 64" $OUT/$1_ENCODING.txt > $OUT/$1_ENCODING_OLD
	
	if [ -s "$OUT/$1_ENCODING_OLD" ]
        then
                echo "*****This base encoding is in ASCII +64"
                echo "*****Fixing the base encoding with seqtk"
                $TOOLS/seqtk-master/seqtk seq -Q64 -V $OUT/$1.fastq > $OUT/$1_sanger.fastq
                echo "*****BWA Alignment to rCRS of fixed encoding fastq files"
                $TOOLS/bwa-0.6.0/bwa aln ../REFs/MitoREFs/rCRS.fa $OUT/$1_sanger.fastq > $OUT/$1.sai
                $TOOLS/bwa-0.6.0/bwa samse ../REFs/MitoREFs/rCRS.fa $OUT/$1.sai $OUT/$1_sanger.fastq > $OUT/$1_realigned.sam
                echo "*****Samtools to set up the sorted and indexed .bam file"
                samtools view -bS $OUT/$1_realigned.sam > $OUT/$1_realigned.bam
                samtools sort $OUT/$1_realigned.bam $OUT/$1_realigned.sorted
                samtools index $OUT/$1_realigned.sorted.bam
        	java -Xmx40g -jar $TOOLS/picard-tools-1.93/DownsampleSam.jar INPUT=$OUT/$1_realigned.sorted.bam OUTPUT=$OUT/$1_realigned-0.5.bam PROBABILITY=0.5 VALIDATION_STRINGENCY=SILENT
                echo .
                echo .
                echo .
        else
                echo "*****This base encoding is in ASCII 33 / Sanger"
		echo "*****BWA Alignment to rCRS of fixed encoding fastq files"
                $TOOLS/bwa-0.6.0/bwa aln ../REFs/MitoREFs/rCRS.fa $OUT/$1.fastq > $OUT/$1.sai
                $TOOLS/bwa-0.6.0/bwa samse ../REFs/MitoREFs/rCRS.fa $OUT/$1.sai $OUT/$1.fastq > $OUT/$1_realigned.sam
                echo "*****Samtools to set up the sorted and indexed .bam file"
                samtools view -bS $OUT/$1_realigned.sam > $OUT/$1_realigned.bam
                samtools sort $OUT/$1_realigned.bam $OUT/$1_realigned.sorted
                samtools index $OUT/$1_realigned.sorted.bam
		
		java -Xmx40g -jar $TOOLS/picard-tools-1.93/DownsampleSam.jar INPUT=$OUT/$1_realigned.sorted.bam OUTPUT=$OUT/$1_realigned-0.5.bam PROBABILITY=0.5 VALIDATION_STRINGENCY=SILENT
        fi	
fi
rm $OUT/$1*fastq
#rm $OUT/$1*bam*
rm $OUT/$1*sam*
rm $OUT/$1*sai
