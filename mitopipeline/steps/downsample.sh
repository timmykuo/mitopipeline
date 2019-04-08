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
#### $3 = directory where you want the fixed ASCII 64 std OUT and tmp realignment files : 
####      /scratch/pbsjobs/sneha/ASCII-scaling-GATK
#### $4 = directory where you want the originally ASCII 33 std OUT files
####      /scratch/pbsjobs/sneha/ASCII-33-safe
#### $5 = mtDNA Copy Number (used in GATK)
#### $6 = File Number.
#### output appended to seqtk1.out / seqtk2.out / etc.,
#######################################

FILE=$1
START=$2
OUT=$3
echo "*****Analyzing File No. $6*****"
echo "*****Patient:$2"
echo .
echo .
echo .
/home/sxg501/local/bam2fastq-1.1.0/bam2fastq -o $3/$2#.fastq $1

if [ -e "$OUT/$2_1.fastq" ]
then
	echo "SAMPLE IS PAIRED-END"
	perl /home/sxg501/Tools/DetermineFastqQualityEncoding.pl $OUT/$2_1.fastq > $OUT/$2_ENCODING.txt
	grep "offset by 64" $OUT/$2_ENCODING.txt > $OUT/$2_ENCODING_OLD

	if [ -s "$OUT/$2_ENCODING_OLD" ]
	then
		echo "*****This base encoding is in ASCII +64"
		echo "*****Fixing the base encoding with seqtk"
		/home/sxg501/local/seqtk-master/seqtk seq -Q64 -V $OUT/$2_1.fastq > $OUT/$2_sanger_1.fastq
		/home/sxg501/local/seqtk-master/seqtk seq -Q64 -V $OUT/$2_2.fastq > $OUT/$2_sanger_2.fastq
		echo "*****BWA Alignment to rCRS of fixed encoding fastq files"
		/home/sxg501/Tools/bwa-0.6.0/bwa aln /home/sxg501/REFs/MitoREFs/rCRS.fa $OUT/$2_sanger_1.fastq > $OUT/$2_1.sai
		/home/sxg501/Tools/bwa-0.6.0/bwa aln /home/sxg501/REFs/MitoREFs/rCRS.fa $OUT/$2_sanger_2.fastq > $OUT/$2_2.sai
		/home/sxg501/Tools/bwa-0.6.0/bwa sampe /home/sxg501/REFs/MitoREFs/rCRS.fa $OUT/$2_1.sai $OUT/$2_2.sai $OUT/$2_sanger_1.fastq $OUT/$2_sanger_2.fastq > $OUT/$2.realigned.sam
		echo "*****Samtools to set up the sorted and indexed .bam file"
		samtools view -bS $OUT/$2.realigned.sam > $OUT/$2.realigned.bam
		samtools sort $OUT/$2.realigned.bam $OUT/$2.realigned.sorted
		samtools index $OUT/$2.realigned.sorted.bam
		
		java -Xmx40g -jar ~/Tools/picard-tools-1.93/DownsampleSam.jar INPUT=$OUT/$2.realigned.sorted.bam OUTPUT=$OUT/$2.realigned-0.5.bam PROBABILITY=0.5 VALIDATION_STRINGENCY=SILENT
		echo "*****Continuing with GATK"
		bash /home/sxg501/GATK-ReAnalyze/script-FINAL.sh $OUT/$2.realigned-0.5.bam /scratch/pbsjobs/sneha/VCFs-fix $2 $5 > $OUT/$2.GATK.OUT 2>&1
		echo .	
		echo .
		echo .
	else
		echo "*****This base encoding is in ASCII 33 / Sanger"
		echo "*****Continuing with GATK"
		echo "*****BWA Alignment to rCRS of fixed encoding fastq files"
                /home/sxg501/Tools/bwa-0.6.0/bwa aln /home/sxg501/REFs/MitoREFs/rCRS.fa $OUT/$2_1.fastq > $OUT/$2_1.sai
                /home/sxg501/Tools/bwa-0.6.0/bwa aln /home/sxg501/REFs/MitoREFs/rCRS.fa $OUT/$2_2.fastq > $OUT/$2_2.sai
                /home/sxg501/Tools/bwa-0.6.0/bwa sampe /home/sxg501/REFs/MitoREFs/rCRS.fa $OUT/$2_1.sai $OUT/$2_2.sai $OUT/$2_1.fastq $OUT/$2_2.fastq > $OUT/$2.realigned.sam
                echo "*****Samtools to set up the sorted and indexed .bam file"
                samtools view -bS $OUT/$2.realigned.sam > $OUT/$2.realigned.bam
                samtools sort $OUT/$2.realigned.bam $OUT/$2.realigned.sorted
                samtools index $OUT/$2.realigned.sorted.bam
		
		java -Xmx40g -jar ~/Tools/picard-tools-1.93/DownsampleSam.jar INPUT=$OUT/$2.realigned.sorted.bam OUTPUT=$OUT/$2.realigned-0.5.bam PROBABILITY=0.5 VALIDATION_STRINGENCY=SILENT
                echo "*****Continuing with GATK"
	        bash /home/sxg501/GATK-ReAnalyze/script-FINAL.sh $OUT/$2.realigned-0.5.bam /scratch/pbsjobs/sneha/VCFs-nofix $2 $5 > $4/$2.GATK.OUT 2>&1
	fi

else
	echo "SAMPLE IS SINGLE-END"
        perl /home/sxg501/Tools/DetermineFastqQualityEncoding.pl $OUT/$2.fastq > $OUT/$2_ENCODING.txt
        grep "offset by 64" $OUT/$2_ENCODING.txt > $OUT/$2_ENCODING_OLD
	
	if [ -s "$OUT/$2_ENCODING_OLD" ]
        then
                echo "*****This base encoding is in ASCII +64"
                echo "*****Fixing the base encoding with seqtk"
                /home/sxg501/local/seqtk-master/seqtk seq -Q64 -V $OUT/$2.fastq > $OUT/$2_sanger.fastq
                echo "*****BWA Alignment to rCRS of fixed encoding fastq files"
                /home/sxg501/Tools/bwa-0.6.0/bwa aln /home/sxg501/REFs/MitoREFs/rCRS.fa $OUT/$2_sanger.fastq > $OUT/$2.sai
                /home/sxg501/Tools/bwa-0.6.0/bwa samse /home/sxg501/REFs/MitoREFs/rCRS.fa $OUT/$2.sai $OUT/$2_sanger.fastq > $OUT/$2.realigned.sam
                echo "*****Samtools to set up the sorted and indexed .bam file"
                samtools view -bS $OUT/$2.realigned.sam > $OUT/$2.realigned.bam
                samtools sort $OUT/$2.realigned.bam $OUT/$2.realigned.sorted
                samtools index $OUT/$2.realigned.sorted.bam
        	java -Xmx40g -jar ~/Tools/picard-tools-1.93/DownsampleSam.jar INPUT=$OUT/$2.realigned.sorted.bam OUTPUT=$OUT/$2.realigned-0.5.bam PROBABILITY=0.5 VALIDATION_STRINGENCY=SILENT
                echo "*****Continuing with GATK"
	        bash /home/sxg501/GATK-ReAnalyze/script-FINAL.sh $OUT/$2.realigned-0.5.bam /scratch/pbsjobs/sneha/VCFs-fix $2 $5 > $OUT/$2.GATK.OUT 2>&1
                echo .
                echo .
                echo .
        else
                echo "*****This base encoding is in ASCII 33 / Sanger"
		echo "*****BWA Alignment to rCRS of fixed encoding fastq files"
                /home/sxg501/Tools/bwa-0.6.0/bwa aln /home/sxg501/REFs/MitoREFs/rCRS.fa $OUT/$2.fastq > $OUT/$2.sai
                /home/sxg501/Tools/bwa-0.6.0/bwa samse /home/sxg501/REFs/MitoREFs/rCRS.fa $OUT/$2.sai $OUT/$2.fastq > $OUT/$2.realigned.sam
                echo "*****Samtools to set up the sorted and indexed .bam file"
                samtools view -bS $OUT/$2.realigned.sam > $OUT/$2.realigned.bam
                samtools sort $OUT/$2.realigned.bam $OUT/$2.realigned.sorted
                samtools index $OUT/$2.realigned.sorted.bam
		
		java -Xmx40g -jar ~/Tools/picard-tools-1.93/DownsampleSam.jar INPUT=$OUT/$2.realigned.sorted.bam OUTPUT=$OUT/$2.realigned-0.5.bam PROBABILITY=0.5 VALIDATION_STRINGENCY=SILENT
                echo "*****Continuing with GATK"
                bash /home/sxg501/GATK-ReAnalyze/script-FINAL.sh $OUT/$2.realigned-0.5.bam /scratch/pbsjobs/sneha/VCFs-nofix $2 $5 > $4/$2.GATK.OUT 2>&1
        fi	
fi
rm $OUT/$2*fastq
rm $OUT/$2*bam*
rm $OUT/$2*sam*
rm $OUT/$2*sai
