#!/usr/bin/
#script.sh
#runs with pipeline.sh and frequencies.sh; .bed file optional
#central file for processing .bams
# *****always fix comments and chromosomes and I/O!!!!***********
###############
#### $1 - bam file to process along with file path
#### $2 - directory to store vcfs
#### $3 - filename unique identifier ex. BLCA_A0C8_01.mito_hg38
#### $4 - mtDNA Copy Number (used in GATK) use 100 if you can't find it on gdc portal
################

#$1 - filename
#$2 - start directory
#$3 - out directory

module load samtools
TMPDIR=$3"/temp"
BAMS=$2"/no_numts_bams"
REFS="./REFs"
#### $2 is the OUT DIR where you want the final.vcf file to be saved.
echo *****$1*****

echo *****AddOrReplaceReadGroups $1*****
java -Xmx8g -jar ../tools/picard-tools-1.93/AddOrReplaceReadGroups.jar \
I=$BAMS/$1.mito_hg38.sorted.mt.bam \
O=$TMPDIR/$1.tcga.sort.2.bam \
SORT_ORDER=coordinate \
RGID=TLL \
RGLB=bar \
RGPL=illumina \
RGSM=foo \
RGPU=bc \
VALIDATION_STRINGENCY=SILENT 
echo *****done*****

echo *****Picard Mark Duplicates and RealignerTargetCreator $1*****
touch $TMPDIR/$1.metricsfile.txt

java -Xmx8g -jar ../tools/picard-tools-1.93/MarkDuplicates.jar \
INPUT=$TMPDIR/$1.tcga.sort.2.bam \
OUTPUT=$TMPDIR/$1.tcga.marked.bam \
METRICS_FILE=$TMPDIR/$1.metricsfile.txt \
CREATE_INDEX=true \
VALIDATION_STRINGENCY=SILENT \
ASSUME_SORTED=true

java -Xmx8g -jar ../tools/gatk/GenomeAnalysisTK.jar \
-T RealignerTargetCreator \
-R $REFS/MitoREFs/rCRS-MT.fa \
-I $TMPDIR/$1.tcga.marked.bam \
-o $TMPDIR/$1.tcga.list
echo *****done*****

echo *****IndelRealigner and Picard Marking $1*****
java -Xmx8g -jar ../tools/gatk/GenomeAnalysisTK.jar \
-T IndelRealigner \
-I $TMPDIR/$1.tcga.marked.bam \
-R $REFS/MitoREFs/rCRS-MT.fa \
--targetIntervals $TMPDIR/$1.tcga.list \
-o $TMPDIR/$1.tcga.marked.realigned.bam

java -Xmx8g -jar ../tools/picard-tools-1.93/FixMateInformation.jar \
INPUT=$TMPDIR/$1.tcga.marked.realigned.bam \
OUTPUT=$TMPDIR/$1.tcga.marked.realigned.fixed.bam \
SO=coordinate \
VALIDATION_STRINGENCY=SILENT \
CREATE_INDEX=true
echo *****done*****

echo *****BaseRecalibrator $1*****
java -Xmx8g -jar ../tools/gatk/GenomeAnalysisTK.jar \
-T BaseRecalibrator \
-I $TMPDIR/$1.tcga.marked.realigned.fixed.bam \
-R $REFS/MitoREFs/rCRS-MT.fa \
-o $TMPDIR/$1.recal_data.grp \
--knownSites ./dbsnp/mtONLY.vcf 
echo *****done*****

echo *****PrintReads $1*****
java -Xmx8g -jar ../tools/gatk/GenomeAnalysisTK.jar \
-T PrintReads \
-R $REFS/MitoREFs/rCRS-MT.fa \
-I $TMPDIR/$1.tcga.marked.realigned.fixed.bam \
-BQSR $TMPDIR/$1.recal_data.grp \
-o $TMPDIR/$1.tcga.marked.realigned.fixed.read.bam
echo *****done*****

echo *****HaplotypeCaller $1*****
java -Xmx10g -jar ../tools/gatk/GenomeAnalysisTK.jar \
-T HaplotypeCaller \
-R $REFS/MitoREFs/rCRS-MT.fa \
-I $TMPDIR/$1.tcga.marked.realigned.fixed.read.bam \
--maxReadsInRegionPerSample 200 \
--sample_ploidy 100 \
-stand_call_conf 50 \
-stand_emit_conf 10 \
--pcr_indel_model HOSTILE \
-minPruning 10 \
-A StrandAlleleCountsBySample \
--dbsnp ./dbsnp/mtONLY.vcf \
-mmq 0 \
-drf DuplicateRead \
-o $TMPDIR/$1.tcga.snps.vcf
echo *****done*****

echo *****VariantFiltration $1******
java -Xmx8g -jar ./.tools/gatk/GenomeAnalysisTK.jar \
-R $REFS/MitoREFs/rCRS-MT.fa \
-T VariantFiltration \
--variant $TMPDIR/$1.tcga.snps.vcf \
-o $TMPDIR/$1.snps.recalibrated.filtered.vcf \
--filterExpression "QD < 2.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0" \
--filterName "my_snp_filter" \

#rm paired*
#rm tcga*
#rm metrics*
#rm snps*
#rm extract*
#rm recal*
#echo *****done*****
mv $TMPDIR/$1.snps.recalibrated.filtered.vcf $3
echo *****done $1*****
