#!/usr/bin/
#script.sh
#runs with pipeline.sh and frequencies.sh; .bed file optional
#central file for processing .bams
# *****always fix comments and chromosomes and I/O!!!!***********

TMPDIR=$2
#### $2 is the OUT DIR where you want the final.vcf file to be saved.

echo *****$1*****

echo *****AddOrReplaceReadGroups $1*****
java -Xmx8g -jar /home/sxg501/Tools/picard-tools-1.93/AddOrReplaceReadGroups.jar \
I=$1 \
O=$TMPDIR/$3.tcga.sort.2.bam \
SORT_ORDER=coordinate \
RGID=TLL \
RGLB=bar \
RGPL=illumina \
RGSM=foo \
RGPU=bc \
VALIDATION_STRINGENCY=SILENT 
echo *****done*****

echo *****Picard Mark Duplicates and RealignerTargetCreator $1*****
touch /home/sxg501/GATK-ReAnalyze/VAF-nofix-tmpfiles/$3.metricsfile.txt

java -Xmx8g -jar /home/sxg501/Tools/picard-tools-1.93/MarkDuplicates.jar \
INPUT=$TMPDIR/$3.tcga.sort.2.bam \
OUTPUT=$TMPDIR/$3.tcga.marked.bam \
METRICS_FILE=$TMPDIR/$3.metricsfile.txt \
CREATE_INDEX=true \
VALIDATION_STRINGENCY=SILENT \
ASSUME_SORTED=true

java -Xmx8g -jar /home/sxg501/Tools/GenomeAnalysisTK.jar \
-T RealignerTargetCreator \
-R /home/sxg501/REFs/MitoREFs/rCRS.fa \
-I $TMPDIR/$3.tcga.marked.bam \
-o $TMPDIR/$3.tcga.list
echo *****done*****

echo *****IndelRealigner and Picard Marking $1*****
java -Xmx8g -jar /home/sxg501/Tools/GenomeAnalysisTK.jar \
-T IndelRealigner \
-I $TMPDIR/$3.tcga.marked.bam \
-R /home/sxg501/REFs/MitoREFs/rCRS.fa \
--targetIntervals $TMPDIR/$3.tcga.list \
-o $TMPDIR/$3.tcga.marked.realigned.bam

java -Xmx8g -jar /home/sxg501/Tools/picard-tools-1.93/FixMateInformation.jar \
INPUT=$TMPDIR/$3.tcga.marked.realigned.bam \
OUTPUT=$TMPDIR/$3.tcga.marked.realigned.fixed.bam \
SO=coordinate \
VALIDATION_STRINGENCY=SILENT \
CREATE_INDEX=true
echo *****done*****

echo *****BaseRecalibrator $1*****
java -Xmx8g -jar /home/sxg501/Tools/GenomeAnalysisTK.jar \
-T BaseRecalibrator \
-I $TMPDIR/$3.tcga.marked.realigned.fixed.bam \
-R /home/sxg501/REFs/MitoREFs/rCRS.fa \
-o $TMPDIR/$3.recal_data.grp \
--knownSites /home/sxg501/GATK-ReAnalyze/dbsnp/mtONLY.vcf 
echo *****done*****

echo *****PrintReads $1*****
java -Xmx8g -jar /home/sxg501/Tools/GenomeAnalysisTK.jar \
-T PrintReads \
-R /home/sxg501/REFs/MitoREFs/rCRS.fa \
-I $TMPDIR/$3.tcga.marked.realigned.fixed.bam \
-BQSR $TMPDIR/$3.recal_data.grp \
-o $TMPDIR/$3.tcga.marked.realigned.fixed.read.bam
echo *****done*****

echo *****HaplotypeCaller $1*****
java -Xmx46g -jar /home/sxg501/Tools/GenomeAnalysisTK.jar \
-T HaplotypeCaller \
-R /home/sxg501/REFs/MitoREFs/rCRS.fa \
-I $TMPDIR/$3.tcga.marked.realigned.fixed.read.bam \
--maxReadsInRegionPerSample 200 \
--sample_ploidy $4 \
-stand_call_conf 50 \
-stand_emit_conf 10 \
--pcr_indel_model HOSTILE \
-minPruning 10 \
-A StrandAlleleCountsBySample \
--dbsnp /home/sxg501/GATK-ReAnalyze/dbsnp/mtONLY.vcf \
-o $TMPDIR/$3.tcga.snps.vcf
echo *****done*****

#echo *****Starting DepthOfCoverage*****
#java -Xmx64g -jar /homeG/LaFramboiseLab/sjm67/GenomeAnalysisTK-2.4-9-g532efad/GenomeAnalysisTK.jar  \
#-R /homeG/LaFramboiseLab/sjm67/TCGA/hg19_karyo_rCRS/hg19_karyo_rCRS.fa \
#-T Coverage \
#-o $1.out.table \
#-I $1.tcga.marked.realigned.fixed.read.bam \
#--outputFormat table
#echo *****done*****

echo *****VariantFiltration $1******
java -Xmx8g -jar /home/sxg501/Tools/GenomeAnalysisTK.jar \
-R /home/sxg501/REFs/MitoREFs/rCRS.fa \
-T VariantFiltration \
--variant $TMPDIR/$3.tcga.snps.vcf \
-o $TMPDIR/$3.snps.recalibrated.filtered.vcf \
--filterExpression "QD < 2.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0" \
--filterName "my_snp_filter" \

#rm paired*
#rm tcga*
#rm metrics*
#rm snps*
#rm extract*
#rm recal*
#echo *****done*****

echo *****Creating Report $1*****
#perl /home/sxg501/GATK-ReAnalyze/frequencies.pl $1.snps.recalibrated.filtered.vcf $1.out.txt
echo *****done*****
#cp $TMPDIR/$3.snps.recalibrated.filtered.vcf $2/
echo *****Creating Report $1*****
#perl /home/sxg501/GATK-ReAnalyze/frequencies.pl $1.snps.recalibrated.filtered.vcf $1.out.txt
echo *****done*****
#cp $TMPDIR/$3.snps.recalibrated.filtered.vcf $2/
#echo ***making pileup file****
#samtools mpileup -f ~/REFs/MitoREFs/rCRS.fa $TMPDIR/$3.tcga.marked.realigned.fixed.read.bam > ~/GATK-ReAnalyze/Pileups/$3.pileup
#perl ~/GATK-ReAnalyze/Pileups/pileupstats2.pl ~/GATK-ReAnalyze/Pileups/$3.pileup ~/GATK-ReAnalyze/Pileups/$3.pOUT
#bash ~/GATK-ReAnalyze/Pileups/pOUT-mod.sh ~/GATK-ReAnalyze/Pileups/$3.pOUT

echo *****done $1*****
