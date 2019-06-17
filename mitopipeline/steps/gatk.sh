#!/bin/bash
#script.sh
###############
################
#$1 - filename
#$2 - start directory
#$3 - out directory
#$4 - tools directory
#$5 is steps directory
#$6 is the refs directory
##############
TMPDIR=$3"/gatk_stor"
BAMS=$2
REFS=$6
TOOLS=$4
STEPS=$5
#### $2 is the OUT DIR where you want the final.vcf file to be saved.
#only set filetype if the value was a step in the pipeline
#last string following / delimeter will be name of the previous job
filetype="_"$(awk -F/ '{print $NF}' <<< "$2" | awk '{print tolower($0)}')
if [ "$filetype" != "_extractmito" ] && [ "$filetype" != "_clipping" ] && [ "$filetype" != "_splitgap" ] && [ "$filetype" != "_removenumts" ];
then
filetype=""
fi
echo *****$1*****

echo *****AddOrReplaceReadGroups $1*****
java -Xmx8g -jar $STEPS/tools/AddOrReplaceReadGroups.jar \
I=$BAMS/$1${filetype}.bam \
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

java -Xmx8g -jar $STEPS/tools/MarkDuplicates.jar \
INPUT=$TMPDIR/$1.tcga.sort.2.bam \
OUTPUT=$TMPDIR/$1.tcga.marked.bam \
METRICS_FILE=$TMPDIR/$1.metricsfile.txt \
CREATE_INDEX=true \
VALIDATION_STRINGENCY=SILENT \
ASSUME_SORTED=true

java -Xmx8g -jar $STEPS/tools/FixMateInformation.jar \
INPUT=$TMPDIR/$1.tcga.marked.bam \
OUTPUT=$TMPDIR/$1.tcga.marked.realigned.fixed.bam \
SO=coordinate \
VALIDATION_STRINGENCY=SILENT \
CREATE_INDEX=true
echo *****done*****

echo *****BaseRecalibrator $1*****
java -Xmx8g -jar $TOOLS/gatk/GenomeAnalysisTK.jar BaseRecalibrator \
-I $TMPDIR/$1.tcga.marked.realigned.fixed.bam \
-R $REFS/rCRS-MT.fa \
-O $TMPDIR/$1.recal_data.grp \
--known-sites $5/dbsnp/mtONLY.vcf 
echo *****done*****

echo *****ApplyBQSR $1*****
java -Xmx8g -jar $TOOLS/gatk/GenomeAnalysisTK.jar ApplyBQSR \
-R $REFS/rCRS-MT.fa \
-I $TMPDIR/$1.tcga.marked.realigned.fixed.bam \
--bqsr-recal-file $TMPDIR/$1.recal_data.grp \
-O $TMPDIR/$1.tcga.marked.realigned.fixed.read.bam
echo *****done*****

echo *****HaplotypeCaller $1*****
java -Xmx10g -jar $TOOLS/gatk/GenomeAnalysisTK.jar HaplotypeCaller \
-R $REFS/rCRS-MT.fa \
-I $TMPDIR/$1.tcga.marked.realigned.fixed.read.bam \
--sample-ploidy 100 \
--stand-call-conf 50 \
--pcr-indel-model HOSTILE \
--min-pruning 10 \
--dbsnp $5/dbsnp/mtONLY.vcf \
-O $TMPDIR/$1.tcga.snps.vcf
echo *****done*****

echo *****VariantFiltration $1******
java -Xmx8g -jar $TOOLS/gatk/GenomeAnalysisTK.jar VariantFiltration \
-R $REFS/rCRS-MT.fa \
--variant $TMPDIR/$1.tcga.snps.vcf \
-O $TMPDIR/$1.snps.recalibrated.filtered.vcf \
--filter-expression "QD < 2.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0" \
--filter-name "my_snp_filter" \

mv $TMPDIR/$1.snps.recalibrated.filtered.vcf $3/$1_gatk.vcf
echo *****done $1*****
