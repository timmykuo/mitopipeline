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

#$1 is filename
#$2 is start directory
#$3 is out directory
module load intel/17
module load openmpi/2.0.1
module load samtools
TMPDIR=$3"/temp"
BAMS=$2"/no_numts_bams"
#### $2 is the OUT DIR where you want the final.vcf file to be saved.

cd /mnt/rds/txl80/LaframboiseLab/tyk3/TARGET_PIPELINE/no_numts_bams
echo *****AddOrReplaceReadGroups $1*****
java -Xmx8g -jar ./tools/picard-tools-1.93/AddOrReplaceReadGroups.jar \
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

java -Xmx8g -jar ./tools/picard-tools-1.93/MarkDuplicates.jar \
INPUT=$TMPDIR/$1.tcga.sort.2.bam \
OUTPUT=$TMPDIR/$1.tcga.marked.bam \
METRICS_FILE=$TMPDIR/$1.metricsfile.txt \
CREATE_INDEX=true \

java -Xmx8g -jar ./tools/picard-tools-1.93/FixMateInformation.jar \
INPUT=$TMPDIR/$1.tcga.marked.bam \
OUTPUT=$TMPDIR/$1.tcga.marked.realigned.fixed.bam \
SO=coordinate \
VALIDATION_STRINGENCY=SILENT \
CREATE_INDEX=true
echo *****done*****

echo *****BaseRecalibrator $1*****
./tools/gatk-4.0.11.0/gatk --java-options "-Xmx8g" BaseRecalibrator \
-I $TMPDIR/$1.tcga.marked.realigned.fixed.bam \
-R /mnt/rds/txl80/LaframboiseLab/tyk3/REFs/MitoREFs/rCRS-MT.fa \
-O $TMPDIR/$1.recal_data.grp \
--known-sites /home/tyk3/GATK-ReAnalyze/dbsnp/mtONLY.vcf 
echo *****done*****

echo *****PrintReads $1*****
./tools/gatk-4.0.11.0/gatk --java-options "-Xmx8g" ApplyBQSR \
-R /mnt/rds/txl80/LaframboiseLab/tyk3/REFs/MitoREFs/rCRS-MT.fa \
-I $TMPDIR/$1.tcga.marked.realigned.fixed.bam \
-bqsr-recal-file $TMPDIR/$1.recal_data.grp \
-O $TMPDIR/$1.tcga.marked.realigned.fixed.read.bam
echo *****done*****

echo *****HaplotypeCaller $1*****
./gatk-4.0.11.0/gatk --java-options "-Xmx10g" HaplotypeCaller \
-R /mnt/rds/txl80/LaframboiseLab/tyk3/REFs/MitoREFs/rCRS-MT.fa \
-I $TMPDIR/$1.tcga.marked.realigned.fixed.read.bam \
--max-reads-per-alignment-start 200 \
--sample-ploidy 100 \
--pcr-indel-model HOSTILE \
-O $TMPDIR/$1.tcga.snps.vcf \
--min-pruning 10 \
-A DepthPerAlleleBySample \
--dbsnp /home/tyk3/GATK-ReAnalyze/dbsnp/mtONLY.vcf \
-DF NotDuplicateReadFilter \
--minimum-mapping-quality 0 

#echo *****Starting DepthOfCoverage*****
#java -Xmx64g -jar /homeG/LaFramboiseLab/sjm67/GenomeAnalysisTK-2.4-9-g532efad/GenomeAnalysisTK.jar  \
#-R /homeG/LaFramboiseLab/sjm67/TCGA/hg19_karyo_rCRS/hg19_karyo_rCRS-MT.fa \
#-T Coverage \
#-o $1.out.table \
#-I $1.tcga.marked.realigned.fixed.read.bam \
#--outputFormat table
#echo *****done*****

echo *****VariantFiltration $1******
./gatk-4.0.11.0/gatk --java-options "-Xmx8g" VariantFiltration \
-R ./REFs/MitoREFs/rCRS-MT.fa \
-V $TMPDIR/$1.tcga.snps.vcf \
-O $3/$1.snps.recalibrated.filtered.vcf \
--filter-expression "QD < 2.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0" \
--filter-name "my_snp_filter" 

#echo *****done*****

echo *****Creating Report $1*****
#perl /home/sxg501/GATK-ReAnalyze/frequencies.pl $1.snps.recalibrated.filtered.vcf $1.out.txt
echo *****done*****
#cp $STOR/$3.snps.recalibrated.filtered.vcf $2/
echo *****Creating Report $1*****
#perl /home/sxg501/GATK-ReAnalyze/frequencies.pl $1.snps.recalibrated.filtered.vcf $1.out.txt
echo *****done*****
#cp $STOR/$3.snps.recalibrated.filtered.vcf $2/
#echo ***making pileup file****
#samtools mpileup -f ~/REFs/MitoREFs/rCRS-MT.fa $STOR/$3.tcga.marked.realigned.fixed.read.bam > ~/GATK-ReAnalyze/Pileups/$3.pileup
#perl ~/GATK-ReAnalyze/Pileups/pileupstats2.pl ~/GATK-ReAnalyze/Pileups/$3.pileup ~/GATK-ReAnalyze/Pileups/$3.pOUT
#bash ~/GATK-ReAnalyze/Pileups/pOUT-mod.sh ~/GATK-ReAnalyze/Pileups/$3.pOUT

echo *****done $1*****
