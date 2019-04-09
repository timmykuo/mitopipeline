# $1 filename
# $2 path to snpeff tool
# $3 path to store

# $1 is filename
# $2 is start directory
# $3 is out directory
VCFS=$2
#vcfs was "/mnt/rds/txl80/LaframboiseLab/tyk3/scripts/gatk3_vcfs"
SNPEFF= $3
#snpeff was "/mnt/rds/txl80/LaframboiseLab/tyk3/TARGET_PIPELINE/snpeffs"
cd ./tools/snpEff
module load intel/17
module load openmpi/2.0.1

java -Xmx4g -jar snpEff.jar GRCh38.86 $VCFS/$1.snps.recalibrated.filtered.vcf > $SNPEFF/$1_snpEff.vcf

java -Xmx4g -jar snpEff.jar GRCh38.86 -classic $VCFS/$1.snps.recalibrated.filtered.vcf > $SNPEFF/$1_formatEff.vcf
