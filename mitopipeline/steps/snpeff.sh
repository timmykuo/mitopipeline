# $1 filename
# $2 path to snpeff tool
# $3 path to store 
VCFS="/mnt/rds/txl80/LaframboiseLab/tyk3/scripts/gatk3_vcfs"
SNPEFF="/mnt/rds/txl80/LaframboiseLab/tyk3/TARGET_PIPELINE/snpeffs"
cd /mnt/rds/txl80/LaframboiseLab/tyk3/Tools/snpEff/snpEff
module load intel/17
module load openmpi/2.0.1

java -Xmx4g -jar snpEff.jar GRCh38.86 $VCFS/$1.snps.recalibrated.filtered.vcf > $SNPEFF/$1_snpEff.vcf

java -Xmx4g -jar snpEff.jar GRCh38.86 -classic $VCFS/$1.snps.recalibrated.filtered.vcf > $SNPEFF/$1_formatEff.vcf
