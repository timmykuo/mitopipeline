#!/bin/bash

#$1 is the file id

TOOLS="/mnt/rds/txl80/LaframboiseLab/tyk3/Tools/annovar"
VCFS="/mnt/rds/txl80/LaframboiseLab/tyk3/scripts/gatk3_vcfs"
ANNOVAR="/mnt/rds/txl80/LaframboiseLab/tyk3/TARGET_PIPELINE/annovar"
PILEUPS="/mnt/rds/txl80/LaframboiseLab/tyk3/TARGET_PIPELINE/pileups"
cd $ANNOVAR

sed 's/MT/M/g' $VCFS/$1.snps.recalibrated.filtered.vcf > $ANNOVAR/$1.vcf

#convert vcf file to avinput file
perl $TOOLS/convert2annovar.pl -format vcf4 $ANNOVAR/$1.vcf  > $1.avinput

#perl $TOOLS/convert2annovar.pl -format vcf4 $VCFS/$1.snps.recalibrated.filtered.vcf  -outfile $1.avinput -includeinfo -comment -allsample 

perl $TOOLS/table_annovar.pl $1.avinput $TOOLS/humandb/ -remove -protocol dbnsfp33a -operation f -build hg38 -nastring .

rm $ANNOVAR/$1.vcf

