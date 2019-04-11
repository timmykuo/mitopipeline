#!/bin/bash

#$1 is the fileNAME
#$2 is the start directory
#$3 is the out directory
#$4 is toosl directory 

TOOLS=$4
VCFS=$2"/gatk3_vcfs"
ANNOVAR=$3"/annovar"
PILEUPS=$3"/pileups"

sed 's/MT/M/g' $VCFS/$1.snps.recalibrated.filtered.vcf > $ANNOVAR/$1.vcf

#convert vcf file to avinput file
perl $TOOLS/convert2annovar.pl -format vcf4 $ANNOVAR/$1.vcf  > $ANNOVAR/$1.avinput

#perl $TOOLS/convert2annovar.pl -format vcf4 $VCFS/$1.snps.recalibrated.filtered.vcf  -outfile $1.avinput -includeinfo -comment -allsample 

perl $TOOLS/table_annovar.pl $ANNOVAR/$1.avinput $TOOLS/humandb/ -remove -protocol dbnsfp33a -operation f -build hg38 -nastring . > $3/$1.avoutput

rm $ANNOVAR/$1.vcf

