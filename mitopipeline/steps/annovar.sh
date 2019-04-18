#!/bin/bash

#$1 is the fileNAME
#$2 is the start directory
#$3 is the out directory
#$4 is toosl directory 
#$5 is the steps directory
#$6 is the refs directory

TOOLS=$4
VCFS=$2
ANNOVAR=$3"/annovar"
PILEUPS=$3"/pileups"
#last string following / delimeter will be name of the previous job
filetype=$(awk -F/ '{print $NF}' <<< "$2" | awk '{print tolower($0)}')

#sed 's/MT/M/g' $VCFS/$1_$filetype.vcf > $ANNOVAR/$1.vcf

#convert vcf file to avinput file
perl $TOOLS/annovar/convert2annovar.pl -format vcf4 $VCFS/$1_$filetype.vcf  > $ANNOVAR/$1.avinput

perl $TOOLS/annovar/table_annovar.pl $ANNOVAR/$1.avinput $TOOLS/humandb/ -remove -protocol dbnsfp33a -operation f -build hg38 -nastring . > $3/$1.avoutput

#rm $ANNOVAR/$1.vcf

