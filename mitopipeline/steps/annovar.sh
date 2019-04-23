#!/bin/bash
#$1 is the fileNAME
#$2 is the start directory
#$3 is the out directory
#$4 is tools directory 
#$5 is the steps directory
#$6 is the refs directory

TOOLS=$4"/annovar"
VCFS=$2
ANNOVAR=$3

#last string following / delimeter will be name of the previous job
#only set filetype if the value was a step in the pipeline
filetype="_"$(awk -F/ '{print $NF}' <<< "$2" | awk '{print tolower($0)}')
if [ "$filetype" != "gatk" ];
then
filename=""
fi
sed 's/MT/M/g' $VCFS/$1$filetype.vcf > $ANNOVAR/$1.vcf

#convert vcf file to avinput file
perl $TOOLS/convert2annovar.pl -format vcf4 $ANNOVAR/$1.vcf  > $ANNOVAR/$1.avinput

perl $TOOLS/table_annovar.pl $ANNOVAR/$1.avinput $TOOLS/humandb/ -remove -protocol dbnsfp33a -operation f -build hg38 -nastring . > $3/$1_annovar.avoutput
rm $ANNOVAR/$1.vcf
rm $ANNOVAR/$1.avinput
rm $ANNOVAR/$1_annovar.avoutput
mv $ANNOVAR/$1.avinput.hg38_multianno.txt $ANNOVAR/$1_avinput.hg38_multianno.txt

