#!/bin/bash
# $1 is filename
# $2 is start directory
# $3 is out directory
# $4 is the tools directory
# $5 is the steps directory
# $6 is the refs directory
VCFS=$2
SNPEFF=$3
TOOLS=$4
#last string following / delimeter will be name of the previous job
#only set filetype if the value was a step in the pipeline
filetype="_"$(awk -F/ '{print $NF}' <<< "$2" | awk '{print tolower($0)}')
if [ "$filetype" != "gatk" ];
then
filename=""
fi

java -Xmx4g -jar $TOOLS/snpEff/snpEff.jar -c $TOOLS/snpEff/snpEff.config -v GRCh38.86 $VCFS/$1$filetype.vcf > $SNPEFF/$1_snpeff.vcf