#!/usr/bin
# $1 is filename
# $2 is start directory
# $3 is out directory
# $4 is the tools directory
# $5 is the steps directory
# $6 is the refs directory

VCFS=$2
SNPEFF=$3"/snpEff"
TOOLS=$4
#last string following / delimeter will be name of the previous job
filetype=$(awk -F/ '{print $NF}' <<< "$2" | awk '{print tolower($0)}')

java -Xmx4g -jar $TOOLS/snpEff.jar GRCh38.86 $VCFS/$1_$filetype.vcf > $SNPEFF/$1_snpEff.vcf