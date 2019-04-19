#!/bin/bash
#$1 is filename
#$2 is script name
#$3 is OUT/slurm directory
#$4 is start directory
#$5 is out directory
#$6 is the tools directory 
# $7 is steps directory
# $8 is the refs directory 
if [ ! -f $3/slurm/$2_$1.slurm ];
then
FILE="$3/slurm/$2_$1.slurm"
echo '#!/bin/bash' >> $FILE
#echo '#SBATCH --mail-user=tyk3@case.edu' >> $FILE
#echo '#SBATCH --mail-type=ALL' >> $FILE
echo '#SBATCH -n 1' >> $FILE
echo '#SBATCH -N 1' >> $FILE
echo '#SBATCH --time=12:00:00' >> $FILE
echo '#SBATCH --mem=20GB' >> $FILE
echo '#SBATCH -J '"$2_$1" >> $FILE
#echo '#SBATCH -A txl80' >> $FILE
echo 'bash '"$7"'/steps/'"$2"'.sh '"$1 $4 $5 $6 $7 $8"' >> '"$3"'/STDOUT/'"$2"'_'"$1"'.out 2>&1' >> $FILE
sleep 1
fi

cd $3/slurm/

sleep 1
#run slurm job
N=4
batchId=`sbatch $3/slurm/$2_$1.slurm | awk -v N=$N '{print $N}'`

queue=$(squeue -u tyk3)
inQueue=$(echo "$queue" | grep $batchId)
#if the batchid is still within the squeue output
while [ ! -z "$inQueue" ]
do
sleep 1m
queue=$(squeue -u tyk3)
inQueue=$(echo "$queue" | grep $batchId)
done

rm $2_$1.slurm
rm slurm-$batchId.out
