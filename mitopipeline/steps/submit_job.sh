#!/bin/bash
#$1 is filename
#$2 is script name
#$3 is slurm directory
#$4 is start directory
#$5 is out directory
#$6 is the tools directory 
# $7 is steps directory
# $8 is the refs directory 
# $9 is email address
if [ ! -f $3/$2_$1.slurm ];
then
FILE="$3/$2_$1.slurm"
touch $3/$2_$1.slurm
echo '#!/bin/bash' >> $FILE
echo '#SBATCH -n 1' >> $FILE
echo '#SBATCH -N 1' >> $FILE
echo '#SBATCH --time=12:00:00' >> $FILE
echo '#SBATCH -J '"$2_$1" >> $FILE
echo '#SBATCH --mail-user='"$9" >> $FILE
echo '#SBATCH --mail-type=ALL' >> $FILE
echo '#SBATCH --mem='"50"'GB' >> $FILE
echo 'bash '"$7"'/'"$2" "$1 $4 $5 $6 $7 $8"' >> '"$3"'/STDOUT'"$2"'_'"$1"'.out 2>&1' >> $FILE
sleep 1
fi

sleep 1
#run slurm job
N=4
batchId=`sbatch $3/$2_$1.slurm | awk -v N=$N '{print $N}'`

queue=$(squeue -u tyk3)
inQueue=$(echo "$queue" | grep $batchId)
#if the batchid is still within the squeue output
while [ ! -z "$inQueue" ]
do
sleep 1m
queue=$(squeue -u tyk3)
inQueue=$(echo "$queue" | grep $batchId)
done

rm slurm-$batchId