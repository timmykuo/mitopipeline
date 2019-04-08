#!/bin/bash
if [ ! -f /mnt/rds/txl80/LaframboiseLab/tyk3/scripts/slurm/gatk_test_$1.slurm ];
then
FILE="/mnt/rds/txl80/LaframboiseLab/tyk3/scripts/slurm/gatk_test_$1.slurm"
echo '#!/bin/bash' >> $FILE
echo '#SBATCH --mail-user=tyk3@case.edu' >> $FILE
echo '#SBATCH --mail-type=ALL' >> $FILE
echo '#SBATCH -n 1' >> $FILE
echo '#SBATCH -N 1' >> $FILE
echo '#SBATCH --time=24:00:00' >> $FILE
echo '#SBATCH --mem=25GB' >> $FILE
echo '#SBATCH -J gatk_test_'"$1" >> $FILE
echo '#SBATCH -A txl80' >> $FILE
echo 'bash /mnt/rds/txl80/LaframboiseLab/tyk3/scripts/gatktest.sh '"$1"' >> /mnt/rds/txl80/LaframboiseLab/tyk3/TARGET_PIPELINE/STDOUT/gatk_test_'"$1"'.out 2>&1' >> $FILE
sleep 1
fi

cd /mnt/rds/txl80/LaframboiseLab/tyk3/scripts/slurm/

sleep 1
#run slurm job
N=4
batchId=`sbatch /mnt/rds/txl80/LaframboiseLab/tyk3/scripts/slurm/gatk_test_$1.slurm | awk -v N=$N '{print $N}'`

queue=$(squeue -u tyk3)
inQueue=$(echo "$queue" | grep $batchId)
#if the batchid is still within the squeue output
while [ ! -z "$inQueue" ]
do
sleep 1m
queue=$(squeue -u tyk3)
inQueue=$(echo "$queue" | grep $batchId)
done

if grep -Fq "InternalError" /mnt/rds/txl80/LaframboiseLab/tyk3/TARGET_PIPELINE/STDOUT/gatk_test_${1}.out;
then
echo $1
else
rm gatk_test_$1.slurm
rm slurm-$batchId.out
fi
