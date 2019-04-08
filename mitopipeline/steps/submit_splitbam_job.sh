#!/bin/bash
if [ ! -f /mnt/rds/txl80/LaframboiseLab/tyk3/scripts/slurm/split_bam_$1.slurm ];
then
FILE="/mnt/rds/txl80/LaframboiseLab/tyk3/scripts/slurm/split_bam_$1.slurm"
echo '#!/bin/bash' >> $FILE
echo '#SBATCH --mail-user=tyk3@case.edu' >> $FILE
echo '#SBATCH --mail-type=ALL' >> $FILE
echo '#SBATCH -n 1' >> $FILE
echo '#SBATCH -N 1' >> $FILE
echo '#SBATCH --time=24:00:00' >> $FILE
echo '#SBATCH --mem=30GB' >> $FILE
echo '#SBATCH -J split_bam_'"$1" >> $FILE
echo '#SBATCH -A txl80' >> $FILE

echo 'bash /mnt/rds/txl80/LaframboiseLab/tyk3/scripts/split_bams.sh ' "$1 $2" ' >> /mnt/rds/txl80/LaframboiseLab/tyk3/TARGET_PIPELINE/STDOUT/'"$1"'_submit_bam.out 2>&1' >> $FILE
sleep 1
fi

cd /mnt/rds/txl80/LaframboiseLab/tyk3/scripts/slurm

N=4
batchId=`sbatch /mnt/rds/txl80/LaframboiseLab/tyk3/scripts/slurm/split_bam_$1.slurm | awk -v N=$N '{print $N}'`

queue=$(squeue -u tyk3)
inQueue=$(echo "$queue" | grep $batchId)
#if the batchid is still within the squeue output
while [ ! -z "$inQueue" ]
do
sleep 1m
queue=$(squeue -u tyk3)
inQueue=$(echo "$queue" | grep $batchId)
done

rm split_bam_$1.slurm
rm slurm-$batchId.out
