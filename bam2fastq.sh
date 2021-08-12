#!/bin/bash
#SBATCH --job-name=bam2fastq
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=6-00:00:00
#SBATCH --mem=10G
#SBATCH --mail-type=ALL ##BEGIN,FAIL,END
#SBATCH --mail-user=slancast@stanford.edu
#SBATCH --workdir=/srv/gsfs0/projects/snyder/fiber-ipop/iPOP1
#SBATCH --account=mpsnyder

module load bedtools

######## to retrieve info before run:
date
echo "I ran on host: $(hostname -s)"

echo "SLURM Environment is:"
env | grep "SLURM" | sort

echo "My limits are:"
ulimit -a

for i in *genome.bam
do
sequence_run=$(echo $i | head -c -17)
output="/srv/gsfs0/projects/snyder/fiber-ipop/iPOP1/fastqs/"$sequence_run".fastq"
echo $output
bamToFastq -i $i -fq $output
done

echo "This job finished on: "
date
