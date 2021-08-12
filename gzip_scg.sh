#!/bin/bash
#SBATCH --job-name=DE
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=6-00:00:00
#SBATCH --mem=2G
#SBATCH --mail-type=ALL ##BEGIN,FAIL,END
#SBATCH --mail-user=slancast@stanford.edu
#SBATCH --chdir=/labs/mpsnyder/fiber-ipop/iPOP1/fastqs
#SBATCH --account=mpsnyder

#it takes long enough to gzip files that it needs to be submitted through sbatch

######## to retrieve info before run:
date
echo "I ran on host: $(hostname -s)"

echo "SLURM Environment is:"
env | grep "SLURM" | sort

echo "My limits are:"
ulimit -a

gzip *.fastq

echo "This job finished on: "
date


