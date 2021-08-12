#!/bin/bash
#SBATCH --job-name=DE
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=6-00:00:00
#SBATCH --mem=1G
#SBATCH --mail-type=ALL ##BEGIN,FAIL,END
#SBATCH --mail-user=slancast@stanford.edu
#SBATCH --chdir=/labs/mpsnyder/fiber-ipop/
#SBATCH --account=mpsnyder

#This is going to rerun the findings I made before on scg.

######## to retrieve info before run:
date
echo "I ran on host: $(hostname -s)"

echo "SLURM Environment is:"
env | grep "SLURM" | sort

echo "My limits are:"
ulimit -a

tar -czvf second-metagenome-analyses.tar.gz /labs/mpsnyder/fiber-ipop/second-metagenome-analyses

echo "This job finished on: "
date


