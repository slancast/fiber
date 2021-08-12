#!/bin/bash
#SBATCH --job-name=DE
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=5
#SBATCH --time=6-00:00:00
#SBATCH --mem=20G
#SBATCH --mail-type=ALL ##BEGIN,FAIL,END
#SBATCH --mail-user=slancast@stanford.edu
#SBATCH --workdir=/srv/gsfs0/projects/snyder/slancast/fiber/rnaseq/fall2019
#SBATCH --account=mpsnyder


######## to retrieve info before run:
date
echo "I ran on host: $(hostname -s)"

echo "SLURM Environment is:"
env | grep "SLURM" | sort

echo "My limits are:"
ulimit -a

Rscript ./deseq_baseline_scg.R

echo "This job finished on: "
date
