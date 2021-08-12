#!/bin/bash
#SBATCH --job-name=batch_effects
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=6-00:00:00
#SBATCH --mem=100G
#SBATCH --mail-type=ALL ##BEGIN,FAIL,END
#SBATCH --mail-user=slancast@stanford.edu
#SBATCH --workdir=/srv/gsfs0/projects/snyder/fiber-ipop/batch_effects/
#SBATCH --account=mpsnyder

######## to retrieve info before run:
date
echo "I ran on host: $(hostname -s)"

echo "SLURM Environment is:"
env | grep "SLURM" | sort

echo "My limits are:"
ulimit -a

Rscript /srv/gsfs0/projects/snyder/fiber-ipop/batch_effects/Batch_Effects_gf_scg.R

echo "This job finished on: "
date
