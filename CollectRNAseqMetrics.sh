#!/bin/bash
#SBATCH --job-name=DE
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --time=6-00:00:00
#SBATCH --mem=20G
#SBATCH --mail-type=ALL ##BEGIN,FAIL,END
#SBATCH --mail-user=slancast@stanford.edu
#SBATCH --chdir=/labs/mpsnyder/fiber-ipop/iPOP1
#SBATCH --account=mpsnyder

#This is going to rerun the findings I made before on scg.

######## to retrieve info before run:
date
echo "I ran on host: $(hostname -s)"

echo "SLURM Environment is:"
env | grep "SLURM" | sort

echo "My limits are:"
ulimit -a

cat CollectRNAseqMetrics.sh

for i in *genome.bam
do
echo $i
output=$(echo $i | head -c -4)
output=$output"RNA_Metrics"
echo $output
java -Xmx8g -jar /srv/gsfs0/software/picard-tools/2.8.0/picard.jar CollectRnaSeqMetrics I=$i OUTPUT=$output REF_FLAT=/gssc/Users/weie/refs/GRCh38_GENCODE/annotation/refFlat.txt.gz STRAND=NONE
done

echo "This job finished on: "
date


