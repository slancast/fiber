#slancast@scg4.stanford.edu:/srv/gsfs0/projects/snyder/slancast/fiber/rnaseq
#qsub -cwd -l h_vmem=8g -l h_rt=168:00:00 -l s_rt=168:00:00 DESeq2pipe.sh

module load r/3.4.1
cat DESeq2_fiber_week_B.r
Rscript DESeq2_fiber_week_B.r