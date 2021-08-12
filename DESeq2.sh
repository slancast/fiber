##!/usr/bin/bash

#This will run deseq2 on gcloud and stop the instance once it's done to save some ducats

Rscript ~/rnaseq-storage/rsem-results/DESeq2_gcloud_baseline.R  > out.txt

sudo poweroff

R

sudo poweroff
