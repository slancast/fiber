#!/usr/bin/bash
#Firstly I will generate the reference files
#Making the referenece in star is more informative than doing it in RSEM
/home/slancast/STAR-master/bin/Linux_x86_64/STAR --runThreadN 15 --runMode genomeGenerate --genomeDir ~/reference_genomes/ --genomeFastaFiles ~/reference_genomes/GCF_000001405.31_GRCh38.p5_genomic.fna --sjdbGTFfile ~/reference_genomes/human_refseq.gtf

#Now time to run star, it looks like those commands are pretty basic
/home/slancast/STAR-master/bin/Linux_x86_64/STAR --runThreadN 15 --genomeDir ~/reference_genomes/ --readFilesIn /home/slancast/Mike/SCGPM_062019-HEPG2-SCFA-RNA-seq-pool-MN_H2HWW_L4_ATCACG_R1.fastq ~/Mike/SCGPM_062019-HEPG2-SCFA-RNA-seq-pool-MN_H2HWW_L4_ATCACG_R2.fastq
#OK STAR Works on to RSEM

#Preparing RSEM Reference
#Need to run both this and STAR
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/vertebrate_mammalian/Homo_sapiens/all_assembly_versions/GCF_000001405.31_GRCh38.p5/GCF_000001405.31_GRCh38.p5_genomic.gff.gz ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/vertebrate_mammalian/Homo_sapiens/all_assembly_versions/GCF_000001405.31_GRCh38.p5/GCF_000001405.31_GRCh38.p5_genomic.fna.gz
~/RSEM-1.3.1/rsem-prepare-reference -p 15 --gtf Homo_sapiens.GRCh38.83.gtf --star --star-path /home/slancast/STAR-master/bin/Linux_x86_64 Homo_sapiens.GRCh38.dna.primary_assembly.fa ~/reference_genomes/human_ensembl

#Running RSEM
~/RSEM-1.3.1/rsem-calculate-expression -p 15 --paired-end --star --star-path /home/slancast/STAR-master/bin/Linux_x86_64/ ~/Mike/SCGPM_062019-HEPG2-SCFA-RNA-seq-pool-MN_H2HWW_L4_ATCACG_R1.fastq ~/Mike/SCGPM_062019-HEPG2-SCFA-RNA-seq-pool-MN_H2HWW_L4_ATCACG_R2.fastq ~/reference_genomes/human_ensembl ~/RSEM_out/ATCACG

#I began writing this trimmomatic code, but I'll move it to "RSEM4Mike.sh"
trimmomatic PE SCGPM_MN-SCFA-CT26-RNA-seq-pool-10232019_H7WYN_L8_ACAGTG_R1.fastq.gz SCGPM_MN-SCFA-CT26-RNA-seq-pool-10232019_H7WYN_L8_ACAGTG_R2.fastq.gz output_forward_paired.fq.gz output_forward_unpaired.fq.gz output_reverse_paired.fq.gz output_reverse_unpaired.fq.gz CROP:25

#I got the mouse reference from here: https://www.gencodegenes.org/mouse/
#preparing mouse reference. RSEM prepare reference didn't generate this file: genomeParameters.txt
/home/slancast/STAR-master/bin/Linux_x86_64/STAR --runThreadN 15 --runMode genomeGenerate --genomeDir ~/mouse_reference/ --genomeFastaFiles ~/mouse_reference/GRCm38.primary_assembly.genome.fa --sjdbGTFfile ~/mouse_reference/gencode.vM23.annotation.gtf
~/RSEM-1.3.1/rsem-prepare-reference -p 15 --gtf gencode.vM23.annotation.gtf --star --star-path /home/slancast/STAR-master/bin/Linux_x86_64 GRCm38.primary_assembly.genome.fa ~/mouse_reference/mouse_references
