#!/usr/bin/bash
#RSEM for Mike
#16 standard computers

mkdir trimmomatic

#The following will run trimmomatic
#Mike's second sequencing required trimming
for i in *R1.fastq.gz
do 
barcode=$(echo $i | head -c -13)
read2=$barcode"_R2.fastq.gz"
output_forward_unpaired="./trimmomatic/"$barcode"output_forward_unpaired.fq.gz"
output_reverse_unpaired="./trimmomatic/"$barcode"output_reverse_unpaired.fq.gz"
output_forward_paired="./trimmomatic/"$barcode"output_forward_paired.fq.gz"
output_reverse_paired="./trimmomatic/"$barcode"output_reverse_paired.fq.gz"
trimmomatic PE $i $read2 $output_forward_paired $output_forward_unpaired $output_reverse_paired $output_reverse_unpaired CROP:75
done

#This snippet is running the code for the first rsem
gunzip *.gz #must be gunzipped
for i in *R1.fastq
do 
barcode=$(echo $i | head -c -10)
read2=$barcode"_R2.fastq"
~/RSEM-1.3.1/rsem-calculate-expression -p 15 --paired-end --star --star-path /home/slancast/STAR-master/bin/Linux_x86_64/ $i $read2 ~/reference_genomes/human_ensembl $barcode
done
gzip *.fastq

#Running the second resem on the mouse genomes
for i in ./trimmomatic/*output_forward_paired.fq
do 
barcode=$(echo $i | head -c -25)
read2=$barcode"output_reverse_paired.fq"
~/RSEM-1.3.1/rsem-calculate-expression -p 15 --paired-end --star --star-path /home/slancast/STAR-master/bin/Linux_x86_64/ $i $read2 ~/mouse_reference/mouse_references $barcode
done



gzip ./trimmomatic/*.fq
sudo poweroff


