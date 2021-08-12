
#gunzip *.gz
for i in *R1.fastq.gz
do 
barcode=$(echo $i | head -c -10)
read2=$barcode"_R2.fastq"
~/RSEM-1.3.1/rsem-calculate-expression -p 15 --calc-ci --calc-pme --paired-end --star --star-path /home/slancast/STAR-master/bin/Linux_x86_64 $i $read2 ~/reference_genomes/human_ensembl $barcode
done




