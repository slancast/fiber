#!/usr/bin/bash

cd ~/library-2-combined
for i in ~/Human-stool-library-2_L3/*_R1.fastq.gz
do
barcode=$(echo $i | head -c -13)
barcode=$(echo $barcode | tail -c -9)
echo $barcode
a=~/Human-stool-library-2_L3/*$barcode*
b=~/Human-stool-library-2_L4/*$barcode*
cat $a $b > $barcode".fastq.gz"
done

