#!/usr/bin/bash

cd ~/Human-stool-library
for i in *_R1.fastq.gz
do
barcode=$(echo $i | head -c -13)
barcode=$(echo $barcode | tail -c -9)
echo $barcode
a=~/Human-stool-library/*$barcode*
cat $a > $barcode".fastq.gz"
done

