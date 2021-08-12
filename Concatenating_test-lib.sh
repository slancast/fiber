#!/usr/bin/bash

cd ~/test-lib-combined
for i in ~/test-lib2/*_R1.fastq.gz
do
barcode=$(echo $i | head -c -13)
barcode=$(echo $barcode | tail -c -9)
echo $barcode
a=~/test-lib2/*$barcode*
b=~/test-lib1/*$barcode*
cat $a $b > $barcode".fastq.gz"
done

