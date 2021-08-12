#!/usr/bin/bash
 
#finding a folder with the fastqs
cd /home/slancast/RNA/180607_COOPER_0203_AHV27JBBXX_L1_Fiber-iPOP-3/stage0_bcl2fastq/fastqs

#making sure no arrays are already made, making array, reseting counter
unset array
counter=0
declare -a array

#looping over files to pull out barcode info
for i in *R2.fastq.gz
do
#cutting barcode info from file name
barcode_temp=$(echo $i | head -c -13)
barcode=$(echo $barcode_temp | tail -c -17)
echo $barcode
#building array and counting
array[counter]=$barcode
counter=$((counter+1))
done

echo ${array[*]}

cd ../../../

#finding all the R2s and R1s of a barcode and concatenating them together.
for i in ${array[*]}
do
#finding all the R2s with that barcode
file=$i"*_R2.fastq.gz"
file="./*Fiber-iPOP-3/stage0_bcl2fastq/fastqs/*"$file
#creating output location
output="./Fiber-iPOP-3_concatenated/"$i
output=$output"_R2_concatenated.fastq.gz"
cat $file > $output 
file=$i"*_R1.fastq.gz"
file="./*Fiber-iPOP-3/stage0_bcl2fastq/fastqs/*"$file
output="./Fiber-iPOP-3_concatenated/"$i
output=$output"_R1_concatenated.fastq.gz"
cat $file > $output
done


