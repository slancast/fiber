!#/usr/bin/bash
 
cd /home/slancast/RNA/180607_COOPER_0203_AHV27JBBXX_L7_Fiber-iPOP-4/stage0_bcl2fastq/fastqs

unset array
counter=0
declare -a array

for i in *R2.fastq.gz
do
barcode_temp=$(echo $i | head -c -13)
barcode=$(echo $barcode_temp | tail -c -17)
echo $barcode
array[counter]=$barcode
counter=$((counter+1))
done

echo ${array[*]}

cd ../../../

for i in ${array[*]}
do
file=$i"*_R2.fastq.gz"
file="./*Fiber-iPOP-4/stage0_bcl2fastq/fastqs/*"$file
output="./concatenated/"$i
output=$output"_R2_concatenated.fastq.gz"
cat $file > $output 
file=$i"*_R1.fastq.gz"
file="./*Fiber-iPOP-4/stage0_bcl2fastq/fastqs/*"$file
output="./concatenated/"$i
output=$output"_R1_concatenated.fastq.gz"
cat $file > $output
done


