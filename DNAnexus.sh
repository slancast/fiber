#!usr/bin/bash

unset array
counter=0
declare -a array

for i in Fiber-iPOP-4_concatenated/*R1_concatenated.fastq.gz
do
echo $i
#cutting barcode info from file name
barcode_temp=$(echo $i | head -c -26)
echo $barcode_temp
barcode=$(echo $barcode_temp | tail -c -17)
echo $barcode
#building array and counting
array[counter]=$barcode
counter=$((counter+1))
echo $counter
done

echo ${array[*]}
echo ${#array[@]}

source dx-toolkit/environment

for i in ${array[*]}
do
dx describe "Fiber-iPOP-4_concatenated/"${i}"_R1_concatenated.fastq.gz"
yes | dx run workflows/EncodePipelineiPOP4  -i0.reads1="Fiber-iPOP-4_concatenated/"${i}"_R1_concatenated.fastq.gz" -i0.reads2="Fiber-iPOP-4_concatenated/"${i}"_R2_concatenated.fastq.gz"
done

#Saving code potentially for later
#printf "Fiber-iPOP_concatenated/"${i}"_R1_concatenated.fastq.gz\nFiber-iPOP_concatenated/"${i}"_R2_concatenated.fastq.gz\n\nY\n" | dx run workflows/EncodePipelineiPOP3 


unset array
counter=0
declare -a array

for i in Fiber-iPOP-3_concatenated/*R1_concatenated.fastq.gz
do
echo $i
#cutting barcode info from file name
barcode_temp=$(echo $i | head -c -26)
echo $barcode_temp
barcode=$(echo $barcode_temp | tail -c -17)
echo $barcode
#building array and counting
array[counter]=$barcode
counter=$((counter+1))
echo $counter
done

echo ${array[*]}
echo ${#array[@]}

source dx-toolkit/environment

for i in ${array[*]}
do
dx describe "Fiber-iPOP-3_concatenated/"${i}"_R1_concatenated.fastq.gz"
yes | dx run workflows/EncodePipelineiPOP3  -i0.reads1="Fiber-iPOP-3_concatenated/"${i}"_R1_concatenated.fastq.gz" -i0.reads2="Fiber-iPOP-3_concatenated/"${i}"_R2_concatenated.fastq.gz"
done

