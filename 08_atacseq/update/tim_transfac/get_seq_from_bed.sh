#!/bin/bash
module load bedtools
SP=$1
BED=$2
cp $BED temp.bed
perl -p -i -e "s/SPECIES/$SP/" temp.bed 
bedtools getfasta -fi $SP.fa -bed temp.bed -fo - |
   perl -p -i -e 's/^>.*$//' |
   perl -p -i -e 's/^$/NNNNNNNNNN/' |
   tr -d '\n' > ${SP}_temp.fa
echo ">$SP" | cat - ${SP}_temp.fa > ${SP}_${BED}.fa
rm ${SP}_temp.fa
rm temp.bed
