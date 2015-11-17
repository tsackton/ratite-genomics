#!/bin/sh

for FILE in $(ls *.txt)
do
python translation_tester_best.py $FILE
done
