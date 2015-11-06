#!/bin/sh

for FILE in $(ls *.txt)
do
python translation_tester.py $FILE
done
