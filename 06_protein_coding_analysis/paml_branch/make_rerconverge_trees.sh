grep "tree1" $1.out | grep -v "_" |
    grep -v "^tree1[[:space:]]*$" | awk 'BEGIN {OFS = "\t"; FS="\t"} {if (NF == 3) print "HOG"$2, $3}' |
    grep -v "^HOG[[:digit:]]*[[:space:]]*$" > ${1}_RERconverge.txt
