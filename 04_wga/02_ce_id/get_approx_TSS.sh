awk 'BEGIN { OFS = "\t" } {if ($6 == "+") {print $1,$2,$2+1,$4,$5,$6} else print $1,$3-1,$3,$4,$5,$6}' $1 | sort -k1,1V -k2,2n 
