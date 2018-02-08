#ratite removed CNEE analysis

#get ratite-only CNEEs
bedtools intersect -a final_ces_noratite.tree2.bed -b final_ces.tree2.bed -v -wa > final_ces_rronly.tree2.bed

#filter to remove those shorter than 50 bp
awk '($3-$2) >= 50' final_ces_rronly.tree2.bed > final_ces_rronly_filtered.tree2.bed

