for GROUP in ratite vl rnd
do
	perl -p -e 's/\// /g' ${GROUP}_pval | perl -p -e 's/:/ /' | perl -p -e 's/\[RELAX\] Likelihood ratio test for relaxation on Test branches, p = //' | awk -v group="$GROUP" 'BEGIN{OFS=" "} {print group, $0}' >> relax_pval.all
	perl -p -e 's/\// /g' ${GROUP}_K | perl -p -e 's/:/ /' | perl -p -e 's/\[RELAX\] Log\(L\) = \S+ Relaxation parameter K = //' | awk -v group="$GROUP" 'BEGIN{OFS=" "} {print group, $0}' >> relax_K.all   
done
