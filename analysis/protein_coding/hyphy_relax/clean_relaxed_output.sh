rm relax_pval.all
rm relax_K.all
for GROUP in Ratite VL RND
do
	perl -p -e 's/\// /g' relax_parsed_Pval_${GROUP}_2017-11-20 | perl -p -e 's/:/ /' | perl -p -e 's/\[RELAX\] Likelihood ratio test for relaxation on Test branches, p = //' | awk -v group="$GROUP" 'BEGIN{OFS=" "} {print group, $0}' >> relax_pval.all
	perl -p -e 's/\// /g' relax_parsed_K_${GROUP}_2017-11-20 | perl -p -e 's/:/ /' | perl -p -e 's/\[RELAX\] Log\(L\) = \S+ Relaxation parameter K = //' | awk -v group="$GROUP" 'BEGIN{OFS=" "} {print group, $0}' >> relax_K.all   
done
