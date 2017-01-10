perl -p -e 's/\// /g' relax_parsed | perl -p -e 's/\[RELAX\] Likelihood ratio test for relaxation on Test branches, p = //' > relax_parsed.clean
perl -p -e 's/\// /g' relax_parsed_K | perl -p -e 's/\[RELAX\] Log\(L\) = \S+ Relaxation parameter K = //' > relax_parsed_K.clean
