Analysis of conserved non-coding elements
===============

This directory contains a large number of mostly R scripts used in the analysis of non-coding elements. These scripts do the following:

Processing CNEE output from PhyloAcc and phyloP
----
analyze_CNEEs_combined.R: produces parsed output from raw phyloP and PhyloAcc data by calling the two scripts below and calculating columns

sourced by analyze_CNEEs_combined to do each part of the analysis
analyze_phyloAcc.R
analyze_phyloP.R


Running permutations
----
The following scripts run permutations to test for convergence, enrichment associated with particular genes, or enrichment associated with GO categories:

run_convergence_perms_sbatch.sh
run_gene_perms_sbatch.sh
run_go_perms_sbatch.sh
run_convergence_perms.R
run_enrichment_perms.R
run_gene_perms.R

After permutations are done, these scripts compute ecdfs (for data compactness):
run_go_ecdf_calculation.R
run_gene_ecdf_calculation.R

The following scripts analyze permutation results using either raw permutation output (convergence perms) or ecdfs:
analyze_convergence_permutations.R
analyze_gene_perms.R
analyze_go_permutations.R


Spatial enrichment 
----
The following scripts were used to analyze patterns of spatial enrichment of accelerated CNEEs

analyze_spatial_enrichment.R
plot_spatial_enrichment.R

Misc
---

The following are miscellaneous plotting scripts for supplemental figures:

analyze_ILS.R
Fig13.R
Figure2.R
FigureS15.R
GC_content_anno.R
make_transfac_beds.R
plot_phyloacc_res.R
REVIGO.r


