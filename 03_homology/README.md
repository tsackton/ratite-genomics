Homology Inference
===============

Homology inference and protein alignment in 4 steps:

1. Extract longest protein translation for each gene from GFF files (00_gff_parsing)
2. Run OMA to infer homologous groups among them (01_oma)
3. Extract proteins not assigned to a group, and use HMMER to screen for possible merges, producing new homologous groups (02_hmm)
4. Run alignment and filtering on those groups to produce final analysis set (03_protein_coding_alignment)

Files in the main directory are largely useful output files and code.

alignment_summary_stats: summary statistics from initial alignments used for error checking homology assignments and filtering
all_hog_info.tsv: primary data table of homologous groups for filtering and processing; columns with species codes (from oma_species_list.csv) represent number of proteins assigned to that hog for that species (NA = 0)
good_PAML_HOG_protein_transcript_info: transcript info for HOGs that passed filtering for PAML and related analysis
HOG_final_alignment_seqids: seqids for further analysis
new_hog_list.txt: hog assignments in long format
oma_species_list.csv: list of species used for OMA analysis with data sources
sp_list.txt: list of all species abbreviations that should be present in output files
updated_hog_matrix.txt: hog assignments as a matrix (wide format)

Code
----
hog_qc.R
check_homology_matrix.R
make_protein_info_table.R

Note
------
updated_hog_matrix.txt is the homologous group assignments after the HMM step. The original is the in the OMA directory.

