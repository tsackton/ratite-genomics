Homology Inference
===============

Homology inference and protein alignment in 4 steps:

1. Extract longest protein translation for each gene from GFF files (00_gff_parsing)
2. Run OMA to infer homologous groups among them (01_oma)
3. Extract proteins not assigned to a group, and use HMMER to screen for possible merges, producing new homologous groups (02_hmm)
4. Run alignment and filtering on those groups to produce final analysis set (03_protein_coding_alignment)

new_hog_list.txt and updated_hog_matrix.txt are the homologous group assignments and count matrix after the HMM step
The originals of these files are in the OMA directory

