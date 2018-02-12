README: OMA HOMOLOGY - PROTEIN-CODING ALIGNMENTS
===============

PURPOSE: compile datasets to run PAML inference of positive selection on HOGs [hierarchical orthologous groups]

Note that scripts were written for the Odyssey Cluster at Harvard University, and therefore might include code
to account for program stops/restarts on Odyssey's serial_requeue partition, and/or hard-coded paths
to e.g. software installations and user directories.

---------------

Step 1: Alignment and filtering of HOG protein sequences to determine which sequences to retain for CDS alignments

File 'new_hog_list.txt' obtained from Tim Sackton lists protein & gene accessions for all sequences assigned to HOGs

Script: write_hog_fastas.pl used to output protein sequences in fasta format for HOGs that have at least 4 species present.
Protein sequences were retrieved from file 'all_proteins.fa' obtained from Tim Sackton.
Total HOGs: 45,367
< 4 species (omitted): 29,216
Output HOGs: 16,151

Script: mafft_align_template.pl used to align protein sequences with mafft v. 7.245 and default options for high sensitivity alignment.
Output to: HOG_fastas_aligned/ (with 4 subdirectory 'batches')
NOTE: the 2 largest HOGs (HOG2_1 from batch1 & HOG2_2 from batch2) failed to finish alignment within a reasonable timeframe and were omitted

Following alignment, 3 rounds of filtering were used to retain only sequences that align reasonably well within each HOG
Round1 script: oma_filter1.pl
-filters alignments by removing entire columns.  Columns are retained only if one or more of the following applies:
(1) minimum 30% of sequences in alignment have a (non-gap) base
(2) minimum 10 sequences total have a (non-gap) base
(3) there is at least 1 sequence with a (non-gap) base at that position from at least 2 of the 3 major taxonomic groups
    [non-avian outgroups, palaeognaths, neognaths]

Round2 filtering: used the Jarvis et al. (Science 346: 1320-1331) Avian Phylogenomics Project scripts for amino acid alignments, which
masks over poorly aligning regions of individual sequences (rather than omitting entire alignment columns)
Scripts/files:
spotProblematicSeqsBase-W12S4.py
spotProblematicSeqsModules.py
spotPromlematicSeqsModules.pyc
blosum62.txt
(all accessed from: ftp://climb.genomics.cn/pub/10.5524.101001_102000/10/041/ on Sept. 30, 2015)

Filtering scripts were implemented with perl wrapper: oma_filter2.pl, and output from Jarvis et al. script
is provided in directory: HOG_jarvis_filter_output/
Script: apply_oma_filter2_template.pl was then used to perform the actual masking of sequences using output text file from Jarvis et al. scripts

Round3: a redo of round1 filtering to remove gappy/low-representation columns from the alignment matrix following round2 masking of poorly aligning regions

Final filtered output is provided in directory: HOG_fastas_aligned_filtered/

---------------

Step 2: Calculate alignment statistics & use these to choose 'good' sequences and HOGs to retain for PAML analyses

Scripts: write_alignment_summary_stats_per_locus.pl & write_alignment_stats_per_species.pl
Output is provided for both the intial alignments and the post-filtering alignments in: HOG_alignment_stats/

Script: write_PAML_groups.pl & then add_PAML_groups_transcript_ids.pl [using 'all_info.summary' input file provided by Tim Sackton]
Outputs 'good_PAML_HOG_protein_transcript_info' text file with list of HOGs to retain, and what sequences within each HOG are to be used

To retain a HOG, we:
(1) omitted non-avian outgroup sequences
(2) allowed a maximum of 3 'good' sequences per species, and the total number of sequences in the alignment not exceed 1.5X the total number of species
(3) required the presence of at least 50% of all avian species in the alignment, AND of 50% of all palaeognath species, in each case with 'good' sequences
'Good' sequences:
- if the original (unfiltered) sequence length is at least 50% of the average input length for all unfiltered sequences for that locus
- if the filtered sequence length is at least 50% of the original input length for that individual sequence
- AND if there is < 1 gap per bp aligned sequence in the filtered sequence
Failure to meet these criteria resulted in the entire sequence being removed from the HOG

These criteria resulted in the retention of 11,274 HOGs for PAML analyses

---------------

Step 3: Compile CDS fastas for HOGs & perform codon-aware sequence alignment

#####FOR PUBLICLY AVAILABLE GENOMES:
Retrieve NCBI rna GenBank format annotations corresponding to versions used by Tim Sackton in running OMA to generate HOGs
Script: get_oma_rna_genbank.sh

Script: parse_cds_from_rna_genbank.pl
Parse GenBank rna annotation file for each species & output all CDS sequences, with ends padded so all sequences begin in phase0 and end as a multiple of 3.

#####FOR NEW PALAEOGNATH GENOMES:
Use Cufflinks v. 2.2.1 gffread utility, together with new palaeognath genomes from Sackton et al. (2018) & maker annotation GFFs to output all CDS for each species.
e.g. gffread input.gff -g /path/to/genome/fasta -x output_cds.fa

#####ALL GENOMES:
Output cds sequences are provided in: cds_fastas/

Assembling multi-species input fastas for HOG cds:
Script: write_prank_input_fastas.pl to compile the nucleotide CDS sequence fastas.
Split into 12 batches: listed in text file HOG_cds_batches

Then, run prank v. 150803 with codon model to align CDS sequences for each HOG.
Script: prank_cds_alignment_template.pl
Output to: PAML_aligned_PRANK_fastas/
NOTE: batch4 HOG2_2081, batch9 HOG2_557, & batch9 HOG2_58 all failed to complete in a reasonable timeframe & were excluded

HOG CDS alignments were filtered with the Jarvis et al. 2014 DNA filtering script to mask over 'single' (e.g. species-specific regions) and low-identity regions.
Script: filter_alignment_fasta_v1.3B.pl (acessed from ftp://climb.genomics.cn/pub/10.5524.101001_102000/10/041/ on Sept. 30, 2015).
Output to: PAML_filtered_PRANK_fastas/ (contains the filtered alignments, plus the output 'stats' file)

An additional step was added to remove gap-only (and non-ACGT only) codons prior to tree-building.
(i.e. removal of alignment columns for codons with no ACGT bases, retaining reading frame)
Output to: PAML_final_PRANK_fastas/

Both the Jarvis et al. filtering & final column removal were performed with the perl wrapper script: filter_prank_alignments.pl

---------------

Step 4: Build phylogenetic tree for each HOG, to be used as guide tree in PAML analysis

Trees were built with RAxML v. 8.1.4, with 200 rapid bootstraps followed by a thorough ML tree search.  All loci were run with a GTR+GAMMA model, and with
alignments partitioned by codon1+codon2, and codon3 (e.g. 3rd codon position specified by a separate partition).

Script: oma_RAxML_template.pl (perl wrapper for RAxML runs)
Output: RAxML_output/
For each locus:
-partition file for running RAxML
-RAxML_bestTree (best ML tree inference)
-RAxML_bipartitions (bootstrap support values drawn on best ML tree)
-RAxML_bipartitionsBranchLabels (as above, but with values migrated to branches)
-RAxML_bootstrap (output from 200 rapid bootstrap replicates)
-RAxML_info
