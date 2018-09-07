README: MOA CONSERVED ELEMENTS
===============

PURPOSE: extract sequence from the little bush moa reference-based genome assembly for a subset of 
conserved non-exonic elements (CNEEs) identified from the palaeognath project whole-genome alignment (wga), 
and align with sequences from all other speces in wga.
Additionally, compile an alignment of fourfould degenerate (4d) sites for all species including moa.

Note that scripts were written for the Odyssey Cluster at Harvard University, and therefore might include code
to account for program stops/restarts on Odyssey's serial_requeue partition, and/or hard-coded paths
to e.g. software installations and user directories.

===============

CNEEs:
Step 1: Lift over coordinates from chicken (galGal) reference to each other species in wga using 
halLiftover v2.1 (from HAL Tools Hickey et al. 2013 Bioinformatics 29:1341-1342, https://github.com/glennhickey/hal).

Input file of galGal CNEE coordinates: final_cnees_long.bed from Tim Sackton lists coordinates for 284,001 cnees.

Script: preprocess_beds_for_moa.sh (by Tim Sackton, archived at:  https://github.com/tsackton/ratite-genomics/blob/master/07_cnee_analysis/01_phyloP/)

Output: .psl format file for each species (e.g. droNov.psl)

---------------

Step 2: Parse liftover .psl files.  Omit CNEE from target species if there are multiple liftover target regions (e.g. a single chicken
reference sequence lifts over to multiple regions in target species), or if the target region is > 2X the reference chicken length.  
Also note if a particular CNEE is 'missing' (no target liftover region), or is 'short' (< 50% of chicken referece length), but do not 
omit short CNEEs.

Script: parse_cnee_halLiftover.pl
Output: a .bed format file for each species, (e.g. droNov_cnees_parsed_liftover.bed)
Additionally, a logfile (final_cnees_long_liftover_parsing_log) that records missing/multiple/long/short CNEEs in each target species

---------------

Step 3: Compile .fasta file for each CNEE (e.g. sequence for that CNEE from all species).

Script: write_de_novo_align_fastas.pl
Output: input_fastas_allspecies_cnees.tar.gz (split into 57 subdirectory 'batches' of ~5000 loci; text file 'cnee_alignment_batches' lists CNEEs
contained in each batch)

---------------

Step 4: Add sequence for little bush moa (Anomalopteryx didiformis [anoDid], a reference-based genome assembly mapping reads to the emu [droNov] genome).
   
This approach used the emu (droNov) CNEE coordinates lifted over from the chicken (galGal) 
reference, and mapped them to their corresponding moa coordinates using whole-scaffold emu-moa
alignments produced during the creation of the moa genome assembly.  The moa genome assembly is 
described in Cloutier et al. 2018 (Bioarchive doi:https://doi.org/10.1101/262816 ), with the nuclear 
genome deposited under NCBI accession (pending), and accompanying data (including emu-moa scaffold 
alignments) available under Dryad digital repository DOI:(pending).  

---------------

Step 5: Align individual CNEEs

Script: mafft_cnee_align_template.pl
Align with mafft v. 7.245 (Katoh and Standley 2013, Mol Biol Evol 30:772-780) and default options 
for high sensitivity alignment with global pairwise alignment during the consistency scoring step ('ginsi' option)

Output: input_fastas_allspecies_cnees_aligned.tar.gz [subdirectory 'batches' as listed in cnee_alignment_batches]

---------------

Step 6: Remove alignment columns with gap in reference chicken (galGal) sequence

This step was done to make the resulting downstream analyses of CNEEs comparable to those done without moa
(e.g. done from wga with galGal as the reference).

Script: remove_galGal_gaps_from_cnee_alignments.pl
Output: input_fastas_allspecies_cnees_aligned_no-galGal-gaps.tar.gz [subdirectory 'batches' as listed in cnee_alignment_batches]

---------------

Step 7: Concatenate CNEEs into a single file

Script: write_cnee_spp_concatenated_seq.pl
Concatenates all CNEE sequences for a given species into a single fasta.  Concatenates loci in consistent manner across species,
and adds gaps (-) for missing loci within a species, and Ns for omitted (multiple/too long) loci to allow single-species fastas to 
be combined into final output (see below).

Script: write_final_fasta_allspecies.pl
Combines all single-species fastas into a single file & outputs in both fasta and phylip format
Output: allspecies_cnee_concatenated.fasta & allspecies_cnee_concatenated.phy

Script: write_cnee_concat_partitionfile.pl
Outputs a text file listing the partitions in the concatenated CNEE file, in both RAxML-style and .bed format
Output: allspecies_cnee_concat_partitions & allspecies_cnee_concat_partitions.bed

---------------

Step 8: Per-partition branch length estimation

Estimate branch lengths for each CNEE, constrained to each of 3 topologies (allowing placement of rheas to differ).
[Constraint trees: ratiteTree_WithMoa.ver1.nh, ratiteTree_WithMoa.ver2.nh, ratiteTree_WithMoa.ver3.nh from Tim Sackton]

Scripts: write_branchlength_concat_alignment.pl & write_branchlength_concat_partitionfile.pl
Concatenates CNEE loci in sets of ~126 loci [because max. for default RAxML compilation per-partition branchlength
estimation is 128].

Script: RAxML_branchlengths_template.pl
Perl wrapper for RAxML, run for each of the 3 input tree topologies.  Outputs files to (non-archived) directory structure,
with a subdirectory for each RAxML run of 126 loci.

Script: parse_RAxML_branchlength_output.pl
Output: cnee_branchlengths_ver1.tar.gz, cnee_branchlengths_ver2.tar.gz, cnee_branchlengths_ver3.tar.gz
Top-level directory contains a subdirectory for each batch (corresponding to cnee_alignment_batches), with
RAxML output files renamed to include locus ID (instead of e.g. PARTITION.0).  Also included is e.g. RAxML_info.cnee_branchlengths_ver1
with all model parameters across runs compiled into a single file.

===============

FOURFOLD DEGENERATE (4d) SITES:
Used a similar approach to compilation of CNEE dataset above, but no alignment step was necessary.

---------------

Step 1: 
Input file of galGal 4d coordinates: galGal4_4d_bed6_withid.bed from Tim Sackton lists coordinates for 5,739,749 4d sites

Script: preprocess_beds_for_moa.sh (by Tim Sackton, archived at:  https://github.com/tsackton/ratite-genomics/blob/master/07_cnee_analysis/01_phyloP/)

Output: .bed format file for each species (e.g. droNov.bed)

---------------

Step 2: Parse liftover .bed files

Script: parse_neutMods_4d_bed.pl
Outputs a .bed file for each (non-moa) target species, with an entry for each chicken reference base.  
Placeholder lines are output for missing or multiple liftover sites.

---------------

Step 3: Concatenate 4d sites into a single file

Script: write_4d_spp_concatenated_seq.pl
Uses parsed .bed annotation files to retrieve 4d sites, and concatenates sites for a given species into a single fasta.  
Concatenates sites in consistent manner across species, and adds gaps (-) for missing sites within a species, and Ns for 
omitted sites (liftover to multiple bases) to allow single-species fastas to be combined into final output (see below).

Script: write_final_fasta_allspecies.pl
Combines all single-species fastas into a single file & outputs in both fasta and phylip format
Output: allspecies_4d_concatenated.fasta & allspecies_4d_concatenated.phy

