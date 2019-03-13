- RARs were associated to their nearest gene using bedtools closest with default settings
- The galGal4 gene locations and names were pulled from the GFF available here: 
http://mirrors.vbi.vt.edu/mirrors/ftp.ncbi.nih.gov/genomes/refseq/vertebrate_other/Gallus_gallus/all_assembly_versions/GCF_000002315.3_Gallus_gallus-4.0/
- that gcf file was grep'd with the following to only link coding genes to the elements: grep "gene_biotype=protein_coding"
- Given that there were only 54 elements, they were all verified as the closest coding genes in the genome browser. 
