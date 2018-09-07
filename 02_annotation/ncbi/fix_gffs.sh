#make coordinate masks to update GFF files
#need two things: 
#first a bed file of all scaffolds to remove entirely

#make fai index file
cp ../../01_assembly/assemblies_ver1/*.fa.gz .
gunzip *.fa.gz
for FASTA in $(ls *.fa);
do
 samtools faidx $FASTA
done

#second a way to adjust coordinates of scaffolds that are not removed but have start trimmed
#finally remove any genes with features entirely contained in a masked region


#this will need to be a python script unfortunately. pseudocode:

function get_lengths:
    read in fai file
    parse to get scaffold id and length

function get_masked_regions:
    read in bed file
    for each entry get scaffold start stop length

function remove_scaffolds:
    take length and masked_region dicts
    for each masked_region, ask if length == scaffold length
    if yes add to list of scaffolds to remove
    if no next

function coordinate_shift:
   take masked_region dict
   if start == 0 that means begins at beginning of contig
   calculate coordinate shift, which will be equal to length of masked region with 0 start
   (e.g., if masked region is 0 93, length is 93, which means position 94 (1 base) is now 1
   make coord shift dict that gives value to subtract for each scaffold

function parse_gff:
    read in gff file
    for each line, keep most fields, however need to extract: scaffold, start, end, associated gene
    add to gff dict
    masked_gene_dict = masked_gene_check(line, masked_regions)

function masked_gene_check:
    take gff, masked regions as input
    if input type is CDS or exon, and input start and stop are entirely within a masked region
    add gene id to masked_gene_dict
    otherwise nothing

function check_gff:
    take a gff line from the gff dict
    first check to see if scaffold is in sacffold_remove
    if yes, skip
    if no, next check to see if coordinate needs to be shifted
       if yes, update coordinate and continue
       if no, leave coordinate and continue
    finally check to see if feature ID (parent or ID) is in the 'bad features to remove list'
    if yes, skip
    if no, print out line (potentially with updated coordinates)
    
function main:
    get species from command line
    lengths = get_lengths(species)
    masked = get_masked_region(species)
    remove = remove_scaffolds(masked, lengths)
    coordinate_shift = coordinate_shift(masked)
    gff_dict,masked_genes = parse_gff(species)
    for line in gff_dict:
        check_gff(line)
    
    
    
