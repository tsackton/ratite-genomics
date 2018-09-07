Final code to prepare data for NCBI submission. Masks are already made, but starts from masks, gff files, and original fasta files.

Working directory: /n/holylfs/LABS/edwards_lab/tsackton/RATITE_PAPER_DATA_FREEZE/FINAL_RELEASE/02_annotation/final_ncbi_submission
Environment: source ~/py2_env.sh 

Copy files:

```bash
cp -r ../gag/masks .
cp -r ../../01_assembly/assemblies_ver1/*.fa.gz .
gunzip *.gz
mkdir fastas
mv *.fa fastas
cp -r ../gag/gffs/ .
```

Copy code and metadata:

```bash
cp ../gag/biosamples .
cp ../gag/generic.sbt .
cp ../gag/locus_tags .
cp ../gag/table2asn .
```

Check to make sure there are no scaffolds in GFFs that are missing from genomes:

```bash
for GENOME in aptHaa aptOwe eudEle rheAme aptRow casCas cryCin droNov notPer rhePen;
do
   grep -v "^#" gffs/$GENOME.genome.gff | cut -f1,1 | sort | uniq > $GENOME.gff.scaffolds
   grep ">" fastas/$GENOME.fa | perl -p -e 's/>//' | sort > $GENOME.fasta.scaffolds
   comm -23 $GENOME.gff.scaffolds $GENOME.fasta.scaffolds > $GENOME.missing
done
```

[tsackton@bioinf01 final_ncbi_submission]$ head *.missing
==> aptHaa.missing <==

==> aptOwe.missing <==

==> aptRow.missing <==

==> casCas.missing <==

==> cryCin.missing <==

==> droNov.missing <==

==> eudEle.missing <==

==> notPer.missing <==

==> rheAme.missing <==

==> rhePen.missing <==

Looks good! Remove stuff we don't need.

```bash
rm *.missing
rm *.fasta.scaffolds
rm *.gff.scaffolds
```

Now, add mtDNA to masks since I want to just remove mtDNA from submissions for simplicity.

```bash
for GENOME in aptHaa aptOwe eudEle rheAme;
do
  bioawk -c fastx '{print $name, 0, length($seq), "mtDNA"}' ../../01_assembly/genbank_submission/processing/${GENOME}_mt.fa >> masks/${GENOME}_ncbi_mask.bed
done
```

Now run bedtools maskfasta to mask regions in bed files and produce masked genome file.

```bash
for GENOME in aptHaa aptOwe eudEle rheAme aptRow casCas cryCin droNov notPer rhePen;
do
   bedtools maskfasta -fi fastas/$GENOME.fa \
   -bed masks/${GENOME}_ncbi_mask.bed \
   -fo fastas/${GENOME}_masked.fa &
done
```

Remove sequences that are all Ns with bioAwk:

```bash
for GENOME in aptHaa aptOwe eudEle rheAme aptRow casCas cryCin droNov notPer rhePen;
do
bioawk -c fastx '{sub(/^N+$/, "", $seq); if(length($seq) > 0)print ">"$name"\n"$seq}' fastas/${GENOME}_masked.fa > fastas/${GENOME}.final.fa
done
```

Check to make sure nothing inappropriate is removed:

```bash
for GENOME in aptHaa aptOwe eudEle rheAme aptRow casCas cryCin droNov notPer rhePen;
do
	grep -v "contaminant" masks/${GENOME}_ncbi_mask.bed | cut -f1,1 | sort | uniq > $GENOME.removed.scaffolds
	grep "^>" fastas/$GENOME.fa | perl -p -e 's/>//' | sort > $GENOME.orig.fasta.scaffolds
	grep "^>" fastas/$GENOME.final.fa | perl -p -e 's/>//' | sort > $GENOME.final.fasta.scaffolds
	comm -23 $GENOME.orig.fasta.scaffolds $GENOME.final.fasta.scaffolds > $GENOME.missing.scaffolds
	comm -23 $GENOME.removed.scaffolds $GENOME.missing.scaffolds > $GENOME.check
done
head -n 1 *.check
```

Check files should be empty, and they are.

Remove stuff not needed.

```bash
rm *.check
rm *.scaffolds
```

Now remove lines from gff with scaffolds we don't want to keep, using the mask bed files

Note added 5/9: wait a minute this doesn't quite work...adding back skip empty scaffolds line

```bash
#NOT RUN
for GENOME in aptHaa aptOwe eudEle rheAme aptRow casCas cryCin droNov notPer rhePen;
do
	grep -v "contaminant" masks/${GENOME}_ncbi_mask.bed > $GENOME.removed.scaffolds.bed
	bedtools subtract -a gffs/$GENOME.genome.gff -b $GENOME.removed.scaffolds.bed -header -A > gffs/$GENOME.final.gff
done
mv *.bed masks
```

Now should be ready to run GAG: fastas/*.final.fa is the masked, final genome versions; gffs/*.final.gff is the masked, final GFF versions

```bash
for GENOME in aptHaa aptOwe eudEle rheAme aptRow casCas cryCin droNov notPer rhePen;
do
    python ~/sw/progs/GAG/gag.py -f fastas/$GENOME.final.fa -g gffs/$GENOME.final.gff -o gag_$GENOME --fix_terminal_ns -ses
done
```

Next, use gffread to discard any mRNAs with in-frame stops and rerun table2asn:

```bash
for GENOME in aptHaa aptOwe eudEle rheAme aptRow casCas cryCin droNov notPer rhePen;
do
/n/home12/tsackton/sw/progs/gffread/gffread gag_$GENOME/genome.gff -g gag_$GENOME/genome.fasta -O -F -V -H -o gag_$GENOME/genome.fixed.gff
done
```

Finally, check to make sure there are no lines in GFF that are missing from genomes:

```bash
for GENOME in aptHaa aptOwe eudEle rheAme aptRow casCas cryCin droNov notPer rhePen;
do
   grep -v "^#" gag_$GENOME/genome.fixed.gff | cut -f1,1 | sort | uniq > $GENOME.gff.scaffolds
   grep ">" gag_$GENOME/genome.fasta | perl -p -e 's/>//' | sort > $GENOME.fasta.scaffolds
   comm -23 $GENOME.gff.scaffolds $GENOME.fasta.scaffolds > $GENOME.missing
done
```

head *.missing is clean

Looks good! Remove stuff we don't need.

```bash
rm *.missing
rm *.fasta.scaffolds
rm *.gff.scaffolds
```

Now have the final fasta files (NCBI) and GFF files (Dryad, etc) to upload.

Make a submission directory and rename files
```bash
mkdir -p final_submission_fastas
mkdir -p final_submission_gffs
for GENOME in aptHaa aptOwe eudEle rheAme aptRow casCas cryCin droNov notPer rhePen;
do
  cp gag_$GENOME/genome.fasta final_submission_fastas/$GENOME.fa
  cp gag_$GENOME/genome.fixed.gff final_submission_gffs/$GENOME.gff
done
cd final_submission_fastas
for FILE in $(ls *.fa); do gzip $FILE & done
```

Upload:

```bash
 ~/.aspera/connect/bin/ascp -i /n/holylfs/LABS/edwards_lab/tsackton/RATITE_PAPER_DATA_FREEZE/NCBI_SRA/aspera.openssh.txt -QT -k1 -d final_submission_fastas subasp@upload.ncbi.nlm.nih.gov:uploads/tsackton@oeb.harvard.edu_kDOcvVqP
```

