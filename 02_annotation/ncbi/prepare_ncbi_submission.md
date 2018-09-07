Goal is to use GAG [http://genomeannotation.github.io/GAG/] to generate inputs to NCBI.

```bash
source ~/py2_env.sh
```

First, copy masks to masks subdir:

```bash
cd masks
cp ../../../01_assembly/genbank_submission/ver1_masks/*.bed .
cd ..
```

Next, use bedtools maskfasta to replace masks with Ns in mt-labeled but otherwise unmodified fastas

```bash
for GENOME in aptHaa aptOwe eudEle rheAme aptRow casCas cryCin droNov notPer rhePen;
do
   bedtools maskfasta -fi /n/holylfs/LABS/edwards_lab/tsackton/RATITE_PAPER_DATA_FREEZE/FINAL_RELEASE/01_assembly/genbank_submission/processing/${GENOME}_mt.fa \
   -bed masks/${GENOME}_ncbi_mask.bed -fo fastas/${GENOME}_mt.fa &
   bedtools maskfasta -fi /n/holylfs/LABS/edwards_lab/tsackton/RATITE_PAPER_DATA_FREEZE/FINAL_RELEASE/01_assembly/genbank_submission/processing/${GENOME}_nt.fa \
   -bed masks/${GENOME}_ncbi_mask.bed -fo fastas/${GENOME}_nt.fa &
done
```

Merge nt and mt
```bash
for GENOME in aptHaa aptOwe eudEle rheAme aptRow casCas cryCin droNov notPer rhePen;
do
   cat fastas/$GENOME*.fa > fastas/$GENOME.masked
done
```

Now copy gffs to gff subdir:
```bash
cp ../final_gffs/* gffs
```

Now test gag on a genome. Goal is to use GAG to produce fasta / GFF than can then be fed to NCBI 

```bash
for GENOME in aptHaa;
do
    python ~/sw/progs/GAG/gag.py -f fastas/$GENOME.masked -g gffs/$GENOME.genome.gff -o out_$GENOME --fix_terminal_ns -ses
done
```

GAG fails on completely masked scaffolds, so remove things that are all Ns. Use bioawk

```bash
for GENOME in aptHaa aptOwe eudEle rheAme aptRow casCas cryCin droNov notPer rhePen;
do
bioawk -c fastx '{sub(/^N+$/, "", $seq); if(length($seq) > 999)print ">"$name"\n"$seq}' fastas/${GENOME}.masked > fastas/${GENOME}.final.fa
done
```

Retry GAG:

```bash
for GENOME in aptHaa;
do
    python ~/sw/progs/GAG/gag.py -f fastas/$GENOME.final.fa -g gffs/$GENOME.genome.gff -o out_$GENOME --fix_terminal_ns -ses
done
```

Run is successful, now need to test feeding into tbl2asn with GFF to see what errors I get that need to be fixed.

Make template file on NCBI website, then run table2asn with this command as a test:
```bash
 ../table2asn -M n -c wsf -euk -t aptHaa.sbt -a r20k -l paired-ends -i genome.fasta -f genome.gff -o aptHaa.sqn -Z aptHaa.dr -locus-tag-prefix C5126
 ```

Some confusing stuff, trying again with a few additional options

```bash
 ../table2asn -M n -J -j "[organism=Aptex haasti]" -c wsf -euk -t aptHaa.sbt -gaps-min 20 -l paired-ends -i genome.fasta -f genome.gff -o aptHaa.sqn -Z aptHaa.dr -locus-tag-prefix C5126
 ```

Need to deal with:
14 ERROR-level messages exist

SEQ_INST.StopInProtein:	7
SEQ_FEAT.InternalStop:	7

Try GAG again with --fix_start_stop flag

```bash
for GENOME in aptHaa;
do
    python ~/sw/progs/GAG/gag.py -f fastas/$GENOME.final.fa -g gffs/$GENOME.genome.gff -o out_$GENOME --fix_terminal_ns -ses --fix_start_stop
done
```

Now run table2asn again:

```bash
 ../table2asn -M n -J -j "[organism=Aptex haasti]" -c wsf -euk -t aptHaa.sbt -gaps-min 20 -l paired-ends -i genome.fasta -f genome.gff -o aptHaa.sqn -Z aptHaa.dr -locus-tag-prefix C5126
```

Better without the --fix_start_stop

Now run GAG on all genomes:

```bash
for GENOME in aptHaa aptOwe eudEle rheAme aptRow casCas cryCin droNov notPer rhePen;
do
    python ~/sw/progs/GAG/gag.py -f fastas/$GENOME.final.fa -g gffs/$GENOME.genome.gff -o out_$GENOME --fix_terminal_ns -ses
done
```

Finally, need do three things before running table2asn:
1. Make template file (update biosample)
2. Assign locus-tag-prefix
3. Assign organism tag

This will make initial versions of everything, then need to look at and deal with errors

```bash
for GENOME in aptHaa aptOwe eudEle rheAme aptRow casCas cryCin droNov notPer rhePen;
do
cd out_$GENOME
BIOSAM=$(grep "$GENOME" ../biosamples | cut -f1,1)
cp ../generic.sbt $GENOME.sbt
perl -p -i -e "s/BIOSAMPLE/$BIOSAM/" $GENOME.sbt
PREFIX=$(grep "$BIOSAM" ../locus_tags | cut -f1,1)
ORG=$(grep "$GENOME" ../biosamples | cut -f3,3)
../table2asn -M n -J -j "[organism=$ORG]" -c wsf -euk -t $GENOME.sbt -gaps-min 20 -l paired-ends -i genome.fasta -f genome.gff -o $GENOME.sqn -Z $GENOME.dr -locus-tag-prefix $PREFIX &> $GENOME.table2asn.log &
cd ..
done
```

Some errors, lets see if there are differences with tbl2asn instead of table2asn using GFF. 

Testing with aptHaa:

```bash
tbl2asn -p . -t aptHaa.sbt -M n -Z discrep -a r20k -l paired-ends -j "[organism=Apteryx haastii]"
```

Slow and no locus_tag option so won't use.

Next, use gffread to discard any mRNAs with in-frame stops and rerun table2asn:

```bash
for GENOME in aptHaa aptOwe eudEle rheAme aptRow casCas cryCin droNov notPer rhePen;
do
/n/home12/tsackton/sw/progs/gffread/gffread out_$GENOME/genome.gff -g out_$GENOME/genome.fasta -O -F -V -H -o out_$GENOME/genome.fixed.gff
done
```

Try table2asn again:

```bash
for GENOME in aptHaa aptOwe eudEle rheAme aptRow casCas cryCin droNov notPer rhePen;
do
cd out_$GENOME
BIOSAM=$(grep "$GENOME" ../biosamples | cut -f1,1)
cp ../generic.sbt $GENOME.sbt
perl -p -i -e "s/BIOSAMPLE/$BIOSAM/" $GENOME.sbt
PREFIX=$(grep "$BIOSAM" ../locus_tags | cut -f1,1)
ORG=$(grep "$GENOME" ../biosamples | cut -f3,3)
mv $GENOME.sqn $GENOME.run1.sqn 
mv $GENOME.dr $GENOME.run1.dr 
mv $GENOME.stats $GENOME.run1.stats
mv $GENOME.val $GENOME.run1.val
../table2asn -M n -J -j "[organism=$ORG]" -c wsf -euk -t $GENOME.sbt -gaps-min 20 -l paired-ends -i genome.fasta -f genome.fixed.gff -o $GENOME.sqn -Z $GENOME.dr -locus-tag-prefix $PREFIX &> $GENOME.table2asn.2.log &
cd ..
done
```

Something is messed up, there are scaffolds that are dropping out from at least some genomes which is then I think causing problems with table2asn (although I'm not sure why).

Think I will have to go back to the beginning? 

Argh. Will tackle next week. Tasks added to Things.
