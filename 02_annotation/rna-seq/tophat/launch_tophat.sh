#!/bin/bash

for KIWI in aptHaa aptRow aptOwe; do sbatch run_tophat.sh $KIWI reads/kiwi/kiwi.R1.pair.fastq.gz,reads/kiwi/kiwi.R1.single.fastq.gz,reads/kiwi/kiwi.R2.single.fastq.gz reads/kiwi/kiwi.R2.pair.fastq.gz; sbatch run_tophat.sh $KIWI reads/dronov/dronov.R1.pair.fastq.gz,reads/dronov/dronov.R1.single.fastq.gz,reads/dronov/dronov.R2.single.fastq.gz reads/dronov/dronov.R2.pair.fastq.gz; done
for EMU in droNov casCas; do sbatch run_tophat.sh $EMU reads/dronov/dronov.R1.pair.fastq.gz,reads/dronov/dronov.R1.single.fastq.gz,reads/dronov/dronov.R2.single.fastq.gz reads/dronov/dronov.R2.pair.fastq.gz; sbatch run_tophat.sh $EMU reads/kiwi/kiwi.R1.pair.fastq.gz,reads/kiwi/kiwi.R1.single.fastq.gz,reads/kiwi/kiwi.R2.single.fastq.gz reads/kiwi/kiwi.R2.pair.fastq.gz; done
for RHEA in rheAme rhePen; do sbatch run_tophat.sh $RHEA reads/dronov/dronov.R1.pair.fastq.gz,reads/dronov/dronov.R1.single.fastq.gz,reads/dronov/dronov.R2.single.fastq.gz reads/dronov/dronov.R2.pair.fastq.gz; sbatch run_tophat.sh $RHEA reads/kiwi/kiwi.R1.pair.fastq.gz,reads/kiwi/kiwi.R1.single.fastq.gz,reads/kiwi/kiwi.R2.single.fastq.gz reads/kiwi/kiwi.R2.pair.fastq.gz; done
for TINA in notPer cryCin eudEle; do sbatch run_tophat.sh $TINA reads/notper/notper.R1.pair.fastq.gz,reads/notper/notper.R1.single.fastq.gz,reads/notper/notper.R2.single.fastq.gz reads/notper/notper.R2.pair.fastq.gz; done

