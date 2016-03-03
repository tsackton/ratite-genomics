#!/bin/bash

NEUTVER=$1
NEUTMOD=/n/regal/edwards_lab/ratites/wga/phast/moa/phyloFit/${NEUTVER}_final.named.mod
mkdir -p $NEUTVER
phyloP --method LRT --features cnees.bed --mode ACC --branch rheAme,rhePen,rheAme-rhePen,strCam,casCas,droNov,casCas-droNov,aptHaa-casCas,aptHaa-aptOwe,aptRow,aptHaa,aptOwe,aptHaa-aptRow,anoDid $NEUTMOD allspecies_cnee_concatenated.fasta > $NEUTVER/allRatite.out &
phyloP --method LRT --features cnees.bed --mode ACC --branch rheAme $NEUTMOD allspecies_cnee_concatenated.fasta > $NEUTVER/rheAme.out &
phyloP --method LRT --features cnees.bed --mode ACC --branch rhePen $NEUTMOD allspecies_cnee_concatenated.fasta > $NEUTVER/rhePen.out &
phyloP --method LRT --features cnees.bed --mode ACC --branch strCam $NEUTMOD allspecies_cnee_concatenated.fasta > $NEUTVER/strCam.out &
phyloP --method LRT --features cnees.bed --mode ACC --branch casCas $NEUTMOD allspecies_cnee_concatenated.fasta > $NEUTVER/casCas.out &
phyloP --method LRT --features cnees.bed --mode ACC --branch droNov $NEUTMOD allspecies_cnee_concatenated.fasta > $NEUTVER/droNov.out &
phyloP --method LRT --features cnees.bed --mode ACC --branch aptRow $NEUTMOD allspecies_cnee_concatenated.fasta > $NEUTVER/aptRow.out &
phyloP --method LRT --features cnees.bed --mode ACC --branch aptHaa $NEUTMOD allspecies_cnee_concatenated.fasta > $NEUTVER/aptHaa.out &
phyloP --method LRT --features cnees.bed --mode ACC --branch aptOwe $NEUTMOD allspecies_cnee_concatenated.fasta > $NEUTVER/aptOwe.out &
phyloP --method LRT --features cnees.bed --mode ACC --branch anoDid $NEUTMOD allspecies_cnee_concatenated.fasta > $NEUTVER/anoDid.out &
phyloP --method LRT --features cnees.bed --mode ACC --branch rheAme,rhePen,rheAme-rhePen $NEUTMOD allspecies_cnee_concatenated.fasta > $NEUTVER/Rhea.out &
phyloP --method LRT --features cnees.bed --mode ACC --branch casCas,droNov,casCas-droNov $NEUTMOD allspecies_cnee_concatenated.fasta > $NEUTVER/Casuar.out &
phyloP --method LRT --features cnees.bed --mode ACC --branch aptHaa-casCas,aptHaa-aptOwe,aptRow,aptHaa,aptOwe,aptHaa-aptRow $NEUTMOD allspecies_cnee_concatenated.fasta > $NEUTVER/Kiwi.out &
phyloP --method LRT --features cnees.bed --mode ACC --branch cryCin,tinGut,cryCin-tinGut,eudEle,notPer,eudEle-notPer,cryCin-eudEle $NEUTMOD allspecies_cnee_concatenated.fasta > $NEUTVER/Tinamou.out &

