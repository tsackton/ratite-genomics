#!/bin/bash

NEUTVER=$1
INPUT=$2
NEUTMOD=/n/regal/edwards_lab/ratites/wga/phast/moa/phyloFit/${NEUTVER}_final.named.mod
mkdir -p $NEUTVER
phyloP --method LRT --features $INPUT.gff --mode ACC --branch rheAme,rhePen,rheAme-rhePen,strCam,casCas,droNov,casCas-droNov,aptHaa-casCas,aptHaa-aptOwe,aptRow,aptHaa,aptOwe,aptHaa-aptRow,anoDid $NEUTMOD $INPUT.fa > $NEUTVER/allRatite.out &
phyloP --method LRT --features $INPUT.gff --mode ACC --branch rheAme $NEUTMOD $INPUT.fa > $NEUTVER/rheAme.out &
phyloP --method LRT --features $INPUT.gff --mode ACC --branch rhePen $NEUTMOD $INPUT.fa > $NEUTVER/rhePen.out &
phyloP --method LRT --features $INPUT.gff --mode ACC --branch strCam $NEUTMOD $INPUT.fa > $NEUTVER/strCam.out &
phyloP --method LRT --features $INPUT.gff --mode ACC --branch casCas $NEUTMOD $INPUT.fa > $NEUTVER/casCas.out &
phyloP --method LRT --features $INPUT.gff --mode ACC --branch droNov $NEUTMOD $INPUT.fa > $NEUTVER/droNov.out &
phyloP --method LRT --features $INPUT.gff --mode ACC --branch aptRow $NEUTMOD $INPUT.fa > $NEUTVER/aptRow.out &
phyloP --method LRT --features $INPUT.gff --mode ACC --branch aptHaa $NEUTMOD $INPUT.fa > $NEUTVER/aptHaa.out &
phyloP --method LRT --features $INPUT.gff --mode ACC --branch aptOwe $NEUTMOD $INPUT.fa > $NEUTVER/aptOwe.out &
phyloP --method LRT --features $INPUT.gff --mode ACC --branch anoDid $NEUTMOD $INPUT.fa > $NEUTVER/anoDid.out &
phyloP --method LRT --features $INPUT.gff --mode ACC --branch rheAme,rhePen,rheAme-rhePen $NEUTMOD $INPUT.fa > $NEUTVER/Rhea.out &
phyloP --method LRT --features $INPUT.gff --mode ACC --branch casCas,droNov,casCas-droNov $NEUTMOD $INPUT.fa > $NEUTVER/Casuar.out &
phyloP --method LRT --features $INPUT.gff --mode ACC --branch aptHaa-casCas,aptHaa-aptOwe,aptRow,aptHaa,aptOwe,aptHaa-aptRow $NEUTMOD $INPUT.fa > $NEUTVER/Kiwi.out &
phyloP --method LRT --features $INPUT.gff --mode ACC --branch cryCin,tinGut,cryCin-tinGut,eudEle,notPer,eudEle-notPer,cryCin-eudEle $NEUTMOD $INPUT.fa > $NEUTVER/Tinamou.out &
phyloP --method LRT --features $INPUT.gff --mode ACC --branch aptHaa-strCam $NEUTMOD $INPUT.fa > $NEUTVER/basalPaleo.out &


