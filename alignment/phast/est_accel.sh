#!/bin/bash

#SBATCH -t 1-00:00
#SBATCH --mem 400
#SBATCH -p serial_requeue
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -J phyloP

NEUTMOD=$1 
CHR=$2
ALIGN=$3
phyloP --method LRT --features $CHR.temp.bed --mode ACC --branch rheAme,rhePen,rheAme-rhePen,strCam,casCas,droNov,casCas-droNov,aptHaa-casCas,aptHaa-aptOwe,aptRow,aptHaa,aptOwe,aptHaa-aptRow $NEUTMOD $ALIGN > ratite/$CHR.Ratite.out
phyloP --method LRT --features $CHR.temp.bed --mode ACC --branch rheAme $NEUTMOD $ALIGN > tips/$CHR.rheAme.out
phyloP --method LRT --features $CHR.temp.bed --mode ACC --branch rhePen $NEUTMOD $ALIGN > tips/$CHR.rhePen.out
phyloP --method LRT --features $CHR.temp.bed --mode ACC --branch strCam $NEUTMOD $ALIGN > tips/$CHR.strCam.out
phyloP --method LRT --features $CHR.temp.bed --mode ACC --branch casCas $NEUTMOD $ALIGN > tips/$CHR.casCas.out
phyloP --method LRT --features $CHR.temp.bed --mode ACC --branch droNov $NEUTMOD $ALIGN > tips/$CHR.droNov.out
phyloP --method LRT --features $CHR.temp.bed --mode ACC --branch aptRow $NEUTMOD $ALIGN > tips/$CHR.aptRow.out
phyloP --method LRT --features $CHR.temp.bed --mode ACC --branch aptHaa $NEUTMOD $ALIGN > tips/$CHR.aptHaa.out
phyloP --method LRT --features $CHR.temp.bed --mode ACC --branch aptOwe $NEUTMOD $ALIGN > tips/$CHR.aptOwe.out
phyloP --method LRT --features $CHR.temp.bed --mode ACC --branch rheAme,rhePen,rheAme-rhePen $NEUTMOD $ALIGN > clades/$CHR.Rhea.out
phyloP --method LRT --features $CHR.temp.bed --mode ACC --branch casCas,droNov,casCas-droNov $NEUTMOD $ALIGN > clades/$CHR.Casuar.out
phyloP --method LRT --features $CHR.temp.bed --mode ACC --branch aptHaa-casCas,aptHaa-aptOwe,aptRow,aptHaa,aptOwe,aptHaa-aptRow $NEUTMOD $ALIGN > clades/$CHR.Kiwi.out
phyloP --method LRT --features $CHR.temp.bed --mode ACC --branch cryCin,tinGut,cryCin-tinGut,eudEle,notPer,eudEle-notPer,cryCin-eudEle $NEUTMOD $ALIGN > tinamou/$CHR.Tinamou.out		
