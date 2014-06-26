#!/bin/bash

#script to fetch existing data from SRA for selected initial species
#written for the Harvard Odyssey cluster module management system and SLURM

#get aspera
module load hpc/aspera

#set up directories
mkdir -p /n/regal/edwards_lab/ratites/data/sra/strCam
mkdir -p /n/regal/edwards_lab/ratites/data/sra/tinGut
mkdir -p /n/regal/edwards_lab/ratites/data/sra/melUnd

#fetch ostrich
ascp -i /n/sw/aspera/etc/asperaweb_id_dsa.putty -QTr -l 250M anonftp@ftp-private.ncbi.nlm.nih.gov:/sra/sra-instant/reads/ByStudy/sra/SRP/SRP028/SRP028745/ /n/regal/edwards_lab/ratites/data/sra/strCam

#fetch tinamou
ascp -i /n/sw/aspera/etc/asperaweb_id_dsa.putty -QTr -l 250M anonftp@ftp-private.ncbi.nlm.nih.gov:/sra/sra-instant/reads/ByStudy/sra/SRP/SRP028/SRP028753/ /n/regal/edwards_lab/ratites/data/sra/tinGut

#fetch budgie
ascp -i /n/sw/aspera/etc/asperaweb_id_dsa.putty -QTr -l 250M anonftp@ftp-private.ncbi.nlm.nih.gov:/sra/sra-instant/reads/ByStudy/sra/ERP/ERP002/ERP002324/ /n/regal/edwards_lab/ratites/data/sra/melUnd
