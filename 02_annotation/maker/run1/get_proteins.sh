#!/bin/bash

wget -O musmus.fa.gz ftp://ftp.ncbi.nlm.nih.gov/genomes/M_musculus/protein/protein.fa.gz
wget -O homsap.fa.gz ftp://ftp.ncbi.nlm.nih.gov/genomes/H_sapiens/protein/protein.fa.gz
wget -O danrer.fa.gz ftp://ftp.ncbi.nlm.nih.gov/genomes/D_rerio/protein/protein.fa.gz
wget -O taegut.fa.gz ftp://ftp.ncbi.nlm.nih.gov/genomes/Taeniopygia_guttata/protein/protein.fa.gz
wget -O galgal.fa.gz ftp://ftp.ncbi.nlm.nih.gov/genomes/Gallus_gallus/protein/protein.fa.gz
wget -O melgal.fa.gz ftp://ftp.ncbi.nlm.nih.gov/genomes/Meleagris_gallopavo/protein/protein.fa.gz
wget -O ficalb.fa.gz ftp://ftp.ncbi.nlm.nih.gov/genomes/Ficedula_albicollis/protein/protein.fa.gz
wget -O strcam.fa.gz ftp://ftp.ncbi.nlm.nih.gov/genomes/Struthio_camelus_australis/protein/protein.fa.gz
wget -O tingut.fa.gz ftp://ftp.ncbi.nlm.nih.gov/genomes/Tinamus_guttatus/protein/protein.fa.gz
wget -O anocar.fa.gz ftp://ftp.ncbi.nlm.nih.gov/genomes/Anolis_carolinensis/protein/protein.fa.gz
cat *.fa.gz > maker_protevi_othersp.fa.gz
gunzip maker_protevi_othersp.fa.gz 	
rm *.fa.gz

