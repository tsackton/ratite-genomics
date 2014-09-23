#!/bin/bash

#cd /scratch/tsackton/assemblies/tinGut/cegma
#cegma -g ../BGI/GCA_000705375.2_ASM70537v2_genomic.fna -o tinGut_BGI -T 8 --vrt &> ../logs/cegma_20140905.log
cd /scratch/tsackton/assemblies/strCam/cegma
cegma -g ../Zhou/Struthio_camelus.20130116.OM.fa -o strCam_Zhou -T 8 --vrt &> ../logs/cegma_2010908.log
#cd /scratch/tsackton/assemblies/galGal/cegma
#cegma -g ../v4/GCF_000002315.3_Gallus_gallus-4.0_genomic.fna -o galGal_v4 -T 8 --vrt &> ../logs/cegma_20140905.log

