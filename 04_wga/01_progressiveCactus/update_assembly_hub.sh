module purge
source /n/sw/progressiveCactus/environment
hal2assemblyHub.py /n/holylfs/INTERNAL_REPOS/RATITES/ratite_final_20150627/ratiteAlign.hal ratiteHub8 \
	--hub=ratiteHub \
        --lod \
	--shortLabel="RatiteHub_ver8" \
	--longLabel="Comparative assembly hub for ratite genomes" \
	--email=tsackton@oeb.harvard.edu \
	--genomes=genomelist.txt \
	--rename=/n/holylfs/INTERNAL_REPOS/RATITES/assemblyHub/inputs/common_names.txt \
	--lodDir=/n/holylfs/INTERNAL_REPOS/RATITES/assemblyHub/inputs/lodFinal \
	--lodTxtFile=/n/holylfs/INTERNAL_REPOS/RATITES/assemblyHub/inputs/lod.txt \
	--twobitdir=twoBit \
	--finalBigwigDirs=GC,alignability,phyloP,GG-FL-1,GG-FL-2,GG-FL-3,GG-HL-1,GG-HL-2,GG-HL-3,RA-FL-1,RA-FL-2,RA-HL-1,RA-HL-2,ATAC_Pooled_FL,ATAC_Pooled_HL,ATAC_Pooled_KeelS1,ATAC_Pooled_KeelS10,ATAC_Pooled_PecMusS12,ATAC_Pooled_PecMusS3,ATAC_Pooled_SternumS11,ATAC_Pooled_SternumS2 \
	--finalBigBedDirs=ATACseqHindlimb,ATACseqForelimb,final_beds/MAKER,final_beds/RNA,final_beds/CE,final_beds/Ensembl,final_beds/Forelimb,final_beds/Hindlimb,final_beds/HOXD13,final_beds/NCBI,final_beds/PITX1 \
	--maxThreads 32 


