# from Tim's email:
# The protocol:
# 1) Run the R code in the repo (https://github.com/tsackton/ratite-genomics/tree/master/07_cnee_analysis/transfac_input) to generate an input bed file. Currently makes a file with just elements with score != "none" but this is easy to modify.

# update to next step: 
# 2) Run get_seq_from_bed.sh (here, Tim's version needed a couple changes) with the first argument as the species you want to create a fasta for and the 2nd argument as the bed file you want to use (e.g. the one from the script above -- crucial though to use one generated with the script in 1 because it does the coordinate conversation between genome and CNEE alignment)

# to go from the fasta files for each species to the files in this directory, the following was called:
# 3) for SP in allMis allSin anaPla anoCar anoDid aptFor aptHaa aptOwe aptRow balReg calAnn casCas chaPel chaVoc cheMyd chrPic colLiv corBra croPor cryCin cucCan droNov eudEle falPer ficAlb fulGla galGal gavGan halLeu lepDis melGal melUnd mesUni nanAur nanBra nanHar nipNip notPer picPub pseHum pygAde rheAme rhePen strCam taeGut tinGut uriPel; do ./get_seq_from_bed.sh $SP cnee_part_acc.bed; perl -p -i -e 's/-//g' ${SP}_cnee_part_acc.bed.fa; done

# put them into < 10MB batches:
# make sure to remove all previous batch* files before running this or it appends
# 4) for f in allMis_cnee_part_acc.bed.fa allSin_cnee_part_acc.bed.fa anaPla_cnee_part_acc.bed.fa anoCar_cnee_part_acc.bed.fa anoDid_cnee_part_acc.bed.fa aptFor_cnee_part_acc.bed.fa aptHaa_cnee_part_acc.bed.fa aptOwe_cnee_part_acc.bed.fa aptRow_cnee_part_acc.bed.fa balReg_cnee_part_acc.bed.fa calAnn_cnee_part_acc.bed.fa casCas_cnee_part_acc.bed.fa chaPel_cnee_part_acc.bed.fa chaVoc_cnee_part_acc.bed.fa cheMyd_cnee_part_acc.bed.fa chrPic_cnee_part_acc.bed.fa; do (cat "${f}"; echo) >> batch1.fa; done
# 5) for f in colLiv_cnee_part_acc.bed.fa corBra_cnee_part_acc.bed.fa croPor_cnee_part_acc.bed.fa cryCin_cnee_part_acc.bed.fa cucCan_cnee_part_acc.bed.fa droNov_cnee_part_acc.bed.fa eudEle_cnee_part_acc.bed.fa falPer_cnee_part_acc.bed.fa ficAlb_cnee_part_acc.bed.fa fulGla_cnee_part_acc.bed.fa galGal_cnee_part_acc.bed.fa gavGan_cnee_part_acc.bed.fa halLeu_cnee_part_acc.bed.fa lepDis_cnee_part_acc.bed.fa melGal_cnee_part_acc.bed.fa melUnd_cnee_part_acc.bed.fa; do (cat "${f}"; echo) >> batch2.fa; done 
# 6) for f in mesUni_cnee_part_acc.bed.fa nanAur_cnee_part_acc.bed.fa nanBra_cnee_part_acc.bed.fa nanHar_cnee_part_acc.bed.fa nipNip_cnee_part_acc.bed.fa notPer_cnee_part_acc.bed.fa picPub_cnee_part_acc.bed.fa pseHum_cnee_part_acc.bed.fa pygAde_cnee_part_acc.bed.fa rheAme_cnee_part_acc.bed.fa rhePen_cnee_part_acc.bed.fa strCam_cnee_part_acc.bed.fa taeGut_cnee_part_acc.bed.fa tinGut_cnee_part_acc.bed.fa uriPel_cnee_part_acc.bed.fa; do (cat "${f}"; echo) >> batch3.fa; done

