#clean up files

zcat 140808_D00365_0310_AC5KY1ANXX/BclToFastq_Lane1_Indexlength6_Run1/Project_C5KY1ANXX/Sample_emu_220_CTCTAC/*.fastq.gz | grep -A 3 '^@.* [^:]*:N:[^:]*:' | grep -v "^--$" > droNov_R2_220.fastq
zcat 140808_D00365_0310_AC5KY1ANXX/BclToFastq_Lane1_Indexlength6_Run1/Project_C5KY1ANXX/Sample_emu_3kb_ACAGTG/*.fastq.gz | grep -A 3 '^@.* [^:]*:N:[^:]*:' | grep -v "^--$" > droNov_R2_3kb.fastq
zcat 140808_D00365_0310_AC5KY1ANXX/BclToFastq_Lane1_Indexlength6_Run1/Project_C5KY1ANXX/Sample_g_rhea_220_ATCAGT/*.fastq.gz | grep -A 3 '^@.* [^:]*:N:[^:]*:' | grep -v "^--$" > rheAme_R2_220.fastq
zcat 140808_D00365_0310_AC5KY1ANXX/BclToFastq_Lane1_Indexlength6_Run1/Project_C5KY1ANXX/Sample_g_rhea_3kb_GTGAAA/*.fastq.gz | grep -A 3 '^@.* [^:]*:N:[^:]*:' | grep -v "^--$" > rheAme_R2_3kb.fastq
zcat 140808_D00365_0310_AC5KY1ANXX/BclToFastq_Lane1_Indexlength6_Run1/Project_C5KY1ANXX/Sample_ls_kiwi_220_TAGTTG/*.fastq.gz | grep -A 3 '^@.* [^:]*:N:[^:]*:' | grep -v "^--$" > aptOwe_R2_220.fastq
zcat 140808_D00365_0310_AC5KY1ANXX/BclToFastq_Lane1_Indexlength6_Run1/Project_C5KY1ANXX/Sample_ls_kiwi_3kb_CTTGTA/*.fastq.gz | grep -A 3 '^@.* [^:]*:N:[^:]*:' | grep -v "^--$" > aptOwe_R2_3kb.fastq
zcat 140808_D00365_0310_AC5KY1ANXX/BclToFastq_Lane1_Indexlength6_Run1/Project_C5KY1ANXX/Sample_s_cassowary_220_GCTACC/*.fastq.gz | grep -A 3 '^@.* [^:]*:N:[^:]*:' | grep -v "^--$" > casCas_R2_220.fastq
zcat 140808_D00365_0310_AC5KY1ANXX/BclToFastq_Lane1_Indexlength6_Run1/Project_C5KY1ANXX/Sample_s_cassowary_3kb_GCCAAT/*.fastq.gz | grep -A 3 '^@.* [^:]*:N:[^:]*:' | grep -v "^--$" > casCas_R2_3kb.fastq

zcat 140808_D00365_0310_AC5KY1ANXX/BclToFastq_Lane2_Indexlength6_Run1/Project_C5KY1ANXX/Sample_emu_220_CTCTAC/*.fastq.gz | grep -A 3 '^@.* [^:]*:N:[^:]*:' | grep -v "^--$" >> droNov_R2_220.fastq
zcat 140808_D00365_0310_AC5KY1ANXX/BclToFastq_Lane2_Indexlength6_Run1/Project_C5KY1ANXX/Sample_emu_3kb_ACAGTG/*.fastq.gz | grep -A 3 '^@.* [^:]*:N:[^:]*:' | grep -v "^--$" >> droNov_R2_3kb.fastq
zcat 140808_D00365_0310_AC5KY1ANXX/BclToFastq_Lane2_Indexlength6_Run1/Project_C5KY1ANXX/Sample_g_rhea_220_ATCAGT/*.fastq.gz | grep -A 3 '^@.* [^:]*:N:[^:]*:' | grep -v "^--$" >> rheAme_R2_220.fastq
zcat 140808_D00365_0310_AC5KY1ANXX/BclToFastq_Lane2_Indexlength6_Run1/Project_C5KY1ANXX/Sample_g_rhea_3kb_GTGAAA/*.fastq.gz | grep -A 3 '^@.* [^:]*:N:[^:]*:' | grep -v "^--$" >> rheAme_R2_3kb.fastq
zcat 140808_D00365_0310_AC5KY1ANXX/BclToFastq_Lane2_Indexlength6_Run1/Project_C5KY1ANXX/Sample_ls_kiwi_220_TAGTTG/*.fastq.gz | grep -A 3 '^@.* [^:]*:N:[^:]*:' | grep -v "^--$" >> aptOwe_R2_220.fastq
zcat 140808_D00365_0310_AC5KY1ANXX/BclToFastq_Lane2_Indexlength6_Run1/Project_C5KY1ANXX/Sample_ls_kiwi_3kb_CTTGTA/*.fastq.gz | grep -A 3 '^@.* [^:]*:N:[^:]*:' | grep -v "^--$" >> aptOwe_R2_3kb.fastq
zcat 140808_D00365_0310_AC5KY1ANXX/BclToFastq_Lane2_Indexlength6_Run1/Project_C5KY1ANXX/Sample_s_cassowary_220_GCTACC/*.fastq.gz | grep -A 3 '^@.* [^:]*:N:[^:]*:' | grep -v "^--$" >> casCas_R2_220.fastq
zcat 140808_D00365_0310_AC5KY1ANXX/BclToFastq_Lane2_Indexlength6_Run1/Project_C5KY1ANXX/Sample_s_cassowary_3kb_GCCAAT/*.fastq.gz | grep -A 3 '^@.* [^:]*:N:[^:]*:' | grep -v "^--$" >> casCas_R2_3kb.fastq

zcat 140808_D00365_0310_AC5KY1ANXX/BclToFastq_Lane3_Indexlength6_Run1/Project_C5KY1ANXX/Sample_emu_220_CTCTAC/*.fastq.gz | grep -A 3 '^@.* [^:]*:N:[^:]*:' | grep -v "^--$" >> droNov_R2_220.fastq
zcat 140808_D00365_0310_AC5KY1ANXX/BclToFastq_Lane3_Indexlength6_Run1/Project_C5KY1ANXX/Sample_emu_3kb_ACAGTG/*.fastq.gz | grep -A 3 '^@.* [^:]*:N:[^:]*:' | grep -v "^--$" >> droNov_R2_3kb.fastq
zcat 140808_D00365_0310_AC5KY1ANXX/BclToFastq_Lane3_Indexlength6_Run1/Project_C5KY1ANXX/Sample_g_rhea_220_ATCAGT/*.fastq.gz | grep -A 3 '^@.* [^:]*:N:[^:]*:' | grep -v "^--$" >> rheAme_R2_220.fastq
zcat 140808_D00365_0310_AC5KY1ANXX/BclToFastq_Lane3_Indexlength6_Run1/Project_C5KY1ANXX/Sample_g_rhea_3kb_GTGAAA/*.fastq.gz | grep -A 3 '^@.* [^:]*:N:[^:]*:' | grep -v "^--$" >> rheAme_R2_3kb.fastq
zcat 140808_D00365_0310_AC5KY1ANXX/BclToFastq_Lane3_Indexlength6_Run1/Project_C5KY1ANXX/Sample_ls_kiwi_220_TAGTTG/*.fastq.gz | grep -A 3 '^@.* [^:]*:N:[^:]*:' | grep -v "^--$" >> aptOwe_R2_220.fastq
zcat 140808_D00365_0310_AC5KY1ANXX/BclToFastq_Lane3_Indexlength6_Run1/Project_C5KY1ANXX/Sample_ls_kiwi_3kb_CTTGTA/*.fastq.gz | grep -A 3 '^@.* [^:]*:N:[^:]*:' | grep -v "^--$" >> aptOwe_R2_3kb.fastq
zcat 140808_D00365_0310_AC5KY1ANXX/BclToFastq_Lane3_Indexlength6_Run1/Project_C5KY1ANXX/Sample_s_cassowary_220_GCTACC/*.fastq.gz | grep -A 3 '^@.* [^:]*:N:[^:]*:' | grep -v "^--$" >> casCas_R2_220.fastq
zcat 140808_D00365_0310_AC5KY1ANXX/BclToFastq_Lane3_Indexlength6_Run1/Project_C5KY1ANXX/Sample_s_cassowary_3kb_GCCAAT/*.fastq.gz | grep -A 3 '^@.* [^:]*:N:[^:]*:' | grep -v "^--$" >> casCas_R2_3kb.fastq

zcat 140808_D00365_0310_AC5KY1ANXX/BclToFastq_Lane4_Indexlength6_Run1/Project_C5KY1ANXX/Sample_emu_220_CTCTAC/*.fastq.gz | grep -A 3 '^@.* [^:]*:N:[^:]*:' | grep -v "^--$" >> droNov_R2_220.fastq
zcat 140808_D00365_0310_AC5KY1ANXX/BclToFastq_Lane4_Indexlength6_Run1/Project_C5KY1ANXX/Sample_emu_3kb_ACAGTG/*.fastq.gz | grep -A 3 '^@.* [^:]*:N:[^:]*:' | grep -v "^--$" >> droNov_R2_3kb.fastq
zcat 140808_D00365_0310_AC5KY1ANXX/BclToFastq_Lane4_Indexlength6_Run1/Project_C5KY1ANXX/Sample_g_rhea_220_ATCAGT/*.fastq.gz | grep -A 3 '^@.* [^:]*:N:[^:]*:' | grep -v "^--$" >> rheAme_R2_220.fastq
zcat 140808_D00365_0310_AC5KY1ANXX/BclToFastq_Lane4_Indexlength6_Run1/Project_C5KY1ANXX/Sample_g_rhea_3kb_GTGAAA/*.fastq.gz | grep -A 3 '^@.* [^:]*:N:[^:]*:' | grep -v "^--$" >> rheAme_R2_3kb.fastq
zcat 140808_D00365_0310_AC5KY1ANXX/BclToFastq_Lane4_Indexlength6_Run1/Project_C5KY1ANXX/Sample_ls_kiwi_220_TAGTTG/*.fastq.gz | grep -A 3 '^@.* [^:]*:N:[^:]*:' | grep -v "^--$" >> aptOwe_R2_220.fastq
zcat 140808_D00365_0310_AC5KY1ANXX/BclToFastq_Lane4_Indexlength6_Run1/Project_C5KY1ANXX/Sample_ls_kiwi_3kb_CTTGTA/*.fastq.gz | grep -A 3 '^@.* [^:]*:N:[^:]*:' | grep -v "^--$" >> aptOwe_R2_3kb.fastq
zcat 140808_D00365_0310_AC5KY1ANXX/BclToFastq_Lane4_Indexlength6_Run1/Project_C5KY1ANXX/Sample_s_cassowary_220_GCTACC/*.fastq.gz | grep -A 3 '^@.* [^:]*:N:[^:]*:' | grep -v "^--$" >> casCas_R2_220.fastq
zcat 140808_D00365_0310_AC5KY1ANXX/BclToFastq_Lane4_Indexlength6_Run1/Project_C5KY1ANXX/Sample_s_cassowary_3kb_GCCAAT/*.fastq.gz | grep -A 3 '^@.* [^:]*:N:[^:]*:' | grep -v "^--$" >> casCas_R2_3kb.fastq

zcat 140808_D00365_0310_AC5KY1ANXX/BclToFastq_Lane5_Indexlength6_Run1/Project_C5KY1ANXX/Sample_emu_220_CTCTAC/*.fastq.gz | grep -A 3 '^@.* [^:]*:N:[^:]*:' | grep -v "^--$" >> droNov_R2_220.fastq
zcat 140808_D00365_0310_AC5KY1ANXX/BclToFastq_Lane5_Indexlength6_Run1/Project_C5KY1ANXX/Sample_emu_3kb_ACAGTG/*.fastq.gz | grep -A 3 '^@.* [^:]*:N:[^:]*:' | grep -v "^--$" >> droNov_R2_3kb.fastq
zcat 140808_D00365_0310_AC5KY1ANXX/BclToFastq_Lane5_Indexlength6_Run1/Project_C5KY1ANXX/Sample_g_rhea_220_ATCAGT/*.fastq.gz | grep -A 3 '^@.* [^:]*:N:[^:]*:' | grep -v "^--$" >> rheAme_R2_220.fastq
zcat 140808_D00365_0310_AC5KY1ANXX/BclToFastq_Lane5_Indexlength6_Run1/Project_C5KY1ANXX/Sample_g_rhea_3kb_GTGAAA/*.fastq.gz | grep -A 3 '^@.* [^:]*:N:[^:]*:' | grep -v "^--$" >> rheAme_R2_3kb.fastq
zcat 140808_D00365_0310_AC5KY1ANXX/BclToFastq_Lane5_Indexlength6_Run1/Project_C5KY1ANXX/Sample_ls_kiwi_220_TAGTTG/*.fastq.gz | grep -A 3 '^@.* [^:]*:N:[^:]*:' | grep -v "^--$" >> aptOwe_R2_220.fastq
zcat 140808_D00365_0310_AC5KY1ANXX/BclToFastq_Lane5_Indexlength6_Run1/Project_C5KY1ANXX/Sample_ls_kiwi_3kb_CTTGTA/*.fastq.gz | grep -A 3 '^@.* [^:]*:N:[^:]*:' | grep -v "^--$" >> aptOwe_R2_3kb.fastq
zcat 140808_D00365_0310_AC5KY1ANXX/BclToFastq_Lane5_Indexlength6_Run1/Project_C5KY1ANXX/Sample_s_cassowary_220_GCTACC/*.fastq.gz | grep -A 3 '^@.* [^:]*:N:[^:]*:' | grep -v "^--$" >> casCas_R2_220.fastq
zcat 140808_D00365_0310_AC5KY1ANXX/BclToFastq_Lane5_Indexlength6_Run1/Project_C5KY1ANXX/Sample_s_cassowary_3kb_GCCAAT/*.fastq.gz | grep -A 3 '^@.* [^:]*:N:[^:]*:' | grep -v "^--$" >> casCas_R2_3kb.fastq

zcat 140808_D00365_0310_AC5KY1ANXX/BclToFastq_Lane6_Indexlength6_Run1/Project_C5KY1ANXX/Sample_emu_220_CTCTAC/*.fastq.gz | grep -A 3 '^@.* [^:]*:N:[^:]*:' | grep -v "^--$" >> droNov_R2_220.fastq
zcat 140808_D00365_0310_AC5KY1ANXX/BclToFastq_Lane6_Indexlength6_Run1/Project_C5KY1ANXX/Sample_emu_3kb_ACAGTG/*.fastq.gz | grep -A 3 '^@.* [^:]*:N:[^:]*:' | grep -v "^--$" >> droNov_R2_3kb.fastq
zcat 140808_D00365_0310_AC5KY1ANXX/BclToFastq_Lane6_Indexlength6_Run1/Project_C5KY1ANXX/Sample_g_rhea_220_ATCAGT/*.fastq.gz | grep -A 3 '^@.* [^:]*:N:[^:]*:' | grep -v "^--$" >> rheAme_R2_220.fastq
zcat 140808_D00365_0310_AC5KY1ANXX/BclToFastq_Lane6_Indexlength6_Run1/Project_C5KY1ANXX/Sample_g_rhea_3kb_GTGAAA/*.fastq.gz | grep -A 3 '^@.* [^:]*:N:[^:]*:' | grep -v "^--$" >> rheAme_R2_3kb.fastq
zcat 140808_D00365_0310_AC5KY1ANXX/BclToFastq_Lane6_Indexlength6_Run1/Project_C5KY1ANXX/Sample_ls_kiwi_220_TAGTTG/*.fastq.gz | grep -A 3 '^@.* [^:]*:N:[^:]*:' | grep -v "^--$" >> aptOwe_R2_220.fastq
zcat 140808_D00365_0310_AC5KY1ANXX/BclToFastq_Lane6_Indexlength6_Run1/Project_C5KY1ANXX/Sample_ls_kiwi_3kb_CTTGTA/*.fastq.gz | grep -A 3 '^@.* [^:]*:N:[^:]*:' | grep -v "^--$" >> aptOwe_R2_3kb.fastq
zcat 140808_D00365_0310_AC5KY1ANXX/BclToFastq_Lane6_Indexlength6_Run1/Project_C5KY1ANXX/Sample_s_cassowary_220_GCTACC/*.fastq.gz | grep -A 3 '^@.* [^:]*:N:[^:]*:' | grep -v "^--$" >> casCas_R2_220.fastq
zcat 140808_D00365_0310_AC5KY1ANXX/BclToFastq_Lane6_Indexlength6_Run1/Project_C5KY1ANXX/Sample_s_cassowary_3kb_GCCAAT/*.fastq.gz | grep -A 3 '^@.* [^:]*:N:[^:]*:' | grep -v "^--$" >> casCas_R2_3kb.fastq

zcat 140808_D00365_0310_AC5KY1ANXX/BclToFastq_Lane7_Indexlength6_Run1/Project_C5KY1ANXX/Sample_emu_220_CTCTAC/*.fastq.gz | grep -A 3 '^@.* [^:]*:N:[^:]*:' | grep -v "^--$" >> droNov_R2_220.fastq
zcat 140808_D00365_0310_AC5KY1ANXX/BclToFastq_Lane7_Indexlength6_Run1/Project_C5KY1ANXX/Sample_emu_3kb_ACAGTG/*.fastq.gz | grep -A 3 '^@.* [^:]*:N:[^:]*:' | grep -v "^--$" >> droNov_R2_3kb.fastq
zcat 140808_D00365_0310_AC5KY1ANXX/BclToFastq_Lane7_Indexlength6_Run1/Project_C5KY1ANXX/Sample_g_rhea_220_ATCAGT/*.fastq.gz | grep -A 3 '^@.* [^:]*:N:[^:]*:' | grep -v "^--$" >> rheAme_R2_220.fastq
zcat 140808_D00365_0310_AC5KY1ANXX/BclToFastq_Lane7_Indexlength6_Run1/Project_C5KY1ANXX/Sample_g_rhea_3kb_GTGAAA/*.fastq.gz | grep -A 3 '^@.* [^:]*:N:[^:]*:' | grep -v "^--$" >> rheAme_R2_3kb.fastq
zcat 140808_D00365_0310_AC5KY1ANXX/BclToFastq_Lane7_Indexlength6_Run1/Project_C5KY1ANXX/Sample_ls_kiwi_220_TAGTTG/*.fastq.gz | grep -A 3 '^@.* [^:]*:N:[^:]*:' | grep -v "^--$" >> aptOwe_R2_220.fastq
zcat 140808_D00365_0310_AC5KY1ANXX/BclToFastq_Lane7_Indexlength6_Run1/Project_C5KY1ANXX/Sample_ls_kiwi_3kb_CTTGTA/*.fastq.gz | grep -A 3 '^@.* [^:]*:N:[^:]*:' | grep -v "^--$" >> aptOwe_R2_3kb.fastq
zcat 140808_D00365_0310_AC5KY1ANXX/BclToFastq_Lane7_Indexlength6_Run1/Project_C5KY1ANXX/Sample_s_cassowary_220_GCTACC/*.fastq.gz | grep -A 3 '^@.* [^:]*:N:[^:]*:' | grep -v "^--$" >> casCas_R2_220.fastq
zcat 140808_D00365_0310_AC5KY1ANXX/BclToFastq_Lane7_Indexlength6_Run1/Project_C5KY1ANXX/Sample_s_cassowary_3kb_GCCAAT/*.fastq.gz | grep -A 3 '^@.* [^:]*:N:[^:]*:' | grep -v "^--$" >> casCas_R2_3kb.fastq

zcat 140808_D00365_0310_AC5KY1ANXX/BclToFastq_Lane8_Indexlength6_Run1/Project_C5KY1ANXX/Sample_emu_220_CTCTAC/*.fastq.gz | grep -A 3 '^@.* [^:]*:N:[^:]*:' | grep -v "^--$" >> droNov_R2_220.fastq
zcat 140808_D00365_0310_AC5KY1ANXX/BclToFastq_Lane8_Indexlength6_Run1/Project_C5KY1ANXX/Sample_emu_3kb_ACAGTG/*.fastq.gz | grep -A 3 '^@.* [^:]*:N:[^:]*:' | grep -v "^--$" >> droNov_R2_3kb.fastq
zcat 140808_D00365_0310_AC5KY1ANXX/BclToFastq_Lane8_Indexlength6_Run1/Project_C5KY1ANXX/Sample_g_rhea_220_ATCAGT/*.fastq.gz | grep -A 3 '^@.* [^:]*:N:[^:]*:' | grep -v "^--$" >> rheAme_R2_220.fastq
zcat 140808_D00365_0310_AC5KY1ANXX/BclToFastq_Lane8_Indexlength6_Run1/Project_C5KY1ANXX/Sample_g_rhea_3kb_GTGAAA/*.fastq.gz | grep -A 3 '^@.* [^:]*:N:[^:]*:' | grep -v "^--$" >> rheAme_R2_3kb.fastq
zcat 140808_D00365_0310_AC5KY1ANXX/BclToFastq_Lane8_Indexlength6_Run1/Project_C5KY1ANXX/Sample_ls_kiwi_220_TAGTTG/*.fastq.gz | grep -A 3 '^@.* [^:]*:N:[^:]*:' | grep -v "^--$" >> aptOwe_R2_220.fastq
zcat 140808_D00365_0310_AC5KY1ANXX/BclToFastq_Lane8_Indexlength6_Run1/Project_C5KY1ANXX/Sample_ls_kiwi_3kb_CTTGTA/*.fastq.gz | grep -A 3 '^@.* [^:]*:N:[^:]*:' | grep -v "^--$" >> aptOwe_R2_3kb.fastq
zcat 140808_D00365_0310_AC5KY1ANXX/BclToFastq_Lane8_Indexlength6_Run1/Project_C5KY1ANXX/Sample_s_cassowary_220_GCTACC/*.fastq.gz | grep -A 3 '^@.* [^:]*:N:[^:]*:' | grep -v "^--$" >> casCas_R2_220.fastq
zcat 140808_D00365_0310_AC5KY1ANXX/BclToFastq_Lane8_Indexlength6_Run1/Project_C5KY1ANXX/Sample_s_cassowary_3kb_GCCAAT/*.fastq.gz | grep -A 3 '^@.* [^:]*:N:[^:]*:' | grep -v "^--$" >> casCas_R2_3kb.fastq
