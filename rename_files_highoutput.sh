#clean up files

zcat 140808_D00365_0310_AC5KY1ANXX/BclToFastq_Lane1_Indexlength6_Run1/Project_C5KY1ANXX/Sample_emu_220_CTCTAC.R1.fastq.gz | grep -A 3 '^@.* [^:]*:N:[^:]*:' | grep -v "^--$" | gzip -c > droNov_R1_220.fastq.gz
zcat 140808_D00365_0310_AC5KY1ANXX/BclToFastq_Lane1_Indexlength6_Run1/Project_C5KY1ANXX/Sample_emu_3kb_ACAGTG.R1.fastq.gz | grep -A 3 '^@.* [^:]*:N:[^:]*:' | grep -v "^--$" | gzip -c > droNov_R1_3kb.fastq.gz
zcat 140808_D00365_0310_AC5KY1ANXX/BclToFastq_Lane1_Indexlength6_Run1/Project_C5KY1ANXX/Sample_g_rhea_220_ATCAGT.R1.fastq.gz | grep -A 3 '^@.* [^:]*:N:[^:]*:' | grep -v "^--$" | gzip -c > rheAme_R1_220.fastq.gz
zcat 140808_D00365_0310_AC5KY1ANXX/BclToFastq_Lane1_Indexlength6_Run1/Project_C5KY1ANXX/Sample_g_rhea_3kb_GTGAAA.R1.fastq.gz | grep -A 3 '^@.* [^:]*:N:[^:]*:' | grep -v "^--$" | gzip -c > rheAme_R1_3kb.fastq.gz
zcat 140808_D00365_0310_AC5KY1ANXX/BclToFastq_Lane1_Indexlength6_Run1/Project_C5KY1ANXX/Sample_ls_kiwi_220_TAGTTG.R1.fastq.gz | grep -A 3 '^@.* [^:]*:N:[^:]*:' | grep -v "^--$" | gzip -c > aptOwe_R1_220.fastq.gz
zcat 140808_D00365_0310_AC5KY1ANXX/BclToFastq_Lane1_Indexlength6_Run1/Project_C5KY1ANXX/Sample_ls_kiwi_3kb_CTTGTA.R1.fastq.gz | grep -A 3 '^@.* [^:]*:N:[^:]*:' | grep -v "^--$" | gzip -c > aptOwe_R1_3kb.fastq.gz
zcat 140808_D00365_0310_AC5KY1ANXX/BclToFastq_Lane1_Indexlength6_Run1/Project_C5KY1ANXX/Sample_s_cassowary_220_GCTACC.R1.fastq.gz | grep -A 3 '^@.* [^:]*:N:[^:]*:' | grep -v "^--$" | gzip -c > casCas_R1_220.fastq.gz
zcat 140808_D00365_0310_AC5KY1ANXX/BclToFastq_Lane1_Indexlength6_Run1/Project_C5KY1ANXX/Sample_s_cassowary_3kb_GCCAAT.R1.fastq.gz | grep -A 3 '^@.* [^:]*:N:[^:]*:' | grep -v "^--$" | gzip -c > casCas_R1_3kb.fastq.gz

zcat 140808_D00365_0310_AC5KY1ANXX/BclToFastq_Lane2_Indexlength6_Run1/Project_C5KY1ANXX/Sample_emu_220_CTCTAC.R1.fastq.gz | grep -A 3 '^@.* [^:]*:N:[^:]*:' | grep -v "^--$" | gzip -c >> droNov_R1_220.fastq.gz
zcat 140808_D00365_0310_AC5KY1ANXX/BclToFastq_Lane2_Indexlength6_Run1/Project_C5KY1ANXX/Sample_emu_3kb_ACAGTG.R1.fastq.gz | grep -A 3 '^@.* [^:]*:N:[^:]*:' | grep -v "^--$" | gzip -c >> droNov_R1_3kb.fastq.gz
zcat 140808_D00365_0310_AC5KY1ANXX/BclToFastq_Lane2_Indexlength6_Run1/Project_C5KY1ANXX/Sample_g_rhea_220_ATCAGT.R1.fastq.gz | grep -A 3 '^@.* [^:]*:N:[^:]*:' | grep -v "^--$" | gzip -c >> rheAme_R1_220.fastq.gz
zcat 140808_D00365_0310_AC5KY1ANXX/BclToFastq_Lane2_Indexlength6_Run1/Project_C5KY1ANXX/Sample_g_rhea_3kb_GTGAAA.R1.fastq.gz | grep -A 3 '^@.* [^:]*:N:[^:]*:' | grep -v "^--$" | gzip -c >> rheAme_R1_3kb.fastq.gz
zcat 140808_D00365_0310_AC5KY1ANXX/BclToFastq_Lane2_Indexlength6_Run1/Project_C5KY1ANXX/Sample_ls_kiwi_220_TAGTTG.R1.fastq.gz | grep -A 3 '^@.* [^:]*:N:[^:]*:' | grep -v "^--$" | gzip -c >> aptOwe_R1_220.fastq.gz
zcat 140808_D00365_0310_AC5KY1ANXX/BclToFastq_Lane2_Indexlength6_Run1/Project_C5KY1ANXX/Sample_ls_kiwi_3kb_CTTGTA.R1.fastq.gz | grep -A 3 '^@.* [^:]*:N:[^:]*:' | grep -v "^--$" | gzip -c >> aptOwe_R1_3kb.fastq.gz
zcat 140808_D00365_0310_AC5KY1ANXX/BclToFastq_Lane2_Indexlength6_Run1/Project_C5KY1ANXX/Sample_s_cassowary_220_GCTACC.R1.fastq.gz | grep -A 3 '^@.* [^:]*:N:[^:]*:' | grep -v "^--$" | gzip -c >> casCas_R1_220.fastq.gz
zcat 140808_D00365_0310_AC5KY1ANXX/BclToFastq_Lane2_Indexlength6_Run1/Project_C5KY1ANXX/Sample_s_cassowary_3kb_GCCAAT.R1.fastq.gz | grep -A 3 '^@.* [^:]*:N:[^:]*:' | grep -v "^--$" | gzip -c >> casCas_R1_3kb.fastq.gz

zcat 140808_D00365_0310_AC5KY1ANXX/BclToFastq_Lane3_Indexlength6_Run1/Project_C5KY1ANXX/Sample_emu_220_CTCTAC.R1.fastq.gz | grep -A 3 '^@.* [^:]*:N:[^:]*:' | grep -v "^--$" | gzip -c >> droNov_R1_220.fastq.gz
zcat 140808_D00365_0310_AC5KY1ANXX/BclToFastq_Lane3_Indexlength6_Run1/Project_C5KY1ANXX/Sample_emu_3kb_ACAGTG.R1.fastq.gz | grep -A 3 '^@.* [^:]*:N:[^:]*:' | grep -v "^--$" | gzip -c >> droNov_R1_3kb.fastq.gz
zcat 140808_D00365_0310_AC5KY1ANXX/BclToFastq_Lane3_Indexlength6_Run1/Project_C5KY1ANXX/Sample_g_rhea_220_ATCAGT.R1.fastq.gz | grep -A 3 '^@.* [^:]*:N:[^:]*:' | grep -v "^--$" | gzip -c >> rheAme_R1_220.fastq.gz
zcat 140808_D00365_0310_AC5KY1ANXX/BclToFastq_Lane3_Indexlength6_Run1/Project_C5KY1ANXX/Sample_g_rhea_3kb_GTGAAA.R1.fastq.gz | grep -A 3 '^@.* [^:]*:N:[^:]*:' | grep -v "^--$" | gzip -c >> rheAme_R1_3kb.fastq.gz
zcat 140808_D00365_0310_AC5KY1ANXX/BclToFastq_Lane3_Indexlength6_Run1/Project_C5KY1ANXX/Sample_ls_kiwi_220_TAGTTG.R1.fastq.gz | grep -A 3 '^@.* [^:]*:N:[^:]*:' | grep -v "^--$" | gzip -c >> aptOwe_R1_220.fastq.gz
zcat 140808_D00365_0310_AC5KY1ANXX/BclToFastq_Lane3_Indexlength6_Run1/Project_C5KY1ANXX/Sample_ls_kiwi_3kb_CTTGTA.R1.fastq.gz | grep -A 3 '^@.* [^:]*:N:[^:]*:' | grep -v "^--$" | gzip -c >> aptOwe_R1_3kb.fastq.gz
zcat 140808_D00365_0310_AC5KY1ANXX/BclToFastq_Lane3_Indexlength6_Run1/Project_C5KY1ANXX/Sample_s_cassowary_220_GCTACC.R1.fastq.gz | grep -A 3 '^@.* [^:]*:N:[^:]*:' | grep -v "^--$" | gzip -c >> casCas_R1_220.fastq.gz
zcat 140808_D00365_0310_AC5KY1ANXX/BclToFastq_Lane3_Indexlength6_Run1/Project_C5KY1ANXX/Sample_s_cassowary_3kb_GCCAAT.R1.fastq.gz | grep -A 3 '^@.* [^:]*:N:[^:]*:' | grep -v "^--$" | gzip -c >> casCas_R1_3kb.fastq.gz

zcat 140808_D00365_0310_AC5KY1ANXX/BclToFastq_Lane4_Indexlength6_Run1/Project_C5KY1ANXX/Sample_emu_220_CTCTAC.R1.fastq.gz | grep -A 3 '^@.* [^:]*:N:[^:]*:' | grep -v "^--$" | gzip -c >> droNov_R1_220.fastq.gz
zcat 140808_D00365_0310_AC5KY1ANXX/BclToFastq_Lane4_Indexlength6_Run1/Project_C5KY1ANXX/Sample_emu_3kb_ACAGTG.R1.fastq.gz | grep -A 3 '^@.* [^:]*:N:[^:]*:' | grep -v "^--$" | gzip -c >> droNov_R1_3kb.fastq.gz
zcat 140808_D00365_0310_AC5KY1ANXX/BclToFastq_Lane4_Indexlength6_Run1/Project_C5KY1ANXX/Sample_g_rhea_220_ATCAGT.R1.fastq.gz | grep -A 3 '^@.* [^:]*:N:[^:]*:' | grep -v "^--$" | gzip -c >> rheAme_R1_220.fastq.gz
zcat 140808_D00365_0310_AC5KY1ANXX/BclToFastq_Lane4_Indexlength6_Run1/Project_C5KY1ANXX/Sample_g_rhea_3kb_GTGAAA.R1.fastq.gz | grep -A 3 '^@.* [^:]*:N:[^:]*:' | grep -v "^--$" | gzip -c >> rheAme_R1_3kb.fastq.gz
zcat 140808_D00365_0310_AC5KY1ANXX/BclToFastq_Lane4_Indexlength6_Run1/Project_C5KY1ANXX/Sample_ls_kiwi_220_TAGTTG.R1.fastq.gz | grep -A 3 '^@.* [^:]*:N:[^:]*:' | grep -v "^--$" | gzip -c >> aptOwe_R1_220.fastq.gz
zcat 140808_D00365_0310_AC5KY1ANXX/BclToFastq_Lane4_Indexlength6_Run1/Project_C5KY1ANXX/Sample_ls_kiwi_3kb_CTTGTA.R1.fastq.gz | grep -A 3 '^@.* [^:]*:N:[^:]*:' | grep -v "^--$" | gzip -c >> aptOwe_R1_3kb.fastq.gz
zcat 140808_D00365_0310_AC5KY1ANXX/BclToFastq_Lane4_Indexlength6_Run1/Project_C5KY1ANXX/Sample_s_cassowary_220_GCTACC.R1.fastq.gz | grep -A 3 '^@.* [^:]*:N:[^:]*:' | grep -v "^--$" | gzip -c >> casCas_R1_220.fastq.gz
zcat 140808_D00365_0310_AC5KY1ANXX/BclToFastq_Lane4_Indexlength6_Run1/Project_C5KY1ANXX/Sample_s_cassowary_3kb_GCCAAT.R1.fastq.gz | grep -A 3 '^@.* [^:]*:N:[^:]*:' | grep -v "^--$" | gzip -c >> casCas_R1_3kb.fastq.gz

zcat 140808_D00365_0310_AC5KY1ANXX/BclToFastq_Lane5_Indexlength6_Run1/Project_C5KY1ANXX/Sample_emu_220_CTCTAC.R1.fastq.gz | grep -A 3 '^@.* [^:]*:N:[^:]*:' | grep -v "^--$" | gzip -c >> droNov_R1_220.fastq.gz
zcat 140808_D00365_0310_AC5KY1ANXX/BclToFastq_Lane5_Indexlength6_Run1/Project_C5KY1ANXX/Sample_emu_3kb_ACAGTG.R1.fastq.gz | grep -A 3 '^@.* [^:]*:N:[^:]*:' | grep -v "^--$" | gzip -c >> droNov_R1_3kb.fastq.gz
zcat 140808_D00365_0310_AC5KY1ANXX/BclToFastq_Lane5_Indexlength6_Run1/Project_C5KY1ANXX/Sample_g_rhea_220_ATCAGT.R1.fastq.gz | grep -A 3 '^@.* [^:]*:N:[^:]*:' | grep -v "^--$" | gzip -c >> rheAme_R1_220.fastq.gz
zcat 140808_D00365_0310_AC5KY1ANXX/BclToFastq_Lane5_Indexlength6_Run1/Project_C5KY1ANXX/Sample_g_rhea_3kb_GTGAAA.R1.fastq.gz | grep -A 3 '^@.* [^:]*:N:[^:]*:' | grep -v "^--$" | gzip -c >> rheAme_R1_3kb.fastq.gz
zcat 140808_D00365_0310_AC5KY1ANXX/BclToFastq_Lane5_Indexlength6_Run1/Project_C5KY1ANXX/Sample_ls_kiwi_220_TAGTTG.R1.fastq.gz | grep -A 3 '^@.* [^:]*:N:[^:]*:' | grep -v "^--$" | gzip -c >> aptOwe_R1_220.fastq.gz
zcat 140808_D00365_0310_AC5KY1ANXX/BclToFastq_Lane5_Indexlength6_Run1/Project_C5KY1ANXX/Sample_ls_kiwi_3kb_CTTGTA.R1.fastq.gz | grep -A 3 '^@.* [^:]*:N:[^:]*:' | grep -v "^--$" | gzip -c >> aptOwe_R1_3kb.fastq.gz
zcat 140808_D00365_0310_AC5KY1ANXX/BclToFastq_Lane5_Indexlength6_Run1/Project_C5KY1ANXX/Sample_s_cassowary_220_GCTACC.R1.fastq.gz | grep -A 3 '^@.* [^:]*:N:[^:]*:' | grep -v "^--$" | gzip -c >> casCas_R1_220.fastq.gz
zcat 140808_D00365_0310_AC5KY1ANXX/BclToFastq_Lane5_Indexlength6_Run1/Project_C5KY1ANXX/Sample_s_cassowary_3kb_GCCAAT.R1.fastq.gz | grep -A 3 '^@.* [^:]*:N:[^:]*:' | grep -v "^--$" | gzip -c >> casCas_R1_3kb.fastq.gz

zcat 140808_D00365_0310_AC5KY1ANXX/BclToFastq_Lane6_Indexlength6_Run1/Project_C5KY1ANXX/Sample_emu_220_CTCTAC.R1.fastq.gz | grep -A 3 '^@.* [^:]*:N:[^:]*:' | grep -v "^--$" | gzip -c >> droNov_R1_220.fastq.gz
zcat 140808_D00365_0310_AC5KY1ANXX/BclToFastq_Lane6_Indexlength6_Run1/Project_C5KY1ANXX/Sample_emu_3kb_ACAGTG.R1.fastq.gz | grep -A 3 '^@.* [^:]*:N:[^:]*:' | grep -v "^--$" | gzip -c >> droNov_R1_3kb.fastq.gz
zcat 140808_D00365_0310_AC5KY1ANXX/BclToFastq_Lane6_Indexlength6_Run1/Project_C5KY1ANXX/Sample_g_rhea_220_ATCAGT.R1.fastq.gz | grep -A 3 '^@.* [^:]*:N:[^:]*:' | grep -v "^--$" | gzip -c >> rheAme_R1_220.fastq.gz
zcat 140808_D00365_0310_AC5KY1ANXX/BclToFastq_Lane6_Indexlength6_Run1/Project_C5KY1ANXX/Sample_g_rhea_3kb_GTGAAA.R1.fastq.gz | grep -A 3 '^@.* [^:]*:N:[^:]*:' | grep -v "^--$" | gzip -c >> rheAme_R1_3kb.fastq.gz
zcat 140808_D00365_0310_AC5KY1ANXX/BclToFastq_Lane6_Indexlength6_Run1/Project_C5KY1ANXX/Sample_ls_kiwi_220_TAGTTG.R1.fastq.gz | grep -A 3 '^@.* [^:]*:N:[^:]*:' | grep -v "^--$" | gzip -c >> aptOwe_R1_220.fastq.gz
zcat 140808_D00365_0310_AC5KY1ANXX/BclToFastq_Lane6_Indexlength6_Run1/Project_C5KY1ANXX/Sample_ls_kiwi_3kb_CTTGTA.R1.fastq.gz | grep -A 3 '^@.* [^:]*:N:[^:]*:' | grep -v "^--$" | gzip -c >> aptOwe_R1_3kb.fastq.gz
zcat 140808_D00365_0310_AC5KY1ANXX/BclToFastq_Lane6_Indexlength6_Run1/Project_C5KY1ANXX/Sample_s_cassowary_220_GCTACC.R1.fastq.gz | grep -A 3 '^@.* [^:]*:N:[^:]*:' | grep -v "^--$" | gzip -c >> casCas_R1_220.fastq.gz
zcat 140808_D00365_0310_AC5KY1ANXX/BclToFastq_Lane6_Indexlength6_Run1/Project_C5KY1ANXX/Sample_s_cassowary_3kb_GCCAAT.R1.fastq.gz | grep -A 3 '^@.* [^:]*:N:[^:]*:' | grep -v "^--$" | gzip -c >> casCas_R1_3kb.fastq.gz

zcat 140808_D00365_0310_AC5KY1ANXX/BclToFastq_Lane7_Indexlength6_Run1/Project_C5KY1ANXX/Sample_emu_220_CTCTAC.R1.fastq.gz | grep -A 3 '^@.* [^:]*:N:[^:]*:' | grep -v "^--$" | gzip -c >> droNov_R1_220.fastq.gz
zcat 140808_D00365_0310_AC5KY1ANXX/BclToFastq_Lane7_Indexlength6_Run1/Project_C5KY1ANXX/Sample_emu_3kb_ACAGTG.R1.fastq.gz | grep -A 3 '^@.* [^:]*:N:[^:]*:' | grep -v "^--$" | gzip -c >> droNov_R1_3kb.fastq.gz
zcat 140808_D00365_0310_AC5KY1ANXX/BclToFastq_Lane7_Indexlength6_Run1/Project_C5KY1ANXX/Sample_g_rhea_220_ATCAGT.R1.fastq.gz | grep -A 3 '^@.* [^:]*:N:[^:]*:' | grep -v "^--$" | gzip -c >> rheAme_R1_220.fastq.gz
zcat 140808_D00365_0310_AC5KY1ANXX/BclToFastq_Lane7_Indexlength6_Run1/Project_C5KY1ANXX/Sample_g_rhea_3kb_GTGAAA.R1.fastq.gz | grep -A 3 '^@.* [^:]*:N:[^:]*:' | grep -v "^--$" | gzip -c >> rheAme_R1_3kb.fastq.gz
zcat 140808_D00365_0310_AC5KY1ANXX/BclToFastq_Lane7_Indexlength6_Run1/Project_C5KY1ANXX/Sample_ls_kiwi_220_TAGTTG.R1.fastq.gz | grep -A 3 '^@.* [^:]*:N:[^:]*:' | grep -v "^--$" | gzip -c >> aptOwe_R1_220.fastq.gz
zcat 140808_D00365_0310_AC5KY1ANXX/BclToFastq_Lane7_Indexlength6_Run1/Project_C5KY1ANXX/Sample_ls_kiwi_3kb_CTTGTA.R1.fastq.gz | grep -A 3 '^@.* [^:]*:N:[^:]*:' | grep -v "^--$" | gzip -c >> aptOwe_R1_3kb.fastq.gz
zcat 140808_D00365_0310_AC5KY1ANXX/BclToFastq_Lane7_Indexlength6_Run1/Project_C5KY1ANXX/Sample_s_cassowary_220_GCTACC.R1.fastq.gz | grep -A 3 '^@.* [^:]*:N:[^:]*:' | grep -v "^--$" | gzip -c >> casCas_R1_220.fastq.gz
zcat 140808_D00365_0310_AC5KY1ANXX/BclToFastq_Lane7_Indexlength6_Run1/Project_C5KY1ANXX/Sample_s_cassowary_3kb_GCCAAT.R1.fastq.gz | grep -A 3 '^@.* [^:]*:N:[^:]*:' | grep -v "^--$" | gzip -c >> casCas_R1_3kb.fastq.gz

zcat 140808_D00365_0310_AC5KY1ANXX/BclToFastq_Lane8_Indexlength6_Run1/Project_C5KY1ANXX/Sample_emu_220_CTCTAC.R1.fastq.gz | grep -A 3 '^@.* [^:]*:N:[^:]*:' | grep -v "^--$" | gzip -c >> droNov_R1_220.fastq.gz
zcat 140808_D00365_0310_AC5KY1ANXX/BclToFastq_Lane8_Indexlength6_Run1/Project_C5KY1ANXX/Sample_emu_3kb_ACAGTG.R1.fastq.gz | grep -A 3 '^@.* [^:]*:N:[^:]*:' | grep -v "^--$" | gzip -c >> droNov_R1_3kb.fastq.gz
zcat 140808_D00365_0310_AC5KY1ANXX/BclToFastq_Lane8_Indexlength6_Run1/Project_C5KY1ANXX/Sample_g_rhea_220_ATCAGT.R1.fastq.gz | grep -A 3 '^@.* [^:]*:N:[^:]*:' | grep -v "^--$" | gzip -c >> rheAme_R1_220.fastq.gz
zcat 140808_D00365_0310_AC5KY1ANXX/BclToFastq_Lane8_Indexlength6_Run1/Project_C5KY1ANXX/Sample_g_rhea_3kb_GTGAAA.R1.fastq.gz | grep -A 3 '^@.* [^:]*:N:[^:]*:' | grep -v "^--$" | gzip -c >> rheAme_R1_3kb.fastq.gz
zcat 140808_D00365_0310_AC5KY1ANXX/BclToFastq_Lane8_Indexlength6_Run1/Project_C5KY1ANXX/Sample_ls_kiwi_220_TAGTTG.R1.fastq.gz | grep -A 3 '^@.* [^:]*:N:[^:]*:' | grep -v "^--$" | gzip -c >> aptOwe_R1_220.fastq.gz
zcat 140808_D00365_0310_AC5KY1ANXX/BclToFastq_Lane8_Indexlength6_Run1/Project_C5KY1ANXX/Sample_ls_kiwi_3kb_CTTGTA.R1.fastq.gz | grep -A 3 '^@.* [^:]*:N:[^:]*:' | grep -v "^--$" | gzip -c >> aptOwe_R1_3kb.fastq.gz
zcat 140808_D00365_0310_AC5KY1ANXX/BclToFastq_Lane8_Indexlength6_Run1/Project_C5KY1ANXX/Sample_s_cassowary_220_GCTACC.R1.fastq.gz | grep -A 3 '^@.* [^:]*:N:[^:]*:' | grep -v "^--$" | gzip -c >> casCas_R1_220.fastq.gz
zcat 140808_D00365_0310_AC5KY1ANXX/BclToFastq_Lane8_Indexlength6_Run1/Project_C5KY1ANXX/Sample_s_cassowary_3kb_GCCAAT.R1.fastq.gz | grep -A 3 '^@.* [^:]*:N:[^:]*:' | grep -v "^--$" | gzip -c >> casCas_R1_3kb.fastq.gz
