##code for running the hmm search on protein set
##first run gff parsing code and oma parsing code in gff_parsing directory (download GFFs, point at HOGFasta)
##this code picks up from oma_parse.unassigned list

cut -f3,3 oma_parse.unassigned > unassigned_prot_list
perl -p -i -e 's/(\w+)\.\d+/$1/' unassigned_prot_list

#make input protein db
cat ~/ratite_store/annotation/maker_final/*.proteins.fasta > new_prots.fa
perl -p -i -e 's/(>\S+)\s*.*$/$1/' new_prots.fa
./get_ncbi_proteins.sh
cat ??????.fa > ncbi_prots.fa
rm ??????.fa
perl -p -i -e 's/>\w*\|\w*\|\w*\|(\w+)\.\d+\|.*$/>$1/' ncbi_prots.fa
cat ncbi_prots.fa new_prots.fa > all_input_prot.fa

#get unassigned
seqtk subseq all_input_prot.fa unassigned_prot_list > unassigned_prot_seq.fa

#verify
grep ">" unassigned_prot_seq.fa > unassigned_deflines
perl -p -i -e 's/>//' unassigned_deflines
diff <(sort unassigned_deflines) <(sort unassigned_prot_list) > missing_cases

#get assigned proteins in one big lump
cat ../OutputFolder/HOGFasta/*.fa > assigned_prot_seq.fa

#clean up deflines
cat *_prot_seqs.fa > seqs_for_hmm.fa
perl -p -i -e 's/>\w+\|\S*\|(\S+)\s+.*$/>$1/' seqs_for_hmm.fa
perl -p -i -e 's/(>\w+)\.\d+.*$/$1/' seqs_for_hmm.fa

#prep hmm database
cat ../hmms/*.hmm > hog_hmms

#should now have a directory with the input file, seqs_for_hmm.fa and the database to search, hog_hmms
#note that one test profile with hmmsearch takes about 15 seconds
#there are 27873 profiles to search against 687187 sequences
#split into 100 array jobs with ~278 profiles each
#each array gets 8 cpus
#use a command line this:
#hmmsearch -E <eval> --cpu 8 --tblout HOG<num>.tabout --noali ../hmms/<INPUT>.hmm seqs_for_hmm.fa 

#first prep file parts
ls hmms > hmm.list
perl -p -i -e 's/\.hmm//' hmm.list 
split -a 3 -d -l 250 hmm.list hmmpart. #make file parts

