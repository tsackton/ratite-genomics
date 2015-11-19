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