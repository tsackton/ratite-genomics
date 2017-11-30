#!/bin/bash

SP=$1
GENE=$2
cd annotation/$SP
fasta_merge -d ${SP}.maker.output/${SP}_master_datastore_index.log -o ${SP}.genome
gff3_merge -d ${SP}.maker.output/${SP}_master_datastore_index.log -o ${SP}.genome.gff
maker_map_ids --prefix ${GENE} --justify 6 ${SP}.genome.gff > ${SP}.id.map
cp ${SP}.genome.gff ${SP}.genome.orig.gff
map_fasta_ids ${SP}.id.map ${SP}.genome.all.maker.proteins.fasta &
map_fasta_ids ${SP}.id.map ${SP}.genome.all.maker.transcripts.fasta &
map_gff_ids ${SP}.id.map ${SP}.genome.gff &
