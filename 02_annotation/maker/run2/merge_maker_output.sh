#!/bin/bash

SP=$1
cd annotation/$SP
fasta_merge -d ${SP}.maker.output/${SP}_master_datastore_index.log -o ${SP}.genome
gff3_merge -d ${SP}.maker.output/${SP}_master_datastore_index.log -o ${SP}.genome.gff
