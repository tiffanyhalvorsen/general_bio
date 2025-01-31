#!/bin/bash

## This script takes a <txt> file as input and uses the NCBI datasets tool
## to download the genome assembly FASTA file, protein sequence genbank file, and
## ORF file for each accession. Files are saved to a folder with the accession
## number in the current working directory.


## USAGE: download_genomes_from_accessions.sh <accessions_list.txt>

# Input accession file
accessions=$1

while read line; do
	datasets download genome accession ${line} --include protein,cds,genome --filename ${line}.zip
	unzip ${line}.zip
	rm *.zip
	rm ./ncbi_dataset/data/assembly_data_report.jsonl ./ncbi_dataset/data/dataset_catalog.json
	rm README.md
	mv ./ncbi_dataset/data/${line} ./
	rm -r ncbi_dataset/
done < ${accessions}

