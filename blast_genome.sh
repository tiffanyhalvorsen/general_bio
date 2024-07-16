#!/bin/bash

### This script takes a FASTA genome file and a query gene file as input and does a nucleotide
### blast search for the query. Must provide full path to both files. All output is written
### to current working directory. 

### Dependencies: prodigal, blast

## User-provided genome and query. 
## Extract names for making new files.

GENOME=$1
QUERY=$2
BASENAME=$(basename -- $GENOME)
GENOME_NAME="${BASENAME%.*}"
BASENAME2=$(basename -- $QUERY)
QUERY_NAME="${BASENAME2%.*}"

#if [[ ! $GENOME =~ ^*+.f[n]?[ast]?a$ ]]
#	echo "Invalid genome file format"
#	exit 1
OPT1="^.+\.fna$"
OPT2="^.+\.fa$"
OPT3="^.+\.fasta$"

until [[ $GENOME =~ ($OPT1|$OPT2|$OPT3) ]]; do
  printf >&2 '\nInvalid extension: %s\n' "$GENOME"
  exit 1
done
## Find ORFs and make gene and protein file with prodigal, if not done already

TEST1="${GENOME_NAME}_genes.fna"
TEST2="${GENOME_NAME}_db.fna.ndb"

if [ ! -s $TEST1 ]; then

	echo -e "\nNo ORF file. Processing genome with prodigal..."

	prodigal -i $GENOME -d ${GENOME_NAME}_genes.fna -a ${GENOME_NAME}_prots.faa \
		&> ${GENOME_NAME}_prodigal.log

	echo -e "\nDone."
fi


if [ ! -s $TEST2 ]; then
	
	echo -e "\nNo database yet. Processing genome with makeblastdb..."

	makeblastdb -in ${GENOME_NAME}_genes.fna -dbtype nucl \
		-out ${GENOME_NAME}_db.fna &> ${GENOME_NAME}_makeblastdb.log

	echo -e "\nDone."
fi

echo -e "\nRunning blastn with $QUERY..."

blastn -query $QUERY -db ${GENOME_NAME}_db.fna -outfmt 7 \
	-out ${GENOME_NAME}_${QUERY_NAME}_blastn.out

echo -e "\nDone.\n"

