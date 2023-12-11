#!/usr/bin/env python3

import sys
import os
from Bio import Entrez
from Bio import SeqIO
from Bio import GenBank

## This script retrieves protein sequences from a list of accessions.
## Accession list should be provided as a <txt> file with line-
## separated protein accession version numbers.

## Usage: /path/to/script acc_list.txt

def entrez_get_prots_from_acc(acc_file):
	
	name, ext = os.path.splitext(acc_file)
	name = name.replace("_protIDs", "")
	output_file = f'{name}.faa'
	Entrez.email = "halvorsen1@llnl.gov"


	## Read in accessions	
	accs = list()
	with open(acc_file, "r") as acc_list:
		for line in acc_list:
			line = line.rstrip()
			if not line.startswith('#'):
				accs.append(line)

	## Create string of all accessions
	accs_str = ",".join(str(acc) for acc in accs)

	## Access sequences and write to output file
	with open(output_file, "w") as output:
		handle = Entrez.efetch(db='protein', id=accs_str, rettype = "fasta", retmode="text")
		for seq_record in SeqIO.parse(handle, "fasta"):
			output.write(f'>{seq_record.description}\n{seq_record.seq}\n')
		handle.close()
	
	print(f'Done writing protein sequences to {output_file}')	
	return output


def main():

	usage = f'\n\n\tusage: {sys.argv[0]} "acc_list.txt"\n\n'

	if len(sys.argv) < 2:
		sys.stderr.write(usage)
		sys.exit(1)

	acc_file = sys.argv[1]

	entrez_get_prots_from_acc(acc_file)

if __name__ == '__main__':
	main()
