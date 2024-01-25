#!/usr/bin/env python3

import sys
from Bio import SeqIO
import re
import os

## This script retrieves protein IDs from the header ID line of 
## a fasta file using regex. The fasta file can be nucleotide or protein
## sequences.

## usage: /path/to/script input.fasta

def get_prot_IDs(input_fasta):
	name, ext = os.path.splitext(input_fasta)
	output_file = f'{name}_protIDs.txt'
	
	with open(output_file, "w") as output:
		for record in SeqIO.parse(input_fasta, "fasta"):
			match = re.search(r'.\d_cds_(.+.\d_\d+).+$', record.id)
			prot_ID = match.group(1)
			output.write(f'{prot_ID}\n')	

	return output


def main():

	
	usage = f'\n\n\tusage: {sys.argv[0]} input.fasta.\n\n\tReturns an output <[input-filename]_protIDs.txt> file to current directory.\n\n'

	if len(sys.argv) < 2:
		sys.stderr.write(usage)
		sys.exit(1)
	
	input_fasta = sys.argv[1]
	get_prot_IDs(input_fasta)

if __name__ == '__main__':
	main()


