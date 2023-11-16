#!/usr/bin/env python3

import sys
import re

def fasta_seq_select(input_file, file_origin):
	summary = {}
	whole_seq = ''
	start = int(sys.argv[2]) - 1
	end = int(sys.argv[3]) + 1
	with open(input_file, "r") as file:

		if file_origin == "atcc":
			for line in file:
				line = line.rstrip()

				## Retrieve information from header

				if line.startswith('>'):
					match = re.search(r'^>.+\s(assembly_id)="(.+)"\s(genome_id)="(.+)"\s(atcc_catalog_number)="(.+)"\s(species)="(.+)"\s(contig_number)="(.+)"\s(topology)="(.+)"', line)
				
					## Make a dictionary of information from FASTA header.
					## If a match group was assigned an odd number, make it a key. 
					## then, assign the next group to the value.

					for i in range(13):
						if i % 2 == 1:
							summary[match.group(i)] = match.group(i+1)
				else:

					## collect sequence into one string

					line = line.rstrip()
					whole_seq += line
					coordinates = range(start+1,end-1)

			## Select sequence based on user input
				
			selection = whole_seq[start:end]

			print(selection)
			
		elif file_origin == "ncbi":
			for line in file:
				line = line.rstrip()
				if line.startswith('>'):
					match = re.search(r'^>(\S+)\s(\S+\s\S+\s\S+),\s(\S+\s\S+)', line)
					
					summary['accession'] = match.group(1)
					summary['species'] = match.group(2)
					summary['seq description'] = match.group(3)

				else:
					line = line.rstrip()
					whole_seq += line
					coordinates = range(start+1,end-1)

			selection = whole_seq[start:end]
				
			print(f'{selection}')

		else:
			sys.stderr.write(usage)
			sys.exit(1)

	return summary, coordinates, selection
				

def main():
	

	usage = f'\n\n\t{sys.argv[0]} input_file.fa seq_start seq_end file_origin (ex: "ncbi" or "atcc")\n\n'
	
	if len(sys.argv) < 5:
		sys.stderr.write(usage)
		sys.exit(1)
	
	input_file = sys.argv[1]
	file_origin = sys.argv[4]
	fasta_seq_select(input_file, file_origin)

if __name__ == '__main__':
	main()

