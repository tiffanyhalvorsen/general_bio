#!/usr/bin/env python3

import sys
import re
import os

## This program allows you to extract specific sequences from a whole genome fasta, using coordinates.
## usage:
## (1) input fasta file
## (2) sequence start
## (3) sequence end
## (4) origin of fasta (ncbi or ATCC) 

class IncompleteAssembly(Exception):
	pass

def rev_comp(dna):
    dna = dna.replace('A', 't')
    dna = dna.replace('T', 'a')
    dna = dna.replace('G', 'c')
    dna = dna.replace('C', 'g')
    dna = dna.upper()
    dna = dna[::-1]
    return dna

def fasta_seq_select(input_file, start, end, strand, file_origin):
	summary = {}
	whole_seq = ''
	start = int(start) - 1
	end = int(end) + 1
	strand = str(strand)

	try:
		file_extension = os.path.splitext(input_file)[1]
		options = ['.fa', '.fna', '.nt', '.fasta']
		if file_extension not in options:
			raise Exception("Filename does not end with either .fa, .fna, .nt, or .fasta")

	except Exception as error:
		print("Error: " + str(error))
		sys.exit(1)	

	try:
		strand = strand
		strands = ['reverse', 'Reverse', 'forward', 'Forward']
		if strand not in strands:
			raise Exception("Strand must be designated. Include 'forward' or 'reverse' as 4th parameter")
	
	except Exception as error:
		print("Error: " + str(error))
		sys.exit(1)	

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
				
				## collect sequence into one string
				else:
					line = line.rstrip()
					whole_seq += line
					coordinates = range(start+1,end-1)
					selection = whole_seq[start:end]
					
					## Select sequence based on user input
					if strand == "forward" or strand == "Forward":
						selection = selection
					elif strand == "reverse" or strand == "Reverse":
						selection = rev_comp(selection)			
		
			print(f'{selection}')
			
		elif file_origin == "ncbi":
			contig_count = 0
			for line in file:
				line = line.rstrip()

				if line.startswith('>'):
					contig_count += 1
					match = re.search(r'^>(\S+)\s(.+),\s(.+)$', line)
					
					summary['accession'] = match.group(1)
					summary['species'] = match.group(2)
					summary['seq description'] = match.group(3)
					
				else:
					line = line.rstrip()
					whole_seq += line
					coordinates = range(start+1,end-1)
					selection = whole_seq[start:end]
					
					## Select sequence based on user input
					if strand == "forward" or strand == "Forward":
						selection = selection
					elif strand == "reverse" or strand == "Reverse":
						line = line.rstrip()
						whole_seq += line

			
			if contig_count > 1: 
				raise IncompleteAssembly(f'\n\nMust provide assembled genome to find loci - this FASTA file has {contig_count} contigs.\n\n')	
			else:
				print(f'{selection}')

		else:
			sys.stderr.write(usage)
			sys.exit(1)

	return summary, coordinates, selection
				

def main():
	

	usage = f'\n\n\t{sys.argv[0]} [input_file.fa] [seq_start] [seq_end] [strand] (forward or reverse) [file_origin] (ex: "ncbi" or "atcc")\n\n'
	
	if len(sys.argv) < 5:
		sys.stderr.write(usage)
		sys.exit(1)
	
	input_file = sys.argv[1]
	start = sys.argv[2]
	end = sys.argv[3]
	strand = sys.argv[4]
	file_origin = sys.argv[5]
	fasta_seq_select(input_file, start, end, strand, file_origin)

if __name__ == '__main__':
	main()

