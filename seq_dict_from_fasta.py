#!/usr/bin/env python3

import sys
import re


	# This program takes a file as input, opens it for reading,
	# removes all line end characters	 
	# and uses a regex to search for headers, then extracts 
	# the seq_id, description, and sequence from the fasta 
	# records into a nested dictionary.

def seq_dict_from_fasta(input_file):
	input_file = sys.argv[1]
	
	usage = f'\n\n\tusage: {sys.argv[0]} filename.fasta'
	
	if len(sys.argv) < 2:
		sys.stderr.write(usage)
		sys.exit(1)

	fastaDict = {}
	with open(input_file, "r") as fasta_all:
		
		# remove end of line characters
		for line in fasta_all:
			line = line.rstrip()
			
			# find headers
			if line.startswith('>'):
				
				# capture header information
				match = re.search(r'^>(\S+)\s+(.+)$', line)
				seq_id = match.group(1)
				desc = match.group(2)
		
				# add header information to dictionary
				fastaDict[seq_id] = {'desc':desc,
													 'seq':""}

			# find sequences and add them to dictionary
			else:
				fastaDict[seq_id]['seq'] += line
		
	print(fastaDict)

	if __name__ = '__main__':
		main()
