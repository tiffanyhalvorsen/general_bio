#!/usr/bin/env python3

import sys
import re

fasta = sys.argv[1]

contigs = 0
nt_comp = {}

contig_lengths = {}

nt_bycontig  = {}
			
with open(fasta, "r") as file:
	for line in file:
		line = line.rstrip()
		if line.startswith('>'):
			match = re.search(r'^>(\S+)\s\D{4}(\d+)\s(.+)$', line)
			contigs += 1
			contig_lengths[match.group(1)] = match.group(2)
			contig_name=match.group(1)
			nt_bycontig[contig_name] = ''

		else:
			
			for nt in line:
				if nt in nt_comp:
					nt_comp[nt] += 1
				else:
					nt_comp[nt] = 1
		
			for nt in line:
				print(nt_bycontig[contig_name])
				if nt in nt_bycontig[contig_name]:
					print(nt_bycontig[contig_name])
					nt_bycontig[contig_name][nt] += 1
				
				else:
					print('in the else')
					nt_bycontig[contig_name][nt] = 1

#print(nt_bycontig)
#print(nt_comp)

print(contig_lengths)
lengths = []
for key, value in contig_lengths.items():
	lengths.append(int(value))

print(f'smallest contig is: {min(lengths)}\tlongest contig is: {max(lengths)}')
