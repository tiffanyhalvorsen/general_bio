#!/usr/bin/env python3

import sys
import re
import json
from Bio import SeqIO


## collect command line arguments as input

fastafile = sys.argv[1]
jsonfile = sys.argv[2]

## create dictionary from fasta records with sequences as keys
## interproscan uses the sequence as a unique identifier (not seq id)


## fasta input file source: ensembl biomart

fasta_dict = {}

with open(fastafile, "r") as file:
    for line in file:
        line = line.rstrip('\n*')
			
        # find headers
        if line.startswith('>'):
				
            # capture header information
            match = re.search(r'^>(\S+)$', line)
            seq_id = match.group(1)
            #desc = match.group(2)
		
            # add header information to dictionary
            fasta_dict[seq_id] = ""

	    # find sequences and add them to dictionary
        else:
            fasta_dict[seq_id] += line

# reverse dictionary to get sequences as keys
fasta_dict = dict((v, k) for k, v in fasta_dict.items())	
## create a dict for counting PANTHER2GO terms
## create a dict for counting InterPro2GO terms
## create a dict for counting PFAM domains

panther2go = {}
interpro2go = {}
pfam = {}

# parse json file into dictionary

with open(jsonfile, "r") as file:

    ## open file in its entirety to work with json.loads()
    data = file.read()

## load in file with json module
iprscan = json.loads(data)

## iprscan data structure = { 'results' : [ {} , {} ] } 
## extract GO term info (id and name) from PANTHER matches

for result in iprscan['results']:
    sequence = result['sequence']
    seq_id = fasta_dict[sequence]

    ## print a line to separate records
    print("==================================")
    print(f'{seq_id}\n{sequence}')

    #for goterm in oprscan['results']['matches']['goXRefs']:
     #   print(goterm)
