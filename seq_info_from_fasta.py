#!/usr/bin/env python3

import sys
import math
from Bio import SeqIO

	# This program parses a fasta file using biopython and returns
    # 1) total number of records, 2) avg sequence length, 3) smallest
    # and largest sequences in file, 4) average, highest and lowest
    # GC content across sequences.

def seq_info_from_fasta(input_file, output_file):
    

    ## set empty values

    seq_records = {}
    record_num = 0
    total_nt = 0
    seq_lengths = []
    GC_contents = []
	
    ## read in file with biopython

    for seq_record in SeqIO.parse(str(input_file), "fasta"):
        record_num += 1

        ## calculate total nt content

        total_nt += len(seq_record.seq)

        ## set variables

        seq_id = seq_record.id
        seq = seq_record.seq

        ## calculate sequence length for each record, add to list

        seq_length = len(seq_record.seq)
        seq_lengths.append(int(seq_length))

        ## calculate GC content

        GC = (seq.count('G') + seq.count('C'))/len(seq)

        ## create entries in dictionary for each record with above information

        seq_records[seq_id] = [seq,seq_length,GC]

        ## use GC_contents list for calculating global GC content values

        GC_contents.append(int(GC*100))

    avg_record_len = total_nt/record_num

    with open(output_file, 'w') as file:
        file.write(f'Total number of sequences: {record_num:.2f}\nTotal number of nucleotides in file: {total_nt}.\nAverage sequence length: {avg_record_len:.2f}.\nSmallest sequence: {min(seq_lengths)} nts\nLongest sequence: {max(seq_lengths)} nts.\nAverage GC content of records in file: {sum(GC_contents)/len(GC_contents):.2f}%.\nLowest GC content in a single sequence: {min(GC_contents)}%.\nHighest GC content in a single sequence: {max(GC_contents)}%')

def main():

    ## gather commandline args

    input_file = sys.argv[1]
    output_file = sys.argv[2]
    
    usage = f'\n\n\tusage: {sys.argv[0]} input_filename.fasta output_filename.txt'

    if len(sys.argv) < 2:
        sys.stderr.write(usage)
        sys.exit(1) 
    
    seq_info_from_fasta(input_file, output_file)

if __name__ == '__main__':
    main()
