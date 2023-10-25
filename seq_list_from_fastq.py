#!/usr/bin/env python

import os, sys, re
from Bio import SeqIO



## method: seq_list_from_fastq_file(fastq_filename)
##
##  Extracts the sequence lines from a fastq file and returns a list
##  of the sequence lines
##
##  input parameters:
##
##  fastq_filename :  name of the fastq file (type: string)
##
##  returns seq_list : list of read sequences.
##                    ie.  ["GATCGCATAG", "CGATGCAG", ...]
    
def seq_list_from_fastq(fastq_filename):

    seq_list = list()
    for seq_record in SeqIO.parse(fastq_filename, "fastq"):
        seq_list.append(seq_record.seq)
    
    #with open(fastq_filename, "r") as file:
    ## begin your code
        #for line in file:
            #for match in re.findall(r"^@\.+\n(\.+)\n[+]\n(\w+)", line):
            #    sequence = match.group(1)
            #    QS = match.group[2]
            #    seq_list.append(sequence)

    return seq_list
    

def main():

    progname = sys.argv[0]
    
    usage = "\n\n\tusage: {} filename.fastq num_seqs_show\n\n\n".format(progname)
    
    if len(sys.argv) < 3:
        sys.stderr.write(usage)
        sys.exit(1)

    # capture command-line arguments
    fastq_filename = sys.argv[1]
    num_seqs_show = int(sys.argv[2])

    seq_list = seq_list_from_fastq(fastq_filename)

    print(seq_list[0:num_seqs_show])

    sys.exit(0)  # always good practice to indicate worked ok!



if __name__ == '__main__':
    main()
    
