#!/usr/bin/env python

import os, sys



## method: sequence_to_kmer_list(sequence, kmer_length)
##
##  Extracts all kmers of a specified length from a sequence
##
##  ie.  sequence: GATCGATCGATCGA
##   and given kmer_length = 4
##   would yield
##                 GATC
##                  ATCG
##                   TCGA
##                    .... and so forth
##                       
##  input parameters:
##
##  seuqence : nucleotide sequence (type: string)
##
##  returns kmer_list : list of kmer sequences.
##                    ie.  ["GATC", "ATCG", ...]
    
def seq_to_kmer_list(sequence, kmer_length):

    kmers_list = list()

    ## iterate over length of sequence, grabbing sequence
    ## of length kmer + 1 at each iteration, add to list of kmers.
    ## filter out the tail end of sequence that results in kmers 
    ## shorter than desired length.

    for position in range(len(sequence)):
        kmer = sequence[position:(position+kmer_length)]
        if len(kmer) == kmer_length:
            kmers_list.append(kmer)
            position += 1

    return kmers_list



def main():

    progname = sys.argv[0]
    
    usage = "\n\n\tusage: {} sequence kmer_length\n\n\n".format(progname)
    
    if len(sys.argv) < 3:
        sys.stderr.write(usage)
        sys.exit(1)

    # capture command-line arguments
    sequence = sys.argv[1]
    kmer_length = int(sys.argv[2])

    kmers  = seq_to_kmer_list(sequence, kmer_length)

    print(kmers)

    sys.exit(0)  # always good practice to indicate worked ok!



if __name__ == '__main__':
    main()
