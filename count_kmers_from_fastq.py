#!/usr/bin/env python

import os, sys

from seq_to_kmer_list import *
from seq_list_from_fastq import *


## method: count_kmers(kmer_list)
##
##  Counts the frequency of each kmer in the given list of kmers
##
##  input parameters:
##
##  kmer_list : list of kmers (type: list)
##               ie.  ["GATC", "TCGA", "GATC", ...]
##
##
##  returns kmer_counts_dict : dict containing ( kmer : count )
##                    ie.  {  "GATC" : 2,
##                            "TCGA" : 1,
##                             ...       }


def count_kmers(kmer_list):

    kmer_count_dict = {}
    for kmer in kmer_list:
        if kmer not in kmer_count_dict:
            kmer_count_dict[kmer] = 1
        else:
            kmer_count_dict[kmer] += 1


    return kmer_count_dict


def main():

    progname = sys.argv[0]

    usage = "\n\n\tusage: {} filename.fastq kmer_length num_top_kmers_show\n\n\n".format(
        progname
    )

    if len(sys.argv) < 4:
        sys.stderr.write(usage)
        sys.exit(1)

    # capture command-line arguments
    fastq_filename = sys.argv[1]
    kmer_length = int(sys.argv[2])
    num_top_kmers_show = int(sys.argv[3])

    # call on previous functions to create seq list and kmer list
    seq_list = seq_list_from_fastq(fastq_filename)

    #######################
    ## Step 1:
    ## begin your code, populate 'all_kmers' list with the
    ## collection of kmers from all sequences
    
    all_kmers = list()
    for seq in seq_list:
       all_kmers.extend(seq_to_kmer_list(seq, 4)) 

    ## create dictionary of counts using list of identified kmers
    ## in provided sequence

    kmer_count_dict = count_kmers(all_kmers)

    ## list of unique kmers, without counts

    unique_kmers = list(kmer_count_dict.keys())

    #########################
    ## Step 3: sort unique_kmers by abundance descendingly
    ## (Note, you can run and test without first implementing Step 3)
    ## begin your code       hint: see the built-in 'sorted' method documentation

    ## sort kmer dictionary

    sorted_kmers = sorted(((value, key) for (key, value) in kmer_count_dict.items()), reverse=True)

    sorted_kmers = dict((k,v) for v,k in sorted_kmers)

    ## printing the num top kmers to show
    top_kmers_show = unique_kmers[0:num_top_kmers_show]

    for kmer in top_kmers_show:
        print("{}: {}".format(kmer, kmer_count_dict[kmer]))

    sys.exit(0)  # always good practice to indicate worked ok!

if __name__ == "__main__":
    main()
