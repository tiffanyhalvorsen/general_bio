#!/usr/bin/env python3

import sys
import re

    ## This program takes a dna sequence and re-formats it to the desired line width

def format_dna(dna, width):
    dna = re.sub(r"\n",r"",dna)
    regex = r"\w" + '{10}'
    new_dna = re.sub(regex,r"\1\n",dna, re.IGNORECASE) 		
    return new_dna

def main():

    ## usage message
     
    usage = f'\n\n\tusage: {sys.argv[0]} DNA_sequence desired_width'

    if len(sys.argv) < 2:
        sys.stderr.write(usage)
        sys.exit(1)

    ## collect command line arguments

    dna = sys.argv[1]
    width = sys.argv[2]

    formatted_dna = format_dna(dna, width)

    print(formatted_dna)

if __name__ == "__main__":
   main()


