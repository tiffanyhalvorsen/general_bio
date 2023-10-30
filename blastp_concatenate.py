#!/usr/bin/env python3

import sys
from pathlib import Path


def blastp_concatenate(input_file_path, input_file_pattern, output_file):

    blastp_dict = {}
    headers = set()
    files = Path(input_file_path).glob(input_file_pattern)
    header_flag = 0

    with open(output_file, "w") as output:
        for file in files:
            with open(file, "r") as input_file:
                for line in input_file:
                    if line.startswith('# Field:'):
                        line.lstrip('# Field:')
                        line.rstrip()
                        headers = line.split(', ')
                        if header_flag == 0:
                            fields = headers.join('\t')
                            output.write(f'{fields}\n')
                            header_flag = 1
                    if not line.startswith('#') and not line.startswith('query_id'):
                        #print(line)
                        output.write(line)
                        line = line.rstrip()
                        data = line.split('\t')
                        blastp_dict[data[1]] = dict(zip(headers, data))
    return blastp_dict

def main():
    
    input_file_path = sys.argv[1]
    input_file_pattern = sys.argv[2]
    output_file = sys.argv[3]

    usage = f'\n\n\tusage: {sys.argv[0]} filename.fasta'

    if len(sys.argv) <3:
        sys.stderr.write(usage)
        sys.exit(1)

    blastp_concatenate(input_file_path, input_file_pattern, output_file)

if __name__ == '__main__':
    main()
       
