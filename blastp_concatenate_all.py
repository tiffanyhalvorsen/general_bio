#!/usr/bin/env python3

import sys
from pathlib import Path
import re
import os
import fnmatch

def blastp_concatenate(input_file_path, input_file_pattern):

    blastp_dict = {}
    headers = set()
#    input_file_path = str(input_file_path)
    if os.path.exits(Path(input_file_path).glob(input_file_pattern)):
        files = Path(input_file_path).glob(input_file_pattern)
    else:
        sys.exit(f'Files at {input_file_path} do not exist.')

    header_flag = 0

    if input_file_path.startswith('GC'):
        output_file = os.path.join(input_file_path, f'{input_file_path}_dmnd_results.fa')
    else:    
        output_file = os.path.join(input_file_path, 'combined_dmnd_results.fa')

    with open(output_file, "w") as output:
        for file in files:
            with open(file, "r") as input_file:
                for line in input_file:
                    if line.startswith('# Field'):
                        if header_flag == 0:
                            #line.lstrip('# Field')
                            line = line.rstrip()
                            newline = re.sub(r'^(# Field.+?\S)', r'', line)
                       	    headers = newline.split(', ')
                            fields = "\t".join(headers)
                            output.write(f'{fields}\n')
                            header_flag = 1
                        else:
                            continue
                    if not line.startswith('#') and not line.startswith('query_id'):
                        output.write(line)
                        line = line.rstrip()
                        data = line.split('\t')
                        #print(data[1])
                        #blastp_dict[data[1]] = dict(zip(headers, data))
        
    output_check = os.path.getsize(output_file)
    if output_check > 0:
        print(f'{output_file} successfully written')
        return output
    else:
        print('Output file failed to be written')

def blastp_concatenate_between_directories(pattern, output_file):
	
    blastp_dict = {}
    headers = set()
#   input_file_path = str(input_file_path)
    header_flag = 0
    with open(output_file, "w") as output:
       for root, dirnames, filenames in os.walk('.'):
           matches = fnmatch.filter(filenames, pattern)
           #dir_list = [d for d in dirnames if d.startswith('GC')]
           
           ## Open only files matching user-provided pattern in regex format
           ## Commented out line can be used to search instead for a file extension, not in regex format

           #for filename in [f for f in filenames if f.endswith(pattern)]:
           for filename in [f for f in filenames if f in matches]:
               file = os.path.join(root, filename)
               print(file)
               with open(file, "r") as input_file:
                   for line in input_file:
                       if line.startswith('# Field'):
                            if header_flag == 0:
                                #line.lstrip('# Field')
                                line = line.rstrip()
                                newline = re.sub(r'^(# Field.+?\S)', r'', line)
                                headers = newline.split(', ')
                                fields = "\t".join(headers)
                                output.write(f'{fields}\n')
                                header_flag = 1
                            else:
                                continue
                       if not line.startswith('#') and not line.startswith('query_id'):
                           output.write(line)
                           line = line.rstrip()
                           data = line.split('\t')
                           #print(data[1])
                           #blastp_dict[data[1]] = dict(zip(headers, data))
        
    output_check = os.path.getsize(output_file)
    if output_check > 0:
        print(f'{output_file} successfully written')
        return output
    else:
        print('No output created :(')

def main():
    
    usage = f'\n\n\tThis program can combine blastp results either in a provided directory using a provided file pattern (usage #1) OR it can combine blastp results between subdirectories using a file extension. But usage #2 requires a user-provided filename\n\n\tusage #1: {sys.argv[0]} path/to/files/ pattern\n\n\tusage #2: {sys.argv[0]} "file_pattern" "output_filename"\n\n'
    if len(sys.argv) <2:
        sys.stderr.write(usage)
        sys.exit(1)

    elif len(sys.argv) == 3:

        if sys.argv[1].endswith('/'):
            input_file_pattern = sys.argv[1]
            output_file = sys.argv[2] 
            blastp_concatenate(input_file_path, input_file_pattern)
        else: 
            pattern = sys.argv[1]
            output_file = sys.argv[2]
            blastp_concatenate_between_directories(pattern, output_file)

       # else:
       #     sys.stderr.write(usage)
       #     sys.exit(1)

    else:
        sys.stderr.write(usage)
        sys.exit(1)

if __name__ == '__main__':
    main()
       
