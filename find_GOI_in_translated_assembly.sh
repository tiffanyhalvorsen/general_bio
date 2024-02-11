#!/bin/bash

## This script finds all subdirectories in the working directory (meant to be 'assemblies/[all assemlies]/' 
## and looks for diamond database files in those subdirectories, then aligns all query sequences
## (located in a 'assemblies/queries/' directory) to each assembly and writes the results to an output
## file in each assembly directory.

## usage: find_GOI_in_translated_assembly.sh [optional sensitivity parameter: --sensitive or --very-sensitive] [optional: --current-dir to save output files in current directory. If left blank, they are saved in each subdirectory where the assemblies are located]


sensitivity=$1
output_location=$2


find . -mindepth 1 -maxdepth 1 -type d | while IFS= read -r dir; do
	
	for db_file in $(find $dir -maxdepth 1 -type f -name "*.dmnd"); do

		for query in $(find ./queries -type f -name "*.faa"); do
			queryname=$(basename $query .faa)
			dbname=$(basename $db_file .dmnd)
			if [[ $output_location == --current-dir ]]; then					
				diamond blastp --query $query --db $db_file $sensitivity --outfmt 6 --header --out ${queryname}_blastp_results_${dbname}.out
			else
				diamond blastp --query $query --db $db_file $sensitivity --outfmt 6 --header --out $dir/$(basename -- $query)_blastp_results_$(basename -- $db_file).out
			fi
		done		
	done
done
