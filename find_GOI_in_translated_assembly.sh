#!/bin/bash

## This script finds all subdirectories in the working directory (meant to be 'assemblies/[all assemlies]/' 
## and looks for diamond database files in those subdirectories, then aligns all query sequences
## (located in a 'assemblies/queries/' directory) to each assembly and writes the results to an output
## file in each assembly directory.

## usage: find_GOI_in_translated_assembly.sh [optional sensitivity parameter: --sensitive or --very-sensitive] [optional: --current-dir to save output files in current directory. If left blank, they are saved in each subdirectory where the assemblies are located]


sensitivity=$1
output_location=$2

# Search for all subdirectories in current directory

find . -mindepth 1 -maxdepth 1 -type d | while IFS= read -r dir; do
	
	# Loop through each genome diamond database file

	for db_file in $(find $dir -maxdepth 1 -type f -name "GC*.dmnd"); do

		# Also loop through each query of interest in ./queries subfolder

		# Blastp each .faa query file against each .dmnd genome database file
		
		for query in $(find ./queries -type f -name "*.faa"); do
			queryname=$(basename $query .faa)
			dbname=$(basename $db_file .dmnd)

			# Check that diamond has not alreadyu been run for the queries in ./queries subfolder

            if [ ! -s ${dir}/${queryname}_blastp_results_${dbname}.out ]; then 
				echo -e "${queryname}_blastp_results_${dbname}.out -- blasting off now \n\n"

				# User can specify that output gets placed in the current working directory

				if [[ $output_location == --current-dir ]]; then					
					diamond blastp --query $query --db $db_file $sensitiviity --outfmt 6 --header --out ${queryname}_blastp_results_${dbname}.out 2> diamond.log

				# If current directory is not specified, output is placed in subdirectories with genome assembly folders

				else
					diamond blastp --query $query --db $db_file $sensitivity --outfmt 6 --header --out ${dir}/${queryname}_blastp_results_${dbname}.out 2> diamond.log
				fi
				
				assembly=$(basename $dir)

			# Add "Assembly_ID" and "Query_name"to the end of the Fields line and the
			# corresponding variables to the end of each blastp hit in which it was found.

			# Note: sed delimiter can change. Using comma here to avoid pattern matching 
			# conflicts. Two sed commands, separated by ";". Must surround sed command 

			# in double-quotes to expand the $assembly variable.

				sed -i '' -e "/# Fields/ s,$,\tAssembly_ID,;/^# /!s,$,\t$assembly,;/# Fields/ s,$,\tQuery_name,;/^# /! s,$,\t$queryname," ${dir}/${queryname}_blastp_results_${dbname}.out
			
				# If database has already been queried, don't run again

				else
					echo -e  "${queryname} has already been queried against ${dbname}. Skipping.\n\n"
		

			fi		
		
			done		
		done
done	
