#!/bin/bash

file_type="GC*genomic.fna"

find . -mindepth 1 -maxdepth 2 -type d | while IFS= read -r dir; do
		
	for file in ${dir}/GC*genomic.fna; do
		outname="${dir}"
		testfile="${dir}_prots.faa"
		assembly="$(basename -- $file)"
		#cd ${dir}/

		echo -e "\nLooking for files matching: $file_type in $dir\n"

		if [[ $file =~ ^.+_MAG_genomic.fna$ ]]; then
			cd ${dir}

			if [ ! -s ${testfile} ]; then
				echo -e "\n\nProcessing metagenome $file\n\n"
				prodigal -i $assembly -d ${outname}_genes.fna -a ${outname}_prots.faa -o ${outname}.out -p meta  2> prodigal.log
				grep -c '>' ${outname}_prots.faa > total_prots.txt
			else
				echo -e "\n\n$assembly has already been processed by prodigal\n\n"
			fi 


			if [ ! -f *dmnd ]; then
				echo -e "\n\nMaking diamond database file for ${outname}\n\n"
				diamond makedb --in ${outname}_prots.faa --db ${outname} 2> diamond.log
				echo "done"
			else
				echo -e "\n\n$assembly has already been processed by diamond\n\n"

			fi
			cd ../

		elif [[ $file =~ ^.+_genomic.fna$ ]]; then
			cd ${dir}
			if [ ! -s ${testfile} ]; then
				echo -e "Processing genome $file\n\n"	
				prodigal -i $assembly -d ${outname}_genes.fna -a ${outname}_prots.faa -o ${outname}.out 2> prodigal.log
				cat prodigal.log
				grep -c '>' ${outname}_prots.faa > total_prots.txt
			else
				echo -e "\n$assembly has already been processed by prodigal\n"
			fi


			if [ ! -f *dmnd ]; then
				echo -e "\n\nMaking diamond database file for ${outname}\n\n"
				diamond makedb --in ${outname}_prots.faa --db ${outname} 2> diamond.log
				cat diamond.log
			else
				echo -e "\n\n$assembly has already been processed by diamond\n\n"
			fi

			cd ../
		else
			continue
		fi
	done	
		#cd - || exit
	#fi
done

