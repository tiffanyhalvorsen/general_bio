#!/bin/bash

file_type="GC*genomic.fna"

find . -mindepth 1 -maxdepth 2 -type d | while IFS= read -r dir; do
		
	for file in ${dir}/GC*genomic.fna; do
		outname="${dir}"
		testfile="${dir}_prots.faa"
		assembly="$(basename -- $file)"
		#cd ${dir}/
		if [[ $file =~ ^.+_MAG_genomic.fna$ ]]; then
			cd ${dir}

			if ! [ -f ${testfile} ]; then
				echo "\nProcessing metagenome $file"
				prodigal -i $assembly -d ${outname}_genes.fna -a ${outname}_prots.faa -p meta -t ${outname}_training.txt 2> prodigal.log
				echo "done"
				grep -c '>' ${outname}_prots.faa > total_prots.txt
			else
				echo -e "\n$assembly has already been processed by prodigal\n"
			fi 


			if ! [[ -f *dmnd ]]; then
				echo "Making diamond database file for ${outname}"
				diamond makedb --in ${outname}_prots.faa --db ${outname} 2> diamond.log
				echo "done"
			else
				echo -e "$assembly has already been processed by diamond\n\n"

			fi
			cd ../

		elif [[ $file =~ ^.+_genomic.fna$ ]]; then
			cd ${dir}
			if ! [ -f ${testfile} ]; then
				echo -e "Processing genome $file\n\n"	
				prodigal -i $assembly -d ${outname}_genes.fna -a ${outname}_prots.faa -o ${outname}.out -t ${outname}_training.txt 2> prodigal.log
				cat prodigal.log
				grep -c '>' ${outname}_prots.faa > total_prots.txt
			else
				echo -e "\n$assembly has already been processed by prodigal\n"
			fi


			if ! [[ -f *dmnd ]]; then
				echo "Making diamond database file for ${outname}"
				diamond makedb --in ${outname}_prots.faa --db ${outname} 2> diamond.log
				cat diamond.log
			else
				echo -e "$assembly has already been processed by diamond\n\n"
			fi

			cd ../
		else
			continue
		fi
	done	
		#cd - || exit
	#fi
done

