#!/usr/bin/env python3

from Bio import SeqIO
from Bio import Entrez
from pathlib import Path
import subprocess
import sys
import os
import urllib
import json
import pandas as pd

def retrieve_esummary(id, db):

	## Get esummary for an entrez id

	esummary_handle = Entrez.esummary(db=db, id=id, report="full")
	esummary_record = Entrez.read(esummary_handle)
	esummary_handle.close()
	return esummary_record


def get_assemblies(email, term, output_prefix, retmax, download=False, path='assemblies'):

	## Download genbank assemblies for a search term
	## term: usually organism name
	## download: default is True
	## path: where to save assemblies

	Entrez.email = email
	
	## Esearch returns uid numbers (unique identifiers). Use the uid to retrieve
	## accession numbers with esummary.
	
	handle = Entrez.esearch(db="assembly", term=term, retmax=retmax)
	record = Entrez.read(handle)
	total_hits = record['Count']
	print(f'\nTotal hits in genome assembly database: {total_hits}\n')
	uids = record['IdList']
	total_records = len(uids)
		
	## ALTERNATIVE: RETRIEVE ACCESSIONS FROM A DIFFERENT DATABASE WITH ELINK

	# genome_handle = Entrez.elink(dbfrom="genome", db="nucleotide", from_uid=uids)
	# genome_record = Entrez.read(genome_handle)

	## Extract nucleotide database accessions from the Entrez.elink() object
		
	# genome_accessions = genome_record[0]['LinkSetDb'][0]['Link']
	# ids = [acc['Id'] for acc in genome_accessions]
		
	# genome_handle.close()

	# for acc in uids:
		
		## Retrieve record                                          
		
		# acc_handle = Entrez.efetch(db='nucleotide', id=acc, retmode="xml")   
		# acc_record = Entrez.read(acc_handle)
		#print(acc_record[0]['GBSeq_primary-accession'])
		
		##  Make list of primary accessions                          
		
		# accessions.append(acc_record[0]['GBSeq_primary-accession'])
		# acc_handle.close()                                         
	accessions = list()
	links = list()
	info = dict()

	## Make a new directory based on provided path name
	directory = path
	if not os.path.exists(directory):
		os.mkdir(directory)
	
	## Move into new directory
	new_folder = Path(f"./{directory}/")
	os.chdir(new_folder)
	

	for uid in uids:
		## Get summary
		summary = retrieve_esummary(id=uid, db="assembly")
		## Add accession to list
		accessions.append(summary['DocumentSummarySet']['DocumentSummary'][0]['AssemblyAccession'])
		
		# Populate information dict
		accession = summary['DocumentSummarySet']['DocumentSummary'][0]['AssemblyAccession']
		organism = ('organism' , summary['DocumentSummarySet']['DocumentSummary'][0]['Organism'])
		taxid = ('taxid' , summary['DocumentSummarySet']['DocumentSummary'][0]['SpeciesTaxid'])
		## Make a variables for obtaining strain information below
		extra_info = summary['DocumentSummarySet']['DocumentSummary'][0]['ExclFromRefSeq']
		strain_info = summary['DocumentSummarySet']['DocumentSummary'][0]['Biosource']['InfraspeciesList']
		info[accession] = [organism, taxid]
		#print(summary['DocumentSummarySet']['DocumentSummary'][0])
		#print(f'\n\n')
		
		## Populate dictionary with strain information.
		## Completed assemblies will have a strain code as 'Sub_value', and a RefSeq accession
		## Assemblies from metagenome projects (or else, incomplete) will not have a RefSeq
		## accession and the dictionary is populated with information on why they are 
		## exluded from the RefSeq database (taken from "extra_info" variable, which mapes to the 
		## 'ExclFromRefSeq' information).

		if len(strain_info) > 0:
			strain = ('strain' , strain_info[0]['Sub_value'])
		else:
			if len(extra_info) > 1:
				str1 = extra_info[0]
				str2 = extra_info[1]
				str3 = summary['DocumentSummarySet']['DocumentSummary'][0]['Biosource']['Isolate']
				strain = ('strain' , f'{str1}; {str2}; {str3}')
			elif len(extra_info) == 1:
				str1 = extra_info[0]
				str3 = summary['DocumentSummarySet']['DocumentSummary'][0]['Biosource']['Isolate']
				strain = ('strain' , f'{str1}; {str3}')
			else:
				str3 = summary['DocumentSummarySet']['DocumentSummary'][0]['Biosource']['Isolate']
				strain = ('strain' , f'{str3}')
		info[accession].append(strain)
	
		if download==True:
			print("Downloading...\n\n")
			## Get ftp link
			url = summary['DocumentSummarySet']['DocumentSummary'][0]['FtpPath_RefSeq']
			if url == '':
				url = summary['DocumentSummarySet']['DocumentSummary'][0]['FtpPath_GenBank']
				print(f'\tNo RefSeq file available for:\n\n\t {info[accession]}.\n\t Downloading GenBank file.\n\t Check metadata file for info.\n')
	
			## Get FASTA link
			label = os.path.basename(url)
			link = os.path.join(url,label+'_genomic.fna.gz')
			print(f'{link}\n')
			links.append(link)
			urllib.request.urlretrieve(link, f'{label}.fna.gz')
	
	## Write accession numbers to output <txt> file
	prefix = output_prefix
	accessions_file = prefix + "_accessions.txt"
 	
	## Create a dataframe from info dict to make csv of strain informatio
	rows = []
	for row_key, inner_list in info.items():
		for col_key, value in inner_list:
			rows.append({'Row': row_key, 'Column': col_key, 'Value': value})
	df = pd.DataFrame(rows)
	df_pivoted = df.pivot(index='Row', columns='Column', values='Value')
	df_pivoted.reset_index(inplace=True)
	df_pivoted.index.name = None
	df_pivoted.rename(columns={'Row': 'accession'}, inplace=True)
	print(f'Search results (saved to {output_prefix}_assembly_info.csv)\n\n {df_pivoted}') 
	
	# Write the DataFrame to a CSV file
	df_pivoted.to_csv(f'{output_prefix}_assembly_info.csv', index=False)
	
	with open(accessions_file, "w") as acc_output:
		accessions_delim = '\n'.join(accessions)
		#output.write(f'# Total number of genomes retrieved with esearch term {term}: {total_records}.\n')
		acc_output.write(accessions_delim)

	print(f'\n{total_records} total UIDs accessed from "assembly" database using search term {term}. Data for {len(accessions)} records written to ./accessions/\n')
	
	handle.close()
	return links, accessions, accessions_file

def main():
    
    usage = f'\n\tThis script can be used to access Entrez Assembly database and query it for strains of interest.\nResults can be used to create a list of accession numbers, written to a <txt> file using input prefix name.\nIf the last parameter provided is "Yes", FASTA files for each accession will be downloaded into "assemblies" in current directory.\n\n\tNOTE: first 3 parameters are required.\n\n\tusage: {sys.argv[0]} \n\n\tParameters:\n\t(1) Email address (e.g, "bionerd@hotmail.com") for Entrez database query.\n\t(2) Search term (e.g., "Pseudomonas[Orgn] AND luxR[Gene]"\n\t(3) Output file prefix (e.g., "pseudomonas-luxR")\n\t(4) Max # of results to return (optional, default is 20)\n\t(5) Do you want to download genome assemblies? Type "Yes" as 5th parameter.\n\nNOTE: accessions.txt file will be over-written if program is re-run in same directory with same prefix.\n\n\n'

    
    if len(sys.argv) < 4:
        sys.stderr.write(usage)
        sys.exit(1)

    elif len(sys.argv) == 4:
        retmax = 20
        email = sys.argv[1]
        term = sys.argv[2]
        output_prefix = sys.argv[3]
        get_assemblies(email, term, output_prefix, retmax)
    
    elif len(sys.argv) == 5:
        email = sys.argv[1]
        term = sys.argv[2]
        output_prefix = sys.argv[3]
        retmax = sys.argv[4]
        get_assemblies(email, term, output_prefix, retmax)
    
    elif len(sys.argv) == 5:
        email = sys.argv[1]
        term = sys.argv[2]
        output_prefix = sys.argv[3]
        retmax = sys.argv[4]
        if sys.argv[5] == "Yes":    
            get_assemblies(email, term, output_prefix, retmax, download=True)
        else:
            sys.stderr.write("Last parameter should be 'Yes' to intitate assembly FASTA download.")
            sys.exit(1)
	
	
if __name__ == '__main__':
	main()
