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
import argparse

def retrieve_esummary(id, db):

	## Get esummary for an entrez id

	esummary_handle = Entrez.esummary(db=db, id=id, report="full")
	esummary_record = Entrez.read(esummary_handle)
	esummary_handle.close()
	return esummary_record


def get_assemblies(email, term, prefix, source, coverage_threshold, retmax, download, path='assemblies'):

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
    filtered_accessions = list()
    removed_accessions = list()
    links = list()
    info = dict()

	## Make a new directory based on provided path name
    directory = path
    if not os.path.exists(directory):
        os.mkdir(directory)
        print("Director does not exist")
    else:
        print("Directory exists")
	
	## Move into new directory
    new_folder = Path(f"./{directory}/")
    os.chdir(new_folder)
	

    for uid in uids:
	    ## Get summary
        summary = retrieve_esummary(id=uid, db="assembly")
        summary.str = json.dumps(summary, indent=4)
        #print(json.dumps(summary['DocumentSummarySet']['DocumentSummary'][0]['ExclFromRefSeq']))
        summary_info = summary['DocumentSummarySet']['DocumentSummary'][0]
        partial_genome = summary['DocumentSummarySet']['DocumentSummary'][0]['PartialGenomeRepresentation']
        accession = summary_info['AssemblyAccession']
        coverage = float(summary['DocumentSummarySet']['DocumentSummary'][0]['Coverage'])
        #print(summary['DocumentSummarySet']['DocumentSummary'][0])

		## Create separate lists for accessions from high quality assemblies and those that are not
        if coverage >= float(coverage_threshold):
            filtered_accessions.append(summary_info['AssemblyAccession'])
        else:
            removed_accessions.append(summary_info['AssemblyAccession'])
	
		# Populate information dict
        assembly_status = ('assemby_status', summary_info['AssemblyStatus'])
        organism = ('organism' , summary_info['Organism'])
        assembly_name = ('assembly_name', summary_info['AssemblyName'])
        coverage = ('coverage', summary_info['Coverage'])
        taxid = ('taxid' , summary_info['SpeciesTaxid'])	
		## Make a variables for obtaining strain information below
        extra_info = summary_info['ExclFromRefSeq']
        strain_info = summary_info['Biosource']['InfraspeciesList']
        info[accession] = [organism, taxid, assembly_name, coverage, assembly_status]
		
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
	
        excluded=summary['DocumentSummarySet']['DocumentSummary'][0]['ExclFromRefSeq']

        ## Retrieve asseblies RefSeq record
        if source == 'refseq':

            # Get ftp link
            if (download == 'nonfiltered'  and float(summary_info['Coverage']) >= float(coverage_threshold)):
                if excluded != []:
                    url = summary['DocumentSummarySet']['DocumentSummary'][0]['FtpPath_RefSeq']
                    label = os.path.basename(url)
                    link = os.path.join(url,label+'_genomic.fna.gz')
                    links.append(link)
                    urllib.request.urlretrieve(link, f'{label}.fna.gz')
                    print(f'Downloading RefSeq file....\n{link}\n')
                else:
                    print(f'No RefSeq file for {accession}, due to {excluded}. Skipping')

            elif (download == 'filtered' and float(summary_info['Coverage']) >= float(coverage_threshold)):
                if (excluded != []):
                    print(f'No RefSeq URL available for {accession}\nInfo:{excluded}\nDownloading Genbank file because of -d nonfiltered flag.')
                    url = summary['DocumentSummarySet']['DocumentSummary'][0]['FtpPath_GenBank']
                    label = os.path.basename(url)
                    link = os.path.join(url,label+'_genomic.fna.gz')
                    links.append(link)
                    urllib.request.urlretrieve(link, f'{label}.fna.gz')
                else:
                    print(f'No RefSeq available for {accession} due to {excluded}. Skipping.')

        ## Check for bad quality assemblies, retrieve assembly if desired using GenBank record
        elif source == 'genbank':
            
            # Get ftp link
            if (download == 'nonfiltered' and float(summary_info['Coverage']) >= float(coverage_threshold)):
                if (excluded != []):
                    print(f'FYI {accession} was flagged by NCBI due to:{excluded}\nDownloading Genbank file anyway because of nonfiltered flag.')
                    url = summary['DocumentSummarySet']['DocumentSummary'][0]['FtpPath_GenBank']
                    label = os.path.basename(url)
                    link = os.path.join(url,label+'_genomic.fna.gz')
                    links.append(link)
                    urllib.request.urlretrieve(link, f'{label}.fna.gz')                
                else:
                    url = summary['DocumentSummarySet']['DocumentSummary'][0]['FtpPath_GenBank']
                    label = os.path.basename(url)
                    link = os.path.join(url,label+'_genomic.fna.gz')
                    links.append(link)
                    urllib.request.urlretrieve(link, f'{label}.fna.gz')                
                    print(f'Downloading GenkBank file....\n{link}\n')

            elif (download == 'filtered' and float(summary_info['Coverage']) >= float(coverage_threshold)):
                if (excluded != []):
                    print(f'{accession} was flagged by NCBI due to:{excluded}\nSkipping download because of -d filtered flag.')
                else:
                    url = summary['DocumentSummarySet']['DocumentSummary'][0]['FtpPath_GenBank']
                    label = os.path.basename(url)
                    link = os.path.join(url,label+'_genomic.fna.gz')
                    links.append(link)
                    urllib.request.urlretrieve(link, f'{label}.fna.gz')                
                    print(f'Downloading GenkBank file....\n{link}\n')
           # else:
           #     print(f'\tNo RefSeq file available for:\n\n\t {accession}\n\tProbably from metagenome?\n\tAdditional info from NCBI: {extra_info} \n\t{coverage} \n\t{assembly_status}.\n\t Downloading GenBank file if coverage >{coverage_threshold}.\n\tCheck metadata file for info.\n')

            #elif url == '' and download == 'filtered':
            #    print(f'No RefSeq file available for:\n\n\t{accession}\n\tAdditional info from NCBI:{extra_info}\n\t{coverage}\n\t{assembly_status}\n\tTo download anyway, rerun and set download to --nonfiltered.\n')

#            ## Get FASTA link
#           label = os.path.basename(url)
#           link = os.path.join(url,label+'_genomic.fna.gz')
            
#           ## Download FASTA
#           print(f'{link}\n\n')
#           links.append(link)
#           urllib.request.urlretrieve(link, f'{label}.fna.gz')

	## Write accession numbers to output <txt> file
    accessions_file = prefix + "_accessions.txt"
	
	## Create a dataframe from info dict to make csv of strain information
    rows = []
    for row_key, inner_list in info.items():
        for col_key, value in inner_list:
            rows.append({'Row': row_key, 'Column': col_key, 'Value': value})
    df = pd.DataFrame(rows)
    df_pivoted = df.pivot(index='Row', columns='Column', values='Value')
    df_pivoted.reset_index(inplace=True)
    df_pivoted.index.name = None
    df_pivoted.rename(columns={'Row': 'accession'}, inplace=True)
    print(f'Search results saved to {prefix}_assembly_info.csv\n\n {df_pivoted})')

	# Write the DataFrame to a CSV file
    df_pivoted.to_csv(f'{prefix}_assembly_info.csv', index=False)

    with open(accessions_file, "w") as acc_output:
        if len(removed_accessions) > 0:
            print(f'Some accessions were of low quality and removed from final list (more info in CSV file):\n\n')
            accessions_delim = '\n'.join(filtered_accessions)
            print('\n'.join(removed_accessions))
			#output.write(f'# Total number of genomes retrieved with esearch term {term}: {total_records}.\n')
            acc_output.write(accessions_delim)
        else:	
            accessions_delim = '\n'.join(filtered_accessions)
			#output.write(f'# Total number of genomes retrieved with esearch term {term}: {total_records}.\n')
            acc_output.write(accessions_delim)

    print(f'\n{total_records} total UIDs accessed from "assembly" database using search term {term}.\nData for {len(filtered_accessions)} records written to ./accessions/\n')
	
    handle.close()
    return links, filtered_accessions, accessions_file

def main():
	
    parser = argparse.ArgumentParser(description="Script for downloading genomes and their metadata from NCBI Entrez browser.")
    parser.add_argument('-e', '--email', type=str, required=True, help='Email to provide server.')
    parser.add_argument('-t', '--term', type=str, required=True, help='A search term for searching NCBI assembly database. Use single quotes with spaces.')
    parser.add_argument('-p', '--prefix', type=str, required=True, help='Prefix appended to output files.')
    parser.add_argument('-s', '--source', choices=['genbank', 'refseq'], nargs='?', default='genbank', const='genbank', required=False, help='Choose database for retrieval. Default is GenBank.')
    parser.add_argument('-c', '--coverage', type=str, required=False, default=5, help='Set coverage threshold for assemblies to download.')
    parser.add_argument('-m', '--retmax', type=int, required=False, default=20, help='Set max number of assemblies to download. Default is 20. Can run script without downloading to see total number of available assemblies in database then set --retmax to that value if downloading all assemblies is desired.')
#    parser.add_argument('-d', '--download', action=argparse.BooleanOptionalAction, required=False, default=False, help='Should assemblies be downloaded into current directory? Defaults to --no-download.')
    parser.add_argument('-d', '--download', choices=['none', 'nonfiltered', 'filtered'], nargs='?', default='none', const='none', required=False, help='Choose whether you want any and all assemblies (nonfiltered) or only those that are not flagged as bad quality by NCBI (filtered).')
    args = parser.parse_args()
    print(args)
    get_assemblies(args.email, args.term, args.prefix, args.source, args.coverage, args.retmax, args.download)	
	
if __name__ == '__main__':
    main()
