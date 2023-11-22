#!/usr/bin/env python3

from Bio import SeqIO
from Bio import Entrez
import subprocess
import sys


def entrez_retrieve_accessions(db, term, output_file, retmax):

	## retrieves records from Entrez database

	Entrez.email = "halvorsen1@llnl.gov"
	handle = Entrez.esearch(idtype="acc", db=db, term=term, retmax=retmax)
	record = Entrez.read(handle)
	
	## esearch returns uid numbers (unique identifiers). Use the uid to retrieve
	## accession numbers for each genome with efetch, or from a different database
	## with elink
	
	uids = record['IdList']
	total_records = len(uids)
	if db == "genome":
		
		# obtain nucleotide records from genome search

		genome_handle = Entrez.elink(dbfrom="genome", db="nucleotide", from_uid=uids)
		genome_record = Entrez.read(genome_handle)

		## extract nucleotide database accessions from the Entrez.elink() object
		genome_accessions = genome_record[0]['LinkSetDb'][0]['Link']
		ids = [acc['Id'] for acc in genome_accessions]
		
		genome_handle.close()

		accessions = list()
		

		## use efetch to get primary accession numbers from uids

		for acc in ids:
			# retrieve record                                          
			acc_handle = Entrez.efetch(db='nucleotide', id=acc, retmode="xml")   
			acc_record = Entrez.read(acc_handle)
			#print(acc_record[0]['GBSeq_primary-accession'])
			# make list of primary accessions                          
			accessions.append(acc_record[0]['GBSeq_primary-accession'])
			acc_handle.close()                                         
        

		## write accession numbers to output <txt> file

		with open(output_file, "w") as output:
			accessions_delim = '\n'.join(accessions)
			#output.write(f'# Total number of genomes retrieved with esearch term {term}: {total_records}.\n')
			output.write(accessions_delim)
		 

		print(f'{total_records} total UIDs accessed from nucleotide db  and {len(accessions)} resulting accessions written to {output_file}')
		
		return accessions, output_file

	elif db == "protein":
		
		## make list of uids from elink for retrieving accessions from another database

		#prot_handle = Entrez.elink(dbfrom="protein", db=db, from_uid=uids)
		#prot_record = Entrez.read(prot_handle)
		#print(prot_record)
		#prot_accessions = prot_handle[1]['LinkSetDb'][0]['Link']
		#accessions = [acc['Id'] for acc in prot_accessions]
		#total_prot = len(prot_accessions)

		#protein_handle.close()

		## Obtain accession numbers for nonredundant protein records with efetch

		accessions = list()

		for uid in uids:

			# retrieve record

			acc_handle = Entrez.efetch(db=db, id=uid, retmode="xml")
			acc_record = Entrez.read(acc_handle)
			
			# make list of primary accessions

			accessions.append(acc_record[0]['GBSeq_primary-accession'])
			acc_handle.close()

		## for now, return the accessions from Entrez.esearch() simply using the user-
		## provided database. can change with the commented out code above.

		with open(output_file, "w") as output:
			accessions_delim = '\n'.join(accessions)
			output.write(accessions_delim)
		
		
		print(f'{total_records} total UIDs accessed from {db} using esearch term {term} and {len(accessions)} resulting accessions written to {output_file}')	
	
		return accessions, output_file

	else:
		sys.stderr.write(usage)
		sys.exit(1)
	
	handle.close()

#def blastx_on_accession_list(accession_file):
#	with open(accession_file, "r") as accession_file_obj:
#		for line in accession_file_obj:

			## create list of accessions from input file
#			accessions.append()

#		for accession in accessions:
			## pull from database with accession handle
#			retrieve = 
#			subprocess.run(retrieve, shell=True)

			## create local blast database from retrieved genome
#			makedb = 

			## blast the translated query against translated ORFs in genome
#			blastx = f'blastx -query {} -db {} -out{}'
		
def main():
	
	usage = f'\n\n\tusage: {sys.argv[0]} "database type" "search term"  "output_file.txt" max # of results to return (optional, default is 20)\n\n\texamples:\n\n\tdatabase type: "genome" or "protein"\n\tsearch term: "Pseudomonas[Orgn] AND luxR[Gene]"\n\n\t'


	if len(sys.argv) < 4:
		sys.stderr.write(usage)
		sys.exit(1)

	if len(sys.argv) == 5:	
		retmax = sys.argv[4]
	else:
		retmax = 20

	db = sys.argv[1]
	term = sys.argv[2]
	output_file = sys.argv[3]

	entrez_retrieve_accessions(db, term, output_file, retmax)

if __name__ == '__main__':
	main()
