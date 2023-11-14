#!/usr/bin/env python3

from Bio import SeqIO
from Bio import Entrez
import subprocess
import sys


def entrez_retrieve_accessions(db, term, return_max):

	## retrieves records from Entrez database

	Entrez.email = "halvorsen1@llnl.gov"
	handle = Entrez.esearch(db=db, term=term, retmax=return_max)
	record = Entrez.read(handle)
	
	## esearch returns uid numbers (unique identifiers). Use the uid to retrieve
	## accession numbers for each genome

	uids = record['IdList']
	print(uids)
	if db == "genome":
		
		# obtain nucleotide records from genome search

		genome_handle = Entrez.elink(dbfrom="genome", db="nuccore", from_uid=uids)
		genome_record = Entrez.read(genome_handle)

		# print('\n\n',genome_record)
		
		## extract nucleotide database accessions from the Entrez.elink() object
		
		genome_accessions = genome_record[1]['LinkSetDb'][0]['Link']
		ids = [acc['Id'] for acc in genome_accessions]
		print('\n',genome_accessions,'\n\n',ids,'\n\n',len(ids))
		total_genomes = len(ids)
		genome_handle.close()

		accessions = list()
		

		## figure out how to get primary accessions with efetch for genomes
	
		for acc in ids:
			# retrieve record                                          
			acc_handle = Entrez.efetch(db=db, id=acc, retmode="xml")   
			acc_record = Entrez.read(acc_handle)
			print(acc_record)                                                   
			# make list of primary accessions                          
			#accessions.append(acc_record[0]['GBSeq_primary-accession'])
			acc_handle.close()                                         
                                                               
		print(accessions)                                              

		return accessions, total_genomes

	elif db == "protein":
		

		#prot_handle = Entrez.elink(dbfrom="protein", db=db, from_uid=uids)
		#prot_record = Entrez.read(prot_handle)
		#print(prot_record)
		#prot_accessions = prot_handle[1]['LinkSetDb'][0]['Link']
		#accessions = [acc['Id'] for acc in prot_accessions]
		#total_prot = len(prot_accessions)

		#print('\n',record, '\n\n', total_prot, '\n\n', accessions)

		#protein_handle.close()

		## Obtain accession numbers for nonredundant protein records with efetch
		## This chunk works

		accessions = list()

		for uid in uids:
			# retrieve record
			acc_handle = Entrez.efetch(db=db, id=uid, retmode="xml")
			acc_record = Entrez.read(acc_handle)

			# make list of primary accessions
			accessions.append(acc_record[0]['GBSeq_primary-accession'])
			acc_handle.close()

		print(accessions)
		
		## for now, return the accessions from Entrez.esearch() simply using the user-
		## provided database

		return accessions
		
	handle.close()

#def entrez_retrieve_seqs(acc_list):

	

accessions = list()

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
	
	usage = f'\n\n\tusage: {sys.argv[0]} "database type" "search term"  "ID type"\n\n\texamples:\n\n\tdatabase type: "pubmed" or "nucleotide"\n\tsearch term: "Pseudomonas[Orgn] AND luxR[Gene]"\n\n\toptional 3rd parameter: max # of results to return (default is 20)\n\n'


	if len(sys.argv) < 3:
		sys.stderr.write(usage)
		sys.exit(1)

	if len(sys.argv) == 4:	
		return_max = sys.argv[3]
	else:
		return_max = 20

	db = sys.argv[1]
	term = sys.argv[2]

	entrez_retrieve_accessions(db, term, return_max)	


if __name__ == '__main__':
	main()
