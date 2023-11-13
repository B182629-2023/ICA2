#!/usr/bin/python3

import os
import sys
import pandas as pd

#os.system("wget ftp://ftp.expasy.org/databases/prosite/prosite.dat")

#sh -c "$(wget -q ftp://ftp.ncbi.nlm.nih.gov/entrez/entrezdirect/install-edirect.sh -O -)"

tax_group = input("Enter taxonomy group\n")


# Asks the user if their taxonomy group input is correct 
#while True:
#	while True:
#		tax_group = input("Enter taxonomy group\n")
#		tax_group_correct = input(f"You entered '{tax_group}'. Is this correct? (yes/no)\n")
#		if tax_group_correct == "yes":
#			print(f"\nQuerying taxonomic group '{tax_group}'...\n")
#			break
#		elif tax_group_correct == "no":
#			print(f"\nResetting taxonomic group...\n")
#		else:
#			print("\nInvalid input. Please enter yes or no.\n")
#
#	query = f"'{tax_group}[Organism:exp] AND swissprot[filter]'"
#	os.system("esearch -db protein -query " + query + " | efetch -format acc > accessions")
#	accession_id = open('accessions').read().rsplit()
#	print("There are " + str(len(accession_id)) + " protein entries for " + tax_group)
#	accession_no_continue = input(f"Would you like to continue querying {tax_group}? (yes/no)\n")
#	if accession_no_continue == "yes":
#		print(f"Fetching protein families/superfamilies... (this may take several minutes)")
#		os.system("esearch -db protein -query " + query + " | efetch -format gp | grep -A 1 \"SIMILARITY\" | awk '{sub(/^[ \\t]+/,\"\");}1' > pf_raw.txt")
#		break
#	elif accession_no_continue == "no":
#		print(f"\nResetting taxonomic group...\n")
#	else:
#		print("\nInvalid input. Please enter yes or no.\n")



# Create a clean file including protein family names for taxonomic group
#pf_clean_open = open("pf_clean.txt", "w")
#pf_replaced_open = open("pf_raw.txt").read().replace("\n"," ").replace("{ECO","\n").replace("[SIMILARITY]","\n").replace(" Belongs to the ","").replace("--","").replace(";","\n") 
#pf_clean_open.write(pf_replaced_open)
#pf_clean_open.close()


# Writes lines containing string 'family' or 'superfamily' to new file
#string = ' family'
#with open("pf_clean.txt", 'r') as input, open("protein_families.txt", 'w') as output:
#	pf_clean_read = input.read()
#	rows = pf_clean_read.split('\n')
#	cleaned_rows = [row for row in rows if string in row]
#	new_rows = '\n'.join(cleaned_rows)
#	output.write(new_rows + '\n')


# Separates protein family name and protein superfamily name and writes lines to a new file
#separated_open = open("separated.txt", 'w')
#separated_replaced_open = open("protein_families.txt").read().replace(". ","\n").replace(" belongs to the ", "")
#separated_open.write(separated_replaced_open)
#separated_open.close()

#Writes lines containing protein family name into a 'family' file
#with open("separated.txt", 'r') as input, open("families.txt", 'w') as output:
    #    pf_count_read = input.read()
   #     rows = pf_count_read.split('\n')
  #      cleaned_rows = [row for row in rows if string in row]
 #      new_rows = '\n'.join(cleaned_rows)
 #       output.write(new_rows + '\n')

#Counts the number of lines containing each protein family name i.e. counts the number of proteins belonging to each family
#with open("families.txt", 'r') as input, open("families_count.txt", 'w') as output:
#        rows = [row.rstrip() for row in input]
#        identical_rows = set(rows)
#        for row in identical_rows:
#                count = rows.count(row)
#                output.write(row + '\t' +  str(count) + '\n')
#print("Family protein names collated.")

#Writes lines containing protein superfamily name into a 'superfamily' file
#string2 = 'superfamily'
#with open("separated.txt", 'r') as input, open("superfamilies.txt", 'w') as output:
#	pf_count_read = input.read()
#	rows = pf_count_read.split('\n')
#	cleaned_rows = [row for row in rows if string2 in row]
#	new_rows = '\n'.join(cleaned_rows)
#	output.write(new_rows + '\n')

#Counts the number of lines containing each protein superfamily name (if present)
#with open("superfamilies.txt", 'r') as input, open("superfamilies_count.txt", 'w') as output:
#        rows = [row.rstrip() for row in input]
#        identical_rows = set(rows)
#        for row in identical_rows:
#                count = rows.count(row)
#                output.write(row + '\t' +  str(count) + '\n')

#print("Superfamily names collated.")

family_df = pd.read_csv("families_count.txt",sep='\t')

family_df.columns = ['Protein Family', 'Protein Count']
row_number = len(family_df)
print("There are " + str(row_number) + " protein families in " + tax_group + ".")

while True:
	table = input("Would you like to view the protein families present in " + tax_group + "? (yes/no)\n")
	if table == "yes":
		while True:
			table_size = input("Enter the number of protein families with the most proteins you would like to view. (e.g. 5, 10, 50 etc.)\n")
			if int(table_size) <= int(row_number):
				table_size_df = family_df.sort_values('Protein Count', ascending=False).head(int(table_size))
				print(table_size_df)
				break
			elif int(table_size) > int(row_number):
				print("There are " + str(row_number) + " protein families in " + tax_group + ". Please enter a smaller number.")
			else:
				print("Invalid input.")
		break
	elif table == "no":
		break
	else:
		print("Invalid input.")

while True:
	pfam_interest = input("Enter protein family of interest")
	pfam = pfam_interest + 
	if pfam_interest in set(family_df['Protein Family']):
		pfam_interest_size = family_df[family_df['#Protein Family'].str.contains(pfam_interest)]['Protein Count']
		print("There are " + str(pfam_interest_size) + " proteins in the " + pfam_interest + ".") 
		print("Searching for proteins in the " + pfam_interest + " in " + tax_group + "...")
		break
	else pfam_interest not in set(family_df['Protein Family']): 
		print("There are 0 protein families with the name " + pfam_interest + " in " + tax_group + ".")



#protein_families_count_open = open("protein_families_count.txt").read()
#print(protein_families_count_open)


#while True:
#	taxonomy_group = input("Enter taxonomy group\n")
#	protein_family = input("Enter protein family\n")
#	query = "{taxonomy_group}[Organism] AND {protein_family}[Protein Name]"
#	os.system("esearch -db protein -query {query} | efetch -format acc >  accessions")
#	accession_id = open('accessions').read().rsplit()
#	print(f"There are {len(accession_id)}results for {protein_family} in {taxonomy_group}")
#	os.system("esearch -db protein -query {query} | efetch -format gp | grep Accession\|Description > gps")
#	descriptions = open('gps').read()
#	view = input("Would you like to view the proteins and organisms generated from query? (yes/no)")
#	if view == "yes":
#		print({descirptions})
#	elif view == "no":
#		print(f"")
#	else:
#		print(f"Invalid input. Please type yes or no.")
#	proceed = input("Would you like to proceed with this query? (yes/no)")
#	if proceed == "yes":
#		print(f"Proceeding with results from  query")
#		break
#	elif proceed == "no":
#		print(f"resetting variables...")
#	else:
#		print(f"Invalid input.")

#os.system("esearch -db protein -query {query} | efetch -format fasta | grep -v '>' | clustalo - - --outfmt=clustal")


