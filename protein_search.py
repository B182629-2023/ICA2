#!/usr/bin/python3

import subprocess
import sys
import pandas as pd

#os.system("wget ftp://ftp.expasy.org/databases/prosite/prosite.dat")

#sh -c "$(wget -q ftp://ftp.ncbi.nlm.nih.gov/entrez/entrezdirect/install-edirect.sh -O -)"


# Taxonomic group: Asks the user for desired taxonomic group and verifies if input is correct
while True:
	while True:
		tax_group = input("Enter taxonomy group\n")
		tax_group_correct = input(f"\nYou entered '{tax_group}'. Is this correct? (yes/no)\n")
		if tax_group_correct == "yes":
			print(f"\n\nQuerying taxonomic group '{tax_group}'...\n")
			break
		elif tax_group_correct == "no":
			print(f"\nResetting taxonomic group...\n")
		else:
			print("\nInvalid input. Please enter yes or no.\n")

	# Searches for entries with taxnomic group with the swissprot filter and collects accession IDs
	query = f"'{tax_group}[Organism:exp] AND swissprot[filter]'"
	subprocess.call("esearch -db protein -query " + query + " | efetch -format acc > accessions.txt", shell=True)

	# Shows user the number of protein entries for taxonomic group and asks if they want to proceed
	with open('accessions.txt', 'r') as accessions:
		accs = accessions.read().rsplit('\n')
		print("\nThere are " + str(len(accs)) + " protein entries for " + tax_group)

	accession_no_continue = input(f"\nWould you like to continue querying {tax_group}? (yes/no)\n")
	if accession_no_continue == "yes":
		# collects protein family and superfamily information from the GenPept summary and removes unwated whitespace for each protein entry then exits loop
                # sometimes this esearch returns 500 internal server error. If this occurs, run the code again as it always works the second/third time.
		print(f"\n\nFetching protein families... (this may take several minutes)\n")
		subprocess.call("esearch -db protein -query " + query + " | efetch -format gp | grep -A 1 \"SIMILARITY\" | awk '{sub(/^[ \\t]+/,\"\");}1' > pf_raw.txt", shell=True)
		break

	elif accession_no_continue == "no":
	        # resets tax_group variable and returns to the start of the loop
		print(f"\nResetting taxonomic group...\n")
	else:
		print("\nInvalid input. Please enter yes or no.\n")




# Cleans the esearch generated file, removing everything but the name of the protein family and outputting content into a new file.

with open("pf_raw.txt", "r") as input_file, open("pf_clean.txt", "w") as output_file:
	string = " family"
	# Removes everything but the name of the protein family, superfamily and subfamily and splits all these onto different lines using '.' as a separator.
	clean_file = input_file.read().replace("\n"," ").replace(". ","\n").replace("[SIMILARITY] Belongs to the ","\n").replace("--","").replace(" belongs to the ", "").split('\n')
	# stores rows in input file which contain the string 'family'
	for row in clean_file:
		if string in row:
			print(row)
			output_file.write(row + "\n")


print("\nCounting protein families for " + tax_group + "...\n")




# Creates a pandas dataframe with separator as \n.
fam_df = pd.read_csv("pf_clean.txt")
fam_df.columns = ['Family']
# Creates a new column removing the string 'family' from all protein family names
fam_df['Protein Family'] = fam_df.apply(lambda x : x['Family'].split(' family')[0], axis=1)
family_df = fam_df[['Protein Family']]
print(str(len(fam_df)) + " out of " + str(len(accs)) + " protein entries contain protein family information.")

# Creates a new dataframe counting the number of unique protein families
protein_family_df = family_df['Protein Family'].value_counts()
print(str(len(protein_family_df)))
row_number = len(protein_family_df)
# Shows the user how many protein families are present in their taxonomic group
print("\nThere are " + str(row_number) + " protein families in " + tax_group + ".\n")




# The point of this section of code is for the user to see what different protein families there are for their chosen tax group.
# From this they can see if their desired protein family is present in the tax group.
# If so, they can use the exact name of the protein family in their query to obtain more accurate search results.

while True:
	# Asks the user if they want to see the different protein families in their taxonomic group
	pfam_table = input("Would you like to view the protein families present in " + tax_group + "? (yes/no)\n")
	if pfam_table == "yes":
		while True:
			while True:
				# Asks how many top protein families they would like to view
				table_size = input("\nEnter the number of protein families with the most proteins you would like to view. (e.g. 5, 10, 50 etc.)\n")
				# If their requested number of protein families (X) is <= the total number of protein families
				# print the top X protein families in the protein family count dataframe, exit loop
				if int(table_size) <= int(row_number):
					table_size_df = protein_family_df.head(int(table_size))
					print(table_size_df)
					break
				#If the requested number is > than the total number of protein families reset loop
				elif int(table_size) > int(row_number):
					print("There are " + str(row_number) + " protein families in " + tax_group + ". Please enter a smaller number.")
				else:
					print("Invalid input.")

			# Asks if they want to view another number of top protein families
			table2 = input("\nWould you like to view another number of protein families? (yes/no)\n")
			# If yes, reset loop
			if table2 == "yes":
				print("\n")
			# If no, exit nested loop
			elif table2 == "no":
				break
			else:
				print("Invalid input.")
		break
	# If the user does not want to see protein families, exit nested loop
	elif pfam_table == "no":
		break
	else:
		print("Invalid input.")





# In this section, the user can search for a specific protein family and obtain the number of proteins in the protein family.
# If they search for a superfamily e.g. kinase, all protein families with name kinase will be returned.

while True:
	while True:
		pfam_interest = input("\nEnter protein family of interest\n")

		# Finds the rows of the dataframe containing the protein family string and counts number of rows. This is case insensitive.
		interest = family_df.loc[family_df['Protein Family'].str.contains(pfam_interest, case=False)].value_counts()

		if len(interest) == 0:
	                # If there are no results with sepcified name, they are asked to enter a different protein family namme.
			print("\nThere are no proteins from the " + pfam_interest  + " family. Please enter a different search.\n")

		else:
			# Prints the number of protein families returned in their search
			print("\nThere are " + str(len(interest)) + " protein families containing the name " + pfam_interest + " in " + tax_group)

		while True:
			# Asks if the user wants to view the protein family names and their counts.
			interest_view = input("\nWould you like to view the protein family search results? (yes/no)\n")
			if interest_view == "yes":
				print(interest)
				break
			elif interest_view == "no":
				break
			else:
				print("Invalid input.")
		break
	# Asks if they want to search for another protein family.
	pfam_interest_2 = input("\nWould you like to search for another protein family? (yes/no)\n")
	if pfam_interest_2 == "yes":
		print("\n")
	elif pfam_interest_2 == "no":
		break
	else:
		print("Invalid input.")




# The user is now asked for the protein family they would like to query, sequences for their full query are generated.


accessionid = []
proteinfam = []

while True:
	while True:
		while True:
			prot_family = input("\nPlease enter a protein family name to query.\n")
			prot_family_correct = input("\nYou entered " + prot_family + ". Is this correct? (yes/no)\n")
			if prot_family_correct == "yes":
				print("\nQuerying taxonomic group " + tax_group + " and protein family " + prot_family + " ...\n")
				break
			elif prot_family_correct == "no":
				print(f"\nResetting protein family...\n")
			else:
				print("Invalid input.")

		# Creates a variable for the query to be used in the esearch.
		pfam_query_df = family_df.loc[family_df['Protein Family'].str.contains(prot_family, case=False)]
		new_query = f"'{tax_group}[Organism:exp] AND swissprot[filter] AND {prot_family}'"
		subprocess.call("esearch -db protein -query " + new_query + " | efetch -format acc > prot_accessions.txt", shell=True)
		print("Progress:")
		with open('prot_accessions.txt', 'r') as input_file, open('new_prot_accessions.txt', 'w') as output_file:
			prot_accessions = input_file.read().rsplit('\n')
			for accession in prot_accessions:
				p_family = subprocess.run("esearch -db protein -query " + accession + "[accession] | efetch -format gp | grep -A 1 \"SIMILARITY\" | awk '{sub(/^[ \\t]+/,\"\");}1'", shell=True, capture_output=True, text=True)
				p_return = p_family.stdout.strip()
				p_replaced = p_return.replace("\n", " ").replace(". ","\n").replace("[SIMILARITY] Belongs to the ","")
				string_f = " family"
				for row in p_replaced.split("\n"):
					if string_f in row:
						if prot_family.lower() in row.lower():
							accessionid.append(accession)
							proteinfam.append(row)
							print(str(len(accessionid)) + "/" + str(len(pfam_query_df)) + "accessions queried.")
							output_file.write(accession + "\t" + row + "\n")
				if len(accessionid) < len(pfam_query_df):
					pass
				else:
					print("\nQuery complete!\n")
					print(str(len(accessionid)) + " " +  prot_family + " family entries identified for " + tax_group)
					break

		# Shows the user the number of protein entries for their desired taxonomic group and protein family
		pfam_query_df2 = family_df.loc[family_df['Protein Family'].str.contains(prot_family, case=False)].value_counts()
		if len(pfam_query_df2) > 1:
			print("\n\nHere are the unique protein families which contain the name " + prot_family + " and their protein counts.\n")
			print(pfam_query_df2)
			# If the user decides to choose 1 protein family, this will reduce the number of sequences used in the conservation analysis.
			# This may also increase accuracy as only closely related proteins from one family are used.
			print("\nYour query contains " + str(len(pfam_query_df2)) + " protein families.\nFor more accurate protein conservation results, it is recommended to query 1 protein family.\n")
			break
		else:
			break

	accession2_no_continue = input("\nWould you like to continue with your query? (yes/no)\n")
	if accession2_no_continue == "yes":
		# Collects the protein sequences in fasta format and stores in a fasta file
		print("\nFetching protein sequences... (this may take several minutes)\n")
		with open("sequences.fa", "w") as output_file:
			for accession in accessionid:
				seqs = subprocess.run("esearch -db protein -query " + accession + "[accession] | efetch -format fasta", shell=True, capture_output=True, text=True)
				seqs_return = seqs.stdout.strip()
				output_file.write(seqs_return + "\n")
				print(accession)
		print("\nProtein sequences generated!")
		break

	elif accession2_no_continue == "no":
		print(f"\nResetting protein family...\n")

	else:
		print("\nInvalid input. Please enter yes or no.\n")

