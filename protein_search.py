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
		confirm_tax_group = input("\nYou entered " + tax_group + ". Is this correct? (yes/no)\n")
		if confirm_tax_group == "yes":
			print(f"\n\nQuerying taxonomic group " + tax_group + "...\n")
			break
		elif confirm_tax_group == "no":
			print(f"\nResetting taxonomic group...\n")
		else:
			print("\nInvalid input. Please enter yes or no.\n")

	# Searches for entries with taxnomic group with the swissprot filter and collects accession IDs
	tax_query = f"'{tax_group}[Organism:exp] AND swissprot[filter]'"
	subprocess.call("esearch -db protein -query " + tax_query + " | efetch -format acc > tax_accessions.txt", shell=True)

	# Shows user the number of protein entries for taxonomic group and asks if they want to proceed
	with open('tax_accessions.txt', 'r') as tax_accessions:
		tax_accs = tax_accessions.read().rsplit('\n')
		print("\nThere are " + str(len(tax_accs)) + " protein entries for " + tax_group)

	continue_tax_group = input(f"\nWould you like to continue querying {tax_group}? (yes/no)\n")
	if continue_tax_group == "yes":
		# collects protein family and superfamily information from the GenPept summary and removes unwated whitespace for each protein entry then exits loop
                # sometimes this esearch returns 500 internal server error. If this occurs, run the code again as it always works the second/third time.
		print(f"\n\nFetching protein families... (this may take several minutes)\n")
		subprocess.call("esearch -db protein -query " + tax_query + " | efetch -format gp | grep -A 1 \"SIMILARITY\" | awk '{sub(/^[ \\t]+/,\"\");}1' > pf_raw.txt", shell=True)
		break

	elif continue_tax_group == "no":
	        # resets tax_group variable and returns to the start of the loop
		print(f"\nResetting taxonomic group...\n")
	else:
		print("\nInvalid input. Please enter yes or no.\n")


string = " family"

# Cleans the esearch generated file, removing everything but the name of the protein family and outputting content into a new file.

with open("pf_raw.txt", "r") as input_file, open("pf_clean.txt", "w") as output_file:
	# Removes everything but the name of the protein family, superfamily and subfamily and splits all these onto different lines using '.' as a separator.
	clean_pfam = input_file.read().replace("\n"," ").replace(". ","\n").replace("[SIMILARITY] Belongs to the ","\n").replace("--","").replace(" belongs to the ", "").split('\n')
	# stores rows in input file which contain the string 'family'
	for row in clean_pfam:
		if string in row:
			print(row)
			output_file.write(row + "\n")


print("\nCounting protein families for " + tax_group + "...\n")




# Creates a pandas dataframe with separator as \n.
df = pd.read_csv("pf_clean.txt")
df.columns = ['Family']
# Creates a new column removing the string 'family' from all protein family names
df['Protein Family'] = df.apply(lambda x : x['Family'].split(' family')[0], axis=1)
family_df = df[['Protein Family']]
print(str(len(family_df)) + " out of " + str(len(tax_accs)) + " protein entries contain protein family information.")

# Creates a new dataframe counting the number of unique protein families
unique_family_df = family_df['Protein Family'].value_counts()
row_number = len(unique_family_df)
# Shows the user how many protein families are present in their taxonomic group
print("\nThere are " + str(row_number) + " protein families in " + tax_group + ".\n")




# The point of this section of code is for the user to see what different protein families there are for their chosen tax group.
# From this they can see if their desired protein family is present in the tax group.
# If so, they can use the exact name of the protein family in their query to obtain more accurate search results.

while True:
	# Asks the user if they want to see the different protein families in their taxonomic group
	df_view = input("Would you like to view the protein families present in " + tax_group + "? (yes/no)\n")
	if df_view == "yes":
		while True:
			while True:
				# Asks how many top protein families they would like to view
				df_size = input("\nEnter the number of top protein families you would like to view. (e.g. 5, 10, 50 etc.)\n")
				# If their requested number of protein families (X) is <= the total number of protein families
				# print the top X protein families in the protein family count dataframe, exit loop
				if int(df_size) <= int(row_number):
					df_view_size = unique_family_df.head(int(df_size))
					print(df_view_size)
					break
				#If the requested number is > than the total number of protein families reset loop
				elif int(df_size) > int(row_number):
					print("There are " + str(row_number) + " protein families in " + tax_group + ". Please enter a smaller number.")
				else:
					print("Invalid input.")

			# Asks if they want to view another number of top protein families
			df_view_2 = input("\nWould you like to view another number of protein families? (yes/no)\n")
			# If yes, reset loop
			if df_view_2 == "yes":
				print("\n")
			# If no, exit nested loop
			elif df_view_2 == "no":
				break
			else:
				print("Invalid input.")
		break
	# If the user does not want to see protein families, exit nested loop
	elif df_view == "no":
		break
	else:
		print("Invalid input.")





# In this section, the user can search for a specific protein family and obtain the number of proteins in the protein family.
# If they search for a superfamily e.g. kinase, all protein families with name kinase will be returned.

while True:
	pfam_df_search = input("\nWould you like to search for a specific protein family? (yes/no)\n")
	if pfam_df_search == "yes":
		while True:
			while True:
				pfam_interest = input("\nEnter protein family of interest\n")

				# Finds the rows of the dataframe containing the protein family string and counts number of rows. This is case insensitive.
				find_interest = family_df.loc[family_df['Protein Family'].str.contains(pfam_interest, case=False)].value_counts()

				if len(find_interest) == 0:
	                	# If there are no results with sepcified name, they are asked to enter a different protein family namme.
					print("\nThere are no proteins from the " + pfam_interest  + " family. Please enter a different search.\n")
				elif len(find_interest) > 0:
					# Prints the number of protein families returned in their search
					print("\nThere are " + str(len(find_interest)) + " protein families containing the name " + pfam_interest + " in " + tax_group)
					print(find_interest)
					break
				else:
					print("Something went wrong. Please try again or enter another protein family.")
			pfam_df_search_2 = input("\nWould you like to search for another protein family? (yes/no)\n")
			if pfam_df_search_2 == "yes":
				pass
			elif pfam_df_search_2 == "no":
				break
			else:
				print("Invalid input.")
		break
	elif pfam_df_search == "no":
		break
	else:
		print("Invalid input. Please enter yes/no.")

# The user is now asked for the protein family they would like to query, sequences for their full query are generated.


accessionid = []
proteinfam = []

while True:
	while True:
		while True:
			protein_family = input("\nPlease enter a protein family name to query.\n")
			confrim_protein_family = input("\nYou entered " + protein_family + ". Is this correct? (yes/no)\n")
			if confrim_protein_family == "yes":
				print("\nQuerying taxonomic group " + tax_group + " and protein family " + protein_family + " ...\n")
				break
			elif confrim_protein_family == "no":
				print(f"\nResetting protein family...\n")
			else:
				print("Invalid input.")

		# Creates a variable for the query to be used in the esearch.
		protein_family_df = family_df.loc[family_df['Protein Family'].str.contains(protein_family, case=False)]
		query = f"'{tax_group}[Organism:exp] AND swissprot[filter] AND {protein_family}'"
		subprocess.call("esearch -db protein -query " + query + " | efetch -format acc > prot_accessions.txt", shell=True)
		with open('prot_accessions.txt', 'r') as input_file2, open('new_prot_accessions.txt', 'w') as output_file2:
			prot_accessions = input_file2.read().rsplit('\n')
			if len(prot_accessions) <= 1:
				print("No proteins with name " + protein_family + " were found.")
				break
			elif len(prot_accessions) > 2:
				for accession in prot_accessions:
					acc_esearch = subprocess.run("esearch -db protein -query " + accession + "[accession] | efetch -format gp | grep -A 1 \"SIMILARITY\" | awk '{sub(/^[ \\t]+/,\"\");}1'", shell=True, capture_output=True, text=True)
					esearch_result = acc_esearch.stdout.strip()
					clean_result = esearch_result.replace("\n", " ").replace(". ","\n").replace("[SIMILARITY] Belongs to the ","")
					for row in clean_result.split("\n"):
						if string in row:
							if protein_family.lower() in row.lower():
								accessionid.append(accession)
								proteinfam.append(row)
								print(protein_family + ": " +str(len(accessionid)) + "/" + str(len(protein_family_df)) +  " identified.")
								output_file2.write(accession + "\t" + row + "\n")
					if len(accessionid) < len(protein_family_df):
						pass
					else:
						print("\nAll proteins found!\n")
						print(str(len(accessionid)) + " " +  protein_family + " family entries identified for " + tax_group)
						break
				break
			else:
				print("Something went wrong. Please try again.")
		# Shows the user the number of protein entries for their desired taxonomic group and protein family
		protein_family_df_counts = family_df.loc[family_df['Protein Family'].str.contains(protein_family, case=False)].value_counts()
		if len(protein_family_df_counts) > 1:
			print("\n\nHere are the unique protein families which contain the name " + protein_family + " and their protein counts.\n")
			print(protein_family_df_counts)
			# If the user decides to choose 1 protein family, this will reduce the number of sequences used in the conservation analysis.
			# This may also increase accuracy as only closely related proteins from one family are used.
			print("\nYour query contains " + str(len(protein_family_df_counts)) + " protein families.\nFor more accurate protein conservation results, it is recommended to query 1 protein family.\n")
			break
		else:
			break

	continue_protein_family = input("\nWould you like to continue with your query? (yes/no)\n")
	if continue_protein_family == "yes":
		# Collects the protein sequences in fasta format and stores in a fasta file
		print("\nFetching protein sequences... (this may take several minutes)\n")
		with open("sequences.fa", "w") as output_file:
			for accession in accessionid:
				seq_esearch = subprocess.run("esearch -db protein -query " + accession + "[accession] | efetch -format fasta", shell=True, capture_output=True, text=True)
				seq_result = seq_esearch.stdout.strip()
				output_file.write(seq_result + "\n")
				print(accession)
		print("\nProtein sequences generated! - see 'sequences.fa'.")
		break

	elif continue_protein_family == "no":
		print(f"\nResetting protein family...\n")

	else:
		print("\nInvalid input. Please enter yes or no.\n")

