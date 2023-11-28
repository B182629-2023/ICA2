#!/usr/bin/python3

import subprocess
import pandas as pd
import os
import sys


print("\n\nThis is the second part of the sequence conservation analysis program.")
print("Here, the sequences in the FASTA file generated from the Entrez esearch query is aligned.")
print("Conservation statistics will be generated for each protein sequence of interest.\n\n")



# Asks the user if they want to run the program, if not, exit
while True:
        cont = input("Do you wish to continue? (yes/no)\n")
        if cont == "yes":
                print("\n\nStarting program...\n")
                break
        if cont == "no":
                print("\nExiting program...")
                sys.exit()
        else:
                print("Invalid input.") 


# Enures the user is still in the analysis directory
os.chdir(os.environ['HOME'])
os.chdir("conservation_analysis")


# Runs a multiple sequence alingment of chosen sequence set in clustal format
subprocess.call("clustalo -i sequences.fa -o aligned_seqs.msf --outfmt=msf", shell=True)
print("\nClustering protein sequences...\n")


# Uses multiple sequence alignment file to create a conservation plot and outputs it to screen
subprocess.call("plotcon -sequence aligned_seqs.msf -graph png -win 10 -goutfile conservation_plot -scorefile EBLOSUM62", shell=True)
print("\nGenerating conservation plot...\n")
subprocess.call("eog *.png", shell=True)


# Uses infoalign to generate multiple sequence alignment data.
subprocess.call("infoalign aligned_seqs.msf -nousa -noname -noweight -nodescription -outfile conservation_score.infoalign", shell=True)
subprocess.call("infoalign aligned_seqs.msf -only -name -outfile conservation_score_name.infoalign", shell=True)


# Create a dataframe to visualise multiple sequence alignment statistics data. Original DF has one big column for USA, NAME and SEQLENGTH.
# To fix this, dataframe containing only identifier for each seq 'names' is created. 
name_df = pd.read_csv("conservation_score_name.infoalign", sep="_")
# To extract Accession ID only, use separator delimiter as _ and extract accession ID
name_df.columns = ['SP', 'Accession ID', 'A', 'B', 'C']
name_df_new = name_df[['Accession ID']]
# Create dataframe with all other columns and join both dataframes together horizontally
score_df = pd.read_csv("conservation_score.infoalign", sep="\t")
aln_score_df = pd.concat([name_df_new, score_df], axis=1)
select_aln_score_df = aln_score_df.iloc[:, 0:9].dropna() # Random 'Unamed' column containing only NaNs removed using iloc to select the columns we want




# Asks the user if they would like to view the alignment statistics
while True:
	view_score = input("\nWould you like to view the sequence conservation statistics for these sequences? (yes/no)\n") 
	if view_score == "yes":
		print(select_aln_score_df)
		break
	elif view_score == "no":
		break
	else:
		print("Invalid Input.")




# Shows the user the distance matrix sores for aligned sequence set

# Break out flag created to exit outer loop when flag=False
break_out_flag1 = True

while break_out_flag1:

	aa_plot = input("\nWould you like to view a distance matrix for the aligned sequences (yes/no)\n")
	if aa_plot == "yes":
		while True:
			# asks user for their preffered scoring method
			protmethod = input("\nWhich scoring matrix method would you like to use?\n(0 = Uncorrected; 1 = Jukes-Cantor; 2 = Kimura Protein)\n")
			if int(protmethod) == 0:
				subprocess.call("distmat aligned_seqs.msf -outfile aligned_seqs.distmat -protmethod 0", shell=True)
				break
			elif int(protmethod) == 1:
				subprocess.call("distmat aligned_seqs.msf -outfile aligned_seqs.distmat -protmethod 1", shell=True)
				break
			elif int(protmethod) == 2:
				subprocess.call("distmat aligned_seqs.msf -outfile aligned_seqs.distmat -protmethod 2", shell=True)
				break
			else:
				print("Invalid input.")

		with open("aligned_seqs.distmat", "r") as matrix:
			print(matrix.read())

		while True:
			plot_seq_2 = input("\nWould you like to view the matrix with another matrix scoring method? (yes/no)\n")
			if plot_seq_2 == "yes":
				# removes existing file to avoid errors/incorrect overwriting
				os.remove("aligned_seqs.distmat")
				# Breaks out of current loop and reses to top of the outer loop
				break
			if plot_seq_2 == "no":
				# Breaks out of outer loop
				break_out_flag1 = False
				break
			else:
				print("Invalid input.")
	elif aa_plot == "no":
		break
	else:
		print("Invalid input.")


# Imports PROSITE database  (I realised there is a way to scan sequences against PROSITE using an emboss program but it refused to work! Had to resort to importing the database... sorry!)
print("\n Importing the PROSITE database...\n")
subprocess.call("wget ftp://ftp.expasy.org/databases/prosite/prosite.dat", shell=True)


# Creates a function which splits all entries in the prosite database and creates a separate dictionaries for all entries
def split_entry(entry):
	rows = entry.split("\n")
	prosite_entry = {}
	for row in rows:
		# Checks that each row has 2 values in order to create a key and value pair
		linesplit = row.split(None, 1)
		if len(linesplit) > 1:
			# Creates key and value pairs
			key, value = linesplit
			# As DR key contains values running over multiple lines, this checks if a key for DR already exists when iterating over the rows
			if key == "DR": 
				if key in prosite_entry:
					# If a DR key exists, the value (i.e. row) is appended to the DR key in the relevant prosite_entry dictionary
					prosite_entry[key].append(value.strip())
				else:
					# If a DR key does not yet exist, it creates a key with the respecive value (i.e. row)
					prosite_entry[key] = [value.strip()]
			else:
				# If the row does not contain DR key, continue as normal with other keys
				prosite_entry[key] = value.strip()
	return prosite_entry





# This dictionary contains all prosite_entry sub-dictionaries with the PROSITE accession as the key of each prosite_entry dictionary.
prosite_dict = {}

with open("prosite.dat", "r") as prosite_file:
	# The first 26 lines are introductory lines, so they are skipped
	content = prosite_file.read()[26:]
	# Each entry is split by // characters
	entries = content.split("//\n")
	for entry in entries:
		rows = entry.split("\n")
		entry_ac = None

		# extract the AC (prosite accession) in entry
		for row in rows:
			if "AC   " in row:
				entry_ac = row.split()[1]
		# call function to create a nested dictionary with prosite accession as key and prosite_entry dictionaries as values.
		if entry_ac:
			prosite_dict[entry_ac] = split_entry(entry)



print("\nScanning PROSITE database for matching motifs...\n")

# Extract chosen protein sequence set accessions to query against prosite dictionary
with open("new_prot_accessions.txt", "r") as prot_accessions, open("protein_accessions.txt", "w") as new_acc_file:
	content = prot_accessions.read().split("\n")
	for row in content:
		# extracts the accession ID only and then gets rid of the version ID
		versions = row.split("\t")[0]
		accessions = versions.split(".")[0]
		new_acc_file.write("\n" + accessions)


# Create empty lists to append values in for loop
prositeid = []
acc = []
prositedesc = []

print("\nFetching matches...\n")

# For each accession, scans the dictionary for a match. If found, extracts prosite ID and motif description and appends to list. 
with open("protein_accessions.txt", "r") as prot_accs:
	content_2 = prot_accs.read().strip()
	accs = content_2.split("\n")
	for accession_2 in accs:
		for key, value in prosite_dict.items():
			# Finds DR values in each prosite entry directory in nested dictionary AND finds accession string in DR values.
			if 'DR' in value and any(accession_2 in DR_value for DR_value in value['DR']):
				# Appends prosite accession ID
				prositeid.append(value.get('AC'))
				# Appends protein accession ID
				acc.append(accession_2)
				# Appends motif description
				prositedesc.append(value.get('DE'))


print("\nSearch complete!\n")





# Create a dataset and a pandas dataframe from prosite data matches.
# Uses lists with values appended in for loop to create DF
prosite_data = {'Protein Count': acc, 'Prosite Accession' : prositeid, 'Motif': prositedesc}
prosite_df = pd.DataFrame(prosite_data)

# Counts the unique prosite accessions found in the protein sequence set and sorts motifs in order of most to least number of sequences containing motif
prosite_df_counts = prosite_df.groupby(['Prosite Accession', 'Motif'])['Protein Count'].nunique().reset_index()
unique_motifs = len(prosite_df_counts)
prosite_df_counts_desc = prosite_df_counts.sort_values('Protein Count', ascending=False)
print(str(unique_motifs) + " prosite motifs identified for chosen sequence set.\n")




# Asks the user how many prosite motifs found in sequence set they want to view

break_out_flag_3 = True

while break_out_flag_3:
	while True:
		view_prosite_df = input("\nHow many top PROSITE motifs would you like to view? (e.g. 5, 10, 50 etc.)\n")
		# If the requested number is <=  the total number of unique motifs, print dataset.
		if int(view_prosite_df) <= int(unique_motifs):
			prosite_df_count = prosite_df_counts_desc.head(int(view_prosite_df))
			print(prosite_df_count)
			break
		# If requested number is smaller, ask to enter a smaller number.
		elif int(view_prosite_df) > int(unique_motifs):
			print(str(unique_motifs) + " prosite motifs identified for chosen sequence set. Please enter a smaller number.\n")
		else:
			print("Invalid input")
	# Asks if  they wish to view another number
	while True:
		view_prosite_data2 = input("\nWould you like to view another number of prosite motifs? (yes/no)\n")
		if view_prosite_data2 == "yes":
			break
		elif view_prosite_data2 == "no":
			break_out_flag_3 = False
			break
		else:
			print("Invalid input.")



print("\nProgram complete!\n")


