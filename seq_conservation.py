#!/usr/bin/python3

import subprocess
import pandas as pd

# Imports PROSITE database
#subprocess.call("wget ftp://ftp.expasy.org/databases/prosite/prosite.dat", shell=True)



# Runs a multiple sequence alingment of chosen sequence set in clustal format
subprocess.call("clustalo -i sequences.fa -o aligned_seqs.msf --outfmt=msf", shell=True)
print("\nClustering protein sequences...\n")


# Uses multiple sequence alignment file to create a conservation plot and outputs it to screen
#subprocess.call("plotcon -sequence aligned_seqs.msf -graph png -win 10 -goutfile conservation_plot -scorefile EBLOSUM62", shell=True)
#subprocess.call("eog *.png", shell=True)


# Uses infoalign to generate multiple sequence alignment data.
subprocess.call("infoalign aligned_seqs.msf -outfile conservation_score.infoalign", shell=True)


# Create a datafile to visualise multiple sequence alignment data
score_df = pd.read_csv("conservation_score.infoalign", sep="\t")



# Asks the user if they would like to view the alignment output
while True:
	view_score = input("Would you like to view the level of conservation between these sequences? (yes/no)\n") 
	if view_score == "yes":
		print(score_df)
		break
	elif view_score == "no":
		break
	else:
		print("Invalid Input.")

subprocess.call("fprotdist aligned_seqs.msf -outfile distance_comparisons", shell=True)

subprocess.call("fdrawtree 


with open("distance_comparisons.txt", "r") as distance_file:
	print(distance_file.read())

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





# Extract chosen protein sequence set accessions to query against prosite dictionary
with open("new_prot_accessions.txt", "r") as prot_accessions, open("protein_accessions.txt", "w") as new_acc_file:
	content = prot_accessions.read().split("\n")
	for row in content:
		versions = row.split("\t")[0]
		accessions = versions.split(".")[0]
		new_acc_file.write("\n" + accessions)



prositeid = []
acc = []
prositedesc = []


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


print("Query complete!")





# Create a dataset and a pandas dataframe from prosite data matches.
prosite_data = {'Protein Accession': acc, 'Prosite Accession': prositeid, 'Motif': prositedesc}
prosite_df = pd.DataFrame(prosite_data)

# Counts the unique prosite accessions found in the protein sequence set and sorts motifs in order of most to least number of sequences containing motif
prosite_df_counts = prosite_df.groupby(['Prosite Accession', 'Motif'])['Protein Accession'].nunique().reset_index()
unique_motifs = len(prosite_df_counts)
prosite_df_counts_desc = prosite_df_counts.sort_values('Protein Accession', ascending=False)
print(str(unique_motifs) + " prosite motifs identified for chosen sequence set.")




# Asks the user how many prosite motifs found in sequence set they want to view
while True:
	while True:
		view_prosite_df = input("How many top prosite motifs would you like to view? (e.g. 5, 10, 50 etc.)\n")
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
	view_prosite_data2 = input("Would you like to view another number of prosite motifs? (yes/no)\n")
	if view_prosite_data2 == "yes":
		pass
	elif view_prosite_data2 == "no":
		break
	else:
		print("Invalid input.")



