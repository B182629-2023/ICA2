#!/usr/bin/python3

import os
import sys
import pandas as pd

import subprocess

tax_group = input("tax")



while True:
        while True:
                while True:
                        prot_family = input("\nPlease enter a protein family name to query.\n")
                        prot_family_correct = input("\nYou entered " + prot_family + ". Is this correct? (yes/no)\n")
                        if prot_family_correct == "yes":
                                print("\nQuerying taxonomic group " + tax_group + " and protein family " + prot_family + " ...\n")
                                new_query = f"'{tax_group}[Organism:exp] AND swissprot[filter] AND {prot_family}'"
                                pfam_query_df = family_df.loc[family_df['Protein Family'].str.contains(prot_family, case=False)]
                                # Collects accession IDs for protein family and tax group
                                subprocess.call("esearch -db protein -query " + new_query + " | efetch -format acc >
 prot_accessions.txt", shell=True)
                                # Reads each accession in accession file and conducts an esearch
                                with open('prot_accessions.txt', 'r') as input_file, open('new_prot_accessions.txt', 'w') as output_file:
                                        prot_accessions = input_file.read().rsplit('\n')
                                        for accession in prot_accessions:
                                                p_family = subprocess.run("esearch -db protein -query " + accession + "[accession]
 | efetch -format gp | grep -A 1 \"SIMILARITY\" | awk '{sub(/^[ \\t]+/,\"\");}1'", shell=True, capture_output=True, text=True)
                                                p_return = p_family.stdout.strip()
                                                p_replaced = p_return.replace("\n", " ").replace(". ","\n")
.replace("[SIMILARITY] Belongs to the ","")
                                                string_f = " family"
                                                for row in p_replaced.split("\n"):
                                                        if string_f in row:
                                                                if prot_family.lower() in row.lower():
                                                                        output_file.write(accession + "\t" + row + "\n")
                                                                        print(accession + "\t" + row + "\n")
                                                                        if len(output_file) == len(pfam_query_df): 
                                                                                break
                                                                        else:
                                                                                pass

                                print("done")
                                break
                        elif prot_family_correct == "no":
                                print(f"\nResetting protein family...\n")
                        else:
                                print("Invalid input.")


                # Shows the user the number of protein entries for their desired taxonomic group and protein family
                prot_accession_id = open('new_prot_accessions.txt').read().rsplit("\n")
                print("\nThere are " + str(len(prot_accession_id)) + " " +  prot_family + " family entries for " + tax_group)
                pfam_query_df = family_df.loc[family_df['Protein Family'].str.contains(prot_family, case=False)].value_counts()
                if len(pfam_query_df) > 1:
                        print("Here are the unique protein families which contain the name " + prot_family + " and their protein counts.\n")
                        print(pfam_query_df)
                        # If the user decides to choose 1 protein family, this will reduce the number of sequences used in the conservation analysis.
                        # This may also increase accuracy as only closely related proteins from one family are used.
                        print("Your query contains " + len(pfam_query_df) + " protein families. For more accurate protein conservation results, it is recommended to query 1 protein family.\n")
                        break
                else:
                        break

        accession2_no_continue = input(f"\nWould you like to continue with your query? (yes/no)\n")
        if accession2_no_continue == "yes":

                # Collects the protein sequences in fasta format and stores in a fasta file
                print("\nFetching protein sequences... (this may take several minutes)\n")
                os.system("esearch -db protein -query " + new_query + " | efetch -format fasta > sequences.fa")
                print("\nProtein sequences generated!")
                break

        elif accession2_no_continue == "no":

                print(f"\nResetting protein family...\n")

        else:
                print("\nInvalid input. Please enter yes or no.\n")

