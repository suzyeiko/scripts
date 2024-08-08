import sys
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from fasta_parser import read_fasta
from matplotlib_venn import venn2, venn3
import random

assembly_metadata_file = str(sys.argv[1])
fastafile = str(sys.argv[2])


assembly_metadata = pd.read_csv(assembly_metadata_file, sep="\t", index_col=False)

##### Genus and species #####
# Function to split the column "scientific_name" into genus and species
def split_scientific_name(name):
    parts = name.split()
    if parts[0] == "Candidatus":
        genus = " ".join(parts[:2])
        species = parts[0]+" "+parts[1]+" "+parts[2]
    else:
        genus = parts[0]
        species = parts[0]+" "+parts[1]
    return pd.Series([genus, species])

# Apply the function and create columns genus and species
assembly_metadata[['genus', 'species']] = assembly_metadata['scientific_name'].apply(split_scientific_name)
assembly_metadata

##### Sequence length #####
# Calculate sequence length
seq_length = {}
for record in read_fasta(fastafile):
   seq_length[record.accession] = len(record)

#Remove .number after accession ID
for n in range(10):
    id = "."+str(n)
    #print(id)
    seq_length = {k.replace(id,''):v for k,v in seq_length.items()}

# Convert seq_length dictionary to DataFrame
seq_length_df = pd.DataFrame(list(seq_length.items()), columns=['sequence_accession', 'contig_length'])

# Merge the DataFrame with seq_length DataFrame
assembly_metadata = pd.merge(assembly_metadata, seq_length_df, on='sequence_accession', how='left')
###########################

# Define parameters for filtering sequences based on assembly information
l50 = 2
n_contigs = 24
contiglength = 1000

def filter_seqs_based_on_assembly_info(assembly_metadata, l50, n_contigs, contiglength):
    filtered_assembly_metadata = assembly_metadata[(assembly_metadata['contig_l50'] <= l50) & (assembly_metadata['number_of_contigs'] <= n_contigs) & (assembly_metadata['contig_length'] >= contiglength)]
    return filtered_assembly_metadata

# Filter sequences
filtered_assembly_metadata = filter_seqs_based_on_assembly_info(assembly_metadata)

print("From", len(assembly_metadata), "initial sequences,", len(filtered_assembly_metadata), "were filtered based on assembly quality.")
print("Filtering criteria:")
print("Assemblies with L50 <=", l50)
print("Assemblies with number of contigs <=", n_contigs)
print("Sequences with length <=", contiglength)

# Write assembly metadata of filtered sequences
filtered_assembly_metadata.to_csv('analysis/1.plasmid_metadata/filtered_contigs_assemblystats.tsv', sep="\t", index=False)

# Write list of accession IDs
with open('analysis/1.plasmid_metadata/accession_filtered_contigs_using_assembly.txt', 'w') as f:
    f.write(filtered_assembly_metadata['sequence_accession'].str.cat(sep='\n'))


##### Write one assembly ID by species #####
seq_count_by_assembly_filtered = filtered_assembly_metadata.assembly_accession.value_counts()
filtered_seqs_by_assembly = pd.merge(seq_count_by_assembly_filtered, filtered_assembly_metadata.groupby('assembly_accession').first(), how = "inner", on=["assembly_accession"])

# Group by species
grouped = filtered_seqs_by_assembly.groupby('species')

# Calculate statistics
stats_species = grouped.agg({
    'count': ['mean', 'std', 'min', 'max', 'median', 'size']
})

# Rename the columns
stats_species.columns = ['mean', 'std', 'min', 'max', 'median', 'count']

# Create a column with list of assembly_IDs
assembly_ids = grouped.apply(lambda x: list(x.index))
stats_species['assembly_IDs'] = assembly_ids

#stats_species

random.seed(123)

with open('analysis/1.plasmid_metadata/assemblyaccession_list_by_species.txt', 'w') as f:
    for species in stats_species["assembly_IDs"]:
        #print(random.choice(species))
        f.write(random.choice(species)+"\n")
