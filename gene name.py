## Written by Susan Xu Summer 2021; edited by Anna October 2024
## Takes the RAREdar_REsults.txt and pulls down the common name.
## The final file (gene_names_with_counts.txt) has three tab-delimited columns:
## <RefSeq ID> <# RAREs> <Gene Name>

#import required packages
from Bio import Entrez
import pandas as pd
import time
import ssl
import sys

def get_gene_name_from_refseq(accs):
    '''
    Input: a list of refseq accession numbers.
    Returns: a list of the gene names, one for each accession number.

    Queries the Entrez database and gets a batch of names.
    '''
    # note this uses Susan Xu's reed email.
    Entrez.email = "xususan@reed.edu" # Always tell NCBI who you are. this is required for accessing their server. 
    ssl_context = ssl._create_unverified_context()
    handle = Entrez.efetch(db="nucleotide", id=','.join(accs), retmode="xml", context=ssl_context)
    records = Entrez.read(handle)
    return [records[i]['GBSeq_definition'] for i in range(len(records))]

# Load your data
data = pd.read_csv('RAREdar_Results.txt', sep='\t')  # adjust the file name and separator to match your file

# Extract the column with RefSeq IDs. Replace 'column_name' with the name of your column
raw_data = data['Gene']

# Split the strings by underscore and keep only the last two parts (the identifier)
refseq_ids = ["_".join(s.split('_')[-2:]) for s in raw_data]
print('%d refseq ids' % (len(refseq_ids)))

# The refseq ids contain duplicates. Make an id_dict
# that stores the number of times each gene is present in the
# original file (one line per RARE instance). 
id_dict = {} # (refseqID, # of times occurs) pairs
for x in refseq_ids:
    if x not in id_dict:
        id_dict[x] = 0
    id_dict[x]+=1
print('%d unique refseq ids' % (len(id_dict)))

# Initialize a dictionary to hold RefSeq ID - gene name pairs
out = open('gene_names_with_counts.txt','w')
out.write('RefSeqID\t# RAREs\tGene Name\n')
accs = []
for i, refseq_id in enumerate(id_dict.keys()):
    accs.append(refseq_id)
    # when the list of accession numbers has 500 ids or is at the end of the full list,
    # map the IDs to gene names.
    if i!= 0 and (i % 500==0 or i == len(id_dict.keys())-1):
        print('#### %d of %d (%.2f)' % (i,len(id_dict.keys()),i/len(id_dict.keys())))
        try:
            gene_names = get_gene_name_from_refseq(accs)

            # there should be exactly one name per accession.
            assert len(gene_names) == len(accs) 
            
            # write the names/ids to a file and print them to screen.
            for j in range(len(accs)): 
                print(f"RefSeq ID: {accs[j]}, Gene Name: {gene_names[j]}")
                out.write('%s\t%d\t%s\n' % (accs[j],id_dict[accs[j]],gene_names[j]))
            out.flush()
            time.sleep(1)  # a pause of half a second between requests to avoid overloading NCBI servers
        except Exception as e:
            print(f"Failed to get gene name for RefSeq ID {refseq_id}: {e}")
            continue
        
        # reset the accession IDs list for the next batch.
        accs = [] 
out.close()