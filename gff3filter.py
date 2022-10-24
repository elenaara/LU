# -*- coding: utf-8 -*-
"""
Created on Mon Oct 24 10:20:30 2022

@author: Elena Aramendía

Transcript filter, the script takes a gff3 file and removes transcripts 
according to the following filters:
    - Non protein coding genes/transcripts: rows with gene_type or 
    transcript_type != protein_coding are removed
    - Transcript support level: transcript with bad support (>= 4) or "NA" are removed
    - Transcripts with the following tags are removed:
        "cds_start_NF","cds_end_NF","fragmented_locus",
        "inferred_exon_combination","low_sequence_quality",
        "non_canonical_genome_sequence_error","non_canonical_TEC",
        "not_best_in_genome_evidence","not_organism_supported"
    - Mithocondrial genes: gene names starting with MT- o MTRNR2
    - Histone genes: gene names starting with H2A,H2B,H3 or H4
    
The returns a text file with the unique transcript IDs.

Usage:
    gff3filter.py -i input_file -o1 [output_table] -o2 [output_tr]

Filtering based on
        Zaitsev, Aleksandr et al. 
        “Precise reconstruction of the TME using bulk RNA-seq and a machine learning algorithm trained on artificial transcriptomes.” 
        Cancer cell vol. 40,8 (2022): 879-894.e16. doi:10.1016/j.ccell.2022.07.006
"""

#%% Modules
import argparse 
#%% Argument parser
parser = argparse.ArgumentParser()
usage = 'This program takes a gff3 file, filters it and generates a new gff3 file with info on last column re-formatted as columns for each field and text file with transcript ids of the transcripts kept'
parser.add_argument(                                                  # distance matrix
    '-i',
    type = argparse.FileType('r'),
    dest = 'gff_input',
    required = True,
    help = 'Input gff3 file'
    )

parser.add_argument(                                                  # output file
    '-o1',
    dest = 'output_table',
    type = argparse.FileType('w'),
    required = False,
    default = 'output_table',
    help = 'Output gff3 table, with info on last column re-formatted as columns for each field'
    )

parser.add_argument(                                                  # output file
    '-o2',
    dest = 'output_tr',
    type = argparse.FileType('w'),
    required = False,
    default = 'output_tr',
    help = 'Output text file with the unique transcript IDs '
    )


args = parser.parse_args() 
#%% Main code
gff = args.gff_input
out_tr = args.output_tr
out_table = args.output_table

#%%
# Test files
# gff = open("C:/Users/Admin/Desktop/MT_2223/RNAseq/ref/gencode.v39.annotation.gff3", 'r') 
# out_table = open("C:/Users/Admin/Desktop/MT_2223/RNAseq/ref/test_table.tsv", 'w') 
# out_tr = open("C:/Users/Admin/Desktop/MT_2223/RNAseq/ref/test_tr_to_keep.txt", 'w') 

# Get fields in the table
header_list = list()


# Dictionary storing everything
tr_dict = dict()
# Read through the gff3 file
# While reading filter for everything
# We want to keep:
    # gene_type="protein_coding"
    # transcript_type="protein_coding"
biotypes = ("IG_V_gene", "IG_D_gene", "IG_J_gene", "TR_V_gene", "TR_D_gene", "TR_J_gene","pseudogene", "polymorphic_pseudogene", "processed_pseudogene", "transcribed_processed_pseudogene", "unitary_pseudogene", "unprocessed_pseudogene", "transcribed_unitary_pseudogene", "IG_C_pseudogene", "IG_J_pseudogene", "IG_V_pseudogene", "transcribed_unprocessed_pseudogene", "translated_unprocessed_pseudogene", "TR_J_pseudogene", "TR_V_pseudogene", "snRNA", "snoRNA", "miRNA", "ribozyme", "rRNA", "Mt_tRNA", "Mt_rRNA", "scaRNA", "retained_intron", "sense_intronic", "sense_overlapping", "nonsense_mediated_decay", "non_stop_decay", "antisense", "lincRNA", "macro_lncRNA", "processed_transcript", "3prime_overlapping_ncrna", "sRNA", "misc_RNA", "vaultRNA", "TEC")
# And remove:
    # Mit genes (start with MTRNR2)
    # trs with TSL<4 or NA
    # these tags
tags_to_remove = ("cds_start_NF","cds_end_NF","fragmented_locus","inferred_exon_combination","low_sequence_quality","non_canonical_genome_sequence_error","non_canonical_TEC","not_best_in_genome_evidence","not_organism_supported")


for line in gff: # read each line
    # Boolean for filtering
    add = True
    if not line.startswith('#'): # Find info lines
    
        # get last column
        line = line.split('\t')
        seqid = line[0]
        source = line[1]
        typ = line[2]
        start = line[3]
        end = line[4]
        score = line[5]
        strand = line[6]
        phase = line[7]
        
        data = line[8]
        data = data.strip()
        
        # Split fields
        data = data.split(';')
        for field in data: # for each field in that column
            if add == True:
                field = field.split("=") # split 
                header = field[0]
                value = field[1]
                
                if header == "ID": # if it is first entry for this transcript
                    ID = value
                    # create sub dict and add ID and all first columns
                    tr_dict[ID] = {header: ID}
                    tr_dict[ID]["seqid"] = seqid
                    tr_dict[ID]["source"] = source
                    tr_dict[ID]["type"] = typ
                    tr_dict[ID]["start"] = start
                    tr_dict[ID]["end"] = end
                    tr_dict[ID]["score"] = score
                    tr_dict[ID]["strand"] = strand
                    tr_dict[ID]["phase"] = phase
                    
                # if the field name is not stored, add it to the list
                if header not in header_list:
                    header_list.append(header)
    
                # if it not is the first entry
                if ID in tr_dict:
                    # Add next entry
                    tr_dict[ID][header] = value
                    
                 # FILTER 
                 # Keep protein coding
                 
                if header == "transcript_type":
                    if value != "protein_coding":
                        add = False
                    # if value in biotypes:
                    #     add = False
                  # Keep TSL > 4
                elif header == "transcript_support_level":
                    if value == "NA":
                        add = False
                    else:
                        if int(value) >= 4:
                            add = False
                         
                  # Remove tags
                elif header == "tag":
                    tags = value.split(",")
                    for tag in tags:                    
                        if tag in tags_to_remove:
                            add = False  
                            break
                        
                # Remove MT genes
                elif header == "gene_name":
                      if value.startswith("MTRNR2") or value.startswith("MT-"):
                        add = False
                        
                # Remove histone genes
                elif header == "gene_name":
                      if value.startswith("H2A") or value.startswith("H2B") or value.startswith("H3") or value.startswith("H4"):
                        add = False
                
                if ID in tr_dict and add == False:
                    del tr_dict[ID]
                
                

# Print out transcript to keep
# Set to avoid repeating transcripts
printed_trs = set()
for tr in tr_dict:
    for field in header_list: 
        if field in tr_dict[tr] and field == "transcript_id":
            # Search for transcript field and print tr_id
            if tr_dict[tr][field] not in printed_trs:    
                printed_trs.add(tr_dict[tr][field])
                print(tr_dict[tr][field], file = out_tr)

#%%# Now that we have everything, create table

# Print header
gff3_fields = ("seqid","source","type","start","end","score","strand","phase")
full_header = gff3_fields + tuple(header_list)
header_print = ""
for field in full_header:
    if field == full_header[0]:
        header_print += field
    else:
        header_print += "\t" + field
        
print(header_print,file=out_table)

# Print all the rows
# for each transcript
for tr in tr_dict:
    row_print = ""
    # For each column
    for field in full_header: 
        # if there is a value for the column print that
        if field in tr_dict[tr]:
            if field == full_header[0]: # first value
                row_print += tr_dict[tr][field]
            else:
                row_print += "\t" + tr_dict[tr][field] 
                
        elif field not in tr_dict[tr]: 
        # If there is no value for that column print .
            row_print += "\t."
        # First field (ID) should always have a value      
            
    print(row_print,file=out_table)
