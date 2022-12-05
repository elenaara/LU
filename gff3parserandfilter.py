# -*- coding: utf-8 -*-
"""
Created on Thu Nov 17 14:23:55 2022


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
    
The returns two transcript tables, one with all the transcripts and one with just
the filtered transcripts.

The table includes:
    
    - seqid (chr)
    - star position
    - end position
    - transcript id
    - gene id
    - gene type
    - gene name (symbol)
    - hgnc id
    - havana name
    - havana transcript
But all the tags are saved when the file is read so more fields could be printed 
including them in the "header" list in the printing part of the script.

"""

#%%

gff = open("D:/ref/gencode.v39.annotation.gff3", 'r') 
out_table1 = open("D:/ref/gencode39_transcript_table.tsv", 'w') 
out_table2 = open("D:/ref/gencode39_FILTEREDtranscript_table.tsv", 'w') 


## Make 2 files, one with everything and one filetered


#%% Unfiltered table

# Get fields in the table
header_list = list()


# Dictionary storing everything
tr_dict = dict()


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
                 
                if typ != "transcript":
                    add = False
                
                if ID in tr_dict and add == False:
                    del tr_dict[ID]
gff.close()
#%%# PRINT UNFILTERED TABLE

# Print header
header = ("seqid","start","end",'transcript_id','gene_id','gene_type','gene_name',
               'transcript_type','transcript_name','hgnc_id','havana_gene',
               'havana_transcript')
header_print = ""
for field in header:
    if field == header[0]:
        header_print += field
    else:
        header_print += "\t" + field
        
print(header_print,file=out_table1)

# Print all the rows
# for each transcript
for tr in tr_dict:
    row_print = ""
    # For each column
    for field in header: 
        # if there is a value for the column print that
        if field in tr_dict[tr]:
            if field == header[0]: # first value
                row_print += tr_dict[tr][field]
            else:
                row_print += "\t" + tr_dict[tr][field] 
                
        elif field not in tr_dict[tr]: 
        # If there is no value for that column print .
            row_print += "\t."
        # First field (ID) should always have a value      
            
    print(row_print,file=out_table1)


#%% Filtered table
gff = open("D:/ref/gencode.v39.annotation.gff3", 'r')
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
                    
                
                if typ != "transcript":
                    add = False
                elif typ == "transcript":
                    # FILTER
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


#%% PRINT FILTERED

# Print header
header = ("seqid","start","end",'transcript_id','gene_id','gene_type','gene_name',
               'transcript_type','transcript_name','hgnc_id','havana_gene',
               'havana_transcript')
header_print = ""
for field in header:
    if field == header[0]:
        header_print += field
    else:
        header_print += "\t" + field
        
print(header_print,file=out_table2)

# Print all the rows
# for each transcript
for tr in tr_dict:
    row_print = ""
    # For each column
    for field in header: 
        # if there is a value for the column print that
        if field in tr_dict[tr]:
            if field == header[0]: # first value
                row_print += tr_dict[tr][field]
            else:
                row_print += "\t" + tr_dict[tr][field] 
                
        elif field not in tr_dict[tr]: 
        # If there is no value for that column print .
            row_print += "\t."
        # First field (transcript_id) should always have a value      
            
    print(row_print,file=out_table2)

# Close files
gff.close()
out_table1.close()
out_table2.close()