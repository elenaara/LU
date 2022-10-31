# Script to get number of reads in each sample from kallisto/salmon json files

library(stringr)

## Kallisto version

# Get all the sample names/file paths (theyre in 4 different directories)
setwd("C:/Users/Admin/Desktop/MT_2223/RNAseq/Kallisto_GRCh38_Gencode39/")
samples1 <- list.files(path="KallistoOut_Gencode39/")
samples1 <- samples1[1:(length(samples1)-4)]
samples1_path <- paste0("KallistoOut_Gencode39/",samples1)

samples2 <- list.files(path="KallistoOut_Gencode39_batch2/")
samples2 <- samples2[1:(length(samples2)-4)]
samples2_path <- paste0("KallistoOut_Gencode39_batch2/",samples2)

samples3 <- list.files(path="KallistoOut_Gencode39_CTG/")
samples3 <- samples3[1:(length(samples3)-4)]
samples3_path <- paste0("KallistoOut_Gencode39_CTG/",samples3)

samples4 <- list.files(path="KallistoOut_Gencode39_CTG_168_026_064/")
samples4 <- samples4[1:(length(samples4)-6)]
samples4_path <- paste0("KallistoOut_Gencode39_CTG_168_026_064/",samples4)

# get all sample names together
samples <- c(samples1,samples2,samples3,samples4)
samples_path <- c(samples1_path,samples2_path,samples3_path,samples4_path)

# Empty dataframe
reads.df <- data.frame(n_processed=rep(0,length(samples)),n_pseudoaligned=rep(0,length(samples)),row.names = samples)

# File paths for the json files
files <- paste0("./", samples_path, "/run_info.json")
# Go through each file
for (x in 1:length(files)) {
  # Get sample name and open the json file
  sample <- files[x]
  sample_name <- samples[x]
  open_file <- file(sample,open="r")
  jsonfile <- readLines(open_file)
  
  # For each line in the file 
  for (i in 1:length(jsonfile)) {
    line = jsonfile[i]
    
    # Search for the n_processed line
    # Get number of reads and store in dataframe
    if (grepl("n_processed", line) == TRUE) {
      line = str_split(line,":")
      n_proc = line[[1]][1]
      n_proc = substr(n_proc,3,nchar(n_proc)-1)
      reads = line[[1]][2]
      reads <- str_remove_all(reads, ",")
      reads <- str_remove_all(reads, "\\s")
      
      reads.df[sample_name,"n_processed"] <- reads
    }
    
    # Search for the n_pseudoaligned line
    # Get number of reads and store in dataframe
    else if (grepl("n_pseudoaligned", line) == TRUE) {
      line = str_split(line,":")
      n_pseu = line[[1]][1]
      n_pseu = substr(n_pseu,3,nchar(n_pseu)-1)
      reads = line[[1]][2]
      reads <- str_remove_all(reads, ",")
      reads <- str_remove_all(reads, "\\s")
        
      reads.df[sample_name,"n_pseudoaligned"] <- reads
    }
  }
  # Close open file
  close(open_file)
}

# Save to table
write.table(reads.df,
            file = "n_processed_kallisto_uroscanseq.txt",
            sep="\t",quote = FALSE)


## Salmon version
# Get all the sample names/file paths 
setwd("C:/Users/Admin/Desktop/MT_2223/RNAseq/UROSCANSEQ_salmon_Gencode39_Elena/")
samples <- list.files(path=".")
samples <- samples[1:(length(samples)-1)]
files <- paste0("./", samples, "/lib_format_counts.json")
samples <- substr(samples,1,nchar(samples)-6)
names(files) <- samples


#num_compatible_fragments

# Empty dataframe
reads.df <- data.frame(n_compatible_fragments=rep(0,length(samples)),row.names = samples)

# Go through each file
for (x in 1:length(files)) {
  # Get sample name and open the json file
  sample <- files[x]
  sample_name <- samples[x]
  open_file <- file(sample,open="r")
  jsonfile <- readLines(open_file)
  
  # For each line in the file 
  for (i in 1:length(jsonfile)) {
    line = jsonfile[i]
    
    # Search for the n_processed line
    # Get number of reads and store in dataframe
    if (grepl("num_compatible_fragments", line) == TRUE) {
      line = str_split(line,":")
      n_proc = line[[1]][1]
      n_proc = substr(n_proc,3,nchar(n_proc)-1)
      reads = line[[1]][2]
      reads <- str_remove_all(reads, ",")
      reads <- str_remove_all(reads, "\\s")
      
      reads.df[sample_name,"n_compatible_fragments"] <- reads
    }
  }
  # Close open file
  close(open_file)
}

# Save to table
write.table(reads.df,
            file = "n_compatible_fragments_salmon_uroscanseq.txt",
            sep="\t",quote = FALSE)
