###############################################################################################################
##################### Script title: Martineau_AlexARMS_bioinformatics_cleantable_seq1 #########################
#####################                      Project: AlexARMS                          #########################
#####################         Author: Gabrielle Martineau (mg62@hawaii.edu)           #########################
#####################                    Updated on: 2024-12-09                       #########################
###############################################################################################################



#####################################
########  A. Load libraries  ########
#####################################

# Install necessary packages if needed
necessarypackages = c("tidyverse", "data.table","seqinr", "here", "knitr", "remotes", "Biostrings", "GUniFrac")
# Install packages for that function
for (package in necessarypackages) {
  if (!require(package, character.only = TRUE)) {
    install.packages(package)
  }
}
# Load necessary libraries
lapply(necessarypackages, library, character.only = TRUE)

# Install specific taxize package and dependency bold
remotes::install_github("ropensci/bold")
library('bold')

remotes::install_github("ropensci/taxize")
library('taxize')



#####################################
########  B. Load functions  ########
#####################################

## B.1 Read in fasta files into data frame function
read_fasta <- function (fasta_file, name, path_to) {
  
  # Define path
  path_to = path_to
  
  # Construct the full file path
  full_path <- here(do.call(file.path, as.list(path_to)),  fasta_file)
  
  # Read in the DNA sequences from the FASTA file
  sequences_df <- readDNAStringSet(full_path)
  
  # Make a dataframe for sequence names and sequences
  seq_names <- names(sequences_df)
  sequences <- as.character(sequences_df)  # Convert DNAStringSet to character
  
  # Combine into a dataframe
  sequences_df_final <- data.frame(ID = seq_names, sequence = sequences, stringsAsFactors = FALSE)
  
  # Assign the result to a dynamically named object in the Global environment
  rownames(sequences_df_final) <- NULL
  
  assign(paste0("ASV_sequence_df_", name), sequences_df_final, envir = .GlobalEnv)
}


## B.2 Delete singletons and remove OTU/ASV that have an abundance of zero in all samples
singleton_removal <- function(abundance_table, treshold = 1, filename) {
  
  abundance_no_singleton = abundance_table
  
  # Extract numerical columns
  numerical_columns <- abundance_no_singleton %>%
    dplyr::select(where(is.numeric)) %>%
    colnames()
  
  # Make all singletons equal to zero 
  for (i in numerical_columns) {
    
    abundance_no_singleton[[i]] =   ifelse(abundance_no_singleton[[i]]  == treshold,  
                                           0,
                                           abundance_no_singleton[[i]])
  }
  
  # Delete ASV that have zero sum
  abundance_no_singleton <- abundance_no_singleton %>%
    dplyr:: mutate(sumASV = rowSums (across(where(is.numeric)))) %>%
    dplyr::filter(sumASV > 0) %>%
    dplyr::select(- (sumASV)) %>%
    dplyr::select(order(names(.))) %>%
    select(where(is.character), where(is.numeric))
  
  # Save abundance no singleton to global environment
  assign(paste0(filename, "_sequence_df_merged_nosg"), abundance_no_singleton, envir = .GlobalEnv)
}


## B.3 Merge abundance tables and save to csv
merge_abundance <- function (abundanceA, abundanceB, filename, sequenceA, sequenceB) {
  
  # Read in the abundance files
  path_to_abundance_and_fastq_A <- c("eDNAflow", "outputs", "zotus_052024")
  path_to_abundance_and_fastq_B <- c("eDNAflow", "outputs", "zotus_070224")
  
  
  abundanceA <- as.data.frame(
    fread(
      here(do.call(file.path, as.list(path_to_abundance_and_fastq_A)),  abundanceA),
      header = TRUE
    )
  )
  names(abundanceA)[1] <- "ID"
  
  abundanceB <- as.data.frame(
    fread(
      here(do.call(file.path, as.list(path_to_abundance_and_fastq_B)),  abundanceB),
      header = TRUE
    )
  )
  names(abundanceB)[1] <- "ID"
  
  
  # Rename the second abundance table
  # Extract the Zotu value of both tables
  numberA <- max(as.numeric(str_extract(abundanceA$ID, "\\d+")))
  numberB <- max(as.numeric(str_extract(abundanceB$ID, "\\d+")))
  
  extract_smallest <- function(numberA, numberB) {
    result <- list()
    
    if (numberA > numberB) {
      result$smallesttable <- "abundanceB"
      result$maxASV <- numberA
      result$lenght <- nrow(abundanceB)
    } else if (numberB > numberA) {
      result$smallesttable <- "abundanceA"
      result$maxASV <- numberB
      result$lenght <- nrow(abundanceA)
    } else {
      return("numberA and numberB are equal")
    }
    return(result)
  }
  
  # Extract results
  largest <- extract_smallest(numberA, numberB)
  rownames_new <-  paste("Zotu",  seq(from = largest$maxASV + 1, to = largest$lenght + largest$maxASV), sep ="")
  
  if (largest$smallesttable == "abundanceA") {
    rownames(abundanceA) <-  rownames_new 
  } else {
    rownames(abundanceB) <- rownames_new  
  }
  
  # Create index data frame to map back the old ASV name
  if (largest$smallesttable == "abundanceA") {
    index_df <- data.frame(original = abundanceA$ID, new = row.names(abundanceA))
    abundanceA$ID <- NULL
    abundanceA <- abundanceA %>%
      rownames_to_column(var = "ID")
  } else {
    index_df <- data.frame(original = abundanceB$ID, new = row.names(abundanceB))
    abundanceB$ID <- NULL
    abundanceB <- abundanceB %>%
      rownames_to_column(var = "ID")
  }
  
  
  # merge data frames
  merged_df <- full_join(abundanceA, abundanceB, by = intersect(names(abundanceA), names(abundanceB)))%>%
    mutate(across(everything(), ~ replace_na(., 0)))
  
  # add sequence info
  fastaB <- read_fasta(fasta_file = paste0("zotus_070224/", sequenceB), name = "B", path_to = c("eDNAflow", "outputs"))
  fastaA <- read_fasta(fasta_file = paste0("zotus_052024/", sequenceA), name = "A", path_to = c("eDNAflow", "outputs"))
  
  
  # Change the ID name for fastaB
  if (largest$smallesttable == "abundanceA") {
    correct_sequences <-  merge(fastaA, index_df, by.x = "ID", by.y = "original") %>%
      dplyr::select(c(sequence, new)) %>%
      dplyr::rename(ID = new)
  } else {
    correct_sequences <-  merge(fastaB, index_df, by.x = "ID", by.y = "original") %>%
      dplyr::select(c(sequence, new)) %>%
      dplyr::rename(ID = new)
  }
  
  
  # Row bind the two fasta files
  if (largest$smallesttable == "abundanceA") {
    combined_sequences <- bind_rows(fastaB, correct_sequences)
  } else {
    combined_sequences <- bind_rows(fastaA, correct_sequences)
  }
  
  
  # Merge the abundance value with sequences ID
  merged_df_combined <- merge(merged_df, combined_sequences, by = "ID")
  
  
  # Aggregate results
  aggregate.data.frame <- aggregate(x =  merged_df_combined[, !names(merged_df_combined) %in% c("ID","sequence")],
                                    by = list(sequence = merged_df_combined$sequence),
                                    sum) %>%
    dplyr::select(order(names(.))) %>%
    dplyr::arrange(sequence)
  
  # VERIFY: make a list of ID of same sequence
  # find unique sequences and check if they are the same
  sequences_left <- merged_df_combined %>%
    distinct(sequence, .keep_all = TRUE) %>%
    select(-c(ID))
  
  # Chedk if they are the same, regardless of the order:
  test_1 <- setequal(aggregate.data.frame$sequence, merged_df_combined$sequence)
  
  
  # Find duplicated sequence and sum them out and add to unique sequence
  # and check if it is equal to the aggragted data frame
  verification_df <- merged_df_combined %>%
    filter(sequence %in% sequence[duplicated(sequence)]) %>% # all duplicate instances
    group_by (sequence) %>%
    dplyr::summarise(across(where(is.numeric), sum), .groups = "keep")
  
  seuquences_non_duplicated <- merged_df_combined %>%
    filter(!(sequence %in% sequence[duplicated(sequence)]))%>%
    select(-c(ID)) %>%
    bind_rows(verification_df) %>% # add them up
    dplyr::select(order(names(.))) %>%
    dplyr::arrange(sequence)
  
  # Check if all duplicates were identical
  test_2 <- identical (aggregate.data.frame[, sapply(aggregate.data.frame, is.numeric)],
                       seuquences_non_duplicated[, sapply(seuquences_non_duplicated, is.numeric)])
  
  if (test_1 == TRUE & test_2 == TRUE) {
    message("Merging abundance tables was succesful")
  } else {
    message("Merging abundance tables failed")
  }
  
  # Add ID back to the aggregated data
  aggregate_with_ID <- merged_df_combined %>%
    select(c(ID, sequence)) %>%
    distinct(sequence, .keep_all = TRUE) %>%
    inner_join(aggregate.data.frame) %>%
    dplyr::select(order(names(.))) %>%
    select(where(is.character), where(is.numeric))
  
  
  # Make an output directory
  subdirectory_name <- "outputs"
  relative_path <- c("R_clean_tables", "ZoTU_Taxonomy")
  
  subdirectory_path <- here(do.call(file.path, as.list(relative_path)), subdirectory_name)
  
  if (!file.exists(subdirectory_path)) {
    # If it doesn't exist, create the subdirectory
    dir.create(subdirectory_path)
    message("Creating subdirectory: ", subdirectory_path)
  } else {
  }
  
  # save final data frame
  output_aggregated_table <- here(do.call(file.path, as.list(relative_path)), subdirectory_name, paste0(filename, "_merged_abundance", ".csv"))
  write.csv(aggregate_with_ID, output_aggregated_table, row.names = FALSE)
  
  # Return index files and abundance tables
  assign(paste0(filename, "_AbundanceA"),abundanceA, envir = .GlobalEnv)
  assign(paste0(filename, "_AbundanceB"), abundanceB, envir = .GlobalEnv)
  assign(paste0(filename, "_Indexes"), index_df, envir = .GlobalEnv)
  assign(paste0(filename, "_sequences_df_full"), merged_df_combined, envir = .GlobalEnv)
  assign(paste0(filename, "_sequence_df_merged"), aggregate_with_ID, envir = .GlobalEnv)
  
}


## B.4 Load blastn results, format in table and save to csv
load_blast <- function (blast_tab_file, db_name){
  
  # Construct the full path to the blastn file
  full_file_path <- here("Blastn", "output", blast_tab_file)
  
  
  # Lod blastn file in the tab format
  blastn_file <- fread(full_file_path)
  
  # Rename the columns
  names(blastn_file) <- c("query.id", 
                          "identity",   
                          "subject.id",
                          "scientific.name",
                          "common.name",
                          "super.kingdom",
                          "percent.identity", 
                          "alignment.length",
                          "query.length",
                          "sequence.length",
                          "mismatches", 
                          "gap.opens",
                          "gaps",
                          "qstart", 
                          "sequence.end",
                          "start.of.alignment",
                          "qend",
                          "subject.title",
                          "evalue", 
                          "bitscore",
                          "query.coverage",
                          "query.coverage.HSP")
  # Save the raw output to csv in the same folder directory recursive directory: create it
  subdirectory_name <- "outputs"
  
  relative_path <- c("R_clean_tables", "ZoTU_Taxonomy")
  
  subdirectory_path <- here(do.call(file.path, as.list(relative_path)), subdirectory_name)
  
  if (!file.exists(subdirectory_path)) {
    # If it doesn't exist, create the subdirectory
    dir.create(subdirectory_path)
  } else {
    message("The subdirectory already exists: ", subdirectory_path)
  }
  output_csv_path <- here(do.call(file.path, as.list(relative_path)), subdirectory_name, paste0(db_name, "_raw_blastn_output.csv"))
  write.csv(blastn_file, output_csv_path, row.names = FALSE)
  
  # Save that table to global.env (rename)
  blastn_var_name <- paste0(db_name, "_blastn_file_raw")
  assign(blastn_var_name, blastn_file, envir = .GlobalEnv)
  
  message("File saved as: ", output_csv_path)
}


## B.5 Filter blastn results based on different parameters and return a fasta file of resuts that didn't pass the test
filter_blast <- function (raw_blastn, min_perc_id, min_query_cov, abundance_table, db_name, return_fasta_and_remaining) {
  
  # Condition based on return_fasta_remaining argument)     
  if (return_fasta_and_remaining == TRUE){
    
    #Filter the blastn results based on query ID and % identity
    top_blastn <- raw_blastn %>%
      dplyr::filter(percent.identity  >= min_perc_id & query.coverage >= min_query_cov) %>%
      dplyr::group_by (query.id) %>%
      # Among them, keep the ones with top % identity
      filter(percent.identity == max (percent.identity, na.rm = TRUE)) %>%
      # Delete duplicates by scientific name
      dplyr::group_by (query.id) %>%
      filter(!(duplicated(scientific.name)))
    
    # Save object to environment
    assign(paste0(db_name, "_filtered_blast"), top_blastn ,envir = .GlobalEnv)
    
    # Import abundance 
    abundance_table <- abundance_table
    
    # Find the ASVs that did not make the cut out of the total ASV from beginning
    all_ASV <- abundance_table[, colnames(abundance_table) %in% "ID"]
    ASV_not_retained <- all_ASV[which(!(all_ASV %in% (unique(top_blastn$query.id))))]
    
    # Import sequences/ Extract sequences
    # sequence_table <- as.data.frame(data.table::fread (sequences_abundance_path), header = TRUE)
    sequences <- abundance_table[, colnames(abundance_table) %in% "sequence"]
    
    # Make a df that has ASV names and sequences
    # sequences_and_ASVs <- as.data.frame(bind_cols(all_ASV, sequences)) 
    sequences_and_ASVs <- abundance_table %>%
      select(c("ID", "sequence")) %>%
      dplyr::rename(all_ASV = ID, sequences = sequence) %>%
      dplyr::filter(!(all_ASV %in% (unique(top_blastn$query.id))))
    
    # Export that leftover ASV to environment
    assign(paste0(db_name, "_remaining_asvs"), sequences_and_ASVs ,envir = .GlobalEnv)
    
    # Make an output directory
    subdirectory_name <- "outputs"
    relative_path <- c("R_clean_tables", "ZoTU_Taxonomy")
    
    subdirectory_path <- here(do.call(file.path, as.list(relative_path)), subdirectory_name)
    
    if (!file.exists(subdirectory_path)) {
      # If it doesn't exist, create the subdirectory
      dir.create(subdirectory_path)
      message("Creating subdirectory: ", subdirectory_path)
    } else {
    }
    
    # Return a fasta file with only unique 
    uniqueseqs <- as.list(sequences_and_ASVs$sequences)
    seqnum <- sequences_and_ASVs$all_ASV
    
    # save fasta and export fasta file
    output_fasta <- here(do.call(file.path, as.list(relative_path)), subdirectory_name, paste0(db_name, "_sequences_not_retained.fasta"))
    write.fasta(uniqueseqs, seqnum, output_fasta)
    
    
    # Verify if the names of sequences in fasta files are the opposite of top blast results
    verify_test <- c(seqnum, unique(top_blastn$query.id))
    if (identical(sort(verify_test), sort(all_ASV))) {
      message("initial filtering successful")
    } else {
      message("initial filtering corrupted")
    }
    
    
    # Condition based on return_fasta_remaining arguemnt)          
  } else {
    
    #Filter the blastn results based on query ID and % identity
    top_blastn <- raw_blastn %>%
      dplyr::filter(percent.identity  >= min_perc_id & query.coverage >= min_query_cov) %>%
      dplyr::group_by (query.id) %>%
      # Among them, keep the ones with top % identity
      filter(percent.identity == max (percent.identity, na.rm = TRUE)) %>%
      # Delete duplicates by scientific name
      dplyr::group_by (query.id) %>%
      filter(!(duplicated(scientific.name)))
    
    # Save object to environment
    assign(paste0(db_name, "_filtered_blast"), top_blastn ,envir = .GlobalEnv)
    
  }
}


## B.6 Delete duplicated taxonomic ID in BLASTn search
delete_duplicated_ID <- function (filtered_blastn, db_name){
  
  # Delete duplicates
  no.duplicated.blast.searches <- filtered_blastn %>%
    group_by(query.id) %>%
    filter(!(duplicated(query.id)))
  
  # Verify if it worked
  if (length(unique(filtered_blastn$query.id)) == length(no.duplicated.blast.searches$query.id)) {
    print("remove duplicates successful")
  } else {
    print("remove duplicates failure")
  }
  # Save deduplicated df to environment
  assign(paste0(db_name, "_blastn_no_duplicates"), no.duplicated.blast.searches, envir = .GlobalEnv)
}


## B.7 Classify BLASTn results
classify_blast <- function (blastn_result, taxonomy_vector, db_name) {
  
  # Change column name 'query.id' to 'ASV'
  blastn_result <-blastn_result %>%
    dplyr::rename(ASV = query.id)
  
  # Separate subject.id and remove junk columns (sometimes the taxonomic ID is duplicated in the column)
  NAs.IDnumber <-blastn_result %>%
    separate(subject.id, into = c("subject.id", "subject.id.junk"), sep = ';') %>%
    dplyr::rename("subject.id" = "subject.id") %>%
    dplyr::select(-subject.id.junk)
  
  # Extract ID numbers
  all_ID_numbers <- as.numeric(NAs.IDnumber$subject.id)
  
  # Split in intervals of 500 because sometimes the classify functions lags and store them into a list
  split_into_intervals <- function(char_vector, interval_size) {
    num_intervals <- ceiling(length(char_vector) / interval_size)
    interval_list <- split(char_vector, ceiling(seq_along(char_vector) / interval_size))
    return(interval_list)
  }
  
  all_ID_numbers_list <- split_into_intervals(all_ID_numbers,500)
  
  # Write classification function
  classify_vector <- function(id_vector) {
    classification(sci_id = id_vector, callopts = list(), db = "ncbi", return_id = TRUE)
  }
  
  classification_test <- lapply(all_ID_numbers_list, classify_vector)
  classified.asv <-  do.call("c", classification_test)
  
  # Quick check to verify if the right numbers of IDs have been researched
  taxonomic_assignment = classified.asv 
  length(taxonomic_assignment) == length(all_ID_numbers)
  
  #Store ASV names into vector (add name for ASV that are non-unique in case the data frame still has duplicates)
  NAs.IDnumber$ASV.unique <- make.unique (NAs.IDnumber$ASV)
  ASVs_taxonomy <- NAs.IDnumber$ASV.unique
  
  # If there are no duplicates, ASV and ASV_unique should be the same
  NAs.IDnumber$ASV.unique == NAs.IDnumber$ASV
  
  # Extract each database for each ASV and store them in a list
  df_list <- lapply(taxonomic_assignment, \(x) as.data.frame(t(x))) |>
    setNames(paste0('database', ASVs_taxonomy))
  
  # Initialize a list to store the processed dataframes
  df_list_2 <- vector("list", length(all_ID_numbers))
  
  # Perform data wrangling for each table to clean up the data frame
  for (i in seq_along(df_list)) {
    # Extract the dataframe and transpose
    df <- df_list[[i]]
    # Set colnames as rank
    colnames(df) = c(df[which(rownames(df) == "rank"), ])
    # Remove the rank row
    df <- df[-which(rownames(df) == "rank"), ]
    # Remove unnecessary columns
    drop <- c("no rank", "clade")
    df <- df[, !(names(df) %in% drop)]
    # Bind the ASV name
    df <- cbind(df, ASV = ASVs_taxonomy[i])
    # Rename the last column to "taxonomic.ID"
    names(df)[length(names(df))] <- "taxonomic.ID"
    # Remove the column that says ID
    df <- df[-which(rownames(df) == "id"), ]
    # Store the dataframe in the list with a unique name
    df_list_2[[i]] <- df
  }
  
  # Merge taxonomic information for each ASV in the taxonomy matrix
  amerged.taxo.df.2 <- rbindlist(df_list_2, fill = TRUE) 
  
  # Ensure that amerged.taxo.df.2 has one row for each element of df_list_2
  nrow(amerged.taxo.df.2) == length(df_list_2) 
  
  
  # Get the common columns between what I wanna keep and what is found in the data frame
  common_cols <- intersect(colnames(amerged.taxo.df.2), taxonomy_vector)
  
  # Keep only columns in vector taxa
  amerged.taxo.df.3 <- amerged.taxo.df.2 %>%
    dplyr::select(c("taxonomic.ID", common_cols))
  
  # Change the NA as unclassified
  amerged.taxo.df.3[is.na(amerged.taxo.df.3)] <- "unclassified"
  amerged.taxo.df.3 = as.data.frame(amerged.taxo.df.3)
  
  # Change the ASV names in case there were some duplicates
  amerged.taxo.df.4 <- amerged.taxo.df.3 %>%
    dplyr::rename(ASV.unique = taxonomic.ID)
  
  amerged.taxo.df.4 <- left_join(amerged.taxo.df.4, NAs.IDnumber, by = "ASV.unique")
  
  amerged.taxo.df.5 <- amerged.taxo.df.4 %>%
    dplyr::select(c(common_cols, ID = ASV.unique)) 
  
  #Final tax table
  rownames(amerged.taxo.df.5) <- amerged.taxo.df.5$ID
  
  # Verify if nothing in the order has been messed up
  input_classification <- as.data.frame(blastn_result) %>%
    dplyr::select(ID = ASV, species = scientific.name)%>%
    separate(species, into = c("species1", "species2"), sep = ' ') %>%
    dplyr::select(c(ID, species1, species2))
  
  output_classification <- as.data.frame(amerged.taxo.df.5)%>%
    dplyr::select(ID, species) %>%
    separate(species, into = c("species1", "species2"), sep = ' ') %>%
    dplyr::select(c(ID, species1, species2))

  # Check if columns in both data frames are identical
  all_successful <- TRUE
  for (col in colnames(input_classification)) {
    if (!identical(input_classification[[col]], output_classification[[col]])) {
      all_successful <- FALSE
      print(paste("Column", col, "comparison: corrupted"))
    }
  }
  
  if(all_successful) {
    print("Classification successful")
  }
  
  # Save cleaned ASV table to environment
  assign(paste0(db_name, "_taxonomy_table"),amerged.taxo.df.5, envir = .GlobalEnv)
  
  amerged.taxo.df.5 <- amerged.taxo.df.5 %>%
    select(ID, everything())
  
  # Make an output directory
  subdirectory_name <- "outputs"
  relative_path <- c("R_clean_tables", "ZoTU_Taxonomy")
  
  subdirectory_path <- here(do.call(file.path, as.list(relative_path)), subdirectory_name)
  
  if (!file.exists(subdirectory_path)) {
    # If it doesn't exist, create the subdirectory
    dir.create(subdirectory_path)
    message("Creating subdirectory: ", subdirectory_path)
  } else {
  }
  
  # save final data frame
  output_classfied_table <- here(do.call(file.path, as.list(relative_path)), subdirectory_name, paste0(db_name, "_taxonomy_table.csv"))
  write.csv(amerged.taxo.df.5, output_classfied_table, row.names = FALSE)
}



## B.8 Merge tables function
merge_taxo_abundance <- function (taxonomyA, taxonomyB, abundanceA,abundanceB,abundance_global, filename){
  
  path_to_abundance_and_taxonomy <- c("R_clean_tables", "ZoTU_Taxonomy", "outputs")
  
  
  if (missing(abundanceA)){
    path_to_abundance_and_taxonomy <- c("R_clean_tables", "ZoTU_Taxonomy", "outputs")
    
    taxonomyA <-  as.data.frame(fread(here(do.call(file.path, as.list(path_to_abundance_and_taxonomy)),  taxonomyA),
        header = TRUE))
    names(taxonomyA)[1] = "ID"
    
    taxonomyB <- as.data.frame(fread(here(do.call(file.path, as.list(path_to_abundance_and_taxonomy)),  taxonomyB),
                                     header = TRUE))
    names(taxonomyB)[1] = "ID"
    
    abundance_global <- as.data.frame(fread(here(do.call(file.path, as.list(path_to_abundance_and_taxonomy)),  abundance_global),
                                            header = TRUE))
    names(abundance_global)[1] <- "ID"
    
    # Get the column names of taxonomyA
    taxonomy_colnames <- colnames(taxonomyA)
    # Remove "ID" from the column names
    taxonomy_colnames <- setdiff(taxonomy_colnames, "ID")
    
    ## Join two taxonomy tables into total.taxonomy
    total.taxonomy <- rbind(taxonomyA, taxonomyB)
    row.names(total.taxonomy) = total.taxonomy$ID
    
    ## Verify if there are duplicates
    nrow(total.taxonomy) == length(unique(row.names(total.taxonomy)))
    
    ## Join combined tax table and abundance_global tables together
    a.combined_table <- inner_join(abundance_global, total.taxonomy, by ="ID")
    
    ## Save abundance table: that one contains only the IDs that were retained in the pipeline (in abundance)
    a.final.abundance.table <-  a.combined_table %>%
      dplyr::select(colnames(abundance_global))
    row.names(a.final.abundance.table)<- a.final.abundance.table$ID
    
    
    ## G.2 Delete duplicates in both the abundance and tax tables
    # Find the number of unique taxa
    rownames(a.combined_table) = a.combined_table$ID
    unique_taxa <- a.combined_table %>%
      distinct(across(all_of(taxonomy_colnames)))
    
    
    number_unique_taxa <- nrow(unique_taxa)
    # Find the unique taxa and their mapping ASV
    unique_taxa$mapping_ASV <- rownames(unique_taxa)
    mapping_taxa <- unique_taxa %>%
      dplyr::select(c(all_of(taxonomy_colnames)), mapping_ASV)
    
    # Remove duplicates from total.taxonomy
    total.taxonomy$ID <- rownames(total.taxonomy)
    
    duplicated_tax_table <- total.taxonomy %>%
      add_count(across(all_of(taxonomy_colnames))) %>%
      pivot_wider(names_from = ID, values_from = n) 
    
    # Add mapping taxa and replace taxa name by mapping taxa ID in duplicated abundance table
    duplicated_tax_table_2 = right_join(mapping_taxa, duplicated_tax_table)
    rownames(duplicated_tax_table_2) <- duplicated_tax_table_2$mapping_ASV
    duplicated_tax_table_2 = as.data.frame(t(duplicated_tax_table_2[,!(colnames(duplicated_tax_table_2) %in% c("mapping_ASV", taxonomy_colnames))]))
    
    duplicated_test <- which(!is.na(duplicated_tax_table_2), arr.ind=TRUE)
    duplicated_tax_table_2[duplicated_test] <- rownames(duplicated_test)
    
    # Try to clean it up and remove NAs
    #define a function
    RemoveNAs <- function(x, size) {
      temp <- x[!is.na(x)]
      c(temp, rep("", size - length(temp)))
    }
    
    # Calculate longest column size in df (non-NA)
    Max <- max(colSums(!is.na(duplicated_tax_table_2)))
    
    # Remove NAs from duplicated_tax_table_2
    Indexes <- setDT(duplicated_tax_table_2)[, lapply(.SD, RemoveNAs, Max)]
    
    # Some verification:
    # Check if the most abundant unique species has the right number of columns
    max.count = duplicated_tax_table <- total.taxonomy%>%
      add_count(across(all_of(taxonomy_colnames)))
    max.count = max(max.count$n)
    
    if (max.count == Max & max.count == nrow(Indexes)){
      print ("V1 good")
    } else {
      print ("V1 is wrong")
    }
    
    
    # Check if the number of unique ASV is good
    if (ncol(Indexes) == number_unique_taxa){
      print( "V2 good" )
    } else {
      print ("V2 is wrong")
    }
    
    # Delete duplicated rows on tax table and abundace tables
    duplicated_tax_table_1 <- total.taxonomy %>%
      add_count(across(all_of(taxonomy_colnames))) %>%
      pivot_wider(names_from = ID, values_from = n)
    
    tax_table_no_duplicate = right_join (mapping_taxa, duplicated_tax_table_1) %>%
      dplyr::select(c(all_of(taxonomy_colnames)), mapping_ASV) 
    rownames(tax_table_no_duplicate) <- tax_table_no_duplicate$mapping_ASV
    
    
    # Delete duplicates on abundance table
    # Add mapping ASV column to the original tax table
    index_taxonomy <- merge(total.taxonomy, mapping_taxa) 
    
    # Check if all the mapping taxa are present
    if (length(unique(mapping_taxa$mapping_ASV)) == length(unique(index_taxonomy$mapping_ASV)))  {
      print( "V3 good" )
    } else {
      print ("V3 is wrong")
    }
    
    abundance_no_duplicates <- inner_join(abundance_global, index_taxonomy)%>%
      dplyr::select(-c(all_of(taxonomy_colnames)))
    
    abundance_no_duplicates_1 <- abundance_no_duplicates %>%
      dplyr::select(-c(ID))
    
    abundance_no_duplicates_2 = aggregate (. ~ mapping_ASV, abundance_no_duplicates_1, sum)
    ## G.3 Verification steps: does the total equals to the sum of its part for different indexes?
    # Create data frames to store the results
    indexes_test_df <-  data.frame(matrix(nrow = nrow(Indexes), ncol = ncol(Indexes)))
    colnames(indexes_test_df) <- colnames(Indexes)
    indexes_test_df[] <- NA
    caculation_test_df <-  data.frame(matrix(nrow = ncol(abundance_no_duplicates_2) -1, ncol = ncol(Indexes)))
    colnames(caculation_test_df) <- colnames(Indexes)
    
    for (index in colnames(Indexes)){
      I <- Indexes[[index]]
      V <- abundance_no_duplicates %>%
        dplyr::filter(mapping_ASV %in% index)%>%
        as.data.frame()
      
      v.vector <-  V$ID %in% I
      replacement_length <- length(v.vector)
      indexes_test_df[, index] <- c(as.character(v.vector), rep(NA, nrow(indexes_test_df) - replacement_length))
      
      columns_to_exclude <- c("mapping_ASV", "ID")
      V.1 <- colSums(V[ , !(names(V) %in% columns_to_exclude)])
      caculation_test_df[,index] <- as.character(V.1 ==  abundance_no_duplicates_2[which(abundance_no_duplicates_2$mapping_ASV == index), !colnames(abundance_no_duplicates_2) %in% columns_to_exclude])
    }
    if (any(is.na(indexes_test_df)) || any(is.na(caculation_test_df))){
      print ("V4 is good")
    } else {
      print ("V4 wrong")
    }
    
    
    # Check if there are still duplicated species
    if (length(unique(abundance_no_duplicates$mapping_ASV)) == length(unique(index_taxonomy$species))) {
      print( "V5 good" )
      
      rownames(abundance_no_duplicates_2) = abundance_no_duplicates_2$mapping_ASV
      mapping_ASV <- "mapping_ASV"
      abundance_no_duplicates_2 <- abundance_no_duplicates_2[, c(mapping_ASV, setdiff(names(abundance_no_duplicates_2), mapping_ASV))]
      colnames(abundance_no_duplicates_2)[1] = "ID"
      
      rownames(tax_table_no_duplicate) = tax_table_no_duplicate$mapping_ASV
      tax_table_no_duplicate <- tax_table_no_duplicate[, c(mapping_ASV, setdiff(names(tax_table_no_duplicate), mapping_ASV))]
      colnames(tax_table_no_duplicate)[1] = "ID"
      
      # Save cleaned ASV table to environment
      assign(paste0(filename, "_taxonomy_table"),tax_table_no_duplicate , envir = .GlobalEnv)
      assign(paste0(filename, "_abundance_table"),abundance_no_duplicates_2 , envir = .GlobalEnv)
      
      # Make output directory if not already made
      subdirectory_name <- "outputs"
      relative_path <- c("R_clean_tables", "ZoTU_Taxonomy")
      
      subdirectory_path <- here(do.call(file.path, as.list(relative_path)), subdirectory_name)
      
      if (!file.exists(subdirectory_path)) {
        # If it doesn't exist, create the subdirectory
        dir.create(subdirectory_path)
        message("Creating subdirectory: ", subdirectory_path)
      } else {
      }
      
      # Save abundance 
      # save final data frame
      output_abundance <- here(do.call(file.path, as.list(relative_path)), subdirectory_name, paste0(filename, "_abundance_table.csv"))
      write.csv(abundance_no_duplicates_2, output_abundance, row.names = FALSE)
      
      output_taxonomy <- here(do.call(file.path, as.list(relative_path)), subdirectory_name, paste0(filename, "_taxonomy_table.csv"))
      write.csv(tax_table_no_duplicate, output_taxonomy, row.names = FALSE)
      
      
    } else {
      rownames(abundance_no_duplicates_2) = abundance_no_duplicates_2$mapping_ASV
      mapping_ASV <- "mapping_ASV"
      abundance_no_duplicates_2 <- abundance_no_duplicates_2[, c(mapping_ASV, setdiff(names(abundance_no_duplicates_2), mapping_ASV))]
      colnames(abundance_no_duplicates_2)[1] = "ID"
      
      rownames(tax_table_no_duplicate) = tax_table_no_duplicate$mapping_ASV
      tax_table_no_duplicate <- tax_table_no_duplicate[, c(mapping_ASV, setdiff(names(tax_table_no_duplicate), mapping_ASV))]
      colnames(tax_table_no_duplicate)[1] = "ID"
      
      duplicated_species(taxonomy = tax_table_no_duplicate, abundance = abundance_no_duplicates_2, filename = "test" )
      if   (length(unique(test_abundance_table$ID)) == length(unique(index_taxonomy$species))) {
        print( "V5 good" )
      } else {
        print ("V5 is wrong")
      }
      # Save cleaned ASV table to environment
      assign(paste0(filename, "_taxonomy_table"),test_taxonomy_table, envir = .GlobalEnv)
      assign(paste0(filename, "_abundance_table"),test_abundance_table, envir = .GlobalEnv)
      
      # Make output directory if not already made
      subdirectory_name <- "outputs"
      relative_path <- c("R_clean_tables", "ZoTU_Taxonomy")
      
      subdirectory_path <- here(do.call(file.path, as.list(relative_path)), subdirectory_name)
      
      if (!file.exists(subdirectory_path)) {
        # If it doesn't exist, create the subdirectory
        dir.create(subdirectory_path)
        message("Creating subdirectory: ", subdirectory_path)
      } else {
      }
      
      
      output_abundance <- here(do.call(file.path, as.list(relative_path)), subdirectory_name, paste0(filename, "_abundance_table.csv"))
      write.csv(test_abundance_table, output_abundance, row.names = FALSE)
      
      output_taxonomy <- here(do.call(file.path, as.list(relative_path)), subdirectory_name, paste0(filename, "_taxonomy_table.csv"))
      write.csv(test_taxonomy_table, output_taxonomy, row.names = FALSE)
    
      
      
    }
    
    
    # Condition based on abundance A argument missing or not   
  } else {
    
    path_to_abundance_and_taxonomy <- c("R_clean_tables", "ZoTU_Taxonomy", "outputs")
    taxonomyA <-  as.data.frame(fread(here(do.call(file.path, as.list(path_to_abundance_and_taxonomy)),  taxonomyA),
                                      header = TRUE))
    names(taxonomyA)[1] = "ID"
    
    taxonomyB <- as.data.frame(fread(here(do.call(file.path, as.list(path_to_abundance_and_taxonomy)),  taxonomyB),
                                     header = TRUE))
    names(taxonomyB)[1] = "ID"
    
    abundanceA <- as.data.frame(fread(here(do.call(file.path, as.list(path_to_abundance_and_taxonomy)),  abundanceA ),
                                           header = TRUE))
    names(abundanceA)[1] <- "ID"
    
    abundanceB <- as.data.frame(fread(here(do.call(file.path, as.list(path_to_abundance_and_taxonomy)),  abundanceB ),
                                      header = TRUE))
    names(abundanceB)[1] <- "ID"
    
    # Change the names for abundance and taxonomy tables B (just to make sure they are not duplicated)
    combined_table_B <- merge(abundanceB, taxonomyB, by = "ID")
    
    # Rename the ASV
    rownames(combined_table_B)=   ASV.sequences = paste("ASV",  seq(from = 10000,to = (9999  + nrow(combined_table_B))), sep ="")
    combined_table_B$mapping = row.names(combined_table_B)
    
    # Save the mapping in another data frame
    Mapping_ASV <- as.data.frame(cbind(combined_table_B$ID, combined_table_B$mapping))
    
    # Resplit them
    combined_table_B$ID <- NULL
    combined_table_B$mapping <- NULL
    
    # Find out how many taxa are in taxonomy table
    taxonomy_colnames_B <- colnames(taxonomyB)
    # Remove "ID" from the column names
    taxonomy_colnames_B <- setdiff(taxonomy_colnames_B, "ID")
    
    taxonomy.B.rename <- combined_table_B[, colnames(combined_table_B) %in% taxonomy_colnames_B]
    taxonomy.B.rename$ID = row.names(taxonomy.B.rename)
    
    abundance.B.rename <- combined_table_B[, !(colnames(combined_table_B) %in% taxonomy_colnames_B)]
    abundance.B.rename$ID <- row.names(abundance.B.rename)
    
    # Get the common taxonomy colnames
    taxonomy.A.colnames <- colnames(taxonomyA)
    
    # Take the intersection of taxonomy A and B colnames
    taxonomy_colnames <- intersect(taxonomy.A.colnames, taxonomy_colnames_B)
    
    # Remove "ID" from the column names
    taxonomy_colnames <- setdiff(taxonomy_colnames, "ID")
    
    ## Join two taxonomy tables into total.taxonomy
    total.taxonomy <- rbind(taxonomyA, taxonomy.B.rename)
    row.names(total.taxonomy) = total.taxonomy$ID
    
    ## Verify if there are duplicates
    nrow(total.taxonomy) == length(unique(row.names(total.taxonomy)))
    
    # Combine abundance table
    abundance_combined <- full_join(abundanceA, abundance.B.rename)
    abundance_combined[is.na(abundance_combined)] <- 0
    
    
    ## Join combined tax table and  abundance_combined tables together
    a.combined_table <- inner_join( abundance_combined, total.taxonomy, by ="ID")
    
    ## Save abundance table: that one contains only the IDs that were retained in the pipeline (in abundance)
    a.final.abundance.table <-  a.combined_table %>%
      dplyr::select(colnames( abundance_combined))
    row.names(a.final.abundance.table)<- a.final.abundance.table$ID
    
    
    ## G.2 Delete duplicates in both the abundance and tax tables
    # Find the number of unique taxa
    rownames(a.combined_table) = a.combined_table$ID
    unique_taxa <- a.combined_table %>%
      distinct(across(all_of(taxonomy_colnames)))
    
    
    number_unique_taxa <- nrow(unique_taxa)
    # Find the unique taxa and their mapping ASV
    unique_taxa$mapping_ASV <- rownames(unique_taxa)
    mapping_taxa <- unique_taxa %>%
      dplyr::select(c(all_of(taxonomy_colnames)), mapping_ASV)
    
    # Remove duplicates from total.taxonomy
    total.taxonomy$ID <- rownames(total.taxonomy)
    
    duplicated_tax_table <- total.taxonomy %>%
      add_count(across(all_of(taxonomy_colnames))) %>%
      pivot_wider(names_from = ID, values_from = n) 
    
    # Add mapping taxa and replace taxa name by mapping taxa ID in duplicated abundance table
    duplicated_tax_table_2 = right_join(mapping_taxa, duplicated_tax_table)
    rownames(duplicated_tax_table_2) <- duplicated_tax_table_2$mapping_ASV
    duplicated_tax_table_2 = as.data.frame(t(duplicated_tax_table_2[,!(colnames(duplicated_tax_table_2) %in% c("mapping_ASV", taxonomy_colnames))]))
    
    duplicated_test <- which(!is.na(duplicated_tax_table_2), arr.ind=TRUE)
    duplicated_tax_table_2[duplicated_test] <- rownames(duplicated_test)
    
    # Try to clean it up and remove NAs
    #define a function
    RemoveNAs <- function(x, size) {
      temp <- x[!is.na(x)]
      c(temp, rep("", size - length(temp)))
    }
    
    # Calculate longest column size in df (non-NA)
    Max <- max(colSums(!is.na(duplicated_tax_table_2)))
    
    # Remove NAs from duplicated_tax_table_2
    Indexes <- setDT(duplicated_tax_table_2)[, lapply(.SD, RemoveNAs, Max)]
    
    # Some verification:
    # Check if the most abundant unique species has the right number of columns
    max.count = duplicated_tax_table <- total.taxonomy%>%
      add_count(across(all_of(taxonomy_colnames)))
    max.count = max(max.count$n)
    
    if (max.count == Max & max.count == nrow(Indexes)){
      print ("V1 good")
    } else {
      print ("V1 is wrong")
    }
    
    
    # Check if the number of unique ASV is good
    if (ncol(Indexes) == number_unique_taxa){
      print( "V2 good" )
    } else {
      print ("V2 is wrong")
    }
    
    # Delete duplicated rows on tax table and abundace tables
    duplicated_tax_table_1 <- total.taxonomy %>%
      add_count(across(all_of(taxonomy_colnames))) %>%
      pivot_wider(names_from = ID, values_from = n)
    
    tax_table_no_duplicate = right_join (mapping_taxa, duplicated_tax_table_1) %>%
      dplyr::select(c(all_of(taxonomy_colnames)), mapping_ASV) 
    rownames(tax_table_no_duplicate) <- tax_table_no_duplicate$mapping_ASV
    
    # Delete duplicates on abundance table
    # Add mapping ASV column to the original tax table
    index_taxonomy <- merge(total.taxonomy, mapping_taxa) 
    
    # Check if all the mapping taxa are present
    if (length(unique(mapping_taxa$mapping_ASV)) == length(unique(index_taxonomy$mapping_ASV)))  {
      print( "V3 good" )
    } else {
      print ("V3 is wrong")
    }
    
    abundance_no_duplicates <- inner_join( abundance_combined, index_taxonomy)%>%
      dplyr::select(-c(all_of(taxonomy_colnames)))
    
    abundance_no_duplicates_1 <- abundance_no_duplicates %>%
      dplyr::select(-c(ID))
    
    abundance_no_duplicates_2 = aggregate (. ~ mapping_ASV, abundance_no_duplicates_1, sum)
    
    ## G.3 Verification steps: does the total equals to the sum of its part for different indexes?
    # Create data frames to store the results
    indexes_test_df <-  data.frame(matrix(nrow = nrow(Indexes), ncol = ncol(Indexes)))
    colnames(indexes_test_df) <- colnames(Indexes)
    indexes_test_df[] <- NA
    caculation_test_df <-  data.frame(matrix(nrow = ncol(abundance_no_duplicates_2) -1, ncol = ncol(Indexes)))
    colnames(caculation_test_df) <- colnames(Indexes)
    
    for (index in colnames(Indexes)){
      I <- Indexes[[index]]
      V <- abundance_no_duplicates %>%
        dplyr::filter(mapping_ASV %in% index)%>%
        as.data.frame()
      
      v.vector <-  V$ID %in% I
      replacement_length <- length(v.vector)
      indexes_test_df[, index] <- c(as.character(v.vector), rep(NA, nrow(indexes_test_df) - replacement_length))
      
      columns_to_exclude <- c("mapping_ASV", "ID")
      V.1 <- colSums(V[ , !(names(V) %in% columns_to_exclude)])
      caculation_test_df[,index] <- as.character(V.1 ==  abundance_no_duplicates_2[which(abundance_no_duplicates_2$mapping_ASV == index), !colnames(abundance_no_duplicates_2) %in% columns_to_exclude])
    }
    if (any(is.na(indexes_test_df)) || any(is.na(caculation_test_df))){
      print ("V4 is good")
    } else {
      print ("V4 wrong")
    }
    
    # Check if there are still duplicated species
    if (length(unique(abundance_no_duplicates$mapping_ASV)) == length(unique(index_taxonomy$species))) {
      print( "V5 good" )
      
      rownames(abundance_no_duplicates_2) = abundance_no_duplicates_2$mapping_ASV
      mapping_ASV <- "mapping_ASV"
      abundance_no_duplicates_2 <- abundance_no_duplicates_2[, c(mapping_ASV, setdiff(names(abundance_no_duplicates_2), mapping_ASV))]
      colnames(abundance_no_duplicates_2)[1] = "ID"
      
      rownames(tax_table_no_duplicate) = tax_table_no_duplicate$mapping_ASV
      tax_table_no_duplicate <- tax_table_no_duplicate[, c(mapping_ASV, setdiff(names(tax_table_no_duplicate), mapping_ASV))]
      colnames(tax_table_no_duplicate)[1] = "ID"
      
      # Save cleaned ASV table to environment
      assign(paste0(filename, "_taxonomy_table"),tax_table_no_duplicate, envir = .GlobalEnv)
      assign(paste0(filename, "_abundance_table"),abundance_no_duplicates_2, envir = .GlobalEnv)
      
      # Make output directory if not already made
      subdirectory_name <- "outputs"
      relative_path <- c("R_clean_tables", "ZoTU_Taxonomy")
      
      subdirectory_path <- here(do.call(file.path, as.list(relative_path)), subdirectory_name)
      
      if (!file.exists(subdirectory_path)) {
        # If it doesn't exist, create the subdirectory
        dir.create(subdirectory_path)
        message("Creating subdirectory: ", subdirectory_path)
      } else {
      }
      
      
      output_abundance <- here(do.call(file.path, as.list(relative_path)), subdirectory_name, paste0(filename, "_abundance_table.csv"))
      write.csv(abundance_no_duplicates_2, output_abundance, row.names = FALSE)
      
      output_taxonomy <- here(do.call(file.path, as.list(relative_path)), subdirectory_name, paste0(filename, "_taxonomy_table.csv"))
      write.csv(tax_table_no_duplicate, output_taxonomy, row.names = FALSE)
      
      
      
    } else {
      rownames(abundance_no_duplicates_2) = abundance_no_duplicates_2$mapping_ASV
      mapping_ASV <- "mapping_ASV"
      abundance_no_duplicates_2 <- abundance_no_duplicates_2[, c(mapping_ASV, setdiff(names(abundance_no_duplicates_2), mapping_ASV))]
      colnames(abundance_no_duplicates_2)[1] = "ID"
      
      rownames(tax_table_no_duplicate) = tax_table_no_duplicate$mapping_ASV
      tax_table_no_duplicate <- tax_table_no_duplicate[, c(mapping_ASV, setdiff(names(tax_table_no_duplicate), mapping_ASV))]
      colnames(tax_table_no_duplicate)[1] = "ID"
      
      duplicated_species(taxonomy = tax_table_no_duplicate, abundance = abundance_no_duplicates_2, filename = "test" )
      if   (length(unique(test_abundance_table$ID)) == length(unique(index_taxonomy$species))) {
        print( "V5 good" )
      } else {
        print ("V5 is wrong")
      }
      # Save cleaned ASV table to environment
      assign(paste0(filename, "_taxonomy_table"),test_taxonomy_table, envir = .GlobalEnv)
      assign(paste0(filename, "_abundance_table"),test_abundance_table, envir = .GlobalEnv)
      
      # Make output directory if not already made
      subdirectory_name <- "outputs"
      relative_path <- c("R_clean_tables", "ZoTU_Taxonomy")
      
      subdirectory_path <- here(do.call(file.path, as.list(relative_path)), subdirectory_name)
      
      if (!file.exists(subdirectory_path)) {
        # If it doesn't exist, create the subdirectory
        dir.create(subdirectory_path)
        message("Creating subdirectory: ", subdirectory_path)
      } else {
      }
      
      
      output_abundance <- here(do.call(file.path, as.list(relative_path)), subdirectory_name, paste0(filename, "_abundance_table.csv"))
      write.csv(test_abundance_table, output_abundance, row.names = FALSE)
      
      output_taxonomy <- here(do.call(file.path, as.list(relative_path)), subdirectory_name, paste0(filename, "_taxonomy_table.csv"))
      write.csv(test_taxonomy_table, output_taxonomy, row.names = FALSE)
      
      
      
    }
    
  }
}




#####################################
########  C. Data Analyses  #########
#####################################

## C.1 Adjust abundance tables so the column names match for each sample in both sequencing runs 
ASV_table_A <- as.data.frame(fread (here(do.call(file.path, as.list(c("eDNAflow", "outputs", "zotus_052024"))),  "zotu_table.tsv") , header = TRUE))
ASV_table_A_modified <-ASV_table_A %>%
  rename_with(~ gsub("^(([^_]*_){2}).*", "\\1", .)) %>%
  rename_with(~ gsub("_$", "", .)) %>%
  write_csv(here("eDNAflow", "outputs", "zotus_052024", "zotu_table_names_modified.csv"))


ASV_table_B <- as.data.frame(fread (here(do.call(file.path, as.list(c("eDNAflow", "outputs", "zotus_070224"))),  "zotu_table.tsv") , header = TRUE))
ASV_table_B_modified <-ASV_table_B %>%
  rename_with(~ gsub("^(([^_]*_){2}).*", "\\1", .)) %>%
  rename_with(~ gsub("_$", "", .)) %>%
  write_csv(here("eDNAflow", "outputs", "zotus_070224", "zotu_table_names_modified.csv"))


## C.2 Merge abundance tables
merge_abundance(abundanceA = "zotu_table_names_modified.csv", 
                abundanceB = "zotu_table_names_modified.csv", 
                sequenceA = "eNAflowAlex_zotus.fasta",
                sequenceB = "eDNAflowAlex_070624_zotus.fasta",
                filename = "both_seqruns")


## C.3 Make histogram of sequencing depth per sample and rarefy
# Rarefy data frame
only_numeric <- both_seqruns_sequence_df_merged %>%
  select(where(is.numeric))

# Get the sum of reads per sample
sum_per_sample <- colSums(only_numeric)
sum_per_sample_df <- data.frame(
  sample = names(sum_per_sample),
  sum_sample = as.numeric(sum_per_sample)) %>%
  mutate(sample = as.character(sample))

# Include island information with meta data file
meta_data <- read_csv(here("R_clean_tables", "MetaData", "Alex_metadata_040824.csv"), show_col_types = FALSE) 
meta_data <- meta_data %>%
  mutate(sample_name_new = sub("^(([^_]*_){2}).*", "\\1", meta_data$Sample)) %>%
  mutate(sample_name_new =  as.character(sub("_$", "", sample_name_new)) )

# Merge sample information with sum of reads
hist_meta <- merge (meta_data,sum_per_sample_df, by.x = "sample_name_new",  by.y = "sample") %>% 
  dplyr::select(c(sample_name_new, sum_sample, IslandName)) %>%
  dplyr::group_by(IslandName) %>%
  dplyr::arrange(sum_sample)

# Identify the maximum x-position of samples with sum_sample smaller than 10000 and 20000
max_x_position <- max(which(hist_meta$sum_sample < 10000)) # red line
max_x_position_2 <- max(which(hist_meta$sum_sample < 20000)) # blue line

# plot histogram
histogram <- ggplot(hist_meta, aes(x = sample_name_new, y = sum_sample)) +
  geom_bar (aes(fill = IslandName),stat= "identity",  width = 0.5) +
  scale_x_discrete (limits = hist_meta$sample_name_new)+
  theme_classic() +
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust = 0.5)) +
  ylim(0, max(hist_meta$sum_sample)) +
  geom_vline(xintercept = max_x_position + 0.5, col = "red", size = 1) +
  geom_vline(xintercept = max_x_position_2 + 0.5, col = "blue", size = 1) +
  ylab("Read abundance") +
  xlab("Sample") +
  labs(fill = "Island") 

histogram

# Export histogram in outputs
ggsave(filename = here("R_clean_tables", "ZoTU_Taxonomy", "outputs", "bothseqruns_histogram_raw_read_abundance.png"),
       plot = histogram, 
       width = 12,
       height = 7)

# Rarefy only for samples below above 20K
samples_below <- sum_per_sample_df[which(sum_per_sample_df$sum_sample < 20000), "sample"]
lenght <- length(samples_below)  # 8 samples

# Flip table for rarefaction
rarefaction_df <- as.data.frame(t(both_seqruns_sequence_df_merged)) %>%
  rename_with(~ both_seqruns_sequence_df_merged[,1 ]) 

# Index necessary columns
rarefaction_df_to_rarefy <- rarefaction_df %>%
  filter(!rownames(.) %in% c("ID", "sequence", samples_below)) %>%
  mutate(across(everything(), as.numeric)) %>%  
  Rarefy(20000) %>%  # Rarefy to 20000
  .$otu.tab.rff  %>% 
  t() %>%
  as.data.frame() %>%
  rownames_to_column("ID")

# Add samples non-rarefied back
rarefaction_combined <- both_seqruns_sequence_df_merged %>%
  dplyr::select(c(samples_below, "ID", "sequence")) %>% 
  inner_join(rarefaction_df_to_rarefy) %>%
  dplyr::select(order(names(.))) %>%
  select(where(is.character), where(is.numeric))

# Validation: compare sample sums in new table vs. after
colSum_rarefied <- colSums(rarefaction_combined[, !(colnames(rarefaction_combined) %in% c(samples_below, "ID", "sequence"))])

colSum_non_rarefied <- colSums(rarefaction_combined[, colnames(rarefaction_combined) %in% c(samples_below)])
colSum_non_rarefied_original <- colSums(both_seqruns_sequence_df_merged[, colnames(both_seqruns_sequence_df_merged) %in% c(samples_below)])

rarefaction_treshold = 20000

if (length(colSum_rarefied)== ncol(rarefaction_df_to_rarefy[, sapply(rarefaction_df_to_rarefy, is.numeric)])
 & all(colSum_non_rarefied == colSum_non_rarefied_original)
 & all(colSum_rarefied == rarefaction_treshold)) {
  message("Rarefaction successful")
} else {
  message("Rarefaction failed")
}


## C.4 Remove singleton
singleton_removal(abundance_table = rarefaction_combined, treshold = 1, filename = "both_seqruns")


## C.5 Write fasta file and save rarefied (no sing) abundance table and corresponding sequences
# Create fasta
uniqueseqs <- as.list(both_seqruns_sequence_df_merged_nosg$sequence)
seqnum <- both_seqruns_sequence_df_merged_nosg$ID

# Make an output directory
subdirectory_name <- "outputs"
relative_path <- c("R_clean_tables", "ZoTU_Taxonomy")

subdirectory_path <- here(do.call(file.path, as.list(relative_path)), subdirectory_name)

if (!file.exists(subdirectory_path)) {
  # If it doesn't exist, create the subdirectory
  dir.create(subdirectory_path)
  message("Creating subdirectory: ", subdirectory_path)
} else {
  message("subdirectory already created")
}

# save fasta and export fasta file
output_fasta <- here(do.call(file.path, as.list(relative_path)), subdirectory_name, paste0("both_seqruns_","rarefied_nosg", ".fasta"))
# write abundance and fasta files
write.fasta(uniqueseqs, seqnum, output_fasta)

# Save rarefied table
output_rarefied_table_nosg <- here(do.call(file.path, as.list(relative_path)), subdirectory_name, paste0("both_seqruns_","merged_abundance_rarefied_nosg",  ".csv"))
write.csv(both_seqruns_sequence_df_merged_nosg, output_rarefied_table_nosg, row.names = FALSE)


# Save rarefied table without sequences
both_seqruns_sequence_df_merged_nosg_no_sequence <- both_seqruns_sequence_df_merged_nosg 
both_seqruns_sequence_df_merged_nosg_no_sequence$sequence <- NULL

output_rarefied_table_nosg_noseq <- here(do.call(file.path, as.list(relative_path)), subdirectory_name, paste0("both_seqruns_","merged_abundance_rarefied_nosg_noseq",  ".csv"))
write.csv(both_seqruns_sequence_df_merged_nosg_no_sequence, output_rarefied_table_nosg_noseq, row.names = FALSE)



###############################################################
########  D. BLASTn taxonomy outputs: sponge data base ########
###############################################################

## D.1 Load blastn results, format in table and save to csv     
load_blast("blast_results_spongedb_120924.tab", db_name = "sponge_both_seqruns")


## D.2 Filter sponge db blast results
filter_blast (raw_blastn = sponge_both_seqruns_blastn_file_raw,
              min_perc_id = 95,
              min_query_cov = 95,
              abundance_table = both_seqruns_sequence_df_merged_nosg,
              db_name = "sponge_both_seqruns",
              return_fasta_and_remaining = TRUE)


## D.3 Delete duplicates in filtered BLASTn results from sponge data base
delete_duplicated_ID(filtered_blastn = sponge_both_seqruns_filtered_blast, db_name = "sponge_both_seqruns")


## D.4 Classify filtered BLASTn results from sponge data base
classify_blast(blastn_result = sponge_both_seqruns_blastn_no_duplicates, 
               taxonomy_vector = c("kingdom", "phylum", "class", "order", "family", "genus", "species"), 
               db_name = "sponge_both_seqruns")



###########################################################
########  E. BLASTn taxonomy outputs: nt data base ########
###########################################################

## E.1 Load blastn results, format in table and save to csv   
load_blast("blast_results_full_nt_120924.tab", db_name = "nt_both_seqruns")


## E.2 Filter nt db blast results
filter_blast (raw_blastn = nt_both_seqruns_blastn_file_raw,
              min_perc_id = 95,
              min_query_cov = 95,
              abundance_table = both_seqruns_sequence_df_merged_nosg,
              db_name = "nt_both_seqruns",
              return_fasta_and_remaining = FALSE)


## E.3 Delete duplicates in filtered BLASTn results from nt data base
delete_duplicated_ID(filtered_blastn = nt_both_seqruns_filtered_blast, 
                     db_name = "nt_both_seqruns")


## E.4 Classify filtered BLASTn results from sponge data base
classify_blast(blastn_result = nt_both_seqruns_blastn_no_duplicates, 
               taxonomy_vector = c("kingdom", "phylum", "class", "order", "family", "genus", "species"), 
               db_name = "nt_both_seqruns")



########################################################################
########  F. Merge 2 BLASTn results together (sponge and nt db) ########
########################################################################

merge_taxo_abundance(taxonomyA = "sponge_both_seqruns_taxonomy_table.csv",
                     taxonomyB = "nt_both_seqruns_taxonomy_table.csv",
                     abundance_global =  "both_seqruns_merged_abundance_rarefied_nosg_noseq.csv",
                     filename = "sponge_nt")

