library(vegan)
library(tidyverse)

# For Kraken2 output, obtaining number of reads
# Note: Raw Kraken2 reads can't/shouldn't be used for estimates of abundance...
get_kraken2_taxids_and_reads <- function(folder_path, taxrank="S") {
  # Get a list of filenames
  filenames <- Sys.glob(paste0(folder_path, "/*.txt"))
  
  # Create output tibble
  all_reports <- tibble()
  
  # Iterate through the reports in the specified folder_path...
  for (f in filenames) {
    # Open the report file and set the colnames
    kraken_report_file <- readr::read_delim(f, col_names=FALSE)
    names(kraken_report_file) <- c("Reads Contained In Clade (%)", 
                                   "Reads Contained In Clade (#)", 
                                   "Reads At This Taxon (#)", 
                                   "Taxonomic Rank", 
                                   "NCBI TaxID", 
                                   "Scientific Name"
    )
    
    # Only looking at this taxonomic rank (usually species)
    filtered_file <- filter(kraken_report_file, `Taxonomic Rank` == taxrank)
    
    # If there's already at least one row in the report... (that means we don't need to add the taxids column)
    if (nrow(all_reports) != 0L) {
      # We only need these two columns from each of the report files
      subset_file <- filtered_file[, c("Reads At This Taxon (#)", "NCBI TaxID")]
      
      # Using regex to obtain the relevant sample/accession part of the name (i.e., it's from SRA)
      names(subset_file) <- c(sprintf("Reads (%s)", 
                                      str_extract(f, "SAMN[A-Za-z0-9]+(_[A-Za-z0-9]+)")), 
                              "NCBI TaxID"
      )
      
      # Join this to the main file by NCBI TaxID
      all_reports <- full_join(all_reports, subset_file, by="NCBI TaxID")
      
    } else { # No rows in all_reports yet!
      # We only need these two columns, and we're adding both straight into all_reports, since it's the first time. No need to join.
      all_reports <- filtered_file[, c("NCBI TaxID", "Reads At This Taxon (#)")]
      names(all_reports) <- c("NCBI TaxID", 
                              sprintf("Reads (%s)", 
                                      str_extract(f, "SAMN[A-Za-z0-9]+(_[A-Za-z0-9]+)")
                              )
      )
    }
    
  }
  return(all_reports)
}

# Exactly the same as the get_kraken2_taxids_and_reads, but for centrifuge
get_centrifuge_taxids_and_reads <- function(folder_path, taxrank="species", uniquereads=FALSE) {
  filenames <- Sys.glob(paste0(folder_path, "/*.report.txt"))
  
  all_reports <- tibble()
  
  for (f in filenames) {
    centrifuge_report_file <- readr::read_delim(f, col_names=TRUE)
    
    filtered_file <- filter(centrifuge_report_file, `taxRank` == taxrank)
    
    if (nrow(all_reports) != 0L) {
      subset_file <- filtered_file[, c(if(uniquereads) "numUniqueReads" else "numReads", "taxID")]
      
      names(subset_file) <- c(sprintf("Reads (%s)", 
                                      str_extract(f, "SAMN[A-Za-z0-9]+(_[A-Za-z0-9]+)")), 
                              "taxID"
      )
      
      all_reports <- full_join(all_reports, subset_file, by="taxID")
      
    } else {
      all_reports <- filtered_file[, c("taxID", if(uniquereads) "numUniqueReads" else "numReads")]
      names(all_reports) <- c("taxID", 
                              sprintf("Reads (%s)", 
                                      str_extract(f, "SAMN[A-Za-z0-9]+(_[A-Za-z0-9]+)")
                              )
      )
    }
  }
  return(all_reports)
}

# Same as above, but getting abundances instead of reads 
get_bracken_taxids_and_abundances <- function(folder_path, taxrank="S", all_taxa=FALSE, tsv=FALSE) {
  filenames <- Sys.glob(paste0(folder_path, (if (tsv) "/*.tsv" else "/*.txt")))
  
  all_reports <- tibble()
  
  for (f in filenames) {
    bracken_report_file <- readr::read_delim(f, col_names=FALSE)
    names(bracken_report_file) <- c("Name", 
                                   "NCBI TaxID", 
                                   "Taxonomic Rank", 
                                   "Kraken Assigned Reads", 
                                   "Bracken Added Reads", 
                                   "Bracken Total Reads",
                                   "Fraction of Total Reads"
    )
    
    if (!all_taxa) {
      filtered_file <- filter(bracken_report_file, `Taxonomic Rank` == taxrank)
    } else {
      filtered_file <- bracken_report_file
    }
    
    if (nrow(all_reports) != 0L) {
      subset_file <- filtered_file[, c("Fraction of Total Reads", "NCBI TaxID")]
      
      # Using regex to extract names of Raman samples (i.e., RAMAN_R003)
      names(subset_file) <- c(sprintf("Abundance (%s)", 
                                      str_extract(f, "R[A-Za-z0-9]+[_][A-Za-z0-9]+")), 
                              "NCBI TaxID"
      )
      
      all_reports <- full_join(all_reports, subset_file, by="NCBI TaxID")
      
    } else {
      all_reports <- filtered_file[, c("NCBI TaxID", "Fraction of Total Reads")]
      names(all_reports) <- c("NCBI TaxID", 
                              sprintf("Abundance (%s)", 
                                      str_extract(f, "R[A-Za-z0-9]+[_][A-Za-z0-9]+")
                              )
      )
    }
    
  }
  return(all_reports)
}

# Same as above, but getting abundances instead of reads
get_centrifuge_taxids_and_abundances <- function(folder_path, taxrank="species", all_taxa=FALSE) {
  filenames <- Sys.glob(paste0(folder_path, "/*.txt"))
  
  all_reports <- tibble()
  
  for (f in filenames) {
    centrifuge_report_file <- readr::read_delim(f, col_names=TRUE)
    
    if (!all_taxa) {
      filtered_file <- filter(centrifuge_report_file, `taxRank` == taxrank)
    } else {
      filtered_file <- centrifuge_report_file
    }
    
    if (nrow(all_reports) != 0L) {
      subset_file <- filtered_file[, c("abundance", "taxID")]
      
      # Using regex to extract names of Raman samples (i.e., RAMAN_R003)
      names(subset_file) <- c(sprintf("Abundance (%s)", 
                                      str_extract(f, "R[A-Za-z0-9]+[_][A-Za-z0-9]+")), 
                              "taxID"
      )
      
      all_reports <- full_join(all_reports, subset_file, by="taxID")
      
    } else {
      all_reports <- filtered_file[, c("taxID", "abundance")]
      names(all_reports) <- c("taxID", 
                              sprintf("Abundance (%s)", 
                                      str_extract(f, "R[A-Za-z0-9]+[_][A-Za-z0-9]+")
                              )
      )
    }
  }
  return(all_reports)
}

get_kreport_taxids_and_abundances <- function(folder_path, taxrank="S", all_taxa=FALSE) {
  filenames <- Sys.glob(paste0(folder_path, "/*.txt"))
  
  all_reports <- tibble()
  
  for (f in filenames) {
    kreport_file <- readr::read_delim(f, col_names=TRUE)
    
    if (!all_taxa) {
      filtered_file <- filter(kreport_file, kreport_file[c(4)] == taxrank)
    } else {
      filtered_file <- kreport_file
    }
    
    if (nrow(all_reports) != 0L) {
      subset_file <- filtered_file[, c(1, 5)]
      subset_file[, c(1)] = as.numeric(unlist(subset_file[, c(1)], use.names=FALSE)) / 100.0
      
      # Using regex to extract names of Raman samples (i.e., RAMAN_R003)
      names(subset_file) <- c(sprintf("Abundance (%s)", 
                                      str_extract(f, "R[A-Za-z0-9]+[_][A-Za-z0-9]+")), 
                              "taxID"
      )
      
      all_reports <- full_join(all_reports, subset_file, by="taxID")
      
    } else {
      all_reports <- filtered_file[, c(5, 1)]
      all_reports[, c(2)] = as.numeric(unlist(all_reports[, c(2)], use.names=FALSE)) / 100.0
      
      names(all_reports) <- c("taxID", 
                              sprintf("Abundance (%s)", 
                                      str_extract(f, "R[A-Za-z0-9]+[_][A-Za-z0-9]+")
                              )
      )
    }
  }
  return(all_reports)
}

# Does this need to be declared on its own? ...
gt_zero <- function(x) { return(x > 0 & !is.na(x)) }

# When joining from multiple report files and only including the species level, sub-species assignments can cause 0 abundance TaxIDs to be included
# We need to remove these! 
remove_absent_taxa <- function(dataset, abundance=TRUE) {
  filtered <- dataset %>% mutate(across(-c(1), as.numeric)) # Decimals will be interpreted as character values (?), make them numeric, except for TaxIDs
  filtered <- filter(filtered, if_any(contains(if (abundance) "Abundance" else "Reads"), gt_zero))
  filtered[is.na(filtered)] <- 0 # replacing NA with 0
  
  return(filtered)
}

# From one of these report structures returned by a get_x_taxids_and_y function, (and using remove_absent_taxa), 
# we can get a table of alpha diversity values per each sample (column) with abundance/counts(reads)
get_alpha_diversity_from_reports <- function(reports, abundance=TRUE) {
  # Create output tibble
  alpha_indices <- tibble(SampleName = character(), 
                          Shannon = numeric(), 
                          InverseSimpson = numeric())
  
  # Should I be using a loop here? ...
  for (col in colnames(reports)) {
    if (grepl(if (abundance) "Abundance" else "Reads", col)) { # Check that this isn't the "taxids" column
      shannon_index <- vegan::diversity(na.omit(reports[col]), index = "shannon")
      inverse_simpson_index <- vegan::diversity(na.omit(reports[col]), index = "invsimpson")
      
      alpha_indices <- add_row(alpha_indices, # Add a row to alpha_indices, containing SampleName, Shannon & InverseSimpson
                               SampleName=col, 
                               Shannon=shannon_index, 
                               InverseSimpson=inverse_simpson_index
      )
    }
  }
  
  return(alpha_indices)
  
}

# This function also adapts reports into BCI-style with t(select(filtered, -c(1))) for use with vegdist
get_rarefied_dataset <- function(reports) { 
  raremax <- min(rowSums(t(select(reports, -c(1))))) 
  rarefied <- rrarefy(t(select(reports, -c(1))), raremax) # This only allows count data... *what is rarefaction again?*
  return(rarefied)
}

# Short-hand, unused
get_kraken2_dataset <- function(filepath) {
  get_rarefied_dataset(
    remove_absent_taxa(
      get_kraken2_taxids_and_reads(filepath)
    )
  )
}

# Short-hand, unused
get_centrifuge_dataset <- function(filepath) {
  get_rarefied_dataset(
    remove_absent_taxa(
      get_centrifuge_taxids_and_reads(filepath)
    )
  )
}