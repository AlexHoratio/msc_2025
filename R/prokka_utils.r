get_product_counts_from_prokka_tsv <- function(filepath) {
  raw <- read_tsv(filepath)
  
  product_frequencies <- tibble(dplyr::count(raw, product))
  colnames(product_frequencies) <- c("product", paste0(str_extract(sapply(strsplit(filepath, "/"), tail, 1), "[R]+[0-9]+"), 
                                                       " ",
                                                       str_extract(sapply(strsplit(filepath, "/"), tail, 1), "(?<=_000)[0-9]+")))
  
  return(product_frequencies)
  
}

get_gene_counts_from_all_bins <- function(filepath) {
  
  bin_folders <- Filter(function(x) !any(grepl("unbinned", x)), list.dirs(filepath, recursive = FALSE))
  gene_counts <- tibble(product = character())
  
  for (bin_folder in bin_folders) {
    gene_counts <- full_join(gene_counts, 
                         get_product_counts_from_prokka_tsv(paste0(bin_folder, "/", sapply(strsplit(bin_folder, "/"), tail, 1), ".tsv")),
                         by = "product"
    )
  }
  
  gene_counts <- gene_counts[!is.na(gene_counts$product), ]
  gene_counts[is.na(gene_counts)] <- 0
  return(gene_counts)
  
}

get_cog_counts_from_prokka_tsv <- function(filepath) {
  raw <- read_tsv(filepath)
  
  cog_frequencies <- tibble(dplyr::count(raw, COG))
  colnames(cog_frequencies) <- c("COG", paste0(str_extract(sapply(strsplit(filepath, "/"), tail, 1), "[R]+[0-9]+"), 
                                 " ",
                                 str_extract(sapply(strsplit(filepath, "/"), tail, 1), "(?<=_000)[0-9]+")))
  
  return(cog_frequencies)
  
}

get_cog_counts_from_all_bins <- function(filepath) {
  
  bin_folders <- Filter(function(x) !any(grepl("unbinned", x)), list.dirs(filepath, recursive = FALSE))
  cog_counts <- tibble(COG = character())
  
  for (bin_folder in bin_folders) {
    cog_counts <- full_join(cog_counts, 
                            get_cog_counts_from_prokka_tsv(paste0(bin_folder, "/", sapply(strsplit(bin_folder, "/"), tail, 1), ".tsv")),
                             by = "COG"
    )
  }
  
  cog_counts <- cog_counts[!is.na(cog_counts$COG), ]
  cog_counts[is.na(cog_counts)] <- 0
  return(cog_counts)
  
}

get_kegg_counts_from_prokka_tsv <- function(filepath) {
  raw <- read_tsv(filepath, col_names = FALSE)
  
  kegg_frequencies <- tibble(dplyr::count(raw, X2))
  colnames(kegg_frequencies) <- c("KO", paste0(str_extract(sapply(strsplit(filepath, "/"), tail, 1), "[R]+[0-9]+"), 
                                               " ",
                                               str_extract(sapply(strsplit(filepath, "/"), tail, 1), "(?<=_000)[0-9]+")))
  
  return(kegg_frequencies)
  
}

get_kegg_counts_from_all_bins <- function(filepath) {
  
  bin_folders <- Filter(function(x) !any(grepl("unbinned", x)), list.dirs(filepath, recursive = FALSE))
  kegg_counts <- tibble(`KO` = character())
  
  for (bin_folder in bin_folders) {
    kegg_counts <- full_join(kegg_counts,
                             get_kegg_counts_from_prokka_tsv(paste0(bin_folder, "/", sapply(strsplit(bin_folder, "/"), tail, 1), ".gbk.ko.out")),
                             by = "KO")
  }

  kegg_counts <- kegg_counts[!is.na(kegg_counts$`KO`), ]
  kegg_counts[is.na(kegg_counts)] <- 0
  return(kegg_counts)
  
}