library(tidyverse)
library(this.path)
library(reticulate)

get_list_of_run_accessions <- function() {
  run_accessions <- list.dirs(recursive = FALSE) %>%
    str_replace("./", "")
    
  return(run_accessions)
}

get_read_statistics <- function(run_accessions, run_name) {
  reticulate::source_python(paste0(this.dir(), "/../python/generate_tsv_from_readqc_summaries.py"))
  generate_tsv_from_readqc_summaries(run_accessions, paste0(this.dir(), paste0("/data/loma/", run_name)))
  
}

get_tsv_from_run_folder <- function(run_name) {
  setwd(paste0(this.dir(), paste0("/data/loma/", run_name)))
  read_tsv(get_read_statistics(get_list_of_run_accessions(), run_name))
}
    
make_ggplot_scatterplot <- function(data, xx, yy) {
  ggplot(data, aes(x = xx, y = yy)) + 
    geom_point()
}