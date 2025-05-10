library(tidyverse)
library(this.path)
library(reticulate)

get_list_of_run_accessions <- function() {
  run_accessions <- list.dirs(recursive = FALSE) %>%
    str_replace("./", "")
    
  return(run_accessions)
}

get_read_statistics <- function(run_accessions, folder_name) {
  reticulate::source_python(paste0(this.dir(), "/../python/generate_tsv_from_readqc_summaries.py"))
  generate_tsv_from_readqc_summaries(run_accessions, paste0(this.dir(), folder_name))
  
}

get_tsv_from_run_folder <- function(folder_name) {
  setwd(paste0(this.dir(), folder_name))
  read_tsv(get_read_statistics(get_list_of_run_accessions(), folder_name))
}
    
make_ggplot_scatterplot <- function(data, xx, yy) {
  ggplot(data, aes(x = xx, y = yy)) + 
    geom_point()
}



format_si <- function(...) {
  # Based on code by Ben Tupper
  # https://stat.ethz.ch/pipermail/r-help/2012-January/299804.html
  
  function(x) {
    limits <- c(1e-24, 1e-21, 1e-18, 1e-15, 1e-12,
                1e-9,  1e-6,  1e-3,  1e0,   1e3,
                1e6,   1e9,   1e12,  1e15,  1e18,
                1e21,  1e24)
    prefix <- c("y",   "z",   "a",   "f",   "p",
                "n",   "µ",   "m",   " ",   "k",
                "M",   "G",   "T",   "P",   "E",
                "Z",   "Y")
    
    # Vector with array indices according to position in intervals
    i <- findInterval(abs(x), limits)
    
    # Set prefix to " " for very small values < 1e-24
    i <- ifelse(i==0, which(limits == 1e0), i)
    
    paste(format(round(x/limits[i], 1),
                 trim=TRUE, scientific=FALSE, ...),
          prefix[i])
  }
}