#!/usr/bin/env Rscript

rm(list = ls())

## ------------------------------------------------------------------------

require("optparse")

## ------------------------------------------------------------------------
suppressPackageStartupMessages(library("optparse"))


## ------------------------------------------------------------------------

prepare_faa_PathoFact_orf <- function(aa_file_path = NULL,
                                      output_dir = "~/Desktop/",
                                      length_fil = 250,
                                      pat = "[[:space:]].*"){
  
  #######-------------------------- Loading required Packges
  require(Biostrings)
  require(tidyverse)
  require(Biostrings)
  
  #######-------------------------- Read faa
  
  aa_file_path %>% 
    Biostrings::readAAStringSet() -> sequences
  
  #######-------------------------- Simplify headers
  
  names(sequences) <- sub(pattern = pat, "", names(sequences))  # https://stackoverflow.com/questions/9319242/remove-everything-after-space-in-string
  
  #######--------------------------
  
  sequences[width(sequences) >= length_fil,] -> sequences
  
  #######--------------------------
  
  data.frame(names(sequences)) -> df 
  rownames(df) <- NULL
  
  #######--------------------------
  
  df %>% 
    mutate(CONTIG = sub("(.*?_.*?)_.*", "\\1", names.sequences.)) %>% 
    rename("ORF" = "names.sequences.") %>% 
    rownames_to_column("id") %>% 
    select(CONTIG, ORF) -> df
  
  colnames(df) <- NULL
  
  #######--------------------------
  
  
  write_tsv(x = df , 
            col_names = FALSE,
            file = paste0(output_dir,"/mapping.contig"))
  
  sequences %>% 
    Biostrings::writeXStringSet(paste0(output_dir,"/orf.faa"), append=FALSE,
                                compress=FALSE, compression_level=NA, format="fasta")
  
  #######--------------------------
  
  # system2("gsed ",
  #         args = c("-i ",paste(c(  's/\*//g' ))),
  #         args = c(" ",paste0(output_dir,"/orf.faa")),
  #         
  #         )
  
          system2("gsed",
                  args = c("-i 's/*//g' ", paste0(output_dir,"/orf.faa"))
  )

}



# prepare_faa_PathoFact_orf(aa_file_path = "/Volumes/Pegasus3 R4/NRP72/chicken/SQM/results/03.chicken.fna",
#                           length_fil = 10000)


## ------------------------------------------------------------------------

option_list = list(
  make_option(c("-f", "--aa_file_path"), type="character", default = NULL, 
              help="Path of the input aminoacid fasta of the orfs", metavar="character"),
  
  make_option(c("-o","--output_dir"), type="character", default = "~/Desktop/", 
              help="Name of the output directory -  created before running the analysis", metavar="character"),
  
  make_option(c("-l","--length_fil"), type="numeric", default= 250, 
              help="Length filtering of orf AA sequence", metavar="character"),
  
  make_option(c("--pat"), type="character", default= "[[:space:]].*", 
              help="Default regex pattern for faa header cleangin", metavar="character")
)

opt_parser = OptionParser(option_list = option_list);
opt = parse_args(opt_parser)

## ------------------------------------------------------------------------
# parse_args(opt_parser, args = c("--help"))

## ------------------------------------------------------------------------



prepare_faa_PathoFact_orf(aa_file_path = opt$aa_file_path,
                          output_dir = opt$output_dir,
                          length_fil = opt$length_fil,
                          pat = opt$pat)

## ------------------------------------------------------------------------


print(sessionInfo())

## ------------------------------------------------------------------------

print(Sys.time())

## ------------------------------------------------------------------------



