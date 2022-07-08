---
title: "Merge chicken1 gene catalogue + quantification "
author: "Hannah Li Hägi & Florentin Constancias "
date: " July 04, 2022 "
output: 
  html_document: 
    toc: yes
    keep_md: yes
---





```r
rm(list= ls())
library('dplyr')
```

```
## 
## Attaching package: 'dplyr'
```

```
## The following objects are masked from 'package:stats':
## 
##     filter, lag
```

```
## The following objects are masked from 'package:base':
## 
##     intersect, setdiff, setequal, union
```

```r
library("tidyverse")
```

```
## ── Attaching packages ────────────────────────────────── tidyverse 1.3.1.9000 ──
```

```
## ✔ ggplot2 3.3.6     ✔ purrr   0.3.4
## ✔ tibble  3.1.7     ✔ stringr 1.4.0
## ✔ tidyr   1.2.0     ✔ forcats 0.5.1
## ✔ readr   2.1.2
```

```
## ── Conflicts ────────────────────────────────────────── tidyverse_conflicts() ──
## ✖ dplyr::filter() masks stats::filter()
## ✖ dplyr::lag()    masks stats::lag()
```

```r
library('here')
```

```
## here() starts at /Users/fconstan/Documents/GitHub/NRP72-FBT
```

```r
# library("SQMtools")
library("phyloseq")
library("readxl")


'%!in%' <- function(x,y)!('%in%'(x,y))

#source("https://raw.githubusercontent.com/fconstancias/DivComAnalyses/master/R/phyloseq_taxa_tests.R")
# source("https://raw.githubusercontent.com/fconstancias/DivComAnalyses/master/R/phyloseq_normalisation.R")
# source("https://raw.githubusercontent.com/fconstancias/DivComAnalyses/master/R/phyloseq_alpha.R")
# source("https://raw.githubusercontent.com/fconstancias/DivComAnalyses/master/R/phyloseq_beta.R")
# source("https://raw.githubusercontent.com/fconstancias/DivComAnalyses/master/R/phyloseq_heatmap.R")
```

# preliminary:

anvio no bin:


```r
# cd /Users/fconstan/Documents/GitHub/NRP72-FBT/data/processed/anvio/chicken1
# anvi-get-split-coverages -p PROFILE.db  -c CONTIGS.db --list-splits  > all_splits.tsv
# "/Users/fconstan/Documents/GitHub/NRP72-FBT/data/processed/anvio/chicken1/all_splits.tsv" %>%
#   read_tsv(col_names = c("split", "bin")) %>%
#   mutate(bin = "all") %>%
#   data.frame() ->  Manual_bins_tmp
# 
# Manual_bins_tmp %>%
#   write_tsv(here::here("data/processed/anvio/chicken1/collection_all_splits.tsv"), col_names = FALSE)


# anvi-import-collection -p PROFILE.db  -c CONTIGS.db  -C all_splits collection_all_splits.tsv
# anvi-summarize  -p PROFILE.db  -c CONTIGS.db  --list-collections
# anvi-summarize  -p PROFILE.db  -c CONTIGS.db  -C all_splits -o summary_no_bins

#anvi-summarize  -p PROFILE.db  -c CONTIGS.db  -C concoct_10_raw -o summary_no_bins # **actually this is al the splits : NOT SURE**
```

* FYI: A subset of split sequences are being initialized (to be precise, only
118590 of 665668 splits the contigs database knows about). Nothing to worry
about. Probably.


# Define inputs:

## tables/ data:


```r
SQMR <- "data/raw/SQM/chicken1_full_SQM.rds"

metadata <- "data/raw/25.11.2021_metadata_updated.tsv"

RGI_CARD <- "data/raw/CARD/sqm_rgi_chicken1.txt"
Viralverify <- "data/raw/mobilome/viralverify/01.chicken1_result_table.csv"
isescan <- "data/raw/mobilome/isescan/01.SqueezeChicken1.fasta.tsv"
PathoFact_AMR <- "data/processed/PathoFact/chicken1/PathoFact_chicken1_AMR/PathoFact_results_orf/PathoFact_report/AMR_MGE_prediction_chicken_l1000_report.tsv"
PathoFact_Tox <- "data/processed/PathoFact/chicken1/PathoFact_chicken1_Tox/PathoFact_results_orf/PathoFact_report/Toxin_prediction_chicken_l1000_report.tsv"
PathoFact_Tox_fam <- "data/processed/PathoFact/chicken1/PathoFact_chicken1_Tox/PathoFact_results_orf/PathoFact_report/Toxin_gene_library_chicken_l1000_report.tsv" 

anvio_nobin_summary <- "data/processed/anvio/chicken1/summary_no_bins/bin_by_bin/all/all-gene_calls.txt"

DAS_collection <- "data/processed/anvio/chicken1/DAS.txt"
Manual_bins <- "data/processed/anvio/chicken1/MAGS/concoct_20_collection_refined_2.txt"

ANVIO_SUMMARY_DAS  <- "data/processed/anvio/chicken1/MAGS/SUMMARY/DAS/bins_summary.txt"
ANVIO_SUMMARY_CONCOCT <- "data/processed/anvio/chicken1/MAGS/SUMMARY/concoct_20_refined_2/bins_summary.txt"

GTDB_tax <- "data/processed/gtdb/GTDB_all.xlsx"
SQM_DAS <- "data/raw/SQM/19.chicken1.bintable"

plasx <- "data/processed/mobilome/chicken1_contig_plasx-scores.txt"

# contigs.lgt.tsv <- "~/Documents/MSc_Thesis-NRP72/data/HGT/human.contigs.lgt.tsv"
# contigs.no_lgt.tsv <- "~/Documents/MSc_Thesis-NRP72/data/HGT/human.contigs.no_lgt.tsv"
# contigs.unclassified.tsv <- "~/Documents/MSc_Thesis-NRP72/data/HGT/human.contigs.unclassified.tsv"
```


```r
gc()
```

```
##           used  (Mb) gc trigger  (Mb) limit (Mb) max used  (Mb)
## Ncells 5058379 270.2    9162152 489.4         NA  6824151 364.5
## Vcells 8521871  65.1   14786712 112.9    1024000 10536360  80.4
```

## parameters:


```r
Model_type <- c("protein homolog model", "protein variant model", "protein overexpression model", "rRNA gene variant model")
Cut_Off <- c("Strict", "Perfect")
besid = 95
Percent_length = 0.80
```

# Gene/annotation - contig classification:

## Load SQM:

```r
SQMR %>% 
  here::here() %>% 
  readRDS() -> human
```

## Load metadata:


```r
metadata %>% 
  here::here() %>% 
  read_tsv() %>% 
  data.frame() %>%
  dplyr::filter(!is.na(metagenomic_sample_name), Model == "Chicken") -> human_meta
```

```
## Rows: 604 Columns: 61
## ── Column specification ────────────────────────────────────────────────────────
## Delimiter: "\t"
## chr (19): sample, Sample_description, I7_Index_ID, index, I5_Index_ID, index...
## dbl (42): input, filtered, denoisedF, denoisedR, merged, tabled, filtered_pc...
## 
## ℹ Use `spec()` to retrieve the full column specification for this data.
## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.
```

Modify metadata: **this should be done once on the metadata**


```r
# create Period column
human_meta %>%
  mutate(Reactor_Treatment = base::replace(Reactor_Treatment, Reactor_Treatment == "TR_2CTX", "TR2_CTX")) %>%
  mutate(Reactor_Treatment_Dose = if_else(!is.na(`Antibiotic_mg.mL`), paste0(Reactor_Treatment, `Antibiotic_mg.mL`), Reactor_Treatment)) %>% 
  mutate(Treatment_Dose = if_else(!is.na(`Antibiotic_mg.mL`), paste0(Treatment, `Antibiotic_mg.mL`), Treatment)) %>%
  mutate(
    Period = case_when(
      Day_of_Treatment <= 0  ~ "pret",
      Day_of_Treatment > 0 ~ as.character(Day_of_Treatment)
    )) %>%
  mutate(Period = base::replace(Period, Treatment == "DONOR", "pret")) %>%
  mutate(Day_of_Treatment = addNA(Day_of_Treatment)) %>%
  column_to_rownames("metagenomic_sample_name") -> human_meta

#human_meta$Treatment_Dose <- factor(human_meta$Treatment_Dose, ordered = TRUE, 
#                                 levels = c("DONOR", "UNTREATED", "CTX+HV292.120", "CTX+HV292.1200",
#                                            "CTX20", "CTX200", "HV292.1", "VAN+CCUG5916890", "VAN+CCUG59168600",
#                                            "VAN90", "VAN600", "CCUG59168"))
```


## PathoFact:

### PathoFact_AMR:


```r
PathoFact_AMR %>% 
  here::here() %>% 
  read_tsv() %>% 
  select(-Contig, -Contig_ID, -ORF_ID) %>% 
  na_if("-") %>% 
  na_if("n/a") -> PathoFact_AMR
```

```
## Rows: 1021675 Columns: 11
## ── Column specification ────────────────────────────────────────────────────────
## Delimiter: "\t"
## chr (11): Contig, Contig_ID, ORF, ORF_ID, ARG, ARG_SNPs, AMR_category, AMR_s...
## 
## ℹ Use `spec()` to retrieve the full column specification for this data.
## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.
```

```r
colnames(PathoFact_AMR) <- paste0("PathoFact_AMR_", colnames(PathoFact_AMR))

PathoFact_AMR %>% 
  rename(ORF_ID = PathoFact_AMR_ORF) -> PathoFact_AMR

PathoFact_AMR %>% 
  dim()
```

```
## [1] 1021675       8
```

### PathoFact_Tox:


```r
PathoFact_Tox %>% 
  here::here() %>% 
  read_tsv() %>% 
  select(-ORF_ID, -Number_of_hits, -Signal_peptide, -Toxin_prediction) %>% 
  na_if("-") %>% 
  na_if("n/a") %>% 
  filter(!is.na(Toxin_confidence_level)) -> PathoFact_Tox
```

```
## Rows: 1021675 Columns: 6
## ── Column specification ────────────────────────────────────────────────────────
## Delimiter: "\t"
## chr (5): ORF, ORF_ID, Toxin_prediction, Signal_peptide, Toxin_confidence_level
## dbl (1): Number_of_hits
## 
## ℹ Use `spec()` to retrieve the full column specification for this data.
## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.
```

```r
PathoFact_Tox %>% 
  dim()
```

```
## [1] 13796     2
```


```r
PathoFact_Tox_fam %>% 
  here::here() %>% 
  read_tsv() %>% 
  # add_count(ORF, name = "count_ORF") %>% arrange(-count_ORF)
  group_by(ORF) %>% 
  top_n(1, wt = `Significance_evalue`) %>% 
  top_n(1, wt = NAME) %>% #add_count(ORF, name = "count_ORF") %>% arrange(-count_ORF)
  select(-ORF_ID, -Database, -Alternative_name) -> PathoFact_Tox_fam
```

```
## Rows: 16057 Columns: 9
## ── Column specification ────────────────────────────────────────────────────────
## Delimiter: "\t"
## chr (7): ORF_ID, ORF, HMM_Name, NAME, Alternative_name, Database, Description
## dbl (2): Score, Significance_evalue
## 
## ℹ Use `spec()` to retrieve the full column specification for this data.
## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.
```

```r
PathoFact_Tox_fam %>% 
  dim()
```

```
## [1] 13796     6
```


```r
PathoFact_Tox %>% 
  left_join(PathoFact_Tox_fam,
            by = c("ORF" = "ORF")) -> PathoFact_Tox


PathoFact_Tox %>% 
  dim()
```

```
## [1] 13796     7
```



```r
colnames(PathoFact_Tox) <- paste0("PathoFact_Tox_", colnames(PathoFact_Tox))

PathoFact_Tox %>% 
  rename(PathoFact_Tox_Toxin_classification = PathoFact_Tox_Toxin_confidence_level,
         ORF_ID = PathoFact_Tox_ORF) -> PathoFact_Tox


PathoFact_Tox %>% 
  head() %>% 
  DT::datatable()
```

```{=html}
<div id="htmlwidget-547ddae7bf76666607c9" style="width:100%;height:auto;" class="datatables html-widget"></div>
<script type="application/json" data-for="htmlwidget-547ddae7bf76666607c9">{"x":{"filter":"none","vertical":false,"data":[["1","2","3","4","5","6"],["megahit_397_2-1306","megahit_409_917-1549","megahit_513_1584-2006","megahit_717_1252-1569","megahit_717_11990-12355","megahit_931_11981-12394"],["2: Non-secreted Toxin","2: Non-secreted Toxin","2: Non-secreted Toxin","2: Non-secreted Toxin","2: Non-secreted Toxin","2: Non-secreted Toxin"],["PF01636.18","PF00717.18","TIGR01593","PF08845.5","PF01325.14","TIGR01593"],[45,53,82.3,42.3,43.8,93.5],[7.1e-12,1.3e-14,1.5e-23,3.1e-11,1.2e-11,4.9e-27],["APH","Peptidase_S24","holin_tox_secr","SymE_toxin","Fe_dep_repress","holin_tox_secr"],["Phosphotransferase enzyme family","Peptidase S24-like","toxin secretion/phage lysis holin","Toxin SymE, type I toxin-antitoxin system","Iron dependent repressor, N-terminal DNA binding domain","toxin secretion/phage lysis holin"]],"container":"<table class=\"display\">\n  <thead>\n    <tr>\n      <th> <\/th>\n      <th>ORF_ID<\/th>\n      <th>PathoFact_Tox_Toxin_classification<\/th>\n      <th>PathoFact_Tox_HMM_Name<\/th>\n      <th>PathoFact_Tox_Score<\/th>\n      <th>PathoFact_Tox_Significance_evalue<\/th>\n      <th>PathoFact_Tox_NAME<\/th>\n      <th>PathoFact_Tox_Description<\/th>\n    <\/tr>\n  <\/thead>\n<\/table>","options":{"columnDefs":[{"className":"dt-right","targets":[4,5]},{"orderable":false,"targets":0}],"order":[],"autoWidth":false,"orderClasses":false}},"evals":[],"jsHooks":[]}</script>
```
## rgi CARD:


```r
# Resistome data:
RGI_CARD %>% 
  here::here() %>% 
  read_tsv() %>% 
  # read.delim() %>%
  select(-`Contig ID`:-Hits) %>%
  dplyr::filter(Model_type %in% Model_type ,
                Cut_Off %in% Cut_Off,
                Nudged %!in% c("TRUE"),
                Best_Identities > besid, # plot distribution of bestID vs percentage length ref and color bitscore and symbol class AMR
                `Percentage.Length.of.Reference.Sequence` > Percent_length) %>%
  separate(ORF_ID, c("Megahit", "Contig", "Split"), sep = "_", fill = "right", remove = FALSE) %>%
  # mutate(ORF_ID = paste(Megahit, Contig, paste(Start, Stop, sep = "-"), sep = "_")) %>%
  #mutate(Contig = paste(Megahit, Contig, sep = "_")) %>%
  select(-c(Contig, Split, CARD_Protein_Sequence, Predicted_DNA, Predicted_Protein, Contig)) %>% 
  na_if("NA") %>% 
  na_if("n/a") %>% 
  na_if("<NA>") %>% 
  na_if("") %>% 
  select(-Start, -Stop, -Orientation) %>% 
  select(ORF_ID, everything()) -> card
```

```
## Rows: 154 Columns: 212
## ── Column specification ────────────────────────────────────────────────────────
## Delimiter: "\t"
## chr  (25): ORF_ID, Cut_Off, Best_Hit_ARO, Model_type, SNPs_in_Best_Hit_ARO, ...
## dbl (185): Start, Stop, Orientation, Pass_Bitscore, Best_Hit_Bitscore, Best_...
## lgl   (2): Predicted_DNA, Nudged
## 
## ℹ Use `spec()` to retrieve the full column specification for this data.
## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.
```

```r
colnames(card) <- paste0("rgi_CARD_", colnames(card))

card %>% 
  rename(ORF_ID = rgi_CARD_ORF_ID) -> card

# card %>% 
#   # head() %>% 
#   DT::datatable()
```

## Waffle HGT: 

Not here since as to be on a per sample basis


```r
# TODO move to more appropriate section not 

# # HGT events
# contigs.lgt.tsv %>% 
#   read_tsv() %>% 
#   select(CONTIG_NAME, CALL) -> hgt
# 
# 
# contigs.no_lgt.tsv %>% 
#   read_tsv() %>% 
#   select(CONTIG_NAME, CALL) -> no_hgt
# 
# contigs.unclassified.tsv %>% 
#   read_tsv() %>% 
#   select(CONTIG_NAME, CALL) -> unclassified_hgt
# 
# hgt_events <- rbind(hgt, no_hgt, unclassified_hgt)
```

## Viralverify:


```r
# Viralverify 
Viralverify  %>%
  read_csv() %>% 
  rename(Contig_ID = `Contig name`,
         viralverify_Prediction = Prediction) %>% 
  select(Contig_ID, viralverify_Prediction) -> viralverify
```

```
## Warning: One or more parsing issues, see `problems()` for details
```

```
## Rows: 655896 Columns: 6
## ── Column specification ────────────────────────────────────────────────────────
## Delimiter: ","
## chr (5): Contig name, Prediction, Circular, Score, Pfam hits
## dbl (1): Length
## 
## ℹ Use `spec()` to retrieve the full column specification for this data.
## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.
```


# isescan:


```r
isescan %>%
  read_tsv() %>% 
  group_by(seqID) %>%
  add_count(cluster, name = "count_IS_contig") %>% 
  top_n(1, wt = `E-value4copy`) %>%
  top_n(1, wt = isBegin) %>%
  mutate(IS = "yes") %>%
  select(seqID:cluster, type, count_IS_contig) -> top_isescan
```

```
## Rows: 12361 Columns: 24
## ── Column specification ────────────────────────────────────────────────────────
## Delimiter: "\t"
## chr  (6): seqID, family, cluster, strand, type, tir
## dbl (18): isBegin, isEnd, isLen, ncopy4is, start1, end1, start2, end2, score...
## 
## ℹ Use `spec()` to retrieve the full column specification for this data.
## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.
```

```r
colnames(top_isescan) <- paste0("isescan_", colnames(top_isescan))

top_isescan %>% 
  rename(Contig_ID = isescan_seqID) -> top_isescan


top_isescan %>% 
  head() %>% 
  DT::datatable()
```

```{=html}
<div id="htmlwidget-9819bb7c456398aa30a6" style="width:100%;height:auto;" class="datatables html-widget"></div>
<script type="application/json" data-for="htmlwidget-9819bb7c456398aa30a6">{"x":{"filter":"none","vertical":false,"data":[["1","2","3","4","5","6"],["megahit_100015","megahit_100038","megahit_100108","megahit_100161","megahit_100224","megahit_100245"],["IS256","IS66","ISL3","IS30","IS21","IS3"],["IS256_128","IS66_43","ISL3_126|ISL3||protein:plasmid:139408","IS30_173","IS21_35","IS3_255"],["p","p","p","p","p","c"],[1,1,1,1,1,1]],"container":"<table class=\"display\">\n  <thead>\n    <tr>\n      <th> <\/th>\n      <th>Contig_ID<\/th>\n      <th>isescan_family<\/th>\n      <th>isescan_cluster<\/th>\n      <th>isescan_type<\/th>\n      <th>isescan_count_IS_contig<\/th>\n    <\/tr>\n  <\/thead>\n<\/table>","options":{"columnDefs":[{"className":"dt-right","targets":5},{"orderable":false,"targets":0}],"order":[],"autoWidth":false,"orderClasses":false}},"evals":[],"jsHooks":[]}</script>
```

## plasx:



```r
plasx %>% 
  here::here() %>% 
  read_tsv(col_names = c("Contig_ID", "plasX_score")) -> pasx_score
```

```
## Rows: 304575 Columns: 2
## ── Column specification ────────────────────────────────────────────────────────
## Delimiter: "\t"
## chr (2): Contig_ID, plasX_score
## 
## ℹ Use `spec()` to retrieve the full column specification for this data.
## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.
```


## Combine all:


```r
human$orfs$table %>%
  rownames_to_column("ORF_ID") %>% 
  data.frame() %>% 
  na_if("") %>% 
  # left_join(
  # human$orfs$tax %>% 
  #   as.data.frame() %>% 
  #   rownames_to_column("ORF_ID") %>%
  #   separate(ORF_ID, c("Megahit", "Contig", "Range"), sep = "_", fill = "right", remove = FALSE) %>%
  #   # mutate(Contig = as.numeric(Contig)) %>%
  #   mutate(Contig = paste(Megahit, Contig, sep = "_")) %>% 
  #   left_join(human$orfs$seqs %>%
  #     as.data.frame() %>%
  #   rownames_to_column("ORF_ID")) %>% 
  #   select(-Megahit) %>% 
#   rename(seq = ".",
#          Contig_ID = Contig,
#          superkingdom_ORF = superkingdom,
#          phylum_ORF = phylum,
#          class_ORF= class,
#          order_ORF = order,
#          family_orf = family,
#          genus_orf = genus,
#          species_orf = species)
# ) %>% 
select(-Method, -Molecule, -Hits) %>% 
  rename(Contig_ID = `Contig.ID`,
         Tax_orf = Tax) -> orf_tmp


orf_tmp %>% 
  select(-(starts_with(c("TPM.","Raw.base.count.","Coverage.", "Raw.read.count")))) -> orf
# starts_with(c("Petal", "Sepal"))
# orf %>% 
#   write_tsv("~/Desktop/unitedhuman1.tsv")
colnames(orf) <- paste0("SQM_", colnames(orf))

orf %>% 
  rename(Contig_ID = SQM_Contig_ID,
         ORF_ID = SQM_ORF_ID) -> orf

# orf %>%
#   head() %>%
#   DT::datatable()
```


```r
human$contigs$table %>%
  rownames_to_column("Contig_ID") %>% 
  data.frame() %>% 
  rename(Tax_contig = `Tax`,
         Disparity_contig = Disparity,
         GC_perc_contig = GC.perc,
         Length_contig = Length,
         Num_genes_contig = Num.genes) %>% 
  select(-(starts_with(c("TPM.","Raw.base.count.","Coverage.", "Raw.read.count")))) %>% 
  select(-Bin.ID) -> contigs


colnames(contigs) <- paste0("SQM_", colnames(contigs))

contigs %>% 
  rename(Contig_ID = SQM_Contig_ID) -> contigs

contigs %>% 
  head() %>% 
  DT::datatable()
```

```{=html}
<div id="htmlwidget-c612669f3aebf3d57f52" style="width:100%;height:auto;" class="datatables html-widget"></div>
<script type="application/json" data-for="htmlwidget-c612669f3aebf3d57f52">{"x":{"filter":"none","vertical":false,"data":[["1","2","3","4","5","6"],["megahit_1","megahit_2","megahit_3","megahit_4","megahit_5","megahit_6"],["k_Bacteria;n_Terrabacteria group;p_Firmicutes;c_Clostridia;o_Eubacteriales","","k_Bacteria;n_Terrabacteria group;p_Firmicutes;c_Clostridia;o_Eubacteriales;n_Eubacteriales incertae sedis;g_Intestinimonas;s_Intestinimonas timonensis","k_Bacteria;n_Terrabacteria group;p_Firmicutes;c_Clostridia;o_Eubacteriales","k_Bacteria;n_Terrabacteria group;p_Firmicutes;c_Clostridia;o_Eubacteriales;f_Oscillospiraceae",""],[0,0,0,0,0,0],[61.17,57.69,66.67,45.79,65.56,35.22],[703,234,234,214,270,318],[1,1,1,1,1,1]],"container":"<table class=\"display\">\n  <thead>\n    <tr>\n      <th> <\/th>\n      <th>Contig_ID<\/th>\n      <th>SQM_Tax_contig<\/th>\n      <th>SQM_Disparity_contig<\/th>\n      <th>SQM_GC_perc_contig<\/th>\n      <th>SQM_Length_contig<\/th>\n      <th>SQM_Num_genes_contig<\/th>\n    <\/tr>\n  <\/thead>\n<\/table>","options":{"columnDefs":[{"className":"dt-right","targets":[3,4,5,6]},{"orderable":false,"targets":0}],"order":[],"autoWidth":false,"orderClasses":false}},"evals":[],"jsHooks":[]}</script>
```


```r
orf %>% 
  dim()
```

```
## [1] 1634985      14
```


```r
gc()
```

```
##              used   (Mb) gc trigger    (Mb) limit (Mb)   max used    (Mb)
## Ncells   10151326  542.2   18356477   980.4         NA   14294612   763.5
## Vcells 1303406573 9944.3 2199626536 16781.9    1024000 1832478970 13980.8
```


```r
orf %>% 
  left_join(card,
            by = c("ORF_ID" = "ORF_ID")) %>% 
  left_join(PathoFact_AMR,
            by = c("ORF_ID" = "ORF_ID")) %>% 
  left_join(pasx_score,
            by = c("Contig_ID" = "Contig_ID")) %>% 
  left_join(PathoFact_Tox,
            by = c("ORF_ID" = "ORF_ID")) %>% 
  left_join(viralverify,
            by = c("Contig_ID" = "Contig_ID")) %>% 
  left_join(top_isescan,
            by = c("Contig_ID" = "Contig_ID")) %>% 
  left_join(contigs,
            by = c("Contig_ID" = "Contig_ID")) -> orf_contig

orf_contig %>% 
  dim()
```

```
## [1] 1634985      56
```

```r
orf_contig %>% 
  head() %>% 
  DT::datatable()
```

```{=html}
<div id="htmlwidget-be4cba94955b298b4a6b" style="width:100%;height:auto;" class="datatables html-widget"></div>
<script type="application/json" data-for="htmlwidget-be4cba94955b298b4a6b">{"x":{"filter":"none","vertical":false,"data":[["1","2","3","4","5","6"],["megahit_1_2-643","megahit_2_2-232","megahit_3_3-233","megahit_4_2-214","megahit_5_1-270","megahit_6_3-317"],["megahit_1","megahit_2","megahit_3","megahit_4","megahit_5","megahit_6"],[642,231,231,213,270,315],[214,77,77,71,90,105],[62.46,57.58,67.1,45.54,65.56,35.24],[null,null,"ABC.CD.P","ABCB-BAC",null,null],["k_Bacteria;n_Terrabacteria group;p_Firmicutes;c_Clostridia;o_Eubacteriales",null,"k_Bacteria;n_Terrabacteria group;p_Firmicutes;c_Clostridia;o_Eubacteriales;n_Eubacteriales incertae sedis;g_Intestinimonas;s_Intestinimonas timonensis","k_Bacteria;n_Terrabacteria group;p_Firmicutes;c_Clostridia;o_Eubacteriales","k_Bacteria;n_Terrabacteria group;p_Firmicutes;c_Clostridia;o_Eubacteriales;f_Oscillospiraceae",null],[null,null,"K02004*","K06147*",null,null],[null,null,"putative ABC transport system permease protein","ATP-binding cassette, subfamily B, bacterial",null,null],[null,null,"Brite Hierarchies; Protein families: signaling and cellular processes; Transporters","Brite Hierarchies; Protein families: signaling and cellular processes; Transporters",null,null],["ENOG4111FY8*",null,"COG0577*","COG1132*","ENOG4111FI4",null],[null,null,"ABC-type antimicrobial peptide transport system, permease component","ABC-type multidrug transport system, ATPase and permease components","domain protein",null],["Function unknown",null,null,null,"Function unknown",null],[null,null,null,null,null,null],[null,null,null,null,null,null],[null,null,null,null,null,null],[null,null,null,null,null,null],[null,null,null,null,null,null],[null,null,null,null,null,null],[null,null,null,null,null,null],[null,null,null,null,null,null],[null,null,null,null,null,null],[null,null,null,null,null,null],[null,null,null,null,null,null],[null,null,null,null,null,null],[null,null,null,null,null,null],[null,null,null,null,null,null],[null,null,null,null,null,null],[null,null,null,null,null,null],[null,null,null,null,null,null],[null,null,null,null,null,null],[null,null,null,null,null,null],[null,null,null,null,null,null],[null,null,null,null,null,null],[null,null,null,null,null,null],[null,null,null,null,null,null],[null,null,null,null,null,null],[null,null,null,null,null,null],[null,null,null,null,null,null],[null,null,null,null,null,null],[null,null,null,null,null,null],[null,null,null,null,null,null],[null,null,null,null,null,null],[null,null,null,null,null,null],[null,null,null,null,null,null],[null,null,null,null,null,null],["Uncertain - too short","Uncertain - too short","Uncertain - too short","Uncertain - too short","Uncertain - too short","Uncertain - too short"],[null,null,null,null,null,null],[null,null,null,null,null,null],[null,null,null,null,null,null],[null,null,null,null,null,null],["k_Bacteria;n_Terrabacteria group;p_Firmicutes;c_Clostridia;o_Eubacteriales","","k_Bacteria;n_Terrabacteria group;p_Firmicutes;c_Clostridia;o_Eubacteriales;n_Eubacteriales incertae sedis;g_Intestinimonas;s_Intestinimonas timonensis","k_Bacteria;n_Terrabacteria group;p_Firmicutes;c_Clostridia;o_Eubacteriales","k_Bacteria;n_Terrabacteria group;p_Firmicutes;c_Clostridia;o_Eubacteriales;f_Oscillospiraceae",""],[0,0,0,0,0,0],[61.17,57.69,66.67,45.79,65.56,35.22],[703,234,234,214,270,318],[1,1,1,1,1,1]],"container":"<table class=\"display\">\n  <thead>\n    <tr>\n      <th> <\/th>\n      <th>ORF_ID<\/th>\n      <th>Contig_ID<\/th>\n      <th>SQM_Length.NT<\/th>\n      <th>SQM_Length.AA<\/th>\n      <th>SQM_GC.perc<\/th>\n      <th>SQM_Gene.name<\/th>\n      <th>SQM_Tax_orf<\/th>\n      <th>SQM_KEGG.ID<\/th>\n      <th>SQM_KEGGFUN<\/th>\n      <th>SQM_KEGGPATH<\/th>\n      <th>SQM_COG.ID<\/th>\n      <th>SQM_COGFUN<\/th>\n      <th>SQM_COGPATH<\/th>\n      <th>SQM_PFAM<\/th>\n      <th>rgi_CARD_Megahit<\/th>\n      <th>rgi_CARD_Cut_Off<\/th>\n      <th>rgi_CARD_Pass_Bitscore<\/th>\n      <th>rgi_CARD_Best_Hit_Bitscore<\/th>\n      <th>rgi_CARD_Best_Hit_ARO<\/th>\n      <th>rgi_CARD_Best_Identities<\/th>\n      <th>rgi_CARD_ARO<\/th>\n      <th>rgi_CARD_Model_type<\/th>\n      <th>rgi_CARD_SNPs_in_Best_Hit_ARO<\/th>\n      <th>rgi_CARD_Other_SNPs<\/th>\n      <th>rgi_CARD_Drug.Class<\/th>\n      <th>rgi_CARD_Resistance.Mechanism<\/th>\n      <th>rgi_CARD_AMR.Gene.Family<\/th>\n      <th>rgi_CARD_Percentage.Length.of.Reference.Sequence<\/th>\n      <th>rgi_CARD_ID<\/th>\n      <th>rgi_CARD_Model_ID<\/th>\n      <th>rgi_CARD_Nudged<\/th>\n      <th>rgi_CARD_Note<\/th>\n      <th>PathoFact_AMR_ARG<\/th>\n      <th>PathoFact_AMR_ARG_SNPs<\/th>\n      <th>PathoFact_AMR_AMR_category<\/th>\n      <th>PathoFact_AMR_AMR_sub_class<\/th>\n      <th>PathoFact_AMR_Resistance_mechanism<\/th>\n      <th>PathoFact_AMR_Database<\/th>\n      <th>PathoFact_AMR_MGE_prediction<\/th>\n      <th>plasX_score<\/th>\n      <th>PathoFact_Tox_Toxin_classification<\/th>\n      <th>PathoFact_Tox_HMM_Name<\/th>\n      <th>PathoFact_Tox_Score<\/th>\n      <th>PathoFact_Tox_Significance_evalue<\/th>\n      <th>PathoFact_Tox_NAME<\/th>\n      <th>PathoFact_Tox_Description<\/th>\n      <th>viralverify_Prediction<\/th>\n      <th>isescan_family<\/th>\n      <th>isescan_cluster<\/th>\n      <th>isescan_type<\/th>\n      <th>isescan_count_IS_contig<\/th>\n      <th>SQM_Tax_contig<\/th>\n      <th>SQM_Disparity_contig<\/th>\n      <th>SQM_GC_perc_contig<\/th>\n      <th>SQM_Length_contig<\/th>\n      <th>SQM_Num_genes_contig<\/th>\n    <\/tr>\n  <\/thead>\n<\/table>","options":{"columnDefs":[{"className":"dt-right","targets":[3,4,5,17,18,20,21,28,30,43,44,51,53,54,55,56]},{"orderable":false,"targets":0}],"order":[],"autoWidth":false,"orderClasses":false}},"evals":[],"jsHooks":[]}</script>
```



```r
orf %>% 
filter(ORF_ID == "c1_megahit_288931_1-1971")
```

```
##  [1] ORF_ID        Contig_ID     SQM_Length.NT SQM_Length.AA SQM_GC.perc  
##  [6] SQM_Gene.name SQM_Tax_orf   SQM_KEGG.ID   SQM_KEGGFUN   SQM_KEGGPATH 
## [11] SQM_COG.ID    SQM_COGFUN    SQM_COGPATH   SQM_PFAM     
## <0 rows> (or 0-length row.names)
```



```r
orf_contig$ORF_ID %>%  unique() %>%  length()
```

```
## [1] 1634985
```


## Test anvio no bin:

** TODO: make sure all splits are there... **


```r
anvio_nobin_summary %>% 
  here::here() %>%  
  read_tsv() -> anvi_no_bin_annot
```

```
## Rows: 922861 Columns: 22
## ── Column specification ────────────────────────────────────────────────────────
## Delimiter: "\t"
## chr (19): contig, direction, COGPATH, COGPATH (ACCESSION), COG, COG (ACCESSI...
## dbl  (3): gene_callers_id, start, stop
## 
## ℹ Use `spec()` to retrieve the full column specification for this data.
## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.
```

Link anvio gene id with SQM gene id:


```r
anvi_no_bin_annot %>% 
  mutate(orf_id = paste0(contig, "_", start+1, "-", stop)) %>% 
  select(-gene_callers_id) -> anvi_no_bin_annot

colnames(anvi_no_bin_annot) <- paste0("anvio_", colnames(anvi_no_bin_annot))

anvi_no_bin_annot %>% 
  select(anvio_orf_id, everything()) %>% 
  select(-anvio_contig, -anvio_start, -anvio_stop, -anvio_direction ,-anvio_dna_sequence) -> anvi_no_bin_annot
```



```r
orf_contig %>% 
  left_join(anvi_no_bin_annot,
            by = c("ORF_ID" = "anvio_orf_id")) -> orf_contig
```



```r
# orf_contig %>%
#   filter(anvio_contig == "megahit_37")
```


```r
orf_contig %>% 
  dim()
```

```
## [1] 1634986      72
```

## MAG informations:

### Collections:

#### SQM DAS collection:


```r
DAS_collection %>% 
  here::here() %>% 
  read_tsv(col_names = c("split", "DAS_bin")) %>% 
  separate(split, into = c("cont", "cont_num", "split", "split_id"), sep = "_", remove = FALSE) %>% 
  mutate(Contig_ID = paste0(cont,"_", cont_num)) %>% 
  distinct(Contig_ID, .keep_all = TRUE) %>% 
  select(Contig_ID, DAS_bin) -> DAS
```

```
## Rows: 70127 Columns: 2
## ── Column specification ────────────────────────────────────────────────────────
## Delimiter: "\t"
## chr (2): split, DAS_bin
## 
## ℹ Use `spec()` to retrieve the full column specification for this data.
## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.
```

```r
DAS %>% 
  head() %>% 
  DT::datatable()
```

```{=html}
<div id="htmlwidget-e24e09e5f9c714ceb242" style="width:100%;height:auto;" class="datatables html-widget"></div>
<script type="application/json" data-for="htmlwidget-e24e09e5f9c714ceb242">{"x":{"filter":"none","vertical":false,"data":[["1","2","3","4","5","6"],["megahit_15787","megahit_378627","megahit_127613","megahit_292113","megahit_193238","megahit_132560"],["maxbin_006_fasta_contigs","maxbin_006_fasta_contigs","maxbin_006_fasta_contigs","maxbin_006_fasta_contigs","maxbin_006_fasta_contigs","maxbin_006_fasta_contigs"]],"container":"<table class=\"display\">\n  <thead>\n    <tr>\n      <th> <\/th>\n      <th>Contig_ID<\/th>\n      <th>DAS_bin<\/th>\n    <\/tr>\n  <\/thead>\n<\/table>","options":{"columnDefs":[{"orderable":false,"targets":0}],"order":[],"autoWidth":false,"orderClasses":false}},"evals":[],"jsHooks":[]}</script>
```

#### Manual MAGs CONCOCT:


```r
Manual_bins %>% 
  read_tsv(col_names = c("split", "CONCOCT_bin")) %>% 
  separate(split, into = c("cont", "cont_num", "split", "split_id"), sep = "_", remove = FALSE) %>% 
  mutate(Contig_ID = paste0(cont,"_", cont_num)) %>% 
  distinct(Contig_ID, .keep_all = TRUE) %>% 
  select(Contig_ID, CONCOCT_bin) -> CONCOCT
```

```
## Rows: 41235 Columns: 2
## ── Column specification ────────────────────────────────────────────────────────
## Delimiter: "\t"
## chr (2): split, CONCOCT_bin
## 
## ℹ Use `spec()` to retrieve the full column specification for this data.
## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.
```

```r
CONCOCT  %>% 
  head() %>% 
  DT::datatable()
```

```{=html}
<div id="htmlwidget-3eb8429d48a8b7ccda82" style="width:100%;height:auto;" class="datatables html-widget"></div>
<script type="application/json" data-for="htmlwidget-3eb8429d48a8b7ccda82">{"x":{"filter":"none","vertical":false,"data":[["1","2","3","4","5","6"],["megahit_419123","megahit_10003","megahit_616572","megahit_241507","megahit_329011","megahit_419949"],["Bin_19","Bin_19","Bin_19","Bin_19","Bin_19","Bin_19"]],"container":"<table class=\"display\">\n  <thead>\n    <tr>\n      <th> <\/th>\n      <th>Contig_ID<\/th>\n      <th>CONCOCT_bin<\/th>\n    <\/tr>\n  <\/thead>\n<\/table>","options":{"columnDefs":[{"orderable":false,"targets":0}],"order":[],"autoWidth":false,"orderClasses":false}},"evals":[],"jsHooks":[]}</script>
```

Generate fake bin for all contigs -> export kegg annotations from anvio?


```r
# Manual_bins %>%
#   read_tsv(col_names = c("split", "bin")) %>%
#   mutate(bin = "all") %>% 
#   data.frame() ->  Manual_bins_tmp
# 
# # names(Manual_bins_tmp) <- NULL
# 
# Manual_bins_tmp %>% 
#   write_tsv(here::here("data/processed/anvio/human1/collection_all_splits.tsv"), col_names = FALSE)
```

Let's add the contig <-> bin to the merged data:


```r
orf_contig %>% 
  left_join(CONCOCT) %>% 
  left_join(DAS) -> orf_contig_2
```

```
## Joining, by = "Contig_ID"
## Joining, by = "Contig_ID"
```



```r
orf_contig_2 %>% dim()
```

```
## [1] 1634986      74
```

### MAGs summary:

#### MAGs GDTB taxonomy:


```r
GTDB_tax %>%
  read_excel() %>%
  filter(dataset == "chicken1") %>%
  mutate(user_genome = str_remove_all(user_genome, c("chicken1_|-contigs|refined_2_"))) %>%
  filter(str_detect(user_genome, c("DAS|concoct"))) %>%
  select(user_genome, classification) %>%
  # separate(classification, 
  #          into = c("domain", "phylum", "class", "order", "family", "genus", "species"), 
  #          sep = ";", 
  #          fill = "right") %>% 
  rename(GTDB_tax = classification) %>% 
  mutate(user_genome = str_remove_all(user_genome, "concoct_20_")) %>% 
  mutate(user_genome = str_remove_all(user_genome, "DAS_")) -> GTDB_tax
```

```
## New names:
## • `` -> `...23`
## • `` -> `...24`
## • `` -> `...25`
## • `` -> `...26`
## • `` -> `...27`
## • `` -> `...28`
```

```r
# TODO clean taxonomy and add highest resolution.

GTDB_tax %>% 
  head() %>% 
  DT::datatable()
```

```{=html}
<div id="htmlwidget-d50616a2b63e9727dfcc" style="width:100%;height:auto;" class="datatables html-widget"></div>
<script type="application/json" data-for="htmlwidget-d50616a2b63e9727dfcc">{"x":{"filter":"none","vertical":false,"data":[["1","2","3","4","5","6"],["maxbin_002_fasta_contigs","maxbin_003_fasta_contigs","maxbin_004_fasta_contigs","maxbin_006_fasta_contigs","maxbin_007_fasta_contigs","maxbin_008_fasta_contigs"],["d__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;o__Enterobacterales;f__Enterobacteriaceae;g__Proteus;s__Proteus mirabilis","d__Bacteria;p__Firmicutes_A;c__Clostridia;o__Oscillospirales;f__Ruminococcaceae;g__Gemmiger;s__","d__Bacteria;p__Firmicutes_A;c__Clostridia;o__Oscillospirales;f__Acutalibacteraceae;g__;s__","d__Bacteria;p__Firmicutes_A;c__Clostridia;o__Lachnospirales;f__Lachnospiraceae;g__Mediterraneibacter;s__Mediterraneibacter surreyensis","d__Bacteria;p__Firmicutes_A;c__Clostridia;o__Lachnospirales;f__Lachnospiraceae;g__Eisenbergiella;s__Eisenbergiella sp904392525","d__Bacteria;p__Firmicutes_A;c__Clostridia;o__Oscillospirales;f__Oscillospiraceae;g__Flavonifractor;s__Flavonifractor avistercoris"]],"container":"<table class=\"display\">\n  <thead>\n    <tr>\n      <th> <\/th>\n      <th>user_genome<\/th>\n      <th>GTDB_tax<\/th>\n    <\/tr>\n  <\/thead>\n<\/table>","options":{"columnDefs":[{"orderable":false,"targets":0}],"order":[],"autoWidth":false,"orderClasses":false}},"evals":[],"jsHooks":[]}</script>
```


```r
orf_contig_2 %>% 
  left_join(GTDB_tax,
            by = c("CONCOCT_bin" = "user_genome")) %>%  #,
  # suffix = c("", "CONCOCT")) %>% 
  rename(GTDB_tax_CONCOCT = GTDB_tax) %>% 
  left_join(GTDB_tax,
            by = c("DAS_bin" = "user_genome")) %>% #,
  # suffix = c("", "_DAS")) %>% 
  rename(GTDB_tax_DAS = GTDB_tax) -> orf_contig_3
```




```r
orf_contig_3 %>% 
  dim()
```

```
## [1] 1634986      76
```


```r
## CONCOCT tax info

# GTDB_tax %>%
#   filter(str_detect(user_genome, "concoct")) %>%
#   mutate(user_genome = str_remove_all(user_genome, "concoct_10_")) -> CONCOCT_tax
# 
# colnames(CONCOCT_tax) <- paste0("CONCOCT_", colnames(CONCOCT_tax))
```

#### SQM DAS info:


```r
SQM_DAS %>% 
  read_tsv(skip = 1) -> SQM_DAS
```

```
## Rows: 237 Columns: 109
## ── Column specification ────────────────────────────────────────────────────────
## Delimiter: "\t"
## chr   (4): Bin ID, Method, Tax, Tax 16S
## dbl (105): Length, GC perc, Num contigs, Disparity, Completeness, Contaminat...
## 
## ℹ Use `spec()` to retrieve the full column specification for this data.
## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.
```

```r
colnames(SQM_DAS) <- paste0("SQM_DAS_", colnames(SQM_DAS))

SQM_DAS %>%
  mutate(SQM_DAS_HQ = ifelse(SQM_DAS_Completeness >= 80 & SQM_DAS_Contamination <= 10, "HQ", NA)) %>% 
  select(-SQM_DAS_Method, -`SQM_DAS_Coverage C1`:-`SQM_DAS_TPM D9`) -> SQM_DAS_HQ

# filter(SQM_DAS_Completeness >= 80,
#        SQM_DAS_Contamination <= 10) %>% 
# mutate(SQM_DAS_bins_HQ = `SQM_DAS_Bin ID`) -> SQM_DAS_HQ

SQM_DAS_HQ %>% 
  head() %>% 
  DT::datatable()
```

```{=html}
<div id="htmlwidget-c4934cb9d78eded89e50" style="width:100%;height:auto;" class="datatables html-widget"></div>
<script type="application/json" data-for="htmlwidget-c4934cb9d78eded89e50">{"x":{"filter":"none","vertical":false,"data":[["1","2","3","4","5","6"],["maxbin.070.fasta.contigs","maxbin.148.fasta.contigs","maxbin.127.fasta.contigs","maxbin.073.fasta.contigs","maxbin.053.fasta.contigs","maxbin.159.fasta.contigs"],["k_Bacteria;p_Proteobacteria;c_Gammaproteobacteria;o_Pseudomonadales;f_Moraxellaceae;g_Acinetobacter","k_Bacteria;n_Terrabacteria group;p_Firmicutes;c_Clostridia;o_Eubacteriales;f_Lachnospiraceae","k_Bacteria;n_Terrabacteria group;p_Firmicutes;c_Bacilli;o_Bacillales;f_Bacillaceae","k_Bacteria;n_Terrabacteria group;p_Firmicutes;c_Clostridia;o_Eubacteriales","k_Bacteria;n_Terrabacteria group;p_Firmicutes;c_Bacilli;o_Lactobacillales;f_Enterococcaceae;g_Enterococcus","k_Bacteria;n_Terrabacteria group;p_Firmicutes;c_Clostridia;o_Eubacteriales"],[null,null,null,null,"superkingdom:Bacteria;clade:Terrabacteria group;phylum:Firmicutes;class:Bacilli;order:Lactobacillales;family:Enterococcaceae;genus:Enterococcus","superkingdom:Bacteria;clade:Terrabacteria group;phylum:Firmicutes;class:Erysipelotrichia;order:Erysipelotrichales;family:Erysipelotrichaceae;genus:Holdemania"],[3921451,2552034,5170540,3847753,3875630,3762815],[38.85,43.69,37.5,41.55,41.72,44.37],[81,95,41,338,108,268],[0,0.474,0,0.065,0,0.495],[100,99.38,99.35,99.09,98.46,98.18],[0,0.91,0.32,2.4,1.24,30.45],[0,0,0,0,0,0],["HQ","HQ","HQ","HQ","HQ",null]],"container":"<table class=\"display\">\n  <thead>\n    <tr>\n      <th> <\/th>\n      <th>SQM_DAS_Bin ID<\/th>\n      <th>SQM_DAS_Tax<\/th>\n      <th>SQM_DAS_Tax 16S<\/th>\n      <th>SQM_DAS_Length<\/th>\n      <th>SQM_DAS_GC perc<\/th>\n      <th>SQM_DAS_Num contigs<\/th>\n      <th>SQM_DAS_Disparity<\/th>\n      <th>SQM_DAS_Completeness<\/th>\n      <th>SQM_DAS_Contamination<\/th>\n      <th>SQM_DAS_Strain het<\/th>\n      <th>SQM_DAS_HQ<\/th>\n    <\/tr>\n  <\/thead>\n<\/table>","options":{"columnDefs":[{"className":"dt-right","targets":[4,5,6,7,8,9,10]},{"orderable":false,"targets":0}],"order":[],"autoWidth":false,"orderClasses":false}},"evals":[],"jsHooks":[]}</script>
```


```r
orf_contig_3 %>% 
  left_join(SQM_DAS_HQ,
            by = c("DAS_bin" = "SQM_DAS_Bin ID")) -> orf_contig_3
```




```r
orf_contig_3 %>% 
  dim()
```

```
## [1] 1634986      86
```


#### ANVIO DAS info:


```r
ANVIO_SUMMARY_DAS %>% 
  here::here() %>% 
  read_tsv() %>% 
  select(-t_domain:-t_species) -> ANVIO_SUMMARY_DAS
```

```
## Rows: 239 Columns: 14
## ── Column specification ────────────────────────────────────────────────────────
## Delimiter: "\t"
## chr (8): bins, t_domain, t_phylum, t_class, t_order, t_family, t_genus, t_sp...
## dbl (6): total_length, num_contigs, N50, GC_content, percent_completion, per...
## 
## ℹ Use `spec()` to retrieve the full column specification for this data.
## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.
```

```r
colnames(ANVIO_SUMMARY_DAS) <- paste0("ANVIO_SUMMARY_DAS_", colnames(ANVIO_SUMMARY_DAS))

ANVIO_SUMMARY_DAS %>%
  # filter(ANVIO_SUMMARY_DAS_percent_completion >= 80,
  #        ANVIO_SUMMARY_DAS_percent_redundancy <= 10) %>% 
  mutate(ANVIO_DAS_HQ = ifelse(ANVIO_SUMMARY_DAS_percent_completion >= 80 & ANVIO_SUMMARY_DAS_percent_redundancy <= 10, "HQ", NA)) -> ANVIO_SUMMARY_DAS

ANVIO_SUMMARY_DAS %>% 
  head() %>% 
  DT::datatable()
```

```{=html}
<div id="htmlwidget-cc84c92f3aa619bf0b2c" style="width:100%;height:auto;" class="datatables html-widget"></div>
<script type="application/json" data-for="htmlwidget-cc84c92f3aa619bf0b2c">{"x":{"filter":"none","vertical":false,"data":[["1","2","3","4","5","6"],["maxbin_002_fasta_contigs","maxbin_003_fasta_contigs","maxbin_004_fasta_contigs","maxbin_006_fasta_contigs","maxbin_007_fasta_contigs","maxbin_008_fasta_contigs"],[3814813,1841718,1514456,3791763,4867209,2589054],[40,98,108,295,1047,343],[219559,48209,47857,24520,5939,11944],[38.8077037673219,57.8443069423649,55.2371287942733,45.3010169138919,52.0581629153053,63.2587714220528],[95.7746478873239,64.7887323943662,53.5211267605634,88.7323943661972,83.0985915492958,64.7887323943662],[0,2.8169014084507,25.3521126760563,22.5352112676056,11.2676056338028,7.04225352112676],["HQ",null,null,null,null,null]],"container":"<table class=\"display\">\n  <thead>\n    <tr>\n      <th> <\/th>\n      <th>ANVIO_SUMMARY_DAS_bins<\/th>\n      <th>ANVIO_SUMMARY_DAS_total_length<\/th>\n      <th>ANVIO_SUMMARY_DAS_num_contigs<\/th>\n      <th>ANVIO_SUMMARY_DAS_N50<\/th>\n      <th>ANVIO_SUMMARY_DAS_GC_content<\/th>\n      <th>ANVIO_SUMMARY_DAS_percent_completion<\/th>\n      <th>ANVIO_SUMMARY_DAS_percent_redundancy<\/th>\n      <th>ANVIO_DAS_HQ<\/th>\n    <\/tr>\n  <\/thead>\n<\/table>","options":{"columnDefs":[{"className":"dt-right","targets":[2,3,4,5,6,7]},{"orderable":false,"targets":0}],"order":[],"autoWidth":false,"orderClasses":false}},"evals":[],"jsHooks":[]}</script>
```


```r
orf_contig_3 %>% 
  left_join(ANVIO_SUMMARY_DAS,
            by = c("DAS_bin" = "ANVIO_SUMMARY_DAS_bins")) -> orf_contig_3
```



```r
orf_contig_3 %>% 
  dim()
```

```
## [1] 1634986      93
```

#### ANVIO CONCOCT info:


```r
ANVIO_SUMMARY_CONCOCT %>% 
  here::here() %>% 
  read_tsv() %>% 
  select(-t_domain:-t_species) -> ANVIO_SUMMARY_CONCOCT
```

```
## Rows: 180 Columns: 14
## ── Column specification ────────────────────────────────────────────────────────
## Delimiter: "\t"
## chr (8): bins, t_domain, t_phylum, t_class, t_order, t_family, t_genus, t_sp...
## dbl (6): total_length, num_contigs, N50, GC_content, percent_completion, per...
## 
## ℹ Use `spec()` to retrieve the full column specification for this data.
## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.
```

```r
colnames(ANVIO_SUMMARY_CONCOCT) <- paste0("ANVIO_SUMMARY_CONCOCT_", colnames(ANVIO_SUMMARY_CONCOCT))

ANVIO_SUMMARY_CONCOCT %>% 
  mutate(ANVIO_CONCOCT_HQ = ifelse(ANVIO_SUMMARY_CONCOCT_percent_completion >= 80 & ANVIO_SUMMARY_CONCOCT_percent_redundancy <= 10, "HQ", NA)) -> ANVIO_SUMMARY_CONCOCT

ANVIO_SUMMARY_CONCOCT %>% 
  head() %>% 
  DT::datatable()
```

```{=html}
<div id="htmlwidget-279f104a2c5a187bbfef" style="width:100%;height:auto;" class="datatables html-widget"></div>
<script type="application/json" data-for="htmlwidget-279f104a2c5a187bbfef">{"x":{"filter":"none","vertical":false,"data":[["1","2","3","4","5","6"],["Bin_0_1","Bin_0_10","Bin_0_11_1","Bin_0_12","Bin_0_13_1","Bin_0_14"],[1520044,1911713,1460748,1595220,1358898,850979],[301,186,352,145,293,15],[6217,20284,4881,19367,5999,83079],[46.8891267269538,40.3301233911154,40.1860333299844,49.0252164980562,51.2929135824418,46.2562250852464],[57.7464788732394,64.7887323943662,40.8450704225352,69.0140845070423,80.2816901408451,73.2394366197183],[0,2.8169014084507,0,4.22535211267606,5.63380281690141,0],[null,null,null,null,"HQ",null]],"container":"<table class=\"display\">\n  <thead>\n    <tr>\n      <th> <\/th>\n      <th>ANVIO_SUMMARY_CONCOCT_bins<\/th>\n      <th>ANVIO_SUMMARY_CONCOCT_total_length<\/th>\n      <th>ANVIO_SUMMARY_CONCOCT_num_contigs<\/th>\n      <th>ANVIO_SUMMARY_CONCOCT_N50<\/th>\n      <th>ANVIO_SUMMARY_CONCOCT_GC_content<\/th>\n      <th>ANVIO_SUMMARY_CONCOCT_percent_completion<\/th>\n      <th>ANVIO_SUMMARY_CONCOCT_percent_redundancy<\/th>\n      <th>ANVIO_CONCOCT_HQ<\/th>\n    <\/tr>\n  <\/thead>\n<\/table>","options":{"columnDefs":[{"className":"dt-right","targets":[2,3,4,5,6,7]},{"orderable":false,"targets":0}],"order":[],"autoWidth":false,"orderClasses":false}},"evals":[],"jsHooks":[]}</script>
```


```r
orf_contig_3 %>% 
  left_join(ANVIO_SUMMARY_CONCOCT,
            by = c("CONCOCT_bin" = "ANVIO_SUMMARY_CONCOCT_bins")) -> orf_contig_3
```



```r
orf_contig_3 %>% 
  dim()
```

```
## [1] 1634986     100
```


Rename ORF and contig ID:

```r
orf_contig_3 %>% 
  mutate(ORF_ID = paste0("c1_", ORF_ID),
         Contig_ID = paste0("c1_", Contig_ID)) -> full_gene_catalogue
```


```r
full_gene_catalogue %>% 
  dim()
```

```
## [1] 1634986     100
```


```r
full_gene_catalogue %>% 
  filter(PathoFact_AMR_MGE_prediction %in% c("plasmid", "ambiguous (plasmid/phage)") &
         viralverify_Prediction %in% c("Uncertain - too short", "Plamid", "Uncertain - plasmid or chromosomal", "Plasmid")  ) %>% 
  filter(!is.na(PathoFact_AMR_AMR_category)) %>% 
  # filter(PathoFact_AMR_MGE_prediction %in% c("plasmid", "ambiguous (plasmid/phage)", "phage","ambiguous (phage/chromosome)" )) %>% 
  select(ORF_ID, PathoFact_AMR_ARG, PathoFact_AMR_AMR_category, viralverify_Prediction,PathoFact_AMR_MGE_prediction,  DAS_bin, CONCOCT_bin) %>% 
  DT::datatable()
```

```{=html}
<div id="htmlwidget-00a255519b448206e59c" style="width:100%;height:auto;" class="datatables html-widget"></div>
<script type="application/json" data-for="htmlwidget-00a255519b448206e59c">{"x":{"filter":"none","vertical":false,"data":[["1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","23","24","25","26","27","28","29","30","31","32","33","34","35","36","37","38","39","40","41","42","43","44","45","46","47","48","49","50","51","52","53","54","55","56","57","58","59","60","61","62","63","64","65","66","67","68","69","70","71","72","73","74","75","76","77","78","79","80","81","82","83","84","85","86","87","88","89","90","91","92","93","94","95","96","97","98","99","100","101","102","103","104","105","106","107","108","109","110","111","112","113","114","115","116","117","118","119","120","121","122","123","124","125","126","127","128","129","130","131","132","133","134","135","136","137","138","139","140","141","142","143","144","145","146","147","148","149","150","151","152","153","154","155","156","157","158","159","160","161","162","163","164","165","166","167","168","169","170","171","172","173","174"],["c1_megahit_8430_2303-3106","c1_megahit_8430_3106-3942","c1_megahit_8430_5473-5703","c1_megahit_8588_1391-1876","c1_megahit_8588_2029-2940","c1_megahit_8588_3368-3976","c1_megahit_8588_3982-5031","c1_megahit_8588_5006-5974","c1_megahit_8588_6189-7283","c1_megahit_8588_7321-8016","c1_megahit_10204_188-2107","c1_megahit_11002_71-394","c1_megahit_12280_2617-4068","c1_megahit_14914_2020-2133","c1_megahit_17599_12821-14053","c1_megahit_18358_8285-8398","c1_megahit_18803_1071-1277","c1_megahit_18868_1040-1723","c1_megahit_20276_139-834","c1_megahit_22207_5667-5885","c1_megahit_30555_3989-4132","c1_megahit_37545_1-171","c1_megahit_40932_3-611","c1_megahit_41600_1-918","c1_megahit_43752_2-91","c1_megahit_45556_1348-1455","c1_megahit_54590_6068-6703","c1_megahit_58580_1394-2353","c1_megahit_60380_10265-11125","c1_megahit_60380_11298-11975","c1_megahit_60380_12054-13253","c1_megahit_60643_2-256","c1_megahit_64732_1474-1578","c1_megahit_67476_1459-2223","c1_megahit_70146_139-366","c1_megahit_74379_105-761","c1_megahit_79788_808-1362","c1_megahit_80088_1085-1951","c1_megahit_87279_500-1603","c1_megahit_91130_225-917","c1_megahit_98297_3-179","c1_megahit_98933_90-1013","c1_megahit_101532_660-1295","c1_megahit_102253_554-2023","c1_megahit_104551_1-1299","c1_megahit_105454_7182-7289","c1_megahit_116776_14411-15577","c1_megahit_117067_410-1669","c1_megahit_119989_10581-11762","c1_megahit_122070_812-1024","c1_megahit_124494_13049-14317","c1_megahit_135054_134-1642","c1_megahit_136014_24218-25384","c1_megahit_137297_14786-15517","c1_megahit_137476_222-1145","c1_megahit_141381_2475-2690","c1_megahit_151168_18269-18499","c1_megahit_153155_1-333","c1_megahit_161298_3879-5078","c1_megahit_163311_1988-2239","c1_megahit_167979_2432-2644","c1_megahit_174128_118-1413","c1_megahit_175927_21699-21911","c1_megahit_176772_1460-1675","c1_megahit_181505_5667-6833","c1_megahit_184034_1011-1679","c1_megahit_184775_414-2333","c1_megahit_186553_1074-1151","c1_megahit_186897_7208-7681","c1_megahit_187286_3-89","c1_megahit_191247_1-225","c1_megahit_200600_3930-4745","c1_megahit_201652_3-74","c1_megahit_202645_7984-9882","c1_megahit_205635_1046-1969","c1_megahit_216648_16117-17367","c1_megahit_219391_905-1780","c1_megahit_224719_19033-19662","c1_megahit_225943_1162-1254","c1_megahit_235183_1-606","c1_megahit_235608_949-1041","c1_megahit_237605_537-950","c1_megahit_238278_3-68","c1_megahit_238953_2-409","c1_megahit_238953_479-1258","c1_megahit_238953_1763-2602","c1_megahit_245533_429-1904","c1_megahit_251792_261-1283","c1_megahit_254009_1156-1311","c1_megahit_254931_1-138","c1_megahit_256419_471-620","c1_megahit_264130_473-1396","c1_megahit_267234_414-827","c1_megahit_271698_426-1124","c1_megahit_273399_2-97","c1_megahit_276820_1280-1684","c1_megahit_281638_1086-1190","c1_megahit_294721_460-663","c1_megahit_301431_23904-24581","c1_megahit_305865_410-1033","c1_megahit_309104_1-1308","c1_megahit_315556_1176-1364","c1_megahit_316981_326-1198","c1_megahit_319106_20147-21328","c1_megahit_322550_2-109","c1_megahit_323384_26-1948","c1_megahit_323713_1-117","c1_megahit_324420_869-1606","c1_megahit_327597_5086-5745","c1_megahit_330281_173-631","c1_megahit_334788_592-1008","c1_megahit_340072_1951-2598","c1_megahit_350807_734-2110","c1_megahit_352506_6849-7061","c1_megahit_355453_32716-32823","c1_megahit_358432_204-842","c1_megahit_365035_3428-4111","c1_megahit_371143_410-1141","c1_megahit_374902_2-82","c1_megahit_379002_77-982","c1_megahit_381648_370-579","c1_megahit_382689_267-1139","c1_megahit_382853_3090-3299","c1_megahit_389506_582-794","c1_megahit_399404_12215-13465","c1_megahit_403241_107-370","c1_megahit_406489_987-1682","c1_megahit_409032_1159-1254","c1_megahit_416398_1337-2050","c1_megahit_425747_1-573","c1_megahit_428405_1443-2651","c1_megahit_436738_1072-1992","c1_megahit_439738_3678-4097","c1_megahit_439856_94-357","c1_megahit_440044_760-2718","c1_megahit_440044_2702-3964","c1_megahit_452389_1928-2140","c1_megahit_479601_6781-8250","c1_megahit_482790_259-357","c1_megahit_483940_298-1086","c1_megahit_483940_1144-1668","c1_megahit_483940_1763-2236","c1_megahit_488803_844-1473","c1_megahit_492570_417-635","c1_megahit_502839_198-1652","c1_megahit_503848_262-1956","c1_megahit_509017_678-1250","c1_megahit_514616_876-1562","c1_megahit_523767_1947-2039","c1_megahit_527694_1-633","c1_megahit_527694_764-1552","c1_megahit_531380_87-1334","c1_megahit_531977_1139-1348","c1_megahit_534004_1484-1597","c1_megahit_535112_1742-1972","c1_megahit_539047_2-190","c1_megahit_555380_1689-1985","c1_megahit_557529_223-1023","c1_megahit_563050_740-1366","c1_megahit_585117_727-1674","c1_megahit_589255_79-1506","c1_megahit_589426_79-1506","c1_megahit_592620_346-1029","c1_megahit_602341_4969-5625","c1_megahit_603840_278-1075","c1_megahit_603840_1823-2752","c1_megahit_615480_443-1120","c1_megahit_615646_1594-1737","c1_megahit_623315_1-972","c1_megahit_634750_948-1166","c1_megahit_636435_294-539","c1_megahit_646934_1623-1748","c1_megahit_649618_662-1300","c1_megahit_655046_4400-4495"],["APH(3'')-IB","APH(6)-ID","VGAC","VANZA","VANYA","VANXA","VANA","VANHA","VANSA","VANRA","TET32","QACG","BIFUNCTIONAL_AMINOGLYCOSIDE_N-ACETYLTRANSFERASE_AND_AMINOGLYCOSIDE_PHOSPHOTRANSFERASE","TET(44)","UGD","TET(44)","VANU","VANR","VANR","VANU","OQXB","ADEC","VANR","BCRA","PSEUDOMONAS AERUGINOSA CATB7","RPHB","VATB","UGD","TEM-1","TETR","TET(A)","OXA-97","OXA-45","ANT(4')-IA","PDC-91","QNRS1","VATB","LLMA_23S_RIBOSOMAL_RNA_METHYLTRANSFERASE","VANI","VANR","RGT1438","BCRA","MEXL","POXTA","LSA","CTX-M-12","UGD","UGD","MACA","VANU","UGD","LSAE","UGD","ERMA","BCRA","VANU","VGAC","TRUNCATED_PUTATIVE_RESPONSE_REGULATOR_ARLR","UGD","VANU","VANU","MAJOR_FACILITATOR_SUPERFAMILY_TRANSPORTER","VANU","VANU","UGD","VATB","TET(W/N/W)","OXA-297","DFRF","TETB(58)","AAD(6)","SUL2","APH(9)-IB","TET(W/N/W)","BCRA","UGD","CTX-M-1","TETM","QNRS3","VATB","MEXN","FOSX","ERMA","ANT(2'')-IA","AADA2","SUL1","LSAE","UGD","BACA","LNUC","ADEK","BCRA","VANS","VANR","AXYY","DNA-BINDING_PROTEIN_H-NS","ANT(6)-IA","VANU","TRUNCATED_PUTATIVE_RESPONSE_REGULATOR_ARLR","VATE","EPTA","APH(6)-ID","LLMA_23S_RIBOSOMAL_RNA_METHYLTRANSFERASE","UGD","STAPHYLOCOCCUS MUPA CONFERRING RESISTANCE TO MUPIROCIN","TET(44)","TET(44)","ERMB","VATB","DFRA5","ADP-RIBOSYLATING_TRANSFERASE_ARR","VATE","TET(L)","VANU","RPHB","VATE","VANR","LLMA_23S_RIBOSOMAL_RNA_METHYLTRANSFERASE","MTRD","MPHB","VANU","AADE","VANU","VANU","UGD","VANU","VANR","BRUCELLA SUIS MPRF","VANR","TRUNCATED_PUTATIVE_RESPONSE_REGULATOR_ARLR","BIFUNCTIONAL_AMINOGLYCOSIDE_N-ACETYLTRANSFERASE_AND_AMINOGLYCOSIDE_PHOSPHOTRANSFERASE","BCRA","FOSB","VANU","TETB(P)","TETP","VANU","LSA","MCR-7.1","ANT(3'')-IIA","SAT-2","DFRA1","VATB","VANU","LSAE","POXTA","VATB","VANR","MOX-3","DFRA17","AADA5","UGD","VANU","VANRG","VGAC","THIN-B","YKKD","ERMZ","VATB","CLASS_A","TET(W/N/W)","TET(W/N/W)","TRUNCATED_PUTATIVE_RESPONSE_REGULATOR_ARLR","VATB","BACA","BCRA","TRUNCATED_PUTATIVE_RESPONSE_REGULATOR_ARLR","TETA(58)","MDTM","VANU","VANU","OCH-3","VATE","APH(3')-IVA"],["aminoglycoside","aminoglycoside","multidrug","glycopeptide","glycopeptide","glycopeptide","glycopeptide","glycopeptide","glycopeptide","glycopeptide","tetracycline","multidrug","aminoglycoside","tetracycline","peptide","tetracycline","glycopeptide","glycopeptide","glycopeptide","glycopeptide","multidrug","multidrug","glycopeptide","bacitracin","phenicol","rifamycin","MLS","peptide","beta-lactam","tetracycline","tetracycline","beta-lactam","beta-lactam","aminoglycoside","beta-lactam","fluoroquinolone","MLS","MLS","glycopeptide","glycopeptide","rifamycin","bacitracin","multidrug","multidrug","multidrug","beta-lactam","peptide","peptide","MLS","glycopeptide","peptide","multidrug","peptide","MLS","bacitracin","glycopeptide","multidrug","unclassified","peptide","glycopeptide","glycopeptide","multidrug","glycopeptide","glycopeptide","peptide","MLS","tetracycline","beta-lactam","diaminopyrimidine","tetracycline","aminoglycoside","sulfonamide","aminoglycoside","tetracycline","bacitracin","peptide","beta-lactam","tetracycline","fluoroquinolone","MLS","phenicol","fosfomycin","MLS","aminoglycoside","aminoglycoside","sulfonamide","multidrug","peptide","bacitracin","MLS","multidrug","bacitracin","glycopeptide","glycopeptide","multidrug","unclassified","aminoglycoside","glycopeptide","unclassified","MLS","peptide","aminoglycoside","MLS","peptide","mupirocin","tetracycline","tetracycline","MLS","MLS","diaminopyrimidine","rifamycin","MLS","tetracycline","glycopeptide","rifamycin","MLS","glycopeptide","MLS","multidrug","MLS","glycopeptide","aminoglycoside","glycopeptide","glycopeptide","peptide","glycopeptide","glycopeptide","peptide","glycopeptide","unclassified","aminoglycoside","bacitracin","fosfomycin","glycopeptide","tetracycline","tetracycline","glycopeptide","multidrug","peptide","aminoglycoside","nucleoside","diaminopyrimidine","MLS","glycopeptide","multidrug","multidrug","MLS","glycopeptide","beta-lactam","diaminopyrimidine","aminoglycoside","peptide","glycopeptide","glycopeptide","multidrug","beta-lactam","multidrug","MLS","MLS","beta-lactam","tetracycline","tetracycline","unclassified","MLS","bacitracin","bacitracin","unclassified","tetracycline","multidrug","glycopeptide","glycopeptide","beta-lactam","MLS","aminoglycoside"],["Plasmid","Plasmid","Plasmid","Plasmid","Plasmid","Plasmid","Plasmid","Plasmid","Plasmid","Plasmid","Uncertain - plasmid or chromosomal","Uncertain - plasmid or chromosomal","Uncertain - plasmid or chromosomal","Uncertain - too short","Uncertain - plasmid or chromosomal","Uncertain - plasmid or chromosomal","Uncertain - too short","Uncertain - plasmid or chromosomal","Uncertain - plasmid or chromosomal","Uncertain - plasmid or chromosomal","Plasmid","Uncertain - plasmid or chromosomal","Uncertain - too short","Uncertain - plasmid or chromosomal","Uncertain - plasmid or chromosomal","Uncertain - too short","Uncertain - plasmid or chromosomal","Uncertain - plasmid or chromosomal","Plasmid","Plasmid","Plasmid","Uncertain - too short","Uncertain - too short","Uncertain - plasmid or chromosomal","Uncertain - too short","Plasmid","Uncertain - too short","Uncertain - too short","Uncertain - plasmid or chromosomal","Uncertain - plasmid or chromosomal","Uncertain - plasmid or chromosomal","Uncertain - plasmid or chromosomal","Uncertain - plasmid or chromosomal","Uncertain - too short","Uncertain - too short","Uncertain - plasmid or chromosomal","Uncertain - plasmid or chromosomal","Uncertain - plasmid or chromosomal","Uncertain - plasmid or chromosomal","Uncertain - too short","Uncertain - plasmid or chromosomal","Uncertain - too short","Uncertain - plasmid or chromosomal","Plasmid","Uncertain - too short","Uncertain - plasmid or chromosomal","Plasmid","Uncertain - too short","Uncertain - plasmid or chromosomal","Uncertain - too short","Plasmid","Uncertain - too short","Plasmid","Uncertain - too short","Uncertain - plasmid or chromosomal","Uncertain - too short","Uncertain - plasmid or chromosomal","Uncertain - too short","Plasmid","Uncertain - too short","Uncertain - too short","Uncertain - plasmid or chromosomal","Uncertain - too short","Uncertain - plasmid or chromosomal","Uncertain - plasmid or chromosomal","Uncertain - plasmid or chromosomal","Plasmid","Plasmid","Uncertain - too short","Uncertain - too short","Uncertain - too short","Uncertain - plasmid or chromosomal","Uncertain - plasmid or chromosomal","Plasmid","Plasmid","Plasmid","Uncertain - too short","Uncertain - too short","Uncertain - plasmid or chromosomal","Uncertain - too short","Uncertain - too short","Uncertain - too short","Uncertain - too short","Uncertain - plasmid or chromosomal","Uncertain - too short","Uncertain - plasmid or chromosomal","Uncertain - too short","Uncertain - too short","Uncertain - plasmid or chromosomal","Uncertain - too short","Uncertain - too short","Uncertain - too short","Uncertain - too short","Uncertain - plasmid or chromosomal","Uncertain - too short","Uncertain - plasmid or chromosomal","Uncertain - too short","Uncertain - plasmid or chromosomal","Uncertain - plasmid or chromosomal","Uncertain - too short","Uncertain - too short","Uncertain - plasmid or chromosomal","Uncertain - plasmid or chromosomal","Uncertain - plasmid or chromosomal","Uncertain - plasmid or chromosomal","Uncertain - too short","Uncertain - plasmid or chromosomal","Uncertain - too short","Uncertain - too short","Uncertain - too short","Uncertain - too short","Uncertain - plasmid or chromosomal","Uncertain - plasmid or chromosomal","Uncertain - too short","Uncertain - plasmid or chromosomal","Uncertain - too short","Uncertain - plasmid or chromosomal","Uncertain - too short","Uncertain - plasmid or chromosomal","Uncertain - plasmid or chromosomal","Uncertain - too short","Uncertain - plasmid or chromosomal","Uncertain - plasmid or chromosomal","Uncertain - too short","Uncertain - plasmid or chromosomal","Uncertain - plasmid or chromosomal","Uncertain - plasmid or chromosomal","Uncertain - plasmid or chromosomal","Uncertain - too short","Plasmid","Plasmid","Plasmid","Uncertain - plasmid or chromosomal","Uncertain - too short","Uncertain - too short","Uncertain - plasmid or chromosomal","Uncertain - too short","Uncertain - plasmid or chromosomal","Uncertain - too short","Uncertain - plasmid or chromosomal","Uncertain - plasmid or chromosomal","Uncertain - plasmid or chromosomal","Uncertain - too short","Uncertain - too short","Uncertain - plasmid or chromosomal","Uncertain - too short","Uncertain - too short","Uncertain - too short","Uncertain - plasmid or chromosomal","Uncertain - plasmid or chromosomal","Uncertain - plasmid or chromosomal","Uncertain - plasmid or chromosomal","Uncertain - plasmid or chromosomal","Uncertain - plasmid or chromosomal","Uncertain - plasmid or chromosomal","Uncertain - plasmid or chromosomal","Uncertain - plasmid or chromosomal","Uncertain - plasmid or chromosomal","Uncertain - too short","Uncertain - too short","Uncertain - too short","Uncertain - too short","Uncertain - too short","Plasmid"],["plasmid","plasmid","plasmid","plasmid","plasmid","plasmid","plasmid","plasmid","plasmid","plasmid","plasmid","plasmid","ambiguous (plasmid/phage)","plasmid","plasmid","ambiguous (plasmid/phage)","plasmid","ambiguous (plasmid/phage)","plasmid","ambiguous (plasmid/phage)","plasmid","plasmid","plasmid","plasmid","plasmid","plasmid","plasmid","plasmid","plasmid","plasmid","plasmid","plasmid","plasmid","plasmid","plasmid","plasmid","plasmid","plasmid","plasmid","plasmid","ambiguous (plasmid/phage)","plasmid","plasmid","plasmid","plasmid","plasmid","ambiguous (plasmid/phage)","plasmid","plasmid","ambiguous (plasmid/phage)","plasmid","plasmid","plasmid","ambiguous (plasmid/phage)","plasmid","plasmid","plasmid","plasmid","plasmid","plasmid","plasmid","plasmid","plasmid","plasmid","plasmid","plasmid","plasmid","plasmid","ambiguous (plasmid/phage)","plasmid","ambiguous (plasmid/phage)","plasmid","plasmid","plasmid","plasmid","ambiguous (plasmid/phage)","plasmid","ambiguous (plasmid/phage)","plasmid","plasmid","plasmid","plasmid","plasmid","plasmid","plasmid","plasmid","plasmid","plasmid","plasmid","plasmid","ambiguous (plasmid/phage)","plasmid","plasmid","plasmid","plasmid","plasmid","plasmid","ambiguous (plasmid/phage)","plasmid","plasmid","plasmid","plasmid","plasmid","plasmid","ambiguous (plasmid/phage)","plasmid","plasmid","plasmid","plasmid","plasmid","plasmid","plasmid","plasmid","plasmid","plasmid","plasmid","plasmid","plasmid","plasmid","plasmid","ambiguous (plasmid/phage)","ambiguous (plasmid/phage)","plasmid","plasmid","ambiguous (plasmid/phage)","plasmid","plasmid","plasmid","plasmid","plasmid","ambiguous (plasmid/phage)","plasmid","plasmid","plasmid","plasmid","plasmid","ambiguous (plasmid/phage)","plasmid","plasmid","plasmid","plasmid","plasmid","plasmid","plasmid","plasmid","plasmid","plasmid","plasmid","ambiguous (plasmid/phage)","plasmid","plasmid","plasmid","ambiguous (plasmid/phage)","plasmid","plasmid","plasmid","plasmid","plasmid","plasmid","plasmid","plasmid","plasmid","plasmid","plasmid","plasmid","plasmid","plasmid","plasmid","plasmid","plasmid","plasmid","plasmid","ambiguous (plasmid/phage)","plasmid"],[null,null,null,null,null,null,null,null,null,null,null,null,null,"maxbin_062_fasta_contigs",null,"maxbin_019_fasta_contigs","maxbin_122_fasta_contigs","maxbin_255_fasta_contigs","maxbin_163_fasta_contigs","maxbin_028_fasta_contigs","maxbin_045_fasta_contigs",null,null,"maxbin_110_fasta_contigs","maxbin_122_fasta_contigs",null,"maxbin_301_fasta_contigs",null,null,null,null,null,"maxbin_069_fasta_sub_contigs","maxbin_317_fasta_contigs",null,null,null,"maxbin_184_fasta_contigs","maxbin_266_fasta_contigs","maxbin_232_fasta_contigs",null,"maxbin_156_fasta_contigs",null,null,"maxbin_006_fasta_contigs","maxbin_047_fasta_contigs","maxbin_023_fasta_contigs",null,null,"maxbin_102_fasta_contigs","maxbin_062_fasta_contigs",null,"maxbin_270_fasta_contigs",null,null,null,null,null,null,"maxbin_017_fasta_contigs","maxbin_211_fasta_contigs",null,"maxbin_166_fasta_contigs","maxbin_016_fasta_contigs","maxbin_062_fasta_contigs",null,"maxbin_051_fasta_sub_contigs",null,"maxbin_015_fasta_contigs",null,null,null,null,null,null,"maxbin_219_fasta_contigs",null,"maxbin_028_fasta_contigs",null,null,null,"maxbin_266_fasta_contigs",null,null,null,null,"maxbin_184_fasta_contigs",null,null,null,null,"maxbin_266_fasta_contigs",null,"maxbin_069_fasta_sub_contigs",null,"maxbin_160_fasta_contigs",null,null,"maxbin_284_fasta_contigs",null,null,null,null,null,null,"maxbin_016_fasta_contigs","maxbin_255_fasta_contigs",null,"maxbin_241_fasta_contigs",null,null,"maxbin_074_fasta_sub_contigs","maxbin_280_fasta_contigs","maxbin_034_fasta_contigs","maxbin_044_fasta_contigs",null,"maxbin_255_fasta_contigs",null,null,null,null,"maxbin_035_fasta_contigs",null,null,"maxbin_112_fasta_contigs",null,null,null,null,null,null,"maxbin_317_fasta_contigs",null,null,null,null,"maxbin_185_fasta_contigs","maxbin_301_fasta_contigs","maxbin_165_fasta_contigs",null,null,null,"maxbin_051_fasta_sub_contigs",null,null,"metabat2_413_fa_sub_contigs",null,null,"maxbin_236_fasta_contigs",null,null,"maxbin_156_fasta_contigs",null,"maxbin_099_fasta_contigs",null,null,"maxbin_315_fasta_contigs",null,"maxbin_163_fasta_contigs",null,"maxbin_028_fasta_contigs",null,null,null,null,null,null,"maxbin_036_fasta_contigs",null,null,null,null,null,"maxbin_239_fasta_contigs"],[null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,"Bin_9_6_1",null,null,null,null,null,"Bin_13_1",null,null,null,null,"Bin_0_13_1","Bin_19",null,null,null,null,null,"Bin_10_3",null,null,null,null,"Bin_19",null,"Bin_19",null,null,null,null,"Bin_13_8","Bin_8_1_1",null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,"Bin_6_25",null,null,null,null,null,null,null,null,null,null,null,null,"Bin_5_8_CTXM",null,null,null,null,"Bin_19","Bin_19","Bin_2_1_1","Bin_2_1_1","Bin_2_1_1",null,null,null,null,null,"Bin_19",null,"Bin_9_7",null,null,null,null,"Bin_14_8",null,null,null,null,null,null,null,"Bin_19",null,null,null,null,null,null,null,"Bin_13_8",null,"Bin_19",null,null,null,null,null,null,null,null,null,null,null,null,null,null,"Bin_10_3","Bin_19",null,null,null,null,"Bin_0_13_1",null,null,null,null,null,null,null,null,null,"Bin_19",null,null,null,null,null,null,null,null,null,null,null,"Bin_19","Bin_7_6_1",null,null,null,null,null,null,"Bin_13_10",null,null,null,null,null,null]],"container":"<table class=\"display\">\n  <thead>\n    <tr>\n      <th> <\/th>\n      <th>ORF_ID<\/th>\n      <th>PathoFact_AMR_ARG<\/th>\n      <th>PathoFact_AMR_AMR_category<\/th>\n      <th>viralverify_Prediction<\/th>\n      <th>PathoFact_AMR_MGE_prediction<\/th>\n      <th>DAS_bin<\/th>\n      <th>CONCOCT_bin<\/th>\n    <\/tr>\n  <\/thead>\n<\/table>","options":{"columnDefs":[{"orderable":false,"targets":0}],"order":[],"autoWidth":false,"orderClasses":false}},"evals":[],"jsHooks":[]}</script>
```

**We have to clarify some cases before saving: e.g.: plasmid but in MAG... it should be plasmid and not MAG -> need to update anvio collection and summarise again to get genomes without plamid contigs... Would like to include plasX as well: https://github.com/michaelkyu/PlasX and build consensus prediction, plasmids are key here.**

**If per sample assembly is retained it is gonna be a huge amount of work. Let's keep all this in mind**



```r
# full_gene_catalogue %>% 
#   mutate(Classification_mobilome_DAS = case_when(
#     PathoFact_AMR_MGE_prediction == "Plasmid" ~ "Plasmid",
#     !is.na(DAS_bin) ~ "MAG", Would like to include
#     TRUE ~ Prediction)) %>%
#   mutate(Classification_mobilome_DAS_HQ = case_when(
#     Prediction == "Plasmid" ~ "Plasmid",
#     !is.na(DAS_bins_HQ) ~ "MAG",
#     TRUE ~ Prediction)) %>%
#   mutate(Classification_mobilome_CONCOCT = case_when(
#     Prediction == "Plasmid" ~ "Plasmid",
#     !is.na(CONCOCT_bin) ~ "MAG",
#     TRUE ~ Prediction)) %>%
#   mutate(Classification_mobilome_CONCOCT_HQ = case_when(
#     Prediction == "Plasmid" ~ "Plasmid",
#     !is.na(CONCOCT_bins_HQ) ~ "MAG",
#     TRUE ~ Prediction)) %>%
#   as.matrix(rownames = TRUE)
```



```r
# 
# 
# left_join(card, by = "ORF_ID") %>%
#   left_join(hgt_events, by = c("Contig" = "CONTIG_NAME")) %>%
#   left_join(viralverify, by = c("Contig" = "Contig name")) %>%
#   left_join(top_isescan, by = c("Contig" = "seqID")) %>%
#   left_join(DAS, by = c("Contig" = "cont_id")) %>%
#   left_join(DAS_tax, by = c("DAS_bin" ="DAS_user_genome")) %>%
#   left_join(DAS_HQ, by = c("DAS_bin" = "DAS_bins_HQ"), keep = TRUE) %>%
#   left_join(CONCOCT, by = c("Contig" = "cont_id")) %>%
#   left_join(CONCOCT_tax, by = c("CONCOCT_bin" ="CONCOCT_user_genome")) %>%
#   left_join(CONCOCT_HQ, by = c("CONCOCT_bin" = "CONCOCT_bins_HQ"), keep = TRUE) %>%
#   select(-c(Megahit, Range)) %>%
#   column_to_rownames("ORF_ID") %>%
#   mutate(Classification_mobilome_DAS = case_when(
#     Prediction == "Plasmid" ~ "Plasmid",
#     !is.na(DAS_bin) ~ "MAG",
#     TRUE ~ Prediction)) %>%
#   mutate(Classification_mobilome_DAS_HQ = case_when(
#     Prediction == "Plasmid" ~ "Plasmid",
#     !is.na(DAS_bins_HQ) ~ "MAG",
#     TRUE ~ Prediction)) %>%
#   mutate(Classification_mobilome_CONCOCT = case_when(
#     Prediction == "Plasmid" ~ "Plasmid",
#     !is.na(CONCOCT_bin) ~ "MAG",
#     TRUE ~ Prediction)) %>%
#   mutate(Classification_mobilome_CONCOCT_HQ = case_when(
#     Prediction == "Plasmid" ~ "Plasmid",
#     !is.na(CONCOCT_bins_HQ) ~ "MAG",
#     TRUE ~ Prediction)) %>%
#   as.matrix(rownames = TRUE) -> human_tax
# 
# 
# rownames(human_tax) <- paste0("h_",rownames(human_tax))
# 
# # Check for hgt events in resistome
# human_tax %>% as.data.frame() %>% filter(!is.na(ARO), CALL == "lgt")
```




##  Extract quantitative information from SQM:

### ORFs:



```r
human$orfs$table %>%
  rownames_to_column("ORF_ID") %>% 
  data.frame() %>% 
  select(ORF_ID, starts_with(c("TPM", "Coverage", "Raw.base.count", "Raw.read.count."))) %>% 
  # select(starts_with("TPM"), "ORF_ID") %>% 
  column_to_rownames("ORF_ID") -> quant_orfs

# colnames(quant) <- colnames(quant) %>% 
#   str_replace_all("//.", "//_")

rownames(quant_orfs) <- paste0("c1_",
                               rownames(quant_orfs))

quant_orfs %>% 
  rownames_to_column("ORF_ID") -> quant_orfs

quant_orfs %>% 
  head() %>% 
  DT::datatable()
```

```{=html}
<div id="htmlwidget-b45f6184cb5d91f8f508" style="width:100%;height:auto;" class="datatables html-widget"></div>
<script type="application/json" data-for="htmlwidget-b45f6184cb5d91f8f508">{"x":{"filter":"none","vertical":false,"data":[["1","2","3","4","5","6"],["c1_megahit_1_2-643","c1_megahit_2_2-232","c1_megahit_3_3-233","c1_megahit_4_2-214","c1_megahit_5_1-270","c1_megahit_6_3-317"],[0.418,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0.636,0,0,0,0,0],[0.341,0,0,0,0,0],[0,0,0,0,0,0],[0.398,0,0,0,0.177,0],[0.131,0,0,0,0,0],[0.287,0,0,0,0,0],[0.746,0,0,0,0,0],[0,0,0,0,0,0],[0.891,0,0,0,0,0],[1.224,0,0,0,0,0],[0.282,0,0,0,0,0],[0,0,0,0,0,0],[0.388,0,0,0,0,0],[0.171,0,0,0,0,0],[0.609,0,0,0,0.121,0],[0,0,0,0,0,0],[1.184,0,0,0,0.251,0],[0.93,0,0,0,0,0],[0.098,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0.644,0,0,0,0,0],[0.047,0,0,0,0,0],[0.533,0,0,0,0,0],[0.291,0,0,0,0,0],[0.503,0,0,0,0,0],[0.718,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[1.347,0,0,0,0,0],[0.26,0,0,0,0,0],[0.348,0,0,0,0,0],[0.158,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0.453,0,0,0,0.045,0],[0.688,0,0,0.059,0,0],[0.438,0,0,0,0,0],[2.444,0,0,0,0.039,0.201],[0.427,0,0,0,0,0],[0.548,0,0,0,0,0],[0.825,0,0,0,0,0],[0.493,0,0,0,0,0],[0.526,0,0,0,0.104,0],[0.3,0,0,0,0,0],[3.863,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[5.037,0,0,0,0,0],[2.804,0,0,0,0,0],[0,0,0,0,0,0],[3.603,0,0,0,1.667,0],[1.402,0,0,0,0,0],[2.539,0,0,0,0,0],[7.579,0,0,0,0,0],[0,0,0,0,0,0],[8.41,0,0,0,0,0],[10.734,0,0,0,0,0],[2.731,0,0,0,0,0],[0,0,0,0,0,0],[4.112,0,0,0,0,0],[1.866,0,0,0,0,0],[5.533,0,0,0,1.107,0],[0,0,0,0,0,0],[13.012,0,0,0,2.778,0],[10.505,0,0,0,0,0],[0.935,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[8.769,0,0,0,0,0],[0.467,0,0,0,0,0],[7.002,0,0,0,0,0],[3.442,0,0,0,0,0],[4.903,0,0,0,0,0],[7.783,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[15.101,0,0,0,0,0],[2.737,0,0,0,0,0],[2.942,0,0,0,0,0],[1.636,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[5.361,0,0,0,0.556,0],[8.081,0,0,0.704,0,0],[4.104,0,0,0,0,0],[34.022,0,0,0,0.556,2.841],[4.008,0,0,0,0,0],[5.994,0,0,0,0,0],[7.603,0,0,0,0,0],[4.693,0,0,0,0,0],[5.268,0,0,0,1.111,0],[2.607,0,0,0,0,0],[2480,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[3234,0,0,0,0,0],[1800,0,0,0,0,0],[0,0,0,0,0,0],[2313,0,0,0,450,0],[900,0,0,0,0,0],[1630,0,0,0,0,0],[4866,0,0,0,0,0],[0,0,0,0,0,0],[5399,0,0,0,0,0],[6891,0,0,0,0,0],[1753,0,0,0,0,0],[0,0,0,0,0,0],[2640,0,0,0,0,0],[1198,0,0,0,0,0],[3552,0,0,0,299,0],[0,0,0,0,0,0],[8354,0,0,0,750,0],[6744,0,0,0,0,0],[600,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[5630,0,0,0,0,0],[300,0,0,0,0,0],[4495,0,0,0,0,0],[2210,0,0,0,0,0],[3148,0,0,0,0,0],[4997,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[9695,0,0,0,0,0],[1757,0,0,0,0,0],[1889,0,0,0,0,0],[1050,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[3442,0,0,0,150,0],[5188,0,0,150,0,0],[2635,0,0,0,0,0],[21842,0,0,0,150,895],[2573,0,0,0,0,0],[3848,0,0,0,0,0],[4881,0,0,0,0,0],[3013,0,0,0,0,0],[3382,0,0,0,300,0],[1674,0,0,0,0,0],[17,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[22,0,0,0,0,0],[12,0,0,0,0,0],[0,0,0,0,0,0],[16,0,0,0,3,0],[6,0,0,0,0,0],[11,0,0,0,0,0],[33,0,0,0,0,0],[0,0,0,0,0,0],[36,0,0,0,0,0],[46,0,0,0,0,0],[12,0,0,0,0,0],[0,0,0,0,0,0],[18,0,0,0,0,0],[8,0,0,0,0,0],[24,0,0,0,2,0],[0,0,0,0,0,0],[56,0,0,0,5,0],[45,0,0,0,0,0],[4,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[38,0,0,0,0,0],[2,0,0,0,0,0],[31,0,0,0,0,0],[15,0,0,0,0,0],[21,0,0,0,0,0],[35,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[66,0,0,0,0,0],[12,0,0,0,0,0],[14,0,0,0,0,0],[7,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[24,0,0,0,1,0],[35,0,0,1,0,0],[18,0,0,0,0,0],[149,0,0,0,1,6],[18,0,0,0,0,0],[26,0,0,0,0,0],[36,0,0,0,0,0],[22,0,0,0,0,0],[24,0,0,0,2,0],[12,0,0,0,0,0]],"container":"<table class=\"display\">\n  <thead>\n    <tr>\n      <th> <\/th>\n      <th>ORF_ID<\/th>\n      <th>TPM.C1<\/th>\n      <th>TPM.C2<\/th>\n      <th>TPM.C3<\/th>\n      <th>TPM.C4<\/th>\n      <th>TPM.C5<\/th>\n      <th>TPM.C6<\/th>\n      <th>TPM.D10<\/th>\n      <th>TPM.D11<\/th>\n      <th>TPM.D12<\/th>\n      <th>TPM.D13<\/th>\n      <th>TPM.D14<\/th>\n      <th>TPM.D15<\/th>\n      <th>TPM.D16<\/th>\n      <th>TPM.D17<\/th>\n      <th>TPM.D18<\/th>\n      <th>TPM.D19<\/th>\n      <th>TPM.D1<\/th>\n      <th>TPM.D20<\/th>\n      <th>TPM.D21<\/th>\n      <th>TPM.D22<\/th>\n      <th>TPM.D23<\/th>\n      <th>TPM.D24<\/th>\n      <th>TPM.D25<\/th>\n      <th>TPM.D26<\/th>\n      <th>TPM.D27<\/th>\n      <th>TPM.D28<\/th>\n      <th>TPM.D29<\/th>\n      <th>TPM.D2<\/th>\n      <th>TPM.D30<\/th>\n      <th>TPM.D31<\/th>\n      <th>TPM.D32<\/th>\n      <th>TPM.D33<\/th>\n      <th>TPM.D34<\/th>\n      <th>TPM.D35<\/th>\n      <th>TPM.D36<\/th>\n      <th>TPM.D37<\/th>\n      <th>TPM.D38<\/th>\n      <th>TPM.D39<\/th>\n      <th>TPM.D3<\/th>\n      <th>TPM.D40<\/th>\n      <th>TPM.D41<\/th>\n      <th>TPM.D42<\/th>\n      <th>TPM.D43<\/th>\n      <th>TPM.D4<\/th>\n      <th>TPM.D5<\/th>\n      <th>TPM.D6<\/th>\n      <th>TPM.D7<\/th>\n      <th>TPM.D8<\/th>\n      <th>TPM.D9<\/th>\n      <th>Coverage.C1<\/th>\n      <th>Coverage.C2<\/th>\n      <th>Coverage.C3<\/th>\n      <th>Coverage.C4<\/th>\n      <th>Coverage.C5<\/th>\n      <th>Coverage.C6<\/th>\n      <th>Coverage.D10<\/th>\n      <th>Coverage.D11<\/th>\n      <th>Coverage.D12<\/th>\n      <th>Coverage.D13<\/th>\n      <th>Coverage.D14<\/th>\n      <th>Coverage.D15<\/th>\n      <th>Coverage.D16<\/th>\n      <th>Coverage.D17<\/th>\n      <th>Coverage.D18<\/th>\n      <th>Coverage.D19<\/th>\n      <th>Coverage.D1<\/th>\n      <th>Coverage.D20<\/th>\n      <th>Coverage.D21<\/th>\n      <th>Coverage.D22<\/th>\n      <th>Coverage.D23<\/th>\n      <th>Coverage.D24<\/th>\n      <th>Coverage.D25<\/th>\n      <th>Coverage.D26<\/th>\n      <th>Coverage.D27<\/th>\n      <th>Coverage.D28<\/th>\n      <th>Coverage.D29<\/th>\n      <th>Coverage.D2<\/th>\n      <th>Coverage.D30<\/th>\n      <th>Coverage.D31<\/th>\n      <th>Coverage.D32<\/th>\n      <th>Coverage.D33<\/th>\n      <th>Coverage.D34<\/th>\n      <th>Coverage.D35<\/th>\n      <th>Coverage.D36<\/th>\n      <th>Coverage.D37<\/th>\n      <th>Coverage.D38<\/th>\n      <th>Coverage.D39<\/th>\n      <th>Coverage.D3<\/th>\n      <th>Coverage.D40<\/th>\n      <th>Coverage.D41<\/th>\n      <th>Coverage.D42<\/th>\n      <th>Coverage.D43<\/th>\n      <th>Coverage.D4<\/th>\n      <th>Coverage.D5<\/th>\n      <th>Coverage.D6<\/th>\n      <th>Coverage.D7<\/th>\n      <th>Coverage.D8<\/th>\n      <th>Coverage.D9<\/th>\n      <th>Raw.base.count.C1<\/th>\n      <th>Raw.base.count.C2<\/th>\n      <th>Raw.base.count.C3<\/th>\n      <th>Raw.base.count.C4<\/th>\n      <th>Raw.base.count.C5<\/th>\n      <th>Raw.base.count.C6<\/th>\n      <th>Raw.base.count.D10<\/th>\n      <th>Raw.base.count.D11<\/th>\n      <th>Raw.base.count.D12<\/th>\n      <th>Raw.base.count.D13<\/th>\n      <th>Raw.base.count.D14<\/th>\n      <th>Raw.base.count.D15<\/th>\n      <th>Raw.base.count.D16<\/th>\n      <th>Raw.base.count.D17<\/th>\n      <th>Raw.base.count.D18<\/th>\n      <th>Raw.base.count.D19<\/th>\n      <th>Raw.base.count.D1<\/th>\n      <th>Raw.base.count.D20<\/th>\n      <th>Raw.base.count.D21<\/th>\n      <th>Raw.base.count.D22<\/th>\n      <th>Raw.base.count.D23<\/th>\n      <th>Raw.base.count.D24<\/th>\n      <th>Raw.base.count.D25<\/th>\n      <th>Raw.base.count.D26<\/th>\n      <th>Raw.base.count.D27<\/th>\n      <th>Raw.base.count.D28<\/th>\n      <th>Raw.base.count.D29<\/th>\n      <th>Raw.base.count.D2<\/th>\n      <th>Raw.base.count.D30<\/th>\n      <th>Raw.base.count.D31<\/th>\n      <th>Raw.base.count.D32<\/th>\n      <th>Raw.base.count.D33<\/th>\n      <th>Raw.base.count.D34<\/th>\n      <th>Raw.base.count.D35<\/th>\n      <th>Raw.base.count.D36<\/th>\n      <th>Raw.base.count.D37<\/th>\n      <th>Raw.base.count.D38<\/th>\n      <th>Raw.base.count.D39<\/th>\n      <th>Raw.base.count.D3<\/th>\n      <th>Raw.base.count.D40<\/th>\n      <th>Raw.base.count.D41<\/th>\n      <th>Raw.base.count.D42<\/th>\n      <th>Raw.base.count.D43<\/th>\n      <th>Raw.base.count.D4<\/th>\n      <th>Raw.base.count.D5<\/th>\n      <th>Raw.base.count.D6<\/th>\n      <th>Raw.base.count.D7<\/th>\n      <th>Raw.base.count.D8<\/th>\n      <th>Raw.base.count.D9<\/th>\n      <th>Raw.read.count.C1<\/th>\n      <th>Raw.read.count.C2<\/th>\n      <th>Raw.read.count.C3<\/th>\n      <th>Raw.read.count.C4<\/th>\n      <th>Raw.read.count.C5<\/th>\n      <th>Raw.read.count.C6<\/th>\n      <th>Raw.read.count.D10<\/th>\n      <th>Raw.read.count.D11<\/th>\n      <th>Raw.read.count.D12<\/th>\n      <th>Raw.read.count.D13<\/th>\n      <th>Raw.read.count.D14<\/th>\n      <th>Raw.read.count.D15<\/th>\n      <th>Raw.read.count.D16<\/th>\n      <th>Raw.read.count.D17<\/th>\n      <th>Raw.read.count.D18<\/th>\n      <th>Raw.read.count.D19<\/th>\n      <th>Raw.read.count.D1<\/th>\n      <th>Raw.read.count.D20<\/th>\n      <th>Raw.read.count.D21<\/th>\n      <th>Raw.read.count.D22<\/th>\n      <th>Raw.read.count.D23<\/th>\n      <th>Raw.read.count.D24<\/th>\n      <th>Raw.read.count.D25<\/th>\n      <th>Raw.read.count.D26<\/th>\n      <th>Raw.read.count.D27<\/th>\n      <th>Raw.read.count.D28<\/th>\n      <th>Raw.read.count.D29<\/th>\n      <th>Raw.read.count.D2<\/th>\n      <th>Raw.read.count.D30<\/th>\n      <th>Raw.read.count.D31<\/th>\n      <th>Raw.read.count.D32<\/th>\n      <th>Raw.read.count.D33<\/th>\n      <th>Raw.read.count.D34<\/th>\n      <th>Raw.read.count.D35<\/th>\n      <th>Raw.read.count.D36<\/th>\n      <th>Raw.read.count.D37<\/th>\n      <th>Raw.read.count.D38<\/th>\n      <th>Raw.read.count.D39<\/th>\n      <th>Raw.read.count.D3<\/th>\n      <th>Raw.read.count.D40<\/th>\n      <th>Raw.read.count.D41<\/th>\n      <th>Raw.read.count.D42<\/th>\n      <th>Raw.read.count.D43<\/th>\n      <th>Raw.read.count.D4<\/th>\n      <th>Raw.read.count.D5<\/th>\n      <th>Raw.read.count.D6<\/th>\n      <th>Raw.read.count.D7<\/th>\n      <th>Raw.read.count.D8<\/th>\n      <th>Raw.read.count.D9<\/th>\n    <\/tr>\n  <\/thead>\n<\/table>","options":{"columnDefs":[{"className":"dt-right","targets":[2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58,59,60,61,62,63,64,65,66,67,68,69,70,71,72,73,74,75,76,77,78,79,80,81,82,83,84,85,86,87,88,89,90,91,92,93,94,95,96,97,98,99,100,101,102,103,104,105,106,107,108,109,110,111,112,113,114,115,116,117,118,119,120,121,122,123,124,125,126,127,128,129,130,131,132,133,134,135,136,137,138,139,140,141,142,143,144,145,146,147,148,149,150,151,152,153,154,155,156,157,158,159,160,161,162,163,164,165,166,167,168,169,170,171,172,173,174,175,176,177,178,179,180,181,182,183,184,185,186,187,188,189,190,191,192,193,194,195,196,197]},{"orderable":false,"targets":0}],"order":[],"autoWidth":false,"orderClasses":false}},"evals":[],"jsHooks":[]}</script>
```

```r
# colnames(human_otu_tpm) <- colnames(human_otu_tpm) %>% str_remove_all("TPM ") %>% str_replace_all("-", "_")
# rownames(human_otu_tpm) <- paste0("h_", rownames(human_otu_tpm))
# 
# # Coverage
# human$orfs$table %>% as.data.frame() %>% select(starts_with("Coverage")) -> human_otu_cov
# 
# colnames(human_otu_cov) <- colnames(human_otu_cov) %>% str_remove_all("Coverage ") %>% str_replace_all("-", "_")
# rownames(human_otu_cov) <- paste0("h_", rownames(human_otu_cov))
# 
# # Raw reads
# human$orfs$table %>% as.data.frame() %>% select(starts_with("Raw read count")) -> human_otu_rr
# 
# colnames(human_otu_rr) <- colnames(human_otu_rr) %>% str_remove_all("Raw read count ") %>% str_replace_all("-", "_")
# rownames(human_otu_rr) <- paste0("h_", rownames(human_otu_rr))
```

### contigs:


```r
human$contigs$table %>%
  rownames_to_column("Contig_ID") %>% 
  data.frame() %>% 
    select(Contig_ID, starts_with(c("TPM", "Coverage", "Raw.base.count", "Raw.read.count."))) %>% 
  # select(Contig_ID, Coverage.C1:Raw.read.count.D9) %>% 
  # select(starts_with("TPM"), "ORF_ID") %>% 
  column_to_rownames("Contig_ID") -> quant_contigs

# colnames(quant_contigs) <- colnames(quant_contigs) %>% 
#   str_replace_all("//.", "//_")

rownames(quant_contigs) <- paste0("c1_",
                                  rownames(quant_contigs))

quant_contigs %>% 
  rownames_to_column("Contig_ID") -> quant_contigs

# colnames(human_otu_tpm) <- colnames(human_otu_tpm) %>% str_remove_all("TPM ") %>% str_replace_all("-", "_")
# rownames(human_otu_tpm) <- paste0("h_", rownames(human_otu_tpm))
# 
# # Coverage
# human$orfs$table %>% as.data.frame() %>% select(starts_with("Coverage")) -> human_otu_cov
# 
# colnames(human_otu_cov) <- colnames(human_otu_cov) %>% str_remove_all("Coverage ") %>% str_replace_all("-", "_")
# rownames(human_otu_cov) <- paste0("h_", rownames(human_otu_cov))
# 
# # Raw reads
# human$orfs$table %>% as.data.frame() %>% select(starts_with("Raw read count")) -> human_otu_rr
# 
# colnames(human_otu_rr) <- colnames(human_otu_rr) %>% str_remove_all("Raw read count ") %>% str_replace_all("-", "_")
# rownames(human_otu_rr) <- paste0("h_", rownames(human_otu_rr))


quant_contigs %>% 
  head() %>% 
  DT::datatable()
```

```{=html}
<div id="htmlwidget-654fd51560b287cbb9c7" style="width:100%;height:auto;" class="datatables html-widget"></div>
<script type="application/json" data-for="htmlwidget-654fd51560b287cbb9c7">{"x":{"filter":"none","vertical":false,"data":[["1","2","3","4","5","6"],["c1_megahit_1","c1_megahit_2","c1_megahit_3","c1_megahit_4","c1_megahit_5","c1_megahit_6"],[3.1,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[5,0,0,0,0,0],[2.6,0,0,0,0,0],[0,0,0,0,0,0],[3,0,0,0,1.5,0],[1.5,0,0,0,0,0],[2.9,0,0,0,0,0],[6,0,0,0,0,0],[0,0,0,0,0,0],[7.2,0,0,0,0,0],[12,0,0,0,0,0],[2.3,0,0,0,0,0],[0,0,0,0,0,0],[3,0,0,0,0,0],[1.8,0,0,0,0,0],[5.3,0,0,0,1.7,0],[0,0,0,0,0,0],[9.2,0,0,0,2.1,0],[7.3,0,0,0,0,0],[0.9,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[6.6,0,0,0,0,0],[0.4,0,0,0,0,0],[4.1,0,0,0,0,0],[2.4,0,0,0,0,0],[3.9,0,0,0,0,0],[5.7,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[9.9,0,0,0,0,0],[2.4,0,0,0,0,0],[3.1,0,0,0,0,0],[1.4,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[3.6,0,0,0,0.8,0],[5.6,0,0,0.5,0,0],[3.3,0,0,0,0,0],[10.7,0,0,0,0.2,1.1],[3.9,0,0,0,0,0],[4.4,0,0,0,0,0],[7.4,0,0,0,0,0],[4.3,0,0,0,0,0],[3.7,0,0,0,0.8,0],[2.5,0,0,0,0,0],[3.62,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[4.69,0,0,0,0,0],[2.56,0,0,0,0,0],[0,0,0,0,0,0],[3.41,0,0,0,1.67,0],[1.28,0,0,0,0,0],[2.35,0,0,0,0,0],[7.04,0,0,0,0,0],[0,0,0,0,0,0],[7.68,0,0,0,0,0],[9.82,0,0,0,0,0],[2.56,0,0,0,0,0],[0,0,0,0,0,0],[3.84,0,0,0,0,0],[1.71,0,0,0,0,0],[5.12,0,0,0,1.66,0],[0,0,0,0,0,0],[11.95,0,0,0,2.78,0],[9.6,0,0,0,0,0],[0.85,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[8.11,0,0,0,0,0],[0.43,0,0,0,0,0],[6.61,0,0,0,0,0],[3.2,0,0,0,0,0],[4.48,0,0,0,0,0],[7.47,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[14.08,0,0,0,0,0],[2.56,0,0,0,0,0],[2.99,0,0,0,0,0],[1.49,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[5.12,0,0,0,1.11,0],[7.47,0,0,0.7,0,0],[3.82,0,0,0,0,0],[31.78,0,0,0,0.56,3.3],[3.83,0,0,0,0,0],[5.55,0,0,0,0,0],[7.68,0,0,0,0,0],[4.69,0,0,0,0,0],[5.12,0,0,0,1.11,0],[2.56,0,0,0,0,0],[17,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[22,0,0,0,0,0],[12,0,0,0,0,0],[0,0,0,0,0,0],[16,0,0,0,3,0],[6,0,0,0,0,0],[11,0,0,0,0,0],[33,0,0,0,0,0],[0,0,0,0,0,0],[36,0,0,0,0,0],[46,0,0,0,0,0],[12,0,0,0,0,0],[0,0,0,0,0,0],[18,0,0,0,0,0],[8,0,0,0,0,0],[24,0,0,0,3,0],[0,0,0,0,0,0],[56,0,0,0,5,0],[45,0,0,0,0,0],[4,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[38,0,0,0,0,0],[2,0,0,0,0,0],[31,0,0,0,0,0],[15,0,0,0,0,0],[21,0,0,0,0,0],[35,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[66,0,0,0,0,0],[12,0,0,0,0,0],[14,0,0,0,0,0],[7,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[24,0,0,0,2,0],[35,0,0,1,0,0],[18,0,0,0,0,0],[149,0,0,0,1,7],[18,0,0,0,0,0],[26,0,0,0,0,0],[36,0,0,0,0,0],[22,0,0,0,0,0],[24,0,0,0,2,0],[12,0,0,0,0,0]],"container":"<table class=\"display\">\n  <thead>\n    <tr>\n      <th> <\/th>\n      <th>Contig_ID<\/th>\n      <th>TPM.C1<\/th>\n      <th>TPM.C2<\/th>\n      <th>TPM.C3<\/th>\n      <th>TPM.C4<\/th>\n      <th>TPM.C5<\/th>\n      <th>TPM.C6<\/th>\n      <th>TPM.D10<\/th>\n      <th>TPM.D11<\/th>\n      <th>TPM.D12<\/th>\n      <th>TPM.D13<\/th>\n      <th>TPM.D14<\/th>\n      <th>TPM.D15<\/th>\n      <th>TPM.D16<\/th>\n      <th>TPM.D17<\/th>\n      <th>TPM.D18<\/th>\n      <th>TPM.D19<\/th>\n      <th>TPM.D1<\/th>\n      <th>TPM.D20<\/th>\n      <th>TPM.D21<\/th>\n      <th>TPM.D22<\/th>\n      <th>TPM.D23<\/th>\n      <th>TPM.D24<\/th>\n      <th>TPM.D25<\/th>\n      <th>TPM.D26<\/th>\n      <th>TPM.D27<\/th>\n      <th>TPM.D28<\/th>\n      <th>TPM.D29<\/th>\n      <th>TPM.D2<\/th>\n      <th>TPM.D30<\/th>\n      <th>TPM.D31<\/th>\n      <th>TPM.D32<\/th>\n      <th>TPM.D33<\/th>\n      <th>TPM.D34<\/th>\n      <th>TPM.D35<\/th>\n      <th>TPM.D36<\/th>\n      <th>TPM.D37<\/th>\n      <th>TPM.D38<\/th>\n      <th>TPM.D39<\/th>\n      <th>TPM.D3<\/th>\n      <th>TPM.D40<\/th>\n      <th>TPM.D41<\/th>\n      <th>TPM.D42<\/th>\n      <th>TPM.D43<\/th>\n      <th>TPM.D4<\/th>\n      <th>TPM.D5<\/th>\n      <th>TPM.D6<\/th>\n      <th>TPM.D7<\/th>\n      <th>TPM.D8<\/th>\n      <th>TPM.D9<\/th>\n      <th>Coverage.C1<\/th>\n      <th>Coverage.C2<\/th>\n      <th>Coverage.C3<\/th>\n      <th>Coverage.C4<\/th>\n      <th>Coverage.C5<\/th>\n      <th>Coverage.C6<\/th>\n      <th>Coverage.D10<\/th>\n      <th>Coverage.D11<\/th>\n      <th>Coverage.D12<\/th>\n      <th>Coverage.D13<\/th>\n      <th>Coverage.D14<\/th>\n      <th>Coverage.D15<\/th>\n      <th>Coverage.D16<\/th>\n      <th>Coverage.D17<\/th>\n      <th>Coverage.D18<\/th>\n      <th>Coverage.D19<\/th>\n      <th>Coverage.D1<\/th>\n      <th>Coverage.D20<\/th>\n      <th>Coverage.D21<\/th>\n      <th>Coverage.D22<\/th>\n      <th>Coverage.D23<\/th>\n      <th>Coverage.D24<\/th>\n      <th>Coverage.D25<\/th>\n      <th>Coverage.D26<\/th>\n      <th>Coverage.D27<\/th>\n      <th>Coverage.D28<\/th>\n      <th>Coverage.D29<\/th>\n      <th>Coverage.D2<\/th>\n      <th>Coverage.D30<\/th>\n      <th>Coverage.D31<\/th>\n      <th>Coverage.D32<\/th>\n      <th>Coverage.D33<\/th>\n      <th>Coverage.D34<\/th>\n      <th>Coverage.D35<\/th>\n      <th>Coverage.D36<\/th>\n      <th>Coverage.D37<\/th>\n      <th>Coverage.D38<\/th>\n      <th>Coverage.D39<\/th>\n      <th>Coverage.D3<\/th>\n      <th>Coverage.D40<\/th>\n      <th>Coverage.D41<\/th>\n      <th>Coverage.D42<\/th>\n      <th>Coverage.D43<\/th>\n      <th>Coverage.D4<\/th>\n      <th>Coverage.D5<\/th>\n      <th>Coverage.D6<\/th>\n      <th>Coverage.D7<\/th>\n      <th>Coverage.D8<\/th>\n      <th>Coverage.D9<\/th>\n      <th>Raw.read.count.C1<\/th>\n      <th>Raw.read.count.C2<\/th>\n      <th>Raw.read.count.C3<\/th>\n      <th>Raw.read.count.C4<\/th>\n      <th>Raw.read.count.C5<\/th>\n      <th>Raw.read.count.C6<\/th>\n      <th>Raw.read.count.D10<\/th>\n      <th>Raw.read.count.D11<\/th>\n      <th>Raw.read.count.D12<\/th>\n      <th>Raw.read.count.D13<\/th>\n      <th>Raw.read.count.D14<\/th>\n      <th>Raw.read.count.D15<\/th>\n      <th>Raw.read.count.D16<\/th>\n      <th>Raw.read.count.D17<\/th>\n      <th>Raw.read.count.D18<\/th>\n      <th>Raw.read.count.D19<\/th>\n      <th>Raw.read.count.D1<\/th>\n      <th>Raw.read.count.D20<\/th>\n      <th>Raw.read.count.D21<\/th>\n      <th>Raw.read.count.D22<\/th>\n      <th>Raw.read.count.D23<\/th>\n      <th>Raw.read.count.D24<\/th>\n      <th>Raw.read.count.D25<\/th>\n      <th>Raw.read.count.D26<\/th>\n      <th>Raw.read.count.D27<\/th>\n      <th>Raw.read.count.D28<\/th>\n      <th>Raw.read.count.D29<\/th>\n      <th>Raw.read.count.D2<\/th>\n      <th>Raw.read.count.D30<\/th>\n      <th>Raw.read.count.D31<\/th>\n      <th>Raw.read.count.D32<\/th>\n      <th>Raw.read.count.D33<\/th>\n      <th>Raw.read.count.D34<\/th>\n      <th>Raw.read.count.D35<\/th>\n      <th>Raw.read.count.D36<\/th>\n      <th>Raw.read.count.D37<\/th>\n      <th>Raw.read.count.D38<\/th>\n      <th>Raw.read.count.D39<\/th>\n      <th>Raw.read.count.D3<\/th>\n      <th>Raw.read.count.D40<\/th>\n      <th>Raw.read.count.D41<\/th>\n      <th>Raw.read.count.D42<\/th>\n      <th>Raw.read.count.D43<\/th>\n      <th>Raw.read.count.D4<\/th>\n      <th>Raw.read.count.D5<\/th>\n      <th>Raw.read.count.D6<\/th>\n      <th>Raw.read.count.D7<\/th>\n      <th>Raw.read.count.D8<\/th>\n      <th>Raw.read.count.D9<\/th>\n    <\/tr>\n  <\/thead>\n<\/table>","options":{"columnDefs":[{"className":"dt-right","targets":[2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58,59,60,61,62,63,64,65,66,67,68,69,70,71,72,73,74,75,76,77,78,79,80,81,82,83,84,85,86,87,88,89,90,91,92,93,94,95,96,97,98,99,100,101,102,103,104,105,106,107,108,109,110,111,112,113,114,115,116,117,118,119,120,121,122,123,124,125,126,127,128,129,130,131,132,133,134,135,136,137,138,139,140,141,142,143,144,145,146,147,148]},{"orderable":false,"targets":0}],"order":[],"autoWidth":false,"orderClasses":false}},"evals":[],"jsHooks":[]}</script>
```

Combine quantification data with annotations:


```r
colnames(quant_orfs) <- paste0("ORF_" ,colnames(quant_orfs))

colnames(quant_contigs) <- paste0("contig_" ,colnames(quant_contigs))
```



```r
full_gene_catalogue %>% 
    distinct(ORF_ID, .keep_all = TRUE) %>% # see chucnk below
  left_join(quant_orfs,
            by = c("ORF_ID" = "ORF_ORF_ID")) -> full_gene_catalogue_orf_quant

full_gene_catalogue_orf_quant %>% 
    distinct(ORF_ID, .keep_all = TRUE) %>% # see chucnk below

  left_join(quant_contigs,
            by = c("Contig_ID" = "contig_Contig_ID")) -> full_quant_full_gene_catalogue
```



```r
rm(list=setdiff(ls(), c("full_gene_catalogue_orf_quant", "human_meta")))
```



```r
gc()
```

```
##             used   (Mb) gc trigger    (Mb) limit (Mb)   max used    (Mb)
## Ncells   7873527  420.5   18356477   980.4         NA   18356477   980.4
## Vcells 426805988 3256.3 3043161888 23217.5    1024000 3445425436 26286.6
```



```r
full_gene_catalogue_orf_quant %>% 
  distinct(ORF_ID, .keep_all = TRUE) %>% # see chucnk below
  dplyr::select(ORF_ID, SQM_Length.NT, starts_with("ORF_Raw.read.count.")) %>% 
  pivot_longer(cols = starts_with("ORF_Raw.read.count."),
               values_to = "read_counts", 
               names_to = "sample")  -> test

# dim(test)
```




```r
# full_gene_catalogue_orf_quant %>% 
#   dplyr::select(ORF_ID, SQM_Length.NT, starts_with("ORF_Raw.read.count.")) %>% group_by(ORF_ID, SQM_Length.NT) %>% dplyr::summarise(n = dplyr::n(), .groups = "drop") %>%  dplyr::filter(n > 1L)   
#   ORF_ID                   SQM_Length.NT     n
#   <chr>                            <int> <int>
# 1 c1_megahit_288931_1-1971          1971     2
```


```r
# test %>% 
#   head() 
```


```r
# test %>% 
#   group_by(sample, SQM_Length.NT, ORF_ID) %>% 
#   dplyr::summarise(n = dplyr::n(), .groups = "drop") %>% 
#   dplyr::filter(n > 1L) -> test_test 
```



```r
gc()
```

```
##             used   (Mb) gc trigger    (Mb) limit (Mb)   max used    (Mb)
## Ncells   7876925  420.7   18356477   980.4         NA   18356477   980.4
## Vcells 667155921 5090.0 2434529511 18574.0    1024000 3445425436 26286.6
```



```r
test %>% # -> test  test %>% %>%  arrange(count)
  ## group_by(sample) %>% 
  mutate(Num_Gi = read_counts / SQM_Length.NT) %>% 
  ## ungroup() %>% 
  select(-read_counts) %>%
  #  dplyr::group_by(ORF_ID, SQM_Length.NT, sample) %>%
  # dplyr::summarise(n = dplyr::n(), .groups = "drop") %>%
  # dplyr::filter(n > 1L)
  pivot_wider(names_from =  sample,
              values_from = Num_Gi) %>% 
  replace(is.na(.), 0) %>% 
  rename_at(
    # select all variables with "factor.sample" in the name
    vars(contains("ORF_Raw.read.count."))
    # use stringr::str_replace to remove factor.sample.
    #   you could do the same with base::gsub()
    , funs(str_replace(., "ORF_Raw.read.count.", "ORF_Num_Gi."))) %>% 
  select(-SQM_Length.NT) -> full_gene_catalogue_orf_Num_Gi
```

```
## Warning: `funs()` was deprecated in dplyr 0.8.0.
## Please use a list of either functions or lambdas: 
## 
##   # Simple named list: 
##   list(mean = mean, median = median)
## 
##   # Auto named with `tibble::lst()`: 
##   tibble::lst(mean, median)
## 
##   # Using lambdas
##   list(~ mean(., trim = .2), ~ median(., na.rm = TRUE))
## This warning is displayed once every 8 hours.
## Call `lifecycle::last_lifecycle_warnings()` to see where this warning was generated.
```

```r
full_gene_catalogue_orf_Num_Gi %>%  
  head() %>% 
  DT::datatable()
```

```{=html}
<div id="htmlwidget-90e76ebd62df295d5826" style="width:100%;height:auto;" class="datatables html-widget"></div>
<script type="application/json" data-for="htmlwidget-90e76ebd62df295d5826">{"x":{"filter":"none","vertical":false,"data":[["1","2","3","4","5","6"],["c1_megahit_1_2-643","c1_megahit_2_2-232","c1_megahit_3_3-233","c1_megahit_4_2-214","c1_megahit_5_1-270","c1_megahit_6_3-317"],[0.0264797507788162,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0.0342679127725857,0,0,0,0,0],[0.0186915887850467,0,0,0,0,0],[0,0,0,0,0,0],[0.0249221183800623,0,0,0,0.0111111111111111,0],[0.00934579439252336,0,0,0,0,0],[0.0171339563862928,0,0,0,0,0],[0.0514018691588785,0,0,0,0,0],[0,0,0,0,0,0],[0.0560747663551402,0,0,0,0,0],[0.0716510903426791,0,0,0,0,0],[0.0186915887850467,0,0,0,0,0],[0,0,0,0,0,0],[0.0280373831775701,0,0,0,0,0],[0.0124610591900312,0,0,0,0,0],[0.0373831775700935,0,0,0,0.00740740740740741,0],[0,0,0,0,0,0],[0.0872274143302181,0,0,0,0.0185185185185185,0],[0.0700934579439252,0,0,0,0,0],[0.00623052959501558,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0.059190031152648,0,0,0,0,0],[0.00311526479750779,0,0,0,0,0],[0.0482866043613707,0,0,0,0,0],[0.0233644859813084,0,0,0,0,0],[0.0327102803738318,0,0,0,0,0],[0.0545171339563863,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0.102803738317757,0,0,0,0,0],[0.0186915887850467,0,0,0,0,0],[0.0218068535825545,0,0,0,0,0],[0.0109034267912773,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0.0373831775700935,0,0,0,0.0037037037037037,0],[0.0545171339563863,0,0,0.00469483568075117,0,0],[0.0280373831775701,0,0,0,0,0],[0.23208722741433,0,0,0,0.0037037037037037,0.019047619047619],[0.0280373831775701,0,0,0,0,0],[0.0404984423676012,0,0,0,0,0],[0.0560747663551402,0,0,0,0,0],[0.0342679127725857,0,0,0,0,0],[0.0373831775700935,0,0,0,0.00740740740740741,0],[0.0186915887850467,0,0,0,0,0]],"container":"<table class=\"display\">\n  <thead>\n    <tr>\n      <th> <\/th>\n      <th>ORF_ID<\/th>\n      <th>ORF_Num_Gi.C1<\/th>\n      <th>ORF_Num_Gi.C2<\/th>\n      <th>ORF_Num_Gi.C3<\/th>\n      <th>ORF_Num_Gi.C4<\/th>\n      <th>ORF_Num_Gi.C5<\/th>\n      <th>ORF_Num_Gi.C6<\/th>\n      <th>ORF_Num_Gi.D10<\/th>\n      <th>ORF_Num_Gi.D11<\/th>\n      <th>ORF_Num_Gi.D12<\/th>\n      <th>ORF_Num_Gi.D13<\/th>\n      <th>ORF_Num_Gi.D14<\/th>\n      <th>ORF_Num_Gi.D15<\/th>\n      <th>ORF_Num_Gi.D16<\/th>\n      <th>ORF_Num_Gi.D17<\/th>\n      <th>ORF_Num_Gi.D18<\/th>\n      <th>ORF_Num_Gi.D19<\/th>\n      <th>ORF_Num_Gi.D1<\/th>\n      <th>ORF_Num_Gi.D20<\/th>\n      <th>ORF_Num_Gi.D21<\/th>\n      <th>ORF_Num_Gi.D22<\/th>\n      <th>ORF_Num_Gi.D23<\/th>\n      <th>ORF_Num_Gi.D24<\/th>\n      <th>ORF_Num_Gi.D25<\/th>\n      <th>ORF_Num_Gi.D26<\/th>\n      <th>ORF_Num_Gi.D27<\/th>\n      <th>ORF_Num_Gi.D28<\/th>\n      <th>ORF_Num_Gi.D29<\/th>\n      <th>ORF_Num_Gi.D2<\/th>\n      <th>ORF_Num_Gi.D30<\/th>\n      <th>ORF_Num_Gi.D31<\/th>\n      <th>ORF_Num_Gi.D32<\/th>\n      <th>ORF_Num_Gi.D33<\/th>\n      <th>ORF_Num_Gi.D34<\/th>\n      <th>ORF_Num_Gi.D35<\/th>\n      <th>ORF_Num_Gi.D36<\/th>\n      <th>ORF_Num_Gi.D37<\/th>\n      <th>ORF_Num_Gi.D38<\/th>\n      <th>ORF_Num_Gi.D39<\/th>\n      <th>ORF_Num_Gi.D3<\/th>\n      <th>ORF_Num_Gi.D40<\/th>\n      <th>ORF_Num_Gi.D41<\/th>\n      <th>ORF_Num_Gi.D42<\/th>\n      <th>ORF_Num_Gi.D43<\/th>\n      <th>ORF_Num_Gi.D4<\/th>\n      <th>ORF_Num_Gi.D5<\/th>\n      <th>ORF_Num_Gi.D6<\/th>\n      <th>ORF_Num_Gi.D7<\/th>\n      <th>ORF_Num_Gi.D8<\/th>\n      <th>ORF_Num_Gi.D9<\/th>\n    <\/tr>\n  <\/thead>\n<\/table>","options":{"columnDefs":[{"className":"dt-right","targets":[2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50]},{"orderable":false,"targets":0}],"order":[],"autoWidth":false,"orderClasses":false}},"evals":[],"jsHooks":[]}</script>
```


```r
full_gene_catalogue_orf_Num_Gi %>% 
  pivot_longer(cols = starts_with("ORF_Num_Gi"),
               values_to = "Num_Gi", 
               names_to = "sample") %>% 
  group_by(sample) %>% 
  mutate(Num_Gi_pc = Num_Gi / sum(Num_Gi)) %>% 
  ungroup() %>% 
  select(-Num_Gi) %>% 
  pivot_wider(names_from =  sample,
              values_from = Num_Gi_pc) %>% 
  replace(is.na(.), 0) %>% 
  rename_at(
    # select all variables with "factor.sample" in the name
    vars(contains("ORF_Num_Gi."))
    # use stringr::str_replace to remove factor.sample.
    #   you could do the same with base::gsub()
    , funs(str_replace(., "ORF_Num_Gi.", "ORF_Num_Gi_pc."))) -> full_gene_catalogue_orf_Num_Gi_pc

full_gene_catalogue_orf_Num_Gi_pc %>%  
  head() %>% 
  DT::datatable()
```

```{=html}
<div id="htmlwidget-945a9f22442a68f7a657" style="width:100%;height:auto;" class="datatables html-widget"></div>
<script type="application/json" data-for="htmlwidget-945a9f22442a68f7a657">{"x":{"filter":"none","vertical":false,"data":[["1","2","3","4","5","6"],["c1_megahit_1_2-643","c1_megahit_2_2-232","c1_megahit_3_3-233","c1_megahit_4_2-214","c1_megahit_5_1-270","c1_megahit_6_3-317"],[4.18218387552105e-07,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[6.35991340521448e-07,0,0,0,0,0],[3.40789472663808e-07,0,0,0,0,0],[0,0,0,0,0,0],[3.9791607268332e-07,0,0,0,1.77404249071313e-07,0],[1.30976085270764e-07,0,0,0,0,0],[2.87267978893692e-07,0,0,0,0,0],[7.46125885846709e-07,0,0,0,0,0],[0,0,0,0,0,0],[8.91265562494407e-07,0,0,0,0,0],[1.22448009884397e-06,0,0,0,0,0],[2.81850780207514e-07,0,0,0,0,0],[0,0,0,0,0,0],[3.88095694274393e-07,0,0,0,0,0],[1.70511587017238e-07,0,0,0,0,0],[6.09281981569173e-07,0,0,0,1.20728096347966e-07,0],[0,0,0,0,0,0],[1.1839215853555e-06,0,0,0,2.51348431811584e-07,0],[9.30275767313403e-07,0,0,0,0,0],[9.78506390165811e-08,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[6.43788672958596e-07,0,0,0,0,0],[4.68366776038295e-08,0,0,0,0,0],[5.32701057606594e-07,0,0,0,0,0],[2.91270220674745e-07,0,0,0,0,0],[5.03379670304167e-07,0,0,0,0,0],[7.17987661565603e-07,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[1.34722104603446e-06,0,0,0,0,0],[2.59902315936775e-07,0,0,0,0,0],[3.48159907291395e-07,0,0,0,0,0],[1.58145956604885e-07,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[4.52874987489062e-07,0,0,0,4.48681700567867e-08,0],[6.8844762320545e-07,0,0,5.92868375718175e-08,0,0],[4.37820837600919e-07,0,0,0,0,0],[2.44427886107983e-06,0,0,0,3.90063889836751e-08,2.00604286201758e-07],[4.2723042733103e-07,0,0,0,0,0],[5.48031839988033e-07,0,0,0,0,0],[8.24595434632126e-07,0,0,0,0,0],[4.92855349788805e-07,0,0,0,0,0],[5.26258261511937e-07,0,0,0,1.04277099966254e-07,0],[3.00454253838188e-07,0,0,0,0,0]],"container":"<table class=\"display\">\n  <thead>\n    <tr>\n      <th> <\/th>\n      <th>ORF_ID<\/th>\n      <th>ORF_Num_Gi_pc.C1<\/th>\n      <th>ORF_Num_Gi_pc.C2<\/th>\n      <th>ORF_Num_Gi_pc.C3<\/th>\n      <th>ORF_Num_Gi_pc.C4<\/th>\n      <th>ORF_Num_Gi_pc.C5<\/th>\n      <th>ORF_Num_Gi_pc.C6<\/th>\n      <th>ORF_Num_Gi_pc.D10<\/th>\n      <th>ORF_Num_Gi_pc.D11<\/th>\n      <th>ORF_Num_Gi_pc.D12<\/th>\n      <th>ORF_Num_Gi_pc.D13<\/th>\n      <th>ORF_Num_Gi_pc.D14<\/th>\n      <th>ORF_Num_Gi_pc.D15<\/th>\n      <th>ORF_Num_Gi_pc.D16<\/th>\n      <th>ORF_Num_Gi_pc.D17<\/th>\n      <th>ORF_Num_Gi_pc.D18<\/th>\n      <th>ORF_Num_Gi_pc.D19<\/th>\n      <th>ORF_Num_Gi_pc.D1<\/th>\n      <th>ORF_Num_Gi_pc.D20<\/th>\n      <th>ORF_Num_Gi_pc.D21<\/th>\n      <th>ORF_Num_Gi_pc.D22<\/th>\n      <th>ORF_Num_Gi_pc.D23<\/th>\n      <th>ORF_Num_Gi_pc.D24<\/th>\n      <th>ORF_Num_Gi_pc.D25<\/th>\n      <th>ORF_Num_Gi_pc.D26<\/th>\n      <th>ORF_Num_Gi_pc.D27<\/th>\n      <th>ORF_Num_Gi_pc.D28<\/th>\n      <th>ORF_Num_Gi_pc.D29<\/th>\n      <th>ORF_Num_Gi_pc.D2<\/th>\n      <th>ORF_Num_Gi_pc.D30<\/th>\n      <th>ORF_Num_Gi_pc.D31<\/th>\n      <th>ORF_Num_Gi_pc.D32<\/th>\n      <th>ORF_Num_Gi_pc.D33<\/th>\n      <th>ORF_Num_Gi_pc.D34<\/th>\n      <th>ORF_Num_Gi_pc.D35<\/th>\n      <th>ORF_Num_Gi_pc.D36<\/th>\n      <th>ORF_Num_Gi_pc.D37<\/th>\n      <th>ORF_Num_Gi_pc.D38<\/th>\n      <th>ORF_Num_Gi_pc.D39<\/th>\n      <th>ORF_Num_Gi_pc.D3<\/th>\n      <th>ORF_Num_Gi_pc.D40<\/th>\n      <th>ORF_Num_Gi_pc.D41<\/th>\n      <th>ORF_Num_Gi_pc.D42<\/th>\n      <th>ORF_Num_Gi_pc.D43<\/th>\n      <th>ORF_Num_Gi_pc.D4<\/th>\n      <th>ORF_Num_Gi_pc.D5<\/th>\n      <th>ORF_Num_Gi_pc.D6<\/th>\n      <th>ORF_Num_Gi_pc.D7<\/th>\n      <th>ORF_Num_Gi_pc.D8<\/th>\n      <th>ORF_Num_Gi_pc.D9<\/th>\n    <\/tr>\n  <\/thead>\n<\/table>","options":{"columnDefs":[{"className":"dt-right","targets":[2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50]},{"orderable":false,"targets":0}],"order":[],"autoWidth":false,"orderClasses":false}},"evals":[],"jsHooks":[]}</script>
```


```r
full_gene_catalogue_orf_quant %>% 
  left_join(full_gene_catalogue_orf_Num_Gi_pc,
            by = c("ORF_ID" = "ORF_ID")) %>% 
  left_join(full_gene_catalogue_orf_Num_Gi,
            by = c("ORF_ID" = "ORF_ID")) -> full_gene_catalogue_orf_quant_all_metrics
```


Generate phyloseq objects:


sample_data:


```r
human_meta
```

```
##          sample input filtered denoisedF denoisedR merged tabled filtered_pc
## D1   CR-15-S292 23558    23449     23410     23355  22970  22970       0.995
## D14   CR-17-S97 35320    35021     34992     34870  34419  34419       0.992
## D21  CR-19-S341 17616    17528     17501     17477  17246  17246       0.995
## D28  CR-21-S126  5791     5764      5753      5734   5585   5585       0.995
## D35  CR-46-S191 24284    24178     24128     24062  23531  23531       0.996
## D36  CR-64-S324 10847    10801     10784     10726  10465  10465       0.996
## D43    D-0-S154 11484    11428     11279     11256  10248  10248       0.995
## D7  TR1-15-S168   153      151       144       115     79     79       0.987
## D8  TR1-17-S246 16895    16833     16807     16749  16429  16429       0.996
## D15 TR1-19-S305  9274     9256      9239      9203   8959   8959       0.998
## D22 TR1-21-S207 20083    20061     20043     19956  19419  19419       0.999
## D29 TR1-46-S170 32752    32656     32579     32459  31595  31595       0.997
## D42 TR1-64-S263 31538    31375     31337     31248  30747  30747       0.995
## D6  TR2-15-S288 31566    31486     31385     31307  30381  30381       0.997
## D9  TR2-17-S248 31490    31404     31369     31205  30516  30516       0.997
## D16 TR2-19-S318 10065    10054     10040     10015   9814   9814       0.999
## D23 TR2-21-S228 15593    15567     15543     15456  14959  14959       0.998
## D30 TR2-46-S239 25929    25888     25862     25757  25098  25098       0.998
## D41 TR2-64-S216 21707    21663     21627     21526  20850  20850       0.998
## D5  TR3-15-S338 17747    17682     17655     17595  17213  17213       0.996
## D10 TR3-17-S222 14735    14652     14631     14571  14285  14285       0.994
## D17 TR3-19-S266 32415    32239     32191     32099  31521  31521       0.995
## D24 TR3-21-S234 23216    23054     23032     22946  22495  22495       0.993
## D31 TR3-46-S224 14862    14757     14738     14685  14391  14391       0.993
## D40 TR3-64-S190 28368    28173     28084     28039  27343  27343       0.993
## D4  TR4-15-S218 21839    21709     21687     21638  21317  21317       0.994
## D11 TR4-17-S114 28993    28752     28709     28624  28097  28097       0.992
## D18 TR4-19-S270 28080    27827     27811     27758  27315  27315       0.991
## D25 TR4-21-S254 21868    21735     21727     21684  21405  21405       0.994
## D32 TR4-46-S215 25294    25090     25053     24965  24729  24729       0.992
## D39 TR4-64-S220 15761    15686     15672     15607  15483  15483       0.995
## C_1 TR5-13-S103 25912    25816     25770     25681  24848  24848       0.996
## D3  TR5-15-S358 22887    22773     22768     22682  22378  22378       0.995
## D12  TR5-17-S82 81610    81154     81069     80828  79535  79535       0.994
## D19  TR5-19-S87 35069    34934     34910     34783  34135  34135       0.996
## C_6 TR5-20-S208 24461    24365     24357     24197  23598  23598       0.996
## D26 TR5-21-S185 35805    35626     35602     35513  34882  34882       0.995
## D33 TR5-46-S310 20221    20118     20093     19985  19794  19794       0.995
## D38 TR5-64-S141 22049    21937     21916     21817  21659  21659       0.995
## D2  TR6-15-S149 21435    21335     21298     21224  20816  20816       0.995
## D13 TR6-17-S138 20193    20062     20030     19961  19569  19569       0.994
## D20 TR6-19-S137 25289    25165     25122     24996  24337  24337       0.995
## D27 TR6-21-S107  4866     4830      4824      4806   4681   4681       0.993
## D34 TR6-46-S176 49716    49525     49414     49275  48250  48250       0.996
## D37 TR6-64-S146 24475    24302     24256     24192  23730  23730       0.993
## C_2      TR5-63    NA       NA        NA        NA     NA     NA          NA
## C_3      TR4-63    NA       NA        NA        NA     NA     NA          NA
## C_4      TR3-63    NA       NA        NA        NA     NA     NA          NA
## C_5      TR2-63    NA       NA        NA        NA     NA     NA          NA
##     denoisedF_pc denoisedR_pc merged_pc filtered_merged_pc input_merged_pc
## D1         0.998        0.996     0.981              0.980           0.975
## D14        0.999        0.996     0.984              0.983           0.974
## D21        0.998        0.997     0.985              0.984           0.979
## D28        0.998        0.995     0.971              0.969           0.964
## D35        0.998        0.995     0.975              0.973           0.969
## D36        0.998        0.993     0.970              0.969           0.965
## D43        0.987        0.985     0.909              0.897           0.892
## D7         0.954        0.762     0.549              0.523           0.516
## D8         0.998        0.995     0.978              0.976           0.972
## D15        0.998        0.994     0.970              0.968           0.966
## D22        0.999        0.995     0.969              0.968           0.967
## D29        0.998        0.994     0.970              0.968           0.965
## D42        0.999        0.996     0.981              0.980           0.975
## D6         0.997        0.994     0.968              0.965           0.962
## D9         0.999        0.994     0.973              0.972           0.969
## D16        0.999        0.996     0.977              0.976           0.975
## D23        0.998        0.993     0.962              0.961           0.959
## D30        0.999        0.995     0.970              0.969           0.968
## D41        0.998        0.994     0.964              0.962           0.961
## D5         0.998        0.995     0.975              0.973           0.970
## D10        0.999        0.994     0.976              0.975           0.969
## D17        0.999        0.996     0.979              0.978           0.972
## D24        0.999        0.995     0.977              0.976           0.969
## D31        0.999        0.995     0.976              0.975           0.968
## D40        0.997        0.995     0.974              0.971           0.964
## D4         0.999        0.997     0.983              0.982           0.976
## D11        0.999        0.996     0.979              0.977           0.969
## D18        0.999        0.998     0.982              0.982           0.973
## D25        1.000        0.998     0.985              0.985           0.979
## D32        0.999        0.995     0.987              0.986           0.978
## D39        0.999        0.995     0.988              0.987           0.982
## C_1        0.998        0.995     0.964              0.963           0.959
## D3         1.000        0.996     0.983              0.983           0.978
## D12        0.999        0.996     0.981              0.980           0.975
## D19        0.999        0.996     0.978              0.977           0.973
## C_6        1.000        0.993     0.969              0.969           0.965
## D26        0.999        0.997     0.980              0.979           0.974
## D33        0.999        0.993     0.985              0.984           0.979
## D38        0.999        0.995     0.988              0.987           0.982
## D2         0.998        0.995     0.977              0.976           0.971
## D13        0.998        0.995     0.977              0.975           0.969
## D20        0.998        0.993     0.969              0.967           0.962
## D27        0.999        0.995     0.970              0.969           0.962
## D34        0.998        0.995     0.976              0.974           0.971
## D37        0.998        0.995     0.978              0.976           0.970
## C_2           NA           NA        NA                 NA              NA
## C_3           NA           NA        NA                 NA              NA
## C_4           NA           NA        NA                 NA              NA
## C_5           NA           NA        NA                 NA              NA
##     tabled_joined chimera_out length_filtered tabled_pc chimera_out_pc
## D1          22970       22102           22098         1           0.96
## D14         34419       32366           32350         1           0.94
## D21         17246       16865           16865         1           0.98
## D28          5585        5498            5498         1           0.98
## D35         23531       22823           22821         1           0.97
## D36         10465       10250           10250         1           0.98
## D43         10248       10134           10134         1           0.99
## D7             79          79              79         1           1.00
## D8          16429       16233           16233         1           0.99
## D15          8959        8815            8815         1           0.98
## D22         19419       18970           18970         1           0.98
## D29         31595       30745           30745         1           0.97
## D42         30747       30193           30193         1           0.98
## D6          30381       28832           28832         1           0.95
## D9          30516       30071           30071         1           0.99
## D16          9814        9662            9662         1           0.98
## D23         14959       14630           14630         1           0.98
## D30         25098       24475           24475         1           0.98
## D41         20850       20344           20344         1           0.98
## D5          17213       16817           16817         1           0.98
## D10         14285       14059           14059         1           0.98
## D17         31521       30971           30970         1           0.98
## D24         22495       21882           21876         1           0.97
## D31         14391       14095           14095         1           0.98
## D40         27343       26205           26202         1           0.96
## D4          21317       20936           20936         1           0.98
## D11         28097       27386           27381         1           0.97
## D18         27315       26523           26523         1           0.97
## D25         21405       20811           20805         1           0.97
## D32         24729       23673           23673         1           0.96
## D39         15483       15160           15160         1           0.98
## C_1         24848       24314           24314         1           0.98
## D3          22378       21686           21680         1           0.97
## D12         79535       76689           76678         1           0.96
## D19         34135       33512           33510         1           0.98
## C_6         23598       22590           22589         1           0.96
## D26         34882       32585           32568         1           0.93
## D33         19794       19015           19015         1           0.96
## D38         21659       21080           21080         1           0.97
## D2          20816       20429           20429         1           0.98
## D13         19569       19164           19164         1           0.98
## D20         24337       23661           23659         1           0.97
## D27          4681        4595            4595         1           0.98
## D34         48250       46576           46575         1           0.97
## D37         23730       23096           23096         1           0.97
## C_2            NA          NA              NA        NA             NA
## C_3            NA          NA              NA        NA             NA
## C_4            NA          NA              NA        NA             NA
## C_5            NA          NA              NA        NA             NA
##     length_filtered_pc Sample_description I7_Index_ID    index I5_Index_ID
## D1                   1              CR-15      N716-D ACTCGCTA      S517-D
## D14                  1              CR-17      N701-A TAAGGCGA      S513-D
## D21                  1              CR-19      N723-D TAGCGCTC      S518-D
## D28                  1              CR-21      N704-A TCCTGAGC      S520-D
## D35                  1              CR-46      N715-A ATCTCAGG      S521-D
## D36                  1              CR-64      N721-D TACGCTGC      S517-D
## D43                  1                D-0      N710-A CGAGGCTG      S515-D
## D7                   1             TR1-15      N711-A AAGAGGCA      S522-D
## D8                   1             TR1-17      N723-D TAGCGCTC      S508-A
## D15                  1             TR1-19      N719-D GCGTAGTA      S513-D
## D22                  1             TR1-21      N718-D GGAGCTAC      S510-A
## D29                  1             TR1-46      N712-A GTAGAGGA      S515-D
## D42                  1             TR1-64      N726-D CCTAAGAC      S510-A
## D6                   1             TR2-15      N729-D TCGACGTC      S511-A
## D9                   1             TR2-17      N723-D TAGCGCTC      S511-A
## D16                  1             TR2-19      N720-D CGGAGCCT      S520-D
## D23                  1             TR2-21      N721-D TACGCTGC      S506-A
## D30                  1             TR2-46      N722-D ATGCGCAG      S510-A
## D41                  1             TR2-64      N719-D GCGTAGTA      S511-A
## D5                   1             TR3-15      N723-D TAGCGCTC      S515-D
## D10                  1             TR3-17      N720-D CGGAGCCT      S508-A
## D17                  1             TR3-19      N727-D CGATCAGT      S503-A
## D24                  1             TR3-21      N722-D ATGCGCAG      S503-A
## D31                  1             TR3-46      N720-D CGGAGCCT      S511-A
## D40                  1             TR3-64      N715-A ATCTCAGG      S520-D
## D4                   1             TR4-15      N720-D CGGAGCCT      S503-A
## D11                  1             TR4-17      N703-A AGGCAGAA      S515-D
## D18                  1             TR4-19      N727-D CGATCAGT      S508-A
## D25                  1             TR4-21      N724-D ACTGAGCG      S508-A
## D32                  1             TR4-46      N719-D GCGTAGTA      S510-A
## D39                  1             TR4-64      N720-D CGGAGCCT      S506-A
## C_1                  1             TR5-13      N701-A TAAGGCGA      S521-D
## D3                   1             TR5-15      N726-D CCTAAGAC      S520-D
## D12                  1             TR5-17      N714-A GCTCATGA      S503-A
## D19                  1             TR5-19      N714-A GCTCATGA      S510-A
## C_6                  1             TR5-20      N718-D GGAGCTAC      S511-A
## D26                  1             TR5-21      N715-A ATCTCAGG      S513-D
## D33                  1             TR5-46      N719-D GCGTAGTA      S520-D
## D38                  1             TR5-64      N706-A TAGGCATG      S518-D
## D2                   1             TR6-15      N707-A CTCTCTAC      S518-D
## D13                  1             TR6-17      N706-A TAGGCATG      S515-D
## D20                  1             TR6-19      N706-A TAGGCATG      S513-D
## D27                  1             TR6-21      N702-A CGTACTAG      S516-D
## D34                  1             TR6-46      N712-A GTAGAGGA      S522-D
## D37                  1             TR6-64      N707-A CTCTCTAC      S515-D
## C_2                 NA             TR5-63        <NA>     <NA>        <NA>
## C_3                 NA             TR4-63        <NA>     <NA>        <NA>
## C_4                 NA             TR3-63        <NA>     <NA>        <NA>
## C_5                 NA             TR2-63        <NA>     <NA>        <NA>
##       index2 Description2 Experiment Reactor     Treatment Day_of_Connection
## D1  GCGTAAGA         <NA> Continuous      CR     UNTREATED                15
## D14 TCGACTAG         <NA> Continuous      CR     UNTREATED                17
## D21 CTATTAAG         <NA> Continuous      CR     UNTREATED                19
## D28 AAGGCTAT         <NA> Continuous      CR     UNTREATED                21
## D35 GAGCCTTA         <NA> Continuous      CR     UNTREATED                46
## D36 GCGTAAGA         <NA> Continuous      CR     UNTREATED                64
## D43 TTCTAGCT        DONOR      Cecum   DONOR         DONOR                NA
## D7  TTATGCGA         <NA> Continuous     TR1   CTX+HV292.1                15
## D8  CTAAGCCT         <NA> Continuous     TR1   CTX+HV292.1                17
## D15 TCGACTAG         <NA> Continuous     TR1   CTX+HV292.1                19
## D22 CGTCTAAT         <NA> Continuous     TR1   CTX+HV292.1                21
## D29 TTCTAGCT         <NA> Continuous     TR1   CTX+HV292.1                46
## D42 CGTCTAAT         <NA> Continuous     TR1   CTX+HV292.1                64
## D6  TCTCTCCG         <NA> Continuous     TR2           CTX                15
## D9  TCTCTCCG         <NA> Continuous     TR2           CTX                17
## D16 AAGGCTAT         <NA> Continuous     TR2           CTX                19
## D23 ACTGCATA         <NA> Continuous     TR2           CTX                21
## D30 CGTCTAAT         <NA> Continuous     TR2           CTX                46
## D41 TCTCTCCG         <NA> Continuous     TR2           CTX                64
## D5  TTCTAGCT         <NA> Continuous     TR3       HV292.1                15
## D10 CTAAGCCT         <NA> Continuous     TR3       HV292.1                17
## D17 TATCCTCT         <NA> Continuous     TR3       HV292.1                19
## D24 TATCCTCT         <NA> Continuous     TR3       HV292.1                21
## D31 TCTCTCCG         <NA> Continuous     TR3       HV292.1                46
## D40 AAGGCTAT         <NA> Continuous     TR3       HV292.1                64
## D4  TATCCTCT         <NA> Continuous     TR4           VAN                15
## D11 TTCTAGCT         <NA> Continuous     TR4           VAN                17
## D18 CTAAGCCT         <NA> Continuous     TR4           VAN                19
## D25 CTAAGCCT         <NA> Continuous     TR4           VAN                21
## D32 CGTCTAAT         <NA> Continuous     TR4           VAN                46
## D39 ACTGCATA         <NA> Continuous     TR4           VAN                64
## C_1 GAGCCTTA         <NA> Continuous     TR5 VAN+CCUG59168                13
## D3  AAGGCTAT         <NA> Continuous     TR5 VAN+CCUG59168                15
## D12 TATCCTCT         <NA> Continuous     TR5 VAN+CCUG59168                17
## D19 CGTCTAAT         <NA> Continuous     TR5 VAN+CCUG59168                19
## C_6 TCTCTCCG         <NA> Continuous     TR5 VAN+CCUG59168                20
## D26 TCGACTAG         <NA> Continuous     TR5 VAN+CCUG59168                21
## D33 AAGGCTAT         <NA> Continuous     TR5 VAN+CCUG59168                46
## D38 CTATTAAG         <NA> Continuous     TR5 VAN+CCUG59168                64
## D2  CTATTAAG         <NA> Continuous     TR6     CCUG59168                15
## D13 TTCTAGCT         <NA> Continuous     TR6     CCUG59168                17
## D20 TCGACTAG         <NA> Continuous     TR6     CCUG59168                19
## D27 CCTAGAGT         <NA> Continuous     TR6     CCUG59168                21
## D34 TTATGCGA         <NA> Continuous     TR6     CCUG59168                46
## D37 TTCTAGCT         <NA> Continuous     TR6     CCUG59168                64
## C_2     <NA>         <NA> Continuous     TR5 VAN+CCUG59168                63
## C_3     <NA>         <NA> Continuous     TR4           VAN                63
## C_4     <NA>         <NA> Continuous     TR3       HV292.1                63
## C_5     <NA>         <NA> Continuous     TR2           CTX                63
##     Day_of_Treatment Day_from_Inoculum  Enrichment Phase    Treatment2
## D1                -1                38 NotEnriched  Stab     UNTREATED
## D14                1                40 NotEnriched Treat     UNTREATED
## D21                3                42 NotEnriched Treat     UNTREATED
## D28                5                44 NotEnriched Treat     UNTREATED
## D35               30                69 NotEnriched Treat     UNTREATED
## D36               48                87 NotEnriched Treat     UNTREATED
## D43             <NA>                NA NotEnriched DONOR         DONOR
## D7                -1                38 NotEnriched  Stab    AB+E. coli
## D8                 1                40 NotEnriched Treat    AB+E. coli
## D15                3                42 NotEnriched Treat    AB+E. coli
## D22                5                44 NotEnriched Treat    AB+E. coli
## D29               30                69 NotEnriched Treat    AB+E. coli
## D42               48                87 NotEnriched Treat    AB+E. coli
## D6                -1                38 NotEnriched  Stab            AB
## D9                 1                40 NotEnriched Treat            AB
## D16                3                42 NotEnriched Treat            AB
## D23                5                44 NotEnriched Treat            AB
## D30               30                69 NotEnriched Treat            AB
## D41               48                87 NotEnriched Treat            AB
## D5                -1                38 NotEnriched  Stab       E. coli
## D10                1                40 NotEnriched Treat       E. coli
## D17                3                42 NotEnriched Treat       E. coli
## D24                5                44 NotEnriched Treat       E. coli
## D31               30                69 NotEnriched Treat       E. coli
## D40               48                87 NotEnriched Treat       E. coli
## D4                -1                38 NotEnriched  Stab            AB
## D11                1                40 NotEnriched Treat            AB
## D18                3                42 NotEnriched Treat            AB
## D25                5                44 NotEnriched Treat            AB
## D32               30                69 NotEnriched Treat            AB
## D39               48                87 NotEnriched Treat            AB
## C_1               -3                36 NotEnriched  Stab AB+E. faecium
## D3                -1                38 NotEnriched  Stab AB+E. faecium
## D12                1                40 NotEnriched Treat AB+E. faecium
## D19                3                42 NotEnriched Treat AB+E. faecium
## C_6                4                43 NotEnriched Treat AB+E. faecium
## D26                5                44 NotEnriched Treat AB+E. faecium
## D33               30                69 NotEnriched Treat AB+E. faecium
## D38               48                87 NotEnriched Treat AB+E. faecium
## D2                -1                38 NotEnriched  Stab    E. faecium
## D13                1                40 NotEnriched Treat    E. faecium
## D20                3                42 NotEnriched Treat    E. faecium
## D27                5                44 NotEnriched Treat    E. faecium
## D34               30                69 NotEnriched Treat    E. faecium
## D37               48                87 NotEnriched Treat    E. faecium
## C_2               47                86 NotEnriched Treat AB+E. faecium
## C_3               47                86 NotEnriched Treat            AB
## C_4               47                86 NotEnriched Treat       E. coli
## C_5               47                86 NotEnriched Treat            AB
##                     Date Paul Reactor_Treatment GeneCopyNumberperML
## D1  2020-06-26T00:00:00Z <NA>      CR_UNTREATED            6.75e+10
## D14 2020-06-28T00:00:00Z <NA>      CR_UNTREATED            7.01e+10
## D21 2020-06-30T00:00:00Z <NA>      CR_UNTREATED            9.01e+10
## D28 2020-07-02T00:00:00Z <NA>      CR_UNTREATED            4.39e+10
## D35 2020-07-27T00:00:00Z <NA>      CR_UNTREATED            4.53e+10
## D36 2020-08-14T00:00:00Z <NA>      CR_UNTREATED            3.26e+10
## D43 2020-05-19T00:00:00Z <NA>             DONOR            1.44e+10
## D7  2020-06-26T00:00:00Z <NA>   TR1_CTX+HV292.1            5.08e+10
## D8  2020-06-28T00:00:00Z <NA>   TR1_CTX+HV292.1            3.11e+10
## D15 2020-06-30T00:00:00Z <NA>   TR1_CTX+HV292.1            3.36e+10
## D22 2020-07-02T00:00:00Z <NA>   TR1_CTX+HV292.1            2.22e+10
## D29 2020-07-27T00:00:00Z <NA>   TR1_CTX+HV292.1            6.05e+10
## D42 2020-08-14T00:00:00Z <NA>   TR1_CTX+HV292.1            3.61e+10
## D6  2020-06-26T00:00:00Z <NA>           TR2_CTX            6.23e+10
## D9  2020-06-28T00:00:00Z <NA>           TR2_CTX            2.71e+10
## D16 2020-06-30T00:00:00Z <NA>           TR2_CTX            4.20e+10
## D23 2020-07-02T00:00:00Z <NA>           TR2_CTX            2.89e+10
## D30 2020-07-27T00:00:00Z <NA>           TR2_CTX            5.05e+10
## D41 2020-08-14T00:00:00Z <NA>           TR2_CTX            2.79e+10
## D5  2020-06-26T00:00:00Z <NA>       TR3_HV292.1            4.87e+10
## D10 2020-06-28T00:00:00Z <NA>       TR3_HV292.1            5.24e+10
## D17 2020-06-30T00:00:00Z <NA>       TR3_HV292.1            4.62e+10
## D24 2020-07-02T00:00:00Z <NA>       TR3_HV292.1            4.79e+10
## D31 2020-07-27T00:00:00Z <NA>       TR3_HV292.1            3.84e+10
## D40 2020-08-14T00:00:00Z <NA>       TR3_HV292.1            3.43e+10
## D4  2020-06-26T00:00:00Z <NA>           TR4_VAN            6.68e+10
## D11 2020-06-28T00:00:00Z <NA>           TR4_VAN            3.74e+10
## D18 2020-06-30T00:00:00Z <NA>           TR4_VAN            5.61e+10
## D25 2020-07-02T00:00:00Z <NA>           TR4_VAN            5.82e+10
## D32 2020-07-27T00:00:00Z <NA>           TR4_VAN            6.31e+10
## D39 2020-08-14T00:00:00Z <NA>           TR4_VAN            6.94e+10
## C_1 2020-06-24T00:00:00Z <NA> TR5_VAN+CCUG59168            3.94e+10
## D3  2020-06-26T00:00:00Z <NA> TR5_VAN+CCUG59168            7.15e+10
## D12 2020-06-28T00:00:00Z <NA> TR5_VAN+CCUG59168            2.60e+10
## D19 2020-06-30T00:00:00Z <NA> TR5_VAN+CCUG59168            5.71e+10
## C_6 2020-07-01T00:00:00Z <NA> TR5_VAN+CCUG59168            2.50e+10
## D26 2020-07-02T00:00:00Z <NA> TR5_VAN+CCUG59168            4.89e+10
## D33 2020-07-27T00:00:00Z <NA> TR5_VAN+CCUG59168            6.59e+10
## D38 2020-08-14T00:00:00Z <NA> TR5_VAN+CCUG59168            5.32e+10
## D2  2020-06-26T00:00:00Z <NA>     TR6_CCUG59168            5.30e+10
## D13 2020-06-28T00:00:00Z <NA>     TR6_CCUG59168            3.93e+10
## D20 2020-06-30T00:00:00Z <NA>     TR6_CCUG59168            5.68e+10
## D27 2020-07-02T00:00:00Z <NA>     TR6_CCUG59168            6.93e+10
## D34 2020-07-27T00:00:00Z <NA>     TR6_CCUG59168            2.81e+10
## D37 2020-08-14T00:00:00Z <NA>     TR6_CCUG59168            3.93e+10
## C_2             13.08.20 <NA> TR5_VAN+CCUG59168                  NA
## C_3             13.08.20 <NA>           TR4_VAN                  NA
## C_4             13.08.20 <NA>       TR3_HV292.1                  NA
## C_5             13.08.20 <NA>           TR2_CTX                  NA
##     HV292.1_Copy_Number_permL CCUG59168_Copy_Number_permL CTX_Copy_Number_permL
## D1                         NA                          NA                    NA
## D14                        NA                          NA                    NA
## D21                        NA                          NA                    NA
## D28                        NA                          NA                    NA
## D35                        NA                          NA                    NA
## D36                        NA                          NA                    NA
## D43                        NA                          NA                    NA
## D7                         NA                          NA                    NA
## D8                         NA                          NA                    NA
## D15                    794000                          NA          1.507417e+06
## D22                        NA                          NA                    NA
## D29                  54054532                          NA          8.070878e+07
## D42                 862908750                          NA          1.927668e+09
## D6                         NA                          NA                    NA
## D9                         NA                          NA                    NA
## D16                        NA                          NA                    NA
## D23                        NA                          NA                    NA
## D30                        NA                          NA                    NA
## D41                        NA                          NA                    NA
## D5                         NA                          NA                    NA
## D10                        NA                          NA                    NA
## D17                    887000                          NA          1.670000e+06
## D24                        NA                          NA                    NA
## D31                         0                          NA          4.713665e+03
## D40                         0                          NA          3.666540e+04
## D4                         NA                          NA                    NA
## D11                        NA                          NA                    NA
## D18                        NA                          NA                    NA
## D25                        NA                          NA                    NA
## D32                        NA                          NA                    NA
## D39                        NA                          NA                    NA
## C_1                        NA                          NA                    NA
## D3                         NA                          NA                    NA
## D12                        NA                          NA                    NA
## D19                        NA                  50400000.0                    NA
## C_6                        NA                          NA                    NA
## D26                        NA                          NA                    NA
## D33                        NA                  91984000.0                    NA
## D38                        NA                 539000000.0                    NA
## D2                         NA                          NA                    NA
## D13                        NA                          NA                    NA
## D20                        NA                         0.0                    NA
## D27                        NA                          NA                    NA
## D34                        NA                    268816.8                    NA
## D37                        NA                      2120.0                    NA
## C_2                        NA                          NA                    NA
## C_3                        NA                          NA                    NA
## C_4                        NA                          NA                    NA
## C_5                        NA                          NA                    NA
##     VAN_Copy_Number_permL   Model Antibiotic_mg.mL Fermentation Antibiotic
## D1                     NA Chicken               NA           NA       <NA>
## D14                    NA Chicken               NA           NA       <NA>
## D21                    NA Chicken               NA           NA       <NA>
## D28                    NA Chicken               NA           NA       <NA>
## D35                    NA Chicken               NA           NA       <NA>
## D36                    NA Chicken               NA           NA       <NA>
## D43                    NA Chicken               NA           NA       <NA>
## D7                     NA Chicken               20           NA        CTX
## D8                     NA Chicken               20           NA        CTX
## D15                    NA Chicken               20           NA        CTX
## D22                    NA Chicken               20           NA        CTX
## D29                    NA Chicken               20           NA        CTX
## D42                    NA Chicken               20           NA        CTX
## D6                     NA Chicken               20           NA        CTX
## D9                     NA Chicken               20           NA        CTX
## D16                    NA Chicken               20           NA        CTX
## D23                    NA Chicken               20           NA        CTX
## D30                    NA Chicken               20           NA        CTX
## D41                    NA Chicken               20           NA        CTX
## D5                     NA Chicken               NA           NA       <NA>
## D10                    NA Chicken               NA           NA       <NA>
## D17                    NA Chicken               NA           NA       <NA>
## D24                    NA Chicken               NA           NA       <NA>
## D31                    NA Chicken               NA           NA       <NA>
## D40                    NA Chicken               NA           NA       <NA>
## D4                     NA Chicken               90           NA        VAN
## D11                    NA Chicken               90           NA        VAN
## D18                    NA Chicken               90           NA        VAN
## D25                    NA Chicken               90           NA        VAN
## D32                    NA Chicken               90           NA        VAN
## D39                    NA Chicken               90           NA        VAN
## C_1                    NA Chicken               90           NA        VAN
## D3            17402500000 Chicken               90           NA        VAN
## D12           41564875000 Chicken               90           NA        VAN
## D19           66827950000 Chicken               90           NA        VAN
## C_6           28473350000 Chicken               90           NA        VAN
## D26           57697500000 Chicken               90           NA        VAN
## D33          577000000000 Chicken               90           NA        VAN
## D38                    NA Chicken               90           NA        VAN
## D2              106535000 Chicken               NA           NA       <NA>
## D13                     0 Chicken               NA           NA       <NA>
## D20                     0 Chicken               NA           NA       <NA>
## D27                     0 Chicken               NA           NA       <NA>
## D34                    NA Chicken               NA           NA       <NA>
## D37                    NA Chicken               NA           NA       <NA>
## C_2                    NA Chicken               90           NA        VAN
## C_3                    NA Chicken               90           NA        VAN
## C_4                    NA Chicken               NA           NA        CTX
## C_5                    NA Chicken               20           NA       <NA>
##     Lactose_mM Glucose_mM Galactose_mM Succinat_mM Lactat_mM Formiat_mM
## D1       0.000      0.947        0.000       1.969     0.000      3.179
## D14      0.000      0.000        1.759       4.295     8.223      5.614
## D21      0.261      0.000        0.000       0.818     4.327      0.000
## D28      0.000      0.000        0.000       0.701     0.000      0.000
## D35         NA         NA           NA          NA        NA         NA
## D36      0.000      0.000        0.000       0.546     0.000      0.000
## D43         NA         NA           NA          NA        NA         NA
## D7       0.000      0.000        0.000       1.090     0.000      0.000
## D8       0.000      0.000        0.000       1.201     0.000      0.000
## D15      0.000      0.000        0.000       2.779     0.000      0.000
## D22      0.000      0.000        0.000       3.069     0.000      0.000
## D29      0.000      0.000        0.000       0.527     8.115      0.000
## D42      0.000      0.000        0.000       0.455     1.035      0.000
## D6       0.000      0.000        0.000       0.792     0.000      0.000
## D9       0.000      0.000        0.000       1.116     0.000      0.000
## D16      0.000      0.000        0.000       2.973     0.000      0.000
## D23      0.000      0.000        0.000       3.420     0.000      0.000
## D30      0.000      0.000        0.000       3.776     0.000      0.000
## D41      0.000      0.000        0.000       0.511     0.000      2.219
## D5       0.000      0.000        0.000       0.650     0.000      0.000
## D10      0.000      0.000        0.000       0.751     0.000      0.000
## D17      0.000      0.000        0.000       0.629     0.000      0.000
## D24      0.000      0.000        0.000       0.639     0.000      0.000
## D31      0.000      0.000        0.000       0.453     0.000      0.000
## D40      0.000      0.000        0.000       0.456     0.000      0.000
## D4       0.000      0.000        0.000       0.621     0.000      0.000
## D11      0.000      1.401        1.045       4.792     0.000      0.000
## D18      0.000      0.000        1.380      11.605     0.000      0.000
## D25      0.598      0.000        1.033      10.796     0.000      9.618
## D32      0.000      0.000        0.000       7.444     8.524      0.000
## D39      0.000      0.000        0.000       0.418     3.808      0.000
## C_1      0.000      0.000        0.000       1.525     0.000      0.000
## D3       0.000      0.000        0.000       0.694     0.000      0.000
## D12      0.000      1.584        1.135       6.021     0.000      0.000
## D19      0.489      0.000        1.448       9.526     0.000      0.000
## C_6      0.550      0.437        0.000       7.267     0.000      0.000
## D26      0.564      0.000        0.000       9.386     0.000      9.711
## D33      0.000      0.000        0.000       7.771     0.000      0.000
## D38      0.000      0.000        0.000       3.990     0.000      0.000
## D2       0.000      0.000        0.000       0.734     0.000      0.000
## D13      0.000      0.000        0.000       0.726     0.000      0.000
## D20      0.000      0.000        0.000       0.748     0.000      0.000
## D27      0.000      0.000        0.000       0.729     0.000      0.000
## D34      0.000      0.000        0.000       1.150     0.000      0.000
## D37      0.000      0.000        0.000       0.429     0.000      0.000
## C_2         NA         NA           NA          NA        NA         NA
## C_3         NA         NA           NA          NA        NA         NA
## C_4         NA         NA           NA          NA        NA         NA
## C_5         NA         NA           NA          NA        NA         NA
##     Acetat_mM Propionat_mM Isobutyrat_mM Butyrat_mM Isovalerat_mM Valerat_mM
## D1     49.775       11.560         5.270     22.600         4.267      3.218
## D14    39.530        7.362         0.000      8.315         1.417      0.365
## D21    81.797       10.133         0.000     14.599         0.881      0.000
## D28    86.705       15.767         6.942     27.682         4.870      0.000
## D35        NA           NA            NA         NA            NA         NA
## D36    75.009       11.518         7.914     36.797         4.585      6.461
## D43        NA           NA            NA         NA            NA         NA
## D7     83.366       13.060         8.268     41.333         6.524      7.962
## D8     82.428       12.110         7.901     46.439         6.534      4.832
## D15    77.988        8.582         4.400     49.923         3.876      1.147
## D22    72.326       12.140         0.000     49.287         7.405      0.522
## D29   100.476       13.971         7.021     59.817         8.804      9.447
## D42    69.191       13.707         7.787     36.459         7.614      8.207
## D6     98.070       13.180         8.776     44.333         6.422      7.313
## D9     95.235       12.018         7.794     47.890         6.038      3.410
## D16    90.991        7.423         0.000     44.651         2.036      0.036
## D23    84.906        5.716         0.000     41.240         1.044      0.000
## D30   110.663       12.096         5.428     44.934         6.609      0.000
## D41    57.910       13.311         7.440     35.756         6.925      0.233
## D5     99.568       14.760         7.584     33.721         6.467      8.363
## D10    86.398       11.689         8.204     35.082         6.429      6.960
## D17    85.603       13.331         6.860     29.604         6.109      0.000
## D24    93.635       13.645         6.714     26.799         5.297      0.000
## D31   105.354       14.368         7.536     36.305         6.443      7.805
## D40    83.474       12.933         6.417     23.882         4.188      6.227
## D4    101.441       14.850         7.813     33.872         6.542      8.752
## D11    79.243       17.005         8.207     28.479         6.182      4.307
## D18    48.991       21.902         4.715      6.737         5.970      0.436
## D25    47.637       24.239         0.000      1.522         5.250      0.000
## D32    46.362       19.358         4.941      0.287         7.161      0.000
## D39    41.394       17.022         5.820      2.230         4.591      0.000
## C_1    92.696       11.958         9.009     46.755         6.637      7.112
## D3     94.605       12.469         7.545     42.811         6.110      7.167
## D12    77.564       16.918         7.590     29.690         5.858      3.135
## D19    47.530       23.120         0.000      7.136         5.905      0.337
## C_6    44.220       21.102         0.000      3.984         5.818      0.000
## D26    45.669       24.147         0.000      1.788         5.231      0.000
## D33    41.884       18.412         5.257      1.018         7.101      0.000
## D38    31.082       13.931         6.202      1.974         4.674      0.000
## D2     97.498       11.455         7.408     47.083         5.748      7.114
## D13    87.252       11.533         8.115     38.443         6.024      6.447
## D20    87.310       13.946         0.000     31.669         6.040      0.000
## D27    94.750       13.997         7.214     26.634         5.309      0.000
## D34   100.967       15.176         7.531     37.422         7.673      7.628
## D37    85.485       14.118         6.381     23.281         5.541      6.703
## C_2        NA           NA            NA         NA            NA         NA
## C_3        NA           NA            NA         NA            NA         NA
## C_4        NA           NA            NA         NA            NA         NA
## C_5        NA           NA            NA         NA            NA         NA
##     Total_SCFA_mM raw_metagenomic_pairs Reactor_Treatment_Dose  Treatment_Dose
## D1        101.838              57121568           CR_UNTREATED       UNTREATED
## D14        76.880              44863798           CR_UNTREATED       UNTREATED
## D21       112.555              56304446           CR_UNTREATED       UNTREATED
## D28       142.667              50088120           CR_UNTREATED       UNTREATED
## D35            NA              55285758           CR_UNTREATED       UNTREATED
## D36       142.830              48021514           CR_UNTREATED       UNTREATED
## D43            NA              68282190                  DONOR           DONOR
## D7        161.603              52833332      TR1_CTX+HV292.120   CTX+HV292.120
## D8        161.445              53009028      TR1_CTX+HV292.120   CTX+HV292.120
## D15       148.695              47492202      TR1_CTX+HV292.120   CTX+HV292.120
## D22       144.749              55638128      TR1_CTX+HV292.120   CTX+HV292.120
## D29       208.178              68334190      TR1_CTX+HV292.120   CTX+HV292.120
## D42       144.455              47906582      TR1_CTX+HV292.120   CTX+HV292.120
## D6        178.886              51668602              TR2_CTX20           CTX20
## D9        173.501              47401086              TR2_CTX20           CTX20
## D16       148.110              44839890              TR2_CTX20           CTX20
## D23       136.326              56766854              TR2_CTX20           CTX20
## D30       183.506              49315184              TR2_CTX20           CTX20
## D41       124.305              59842344              TR2_CTX20           CTX20
## D5        171.113              55615446            TR3_HV292.1         HV292.1
## D10       155.513              46820988            TR3_HV292.1         HV292.1
## D17       142.136              49791554            TR3_HV292.1         HV292.1
## D24       146.729              48166738            TR3_HV292.1         HV292.1
## D31       178.264              56873318            TR3_HV292.1         HV292.1
## D40       137.577              61210784            TR3_HV292.1         HV292.1
## D4        173.891              50262020              TR4_VAN90           VAN90
## D11       148.215              55159700              TR4_VAN90           VAN90
## D18       100.356              42389522              TR4_VAN90           VAN90
## D25        99.062              49858296              TR4_VAN90           VAN90
## D32        85.553              50291070              TR4_VAN90           VAN90
## D39        71.475              64514742              TR4_VAN90           VAN90
## C_1       175.692                    NA    TR5_VAN+CCUG5916890 VAN+CCUG5916890
## D3        171.401              54057756    TR5_VAN+CCUG5916890 VAN+CCUG5916890
## D12       146.776              45384782    TR5_VAN+CCUG5916890 VAN+CCUG5916890
## D19        93.554              54226192    TR5_VAN+CCUG5916890 VAN+CCUG5916890
## C_6        82.391                    NA    TR5_VAN+CCUG5916890 VAN+CCUG5916890
## D26        86.221              50932862    TR5_VAN+CCUG5916890 VAN+CCUG5916890
## D33        81.443              59788084    TR5_VAN+CCUG5916890 VAN+CCUG5916890
## D38        61.853              38395332    TR5_VAN+CCUG5916890 VAN+CCUG5916890
## D2        177.040              60148698          TR6_CCUG59168       CCUG59168
## D13       158.540              52085624          TR6_CCUG59168       CCUG59168
## D20       139.713              46655108          TR6_CCUG59168       CCUG59168
## D27       148.633              70694598          TR6_CCUG59168       CCUG59168
## D34       177.547              56546096          TR6_CCUG59168       CCUG59168
## D37       141.938              51876272          TR6_CCUG59168       CCUG59168
## C_2            NA                    NA    TR5_VAN+CCUG5916890 VAN+CCUG5916890
## C_3            NA                    NA              TR4_VAN90           VAN90
## C_4            NA                    NA            TR3_HV292.1         HV292.1
## C_5            NA                    NA              TR2_CTX20           CTX20
##     Period
## D1    pret
## D14      1
## D21      3
## D28      5
## D35     30
## D36     48
## D43   pret
## D7    pret
## D8       1
## D15      3
## D22      5
## D29     30
## D42     48
## D6    pret
## D9       1
## D16      3
## D23      5
## D30     30
## D41     48
## D5    pret
## D10      1
## D17      3
## D24      5
## D31     30
## D40     48
## D4    pret
## D11      1
## D18      3
## D25      5
## D32     30
## D39     48
## C_1   pret
## D3    pret
## D12      1
## D19      3
## C_6      4
## D26      5
## D33     30
## D38     48
## D2    pret
## D13      1
## D20      3
## D27      5
## D34     30
## D37     48
## C_2     47
## C_3     47
## C_4     47
## C_5     47
```

```r
human_meta %>%  rownames()
```

```
##  [1] "D1"  "D14" "D21" "D28" "D35" "D36" "D43" "D7"  "D8"  "D15" "D22" "D29"
## [13] "D42" "D6"  "D9"  "D16" "D23" "D30" "D41" "D5"  "D10" "D17" "D24" "D31"
## [25] "D40" "D4"  "D11" "D18" "D25" "D32" "D39" "C_1" "D3"  "D12" "D19" "C_6"
## [37] "D26" "D33" "D38" "D2"  "D13" "D20" "D27" "D34" "D37" "C_2" "C_3" "C_4"
## [49] "C_5"
```


```r
rownames(human_meta) <- str_replace(rownames(human_meta),
                          "_", "")
```

tax_table:


```r
full_gene_catalogue_orf_quant_all_metrics %>% 
  column_to_rownames("ORF_ID") %>% 
  select(Contig_ID:ANVIO_CONCOCT_HQ) %>% 
  as.matrix() %>% 
  tax_table() -> tax
```


ORFs:

TPM:



```r
full_gene_catalogue_orf_quant_all_metrics %>% 
  column_to_rownames("ORF_ID") %>% 
  select(starts_with("ORF_TPM.")) %>% 
  rename_at(
    # select all variables with "factor.sample" in the name
    vars(contains("ORF_TPM."))
    # use stringr::str_replace to remove factor.sample.
    #   you could do the same with base::gsub()
    , funs(str_replace(., "ORF_TPM.", ""))) %>% 
  # stringr::str_replace_all("\\s","_")
  rename_at(
    # select all variables with "factor.sample" in the name
    vars(contains("."))
    , funs(str_replace_all(., "\\.", "_")))  %>% 
  as.matrix() %>% 
  otu_table(taxa_are_rows = TRUE) -> tpm
```


Num_Gi_pc


```r
full_gene_catalogue_orf_quant_all_metrics %>% 
  column_to_rownames("ORF_ID") %>% 
  select(starts_with("ORF_Num_Gi_pc.")) %>% 
  rename_at(
    # select all variables with "factor.sample" in the name
    vars(contains("ORF_Num_Gi_pc."))
    # use stringr::str_replace to remove factor.sample.
    #   you could do the same with base::gsub()
    , funs(str_replace(., "ORF_Num_Gi_pc.", ""))) %>% 
  # stringr::str_replace_all("\\s","_")
  rename_at(
    # select all variables with "factor.sample" in the name
    vars(contains("."))
    , funs(str_replace_all(., "\\.", "_")))  %>% 
  as.matrix() %>% 
  otu_table(taxa_are_rows = TRUE) -> Num_Gi_pc
```

Num_Gi


```r
full_gene_catalogue_orf_quant_all_metrics %>% 
  column_to_rownames("ORF_ID") %>% 
  select(starts_with("ORF_Num_Gi.")) %>% 
  rename_at(
    # select all variables with "factor.sample" in the name
    vars(contains("ORF_Num_Gi."))
    # use stringr::str_replace to remove factor.sample.
    #   you could do the same with base::gsub()
    , funs(str_replace(., "ORF_Num_Gi.", ""))) %>% 
  # stringr::str_replace_all("\\s","_")
  rename_at(
    # select all variables with "factor.sample" in the name
    vars(contains("."))
    , funs(str_replace_all(., "\\.", "_")))  %>% 
  as.matrix() %>% 
  otu_table(taxa_are_rows = TRUE) -> Num_Gi
```

Coverage


```r
full_gene_catalogue_orf_quant_all_metrics %>% 
  column_to_rownames("ORF_ID") %>% 
  select(starts_with("ORF_Coverage.")) %>% 
  rename_at(
    # select all variables with "factor.sample" in the name
    vars(contains("ORF_Coverage."))
    # use stringr::str_replace to remove factor.sample.
    #   you could do the same with base::gsub()
    , funs(str_replace(., "ORF_Coverage.", ""))) %>% 
  # stringr::str_replace_all("\\s","_")
  rename_at(
    # select all variables with "factor.sample" in the name
    vars(contains("."))
    , funs(str_replace_all(., "\\.", "_")))  %>% 
  as.matrix() %>% 
  otu_table(taxa_are_rows = TRUE) -> cov
```

ORF_Raw.read


```r
full_gene_catalogue_orf_quant_all_metrics %>% 
  column_to_rownames("ORF_ID") %>% 
  select(starts_with("ORF_Raw.read.count.")) %>% 
  rename_at(
    # select all variables with "factor.sample" in the name
    vars(contains("ORF_Raw.read.count."))
    # use stringr::str_replace to remove factor.sample.
    #   you could do the same with base::gsub()
    , funs(str_replace(., "ORF_Raw.read.count.", ""))) %>% 
  # stringr::str_replace_all("\\s","_")
  rename_at(
    # select all variables with "factor.sample" in the name
    vars(contains("."))
    , funs(str_replace_all(., "\\.", "_")))  %>% 
  as.matrix() %>% 
  otu_table(taxa_are_rows = TRUE) -> read_count
```



Create phyloseq objects and tsv files:


```r
dim(tpm)
```

```
## [1] 1634985      49
```

```r
physeq_tpm <- phyloseq(tpm,
                       tax,
                       human_meta %>% sample_data())


physeq_tpm
```

```
## phyloseq-class experiment-level object
## otu_table()   OTU Table:         [ 1634985 taxa and 49 samples ]
## sample_data() Sample Data:       [ 49 samples by 63 sample variables ]
## tax_table()   Taxonomy Table:    [ 1634985 taxa by 99 taxonomic ranks ]
```


```r
dim(cov)
```

```
## [1] 1634985      49
```

```r
physeq_cov <- phyloseq(cov,
                       tax,
                       human_meta %>% sample_data())

physeq_cov
```

```
## phyloseq-class experiment-level object
## otu_table()   OTU Table:         [ 1634985 taxa and 49 samples ]
## sample_data() Sample Data:       [ 49 samples by 63 sample variables ]
## tax_table()   Taxonomy Table:    [ 1634985 taxa by 99 taxonomic ranks ]
```



```r
dim(Num_Gi_pc)
```

```
## [1] 1634985      49
```

```r
physeq_Num_Gi_pc <- phyloseq(Num_Gi_pc,
                             tax,
                             human_meta %>% sample_data())

physeq_Num_Gi_pc
```

```
## phyloseq-class experiment-level object
## otu_table()   OTU Table:         [ 1634985 taxa and 49 samples ]
## sample_data() Sample Data:       [ 49 samples by 63 sample variables ]
## tax_table()   Taxonomy Table:    [ 1634985 taxa by 99 taxonomic ranks ]
```


```r
dim(Num_Gi)
```

```
## [1] 1634985      49
```

```r
physeq_Num_Gi<- phyloseq(Num_Gi,
                         tax,
                         human_meta %>% sample_data())

physeq_Num_Gi
```

```
## phyloseq-class experiment-level object
## otu_table()   OTU Table:         [ 1634985 taxa and 49 samples ]
## sample_data() Sample Data:       [ 49 samples by 63 sample variables ]
## tax_table()   Taxonomy Table:    [ 1634985 taxa by 99 taxonomic ranks ]
```



```r
dim(read_count)
```

```
## [1] 1634985      49
```

```r
physeq_read_count <- phyloseq(read_count,
                              tax,
                              human_meta %>% sample_data())

physeq_read_count
```

```
## phyloseq-class experiment-level object
## otu_table()   OTU Table:         [ 1634985 taxa and 49 samples ]
## sample_data() Sample Data:       [ 49 samples by 63 sample variables ]
## tax_table()   Taxonomy Table:    [ 1634985 taxa by 99 taxonomic ranks ]
```



```r
phyloseq_list_out = list("read_count" = physeq_read_count,
                         "Num_Gi" = physeq_Num_Gi,
                         "Num_Gi_pc" = physeq_Num_Gi_pc,
                         "cov" = physeq_cov,
                         "tpm" = physeq_tpm)
```


Export outputs:


```r
phyloseq_list_out %>% 
  saveRDS(here::here("data/processed/chicken1_full_gene_catalog_phyloseq.RDS"))

full_gene_catalogue_orf_quant_all_metrics %>% 
  saveRDS(here::here("data/processed/chicken1_full_gene_catalog_full_metrics_table.RDS"))

full_gene_catalogue_orf_quant_all_metrics %>% 
  write_tsv(here::here("data/processed/chicken1_full_gene_catalog_full_metrics_table.tsv.gz"))
```


Export only AMR and Tox:


```r
full_gene_catalogue_orf_quant_all_metrics %>% 
  filter(!is.na(PathoFact_AMR_ARG) | !is.na(rgi_CARD_ARO) | !is.na(PathoFact_Tox_Toxin_classification)) %>% 
  write_tsv(here::here("data/processed/chicken1_full_gene_catalog_full_metrics_table_AMR_tox.tsv"))
```




```r
devtools::session_info()
```

```
## ─ Session info ───────────────────────────────────────────────────────────────
##  setting  value
##  version  R version 4.2.0 (2022-04-22)
##  os       macOS Mojave 10.14.6
##  system   x86_64, darwin17.0
##  ui       X11
##  language (EN)
##  collate  en_US.UTF-8
##  ctype    en_US.UTF-8
##  tz       Europe/Zurich
##  date     2022-07-04
##  pandoc   2.14.0.3 @ /Applications/RStudio.app/Contents/MacOS/pandoc/ (via rmarkdown)
## 
## ─ Packages ───────────────────────────────────────────────────────────────────
##  package          * version    date (UTC) lib source
##  ade4               1.7-19     2022-04-19 [1] CRAN (R 4.2.0)
##  ape                5.6-2      2022-03-02 [1] CRAN (R 4.2.0)
##  assertthat         0.2.1      2019-03-21 [1] CRAN (R 4.2.0)
##  backports          1.4.1      2021-12-13 [1] CRAN (R 4.2.0)
##  Biobase            2.56.0     2022-04-26 [1] Bioconductor
##  BiocGenerics       0.42.0     2022-04-26 [1] Bioconductor
##  biomformat         1.24.0     2022-04-26 [1] Bioconductor
##  Biostrings         2.64.0     2022-04-26 [1] Bioconductor
##  bit                4.0.4      2020-08-04 [1] CRAN (R 4.2.0)
##  bit64              4.0.5      2020-08-30 [1] CRAN (R 4.2.0)
##  bitops             1.0-7      2021-04-24 [1] CRAN (R 4.2.0)
##  brio               1.1.3      2021-11-30 [1] CRAN (R 4.2.0)
##  broom              0.8.0      2022-04-13 [1] CRAN (R 4.2.0)
##  bslib              0.3.1      2021-10-06 [1] CRAN (R 4.2.0)
##  cachem             1.0.6      2021-08-19 [1] CRAN (R 4.2.0)
##  callr              3.7.0      2021-04-20 [1] CRAN (R 4.2.0)
##  cellranger         1.1.0      2016-07-27 [1] CRAN (R 4.2.0)
##  cli                3.3.0      2022-04-25 [1] CRAN (R 4.2.0)
##  cluster            2.1.3      2022-03-28 [1] CRAN (R 4.2.0)
##  codetools          0.2-18     2020-11-04 [1] CRAN (R 4.2.0)
##  colorspace         2.0-3      2022-02-21 [1] CRAN (R 4.2.0)
##  crayon             1.5.1      2022-03-26 [1] CRAN (R 4.2.0)
##  crosstalk          1.2.0      2021-11-04 [1] CRAN (R 4.2.0)
##  data.table         1.14.2     2021-09-27 [1] CRAN (R 4.2.0)
##  DBI                1.1.2      2021-12-20 [1] CRAN (R 4.2.0)
##  dbplyr             2.2.0      2022-06-05 [1] CRAN (R 4.2.0)
##  desc               1.4.1      2022-03-06 [1] CRAN (R 4.2.0)
##  devtools           2.4.3      2021-11-30 [1] CRAN (R 4.2.0)
##  digest             0.6.29     2021-12-01 [1] CRAN (R 4.2.0)
##  dplyr            * 1.0.9      2022-04-28 [1] CRAN (R 4.2.0)
##  DT                 0.23       2022-05-10 [1] CRAN (R 4.2.0)
##  dtplyr             1.2.1      2022-01-19 [1] CRAN (R 4.2.0)
##  ellipsis           0.3.2      2021-04-29 [1] CRAN (R 4.2.0)
##  evaluate           0.15       2022-02-18 [1] CRAN (R 4.2.0)
##  fansi              1.0.3      2022-03-24 [1] CRAN (R 4.2.0)
##  fastmap            1.1.0      2021-01-25 [1] CRAN (R 4.2.0)
##  forcats          * 0.5.1      2021-01-27 [1] CRAN (R 4.2.0)
##  foreach            1.5.2      2022-02-02 [1] CRAN (R 4.2.0)
##  fs                 1.5.2      2021-12-08 [1] CRAN (R 4.2.0)
##  gargle             1.2.0      2021-07-02 [1] CRAN (R 4.2.0)
##  generics           0.1.2      2022-01-31 [1] CRAN (R 4.2.0)
##  GenomeInfoDb       1.32.2     2022-05-15 [1] Bioconductor
##  GenomeInfoDbData   1.2.8      2022-06-13 [1] Bioconductor
##  ggplot2          * 3.3.6      2022-05-03 [1] CRAN (R 4.2.0)
##  glue               1.6.2      2022-02-24 [1] CRAN (R 4.2.0)
##  googledrive        2.0.0      2021-07-08 [1] CRAN (R 4.2.0)
##  googlesheets4      1.0.0      2021-07-21 [1] CRAN (R 4.2.0)
##  gtable             0.3.0      2019-03-25 [1] CRAN (R 4.2.0)
##  haven              2.5.0      2022-04-15 [1] CRAN (R 4.2.0)
##  here             * 1.0.1      2020-12-13 [1] CRAN (R 4.2.0)
##  hms                1.1.1      2021-09-26 [1] CRAN (R 4.2.0)
##  htmltools          0.5.2      2021-08-25 [1] CRAN (R 4.2.0)
##  htmlwidgets        1.5.4      2021-09-08 [1] CRAN (R 4.2.0)
##  httr               1.4.3      2022-05-04 [1] CRAN (R 4.2.0)
##  igraph             1.3.1      2022-04-20 [1] CRAN (R 4.2.0)
##  IRanges            2.30.0     2022-04-26 [1] Bioconductor
##  iterators          1.0.14     2022-02-05 [1] CRAN (R 4.2.0)
##  jquerylib          0.1.4      2021-04-26 [1] CRAN (R 4.2.0)
##  jsonlite           1.8.0      2022-02-22 [1] CRAN (R 4.2.0)
##  knitr              1.39       2022-04-26 [1] CRAN (R 4.2.0)
##  lattice            0.20-45    2021-09-22 [1] CRAN (R 4.2.0)
##  lifecycle          1.0.1      2021-09-24 [1] CRAN (R 4.2.0)
##  lubridate          1.8.0      2021-10-07 [1] CRAN (R 4.2.0)
##  magrittr           2.0.3      2022-03-30 [1] CRAN (R 4.2.0)
##  MASS               7.3-57     2022-04-22 [1] CRAN (R 4.2.0)
##  Matrix             1.4-1      2022-03-23 [1] CRAN (R 4.2.0)
##  memoise            2.0.1      2021-11-26 [1] CRAN (R 4.2.0)
##  mgcv               1.8-40     2022-03-29 [1] CRAN (R 4.2.0)
##  modelr             0.1.8      2020-05-19 [1] CRAN (R 4.2.0)
##  multtest           2.52.0     2022-04-26 [1] Bioconductor
##  munsell            0.5.0      2018-06-12 [1] CRAN (R 4.2.0)
##  nlme               3.1-157    2022-03-25 [1] CRAN (R 4.2.0)
##  permute            0.9-7      2022-01-27 [1] CRAN (R 4.2.0)
##  phyloseq         * 1.40.0     2022-04-26 [1] Bioconductor
##  pillar             1.7.0      2022-02-01 [1] CRAN (R 4.2.0)
##  pkgbuild           1.3.1      2021-12-20 [1] CRAN (R 4.2.0)
##  pkgconfig          2.0.3      2019-09-22 [1] CRAN (R 4.2.0)
##  pkgload            1.2.4      2021-11-30 [1] CRAN (R 4.2.0)
##  plyr               1.8.7      2022-03-24 [1] CRAN (R 4.2.0)
##  prettyunits        1.1.1      2020-01-24 [1] CRAN (R 4.2.0)
##  processx           3.6.0      2022-06-10 [1] CRAN (R 4.2.0)
##  ps                 1.7.0      2022-04-23 [1] CRAN (R 4.2.0)
##  purrr            * 0.3.4      2020-04-17 [1] CRAN (R 4.2.0)
##  R6                 2.5.1      2021-08-19 [1] CRAN (R 4.2.0)
##  Rcpp               1.0.8.3    2022-03-17 [1] CRAN (R 4.2.0)
##  RCurl              1.98-1.7   2022-06-09 [1] CRAN (R 4.2.0)
##  readr            * 2.1.2      2022-01-30 [1] CRAN (R 4.2.0)
##  readxl           * 1.4.0      2022-03-28 [1] CRAN (R 4.2.0)
##  remotes            2.4.2      2021-11-30 [1] CRAN (R 4.2.0)
##  reprex             2.0.1      2021-08-05 [1] CRAN (R 4.2.0)
##  reshape2           1.4.4      2020-04-09 [1] CRAN (R 4.2.0)
##  rhdf5              2.40.0     2022-04-26 [1] Bioconductor
##  rhdf5filters       1.8.0      2022-04-26 [1] Bioconductor
##  Rhdf5lib           1.18.2     2022-05-17 [1] Bioconductor
##  rlang              1.0.2      2022-03-04 [1] CRAN (R 4.2.0)
##  rmarkdown          2.14       2022-04-25 [1] CRAN (R 4.2.0)
##  rprojroot          2.0.3      2022-04-02 [1] CRAN (R 4.2.0)
##  rstudioapi         0.13       2020-11-12 [1] CRAN (R 4.2.0)
##  rvest              1.0.2      2021-10-16 [1] CRAN (R 4.2.0)
##  S4Vectors          0.34.0     2022-04-26 [1] Bioconductor
##  sass               0.4.1      2022-03-23 [1] CRAN (R 4.2.0)
##  scales             1.2.0      2022-04-13 [1] CRAN (R 4.2.0)
##  sessioninfo        1.2.2      2021-12-06 [1] CRAN (R 4.2.0)
##  stringi            1.7.6      2021-11-29 [1] CRAN (R 4.2.0)
##  stringr          * 1.4.0      2019-02-10 [1] CRAN (R 4.2.0)
##  survival           3.3-1      2022-03-03 [1] CRAN (R 4.2.0)
##  testthat           3.1.4      2022-04-26 [1] CRAN (R 4.2.0)
##  tibble           * 3.1.7      2022-05-03 [1] CRAN (R 4.2.0)
##  tidyr            * 1.2.0      2022-02-01 [1] CRAN (R 4.2.0)
##  tidyselect         1.1.2      2022-02-21 [1] CRAN (R 4.2.0)
##  tidyverse        * 1.3.1.9000 2022-06-13 [1] Github (tidyverse/tidyverse@6186fbf)
##  tzdb               0.3.0      2022-03-28 [1] CRAN (R 4.2.0)
##  usethis            2.1.6      2022-05-25 [1] CRAN (R 4.2.0)
##  utf8               1.2.2      2021-07-24 [1] CRAN (R 4.2.0)
##  vctrs              0.4.1      2022-04-13 [1] CRAN (R 4.2.0)
##  vegan              2.6-2      2022-04-17 [1] CRAN (R 4.2.0)
##  vroom              1.5.7      2021-11-30 [1] CRAN (R 4.2.0)
##  withr              2.5.0      2022-03-03 [1] CRAN (R 4.2.0)
##  xfun               0.31       2022-05-10 [1] CRAN (R 4.2.0)
##  xml2               1.3.3      2021-11-30 [1] CRAN (R 4.2.0)
##  XVector            0.36.0     2022-04-26 [1] Bioconductor
##  yaml               2.3.5      2022-02-21 [1] CRAN (R 4.2.0)
##  zlibbioc           1.42.0     2022-04-26 [1] Bioconductor
## 
##  [1] /Library/Frameworks/R.framework/Versions/4.2/Resources/library
## 
## ──────────────────────────────────────────────────────────────────────────────
```

