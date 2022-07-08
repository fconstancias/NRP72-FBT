---
title: "Merge human gene catalogue + quantification "
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
# cd /Users/fconstan/Documents/GitHub/NRP72-FBT/data/processed/anvio/human1/
# # anvi-get-split-coverages -p PROFILE.db  -c CONTIGS.db --list-splits  > all_splits.tsv
# "/Users/fconstan/Documents/GitHub/NRP72-FBT/data/processed/anvio/human1/all_splits.tsv" %>%
#   read_tsv(col_names = c("split", "bin")) %>%
#   mutate(bin = "all") %>%
#   data.frame() ->  Manual_bins_tmp
# # # 
# # # 
# Manual_bins_tmp %>%
#   write_tsv(here::here("data/processed/anvio/human1/collection_all_splits.tsv"), col_names = FALSE)
# 
# 
# anvi-import-collection -p PROFILE.db  -c CONTIGS.db  -C all_splits collection_all_splits.tsv
# anvi-summarize  -p PROFILE.db  -c CONTIGS.db  --list-collections
# anvi-summarize  -p PROFILE.db  -c CONTIGS.db -C all_splits -o summary_no_bins 

#anvi-summarize  -p PROFILE.db  -c CONTIGS.db  -C concoct_10_raw -o summary_no_bins # **actually this is al the splits : NOT SURE**
```

* FYI: A subset of split sequences are being initialized (to be precise, only
57125 of 261459 splits the contigs database knows about). Nothing to worry
about. Probably.



# Define inputs:




## tables/ data:


```r
SQMR <- "data/raw/SQM/human1_SQM.rds"

metadata <- "data/raw/25.11.2021_metadata_updated.tsv"


RGI_CARD <- "data/raw/CARD/rgi_contig_human1.txt"
Viralverify <- "data/raw/mobilome/viralverify/01.human1_result_table.csv"
isescan <- "data/raw/mobilome/isescan/01.SqueezeHuman1.fasta.tsv"
PathoFact_AMR <- "data/processed/PathoFact/humann1/AMR/AMR_MGE_prediction_sample_l1000_report.tsv"
PathoFact_Tox <- "data/processed/PathoFact/humann1/Tox/Toxin_prediction_sample_l1000_report.tsv"
PathoFact_Tox_fam <- "data/processed/PathoFact/humann1/Tox/Toxin_gene_library_sample_l1000_report.tsv"


anvio_nobin_summary <- "data/processed/anvio/human1/summary_no_bins/bin_by_bin/all/all-gene_calls.txt"

DAS_collection <- "data/processed/anvio/human1/DAS.txt"
Manual_bins <- "data/processed/anvio/human1/MAGS/concoct_10_collection_refined_2.txt"

ANVIO_SUMMARY_DAS  <- "data/processed/anvio/human1/MAGS/SUMMARY/DAS/bins_summary.txt"
ANVIO_SUMMARY_CONCOCT <- "data/processed/anvio/human1/MAGS/SUMMARY/concoct_10_refined_2/bins_summary.txt"

GTDB_tax <- "data/processed/gtdb/GTDB_all.xlsx"
SQM_DAS <- "data/raw/SQM/19.human1.bintable"

# contigs.lgt.tsv <- "~/Documents/MSc_Thesis-NRP72/data/HGT/human.contigs.lgt.tsv"
# contigs.no_lgt.tsv <- "~/Documents/MSc_Thesis-NRP72/data/HGT/human.contigs.no_lgt.tsv"
# contigs.unclassified.tsv <- "~/Documents/MSc_Thesis-NRP72/data/HGT/human.contigs.unclassified.tsv"
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
  filter(!is.na(metagenomic_sample_name), Model == "Human") -> human_meta
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
## Rows: 563766 Columns: 11
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
## Rows: 563766 Columns: 6
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
## [1] 8834    2
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
## Rows: 10189 Columns: 9
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
## [1] 8834    6
```


```r
PathoFact_Tox %>% 
  left_join(PathoFact_Tox_fam,
            by = c("ORF" = "ORF")) -> PathoFact_Tox


PathoFact_Tox %>% 
  dim()
```

```
## [1] 8834    7
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
<div id="htmlwidget-829c3feb025761a4fb2c" style="width:100%;height:auto;" class="datatables html-widget"></div>
<script type="application/json" data-for="htmlwidget-829c3feb025761a4fb2c">{"x":{"filter":"none","vertical":false,"data":[["1","2","3","4","5","6"],["megahit_117_1617-2000","megahit_137_1282-2034","megahit_170_64-381","megahit_364_413-814","megahit_374_4025-4759","megahit_439_1295-2875"],["2: Non-secreted Toxin","2: Non-secreted Toxin","2: Non-secreted Toxin","1: Secreted Toxin","2: Non-secreted Toxin","1: Secreted Toxin"],["PF00903.20","K11068_1278_95_625","PF01325.14","PF07883.6","PF13950.1","K10815_8"],[62.7,190.2,40.5,53,101.9,130.4],[2.1e-17,1.7e-56,1.2e-10,1.1e-14,8.2e-30,3.3e-38],["Glyoxalase","hlyIII","Fe_dep_repress","Cupin_2","Epimerase_Csub","hcnB"],["Glyoxalase/Bleomycin resistance protein/Dioxygenase superfamily","hemolysin III","Iron dependent repressor, N-terminal DNA binding domain","Cupin domain","UDP-glucose 4-epimerase C-term subunit","hydrogen cyanide synthase HcnB [EC:1.4.99.5]"]],"container":"<table class=\"display\">\n  <thead>\n    <tr>\n      <th> <\/th>\n      <th>ORF_ID<\/th>\n      <th>PathoFact_Tox_Toxin_classification<\/th>\n      <th>PathoFact_Tox_HMM_Name<\/th>\n      <th>PathoFact_Tox_Score<\/th>\n      <th>PathoFact_Tox_Significance_evalue<\/th>\n      <th>PathoFact_Tox_NAME<\/th>\n      <th>PathoFact_Tox_Description<\/th>\n    <\/tr>\n  <\/thead>\n<\/table>","options":{"columnDefs":[{"className":"dt-right","targets":[4,5]},{"orderable":false,"targets":0}],"order":[],"autoWidth":false,"orderClasses":false}},"evals":[],"jsHooks":[]}</script>
```
## rgi CARD:


```r
# Resistome data:
RGI_CARD %>% 
  here::here() %>% 
  read_tsv() %>% 
  # read.delim() %>%
  select(-ORF_ID) %>%
  dplyr::filter(Model_type %in% Model_type ,
                Cut_Off %in% Cut_Off,
                Nudged %!in% c("TRUE"),
                Best_Identities > besid, # plot distribution of bestID vs percentage length ref and color bitscore and symbol class AMR
                `Percentage Length of Reference Sequence` > Percent_length) %>%
  separate(Contig, c("Megahit", "Contig", "Split"), sep = "_", fill = "right") %>%
  mutate(ORF_ID = paste(Megahit, Contig, paste(Start, Stop, sep = "-"), sep = "_")) %>%
  #mutate(Contig = paste(Megahit, Contig, sep = "_")) %>%
  select(-c(Megahit, Split, CARD_Protein_Sequence, Predicted_DNA, Predicted_Protein, Contig)) %>% 
  na_if("NA") %>% 
  na_if("n/a") %>% 
  na_if("<NA>") %>% 
  na_if("") %>% 
  select(-Start, -Stop, -Orientation) %>% 
  select(ORF_ID, everything()) -> card
```

```
## Rows: 596 Columns: 25
## ── Column specification ────────────────────────────────────────────────────────
## Delimiter: "\t"
## chr (16): ORF_ID, Contig, Orientation, Cut_Off, Best_Hit_ARO, Model_type, SN...
## dbl  (8): Start, Stop, Pass_Bitscore, Best_Hit_Bitscore, Best_Identities, AR...
## lgl  (1): Nudged
## 
## ℹ Use `spec()` to retrieve the full column specification for this data.
## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.
```

```r
colnames(card) <- paste0("rgi_CARD_", colnames(card))

card %>% 
  rename(ORF_ID = rgi_CARD_ORF_ID) -> card

card %>% 
  head() %>% 
  DT::datatable()
```

```{=html}
<div id="htmlwidget-2bcea764f9df3988b898" style="width:100%;height:auto;" class="datatables html-widget"></div>
<script type="application/json" data-for="htmlwidget-2bcea764f9df3988b898">{"x":{"filter":"none","vertical":false,"data":[["1","2","3","4","5","6"],["megahit_1512_86-616","megahit_1512_743-1915","megahit_1512_1932-3470","megahit_1935_13232-13930","megahit_4077_533-1408","megahit_4631_5392-6930"],["Perfect","Strict","Perfect","Strict","Perfect","Strict"],[280,675,900,400,500,900],[349,739.6,996.1,458,558.1,963.4],["emrR","emrA","emrB","vanRD","CTX-M-1","emrB"],[100,99.74,100,96.12,100,95.51],[3000516,3000027,3000074,3002923,3001864,3000074],["protein homolog model","protein homolog model","protein homolog model","protein homolog model","protein homolog model","protein homolog model"],[null,null,null,null,null,null],[null,null,null,null,null,null],["fluoroquinolone antibiotic","fluoroquinolone antibiotic","fluoroquinolone antibiotic","glycopeptide antibiotic","cephalosporin","fluoroquinolone antibiotic"],["antibiotic efflux","antibiotic efflux","antibiotic efflux","antibiotic target alteration","antibiotic inactivation","antibiotic efflux"],["major facilitator superfamily (MFS) antibiotic efflux pump","major facilitator superfamily (MFS) antibiotic efflux pump","major facilitator superfamily (MFS) antibiotic efflux pump","glycopeptide resistance gene cluster; vanR","CTX-M beta-lactamase","major facilitator superfamily (MFS) antibiotic efflux pump"],[100,100,100,100,100,100],["gnl|BL_ORD_ID|1254|hsp_num:0","gnl|BL_ORD_ID|1657|hsp_num:0","gnl|BL_ORD_ID|1746|hsp_num:0","gnl|BL_ORD_ID|950|hsp_num:0","gnl|BL_ORD_ID|392|hsp_num:0","gnl|BL_ORD_ID|1746|hsp_num:0"],[1330,1757,1847,1004,420,1847],[null,null,null,null,null,null],[null,null,null,null,null,null]],"container":"<table class=\"display\">\n  <thead>\n    <tr>\n      <th> <\/th>\n      <th>ORF_ID<\/th>\n      <th>rgi_CARD_Cut_Off<\/th>\n      <th>rgi_CARD_Pass_Bitscore<\/th>\n      <th>rgi_CARD_Best_Hit_Bitscore<\/th>\n      <th>rgi_CARD_Best_Hit_ARO<\/th>\n      <th>rgi_CARD_Best_Identities<\/th>\n      <th>rgi_CARD_ARO<\/th>\n      <th>rgi_CARD_Model_type<\/th>\n      <th>rgi_CARD_SNPs_in_Best_Hit_ARO<\/th>\n      <th>rgi_CARD_Other_SNPs<\/th>\n      <th>rgi_CARD_Drug Class<\/th>\n      <th>rgi_CARD_Resistance Mechanism<\/th>\n      <th>rgi_CARD_AMR Gene Family<\/th>\n      <th>rgi_CARD_Percentage Length of Reference Sequence<\/th>\n      <th>rgi_CARD_ID<\/th>\n      <th>rgi_CARD_Model_ID<\/th>\n      <th>rgi_CARD_Nudged<\/th>\n      <th>rgi_CARD_Note<\/th>\n    <\/tr>\n  <\/thead>\n<\/table>","options":{"columnDefs":[{"className":"dt-right","targets":[3,4,6,7,14,16]},{"orderable":false,"targets":0}],"order":[],"autoWidth":false,"orderClasses":false}},"evals":[],"jsHooks":[]}</script>
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
## Rows: 253075 Columns: 6
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
## Rows: 6291 Columns: 24
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
<div id="htmlwidget-ab17de5797c309e0abf9" style="width:100%;height:auto;" class="datatables html-widget"></div>
<script type="application/json" data-for="htmlwidget-ab17de5797c309e0abf9">{"x":{"filter":"none","vertical":false,"data":[["1","2","3","4","5","6"],["megahit_100032","megahit_10009","megahit_100143","megahit_100161","megahit_100166","megahit_100186"],["IS30","ISNCY","IS66","IS200/IS605","IS30","IS66"],["IS30_241","ISNCY_86","IS66_46","IS200/IS605_96","IS30_223","IS66_43"],["p","p","p","p","c","p"],[1,1,1,1,1,2]],"container":"<table class=\"display\">\n  <thead>\n    <tr>\n      <th> <\/th>\n      <th>Contig_ID<\/th>\n      <th>isescan_family<\/th>\n      <th>isescan_cluster<\/th>\n      <th>isescan_type<\/th>\n      <th>isescan_count_IS_contig<\/th>\n    <\/tr>\n  <\/thead>\n<\/table>","options":{"columnDefs":[{"className":"dt-right","targets":5},{"orderable":false,"targets":0}],"order":[],"autoWidth":false,"orderClasses":false}},"evals":[],"jsHooks":[]}</script>
```
## plasx:



```r
# plasx %>% 
#   here::here() %>% 
#   read_tsv(col_names = c("Contig_ID", "plasX_score")) -> pasx_score
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
select(-Method, -Molecule) %>% 
  rename(Contig_ID = `Contig.ID`,
         Tax_orf = Tax) %>% 
  select(-TPM.D.1:-Hits) -> orf_tmp


orf_tmp %>% 
  select(-(starts_with(c("TPM.","Raw.base.count.","Coverage.", "Raw.read.count")))) -> orf
# starts_with(c("Petal", "Sepal"))
# orf %>% 
#   write_tsv("~/Desktop/unitedhuman1.tsv")
# orf %>% 
#   write_tsv("~/Desktop/unitedhuman1.tsv")
colnames(orf) <- paste0("SQM_", colnames(orf))

orf %>% 
  rename(Contig_ID = SQM_Contig_ID,
         ORF_ID = SQM_ORF_ID) -> orf

orf %>% 
  head() %>% 
  DT::datatable()
```

```{=html}
<div id="htmlwidget-c98ac71b298e8735f2c5" style="width:100%;height:auto;" class="datatables html-widget"></div>
<script type="application/json" data-for="htmlwidget-c98ac71b298e8735f2c5">{"x":{"filter":"none","vertical":false,"data":[["1","2","3","4","5","6"],["megahit_1_1-312","megahit_2_3-176","megahit_2_173-304","megahit_3_3-116","megahit_3_82-219","megahit_4_1-72"],["megahit_1","megahit_2","megahit_2","megahit_3","megahit_3","megahit_4"],[312,174,132,114,138,72],[104,58,44,38,46,24],[34.62,38.51,40.91,51.75,47.83,36.11],["RP-S5, MRPS5, rpsE",null,null,null,null,null],["k_Bacteria;n_Terrabacteria group;p_Firmicutes;c_Erysipelotrichia;o_Erysipelotrichales;f_Erysipelotrichaceae;g_Massilimicrobiota","k_Bacteria;n_Terrabacteria group;p_Firmicutes;c_Clostridia;o_Eubacteriales;n_unclassified Eubacteriales;s_Clostridiales bacterium","k_Bacteria;n_Terrabacteria group;p_Firmicutes","k_Bacteria;n_Terrabacteria group;p_Firmicutes","k_Bacteria;n_Terrabacteria group;p_Firmicutes;c_Clostridia;o_Eubacteriales","k_Bacteria;n_Terrabacteria group;p_Firmicutes;c_Bacilli;o_Lactobacillales;f_Enterococcaceae;g_Enterococcus;s_Enterococcus faecium"],["K02988*",null,null,null,null,null],["small subunit ribosomal protein S5",null,null,null,null,null],["Genetic Information Processing; Translation; Ribosome | Brite Hierarchies; Protein families: genetic information processing; Ribosome",null,null,null,null,null],["COG0098*",null,null,"ENOG410ZP7K*","COG3843*",null],["Ribosomal protein S5",null,null,"D repair protein RadA domain","Type IV secretory pathway, VirD2 components (relaxase)",null],["Translation, ribosomal structure and biogenesis",null,null,"Function unknown",null,null],["PF03719 [Ribosomal protein S5, C-terminal domain]",null,null,null,null,null]],"container":"<table class=\"display\">\n  <thead>\n    <tr>\n      <th> <\/th>\n      <th>ORF_ID<\/th>\n      <th>Contig_ID<\/th>\n      <th>SQM_Length.NT<\/th>\n      <th>SQM_Length.AA<\/th>\n      <th>SQM_GC.perc<\/th>\n      <th>SQM_Gene.name<\/th>\n      <th>SQM_Tax_orf<\/th>\n      <th>SQM_KEGG.ID<\/th>\n      <th>SQM_KEGGFUN<\/th>\n      <th>SQM_KEGGPATH<\/th>\n      <th>SQM_COG.ID<\/th>\n      <th>SQM_COGFUN<\/th>\n      <th>SQM_COGPATH<\/th>\n      <th>SQM_PFAM<\/th>\n    <\/tr>\n  <\/thead>\n<\/table>","options":{"columnDefs":[{"className":"dt-right","targets":[3,4,5]},{"orderable":false,"targets":0}],"order":[],"autoWidth":false,"orderClasses":false}},"evals":[],"jsHooks":[]}</script>
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
<div id="htmlwidget-69ac827b3ad27432a89c" style="width:100%;height:auto;" class="datatables html-widget"></div>
<script type="application/json" data-for="htmlwidget-69ac827b3ad27432a89c">{"x":{"filter":"none","vertical":false,"data":[["1","2","3","4","5","6"],["megahit_1","megahit_2","megahit_3","megahit_4","megahit_5","megahit_6"],["k_Bacteria;n_Terrabacteria group;p_Firmicutes;c_Erysipelotrichia;o_Erysipelotrichales;f_Erysipelotrichaceae;g_Massilimicrobiota","k_Bacteria;n_Terrabacteria group;p_Firmicutes;c_Clostridia;o_Eubacteriales;n_unclassified Eubacteriales;s_Clostridiales bacterium","k_Bacteria;n_Terrabacteria group;p_Firmicutes;c_Clostridia;o_Eubacteriales","k_Bacteria;n_Terrabacteria group;p_Firmicutes;c_Bacilli;o_Lactobacillales;f_Enterococcaceae;g_Enterococcus;s_Enterococcus faecium","","k_Bacteria;n_Terrabacteria group;p_Firmicutes;c_Clostridia;o_Eubacteriales"],[0,0,0,0,0,0],[34.5,39.67,51.13,38.89,50.74,45.22],[313,305,221,324,339,230],[1,2,2,2,1,1]],"container":"<table class=\"display\">\n  <thead>\n    <tr>\n      <th> <\/th>\n      <th>Contig_ID<\/th>\n      <th>SQM_Tax_contig<\/th>\n      <th>SQM_Disparity_contig<\/th>\n      <th>SQM_GC_perc_contig<\/th>\n      <th>SQM_Length_contig<\/th>\n      <th>SQM_Num_genes_contig<\/th>\n    <\/tr>\n  <\/thead>\n<\/table>","options":{"columnDefs":[{"className":"dt-right","targets":[3,4,5,6]},{"orderable":false,"targets":0}],"order":[],"autoWidth":false,"orderClasses":false}},"evals":[],"jsHooks":[]}</script>
```



```r
orf %>% 
  left_join(card,
            by = c("ORF_ID" = "ORF_ID")) %>% 
  left_join(PathoFact_AMR,
            by = c("ORF_ID" = "ORF_ID")) %>% 
  left_join(PathoFact_Tox,
            by = c("ORF_ID" = "ORF_ID")) %>% 
  left_join(viralverify,
            by = c("Contig_ID" = "Contig_ID")) %>% 
  left_join(top_isescan,
            by = c("Contig_ID" = "Contig_ID")) %>% 
  left_join(contigs,
            by = c("Contig_ID" = "Contig_ID")) -> orf_contig

orf_contig %>% 
  head() %>% 
  DT::datatable()
```

```{=html}
<div id="htmlwidget-45eaf4164604819f8b59" style="width:100%;height:auto;" class="datatables html-widget"></div>
<script type="application/json" data-for="htmlwidget-45eaf4164604819f8b59">{"x":{"filter":"none","vertical":false,"data":[["1","2","3","4","5","6"],["megahit_1_1-312","megahit_2_3-176","megahit_2_173-304","megahit_3_3-116","megahit_3_82-219","megahit_4_1-72"],["megahit_1","megahit_2","megahit_2","megahit_3","megahit_3","megahit_4"],[312,174,132,114,138,72],[104,58,44,38,46,24],[34.62,38.51,40.91,51.75,47.83,36.11],["RP-S5, MRPS5, rpsE",null,null,null,null,null],["k_Bacteria;n_Terrabacteria group;p_Firmicutes;c_Erysipelotrichia;o_Erysipelotrichales;f_Erysipelotrichaceae;g_Massilimicrobiota","k_Bacteria;n_Terrabacteria group;p_Firmicutes;c_Clostridia;o_Eubacteriales;n_unclassified Eubacteriales;s_Clostridiales bacterium","k_Bacteria;n_Terrabacteria group;p_Firmicutes","k_Bacteria;n_Terrabacteria group;p_Firmicutes","k_Bacteria;n_Terrabacteria group;p_Firmicutes;c_Clostridia;o_Eubacteriales","k_Bacteria;n_Terrabacteria group;p_Firmicutes;c_Bacilli;o_Lactobacillales;f_Enterococcaceae;g_Enterococcus;s_Enterococcus faecium"],["K02988*",null,null,null,null,null],["small subunit ribosomal protein S5",null,null,null,null,null],["Genetic Information Processing; Translation; Ribosome | Brite Hierarchies; Protein families: genetic information processing; Ribosome",null,null,null,null,null],["COG0098*",null,null,"ENOG410ZP7K*","COG3843*",null],["Ribosomal protein S5",null,null,"D repair protein RadA domain","Type IV secretory pathway, VirD2 components (relaxase)",null],["Translation, ribosomal structure and biogenesis",null,null,"Function unknown",null,null],["PF03719 [Ribosomal protein S5, C-terminal domain]",null,null,null,null,null],[null,null,null,null,null,null],[null,null,null,null,null,null],[null,null,null,null,null,null],[null,null,null,null,null,null],[null,null,null,null,null,null],[null,null,null,null,null,null],[null,null,null,null,null,null],[null,null,null,null,null,null],[null,null,null,null,null,null],[null,null,null,null,null,null],[null,null,null,null,null,null],[null,null,null,null,null,null],[null,null,null,null,null,null],[null,null,null,null,null,null],[null,null,null,null,null,null],[null,null,null,null,null,null],[null,null,null,null,null,null],[null,null,null,null,null,null],[null,null,null,null,null,null],[null,null,null,null,null,null],[null,null,null,null,null,null],[null,null,null,null,null,null],[null,null,null,null,null,null],[null,null,null,null,null,null],[null,null,null,null,null,null],[null,null,null,null,null,null],[null,null,null,null,null,null],[null,null,null,null,null,null],[null,null,null,null,null,null],[null,null,null,null,null,null],["Uncertain - too short","Uncertain - too short","Uncertain - too short","Uncertain - too short","Uncertain - too short","Uncertain - too short"],[null,null,null,null,null,null],[null,null,null,null,null,null],[null,null,null,null,null,null],[null,null,null,null,null,null],["k_Bacteria;n_Terrabacteria group;p_Firmicutes;c_Erysipelotrichia;o_Erysipelotrichales;f_Erysipelotrichaceae;g_Massilimicrobiota","k_Bacteria;n_Terrabacteria group;p_Firmicutes;c_Clostridia;o_Eubacteriales;n_unclassified Eubacteriales;s_Clostridiales bacterium","k_Bacteria;n_Terrabacteria group;p_Firmicutes;c_Clostridia;o_Eubacteriales;n_unclassified Eubacteriales;s_Clostridiales bacterium","k_Bacteria;n_Terrabacteria group;p_Firmicutes;c_Clostridia;o_Eubacteriales","k_Bacteria;n_Terrabacteria group;p_Firmicutes;c_Clostridia;o_Eubacteriales","k_Bacteria;n_Terrabacteria group;p_Firmicutes;c_Bacilli;o_Lactobacillales;f_Enterococcaceae;g_Enterococcus;s_Enterococcus faecium"],[0,0,0,0,0,0],[34.5,39.67,39.67,51.13,51.13,38.89],[313,305,305,221,221,324],[1,2,2,2,2,2]],"container":"<table class=\"display\">\n  <thead>\n    <tr>\n      <th> <\/th>\n      <th>ORF_ID<\/th>\n      <th>Contig_ID<\/th>\n      <th>SQM_Length.NT<\/th>\n      <th>SQM_Length.AA<\/th>\n      <th>SQM_GC.perc<\/th>\n      <th>SQM_Gene.name<\/th>\n      <th>SQM_Tax_orf<\/th>\n      <th>SQM_KEGG.ID<\/th>\n      <th>SQM_KEGGFUN<\/th>\n      <th>SQM_KEGGPATH<\/th>\n      <th>SQM_COG.ID<\/th>\n      <th>SQM_COGFUN<\/th>\n      <th>SQM_COGPATH<\/th>\n      <th>SQM_PFAM<\/th>\n      <th>rgi_CARD_Cut_Off<\/th>\n      <th>rgi_CARD_Pass_Bitscore<\/th>\n      <th>rgi_CARD_Best_Hit_Bitscore<\/th>\n      <th>rgi_CARD_Best_Hit_ARO<\/th>\n      <th>rgi_CARD_Best_Identities<\/th>\n      <th>rgi_CARD_ARO<\/th>\n      <th>rgi_CARD_Model_type<\/th>\n      <th>rgi_CARD_SNPs_in_Best_Hit_ARO<\/th>\n      <th>rgi_CARD_Other_SNPs<\/th>\n      <th>rgi_CARD_Drug Class<\/th>\n      <th>rgi_CARD_Resistance Mechanism<\/th>\n      <th>rgi_CARD_AMR Gene Family<\/th>\n      <th>rgi_CARD_Percentage Length of Reference Sequence<\/th>\n      <th>rgi_CARD_ID<\/th>\n      <th>rgi_CARD_Model_ID<\/th>\n      <th>rgi_CARD_Nudged<\/th>\n      <th>rgi_CARD_Note<\/th>\n      <th>PathoFact_AMR_ARG<\/th>\n      <th>PathoFact_AMR_ARG_SNPs<\/th>\n      <th>PathoFact_AMR_AMR_category<\/th>\n      <th>PathoFact_AMR_AMR_sub_class<\/th>\n      <th>PathoFact_AMR_Resistance_mechanism<\/th>\n      <th>PathoFact_AMR_Database<\/th>\n      <th>PathoFact_AMR_MGE_prediction<\/th>\n      <th>PathoFact_Tox_Toxin_classification<\/th>\n      <th>PathoFact_Tox_HMM_Name<\/th>\n      <th>PathoFact_Tox_Score<\/th>\n      <th>PathoFact_Tox_Significance_evalue<\/th>\n      <th>PathoFact_Tox_NAME<\/th>\n      <th>PathoFact_Tox_Description<\/th>\n      <th>viralverify_Prediction<\/th>\n      <th>isescan_family<\/th>\n      <th>isescan_cluster<\/th>\n      <th>isescan_type<\/th>\n      <th>isescan_count_IS_contig<\/th>\n      <th>SQM_Tax_contig<\/th>\n      <th>SQM_Disparity_contig<\/th>\n      <th>SQM_GC_perc_contig<\/th>\n      <th>SQM_Length_contig<\/th>\n      <th>SQM_Num_genes_contig<\/th>\n    <\/tr>\n  <\/thead>\n<\/table>","options":{"columnDefs":[{"className":"dt-right","targets":[3,4,5,16,17,19,20,27,29,41,42,49,51,52,53,54]},{"orderable":false,"targets":0}],"order":[],"autoWidth":false,"orderClasses":false}},"evals":[],"jsHooks":[]}</script>
```
## Test anvio no bin:


```r
anvio_nobin_summary %>% 
  here::here() %>%  
  read_tsv() -> anvi_no_bin_annot
```

```
## Rows: 527822 Columns: 22
## ── Column specification ────────────────────────────────────────────────────────
## Delimiter: "\t"
## chr (19): contig, direction, KEGG, KEGG (ACCESSION), COG, COG (ACCESSION), C...
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
# anvi_no_bin_annot$anvio_orf_id %>% 
#   head()
```


```r
# orf_contig %>%
#   filter(Contig_ID == "megahit_37") %>%
#   distinct(ORF_ID, .keep_all = TRUE)
```


```r
# orf_contig %>%
#   filter(ORF_ID == "megahit_37_124-318")
```

```r
# anvi_no_bin_annot %>% 
#     filter(anvio_contig == "megahit_37") %>%
#   distinct(anvio_orf_id, .keep_all = TRUE)
```


```r
# anvi_no_bin_annot %>% 
#     filter(anvio_orf_id == "megahit_37_124-318")
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
## Rows: 36166 Columns: 2
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
<div id="htmlwidget-bd786fc83ccbac70f1ca" style="width:100%;height:auto;" class="datatables html-widget"></div>
<script type="application/json" data-for="htmlwidget-bd786fc83ccbac70f1ca">{"x":{"filter":"none","vertical":false,"data":[["1","2","3","4","5","6"],["megahit_94475","megahit_135602","megahit_101245","megahit_222020","megahit_229957","megahit_94514"],["maxbin_154_fasta_contigs","maxbin_154_fasta_contigs","maxbin_154_fasta_contigs","maxbin_154_fasta_contigs","maxbin_154_fasta_contigs","maxbin_154_fasta_contigs"]],"container":"<table class=\"display\">\n  <thead>\n    <tr>\n      <th> <\/th>\n      <th>Contig_ID<\/th>\n      <th>DAS_bin<\/th>\n    <\/tr>\n  <\/thead>\n<\/table>","options":{"columnDefs":[{"orderable":false,"targets":0}],"order":[],"autoWidth":false,"orderClasses":false}},"evals":[],"jsHooks":[]}</script>
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
## Rows: 28212 Columns: 2
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
<div id="htmlwidget-f4add39831f97c659c2c" style="width:100%;height:auto;" class="datatables html-widget"></div>
<script type="application/json" data-for="htmlwidget-f4add39831f97c659c2c">{"x":{"filter":"none","vertical":false,"data":[["1","2","3","4","5","6"],["megahit_69635","megahit_218792","megahit_96692","megahit_139453","megahit_225023","megahit_105955"],["Bin_8_1","Bin_8_1","Bin_8_1","Bin_8_1","Bin_8_1","Bin_8_1"]],"container":"<table class=\"display\">\n  <thead>\n    <tr>\n      <th> <\/th>\n      <th>Contig_ID<\/th>\n      <th>CONCOCT_bin<\/th>\n    <\/tr>\n  <\/thead>\n<\/table>","options":{"columnDefs":[{"orderable":false,"targets":0}],"order":[],"autoWidth":false,"orderClasses":false}},"evals":[],"jsHooks":[]}</script>
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

### MAGs summary:

#### MAGs GDTB taxonomy:


```r
GTDB_tax %>%
  read_excel() %>%
  filter(dataset == "human1") %>%
  mutate(user_genome = str_remove_all(user_genome, c("human1_|-contigs|refined_2_"))) %>%
  filter(str_detect(user_genome, c("DAS|concoct"))) %>%
  select(user_genome, classification) %>%
  # separate(classification, 
  #          into = c("domain", "phylum", "class", "order", "family", "genus", "species"), 
  #          sep = ";", 
  #          fill = "right") %>% 
  rename(GTDB_tax = classification) %>% 
  mutate(user_genome = str_remove_all(user_genome, "concoct_10_")) %>% 
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
<div id="htmlwidget-657c5f90f5c3f82438ea" style="width:100%;height:auto;" class="datatables html-widget"></div>
<script type="application/json" data-for="htmlwidget-657c5f90f5c3f82438ea">{"x":{"filter":"none","vertical":false,"data":[["1","2","3","4","5","6"],["maxbin_002_fasta_contigs","maxbin_004_fasta_contigs","maxbin_005_fasta_contigs","maxbin_006_fasta_contigs","maxbin_008_fasta_contigs","maxbin_009_fasta_contigs"],["d__Bacteria;p__Bacteroidota;c__Bacteroidia;o__Bacteroidales;f__Bacteroidaceae;g__Phocaeicola;s__Phocaeicola plebeius_A","d__Bacteria;p__Firmicutes_A;c__Clostridia;o__Oscillospirales;f__Acutalibacteraceae;g__Ruminococcus_E;s__Ruminococcus_E bromii_B","d__Bacteria;p__Firmicutes_A;c__Clostridia;o__Oscillospirales;f__Ruminococcaceae;g__Faecalibacterium;s__Faecalibacterium prausnitzii_C","d__Bacteria;p__Firmicutes_A;c__Clostridia;o__Lachnospirales;f__Lachnospiraceae;g__Agathobacter;s__Agathobacter rectalis","d__Bacteria;p__Firmicutes_A;c__Clostridia;o__Lachnospirales;f__Lachnospiraceae;g__RUG115;s__RUG115 sp900066395","d__Bacteria;p__Firmicutes_A;c__Clostridia;o__Lachnospirales;f__Lachnospiraceae;g__Anaerostipes;s__Anaerostipes hadrus"]],"container":"<table class=\"display\">\n  <thead>\n    <tr>\n      <th> <\/th>\n      <th>user_genome<\/th>\n      <th>GTDB_tax<\/th>\n    <\/tr>\n  <\/thead>\n<\/table>","options":{"columnDefs":[{"orderable":false,"targets":0}],"order":[],"autoWidth":false,"orderClasses":false}},"evals":[],"jsHooks":[]}</script>
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
## Rows: 136 Columns: 109
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
  select(-SQM_DAS_Method, -`SQM_DAS_Coverage D-1`:-`SQM_DAS_TPM HV-24`) -> SQM_DAS_HQ

# filter(SQM_DAS_Completeness >= 80,
#        SQM_DAS_Contamination <= 10) %>% 
# mutate(SQM_DAS_bins_HQ = `SQM_DAS_Bin ID`) -> SQM_DAS_HQ

SQM_DAS_HQ %>% 
  head() %>% 
  DT::datatable()
```

```{=html}
<div id="htmlwidget-082ad8b64b19ec83c166" style="width:100%;height:auto;" class="datatables html-widget"></div>
<script type="application/json" data-for="htmlwidget-082ad8b64b19ec83c166">{"x":{"filter":"none","vertical":false,"data":[["1","2","3","4","5","6"],["maxbin.163.fasta.contigs","maxbin.043.fasta.contigs","maxbin.075.fasta.contigs","maxbin.153.fasta.contigs","maxbin.094.fasta.contigs","maxbin.171.fasta.contigs"],["k_Bacteria;n_Terrabacteria group;p_Firmicutes;c_Clostridia;o_Eubacteriales;f_Peptostreptococcaceae;g_Clostridioides;s_Clostridioides difficile","k_Bacteria;n_Terrabacteria group;p_Firmicutes;c_Bacilli;o_Lactobacillales;f_Enterococcaceae;g_Enterococcus","k_Bacteria;n_Terrabacteria group;p_Firmicutes;c_Erysipelotrichia;o_Erysipelotrichales","k_Bacteria;n_Terrabacteria group;p_Firmicutes;c_Bacilli;o_Lactobacillales;f_Enterococcaceae;g_Enterococcus","k_Bacteria;n_Terrabacteria group;p_Firmicutes;c_Clostridia;o_Eubacteriales","k_Bacteria;n_Terrabacteria group;p_Actinobacteria;c_Coriobacteriia;o_Coriobacteriales"],[null,null,null,"superkingdom:Bacteria;clade:Terrabacteria group;phylum:Firmicutes;class:Bacilli;order:Lactobacillales;family:Enterococcaceae;genus:Enterococcus",null,"superkingdom:Bacteria;clade:Terrabacteria group;phylum:Actinobacteria;class:Coriobacteriia;order:Coriobacteriales;family:Coriobacteriaceae"],[4025395,3586160,3293471,2805040,3182586,1869581],[28.35,42.33,31.31,37.88,55.16,49.32],[47,60,51,119,162,38],[0.094,0,0.08,0,0.013,0.066],[100,99.61,99.37,99.19,99.09,98.71],[0.33,0.62,0,0.78,12.62,0],[0,0,0,81.25,7.69,0],["HQ","HQ","HQ","HQ",null,"HQ"]],"container":"<table class=\"display\">\n  <thead>\n    <tr>\n      <th> <\/th>\n      <th>SQM_DAS_Bin ID<\/th>\n      <th>SQM_DAS_Tax<\/th>\n      <th>SQM_DAS_Tax 16S<\/th>\n      <th>SQM_DAS_Length<\/th>\n      <th>SQM_DAS_GC perc<\/th>\n      <th>SQM_DAS_Num contigs<\/th>\n      <th>SQM_DAS_Disparity<\/th>\n      <th>SQM_DAS_Completeness<\/th>\n      <th>SQM_DAS_Contamination<\/th>\n      <th>SQM_DAS_Strain het<\/th>\n      <th>SQM_DAS_HQ<\/th>\n    <\/tr>\n  <\/thead>\n<\/table>","options":{"columnDefs":[{"className":"dt-right","targets":[4,5,6,7,8,9,10]},{"orderable":false,"targets":0}],"order":[],"autoWidth":false,"orderClasses":false}},"evals":[],"jsHooks":[]}</script>
```


```r
orf_contig_3 %>% 
  left_join(SQM_DAS_HQ,
            by = c("DAS_bin" = "SQM_DAS_Bin ID")) -> orf_contig_3
```

#### ANVIO DAS info:


```r
ANVIO_SUMMARY_DAS %>% 
  here::here() %>% 
  read_tsv() %>% 
  select(-t_domain:-t_species) -> ANVIO_SUMMARY_DAS
```

```
## Rows: 136 Columns: 14
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
<div id="htmlwidget-892a0e24176d99b7b833" style="width:100%;height:auto;" class="datatables html-widget"></div>
<script type="application/json" data-for="htmlwidget-892a0e24176d99b7b833">{"x":{"filter":"none","vertical":false,"data":[["1","2","3","4","5","6"],["maxbin_002_fasta_contigs","maxbin_004_fasta_contigs","maxbin_005_fasta_contigs","maxbin_006_fasta_contigs","maxbin_008_fasta_contigs","maxbin_009_fasta_contigs"],[3373122,2078616,3887172,2170300,2596485,2676589],[99,36,465,288,54,325],[65526,155262,22282,13284,117427,14386],[44.100777643127,40.7792212729088,57.8414339713452,41.1785262144328,40.7337175045291,37.0186970597009],[88.7323943661972,91.5492957746479,80.2816901408451,39.4366197183099,64.7887323943662,95.7746478873239],[2.8169014084507,4.22535211267606,22.5352112676056,0,11.2676056338028,8.45070422535211],["HQ","HQ",null,null,null,"HQ"]],"container":"<table class=\"display\">\n  <thead>\n    <tr>\n      <th> <\/th>\n      <th>ANVIO_SUMMARY_DAS_bins<\/th>\n      <th>ANVIO_SUMMARY_DAS_total_length<\/th>\n      <th>ANVIO_SUMMARY_DAS_num_contigs<\/th>\n      <th>ANVIO_SUMMARY_DAS_N50<\/th>\n      <th>ANVIO_SUMMARY_DAS_GC_content<\/th>\n      <th>ANVIO_SUMMARY_DAS_percent_completion<\/th>\n      <th>ANVIO_SUMMARY_DAS_percent_redundancy<\/th>\n      <th>ANVIO_DAS_HQ<\/th>\n    <\/tr>\n  <\/thead>\n<\/table>","options":{"columnDefs":[{"className":"dt-right","targets":[2,3,4,5,6,7]},{"orderable":false,"targets":0}],"order":[],"autoWidth":false,"orderClasses":false}},"evals":[],"jsHooks":[]}</script>
```


```r
orf_contig_3 %>% 
  left_join(ANVIO_SUMMARY_DAS,
            by = c("DAS_bin" = "ANVIO_SUMMARY_DAS_bins")) -> orf_contig_3
```

#### ANVIO CONCOCT info:


```r
ANVIO_SUMMARY_CONCOCT %>% 
  here::here() %>% 
  read_tsv() %>% 
  select(-t_domain:-t_species) -> ANVIO_SUMMARY_CONCOCT
```

```
## Rows: 129 Columns: 14
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
<div id="htmlwidget-0186194d8b8157e8a0ac" style="width:100%;height:auto;" class="datatables html-widget"></div>
<script type="application/json" data-for="htmlwidget-0186194d8b8157e8a0ac">{"x":{"filter":"none","vertical":false,"data":[["1","2","3","4","5","6"],["Bin_0_1_1","Bin_0_2_1","Bin_0_3_1","Bin_0_4_1","Bin_0_5_1","Bin_0_6_1"],[4997430,2611266,1816009,3210323,3351680,2332080],[93,118,22,74,21,31],[95529,48119,147912,71289,221494,117575],[44.8850463626556,45.7147524370042,44.5251091026598,44.3226862025176,59.2794959337204,43.7523527454812],[98.5915492957746,98.5915492957746,88.7323943661972,85.9154929577465,80.2816901408451,100],[1.40845070422535,0,0,2.8169014084507,1.40845070422535,1.40845070422535],["HQ","HQ","HQ","HQ","HQ","HQ"]],"container":"<table class=\"display\">\n  <thead>\n    <tr>\n      <th> <\/th>\n      <th>ANVIO_SUMMARY_CONCOCT_bins<\/th>\n      <th>ANVIO_SUMMARY_CONCOCT_total_length<\/th>\n      <th>ANVIO_SUMMARY_CONCOCT_num_contigs<\/th>\n      <th>ANVIO_SUMMARY_CONCOCT_N50<\/th>\n      <th>ANVIO_SUMMARY_CONCOCT_GC_content<\/th>\n      <th>ANVIO_SUMMARY_CONCOCT_percent_completion<\/th>\n      <th>ANVIO_SUMMARY_CONCOCT_percent_redundancy<\/th>\n      <th>ANVIO_CONCOCT_HQ<\/th>\n    <\/tr>\n  <\/thead>\n<\/table>","options":{"columnDefs":[{"className":"dt-right","targets":[2,3,4,5,6,7]},{"orderable":false,"targets":0}],"order":[],"autoWidth":false,"orderClasses":false}},"evals":[],"jsHooks":[]}</script>
```


```r
orf_contig_3 %>% 
  left_join(ANVIO_SUMMARY_CONCOCT,
            by = c("CONCOCT_bin" = "ANVIO_SUMMARY_CONCOCT_bins")) -> orf_contig_3
```


Rename ORF and contig ID:

```r
orf_contig_3 %>% 
  mutate(ORF_ID = paste0("h1_", ORF_ID),
         Contig_ID = paste0("h1_", Contig_ID)) -> full_gene_catalogue
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
<div id="htmlwidget-dd6d32dfe4e4c1d0a3a8" style="width:100%;height:auto;" class="datatables html-widget"></div>
<script type="application/json" data-for="htmlwidget-dd6d32dfe4e4c1d0a3a8">{"x":{"filter":"none","vertical":false,"data":[["1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","23","24","25","26","27","28","29","30","31","32","33","34","35","36","37","38","39","40","41","42","43","44","45","46","47","48","49","50","51","52","53","54","55","56","57","58","59","60","61","62","63","64","65","66","67","68","69","70","71","72","73","74","75","76","77","78","79","80","81","82","83","84","85","86","87","88","89","90","91","92","93","94","95","96","97","98","99","100","101","102","103","104","105","106","107","108","109","110","111","112","113","114","115","116","117","118","119","120","121","122","123","124","125","126","127","128","129"],["h1_megahit_3458_1913-2263","h1_megahit_3700_1529-2170","h1_megahit_4077_533-1408","h1_megahit_4383_62-1096","h1_megahit_4383_1114-2592","h1_megahit_5059_1937-2077","h1_megahit_7781_5471-6166","h1_megahit_7781_6204-7298","h1_megahit_7781_7513-8481","h1_megahit_7781_8456-9505","h1_megahit_7781_9511-10119","h1_megahit_7781_10547-11458","h1_megahit_7781_11611-12096","h1_megahit_8637_6630-6800","h1_megahit_10816_3-1409","h1_megahit_12983_96-836","h1_megahit_13043_3-1073","h1_megahit_14818_1-867","h1_megahit_17429_3-1667","h1_megahit_21838_49-1107","h1_megahit_22423_368-1294","h1_megahit_22665_178-801","h1_megahit_24486_104-1030","h1_megahit_26875_13-837","h1_megahit_31365_2-67","h1_megahit_32317_1113-2249","h1_megahit_34365_273-995","h1_megahit_35417_11042-11755","h1_megahit_38276_5712-5939","h1_megahit_43660_352-1113","h1_megahit_43950_3-968","h1_megahit_44613_2883-3806","h1_megahit_45838_3-1151","h1_megahit_47082_3-242","h1_megahit_47706_2293-2853","h1_megahit_48851_391-1086","h1_megahit_50927_1443-1502","h1_megahit_51160_1974-3914","h1_megahit_56909_6180-6836","h1_megahit_59954_1177-1611","h1_megahit_60788_49251-49448","h1_megahit_60864_1-60","h1_megahit_66397_331-525","h1_megahit_69239_1-486","h1_megahit_72350_6904-8370","h1_megahit_76455_2170-3036","h1_megahit_76882_146-1021","h1_megahit_76972_10094-10252","h1_megahit_80040_12135-13358","h1_megahit_80612_1886-2809","h1_megahit_82493_268-1194","h1_megahit_88117_6614-7453","h1_megahit_88117_7958-8041","h1_megahit_88784_1007-1933","h1_megahit_89242_3-416","h1_megahit_89242_832-1179","h1_megahit_90033_128-943","h1_megahit_90033_2422-3282","h1_megahit_92101_2013-2930","h1_megahit_92307_3933-4562","h1_megahit_92462_5562-5726","h1_megahit_96881_1-174","h1_megahit_97941_8380-8610","h1_megahit_97941_13243-13473","h1_megahit_101494_48-836","h1_megahit_104003_1-336","h1_megahit_104003_333-1736","h1_megahit_108574_295-516","h1_megahit_111254_1984-2907","h1_megahit_111254_6015-6257","h1_megahit_112028_3-152","h1_megahit_116885_2-1588","h1_megahit_117104_1356-2321","h1_megahit_119962_28181-28585","h1_megahit_121143_236-1663","h1_megahit_121608_236-1663","h1_megahit_123691_1441-1842","h1_megahit_128451_3-710","h1_megahit_128765_1-681","h1_megahit_128765_739-1416","h1_megahit_128887_313-1563","h1_megahit_130881_1-681","h1_megahit_130881_774-2033","h1_megahit_130881_2295-2378","h1_megahit_134352_5274-5693","h1_megahit_140180_2-130","h1_megahit_146089_3-632","h1_megahit_147156_1-81","h1_megahit_148486_3-350","h1_megahit_159299_357-2270","h1_megahit_167681_1839-2189","h1_megahit_172549_1950-2195","h1_megahit_174172_416-907","h1_megahit_175793_1720-2808","h1_megahit_175793_2795-3493","h1_megahit_185618_1-60","h1_megahit_188016_651-1115","h1_megahit_188016_1201-1524","h1_megahit_191331_3-122","h1_megahit_192793_821-1372","h1_megahit_195133_149-985","h1_megahit_195133_985-1404","h1_megahit_195133_1849-2664","h1_megahit_196698_3-98","h1_megahit_196936_743-973","h1_megahit_203866_1698-1934","h1_megahit_204983_405-1085","h1_megahit_208133_14655-15578","h1_megahit_209006_1-894","h1_megahit_209006_914-2005","h1_megahit_209006_2105-2791","h1_megahit_209753_1462-1866","h1_megahit_210647_1-60","h1_megahit_223849_916-1113","h1_megahit_223975_1-93","h1_megahit_224114_536-1132","h1_megahit_224435_756-1346","h1_megahit_224474_1664-1873","h1_megahit_225321_1-87","h1_megahit_229751_3-158","h1_megahit_237523_493-1011","h1_megahit_242168_461-1084","h1_megahit_242168_1166-2371","h1_megahit_242168_2484-3077","h1_megahit_242168_3165-3581","h1_megahit_245522_2335-2661","h1_megahit_246135_54-1262","h1_megahit_247479_472-1671","h1_megahit_247727_590-796"],["FOSX","VATE","CTX-M-1","CHLORAMPHENICOL_AND_FLORFENICOL_RESISTANCE","LSAE","RLMA(II)","VANRA","VANSA","VANHA","VANA","VANXA","VANYA","VANZA","CHRB","MDTO","OMPR","BICYCLOMYCIN-MULTIDRUG_EFFLUX_PROTEIN_BCR","BCRA","MSBA","OMPF","BCRA","CAMPYLOBACTER COLI CHLORAMPHENICOL ACETYLTRANSFERASE","BCRA","GADX","MYRA","OMPR","BAER","VANR","AADA14","ERMB","UGD","BCRA","ESCHERICHIA COLI AMPH BETA-LACTAMASE","KDPE","BACA","VANR","TMB-2","TETS","VANR","ESCHERICHIA COLI MARR MUTANT CONFERRING ANTIBIOTIC RESISTANCE","CRPP","MDSC","APH(3')-VC","ESCHERICHIA COLI EMRE","EMRB-QACA_FAMILY_MAJOR_FACILITATOR_TRANSPORTER","AADE","ESCHERICHIA COLI AMPC BETA-LACTAMASE","VANRM","VANU","BCRA","BCRA","SUL1","AADA12","BCRA","EMRK","EVGA","APH(3')-IA","TEM-1","BCRA","VATG","SUL2","SMEE","VGAC","VGAC","AADA5","BAER","BAES","VANU","BCRA","VANU","QNRS3","TET(W/N/W)","CFXA2","DNA-BINDING_PROTEIN_H-NS","TET(W/N/W)","TET(W/N/W)","DNA-BINDING_PROTEIN_H-NS","VANR","AADA3","ANT(2'')-IA","CHLORAMPHENICOL_AND_FLORFENICOL_RESISTANCE","AADA15","CMLA6","AADA12","FOSB3","TET(X3)","TETM","OPRM","TRUNCATED_PUTATIVE_RESPONSE_REGULATOR_ARLR","TET(W/N/W)","TETW","PP-FLO","LNUA","VANS","VANR","MDSC","ESCHERICHIA COLI SOXR WITH MUTATION CONFERRING ANTIBIOTIC RESISTANCE","ESCHERICHIA COLI SOXS WITH MUTATION CONFERRING ANTIBIOTIC RESISTANCE","EMRY","APH(6)-IC","APH(6)-ID","APH(3'')-IB","SUL2","OXA-43","VGAC","VANU","VATB","BCRA","VANHD","VANYD","VANSD","DNA-BINDING_PROTEIN_H-NS","MDSC","DFRG","MDTM","VANZ","VANR","VANU","VANN","EREB","LNUE","TETR","TET(B)","TETC","TETD","ESCHERICHIA COLI SOXS WITH MUTATION CONFERRING ANTIBIOTIC RESISTANCE","EMRD","UGD","VANU"],["fosfomycin","MLS","beta-lactam","phenicol","multidrug","MLS","glycopeptide","glycopeptide","glycopeptide","glycopeptide","glycopeptide","glycopeptide","glycopeptide","MLS","multidrug","beta-lactam","multidrug","bacitracin","nitroimidazole","beta-lactam","bacitracin","phenicol","bacitracin","multidrug","MLS","beta-lactam","aminoglycoside:aminocoumarin","glycopeptide","aminoglycoside","MLS","peptide","bacitracin","beta-lactam","aminoglycoside","bacitracin","glycopeptide","beta-lactam","tetracycline","glycopeptide","multidrug","fluoroquinolone","multidrug","aminoglycoside","MLS","fluoroquinolone","aminoglycoside","beta-lactam","glycopeptide","glycopeptide","bacitracin","bacitracin","sulfonamide","aminoglycoside","bacitracin","tetracycline","multidrug","aminoglycoside","beta-lactam","bacitracin","MLS","sulfonamide","multidrug","multidrug","multidrug","aminoglycoside","aminoglycoside:aminocoumarin","aminoglycoside:aminocoumarin","glycopeptide","bacitracin","glycopeptide","fluoroquinolone","tetracycline","beta-lactam","unclassified","tetracycline","tetracycline","unclassified","glycopeptide","aminoglycoside","aminoglycoside","phenicol","aminoglycoside","phenicol","aminoglycoside","fosfomycin","tetracycline","tetracycline","multidrug","unclassified","tetracycline","tetracycline","phenicol","MLS","glycopeptide","glycopeptide","multidrug","multidrug","multidrug","tetracycline","aminoglycoside","aminoglycoside","aminoglycoside","sulfonamide","beta-lactam","multidrug","glycopeptide","MLS","bacitracin","glycopeptide","glycopeptide","glycopeptide","unclassified","multidrug","diaminopyrimidine","multidrug","glycopeptide","glycopeptide","glycopeptide","glycopeptide","MLS","MLS","tetracycline","tetracycline","tetracycline","tetracycline","multidrug","multidrug","peptide","glycopeptide"],["Uncertain - plasmid or chromosomal","Uncertain - plasmid or chromosomal","Uncertain - plasmid or chromosomal","Uncertain - plasmid or chromosomal","Uncertain - plasmid or chromosomal","Uncertain - plasmid or chromosomal","Plasmid","Plasmid","Plasmid","Plasmid","Plasmid","Plasmid","Plasmid","Uncertain - plasmid or chromosomal","Uncertain - too short","Uncertain - plasmid or chromosomal","Uncertain - too short","Uncertain - plasmid or chromosomal","Uncertain - too short","Uncertain - plasmid or chromosomal","Uncertain - plasmid or chromosomal","Uncertain - too short","Uncertain - plasmid or chromosomal","Uncertain - too short","Uncertain - too short","Uncertain - plasmid or chromosomal","Uncertain - plasmid or chromosomal","Uncertain - plasmid or chromosomal","Uncertain - plasmid or chromosomal","Plasmid","Uncertain - too short","Uncertain - plasmid or chromosomal","Uncertain - too short","Uncertain - plasmid or chromosomal","Plasmid","Uncertain - plasmid or chromosomal","Uncertain - plasmid or chromosomal","Uncertain - plasmid or chromosomal","Uncertain - plasmid or chromosomal","Uncertain - plasmid or chromosomal","Plasmid","Uncertain - plasmid or chromosomal","Uncertain - too short","Uncertain - plasmid or chromosomal","Uncertain - plasmid or chromosomal","Uncertain - plasmid or chromosomal","Uncertain - too short","Uncertain - plasmid or chromosomal","Plasmid","Uncertain - plasmid or chromosomal","Uncertain - plasmid or chromosomal","Plasmid","Plasmid","Uncertain - plasmid or chromosomal","Uncertain - too short","Uncertain - too short","Plasmid","Plasmid","Uncertain - plasmid or chromosomal","Plasmid","Uncertain - plasmid or chromosomal","Uncertain - too short","Plasmid","Plasmid","Uncertain - plasmid or chromosomal","Uncertain - plasmid or chromosomal","Uncertain - plasmid or chromosomal","Uncertain - too short","Uncertain - plasmid or chromosomal","Uncertain - plasmid or chromosomal","Uncertain - plasmid or chromosomal","Uncertain - too short","Uncertain - plasmid or chromosomal","Plasmid","Uncertain - plasmid or chromosomal","Uncertain - plasmid or chromosomal","Uncertain - plasmid or chromosomal","Uncertain - plasmid or chromosomal","Uncertain - too short","Uncertain - too short","Uncertain - too short","Uncertain - too short","Uncertain - too short","Uncertain - too short","Uncertain - plasmid or chromosomal","Uncertain - too short","Plasmid","Uncertain - too short","Uncertain - plasmid or chromosomal","Uncertain - too short","Uncertain - too short","Uncertain - plasmid or chromosomal","Uncertain - plasmid or chromosomal","Uncertain - plasmid or chromosomal","Uncertain - plasmid or chromosomal","Uncertain - plasmid or chromosomal","Uncertain - plasmid or chromosomal","Uncertain - plasmid or chromosomal","Uncertain - too short","Uncertain - plasmid or chromosomal","Plasmid","Plasmid","Plasmid","Uncertain - too short","Uncertain - plasmid or chromosomal","Uncertain - too short","Uncertain - too short","Uncertain - plasmid or chromosomal","Uncertain - plasmid or chromosomal","Uncertain - plasmid or chromosomal","Uncertain - plasmid or chromosomal","Uncertain - plasmid or chromosomal","Uncertain - too short","Uncertain - too short","Uncertain - too short","Uncertain - too short","Uncertain - plasmid or chromosomal","Uncertain - plasmid or chromosomal","Uncertain - too short","Uncertain - too short","Uncertain - too short","Plasmid","Plasmid","Plasmid","Plasmid","Uncertain - plasmid or chromosomal","Uncertain - too short","Uncertain - plasmid or chromosomal","Uncertain - too short"],["ambiguous (plasmid/phage)","plasmid","plasmid","plasmid","plasmid","plasmid","plasmid","plasmid","plasmid","plasmid","plasmid","plasmid","plasmid","ambiguous (plasmid/phage)","plasmid","plasmid","plasmid","plasmid","plasmid","plasmid","ambiguous (plasmid/phage)","plasmid","plasmid","plasmid","plasmid","plasmid","plasmid","plasmid","plasmid","plasmid","plasmid","ambiguous (plasmid/phage)","plasmid","plasmid","plasmid","plasmid","plasmid","plasmid","plasmid","plasmid","plasmid","plasmid","plasmid","plasmid","plasmid","ambiguous (plasmid/phage)","plasmid","plasmid","plasmid","plasmid","plasmid","plasmid","plasmid","plasmid","plasmid","plasmid","plasmid","plasmid","plasmid","ambiguous (plasmid/phage)","plasmid","plasmid","plasmid","plasmid","plasmid","plasmid","plasmid","plasmid","ambiguous (plasmid/phage)","ambiguous (plasmid/phage)","plasmid","plasmid","plasmid","ambiguous (plasmid/phage)","plasmid","plasmid","plasmid","ambiguous (plasmid/phage)","plasmid","plasmid","plasmid","plasmid","plasmid","plasmid","ambiguous (plasmid/phage)","plasmid","ambiguous (plasmid/phage)","plasmid","ambiguous (plasmid/phage)","plasmid","plasmid","plasmid","plasmid","plasmid","plasmid","plasmid","plasmid","plasmid","plasmid","plasmid","plasmid","plasmid","plasmid","plasmid","plasmid","plasmid","plasmid","ambiguous (plasmid/phage)","plasmid","plasmid","plasmid","plasmid","plasmid","plasmid","plasmid","plasmid","plasmid","plasmid","plasmid","plasmid","plasmid","plasmid","plasmid","plasmid","plasmid","plasmid","plasmid","plasmid","plasmid"],["maxbin_167_fasta_contigs","maxbin_014_fasta_contigs","maxbin_140_fasta_contigs","maxbin_065_fasta_contigs","maxbin_065_fasta_contigs","maxbin_157_fasta_contigs",null,null,null,null,null,null,null,"maxbin_034_fasta_contigs",null,"maxbin_088_fasta_contigs",null,"maxbin_133_fasta_contigs",null,null,"maxbin_066_fasta_contigs",null,null,null,"maxbin_128_fasta_contigs","maxbin_110_fasta_contigs",null,"maxbin_058_fasta_contigs","maxbin_137_fasta_contigs","maxbin_010_fasta_contigs",null,"maxbin_107_fasta_contigs",null,"maxbin_110_fasta_contigs","maxbin_165_fasta_contigs",null,"maxbin_036_fasta_contigs","maxbin_154_fasta_contigs","maxbin_122_fasta_contigs",null,"maxbin_166_fasta_contigs","maxbin_165_fasta_contigs",null,null,"maxbin_162_fasta_contigs","maxbin_013_fasta_contigs",null,"maxbin_055_fasta_sub_contigs","maxbin_084_fasta_contigs",null,"maxbin_122_fasta_contigs",null,null,"maxbin_006_fasta_contigs",null,null,null,null,"maxbin_136_fasta_contigs",null,null,"maxbin_140_fasta_contigs",null,null,null,null,null,null,"maxbin_164_fasta_contigs","maxbin_164_fasta_contigs","maxbin_074_fasta_contigs","maxbin_013_fasta_contigs","maxbin_017_fasta_contigs",null,"maxbin_143_fasta_contigs","maxbin_167_fasta_contigs",null,null,null,null,null,null,null,null,"maxbin_026_fasta_contigs","maxbin_035_fasta_contigs","maxbin_143_fasta_contigs",null,null,"maxbin_138_fasta_contigs",null,"maxbin_166_fasta_contigs","maxbin_125_fasta_contigs","maxbin_055_fasta_sub_contigs","maxbin_055_fasta_sub_contigs","maxbin_165_fasta_contigs","maxbin_110_fasta_contigs","maxbin_110_fasta_contigs",null,"maxbin_053_fasta_contigs",null,null,null,null,null,"maxbin_152_fasta_contigs",null,"maxbin_164_fasta_contigs","maxbin_125_fasta_contigs","maxbin_125_fasta_contigs","maxbin_125_fasta_contigs",null,"maxbin_165_fasta_contigs",null,"maxbin_135_fasta_contigs","maxbin_013_fasta_contigs",null,null,null,null,null,"maxbin_140_fasta_contigs","maxbin_140_fasta_contigs","maxbin_140_fasta_contigs","maxbin_140_fasta_contigs",null,null,"maxbin_098_fasta_contigs",null],["Bin_0_9_1",null,null,null,null,null,null,null,null,null,null,null,null,"Bin_5_5",null,null,null,null,null,null,null,null,"Bin_2_26",null,null,null,null,null,"Bin_2_21",null,null,"Bin_9_21",null,"Bin_7_4_1",null,null,null,null,null,null,"Bin_7_3","Bin_0_2_1",null,null,"Bin_9_9_1","Bin_5_15_1",null,"Bin_5_7_1","Bin_4_1_1","Bin_6_12",null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,"Bin_6_15",null,"Bin_6_14","Bin_0_9_1",null,null,null,null,"Bin_6_10",null,null,null,"Bin_0_6_1",null,null,null,null,null,null,"Bin_7_3",null,null,null,"Bin_0_2_1","Bin_7_4_1","Bin_7_4_1",null,"Bin_2_1",null,null,null,null,null,"Bin_2_20",null,null,null,null,null,null,"Bin_0_2_1",null,"Bin_2_10",null,null,"Bin_9_15_1",null,null,"Bin_8_4_1",null,null,null,null,null,null,"Bin_6_8",null]],"container":"<table class=\"display\">\n  <thead>\n    <tr>\n      <th> <\/th>\n      <th>ORF_ID<\/th>\n      <th>PathoFact_AMR_ARG<\/th>\n      <th>PathoFact_AMR_AMR_category<\/th>\n      <th>viralverify_Prediction<\/th>\n      <th>PathoFact_AMR_MGE_prediction<\/th>\n      <th>DAS_bin<\/th>\n      <th>CONCOCT_bin<\/th>\n    <\/tr>\n  <\/thead>\n<\/table>","options":{"columnDefs":[{"orderable":false,"targets":0}],"order":[],"autoWidth":false,"orderClasses":false}},"evals":[],"jsHooks":[]}</script>
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
  # select(ORF_ID, TPM.D.1:Raw.base.count.HV.24) %>% 
    select(ORF_ID, starts_with(c("TPM", "Coverage", "Raw.base.count", "Raw.read.count."))) %>% 
  # select(starts_with("TPM"), "ORF_ID") %>% 
  column_to_rownames("ORF_ID") -> quant_orfs

# colnames(quant) <- colnames(quant) %>% 
#   str_replace_all("//.", "//_")

rownames(quant_orfs) <- paste0("h1_",
                               rownames(quant_orfs))

quant_orfs %>% 
  rownames_to_column("ORF_ID") -> quant_orfs

quant_orfs %>% 
  head() %>% 
  DT::datatable()
```

```{=html}
<div id="htmlwidget-8284da1c631777b04da7" style="width:100%;height:auto;" class="datatables html-widget"></div>
<script type="application/json" data-for="htmlwidget-8284da1c631777b04da7">{"x":{"filter":"none","vertical":false,"data":[["1","2","3","4","5","6"],["h1_megahit_1_1-312","h1_megahit_2_3-176","h1_megahit_2_173-304","h1_megahit_3_3-116","h1_megahit_3_82-219","h1_megahit_4_1-72"],[0.165,0,0,0,0,0],[0,0.229,0.151,0,0,0],[0,0,0,0.298,0.246,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0.23,0.304,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0.405],[0,0,0,0.156,0.129,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0.178,0.147,0],[0,0,0,0.141,0.116,0],[0,0,0,0,0,0],[0,0,0,0.154,0.127,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0.553],[0,0,0,0,0,0],[0,0,0,0.165,0.137,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0.434,0.358,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[3.349,0,0,0,0,0],[0,1.149,0.78,0,0,0],[0,0,0,1.684,1.261,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,1.172,0.773,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,1.972],[0,0,0,0.404,0.993,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0.86,0.623,0],[0,0,0,0.746,0.717,0],[0,0,0,0,0,0],[0,0,0,0.561,0.87,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0.056],[0,0,0,0,0,0],[0,0,0,0.096,0.768,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,1.974,2.37,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[1045,0,0,0,0,0],[0,200,103,0,0,0],[0,0,0,192,174,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,204,102,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,142],[0,0,0,46,137,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,98,86,0],[0,0,0,85,99,0],[0,0,0,0,0,0],[0,0,0,64,120,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,4],[0,0,0,0,0,0],[0,0,0,11,106,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,225,327,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[7,0,0,0,0,0],[0,2,1,0,0,0],[0,0,0,2,2,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,2,2,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,2],[0,0,0,1,1,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,1,1,0],[0,0,0,1,1,0],[0,0,0,0,0,0],[0,0,0,1,1,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,2],[0,0,0,0,0,0],[0,0,0,1,1,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,3,3,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0]],"container":"<table class=\"display\">\n  <thead>\n    <tr>\n      <th> <\/th>\n      <th>ORF_ID<\/th>\n      <th>TPM.D.1<\/th>\n      <th>TPM.HC.1<\/th>\n      <th>TPM.HC.2<\/th>\n      <th>TPM.HC.3<\/th>\n      <th>TPM.HC.4<\/th>\n      <th>TPM.HC.5<\/th>\n      <th>TPM.HC.6<\/th>\n      <th>TPM.HC.7<\/th>\n      <th>TPM.HC.8<\/th>\n      <th>TPM.HC.9<\/th>\n      <th>TPM.HC.10<\/th>\n      <th>TPM.HC.11<\/th>\n      <th>TPM.HC.12<\/th>\n      <th>TPM.HC.13<\/th>\n      <th>TPM.HC.14<\/th>\n      <th>TPM.HC.15<\/th>\n      <th>TPM.HC.16<\/th>\n      <th>TPM.HC.17<\/th>\n      <th>TPM.HC.18<\/th>\n      <th>TPM.HC.19<\/th>\n      <th>TPM.HC.20<\/th>\n      <th>TPM.HC.21<\/th>\n      <th>TPM.HC.22<\/th>\n      <th>TPM.HC.23<\/th>\n      <th>TPM.HC.24<\/th>\n      <th>TPM.HV.1<\/th>\n      <th>TPM.HV.2<\/th>\n      <th>TPM.HV.3<\/th>\n      <th>TPM.HV.4<\/th>\n      <th>TPM.HV.5<\/th>\n      <th>TPM.HV.6<\/th>\n      <th>TPM.HV.7<\/th>\n      <th>TPM.HV.8<\/th>\n      <th>TPM.HV.9<\/th>\n      <th>TPM.HV.10<\/th>\n      <th>TPM.HV.11<\/th>\n      <th>TPM.HV.12<\/th>\n      <th>TPM.HV.13<\/th>\n      <th>TPM.HV.14<\/th>\n      <th>TPM.HV.15<\/th>\n      <th>TPM.HV.16<\/th>\n      <th>TPM.HV.17<\/th>\n      <th>TPM.HV.18<\/th>\n      <th>TPM.HV.19<\/th>\n      <th>TPM.HV.20<\/th>\n      <th>TPM.HV.21<\/th>\n      <th>TPM.HV.22<\/th>\n      <th>TPM.HV.23<\/th>\n      <th>TPM.HV.24<\/th>\n      <th>Coverage.D.1<\/th>\n      <th>Coverage.HC.1<\/th>\n      <th>Coverage.HC.2<\/th>\n      <th>Coverage.HC.3<\/th>\n      <th>Coverage.HC.4<\/th>\n      <th>Coverage.HC.5<\/th>\n      <th>Coverage.HC.6<\/th>\n      <th>Coverage.HC.7<\/th>\n      <th>Coverage.HC.8<\/th>\n      <th>Coverage.HC.9<\/th>\n      <th>Coverage.HC.10<\/th>\n      <th>Coverage.HC.11<\/th>\n      <th>Coverage.HC.12<\/th>\n      <th>Coverage.HC.13<\/th>\n      <th>Coverage.HC.14<\/th>\n      <th>Coverage.HC.15<\/th>\n      <th>Coverage.HC.16<\/th>\n      <th>Coverage.HC.17<\/th>\n      <th>Coverage.HC.18<\/th>\n      <th>Coverage.HC.19<\/th>\n      <th>Coverage.HC.20<\/th>\n      <th>Coverage.HC.21<\/th>\n      <th>Coverage.HC.22<\/th>\n      <th>Coverage.HC.23<\/th>\n      <th>Coverage.HC.24<\/th>\n      <th>Coverage.HV.1<\/th>\n      <th>Coverage.HV.2<\/th>\n      <th>Coverage.HV.3<\/th>\n      <th>Coverage.HV.4<\/th>\n      <th>Coverage.HV.5<\/th>\n      <th>Coverage.HV.6<\/th>\n      <th>Coverage.HV.7<\/th>\n      <th>Coverage.HV.8<\/th>\n      <th>Coverage.HV.9<\/th>\n      <th>Coverage.HV.10<\/th>\n      <th>Coverage.HV.11<\/th>\n      <th>Coverage.HV.12<\/th>\n      <th>Coverage.HV.13<\/th>\n      <th>Coverage.HV.14<\/th>\n      <th>Coverage.HV.15<\/th>\n      <th>Coverage.HV.16<\/th>\n      <th>Coverage.HV.17<\/th>\n      <th>Coverage.HV.18<\/th>\n      <th>Coverage.HV.19<\/th>\n      <th>Coverage.HV.20<\/th>\n      <th>Coverage.HV.21<\/th>\n      <th>Coverage.HV.22<\/th>\n      <th>Coverage.HV.23<\/th>\n      <th>Coverage.HV.24<\/th>\n      <th>Raw.base.count.D.1<\/th>\n      <th>Raw.base.count.HC.1<\/th>\n      <th>Raw.base.count.HC.2<\/th>\n      <th>Raw.base.count.HC.3<\/th>\n      <th>Raw.base.count.HC.4<\/th>\n      <th>Raw.base.count.HC.5<\/th>\n      <th>Raw.base.count.HC.6<\/th>\n      <th>Raw.base.count.HC.7<\/th>\n      <th>Raw.base.count.HC.8<\/th>\n      <th>Raw.base.count.HC.9<\/th>\n      <th>Raw.base.count.HC.10<\/th>\n      <th>Raw.base.count.HC.11<\/th>\n      <th>Raw.base.count.HC.12<\/th>\n      <th>Raw.base.count.HC.13<\/th>\n      <th>Raw.base.count.HC.14<\/th>\n      <th>Raw.base.count.HC.15<\/th>\n      <th>Raw.base.count.HC.16<\/th>\n      <th>Raw.base.count.HC.17<\/th>\n      <th>Raw.base.count.HC.18<\/th>\n      <th>Raw.base.count.HC.19<\/th>\n      <th>Raw.base.count.HC.20<\/th>\n      <th>Raw.base.count.HC.21<\/th>\n      <th>Raw.base.count.HC.22<\/th>\n      <th>Raw.base.count.HC.23<\/th>\n      <th>Raw.base.count.HC.24<\/th>\n      <th>Raw.base.count.HV.1<\/th>\n      <th>Raw.base.count.HV.2<\/th>\n      <th>Raw.base.count.HV.3<\/th>\n      <th>Raw.base.count.HV.4<\/th>\n      <th>Raw.base.count.HV.5<\/th>\n      <th>Raw.base.count.HV.6<\/th>\n      <th>Raw.base.count.HV.7<\/th>\n      <th>Raw.base.count.HV.8<\/th>\n      <th>Raw.base.count.HV.9<\/th>\n      <th>Raw.base.count.HV.10<\/th>\n      <th>Raw.base.count.HV.11<\/th>\n      <th>Raw.base.count.HV.12<\/th>\n      <th>Raw.base.count.HV.13<\/th>\n      <th>Raw.base.count.HV.14<\/th>\n      <th>Raw.base.count.HV.15<\/th>\n      <th>Raw.base.count.HV.16<\/th>\n      <th>Raw.base.count.HV.17<\/th>\n      <th>Raw.base.count.HV.18<\/th>\n      <th>Raw.base.count.HV.19<\/th>\n      <th>Raw.base.count.HV.20<\/th>\n      <th>Raw.base.count.HV.21<\/th>\n      <th>Raw.base.count.HV.22<\/th>\n      <th>Raw.base.count.HV.23<\/th>\n      <th>Raw.base.count.HV.24<\/th>\n      <th>Raw.read.count.D.1<\/th>\n      <th>Raw.read.count.HC.1<\/th>\n      <th>Raw.read.count.HC.2<\/th>\n      <th>Raw.read.count.HC.3<\/th>\n      <th>Raw.read.count.HC.4<\/th>\n      <th>Raw.read.count.HC.5<\/th>\n      <th>Raw.read.count.HC.6<\/th>\n      <th>Raw.read.count.HC.7<\/th>\n      <th>Raw.read.count.HC.8<\/th>\n      <th>Raw.read.count.HC.9<\/th>\n      <th>Raw.read.count.HC.10<\/th>\n      <th>Raw.read.count.HC.11<\/th>\n      <th>Raw.read.count.HC.12<\/th>\n      <th>Raw.read.count.HC.13<\/th>\n      <th>Raw.read.count.HC.14<\/th>\n      <th>Raw.read.count.HC.15<\/th>\n      <th>Raw.read.count.HC.16<\/th>\n      <th>Raw.read.count.HC.17<\/th>\n      <th>Raw.read.count.HC.18<\/th>\n      <th>Raw.read.count.HC.19<\/th>\n      <th>Raw.read.count.HC.20<\/th>\n      <th>Raw.read.count.HC.21<\/th>\n      <th>Raw.read.count.HC.22<\/th>\n      <th>Raw.read.count.HC.23<\/th>\n      <th>Raw.read.count.HC.24<\/th>\n      <th>Raw.read.count.HV.1<\/th>\n      <th>Raw.read.count.HV.2<\/th>\n      <th>Raw.read.count.HV.3<\/th>\n      <th>Raw.read.count.HV.4<\/th>\n      <th>Raw.read.count.HV.5<\/th>\n      <th>Raw.read.count.HV.6<\/th>\n      <th>Raw.read.count.HV.7<\/th>\n      <th>Raw.read.count.HV.8<\/th>\n      <th>Raw.read.count.HV.9<\/th>\n      <th>Raw.read.count.HV.10<\/th>\n      <th>Raw.read.count.HV.11<\/th>\n      <th>Raw.read.count.HV.12<\/th>\n      <th>Raw.read.count.HV.13<\/th>\n      <th>Raw.read.count.HV.14<\/th>\n      <th>Raw.read.count.HV.15<\/th>\n      <th>Raw.read.count.HV.16<\/th>\n      <th>Raw.read.count.HV.17<\/th>\n      <th>Raw.read.count.HV.18<\/th>\n      <th>Raw.read.count.HV.19<\/th>\n      <th>Raw.read.count.HV.20<\/th>\n      <th>Raw.read.count.HV.21<\/th>\n      <th>Raw.read.count.HV.22<\/th>\n      <th>Raw.read.count.HV.23<\/th>\n      <th>Raw.read.count.HV.24<\/th>\n    <\/tr>\n  <\/thead>\n<\/table>","options":{"columnDefs":[{"className":"dt-right","targets":[2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58,59,60,61,62,63,64,65,66,67,68,69,70,71,72,73,74,75,76,77,78,79,80,81,82,83,84,85,86,87,88,89,90,91,92,93,94,95,96,97,98,99,100,101,102,103,104,105,106,107,108,109,110,111,112,113,114,115,116,117,118,119,120,121,122,123,124,125,126,127,128,129,130,131,132,133,134,135,136,137,138,139,140,141,142,143,144,145,146,147,148,149,150,151,152,153,154,155,156,157,158,159,160,161,162,163,164,165,166,167,168,169,170,171,172,173,174,175,176,177,178,179,180,181,182,183,184,185,186,187,188,189,190,191,192,193,194,195,196,197]},{"orderable":false,"targets":0}],"order":[],"autoWidth":false,"orderClasses":false}},"evals":[],"jsHooks":[]}</script>
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
  # select(starts_with("TPM"), "ORF_ID") %>% 
  column_to_rownames("Contig_ID") -> quant_contigs

# colnames(quant_contigs) <- colnames(quant_contigs) %>% 
#   str_replace_all("//.", "//_")

rownames(quant_contigs) <- paste0("h1_",
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
<div id="htmlwidget-6908e02842004b8bd79f" style="width:100%;height:auto;" class="datatables html-widget"></div>
<script type="application/json" data-for="htmlwidget-6908e02842004b8bd79f">{"x":{"filter":"none","vertical":false,"data":[["1","2","3","4","5","6"],["h1_megahit_1","h1_megahit_2","h1_megahit_3","h1_megahit_4","h1_megahit_5","h1_megahit_6"],[1.5,0,0,0,0,0],[0,1.5,0,0,0,0],[0,0,1.7,0,1.7,0],[0,0,0,0,1.4,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0.7,0],[0,1.6,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,2.8,0,0],[0,0,1.1,0,0,0],[0,0,0,0,0.6,0],[0,0,0,0,0,0],[0,0,1.4,0,0,0],[0,0,0.9,0,0,0],[0,0,0,0,0,0],[0,0,1,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0.8,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,2.1,0,0],[0,0,0,0,0,0],[0,0,1.6,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,2.2,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[3.83,0,0,0,0,0],[0,0.98,0,0,0,0],[0,0,1.36,0,0.88,0],[0,0,0,0,0.88,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0.13,0],[0,0.98,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,1.85,0,0],[0,0,0.68,0,0,0],[0,0,0,0,0.21,0],[0,0,0,0,0,0],[0,0,0.68,0,0,0],[0,0,0.68,0,0,0],[0,0,0,0,0,0],[0,0,0.68,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0.93,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0.92,0,0],[0,0,0,0,0,0],[0,0,0.92,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,2.04,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[8,0,0,0,0,0],[0,2,0,0,0,0],[0,0,2,0,3,0],[0,0,0,0,2,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,1,0],[0,2,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,4,0,0],[0,0,1,0,0,0],[0,0,0,0,1,0],[0,0,0,0,0,0],[0,0,1,0,0,0],[0,0,1,0,0,0],[0,0,0,0,0,0],[0,0,1,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,2,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,2,0,0],[0,0,0,0,0,0],[0,0,2,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,3,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0]],"container":"<table class=\"display\">\n  <thead>\n    <tr>\n      <th> <\/th>\n      <th>Contig_ID<\/th>\n      <th>TPM.D.1<\/th>\n      <th>TPM.HC.1<\/th>\n      <th>TPM.HC.2<\/th>\n      <th>TPM.HC.3<\/th>\n      <th>TPM.HC.4<\/th>\n      <th>TPM.HC.5<\/th>\n      <th>TPM.HC.6<\/th>\n      <th>TPM.HC.7<\/th>\n      <th>TPM.HC.8<\/th>\n      <th>TPM.HC.9<\/th>\n      <th>TPM.HC.10<\/th>\n      <th>TPM.HC.11<\/th>\n      <th>TPM.HC.12<\/th>\n      <th>TPM.HC.13<\/th>\n      <th>TPM.HC.14<\/th>\n      <th>TPM.HC.15<\/th>\n      <th>TPM.HC.16<\/th>\n      <th>TPM.HC.17<\/th>\n      <th>TPM.HC.18<\/th>\n      <th>TPM.HC.19<\/th>\n      <th>TPM.HC.20<\/th>\n      <th>TPM.HC.21<\/th>\n      <th>TPM.HC.22<\/th>\n      <th>TPM.HC.23<\/th>\n      <th>TPM.HC.24<\/th>\n      <th>TPM.HV.1<\/th>\n      <th>TPM.HV.2<\/th>\n      <th>TPM.HV.3<\/th>\n      <th>TPM.HV.4<\/th>\n      <th>TPM.HV.5<\/th>\n      <th>TPM.HV.6<\/th>\n      <th>TPM.HV.7<\/th>\n      <th>TPM.HV.8<\/th>\n      <th>TPM.HV.9<\/th>\n      <th>TPM.HV.10<\/th>\n      <th>TPM.HV.11<\/th>\n      <th>TPM.HV.12<\/th>\n      <th>TPM.HV.13<\/th>\n      <th>TPM.HV.14<\/th>\n      <th>TPM.HV.15<\/th>\n      <th>TPM.HV.16<\/th>\n      <th>TPM.HV.17<\/th>\n      <th>TPM.HV.18<\/th>\n      <th>TPM.HV.19<\/th>\n      <th>TPM.HV.20<\/th>\n      <th>TPM.HV.21<\/th>\n      <th>TPM.HV.22<\/th>\n      <th>TPM.HV.23<\/th>\n      <th>TPM.HV.24<\/th>\n      <th>Coverage.D.1<\/th>\n      <th>Coverage.HC.1<\/th>\n      <th>Coverage.HC.2<\/th>\n      <th>Coverage.HC.3<\/th>\n      <th>Coverage.HC.4<\/th>\n      <th>Coverage.HC.5<\/th>\n      <th>Coverage.HC.6<\/th>\n      <th>Coverage.HC.7<\/th>\n      <th>Coverage.HC.8<\/th>\n      <th>Coverage.HC.9<\/th>\n      <th>Coverage.HC.10<\/th>\n      <th>Coverage.HC.11<\/th>\n      <th>Coverage.HC.12<\/th>\n      <th>Coverage.HC.13<\/th>\n      <th>Coverage.HC.14<\/th>\n      <th>Coverage.HC.15<\/th>\n      <th>Coverage.HC.16<\/th>\n      <th>Coverage.HC.17<\/th>\n      <th>Coverage.HC.18<\/th>\n      <th>Coverage.HC.19<\/th>\n      <th>Coverage.HC.20<\/th>\n      <th>Coverage.HC.21<\/th>\n      <th>Coverage.HC.22<\/th>\n      <th>Coverage.HC.23<\/th>\n      <th>Coverage.HC.24<\/th>\n      <th>Coverage.HV.1<\/th>\n      <th>Coverage.HV.2<\/th>\n      <th>Coverage.HV.3<\/th>\n      <th>Coverage.HV.4<\/th>\n      <th>Coverage.HV.5<\/th>\n      <th>Coverage.HV.6<\/th>\n      <th>Coverage.HV.7<\/th>\n      <th>Coverage.HV.8<\/th>\n      <th>Coverage.HV.9<\/th>\n      <th>Coverage.HV.10<\/th>\n      <th>Coverage.HV.11<\/th>\n      <th>Coverage.HV.12<\/th>\n      <th>Coverage.HV.13<\/th>\n      <th>Coverage.HV.14<\/th>\n      <th>Coverage.HV.15<\/th>\n      <th>Coverage.HV.16<\/th>\n      <th>Coverage.HV.17<\/th>\n      <th>Coverage.HV.18<\/th>\n      <th>Coverage.HV.19<\/th>\n      <th>Coverage.HV.20<\/th>\n      <th>Coverage.HV.21<\/th>\n      <th>Coverage.HV.22<\/th>\n      <th>Coverage.HV.23<\/th>\n      <th>Coverage.HV.24<\/th>\n      <th>Raw.read.count.D.1<\/th>\n      <th>Raw.read.count.HC.1<\/th>\n      <th>Raw.read.count.HC.2<\/th>\n      <th>Raw.read.count.HC.3<\/th>\n      <th>Raw.read.count.HC.4<\/th>\n      <th>Raw.read.count.HC.5<\/th>\n      <th>Raw.read.count.HC.6<\/th>\n      <th>Raw.read.count.HC.7<\/th>\n      <th>Raw.read.count.HC.8<\/th>\n      <th>Raw.read.count.HC.9<\/th>\n      <th>Raw.read.count.HC.10<\/th>\n      <th>Raw.read.count.HC.11<\/th>\n      <th>Raw.read.count.HC.12<\/th>\n      <th>Raw.read.count.HC.13<\/th>\n      <th>Raw.read.count.HC.14<\/th>\n      <th>Raw.read.count.HC.15<\/th>\n      <th>Raw.read.count.HC.16<\/th>\n      <th>Raw.read.count.HC.17<\/th>\n      <th>Raw.read.count.HC.18<\/th>\n      <th>Raw.read.count.HC.19<\/th>\n      <th>Raw.read.count.HC.20<\/th>\n      <th>Raw.read.count.HC.21<\/th>\n      <th>Raw.read.count.HC.22<\/th>\n      <th>Raw.read.count.HC.23<\/th>\n      <th>Raw.read.count.HC.24<\/th>\n      <th>Raw.read.count.HV.1<\/th>\n      <th>Raw.read.count.HV.2<\/th>\n      <th>Raw.read.count.HV.3<\/th>\n      <th>Raw.read.count.HV.4<\/th>\n      <th>Raw.read.count.HV.5<\/th>\n      <th>Raw.read.count.HV.6<\/th>\n      <th>Raw.read.count.HV.7<\/th>\n      <th>Raw.read.count.HV.8<\/th>\n      <th>Raw.read.count.HV.9<\/th>\n      <th>Raw.read.count.HV.10<\/th>\n      <th>Raw.read.count.HV.11<\/th>\n      <th>Raw.read.count.HV.12<\/th>\n      <th>Raw.read.count.HV.13<\/th>\n      <th>Raw.read.count.HV.14<\/th>\n      <th>Raw.read.count.HV.15<\/th>\n      <th>Raw.read.count.HV.16<\/th>\n      <th>Raw.read.count.HV.17<\/th>\n      <th>Raw.read.count.HV.18<\/th>\n      <th>Raw.read.count.HV.19<\/th>\n      <th>Raw.read.count.HV.20<\/th>\n      <th>Raw.read.count.HV.21<\/th>\n      <th>Raw.read.count.HV.22<\/th>\n      <th>Raw.read.count.HV.23<\/th>\n      <th>Raw.read.count.HV.24<\/th>\n    <\/tr>\n  <\/thead>\n<\/table>","options":{"columnDefs":[{"className":"dt-right","targets":[2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58,59,60,61,62,63,64,65,66,67,68,69,70,71,72,73,74,75,76,77,78,79,80,81,82,83,84,85,86,87,88,89,90,91,92,93,94,95,96,97,98,99,100,101,102,103,104,105,106,107,108,109,110,111,112,113,114,115,116,117,118,119,120,121,122,123,124,125,126,127,128,129,130,131,132,133,134,135,136,137,138,139,140,141,142,143,144,145,146,147,148]},{"orderable":false,"targets":0}],"order":[],"autoWidth":false,"orderClasses":false}},"evals":[],"jsHooks":[]}</script>
```

Combine quantification data with annotations:


```r
colnames(quant_orfs) <- paste0("ORF_" ,colnames(quant_orfs))

colnames(quant_contigs) <- paste0("contig_" ,colnames(quant_contigs))
```



```r
full_gene_catalogue %>% 
  left_join(quant_orfs,
            by = c("ORF_ID" = "ORF_ORF_ID")) -> full_gene_catalogue_orf_quant

full_gene_catalogue_orf_quant %>% 
  left_join(quant_contigs,
            by = c("Contig_ID" = "contig_Contig_ID")) -> full_quant_full_gene_catalogue
```


```r
full_gene_catalogue_orf_quant %>% 
  dplyr::select(ORF_ID, SQM_Length.NT, starts_with("ORF_Raw.read.count.")) %>% 
  pivot_longer(cols = starts_with("ORF_Raw.read.count."),
               values_to = "read_counts", 
               names_to = "sample") %>% # -> test  test %>%  group_by(sample, ORF_ID) %>%  dplyr::summarise(count = n()) -> test test %>%  arrange(count)
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
<div id="htmlwidget-bcf80cde3ed5fe2a5ede" style="width:100%;height:auto;" class="datatables html-widget"></div>
<script type="application/json" data-for="htmlwidget-bcf80cde3ed5fe2a5ede">{"x":{"filter":"none","vertical":false,"data":[["1","2","3","4","5","6"],["h1_megahit_1_1-312","h1_megahit_2_3-176","h1_megahit_2_173-304","h1_megahit_3_3-116","h1_megahit_3_82-219","h1_megahit_4_1-72"],[0.0224358974358974,0,0,0,0,0],[0,0.0114942528735632,0.00757575757575758,0,0,0],[0,0,0,0.0175438596491228,0.0144927536231884,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0.0114942528735632,0.0151515151515152,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0.0277777777777778],[0,0,0,0.0087719298245614,0.0072463768115942,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0.0087719298245614,0.0072463768115942,0],[0,0,0,0.0087719298245614,0.0072463768115942,0],[0,0,0,0,0,0],[0,0,0,0.0087719298245614,0.0072463768115942,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0.0277777777777778],[0,0,0,0,0,0],[0,0,0,0.0087719298245614,0.0072463768115942,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0.0263157894736842,0.0217391304347826,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0]],"container":"<table class=\"display\">\n  <thead>\n    <tr>\n      <th> <\/th>\n      <th>ORF_ID<\/th>\n      <th>ORF_Num_Gi.D.1<\/th>\n      <th>ORF_Num_Gi.HC.1<\/th>\n      <th>ORF_Num_Gi.HC.2<\/th>\n      <th>ORF_Num_Gi.HC.3<\/th>\n      <th>ORF_Num_Gi.HC.4<\/th>\n      <th>ORF_Num_Gi.HC.5<\/th>\n      <th>ORF_Num_Gi.HC.6<\/th>\n      <th>ORF_Num_Gi.HC.7<\/th>\n      <th>ORF_Num_Gi.HC.8<\/th>\n      <th>ORF_Num_Gi.HC.9<\/th>\n      <th>ORF_Num_Gi.HC.10<\/th>\n      <th>ORF_Num_Gi.HC.11<\/th>\n      <th>ORF_Num_Gi.HC.12<\/th>\n      <th>ORF_Num_Gi.HC.13<\/th>\n      <th>ORF_Num_Gi.HC.14<\/th>\n      <th>ORF_Num_Gi.HC.15<\/th>\n      <th>ORF_Num_Gi.HC.16<\/th>\n      <th>ORF_Num_Gi.HC.17<\/th>\n      <th>ORF_Num_Gi.HC.18<\/th>\n      <th>ORF_Num_Gi.HC.19<\/th>\n      <th>ORF_Num_Gi.HC.20<\/th>\n      <th>ORF_Num_Gi.HC.21<\/th>\n      <th>ORF_Num_Gi.HC.22<\/th>\n      <th>ORF_Num_Gi.HC.23<\/th>\n      <th>ORF_Num_Gi.HC.24<\/th>\n      <th>ORF_Num_Gi.HV.1<\/th>\n      <th>ORF_Num_Gi.HV.2<\/th>\n      <th>ORF_Num_Gi.HV.3<\/th>\n      <th>ORF_Num_Gi.HV.4<\/th>\n      <th>ORF_Num_Gi.HV.5<\/th>\n      <th>ORF_Num_Gi.HV.6<\/th>\n      <th>ORF_Num_Gi.HV.7<\/th>\n      <th>ORF_Num_Gi.HV.8<\/th>\n      <th>ORF_Num_Gi.HV.9<\/th>\n      <th>ORF_Num_Gi.HV.10<\/th>\n      <th>ORF_Num_Gi.HV.11<\/th>\n      <th>ORF_Num_Gi.HV.12<\/th>\n      <th>ORF_Num_Gi.HV.13<\/th>\n      <th>ORF_Num_Gi.HV.14<\/th>\n      <th>ORF_Num_Gi.HV.15<\/th>\n      <th>ORF_Num_Gi.HV.16<\/th>\n      <th>ORF_Num_Gi.HV.17<\/th>\n      <th>ORF_Num_Gi.HV.18<\/th>\n      <th>ORF_Num_Gi.HV.19<\/th>\n      <th>ORF_Num_Gi.HV.20<\/th>\n      <th>ORF_Num_Gi.HV.21<\/th>\n      <th>ORF_Num_Gi.HV.22<\/th>\n      <th>ORF_Num_Gi.HV.23<\/th>\n      <th>ORF_Num_Gi.HV.24<\/th>\n    <\/tr>\n  <\/thead>\n<\/table>","options":{"columnDefs":[{"className":"dt-right","targets":[2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50]},{"orderable":false,"targets":0}],"order":[],"autoWidth":false,"orderClasses":false}},"evals":[],"jsHooks":[]}</script>
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
<div id="htmlwidget-e7d0dc9929da45f03498" style="width:100%;height:auto;" class="datatables html-widget"></div>
<script type="application/json" data-for="htmlwidget-e7d0dc9929da45f03498">{"x":{"filter":"none","vertical":false,"data":[["1","2","3","4","5","6"],["h1_megahit_1_1-312","h1_megahit_2_3-176","h1_megahit_2_173-304","h1_megahit_3_3-116","h1_megahit_3_82-219","h1_megahit_4_1-72"],[1.65106252977147e-07,0,0,0,0,0],[0,2.28833590277493e-07,1.50822139046529e-07,0,0,0],[0,0,0,2.97590210493958e-07,2.45835391277617e-07,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,2.30437459025747e-07,3.03758468715757e-07,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,4.04782802479438e-07],[0,0,0,1.56397247893644e-07,1.29197726520837e-07,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,1.78217802413832e-07,1.47223401994035e-07,0],[0,0,0,1.40701542858223e-07,1.16231709317663e-07,0],[0,0,0,0,0,0],[0,0,0,1.5402926080401e-07,1.27241563272878e-07,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,5.52850535761806e-07],[0,0,0,0,0,0],[0,0,0,1.65488750877943e-07,1.36708098551344e-07,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,4.33713887036034e-07,3.58285384942811e-07,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0]],"container":"<table class=\"display\">\n  <thead>\n    <tr>\n      <th> <\/th>\n      <th>ORF_ID<\/th>\n      <th>ORF_Num_Gi_pc.D.1<\/th>\n      <th>ORF_Num_Gi_pc.HC.1<\/th>\n      <th>ORF_Num_Gi_pc.HC.2<\/th>\n      <th>ORF_Num_Gi_pc.HC.3<\/th>\n      <th>ORF_Num_Gi_pc.HC.4<\/th>\n      <th>ORF_Num_Gi_pc.HC.5<\/th>\n      <th>ORF_Num_Gi_pc.HC.6<\/th>\n      <th>ORF_Num_Gi_pc.HC.7<\/th>\n      <th>ORF_Num_Gi_pc.HC.8<\/th>\n      <th>ORF_Num_Gi_pc.HC.9<\/th>\n      <th>ORF_Num_Gi_pc.HC.10<\/th>\n      <th>ORF_Num_Gi_pc.HC.11<\/th>\n      <th>ORF_Num_Gi_pc.HC.12<\/th>\n      <th>ORF_Num_Gi_pc.HC.13<\/th>\n      <th>ORF_Num_Gi_pc.HC.14<\/th>\n      <th>ORF_Num_Gi_pc.HC.15<\/th>\n      <th>ORF_Num_Gi_pc.HC.16<\/th>\n      <th>ORF_Num_Gi_pc.HC.17<\/th>\n      <th>ORF_Num_Gi_pc.HC.18<\/th>\n      <th>ORF_Num_Gi_pc.HC.19<\/th>\n      <th>ORF_Num_Gi_pc.HC.20<\/th>\n      <th>ORF_Num_Gi_pc.HC.21<\/th>\n      <th>ORF_Num_Gi_pc.HC.22<\/th>\n      <th>ORF_Num_Gi_pc.HC.23<\/th>\n      <th>ORF_Num_Gi_pc.HC.24<\/th>\n      <th>ORF_Num_Gi_pc.HV.1<\/th>\n      <th>ORF_Num_Gi_pc.HV.2<\/th>\n      <th>ORF_Num_Gi_pc.HV.3<\/th>\n      <th>ORF_Num_Gi_pc.HV.4<\/th>\n      <th>ORF_Num_Gi_pc.HV.5<\/th>\n      <th>ORF_Num_Gi_pc.HV.6<\/th>\n      <th>ORF_Num_Gi_pc.HV.7<\/th>\n      <th>ORF_Num_Gi_pc.HV.8<\/th>\n      <th>ORF_Num_Gi_pc.HV.9<\/th>\n      <th>ORF_Num_Gi_pc.HV.10<\/th>\n      <th>ORF_Num_Gi_pc.HV.11<\/th>\n      <th>ORF_Num_Gi_pc.HV.12<\/th>\n      <th>ORF_Num_Gi_pc.HV.13<\/th>\n      <th>ORF_Num_Gi_pc.HV.14<\/th>\n      <th>ORF_Num_Gi_pc.HV.15<\/th>\n      <th>ORF_Num_Gi_pc.HV.16<\/th>\n      <th>ORF_Num_Gi_pc.HV.17<\/th>\n      <th>ORF_Num_Gi_pc.HV.18<\/th>\n      <th>ORF_Num_Gi_pc.HV.19<\/th>\n      <th>ORF_Num_Gi_pc.HV.20<\/th>\n      <th>ORF_Num_Gi_pc.HV.21<\/th>\n      <th>ORF_Num_Gi_pc.HV.22<\/th>\n      <th>ORF_Num_Gi_pc.HV.23<\/th>\n      <th>ORF_Num_Gi_pc.HV.24<\/th>\n    <\/tr>\n  <\/thead>\n<\/table>","options":{"columnDefs":[{"className":"dt-right","targets":[2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50]},{"orderable":false,"targets":0}],"order":[],"autoWidth":false,"orderClasses":false}},"evals":[],"jsHooks":[]}</script>
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
human_meta %>% 
  select(-sample:-index2) -> human_meta
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
## [1] 791180     49
```

```r
physeq_tpm <- phyloseq(tpm,
                       tax,
                       human_meta %>% sample_data())


physeq_tpm
```

```
## phyloseq-class experiment-level object
## otu_table()   OTU Table:         [ 791180 taxa and 49 samples ]
## sample_data() Sample Data:       [ 49 samples by 39 sample variables ]
## tax_table()   Taxonomy Table:    [ 791180 taxa by 97 taxonomic ranks ]
```


```r
dim(cov)
```

```
## [1] 791180     49
```

```r
physeq_cov <- phyloseq(cov,
                       tax,
                       human_meta %>% sample_data())

physeq_cov
```

```
## phyloseq-class experiment-level object
## otu_table()   OTU Table:         [ 791180 taxa and 49 samples ]
## sample_data() Sample Data:       [ 49 samples by 39 sample variables ]
## tax_table()   Taxonomy Table:    [ 791180 taxa by 97 taxonomic ranks ]
```



```r
dim(Num_Gi_pc)
```

```
## [1] 791180     49
```

```r
physeq_Num_Gi_pc <- phyloseq(Num_Gi_pc,
                             tax,
                             human_meta %>% sample_data())

physeq_Num_Gi_pc
```

```
## phyloseq-class experiment-level object
## otu_table()   OTU Table:         [ 791180 taxa and 49 samples ]
## sample_data() Sample Data:       [ 49 samples by 39 sample variables ]
## tax_table()   Taxonomy Table:    [ 791180 taxa by 97 taxonomic ranks ]
```


```r
dim(Num_Gi)
```

```
## [1] 791180     49
```

```r
physeq_Num_Gi<- phyloseq(Num_Gi,
                         tax,
                         human_meta %>% sample_data())

physeq_Num_Gi
```

```
## phyloseq-class experiment-level object
## otu_table()   OTU Table:         [ 791180 taxa and 49 samples ]
## sample_data() Sample Data:       [ 49 samples by 39 sample variables ]
## tax_table()   Taxonomy Table:    [ 791180 taxa by 97 taxonomic ranks ]
```



```r
dim(read_count)
```

```
## [1] 791180     49
```

```r
physeq_read_count <- phyloseq(read_count,
                              tax,
                              human_meta %>% sample_data())

physeq_read_count
```

```
## phyloseq-class experiment-level object
## otu_table()   OTU Table:         [ 791180 taxa and 49 samples ]
## sample_data() Sample Data:       [ 49 samples by 39 sample variables ]
## tax_table()   Taxonomy Table:    [ 791180 taxa by 97 taxonomic ranks ]
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
  saveRDS(here::here("data/processed/human1_full_gene_catalog_phyloseq.RDS"))

full_gene_catalogue_orf_quant_all_metrics %>% 
  saveRDS(here::here("data/processed/human1_full_gene_catalog_full_metrics_table.RDS"))

full_gene_catalogue_orf_quant_all_metrics %>% 
  write_tsv(here::here("data/processed/human1_full_gene_catalog_full_metrics_table.tsv.gz"))
```


Export only AMR and Tox:


```r
full_gene_catalogue_orf_quant_all_metrics %>% 
  filter(!is.na(PathoFact_AMR_ARG) | !is.na(rgi_CARD_ARO) | !is.na(PathoFact_Tox_Toxin_classification)) %>% 
  write_tsv(here::here("data/processed/human1_full_gene_catalog_full_metrics_table_AMR_tox.tsv"))
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

