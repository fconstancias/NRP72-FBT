---
title: "Merge human gene catalogue + quantification"
author: "Hannah Li Hägi & Florentin Constancias "
date: " June 27, 2022 "
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
source("https://raw.githubusercontent.com/fconstancias/DivComAnalyses/master/R/phyloseq_normalisation.R")
```

```
## Loading required package: scales
```

```
## 
## Attaching package: 'scales'
```

```
## The following object is masked from 'package:purrr':
## 
##     discard
```

```
## The following object is masked from 'package:readr':
## 
##     col_factor
```

```
## Loading required package: reshape2
```

```
## 
## Attaching package: 'reshape2'
```

```
## The following object is masked from 'package:tidyr':
## 
##     smiths
```

```r
source("https://raw.githubusercontent.com/fconstancias/DivComAnalyses/master/R/phyloseq_alpha.R")
source("https://raw.githubusercontent.com/fconstancias/DivComAnalyses/master/R/phyloseq_beta.R")
source("https://raw.githubusercontent.com/fconstancias/DivComAnalyses/master/R/phyloseq_heatmap.R")
```

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
  na_if("n/a") -> PathoFact_Tox
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
colnames(PathoFact_Tox) <- paste0("PathoFact_Tox_", colnames(PathoFact_Tox))

PathoFact_Tox %>% 
  rename(PathoFact_Tox_Toxin_classification = PathoFact_Tox_Toxin_confidence_level,
         ORF_ID = PathoFact_Tox_ORF) -> PathoFact_Tox
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
<div id="htmlwidget-ed124481f0e82d66595a" style="width:100%;height:auto;" class="datatables html-widget"></div>
<script type="application/json" data-for="htmlwidget-ed124481f0e82d66595a">{"x":{"filter":"none","vertical":false,"data":[["1","2","3","4","5","6"],["megahit_1512_86-616","megahit_1512_743-1915","megahit_1512_1932-3470","megahit_1935_13232-13930","megahit_4077_533-1408","megahit_4631_5392-6930"],["Perfect","Strict","Perfect","Strict","Perfect","Strict"],[280,675,900,400,500,900],[349,739.6,996.1,458,558.1,963.4],["emrR","emrA","emrB","vanRD","CTX-M-1","emrB"],[100,99.74,100,96.12,100,95.51],[3000516,3000027,3000074,3002923,3001864,3000074],["protein homolog model","protein homolog model","protein homolog model","protein homolog model","protein homolog model","protein homolog model"],[null,null,null,null,null,null],[null,null,null,null,null,null],["fluoroquinolone antibiotic","fluoroquinolone antibiotic","fluoroquinolone antibiotic","glycopeptide antibiotic","cephalosporin","fluoroquinolone antibiotic"],["antibiotic efflux","antibiotic efflux","antibiotic efflux","antibiotic target alteration","antibiotic inactivation","antibiotic efflux"],["major facilitator superfamily (MFS) antibiotic efflux pump","major facilitator superfamily (MFS) antibiotic efflux pump","major facilitator superfamily (MFS) antibiotic efflux pump","glycopeptide resistance gene cluster; vanR","CTX-M beta-lactamase","major facilitator superfamily (MFS) antibiotic efflux pump"],[100,100,100,100,100,100],["gnl|BL_ORD_ID|1254|hsp_num:0","gnl|BL_ORD_ID|1657|hsp_num:0","gnl|BL_ORD_ID|1746|hsp_num:0","gnl|BL_ORD_ID|950|hsp_num:0","gnl|BL_ORD_ID|392|hsp_num:0","gnl|BL_ORD_ID|1746|hsp_num:0"],[1330,1757,1847,1004,420,1847],[null,null,null,null,null,null],[null,null,null,null,null,null]],"container":"<table class=\"display\">\n  <thead>\n    <tr>\n      <th> <\/th>\n      <th>ORF_ID<\/th>\n      <th>rgi_CARD_Cut_Off<\/th>\n      <th>rgi_CARD_Pass_Bitscore<\/th>\n      <th>rgi_CARD_Best_Hit_Bitscore<\/th>\n      <th>rgi_CARD_Best_Hit_ARO<\/th>\n      <th>rgi_CARD_Best_Identities<\/th>\n      <th>rgi_CARD_ARO<\/th>\n      <th>rgi_CARD_Model_type<\/th>\n      <th>rgi_CARD_SNPs_in_Best_Hit_ARO<\/th>\n      <th>rgi_CARD_Other_SNPs<\/th>\n      <th>rgi_CARD_Drug Class<\/th>\n      <th>rgi_CARD_Resistance Mechanism<\/th>\n      <th>rgi_CARD_AMR Gene Family<\/th>\n      <th>rgi_CARD_Percentage Length of Reference Sequence<\/th>\n      <th>rgi_CARD_ID<\/th>\n      <th>rgi_CARD_Model_ID<\/th>\n      <th>rgi_CARD_Nudged<\/th>\n      <th>rgi_CARD_Note<\/th>\n    <\/tr>\n  <\/thead>\n<\/table>","options":{"columnDefs":[{"className":"dt-right","targets":[3,4,6,7,14,16]},{"orderable":false,"targets":0}],"order":[],"autoWidth":false,"orderClasses":false}},"evals":[],"jsHooks":[]}</script>
```

## Waffle HGT: 



```r
# TODO move to more appropriate section

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
<div id="htmlwidget-f75c638e93d792a7d8dd" style="width:100%;height:auto;" class="datatables html-widget"></div>
<script type="application/json" data-for="htmlwidget-f75c638e93d792a7d8dd">{"x":{"filter":"none","vertical":false,"data":[["1","2","3","4","5","6"],["megahit_100032","megahit_10009","megahit_100143","megahit_100161","megahit_100166","megahit_100186"],["IS30","ISNCY","IS66","IS200/IS605","IS30","IS66"],["IS30_241","ISNCY_86","IS66_46","IS200/IS605_96","IS30_223","IS66_43"],["p","p","p","p","c","p"],[1,1,1,1,1,2]],"container":"<table class=\"display\">\n  <thead>\n    <tr>\n      <th> <\/th>\n      <th>Contig_ID<\/th>\n      <th>isescan_family<\/th>\n      <th>isescan_cluster<\/th>\n      <th>isescan_type<\/th>\n      <th>isescan_count_IS_contig<\/th>\n    <\/tr>\n  <\/thead>\n<\/table>","options":{"columnDefs":[{"className":"dt-right","targets":5},{"orderable":false,"targets":0}],"order":[],"autoWidth":false,"orderClasses":false}},"evals":[],"jsHooks":[]}</script>
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
  select(-TPM.D.1:-Hits) -> orf

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
<div id="htmlwidget-023fa45fcd5a2088887d" style="width:100%;height:auto;" class="datatables html-widget"></div>
<script type="application/json" data-for="htmlwidget-023fa45fcd5a2088887d">{"x":{"filter":"none","vertical":false,"data":[["1","2","3","4","5","6"],["megahit_1_1-312","megahit_2_3-176","megahit_2_173-304","megahit_3_3-116","megahit_3_82-219","megahit_4_1-72"],["megahit_1","megahit_2","megahit_2","megahit_3","megahit_3","megahit_4"],[312,174,132,114,138,72],[104,58,44,38,46,24],[34.62,38.51,40.91,51.75,47.83,36.11],["RP-S5, MRPS5, rpsE",null,null,null,null,null],["k_Bacteria;n_Terrabacteria group;p_Firmicutes;c_Erysipelotrichia;o_Erysipelotrichales;f_Erysipelotrichaceae;g_Massilimicrobiota","k_Bacteria;n_Terrabacteria group;p_Firmicutes;c_Clostridia;o_Eubacteriales;n_unclassified Eubacteriales;s_Clostridiales bacterium","k_Bacteria;n_Terrabacteria group;p_Firmicutes","k_Bacteria;n_Terrabacteria group;p_Firmicutes","k_Bacteria;n_Terrabacteria group;p_Firmicutes;c_Clostridia;o_Eubacteriales","k_Bacteria;n_Terrabacteria group;p_Firmicutes;c_Bacilli;o_Lactobacillales;f_Enterococcaceae;g_Enterococcus;s_Enterococcus faecium"],["K02988*",null,null,null,null,null],["small subunit ribosomal protein S5",null,null,null,null,null],["Genetic Information Processing; Translation; Ribosome | Brite Hierarchies; Protein families: genetic information processing; Ribosome",null,null,null,null,null],["COG0098*",null,null,"ENOG410ZP7K*","COG3843*",null],["Ribosomal protein S5",null,null,"D repair protein RadA domain","Type IV secretory pathway, VirD2 components (relaxase)",null],["Translation, ribosomal structure and biogenesis",null,null,"Function unknown",null,null],["PF03719 [Ribosomal protein S5, C-terminal domain]",null,null,null,null,null]],"container":"<table class=\"display\">\n  <thead>\n    <tr>\n      <th> <\/th>\n      <th>ORF_ID<\/th>\n      <th>Contig_ID<\/th>\n      <th>SQM_Length.NT<\/th>\n      <th>SQM_Length.AA<\/th>\n      <th>SQM_GC.perc<\/th>\n      <th>SQM_Gene.name<\/th>\n      <th>SQM_Tax_orf<\/th>\n      <th>SQM_KEGG.ID<\/th>\n      <th>SQM_KEGGFUN<\/th>\n      <th>SQM_KEGGPATH<\/th>\n      <th>SQM_COG.ID<\/th>\n      <th>SQM_COGFUN<\/th>\n      <th>SQM_COGPATH<\/th>\n      <th>SQM_PFAM<\/th>\n    <\/tr>\n  <\/thead>\n<\/table>","options":{"columnDefs":[{"className":"dt-right","targets":[3,4,5]},{"orderable":false,"targets":0}],"order":[],"autoWidth":false,"orderClasses":false}},"evals":[],"jsHooks":[]}</script>
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
  select(-Coverage.D.1:-Raw.read.count.HV.24, -Bin.ID) -> contigs


colnames(contigs) <- paste0("SQM_", colnames(contigs))

contigs %>% 
  rename(Contig_ID = SQM_Contig_ID) -> contigs

contigs %>% 
  head() %>% 
  DT::datatable()
```

```{=html}
<div id="htmlwidget-06564437a2a15aead1e3" style="width:100%;height:auto;" class="datatables html-widget"></div>
<script type="application/json" data-for="htmlwidget-06564437a2a15aead1e3">{"x":{"filter":"none","vertical":false,"data":[["1","2","3","4","5","6"],["megahit_1","megahit_2","megahit_3","megahit_4","megahit_5","megahit_6"],["k_Bacteria;n_Terrabacteria group;p_Firmicutes;c_Erysipelotrichia;o_Erysipelotrichales;f_Erysipelotrichaceae;g_Massilimicrobiota","k_Bacteria;n_Terrabacteria group;p_Firmicutes;c_Clostridia;o_Eubacteriales;n_unclassified Eubacteriales;s_Clostridiales bacterium","k_Bacteria;n_Terrabacteria group;p_Firmicutes;c_Clostridia;o_Eubacteriales","k_Bacteria;n_Terrabacteria group;p_Firmicutes;c_Bacilli;o_Lactobacillales;f_Enterococcaceae;g_Enterococcus;s_Enterococcus faecium","","k_Bacteria;n_Terrabacteria group;p_Firmicutes;c_Clostridia;o_Eubacteriales"],[0,0,0,0,0,0],[34.5,39.67,51.13,38.89,50.74,45.22],[313,305,221,324,339,230],[1,2,2,2,1,1]],"container":"<table class=\"display\">\n  <thead>\n    <tr>\n      <th> <\/th>\n      <th>Contig_ID<\/th>\n      <th>SQM_Tax_contig<\/th>\n      <th>SQM_Disparity_contig<\/th>\n      <th>SQM_GC_perc_contig<\/th>\n      <th>SQM_Length_contig<\/th>\n      <th>SQM_Num_genes_contig<\/th>\n    <\/tr>\n  <\/thead>\n<\/table>","options":{"columnDefs":[{"className":"dt-right","targets":[3,4,5,6]},{"orderable":false,"targets":0}],"order":[],"autoWidth":false,"orderClasses":false}},"evals":[],"jsHooks":[]}</script>
```



```r
orf %>% 
  left_join(card,
            by = c("ORF_ID" = "ORF_ID")) %>% 
  left_join(PathoFact_AMR) %>% 
  left_join(PathoFact_Tox) %>% 
  left_join(viralverify,
            by = c("Contig_ID" = "Contig_ID")) %>% 
  left_join(top_isescan) %>% 
  left_join(contigs,
            by = c("Contig_ID" = "Contig_ID")) -> orf_contig
```

```
## Joining, by = "ORF_ID"
## Joining, by = "ORF_ID"
## Joining, by = "Contig_ID"
```

```r
orf_contig %>% 
  head() %>% 
  DT::datatable()
```

```{=html}
<div id="htmlwidget-a5b572cbfa9f8f737332" style="width:100%;height:auto;" class="datatables html-widget"></div>
<script type="application/json" data-for="htmlwidget-a5b572cbfa9f8f737332">{"x":{"filter":"none","vertical":false,"data":[["1","2","3","4","5","6"],["megahit_1_1-312","megahit_2_3-176","megahit_2_173-304","megahit_3_3-116","megahit_3_82-219","megahit_4_1-72"],["megahit_1","megahit_2","megahit_2","megahit_3","megahit_3","megahit_4"],[312,174,132,114,138,72],[104,58,44,38,46,24],[34.62,38.51,40.91,51.75,47.83,36.11],["RP-S5, MRPS5, rpsE",null,null,null,null,null],["k_Bacteria;n_Terrabacteria group;p_Firmicutes;c_Erysipelotrichia;o_Erysipelotrichales;f_Erysipelotrichaceae;g_Massilimicrobiota","k_Bacteria;n_Terrabacteria group;p_Firmicutes;c_Clostridia;o_Eubacteriales;n_unclassified Eubacteriales;s_Clostridiales bacterium","k_Bacteria;n_Terrabacteria group;p_Firmicutes","k_Bacteria;n_Terrabacteria group;p_Firmicutes","k_Bacteria;n_Terrabacteria group;p_Firmicutes;c_Clostridia;o_Eubacteriales","k_Bacteria;n_Terrabacteria group;p_Firmicutes;c_Bacilli;o_Lactobacillales;f_Enterococcaceae;g_Enterococcus;s_Enterococcus faecium"],["K02988*",null,null,null,null,null],["small subunit ribosomal protein S5",null,null,null,null,null],["Genetic Information Processing; Translation; Ribosome | Brite Hierarchies; Protein families: genetic information processing; Ribosome",null,null,null,null,null],["COG0098*",null,null,"ENOG410ZP7K*","COG3843*",null],["Ribosomal protein S5",null,null,"D repair protein RadA domain","Type IV secretory pathway, VirD2 components (relaxase)",null],["Translation, ribosomal structure and biogenesis",null,null,"Function unknown",null,null],["PF03719 [Ribosomal protein S5, C-terminal domain]",null,null,null,null,null],[null,null,null,null,null,null],[null,null,null,null,null,null],[null,null,null,null,null,null],[null,null,null,null,null,null],[null,null,null,null,null,null],[null,null,null,null,null,null],[null,null,null,null,null,null],[null,null,null,null,null,null],[null,null,null,null,null,null],[null,null,null,null,null,null],[null,null,null,null,null,null],[null,null,null,null,null,null],[null,null,null,null,null,null],[null,null,null,null,null,null],[null,null,null,null,null,null],[null,null,null,null,null,null],[null,null,null,null,null,null],[null,null,null,null,null,null],[null,null,null,null,null,null],[null,null,null,null,null,null],[null,null,null,null,null,null],[null,null,null,null,null,null],[null,null,null,null,null,null],[null,null,null,null,null,null],[null,null,null,null,null,null],["Uncertain - too short","Uncertain - too short","Uncertain - too short","Uncertain - too short","Uncertain - too short","Uncertain - too short"],[null,null,null,null,null,null],[null,null,null,null,null,null],[null,null,null,null,null,null],[null,null,null,null,null,null],["k_Bacteria;n_Terrabacteria group;p_Firmicutes;c_Erysipelotrichia;o_Erysipelotrichales;f_Erysipelotrichaceae;g_Massilimicrobiota","k_Bacteria;n_Terrabacteria group;p_Firmicutes;c_Clostridia;o_Eubacteriales;n_unclassified Eubacteriales;s_Clostridiales bacterium","k_Bacteria;n_Terrabacteria group;p_Firmicutes;c_Clostridia;o_Eubacteriales;n_unclassified Eubacteriales;s_Clostridiales bacterium","k_Bacteria;n_Terrabacteria group;p_Firmicutes;c_Clostridia;o_Eubacteriales","k_Bacteria;n_Terrabacteria group;p_Firmicutes;c_Clostridia;o_Eubacteriales","k_Bacteria;n_Terrabacteria group;p_Firmicutes;c_Bacilli;o_Lactobacillales;f_Enterococcaceae;g_Enterococcus;s_Enterococcus faecium"],[0,0,0,0,0,0],[34.5,39.67,39.67,51.13,51.13,38.89],[313,305,305,221,221,324],[1,2,2,2,2,2]],"container":"<table class=\"display\">\n  <thead>\n    <tr>\n      <th> <\/th>\n      <th>ORF_ID<\/th>\n      <th>Contig_ID<\/th>\n      <th>SQM_Length.NT<\/th>\n      <th>SQM_Length.AA<\/th>\n      <th>SQM_GC.perc<\/th>\n      <th>SQM_Gene.name<\/th>\n      <th>SQM_Tax_orf<\/th>\n      <th>SQM_KEGG.ID<\/th>\n      <th>SQM_KEGGFUN<\/th>\n      <th>SQM_KEGGPATH<\/th>\n      <th>SQM_COG.ID<\/th>\n      <th>SQM_COGFUN<\/th>\n      <th>SQM_COGPATH<\/th>\n      <th>SQM_PFAM<\/th>\n      <th>rgi_CARD_Cut_Off<\/th>\n      <th>rgi_CARD_Pass_Bitscore<\/th>\n      <th>rgi_CARD_Best_Hit_Bitscore<\/th>\n      <th>rgi_CARD_Best_Hit_ARO<\/th>\n      <th>rgi_CARD_Best_Identities<\/th>\n      <th>rgi_CARD_ARO<\/th>\n      <th>rgi_CARD_Model_type<\/th>\n      <th>rgi_CARD_SNPs_in_Best_Hit_ARO<\/th>\n      <th>rgi_CARD_Other_SNPs<\/th>\n      <th>rgi_CARD_Drug Class<\/th>\n      <th>rgi_CARD_Resistance Mechanism<\/th>\n      <th>rgi_CARD_AMR Gene Family<\/th>\n      <th>rgi_CARD_Percentage Length of Reference Sequence<\/th>\n      <th>rgi_CARD_ID<\/th>\n      <th>rgi_CARD_Model_ID<\/th>\n      <th>rgi_CARD_Nudged<\/th>\n      <th>rgi_CARD_Note<\/th>\n      <th>PathoFact_AMR_ARG<\/th>\n      <th>PathoFact_AMR_ARG_SNPs<\/th>\n      <th>PathoFact_AMR_AMR_category<\/th>\n      <th>PathoFact_AMR_AMR_sub_class<\/th>\n      <th>PathoFact_AMR_Resistance_mechanism<\/th>\n      <th>PathoFact_AMR_Database<\/th>\n      <th>PathoFact_AMR_MGE_prediction<\/th>\n      <th>PathoFact_Tox_Toxin_classification<\/th>\n      <th>viralverify_Prediction<\/th>\n      <th>isescan_family<\/th>\n      <th>isescan_cluster<\/th>\n      <th>isescan_type<\/th>\n      <th>isescan_count_IS_contig<\/th>\n      <th>SQM_Tax_contig<\/th>\n      <th>SQM_Disparity_contig<\/th>\n      <th>SQM_GC_perc_contig<\/th>\n      <th>SQM_Length_contig<\/th>\n      <th>SQM_Num_genes_contig<\/th>\n    <\/tr>\n  <\/thead>\n<\/table>","options":{"columnDefs":[{"className":"dt-right","targets":[3,4,5,16,17,19,20,27,29,44,46,47,48,49]},{"orderable":false,"targets":0}],"order":[],"autoWidth":false,"orderClasses":false}},"evals":[],"jsHooks":[]}</script>
```
## Test anvio no bin:


```r
anvio_nobin_summary %>% 
  here::here() %>%  
  read_tsv() -> anvi_no_bin_annot
```

```
## Rows: 352828 Columns: 22
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
<div id="htmlwidget-503c59acba863ff13ca1" style="width:100%;height:auto;" class="datatables html-widget"></div>
<script type="application/json" data-for="htmlwidget-503c59acba863ff13ca1">{"x":{"filter":"none","vertical":false,"data":[["1","2","3","4","5","6"],["megahit_94475","megahit_135602","megahit_101245","megahit_222020","megahit_229957","megahit_94514"],["maxbin_154_fasta_contigs","maxbin_154_fasta_contigs","maxbin_154_fasta_contigs","maxbin_154_fasta_contigs","maxbin_154_fasta_contigs","maxbin_154_fasta_contigs"]],"container":"<table class=\"display\">\n  <thead>\n    <tr>\n      <th> <\/th>\n      <th>Contig_ID<\/th>\n      <th>DAS_bin<\/th>\n    <\/tr>\n  <\/thead>\n<\/table>","options":{"columnDefs":[{"orderable":false,"targets":0}],"order":[],"autoWidth":false,"orderClasses":false}},"evals":[],"jsHooks":[]}</script>
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
<div id="htmlwidget-a1a0181f0aa91748e6c8" style="width:100%;height:auto;" class="datatables html-widget"></div>
<script type="application/json" data-for="htmlwidget-a1a0181f0aa91748e6c8">{"x":{"filter":"none","vertical":false,"data":[["1","2","3","4","5","6"],["megahit_69635","megahit_218792","megahit_96692","megahit_139453","megahit_225023","megahit_105955"],["Bin_8_1","Bin_8_1","Bin_8_1","Bin_8_1","Bin_8_1","Bin_8_1"]],"container":"<table class=\"display\">\n  <thead>\n    <tr>\n      <th> <\/th>\n      <th>Contig_ID<\/th>\n      <th>CONCOCT_bin<\/th>\n    <\/tr>\n  <\/thead>\n<\/table>","options":{"columnDefs":[{"orderable":false,"targets":0}],"order":[],"autoWidth":false,"orderClasses":false}},"evals":[],"jsHooks":[]}</script>
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
<div id="htmlwidget-cba7bd33bde8bb3f9afe" style="width:100%;height:auto;" class="datatables html-widget"></div>
<script type="application/json" data-for="htmlwidget-cba7bd33bde8bb3f9afe">{"x":{"filter":"none","vertical":false,"data":[["1","2","3","4","5","6"],["maxbin_002_fasta_contigs","maxbin_004_fasta_contigs","maxbin_005_fasta_contigs","maxbin_006_fasta_contigs","maxbin_008_fasta_contigs","maxbin_009_fasta_contigs"],["d__Bacteria;p__Bacteroidota;c__Bacteroidia;o__Bacteroidales;f__Bacteroidaceae;g__Phocaeicola;s__Phocaeicola plebeius_A","d__Bacteria;p__Firmicutes_A;c__Clostridia;o__Oscillospirales;f__Acutalibacteraceae;g__Ruminococcus_E;s__Ruminococcus_E bromii_B","d__Bacteria;p__Firmicutes_A;c__Clostridia;o__Oscillospirales;f__Ruminococcaceae;g__Faecalibacterium;s__Faecalibacterium prausnitzii_C","d__Bacteria;p__Firmicutes_A;c__Clostridia;o__Lachnospirales;f__Lachnospiraceae;g__Agathobacter;s__Agathobacter rectalis","d__Bacteria;p__Firmicutes_A;c__Clostridia;o__Lachnospirales;f__Lachnospiraceae;g__RUG115;s__RUG115 sp900066395","d__Bacteria;p__Firmicutes_A;c__Clostridia;o__Lachnospirales;f__Lachnospiraceae;g__Anaerostipes;s__Anaerostipes hadrus"]],"container":"<table class=\"display\">\n  <thead>\n    <tr>\n      <th> <\/th>\n      <th>user_genome<\/th>\n      <th>GTDB_tax<\/th>\n    <\/tr>\n  <\/thead>\n<\/table>","options":{"columnDefs":[{"orderable":false,"targets":0}],"order":[],"autoWidth":false,"orderClasses":false}},"evals":[],"jsHooks":[]}</script>
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
<div id="htmlwidget-760f3be48f491c07c2b6" style="width:100%;height:auto;" class="datatables html-widget"></div>
<script type="application/json" data-for="htmlwidget-760f3be48f491c07c2b6">{"x":{"filter":"none","vertical":false,"data":[["1","2","3","4","5","6"],["maxbin.163.fasta.contigs","maxbin.043.fasta.contigs","maxbin.075.fasta.contigs","maxbin.153.fasta.contigs","maxbin.094.fasta.contigs","maxbin.171.fasta.contigs"],["k_Bacteria;n_Terrabacteria group;p_Firmicutes;c_Clostridia;o_Eubacteriales;f_Peptostreptococcaceae;g_Clostridioides;s_Clostridioides difficile","k_Bacteria;n_Terrabacteria group;p_Firmicutes;c_Bacilli;o_Lactobacillales;f_Enterococcaceae;g_Enterococcus","k_Bacteria;n_Terrabacteria group;p_Firmicutes;c_Erysipelotrichia;o_Erysipelotrichales","k_Bacteria;n_Terrabacteria group;p_Firmicutes;c_Bacilli;o_Lactobacillales;f_Enterococcaceae;g_Enterococcus","k_Bacteria;n_Terrabacteria group;p_Firmicutes;c_Clostridia;o_Eubacteriales","k_Bacteria;n_Terrabacteria group;p_Actinobacteria;c_Coriobacteriia;o_Coriobacteriales"],[null,null,null,"superkingdom:Bacteria;clade:Terrabacteria group;phylum:Firmicutes;class:Bacilli;order:Lactobacillales;family:Enterococcaceae;genus:Enterococcus",null,"superkingdom:Bacteria;clade:Terrabacteria group;phylum:Actinobacteria;class:Coriobacteriia;order:Coriobacteriales;family:Coriobacteriaceae"],[4025395,3586160,3293471,2805040,3182586,1869581],[28.35,42.33,31.31,37.88,55.16,49.32],[47,60,51,119,162,38],[0.094,0,0.08,0,0.013,0.066],[100,99.61,99.37,99.19,99.09,98.71],[0.33,0.62,0,0.78,12.62,0],[0,0,0,81.25,7.69,0],["HQ","HQ","HQ","HQ",null,"HQ"]],"container":"<table class=\"display\">\n  <thead>\n    <tr>\n      <th> <\/th>\n      <th>SQM_DAS_Bin ID<\/th>\n      <th>SQM_DAS_Tax<\/th>\n      <th>SQM_DAS_Tax 16S<\/th>\n      <th>SQM_DAS_Length<\/th>\n      <th>SQM_DAS_GC perc<\/th>\n      <th>SQM_DAS_Num contigs<\/th>\n      <th>SQM_DAS_Disparity<\/th>\n      <th>SQM_DAS_Completeness<\/th>\n      <th>SQM_DAS_Contamination<\/th>\n      <th>SQM_DAS_Strain het<\/th>\n      <th>SQM_DAS_HQ<\/th>\n    <\/tr>\n  <\/thead>\n<\/table>","options":{"columnDefs":[{"className":"dt-right","targets":[4,5,6,7,8,9,10]},{"orderable":false,"targets":0}],"order":[],"autoWidth":false,"orderClasses":false}},"evals":[],"jsHooks":[]}</script>
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
<div id="htmlwidget-c86ababd5bdcc7ba92b0" style="width:100%;height:auto;" class="datatables html-widget"></div>
<script type="application/json" data-for="htmlwidget-c86ababd5bdcc7ba92b0">{"x":{"filter":"none","vertical":false,"data":[["1","2","3","4","5","6"],["maxbin_002_fasta_contigs","maxbin_004_fasta_contigs","maxbin_005_fasta_contigs","maxbin_006_fasta_contigs","maxbin_008_fasta_contigs","maxbin_009_fasta_contigs"],[3373122,2078616,3887172,2170300,2596485,2676589],[99,36,465,288,54,325],[65526,155262,22282,13284,117427,14386],[44.100777643127,40.7792212729088,57.8414339713452,41.1785262144328,40.7337175045291,37.0186970597009],[88.7323943661972,91.5492957746479,80.2816901408451,39.4366197183099,64.7887323943662,95.7746478873239],[2.8169014084507,4.22535211267606,22.5352112676056,0,11.2676056338028,8.45070422535211],["HQ","HQ",null,null,null,"HQ"]],"container":"<table class=\"display\">\n  <thead>\n    <tr>\n      <th> <\/th>\n      <th>ANVIO_SUMMARY_DAS_bins<\/th>\n      <th>ANVIO_SUMMARY_DAS_total_length<\/th>\n      <th>ANVIO_SUMMARY_DAS_num_contigs<\/th>\n      <th>ANVIO_SUMMARY_DAS_N50<\/th>\n      <th>ANVIO_SUMMARY_DAS_GC_content<\/th>\n      <th>ANVIO_SUMMARY_DAS_percent_completion<\/th>\n      <th>ANVIO_SUMMARY_DAS_percent_redundancy<\/th>\n      <th>ANVIO_DAS_HQ<\/th>\n    <\/tr>\n  <\/thead>\n<\/table>","options":{"columnDefs":[{"className":"dt-right","targets":[2,3,4,5,6,7]},{"orderable":false,"targets":0}],"order":[],"autoWidth":false,"orderClasses":false}},"evals":[],"jsHooks":[]}</script>
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
<div id="htmlwidget-79ab01b3c558bda84f48" style="width:100%;height:auto;" class="datatables html-widget"></div>
<script type="application/json" data-for="htmlwidget-79ab01b3c558bda84f48">{"x":{"filter":"none","vertical":false,"data":[["1","2","3","4","5","6"],["Bin_0_1_1","Bin_0_2_1","Bin_0_3_1","Bin_0_4_1","Bin_0_5_1","Bin_0_6_1"],[4997430,2611266,1816009,3210323,3351680,2332080],[93,118,22,74,21,31],[95529,48119,147912,71289,221494,117575],[44.8850463626556,45.7147524370042,44.5251091026598,44.3226862025176,59.2794959337204,43.7523527454812],[98.5915492957746,98.5915492957746,88.7323943661972,85.9154929577465,80.2816901408451,100],[1.40845070422535,0,0,2.8169014084507,1.40845070422535,1.40845070422535],["HQ","HQ","HQ","HQ","HQ","HQ"]],"container":"<table class=\"display\">\n  <thead>\n    <tr>\n      <th> <\/th>\n      <th>ANVIO_SUMMARY_CONCOCT_bins<\/th>\n      <th>ANVIO_SUMMARY_CONCOCT_total_length<\/th>\n      <th>ANVIO_SUMMARY_CONCOCT_num_contigs<\/th>\n      <th>ANVIO_SUMMARY_CONCOCT_N50<\/th>\n      <th>ANVIO_SUMMARY_CONCOCT_GC_content<\/th>\n      <th>ANVIO_SUMMARY_CONCOCT_percent_completion<\/th>\n      <th>ANVIO_SUMMARY_CONCOCT_percent_redundancy<\/th>\n      <th>ANVIO_CONCOCT_HQ<\/th>\n    <\/tr>\n  <\/thead>\n<\/table>","options":{"columnDefs":[{"className":"dt-right","targets":[2,3,4,5,6,7]},{"orderable":false,"targets":0}],"order":[],"autoWidth":false,"orderClasses":false}},"evals":[],"jsHooks":[]}</script>
```


```r
orf_contig_3 %>% 
  left_join(ANVIO_SUMMARY_CONCOCT,
            by = c("CONCOCT_bin" = "ANVIO_SUMMARY_CONCOCT_bins")) -> orf_contig_3
```


Rename ORF and contig ID:

```r
orf_contig_3 %>% 
  mutate(ORF_ID = paste0("h_", ORF_ID),
         Contig_ID = paste0("h_", Contig_ID)) -> full_gene_catalogue
```

##  Extract quantitative information from SQM:

### ORFs:



```r
human$orfs$table %>%
  rownames_to_column("ORF_ID") %>% 
  data.frame() %>% 
  select(ORF_ID, TPM.D.1:Raw.base.count.HV.24) %>% 
  # select(starts_with("TPM"), "ORF_ID") %>% 
  column_to_rownames("ORF_ID") -> quant_orfs

# colnames(quant) <- colnames(quant) %>% 
#   str_replace_all("//.", "//_")

rownames(quant_orfs) <- paste0("h_",
                               rownames(quant_orfs))

quant_orfs %>% 
  rownames_to_column("ORF_ID") -> quant_orfs

quant_orfs %>% 
  head() %>% 
  DT::datatable()
```

```{=html}
<div id="htmlwidget-013cf31905d48effde03" style="width:100%;height:auto;" class="datatables html-widget"></div>
<script type="application/json" data-for="htmlwidget-013cf31905d48effde03">{"x":{"filter":"none","vertical":false,"data":[["1","2","3","4","5","6"],["h_megahit_1_1-312","h_megahit_2_3-176","h_megahit_2_173-304","h_megahit_3_3-116","h_megahit_3_82-219","h_megahit_4_1-72"],[0.165,0,0,0,0,0],[0,0.229,0.151,0,0,0],[0,0,0,0.298,0.246,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0.23,0.304,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0.405],[0,0,0,0.156,0.129,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0.178,0.147,0],[0,0,0,0.141,0.116,0],[0,0,0,0,0,0],[0,0,0,0.154,0.127,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0.553],[0,0,0,0,0,0],[0,0,0,0.165,0.137,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0.434,0.358,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[3.349,0,0,0,0,0],[0,1.149,0.78,0,0,0],[0,0,0,1.684,1.261,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,1.172,0.773,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,1.972],[0,0,0,0.404,0.993,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0.86,0.623,0],[0,0,0,0.746,0.717,0],[0,0,0,0,0,0],[0,0,0,0.561,0.87,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0.056],[0,0,0,0,0,0],[0,0,0,0.096,0.768,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,1.974,2.37,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[7,0,0,0,0,0],[0,2,1,0,0,0],[0,0,0,2,2,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,2,2,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,2],[0,0,0,1,1,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,1,1,0],[0,0,0,1,1,0],[0,0,0,0,0,0],[0,0,0,1,1,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,2],[0,0,0,0,0,0],[0,0,0,1,1,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,3,3,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[1045,0,0,0,0,0],[0,200,103,0,0,0],[0,0,0,192,174,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,204,102,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,142],[0,0,0,46,137,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,98,86,0],[0,0,0,85,99,0],[0,0,0,0,0,0],[0,0,0,64,120,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,4],[0,0,0,0,0,0],[0,0,0,11,106,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,225,327,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0]],"container":"<table class=\"display\">\n  <thead>\n    <tr>\n      <th> <\/th>\n      <th>ORF_ID<\/th>\n      <th>TPM.D.1<\/th>\n      <th>TPM.HC.1<\/th>\n      <th>TPM.HC.2<\/th>\n      <th>TPM.HC.3<\/th>\n      <th>TPM.HC.4<\/th>\n      <th>TPM.HC.5<\/th>\n      <th>TPM.HC.6<\/th>\n      <th>TPM.HC.7<\/th>\n      <th>TPM.HC.8<\/th>\n      <th>TPM.HC.9<\/th>\n      <th>TPM.HC.10<\/th>\n      <th>TPM.HC.11<\/th>\n      <th>TPM.HC.12<\/th>\n      <th>TPM.HC.13<\/th>\n      <th>TPM.HC.14<\/th>\n      <th>TPM.HC.15<\/th>\n      <th>TPM.HC.16<\/th>\n      <th>TPM.HC.17<\/th>\n      <th>TPM.HC.18<\/th>\n      <th>TPM.HC.19<\/th>\n      <th>TPM.HC.20<\/th>\n      <th>TPM.HC.21<\/th>\n      <th>TPM.HC.22<\/th>\n      <th>TPM.HC.23<\/th>\n      <th>TPM.HC.24<\/th>\n      <th>TPM.HV.1<\/th>\n      <th>TPM.HV.2<\/th>\n      <th>TPM.HV.3<\/th>\n      <th>TPM.HV.4<\/th>\n      <th>TPM.HV.5<\/th>\n      <th>TPM.HV.6<\/th>\n      <th>TPM.HV.7<\/th>\n      <th>TPM.HV.8<\/th>\n      <th>TPM.HV.9<\/th>\n      <th>TPM.HV.10<\/th>\n      <th>TPM.HV.11<\/th>\n      <th>TPM.HV.12<\/th>\n      <th>TPM.HV.13<\/th>\n      <th>TPM.HV.14<\/th>\n      <th>TPM.HV.15<\/th>\n      <th>TPM.HV.16<\/th>\n      <th>TPM.HV.17<\/th>\n      <th>TPM.HV.18<\/th>\n      <th>TPM.HV.19<\/th>\n      <th>TPM.HV.20<\/th>\n      <th>TPM.HV.21<\/th>\n      <th>TPM.HV.22<\/th>\n      <th>TPM.HV.23<\/th>\n      <th>TPM.HV.24<\/th>\n      <th>Coverage.D.1<\/th>\n      <th>Coverage.HC.1<\/th>\n      <th>Coverage.HC.2<\/th>\n      <th>Coverage.HC.3<\/th>\n      <th>Coverage.HC.4<\/th>\n      <th>Coverage.HC.5<\/th>\n      <th>Coverage.HC.6<\/th>\n      <th>Coverage.HC.7<\/th>\n      <th>Coverage.HC.8<\/th>\n      <th>Coverage.HC.9<\/th>\n      <th>Coverage.HC.10<\/th>\n      <th>Coverage.HC.11<\/th>\n      <th>Coverage.HC.12<\/th>\n      <th>Coverage.HC.13<\/th>\n      <th>Coverage.HC.14<\/th>\n      <th>Coverage.HC.15<\/th>\n      <th>Coverage.HC.16<\/th>\n      <th>Coverage.HC.17<\/th>\n      <th>Coverage.HC.18<\/th>\n      <th>Coverage.HC.19<\/th>\n      <th>Coverage.HC.20<\/th>\n      <th>Coverage.HC.21<\/th>\n      <th>Coverage.HC.22<\/th>\n      <th>Coverage.HC.23<\/th>\n      <th>Coverage.HC.24<\/th>\n      <th>Coverage.HV.1<\/th>\n      <th>Coverage.HV.2<\/th>\n      <th>Coverage.HV.3<\/th>\n      <th>Coverage.HV.4<\/th>\n      <th>Coverage.HV.5<\/th>\n      <th>Coverage.HV.6<\/th>\n      <th>Coverage.HV.7<\/th>\n      <th>Coverage.HV.8<\/th>\n      <th>Coverage.HV.9<\/th>\n      <th>Coverage.HV.10<\/th>\n      <th>Coverage.HV.11<\/th>\n      <th>Coverage.HV.12<\/th>\n      <th>Coverage.HV.13<\/th>\n      <th>Coverage.HV.14<\/th>\n      <th>Coverage.HV.15<\/th>\n      <th>Coverage.HV.16<\/th>\n      <th>Coverage.HV.17<\/th>\n      <th>Coverage.HV.18<\/th>\n      <th>Coverage.HV.19<\/th>\n      <th>Coverage.HV.20<\/th>\n      <th>Coverage.HV.21<\/th>\n      <th>Coverage.HV.22<\/th>\n      <th>Coverage.HV.23<\/th>\n      <th>Coverage.HV.24<\/th>\n      <th>Raw.read.count.D.1<\/th>\n      <th>Raw.read.count.HC.1<\/th>\n      <th>Raw.read.count.HC.2<\/th>\n      <th>Raw.read.count.HC.3<\/th>\n      <th>Raw.read.count.HC.4<\/th>\n      <th>Raw.read.count.HC.5<\/th>\n      <th>Raw.read.count.HC.6<\/th>\n      <th>Raw.read.count.HC.7<\/th>\n      <th>Raw.read.count.HC.8<\/th>\n      <th>Raw.read.count.HC.9<\/th>\n      <th>Raw.read.count.HC.10<\/th>\n      <th>Raw.read.count.HC.11<\/th>\n      <th>Raw.read.count.HC.12<\/th>\n      <th>Raw.read.count.HC.13<\/th>\n      <th>Raw.read.count.HC.14<\/th>\n      <th>Raw.read.count.HC.15<\/th>\n      <th>Raw.read.count.HC.16<\/th>\n      <th>Raw.read.count.HC.17<\/th>\n      <th>Raw.read.count.HC.18<\/th>\n      <th>Raw.read.count.HC.19<\/th>\n      <th>Raw.read.count.HC.20<\/th>\n      <th>Raw.read.count.HC.21<\/th>\n      <th>Raw.read.count.HC.22<\/th>\n      <th>Raw.read.count.HC.23<\/th>\n      <th>Raw.read.count.HC.24<\/th>\n      <th>Raw.read.count.HV.1<\/th>\n      <th>Raw.read.count.HV.2<\/th>\n      <th>Raw.read.count.HV.3<\/th>\n      <th>Raw.read.count.HV.4<\/th>\n      <th>Raw.read.count.HV.5<\/th>\n      <th>Raw.read.count.HV.6<\/th>\n      <th>Raw.read.count.HV.7<\/th>\n      <th>Raw.read.count.HV.8<\/th>\n      <th>Raw.read.count.HV.9<\/th>\n      <th>Raw.read.count.HV.10<\/th>\n      <th>Raw.read.count.HV.11<\/th>\n      <th>Raw.read.count.HV.12<\/th>\n      <th>Raw.read.count.HV.13<\/th>\n      <th>Raw.read.count.HV.14<\/th>\n      <th>Raw.read.count.HV.15<\/th>\n      <th>Raw.read.count.HV.16<\/th>\n      <th>Raw.read.count.HV.17<\/th>\n      <th>Raw.read.count.HV.18<\/th>\n      <th>Raw.read.count.HV.19<\/th>\n      <th>Raw.read.count.HV.20<\/th>\n      <th>Raw.read.count.HV.21<\/th>\n      <th>Raw.read.count.HV.22<\/th>\n      <th>Raw.read.count.HV.23<\/th>\n      <th>Raw.read.count.HV.24<\/th>\n      <th>Raw.base.count.D.1<\/th>\n      <th>Raw.base.count.HC.1<\/th>\n      <th>Raw.base.count.HC.2<\/th>\n      <th>Raw.base.count.HC.3<\/th>\n      <th>Raw.base.count.HC.4<\/th>\n      <th>Raw.base.count.HC.5<\/th>\n      <th>Raw.base.count.HC.6<\/th>\n      <th>Raw.base.count.HC.7<\/th>\n      <th>Raw.base.count.HC.8<\/th>\n      <th>Raw.base.count.HC.9<\/th>\n      <th>Raw.base.count.HC.10<\/th>\n      <th>Raw.base.count.HC.11<\/th>\n      <th>Raw.base.count.HC.12<\/th>\n      <th>Raw.base.count.HC.13<\/th>\n      <th>Raw.base.count.HC.14<\/th>\n      <th>Raw.base.count.HC.15<\/th>\n      <th>Raw.base.count.HC.16<\/th>\n      <th>Raw.base.count.HC.17<\/th>\n      <th>Raw.base.count.HC.18<\/th>\n      <th>Raw.base.count.HC.19<\/th>\n      <th>Raw.base.count.HC.20<\/th>\n      <th>Raw.base.count.HC.21<\/th>\n      <th>Raw.base.count.HC.22<\/th>\n      <th>Raw.base.count.HC.23<\/th>\n      <th>Raw.base.count.HC.24<\/th>\n      <th>Raw.base.count.HV.1<\/th>\n      <th>Raw.base.count.HV.2<\/th>\n      <th>Raw.base.count.HV.3<\/th>\n      <th>Raw.base.count.HV.4<\/th>\n      <th>Raw.base.count.HV.5<\/th>\n      <th>Raw.base.count.HV.6<\/th>\n      <th>Raw.base.count.HV.7<\/th>\n      <th>Raw.base.count.HV.8<\/th>\n      <th>Raw.base.count.HV.9<\/th>\n      <th>Raw.base.count.HV.10<\/th>\n      <th>Raw.base.count.HV.11<\/th>\n      <th>Raw.base.count.HV.12<\/th>\n      <th>Raw.base.count.HV.13<\/th>\n      <th>Raw.base.count.HV.14<\/th>\n      <th>Raw.base.count.HV.15<\/th>\n      <th>Raw.base.count.HV.16<\/th>\n      <th>Raw.base.count.HV.17<\/th>\n      <th>Raw.base.count.HV.18<\/th>\n      <th>Raw.base.count.HV.19<\/th>\n      <th>Raw.base.count.HV.20<\/th>\n      <th>Raw.base.count.HV.21<\/th>\n      <th>Raw.base.count.HV.22<\/th>\n      <th>Raw.base.count.HV.23<\/th>\n      <th>Raw.base.count.HV.24<\/th>\n    <\/tr>\n  <\/thead>\n<\/table>","options":{"columnDefs":[{"className":"dt-right","targets":[2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58,59,60,61,62,63,64,65,66,67,68,69,70,71,72,73,74,75,76,77,78,79,80,81,82,83,84,85,86,87,88,89,90,91,92,93,94,95,96,97,98,99,100,101,102,103,104,105,106,107,108,109,110,111,112,113,114,115,116,117,118,119,120,121,122,123,124,125,126,127,128,129,130,131,132,133,134,135,136,137,138,139,140,141,142,143,144,145,146,147,148,149,150,151,152,153,154,155,156,157,158,159,160,161,162,163,164,165,166,167,168,169,170,171,172,173,174,175,176,177,178,179,180,181,182,183,184,185,186,187,188,189,190,191,192,193,194,195,196,197]},{"orderable":false,"targets":0}],"order":[],"autoWidth":false,"orderClasses":false}},"evals":[],"jsHooks":[]}</script>
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
  select(Contig_ID, Coverage.D.1:Raw.read.count.HV.24) %>% 
  # select(starts_with("TPM"), "ORF_ID") %>% 
  column_to_rownames("Contig_ID") -> quant_contigs

# colnames(quant_contigs) <- colnames(quant_contigs) %>% 
#   str_replace_all("//.", "//_")

rownames(quant_contigs) <- paste0("h_",
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
<div id="htmlwidget-340e456925d71a8e31a1" style="width:100%;height:auto;" class="datatables html-widget"></div>
<script type="application/json" data-for="htmlwidget-340e456925d71a8e31a1">{"x":{"filter":"none","vertical":false,"data":[["1","2","3","4","5","6"],["h_megahit_1","h_megahit_2","h_megahit_3","h_megahit_4","h_megahit_5","h_megahit_6"],[3.83,0,0,0,0,0],[1.5,0,0,0,0,0],[8,0,0,0,0,0],[0,0.98,0,0,0,0],[0,1.5,0,0,0,0],[0,2,0,0,0,0],[0,0,1.36,0,0.88,0],[0,0,1.7,0,1.7,0],[0,0,2,0,3,0],[0,0,0,0,0.88,0],[0,0,0,0,1.4,0],[0,0,0,0,2,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0.13,0],[0,0,0,0,0.7,0],[0,0,0,0,1,0],[0,0.98,0,0,0,0],[0,1.6,0,0,0,0],[0,2,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,1.85,0,0],[0,0,0,2.8,0,0],[0,0,0,4,0,0],[0,0,0.68,0,0,0],[0,0,1.1,0,0,0],[0,0,1,0,0,0],[0,0,0,0,0.21,0],[0,0,0,0,0.6,0],[0,0,0,0,1,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0.68,0,0,0],[0,0,1.4,0,0,0],[0,0,1,0,0,0],[0,0,0.68,0,0,0],[0,0,0.9,0,0,0],[0,0,1,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0.68,0,0,0],[0,0,1,0,0,0],[0,0,1,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0.93,0,0],[0,0,0,0.8,0,0],[0,0,0,2,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0.92,0,0],[0,0,0,2.1,0,0],[0,0,0,2,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0.92,0,0,0],[0,0,1.6,0,0,0],[0,0,2,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,2.04,0,0,0],[0,0,2.2,0,0,0],[0,0,3,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0]],"container":"<table class=\"display\">\n  <thead>\n    <tr>\n      <th> <\/th>\n      <th>Contig_ID<\/th>\n      <th>Coverage.D.1<\/th>\n      <th>TPM.D.1<\/th>\n      <th>Raw.read.count.D.1<\/th>\n      <th>Coverage.HC.1<\/th>\n      <th>TPM.HC.1<\/th>\n      <th>Raw.read.count.HC.1<\/th>\n      <th>Coverage.HC.2<\/th>\n      <th>TPM.HC.2<\/th>\n      <th>Raw.read.count.HC.2<\/th>\n      <th>Coverage.HC.3<\/th>\n      <th>TPM.HC.3<\/th>\n      <th>Raw.read.count.HC.3<\/th>\n      <th>Coverage.HC.4<\/th>\n      <th>TPM.HC.4<\/th>\n      <th>Raw.read.count.HC.4<\/th>\n      <th>Coverage.HC.5<\/th>\n      <th>TPM.HC.5<\/th>\n      <th>Raw.read.count.HC.5<\/th>\n      <th>Coverage.HC.6<\/th>\n      <th>TPM.HC.6<\/th>\n      <th>Raw.read.count.HC.6<\/th>\n      <th>Coverage.HC.7<\/th>\n      <th>TPM.HC.7<\/th>\n      <th>Raw.read.count.HC.7<\/th>\n      <th>Coverage.HC.8<\/th>\n      <th>TPM.HC.8<\/th>\n      <th>Raw.read.count.HC.8<\/th>\n      <th>Coverage.HC.9<\/th>\n      <th>TPM.HC.9<\/th>\n      <th>Raw.read.count.HC.9<\/th>\n      <th>Coverage.HC.10<\/th>\n      <th>TPM.HC.10<\/th>\n      <th>Raw.read.count.HC.10<\/th>\n      <th>Coverage.HC.11<\/th>\n      <th>TPM.HC.11<\/th>\n      <th>Raw.read.count.HC.11<\/th>\n      <th>Coverage.HC.12<\/th>\n      <th>TPM.HC.12<\/th>\n      <th>Raw.read.count.HC.12<\/th>\n      <th>Coverage.HC.13<\/th>\n      <th>TPM.HC.13<\/th>\n      <th>Raw.read.count.HC.13<\/th>\n      <th>Coverage.HC.14<\/th>\n      <th>TPM.HC.14<\/th>\n      <th>Raw.read.count.HC.14<\/th>\n      <th>Coverage.HC.15<\/th>\n      <th>TPM.HC.15<\/th>\n      <th>Raw.read.count.HC.15<\/th>\n      <th>Coverage.HC.16<\/th>\n      <th>TPM.HC.16<\/th>\n      <th>Raw.read.count.HC.16<\/th>\n      <th>Coverage.HC.17<\/th>\n      <th>TPM.HC.17<\/th>\n      <th>Raw.read.count.HC.17<\/th>\n      <th>Coverage.HC.18<\/th>\n      <th>TPM.HC.18<\/th>\n      <th>Raw.read.count.HC.18<\/th>\n      <th>Coverage.HC.19<\/th>\n      <th>TPM.HC.19<\/th>\n      <th>Raw.read.count.HC.19<\/th>\n      <th>Coverage.HC.20<\/th>\n      <th>TPM.HC.20<\/th>\n      <th>Raw.read.count.HC.20<\/th>\n      <th>Coverage.HC.21<\/th>\n      <th>TPM.HC.21<\/th>\n      <th>Raw.read.count.HC.21<\/th>\n      <th>Coverage.HC.22<\/th>\n      <th>TPM.HC.22<\/th>\n      <th>Raw.read.count.HC.22<\/th>\n      <th>Coverage.HC.23<\/th>\n      <th>TPM.HC.23<\/th>\n      <th>Raw.read.count.HC.23<\/th>\n      <th>Coverage.HC.24<\/th>\n      <th>TPM.HC.24<\/th>\n      <th>Raw.read.count.HC.24<\/th>\n      <th>Coverage.HV.1<\/th>\n      <th>TPM.HV.1<\/th>\n      <th>Raw.read.count.HV.1<\/th>\n      <th>Coverage.HV.2<\/th>\n      <th>TPM.HV.2<\/th>\n      <th>Raw.read.count.HV.2<\/th>\n      <th>Coverage.HV.3<\/th>\n      <th>TPM.HV.3<\/th>\n      <th>Raw.read.count.HV.3<\/th>\n      <th>Coverage.HV.4<\/th>\n      <th>TPM.HV.4<\/th>\n      <th>Raw.read.count.HV.4<\/th>\n      <th>Coverage.HV.5<\/th>\n      <th>TPM.HV.5<\/th>\n      <th>Raw.read.count.HV.5<\/th>\n      <th>Coverage.HV.6<\/th>\n      <th>TPM.HV.6<\/th>\n      <th>Raw.read.count.HV.6<\/th>\n      <th>Coverage.HV.7<\/th>\n      <th>TPM.HV.7<\/th>\n      <th>Raw.read.count.HV.7<\/th>\n      <th>Coverage.HV.8<\/th>\n      <th>TPM.HV.8<\/th>\n      <th>Raw.read.count.HV.8<\/th>\n      <th>Coverage.HV.9<\/th>\n      <th>TPM.HV.9<\/th>\n      <th>Raw.read.count.HV.9<\/th>\n      <th>Coverage.HV.10<\/th>\n      <th>TPM.HV.10<\/th>\n      <th>Raw.read.count.HV.10<\/th>\n      <th>Coverage.HV.11<\/th>\n      <th>TPM.HV.11<\/th>\n      <th>Raw.read.count.HV.11<\/th>\n      <th>Coverage.HV.12<\/th>\n      <th>TPM.HV.12<\/th>\n      <th>Raw.read.count.HV.12<\/th>\n      <th>Coverage.HV.13<\/th>\n      <th>TPM.HV.13<\/th>\n      <th>Raw.read.count.HV.13<\/th>\n      <th>Coverage.HV.14<\/th>\n      <th>TPM.HV.14<\/th>\n      <th>Raw.read.count.HV.14<\/th>\n      <th>Coverage.HV.15<\/th>\n      <th>TPM.HV.15<\/th>\n      <th>Raw.read.count.HV.15<\/th>\n      <th>Coverage.HV.16<\/th>\n      <th>TPM.HV.16<\/th>\n      <th>Raw.read.count.HV.16<\/th>\n      <th>Coverage.HV.17<\/th>\n      <th>TPM.HV.17<\/th>\n      <th>Raw.read.count.HV.17<\/th>\n      <th>Coverage.HV.18<\/th>\n      <th>TPM.HV.18<\/th>\n      <th>Raw.read.count.HV.18<\/th>\n      <th>Coverage.HV.19<\/th>\n      <th>TPM.HV.19<\/th>\n      <th>Raw.read.count.HV.19<\/th>\n      <th>Coverage.HV.20<\/th>\n      <th>TPM.HV.20<\/th>\n      <th>Raw.read.count.HV.20<\/th>\n      <th>Coverage.HV.21<\/th>\n      <th>TPM.HV.21<\/th>\n      <th>Raw.read.count.HV.21<\/th>\n      <th>Coverage.HV.22<\/th>\n      <th>TPM.HV.22<\/th>\n      <th>Raw.read.count.HV.22<\/th>\n      <th>Coverage.HV.23<\/th>\n      <th>TPM.HV.23<\/th>\n      <th>Raw.read.count.HV.23<\/th>\n      <th>Coverage.HV.24<\/th>\n      <th>TPM.HV.24<\/th>\n      <th>Raw.read.count.HV.24<\/th>\n    <\/tr>\n  <\/thead>\n<\/table>","options":{"columnDefs":[{"className":"dt-right","targets":[2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58,59,60,61,62,63,64,65,66,67,68,69,70,71,72,73,74,75,76,77,78,79,80,81,82,83,84,85,86,87,88,89,90,91,92,93,94,95,96,97,98,99,100,101,102,103,104,105,106,107,108,109,110,111,112,113,114,115,116,117,118,119,120,121,122,123,124,125,126,127,128,129,130,131,132,133,134,135,136,137,138,139,140,141,142,143,144,145,146,147,148]},{"orderable":false,"targets":0}],"order":[],"autoWidth":false,"orderClasses":false}},"evals":[],"jsHooks":[]}</script>
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
  select(ORF_ID, SQM_Length.NT, starts_with("ORF_Raw.read.count.")) %>% 
  pivot_longer(cols = starts_with("ORF_Raw.read.count."),
               values_to = "read_counts", 
               names_to = "sample") %>% 
  mutate(Num_Gi = read_counts / SQM_Length.NT) %>% 
  select(-read_counts) %>% 
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
<div id="htmlwidget-713a47829faa130e9c47" style="width:100%;height:auto;" class="datatables html-widget"></div>
<script type="application/json" data-for="htmlwidget-713a47829faa130e9c47">{"x":{"filter":"none","vertical":false,"data":[["1","2","3","4","5","6"],["h_megahit_1_1-312","h_megahit_2_3-176","h_megahit_2_173-304","h_megahit_3_3-116","h_megahit_3_82-219","h_megahit_4_1-72"],[0.0224358974358974,0,0,0,0,0],[0,0.0114942528735632,0.00757575757575758,0,0,0],[0,0,0,0.0175438596491228,0.0144927536231884,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0.0114942528735632,0.0151515151515152,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0.0277777777777778],[0,0,0,0.0087719298245614,0.0072463768115942,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0.0087719298245614,0.0072463768115942,0],[0,0,0,0.0087719298245614,0.0072463768115942,0],[0,0,0,0,0,0],[0,0,0,0.0087719298245614,0.0072463768115942,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0.0277777777777778],[0,0,0,0,0,0],[0,0,0,0.0087719298245614,0.0072463768115942,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0.0263157894736842,0.0217391304347826,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0]],"container":"<table class=\"display\">\n  <thead>\n    <tr>\n      <th> <\/th>\n      <th>ORF_ID<\/th>\n      <th>ORF_Num_Gi.D.1<\/th>\n      <th>ORF_Num_Gi.HC.1<\/th>\n      <th>ORF_Num_Gi.HC.2<\/th>\n      <th>ORF_Num_Gi.HC.3<\/th>\n      <th>ORF_Num_Gi.HC.4<\/th>\n      <th>ORF_Num_Gi.HC.5<\/th>\n      <th>ORF_Num_Gi.HC.6<\/th>\n      <th>ORF_Num_Gi.HC.7<\/th>\n      <th>ORF_Num_Gi.HC.8<\/th>\n      <th>ORF_Num_Gi.HC.9<\/th>\n      <th>ORF_Num_Gi.HC.10<\/th>\n      <th>ORF_Num_Gi.HC.11<\/th>\n      <th>ORF_Num_Gi.HC.12<\/th>\n      <th>ORF_Num_Gi.HC.13<\/th>\n      <th>ORF_Num_Gi.HC.14<\/th>\n      <th>ORF_Num_Gi.HC.15<\/th>\n      <th>ORF_Num_Gi.HC.16<\/th>\n      <th>ORF_Num_Gi.HC.17<\/th>\n      <th>ORF_Num_Gi.HC.18<\/th>\n      <th>ORF_Num_Gi.HC.19<\/th>\n      <th>ORF_Num_Gi.HC.20<\/th>\n      <th>ORF_Num_Gi.HC.21<\/th>\n      <th>ORF_Num_Gi.HC.22<\/th>\n      <th>ORF_Num_Gi.HC.23<\/th>\n      <th>ORF_Num_Gi.HC.24<\/th>\n      <th>ORF_Num_Gi.HV.1<\/th>\n      <th>ORF_Num_Gi.HV.2<\/th>\n      <th>ORF_Num_Gi.HV.3<\/th>\n      <th>ORF_Num_Gi.HV.4<\/th>\n      <th>ORF_Num_Gi.HV.5<\/th>\n      <th>ORF_Num_Gi.HV.6<\/th>\n      <th>ORF_Num_Gi.HV.7<\/th>\n      <th>ORF_Num_Gi.HV.8<\/th>\n      <th>ORF_Num_Gi.HV.9<\/th>\n      <th>ORF_Num_Gi.HV.10<\/th>\n      <th>ORF_Num_Gi.HV.11<\/th>\n      <th>ORF_Num_Gi.HV.12<\/th>\n      <th>ORF_Num_Gi.HV.13<\/th>\n      <th>ORF_Num_Gi.HV.14<\/th>\n      <th>ORF_Num_Gi.HV.15<\/th>\n      <th>ORF_Num_Gi.HV.16<\/th>\n      <th>ORF_Num_Gi.HV.17<\/th>\n      <th>ORF_Num_Gi.HV.18<\/th>\n      <th>ORF_Num_Gi.HV.19<\/th>\n      <th>ORF_Num_Gi.HV.20<\/th>\n      <th>ORF_Num_Gi.HV.21<\/th>\n      <th>ORF_Num_Gi.HV.22<\/th>\n      <th>ORF_Num_Gi.HV.23<\/th>\n      <th>ORF_Num_Gi.HV.24<\/th>\n    <\/tr>\n  <\/thead>\n<\/table>","options":{"columnDefs":[{"className":"dt-right","targets":[2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50]},{"orderable":false,"targets":0}],"order":[],"autoWidth":false,"orderClasses":false}},"evals":[],"jsHooks":[]}</script>
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
<div id="htmlwidget-3fe68c7ace0ec6d18dc9" style="width:100%;height:auto;" class="datatables html-widget"></div>
<script type="application/json" data-for="htmlwidget-3fe68c7ace0ec6d18dc9">{"x":{"filter":"none","vertical":false,"data":[["1","2","3","4","5","6"],["h_megahit_1_1-312","h_megahit_2_3-176","h_megahit_2_173-304","h_megahit_3_3-116","h_megahit_3_82-219","h_megahit_4_1-72"],[1.65106252977147e-07,0,0,0,0,0],[0,2.28833590277493e-07,1.50822139046529e-07,0,0,0],[0,0,0,2.97590210493958e-07,2.45835391277617e-07,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,2.30437459025747e-07,3.03758468715757e-07,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,4.04782802479438e-07],[0,0,0,1.56397247893644e-07,1.29197726520837e-07,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,1.78217802413832e-07,1.47223401994035e-07,0],[0,0,0,1.40701542858223e-07,1.16231709317663e-07,0],[0,0,0,0,0,0],[0,0,0,1.5402926080401e-07,1.27241563272878e-07,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,5.52850535761806e-07],[0,0,0,0,0,0],[0,0,0,1.65488750877943e-07,1.36708098551344e-07,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,4.33713887036034e-07,3.58285384942811e-07,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0]],"container":"<table class=\"display\">\n  <thead>\n    <tr>\n      <th> <\/th>\n      <th>ORF_ID<\/th>\n      <th>ORF_Num_Gi_pc.D.1<\/th>\n      <th>ORF_Num_Gi_pc.HC.1<\/th>\n      <th>ORF_Num_Gi_pc.HC.2<\/th>\n      <th>ORF_Num_Gi_pc.HC.3<\/th>\n      <th>ORF_Num_Gi_pc.HC.4<\/th>\n      <th>ORF_Num_Gi_pc.HC.5<\/th>\n      <th>ORF_Num_Gi_pc.HC.6<\/th>\n      <th>ORF_Num_Gi_pc.HC.7<\/th>\n      <th>ORF_Num_Gi_pc.HC.8<\/th>\n      <th>ORF_Num_Gi_pc.HC.9<\/th>\n      <th>ORF_Num_Gi_pc.HC.10<\/th>\n      <th>ORF_Num_Gi_pc.HC.11<\/th>\n      <th>ORF_Num_Gi_pc.HC.12<\/th>\n      <th>ORF_Num_Gi_pc.HC.13<\/th>\n      <th>ORF_Num_Gi_pc.HC.14<\/th>\n      <th>ORF_Num_Gi_pc.HC.15<\/th>\n      <th>ORF_Num_Gi_pc.HC.16<\/th>\n      <th>ORF_Num_Gi_pc.HC.17<\/th>\n      <th>ORF_Num_Gi_pc.HC.18<\/th>\n      <th>ORF_Num_Gi_pc.HC.19<\/th>\n      <th>ORF_Num_Gi_pc.HC.20<\/th>\n      <th>ORF_Num_Gi_pc.HC.21<\/th>\n      <th>ORF_Num_Gi_pc.HC.22<\/th>\n      <th>ORF_Num_Gi_pc.HC.23<\/th>\n      <th>ORF_Num_Gi_pc.HC.24<\/th>\n      <th>ORF_Num_Gi_pc.HV.1<\/th>\n      <th>ORF_Num_Gi_pc.HV.2<\/th>\n      <th>ORF_Num_Gi_pc.HV.3<\/th>\n      <th>ORF_Num_Gi_pc.HV.4<\/th>\n      <th>ORF_Num_Gi_pc.HV.5<\/th>\n      <th>ORF_Num_Gi_pc.HV.6<\/th>\n      <th>ORF_Num_Gi_pc.HV.7<\/th>\n      <th>ORF_Num_Gi_pc.HV.8<\/th>\n      <th>ORF_Num_Gi_pc.HV.9<\/th>\n      <th>ORF_Num_Gi_pc.HV.10<\/th>\n      <th>ORF_Num_Gi_pc.HV.11<\/th>\n      <th>ORF_Num_Gi_pc.HV.12<\/th>\n      <th>ORF_Num_Gi_pc.HV.13<\/th>\n      <th>ORF_Num_Gi_pc.HV.14<\/th>\n      <th>ORF_Num_Gi_pc.HV.15<\/th>\n      <th>ORF_Num_Gi_pc.HV.16<\/th>\n      <th>ORF_Num_Gi_pc.HV.17<\/th>\n      <th>ORF_Num_Gi_pc.HV.18<\/th>\n      <th>ORF_Num_Gi_pc.HV.19<\/th>\n      <th>ORF_Num_Gi_pc.HV.20<\/th>\n      <th>ORF_Num_Gi_pc.HV.21<\/th>\n      <th>ORF_Num_Gi_pc.HV.22<\/th>\n      <th>ORF_Num_Gi_pc.HV.23<\/th>\n      <th>ORF_Num_Gi_pc.HV.24<\/th>\n    <\/tr>\n  <\/thead>\n<\/table>","options":{"columnDefs":[{"className":"dt-right","targets":[2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50]},{"orderable":false,"targets":0}],"order":[],"autoWidth":false,"orderClasses":false}},"evals":[],"jsHooks":[]}</script>
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
##            sample  input filtered denoisedF denoisedR merged tabled filtered_pc
## HC_1   TR1-61-S47 238630   238248    237846    237395 231166 231166       0.998
## HC_2  TR2-61-S196  24392    24363     24296     24257  22096  22096       0.999
## HC_3   TR3-61-S45  49196    49111     49012     48940  47064  47064       0.998
## HC_4   TR4-61-S11  28442    28406     28348     28277  27014  27014       0.999
## HC_5  TR5-61-S142    245      244       240       237    219    219       0.996
## HC_6  TR6-61-S186 252856   252414    252020    251515 245450 245450       0.998
## HC_7  TR1-63-S100     32       31        31        30     26     26       0.969
## HC_8  TR2-63-S160  36811    36755     36712     36601  35561  35561       0.998
## HC_9  TR3-63-S206  26658    26606     26549     26501  25740  25740       0.998
## HC_10  TR4-63-S74  41734    41656     41517     41407  39936  39936       0.998
## HC_11  TR5-63-S82  43222    43149     43078     42988  41415  41415       0.998
## HC_12 TR6-63-S195  34779    34728     34672     34605  33178  33178       0.999
## HC_13 TR1-65-S128     60       59        53        51     36     36       0.983
## HC_14 TR2-65-S158  43440    43349     43294     43167  41820  41820       0.998
## HC_15  TR3-65-S92  32513    32447     32308     32211  29961  29961       0.998
## HC_16 TR4-65-S155     27       26        25        24     15     15       0.963
## HC_17 TR5-65-S187     35       35        32        32     25     25       1.000
## HC_18  TR6-65-S94  33692    33627     33496     33425  30776  30776       0.998
## HC_19 TR1-69-S145  40431    40343     40200     40139  38039  38039       0.998
## HC_20  TR2-69-S57 178654   178284    177932    177646 174558 174558       0.998
## HC_21 TR3-69-S139  30574    30520     30460     30405  29510  29510       0.998
## HC_22 TR4-69-S197  25982    25946     25869     25853  24646  24646       0.999
## HC_23 TR5-69-S107     23       23        22        22     19     19       1.000
## HC_24 TR6-69-S174  44451    44372     44282     44188  42841  42841       0.998
## HV_1   TR1-28-S35  34188    34109     34047     34001  32656  32656       0.998
## HV_2   TR2-28-S64  25719    25668     25602     25529  23955  23955       0.998
## HV_3   TR3-28-S38 833850   832269    830661    829156 814071 814071       0.998
## HV_4  TR4-28-S138    119      119       112       111     94     94       1.000
## HV_5    TR5-28-S5  31023    30972     30868     30821  29212  29212       0.998
## HV_6   TR6-28-S80  28872    28808     28729     28664  27021  27021       0.998
## HV_7   TR1-30-S76  36509    36431     36343     36255  35484  35484       0.998
## HV_8   TR2-30-S73  38226    38119     37993     37920  36739  36739       0.997
## HV_9  TR3-30-S190  35109    35014     34904     34824  33704  33704       0.997
## HV_10  TR4-30-S48  26568    26495     26335     26274  24097  24097       0.997
## HV_11  TR5-30-S62  33028    32966     32872     32818  31102  31102       0.998
## HV_12 TR6-30-S144  31410    31362     31309     31213  29786  29786       0.998
## HV_13 TR1-32-S111     32       32        30        27     23     23       1.000
## HV_14  TR2-32-S90  33510    33387     33329     33287  32683  32683       0.996
## HV_15 TR3-32-S211  24657    24582     24552     24494  24072  24072       0.997
## HV_16  TR4-32-S17  44815    44659     44563     44486  43879  43879       0.997
## HV_17  TR5-32-S56  23120    23079     23027     22987  21785  21785       0.998
## HV_18  TR6-32-S39  30398    30364     30310     30240  28866  28866       0.999
## D_1      D-2-S164  29089    29028     28917     28759  27917  27917       0.998
## HV_19 TR1-36-S198  34485    34334     34287     33498  28774  28774       0.996
## HV_20 TR2-36-S183  46398    46202     46139     45199  39386  39386       0.996
## HV_21   TR3-36-S8  19943    19861     19755     19696  18108  18108       0.996
## HV_22  TR4-36-S12  31375    31261     31220     31194  30717  30717       0.996
## HV_23  TR5-36-S85  33570    33504     33364     33301  31257  31257       0.998
## HV_24 TR6-36-S124     37       37        34        32     27     27       1.000
##       denoisedF_pc denoisedR_pc merged_pc filtered_merged_pc input_merged_pc
## HC_1         0.998        0.996     0.972              0.970           0.969
## HC_2         0.997        0.996     0.909              0.907           0.906
## HC_3         0.998        0.997     0.960              0.958           0.957
## HC_4         0.998        0.995     0.953              0.951           0.950
## HC_5         0.984        0.971     0.912              0.898           0.894
## HC_6         0.998        0.996     0.974              0.972           0.971
## HC_7         1.000        0.968     0.839              0.839           0.812
## HC_8         0.999        0.996     0.969              0.968           0.966
## HC_9         0.998        0.996     0.970              0.967           0.966
## HC_10        0.997        0.994     0.962              0.959           0.957
## HC_11        0.998        0.996     0.961              0.960           0.958
## HC_12        0.998        0.996     0.957              0.955           0.954
## HC_13        0.898        0.864     0.679              0.610           0.600
## HC_14        0.999        0.996     0.966              0.965           0.963
## HC_15        0.996        0.993     0.927              0.923           0.922
## HC_16        0.962        0.923     0.600              0.577           0.556
## HC_17        0.914        0.914     0.781              0.714           0.714
## HC_18        0.996        0.994     0.919              0.915           0.913
## HC_19        0.996        0.995     0.946              0.943           0.941
## HC_20        0.998        0.996     0.981              0.979           0.977
## HC_21        0.998        0.996     0.969              0.967           0.965
## HC_22        0.997        0.996     0.953              0.950           0.949
## HC_23        0.957        0.957     0.864              0.826           0.826
## HC_24        0.998        0.996     0.967              0.965           0.964
## HV_1         0.998        0.997     0.959              0.957           0.955
## HV_2         0.997        0.995     0.936              0.933           0.931
## HV_3         0.998        0.996     0.980              0.978           0.976
## HV_4         0.941        0.933     0.839              0.790           0.790
## HV_5         0.997        0.995     0.946              0.943           0.942
## HV_6         0.997        0.995     0.941              0.938           0.936
## HV_7         0.998        0.995     0.976              0.974           0.972
## HV_8         0.997        0.995     0.967              0.964           0.961
## HV_9         0.997        0.995     0.966              0.963           0.960
## HV_10        0.994        0.992     0.915              0.909           0.907
## HV_11        0.997        0.996     0.946              0.943           0.942
## HV_12        0.998        0.995     0.951              0.950           0.948
## HV_13        0.938        0.844     0.767              0.719           0.719
## HV_14        0.998        0.997     0.981              0.979           0.975
## HV_15        0.999        0.996     0.980              0.979           0.976
## HV_16        0.998        0.996     0.985              0.983           0.979
## HV_17        0.998        0.996     0.946              0.944           0.942
## HV_18        0.998        0.996     0.952              0.951           0.950
## D_1          0.996        0.991     0.965              0.962           0.960
## HV_19        0.999        0.976     0.839              0.838           0.834
## HV_20        0.999        0.978     0.854              0.852           0.849
## HV_21        0.995        0.992     0.917              0.912           0.908
## HV_22        0.999        0.998     0.984              0.983           0.979
## HV_23        0.996        0.994     0.937              0.933           0.931
## HV_24        0.919        0.865     0.794              0.730           0.730
##       tabled_joined chimera_out length_filtered tabled_pc chimera_out_pc
## HC_1         231166      188953          188948         1           0.82
## HC_2          22096       10235           10235         1           0.46
## HC_3          47064       36864           36862         1           0.78
## HC_4          27014       19750           19745         1           0.73
## HC_5            219         183             183         1           0.84
## HC_6         245450      201757          201745         1           0.82
## HC_7             26          24              24         1           0.92
## HC_8          35561       28863           28862         1           0.81
## HC_9          25740       20454           20454         1           0.79
## HC_10         39936       32773           32768         1           0.82
## HC_11         41415       31488           31486         1           0.76
## HC_12         33178       23813           23807         1           0.72
## HC_13            36          29              29         1           0.81
## HC_14         41820       32880           32880         1           0.79
## HC_15         29961       17134           17133         1           0.57
## HC_16            15          14              14         1           0.93
## HC_17            25          23              23         1           0.92
## HC_18         30776       15573           15573         1           0.51
## HC_19         38039       24528           24524         1           0.64
## HC_20        174558      149775          149775         1           0.86
## HC_21         29510       23544           23544         1           0.80
## HC_22         24646       15361           15359         1           0.62
## HC_23            19          17              17         1           0.89
## HC_24         42841       35074           35074         1           0.82
## HV_1          32656       24739           24739         1           0.76
## HV_2          23955       14713           14713         1           0.61
## HV_3         814071      695560          695530         1           0.85
## HV_4             94          83              83         1           0.88
## HV_5          29212       18666           18666         1           0.64
## HV_6          27021       17152           17151         1           0.63
## HV_7          35484       28461           28461         1           0.80
## HV_8          36739       25840           25840         1           0.70
## HV_9          33704       26093           26093         1           0.77
## HV_10         24097       11843           11843         1           0.49
## HV_11         31102       21287           21287         1           0.68
## HV_12         29786       21835           21832         1           0.73
## HV_13            23          18              18         1           0.78
## HV_14         32683       25240           25240         1           0.77
## HV_15         24072       19576           19576         1           0.81
## HV_16         43879       34419           34419         1           0.78
## HV_17         21785       15674           15673         1           0.72
## HV_18         28866       21084           21084         1           0.73
## D_1           27917       26104           26104         1           0.94
## HV_19         28774       22660           22660         1           0.79
## HV_20         39386       33722           33722         1           0.86
## HV_21         18108        8748            8748         1           0.48
## HV_22         30717       23394           23394         1           0.76
## HV_23         31257       19614           19614         1           0.63
## HV_24            27          20              20         1           0.74
##       length_filtered_pc Sample_description I7_Index_ID index I5_Index_ID
## HC_1                   1             TR1-61      A-N706  <NA>      A-S510
## HC_2                   1             TR2-61      D-N716  <NA>      A-S506
## HC_3                   1             TR3-61      A-N706  <NA>      A-S507
## HC_4                   1             TR4-61      A-N702  <NA>      A-S505
## HC_5                   1             TR5-61      A-N706  <NA>      D-S520
## HC_6                   1             TR6-61      A-N715  <NA>      D-S515
## HC_7                   1             TR1-63      A-N701  <NA>      D-S517
## HC_8                   1             TR2-63      A-N710  <NA>      D-S522
## HC_9                   1             TR3-63      D-N718  <NA>      A-S508
## HC_10                  1             TR4-63      A-N712  <NA>      A-S503
## HC_11                  1             TR5-63      A-N714  <NA>      A-S503
## HC_12                  1             TR6-63      D-N716  <NA>      A-S505
## HC_13                  1             TR1-65      A-N704  <NA>      D-S522
## HC_14                  1             TR2-65      A-N710  <NA>      D-S520
## HC_15                  1             TR3-65      A-N715  <NA>      A-S506
## HC_16                  1             TR4-65      A-N710  <NA>      D-S516
## HC_17                  1             TR5-65      A-N715  <NA>      D-S516
## HC_18                  1             TR6-65      A-N715  <NA>      A-S508
## HC_19                  1             TR1-69      A-N707  <NA>      D-S513
## HC_20                  1             TR2-69      A-N710  <NA>      A-S502
## HC_21                  1             TR3-69      A-N706  <NA>      D-S516
## HC_22                  1             TR4-69      D-N716  <NA>      A-S507
## HC_23                  1             TR5-69      A-N702  <NA>      D-S516
## HC_24                  1             TR6-69      A-N712  <NA>      D-S520
## HV_1                   1             TR1-28      A-N705  <NA>      A-S505
## HV_2                   1             TR2-28      A-N710  <NA>      A-S511
## HV_3                   1             TR3-28      A-N705  <NA>      A-S508
## HV_4                   1             TR4-28      A-N706  <NA>      D-S515
## HV_5                   1             TR5-28      A-N701  <NA>      A-S507
## HV_6                   1             TR6-28      A-N712  <NA>      A-S511
## HV_7                   1             TR1-30      A-N712  <NA>      A-S506
## HV_8                   1             TR2-30      A-N712  <NA>      A-S502
## HV_9                   1             TR3-30      A-N715  <NA>      D-S520
## HV_10                  1             TR4-30      A-N706  <NA>      A-S511
## HV_11                  1             TR5-30      A-N710  <NA>      A-S508
## HV_12                  1             TR6-30      A-N706  <NA>      D-S522
## HV_13                  1             TR1-32      A-N702  <NA>      D-S521
## HV_14                  1             TR2-32      A-N715  <NA>      A-S503
## HV_15                  1             TR3-32      D-N719  <NA>      A-S505
## HV_16                  1             TR4-32      A-N703  <NA>      A-S502
## HV_17                  1             TR5-32      A-N707  <NA>      A-S511
## HV_18                  1             TR6-32      A-N705  <NA>      A-S510
## D_1                    1                D-2      A-N711  <NA>      D-S517
## HV_19                  1             TR1-36      D-N716  <NA>      A-S508
## HV_20                  1             TR2-36      A-N714  <NA>      D-S521
## HV_21                  1             TR3-36      A-N701  <NA>      A-S511
## HV_22                  1             TR4-36      A-N702  <NA>      A-S506
## HV_23                  1             TR5-36      A-N714  <NA>      A-S507
## HV_24                  1             TR6-36      A-N704  <NA>      D-S517
##       index2 Description2 Experiment Reactor     Treatment Day_of_Connection
## HC_1    <NA>         <NA> Continuous     TR1   CTX+HV292.1                 8
## HC_2    <NA>         <NA> Continuous     TR2           CTX                 8
## HC_3    <NA>         <NA> Continuous     TR3   CTX+HV292.1                 8
## HC_4    <NA>         <NA> Continuous     TR4           CTX                 8
## HC_5    <NA>         <NA> Continuous     TR5       HV292.1                 8
## HC_6    <NA>         <NA> Continuous      CR     UNTREATED                 8
## HC_7    <NA>         <NA> Continuous     TR1   CTX+HV292.1                10
## HC_8    <NA>         <NA> Continuous     TR2           CTX                10
## HC_9    <NA>         <NA> Continuous     TR3   CTX+HV292.1                10
## HC_10   <NA>         <NA> Continuous     TR4           CTX                10
## HC_11   <NA>         <NA> Continuous     TR5       HV292.1                10
## HC_12   <NA>         <NA> Continuous      CR     UNTREATED                10
## HC_13   <NA>         <NA> Continuous     TR1   CTX+HV292.1                12
## HC_14   <NA>         <NA> Continuous     TR2           CTX                12
## HC_15   <NA>         <NA> Continuous     TR3   CTX+HV292.1                12
## HC_16   <NA>         <NA> Continuous     TR4           CTX                12
## HC_17   <NA>         <NA> Continuous     TR5       HV292.1                12
## HC_18   <NA>         <NA> Continuous      CR     UNTREATED                12
## HC_19   <NA>         <NA> Continuous     TR1   CTX+HV292.1                16
## HC_20   <NA>         <NA> Continuous     TR2           CTX                16
## HC_21   <NA>         <NA> Continuous     TR3   CTX+HV292.1                16
## HC_22   <NA>         <NA> Continuous     TR4           CTX                16
## HC_23   <NA>         <NA> Continuous     TR5       HV292.1                16
## HC_24   <NA>         <NA> Continuous      CR     UNTREATED                16
## HV_1    <NA>         <NA> Continuous     TR1 VAN+CCUG59168                28
## HV_2    <NA>         <NA> Continuous     TR2           VAN                28
## HV_3    <NA>         <NA> Continuous     TR3 VAN+CCUG59168                28
## HV_4    <NA>         <NA> Continuous     TR4           VAN                28
## HV_5    <NA>         <NA> Continuous     TR5     CCUG59168                28
## HV_6    <NA>         <NA> Continuous      CR     UNTREATED                28
## HV_7    <NA>         <NA> Continuous     TR1 VAN+CCUG59168                30
## HV_8    <NA>         <NA> Continuous     TR2           VAN                30
## HV_9    <NA>         <NA> Continuous     TR3 VAN+CCUG59168                30
## HV_10   <NA>         <NA> Continuous     TR4           VAN                30
## HV_11   <NA>         <NA> Continuous     TR5     CCUG59168                30
## HV_12   <NA>         <NA> Continuous      CR     UNTREATED                30
## HV_13   <NA>         <NA> Continuous     TR1 VAN+CCUG59168                32
## HV_14   <NA>         <NA> Continuous     TR2           VAN                32
## HV_15   <NA>         <NA> Continuous     TR3 VAN+CCUG59168                32
## HV_16   <NA>         <NA> Continuous     TR4           VAN                32
## HV_17   <NA>         <NA> Continuous     TR5     CCUG59168                32
## HV_18   <NA>         <NA> Continuous      CR     UNTREATED                32
## D_1     <NA>        DONOR      Cecum   DONOR         DONOR                NA
## HV_19   <NA>         <NA> Continuous     TR1 VAN+CCUG59168                36
## HV_20   <NA>         <NA> Continuous     TR2           VAN                36
## HV_21   <NA>         <NA> Continuous     TR3 VAN+CCUG59168                36
## HV_22   <NA>         <NA> Continuous     TR4           VAN                36
## HV_23   <NA>         <NA> Continuous     TR5     CCUG59168                36
## HV_24   <NA>         <NA> Continuous      CR     UNTREATED                36
##       Day_of_Treatment Day_from_Inoculum  Enrichment Phase    Treatment2 Date
## HC_1                -1                76 NotEnriched  Stab    AB+E. coli <NA>
## HC_2                -1                76 NotEnriched  Stab            AB <NA>
## HC_3                -1                76 NotEnriched  Stab    AB+E. coli <NA>
## HC_4                -1                76 NotEnriched  Stab            AB <NA>
## HC_5                -1                76 NotEnriched  Stab       E. coli <NA>
## HC_6                -1                76 NotEnriched  Stab     UNTREATED <NA>
## HC_7                 1                78 NotEnriched Treat    AB+E. coli <NA>
## HC_8                 1                78 NotEnriched Treat            AB <NA>
## HC_9                 1                78 NotEnriched Treat    AB+E. coli <NA>
## HC_10                1                78 NotEnriched Treat            AB <NA>
## HC_11                1                78 NotEnriched Treat       E. coli <NA>
## HC_12                1                78 NotEnriched Treat     UNTREATED <NA>
## HC_13                3                80 NotEnriched Treat    AB+E. coli <NA>
## HC_14                3                80 NotEnriched Treat            AB <NA>
## HC_15                3                80 NotEnriched Treat    AB+E. coli <NA>
## HC_16                3                80 NotEnriched Treat            AB <NA>
## HC_17                3                80 NotEnriched Treat       E. coli <NA>
## HC_18                3                80 NotEnriched Treat     UNTREATED <NA>
## HC_19                7                84 NotEnriched Treat    AB+E. coli <NA>
## HC_20                7                84 NotEnriched Treat            AB <NA>
## HC_21                7                84 NotEnriched Treat    AB+E. coli <NA>
## HC_22                7                84 NotEnriched Treat            AB <NA>
## HC_23                7                84 NotEnriched Treat       E. coli <NA>
## HC_24                7                84 NotEnriched Treat     UNTREATED <NA>
## HV_1                -1                43 NotEnriched  Stab AB+E. faecium <NA>
## HV_2                -1                43 NotEnriched  Stab            AB <NA>
## HV_3                -1                43 NotEnriched  Stab AB+E. faecium <NA>
## HV_4                -1                43 NotEnriched  Stab            AB <NA>
## HV_5                -1                43 NotEnriched  Stab    E. faecium <NA>
## HV_6                -1                43 NotEnriched  Stab     UNTREATED <NA>
## HV_7                 1                45 NotEnriched Treat AB+E. faecium <NA>
## HV_8                 1                45 NotEnriched Treat            AB <NA>
## HV_9                 1                45 NotEnriched Treat AB+E. faecium <NA>
## HV_10                1                45 NotEnriched Treat            AB <NA>
## HV_11                1                45 NotEnriched Treat    E. faecium <NA>
## HV_12                1                45 NotEnriched Treat     UNTREATED <NA>
## HV_13                3                47 NotEnriched Treat AB+E. faecium <NA>
## HV_14                3                47 NotEnriched Treat            AB <NA>
## HV_15                3                47 NotEnriched Treat AB+E. faecium <NA>
## HV_16                3                47 NotEnriched Treat            AB <NA>
## HV_17                3                47 NotEnriched Treat    E. faecium <NA>
## HV_18                3                47 NotEnriched Treat     UNTREATED <NA>
## D_1               <NA>                NA NotEnriched DONOR         DONOR <NA>
## HV_19                7                51 NotEnriched Treat AB+E. faecium <NA>
## HV_20                7                51 NotEnriched Treat            AB <NA>
## HV_21                7                51 NotEnriched Treat AB+E. faecium <NA>
## HV_22                7                51 NotEnriched Treat            AB <NA>
## HV_23                7                51 NotEnriched Treat    E. faecium <NA>
## HV_24                7                51 NotEnriched Treat     UNTREATED <NA>
##       Paul Reactor_Treatment GeneCopyNumberperML HV292.1_Copy_Number_permL
## HC_1  <NA>   TR1_CTX+HV292.1            2.14e+11                        NA
## HC_2  <NA>           TR2_CTX            2.07e+11                        NA
## HC_3  <NA>   TR3_CTX+HV292.1            1.98e+11                        NA
## HC_4  <NA>           TR4_CTX            2.05e+11                        NA
## HC_5  <NA>       TR5_HV292.1            2.43e+11                        NA
## HC_6  <NA>      CR_UNTREATED            2.37e+11                        NA
## HC_7  <NA>   TR1_CTX+HV292.1            5.59e+11                   5520000
## HC_8  <NA>           TR2_CTX            2.66e+11                         0
## HC_9  <NA>   TR3_CTX+HV292.1            3.83e+11                  30487500
## HC_10 <NA>           TR4_CTX            3.22e+11                         0
## HC_11 <NA>       TR5_HV292.1            4.52e+11                  20600000
## HC_12 <NA>      CR_UNTREATED            2.39e+11                         0
## HC_13 <NA>   TR1_CTX+HV292.1            1.36e+11                         0
## HC_14 <NA>           TR2_CTX            2.85e+11                         0
## HC_15 <NA>   TR3_CTX+HV292.1            2.04e+11                         0
## HC_16 <NA>           TR4_CTX            2.39e+11                         0
## HC_17 <NA>       TR5_HV292.1            2.54e+11                         0
## HC_18 <NA>      CR_UNTREATED            1.26e+11                         0
## HC_19 <NA>   TR1_CTX+HV292.1            2.10e+11                         0
## HC_20 <NA>           TR2_CTX            1.02e+11                         0
## HC_21 <NA>   TR3_CTX+HV292.1            2.41e+11                         0
## HC_22 <NA>           TR4_CTX            2.42e+11                         0
## HC_23 <NA>       TR5_HV292.1            3.93e+11                         0
## HC_24 <NA>      CR_UNTREATED            2.60e+11                         0
## HV_1  <NA> TR1_VAN+CCUG59168            1.11e+11                        NA
## HV_2  <NA>           TR2_VAN            2.50e+11                        NA
## HV_3  <NA> TR3_VAN+CCUG59168            3.13e+11                        NA
## HV_4  <NA>           TR4_VAN            2.99e+10                        NA
## HV_5  <NA>     TR5_CCUG59168            1.92e+11                        NA
## HV_6  <NA>      CR_UNTREATED            2.58e+11                        NA
## HV_7  <NA> TR1_VAN+CCUG59168            3.03e+11                        NA
## HV_8  <NA>           TR2_VAN            7.45e+10                        NA
## HV_9  <NA> TR3_VAN+CCUG59168            1.21e+11                        NA
## HV_10 <NA>           TR4_VAN            1.32e+11                        NA
## HV_11 <NA>     TR5_CCUG59168            2.77e+11                        NA
## HV_12 <NA>      CR_UNTREATED            4.92e+10                        NA
## HV_13 <NA> TR1_VAN+CCUG59168            1.94e+11                        NA
## HV_14 <NA>           TR2_VAN            2.43e+11                        NA
## HV_15 <NA> TR3_VAN+CCUG59168            2.00e+11                        NA
## HV_16 <NA>           TR4_VAN            9.63e+10                        NA
## HV_17 <NA>     TR5_CCUG59168            2.00e+11                        NA
## HV_18 <NA>      CR_UNTREATED            3.02e+11                        NA
## D_1   <NA>             DONOR            1.90e+10                        NA
## HV_19 <NA> TR1_VAN+CCUG59168            4.92e+10                        NA
## HV_20 <NA>           TR2_VAN            3.21e+11                        NA
## HV_21 <NA> TR3_VAN+CCUG59168            2.91e+11                        NA
## HV_22 <NA>           TR4_VAN            1.17e+11                        NA
## HV_23 <NA>     TR5_CCUG59168            1.68e+11                        NA
## HV_24 <NA>      CR_UNTREATED            4.39e+11                        NA
##       CCUG59168_Copy_Number_permL CTX_Copy_Number_permL VAN_Copy_Number_permL
## HC_1                           NA                    NA                    NA
## HC_2                           NA                    NA                    NA
## HC_3                           NA                    NA                    NA
## HC_4                           NA                    NA                    NA
## HC_5                           NA                    NA                    NA
## HC_6                           NA                    NA                    NA
## HC_7                           NA              11025000                    NA
## HC_8                           NA                     0                    NA
## HC_9                           NA              53737500                    NA
## HC_10                          NA                     0                    NA
## HC_11                          NA              41600000                    NA
## HC_12                          NA                     0                    NA
## HC_13                          NA               2160000                    NA
## HC_14                          NA                     0                    NA
## HC_15                          NA               3550000                    NA
## HC_16                          NA                     0                    NA
## HC_17                          NA                210000                    NA
## HC_18                          NA                     0                    NA
## HC_19                          NA                     0                    NA
## HC_20                          NA                     0                    NA
## HC_21                          NA               1117500                    NA
## HC_22                          NA                     0                    NA
## HC_23                          NA                    NA                    NA
## HC_24                          NA                     0                    NA
## HV_1                           NA                    NA                    NA
## HV_2                           NA                    NA                    NA
## HV_3                           NA                    NA                    NA
## HV_4                           NA                    NA                    NA
## HV_5                           NA                    NA                    NA
## HV_6                           NA                    NA                    NA
## HV_7                    1.380e+10                    NA             9.600e+12
## HV_8                    0.000e+00                    NA             0.000e+00
## HV_9                    2.080e+09                    NA             8.500e+11
## HV_10                   0.000e+00                    NA             0.000e+00
## HV_11                   1.570e+07                    NA             6.125e+09
## HV_12                   0.000e+00                    NA             0.000e+00
## HV_13                   8.450e+09                    NA             8.400e+12
## HV_14                   0.000e+00                    NA             0.000e+00
## HV_15                   2.690e+09                    NA             1.190e+12
## HV_16                   0.000e+00                    NA             0.000e+00
## HV_17                   8.830e+05                    NA             3.470e+08
## HV_18                   0.000e+00                    NA             0.000e+00
## D_1                            NA                    NA                    NA
## HV_19                   5.165e+09                    NA             1.300e+12
## HV_20                   0.000e+00                    NA             0.000e+00
## HV_21                   1.670e+09                    NA             4.800e+11
## HV_22                   0.000e+00                    NA             0.000e+00
## HV_23                   1.290e+05                    NA             0.000e+00
## HV_24                   0.000e+00                    NA             0.000e+00
##       Model Antibiotic_mg.mL Fermentation Antibiotic Lactose_mM Glucose_mM
## HC_1  Human               20            2        CTX      0.000      0.000
## HC_2  Human               20            2        CTX      0.000      0.000
## HC_3  Human              200            2        CTX      0.000      0.000
## HC_4  Human              200            2        CTX      0.000      0.000
## HC_5  Human               NA            2        CTX      0.000      0.000
## HC_6  Human               NA            2        CTX      0.000      0.000
## HC_7  Human               20            2        CTX      0.000      0.000
## HC_8  Human               20            2        CTX      0.000      0.000
## HC_9  Human              200            2        CTX      0.000      0.305
## HC_10 Human              200            2        CTX      0.000      0.309
## HC_11 Human               NA            2        CTX      0.000      0.000
## HC_12 Human               NA            2        CTX      0.000      0.000
## HC_13 Human               20            2        CTX      0.000      0.000
## HC_14 Human               20            2        CTX      0.000      0.000
## HC_15 Human              200            2        CTX      0.000      0.000
## HC_16 Human              200            2        CTX      0.000      0.000
## HC_17 Human               NA            2        CTX      0.000      0.000
## HC_18 Human               NA            2        CTX      0.000      0.000
## HC_19 Human               20            2        CTX      0.000      0.000
## HC_20 Human               20            2        CTX      0.000      0.000
## HC_21 Human              200            2        CTX      0.000      0.000
## HC_22 Human              200            2        CTX      0.000      0.000
## HC_23 Human               NA            2        CTX      0.000      0.000
## HC_24 Human               NA            2        CTX      0.000      0.000
## HV_1  Human               90            1        VAN      0.053      0.173
## HV_2  Human               90            1        VAN      0.000      0.228
## HV_3  Human              600            1        VAN      0.069      0.000
## HV_4  Human              600            1        VAN      0.000      0.000
## HV_5  Human               NA            1        VAN      0.037      0.284
## HV_6  Human               NA            1        VAN      0.046      0.000
## HV_7  Human               90            1        VAN      0.000      3.207
## HV_8  Human               90            1        VAN      0.000      9.781
## HV_9  Human              600            1        VAN      0.000      0.000
## HV_10 Human              600            1        VAN      0.000      6.682
## HV_11 Human               NA            1        VAN      0.000      0.000
## HV_12 Human               NA            1        VAN      0.000      0.000
## HV_13 Human               90            1        VAN      1.580      0.767
## HV_14 Human               90            1        VAN      1.238      1.616
## HV_15 Human              600            1        VAN      0.904      0.000
## HV_16 Human              600            1        VAN      1.663      0.000
## HV_17 Human               NA            1        VAN      0.000      0.000
## HV_18 Human               NA            1        VAN      0.000      0.000
## D_1   Human               NA           NA       <NA>         NA         NA
## HV_19 Human               90            1        VAN      0.777      0.244
## HV_20 Human               90            1        VAN      0.731      0.000
## HV_21 Human              600            1        VAN      1.220      0.000
## HV_22 Human              600            1        VAN      1.579      0.000
## HV_23 Human               NA            1        VAN      0.000      0.000
## HV_24 Human               NA            1        VAN      0.000      0.000
##       Galactose_mM Succinat_mM Lactat_mM Formiat_mM Acetat_mM Propionat_mM
## HC_1         0.536       1.519     1.338      5.898    68.571       27.257
## HC_2         0.588       2.608     1.287      8.273    71.535       26.884
## HC_3         0.655       1.531     1.310      7.826    76.633       27.355
## HC_4         0.622       1.513     1.238      4.607    70.801       25.104
## HC_5         0.602       1.448     1.321      7.221    73.189       36.735
## HC_6         0.703       1.356     1.423      6.459    52.167       24.281
## HC_7         0.314       1.558     1.613      6.574    75.484       32.194
## HC_8         0.703       1.497     1.750      9.178    71.548       29.822
## HC_9         0.989       1.200     3.095      6.354    69.382       34.611
## HC_10        0.958       1.300     3.283      0.000    64.353       34.868
## HC_11        0.680       1.435     1.293      6.330    70.102       38.008
## HC_12        0.464       1.390     1.453      5.867    52.980       25.113
## HC_13        0.367       1.370     1.552      3.816    70.605       26.932
## HC_14        0.000       1.085     1.505      5.782    68.421       23.025
## HC_15        0.391       1.100     2.160      5.931    65.980       27.327
## HC_16        0.397       1.051     2.611      0.000    64.822       36.192
## HC_17        0.000       1.110     1.029      5.727    67.573       39.064
## HC_18        0.000       1.258     0.905      5.953    54.688       26.113
## HC_19        0.229       1.194     1.038      3.661    69.252       25.772
## HC_20        0.759       1.121     2.450      0.000    59.247       25.199
## HC_21        0.000       1.739     1.127      5.850    66.606       24.514
## HC_22        0.000       1.796     1.322      0.000    58.881       33.442
## HC_23        0.000       1.730     1.331     10.916    67.988       35.881
## HC_24        0.000       1.784     0.853      4.163    61.137       27.614
## HV_1         0.000       0.901     0.534      5.278    77.645       23.407
## HV_2         0.000       0.886     0.734      6.244    78.499       27.767
## HV_3         0.000       0.844     0.679      7.344    57.376       36.797
## HV_4         0.000       0.932     0.449      7.921    69.627       22.762
## HV_5         0.000       1.135     0.364      3.277    76.749       19.707
## HV_6         0.000       0.819     0.000      8.203    64.866       17.193
## HV_7         4.096       1.351    29.016     10.029    23.392       19.437
## HV_8         4.362       1.558     5.440      2.891    21.269       19.543
## HV_9         2.502       1.962    31.934     10.365    19.888       17.952
## HV_10        4.270       1.587     5.702      2.731    29.523       14.073
## HV_11        0.000       1.239     1.215      5.294    78.115       17.219
## HV_12        0.000       1.246     0.994      9.137    71.280       19.527
## HV_13        2.044       1.473     9.844     13.503    11.513        7.533
## HV_14        1.323       1.637     5.033     12.418    20.830       25.674
## HV_15        2.563       1.559    15.308      8.918    11.689       13.911
## HV_16        2.668       1.479     7.130     14.424    12.683        6.561
## HV_17        0.000       1.233     1.269      6.174    73.792       17.652
## HV_18        0.000       1.251     1.113     10.938    74.408       23.448
## D_1             NA          NA        NA         NA        NA           NA
## HV_19        1.325       1.417     5.496      8.836    18.665       41.793
## HV_20        1.144       1.378     4.138      8.562    19.052       45.481
## HV_21        0.784       1.435    13.537     10.322    10.928       11.060
## HV_22        0.750       1.401     4.294     12.778    15.865       10.169
## HV_23        0.266       1.208     1.372      6.666    72.238       14.474
## HV_24        0.000       1.187     1.171      9.060    64.405       19.489
##       Isobutyrat_mM Butyrat_mM Isovalerat_mM Valerat_mM Total_SCFA_mM
## HC_1          4.937     48.394         3.842      3.761       166.053
## HC_2          4.776     45.421         3.394      3.482       168.248
## HC_3          4.521     44.309         3.614      1.935       169.689
## HC_4          4.911     50.198         4.062      1.715       164.771
## HC_5          3.001     33.900         2.809      1.264       161.490
## HC_6          3.506     53.009         3.828      1.687       148.419
## HC_7          4.251     41.199         3.894      4.031       171.112
## HC_8          4.000     43.675         2.980      3.486       168.639
## HC_9          3.285     36.085         3.394      1.704       160.404
## HC_10         3.117     39.773         3.610      1.147       152.718
## HC_11         3.627     35.459         3.258      1.282       161.474
## HC_12         3.912     52.049         3.952      1.797       148.977
## HC_13         4.857     39.501         4.427      5.098       158.525
## HC_14         4.640     42.896         3.259      4.630       155.243
## HC_15         3.912     40.000         5.518      1.598       153.917
## HC_16         3.756     38.253         5.307      1.038       153.427
## HC_17         3.187     35.241         3.844      1.038       157.813
## HC_18         3.488     51.621         4.236      1.783       150.045
## HC_19         5.347     47.850         5.290      6.032       165.665
## HC_20         5.255     49.825         6.343      4.494       154.693
## HC_21         5.216     48.333         4.836      6.599       164.820
## HC_22         5.939     48.455         7.871      1.967       159.673
## HC_23         0.000     38.423         4.747      1.182       162.198
## HC_24         4.315     48.676         5.126      1.616       155.284
## HV_1          3.886     40.871         4.754      1.239       158.741
## HV_2          4.418     37.598         5.287      1.303       162.964
## HV_3          4.541     33.022         5.700      0.000       146.372
## HV_4          3.449     43.182         2.876      0.475       151.673
## HV_5          4.269     42.214         4.334      0.284       152.654
## HV_6          2.603     44.874         3.593      0.390       142.587
## HV_7          2.571     14.661         4.272      1.059       113.091
## HV_8          3.698     16.553         6.617      1.054        92.766
## HV_9          3.586     18.726         5.648      1.069       113.632
## HV_10         2.032     12.535         2.750      1.147        83.032
## HV_11         4.725     46.529         4.203      1.259       159.798
## HV_12         2.548     44.810         3.415      1.293       154.250
## HV_13         3.121      9.105         5.593      1.046        67.122
## HV_14         3.302     12.666         5.969      1.079        92.785
## HV_15         3.508      7.647         5.601      1.023        72.631
## HV_16         2.260      5.816         4.042      0.958        59.684
## HV_17         5.184     48.061         4.140      1.393       158.898
## HV_18         2.422     42.201         2.336      1.288       159.405
## D_1              NA         NA            NA         NA            NA
## HV_19         4.202     13.490         6.903      1.239       104.387
## HV_20         4.842     12.833         7.904      1.243       107.308
## HV_21         2.549      5.173         3.751      0.973        61.732
## HV_22         2.208      4.905         4.026      0.971        58.946
## HV_23         4.565     48.449         3.707      1.426       154.371
## HV_24         2.321     42.328         2.300      1.239       143.500
##       raw_metagenomic_pairs Reactor_Treatment_Dose   Treatment_Dose Period
## HC_1                     NA      TR1_CTX+HV292.120    CTX+HV292.120   pret
## HC_2                     NA              TR2_CTX20            CTX20   pret
## HC_3                     NA     TR3_CTX+HV292.1200   CTX+HV292.1200   pret
## HC_4                     NA             TR4_CTX200           CTX200   pret
## HC_5                     NA            TR5_HV292.1          HV292.1   pret
## HC_6                     NA           CR_UNTREATED        UNTREATED   pret
## HC_7                     NA      TR1_CTX+HV292.120    CTX+HV292.120      1
## HC_8                     NA              TR2_CTX20            CTX20      1
## HC_9                     NA     TR3_CTX+HV292.1200   CTX+HV292.1200      1
## HC_10                    NA             TR4_CTX200           CTX200      1
## HC_11                    NA            TR5_HV292.1          HV292.1      1
## HC_12                    NA           CR_UNTREATED        UNTREATED      1
## HC_13                    NA      TR1_CTX+HV292.120    CTX+HV292.120      3
## HC_14                    NA              TR2_CTX20            CTX20      3
## HC_15                    NA     TR3_CTX+HV292.1200   CTX+HV292.1200      3
## HC_16                    NA             TR4_CTX200           CTX200      3
## HC_17                    NA            TR5_HV292.1          HV292.1      3
## HC_18                    NA           CR_UNTREATED        UNTREATED      3
## HC_19                    NA      TR1_CTX+HV292.120    CTX+HV292.120      7
## HC_20                    NA              TR2_CTX20            CTX20      7
## HC_21                    NA     TR3_CTX+HV292.1200   CTX+HV292.1200      7
## HC_22                    NA             TR4_CTX200           CTX200      7
## HC_23                    NA            TR5_HV292.1          HV292.1      7
## HC_24                    NA           CR_UNTREATED        UNTREATED      7
## HV_1                     NA    TR1_VAN+CCUG5916890  VAN+CCUG5916890   pret
## HV_2                     NA              TR2_VAN90            VAN90   pret
## HV_3                     NA   TR3_VAN+CCUG59168600 VAN+CCUG59168600   pret
## HV_4                     NA             TR4_VAN600           VAN600   pret
## HV_5                     NA          TR5_CCUG59168        CCUG59168   pret
## HV_6                     NA           CR_UNTREATED        UNTREATED   pret
## HV_7                     NA    TR1_VAN+CCUG5916890  VAN+CCUG5916890      1
## HV_8                     NA              TR2_VAN90            VAN90      1
## HV_9                     NA   TR3_VAN+CCUG59168600 VAN+CCUG59168600      1
## HV_10                    NA             TR4_VAN600           VAN600      1
## HV_11                    NA          TR5_CCUG59168        CCUG59168      1
## HV_12                    NA           CR_UNTREATED        UNTREATED      1
## HV_13                    NA    TR1_VAN+CCUG5916890  VAN+CCUG5916890      3
## HV_14                    NA              TR2_VAN90            VAN90      3
## HV_15                    NA   TR3_VAN+CCUG59168600 VAN+CCUG59168600      3
## HV_16                    NA             TR4_VAN600           VAN600      3
## HV_17                    NA          TR5_CCUG59168        CCUG59168      3
## HV_18                    NA           CR_UNTREATED        UNTREATED      3
## D_1                      NA                  DONOR            DONOR   pret
## HV_19                    NA    TR1_VAN+CCUG5916890  VAN+CCUG5916890      7
## HV_20                    NA              TR2_VAN90            VAN90      7
## HV_21                    NA   TR3_VAN+CCUG59168600 VAN+CCUG59168600      7
## HV_22                    NA             TR4_VAN600           VAN600      7
## HV_23                    NA          TR5_CCUG59168        CCUG59168      7
## HV_24                    NA           CR_UNTREATED        UNTREATED      7
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
  select(ORF_TPM.D.1:ORF_TPM.HV.24) %>% 
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
## sample_data() Sample Data:       [ 49 samples by 63 sample variables ]
## tax_table()   Taxonomy Table:    [ 791180 taxa by 92 taxonomic ranks ]
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
## sample_data() Sample Data:       [ 49 samples by 63 sample variables ]
## tax_table()   Taxonomy Table:    [ 791180 taxa by 92 taxonomic ranks ]
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
## sample_data() Sample Data:       [ 49 samples by 63 sample variables ]
## tax_table()   Taxonomy Table:    [ 791180 taxa by 92 taxonomic ranks ]
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
## sample_data() Sample Data:       [ 49 samples by 63 sample variables ]
## tax_table()   Taxonomy Table:    [ 791180 taxa by 92 taxonomic ranks ]
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
## sample_data() Sample Data:       [ 49 samples by 63 sample variables ]
## tax_table()   Taxonomy Table:    [ 791180 taxa by 92 taxonomic ranks ]
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
  saveRDS(here::here("data/processed/human_full_gene_catalog_phyloseq.RDS"))

full_gene_catalogue_orf_quant_all_metrics %>% 
  saveRDS(here::here("data/processed/human_full_gene_catalog_full_metrics_table.RDS"))

full_gene_catalogue_orf_quant_all_metrics %>% 
  write_tsv(here::here("data/processed/human_full_gene_catalog_full_metrics_table.tsv.gz"))
```


Export only AMR and Tox:


```r
full_gene_catalogue_orf_quant_all_metrics %>% 
  filter(!is.na(PathoFact_AMR_ARG) | !is.na(rgi_CARD_ARO) | !is.na(PathoFact_Tox_Toxin_classification)) %>% 
  write_tsv(here::here("data/processed/human_full_gene_catalog_full_metrics_table_AMR_tox.tsv"))
```



We might need to clarify some cases before saving: e.g.: plasmid but in MAG...



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
##  date     2022-06-27
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
##  reshape2         * 1.4.4      2020-04-09 [1] CRAN (R 4.2.0)
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
##  scales           * 1.2.0      2022-04-13 [1] CRAN (R 4.2.0)
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

