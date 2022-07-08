---
title: "NRP72 -metadata - humann - chiken "
author: "Florentin Constancias"
date: "July 06, 2022"
output: 
  html_document: 
    toc: yes
    keep_md: yes

---



#### Load required packages


```r
library(tidyverse)
library(phyloseq)
library(microbiome)
library(here)
options(getClass.msg=FALSE) # https://github.com/epurdom/clusterExperiment/issues/66
#this fixes an error message that pops up because the class 'Annotated' is defined in two different packages
source("https://raw.githubusercontent.com/fconstancias/metabarcodingRpipeline/dev/scripts/functions_export_simplified.R")
'%!in%' <- function(x,y)!('%in%'(x,y))
```


# Import phyloseq object


```r
ps = "data/processed/16S/1/ps_silva_dada2_human_chicken_meta.RDS"


ps %>% 
  here::here() %>%
  readRDS()  -> physeq

# physeq$physeq -> physeq

physeq
```

```
## phyloseq-class experiment-level object
## otu_table()   OTU Table:         [ 1155 taxa and 600 samples ]
## sample_data() Sample Data:       [ 600 samples by 63 sample variables ]
## tax_table()   Taxonomy Table:    [ 1155 taxa by 7 taxonomic ranks ]
## phy_tree()    Phylogenetic Tree: [ 1155 tips and 1154 internal nodes ]
## refseq()      DNAStringSet:      [ 1155 reference sequences ]
```


```r
physeq %>% 
  sample_names() %>% 
  sort() %>% 
  head()
```

```
## [1] "AB24-1-S55"  "AB24-1E-S38" "AB24-2-S48"  "AB24-2E-S23" "AB24-3-S71" 
## [6] "AB24-3E-S25"
```



```r
# here::here("data/raw/metabarcoding/merged_chicken_human_04.08.2021.tsv") %>%
here::here("data/raw/25.11.2021_metadata_updated.tsv") %>%
  readr::read_tsv() %>% 
  pull("sample") %>% 
  sort() %>% 
  head()
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

```
## [1] "AB24-1-S55"  "AB24-1E-S38" "AB24-2-S48"  "AB24-2E-S23" "AB24-3-S71" 
## [6] "AB24-3E-S25"
```


```r
physeq %>%
  physeq_add_metadata(physeq = .,
                      metadata = here::here("data/raw/25.11.2021_metadata_updated.tsv") %>%
                        readr::read_tsv(), sample_column = "sample") -> physeq_meta
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

```r
physeq;physeq_meta
```

```
## phyloseq-class experiment-level object
## otu_table()   OTU Table:         [ 1155 taxa and 600 samples ]
## sample_data() Sample Data:       [ 600 samples by 63 sample variables ]
## tax_table()   Taxonomy Table:    [ 1155 taxa by 7 taxonomic ranks ]
## phy_tree()    Phylogenetic Tree: [ 1155 tips and 1154 internal nodes ]
## refseq()      DNAStringSet:      [ 1155 reference sequences ]
```

```
## phyloseq-class experiment-level object
## otu_table()   OTU Table:         [ 1155 taxa and 600 samples ]
## sample_data() Sample Data:       [ 600 samples by 61 sample variables ]
## tax_table()   Taxonomy Table:    [ 1155 taxa by 7 taxonomic ranks ]
## phy_tree()    Phylogenetic Tree: [ 1155 tips and 1154 internal nodes ]
## refseq()      DNAStringSet:      [ 1155 reference sequences ]
```




```r
physeq %>%
  physeq_add_metadata(physeq = .,
                      metadata = here::here("data/raw/25.11.2021_metadata_updated.tsv") %>%
                        readr::read_tsv() %>% 
                        mutate(Reactor_Treatment_Dose = ifelse(!is.na(`Antibiotic_mg/mL`),
                                                               paste0(Reactor_Treatment, `Antibiotic_mg/mL`), Reactor_Treatment),
                               Treatment_Dose = ifelse(!is.na(`Antibiotic_mg/mL`),
                                                       paste0(Treatment, `Antibiotic_mg/mL`), Reactor_Treatment))
                      , sample_column = "sample") -> physeq_meta
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

```r
# 
# meta %>% 
#     mutate(Reactor_Treatment_Dose = ifelse(!is.na(`Antibiotic_mg/mL`), 
#        paste0(Reactor_Treatment, `Antibiotic_mg/mL`), Reactor_Treatment)
#     Treatment_Dose = ifelse(!is.na(`Antibiotic_mg/mL`), 
#        paste0(Treatment, `Antibiotic_mg/mL`), Reactor_Treatment))

# ifelse(!is.na(sample_data(physeq_meta)$Antibiotic_mg.mL), 
#        sample_data(physeq_meta)$Antibiotic_mg.mL, "") %>% 
# mutate(gradebook, Pass.Fail = ifelse(grade > 60, "Pass", "Fail"))



physeq;physeq_meta
```

```
## phyloseq-class experiment-level object
## otu_table()   OTU Table:         [ 1155 taxa and 600 samples ]
## sample_data() Sample Data:       [ 600 samples by 63 sample variables ]
## tax_table()   Taxonomy Table:    [ 1155 taxa by 7 taxonomic ranks ]
## phy_tree()    Phylogenetic Tree: [ 1155 tips and 1154 internal nodes ]
## refseq()      DNAStringSet:      [ 1155 reference sequences ]
```

```
## phyloseq-class experiment-level object
## otu_table()   OTU Table:         [ 1155 taxa and 600 samples ]
## sample_data() Sample Data:       [ 600 samples by 63 sample variables ]
## tax_table()   Taxonomy Table:    [ 1155 taxa by 7 taxonomic ranks ]
## phy_tree()    Phylogenetic Tree: [ 1155 tips and 1154 internal nodes ]
## refseq()      DNAStringSet:      [ 1155 reference sequences ]
```



600 samples initially and 471 with metadata ????


```r
intersect(
  sample_names(physeq_meta),
  sample_names(physeq)) %>% 
  length()
```

```
## [1] 600
```


```r
difference <- function(x, y) {
  c(setdiff(x, y), setdiff(y, x))
}

difference(
  sample_names(physeq),
  sample_names(physeq_meta))
```

```
## character(0)
```

Above missing from the metadata but present in the initiall phyloseq object.

Update factor levels:


```r
physeq_meta %>% 
  sample_data() %>% 
  data.frame() %>% 
  select(-input:-length_filtered_pc,-I7_Index_ID:-index2) %>% 
  mutate(Treatment2 = factor(Treatment2, levels = c("DONOR", "UNTREATED", "AB+E. coli", "AB", "E. coli", "AB+E. faecium", "E. faecium")),
         Reactor_Treatment = base::replace(Reactor_Treatment, Reactor_Treatment == "TR_2CTX", "TR2_CTX"),
         Phase = factor(Phase, levels = c("Stab", "Treat")),
         Model = factor(Model, levels = c("Chicken", "Human")),
         Fermentation = factor(Fermentation, levels = c("1", "2")),
         Period = case_when( # if chicken we have period if human
           Day_of_Treatment == -2 | Day_of_Treatment == -3 | Day_of_Treatment == -6 | Day_of_Treatment == -7 ~ "pret",
           Day_of_Treatment < 10 & Day_of_Treatment >= 0 ~ "t1",
           Day_of_Treatment < 20 & Day_of_Treatment >= 10 ~ "t2",
           Day_of_Treatment < 30 & Day_of_Treatment >= 20 ~ "t3",
           Day_of_Treatment < 40 & Day_of_Treatment >= 30 ~ "t4",
           Day_of_Treatment < 50 & Day_of_Treatment >= 40 ~ "t5"),
         Period = factor(Period, levels = c("pret", "t1", "t2", "t3", "t4", "t5")),
         Reactor_Treatment_Dose = if_else(!is.na(`Antibiotic_mg.mL`), paste0(Reactor_Treatment, `Antibiotic_mg.mL`), Reactor_Treatment),
         Treatment_Dose = if_else(!is.na(`Antibiotic_mg.mL`), paste0(Treatment, `Antibiotic_mg.mL`), Treatment),
         Treatment = base::replace(Treatment, Experiment == "Cecum", "DONOR"),
         Treatment = factor(Treatment, levels = c("negative","positive","STRAIN", "DONOR", "UNTREATED", "CTX", "CTX+HV292.1", "HV292.1", "VAN", "VAN+CCUG59168", "CCUG59168")),
         # Reactor = factor(Reactor, levels = c("negative","positive","STRAIN", "DONOR", "IR1", "IR2", "CR", "TR1", "TR2", "TR3", "TR4", "TR5", "TR6","Batch3","Batch4")),
         Experiment = factor(Experiment, levels = c( "NTC","Mock", "HV292.1", "CCUG59168", "Cecum","Continuous","Batch")),
         
         # mutate( Treatment  = case_when(Experiment == "Cecum" ~ "DONOR",
         #                               TRUE ~ Treatment)) %>% 
         # mutate(Treatment = ifelse(Experiment == "Cecum", "Donor", Treatment)) %>% 
         
  ) -> sample_data(physeq_meta)
# mutate(Reactor_Treatment = base::replace(Reactor_Treatment, Reactor_Treatment == "TR_2CTX", "TR2_CTX")) %>%
#   mutate(Reactor_Treatment_Dose = if_else(!is.na(`Antibiotic_mg.mL`), paste0(Reactor_Treatment, `Antibiotic_mg.mL`), Reactor_Treatment)) %>% 
#   mutate(Treatment_Dose = if_else(!is.na(`Antibiotic_mg.mL`), paste0(Treatment, `Antibiotic_mg.mL`), Treatment)) %>%
#   mutate(
#     Period = case_when(
#       Day_of_Treatment <= 0  ~ "pret",
#       Day_of_Treatment > 0 ~ as.character(Day_of_Treatment)
#     )) %>%
#   mutate(Period = base::replace(Period, Treatment == "DONOR", "pret")) %>%
#   mutate(Day_of_Treatment = addNA(Day_of_Treatment)) %>%
```


```r
# physeq_meta %>% 
#   add_phylogeny_to_phyloseq(export = FALSE) %>% 
#   saveRDS(here::here("data/raw/metabarcoding/ps_silva_dada2_humanonly_phylo_meta.RDS"))
```


```r
physeq_meta %>%
  saveRDS(here::here("data/processed/16S/1/ps_silva_dada2_human_chicken_meta_fact.RDS"))
```


```r
sessionInfo()
```

```
## R version 4.2.0 (2022-04-22)
## Platform: x86_64-apple-darwin17.0 (64-bit)
## Running under: macOS Mojave 10.14.6
## 
## Matrix products: default
## BLAS:   /Library/Frameworks/R.framework/Versions/4.2/Resources/lib/libRblas.0.dylib
## LAPACK: /Library/Frameworks/R.framework/Versions/4.2/Resources/lib/libRlapack.dylib
## 
## locale:
## [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
## 
## attached base packages:
## [1] stats     graphics  grDevices utils     datasets  methods   base     
## 
## other attached packages:
##  [1] here_1.0.1           microbiome_1.18.0    phyloseq_1.40.0     
##  [4] forcats_0.5.1        stringr_1.4.0        dplyr_1.0.9         
##  [7] purrr_0.3.4          readr_2.1.2          tidyr_1.2.0         
## [10] tibble_3.1.7         ggplot2_3.3.6        tidyverse_1.3.1.9000
## 
## loaded via a namespace (and not attached):
##   [1] nlme_3.1-157           bitops_1.0-7           fs_1.5.2              
##   [4] bit64_4.0.5            lubridate_1.8.0        httr_1.4.3            
##   [7] rprojroot_2.0.3        GenomeInfoDb_1.32.2    tools_4.2.0           
##  [10] backports_1.4.1        bslib_0.3.1            vegan_2.6-2           
##  [13] utf8_1.2.2             R6_2.5.1               mgcv_1.8-40           
##  [16] DBI_1.1.2              BiocGenerics_0.42.0    colorspace_2.0-3      
##  [19] permute_0.9-7          rhdf5filters_1.8.0     ade4_1.7-19           
##  [22] withr_2.5.0            tidyselect_1.1.2       bit_4.0.4             
##  [25] compiler_4.2.0         cli_3.3.0              rvest_1.0.2           
##  [28] Biobase_2.56.0         xml2_1.3.3             sass_0.4.1            
##  [31] scales_1.2.0           digest_0.6.29          rmarkdown_2.14        
##  [34] XVector_0.36.0         pkgconfig_2.0.3        htmltools_0.5.2       
##  [37] dbplyr_2.2.0           fastmap_1.1.0          rlang_1.0.3           
##  [40] readxl_1.4.0           rstudioapi_0.13        jquerylib_0.1.4       
##  [43] generics_0.1.2         jsonlite_1.8.0         vroom_1.5.7           
##  [46] googlesheets4_1.0.0    RCurl_1.98-1.7         magrittr_2.0.3        
##  [49] GenomeInfoDbData_1.2.8 biomformat_1.24.0      Matrix_1.4-1          
##  [52] Rhdf5lib_1.18.2        Rcpp_1.0.8.3           munsell_0.5.0         
##  [55] S4Vectors_0.34.0       fansi_1.0.3            ape_5.6-2             
##  [58] lifecycle_1.0.1        stringi_1.7.6          yaml_2.3.5            
##  [61] MASS_7.3-57            zlibbioc_1.42.0        Rtsne_0.16            
##  [64] rhdf5_2.40.0           plyr_1.8.7             grid_4.2.0            
##  [67] parallel_4.2.0         crayon_1.5.1           lattice_0.20-45       
##  [70] splines_4.2.0          Biostrings_2.64.0      haven_2.5.0           
##  [73] multtest_2.52.0        hms_1.1.1              knitr_1.39            
##  [76] pillar_1.7.0           igraph_1.3.1           dtplyr_1.2.1          
##  [79] reshape2_1.4.4         codetools_0.2-18       stats4_4.2.0          
##  [82] reprex_2.0.1           glue_1.6.2             evaluate_0.15         
##  [85] data.table_1.14.2      modelr_0.1.8           vctrs_0.4.1           
##  [88] tzdb_0.3.0             foreach_1.5.2          cellranger_1.1.0      
##  [91] gtable_0.3.0           assertthat_0.2.1       xfun_0.31             
##  [94] broom_1.0.0            survival_3.3-1         googledrive_2.0.0     
##  [97] gargle_1.2.0           iterators_1.0.14       IRanges_2.30.0        
## [100] cluster_2.1.3          ellipsis_0.3.2
```

