---
title: "Chicken Resistome Visualization"
author: "Hannah Li Hägi & FC"
date: "January 06, 2022"
output: 
  html_document: 
    toc: yes
    keep_md: yes
---



#### Load required packages


```r
library(tidyverse)
library(phyloseq)
library(speedyseq)
library(ggrepel)
library(here)
library(gridExtra)
library(readxl)
options(getClass.msg=FALSE) # https://github.com/epurdom/clusterExperiment/issues/66
#this fixes an error message that pops up because the class 'Annotated' is defined in two different packages
```

#### Source required functions


```r
rm(list = ls())

'%!in%' <- function(x,y)!('%in%'(x,y))

source("https://raw.githubusercontent.com/fconstancias/DivComAnalyses/master/R/phyloseq_biobakery_functions.R")
```

# Import data:


```r
chick = "~/Projects/ETH/Alessia/mobilome/chicken/chicken_resistome_phyloseq.rds"
human = "~/Projects/ETH/Alessia/mobilome/human/human_resistome_phyloseq.rds"

chick %>% 
  readRDS() -> ps_chick


human %>% 
  readRDS() -> ps_human
```




```r
taxa_names(ps_chick) <- paste0("chick_", taxa_names(ps_chick))
taxa_names(ps_human) <- paste0("human_", taxa_names(ps_human))
```





```r
as(tax_table(ps_chick), "matrix") %>% 
  as.data.frame() %>%
  rownames_to_column('ASV') -> tax_table_chick

as(otu_table(ps_chick), "matrix") %>% 
  as.data.frame() %>%
  rownames_to_column('ASV') -> otu_table_chick

ps_chick %>% 
  sample_data() %>% 
  data.frame() %>% 
  rownames_to_column('sample_id') ->  meta_chick


as(tax_table(ps_human), "matrix") %>% 
  as.data.frame() %>%
  rownames_to_column('ASV') -> tax_table_human

as(otu_table(ps_human), "matrix") %>% 
  as.data.frame() %>%
  rownames_to_column('ASV') -> otu_table_human

ps_human %>% 
  sample_data() %>% 
  data.frame() %>% 
  rownames_to_column('sample_id') ->  meta_human
```

We will combine everything at this stage so we can later agglomerate at the Best_Hit_ARO:

```r
otu_table_human %>%
  # left_join(refseq_df, by = 'ASV') %>%
  left_join(tax_table_human, by = c("ASV" = "ASV")) %>%
  dplyr::select(ASV, everything()) -> merged_human

otu_table_chick %>%
  # left_join(refseq_df, by = 'ASV') %>%
  left_join(tax_table_chick, by = c("ASV" = "ASV")) %>%
  dplyr::select(ASV, everything()) -> merged_chick
```


```r
full_join(merged_human,
          merged_chick) -> full_merged
```

```
## Joining, by = c("ASV", "Best_Hit_ARO", "Model_type", "Drug_Class", "Resistance_Mechanism", "AMR_Gene Family", "Best_Identities", "Percentage_Length of Reference Sequence", "Note", "Predicted_DNA", "Length_AA", "Gene_name", "Contig_ID", "Tax", "KEGG_ID", "KEGGFUN", "KEGGPATH", "COG_ID", "COGFUN", "COGPATH", "PFAM", "Drug_Class_multi")
```

```r
full_merged %>% 
  write_tsv("~/Projects/ETH/Alessia/mobilome/meged.tsv")
```

Seems perfect. LEt's do the same for the metadata.



```r
full_join(meta_human,
          meta_chick) -> meta_merged
```

```
## Joining, by = c("sample_id", "sample", "input", "filtered", "denoisedF", "denoisedR", "merged", "tabled", "filtered_pc", "denoisedF_pc", "denoisedR_pc", "merged_pc", "filtered_merged_pc", "input_merged_pc", "tabled_joined", "chimera_out", "length_filtered", "tabled_pc", "chimera_out_pc", "length_filtered_pc", "Sample_description", "I7_Index_ID", "index", "I5_Index_ID", "index2", "Description2", "Experiment", "Reactor", "Treatment", "Day_of_Connection", "Day_of_Treatment", "Day_from_Inoculum", "Enrichment", "Phase", "Treatment2", "Date", "Paul", "Reactor_Treatment", "GeneCopyNumberperML", "HV292.1_Copy_Number_permL", "CCUG59168_Copy_Number_permL", "CTX_Copy_Number_permL", "VAN_Copy_Number_permL", "Model", "Antibiotic_mg.mL", "Fermentation", "Antibiotic", "Lactose_mM", "Glucose_mM", "Galactose_mM", "Succinat_mM", "Lactat_mM", "Formiat_mM", "Acetat_mM", "Propionat_mM", "Isobutyrat_mM", "Butyrat_mM", "Isovalerat_mM", "Valerat_mM", "Total_SCFA_mM", "raw_metagenomic_pairs")
```

```r
meta_merged %>% 
  write_tsv("~/Projects/ETH/Alessia/mobilome/meta_meged.tsv")
```

Create merged phyloseq object:


```r
full_merged %>% 
  select(ASV, 
         sample_names(ps_chick), sample_names(ps_human)) %>% 
  column_to_rownames("ASV") %>% 
  as.matrix() %>% 
  otu_table(taxa_are_rows = TRUE) -> otu_combined

full_merged %>% 
  select(ASV, 
         colnames(tax_table(ps_chick))) %>% 
  column_to_rownames("ASV") %>% 
  as.matrix() %>% 
  tax_table() -> tax_combined


phyloseq(otu_combined,
         tax_combined,
         meta_merged %>%
           column_to_rownames("sample_id") %>% 
           sample_data()) -> ps_combined


nsamples(ps_chick) + nsamples(ps_human); nsamples(ps_combined); print("looks good !")
```

```
## [1] 92
```

```
## [1] 92
```

```
## [1] "looks good !"
```


```r
ps_combined
```

```
## phyloseq-class experiment-level object
## otu_table()   OTU Table:          [ 244 taxa and 92 samples ]:
## sample_data() Sample Data:        [ 92 samples by 60 sample variables ]:
## tax_table()   Taxonomy Table:     [ 244 taxa by 21 taxonomic ranks ]:
## taxa are rows
```




```r
ps_combined %>% 
  sample_data() %>% 
  data.frame() %>% 
  mutate_at(vars(Day_of_Treatment), ~ replace(., is.na(.), -1)) -> sample_data(ps_combined)
```



```r
as(tax_table(ps_combined), "matrix") %>% 
  as.data.frame() %>%
  rownames_to_column('ASV') -> tax_table_ps_combined

as(otu_table(ps_combined), "matrix") %>% 
  as.data.frame() %>%
  rownames_to_column('ASV') -> otu_table_ps_combined

ps_combined %>% 
  sample_data() %>% 
  data.frame() %>% 
  rownames_to_column('sample_id') ->  meta_ps_combined

otu_table_ps_combined %>%
  # left_join(refseq_df, by = 'ASV') %>%
  left_join(tax_table_ps_combined, by = c("ASV" = "ASV")) %>%
  dplyr::select(ASV, everything()) -> merged_ps_combined


merged_ps_combined %>% 
  write_tsv("~/Projects/ETH/Alessia/mobilome/merged_ps_combined.tsv")

meta_ps_combined %>% 
  write_tsv("~/Projects/ETH/Alessia/mobilome/merged_meta_ps_combined.tsv")

ps_combined %>% 
  saveRDS("~/Projects/ETH/Alessia/mobilome/ps_combined.rds")
```

<!-- Create phyloseq object: -->
<!-- ```{r, fig.height=18, fig.width=16} -->
<!-- fil_df %>%  -->
<!--   # select(ORF_ID, Best_Hit_ARO, Model_type, `Drug Class`, `Resistance Mechanism`, `AMR Gene Family`, -->
<!--   #        starts_with("TPM "), Predicted_DNA) %>%  -->
<!--   mutate(`Drug Class` = str_remove_all(`Drug Class`," antibiotic")) %>%  -->
<!--   mutate(Drug_Class_multi = `Drug Class`) %>%  -->

<!--   column_to_rownames("ORF_ID") -> ps_ready -->

<!-- ps_ready$Drug_Class_multi[grepl(";", ps_ready$Drug_Class_multi)] <- "multi-drug" -->

<!-- ps_ready %>%  -->
<!--   select(-starts_with("TPM ")) %>%  -->
<!--   as.matrix() %>%  -->
<!--   tax_table() -> ps_tax -->

<!-- ps_ready %>%  -->
<!--   select(starts_with("TPM ")) %>%  -->
<!--   as.matrix() %>%  -->
<!--   otu_table(taxa_are_rows = TRUE) -> ps_otu -->


<!-- ps_ready$Predicted_DNA -> DNA -->

<!-- names(DNA) <- ps_ready$ORF_ID -->

<!-- physeq = phyloseq(ps_tax,  -->
<!--                   ps_otu) %>%  -->
<!--   clean_phyloseq_sample_names(sub_pat =  "TPM ") %>%  -->
<!--   clean_phyloseq_sample_names(sub_pat =  "-", str_replace = "_") -->

<!-- # physeq@refseq <- Biostrings::DNAStringSet(DNA) -->

<!-- merge_phyloseq(physeq, -->
<!--                meta %>% column_to_rownames("metagenomic_sample_name") %>% sample_data()) -> physeq -->

<!-- colnames(tax_table(physeq)) <- str_replace(colnames(tax_table(physeq)), " ", "_") -->





<!-- sample_data(physeq)$Reactor_Treatment <- fct_relevel(sample_data(physeq)$Reactor_Treatment, "DONOR", "CR_UNTREATED", "TR1_CTX+HV292.1", -->
<!--                                                      "TR2_CTX", "TR3_HV292.1", "TR5_VAN+CCUG59168", "TR4_VAN", "TR6_CCUG59168")  -->

<!-- sample_data(physeq)$Treatment <- fct_relevel(sample_data(physeq)$Treatment, "DONOR", "UNTREATED",  "CTX+HV292.1", "CTX","HV292.1","VAN+CCUG59168", "VAN",  "CCUG59168") -->


<!-- # source("https://raw.githubusercontent.com/fconstancias/metabarcodingRpipeline/dev/scripts/functions_export_simplified.R") -->
<!-- #  -->
<!-- # physeq %>%  -->
<!-- #   phloseq_export_otu_tax() %>%  -->
<!-- #   write_tsv("human_resistome_phyloseq_rds.tsv") -->

<!-- as(tax_table(physeq), "matrix") %>%  -->
<!--   as.data.frame() %>% -->
<!--   rownames_to_column('ASV') -> tax_table -->

<!-- as(otu_table(physeq), "matrix") %>%  -->
<!--   as.data.frame() %>% -->
<!--   rownames_to_column('ASV') -> otu_table -->

<!-- otu_table %>% -->
<!--   # left_join(refseq_df, by = 'ASV') %>% -->
<!--   left_join(tax_table, by = c("ASV" = "ASV")) %>% -->
<!--   dplyr::select(ASV, everything()) %>%  -->
<!--   write_tsv("~/Projects/ETH/Alessia/mobilome/chicken/chicken_resistome_phyloseq_rds.tsv") -->

<!-- saveRDS(physeq, file = "~/Projects/ETH/Alessia/mobilome/chicken/chicken_resistome_phyloseq.rds") -->

<!-- #  -->
<!-- #  -->
<!-- # colnames(ps_ready) <- sub("TPM ", "", colnames(ps_ready)) -->
<!-- #  -->
<!-- # agg_df <- aggregate(fil_df[c(paste0("TPM D", 1:43))], by = fil_df["Best_Hit_ARO"] , sum, drop = F)  -->
<!-- #  -->
<!-- # agg_df %>%  -->
<!-- #   left_join(., fil_df[c("Best_Hit_ARO", "Drug.Class")], by ="Best_Hit_ARO") %>% -->
<!-- #   unique %>% -->
<!-- #   as.data.frame() %>% -->
<!-- #   remove_rownames() %>% -->
<!-- #   column_to_rownames(., "Best_Hit_ARO") -> df  -->
<!-- #  -->
<!-- #  -->
<!-- # colnames(df) <- sub("TPM ", "", colnames(df)) -->
<!-- # df$Drug.Class <- str_remove_all(df$Drug.Class," antibiotic") -->
<!-- #  -->
<!-- # # Rename the entries containing more than 1 drug resistance -->
<!-- # df$Drug.Class[grepl(";", df$Drug.Class)] <- "multi-drug" -->
<!-- #  -->
<!-- # #create tax_table -->
<!-- # df %>% -->
<!-- #   dplyr::select(c("Drug.Class")) %>% -->
<!-- #   as.matrix() %>% -->
<!-- #   tax_table()-> TAX -->
<!-- #  -->
<!-- # #create otu table containing the abundances of each gen -->
<!-- # df %>%  -->
<!-- #   dplyr::select(c(paste0("D", 1:43))) %>%  -->
<!-- #   otu_table(taxa_are_rows = TRUE)-> OTU -->
<!-- #  -->
<!-- # #create sample_info table containing meta information -->
<!-- # meta$Treatment <- factor(meta$Treatment,  -->
<!-- #                          ordered = T,  -->
<!-- #                          levels = c("DONOR", "UNTREATED", "CCUG59168", "VAN",  -->
<!-- #                                     "VAN+CCUG59168", "HV292.1", "CTX", "CTX+HV292.1")) -->
<!-- #  -->
<!-- # meta$Day_of_Treatment <- factor(meta$Day_of_Treatment,  -->
<!-- #                                 levels = c("-1", "1", "3", "5", "30", "48")) -->
<!-- #  -->
<!-- # meta %>% -->
<!-- #   as.data.frame() %>%  -->
<!-- #   dplyr::select(. , raw_metagenomic_pairs, metagenomic_sample_name, -->
<!-- #                 Treatment, Day_of_Treatment, Fermentation, Reactor_Treatment) %>%  -->
<!-- #   column_to_rownames(., "metagenomic_sample_name") %>% -->
<!-- #   mutate(Day_of_Treatment = as.factor(Day_of_Treatment)) %>% -->
<!-- #   sample_data()-> sampleinfo -->
<!-- #  -->
<!-- #  -->
<!-- # physeq = phyloseq(OTU, TAX, sampleinfo) -->
<!-- # saveRDS(physeq, file = "resistome_phyloseq.rds") -->

<!-- ``` -->




```r
sessionInfo()
```

```
## R version 4.1.2 (2021-11-01)
## Platform: x86_64-apple-darwin17.0 (64-bit)
## Running under: macOS Mojave 10.14.6
## 
## Matrix products: default
## BLAS:   /Library/Frameworks/R.framework/Versions/4.1/Resources/lib/libRblas.0.dylib
## LAPACK: /Library/Frameworks/R.framework/Versions/4.1/Resources/lib/libRlapack.dylib
## 
## locale:
## [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
## 
## attached base packages:
## [1] stats     graphics  grDevices utils     datasets  methods   base     
## 
## other attached packages:
##  [1] readxl_1.3.1         gridExtra_2.3        here_1.0.1          
##  [4] ggrepel_0.9.1        speedyseq_0.5.3.9018 phyloseq_1.36.0     
##  [7] forcats_0.5.1        stringr_1.4.0        dplyr_1.0.7         
## [10] purrr_0.3.4          readr_2.1.0          tidyr_1.1.4         
## [13] tibble_3.1.6         ggplot2_3.3.5        tidyverse_1.3.1.9000
## 
## loaded via a namespace (and not attached):
##  [1] nlme_3.1-153           bitops_1.0-7           fs_1.5.2              
##  [4] bit64_4.0.5            lubridate_1.8.0        httr_1.4.2            
##  [7] rprojroot_2.0.2        GenomeInfoDb_1.28.4    tools_4.1.2           
## [10] backports_1.4.1        bslib_0.3.1            vegan_2.5-7           
## [13] utf8_1.2.2             R6_2.5.1               mgcv_1.8-38           
## [16] DBI_1.1.1              BiocGenerics_0.38.0    colorspace_2.0-2      
## [19] permute_0.9-5          rhdf5filters_1.4.0     ade4_1.7-18           
## [22] withr_2.4.3            tidyselect_1.1.1       bit_4.0.4             
## [25] compiler_4.1.2         cli_3.1.0              rvest_1.0.2           
## [28] Biobase_2.52.0         xml2_1.3.2             sass_0.4.0            
## [31] scales_1.1.1           digest_0.6.29          rmarkdown_2.11        
## [34] XVector_0.32.0         pkgconfig_2.0.3        htmltools_0.5.2       
## [37] dbplyr_2.1.1           fastmap_1.1.0          rlang_0.4.12          
## [40] rstudioapi_0.13        jquerylib_0.1.4        generics_0.1.1        
## [43] jsonlite_1.7.2         vroom_1.5.6            RCurl_1.98-1.5        
## [46] magrittr_2.0.1         GenomeInfoDbData_1.2.6 biomformat_1.20.0     
## [49] Matrix_1.3-4           Rcpp_1.0.7             munsell_0.5.0         
## [52] S4Vectors_0.30.2       Rhdf5lib_1.14.2        fansi_0.5.0           
## [55] ape_5.6                lifecycle_1.0.1        stringi_1.7.6         
## [58] yaml_2.2.1             MASS_7.3-54            zlibbioc_1.38.0       
## [61] rhdf5_2.36.0           plyr_1.8.6             grid_4.1.2            
## [64] parallel_4.1.2         crayon_1.4.2           lattice_0.20-45       
## [67] splines_4.1.2          Biostrings_2.60.2      haven_2.4.3           
## [70] multtest_2.48.0        hms_1.1.1              knitr_1.36            
## [73] pillar_1.6.4           igraph_1.2.10          reshape2_1.4.4        
## [76] codetools_0.2-18       stats4_4.1.2           reprex_2.0.1          
## [79] glue_1.6.0             evaluate_0.14          data.table_1.14.2     
## [82] modelr_0.1.8           vctrs_0.3.8            tzdb_0.2.0            
## [85] foreach_1.5.1          cellranger_1.1.0       gtable_0.3.0          
## [88] assertthat_0.2.1       xfun_0.28              broom_0.7.11          
## [91] survival_3.2-13        iterators_1.0.13       IRanges_2.26.0        
## [94] cluster_2.1.2          ellipsis_0.3.2
```
