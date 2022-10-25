---
title: "NTP72 - data processing"
author: "Florentin Constancias"
date: "October 18, 2022"
output: 
  html_document: 
    toc: yes
    keep_md: yes
---




```r
#install.packages("tidyverse")
require(tidyverse); packageVersion("tidyverse")
```

```
## Loading required package: tidyverse
```

```
## ── Attaching packages ─────────────────────────────────────── tidyverse 1.3.2 ──
## ✔ ggplot2 3.3.6      ✔ purrr   0.3.5 
## ✔ tibble  3.1.8      ✔ dplyr   1.0.10
## ✔ tidyr   1.2.1      ✔ stringr 1.4.1 
## ✔ readr   2.1.3      ✔ forcats 0.5.2 
## ── Conflicts ────────────────────────────────────────── tidyverse_conflicts() ──
## ✖ dplyr::filter() masks stats::filter()
## ✖ dplyr::lag()    masks stats::lag()
```

```r
require(phyloseq); packageVersion("phyloseq")
```

```
## Loading required package: phyloseq
```


```r
source("https://raw.githubusercontent.com/fconstancias/DivComAnalyses/master/R/phyloseq_taxa_tests.R") 
```

```
## 
## Attaching package: 'nlme'
```

```
## The following object is masked from 'package:dplyr':
## 
##     collapse
```

```
## Welcome to compositions, a package for compositional data analysis.
## Find an intro with "? compositions"
```

```
## 
## Attaching package: 'compositions'
```

```
## The following objects are masked from 'package:stats':
## 
##     anova, cor, cov, dist, var
```

```
## The following objects are masked from 'package:base':
## 
##     %*%, norm, scale, scale.default
```

```r
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


```r
'%!in%' <- function(x,y)!('%in%'(x,y))
out_pptx = "~/Desktop/16S-human.pptx"
```


```r
load(url("https://github.com/fconstancias/NRP72-FBT/blob/master/Figures/Rastetics.Rdata?raw=true"))
```


## load data:



```r
(url("https://github.com/fconstancias/NRP72-FBT/blob/master/data/processed/16S/phyloseq_phylo_aline_human2_donor.RDS?raw=true" )) %>% 
  readRDS() %>% 
  phyloseq_get_strains() %>% 
  subset_samples(Sample_type == "Fecal Sample") -> ps_aline
```

```
## Joining, by = "ASV"
```

```r
sample_names(ps_aline) <- "S-238-S238"
```



```r
source("https://raw.githubusercontent.com/fconstancias/metabaRpipe-source/master/Rscripts/functions.R")
source("https://raw.githubusercontent.com/fconstancias/DivComAnalyses/master/R/phyloseq_normalisation.R") 

(url("https://github.com/fconstancias/NRP72-FBT/blob/master/data/processed/16S/ps_silva_dada2_human_chicken_all_fct_polyferms_051022.RDS?raw=true" )) %>% 
  readRDS() %>% 
  phyloseq_get_strains() %>% 
  # physeq_add_metadata(physeq = .,
  # metadata = ("~/Downloads/16S-NRP72/27.07.2022_metadata_updated_chicken.xlsx") %>%
  # readxl::read_excel(sheet = 1), sample_column = "sample_name") %>% 
  subset_samples(Sample_description %!in% c("TR5-15", "TR4-1")) -> ps_init
```

```
## Joining, by = "ASV"
```


```r
ps_init %>%  
  sample_data() %>% 
  data.frame() %>% 
  rownames_to_column("sample_names") -> metadata
```


```r
phyloseq_combine_objects(ps_init %>% 
                           subset_samples(sample_name != "S-238-S238") ,  
                         ps_aline %>%   filter_taxa(function(x) sum(x > 0) > 0, TRUE), 
                         merge_metada = FALSE) -> ps_all
```

```
## Loading required package: speedyseq
```

```
## 
## Attaching package: 'speedyseq'
```

```
## The following objects are masked from 'package:phyloseq':
## 
##     filter_taxa, plot_bar, plot_heatmap, plot_tree, psmelt, tax_glom,
##     tip_glom, transform_sample_counts
```

```
## Loading required package: DECIPHER
```

```
## Loading required package: Biostrings
```

```
## Loading required package: BiocGenerics
```

```
## 
## Attaching package: 'BiocGenerics'
```

```
## The following objects are masked from 'package:compositions':
## 
##     normalize, var
```

```
## The following objects are masked from 'package:dplyr':
## 
##     combine, intersect, setdiff, union
```

```
## The following objects are masked from 'package:stats':
## 
##     IQR, mad, sd, var, xtabs
```

```
## The following objects are masked from 'package:base':
## 
##     anyDuplicated, append, as.data.frame, basename, cbind, colnames,
##     dirname, do.call, duplicated, eval, evalq, Filter, Find, get, grep,
##     grepl, intersect, is.unsorted, lapply, Map, mapply, match, mget,
##     order, paste, pmax, pmax.int, pmin, pmin.int, Position, rank,
##     rbind, Reduce, rownames, sapply, setdiff, sort, table, tapply,
##     union, unique, unsplit, which.max, which.min
```

```
## Loading required package: S4Vectors
```

```
## Loading required package: stats4
```

```
## 
## Attaching package: 'S4Vectors'
```

```
## The following objects are masked from 'package:dplyr':
## 
##     first, rename
```

```
## The following object is masked from 'package:tidyr':
## 
##     expand
```

```
## The following objects are masked from 'package:base':
## 
##     expand.grid, I, unname
```

```
## Loading required package: IRanges
```

```
## 
## Attaching package: 'IRanges'
```

```
## The following object is masked from 'package:nlme':
## 
##     collapse
```

```
## The following object is masked from 'package:phyloseq':
## 
##     distance
```

```
## The following objects are masked from 'package:dplyr':
## 
##     collapse, desc, slice
```

```
## The following object is masked from 'package:purrr':
## 
##     reduce
```

```
## Loading required package: XVector
```

```
## 
## Attaching package: 'XVector'
```

```
## The following object is masked from 'package:purrr':
## 
##     compact
```

```
## Loading required package: GenomeInfoDb
```

```
## 
## Attaching package: 'Biostrings'
```

```
## The following object is masked from 'package:base':
## 
##     strsplit
```

```
## Loading required package: RSQLite
```

```
## Loading required package: parallel
```

```
## Joining, by = "ASV"
```

```
## [1] "Number of sequences clustered = 0"
```


```r
ps_all %>% 
  physeq_add_metadata(physeq = .,
                      metadata = metadata, 
                      sample_column = "sample_name") %>% 
  phyloseq_get_strains() -> ps_temp
```

```
## Joining, by = "ASV"
```

```r
# sample_data(ps_temp) %>% 
#   data.frame()
```


```r
# tax_table(ps_init)[,"Strain"]
```



```r
# tax_table(ps_temp)[,"Strain"]
```


```r
ps_temp %>% 
  # microViz::ps_mutate(Antibiotic_mg.L = ifelse(Phase == "Stab", NA, Antibiotic_mg.L))  %>%  sample_data() %>%  data.frame() %>%  select(Antibiotic_mg.L)
  microViz::ps_mutate(Day_of_Treatment = as.double(Day_of_Treatment)) %>% 
  # microViz::ps_mutate(Antibiotic_mg.L = factor(Antibiotic_mg.L, levels = c(NA, 20, 200,90, 600))) %>% 
  microViz::ps_mutate(Treatment = factor(Treatment, levels = c("DONOR","UNTREATED","CTX","CTX+HV292.1","HV292.1","VAN","VAN+CCUG59168", "CCUG59168"
  ))) %>% 
  microViz::ps_mutate(Reactor_Treatment_Dose = factor(Reactor_Treatment_Dose, levels = c("DONOR",  
                                                                                         "IR_UNTREATED",
                                                                                         "CR_UNTREATED",
                                                                                         "TR2_CTX20",
                                                                                         "TR3_CTX20",
                                                                                         "TR3_CTX200",
                                                                                         "TR4_CTX200",
                                                                                         "TR1_CTX+HV292.120",
                                                                                         "TR2_CTX+HV292.120",
                                                                                         "TR3_CTX+HV292.1200",
                                                                                         "TR5_CTX+HV292.1200",
                                                                                         "TR3_HV292.1",
                                                                                         "TR5_HV292.1",
                                                                                         "TR7_HV292.1",
                                                                                         "TR4_VAN90",
                                                                                         "TR2_VAN90",
                                                                                         "TR6_VAN90",
                                                                                         "TR4_VAN600",
                                                                                         "TR2_VAN600",
                                                                                         "TR5_VAN+CCUG5916890",
                                                                                         "TR1_VAN+CCUG5916890",
                                                                                         "TR7_VAN+CCUG5916890",
                                                                                         "TR3_VAN+CCUG59168600",
                                                                                         "TR4_VAN+CCUG59168600",
                                                                                         "TR6_CCUG59168",
                                                                                         "TR5_CCUG59168",
                                                                                         "TR4_CCUG59168"))) %>% 
  microViz::ps_mutate(Reactor_Treatment_Dose_fermentation = paste0(Reactor_Treatment_Dose, "_", Fermentation)) -> ps_temp

ps_temp %>%
  subset_samples(Enrichment == "NotEnriched") %>%
  subset_samples(Experiment %in% c("Cecum", "Continuous")) %>%
  subset_samples(Day_of_Treatment > -6 & Day_of_Treatment <=30 | Reactor == "DONOR") %>% 
  # subset_samples(Model == "Human") %>% 
  # subset_samples(Reactor != "CR" & Model2 == "Chicken1" | Model2 == "Chicken2" &   Reactor %in% c("TR2", "TR3", "TR4", "TR5", "TR6", "TR7", "CR", "DONOR")) %>% 
  subset_samples(Reactor != "IR") %>% 
   filter_taxa(function(x) sum(x > 0) > 0, TRUE)  -> ps_filtered
```


```r
ps_filtered %>% 
  saveRDS(here::here("data/processed/16S/16S_working_phyloseq.RDS"))
```


```r
sessionInfo()
```

```
## R version 4.2.1 (2022-06-23)
## Platform: x86_64-apple-darwin17.0 (64-bit)
## Running under: macOS Big Sur ... 10.16
## 
## Matrix products: default
## BLAS:   /Library/Frameworks/R.framework/Versions/4.2/Resources/lib/libRblas.0.dylib
## LAPACK: /Library/Frameworks/R.framework/Versions/4.2/Resources/lib/libRlapack.dylib
## 
## locale:
## [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
## 
## attached base packages:
## [1] parallel  stats4    stats     graphics  grDevices utils     datasets 
## [8] methods   base     
## 
## other attached packages:
##  [1] DECIPHER_2.24.0      RSQLite_2.2.18       Biostrings_2.64.1   
##  [4] GenomeInfoDb_1.32.4  XVector_0.36.0       IRanges_2.30.1      
##  [7] S4Vectors_0.34.0     BiocGenerics_0.42.0  speedyseq_0.5.3.9018
## [10] reshape2_1.4.4       scales_1.2.1         compositions_2.0-4  
## [13] nlme_3.1-159         phyloseq_1.40.0      forcats_0.5.2       
## [16] stringr_1.4.1        dplyr_1.0.10         purrr_0.3.5         
## [19] readr_2.1.3          tidyr_1.2.1          tibble_3.1.8        
## [22] ggplot2_3.3.6        tidyverse_1.3.2     
## 
## loaded via a namespace (and not attached):
##  [1] googledrive_2.0.0      colorspace_2.0-3       ellipsis_0.3.2        
##  [4] rprojroot_2.0.3        fs_1.5.2               rstudioapi_0.14       
##  [7] bit64_4.0.5            fansi_1.0.3            lubridate_1.8.0       
## [10] xml2_1.3.3             codetools_0.2-18       splines_4.2.1         
## [13] cachem_1.0.6           robustbase_0.95-0      knitr_1.40            
## [16] ade4_1.7-19            jsonlite_1.8.2         broom_1.0.1           
## [19] cluster_2.1.4          dbplyr_2.2.1           compiler_4.2.1        
## [22] httr_1.4.4             backports_1.4.1        assertthat_0.2.1      
## [25] Matrix_1.5-1           fastmap_1.1.0          gargle_1.2.1          
## [28] cli_3.4.1              htmltools_0.5.3        tools_4.2.1           
## [31] igraph_1.3.5           gtable_0.3.1           glue_1.6.2            
## [34] GenomeInfoDbData_1.2.8 Rcpp_1.0.9             Biobase_2.56.0        
## [37] cellranger_1.1.0       jquerylib_0.1.4        vctrs_0.4.2           
## [40] rhdf5filters_1.8.0     multtest_2.52.0        ape_5.6-2             
## [43] iterators_1.0.14       tensorA_0.36.2         xfun_0.33             
## [46] rvest_1.0.3            lifecycle_1.0.3        googlesheets4_1.0.1   
## [49] DEoptimR_1.0-11        zlibbioc_1.42.0        MASS_7.3-58.1         
## [52] hms_1.1.2              biomformat_1.24.0      rhdf5_2.40.0          
## [55] yaml_2.3.5             memoise_2.0.1          sass_0.4.2            
## [58] stringi_1.7.8          foreach_1.5.2          permute_0.9-7         
## [61] rlang_1.0.6            pkgconfig_2.0.3        bitops_1.0-7          
## [64] evaluate_0.17          lattice_0.20-45        Rhdf5lib_1.18.2       
## [67] bit_4.0.4              tidyselect_1.2.0       here_1.0.1            
## [70] plyr_1.8.7             magrittr_2.0.3         R6_2.5.1              
## [73] generics_0.1.3         DBI_1.1.3              pillar_1.8.1          
## [76] haven_2.5.1            withr_2.5.0            mgcv_1.8-40           
## [79] survival_3.4-0         RCurl_1.98-1.9         bayesm_3.1-4          
## [82] modelr_0.1.9           crayon_1.5.2           utf8_1.2.2            
## [85] tzdb_0.3.0             rmarkdown_2.16         grid_4.2.1            
## [88] readxl_1.4.1           data.table_1.14.2      blob_1.2.3            
## [91] vegan_2.6-2            reprex_2.0.2           digest_0.6.29         
## [94] munsell_0.5.0          microViz_0.9.2         bslib_0.4.0
```

