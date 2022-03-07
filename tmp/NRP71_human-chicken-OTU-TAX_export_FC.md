---
title: "NRP72 -metadata - humann - chiken "
author: "Florentin Constancias"
date: "March 07, 2022"
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
```


# Import phyloseq object


```r
ps = "data/processed/16S/ps_silva_dada2_human_chicken_meta.RDS"

ps %>% 
  here::here() %>%
  readRDS()  -> physeq

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
source("https://raw.githubusercontent.com/fconstancias/metabaRpipe/master/Rscripts/functions.R")
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
## The following object is masked from 'package:microbiome':
## 
##     alpha
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
physeq %>% 
  phyloseq_get_strains_fast() %>% 
  phloseq_export_otu_tax() %>% 
  write_tsv( here::here("data/processed/16S/ps_silva_dada2_human_chicken_meta.tsv")) 
```

```
## Joining, by = "ASV"
```


```r
sessionInfo()
```

```
## R version 4.0.2 (2020-06-22)
## Platform: x86_64-apple-darwin17.0 (64-bit)
## Running under: macOS High Sierra 10.13.6
## 
## Matrix products: default
## BLAS:   /Library/Frameworks/R.framework/Versions/4.0/Resources/lib/libRblas.dylib
## LAPACK: /Library/Frameworks/R.framework/Versions/4.0/Resources/lib/libRlapack.dylib
## 
## locale:
## [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
## 
## attached base packages:
## [1] stats     graphics  grDevices utils     datasets  methods   base     
## 
## other attached packages:
##  [1] reshape2_1.4.4    scales_1.1.1      here_1.0.1        microbiome_1.10.0
##  [5] phyloseq_1.34.0   forcats_0.5.1     stringr_1.4.0     dplyr_1.0.7      
##  [9] purrr_0.3.4       readr_2.0.1       tidyr_1.1.4       tibble_3.1.6     
## [13] ggplot2_3.3.5     tidyverse_1.3.1  
## 
## loaded via a namespace (and not attached):
##  [1] nlme_3.1-152        fs_1.5.0            bit64_4.0.5        
##  [4] lubridate_1.8.0     httr_1.4.2          rprojroot_2.0.2    
##  [7] tools_4.0.2         backports_1.4.1     bslib_0.2.5.1      
## [10] vegan_2.5-7         utf8_1.2.2          R6_2.5.1           
## [13] mgcv_1.8-36         DBI_1.1.1           BiocGenerics_0.34.0
## [16] colorspace_2.0-2    permute_0.9-5       ade4_1.7-17        
## [19] withr_2.4.3         tidyselect_1.1.1    bit_4.0.4          
## [22] compiler_4.0.2      cli_3.1.0           rvest_1.0.1        
## [25] Biobase_2.50.0      xml2_1.3.2          sass_0.4.0         
## [28] digest_0.6.29       rmarkdown_2.10      XVector_0.28.0     
## [31] pkgconfig_2.0.3     htmltools_0.5.2     dbplyr_2.1.1       
## [34] fastmap_1.1.0       rlang_0.4.12        readxl_1.3.1       
## [37] rstudioapi_0.13     jquerylib_0.1.4     generics_0.1.1     
## [40] jsonlite_1.7.2      vroom_1.5.4         magrittr_2.0.1     
## [43] biomformat_1.16.0   Matrix_1.3-4        Rcpp_1.0.7         
## [46] munsell_0.5.0       S4Vectors_0.26.1    Rhdf5lib_1.10.1    
## [49] fansi_0.5.0         ape_5.5             lifecycle_1.0.1    
## [52] stringi_1.7.6       yaml_2.2.1          MASS_7.3-54        
## [55] zlibbioc_1.34.0     Rtsne_0.15          rhdf5_2.32.4       
## [58] plyr_1.8.6          grid_4.0.2          parallel_4.0.2     
## [61] crayon_1.4.2        lattice_0.20-44     splines_4.0.2      
## [64] Biostrings_2.56.0   haven_2.4.3         multtest_2.44.0    
## [67] hms_1.1.0           knitr_1.37          pillar_1.6.4       
## [70] igraph_1.2.6        codetools_0.2-18    stats4_4.0.2       
## [73] reprex_2.0.1        glue_1.6.0          evaluate_0.14      
## [76] data.table_1.14.2   modelr_0.1.8        vctrs_0.3.8        
## [79] tzdb_0.1.2          foreach_1.5.1       cellranger_1.1.0   
## [82] gtable_0.3.0        assertthat_0.2.1    xfun_0.29          
## [85] broom_0.7.11        survival_3.2-13     iterators_1.0.13   
## [88] IRanges_2.22.2      cluster_2.1.2       ellipsis_0.3.2
```

