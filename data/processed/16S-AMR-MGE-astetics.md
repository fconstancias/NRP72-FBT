---
title: "Plotting preparation"
author: "Florentin Constancias "
date: " July 15, 2022 "
output: 
  html_document: 
    toc: yes
    keep_md: yes
---




```r
rm(list= ls())
# gc()
require("tidyverse")
```

```
## Loading required package: tidyverse
```

```
## ── Attaching packages ────────────────────────────────── tidyverse 1.3.1.9000 ──
```

```
## ✔ ggplot2 3.3.6     ✔ purrr   0.3.4
## ✔ tibble  3.1.7     ✔ dplyr   1.0.9
## ✔ tidyr   1.2.0     ✔ stringr 1.4.0
## ✔ readr   2.1.2     ✔ forcats 0.5.1
```

```
## ── Conflicts ────────────────────────────────────────── tidyverse_conflicts() ──
## ✖ dplyr::filter() masks stats::filter()
## ✖ dplyr::lag()    masks stats::lag()
```

```r
require('speedyseq')
```

```
## Loading required package: speedyseq
```

```
## Loading required package: phyloseq
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

```r
'%!in%' <- function(x,y)!('%in%'(x,y))
source("https://raw.githubusercontent.com/fconstancias/DivComAnalyses/master/R/phyloseq_varia.R")  # source the functions from github.
```



```r
export_ppt_width = NULL
export_ppt_height = NULL
out_pptx = here::here("output/plots_2.pptx")
```

# Load data:

## 16S:


```r
"data/processed/16S/ps_silva_dada2_human_chicken_all_fct_polyferms.RDS" %>% 
  here::here() %>% 
  readRDS() %>% 
  filter_taxa(function(x) sum(x > 0) > 0, TRUE) -> ps_16S
```

## Export metadata:

To be merged with phyloseq data. TODO: Would be better to do that before in the gene phyloseq data generation.


```r
ps_16S %>% 
  sample_data() %>% 
  data.frame() %>% 
  filter(!is.na(metagenomic_sample_name)) %>% 
  rownames_to_column('tmp') %>% 
  select(-tmp) -> metadata
```

## phyloseq object of SQM + Pathofact + bins + ... with Num_Gi_pc quantif:


```r
# "data/processed/exp1_full_gene_catalog_phyloseq.RDS" %>% 
#   here::here() %>% 
#   readRDS() %>% 
#   .[["Num_Gi_pc"]] -> ps_coass
# 
# ps_coass
```

## extract on only AMR gene catalogue - qauntification on all gene catalogue:


```r
# ps_coass %>% 
#   subset_taxa(!is.na(PathoFact_AMR_ARG)) %>% 
#   subset_taxa(PathoFact_AMR_AMR_category != "unclassified") -> ps_AMR
```


```r
"data/processed/exp1_ps_AMR.RDS" %>% 
  here::here() %>% 
  readRDS() -> ps_AMR
```


```r
# ps_AMR %>% 
#   tax_table() %>% 
#   data.frame() %>% 
#   glimpse()
```

I hvae to do that otherwise it is not stored as proper chr when filtering... no clue what is happening here.



```r
ps_AMR %>%
  tax_table() %>%
  data.frame() %>%
  mutate(PathoFact_AMR_Resistance_mechanism_multi = ifelse(grepl(";", PathoFact_AMR_Resistance_mechanism), "multi-mech", PathoFact_AMR_Resistance_mechanism)) %>%
  mutate(PathoFact_AMR_AMR_sub_class_multi = ifelse(grepl(";", PathoFact_AMR_AMR_sub_class), "multidrug", PathoFact_AMR_AMR_sub_class)) %>%
  mutate(ANVIO_CONCOCT_HQ_clean = ifelse(PathoFact_AMR_MGE_prediction == "plasmid", NA, ANVIO_CONCOCT_HQ)) %>%
  mutate(GTDB_tax_CONCOCT_clean = ifelse(PathoFact_AMR_MGE_prediction == "plasmid", NA, GTDB_tax_CONCOCT)) %>%
  mutate(PathoFact_AMR_MGE_prediction_clean = ifelse(!is.na(ANVIO_CONCOCT_HQ), "MAG", PathoFact_AMR_MGE_prediction)) %>%
  mutate(PathoFact_AMR_MGE_prediction_clean_2 = case_when(grepl("ambiguous", PathoFact_AMR_MGE_prediction_clean) ~ "ambiguous", TRUE ~ PathoFact_AMR_MGE_prediction_clean)) %>%
  mutate(PathoFact_AMR_MGE_prediction_clean_2 = fct_relevel(PathoFact_AMR_MGE_prediction_clean_2, c("MAG", "chromosome", "plasmid", "phage", "ambiguous", "unclassified"))) %>%
  # mutate(Treatment = factor(Treatment, levels = c)) %>%
  # mutate(Day_of_Treatment = fct_relevel(Day_of_Treatment, c(NA, "-3", "-1", "1" , "3", "4", "5","7", "30", "47", "48")))  %>% 
  
  as.matrix() -> tax_table(ps_AMR)

ps_AMR %>% 
  saveRDS(., file = here::here("data/processed/exp1_ps_AMR.RDS") )
```




```r
# ps_AMR %>% 
#   tax_table() %>% 
#   colnames()
```

# Define astetics:

## Treatment colors:

<https://colorbrewer2.org/#type=diverging&scheme=RdYlBu&n=6>

```r
treat_col <- c("#3CB371", "#7C8181",'#d73027','#fc8d59','#fee090','#4575b4','#91bfdb', '#e0f3f8')

names(treat_col) <- get_variable(ps_AMR, "Treatment") %>%  
  levels()

treat_col %>% 
  scales::show_col(cex_label = 0.5)
```

<img src="16S-AMR-MGE-astetics_files/figure-html/unnamed-chunk-11-1.png" width="50%" style="display: block; margin: auto;" />

## Drug classes colors:


```r
ps_AMR %>% 
  generate_color_palette(var = "PathoFact_AMR_AMR_sub_class_multi", seed = 280386) -> drug_classe_multi
```

```
## Loading required package: randomcoloR
```

<img src="16S-AMR-MGE-astetics_files/figure-html/unnamed-chunk-12-1.png" width="50%" style="display: block; margin: auto;" />



```r
# ps_AMR %>% 
#   generate_color_palette(var = "Model_type", pal = "jco") -> model_type
```


```r
ps_AMR %>% 
  generate_color_palette(var = "PathoFact_AMR_Resistance_mechanism_multi",  pal = "npg") -> resistance_type
```

```
## Loading required package: ggpubr
```

<img src="16S-AMR-MGE-astetics_files/figure-html/unnamed-chunk-14-1.png" width="50%" style="display: block; margin: auto;" />


```r
ps_AMR %>% 
  generate_color_palette(var = "PathoFact_AMR_MGE_prediction_clean_2", seed = 72) -> MGE
```

```
## Loading required package: randomcoloR
```

<img src="16S-AMR-MGE-astetics_files/figure-html/unnamed-chunk-15-1.png" width="50%" style="display: block; margin: auto;" />

## Others:

<https://r-graphics.org/recipe-scatter-shapes>

```r
fermentaion_shape <- c(22, 21)
names(fermentaion_shape) <- sample_data(ps_AMR)$Fermentation %>%  levels()

model_shape <- c(22, 21)
names(model_shape) <- sample_data(ps_AMR)$Model %>%  levels()

antibio_shape <- c(22, 21)
names(antibio_shape) <- sample_data(ps_AMR)$Antibiotic %>%  levels()


conc_stroke <- c(1, 2, 3, 4, 5)
names(conc_stroke) <- sample_data(ps_AMR)$Antibiotic_mg.mL %>%  levels()
# scale_discrete_manual(
#   aesthetics = "stroke",
#   values = c(`A` = 2, `B` = 1, `C` = 2)
# )

theme_set(theme_classic() + theme(legend.position = 'bottom'))
```

save all ascetics parameters as Rdata which can be loaded afterwards:


```r
save(treat_col, 
     antibio_shape,
     drug_classe_multi,
     conc_stroke,
     MGE,
     fermentaion_shape,
     resistance_type,
     model_shape, file = here::here("Figures/Rastetics.Rdata") )
# 
# Then you can
# 
# load("location.filename.RData") 
```

<https://r4ds.had.co.nz/graphics-for-communication.html#figure-sizing>

<https://www.tidyverse.org/blog/2020/08/taking-control-of-plot-scaling/>




```r
devtools::session_info()
```

```
## ─ Session info ───────────────────────────────────────────────────────────────
##  setting  value
##  version  R version 4.2.1 (2022-06-23)
##  os       macOS Big Sur ... 10.16
##  system   x86_64, darwin17.0
##  ui       X11
##  language (EN)
##  collate  en_US.UTF-8
##  ctype    en_US.UTF-8
##  tz       Europe/Zurich
##  date     2022-07-15
##  pandoc   2.18 @ /Applications/RStudio.app/Contents/MacOS/quarto/bin/tools/ (via rmarkdown)
## 
## ─ Packages ───────────────────────────────────────────────────────────────────
##  package          * version    date (UTC) lib source
##  abind              1.4-5      2016-07-21 [1] CRAN (R 4.2.0)
##  ade4               1.7-19     2022-04-19 [1] CRAN (R 4.2.0)
##  ape                5.6-2      2022-03-02 [1] CRAN (R 4.2.0)
##  assertthat         0.2.1      2019-03-21 [1] CRAN (R 4.2.0)
##  backports          1.4.1      2021-12-13 [1] CRAN (R 4.2.0)
##  Biobase            2.56.0     2022-04-26 [1] Bioconductor
##  BiocGenerics       0.42.0     2022-04-26 [1] Bioconductor
##  biomformat         1.24.0     2022-04-26 [1] Bioconductor
##  Biostrings         2.64.0     2022-04-26 [1] Bioconductor
##  bitops             1.0-7      2021-04-24 [1] CRAN (R 4.2.0)
##  broom              1.0.0      2022-07-01 [1] CRAN (R 4.2.0)
##  bslib              0.3.1      2021-10-06 [1] CRAN (R 4.2.0)
##  cachem             1.0.6      2021-08-19 [1] CRAN (R 4.2.0)
##  callr              3.7.1      2022-07-13 [1] CRAN (R 4.2.1)
##  car                3.1-0      2022-06-15 [1] CRAN (R 4.2.0)
##  carData            3.0-5      2022-01-06 [1] CRAN (R 4.2.0)
##  cellranger         1.1.0      2016-07-27 [1] CRAN (R 4.2.0)
##  cli                3.3.0      2022-04-25 [1] CRAN (R 4.2.0)
##  cluster            2.1.3      2022-03-28 [1] CRAN (R 4.2.1)
##  codetools          0.2-18     2020-11-04 [1] CRAN (R 4.2.1)
##  colorspace         2.0-3      2022-02-21 [1] CRAN (R 4.2.0)
##  crayon             1.5.1      2022-03-26 [1] CRAN (R 4.2.0)
##  curl               4.3.2      2021-06-23 [1] CRAN (R 4.2.0)
##  data.table         1.14.2     2021-09-27 [1] CRAN (R 4.2.0)
##  DBI                1.1.3      2022-06-18 [1] CRAN (R 4.2.0)
##  dbplyr             2.2.1      2022-06-27 [1] CRAN (R 4.2.0)
##  devtools           2.4.3      2021-11-30 [1] CRAN (R 4.2.0)
##  digest             0.6.29     2021-12-01 [1] CRAN (R 4.2.0)
##  dplyr            * 1.0.9      2022-04-28 [1] CRAN (R 4.2.0)
##  dtplyr             1.2.1      2022-01-19 [1] CRAN (R 4.2.0)
##  ellipsis           0.3.2      2021-04-29 [1] CRAN (R 4.2.0)
##  evaluate           0.15       2022-02-18 [1] CRAN (R 4.2.0)
##  fansi              1.0.3      2022-03-24 [1] CRAN (R 4.2.0)
##  farver             2.1.1      2022-07-06 [1] CRAN (R 4.2.0)
##  fastmap            1.1.0      2021-01-25 [1] CRAN (R 4.2.0)
##  forcats          * 0.5.1      2021-01-27 [1] CRAN (R 4.2.0)
##  foreach            1.5.2      2022-02-02 [1] CRAN (R 4.2.0)
##  fs                 1.5.2      2021-12-08 [1] CRAN (R 4.2.0)
##  gargle             1.2.0      2021-07-02 [1] CRAN (R 4.2.0)
##  generics           0.1.3      2022-07-05 [1] CRAN (R 4.2.0)
##  GenomeInfoDb       1.32.2     2022-05-15 [1] Bioconductor
##  GenomeInfoDbData   1.2.8      2022-06-13 [1] Bioconductor
##  ggplot2          * 3.3.6      2022-05-03 [1] CRAN (R 4.2.0)
##  ggsci              2.9        2018-05-14 [1] CRAN (R 4.2.0)
##  ggsignif           0.6.3      2021-09-09 [1] CRAN (R 4.2.0)
##  glue               1.6.2      2022-02-24 [1] CRAN (R 4.2.0)
##  googledrive        2.0.0      2021-07-08 [1] CRAN (R 4.2.0)
##  googlesheets4      1.0.0      2021-07-21 [1] CRAN (R 4.2.0)
##  gtable             0.3.0      2019-03-25 [1] CRAN (R 4.2.0)
##  haven              2.5.0      2022-04-15 [1] CRAN (R 4.2.0)
##  here               1.0.1      2020-12-13 [1] CRAN (R 4.2.0)
##  highr              0.9        2021-04-16 [1] CRAN (R 4.2.0)
##  hms                1.1.1      2021-09-26 [1] CRAN (R 4.2.0)
##  htmltools          0.5.2      2021-08-25 [1] CRAN (R 4.2.0)
##  httr               1.4.3      2022-05-04 [1] CRAN (R 4.2.0)
##  igraph             1.3.2      2022-06-13 [1] CRAN (R 4.2.0)
##  IRanges            2.30.0     2022-04-26 [1] Bioconductor
##  iterators          1.0.14     2022-02-05 [1] CRAN (R 4.2.0)
##  jquerylib          0.1.4      2021-04-26 [1] CRAN (R 4.2.0)
##  jsonlite           1.8.0      2022-02-22 [1] CRAN (R 4.2.0)
##  knitr              1.39       2022-04-26 [1] CRAN (R 4.2.0)
##  lattice            0.20-45    2021-09-22 [1] CRAN (R 4.2.1)
##  lifecycle          1.0.1      2021-09-24 [1] CRAN (R 4.2.0)
##  lubridate          1.8.0      2021-10-07 [1] CRAN (R 4.2.0)
##  magrittr           2.0.3      2022-03-30 [1] CRAN (R 4.2.0)
##  MASS               7.3-58     2022-07-14 [1] CRAN (R 4.2.1)
##  Matrix             1.4-1      2022-03-23 [1] CRAN (R 4.2.1)
##  memoise            2.0.1      2021-11-26 [1] CRAN (R 4.2.0)
##  mgcv               1.8-40     2022-03-29 [1] CRAN (R 4.2.1)
##  modelr             0.1.8      2020-05-19 [1] CRAN (R 4.2.0)
##  multtest           2.52.0     2022-04-26 [1] Bioconductor
##  munsell            0.5.0      2018-06-12 [1] CRAN (R 4.2.0)
##  nlme               3.1-158    2022-06-15 [1] CRAN (R 4.2.0)
##  permute            0.9-7      2022-01-27 [1] CRAN (R 4.2.0)
##  phyloseq         * 1.40.0     2022-04-26 [1] Bioconductor
##  pillar             1.7.0      2022-02-01 [1] CRAN (R 4.2.0)
##  pkgbuild           1.3.1      2021-12-20 [1] CRAN (R 4.2.0)
##  pkgconfig          2.0.3      2019-09-22 [1] CRAN (R 4.2.0)
##  pkgload            1.3.0      2022-06-27 [1] CRAN (R 4.2.0)
##  plyr               1.8.7      2022-03-24 [1] CRAN (R 4.2.0)
##  prettyunits        1.1.1      2020-01-24 [1] CRAN (R 4.2.0)
##  processx           3.7.0      2022-07-07 [1] CRAN (R 4.2.0)
##  ps                 1.7.1      2022-06-18 [1] CRAN (R 4.2.0)
##  purrr            * 0.3.4      2020-04-17 [1] CRAN (R 4.2.0)
##  R6                 2.5.1      2021-08-19 [1] CRAN (R 4.2.0)
##  Rcpp               1.0.9      2022-07-08 [1] CRAN (R 4.2.0)
##  RCurl              1.98-1.7   2022-06-09 [1] CRAN (R 4.2.0)
##  readr            * 2.1.2      2022-01-30 [1] CRAN (R 4.2.0)
##  readxl             1.4.0      2022-03-28 [1] CRAN (R 4.2.0)
##  remotes            2.4.2      2021-11-30 [1] CRAN (R 4.2.0)
##  reprex             2.0.1      2021-08-05 [1] CRAN (R 4.2.0)
##  reshape2           1.4.4      2020-04-09 [1] CRAN (R 4.2.0)
##  rhdf5              2.40.0     2022-04-26 [1] Bioconductor
##  rhdf5filters       1.8.0      2022-04-26 [1] Bioconductor
##  Rhdf5lib           1.18.2     2022-05-17 [1] Bioconductor
##  rlang              1.0.4      2022-07-12 [1] CRAN (R 4.2.0)
##  rmarkdown          2.14       2022-04-25 [1] CRAN (R 4.2.0)
##  rprojroot          2.0.3      2022-04-02 [1] CRAN (R 4.2.0)
##  rstatix            0.7.0      2021-02-13 [1] CRAN (R 4.2.0)
##  rstudioapi         0.13       2020-11-12 [1] CRAN (R 4.2.0)
##  Rtsne              0.16       2022-04-17 [1] CRAN (R 4.2.0)
##  rvest              1.0.2      2021-10-16 [1] CRAN (R 4.2.0)
##  S4Vectors          0.34.0     2022-04-26 [1] Bioconductor
##  sass               0.4.1      2022-03-23 [1] CRAN (R 4.2.0)
##  scales             1.2.0      2022-04-13 [1] CRAN (R 4.2.0)
##  sessioninfo        1.2.2      2021-12-06 [1] CRAN (R 4.2.0)
##  speedyseq        * 0.5.3.9018 2022-06-13 [1] Github (mikemc/speedyseq@ceb941f)
##  stringi            1.7.8      2022-07-11 [1] CRAN (R 4.2.0)
##  stringr          * 1.4.0      2019-02-10 [1] CRAN (R 4.2.0)
##  survival           3.3-1      2022-03-03 [1] CRAN (R 4.2.1)
##  tibble           * 3.1.7      2022-05-03 [1] CRAN (R 4.2.0)
##  tidyr            * 1.2.0      2022-02-01 [1] CRAN (R 4.2.0)
##  tidyselect         1.1.2      2022-02-21 [1] CRAN (R 4.2.0)
##  tidyverse        * 1.3.1.9000 2022-06-13 [1] Github (tidyverse/tidyverse@6186fbf)
##  tzdb               0.3.0      2022-03-28 [1] CRAN (R 4.2.0)
##  usethis            2.1.6      2022-05-25 [1] CRAN (R 4.2.0)
##  utf8               1.2.2      2021-07-24 [1] CRAN (R 4.2.0)
##  V8                 4.2.0      2022-05-14 [1] CRAN (R 4.2.0)
##  vctrs              0.4.1      2022-04-13 [1] CRAN (R 4.2.0)
##  vegan              2.6-2      2022-04-17 [1] CRAN (R 4.2.0)
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

