---
title: "NTP72 - 16S - taxa compostiion - chicken draft "
author: "Florentin Constancias"
date: "December 21, 2022"
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
## ✔ ggplot2 3.4.0      ✔ purrr   0.3.5 
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
# source("https://raw.githubusercontent.com/fconstancias/DivComAnalyses/master/R/phyloseq_taxa_tests.R") 
source("https://raw.githubusercontent.com/fconstancias/DivComAnalyses/master/R/phyloseq_alpha.R") 
source("https://raw.githubusercontent.com/fconstancias/DivComAnalyses/master/R/phyloseq_beta.R") 
source("https://raw.githubusercontent.com/fconstancias/DivComAnalyses/master/R/phyloseq_heatmap.R")
# source("https://raw.githubusercontent.com/fconstancias/metabaRpipe-source/master/Rscripts/functions.R")
source("https://raw.githubusercontent.com/fconstancias/DivComAnalyses/master/R/phyloseq_normalisation.R")
source("https://raw.githubusercontent.com/fconstancias/DivComAnalyses/master/R/phyloseq_varia.R")


'%!in%' <- function(x,y)!('%in%'(x,y))
```



```r
# out_pptx = "~/Desktop/16S-chick.pptx"
# out_xlsx = "~/Desktop/16S-chick.xlsx"
```


# Load data:

aesthetics:


```r
load(url("https://github.com/fconstancias/NRP72-FBT/blob/master/Figures/Rastetics.Rdata?raw=true"))
```



```r
"https://github.com/fconstancias/NRP72-FBT/blob/master/data/processed/16S/16S_working_phyloseq.RDS?raw=true" %>% 
  url() %>% 
  readRDS()  %>% 
  subset_samples(Enrichment == "NotEnriched") %>%
  subset_samples(Experiment %in% c("Continuous", "Cecum")) %>%
  subset_samples(Day_of_Treatment > -4 & Day_of_Treatment <= 30 | Reactor == "DONOR") %>% 
  subset_samples(Model == "Chicken") %>% 
  subset_samples(Reactor != "CR" & Model2 == "Chicken1" | Model2 == "Chicken2" & Reactor %in% c("TR2", "TR3", "TR4", "TR5", "TR6", "TR7", "CR", "DONOR")) %>% 
  subset_samples(Reactor != "IR") %>% 
  microViz::ps_mutate("Model2_reactor_treatment_dose" = paste0(Model2,"_",Reactor_Treatment_Dose))-> ps_filtered
```

# Filtering:


```r
min_sample_size = 3738

ps_filtered %>%
  rarefy_even_depth(sample.size = min_sample_size, rngseed = 123)  -> ps_rare
```

```
## `set.seed(123)` was used to initialize repeatable random subsampling.
```

```
## Please record this for your records so others can reproduce.
```

```
## Try `set.seed(123); .Random.seed` for the full vector
```

```
## ...
```

```
## 9 samples removedbecause they contained fewer reads than `sample.size`.
```

```
## Up to first five removed samples are:
```

```
## TR1-15-S168TR1-16-S199TR2-34-S302TR3-14-S136TR3-28-S162
```

```
## ...
```

```
## 658OTUs were removed because they are no longer 
## present in any sample after random subsampling
```

```
## ...
```

```r
ps_filtered %>%
  prune_samples(sample_sums(.)>= min_sample_size, .) -> ps_fil
```

# Agglomeration at Genus level:

## phyloseq:


```r
ps_rare %>% 
  # prune_samples(get_variable(chick2, var) %in% group,
  #               .)  %>%
  filter_taxa(function(x) sum(x > 0) > 0, TRUE) %>% 
  transform_sample_counts(function(x) x/sum(x) * 100) %>% 
  physeq_glom_rename(taxrank = "Genus", rename_ASV = "Genus") -> ps_filt
```


## long dataframe:


```r
ps_filt %>% 
  transform_sample_counts(function(x) x/sum(x) * 100) %>% 
  speedyseq::psmelt() %>% 
  select(OTU, Abundance, Model2, Day_of_Treatment, Model2_reactor_treatment_dose ,Day_of_Treatment, Treatment, Reactor_Treatment_Dose) %>% 
  arrange(Day_of_Treatment)  -> ps_filt2
```


# compute count/ mean per reactors:


```r
ps_filt2 %>% 
  group_by(OTU, Model2_reactor_treatment_dose) %>% 
  add_count() %>% 
  select(-Abundance) %>% 
  distinct(OTU, Model2_reactor_treatment_dose, n) %>% 
  mutate(OTU_Reactor_Treatment_Dose = paste0(OTU,"_",Model2_reactor_treatment_dose)) -> count_OTU_reactors


ps_filt2 %>% 
  group_by(OTU, Model2_reactor_treatment_dose) %>% 
  summarise(mean_ab = mean(Abundance)) %>% 
  mutate(OTU_Reactor_Treatment_Dose = paste0(OTU,"_",Model2_reactor_treatment_dose))  -> mean_ab_OTU_reactors
```

```
## `summarise()` has grouped output by 'OTU'. You can override using the `.groups`
## argument.
```


# define treatment phases:


```r
#define_groups:
ps_filt2 %>% 
  # select(Day_of_Treatment) %>% 
  distinct(Day_of_Treatment) %>% 
  mutate(period = case_when(Day_of_Treatment <= 0 ~ "pret",
                            Day_of_Treatment > 0 &  Day_of_Treatment <= 5 ~ "t1" ,
                            Day_of_Treatment > 5 &  Day_of_Treatment <= 10 ~ "t2",
                            Day_of_Treatment > 10  ~ "t3" )) -> day_period
```


# perform classification of OTU (Genera) per reactors:


```r
ps_filt2 %>% 
  mutate(OTU_Reactor_Treatment_Dose = paste0(OTU,"_",Model2_reactor_treatment_dose)) %>% 
  left_join(count_OTU_reactors %>%  ungroup() %>%  select(-OTU, -Model2_reactor_treatment_dose),
            by = c("OTU_Reactor_Treatment_Dose" = "OTU_Reactor_Treatment_Dose"),
            keep = FALSE) %>% 
  left_join(mean_ab_OTU_reactors %>%  ungroup() %>%  select(-OTU, -Model2_reactor_treatment_dose),
            by = c("OTU_Reactor_Treatment_Dose" = "OTU_Reactor_Treatment_Dose"),
            keep = FALSE) %>% 
  # filter(mean_ab > 1 & n > 8)  %>%
  mutate(Day_of_Treatment = as.double(Day_of_Treatment)) %>%
  left_join(day_period,
            by = c("Day_of_Treatment" = "Day_of_Treatment")) %>% 
  group_by(Model2_reactor_treatment_dose, OTU, period) %>% 
  summarise(mean = base::mean(Abundance, na.rm = TRUE)) %>% 
  group_by(Model2_reactor_treatment_dose, OTU) %>% 
  pivot_wider(names_from = period, 
              values_from = mean) %>% 
summarise(category = 
            case_when(
              # pret == 0 & (t1 >0 & t2>0 & t3 >0) ~ "new",
              (t1 > pret  &
                 t2 > pret &
                 t3 > pret) & (t3 > t2 | t3 > t1 ) &
                (t3/(pret + 0.00000000001) > 1.5 & t2/(pret + 0.00000000001) > 1.5)  ~ "prol",
              (pret + 0.00000000001)/t1 > 1.5  &
                t1/t2 > 1.5 &
                t1/t3 > 1.5 ~ "dis",
              pret/t1 > 1.5 &
                t3/t1 > 1.25  ~ "recovering",
              (t1/(pret + 0.00000000001) > 1.5 | t2/(pret + 0.00000000001) > 1.5) &
                (t1/t3 > 1.5 | t2/t3 > 1.5) ~ "mid",
              # pret == 0 & t1>0 &t2>0 &t3>0 ~ "new",
              TRUE ~ "normal")) %>% 
  arrange(Model2_reactor_treatment_dose,category) %>% 
    rename("cat_period" = "category") %>% 
  mutate(OTU_Reactor_Treatment_Dose = paste0(OTU,"_",Model2_reactor_treatment_dose)) -> tmp
```

```
## `summarise()` has grouped output by 'Model2_reactor_treatment_dose', 'OTU'. You
## can override using the `.groups` argument.
## `summarise()` has grouped output by 'Model2_reactor_treatment_dose'. You can
## override using the `.groups` argument.
```

# Combine info:


```r
ps_filt2 %>%
  select(OTU, Abundance, Day_of_Treatment, Model2 ,Model2_reactor_treatment_dose) %>% 
  left_join(tmp,
            by = c("Model2_reactor_treatment_dose", "OTU")) %>%
  left_join(mean_ab_OTU_reactors,
            by = c("Model2_reactor_treatment_dose", "OTU")) %>%
  # filter(mean_ab <= 5) %>%
  filter(cat_period != "normal") %>%
  # Reactor_Treatment_Dose %in% ("TR4_CCUG59168")) %>% 
  arrange(Model2_reactor_treatment_dose) -> pre_plot
```

# generate colors:


```r
col <- c(Lachnospira = "#E6BD7F", Negativibacillus = "#3AE6C1", Tyzzerella = "#E071F6", 
         Blautia = "#568243", `Escherichia-Shigella` = "#B97CE5", `Defluviitaleaceae UCG-011` = "#C7FDC0", 
         `NK4A214 group` = "#AAEAFD", `Lachnospiraceae UCG-010` = "#F9A5FF", 
         `[Eubacterium] hallii group` = "#B76D70", Enterococcus = "#8C73EE", 
         `GCA-900066575` = "#5FE672", Monoglobus = "#99AEE0", Shuttleworthia = "#F8E1C5", 
         Butyricicoccus = "#B2ED9A", `CAG-56` = "#A7D2FB", Holdemania = "#7295FA", 
         Gordonibacter = "#EBC5E9", `Lachnospiraceae FCS020 group` = "#CD27CD", 
         Lactobacillus = "#D6DB8F", Ligilactobacillus = "#BFA6F1", Marvinbryantia = "#A8BF92", 
         Paludicola = "#EB743D", Pseudoflavonifractor = "#CCAFE1", Pygmaiobacter = "#90F3C8", 
         UBA1819 = "#D597D2", `UCG-005` = "#58EE4F", Eisenbergiella = "#C7EDD2", 
         Intestinimonas = "#8848FD", Oscillibacter = "#BBF47B", Erysipelatoclostridium = "#3A8477", 
         Oscillospira = "#F3F7F7", `[Eubacterium] nodatum group` = "#958CE9", 
         Pseudomonas = "#F5A75A", Anaerostipes = "#D67764", DTU089 = "#4AA2E4", 
         Ruminococcus = "#E281E6", Stenotrophomonas = "#EB96A7", Subdoligranulum = "#D844C3", 
         `[Eubacterium] brachy group` = "#3F88DA", Harryflintia = "#BE4EF0", 
         Anaerofilum = "#6ECCFD", Bacillus = "#CAEC53", CHKCI001 = "#9EB0FA", 
         `UC5-1-2E3` = "#B0CEE0", Lachnoclostridium = "#CDFDF1", `[Ruminococcus] torques group` = "#81D7B3", 
         Sellimonas = "#EBBAC0", Flavonifractor = "#D9C4C4", Colidextribacter = "#D0A893", 
         CHKCI002 = "#F490B6", Acinetobacter = "#E1D14D", `Family XIII AD3011 group` = "#9D48CC", 
         `Candidatus Soleaferrea` = "#B89C9B", Aerococcus = "#1CBA6F", 
         Merdibacter = "#E9BD9C", Parabacteroides = "#A0F789", Streptococcus = "#EF6FC8", 
         Prevotella_9 = "#5F5AED", Angelakisella = "#C5CBDB", Campylobacter = "#C6CCB8", 
         Catenibacterium = "#55B7A2", `Christensenellaceae R-7 group` = "#913298", 
         Cosenzaea = "#07DFA0", `Lachnospiraceae NK3A20 group` = "#EA585E", 
         Lactiplantibacillus = "#E7B0E1", Latilactobacillus = "#F8E9F8", 
         Limosilactobacillus = "#D5DE73", Megasphaera = "#F8EBAE", Sutterella = "#81D184", 
         Proteus = "#63F5F6", Lactococcus = "#FBAE9A", Oceanobacillus = "#8C7085", 
         Pediococcus = "#EDD484", Staphylococcus = "#60D8E4", Veillonella = "#E34172", 
         Tuzzerella = "#93EA46", Parasutterella = "#872253", Alistipes = "#CCE98C", 
         Anaerococcus = "#ACB54A", Collinsella = "#C6F593", Dialister = "#D8DF53", 
         Dorea = "#83B92D", Leuconostoc = "#7B6392", Roseburia = "#F5F6E3", 
         Fournierella = "#54B0B8", Anaerofustis = "#5985A1", Coprococcus = "#6F80A4", 
         Intestinibacter = "#9CE2E3", Peptococcus = "#E0F2C9", Frisingicoccus = "#B8BFE4", 
         HT002 = "#82E5DD", Methanosphaera = "#CF79BA", Prevotella = "#58F0E3", 
         `UCG-002` = "#DAE61B", Clostridioides = "#99CEB2", `Clostridium sensu stricto 1` = "#C4B6CC", 
         Listeria = "#88B9B3", Mogibacterium = "#CAAA33", Salmonella = "#B47EC8", 
         Megamonas = "#574695", Sarcina = "#C697B4")
```

# visualize - everything but normal:


```r
pre_plot %>%
  ggplot(data = ., aes_string(color = "OTU", fill = "OTU")) +
  geom_area(aes_string(x = "Day_of_Treatment", y = "Abundance")) + 
  facet_grid(cat_period ~  Model2 + Model2_reactor_treatment_dose, scales = "free") + 
  scale_fill_manual(values =  col) + 
  scale_color_manual(values = col) -> p_cat
```

```
## Warning: `aes_string()` was deprecated in ggplot2 3.0.0.
## ℹ Please use tidy evaluation ideoms with `aes()`
```

```r
p_cat + theme(legend.position = "none")
```

![](16S-taxa-classification_final_files/figure-html/unnamed-chunk-12-1.png)<!-- -->


# plot per period/ per treatment catagories:



```r
ps_filt2 %>%
   mutate(OTU_Reactor_Treatment_Dose = paste0(OTU,"_",Model2_reactor_treatment_dose)) %>% 
  select(OTU, Abundance, Day_of_Treatment, Model2 ,Treatment, Model2_reactor_treatment_dose, Reactor_Treatment_Dose, OTU_Reactor_Treatment_Dose) %>% 
  left_join(tmp,
            by = c("Model2_reactor_treatment_dose", "OTU")) %>%
  left_join(mean_ab_OTU_reactors,
            by = c("Model2_reactor_treatment_dose", "OTU")) %>% 
    left_join(day_period,
            by = c("Day_of_Treatment" = "Day_of_Treatment")) %>% 
  select(cat_period, Abundance, period, Treatment,Model2,  OTU_Reactor_Treatment_Dose, Model2_reactor_treatment_dose, Reactor_Treatment_Dose, OTU, Day_of_Treatment) %>% 
  group_by(Model2, Model2_reactor_treatment_dose,Treatment, period, OTU_Reactor_Treatment_Dose, cat_period) %>% 
  arrange(Model2, Model2_reactor_treatment_dose, Reactor_Treatment_Dose) %>% 
  summarise(mean_period = mean(Abundance)) %>% 
  group_by(Model2, Model2_reactor_treatment_dose, Treatment, period, cat_period) %>% 
  summarise(sum_cat = sum(mean_period)) %>% 
  ggplot(data = ., aes_string(x = "period", y = "sum_cat",  fill = "cat_period")) + #, color = "grey90") +
  geom_bar(position="fill", stat="identity") + 
  facet_grid(Model2 ~ Treatment, scales = "free", drop = TRUE, space = "free") + 
  theme(legend.position = "bottom") + ggpubr::rotate_x_text(45) + xlab(NULL) + ylab(NULL) +
  scale_fill_manual(values = ggpubr::get_palette(palette = "jco", k = 5)) -> p_cat2
```

```
## `summarise()` has grouped output by 'Model2', 'Model2_reactor_treatment_dose',
## 'Treatment', 'period', 'OTU_Reactor_Treatment_Dose'. You can override using the
## `.groups` argument.
## `summarise()` has grouped output by 'Model2', 'Model2_reactor_treatment_dose',
## 'Treatment', 'period'. You can override using the `.groups` argument.
```

```r
p_cat2
```

![](16S-taxa-classification_final_files/figure-html/unnamed-chunk-13-1.png)<!-- -->

# save:




```r
ps_filt2 %>%
   mutate(OTU_Reactor_Treatment_Dose = paste0(OTU,"_",Model2_reactor_treatment_dose)) %>% 
  select(OTU, Abundance, Day_of_Treatment, Model2 ,Treatment, Model2_reactor_treatment_dose, Reactor_Treatment_Dose, OTU_Reactor_Treatment_Dose) %>% 
  left_join(tmp,
            by = c("Model2_reactor_treatment_dose", "OTU")) %>%
  left_join(mean_ab_OTU_reactors,
            by = c("Model2_reactor_treatment_dose", "OTU")) %>% 
    left_join(day_period,
            by = c("Day_of_Treatment" = "Day_of_Treatment")) %>% 
  select(cat_period, Abundance, period, Treatment,Model2,  OTU_Reactor_Treatment_Dose, Model2_reactor_treatment_dose, Reactor_Treatment_Dose, OTU, Day_of_Treatment) %>% 
  group_by(Model2, Model2_reactor_treatment_dose,Treatment, period, OTU_Reactor_Treatment_Dose, cat_period) %>% 
  arrange(Model2, Model2_reactor_treatment_dose, Reactor_Treatment_Dose) %>% 
  distinct(cat_period, Abundance, period, Treatment, Model2, OTU_Reactor_Treatment_Dose, Model2_reactor_treatment_dose, Reactor_Treatment_Dose, OTU) %>% 
  write_tsv(file = "~/Desktop/16S-chick_classif.tsv")
```



```r
p_cat2$data %>% 
  write_tsv(file = "~/Desktop/16S-chick_classif_2.tsv")
```


```r
sessionInfo()
```

```
## R version 4.2.2 (2022-10-31)
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
## [1] stats     graphics  grDevices utils     datasets  methods   base     
## 
## other attached packages:
##  [1] reshape2_1.4.4  scales_1.2.1    phyloseq_1.42.0 forcats_0.5.2  
##  [5] stringr_1.4.1   dplyr_1.0.10    purrr_0.3.5     readr_2.1.3    
##  [9] tidyr_1.2.1     tibble_3.1.8    ggplot2_3.4.0   tidyverse_1.3.2
## 
## loaded via a namespace (and not attached):
##   [1] googledrive_2.0.0      colorspace_2.0-3       ggsignif_0.6.4        
##   [4] ellipsis_0.3.2         XVector_0.38.0         fs_1.5.2              
##   [7] rstudioapi_0.14        ggpubr_0.4.0           farver_2.1.1          
##  [10] bit64_4.0.5            fansi_1.0.3            lubridate_1.9.0       
##  [13] xml2_1.3.3             codetools_0.2-18       splines_4.2.2         
##  [16] cachem_1.0.6           knitr_1.40             ade4_1.7-20           
##  [19] jsonlite_1.8.3         broom_1.0.1            cluster_2.1.4         
##  [22] dbplyr_2.2.1           compiler_4.2.2         httr_1.4.4            
##  [25] backports_1.4.1        assertthat_0.2.1       Matrix_1.5-1          
##  [28] fastmap_1.1.0          gargle_1.2.1           cli_3.4.1             
##  [31] htmltools_0.5.3        tools_4.2.2            igraph_1.3.5          
##  [34] gtable_0.3.1           glue_1.6.2             GenomeInfoDbData_1.2.9
##  [37] Rcpp_1.0.9             carData_3.0-5          Biobase_2.56.0        
##  [40] cellranger_1.1.0       jquerylib_0.1.4        vctrs_0.5.0           
##  [43] Biostrings_2.66.0      rhdf5filters_1.10.0    multtest_2.54.0       
##  [46] ape_5.6-2              nlme_3.1-160           iterators_1.0.14      
##  [49] xfun_0.34              rvest_1.0.3            timechange_0.1.1      
##  [52] lifecycle_1.0.3        rstatix_0.7.1          googlesheets4_1.0.1   
##  [55] zlibbioc_1.44.0        MASS_7.3-58.1          vroom_1.6.0           
##  [58] ragg_1.2.4             hms_1.1.2              parallel_4.2.2        
##  [61] biomformat_1.26.0      rhdf5_2.42.0           yaml_2.3.6            
##  [64] speedyseq_0.5.3.9018   sass_0.4.2             stringi_1.7.8         
##  [67] highr_0.9              S4Vectors_0.36.0       foreach_1.5.2         
##  [70] permute_0.9-7          BiocGenerics_0.44.0    GenomeInfoDb_1.34.3   
##  [73] rlang_1.0.6            pkgconfig_2.0.3        systemfonts_1.0.4     
##  [76] bitops_1.0-7           evaluate_0.18          lattice_0.20-45       
##  [79] Rhdf5lib_1.20.0        labeling_0.4.2         bit_4.0.4             
##  [82] tidyselect_1.2.0       ggsci_2.9              plyr_1.8.8            
##  [85] magrittr_2.0.3         R6_2.5.1               IRanges_2.32.0        
##  [88] generics_0.1.3         DBI_1.1.3              pillar_1.8.1          
##  [91] haven_2.5.1            withr_2.5.0            mgcv_1.8-41           
##  [94] abind_1.4-5            survival_3.4-0         RCurl_1.98-1.9        
##  [97] car_3.1-1              modelr_0.1.9           crayon_1.5.2          
## [100] utf8_1.2.2             tzdb_0.3.0             rmarkdown_2.18        
## [103] grid_4.2.2             readxl_1.4.1           data.table_1.14.4     
## [106] vegan_2.6-4            reprex_2.0.2           digest_0.6.30         
## [109] textshaping_0.3.6      stats4_4.2.2           munsell_0.5.0         
## [112] microViz_0.9.7         bslib_0.4.1
```
