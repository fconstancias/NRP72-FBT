---
title: "Resistome data chicken1 - draft"
author: "Florentin Constancias "
date: " September 28, 2022 "
output: 
  html_document: 
    toc: yes
    keep_md: yes
---




```r
rm(list= ls())

gc()
```

```
##          used (Mb) gc trigger (Mb) limit (Mb) max used (Mb)
## Ncells 525950 28.1    1168072 62.4         NA   669282 35.8
## Vcells 991163  7.6    8388608 64.0      16384  1839954 14.1
```

```r
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
## ✔ tibble  3.1.8     ✔ dplyr   1.0.9
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
require(microViz)
```

```
## Loading required package: microViz
```

```
## 
## microViz version 0.9.2 - Copyright (C) 2021 David Barnett
## * Website: https://david-barnett.github.io/microViz/
## * Useful? For citation info, run: citation('microViz')
## * Silence: suppressPackageStartupMessages(library(microViz))
```

```r
'%!in%' <- function(x,y)!('%in%'(x,y))

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
source("https://raw.githubusercontent.com/fconstancias/DivComAnalyses/master/R/phyloseq_varia.R")
source("https://raw.githubusercontent.com/fconstancias/metabaRpipe-source/master/Rscripts/functions.R")
```



```r
export_ppt_width = NULL
export_ppt_height = NULL
# out_pptx = here::here("output/plots_ISME.pptx")
out_pptx = "~/Desktop/plots_manuscript_chicken1.pptx"
```

# Load data:

## 16S:


```r
# "data/processed/16S/ps_silva_dada2_human_chicken_all_fct_polyferms.RDS" %>% 
#   here::here() %>% 
#   readRDS() %>% 
#   filter_taxa(function(x) sum(x > 0) > 0, TRUE) -> ps_16S
```

## resistome:


```r
"data/processed/exp1_ps_AMR.RDS" %>% 
  here::here() %>% 
  readRDS() -> ps_AMR
```



```r
# source("https://raw.githubusercontent.com/fconstancias/metabaRpipe-source/master/Rscripts/functions.R")
# 
#   as(tax_table(ps_AMR), "matrix") %>% 
#     as.data.frame() %>%
#     rownames_to_column('ASV') -> tax_table
#   
#   as(otu_table(ps_AMR), "matrix") %>% 
#     as.data.frame() %>%
#     rownames_to_column('ASV') -> otu_table
#   
#   otu_table %>%
#     left_join(tax_table, by = c("ASV" = "ASV")) %>%
#     write_tsv("~/Desktop/TEMP/NRP72/ps_AMR.tsv")
```

## Aestetics:


```r
load( here::here("Figures/Rastetics.Rdata"))

# remove donor from legends - we better not show this guy:
treat_col <- treat_col[! treat_col %in% c('#3CB371')]
```


```r
# ps_AMR %>% 
#   tax_table() %>% 
#   data.frame()
```

## Others:


```r
ps_AMR %>% 
  tax_table() %>% 
  data.frame() %>% 
  rownames_to_column("ORF_ID") -> tax_mapping
```



```r
# ps_AMR %>% 
#   sample_data() %>% 
#   data.frame() %>% 
#   rownames_to_column("sample_id") %>% 
#     # filter(Model == "Chicken") %>% 
#   mutate(Day_of_Treatment_num = as.double(Day_of_Treatment)) %>% 
#   # mutate(Day_of_Treatment = fct_relevel(Day_of_Treatment, c(NA, "-3", "-1", "1" , "3", "4", "5","7", "30", "47", "48")))  %>% 
#   mutate(Antibiotic_mg.L = factor(Antibiotic_mg.mL, levels = c(NA, 20, 200,90, 600))) -> metadata_mapping

# metadata_mapping %>% 
#   write_tsv("~/Desktop/metadata_mapping.tsv")

"data/raw/metadata_mapping_AMR.tsv" %>% 
  here::here() %>% 
  read_tsv() %>% 
  mutate(Antibiotic_mg.L = ifelse(Phase == "Stab", NA, Antibiotic_mg.L)) %>% 
  mutate(Antibiotic_mg.L = factor(Antibiotic_mg.L, levels = c(NA, 20, 200,90, 600))) -> metadata_mapping
```

```
## Rows: 98 Columns: 42
## ── Column specification ────────────────────────────────────────────────────────
## Delimiter: "\t"
## chr (15): sample_id, Description2, Experiment, Reactor, Treatment, Enrichmen...
## dbl (26): Day_of_Connection, Day_of_Treatment, Day_from_Inoculum, GeneCopyNu...
## lgl  (1): Paul
## 
## ℹ Use `spec()` to retrieve the full column specification for this data.
## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.
```

No need to melt for this one since every otu is a gene with AMR hit and MGE / MAG ...


```r
ps_AMR %>%
  physeq_sel_tax_table(c("PathoFact_AMR_AMR_sub_class", "PathoFact_AMR_AMR_sub_class_multi",
                         "PathoFact_AMR_Resistance_mechanism","PathoFact_AMR_Resistance_mechanism_multi",
                         "PathoFact_AMR_ARG", "PathoFact_AMR_Database", "PathoFact_AMR_MGE_prediction_clean_2",
                         "viralverify_Prediction", "isescan_type", "isescan_count_IS_contig" , "ANVIO_CONCOCT_HQ_clean", "CONCOCT_bin",
                         "GTDB_tax_CONCOCT_clean", "PathoFact_AMR_MGE_prediction_clean_2")) %>% 
  subset_samples(Model == "Chicken") -> ps_AMR_all
```

```
## Note: Using an external vector in selections is ambiguous.
## ℹ Use `all_of(tax_sel)` instead of `tax_sel` to silence this message.
## ℹ See <https://tidyselect.r-lib.org/reference/faq-external-vector.html>.
## This message is displayed once per session.
```

```r
ps_AMR_all
```

```
## phyloseq-class experiment-level object
## otu_table()   OTU Table:          [ 3698 taxa and 49 samples ]:
## sample_data() Sample Data:        [ 49 samples by 39 sample variables ]:
## tax_table()   Taxonomy Table:     [ 3698 taxa by 13 taxonomic ranks ]:
## taxa are rows
```

Looks good, let's go.

# Total Resistome :

## Abundance - sum Rnum_Gi:


```r
ps_AMR_all %>%
  # physeq_sel_tax_table(c("PathoFact_AMR_AMR_sub_class", "PathoFact_AMR_AMR_sub_class_multi", "PathoFact_AMR_ARG"))  %>% 
  physeq_sel_tax_table(c("PathoFact_AMR_AMR_sub_class", "PathoFact_AMR_AMR_sub_class_multi", "PathoFact_AMR_ARG"))  %>% 
  speedyseq::tax_glom(., taxrank = "PathoFact_AMR_ARG")  %>% 
  filter_taxa(function(x) sum(x > 0) > 0, TRUE) %>% 
  speedyseq::psmelt() %>% 
  select(-Lactose_mM:-Total_SCFA_mM) %>% 
  filter(Day_of_Treatment < 34 &
           !Sample %in% c("D3", "D19")) %>% 
  filter(Reactor != "CR" & Model == "Chicken" | Model == "Human" & Reactor %in% c("TR1", "TR2", "TR3", "TR4", "TR5", "TR6", "CR"))  %>% 
  group_by(Sample) %>% 
  summarise(sum_abundance = sum(Abundance), .groups = "drop") %>% 
  mutate(sum_abundance = na_if(sum_abundance, 0)) %>%
  left_join(metadata_mapping,
            by = c("Sample" = "sample_id")) %>% 
  # mutate(Antibiotic_mg.mL = factor(Antibiotic_mg.mL, levels = c(NA, 20, 200,90, 600))) %>%
  ggplot(aes_string(x = "Day_of_Treatment", y = "sum_abundance", color = "Treatment", fill  = "Treatment", shape = "Antibiotic_mg.L")) + #shape = Fermentation
  geom_point(size = 2.5) + 
  geom_line(linetype = 2,  size = 0.5,  aes_string(group = "interaction(Model,  Reactor_Treatment_Dose)")) +
  labs(y = "sum Rnum_Gi", x = "Day_of_Treatment") +
  scale_color_manual(name = "", values = treat_col,
                     na.value = "black") +
  scale_fill_manual(name = "", values = treat_col,
                    na.value = "black") +
  scale_shape_manual(name = "" ,values = c(15,16,18,19), na.value =  17) +
  # scale_alpha_continuous(name = "", range=c(0.6, 1), na.value =  1) +
  geom_vline(xintercept = 0, size = 0.5, linetype = 1, color = "grey60") +
  theme_light() + 
  scale_y_continuous(labels=scientific) -> p1 # scale_y_continuous(labels = function(x) format(x, scientific = TRUE))
```



```r
p1 + facet_grid(as.formula(paste0("Antibiotic  ~  Model  ")), scales = "free", space = "fixed", drop = TRUE) -> ptmp

ptmp
```

<img src="chicken1_resistome_plots_draft_files/figure-html/unnamed-chunk-12-1.png" style="display: block; margin: auto;" />

```r
ptmp %>% 
  ggsave(plot = .,bg = "transparent",scale = 2,
         units = "mm",
         dpi = 600,
         height = 0.618 * 84,
         width = 84,
         filename = "~/Desktop/Resistome_abundance.eps")

ptmp %>%
  export::graph2ppt(append = TRUE, width = 317.48031496 , height = 0.618 * 317.48031496, paper = "A4",  scaling = 2,
                    file = out_pptx)
```

```
## Exported graph as ~/Desktop/plots_manuscript_chicken1.pptx
```

Also retrieve the samples we now keep for the final analysis:


```r
p1$data %>% distinct(Sample)  %>% pull %>%  as.character() -> samples
```

Export values displayed:


```r
p1$data %>% 
  select(Sample, sum_abundance, Day_of_Treatment, Treatment, Model) %>% 
  arrange(Day_of_Treatment) %>% 
  xlsx::write.xlsx(x = . ,
                   file = "~/Desktop/data_plots.xlsx",
                   sheetName = "resistome_abundance")
```



## / Baseline reference :


```r
#https://stackoverflow.com/questions/52581469/dynamically-normalize-all-rows-with-first-element-within-a-group
#<https://blog.exploratory.io/introducing-time-series-analysis-with-dplyr-60683587cf8a>

ps_AMR_all %>% 
  physeq_sel_tax_table(c("PathoFact_AMR_AMR_sub_class", "PathoFact_AMR_AMR_sub_class_multi", "PathoFact_AMR_ARG"))  %>% 
  speedyseq::tax_glom(., taxrank = "PathoFact_AMR_ARG")  %>% 
  prune_samples(samples, .) %>% 
  filter_taxa(function(x) sum(x > 0) > 0, TRUE) %>% 
  speedyseq::psmelt() %>% 
  select(-Lactose_mM:-Total_SCFA_mM) %>% 
  group_by(Sample) %>% 
  summarise(sum_abundance = sum(Abundance), .groups = "drop") %>% 
  mutate(sum_abundance = na_if(sum_abundance, 0)) %>%
  left_join(metadata_mapping,
            by = c("Sample" = "sample_id")) %>% 
  ungroup() %>% 
  group_by(Model, Reactor_Treatment_Dose) %>% 
  arrange(Day_of_Treatment) %>% 
  mutate(per_cent_change = sum_abundance / first(sum_abundance)) %>% 
  ggplot(aes_string(x = "Day_of_Treatment", y = "per_cent_change", color = "Treatment", fill  = "Treatment", shape = "Antibiotic_mg.L")) +
  geom_point(size = 2.5) + 
  geom_line(linetype = 2,  size = 0.5,  aes_string(group = "interaction(Model, Reactor_Treatment_Dose)")) +
  labs(y = "sum Rnum_Gi / baseline", x = "Day_of_Treatment") +
  scale_color_manual(name = "", values = treat_col,
                     na.value = "black") +
  scale_fill_manual(name = "", values = treat_col,
                    na.value = "black") +
  scale_shape_manual(name = "" ,values = c(15,16,18,19), na.value =  17) +
  geom_vline(xintercept = 0, size = 0.5, linetype = 1, color = "grey60") +
  geom_hline(yintercept = 1, size = 0.5, linetype = 2, color = "grey60") +
  theme_light() + facet_grid(as.formula(paste0("Antibiotic ~  Model  ")), 
                             scales = "free_x", space = "fixed", drop = TRUE) -> p1

p1
```

<img src="chicken1_resistome_plots_draft_files/figure-html/unnamed-chunk-15-1.png" style="display: block; margin: auto;" />

```r
p1 %>% 
  ggsave(plot = .,bg = "transparent",scale = 2,
         units = "mm",
         dpi = 600,
         height = 0.618 * 84,
         width = 84,
         filename = "~/Desktop/Resistome_baseline.eps")

p1 %>%
  export::graph2ppt(append = TRUE, width = 317.48031496 / 1.25, height = 0.618 * 317.48031496  / 1.25, paper = "A4",  scaling = 2,
                    file = out_pptx)
```

```
## Exported graph as ~/Desktop/plots_manuscript_chicken1.pptx
```

Export values displayed:



```r
p1$data %>% 
  data.frame() %>% 
  select(Sample, per_cent_change, Day_of_Treatment, Treatment, Model) %>% 
  # arrange(Day_of_Treatment) %>% 
  xlsx::write.xlsx(x = . ,
                   append = TRUE,
                   file = "~/Desktop/data_plots.xlsx",
                   sheetName = "resistomebaseline")
```


# AMR drug classes


```r
ps_AMR_all %>% 
  physeq_sel_tax_table(c("PathoFact_AMR_AMR_sub_class", "PathoFact_AMR_AMR_sub_class_multi", "PathoFact_AMR_ARG"))  %>%
  speedyseq::tax_glom(., taxrank = "PathoFact_AMR_AMR_sub_class_multi")  %>% 
  physeq_sel_tax_table(c("PathoFact_AMR_AMR_sub_class_multi")) %>% 
  prune_samples(samples, .) %>% 
  filter_taxa(function(x) sum(x > 0) > 0, TRUE) %>% 
  physeq_add_metadata(metadata_mapping,
                      sample_column = "sample_id") %>% 
  transform_sample_counts(function(x) x/sum(x) * 100) %>% 
  subset_taxa(PathoFact_AMR_AMR_sub_class_multi %!in% c("Unclassified", "multidrug (unclassified)")) %>% 
  # ps_mutate(Day_of_Treatment = as.numeric(Day_of_Treatment)) %>% 
  ps_arrange(Day_of_Treatment) %>% 
  comp_barplot(
    palette = distinct_palette(12, pal = "kelly"),
    tax_level = "PathoFact_AMR_AMR_sub_class_multi",
    label = "Day_of_Treatment", # name an alternative variable to label axes
    n_taxa = 12, # give more taxa unique colours
    merge_other = FALSE, # split the "other" category to display alpha diversity
    bar_width = 0.9, # reduce the bar width to 70% of one row
    bar_outline_colour = "grey5",
    sample_order = "default") +#,# is the default (use NA to remove outlines)
  # facet_by = c("Model", "Antibiotic", "Antibiotic_mg.mL", "Treatment2")) +
  # group_by = c("Model", "Antibiotic", "Antibiotic_mg.mL", "Treatment2"))
  # coord_flip() + 
  theme_light() + 
  ylab("Proportion - %") +
  scale_y_continuous(labels=percent) + 
  facet_grid(. ~ Antibiotic + Antibiotic_mg.mL + Treatment2, scales = "free", space = "free") -> p2
```

```
## Warning in ps_counts(data, warn = TRUE): otu_table of counts is NOT available!
## Available otu_table contains non-zero values that are less than 1
```

```
## Warning in ps_counts(data = data): otu_table of counts is NOT available!
## Available otu_table contains non-zero values that are less than 1
```

```r
# facet_wrap(Model ~ Antibiotic + Antibiotic_mg.mL + Treatment2, scales = "free", ncol = 6, nrow = 4) -> p2

# # Plot them side by side with the patchwork package.
# patch <- patchwork::wrap_plots(plot_list, nrow = 1, guides = "collect")
# patch & coord_flip() # make all plots horizontal (note: use & instead of +)
```



```r
plot_AMR_multi <- p2


p2 
```

<img src="chicken1_resistome_plots_draft_files/figure-html/unnamed-chunk-18-1.png" style="display: block; margin: auto;" />

```r
p2 %>% 
  ggsave(plot = .,bg = "transparent",scale = 2,
         units = "mm",
         dpi = 600,
         height = 0.618 * 84 * 2,
         width = 84 * 2,
         filename = "~/Desktop/Resistome_sub_class_multi_prop.eps")

p2 %>%
  export::graph2ppt(append = TRUE, width = 317.48031496 * 1.5  , height = 0.618 * 317.48031496  , paper = "A4",  scaling = 2,
                    file = out_pptx)
```

```
## Exported graph as ~/Desktop/plots_manuscript_chicken1.pptx
```


Export values displayed:



```r
p2$data %>% 
  data.frame() %>% 
  select(Sample, Abundance, Day_of_Treatment, Treatment, Model, OTU) %>% 
  # arrange(Day_of_Treatment) %>% 
  xlsx::write.xlsx(x = . ,
                   append = TRUE,
                   file = "~/Desktop/data_plots.xlsx",
                   sheetName = "drug_class")
```

same but heatmap:


```r
ps_AMR_all %>% 
  physeq_sel_tax_table(c("PathoFact_AMR_AMR_sub_class_multi"))  %>%
  speedyseq::tax_glom(., taxrank = "PathoFact_AMR_AMR_sub_class_multi")  %>% 
  physeq_sel_tax_table(c("PathoFact_AMR_AMR_sub_class_multi")) %>% 
  prune_samples(samples, .) %>% 
  filter_taxa(function(x) sum(x > 0) > 0, TRUE) %>% 
  physeq_add_metadata(metadata_mapping,
                      sample_column = "sample_id") %>% 
  transform_sample_counts(function(x) x/sum(x) * 100) %>% 
  subset_taxa(PathoFact_AMR_AMR_sub_class_multi %!in% c("Unclassified", "multidrug (unclassified)", "MLS (unclassified)")) %>% 
  # ps_mutate(Day_of_Treatment = as.numeric(Day_of_Treatment)) %>% 
  # ps_arrange(-Day_of_Treatment) %>% 
  speedyseq::psmelt() %>%
  # left_join(tax_mapping,
  #           by = c("Best_Hit_ARO" = "Best_Hit_ARO"),
  #           suffix = c("_x", "")) %>%
  # mutate(Sample = fct_reorder(Sample, Model,Treatment,Fermentation,Antibiotic_mg.mL, Day_of_Treatment )) %>%
  mutate(Abundance = na_if(Abundance, 0)) %>%
  mutate(logTPM = log10(Abundance + 1)) -> pheatmap_df


pheatmap_df %>%
  # replace(is.na(.), 0)  %>%
  mutate(log10Abundance = log10(Abundance)) %>% 
  select(Sample, PathoFact_AMR_AMR_sub_class_multi, Abundance) %>%
  group_by(Sample, PathoFact_AMR_AMR_sub_class_multi) %>%
  pivot_wider(names_from =  Sample, PathoFact_AMR_AMR_sub_class_multi,
              values_from = Abundance) %>%
  column_to_rownames("PathoFact_AMR_AMR_sub_class_multi") -> pheatmap_samples


pheatmap_df %>%
  # replace(is.na(.), 0)  %>%
  mutate(log10Abundance = log10(Abundance)) %>% 
  select(Sample, PathoFact_AMR_AMR_sub_class_multi, Abundance) %>%
  group_by(PathoFact_AMR_AMR_sub_class_multi) %>% 
  summarise(meanprop = mean(Abundance, na.rm = TRUE)) -> mean_ab


# we want to knownumber of AMR gene behind each category:

ps_AMR_all %>% 
  physeq_sel_tax_table(c("PathoFact_AMR_AMR_sub_class_multi", "PathoFact_AMR_ARG"))  %>%
  speedyseq::tax_glom(., taxrank = "PathoFact_AMR_ARG")  %>% 
  transform_sample_counts(function(x) x/sum(x) * 100) %>% 
  subset_taxa(PathoFact_AMR_AMR_sub_class_multi %!in% c("Unclassified", "multidrug (unclassified)", "MLS (unclassified)")) %>% 
  # physeq_sel_tax_table(c("PathoFact_AMR_AMR_sub_class_multi")) %>% 
  prune_samples(samples, .) %>% 
  filter_taxa(function(x) sum(x > 0) > 0, TRUE) %>% 
  physeq_add_metadata(metadata_mapping,
                      sample_column = "sample_id") %>% 
  # ps_mutate(Day_of_Treatment = as.numeric(Day_of_Treatment)) %>% 
  # ps_arrange(-Day_of_Treatment) %>% 
  speedyseq::psmelt() %>% 
  select(Sample, PathoFact_AMR_AMR_sub_class_multi, PathoFact_AMR_ARG, Abundance) %>%
  filter(Abundance > 0) %>% 
  select(-Sample, -Abundance) %>% 
  distinct(PathoFact_AMR_ARG, .keep_all = TRUE) %>% 
  group_by(PathoFact_AMR_AMR_sub_class_multi) %>%
  summarise(nAMR = length(PathoFact_AMR_ARG)) -> num_AMR

# we want to knownumber of gene behind each AMR gene name:

ps_AMR_all %>% 
  physeq_sel_tax_table(c("PathoFact_AMR_AMR_sub_class_multi", "PathoFact_AMR_ARG"))  %>%
  # speedyseq::tax_glom(., taxrank = "PathoFact_AMR_ARG")  %>% 
  transform_sample_counts(function(x) x/sum(x) * 100) %>% 
  subset_taxa(PathoFact_AMR_AMR_sub_class_multi %!in% c("Unclassified", "multidrug (unclassified)", "MLS (unclassified)")) %>% 
  # physeq_sel_tax_table(c("PathoFact_AMR_AMR_sub_class_multi")) %>% 
  prune_samples(samples, .) %>% 
  filter_taxa(function(x) sum(x > 0) > 0, TRUE) %>% 
  physeq_add_metadata(metadata_mapping,
                      sample_column = "sample_id") %>% 
  # ps_mutate(Day_of_Treatment = as.numeric(Day_of_Treatment)) %>% 
  # ps_arrange(-Day_of_Treatment) %>% 
  speedyseq::psmelt() %>% 
  select(Sample, PathoFact_AMR_AMR_sub_class_multi, OTU, Abundance) %>%
  filter(Abundance > 0) %>% 
  select(-Sample, -Abundance) %>% 
  distinct(OTU, .keep_all = TRUE) %>% 
  group_by(PathoFact_AMR_AMR_sub_class_multi) %>%
  summarise(ngn = length(OTU)) -> num_gn


cbind(mean_ab,
      nGN = num_gn$ngn,
      nAMR = num_AMR$nAMR) -> row_anot
```


dendrogram gN:


```r
ps_AMR_all %>% 
  physeq_sel_tax_table(c("PathoFact_AMR_AMR_sub_class_multi", "PathoFact_AMR_ARG"))  %>%
  # speedyseq::tax_glom(., taxrank = "PathoFact_AMR_ARG")  %>% 
  transform_sample_counts(function(x) x/sum(x) * 100) %>% 
  subset_taxa(PathoFact_AMR_AMR_sub_class_multi %!in% c("Unclassified", "multidrug (unclassified)", "MLS (unclassified)")) %>% 
  # physeq_sel_tax_table(c("PathoFact_AMR_AMR_sub_class_multi")) %>% 
  prune_samples(samples, .) %>% 
  filter_taxa(function(x) sum(x > 0) > 0, TRUE) %>% 
  physeq_add_metadata(metadata_mapping,
                      sample_column = "sample_id") %>% 
  # ps_mutate(Day_of_Treatment = as.numeric(Day_of_Treatment)) %>% 
  # ps_arrange(-Day_of_Treatment) %>% 
  speedyseq::psmelt() %>%
  # left_join(tax_mapping,
  #           by = c("Best_Hit_ARO" = "Best_Hit_ARO"),
  #           suffix = c("_x", "")) %>%
  # mutate(Sample = fct_reorder(Sample, Model,Treatment,Fermentation,Antibiotic_mg.mL, Day_of_Treatment )) %>%
  mutate(Abundance = na_if(Abundance, 0)) %>%
  # mutate(logTPM = log10(Abundance + 1)) %>% 
  # mutate(log10Abundance = log10(Abundance)) %>% 
  select(Sample, OTU, Abundance) %>%
  group_by(Sample, OTU) %>%
  pivot_wider(names_from =  Sample, OTU,
              values_from = Abundance) %>%
  column_to_rownames("OTU") %>% 
  replace(is.na(.), 0)  %>%
  t() %>%
  vegan::vegdist(na.rm = TRUE) %>%
  hclust(method = "ward.D2") -> clust_t_OTU
```

dendrogram AMRgN:


```r
ps_AMR_all %>% 
  physeq_sel_tax_table(c("PathoFact_AMR_AMR_sub_class_multi", "PathoFact_AMR_ARG"))  %>%
  speedyseq::tax_glom(., taxrank = "PathoFact_AMR_ARG")  %>% 
  transform_sample_counts(function(x) x/sum(x) * 100) %>% 
  subset_taxa(PathoFact_AMR_AMR_sub_class_multi %!in% c("Unclassified", "multidrug (unclassified)", "MLS (unclassified)")) %>% 
  # physeq_sel_tax_table(c("PathoFact_AMR_AMR_sub_class_multi")) %>% 
  prune_samples(samples, .) %>% 
  filter_taxa(function(x) sum(x > 0) > 0, TRUE) %>% 
  physeq_add_metadata(metadata_mapping,
                      sample_column = "sample_id") %>% 
  # ps_mutate(Day_of_Treatment = as.numeric(Day_of_Treatment)) %>% 
  # ps_arrange(-Day_of_Treatment) %>% 
  speedyseq::psmelt() %>%
  # left_join(tax_mapping,
  #           by = c("Best_Hit_ARO" = "Best_Hit_ARO"),
  #           suffix = c("_x", "")) %>%
  # mutate(Sample = fct_reorder(Sample, Model,Treatment,Fermentation,Antibiotic_mg.mL, Day_of_Treatment )) %>%
  mutate(Abundance = na_if(Abundance, 0)) %>%
  # mutate(logTPM = log10(Abundance + 1)) %>% 
  # mutate(log10Abundance = log10(Abundance)) %>% 
  select(Sample, PathoFact_AMR_ARG, Abundance) %>%
  group_by(Sample, PathoFact_AMR_ARG) %>%
  pivot_wider(names_from =  Sample, PathoFact_AMR_ARG,
              values_from = Abundance) %>%
  column_to_rownames("PathoFact_AMR_ARG") %>% 
  replace(is.na(.), 0)  %>%
  t() %>%
  # sqrt() %>% 
  vegan::vegdist(na.rm = TRUE, method = "robust.aitchison") %>%
  hclust(method = "ward.D2") -> clust_t_AMR

# NO : median, centroid
# "ward.D", "ward.D2", "single", "complete", "average" (= UPGMA), "mcquitty" (= WPGMA), "median" (= WPGMC) or "centroid" (= UPGMC)
```



```r
pheatmap_samples %>%
  replace(is.na(.), 0)  %>%
  t() %>%
  vegan::vegdist(na.rm = TRUE, method = "robust.aitchison") %>%
  hclust(method = "ward.D2") -> clust_t


pheatmap_samples %>%
  replace(is.na(.), 0)  %>%
  vegan::vegdist(na.rm = TRUE, method = "robust.aitchison") %>%
  # cor() %>%
  hclust(method = "ward.D2") -> clust
```


Pheatmap:

dendro class


```r
pheatmap_samples %>%
  pheatmap::pheatmap(cluster_rows = FALSE ,
                     scale = "row",
                     cluster_cols = clust_t,
                     na_col = "transparent",
                     annotation_row = tax_mapping %>% rownames_to_column('tmp_col') %>%
                       filter(PathoFact_AMR_AMR_sub_class_multi %in% rownames(pheatmap_samples)) %>%
                       distinct(PathoFact_AMR_AMR_sub_class_multi, .keep_all = TRUE) %>%
                       left_join(., 
                                 row_anot,
                                 by = c("PathoFact_AMR_AMR_sub_class_multi" = "PathoFact_AMR_AMR_sub_class_multi")) %>% 
                       column_to_rownames("PathoFact_AMR_AMR_sub_class_multi") %>%
                       select(PathoFact_AMR_Resistance_mechanism_multi, meanprop, nGN, nAMR) %>% 
                       rename(mechanism = PathoFact_AMR_Resistance_mechanism_multi),
                     annotation_colors = list("mechanism" = resistance_type,
                                              "Treatment" = treat_col),
                     annotation_col = pheatmap_df %>%
                       select(Day_of_Treatment, Antibiotic_mg.mL, Treatment, Sample) %>% 
                       mutate(Day_of_Treatment = ifelse(Day_of_Treatment < 0, NA,Day_of_Treatment)) %>% 
                       distinct(Sample, .keep_all = TRUE) %>% 
                       column_to_rownames("Sample"),
                     annotation_legend = TRUE)  -> p_pheat

p_pheat
```

<img src="chicken1_resistome_plots_draft_files/figure-html/unnamed-chunk-24-1.png" style="display: block; margin: auto;" />

```r
p_pheat %>% 
  ggsave(plot = .,bg = "transparent",scale = 2,
         units = "mm",
         dpi = 600,
         height = 0.618 * 84 * 2,
         width = 84  * 1.5,
         filename = "~/Desktop/Resistome_heatmap_class.eps")

p_pheat %>%
  export::graph2ppt(append = TRUE, width = 317.48031496 * 8, height = 0.618 * 317.48031496 * 10 , paper = "A4",  scaling = 1,
                    file = out_pptx)
```

```
## Exported graph as ~/Desktop/plots_manuscript_chicken1.pptx
```

dendro AMR


```r
pheatmap_samples %>%
  pheatmap::pheatmap(cluster_rows = FALSE ,
                     scale = "row",
                     cluster_cols = clust_t_AMR,
                     na_col = "transparent",
                     annotation_row = tax_mapping %>% rownames_to_column('tmp_col') %>%
                       filter(PathoFact_AMR_AMR_sub_class_multi %in% rownames(pheatmap_samples)) %>%
                       distinct(PathoFact_AMR_AMR_sub_class_multi, .keep_all = TRUE) %>%
                       left_join(., 
                                 row_anot,
                                 by = c("PathoFact_AMR_AMR_sub_class_multi" = "PathoFact_AMR_AMR_sub_class_multi")) %>% 
                       column_to_rownames("PathoFact_AMR_AMR_sub_class_multi") %>%
                       select(PathoFact_AMR_Resistance_mechanism_multi, meanprop, nGN, nAMR) %>% 
                       rename(mechanism = PathoFact_AMR_Resistance_mechanism_multi),
                     annotation_colors = list("mechanism" = resistance_type,
                                              "Treatment" = treat_col),
                     annotation_col = pheatmap_df %>%
                       select(Day_of_Treatment, Antibiotic_mg.mL, Treatment, Sample) %>% 
                       mutate(Day_of_Treatment = ifelse(Day_of_Treatment < 0, NA,Day_of_Treatment)) %>% 
                       distinct(Sample, .keep_all = TRUE) %>% 
                       column_to_rownames("Sample"),
                     annotation_legend = TRUE)  -> p_pheat

p_pheat
```

<img src="chicken1_resistome_plots_draft_files/figure-html/unnamed-chunk-25-1.png" style="display: block; margin: auto;" />

```r
p_pheat %>% 
  ggsave(plot = .,bg = "transparent",scale = 2,
         units = "mm",
         dpi = 600,
         height = 0.618 * 84 * 2,
         width = 84  * 1.5,
         filename = "~/Desktop/Resistome_heatmap_class.eps")

p_pheat %>%
  export::graph2ppt(append = TRUE, width = 317.48031496 * 8, height = 0.618 * 317.48031496 * 10 , paper = "A4",  scaling = 1,
                    file = out_pptx)
```

```
## Exported graph as ~/Desktop/plots_manuscript_chicken1.pptx
```


dendro gN


```r
pheatmap_samples %>%
  pheatmap::pheatmap(cluster_rows = FALSE ,
                     scale = "row",
                     cluster_cols = clust_t_OTU,
                     na_col = "transparent",
                     annotation_row = tax_mapping %>% rownames_to_column('tmp_col') %>%
                       filter(PathoFact_AMR_AMR_sub_class_multi %in% rownames(pheatmap_samples)) %>%
                       distinct(PathoFact_AMR_AMR_sub_class_multi, .keep_all = TRUE) %>%
                       left_join(., 
                                 row_anot,
                                 by = c("PathoFact_AMR_AMR_sub_class_multi" = "PathoFact_AMR_AMR_sub_class_multi")) %>% 
                       column_to_rownames("PathoFact_AMR_AMR_sub_class_multi") %>%
                       select(PathoFact_AMR_Resistance_mechanism_multi, meanprop, nGN, nAMR) %>% 
                       rename(mechanism = PathoFact_AMR_Resistance_mechanism_multi),
                     annotation_colors = list("mechanism" = resistance_type,
                                              "Treatment" = treat_col),
                     annotation_col = pheatmap_df %>%
                       select(Day_of_Treatment, Antibiotic_mg.mL, Treatment, Sample) %>% 
                       mutate(Day_of_Treatment = ifelse(Day_of_Treatment < 0, NA,Day_of_Treatment)) %>% 
                       distinct(Sample, .keep_all = TRUE) %>% 
                       column_to_rownames("Sample"),
                     annotation_legend = TRUE)  -> p_pheat

p_pheat
```

<img src="chicken1_resistome_plots_draft_files/figure-html/unnamed-chunk-26-1.png" style="display: block; margin: auto;" />

```r
p_pheat %>% 
  ggsave(plot = .,bg = "transparent",scale = 2,
         units = "mm",
         dpi = 600,
         height = 0.618 * 84 * 2,
         width = 84  * 1.5,
         filename = "~/Desktop/Resistome_heatmap_class.eps")

p_pheat %>%
  export::graph2ppt(append = TRUE, width = 317.48031496 * 8, height = 0.618 * 317.48031496 * 10 , paper = "A4",  scaling = 1,
                    file = out_pptx)
```

```
## Exported graph as ~/Desktop/plots_manuscript_chicken1.pptx
```

export underlying data:


```r
tax_mapping %>% rownames_to_column('tmp_col') %>%
  filter(PathoFact_AMR_AMR_sub_class_multi %in% rownames(pheatmap_samples)) %>%
  distinct(PathoFact_AMR_AMR_sub_class_multi, .keep_all = TRUE) %>%
  left_join(., 
            row_anot,
            by = c("PathoFact_AMR_AMR_sub_class_multi" = "PathoFact_AMR_AMR_sub_class_multi")) %>% 
  column_to_rownames("PathoFact_AMR_AMR_sub_class_multi") %>%
  select(PathoFact_AMR_Resistance_mechanism_multi, meanprop, nGN, nAMR) %>% 
  rename(mechanism = PathoFact_AMR_Resistance_mechanism_multi) %>% 
  data.frame() %>% 
  # arrange(Day_of_Treatment) %>% 
  xlsx::write.xlsx(x = . ,
                   append = TRUE,
                   file = "~/Desktop/data_plots.xlsx",
                   sheetName = "drug_class_heatmap")
```

only cephalosporin:


```r
ps_AMR_all %>% 
  physeq_sel_tax_table(c("PathoFact_AMR_AMR_sub_class", "PathoFact_AMR_AMR_sub_class_multi", "PathoFact_AMR_ARG")) %>% 
  subset_taxa(PathoFact_AMR_AMR_sub_class_multi == "cephalosporin") %>% 
  speedyseq::tax_glom(., taxrank = "PathoFact_AMR_ARG")  %>% 
  physeq_sel_tax_table(c("PathoFact_AMR_ARG")) %>% 
  prune_samples(samples, .) %>% 
  physeq_add_metadata(metadata_mapping,
                      sample_column = "sample_id") %>% 
  transform_sample_counts(function(x) x/sum(x) * 100) %>% 
  filter_taxa(function(x) sum(x > 0) > 0, TRUE) %>% 
  # subset_taxa(PathoFact_AMR_AMR_sub_class_multi %!in% c("Unclassified", "multidrug (unclassified)", "MLS (unclassified)")) %>% 
  # ps_mutate(Day_of_Treatment = as.numeric(Day_of_Treatment)) %>% 
  ps_arrange(-Day_of_Treatment) %>% 
  speedyseq::psmelt() %>%
  # left_join(tax_mapping,
  #           by = c("Best_Hit_ARO" = "Best_Hit_ARO"),
  #           suffix = c("_x", "")) %>%
  # mutate(Sample = fct_reorder(Sample, Model,Treatment,Fermentation,Antibiotic_mg.mL, Day_of_Treatment )) %>%
  mutate(Abundance = na_if(Abundance, 0)) %>%
  mutate(logAbundance = log10(Abundance + 1)) %>% 
  mutate(PathoFact_AMR_ARG = PathoFact_AMR_ARG %>% stringr::str_trunc(45,  side ="center")) -> pheatmap_df


pheatmap_df %>%
  # replace(is.na(.), 0)  %>%
  mutate(log10Abundance = log10(Abundance)) %>% 
  select(Sample, PathoFact_AMR_ARG, Abundance) %>%
  filter(Abundance > 0) %>%
  # filter(Abundance > 1) %>%
  group_by(Sample, PathoFact_AMR_ARG) %>%
  pivot_wider(names_from =  Sample, PathoFact_AMR_ARG,
              values_from = Abundance) %>%
  # mutate(PathoFact_AMR_ARG = PathoFact_AMR_ARG %>%  str_to_lower(.)) %>% 
  # %>% str_to_lower(., locale = "en")
  column_to_rownames("PathoFact_AMR_ARG") -> pheatmap_samples

pheatmap_df %>%
  # replace(is.na(.), 0)  %>%
  mutate(log10Abundance = log10(Abundance)) %>% 
  select(Sample, PathoFact_AMR_ARG, Abundance) %>%
  group_by(PathoFact_AMR_ARG) %>% 
  summarise(meanAbundance = mean(Abundance, na.rm = TRUE)) -> mean_ab


pheatmap_samples %>%
  replace(is.na(.), 0)  %>%
  t() %>%
  vegan::vegdist(na.rm = TRUE) %>%
  hclust() -> clust_t

pheatmap_samples %>%
  replace(is.na(.), 0)  %>%
  vegan::vegdist(na.rm = TRUE, method = "euc") %>%
  # cor() %>%
  hclust() -> clust

pheatmap_samples %>%
  pheatmap::pheatmap(fontsize_col = 6,
                     fontsize = 6,
                     cluster_rows = FALSE ,
                     scale = "row",
                     cluster_cols = clust_t,
                     na_col = "transparent",
                     annotation_row = tax_mapping %>% rownames_to_column('tmp_col') %>%
                       filter(PathoFact_AMR_ARG %in% rownames(pheatmap_samples)) %>%
                       distinct(PathoFact_AMR_ARG, .keep_all = TRUE) %>%
                       left_join(., 
                                 mean_ab,
                                 by = c("PathoFact_AMR_ARG" = "PathoFact_AMR_ARG")) %>% 
                       column_to_rownames("PathoFact_AMR_ARG") %>%
                       select(PathoFact_AMR_AMR_sub_class_multi, meanAbundance),
                     annotation_colors = list("PathoFact_AMR_AMR_sub_class_multi" = drug_classe_multi,
                                              "Treatment" = treat_col),
                     annotation_col = pheatmap_df %>%
                       select(Antibiotic_mg.mL,Day_of_Treatment, Treatment, Model, Sample) %>% 
                       distinct(Sample, .keep_all = TRUE) %>% 
                       column_to_rownames("Sample"))  -> p_pheat


p_pheat
```

<img src="chicken1_resistome_plots_draft_files/figure-html/unnamed-chunk-28-1.png" style="display: block; margin: auto;" />

```r
p_pheat %>% 
  ggsave(plot = .,bg = "transparent",scale = 2,
         units = "mm",
         dpi = 600,
         height = 0.618 * 84,
         width = 84  ,
         filename = "~/Desktop/Resistome_heatmap_CTX.eps")


p_pheat %>%
  export::graph2ppt(append = TRUE, width = 317.48031496 * 2, height = 0.618 * 317.48031496 * 2 , paper = "A4",  scaling = 1,
                    file = out_pptx)
```

```
## Exported graph as ~/Desktop/plots_manuscript_chicken1.pptx
```


same but only ....:


```r
ps_AMR_all %>% 
  physeq_sel_tax_table(c("PathoFact_AMR_AMR_sub_class", "PathoFact_AMR_AMR_sub_class_multi", "PathoFact_AMR_ARG")) %>% 
  subset_taxa(PathoFact_AMR_AMR_sub_class_multi == "glycopeptide antibiotic") %>% 
  speedyseq::tax_glom(., taxrank = "PathoFact_AMR_ARG")  %>% 
  physeq_sel_tax_table(c("PathoFact_AMR_ARG")) %>% 
  prune_samples(samples, .) %>% 
  physeq_add_metadata(metadata_mapping,
                      sample_column = "sample_id") %>% 
  transform_sample_counts(function(x) x/sum(x) * 100) %>% 
  filter_taxa(function(x) sum(x > 0) > 0, TRUE) %>% 
  # subset_taxa(PathoFact_AMR_AMR_sub_class_multi %!in% c("Unclassified", "multidrug (unclassified)", "MLS (unclassified)")) %>% 
  # ps_mutate(Day_of_Treatment = as.numeric(Day_of_Treatment)) %>% 
  ps_arrange(-Day_of_Treatment) %>% 
  speedyseq::psmelt() %>%
  # left_join(tax_mapping,
  #           by = c("Best_Hit_ARO" = "Best_Hit_ARO"),
  #           suffix = c("_x", "")) %>%
  # mutate(Sample = fct_reorder(Sample, Model,Treatment,Fermentation,Antibiotic_mg.mL, Day_of_Treatment )) %>%
  mutate(Abundance = na_if(Abundance, 0)) %>%
  mutate(logAbundance = log10(Abundance + 1)) %>% 
  mutate(PathoFact_AMR_ARG = PathoFact_AMR_ARG %>% stringr::str_trunc(45,  side ="center")) -> pheatmap_df


pheatmap_df %>%
  # replace(is.na(.), 0)  %>%
  mutate(log10Abundance = log10(Abundance)) %>% 
  select(Sample, PathoFact_AMR_ARG, Abundance) %>%
  filter(Abundance > 0) %>%
  # filter(Abundance > 1) %>%
  group_by(Sample, PathoFact_AMR_ARG) %>%
  pivot_wider(names_from =  Sample, PathoFact_AMR_ARG,
              values_from = Abundance) %>%
  # mutate(PathoFact_AMR_ARG = PathoFact_AMR_ARG %>%  str_to_lower(.)) %>% 
  # %>% str_to_lower(., locale = "en")
  column_to_rownames("PathoFact_AMR_ARG") -> pheatmap_samples

pheatmap_df %>%
  # replace(is.na(.), 0)  %>%
  mutate(log10Abundance = log10(Abundance)) %>% 
  select(Sample, PathoFact_AMR_ARG, Abundance) %>%
  group_by(PathoFact_AMR_ARG) %>% 
  summarise(meanAbundance = mean(Abundance, na.rm = TRUE)) -> mean_ab


pheatmap_samples %>%
  replace(is.na(.), 0)  %>%
  t() %>%
  vegan::vegdist(na.rm = TRUE) %>%
  hclust() -> clust_t

pheatmap_samples %>%
  replace(is.na(.), 0)  %>%
  vegan::vegdist(na.rm = TRUE, method = "euc") %>%
  # cor() %>%
  hclust() -> clust

pheatmap_samples %>%
  pheatmap::pheatmap(fontsize_col = 6,
                     fontsize = 6,
                     cluster_rows = FALSE ,
                     scale = "row",
                     cluster_cols = clust_t,
                     na_col = "transparent",
                     annotation_row = tax_mapping %>% rownames_to_column('tmp_col') %>%
                       filter(PathoFact_AMR_ARG %in% rownames(pheatmap_samples)) %>%
                       distinct(PathoFact_AMR_ARG, .keep_all = TRUE) %>%
                       left_join(., 
                                 mean_ab,
                                 by = c("PathoFact_AMR_ARG" = "PathoFact_AMR_ARG")) %>% 
                       column_to_rownames("PathoFact_AMR_ARG") %>%
                       select(PathoFact_AMR_AMR_sub_class_multi, meanAbundance),
                     annotation_colors = list("PathoFact_AMR_AMR_sub_class_multi" = drug_classe_multi,
                                              "Treatment" = treat_col),
                     annotation_col = pheatmap_df %>%
                       select(Antibiotic_mg.mL,Day_of_Treatment, Treatment, Model, Sample) %>% 
                       distinct(Sample, .keep_all = TRUE) %>% 
                       column_to_rownames("Sample"))  -> p_pheat


p_pheat
```

<img src="chicken1_resistome_plots_draft_files/figure-html/unnamed-chunk-29-1.png" style="display: block; margin: auto;" />

```r
p_pheat %>% 
  ggsave(plot = .,bg = "transparent",scale = 2,
         units = "mm",
         dpi = 600,
         height = 0.618 * 84 * 1.5,
         width = 84  ,
         filename = "~/Desktop/Resistome_heatmap_VAN.eps")


p_pheat %>%
  export::graph2ppt(append = TRUE, width = 317.48031496 * 2, height = 0.618 * 317.48031496 * 3 , paper = "A4",  scaling = 1,
                    file = out_pptx)
```

```
## Exported graph as ~/Desktop/plots_manuscript_chicken1.pptx
```




all:



```r
ps_AMR_all %>% 
  physeq_sel_tax_table(c("PathoFact_AMR_AMR_sub_class", "PathoFact_AMR_AMR_sub_class_multi", "PathoFact_AMR_ARG")) %>% 
  speedyseq::tax_glom(., taxrank = "PathoFact_AMR_ARG")  %>% 
  physeq_sel_tax_table(c("PathoFact_AMR_ARG")) %>% 
  prune_samples(samples, .) %>% 
  physeq_add_metadata(metadata_mapping,
                      sample_column = "sample_id") %>% 
  transform_sample_counts(function(x) x/sum(x) * 100) %>% 
  filter_taxa(function(x) sum(x > 0) > 0, TRUE) %>% 
  # subset_taxa(PathoFact_AMR_AMR_sub_class_multi %!in% c("Unclassified", "multidrug (unclassified)", "MLS (unclassified)")) %>% 
  # ps_mutate(Day_of_Treatment = as.numeric(Day_of_Treatment)) %>% 
  ps_arrange(-Day_of_Treatment) %>% 
  speedyseq::psmelt() %>%
  # left_join(tax_mapping,
  #           by = c("Best_Hit_ARO" = "Best_Hit_ARO"),
  #           suffix = c("_x", "")) %>%
  # mutate(Sample = fct_reorder(Sample, Model,Treatment,Fermentation,Antibiotic_mg.mL, Day_of_Treatment )) %>%
  mutate(Abundance = na_if(Abundance, 0)) %>%
  mutate(logAbundance = log10(Abundance + 1)) %>% 
  mutate(PathoFact_AMR_ARG = PathoFact_AMR_ARG %>% stringr::str_trunc(45,  side ="center")) -> pheatmap_df


pheatmap_df %>%
  # replace(is.na(.), 0)  %>%
  mutate(log10Abundance = log10(Abundance)) %>% 
  select(Sample, PathoFact_AMR_ARG, Abundance) %>%
  filter(Abundance > 0) %>% 
  # filter(Abundance > 1) %>% 
  group_by(Sample, PathoFact_AMR_ARG) %>%
  pivot_wider(names_from =  Sample, PathoFact_AMR_ARG,
              values_from = Abundance) %>%
  # mutate(PathoFact_AMR_ARG = PathoFact_AMR_ARG %>%  str_to_lower(.)) %>% 
  # %>% str_to_lower(., locale = "en")
  column_to_rownames("PathoFact_AMR_ARG") -> pheatmap_samples

pheatmap_df %>%
  # replace(is.na(.), 0)  %>%
  mutate(log10Abundance = log10(Abundance)) %>% 
  select(Sample, PathoFact_AMR_ARG, Abundance) %>%
  group_by(PathoFact_AMR_ARG) %>% 
  summarise(meanAbundance = mean(Abundance, na.rm = TRUE)) -> mean_ab


pheatmap_samples %>%
  replace(is.na(.), 0)  %>%
  t() %>%
  vegan::vegdist(na.rm = TRUE) %>%
  hclust() -> clust_t

pheatmap_samples %>% 
  replace(is.na(.), 0)  %>%
  vegan::vegdist(na.rm = TRUE, method = "robust.aitchison") %>%
  # cor() %>%
  hclust(method = "ward.D2") -> clust

pheatmap_samples %>%
  pheatmap::pheatmap(fontsize_col = 2,
                     fontsize = 2,
                     cluster_rows = FALSE ,
                     scale = "row",
                     cluster_cols = clust_t,
                     na_col = "transparent",
                     annotation_row = tax_mapping %>% rownames_to_column('tmp_col') %>%
                       filter(PathoFact_AMR_ARG %in% rownames(pheatmap_samples)) %>%
                       distinct(PathoFact_AMR_ARG, .keep_all = TRUE) %>%
                       left_join(., 
                                 mean_ab,
                                 by = c("PathoFact_AMR_ARG" = "PathoFact_AMR_ARG")) %>% 
                       column_to_rownames("PathoFact_AMR_ARG") %>%
                       select(PathoFact_AMR_AMR_sub_class_multi, meanAbundance),
                     annotation_colors = list("PathoFact_AMR_AMR_sub_class_multi" = drug_classe_multi,
                                              "Treatment" = treat_col),
                     annotation_col = pheatmap_df %>%
                       select(Antibiotic_mg.mL,Day_of_Treatment, Treatment, Model, Sample) %>% 
                       distinct(Sample, .keep_all = TRUE) %>% 
                       column_to_rownames("Sample"))  -> p_pheat


p_pheat
```

<img src="chicken1_resistome_plots_draft_files/figure-html/unnamed-chunk-30-1.png" style="display: block; margin: auto;" />

```r
p_pheat %>% 
  ggsave(plot = .,bg = "transparent",scale = 2,
         units = "mm",
         dpi = 600,
         height = 0.618 * 84 * 2.5,
         width = 84  * 0.5,
         filename = "~/Desktop/Resistome_heatmap_all_gn.eps")


p_pheat %>%
  export::graph2ppt(append = TRUE, width = 317.48031496 * 60, height = 0.618 * 317.48031496 * 260 , paper = "A4",  scaling = 1,
                    file = out_pptx)
```

```
## Exported graph as ~/Desktop/plots_manuscript_chicken1.pptx
```

# Beta-div resistome:

Generate a temporary phyloseq object first - containing PathoFact_AMR_ARG - which will be used for the 2 models:
*TODO*: filter based on BP's expertise? Unclassified... ?

## AMR gNs:


```r
ps_AMR_all %>% 
  physeq_sel_tax_table(c("PathoFact_AMR_AMR_sub_class", "PathoFact_AMR_AMR_sub_class_multi", "PathoFact_AMR_ARG"))  %>%
  speedyseq::tax_glom(., taxrank = "PathoFact_AMR_ARG")  %>%
  filter_taxa(function(x) sum(x > 0) > 0, TRUE) %>% 
  prune_samples(samples, .) %>% 
  physeq_add_metadata(metadata_mapping,
                      sample_column = "sample_id") -> ps_tmp

ps_tmp
```

```
## phyloseq-class experiment-level object
## otu_table()   OTU Table:          [ 446 taxa and 30 samples ]:
## sample_data() Sample Data:        [ 30 samples by 42 sample variables ]:
## tax_table()   Taxonomy Table:     [ 446 taxa by 3 taxonomic ranks ]:
## taxa are rows
```


### Generate Aitchinson PCA:


```r
ps_tmp %>%
  subset_samples(Model == "Chicken") %>% 
  filter_taxa(function(x) sum(x > 0) > 0, TRUE) %>%
  # microbiome::transform(transform = "compositional") %>% 
  phyloseq_plot_bdiv(dlist = NULL,
                     m = "CoDa",
                     seed = 123,
                     axis1 = 1,
                     axis2 = 2) -> coda_resistome


ps_tmp %>%
  subset_samples(Model == "Chicken") %>%
  filter_taxa(function(x) sum(x > 0) > 0, TRUE) %>%
  microbiome::transform(transform = "compositional") %>%
  speedyseq::psmelt() %>%
    pivot_wider(names_from =  Sample, PathoFact_AMR_ARG,
              values_from = Abundance) %>%
  # mutate(PathoFact_AMR_ARG = PathoFact_AMR_ARG %>%  str_to_lower(.)) %>%
  # %>% str_to_lower(., locale = "en")
  column_to_rownames("PathoFact_AMR_ARG") %>%
  t() %>%
  replace(is.na(.), 0)  %>%
  vegan::vegdist(na.rm = TRUE, method = "robust.aitchison") %>% 
  magrittr::divide_by(100) -> robust.aitchison
#   # cor() %>%
#   # hclust(method = "ward.D2") -> clust
# 
# 
plot(
  robust.aitchison ,
  as.matrix(coda_resistome$aidist)[colnames(robust.aitchison %>% as.matrix()), rownames(robust.aitchison %>% as.matrix())] %>% as.dist()
)
```

<img src="chicken1_resistome_plots_draft_files/figure-html/unnamed-chunk-32-1.png" style="display: block; margin: auto;" />



```r
list = NULL

list$robust.aitchison <- robust.aitchison

ps_tmp %>%
  subset_samples(Model == "Chicken") %>% 
  filter_taxa(function(x) sum(x > 0) > 0, TRUE) %>%
  # microbiome::transform(transform = "compositional") %>% 
  phyloseq_plot_bdiv(dlist = list,
                     m = "PCoA",
                     seed = 123,
                     axis1 = 1,
                     axis2 = 2) -> robust.aitchison.pcoa
```

```
## [1] "robust.aitchison"
```


### Initial PCA:


```r
coda_resistome$PCA$layers[[1]] = NULL

coda_resistome$PCA + geom_point(size= 3.5,
                                aes(colour = Treatment ,
                                    shape = Antibiotic_mg.L,
                                    alpha = NULL)) +
  geom_path(data = coda_resistome$PCA$data %>%
              arrange(Day_of_Treatment) ,
            aes(colour = Treatment, group = interaction(Model, Antibiotic, Reactor_Treatment_Dose)),
            arrow = arrow(
              angle = 30, length = unit(0.15, "inches"),
              ends = "last", type = "open"
            ), linetype = "longdash", size = 0.1) +
  theme_light() +
  scale_color_manual(name = "", values = treat_col,
                     na.value = "black") +
  scale_fill_manual(name = "", values = treat_col,
                    na.value = "black") +
  scale_shape_manual(name = "" ,values = c(15,16,18,19), na.value =  17) +
  theme(legend.position = "bottom") +
  facet_grid(. ~ Model) -> pca_treat


pca_treat
```

<img src="chicken1_resistome_plots_draft_files/figure-html/unnamed-chunk-34-1.png" style="display: block; margin: auto;" />


```r
robust.aitchison.pcoa$robust.aitchison$layers[[1]] = NULL
robust.aitchison.pcoa$robust.aitchison$layers[[1]] = NULL
robust.aitchison.pcoa$robust.aitchison$layers[[2]] = NULL

robust.aitchison.pcoa$robust.aitchison + geom_point(size = 3.5,
                                aes(colour = Treatment ,
                                    shape = Antibiotic_mg.L,
                                    alpha = NULL)) +
  geom_path(data = robust.aitchison.pcoa$robust.aitchison$data %>% 
              arrange(Day_of_Treatment) ,
            aes(colour = Treatment, group = interaction(Model, Antibiotic, Reactor_Treatment_Dose)),         
            arrow = arrow(
              angle = 30, length = unit(0.15, "inches"),
              ends = "last", type = "open"
            ), linetype = "longdash", size = 0.1) +
  theme_light() +
  scale_color_manual(name = "", values = treat_col,
                     na.value = "black") +
  scale_fill_manual(name = "", values = treat_col,
                    na.value = "black") +
  scale_shape_manual(name = "" ,values = c(15,16,18,19), na.value =  17) +
  theme(legend.position = "bottom") + 
  facet_grid(. ~ Model) -> pca_treat

pca_treat
```

<img src="chicken1_resistome_plots_draft_files/figure-html/unnamed-chunk-35-1.png" style="display: block; margin: auto;" />


```r
pca_treat%>% 
  ggsave(plot = .,bg = "transparent",scale = 2,
         units = "mm",
         dpi = 600,
         height = 0.618 * 84,
         width = 84,
         filename = "~/Desktop/Resistome_beta_chick.eps")

pca_treat%>%
  export::graph2ppt(append = TRUE, width = 317.48031496  , height = 0.618 * 317.48031496  , paper = "A4",  scaling = 2,
                    file = out_pptx)
```

```
## Exported graph as ~/Desktop/plots_manuscript_chicken1.pptx
```

We generate a empty plot for vectors: 


```r
pca_treat + 
  scale_fill_manual(values = c("transparent")) + 
  scale_color_manual(values = c(rep("transparent", 6))) + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) -> empty_plot_tmp
```

```
## Scale for 'fill' is already present. Adding another scale for 'fill', which
## will replace the existing scale.
```

```
## Scale for 'colour' is already present. Adding another scale for 'colour',
## which will replace the existing scale.
```

### Next, we add AM drug classes info as envfit vectors on the generated 'empty' plot:


```r
phyloseq_add_taxa_vector_fix(dist = list$robust.aitchison,
                             phyloseq = ps_tmp %>%  
                               physeq_sel_tax_table(c("PathoFact_AMR_AMR_sub_class_multi")) %>% 
                               subset_samples(Model == "Chicken") %>% 
                               filter_taxa(function(x) sum(x > 0) > 0, TRUE),
                             figure_ord = empty_plot_tmp,
                             taxrank_glom = "PathoFact_AMR_AMR_sub_class_multi") -> envfit_class
```

```
## Loading required package: vegan
```

```
## Loading required package: permute
```

```
## Loading required package: lattice
```

```
## This is vegan 2.6-2
```

```r
envfit_class$signenvfit %>% 
  mutate_if(is.numeric, ~round(., 4)) %>% 
  DT::datatable()
```


```r
envfit_class$plot
```

<img src="chicken1_resistome_plots_draft_files/figure-html/unnamed-chunk-39-1.png" style="display: block; margin: auto;" />

```r
# ppp$plot + 
#   ggrepel::geom_text_repel(cex=2,
#                            aes(label= Phase),
#                            segment.color = 'black',
#                            segment.size = 0.5,
#                            ppp$plot$data %>% filter(Phase %in% c("Stab")) %>% unique()
#   ) -> pbeta

envfit_class$plot %>%
  export::graph2ppt(append = TRUE, width = 317.48031496, height = 0.618 * 317.48031496, paper = "A4",  scaling = 2,
                    file = out_pptx)
```

```
## Exported graph as ~/Desktop/plots_manuscript_chicken1.pptx
```

```r
envfit_class$plot %>% 
  ggsave(plot = .,bg = "transparent",scale = 2,
         units = "mm",
         dpi = 600,
         height = 0.618 * 84,
         width = 84,
         filename = "~/Desktop/Resistome_beta_chick_class.eps")
```

### Next, we add MGE prediction genes abundance as envfit vectors on the generated 'empty' plot:


```r
ps_AMR_all %>% 
  physeq_sel_tax_table(c("PathoFact_AMR_MGE_prediction_clean_2"))  %>%
  speedyseq::tax_glom(., taxrank = "PathoFact_AMR_MGE_prediction_clean_2") %>% 
  # physeq_sel_tax_table(c("PathoFact_AMR_AMR_sub_class_multi")) %>% 
  prune_samples(samples, .) %>% 
  subset_samples(Model == "Chicken") %>% 
  filter_taxa(function(x) sum(x > 0) > 0, TRUE) %>% 
  physeq_add_metadata(metadata_mapping,
                      sample_column = "sample_id") -> ps_tmp_envfit

# custom PathoFact_AMR_MGE_prediction cleaning- should be done once at the begining...

ps_tmp_envfit %>% 
  tax_table() %>% 
  data.frame() %>% 
  mutate(PathoFact_AMR_MGE_prediction_clean_2 = case_when(grepl("ambiguous", PathoFact_AMR_MGE_prediction_clean_2) ~ "ambiguous", TRUE ~ PathoFact_AMR_MGE_prediction_clean_2)) %>% 
  mutate(PathoFact_AMR_MGE_prediction_clean_2 = case_when(PathoFact_AMR_MGE_prediction_clean_2 %in% c("MAG", "chromosome")  ~ "MAG - chromosome", TRUE ~ PathoFact_AMR_MGE_prediction_clean_2)) %>% 
  as.matrix() -> tax_table(ps_tmp_envfit)

ps_tmp_envfit %>% 
  phyloseq_add_taxa_vector_fix(dist = coda_resistome$aidist,
                               phyloseq = . ,
                               taxnames_rm = c("unclassified", "ambiguous"),
                               vector_color = "tomato",
                               figure_ord = empty_plot_tmp,
                               taxrank_glom = "PathoFact_AMR_MGE_prediction_clean_2") -> envfit_MGE



envfit_MGE$signenvfit %>% 
  mutate_if(is.numeric, ~round(., 4)) %>% 
  DT::datatable()
```

```{=html}
<div id="htmlwidget-f640cd029f37c0438e19" style="width:100%;height:auto;" class="datatables html-widget"></div>
<script type="application/json" data-for="htmlwidget-f640cd029f37c0438e19">{"x":{"filter":"none","vertical":false,"data":[["c1_megahit_367038_90781.91203","c1_megahit_148336_31385.31603","c1_megahit_107010_14.880"],[-0.6636,0.5313,0.497],[-0.5241,0.5118,0.4279],[0.001,0.001,0.001],[0.7149,0.5443,0.4302],["MAG - chromosome","plasmid","phage"],[0.0013,0.0013,0.0013]],"container":"<table class=\"display\">\n  <thead>\n    <tr>\n      <th> <\/th>\n      <th>Axis.1<\/th>\n      <th>Axis.2<\/th>\n      <th>pval<\/th>\n      <th>r<\/th>\n      <th>tax_rank_plot<\/th>\n      <th>pval.adj<\/th>\n    <\/tr>\n  <\/thead>\n<\/table>","options":{"columnDefs":[{"className":"dt-right","targets":[1,2,3,4,6]},{"orderable":false,"targets":0}],"order":[],"autoWidth":false,"orderClasses":false}},"evals":[],"jsHooks":[]}</script>
```


```r
envfit_MGE$plot 
```

<img src="chicken1_resistome_plots_draft_files/figure-html/unnamed-chunk-41-1.png" style="display: block; margin: auto;" />

```r
# ppp$plot + 
#   ggrepel::geom_text_repel(cex=2,
#                            aes(label= Phase),
#                            segment.color = 'black',
#                            segment.size = 0.5,
#                            ppp$plot$data %>% filter(Phase %in% c("Stab")) %>% unique()
#   ) -> pbeta

envfit_MGE$plot %>% 
  ggsave(plot = .,bg = "transparent",scale = 2,
         units = "mm",
         dpi = 600,
         height = 0.618 * 84,
         width = 84,
         filename = "~/Desktop/Resistome_beta_chick_MGE.eps")

envfit_MGE$plot %>%
  export::graph2ppt(append = TRUE, width = 317.48031496, height = 0.618 * 317.48031496, paper = "A4",  scaling = 2,
                    file = out_pptx)
```

```
## Exported graph as ~/Desktop/plots_manuscript_chicken1.pptx
```

### Finally, we add AMR genes abunbundance as envfit vectors - similarly as before:


```r
# phyloseq_add_taxa_vector_fix(dist = coda_resistome$aidist, 
#                              ggrepel_force = 30,
#                              phyloseq = ps_tmp %>% 
#                                physeq_sel_tax_table(c("PathoFact_AMR_ARG")) %>% 
#                                subset_samples(Model == "Chicken") %>% 
#                                filter_taxa(function(x) sum(x > 0) > 0, TRUE) %>% 
#                                filter_taxa(., function(x){sum(x > 0.0000001) > 3}, prune = TRUE), # only keeping dominant and occurent ones
#                              figure_ord = empty_plot_tmp, 
#                              top_r = 20,
#                              vector_color = "chartreuse",
#                              taxrank_glom = "PathoFact_AMR_ARG") -> pca_env_ARG
# 
# pca_env_ARG$signenvfit %>% 
#   mutate_if(is.numeric, ~round(., 4)) %>% 
#   DT::datatable()  
```


```r
# pca_env_ARG$plot
```


```r
# pca_env_ARG$vectors
# 
# pca_env_ARG$vectors %>% 
#   ggsave(plot = .,bg = "transparent",scale = 2,
#          units = "mm",
#          dpi = 600,
#          height = 0.618 * 84,
#          width = 84,
#          filename = "~/Desktop/Resistome_beta_chick_ARG.eps")
# 
# pca_env_ARG$vectors  %>% 
#   export::graph2ppt(append = TRUE, width = 317.48031496  , height = 0.618 * 317.48031496  , paper = "A4",  scaling = 2,
#                     file = out_pptx)
```

## Gn


```r
ps_AMR_all %>% 
  physeq_sel_tax_table(c("PathoFact_AMR_AMR_sub_class", "PathoFact_AMR_AMR_sub_class_multi", "PathoFact_AMR_ARG"))  %>%
  # speedyseq::tax_glom(., taxrank = "PathoFact_AMR_ARG")  %>%
  filter_taxa(function(x) sum(x > 0) > 0, TRUE) %>% 
  prune_samples(samples, .) %>% 
  physeq_add_metadata(metadata_mapping,
                      sample_column = "sample_id") -> ps_tmp

ps_tmp
```

```
## phyloseq-class experiment-level object
## otu_table()   OTU Table:          [ 2322 taxa and 30 samples ]:
## sample_data() Sample Data:        [ 30 samples by 42 sample variables ]:
## tax_table()   Taxonomy Table:     [ 2322 taxa by 3 taxonomic ranks ]:
## taxa are rows
```


### Generate Aitchinson PCA:


```r
ps_tmp %>%
  subset_samples(Model == "Chicken") %>% 
  filter_taxa(function(x) sum(x > 0) > 0, TRUE) %>%
  # microbiome::transform(transform = "compositional") %>% 
  phyloseq_plot_bdiv(dlist = NULL,
                     m = "CoDa",
                     seed = 123,
                     axis1 = 1,
                     axis2 = 2) -> coda_resistome


ps_tmp %>%
  subset_samples(Model == "Chicken") %>%
  filter_taxa(function(x) sum(x > 0) > 0, TRUE) %>%
  microbiome::transform(transform = "compositional") %>%
  speedyseq::psmelt() %>%
    pivot_wider(names_from =  Sample, OTU,
              values_from = Abundance) %>%
  # mutate(PathoFact_AMR_ARG = PathoFact_AMR_ARG %>%  str_to_lower(.)) %>%
  # %>% str_to_lower(., locale = "en")
  column_to_rownames("OTU") %>%
  t() %>%
  replace(is.na(.), 0)  %>%
  vegan::vegdist(na.rm = TRUE, method = "robust.aitchison") %>% 
  magrittr::divide_by(100) -> robust.aitchison
#   # cor() %>%
#   # hclust(method = "ward.D2") -> clust
# 
# 
plot(
  robust.aitchison ,
  as.matrix(coda_resistome$aidist)[colnames(robust.aitchison %>% as.matrix()), rownames(robust.aitchison %>% as.matrix())] %>% as.dist()
)
```

<img src="chicken1_resistome_plots_draft_files/figure-html/unnamed-chunk-46-1.png" style="display: block; margin: auto;" />



```r
list = NULL

list$robust.aitchison <- robust.aitchison

ps_tmp %>%
  subset_samples(Model == "Chicken") %>% 
  filter_taxa(function(x) sum(x > 0) > 0, TRUE) %>%
  # microbiome::transform(transform = "compositional") %>% 
  phyloseq_plot_bdiv(dlist = list,
                     m = "PCoA",
                     seed = 123,
                     axis1 = 1,
                     axis2 = 2) -> robust.aitchison.pcoa
```

```
## [1] "robust.aitchison"
```


### Initial PCA:


```r
coda_resistome$PCA$layers[[1]] = NULL

coda_resistome$PCA + geom_point(size= 3.5,
                                aes(colour = Treatment ,
                                    shape = Antibiotic_mg.L,
                                    alpha = NULL)) +
  geom_path(data = coda_resistome$PCA$data %>%
              arrange(Day_of_Treatment) ,
            aes(colour = Treatment, group = interaction(Model, Antibiotic, Reactor_Treatment_Dose)),
            arrow = arrow(
              angle = 30, length = unit(0.15, "inches"),
              ends = "last", type = "open"
            ), linetype = "longdash", size = 0.1) +
  theme_light() +
  scale_color_manual(name = "", values = treat_col,
                     na.value = "black") +
  scale_fill_manual(name = "", values = treat_col,
                    na.value = "black") +
  scale_shape_manual(name = "" ,values = c(15,16,18,19), na.value =  17) +
  theme(legend.position = "bottom") +
  facet_grid(. ~ Model) -> pca_treat


pca_treat
```

<img src="chicken1_resistome_plots_draft_files/figure-html/unnamed-chunk-48-1.png" style="display: block; margin: auto;" />


```r
robust.aitchison.pcoa$robust.aitchison$layers[[1]] = NULL
robust.aitchison.pcoa$robust.aitchison$layers[[1]] = NULL
robust.aitchison.pcoa$robust.aitchison$layers[[2]] = NULL

robust.aitchison.pcoa$robust.aitchison + geom_point(size = 3.5,
                                aes(colour = Treatment ,
                                    shape = Antibiotic_mg.L,
                                    alpha = NULL)) +
  geom_path(data = robust.aitchison.pcoa$robust.aitchison$data %>% 
              arrange(Day_of_Treatment) ,
            aes(colour = Treatment, group = interaction(Model, Antibiotic, Reactor_Treatment_Dose)),         
            arrow = arrow(
              angle = 30, length = unit(0.15, "inches"),
              ends = "last", type = "open"
            ), linetype = "longdash", size = 0.1) +
  theme_light() +
  scale_color_manual(name = "", values = treat_col,
                     na.value = "black") +
  scale_fill_manual(name = "", values = treat_col,
                    na.value = "black") +
  scale_shape_manual(name = "" ,values = c(15,16,18,19), na.value =  17) +
  theme(legend.position = "bottom") + 
  facet_grid(. ~ Model) -> pca_treat

pca_treat
```

<img src="chicken1_resistome_plots_draft_files/figure-html/unnamed-chunk-49-1.png" style="display: block; margin: auto;" />


```r
pca_treat%>% 
  ggsave(plot = .,bg = "transparent",scale = 2,
         units = "mm",
         dpi = 600,
         height = 0.618 * 84,
         width = 84,
         filename = "~/Desktop/Resistome_beta_chick_gn.eps")

pca_treat%>%
  export::graph2ppt(append = TRUE, width = 317.48031496  , height = 0.618 * 317.48031496  , paper = "A4",  scaling = 2,
                    file = out_pptx)
```

```
## Exported graph as ~/Desktop/plots_manuscript_chicken1.pptx
```

We generate a empty plot for vectors: 


```r
pca_treat + 
  scale_fill_manual(values = c("transparent")) + 
  scale_color_manual(values = c(rep("transparent", 6))) + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) -> empty_plot_tmp
```

```
## Scale for 'fill' is already present. Adding another scale for 'fill', which
## will replace the existing scale.
```

```
## Scale for 'colour' is already present. Adding another scale for 'colour',
## which will replace the existing scale.
```

### Next, we add AM drug classes info as envfit vectors on the generated 'empty' plot:


```r
phyloseq_add_taxa_vector_fix(dist = list$robust.aitchison,
                             phyloseq = ps_tmp %>%  
                               physeq_sel_tax_table(c("PathoFact_AMR_AMR_sub_class_multi")) %>% 
                               subset_samples(Model == "Chicken") %>% 
                               filter_taxa(function(x) sum(x > 0) > 0, TRUE),
                             figure_ord = empty_plot_tmp,
                             taxrank_glom = "PathoFact_AMR_AMR_sub_class_multi") -> envfit_class

envfit_class$signenvfit %>% 
  mutate_if(is.numeric, ~round(., 4)) %>% 
  DT::datatable()
```


```r
envfit_class$plot
```

<img src="chicken1_resistome_plots_draft_files/figure-html/unnamed-chunk-53-1.png" style="display: block; margin: auto;" />

```r
# ppp$plot + 
#   ggrepel::geom_text_repel(cex=2,
#                            aes(label= Phase),
#                            segment.color = 'black',
#                            segment.size = 0.5,
#                            ppp$plot$data %>% filter(Phase %in% c("Stab")) %>% unique()
#   ) -> pbeta

envfit_class$plot %>%
  export::graph2ppt(append = TRUE, width = 317.48031496, height = 0.618 * 317.48031496, paper = "A4",  scaling = 2,
                    file = out_pptx)
```

```
## Exported graph as ~/Desktop/plots_manuscript_chicken1.pptx
```

```r
envfit_class$plot %>% 
  ggsave(plot = .,bg = "transparent",scale = 2,
         units = "mm",
         dpi = 600,
         height = 0.618 * 84,
         width = 84,
         filename = "~/Desktop/Resistome_beta_chick_class_gn.eps")
```

### Next, we add MGE prediction genes abundance as envfit vectors on the generated 'empty' plot:


```r
ps_AMR_all %>% 
  physeq_sel_tax_table(c("PathoFact_AMR_MGE_prediction_clean_2"))  %>%
  speedyseq::tax_glom(., taxrank = "PathoFact_AMR_MGE_prediction_clean_2") %>% 
  # physeq_sel_tax_table(c("PathoFact_AMR_AMR_sub_class_multi")) %>% 
  prune_samples(samples, .) %>% 
  subset_samples(Model == "Chicken") %>% 
  filter_taxa(function(x) sum(x > 0) > 0, TRUE) %>% 
  physeq_add_metadata(metadata_mapping,
                      sample_column = "sample_id") -> ps_tmp_envfit

# custom PathoFact_AMR_MGE_prediction cleaning- should be done once at the begining...

ps_tmp_envfit %>% 
  tax_table() %>% 
  data.frame() %>% 
  mutate(PathoFact_AMR_MGE_prediction_clean_2 = case_when(grepl("ambiguous", PathoFact_AMR_MGE_prediction_clean_2) ~ "ambiguous", TRUE ~ PathoFact_AMR_MGE_prediction_clean_2)) %>% 
  mutate(PathoFact_AMR_MGE_prediction_clean_2 = case_when(PathoFact_AMR_MGE_prediction_clean_2 %in% c("MAG", "chromosome")  ~ "MAG - chromosome", TRUE ~ PathoFact_AMR_MGE_prediction_clean_2)) %>% 
  as.matrix() -> tax_table(ps_tmp_envfit)

ps_tmp_envfit %>% 
  phyloseq_add_taxa_vector_fix(dist = coda_resistome$aidist,
                               phyloseq = . ,
                               taxnames_rm = c("unclassified", "ambiguous"),
                               vector_color = "tomato",
                               figure_ord = empty_plot_tmp,
                               taxrank_glom = "PathoFact_AMR_MGE_prediction_clean_2") -> envfit_MGE



envfit_MGE$signenvfit %>% 
  mutate_if(is.numeric, ~round(., 4)) %>% 
  DT::datatable()
```

```{=html}
<div id="htmlwidget-ef00aba85768a5b9bc73" style="width:100%;height:auto;" class="datatables html-widget"></div>
<script type="application/json" data-for="htmlwidget-ef00aba85768a5b9bc73">{"x":{"filter":"none","vertical":false,"data":[["c1_megahit_367038_90781.91203","c1_megahit_148336_31385.31603","c1_megahit_107010_14.880"],[0.7001,-0.5868,-0.5377],[0.3955,-0.3794,-0.3141],[0.001,0.001,0.004],[0.6466,0.4883,0.3878],["MAG - chromosome","plasmid","phage"],[0.0017,0.0017,0.005]],"container":"<table class=\"display\">\n  <thead>\n    <tr>\n      <th> <\/th>\n      <th>Axis.1<\/th>\n      <th>Axis.2<\/th>\n      <th>pval<\/th>\n      <th>r<\/th>\n      <th>tax_rank_plot<\/th>\n      <th>pval.adj<\/th>\n    <\/tr>\n  <\/thead>\n<\/table>","options":{"columnDefs":[{"className":"dt-right","targets":[1,2,3,4,6]},{"orderable":false,"targets":0}],"order":[],"autoWidth":false,"orderClasses":false}},"evals":[],"jsHooks":[]}</script>
```


```r
envfit_MGE$plot 
```

<img src="chicken1_resistome_plots_draft_files/figure-html/unnamed-chunk-55-1.png" style="display: block; margin: auto;" />

```r
# ppp$plot + 
#   ggrepel::geom_text_repel(cex=2,
#                            aes(label= Phase),
#                            segment.color = 'black',
#                            segment.size = 0.5,
#                            ppp$plot$data %>% filter(Phase %in% c("Stab")) %>% unique()
#   ) -> pbeta

envfit_MGE$plot %>% 
  ggsave(plot = .,bg = "transparent",scale = 2,
         units = "mm",
         dpi = 600,
         height = 0.618 * 84,
         width = 84,
         filename = "~/Desktop/Resistome_beta_chick_MGE_gn.eps")

envfit_MGE$plot %>%
  export::graph2ppt(append = TRUE, width = 317.48031496, height = 0.618 * 317.48031496, paper = "A4",  scaling = 2,
                    file = out_pptx)
```

```
## Exported graph as ~/Desktop/plots_manuscript_chicken1.pptx
```

### Finally, we add AMR genes abunbundance as envfit vectors - similarly as before:


```r
# phyloseq_add_taxa_vector_fix(dist = coda_resistome$aidist, 
#                              ggrepel_force = 30,
#                              phyloseq = ps_tmp %>% 
#                                physeq_sel_tax_table(c("PathoFact_AMR_ARG")) %>% 
#                                subset_samples(Model == "Chicken") %>% 
#                                filter_taxa(function(x) sum(x > 0) > 0, TRUE) %>% 
#                                filter_taxa(., function(x){sum(x > 0.0000001) > 3}, prune = TRUE), # only keeping dominant and occurent ones
#                              figure_ord = empty_plot_tmp, 
#                              top_r = 20,
#                              vector_color = "chartreuse",
#                              taxrank_glom = "PathoFact_AMR_ARG") -> pca_env_ARG
# 
# pca_env_ARG$signenvfit %>% 
#   mutate_if(is.numeric, ~round(., 4)) %>% 
#   DT::datatable()  
```


```r
# pca_env_ARG$plot
```


```r
# pca_env_ARG$vectors
# 
# pca_env_ARG$vectors %>% 
#   ggsave(plot = .,bg = "transparent",scale = 2,
#          units = "mm",
#          dpi = 600,
#          height = 0.618 * 84,
#          width = 84,
#          filename = "~/Desktop/Resistome_beta_chick_ARG.eps")
# 
# pca_env_ARG$vectors  %>% 
#   export::graph2ppt(append = TRUE, width = 317.48031496  , height = 0.618 * 317.48031496  , paper = "A4",  scaling = 2,
#                     file = out_pptx)
```


# Alpha-diversity resistome:

##  gN:


```r
ps_AMR_all %>% 
  physeq_sel_tax_table(c("PathoFact_AMR_AMR_sub_class", "PathoFact_AMR_AMR_sub_class_multi", "PathoFact_AMR_ARG"))  %>%
  # speedyseq::tax_glom(., taxrank = "PathoFact_AMR_ARG")  %>%
  filter_taxa(function(x) sum(x > 0) > 0, TRUE) %>% 
  prune_samples(samples, .) %>% 
  physeq_add_metadata(metadata_mapping,
                      sample_column = "sample_id") -> ps_tmp

ps_tmp
```

```
## phyloseq-class experiment-level object
## otu_table()   OTU Table:          [ 2322 taxa and 30 samples ]:
## sample_data() Sample Data:        [ 30 samples by 42 sample variables ]:
## tax_table()   Taxonomy Table:     [ 2322 taxa by 3 taxonomic ranks ]:
## taxa are rows
```



```r
sample_data(ps_tmp)$Read_nrom = 1000000

ps_tmp %>% 
  phyloseq_density_normalize(physeq = .,
                             # sample_idx = "Reactor",
                             value_idx = "Read_nrom") %>% 
  # transform_sample_counts(function(x) (x*sum(x)) *100) %>% 
  phyloseq_alphas(phylo = FALSE) -> resistome_alpha
```

```
## Loading required package: metagMisc
```

```
## 
## Attaching package: 'metagMisc'
```

```
## The following object is masked from 'package:purrr':
## 
##     some
```

```
## Loading required package: microbiome
```

```
## 
## microbiome R package (microbiome.github.com)
##     
## 
## 
##  Copyright (C) 2011-2022 Leo Lahti, 
##     Sudarshan Shetty et al. <microbiome.github.io>
```

```
## 
## Attaching package: 'microbiome'
```

```
## The following object is masked from 'package:vegan':
## 
##     diversity
```

```
## The following object is masked from 'package:scales':
## 
##     alpha
```

```
## The following object is masked from 'package:ggplot2':
## 
##     alpha
```

```
## The following object is masked from 'package:base':
## 
##     transform
```

```
## Observed richness
```

```
## Other forms of richness
```

```
## Diversity
```

```
## Evenness
```

```
## Dominance
```

```
## Rarity
```

```
## New names:
## • `id` -> `id...1`
## • `id` -> `id...45`
## • `id` -> `id...68`
```

```r
resistome_alpha %>% 
  glimpse()
```

```
## Rows: 30
## Columns: 73
## $ id...1                      <chr> "C1", "C6", "D10", "D11", "D12", "D13", "D…
## $ sample_id                   <chr> "C1", "C6", "D10", "D11", "D12", "D13", "D…
## $ Description2                <chr> NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA…
## $ Experiment                  <chr> "Continuous", "Continuous", "Continuous", …
## $ Reactor                     <chr> "TR5", "TR5", "TR3", "TR4", "TR5", "TR6", …
## $ Treatment                   <chr> "VAN+CCUG59168", "VAN+CCUG59168", "HV292.1…
## $ Day_of_Connection           <dbl> 13, 20, 17, 17, 17, 17, 19, 19, 19, 19, 19…
## $ Day_of_Treatment            <dbl> -3, 4, 1, 1, 1, 1, 3, 3, 3, 3, 3, 5, 5, 5,…
## $ Day_from_Inoculum           <dbl> 36, 43, 40, 40, 40, 40, 42, 42, 42, 42, 42…
## $ Enrichment                  <chr> "NotEnriched", "NotEnriched", "NotEnriched…
## $ Phase                       <chr> "Stab", "Treat", "Treat", "Treat", "Treat"…
## $ Treatment2                  <chr> "AB+E. faecium", "AB+E. faecium", "E. coli…
## $ Date                        <chr> "2020-06-24T00:00:00Z", "2020-07-01T00:00:…
## $ Paul                        <lgl> NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA…
## $ Reactor_Treatment           <chr> "TR5_VAN+CCUG59168", "TR5_VAN+CCUG59168", …
## $ GeneCopyNumberperML         <dbl> 3.94e+10, 2.50e+10, 5.24e+10, 3.74e+10, 2.…
## $ HV292.1_Copy_Number_permL   <dbl> NA, NA, NA, NA, NA, NA, 794000, NA, 887000…
## $ CCUG59168_Copy_Number_permL <dbl> NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, 0.…
## $ CTX_Copy_Number_permL       <dbl> NA, NA, NA, NA, NA, NA, 1507417.000, NA, 1…
## $ VAN_Copy_Number_permL       <dbl> NA, 28473350000, NA, NA, 41564875000, 0, N…
## $ Model                       <chr> "Chicken", "Chicken", "Chicken", "Chicken"…
## $ Antibiotic_mg.mL            <dbl> 90, 90, NA, 90, 90, NA, 20, 20, NA, 90, NA…
## $ Fermentation                <dbl> NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA…
## $ Antibiotic                  <chr> "VAN", "VAN", "CTX", "VAN", "VAN", "VAN", …
## $ Lactose_mM                  <dbl> 0.000, 0.550, 0.000, 0.000, 0.000, 0.000, …
## $ Glucose_mM                  <dbl> 0.000, 0.437, 0.000, 1.401, 1.584, 0.000, …
## $ Galactose_mM                <dbl> 0.000, 0.000, 0.000, 1.045, 1.135, 0.000, …
## $ Succinat_mM                 <dbl> 1.525, 7.267, 0.751, 4.792, 6.021, 0.726, …
## $ Lactat_mM                   <dbl> 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, …
## $ Formiat_mM                  <dbl> 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, …
## $ Acetat_mM                   <dbl> 92.696, 44.220, 86.398, 79.243, 77.564, 87…
## $ Propionat_mM                <dbl> 11.958, 21.102, 11.689, 17.005, 16.918, 11…
## $ Isobutyrat_mM               <dbl> 9.009, 0.000, 8.204, 8.207, 7.590, 8.115, …
## $ Butyrat_mM                  <dbl> 46.755, 3.984, 35.082, 28.479, 29.690, 38.…
## $ Isovalerat_mM               <dbl> 6.637, 5.818, 6.429, 6.182, 5.858, 6.024, …
## $ Valerat_mM                  <dbl> 7.112, 0.000, 6.960, 4.307, 3.135, 6.447, …
## $ Total_SCFA_mM               <dbl> 175.692, 82.391, 155.513, 148.215, 146.776…
## $ raw_metagenomic_pairs       <dbl> NA, NA, 46820988, 55159700, 45384782, 5208…
## $ Period                      <chr> "pret", "t1", "t1", "t1", "t1", "t1", "t1"…
## $ Reactor_Treatment_Dose      <chr> "TR5_VAN+CCUG5916890", "TR5_VAN+CCUG591689…
## $ Treatment_Dose              <chr> "VAN+CCUG5916890", "VAN+CCUG5916890", "HV2…
## $ Day_of_Treatment_num        <dbl> -3, 4, 1, 1, 1, 1, 3, 3, 3, 3, 3, 5, 5, 5,…
## $ Antibiotic_mg.L             <fct> NA, 90, NA, 90, 90, NA, 20, 20, NA, 90, NA…
## $ Read_nrom                   <dbl> 1e+06, 1e+06, 1e+06, 1e+06, 1e+06, 1e+06, …
## $ id...45                     <chr> "C1", "C6", "D10", "D11", "D12", "D13", "D…
## $ observed                    <dbl> 1325, 474, 1261, 1179, 1183, 1334, 1101, 9…
## $ chao1                       <dbl> 1325.0000, 477.9286, 1261.0000, 1181.5610,…
## $ diversity_inverse_simpson   <dbl> 187.43610, 80.19198, 229.15411, 99.36983, …
## $ diversity_gini_simpson      <dbl> 0.9946648, 0.9875299, 0.9956361, 0.9899366…
## $ diversity_shannon           <dbl> 5.855378, 4.759824, 5.881614, 5.136870, 5.…
## $ diversity_fisher            <dbl> 150.54432, 47.62789, 142.36969, 131.98324,…
## $ diversity_coverage          <dbl> 63, 28, 88, 32, 51, 80, 31, 34, 91, 25, 70…
## $ evenness_camargo            <dbl> 0.3874239, 0.8434157, 0.3852392, 0.4002576…
## $ evenness_pielou             <dbl> 0.8144723, 0.7725472, 0.8237947, 0.7263241…
## $ evenness_simpson            <dbl> 0.14146121, 0.16918140, 0.18172412, 0.0842…
## $ evenness_evar               <dbl> 0.17342497, 0.07843821, 0.14837578, 0.1206…
## $ evenness_bulla              <dbl> 0.3653134, 0.3509595, 0.3468210, 0.2421125…
## $ dominance_dbp               <dbl> 0.02474820, 0.02323370, 0.02002904, 0.0207…
## $ dominance_dmn               <dbl> 0.04091433, 0.04500841, 0.03195306, 0.0390…
## $ dominance_absolute          <dbl> 24748, 23234, 20029, 20701, 15061, 23553, …
## $ dominance_relative          <dbl> 0.02474820, 0.02323370, 0.02002904, 0.0207…
## $ dominance_simpson           <dbl> 0.005335152, 0.012470074, 0.004363875, 0.0…
## $ dominance_core_abundance    <dbl> 0.6009738, 0.8209173, 0.6203842, 0.8179441…
## $ dominance_gini              <dbl> 0.8786363, 0.9612000, 0.8829345, 0.9399624…
## $ rarity_log_modulo_skewness  <dbl> 2.043266, 2.001606, 2.048985, 2.052815, 2.…
## $ rarity_low_abundance        <dbl> 0.38554708, 0.08903784, 0.26037552, 0.1396…
## $ rarity_rare_abundance       <dbl> 0.31536252, 0.10907858, 0.26547153, 0.1056…
## $ id...68                     <chr> "C1", "C6", "D10", "D11", "D12", "D13", "D…
## $ Observed                    <dbl> 1325, 474, 1261, 1179, 1183, 1334, 1101, 9…
## $ Chao1                       <dbl> 1325.0000, 477.9286, 1261.0000, 1181.5610,…
## $ se.chao1                    <dbl> 0.04165094, 3.33748175, 0.03332011, 2.1431…
## $ ACE                         <dbl> 1325.1559, 477.5355, 1261.1713, 1182.8912,…
## $ se.ACE                      <dbl> 10.355195, 9.191495, 11.317160, 15.034705,…
```


```r
plot_alpha_div_NRP72 <- function(df = resistome_alpha, x = "Day_of_Treatment", y = "value",  point_size = 2.5, color = "Treatment", fill  = "Treatment", shape = "Antibiotic", alpha = "Antibiotic_mg.mL", facet_formula = paste0("alphadiversiy ~  Model "), ylab = "Resistome alpha-diversity", xlab = "Days (Treatment)", measures = c("Observed", "diversity_shannon"), path_group = "interaction(Model, Reactor_Treatment_Dose)"){
  
  df %>% 
    pivot_longer(cols = all_of(measures), values_to = "value", names_to = 'alphadiversiy', values_drop_na  = TRUE) %>%
    mutate(alphadiversiy = fct_relevel(alphadiversiy, measures)) -> df_ready
  
  df_ready %>% 
    ggplot(aes_string(x = x, y = y, color = color, fill = fill, shape = shape, alpha = alpha)) + #shape = Fermentation
    geom_point(size = point_size) + 
    geom_line(linetype = 2,  size = 0.5,  aes_string(group = path_group)) +
    labs(y = ylab, x = xlab) +
    # scale_color_manual(name = "", values = treat_col,
    #                    na.value = "black") +
    # scale_fill_manual(name = "", values = treat_col,
    #                   na.value = "black") +
    facet_grid(as.formula(facet_formula), scales = "free", space = "fixed", drop = TRUE) +
    # scale_shape_manual(name = "" ,values = antibio_shape, na.value =  17) +
    # scale_alpha_discrete(name = "", range=c(0.6, 1), na.value =  0.6) + 
    geom_vline(xintercept = 0, size = 0.5, linetype = 1, color = "grey60") -> p_alpha_both
  
  return(p_alpha_both)
}
```



```r
resistome_alpha %>% 
  plot_alpha_div_NRP72( shape = "as.factor(Antibiotic_mg.mL)", alpha = NULL) -> p_alpha 

p_alpha +  scale_color_manual(name = "", values = treat_col,
                              na.value = "black") +
  scale_fill_manual(name = "", values = treat_col,
                    na.value = "black") +
  scale_shape_manual(name = "" ,values = c(15,16,18,19), na.value =  17) + 
  theme_light() + facet_null() + facet_grid(Antibiotic + alphadiversiy  ~ Model, scales = "free", drop = TRUE) -> p_alpha #+
# ggh4x::facetted_pos_scales(
#   y = list(
#     alphadiversiy == "Observed" ~ scale_y_continuous(limits = c(120,270)),
#     alphadiversiy == "evenness_pielou" ~ scale_y_continuous(limits = c(0.6, 1))
#   )
# ) -> p_alpha

p_alpha
```

<img src="chicken1_resistome_plots_draft_files/figure-html/unnamed-chunk-62-1.png" style="display: block; margin: auto;" />

```r
p_alpha %>% 
  ggsave(plot = .,bg = "transparent",scale = 2,
         units = "mm",
         dpi = 600,
         height = 0.618 * 84,
         width = 84,
         filename = "~/Desktop/alpha_gn.eps")

p_alpha %>%
  export::graph2ppt(append = TRUE, width = 317.48031496  , height = 0.618 * 317.48031496  , paper = "A4",  scaling = 2,
                    file = out_pptx)
```

```
## Exported graph as ~/Desktop/plots_manuscript_chicken1.pptx
```

## AMR gNs:


```r
ps_AMR_all %>% 
  physeq_sel_tax_table(c("PathoFact_AMR_AMR_sub_class", "PathoFact_AMR_AMR_sub_class_multi", "PathoFact_AMR_ARG"))  %>%
  speedyseq::tax_glom(., taxrank = "PathoFact_AMR_ARG")  %>%
  filter_taxa(function(x) sum(x > 0) > 0, TRUE) %>% 
  prune_samples(samples, .) %>% 
  physeq_add_metadata(metadata_mapping,
                      sample_column = "sample_id") -> ps_tmp

ps_tmp
```

```
## phyloseq-class experiment-level object
## otu_table()   OTU Table:          [ 446 taxa and 30 samples ]:
## sample_data() Sample Data:        [ 30 samples by 42 sample variables ]:
## tax_table()   Taxonomy Table:     [ 446 taxa by 3 taxonomic ranks ]:
## taxa are rows
```



```r
sample_data(ps_tmp)$Read_nrom = 1000000

ps_tmp %>% 
  phyloseq_density_normalize(physeq = .,
                             # sample_idx = "Reactor",
                             value_idx = "Read_nrom") %>% 
  # transform_sample_counts(function(x) (x*sum(x)) *100) %>% 
  phyloseq_alphas(phylo = FALSE) -> resistome_alpha
```

```
## Observed richness
```

```
## Other forms of richness
```

```
## Diversity
```

```
## Evenness
```

```
## Dominance
```

```
## Rarity
```

```
## New names:
## • `id` -> `id...1`
## • `id` -> `id...45`
## • `id` -> `id...68`
```


```r
resistome_alpha %>% 
  plot_alpha_div_NRP72( shape = "as.factor(Antibiotic_mg.mL)", alpha = NULL) -> p_alpha 

p_alpha +  scale_color_manual(name = "", values = treat_col,
                              na.value = "black") +
  scale_fill_manual(name = "", values = treat_col,
                    na.value = "black") +
  scale_shape_manual(name = "" ,values = c(15,16,18,19), na.value =  17) + 
  theme_light() + facet_null() + facet_grid(Antibiotic + alphadiversiy  ~ Model, scales = "free", drop = TRUE) -> p_alpha #+
# ggh4x::facetted_pos_scales(
#   y = list(
#     alphadiversiy == "Observed" ~ scale_y_continuous(limits = c(120,270)),
#     alphadiversiy == "evenness_pielou" ~ scale_y_continuous(limits = c(0.6, 1))
#   )
# ) -> p_alpha

p_alpha
```

<img src="chicken1_resistome_plots_draft_files/figure-html/unnamed-chunk-65-1.png" style="display: block; margin: auto;" />

```r
p_alpha %>% 
  ggsave(plot = .,bg = "transparent",scale = 2,
         units = "mm",
         dpi = 600,
         height = 0.618 * 84,
         width = 84,
         filename = "~/Desktop/alpha.eps")

p_alpha %>%
  export::graph2ppt(append = TRUE, width = 317.48031496  , height = 0.618 * 317.48031496  , paper = "A4",  scaling = 2,
                    file = out_pptx)
```

```
## Exported graph as ~/Desktop/plots_manuscript_chicken1.pptx
```


# MGE :


```r
ps_AMR_all %>% 
  physeq_sel_tax_table(c("PathoFact_AMR_AMR_sub_class", "PathoFact_AMR_AMR_sub_class_multi", "PathoFact_AMR_ARG", "PathoFact_AMR_MGE_prediction_clean_2"))  %>%
  speedyseq::tax_glom(., taxrank = "PathoFact_AMR_MGE_prediction_clean_2")  %>% 
  physeq_sel_tax_table(c("PathoFact_AMR_MGE_prediction_clean_2")) %>% 
  # physeq_sel_tax_table(c("PathoFact_AMR_AMR_sub_class_multi")) %>% 
  prune_samples(samples, .) %>% 
  filter_taxa(function(x) sum(x > 0) > 0, TRUE) %>% 
  physeq_add_metadata(metadata_mapping,
                      sample_column = "sample_id") -> ps_tmp

ps_tmp %>% 
  tax_table() %>% 
  data.frame() %>% 
  mutate(PathoFact_AMR_MGE_prediction_clean_2 = case_when(grepl("ambiguous", PathoFact_AMR_MGE_prediction_clean_2) ~ "ambiguous", TRUE ~ PathoFact_AMR_MGE_prediction_clean_2)) %>% 
  mutate(PathoFact_AMR_MGE_prediction_clean_2 = case_when(PathoFact_AMR_MGE_prediction_clean_2 %in% c("MAG", "chromosome")  ~ "MAG - chromosome", TRUE ~ PathoFact_AMR_MGE_prediction_clean_2)) %>% 
  mutate(PathoFact_AMR_MGE_prediction_clean_2 = fct_relevel(PathoFact_AMR_MGE_prediction_clean_2, c("MAG", "chromosome", "plasmid", "phage", "ambiguous", "unclassified"))) %>% 
  as.matrix() -> tax_table(ps_tmp)
```

```
## Warning: Unknown levels in `f`: MAG, chromosome
```

```r
ps_tmp %>% 
  ps_arrange(-Day_of_Treatment) %>% 
  comp_barplot(
    tax_level = "PathoFact_AMR_MGE_prediction_clean_2",
    label = "Day_of_Treatment", # name an alternative variable to label axes
    n_taxa = 20, # give more taxa unique colours
    merge_other = FALSE, # split the "other" category to display alpha diversity
    bar_width = 0.9, # reduce the bar width to 70% of one row
    bar_outline_colour = "grey5",
    sample_order = "default") +#,# is the default (use NA to remove outlines)
  # facet_by = c("Model", "Antibiotic", "Antibiotic_mg.mL", "Treatment2")) +
  # group_by = c("Model", "Antibiotic", "Antibiotic_mg.mL", "Treatment2"))
  # coord_flip() + 
  theme_light() +
  scale_y_continuous(labels=percent) + 
  facet_wrap(Model ~ Antibiotic + Antibiotic_mg.mL + Treatment2, scales = "free", ncol = 6, nrow = 4) -> p2
```

```
## Warning in ps_counts(data, warn = TRUE): otu_table of counts is NOT available!
## Available otu_table contains non-zero values that are less than 1
```

```
## Warning in ps_counts(data = data): otu_table of counts is NOT available!
## Available otu_table contains non-zero values that are less than 1
```

```r
# 
# # Plot them side by side with the patchwork package.
# patch <- patchwork::wrap_plots(plot_list, nrow = 1, guides = "collect")
# patch & coord_flip() # make all plots horizontal (note: use & instead of +)
#   

p2 %>%
  ggpubr::set_palette("npg") -> p2
```

```
## Scale for 'fill' is already present. Adding another scale for 'fill', which
## will replace the existing scale.
```

```r
p2
```

<img src="chicken1_resistome_plots_draft_files/figure-html/unnamed-chunk-66-1.png" style="display: block; margin: auto;" />

```r
p2 %>% 
  ggsave(plot = .,bg = "transparent",scale = 2,
         units = "mm",
         dpi = 600,
         height = 0.618 * 84 * 2,
         width = 84 * 2,
         filename = "~/Desktop/AMR_carrier.eps")

p2 %>%
  export::graph2ppt(append = TRUE, width = 317.48031496 * 1.5  , height = 0.618 * 317.48031496  , paper = "A4",  scaling = 2,
                    file = out_pptx)
```

```
## Exported graph as ~/Desktop/plots_manuscript_chicken1.pptx
```

using the data stored in p2 object:


```r
p2$data %>% 
  mutate(top = fct_relevel(top, c("MAG", "chromosome", "plasmid", "phage", "ambiguous", "unclassified"))) %>% 
  mutate(Antibiotic_mg.mL = factor(Antibiotic_mg.mL, levels = c(NA, 20, 200,90, 600))) %>%
  ggplot(aes_string(x = "Day_of_Treatment", y = "Abundance", color = "Treatment", fill  = "Treatment", shape = "Antibiotic_mg.mL")) + #shape = Fermentation
  geom_point(size = 2.5) + 
  geom_line(linetype = 2,  size = 0.5,  aes_string(group = "interaction(Model, Reactor_Treatment_Dose)")) +
  labs(y = paste0("sum Rnum_Gi ") , x = "Day_of_Treatment") +
  scale_color_manual(name = "", values = treat_col,
                     na.value = "black") +
  scale_fill_manual(name = "", values = treat_col,
                    na.value = "black") +
  scale_shape_manual(name = "" ,values = c(15,16,18,19), na.value =  17) +
  # scale_alpha_continuous(name = "", range=c(0.6, 1), na.value =  1) +
  geom_vline(xintercept = 0, size = 0.5, linetype = 1, color = "grey60") +
  theme_light() + 
  facet_grid(Antibiotic + top ~ Model , scales = "free", drop = TRUE) +  scale_y_continuous(labels=percent) -> p1
```

```
## Warning: Unknown levels in `f`: MAG, chromosome
```

```r
p1
```

<img src="chicken1_resistome_plots_draft_files/figure-html/unnamed-chunk-67-1.png" style="display: block; margin: auto;" />

```r
# p1 %>%
#   ggpubr::set_palette("npg") -> p1

p1 %>% 
  ggsave(plot = .,bg = "transparent",scale = 2,
         units = "mm",
         dpi = 600,
         height = 0.618 * 84 * 3,
         width = 84 * 2,
         filename = "~/Desktop/AMRcarrier_2.eps")

p1 %>%
  export::graph2ppt(append = TRUE, width = 317.48031496 * 2  , height = 0.618 * 317.48031496 * 3 , paper = "A4",  scaling = 2,
                    file = out_pptx)
```

```
## Exported graph as ~/Desktop/plots_manuscript_chicken1.pptx
```

```r
# ps_AMR_all %>% 
#   physeq_sel_tax_table(c("PathoFact_AMR_AMR_sub_class", "PathoFact_AMR_AMR_sub_class_multi", "PathoFact_AMR_ARG", "PathoFact_AMR_MGE_prediction_clean_2"))  %>%
#   speedyseq::tax_glom(., taxrank = "PathoFact_AMR_MGE_prediction_clean_2")  %>% 
#   physeq_sel_tax_table(c("PathoFact_AMR_MGE_prediction_clean_2")) %>% 
#   # physeq_sel_tax_table(c("PathoFact_AMR_AMR_sub_class_multi")) %>% 
#   prune_samples(samples, .) %>% 
#   filter_taxa(function(x) sum(x > 0) > 0, TRUE) %>% 
#   physeq_add_metadata(metadata_mapping,
#                       sample_column = "sample_id") -> ps_tmp
# 
# ps_tmp %>% 
#   tax_table() %>% 
#   data.frame() %>% 
#   mutate(PathoFact_AMR_MGE_prediction_clean_2 = case_when(grepl("ambiguous", PathoFact_AMR_MGE_prediction_clean_2) ~ "ambiguous", TRUE ~ PathoFact_AMR_MGE_prediction_clean_2)) %>% 
#   # mutate(PathoFact_AMR_MGE_prediction_clean_2 = case_when(PathoFact_AMR_MGE_prediction_clean_2 %in% c("MAG", "chromosome")  ~ "MAG - chromosome", TRUE ~ PathoFact_AMR_MGE_prediction_clean_2)) %>% 
#   mutate(PathoFact_AMR_MGE_prediction_clean_2 = fct_relevel(PathoFact_AMR_MGE_prediction_clean_2, c("MAG", "chromosome", "plasmid", "phage", "ambiguous", "unclassified"))) %>% 
#   as.matrix() -> tax_table(ps_tmp)
# 
# ps_tmp %>% 
#   ps_arrange(-Day_of_Treatment) %>% 
#   comp_barplot(
#     tax_level = "PathoFact_AMR_MGE_prediction_clean_2",
#     label = "Day_of_Treatment", # name an alternative variable to label axes
#     n_taxa = 20, # give more taxa unique colours
#     merge_other = FALSE, # split the "other" category to display alpha diversity
#     bar_width = 0.9, # reduce the bar width to 70% of one row
#     bar_outline_colour = "grey5",
#     sample_order = "default") +#,# is the default (use NA to remove outlines)
#   # facet_by = c("Model", "Antibiotic", "Antibiotic_mg.mL", "Treatment2")) +
#   # group_by = c("Model", "Antibiotic", "Antibiotic_mg.mL", "Treatment2"))
#   coord_flip() +  theme_light() +
#   scale_y_continuous(labels=percent) + 
#   facet_wrap(Model ~ Antibiotic + Antibiotic_mg.mL + Treatment2, scales = "free", ncol = 6, nrow = 4) -> p2
# 
# p_carrier_all <- p2
# # # Plot them side by side with the patchwork package.
# # patch <- patchwork::wrap_plots(plot_list, nrow = 1, guides = "collect")
# # patch & coord_flip() # make all plots horizontal (note: use & instead of +)
# #   
# 
# p2 %>%
#   ggpubr::set_palette("npg") -> p2
# 
# p2
# 
# p2 %>% 
#   ggsave(plot = .,bg = "transparent",scale = 2,
#          units = "mm",
#          dpi = 600,
#          height = 0.618 * 84 * 2,
#          width = 84 * 2,
#          filename = "~/Desktop/AMR_carrier_3.eps")
# 
# p2 %>%
#   export::graph2ppt(append = TRUE, width = 317.48031496 * 2  , height = 0.618 * 317.48031496 * 0.75 , paper = "A4",  scaling = 2,
#                     file = out_pptx)
```

### Fig. 2b


```r
ps_AMR_all %>% 
  physeq_sel_tax_table(c("PathoFact_AMR_AMR_sub_class_multi",
                         "PathoFact_AMR_ARG"))  %>%
  # physeq_sel_tax_table(c("PathoFact_AMR_ARG"))  %>% 
  speedyseq::tax_glom(., taxrank = "PathoFact_AMR_ARG")  %>% 
  prune_samples(samples, .) %>% 
  physeq_add_metadata(metadata_mapping,
                      sample_column = "sample_id") %>% 
  filter_taxa(function(x) sum(x > 0) > 0, TRUE) %>% 
  speedyseq::psmelt()  %>% 
  filter(! grepl('nclassified', PathoFact_AMR_AMR_sub_class_multi)) -> df

results <- vector("list", length(df$PathoFact_AMR_AMR_sub_class_multi %>%  unique()))
names(results) = c(df$PathoFact_AMR_AMR_sub_class_multi %>%  unique())

for (dd in df$PathoFact_AMR_AMR_sub_class_multi %>%  unique())
{
  print(dd)
  df %>% 
    filter(PathoFact_AMR_AMR_sub_class_multi == !!dd) %>% 
    select(-Lactose_mM:-Total_SCFA_mM) %>% 
    group_by(Sample, PathoFact_AMR_AMR_sub_class_multi) %>% 
    summarise(sum_abundance = sum(Abundance), .groups = "drop") %>% 
    mutate(sum_abundance = na_if(sum_abundance, 0)) %>%
    left_join(metadata_mapping,
              by = c("Sample" = "sample_id")) %>% 
    # filter(Model == "Chicken") %>% 
    # mutate(Day_of_Treatment_num = as.double(Day_of_Treatment)) %>% 
    # mutate(Day_of_Treatment = fct_relevel(Day_of_Treatment, c(NA, "-3", "-1", "1" , "3", "4", "5","7", "30", "47", "48")))  %>% 
    mutate(Antibiotic_mg.mL = factor(Antibiotic_mg.mL, levels = c(NA, 20, 200,90, 600))) %>%
    ggplot(aes_string(x = "Day_of_Treatment", y = "sum_abundance", color = "Treatment", fill  = "Treatment", shape = "Antibiotic_mg.mL")) + #shape = Fermentation
    geom_point(size = 2.5) + 
    geom_line(linetype = 2,  size = 0.5,  aes_string(group = "interaction(Model, Reactor_Treatment_Dose)")) +
    labs(y = paste0("sum Rnum_Gi ") , x = "Day_of_Treatment") +
    scale_color_manual(name = "", values = treat_col,
                       na.value = "black") +
    scale_fill_manual(name = "", values = treat_col,
                      na.value = "black") +
    scale_shape_manual(name = "" ,values = c(15,16,18,19), na.value =  17) +
    # scale_alpha_continuous(name = "", range=c(0.6, 1), na.value =  1) +
    geom_vline(xintercept = 0, size = 0.5, linetype = 1, color = "grey60") +
    ggtitle(dd) +
    theme_light() + 
    facet_grid(Antibiotic  ~ Model , scales = "free", drop = TRUE) +  scale_y_continuous(labels=scientific) -> p1
  # scale_y_continuous(labels = function(x) format(x, scientific = TRUE)) -> p1
  
  # p1 + facet_grid(as.formula(paste0("Antibiotic_mg.mL + PathoFact_AMR_AMR_sub_class_multi ~  Model  ")), scales = "fixed", space = "fixed", drop = TRUE)
  # 
  # p1 + facet_wrap(as.formula(paste0("PathoFact_AMR_AMR_sub_class_multi ~  Model  ")),scales = "free"
  #)
  results[[dd]] <- p1
  
  
  results[[dd]]  %>% 
    ggsave(plot = .,bg = "transparent",scale = 2,
           units = "mm",
           dpi = 600,
           height = 0.618 * 84,
           width = 84,
           filename = paste0("~/Desktop/AMR_class_", dd, ".eps"))
  
  results[[dd]]  %>%
    export::graph2ppt(append = TRUE, width = 317.48031496 , height = 0.618 * 317.48031496, paper = "A4",  scaling = 2,
                      file = out_pptx)
  
}
```

```
## [1] "peptide antibiotic"
```

```
## Exported graph as ~/Desktop/plots_manuscript_chicken1.pptx
```

```
## [1] "glycopeptide antibiotic"
```

```
## Exported graph as ~/Desktop/plots_manuscript_chicken1.pptx
```

```
## [1] "lincosamide antibiotic"
```

```
## Exported graph as ~/Desktop/plots_manuscript_chicken1.pptx
```

```
## [1] "tetracycline antibiotic"
```

```
## Exported graph as ~/Desktop/plots_manuscript_chicken1.pptx
```

```
## [1] "multidrug"
```

```
## Exported graph as ~/Desktop/plots_manuscript_chicken1.pptx
```

```
## [1] "macrolide antibiotic"
```

```
## Exported graph as ~/Desktop/plots_manuscript_chicken1.pptx
```

```
## [1] "aminoglycoside antibiotic"
```

```
## Exported graph as ~/Desktop/plots_manuscript_chicken1.pptx
```

```
## [1] "fluoroquinolone antibiotic"
```

```
## Exported graph as ~/Desktop/plots_manuscript_chicken1.pptx
```

```
## [1] "penam"
```

```
## Exported graph as ~/Desktop/plots_manuscript_chicken1.pptx
```

```
## [1] "fosfomycin"
```

```
## Exported graph as ~/Desktop/plots_manuscript_chicken1.pptx
```

```
## [1] "phenicol antibiotic"
```

```
## Exported graph as ~/Desktop/plots_manuscript_chicken1.pptx
```

```
## [1] "bicyclomycin"
```

```
## Warning: Removed 1 rows containing missing values (geom_point).
```

```
## Warning: Removed 1 row(s) containing missing values (geom_path).
```

```
## Warning: Removed 1 rows containing missing values (geom_point).
```

```
## Warning: Removed 1 row(s) containing missing values (geom_path).
```

```
## Exported graph as ~/Desktop/plots_manuscript_chicken1.pptx
```

```
## [1] "aminocoumarin antibiotic"
```

```
## Exported graph as ~/Desktop/plots_manuscript_chicken1.pptx
```

```
## [1] "streptogramin antibiotic"
```

```
## Exported graph as ~/Desktop/plots_manuscript_chicken1.pptx
```

```
## [1] "nitroimidazole antibiotic"
```

```
## Warning: Removed 1 rows containing missing values (geom_point).
```

```
## Warning: Removed 1 rows containing missing values (geom_point).
```

```
## Exported graph as ~/Desktop/plots_manuscript_chicken1.pptx
```

```
## [1] "diaminopyrimidine antibiotic"
```

```
## Exported graph as ~/Desktop/plots_manuscript_chicken1.pptx
```

```
## [1] "cephamycin"
```

```
## Warning: Removed 18 rows containing missing values (geom_point).
```

```
## Warning: Removed 18 row(s) containing missing values (geom_path).
```

```
## Warning: Removed 18 rows containing missing values (geom_point).
```

```
## Warning: Removed 18 row(s) containing missing values (geom_path).
```

```
## Exported graph as ~/Desktop/plots_manuscript_chicken1.pptx
```

```
## [1] "rifamycin antibiotic"
```

```
## Exported graph as ~/Desktop/plots_manuscript_chicken1.pptx
```

```
## [1] "beta-lactam antibiotics"
```

```
## Warning: Removed 3 rows containing missing values (geom_point).
```

```
## Warning: Removed 3 row(s) containing missing values (geom_path).
```

```
## Warning: Removed 3 rows containing missing values (geom_point).
```

```
## Warning: Removed 3 row(s) containing missing values (geom_path).
```

```
## Exported graph as ~/Desktop/plots_manuscript_chicken1.pptx
```

```
## [1] "nucleoside antibiotic"
```

```
## Warning: Removed 2 rows containing missing values (geom_point).
```

```
## Warning: Removed 2 rows containing missing values (geom_point).
```

```
## Exported graph as ~/Desktop/plots_manuscript_chicken1.pptx
```

```
## [1] "mupirocin"
```

```
## Exported graph as ~/Desktop/plots_manuscript_chicken1.pptx
```

```
## [1] "carbapenem"
```

```
## Warning: Removed 1 rows containing missing values (geom_point).
```

```
## Warning: Removed 1 row(s) containing missing values (geom_path).
```

```
## Warning: Removed 1 rows containing missing values (geom_point).
```

```
## Warning: Removed 1 row(s) containing missing values (geom_path).
```

```
## Exported graph as ~/Desktop/plots_manuscript_chicken1.pptx
```

```
## [1] "cephalosporin"
```

```
## Exported graph as ~/Desktop/plots_manuscript_chicken1.pptx
```

```
## [1] "triclosan"
```

```
## Warning: Removed 2 rows containing missing values (geom_point).
```

```
## Warning: Removed 2 row(s) containing missing values (geom_path).
```

```
## Warning: Removed 2 rows containing missing values (geom_point).
```

```
## Warning: Removed 2 row(s) containing missing values (geom_path).
```

```
## Exported graph as ~/Desktop/plots_manuscript_chicken1.pptx
```

```
## [1] "acridine dye"
```

```
## Warning: Removed 2 rows containing missing values (geom_point).
## Removed 2 row(s) containing missing values (geom_path).
```

```
## Warning: Removed 2 rows containing missing values (geom_point).
```

```
## Warning: Removed 2 row(s) containing missing values (geom_path).
```

```
## Exported graph as ~/Desktop/plots_manuscript_chicken1.pptx
```

```
## [1] "macrolide"
```

```
## Warning: Removed 19 rows containing missing values (geom_point).
```

```
## Warning: Removed 18 row(s) containing missing values (geom_path).
```

```
## Warning: Removed 19 rows containing missing values (geom_point).
```

```
## Warning: Removed 18 row(s) containing missing values (geom_path).
```

```
## Exported graph as ~/Desktop/plots_manuscript_chicken1.pptx
```

```
## [1] "fusidic acid"
```

```
## Warning: Removed 11 rows containing missing values (geom_point).
```

```
## Warning: Removed 9 row(s) containing missing values (geom_path).
```

```
## Warning: Removed 11 rows containing missing values (geom_point).
```

```
## Warning: Removed 9 row(s) containing missing values (geom_path).
```

```
## Exported graph as ~/Desktop/plots_manuscript_chicken1.pptx
```

```
## [1] "para-aminosalicylic acid"
```

```
## Warning: Removed 29 rows containing missing values (geom_point).
```

```
## Warning: Removed 29 row(s) containing missing values (geom_path).
```

```
## geom_path: Each group consists of only one observation. Do you need to adjust
## the group aesthetic?
```

```
## Warning: Removed 29 rows containing missing values (geom_point).
## Removed 29 row(s) containing missing values (geom_path).
```

```
## geom_path: Each group consists of only one observation. Do you need to adjust
## the group aesthetic?
## Exported graph as ~/Desktop/plots_manuscript_chicken1.pptx
```

```
## [1] "pleuromutilin antibiotic"
```

```
## Warning: Removed 28 rows containing missing values (geom_point).
```

```
## Warning: Removed 28 row(s) containing missing values (geom_path).
```

```
## geom_path: Each group consists of only one observation. Do you need to adjust
## the group aesthetic?
```

```
## Warning: Removed 28 rows containing missing values (geom_point).
## Removed 28 row(s) containing missing values (geom_path).
```

```
## geom_path: Each group consists of only one observation. Do you need to adjust
## the group aesthetic?
## Exported graph as ~/Desktop/plots_manuscript_chicken1.pptx
```

Vancomycin is a glycopeptide antibiotic used in the prophylaxis and treatment of infections caused by Gram-positive bacteria. Vancomycin inhibits the synthesis of peptidoglycan, the major component of the cell wall of gram-positive bacteria. Its mechanism of action is unusual in that it acts by binding precursors of peptidoglycan, rather than by interacting with an enzyme. -> glycopeptide antibiotic


```r
results$`glycopeptide antibiotic`
```

<img src="chicken1_resistome_plots_draft_files/figure-html/unnamed-chunk-70-1.png" style="display: block; margin: auto;" />

```r
# results$`glycopeptide antibiotic` %>% 
#   export::graph2ppt(append = TRUE, width = 8, height = 4,
#                     file = out_pptx)
```

Cefotaxime is a semisynthetic cephalosporin taken parenterally. It is resistant to most beta-lactamases and active against Gram-negative rods and cocci due to its aminothiazoyl and methoximino functional groups. -> cephalosporin


```r
results$cephalosporin 
```

<img src="chicken1_resistome_plots_draft_files/figure-html/unnamed-chunk-71-1.png" style="display: block; margin: auto;" />

```r
# results$cephalosporin %>% 
#   export::graph2ppt(append = TRUE, width = 8, height = 4,
#                     file = out_pptx)
```

same but normalized by baseline:


```r
# ps_AMR_all %>% 
#   physeq_sel_tax_table(c("PathoFact_AMR_AMR_sub_class_multi",
#                          "PathoFact_AMR_ARG"))  %>%
#   # physeq_sel_tax_table(c("PathoFact_AMR_ARG"))  %>% 
#   speedyseq::tax_glom(., taxrank = "PathoFact_AMR_ARG")  %>% 
#   prune_samples(samples, .) %>% 
#   physeq_add_metadata(metadata_mapping,
#                       sample_column = "sample_id") %>% 
#   filter_taxa(function(x) sum(x > 0) > 0, TRUE) %>% 
#   speedyseq::psmelt() %>% 
#   filter(! grepl('nclassified', PathoFact_AMR_AMR_sub_class_multi)) -> df
# 
# results <- vector("list", length(df$PathoFact_AMR_AMR_sub_class_multi %>%  unique()))
# names(results) = c(df$PathoFact_AMR_AMR_sub_class_multi %>%  unique())
# 
# for (dd in df$PathoFact_AMR_AMR_sub_class_multi %>%  unique())
# {
#   print(dd)
#   
#   df %>% 
#     filter(PathoFact_AMR_AMR_sub_class_multi == !!dd) %>% 
#     select(-Lactose_mM:-Total_SCFA_mM) %>% 
#     group_by(Sample, PathoFact_AMR_AMR_sub_class_multi) %>% 
#     summarise(sum_abundance = sum(Abundance), .groups = "drop") %>% 
#     mutate(sum_abundance = na_if(sum_abundance, 0)) %>%
#     left_join(metadata_mapping,
#               by = c("Sample" = "sample_id")) %>% 
#     ungroup() %>% 
#     # mutate(Reactor_Treatment_Dose = paste0(Period, Reactor_Treatment_Dose)) %>% 
#     group_by(Model, Fermentation, 
#              Reactor_Treatment_Dose) %>% 
#     arrange(Day_of_Treatment) %>% 
#     mutate(per_cent_change = sum_abundance / first(sum_abundance)) %>% 
#     # filter(Model == "Chicken") %>% 
#     # mutate(Day_of_Treatment_num = as.double(Day_of_Treatment)) %>% 
#     # mutate(Day_of_Treatment = fct_relevel(Day_of_Treatment, c(NA, "-3", "-1", "1" , "3", "4", "5","7", "30", "47", "48")))  %>% 
#     mutate(Antibiotic_mg.mL = factor(Antibiotic_mg.mL, levels = c(NA, 20, 200,90, 600))) %>%
#     ggplot(aes_string(x = "Day_of_Treatment", y = "per_cent_change", color = "Treatment", fill  = "Treatment", shape = "Antibiotic_mg.mL")) + #shape = Fermentation
#     geom_point(size = 2.5) + 
#     geom_line(linetype = 2,  size = 0.5,  aes_string(group = "interaction(Model, Reactor_Treatment_Dose)")) +
#     labs(y = paste0("sum Rnum_Gi / baseline ") , x = "Day_of_Treatment") +
#     scale_color_manual(name = "", values = treat_col,
#                        na.value = "black") +
#     scale_fill_manual(name = "", values = treat_col,
#                       na.value = "black") +
#     scale_shape_manual(name = "" ,values = c(15,16,18,19), na.value =  17) +
#     # scale_alpha_continuous(name = "", range=c(0.6, 1), na.value =  1) +
#     geom_vline(xintercept = 0, size = 0.5, linetype = 1, color = "grey60") +
#     ggtitle(dd) +
#     theme_light() + 
#     facet_grid(Antibiotic  ~ Model , scales = "free_x", drop = TRUE) -> p1
#   # scale_y_continuous(labels = function(x) format(x, scientific = TRUE)) -> p1
#   
#   # p1 + facet_grid(as.formula(paste0("Antibiotic_mg.mL + PathoFact_AMR_AMR_sub_class_multi ~  Model  ")), scales = "fixed", space = "fixed", drop = TRUE)
#   # 
#   # p1 + facet_wrap(as.formula(paste0("PathoFact_AMR_AMR_sub_class_multi ~  Model  ")),scales = "free"
#   #)
#   results[[dd]] <- p1
#   
#   
#   results[[dd]]  %>% 
#     ggsave(plot = .,bg = "transparent",scale = 2,
#            units = "mm",
#            dpi = 600,
#            height = 0.618 * 84,
#            width = 84,
#            filename = paste0("~/Desktop/AMR_class_", dd, "_baseline.eps"))
#   
#   results[[dd]]  %>%
#     export::graph2ppt(append = TRUE, width = 317.48031496 , height = 0.618 * 317.48031496, paper = "A4",  scaling = 2,
#                       file = out_pptx)
#   
# }
```



Vancomycin is a glycopeptide antibiotic used in the prophylaxis and treatment of infections caused by Gram-positive bacteria. Vancomycin inhibits the synthesis of peptidoglycan, the major component of the cell wall of gram-positive bacteria. Its mechanism of action is unusual in that it acts by binding precursors of peptidoglycan, rather than by interacting with an enzyme. -> glycopeptide antibiotic


```r
# results$`glycopeptide antibiotic`

# results$`glycopeptide antibiotic` %>% 
#   export::graph2ppt(append = TRUE, width = 8, height = 4,
#                     file = out_pptx)
```

Cefotaxime is a semisynthetic cephalosporin taken parenterally. It is resistant to most beta-lactamases and active against Gram-negative rods and cocci due to its aminothiazoyl and methoximino functional groups. -> cephalosporin


```r
# results$cephalosporin 

# results$cephalosporin %>% 
#   export::graph2ppt(append = TRUE, width = 8, height = 4,
#                     file = out_pptx)
```

- line plot MAG-Chrosmosome - Plasmid - phage - unknown facet | quant & % 

-- quant


```r
# ps_AMR_all %>% 
#   physeq_sel_tax_table(c("PathoFact_AMR_AMR_sub_class_multi",
#                          "PathoFact_AMR_MGE_prediction_clean_2",
#                          "PathoFact_AMR_ARG"))  %>%
#   # physeq_sel_tax_table(c("PathoFact_AMR_ARG"))  %>% 
#   speedyseq::tax_glom(., taxrank = "PathoFact_AMR_ARG")  %>% 
#   physeq_sel_tax_table(c("PathoFact_AMR_AMR_sub_class_multi",
#                          "PathoFact_AMR_ARG",
#                          "PathoFact_AMR_MGE_prediction_clean_2"))  %>%
#   speedyseq::tax_glom(., taxrank = "PathoFact_AMR_MGE_prediction_clean_2")  %>% 
#   prune_samples(samples, .) %>% 
#   physeq_add_metadata(metadata_mapping,
#                       sample_column = "sample_id") %>% 
#   filter_taxa(function(x) sum(x > 0) > 0, TRUE) %>% 
#   speedyseq::psmelt() %>% 
#   filter(!grepl('nclassified', PathoFact_AMR_AMR_sub_class_multi))  -> df
# 
# results <- vector("list", length(df$PathoFact_AMR_AMR_sub_class_multi %>%  unique()))
# names(results) = c(df$PathoFact_AMR_AMR_sub_class_multi %>%  unique())
# 
# for (dd in df %>% filter(PathoFact_AMR_AMR_sub_class_multi != "pleuromutilin antibiotic") %>% pull(PathoFact_AMR_AMR_sub_class_multi) %>%  unique()
# ){
#   
#   print(dd)
#   
#   df %>% 
#     filter(PathoFact_AMR_AMR_sub_class_multi == !!dd) %>% 
#     select(-Lactose_mM:-Total_SCFA_mM) %>% 
#     mutate(PathoFact_AMR_MGE_prediction_clean_2 = case_when(PathoFact_AMR_MGE_prediction_clean_2 %in% c("MAG", "chromosome")  ~ "MAG - chromosome", TRUE ~ PathoFact_AMR_MGE_prediction_clean_2)) %>% 
#     mutate(PathoFact_AMR_MGE_prediction_clean_2 = fct_relevel(PathoFact_AMR_MGE_prediction_clean_2, c("MAG - chromosome", "plasmid", "phage", "ambiguous", "unclassified"))) %>% 
#     group_by(Sample, PathoFact_AMR_MGE_prediction_clean_2, PathoFact_AMR_AMR_sub_class_multi) %>% 
#     summarise(sum_abundance = sum(Abundance), .groups = "drop") %>% 
#     mutate(sum_abundance = na_if(sum_abundance, 0)) %>%
#     left_join(metadata_mapping,
#               by = c("Sample" = "sample_id")) %>% 
#     filter(! PathoFact_AMR_MGE_prediction_clean_2  %in% c("unclassified", "ambiguous")) %>% 
#     # filter(Model == "Chicken") %>% 
#     # mutate(Day_of_Treatment_num = as.double(Day_of_Treatment)) %>% 
#     # mutate(Day_of_Treatment = fct_relevel(Day_of_Treatment, c(NA, "-3", "-1", "1" , "3", "4", "5","7", "30", "47", "48")))  %>% 
#     mutate(Antibiotic_mg.mL = factor(Antibiotic_mg.mL, levels = c(NA, 20, 200,90, 600))) %>%
#     ggplot(aes_string(x = "Day_of_Treatment", y = "sum_abundance", color = "Treatment", fill  = "Treatment", shape = "Antibiotic_mg.mL")) + #shape = Fermentation
#     geom_point(size = 2.5) + 
#     geom_line(linetype = 2,  size = 0.5,  aes_string(group = "interaction(Model, Reactor_Treatment_Dose)")) +
#     labs(y = paste0("sum Rnum_Gi ") , x = "Day_of_Treatment") +
#     scale_color_manual(name = "", values = treat_col,
#                        na.value = "black") +
#     scale_fill_manual(name = "", values = treat_col,
#                       na.value = "black") +
#     scale_shape_manual(name = "" ,values = c(15,16,18,19), na.value =  17) +
#     # scale_alpha_continuous(name = "", range=c(0.6, 1), na.value =  1) +
#     geom_vline(xintercept = 0, size = 0.5, linetype = 1, color = "grey60") +
#     ggtitle(dd) +
#     theme_light() + 
#     facet_grid(Antibiotic + PathoFact_AMR_MGE_prediction_clean_2 ~ Model, scales = "free", drop = TRUE) +  scale_y_continuous(labels=scientific) -> p1
#   
#   # p1 + facet_grid(as.formula(paste0("Antibiotic_mg.mL + PathoFact_AMR_AMR_sub_class_multi ~  Model  ")), scales = "fixed", space = "fixed", drop = TRUE)
#   # 
#   # p1 + facet_wrap(as.formula(paste0("PathoFact_AMR_AMR_sub_class_multi ~  Model  ")),scales = "free"
#   #)
#   results[[dd]] <- p1
#   
#   results[[dd]]  %>% 
#     ggsave(plot = .,bg = "transparent",scale = 2,
#            units = "mm",
#            dpi = 600,
#            height = 0.618 * 84 * 1.5,
#            width = 84,
#            filename = paste0("~/Desktop/lineplot_", dd, ".eps"))
#   
#   results[[dd]]  %>%
#     export::graph2ppt(append = TRUE, width = 317.48031496 , height = 0.618 * 317.48031496, paper = "A4",  scaling = 2,
#                       file = out_pptx)
#   
# }
```

-- then %


```r
# ps_AMR_all %>% 
#   physeq_sel_tax_table(c("PathoFact_AMR_AMR_sub_class_multi",
#                          "PathoFact_AMR_MGE_prediction_clean_2",
#                          "PathoFact_AMR_ARG"))  %>%
#   # physeq_sel_tax_table(c("PathoFact_AMR_ARG"))  %>% 
#   speedyseq::tax_glom(., taxrank = "PathoFact_AMR_ARG")  %>% 
#   physeq_sel_tax_table(c("PathoFact_AMR_AMR_sub_class_multi",
#                          "PathoFact_AMR_ARG",
#                          "PathoFact_AMR_MGE_prediction_clean_2"))  %>%
#   speedyseq::tax_glom(., taxrank = "PathoFact_AMR_MGE_prediction_clean_2")  %>% 
#   prune_samples(samples, .) %>% 
#   physeq_add_metadata(metadata_mapping,
#                       sample_column = "sample_id") %>% 
#   filter_taxa(function(x) sum(x > 0) > 0, TRUE) %>% 
#   speedyseq::psmelt() %>% 
#   select(-Lactose_mM:-Total_SCFA_mM) %>% 
#   mutate(PathoFact_AMR_MGE_prediction_clean_2 = case_when(PathoFact_AMR_MGE_prediction_clean_2 %in% c("MAG", "chromosome")  ~ "MAG - chromosome", TRUE ~ PathoFact_AMR_MGE_prediction_clean_2)) %>% 
#   mutate(PathoFact_AMR_MGE_prediction_clean_2 = fct_relevel(PathoFact_AMR_MGE_prediction_clean_2, c("MAG - chromosome", "plasmid", "phage", "ambiguous", "unclassified"))) -> df  #%>% 
# # filter(grepl('unclassified', PathoFact_AMR_AMR_sub_class_multi))  -> df
# 
# results <- vector("list", length(df$PathoFact_AMR_AMR_sub_class_multi %>%  unique()))
# names(results) = c(df$PathoFact_AMR_AMR_sub_class_multi %>%  unique())
# 
# for (dd in df %>% filter(PathoFact_AMR_AMR_sub_class_multi %!in% c("pleuromutilin antibiotic", "MLS (unclassified)", "Unclassified", "multidrug (unclassified)")) %>% pull(PathoFact_AMR_AMR_sub_class_multi) %>%  unique()
# )
# {
#   print(dd)
#   
#   df %>% 
#     filter(PathoFact_AMR_AMR_sub_class_multi == !!dd) %>% 
#     select(Sample, PathoFact_AMR_AMR_sub_class_multi, PathoFact_AMR_MGE_prediction_clean_2, Abundance) %>%
#     group_by(PathoFact_AMR_AMR_sub_class_multi, PathoFact_AMR_MGE_prediction_clean_2, Sample) -> df_tmp
#   
#   df_tmp %>% 
#     summarise(sum_MGE = sum(Abundance), .groups = "drop") %>%
#     mutate(sum_MGE = na_if(sum_MGE, 0)) -> df_MGE_sum
#   
#   df_tmp %>% 
#     group_by(Sample) %>% 
#     summarise(sample_sum = sum(Abundance), .groups = "drop") %>% 
#     mutate(sample_sum = na_if(sample_sum, 0)) -> df_sample_sum 
#   
#   df_MGE_sum %>% 
#     left_join(df_sample_sum,
#               by = c("Sample" = "Sample")) %>% 
#     mutate(percent_MGE = sum_MGE / sample_sum) %>% 
#     # left_join()
#     left_join(metadata_mapping,
#               by = c("Sample" = "sample_id")) %>% 
#     filter(! PathoFact_AMR_MGE_prediction_clean_2  %in% c("unclassified", "ambiguous")) %>%
#     # filter(Model == "Chicken") %>% 
#     # mutate(Day_of_Treatment_num = as.double(Day_of_Treatment)) %>% 
#     # mutate(Day_of_Treatment = fct_relevel(Day_of_Treatment, c(NA, "-3", "-1", "1" , "3", "4", "5","7", "30", "47", "48")))  %>% 
#     mutate(Antibiotic_mg.mL = factor(Antibiotic_mg.mL, levels = c(NA, 20, 200,90, 600))) %>%
#     ggplot(aes_string(x = "Day_of_Treatment", y = "percent_MGE", color = "Treatment", fill  = "Treatment", shape = "Antibiotic_mg.mL")) + #shape = Fermentation
#     geom_point(size = 2.5) + 
#     geom_line(linetype = 2,  size = 0.5,  aes_string(group = "interaction(Model, Reactor_Treatment_Dose)")) +
#     labs(y = paste0(" proportion % ") , x = "Day_of_Treatment") +
#     scale_color_manual(name = "", values = treat_col,
#                        na.value = "black") +
#     scale_fill_manual(name = "", values = treat_col,
#                       na.value = "black") +
#     scale_shape_manual(name = "" ,values = c(15,16,18,19), na.value =  17) +
#     # scale_alpha_continuous(name = "", range=c(0.6, 1), na.value =  1) +
#     geom_vline(xintercept = 0, size = 0.5, linetype = 1, color = "grey60") +
#     ggtitle(dd) +
#     theme_light() + 
#     facet_grid(Antibiotic + PathoFact_AMR_MGE_prediction_clean_2 ~ Model, scales = "free", drop = TRUE) +  scale_y_continuous(labels=percent) -> p1
#   # scale_y_continuous(labels = function(x) format(x, scientific = TRUE)) -> p1
#   
#   # p1 + facet_grid(as.formula(paste0("Antibiotic_mg.mL + PathoFact_AMR_AMR_sub_class_multi ~  Model  ")), scales = "fixed", space = "fixed", drop = TRUE)
#   # 
#   # p1 + facet_wrap(as.formula(paste0("PathoFact_AMR_AMR_sub_class_multi ~  Model  ")),scales = "free"
#   #)
#   results[[dd]] <- p1
#   
#   
#   results[[dd]]  %>% 
#     ggsave(plot = .,bg = "transparent",scale = 2,
#            units = "mm",
#            dpi = 600,
#            height = 0.618 * 84 * 1.5,
#            width = 84,
#            filename = paste0("~/Desktop/lineplot_", dd, "_percent.eps"))
#   
#   results[[dd]]  %>%
#     export::graph2ppt(append = TRUE, width = 317.48031496 , height = 0.618 * 317.48031496, paper = "A4",  scaling = 2,
#                       file = out_pptx)
#   
# }
```

and then baseline ratio.


```r
# ps_AMR_all %>% 
#   physeq_sel_tax_table(c("PathoFact_AMR_AMR_sub_class_multi",
#                          "PathoFact_AMR_MGE_prediction_clean_2",
#                          "PathoFact_AMR_ARG"))  %>%
#   # physeq_sel_tax_table(c("PathoFact_AMR_ARG"))  %>% 
#   speedyseq::tax_glom(., taxrank = "PathoFact_AMR_ARG")  %>% 
#   physeq_sel_tax_table(c("PathoFact_AMR_AMR_sub_class_multi",
#                          "PathoFact_AMR_ARG",
#                          "PathoFact_AMR_MGE_prediction_clean_2"))  %>%
#   speedyseq::tax_glom(., taxrank = "PathoFact_AMR_MGE_prediction_clean_2")  %>% 
#   prune_samples(samples, .) %>% 
#   physeq_add_metadata(metadata_mapping,
#                       sample_column = "sample_id") %>% 
#   filter_taxa(function(x) sum(x > 0) > 0, TRUE) %>% 
#   speedyseq::psmelt() %>% 
#   filter(! grepl('nclassified', PathoFact_AMR_AMR_sub_class_multi))  -> df
# 
# results <- vector("list", length(df$PathoFact_AMR_AMR_sub_class_multi %>%  unique()))
# names(results) = c(df$PathoFact_AMR_AMR_sub_class_multi %>%  unique())
# 
# for (dd in df %>% filter(PathoFact_AMR_AMR_sub_class_multi != "pleuromutilin antibiotic") %>% pull(PathoFact_AMR_AMR_sub_class_multi) %>%  unique()
# )
# {
#   print(dd)
#   
#   df %>% 
#     filter(PathoFact_AMR_AMR_sub_class_multi == !!dd) %>% 
#     select(-Lactose_mM:-Total_SCFA_mM) %>% 
#     mutate(PathoFact_AMR_MGE_prediction_clean_2 = case_when(PathoFact_AMR_MGE_prediction_clean_2 %in% c("MAG", "chromosome")  ~ "MAG - chromosome", TRUE ~ PathoFact_AMR_MGE_prediction_clean_2)) %>% 
#     mutate(PathoFact_AMR_MGE_prediction_clean_2 = fct_relevel(PathoFact_AMR_MGE_prediction_clean_2, c("MAG - chromosome", "plasmid", "phage", "ambiguous", "unclassified"))) %>% 
#     group_by(Sample, PathoFact_AMR_AMR_sub_class_multi, PathoFact_AMR_MGE_prediction_clean_2) %>% 
#     summarise(sum_abundance = sum(Abundance), .groups = "drop") %>% 
#     mutate(sum_abundance = na_if(sum_abundance, 0)) %>%
#     ungroup() %>% 
#     left_join(metadata_mapping,
#               by = c("Sample" = "sample_id")) %>% 
#     filter(! PathoFact_AMR_MGE_prediction_clean_2  %in% c("unclassified", "ambiguous")) %>%
#     group_by(PathoFact_AMR_AMR_sub_class_multi, PathoFact_AMR_MGE_prediction_clean_2, Model, Fermentation, Reactor_Treatment_Dose) %>% 
#     arrange(Day_of_Treatment) %>% 
#     mutate(per_cent_change = sum_abundance / first(sum_abundance)) %>% 
#     ggplot(aes_string(x = "Day_of_Treatment", y = "per_cent_change", color = "Treatment", fill  = "Treatment", shape = "Antibiotic_mg.L")) + #shape = Fermentation
#     geom_point(size = 2.5) + 
#     geom_line(linetype = 2,  size = 0.5,  aes_string(group = "interaction(Model, Reactor_Treatment_Dose)")) +
#     labs(y = "sum Rnum_Gi / baseline", x = "Day_of_Treatment") +
#     scale_color_manual(name = "", values = treat_col,
#                        na.value = "black") +
#     scale_fill_manual(name = "", values = treat_col,
#                       na.value = "black") +
#     scale_shape_manual(name = "" ,values = c(15,16,18,19), na.value =  17) +
#     # scale_alpha_continuous(name = "", range=c(0.6, 1), na.value =  1) +
#     geom_vline(xintercept = 0, size = 0.5, linetype = 1, color = "grey60") +
#     geom_hline(yintercept = 1, size = 0.5, linetype = 2, color = "grey60") +
#     ggtitle(dd) +
#     theme_light() + facet_grid(Antibiotic + PathoFact_AMR_MGE_prediction_clean_2 ~ Model, scales = "free", drop = TRUE) -> p1
#   
#   p1
#   
#   results[[dd]] <- p1
#   
#   
#   results[[dd]]  %>% 
#     ggsave(plot = .,bg = "transparent",scale = 2,
#            units = "mm",
#            dpi = 600,
#            height = 0.618 * 84 * 1.5,
#            width = 84,
#            filename = paste0("~/Desktop/lineplot_", dd, "_baseline_ratio.eps"))
#   
#   results[[dd]]  %>%
#     export::graph2ppt(append = TRUE, width = 317.48031496 , height = 0.618 * 317.48031496, paper = "A4",  scaling = 2,
#                       file = out_pptx)
#   
# }
```

- glycopeptide antibiotic:


```r
ps_AMR_all %>% 
  physeq_sel_tax_table(c("PathoFact_AMR_AMR_sub_class", "PathoFact_AMR_AMR_sub_class_multi", "PathoFact_AMR_ARG", "PathoFact_AMR_MGE_prediction_clean_2"))  %>%
  prune_samples(samples, .) %>% 
  physeq_add_metadata(metadata_mapping,
                      sample_column = "sample_id") %>% 
  subset_taxa(PathoFact_AMR_AMR_sub_class_multi %in% c("glycopeptide antibiotic")) %>% 
  filter_taxa(function(x) sum(x > 0) > 0, TRUE) -> ps_tmp

ps_tmp %>% 
  tax_table() %>% 
  data.frame() %>% 
  mutate(PathoFact_AMR_MGE_prediction_clean_2 = case_when(grepl("ambiguous", PathoFact_AMR_MGE_prediction_clean_2) ~ "ambiguous", TRUE ~ PathoFact_AMR_MGE_prediction_clean_2)) %>% 
  # mutate(PathoFact_AMR_MGE_prediction_clean_2 = case_when(PathoFact_AMR_MGE_prediction_clean_2 %in% c("MAG", "chromosome")  ~ "MAG - chromosome", TRUE ~ PathoFact_AMR_MGE_prediction_clean_2)) %>% 
  mutate(PathoFact_AMR_MGE_prediction_clean_2 = fct_relevel(PathoFact_AMR_MGE_prediction_clean_2, c("MAG", "chromosome", "plasmid", "phage", "ambiguous", "unclassified"))) %>% 
  as.matrix() -> tax_table(ps_tmp)

ps_tmp %>% 
  ps_arrange(-Day_of_Treatment) %>% 
  comp_barplot(
    tax_transform_for_plot = "identity",
    tax_level = "PathoFact_AMR_ARG",
    label = "Day_of_Treatment", # name an alternative variable to label axes
    n_taxa = 20, # give more taxa unique colours
    merge_other = FALSE, # split the "other" category to display alpha diversity
    bar_width = 0.9, # reduce the bar width to 70% of one row
    bar_outline_colour = "grey5",
    sample_order = "default") +#,# is the default (use NA to remove outlines)
  # facet_by = c("Model", "Antibiotic", "Antibiotic_mg.mL", "Treatment2")) +
  # group_by = c("Model", "Antibiotic", "Antibiotic_mg.mL", "Treatment2"))
  coord_flip() +  theme_light() +
  scale_y_continuous(labels=scientific) -> p2
```

```
## Short values detected in phyloseq tax_table (nchar<4) :
## Consider using tax_fix() to make taxa uniquely identifiable
```

```
## Warning in ps_counts(data, warn = TRUE): otu_table of counts is NOT available!
## Available otu_table contains non-zero values that are less than 1
```

```
## Warning in ps_counts(data = data): otu_table of counts is NOT available!
## Available otu_table contains non-zero values that are less than 1
```

```r
# 
# # Plot them side by side with the patchwork package.
# patch <- patchwork::wrap_plots(plot_list, nrow = 1, guides = "collect")
# patch & coord_flip() # make all plots horizontal (note: use & instead of +)
#   



p2 + facet_wrap(Model ~ Antibiotic + Antibiotic_mg.mL + Treatment2, scales = "free", ncol = 6, nrow = 4) -> p2_1

p2_1
```

<img src="chicken1_resistome_plots_draft_files/figure-html/unnamed-chunk-78-1.png" style="display: block; margin: auto;" />

```r
p2 + facet_wrap(Model ~ Antibiotic + Antibiotic_mg.mL + Treatment2, scales = "free_y", ncol = 6, nrow = 4) -> p2_2

p2_2
```

<img src="chicken1_resistome_plots_draft_files/figure-html/unnamed-chunk-78-2.png" style="display: block; margin: auto;" />

```r
p2_1 %>% 
  ggsave(plot = .,bg = "transparent",scale = 2,
         units = "mm",
         dpi = 600,
         height = 0.618 * 84 * 2,
         width = 84 * 2,
         filename = "~/Desktop/glyco_1-1.eps")

p2_2 %>% 
  ggsave(plot = .,bg = "transparent",scale = 2,
         units = "mm",
         dpi = 600,
         height = 0.618 * 84 * 2,
         width = 84 * 2,
         filename = "~/Desktop/glyco_1-2.eps")


p2_1 %>%
  export::graph2ppt(append = TRUE, width = 317.48031496 * 2  , height = 0.618 * 317.48031496 * 0.75 , paper = "A4",  scaling = 2,
                    file = out_pptx)
```

```
## Exported graph as ~/Desktop/plots_manuscript_chicken1.pptx
```

```r
p2_2 %>%
  export::graph2ppt(append = TRUE, width = 317.48031496 * 2  , height = 0.618 * 317.48031496 * 0.75 , paper = "A4",  scaling = 2,
                    file = out_pptx)
```

```
## Exported graph as ~/Desktop/plots_manuscript_chicken1.pptx
```

- cephalosporin:


```r
ps_AMR_all %>% 
  physeq_sel_tax_table(c("PathoFact_AMR_AMR_sub_class", "PathoFact_AMR_AMR_sub_class_multi", "PathoFact_AMR_ARG", "PathoFact_AMR_MGE_prediction_clean_2"))  %>%
  prune_samples(samples, .) %>% 
  physeq_add_metadata(metadata_mapping,
                      sample_column = "sample_id") %>% 
  subset_taxa(PathoFact_AMR_AMR_sub_class_multi %in% c("cephalosporin")) %>% 
  filter_taxa(function(x) sum(x > 0) > 0, TRUE) -> ps_tmp

ps_tmp %>% 
  tax_table() %>% 
  data.frame() %>% 
  mutate(PathoFact_AMR_MGE_prediction_clean_2 = case_when(grepl("ambiguous", PathoFact_AMR_MGE_prediction_clean_2) ~ "ambiguous", TRUE ~ PathoFact_AMR_MGE_prediction_clean_2)) %>% 
  # mutate(PathoFact_AMR_MGE_prediction_clean_2 = case_when(PathoFact_AMR_MGE_prediction_clean_2 %in% c("MAG", "chromosome")  ~ "MAG - chromosome", TRUE ~ PathoFact_AMR_MGE_prediction_clean_2)) %>% 
  mutate(PathoFact_AMR_MGE_prediction_clean_2 = fct_relevel(PathoFact_AMR_MGE_prediction_clean_2, c("MAG", "chromosome", "plasmid", "phage", "ambiguous", "unclassified"))) %>% 
  as.matrix() -> tax_table(ps_tmp)
```

```
## Warning: Unknown levels in `f`: phage, ambiguous, unclassified
```

```r
ps_tmp %>% 
  ps_arrange(-Day_of_Treatment) %>% 
  comp_barplot(
    tax_transform_for_plot = "identity",
    tax_level = "PathoFact_AMR_ARG",
    label = "Day_of_Treatment", # name an alternative variable to label axes
    n_taxa = 20, # give more taxa unique colours
    merge_other = FALSE, # split the "other" category to display alpha diversity
    bar_width = 0.9, # reduce the bar width to 70% of one row
    bar_outline_colour = "grey5",
    sample_order = "default") +#,# is the default (use NA to remove outlines)
  # facet_by = c("Model", "Antibiotic", "Antibiotic_mg.mL", "Treatment2")) +
  # group_by = c("Model", "Antibiotic", "Antibiotic_mg.mL", "Treatment2"))
  coord_flip() +  theme_light() +
  scale_y_continuous(labels=scientific) -> p2
```

```
## Short values detected in phyloseq tax_table (nchar<4) :
## Consider using tax_fix() to make taxa uniquely identifiable
```

```
## Warning in ps_counts(data, warn = TRUE): otu_table of counts is NOT available!
## Available otu_table contains non-zero values that are less than 1
```

```
## Warning in ps_counts(data = data): otu_table of counts is NOT available!
## Available otu_table contains non-zero values that are less than 1
```

```r
# 
# # Plot them side by side with the patchwork package.
# patch <- patchwork::wrap_plots(plot_list, nrow = 1, guides = "collect")
# patch & coord_flip() # make all plots horizontal (note: use & instead of +)
#   



p2 + facet_wrap(Model ~ Antibiotic + Antibiotic_mg.mL + Treatment2, scales = "free", ncol = 6, nrow = 4) -> p2_1

p2_1
```

<img src="chicken1_resistome_plots_draft_files/figure-html/unnamed-chunk-79-1.png" style="display: block; margin: auto;" />

```r
p2 + facet_wrap(Model ~ Antibiotic + Antibiotic_mg.mL + Treatment2, scales = "free_y", ncol = 6, nrow = 4) -> p2_2

p2_2
```

<img src="chicken1_resistome_plots_draft_files/figure-html/unnamed-chunk-79-2.png" style="display: block; margin: auto;" />

```r
p2_1 %>% 
  ggsave(plot = .,bg = "transparent",scale = 2,
         units = "mm",
         dpi = 600,
         height = 0.618 * 84 * 2,
         width = 84 * 2,
         filename = "~/Desktop/cephalosporin_1-1.eps")

p2_2 %>% 
  ggsave(plot = .,bg = "transparent",scale = 2,
         units = "mm",
         dpi = 600,
         height = 0.618 * 84 * 2,
         width = 84 * 2,
         filename = "~/Desktop/cephalosporin_1-2.eps")


p2_1 %>%
  export::graph2ppt(append = TRUE, width = 317.48031496 * 2  , height = 0.618 * 317.48031496 * 0.75 , paper = "A4",  scaling = 2,
                    file = out_pptx)
```

```
## Exported graph as ~/Desktop/plots_manuscript_chicken1.pptx
```

```r
p2_2 %>%
  export::graph2ppt(append = TRUE, width = 317.48031496 * 2  , height = 0.618 * 317.48031496 * 0.75 , paper = "A4",  scaling = 2,
                    file = out_pptx)
```

```
## Exported graph as ~/Desktop/plots_manuscript_chicken1.pptx
```


- VANA


```r
ps_AMR_all %>% 
  physeq_sel_tax_table(c("PathoFact_AMR_AMR_sub_class", "PathoFact_AMR_AMR_sub_class_multi", "PathoFact_AMR_MGE_prediction_clean_2", "PathoFact_AMR_ARG"))  %>%
  prune_samples(samples, .) %>% 
  physeq_add_metadata(metadata_mapping,
                      sample_column = "sample_id") %>% 
  subset_taxa(PathoFact_AMR_ARG %in% c("VANA")) %>% 
  filter_taxa(function(x) sum(x > 0) > 0, TRUE) -> ps_tmp

ps_tmp %>% 
  tax_table() %>% 
  data.frame() %>% 
  rownames_to_column('GN_id') %>% 
  mutate(PathoFact_AMR_MGE_prediction_clean_2 = case_when(grepl("ambiguous", PathoFact_AMR_MGE_prediction_clean_2) ~ "ambiguous", TRUE ~ PathoFact_AMR_MGE_prediction_clean_2)) %>% 
  # mutate(PathoFact_AMR_MGE_prediction_clean_2 = case_when(PathoFact_AMR_MGE_prediction_clean_2 %in% c("MAG", "chromosome")  ~ "MAG - chromosome", TRUE ~ PathoFact_AMR_MGE_prediction_clean_2)) %>% 
  mutate(PathoFact_AMR_MGE_prediction_clean_2 = fct_relevel(PathoFact_AMR_MGE_prediction_clean_2, c("MAG", "chromosome", "plasmid", "phage", "ambiguous", "unclassified"))) %>% 
  mutate(GN = paste(GN_id, PathoFact_AMR_MGE_prediction_clean_2, PathoFact_AMR_ARG, sep  = "_")) %>% 
  column_to_rownames("GN_id") %>% 
  as.matrix() -> tax_table(ps_tmp)

ps_tmp %>% 
  tax_filter(undetected = 0, use_counts = FALSE) %>% 
  # microViz::phyloseq_validate(use_counts = FALSE) %>% 
  microViz::ps_arrange(-Day_of_Treatment) %>% 
  microViz::comp_barplot(
    tax_transform_for_plot = "identity",
    tax_level = "GN",
    label = "Day_of_Treatment", # name an alternative variable to label axes
    n_taxa = 15, # give more taxa unique colours
    merge_other = FALSE, # split the "other" category to display alpha diversity
    bar_width = 0.9, # reduce the bar width to 70% of one row
    bar_outline_colour = "grey5",
    sample_order = "default") +#,# is the default (use NA to remove outlines)
  # facet_by = c("Model", "Antibiotic", "Antibiotic_mg.mL", "Treatment2")) +
  # group_by = c("Model", "Antibiotic", "Antibiotic_mg.mL", "Treatment2"))
  coord_flip() +  theme_light() +
  scale_y_continuous(labels=percent) -> p2


p2 + facet_wrap(Model ~ Antibiotic + Antibiotic_mg.mL + Treatment2, scales = "free", ncol = 6, nrow = 4) -> p2_1

p2_1

p2 + facet_wrap(Model ~ Antibiotic + Antibiotic_mg.mL + Treatment2, scales = "free_y", ncol = 6, nrow = 4) -> p2_2

p2_2

p2_1 %>% 
  ggsave(plot = .,bg = "transparent",scale = 2,
         units = "mm",
         dpi = 600,
         height = 0.618 * 84 * 2,
         width = 84 * 2,
         filename = "~/Desktop/vana_2-1.eps")

p2_2 %>% 
  ggsave(plot = .,bg = "transparent",scale = 2,
         units = "mm",
         dpi = 600,
         height = 0.618 * 84 * 2,
         width = 84 * 2,
         filename = "~/Desktop/vana_2-2.eps")


p2_1 %>%
  export::graph2ppt(append = TRUE, width = 317.48031496 * 2  , height = 0.618 * 317.48031496 * 0.75 , paper = "A4",  scaling = 2,
                    file = out_pptx)

p2_2 %>%
  export::graph2ppt(append = TRUE, width = 317.48031496 * 2  , height = 0.618 * 317.48031496 * 0.75 , paper = "A4",  scaling = 2,
                    file = out_pptx)
```

- CTX-M-1


```r
# ps_AMR_all %>% 
#   physeq_sel_tax_table(c("PathoFact_AMR_AMR_sub_class", "PathoFact_AMR_AMR_sub_class_multi", "PathoFact_AMR_MGE_prediction_clean_2", "PathoFact_AMR_ARG"))  %>%
#   prune_samples(samples, .) %>% 
#   physeq_add_metadata(metadata_mapping,
#                       sample_column = "sample_id") %>% 
#   subset_taxa(PathoFact_AMR_ARG %in% c("CTX-M-1")) %>% 
#   filter_taxa(function(x) sum(x > 0) > 0, TRUE) -> ps_tmp
# 
# ps_tmp %>% 
#   tax_table() %>% 
#   data.frame() %>% 
#   rownames_to_column('GN_id') %>% 
#   mutate(PathoFact_AMR_MGE_prediction_clean_2 = case_when(grepl("ambiguous", PathoFact_AMR_MGE_prediction_clean_2) ~ "ambiguous", TRUE ~ PathoFact_AMR_MGE_prediction_clean_2)) %>% 
#   # mutate(PathoFact_AMR_MGE_prediction_clean_2 = case_when(PathoFact_AMR_MGE_prediction_clean_2 %in% c("MAG", "chromosome")  ~ "MAG - chromosome", TRUE ~ PathoFact_AMR_MGE_prediction_clean_2)) %>% 
#   mutate(PathoFact_AMR_MGE_prediction_clean_2 = fct_relevel(PathoFact_AMR_MGE_prediction_clean_2, c("MAG", "chromosome", "plasmid", "phage", "ambiguous", "unclassified"))) %>% 
#   mutate(GN = paste(GN_id, PathoFact_AMR_MGE_prediction_clean_2, PathoFact_AMR_ARG, sep  = "_")) %>% 
#   column_to_rownames("GN_id") %>% 
#   as.matrix() -> tax_table(ps_tmp)
# 
# ps_tmp %>% 
#   ps_arrange(-Day_of_Treatment) %>% 
#   comp_barplot(
#     tax_transform_for_plot = "identity",
#     tax_level = "GN",
#     label = "Day_of_Treatment", # name an alternative variable to label axes
#     n_taxa = 18, # give more taxa unique colours
#     merge_other = FALSE, # split the "other" category to display alpha diversity
#     bar_width = 0.9, # reduce the bar width to 70% of one row
#     bar_outline_colour = "grey5",
#     sample_order = "default") +#,# is the default (use NA to remove outlines)
#   # facet_by = c("Model", "Antibiotic", "Antibiotic_mg.mL", "Treatment2")) +
#   # group_by = c("Model", "Antibiotic", "Antibiotic_mg.mL", "Treatment2"))
#   coord_flip() +  theme_light() +
#   scale_y_continuous(labels=scientific) -> p2
# # 
# # # Plot them side by side with the patchwork package.
# # patch <- patchwork::wrap_plots(plot_list, nrow = 1, guides = "collect")
# # patch & coord_flip() # make all plots horizontal (note: use & instead of +)
# #   
# 
# p2 + facet_wrap(Model ~ Antibiotic + Antibiotic_mg.mL + Treatment2, scales = "free", ncol = 6, nrow = 4) -> p2_1
# 
# p2_1
# 
# p2 + facet_wrap(Model ~ Antibiotic + Antibiotic_mg.mL + Treatment2, scales = "free_y", ncol = 6, nrow = 4) -> p2_2
# 
# p2_2
# 
# p2_1 %>% 
#   ggsave(plot = .,bg = "transparent",scale = 2,
#          units = "mm",
#          dpi = 600,
#          height = 0.618 * 84 * 2,
#          width = 84 * 2,
#          filename = "~/Desktop/CTX_1.eps")
# 
# p2_2 %>% 
#   ggsave(plot = .,bg = "transparent",scale = 2,
#          units = "mm",
#          dpi = 600,
#          height = 0.618 * 84 * 2,
#          width = 84 * 2,
#          filename = "~/Desktop/CTX_2.eps")
# 
# 
# p2_1 %>%
#   export::graph2ppt(append = TRUE, width = 317.48031496 * 2  , height = 0.618 * 317.48031496 * 0.75 , paper = "A4",  scaling = 2,
#                     file = out_pptx)
# 
# p2_2 %>%
#   export::graph2ppt(append = TRUE, width = 317.48031496 * 2  , height = 0.618 * 317.48031496 * 0.75 , paper = "A4",  scaling = 2,
#                     file = out_pptx)
```

HQ MAG contributors:


```r
ps_AMR_all %>% 
  physeq_sel_tax_table(c("PathoFact_AMR_AMR_sub_class", "PathoFact_AMR_AMR_sub_class_multi", "PathoFact_AMR_MGE_prediction_clean_2", "ANVIO_CONCOCT_HQ_clean", "CONCOCT_bin",
                         "GTDB_tax_CONCOCT_clean" ,"PathoFact_AMR_ARG")) %>% 
  prune_samples(samples, .) %>% 
  physeq_add_metadata(metadata_mapping,
                      sample_column = "sample_id") %>% 
  # subset_taxa(!is.na("GTDB_tax_CONCOCT_clean")) %>% 
  # physeq_sel_tax_table(c("PathoFact_AMR_ARG"))  %>% 
  speedyseq::tax_glom(., taxrank = "GTDB_tax_CONCOCT_clean")  %>% 
  filter_taxa(function(x) sum(x > 0) > 0, TRUE) %>% 
  microbiome::transform("compositional") %>%
  speedyseq::psmelt() %>% 
  separate(GTDB_tax_CONCOCT_clean, into = c("domain", "phylum", "class", "order", "family", "genus", "species"), sep  = ";" , remove = FALSE) -> df
```

```
## Warning: Expected 7 pieces. Missing pieces filled with `NA` in 510 rows [3618,
## 3765, 3781, 3830, 4292, 4364, 5258, 6054, 6406, 6758, 6922, 7492, 7494, 7498,
## 8076, 8426, 8739, 8971, 9028, 9215, ...].
```



```r
df %>% 
  filter(genus %in% c("g__Enterococcus_B",
                      "g__Escherichia")) %>%
  mutate(Abundance = na_if(Abundance, 0)) %>%
  mutate(logTPM = log10(Abundance + 1))  %>% 
  # separate_rows(., Resistance_Mechanism,sep = '; ',convert = TRUE) %>% 
  # ggplot(aes(x = as.factor(Day_of_Treatment), y = Abundance, fill = fct_reorder(Strain, logTPM))) +   
  ggplot(aes(x = as.factor(Day_of_Treatment), y = Abundance, fill = species)) +#reorder(genus, Abundance, FUN=median, na.rm = TRUE))) +
  geom_bar( stat = "identity", colour="black", size=0.025) +
  facet_wrap(Model ~ Antibiotic + Antibiotic_mg.mL + Treatment2, scales = "fixed", ncol = 6, nrow = 4) +
  # facet_grid(Model ~ Treatment + Fermentation + Antibiotic_mg.mL , scales = "free", space = "free_x", drop = TRUE) +
  ggtitle("ARG Type Abundance per Sample") +
  xlab("Day (Treatment)")  +
  # ylab("Proportion") + 
  theme_light() +
  scale_y_continuous(labels=percent) + 
  # scale_fill_manual(name = "", values = col_strains) +
  theme(legend.position = "bottom") +
  ggpubr::rotate_x_text(60) +
  theme(legend.position = ("bottom")) -> p_mec

# p_mec %>% 
#   ggpubr::set_palette("jco") -> p_mec

p_mec
```

```
## Warning: Removed 99 rows containing missing values (position_stack).
```

<img src="chicken1_resistome_plots_draft_files/figure-html/unnamed-chunk-83-1.png" style="display: block; margin: auto;" />

```r
p_mec %>% 
  export::graph2ppt(append = TRUE, width = 317.48031496 * 2  , height = 0.618 * 317.48031496 * 0.75 , paper = "A4",  scaling = 2,
                    file = out_pptx)
```

```
## Warning: Removed 99 rows containing missing values (position_stack).
```

```
## Exported graph as ~/Desktop/plots_manuscript_chicken1.pptx
```

Need to clarify resistome contributors since some gene are assigned to bins but are lacking GTDB taxonomy. updated collection?

contigs of interest: 
c1_megahit_8588
c1_megahit_219391



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
##  date     2022-09-28
##  pandoc   2.18 @ /Applications/RStudio.app/Contents/MacOS/quarto/bin/tools/ (via rmarkdown)
## 
## ─ Packages ───────────────────────────────────────────────────────────────────
##  package          * version    date (UTC) lib source
##  abind              1.4-5      2016-07-21 [1] CRAN (R 4.2.0)
##  ade4               1.7-19     2022-04-19 [1] CRAN (R 4.2.0)
##  ape                5.6-2      2022-03-02 [1] CRAN (R 4.2.0)
##  assertthat         0.2.1      2019-03-21 [1] CRAN (R 4.2.0)
##  backports          1.4.1      2021-12-13 [1] CRAN (R 4.2.0)
##  base64enc          0.1-3      2015-07-28 [1] CRAN (R 4.2.0)
##  bayesm             3.1-4      2019-10-15 [1] CRAN (R 4.2.0)
##  Biobase            2.56.0     2022-04-26 [1] Bioconductor
##  BiocGenerics       0.42.0     2022-04-26 [1] Bioconductor
##  biomformat         1.24.0     2022-04-26 [1] Bioconductor
##  Biostrings         2.64.0     2022-04-26 [1] Bioconductor
##  bit                4.0.4      2020-08-04 [1] CRAN (R 4.2.0)
##  bit64              4.0.5      2020-08-30 [1] CRAN (R 4.2.0)
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
##  compositions     * 2.0-4      2022-01-05 [1] CRAN (R 4.2.0)
##  crayon             1.5.1      2022-03-26 [1] CRAN (R 4.2.0)
##  crosstalk          1.2.0      2021-11-04 [1] CRAN (R 4.2.0)
##  data.table         1.14.2     2021-09-27 [1] CRAN (R 4.2.0)
##  DBI                1.1.3      2022-06-18 [1] CRAN (R 4.2.0)
##  dbplyr             2.2.1      2022-06-27 [1] CRAN (R 4.2.0)
##  DEoptimR           1.0-11     2022-04-03 [1] CRAN (R 4.2.0)
##  devEMF             4.1        2022-06-24 [1] CRAN (R 4.2.0)
##  devtools           2.4.3      2021-11-30 [1] CRAN (R 4.2.0)
##  digest             0.6.29     2021-12-01 [1] CRAN (R 4.2.0)
##  dplyr            * 1.0.9      2022-04-28 [1] CRAN (R 4.2.0)
##  DT                 0.23       2022-05-10 [1] CRAN (R 4.2.0)
##  dtplyr             1.2.1      2022-01-19 [1] CRAN (R 4.2.0)
##  ellipsis           0.3.2      2021-04-29 [1] CRAN (R 4.2.0)
##  evaluate           0.15       2022-02-18 [1] CRAN (R 4.2.0)
##  export             0.3.0      2022-06-13 [1] Github (tomwenseleers/export@1afc8e2)
##  extrafont          0.18       2022-04-12 [1] CRAN (R 4.2.0)
##  extrafontdb        1.0        2012-06-11 [1] CRAN (R 4.2.0)
##  fansi              1.0.3      2022-03-24 [1] CRAN (R 4.2.0)
##  farver             2.1.1      2022-07-06 [1] CRAN (R 4.2.0)
##  fastmap            1.1.0      2021-01-25 [1] CRAN (R 4.2.0)
##  flextable          0.7.2      2022-06-12 [1] CRAN (R 4.2.0)
##  forcats          * 0.5.1      2021-01-27 [1] CRAN (R 4.2.0)
##  foreach            1.5.2      2022-02-02 [1] CRAN (R 4.2.0)
##  fs                 1.5.2      2021-12-08 [1] CRAN (R 4.2.0)
##  gargle             1.2.0      2021-07-02 [1] CRAN (R 4.2.0)
##  gdtools          * 0.2.4      2022-02-14 [1] CRAN (R 4.2.0)
##  generics           0.1.3      2022-07-05 [1] CRAN (R 4.2.0)
##  GenomeInfoDb       1.32.2     2022-05-15 [1] Bioconductor
##  GenomeInfoDbData   1.2.8      2022-06-13 [1] Bioconductor
##  ggplot2          * 3.3.6      2022-05-03 [1] CRAN (R 4.2.0)
##  ggpubr             0.4.0      2020-06-27 [1] CRAN (R 4.2.0)
##  ggrepel            0.9.1      2021-01-15 [1] CRAN (R 4.2.0)
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
##  htmltools          0.5.3      2022-07-18 [1] CRAN (R 4.2.0)
##  htmlwidgets        1.5.4      2021-09-08 [1] CRAN (R 4.2.0)
##  httr               1.4.3      2022-05-04 [1] CRAN (R 4.2.0)
##  igraph             1.3.4      2022-07-19 [1] CRAN (R 4.2.0)
##  IRanges            2.30.0     2022-04-26 [1] Bioconductor
##  iterators          1.0.14     2022-02-05 [1] CRAN (R 4.2.0)
##  jquerylib          0.1.4      2021-04-26 [1] CRAN (R 4.2.0)
##  jsonlite           1.8.0      2022-02-22 [1] CRAN (R 4.2.0)
##  knitr              1.39       2022-04-26 [1] CRAN (R 4.2.0)
##  labeling           0.4.2      2020-10-20 [1] CRAN (R 4.2.0)
##  lattice          * 0.20-45    2021-09-22 [1] CRAN (R 4.2.1)
##  lifecycle          1.0.1      2021-09-24 [1] CRAN (R 4.2.0)
##  lubridate          1.8.0      2021-10-07 [1] CRAN (R 4.2.0)
##  magrittr           2.0.3      2022-03-30 [1] CRAN (R 4.2.0)
##  MASS               7.3-58     2022-07-14 [1] CRAN (R 4.2.1)
##  Matrix             1.4-1      2022-03-23 [1] CRAN (R 4.2.1)
##  memoise            2.0.1      2021-11-26 [1] CRAN (R 4.2.0)
##  metagMisc        * 0.0.4      2022-06-14 [1] Github (vmikk/metagMisc@9cb3131)
##  mgcv               1.8-40     2022-03-29 [1] CRAN (R 4.2.1)
##  microbiome       * 1.19.1     2022-08-03 [1] Github (microbiome/microbiome@f1a7671)
##  microViz         * 0.9.2      2022-07-08 [1] Github (david-barnett/microViz@4d5f969)
##  modelr             0.1.8      2020-05-19 [1] CRAN (R 4.2.0)
##  multtest           2.52.0     2022-04-26 [1] Bioconductor
##  munsell            0.5.0      2018-06-12 [1] CRAN (R 4.2.0)
##  nlme             * 3.1-158    2022-06-15 [1] CRAN (R 4.2.0)
##  officer            0.4.3      2022-06-12 [1] CRAN (R 4.2.0)
##  openxlsx           4.2.5      2021-12-14 [1] CRAN (R 4.2.0)
##  permute          * 0.9-7      2022-01-27 [1] CRAN (R 4.2.0)
##  pheatmap           1.0.12     2019-01-04 [1] CRAN (R 4.2.0)
##  phyloseq         * 1.39.1     2022-08-03 [1] Github (joey711/phyloseq@9b211d9)
##  pillar             1.8.0      2022-07-18 [1] CRAN (R 4.2.0)
##  pkgbuild           1.3.1      2021-12-20 [1] CRAN (R 4.2.0)
##  pkgconfig          2.0.3      2019-09-22 [1] CRAN (R 4.2.0)
##  pkgdown            2.0.5      2022-06-23 [1] CRAN (R 4.2.0)
##  pkgload            1.3.0      2022-06-27 [1] CRAN (R 4.2.0)
##  plyr               1.8.7      2022-03-24 [1] CRAN (R 4.2.0)
##  prettyunits        1.1.1      2020-01-24 [1] CRAN (R 4.2.0)
##  processx           3.7.0      2022-07-07 [1] CRAN (R 4.2.0)
##  ps                 1.7.1      2022-06-18 [1] CRAN (R 4.2.0)
##  purrr            * 0.3.4      2020-04-17 [1] CRAN (R 4.2.0)
##  R6                 2.5.1      2021-08-19 [1] CRAN (R 4.2.0)
##  ragg               1.2.2      2022-02-21 [1] CRAN (R 4.2.0)
##  RColorBrewer       1.1-3      2022-04-03 [1] CRAN (R 4.2.0)
##  Rcpp               1.0.9      2022-07-08 [1] CRAN (R 4.2.0)
##  RCurl              1.98-1.8   2022-07-30 [1] CRAN (R 4.2.0)
##  readr            * 2.1.2      2022-01-30 [1] CRAN (R 4.2.0)
##  readxl             1.4.0      2022-03-28 [1] CRAN (R 4.2.0)
##  remotes            2.4.2      2021-11-30 [1] CRAN (R 4.2.0)
##  reprex             2.0.1      2021-08-05 [1] CRAN (R 4.2.0)
##  reshape2         * 1.4.4      2020-04-09 [1] CRAN (R 4.2.0)
##  rgl                0.108.3    2021-11-21 [1] CRAN (R 4.2.0)
##  rhdf5              2.40.0     2022-04-26 [1] Bioconductor
##  rhdf5filters       1.8.0      2022-04-26 [1] Bioconductor
##  Rhdf5lib           1.18.2     2022-05-17 [1] Bioconductor
##  rJava              1.0-6      2021-12-10 [1] CRAN (R 4.2.0)
##  rlang              1.0.4      2022-07-12 [1] CRAN (R 4.2.0)
##  rmarkdown          2.14       2022-04-25 [1] CRAN (R 4.2.0)
##  robustbase         0.95-0     2022-04-02 [1] CRAN (R 4.2.0)
##  rprojroot          2.0.3      2022-04-02 [1] CRAN (R 4.2.0)
##  rstatix            0.7.0      2021-02-13 [1] CRAN (R 4.2.0)
##  rstudioapi         0.13       2020-11-12 [1] CRAN (R 4.2.0)
##  Rtsne              0.16       2022-04-17 [1] CRAN (R 4.2.0)
##  Rttf2pt1           1.3.10     2022-02-07 [1] CRAN (R 4.2.0)
##  rvest              1.0.2      2021-10-16 [1] CRAN (R 4.2.0)
##  rvg                0.2.5      2020-06-30 [1] CRAN (R 4.2.0)
##  S4Vectors          0.34.0     2022-04-26 [1] Bioconductor
##  sass               0.4.1      2022-03-23 [1] CRAN (R 4.2.0)
##  scales           * 1.2.0      2022-04-13 [1] CRAN (R 4.2.0)
##  sessioninfo        1.2.2      2021-12-06 [1] CRAN (R 4.2.0)
##  speedyseq        * 0.5.3.9018 2022-06-13 [1] Github (mikemc/speedyseq@ceb941f)
##  stargazer          5.2.3      2022-03-04 [1] CRAN (R 4.2.0)
##  stringi            1.7.8      2022-07-11 [1] CRAN (R 4.2.0)
##  stringr          * 1.4.0      2019-02-10 [1] CRAN (R 4.2.0)
##  survival           3.3-1      2022-03-03 [1] CRAN (R 4.2.1)
##  systemfonts        1.0.4      2022-02-11 [1] CRAN (R 4.2.0)
##  tensorA            0.36.2     2020-11-19 [1] CRAN (R 4.2.0)
##  textshaping        0.3.6      2021-10-13 [1] CRAN (R 4.2.0)
##  tibble           * 3.1.8      2022-07-22 [1] CRAN (R 4.2.0)
##  tidyr            * 1.2.0      2022-02-01 [1] CRAN (R 4.2.0)
##  tidyselect         1.1.2      2022-02-21 [1] CRAN (R 4.2.0)
##  tidyverse        * 1.3.1.9000 2022-06-13 [1] Github (tidyverse/tidyverse@6186fbf)
##  tzdb               0.3.0      2022-03-28 [1] CRAN (R 4.2.0)
##  usethis            2.1.6      2022-05-25 [1] CRAN (R 4.2.0)
##  utf8               1.2.2      2021-07-24 [1] CRAN (R 4.2.0)
##  uuid               1.1-0      2022-04-19 [1] CRAN (R 4.2.0)
##  vctrs              0.4.1      2022-04-13 [1] CRAN (R 4.2.0)
##  vegan            * 2.6-2      2022-04-17 [1] CRAN (R 4.2.0)
##  vroom              1.5.7      2021-11-30 [1] CRAN (R 4.2.0)
##  withr              2.5.0      2022-03-03 [1] CRAN (R 4.2.0)
##  xfun               0.31       2022-05-10 [1] CRAN (R 4.2.0)
##  xlsx               0.6.5      2020-11-10 [1] CRAN (R 4.2.0)
##  xlsxjars           0.6.1      2014-08-22 [1] CRAN (R 4.2.0)
##  xml2               1.3.3      2021-11-30 [1] CRAN (R 4.2.0)
##  xtable             1.8-4      2019-04-21 [1] CRAN (R 4.2.0)
##  XVector            0.36.0     2022-04-26 [1] Bioconductor
##  yaml               2.3.5      2022-02-21 [1] CRAN (R 4.2.0)
##  zip                2.2.0      2021-05-31 [1] CRAN (R 4.2.0)
##  zlibbioc           1.42.0     2022-04-26 [1] Bioconductor
## 
##  [1] /Library/Frameworks/R.framework/Versions/4.2/Resources/library
## 
## ──────────────────────────────────────────────────────────────────────────────
```



- the top gene -> where line plot facet localisation? heatmap?
- heatmap AMR class + genetic material + dendro
- heatmap AMR GN + dendro


Resistome Proportion
AMR categories
Gene heatmap
Diversity: alpha/beta
Contributing genetic material (MAGs/Chromosomes/plasmids/…)
Number of AMR/ contigs / plasmids – length
Number of unmapped reads / sample
Prop/coverage or whatsoever per gene per mag per sample (rector/time)
Visualise contigs with ARG in anvo -> coverage variability. 

proof of concept with WAFLE METACHIP
proof of concept ANVIO -> population variability

MAG resistome?
donor strain chromosomal resistome?
