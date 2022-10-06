---
title: "Resistome data chicken1 - draft"
author: "Florentin Constancias "
date: " September 27, 2022 "
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
## Ncells 525695 28.1    1167338 62.4         NA   669294 35.8
## Vcells 990944  7.6    8388608 64.0      16384  1840039 14.1
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
out_pptx = "~/Desktop/plots_draft_chicken1.pptx"
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

## Aestetics:


```r
load( here::here("Figures/Rastetics.Rdata"))
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
  mutate(Antibiotic_mg.L = factor(Antibiotic_mg.mL, levels = c(NA, 20, 200,90, 600))) -> metadata_mapping
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
  mutate(Antibiotic_mg.mL = factor(Antibiotic_mg.mL, levels = c(NA, 20, 200,90, 600))) %>%
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

<img src="combine_chicken1_resistome_plots_draft_files/figure-html/unnamed-chunk-11-1.png" style="display: block; margin: auto;" />

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
## Exported graph as ~/Desktop/plots_draft_chicken1.pptx
```

Also retrieve the samples we now keep for the final analysis:


```r
p1$data %>% distinct(Sample)  %>% pull %>%  as.character() -> samples
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

<img src="combine_chicken1_resistome_plots_draft_files/figure-html/unnamed-chunk-13-1.png" style="display: block; margin: auto;" />

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
## Exported graph as ~/Desktop/plots_draft_chicken1.pptx
```


```r
p1 + facet_null() + facet_grid(as.formula(paste0("Antibiotic ~  Model  ")), 
                               scales = "free", space = "fixed", drop = TRUE) -> p1_2

p1_2
```

<img src="combine_chicken1_resistome_plots_draft_files/figure-html/unnamed-chunk-14-1.png" style="display: block; margin: auto;" />

```r
p1_2 %>% 
  ggsave(plot = .,bg = "transparent",scale = 2,
         units = "mm",
         dpi = 600,
         height = 0.618 * 84,
         width = 84,
         filename = "~/Desktop/Resistome_baseline2.eps")

p1_2 %>%
  export::graph2ppt(append = TRUE, width = 317.48031496 / 1.25 , height = 0.618 * 317.48031496 / 1.25 , paper = "A4",  scaling = 2,
                    file = out_pptx)
```

```
## Exported graph as ~/Desktop/plots_draft_chicken1.pptx
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
  ps_arrange(-Day_of_Treatment) %>% 
  comp_barplot(
    tax_level = "PathoFact_AMR_AMR_sub_class_multi",
    label = "Day_of_Treatment", # name an alternative variable to label axes
    n_taxa = 15, # give more taxa unique colours
    merge_other = FALSE, # split the "other" category to display alpha diversity
    bar_width = 0.9, # reduce the bar width to 70% of one row
    bar_outline_colour = "grey5",
    sample_order = "default") +#,# is the default (use NA to remove outlines)
  # facet_by = c("Model", "Antibiotic", "Antibiotic_mg.mL", "Treatment2")) +
  # group_by = c("Model", "Antibiotic", "Antibiotic_mg.mL", "Treatment2"))
  coord_flip() +  theme_light() + 
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
# # Plot them side by side with the patchwork package.
# patch <- patchwork::wrap_plots(plot_list, nrow = 1, guides = "collect")
# patch & coord_flip() # make all plots horizontal (note: use & instead of +)
```



```r
plot_AMR_multi <- p2


p2
```

<img src="combine_chicken1_resistome_plots_draft_files/figure-html/unnamed-chunk-16-1.png" style="display: block; margin: auto;" />

```r
p2 %>% 
  ggsave(plot = .,bg = "transparent",scale = 2,
         units = "mm",
         dpi = 600,
         height = 0.618 * 84 * 2,
         width = 84 * 2,
         filename = "~/Desktop/Resistome_sub_class_multi_prop.eps")

p2 %>%
  export::graph2ppt(append = TRUE, width = 317.48031496 * 2  , height = 0.618 * 317.48031496 * 1.25 , paper = "A4",  scaling = 2,
                    file = out_pptx)
```

```
## Exported graph as ~/Desktop/plots_draft_chicken1.pptx
```


```r
ps_AMR_all %>%
  physeq_sel_tax_table(c("PathoFact_AMR_AMR_sub_class", "PathoFact_AMR_AMR_sub_class_multi", "PathoFact_AMR_ARG")) %>%
  speedyseq::tax_glom(., taxrank = "PathoFact_AMR_AMR_sub_class_multi") %>% 
  physeq_sel_tax_table(c("PathoFact_AMR_AMR_sub_class_multi")) %>% 
  prune_samples(samples, .) %>% 
  filter_taxa(function(x) sum(x > 0) > 0, TRUE) %>% 
  physeq_add_metadata(metadata_mapping,
                      sample_column = "sample_id") %>% 
  subset_taxa(PathoFact_AMR_AMR_sub_class_multi %!in% c("Unclassified", "multidrug (unclassified)")) %>% 
  # ps_mutate(Day_of_Treatment = as.numeric(Day_of_Treatment)) %>% 
  ps_arrange(-Day_of_Treatment) %>% 
  comp_barplot(
    tax_transform_for_plot = "identity",
    tax_level = "PathoFact_AMR_AMR_sub_class_multi",
    label = "Day_of_Treatment", # name an alternative variable to label axes
    n_taxa = 15, # give more taxa unique colours
    merge_other = FALSE, # split the "other" category to display alpha diversity
    bar_width = 0.9, # reduce the bar width to 70% of one row
    bar_outline_colour = "grey5",
    sample_order = "default") +#,# is the default (use NA to remove outlines)
  # facet_by = c("Model", "Antibiotic", "Antibiotic_mg.mL", "Treatment2")) +
  # group_by = c("Model", "Antibiotic", "Antibiotic_mg.mL", "Treatment2"))
  coord_flip() +  theme_light() + 
  scale_y_continuous(labels=scientific) + 
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
p2
```

<img src="combine_chicken1_resistome_plots_draft_files/figure-html/unnamed-chunk-17-1.png" style="display: block; margin: auto;" />

```r
p2 %>% 
  ggsave(plot = .,bg = "transparent",scale = 2,
         units = "mm",
         dpi = 600,
         height = 0.618 * 84 * 2,
         width = 84 * 2,
         filename = "~/Desktop/Resistome_sub_class_multi.eps")

p2 %>%
  export::graph2ppt(append = TRUE, width = 317.48031496 * 2  , height = 0.618 * 317.48031496 * 1.25 , paper = "A4",  scaling = 2,
                    file = out_pptx)
```

```
## Exported graph as ~/Desktop/plots_draft_chicken1.pptx
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
  ps_arrange(-Day_of_Treatment) %>% 
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
  summarise(meanAbundance = mean(Abundance, na.rm = TRUE)) -> mean_ab


pheatmap_samples %>%
  replace(is.na(.), 0)  %>%
  t() %>%
  vegan::vegdist(na.rm = TRUE) %>%
  hclust() -> clust_t

pheatmap_samples %>%
  replace(is.na(.), 0)  %>%
  vegan::vegdist(na.rm = TRUE) %>%
  # cor() %>%
  hclust() -> clust

pheatmap_samples %>%
  pheatmap::pheatmap(cluster_rows = FALSE ,
                     scale = "row",
                     cluster_cols = clust_t,
                     na_col = "transparent",
                     annotation_row = tax_mapping %>% rownames_to_column('tmp_col') %>%
                       filter(PathoFact_AMR_AMR_sub_class_multi %in% rownames(pheatmap_samples)) %>%
                       distinct(PathoFact_AMR_AMR_sub_class_multi, .keep_all = TRUE) %>%
                       left_join(., 
                                 mean_ab,
                                 by = c("PathoFact_AMR_AMR_sub_class_multi" = "PathoFact_AMR_AMR_sub_class_multi")) %>% 
                       column_to_rownames("PathoFact_AMR_AMR_sub_class_multi") %>%
                       select(PathoFact_AMR_Resistance_mechanism_multi, meanAbundance),
                     annotation_colors = list("PathoFact_AMR_AMR_sub_class_multi" = resistance_type,
                                              "Treatment" = treat_col),
                     annotation_col = pheatmap_df %>%
                       select(Day_of_Treatment, Antibiotic_mg.mL, Treatment, Model, Sample) %>% 
                       distinct(Sample, .keep_all = TRUE) %>% 
                       column_to_rownames("Sample"))  -> p_pheat


p_pheat
```

<img src="combine_chicken1_resistome_plots_draft_files/figure-html/unnamed-chunk-18-1.png" style="display: block; margin: auto;" />

```r
p_pheat %>% 
  ggsave(plot = .,bg = "transparent",scale = 2,
         units = "mm",
         dpi = 600,
         height = 0.618 * 84 * 2,
         width = 84  * 2,
         filename = "~/Desktop/Resistome_heatmap_class.eps")

p_pheat %>%
  export::graph2ppt(append = TRUE, width = 317.48031496 * 8, height = 0.618 * 317.48031496 * 10 , paper = "A4",  scaling = 1,
                    file = out_pptx)
```

```
## Exported graph as ~/Desktop/plots_draft_chicken1.pptx
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

<img src="combine_chicken1_resistome_plots_draft_files/figure-html/unnamed-chunk-19-1.png" style="display: block; margin: auto;" />

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
## Exported graph as ~/Desktop/plots_draft_chicken1.pptx
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

<img src="combine_chicken1_resistome_plots_draft_files/figure-html/unnamed-chunk-20-1.png" style="display: block; margin: auto;" />

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
## Exported graph as ~/Desktop/plots_draft_chicken1.pptx
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
  vegan::vegdist(na.rm = TRUE, method = "euc") %>%
  # cor() %>%
  hclust() -> clust

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

<img src="combine_chicken1_resistome_plots_draft_files/figure-html/unnamed-chunk-21-1.png" style="display: block; margin: auto;" />

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
## Exported graph as ~/Desktop/plots_draft_chicken1.pptx
```



# Beta-div resistome:

Generate a temporary phyloseq object first - containing PathoFact_AMR_ARG - which will be used for the 2 models:
*TODO*: filter based on BP's expertise? Unclassified... ?


```r
ps_AMR_all %>% 
  physeq_sel_tax_table(c("PathoFact_AMR_AMR_sub_class", "PathoFact_AMR_AMR_sub_class_multi", "PathoFact_AMR_ARG"))  %>%
  speedyseq::tax_glom(., taxrank = "PathoFact_AMR_ARG")  %>% 
  filter_taxa(function(x) sum(x > 0) > 0, TRUE) %>% 
  prune_samples(samples, .) %>% 
  physeq_add_metadata(metadata_mapping,
                      sample_column = "sample_id") -> ps_tmp
```

## Chicken:

### Generate Aitchinson PCA:


```r
ps_tmp %>%
  subset_samples(Model == "Chicken") %>% 
  filter_taxa(function(x) sum(x > 0) > 0, TRUE) %>%
  phyloseq_plot_bdiv(dlist = NULL,
                     m = "CoDa",
                     seed = 123,
                     axis1 = 1,
                     axis2 = 2) -> coda_resistome
```

### Initial PCA:


```r
coda_resistome$PCA$layers[[1]] = NULL

coda_resistome$PCA + geom_point(size=2,
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

<img src="combine_chicken1_resistome_plots_draft_files/figure-html/unnamed-chunk-24-1.png" style="display: block; margin: auto;" />



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
## Exported graph as ~/Desktop/plots_draft_chicken1.pptx
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
phyloseq_add_taxa_vector_fix(dist = coda_resistome$aidist,
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

```
## Warning: ggrepel: 1 unlabeled data points (too many overlaps). Consider
## increasing max.overlaps
```

<img src="combine_chicken1_resistome_plots_draft_files/figure-html/unnamed-chunk-28-1.png" style="display: block; margin: auto;" />

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
## Exported graph as ~/Desktop/plots_draft_chicken1.pptx
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

<img src="combine_chicken1_resistome_plots_draft_files/figure-html/unnamed-chunk-30-1.png" style="display: block; margin: auto;" />

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
## Exported graph as ~/Desktop/plots_draft_chicken1.pptx
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
## $ Antibiotic_mg.L             <fct> 90, 90, NA, 90, 90, NA, 20, 20, NA, 90, NA…
## $ Read_nrom                   <dbl> 1e+06, 1e+06, 1e+06, 1e+06, 1e+06, 1e+06, …
## $ id...45                     <chr> "C1", "C6", "D10", "D11", "D12", "D13", "D…
## $ observed                    <dbl> 287, 219, 264, 265, 287, 296, 241, 201, 25…
## $ chao1                       <dbl> 287.00, 219.00, 264.00, 265.00, 287.00, 29…
## $ diversity_inverse_simpson   <dbl> 42.42549, 65.93051, 53.57182, 69.64815, 83…
## $ diversity_gini_simpson      <dbl> 0.9764293, 0.9848325, 0.9813335, 0.9856421…
## $ diversity_shannon           <dbl> 4.392901, 4.450068, 4.553329, 4.501119, 4.…
## $ diversity_fisher            <dbl> 27.31195, 20.26535, 24.90449, 25.00868, 27…
## $ diversity_coverage          <dbl> 18, 25, 25, 27, 31, 24, 6, 7, 29, 23, 29, …
## $ evenness_camargo            <dbl> 0.2395164, 0.3723076, 0.2938254, 0.2888339…
## $ evenness_pielou             <dbl> 0.7762018, 0.8257579, 0.8166017, 0.8066912…
## $ evenness_simpson            <dbl> 0.14782400, 0.30105254, 0.20292356, 0.2628…
## $ evenness_evar               <dbl> 0.16246759, 0.09178509, 0.13392919, 0.0788…
## $ evenness_bulla              <dbl> 0.3787306, 0.4134235, 0.4746145, 0.3892407…
## $ dominance_dbp               <dbl> 0.06587407, 0.02455990, 0.05591778, 0.0286…
## $ dominance_dmn               <dbl> 0.12368312, 0.04781481, 0.10536258, 0.0515…
## $ dominance_absolute          <dbl> 65874, 24560, 55918, 28700, 31274, 54079, …
## $ dominance_relative          <dbl> 0.06587407, 0.02455990, 0.05591778, 0.0286…
## $ dominance_simpson           <dbl> 0.023570737, 0.015167485, 0.018666530, 0.0…
## $ dominance_core_abundance    <dbl> 0.9114019, 0.9138453, 0.9091184, 0.9445082…
## $ dominance_gini              <dbl> 0.8509026, 0.8538181, 0.8262136, 0.8458745…
## $ rarity_log_modulo_skewness  <dbl> 2.033017, 2.005974, 1.863853, 1.978196, 2.…
## $ rarity_low_abundance        <dbl> 0.11715412, 0.04722881, 0.04845981, 0.0375…
## $ rarity_rare_abundance       <dbl> 0.04817305, 0.04063284, 0.03826785, 0.0181…
## $ id...68                     <chr> "C1", "C6", "D10", "D11", "D12", "D13", "D…
## $ Observed                    <dbl> 287, 219, 264, 265, 287, 296, 241, 201, 25…
## $ Chao1                       <dbl> 287.00, 219.00, 264.00, 265.00, 287.00, 29…
## $ se.chao1                    <dbl> 0.0000000, 0.0000000, 0.0000000, 0.0000000…
## $ ACE                         <dbl> 287.0000, 219.0000, 264.0000, 265.0000, 28…
## $ se.ACE                      <dbl> 1.409277, 4.165525, 2.421495, 5.767181, 3.…
```


```r
plot_alpha_div_NRP72 <- function(df = resistome_alpha, x = "Day_of_Treatment", y = "value",  point_size = 2.5, color = "Treatment", fill  = "Treatment", shape = "Antibiotic", alpha = "Antibiotic_mg.mL", facet_formula = paste0("alphadiversiy ~  Model "), ylab = "Resistome alpha-diversity", xlab = "Days (Treatment)", measures = c("Observed", "evenness_pielou"), path_group = "interaction(Model, Reactor_Treatment_Dose)"){
  
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
  theme_light() + facet_null() + facet_grid(Antibiotic + alphadiversiy  ~ Model, scales = "free", drop = TRUE) +
  ggh4x::facetted_pos_scales(
    y = list(
      alphadiversiy == "Observed" ~ scale_y_continuous(limits = c(120,270)),
      alphadiversiy == "evenness_pielou" ~ scale_y_continuous(limits = c(0.6, 1))
    )
  ) -> p_alpha

p_alpha
```

<img src="combine_chicken1_resistome_plots_draft_files/figure-html/unnamed-chunk-36-1.png" style="display: block; margin: auto;" />

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
## Exported graph as ~/Desktop/plots_draft_chicken1.pptx
```

Same with subtraction from baseline.


```r
x = "Day_of_Treatment"
y = "value"
point_size = 2.5
color = "Treatment"
fill  = "Treatment"
shape = "as.factor(Antibiotic_mg.mL)"
alpha = NULL
facet_formula = paste0("alphadiversiy ~  Model ")
ylab = "Resistome baseline delta"
xlab = "Days (Treatment)"
measures = c("Observed", "evenness_pielou")
path_group = "interaction(Model, Reactor_Treatment_Dose)"


resistome_alpha %>% 
  pivot_longer(cols = all_of(measures), values_to = "value", names_to = 'alphadiversiy', values_drop_na  = TRUE) %>%
  group_by(alphadiversiy, Model,Fermentation,  Reactor_Treatment_Dose) %>% 
  arrange(Day_of_Treatment) %>% 
  mutate(delta_change = value - first(value),
         ratio_change =  value / first(value)) %>% 
  ggplot(aes_string(x = x, y = "delta_change", color = color, fill = fill, shape = shape, alpha = alpha)) + #shape = Fermentation
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
  geom_vline(xintercept = 0, size = 0.5, linetype = 1, color = "grey60")+
  scale_color_manual(name = "", values = treat_col,
                     na.value = "black") +
  scale_fill_manual(name = "", values = treat_col,
                    na.value = "black") +
  scale_shape_manual(name = "" ,values = c(15,16,18,19), na.value =  17) + 
  theme_light() + facet_null() +  facet_grid(Antibiotic + alphadiversiy  ~ Model , scales = "free", drop = TRUE) +
  ggh4x::facetted_pos_scales(
    y = list(
      alphadiversiy == "Observed" ~ scale_y_continuous(limits = c(-125,50)),
      alphadiversiy == "evenness_pielou" ~ scale_y_continuous(limits = c(-0.2, 0.2))
    )
  ) -> p_alpha_delta


p_alpha_delta
```

<img src="combine_chicken1_resistome_plots_draft_files/figure-html/unnamed-chunk-37-1.png" style="display: block; margin: auto;" />

```r
p_alpha_delta %>% 
  ggsave(plot = .,bg = "transparent",scale = 2,
         units = "mm",
         dpi = 600,
         height = 0.618 * 84,
         width = 84,
         filename = "~/Desktop/palpha_d.eps")

p_alpha_delta %>%
  export::graph2ppt(append = TRUE, width = 317.48031496  , height = 0.618 * 317.48031496  , paper = "A4",  scaling = 2,
                    file = out_pptx)
```

```
## Exported graph as ~/Desktop/plots_draft_chicken1.pptx
```

Same with ratio from baseline.


```r
x = "Day_of_Treatment"
y = "value"
point_size = 2.5
color = "Treatment"
fill  = "Treatment"
shape = "as.factor(Antibiotic_mg.mL)"
alpha = NULL
facet_formula = paste0("alphadiversiy ~  Model ")
ylab = "Resistome baseline ratio"
xlab = "Days (Treatment)"
measures = c("Observed", "evenness_pielou")
path_group = "interaction(Model, Reactor_Treatment_Dose)"


resistome_alpha %>% 
  pivot_longer(cols = all_of(measures), values_to = "value", names_to = 'alphadiversiy', values_drop_na  = TRUE) %>%
  group_by(alphadiversiy, Model,Fermentation,  Reactor_Treatment_Dose) %>% 
  arrange(Day_of_Treatment) %>% 
  mutate(delta_change = value - first(value),
         ratio_change =  value / first(value)) %>% 
  ggplot(aes_string(x = x, y = "ratio_change", color = color, fill = fill, shape = shape, alpha = alpha)) + #shape = Fermentation
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
  geom_vline(xintercept = 0, size = 0.5, linetype = 1, color = "grey60")+
  scale_color_manual(name = "", values = treat_col,
                     na.value = "black") +
  scale_fill_manual(name = "", values = treat_col,
                    na.value = "black") +
  scale_shape_manual(name = "" ,values = c(15,16,18,19), na.value =  17) + 
  theme_light() + facet_null() +  facet_grid(Antibiotic + alphadiversiy  ~ Model , scales = "free", drop = TRUE) +
  ggh4x::facetted_pos_scales(
    y = list(
      alphadiversiy == "Observed" ~ scale_y_continuous(limits = c(0.5,1.2)),
      alphadiversiy == "evenness_pielou" ~ scale_y_continuous(limits = c(0.5, 1.4))
    )
  )-> p_alpha_delta


p_alpha_delta
```

<img src="combine_chicken1_resistome_plots_draft_files/figure-html/unnamed-chunk-38-1.png" style="display: block; margin: auto;" />

```r
p_alpha_delta %>% 
  ggsave(plot = .,bg = "transparent",scale = 2,
         units = "mm",
         dpi = 600,
         height = 0.618 * 84,
         width = 84,
         filename = "~/Desktop/palpha_r.eps")

p_alpha_delta %>%
  export::graph2ppt(append = TRUE, width = 317.48031496  , height = 0.618 * 317.48031496  , paper = "A4",  scaling = 2,
                    file = out_pptx)
```

```
## Exported graph as ~/Desktop/plots_draft_chicken1.pptx
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
  coord_flip() +  theme_light() +
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

<img src="combine_chicken1_resistome_plots_draft_files/figure-html/unnamed-chunk-39-1.png" style="display: block; margin: auto;" />

```r
p2 %>% 
  ggsave(plot = .,bg = "transparent",scale = 2,
         units = "mm",
         dpi = 600,
         height = 0.618 * 84 * 2,
         width = 84 * 2,
         filename = "~/Desktop/AMR_carrier.eps")

p2 %>%
  export::graph2ppt(append = TRUE, width = 317.48031496 * 2  , height = 0.618 * 317.48031496 * 0.75 , paper = "A4",  scaling = 2,
                    file = out_pptx)
```

```
## Exported graph as ~/Desktop/plots_draft_chicken1.pptx
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

<img src="combine_chicken1_resistome_plots_draft_files/figure-html/unnamed-chunk-40-1.png" style="display: block; margin: auto;" />

```r
p1 %>%
  ggpubr::set_palette("npg") -> p1
```

```
## Scale for 'colour' is already present. Adding another scale for 'colour',
## which will replace the existing scale.
```

```
## Scale for 'fill' is already present. Adding another scale for 'fill', which
## will replace the existing scale.
```

```r
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
## Exported graph as ~/Desktop/plots_draft_chicken1.pptx
```



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
  # mutate(PathoFact_AMR_MGE_prediction_clean_2 = case_when(PathoFact_AMR_MGE_prediction_clean_2 %in% c("MAG", "chromosome")  ~ "MAG - chromosome", TRUE ~ PathoFact_AMR_MGE_prediction_clean_2)) %>% 
  mutate(PathoFact_AMR_MGE_prediction_clean_2 = fct_relevel(PathoFact_AMR_MGE_prediction_clean_2, c("MAG", "chromosome", "plasmid", "phage", "ambiguous", "unclassified"))) %>% 
  as.matrix() -> tax_table(ps_tmp)

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
  coord_flip() +  theme_light() +
  scale_y_continuous(labels=percent) + 
  facet_wrap(Model ~ Antibiotic + Antibiotic_mg.mL + Treatment2, scales = "free", ncol = 6, nrow = 4) -> p2
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
## Short values detected in phyloseq tax_table (nchar<4) :
## Consider using tax_fix() to make taxa uniquely identifiable
```

```
## Warning in ps_counts(data = data): otu_table of counts is NOT available!
## Available otu_table contains non-zero values that are less than 1
```

```r
p_carrier_all <- p2
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

<img src="combine_chicken1_resistome_plots_draft_files/figure-html/unnamed-chunk-41-1.png" style="display: block; margin: auto;" />

```r
p2 %>% 
  ggsave(plot = .,bg = "transparent",scale = 2,
         units = "mm",
         dpi = 600,
         height = 0.618 * 84 * 2,
         width = 84 * 2,
         filename = "~/Desktop/AMR_carrier_3.eps")

p2 %>%
  export::graph2ppt(append = TRUE, width = 317.48031496 * 2  , height = 0.618 * 317.48031496 * 0.75 , paper = "A4",  scaling = 2,
                    file = out_pptx)
```

```
## Exported graph as ~/Desktop/plots_draft_chicken1.pptx
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
## Exported graph as ~/Desktop/plots_draft_chicken1.pptx
```

```
## [1] "glycopeptide antibiotic"
```

```
## Exported graph as ~/Desktop/plots_draft_chicken1.pptx
```

```
## [1] "lincosamide antibiotic"
```

```
## Exported graph as ~/Desktop/plots_draft_chicken1.pptx
```

```
## [1] "tetracycline antibiotic"
```

```
## Exported graph as ~/Desktop/plots_draft_chicken1.pptx
```

```
## [1] "multidrug"
```

```
## Exported graph as ~/Desktop/plots_draft_chicken1.pptx
```

```
## [1] "macrolide antibiotic"
```

```
## Exported graph as ~/Desktop/plots_draft_chicken1.pptx
```

```
## [1] "aminoglycoside antibiotic"
```

```
## Exported graph as ~/Desktop/plots_draft_chicken1.pptx
```

```
## [1] "fluoroquinolone antibiotic"
```

```
## Exported graph as ~/Desktop/plots_draft_chicken1.pptx
```

```
## [1] "penam"
```

```
## Exported graph as ~/Desktop/plots_draft_chicken1.pptx
```

```
## [1] "fosfomycin"
```

```
## Exported graph as ~/Desktop/plots_draft_chicken1.pptx
```

```
## [1] "phenicol antibiotic"
```

```
## Exported graph as ~/Desktop/plots_draft_chicken1.pptx
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
## Exported graph as ~/Desktop/plots_draft_chicken1.pptx
```

```
## [1] "aminocoumarin antibiotic"
```

```
## Exported graph as ~/Desktop/plots_draft_chicken1.pptx
```

```
## [1] "streptogramin antibiotic"
```

```
## Exported graph as ~/Desktop/plots_draft_chicken1.pptx
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
## Exported graph as ~/Desktop/plots_draft_chicken1.pptx
```

```
## [1] "diaminopyrimidine antibiotic"
```

```
## Exported graph as ~/Desktop/plots_draft_chicken1.pptx
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
## Exported graph as ~/Desktop/plots_draft_chicken1.pptx
```

```
## [1] "rifamycin antibiotic"
```

```
## Exported graph as ~/Desktop/plots_draft_chicken1.pptx
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
## Exported graph as ~/Desktop/plots_draft_chicken1.pptx
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
## Exported graph as ~/Desktop/plots_draft_chicken1.pptx
```

```
## [1] "mupirocin"
```

```
## Exported graph as ~/Desktop/plots_draft_chicken1.pptx
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
## Exported graph as ~/Desktop/plots_draft_chicken1.pptx
```

```
## [1] "cephalosporin"
```

```
## Exported graph as ~/Desktop/plots_draft_chicken1.pptx
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
## Exported graph as ~/Desktop/plots_draft_chicken1.pptx
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
## Exported graph as ~/Desktop/plots_draft_chicken1.pptx
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
## Exported graph as ~/Desktop/plots_draft_chicken1.pptx
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
## Exported graph as ~/Desktop/plots_draft_chicken1.pptx
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
## Exported graph as ~/Desktop/plots_draft_chicken1.pptx
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
## Exported graph as ~/Desktop/plots_draft_chicken1.pptx
```

Vancomycin is a glycopeptide antibiotic used in the prophylaxis and treatment of infections caused by Gram-positive bacteria. Vancomycin inhibits the synthesis of peptidoglycan, the major component of the cell wall of gram-positive bacteria. Its mechanism of action is unusual in that it acts by binding precursors of peptidoglycan, rather than by interacting with an enzyme. -> glycopeptide antibiotic


```r
results$`glycopeptide antibiotic`
```

<img src="combine_chicken1_resistome_plots_draft_files/figure-html/unnamed-chunk-43-1.png" style="display: block; margin: auto;" />

```r
# results$`glycopeptide antibiotic` %>% 
#   export::graph2ppt(append = TRUE, width = 8, height = 4,
#                     file = out_pptx)
```

Cefotaxime is a semisynthetic cephalosporin taken parenterally. It is resistant to most beta-lactamases and active against Gram-negative rods and cocci due to its aminothiazoyl and methoximino functional groups. -> cephalosporin


```r
results$cephalosporin 
```

<img src="combine_chicken1_resistome_plots_draft_files/figure-html/unnamed-chunk-44-1.png" style="display: block; margin: auto;" />

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

<img src="combine_chicken1_resistome_plots_draft_files/figure-html/unnamed-chunk-51-1.png" style="display: block; margin: auto;" />

```r
p2 + facet_wrap(Model ~ Antibiotic + Antibiotic_mg.mL + Treatment2, scales = "free_y", ncol = 6, nrow = 4) -> p2_2

p2_2
```

<img src="combine_chicken1_resistome_plots_draft_files/figure-html/unnamed-chunk-51-2.png" style="display: block; margin: auto;" />

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
## Exported graph as ~/Desktop/plots_draft_chicken1.pptx
```

```r
p2_2 %>%
  export::graph2ppt(append = TRUE, width = 317.48031496 * 2  , height = 0.618 * 317.48031496 * 0.75 , paper = "A4",  scaling = 2,
                    file = out_pptx)
```

```
## Exported graph as ~/Desktop/plots_draft_chicken1.pptx
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

<img src="combine_chicken1_resistome_plots_draft_files/figure-html/unnamed-chunk-52-1.png" style="display: block; margin: auto;" />

```r
p2 + facet_wrap(Model ~ Antibiotic + Antibiotic_mg.mL + Treatment2, scales = "free_y", ncol = 6, nrow = 4) -> p2_2

p2_2
```

<img src="combine_chicken1_resistome_plots_draft_files/figure-html/unnamed-chunk-52-2.png" style="display: block; margin: auto;" />

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
## Exported graph as ~/Desktop/plots_draft_chicken1.pptx
```

```r
p2_2 %>%
  export::graph2ppt(append = TRUE, width = 317.48031496 * 2  , height = 0.618 * 317.48031496 * 0.75 , paper = "A4",  scaling = 2,
                    file = out_pptx)
```

```
## Exported graph as ~/Desktop/plots_draft_chicken1.pptx
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
  filter(genus %in% c("g__Enterococcus",
                      "g__Escherichia")) %>%
  mutate(Abundance = na_if(Abundance, 0)) %>%
  mutate(logTPM = log10(Abundance + 1))  %>% 
  # separate_rows(., Resistance_Mechanism,sep = '; ',convert = TRUE) %>% 
  # ggplot(aes(x = as.factor(Day_of_Treatment), y = Abundance, fill = fct_reorder(Strain, logTPM))) +   
  ggplot(aes(x = as.factor(Day_of_Treatment), y = Abundance, fill = GTDB_tax_CONCOCT_clean)) +#reorder(genus, Abundance, FUN=median, na.rm = TRUE))) +
  geom_bar( stat = "identity", colour="black", size=0.025) +
  facet_wrap(Model ~ Antibiotic + Antibiotic_mg.mL + Treatment2, scales = "free", ncol = 6, nrow = 4) +
  # facet_grid(Model ~ Treatment + Fermentation + Antibiotic_mg.mL , scales = "free", space = "free_x", drop = TRUE) +
  ggtitle("ARG Type Abundance per Sample") +
  xlab("Day (Treatment)")  +
  # ylab("Proportion") + 
  theme_light() +
  scale_y_continuous(labels=percent) + 
  # scale_fill_manual(name = "", values = col_strains) +
  theme(legend.position = "bottom") +
  ggpubr::rotate_x_text(60) +
  theme(legend.position = ("none")) -> p_mec

# p_mec %>% 
#   ggpubr::set_palette("jco") -> p_mec

p_mec
```

```
## Warning: Removed 46 rows containing missing values (position_stack).
```

<img src="combine_chicken1_resistome_plots_draft_files/figure-html/unnamed-chunk-56-1.png" style="display: block; margin: auto;" />

```r
# p_mec %>%
#   export::graph2ppt(append = TRUE, width = 12, height = 8,
#                     file = out_pptx)
```




```r
ps_AMR_all %>%
  physeq_sel_tax_table(c("PathoFact_AMR_AMR_sub_class", "PathoFact_AMR_AMR_sub_class_multi", "PathoFact_AMR_MGE_prediction_clean_2", "ANVIO_CONCOCT_HQ_clean", "CONCOCT_bin",
                         "GTDB_tax_CONCOCT_clean" ,"PathoFact_AMR_ARG")) %>% 
  prune_samples(samples, .) %>% 
  physeq_add_metadata(metadata_mapping,
                      sample_column = "sample_id") %>% 
  # subset_taxa(!is.na("GTDB_tax_CONCOCT_clean")) %>% 
  # physeq_sel_tax_table(c("PathoFact_AMR_ARG"))  %>% 
  # speedyseq::tax_glom(., taxrank = "GTDB_tax_CONCOCT_clean")  %>% 
  filter_taxa(function(x) sum(x > 0) > 0, TRUE) -> ps_tmp


ps_tmp %>%
  tax_table() %>% 
  data.frame() %>% 
  rownames_to_column('GN_id') %>%
  mutate(PathoFact_AMR_MGE_prediction_clean_2 = case_when(grepl("ambiguous", PathoFact_AMR_MGE_prediction_clean_2) ~ "ambiguous", TRUE ~ PathoFact_AMR_MGE_prediction_clean_2)) %>% 
  # mutate(PathoFact_AMR_MGE_prediction_clean_2 = case_when(PathoFact_AMR_MGE_prediction_clean_2 %in% c("MAG", "chromosome")  ~ "MAG - chromosome", TRUE ~ PathoFact_AMR_MGE_prediction_clean_2)) %>% 
  # mutate(PathoFact_AMR_MGE_prediction_clean_2 = fct_relevel(PathoFact_AMR_MGE_prediction_clean_2, c("MAG", "chromosome", "plasmid", "phage", "ambiguous", "unclassified"))) %>% 
  mutate(GN = paste(GN_id, PathoFact_AMR_MGE_prediction_clean_2, PathoFact_AMR_ARG, sep  = "_")) %>% 
  separate(GTDB_tax_CONCOCT_clean, into = c("domain", "phylum", "class", "order", "family", "genus", "species"), sep  = ";" , remove = FALSE) %>% 
  mutate(binid_phylum_species = paste(CONCOCT_bin, phylum, species, sep = ";")) %>% 
  mutate(MAG_plamid = case_when(PathoFact_AMR_MGE_prediction_clean_2 %in% c("plasmid") ~ PathoFact_AMR_MGE_prediction_clean_2, TRUE ~ binid_phylum_species)) %>%
  # mutate(MAG_plamid_gn = paste0(GN_id, MAG_plamid)) %>%
  column_to_rownames("GN_id") %>% 
  as.matrix() -> tax_table(ps_tmp)
```

```
## Warning: Expected 7 pieces. Missing pieces filled with `NA` in 28 rows [283,
## 468, 518, 528, 534, 543, 678, 704, 796, 839, 1023, 1237, 1312, 1440, 1530, 1550,
## 1635, 1728, 1787, 1810, ...].
```

```r
ps_tmp %>% 
  speedyseq::tax_glom("MAG_plamid") %>%
  physeq_sel_tax_table("MAG_plamid") %>%  
  transform_sample_counts(function(x) x/sum(x)) %>% 
  # subset_taxa(PathoFact_AMR_AMR_sub_class_multi %!in% c("Unclassified", "multidrug (unclassified)")) %>% 
  # ps_mutate(Day_of_Treatment = as.numeric(Day_of_Treatment)) %>% 
  ps_arrange(-Day_of_Treatment) %>% 
  comp_barplot(
    tax_transform_for_plot = "identity",
    tax_level = "MAG_plamid",
    label = "Day_of_Treatment", # name an alternative variable to label axes
    n_taxa = 20, # give more taxa unique colours
    merge_other = FALSE, # split the "other" category to display alpha diversity
    bar_width = 0.9, # reduce the bar width to 70% of one row
    bar_outline_colour = "grey5",
    sample_order = "default") +#,# is the default (use NA to remove outlines)
  # facet_by = c("Model", "Antibiotic", "Antibiotic_mg.mL", "Treatment2")) +
  # group_by = c("Model", "Antibiotic", "Antibiotic_mg.mL", "Treatment2"))
  coord_flip() + theme_light() + 
  scale_y_continuous(labels=percent) + 
  facet_wrap(Model ~ Antibiotic + Antibiotic_mg.mL + Treatment2, scales = "free", ncol = 6, nrow = 4) +
  guides(fill=guide_legend(ncol = 1)) -> p2
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
p2 +
  theme(legend.position = "none") -> p2_noleg

p2_noleg
```

<img src="combine_chicken1_resistome_plots_draft_files/figure-html/unnamed-chunk-58-1.png" style="display: block; margin: auto;" />

```r
p2_noleg %>%
  export::graph2ppt(append = TRUE, width = 317.48031496 * 2  , height = 0.618 * 317.48031496 *0.75, paper = "A4",  scaling = 2,
                    file = out_pptx)
```

```
## Exported graph as ~/Desktop/plots_draft_chicken1.pptx
```


```r
p2 %>% 
  ggpubr::get_legend() %>% 
  ggpubr::as_ggplot() -> p2_leg

p2_leg
```

<img src="combine_chicken1_resistome_plots_draft_files/figure-html/unnamed-chunk-59-1.png" style="display: block; margin: auto;" />

```r
p2_leg %>%
  export::graph2ppt(append = TRUE, width = 317.48031496 *2  , height = 0.618 * 317.48031496 *2 , paper = "A4",  scaling = 2,
                    file = out_pptx)
```

```
## Exported graph as ~/Desktop/plots_draft_chicken1.pptx
```
Lokking only at only E coli and E facium MAGs: s_Escherichia coli s_Enterococcus faecalis



```r
ps_AMR_all %>%
  physeq_sel_tax_table(c("PathoFact_AMR_AMR_sub_class", "PathoFact_AMR_AMR_sub_class_multi", "PathoFact_AMR_MGE_prediction_clean_2", "ANVIO_CONCOCT_HQ_clean", "CONCOCT_bin",
                         "GTDB_tax_CONCOCT_clean" ,"PathoFact_AMR_ARG")) %>% 
  prune_samples(samples, .) %>% 
  physeq_add_metadata(metadata_mapping,
                      sample_column = "sample_id") %>% 
  # subset_taxa(!is.na("GTDB_tax_CONCOCT_clean")) %>% 
  # physeq_sel_tax_table(c("PathoFact_AMR_ARG"))  %>% 
  # speedyseq::tax_glom(., taxrank = "GTDB_tax_CONCOCT_clean")  %>% 
  filter_taxa(function(x) sum(x > 0) > 0, TRUE) -> ps_tmp


ps_tmp %>%
  tax_table() %>% 
  data.frame() %>% 
  rownames_to_column('GN_id') %>%
  mutate(PathoFact_AMR_MGE_prediction_clean_2 = case_when(grepl("ambiguous", PathoFact_AMR_MGE_prediction_clean_2) ~ "ambiguous", TRUE ~ PathoFact_AMR_MGE_prediction_clean_2)) %>% 
  # mutate(PathoFact_AMR_MGE_prediction_clean_2 = case_when(PathoFact_AMR_MGE_prediction_clean_2 %in% c("MAG", "chromosome")  ~ "MAG - chromosome", TRUE ~ PathoFact_AMR_MGE_prediction_clean_2)) %>% 
  # mutate(PathoFact_AMR_MGE_prediction_clean_2 = fct_relevel(PathoFact_AMR_MGE_prediction_clean_2, c("MAG", "chromosome", "plasmid", "phage", "ambiguous", "unclassified"))) %>% 
  mutate(GN = paste(GN_id, PathoFact_AMR_MGE_prediction_clean_2, PathoFact_AMR_ARG, sep  = "_")) %>% 
  separate(GTDB_tax_CONCOCT_clean, into = c("domain", "phylum", "class", "order", "family", "genus", "species"), sep  = ";" , remove = FALSE) %>% 
  mutate(binid_phylum_species = paste(CONCOCT_bin, phylum, species, sep = ";")) %>% 
  mutate(MAG_plamid = case_when(PathoFact_AMR_MGE_prediction_clean_2 %in% c("plasmid") ~ PathoFact_AMR_MGE_prediction_clean_2, TRUE ~ binid_phylum_species)) %>%
  # mutate(MAG_plamid_gn = paste0(GN_id, MAG_plamid)) %>%
  column_to_rownames("GN_id") %>% 
  as.matrix() -> tax_table(ps_tmp)
```

```
## Warning: Expected 7 pieces. Missing pieces filled with `NA` in 28 rows [283,
## 468, 518, 528, 534, 543, 678, 704, 796, 839, 1023, 1237, 1312, 1440, 1530, 1550,
## 1635, 1728, 1787, 1810, ...].
```

```r
ps_tmp %>% 
  physeq_sel_tax_table("MAG_plamid") %>%  
  speedyseq::tax_glom("MAG_plamid") %>%
  transform_sample_counts(function(x) x/sum(x)) %>% 
  subset_taxa(grepl("coli|plasmid|Enterococcus_B faecium",  MAG_plamid)) %>%  # %!in% c("Unclassified", "multidrug (unclassified)")) %>%
  # ps_mutate(Day_of_Treatment = as.numeric(Day_of_Treatment)) %>% 
  ps_arrange(-Day_of_Treatment) %>% 
  comp_barplot(
    tax_transform_for_plot = "identity",
    tax_level = "MAG_plamid",
    label = "Day_of_Treatment", # name an alternative variable to label axes
    n_taxa = 30, # give more taxa unique colours
    merge_other = FALSE, # split the "other" category to display alpha diversity
    bar_width = 0.9, # reduce the bar width to 70% of one row
    bar_outline_colour = "grey5",
    sample_order = "default") +#,# is the default (use NA to remove outlines)
  # facet_by = c("Model", "Antibiotic", "Antibiotic_mg.mL", "Treatment2")) +
  # group_by = c("Model", "Antibiotic", "Antibiotic_mg.mL", "Treatment2"))
  coord_flip() + theme_light() + 
  scale_y_continuous(labels=percent) + 
  facet_wrap(Model ~ Antibiotic + Antibiotic_mg.mL + Treatment2, scales = "free_y", ncol = 6, nrow = 4) +
  guides(fill=guide_legend(ncol = 1)) -> p2
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
p2
```

<img src="combine_chicken1_resistome_plots_draft_files/figure-html/unnamed-chunk-61-1.png" style="display: block; margin: auto;" />

```r
p2 +
  theme(legend.position = "none") -> p2_noleg

p2_noleg
```

<img src="combine_chicken1_resistome_plots_draft_files/figure-html/unnamed-chunk-61-2.png" style="display: block; margin: auto;" />

```r
p2_noleg %>%
  export::graph2ppt(append = TRUE, width = 317.48031496 * 2  , height = 0.618 * 317.48031496 * 0.75 , paper = "A4",  scaling = 2,
                    file = out_pptx)
```

```
## Exported graph as ~/Desktop/plots_draft_chicken1.pptx
```

Prepare ANVIO stuff:

Bin_name
GTDB
phylum
Genera
Species

number AMRG



```r
ps_AMR_all %>%
  physeq_sel_tax_table(c("PathoFact_AMR_AMR_sub_class", "PathoFact_AMR_AMR_sub_class_multi", "PathoFact_AMR_MGE_prediction_clean_2", "ANVIO_CONCOCT_HQ_clean", "CONCOCT_bin",
                         "GTDB_tax_CONCOCT_clean" ,"PathoFact_AMR_ARG")) %>% 
  prune_samples(samples, .) %>% 
  physeq_add_metadata(metadata_mapping,
                      sample_column = "sample_id") %>% 
  # subset_taxa(!is.na("GTDB_tax_CONCOCT_clean")) %>% 
  # physeq_sel_tax_table(c("PathoFact_AMR_ARG"))  %>% 
  # speedyseq::tax_glom(., taxrank = "GTDB_tax_CONCOCT_clean")  %>% 
  filter_taxa(function(x) sum(x > 0) > 0, TRUE) -> ps_tmp


ps_tmp %>%
  tax_table() %>% 
  data.frame() %>% 
  rownames_to_column('GN_id') %>%
  mutate(PathoFact_AMR_MGE_prediction_clean_2 = case_when(grepl("ambiguous", PathoFact_AMR_MGE_prediction_clean_2) ~ "ambiguous", TRUE ~ PathoFact_AMR_MGE_prediction_clean_2)) %>% 
  # mutate(PathoFact_AMR_MGE_prediction_clean_2 = case_when(PathoFact_AMR_MGE_prediction_clean_2 %in% c("MAG", "chromosome")  ~ "MAG - chromosome", TRUE ~ PathoFact_AMR_MGE_prediction_clean_2)) %>% 
  # mutate(PathoFact_AMR_MGE_prediction_clean_2 = fct_relevel(PathoFact_AMR_MGE_prediction_clean_2, c("MAG", "chromosome", "plasmid", "phage", "ambiguous", "unclassified"))) %>% 
  mutate(GN = paste(GN_id, PathoFact_AMR_MGE_prediction_clean_2, PathoFact_AMR_ARG, sep  = "_")) %>% 
  separate(GTDB_tax_CONCOCT_clean, into = c("domain", "phylum", "class", "order", "family", "genus", "species"), sep  = ";" , remove = FALSE) %>% 
  mutate(binid_phylum_species = paste(CONCOCT_bin, phylum, species, sep = ";")) %>% 
  mutate(MAG_plamid = case_when(PathoFact_AMR_MGE_prediction_clean_2 %in% c("plasmid") ~ PathoFact_AMR_MGE_prediction_clean_2, TRUE ~ binid_phylum_species)) %>%
  # mutate(MAG_plamid_gn = paste0(GN_id, MAG_plamid)) %>%
  column_to_rownames("GN_id") %>% 
  as.matrix() -> tax_table(ps_tmp)
```

```
## Warning: Expected 7 pieces. Missing pieces filled with `NA` in 28 rows [283,
## 468, 518, 528, 534, 543, 678, 704, 796, 839, 1023, 1237, 1312, 1440, 1530, 1550,
## 1635, 1728, 1787, 1810, ...].
```

```r
ps_tmp %>%
  tax_table() %>% 
  data.frame() %>% 
  rownames_to_column('GN_id') %>%
  mutate(PathoFact_AMR_MGE_prediction_clean_2 = case_when(grepl("ambiguous", PathoFact_AMR_MGE_prediction_clean_2) ~ "ambiguous", TRUE ~ PathoFact_AMR_MGE_prediction_clean_2)) %>% 
  # mutate(PathoFact_AMR_MGE_prediction_clean_2 = case_when(PathoFact_AMR_MGE_prediction_clean_2 %in% c("MAG", "chromosome")  ~ "MAG - chromosome", TRUE ~ PathoFact_AMR_MGE_prediction_clean_2)) %>% 
  # mutate(PathoFact_AMR_MGE_prediction_clean_2 = fct_relevel(PathoFact_AMR_MGE_prediction_clean_2, c("MAG", "chromosome", "plasmid", "phage", "ambiguous", "unclassified"))) %>% 
  mutate(GN = paste(GN_id, PathoFact_AMR_MGE_prediction_clean_2, PathoFact_AMR_ARG, sep  = "_")) %>% 
  separate(GTDB_tax_CONCOCT_clean, into = c("domain", "phylum", "class", "order", "family", "genus", "species"), sep  = ";" , remove = FALSE) %>% 
  mutate(binid_phylum_species = paste(CONCOCT_bin, phylum, species, sep = ";")) %>% 
  mutate(MAG_plamid = case_when(PathoFact_AMR_MGE_prediction_clean_2 %in% c("plasmid") ~ PathoFact_AMR_MGE_prediction_clean_2, TRUE ~ binid_phylum_species)) %>%
  # mutate(MAG_plamid_gn = paste0(GN_id, MAG_plamid)) %>%
  column_to_rownames("GN_id") %>% 
  as.matrix() -> tax_table(ps_tmp)
```

```
## Warning: Expected 7 pieces. Missing pieces filled with `NA` in 28 rows [283,
## 468, 518, 528, 534, 543, 678, 704, 796, 839, 1023, 1237, 1312, 1440, 1530, 1550,
## 1635, 1728, 1787, 1810, ...].
```

```r
# tax_table(ps_tmp) %>% 
#   data.frame() %>% 
#   write_tsv("~/Desktop/MAG_taxa_ARG_h1_2_work_inprogress_first.tsv")
```


Count ARG:


```r
tax_table(ps_tmp) %>% 
  data.frame() %>% 
  rownames_to_column('rowname') %>% 
  filter(! grepl('nclassified', PathoFact_AMR_AMR_sub_class_multi)) %>% 
  separate(rowname,into = c("model", "mega", "contigid", "start", "stop"), sep = "_") %>% 
  filter(!is.na(CONCOCT_bin)) %>% 
  # filter(ANVIO_CONCOCT_HQ_clean == "HQ") %>% 
  select(GN, model, binid_phylum_species, PathoFact_AMR_ARG, species, genus, phylum, CONCOCT_bin, PathoFact_AMR_AMR_sub_class, PathoFact_AMR_AMR_sub_class_multi, GTDB_tax_CONCOCT_clean) %>% 
  select(model, CONCOCT_bin, binid_phylum_species, PathoFact_AMR_ARG, GTDB_tax_CONCOCT_clean) %>% 
  # separate_rows(PathoFact_AMR_AMR_sub_class, sep = ";") %>% 
  group_by(model, CONCOCT_bin, binid_phylum_species, GTDB_tax_CONCOCT_clean) %>% 
  summarise(num_ARG = length(PathoFact_AMR_ARG)) -> MAGs_count_ARG
```

```
## Warning: Expected 5 pieces. Missing pieces filled with `NA` in 1929 rows [1, 2,
## 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, ...].
```

```
## `summarise()` has grouped output by 'model', 'CONCOCT_bin',
## 'binid_phylum_species'. You can override using the `.groups` argument.
```

```r
MAGs_count_ARG
```

```
## # A tibble: 177 × 5
## # Groups:   model, CONCOCT_bin, binid_phylum_species [176]
##    model CONCOCT_bin binid_phylum_species                        GTDB_…¹ num_ARG
##    <chr> <chr>       <chr>                                       <chr>     <int>
##  1 c1    Bin_1_1     Bin_1_1;p__Actinobacteriota;s__Rhodococcus… d__Bac…      12
##  2 c1    Bin_1_2     Bin_1_2;p__Firmicutes_A;s__                 d__Bac…       4
##  3 c1    Bin_10_1    Bin_10_1;p__Firmicutes_A;s__                d__Bac…       4
##  4 c1    Bin_10_2    Bin_10_2;p__Firmicutes_A;s__                d__Bac…       4
##  5 c1    Bin_11_1_1  Bin_11_1_1;p__Firmicutes_A;s__Sellimonas a… d__Bac…       8
##  6 c1    Bin_11_10_1 Bin_11_10_1;p__Firmicutes_A;s__Mediterrane… d__Bac…       5
##  7 c1    Bin_11_11   Bin_11_11;p__Firmicutes_A;s__Enterocloster… d__Bac…       6
##  8 c1    Bin_11_12   Bin_11_12;p__Firmicutes_A;s__Lachnoclostri… d__Bac…       5
##  9 c1    Bin_11_13   Bin_11_13;p__Firmicutes_A;s__Ruminococcus_… d__Bac…       4
## 10 c1    Bin_11_14   Bin_11_14;p__Firmicutes_A;s__Blautia pulli… d__Bac…       5
## # … with 167 more rows, and abbreviated variable name ¹​GTDB_tax_CONCOCT_clean
## # ℹ Use `print(n = ...)` to see more rows
```

```r
MAGs_count_ARG %>%  filter(model == "h1") %>%  add_count(CONCOCT_bin) %>%  arrange(desc(n))
```

```
## # A tibble: 0 × 6
## # Groups:   model, CONCOCT_bin, binid_phylum_species [0]
## # … with 6 variables: model <chr>, CONCOCT_bin <chr>,
## #   binid_phylum_species <chr>, GTDB_tax_CONCOCT_clean <chr>, num_ARG <int>,
## #   n <int>
## # ℹ Use `colnames()` to see all variable names
```

```r
# MAGs_count_ARG %>% 
#   DT::datatable()
```


drug class:


```r
# tax_table(ps_tmp) %>% 
#   data.frame() %>% 
#   rownames_to_column('rowname') %>% 
#   separate(rowname,into = c("model", "mega", "contigid", "start", "stop"), sep = "_") %>% 
#     filter(!is.na(CONCOCT_bin)) %>% 
#   filter(ANVIO_CONCOCT_HQ_clean == "HQ") %>% 
#     filter(! grepl('nclassified', PathoFact_AMR_AMR_sub_class_multi)) %>% 
#   select(GN, model, binid_phylum_species, PathoFact_AMR_ARG, species, genus, phylum, CONCOCT_bin, PathoFact_AMR_AMR_sub_class, PathoFact_AMR_AMR_sub_class_multi, GTDB_tax_CONCOCT_clean) %>% 
#   select(model, CONCOCT_bin, binid_phylum_species, CONCOCT_bin, PathoFact_AMR_AMR_sub_class_multi) %>% 
#   mutate(mode_binid_phylum_species = paste0(model, "_", binid_phylum_species)) %>% 
#   mutate(unic_id = paste0(mode_binid_phylum_species, PathoFact_AMR_AMR_sub_class_multi)) %>% 
#   # separate_rows(PathoFact_AMR_AMR_sub_class, sep = ";") %>% 
#   distinct(unic_id, .keep_all = TRUE) %>% 
#   select(-unic_id) %>% 
#   # group_by(model, CONCOCT_bin, binid_phylum_species, PathoFact_AMR_AMR_sub_class_multi) %>% 
#     # summarise(num_ARG = length(PathoFact_AMR_AMR_sub_class_multi)) %>% 
#     pivot_wider(names_from =  PathoFact_AMR_AMR_sub_class_multi,
#               values_from = mode_binid_phylum_species ) 
```



```r
# tax_table(ps_tmp) %>% 
#   data.frame() %>% 
#   rownames_to_column('rowname') %>% 
#   separate(rowname,into = c("model", "mega", "contigid", "start", "stop"), sep = "_") %>% 
#   filter(!is.na(CONCOCT_bin)) %>% 
#   filter(ANVIO_CONCOCT_HQ_clean == "HQ") %>% 
#   filter(! grepl('nclassified', PathoFact_AMR_AMR_sub_class_multi)) %>% 
#   separate_rows(PathoFact_AMR_AMR_sub_class, sep = ";", convert = TRUE) %>%
#   select(GN, model, binid_phylum_species, PathoFact_AMR_ARG, species, genus, phylum, CONCOCT_bin, PathoFact_AMR_AMR_sub_class, PathoFact_AMR_AMR_sub_class_multi, GTDB_tax_CONCOCT_clean) %>% 
#   select(model, CONCOCT_bin, binid_phylum_species, CONCOCT_bin, PathoFact_AMR_AMR_sub_class) %>% 
#   mutate(mode_binid_phylum_species = paste0(model, "_", binid_phylum_species)) %>% 
#   # separate_rows(PathoFact_AMR_AMR_sub_class, sep = ";") %>% 
#   group_by(model, CONCOCT_bin, binid_phylum_species, PathoFact_AMR_AMR_sub_class) %>%
#   summarise(num_ARG = length(PathoFact_AMR_AMR_sub_class)) %>%
#   pivot_wider(names_from =  PathoFact_AMR_AMR_sub_class,
#               values_from = num_ARG ) -> MAGs_count_class
# 
# MAGs_count_class
# 
# MAGs_count_class %>%  filter(model == "h1") %>%  add_count(CONCOCT_bin) %>%  arrange(desc(n))
```



```r
# tax_table(ps_tmp) %>% 
#   data.frame() %>% 
#   rownames_to_column('rowname') %>% 
#   separate(rowname,into = c("model", "mega", "contigid", "start", "stop"), sep = "_") %>% 
#   filter(!is.na(CONCOCT_bin)) %>% 
#   # filter(ANVIO_CONCOCT_HQ_clean == "HQ") %>% 
#   filter(! grepl('nclassified', PathoFact_AMR_AMR_sub_class_multi)) %>% 
#   select(GN, model, binid_phylum_species, PathoFact_AMR_ARG, species, genus, phylum, CONCOCT_bin, PathoFact_AMR_AMR_sub_class_multi, GTDB_tax_CONCOCT_clean) %>% 
#   select(model, CONCOCT_bin, binid_phylum_species, CONCOCT_bin, PathoFact_AMR_AMR_sub_class_multi) %>% 
#   mutate(mode_binid_phylum_species = paste0(model, "_", binid_phylum_species)) %>% 
#   # separate_rows(PathoFact_AMR_AMR_sub_class, sep = ";") %>% 
#   group_by(model, CONCOCT_bin, binid_phylum_species, PathoFact_AMR_AMR_sub_class_multi) %>%
#   summarise(num_ARG = length(PathoFact_AMR_AMR_sub_class_multi)) %>%
#   pivot_wider(names_from =  PathoFact_AMR_AMR_sub_class_multi,
#               values_from = num_ARG ) -> MAGs_count_class_multi
# 
# MAGs_count_class_multi
# 
# 
# MAGs_count_class_multi %>%  filter(model == "h1") %>%  add_count(CONCOCT_bin) %>%  arrange(desc(n))

# MAGs_count_class_multi %>%
#   DT::datatable()
```


Combined results:


```r
# MAGs_count_ARG %>% 
#   # left_join(.,
#   #           MAGs_count_class,
#   #           by = c("binid_phylum_species" = "binid_phylum_species")) %>% 
#   left_join(.,
#             MAGs_count_class_multi,
#             by = c("binid_phylum_species" = "binid_phylum_species"),
#             suffix = c("", "_multi")) %>% 
#   # select(model,CONCOCT_bin, binid_phylum_species,num_ARG) %>% 
#   # rename_at(vars(everything()), ~str_replace_all(., "\\s+", "")) %>% 
#   # replace(is.na(.), 0) %>% 
#   # left_join(ps_tmp %>%
#   #             tax_table() %>%
#   #             data.frame() %>% 
#   #           select(binid_phylum_species, GTDB_tax_CONCOCT_clean, domain:species) %>%  distinct(binid_phylum_species, .keep_all = TRUE),
#   #                       by = c("binid_phylum_species" = "binid_phylum_species")) %>%
#   select(CONCOCT_bin, everything()) %>% 
#   filter(model == "h1" & model_multi == "h1") %>% 
#   separate(GTDB_tax_CONCOCT_clean, into = c("domain", "phylum", "class", "order", "family", "genus", "species"), sep  = ";" , remove = FALSE) %>% 
#   write_tsv("~/Desktop/MAG_taxa_ARG_h1_2_work_inprogress.tsv")
```

sample_data_stuff:

sample carrier stuff:


```r
# p_carrier_all$data %>% 
#   select(OTU, Sample, Abundance) %>% 
#   dplyr::group_by(Sample) %>% 
#   pivot_wider(values_from = Abundance,
#               names_from = OTU) -> sample_carrier_table
```


Class_multi:


```r
# plot_AMR_multi$data %>% 
#   filter(! grepl('nclassified', OTU)) %>% 
#   select(OTU, Sample, Abundance) %>% 
#   dplyr::group_by(Sample) %>% 
#   pivot_wider(values_from = Abundance,
#               names_from = OTU) -> sample_AMR_table
```




```r
ps_AMR %>% 
  physeq_add_metadata(metadata_mapping,
                      sample_column = "sample_id") %>% 
  sample_data() %>% 
  data.frame() %>% 
  filter(Model == "Human") #%>% 
```

```
##       sample_id Description2 Experiment Reactor     Treatment Day_of_Connection
## D.1         D.1        DONOR      Cecum   DONOR         DONOR                NA
## HC.1       HC.1         <NA> Continuous     TR1   CTX+HV292.1                 8
## HC.2       HC.2         <NA> Continuous     TR2           CTX                 8
## HC.3       HC.3         <NA> Continuous     TR3   CTX+HV292.1                 8
## HC.4       HC.4         <NA> Continuous     TR4           CTX                 8
## HC.5       HC.5         <NA> Continuous     TR5       HV292.1                 8
## HC.6       HC.6         <NA> Continuous      CR     UNTREATED                 8
## HC.7       HC.7         <NA> Continuous     TR1   CTX+HV292.1                10
## HC.8       HC.8         <NA> Continuous     TR2           CTX                10
## HC.9       HC.9         <NA> Continuous     TR3   CTX+HV292.1                10
## HC.10     HC.10         <NA> Continuous     TR4           CTX                10
## HC.11     HC.11         <NA> Continuous     TR5       HV292.1                10
## HC.12     HC.12         <NA> Continuous      CR     UNTREATED                10
## HC.13     HC.13         <NA> Continuous     TR1   CTX+HV292.1                12
## HC.14     HC.14         <NA> Continuous     TR2           CTX                12
## HC.15     HC.15         <NA> Continuous     TR3   CTX+HV292.1                12
## HC.16     HC.16         <NA> Continuous     TR4           CTX                12
## HC.17     HC.17         <NA> Continuous     TR5       HV292.1                12
## HC.18     HC.18         <NA> Continuous      CR     UNTREATED                12
## HC.19     HC.19         <NA> Continuous     TR1   CTX+HV292.1                16
## HC.20     HC.20         <NA> Continuous     TR2           CTX                16
## HC.21     HC.21         <NA> Continuous     TR3   CTX+HV292.1                16
## HC.22     HC.22         <NA> Continuous     TR4           CTX                16
## HC.23     HC.23         <NA> Continuous     TR5       HV292.1                16
## HC.24     HC.24         <NA> Continuous      CR     UNTREATED                16
## HV.1       HV.1         <NA> Continuous     TR1 VAN+CCUG59168                28
## HV.2       HV.2         <NA> Continuous     TR2           VAN                28
## HV.3       HV.3         <NA> Continuous     TR3 VAN+CCUG59168                28
## HV.4       HV.4         <NA> Continuous     TR4           VAN                28
## HV.5       HV.5         <NA> Continuous     TR5     CCUG59168                28
## HV.6       HV.6         <NA> Continuous      CR     UNTREATED                28
## HV.7       HV.7         <NA> Continuous     TR1 VAN+CCUG59168                30
## HV.8       HV.8         <NA> Continuous     TR2           VAN                30
## HV.9       HV.9         <NA> Continuous     TR3 VAN+CCUG59168                30
## HV.10     HV.10         <NA> Continuous     TR4           VAN                30
## HV.11     HV.11         <NA> Continuous     TR5     CCUG59168                30
## HV.12     HV.12         <NA> Continuous      CR     UNTREATED                30
## HV.13     HV.13         <NA> Continuous     TR1 VAN+CCUG59168                32
## HV.14     HV.14         <NA> Continuous     TR2           VAN                32
## HV.15     HV.15         <NA> Continuous     TR3 VAN+CCUG59168                32
## HV.16     HV.16         <NA> Continuous     TR4           VAN                32
## HV.17     HV.17         <NA> Continuous     TR5     CCUG59168                32
## HV.18     HV.18         <NA> Continuous      CR     UNTREATED                32
## HV.19     HV.19         <NA> Continuous     TR1 VAN+CCUG59168                36
## HV.20     HV.20         <NA> Continuous     TR2           VAN                36
## HV.21     HV.21         <NA> Continuous     TR3 VAN+CCUG59168                36
## HV.22     HV.22         <NA> Continuous     TR4           VAN                36
## HV.23     HV.23         <NA> Continuous     TR5     CCUG59168                36
## HV.24     HV.24         <NA> Continuous      CR     UNTREATED                36
##       Day_of_Treatment Day_from_Inoculum  Enrichment Phase    Treatment2 Date
## D.1                 NA                NA NotEnriched  <NA>         DONOR <NA>
## HC.1                -1                76 NotEnriched  Stab    AB+E. coli <NA>
## HC.2                -1                76 NotEnriched  Stab            AB <NA>
## HC.3                -1                76 NotEnriched  Stab    AB+E. coli <NA>
## HC.4                -1                76 NotEnriched  Stab            AB <NA>
## HC.5                -1                76 NotEnriched  Stab       E. coli <NA>
## HC.6                -1                76 NotEnriched  Stab     UNTREATED <NA>
## HC.7                 1                78 NotEnriched Treat    AB+E. coli <NA>
## HC.8                 1                78 NotEnriched Treat            AB <NA>
## HC.9                 1                78 NotEnriched Treat    AB+E. coli <NA>
## HC.10                1                78 NotEnriched Treat            AB <NA>
## HC.11                1                78 NotEnriched Treat       E. coli <NA>
## HC.12                1                78 NotEnriched Treat     UNTREATED <NA>
## HC.13                3                80 NotEnriched Treat    AB+E. coli <NA>
## HC.14                3                80 NotEnriched Treat            AB <NA>
## HC.15                3                80 NotEnriched Treat    AB+E. coli <NA>
## HC.16                3                80 NotEnriched Treat            AB <NA>
## HC.17                3                80 NotEnriched Treat       E. coli <NA>
## HC.18                3                80 NotEnriched Treat     UNTREATED <NA>
## HC.19                7                84 NotEnriched Treat    AB+E. coli <NA>
## HC.20                7                84 NotEnriched Treat            AB <NA>
## HC.21                7                84 NotEnriched Treat    AB+E. coli <NA>
## HC.22                7                84 NotEnriched Treat            AB <NA>
## HC.23                7                84 NotEnriched Treat       E. coli <NA>
## HC.24                7                84 NotEnriched Treat     UNTREATED <NA>
## HV.1                -1                43 NotEnriched  Stab AB+E. faecium <NA>
## HV.2                -1                43 NotEnriched  Stab            AB <NA>
## HV.3                -1                43 NotEnriched  Stab AB+E. faecium <NA>
## HV.4                -1                43 NotEnriched  Stab            AB <NA>
## HV.5                -1                43 NotEnriched  Stab    E. faecium <NA>
## HV.6                -1                43 NotEnriched  Stab     UNTREATED <NA>
## HV.7                 1                45 NotEnriched Treat AB+E. faecium <NA>
## HV.8                 1                45 NotEnriched Treat            AB <NA>
## HV.9                 1                45 NotEnriched Treat AB+E. faecium <NA>
## HV.10                1                45 NotEnriched Treat            AB <NA>
## HV.11                1                45 NotEnriched Treat    E. faecium <NA>
## HV.12                1                45 NotEnriched Treat     UNTREATED <NA>
## HV.13                3                47 NotEnriched Treat AB+E. faecium <NA>
## HV.14                3                47 NotEnriched Treat            AB <NA>
## HV.15                3                47 NotEnriched Treat AB+E. faecium <NA>
## HV.16                3                47 NotEnriched Treat            AB <NA>
## HV.17                3                47 NotEnriched Treat    E. faecium <NA>
## HV.18                3                47 NotEnriched Treat     UNTREATED <NA>
## HV.19                7                51 NotEnriched Treat AB+E. faecium <NA>
## HV.20                7                51 NotEnriched Treat            AB <NA>
## HV.21                7                51 NotEnriched Treat AB+E. faecium <NA>
## HV.22                7                51 NotEnriched Treat            AB <NA>
## HV.23                7                51 NotEnriched Treat    E. faecium <NA>
## HV.24                7                51 NotEnriched Treat     UNTREATED <NA>
##       Paul Reactor_Treatment GeneCopyNumberperML HV292.1_Copy_Number_permL
## D.1     NA             DONOR            1.90e+10                        NA
## HC.1    NA   TR1_CTX+HV292.1            2.14e+11                        NA
## HC.2    NA           TR2_CTX            2.07e+11                        NA
## HC.3    NA   TR3_CTX+HV292.1            1.98e+11                        NA
## HC.4    NA           TR4_CTX            2.05e+11                        NA
## HC.5    NA       TR5_HV292.1            2.43e+11                        NA
## HC.6    NA      CR_UNTREATED            2.37e+11                        NA
## HC.7    NA   TR1_CTX+HV292.1            5.59e+11                   5520000
## HC.8    NA           TR2_CTX            2.66e+11                         0
## HC.9    NA   TR3_CTX+HV292.1            3.83e+11                  30487500
## HC.10   NA           TR4_CTX            3.22e+11                         0
## HC.11   NA       TR5_HV292.1            4.52e+11                  20600000
## HC.12   NA      CR_UNTREATED            2.39e+11                         0
## HC.13   NA   TR1_CTX+HV292.1            1.36e+11                         0
## HC.14   NA           TR2_CTX            2.85e+11                         0
## HC.15   NA   TR3_CTX+HV292.1            2.04e+11                         0
## HC.16   NA           TR4_CTX            2.39e+11                         0
## HC.17   NA       TR5_HV292.1            2.54e+11                         0
## HC.18   NA      CR_UNTREATED            1.26e+11                         0
## HC.19   NA   TR1_CTX+HV292.1            2.10e+11                         0
## HC.20   NA           TR2_CTX            1.02e+11                         0
## HC.21   NA   TR3_CTX+HV292.1            2.41e+11                         0
## HC.22   NA           TR4_CTX            2.42e+11                         0
## HC.23   NA       TR5_HV292.1            3.93e+11                         0
## HC.24   NA      CR_UNTREATED            2.60e+11                         0
## HV.1    NA TR1_VAN+CCUG59168            1.11e+11                        NA
## HV.2    NA           TR2_VAN            2.50e+11                        NA
## HV.3    NA TR3_VAN+CCUG59168            3.13e+11                        NA
## HV.4    NA           TR4_VAN            2.99e+10                        NA
## HV.5    NA     TR5_CCUG59168            1.92e+11                        NA
## HV.6    NA      CR_UNTREATED            2.58e+11                        NA
## HV.7    NA TR1_VAN+CCUG59168            3.03e+11                        NA
## HV.8    NA           TR2_VAN            7.45e+10                        NA
## HV.9    NA TR3_VAN+CCUG59168            1.21e+11                        NA
## HV.10   NA           TR4_VAN            1.32e+11                        NA
## HV.11   NA     TR5_CCUG59168            2.77e+11                        NA
## HV.12   NA      CR_UNTREATED            4.92e+10                        NA
## HV.13   NA TR1_VAN+CCUG59168            1.94e+11                        NA
## HV.14   NA           TR2_VAN            2.43e+11                        NA
## HV.15   NA TR3_VAN+CCUG59168            2.00e+11                        NA
## HV.16   NA           TR4_VAN            9.63e+10                        NA
## HV.17   NA     TR5_CCUG59168            2.00e+11                        NA
## HV.18   NA      CR_UNTREATED            3.02e+11                        NA
## HV.19   NA TR1_VAN+CCUG59168            4.92e+10                        NA
## HV.20   NA           TR2_VAN            3.21e+11                        NA
## HV.21   NA TR3_VAN+CCUG59168            2.91e+11                        NA
## HV.22   NA           TR4_VAN            1.17e+11                        NA
## HV.23   NA     TR5_CCUG59168            1.68e+11                        NA
## HV.24   NA      CR_UNTREATED            4.39e+11                        NA
##       CCUG59168_Copy_Number_permL CTX_Copy_Number_permL VAN_Copy_Number_permL
## D.1                            NA                    NA                    NA
## HC.1                           NA                    NA                    NA
## HC.2                           NA                    NA                    NA
## HC.3                           NA                    NA                    NA
## HC.4                           NA                    NA                    NA
## HC.5                           NA                    NA                    NA
## HC.6                           NA                    NA                    NA
## HC.7                           NA              11025000                    NA
## HC.8                           NA                     0                    NA
## HC.9                           NA              53737500                    NA
## HC.10                          NA                     0                    NA
## HC.11                          NA              41600000                    NA
## HC.12                          NA                     0                    NA
## HC.13                          NA               2160000                    NA
## HC.14                          NA                     0                    NA
## HC.15                          NA               3550000                    NA
## HC.16                          NA                     0                    NA
## HC.17                          NA                210000                    NA
## HC.18                          NA                     0                    NA
## HC.19                          NA                     0                    NA
## HC.20                          NA                     0                    NA
## HC.21                          NA               1117500                    NA
## HC.22                          NA                     0                    NA
## HC.23                          NA                    NA                    NA
## HC.24                          NA                     0                    NA
## HV.1                           NA                    NA                    NA
## HV.2                           NA                    NA                    NA
## HV.3                           NA                    NA                    NA
## HV.4                           NA                    NA                    NA
## HV.5                           NA                    NA                    NA
## HV.6                           NA                    NA                    NA
## HV.7                    1.380e+10                    NA             9.600e+12
## HV.8                    0.000e+00                    NA             0.000e+00
## HV.9                    2.080e+09                    NA             8.500e+11
## HV.10                   0.000e+00                    NA             0.000e+00
## HV.11                   1.570e+07                    NA             6.125e+09
## HV.12                   0.000e+00                    NA             0.000e+00
## HV.13                   8.450e+09                    NA             8.400e+12
## HV.14                   0.000e+00                    NA             0.000e+00
## HV.15                   2.690e+09                    NA             1.190e+12
## HV.16                   0.000e+00                    NA             0.000e+00
## HV.17                   8.830e+05                    NA             3.470e+08
## HV.18                   0.000e+00                    NA             0.000e+00
## HV.19                   5.165e+09                    NA             1.300e+12
## HV.20                   0.000e+00                    NA             0.000e+00
## HV.21                   1.670e+09                    NA             4.800e+11
## HV.22                   0.000e+00                    NA             0.000e+00
## HV.23                   1.290e+05                    NA             0.000e+00
## HV.24                   0.000e+00                    NA             0.000e+00
##       Model Antibiotic_mg.mL Fermentation Antibiotic Lactose_mM Glucose_mM
## D.1   Human               NA           NA       <NA>         NA         NA
## HC.1  Human               20            2        CTX      0.000      0.000
## HC.2  Human               20            2        CTX      0.000      0.000
## HC.3  Human              200            2        CTX      0.000      0.000
## HC.4  Human              200            2        CTX      0.000      0.000
## HC.5  Human               NA            2        CTX      0.000      0.000
## HC.6  Human               NA            2        CTX      0.000      0.000
## HC.7  Human               20            2        CTX      0.000      0.000
## HC.8  Human               20            2        CTX      0.000      0.000
## HC.9  Human              200            2        CTX      0.000      0.305
## HC.10 Human              200            2        CTX      0.000      0.309
## HC.11 Human               NA            2        CTX      0.000      0.000
## HC.12 Human               NA            2        CTX      0.000      0.000
## HC.13 Human               20            2        CTX      0.000      0.000
## HC.14 Human               20            2        CTX      0.000      0.000
## HC.15 Human              200            2        CTX      0.000      0.000
## HC.16 Human              200            2        CTX      0.000      0.000
## HC.17 Human               NA            2        CTX      0.000      0.000
## HC.18 Human               NA            2        CTX      0.000      0.000
## HC.19 Human               20            2        CTX      0.000      0.000
## HC.20 Human               20            2        CTX      0.000      0.000
## HC.21 Human              200            2        CTX      0.000      0.000
## HC.22 Human              200            2        CTX      0.000      0.000
## HC.23 Human               NA            2        CTX      0.000      0.000
## HC.24 Human               NA            2        CTX      0.000      0.000
## HV.1  Human               90            1        VAN      0.053      0.173
## HV.2  Human               90            1        VAN      0.000      0.228
## HV.3  Human              600            1        VAN      0.069      0.000
## HV.4  Human              600            1        VAN      0.000      0.000
## HV.5  Human               NA            1        VAN      0.037      0.284
## HV.6  Human               NA            1        VAN      0.046      0.000
## HV.7  Human               90            1        VAN      0.000      3.207
## HV.8  Human               90            1        VAN      0.000      9.781
## HV.9  Human              600            1        VAN      0.000      0.000
## HV.10 Human              600            1        VAN      0.000      6.682
## HV.11 Human               NA            1        VAN      0.000      0.000
## HV.12 Human               NA            1        VAN      0.000      0.000
## HV.13 Human               90            1        VAN      1.580      0.767
## HV.14 Human               90            1        VAN      1.238      1.616
## HV.15 Human              600            1        VAN      0.904      0.000
## HV.16 Human              600            1        VAN      1.663      0.000
## HV.17 Human               NA            1        VAN      0.000      0.000
## HV.18 Human               NA            1        VAN      0.000      0.000
## HV.19 Human               90            1        VAN      0.777      0.244
## HV.20 Human               90            1        VAN      0.731      0.000
## HV.21 Human              600            1        VAN      1.220      0.000
## HV.22 Human              600            1        VAN      1.579      0.000
## HV.23 Human               NA            1        VAN      0.000      0.000
## HV.24 Human               NA            1        VAN      0.000      0.000
##       Galactose_mM Succinat_mM Lactat_mM Formiat_mM Acetat_mM Propionat_mM
## D.1             NA          NA        NA         NA        NA           NA
## HC.1         0.536       1.519     1.338      5.898    68.571       27.257
## HC.2         0.588       2.608     1.287      8.273    71.535       26.884
## HC.3         0.655       1.531     1.310      7.826    76.633       27.355
## HC.4         0.622       1.513     1.238      4.607    70.801       25.104
## HC.5         0.602       1.448     1.321      7.221    73.189       36.735
## HC.6         0.703       1.356     1.423      6.459    52.167       24.281
## HC.7         0.314       1.558     1.613      6.574    75.484       32.194
## HC.8         0.703       1.497     1.750      9.178    71.548       29.822
## HC.9         0.989       1.200     3.095      6.354    69.382       34.611
## HC.10        0.958       1.300     3.283      0.000    64.353       34.868
## HC.11        0.680       1.435     1.293      6.330    70.102       38.008
## HC.12        0.464       1.390     1.453      5.867    52.980       25.113
## HC.13        0.367       1.370     1.552      3.816    70.605       26.932
## HC.14        0.000       1.085     1.505      5.782    68.421       23.025
## HC.15        0.391       1.100     2.160      5.931    65.980       27.327
## HC.16        0.397       1.051     2.611      0.000    64.822       36.192
## HC.17        0.000       1.110     1.029      5.727    67.573       39.064
## HC.18        0.000       1.258     0.905      5.953    54.688       26.113
## HC.19        0.229       1.194     1.038      3.661    69.252       25.772
## HC.20        0.759       1.121     2.450      0.000    59.247       25.199
## HC.21        0.000       1.739     1.127      5.850    66.606       24.514
## HC.22        0.000       1.796     1.322      0.000    58.881       33.442
## HC.23        0.000       1.730     1.331     10.916    67.988       35.881
## HC.24        0.000       1.784     0.853      4.163    61.137       27.614
## HV.1         0.000       0.901     0.534      5.278    77.645       23.407
## HV.2         0.000       0.886     0.734      6.244    78.499       27.767
## HV.3         0.000       0.844     0.679      7.344    57.376       36.797
## HV.4         0.000       0.932     0.449      7.921    69.627       22.762
## HV.5         0.000       1.135     0.364      3.277    76.749       19.707
## HV.6         0.000       0.819     0.000      8.203    64.866       17.193
## HV.7         4.096       1.351    29.016     10.029    23.392       19.437
## HV.8         4.362       1.558     5.440      2.891    21.269       19.543
## HV.9         2.502       1.962    31.934     10.365    19.888       17.952
## HV.10        4.270       1.587     5.702      2.731    29.523       14.073
## HV.11        0.000       1.239     1.215      5.294    78.115       17.219
## HV.12        0.000       1.246     0.994      9.137    71.280       19.527
## HV.13        2.044       1.473     9.844     13.503    11.513        7.533
## HV.14        1.323       1.637     5.033     12.418    20.830       25.674
## HV.15        2.563       1.559    15.308      8.918    11.689       13.911
## HV.16        2.668       1.479     7.130     14.424    12.683        6.561
## HV.17        0.000       1.233     1.269      6.174    73.792       17.652
## HV.18        0.000       1.251     1.113     10.938    74.408       23.448
## HV.19        1.325       1.417     5.496      8.836    18.665       41.793
## HV.20        1.144       1.378     4.138      8.562    19.052       45.481
## HV.21        0.784       1.435    13.537     10.322    10.928       11.060
## HV.22        0.750       1.401     4.294     12.778    15.865       10.169
## HV.23        0.266       1.208     1.372      6.666    72.238       14.474
## HV.24        0.000       1.187     1.171      9.060    64.405       19.489
##       Isobutyrat_mM Butyrat_mM Isovalerat_mM Valerat_mM Total_SCFA_mM
## D.1              NA         NA            NA         NA            NA
## HC.1          4.937     48.394         3.842      3.761       166.053
## HC.2          4.776     45.421         3.394      3.482       168.248
## HC.3          4.521     44.309         3.614      1.935       169.689
## HC.4          4.911     50.198         4.062      1.715       164.771
## HC.5          3.001     33.900         2.809      1.264       161.490
## HC.6          3.506     53.009         3.828      1.687       148.419
## HC.7          4.251     41.199         3.894      4.031       171.112
## HC.8          4.000     43.675         2.980      3.486       168.639
## HC.9          3.285     36.085         3.394      1.704       160.404
## HC.10         3.117     39.773         3.610      1.147       152.718
## HC.11         3.627     35.459         3.258      1.282       161.474
## HC.12         3.912     52.049         3.952      1.797       148.977
## HC.13         4.857     39.501         4.427      5.098       158.525
## HC.14         4.640     42.896         3.259      4.630       155.243
## HC.15         3.912     40.000         5.518      1.598       153.917
## HC.16         3.756     38.253         5.307      1.038       153.427
## HC.17         3.187     35.241         3.844      1.038       157.813
## HC.18         3.488     51.621         4.236      1.783       150.045
## HC.19         5.347     47.850         5.290      6.032       165.665
## HC.20         5.255     49.825         6.343      4.494       154.693
## HC.21         5.216     48.333         4.836      6.599       164.820
## HC.22         5.939     48.455         7.871      1.967       159.673
## HC.23         0.000     38.423         4.747      1.182       162.198
## HC.24         4.315     48.676         5.126      1.616       155.284
## HV.1          3.886     40.871         4.754      1.239       158.741
## HV.2          4.418     37.598         5.287      1.303       162.964
## HV.3          4.541     33.022         5.700      0.000       146.372
## HV.4          3.449     43.182         2.876      0.475       151.673
## HV.5          4.269     42.214         4.334      0.284       152.654
## HV.6          2.603     44.874         3.593      0.390       142.587
## HV.7          2.571     14.661         4.272      1.059       113.091
## HV.8          3.698     16.553         6.617      1.054        92.766
## HV.9          3.586     18.726         5.648      1.069       113.632
## HV.10         2.032     12.535         2.750      1.147        83.032
## HV.11         4.725     46.529         4.203      1.259       159.798
## HV.12         2.548     44.810         3.415      1.293       154.250
## HV.13         3.121      9.105         5.593      1.046        67.122
## HV.14         3.302     12.666         5.969      1.079        92.785
## HV.15         3.508      7.647         5.601      1.023        72.631
## HV.16         2.260      5.816         4.042      0.958        59.684
## HV.17         5.184     48.061         4.140      1.393       158.898
## HV.18         2.422     42.201         2.336      1.288       159.405
## HV.19         4.202     13.490         6.903      1.239       104.387
## HV.20         4.842     12.833         7.904      1.243       107.308
## HV.21         2.549      5.173         3.751      0.973        61.732
## HV.22         2.208      4.905         4.026      0.971        58.946
## HV.23         4.565     48.449         3.707      1.426       154.371
## HV.24         2.321     42.328         2.300      1.239       143.500
##       raw_metagenomic_pairs Period Reactor_Treatment_Dose   Treatment_Dose
## D.1                      NA   <NA>                  DONOR            DONOR
## HC.1                     NA   <NA>      TR1_CTX+HV292.120    CTX+HV292.120
## HC.2                     NA   <NA>              TR2_CTX20            CTX20
## HC.3                     NA   <NA>     TR3_CTX+HV292.1200   CTX+HV292.1200
## HC.4                     NA   <NA>             TR4_CTX200           CTX200
## HC.5                     NA   <NA>            TR5_HV292.1          HV292.1
## HC.6                     NA   <NA>           CR_UNTREATED        UNTREATED
## HC.7                     NA     t1      TR1_CTX+HV292.120    CTX+HV292.120
## HC.8                     NA     t1              TR2_CTX20            CTX20
## HC.9                     NA     t1     TR3_CTX+HV292.1200   CTX+HV292.1200
## HC.10                    NA     t1             TR4_CTX200           CTX200
## HC.11                    NA     t1            TR5_HV292.1          HV292.1
## HC.12                    NA     t1           CR_UNTREATED        UNTREATED
## HC.13                    NA     t1      TR1_CTX+HV292.120    CTX+HV292.120
## HC.14                    NA     t1              TR2_CTX20            CTX20
## HC.15                    NA     t1     TR3_CTX+HV292.1200   CTX+HV292.1200
## HC.16                    NA     t1             TR4_CTX200           CTX200
## HC.17                    NA     t1            TR5_HV292.1          HV292.1
## HC.18                    NA     t1           CR_UNTREATED        UNTREATED
## HC.19                    NA     t1      TR1_CTX+HV292.120    CTX+HV292.120
## HC.20                    NA     t1              TR2_CTX20            CTX20
## HC.21                    NA     t1     TR3_CTX+HV292.1200   CTX+HV292.1200
## HC.22                    NA     t1             TR4_CTX200           CTX200
## HC.23                    NA     t1            TR5_HV292.1          HV292.1
## HC.24                    NA     t1           CR_UNTREATED        UNTREATED
## HV.1                     NA   <NA>    TR1_VAN+CCUG5916890  VAN+CCUG5916890
## HV.2                     NA   <NA>              TR2_VAN90            VAN90
## HV.3                     NA   <NA>   TR3_VAN+CCUG59168600 VAN+CCUG59168600
## HV.4                     NA   <NA>             TR4_VAN600           VAN600
## HV.5                     NA   <NA>          TR5_CCUG59168        CCUG59168
## HV.6                     NA   <NA>           CR_UNTREATED        UNTREATED
## HV.7                     NA     t1    TR1_VAN+CCUG5916890  VAN+CCUG5916890
## HV.8                     NA     t1              TR2_VAN90            VAN90
## HV.9                     NA     t1   TR3_VAN+CCUG59168600 VAN+CCUG59168600
## HV.10                    NA     t1             TR4_VAN600           VAN600
## HV.11                    NA     t1          TR5_CCUG59168        CCUG59168
## HV.12                    NA     t1           CR_UNTREATED        UNTREATED
## HV.13                    NA     t1    TR1_VAN+CCUG5916890  VAN+CCUG5916890
## HV.14                    NA     t1              TR2_VAN90            VAN90
## HV.15                    NA     t1   TR3_VAN+CCUG59168600 VAN+CCUG59168600
## HV.16                    NA     t1             TR4_VAN600           VAN600
## HV.17                    NA     t1          TR5_CCUG59168        CCUG59168
## HV.18                    NA     t1           CR_UNTREATED        UNTREATED
## HV.19                    NA     t1    TR1_VAN+CCUG5916890  VAN+CCUG5916890
## HV.20                    NA     t1              TR2_VAN90            VAN90
## HV.21                    NA     t1   TR3_VAN+CCUG59168600 VAN+CCUG59168600
## HV.22                    NA     t1             TR4_VAN600           VAN600
## HV.23                    NA     t1          TR5_CCUG59168        CCUG59168
## HV.24                    NA     t1           CR_UNTREATED        UNTREATED
##       Day_of_Treatment_num Antibiotic_mg.L
## D.1                     NA            <NA>
## HC.1                    -1              20
## HC.2                    -1              20
## HC.3                    -1             200
## HC.4                    -1             200
## HC.5                    -1            <NA>
## HC.6                    -1            <NA>
## HC.7                     1              20
## HC.8                     1              20
## HC.9                     1             200
## HC.10                    1             200
## HC.11                    1            <NA>
## HC.12                    1            <NA>
## HC.13                    3              20
## HC.14                    3              20
## HC.15                    3             200
## HC.16                    3             200
## HC.17                    3            <NA>
## HC.18                    3            <NA>
## HC.19                    7              20
## HC.20                    7              20
## HC.21                    7             200
## HC.22                    7             200
## HC.23                    7            <NA>
## HC.24                    7            <NA>
## HV.1                    -1              90
## HV.2                    -1              90
## HV.3                    -1             600
## HV.4                    -1             600
## HV.5                    -1            <NA>
## HV.6                    -1            <NA>
## HV.7                     1              90
## HV.8                     1              90
## HV.9                     1             600
## HV.10                    1             600
## HV.11                    1            <NA>
## HV.12                    1            <NA>
## HV.13                    3              90
## HV.14                    3              90
## HV.15                    3             600
## HV.16                    3             600
## HV.17                    3            <NA>
## HV.18                    3            <NA>
## HV.19                    7              90
## HV.20                    7              90
## HV.21                    7             600
## HV.22                    7             600
## HV.23                    7            <NA>
## HV.24                    7            <NA>
```

```r
 # write_tsv("~/Desktop/samples_data_h1.tsv")
```

combine:


```r
# ps_AMR %>% 
#   physeq_add_metadata(metadata_mapping,
#                       sample_column = "sample_id") %>% 
#   sample_data() %>% 
#   data.frame() %>% 
#   filter(Model == "Human") %>% 
#   left_join(.,
#             sample_AMR_table,
#             by = c("sample_id" = "Sample")) %>% 
#   left_join(sample_carrier_table,
#             by = c("sample_id" = "Sample")) %>% 
#   select(-Description2, -Experiment, -Day_of_Connection, -Day_from_Inoculum, -Phase, -Date:-VAN_Copy_Number_permL,
#          -raw_metagenomic_pairs, -Period, -Day_of_Treatment_num) %>% 
#   # mutate(Others = sum(`pleuromutilin antibiotic`, `para-aminosalicylic acid`, macrolide, `fusidic acid`, carbapenem, `acridine dye` ,`mupirocin`, na.rm = TRUE)) %>% 
#   select(-`pleuromutilin antibiotic`:-`mupirocin`) %>% 
#   write_tsv("~/Desktop/samples_data_h1_end.tsv")
```



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
##  date     2022-09-27
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
##  cowplot            1.1.1      2020-12-30 [1] CRAN (R 4.2.0)
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
##  ggh4x              0.2.1.9000 2022-08-08 [1] Github (teunbrand/ggh4x@71446d4)
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
