---
title: " Fig. 1 - S1 - chicken manuscript "
author: "Florentin Constancias"
date: "December 23, 2022"
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
## ✔ tidyr   1.2.1      ✔ stringr 1.5.0 
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
require(microViz); packageVersion("microViz")
```

```
## Loading required package: microViz
## 
## microViz version 0.9.7 - Copyright (C) 2022 David Barnett
## * Website: https://david-barnett.github.io/microViz/
## * Useful? For citation info, run: citation('microViz')
## * Silence: suppressPackageStartupMessages(library(microViz))
```


```r
source("https://raw.githubusercontent.com/fconstancias/DivComAnalyses/master/R/phyloseq_taxa_tests.R") 
source("https://raw.githubusercontent.com/fconstancias/DivComAnalyses/master/R/phyloseq_alpha.R") 
source("https://raw.githubusercontent.com/fconstancias/DivComAnalyses/master/R/phyloseq_beta.R") 
source("https://raw.githubusercontent.com/fconstancias/DivComAnalyses/master/R/phyloseq_heatmap.R")
# source("https://raw.githubusercontent.com/fconstancias/metabaRpipe-source/master/Rscripts/functions.R")
source("https://raw.githubusercontent.com/fconstancias/DivComAnalyses/master/R/phyloseq_normalisation.R") 
source("https://raw.githubusercontent.com/fconstancias/DivComAnalyses/master/R/phyloseq_varia.R")


'%!in%' <- function(x,y)!('%in%'(x,y))
```


```r
out_pptx = "~/Desktop/16S-human_S1.pptx"
out_xlsx = "~/Desktop/16S-human_S1.xlsx"
```


```r
load(url("https://github.com/fconstancias/NRP72-FBT/blob/master/Figures/Rastetics.Rdata?raw=true"))
```

# load data:


```r
"data/processed/16S/16S_working_phyloseq.RDS" %>% 
  here::here() %>% 
  readRDS()  %>% 
  # subset_samples(Enrichment == "NotEnriched") %>%
  # subset_samples(Experiment %in% c("Continuous", "Cecum")) %>%
  subset_samples(Day_of_Treatment > -6 & Day_of_Treatment <= 30 | Reactor == "DONOR") %>% 
  subset_samples(Model == "Human") %>% 
  # subset_samples(Reactor != "CR" & Model2 == "Chicken1" | Model2 == "Chicken2" & Reactor %in% c("TR2", "TR3", "TR4", "TR5", "TR6", "TR7", "CR", "DONOR")) %>% 
  subset_samples(Reactor != "IR") -> ps_filtered
```

# Rarefraction


```r
min_sample_size = 3738

ps_filtered %>%
  filter_taxa(function(x) sum(x > 0) > 0, TRUE) %>% 
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
## 35 samples removedbecause they contained fewer reads than `sample.size`.
```

```
## Up to first five removed samples are:
```

```
## D-1-S49TR1-31-S109TR1-32-S111TR1-60-S114TR1-63-S100
```

```
## ...
```

```
## 135OTUs were removed because they are no longer 
## present in any sample after random subsampling
```

```
## ...
```

```r
ps_filtered %>%
  filter_taxa(function(x) sum(x > 0) > 0, TRUE) %>% 
  prune_samples(sample_sums(.)>= min_sample_size, .) -> ps_fil
```

## Before Treatment:

Obj: compare donor microbiota and pre-treatment / CR - heatmap
day -2, -1, 0
top feces

Identify last 3 samples including 0:


```r
ps_rare %>%
  filter_taxa(function(x) sum(x > 0) > 0, TRUE) %>% 
  subset_samples(Day_of_Treatment <= 0 & Day_of_Treatment > -4 | Reactor %in% c("DONOR")) -> before_treat # %>% 

before_treat %>% 
  sample_data() %>% 
  data.frame() %>% 
  group_by(Model, Model2, Reactor_Treatment_Dose_fermentation) %>% 
  arrange(Day_of_Treatment) %>% 
  slice(tail(row_number(), 4)) %>% 
  pull(sample_names) %>% 
  as.character() -> before_treat_samples
```


```r
ps_rare %>%
  filter_taxa(function(x) sum(x > 0) > 0, TRUE) %>% 
  subset_samples(sample_name %in% c(before_treat_samples)) -> before_treat
# prune_samples(x = ., samples  %in% c(before_treat_samples)) -> before_treat #%>% 
# rarefy_even_depth(sample.size = min_sample_size, rngseed = 123) %>%
# subset_samples(Day_of_Treatment %in% c(-2, -1, 0) | Reactor %in% c("DONOR")) -> before_treat # %>% 
# subset_samples(Reactor %!in%c("IR")) %>% 
# subset_samples() -> before_treat
```

### ASV:

top ASV  per fecal samples in each donors:


```r
ntax = 20 

before_treat %>%
  subset_samples(Reactor == "DONOR") %>% 
  physeq_most_abundant(group_var = "Model2",
                       ntax = ntax,
                       tax_level = "Strain") -> top_ASV_donor_feces

# top_ASV_donor_feces

before_treat %>%
  physeq_most_abundant(group_var = "Reactor",
                       ntax = ntax,
                       tax_level = "Strain") -> top_ASV_group

# top_ASV_group
```

#### top feces: 


```r
before_treat %>%
  filter_taxa(function(x) sum(x > 0) > 0, TRUE) %>% 
  transform_sample_counts(function(x) x/sum(x) * 100) %>%
  subset_taxa(Strain %in% top_ASV_donor_feces) %>% # extract only the taxa to display - after percentage normalisation
  phyloseq_ampvis_heatmap(physeq = ., tax_aggregate = "Species",
                          transform = FALSE,
                          plot_values = FALSE,
                          facet_by = c("Model2", "Model", "Fermentation","Reactor_Treatment_Dose"),
                          group_by = "Day_of_Treatment",#"sample_name",
                          ntax = Inf
  ) -> pre_speces_feces
```

```
## Loading required package: ampvis2
```

```
## Warning: There are only 33 taxa, showing all
```

```r
pre_speces_feces$data %>% 
  mutate(Abundance = na_if(Abundance, 0)) %>% 
  mutate(Fermentation = Fermentation %>%  as.character() %>%  replace_na(., "NA")) %>% 
  mutate(Fermentation = factor(Fermentation, levels = c("NA", 1, 2))) -> pre_speces_feces$data 

pre_speces_feces + 
  facet_grid(Model ~  Model2 + Fermentation + Reactor_Treatment_Dose, scales = "free", space = "free") +
  scale_color_manual(values = c("#31a354", "orange", 
                                "#f03b20"), drop = FALSE, na.value = 'transparent') -> pre_speces_feces
# scale_fill_viridis_c(breaks = c(0,  0.01, 1, 10, 25, 50, 75, 100),
# labels = c(0,  0.01, 1, 10, 25,  50, 75, 100),
# trans = scales::pseudo_log_trans(sigma = 0.9),
# na.value = 'transparent') -> pre_speces

pre_speces_feces
```

<img src="manuscript_Figure1-S1_files/figure-html/unnamed-chunk-8-1.png" style="display: block; margin: auto;" />

```r
# pre_speces %>% 
#   ggsave(file="~/Desktop/16S-chick_S1.svg", plot=., width=10, height=8)
```

Representative of those:


```r
# before_treat %>%
#   filter_taxa(function(x) sum(x > 0) > 0, TRUE) %>% 
#   transform_sample_counts(function(x) x/sum(x) * 100) %>%
#   subset_taxa(Strain %in% top_ASV_donor_feces) %>%
#   ps_arrange(Day_of_Treatment) %>% 
#   speedyseq::psmelt() %>%
#   # group_by(Sample, Day_of_Treatment, Fermentation, Model, Model2, Reactor_Treatment_Dose) %>%
#   group_by(sample_names, Model, Model2, Day_of_Treatment, Fermentation , Reactor_Treatment_Dose) %>% 
#   summarise(coverage = sum(Abundance)) %>%
#   ungroup() %>%
#   mutate(Day_of_Treatment = as.double(Day_of_Treatment)) %>%
#   mutate_at(vars(Day_of_Treatment), ~replace_na(., 0)) %>% 
#       # arrange(Day_of_Treatment) %>%
#   # group_by(sample_names, Model, Model2, Day_of_Treatment, Reactor_Treatment_Dose_fermentation) %>% 
#   #      summarise(count = length(coverage)) %>% 
#   # ggpubr::ggdotchart(., x = "Day_of_Treatment", y = "coverage",
#   #                    fill = "grey20", color = "grey5", width = 0.8, add = "segments") + #,
#   # ylab = "gdh copies / 10^10 cells \n", xlab = "") +
#   # facet_grid(variable ~ ., scales = "free") +
#   ggplot(., aes(x=Day_of_Treatment, y=coverage)) +
#   geom_point() + ggpubr::rotate_x_text(60)  +
#   theme_light() + 
#   # facet_grid(Model ~  Model2 + Reactor_Treatment_Dose_fermentation,  scales = "fixed", space = "fixed", drop = TRUE) -> toto
#   facet_grid(Model ~  Model2 + Fermentation + Reactor_Treatment_Dose, scales = "free", space = "free", drop = FALSE) -> toto
# 
# toto



before_treat %>%
  filter_taxa(function(x) sum(x > 0) > 0, TRUE) %>% 
  transform_sample_counts(function(x) x/sum(x) * 100) %>%
  subset_taxa(Strain %in% top_ASV_donor_feces) %>%
  ps_arrange(Day_of_Treatment) %>% 
  speedyseq::psmelt() %>%
  # group_by(Sample, Day_of_Treatment, Fermentation, Model, Model2, Reactor_Treatment_Dose) %>%
  group_by(sample_names, Model, Model2, Day_of_Treatment, Fermentation , Reactor_Treatment_Dose) %>% 
  summarise(coverage = sum(Abundance)) %>%
  ungroup() %>%
  mutate(Day_of_Treatment = as.double(Day_of_Treatment)) %>%
  mutate_at(vars(Day_of_Treatment), ~replace_na(., 0)) %>%
  arrange(Day_of_Treatment) %>%
  # group_by(sample_names, Model, Model2, Day_of_Treatment, Reactor_Treatment_Dose_fermentation) %>% 
  #      summarise(count = length(coverage)) %>% 
  ggpubr::ggbarplot(., x = "Day_of_Treatment", y = "coverage",
                    fill = "grey20", color = "grey5", width = 0.8, add = "segments") + #,
  # ylab = "gdh copies / 10^10 cells \n", xlab = "") +
  # facet_grid(variable ~ ., scales = "free") +
  # ggplot(., aes(x=Day_of_Treatment, y=coverage)) +
  # geom_point() +
  ggpubr::rotate_x_text(60)  +
  theme_light() + 
  facet_grid(Model ~  Model2 + Fermentation + Reactor_Treatment_Dose, scales = "free", space = "free") -> repr_top_ASV_donor_feces
```

```
## `summarise()` has grouped output by 'sample_names', 'Model', 'Model2',
## 'Day_of_Treatment', 'Fermentation'. You can override using the `.groups`
## argument.
```

```r
# facet_grid(Model ~  Model2 + Reactor_Treatment_Dose_fermentation,  scales = "fixed", space = "fixed", drop = TRUE) -> toto
# facet_grid(Model ~  Model2 + Fermentation + Reactor_Treatment_Dose, scales = "free", space = "free", drop = FALSE) -> toto
# facet_wrap()

repr_top_ASV_donor_feces
```

<img src="manuscript_Figure1-S1_files/figure-html/unnamed-chunk-9-1.png" style="display: block; margin: auto;" />

#### top top_ASV_group: 


```r
before_treat %>%
  filter_taxa(function(x) sum(x > 0) > 0, TRUE) %>% 
  transform_sample_counts(function(x) x/sum(x) * 100) %>%
  subset_taxa(Strain %in% top_ASV_group) %>% # extract only the taxa to display - after percentage normalisation
  phyloseq_ampvis_heatmap(physeq = ., tax_aggregate = "Species",
                          transform = FALSE,
                          plot_values = FALSE,
                          facet_by = c("Model2", "Model", "Fermentation","Reactor_Treatment_Dose"),
                          group_by = "Day_of_Treatment",#"sample_name",
                          ntax = Inf
  ) -> pre_speces_group
```

```
## Warning: There are only 53 taxa, showing all
```

```r
pre_speces_group$data %>% 
  mutate(Abundance = na_if(Abundance, 0)) %>% 
  mutate(Fermentation = Fermentation %>%  as.character() %>%  replace_na(., "NA")) %>% 
  mutate(Fermentation = factor(Fermentation, levels = c("NA", 1, 2))) -> pre_speces_group$data 

pre_speces_group + 
  facet_grid(Model ~   Model2 + Fermentation + Reactor_Treatment_Dose, scales = "free", space = "free") +
  scale_color_manual(values = c("#31a354", "orange", 
                                "#f03b20"), drop = FALSE, na.value = 'transparent') -> pre_speces_group
# scale_fill_viridis_c(breaks = c(0,  0.01, 1, 10, 25, 50, 75, 100),
# labels = c(0,  0.01, 1, 10, 25,  50, 75, 100),
# trans = scales::pseudo_log_trans(sigma = 0.9),
# na.value = 'transparent') -> pre_speces

pre_speces_group
```

<img src="manuscript_Figure1-S1_files/figure-html/unnamed-chunk-10-1.png" style="display: block; margin: auto;" />

```r
# pre_speces %>% 
#   ggsave(file="~/Desktop/16S-chick_S1.svg", plot=., width=10, height=8)
```

Representative of those:


```r
# before_treat %>%
#   filter_taxa(function(x) sum(x > 0) > 0, TRUE) %>% 
#   transform_sample_counts(function(x) x/sum(x) * 100) %>%
#   subset_taxa(Strain %in% top_ASV_donor_feces) %>%
#   ps_arrange(Day_of_Treatment) %>% 
#   speedyseq::psmelt() %>%
#   # group_by(Sample, Day_of_Treatment, Fermentation, Model, Model2, Reactor_Treatment_Dose) %>%
#   group_by(sample_names, Model, Model2, Day_of_Treatment, Fermentation , Reactor_Treatment_Dose) %>% 
#   summarise(coverage = sum(Abundance)) %>%
#   ungroup() %>%
#   mutate(Day_of_Treatment = as.double(Day_of_Treatment)) %>%
#   mutate_at(vars(Day_of_Treatment), ~replace_na(., 0)) %>% 
#       # arrange(Day_of_Treatment) %>%
#   # group_by(sample_names, Model, Model2, Day_of_Treatment, Reactor_Treatment_Dose_fermentation) %>% 
#   #      summarise(count = length(coverage)) %>% 
#   # ggpubr::ggdotchart(., x = "Day_of_Treatment", y = "coverage",
#   #                    fill = "grey20", color = "grey5", width = 0.8, add = "segments") + #,
#   # ylab = "gdh copies / 10^10 cells \n", xlab = "") +
#   # facet_grid(variable ~ ., scales = "free") +
#   ggplot(., aes(x=Day_of_Treatment, y=coverage)) +
#   geom_point() + ggpubr::rotate_x_text(60)  +
#   theme_light() + 
#   # facet_grid(Model ~  Model2 + Reactor_Treatment_Dose_fermentation,  scales = "fixed", space = "fixed", drop = TRUE) -> toto
#   facet_grid(Model ~  Model2 + Fermentation + Reactor_Treatment_Dose, scales = "free", space = "free", drop = FALSE) -> toto
# 
# toto



before_treat %>%
  filter_taxa(function(x) sum(x > 0) > 0, TRUE) %>% 
  transform_sample_counts(function(x) x/sum(x) * 100) %>%
  subset_taxa(Strain %in% top_ASV_group) %>%
  ps_arrange(Day_of_Treatment) %>% 
  speedyseq::psmelt() %>%
  # group_by(Sample, Day_of_Treatment, Fermentation, Model, Model2, Reactor_Treatment_Dose) %>%
  group_by(sample_names, Model, Model2, Day_of_Treatment, Fermentation , Reactor_Treatment_Dose) %>% 
  summarise(coverage = sum(Abundance)) %>%
  ungroup() %>%
  mutate(Day_of_Treatment = as.double(Day_of_Treatment)) %>%
  # mutate_at(vars(Day_of_Treatment), ~replace_na(., 0)) %>%
  arrange(Day_of_Treatment) %>%
  # group_by(sample_names, Model, Model2, Day_of_Treatment, Reactor_Treatment_Dose_fermentation) %>% 
  #      summarise(count = length(coverage)) %>% 
  ggpubr::ggbarplot(., x = "Day_of_Treatment", y = "coverage",
                    fill = "grey20", color = "grey5", width = 0.8, add = "segments") + #,
  # ylab = "gdh copies / 10^10 cells \n", xlab = "") +
  # facet_grid(variable ~ ., scales = "free") +
  # ggplot(., aes(x=Day_of_Treatment, y=coverage)) +
  # geom_point() +
  ggpubr::rotate_x_text(60)  +
  theme_light() + 
  facet_grid(Model ~  Model2 + Fermentation + Reactor_Treatment_Dose, scales = "free", space = "free") -> repr_top_ASV_pre_speces_group
```

```
## `summarise()` has grouped output by 'sample_names', 'Model', 'Model2',
## 'Day_of_Treatment', 'Fermentation'. You can override using the `.groups`
## argument.
```

```r
# facet_grid(Model ~  Model2 + Reactor_Treatment_Dose_fermentation,  scales = "fixed", space = "fixed", drop = TRUE) -> toto
# facet_grid(Model ~  Model2 + Fermentation + Reactor_Treatment_Dose, scales = "free", space = "free", drop = FALSE) -> toto
# facet_wrap()

repr_top_ASV_pre_speces_group
```

<img src="manuscript_Figure1-S1_files/figure-html/unnamed-chunk-11-1.png" style="display: block; margin: auto;" />
### Family:

top   per fecal samples in each donors:


```r
ntax = 6 

before_treat %>%
  subset_samples(Reactor == "DONOR") %>% 
  physeq_most_abundant(group_var = "Model2",
                       ntax = ntax,
                       tax_level = "Family") -> top_fam_donor_feces

# top_ASV_donor_feces

before_treat %>%
  physeq_most_abundant(group_var = "Reactor",
                       ntax = ntax,
                       tax_level = "Family") -> top_fam_group

# top_ASV_group
```



```r
before_treat %>%
  filter_taxa(function(x) sum(x > 0) > 0, TRUE) %>% 
  # subset_samples(Model2 == "Chicken1") %>% 
  # subset_samples(Reactor == "DONOR") %>%  sample_data() %>%  data.frame()
  speedyseq::tax_glom(taxrank = "Family") %>%
  transform_sample_counts(function(x) x/sum(x) * 100) %>%
  subset_taxa(Family %in% top_fam_donor_feces) %>% # extract only the taxa to display - after percentage normalisation
  
  phyloseq_ampvis_heatmap(physeq = ., tax_aggregate = "Family",
                          transform = FALSE,
                          plot_values = FALSE,
                          facet_by = c("Model2", "Model", "Fermentation","Reactor_Treatment_Dose"),
                          group_by = "Day_of_Treatment",#"sample_name",
                          ntax = Inf
  ) -> pre_Family
```

```
## Warning: There are only 4 taxa, showing all
```

```r
pre_Family$data %>% 
  mutate(Abundance = na_if(Abundance, 0)) %>% 
  mutate(Fermentation = Fermentation %>%  as.character() %>%  replace_na(., "NA")) %>% 
  mutate(Fermentation = factor(Fermentation, levels = c("NA", 1, 2))) -> pre_Family$data 


pre_Family + 
  facet_grid(Model ~ Model2 + Fermentation + Reactor_Treatment_Dose, scales = "free", space = "free") +
  scale_color_manual(values = c("#31a354", "orange", 
                                "#f03b20"), drop = FALSE, na.value = 'transparent') -> pre_Family


pre_Family
```

<img src="manuscript_Figure1-S1_files/figure-html/unnamed-chunk-13-1.png" style="display: block; margin: auto;" />

```r
before_treat %>%
  filter_taxa(function(x) sum(x > 0) > 0, TRUE) %>% 
  transform_sample_counts(function(x) x/sum(x) * 100) %>%
  subset_taxa(Family %in% top_fam_donor_feces) %>% 
  ps_arrange(Day_of_Treatment) %>% 
  speedyseq::psmelt() %>%
  # group_by(Sample, Day_of_Treatment, Fermentation, Model, Model2, Reactor_Treatment_Dose) %>%
  group_by(sample_names, Model, Model2, Day_of_Treatment, Fermentation , Reactor_Treatment_Dose) %>% 
  summarise(coverage = sum(Abundance)) %>%
  ungroup() %>%
  mutate(Day_of_Treatment = as.double(Day_of_Treatment)) %>%
  # mutate_at(vars(Day_of_Treatment), ~replace_na(., 0)) %>%
  arrange(Day_of_Treatment) %>%
  # group_by(sample_names, Model, Model2, Day_of_Treatment, Reactor_Treatment_Dose_fermentation) %>% 
  #      summarise(count = length(coverage)) %>% 
  ggpubr::ggbarplot(., x = "Day_of_Treatment", y = "coverage",
                    fill = "grey20", color = "grey5", width = 0.8, add = "segments") + #,
  # ylab = "gdh copies / 10^10 cells \n", xlab = "") +
  # facet_grid(variable ~ ., scales = "free") +
  # ggplot(., aes(x=Day_of_Treatment, y=coverage)) +
  # geom_point() +
  ggpubr::rotate_x_text(60)  +
  theme_light() + 
  facet_grid(Model ~   Model2 + Fermentation + Reactor_Treatment_Dose, scales = "free", space = "free") -> repr_top_Fam_donor_feces
```

```
## `summarise()` has grouped output by 'sample_names', 'Model', 'Model2',
## 'Day_of_Treatment', 'Fermentation'. You can override using the `.groups`
## argument.
```

```r
# facet_grid(Model ~  Model2 + Reactor_Treatment_Dose_fermentation,  scales = "fixed", space = "fixed", drop = TRUE) -> toto
# facet_grid(Model ~  Model2 + Fermentation + Reactor_Treatment_Dose, scales = "free", space = "free", drop = FALSE) -> toto
# facet_wrap()

repr_top_Fam_donor_feces
```

<img src="manuscript_Figure1-S1_files/figure-html/unnamed-chunk-14-1.png" style="display: block; margin: auto;" />



```r
before_treat %>%
  filter_taxa(function(x) sum(x > 0) > 0, TRUE) %>% 
  # subset_samples(Model2 == "Chicken1") %>% 
  # subset_samples(Reactor == "DONOR") %>%  sample_data() %>%  data.frame()
  speedyseq::tax_glom(taxrank = "Family") %>%
  transform_sample_counts(function(x) x/sum(x) * 100) %>%
  subset_taxa(Family %in% top_fam_group) %>% # extract only the taxa to display - after percentage normalisation
  
  phyloseq_ampvis_heatmap(physeq = ., tax_aggregate = "Family",
                          transform = FALSE,
                          plot_values = FALSE,
                          facet_by = c("Model2", "Model", "Fermentation","Reactor_Treatment_Dose"),
                          group_by = "Day_of_Treatment",#"sample_name",
                          ntax = Inf
  ) -> pre_Family_grp
```

```
## Warning: There are only 5 taxa, showing all
```

```r
pre_Family_grp$data %>% 
  mutate(Abundance = na_if(Abundance, 0)) %>% 
  mutate(Fermentation = Fermentation %>%  as.character() %>%  replace_na(., "NA")) %>% 
  mutate(Fermentation = factor(Fermentation, levels = c("NA", 1, 2))) -> pre_Family_grp$data 


pre_Family_grp + 
  facet_grid(Model ~ Model2 + Fermentation + Reactor_Treatment_Dose, scales = "free", space = "free") +
  scale_color_manual(values = c("#31a354", "orange", 
                                "#f03b20"), drop = FALSE, na.value = 'transparent') -> pre_Family_grp


pre_Family_grp
```

<img src="manuscript_Figure1-S1_files/figure-html/unnamed-chunk-15-1.png" style="display: block; margin: auto;" />

```r
before_treat %>%
  filter_taxa(function(x) sum(x > 0) > 0, TRUE) %>% 
  transform_sample_counts(function(x) x/sum(x) * 100) %>%
  subset_taxa(Family %in% top_fam_group) %>% 
  ps_arrange(Day_of_Treatment) %>% 
  speedyseq::psmelt() %>%
  # group_by(Sample, Day_of_Treatment, Fermentation, Model, Model2, Reactor_Treatment_Dose) %>%
  group_by(sample_names, Model, Model2, Day_of_Treatment, Fermentation , Reactor_Treatment_Dose) %>% 
  summarise(coverage = sum(Abundance)) %>%
  ungroup() %>%
  mutate(Day_of_Treatment = as.double(Day_of_Treatment)) %>%
  # mutate_at(vars(Day_of_Treatment), ~replace_na(., 0)) %>%
  arrange(Day_of_Treatment) %>%
  # group_by(sample_names, Model, Model2, Day_of_Treatment, Reactor_Treatment_Dose_fermentation) %>% 
  #      summarise(count = length(coverage)) %>% 
  ggpubr::ggbarplot(., x = "Day_of_Treatment", y = "coverage",
                    fill = "grey20", color = "grey5", width = 0.8, add = "segments") + #,
  # ylab = "gdh copies / 10^10 cells \n", xlab = "") +
  # facet_grid(variable ~ ., scales = "free") +
  # ggplot(., aes(x=Day_of_Treatment, y=coverage)) +
  # geom_point() +
  ggpubr::rotate_x_text(60)  +
  theme_light() + 
  facet_grid(Model ~   Model2+ Fermentation + Reactor_Treatment_Dose, scales = "free", space = "free") -> repr_top_Fam_donor_grp
```

```
## `summarise()` has grouped output by 'sample_names', 'Model', 'Model2',
## 'Day_of_Treatment', 'Fermentation'. You can override using the `.groups`
## argument.
```

```r
# facet_grid(Model ~  Model2 + Reactor_Treatment_Dose_fermentation,  scales = "fixed", space = "fixed", drop = TRUE) -> toto
# facet_grid(Model ~  Model2 + Fermentation + Reactor_Treatment_Dose, scales = "free", space = "free", drop = FALSE) -> toto
# facet_wrap()

repr_top_Fam_donor_grp
```

<img src="manuscript_Figure1-S1_files/figure-html/unnamed-chunk-16-1.png" style="display: block; margin: auto;" />

### Genus:

top  per fecal samples in each donors:


```r
ntax = 10

before_treat %>%
  subset_samples(Reactor == "DONOR") %>% 
  physeq_most_abundant(group_var = "Model2",
                       ntax = ntax,
                       tax_level = "Genus") -> top_gen_donor_feces

# top_ASV_donor_feces

before_treat %>%
  physeq_most_abundant(group_var = "Reactor",
                       ntax = ntax,
                       tax_level = "Genus") -> top_gen_group

# top_ASV_group
```



```r
before_treat %>%
  filter_taxa(function(x) sum(x > 0) > 0, TRUE) %>% 
  # subset_samples(Model2 == "Chicken1") %>% 
  # subset_samples(Reactor == "DONOR") %>%  sample_data() %>%  data.frame()
  speedyseq::tax_glom(taxrank = "Genus") %>%
  transform_sample_counts(function(x) x/sum(x) * 100) %>%
  subset_taxa(Genus %in% top_gen_donor_feces) %>% # extract only the taxa to display - after percentage normalisation
  
  phyloseq_ampvis_heatmap(physeq = ., tax_aggregate = "Genus",
                          transform = FALSE,
                          plot_values = FALSE,
                          facet_by = c("Model2", "Model", "Fermentation","Reactor_Treatment_Dose"),
                          group_by = "Day_of_Treatment",#"sample_name",
                          ntax = Inf
  ) -> pre_gen
```

```
## Warning: There are only 7 taxa, showing all
```

```r
pre_gen$data %>% 
  mutate(Abundance = na_if(Abundance, 0)) %>% 
  mutate(Fermentation = Fermentation %>%  as.character() %>%  replace_na(., "NA")) %>% 
  mutate(Fermentation = factor(Fermentation, levels = c("NA", 1, 2))) -> pre_gen$data 


pre_gen + 
  facet_grid(Model ~  Model2 + Fermentation + Reactor_Treatment_Dose, scales = "free", space = "free") +
  scale_color_manual(values = c("#31a354", "orange", 
                                "#f03b20"), drop = FALSE, na.value = 'transparent') -> pre_gen


pre_gen
```

<img src="manuscript_Figure1-S1_files/figure-html/unnamed-chunk-18-1.png" style="display: block; margin: auto;" />

```r
before_treat %>%
  filter_taxa(function(x) sum(x > 0) > 0, TRUE) %>% 
  transform_sample_counts(function(x) x/sum(x) * 100) %>%
  subset_taxa(Genus %in% top_gen_donor_feces) %>% 
  ps_arrange(Day_of_Treatment) %>% 
  speedyseq::psmelt() %>%
  # group_by(Sample, Day_of_Treatment, Fermentation, Model, Model2, Reactor_Treatment_Dose) %>%
  group_by(sample_names, Model, Model2, Day_of_Treatment, Fermentation , Reactor_Treatment_Dose) %>% 
  summarise(coverage = sum(Abundance)) %>%
  ungroup() %>%
  mutate(Day_of_Treatment = as.double(Day_of_Treatment)) %>%
  # mutate_at(vars(Day_of_Treatment), ~replace_na(., 0)) %>%
  arrange(Day_of_Treatment) %>%
  # group_by(sample_names, Model, Model2, Day_of_Treatment, Reactor_Treatment_Dose_fermentation) %>% 
  #      summarise(count = length(coverage)) %>% 
  ggpubr::ggbarplot(., x = "Day_of_Treatment", y = "coverage",
                    fill = "grey20", color = "grey5", width = 0.8, add = "segments") + #,
  # ylab = "gdh copies / 10^10 cells \n", xlab = "") +
  # facet_grid(variable ~ ., scales = "free") +
  # ggplot(., aes(x=Day_of_Treatment, y=coverage)) +
  # geom_point() +
  ggpubr::rotate_x_text(60)  +
  theme_light() + 
  facet_grid(Model ~  Model2 + Fermentation + Reactor_Treatment_Dose, scales = "free", space = "free") -> repr_top_gen_donor_feces
```

```
## `summarise()` has grouped output by 'sample_names', 'Model', 'Model2',
## 'Day_of_Treatment', 'Fermentation'. You can override using the `.groups`
## argument.
```

```r
# facet_grid(Model ~  Model2 + Reactor_Treatment_Dose_fermentation,  scales = "fixed", space = "fixed", drop = TRUE) -> toto
# facet_grid(Model ~  Model2 + Fermentation + Reactor_Treatment_Dose, scales = "free", space = "free", drop = FALSE) -> toto
# facet_wrap()

repr_top_gen_donor_feces
```

<img src="manuscript_Figure1-S1_files/figure-html/unnamed-chunk-19-1.png" style="display: block; margin: auto;" />



```r
before_treat %>%
  filter_taxa(function(x) sum(x > 0) > 0, TRUE) %>% 
  # subset_samples(Model2 == "Chicken1") %>% 
  # subset_samples(Reactor == "DONOR") %>%  sample_data() %>%  data.frame()
  speedyseq::tax_glom(taxrank = "Genus") %>%
  transform_sample_counts(function(x) x/sum(x) * 100) %>%
  subset_taxa(Genus %in% top_gen_group) %>% # extract only the taxa to display - after percentage normalisation
  
  phyloseq_ampvis_heatmap(physeq = ., tax_aggregate = "Genus",
                          transform = FALSE,
                          plot_values = FALSE,
                          facet_by = c("Model2", "Model", "Fermentation","Reactor_Treatment_Dose"),
                          group_by = "Day_of_Treatment",#"sample_name",
                          ntax = Inf
  ) -> pre_gen_grp
```

```
## Warning: There are only 14 taxa, showing all
```

```r
pre_gen_grp$data %>% 
  mutate(Abundance = na_if(Abundance, 0)) %>% 
  mutate(Fermentation = Fermentation %>%  as.character() %>%  replace_na(., "NA")) %>% 
  mutate(Fermentation = factor(Fermentation, levels = c("NA", 1, 2))) -> pre_gen_grp$data 


pre_gen_grp + 
  facet_grid(Model ~  Model2 + Fermentation + Reactor_Treatment_Dose, scales = "free", space = "free") +
  scale_color_manual(values = c("#31a354", "orange", 
                                "#f03b20"), drop = FALSE, na.value = 'transparent') -> pre_gen_grp


pre_gen_grp
```

<img src="manuscript_Figure1-S1_files/figure-html/unnamed-chunk-20-1.png" style="display: block; margin: auto;" />

```r
before_treat %>%
  filter_taxa(function(x) sum(x > 0) > 0, TRUE) %>% 
  transform_sample_counts(function(x) x/sum(x) * 100) %>%
  subset_taxa(Genus %in% top_gen_group) %>% 
  ps_arrange(Day_of_Treatment) %>% 
  speedyseq::psmelt() %>%
  # group_by(Sample, Day_of_Treatment, Fermentation, Model, Model2, Reactor_Treatment_Dose) %>%
  group_by(sample_names, Model, Model2, Day_of_Treatment , Fermentation,  Reactor_Treatment_Dose) %>% 
  summarise(coverage = sum(Abundance)) %>%
  ungroup() %>%
  mutate(Day_of_Treatment = as.double(Day_of_Treatment)) %>%
  # mutate_at(vars(Day_of_Treatment), ~replace_na(., 0)) %>%
  arrange(Day_of_Treatment) %>%
  # group_by(sample_names, Model, Model2, Day_of_Treatment, Reactor_Treatment_Dose_fermentation) %>% 
  #      summarise(count = length(coverage)) %>% 
  ggpubr::ggbarplot(., x = "Day_of_Treatment", y = "coverage",
                    fill = "grey20", color = "grey5", width = 0.8, add = "segments") + #,
  # ylab = "gdh copies / 10^10 cells \n", xlab = "") +
  # facet_grid(variable ~ ., scales = "free") +
  # ggplot(., aes(x=Day_of_Treatment, y=coverage)) +
  # geom_point() +
  ggpubr::rotate_x_text(60)  +
  theme_light() + 
  facet_grid(Model ~  Model2 + Fermentation + Reactor_Treatment_Dose, scales = "free", space = "free") -> repr_top_gen_donor_grp
```

```
## `summarise()` has grouped output by 'sample_names', 'Model', 'Model2',
## 'Day_of_Treatment', 'Fermentation'. You can override using the `.groups`
## argument.
```

```r
# facet_grid(Model ~  Model2 + Reactor_Treatment_Dose_fermentation,  scales = "fixed", space = "fixed", drop = TRUE) -> toto
# facet_grid(Model ~  Model2 + Fermentation + Reactor_Treatment_Dose, scales = "free", space = "free", drop = FALSE) -> toto
# facet_wrap()

repr_top_gen_donor_grp
```

<img src="manuscript_Figure1-S1_files/figure-html/unnamed-chunk-21-1.png" style="display: block; margin: auto;" />
## gram plots:

histogram of gram+ gram- & intrinsicly resistant or not an annotation color in a heatmap would be cool as well.


```r
"data/processed/16S/16S_based_gram.RDS" %>% 
  here::here() %>% 
  readRDS() -> gram

gram %>% 
  summarise_all(~ sum(is.na(.)))
```

```
## # A tibble: 1 × 9
##     ASV Family Class Phylum Genus Species gram_neg_genus gram_neg_phylum gram_…¹
##   <int>  <int> <int>  <int> <int>   <int>          <int>           <int>   <int>
## 1     0    146    15     11   426    1206            583             188      23
## # … with abbreviated variable name ¹​gram_neg_class
```


```r
dim(gram)
```

```
## [1] 1321    9
```

```r
intersect(gram$ASV,
          tax_table(ps_filtered)[,"Strain"] %>%  as.character() ) %>%  length() 
```

```
## [1] 1321
```

Let's take the Class info here.

## intrinsic:


```r
"data/processed/16S/16S_based_intrinsinc_AMR.RDS" %>% 
  here::here() %>% 
  readRDS() -> intrinsic_AMR
```

Wich level to chose? the one with less NA:


```r
intrinsic_AMR %>% 
  select(starts_with("intrinsic")) %>% 
  summary()
```

```
##  intrinsic_strain_ctx intrinsic_strain_van
##  Mode :logical        Mode :logical       
##  FALSE:588            FALSE:554           
##  TRUE :1              TRUE :35
```

Let's take the Strain info here.


```r
ps_rare_gram_int <- before_treat

# We want sample, reactor..., day_0f_treatment 
ps_rare_gram_int %>% 
  tax_table() %>% 
  data.frame() %>% 
  rownames_to_column("ASVid") %>% 
  left_join(gram %>% dplyr::rename("Strain" = "ASV") %>%  
              select(Strain , gram_neg_class)
  ) %>%  #,
  # by = c("Strain" = "ASV"),
  # suffix = c("", ".y")) %>% 
  select(ASVid:Strain, gram_neg_class) %>% 
  dplyr::rename('gram_stain' = 'gram_neg_class') %>% 
  mutate(gram_stain = replace_na(gram_stain, "unknown-gram")) %>% 
  left_join(., intrinsic_AMR,
            suffix = c("", ".y"),
            by = c("Strain" = "ASV")) %>%
  mutate(intrinsic = ifelse(intrinsic_strain_van == TRUE & intrinsic_strain_ctx == TRUE, "intrinsic-VAN-CTX",  ifelse(intrinsic_strain_van == TRUE, "intrinsic-VAN", ifelse(intrinsic_strain_ctx == TRUE, "intrinsic-CTX", "No-intrinsic")))) %>% 
  select(ASVid:gram_stain,intrinsic) %>%
  mutate(gram_intrinsic = paste0(gram_stain, 
                                 ifelse(is.na(intrinsic), "", paste0("_",intrinsic)))) %>% 
  # dplyr::rename('intrinsic_CTX' = 'intrinsic_strain_ctx',
  # 'intrinsic_VAN' = 'intrinsic_strain_van') %>% 
  column_to_rownames('ASVid') %>% 
  as.matrix() -> tax_table(ps_rare_gram_int)
```

```
## Joining, by = "Strain"
```

```r
# ps_rare_gram_int %>% 
#   phloseq_export_otu_tax() %>% 
#   select(-ASV_sequence) %>% 
#     xlsx::write.xlsx(file = out_xlsx,
#                    append = TRUE,
#                    sheetName = "gram_classification",
#                    showNA = TRUE)
```



```r
require(microViz)

ps_rare_gram_int %>% 
  subset_samples(sample_name %in% c(before_treat_samples)) %>% 
  # subset_samples(Model2 == "Chicken2" & Treatment != "DONOR" & Day_of_Treatment >= -1) %>% 
  # microViz::ps_mutate(Reactor_Treatment_Dose = paste0(Reactor_Treatment_Dose, "_", Fermentation)) %>% 
  # microViz::ps_mutate(Treatment = fct_relevel(Treatment, sample_data(ps_rare)$Treatment %>%  levels())) %>% 
  transform_sample_counts(function(x) x/sum(x) * 100) %>%
  physeq_sel_tax_table(c("gram_stain")) %>% 
  speedyseq::tax_glom("gram_stain") %>% 
  speedyseq::psmelt() %>% 
  arrange(-Day_of_Treatment) %>% 
  ggpubr::ggbarplot(., x = "Day_of_Treatment", y = "Abundance", 
                    fill = "gram_stain", width = 0.8) + #,
  # ylab = "gdh copies / 10^10 cells \n", xlab = "") +
  # facet_grid(variable ~ ., scales = "free") + 
  ggpubr::rotate_x_text(60)  + 
  scale_color_manual(values  =   microViz::distinct_palette(3, pal = "greenArmytage")) +
  scale_fill_manual(values  =   microViz::distinct_palette(3, pal = "greenArmytage")) +
  theme_light() + facet_grid(Model ~  Model2 + Fermentation + Reactor_Treatment_Dose, scales = "free", space = "free_x", drop = TRUE) 
```

```
## Warning: Using an external vector in selections was deprecated in tidyselect 1.1.0.
## ℹ Please use `all_of()` or `any_of()` instead.
##   # Was:
##   data %>% select(tax_sel)
## 
##   # Now:
##   data %>% select(all_of(tax_sel))
## 
## See <https://tidyselect.r-lib.org/reference/faq-external-vector.html>.
```

<img src="manuscript_Figure1-S1_files/figure-html/unnamed-chunk-27-1.png" style="display: block; margin: auto;" />

```r
# scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
#               labels = trans_format("log10", math_format(10^.x))) 
```

## alpha div plots:


```r
before_treat %>% 
  phyloseq_alphas(physeq = ., phylo = FALSE, compute_NTI_PD = TRUE, export_hill = TRUE, return = "long") -> a_div
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

```r
a_div %>% 
  filter(! alphadiversiy %in% c("diversity_coverage", "evenness_pielou", "Hill-d1",  "Hill-d2") ) %>% 
  mutate(alphadiversiy = fct_relevel(alphadiversiy, c("Hill-d0"))) %>%
  ggpubr::ggbarplot(., x = "Day_of_Treatment", y = "value",
                    fill = "grey20", color = "grey5", width = 0.8, add = "segments") + #,
  # ylab = "gdh copies / 10^10 cells \n", xlab = "") +
  # facet_grid(variable ~ ., scales = "free") +
  # ggplot(., aes(x=Day_of_Treatment, y=coverage)) +
  # geom_point() +
  ggpubr::rotate_x_text(60)  +
  theme_light() + 
  facet_grid(alphadiversiy ~  Model2  +  Reactor_Treatment_Dose, scales = "free", space = "free_x") -> alpha_hist
# mutate(Day_of_Treatment = sort(Day_of_Treatment)) %>% 
# ggplot(data = .,
#        mapping = aes(x = "Day_of_Treatment", y = "value")) +
# ggpubr::ggdotchart(., x = "Day_of_Treatment", y = "value", 
#                    fill = "grey20", color = "grey5", width = 0.8, add = "segments") + 

# ggplot(aes(x = Day_of_Treatment, y = value)) +
# geom_point(color = "grey20") +
# ggpubr::rotate_x_text(60)  + 
# theme_light() + facet_grid(alphadiversiy ~  Model2 + Fermentation + Reactor_Treatment_Dose, scales = "free", space = "free_x", drop = TRUE) -> alpha_hist
# scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
#               labels = trans_format("log10", math_format(10^.x))) 

alpha_hist
```

<img src="manuscript_Figure1-S1_files/figure-html/unnamed-chunk-28-1.png" style="display: block; margin: auto;" />

```r
# alpha_hist %>% 
#   export::graph2ppt(append = TRUE, 
#                     width = 317.48031496 * 2.5, 
#                     height = 0.618 * 317.48031496 * 1 , paper = "A4",  scaling = 2,
#                     file = out_pptx)
```

# Export - plots:

## vs per rector


```r
ggpubr::ggarrange(alpha_hist +
                    theme(legend.position = "null") +
                    theme(axis.title.x=element_blank(),
                          axis.text.x=element_blank(),
                          axis.ticks.x=element_blank()),
                  pre_Family_grp + theme(
                    strip.background = element_blank(),
                    strip.text.x = element_blank(),
                    legend.position = "null",
                    axis.title.x=element_blank(),
                    axis.text.x=element_blank(),
                    axis.ticks.x=element_blank()),
                  pre_gen_grp + theme(
                    strip.background = element_blank(),
                    strip.text.x = element_blank(),
                    legend.position = "null",
                    axis.title.x=element_blank(),
                    axis.text.x=element_blank(),
                    axis.ticks.x=element_blank()),
                  repr_top_gen_donor_grp + theme(
                    strip.background = element_blank(),
                    strip.text.x = element_blank(),
                    legend.position = "null",
                    axis.title.x=element_blank(),
                    axis.text.x=element_blank(),
                    axis.ticks.x=element_blank()),
                  pre_speces_group + theme(
                    strip.background = element_blank(),
                    strip.text.x = element_blank(),
                    legend.position = "null",
                    axis.title.x=element_blank(),
                    axis.text.x=element_blank(),
                    axis.ticks.x=element_blank()),
                  repr_top_ASV_pre_speces_group + theme(
                    strip.background = element_blank(),
                    strip.text.x = element_blank()
                  )  + theme(legend.position = "null"),
                  align = "v",
                  ncol = 1,
                  heights = c(0.75, 1, 1.5 ,0.3 , 3.75, 0.6),
                  common.legend = FALSE) -> pre_plots_grp
```

```
## Warning: Graphs cannot be vertically aligned unless the axis parameter is set.
## Placing graphs unaligned.
```

```r
pre_plots_grp
```

<img src="manuscript_Figure1-S1_files/figure-html/unnamed-chunk-29-1.png" style="display: block; margin: auto;" />

```r
pre_plots_grp %>%
  ggsave(file="~/Desktop/16S-chick_Fig1b.svg", plot=., width=8, height=10)

# pre_plots_grp %>% 
#   export::graph2ppt(append = TRUE, 
#                     width = 317.48031496 * 6, 
#                     height = 0.618 * 317.48031496 * 8 , paper = "A2",  scaling = 2,
#                     file = out_pptx)
# 

# alpha_hist$data %>%
#   select(alphadiversiy, sample_name, value,Reactor_Treatment_Dose) %>%
#   pivot_wider(names_from = alphadiversiy, values_from = value) %>%
#   xlsx::write.xlsx(file = out_xlsx,
#                    append = TRUE,
#                    sheetName = "pre_alpha",
#                    showNA = TRUE)

pre_Family_grp$data %>%
  select(Display, Sample, Abundance,Reactor_Treatment_Dose) %>%
  pivot_wider(names_from = Display, values_from = Abundance) %>%
  xlsx::write.xlsx(file = out_xlsx,
                   append = TRUE,
                   sheetName = "pre_family_grp",
                   showNA = TRUE)

pre_gen_grp$data %>%
  select(Display, Sample, Abundance,Reactor_Treatment_Dose) %>%
  pivot_wider(names_from = Display, values_from = Abundance) %>%
  xlsx::write.xlsx(file = out_xlsx,
                   append = TRUE,
                   sheetName = "pre_gen_grp",
                   showNA = TRUE)

repr_top_gen_donor_grp$data %>%
  select(sample_names, coverage,Reactor_Treatment_Dose) %>%
  # pivot_wider(names_from = sample_names, values_from = coverage) %>%
  xlsx::write.xlsx(file = out_xlsx,
                   append = TRUE,
                   sheetName = "pre_gen_rep_grp",
                   showNA = TRUE)


pre_speces_group$data %>%
  select(Display, Sample, Abundance,Reactor_Treatment_Dose) %>%
  pivot_wider(names_from = Display, values_from = Abundance) %>%
  xlsx::write.xlsx(file = out_xlsx,
                   append = TRUE,
                   sheetName = "pre_ASV_grp",
                   showNA = TRUE)


repr_top_ASV_pre_speces_group$data %>%
  select(sample_names, coverage,Reactor_Treatment_Dose) %>%
  # pivot_wider(names_from = sample_names, values_from = coverage) %>%
  xlsx::write.xlsx(file = out_xlsx,
                   append = TRUE,
                   sheetName = "pre_rep_ASV_grp",
                   showNA = TRUE)
```


## top vs feces:


```r
ggpubr::ggarrange(alpha_hist + ylab(NULL) + scale_y_continuous(breaks = c(0, 75, 150),  labels = c(0, 75, 150)) + 
                    theme(legend.position = "null") + 
                    theme(axis.title.x=element_blank(),
                          axis.text.x=element_blank(),
                          axis.ticks.x=element_blank()),
                  pre_Family + ylab(NULL) + theme(
                    strip.background = element_blank(),
                    strip.text.x = element_blank(),
                    legend.position = "null",
                    axis.title.x=element_blank(),
                    axis.text.x=element_blank(),
                    axis.ticks.x=element_blank()),
                  pre_gen + ylab(NULL) + theme(
                    strip.background = element_blank(),
                    strip.text.x = element_blank(),
                    legend.position = "null",
                    axis.title.x=element_blank(),
                    axis.text.x=element_blank(),
                    axis.ticks.x=element_blank()),
                  repr_top_gen_donor_feces + ylab(NULL) + scale_y_continuous(breaks = c(0, 50, 100),  labels = c(0, 50, 100)) + theme(
                    strip.background = element_blank(),
                    strip.text.x = element_blank(),
                    legend.position = "null",
                    axis.title.x=element_blank(),
                    axis.text.x=element_blank(),
                    axis.ticks.x=element_blank()),
                  pre_speces_feces + ylab(NULL) + theme(
                    strip.background = element_blank(),
                    strip.text.x = element_blank(),
                    legend.position = "null",
                    axis.title.x=element_blank(),
                    axis.text.x=element_blank(),
                    axis.ticks.x=element_blank()),
                  repr_top_ASV_donor_feces  + ylab(NULL) + scale_y_continuous(breaks = c(0, 50, 100),  labels = c(0, 50, 100)) + theme(
                    strip.background = element_blank(),
                    strip.text.x = element_blank()
                  )  + theme(legend.position = "null"),
                  align = "v",
                  ncol = 1, 
                  nrow = 6,
                  # heights = c(1, 1, 1 ,0.25 , 2.5, 0.5),
                  heights = c(2, 2.75, 3 , 1 , 6.25, 1.5),
                  common.legend = TRUE) -> pre_plots
```

```
## Warning: Graphs cannot be vertically aligned unless the axis parameter is set.
## Placing graphs unaligned.
```

```r
pre_plots
```

<img src="manuscript_Figure1-S1_files/figure-html/unnamed-chunk-30-1.png" style="display: block; margin: auto;" />

```r
pre_plots %>% 
  ggsave(file="~/Desktop/16S-chick_Fig1.svg", plot=., width=10, height=8)

pre_plots %>% 
  export::graph2ppt(append = TRUE, 
                    width = 317.48031496 * 6, 
                    height = 0.618 * 317.48031496 * 8 , paper = "A3",  scaling = 2,
                    file = out_pptx)
```

```
## Exported graph as ~/Desktop/16S-human_S1.pptx
```

```r
# alpha_hist$data %>% 
#   select(alphadiversiy, sample_name, value,Reactor_Treatment_Dose) %>% 
#   pivot_wider(names_from = alphadiversiy, values_from = value) %>% 
#   xlsx::write.xlsx(file = out_xlsx,
#                    append = TRUE,
#                    sheetName = "pre_alpha",
#                    showNA = TRUE)

pre_Family$data %>% 
  select(Display, Sample, Abundance,Reactor_Treatment_Dose) %>% 
  pivot_wider(names_from = Display, values_from = Abundance) %>% 
  xlsx::write.xlsx(file = out_xlsx,
                   append = TRUE,
                   sheetName = "pre_family",
                   showNA = TRUE)

pre_gen$data %>% 
  select(Display, Sample, Abundance,Reactor_Treatment_Dose) %>% 
  pivot_wider(names_from = Display, values_from = Abundance) %>% 
  xlsx::write.xlsx(file = out_xlsx,
                   append = TRUE,
                   sheetName = "pre_gen",
                   showNA = TRUE)

repr_top_gen_donor_feces$data %>% 
  select(sample_names, coverage,Reactor_Treatment_Dose) %>% 
  # pivot_wider(names_from = sample_names, values_from = coverage) %>% 
  xlsx::write.xlsx(file = out_xlsx,
                   append = TRUE,
                   sheetName = "pre_gen_rep",
                   showNA = TRUE)


pre_speces_feces$data %>% 
  select(Display, Sample, Abundance,Reactor_Treatment_Dose) %>% 
  pivot_wider(names_from = Display, values_from = Abundance) %>% 
  xlsx::write.xlsx(file = out_xlsx,
                   append = TRUE,
                   sheetName = "pre_ASV",
                   showNA = TRUE)


repr_top_ASV_donor_feces$data %>% 
  select(sample_names, coverage,Reactor_Treatment_Dose) %>% 
  # pivot_wider(names_from = sample_names, values_from = coverage) %>% 
  xlsx::write.xlsx(file = out_xlsx,
                   append = TRUE,
                   sheetName = "pre_rep_ASV",
                   showNA = TRUE)
```
two feces plot:




```r
before_treat %>% 
  subset_samples(Treatment == "DONOR" ) %>% 
  filter_taxa(function(x) sum(x > 0) > 0, TRUE) %>% 
  # microViz::ps_mutate(Reactor_Treatment_Dose = paste0(Reactor_Treatment_Dose, "_", Fermentation)) %>% 
  transform_sample_counts(function(x) x/sum(x) * 100) %>%
  microViz::phyloseq_validate() %>%
  microViz::tax_fix() %>% 
  microViz::ps_arrange(Day_of_Treatment) %>% 
  microViz::comp_barplot(
    tax_transform_for_plot = "identity",
    tax_level = "Family",
    label = "Model2", # name an alternative variable to label axes
    n_taxa = 20, # give more taxa unique colours
    palette = microViz::distinct_palette(20, pal = "kelly"),
    merge_other = FALSE, # split the "other" category to display alpha diversity
    bar_width = 0.9, # reduce the bar width to 70% of one row
    bar_outline_colour = "grey5",
    sample_order = "default") +
  ylab("Proportion - %") + xlab( " Donor ") -> p_hist
```

```
## NAs detected in phyloseq tax_table:
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
p_hist + facet_grid(. ~ Reactor_Treatment_Dose, scales = "free", space = "free_x", drop = TRUE) -> p_hist_human1_fam

p_hist_human1_fam  +  theme_light() + theme(legend.position = "none") -> p_hist_human1_fam_no_leg

p_hist_human1_fam_no_leg
```

<img src="manuscript_Figure1-S1_files/figure-html/unnamed-chunk-31-1.png" style="display: block; margin: auto;" />

```r
p_hist_human1_fam %>%  
  ggpubr::get_legend() %>% 
  ggpubr::as_ggplot() -> p_hist_human1_fam_leg
```


```r
hueRank <- "Phylum"
hueRankPlural <- "Phyla"
shadeRank <- "Family"

# Sort phyloseq at lower, and then higher ranks
pseq2 <- 
  before_treat %>% 
  subset_samples(Treatment == "DONOR" ) %>% 
  filter_taxa(function(x) sum(x > 0) > 0, TRUE) %>% 
  # microViz::ps_mutate(Reactor_Treatment_Dose = paste0(Reactor_Treatment_Dose, "_", Fermentation)) %>% 
  transform_sample_counts(function(x) x/sum(x) * 100) %>%
  microViz::phyloseq_validate() %>%
  microViz::tax_fix() %>%
  tax_sort(by = sum, at = shadeRank) %>%
  tax_sort(by = sum, at = hueRank) %>%
  tax_agg(rank = shadeRank)
```

```
## NAs detected in phyloseq tax_table:
## Consider using tax_fix() to make taxa uniquely identifiable
```

```
## Warning in ps_counts(data, warn = TRUE): otu_table of counts is NOT available!
## Available otu_table contains non-zero values that are less than 1

## Warning in ps_counts(data, warn = TRUE): otu_table of counts is NOT available!
## Available otu_table contains non-zero values that are less than 1

## Warning in ps_counts(data, warn = TRUE): otu_table of counts is NOT available!
## Available otu_table contains non-zero values that are less than 1

## Warning in ps_counts(data, warn = TRUE): otu_table of counts is NOT available!
## Available otu_table contains non-zero values that are less than 1
```

```r
# Specify number of hues and shades desired
nHues <- 3 # "Other" phyla will be shades of grey
nShades <- 4 # "Other" families will be the lightest shade of each hue

hierarchicalPalInfo <- data.frame(
  hue = as.vector(tt_get(pseq2)[, hueRank]),
  shade = as.vector(tt_get(pseq2)[, shadeRank]),
  counts = taxa_sums(otu_get(pseq2))
)

hierarchicalPalInfo <- hierarchicalPalInfo %>%
  dplyr::mutate(
    hue = forcats::fct_other(
      f = hue, keep = unique(hue)[seq_len(nHues)],
      other_level = paste("Other", hueRankPlural)
    ),
    nChrHue = nchar(as.character(hue)), padHue = max(nChrHue) - nChrHue
  ) %>%
  dplyr::group_by(hue) %>%
  dplyr::mutate(
    shade = forcats::fct_other(
      f = shade, keep = unique(shade)[seq_len(nShades - 1)],
      other_level = "Other"
    )
  ) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(
    nChrShade = nchar(as.character(shade)), padShade = max(nChrShade) - nChrShade,
    Taxa = paste0(hue, ": ", strrep(" ", padHue), shade, strrep(" ", padShade))
  )


hierarchicalPalMatrix <- matrix(
  data = sapply(
    X = seq(from = 30, to = 75, length.out = nShades),
    FUN = function(l) scales::hue_pal(l = l, h.start = 30)(n = nHues)
  ),
  byrow = TRUE, ncol = nHues
)
hierarchicalPalMatrix <- cbind(hierarchicalPalMatrix, grey.colors(n = nShades))

hierarchicalPal <- hierarchicalPalMatrix %>%
  as.vector() %>%
  setNames(unique(hierarchicalPalInfo$Taxa))

tax_palette_plot(hierarchicalPal) +
  theme(axis.text.y.left = element_text(family = "mono"))
```

<img src="manuscript_Figure1-S1_files/figure-html/unnamed-chunk-32-1.png" style="display: block; margin: auto;" />

```r
pseq2 %>%
  ps_get() %>%
  tax_mutate("Phylum: Family" = hierarchicalPalInfo$Taxa, .keep = "none") %>%
  comp_barplot(
    label = "Model2",     merge_other = FALSE,     tax_transform_for_plot = "identity",
    tax_level = "Phylum: Family", n_taxa = length(hierarchicalPal),
    tax_order = "asis", palette = hierarchicalPal, bar_width = 0.975
  ) +
  # coord_flip() +
  theme(legend.text = element_text(family = "mono")) + ylab("Proportion") -> plot_tmp# for text alignment
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
plot_tmp
```

<img src="manuscript_Figure1-S1_files/figure-html/unnamed-chunk-33-1.png" style="display: block; margin: auto;" />

```r
# plot_tmp + facet_grid(Phylum ~ ., scales = "free", space = "free_x", drop = TRUE) 

plot_tmp  +  theme_light() + theme(legend.position = "none") -> plot_tmp_noleg

plot_tmp_noleg %>% 
  ggsave(file="~/Desktop/16S-chick_FigS1.svg", plot=., width=2, height=3)

plot_tmp_noleg %>% 
  export::graph2ppt(append = TRUE, 
                    width = 317.48031496  / 3, 
                    height = 0.618 * 317.48031496 , paper = "A4",  scaling = 2,
                    file = out_pptx)
```

```
## Exported graph as ~/Desktop/16S-human_S1.pptx
```

```r
plot_tmp %>%  
  ggpubr::get_legend() %>% 
  ggpubr::as_ggplot() -> plot_tmp_leg

plot_tmp_leg %>% 
  export::graph2ppt(append = TRUE, 
                    width = 317.48031496 , 
                    height = 0.618 * 317.48031496  , paper = "A4",  scaling = 2,
                    file = out_pptx)
```

```
## Exported graph as ~/Desktop/16S-human_S1.pptx
```



```r
pseq2 %>%
  ps_get() %>%
  tax_mutate("Phylum: Family" = hierarchicalPalInfo$Taxa, .keep = "none") %>%
  comp_barplot(
    label = "Model2",     merge_other = FALSE,     tax_transform_for_plot = "identity",
    tax_level = "Phylum: Family", n_taxa = length(hierarchicalPal),
    tax_order = "asis", palette = hierarchicalPal, bar_width = 0.975
  ) +
  # coord_flip() +
  theme(legend.text = element_text(family = "mono")) + ylab("Proportion") -> plot_tmp# for text alignment
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
plot_tmp
```

<img src="manuscript_Figure1-S1_files/figure-html/unnamed-chunk-34-1.png" style="display: block; margin: auto;" />

```r
# plot_tmp + facet_grid(Phylum ~ ., scales = "free", space = "free_x", drop = TRUE) 

plot_tmp  +  theme_light() + theme(legend.position = "none") -> plot_tmp_noleg

plot_tmp_noleg %>% 
  ggsave(file="~/Desktop/16S-chick_FigS1.svg", plot=., width=2, height=3)

plot_tmp_noleg %>% 
  export::graph2ppt(append = TRUE, 
                    width = 317.48031496  / 3, 
                    height = 0.618 * 317.48031496 , paper = "A4",  scaling = 2,
                    file = out_pptx)
```

```
## Exported graph as ~/Desktop/16S-human_S1.pptx
```

```r
plot_tmp %>%  
  ggpubr::get_legend() %>% 
  ggpubr::as_ggplot() -> plot_tmp_leg

plot_tmp_leg %>% 
  export::graph2ppt(append = TRUE, 
                    width = 317.48031496 , 
                    height = 0.618 * 317.48031496  , paper = "A4",  scaling = 2,
                    file = out_pptx)
```

```
## Exported graph as ~/Desktop/16S-human_S1.pptx
```



```r
require(tidyverse); require(phyloseq)
min_sample_size = 3738
taxrank = "Class"

url("https://github.com/fconstancias/NRP72-FBT/blob/master/data/processed/16S/16S_working_phyloseq.RDS?raw=true")%>% 
  readRDS()  %>% 
  subset_samples(Enrichment == "NotEnriched") %>%
  subset_samples(Experiment %in% c("Continuous", "Cecum")) %>%
  subset_samples(Day_of_Treatment > -6 & Day_of_Treatment <= 30 | Reactor == "DONOR") %>% 
  subset_samples(Model == "Chicken") %>% 
  subset_samples(Reactor != "CR" & Model2 == "Chicken1" | Model2 == "Chicken2" & Reactor %in% c("TR2", "TR3", "TR4", "TR5", "TR6", "TR7", "CR", "DONOR")) %>% 
  subset_samples(Reactor != "IR") %>%
  filter_taxa(function(x) sum(x > 0) > 0, TRUE) %>% 
  rarefy_even_depth(sample.size = min_sample_size, rngseed = 123) %>% 
  # subset_samples(Reactor == "DONOR" & Model == "Chicken") %>% 
  tax_glom(.,  taxrank = eval(taxrank)) %>% 
  # tax_glom(taxrank = taxrank) %>% 
  transform_sample_counts(function(x) x/sum(x) * 100) -> ps_tmp
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
## 158OTUs were removed because they are no longer 
## present in any sample after random subsampling
```

```
## ...
```

```r
as(tax_table(ps_tmp), "matrix") %>% 
  as.data.frame() %>%
  rownames_to_column('ASV') -> tax

as(otu_table(ps_tmp), "matrix") %>% 
  as.data.frame() %>%
  rownames_to_column('ASV') -> otu

otu %>% 
  left_join(tax) %>% 
  # physeq_glom_rename(taxrank = "Strain", rename_ASV = FALSE) %>% 
  write_tsv(paste0("~/Desktop/",taxrank,"_chicken_donors.tsv"))
```

```
## Joining, by = "ASV"
```

```r
ps_tmp %>% 
  sample_data() %>% 
  data.frame() %>% 
    write_tsv(paste0("~/Desktop/metadata_chicken.tsv"))
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
##  [1] gdtools_0.2.4      microbiome_1.19.1  metagMisc_0.0.4    ampvis2_2.7.29    
##  [5] reshape2_1.4.4     scales_1.2.1       compositions_2.0-4 nlme_3.1-161      
##  [9] microViz_0.9.7     phyloseq_1.42.0    forcats_0.5.2      stringr_1.5.0     
## [13] dplyr_1.0.10       purrr_0.3.5        readr_2.1.3        tidyr_1.2.1       
## [17] tibble_3.1.8       ggplot2_3.4.0      tidyverse_1.3.2   
## 
## loaded via a namespace (and not attached):
##   [1] uuid_1.1-0             readxl_1.4.1           backports_1.4.1       
##   [4] systemfonts_1.0.4      plyr_1.8.8             igraph_1.3.5          
##   [7] lazyeval_0.2.2         splines_4.2.2          fantaxtic_0.1.0       
##  [10] GenomeInfoDb_1.34.3    digest_0.6.31          foreach_1.5.2         
##  [13] htmltools_0.5.4        fansi_1.0.3            magrittr_2.0.3        
##  [16] xlsx_0.6.5             googlesheets4_1.0.1    cluster_2.1.4         
##  [19] openxlsx_4.2.5.1       tzdb_0.3.0             Biostrings_2.66.0     
##  [22] extrafont_0.18         modelr_0.1.10          bayesm_3.1-5          
##  [25] vroom_1.6.0            officer_0.5.0          extrafontdb_1.0       
##  [28] svglite_2.1.0          timechange_0.1.1       colorspace_2.0-3      
##  [31] rvest_1.0.3            ggrepel_0.9.2          textshaping_0.3.6     
##  [34] haven_2.5.1            xfun_0.35              crayon_1.5.2          
##  [37] RCurl_1.98-1.9         jsonlite_1.8.4         survival_3.4-0        
##  [40] iterators_1.0.14       ape_5.6-2              glue_1.6.2            
##  [43] rvg_0.3.0              gtable_0.3.1           gargle_1.2.1          
##  [46] zlibbioc_1.44.0        XVector_0.38.0         car_3.1-1             
##  [49] Rttf2pt1_1.3.11        Rhdf5lib_1.20.0        BiocGenerics_0.44.0   
##  [52] DEoptimR_1.0-11        abind_1.4-5            DBI_1.1.3             
##  [55] rstatix_0.7.1          Rcpp_1.0.9             xtable_1.8-4          
##  [58] viridisLite_0.4.1      bit_4.0.5              stats4_4.2.2          
##  [61] htmlwidgets_1.6.0      httr_1.4.4             RColorBrewer_1.1-3    
##  [64] ellipsis_0.3.2         rJava_1.0-6            pkgconfig_2.0.3       
##  [67] farver_2.1.1           sass_0.4.4             dbplyr_2.2.1          
##  [70] utf8_1.2.2             here_1.0.1             tidyselect_1.2.0      
##  [73] labeling_0.4.2         rlang_1.0.6            munsell_0.5.0         
##  [76] cellranger_1.1.0       tools_4.2.2            cachem_1.0.6          
##  [79] cli_3.4.1              generics_0.1.3         devEMF_4.1-1          
##  [82] ade4_1.7-20            export_0.3.0           broom_1.0.2           
##  [85] evaluate_0.19          biomformat_1.26.0      fastmap_1.1.0         
##  [88] yaml_2.3.6             ragg_1.2.4             bit64_4.0.5           
##  [91] knitr_1.41             fs_1.5.2               zip_2.2.2             
##  [94] robustbase_0.95-0      rgl_0.110.2            xml2_1.3.3            
##  [97] compiler_4.2.2         rstudioapi_0.14        plotly_4.10.1         
## [100] ggsignif_0.6.4         reprex_2.0.2           bslib_0.4.2           
## [103] stringi_1.7.8          highr_0.9              stargazer_5.2.3       
## [106] lattice_0.20-45        Matrix_1.5-3           vegan_2.6-4           
## [109] permute_0.9-7          tensorA_0.36.2         multtest_2.54.0       
## [112] vctrs_0.5.1            pillar_1.8.1           lifecycle_1.0.3       
## [115] rhdf5filters_1.10.0    jquerylib_0.1.4        flextable_0.8.3       
## [118] data.table_1.14.6      cowplot_1.1.1          bitops_1.0-7          
## [121] R6_2.5.1               gridExtra_2.3          IRanges_2.32.0        
## [124] codetools_0.2-18       MASS_7.3-58.1          assertthat_0.2.1      
## [127] xlsxjars_0.6.1         rhdf5_2.42.0           rprojroot_2.0.3       
## [130] withr_2.5.0            S4Vectors_0.36.0       GenomeInfoDbData_1.2.9
## [133] mgcv_1.8-41            parallel_4.2.2         hms_1.1.2             
## [136] grid_4.2.2             speedyseq_0.5.3.9018   rmarkdown_2.19        
## [139] carData_3.0-5          googledrive_2.0.0      Rtsne_0.16            
## [142] ggpubr_0.5.0           base64enc_0.1-3        Biobase_2.56.0        
## [145] lubridate_1.9.0
```
