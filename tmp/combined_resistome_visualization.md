---
title: "Resistome_visualization"
author: "Florentin Constancias"
date: "January 07, 2022"
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
library(microViz)
library(RColorBrewer)
library(vegan)
library(randomcoloR)
options(getClass.msg=FALSE) # https://github.com/epurdom/clusterExperiment/issues/66
#this fixes an error message that pops up because the class 'Annotated' is defined in two different packages
```

#### Source required functions


```r
rm(list = ls())

'%!in%' <- function(x,y)!('%in%'(x,y))

source("https://raw.githubusercontent.com/fconstancias/DivComAnalyses/master/R/phyloseq_taxa_tests.R")
source("https://raw.githubusercontent.com/fconstancias/DivComAnalyses/master/R/phyloseq_normalisation.R")
source("https://raw.githubusercontent.com/fconstancias/DivComAnalyses/master/R/phyloseq_alpha.R")
source("https://raw.githubusercontent.com/fconstancias/DivComAnalyses/master/R/phyloseq_beta.R")
source("https://raw.githubusercontent.com/fconstancias/DivComAnalyses/master/R/phyloseq_heatmap.R")
```

# Import data and define variables:


```r
physeq = readRDS("~/Projects/ETH/Alessia/mobilome/ps_combined.rds")
out_pptx = "~/Projects/ETH/Alessia/mobilome/graph.pptx"


physeq %>%
  tax_table() %>%
  data.frame() %>%
  mutate(Resistance_Mechanism_multi = ifelse(grepl(";", Resistance_Mechanism), "multi-mech", Resistance_Mechanism)) %>% 
  as.matrix() -> tax_table(physeq)


physeq %>% 
  tax_table() %>% 
  data.frame() -> tax_mapping

# physeq %>% 
#   tax_table() %>% 
#   data.frame() %>% 
#   mutate_if(is.character, as.factor) %>% 
#   as.matrix() -> tax_table(physeq)

# 
# physeq %>%
#   tax_table() %>%
#   data.frame() %>%
#   mutate_all(funs(str_replace(., "[^[:alnum:]]", "_"))) %>%
#   mutate_all(funs(str_replace(., "[[:punct:]]", "_"))) %>%
#   as.matrix() -> tax_table(physeq)

physeq %>% 
  otu_table() %>% 
  data.frame() %>% 
  replace(is.na(.), 0)  %>% 
  # mutate_all(funs(str_replace(., "[^[:alnum:]]", ""))) %>%
  mutate_all(as.numeric) %>% 
  as.matrix() %>% 
  otu_table(taxa_are_rows = TRUE) -> otu_table(physeq) 
# df[,'values'] <- as.numeric(df[,'values'])
# otu_table(physeq) <- tmp
```



```r
all.equal(taxa_names(physeq %>%  tax_table()), taxa_names(physeq %>%  otu_table()))
```

```
## [1] TRUE
```

```r
all.equal(taxa_names(physeq), taxa_names(physeq %>%  otu_table()))
```

```
## [1] TRUE
```



```r
physeq_tmp <- physeq
tax_table(physeq_tmp) <- tax_table(physeq_tmp)[,c("Resistance_Mechanism_multi", "Best_Hit_ARO")]
# tax_table(physeq_tmp) <- tax_table(physeq_tmp)[,c("Best_Hit_ARO", "Drug_Class")]

physeq_tmp %>% 
  tax_glom(., taxrank = "Best_Hit_ARO") %>% 
  tax_glom(., taxrank = "Resistance_Mechanism_multi") -> physeq_res_multi

taxa_names(physeq_res_multi) <- tax_table(physeq_res_multi)[,"Resistance_Mechanism_multi"]

physeq_res_multi
```

```
## phyloseq-class experiment-level object
## otu_table()   OTU Table:          [ 6 taxa and 92 samples ]:
## sample_data() Sample Data:        [ 92 samples by 60 sample variables ]:
## tax_table()   Taxonomy Table:     [ 6 taxa by 2 taxonomic ranks ]:
## taxa are rows
```


```r
physeq_tmp <- physeq
tax_table(physeq_tmp) <- tax_table(physeq_tmp)[,c("Drug_Class_multi", "Best_Hit_ARO")]
# tax_table(physeq_tmp) <- tax_table(physeq_tmp)[,c("Best_Hit_ARO", "Drug_Class")]

physeq_tmp %>% 
  tax_glom(., taxrank = "Best_Hit_ARO") %>% 
  tax_glom(., taxrank = "Drug_Class_multi") -> physeq_drug_multi

taxa_names(physeq_drug_multi) <- tax_table(physeq_drug_multi)[,"Drug_Class_multi"]

physeq_drug_multi
```

```
## phyloseq-class experiment-level object
## otu_table()   OTU Table:          [ 19 taxa and 92 samples ]:
## sample_data() Sample Data:        [ 92 samples by 60 sample variables ]:
## tax_table()   Taxonomy Table:     [ 19 taxa by 2 taxonomic ranks ]:
## taxa are rows
```



```r
physeq_tmp <- physeq
tax_table(physeq_tmp) <- tax_table(physeq_tmp)[,c("Drug_Class", "Best_Hit_ARO")]
# tax_table(physeq_tmp) <- tax_table(physeq_tmp)[,c("Best_Hit_ARO", "Drug_Class")]

physeq_tmp %>% 
  tax_glom(., taxrank = "Best_Hit_ARO") %>% 
  tax_glom(., taxrank = "Drug_Class") -> physeq_drug

taxa_names(physeq_drug) <- tax_table(physeq_drug)[,"Drug_Class"]

physeq_drug
```

```
## phyloseq-class experiment-level object
## otu_table()   OTU Table:          [ 43 taxa and 92 samples ]:
## sample_data() Sample Data:        [ 92 samples by 60 sample variables ]:
## tax_table()   Taxonomy Table:     [ 43 taxa by 2 taxonomic ranks ]:
## taxa are rows
```


```r
physeq_tmp <- physeq

tax_table(physeq_tmp) <- tax_table(physeq_tmp)[,c("Best_Hit_ARO", "Drug_Class")]

physeq_tmp %>% 
  tax_glom(., taxrank = "Best_Hit_ARO") -> physeq_AMRgn

taxa_names(physeq_AMRgn) <- tax_table(physeq_AMRgn)[,"Best_Hit_ARO"]

physeq_AMRgn
```

```
## phyloseq-class experiment-level object
## otu_table()   OTU Table:          [ 142 taxa and 92 samples ]:
## sample_data() Sample Data:        [ 92 samples by 60 sample variables ]:
## tax_table()   Taxonomy Table:     [ 142 taxa by 2 taxonomic ranks ]:
## taxa are rows
```


<!-- Since we somehow can't use tax_glom, we have to do it manually: -->

<!-- ```{r} -->
<!-- as(tax_table(physeq), "matrix") %>%  -->
<!--   as.data.frame() %>% -->
<!--   rownames_to_column('ASV') -> tax_table_ps_combined -->

<!-- as(otu_table(physeq), "matrix") %>%  -->
<!--   as.data.frame() %>% -->
<!--   rownames_to_column('ASV') -> otu_table_ps_combined -->

<!-- physeq %>%  -->
<!--   sample_data() %>%  -->
<!--   data.frame() %>%  -->
<!--   rownames_to_column('sample_id') ->  meta_ps_combined -->

<!-- otu_table_ps_combined %>% -->
<!--   # left_join(refseq_df, by = 'ASV') %>% -->
<!--   left_join(tax_table_ps_combined, by = c("ASV" = "ASV")) %>% -->
<!--   dplyr::select(ASV, everything()) -> merged_ps_combined -->
<!-- ``` -->


<!-- ```{r} -->
<!-- merged_ps_combined %>%  -->
<!--   pivot_longer(cols = sample_names(physeq), -->
<!--                names_to = "sample") %>%  -->
<!--   group_by(sample, Best_Hit_ARO) %>%  -->
<!--   summarise(Best_Hit_ARO_value = sum(value)) %>%  -->
<!--   replace(is.na(.), 0) %>%  -->
<!--   dplyr::ungroup() %>%  -->
<!--   pivot_wider(names_from = "sample", -->
<!--               values_from =  Best_Hit_ARO_value ) -> gene_df -->

<!-- gene_df -->
<!-- ``` -->


Define colors for Drug classes:


```r
set.seed(123)
treat_col <- ggpubr::get_palette(palette = "lancet",
                           k = get_variable(physeq_drug_multi, "Treatment") %>%  
                              levels() %>% 
                             length())

names(treat_col) <- get_variable(physeq_drug_multi, "Treatment") %>%  
                              levels()
```


Define colors for Treatment - no matter concentration / fermentation / reactors:


```r
myCol <- viridis::viridis(n = get_variable(physeq_drug_multi, "Treatment") %>%  
                              levels() %>% 
                             length() + 20)
myCol
```

```
##  [1] "#440154FF" "#470F62FF" "#481C6EFF" "#482878FF" "#463480FF" "#423F85FF"
##  [7] "#3E4A89FF" "#3A548CFF" "#355F8DFF" "#31688EFF" "#2D718EFF" "#297A8EFF"
## [13] "#26828EFF" "#228C8DFF" "#1F948CFF" "#1F9E89FF" "#21A685FF" "#29AF7FFF"
## [19] "#35B779FF" "#45BF70FF" "#58C765FF" "#6DCD59FF" "#83D44CFF" "#9BD93CFF"
## [25] "#B4DE2CFF" "#CDE11DFF" "#E6E419FF" "#FDE725FF"
```

```r
scales::show_col(myCol,
                 cex_label = 0.5)
```

![](combined_resistome_visualization_files/figure-html/unnamed-chunk-9-1.png)<!-- -->


```r
myCol <- ggpubr::get_palette(palette = "npg",
                             k = get_variable(physeq_drug_multi, "Treatment") %>%  
                              levels() %>% 
                             length() + 3)
myCol
```

```
##  [1] "#E64B35" "#5CAFC5" "#0FA596" "#2A6A87" "#A97E82" "#BB9599" "#89AAB9"
##  [8] "#A79287" "#C9130E" "#83664E" "#B09C85"
```

```r
scales::show_col(myCol,
                 cex_label = 0.5)
```

![](combined_resistome_visualization_files/figure-html/unnamed-chunk-10-1.png)<!-- -->



```r
myCol <- ggpubr::get_palette(palette = "jama",
                             k = get_variable(physeq_drug_multi, "Treatment") %>%  
                              levels() %>% 
                             length() + 6)
myCol
```

```
##  [1] "#374E55" "#846C4D" "#D28A45" "#89957B" "#229EBE" "#3685A8" "#885B66"
##  [8] "#A45F57" "#8A8F7D" "#76A397" "#6F8198" "#6B6695" "#756F80" "#80796B"
```

```r
scales::show_col(myCol,
                 cex_label = 0.5)
```

![](combined_resistome_visualization_files/figure-html/unnamed-chunk-11-1.png)<!-- -->


```r
myCol <- ggpubr::get_palette(palette = "lancet",
                             k = get_variable(physeq_drug_multi, "Treatment") %>%  
                              levels() %>% 
                             length() + 5)
myCol
```

```
##  [1] "#00468B" "#9E172E" "#B43C15" "#42B540" "#16A28D" "#3085AD" "#925E9F"
##  [8] "#D99395" "#E2746E" "#AD002A" "#AD7987" "#7C8181" "#1B1919"
```

```r
scales::show_col(myCol,
                 cex_label = 0.5)
```

![](combined_resistome_visualization_files/figure-html/unnamed-chunk-12-1.png)<!-- -->


```r
treat_col = get_variable(physeq_drug_multi, "Treatment") %>%  
                              levels()


# treat_col <- c("#1B1919", "#7C8181", "#00468B", "#2589AE","#029AAF", "#463161", "#7C66A2", "#ad99cf")

treat_col <- c("#1B1919", "#7C8181", "#00468B", "#3085AD", "#42B540", "#9E172E", "#B43C15", "#AD7987")


names(treat_col) <- get_variable(physeq_drug_multi, "Treatment") %>%  
                              levels()

scales::show_col(treat_col,
                 cex_label = 0.5)
```

![](combined_resistome_visualization_files/figure-html/unnamed-chunk-13-1.png)<!-- -->


```r
set.seed(123)
col <- distinctColorPalette(ntaxa(physeq_drug_multi))

names(col) <- taxa_names(physeq_drug_multi)
```



```r
physeq_AMRgn %>%
  microbiome::transform(transform = "log10p") %>% 
  speedyseq::psmelt() %>% 
  left_join(tax_mapping,
            by = c("Best_Hit_ARO" = "Best_Hit_ARO"),
            suffix = c("_x", "")) %>% 
  mutate(Abundance = na_if(Abundance, 0)) %>%
  mutate(logTPM = log10(Abundance + 1)) %>%
  separate_rows(., Drug_Class,sep = '; ',convert = TRUE) %>% 
  ggplot(mapping = aes(x = as.factor(Day_of_Treatment),
                       y = Best_Hit_ARO ,
                       fill = logTPM,
  )) +
  facet_grid(Drug_Class ~ Model + Treatment + Fermentation + Antibiotic_mg.mL , scales = "free", space = "free") + #,  labeller = labeller(Drug_Class_multi = label_wrap_gen(width = 15)))  +
  # facet_grid(. ~ Reactor_Treatment, scales = "free", space = "free", drop = TRUE,labeller = labeller(Drug.Class = label_wrap_gen(width = 15))) +
  #we adjust the vertical labels such that it is printed correctly. 
  geom_tile() +
  theme_light()  +
  # scale_fill_viridis_c(name = "TPM",
  #                      na.value = "transparent", trans = scales::pseudo_log_trans(sigma = 0.1),
  #                      # trans =  scales::pseudo_log_trans(sigma = 1),
  #                      breaks = c(1, 10, 50, 100, 400), labels = c(1, 10, 50, 100, 400),
  #                      limits = c(0, 500)) + #need to make sure that we include max value of tpm as upper limit
  scale_fill_viridis_c(na.value = "transparent") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        strip.text.x = element_text(size=8),
        #change the size of the label and the angle so the text is printed horizontally
        strip.text.y = element_text(size = 7,angle=0)) -> perma_bbuble

perma_bbuble +
  theme(legend.position = "bottom") -> perma_bbuble 

perma_bbuble +
  scale_y_discrete(label = function(x) stringr::str_trunc(x, 34, side = "center")) ->perma_bbuble

perma_bbuble
```

![](combined_resistome_visualization_files/figure-html/unnamed-chunk-15-1.png)<!-- -->

```r
perma_bbuble %>%
  export::graph2ppt(append = TRUE, width = 20, height = 50,
                    file = out_pptx)
```

```
## Exported graph as ~/Projects/ETH/Alessia/mobilome/graph.pptx
```



```r
physeq_AMRgn %>%
  microbiome::transform(transform = "log10p") %>% 
  speedyseq::psmelt() %>% 
  left_join(tax_mapping,
            by = c("Best_Hit_ARO" = "Best_Hit_ARO"),
            suffix = c("_x", "")) %>% 
  mutate(Best_Hit_ARO = fct_reorder(Best_Hit_ARO, Abundance, .desc=FALSE)) %>% 
  mutate(Abundance = na_if(Abundance, 0)) %>%
  mutate(logTPM = log10(Abundance + 1))  %>% 
  # df$Item = factor(df$Item, levels = a_lot$Item[order(a_lot$percent)])

# tmp$variable <- factor(melted$variable, levels=rev(levels(melted$variable)))
  # mutate(Best_Hit_ARO =  fct_reorder(Best_Hit_ARO, logTPM)) %>% 
  # mutate(Best_Hit_ARO =  fct_relevel(Best_Hit_ARO, logTPM)) %>% 
  # separate_rows(., Drug_Class,sep = '; ',convert = TRUE) %>% 
  ggplot(mapping = aes(x = as.factor(Day_of_Treatment),
                       y = Best_Hit_ARO,
                       fill = logTPM,
  )) +
  facet_grid(.  ~ Model + Treatment + Fermentation + Antibiotic_mg.mL , scales = "free", space = "free") + #,  labeller = labeller(Drug_Class_multi = label_wrap_gen(width = 15)))  +
  # facet_grid(. ~ Reactor_Treatment, scales = "free", space = "free", drop = TRUE,labeller = labeller(Drug.Class = label_wrap_gen(width = 15))) +
  #we adjust the vertical labels such that it is printed correctly. 
  geom_tile() +
  theme_light()  +
  # scale_fill_viridis_c(name = "TPM",
  #                      na.value = "transparent", trans = scales::pseudo_log_trans(sigma = 0.1),
  #                      # trans =  scales::pseudo_log_trans(sigma = 1),
  #                      breaks = c(1, 10, 50, 100, 400), labels = c(1, 10, 50, 100, 400),
  #                      limits = c(0, 500)) + #need to make sure that we include max value of tpm as upper limit
  scale_fill_viridis_c(na.value = "transparent") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        strip.text.x = element_text(size=8),
        #change the size of the label and the angle so the text is printed horizontally
        strip.text.y = element_text(size = 7,angle=0)) -> perma_bbuble

perma_bbuble +
  theme(legend.position = "bottom") +
  scale_y_discrete(label = function(x) stringr::str_trunc(x, 34, side = "center"))  -> perma_bbuble

perma_bbuble
```

![](combined_resistome_visualization_files/figure-html/unnamed-chunk-16-1.png)<!-- -->

```r
perma_bbuble %>%
  export::graph2ppt(append = TRUE, width = 20, height = 28,
                    file = out_pptx)
```

```
## Exported graph as ~/Projects/ETH/Alessia/mobilome/graph.pptx
```



```r
physeq_AMRgn %>%
  microbiome::transform(transform = "log10p") %>% 
  speedyseq::psmelt() %>% 
  left_join(tax_mapping,
            by = c("Best_Hit_ARO" = "Best_Hit_ARO"),
            suffix = c("_x", "")) %>% 
  mutate(Best_Hit_ARO = fct_reorder(Best_Hit_ARO, Abundance, .desc=FALSE)) %>%
  mutate(Abundance = na_if(Abundance, 0)) %>%
  mutate(logTPM = log10(Abundance + 1))  %>% 
  # df$Item = factor(df$Item, levels = a_lot$Item[order(a_lot$percent)])

# tmp$variable <- factor(melted$variable, levels=rev(levels(melted$variable)))
  # mutate(Best_Hit_ARO =  fct_reorder(Best_Hit_ARO, logTPM)) %>% 
  # mutate(Best_Hit_ARO =  fct_relevel(Best_Hit_ARO, logTPM)) %>% 
  # separate_rows(., Drug_Class,sep = '; ',convert = TRUE) %>% 
  ggplot(mapping = aes(x = as.factor(Day_of_Treatment),
                       y = Best_Hit_ARO,
                       fill = logTPM,
  )) +
  facet_grid(Drug_Class_multi  ~  Model + Treatment + Fermentation + Antibiotic_mg.mL , scales = "free", space = "free") + #,  labeller = labeller(Drug_Class_multi = label_wrap_gen(width = 15)))  +
  # facet_grid(. ~ Reactor_Treatment, scales = "free", space = "free", drop = TRUE,labeller = labeller(Drug.Class = label_wrap_gen(width = 15))) +
  #we adjust the vertical labels such that it is printed correctly. 
  geom_tile() +
  theme_light()  +
  # scale_fill_viridis_c(name = "TPM",
  #                      na.value = "transparent", trans = scales::pseudo_log_trans(sigma = 0.1),
  #                      # trans =  scales::pseudo_log_trans(sigma = 1),
  #                      breaks = c(1, 10, 50, 100, 400), labels = c(1, 10, 50, 100, 400),
  #                      limits = c(0, 500)) + #need to make sure that we include max value of tpm as upper limit
  scale_fill_viridis_c(na.value = "transparent") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        strip.text.x = element_text(size=8),
        #change the size of the label and the angle so the text is printed horizontally
        strip.text.y = element_text(size = 7,angle=0)) -> perma_bbuble

perma_bbuble +
  theme(legend.position = "bottom") +
  scale_y_discrete(label = function(x) stringr::str_trunc(x, 34, side = "center"))  -> perma_bbuble

perma_bbuble
```

![](combined_resistome_visualization_files/figure-html/unnamed-chunk-17-1.png)<!-- -->

```r
perma_bbuble %>%
  export::graph2ppt(append = TRUE, width = 20, height = 28,
                    file = out_pptx)
```

```
## Exported graph as ~/Projects/ETH/Alessia/mobilome/graph.pptx
```

```r
# perma_bbuble + facet_null() +  facet_grid(.  ~  Model + Treatment + Fermentation + Antibiotic_mg.mL , scales = "free", space = "free") +
#  ggside::geom_ysidetile(aes(x = 0, yfill = col[perma_bbuble$data$Drug_Class_multi]), show.legend = TRUE) -> test_annot
# 
# 
# test_annot
```

pheatmap test:

all

```r
physeq_AMRgn %>%
  speedyseq::tax_glom("Best_Hit_ARO") %>% 
  microbiome::transform(transform = "log10p") %>% 
  speedyseq::psmelt() %>% 
  # left_join(tax_mapping,
  #           by = c("Best_Hit_ARO" = "Best_Hit_ARO"),
  #           suffix = c("_x", "")) %>% 
  # mutate(Sample = fct_reorder(Sample, Model,Treatment,Fermentation,Antibiotic_mg.mL, Day_of_Treatment )) %>%
  mutate(Abundance = na_if(Abundance, 0)) %>%
  mutate(logTPM = log10(Abundance + 1)) -> pheatmap_df


pheatmap_df %>% 
  # replace(is.na(.), 0)  %>% 
  # select(Sample, Best_Hit_ARO, logTPM) %>% 
  group_by(Sample, Best_Hit_ARO) %>%
  pivot_wider(names_from =  Sample, Best_Hit_ARO,
              values_from = logTPM) %>% 
  column_to_rownames("Best_Hit_ARO") -> pheatmap_samples
  
pheatmap_samples %>% 
  replace(is.na(.), 0)  %>% 
  t() %>% 
  vegdist(na.rm = TRUE) %>% 
  hclust() -> clust_t

pheatmap_samples %>% 
    replace(is.na(.), 0)  %>%
    vegdist(na.rm = TRUE) %>% 
  # cor() %>% 
  hclust() -> clust

 pheatmap_samples %>% 
  pheatmap::pheatmap(cluster_rows = clust ,
                     cluster_cols = clust_t,
                     na_col = "transparent",
                     annotation_row = tax_mapping %>% rownames_to_column('tmp_col') %>% 
                       filter(Best_Hit_ARO %in% taxa_names(physeq_AMRgn)) %>% 
                       distinct(Best_Hit_ARO, .keep_all = TRUE) %>% 
                       column_to_rownames("Best_Hit_ARO") %>% 
                       select(Drug_Class_multi, Resistance_Mechanism_multi),
                     annotation_colors = list(Drug_Class_multi = col,
                                              Treatment = treat_col),
                     annotation_col = physeq_AMRgn %>% 
                       sample_data() %>% 
                       data.frame() %>% 
                       select(Model,Treatment,Fermentation,Antibiotic_mg.mL, Day_of_Treatment)) -> p_pheat
 
 
 p_pheat %>% 
   print()
```

![](combined_resistome_visualization_files/figure-html/unnamed-chunk-18-1.png)<!-- -->

```r
 p_pheat %>%
  export::graph2ppt(append = TRUE, width = 30, height = 30,
                    file = out_pptx)
```

```
## Exported graph as ~/Projects/ETH/Alessia/mobilome/graph.pptx
```

chick

```r
physeq_AMRgn %>%
  speedyseq::tax_glom("Best_Hit_ARO") %>% 
  subset_samples(Model == "Chicken") %>% 
  filter_taxa(function(x) sum(x > 0) > 0, TRUE) -> ps_tmp

ps_tmp %>% 
  microbiome::transform(transform = "log10p") %>% 
  speedyseq::psmelt() %>% 
  # left_join(tax_mapping,
  #           by = c("Best_Hit_ARO" = "Best_Hit_ARO"),
  #           suffix = c("_x", "")) %>% 
  # mutate(Sample = fct_reorder(Sample, Model,Treatment,Fermentation,Antibiotic_mg.mL, Day_of_Treatment )) %>%
  mutate(Abundance = na_if(Abundance, 0)) %>%
  mutate(logTPM = log10(Abundance + 1)) -> pheatmap_df


pheatmap_df %>% 
  # replace(is.na(.), 0)  %>% 
  # select(Sample, Best_Hit_ARO, logTPM) %>% 
  group_by(Sample, Best_Hit_ARO) %>%
  pivot_wider(names_from =  Sample, Best_Hit_ARO,
              values_from = logTPM) %>% 
  column_to_rownames("Best_Hit_ARO") -> pheatmap_samples
  
pheatmap_samples %>% 
  replace(is.na(.), 0)  %>% 
  t() %>% 
  vegdist(na.rm = TRUE) %>% 
  hclust() -> clust_t

# pheatmap_samples %>% 
#     replace(is.na(.), 0)  %>%
#     vegdist(na.rm = TRUE) %>% 
#   # cor() %>% 
#   hclust() -> clust

 pheatmap_samples %>% 
  pheatmap::pheatmap(cluster_rows = FALSE ,
                     cluster_cols = clust_t,
                     na_col = "transparent",
                     annotation_row = tax_mapping %>% rownames_to_column('tmp_col') %>% 
                       filter(Best_Hit_ARO %in% taxa_names(physeq_AMRgn)) %>% 
                       distinct(Best_Hit_ARO, .keep_all = TRUE) %>% 
                       column_to_rownames("Best_Hit_ARO") %>% 
                       select(Drug_Class_multi, Resistance_Mechanism_multi),
                     annotation_colors = list(Drug_Class_multi = col)) -> p_pheat #,
                     # annotation_col = ps_tmp %>% 
                     #   sample_data() %>% 
                     #   data.frame() %>% 
                     #   select(Model,Treatment,Fermentation,Antibiotic_mg.mL, Day_of_Treatment)) -> p_pheat
 
 
 p_pheat %>% 
   print()
```

![](combined_resistome_visualization_files/figure-html/unnamed-chunk-19-1.png)<!-- -->

```r
  p_pheat %>%
  export::graph2ppt(append = TRUE, width = 20, height = 20,
                    file = out_pptx)
```

```
## Exported graph as ~/Projects/ETH/Alessia/mobilome/graph.pptx
```

Human

```r
physeq_AMRgn %>%
  speedyseq::tax_glom("Best_Hit_ARO") %>% 
  subset_samples(Model == "Human") %>% 
  filter_taxa(function(x) sum(x > 0) > 0, TRUE) -> ps_tmp

ps_tmp %>% 
  microbiome::transform(transform = "log10p") %>% 
  speedyseq::psmelt() %>% 
  # left_join(tax_mapping,
  #           by = c("Best_Hit_ARO" = "Best_Hit_ARO"),
  #           suffix = c("_x", "")) %>% 
  # mutate(Sample = fct_reorder(Sample, Model,Treatment,Fermentation,Antibiotic_mg.mL, Day_of_Treatment )) %>%
  mutate(Abundance = na_if(Abundance, 0)) %>%
  mutate(logTPM = log10(Abundance + 1)) -> pheatmap_df


pheatmap_df %>% 
  # replace(is.na(.), 0)  %>% 
  # select(Sample, Best_Hit_ARO, logTPM) %>% 
  group_by(Sample, Best_Hit_ARO) %>%
  pivot_wider(names_from =  Sample, Best_Hit_ARO,
              values_from = logTPM) %>% 
  column_to_rownames("Best_Hit_ARO") -> pheatmap_samples
  
pheatmap_samples %>% 
  replace(is.na(.), 0)  %>% 
  t() %>% 
  vegdist(na.rm = TRUE) %>% 
  hclust() -> clust_t

pheatmap_samples %>%
    replace(is.na(.), 0)  %>%
    vegdist(na.rm = TRUE) %>%
  # cor() %>%
  hclust() -> clust

 pheatmap_samples %>% 
  pheatmap::pheatmap(cluster_rows = clust ,
                     cluster_cols = clust_t,
                     na_col = "transparent",
                     annotation_row = tax_mapping %>% rownames_to_column('tmp_col') %>% 
                       filter(Best_Hit_ARO %in% taxa_names(physeq_AMRgn)) %>% 
                       distinct(Best_Hit_ARO, .keep_all = TRUE) %>% 
                       column_to_rownames("Best_Hit_ARO") %>% 
                       select(Drug_Class_multi, Resistance_Mechanism_multi),
                     annotation_colors = list(Drug_Class_multi = col,
                                              Treatment = treat_col),
                     annotation_col = ps_tmp %>%
                       sample_data() %>%
                       data.frame() %>%
                       select(Model,Treatment,Fermentation,Antibiotic_mg.mL, Day_of_Treatment)) -> p_pheat
 
 
 p_pheat %>% 
   print()
```

![](combined_resistome_visualization_files/figure-html/unnamed-chunk-20-1.png)<!-- -->

```r
   p_pheat %>%
  export::graph2ppt(append = TRUE, width = 20, height = 20,
                    file = out_pptx)
```

```
## Exported graph as ~/Projects/ETH/Alessia/mobilome/graph.pptx
```



```r
# see Sneha's & Hannah's code for "plot.df.sep <- separate_rows(plot.df, Drug.Class,sep = '; ',convert = TRUE)"

physeq_AMRgn %>%
  microbiome::transform(transform = "log10p") %>% 
  speedyseq::psmelt() %>% 
    left_join(tax_mapping,
            by = c("Best_Hit_ARO" = "Best_Hit_ARO"),
            suffix = c("_x", "")) %>%
    mutate(Abundance = na_if(Abundance, 0)) %>%
  mutate(logTPM = log10(Abundance + 1))  %>% 
  separate_rows(., Resistance_Mechanism,sep = '; ',convert = TRUE) %>% 
  ggplot(aes(x = as.factor(Day_of_Treatment), y = Abundance, fill = Resistance_Mechanism)) +   
  geom_bar( stat = "identity", colour="black", size=0.025) +
  facet_grid(Drug_Class_multi  ~ Model +  Treatment + Fermentation + Antibiotic_mg.mL , scales = "free", space = "free_x") +
  ggtitle("ARG Type Abundance per Sample") +
  xlab("Day (Treatment)")  +
  ylab("log10 TPM") + 
  theme_light() -> p_mec

p_mec %>% 
  ggpubr::set_palette("jco") -> p_mec

p_mec
```

```
## Warning: Removed 9564 rows containing missing values (position_stack).
```

![](combined_resistome_visualization_files/figure-html/unnamed-chunk-21-1.png)<!-- -->

```r
p_mec %>%
  export::graph2ppt(append = TRUE, width = 20, height = 18,
                    file = out_pptx)
```

```
## Warning: Removed 9564 rows containing missing values (position_stack).
```

```
## Exported graph as ~/Projects/ETH/Alessia/mobilome/graph.pptx
```


```r
# see Sneha's & Hannah's code for "plot.df.sep <- separate_rows(plot.df, Drug.Class,sep = '; ',convert = TRUE)"

physeq_AMRgn %>%
  microbiome::transform(transform = "log10p") %>% 
  speedyseq::psmelt() %>% 
    left_join(tax_mapping,
            by = c("Best_Hit_ARO" = "Best_Hit_ARO"),
            suffix = c("_x", "")) %>%
    mutate(Abundance = na_if(Abundance, 0)) %>%
  mutate(logTPM = log10(Abundance + 1))  %>% 
  separate_rows(., Resistance_Mechanism,sep = '; ',convert = TRUE) %>% 
  ggplot(aes(x = as.factor(Day_of_Treatment), y = Abundance, fill = Resistance_Mechanism_multi)) +   
  geom_bar( stat = "identity", colour="black", size=0.025) +
  facet_grid(.  ~ Model +  Treatment + Fermentation + Antibiotic_mg.mL , scales = "free", space = "free_x") +
  ggtitle("ARG Type Abundance per Sample") +
  xlab("Day (Treatment)")  +
  ylab("log10 TPM") + 
  theme_light() -> p_mec

p_mec %>% 
  ggpubr::set_palette("jco") -> p_mec

p_mec
```

```
## Warning: Removed 9564 rows containing missing values (position_stack).
```

![](combined_resistome_visualization_files/figure-html/unnamed-chunk-22-1.png)<!-- -->

```r
p_mec %>%
  export::graph2ppt(append = TRUE, width = 18, height = 6,
                    file = out_pptx)
```

```
## Warning: Removed 9564 rows containing missing values (position_stack).
```

```
## Exported graph as ~/Projects/ETH/Alessia/mobilome/graph.pptx
```



```r
# see Sneha's & Hannah's code for "plot.df.sep <- separate_rows(plot.df, Drug.Class,sep = '; ',convert = TRUE)"

physeq_AMRgn %>%
  microbiome::transform(transform = "log10p") %>% 
  speedyseq::psmelt() %>% 
    left_join(tax_mapping,
            by = c("Best_Hit_ARO" = "Best_Hit_ARO"),
            suffix = c("_x", "")) %>%
    mutate(Abundance = na_if(Abundance, 0)) %>%
  mutate(logTPM = log10(Abundance + 1))  %>% 
  # separate_rows(., Resistance_Mechanism,sep = '; ',convert = TRUE) %>% 
  ggplot(aes(x = as.factor(Day_of_Treatment), y = Abundance, fill = Resistance_Mechanism_multi)) +   
  geom_bar( stat = "identity", colour="black", size=0.025) +
  facet_grid(.  ~ Model +  Treatment + Fermentation + Antibiotic_mg.mL , scales = "free", space = "free_x") +
  ggtitle("ARG Type Abundance per Sample") +
  xlab("Day (Treatment)")  +
  ylab("log10 TPM") + 
  theme_light() -> p_mec

p_mec %>% 
  ggpubr::set_palette("jco") -> p_mec

p_mec
```

```
## Warning: Removed 9108 rows containing missing values (position_stack).
```

![](combined_resistome_visualization_files/figure-html/unnamed-chunk-23-1.png)<!-- -->

```r
p_mec %>%
  export::graph2ppt(append = TRUE, width = 18, height = 6,
                    file = out_pptx)
```

```
## Warning: Removed 9108 rows containing missing values (position_stack).
```

```
## Exported graph as ~/Projects/ETH/Alessia/mobilome/graph.pptx
```



```r
# see Sneha's & Hannah's code for "plot.df.sep <- separate_rows(plot.df, Drug.Class,sep = '; ',convert = TRUE)"

physeq_AMRgn %>%
  microbiome::transform(transform = "log10p") %>% 
  speedyseq::psmelt() %>% 
    left_join(tax_mapping,
            by = c("Best_Hit_ARO" = "Best_Hit_ARO"),
            suffix = c("_x", "")) %>%
    mutate(Abundance = na_if(Abundance, 0)) %>%
  mutate(logTPM = log10(Abundance + 1))  %>% 
  separate_rows(., Resistance_Mechanism,sep = '; ',convert = TRUE) %>%
  ggplot(aes(x = as.factor(Day_of_Treatment), y = Abundance, fill = Resistance_Mechanism)) +   
  geom_bar( stat = "identity", position = "fill", colour="black", size=0.025) +
  facet_grid(.  ~ Model +  Treatment + Fermentation + Antibiotic_mg.mL , scales = "free", space = "free_x") +
  ggtitle("ARG Type Abundance per Sample") +
  xlab("Day (Treatment)")  +
  ylab("Proportion") + 
  theme_light() -> p_mec

p_mec %>% 
  ggpubr::set_palette("jco") -> p_mec

p_mec
```

```
## Warning: Removed 9564 rows containing missing values (position_stack).
```

![](combined_resistome_visualization_files/figure-html/unnamed-chunk-24-1.png)<!-- -->

```r
p_mec %>%
  export::graph2ppt(append = TRUE, width = 18, height = 6,
                    file = out_pptx)
```

```
## Warning: Removed 9564 rows containing missing values (position_stack).
```

```
## Exported graph as ~/Projects/ETH/Alessia/mobilome/graph.pptx
```


```r
# https://newbedev.com/how-to-control-ordering-of-stacked-bar-chart-using-identity-on-ggplot2
physeq_AMRgn %>%
  microbiome::transform(transform = "log10p") %>% 
  speedyseq::psmelt() %>% 
    left_join(tax_mapping,
            by = c("Best_Hit_ARO" = "Best_Hit_ARO"),
            suffix = c("_x", "")) %>%
    mutate(Abundance = na_if(Abundance, 0)) %>%
  mutate(logTPM = log10(Abundance + 1))  %>% 
  # separate_rows(., Resistance_Mechanism,sep = '; ',convert = TRUE) %>% 
  ggplot(aes(x = as.factor(Day_of_Treatment), y = Abundance, fill = Resistance_Mechanism_multi)) +   
  geom_bar( stat = "identity", position = "fill", colour="black", size=0.00125) +
  facet_grid(.  ~ Model +  Treatment + Fermentation + Antibiotic_mg.mL , scales = "free", space = "free_x") +
  ggtitle("ARG Type Abundance per Sample") +
  xlab("Day (Treatment)")  +
  ylab("Proportion") + 
  theme_light() -> p_mec

p_mec %>% 
  ggpubr::set_palette("jco") -> p_mec

p_mec
```

```
## Warning: Removed 9108 rows containing missing values (position_stack).
```

![](combined_resistome_visualization_files/figure-html/unnamed-chunk-25-1.png)<!-- -->

```r
p_mec %>%
  export::graph2ppt(append = TRUE, width = 18, height = 6,
                    file = out_pptx)
```

```
## Warning: Removed 9108 rows containing missing values (position_stack).
```

```
## Exported graph as ~/Projects/ETH/Alessia/mobilome/graph.pptx
```



```r
physeq_drug_multi %>%
  microbiome::transform(transform = "log10p") %>% 
  speedyseq::psmelt() %>% 
  ggplot(aes(x = as.factor(Day_of_Treatment), y = Abundance, fill = Drug_Class_multi)) +   
  geom_bar( stat = "identity", colour="black", size=0.025) +
  scale_fill_manual(values=col) +
  #geom_bar(stat="identity",colour=NA,size=0) +
  facet_grid( ~ Model +  Treatment + Fermentation + Antibiotic_mg.mL , scales = "free", space = "free_x") +
  ggtitle("ARG Type Abundance per Sample") +
  xlab("Day (Treatment)")  +
  ylab("log10 TPM") + 
  theme_light() -> p_drug

p_drug
```

![](combined_resistome_visualization_files/figure-html/unnamed-chunk-26-1.png)<!-- -->

```r
p_drug %>%
  export::graph2ppt(append = TRUE, width = 18, height = 6,
                    file = out_pptx)
```

```
## Exported graph as ~/Projects/ETH/Alessia/mobilome/graph.pptx
```


```r
physeq_drug_multi %>%
  microbiome::transform(transform = "log10p") %>% 
  speedyseq::psmelt() %>% 
  ggplot(aes(x = as.factor(Day_of_Treatment), y = Abundance, fill = Drug_Class_multi)) +   
  geom_bar( stat = "identity", position="fill", colour="black", size=0.025) +
  scale_fill_manual(values=col) +
  #geom_bar(stat="identity",colour=NA,size=0) +
  facet_grid(. ~ Model + Treatment + Fermentation + Antibiotic_mg.mL , scales = "free", space = "free_x", drop = TRUE) +
  ggtitle("ARG Type Abundance per Sample") +
  xlab("Day (Treatment)")  +
  ylab("Proportion") +
  theme_light() -> p_drug_prop

p_drug_prop
```

![](combined_resistome_visualization_files/figure-html/unnamed-chunk-27-1.png)<!-- -->

```r
p_drug_prop %>%
  export::graph2ppt(append = TRUE, width = 18, height = 6,
                    file = out_pptx)
```

```
## Exported graph as ~/Projects/ETH/Alessia/mobilome/graph.pptx
```





```r
physeq_AMRgn_round <- physeq_AMRgn
round(otu_table(physeq_AMRgn_round)) -> otu_table(physeq_AMRgn_round)


physeq_AMRgn_round %>% 
  phyloseq_alphas(phylo = FALSE) -> alpha_df
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
##  Copyright (C) 2011-2020 Leo Lahti, 
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
## The following object is masked from 'package:vegan':
## 
##     diversity
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
## * id -> id...1
## * id -> id...62
## * id -> id...85
```



```r
# plot panels for each treatment
alpha_df %>% 
    pivot_longer(cols = all_of(c("Observed", "diversity_shannon")), values_to = "value", names_to = 'alphadiversiy', values_drop_na  = TRUE) %>%
    mutate(alphadiversiy = fct_relevel(alphadiversiy, c("Observed", "diversity_shannon"))) %>% 
  # filter(!is.na(Fermentation)) %>% 
  # filter(Treatment != "DONOR") %>%
  ggplot(aes(x = Day_of_Treatment, y = value, color = Treatment)) +
  geom_point() + 
  geom_line(linetype = "dashed",  size = 0.25) +
  labs(title = "AMR Gene Richness",
       y = "Richness", x = "Days of Treatment") +
  scale_color_manual(values = treat_col,
                     na.value = "black") +
  facet_grid(alphadiversiy  ~  Model + Treatment + Fermentation + Antibiotic_mg.mL , scales = "free", space = "fixed", drop = TRUE) +
  theme_light() -> p_alpha

p_alpha
```

```
## geom_path: Each group consists of only one observation. Do you need to adjust
## the group aesthetic?
## geom_path: Each group consists of only one observation. Do you need to adjust
## the group aesthetic?
## geom_path: Each group consists of only one observation. Do you need to adjust
## the group aesthetic?
## geom_path: Each group consists of only one observation. Do you need to adjust
## the group aesthetic?
```

![](combined_resistome_visualization_files/figure-html/unnamed-chunk-30-1.png)<!-- -->

```r
p_alpha %>%
  export::graph2ppt(append = TRUE, width = 20, height = 6,
                    file = out_pptx)
```

```
## geom_path: Each group consists of only one observation. Do you need to adjust
## the group aesthetic?
## geom_path: Each group consists of only one observation. Do you need to adjust
## the group aesthetic?
## geom_path: Each group consists of only one observation. Do you need to adjust
## the group aesthetic?
## geom_path: Each group consists of only one observation. Do you need to adjust
## the group aesthetic?
```

```
## Exported graph as ~/Projects/ETH/Alessia/mobilome/graph.pptx
```

Beta Diversity:


```r
# sample_data(physeq)$Reactor_Treatment <- fct_relevel(sample_data(physeq)$Reactor_Treatment, "DONOR", "CR_UNTREATED", "TR1_CTX+HV292.1",
#                                                      "TR2_CTX","TR3_HV292.1", "TR4_VAN", "TR5_VAN+CCUG59168", "TR6_CCUG59168") 
# 
# sample_data(physeq)$Treatment <- fct_relevel(sample_data(physeq)$Treatment, "DONOR", "UNTREATED",  "CTX+HV292.1", "CTX","HV292.1","VAN+CCUG59168", "VAN",  "CCUG59168") 


# physeq %>% 
#   rarefy_even_depth(sample.size = 43,
#                     rngseed = 123) -> phyloseq_rare

physeq_AMRgn %>%
  phyloseq_compute_bdiv(norm = "pc",
                        phylo = FALSE,
                        seed = 123) -> bdiv_list
```

```
## Loading required package: ape
```

```
## Loading required package: GUniFrac
```

```
## Registered S3 method overwritten by 'rmutil':
##   method         from
##   print.response httr
```

```r
physeq_AMRgn  %>%
  # subset_samples(Treatment != "DONOR") %>%
  phyloseq_plot_bdiv(dlist = bdiv_list,
                     seed = 123,
                     axis1 = 1,
                     axis2 = 2)  -> pcoa
```

```
## [1] "bray"
## [1] "sorensen"
## [1] "bjaccard"
## [1] "wjaccard"
```

```r
# phyloseq_plot_bdiv(bdiv_list,
#                    # m = "CoDa",
#                    seed = 123) -> coda
# 
pcoa$wjaccard$layers = NULL

pcoa$wjaccard + geom_point(size=3,
                           aes(color = Treatment, 
                               fill = NULL,
                               shape = Fermentation %>%  as.factor(),
                               alpha = Day_of_Treatment)) + 
  theme_light() +
  geom_path(arrow = arrow(type = "open", angle = 30, length = unit(0.15, "inches")),
            size = 0.05, linetype = "dashed", inherit.aes = TRUE, aes(group=Reactor_Treatment), show.legend = FALSE) +
  scale_alpha_continuous(range=c( 0.9, 0.3)) +
  scale_fill_manual(values = treat_col,
                     na.value = "black") +
  scale_color_manual(values = treat_col,
                     na.value = "black") +
  scale_shape_manual(values = c(15, 19), na.value =  17) + 
  theme_classic() -> p1 # +
  # labs(col=NULL, fill = NULL, shape = NULL) + guides(shape=TRUE) -> p1

p1 + facet_grid(Model ~.)
```

![](combined_resistome_visualization_files/figure-html/unnamed-chunk-31-1.png)<!-- -->

```r
p1 %>%
  export::graph2ppt(append = TRUE, width = 8, height = 6,
                    file = out_pptx)
```

```
## Exported graph as ~/Projects/ETH/Alessia/mobilome/graph.pptx
```


```r
physeq_AMRgn %>% 
  tax_table() %>% 
  data.frame() %>% 
  rownames_to_column("id") %>% 
  mutate(gene = id) %>% 
  column_to_rownames("id") %>% 
  as.matrix %>% 
  tax_table() -> tax_table(physeq_AMRgn)



physeq_AMRgn  %>%
  # subset_samples(Treatment != "DONOR") %>% 
  phyloseq_add_taxa_vector(dist = bdiv_list$wjaccard,
                           phyloseq = .,
                           figure_ord = p1,
                           tax_rank_plot = "Best_Hit_ARO", taxrank_glom = "Best_Hit_ARO",
                           top_r = 10, fact = 0.6) -> pco_env

# pco_env$plot %>% 
#   export::graph2ppt(append = TRUE, width = 8, height = 6,
#                     file = out_pptx)

pco_env$signenvfit %>% 
  DT::datatable()
```

```{=html}
<div id="htmlwidget-73e0abf6f248952da8e3" style="width:100%;height:auto;" class="datatables html-widget"></div>
<script type="application/json" data-for="htmlwidget-73e0abf6f248952da8e3">{"x":{"filter":"none","vertical":false,"data":[["1","2","3","4","5","6","7","8","9","10"],[-0.764875563622303,-0.749888192247643,-0.76645666408323,-0.764226283136953,-0.750306639367613,-0.773135872269911,-0.760002320793731,-0.76474219154319,-0.761048768322285,-0.750050085633061],[-0.620078572342329,-0.633661119273334,-0.616970541338826,-0.615824402646257,-0.633977504549405,-0.608996040449816,-0.624714554287498,-0.615627075927356,-0.621106712346566,-0.633250441018094],["emrR","emrB","evgS","evgA","Escherichia.coli.mdfA","mdtP","mdtF","Escherichia.coli.acrA","acrB","msbA"],[0.001,0.001,0.001,0.001,0.001,0.001,0.001,0.001,0.001,0.001],[0.969532063704638,0.963858714951172,0.968108466797518,0.963281506731942,0.964887529353813,0.96861525427411,0.967871801950485,0.963827316141148,0.964968775886826,0.963581252007772],["emrR","emrB","evgS","evgA",null,"mdtP","mdtF",null,"acrB","msbA"],[null,null,null,null,null,null,null,null,null,null],["emrR","emrB","evgS","evgA",null,"mdtP","mdtF",null,"acrB","msbA"]],"container":"<table class=\"display\">\n  <thead>\n    <tr>\n      <th> <\/th>\n      <th>Axis.1<\/th>\n      <th>Axis.2<\/th>\n      <th>id<\/th>\n      <th>pval<\/th>\n      <th>r<\/th>\n      <th>tax_rank_plot<\/th>\n      <th>Drug_Class<\/th>\n      <th>gene<\/th>\n    <\/tr>\n  <\/thead>\n<\/table>","options":{"columnDefs":[{"className":"dt-right","targets":[1,2,4,5]},{"orderable":false,"targets":0}],"order":[],"autoWidth":false,"orderClasses":false}},"evals":[],"jsHooks":[]}</script>
```


```r
physeq_drug_multi  %>%
  # subset_samples(Treatment != "DONOR") %>% 
  phyloseq_add_taxa_vector(dist = bdiv_list$wjaccard,
                           phyloseq = .,
                           figure_ord = p1,
                           tax_rank_plot = "Drug_Class_multi", taxrank_glom = "Drug_Class_multi",
                           top_r = 10, fact = 0.6) -> pco_env


pco_env$plot + facet_grid(Model ~ .) -> pco_env$plot

pco_env$plot
```

```
## Warning: Removed 2 rows containing missing values (geom_text_repel).
```

![](combined_resistome_visualization_files/figure-html/unnamed-chunk-33-1.png)<!-- -->

```r
pco_env$plot  %>%
  export::graph2ppt(append = TRUE, width = 8, height = 6,
                    file = out_pptx)
```

```
## Warning: Removed 2 rows containing missing values (geom_text_repel).
```

```
## Exported graph as ~/Projects/ETH/Alessia/mobilome/graph.pptx
```

```r
pco_env$signenvfit %>% 
  DT::datatable()
```

```{=html}
<div id="htmlwidget-1516fd1423c88eb68dbe" style="width:100%;height:auto;" class="datatables html-widget"></div>
<script type="application/json" data-for="htmlwidget-1516fd1423c88eb68dbe">{"x":{"filter":"none","vertical":false,"data":[["1","2","3","4","5","6","7","8","9","10"],[-0.682540917380305,0.516701010135203,-0.568842203237661,0.665502841815607,-0.770952454498952,0.52759244039967,-0.684742557423979,-0.169784633761131,0.704272179224744,-0.750050085633061],[-0.577807956861462,-0.53980340566923,-0.680701800277176,-0.599367457023382,-0.619240649902978,0.519296562294092,-0.644431871785898,0.697201389616232,-0.612133109214062,-0.633250441018094],["fluoroquinolone","cephalosporin","multi.drug","macrolide","peptide","tetracycline","aminocoumarin","lincosamide","cephamycin","nitroimidazole"],[0.001,0.001,0.001,0.001,0.001,0.001,0.001,0.001,0.001,0.001],[0.799724138910765,0.558367650646838,0.786936393084865,0.802135381003324,0.977826669590222,0.548022702777341,0.884164807321006,0.514916599544206,0.870706245826045,0.963581252007772],["fluoroquinolone","cephalosporin",null,"macrolide","peptide","tetracycline","aminocoumarin","lincosamide","cephamycin","nitroimidazole"],[null,null,null,null,null,null,null,null,null,null]],"container":"<table class=\"display\">\n  <thead>\n    <tr>\n      <th> <\/th>\n      <th>Axis.1<\/th>\n      <th>Axis.2<\/th>\n      <th>id<\/th>\n      <th>pval<\/th>\n      <th>r<\/th>\n      <th>tax_rank_plot<\/th>\n      <th>Best_Hit_ARO<\/th>\n    <\/tr>\n  <\/thead>\n<\/table>","options":{"columnDefs":[{"className":"dt-right","targets":[1,2,4,5]},{"orderable":false,"targets":0}],"order":[],"autoWidth":false,"orderClasses":false}},"evals":[],"jsHooks":[]}</script>
```

Make the heatmap pretty 

```r
physeq %>% 
  # speedyseq::tax_glom("Best_Hit_ARO") %>% 
  speedyseq::psmelt() %>% 
  rename(TPM = Abundance) %>% 
  mutate(TPM = na_if(TPM, 0)) %>%
  # mutate(logTPM = log10(TPM + 1)) %>% 
  # filter(Drug_Class %in% c("cephalosporin", "glycopeptide")) %>% 
  # mutate(Drug_Class = fct_relevel(Drug_Class, c("cephalosporin", "glycopeptide"))) %>% 
  filter(Best_Hit_ARO %in% c("vanZA", "vanYA", "vanXA", "vanA", "vanHA", "vanRA", "CTX-M-1")) %>% 
  # mutate(Classification = fct_relevel(Classification, c("MAG", "Chromosome", "Plasmid", "Uncertain - plasmid or chromosomal"))) %>% 
  ggplot(mapping = aes(x = Day_of_Treatment ,
                       y = Best_Hit_ARO ,
                       fill = TPM,
  )) +
  facet_grid(Drug_Class ~  Model + Treatment + Fermentation + Antibiotic_mg.mL, scales = "free",  space = "fixed", drop = TRUE) +
  # facet_grid(. ~ Treatment,  scales = "free", space = "free", drop = TRUE) +
  # facet_grid(. ~ Reactor_Treatment, scales = "free", space = "free", drop = TRUE,labeller = labeller(Drug_Class = label_wrap_gen(width = 15))) +
  #we adjust the vertical labels such that it is printed correctly. 
  geom_tile() +
  theme_bw()  +
  # scale_fill_viridis_c(name = "TPM",
  #                      na.value = "transparent", trans = scales::pseudo_log_trans(sigma = 0.1),
  #                      # trans =  scales::pseudo_log_trans(sigma = 1),
  #                      breaks = c(1, 10, 50, 100, 400), labels = c(1, 10, 50, 100, 400),
  #                      limits = c(0, 500)) + #need to make sure that we include max value of tpm as upper limit
  scale_fill_viridis_c(na.value = "transparent") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        strip.text.x = element_text(size=8),
        #change the size of the label and the angle so the text is printed horizontally
        strip.text.y = element_text(size = 7,angle=0)) -> perma_bbuble

perma_bbuble +
  theme(legend.position = "bottom") -> perma_bbuble

perma_bbuble
```

![](combined_resistome_visualization_files/figure-html/unnamed-chunk-34-1.png)<!-- -->

```r
perma_bbuble %>%
  export::graph2ppt(append = TRUE, width = 16, height = 4,
                    file = out_pptx)
```

```
## Exported graph as ~/Projects/ETH/Alessia/mobilome/graph.pptx
```

What's next? 
- MAG, plasmid, Chromosome ... proportions.
- SQM contig taxonomy coloring/facetting.



```r
source("https://raw.githubusercontent.com/fconstancias/DivComAnalyses/master/R/phyloseq_normalisation.R")


ps_tmp_tax <- physeq


ps_tmp_tax %>% 
  tax_table() %>% 
  data.frame() %>% 
  # select(Tax) %>% 
  separate(Tax, into = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"), sep = ";") %>% 
  mutate(across(c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"), gsub, pattern = "[a-z]_", replacement = "")) %>% 
  as.matrix() -> tax_table(ps_tmp_tax)
```

```
## Warning: Expected 7 pieces. Additional pieces discarded in 23 rows [8, 10,
## 11, 13, 41, 44, 58, 79, 83, 95, 110, 135, 138, 139, 166, 181, 189, 207, 210,
## 220, ...].
```

```
## Warning: Expected 7 pieces. Missing pieces filled with `NA` in 111 rows [4, 6,
## 9, 16, 17, 21, 23, 26, 29, 30, 33, 36, 45, 48, 51, 52, 53, 55, 56, 57, ...].
```

```r
ps_tmp_tax_2 <- ps_tmp_tax

ps_tmp_tax %>% 
   tax_table() %>% 
  data.frame() %>% 
  select(c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"))  %>% 
    as.matrix() -> tax_table(ps_tmp_tax_2)

ps_tmp_tax_2 %>%   
  phyloseq_get_strains_fast(species = FALSE) %>% 
  tax_table() %>% 
  data.frame()  %>% 
    mutate(Strain =  gsub(" human_megahit.*","",Strain)) %>% 
    mutate(Strain =  gsub(" chick_megahit.*","",Strain)) -> full_tax
```

```
## Joining, by = "ASV"
```

```r
ps_tmp_tax %>% 
  tax_table() %>% 
  data.frame() %>% 
  rownames_to_column("ASV") %>% 
  # select(-Tax) %>% 
  left_join(full_tax %>% 
            select(Strain) %>% rownames_to_column("ASV"),
             by = c("ASV" = "ASV")) %>% 
  select(-Predicted_DNA) %>% 
  column_to_rownames("ASV") %>% 
  as.matrix() -> tax_table(ps_tmp_tax)

  
  
  
tax_table(ps_tmp_tax) %>% 
    data.frame() %>% 
  select(c("Best_Hit_ARO", "Strain")) %>% 
  head()
```

```
##                                Best_Hit_ARO                      Strain
## human_megahit_1512_86-616              emrR            Escherichia coli
## human_megahit_1512_743-1915            emrA            Escherichia coli
## human_megahit_1512_1932-3470           emrB            Escherichia coli
## human_megahit_1935_13232-13930        vanRD  unknown Firmicutes (Class)
## human_megahit_4077_533-1408         CTX-M-1            Escherichia coli
## human_megahit_4631_5392-6930           emrB unknown Citrobacter (Genus)
```

```r
ps_tmp_tax %>% 
  saveRDS("~/Projects/ETH/Alessia/mobilome/ps_combined_tax_clean.rds")


as(tax_table(ps_tmp_tax), "matrix") %>% 
  as.data.frame() %>%
  rownames_to_column('ASV') -> tax_table_ps_combined

as(otu_table(ps_tmp_tax), "matrix") %>% 
  as.data.frame() %>%
  rownames_to_column('ASV') -> otu_table_ps_combined

ps_tmp_tax %>% 
  sample_data() %>% 
  data.frame() %>% 
  rownames_to_column('sample_id') ->  meta_ps_combined

otu_table_ps_combined %>%
  # left_join(refseq_df, by = 'ASV') %>%
  left_join(tax_table_ps_combined, by = c("ASV" = "ASV")) %>%
  dplyr::select(ASV, everything()) -> merged_ps_combined


merged_ps_combined %>% 
  write_tsv("~/Projects/ETH/Alessia/mobilome/ps_combined_tax_clean_rds.tsv")
```


```r
tax_table(ps_tmp_tax) %>% 
    data.frame() %>% 
  select(c("Best_Hit_ARO", "Strain")) %>% 
  group_by(Strain) %>% 
  add_count() %>% 
  arrange(-n)
```

```
## # A tibble: 244 × 3
## # Groups:   Strain [29]
##    Best_Hit_ARO          Strain               n
##    <chr>                 <chr>            <int>
##  1 emrR                  Escherichia coli    99
##  2 emrA                  Escherichia coli    99
##  3 emrB                  Escherichia coli    99
##  4 CTX-M-1               Escherichia coli    99
##  5 baeS                  Escherichia coli    99
##  6 CRP                   Escherichia coli    99
##  7 YojI                  Escherichia coli    99
##  8 evgS                  Escherichia coli    99
##  9 evgA                  Escherichia coli    99
## 10 Escherichia coli mdfA Escherichia coli    99
## # … with 234 more rows
```

Define colors for Taxonomic information:


```r
ps_tmp_tax %>%  speedyseq::psmelt() %>%  pull(Strain) %>%  unique() -> strains
set.seed(456)
col_strains <- distinctColorPalette(length(strains),
                               runTsne = FALSE)

names(col_strains) <- strains
```


```r
# see Sneha's & Hannah's code for "plot.df.sep <- separate_rows(plot.df, Drug.Class,sep = '; ',convert = TRUE)"

ps_tmp_tax %>%
  microbiome::transform(transform = "log10p") %>% 
  speedyseq::psmelt() %>% 
    # left_join(tax_mapping,
    #         by = c("Best_Hit_ARO" = "Best_Hit_ARO"),
    #         suffix = c("_x", "")) %>%
    mutate(Abundance = na_if(Abundance, 0)) %>%
  mutate(logTPM = log10(Abundance + 1))  %>% 
  # separate_rows(., Resistance_Mechanism,sep = '; ',convert = TRUE) %>% 
  ggplot(aes(x = as.factor(Day_of_Treatment), y = Abundance, fill = Strain)) +  
  geom_bar( stat = "identity", colour="black", size=0.025) +
  facet_grid(Drug_Class_multi  ~ Model +  Treatment + Fermentation + Antibiotic_mg.mL , scales = "free", space = "free_x") +
  ggtitle("ARG Type Abundance per Sample") +
  xlab("Day (Treatment)")  +
  ylab("log10 TPM") + 
  theme_light() +
  scale_fill_manual(name = "", values = col_strains) +
  theme(legend.position = "bottom") -> p_mec

# p_mec %>% 
#   ggpubr::set_palette("jco") -> p_mec

p_mec
```

```
## Warning: Removed 15421 rows containing missing values (position_stack).
```

![](combined_resistome_visualization_files/figure-html/unnamed-chunk-38-1.png)<!-- -->

```r
p_mec %>%
  export::graph2ppt(append = TRUE, width = 20, height = 18,
                    file = out_pptx)
```

```
## Warning: Removed 15421 rows containing missing values (position_stack).
```

```
## Exported graph as ~/Projects/ETH/Alessia/mobilome/graph.pptx
```


```r
ps_tmp_tax %>%
  microbiome::transform(transform = "log10p") %>% 
  speedyseq::psmelt() %>% 
    # left_join(tax_mapping,
    #         by = c("Best_Hit_ARO" = "Best_Hit_ARO"),
    #         suffix = c("_x", "")) %>%
    mutate(Abundance = na_if(Abundance, 0)) %>%
  mutate(logTPM = log10(Abundance + 1))  %>% 
  # separate_rows(., Resistance_Mechanism,sep = '; ',convert = TRUE) %>% 
  ggplot(aes(x = as.factor(Day_of_Treatment), y = Abundance, fill = Drug_Class_multi)) +  
  geom_bar( stat = "identity", colour="black", size=0.025) +
  facet_grid( Strain   ~ Model +  Treatment + Fermentation + Antibiotic_mg.mL , scales = "free", space = "free_x") +
  ggtitle("ARG Type Abundance per Sample") +
  xlab("Day (Treatment)")  +
  ylab("log10 TPM") + 
  theme_light() +
  scale_fill_manual(name = "", values = col) +
  theme(legend.position = "bottom") -> p_mec

# p_mec %>% 
#   ggpubr::set_palette("jco") -> p_mec

p_mec
```

```
## Warning: Removed 15421 rows containing missing values (position_stack).
```

![](combined_resistome_visualization_files/figure-html/unnamed-chunk-39-1.png)<!-- -->

```r
p_mec %>%
  export::graph2ppt(append = TRUE, width = 20, height = 18,
                    file = out_pptx)
```

```
## Warning: Removed 15421 rows containing missing values (position_stack).
```

```
## Exported graph as ~/Projects/ETH/Alessia/mobilome/graph.pptx
```


```r
# see Sneha's & Hannah's code for "plot.df.sep <- separate_rows(plot.df, Drug.Class,sep = '; ',convert = TRUE)"

ps_tmp_tax %>%
  microbiome::transform(transform = "log10p") %>% 
  speedyseq::psmelt() %>% 
    # left_join(tax_mapping,
    #         by = c("Best_Hit_ARO" = "Best_Hit_ARO"),
    #         suffix = c("_x", "")) %>%
    mutate(Abundance = na_if(Abundance, 0)) %>%
  mutate(logTPM = log10(Abundance + 1))  %>% 
  # separate_rows(., Resistance_Mechanism,sep = '; ',convert = TRUE) %>% 
  ggplot(aes(x = as.factor(Day_of_Treatment), y = Abundance, fill = Strain)) +  
  geom_bar( stat = "identity", colour="black", position = "fill", size=0.05) +
  facet_grid(Drug_Class_multi  ~ Model +  Treatment + Fermentation + Antibiotic_mg.mL , scales = "free", space = "free_x") +
  ggtitle("ARG Type Abundance per Sample") +
  xlab("Day (Treatment)")  +
  ylab("Proportion") + 
  theme_light() +
  scale_fill_manual(name = "", values = col_strains) +
  theme(legend.position = "bottom") -> p_mec

# p_mec %>% 
#   ggpubr::set_palette("jco") -> p_mec

p_mec
```

```
## Warning: Removed 15421 rows containing missing values (position_stack).
```

![](combined_resistome_visualization_files/figure-html/unnamed-chunk-40-1.png)<!-- -->

```r
p_mec %>%
  export::graph2ppt(append = TRUE, width = 20, height = 18,
                    file = out_pptx)
```

```
## Warning: Removed 15421 rows containing missing values (position_stack).
```

```
## Exported graph as ~/Projects/ETH/Alessia/mobilome/graph.pptx
```


```r
# see Sneha's & Hannah's code for "plot.df.sep <- separate_rows(plot.df, Drug.Class,sep = '; ',convert = TRUE)"

ps_tmp_tax %>%
  microbiome::transform(transform = "log10p") %>% 
  speedyseq::psmelt() %>% 
    # left_join(tax_mapping,
    #         by = c("Best_Hit_ARO" = "Best_Hit_ARO"),
    #         suffix = c("_x", "")) %>%
    mutate(Abundance = na_if(Abundance, 0)) %>%
  mutate(logTPM = log10(Abundance + 1))  %>% 
  # separate_rows(., Resistance_Mechanism,sep = '; ',convert = TRUE) %>% 
  ggplot(aes(x = as.factor(Day_of_Treatment), y = Abundance, fill = Strain)) +  
  geom_bar( stat = "identity", colour="black", size=0.025) +
  facet_grid(.  ~ Model +  Treatment + Fermentation + Antibiotic_mg.mL , scales = "free", space = "free_x") +
  ggtitle("ARG Type Abundance per Sample") +
  xlab("Day (Treatment)")  +
  ylab("log10 TPM") + 
  theme_light() +
  scale_fill_manual(name = "", values = col_strains) +
  theme(legend.position = "bottom") -> p_mec

# p_mec %>% 
#   ggpubr::set_palette("jco") -> p_mec

p_mec
```

```
## Warning: Removed 15421 rows containing missing values (position_stack).
```

![](combined_resistome_visualization_files/figure-html/unnamed-chunk-41-1.png)<!-- -->

```r
p_mec %>%
  export::graph2ppt(append = TRUE, width = 18, height = 6,
                    file = out_pptx)
```

```
## Warning: Removed 15421 rows containing missing values (position_stack).
```

```
## Exported graph as ~/Projects/ETH/Alessia/mobilome/graph.pptx
```


```r
# see Sneha's & Hannah's code for "plot.df.sep <- separate_rows(plot.df, Drug.Class,sep = '; ',convert = TRUE)"

ps_tmp_tax %>%
  microbiome::transform(transform = "log10p") %>% 
  speedyseq::psmelt() %>% 
    # left_join(tax_mapping,
    #         by = c("Best_Hit_ARO" = "Best_Hit_ARO"),
    #         suffix = c("_x", "")) %>%
    mutate(Abundance = na_if(Abundance, 0)) %>%
  mutate(logTPM = log10(Abundance + 1))  %>% 
  # separate_rows(., Resistance_Mechanism,sep = '; ',convert = TRUE) %>% 
  # ggplot(aes(x = as.factor(Day_of_Treatment), y = Abundance, fill = fct_reorder(Strain, logTPM))) +   
  ggplot(aes(x = as.factor(Day_of_Treatment), y = Abundance, fill = reorder(Strain, logTPM, FUN=median, na.rm = TRUE))) +
  geom_bar( stat = "identity", position="fill", colour="black", size=0.025) +
  facet_grid(.  ~ Model +  Treatment + Fermentation + Antibiotic_mg.mL , scales = "free", space = "free_x") +
  ggtitle("ARG Type Abundance per Sample") +
  xlab("Day (Treatment)")  +
  ylab("Proportion") + 
  theme_light() +
  scale_fill_manual(name = "", values = col_strains) +
  theme(legend.position = "bottom") -> p_mec

# p_mec %>% 
#   ggpubr::set_palette("jco") -> p_mec

p_mec
```

```
## Warning: Removed 15421 rows containing missing values (position_stack).
```

![](combined_resistome_visualization_files/figure-html/unnamed-chunk-42-1.png)<!-- -->

```r
p_mec %>%
  export::graph2ppt(append = TRUE, width = 18, height = 6,
                    file = out_pptx)
```

```
## Warning: Removed 15421 rows containing missing values (position_stack).
```

```
## Exported graph as ~/Projects/ETH/Alessia/mobilome/graph.pptx
```


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
##  [1] GUniFrac_1.4         ape_5.6              microbiome_1.14.0   
##  [4] metagMisc_0.0.4      gdtools_0.2.3        reshape2_1.4.4      
##  [7] scales_1.1.1         randomcoloR_1.1.0.1  vegan_2.5-7         
## [10] lattice_0.20-45      permute_0.9-5        RColorBrewer_1.1-2  
## [13] microViz_0.9.0       here_1.0.1           ggrepel_0.9.1       
## [16] speedyseq_0.5.3.9018 phyloseq_1.36.0      forcats_0.5.1       
## [19] stringr_1.4.0        dplyr_1.0.7          purrr_0.3.4         
## [22] readr_2.1.0          tidyr_1.1.4          tibble_3.1.6        
## [25] ggplot2_3.3.5        tidyverse_1.3.1.9000
## 
## loaded via a namespace (and not attached):
##   [1] utf8_1.2.2             tidyselect_1.1.1       lme4_1.1-27.1         
##   [4] htmlwidgets_1.5.4      grid_4.1.2             Rtsne_0.15            
##   [7] modeest_2.4.0          munsell_0.5.0          codetools_0.2-18      
##  [10] statmod_1.4.36         DT_0.20                withr_2.4.3           
##  [13] colorspace_2.0-2       Biobase_2.52.0         highr_0.9             
##  [16] knitr_1.36             uuid_1.0-3             rstudioapi_0.13       
##  [19] stats4_4.1.2           ggsignif_0.6.3         officer_0.4.0         
##  [22] labeling_0.4.2         GenomeInfoDbData_1.2.6 bit64_4.0.5           
##  [25] farver_2.1.0           pheatmap_1.0.12        rhdf5_2.36.0          
##  [28] fBasics_3042.89.1      rprojroot_2.0.2        vctrs_0.3.8           
##  [31] generics_0.1.1         xfun_0.28              R6_2.5.1              
##  [34] GenomeInfoDb_1.28.4    clue_0.3-60            bitops_1.0-7          
##  [37] rhdf5filters_1.4.0     assertthat_0.2.1       vroom_1.5.6           
##  [40] gtable_0.3.0           spatial_7.3-14         timeDate_3043.102     
##  [43] rlang_0.4.12           systemfonts_1.0.2      splines_4.1.2         
##  [46] rstatix_0.7.0          broom_0.7.11           rgl_0.107.14          
##  [49] yaml_2.2.1             abind_1.4-5            modelr_0.1.8          
##  [52] crosstalk_1.2.0        backports_1.4.1        tools_4.1.2           
##  [55] ellipsis_0.3.2         jquerylib_0.1.4        biomformat_1.20.0     
##  [58] BiocGenerics_0.38.0    stargazer_5.2.2        stabledist_0.7-1      
##  [61] devEMF_4.0-2           Rcpp_1.0.7             plyr_1.8.6            
##  [64] base64enc_0.1-3        zlibbioc_1.38.0        RCurl_1.98-1.5        
##  [67] ggpubr_0.4.0           rpart_4.1-15           viridis_0.6.2         
##  [70] statip_0.2.3           S4Vectors_0.30.2       haven_2.4.3           
##  [73] cluster_2.1.2          fs_1.5.2               magrittr_2.0.1        
##  [76] data.table_1.14.2      timeSeries_3062.100    openxlsx_4.2.4        
##  [79] lmerTest_3.1-3         flextable_0.6.9        reprex_2.0.1          
##  [82] matrixStats_0.61.0     hms_1.1.1              evaluate_0.14         
##  [85] xtable_1.8-4           rio_0.5.27             readxl_1.3.1          
##  [88] IRanges_2.26.0         gridExtra_2.3          compiler_4.1.2        
##  [91] V8_3.6.0               crayon_1.4.2           minqa_1.2.4           
##  [94] htmltools_0.5.2        mgcv_1.8-38            tzdb_0.2.0            
##  [97] lubridate_1.8.0        DBI_1.1.1              rmutil_1.1.5          
## [100] dbplyr_2.1.1           MASS_7.3-54            boot_1.3-28           
## [103] Matrix_1.3-4           ade4_1.7-18            car_3.0-11            
## [106] cli_3.1.0              parallel_4.1.2         igraph_1.2.10         
## [109] pkgconfig_2.0.3        numDeriv_2016.8-1.1    foreign_0.8-81        
## [112] xml2_1.3.2             foreach_1.5.1          bslib_0.3.1           
## [115] multtest_2.48.0        XVector_0.32.0         rvg_0.2.5             
## [118] rvest_1.0.2            digest_0.6.29          Biostrings_2.60.2     
## [121] rmarkdown_2.11         cellranger_1.1.0       curl_4.3.2            
## [124] export_0.3.0           nloptr_1.2.2.2         lifecycle_1.0.1       
## [127] nlme_3.1-153           jsonlite_1.7.2         Rhdf5lib_1.14.2       
## [130] carData_3.0-4          viridisLite_0.4.0      fansi_0.5.0           
## [133] pillar_1.6.4           ggsci_2.9              fastmap_1.1.0         
## [136] httr_1.4.2             survival_3.2-13        glue_1.6.0            
## [139] zip_2.2.0              iterators_1.0.13       bit_4.0.4             
## [142] stringi_1.7.6          sass_0.4.0             stable_1.1.4
```

