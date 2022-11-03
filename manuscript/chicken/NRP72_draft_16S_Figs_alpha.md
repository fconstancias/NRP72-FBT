---
title: "NTP72 - 16S - alpha - chicken draft "
author: "Florentin Constancias"
date: "October 27, 2022"
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
source("https://raw.githubusercontent.com/fconstancias/DivComAnalyses/master/R/phyloseq_alpha.R") 
source("https://raw.githubusercontent.com/fconstancias/DivComAnalyses/master/R/phyloseq_beta.R") 
source("https://raw.githubusercontent.com/fconstancias/DivComAnalyses/master/R/phyloseq_heatmap.R")
source("https://raw.githubusercontent.com/fconstancias/metabaRpipe-source/master/Rscripts/functions.R")
source("https://raw.githubusercontent.com/fconstancias/DivComAnalyses/master/R/phyloseq_normalisation.R") 
source("https://raw.githubusercontent.com/fconstancias/DivComAnalyses/master/R/phyloseq_varia.R")


'%!in%' <- function(x,y)!('%in%'(x,y))
```


```r
out_pptx = "~/Desktop/16S-chicken-alpha.pptx"
out_xlsx = "~/Desktop/16S-chicken-alpha.xlsx"
```


```r
load(url("https://github.com/fconstancias/NRP72-FBT/blob/master/Figures/Rastetics.Rdata?raw=true"))
```

## load data:


```r
"data/processed/16S/16S_working_phyloseq.RDS" %>% 
  here::here() %>% 
  readRDS()  %>% 
  subset_samples(Enrichment == "NotEnriched") %>%
  subset_samples(Experiment %in% c("Continuous", "Cecum")) %>%
  subset_samples(Day_of_Treatment > -6 & Day_of_Treatment <= 30 | Reactor == "DONOR") %>% 
  subset_samples(Model == "Chicken") %>% 
  subset_samples(Reactor != "CR" & Model2 == "Chicken1" | Model2 == "Chicken2" & Reactor %in% c("TR2", "TR3", "TR4", "TR5", "TR6", "TR7", "CR", "DONOR")) %>% 
  subset_samples(Reactor != "IR")  -> ps_filtered
```

## rarefraction


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
## 639OTUs were removed because they are no longer 
## present in any sample after random subsampling
```

```
## ...
```

```r
ps_filtered %>%
  prune_samples(sample_sums(.)>= min_sample_size, .) -> ps_fil
```

Alpha Diversity:


```r
ps_rare %>% 
  subset_samples(Treatment != "DONOR") %>% 
  plot_richness(measures = c("Observed", "Shannon", "InvSimpson"),
                shape = "Model2") -> alpha_plot

ps_rare %>% 
  subset_samples(Treatment != "DONOR") %>% 
  microbiome::alpha(index = "all") %>% 
  rownames_to_column("sample_id") %>% 
  pivot_longer(cols = -sample_id, names_to = "variable") %>% 
  mutate(variable = as.factor(variable)) %>% 
  left_join(ps_rare %>%
              sample_data() %>% 
              data.frame() %>%  
              rownames_to_column("sample_id"),
            by = c("sample_id" = "sample_id")
  ) -> alpha
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
# ps_rare %>% 
#   add_phylogeny_to_phyloseq() %>% 
#   metagMisc::phyloseq_phylo_div() -> phylo_alpha

alpha %>% 
  mutate(variable = recode(variable, 
                           observed="Hill-d0",
                           diversity_shannon="Hill-d1",
                           diversity_inverse_simpson = "Hill-d2")) %>% 
  mutate(value = ifelse(variable == "Hill-d1", exp(value), value)) -> hill_df
```


```r
hill_df %>%
  filter(variable %in% c("Hill-d0", "Hill-d1", "Hill-d2", "evenness_pielou")) %>% 
  mutate(variable = fct_relevel(variable, c("Hill-d0", "Hill-d1", "Hill-d2", "evenness_pielou"))) %>% 
  mutate(Day_of_Treatment = as.double(Day_of_Treatment)) %>% 
  filter(!is.na(Day_of_Treatment)) %>% 
  # filter(Day_of_Treatment >= -1) %>% 
  ##filter(samples != "B3","B4") %>%
  ggplot(aes(Day_of_Treatment,
             value,
             colour = Treatment,
             fill = Treatment)) + #,
  # group = interaction(Model2, Reactor, Treatment, Treatment_Dose, Fermentation))) +
  facet_grid(variable ~ Model2,scale="free",space="fixed") +
  # geom_boxplot(outlier.colour = NA, alpha = 0.2) +
  # ggtitle("Alpha_Div_chicken_2 ") +
  geom_point(size=2,position=position_jitterdodge(dodge.width = 0.2), 
             aes(shape = Antibiotic_mg.L)) +
  ylab("Alpha-diversity index") + xlab("Days (Treatment)") + theme_bw() +
  geom_smooth(alpha = 0.08) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  scale_color_manual(name = "", values = treat_col[! treat_col %in% c("#3CB371")],
                     na.value = "black") +
  scale_fill_manual(name = "", values = treat_col[! treat_col %in% c("#3CB371")],
                    na.value = "black") +
  scale_shape_manual(name = "" ,values = c(15,16,18,19), na.value =  17) -> alpha_plot

alpha_plot
```

![](NRP72_draft_16S_Figs_alpha_files/figure-html/unnamed-chunk-6-1.png)<!-- -->



```r
# alpha_plot + facet_null() + facet_grid(Model2 ~ variable,scale="free",space="free")
```


```r
alpha_plot + theme(legend.position = "none") -> alpha_plot_noleg

alpha_plot  %>%  ggpubr::get_legend() %>% ggpubr::as_ggplot() -> alpha_plot_leg

alpha_plot_noleg
```

![](NRP72_draft_16S_Figs_alpha_files/figure-html/unnamed-chunk-8-1.png)<!-- -->

```r
alpha_plot_noleg %>%
  export::graph2ppt(append = TRUE,
                    width = 317.48031496 * 1,
                    height = 0.618 * 317.48031496 * 1.25 , paper = "A4",  scaling = 2,
                    file = out_pptx)


alpha_plot_leg %>%
  export::graph2ppt(append = TRUE,
                    width = 317.48031496 * 1,
                    height = 0.618 * 317.48031496 * 1.25 , paper = "A4",  scaling = 2,
                    file = out_pptx)
```

```r
alpha_plot_noleg + 
  facet_null() + 
  facet_wrap(Model2  ~ variable, scale="free_y", nrow = 2) -> alpha_plot_noleg_2

alpha_plot_noleg_2
```

![](NRP72_draft_16S_Figs_alpha_files/figure-html/unnamed-chunk-9-1.png)<!-- -->

```r
alpha_plot_noleg_2 %>% 
  export::graph2ppt(append = TRUE,
                    width = 317.48031496 * 1,
                    height = 0.618 * 317.48031496 * 1.25 , paper = "A4",  scaling = 2,
                    file = out_pptx)
```


```r
alpha_plot_noleg$data %>% 
  dplyr::select(sample_name, value, variable, Reactor_Treatment, Day_of_Treatment) %>% 
  xlsx::write.xlsx(file = out_xlsx,
                   append = TRUE,
                   # sheetName = "p_gram_h2",
                   showNA = TRUE)
```

# statistical evaluation:

some periods: 
- pret - Day_of_Treatment <= 0
- treatment - Day_of_Treatment >= 1
- 0-5 5-end

pret:


```r
alpha_plot_noleg$data %>% 
  filter(Day_of_Treatment <= 0) %>% 
  # group_by(Model2) %>% 
  # group_by(Model2, variable) %>%
  # rstatix::wilcox_test(log10_ratio ~ Treatment_A,
  #                      # comparisons = list(c("Inulin_3", "Control")),
  #                      # ref.group = "Control",
  #                      data = . ) %>%  #,
  ggpubr::compare_means(value ~ Reactor_Treatment,
                        data = .,
                        method = "wilcox.test",
                        group.by = c("Model2", "variable"),
                        p.adjust.method = "none") %>% 
  arrange(p) %>% 
  ungroup() %>% 
  group_by(Model2) %>% 
  rstatix::adjust_pvalue(method = "fdr",p.col = "p", output.col = "p.adj") %>% 
  arrange(p.adj) -> p_tmp

p_tmp %>% 
  xlsx::write.xlsx(file = out_xlsx,
                   append = TRUE,
                   sheetName = "alpha_div_comp_pret",
                   showNA = TRUE)

p_tmp %>% 
  filter(p.adj <= 0.05) %>% 
  DT::datatable()
```

```{=html}
<div id="htmlwidget-0f3c9d7bb8389a54a0ab" style="width:100%;height:auto;" class="datatables html-widget"></div>
<script type="application/json" data-for="htmlwidget-0f3c9d7bb8389a54a0ab">{"x":{"filter":"none","vertical":false,"data":[[],[],[],[],[],[],[],[],[],[],[]],"container":"<table class=\"display\">\n  <thead>\n    <tr>\n      <th> <\/th>\n      <th>Model2<\/th>\n      <th>variable<\/th>\n      <th>.y.<\/th>\n      <th>group1<\/th>\n      <th>group2<\/th>\n      <th>p<\/th>\n      <th>p.adj<\/th>\n      <th>p.format<\/th>\n      <th>p.signif<\/th>\n      <th>method<\/th>\n    <\/tr>\n  <\/thead>\n<\/table>","options":{"columnDefs":[{"className":"dt-right","targets":[6,7]},{"orderable":false,"targets":0}],"order":[],"autoWidth":false,"orderClasses":false}},"evals":[],"jsHooks":[]}</script>
```
pretreatment: nothing significantly different.



```r
alpha_plot_noleg$data %>% 
  filter(Day_of_Treatment >= 1) %>% 
  # group_by(Model2) %>% 
  # group_by(Model2, variable) %>%
  # rstatix::wilcox_test(log10_ratio ~ Treatment_A,
  #                      # comparisons = list(c("Inulin_3", "Control")),
  #                      # ref.group = "Control",
  #                      data = . ) %>%  #,
  ggpubr::compare_means(value ~ Reactor_Treatment,
                        data = .,
                        method = "wilcox.test",
                        group.by = c("Model2", "variable"),
                        p.adjust.method = "none") %>% 
  arrange(p) %>% 
  ungroup() %>% 
  group_by(Model2) %>% 
  rstatix::adjust_pvalue(method = "fdr",p.col = "p", output.col = "p.adj") %>% 
  arrange(p.adj) -> p_tmp

p_tmp %>% 
  xlsx::write.xlsx(file = out_xlsx,
                   append = TRUE,
                   sheetName = "alpha_div_comp_treat",
                   showNA = TRUE)


p_tmp %>% 
  filter(p.adj <= 0.05) %>% 
  DT::datatable()
```

```{=html}
<div id="htmlwidget-22b5dc01dcc9bcca5343" style="width:100%;height:auto;" class="datatables html-widget"></div>
<script type="application/json" data-for="htmlwidget-22b5dc01dcc9bcca5343">{"x":{"filter":"none","vertical":false,"data":[["1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","23","24","25","26","27","28","29","30","31","32","33","34","35","36","37","38","39","40","41","42","43","44","45","46","47","48","49","50","51","52","53","54","55","56","57","58","59","60","61","62","63","64","65","66","67","68","69","70","71","72","73","74","75","76","77","78","79","80","81","82","83","84","85","86","87","88","89","90","91","92","93","94","95","96"],["Chicken2","Chicken2","Chicken2","Chicken2","Chicken2","Chicken2","Chicken2","Chicken2","Chicken2","Chicken2","Chicken2","Chicken2","Chicken2","Chicken2","Chicken2","Chicken2","Chicken2","Chicken2","Chicken2","Chicken2","Chicken2","Chicken2","Chicken2","Chicken2","Chicken2","Chicken2","Chicken2","Chicken2","Chicken2","Chicken2","Chicken2","Chicken2","Chicken2","Chicken2","Chicken2","Chicken2","Chicken2","Chicken2","Chicken1","Chicken1","Chicken1","Chicken1","Chicken1","Chicken1","Chicken1","Chicken2","Chicken2","Chicken2","Chicken2","Chicken2","Chicken2","Chicken1","Chicken1","Chicken2","Chicken2","Chicken2","Chicken2","Chicken1","Chicken2","Chicken2","Chicken1","Chicken1","Chicken1","Chicken1","Chicken1","Chicken1","Chicken1","Chicken2","Chicken1","Chicken2","Chicken1","Chicken1","Chicken2","Chicken1","Chicken1","Chicken1","Chicken1","Chicken2","Chicken1","Chicken2","Chicken1","Chicken2","Chicken1","Chicken2","Chicken2","Chicken1","Chicken2","Chicken2","Chicken2","Chicken2","Chicken2","Chicken2","Chicken2","Chicken2","Chicken2","Chicken2"],["Hill-d2","Hill-d2","Hill-d2","Hill-d2","Hill-d2","Hill-d2","Hill-d1","Hill-d1","Hill-d1","Hill-d1","Hill-d1","Hill-d1","Hill-d1","evenness_pielou","evenness_pielou","evenness_pielou","evenness_pielou","evenness_pielou","evenness_pielou","Hill-d2","Hill-d2","Hill-d2","Hill-d2","Hill-d2","Hill-d2","Hill-d1","Hill-d1","Hill-d1","Hill-d1","Hill-d1","evenness_pielou","evenness_pielou","evenness_pielou","evenness_pielou","evenness_pielou","Hill-d1","evenness_pielou","Hill-d1","Hill-d2","Hill-d1","Hill-d1","Hill-d2","Hill-d1","Hill-d1","Hill-d2","Hill-d0","Hill-d0","Hill-d0","Hill-d0","Hill-d0","Hill-d0","Hill-d2","Hill-d2","Hill-d0","Hill-d0","Hill-d0","Hill-d0","Hill-d2","Hill-d2","evenness_pielou","Hill-d1","Hill-d1","Hill-d1","Hill-d0","Hill-d2","Hill-d2","Hill-d0","evenness_pielou","Hill-d0","Hill-d0","Hill-d1","Hill-d0","Hill-d2","Hill-d0","Hill-d0","Hill-d0","Hill-d0","evenness_pielou","Hill-d1","Hill-d1","Hill-d0","Hill-d2","Hill-d0","evenness_pielou","Hill-d1","Hill-d0","evenness_pielou","Hill-d2","evenness_pielou","Hill-d0","Hill-d1","Hill-d0","Hill-d2","Hill-d0","Hill-d0","Hill-d1"],["value","value","value","value","value","value","value","value","value","value","value","value","value","value","value","value","value","value","value","value","value","value","value","value","value","value","value","value","value","value","value","value","value","value","value","value","value","value","value","value","value","value","value","value","value","value","value","value","value","value","value","value","value","value","value","value","value","value","value","value","value","value","value","value","value","value","value","value","value","value","value","value","value","value","value","value","value","value","value","value","value","value","value","value","value","value","value","value","value","value","value","value","value","value","value","value"],["TR3_CTX","TR3_CTX","TR6_VAN","TR6_VAN","TR7_VAN+CCUG59168","TR7_VAN+CCUG59168","TR3_CTX","TR3_CTX","TR3_CTX","TR6_VAN","TR6_VAN","TR7_VAN+CCUG59168","TR7_VAN+CCUG59168","TR3_CTX","TR3_CTX","TR6_VAN","TR6_VAN","TR7_VAN+CCUG59168","TR7_VAN+CCUG59168","TR3_CTX","TR3_CTX","TR4_CCUG59168","TR4_CCUG59168","TR5_HV292.1","TR5_HV292.1","TR3_CTX","TR4_CCUG59168","TR4_CCUG59168","TR5_HV292.1","TR5_HV292.1","TR3_CTX","TR3_CTX","TR4_CCUG59168","TR5_HV292.1","TR5_HV292.1","TR3_CTX","TR4_CCUG59168","TR3_CTX","TR1_CTX+HV292.1","TR1_CTX+HV292.1","TR4_VAN","TR3_HV292.1","TR2_CTX","TR3_HV292.1","TR3_HV292.1","TR3_CTX","TR6_VAN","TR6_VAN","TR3_CTX","TR7_VAN+CCUG59168","TR7_VAN+CCUG59168","TR1_CTX+HV292.1","TR4_VAN","TR4_CCUG59168","TR5_HV292.1","TR4_CCUG59168","TR5_HV292.1","TR2_CTX","TR3_CTX","TR4_CCUG59168","TR1_CTX+HV292.1","TR3_HV292.1","TR5_VAN+CCUG59168","TR4_VAN","TR2_CTX","TR5_VAN+CCUG59168","TR3_HV292.1","TR3_CTX","TR1_CTX+HV292.1","TR3_CTX","TR2_CTX","TR2_CTX","TR4_CCUG59168","TR3_HV292.1","TR2_CTX","TR1_CTX+HV292.1","TR5_VAN+CCUG59168","TR6_VAN","TR2_CTX","TR4_CCUG59168","TR2_CTX","TR4_CCUG59168","TR2_CTX","CR_UNTREATED","TR6_VAN","TR1_CTX+HV292.1","TR5_HV292.1","TR6_VAN","TR4_CCUG59168","TR5_HV292.1","TR4_CCUG59168","TR4_CCUG59168","TR3_CTX","TR3_CTX","TR5_HV292.1","TR4_CCUG59168"],["CR_UNTREATED","TR2_CTX+HV292.1","CR_UNTREATED","TR2_CTX+HV292.1","CR_UNTREATED","TR2_CTX+HV292.1","TR6_VAN","CR_UNTREATED","TR2_CTX+HV292.1","CR_UNTREATED","TR2_CTX+HV292.1","CR_UNTREATED","TR2_CTX+HV292.1","CR_UNTREATED","TR2_CTX+HV292.1","CR_UNTREATED","TR2_CTX+HV292.1","CR_UNTREATED","TR2_CTX+HV292.1","TR4_CCUG59168","TR5_HV292.1","TR6_VAN","TR7_VAN+CCUG59168","TR6_VAN","TR7_VAN+CCUG59168","TR7_VAN+CCUG59168","TR6_VAN","TR7_VAN+CCUG59168","TR6_VAN","TR7_VAN+CCUG59168","TR4_CCUG59168","TR5_HV292.1","TR6_VAN","TR6_VAN","TR7_VAN+CCUG59168","TR5_HV292.1","TR7_VAN+CCUG59168","TR4_CCUG59168","TR4_VAN","TR4_VAN","TR6_CCUG59168","TR4_VAN","TR4_VAN","TR4_VAN","TR5_VAN+CCUG59168","TR6_VAN","TR2_CTX+HV292.1","CR_UNTREATED","TR7_VAN+CCUG59168","TR2_CTX+HV292.1","CR_UNTREATED","TR5_VAN+CCUG59168","TR6_CCUG59168","TR6_VAN","TR6_VAN","TR7_VAN+CCUG59168","TR7_VAN+CCUG59168","TR4_VAN","TR6_VAN","CR_UNTREATED","TR5_VAN+CCUG59168","TR5_VAN+CCUG59168","TR6_CCUG59168","TR6_CCUG59168","TR5_VAN+CCUG59168","TR6_CCUG59168","TR4_VAN","TR6_VAN","TR4_VAN","TR5_HV292.1","TR5_VAN+CCUG59168","TR4_VAN","CR_UNTREATED","TR5_VAN+CCUG59168","TR3_HV292.1","TR5_VAN+CCUG59168","TR6_CCUG59168","TR7_VAN+CCUG59168","TR3_HV292.1","CR_UNTREATED","TR5_VAN+CCUG59168","TR2_CTX+HV292.1","TR6_CCUG59168","TR2_CTX+HV292.1","TR7_VAN+CCUG59168","TR2_CTX","CR_UNTREATED","TR7_VAN+CCUG59168","TR2_CTX+HV292.1","TR2_CTX+HV292.1","TR2_CTX+HV292.1","TR5_HV292.1","TR7_VAN+CCUG59168","TR2_CTX+HV292.1","CR_UNTREATED","TR5_HV292.1"],[4.98546736263791e-08,4.98546736263791e-08,4.98546736263791e-08,4.98546736263791e-08,4.98546736263791e-08,4.98546736263791e-08,4.98546736263791e-08,4.98546736263791e-08,4.98546736263791e-08,4.98546736263791e-08,4.98546736263791e-08,4.98546736263791e-08,4.98546736263791e-08,4.98546736263791e-08,4.98546736263791e-08,4.98546736263791e-08,4.98546736263791e-08,4.98546736263791e-08,4.98546736263791e-08,9.97093472527582e-08,9.97093472527582e-08,9.97093472527582e-08,9.97093472527582e-08,9.97093472527582e-08,9.97093472527582e-08,9.97093472527582e-08,9.97093472527582e-08,9.97093472527582e-08,9.97093472527582e-08,9.97093472527582e-08,9.97093472527582e-08,9.97093472527582e-08,9.97093472527582e-08,9.97093472527582e-08,9.97093472527582e-08,1.99418694505516e-07,1.99418694505516e-07,3.98837389011033e-07,2.0708864429419e-07,2.0708864429419e-07,3.84593196546353e-07,7.39602301050679e-07,7.39602301050679e-07,7.39602301050679e-07,1.53837278618541e-06,7.27683527058832e-06,7.29786000337266e-06,7.31893413685682e-06,7.36123095205e-06,7.38245380704547e-06,7.40372640931371e-06,1.89447759780241e-06,2.69215237582447e-06,1.09423079039007e-05,1.09423079039007e-05,1.10786751740682e-05,1.10786751740682e-05,5.17721610735475e-06,2.53261742022006e-05,3.71915865252788e-05,1.38595992681334e-05,2.57677441686057e-05,2.67292271599715e-05,3.42152759318597e-05,3.73055400649963e-05,3.74978366632694e-05,4.48463390540862e-05,0.000101454260829681,6.52421937798617e-05,0.00035854873936163,0.000193834971059362,0.000296109076169129,0.00066117268163304,0.000657374375519621,0.000715973111844374,0.000798174571912088,0.00101018005481446,0.00248630242842115,0.00182977609279938,0.00395407387465538,0.00296567142042042,0.00471395881006865,0.00326371690167024,0.00674339300937766,0.0106001505611144,0.00802203260033275,0.0124673576524431,0.0162031179112886,0.0193253665564878,0.0210095897753113,0.0222209260007079,0.0238040467722237,0.0241243774397631,0.0324429700968352,0.0346435329606211,0.0386944983943234],[2.20410136032413e-07,2.20410136032413e-07,2.20410136032413e-07,2.20410136032413e-07,2.20410136032413e-07,2.20410136032413e-07,2.20410136032413e-07,2.20410136032413e-07,2.20410136032413e-07,2.20410136032413e-07,2.20410136032413e-07,2.20410136032413e-07,2.20410136032413e-07,2.20410136032413e-07,2.20410136032413e-07,2.20410136032413e-07,2.20410136032413e-07,2.20410136032413e-07,2.20410136032413e-07,2.3930243340662e-07,2.3930243340662e-07,2.3930243340662e-07,2.3930243340662e-07,2.3930243340662e-07,2.3930243340662e-07,2.3930243340662e-07,2.3930243340662e-07,2.3930243340662e-07,2.3930243340662e-07,2.3930243340662e-07,2.3930243340662e-07,2.3930243340662e-07,2.3930243340662e-07,2.3930243340662e-07,2.3930243340662e-07,4.52734333471983e-07,4.52734333471983e-07,8.81640544129652e-07,6.2126593288257e-06,6.2126593288257e-06,7.39602301050679e-06,7.39602301050679e-06,7.39602301050679e-06,7.39602301050679e-06,1.31860524530178e-05,1.41343867814171e-05,1.41343867814171e-05,1.41343867814171e-05,1.41343867814171e-05,1.41343867814171e-05,1.41343867814171e-05,1.4208581983518e-05,1.79476825054965e-05,1.93876815546194e-05,1.93876815546194e-05,1.93876815546194e-05,1.93876815546194e-05,3.10632966441285e-05,4.34162986323439e-05,6.24818653624684e-05,7.55978141898185e-05,0.000123365663815253,0.000123365663815253,0.00014061688748726,0.00014061688748726,0.00014061688748726,0.000158281196661481,0.000167101135484181,0.000217473979266206,0.000579194117430326,0.000612110434924301,0.000888327228507388,0.00104789632560708,0.00187821250148463,0.00195265394139375,0.00208219453542284,0.00252545013703616,0.00386758155532179,0.00439146262271851,0.00603894919038276,0.00684385712404711,0.00707093821510298,0.00725270422593386,0.0099376318032934,0.0153519421919587,0.0171900698578559,0.0177501363187326,0.022684365075804,0.0266119801761472,0.0284646055020347,0.0296279013342772,0.03117611853754,0.03117611853754,0.0412910528505175,0.0434336831148085,0.0477990862518113],["5.0e-08","5.0e-08","5.0e-08","5.0e-08","5.0e-08","5.0e-08","5.0e-08","5.0e-08","5.0e-08","5.0e-08","5.0e-08","5.0e-08","5.0e-08","5.0e-08","5.0e-08","5.0e-08","5.0e-08","5.0e-08","5.0e-08","1.0e-07","1.0e-07","1.0e-07","1.0e-07","1.0e-07","1.0e-07","1.0e-07","1.0e-07","1.0e-07","1.0e-07","1.0e-07","1.0e-07","1.0e-07","1.0e-07","1.0e-07","1.0e-07","2.0e-07","2.0e-07","4.0e-07","2.1e-07","2.1e-07","3.8e-07","7.4e-07","7.4e-07","7.4e-07","1.5e-06","7.3e-06","7.3e-06","7.3e-06","7.4e-06","7.4e-06","7.4e-06","1.9e-06","2.7e-06","1.1e-05","1.1e-05","1.1e-05","1.1e-05","5.2e-06","2.5e-05","3.7e-05","1.4e-05","2.6e-05","2.7e-05","3.4e-05","3.7e-05","3.7e-05","4.5e-05","0.00010","6.5e-05","0.00036","0.00019","0.00030","0.00066","0.00066","0.00072","0.00080","0.00101","0.00249","0.00183","0.00395","0.00297","0.00471","0.00326","0.00674","0.01060","0.00802","0.01247","0.01620","0.01933","0.02101","0.02222","0.02380","0.02412","0.03244","0.03464","0.03869"],["****","****","****","****","****","****","****","****","****","****","****","****","****","****","****","****","****","****","****","****","****","****","****","****","****","****","****","****","****","****","****","****","****","****","****","****","****","****","****","****","****","****","****","****","****","****","****","****","****","****","****","****","****","****","****","****","****","****","****","****","****","****","****","****","****","****","****","***","****","***","***","***","***","***","***","***","**","**","**","**","**","**","**","**","*","**","*","*","*","*","*","*","*","*","*","*"],["Wilcoxon","Wilcoxon","Wilcoxon","Wilcoxon","Wilcoxon","Wilcoxon","Wilcoxon","Wilcoxon","Wilcoxon","Wilcoxon","Wilcoxon","Wilcoxon","Wilcoxon","Wilcoxon","Wilcoxon","Wilcoxon","Wilcoxon","Wilcoxon","Wilcoxon","Wilcoxon","Wilcoxon","Wilcoxon","Wilcoxon","Wilcoxon","Wilcoxon","Wilcoxon","Wilcoxon","Wilcoxon","Wilcoxon","Wilcoxon","Wilcoxon","Wilcoxon","Wilcoxon","Wilcoxon","Wilcoxon","Wilcoxon","Wilcoxon","Wilcoxon","Wilcoxon","Wilcoxon","Wilcoxon","Wilcoxon","Wilcoxon","Wilcoxon","Wilcoxon","Wilcoxon","Wilcoxon","Wilcoxon","Wilcoxon","Wilcoxon","Wilcoxon","Wilcoxon","Wilcoxon","Wilcoxon","Wilcoxon","Wilcoxon","Wilcoxon","Wilcoxon","Wilcoxon","Wilcoxon","Wilcoxon","Wilcoxon","Wilcoxon","Wilcoxon","Wilcoxon","Wilcoxon","Wilcoxon","Wilcoxon","Wilcoxon","Wilcoxon","Wilcoxon","Wilcoxon","Wilcoxon","Wilcoxon","Wilcoxon","Wilcoxon","Wilcoxon","Wilcoxon","Wilcoxon","Wilcoxon","Wilcoxon","Wilcoxon","Wilcoxon","Wilcoxon","Wilcoxon","Wilcoxon","Wilcoxon","Wilcoxon","Wilcoxon","Wilcoxon","Wilcoxon","Wilcoxon","Wilcoxon","Wilcoxon","Wilcoxon","Wilcoxon"]],"container":"<table class=\"display\">\n  <thead>\n    <tr>\n      <th> <\/th>\n      <th>Model2<\/th>\n      <th>variable<\/th>\n      <th>.y.<\/th>\n      <th>group1<\/th>\n      <th>group2<\/th>\n      <th>p<\/th>\n      <th>p.adj<\/th>\n      <th>p.format<\/th>\n      <th>p.signif<\/th>\n      <th>method<\/th>\n    <\/tr>\n  <\/thead>\n<\/table>","options":{"columnDefs":[{"className":"dt-right","targets":[6,7]},{"orderable":false,"targets":0}],"order":[],"autoWidth":false,"orderClasses":false}},"evals":[],"jsHooks":[]}</script>
```

treatment: significant differences. ToConsider: pariwise heatmap visualisation?


```r
alpha_plot_noleg$data %>% 
  filter(Day_of_Treatment >= 1 & Day_of_Treatment <= 6) %>% 
  # group_by(Model2) %>% 
  # group_by(Model2, variable) %>%
  # rstatix::wilcox_test(log10_ratio ~ Treatment_A,
  #                      # comparisons = list(c("Inulin_3", "Control")),
  #                      # ref.group = "Control",
  #                      data = . ) %>%  #,
  ggpubr::compare_means(value ~ Reactor_Treatment_Dose_fermentation,
                        data = .,
                        method = "wilcox.test",
                        group.by = c("Model2", "variable"),
                        p.adjust.method = "none") %>% 
  ungroup() %>% 
  group_by(Model2) %>% 
  rstatix::adjust_pvalue(method = "fdr",p.col = "p", output.col = "p.adj") %>% 
  arrange(p.adj) -> p_tmp

p_tmp %>% 
  xlsx::write.xlsx(file = out_xlsx,
                   append = TRUE,
                   sheetName = "alpha_div_comp_d1-d6",
                   showNA = TRUE)


p_tmp %>% 
  filter(p.adj <= 0.05) %>% 
  DT::datatable()
```

```{=html}
<div id="htmlwidget-e9bc75b139554114a9d5" style="width:100%;height:auto;" class="datatables html-widget"></div>
<script type="application/json" data-for="htmlwidget-e9bc75b139554114a9d5">{"x":{"filter":"none","vertical":false,"data":[["1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","23","24","25","26","27","28","29","30","31","32","33","34","35","36","37","38","39","40","41","42","43","44","45","46","47","48","49","50","51","52","53","54","55","56","57","58","59","60","61","62","63","64","65","66","67","68","69","70","71","72","73","74","75","76","77","78"],["Chicken2","Chicken2","Chicken2","Chicken2","Chicken2","Chicken2","Chicken2","Chicken2","Chicken2","Chicken2","Chicken2","Chicken2","Chicken2","Chicken2","Chicken2","Chicken2","Chicken2","Chicken2","Chicken2","Chicken2","Chicken2","Chicken2","Chicken2","Chicken2","Chicken2","Chicken2","Chicken2","Chicken2","Chicken2","Chicken2","Chicken2","Chicken2","Chicken2","Chicken2","Chicken2","Chicken2","Chicken2","Chicken2","Chicken2","Chicken2","Chicken2","Chicken2","Chicken2","Chicken2","Chicken2","Chicken2","Chicken2","Chicken2","Chicken2","Chicken2","Chicken2","Chicken2","Chicken2","Chicken2","Chicken2","Chicken2","Chicken2","Chicken2","Chicken2","Chicken2","Chicken1","Chicken1","Chicken1","Chicken1","Chicken1","Chicken1","Chicken1","Chicken1","Chicken1","Chicken1","Chicken1","Chicken1","Chicken1","Chicken1","Chicken1","Chicken1","Chicken1","Chicken2"],["Hill-d2","Hill-d2","Hill-d2","Hill-d2","Hill-d2","Hill-d2","Hill-d2","Hill-d2","Hill-d2","Hill-d1","Hill-d1","Hill-d1","Hill-d1","Hill-d1","Hill-d1","Hill-d1","Hill-d1","Hill-d1","Hill-d1","evenness_pielou","evenness_pielou","evenness_pielou","evenness_pielou","evenness_pielou","evenness_pielou","evenness_pielou","evenness_pielou","evenness_pielou","Hill-d0","Hill-d0","Hill-d0","Hill-d0","Hill-d0","Hill-d0","Hill-d0","Hill-d0","Hill-d2","Hill-d2","Hill-d2","Hill-d2","Hill-d2","Hill-d1","Hill-d1","Hill-d1","Hill-d1","Hill-d1","Hill-d1","evenness_pielou","evenness_pielou","evenness_pielou","evenness_pielou","evenness_pielou","Hill-d0","Hill-d0","Hill-d0","Hill-d0","Hill-d0","Hill-d2","evenness_pielou","evenness_pielou","Hill-d0","Hill-d2","Hill-d2","Hill-d2","Hill-d2","Hill-d2","Hill-d1","Hill-d1","Hill-d1","Hill-d1","evenness_pielou","evenness_pielou","evenness_pielou","evenness_pielou","evenness_pielou","Hill-d0","evenness_pielou","evenness_pielou"],["value","value","value","value","value","value","value","value","value","value","value","value","value","value","value","value","value","value","value","value","value","value","value","value","value","value","value","value","value","value","value","value","value","value","value","value","value","value","value","value","value","value","value","value","value","value","value","value","value","value","value","value","value","value","value","value","value","value","value","value","value","value","value","value","value","value","value","value","value","value","value","value","value","value","value","value","value","value"],["TR3_CTX20_NA","TR3_CTX20_NA","TR3_CTX20_NA","TR4_CCUG59168_NA","TR4_CCUG59168_NA","TR6_VAN90_NA","TR6_VAN90_NA","TR7_VAN+CCUG5916890_NA","TR7_VAN+CCUG5916890_NA","TR3_CTX20_NA","TR3_CTX20_NA","TR3_CTX20_NA","TR3_CTX20_NA","TR4_CCUG59168_NA","TR4_CCUG59168_NA","TR6_VAN90_NA","TR6_VAN90_NA","TR7_VAN+CCUG5916890_NA","TR7_VAN+CCUG5916890_NA","TR3_CTX20_NA","TR3_CTX20_NA","TR3_CTX20_NA","TR4_CCUG59168_NA","TR4_CCUG59168_NA","TR6_VAN90_NA","TR6_VAN90_NA","TR7_VAN+CCUG5916890_NA","TR7_VAN+CCUG5916890_NA","TR3_CTX20_NA","TR3_CTX20_NA","TR4_CCUG59168_NA","TR4_CCUG59168_NA","TR6_VAN90_NA","TR6_VAN90_NA","TR7_VAN+CCUG5916890_NA","TR7_VAN+CCUG5916890_NA","TR3_CTX20_NA","TR4_CCUG59168_NA","TR5_HV292.1_NA","TR5_HV292.1_NA","TR5_HV292.1_NA","TR3_CTX20_NA","TR3_CTX20_NA","TR4_CCUG59168_NA","TR5_HV292.1_NA","TR5_HV292.1_NA","TR5_HV292.1_NA","TR3_CTX20_NA","TR5_HV292.1_NA","TR5_HV292.1_NA","TR5_HV292.1_NA","CR_UNTREATED_NA","TR3_CTX20_NA","TR4_CCUG59168_NA","TR5_HV292.1_NA","TR5_HV292.1_NA","TR5_HV292.1_NA","TR4_CCUG59168_NA","TR3_CTX20_NA","TR4_CCUG59168_NA","TR3_HV292.1_NA","TR1_CTX+HV292.120_NA","TR2_CTX20_NA","TR3_HV292.1_NA","TR3_HV292.1_NA","TR4_VAN90_NA","TR1_CTX+HV292.120_NA","TR2_CTX20_NA","TR3_HV292.1_NA","TR4_VAN90_NA","TR1_CTX+HV292.120_NA","TR1_CTX+HV292.120_NA","TR3_HV292.1_NA","TR3_HV292.1_NA","TR4_VAN90_NA","TR4_VAN90_NA","TR5_VAN+CCUG5916890_NA","TR4_CCUG59168_NA"],["TR4_CCUG59168_NA","CR_UNTREATED_NA","TR2_CTX+HV292.120_NA","TR6_VAN90_NA","TR7_VAN+CCUG5916890_NA","CR_UNTREATED_NA","TR2_CTX+HV292.120_NA","CR_UNTREATED_NA","TR2_CTX+HV292.120_NA","TR4_CCUG59168_NA","TR6_VAN90_NA","CR_UNTREATED_NA","TR2_CTX+HV292.120_NA","TR6_VAN90_NA","TR7_VAN+CCUG5916890_NA","CR_UNTREATED_NA","TR2_CTX+HV292.120_NA","CR_UNTREATED_NA","TR2_CTX+HV292.120_NA","TR4_CCUG59168_NA","CR_UNTREATED_NA","TR2_CTX+HV292.120_NA","TR6_VAN90_NA","TR7_VAN+CCUG5916890_NA","CR_UNTREATED_NA","TR2_CTX+HV292.120_NA","CR_UNTREATED_NA","TR2_CTX+HV292.120_NA","TR6_VAN90_NA","TR7_VAN+CCUG5916890_NA","TR6_VAN90_NA","TR7_VAN+CCUG5916890_NA","CR_UNTREATED_NA","TR2_CTX+HV292.120_NA","CR_UNTREATED_NA","TR2_CTX+HV292.120_NA","TR5_HV292.1_NA","TR5_HV292.1_NA","TR6_VAN90_NA","TR7_VAN+CCUG5916890_NA","TR2_CTX+HV292.120_NA","TR5_HV292.1_NA","TR7_VAN+CCUG5916890_NA","TR5_HV292.1_NA","TR6_VAN90_NA","TR7_VAN+CCUG5916890_NA","TR2_CTX+HV292.120_NA","TR5_HV292.1_NA","TR6_VAN90_NA","TR7_VAN+CCUG5916890_NA","TR2_CTX+HV292.120_NA","TR2_CTX+HV292.120_NA","TR5_HV292.1_NA","TR5_HV292.1_NA","TR6_VAN90_NA","TR7_VAN+CCUG5916890_NA","CR_UNTREATED_NA","CR_UNTREATED_NA","TR6_VAN90_NA","TR5_HV292.1_NA","TR4_VAN90_NA","TR4_VAN90_NA","TR4_VAN90_NA","TR4_VAN90_NA","TR5_VAN+CCUG5916890_NA","TR6_CCUG59168_NA","TR4_VAN90_NA","TR4_VAN90_NA","TR4_VAN90_NA","TR6_CCUG59168_NA","TR4_VAN90_NA","TR5_VAN+CCUG5916890_NA","TR4_VAN90_NA","TR5_VAN+CCUG5916890_NA","TR6_CCUG59168_NA","TR6_CCUG59168_NA","TR6_CCUG59168_NA","CR_UNTREATED_NA"],[0.00216450216450216,0.00216450216450216,0.00216450216450216,0.00216450216450216,0.00216450216450216,0.00216450216450216,0.00216450216450216,0.00216450216450216,0.00216450216450216,0.00216450216450216,0.00216450216450216,0.00216450216450216,0.00216450216450216,0.00216450216450216,0.00216450216450216,0.00216450216450216,0.00216450216450216,0.00216450216450216,0.00216450216450216,0.00216450216450216,0.00216450216450216,0.00216450216450216,0.00216450216450216,0.00216450216450216,0.00216450216450216,0.00216450216450216,0.00216450216450216,0.00216450216450216,0.00492203567532315,0.00492203567532315,0.00499812476508246,0.00499812476508246,0.00499812476508246,0.00499812476508246,0.00499812476508246,0.00499812476508246,0.00432900432900433,0.00432900432900433,0.00432900432900433,0.00432900432900433,0.00432900432900433,0.00432900432900433,0.00432900432900433,0.00432900432900433,0.00432900432900433,0.00432900432900433,0.00432900432900433,0.00432900432900433,0.00432900432900433,0.00432900432900433,0.00432900432900433,0.00432900432900433,0.00796941302412179,0.00865800865800866,0.00796941302412179,0.00796941302412179,0.00865800865800866,0.00865800865800866,0.00865800865800866,0.00865800865800866,0.00432900432900433,0.00432900432900433,0.00432900432900433,0.00432900432900433,0.00432900432900433,0.00432900432900433,0.00432900432900433,0.00432900432900433,0.00432900432900433,0.00432900432900433,0.00432900432900433,0.00432900432900433,0.00432900432900433,0.00432900432900433,0.00432900432900433,0.00865800865800866,0.00865800865800866,0.025974025974026],[0.00649350649350649,0.00649350649350649,0.00649350649350649,0.00649350649350649,0.00649350649350649,0.00649350649350649,0.00649350649350649,0.00649350649350649,0.00649350649350649,0.00649350649350649,0.00649350649350649,0.00649350649350649,0.00649350649350649,0.00649350649350649,0.00649350649350649,0.00649350649350649,0.00649350649350649,0.00649350649350649,0.00649350649350649,0.00649350649350649,0.00649350649350649,0.00649350649350649,0.00649350649350649,0.00649350649350649,0.00649350649350649,0.00649350649350649,0.00649350649350649,0.00649350649350649,0.00807389385128705,0.00807389385128705,0.00807389385128705,0.00807389385128705,0.00807389385128705,0.00807389385128705,0.00807389385128705,0.00807389385128705,0.00807389385128705,0.00807389385128705,0.00807389385128705,0.00807389385128705,0.00807389385128705,0.00807389385128705,0.00807389385128705,0.00807389385128705,0.00807389385128705,0.00807389385128705,0.00807389385128705,0.00807389385128705,0.00807389385128705,0.00807389385128705,0.00807389385128705,0.00807389385128705,0.0121212121212121,0.0121212121212121,0.0121212121212121,0.0121212121212121,0.0121212121212121,0.0121212121212121,0.0121212121212121,0.0121212121212121,0.0173160173160173,0.0173160173160173,0.0173160173160173,0.0173160173160173,0.0173160173160173,0.0173160173160173,0.0173160173160173,0.0173160173160173,0.0173160173160173,0.0173160173160173,0.0173160173160173,0.0173160173160173,0.0173160173160173,0.0173160173160173,0.0173160173160173,0.0305576776165011,0.0305576776165011,0.0357675111773472],["0.0022","0.0022","0.0022","0.0022","0.0022","0.0022","0.0022","0.0022","0.0022","0.0022","0.0022","0.0022","0.0022","0.0022","0.0022","0.0022","0.0022","0.0022","0.0022","0.0022","0.0022","0.0022","0.0022","0.0022","0.0022","0.0022","0.0022","0.0022","0.0049","0.0049","0.0050","0.0050","0.0050","0.0050","0.0050","0.0050","0.0043","0.0043","0.0043","0.0043","0.0043","0.0043","0.0043","0.0043","0.0043","0.0043","0.0043","0.0043","0.0043","0.0043","0.0043","0.0043","0.0080","0.0087","0.0080","0.0080","0.0087","0.0087","0.0087","0.0087","0.0043","0.0043","0.0043","0.0043","0.0043","0.0043","0.0043","0.0043","0.0043","0.0043","0.0043","0.0043","0.0043","0.0043","0.0043","0.0087","0.0087","0.0260"],["**","**","**","**","**","**","**","**","**","**","**","**","**","**","**","**","**","**","**","**","**","**","**","**","**","**","**","**","**","**","**","**","**","**","**","**","**","**","**","**","**","**","**","**","**","**","**","**","**","**","**","**","**","**","**","**","**","**","**","**","**","**","**","**","**","**","**","**","**","**","**","**","**","**","**","**","**","*"],["Wilcoxon","Wilcoxon","Wilcoxon","Wilcoxon","Wilcoxon","Wilcoxon","Wilcoxon","Wilcoxon","Wilcoxon","Wilcoxon","Wilcoxon","Wilcoxon","Wilcoxon","Wilcoxon","Wilcoxon","Wilcoxon","Wilcoxon","Wilcoxon","Wilcoxon","Wilcoxon","Wilcoxon","Wilcoxon","Wilcoxon","Wilcoxon","Wilcoxon","Wilcoxon","Wilcoxon","Wilcoxon","Wilcoxon","Wilcoxon","Wilcoxon","Wilcoxon","Wilcoxon","Wilcoxon","Wilcoxon","Wilcoxon","Wilcoxon","Wilcoxon","Wilcoxon","Wilcoxon","Wilcoxon","Wilcoxon","Wilcoxon","Wilcoxon","Wilcoxon","Wilcoxon","Wilcoxon","Wilcoxon","Wilcoxon","Wilcoxon","Wilcoxon","Wilcoxon","Wilcoxon","Wilcoxon","Wilcoxon","Wilcoxon","Wilcoxon","Wilcoxon","Wilcoxon","Wilcoxon","Wilcoxon","Wilcoxon","Wilcoxon","Wilcoxon","Wilcoxon","Wilcoxon","Wilcoxon","Wilcoxon","Wilcoxon","Wilcoxon","Wilcoxon","Wilcoxon","Wilcoxon","Wilcoxon","Wilcoxon","Wilcoxon","Wilcoxon","Wilcoxon"]],"container":"<table class=\"display\">\n  <thead>\n    <tr>\n      <th> <\/th>\n      <th>Model2<\/th>\n      <th>variable<\/th>\n      <th>.y.<\/th>\n      <th>group1<\/th>\n      <th>group2<\/th>\n      <th>p<\/th>\n      <th>p.adj<\/th>\n      <th>p.format<\/th>\n      <th>p.signif<\/th>\n      <th>method<\/th>\n    <\/tr>\n  <\/thead>\n<\/table>","options":{"columnDefs":[{"className":"dt-right","targets":[6,7]},{"orderable":false,"targets":0}],"order":[],"autoWidth":false,"orderClasses":false}},"evals":[],"jsHooks":[]}</script>
```


```r
alpha_plot_noleg$data %>% 
  filter(Day_of_Treatment > 7) %>% 
  # group_by(Model2) %>% 
  # group_by(Model2, variable) %>%
  # rstatix::wilcox_test(log10_ratio ~ Treatment_A,
  #                      # comparisons = list(c("Inulin_3", "Control")),
  #                      # ref.group = "Control",
  #                      data = . ) %>%  #,
  ggpubr::compare_means(value ~ Reactor_Treatment_Dose_fermentation,
                        data = .,
                        method = "wilcox.test",
                        group.by = c("Model2", "variable"),
                        p.adjust.method = "none") %>% 
  arrange(p) %>% 
  ungroup() %>% 
  group_by(Model2) %>% 
  rstatix::adjust_pvalue(method = "fdr",p.col = "p", output.col = "p.adj") %>% 
  arrange(p.adj) -> p_tmp

p_tmp %>% 
  xlsx::write.xlsx(file = out_xlsx,
                   append = TRUE,
                   sheetName = "alpha_div_comp_d7-end",
                   showNA = TRUE)


p_tmp %>% 
  filter(p.adj <= 0.05) %>% 
  DT::datatable()
```

```{=html}
<div id="htmlwidget-ec72d682becb25d7f71c" style="width:100%;height:auto;" class="datatables html-widget"></div>
<script type="application/json" data-for="htmlwidget-ec72d682becb25d7f71c">{"x":{"filter":"none","vertical":false,"data":[["1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","23","24","25","26","27","28","29","30","31","32","33","34","35","36","37","38","39","40","41","42","43","44","45","46","47","48","49","50","51","52","53","54","55","56","57","58","59","60","61","62","63","64","65","66","67","68","69","70","71","72","73","74","75","76","77","78","79","80","81","82","83","84","85","86"],["Chicken2","Chicken2","Chicken2","Chicken2","Chicken2","Chicken2","Chicken2","Chicken2","Chicken2","Chicken2","Chicken2","Chicken2","Chicken2","Chicken2","Chicken2","Chicken2","Chicken2","Chicken2","Chicken2","Chicken2","Chicken2","Chicken2","Chicken2","Chicken2","Chicken2","Chicken2","Chicken2","Chicken2","Chicken2","Chicken2","Chicken2","Chicken2","Chicken2","Chicken2","Chicken2","Chicken2","Chicken2","Chicken2","Chicken2","Chicken2","Chicken2","Chicken2","Chicken2","Chicken2","Chicken2","Chicken2","Chicken2","Chicken2","Chicken2","Chicken2","Chicken1","Chicken1","Chicken1","Chicken1","Chicken1","Chicken1","Chicken1","Chicken1","Chicken1","Chicken1","Chicken1","Chicken1","Chicken1","Chicken1","Chicken1","Chicken1","Chicken1","Chicken1","Chicken2","Chicken2","Chicken1","Chicken1","Chicken2","Chicken1","Chicken1","Chicken1","Chicken1","Chicken2","Chicken2","Chicken2","Chicken1","Chicken2","Chicken2","Chicken1","Chicken1","Chicken1"],["Hill-d2","Hill-d2","Hill-d2","Hill-d2","Hill-d2","Hill-d2","Hill-d2","Hill-d2","Hill-d2","Hill-d2","Hill-d1","Hill-d1","Hill-d1","Hill-d1","Hill-d1","Hill-d1","Hill-d1","Hill-d1","Hill-d1","Hill-d1","evenness_pielou","evenness_pielou","evenness_pielou","evenness_pielou","evenness_pielou","evenness_pielou","evenness_pielou","evenness_pielou","evenness_pielou","evenness_pielou","Hill-d2","Hill-d2","Hill-d2","Hill-d2","Hill-d1","Hill-d1","Hill-d1","Hill-d1","evenness_pielou","evenness_pielou","Hill-d0","Hill-d0","Hill-d0","Hill-d0","Hill-d0","Hill-d0","Hill-d0","Hill-d0","evenness_pielou","evenness_pielou","Hill-d2","Hill-d2","Hill-d1","Hill-d1","Hill-d2","Hill-d2","Hill-d1","Hill-d1","Hill-d2","Hill-d2","Hill-d2","Hill-d2","Hill-d1","Hill-d1","Hill-d1","Hill-d1","Hill-d0","Hill-d0","Hill-d0","Hill-d0","Hill-d0","Hill-d0","evenness_pielou","Hill-d0","Hill-d0","Hill-d0","Hill-d0","Hill-d2","evenness_pielou","Hill-d1","Hill-d0","Hill-d0","evenness_pielou","Hill-d0","Hill-d0","Hill-d1"],["value","value","value","value","value","value","value","value","value","value","value","value","value","value","value","value","value","value","value","value","value","value","value","value","value","value","value","value","value","value","value","value","value","value","value","value","value","value","value","value","value","value","value","value","value","value","value","value","value","value","value","value","value","value","value","value","value","value","value","value","value","value","value","value","value","value","value","value","value","value","value","value","value","value","value","value","value","value","value","value","value","value","value","value","value","value"],["TR5_HV292.1_NA","TR5_HV292.1_NA","TR5_HV292.1_NA","TR6_VAN90_NA","TR6_VAN90_NA","TR6_VAN90_NA","TR7_VAN+CCUG5916890_NA","TR7_VAN+CCUG5916890_NA","CR_UNTREATED_NA","TR2_CTX+HV292.120_NA","TR5_HV292.1_NA","TR5_HV292.1_NA","TR6_VAN90_NA","TR6_VAN90_NA","TR6_VAN90_NA","TR7_VAN+CCUG5916890_NA","TR7_VAN+CCUG5916890_NA","TR7_VAN+CCUG5916890_NA","CR_UNTREATED_NA","TR2_CTX+HV292.120_NA","TR5_HV292.1_NA","TR5_HV292.1_NA","TR5_HV292.1_NA","TR6_VAN90_NA","TR6_VAN90_NA","TR6_VAN90_NA","TR7_VAN+CCUG5916890_NA","TR7_VAN+CCUG5916890_NA","CR_UNTREATED_NA","TR2_CTX+HV292.120_NA","TR6_VAN90_NA","TR7_VAN+CCUG5916890_NA","TR7_VAN+CCUG5916890_NA","TR3_CTX20_NA","TR5_HV292.1_NA","TR6_VAN90_NA","TR7_VAN+CCUG5916890_NA","TR3_CTX20_NA","TR6_VAN90_NA","TR3_CTX20_NA","TR6_VAN90_NA","TR7_VAN+CCUG5916890_NA","TR5_HV292.1_NA","TR5_HV292.1_NA","TR6_VAN90_NA","TR6_VAN90_NA","TR7_VAN+CCUG5916890_NA","TR7_VAN+CCUG5916890_NA","TR5_HV292.1_NA","TR7_VAN+CCUG5916890_NA","TR1_CTX+HV292.120_NA","TR1_CTX+HV292.120_NA","TR1_CTX+HV292.120_NA","TR1_CTX+HV292.120_NA","TR4_VAN90_NA","TR5_VAN+CCUG5916890_NA","TR4_VAN90_NA","TR5_VAN+CCUG5916890_NA","TR2_CTX20_NA","TR2_CTX20_NA","TR3_HV292.1_NA","TR3_HV292.1_NA","TR2_CTX20_NA","TR2_CTX20_NA","TR3_HV292.1_NA","TR3_HV292.1_NA","TR1_CTX+HV292.120_NA","TR1_CTX+HV292.120_NA","TR6_VAN90_NA","TR7_VAN+CCUG5916890_NA","TR5_VAN+CCUG5916890_NA","TR4_VAN90_NA","CR_UNTREATED_NA","TR2_CTX20_NA","TR3_HV292.1_NA","TR2_CTX20_NA","TR3_HV292.1_NA","TR6_VAN90_NA","TR6_VAN90_NA","TR6_VAN90_NA","TR1_CTX+HV292.120_NA","TR5_HV292.1_NA","TR2_CTX+HV292.120_NA","TR2_CTX20_NA","TR2_CTX20_NA","TR1_CTX+HV292.120_NA"],["TR6_VAN90_NA","TR7_VAN+CCUG5916890_NA","TR3_CTX20_NA","CR_UNTREATED_NA","TR2_CTX+HV292.120_NA","TR3_CTX20_NA","CR_UNTREATED_NA","TR2_CTX+HV292.120_NA","TR3_CTX20_NA","TR3_CTX20_NA","TR6_VAN90_NA","TR7_VAN+CCUG5916890_NA","CR_UNTREATED_NA","TR2_CTX+HV292.120_NA","TR3_CTX20_NA","CR_UNTREATED_NA","TR2_CTX+HV292.120_NA","TR3_CTX20_NA","TR3_CTX20_NA","TR3_CTX20_NA","TR6_VAN90_NA","TR7_VAN+CCUG5916890_NA","TR3_CTX20_NA","CR_UNTREATED_NA","TR2_CTX+HV292.120_NA","TR3_CTX20_NA","CR_UNTREATED_NA","TR2_CTX+HV292.120_NA","TR3_CTX20_NA","TR3_CTX20_NA","TR4_CCUG59168_NA","TR3_CTX20_NA","TR4_CCUG59168_NA","TR4_CCUG59168_NA","TR3_CTX20_NA","TR4_CCUG59168_NA","TR4_CCUG59168_NA","TR4_CCUG59168_NA","TR4_CCUG59168_NA","TR4_CCUG59168_NA","TR3_CTX20_NA","TR3_CTX20_NA","TR6_VAN90_NA","TR7_VAN+CCUG5916890_NA","CR_UNTREATED_NA","TR2_CTX+HV292.120_NA","CR_UNTREATED_NA","TR2_CTX+HV292.120_NA","CR_UNTREATED_NA","TR4_CCUG59168_NA","TR4_VAN90_NA","TR5_VAN+CCUG5916890_NA","TR4_VAN90_NA","TR5_VAN+CCUG5916890_NA","TR6_CCUG59168_NA","TR6_CCUG59168_NA","TR6_CCUG59168_NA","TR6_CCUG59168_NA","TR4_VAN90_NA","TR5_VAN+CCUG5916890_NA","TR4_VAN90_NA","TR5_VAN+CCUG5916890_NA","TR4_VAN90_NA","TR5_VAN+CCUG5916890_NA","TR4_VAN90_NA","TR5_VAN+CCUG5916890_NA","TR5_VAN+CCUG5916890_NA","TR4_VAN90_NA","TR4_CCUG59168_NA","TR4_CCUG59168_NA","TR6_CCUG59168_NA","TR6_CCUG59168_NA","TR4_CCUG59168_NA","TR5_VAN+CCUG5916890_NA","TR5_VAN+CCUG5916890_NA","TR4_VAN90_NA","TR4_VAN90_NA","TR7_VAN+CCUG5916890_NA","TR7_VAN+CCUG5916890_NA","TR7_VAN+CCUG5916890_NA","TR2_CTX20_NA","TR3_CTX20_NA","TR4_CCUG59168_NA","TR6_CCUG59168_NA","TR3_HV292.1_NA","TR2_CTX20_NA"],[0.000582750582750583,0.000582750582750583,0.000582750582750583,0.000582750582750583,0.000582750582750583,0.000582750582750583,0.000582750582750583,0.000582750582750583,0.000582750582750583,0.000582750582750583,0.000582750582750583,0.000582750582750583,0.000582750582750583,0.000582750582750583,0.000582750582750583,0.000582750582750583,0.000582750582750583,0.000582750582750583,0.000582750582750583,0.000582750582750583,0.000582750582750583,0.000582750582750583,0.000582750582750583,0.000582750582750583,0.000582750582750583,0.000582750582750583,0.000582750582750583,0.000582750582750583,0.000582750582750583,0.000582750582750583,0.00116550116550117,0.00116550116550117,0.00116550116550117,0.00116550116550117,0.00116550116550117,0.00116550116550117,0.00116550116550117,0.00116550116550117,0.00116550116550117,0.00116550116550117,0.00211656160257058,0.00211656160257058,0.00214070973158512,0.00214070973158512,0.00214070973158512,0.00214070973158512,0.00214070973158512,0.00214070973158512,0.00233100233100233,0.00233100233100233,0.000310800310800311,0.000310800310800311,0.000310800310800311,0.000310800310800311,0.000582750582750583,0.000582750582750583,0.000582750582750583,0.000582750582750583,0.00116550116550117,0.00116550116550117,0.00116550116550117,0.00116550116550117,0.00116550116550117,0.00116550116550117,0.00116550116550117,0.00116550116550117,0.00138956391640281,0.00141749107279779,0.00331778042978762,0.00331778042978762,0.00204514318801797,0.00209258466211306,0.00466200466200466,0.00314691200925023,0.00318912451057108,0.00323167322080157,0.00327455843386829,0.00699300699300699,0.00699300699300699,0.0110722610722611,0.01172847850548,0.0213077106716282,0.0221445221445221,0.0182558605874593,0.0193733848500302,0.01998001998002],[0.00163170163170163,0.00163170163170163,0.00163170163170163,0.00163170163170163,0.00163170163170163,0.00163170163170163,0.00163170163170163,0.00163170163170163,0.00163170163170163,0.00163170163170163,0.00163170163170163,0.00163170163170163,0.00163170163170163,0.00163170163170163,0.00163170163170163,0.00163170163170163,0.00163170163170163,0.00163170163170163,0.00163170163170163,0.00163170163170163,0.00163170163170163,0.00163170163170163,0.00163170163170163,0.00163170163170163,0.00163170163170163,0.00163170163170163,0.00163170163170163,0.00163170163170163,0.00163170163170163,0.00163170163170163,0.00244755244755245,0.00244755244755245,0.00244755244755245,0.00244755244755245,0.00244755244755245,0.00244755244755245,0.00244755244755245,0.00244755244755245,0.00244755244755245,0.00244755244755245,0.00374624203027396,0.00374624203027396,0.00374624203027396,0.00374624203027396,0.00374624203027396,0.00374624203027396,0.00374624203027396,0.00374624203027396,0.00391608391608392,0.00391608391608392,0.00437062937062937,0.00437062937062937,0.00437062937062937,0.00437062937062937,0.00437062937062937,0.00437062937062937,0.00437062937062937,0.00437062937062937,0.00437062937062937,0.00437062937062937,0.00437062937062937,0.00437062937062937,0.00437062937062937,0.00437062937062937,0.00437062937062937,0.00437062937062937,0.00472497024265931,0.00472497024265931,0.00535949146350309,0.00535949146350309,0.00627775398633918,0.00627775398633918,0.00738883757751682,0.00818639608467073,0.00818639608467073,0.00818639608467073,0.00818639608467073,0.010680228862047,0.010680228862047,0.0166083916083916,0.028148348413152,0.0314008367792416,0.0320713768989631,0.042128909047983,0.0428143285286142,0.0428143285286142],["0.00058","0.00058","0.00058","0.00058","0.00058","0.00058","0.00058","0.00058","0.00058","0.00058","0.00058","0.00058","0.00058","0.00058","0.00058","0.00058","0.00058","0.00058","0.00058","0.00058","0.00058","0.00058","0.00058","0.00058","0.00058","0.00058","0.00058","0.00058","0.00058","0.00058","0.00117","0.00117","0.00117","0.00117","0.00117","0.00117","0.00117","0.00117","0.00117","0.00117","0.00212","0.00212","0.00214","0.00214","0.00214","0.00214","0.00214","0.00214","0.00233","0.00233","0.00031","0.00031","0.00031","0.00031","0.00058","0.00058","0.00058","0.00058","0.00117","0.00117","0.00117","0.00117","0.00117","0.00117","0.00117","0.00117","0.00139","0.00142","0.00332","0.00332","0.00205","0.00209","0.00466","0.00315","0.00319","0.00323","0.00327","0.00699","0.00699","0.01107","0.01173","0.02131","0.02214","0.01826","0.01937","0.01998"],["***","***","***","***","***","***","***","***","***","***","***","***","***","***","***","***","***","***","***","***","***","***","***","***","***","***","***","***","***","***","**","**","**","**","**","**","**","**","**","**","**","**","**","**","**","**","**","**","**","**","***","***","***","***","***","***","***","***","**","**","**","**","**","**","**","**","**","**","**","**","**","**","**","**","**","**","**","**","**","*","*","*","*","*","*","*"],["Wilcoxon","Wilcoxon","Wilcoxon","Wilcoxon","Wilcoxon","Wilcoxon","Wilcoxon","Wilcoxon","Wilcoxon","Wilcoxon","Wilcoxon","Wilcoxon","Wilcoxon","Wilcoxon","Wilcoxon","Wilcoxon","Wilcoxon","Wilcoxon","Wilcoxon","Wilcoxon","Wilcoxon","Wilcoxon","Wilcoxon","Wilcoxon","Wilcoxon","Wilcoxon","Wilcoxon","Wilcoxon","Wilcoxon","Wilcoxon","Wilcoxon","Wilcoxon","Wilcoxon","Wilcoxon","Wilcoxon","Wilcoxon","Wilcoxon","Wilcoxon","Wilcoxon","Wilcoxon","Wilcoxon","Wilcoxon","Wilcoxon","Wilcoxon","Wilcoxon","Wilcoxon","Wilcoxon","Wilcoxon","Wilcoxon","Wilcoxon","Wilcoxon","Wilcoxon","Wilcoxon","Wilcoxon","Wilcoxon","Wilcoxon","Wilcoxon","Wilcoxon","Wilcoxon","Wilcoxon","Wilcoxon","Wilcoxon","Wilcoxon","Wilcoxon","Wilcoxon","Wilcoxon","Wilcoxon","Wilcoxon","Wilcoxon","Wilcoxon","Wilcoxon","Wilcoxon","Wilcoxon","Wilcoxon","Wilcoxon","Wilcoxon","Wilcoxon","Wilcoxon","Wilcoxon","Wilcoxon","Wilcoxon","Wilcoxon","Wilcoxon","Wilcoxon","Wilcoxon","Wilcoxon"]],"container":"<table class=\"display\">\n  <thead>\n    <tr>\n      <th> <\/th>\n      <th>Model2<\/th>\n      <th>variable<\/th>\n      <th>.y.<\/th>\n      <th>group1<\/th>\n      <th>group2<\/th>\n      <th>p<\/th>\n      <th>p.adj<\/th>\n      <th>p.format<\/th>\n      <th>p.signif<\/th>\n      <th>method<\/th>\n    <\/tr>\n  <\/thead>\n<\/table>","options":{"columnDefs":[{"className":"dt-right","targets":[6,7]},{"orderable":false,"targets":0}],"order":[],"autoWidth":false,"orderClasses":false}},"evals":[],"jsHooks":[]}</script>
```


```r
# alpha_plot_2
# 
# 
# library(tidyverse)
# library(broom)
# 
# # original graph with smoother
# ggplot(data=mtcars, aes(hp,wt)) + 
#   stat_smooth(method = "loess", span = 0.75)
# 
# # Create model that will do the same thing as under the hood in ggplot2
# model <- loess(wt ~ hp, data = mtcars, span = 0.75)
# 
# # Add predicted values from model to original dataset using broom library
# mtcars2 <- augment(model, mtcars)
# 
# # Plot both lines 
# ggplot(data=mtcars2, aes(hp,wt)) + 
#   geom_line(aes(hp, .fitted), color = "red") +
#   stat_smooth(method = "loess", span = 0.75)
```


```r
# library(ggplot2)
# # Make the plot
# ggplot(aes(x = speed, y = dist), data = cars) + geom_point() +
#   stat_smooth(method = "loess")
# # Get the values
# smooth_vals = predict(loess(dist~speed,cars), cars$speed)
# 
# 
# plx<-predict(loess(cars$dist ~ cars$speed), se=T)
# 
# lines(cars$speed,plx$fit)
# lines(cars$speed,plx$fit - qt(0.975,plx$df)*plx$se, lty=2)
# lines(cars$speed,plx$fit + qt(0.975,plx$df)*plx$se, lty=2)
```


```r
# alpha_plot_2 %>% 
#   export::graph2ppt(file = "~/Downloads/16S-NRP72/pptx.pptx", append = TRUE,  width = 317.48031496 * 1.5, height = 0.618 * 317.48031496 * 1.075, paper = "A4",  scaling = 2)
```


```r
# ggsave("Alpha_chicken_2_Days_new.pdf", 
#        plot = alpha_plot_2, width = 10, height = 10, units = "cm", bg = "transparent", scale = 2, dpi = 600)
```


```r
# alpha_plot_2 %>% 
#   saveRDS("C:/Users/domga/Desktop/R stuff/Mapfiles/alpha_plot_2.RDS")
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
## [1] stats     graphics  grDevices utils     datasets  methods   base     
## 
## other attached packages:
##  [1] gdtools_0.2.4      reshape2_1.4.4     scales_1.2.1       compositions_2.0-4
##  [5] nlme_3.1-159       phyloseq_1.40.0    forcats_0.5.2      stringr_1.4.1     
##  [9] dplyr_1.0.10       purrr_0.3.5        readr_2.1.3        tidyr_1.2.1       
## [13] tibble_3.1.8       ggplot2_3.3.6      tidyverse_1.3.2   
## 
## loaded via a namespace (and not attached):
##   [1] readxl_1.4.1           uuid_1.1-0             backports_1.4.1       
##   [4] systemfonts_1.0.4      plyr_1.8.7             igraph_1.3.5          
##   [7] splines_4.2.1          crosstalk_1.2.0        GenomeInfoDb_1.32.4   
##  [10] digest_0.6.30          foreach_1.5.2          htmltools_0.5.3       
##  [13] fansi_1.0.3            magrittr_2.0.3         xlsx_0.6.5            
##  [16] googlesheets4_1.0.1    cluster_2.1.4          tzdb_0.3.0            
##  [19] openxlsx_4.2.5         Biostrings_2.64.1      extrafont_0.18        
##  [22] modelr_0.1.9           bayesm_3.1-4           officer_0.4.4         
##  [25] extrafontdb_1.0        colorspace_2.0-3       rvest_1.0.3           
##  [28] textshaping_0.3.6      haven_2.5.1            xfun_0.34             
##  [31] crayon_1.5.2           RCurl_1.98-1.9         jsonlite_1.8.3        
##  [34] survival_3.4-0         iterators_1.0.14       ape_5.6-2             
##  [37] glue_1.6.2             rvg_0.2.5              gtable_0.3.1          
##  [40] gargle_1.2.1           zlibbioc_1.42.0        XVector_0.36.0        
##  [43] car_3.1-0              Rttf2pt1_1.3.10        Rhdf5lib_1.18.2       
##  [46] BiocGenerics_0.42.0    DEoptimR_1.0-11        abind_1.4-5           
##  [49] DBI_1.1.3              rstatix_0.7.0          Rcpp_1.0.9            
##  [52] xtable_1.8-4           DT_0.25                stats4_4.2.1          
##  [55] htmlwidgets_1.5.4      httr_1.4.4             ellipsis_0.3.2        
##  [58] rJava_1.0-6            pkgconfig_2.0.3        farver_2.1.1          
##  [61] sass_0.4.2             dbplyr_2.2.1           utf8_1.2.2            
##  [64] here_1.0.1             tidyselect_1.2.0       labeling_0.4.2        
##  [67] rlang_1.0.6            munsell_0.5.0          cellranger_1.1.0      
##  [70] tools_4.2.1            cachem_1.0.6           cli_3.4.1             
##  [73] generics_0.1.3         devEMF_4.1             ade4_1.7-19           
##  [76] export_0.3.0           broom_1.0.1            evaluate_0.17         
##  [79] biomformat_1.24.0      fastmap_1.1.0          yaml_2.3.6            
##  [82] ragg_1.2.3             knitr_1.40             fs_1.5.2              
##  [85] zip_2.2.1              robustbase_0.95-0      rgl_0.110.2           
##  [88] xml2_1.3.3             compiler_4.2.1         rstudioapi_0.14       
##  [91] ggsignif_0.6.4         reprex_2.0.2           bslib_0.4.0           
##  [94] stringi_1.7.8          highr_0.9              stargazer_5.2.3       
##  [97] lattice_0.20-45        Matrix_1.5-1           vegan_2.6-4           
## [100] microbiome_1.19.1      permute_0.9-7          tensorA_0.36.2        
## [103] multtest_2.52.0        vctrs_0.4.2            pillar_1.8.1          
## [106] lifecycle_1.0.3        rhdf5filters_1.8.0     jquerylib_0.1.4       
## [109] data.table_1.14.4      cowplot_1.1.1          bitops_1.0-7          
## [112] flextable_0.8.2        R6_2.5.1               IRanges_2.30.1        
## [115] codetools_0.2-18       MASS_7.3-58.1          assertthat_0.2.1      
## [118] xlsxjars_0.6.1         rhdf5_2.40.0           rprojroot_2.0.3       
## [121] withr_2.5.0            S4Vectors_0.34.0       GenomeInfoDbData_1.2.8
## [124] mgcv_1.8-40            parallel_4.2.1         hms_1.1.2             
## [127] grid_4.2.1             rmarkdown_2.16         carData_3.0-5         
## [130] googledrive_2.0.0      Rtsne_0.16             ggpubr_0.4.0          
## [133] Biobase_2.56.0         lubridate_1.8.0        base64enc_0.1-3
```
