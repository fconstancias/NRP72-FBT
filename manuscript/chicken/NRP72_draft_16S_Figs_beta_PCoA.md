---
title: "NTP72 - 16S - beta - chicken draft "
author: "Florentin Constancias"
date: "November 01, 2022"
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
out_pptx = "~/Desktop/16S-chicken-beta_PCoA.pptx"
out_xlsx = "~/Desktop/16S-chicken-beta_PCoA.xlsx"
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
  subset_samples(Reactor != "IR") -> ps_filtered
```
a
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

## Beta diversits plots:

Ideal stuff: we want NMDS or T-SNE to show off here + sheppard plots.


```r
# depth = 8273  #or 18017*0.85 #9182 mIMT1_Scharl_alpha #8256 mIMT1_Borsig_alpha.html #9182 mIMT2_engraftment.html 9182 mIMT2_Cecum_alpha.html  9182 mIMT2_20200204_beta_trajectories.html
# 
# physeq %>%
#   rarefy_even_depth(sample.size = depth,
#                     rngseed = 123) -> physeq_rare
#
# physeq %>%
#   prune_samples(sample_sums(.)>= depth, .) -> physeq_fil
```


```r
ps_rare %>%
  # microbiome::transform(transform = "hellinger") %>% 
  phyloseq_compute_bdiv(phylo = FALSE) -> dist

dist$bray = NULL
dist$sorensen = NULL
dist$wjaccard = NULL
# dist$bjaccard = NULL

ps_fil %>%
  phyloseq_plot_bdiv(ps_rare = .,
                     m = "CoDa") -> ps_coda

dist$aichinson <- ps_coda$aidist %>% magrittr::divide_by(100)
# ps_fil  %>%
#   filter_taxa(function(x) sum(x > 0) > 0, TRUE) %>%
#   microbiome::transform(transform = "compositional") %>%
#   speedyseq::psmelt() %>%
#   pivot_wider(names_from =  Sample, OTU,
#               values_from = Abundance) %>%
#   # mutate(PathoFact_AMR_ARG = PathoFact_AMR_ARG %>%  str_to_lower(.)) %>%
#   # %>% str_to_lower(., locale = "en")
#   column_to_rownames("OTU") %>%
#   t() %>%
#   replace(is.na(.), 0)  %>%
#   vegan::vegdist(na.rm = TRUE, method = "robust.aitchison") %>% 
#   magrittr::divide_by(100) -> dist$robust.aitchison
#   # cor() %>%
#   # hclust(method = "ward.D2") -> clust
#
```



```r
# plot(dist$wjaccard, 
#      dist$aichinson)
# 
plot(dist$bjaccard,
     dist$aichinson)
```

![](NRP72_draft_16S_Figs_beta_PCoA_files/figure-html/unnamed-chunk-7-1.png)<!-- -->

### F1 & F2:


```r
ps_rare %>%
  # microbiome::transform(transform = "compositional") %>% 
  phyloseq_plot_bdiv(dlist = dist,
                     m = "PCoA",
                     seed = 123,
                     axis1 = 1,
                     axis2 = 2) -> plots_hall_humans
```

```
## [1] "bjaccard"
## [1] "aichinson"
```

```r
phyloseq_plot_ordinations_facet(plot_list = plots_hall_humans,
                                shape_group = "Antibiotic_mg.L",
                                color_group = "Treatment",
                                alpha = NULL
) +
  scale_color_manual(name = "", values = treat_col,
                     na.value = "black") +
  scale_fill_manual(name = "", values = treat_col,
                    na.value = "black") +
  scale_shape_manual(name = "" ,values = c(15,16,18,19), na.value =  17) +
  theme(legend.position = "bottom") +
  facet_wrap(distance ~ ., scales = "free")
```

![](NRP72_draft_16S_Figs_beta_PCoA_files/figure-html/unnamed-chunk-8-1.png)<!-- -->


```r
# plots_hall_humans$bjaccard = NULL
plots_hall_humans$wjaccard = NULL

phyloseq_plot_ordinations_facet(plot_list = plots_hall_humans,
                                shape_group = "Antibiotic_mg.L",
                                color_group = "Treatment",
                                alpha = NULL
) +
  geom_path(data = rbind(plots_hall_humans$bjaccard$data %>% mutate(distance = "bjaccard"), plots_hall_humans$aichinson$data %>% mutate(distance = "aichinson")) %>%
              arrange(Day_of_Treatment) ,
            # aes(colour = Treatment, group = interaction(Model, Model2, Antibiotic, Treatment, Fermentation, Reactor,Antibiotic_mg.L)),
            aes(colour = Treatment, group = interaction(Model2, Reactor, Treatment, Treatment_Dose)),
            
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
  theme(legend.position = "right") +
  ggtitle(NULL) +
  facet_grid(distance ~ ., scales = "free") -> pcoa_all

pcoa_all + theme(legend.position = "none") -> pcoa_all_noleg

pcoa_all  %>%  ggpubr::get_legend() %>% ggpubr::as_ggplot() -> pcoa_all_leg

pcoa_all_noleg
```

![](NRP72_draft_16S_Figs_beta_PCoA_files/figure-html/unnamed-chunk-9-1.png)<!-- -->



```r
pcoa_all_noleg
```

![](NRP72_draft_16S_Figs_beta_PCoA_files/figure-html/unnamed-chunk-10-1.png)<!-- -->


```r
pcoa_all_noleg %>%
  export::graph2ppt(append = TRUE,
                    width = 317.48031496 * 1,
                    height = 0.618 * 317.48031496 * 1 , paper = "A4",  scaling = 2,
                    file = out_pptx)
```

```
## Exported graph as ~/Desktop/16S-chicken-beta_PCoA.pptx
```

```r
pcoa_all_leg %>%
  export::graph2ppt(append = TRUE,
                    width = 317.48031496 * 1,
                    height = 0.618 * 317.48031496 * 1 , paper = "A4",  scaling = 2,
                    file = out_pptx)
```

```
## Exported graph as ~/Desktop/16S-chicken-beta_PCoA.pptx
```



```r
plots_hall_humans %>% 
  phyloseq_ordinations_expl_var(.) -> expl_var
```

```
## New names:
## • `.id` -> `.id...1`
## • `V1` -> `V1...2`
## • `.id` -> `.id...3`
## • `V1` -> `V1...4`
```

```r
expl_var %>% 
  DT::datatable()
```

```{=html}
<div id="htmlwidget-9eceb20dc2a931752ae4" style="width:100%;height:auto;" class="datatables html-widget"></div>
<script type="application/json" data-for="htmlwidget-9eceb20dc2a931752ae4">{"x":{"filter":"none","vertical":false,"data":[["1","2"],["bjaccard","aichinson"],["Axis.1   [25.7%]","Axis.1   [27.1%]"],["bjaccard","aichinson"],["Axis.2   [20.6%]","Axis.2   [19.6%]"]],"container":"<table class=\"display\">\n  <thead>\n    <tr>\n      <th> <\/th>\n      <th>.id...1<\/th>\n      <th>V1...2<\/th>\n      <th>.id...3<\/th>\n      <th>V1...4<\/th>\n    <\/tr>\n  <\/thead>\n<\/table>","options":{"columnDefs":[{"orderable":false,"targets":0}],"order":[],"autoWidth":false,"orderClasses":false}},"evals":[],"jsHooks":[]}</script>
```

```r
expl_var %>%
  data.frame() %>% 
  export::table2ppt(append = TRUE,
                    file = out_pptx)
```

```
## Exported table as ~/Desktop/16S-chicken-beta_PCoA.pptx
```

```{=html}
<template id="43cb2120-5936-4f06-a678-3958ca2b40f8"><style>
.tabwid table{
  border-spacing:0px !important;
  border-collapse:collapse;
  line-height:1;
  margin-left:auto;
  margin-right:auto;
  border-width: 0;
  display: table;
  border-color: transparent;
  caption-side: top;
}
.tabwid-caption-bottom table{
  caption-side: bottom;
}
.tabwid_left table{
  margin-left:0;
}
.tabwid_right table{
  margin-right:0;
}
.tabwid td {
    padding: 0;
}
.tabwid a {
  text-decoration: none;
}
.tabwid thead {
    background-color: transparent;
}
.tabwid tfoot {
    background-color: transparent;
}
.tabwid table tr {
background-color: transparent;
}
.katex-display {
    margin: 0 0 !important;
}
</style><div class="tabwid"><style>.cl-a1a26940{}.cl-a19ccf6c{font-family:'Helvetica';font-size:12pt;font-weight:bold;font-style:normal;text-decoration:none;color:rgba(0, 0, 0, 1.00);background-color:transparent;}.cl-a19ccf80{font-family:'Helvetica';font-size:12pt;font-weight:normal;font-style:normal;text-decoration:none;color:rgba(0, 0, 0, 1.00);background-color:transparent;}.cl-a19f4df0{margin:0;text-align:left;border-bottom: 0 solid rgba(0, 0, 0, 1.00);border-top: 0 solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);padding-bottom:5pt;padding-top:5pt;padding-left:5pt;padding-right:5pt;line-height: 1;background-color:transparent;}.cl-a19f5a5c{width:2.14in;background-color:transparent;vertical-align: middle;border-bottom: 1pt solid rgba(0, 0, 0, 1.00);border-top: 1pt solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-a19f5a66{width:2.14in;background-color:transparent;vertical-align: middle;border-bottom: 0 solid rgba(0, 0, 0, 1.00);border-top: 0 solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-a19f5a70{width:2.14in;background-color:transparent;vertical-align: middle;border-bottom: 1pt solid rgba(0, 0, 0, 1.00);border-top: 0 solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}</style><table class='cl-a1a26940'><thead><tr style="overflow-wrap:break-word;"><td class="cl-a19f5a5c"><p class="cl-a19f4df0"><span class="cl-a19ccf6c">.id...1</span></p></td><td class="cl-a19f5a5c"><p class="cl-a19f4df0"><span class="cl-a19ccf6c">V1...2</span></p></td><td class="cl-a19f5a5c"><p class="cl-a19f4df0"><span class="cl-a19ccf6c">.id...3</span></p></td><td class="cl-a19f5a5c"><p class="cl-a19f4df0"><span class="cl-a19ccf6c">V1...4</span></p></td></tr></thead><tbody><tr style="overflow-wrap:break-word;"><td class="cl-a19f5a66"><p class="cl-a19f4df0"><span class="cl-a19ccf80">bjaccard</span></p></td><td class="cl-a19f5a66"><p class="cl-a19f4df0"><span class="cl-a19ccf80">Axis.1   [25.7%]</span></p></td><td class="cl-a19f5a66"><p class="cl-a19f4df0"><span class="cl-a19ccf80">bjaccard</span></p></td><td class="cl-a19f5a66"><p class="cl-a19f4df0"><span class="cl-a19ccf80">Axis.2   [20.6%]</span></p></td></tr><tr style="overflow-wrap:break-word;"><td class="cl-a19f5a70"><p class="cl-a19f4df0"><span class="cl-a19ccf80">aichinson</span></p></td><td class="cl-a19f5a70"><p class="cl-a19f4df0"><span class="cl-a19ccf80">Axis.1   [27.1%]</span></p></td><td class="cl-a19f5a70"><p class="cl-a19f4df0"><span class="cl-a19ccf80">aichinson</span></p></td><td class="cl-a19f5a70"><p class="cl-a19f4df0"><span class="cl-a19ccf80">Axis.2   [19.6%]</span></p></td></tr></tbody></table></div></template>
<div class="flextable-shadow-host" id="f22a9731-2a81-4e82-97c0-d2aaf2b26528"></div>
<script>
var dest = document.getElementById("f22a9731-2a81-4e82-97c0-d2aaf2b26528");
var template = document.getElementById("43cb2120-5936-4f06-a678-3958ca2b40f8");
var caption = template.content.querySelector("caption");
var fantome = dest.attachShadow({mode: 'open'});
var templateContent = template.content;
fantome.appendChild(templateContent);
</script>

```


Stats evaluation:

Treatment vs host impact?


```r
lapply(
  dist,
  FUN = phyloseq_adonis2,
  physeq = ps_rare %>%
    subset_samples(Treatment != "DONOR" & Day_of_Treatment >= 1),
  formula = paste0(c("Model2", "Treatment"), collapse=" + "),
  nrep = 999,
  strata = "none",
  terms_margins = "margin"
) %>%
  bind_rows(.id = "Distance") %>%
  mutate_if(is.numeric, round, 3) %>%
  DT::datatable()
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
## This is vegan 2.6-4
```

```{=html}
<div id="htmlwidget-f30d7598595affc4f25a" style="width:100%;height:auto;" class="datatables html-widget"></div>
<script type="application/json" data-for="htmlwidget-f30d7598595affc4f25a">{"x":{"filter":"none","vertical":false,"data":[["1","2","3","4","5","6","7","8"],["bjaccard","bjaccard","bjaccard","bjaccard","aichinson","aichinson","aichinson","aichinson"],["Model2","Treatment","Residual","Total","Model2","Treatment","Residual","Total"],[1,6,164,171,1,6,164,171],[8.538,15.709,24.526,49.637,4.155,6.697,10.081,21.51],[0.172,0.316,0.494,1,0.193,0.311,0.469,1],[57.091,17.507,null,null,67.596,18.16,null,null],[0.001,0.001,null,null,0.001,0.001,null,null]],"container":"<table class=\"display\">\n  <thead>\n    <tr>\n      <th> <\/th>\n      <th>Distance<\/th>\n      <th>terms<\/th>\n      <th>Df<\/th>\n      <th>SumOfSqs<\/th>\n      <th>R2<\/th>\n      <th>F<\/th>\n      <th>Pr(&gt;F)<\/th>\n    <\/tr>\n  <\/thead>\n<\/table>","options":{"columnDefs":[{"className":"dt-right","targets":[3,4,5,6,7]},{"orderable":false,"targets":0}],"order":[],"autoWidth":false,"orderClasses":false}},"evals":[],"jsHooks":[]}</script>
```



### Chicken1:


```r
ps_rare %>%
  subset_samples(Model2 == "Chicken1" & Treatment != "DONOR") -> ps_tmp_forbdiv
```

#### PCoA:


```r
ps_tmp_forbdiv %>%
  # microbiome::transform(transform = "compositional") %>% 
  phyloseq_plot_bdiv(dlist = dist,
                     m = "PCoA",
                     seed = 123,
                     axis1 = 1,
                     axis2 = 2) -> PCoA_tmp
```

```
## [1] "bjaccard"
## [1] "aichinson"
```

#### aichinson


```r
PCoA_tmp$aichinson$layers[[1]] = NULL; PCoA_tmp$aichinson$layers[[1]] = NULL
PCoA_tmp$aichinson$layers[[2]] = NULL; PCoA_tmp$aichinson$layers[[1]] = NULL

# plots_hall_humans$aichinson$layers[[1]] = NULL;plots_hall_humans$aichinson$layers[[1]] = NULL
# plots_hall_humans$aichinson$layers[[2]] = NULL;plots_hall_humans$aichinson$layers[[2]] = NULL


PCoA_tmp$aichinson + geom_point(size = 3,
                                aes(colour = Treatment ,
                                    shape = Antibiotic_mg.L,
                                    alpha = NULL)) +
  geom_path(data = PCoA_tmp$aichinson$data %>%
              arrange(Day_of_Treatment) ,
            # aes(colour = Treatment, group = interaction(Model, Model2, Antibiotic, Treatment, Fermentation, Reactor,Antibiotic_mg.L)),
            aes(colour = Treatment, group = interaction(Model2, Reactor, Treatment, Treatment_Dose)),
            
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
  theme(legend.position = "right") +
  facet_grid(. ~ Model2, scales = "free") -> PCoA_tmp2

PCoA_tmp2
```

![](NRP72_draft_16S_Figs_beta_PCoA_files/figure-html/unnamed-chunk-16-1.png)<!-- -->


```r
PCoA_tmp2 + 
  scale_fill_manual(values = c("transparent")) + 
  scale_color_manual(values = c(rep("transparent", 7))) + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) -> empty_plot_tmp


ps_tmp_forbdiv %>%
  phyloseq_add_taxa_vector_fix(phyloseq = .,
                               dist = dist$aichinson,
                               taxrank_glom = "Family",
                               figure_ord = empty_plot_tmp, 
                               fact = 0.2, pval_cutoff = 0.05,
                               top_r = 8) -> out #top 10 correlated genera

out$signenvfit %>% 
  DT::datatable()
```

```{=html}
<div id="htmlwidget-b85bf888dc62471906a9" style="width:100%;height:auto;" class="datatables html-widget"></div>
<script type="application/json" data-for="htmlwidget-b85bf888dc62471906a9">{"x":{"filter":"none","vertical":false,"data":[["ASV0001","ASV0008","ASV0006","ASV0134","ASV0125","ASV0165","ASV0187","ASV0278"],[0.179126285157017,-0.417317172405722,-0.520111335186489,-0.469901977442166,-0.568941471660472,0.549879695939574,-0.617503745304176,-0.620583212859231],[0.760825246934185,-0.768280896635797,-0.722683637195785,-0.737611838158516,0.420534229374111,0.462664680607074,0.42291115640916,0.417913549423581],[0.001,0.001,0.001,0.001,0.001,0.001,0.001,0.001],[0.610941282406617,0.764409158520211,0.792787440460002,0.764879092195644,0.500543436250461,0.516426286687844,0.560164721680017,0.559775258874501],["Bacteria","Bacteria","Bacteria","Bacteria","Bacteria","Bacteria","Bacteria","Bacteria"],["Proteobacteria","Firmicutes","Firmicutes","Firmicutes","Firmicutes","Firmicutes","Firmicutes","Firmicutes"],["Gammaproteobacteria","Clostridia","Clostridia","Clostridia","Bacilli","Clostridia","Clostridia","Clostridia"],["Enterobacterales","Oscillospirales","Lachnospirales","Oscillospirales","Erysipelotrichales","Peptostreptococcales-Tissierellales","Christensenellales","Monoglobales"],["Morganellaceae","Ruminococcaceae","Lachnospiraceae","Oscillospiraceae","Erysipelatoclostridiaceae","Anaerovoracaceae","Christensenellaceae","Monoglobaceae"],[null,null,null,null,null,null,null,null],[null,null,null,null,null,null,null,null],[null,null,null,null,null,null,null,null],[0.00263157894736842,0.00263157894736842,0.00263157894736842,0.00263157894736842,0.00263157894736842,0.00263157894736842,0.00263157894736842,0.00263157894736842]],"container":"<table class=\"display\">\n  <thead>\n    <tr>\n      <th> <\/th>\n      <th>Axis.1<\/th>\n      <th>Axis.2<\/th>\n      <th>pval<\/th>\n      <th>r<\/th>\n      <th>Kingdom<\/th>\n      <th>Phylum<\/th>\n      <th>Class<\/th>\n      <th>Order<\/th>\n      <th>tax_rank_plot<\/th>\n      <th>Genus<\/th>\n      <th>Species<\/th>\n      <th>Strain<\/th>\n      <th>pval.adj<\/th>\n    <\/tr>\n  <\/thead>\n<\/table>","options":{"columnDefs":[{"className":"dt-right","targets":[1,2,3,4,13]},{"orderable":false,"targets":0}],"order":[],"autoWidth":false,"orderClasses":false}},"evals":[],"jsHooks":[]}</script>
```


```r
out$plot -> PCoA_tmp2_fam

PCoA_tmp2 + 
  theme(legend.position = "none") -> PCoA_tmp2_noleg

PCoA_tmp2_noleg
```

![](NRP72_draft_16S_Figs_beta_PCoA_files/figure-html/unnamed-chunk-18-1.png)<!-- -->

```r
PCoA_tmp2  %>%  ggpubr::get_legend() %>% ggpubr::as_ggplot() -> PCoA_tmp2_leg

PCoA_tmp2_noleg %>% 
  export::graph2ppt(append = TRUE,
                    width = 317.48031496 * 1,
                    height = 0.618 * 317.48031496 * 1 , paper = "A4",  scaling = 2,
                    file = out_pptx)

PCoA_tmp2_leg %>% 
  export::graph2ppt(append = TRUE,
                    width = 317.48031496 * 1,
                    height = 0.618 * 317.48031496 * 1 , paper = "A4",  scaling = 2,
                    file = out_pptx)

PCoA_tmp2_fam %>% 
  export::graph2ppt(append = TRUE,
                    width = 317.48031496 * 1,
                    height = 0.618 * 317.48031496 * 1 , paper = "A4",  scaling = 2,
                    file = out_pptx)
```


```r
ps_tmp_forbdiv %>% 
  microViz::ps_mutate(Lactose_mM = as.double(Lactose_mM), 
                      Glucose_mM = as.double(Glucose_mM),
                      Galactose_mM = as.double(Galactose_mM),
                      Succinat_mM = as.double(Succinat_mM),
                      Lactat_mM = as.double(Lactat_mM),
                      Formiat_mM = as.double(Formiat_mM),
                      Acetat_mM = as.double(Acetat_mM),
                      Propionat_mM = as.double(Propionat_mM),
                      Butyrat_mM = as.double(Butyrat_mM),
                      Valerat_mM = as.double(Valerat_mM),
                      Isobutyrat_mM = as.double(Isobutyrat_mM),
                      Isovalerat_mM = as.double(Isovalerat_mM),
                      Total_SCFA_mM = as.double(Total_SCFA_mM))  %>% # sample_data() %>% data.frame()
  subset_samples(Total_SCFA_mM > 0) %>% 
  phyloseq_add_metadata_vector(dist = dist$aichinson,
                               phyloseq = .,
                               figure_ord = empty_plot_tmp,
                               metadata_sel = c("Succinat_mM", "Formiat_mM",  "Acetat_mM" ,
                                                "Propionat_mM", "Butyrat_mM",  "Valerat_mM", "Total_SCFA_mM"),
                               top_r = 6, fact = 0.3, na.rm = TRUE,
                               pval_cutoff = 0.05) -> env_out


env_out$signenvfit %>% 
  DT::datatable()
```

```{=html}
<div id="htmlwidget-88fdd2a34d502a66744b" style="width:100%;height:auto;" class="datatables html-widget"></div>
<script type="application/json" data-for="htmlwidget-88fdd2a34d502a66744b">{"x":{"filter":"none","vertical":false,"data":[["1","2","3","4","5","6"],["Succinat_mM","Acetat_mM","Propionat_mM","Butyrat_mM","Valerat_mM","Total_SCFA_mM"],[0.7987843132071,-0.763822131606401,0.27159917848853,-0.61519382446119,-0.639755304442811,-0.811340602560074],[0.178331284674618,-0.400906971626607,0.683522424509029,-0.709613437960205,0.412493789226407,-0.389815834112515],["Succinat_mM","Acetat_mM","Propionat_mM","Butyrat_mM","Valerat_mM","Total_SCFA_mM"],[0.001,0.001,0.001,0.001,0.001,0.001],[0.669858426119437,0.744150648630563,0.540969018562345,0.882014672988888,0.579437975713073,0.810229957887379]],"container":"<table class=\"display\">\n  <thead>\n    <tr>\n      <th> <\/th>\n      <th>label<\/th>\n      <th>Axis.1<\/th>\n      <th>Axis.2<\/th>\n      <th>id<\/th>\n      <th>pval<\/th>\n      <th>r<\/th>\n    <\/tr>\n  <\/thead>\n<\/table>","options":{"columnDefs":[{"className":"dt-right","targets":[2,3,5,6]},{"orderable":false,"targets":0}],"order":[],"autoWidth":false,"orderClasses":false}},"evals":[],"jsHooks":[]}</script>
```


```r
env_out$plot %>% 
  export::graph2ppt(append = TRUE,
                    width = 317.48031496 * 1,
                    height = 0.618 * 317.48031496 * 1 , paper = "A4",  scaling = 2,
                    file = out_pptx)
```

#### bjaccard


```r
PCoA_tmp$bjaccard$layers[[1]] = NULL; PCoA_tmp$bjaccard$layers[[1]] = NULL
PCoA_tmp$bjaccard$layers[[2]] = NULL; PCoA_tmp$bjaccard$layers[[1]] = NULL

# plots_hall_humans$aichinson$layers[[1]] = NULL;plots_hall_humans$aichinson$layers[[1]] = NULL
# plots_hall_humans$aichinson$layers[[2]] = NULL;plots_hall_humans$aichinson$layers[[2]] = NULL


PCoA_tmp$bjaccard + geom_point(size = 3,
                               aes(colour = Treatment ,
                                   shape = Antibiotic_mg.L,
                                   alpha = NULL)) +
  geom_path(data = PCoA_tmp$bjaccard$data %>%
              arrange(Day_of_Treatment) ,
            # aes(colour = Treatment, group = interaction(Model, Model2, Antibiotic, Treatment, Fermentation, Reactor,Antibiotic_mg.L)),
            aes(colour = Treatment, group = interaction(Model2, Reactor, Treatment, Treatment_Dose)),
            
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
  theme(legend.position = "right") +
  facet_grid(. ~ Model2, scales = "free") -> PCoA_tmp2

PCoA_tmp2
```

![](NRP72_draft_16S_Figs_beta_PCoA_files/figure-html/unnamed-chunk-21-1.png)<!-- -->


```r
PCoA_tmp2 + 
  scale_fill_manual(values = c("transparent")) + 
  scale_color_manual(values = c(rep("transparent", 7))) + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) -> empty_plot_tmp


ps_tmp_forbdiv %>%
  phyloseq_add_taxa_vector_fix(phyloseq = .,
                               dist = dist$bjaccard,
                               taxrank_glom = "Family",
                               figure_ord = empty_plot_tmp, 
                               fact = 0.2, pval_cutoff = 0.05,
                               top_r = 8) -> out #top 10 correlated genera

out$signenvfit %>% 
  DT::datatable()
```

```{=html}
<div id="htmlwidget-8e1fe96c0688a849169a" style="width:100%;height:auto;" class="datatables html-widget"></div>
<script type="application/json" data-for="htmlwidget-8e1fe96c0688a849169a">{"x":{"filter":"none","vertical":false,"data":[["ASV0001","ASV0008","ASV0006","ASV0134","ASV0109","ASV0165","ASV0187","ASV0278"],[0.231618167499601,-0.463549587310127,-0.576476302781683,-0.525763548777435,-0.344570047491009,0.649893612073232,-0.473121554133513,-0.427436772742712],[-0.715270448514504,0.710032415965814,0.671084322176517,0.649010202443495,-0.589108997800361,-0.375615881402725,-0.610050545663667,-0.606103099393499],[0.001,0.001,0.001,0.001,0.001,0.001,0.001,0.001],[0.565258790034014,0.71902425161764,0.782679095139955,0.697641552098789,0.465777928917301,0.563448997375539,0.596005673250249,0.55006316178711],["Bacteria","Bacteria","Bacteria","Bacteria","Bacteria","Bacteria","Bacteria","Bacteria"],["Proteobacteria","Firmicutes","Firmicutes","Firmicutes","Firmicutes","Firmicutes","Firmicutes","Firmicutes"],["Gammaproteobacteria","Clostridia","Clostridia","Clostridia","Clostridia","Clostridia","Clostridia","Clostridia"],["Enterobacterales","Oscillospirales","Lachnospirales","Oscillospirales","Lachnospirales","Peptostreptococcales-Tissierellales","Christensenellales","Monoglobales"],["Morganellaceae","Ruminococcaceae","Lachnospiraceae","Oscillospiraceae","Defluviitaleaceae","Anaerovoracaceae","Christensenellaceae","Monoglobaceae"],[null,null,null,null,null,null,null,null],[null,null,null,null,null,null,null,null],[null,null,null,null,null,null,null,null],[0.00263157894736842,0.00263157894736842,0.00263157894736842,0.00263157894736842,0.00263157894736842,0.00263157894736842,0.00263157894736842,0.00263157894736842]],"container":"<table class=\"display\">\n  <thead>\n    <tr>\n      <th> <\/th>\n      <th>Axis.1<\/th>\n      <th>Axis.2<\/th>\n      <th>pval<\/th>\n      <th>r<\/th>\n      <th>Kingdom<\/th>\n      <th>Phylum<\/th>\n      <th>Class<\/th>\n      <th>Order<\/th>\n      <th>tax_rank_plot<\/th>\n      <th>Genus<\/th>\n      <th>Species<\/th>\n      <th>Strain<\/th>\n      <th>pval.adj<\/th>\n    <\/tr>\n  <\/thead>\n<\/table>","options":{"columnDefs":[{"className":"dt-right","targets":[1,2,3,4,13]},{"orderable":false,"targets":0}],"order":[],"autoWidth":false,"orderClasses":false}},"evals":[],"jsHooks":[]}</script>
```


```r
out$plot -> PCoA_tmp2_fam

PCoA_tmp2 + 
  theme(legend.position = "none") -> PCoA_tmp2_noleg

PCoA_tmp2_noleg
```

![](NRP72_draft_16S_Figs_beta_PCoA_files/figure-html/unnamed-chunk-23-1.png)<!-- -->

```r
PCoA_tmp2  %>%  ggpubr::get_legend() %>% ggpubr::as_ggplot() -> PCoA_tmp2_leg

PCoA_tmp2_noleg %>% 
  export::graph2ppt(append = TRUE,
                    width = 317.48031496 * 1,
                    height = 0.618 * 317.48031496 * 1 , paper = "A4",  scaling = 2,
                    file = out_pptx)

PCoA_tmp2_leg %>% 
  export::graph2ppt(append = TRUE,
                    width = 317.48031496 * 1,
                    height = 0.618 * 317.48031496 * 1 , paper = "A4",  scaling = 2,
                    file = out_pptx)

PCoA_tmp2_fam %>% 
  export::graph2ppt(append = TRUE,
                    width = 317.48031496 * 1,
                    height = 0.618 * 317.48031496 * 1 , paper = "A4",  scaling = 2,
                    file = out_pptx)
```


```r
ps_tmp_forbdiv %>% 
  microViz::ps_mutate(Lactose_mM = as.double(Lactose_mM), 
                      Glucose_mM = as.double(Glucose_mM),
                      Galactose_mM = as.double(Galactose_mM),
                      Succinat_mM = as.double(Succinat_mM),
                      Lactat_mM = as.double(Lactat_mM),
                      Formiat_mM = as.double(Formiat_mM),
                      Acetat_mM = as.double(Acetat_mM),
                      Propionat_mM = as.double(Propionat_mM),
                      Butyrat_mM = as.double(Butyrat_mM),
                      Valerat_mM = as.double(Valerat_mM),
                      Isobutyrat_mM = as.double(Isobutyrat_mM),
                      Isovalerat_mM = as.double(Isovalerat_mM),
                      Total_SCFA_mM = as.double(Total_SCFA_mM))  %>% # sample_data() %>% data.frame()
  subset_samples(Total_SCFA_mM > 0) %>% 
  phyloseq_add_metadata_vector(dist = dist$bjaccard,
                               phyloseq = .,
                               figure_ord = empty_plot_tmp,
                               metadata_sel = c("Succinat_mM", "Formiat_mM",  "Acetat_mM" ,
                                                "Propionat_mM", "Butyrat_mM",  "Valerat_mM", "Total_SCFA_mM"),
                               top_r = 6, fact = 0.3, na.rm = TRUE,
                               pval_cutoff = 0.05) -> env_out


env_out$signenvfit %>% 
  DT::datatable()
```

```{=html}
<div id="htmlwidget-bb3486349adfa5d821dc" style="width:100%;height:auto;" class="datatables html-widget"></div>
<script type="application/json" data-for="htmlwidget-bb3486349adfa5d821dc">{"x":{"filter":"none","vertical":false,"data":[["1","2","3","4","5","6"],["Succinat_mM","Acetat_mM","Propionat_mM","Butyrat_mM","Valerat_mM","Total_SCFA_mM"],[0.805376707136237,-0.862698523914289,0.324337009107315,-0.717381936893163,-0.465080250779239,-0.90088101316165],[0.13593975442181,0.201311834400704,-0.660808679152135,0.590096643093981,-0.607635844076233,0.135182677818885],["Succinat_mM","Acetat_mM","Propionat_mM","Butyrat_mM","Valerat_mM","Total_SCFA_mM"],[0.001,0.001,0.001,0.001,0.001,0.001],[0.66711125722987,0.784775197833669,0.541862605919468,0.862850891571372,0.585520958671116,0.829860956257444]],"container":"<table class=\"display\">\n  <thead>\n    <tr>\n      <th> <\/th>\n      <th>label<\/th>\n      <th>Axis.1<\/th>\n      <th>Axis.2<\/th>\n      <th>id<\/th>\n      <th>pval<\/th>\n      <th>r<\/th>\n    <\/tr>\n  <\/thead>\n<\/table>","options":{"columnDefs":[{"className":"dt-right","targets":[2,3,5,6]},{"orderable":false,"targets":0}],"order":[],"autoWidth":false,"orderClasses":false}},"evals":[],"jsHooks":[]}</script>
```


```r
env_out$plot %>% 
  export::graph2ppt(append = TRUE,
                    width = 317.48031496 * 1,
                    height = 0.618 * 317.48031496 * 1 , paper = "A4",  scaling = 2,
                    file = out_pptx)
```

Statistical assessment:

- pret - Day_of_Treatment <= 0:


```r
lapply(
  dist,
  FUN = phyloseq_adonis2,
  physeq = ps_tmp_forbdiv %>%
    subset_samples(Day_of_Treatment <= 0),
  formula = paste0(c("Treatment"), collapse=" + "),
  nrep = 999,
  strata = "none",
  terms_margins = "margin"
) %>%
  bind_rows(.id = "Distance")  %>% 
  mutate_if(is.numeric, round, 3) -> tmp_df


tmp_df %>% 
  DT::datatable()
```

```{=html}
<div id="htmlwidget-44474b3e8cd160f85a98" style="width:100%;height:auto;" class="datatables html-widget"></div>
<script type="application/json" data-for="htmlwidget-44474b3e8cd160f85a98">{"x":{"filter":"none","vertical":false,"data":[["1","2","3","4","5","6"],["bjaccard","bjaccard","bjaccard","aichinson","aichinson","aichinson"],["Treatment","Residual","Total","Treatment","Residual","Total"],[5,13,18,5,13,18],[0.694,0.446,1.14,0.471,0.306,0.776],[0.609,0.391,1,0.606,0.394,1],[4.045,null,null,4.004,null,null],[0.001,null,null,0.001,null,null]],"container":"<table class=\"display\">\n  <thead>\n    <tr>\n      <th> <\/th>\n      <th>Distance<\/th>\n      <th>terms<\/th>\n      <th>Df<\/th>\n      <th>SumOfSqs<\/th>\n      <th>R2<\/th>\n      <th>F<\/th>\n      <th>Pr(&gt;F)<\/th>\n    <\/tr>\n  <\/thead>\n<\/table>","options":{"columnDefs":[{"className":"dt-right","targets":[3,4,5,6,7]},{"orderable":false,"targets":0}],"order":[],"autoWidth":false,"orderClasses":false}},"evals":[],"jsHooks":[]}</script>
```

```r
tmp_df %>% 
  xlsx::write.xlsx(file = out_xlsx,
                   append = TRUE,
                   sheetName = "chicken1-beta-pretreat",
                   showNA = TRUE)
```


```r
lapply(
  dist,
  FUN = phyloseq_adonis2,
  physeq = ps_tmp_forbdiv %>%
    subset_samples(Day_of_Treatment <= 0),
  formula = paste0(c("Treatment", "Day_of_Treatment"), collapse=" + "),
  nrep = 999,
  strata = "none",
  terms_margins = "margin"
) %>%
  bind_rows(.id = "Distance")  %>% 
  mutate_if(is.numeric, round, 3) -> tmp_df


tmp_df %>% 
  DT::datatable()
```

```{=html}
<div id="htmlwidget-2514ce28ae593bde0f73" style="width:100%;height:auto;" class="datatables html-widget"></div>
<script type="application/json" data-for="htmlwidget-2514ce28ae593bde0f73">{"x":{"filter":"none","vertical":false,"data":[["1","2","3","4","5","6","7","8"],["bjaccard","bjaccard","bjaccard","bjaccard","aichinson","aichinson","aichinson","aichinson"],["Treatment","Day_of_Treatment","Residual","Total","Treatment","Day_of_Treatment","Residual","Total"],[5,1,12,18,5,1,12,18],[0.71,0.083,0.364,1.14,0.474,0.046,0.259,0.776],[0.623,0.072,0.319,1,0.61,0.06,0.334,1],[4.69,2.725,null,null,4.382,2.14,null,null],[0.001,0.004,null,null,0.001,0.02,null,null]],"container":"<table class=\"display\">\n  <thead>\n    <tr>\n      <th> <\/th>\n      <th>Distance<\/th>\n      <th>terms<\/th>\n      <th>Df<\/th>\n      <th>SumOfSqs<\/th>\n      <th>R2<\/th>\n      <th>F<\/th>\n      <th>Pr(&gt;F)<\/th>\n    <\/tr>\n  <\/thead>\n<\/table>","options":{"columnDefs":[{"className":"dt-right","targets":[3,4,5,6,7]},{"orderable":false,"targets":0}],"order":[],"autoWidth":false,"orderClasses":false}},"evals":[],"jsHooks":[]}</script>
```

```r
tmp_df %>% 
  xlsx::write.xlsx(file = out_xlsx,
                   append = TRUE,
                   sheetName = "chicken1-beta-pretreat-time-treat",
                   showNA = TRUE)
```


```r
lapply(
  dist,
  FUN = physeq_pairwise_permanovas_adonis2,
  physeq = ps_tmp_forbdiv %>%
    subset_samples(Day_of_Treatment <= 0),
  compare_header = "Treatment",
  n_perm = 999,
  strat = "none"
) %>%
  bind_rows(.id = "Distance") %>%
  mutate_if(is.numeric, round, 3) %>%
  select(-pvalBon)  -> tmp_df
```

```
## 'nperm' >= set of all permutations: complete enumeration.
```

```
## Set of permutations < 'minperm'. Generating entire set.
## Set of permutations < 'minperm'. Generating entire set.
## Set of permutations < 'minperm'. Generating entire set.
## Set of permutations < 'minperm'. Generating entire set.
```

```
## 'nperm' >= set of all permutations: complete enumeration.
```

```
## Set of permutations < 'minperm'. Generating entire set.
```

```
## 'nperm' >= set of all permutations: complete enumeration.
```

```
## Set of permutations < 'minperm'. Generating entire set.
```

```
## 'nperm' >= set of all permutations: complete enumeration.
```

```
## Set of permutations < 'minperm'. Generating entire set.
```

```
## 'nperm' >= set of all permutations: complete enumeration.
```

```
## Set of permutations < 'minperm'. Generating entire set.
## Set of permutations < 'minperm'. Generating entire set.
```

```
## 'nperm' >= set of all permutations: complete enumeration.
```

```
## Set of permutations < 'minperm'. Generating entire set.
```

```
## 'nperm' >= set of all permutations: complete enumeration.
```

```
## Set of permutations < 'minperm'. Generating entire set.
## Set of permutations < 'minperm'. Generating entire set.
## Set of permutations < 'minperm'. Generating entire set.
```

```
## 'nperm' >= set of all permutations: complete enumeration.
```

```
## Set of permutations < 'minperm'. Generating entire set.
```

```
## 'nperm' >= set of all permutations: complete enumeration.
```

```
## Set of permutations < 'minperm'. Generating entire set.
## Set of permutations < 'minperm'. Generating entire set.
## Set of permutations < 'minperm'. Generating entire set.
## Set of permutations < 'minperm'. Generating entire set.
```

```
## 'nperm' >= set of all permutations: complete enumeration.
```

```
## Set of permutations < 'minperm'. Generating entire set.
```

```
## 'nperm' >= set of all permutations: complete enumeration.
```

```
## Set of permutations < 'minperm'. Generating entire set.
```

```
## 'nperm' >= set of all permutations: complete enumeration.
```

```
## Set of permutations < 'minperm'. Generating entire set.
```

```
## 'nperm' >= set of all permutations: complete enumeration.
```

```
## Set of permutations < 'minperm'. Generating entire set.
## Set of permutations < 'minperm'. Generating entire set.
```

```
## 'nperm' >= set of all permutations: complete enumeration.
```

```
## Set of permutations < 'minperm'. Generating entire set.
```

```
## 'nperm' >= set of all permutations: complete enumeration.
```

```
## Set of permutations < 'minperm'. Generating entire set.
## Set of permutations < 'minperm'. Generating entire set.
## Set of permutations < 'minperm'. Generating entire set.
```

```
## 'nperm' >= set of all permutations: complete enumeration.
```

```
## Set of permutations < 'minperm'. Generating entire set.
```

```r
tmp_df %>% 
  DT::datatable()
```

```{=html}
<div id="htmlwidget-5d05eac75555269eefd0" style="width:100%;height:auto;" class="datatables html-widget"></div>
<script type="application/json" data-for="htmlwidget-5d05eac75555269eefd0">{"x":{"filter":"none","vertical":false,"data":[["1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","23","24","25","26","27","28","29","30"],["bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson"],["CTX","CTX","CTX","CTX","CTX","CTX+HV292.1","CTX+HV292.1","CTX+HV292.1","CTX+HV292.1","HV292.1","HV292.1","HV292.1","VAN","VAN","VAN+CCUG59168","CTX","CTX","CTX","CTX","CTX","CTX+HV292.1","CTX+HV292.1","CTX+HV292.1","CTX+HV292.1","HV292.1","HV292.1","HV292.1","VAN","VAN","VAN+CCUG59168"],["CTX+HV292.1","HV292.1","VAN","VAN+CCUG59168","CCUG59168","HV292.1","VAN","VAN+CCUG59168","CCUG59168","VAN","VAN+CCUG59168","CCUG59168","VAN+CCUG59168","CCUG59168","CCUG59168","CTX+HV292.1","HV292.1","VAN","VAN+CCUG59168","CCUG59168","HV292.1","VAN","VAN+CCUG59168","CCUG59168","VAN","VAN+CCUG59168","CCUG59168","VAN+CCUG59168","CCUG59168","CCUG59168"],[0.306,0.362,0.503,0.462,0.457,0.496,0.608,0.546,0.555,0.463,0.486,0.496,0.555,0.508,0.312,0.316,0.362,0.463,0.502,0.482,0.498,0.549,0.577,0.554,0.442,0.496,0.508,0.523,0.489,0.332],[0.067,0.037,0.036,0.039,0.031,0.1,0.067,0.1,0.1,0.031,0.1,0.1,0.032,0.036,0.1,0.067,0.033,0.031,0.025,0.025,0.1,0.067,0.1,0.1,0.03,0.1,0.1,0.029,0.028,0.1],[0.118,0.092,0.12,0.084,0.31,0.12,0.118,0.12,0.12,0.31,0.12,0.12,0.16,0.12,0.12,0.118,0.071,0.078,0.25,0.25,0.12,0.118,0.12,0.12,0.09,0.12,0.12,0.109,0.14,0.12]],"container":"<table class=\"display\">\n  <thead>\n    <tr>\n      <th> <\/th>\n      <th>Distance<\/th>\n      <th>X1<\/th>\n      <th>X2<\/th>\n      <th>R2<\/th>\n      <th>pval<\/th>\n      <th>pvalFDR<\/th>\n    <\/tr>\n  <\/thead>\n<\/table>","options":{"columnDefs":[{"className":"dt-right","targets":[4,5,6]},{"orderable":false,"targets":0}],"order":[],"autoWidth":false,"orderClasses":false}},"evals":[],"jsHooks":[]}</script>
```

```r
tmp_df %>% 
  xlsx::write.xlsx(file = out_xlsx,
                   append = TRUE,
                   sheetName = "chicken1-beta-pretreat-pw",
                   showNA = TRUE)
```


```r
lapply(
  dist,
  FUN = phyloseq_TW,
  physeq =  ps_tmp_forbdiv %>%
    subset_samples(Day_of_Treatment <= 0),
  variable = "Treatment",
  nrep = 999,
  strata = NULL
) %>%
  bind_rows(.id = "Distance") %>%
  mutate_if(is.numeric, round, 3)   -> tmp_df

tmp_df %>% 
  DT::datatable()
```

```{=html}
<div id="htmlwidget-7531f021030c44f57ef2" style="width:100%;height:auto;" class="datatables html-widget"></div>
<script type="application/json" data-for="htmlwidget-7531f021030c44f57ef2">{"x":{"filter":"none","vertical":false,"data":[["1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","23","24","25","26","27","28","29","30"],["bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson"],[["CTX"],["CTX"],["CTX"],["CTX"],["CTX"],["CTX+HV292.1"],["CTX+HV292.1"],["CTX+HV292.1"],["CTX+HV292.1"],["HV292.1"],["HV292.1"],["HV292.1"],["VAN"],["VAN"],["VAN+CCUG59168"],["CTX"],["CTX"],["CTX"],["CTX"],["CTX"],["CTX+HV292.1"],["CTX+HV292.1"],["CTX+HV292.1"],["CTX+HV292.1"],["HV292.1"],["HV292.1"],["HV292.1"],["VAN"],["VAN"],["VAN+CCUG59168"]],[["CTX+HV292.1"],["HV292.1"],["VAN"],["VAN+CCUG59168"],["CCUG59168"],["HV292.1"],["VAN"],["VAN+CCUG59168"],["CCUG59168"],["VAN"],["VAN+CCUG59168"],["CCUG59168"],["VAN+CCUG59168"],["CCUG59168"],["CCUG59168"],["CTX+HV292.1"],["HV292.1"],["VAN"],["VAN+CCUG59168"],["CCUG59168"],["HV292.1"],["VAN"],["VAN+CCUG59168"],["CCUG59168"],["VAN"],["VAN+CCUG59168"],["CCUG59168"],["VAN+CCUG59168"],["CCUG59168"],["CCUG59168"]],[[4],[4],[4],[4],[4],[2],[2],[2],[2],[3],[3],[3],[4],[4],[3],[4],[4],[4],[4],[4],[2],[2],[2],[2],[3],[3],[3],[4],[4],[3]],[[2],[3],[4],[3],[3],[3],[4],[3],[3],[4],[3],[3],[3],[3],[3],[2],[3],[4],[3],[3],[3],[4],[3],[3],[4],[3],[3],[3],[3],[3]],[[0.058],[0.032],[0.035],[0.029],[0.036],[0.102],[0.061],[0.092],[0.111],[0.017],[0.093],[0.097],[0.03],[0.026],[0.105],[0.074],[0.032],[0.032],[0.031],[0.024],[0.096],[0.062],[0.102],[0.088],[0.032],[0.092],[0.099],[0.031],[0.024],[0.099]],[[2.27710074325349],[2.85508075078049],[6.06755951565703],[4.5759945581685],[4.39760069999656],[3.45477460027463],[6.59968339259196],[3.86031897047067],[4.13154705894282],[3.97735521358304],[3.78202389785935],[3.93597762106819],[6.12284632806336],[4.96092083370271],[1.81129161560545],[2.18099365598263],[2.82184326121165],[5.16695101411695],[5.2054202871342],[4.72781040354967],[3.33054021071813],[5.23212896166219],[4.35040649731935],[4.05296298539897],[3.78873679460449],[3.94390839748155],[4.13293292398528],[5.43094829979872],[4.66726089428945],[1.99128438539624]],[[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999]],[0.102,0.077,0.077,0.077,0.077,0.111,0.102,0.111,0.111,0.077,0.111,0.111,0.077,0.077,0.111,0.102,0.069,0.069,0.069,0.069,0.102,0.102,0.102,0.102,0.069,0.102,0.102,0.069,0.069,0.102]],"container":"<table class=\"display\">\n  <thead>\n    <tr>\n      <th> <\/th>\n      <th>Distance<\/th>\n      <th>Level1<\/th>\n      <th>Level2<\/th>\n      <th>N1<\/th>\n      <th>N2<\/th>\n      <th>p.value<\/th>\n      <th>tw2.stat<\/th>\n      <th>nrep<\/th>\n      <th>pvalFDR<\/th>\n    <\/tr>\n  <\/thead>\n<\/table>","options":{"columnDefs":[{"className":"dt-right","targets":9},{"orderable":false,"targets":0}],"order":[],"autoWidth":false,"orderClasses":false}},"evals":[],"jsHooks":[]}</script>
```

```r
tmp_df %>% 
  xlsx::write.xlsx(file = out_xlsx,
                   append = TRUE,
                   sheetName = "chicken1-beta-pretreat-pw2",
                   showNA = TRUE)
```



```r
# lapply(
#   dlist,
#   FUN = betadisper,
#   physeq = physeq,
#   variable = "Health_status"
# ) %>%
#   bind_rows() %>%
#   mutate_if(is.numeric, round, 3) %>%
#   t() %>%
#   data.frame() %>%
#   dplyr::rename(bdisper_pvalue = '.') %>%
#   DT::datatable()
```



- treatment - Day_of_Treatment >= 1



```r
lapply(
  dist,
  FUN = phyloseq_adonis2,
  physeq = ps_tmp_forbdiv %>%
    subset_samples(Day_of_Treatment >= 1),
  formula = paste0(c("Treatment"), collapse=" + "),
  nrep = 999,
  strata = "none",
  terms_margins = "margin"
) %>%
  bind_rows(.id = "Distance")  %>% 
  mutate_if(is.numeric, round, 3) -> tmp_df


tmp_df %>% 
  DT::datatable()
```

```{=html}
<div id="htmlwidget-983315589fed18d622a7" style="width:100%;height:auto;" class="datatables html-widget"></div>
<script type="application/json" data-for="htmlwidget-983315589fed18d622a7">{"x":{"filter":"none","vertical":false,"data":[["1","2","3","4","5","6"],["bjaccard","bjaccard","bjaccard","aichinson","aichinson","aichinson"],["Treatment","Residual","Total","Treatment","Residual","Total"],[5,70,75,5,70,75],[10.589,8.322,18.912,3.827,3.398,7.226],[0.56,0.44,1,0.53,0.47,1],[17.814,null,null,15.767,null,null],[0.001,null,null,0.001,null,null]],"container":"<table class=\"display\">\n  <thead>\n    <tr>\n      <th> <\/th>\n      <th>Distance<\/th>\n      <th>terms<\/th>\n      <th>Df<\/th>\n      <th>SumOfSqs<\/th>\n      <th>R2<\/th>\n      <th>F<\/th>\n      <th>Pr(&gt;F)<\/th>\n    <\/tr>\n  <\/thead>\n<\/table>","options":{"columnDefs":[{"className":"dt-right","targets":[3,4,5,6,7]},{"orderable":false,"targets":0}],"order":[],"autoWidth":false,"orderClasses":false}},"evals":[],"jsHooks":[]}</script>
```

```r
tmp_df %>% 
  xlsx::write.xlsx(file = out_xlsx,
                   append = TRUE,
                   sheetName = "chicken1-beta-treat",
                   showNA = TRUE)
```


```r
lapply(
  dist,
  FUN = phyloseq_adonis2,
  physeq = ps_tmp_forbdiv %>%
    subset_samples(Day_of_Treatment >= 1),
  formula = paste0(c("Treatment", "Day_of_Treatment"), collapse=" + "),
  nrep = 999,
  strata = "none",
  terms_margins = "margin"
) %>%
  bind_rows(.id = "Distance")  %>% 
  mutate_if(is.numeric, round, 3) -> tmp_df


tmp_df %>% 
  DT::datatable()
```

```{=html}
<div id="htmlwidget-6f1d8b9e786d9d3fe28d" style="width:100%;height:auto;" class="datatables html-widget"></div>
<script type="application/json" data-for="htmlwidget-6f1d8b9e786d9d3fe28d">{"x":{"filter":"none","vertical":false,"data":[["1","2","3","4","5","6","7","8"],["bjaccard","bjaccard","bjaccard","bjaccard","aichinson","aichinson","aichinson","aichinson"],["Treatment","Day_of_Treatment","Residual","Total","Treatment","Day_of_Treatment","Residual","Total"],[5,1,69,75,5,1,69,75],[10.542,0.701,7.621,18.912,3.815,0.286,3.112,7.226],[0.557,0.037,0.403,1,0.528,0.04,0.431,1],[19.09,6.347,null,null,16.916,6.347,null,null],[0.001,0.002,null,null,0.001,0.001,null,null]],"container":"<table class=\"display\">\n  <thead>\n    <tr>\n      <th> <\/th>\n      <th>Distance<\/th>\n      <th>terms<\/th>\n      <th>Df<\/th>\n      <th>SumOfSqs<\/th>\n      <th>R2<\/th>\n      <th>F<\/th>\n      <th>Pr(&gt;F)<\/th>\n    <\/tr>\n  <\/thead>\n<\/table>","options":{"columnDefs":[{"className":"dt-right","targets":[3,4,5,6,7]},{"orderable":false,"targets":0}],"order":[],"autoWidth":false,"orderClasses":false}},"evals":[],"jsHooks":[]}</script>
```

```r
tmp_df %>% 
  xlsx::write.xlsx(file = out_xlsx,
                   append = TRUE,
                   sheetName = "chicken1-beta-treat-time-treat",
                   showNA = TRUE)
```



```r
lapply(
  dist,
  FUN = physeq_pairwise_permanovas_adonis2,
  physeq = ps_tmp_forbdiv %>%
    subset_samples(Day_of_Treatment >= 1),
  compare_header = "Treatment",
  n_perm = 999,
  strat = "none"
) %>%
  bind_rows(.id = "Distance") %>%
  mutate_if(is.numeric, round, 3) %>%
  select(-pvalBon)  -> tmp_df


tmp_df %>% 
  DT::datatable()
```

```{=html}
<div id="htmlwidget-dcdfc64296f5794e54a4" style="width:100%;height:auto;" class="datatables html-widget"></div>
<script type="application/json" data-for="htmlwidget-dcdfc64296f5794e54a4">{"x":{"filter":"none","vertical":false,"data":[["1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","23","24","25","26","27","28","29","30"],["bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson"],["CTX","CTX","CTX","CTX","CTX","CTX+HV292.1","CTX+HV292.1","CTX+HV292.1","CTX+HV292.1","HV292.1","HV292.1","HV292.1","VAN","VAN","VAN+CCUG59168","CTX","CTX","CTX","CTX","CTX","CTX+HV292.1","CTX+HV292.1","CTX+HV292.1","CTX+HV292.1","HV292.1","HV292.1","HV292.1","VAN","VAN","VAN+CCUG59168"],["CTX+HV292.1","HV292.1","VAN","VAN+CCUG59168","CCUG59168","HV292.1","VAN","VAN+CCUG59168","CCUG59168","VAN","VAN+CCUG59168","CCUG59168","VAN+CCUG59168","CCUG59168","CCUG59168","CTX+HV292.1","HV292.1","VAN","VAN+CCUG59168","CCUG59168","HV292.1","VAN","VAN+CCUG59168","CCUG59168","VAN","VAN+CCUG59168","CCUG59168","VAN+CCUG59168","CCUG59168","CCUG59168"],[0.19,0.433,0.526,0.448,0.4,0.23,0.532,0.45,0.217,0.633,0.535,0.12,0.083,0.587,0.491,0.178,0.441,0.499,0.433,0.421,0.238,0.473,0.407,0.238,0.623,0.532,0.156,0.078,0.561,0.468],[0.001,0.001,0.001,0.001,0.001,0.001,0.001,0.001,0.001,0.001,0.001,0.001,0.093,0.001,0.001,0.001,0.001,0.001,0.001,0.001,0.001,0.001,0.001,0.001,0.001,0.001,0.001,0.11,0.001,0.001],[0.002,0.002,0.002,0.002,0.002,0.002,0.002,0.002,0.002,0.002,0.002,0.002,0.093,0.002,0.002,0.002,0.002,0.002,0.002,0.002,0.002,0.002,0.002,0.002,0.002,0.002,0.002,0.11,0.002,0.002]],"container":"<table class=\"display\">\n  <thead>\n    <tr>\n      <th> <\/th>\n      <th>Distance<\/th>\n      <th>X1<\/th>\n      <th>X2<\/th>\n      <th>R2<\/th>\n      <th>pval<\/th>\n      <th>pvalFDR<\/th>\n    <\/tr>\n  <\/thead>\n<\/table>","options":{"columnDefs":[{"className":"dt-right","targets":[4,5,6]},{"orderable":false,"targets":0}],"order":[],"autoWidth":false,"orderClasses":false}},"evals":[],"jsHooks":[]}</script>
```

```r
tmp_df %>% 
  xlsx::write.xlsx(file = out_xlsx,
                   append = TRUE,
                   sheetName = "chicken1-beta-treat-pw",
                   showNA = TRUE)
```


```r
lapply(
  dist,
  FUN = phyloseq_TW,
  physeq =  ps_tmp_forbdiv %>%
    subset_samples(Day_of_Treatment >= 1),
  variable = "Treatment",
  nrep = 999,
  strata = NULL
) %>%
  bind_rows(.id = "Distance") %>%
  mutate_if(is.numeric, round, 3) -> tmp_df

tmp_df %>% 
  DT::datatable()
```

```{=html}
<div id="htmlwidget-873249b169eb6149d57e" style="width:100%;height:auto;" class="datatables html-widget"></div>
<script type="application/json" data-for="htmlwidget-873249b169eb6149d57e">{"x":{"filter":"none","vertical":false,"data":[["1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","23","24","25","26","27","28","29","30"],["bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson"],[["CTX"],["CTX"],["CTX"],["CTX"],["CTX"],["CTX+HV292.1"],["CTX+HV292.1"],["CTX+HV292.1"],["CTX+HV292.1"],["HV292.1"],["HV292.1"],["HV292.1"],["VAN"],["VAN"],["VAN+CCUG59168"],["CTX"],["CTX"],["CTX"],["CTX"],["CTX"],["CTX+HV292.1"],["CTX+HV292.1"],["CTX+HV292.1"],["CTX+HV292.1"],["HV292.1"],["HV292.1"],["HV292.1"],["VAN"],["VAN"],["VAN+CCUG59168"]],[["CTX+HV292.1"],["HV292.1"],["VAN"],["VAN+CCUG59168"],["CCUG59168"],["HV292.1"],["VAN"],["VAN+CCUG59168"],["CCUG59168"],["VAN"],["VAN+CCUG59168"],["CCUG59168"],["VAN+CCUG59168"],["CCUG59168"],["CCUG59168"],["CTX+HV292.1"],["HV292.1"],["VAN"],["VAN+CCUG59168"],["CCUG59168"],["HV292.1"],["VAN"],["VAN+CCUG59168"],["CCUG59168"],["VAN"],["VAN+CCUG59168"],["CCUG59168"],["VAN+CCUG59168"],["CCUG59168"],["CCUG59168"]],[[12],[12],[12],[12],[12],[14],[14],[14],[14],[12],[12],[12],[12],[12],[13],[12],[12],[12],[12],[12],[14],[14],[14],[14],[12],[12],[12],[12],[12],[13]],[[14],[12],[12],[13],[13],[12],[12],[13],[13],[12],[13],[13],[13],[13],[13],[14],[12],[12],[13],[13],[12],[12],[13],[13],[12],[13],[13],[13],[13],[13]],[[0.001],[0.001],[0.001],[0.001],[0.001],[0.001],[0.001],[0.001],[0.001],[0.001],[0.001],[0.001],[0.088],[0.001],[0.001],[0.002],[0.001],[0.001],[0.001],[0.001],[0.001],[0.001],[0.001],[0.001],[0.001],[0.001],[0.001],[0.135],[0.001],[0.001]],[[5.62955129634188],[16.7952015586551],[24.3774111449177],[18.9133588910552],[15.0630402306501],[7.56851568773938],[27.238175811617],[20.2311526500657],[7.04472009045441],[37.9409559570189],[27.43358657796],[3.17603236638448],[2.10981436224055],[32.0616652395467],[23.118248852649],[5.27858817351732],[17.3433235461057],[21.9344865173774],[17.3834808735498],[16.5245287840095],[7.83170469880149],[23.1494740237639],[17.4952782190865],[7.93553139200452],[36.2795867910752],[26.2693927092316],[4.28164708409463],[1.99522830485],[29.9675165551828],[21.1414435702442]],[[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999]],[0.001,0.001,0.001,0.001,0.001,0.001,0.001,0.001,0.001,0.001,0.001,0.001,0.088,0.001,0.001,0.002,0.001,0.001,0.001,0.001,0.001,0.001,0.001,0.001,0.001,0.001,0.001,0.135,0.001,0.001]],"container":"<table class=\"display\">\n  <thead>\n    <tr>\n      <th> <\/th>\n      <th>Distance<\/th>\n      <th>Level1<\/th>\n      <th>Level2<\/th>\n      <th>N1<\/th>\n      <th>N2<\/th>\n      <th>p.value<\/th>\n      <th>tw2.stat<\/th>\n      <th>nrep<\/th>\n      <th>pvalFDR<\/th>\n    <\/tr>\n  <\/thead>\n<\/table>","options":{"columnDefs":[{"className":"dt-right","targets":9},{"orderable":false,"targets":0}],"order":[],"autoWidth":false,"orderClasses":false}},"evals":[],"jsHooks":[]}</script>
```

```r
tmp_df %>% 
  xlsx::write.xlsx(file = out_xlsx,
                   append = TRUE,
                   sheetName = "chicken1-beta-treat-pw2",
                   showNA = TRUE)
```

- filter(Day_of_Treatment >= 1 & Day_of_Treatment <= 6) %>% 



```r
lapply(
  dist,
  FUN = phyloseq_adonis2,
  physeq = ps_tmp_forbdiv %>%
    subset_samples(Day_of_Treatment >= 1 & Day_of_Treatment <= 6),
  formula = paste0(c("Treatment"), collapse=" + "),
  nrep = 999,
  strata = "none",
  terms_margins = "margin"
) %>%
  bind_rows(.id = "Distance")  %>% 
  mutate_if(is.numeric, round, 3) -> tmp_df


tmp_df %>% 
  DT::datatable()
```

```{=html}
<div id="htmlwidget-cbd3ad10bd2d35a2e76d" style="width:100%;height:auto;" class="datatables html-widget"></div>
<script type="application/json" data-for="htmlwidget-cbd3ad10bd2d35a2e76d">{"x":{"filter":"none","vertical":false,"data":[["1","2","3","4","5","6"],["bjaccard","bjaccard","bjaccard","aichinson","aichinson","aichinson"],["Treatment","Residual","Total","Treatment","Residual","Total"],[5,29,34,5,29,34],[4.102,3.863,7.965,1.685,1.398,3.082],[0.515,0.485,1,0.547,0.453,1],[6.158,null,null,6.991,null,null],[0.001,null,null,0.001,null,null]],"container":"<table class=\"display\">\n  <thead>\n    <tr>\n      <th> <\/th>\n      <th>Distance<\/th>\n      <th>terms<\/th>\n      <th>Df<\/th>\n      <th>SumOfSqs<\/th>\n      <th>R2<\/th>\n      <th>F<\/th>\n      <th>Pr(&gt;F)<\/th>\n    <\/tr>\n  <\/thead>\n<\/table>","options":{"columnDefs":[{"className":"dt-right","targets":[3,4,5,6,7]},{"orderable":false,"targets":0}],"order":[],"autoWidth":false,"orderClasses":false}},"evals":[],"jsHooks":[]}</script>
```

```r
tmp_df %>% 
  xlsx::write.xlsx(file = out_xlsx,
                   append = TRUE,
                   sheetName = "chicken1-beta-d1-d6",
                   showNA = TRUE)
```


```r
lapply(
  dist,
  FUN = phyloseq_adonis2,
  physeq = ps_tmp_forbdiv %>%
    subset_samples(Day_of_Treatment >= 1 & Day_of_Treatment <= 6),
  formula = paste0(c("Treatment", "Day_of_Treatment"), collapse=" + "),
  nrep = 999,
  strata = "none",
  terms_margins = "margin"
) %>%
  bind_rows(.id = "Distance")  %>% 
  mutate_if(is.numeric, round, 3) -> tmp_df


tmp_df %>% 
  DT::datatable()
```

```{=html}
<div id="htmlwidget-f25d4e1770d66a86dcce" style="width:100%;height:auto;" class="datatables html-widget"></div>
<script type="application/json" data-for="htmlwidget-f25d4e1770d66a86dcce">{"x":{"filter":"none","vertical":false,"data":[["1","2","3","4","5","6","7","8"],["bjaccard","bjaccard","bjaccard","bjaccard","aichinson","aichinson","aichinson","aichinson"],["Treatment","Day_of_Treatment","Residual","Total","Treatment","Day_of_Treatment","Residual","Total"],[5,1,28,34,5,1,28,34],[4.028,0.664,3.199,7.965,1.664,0.228,1.169,3.082],[0.506,0.083,0.402,1,0.54,0.074,0.379,1],[7.051,5.813,null,null,7.971,5.466,null,null],[0.001,0.002,null,null,0.001,0.002,null,null]],"container":"<table class=\"display\">\n  <thead>\n    <tr>\n      <th> <\/th>\n      <th>Distance<\/th>\n      <th>terms<\/th>\n      <th>Df<\/th>\n      <th>SumOfSqs<\/th>\n      <th>R2<\/th>\n      <th>F<\/th>\n      <th>Pr(&gt;F)<\/th>\n    <\/tr>\n  <\/thead>\n<\/table>","options":{"columnDefs":[{"className":"dt-right","targets":[3,4,5,6,7]},{"orderable":false,"targets":0}],"order":[],"autoWidth":false,"orderClasses":false}},"evals":[],"jsHooks":[]}</script>
```

```r
tmp_df %>% 
  xlsx::write.xlsx(file = out_xlsx,
                   append = TRUE,
                   sheetName = "chicken1-beta-d1-d6-time-treat",
                   showNA = TRUE)
```


```r
lapply(
  dist,
  FUN = physeq_pairwise_permanovas_adonis2,
  physeq = ps_tmp_forbdiv %>%
    subset_samples(Day_of_Treatment >= 1 & Day_of_Treatment <= 6),
  compare_header = "Treatment",
  n_perm = 999,
  strat = "none"
) %>%
  bind_rows(.id = "Distance") %>%
  mutate_if(is.numeric, round, 3) %>%
  select(-pvalBon)  -> tmp_df


tmp_df %>% 
  DT::datatable()
```

```{=html}
<div id="htmlwidget-969965edfe9f869aa904" style="width:100%;height:auto;" class="datatables html-widget"></div>
<script type="application/json" data-for="htmlwidget-969965edfe9f869aa904">{"x":{"filter":"none","vertical":false,"data":[["1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","23","24","25","26","27","28","29","30"],["bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson"],["CTX","CTX","CTX","CTX","CTX","CTX+HV292.1","CTX+HV292.1","CTX+HV292.1","CTX+HV292.1","HV292.1","HV292.1","HV292.1","VAN","VAN","VAN+CCUG59168","CTX","CTX","CTX","CTX","CTX","CTX+HV292.1","CTX+HV292.1","CTX+HV292.1","CTX+HV292.1","HV292.1","HV292.1","HV292.1","VAN","VAN","VAN+CCUG59168"],["CTX+HV292.1","HV292.1","VAN","VAN+CCUG59168","CCUG59168","HV292.1","VAN","VAN+CCUG59168","CCUG59168","VAN","VAN+CCUG59168","CCUG59168","VAN+CCUG59168","CCUG59168","CCUG59168","CTX+HV292.1","HV292.1","VAN","VAN+CCUG59168","CCUG59168","HV292.1","VAN","VAN+CCUG59168","CCUG59168","VAN","VAN+CCUG59168","CCUG59168","VAN+CCUG59168","CCUG59168","CCUG59168"],[0.143,0.49,0.46,0.342,0.428,0.441,0.478,0.346,0.374,0.599,0.419,0.258,0.104,0.509,0.338,0.143,0.501,0.453,0.372,0.464,0.449,0.474,0.377,0.414,0.643,0.489,0.296,0.099,0.543,0.382],[0.122,0.002,0.002,0.004,0.001,0.004,0.004,0.002,0.002,0.004,0.006,0.004,0.339,0.006,0.003,0.117,0.004,0.003,0.004,0.004,0.002,0.002,0.004,0.002,0.007,0.001,0.003,0.285,0.004,0.01],[0.131,0.009,0.009,0.007,0.015,0.007,0.007,0.009,0.009,0.007,0.007,0.007,0.339,0.007,0.007,0.125,0.007,0.008,0.007,0.007,0.01,0.01,0.007,0.01,0.009,0.015,0.008,0.285,0.007,0.012]],"container":"<table class=\"display\">\n  <thead>\n    <tr>\n      <th> <\/th>\n      <th>Distance<\/th>\n      <th>X1<\/th>\n      <th>X2<\/th>\n      <th>R2<\/th>\n      <th>pval<\/th>\n      <th>pvalFDR<\/th>\n    <\/tr>\n  <\/thead>\n<\/table>","options":{"columnDefs":[{"className":"dt-right","targets":[4,5,6]},{"orderable":false,"targets":0}],"order":[],"autoWidth":false,"orderClasses":false}},"evals":[],"jsHooks":[]}</script>
```

```r
tmp_df %>% 
  xlsx::write.xlsx(file = out_xlsx,
                   append = TRUE,
                   sheetName = "chicken1-beta-d1-d6-pw",
                   showNA = TRUE)
```


```r
lapply(
  dist,
  FUN = phyloseq_TW,
  physeq =  ps_tmp_forbdiv %>%
    subset_samples(Day_of_Treatment >= 1 & Day_of_Treatment <= 6),
  variable = "Treatment",
  nrep = 999,
  strata = NULL
) %>%
  bind_rows(.id = "Distance") %>%
  mutate_if(is.numeric, round, 3) -> tmp_df

tmp_df %>% 
  DT::datatable()
```

```{=html}
<div id="htmlwidget-450f171f14364aa3ebed" style="width:100%;height:auto;" class="datatables html-widget"></div>
<script type="application/json" data-for="htmlwidget-450f171f14364aa3ebed">{"x":{"filter":"none","vertical":false,"data":[["1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","23","24","25","26","27","28","29","30"],["bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson"],[["CTX"],["CTX"],["CTX"],["CTX"],["CTX"],["CTX+HV292.1"],["CTX+HV292.1"],["CTX+HV292.1"],["CTX+HV292.1"],["HV292.1"],["HV292.1"],["HV292.1"],["VAN"],["VAN"],["VAN+CCUG59168"],["CTX"],["CTX"],["CTX"],["CTX"],["CTX"],["CTX+HV292.1"],["CTX+HV292.1"],["CTX+HV292.1"],["CTX+HV292.1"],["HV292.1"],["HV292.1"],["HV292.1"],["VAN"],["VAN"],["VAN+CCUG59168"]],[["CTX+HV292.1"],["HV292.1"],["VAN"],["VAN+CCUG59168"],["CCUG59168"],["HV292.1"],["VAN"],["VAN+CCUG59168"],["CCUG59168"],["VAN"],["VAN+CCUG59168"],["CCUG59168"],["VAN+CCUG59168"],["CCUG59168"],["CCUG59168"],["CTX+HV292.1"],["HV292.1"],["VAN"],["VAN+CCUG59168"],["CCUG59168"],["HV292.1"],["VAN"],["VAN+CCUG59168"],["CCUG59168"],["VAN"],["VAN+CCUG59168"],["CCUG59168"],["VAN+CCUG59168"],["CCUG59168"],["CCUG59168"]],[[6],[6],[6],[6],[6],[6],[6],[6],[6],[6],[6],[6],[5],[5],[6],[6],[6],[6],[6],[6],[6],[6],[6],[6],[6],[6],[6],[5],[5],[6]],[[6],[6],[5],[6],[6],[6],[5],[6],[6],[5],[6],[6],[6],[6],[6],[6],[6],[5],[6],[6],[6],[5],[6],[6],[5],[6],[6],[6],[6],[6]],[[0.111],[0.001],[0.003],[0.002],[0.002],[0.001],[0.001],[0.002],[0.005],[0.004],[0.007],[0.005],[0.301],[0.003],[0.006],[0.119],[0.003],[0.002],[0.003],[0.004],[0.004],[0.002],[0.002],[0.005],[0.003],[0.003],[0.002],[0.265],[0.005],[0.002]],[[1.66728111648946],[9.61717754270057],[7.40412845088531],[5.20675415571432],[7.49224585593084],[7.88516604784861],[7.77867292840547],[5.30040577727223],[5.98218888117334],[11.7038434115344],[7.19952821467173],[3.48186732199651],[1.07017386961333],[8.53479135974693],[5.10336283534725],[1.67138288099198],[10.0535929548539],[7.78212531000009],[5.91421765055282],[8.67150711918878],[8.16366339685637],[8.40176041478224],[6.04298306297006],[7.0741531263652],[15.6355324798611],[9.55514353702956],[4.21172161673649],[1.03849451729933],[10.6558104321424],[6.18033715331262]],[[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999]],[0.119,0.005,0.006,0.005,0.005,0.005,0.005,0.005,0.007,0.007,0.008,0.007,0.301,0.006,0.007,0.128,0.005,0.005,0.005,0.005,0.005,0.005,0.005,0.006,0.005,0.005,0.005,0.265,0.006,0.005]],"container":"<table class=\"display\">\n  <thead>\n    <tr>\n      <th> <\/th>\n      <th>Distance<\/th>\n      <th>Level1<\/th>\n      <th>Level2<\/th>\n      <th>N1<\/th>\n      <th>N2<\/th>\n      <th>p.value<\/th>\n      <th>tw2.stat<\/th>\n      <th>nrep<\/th>\n      <th>pvalFDR<\/th>\n    <\/tr>\n  <\/thead>\n<\/table>","options":{"columnDefs":[{"className":"dt-right","targets":9},{"orderable":false,"targets":0}],"order":[],"autoWidth":false,"orderClasses":false}},"evals":[],"jsHooks":[]}</script>
```

```r
tmp_df %>% 
  xlsx::write.xlsx(file = out_xlsx,
                   append = TRUE,
                   sheetName = "chicken1-beta-d1-d6-pw2",
                   showNA = TRUE)
```


- 7-end



```r
lapply(
  dist,
  FUN = phyloseq_adonis2,
  physeq = ps_tmp_forbdiv %>%
    subset_samples(Day_of_Treatment > 7),
  formula = paste0(c("Treatment"), collapse=" + "),
  nrep = 999,
  strata = "none",
  terms_margins = "margin"
) %>%
  bind_rows(.id = "Distance")  %>% 
  mutate_if(is.numeric, round, 3) -> tmp_df


tmp_df %>% 
  DT::datatable()
```

```{=html}
<div id="htmlwidget-fbfc52d82635d988c6e1" style="width:100%;height:auto;" class="datatables html-widget"></div>
<script type="application/json" data-for="htmlwidget-fbfc52d82635d988c6e1">{"x":{"filter":"none","vertical":false,"data":[["1","2","3","4","5","6"],["bjaccard","bjaccard","bjaccard","aichinson","aichinson","aichinson"],["Treatment","Residual","Total","Treatment","Residual","Total"],[5,35,40,5,35,40],[8.045,2.094,10.139,2.743,1.103,3.846],[0.793,0.207,1,0.713,0.287,1],[26.892,null,null,17.417,null,null],[0.001,null,null,0.001,null,null]],"container":"<table class=\"display\">\n  <thead>\n    <tr>\n      <th> <\/th>\n      <th>Distance<\/th>\n      <th>terms<\/th>\n      <th>Df<\/th>\n      <th>SumOfSqs<\/th>\n      <th>R2<\/th>\n      <th>F<\/th>\n      <th>Pr(&gt;F)<\/th>\n    <\/tr>\n  <\/thead>\n<\/table>","options":{"columnDefs":[{"className":"dt-right","targets":[3,4,5,6,7]},{"orderable":false,"targets":0}],"order":[],"autoWidth":false,"orderClasses":false}},"evals":[],"jsHooks":[]}</script>
```

```r
tmp_df %>% 
  xlsx::write.xlsx(file = out_xlsx,
                   append = TRUE,
                   sheetName = "chicken1-beta-d7-end",
                   showNA = TRUE)
```




```r
lapply(
  dist,
  FUN = phyloseq_adonis2,
  physeq = ps_tmp_forbdiv %>%
    subset_samples(Day_of_Treatment > 7),
  formula = paste0(c("Treatment", "Day_of_Treatment"), collapse=" + "),
  nrep = 999,
  strata = "none",
  terms_margins = "margin"
) %>%
  bind_rows(.id = "Distance")  %>% 
  mutate_if(is.numeric, round, 3) -> tmp_df


tmp_df %>% 
  DT::datatable()
```

```{=html}
<div id="htmlwidget-2f4c52b7f36e26c300a7" style="width:100%;height:auto;" class="datatables html-widget"></div>
<script type="application/json" data-for="htmlwidget-2f4c52b7f36e26c300a7">{"x":{"filter":"none","vertical":false,"data":[["1","2","3","4","5","6","7","8"],["bjaccard","bjaccard","bjaccard","bjaccard","aichinson","aichinson","aichinson","aichinson"],["Treatment","Day_of_Treatment","Residual","Total","Treatment","Day_of_Treatment","Residual","Total"],[5,1,34,40,5,1,34,40],[8.041,0.129,1.966,10.139,2.743,0.067,1.035,3.846],[0.793,0.013,0.194,1,0.713,0.018,0.269,1],[27.818,2.225,null,null,18.017,2.214,null,null],[0.001,0.088,null,null,0.001,0.06,null,null]],"container":"<table class=\"display\">\n  <thead>\n    <tr>\n      <th> <\/th>\n      <th>Distance<\/th>\n      <th>terms<\/th>\n      <th>Df<\/th>\n      <th>SumOfSqs<\/th>\n      <th>R2<\/th>\n      <th>F<\/th>\n      <th>Pr(&gt;F)<\/th>\n    <\/tr>\n  <\/thead>\n<\/table>","options":{"columnDefs":[{"className":"dt-right","targets":[3,4,5,6,7]},{"orderable":false,"targets":0}],"order":[],"autoWidth":false,"orderClasses":false}},"evals":[],"jsHooks":[]}</script>
```

```r
tmp_df %>% 
  xlsx::write.xlsx(file = out_xlsx,
                   append = TRUE,
                   sheetName = "chicken1-beta-d7-end-time-treat",
                   showNA = TRUE)
```



```r
lapply(
  dist,
  FUN = physeq_pairwise_permanovas_adonis2,
  physeq = ps_tmp_forbdiv %>%
    subset_samples(Day_of_Treatment > 7),
  compare_header = "Treatment",
  n_perm = 999,
  strat = "none"
) %>%
  bind_rows(.id = "Distance") %>%
  mutate_if(is.numeric, round, 3) %>%
  select(-pvalBon)  -> tmp_df


tmp_df %>% 
  DT::datatable()
```

```{=html}
<div id="htmlwidget-1d5df3d80f5735caaf13" style="width:100%;height:auto;" class="datatables html-widget"></div>
<script type="application/json" data-for="htmlwidget-1d5df3d80f5735caaf13">{"x":{"filter":"none","vertical":false,"data":[["1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","23","24","25","26","27","28","29","30"],["bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson"],["CTX","CTX","CTX","CTX","CTX","CTX+HV292.1","CTX+HV292.1","CTX+HV292.1","CTX+HV292.1","HV292.1","HV292.1","HV292.1","VAN","VAN","VAN+CCUG59168","CTX","CTX","CTX","CTX","CTX","CTX+HV292.1","CTX+HV292.1","CTX+HV292.1","CTX+HV292.1","HV292.1","HV292.1","HV292.1","VAN","VAN","VAN+CCUG59168"],["CTX+HV292.1","HV292.1","VAN","VAN+CCUG59168","CCUG59168","HV292.1","VAN","VAN+CCUG59168","CCUG59168","VAN","VAN+CCUG59168","CCUG59168","VAN+CCUG59168","CCUG59168","CCUG59168","CTX+HV292.1","HV292.1","VAN","VAN+CCUG59168","CCUG59168","HV292.1","VAN","VAN+CCUG59168","CCUG59168","VAN","VAN+CCUG59168","CCUG59168","VAN+CCUG59168","CCUG59168","CCUG59168"],[0.414,0.553,0.784,0.794,0.538,0.298,0.779,0.784,0.315,0.842,0.852,0.167,0.492,0.832,0.84,0.388,0.575,0.726,0.734,0.549,0.304,0.662,0.663,0.309,0.774,0.781,0.219,0.41,0.734,0.738],[0.001,0.001,0.002,0.002,0.001,0.003,0.001,0.001,0.001,0.002,0.002,0.03,0.003,0.001,0.001,0.001,0.003,0.001,0.001,0.001,0.005,0.001,0.001,0.001,0.002,0.001,0.002,0.002,0.002,0.001],[0.003,0.003,0.003,0.003,0.003,0.003,0.003,0.003,0.003,0.003,0.003,0.03,0.003,0.003,0.003,0.003,0.003,0.003,0.003,0.003,0.005,0.003,0.003,0.003,0.003,0.003,0.003,0.003,0.003,0.003]],"container":"<table class=\"display\">\n  <thead>\n    <tr>\n      <th> <\/th>\n      <th>Distance<\/th>\n      <th>X1<\/th>\n      <th>X2<\/th>\n      <th>R2<\/th>\n      <th>pval<\/th>\n      <th>pvalFDR<\/th>\n    <\/tr>\n  <\/thead>\n<\/table>","options":{"columnDefs":[{"className":"dt-right","targets":[4,5,6]},{"orderable":false,"targets":0}],"order":[],"autoWidth":false,"orderClasses":false}},"evals":[],"jsHooks":[]}</script>
```

```r
tmp_df %>% 
  xlsx::write.xlsx(file = out_xlsx,
                   append = TRUE,
                   sheetName = "chicken1-beta-d7-end-pw",
                   showNA = TRUE)
```


```r
lapply(
  dist,
  FUN = phyloseq_TW,
  physeq =  ps_tmp_forbdiv %>%
    subset_samples(Day_of_Treatment > 7),
  variable = "Treatment",
  nrep = 999,
  strata = NULL
) %>%
  bind_rows(.id = "Distance") %>%
  mutate_if(is.numeric, round, 3) -> tmp_df

tmp_df %>% 
  DT::datatable()
```

```{=html}
<div id="htmlwidget-67c25ca3336b15b0293a" style="width:100%;height:auto;" class="datatables html-widget"></div>
<script type="application/json" data-for="htmlwidget-67c25ca3336b15b0293a">{"x":{"filter":"none","vertical":false,"data":[["1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","23","24","25","26","27","28","29","30"],["bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson"],[["CTX"],["CTX"],["CTX"],["CTX"],["CTX"],["CTX+HV292.1"],["CTX+HV292.1"],["CTX+HV292.1"],["CTX+HV292.1"],["HV292.1"],["HV292.1"],["HV292.1"],["VAN"],["VAN"],["VAN+CCUG59168"],["CTX"],["CTX"],["CTX"],["CTX"],["CTX"],["CTX+HV292.1"],["CTX+HV292.1"],["CTX+HV292.1"],["CTX+HV292.1"],["HV292.1"],["HV292.1"],["HV292.1"],["VAN"],["VAN"],["VAN+CCUG59168"]],[["CTX+HV292.1"],["HV292.1"],["VAN"],["VAN+CCUG59168"],["CCUG59168"],["HV292.1"],["VAN"],["VAN+CCUG59168"],["CCUG59168"],["VAN"],["VAN+CCUG59168"],["CCUG59168"],["VAN+CCUG59168"],["CCUG59168"],["CCUG59168"],["CTX+HV292.1"],["HV292.1"],["VAN"],["VAN+CCUG59168"],["CCUG59168"],["HV292.1"],["VAN"],["VAN+CCUG59168"],["CCUG59168"],["VAN"],["VAN+CCUG59168"],["CCUG59168"],["VAN+CCUG59168"],["CCUG59168"],["CCUG59168"]],[[6],[6],[6],[6],[6],[8],[8],[8],[8],[6],[6],[6],[7],[7],[7],[6],[6],[6],[6],[6],[8],[8],[8],[8],[6],[6],[6],[7],[7],[7]],[[8],[6],[7],[7],[7],[6],[7],[7],[7],[7],[7],[7],[7],[7],[7],[8],[6],[7],[7],[7],[6],[7],[7],[7],[7],[7],[7],[7],[7],[7]],[[0.001],[0.002],[0.001],[0.002],[0.002],[0.004],[0.001],[0.001],[0.001],[0.001],[0.001],[0.017],[0.001],[0.002],[0.002],[0.001],[0.002],[0.002],[0.002],[0.001],[0.001],[0.001],[0.001],[0.001],[0.001],[0.002],[0.002],[0.001],[0.001],[0.001]],[[8.53313004303612],[12.350357047768],[35.9289205528397],[37.3835575527423],[12.4187469773001],[5.3810240500763],[50.2030745072786],[52.5205242437094],[6.15782728203866],[53.8124158380675],[56.4010095673785],[2.19724259420943],[11.6291161697717],[59.5563850668537],[63.080527847121],[8.11190134376749],[13.5253650356858],[25.7852884172108],[26.5923328335292],[13.3391208653124],[5.61268748642504],[28.5751877124528],[28.8883765402037],[6.00894014083955],[33.4999263783951],[34.382756635146],[3.0953930944048],[8.34136639507506],[33.0319306964896],[33.8498337466393]],[[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999]],[0.002,0.002,0.002,0.002,0.002,0.004,0.002,0.002,0.002,0.002,0.002,0.017,0.002,0.002,0.002,0.002,0.002,0.002,0.002,0.002,0.002,0.002,0.002,0.002,0.002,0.002,0.002,0.002,0.002,0.002]],"container":"<table class=\"display\">\n  <thead>\n    <tr>\n      <th> <\/th>\n      <th>Distance<\/th>\n      <th>Level1<\/th>\n      <th>Level2<\/th>\n      <th>N1<\/th>\n      <th>N2<\/th>\n      <th>p.value<\/th>\n      <th>tw2.stat<\/th>\n      <th>nrep<\/th>\n      <th>pvalFDR<\/th>\n    <\/tr>\n  <\/thead>\n<\/table>","options":{"columnDefs":[{"className":"dt-right","targets":9},{"orderable":false,"targets":0}],"order":[],"autoWidth":false,"orderClasses":false}},"evals":[],"jsHooks":[]}</script>
```

```r
tmp_df %>% 
  xlsx::write.xlsx(file = out_xlsx,
                   append = TRUE,
                   sheetName = "chicken1-beta-d7-end-pw2",
                   showNA = TRUE)
```



### Chicken2:



```r
ps_rare %>%
  subset_samples(Model2 == "Chicken2"  & Treatment != "DONOR" ) -> ps_tmp_forbdiv
```

#### PCoA:


```r
ps_tmp_forbdiv %>%
  # microbiome::transform(transform = "compositional") %>% 
  phyloseq_plot_bdiv(dlist = dist,
                     m = "PCoA",
                     seed = 123,
                     axis1 = 1,
                     axis2 = 2) -> PCoA_tmp
```

```
## [1] "bjaccard"
## [1] "aichinson"
```

#### aichinson


```r
PCoA_tmp$aichinson$layers[[1]] = NULL; PCoA_tmp$aichinson$layers[[1]] = NULL
PCoA_tmp$aichinson$layers[[2]] = NULL; PCoA_tmp$aichinson$layers[[1]] = NULL

# plots_hall_humans$aichinson$layers[[1]] = NULL;plots_hall_humans$aichinson$layers[[1]] = NULL
# plots_hall_humans$aichinson$layers[[2]] = NULL;plots_hall_humans$aichinson$layers[[2]] = NULL


PCoA_tmp$aichinson + geom_point(size = 3,
                                aes(colour = Treatment ,
                                    shape = Antibiotic_mg.L,
                                    alpha = NULL)) +
  geom_path(data = PCoA_tmp$aichinson$data %>%
              arrange(Day_of_Treatment) ,
            # aes(colour = Treatment, group = interaction(Model, Model2, Antibiotic, Treatment, Fermentation, Reactor,Antibiotic_mg.L)),
            aes(colour = Treatment, group = interaction(Model2, Reactor, Treatment, Treatment_Dose)),
            
            arrow = arrow(
              angle = 30, length = unit(0.15, "inches"),
              ends = "last", type = "open"
            ), linetype = "longdash", size = 0.1) +
  theme_light() +
  scale_color_manual(name = "", values = treat_col,
                     na.value = "black") +
  scale_fill_manual(name = "", values = treat_col,
                    na.value = "black") +
  scale_shape_manual(name = "" ,values = c(16,19), na.value =  17) +
  theme(legend.position = "right") +
  facet_grid(. ~ Model2, scales = "free") -> PCoA_tmp2

PCoA_tmp2
```

![](NRP72_draft_16S_Figs_beta_PCoA_files/figure-html/unnamed-chunk-45-1.png)<!-- -->


```r
PCoA_tmp2 + 
  scale_fill_manual(values = c("transparent")) + 
  scale_color_manual(values = c(rep("transparent", 7))) + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) -> empty_plot_tmp


ps_tmp_forbdiv %>%
  phyloseq_add_taxa_vector_fix(phyloseq = .,
                               dist = dist$aichinson,
                               taxrank_glom = "Family",
                               figure_ord = empty_plot_tmp, 
                               fact = 0.2, pval_cutoff = 0.05,
                               top_r = 8) -> out #top 10 correlated genera

out$signenvfit %>% 
  DT::datatable()
```

```{=html}
<div id="htmlwidget-530a78b6872d546cc16a" style="width:100%;height:auto;" class="datatables html-widget"></div>
<script type="application/json" data-for="htmlwidget-530a78b6872d546cc16a">{"x":{"filter":"none","vertical":false,"data":[["ASV0017","ASV0088","ASV0007","ASV0021","ASV0026","ASV0037","ASV0101","ASV0188"],[0.121951517572078,0.3903789246638,-0.574050010181196,0.223621804365789,0.795485555849178,-0.811759510728424,0.524826030206649,0.71786257852125],[-0.784334376671912,-0.727432738348284,0.640754980748867,-0.581167163295317,0.324605021176264,0.42486783158517,0.690527099440638,0.294577706525419],[0.001,0.001,0.001,0.001,0.001,0.001,0.001,0.001],[0.630052587067449,0.681554093642548,0.740100359543512,0.387761983080536,0.738165689337519,0.839466177573934,0.752270037044371,0.602102706822954],["Bacteria","Bacteria","Bacteria","Bacteria","Bacteria","Bacteria","Bacteria","Bacteria"],["Bacteroidota","Firmicutes","Firmicutes","Firmicutes","Firmicutes","Proteobacteria","Firmicutes","Firmicutes"],["Bacteroidia","Clostridia","Bacilli","Bacilli","Clostridia","Gammaproteobacteria","Clostridia","Clostridia"],["Bacteroidales","Oscillospirales","Lactobacillales","Bacillales","Oscillospirales","Burkholderiales","Lachnospirales","Monoglobales"],["Bacteroidaceae","Ruminococcaceae","Lactobacillaceae","Bacillaceae","Oscillospiraceae","Sutterellaceae","Defluviitaleaceae","Monoglobaceae"],[null,null,null,null,null,null,null,null],[null,null,null,null,null,null,null,null],[null,null,null,null,null,null,null,null],[0.00294117647058824,0.00294117647058824,0.00294117647058824,0.00294117647058824,0.00294117647058824,0.00294117647058824,0.00294117647058824,0.00294117647058824]],"container":"<table class=\"display\">\n  <thead>\n    <tr>\n      <th> <\/th>\n      <th>Axis.1<\/th>\n      <th>Axis.2<\/th>\n      <th>pval<\/th>\n      <th>r<\/th>\n      <th>Kingdom<\/th>\n      <th>Phylum<\/th>\n      <th>Class<\/th>\n      <th>Order<\/th>\n      <th>tax_rank_plot<\/th>\n      <th>Genus<\/th>\n      <th>Species<\/th>\n      <th>Strain<\/th>\n      <th>pval.adj<\/th>\n    <\/tr>\n  <\/thead>\n<\/table>","options":{"columnDefs":[{"className":"dt-right","targets":[1,2,3,4,13]},{"orderable":false,"targets":0}],"order":[],"autoWidth":false,"orderClasses":false}},"evals":[],"jsHooks":[]}</script>
```


```r
out$plot -> PCoA_tmp2_fam

PCoA_tmp2 + 
  theme(legend.position = "none") -> PCoA_tmp2_noleg

PCoA_tmp2_noleg
```

![](NRP72_draft_16S_Figs_beta_PCoA_files/figure-html/unnamed-chunk-47-1.png)<!-- -->

```r
PCoA_tmp2  %>%  ggpubr::get_legend() %>% ggpubr::as_ggplot() -> PCoA_tmp2_leg

PCoA_tmp2_noleg %>% 
  export::graph2ppt(append = TRUE,
                    width = 317.48031496 * 1,
                    height = 0.618 * 317.48031496 * 1 , paper = "A4",  scaling = 2,
                    file = out_pptx)

PCoA_tmp2_leg %>% 
  export::graph2ppt(append = TRUE,
                    width = 317.48031496 * 1,
                    height = 0.618 * 317.48031496 * 1 , paper = "A4",  scaling = 2,
                    file = out_pptx)

PCoA_tmp2_fam %>% 
  export::graph2ppt(append = TRUE,
                    width = 317.48031496 * 1,
                    height = 0.618 * 317.48031496 * 1 , paper = "A4",  scaling = 2,
                    file = out_pptx)
```


```r
ps_tmp_forbdiv %>% 
  microViz::ps_mutate(Lactose_mM = as.double(Lactose_mM), 
                      Glucose_mM = as.double(Glucose_mM),
                      Galactose_mM = as.double(Galactose_mM),
                      Succinat_mM = as.double(Succinat_mM),
                      Lactat_mM = as.double(Lactat_mM),
                      Formiat_mM = as.double(Formiat_mM),
                      Acetat_mM = as.double(Acetat_mM),
                      Propionat_mM = as.double(Propionat_mM),
                      Butyrat_mM = as.double(Butyrat_mM),
                      Valerat_mM = as.double(Valerat_mM),
                      Isobutyrat_mM = as.double(Isobutyrat_mM),
                      Isovalerat_mM = as.double(Isovalerat_mM),
                      Total_SCFA_mM = as.double(Total_SCFA_mM))  %>% # sample_data() %>% data.frame()
  subset_samples(Total_SCFA_mM > 0) %>% 
  phyloseq_add_metadata_vector(dist = dist$aichinson,
                               phyloseq = .,
                               figure_ord = empty_plot_tmp,
                               metadata_sel = c("Succinat_mM", "Formiat_mM",  "Acetat_mM" ,
                                                "Propionat_mM", "Butyrat_mM",  "Valerat_mM", "Total_SCFA_mM"),
                               top_r = 6, fact = 0.3, na.rm = TRUE,
                               pval_cutoff = 0.05) -> env_out


env_out$signenvfit %>% 
  DT::datatable()
```

```{=html}
<div id="htmlwidget-4cadc5f8c7ab9ff923ac" style="width:100%;height:auto;" class="datatables html-widget"></div>
<script type="application/json" data-for="htmlwidget-4cadc5f8c7ab9ff923ac">{"x":{"filter":"none","vertical":false,"data":[["1","2","3","4","5","6"],["Succinat_mM","Acetat_mM","Propionat_mM","Butyrat_mM","Valerat_mM","Total_SCFA_mM"],[-0.529147382478684,0.81213679042208,-0.691063349217276,0.850938508204205,0.816344909441183,0.838127471001857],[0.28729710568181,-0.173455788507577,-0.265058981587657,-0.15452103634438,0.327477365938189,-0.0867683544108123],["Succinat_mM","Acetat_mM","Propionat_mM","Butyrat_mM","Valerat_mM","Total_SCFA_mM"],[0.001,0.001,0.001,0.001,0.001,0.001],[0.362536579317188,0.689653076923863,0.547824816351684,0.747973095417739,0.773660436372348,0.70998640497513]],"container":"<table class=\"display\">\n  <thead>\n    <tr>\n      <th> <\/th>\n      <th>label<\/th>\n      <th>Axis.1<\/th>\n      <th>Axis.2<\/th>\n      <th>id<\/th>\n      <th>pval<\/th>\n      <th>r<\/th>\n    <\/tr>\n  <\/thead>\n<\/table>","options":{"columnDefs":[{"className":"dt-right","targets":[2,3,5,6]},{"orderable":false,"targets":0}],"order":[],"autoWidth":false,"orderClasses":false}},"evals":[],"jsHooks":[]}</script>
```


```r
env_out$plot %>% 
  export::graph2ppt(append = TRUE,
                    width = 317.48031496 * 1,
                    height = 0.618 * 317.48031496 * 1 , paper = "A4",  scaling = 2,
                    file = out_pptx)
```

#### bjaccard


```r
PCoA_tmp$bjaccard$layers[[1]] = NULL; PCoA_tmp$bjaccard$layers[[1]] = NULL
PCoA_tmp$bjaccard$layers[[2]] = NULL; PCoA_tmp$bjaccard$layers[[1]] = NULL

# plots_hall_humans$aichinson$layers[[1]] = NULL;plots_hall_humans$aichinson$layers[[1]] = NULL
# plots_hall_humans$aichinson$layers[[2]] = NULL;plots_hall_humans$aichinson$layers[[2]] = NULL


PCoA_tmp$bjaccard + geom_point(size = 3,
                               aes(colour = Treatment ,
                                   shape = Antibiotic_mg.L,
                                   alpha = NULL)) +
  geom_path(data = PCoA_tmp$bjaccard$data %>%
              arrange(Day_of_Treatment) ,
            # aes(colour = Treatment, group = interaction(Model, Model2, Antibiotic, Treatment, Fermentation, Reactor,Antibiotic_mg.L)),
            aes(colour = Treatment, group = interaction(Model2, Reactor, Treatment, Treatment_Dose)),
            
            arrow = arrow(
              angle = 30, length = unit(0.15, "inches"),
              ends = "last", type = "open"
            ), linetype = "longdash", size = 0.1) +
  theme_light() +
  scale_color_manual(name = "", values = treat_col,
                     na.value = "black") +
  scale_fill_manual(name = "", values = treat_col,
                    na.value = "black") +
  scale_shape_manual(name = "" ,values = c(16,19), na.value =  17) +
  theme(legend.position = "right") +
  facet_grid(. ~ Model2, scales = "free") -> PCoA_tmp2

PCoA_tmp2
```

![](NRP72_draft_16S_Figs_beta_PCoA_files/figure-html/unnamed-chunk-50-1.png)<!-- -->


```r
PCoA_tmp2 + 
  scale_fill_manual(values = c("transparent")) + 
  scale_color_manual(values = c(rep("transparent", 7))) + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) -> empty_plot_tmp


ps_tmp_forbdiv %>%
  phyloseq_add_taxa_vector_fix(phyloseq = .,
                               dist = dist$bjaccard,
                               taxrank_glom = "Family",
                               figure_ord = empty_plot_tmp, 
                               fact = 0.2, pval_cutoff = 0.05,
                               top_r = 8) -> out #top 10 correlated genera

out$signenvfit %>% 
  DT::datatable()
```

```{=html}
<div id="htmlwidget-fc8df40fd67c2027a709" style="width:100%;height:auto;" class="datatables html-widget"></div>
<script type="application/json" data-for="htmlwidget-fc8df40fd67c2027a709">{"x":{"filter":"none","vertical":false,"data":[["ASV0017","ASV0088","ASV0007","ASV0021","ASV0026","ASV0037","ASV0101","ASV0188"],[-0.155893725671808,-0.413842315029179,0.607678896355686,-0.233359443007561,-0.623494880289104,0.861704456168324,-0.340307642417437,-0.607019222277422],[-0.651788326623972,-0.717053733649682,0.628652573232729,-0.608858337143351,0.606848572086881,0.356073867053972,0.816685999540255,0.446736709942202],[0.001,0.001,0.001,0.001,0.001,0.001,0.001,0.001],[0.449130876427114,0.685431518649659,0.764477698908196,0.425165104349765,0.75701105519061,0.869323168579117,0.782785313332779,0.56804602422427],["Bacteria","Bacteria","Bacteria","Bacteria","Bacteria","Bacteria","Bacteria","Bacteria"],["Bacteroidota","Firmicutes","Firmicutes","Firmicutes","Firmicutes","Proteobacteria","Firmicutes","Firmicutes"],["Bacteroidia","Clostridia","Bacilli","Bacilli","Clostridia","Gammaproteobacteria","Clostridia","Clostridia"],["Bacteroidales","Oscillospirales","Lactobacillales","Bacillales","Oscillospirales","Burkholderiales","Lachnospirales","Monoglobales"],["Bacteroidaceae","Ruminococcaceae","Lactobacillaceae","Bacillaceae","Oscillospiraceae","Sutterellaceae","Defluviitaleaceae","Monoglobaceae"],[null,null,null,null,null,null,null,null],[null,null,null,null,null,null,null,null],[null,null,null,null,null,null,null,null],[0.00294117647058824,0.00294117647058824,0.00294117647058824,0.00294117647058824,0.00294117647058824,0.00294117647058824,0.00294117647058824,0.00294117647058824]],"container":"<table class=\"display\">\n  <thead>\n    <tr>\n      <th> <\/th>\n      <th>Axis.1<\/th>\n      <th>Axis.2<\/th>\n      <th>pval<\/th>\n      <th>r<\/th>\n      <th>Kingdom<\/th>\n      <th>Phylum<\/th>\n      <th>Class<\/th>\n      <th>Order<\/th>\n      <th>tax_rank_plot<\/th>\n      <th>Genus<\/th>\n      <th>Species<\/th>\n      <th>Strain<\/th>\n      <th>pval.adj<\/th>\n    <\/tr>\n  <\/thead>\n<\/table>","options":{"columnDefs":[{"className":"dt-right","targets":[1,2,3,4,13]},{"orderable":false,"targets":0}],"order":[],"autoWidth":false,"orderClasses":false}},"evals":[],"jsHooks":[]}</script>
```


```r
out$plot -> PCoA_tmp2_fam

PCoA_tmp2 + 
  theme(legend.position = "none") -> PCoA_tmp2_noleg

PCoA_tmp2_noleg
```

![](NRP72_draft_16S_Figs_beta_PCoA_files/figure-html/unnamed-chunk-52-1.png)<!-- -->

```r
PCoA_tmp2  %>%  ggpubr::get_legend() %>% ggpubr::as_ggplot() -> PCoA_tmp2_leg

PCoA_tmp2_noleg %>% 
  export::graph2ppt(append = TRUE,
                    width = 317.48031496 * 1,
                    height = 0.618 * 317.48031496 * 1 , paper = "A4",  scaling = 2,
                    file = out_pptx)

PCoA_tmp2_leg %>% 
  export::graph2ppt(append = TRUE,
                    width = 317.48031496 * 1,
                    height = 0.618 * 317.48031496 * 1 , paper = "A4",  scaling = 2,
                    file = out_pptx)

PCoA_tmp2_fam %>% 
  export::graph2ppt(append = TRUE,
                    width = 317.48031496 * 1,
                    height = 0.618 * 317.48031496 * 1 , paper = "A4",  scaling = 2,
                    file = out_pptx)
```


```r
ps_tmp_forbdiv %>% 
  microViz::ps_mutate(Lactose_mM = as.double(Lactose_mM), 
                      Glucose_mM = as.double(Glucose_mM),
                      Galactose_mM = as.double(Galactose_mM),
                      Succinat_mM = as.double(Succinat_mM),
                      Lactat_mM = as.double(Lactat_mM),
                      Formiat_mM = as.double(Formiat_mM),
                      Acetat_mM = as.double(Acetat_mM),
                      Propionat_mM = as.double(Propionat_mM),
                      Butyrat_mM = as.double(Butyrat_mM),
                      Valerat_mM = as.double(Valerat_mM),
                      Isobutyrat_mM = as.double(Isobutyrat_mM),
                      Isovalerat_mM = as.double(Isovalerat_mM),
                      Total_SCFA_mM = as.double(Total_SCFA_mM))  %>% # sample_data() %>% data.frame()
  subset_samples(Total_SCFA_mM > 0) %>% 
  phyloseq_add_metadata_vector(dist = dist$bjaccard,
                               phyloseq = .,
                               figure_ord = empty_plot_tmp,
                               metadata_sel = c("Succinat_mM", "Formiat_mM",  "Acetat_mM" ,
                                                "Propionat_mM", "Butyrat_mM",  "Valerat_mM", "Total_SCFA_mM"),
                               top_r = 6, fact = 0.3, na.rm = TRUE,
                               pval_cutoff = 0.05) -> env_out


env_out$signenvfit %>% 
  DT::datatable()
```

```{=html}
<div id="htmlwidget-d8c7d17cf4080686b93c" style="width:100%;height:auto;" class="datatables html-widget"></div>
<script type="application/json" data-for="htmlwidget-d8c7d17cf4080686b93c">{"x":{"filter":"none","vertical":false,"data":[["1","2","3","4","5","6"],["Succinat_mM","Acetat_mM","Propionat_mM","Butyrat_mM","Valerat_mM","Total_SCFA_mM"],[0.585052051276282,-0.821564776006437,0.644134692028062,-0.842693140747185,-0.659182695904131,-0.807880122704603],[0.137835067181004,0.112884348508983,-0.309596044903804,0.19370665347383,0.563301066484276,0.247692902410314],["Succinat_mM","Acetat_mM","Propionat_mM","Butyrat_mM","Valerat_mM","Total_SCFA_mM"],[0.001,0.001,0.001,0.001,0.001,0.001],[0.361284408447377,0.687711557312805,0.510759212494164,0.747653997062385,0.751829918081761,0.714022066565649]],"container":"<table class=\"display\">\n  <thead>\n    <tr>\n      <th> <\/th>\n      <th>label<\/th>\n      <th>Axis.1<\/th>\n      <th>Axis.2<\/th>\n      <th>id<\/th>\n      <th>pval<\/th>\n      <th>r<\/th>\n    <\/tr>\n  <\/thead>\n<\/table>","options":{"columnDefs":[{"className":"dt-right","targets":[2,3,5,6]},{"orderable":false,"targets":0}],"order":[],"autoWidth":false,"orderClasses":false}},"evals":[],"jsHooks":[]}</script>
```


```r
env_out$plot %>% 
  export::graph2ppt(append = TRUE,
                    width = 317.48031496 * 1,
                    height = 0.618 * 317.48031496 * 1 , paper = "A4",  scaling = 2,
                    file = out_pptx)
```


Statistical assessment:

- pret - Day_of_Treatment <= 0:


```r
lapply(
  dist,
  FUN = phyloseq_adonis2,
  physeq = ps_tmp_forbdiv %>%
    subset_samples(Day_of_Treatment <= 0),
  formula = paste0(c("Treatment"), collapse=" + "),
  nrep = 999,
  strata = "none",
  terms_margins = "margin"
) %>%
  bind_rows(.id = "Distance")  %>% 
  mutate_if(is.numeric, round, 3) -> tmp_df


tmp_df %>% 
  DT::datatable()
```

```{=html}
<div id="htmlwidget-fbc9973e359619ea20b8" style="width:100%;height:auto;" class="datatables html-widget"></div>
<script type="application/json" data-for="htmlwidget-fbc9973e359619ea20b8">{"x":{"filter":"none","vertical":false,"data":[["1","2","3","4","5","6"],["bjaccard","bjaccard","bjaccard","aichinson","aichinson","aichinson"],["Treatment","Residual","Total","Treatment","Residual","Total"],[6,20,26,6,20,26],[1.139,2.044,3.183,0.671,1.15,1.821],[0.358,0.642,1,0.368,0.632,1],[1.858,null,null,1.945,null,null],[0.001,null,null,0.001,null,null]],"container":"<table class=\"display\">\n  <thead>\n    <tr>\n      <th> <\/th>\n      <th>Distance<\/th>\n      <th>terms<\/th>\n      <th>Df<\/th>\n      <th>SumOfSqs<\/th>\n      <th>R2<\/th>\n      <th>F<\/th>\n      <th>Pr(&gt;F)<\/th>\n    <\/tr>\n  <\/thead>\n<\/table>","options":{"columnDefs":[{"className":"dt-right","targets":[3,4,5,6,7]},{"orderable":false,"targets":0}],"order":[],"autoWidth":false,"orderClasses":false}},"evals":[],"jsHooks":[]}</script>
```

```r
tmp_df %>% 
  xlsx::write.xlsx(file = out_xlsx,
                   append = TRUE,
                   sheetName = "chicken2-beta-pretreat",
                   showNA = TRUE)
```



```r
lapply(
  dist,
  FUN = phyloseq_adonis2,
  physeq = ps_tmp_forbdiv %>%
    subset_samples(Day_of_Treatment <= 0),
  formula = paste0(c("Treatment", "Day_of_Treatment"), collapse=" + "),
  nrep = 999,
  strata = "none",
  terms_margins = "margin"
) %>%
  bind_rows(.id = "Distance")  %>% 
  mutate_if(is.numeric, round, 3) -> tmp_df


tmp_df %>% 
  DT::datatable()
```

```{=html}
<div id="htmlwidget-fc9c72093435cc6a0a18" style="width:100%;height:auto;" class="datatables html-widget"></div>
<script type="application/json" data-for="htmlwidget-fc9c72093435cc6a0a18">{"x":{"filter":"none","vertical":false,"data":[["1","2","3","4","5","6","7","8"],["bjaccard","bjaccard","bjaccard","bjaccard","aichinson","aichinson","aichinson","aichinson"],["Treatment","Day_of_Treatment","Residual","Total","Treatment","Day_of_Treatment","Residual","Total"],[6,1,19,26,6,1,19,26],[1.139,0.533,1.511,3.183,0.671,0.283,0.867,1.821],[0.358,0.167,0.475,1,0.368,0.156,0.476,1],[2.388,6.703,null,null,2.452,6.211,null,null],[0.001,0.001,null,null,0.001,0.001,null,null]],"container":"<table class=\"display\">\n  <thead>\n    <tr>\n      <th> <\/th>\n      <th>Distance<\/th>\n      <th>terms<\/th>\n      <th>Df<\/th>\n      <th>SumOfSqs<\/th>\n      <th>R2<\/th>\n      <th>F<\/th>\n      <th>Pr(&gt;F)<\/th>\n    <\/tr>\n  <\/thead>\n<\/table>","options":{"columnDefs":[{"className":"dt-right","targets":[3,4,5,6,7]},{"orderable":false,"targets":0}],"order":[],"autoWidth":false,"orderClasses":false}},"evals":[],"jsHooks":[]}</script>
```

```r
tmp_df %>% 
  xlsx::write.xlsx(file = out_xlsx,
                   append = TRUE,
                   sheetName = "chicken1-beta-pretreat-time-day",
                   showNA = TRUE)
```


```r
lapply(
  dist,
  FUN = physeq_pairwise_permanovas_adonis2,
  physeq = ps_tmp_forbdiv %>%
    subset_samples(Day_of_Treatment <= 0),
  compare_header = "Treatment",
  n_perm = 999,
  strat = "none"
) %>%
  bind_rows(.id = "Distance") %>%
  mutate_if(is.numeric, round, 3) %>%
  select(-pvalBon)  -> tmp_df
```

```
## Set of permutations < 'minperm'. Generating entire set.
## Set of permutations < 'minperm'. Generating entire set.
## Set of permutations < 'minperm'. Generating entire set.
## Set of permutations < 'minperm'. Generating entire set.
## Set of permutations < 'minperm'. Generating entire set.
## Set of permutations < 'minperm'. Generating entire set.
## Set of permutations < 'minperm'. Generating entire set.
## Set of permutations < 'minperm'. Generating entire set.
## Set of permutations < 'minperm'. Generating entire set.
## Set of permutations < 'minperm'. Generating entire set.
## Set of permutations < 'minperm'. Generating entire set.
## Set of permutations < 'minperm'. Generating entire set.
```

```r
tmp_df %>% 
  DT::datatable()
```

```{=html}
<div id="htmlwidget-33198698d6599fd3fc97" style="width:100%;height:auto;" class="datatables html-widget"></div>
<script type="application/json" data-for="htmlwidget-33198698d6599fd3fc97">{"x":{"filter":"none","vertical":false,"data":[["1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","23","24","25","26","27","28","29","30","31","32","33","34","35","36","37","38","39","40","41","42"],["bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson"],["UNTREATED","UNTREATED","UNTREATED","UNTREATED","UNTREATED","UNTREATED","CTX","CTX","CTX","CTX","CTX","CTX+HV292.1","CTX+HV292.1","CTX+HV292.1","CTX+HV292.1","HV292.1","HV292.1","HV292.1","VAN","VAN","VAN+CCUG59168","UNTREATED","UNTREATED","UNTREATED","UNTREATED","UNTREATED","UNTREATED","CTX","CTX","CTX","CTX","CTX","CTX+HV292.1","CTX+HV292.1","CTX+HV292.1","CTX+HV292.1","HV292.1","HV292.1","HV292.1","VAN","VAN","VAN+CCUG59168"],["CTX","CTX+HV292.1","HV292.1","VAN","VAN+CCUG59168","CCUG59168","CTX+HV292.1","HV292.1","VAN","VAN+CCUG59168","CCUG59168","HV292.1","VAN","VAN+CCUG59168","CCUG59168","VAN","VAN+CCUG59168","CCUG59168","VAN+CCUG59168","CCUG59168","CCUG59168","CTX","CTX+HV292.1","HV292.1","VAN","VAN+CCUG59168","CCUG59168","CTX+HV292.1","HV292.1","VAN","VAN+CCUG59168","CCUG59168","HV292.1","VAN","VAN+CCUG59168","CCUG59168","VAN","VAN+CCUG59168","CCUG59168","VAN+CCUG59168","CCUG59168","CCUG59168"],[0.288,0.279,0.256,0.22,0.225,0.237,0.277,0.33,0.315,0.328,0.297,0.25,0.247,0.222,0.195,0.228,0.199,0.165,0.163,0.212,0.156,0.277,0.27,0.272,0.243,0.243,0.261,0.27,0.322,0.323,0.321,0.308,0.242,0.252,0.221,0.214,0.25,0.194,0.165,0.199,0.232,0.171],[0.03,0.037,0.056,0.081,0.173,0.119,0.066,0.033,0.036,0.028,0.03,0.073,0.021,0.243,0.205,0.033,0.261,0.355,0.42,0.099,0.542,0.03,0.022,0.052,0.053,0.059,0.057,0.02,0.023,0.03,0.028,0.037,0.025,0.026,0.196,0.114,0.03,0.328,0.403,0.313,0.091,0.507],[0.18,0.097,0.131,0.142,0.242,0.178,0.139,0.126,0.108,0.294,0.18,0.139,0.441,0.3,0.269,0.126,0.305,0.392,0.441,0.16,0.542,0.079,0.231,0.099,0.093,0.088,0.092,0.42,0.161,0.079,0.098,0.078,0.131,0.109,0.242,0.15,0.079,0.363,0.423,0.365,0.127,0.507]],"container":"<table class=\"display\">\n  <thead>\n    <tr>\n      <th> <\/th>\n      <th>Distance<\/th>\n      <th>X1<\/th>\n      <th>X2<\/th>\n      <th>R2<\/th>\n      <th>pval<\/th>\n      <th>pvalFDR<\/th>\n    <\/tr>\n  <\/thead>\n<\/table>","options":{"columnDefs":[{"className":"dt-right","targets":[4,5,6]},{"orderable":false,"targets":0}],"order":[],"autoWidth":false,"orderClasses":false}},"evals":[],"jsHooks":[]}</script>
```

```r
tmp_df %>% 
  xlsx::write.xlsx(file = out_xlsx,
                   append = TRUE,
                   sheetName = "chicken2-beta-pretreat-pw",
                   showNA = TRUE)
```


```r
lapply(
  dist,
  FUN = phyloseq_TW,
  physeq =  ps_tmp_forbdiv %>%
    subset_samples(Day_of_Treatment <= 0),
  variable = "Treatment",
  nrep = 999,
  strata = NULL
) %>%
  bind_rows(.id = "Distance") %>%
  mutate_if(is.numeric, round, 3)   -> tmp_df

tmp_df %>% 
  DT::datatable()
```

```{=html}
<div id="htmlwidget-38c86864160d6d855a80" style="width:100%;height:auto;" class="datatables html-widget"></div>
<script type="application/json" data-for="htmlwidget-38c86864160d6d855a80">{"x":{"filter":"none","vertical":false,"data":[["1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","23","24","25","26","27","28","29","30","31","32","33","34","35","36","37","38","39","40","41","42"],["bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson"],[["UNTREATED"],["UNTREATED"],["UNTREATED"],["UNTREATED"],["UNTREATED"],["UNTREATED"],["CTX"],["CTX"],["CTX"],["CTX"],["CTX"],["CTX+HV292.1"],["CTX+HV292.1"],["CTX+HV292.1"],["CTX+HV292.1"],["HV292.1"],["HV292.1"],["HV292.1"],["VAN"],["VAN"],["VAN+CCUG59168"],["UNTREATED"],["UNTREATED"],["UNTREATED"],["UNTREATED"],["UNTREATED"],["UNTREATED"],["CTX"],["CTX"],["CTX"],["CTX"],["CTX"],["CTX+HV292.1"],["CTX+HV292.1"],["CTX+HV292.1"],["CTX+HV292.1"],["HV292.1"],["HV292.1"],["HV292.1"],["VAN"],["VAN"],["VAN+CCUG59168"]],[["CTX"],["CTX+HV292.1"],["HV292.1"],["VAN"],["VAN+CCUG59168"],["CCUG59168"],["CTX+HV292.1"],["HV292.1"],["VAN"],["VAN+CCUG59168"],["CCUG59168"],["HV292.1"],["VAN"],["VAN+CCUG59168"],["CCUG59168"],["VAN"],["VAN+CCUG59168"],["CCUG59168"],["VAN+CCUG59168"],["CCUG59168"],["CCUG59168"],["CTX"],["CTX+HV292.1"],["HV292.1"],["VAN"],["VAN+CCUG59168"],["CCUG59168"],["CTX+HV292.1"],["HV292.1"],["VAN"],["VAN+CCUG59168"],["CCUG59168"],["HV292.1"],["VAN"],["VAN+CCUG59168"],["CCUG59168"],["VAN"],["VAN+CCUG59168"],["CCUG59168"],["VAN+CCUG59168"],["CCUG59168"],["CCUG59168"]],[[4],[4],[4],[4],[4],[4],[4],[4],[4],[4],[4],[4],[4],[4],[4],[4],[4],[4],[4],[4],[3],[4],[4],[4],[4],[4],[4],[4],[4],[4],[4],[4],[4],[4],[4],[4],[4],[4],[4],[4],[4],[3]],[[4],[4],[4],[4],[3],[4],[4],[4],[4],[3],[4],[4],[4],[3],[4],[4],[3],[4],[3],[4],[4],[4],[4],[4],[4],[3],[4],[4],[4],[4],[3],[4],[4],[4],[3],[4],[4],[3],[4],[3],[4],[4]],[[0.036],[0.038],[0.05],[0.092],[0.172],[0.115],[0.074],[0.027],[0.023],[0.027],[0.035],[0.053],[0.027],[0.281],[0.209],[0.028],[0.344],[0.347],[0.453],[0.071],[0.543],[0.02],[0.02],[0.046],[0.05],[0.055],[0.047],[0.025],[0.035],[0.045],[0.028],[0.028],[0.033],[0.041],[0.238],[0.11],[0.028],[0.359],[0.402],[0.338],[0.091],[0.519]],[[2.42768731704433],[2.3176435349535],[2.06407963996001],[1.69044480882103],[1.39825719844339],[1.86138594904126],[2.29910472270134],[2.95195634874219],[2.76129533895811],[2.30248435121415],[2.53744553116146],[2.00080306668269],[1.96454198108462],[1.35512450949592],[1.45337740732689],[1.77480560087452],[1.21224402205629],[1.18315467296777],[0.967554646114228],[1.61886672977221],[0.928915040641932],[2.30203746944132],[2.2178399491206],[2.24197435718492],[1.9262637717644],[1.58697154231561],[2.11841504360238],[2.21973341570745],[2.84820488698467],[2.85743331810361],[2.27527342328822],[2.6648408814652],[1.91417919887335],[2.0206529795731],[1.38312048966709],[1.63038580849202],[1.99984960103855],[1.19246259968342],[1.18282929013113],[1.21301488060513],[1.80861943700909],[1.02809403368817]],[[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999]],[0.1,0.1,0.111,0.149,0.241,0.173,0.13,0.1,0.1,0.1,0.1,0.111,0.1,0.347,0.274,0.1,0.384,0.384,0.476,0.13,0.543,0.081,0.081,0.081,0.081,0.082,0.081,0.081,0.081,0.081,0.081,0.081,0.081,0.081,0.294,0.144,0.081,0.397,0.422,0.394,0.127,0.519]],"container":"<table class=\"display\">\n  <thead>\n    <tr>\n      <th> <\/th>\n      <th>Distance<\/th>\n      <th>Level1<\/th>\n      <th>Level2<\/th>\n      <th>N1<\/th>\n      <th>N2<\/th>\n      <th>p.value<\/th>\n      <th>tw2.stat<\/th>\n      <th>nrep<\/th>\n      <th>pvalFDR<\/th>\n    <\/tr>\n  <\/thead>\n<\/table>","options":{"columnDefs":[{"className":"dt-right","targets":9},{"orderable":false,"targets":0}],"order":[],"autoWidth":false,"orderClasses":false}},"evals":[],"jsHooks":[]}</script>
```

```r
tmp_df %>% 
  xlsx::write.xlsx(file = out_xlsx,
                   append = TRUE,
                   sheetName = "chicken2-beta-pretreat-pw2",
                   showNA = TRUE)
```



```r
# lapply(
#   dlist,
#   FUN = betadisper,
#   physeq = physeq,
#   variable = "Health_status"
# ) %>%
#   bind_rows() %>%
#   mutate_if(is.numeric, round, 3) %>%
#   t() %>%
#   data.frame() %>%
#   dplyr::rename(bdisper_pvalue = '.') %>%
#   DT::datatable()
```



- treatment - Day_of_Treatment >= 1



```r
lapply(
  dist,
  FUN = phyloseq_adonis2,
  physeq = ps_tmp_forbdiv %>%
    subset_samples(Day_of_Treatment >= 1),
  formula = paste0(c("Treatment"), collapse=" + "),
  nrep = 999,
  strata = "none",
  terms_margins = "margin"
) %>%
  bind_rows(.id = "Distance")  %>% 
  mutate_if(is.numeric, round, 3) -> tmp_df


tmp_df %>% 
  DT::datatable()
```

```{=html}
<div id="htmlwidget-d96970b9501c191ad406" style="width:100%;height:auto;" class="datatables html-widget"></div>
<script type="application/json" data-for="htmlwidget-d96970b9501c191ad406">{"x":{"filter":"none","vertical":false,"data":[["1","2","3","4","5","6"],["bjaccard","bjaccard","bjaccard","aichinson","aichinson","aichinson"],["Treatment","Residual","Total","Treatment","Residual","Total"],[6,89,95,6,89,95],[11.193,10.13,21.323,5.417,4.135,9.553],[0.525,0.475,1,0.567,0.433,1],[16.391,null,null,19.432,null,null],[0.001,null,null,0.001,null,null]],"container":"<table class=\"display\">\n  <thead>\n    <tr>\n      <th> <\/th>\n      <th>Distance<\/th>\n      <th>terms<\/th>\n      <th>Df<\/th>\n      <th>SumOfSqs<\/th>\n      <th>R2<\/th>\n      <th>F<\/th>\n      <th>Pr(&gt;F)<\/th>\n    <\/tr>\n  <\/thead>\n<\/table>","options":{"columnDefs":[{"className":"dt-right","targets":[3,4,5,6,7]},{"orderable":false,"targets":0}],"order":[],"autoWidth":false,"orderClasses":false}},"evals":[],"jsHooks":[]}</script>
```

```r
tmp_df %>% 
  xlsx::write.xlsx(file = out_xlsx,
                   append = TRUE,
                   sheetName = "chicken2-beta-treat",
                   showNA = TRUE)
```




```r
lapply(
  dist,
  FUN = phyloseq_adonis2,
  physeq = ps_tmp_forbdiv %>%
    subset_samples(Day_of_Treatment >= 1),
  formula = paste0(c("Treatment", "Day_of_Treatment"), collapse=" + "),
  nrep = 999,
  strata = "none",
  terms_margins = "margin"
) %>%
  bind_rows(.id = "Distance")  %>% 
  mutate_if(is.numeric, round, 3) -> tmp_df


tmp_df %>% 
  DT::datatable()
```

```{=html}
<div id="htmlwidget-bbf49928d385f1dd9e5c" style="width:100%;height:auto;" class="datatables html-widget"></div>
<script type="application/json" data-for="htmlwidget-bbf49928d385f1dd9e5c">{"x":{"filter":"none","vertical":false,"data":[["1","2","3","4","5","6","7","8"],["bjaccard","bjaccard","bjaccard","bjaccard","aichinson","aichinson","aichinson","aichinson"],["Treatment","Day_of_Treatment","Residual","Total","Treatment","Day_of_Treatment","Residual","Total"],[6,1,88,95,6,1,88,95],[11.189,0.701,9.429,21.323,5.415,0.359,3.776,9.553],[0.525,0.033,0.442,1,0.567,0.038,0.395,1],[17.405,6.542,null,null,21.032,8.364,null,null],[0.001,0.001,null,null,0.001,0.001,null,null]],"container":"<table class=\"display\">\n  <thead>\n    <tr>\n      <th> <\/th>\n      <th>Distance<\/th>\n      <th>terms<\/th>\n      <th>Df<\/th>\n      <th>SumOfSqs<\/th>\n      <th>R2<\/th>\n      <th>F<\/th>\n      <th>Pr(&gt;F)<\/th>\n    <\/tr>\n  <\/thead>\n<\/table>","options":{"columnDefs":[{"className":"dt-right","targets":[3,4,5,6,7]},{"orderable":false,"targets":0}],"order":[],"autoWidth":false,"orderClasses":false}},"evals":[],"jsHooks":[]}</script>
```

```r
tmp_df %>% 
  xlsx::write.xlsx(file = out_xlsx,
                   append = TRUE,
                   sheetName = "chicken2-beta-treat-time-treat",
                   showNA = TRUE)
```



```r
lapply(
  dist,
  FUN = physeq_pairwise_permanovas_adonis2,
  physeq = ps_tmp_forbdiv %>%
    subset_samples(Day_of_Treatment >= 1),
  compare_header = "Treatment",
  n_perm = 999,
  strat = "none"
) %>%
  bind_rows(.id = "Distance") %>%
  mutate_if(is.numeric, round, 3) %>%
  select(-pvalBon)  -> tmp_df


tmp_df %>% 
  DT::datatable()
```

```{=html}
<div id="htmlwidget-90bfa43893a2d8eda041" style="width:100%;height:auto;" class="datatables html-widget"></div>
<script type="application/json" data-for="htmlwidget-90bfa43893a2d8eda041">{"x":{"filter":"none","vertical":false,"data":[["1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","23","24","25","26","27","28","29","30","31","32","33","34","35","36","37","38","39","40","41","42"],["bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson"],["UNTREATED","UNTREATED","UNTREATED","UNTREATED","UNTREATED","UNTREATED","CTX","CTX","CTX","CTX","CTX","CTX+HV292.1","CTX+HV292.1","CTX+HV292.1","CTX+HV292.1","HV292.1","HV292.1","HV292.1","VAN","VAN","VAN+CCUG59168","UNTREATED","UNTREATED","UNTREATED","UNTREATED","UNTREATED","UNTREATED","CTX","CTX","CTX","CTX","CTX","CTX+HV292.1","CTX+HV292.1","CTX+HV292.1","CTX+HV292.1","HV292.1","HV292.1","HV292.1","VAN","VAN","VAN+CCUG59168"],["CTX","CTX+HV292.1","HV292.1","VAN","VAN+CCUG59168","CCUG59168","CTX+HV292.1","HV292.1","VAN","VAN+CCUG59168","CCUG59168","HV292.1","VAN","VAN+CCUG59168","CCUG59168","VAN","VAN+CCUG59168","CCUG59168","VAN+CCUG59168","CCUG59168","CCUG59168","CTX","CTX+HV292.1","HV292.1","VAN","VAN+CCUG59168","CCUG59168","CTX+HV292.1","HV292.1","VAN","VAN+CCUG59168","CCUG59168","HV292.1","VAN","VAN+CCUG59168","CCUG59168","VAN","VAN+CCUG59168","CCUG59168","VAN+CCUG59168","CCUG59168","CCUG59168"],[0.321,0.168,0.177,0.48,0.522,0.124,0.284,0.326,0.456,0.485,0.284,0.17,0.476,0.508,0.115,0.468,0.503,0.104,0.085,0.448,0.487,0.367,0.187,0.184,0.587,0.584,0.123,0.345,0.389,0.541,0.536,0.341,0.187,0.586,0.573,0.126,0.589,0.581,0.115,0.116,0.552,0.546],[0.001,0.001,0.001,0.001,0.001,0.001,0.001,0.001,0.001,0.001,0.001,0.001,0.001,0.001,0.001,0.001,0.001,0.005,0.005,0.001,0.001,0.001,0.001,0.001,0.001,0.001,0.001,0.001,0.001,0.001,0.001,0.001,0.001,0.001,0.001,0.001,0.001,0.001,0.001,0.003,0.001,0.001],[0.002,0.002,0.002,0.002,0.002,0.002,0.002,0.002,0.002,0.002,0.002,0.002,0.002,0.002,0.002,0.002,0.002,0.005,0.005,0.002,0.002,0.002,0.002,0.002,0.002,0.002,0.002,0.002,0.002,0.002,0.002,0.002,0.002,0.002,0.002,0.002,0.002,0.002,0.002,0.003,0.002,0.002]],"container":"<table class=\"display\">\n  <thead>\n    <tr>\n      <th> <\/th>\n      <th>Distance<\/th>\n      <th>X1<\/th>\n      <th>X2<\/th>\n      <th>R2<\/th>\n      <th>pval<\/th>\n      <th>pvalFDR<\/th>\n    <\/tr>\n  <\/thead>\n<\/table>","options":{"columnDefs":[{"className":"dt-right","targets":[4,5,6]},{"orderable":false,"targets":0}],"order":[],"autoWidth":false,"orderClasses":false}},"evals":[],"jsHooks":[]}</script>
```

```r
tmp_df %>% 
  xlsx::write.xlsx(file = out_xlsx,
                   append = TRUE,
                   sheetName = "chicken2-beta-treat-pw",
                   showNA = TRUE)
```


```r
lapply(
  dist,
  FUN = phyloseq_TW,
  physeq =  ps_tmp_forbdiv %>%
    subset_samples(Day_of_Treatment >= 1),
  variable = "Treatment",
  nrep = 999,
  strata = NULL
) %>%
  bind_rows(.id = "Distance") %>%
  mutate_if(is.numeric, round, 3) -> tmp_df

tmp_df %>% 
  DT::datatable()
```

```{=html}
<div id="htmlwidget-0349e08474188736d854" style="width:100%;height:auto;" class="datatables html-widget"></div>
<script type="application/json" data-for="htmlwidget-0349e08474188736d854">{"x":{"filter":"none","vertical":false,"data":[["1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","23","24","25","26","27","28","29","30","31","32","33","34","35","36","37","38","39","40","41","42"],["bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson"],[["UNTREATED"],["UNTREATED"],["UNTREATED"],["UNTREATED"],["UNTREATED"],["UNTREATED"],["CTX"],["CTX"],["CTX"],["CTX"],["CTX"],["CTX+HV292.1"],["CTX+HV292.1"],["CTX+HV292.1"],["CTX+HV292.1"],["HV292.1"],["HV292.1"],["HV292.1"],["VAN"],["VAN"],["VAN+CCUG59168"],["UNTREATED"],["UNTREATED"],["UNTREATED"],["UNTREATED"],["UNTREATED"],["UNTREATED"],["CTX"],["CTX"],["CTX"],["CTX"],["CTX"],["CTX+HV292.1"],["CTX+HV292.1"],["CTX+HV292.1"],["CTX+HV292.1"],["HV292.1"],["HV292.1"],["HV292.1"],["VAN"],["VAN"],["VAN+CCUG59168"]],[["CTX"],["CTX+HV292.1"],["HV292.1"],["VAN"],["VAN+CCUG59168"],["CCUG59168"],["CTX+HV292.1"],["HV292.1"],["VAN"],["VAN+CCUG59168"],["CCUG59168"],["HV292.1"],["VAN"],["VAN+CCUG59168"],["CCUG59168"],["VAN"],["VAN+CCUG59168"],["CCUG59168"],["VAN+CCUG59168"],["CCUG59168"],["CCUG59168"],["CTX"],["CTX+HV292.1"],["HV292.1"],["VAN"],["VAN+CCUG59168"],["CCUG59168"],["CTX+HV292.1"],["HV292.1"],["VAN"],["VAN+CCUG59168"],["CCUG59168"],["HV292.1"],["VAN"],["VAN+CCUG59168"],["CCUG59168"],["VAN"],["VAN+CCUG59168"],["CCUG59168"],["VAN+CCUG59168"],["CCUG59168"],["CCUG59168"]],[[14],[14],[14],[14],[14],[14],[14],[14],[14],[14],[14],[14],[14],[14],[14],[13],[13],[13],[14],[14],[14],[14],[14],[14],[14],[14],[14],[14],[14],[14],[14],[14],[14],[14],[14],[14],[13],[13],[13],[14],[14],[14]],[[14],[14],[13],[14],[14],[13],[14],[13],[14],[14],[13],[13],[14],[14],[13],[14],[14],[13],[14],[13],[13],[14],[14],[13],[14],[14],[13],[14],[13],[14],[14],[13],[13],[14],[14],[13],[14],[14],[13],[14],[13],[13]],[[0.001],[0.001],[0.001],[0.001],[0.001],[0.001],[0.001],[0.001],[0.001],[0.001],[0.001],[0.001],[0.001],[0.001],[0.001],[0.001],[0.001],[0.001],[0.005],[0.001],[0.001],[0.001],[0.001],[0.001],[0.001],[0.001],[0.001],[0.001],[0.001],[0.001],[0.001],[0.001],[0.001],[0.001],[0.001],[0.001],[0.001],[0.001],[0.002],[0.004],[0.001],[0.001]],[[12.275013860183],[5.26319741545734],[5.38195079734327],[23.9544929732219],[28.4163454767567],[3.51829800877442],[10.3298741208959],[12.2203218863295],[21.799670369467],[24.5341384782912],[9.9570767067841],[5.14126597086098],[23.6130540593819],[26.8325936087373],[3.25261979211152],[22.7457990487514],[26.029091527308],[2.79584583491889],[2.40034491264051],[20.8808133390709],[24.3468595424238],[15.0816107206083],[5.96313536820619],[5.6539612020765],[36.9195843251591],[36.4956816041149],[3.50486700609353],[13.6930522220193],[15.9925081146084],[30.6832216320199],[29.9762646100757],[12.9400473429492],[5.7752334651681],[36.8189886094972],[34.8312938106115],[3.58823560053639],[35.3443733513806],[34.2562378559437],[3.11550101437793],[3.42610906565254],[30.1466391921396],[29.5970459236577]],[[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999]],[0.001,0.001,0.001,0.001,0.001,0.001,0.001,0.001,0.001,0.001,0.001,0.001,0.001,0.001,0.001,0.001,0.001,0.001,0.005,0.001,0.001,0.001,0.001,0.001,0.001,0.001,0.001,0.001,0.001,0.001,0.001,0.001,0.001,0.001,0.001,0.001,0.001,0.001,0.002,0.004,0.001,0.001]],"container":"<table class=\"display\">\n  <thead>\n    <tr>\n      <th> <\/th>\n      <th>Distance<\/th>\n      <th>Level1<\/th>\n      <th>Level2<\/th>\n      <th>N1<\/th>\n      <th>N2<\/th>\n      <th>p.value<\/th>\n      <th>tw2.stat<\/th>\n      <th>nrep<\/th>\n      <th>pvalFDR<\/th>\n    <\/tr>\n  <\/thead>\n<\/table>","options":{"columnDefs":[{"className":"dt-right","targets":9},{"orderable":false,"targets":0}],"order":[],"autoWidth":false,"orderClasses":false}},"evals":[],"jsHooks":[]}</script>
```

```r
tmp_df %>% 
  xlsx::write.xlsx(file = out_xlsx,
                   append = TRUE,
                   sheetName = "chicken2-beta-treat-pw2",
                   showNA = TRUE)
```

- filter(Day_of_Treatment >= 1 & Day_of_Treatment <= 6) %>% 



```r
lapply(
  dist,
  FUN = phyloseq_adonis2,
  physeq = ps_tmp_forbdiv %>%
    subset_samples(Day_of_Treatment >= 1 & Day_of_Treatment <= 6),
  formula = paste0(c("Treatment"), collapse=" + "),
  nrep = 999,
  strata = "none",
  terms_margins = "margin"
) %>%
  bind_rows(.id = "Distance")  %>% 
  mutate_if(is.numeric, round, 3) -> tmp_df


tmp_df %>% 
  DT::datatable()
```

```{=html}
<div id="htmlwidget-2c91a99e85c6bc68af90" style="width:100%;height:auto;" class="datatables html-widget"></div>
<script type="application/json" data-for="htmlwidget-2c91a99e85c6bc68af90">{"x":{"filter":"none","vertical":false,"data":[["1","2","3","4","5","6"],["bjaccard","bjaccard","bjaccard","aichinson","aichinson","aichinson"],["Treatment","Residual","Total","Treatment","Residual","Total"],[6,34,40,6,34,40],[4.92,3.844,8.764,2.425,1.444,3.869],[0.561,0.439,1,0.627,0.373,1],[7.252,null,null,9.515,null,null],[0.001,null,null,0.001,null,null]],"container":"<table class=\"display\">\n  <thead>\n    <tr>\n      <th> <\/th>\n      <th>Distance<\/th>\n      <th>terms<\/th>\n      <th>Df<\/th>\n      <th>SumOfSqs<\/th>\n      <th>R2<\/th>\n      <th>F<\/th>\n      <th>Pr(&gt;F)<\/th>\n    <\/tr>\n  <\/thead>\n<\/table>","options":{"columnDefs":[{"className":"dt-right","targets":[3,4,5,6,7]},{"orderable":false,"targets":0}],"order":[],"autoWidth":false,"orderClasses":false}},"evals":[],"jsHooks":[]}</script>
```

```r
tmp_df %>% 
  xlsx::write.xlsx(file = out_xlsx,
                   append = TRUE,
                   sheetName = "chicken2-beta-d1-d6",
                   showNA = TRUE)
```


```r
lapply(
  dist,
  FUN = phyloseq_adonis2,
  physeq = ps_tmp_forbdiv %>%
    subset_samples(Day_of_Treatment >= 1 & Day_of_Treatment <= 6),
  formula = paste0(c("Treatment", "Day_of_Treatment"), collapse=" + "),
  nrep = 999,
  strata = "none",
  terms_margins = "margin"
) %>%
  bind_rows(.id = "Distance")  %>% 
  mutate_if(is.numeric, round, 3) -> tmp_df


tmp_df %>% 
  DT::datatable()
```

```{=html}
<div id="htmlwidget-f4acdbbe3467343ec6d3" style="width:100%;height:auto;" class="datatables html-widget"></div>
<script type="application/json" data-for="htmlwidget-f4acdbbe3467343ec6d3">{"x":{"filter":"none","vertical":false,"data":[["1","2","3","4","5","6","7","8"],["bjaccard","bjaccard","bjaccard","bjaccard","aichinson","aichinson","aichinson","aichinson"],["Treatment","Day_of_Treatment","Residual","Total","Treatment","Day_of_Treatment","Residual","Total"],[6,1,33,40,6,1,33,40],[4.917,0.377,3.467,8.764,2.426,0.145,1.299,3.869],[0.561,0.043,0.396,1,0.627,0.037,0.336,1],[7.8,3.593,null,null,10.268,3.674,null,null],[0.001,0.005,null,null,0.001,0.005,null,null]],"container":"<table class=\"display\">\n  <thead>\n    <tr>\n      <th> <\/th>\n      <th>Distance<\/th>\n      <th>terms<\/th>\n      <th>Df<\/th>\n      <th>SumOfSqs<\/th>\n      <th>R2<\/th>\n      <th>F<\/th>\n      <th>Pr(&gt;F)<\/th>\n    <\/tr>\n  <\/thead>\n<\/table>","options":{"columnDefs":[{"className":"dt-right","targets":[3,4,5,6,7]},{"orderable":false,"targets":0}],"order":[],"autoWidth":false,"orderClasses":false}},"evals":[],"jsHooks":[]}</script>
```

```r
tmp_df %>% 
  xlsx::write.xlsx(file = out_xlsx,
                   append = TRUE,
                   sheetName = "chicken2-beta-d1-d6-time-treat",
                   showNA = TRUE)
```



```r
lapply(
  dist,
  FUN = physeq_pairwise_permanovas_adonis2,
  physeq = ps_tmp_forbdiv %>%
    subset_samples(Day_of_Treatment >= 1 & Day_of_Treatment <= 6),
  compare_header = "Treatment",
  n_perm = 999,
  strat = "none"
) %>%
  bind_rows(.id = "Distance") %>%
  mutate_if(is.numeric, round, 3) %>%
  select(-pvalBon)  -> tmp_df


tmp_df %>% 
  DT::datatable()
```

```{=html}
<div id="htmlwidget-d15c90e001b0818c17b0" style="width:100%;height:auto;" class="datatables html-widget"></div>
<script type="application/json" data-for="htmlwidget-d15c90e001b0818c17b0">{"x":{"filter":"none","vertical":false,"data":[["1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","23","24","25","26","27","28","29","30","31","32","33","34","35","36","37","38","39","40","41","42"],["bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson"],["UNTREATED","UNTREATED","UNTREATED","UNTREATED","UNTREATED","UNTREATED","CTX","CTX","CTX","CTX","CTX","CTX+HV292.1","CTX+HV292.1","CTX+HV292.1","CTX+HV292.1","HV292.1","HV292.1","HV292.1","VAN","VAN","VAN+CCUG59168","UNTREATED","UNTREATED","UNTREATED","UNTREATED","UNTREATED","UNTREATED","CTX","CTX","CTX","CTX","CTX","CTX+HV292.1","CTX+HV292.1","CTX+HV292.1","CTX+HV292.1","HV292.1","HV292.1","HV292.1","VAN","VAN","VAN+CCUG59168"],["CTX","CTX+HV292.1","HV292.1","VAN","VAN+CCUG59168","CCUG59168","CTX+HV292.1","HV292.1","VAN","VAN+CCUG59168","CCUG59168","HV292.1","VAN","VAN+CCUG59168","CCUG59168","VAN","VAN+CCUG59168","CCUG59168","VAN+CCUG59168","CCUG59168","CCUG59168","CTX","CTX+HV292.1","HV292.1","VAN","VAN+CCUG59168","CCUG59168","CTX+HV292.1","HV292.1","VAN","VAN+CCUG59168","CCUG59168","HV292.1","VAN","VAN+CCUG59168","CCUG59168","VAN","VAN+CCUG59168","CCUG59168","VAN+CCUG59168","CCUG59168","CCUG59168"],[0.504,0.267,0.322,0.475,0.538,0.234,0.403,0.491,0.458,0.497,0.446,0.204,0.445,0.485,0.18,0.463,0.51,0.252,0.155,0.444,0.495,0.529,0.315,0.351,0.608,0.616,0.246,0.461,0.543,0.558,0.563,0.49,0.256,0.572,0.561,0.209,0.618,0.611,0.257,0.189,0.565,0.56],[0.003,0.003,0.006,0.003,0.002,0.003,0.004,0.006,0.002,0.004,0.002,0.002,0.002,0.002,0.007,0.004,0.005,0.008,0.027,0.004,0.001,0.005,0.003,0.002,0.003,0.001,0.005,0.002,0.002,0.004,0.002,0.003,0.002,0.003,0.005,0.008,0.001,0.002,0.004,0.007,0.002,0.003],[0.007,0.007,0.007,0.007,0.009,0.007,0.006,0.007,0.009,0.006,0.009,0.009,0.009,0.009,0.008,0.006,0.007,0.008,0.027,0.006,0.021,0.006,0.005,0.007,0.005,0.014,0.006,0.007,0.007,0.005,0.007,0.005,0.007,0.005,0.006,0.008,0.014,0.007,0.005,0.007,0.007,0.005]],"container":"<table class=\"display\">\n  <thead>\n    <tr>\n      <th> <\/th>\n      <th>Distance<\/th>\n      <th>X1<\/th>\n      <th>X2<\/th>\n      <th>R2<\/th>\n      <th>pval<\/th>\n      <th>pvalFDR<\/th>\n    <\/tr>\n  <\/thead>\n<\/table>","options":{"columnDefs":[{"className":"dt-right","targets":[4,5,6]},{"orderable":false,"targets":0}],"order":[],"autoWidth":false,"orderClasses":false}},"evals":[],"jsHooks":[]}</script>
```

```r
tmp_df %>% 
  xlsx::write.xlsx(file = out_xlsx,
                   append = TRUE,
                   sheetName = "chicken2-beta-d1-d6-pw",
                   showNA = TRUE)
```


```r
lapply(
  dist,
  FUN = phyloseq_TW,
  physeq =  ps_tmp_forbdiv %>%
    subset_samples(Day_of_Treatment >= 1 & Day_of_Treatment <= 6),
  variable = "Treatment",
  nrep = 999,
  strata = NULL
) %>%
  bind_rows(.id = "Distance") %>%
  mutate_if(is.numeric, round, 3) -> tmp_df

tmp_df %>% 
  DT::datatable()
```

```{=html}
<div id="htmlwidget-fec2d181a4be4b780ca9" style="width:100%;height:auto;" class="datatables html-widget"></div>
<script type="application/json" data-for="htmlwidget-fec2d181a4be4b780ca9">{"x":{"filter":"none","vertical":false,"data":[["1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","23","24","25","26","27","28","29","30","31","32","33","34","35","36","37","38","39","40","41","42"],["bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson"],[["UNTREATED"],["UNTREATED"],["UNTREATED"],["UNTREATED"],["UNTREATED"],["UNTREATED"],["CTX"],["CTX"],["CTX"],["CTX"],["CTX"],["CTX+HV292.1"],["CTX+HV292.1"],["CTX+HV292.1"],["CTX+HV292.1"],["HV292.1"],["HV292.1"],["HV292.1"],["VAN"],["VAN"],["VAN+CCUG59168"],["UNTREATED"],["UNTREATED"],["UNTREATED"],["UNTREATED"],["UNTREATED"],["UNTREATED"],["CTX"],["CTX"],["CTX"],["CTX"],["CTX"],["CTX+HV292.1"],["CTX+HV292.1"],["CTX+HV292.1"],["CTX+HV292.1"],["HV292.1"],["HV292.1"],["HV292.1"],["VAN"],["VAN"],["VAN+CCUG59168"]],[["CTX"],["CTX+HV292.1"],["HV292.1"],["VAN"],["VAN+CCUG59168"],["CCUG59168"],["CTX+HV292.1"],["HV292.1"],["VAN"],["VAN+CCUG59168"],["CCUG59168"],["HV292.1"],["VAN"],["VAN+CCUG59168"],["CCUG59168"],["VAN"],["VAN+CCUG59168"],["CCUG59168"],["VAN+CCUG59168"],["CCUG59168"],["CCUG59168"],["CTX"],["CTX+HV292.1"],["HV292.1"],["VAN"],["VAN+CCUG59168"],["CCUG59168"],["CTX+HV292.1"],["HV292.1"],["VAN"],["VAN+CCUG59168"],["CCUG59168"],["HV292.1"],["VAN"],["VAN+CCUG59168"],["CCUG59168"],["VAN"],["VAN+CCUG59168"],["CCUG59168"],["VAN+CCUG59168"],["CCUG59168"],["CCUG59168"]],[[6],[6],[6],[6],[6],[6],[6],[6],[6],[6],[6],[6],[6],[6],[6],[5],[5],[5],[6],[6],[6],[6],[6],[6],[6],[6],[6],[6],[6],[6],[6],[6],[6],[6],[6],[6],[5],[5],[5],[6],[6],[6]],[[6],[6],[5],[6],[6],[6],[6],[5],[6],[6],[6],[5],[6],[6],[6],[6],[6],[6],[6],[6],[6],[6],[6],[5],[6],[6],[6],[6],[5],[6],[6],[6],[5],[6],[6],[6],[6],[6],[6],[6],[6],[6]],[[0.002],[0.003],[0.003],[0.005],[0.003],[0.007],[0.004],[0.002],[0.005],[0.008],[0.001],[0.006],[0.006],[0.002],[0.01],[0.006],[0.003],[0.006],[0.033],[0.002],[0.003],[0.002],[0.005],[0.003],[0.006],[0.005],[0.003],[0.006],[0.002],[0.005],[0.005],[0.002],[0.002],[0.002],[0.003],[0.011],[0.003],[0.002],[0.003],[0.006],[0.004],[0.003]],[[10.1719190123566],[3.64362318995762],[4.24674862445597],[9.05843831002861],[11.653268081804],[3.05456055657452],[6.75011159591556],[8.95473873783033],[8.44837736244418],[9.86617797591006],[8.0551494048434],[2.38606837136972],[8.03418987284241],[9.43633344874204],[2.196591057439],[8.68501534094451],[10.3104750050889],[3.07840129937027],[1.83969358820631],[7.99219336258711],[9.80623983956641],[11.2424540547971],[4.59792656955416],[4.81526230188922],[15.4997553262797],[16.0584061906422],[3.25775037495296],[8.55196004565511],[10.8430886576117],[12.6394299355037],[12.8673694914918],[9.59501726897082],[3.15047880393915],[13.3643972364844],[12.7966650488171],[2.64736194763329],[14.4369447833394],[13.9217569289494],[3.11262740459403],[2.33173901497172],[12.9888492823664],[12.707724846863]],[[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999]],[0.006,0.006,0.006,0.007,0.006,0.008,0.007,0.006,0.007,0.009,0.006,0.007,0.007,0.006,0.011,0.007,0.006,0.007,0.033,0.006,0.006,0.005,0.006,0.005,0.006,0.006,0.005,0.006,0.005,0.006,0.006,0.005,0.005,0.005,0.005,0.011,0.005,0.005,0.005,0.006,0.006,0.005]],"container":"<table class=\"display\">\n  <thead>\n    <tr>\n      <th> <\/th>\n      <th>Distance<\/th>\n      <th>Level1<\/th>\n      <th>Level2<\/th>\n      <th>N1<\/th>\n      <th>N2<\/th>\n      <th>p.value<\/th>\n      <th>tw2.stat<\/th>\n      <th>nrep<\/th>\n      <th>pvalFDR<\/th>\n    <\/tr>\n  <\/thead>\n<\/table>","options":{"columnDefs":[{"className":"dt-right","targets":9},{"orderable":false,"targets":0}],"order":[],"autoWidth":false,"orderClasses":false}},"evals":[],"jsHooks":[]}</script>
```

```r
tmp_df %>% 
  xlsx::write.xlsx(file = out_xlsx,
                   append = TRUE,
                   sheetName = "chicken2-beta-d1-d6-pw2",
                   showNA = TRUE)
```


- 7-end



```r
lapply(
  dist,
  FUN = phyloseq_adonis2,
  physeq = ps_tmp_forbdiv %>%
    subset_samples(Day_of_Treatment > 7),
  formula = paste0(c("Treatment"), collapse=" + "),
  nrep = 999,
  strata = "none",
  terms_margins = "margin"
) %>%
  bind_rows(.id = "Distance")  %>% 
  mutate_if(is.numeric, round, 3) -> tmp_df


tmp_df %>% 
  DT::datatable()
```

```{=html}
<div id="htmlwidget-bcc569e9a5fb6f3f4fe2" style="width:100%;height:auto;" class="datatables html-widget"></div>
<script type="application/json" data-for="htmlwidget-bcc569e9a5fb6f3f4fe2">{"x":{"filter":"none","vertical":false,"data":[["1","2","3","4","5","6"],["bjaccard","bjaccard","bjaccard","aichinson","aichinson","aichinson"],["Treatment","Residual","Total","Treatment","Residual","Total"],[6,41,47,6,41,47],[7.134,3.257,10.391,3.339,1.32,4.659],[0.687,0.313,1,0.717,0.283,1],[14.969,null,null,17.286,null,null],[0.001,null,null,0.001,null,null]],"container":"<table class=\"display\">\n  <thead>\n    <tr>\n      <th> <\/th>\n      <th>Distance<\/th>\n      <th>terms<\/th>\n      <th>Df<\/th>\n      <th>SumOfSqs<\/th>\n      <th>R2<\/th>\n      <th>F<\/th>\n      <th>Pr(&gt;F)<\/th>\n    <\/tr>\n  <\/thead>\n<\/table>","options":{"columnDefs":[{"className":"dt-right","targets":[3,4,5,6,7]},{"orderable":false,"targets":0}],"order":[],"autoWidth":false,"orderClasses":false}},"evals":[],"jsHooks":[]}</script>
```

```r
tmp_df %>% 
  xlsx::write.xlsx(file = out_xlsx,
                   append = TRUE,
                   sheetName = "chicken2-beta-d7-end",
                   showNA = TRUE)
```


```r
lapply(
  dist,
  FUN = phyloseq_adonis2,
  physeq = ps_tmp_forbdiv %>%
    subset_samples(Day_of_Treatment > 7),
  formula = paste0(c("Treatment", "Day_of_Treatment"), collapse=" + "),
  nrep = 999,
  strata = "none",
  terms_margins = "margin"
) %>%
  bind_rows(.id = "Distance")  %>% 
  mutate_if(is.numeric, round, 3) -> tmp_df


tmp_df %>% 
  DT::datatable()
```

```{=html}
<div id="htmlwidget-23dada22fc5642de2fc1" style="width:100%;height:auto;" class="datatables html-widget"></div>
<script type="application/json" data-for="htmlwidget-23dada22fc5642de2fc1">{"x":{"filter":"none","vertical":false,"data":[["1","2","3","4","5","6","7","8"],["bjaccard","bjaccard","bjaccard","bjaccard","aichinson","aichinson","aichinson","aichinson"],["Treatment","Day_of_Treatment","Residual","Total","Treatment","Day_of_Treatment","Residual","Total"],[6,1,40,47,6,1,40,47],[7.136,0.184,3.073,10.391,3.34,0.103,1.217,4.659],[0.687,0.018,0.296,1,0.717,0.022,0.261,1],[15.48,2.39,null,null,18.303,3.395,null,null],[0.001,0.049,null,null,0.001,0.007,null,null]],"container":"<table class=\"display\">\n  <thead>\n    <tr>\n      <th> <\/th>\n      <th>Distance<\/th>\n      <th>terms<\/th>\n      <th>Df<\/th>\n      <th>SumOfSqs<\/th>\n      <th>R2<\/th>\n      <th>F<\/th>\n      <th>Pr(&gt;F)<\/th>\n    <\/tr>\n  <\/thead>\n<\/table>","options":{"columnDefs":[{"className":"dt-right","targets":[3,4,5,6,7]},{"orderable":false,"targets":0}],"order":[],"autoWidth":false,"orderClasses":false}},"evals":[],"jsHooks":[]}</script>
```

```r
tmp_df %>% 
  xlsx::write.xlsx(file = out_xlsx,
                   append = TRUE,
                   sheetName = "chicken2-beta-d7-end-time-treat",
                   showNA = TRUE)
```



```r
lapply(
  dist,
  FUN = physeq_pairwise_permanovas_adonis2,
  physeq = ps_tmp_forbdiv %>%
    subset_samples(Day_of_Treatment > 7),
  compare_header = "Treatment",
  n_perm = 999,
  strat = "none"
) %>%
  bind_rows(.id = "Distance") %>%
  mutate_if(is.numeric, round, 3) %>%
  select(-pvalBon)  -> tmp_df


tmp_df %>% 
  DT::datatable()
```

```{=html}
<div id="htmlwidget-6e16c2883c26dd1edd48" style="width:100%;height:auto;" class="datatables html-widget"></div>
<script type="application/json" data-for="htmlwidget-6e16c2883c26dd1edd48">{"x":{"filter":"none","vertical":false,"data":[["1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","23","24","25","26","27","28","29","30","31","32","33","34","35","36","37","38","39","40","41","42"],["bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson"],["UNTREATED","UNTREATED","UNTREATED","UNTREATED","UNTREATED","UNTREATED","CTX","CTX","CTX","CTX","CTX","CTX+HV292.1","CTX+HV292.1","CTX+HV292.1","CTX+HV292.1","HV292.1","HV292.1","HV292.1","VAN","VAN","VAN+CCUG59168","UNTREATED","UNTREATED","UNTREATED","UNTREATED","UNTREATED","UNTREATED","CTX","CTX","CTX","CTX","CTX","CTX+HV292.1","CTX+HV292.1","CTX+HV292.1","CTX+HV292.1","HV292.1","HV292.1","HV292.1","VAN","VAN","VAN+CCUG59168"],["CTX","CTX+HV292.1","HV292.1","VAN","VAN+CCUG59168","CCUG59168","CTX+HV292.1","HV292.1","VAN","VAN+CCUG59168","CCUG59168","HV292.1","VAN","VAN+CCUG59168","CCUG59168","VAN","VAN+CCUG59168","CCUG59168","VAN+CCUG59168","CCUG59168","CCUG59168","CTX","CTX+HV292.1","HV292.1","VAN","VAN+CCUG59168","CCUG59168","CTX+HV292.1","HV292.1","VAN","VAN+CCUG59168","CCUG59168","HV292.1","VAN","VAN+CCUG59168","CCUG59168","VAN","VAN+CCUG59168","CCUG59168","VAN+CCUG59168","CCUG59168","CCUG59168"],[0.409,0.263,0.251,0.64,0.682,0.245,0.411,0.407,0.649,0.686,0.434,0.302,0.656,0.698,0.258,0.627,0.667,0.147,0.189,0.63,0.676,0.468,0.268,0.268,0.727,0.722,0.239,0.481,0.487,0.741,0.731,0.49,0.314,0.754,0.744,0.271,0.727,0.723,0.186,0.285,0.742,0.737],[0.002,0.001,0.001,0.002,0.002,0.002,0.002,0.002,0.003,0.001,0.001,0.002,0.001,0.001,0.002,0.002,0.002,0.014,0.001,0.002,0.001,0.001,0.001,0.003,0.005,0.001,0.001,0.001,0.001,0.002,0.001,0.002,0.002,0.002,0.001,0.004,0.001,0.002,0.006,0.001,0.001,0.001],[0.003,0.005,0.005,0.003,0.003,0.003,0.003,0.003,0.003,0.005,0.005,0.003,0.005,0.005,0.003,0.003,0.003,0.014,0.005,0.003,0.005,0.003,0.003,0.004,0.005,0.003,0.003,0.003,0.003,0.003,0.003,0.003,0.003,0.003,0.003,0.004,0.003,0.003,0.006,0.003,0.003,0.003]],"container":"<table class=\"display\">\n  <thead>\n    <tr>\n      <th> <\/th>\n      <th>Distance<\/th>\n      <th>X1<\/th>\n      <th>X2<\/th>\n      <th>R2<\/th>\n      <th>pval<\/th>\n      <th>pvalFDR<\/th>\n    <\/tr>\n  <\/thead>\n<\/table>","options":{"columnDefs":[{"className":"dt-right","targets":[4,5,6]},{"orderable":false,"targets":0}],"order":[],"autoWidth":false,"orderClasses":false}},"evals":[],"jsHooks":[]}</script>
```

```r
tmp_df %>% 
  xlsx::write.xlsx(file = out_xlsx,
                   append = TRUE,
                   sheetName = "chicken2-beta-d7-end-pw",
                   showNA = TRUE)
```


```r
lapply(
  dist,
  FUN = phyloseq_TW,
  physeq =  ps_tmp_forbdiv %>%
    subset_samples(Day_of_Treatment > 7),
  variable = "Treatment",
  nrep = 999,
  strata = NULL
) %>%
  bind_rows(.id = "Distance") %>%
  mutate_if(is.numeric, round, 3) -> tmp_df

tmp_df %>% 
  DT::datatable()
```

```{=html}
<div id="htmlwidget-23e11579ab285b5c7c2a" style="width:100%;height:auto;" class="datatables html-widget"></div>
<script type="application/json" data-for="htmlwidget-23e11579ab285b5c7c2a">{"x":{"filter":"none","vertical":false,"data":[["1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","23","24","25","26","27","28","29","30","31","32","33","34","35","36","37","38","39","40","41","42"],["bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson"],[["UNTREATED"],["UNTREATED"],["UNTREATED"],["UNTREATED"],["UNTREATED"],["UNTREATED"],["CTX"],["CTX"],["CTX"],["CTX"],["CTX"],["CTX+HV292.1"],["CTX+HV292.1"],["CTX+HV292.1"],["CTX+HV292.1"],["HV292.1"],["HV292.1"],["HV292.1"],["VAN"],["VAN"],["VAN+CCUG59168"],["UNTREATED"],["UNTREATED"],["UNTREATED"],["UNTREATED"],["UNTREATED"],["UNTREATED"],["CTX"],["CTX"],["CTX"],["CTX"],["CTX"],["CTX+HV292.1"],["CTX+HV292.1"],["CTX+HV292.1"],["CTX+HV292.1"],["HV292.1"],["HV292.1"],["HV292.1"],["VAN"],["VAN"],["VAN+CCUG59168"]],[["CTX"],["CTX+HV292.1"],["HV292.1"],["VAN"],["VAN+CCUG59168"],["CCUG59168"],["CTX+HV292.1"],["HV292.1"],["VAN"],["VAN+CCUG59168"],["CCUG59168"],["HV292.1"],["VAN"],["VAN+CCUG59168"],["CCUG59168"],["VAN"],["VAN+CCUG59168"],["CCUG59168"],["VAN+CCUG59168"],["CCUG59168"],["CCUG59168"],["CTX"],["CTX+HV292.1"],["HV292.1"],["VAN"],["VAN+CCUG59168"],["CCUG59168"],["CTX+HV292.1"],["HV292.1"],["VAN"],["VAN+CCUG59168"],["CCUG59168"],["HV292.1"],["VAN"],["VAN+CCUG59168"],["CCUG59168"],["VAN"],["VAN+CCUG59168"],["CCUG59168"],["VAN+CCUG59168"],["CCUG59168"],["CCUG59168"]],[[7],[7],[7],[7],[7],[7],[7],[7],[7],[7],[7],[7],[7],[7],[7],[7],[7],[7],[7],[7],[7],[7],[7],[7],[7],[7],[7],[7],[7],[7],[7],[7],[7],[7],[7],[7],[7],[7],[7],[7],[7],[7]],[[7],[7],[7],[7],[7],[6],[7],[7],[7],[7],[6],[7],[7],[7],[6],[7],[7],[6],[7],[6],[6],[7],[7],[7],[7],[7],[6],[7],[7],[7],[7],[6],[7],[7],[7],[6],[7],[7],[6],[7],[6],[6]],[[0.001],[0.002],[0.001],[0.001],[0.002],[0.001],[0.001],[0.001],[0.001],[0.002],[0.001],[0.001],[0.002],[0.001],[0.002],[0.001],[0.002],[0.022],[0.003],[0.001],[0.001],[0.002],[0.001],[0.003],[0.001],[0.001],[0.002],[0.001],[0.001],[0.001],[0.002],[0.001],[0.002],[0.001],[0.001],[0.002],[0.001],[0.003],[0.002],[0.001],[0.001],[0.002]],[[8.31970520375296],[4.28734707422854],[4.02569126040891],[21.3449902843217],[25.6976732520508],[3.57740076541489],[8.37249099212652],[8.25218724070032],[22.2120847293217],[26.2746519110328],[8.474750278158],[5.18857315998673],[22.9266880683215],[27.6793262487296],[3.79341673540311],[20.1524374592084],[24.0615255544875],[1.90983351640976],[2.80408642888525],[19.8432016839777],[23.9285036807009],[10.5640610644187],[4.40166382573709],[4.38453974996799],[31.9622309283787],[31.1197977550124],[3.488750403582],[11.1183181700606],[11.3721611171832],[34.2674000474704],[32.6631745294994],[10.5058975722317],[5.49483078090999],[36.7164447185968],[34.8041035981964],[4.08112906980664],[31.9724772352878],[31.3318116593754],[2.5127271641934],[4.78359450179662],[29.2975848236427],[28.976808561324]],[[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999]],[0.002,0.002,0.002,0.002,0.002,0.002,0.002,0.002,0.002,0.002,0.002,0.002,0.002,0.002,0.002,0.002,0.002,0.022,0.003,0.002,0.002,0.002,0.002,0.003,0.002,0.002,0.002,0.002,0.002,0.002,0.002,0.002,0.002,0.002,0.002,0.002,0.002,0.003,0.002,0.002,0.002,0.002]],"container":"<table class=\"display\">\n  <thead>\n    <tr>\n      <th> <\/th>\n      <th>Distance<\/th>\n      <th>Level1<\/th>\n      <th>Level2<\/th>\n      <th>N1<\/th>\n      <th>N2<\/th>\n      <th>p.value<\/th>\n      <th>tw2.stat<\/th>\n      <th>nrep<\/th>\n      <th>pvalFDR<\/th>\n    <\/tr>\n  <\/thead>\n<\/table>","options":{"columnDefs":[{"className":"dt-right","targets":9},{"orderable":false,"targets":0}],"order":[],"autoWidth":false,"orderClasses":false}},"evals":[],"jsHooks":[]}</script>
```

```r
tmp_df %>% 
  xlsx::write.xlsx(file = out_xlsx,
                   append = TRUE,
                   sheetName = "chicken2-beta-d7-end-pw2",
                   showNA = TRUE)
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
##  [1] vegan_2.6-4        lattice_0.20-45    permute_0.9-7      gdtools_0.2.4     
##  [5] GUniFrac_1.6       ape_5.6-2          reshape2_1.4.4     scales_1.2.1      
##  [9] compositions_2.0-4 nlme_3.1-159       phyloseq_1.40.0    forcats_0.5.2     
## [13] stringr_1.4.1      dplyr_1.0.10       purrr_0.3.5        readr_2.1.3       
## [17] tidyr_1.2.1        tibble_3.1.8       ggplot2_3.3.6      tidyverse_1.3.2   
## 
## loaded via a namespace (and not attached):
##   [1] uuid_1.1-0             readxl_1.4.1           backports_1.4.1       
##   [4] systemfonts_1.0.4      plyr_1.8.7             igraph_1.3.5          
##   [7] splines_4.2.1          crosstalk_1.2.0        GenomeInfoDb_1.32.4   
##  [10] digest_0.6.30          foreach_1.5.2          htmltools_0.5.3       
##  [13] fansi_1.0.3            magrittr_2.0.3         xlsx_0.6.5            
##  [16] googlesheets4_1.0.1    cluster_2.1.4          openxlsx_4.2.5        
##  [19] tzdb_0.3.0             Biostrings_2.64.1      extrafont_0.18        
##  [22] modelr_0.1.9           bayesm_3.1-4           matrixStats_0.62.0    
##  [25] officer_0.4.4          stabledist_0.7-1       extrafontdb_1.0       
##  [28] colorspace_2.0-3       rvest_1.0.3            ggrepel_0.9.1         
##  [31] textshaping_0.3.6      haven_2.5.1            xfun_0.34             
##  [34] crayon_1.5.2           RCurl_1.98-1.9         jsonlite_1.8.3        
##  [37] survival_3.4-0         iterators_1.0.14       glue_1.6.2            
##  [40] rvg_0.2.5              gtable_0.3.1           gargle_1.2.1          
##  [43] zlibbioc_1.42.0        XVector_0.36.0         car_3.1-0             
##  [46] Rttf2pt1_1.3.10        Rhdf5lib_1.18.2        BiocGenerics_0.42.0   
##  [49] DEoptimR_1.0-11        abind_1.4-5            DBI_1.1.3             
##  [52] rstatix_0.7.0          Rcpp_1.0.9             xtable_1.8-4          
##  [55] clue_0.3-62            DT_0.25                stats4_4.2.1          
##  [58] htmlwidgets_1.5.4      timeSeries_4021.104    httr_1.4.4            
##  [61] ellipsis_0.3.2         spatial_7.3-15         rJava_1.0-6           
##  [64] pkgconfig_2.0.3        farver_2.1.1           sass_0.4.2            
##  [67] dbplyr_2.2.1           utf8_1.2.2             here_1.0.1            
##  [70] tidyselect_1.2.0       labeling_0.4.2         rlang_1.0.6           
##  [73] munsell_0.5.0          cellranger_1.1.0       tools_4.2.1           
##  [76] cachem_1.0.6           cli_3.4.1              generics_0.1.3        
##  [79] devEMF_4.1             ade4_1.7-19            export_0.3.0          
##  [82] broom_1.0.1            evaluate_0.17          biomformat_1.24.0     
##  [85] fastmap_1.1.0          yaml_2.3.6             ragg_1.2.3            
##  [88] knitr_1.40             fs_1.5.2               zip_2.2.1             
##  [91] robustbase_0.95-0      rgl_0.110.2            xml2_1.3.3            
##  [94] compiler_4.2.1         rstudioapi_0.14        ggsignif_0.6.4        
##  [97] reprex_2.0.2           statmod_1.4.37         bslib_0.4.0           
## [100] stringi_1.7.8          statip_0.2.3           highr_0.9             
## [103] stargazer_5.2.3        modeest_2.4.0          fBasics_4021.92       
## [106] Matrix_1.5-1           microbiome_1.19.1      tensorA_0.36.2        
## [109] multtest_2.52.0        vctrs_0.4.2            pillar_1.8.1          
## [112] lifecycle_1.0.3        rhdf5filters_1.8.0     microViz_0.9.7        
## [115] jquerylib_0.1.4        flextable_0.8.2        cowplot_1.1.1         
## [118] data.table_1.14.4      bitops_1.0-7           R6_2.5.1              
## [121] stable_1.1.6           IRanges_2.30.1         codetools_0.2-18      
## [124] MASS_7.3-58.1          assertthat_0.2.1       xlsxjars_0.6.1        
## [127] rhdf5_2.40.0           rprojroot_2.0.3        withr_2.5.0           
## [130] S4Vectors_0.34.0       GenomeInfoDbData_1.2.8 mgcv_1.8-40           
## [133] parallel_4.2.1         hms_1.1.2              speedyseq_0.5.3.9018  
## [136] grid_4.2.1             rpart_4.1.16           timeDate_4021.106     
## [139] rmarkdown_2.16         carData_3.0-5          rmutil_1.1.9          
## [142] googledrive_2.0.0      Rtsne_0.16             ggpubr_0.4.0          
## [145] base64enc_0.1-3        Biobase_2.56.0         lubridate_1.8.0
```
