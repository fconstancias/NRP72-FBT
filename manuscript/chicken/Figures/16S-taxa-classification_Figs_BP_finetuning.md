---
title: "NTP72 - 16S - taxa compostiion - chicken draft "
author: "Florentin Constancias"
date: "December 14, 2022"
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
out_pptx = "~/Desktop/16S-chick.pptx"
out_xlsx = "~/Desktop/16S-chick.xlsx"
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
  subset_samples(Day_of_Treatment > -4 & Day_of_Treatment <= 30 | Reactor == "DONOR") %>% 
  subset_samples(Model == "Chicken") %>% 
  subset_samples(Reactor != "CR" & Model2 == "Chicken1" | Model2 == "Chicken2" & Reactor %in% c("TR2", "TR3", "TR4", "TR5", "TR6", "TR7", "CR", "DONOR")) %>% 
  subset_samples(Reactor != "IR") %>% 
  microViz::ps_mutate("Model2_reactor_treatment_dose" = paste0(Model2,"_",Reactor_Treatment_Dose))-> ps_filtered
```


```r
ps_filtered %>% 
  phyloseq_check_lib_size(data_color = "Reactor", data_facet = "Reactor", first_n = 50, nreads_display = 1000) -> out
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
## Warning: `aes_string()` was deprecated in ggplot2 3.0.0.
## ℹ Please use tidy evaluation ideoms with `aes()`
```

```
## Warning: Using `size` aesthetic for lines was deprecated in ggplot2 3.4.0.
## ℹ Please use `linewidth` instead.
```

```r
out$df %>% 
  # rownames_to_column('sample_id') %>% 
  # filter(Treatment == "DONOR")  %>% 
  select(SampleID, Day_of_Treatment,Treatment, Reactor, Model2, Model, Period, LibrarySize) 
```

```
##         SampleID Day_of_Treatment     Treatment Reactor   Model2   Model Period
## 191  p-2-H8-S188               10     CCUG59168     TR4 Chicken2 Chicken     t2
## 6    TR1-16-S199                0   CTX+HV292.1     TR1 Chicken1 Chicken     t1
## 59   TR4-18-S200                2           VAN     TR4 Chicken1 Chicken     t2
## 33   TR2-34-S302               18           CTX     TR2 Chicken1 Chicken     t2
## 48   TR3-28-S162               12       HV292.1     TR3 Chicken1 Chicken     t2
## 5    TR1-15-S168               -1   CTX+HV292.1     TR1 Chicken1 Chicken   <NA>
## 87   TR6-13-S164               -3     CCUG59168     TR6 Chicken1 Chicken   pret
## 197  p-3-A5-S197               -2 VAN+CCUG59168     TR7 Chicken2 Chicken   pret
## 38   TR3-14-S136               -2       HV292.1     TR3 Chicken1 Chicken   pret
## 151 p-2-D11-S143               13     UNTREATED      CR Chicken2 Chicken     t2
## 95   TR6-21-S107                5     CCUG59168     TR6 Chicken1 Chicken     t5
## 215  p-3-G1-S265               14       HV292.1     TR5 Chicken2 Chicken     t2
## 181  p-2-G9-S177               11           CTX     TR3 Chicken2 Chicken     t2
## 194  p-3-A2-S194               14 VAN+CCUG59168     TR7 Chicken2 Chicken     t2
## 51   TR3-37-S221               21       HV292.1     TR3 Chicken1 Chicken     t3
## 9    TR1-19-S305                3   CTX+HV292.1     TR1 Chicken1 Chicken     t3
## 112  p-1-F11-S71                2     UNTREATED      CR Chicken2 Chicken     t2
## 100  TR6-34-S143               18     CCUG59168     TR6 Chicken1 Chicken     t2
## 146  p-2-C7-S127                8 VAN+CCUG59168     TR7 Chicken2 Chicken   <NA>
## 37   TR3-13-S304               -3       HV292.1     TR3 Chicken1 Chicken   pret
## 30   TR2-25-S133                9           CTX     TR2 Chicken1 Chicken   <NA>
## 8    TR1-18-S124                2   CTX+HV292.1     TR1 Chicken1 Chicken     t2
## 26   TR2-19-S318                3           CTX     TR2 Chicken1 Chicken     t3
## 1       D-0-S154             <NA>         DONOR   DONOR Chicken1 Chicken   <NA>
## 204  p-3-D2-S230                0   CTX+HV292.1     TR2 Chicken2 Chicken     t1
## 83   TR5-34-S127               18 VAN+CCUG59168     TR5 Chicken1 Chicken     t2
## 198  p-3-B1-S205               14     UNTREATED      CR Chicken2 Chicken     t2
## 217  p-3-G3-S267               -1       HV292.1     TR5 Chicken2 Chicken   <NA>
## 212  p-3-F2-S254                0     CCUG59168     TR4 Chicken2 Chicken     t1
## 176  p-2-G4-S172                6           CTX     TR3 Chicken2 Chicken   <NA>
## 98   TR6-28-S158               12     CCUG59168     TR6 Chicken1 Chicken     t2
## 208  p-3-E2-S242                0           CTX     TR3 Chicken2 Chicken     t1
## 16   TR1-34-S317               18   CTX+HV292.1     TR1 Chicken1 Chicken     t2
## 205  p-3-D3-S231               -1   CTX+HV292.1     TR2 Chicken2 Chicken   <NA>
## 222  p-3-H4-S280               -2           VAN     TR6 Chicken2 Chicken   pret
## 211  p-3-F1-S253               14     CCUG59168     TR4 Chicken2 Chicken     t2
## 199  p-3-B2-S206                0     UNTREATED      CR Chicken2 Chicken     t1
## 65   TR4-28-S183               12           VAN     TR4 Chicken1 Chicken     t2
## 144  p-2-C5-S125                6 VAN+CCUG59168     TR7 Chicken2 Chicken   <NA>
## 206  p-3-D4-S232               -2   CTX+HV292.1     TR2 Chicken2 Chicken   pret
## 67   TR4-34-S342               18           VAN     TR4 Chicken1 Chicken     t2
## 209  p-3-E3-S243               -1           CTX     TR3 Chicken2 Chicken   <NA>
## 82   TR5-31-S313               15 VAN+CCUG59168     TR5 Chicken1 Chicken     t2
## 41   TR3-17-S222                1       HV292.1     TR3 Chicken1 Chicken     t1
## 53   TR3-46-S224               30       HV292.1     TR3 Chicken1 Chicken     t4
## 102  TR6-40-S160               24     CCUG59168     TR6 Chicken1 Chicken     t3
## 195  p-3-A3-S195                0 VAN+CCUG59168     TR7 Chicken2 Chicken     t1
## 113  p-1-H10-S94                1   CTX+HV292.1     TR2 Chicken2 Chicken     t1
## 64   TR4-25-S333                9           VAN     TR4 Chicken1 Chicken   <NA>
## 128 p-2-B11-S119               12           VAN     TR6 Chicken2 Chicken     t2
## 101  TR6-37-S312               21     CCUG59168     TR6 Chicken1 Chicken     t3
## 28   TR2-21-S228                5           CTX     TR2 Chicken1 Chicken     t5
## 213  p-3-F3-S255               -1     CCUG59168     TR4 Chicken2 Chicken   <NA>
## 185  p-2-H2-S182                4     CCUG59168     TR4 Chicken2 Chicken     t4
## 203  p-3-D1-S229               14   CTX+HV292.1     TR2 Chicken2 Chicken     t2
## 34   TR2-37-S330               21           CTX     TR2 Chicken1 Chicken     t3
## 63   TR4-22-S344                6           VAN     TR4 Chicken1 Chicken   <NA>
## 25   TR2-18-S326                2           CTX     TR2 Chicken1 Chicken     t2
## 216  p-3-G2-S266                0       HV292.1     TR5 Chicken2 Chicken     t1
## 7    TR1-17-S246                1   CTX+HV292.1     TR1 Chicken1 Chicken     t1
## 218  p-3-G4-S268               -2       HV292.1     TR5 Chicken2 Chicken   pret
## 57   TR4-16-S276                0           VAN     TR4 Chicken1 Chicken     t1
## 110  p-1-E11-S59                1 VAN+CCUG59168     TR7 Chicken2 Chicken     t1
## 111  p-1-F10-S70                1     UNTREATED      CR Chicken2 Chicken     t1
## 39   TR3-15-S338               -1       HV292.1     TR3 Chicken1 Chicken   <NA>
## 182  p-2-H1-S181                3     CCUG59168     TR4 Chicken2 Chicken     t3
## 94   TR6-20-S132                4     CCUG59168     TR6 Chicken1 Chicken     t4
## 29   TR2-22-S364                6           CTX     TR2 Chicken1 Chicken   <NA>
## 54   TR4-13-S142               -3           VAN     TR4 Chicken1 Chicken   pret
## 214  p-3-F4-S256               -2     CCUG59168     TR4 Chicken2 Chicken   pret
## 85   TR5-40-S335               24 VAN+CCUG59168     TR5 Chicken1 Chicken     t3
## 114  p-1-H11-S95                2   CTX+HV292.1     TR2 Chicken2 Chicken     t2
## 106  p-1-B11-S23                1     CCUG59168     TR4 Chicken2 Chicken     t1
## 11   TR1-21-S207                5   CTX+HV292.1     TR1 Chicken1 Chicken     t5
## 91   TR6-17-S138                1     CCUG59168     TR6 Chicken1 Chicken     t1
## 50   TR3-34-S354               18       HV292.1     TR3 Chicken1 Chicken     t2
## 220  p-3-H2-S278                0           VAN     TR6 Chicken2 Chicken     t1
## 86   TR5-46-S310               30 VAN+CCUG59168     TR5 Chicken1 Chicken     t4
## 184 p-2-H11-S191               13     CCUG59168     TR4 Chicken2 Chicken     t2
## 188  p-2-H5-S185                7     CCUG59168     TR4 Chicken2 Chicken   <NA>
## 21   TR2-14-S232               -2           CTX     TR2 Chicken1 Chicken   pret
## 61   TR4-20-S269                4           VAN     TR4 Chicken1 Chicken     t4
## 55   TR4-14-S226               -2           VAN     TR4 Chicken1 Chicken   pret
## 44   TR3-20-S265                4       HV292.1     TR3 Chicken1 Chicken     t4
## 207  p-3-E1-S241               14           CTX     TR3 Chicken2 Chicken     t2
## 131  p-2-B3-S111                4           VAN     TR6 Chicken2 Chicken     t4
## 35   TR2-40-S259               24           CTX     TR2 Chicken1 Chicken     t3
## 149  p-2-D1-S133                3     UNTREATED      CR Chicken2 Chicken     t3
## 187  p-2-H4-S184                6     CCUG59168     TR4 Chicken2 Chicken   <NA>
## 196  p-3-A4-S196               -1 VAN+CCUG59168     TR7 Chicken2 Chicken   <NA>
## 17   TR1-37-S306               21   CTX+HV292.1     TR1 Chicken1 Chicken     t3
## 89   TR6-15-S149               -1     CCUG59168     TR6 Chicken1 Chicken   <NA>
## 202  p-3-C7-S223             <NA>         DONOR   DONOR Chicken2 Chicken   <NA>
## 56   TR4-15-S218               -1           VAN     TR4 Chicken1 Chicken   <NA>
## 20   TR2-13-S231               -3           CTX     TR2 Chicken1 Chicken   pret
## 62   TR4-21-S254                5           VAN     TR4 Chicken1 Chicken     t5
## 84   TR5-37-S368               21 VAN+CCUG59168     TR5 Chicken1 Chicken     t3
## 133  p-2-B5-S113                6           VAN     TR6 Chicken2 Chicken   <NA>
## 108  p-1-C11-S35                1       HV292.1     TR5 Chicken2 Chicken     t1
## 221  p-3-H3-S279               -1           VAN     TR6 Chicken2 Chicken   <NA>
## 10   TR1-20-S118                4   CTX+HV292.1     TR1 Chicken1 Chicken     t4
## 45   TR3-21-S234                5       HV292.1     TR3 Chicken1 Chicken     t5
## 80   TR5-25-S240                9 VAN+CCUG59168     TR5 Chicken1 Chicken   <NA>
## 88   TR6-14-S211               -2     CCUG59168     TR6 Chicken1 Chicken   pret
## 134  p-2-B6-S114                7           VAN     TR6 Chicken2 Chicken   <NA>
## 169  p-2-F8-S164               10   CTX+HV292.1     TR2 Chicken2 Chicken     t2
## 68   TR4-37-S188               21           VAN     TR4 Chicken1 Chicken     t3
## 135  p-2-B7-S115                8           VAN     TR6 Chicken2 Chicken   <NA>
## 81   TR5-28-S214               12 VAN+CCUG59168     TR5 Chicken1 Chicken     t2
## 130  p-2-B2-S110                3           VAN     TR6 Chicken2 Chicken     t3
## 77   TR5-20-S208                4 VAN+CCUG59168     TR5 Chicken1 Chicken     t4
## 107  p-1-B12-S24                2     CCUG59168     TR4 Chicken2 Chicken     t2
## 132  p-2-B4-S112                5           VAN     TR6 Chicken2 Chicken     t5
## 138  p-2-C1-S121                2 VAN+CCUG59168     TR7 Chicken2 Chicken     t2
## 13   TR1-25-S117                9   CTX+HV292.1     TR1 Chicken1 Chicken   <NA>
## 168  p-2-F7-S163                9   CTX+HV292.1     TR2 Chicken2 Chicken   <NA>
## 93   TR6-19-S137                3     CCUG59168     TR6 Chicken1 Chicken     t3
## 150 p-2-D10-S142               12     UNTREATED      CR Chicken2 Chicken     t2
## 70   TR4-46-S215               30           VAN     TR4 Chicken1 Chicken     t4
## 90   TR6-16-S371                0     CCUG59168     TR6 Chicken1 Chicken     t1
## 52   TR3-40-S257               24       HV292.1     TR3 Chicken1 Chicken     t3
## 71   TR5-13-S103               -3 VAN+CCUG59168     TR5 Chicken1 Chicken   pret
## 36   TR2-46-S239               30           CTX     TR2 Chicken1 Chicken     t4
## 166  p-2-F5-S161                7   CTX+HV292.1     TR2 Chicken2 Chicken   <NA>
## 32   TR2-31-S182               15           CTX     TR2 Chicken1 Chicken     t2
## 192  p-2-H9-S189               11     CCUG59168     TR4 Chicken2 Chicken     t2
## 99   TR6-31-S272               15     CCUG59168     TR6 Chicken1 Chicken     t2
## 73   TR5-16-S290                0 VAN+CCUG59168     TR5 Chicken1 Chicken     t1
## 79   TR5-22-S145                6 VAN+CCUG59168     TR5 Chicken1 Chicken   <NA>
## 162 p-2-F11-S167               13   CTX+HV292.1     TR2 Chicken2 Chicken     t2
## 47   TR3-25-S102                9       HV292.1     TR3 Chicken1 Chicken   <NA>
## 154  p-2-D4-S136                6     UNTREATED      CR Chicken2 Chicken   <NA>
## 124  p-2-A8-S104                9       HV292.1     TR5 Chicken2 Chicken   <NA>
## 143  p-2-C4-S124                5 VAN+CCUG59168     TR7 Chicken2 Chicken     t5
## 159  p-2-D9-S141               11     UNTREATED      CR Chicken2 Chicken     t2
## 60   TR4-19-S270                3           VAN     TR4 Chicken1 Chicken     t3
## 92   TR6-18-S262                2     CCUG59168     TR6 Chicken1 Chicken     t2
## 105  p-1-A12-S12                2           CTX     TR3 Chicken2 Chicken     t2
## 75   TR5-18-S201                2 VAN+CCUG59168     TR5 Chicken1 Chicken     t2
## 153  p-2-D3-S135                5     UNTREATED      CR Chicken2 Chicken     t5
## 58   TR4-17-S114                1           VAN     TR4 Chicken1 Chicken     t1
## 161 p-2-F10-S166               12   CTX+HV292.1     TR2 Chicken2 Chicken     t2
## 158  p-2-D8-S140               10     UNTREATED      CR Chicken2 Chicken     t2
## 155  p-2-D5-S137                7     UNTREATED      CR Chicken2 Chicken   <NA>
## 219  p-3-H1-S277               14           VAN     TR6 Chicken2 Chicken     t2
## 109  p-1-D11-S47                1           VAN     TR6 Chicken2 Chicken     t1
## 136  p-2-B8-S116                9           VAN     TR6 Chicken2 Chicken   <NA>
## 193  p-3-A1-S193               13 VAN+CCUG59168     TR7 Chicken2 Chicken     t2
## 142  p-2-C3-S123                4 VAN+CCUG59168     TR7 Chicken2 Chicken     t4
## 18   TR1-40-S209               24   CTX+HV292.1     TR1 Chicken1 Chicken     t3
## 140 p-2-C11-S131               12 VAN+CCUG59168     TR7 Chicken2 Chicken     t2
## 170  p-2-F9-S165               11   CTX+HV292.1     TR2 Chicken2 Chicken     t2
## 126  p-2-B1-S109                2           VAN     TR6 Chicken2 Chicken     t2
## 127 p-2-B10-S118               11           VAN     TR6 Chicken2 Chicken     t2
## 46   TR3-22-S235                6       HV292.1     TR3 Chicken1 Chicken   <NA>
## 165  p-2-F4-S160                6   CTX+HV292.1     TR2 Chicken2 Chicken   <NA>
## 3    TR1-13-S171               -3   CTX+HV292.1     TR1 Chicken1 Chicken   pret
## 147  p-2-C8-S128                9 VAN+CCUG59168     TR7 Chicken2 Chicken   <NA>
## 31   TR2-28-S286               12           CTX     TR2 Chicken1 Chicken     t2
## 66   TR4-31-S278               15           VAN     TR4 Chicken1 Chicken     t2
## 2    IR1-35-S370               19   CTX+HV292.1     TR1 Chicken1 Chicken     t2
## 42   TR3-18-S293                2       HV292.1     TR3 Chicken1 Chicken     t2
## 139 p-2-C10-S130               11 VAN+CCUG59168     TR7 Chicken2 Chicken     t2
## 14   TR1-28-S289               12   CTX+HV292.1     TR1 Chicken1 Chicken     t2
## 22   TR2-15-S288               -1           CTX     TR2 Chicken1 Chicken   <NA>
## 167  p-2-F6-S162                8   CTX+HV292.1     TR2 Chicken2 Chicken   <NA>
## 164  p-2-F3-S159                5   CTX+HV292.1     TR2 Chicken2 Chicken     t5
## 156  p-2-D6-S138                8     UNTREATED      CR Chicken2 Chicken   <NA>
## 27   TR2-20-S271                4           CTX     TR2 Chicken1 Chicken     t4
## 24   TR2-17-S248                1           CTX     TR2 Chicken1 Chicken     t1
## 118 p-2-A12-S108               13       HV292.1     TR5 Chicken2 Chicken     t2
## 201  p-3-B4-S208               -2     UNTREATED      CR Chicken2 Chicken   pret
## 160  p-2-F1-S157                3   CTX+HV292.1     TR2 Chicken2 Chicken     t3
## 152  p-2-D2-S134                4     UNTREATED      CR Chicken2 Chicken     t4
## 19   TR1-46-S170               30   CTX+HV292.1     TR1 Chicken1 Chicken     t4
## 43   TR3-19-S266                3       HV292.1     TR3 Chicken1 Chicken     t3
## 72   TR5-14-S267               -2 VAN+CCUG59168     TR5 Chicken1 Chicken   pret
## 122  p-2-A6-S102                7       HV292.1     TR5 Chicken2 Chicken   <NA>
## 120   p-2-A3-S99                4       HV292.1     TR5 Chicken2 Chicken     t4
## 69   TR4-40-S242               24           VAN     TR4 Chicken1 Chicken     t3
## 186  p-2-H3-S183                5     CCUG59168     TR4 Chicken2 Chicken     t5
## 148  p-2-C9-S129               10 VAN+CCUG59168     TR7 Chicken2 Chicken     t2
## 183 p-2-H10-S190               12     CCUG59168     TR4 Chicken2 Chicken     t2
## 137  p-2-B9-S117               10           VAN     TR6 Chicken2 Chicken     t2
## 78   TR5-21-S185                5 VAN+CCUG59168     TR5 Chicken1 Chicken     t5
## 172 p-2-G10-S178               12           CTX     TR3 Chicken2 Chicken     t2
## 119   p-2-A2-S98                3       HV292.1     TR5 Chicken2 Chicken     t3
## 12   TR1-22-S282                6   CTX+HV292.1     TR1 Chicken1 Chicken   <NA>
## 200  p-3-B3-S207               -1     UNTREATED      CR Chicken2 Chicken   <NA>
## 116 p-2-A10-S106               11       HV292.1     TR5 Chicken2 Chicken     t2
## 76    TR5-19-S87                3 VAN+CCUG59168     TR5 Chicken1 Chicken     t3
## 145  p-2-C6-S126                7 VAN+CCUG59168     TR7 Chicken2 Chicken   <NA>
## 173 p-2-G11-S179               13           CTX     TR3 Chicken2 Chicken     t2
## 179  p-2-G7-S175                9           CTX     TR3 Chicken2 Chicken   <NA>
## 163  p-2-F2-S158                4   CTX+HV292.1     TR2 Chicken2 Chicken     t4
## 23   TR2-16-S184                0           CTX     TR2 Chicken1 Chicken     t1
## 177  p-2-G5-S173                7           CTX     TR3 Chicken2 Chicken   <NA>
## 129 p-2-B12-S120               13           VAN     TR6 Chicken2 Chicken     t2
## 141  p-2-C2-S122                3 VAN+CCUG59168     TR7 Chicken2 Chicken     t3
## 123  p-2-A7-S103                8       HV292.1     TR5 Chicken2 Chicken   <NA>
## 174  p-2-G2-S170                4           CTX     TR3 Chicken2 Chicken     t4
## 175  p-2-G3-S171                5           CTX     TR3 Chicken2 Chicken     t5
## 178  p-2-G6-S174                8           CTX     TR3 Chicken2 Chicken   <NA>
## 4    TR1-14-S274               -2   CTX+HV292.1     TR1 Chicken1 Chicken   pret
## 189  p-2-H6-S186                8     CCUG59168     TR4 Chicken2 Chicken   <NA>
## 104  p-1-A11-S11                1           CTX     TR3 Chicken2 Chicken     t1
## 125  p-2-A9-S105               10       HV292.1     TR5 Chicken2 Chicken     t2
## 190  p-2-H7-S187                9     CCUG59168     TR4 Chicken2 Chicken   <NA>
## 97   TR6-25-S378                9     CCUG59168     TR6 Chicken1 Chicken   <NA>
## 180  p-2-G8-S176               10           CTX     TR3 Chicken2 Chicken     t2
## 210  p-3-E4-S244               -2           CTX     TR3 Chicken2 Chicken   pret
## 121  p-2-A4-S100                5       HV292.1     TR5 Chicken2 Chicken     t5
## 117 p-2-A11-S107               12       HV292.1     TR5 Chicken2 Chicken     t2
## 157  p-2-D7-S139                9     UNTREATED      CR Chicken2 Chicken   <NA>
## 96   TR6-22-S152                6     CCUG59168     TR6 Chicken1 Chicken   <NA>
## 171  p-2-G1-S169                3           CTX     TR3 Chicken2 Chicken     t3
## 49   TR3-31-S283               15       HV292.1     TR3 Chicken1 Chicken     t2
## 103  TR6-46-S176               30     CCUG59168     TR6 Chicken1 Chicken     t4
## 15   TR1-31-S153               15   CTX+HV292.1     TR1 Chicken1 Chicken     t2
## 40   TR3-16-S296                0       HV292.1     TR3 Chicken1 Chicken     t1
## 115   p-2-A1-S97                2       HV292.1     TR5 Chicken2 Chicken     t2
## 74    TR5-17-S82                1 VAN+CCUG59168     TR5 Chicken1 Chicken     t1
##     LibrarySize
## 191           6
## 6            39
## 59           41
## 33           50
## 48           67
## 5            69
## 87           69
## 197         131
## 38          683
## 151        3738
## 95         4493
## 215        5758
## 181        6822
## 194        6836
## 51         8291
## 9          8608
## 112        8698
## 100        8828
## 146        8875
## 37         9289
## 30         9469
## 8          9526
## 26         9625
## 1          9959
## 204       10389
## 83        10416
## 198       11172
## 217       11293
## 212       11348
## 176       11486
## 98        11501
## 208       11549
## 16        11620
## 205       11856
## 222       12409
## 211       12815
## 199       13057
## 65        13155
## 144       13160
## 206       13423
## 67        13435
## 209       13766
## 82        13865
## 41        13985
## 53        14010
## 102       14053
## 195       14255
## 113       14320
## 64        14385
## 128       14537
## 101       14568
## 28        14695
## 213       14894
## 185       14980
## 203       15414
## 34        15567
## 63        15575
## 25        15579
## 216       15581
## 7         16039
## 218       16084
## 57        16197
## 110       16438
## 111       16732
## 39        16798
## 182       17000
## 94        17578
## 29        17687
## 54        17927
## 214       18325
## 85        18683
## 114       18707
## 106       18848
## 11        18979
## 91        19186
## 50        19218
## 220       19267
## 86        19284
## 184       19403
## 188       19506
## 21        19593
## 61        19601
## 55        19829
## 44        19878
## 207       19937
## 131       20114
## 35        20231
## 149       20281
## 187       20384
## 196       20395
## 17        20499
## 89        20502
## 202       20612
## 56        20762
## 20        20811
## 62        20932
## 84        21108
## 133       21198
## 108       21367
## 221       21578
## 10        21664
## 45        21840
## 80        21903
## 88        22203
## 134       22433
## 169       22471
## 68        22476
## 135       22531
## 81        22776
## 130       22794
## 77        22881
## 107       22889
## 132       23306
## 138       23454
## 13        23652
## 168       23668
## 93        23734
## 150       23898
## 70        23992
## 90        24099
## 52        24381
## 71        24401
## 36        24536
## 166       24666
## 32        24687
## 192       24915
## 99        24982
## 73        24995
## 79        25035
## 162       25375
## 47        25824
## 154       26118
## 124       26381
## 143       26499
## 159       26625
## 60        26675
## 92        26814
## 105       27026
## 75        27037
## 153       27050
## 58        27279
## 161       27431
## 158       27440
## 155       27453
## 219       27516
## 109       27562
## 136       27700
## 193       27732
## 142       27843
## 18        27896
## 140       28030
## 170       28191
## 126       28291
## 127       28321
## 46        28605
## 165       28617
## 3         28694
## 147       28697
## 31        28760
## 66        28770
## 2         28797
## 42        28808
## 139       28947
## 14        28961
## 22        29079
## 167       29345
## 164       29390
## 156       29425
## 27        29590
## 24        30114
## 118       30184
## 201       30408
## 160       30439
## 152       30446
## 19        30955
## 43        31111
## 72        31225
## 122       31875
## 120       31908
## 69        31977
## 186       31989
## 148       32142
## 183       32147
## 137       32173
## 78        32730
## 172       32958
## 119       33039
## 12        33137
## 200       33197
## 116       33501
## 76        33683
## 145       33948
## 173       33973
## 179       34283
## 163       34361
## 23        34452
## 177       35029
## 129       35161
## 141       35214
## 123       35419
## 174       35527
## 175       36043
## 178       36583
## 4         36812
## 189       37048
## 104       37812
## 125       37879
## 190       37982
## 97        39174
## 180       39478
## 210       39627
## 121       39816
## 117       39889
## 157       41007
## 96        41151
## 171       43153
## 49        43601
## 103       46983
## 15        50961
## 40        51543
## 115       55364
## 74        76823
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


```r
# ps_rare %>%
#   subset_samples(! Reactor %in% c("DONOR")) %>% 
#   subset_samples(Model2 == "Chicken2") %>%
#   subset_samples(Day_of_Treatment <= 15) %>% 
#   filter_taxa(function(x) sum(x > 0) > 0, TRUE) -> chick2 # %>% 
# subset_samples(Reactor %!in%c("IR")) %>% 
# subset_samples() -> before_treat

# ps_rare %>%
#   subset_samples(! Reactor %in% c("DONOR", "CR")) %>% 
#   subset_samples(Model2 = "Chicken1") %>% 
#   filter_taxa(function(x) sum(x > 0) > 0, TRUE) -> chick1 # %>% 


# mat.train.svd.corr <- svd(mat.train.clr.corr)


ps_rare %>%
  subset_samples(! Reactor %in% c("DONOR")) %>% 
  subset_samples(Model2 %in% c("Chicken1","Chicken2")) %>% 
  filter_taxa(function(x) sum(x > 0) > 0, TRUE) -> chicks # %>% 

# rm(list=setdiff(ls(), "chicks"))
```

# linear models: CLR taxa - ASV:

Per reactor:


```r
function_lm_per_reactor <- function(chick2, var = "Model2_reactor_treatment_dose", group = "TR6_VAN90", transform = "log", rank_glom = "Genus", ranks = "Genus",  min_prev = 0.25, adj_pval_cut = 0.05, estimate_cut = 0, min_tot_ab = 0.1, padj_group = "rank", sub_test = TRUE, time_split = 4,  use_cluster = FALSE, max_clust = 9)
{
  #######------------------
  
  require(microViz); require(tidyverse)
  
  #######------------------
  
  # The code for `taxatree_models` is quite similar to tax_model. 
  # However, you might need to run `tax_prepend_ranks` to ensure that each taxon at each rank is always unique. 
  chick2 %>% 
    prune_samples(get_variable(chick2, var) %in% group,
                  .)  %>% 
    filter_taxa(function(x) sum(x > 0) > 0, TRUE) %>% 
    # speedyseq::tax_glom(taxrank = rank_glom) %>%
    physeq_glom_rename(taxrank = rank_glom, rename_ASV = rank_glom) %>% 
    microbiome::transform(transform = transform) %>% 
    # subset_samples(eval(as.name(var)) %in% as.name(group)) %>% 
    # subset_samples(Reactor_Treatment_Dose == "TR6_VAN90") %>% 
    microViz::phyloseq_validate() %>%
    # tax_fix() %>%
    # tax_prepend_ranks() %>%
    # tax_transform("clr", rank = "Strain", keep_counts = TRUE) %>% 
    # tax_transform(transform = transform, rank = rank_glom, keep_counts = TRUE) %>%
    microViz::tax_filter(min_total_abundance = min_tot_ab, min_prevalence = min_prev, undetected = 0, use_counts = FALSE) -> ps_fil
  
  #######------------------
  
  lms <- ps_fil %>% 
    # tax_transform(trans = "log2", chain = TRUE, zero_replace = "halfmin") %>% 
    microViz::taxatree_models(
      type = lm, 
      ranks = ranks,
      variables = c('Day_of_Treatment')
    )
  
  #######------------------
  
  lms %>% 
    taxatree_models2stats(.) %>% 
    taxatree_stats_p_adjust(method = "BH", grouping = padj_group) %>% 
    .$taxatree_stats %>% 
    # mutate(taxon =  stringr::str_remove(taxon, "[PFGCOS]: ")) %>%  #-> tmp #%>% 
    # arrange(starts_with("p.adj")) %>%
    select(taxon,p.adj.BH.rank, everything()) -> tmp
  
  #######------------------
  
  tmp %>% 
    # filter(p.adj.BH.rank <= adj_pval_cut) %>% 
    # filter(rank == "Genus") %>% 
    filter(estimate > estimate_cut) %>% 
    mutate(comp = "pos") -> pos_taxa
  
  tmp %>% 
    # filter(p.adj.BH.rank <= adj_pval_cut) %>% 
    # filter(rank == "Genus") %>% 
    filter(estimate < -estimate_cut) %>% 
    mutate(comp = "neg") -> neg_taxa
  
  #######------------------
  
  rbind(pos_taxa,
        neg_taxa) %>% 
    filter(p.adj.BH.rank <= adj_pval_cut) %>%
    select(comp, everything())  -> overall_taxa_long
  
  overall_taxa_long %>% 
    select(term, rank, taxon, comp, p.adj.BH.rank) %>% 
    pivot_wider(   names_from = comp,
                   values_from = p.adj.BH.rank) %>% 
    rowwise() %>%
    mutate(pos = ifelse("pos" %in% names(.), pos, NA),
           neg = ifelse("neg" %in% names(.), neg, NA)) %>% 
    #        p2_pos = ifelse("p2_pos" %in% names(.), p2_neg, NA),
    #        p2_neg = ifelse("p2_neg" %in% names(.), p2_neg, NA)) %>% 
    mutate(class = case_when(!is.na(pos) ~ "inc",
                             !is.na(neg) ~ "dec",
                             TRUE ~ "unclassified")) %>% 
    mutate(Group = group)-> overall_taxa
  # select(-pos) %>% 
  # add_column(., !!!c("pos", "neg")[setdiff(c("pos", "neg"), names(overall_taxa))])
  
  
  #######------------------
  
  if ( use_cluster ==TRUE) {
    
    
    ps_fil$ps %>% 
      t() %>% 
      otu_table() -> t_otu
    
    t_otu %>% 
      Hmisc::rcorr() -> cor_tax
    
    #######------------------
    
    rows.cor <- cor((t_otu),  method = "spearman")
    
    hclust.row <- hclust(as.dist(1-rows.cor),"ward.D2" )
    
    # factoextra::fviz_nbclust(rows.cor, kmeans, method= 'gap_stat') -> outt
    # 
    # n_clust<-outt$data
    
    #######------------------
    
    # library(fpc)
    # pamk.best <- pamk(rows.cor)
    # cat("number of clusters estimated by optimum average silhouette width:", pamk.best$nc, "\n")
    # plot(cluster::pam(rows.cor, pamk.best$nc))
    # 
    #######------------------
    
    library(mclust)
    # Run the function to see how many clusters
    # it finds to be optimal, set it to search for
    # at least 1 model and up 20.
    d_clust <- Mclust(as.matrix(rows.cor), G=1:max_clust)
    m.best <- dim(d_clust$z)[2]
    cat("model-based optimal number of clusters:", m.best, "\n")
    
    cutree(hclust.row, k = m.best) %>% 
      data.frame() %>% 
      rownames_to_column('taxon') -> taxa_gpes
    # mutate(taxon =  stringr::str_remove(taxon, "[PFCGOS]: ")) -> taxa_gpes
    
    colnames(taxa_gpes)  <- c("taxon", "cluster")
    
    
    out = list(); out_heat = list()
    
    for (cluster_nb in taxa_gpes$cluster %>%  unique()){
      
      print(cluster_nb)
      
      taxa_gpes %>% 
        filter(cluster == cluster_nb) %>% 
        pull(taxon) %>% 
        as.character()-> tax
      
      if(tax %>%  length() > 1 )
      {
        chick2 %>% 
          prune_samples(get_variable(chick2, var) %in% group,
                        .)  %>% 
          filter_taxa(function(x) sum(x > 0) > 0, TRUE) %>% 
          microbiome::transform(., "compositional") %>% 
          subset_taxa(.,  eval(as.name(ranks)) %in% tax) %>% 
          
          # subset_taxa( Genus %in% tax) %>% 
          # subset_taxa(Genus %in% c("Proteus", "Ligilactobacillus")) %>% 
          microbiomeutilities::plot_area(., xvar="Day_of_Treatment", 
                                         level = ranks,
                                         facet.by =ranks,
                                         abund.thres = 0,
                                         prev.thres=0,
                                         fill.colors = microViz::distinct_palette(20, pal = "kelly"),
                                         ncol=6,
                                         nrow=4) -> out[[cluster_nb]]
        
      }
      
      # chick2 %>% 
      #   prune_samples(get_variable(chick2, var) %in% group,
      #                 .)  %>% 
      #   filter_taxa(function(x) sum(x > 0) > 0, TRUE) %>% 
      #   microbiome::transform(., "compositional") %>% 
      #     subset_taxa(.,  eval(as.name(ranks)) %in% tax) %>% 
      #   phyloseq_ampvis_heatmap(physeq = ., tax_aggregate = ifelse(ranks == "Strains", "Species",ranks), 
      #                           transform = FALSE,
      #                           plot_values = FALSE,
      #                           facet_by = c("Model2", "Model", "Fermentation","Reactor_Treatment_Dose"),
      #                           group_by = "Day_of_Treatment",#"sample_name",
      #                           ntax = Inf
      #   ) -> out_heat[[cluster_nb]]
      # 
      
    }
    
    out <- list("ps_fil" = ps_fil,
                "model" = tmp,
                "overall_taxa" = overall_taxa,
                "overall_taxa_long" = overall_taxa_long,
                "taxa_cluster" = taxa_gpes,
                "out_heat" = out_heat,
                "out" = out)
    
    
  }else{
    out <- list("ps_fil" = ps_fil,
                "model" = tmp,
                "overall_taxa" = overall_taxa,
                "overall_taxa_long" = overall_taxa_long)
    
  }
  #######------------------
  # then for the stationary ones, we can test on -4 to 5 and 5 to end.
  
  if(sub_test == TRUE){
    
    #######------------------
    
    ps_fil %>% 
      ps_filter(Day_of_Treatment <= time_split) %>% 
      # tax_transform(trans = "log2", chain = TRUE, zero_replace = "halfmin") %>% 
      taxatree_models(
        type = lm, 
        ranks = ranks,
        variables = c('Day_of_Treatment')
      ) %>% 
      taxatree_models2stats(.) %>% 
      taxatree_stats_p_adjust(method = "BH", grouping = padj_group) %>% 
      .$taxatree_stats %>% 
      # mutate(taxon =  stringr::str_remove(taxon, "[PFCGOS]: ")) %>% 
      arrange(p.adj.BH.rank) %>% 
      select(taxon,p.adj.BH.rank, everything()) -> psfil_1
    
    psfil_1 %>% 
      # filter(p.adj.BH.rank <= adj_pval_cut) %>% 
      # filter(rank == "Genus") %>% 
      filter(estimate > estimate_cut) %>% 
      mutate(comp = "p1_pos") -> fil1_pos_taxa
    
    psfil_1 %>% 
      # filter(p.adj.BH.rank <= adj_pval_cut) %>% 
      # filter(rank == "Genus") %>% 
      filter(estimate < -estimate_cut) %>% 
      mutate(comp = "p1_neg") -> fil1_neg_taxa
    
    #######------------------
    
    ps_fil %>% 
      ps_filter(Day_of_Treatment > time_split) %>% 
      # tax_transform(trans = "log2", chain = TRUE, zero_replace = "halfmin") %>% 
      taxatree_models(
        type = lm, 
        ranks = ranks,
        variables = c('Day_of_Treatment')
      ) %>% 
      taxatree_models2stats(.) %>% 
      taxatree_stats_p_adjust(method = "BH", grouping = padj_group) %>% 
      .$taxatree_stats %>% 
      # mutate(taxon =  stringr::str_remove(taxon, "[PFCGOS]: ")) %>% 
      arrange(p.adj.BH.rank) %>% 
      select(taxon,p.adj.BH.rank, everything()) -> psfil_2
    
    psfil_2 %>% 
      # filter(p.adj.BH.rank <= adj_pval_cut) %>% 
      # filter(rank == "Genus") %>% 
      filter(estimate > estimate_cut) %>% 
      mutate(comp = "p2_pos") -> fil2_pos_taxa
    
    psfil_2 %>% 
      # filter(p.adj.BH.rank <= adj_pval_cut) %>% 
      # filter(rank == "Genus") %>% 
      filter(estimate < -estimate_cut) %>% 
      mutate(comp = "p2_neg") -> fil2_neg_taxa
    
    #######------------------
    
    
    rbind(pos_taxa,
          neg_taxa,
          fil1_pos_taxa,
          fil1_neg_taxa,
          fil2_pos_taxa,
          fil2_neg_taxa) %>% 
      filter(p.adj.BH.rank <= adj_pval_cut) -> overall_taxa_long
    
    overall_taxa_long %>% 
      select(term, rank, taxon, comp, p.adj.BH.rank) %>% 
      pivot_wider(   names_from = comp,
                     values_from = p.adj.BH.rank) %>% 
      rowwise() %>%
      mutate(pos = ifelse("pos" %in% names(.), pos, NA),
             neg = ifelse("neg" %in% names(.), neg, NA)) %>% 
      mutate(p1_pos = ifelse("p1_pos" %in% names(.), p1_pos, NA),
             p1_neg = ifelse("p1_neg" %in% names(.), p1_neg, NA),
             p2_pos = ifelse("p2_pos" %in% names(.), p2_pos, NA),
             p2_neg = ifelse("p2_neg" %in% names(.), p2_neg, NA)) %>% 
      mutate(class = case_when(
        # increasing = signif pos or signif p1pos AND signif p2_pos
        !is.na(pos) | (!is.na(p1_pos) & !is.na(p2_pos)) ~ "inc",
        # secondary sucesional - inV = p1_pos AND p2Pos
        (!is.na(p1_pos) & !is.na(p2_neg)) | (!is.na(neg) & !is.na(p1_pos) & is.na(p1_neg) & is.na(p2_neg) & is.na(p2_pos) | (!is.na(p1_pos) & is.na(pos) &is.na(neg) & is.na(p2_pos) & is.na(p2_neg))) ~ "inV",
        # decrasing: signif neg AND ns p1pos and ns p2pos OR signif p1 neg AND 
        (!is.na(neg)  & (is.na(p1_pos) & is.na(p2_pos))) | (!is.na(p1_neg) & is.na(p1_pos) | !is.na(p2_neg) & is.na(pos) & is.na(neg) & is.na(p2_pos) & is.na(p1_pos) & is.na(p1_neg)) ~ "dec",
        # recovering= signif p1 neg ADN signif p2_pos OR signif neg AND signif p2_pos OR ns pos AND neg and p1pos OR 
        is.na(p1_neg) & !is.na(p2_pos) & is.na(pos) & is.na(p2_neg) ~ "V",
        TRUE ~ "unclassified"))  %>% 
      mutate(Group = group) -> overall_taxa
    
    
    
    
    if ( use_cluster ==TRUE) {
      
      out <- list("ps_fil" = ps_fil,
                  "model" = tmp,
                  "model_split1" = psfil_1,
                  "model_split2" = psfil_2,
                  "overall_taxa_long" = overall_taxa_long,
                  "overall_taxa" = overall_taxa,
                  "taxa_cluster" = taxa_gpes)
    }
    else{
      
      out <- list("ps_fil" = ps_fil,
                  "model" = tmp,
                  "model_split1" = psfil_1,
                  "model_split2" = psfil_2,
                  "overall_taxa_long" = overall_taxa_long,
                  "overall_taxa" = overall_taxa)
    }
  }
  
  #######------------------
  
  
  
  
  #######------------------
  
  return(out)
}
```



```r
# chicks %>%
#   function_lm_per_reactor(transform = "compositional", #compositional, 
#                           time_split = 3, rank_glom = "Genus",
#                           ranks = "Genus") -> out
# 
# out$overall_taxa %>% 
#   filter(class == "inV")
```


```r
setNames(vector("list", 
                length(physeq_get_unique(chicks, "Model2_reactor_treatment_dose"))), 
         physeq_get_unique(chicks, "Model2_reactor_treatment_dose")) -> lis
```


```r
for (grp in names(lis)){
  print(grp)
  
  function_lm_per_reactor(chicks, time_split = 3, rank_glom = "Genus",  
                          transform = "log", min_prev = 0.25, adj_pval_cut = 0.05, estimate_cut = 0, min_tot_ab = 0.01,
                          ranks = "Genus", group = grp) -> lis[[grp]]
}
```

```
## [1] "Chicken1_TR4_VAN90"
```

```
## Loading required package: microViz
```

```
## 
## microViz version 0.9.7 - Copyright (C) 2022 David Barnett
## * Website: https://david-barnett.github.io/microViz/
## * Useful? For citation info, run: citation('microViz')
## * Silence: suppressPackageStartupMessages(library(microViz))
```

```
## NAs detected in phyloseq tax_table:
## Consider using tax_fix() to make taxa uniquely identifiable
```

```
## Proportional min_prevalence given: 0.25 --> min 4/16 samples.
```

```
## Proportional min_total_abundance given: 0.01 --> min 16/1579.54612812283 reads.
```

```
## 2022-12-14 04:15:20 - modelling at rank: Genus
```

```
## NAs detected in phyloseq tax_table:
## Consider using tax_fix() to make taxa uniquely identifiable
```

```
## 2022-12-14 04:15:21 - modelling at rank: Genus
```

```
## NAs detected in phyloseq tax_table:
## Consider using tax_fix() to make taxa uniquely identifiable
```

```
## 2022-12-14 04:15:22 - modelling at rank: Genus
```

```
## NAs detected in phyloseq tax_table:
## Consider using tax_fix() to make taxa uniquely identifiable
```

```
## [1] "Chicken2_TR7_VAN+CCUG5916890"
```

```
## NAs detected in phyloseq tax_table:
## Consider using tax_fix() to make taxa uniquely identifiable
```

```
## Proportional min_prevalence given: 0.25 --> min 4/16 samples.
```

```
## Proportional min_total_abundance given: 0.01 --> min 19/1855.11748323388 reads.
```

```
## 2022-12-14 04:15:22 - modelling at rank: Genus
```

```
## NAs detected in phyloseq tax_table:
## Consider using tax_fix() to make taxa uniquely identifiable
```

```
## 2022-12-14 04:15:23 - modelling at rank: Genus
```

```
## NAs detected in phyloseq tax_table:
## Consider using tax_fix() to make taxa uniquely identifiable
```

```
## 2022-12-14 04:15:23 - modelling at rank: Genus
```

```
## NAs detected in phyloseq tax_table:
## Consider using tax_fix() to make taxa uniquely identifiable
```

```
## [1] "Chicken2_TR6_VAN90"
```

```
## NAs detected in phyloseq tax_table:
## Consider using tax_fix() to make taxa uniquely identifiable
```

```
## Proportional min_prevalence given: 0.25 --> min 5/17 samples.
```

```
## Proportional min_total_abundance given: 0.01 --> min 20/1919.75109904389 reads.
```

```
## 2022-12-14 04:15:24 - modelling at rank: Genus
```

```
## NAs detected in phyloseq tax_table:
## Consider using tax_fix() to make taxa uniquely identifiable
```

```
## 2022-12-14 04:15:25 - modelling at rank: Genus
```

```
## NAs detected in phyloseq tax_table:
## Consider using tax_fix() to make taxa uniquely identifiable
```

```
## 2022-12-14 04:15:25 - modelling at rank: Genus
```

```
## NAs detected in phyloseq tax_table:
## Consider using tax_fix() to make taxa uniquely identifiable
```

```
## [1] "Chicken1_TR5_VAN+CCUG5916890"
```

```
## NAs detected in phyloseq tax_table:
## Consider using tax_fix() to make taxa uniquely identifiable
```

```
## Proportional min_prevalence given: 0.25 --> min 4/16 samples.
```

```
## Proportional min_total_abundance given: 0.01 --> min 19/1816.45221856901 reads.
```

```
## 2022-12-14 04:15:26 - modelling at rank: Genus
```

```
## NAs detected in phyloseq tax_table:
## Consider using tax_fix() to make taxa uniquely identifiable
```

```
## 2022-12-14 04:15:26 - modelling at rank: Genus
```

```
## NAs detected in phyloseq tax_table:
## Consider using tax_fix() to make taxa uniquely identifiable
```

```
## 2022-12-14 04:15:27 - modelling at rank: Genus
```

```
## NAs detected in phyloseq tax_table:
## Consider using tax_fix() to make taxa uniquely identifiable
```

```
## [1] "Chicken2_TR3_CTX20"
```

```
## NAs detected in phyloseq tax_table:
## Consider using tax_fix() to make taxa uniquely identifiable
```

```
## Proportional min_prevalence given: 0.25 --> min 5/17 samples.
```

```
## Proportional min_total_abundance given: 0.01 --> min 35/3406.96384670363 reads.
```

```
## 2022-12-14 04:15:27 - modelling at rank: Genus
```

```
## NAs detected in phyloseq tax_table:
## Consider using tax_fix() to make taxa uniquely identifiable
```

```
## 2022-12-14 04:15:28 - modelling at rank: Genus
```

```
## NAs detected in phyloseq tax_table:
## Consider using tax_fix() to make taxa uniquely identifiable
```

```
## 2022-12-14 04:15:29 - modelling at rank: Genus
```

```
## NAs detected in phyloseq tax_table:
## Consider using tax_fix() to make taxa uniquely identifiable
```

```
## [1] "Chicken1_TR6_CCUG59168"
```

```
## NAs detected in phyloseq tax_table:
## Consider using tax_fix() to make taxa uniquely identifiable
```

```
## Proportional min_prevalence given: 0.25 --> min 4/16 samples.
```

```
## Proportional min_total_abundance given: 0.01 --> min 32/3108.52957172442 reads.
```

```
## 2022-12-14 04:15:30 - modelling at rank: Genus
```

```
## NAs detected in phyloseq tax_table:
## Consider using tax_fix() to make taxa uniquely identifiable
```

```
## 2022-12-14 04:15:30 - modelling at rank: Genus
```

```
## NAs detected in phyloseq tax_table:
## Consider using tax_fix() to make taxa uniquely identifiable
```

```
## 2022-12-14 04:15:31 - modelling at rank: Genus
```

```
## NAs detected in phyloseq tax_table:
## Consider using tax_fix() to make taxa uniquely identifiable
```

```
## [1] "Chicken1_TR2_CTX20"
```

```
## NAs detected in phyloseq tax_table:
## Consider using tax_fix() to make taxa uniquely identifiable
```

```
## Proportional min_prevalence given: 0.25 --> min 4/16 samples.
```

```
## Proportional min_total_abundance given: 0.01 --> min 24/2361.87767430103 reads.
```

```
## 2022-12-14 04:15:31 - modelling at rank: Genus
```

```
## NAs detected in phyloseq tax_table:
## Consider using tax_fix() to make taxa uniquely identifiable
```

```
## 2022-12-14 04:15:32 - modelling at rank: Genus
```

```
## NAs detected in phyloseq tax_table:
## Consider using tax_fix() to make taxa uniquely identifiable
```

```
## 2022-12-14 04:15:33 - modelling at rank: Genus
```

```
## NAs detected in phyloseq tax_table:
## Consider using tax_fix() to make taxa uniquely identifiable
```

```
## [1] "Chicken1_TR1_CTX+HV292.120"
```

```
## NAs detected in phyloseq tax_table:
## Consider using tax_fix() to make taxa uniquely identifiable
```

```
## Proportional min_prevalence given: 0.25 --> min 4/16 samples.
```

```
## Proportional min_total_abundance given: 0.01 --> min 28/2784.47408074295 reads.
```

```
## 2022-12-14 04:15:33 - modelling at rank: Genus
```

```
## NAs detected in phyloseq tax_table:
## Consider using tax_fix() to make taxa uniquely identifiable
```

```
## 2022-12-14 04:15:34 - modelling at rank: Genus
```

```
## NAs detected in phyloseq tax_table:
## Consider using tax_fix() to make taxa uniquely identifiable
```

```
## 2022-12-14 04:15:34 - modelling at rank: Genus
```

```
## NAs detected in phyloseq tax_table:
## Consider using tax_fix() to make taxa uniquely identifiable
```

```
## [1] "Chicken2_TR4_CCUG59168"
```

```
## NAs detected in phyloseq tax_table:
## Consider using tax_fix() to make taxa uniquely identifiable
```

```
## Proportional min_prevalence given: 0.25 --> min 4/16 samples.
```

```
## Proportional min_total_abundance given: 0.01 --> min 37/3621.91681487319 reads.
```

```
## 2022-12-14 04:15:35 - modelling at rank: Genus
```

```
## NAs detected in phyloseq tax_table:
## Consider using tax_fix() to make taxa uniquely identifiable
```

```
## 2022-12-14 04:15:36 - modelling at rank: Genus
```

```
## NAs detected in phyloseq tax_table:
## Consider using tax_fix() to make taxa uniquely identifiable
```

```
## 2022-12-14 04:15:36 - modelling at rank: Genus
```

```
## NAs detected in phyloseq tax_table:
## Consider using tax_fix() to make taxa uniquely identifiable
```

```
## [1] "Chicken1_TR3_HV292.1"
```

```
## NAs detected in phyloseq tax_table:
## Consider using tax_fix() to make taxa uniquely identifiable
```

```
## Proportional min_prevalence given: 0.25 --> min 4/15 samples.
```

```
## Proportional min_total_abundance given: 0.01 --> min 28/2760.40675836608 reads.
```

```
## 2022-12-14 04:15:37 - modelling at rank: Genus
```

```
## NAs detected in phyloseq tax_table:
## Consider using tax_fix() to make taxa uniquely identifiable
```

```
## 2022-12-14 04:15:38 - modelling at rank: Genus
```

```
## NAs detected in phyloseq tax_table:
## Consider using tax_fix() to make taxa uniquely identifiable
```

```
## 2022-12-14 04:15:38 - modelling at rank: Genus
```

```
## NAs detected in phyloseq tax_table:
## Consider using tax_fix() to make taxa uniquely identifiable
```

```
## [1] "Chicken2_TR5_HV292.1"
```

```
## NAs detected in phyloseq tax_table:
## Consider using tax_fix() to make taxa uniquely identifiable
```

```
## Proportional min_prevalence given: 0.25 --> min 4/16 samples.
```

```
## Proportional min_total_abundance given: 0.01 --> min 39/3807.63757397037 reads.
```

```
## 2022-12-14 04:15:39 - modelling at rank: Genus
```

```
## NAs detected in phyloseq tax_table:
## Consider using tax_fix() to make taxa uniquely identifiable
```

```
## 2022-12-14 04:15:40 - modelling at rank: Genus
```

```
## NAs detected in phyloseq tax_table:
## Consider using tax_fix() to make taxa uniquely identifiable
```

```
## 2022-12-14 04:15:40 - modelling at rank: Genus
```

```
## NAs detected in phyloseq tax_table:
## Consider using tax_fix() to make taxa uniquely identifiable
```

```
## [1] "Chicken2_CR_UNTREATED"
```

```
## NAs detected in phyloseq tax_table:
## Consider using tax_fix() to make taxa uniquely identifiable
```

```
## Proportional min_prevalence given: 0.25 --> min 5/17 samples.
```

```
## Proportional min_total_abundance given: 0.01 --> min 41/4068.53193091619 reads.
```

```
## 2022-12-14 04:15:41 - modelling at rank: Genus
```

```
## NAs detected in phyloseq tax_table:
## Consider using tax_fix() to make taxa uniquely identifiable
```

```
## 2022-12-14 04:15:42 - modelling at rank: Genus
```

```
## NAs detected in phyloseq tax_table:
## Consider using tax_fix() to make taxa uniquely identifiable
```

```
## 2022-12-14 04:15:43 - modelling at rank: Genus
```

```
## NAs detected in phyloseq tax_table:
## Consider using tax_fix() to make taxa uniquely identifiable
```

```
## [1] "Chicken2_TR2_CTX+HV292.120"
```

```
## NAs detected in phyloseq tax_table:
## Consider using tax_fix() to make taxa uniquely identifiable
```

```
## Proportional min_prevalence given: 0.25 --> min 5/17 samples.
```

```
## Proportional min_total_abundance given: 0.01 --> min 39/3861.25842229525 reads.
```

```
## 2022-12-14 04:15:43 - modelling at rank: Genus
```

```
## NAs detected in phyloseq tax_table:
## Consider using tax_fix() to make taxa uniquely identifiable
```

```
## 2022-12-14 04:15:44 - modelling at rank: Genus
```

```
## NAs detected in phyloseq tax_table:
## Consider using tax_fix() to make taxa uniquely identifiable
```

```
## 2022-12-14 04:15:45 - modelling at rank: Genus
```

```
## NAs detected in phyloseq tax_table:
## Consider using tax_fix() to make taxa uniquely identifiable
```

```r
do.call(rbind, lapply(lis, "[[", "overall_taxa")) %>% 
  ungroup() %>% 
  # select(taxon, Group, class) %>% 
  rename(OTU = taxon,
         Model2_reactor_treatment_dose = Group,
         category = class)-> map_class

map_class %>% 
  filter(!category %in% c("V","inV", "dec", "inc"))
```

```
## # A tibble: 0 × 11
## # … with 11 variables: term <fct>, rank <fct>, OTU <chr>, pos <dbl>, neg <dbl>,
## #   p2_pos <dbl>, p2_neg <dbl>, p1_pos <dbl>, p1_neg <dbl>, category <chr>,
## #   Model2_reactor_treatment_dose <chr>
```




```r
chicks %>% 
  # prune_samples(get_variable(chick2, var) %in% group,
  #               .)  %>%
  filter_taxa(function(x) sum(x > 0) > 0, TRUE) %>% 
  transform_sample_counts(function(x) x/sum(x) * 100) %>% 
  physeq_glom_rename(taxrank = "Genus", rename_ASV = "Genus") -> ps_filt
```



```r
map_class %>%
  select(category, Model2_reactor_treatment_dose, everything()) %>%
  # filter(category %in% c("V"),
  #        Reactor_Treatment_Dose == "CR_UNTREATED") %>%
  filter(!is.na(pos) & !is.na(p2_pos)) %>% 
  distinct(OTU, Model2_reactor_treatment_dose) -> sel

ps_filt %>% 
  speedyseq::psmelt() %>%
  # filter(OTU %in% sel$OTU,
  # Reactor_Treatment_Dose %in% sel$Reactor_Treatment_Dose) %>%
  # filter(Reactor_Treatment_Dose %in% c("TR6_VAN90") ) %>%
  select(OTU, Abundance, Day_of_Treatment, Model2_reactor_treatment_dose) %>% #!!rank_glom,
  left_join(map_class,
            by = c("Model2_reactor_treatment_dose", "OTU")) %>% 
  arrange(Model2_reactor_treatment_dose) %>% 
  ggplot(data = ., aes_string(color = "OTU", fill = "OTU")) +
  geom_area(aes_string(x = "Day_of_Treatment", y = "Abundance")) + 
  facet_grid(category ~  Model2_reactor_treatment_dose, scales = "free") + theme(legend.position = "none")
```

![](16S-taxa-classification_Figs_BP_finetuning_files/figure-html/unnamed-chunk-12-1.png)<!-- -->

same as Hanna but with long df with wide time abundance:

```r
var = "Model2_reactor_treatment_dose"
transform = "log"
rank_glom = "Genus"
# ranks = "Genus"
min_prev = 0.25
adj_pval_cut = 0.05
padj_group = "rank"
sub_test = TRUE

#######------------------

require(tidyverse)

#######------------------

chicks %>% 
  # prune_samples(get_variable(chick2, var) %in% group,
  #               .)  %>%
  filter_taxa(function(x) sum(x > 0) > 0, TRUE) %>% 
  subset_samples(Day_of_Treatment <= 15) %>% 
  physeq_glom_rename(taxrank = rank_glom, rename_ASV = rank_glom) -> ps_filt
```



```r
ps_filt %>% 
  transform_sample_counts(function(x) x/sum(x) * 100) %>% 
  speedyseq::psmelt() %>% 
  select(OTU, Abundance, Model2, Day_of_Treatment, Model2_reactor_treatment_dose ,Day_of_Treatment, Treatment, Reactor_Treatment_Dose) %>% 
  arrange(Day_of_Treatment)  -> ps_filt2


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

```r
ps_filt %>%
  microbiome::transform(transform = "Z") %>%
  speedyseq::psmelt() %>%
  arrange(Day_of_Treatment) %>%
  group_by(OTU, Model2_reactor_treatment_dose) %>%
  # do(model = tseries::adf.test(.))
  # filter(OTU == "Alloprevotella" & Model2_reactor_treatment_dose == "Chicken1_TR1_CTX+HV292.120") %>% #pull(Abundance) -> Abundance
  
  # summarise(pval_ur.df2 = ifelse(n() <= 4,NA, urca::ur.df(Abundance, type = "none", selectlags = "AIC")$out)) %>%
  # summarise(pval_ur.df2_none = ifelse(n() <= 4,NA, erer::ur.df2(Abundance, type = "none", selectlags = "AIC")$out)  %>%
  summarise(pval_ur.df2_none =  ifelse(sum(Abundance > 0, na.rm = FALSE) <= 4,NA,coef(attributes(urca::ur.df(y = Abundance, type = "none", selectlags = "AIC"))$testreg)[1, 4]),
            pval_ur.df2_drift = ifelse(sum(Abundance > 0, na.rm = FALSE) <= 4, NA,coef(attributes(urca::ur.df(y = Abundance, type = "drift", selectlags = "AIC"))$testreg)[1, 4]),
            pval_ur.df2_trend = ifelse(sum(Abundance > 0, na.rm = FALSE)  <= 4,NA,coef(attributes(urca::ur.df(y = Abundance, type = "trend", selectlags = "AIC"))$testreg)[1, 4]),
            pval_adf.test = tseries::adf.test(Abundance, alternative = "stationary")$p.value,
            pval_kpss_type1 = ifelse(sum(Abundance > 0, na.rm = FALSE) <= 4,NA, aTSA::kpss.test(Abundance, lag.short = TRUE, output =FALSE) %>%  data.frame() %>%  rownames_to_column("type") %>% filter(type == "type 1") %>% pull(p.value)),
            pval_kpss_type2 = ifelse(sum(Abundance > 0, na.rm = FALSE)  <= 4,NA, aTSA::kpss.test(Abundance, lag.short = TRUE, output =FALSE) %>%  data.frame() %>%  rownames_to_column("type") %>% filter(type == "type 2") %>% pull(p.value)),
            pval_kpss_type3 = ifelse(sum(Abundance > 0, na.rm = FALSE)  <= 4,NA, aTSA::kpss.test(Abundance, lag.short = TRUE, output =FALSE) %>%  data.frame() %>%  rownames_to_column("type") %>% filter(type == "type 3") %>% pull(p.value)),
            pval_kpss.test_level = ifelse(sum(Abundance > 0, na.rm = FALSE)  <= 4,NA, tseries::kpss.test(Abundance, null = "Level")$p.value),
            pval_Box.test = ifelse(sum(Abundance > 0, na.rm = FALSE) <= 4,NA, stats::Box.test(Abundance, lag = 1, type = "Ljung-Box")$p.value) #Here we see a p-value much smaller than .01, thus we can reject the null hypothesis, indicating the time series does contain an autocorrelation.
            
            
            # pval_pp.test = ifelse(n() <= 4,NA, aTSA::pp.test (vec, lag.short = TRUE, output =FALSE) %>%  data.frame() %>%  rownames_to_column("type") %>% filter(type == "type 3")
            # pval_aTSAadf.tes =  ifelse(n() <= 4,NA, aTSA::adf.test())
  ) %>%
  # pval_CADFtest = ifelse(n() <= 12,NA,CADFtest::CADFtest(Abundance ~ 1,x = Abundance, type="trend", criterion = "AIC")$p.value)) %>%  # 0.000
  # summarise(pval_ur.df2_none = ifelse(n() <= 4,NA, erer::ur.df2(Abundance, type = "none", selectlags = "AIC")$out)) %>%
  # https://search.r-project.org/CRAN/refmans/erer/html/ur.df2.html
  # summarise(pval = ifelse(n() <= 4,NA, CADFtest::CADFtest(x = Abundance, type="trend", criterion = "AIC")$p.value)) %>%
  # group_by(Reactor_Treatment_Dose) %>%
  # mutate(pval_ur.df2_none = str_replace_all(pval_ur.df2_none, "[^*]", "")) %>%
  # rstatix::adjust_pvalue(p.col = "pval_ur.df2_none", method = "fdr") %>%
  # filter(pval.adj <= 0.05) %>%
  arrange(Model2_reactor_treatment_dose) %>%
  mutate(OTU_Reactor_Treatment_Dose = paste0(OTU,"_",Model2_reactor_treatment_dose)) -> ur.df2
```

```
## Registered S3 method overwritten by 'quantmod':
##   method            from
##   as.zoo.data.frame zoo 
## `summarise()` has grouped output by 'OTU'. You can override using the `.groups` argument.
```

```r
ur.df2
```

```
## # A tibble: 1,508 × 12
## # Groups:   OTU [116]
##    OTU   Model…¹ pval_…² pval_…³ pval_…⁴ pval_…⁵ pval_…⁶ pval_…⁷ pval_…⁸ pval_…⁹
##    <chr> <chr>     <dbl>   <dbl>   <dbl>   <dbl>   <dbl>   <dbl>   <dbl>   <dbl>
##  1 [Eub… Chicke…  NA      NA      NA     NaN        NA      NA      NA      NA  
##  2 [Eub… Chicke…   0.542   0.227   0.409   0.01      0.1     0.1     0.1     0.1
##  3 [Eub… Chicke…   0.231   0.354   0.410   0.600     0.1     0.1     0.1     0.1
##  4 [Rum… Chicke…   0.789   0.133   0.238   0.01      0.1     0.1     0.1     0.1
##  5 Acin… Chicke…  NA      NA      NA     NaN        NA      NA      NA      NA  
##  6 Aero… Chicke…  NA      NA      NA     NaN        NA      NA      NA      NA  
##  7 Alis… Chicke…  NA      NA      NA     NaN        NA      NA      NA      NA  
##  8 Allo… Chicke…  NA      NA      NA     NaN        NA      NA      NA      NA  
##  9 Anae… Chicke…  NA      NA      NA     NaN        NA      NA      NA      NA  
## 10 Anae… Chicke…  NA      NA      NA     NaN        NA      NA      NA      NA  
## # … with 1,498 more rows, 2 more variables: pval_Box.test <dbl>,
## #   OTU_Reactor_Treatment_Dose <chr>, and abbreviated variable names
## #   ¹​Model2_reactor_treatment_dose, ²​pval_ur.df2_none, ³​pval_ur.df2_drift,
## #   ⁴​pval_ur.df2_trend, ⁵​pval_adf.test, ⁶​pval_kpss_type1, ⁷​pval_kpss_type2,
## #   ⁸​pval_kpss_type3, ⁹​pval_kpss.test_level
```


```r
library(dplyr)
#define_groups:
ps_filt2 %>% 
  # select(Day_of_Treatment) %>% 
  distinct(Day_of_Treatment) %>% 
  mutate(period = case_when(Day_of_Treatment <= 0 ~ "pret",
                            Day_of_Treatment > 0 &  Day_of_Treatment <= 5 ~ "t1" ,
                            Day_of_Treatment > 5 &  Day_of_Treatment <= 10 ~ "t2",
                            Day_of_Treatment > 10  ~ "t3" )) -> day_period

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
  # mutate(period = ifelse(Day_of_Treatment <= 0, "pret",
  #                        ifelse(0 <  Day_of_treatment <= 5, "t1"
  #                        ))
  # mutate(period = case_when(Day_of_Treatment <= 0 ~ "pret",
  #                           0 <  Day_of_treatment <= 5 ~ "tt"
  #                           Day_of_treatment > 0 & Day_of_treatment <= 5 ~ "t1",
  #                           Day_of_treatment > 6 & Day_of_treatment <= 10 ~ "t2",
  #                           Day_of_treatment > 10 & Day_of_treatment <= 16 ~ "t3",
  #                           TRUE ~ "tunknown")) %>%
  group_by(Model2_reactor_treatment_dose, OTU, period) %>% 
  summarise(mean = base::mean(Abundance, na.rm = TRUE)) %>% 
  group_by(Model2_reactor_treatment_dose, OTU) %>% 
  pivot_wider(names_from = period, 
              values_from = mean) %>% 
  
  # original
  # summarise(category = 
  #           case_when(
  #             # pret == 0 & (t1 >0 & t2>0 & t3 >0) ~ "new",
  #                     (t1 > pret  &
  #                        t2 > pret &
  #                        t3 > pret) & (t3 > t2 | t3 > t1 ) & 
  #                       (t3/(pret+0.00000000001) > 1.5 & t2/(pret+0.00000000001) > 1.25)  ~ "prol",
  #                     t1 < pret  &
  #                       t3 < t1 ~ "dis",
#                     pret/t1 > 2 &
#                       t3/t1 > 1.25  ~ "recovering",
#                     (t1/(pret+0.00000000001) > 1.25 | t2/(pret+0.00000000001) > 1.25) &
#                       (t1/t3 > 1.25 | t2/t3 > 1.25) ~ "mid",
#                     # pret == 0 & t1>0 &t2>0 &t3>0 ~ "new",
#                     TRUE ~ "normal")) %>% 

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
  arrange(Model2_reactor_treatment_dose,category) -> classif
```

```
## `summarise()` has grouped output by 'Model2_reactor_treatment_dose', 'OTU'. You
## can override using the `.groups` argument.
## `summarise()` has grouped output by 'Model2_reactor_treatment_dose'. You can
## override using the `.groups` argument.
```

```r
classif$category %>%
  unique()
```

```
## [1] "dis"        "mid"        "normal"     "prol"       "recovering"
```



```r
classif %>% 
  rename("cat_period" = "category") %>% 
  mutate(OTU_Reactor_Treatment_Dose = paste0(OTU,"_",Model2_reactor_treatment_dose)) %>% 
  left_join(ur.df2,
            by = c("Model2_reactor_treatment_dose", "OTU"),
            suffix = c("", ".y")) %>%
  left_join(map_class %>%  rename("cat_lm" = "category"),
            by = c("Model2_reactor_treatment_dose", "OTU")) %>% 
  select(-term, -rank, ) -> tmp

tmp
```

```
## # A tibble: 1,508 × 21
## # Groups:   Model2_reactor_treatment_dose [13]
##    Model…¹ OTU   cat_p…² OTU_R…³ pval_…⁴ pval_…⁵ pval_…⁶ pval_…⁷ pval_…⁸ pval_…⁹
##    <chr>   <chr> <chr>   <chr>     <dbl>   <dbl>   <dbl>   <dbl>   <dbl>   <dbl>
##  1 Chicke… Hold… dis     Holdem… NA      NA      NA       0.01   NA         NA  
##  2 Chicke… UBA1… dis     UBA181… NA      NA      NA       0.858  NA         NA  
##  3 Chicke… [Eub… mid     [Eubac…  0.542   0.227   0.409   0.01    0.1        0.1
##  4 Chicke… Baci… mid     Bacill…  0.387   0.283   0.0573  0.712   0.0949     0.1
##  5 Chicke… Blau… mid     Blauti…  0.766   0.241   0.0716  0.804   0.1        0.1
##  6 Chicke… CAG-… mid     CAG-56…  0.0861  0.0638  0.754   0.496   0.1        0.1
##  7 Chicke… Flav… mid     Flavon…  0.910   0.0212  0.0111  0.0559  0.1        0.1
##  8 Chicke… Gord… mid     Gordon… NA      NA      NA       0.959  NA         NA  
##  9 Chicke… Lach… mid     Lachno…  0.912   0.0419  0.0151  0.912   0.1        0.1
## 10 Chicke… Lach… mid     Lachno…  0.547   0.0797  0.153   0.905   0.1        0.1
## # … with 1,498 more rows, 11 more variables: pval_kpss_type3 <dbl>,
## #   pval_kpss.test_level <dbl>, pval_Box.test <dbl>,
## #   OTU_Reactor_Treatment_Dose.y <chr>, pos <dbl>, neg <dbl>, p2_pos <dbl>,
## #   p2_neg <dbl>, p1_pos <dbl>, p1_neg <dbl>, cat_lm <chr>, and abbreviated
## #   variable names ¹​Model2_reactor_treatment_dose, ²​cat_period,
## #   ³​OTU_Reactor_Treatment_Dose, ⁴​pval_ur.df2_none, ⁵​pval_ur.df2_drift,
## #   ⁶​pval_ur.df2_trend, ⁷​pval_adf.test, ⁸​pval_kpss_type1, ⁹​pval_kpss_type2
```



```r
ps_filt2 %>%
  # filter(Reactor_Treatment_Dose %in% c("TR6_VAN90","TR4_CCUG59168", "TR5_HV292.1", "TR1_CTX+HV292.120")) %>%
  select(OTU, Abundance, Day_of_Treatment, Model2 ,Model2_reactor_treatment_dose) %>% 
  left_join(tmp,
            by = c("Model2_reactor_treatment_dose", "OTU")) %>% 
  left_join(mean_ab_OTU_reactors,
            by = c("Model2_reactor_treatment_dose", "OTU")) %>% 
  # filter(mean_ab <= 5) %>%
  filter(cat_period != "normal") %>%
  # Reactor_Treatment_Dose %in% ("TR4_CCUG59168")) %>% 
  arrange(Model2_reactor_treatment_dose) -> pre_plot

# pre_plot %>%
#   pull(OTU) %>%
#   unique() %>%
#   length() %>%
#   randomcoloR::distinctColorPalette(k = ., runTsne = TRUE) -> col


# names(col) <- pre_plot %>%
#   pull(OTU) %>%
#   unique() 


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



pre_plot %>% 
  ggplot(data = ., aes_string(color = "OTU", fill = "OTU")) +
  geom_area(aes_string(x = "Day_of_Treatment", y = "Abundance")) + 
  facet_grid(cat_period ~  Model2 + Model2_reactor_treatment_dose, scales = "free") + 
  scale_fill_manual(values =  col) + 
  scale_color_manual(values = col) -> p_cat


p_cat + theme(legend.position = "none") -> p_cat_no_leg

p_cat_no_leg
```

![](16S-taxa-classification_Figs_BP_finetuning_files/figure-html/unnamed-chunk-17-1.png)<!-- -->

```r
p_cat  %>% 
  ggpubr::get_legend() %>% 
  ggpubr::as_ggplot() -> p_cat_leg

p_cat_leg
```

![](16S-taxa-classification_Figs_BP_finetuning_files/figure-html/unnamed-chunk-17-2.png)<!-- -->

```r
p_cat_no_leg %>% 
  export::graph2ppt(append = TRUE,
                    width = 317.48031496 * 6,
                    height = 0.618 * 317.48031496 * 8 , paper = "A2",  scaling = 2,
                    file = out_pptx)
```

```
## Exported graph as ~/Desktop/16S-chick.pptx
```

```r
p_cat_leg %>% 
  export::graph2ppt(append = TRUE,
                    width = 317.48031496 * 6,
                    height = 0.618 * 317.48031496 * 8 , paper = "A2",  scaling = 2,
                    file = out_pptx)
```

```
## Exported graph as ~/Desktop/16S-chick.pptx
```


```r
ps_filt2 %>% 
  # filter(Reactor_Treatment_Dose %in% c("TR6_VAN90","TR4_CCUG59168", "TR5_HV292.1", "TR1_CTX+HV292.120")) %>% 
  select(OTU, Abundance, Day_of_Treatment, Model2 ,Model2_reactor_treatment_dose, Reactor_Treatment_Dose, Treatment, ) %>% 
  left_join(tmp,
            by = c("Model2_reactor_treatment_dose", "OTU")) %>%  
  mutate(cat_lm = replace_na(cat_lm, "normal")) %>% 
  left_join(day_period)  -> tmp_ready 
```

```
## Joining, by = "Day_of_Treatment"
```
Ljun-Box test is a hypothesis test that checks if a time series contains an autocorrelation. The null Hypothesis H0 is that the residuals are independently distributed. sign -> time pattern
tseries::adf.test(lc, alternative = "stationary") -> ns -> time pattern
urca::ur.df not sure sign -> pattern?
kpss.test() signif : pattern


```r
tmp_ready %>% 
  filter(pval_Box.test < 0.05 |  pval_adf.test > 0.05 | pval_kpss.test_level < 0.05) %>%
  arrange(OTU, Model2_reactor_treatment_dose) %>% 
  select(OTU, Model2_reactor_treatment_dose, cat_period, cat_lm,  pval_Box.test, pval_adf.test, pval_kpss.test_level) %>% 
  distinct(OTU, Model2_reactor_treatment_dose ,.keep_all = TRUE) #%>% 
```

```
##                               OTU Model2_reactor_treatment_dose cat_period
## 1      [Eubacterium] brachy group            Chicken1_TR2_CTX20        dis
## 2      [Eubacterium] brachy group          Chicken1_TR3_HV292.1       prol
## 3      [Eubacterium] brachy group            Chicken1_TR4_VAN90        dis
## 4      [Eubacterium] brachy group        Chicken1_TR6_CCUG59168       prol
## 5      [Eubacterium] brachy group            Chicken2_TR3_CTX20 recovering
## 6      [Eubacterium] brachy group        Chicken2_TR4_CCUG59168        mid
## 7      [Eubacterium] brachy group          Chicken2_TR5_HV292.1        mid
## 8      [Eubacterium] hallii group    Chicken1_TR1_CTX+HV292.120        mid
## 9      [Eubacterium] hallii group            Chicken1_TR2_CTX20        mid
## 10     [Eubacterium] hallii group          Chicken1_TR3_HV292.1     normal
## 11     [Eubacterium] hallii group            Chicken1_TR4_VAN90        dis
## 12     [Eubacterium] hallii group  Chicken1_TR5_VAN+CCUG5916890        dis
## 13     [Eubacterium] hallii group        Chicken1_TR6_CCUG59168     normal
## 14     [Eubacterium] hallii group         Chicken2_CR_UNTREATED recovering
## 15     [Eubacterium] hallii group    Chicken2_TR2_CTX+HV292.120     normal
## 16     [Eubacterium] hallii group            Chicken2_TR3_CTX20     normal
## 17     [Eubacterium] hallii group        Chicken2_TR4_CCUG59168       prol
## 18     [Eubacterium] hallii group          Chicken2_TR5_HV292.1        mid
## 19     [Eubacterium] hallii group            Chicken2_TR6_VAN90        dis
## 20    [Eubacterium] nodatum group    Chicken1_TR1_CTX+HV292.120     normal
## 21    [Eubacterium] nodatum group            Chicken1_TR2_CTX20        mid
## 22    [Eubacterium] nodatum group          Chicken1_TR3_HV292.1     normal
## 23    [Eubacterium] nodatum group            Chicken1_TR4_VAN90       prol
## 24    [Eubacterium] nodatum group  Chicken1_TR5_VAN+CCUG5916890       prol
## 25    [Eubacterium] nodatum group         Chicken2_CR_UNTREATED       prol
## 26    [Eubacterium] nodatum group    Chicken2_TR2_CTX+HV292.120        mid
## 27    [Eubacterium] nodatum group            Chicken2_TR3_CTX20        mid
## 28    [Eubacterium] nodatum group        Chicken2_TR4_CCUG59168        mid
## 29    [Eubacterium] nodatum group          Chicken2_TR5_HV292.1        mid
## 30    [Eubacterium] nodatum group            Chicken2_TR6_VAN90     normal
## 31    [Eubacterium] nodatum group  Chicken2_TR7_VAN+CCUG5916890     normal
## 32   [Ruminococcus] torques group    Chicken1_TR1_CTX+HV292.120     normal
## 33   [Ruminococcus] torques group            Chicken1_TR2_CTX20       prol
## 34   [Ruminococcus] torques group          Chicken1_TR3_HV292.1     normal
## 35   [Ruminococcus] torques group            Chicken1_TR4_VAN90        dis
## 36   [Ruminococcus] torques group  Chicken1_TR5_VAN+CCUG5916890        dis
## 37   [Ruminococcus] torques group        Chicken1_TR6_CCUG59168     normal
## 38   [Ruminococcus] torques group         Chicken2_CR_UNTREATED     normal
## 39   [Ruminococcus] torques group            Chicken2_TR3_CTX20     normal
## 40   [Ruminococcus] torques group        Chicken2_TR4_CCUG59168        dis
## 41   [Ruminococcus] torques group          Chicken2_TR5_HV292.1     normal
## 42   [Ruminococcus] torques group            Chicken2_TR6_VAN90        dis
## 43                  Acinetobacter            Chicken1_TR4_VAN90       prol
## 44                  Acinetobacter  Chicken1_TR5_VAN+CCUG5916890       prol
## 45                  Acinetobacter        Chicken1_TR6_CCUG59168       prol
## 46                  Acinetobacter  Chicken2_TR7_VAN+CCUG5916890        mid
## 47                     Aerococcus  Chicken1_TR5_VAN+CCUG5916890        dis
## 48                      Alistipes            Chicken2_TR3_CTX20        mid
## 49                      Alistipes  Chicken2_TR7_VAN+CCUG5916890     normal
## 50                   Anaerococcus        Chicken2_TR4_CCUG59168 recovering
## 51                    Anaerofilum          Chicken1_TR3_HV292.1        mid
## 52                    Anaerofilum         Chicken2_CR_UNTREATED     normal
## 53                    Anaerofilum    Chicken2_TR2_CTX+HV292.120     normal
## 54                    Anaerofilum            Chicken2_TR3_CTX20        mid
## 55                    Anaerofilum        Chicken2_TR4_CCUG59168     normal
## 56                    Anaerofilum          Chicken2_TR5_HV292.1       prol
## 57                   Anaerofustis    Chicken2_TR2_CTX+HV292.120     normal
## 58                   Anaerofustis        Chicken2_TR4_CCUG59168       prol
## 59                   Anaerofustis          Chicken2_TR5_HV292.1       prol
## 60                   Anaerostipes            Chicken1_TR2_CTX20 recovering
## 61                   Anaerostipes        Chicken1_TR6_CCUG59168        mid
## 62                   Anaerostipes            Chicken2_TR3_CTX20        mid
## 63                   Anaerostipes          Chicken2_TR5_HV292.1        dis
## 64                   Anaerostipes            Chicken2_TR6_VAN90        mid
## 65                  Angelakisella         Chicken2_CR_UNTREATED        mid
## 66                  Angelakisella    Chicken2_TR2_CTX+HV292.120        mid
## 67                  Angelakisella            Chicken2_TR3_CTX20     normal
## 68                  Angelakisella        Chicken2_TR4_CCUG59168 recovering
## 69                  Angelakisella          Chicken2_TR5_HV292.1 recovering
## 70                       Bacillus    Chicken1_TR1_CTX+HV292.120        mid
## 71                       Bacillus            Chicken1_TR2_CTX20       prol
## 72                       Bacillus          Chicken1_TR3_HV292.1        mid
## 73                       Bacillus        Chicken1_TR6_CCUG59168 recovering
## 74                       Bacillus         Chicken2_CR_UNTREATED     normal
## 75                       Bacillus    Chicken2_TR2_CTX+HV292.120     normal
## 76                       Bacillus            Chicken2_TR3_CTX20     normal
## 77                       Bacillus        Chicken2_TR4_CCUG59168     normal
## 78                       Bacillus          Chicken2_TR5_HV292.1        dis
## 79                       Bacillus            Chicken2_TR6_VAN90        dis
## 80                       Bacillus  Chicken2_TR7_VAN+CCUG5916890        dis
## 81                    Bacteroides    Chicken2_TR2_CTX+HV292.120     normal
## 82                    Bacteroides        Chicken2_TR4_CCUG59168     normal
## 83                    Bacteroides          Chicken2_TR5_HV292.1     normal
## 84                Bifidobacterium         Chicken2_CR_UNTREATED     normal
## 85                Bifidobacterium    Chicken2_TR2_CTX+HV292.120     normal
## 86                Bifidobacterium            Chicken2_TR3_CTX20     normal
## 87                Bifidobacterium        Chicken2_TR4_CCUG59168     normal
## 88                Bifidobacterium          Chicken2_TR5_HV292.1     normal
## 89                Bifidobacterium            Chicken2_TR6_VAN90     normal
## 90                Bifidobacterium  Chicken2_TR7_VAN+CCUG5916890     normal
## 91                        Blautia    Chicken1_TR1_CTX+HV292.120        mid
## 92                        Blautia            Chicken1_TR2_CTX20        mid
## 93                        Blautia          Chicken1_TR3_HV292.1     normal
## 94                        Blautia            Chicken1_TR4_VAN90     normal
## 95                        Blautia  Chicken1_TR5_VAN+CCUG5916890        dis
## 96                        Blautia        Chicken1_TR6_CCUG59168     normal
## 97                        Blautia         Chicken2_CR_UNTREATED     normal
## 98                        Blautia    Chicken2_TR2_CTX+HV292.120     normal
## 99                        Blautia            Chicken2_TR3_CTX20     normal
## 100                       Blautia        Chicken2_TR4_CCUG59168     normal
## 101                       Blautia          Chicken2_TR5_HV292.1     normal
## 102                       Blautia            Chicken2_TR6_VAN90        dis
## 103                       Blautia  Chicken2_TR7_VAN+CCUG5916890        dis
## 104                Butyricicoccus    Chicken1_TR1_CTX+HV292.120 recovering
## 105                Butyricicoccus            Chicken1_TR2_CTX20        mid
## 106                Butyricicoccus          Chicken1_TR3_HV292.1 recovering
## 107                Butyricicoccus        Chicken1_TR6_CCUG59168        dis
## 108                Butyricicoccus         Chicken2_CR_UNTREATED        mid
## 109                Butyricicoccus    Chicken2_TR2_CTX+HV292.120        mid
## 110                Butyricicoccus            Chicken2_TR3_CTX20     normal
## 111                Butyricicoccus          Chicken2_TR5_HV292.1        mid
## 112                Butyricicoccus  Chicken2_TR7_VAN+CCUG5916890        mid
## 113                        CAG-56    Chicken1_TR1_CTX+HV292.120        mid
## 114                        CAG-56            Chicken1_TR2_CTX20        mid
## 115                        CAG-56  Chicken1_TR5_VAN+CCUG5916890        dis
## 116                        CAG-56        Chicken1_TR6_CCUG59168        dis
## 117                        CAG-56    Chicken2_TR2_CTX+HV292.120        mid
## 118                        CAG-56            Chicken2_TR3_CTX20        mid
## 119                 Campylobacter         Chicken2_CR_UNTREATED        mid
## 120        Candidatus Soleaferrea    Chicken1_TR1_CTX+HV292.120     normal
## 121        Candidatus Soleaferrea            Chicken1_TR2_CTX20     normal
## 122        Candidatus Soleaferrea        Chicken1_TR6_CCUG59168 recovering
## 123        Candidatus Soleaferrea         Chicken2_CR_UNTREATED        mid
## 124        Candidatus Soleaferrea    Chicken2_TR2_CTX+HV292.120       prol
## 125        Candidatus Soleaferrea            Chicken2_TR3_CTX20     normal
## 126        Candidatus Soleaferrea        Chicken2_TR4_CCUG59168     normal
## 127        Candidatus Soleaferrea          Chicken2_TR5_HV292.1        mid
## 128               Catenibacterium         Chicken2_CR_UNTREATED recovering
## 129                      CHKCI001          Chicken1_TR3_HV292.1        mid
## 130                      CHKCI002    Chicken1_TR1_CTX+HV292.120     normal
## 131                      CHKCI002            Chicken1_TR2_CTX20     normal
## 132                      CHKCI002          Chicken1_TR3_HV292.1     normal
## 133                      CHKCI002            Chicken1_TR4_VAN90        dis
## 134                      CHKCI002  Chicken1_TR5_VAN+CCUG5916890        dis
## 135                      CHKCI002        Chicken1_TR6_CCUG59168     normal
## 136                      CHKCI002         Chicken2_CR_UNTREATED       prol
## 137                      CHKCI002    Chicken2_TR2_CTX+HV292.120       prol
## 138                      CHKCI002            Chicken2_TR3_CTX20     normal
## 139                      CHKCI002        Chicken2_TR4_CCUG59168       prol
## 140                      CHKCI002          Chicken2_TR5_HV292.1       prol
## 141 Christensenellaceae R-7 group         Chicken2_CR_UNTREATED        mid
## 142 Christensenellaceae R-7 group    Chicken2_TR2_CTX+HV292.120     normal
## 143 Christensenellaceae R-7 group        Chicken2_TR4_CCUG59168       prol
## 144 Christensenellaceae R-7 group          Chicken2_TR5_HV292.1 recovering
## 145                Clostridioides            Chicken2_TR6_VAN90        mid
## 146   Clostridium sensu stricto 1            Chicken2_TR6_VAN90        mid
## 147              Colidextribacter    Chicken1_TR1_CTX+HV292.120     normal
## 148              Colidextribacter            Chicken1_TR2_CTX20     normal
## 149              Colidextribacter          Chicken1_TR3_HV292.1     normal
## 150              Colidextribacter            Chicken1_TR4_VAN90        dis
## 151              Colidextribacter  Chicken1_TR5_VAN+CCUG5916890        dis
## 152              Colidextribacter        Chicken1_TR6_CCUG59168     normal
## 153              Colidextribacter         Chicken2_CR_UNTREATED     normal
## 154              Colidextribacter    Chicken2_TR2_CTX+HV292.120       prol
## 155              Colidextribacter            Chicken2_TR3_CTX20 recovering
## 156              Colidextribacter        Chicken2_TR4_CCUG59168       prol
## 157              Colidextribacter          Chicken2_TR5_HV292.1       prol
## 158                   Collinsella    Chicken2_TR2_CTX+HV292.120     normal
## 159                   Collinsella            Chicken2_TR3_CTX20        mid
## 160                   Collinsella        Chicken2_TR4_CCUG59168     normal
## 161                   Collinsella          Chicken2_TR5_HV292.1 recovering
## 162                   Collinsella            Chicken2_TR6_VAN90     normal
## 163                   Collinsella  Chicken2_TR7_VAN+CCUG5916890 recovering
## 164                   Coprococcus          Chicken2_TR5_HV292.1        mid
## 165                   Coprococcus  Chicken2_TR7_VAN+CCUG5916890        mid
## 166                     Cosenzaea            Chicken2_TR3_CTX20     normal
## 167     Defluviitaleaceae UCG-011    Chicken1_TR1_CTX+HV292.120 recovering
## 168     Defluviitaleaceae UCG-011            Chicken1_TR2_CTX20     normal
## 169     Defluviitaleaceae UCG-011          Chicken1_TR3_HV292.1       prol
## 170     Defluviitaleaceae UCG-011            Chicken1_TR4_VAN90        dis
## 171     Defluviitaleaceae UCG-011  Chicken1_TR5_VAN+CCUG5916890        dis
## 172     Defluviitaleaceae UCG-011        Chicken1_TR6_CCUG59168     normal
## 173     Defluviitaleaceae UCG-011         Chicken2_CR_UNTREATED       prol
## 174     Defluviitaleaceae UCG-011    Chicken2_TR2_CTX+HV292.120     normal
## 175     Defluviitaleaceae UCG-011            Chicken2_TR3_CTX20 recovering
## 176     Defluviitaleaceae UCG-011        Chicken2_TR4_CCUG59168       prol
## 177     Defluviitaleaceae UCG-011  Chicken2_TR7_VAN+CCUG5916890     normal
## 178                     Dialister         Chicken2_CR_UNTREATED     normal
## 179                     Dialister    Chicken2_TR2_CTX+HV292.120     normal
## 180                     Dialister            Chicken2_TR3_CTX20        mid
## 181                     Dialister          Chicken2_TR5_HV292.1     normal
## 182                     Dialister            Chicken2_TR6_VAN90 recovering
## 183                     Dialister  Chicken2_TR7_VAN+CCUG5916890        mid
## 184                         Dorea          Chicken2_TR5_HV292.1     normal
## 185                         Dorea            Chicken2_TR6_VAN90        mid
## 186                        DTU089    Chicken1_TR1_CTX+HV292.120     normal
## 187                        DTU089            Chicken1_TR2_CTX20        mid
## 188                        DTU089          Chicken1_TR3_HV292.1     normal
## 189                        DTU089         Chicken2_CR_UNTREATED       prol
## 190                        DTU089    Chicken2_TR2_CTX+HV292.120       prol
## 191                        DTU089        Chicken2_TR4_CCUG59168        mid
## 192                        DTU089          Chicken2_TR5_HV292.1        mid
## 193                Eisenbergiella    Chicken1_TR1_CTX+HV292.120     normal
## 194                Eisenbergiella            Chicken1_TR2_CTX20 recovering
## 195                Eisenbergiella          Chicken1_TR3_HV292.1     normal
## 196                Eisenbergiella            Chicken1_TR4_VAN90        dis
## 197                Eisenbergiella  Chicken1_TR5_VAN+CCUG5916890        dis
## 198                Eisenbergiella        Chicken1_TR6_CCUG59168     normal
## 199                Eisenbergiella         Chicken2_CR_UNTREATED     normal
## 200                Eisenbergiella    Chicken2_TR2_CTX+HV292.120     normal
## 201                Eisenbergiella        Chicken2_TR4_CCUG59168     normal
## 202                Eisenbergiella          Chicken2_TR5_HV292.1        dis
## 203                Eisenbergiella            Chicken2_TR6_VAN90        dis
## 204                  Enterococcus    Chicken1_TR1_CTX+HV292.120       prol
## 205                  Enterococcus            Chicken1_TR2_CTX20       prol
## 206                  Enterococcus  Chicken1_TR5_VAN+CCUG5916890     normal
## 207                  Enterococcus        Chicken1_TR6_CCUG59168     normal
## 208                  Enterococcus         Chicken2_CR_UNTREATED     normal
## 209                  Enterococcus    Chicken2_TR2_CTX+HV292.120       prol
## 210                  Enterococcus        Chicken2_TR4_CCUG59168        mid
## 211                  Enterococcus          Chicken2_TR5_HV292.1        mid
## 212                  Enterococcus            Chicken2_TR6_VAN90        dis
## 213                  Enterococcus  Chicken2_TR7_VAN+CCUG5916890     normal
## 214        Erysipelatoclostridium    Chicken1_TR1_CTX+HV292.120     normal
## 215        Erysipelatoclostridium            Chicken1_TR2_CTX20        dis
## 216        Erysipelatoclostridium          Chicken1_TR3_HV292.1     normal
## 217        Erysipelatoclostridium            Chicken1_TR4_VAN90        dis
## 218        Erysipelatoclostridium  Chicken1_TR5_VAN+CCUG5916890        dis
## 219        Erysipelatoclostridium        Chicken1_TR6_CCUG59168     normal
## 220        Erysipelatoclostridium    Chicken2_TR2_CTX+HV292.120        mid
## 221          Escherichia-Shigella    Chicken1_TR1_CTX+HV292.120       prol
## 222          Escherichia-Shigella            Chicken1_TR2_CTX20        dis
## 223          Escherichia-Shigella          Chicken1_TR3_HV292.1     normal
## 224          Escherichia-Shigella            Chicken1_TR4_VAN90 recovering
## 225          Escherichia-Shigella  Chicken1_TR5_VAN+CCUG5916890     normal
## 226          Escherichia-Shigella        Chicken1_TR6_CCUG59168       prol
## 227          Escherichia-Shigella         Chicken2_CR_UNTREATED        mid
## 228          Escherichia-Shigella    Chicken2_TR2_CTX+HV292.120        dis
## 229          Escherichia-Shigella            Chicken2_TR3_CTX20     normal
## 230          Escherichia-Shigella        Chicken2_TR4_CCUG59168        mid
## 231          Escherichia-Shigella  Chicken2_TR7_VAN+CCUG5916890        mid
## 232      Family XIII AD3011 group          Chicken1_TR3_HV292.1     normal
## 233      Family XIII AD3011 group            Chicken1_TR4_VAN90 recovering
## 234      Family XIII AD3011 group  Chicken1_TR5_VAN+CCUG5916890       prol
## 235      Family XIII AD3011 group        Chicken1_TR6_CCUG59168 recovering
## 236      Family XIII AD3011 group         Chicken2_CR_UNTREATED        mid
## 237      Family XIII AD3011 group    Chicken2_TR2_CTX+HV292.120       prol
## 238      Family XIII AD3011 group        Chicken2_TR4_CCUG59168       prol
## 239      Family XIII AD3011 group          Chicken2_TR5_HV292.1        mid
## 240      Family XIII AD3011 group            Chicken2_TR6_VAN90     normal
## 241      Family XIII AD3011 group  Chicken2_TR7_VAN+CCUG5916890 recovering
## 242                Flavonifractor    Chicken1_TR1_CTX+HV292.120        mid
## 243                Flavonifractor            Chicken1_TR2_CTX20       prol
## 244                Flavonifractor          Chicken1_TR3_HV292.1     normal
## 245                Flavonifractor            Chicken1_TR4_VAN90        dis
## 246                Flavonifractor  Chicken1_TR5_VAN+CCUG5916890        dis
## 247                Flavonifractor        Chicken1_TR6_CCUG59168 recovering
## 248                Flavonifractor         Chicken2_CR_UNTREATED recovering
## 249                Flavonifractor    Chicken2_TR2_CTX+HV292.120       prol
## 250                Flavonifractor            Chicken2_TR3_CTX20       prol
## 251                Flavonifractor        Chicken2_TR4_CCUG59168        mid
## 252                Flavonifractor          Chicken2_TR5_HV292.1 recovering
## 253                  Fournierella         Chicken2_CR_UNTREATED     normal
## 254                  Fournierella    Chicken2_TR2_CTX+HV292.120     normal
## 255                  Fournierella            Chicken2_TR3_CTX20     normal
## 256                  Fournierella        Chicken2_TR4_CCUG59168 recovering
## 257                  Fournierella          Chicken2_TR5_HV292.1        dis
## 258                Frisingicoccus          Chicken2_TR5_HV292.1        mid
## 259                 GCA-900066575    Chicken1_TR1_CTX+HV292.120 recovering
## 260                 GCA-900066575            Chicken1_TR2_CTX20        dis
## 261                 GCA-900066575          Chicken1_TR3_HV292.1     normal
## 262                 GCA-900066575            Chicken1_TR4_VAN90        dis
## 263                 GCA-900066575  Chicken1_TR5_VAN+CCUG5916890        dis
## 264                 GCA-900066575        Chicken1_TR6_CCUG59168     normal
## 265                 GCA-900066575         Chicken2_CR_UNTREATED       prol
## 266                 GCA-900066575    Chicken2_TR2_CTX+HV292.120       prol
## 267                 GCA-900066575            Chicken2_TR3_CTX20 recovering
## 268                 GCA-900066575        Chicken2_TR4_CCUG59168       prol
## 269                 GCA-900066575          Chicken2_TR5_HV292.1     normal
## 270                 GCA-900066575            Chicken2_TR6_VAN90        dis
## 271                 Gordonibacter    Chicken1_TR1_CTX+HV292.120        mid
## 272                 Gordonibacter          Chicken1_TR3_HV292.1     normal
## 273                 Gordonibacter            Chicken1_TR4_VAN90     normal
## 274                 Gordonibacter  Chicken1_TR5_VAN+CCUG5916890        dis
## 275                 Gordonibacter        Chicken1_TR6_CCUG59168     normal
## 276                 Gordonibacter         Chicken2_CR_UNTREATED     normal
## 277                 Gordonibacter        Chicken2_TR4_CCUG59168        dis
## 278                 Gordonibacter          Chicken2_TR5_HV292.1        mid
## 279                 Gordonibacter            Chicken2_TR6_VAN90        mid
## 280                  Harryflintia    Chicken1_TR1_CTX+HV292.120     normal
## 281                  Harryflintia            Chicken1_TR2_CTX20     normal
## 282                  Harryflintia          Chicken1_TR3_HV292.1 recovering
## 283                  Harryflintia  Chicken1_TR5_VAN+CCUG5916890        dis
## 284                  Harryflintia        Chicken1_TR6_CCUG59168     normal
## 285                  Harryflintia         Chicken2_CR_UNTREATED     normal
## 286                  Harryflintia            Chicken2_TR3_CTX20     normal
## 287                  Harryflintia        Chicken2_TR4_CCUG59168        mid
## 288                  Harryflintia          Chicken2_TR5_HV292.1        mid
## 289                    Holdemania          Chicken1_TR3_HV292.1        mid
## 290                    Holdemania            Chicken1_TR4_VAN90       prol
## 291                    Holdemania  Chicken1_TR5_VAN+CCUG5916890 recovering
## 292                    Holdemania        Chicken1_TR6_CCUG59168        mid
## 293                    Holdemania    Chicken2_TR2_CTX+HV292.120        mid
## 294                    Holdemania        Chicken2_TR4_CCUG59168        mid
## 295                         HT002          Chicken2_TR5_HV292.1        mid
## 296               Intestinibacter        Chicken2_TR4_CCUG59168        mid
## 297               Intestinibacter          Chicken2_TR5_HV292.1     normal
## 298                Intestinimonas    Chicken1_TR1_CTX+HV292.120     normal
## 299                Intestinimonas            Chicken1_TR2_CTX20 recovering
## 300                Intestinimonas          Chicken1_TR3_HV292.1        mid
## 301                Intestinimonas            Chicken1_TR4_VAN90        dis
## 302                Intestinimonas  Chicken1_TR5_VAN+CCUG5916890        dis
## 303                Intestinimonas        Chicken1_TR6_CCUG59168     normal
## 304                Intestinimonas         Chicken2_CR_UNTREATED     normal
## 305                Intestinimonas    Chicken2_TR2_CTX+HV292.120     normal
## 306                Intestinimonas            Chicken2_TR3_CTX20     normal
## 307                Intestinimonas        Chicken2_TR4_CCUG59168       prol
## 308                Intestinimonas          Chicken2_TR5_HV292.1     normal
## 309                Intestinimonas  Chicken2_TR7_VAN+CCUG5916890        dis
## 310             Lachnoclostridium    Chicken1_TR1_CTX+HV292.120        mid
## 311             Lachnoclostridium            Chicken1_TR2_CTX20     normal
## 312             Lachnoclostridium          Chicken1_TR3_HV292.1     normal
## 313             Lachnoclostridium            Chicken1_TR4_VAN90        dis
## 314             Lachnoclostridium  Chicken1_TR5_VAN+CCUG5916890        dis
## 315             Lachnoclostridium        Chicken1_TR6_CCUG59168     normal
## 316             Lachnoclostridium         Chicken2_CR_UNTREATED recovering
## 317             Lachnoclostridium    Chicken2_TR2_CTX+HV292.120     normal
## 318             Lachnoclostridium            Chicken2_TR3_CTX20 recovering
## 319             Lachnoclostridium        Chicken2_TR4_CCUG59168     normal
## 320             Lachnoclostridium          Chicken2_TR5_HV292.1       prol
## 321                   Lachnospira    Chicken1_TR1_CTX+HV292.120 recovering
## 322                   Lachnospira            Chicken1_TR2_CTX20        dis
## 323                   Lachnospira          Chicken1_TR3_HV292.1     normal
## 324                   Lachnospira            Chicken1_TR4_VAN90        dis
## 325                   Lachnospira  Chicken1_TR5_VAN+CCUG5916890        dis
## 326                   Lachnospira         Chicken2_CR_UNTREATED     normal
## 327                   Lachnospira    Chicken2_TR2_CTX+HV292.120       prol
## 328                   Lachnospira            Chicken2_TR3_CTX20 recovering
## 329                   Lachnospira        Chicken2_TR4_CCUG59168       prol
## 330                   Lachnospira          Chicken2_TR5_HV292.1       prol
## 331                   Lachnospira            Chicken2_TR6_VAN90        dis
## 332  Lachnospiraceae FCS020 group    Chicken1_TR1_CTX+HV292.120        mid
## 333  Lachnospiraceae FCS020 group            Chicken1_TR2_CTX20        mid
## 334  Lachnospiraceae FCS020 group          Chicken1_TR3_HV292.1        mid
## 335  Lachnospiraceae FCS020 group  Chicken1_TR5_VAN+CCUG5916890        dis
## 336  Lachnospiraceae FCS020 group        Chicken1_TR6_CCUG59168       prol
## 337  Lachnospiraceae NK3A20 group         Chicken2_CR_UNTREATED recovering
## 338  Lachnospiraceae NK3A20 group    Chicken2_TR2_CTX+HV292.120        mid
## 339  Lachnospiraceae NK3A20 group            Chicken2_TR3_CTX20        dis
## 340  Lachnospiraceae NK3A20 group        Chicken2_TR4_CCUG59168        dis
## 341  Lachnospiraceae NK3A20 group          Chicken2_TR5_HV292.1       prol
## 342  Lachnospiraceae NK3A20 group            Chicken2_TR6_VAN90     normal
## 343  Lachnospiraceae NK3A20 group  Chicken2_TR7_VAN+CCUG5916890       prol
## 344 Lachnospiraceae NK4A136 group    Chicken1_TR1_CTX+HV292.120     normal
## 345       Lachnospiraceae UCG-004    Chicken1_TR1_CTX+HV292.120     normal
## 346       Lachnospiraceae UCG-010    Chicken1_TR1_CTX+HV292.120 recovering
## 347       Lachnospiraceae UCG-010            Chicken1_TR2_CTX20        dis
## 348       Lachnospiraceae UCG-010          Chicken1_TR3_HV292.1     normal
## 349       Lachnospiraceae UCG-010            Chicken1_TR4_VAN90        dis
## 350       Lachnospiraceae UCG-010  Chicken1_TR5_VAN+CCUG5916890        dis
## 351       Lachnospiraceae UCG-010        Chicken1_TR6_CCUG59168        mid
## 352       Lachnospiraceae UCG-010         Chicken2_CR_UNTREATED       prol
## 353       Lachnospiraceae UCG-010    Chicken2_TR2_CTX+HV292.120     normal
## 354       Lachnospiraceae UCG-010            Chicken2_TR3_CTX20 recovering
## 355       Lachnospiraceae UCG-010          Chicken2_TR5_HV292.1     normal
## 356            Lacticaseibacillus         Chicken2_CR_UNTREATED       prol
## 357            Lacticaseibacillus    Chicken2_TR2_CTX+HV292.120     normal
## 358            Lacticaseibacillus            Chicken2_TR3_CTX20     normal
## 359            Lacticaseibacillus        Chicken2_TR4_CCUG59168       prol
## 360            Lacticaseibacillus          Chicken2_TR5_HV292.1       prol
## 361            Lacticaseibacillus            Chicken2_TR6_VAN90       prol
## 362            Lacticaseibacillus  Chicken2_TR7_VAN+CCUG5916890       prol
## 363           Lactiplantibacillus         Chicken2_CR_UNTREATED        mid
## 364           Lactiplantibacillus    Chicken2_TR2_CTX+HV292.120        mid
## 365           Lactiplantibacillus            Chicken2_TR3_CTX20        mid
## 366           Lactiplantibacillus  Chicken2_TR7_VAN+CCUG5916890     normal
## 367                 Lactobacillus    Chicken1_TR1_CTX+HV292.120        mid
## 368                 Lactobacillus  Chicken1_TR5_VAN+CCUG5916890        mid
## 369                   Lactococcus            Chicken2_TR3_CTX20     normal
## 370                   Lactococcus        Chicken2_TR4_CCUG59168        dis
## 371                   Lactococcus          Chicken2_TR5_HV292.1        mid
## 372                   Lactococcus            Chicken2_TR6_VAN90     normal
## 373                   Lactococcus  Chicken2_TR7_VAN+CCUG5916890 recovering
## 374             Latilactobacillus    Chicken2_TR2_CTX+HV292.120     normal
## 375             Latilactobacillus          Chicken2_TR5_HV292.1     normal
## 376             Latilactobacillus            Chicken2_TR6_VAN90        mid
## 377                   Leuconostoc         Chicken2_CR_UNTREATED     normal
## 378                   Leuconostoc    Chicken2_TR2_CTX+HV292.120     normal
## 379                   Leuconostoc            Chicken2_TR3_CTX20        mid
## 380                   Leuconostoc        Chicken2_TR4_CCUG59168        dis
## 381                   Leuconostoc          Chicken2_TR5_HV292.1     normal
## 382                   Leuconostoc  Chicken2_TR7_VAN+CCUG5916890     normal
## 383             Ligilactobacillus    Chicken1_TR1_CTX+HV292.120        mid
## 384             Ligilactobacillus            Chicken1_TR2_CTX20     normal
## 385             Ligilactobacillus            Chicken1_TR4_VAN90        mid
## 386             Ligilactobacillus  Chicken1_TR5_VAN+CCUG5916890        mid
## 387             Ligilactobacillus        Chicken1_TR6_CCUG59168     normal
## 388             Ligilactobacillus            Chicken2_TR3_CTX20        mid
## 389             Ligilactobacillus          Chicken2_TR5_HV292.1        mid
## 390             Ligilactobacillus            Chicken2_TR6_VAN90        mid
## 391             Ligilactobacillus  Chicken2_TR7_VAN+CCUG5916890        mid
## 392           Limosilactobacillus         Chicken2_CR_UNTREATED        mid
## 393           Limosilactobacillus  Chicken2_TR7_VAN+CCUG5916890        mid
## 394                      Listeria            Chicken2_TR6_VAN90        mid
## 395                Marvinbryantia    Chicken1_TR1_CTX+HV292.120        mid
## 396                Marvinbryantia            Chicken1_TR2_CTX20        mid
## 397                Marvinbryantia          Chicken1_TR3_HV292.1     normal
## 398                Marvinbryantia  Chicken1_TR5_VAN+CCUG5916890        dis
## 399                Marvinbryantia         Chicken2_CR_UNTREATED        dis
## 400                Marvinbryantia    Chicken2_TR2_CTX+HV292.120     normal
## 401                Marvinbryantia            Chicken2_TR3_CTX20        mid
## 402                Marvinbryantia        Chicken2_TR4_CCUG59168     normal
## 403                Marvinbryantia          Chicken2_TR5_HV292.1     normal
## 404                     Megamonas  Chicken2_TR7_VAN+CCUG5916890        mid
## 405                   Megasphaera        Chicken2_TR4_CCUG59168        mid
## 406                   Megasphaera            Chicken2_TR6_VAN90        mid
## 407                   Merdibacter  Chicken1_TR5_VAN+CCUG5916890        mid
## 408                   Merdibacter        Chicken1_TR6_CCUG59168        mid
## 409                   Merdibacter         Chicken2_CR_UNTREATED recovering
## 410                   Merdibacter    Chicken2_TR2_CTX+HV292.120        mid
## 411                   Merdibacter            Chicken2_TR3_CTX20     normal
## 412                   Merdibacter        Chicken2_TR4_CCUG59168       prol
## 413                   Merdibacter          Chicken2_TR5_HV292.1       prol
## 414                Methanosphaera        Chicken2_TR4_CCUG59168     normal
## 415                Methanosphaera            Chicken2_TR6_VAN90        mid
## 416                 Mogibacterium         Chicken2_CR_UNTREATED     normal
## 417                 Mogibacterium            Chicken2_TR6_VAN90        dis
## 418                    Monoglobus    Chicken1_TR1_CTX+HV292.120 recovering
## 419                    Monoglobus            Chicken1_TR2_CTX20        dis
## 420                    Monoglobus          Chicken1_TR3_HV292.1     normal
## 421                    Monoglobus            Chicken1_TR4_VAN90     normal
## 422                    Monoglobus  Chicken1_TR5_VAN+CCUG5916890        dis
## 423                    Monoglobus        Chicken1_TR6_CCUG59168     normal
## 424                    Monoglobus         Chicken2_CR_UNTREATED       prol
## 425                    Monoglobus            Chicken2_TR3_CTX20 recovering
## 426                    Monoglobus        Chicken2_TR4_CCUG59168     normal
## 427                    Monoglobus          Chicken2_TR5_HV292.1     normal
## 428                    Monoglobus            Chicken2_TR6_VAN90        dis
## 429              Negativibacillus    Chicken1_TR1_CTX+HV292.120 recovering
## 430              Negativibacillus            Chicken1_TR2_CTX20        dis
## 431              Negativibacillus          Chicken1_TR3_HV292.1     normal
## 432              Negativibacillus            Chicken1_TR4_VAN90        dis
## 433              Negativibacillus  Chicken1_TR5_VAN+CCUG5916890        dis
## 434              Negativibacillus        Chicken1_TR6_CCUG59168     normal
## 435              Negativibacillus         Chicken2_CR_UNTREATED     normal
## 436              Negativibacillus    Chicken2_TR2_CTX+HV292.120     normal
## 437              Negativibacillus            Chicken2_TR3_CTX20 recovering
## 438              Negativibacillus        Chicken2_TR4_CCUG59168        mid
## 439              Negativibacillus          Chicken2_TR5_HV292.1     normal
## 440              Negativibacillus  Chicken2_TR7_VAN+CCUG5916890        dis
## 441                 NK4A214 group    Chicken1_TR1_CTX+HV292.120        mid
## 442                 NK4A214 group            Chicken1_TR2_CTX20     normal
## 443                 NK4A214 group          Chicken1_TR3_HV292.1     normal
## 444                 NK4A214 group            Chicken1_TR4_VAN90        dis
## 445                 NK4A214 group  Chicken1_TR5_VAN+CCUG5916890        dis
## 446                 NK4A214 group        Chicken1_TR6_CCUG59168        mid
## 447                 NK4A214 group         Chicken2_CR_UNTREATED       prol
## 448                 NK4A214 group    Chicken2_TR2_CTX+HV292.120        mid
## 449                 NK4A214 group            Chicken2_TR3_CTX20     normal
## 450                 NK4A214 group        Chicken2_TR4_CCUG59168     normal
## 451                 NK4A214 group          Chicken2_TR5_HV292.1     normal
## 452                 NK4A214 group            Chicken2_TR6_VAN90        dis
## 453                Oceanobacillus    Chicken2_TR2_CTX+HV292.120 recovering
## 454                 Oscillibacter    Chicken1_TR1_CTX+HV292.120     normal
## 455                 Oscillibacter            Chicken1_TR2_CTX20 recovering
## 456                 Oscillibacter          Chicken1_TR3_HV292.1     normal
## 457                 Oscillibacter            Chicken1_TR4_VAN90        dis
## 458                 Oscillibacter  Chicken1_TR5_VAN+CCUG5916890        dis
## 459                 Oscillibacter        Chicken1_TR6_CCUG59168     normal
## 460                 Oscillibacter         Chicken2_CR_UNTREATED     normal
## 461                 Oscillibacter    Chicken2_TR2_CTX+HV292.120     normal
## 462                 Oscillibacter            Chicken2_TR3_CTX20        mid
## 463                 Oscillibacter        Chicken2_TR4_CCUG59168     normal
## 464                 Oscillibacter          Chicken2_TR5_HV292.1     normal
## 465                  Oscillospira    Chicken1_TR1_CTX+HV292.120     normal
## 466                  Oscillospira            Chicken1_TR2_CTX20        dis
## 467                  Oscillospira          Chicken1_TR3_HV292.1 recovering
## 468                  Oscillospira            Chicken1_TR4_VAN90     normal
## 469                  Oscillospira  Chicken1_TR5_VAN+CCUG5916890        dis
## 470                  Oscillospira        Chicken1_TR6_CCUG59168 recovering
## 471                  Oscillospira         Chicken2_CR_UNTREATED     normal
## 472                  Oscillospira    Chicken2_TR2_CTX+HV292.120        mid
## 473                  Oscillospira            Chicken2_TR3_CTX20 recovering
## 474                  Oscillospira        Chicken2_TR4_CCUG59168        mid
## 475                  Oscillospira          Chicken2_TR5_HV292.1        mid
## 476                 Paenibacillus    Chicken2_TR2_CTX+HV292.120     normal
## 477                    Paludicola    Chicken1_TR1_CTX+HV292.120        mid
## 478                    Paludicola  Chicken1_TR5_VAN+CCUG5916890        dis
## 479               Parabacteroides         Chicken2_CR_UNTREATED recovering
## 480               Parabacteroides    Chicken2_TR2_CTX+HV292.120 recovering
## 481               Parabacteroides            Chicken2_TR3_CTX20 recovering
## 482               Parabacteroides        Chicken2_TR4_CCUG59168       prol
## 483               Parabacteroides          Chicken2_TR5_HV292.1     normal
## 484                Parasutterella         Chicken2_CR_UNTREATED     normal
## 485                Parasutterella    Chicken2_TR2_CTX+HV292.120     normal
## 486                Parasutterella        Chicken2_TR4_CCUG59168     normal
## 487                Parasutterella          Chicken2_TR5_HV292.1       prol
## 488                Parasutterella            Chicken2_TR6_VAN90       prol
## 489                Parasutterella  Chicken2_TR7_VAN+CCUG5916890       prol
## 490                   Pediococcus    Chicken2_TR2_CTX+HV292.120        mid
## 491                   Pediococcus          Chicken2_TR5_HV292.1        mid
## 492                   Peptococcus            Chicken2_TR6_VAN90        mid
## 493                    Prevotella          Chicken2_TR5_HV292.1        mid
## 494                  Prevotella_9         Chicken2_CR_UNTREATED recovering
## 495                  Prevotella_9    Chicken2_TR2_CTX+HV292.120     normal
## 496                  Prevotella_9            Chicken2_TR3_CTX20        mid
## 497                  Prevotella_9        Chicken2_TR4_CCUG59168 recovering
## 498                  Prevotella_9          Chicken2_TR5_HV292.1 recovering
## 499                  Prevotella_9            Chicken2_TR6_VAN90        dis
## 500                  Prevotella_9  Chicken2_TR7_VAN+CCUG5916890        dis
## 501                       Proteus    Chicken1_TR1_CTX+HV292.120 recovering
## 502                       Proteus            Chicken1_TR2_CTX20        dis
## 503                       Proteus          Chicken1_TR3_HV292.1     normal
## 504                       Proteus            Chicken1_TR4_VAN90     normal
## 505                       Proteus  Chicken1_TR5_VAN+CCUG5916890     normal
## 506                       Proteus        Chicken1_TR6_CCUG59168     normal
## 507                       Proteus         Chicken2_CR_UNTREATED     normal
## 508                       Proteus    Chicken2_TR2_CTX+HV292.120 recovering
## 509                       Proteus        Chicken2_TR4_CCUG59168     normal
## 510                       Proteus          Chicken2_TR5_HV292.1     normal
## 511                       Proteus            Chicken2_TR6_VAN90        mid
## 512                       Proteus  Chicken2_TR7_VAN+CCUG5916890       prol
## 513          Pseudoflavonifractor            Chicken1_TR2_CTX20       prol
## 514          Pseudoflavonifractor          Chicken1_TR3_HV292.1       prol
## 515          Pseudoflavonifractor        Chicken1_TR6_CCUG59168     normal
## 516          Pseudoflavonifractor         Chicken2_CR_UNTREATED recovering
## 517          Pseudoflavonifractor    Chicken2_TR2_CTX+HV292.120 recovering
## 518          Pseudoflavonifractor            Chicken2_TR3_CTX20        mid
## 519          Pseudoflavonifractor        Chicken2_TR4_CCUG59168 recovering
## 520          Pseudoflavonifractor          Chicken2_TR5_HV292.1 recovering
## 521                   Pseudomonas    Chicken1_TR1_CTX+HV292.120     normal
## 522                   Pseudomonas            Chicken1_TR2_CTX20        dis
## 523                   Pseudomonas            Chicken1_TR4_VAN90     normal
## 524                   Pseudomonas  Chicken1_TR5_VAN+CCUG5916890     normal
## 525                   Pseudomonas        Chicken1_TR6_CCUG59168        mid
## 526                   Pseudomonas         Chicken2_CR_UNTREATED        mid
## 527                   Pseudomonas            Chicken2_TR3_CTX20     normal
## 528                   Pseudomonas            Chicken2_TR6_VAN90        mid
## 529                 Pygmaiobacter    Chicken1_TR1_CTX+HV292.120        mid
## 530                 Pygmaiobacter            Chicken1_TR2_CTX20       prol
## 531                 Pygmaiobacter        Chicken1_TR6_CCUG59168        mid
## 532                 Pygmaiobacter         Chicken2_CR_UNTREATED       prol
## 533                 Pygmaiobacter    Chicken2_TR2_CTX+HV292.120       prol
## 534                 Pygmaiobacter            Chicken2_TR3_CTX20        mid
## 535                 Pygmaiobacter        Chicken2_TR4_CCUG59168       prol
## 536                 Pygmaiobacter          Chicken2_TR5_HV292.1     normal
## 537                     Roseburia            Chicken2_TR3_CTX20        mid
## 538                     Roseburia        Chicken2_TR4_CCUG59168        mid
## 539                     Roseburia            Chicken2_TR6_VAN90     normal
## 540                     Roseburia  Chicken2_TR7_VAN+CCUG5916890     normal
## 541                  Ruminococcus    Chicken1_TR1_CTX+HV292.120     normal
## 542                  Ruminococcus            Chicken1_TR2_CTX20        mid
## 543                  Ruminococcus  Chicken1_TR5_VAN+CCUG5916890     normal
## 544                  Ruminococcus        Chicken1_TR6_CCUG59168     normal
## 545                  Ruminococcus         Chicken2_CR_UNTREATED        mid
## 546                  Ruminococcus    Chicken2_TR2_CTX+HV292.120        mid
## 547                  Ruminococcus            Chicken2_TR3_CTX20 recovering
## 548                  Ruminococcus        Chicken2_TR4_CCUG59168        mid
## 549                  Ruminococcus          Chicken2_TR5_HV292.1     normal
## 550                    Salmonella            Chicken2_TR6_VAN90        mid
## 551                       Sarcina  Chicken2_TR7_VAN+CCUG5916890        mid
## 552                    Sellimonas    Chicken1_TR1_CTX+HV292.120        mid
## 553                    Sellimonas            Chicken1_TR2_CTX20        mid
## 554                    Sellimonas          Chicken1_TR3_HV292.1     normal
## 555                    Sellimonas            Chicken1_TR4_VAN90        dis
## 556                    Sellimonas  Chicken1_TR5_VAN+CCUG5916890        dis
## 557                    Sellimonas        Chicken1_TR6_CCUG59168     normal
## 558                    Sellimonas         Chicken2_CR_UNTREATED       prol
## 559                    Sellimonas    Chicken2_TR2_CTX+HV292.120     normal
## 560                    Sellimonas            Chicken2_TR3_CTX20        mid
## 561                    Sellimonas        Chicken2_TR4_CCUG59168     normal
## 562                    Sellimonas          Chicken2_TR5_HV292.1     normal
## 563                Shuttleworthia    Chicken1_TR1_CTX+HV292.120 recovering
## 564                Shuttleworthia            Chicken1_TR2_CTX20        dis
## 565                Shuttleworthia          Chicken1_TR3_HV292.1     normal
## 566                Shuttleworthia            Chicken1_TR4_VAN90     normal
## 567                Shuttleworthia  Chicken1_TR5_VAN+CCUG5916890        dis
## 568                Shuttleworthia        Chicken1_TR6_CCUG59168     normal
## 569                Shuttleworthia         Chicken2_CR_UNTREATED       prol
## 570                Shuttleworthia    Chicken2_TR2_CTX+HV292.120       prol
## 571                Shuttleworthia            Chicken2_TR3_CTX20     normal
## 572                Shuttleworthia        Chicken2_TR4_CCUG59168 recovering
## 573                Shuttleworthia          Chicken2_TR5_HV292.1        mid
## 574                Staphylococcus         Chicken2_CR_UNTREATED        mid
## 575                Staphylococcus    Chicken2_TR2_CTX+HV292.120        mid
## 576                Staphylococcus            Chicken2_TR3_CTX20     normal
## 577                Staphylococcus        Chicken2_TR4_CCUG59168        dis
## 578                Staphylococcus          Chicken2_TR5_HV292.1        dis
## 579                Staphylococcus            Chicken2_TR6_VAN90        dis
## 580                Staphylococcus  Chicken2_TR7_VAN+CCUG5916890        dis
## 581              Stenotrophomonas            Chicken1_TR2_CTX20        mid
## 582              Stenotrophomonas            Chicken1_TR4_VAN90       prol
## 583              Stenotrophomonas  Chicken1_TR5_VAN+CCUG5916890     normal
## 584              Stenotrophomonas        Chicken1_TR6_CCUG59168       prol
## 585              Stenotrophomonas            Chicken2_TR6_VAN90        mid
## 586                 Streptococcus         Chicken2_CR_UNTREATED        mid
## 587                 Streptococcus    Chicken2_TR2_CTX+HV292.120     normal
## 588                 Streptococcus            Chicken2_TR3_CTX20        mid
## 589                 Streptococcus  Chicken2_TR7_VAN+CCUG5916890 recovering
## 590               Subdoligranulum            Chicken1_TR2_CTX20        mid
## 591               Subdoligranulum         Chicken2_CR_UNTREATED        dis
## 592               Subdoligranulum    Chicken2_TR2_CTX+HV292.120        dis
## 593               Subdoligranulum            Chicken2_TR3_CTX20        dis
## 594               Subdoligranulum        Chicken2_TR4_CCUG59168        dis
## 595               Subdoligranulum          Chicken2_TR5_HV292.1        dis
## 596               Subdoligranulum            Chicken2_TR6_VAN90        dis
## 597               Subdoligranulum  Chicken2_TR7_VAN+CCUG5916890        dis
## 598                    Sutterella         Chicken2_CR_UNTREATED        mid
## 599                    Sutterella            Chicken2_TR3_CTX20        mid
## 600                    Tuzzerella            Chicken1_TR2_CTX20     normal
## 601                    Tuzzerella          Chicken1_TR3_HV292.1        mid
## 602                    Tuzzerella            Chicken1_TR4_VAN90       prol
## 603                    Tuzzerella  Chicken1_TR5_VAN+CCUG5916890       prol
## 604                    Tuzzerella        Chicken1_TR6_CCUG59168     normal
## 605                    Tuzzerella         Chicken2_CR_UNTREATED     normal
## 606                    Tuzzerella    Chicken2_TR2_CTX+HV292.120     normal
## 607                    Tuzzerella            Chicken2_TR3_CTX20 recovering
## 608                    Tuzzerella        Chicken2_TR4_CCUG59168       prol
## 609                    Tuzzerella          Chicken2_TR5_HV292.1     normal
## 610                    Tuzzerella            Chicken2_TR6_VAN90       prol
## 611                    Tuzzerella  Chicken2_TR7_VAN+CCUG5916890       prol
## 612                    Tyzzerella    Chicken1_TR1_CTX+HV292.120 recovering
## 613                    Tyzzerella            Chicken1_TR2_CTX20        dis
## 614                    Tyzzerella          Chicken1_TR3_HV292.1     normal
## 615                    Tyzzerella            Chicken1_TR4_VAN90        dis
## 616                    Tyzzerella  Chicken1_TR5_VAN+CCUG5916890        dis
## 617                    Tyzzerella        Chicken1_TR6_CCUG59168     normal
## 618                    Tyzzerella         Chicken2_CR_UNTREATED     normal
## 619                    Tyzzerella    Chicken2_TR2_CTX+HV292.120 recovering
## 620                    Tyzzerella            Chicken2_TR3_CTX20 recovering
## 621                    Tyzzerella        Chicken2_TR4_CCUG59168     normal
## 622                    Tyzzerella          Chicken2_TR5_HV292.1 recovering
## 623                       UBA1819    Chicken1_TR1_CTX+HV292.120        dis
## 624                       UBA1819            Chicken1_TR2_CTX20 recovering
## 625                       UBA1819          Chicken1_TR3_HV292.1     normal
## 626                       UBA1819            Chicken1_TR4_VAN90        dis
## 627                       UBA1819  Chicken1_TR5_VAN+CCUG5916890        dis
## 628                       UBA1819        Chicken1_TR6_CCUG59168        mid
## 629                       UBA1819         Chicken2_CR_UNTREATED recovering
## 630                       UBA1819    Chicken2_TR2_CTX+HV292.120        mid
## 631                       UBA1819          Chicken2_TR5_HV292.1       prol
## 632                     UC5-1-2E3    Chicken1_TR1_CTX+HV292.120     normal
## 633                     UC5-1-2E3          Chicken1_TR3_HV292.1        mid
## 634                     UC5-1-2E3            Chicken1_TR4_VAN90        dis
## 635                     UC5-1-2E3  Chicken1_TR5_VAN+CCUG5916890        dis
## 636                     UC5-1-2E3        Chicken1_TR6_CCUG59168       prol
## 637                     UC5-1-2E3         Chicken2_CR_UNTREATED       prol
## 638                     UC5-1-2E3            Chicken2_TR3_CTX20     normal
## 639                     UC5-1-2E3        Chicken2_TR4_CCUG59168     normal
## 640                     UC5-1-2E3          Chicken2_TR5_HV292.1     normal
## 641                     UC5-1-2E3            Chicken2_TR6_VAN90        dis
## 642                       UCG-002          Chicken2_TR5_HV292.1        mid
## 643                       UCG-002  Chicken2_TR7_VAN+CCUG5916890     normal
## 644                       UCG-005    Chicken1_TR1_CTX+HV292.120        mid
## 645                       UCG-005            Chicken1_TR2_CTX20        mid
## 646                       UCG-005         Chicken2_CR_UNTREATED     normal
## 647                       UCG-005    Chicken2_TR2_CTX+HV292.120     normal
## 648                       UCG-005            Chicken2_TR3_CTX20     normal
## 649                       UCG-005        Chicken2_TR4_CCUG59168     normal
## 650                       UCG-005          Chicken2_TR5_HV292.1     normal
## 651                   Veillonella    Chicken2_TR2_CTX+HV292.120        mid
## 652                   Veillonella        Chicken2_TR4_CCUG59168     normal
## 653                   Veillonella          Chicken2_TR5_HV292.1        mid
## 654                   Veillonella            Chicken2_TR6_VAN90     normal
## 655                   Veillonella  Chicken2_TR7_VAN+CCUG5916890     normal
##     cat_lm pval_Box.test pval_adf.test pval_kpss.test_level
## 1   normal            NA    0.74726657                   NA
## 2   normal  0.2844056783    0.28627206           0.08474662
## 3   normal  0.0038134820    0.95815202           0.05610837
## 4   normal  0.1605252199    0.67663673           0.10000000
## 5   normal  0.2434371857    0.52735403           0.08314014
## 6   normal  0.9412046701    0.74944955           0.10000000
## 7   normal  0.9086710153    0.51114645           0.10000000
## 8   normal  0.0385215075    0.01000000           0.10000000
## 9   normal  0.0042094789    0.33931331           0.10000000
## 10  normal  0.5723941548    0.53374115           0.10000000
## 11     dec            NA    0.84829743                   NA
## 12     dec            NA    0.73758058                   NA
## 13  normal  0.0635628631    0.49092913           0.10000000
## 14     inc  0.0006155906    0.33745332           0.08079847
## 15     inc  0.0083690798    0.80386038           0.04274626
## 16       V  0.0298757880    0.95043687           0.10000000
## 17     inc  0.0005451588    0.81399800           0.03270626
## 18  normal  0.0057728996    0.99000000           0.10000000
## 19     dec            NA    0.39276928                   NA
## 20  normal  0.9107238391    0.59980079           0.10000000
## 21  normal  0.2138219915    0.34250768           0.10000000
## 22  normal  0.3484459962    0.81130794           0.09199948
## 23     inc  0.0315571114    0.63196036           0.05565873
## 24     inc  0.0087981624    0.72387328           0.05988155
## 25  normal  0.2672333212    0.58469584           0.05198515
## 26  normal  0.2506215163    0.42769687           0.10000000
## 27  normal  0.0856001831    0.32031009           0.10000000
## 28  normal            NA    0.99000000                   NA
## 29  normal  0.0418134353    0.06499288           0.10000000
## 30  normal            NA    0.79286321                   NA
## 31     inc  0.0001230407    0.49676631           0.03662666
## 32  normal  0.0147180766    0.01000000           0.10000000
## 33  normal  0.0252978432    0.83927641           0.10000000
## 34  normal  0.7813826929    0.91510334           0.10000000
## 35     dec            NA    0.93043909                   NA
## 36     dec  0.0326048810    0.92132714           0.04959409
## 37     dec  0.5621626819    0.87998350           0.06544864
## 38     inc  0.0005428576    0.55833675           0.02251949
## 39  normal            NA    0.42792266                   NA
## 40     dec  0.0008320392    0.99000000           0.06010030
## 41     dec  0.0024985909    0.34202533           0.05728038
## 42     dec            NA    0.15709453                   NA
## 43     inc  0.0026179384    0.03010863           0.04288902
## 44     inc  0.0634619495    0.09251437           0.04821335
## 45     inc  0.0082899332    0.45526381           0.03855040
## 46  normal            NA    0.52489122                   NA
## 47  normal            NA    0.77087665                   NA
## 48  normal  0.0875877729    0.30142768           0.10000000
## 49     inc            NA    0.85618712                   NA
## 50  normal            NA    0.90743977                   NA
## 51  normal            NA    0.71231436                   NA
## 52  normal  0.0442944392    0.99000000           0.10000000
## 53  normal  0.7912144692    0.30355932           0.10000000
## 54  normal  0.0771828443    0.15762396           0.10000000
## 55     inc  0.1175809662    0.79499591           0.06854252
## 56     inc  0.0007228611    0.75238882           0.08204737
## 57  normal            NA    0.41730258                   NA
## 58  normal            NA    0.08064532                   NA
## 59  normal  0.1757221180    0.52142089           0.10000000
## 60  normal            NA    0.72502608                   NA
## 61  normal            NA    0.21406165                   NA
## 62  normal            NA    0.58566603                   NA
## 63  normal            NA    0.90453152                   NA
## 64  normal            NA    0.52764806                   NA
## 65  normal  0.6576502873    0.62863375           0.10000000
## 66  normal  0.6020459205    0.78574784           0.10000000
## 67  normal  0.4555990455    0.12407804           0.10000000
## 68  normal  0.0207982364    0.95897871           0.10000000
## 69  normal  0.1148909624    0.81213895           0.10000000
## 70  normal  0.7724790541    0.71204840           0.10000000
## 71     inV  0.0084132796    0.94918553           0.09607153
## 72  normal            NA    0.84765347                   NA
## 73  normal            NA    0.36541612                   NA
## 74     dec  0.0164080901    0.43458833           0.07681060
## 75  normal  0.6504699589    0.84713876           0.10000000
## 76  normal  0.7992314589    0.55088024           0.10000000
## 77     dec  0.0007618246    0.69315507           0.05131019
## 78     dec  0.0024263587    0.92608796           0.02060902
## 79     dec            NA    0.75570622                   NA
## 80     dec            NA    0.73375260                   NA
## 81  normal  0.5402158278    0.26224797           0.07539558
## 82     dec  0.0688944095    0.96576244           0.10000000
## 83     dec  0.0163099374    0.44999173           0.04802051
## 84  normal  0.7196356248    0.50858140           0.10000000
## 85  normal  0.8520072568    0.12419839           0.10000000
## 86  normal  0.4981157204    0.73478124           0.10000000
## 87  normal  0.2559003679    0.72300976           0.10000000
## 88  normal  0.7256745856    0.56876052           0.10000000
## 89  normal  0.1423164749    0.18994402           0.10000000
## 90  normal  0.5066827106    0.25655894           0.10000000
## 91     dec  0.4089203158    0.80440447           0.10000000
## 92  normal  0.1276261979    0.99000000           0.10000000
## 93  normal  0.0520009040    0.35037563           0.10000000
## 94     dec  0.0023508983    0.48584235           0.06390673
## 95     dec  0.0069482980    0.62786389           0.05660343
## 96  normal  0.6530385110    0.37938763           0.10000000
## 97     inc  0.9915067764    0.39530246           0.03745193
## 98  normal  0.4068271472    0.08589443           0.10000000
## 99  normal  0.1158769200    0.95064858           0.10000000
## 100 normal  0.0224861923    0.94147449           0.10000000
## 101 normal  0.3263005463    0.31135690           0.10000000
## 102    dec            NA    0.65170548                   NA
## 103    dec            NA    0.06385438                   NA
## 104 normal  0.7143838625    0.46291499           0.10000000
## 105 normal  0.2439257960    0.18151633           0.10000000
## 106 normal  0.3076607602    0.20239971           0.10000000
## 107 normal  0.1870526021    0.69997414           0.05682103
## 108 normal            NA    0.52383811                   NA
## 109 normal            NA    0.31251766                   NA
## 110 normal            NA    0.62076831                   NA
## 111 normal            NA    0.82613656                   NA
## 112 normal            NA    0.60245407                   NA
## 113 normal  0.2190098092    0.49567567           0.10000000
## 114    inV  0.0028984581    0.12514508           0.10000000
## 115 normal            NA    0.13950060                   NA
## 116 normal            NA    0.85166014                   NA
## 117 normal            NA    0.22261415                   NA
## 118 normal            NA    0.45538811                   NA
## 119 normal            NA    0.57000977                   NA
## 120 normal            NA    0.99000000                   NA
## 121 normal  0.9487449321    0.91496178           0.10000000
## 122 normal            NA    0.99000000                   NA
## 123 normal            NA    0.57000977                   NA
## 124 normal  0.6604620941    0.12792332           0.10000000
## 125 normal  0.9539760442    0.47449124           0.09448715
## 126 normal  0.4276019409    0.72123625           0.10000000
## 127 normal  0.0023881091    0.48005350           0.10000000
## 128 normal            NA    0.29256745                   NA
## 129 normal            NA    0.24709685                   NA
## 130 normal  0.2660888518    0.41714418           0.10000000
## 131 normal  0.0593350984    0.82404223           0.10000000
## 132 normal  0.7589740694    0.99000000           0.10000000
## 133    dec  0.0040755841    0.70845247           0.06239239
## 134    dec  0.0067578080    0.98169216           0.04940961
## 135    dec  0.5675417322    0.08004706           0.10000000
## 136 normal            NA    0.74150013                   NA
## 137 normal            NA    0.38829559                   NA
## 138 normal            NA    0.44030280                   NA
## 139 normal            NA    0.62654519                   NA
## 140 normal            NA    0.64006934                   NA
## 141 normal            NA    0.25900801                   NA
## 142 normal            NA    0.23322176                   NA
## 143 normal  0.0010215266    0.98191826           0.05451367
## 144 normal            NA    0.73482225                   NA
## 145 normal            NA    0.52764806                   NA
## 146 normal            NA    0.38116763                   NA
## 147      V  0.0382621275    0.91800672           0.10000000
## 148 normal  0.0035354679    0.82956198           0.10000000
## 149 normal  0.4996465770    0.87705958           0.10000000
## 150    dec  0.0019122679    0.65535042           0.06672995
## 151    dec  0.0048083873    0.87317111           0.05285419
## 152 normal  0.7484914718    0.26201770           0.10000000
## 153 normal  0.5592960613    0.82396353           0.10000000
## 154    inc  0.0037195089    0.20005265           0.03515230
## 155 normal  0.0015500546    0.01000000           0.10000000
## 156    inc  0.0014021528    0.28949016           0.03773141
## 157    inc  0.0169479602    0.91171762           0.08262488
## 158 normal            NA    0.34732802                   NA
## 159 normal            NA    0.13579879                   NA
## 160 normal  0.9744497229    0.22089011           0.10000000
## 161 normal  0.7125823982    0.52949515           0.08549275
## 162 normal            NA    0.48987006                   NA
## 163 normal  0.4556703063    0.31003675           0.10000000
## 164 normal            NA    0.60245407                   NA
## 165 normal            NA    0.26368698                   NA
## 166 normal            NA    0.08659876                   NA
## 167 normal            NA    0.64948036                   NA
## 168 normal            NA    0.89827483                   NA
## 169 normal  0.0133003628    0.39658165           0.10000000
## 170    dec            NA    0.11266526                   NA
## 171 normal            NA    0.23427347                   NA
## 172 normal  0.7194361730    0.95410732           0.10000000
## 173    inc  0.0014603519    0.01000000           0.06086501
## 174 normal  0.0646481595    0.57468532           0.10000000
## 175    inc  0.0003199737    0.04914529           0.03866277
## 176    inc  0.1325135156    0.45196933           0.06343291
## 177 normal            NA    0.15286084                   NA
## 178 normal            NA    0.72462027                   NA
## 179 normal  0.7676779512    0.95329681           0.10000000
## 180 normal            NA    0.40445886                   NA
## 181 normal            NA    0.18528960                   NA
## 182 normal            NA    0.91177472                   NA
## 183 normal            NA    0.06596397                   NA
## 184 normal            NA    0.52199207                   NA
## 185 normal            NA    0.31769191                   NA
## 186 normal            NA    0.99000000                   NA
## 187 normal            NA    0.07221275                   NA
## 188 normal            NA    0.99000000                   NA
## 189 normal  0.0048286562    0.40331345           0.03485981
## 190    inc  0.0179126729    0.77674837           0.05624475
## 191 normal            NA    0.92347238                   NA
## 192 normal  0.2081011292    0.76326755           0.10000000
## 193 normal  0.6565459583    0.92806286           0.10000000
## 194    dec            NA    0.99000000                   NA
## 195 normal  0.2505325930    0.12877145           0.10000000
## 196    dec            NA    0.90556132                   NA
## 197    dec            NA    0.43208326                   NA
## 198 normal  0.0081087982    0.86051613           0.09996393
## 199    dec  0.0005087656    0.69028127           0.02443433
## 200    dec  0.0080057187    0.97239116           0.04243036
## 201    dec  0.0015163134    0.99000000           0.03737530
## 202    dec  0.0016185167    0.89642311           0.03160269
## 203    dec            NA    0.08740537                   NA
## 204    inc  0.1175176930    0.11875268           0.05342615
## 205    inc  0.0017574158    0.93667501           0.03547012
## 206 normal  0.4670845626    0.19982801           0.10000000
## 207 normal  0.0218723390    0.37005947           0.10000000
## 208 normal  0.0490126623    0.01000000           0.10000000
## 209    inc  0.0108659766    0.35276835           0.03782470
## 210    dec            NA    0.56144773                   NA
## 211 normal            NA    0.08630216                   NA
## 212 normal            NA    0.55513963                   NA
## 213    inc  0.3781724214    0.51163887           0.07261957
## 214    dec  0.3575304076    0.94216443           0.10000000
## 215    dec  0.0007596450    0.81290733           0.06163765
## 216    dec  0.1994249914    0.99000000           0.06369676
## 217    dec  0.0073917674    0.84867036           0.05951676
## 218    dec  0.1851018298    0.91351021           0.04773840
## 219    dec  0.0429197513    0.98003209           0.04772115
## 220 normal            NA    0.38116763                   NA
## 221 normal            NA    0.87984844                   NA
## 222 normal            NA    0.92124740                   NA
## 223 normal  0.3991577115    0.78107490           0.10000000
## 224    inc  0.1131256339    0.99000000           0.10000000
## 225    inc  0.1473393723    0.89491993           0.10000000
## 226 normal  0.0097066346    0.01000000           0.10000000
## 227    dec  0.0173321193    0.98683911           0.10000000
## 228    dec            NA    0.88037736                   NA
## 229 normal            NA    0.43903555                   NA
## 230    dec  0.0027215573    0.92514776           0.05744756
## 231    dec  0.0001569278    0.47303054           0.06563750
## 232 normal            NA    0.99000000                   NA
## 233    inc            NA    0.74000625                   NA
## 234    inc  0.0365982863    0.92985556           0.05358204
## 235 normal            NA    0.99000000                   NA
## 236 normal            NA    0.72843115                   NA
## 237 normal  0.3635899924    0.07932894           0.04527051
## 238 normal            NA    0.17575216                   NA
## 239 normal            NA    0.63306415                   NA
## 240    inc  0.0003264437    0.41373200           0.02542534
## 241    inc  0.0005355353    0.28407546           0.03905890
## 242    dec  0.5396595374    0.05593565           0.10000000
## 243 normal  0.0173601258    0.78850147           0.09775636
## 244 normal  0.0755876181    0.67508496           0.10000000
## 245    dec  0.0063772391    0.77406232           0.06474994
## 246    dec  0.0293691569    0.95046873           0.04821448
## 247    dec  0.0331457213    0.52647938           0.10000000
## 248 normal  0.6272968472    0.39289011           0.10000000
## 249    inc  0.0395698423    0.78195080           0.03991153
## 250    inc  0.0210375465    0.01249889           0.02971846
## 251    dec  0.0497458268    0.99000000           0.10000000
## 252 normal  0.1151681635    0.69767273           0.10000000
## 253 normal            NA    0.36287929                   NA
## 254 normal  0.3170918174    0.52617318           0.05395240
## 255 normal  0.0133920542    0.26356706           0.03301932
## 256 normal  0.0308752192    0.78600259           0.10000000
## 257 normal            NA    0.33462099                   NA
## 258 normal            NA    0.49078368                   NA
## 259      V  0.0184919842    0.34313222           0.10000000
## 260 normal  0.0021512042    0.91285794           0.05693978
## 261 normal  0.1656845359    0.68019223           0.10000000
## 262    dec  0.0058405106    0.71545581           0.06712673
## 263    dec  0.0788938231    0.87700921           0.04726550
## 264 normal  0.4653781716    0.71868138           0.10000000
## 265    inc  0.1669354401    0.92808461           0.04429551
## 266 normal  0.2110097308    0.62603556           0.10000000
## 267 normal  0.0058447973    0.01729436           0.10000000
## 268 normal  0.2483305374    0.40914658           0.10000000
## 269 normal  0.1136872384    0.97732437           0.10000000
## 270 normal            NA    0.67086515                   NA
## 271 normal            NA    0.95930132                   NA
## 272 normal            NA    0.99000000                   NA
## 273 normal            NA    0.72160085                   NA
## 274 normal            NA    0.45775801                   NA
## 275 normal            NA    0.40717480                   NA
## 276 normal            NA    0.39615565                   NA
## 277 normal            NA    0.32313774                   NA
## 278 normal            NA    0.29950033                   NA
## 279 normal            NA    0.52764806                   NA
## 280 normal            NA    0.99000000                   NA
## 281 normal            NA    0.66393783                   NA
## 282 normal  0.9051198371    0.97895636           0.10000000
## 283 normal            NA    0.09763187                   NA
## 284 normal            NA    0.35856393                   NA
## 285 normal  0.2935607070    0.46673352           0.06640747
## 286 normal  0.5423549454    0.60632811           0.10000000
## 287 normal  0.5127468775    0.98124008           0.10000000
## 288 normal  0.1762034063    0.78658682           0.10000000
## 289 normal            NA    0.90262696                   NA
## 290    dec  0.0032645769    0.24876073           0.06900032
## 291 normal  0.1152156311    0.98424471           0.10000000
## 292 normal            NA    0.80404305                   NA
## 293 normal            NA    0.47241650                   NA
## 294 normal            NA    0.60321972                   NA
## 295 normal            NA    0.26368698                   NA
## 296 normal            NA    0.60245407                   NA
## 297 normal            NA    0.25512815                   NA
## 298      V  0.3987479528    0.48863933           0.10000000
## 299    dec  0.0028485368    0.06816432           0.10000000
## 300 normal            NA    0.83362821                   NA
## 301    dec            NA    0.89454450                   NA
## 302    dec            NA    0.92717131                   NA
## 303 normal  0.4633765330    0.65854912           0.10000000
## 304    inc  0.0029445971    0.31766155           0.05672871
## 305 normal  0.0110544482    0.55762181           0.10000000
## 306 normal  0.0465462550    0.54267803           0.10000000
## 307    inc  0.0001331704    0.65210639           0.04431071
## 308 normal  0.0122250283    0.58466467           0.10000000
## 309    dec            NA    0.20339279                   NA
## 310 normal  0.1281287660    0.91238937           0.10000000
## 311    inc  0.2528580683    0.61306726           0.07783656
## 312 normal  0.1131867827    0.42121809           0.10000000
## 313    dec  0.0020364027    0.86303829           0.06618312
## 314    dec            NA    0.85091860                   NA
## 315 normal  0.3821457993    0.24464111           0.10000000
## 316      V  0.0106941968    0.38687640           0.10000000
## 317 normal  0.0405706537    0.43115988           0.10000000
## 318 normal  0.1898669681    0.19691877           0.10000000
## 319    inc  0.0038378254    0.13270965           0.06404683
## 320    inc  0.0008364840    0.61056800           0.03191495
## 321      V  0.0258828756    0.99000000           0.10000000
## 322    dec  0.0010189771    0.94985132           0.04938666
## 323    dec  0.0644445022    0.27656632           0.04420811
## 324    dec            NA    0.91665060                   NA
## 325    dec            NA    0.77902822                   NA
## 326 normal  0.0089630181    0.12152284           0.10000000
## 327    inc  0.0148916291    0.44746774           0.03775720
## 328 normal  0.0010214801    0.01000000           0.10000000
## 329    inc  0.0007138145    0.55432518           0.05754693
## 330    inc  0.0021314025    0.80133779           0.09512893
## 331    dec            NA    0.09480437                   NA
## 332    inV  0.1122760684    0.90484818           0.10000000
## 333    inV  0.0044438277    0.77085518           0.08227287
## 334 normal  0.8853709237    0.57187659           0.07229974
## 335 normal            NA    0.73269072                   NA
## 336 normal  0.9404469717    0.79041719           0.04713877
## 337 normal  0.0459960641    0.46364468           0.10000000
## 338 normal            NA    0.79311814                   NA
## 339 normal            NA    0.92005213                   NA
## 340 normal            NA    0.79683768                   NA
## 341 normal  0.3435123087    0.36203976           0.10000000
## 342 normal            NA    0.91576536                   NA
## 343 normal  0.8911853338    0.20985160           0.10000000
## 344 normal            NA    0.99000000                   NA
## 345 normal            NA    0.99000000                   NA
## 346 normal            NA    0.66277120                   NA
## 347 normal  0.0014045034    0.90193006           0.05967102
## 348 normal  0.0939957052    0.93021561           0.10000000
## 349    dec  0.0024509613    0.91563719           0.06348821
## 350    dec            NA    0.06997305                   NA
## 351 normal  0.2885971056    0.59316100           0.10000000
## 352    inc  0.0075964666    0.98201887           0.07893130
## 353 normal  0.4416304702    0.53370633           0.10000000
## 354 normal  0.0240376698    0.37986138           0.10000000
## 355 normal  0.1579472115    0.54267003           0.09255940
## 356    inc  0.0006646399    0.22310699           0.02179212
## 357 normal            NA    0.96468357                   NA
## 358    inc            NA    0.98813973                   NA
## 359    inc  0.0005901106    0.47456474           0.06051211
## 360    inc  0.0028746542    0.09277991           0.02340907
## 361    inc  0.0003555412    0.94052660           0.03338380
## 362    inc  0.0015927371    0.88514666           0.04739188
## 363 normal            NA    0.52383811                   NA
## 364 normal            NA    0.83443579                   NA
## 365 normal            NA    0.56062486                   NA
## 366 normal            NA    0.98552319                   NA
## 367 normal            NA    0.81613446                   NA
## 368 normal            NA    0.62380354                   NA
## 369 normal  0.6872706626    0.77098020           0.10000000
## 370 normal  0.3934034011    0.97037785           0.10000000
## 371 normal  0.5987529047    0.39295217           0.10000000
## 372 normal  0.4692471717    0.61993205           0.10000000
## 373 normal  0.9509374930    0.08371527           0.10000000
## 374 normal  0.0862647649    0.96285696           0.10000000
## 375 normal            NA    0.50718306                   NA
## 376 normal            NA    0.22261415                   NA
## 377 normal            NA    0.99000000                   NA
## 378 normal            NA    0.12989530                   NA
## 379 normal            NA    0.47241650                   NA
## 380 normal            NA    0.58558896                   NA
## 381 normal            NA    0.45867962                   NA
## 382 normal            NA    0.25512815                   NA
## 383 normal            NA    0.81613446                   NA
## 384 normal            NA    0.99000000                   NA
## 385    dec  0.0035087546    0.99000000           0.10000000
## 386    dec  0.2434152051    0.49839726           0.08854714
## 387 normal            NA    0.48188319                   NA
## 388 normal            NA    0.39329974                   NA
## 389 normal            NA    0.42538070                   NA
## 390    inV  0.0014464121    0.01000000           0.10000000
## 391    dec  0.0015866940    0.31815474           0.10000000
## 392 normal            NA    0.52764806                   NA
## 393 normal            NA    0.49078368                   NA
## 394 normal            NA    0.55877630                   NA
## 395 normal  0.0343193834    0.04390497           0.10000000
## 396 normal  0.0011726002    0.74103864           0.08531773
## 397 normal            NA    0.75089912                   NA
## 398 normal            NA    0.91003089                   NA
## 399 normal  0.0726152904    0.70491802           0.03518273
## 400 normal  0.3356819381    0.97895757           0.07318281
## 401    dec  0.0107179255    0.04548491           0.05405563
## 402 normal  0.9117820276    0.42454488           0.10000000
## 403 normal  0.9420982505    0.64584874           0.09973223
## 404 normal            NA    0.60245407                   NA
## 405 normal            NA    0.57370402                   NA
## 406 normal            NA    0.43908569                   NA
## 407 normal            NA    0.71806755                   NA
## 408 normal            NA    0.35406854                   NA
## 409 normal  0.0930053167    0.50440824           0.10000000
## 410 normal  0.6856188132    0.52351564           0.10000000
## 411 normal            NA    0.24504474                   NA
## 412 normal  0.2749239631    0.25684642           0.06246187
## 413 normal  0.7867601086    0.89175410           0.08092063
## 414 normal            NA    0.12785088                   NA
## 415 normal            NA    0.43908569                   NA
## 416 normal            NA    0.59072441                   NA
## 417 normal            NA    0.09350611                   NA
## 418      V            NA    0.93226882                   NA
## 419 normal            NA    0.65367222                   NA
## 420 normal  0.3685845397    0.68849991           0.08534748
## 421 normal            NA    0.87422193                   NA
## 422 normal            NA    0.08563826                   NA
## 423 normal  0.5596368296    0.81273977           0.10000000
## 424 normal  0.3758163679    0.44154724           0.03359540
## 425 normal  0.0042672823    0.01196259           0.10000000
## 426 normal  0.3078562325    0.58799642           0.10000000
## 427 normal  0.5895865831    0.36486506           0.10000000
## 428 normal            NA    0.57723275                   NA
## 429      V  0.0160076883    0.54562400           0.10000000
## 430    dec  0.0006958561    0.94375794           0.06011460
## 431 normal  0.0165530928    0.33885556           0.05181185
## 432    dec  0.0014690307    0.82677147           0.06155167
## 433    dec  0.0147049001    0.98331815           0.04863769
## 434 normal  0.2019306600    0.35505293           0.06852030
## 435    inc  0.0547339782    0.51292766           0.08145168
## 436 normal  0.7766001795    0.52120194           0.10000000
## 437 normal  0.0034969891    0.01709866           0.10000000
## 438 normal  0.0534691675    0.40835797           0.10000000
## 439    dec  0.0869457087    0.74779872           0.10000000
## 440    dec            NA    0.39124215                   NA
## 441    dec  0.0669530115    0.76651755           0.10000000
## 442    inV  0.0091691202    0.96339992           0.08677164
## 443    inc  0.3191703253    0.06692551           0.10000000
## 444    dec  0.0214787104    0.48897342           0.07440780
## 445    dec  0.0013652185    0.69819620           0.05956728
## 446 normal  0.1454648359    0.67509501           0.10000000
## 447 normal  0.4780583101    0.10584256           0.10000000
## 448 normal  0.2480141545    0.11090320           0.10000000
## 449 normal            NA    0.91046173                   NA
## 450 normal            NA    0.95098535                   NA
## 451 normal  0.1044491381    0.66655518           0.10000000
## 452 normal            NA    0.30681259                   NA
## 453 normal  0.2180904736    0.44006404           0.10000000
## 454 normal  0.8082323901    0.40663907           0.06684781
## 455 normal  0.0412548166    0.63430336           0.10000000
## 456 normal  0.5211585917    0.39782850           0.08246672
## 457    dec            NA    0.71933259                   NA
## 458    dec            NA    0.72446319                   NA
## 459 normal  0.6040573216    0.53410129           0.10000000
## 460 normal  0.2916681490    0.39217694           0.10000000
## 461 normal  0.7223028317    0.41948976           0.10000000
## 462 normal  0.1084626224    0.65577800           0.10000000
## 463 normal  0.2582764329    0.29409409           0.10000000
## 464 normal  0.4398509690    0.42888566           0.08003381
## 465 normal  0.2840458845    0.59486695           0.10000000
## 466 normal  0.0014277061    0.90047218           0.05633525
## 467 normal  0.2865395422    0.23705768           0.10000000
## 468 normal            NA    0.18188810                   NA
## 469 normal            NA    0.15976570                   NA
## 470 normal  0.5337469309    0.86810259           0.10000000
## 471 normal  0.0956556712    0.68889273           0.10000000
## 472 normal  0.8313225880    0.34291945           0.10000000
## 473 normal  0.3633999677    0.09199105           0.10000000
## 474    dec  0.0485433537    0.98698322           0.10000000
## 475    dec  0.2618628941    0.65987241           0.10000000
## 476 normal            NA    0.60431289                   NA
## 477 normal            NA    0.24709685                   NA
## 478 normal            NA    0.11391454                   NA
## 479 normal  0.0018725530    0.61688924           0.10000000
## 480 normal  0.0045330028    0.40720457           0.10000000
## 481 normal  0.8434033328    0.10409103           0.10000000
## 482    inc  0.0018871420    0.98708193           0.08005112
## 483    inc  0.0321265468    0.44594106           0.07012492
## 484 normal  0.7953149276    0.05723842           0.10000000
## 485 normal  0.8763425337    0.24789834           0.10000000
## 486    inc  0.1991193999    0.51272661           0.04682302
## 487    inc  0.0979416325    0.09793150           0.03978390
## 488    inc  0.0019954975    0.01000000           0.04564096
## 489    inc  0.0057273991    0.94597146           0.04042116
## 490 normal            NA    0.52764806                   NA
## 491 normal            NA    0.34824249                   NA
## 492 normal            NA    0.55877630                   NA
## 493 normal            NA    0.52489122                   NA
## 494 normal  0.7291055270    0.48507318           0.10000000
## 495 normal  0.8458150994    0.63313997           0.10000000
## 496 normal  0.5653491315    0.08518531           0.10000000
## 497 normal  0.5084441457    0.63773814           0.10000000
## 498 normal  0.8984149333    0.66111585           0.10000000
## 499 normal  0.5573339667    0.52214549           0.10000000
## 500 normal  0.5545493929    0.62816207           0.10000000
## 501 normal            NA    0.60665927                   NA
## 502    dec  0.0021654573    0.97565386           0.04878858
## 503 normal  0.4665900702    0.78994891           0.10000000
## 504    dec  0.0959235190    0.49065772           0.10000000
## 505 normal  0.3430135685    0.66035391           0.09789181
## 506 normal  0.3119923329    0.33270418           0.10000000
## 507    dec  0.0200256037    0.54512278           0.03619353
## 508 normal            NA    0.23900673                   NA
## 509 normal            NA    0.71270207                   NA
## 510 normal            NA    0.77159575                   NA
## 511    inV  0.0005810677    0.19018531           0.10000000
## 512    inc  0.0007552744    0.94282280           0.07597961
## 513    inV  0.0914082395    0.06927284           0.10000000
## 514 normal            NA    0.95519693                   NA
## 515 normal            NA    0.45482321                   NA
## 516 normal  0.1006915932    0.49042799           0.10000000
## 517    inc  0.0024643175    0.86774851           0.06921034
## 518 normal  0.0844906084    0.56017721           0.10000000
## 519 normal  0.0016835336    0.84051979           0.10000000
## 520    inc  0.0005332240    0.79728079           0.06289418
## 521 normal  0.1265574942    0.72409234           0.10000000
## 522 normal  0.8863827808    0.50837078           0.10000000
## 523 normal            NA    0.80176243                   NA
## 524 normal            NA    0.22527170                   NA
## 525 normal            NA    0.53807791                   NA
## 526 normal            NA    0.09749799                   NA
## 527 normal            NA    0.35390916                   NA
## 528 normal  0.0014801490    0.01906714           0.10000000
## 529 normal            NA    0.86993634                   NA
## 530 normal  0.0816654638    0.56378312           0.07791047
## 531 normal            NA    0.12953070                   NA
## 532    inc  0.0006806673    0.61322715           0.02334206
## 533    inc  0.0274937535    0.03369109           0.05197567
## 534 normal  0.0280993081    0.11992455           0.10000000
## 535    inc  0.0758696847    0.03606255           0.03257357
## 536 normal  0.6624045145    0.66319298           0.10000000
## 537 normal            NA    0.22261415                   NA
## 538 normal            NA    0.26368698                   NA
## 539 normal            NA    0.50873429                   NA
## 540 normal            NA    0.31477691                   NA
## 541 normal  0.6443468825    0.70754920           0.10000000
## 542 normal  0.4947495725    0.76396057           0.10000000
## 543 normal  0.0261314420    0.73656604           0.10000000
## 544 normal  0.2793990244    0.58751225           0.10000000
## 545 normal  0.7690627356    0.96519093           0.10000000
## 546 normal  0.5044169481    0.77207076           0.10000000
## 547 normal  0.7785673725    0.25717192           0.10000000
## 548 normal  0.8168667155    0.87150510           0.10000000
## 549 normal  0.3790405702    0.90666167           0.10000000
## 550 normal            NA    0.52764806                   NA
## 551 normal            NA    0.57370402                   NA
## 552 normal  0.0132829024    0.19082062           0.10000000
## 553 normal  0.0099790896    0.73507316           0.10000000
## 554    dec  0.6498383342    0.70027172           0.10000000
## 555    dec  0.0042721579    0.90992523           0.06002974
## 556    dec            NA    0.80552696                   NA
## 557 normal  0.1848954735    0.64366214           0.10000000
## 558    inc  0.3112148974    0.53555825           0.04172003
## 559 normal  0.6241786059    0.15242914           0.10000000
## 560 normal  0.0172005573    0.06952637           0.10000000
## 561 normal  0.7847467788    0.58902054           0.10000000
## 562    dec  0.1992449718    0.28611701           0.03050834
## 563 normal            NA    0.99000000                   NA
## 564 normal  0.0046388893    0.97253357           0.05161509
## 565 normal  0.5385746245    0.51657337           0.10000000
## 566 normal            NA    0.39936002                   NA
## 567 normal            NA    0.07250199                   NA
## 568 normal  0.3203992020    0.38224140           0.10000000
## 569 normal  0.3059935680    0.10636345           0.10000000
## 570 normal  0.2489490589    0.68053002           0.10000000
## 571 normal  0.7276680259    0.05260056           0.10000000
## 572 normal  0.6812251198    0.06077938           0.10000000
## 573 normal  0.9351464167    0.90023023           0.10000000
## 574    dec  0.0053183296    0.61322294           0.10000000
## 575 normal            NA    0.41899107                   NA
## 576 normal            NA    0.71955548                   NA
## 577    dec  0.0076109911    0.23450826           0.04099869
## 578    dec  0.0328134000    0.22246960           0.02064654
## 579    dec  0.0002883666    0.97780190           0.02221730
## 580    dec  0.0018449969    0.80839866           0.02183680
## 581 normal            NA    0.32777259                   NA
## 582    inc  0.0028575954    0.45112271           0.04904614
## 583    inc  0.0398113834    0.90304802           0.10000000
## 584 normal  0.0358984623    0.01139292           0.05841389
## 585 normal            NA    0.38116763                   NA
## 586    dec  0.0457815688    0.92022892           0.10000000
## 587 normal            NA    0.99000000                   NA
## 588 normal            NA    0.70848314                   NA
## 589 normal            NA    0.73482225                   NA
## 590 normal            NA    0.73309993                   NA
## 591    dec  0.0018361610    0.71453262           0.04579886
## 592    dec  0.0004405844    0.99000000           0.03773727
## 593    dec  0.0052574087    0.52086853           0.02351619
## 594    dec  0.0003913773    0.96088603           0.02433817
## 595    dec  0.0033904124    0.53685740           0.02875047
## 596    dec  0.0002770315    0.77298726           0.03274387
## 597    dec  0.0027428222    0.39839012           0.07449671
## 598 normal            NA    0.57000977                   NA
## 599 normal            NA    0.56062486                   NA
## 600 normal  0.0017833874    0.98660107           0.10000000
## 601 normal  0.9871676093    0.95996392           0.10000000
## 602    dec  0.0039318647    0.99000000           0.07599514
## 603    dec  0.0403474373    0.92066462           0.06890420
## 604 normal            NA    0.71918947                   NA
## 605    inc  0.0056107519    0.03041416           0.03896395
## 606    inc  0.2840515321    0.11706973           0.03303093
## 607 normal            NA    0.17866417                   NA
## 608    inc  0.0002385019    0.64782459           0.03670951
## 609    inc  0.0007294969    0.88306131           0.09560347
## 610    inc  0.0061488392    0.42328319           0.02712567
## 611    inc  0.0125088598    0.85347202           0.04416138
## 612 normal            NA    0.96464305                   NA
## 613    dec            NA    0.89893503                   NA
## 614 normal  0.0854281713    0.98641326           0.10000000
## 615    dec            NA    0.93827184                   NA
## 616    dec            NA    0.57893439                   NA
## 617 normal  0.1738945862    0.93101533           0.10000000
## 618      V  0.0968546712    0.88847974           0.10000000
## 619 normal  0.0401276524    0.01000000           0.10000000
## 620 normal  0.0006859132    0.01906505           0.10000000
## 621 normal  0.0724624130    0.69061266           0.10000000
## 622 normal  0.9459388935    0.82299998           0.10000000
## 623 normal            NA    0.85791459                   NA
## 624 normal  0.5348752522    0.93145650           0.10000000
## 625 normal            NA    0.87982130                   NA
## 626    dec  0.0387375020    0.87623040           0.06786037
## 627 normal            NA    0.26520849                   NA
## 628 normal  0.2432665842    0.52763729           0.10000000
## 629 normal  0.0901444004    0.72287410           0.10000000
## 630 normal  0.0338173841    0.59770660           0.10000000
## 631 normal  0.1345114051    0.73970666           0.10000000
## 632    inc            NA    0.99000000                   NA
## 633 normal            NA    0.92755563                   NA
## 634    dec  0.0018488863    0.73650113           0.06475873
## 635 normal            NA    0.66676281                   NA
## 636 normal  0.8662581894    0.49913520           0.08223692
## 637 normal  0.0888301924    0.77270309           0.10000000
## 638    dec  0.0134816931    0.78088659           0.05613797
## 639 normal  0.3728109056    0.20674737           0.07643660
## 640 normal  0.7221874266    0.09341033           0.10000000
## 641    dec            NA    0.26821580                   NA
## 642 normal            NA    0.08664886                   NA
## 643 normal            NA    0.09703558                   NA
## 644 normal            NA    0.38190476                   NA
## 645 normal            NA    0.46145184                   NA
## 646    dec  0.0860139622    0.57190025           0.02429218
## 647    dec  0.9469765242    0.26893662           0.02475036
## 648 normal  0.6735662837    0.37292321           0.10000000
## 649 normal  0.5530819767    0.76736772           0.10000000
## 650 normal  0.5040539554    0.41667096           0.07580032
## 651 normal            NA    0.33449852                   NA
## 652 normal            NA    0.99000000                   NA
## 653 normal            NA    0.43285556                   NA
## 654 normal            NA    0.08659876                   NA
## 655 normal            NA    0.23658541                   NA
```

```r
# mutate(final_class = case_when()) -> tmp_ready
# 
# tmp_ready
```



```r
# based on manual categories

tmp_ready %>% 
  # filter(category %in% c("normal"),
  # filter(category %in% c("normal"),
  # Reactor_Treatment_Dose %in% ("TR4_CCUG59168")) %>% 
  select(cat_period, Abundance, period, Treatment,Model2,  OTU_Reactor_Treatment_Dose, Model2_reactor_treatment_dose, Reactor_Treatment_Dose, OTU, Day_of_Treatment) %>% 
  group_by(Model2, Model2_reactor_treatment_dose,Treatment, period, OTU_Reactor_Treatment_Dose, cat_period) %>% 
  arrange(Model2, Model2_reactor_treatment_dose, Reactor_Treatment_Dose) %>% 
  summarise(mean_period = mean(Abundance)) %>% 
  group_by(Model2, Model2_reactor_treatment_dose, Treatment, period, cat_period) %>% 
  summarise(sum_cat = sum(mean_period)) %>% 
  ggplot(data = ., aes_string(x = "period", y = "sum_cat",  fill = "cat_period")) + #, color = "grey90") +
  geom_bar(position="fill", stat="identity") + 
  # geom_area(aes_string(x = "Day_of_Treatment", y = "Abundance")) + 
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

![](16S-taxa-classification_Figs_BP_finetuning_files/figure-html/unnamed-chunk-20-1.png)<!-- -->



```r
tmp_ready %>% 
  # filter(category %in% c("normal"),
  # filter(category %in% c("normal"),
  # Reactor_Treatment_Dose %in% ("TR4_CCUG59168")) %>% 
  # filter(category %in% c("normal"),
  # filter(category %in% c("normal"),
  # Reactor_Treatment_Dose %in% ("TR4_CCUG59168")) %>% 
  select(cat_period, Abundance, period, Treatment,Model2,  OTU_Reactor_Treatment_Dose, Model2_reactor_treatment_dose, Reactor_Treatment_Dose, OTU, Day_of_Treatment) %>% 
  ggplot(data = ., aes_string(x = "period", y = "Abundance",  fill = "cat_period")) + #, color = "grey90") +
  geom_bar(position="fill", stat="identity") + 
  # geom_area(aes_string(x = "Day_of_Treatment", y = "Abundance")) + 
  facet_grid(Model2 ~ Treatment, scales = "free", drop = TRUE, space = "free") + 
  theme(legend.position = "bottom") + ggpubr::rotate_x_text(45) + xlab(NULL) + ylab(NULL) +
  scale_fill_manual(values = ggpubr::get_palette(palette = "jco", k = 5)) -> p_cat

p_cat
```

![](16S-taxa-classification_Figs_BP_finetuning_files/figure-html/unnamed-chunk-21-1.png)<!-- -->

```r
p_cat %>%
  export::graph2ppt(append = TRUE, 
                    width = 317.48031496 * 1.5, 
                    height = 0.618 * 317.48031496 * 1.25, paper = "A4",  scaling = 2,
                    file = out_pptx)
```

```
## Exported graph as ~/Desktop/16S-chick.pptx
```

```r
tmp_ready %>% 
  # filter(category %in% c("normal"),
  # filter(category %in% c("normal"),
  # Reactor_Treatment_Dose %in% ("TR4_CCUG59168")) %>% 
  select(cat_period, Abundance, period, Treatment,Model2,  OTU_Reactor_Treatment_Dose, Model2_reactor_treatment_dose, Reactor_Treatment_Dose, OTU, Day_of_Treatment) %>% 
  group_by(Model2_reactor_treatment_dose, period, OTU_Reactor_Treatment_Dose, cat_period) %>% 
  arrange(Model2_reactor_treatment_dose, Reactor_Treatment_Dose) %>% 
  summarise(mean_period = mean(Abundance)) %>% 
  group_by(Model2_reactor_treatment_dose, period, cat_period) %>% 
  summarise(sum_cat = sum(mean_period)) %>% 
  write_tsv(file = "~/Desktop/16S_fig_cat.tsv")
```

```
## `summarise()` has grouped output by 'Model2_reactor_treatment_dose', 'period',
## 'OTU_Reactor_Treatment_Dose'. You can override using the `.groups` argument.
```

```
## `summarise()` has grouped output by 'Model2_reactor_treatment_dose', 'period'.
## You can override using the `.groups` argument.
```


```r
# based on cat_lm categories
tmp_ready %>% 
  # filter(category %in% c("normal"),
  # filter(category %in% c("normal"),
  # Reactor_Treatment_Dose %in% ("TR4_CCUG59168")) %>% 
  arrange(Model2_reactor_treatment_dose) %>% 
  ggplot(data = ., aes_string(x = "period", y = "Abundance", fill = "cat_lm")) +
  geom_bar(position="fill", stat="identity") + 
  # geom_area(aes_string(x = "Day_of_Treatment", y = "Abundance")) + 
  facet_grid(Model2 ~ Treatment, scales = "free", drop = TRUE, space = "free") + 
  theme(legend.position = "bottom") 
```

![](16S-taxa-classification_Figs_BP_finetuning_files/figure-html/unnamed-chunk-22-1.png)<!-- -->


```r
ps_filt2 %>% 
  # filter(Reactor_Treatment_Dose %in% c("TR6_VAN90","TR4_CCUG59168", "TR5_HV292.1", "TR1_CTX+HV292.120")) %>% 
  select(OTU, Abundance, Day_of_Treatment, Model2 ,Model2_reactor_treatment_dose, Treatment) %>% 
  left_join(classif,
            by = c("Model2_reactor_treatment_dose", "OTU")) %>% 
  left_join(day_period) %>%   # filter(category %in% c("normal"),
  # filter(category %in% c("normal"),
  # Reactor_Treatment_Dose %in% ("TR4_CCUG59168")) %>% 
  arrange(Model2_reactor_treatment_dose) %>% 
  ggplot(data = ., aes_string(x = "interaction(period, category)", y = "Abundance", color = "OTU", fill = "OTU")) +
  geom_bar(position="fill", stat="identity") + 
  # geom_area(aes_string(x = "Day_of_Treatment", y = "Abundance")) + 
  facet_grid(Model2 ~ Treatment , scales = "free", drop = TRUE, space = "free") + 
  theme(legend.position = "none") + ggpubr::rotate_x_text(90) 
```

```
## Joining, by = "Day_of_Treatment"
```

```
## Warning: Removed 236 rows containing missing values (`geom_bar()`).
```

![](16S-taxa-classification_Figs_BP_finetuning_files/figure-html/unnamed-chunk-23-1.png)<!-- -->


```r
ps_filt2 %>% 
  mutate(OTU_Reactor_Treatment_Dose = paste0(OTU,"_",Model2_reactor_treatment_dose)) %>% 
  group_by(OTU, Model2_reactor_treatment_dose) %>% 
  arrange(Day_of_Treatment) %>% 
  mutate(global_mean = mean(Abundance),
         global_CV = (sd(Abundance) / mean(Abundance)) * 100) %>% 
  mutate(Mean_avg_change = mean(cummean(Abundance))) %>% 
  # mutate(cum_avg_change_pos = sum(Mean_avg_change[Mean_avg_change>global_mean]), 
  #  cum_avg_change_neg = sum(Mean_avg_change[Mean_avg_change<global_mean])) %>% 
  mutate(cum_avg_change_pos = sum(Abundance[Abundance>mean(nth(Abundance, n = 3))]),
         cum_avg_change_neg = sum(Abundance[Abundance<mean(nth(Abundance,  n = 3))])) -> cum_stuff
# mutate(avg_cum_change_pos = mean(Abundance[Abundance>first(Abundance)]), 
# avg_cum_change_neg = mean(Abundance[Abundance<first(Abundance)])) %>% 

# mutate(Mean_avg_inc = ifelse(Abundance < first(Abundance)), cummean(Abundance), 0) %>% 
# arrange(-cum_avg_change_neg) %>%
```



```r
p_cat$data %>% 
  select(-starts_with("pval"), -Day_of_Treatment, -Abundance) %>% 
  arrange(Model2_reactor_treatment_dose, OTU) %>% 
  # left_join(cum_stuff,
  #           by = c("OTU_Reactor_Treatment_Dose" = "OTU_Reactor_Treatment_Dose"),
  #           copy = FALSE, keep = FALSE,
  #           suffix = c("", ".y")) %>%
  # select(-ends_with(".y")) %>%
  # select(-Abundance, -Day_of_Treatment, -pos:-p1_neg, -OTU_Reactor_Treatment_Dose.y, -cat_lm, -period) %>% 
  # select(-pos:-p1_neg, -OTU_Reactor_Treatment_Dose.y, -cat_lm) %>% 
  distinct(OTU_Reactor_Treatment_Dose, .keep_all = TRUE) %>%
  # as_tibble() %>% 
  # pivot_wider(names_from = cat_period, values_from = OTU_Reactor_Treatment_Dose )  %>% 
  write_tsv(file = "~/Desktop/16S-chick_classif.tsv")
# pivot_wider(names_from = OTU_Reactor_Treatment_Dose, values_from = cum_abund) -> ps_merged
# xlsx::write.xlsx(out_xlsx,
# sheetName = "cat_period_data")


p_cat$data %>% 
  # select(-starts_with("pval")) %>% 
  arrange(Model2_reactor_treatment_dose, OTU, Day_of_Treatment) %>% 
  # left_join(cum_stuff,
  #           by = c("OTU_Reactor_Treatment_Dose" = "OTU_Reactor_Treatment_Dose"),
  #           copy = FALSE, keep = FALSE,
  #           suffix = c("", ".y")) %>% 
  # select(-ends_with(".y")) %>% 
  # select(-Abundance, -Day_of_Treatment, -pos:-p1_neg, -OTU_Reactor_Treatment_Dose.y, -cat_lm, -period) %>% 
  # select(-pos:-p1_neg, -OTU_Reactor_Treatment_Dose.y, -cat_lm) %>% 
  # distinct(OTU_Reactor_Treatment_Dose, .keep_all = TRUE) %>% 
  # as_tibble() %>% 
  # pivot_wider(names_from = cat_period, values_from = OTU_Reactor_Treatment_Dose )  %>% 
  write_tsv(file = "~/Desktop/16S-chick_classif_2.tsv")
# pivot_wider(names_from = OTU_Reactor_Treatment_Dose, values_from = cum_abund) -> ps_merged
# xlsx::write.xlsx(out_xlsx,
# sheetName = "cat_period_data")


cum_stuff %>% 
  select(-Abundance, -Day_of_Treatment) %>% 
  distinct(OTU_Reactor_Treatment_Dose, .keep_all = TRUE) %>% 
  # xlsx::write.xlsx(out_xlsx,
  # sheetName = "additional_info")
  write_tsv(file = "~/Desktop/16S-chick_additional_info.tsv")
```

-   explore stability metrics: <https://microbiome.github.io/tutorials/Stability.html> to detect taxa impacted. or filter by coefficient variation. or Hannah's classification

Bimodal population distribution across individuals is often associated with instable intermediate abundances within individuals:


```r
# Use relative abundances
pseq <- microbiome::transform(ps_filt, "compositional")

# Merge rare taxa to speed up examples
pseq <- microbiome::aggregate_rare(pseq, level = "Strain", detection = .1/100, prevalence = 10/100)

# For cross-sectional analysis, include
# only the treated points.

pseq0 <- subset_samples(pseq, Day_of_Treatment >= 0)

# intermediate_stability: within subjects:
pseq %>% 
  microViz::ps_mutate(subject = Model2_reactor_treatment_dose,
                      time = Day_of_Treatment) %>% 
  microbiome::intermediate_stability(., #output = "scores", 
                                     output = "scores") -> intermediate.stability

intermediate.stability %>% 
  head()
```


```r
# Bimodality is better estimated from abundances at log scale (such as CLR)
pseq0.clr <- microbiome::transform(pseq0, "clr")

set.seed(4433)

# bimodality: between subjects
bimodality.score <- microbiome::bimodality(pseq0.clr, method = "potential_analysis",
                                           bs.iter = 999, peak.threshold = 5,
                                           min.density = 5)

bimodality.score %>% 
  head()
```


```r
taxa <-microbiome::taxa(pseq0)
df <- data.frame(group = taxa,
                 intermediate.stability = intermediate.stability[taxa],
                 bimodality = bimodality.score[taxa])
theme_set(theme_bw(20))

library(ggrepel)

p <- ggplot(df,
            aes(x = intermediate.stability, y = bimodality, label = '')) + #Doesn't add labels to points
  geom_text_repel() +
  geom_point()

#Labels only those that have bimodality > 0.4 and intermediate stability < 0, adjusts the placement of labels
p <- p + geom_text(aes(label = ifelse(bimodality > 0.75 | intermediate.stability < -0.25, group, '')), vjust = "inward", hjust = "inward")

print(p)
```


```r
p$data
```


```r
# Log10 abundance for a selected taxonomic group
# Pick the most bimodal taxa as an example
tax  <- names(which.max(bimodality.score))

# Detect tipping points at log10 abundances
x <- microbiome::abundances(microbiome::transform(pseq, "clr"))[tax,]

# Bootstrapped potential analysis to identify potential minima
# in practice, use more bootstrap iterations
potential.minima <- microbiome::potential_analysis(x, bs.iter = 10)$minima

# Same with earlywarnings package (without bootstrap ie. less robust)
# library(earlywarnings)
# res <- livpotential_ews(x)$min.points

# Identify the potential minimum locations as tipping point candidates
tipping.samples <- sapply(potential.minima, function (m) {names(which.min(abs(sort(x) - m)))})
tipping.point <- microbiome::abundances(pseq)[tax, tipping.samples]
print(tipping.point)
```



```r
# Illustrate the detected tipping point
plot(density(x), main = tax)
abline(v = x[tipping.samples])
```



```r
# # Bimodality hotplot:
# # Consider a unique sample from each subject: the baseline time point 
# p <- microbiome::hotplot(pseq0, tax, tipping.point = 0.005)
# print(p)
# 
# # Visualize bimodality
# pv <- plot_tipping(pseq, tax, tipping.point = 0.005)
# print(pv)
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
##  [1] gdtools_0.2.4        microViz_0.9.7       speedyseq_0.5.3.9018
##  [4] reshape2_1.4.4       scales_1.2.1         phyloseq_1.42.0     
##  [7] forcats_0.5.2        stringr_1.4.1        dplyr_1.0.10        
## [10] purrr_0.3.5          readr_2.1.3          tidyr_1.2.1         
## [13] tibble_3.1.8         ggplot2_3.4.0        tidyverse_1.3.2     
## 
## loaded via a namespace (and not attached):
##   [1] uuid_1.1-0             readxl_1.4.1           backports_1.4.1       
##   [4] systemfonts_1.0.4      plyr_1.8.8             igraph_1.3.5          
##   [7] splines_4.2.2          listenv_0.8.0          GenomeInfoDb_1.34.3   
##  [10] digest_0.6.30          aTSA_3.1.2             foreach_1.5.2         
##  [13] htmltools_0.5.3        fansi_1.0.3            magrittr_2.0.3        
##  [16] googlesheets4_1.0.1    cluster_2.1.4          openxlsx_4.2.5.1      
##  [19] tzdb_0.3.0             globals_0.16.1         Biostrings_2.66.0     
##  [22] extrafont_0.18         modelr_0.1.9           vroom_1.6.0           
##  [25] officer_0.4.4          extrafontdb_1.0        xts_0.12.2            
##  [28] timechange_0.1.1       tseries_0.10-52        colorspace_2.0-3      
##  [31] corncob_0.3.0          rvest_1.0.3            ggrepel_0.9.2         
##  [34] textshaping_0.3.6      haven_2.5.1            xfun_0.34             
##  [37] crayon_1.5.2           RCurl_1.98-1.9         jsonlite_1.8.3        
##  [40] survival_3.4-0         zoo_1.8-11             iterators_1.0.14      
##  [43] ape_5.6-2              glue_1.6.2             rvg_0.2.5             
##  [46] gtable_0.3.1           gargle_1.2.1           zlibbioc_1.44.0       
##  [49] XVector_0.38.0         car_3.1-1              Rttf2pt1_1.3.11       
##  [52] Rhdf5lib_1.20.0        future.apply_1.10.0    BiocGenerics_0.44.0   
##  [55] quantmod_0.4.20        abind_1.4-5            DBI_1.1.3             
##  [58] rstatix_0.7.1          Rcpp_1.0.9             xtable_1.8-4          
##  [61] bit_4.0.4              stats4_4.2.2           htmlwidgets_1.5.4     
##  [64] httr_1.4.4             ellipsis_0.3.2         pkgconfig_2.0.3       
##  [67] farver_2.1.1           sass_0.4.2             dbplyr_2.2.1          
##  [70] utf8_1.2.2             here_1.0.1             tidyselect_1.2.0      
##  [73] labeling_0.4.2         rlang_1.0.6            munsell_0.5.0         
##  [76] cellranger_1.1.0       tools_4.2.2            cachem_1.0.6          
##  [79] cli_3.4.1              generics_0.1.3         devEMF_4.1-1          
##  [82] ade4_1.7-20            export_0.3.0           broom_1.0.1           
##  [85] evaluate_0.18          biomformat_1.26.0      fastmap_1.1.0         
##  [88] yaml_2.3.6             ragg_1.2.4             bit64_4.0.5           
##  [91] knitr_1.40             fs_1.5.2               zip_2.2.2             
##  [94] rgl_0.110.2            future_1.29.0          nlme_3.1-160          
##  [97] xml2_1.3.3             compiler_4.2.2         rstudioapi_0.14       
## [100] curl_4.3.3             ggsignif_0.6.4         reprex_2.0.2          
## [103] bslib_0.4.1            stringi_1.7.8          highr_0.9             
## [106] stargazer_5.2.3        lattice_0.20-45        Matrix_1.5-1          
## [109] ggsci_2.9              vegan_2.6-4            microbiome_1.19.1     
## [112] permute_0.9-7          urca_1.3-3             multtest_2.54.0       
## [115] vctrs_0.5.0            pillar_1.8.1           lifecycle_1.0.3       
## [118] rhdf5filters_1.10.0    jquerylib_0.1.4        flextable_0.8.3       
## [121] cowplot_1.1.1          data.table_1.14.4      bitops_1.0-7          
## [124] R6_2.5.1               IRanges_2.32.0         parallelly_1.32.1     
## [127] codetools_0.2-18       MASS_7.3-58.1          assertthat_0.2.1      
## [130] rhdf5_2.42.0           rprojroot_2.0.3        withr_2.5.0           
## [133] S4Vectors_0.36.0       GenomeInfoDbData_1.2.9 mgcv_1.8-41           
## [136] parallel_4.2.2         hms_1.1.2              quadprog_1.5-8        
## [139] grid_4.2.2             rmarkdown_2.18         carData_3.0-5         
## [142] googledrive_2.0.0      Rtsne_0.16             TTR_0.24.3            
## [145] ggpubr_0.4.0           base64enc_0.1-3        Biobase_2.56.0        
## [148] lubridate_1.9.0
```
