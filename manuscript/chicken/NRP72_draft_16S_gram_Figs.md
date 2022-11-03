---
title: "NTP72 - 16S - taxa gram - chicken draft "
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
   subset_samples(Day_of_Treatment > -6 & Day_of_Treatment <= 30 | Reactor == "DONOR") %>% 
  subset_samples(Model == "Chicken") %>% 
   subset_samples(Reactor != "CR" & Model2 == "Chicken1" | Model2 == "Chicken2" & Reactor %in% c("TR2", "TR3", "TR4", "TR5", "TR6", "TR7", "CR", "DONOR")) %>% 
   subset_samples(Reactor != "IR") -> ps_filtered
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
## 220  p-3-G1-S265               14       HV292.1     TR5 Chicken2 Chicken     t2
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
## 229  p-3-H5-S281               -5           VAN     TR6 Chicken2 Chicken   <NA>
## 203  p-3-B5-S209               -5     UNTREATED      CR Chicken2 Chicken   <NA>
## 214  p-3-E5-S245               -5           CTX     TR3 Chicken2 Chicken   <NA>
## 1       D-0-S154             <NA>         DONOR   DONOR Chicken1 Chicken   <NA>
## 206  p-3-D2-S230                0   CTX+HV292.1     TR2 Chicken2 Chicken     t1
## 83   TR5-34-S127               18 VAN+CCUG59168     TR5 Chicken1 Chicken     t2
## 198  p-3-A6-S198               -5 VAN+CCUG59168     TR7 Chicken2 Chicken   <NA>
## 199  p-3-B1-S205               14     UNTREATED      CR Chicken2 Chicken     t2
## 222  p-3-G3-S267               -1       HV292.1     TR5 Chicken2 Chicken   <NA>
## 216  p-3-F2-S254                0     CCUG59168     TR4 Chicken2 Chicken     t1
## 176  p-2-G4-S172                6           CTX     TR3 Chicken2 Chicken   <NA>
## 98   TR6-28-S158               12     CCUG59168     TR6 Chicken1 Chicken     t2
## 211  p-3-E2-S242                0           CTX     TR3 Chicken2 Chicken     t1
## 16   TR1-34-S317               18   CTX+HV292.1     TR1 Chicken1 Chicken     t2
## 207  p-3-D3-S231               -1   CTX+HV292.1     TR2 Chicken2 Chicken   <NA>
## 228  p-3-H4-S280               -2           VAN     TR6 Chicken2 Chicken   pret
## 224  p-3-G5-S269               -5       HV292.1     TR5 Chicken2 Chicken   <NA>
## 215  p-3-F1-S253               14     CCUG59168     TR4 Chicken2 Chicken     t2
## 200  p-3-B2-S206                0     UNTREATED      CR Chicken2 Chicken     t1
## 65   TR4-28-S183               12           VAN     TR4 Chicken1 Chicken     t2
## 144  p-2-C5-S125                6 VAN+CCUG59168     TR7 Chicken2 Chicken   <NA>
## 208  p-3-D4-S232               -2   CTX+HV292.1     TR2 Chicken2 Chicken   pret
## 67   TR4-34-S342               18           VAN     TR4 Chicken1 Chicken     t2
## 209  p-3-D5-S233               -5   CTX+HV292.1     TR2 Chicken2 Chicken   <NA>
## 212  p-3-E3-S243               -1           CTX     TR3 Chicken2 Chicken   <NA>
## 82   TR5-31-S313               15 VAN+CCUG59168     TR5 Chicken1 Chicken     t2
## 219  p-3-F5-S257               -5     CCUG59168     TR4 Chicken2 Chicken   <NA>
## 41   TR3-17-S222                1       HV292.1     TR3 Chicken1 Chicken     t1
## 53   TR3-46-S224               30       HV292.1     TR3 Chicken1 Chicken     t4
## 102  TR6-40-S160               24     CCUG59168     TR6 Chicken1 Chicken     t3
## 195  p-3-A3-S195                0 VAN+CCUG59168     TR7 Chicken2 Chicken     t1
## 113  p-1-H10-S94                1   CTX+HV292.1     TR2 Chicken2 Chicken     t1
## 64   TR4-25-S333                9           VAN     TR4 Chicken1 Chicken   <NA>
## 128 p-2-B11-S119               12           VAN     TR6 Chicken2 Chicken     t2
## 101  TR6-37-S312               21     CCUG59168     TR6 Chicken1 Chicken     t3
## 28   TR2-21-S228                5           CTX     TR2 Chicken1 Chicken     t5
## 217  p-3-F3-S255               -1     CCUG59168     TR4 Chicken2 Chicken   <NA>
## 185  p-2-H2-S182                4     CCUG59168     TR4 Chicken2 Chicken     t4
## 205  p-3-D1-S229               14   CTX+HV292.1     TR2 Chicken2 Chicken     t2
## 34   TR2-37-S330               21           CTX     TR2 Chicken1 Chicken     t3
## 63   TR4-22-S344                6           VAN     TR4 Chicken1 Chicken   <NA>
## 25   TR2-18-S326                2           CTX     TR2 Chicken1 Chicken     t2
## 221  p-3-G2-S266                0       HV292.1     TR5 Chicken2 Chicken     t1
## 7    TR1-17-S246                1   CTX+HV292.1     TR1 Chicken1 Chicken     t1
## 223  p-3-G4-S268               -2       HV292.1     TR5 Chicken2 Chicken   pret
## 57   TR4-16-S276                0           VAN     TR4 Chicken1 Chicken     t1
## 110  p-1-E11-S59                1 VAN+CCUG59168     TR7 Chicken2 Chicken     t1
## 111  p-1-F10-S70                1     UNTREATED      CR Chicken2 Chicken     t1
## 39   TR3-15-S338               -1       HV292.1     TR3 Chicken1 Chicken   <NA>
## 182  p-2-H1-S181                3     CCUG59168     TR4 Chicken2 Chicken     t3
## 94   TR6-20-S132                4     CCUG59168     TR6 Chicken1 Chicken     t4
## 29   TR2-22-S364                6           CTX     TR2 Chicken1 Chicken   <NA>
## 54   TR4-13-S142               -3           VAN     TR4 Chicken1 Chicken   pret
## 218  p-3-F4-S256               -2     CCUG59168     TR4 Chicken2 Chicken   pret
## 85   TR5-40-S335               24 VAN+CCUG59168     TR5 Chicken1 Chicken     t3
## 114  p-1-H11-S95                2   CTX+HV292.1     TR2 Chicken2 Chicken     t2
## 106  p-1-B11-S23                1     CCUG59168     TR4 Chicken2 Chicken     t1
## 11   TR1-21-S207                5   CTX+HV292.1     TR1 Chicken1 Chicken     t5
## 91   TR6-17-S138                1     CCUG59168     TR6 Chicken1 Chicken     t1
## 50   TR3-34-S354               18       HV292.1     TR3 Chicken1 Chicken     t2
## 226  p-3-H2-S278                0           VAN     TR6 Chicken2 Chicken     t1
## 86   TR5-46-S310               30 VAN+CCUG59168     TR5 Chicken1 Chicken     t4
## 184 p-2-H11-S191               13     CCUG59168     TR4 Chicken2 Chicken     t2
## 188  p-2-H5-S185                7     CCUG59168     TR4 Chicken2 Chicken   <NA>
## 21   TR2-14-S232               -2           CTX     TR2 Chicken1 Chicken   pret
## 61   TR4-20-S269                4           VAN     TR4 Chicken1 Chicken     t4
## 55   TR4-14-S226               -2           VAN     TR4 Chicken1 Chicken   pret
## 44   TR3-20-S265                4       HV292.1     TR3 Chicken1 Chicken     t4
## 210  p-3-E1-S241               14           CTX     TR3 Chicken2 Chicken     t2
## 131  p-2-B3-S111                4           VAN     TR6 Chicken2 Chicken     t4
## 35   TR2-40-S259               24           CTX     TR2 Chicken1 Chicken     t3
## 149  p-2-D1-S133                3     UNTREATED      CR Chicken2 Chicken     t3
## 187  p-2-H4-S184                6     CCUG59168     TR4 Chicken2 Chicken   <NA>
## 196  p-3-A4-S196               -1 VAN+CCUG59168     TR7 Chicken2 Chicken   <NA>
## 17   TR1-37-S306               21   CTX+HV292.1     TR1 Chicken1 Chicken     t3
## 89   TR6-15-S149               -1     CCUG59168     TR6 Chicken1 Chicken   <NA>
## 204  p-3-C7-S223             <NA>         DONOR   DONOR Chicken2 Chicken   <NA>
## 56   TR4-15-S218               -1           VAN     TR4 Chicken1 Chicken   <NA>
## 20   TR2-13-S231               -3           CTX     TR2 Chicken1 Chicken   pret
## 62   TR4-21-S254                5           VAN     TR4 Chicken1 Chicken     t5
## 84   TR5-37-S368               21 VAN+CCUG59168     TR5 Chicken1 Chicken     t3
## 133  p-2-B5-S113                6           VAN     TR6 Chicken2 Chicken   <NA>
## 108  p-1-C11-S35                1       HV292.1     TR5 Chicken2 Chicken     t1
## 227  p-3-H3-S279               -1           VAN     TR6 Chicken2 Chicken   <NA>
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
## 225  p-3-H1-S277               14           VAN     TR6 Chicken2 Chicken     t2
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
## 202  p-3-B4-S208               -2     UNTREATED      CR Chicken2 Chicken   pret
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
## 201  p-3-B3-S207               -1     UNTREATED      CR Chicken2 Chicken   <NA>
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
## 213  p-3-E4-S244               -2           CTX     TR3 Chicken2 Chicken   pret
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
## 220        5758
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
## 229        9654
## 203        9849
## 214        9899
## 1          9959
## 206       10389
## 83        10416
## 198       11147
## 199       11172
## 222       11293
## 216       11348
## 176       11486
## 98        11501
## 211       11549
## 16        11620
## 207       11856
## 228       12409
## 224       12689
## 215       12815
## 200       13057
## 65        13155
## 144       13160
## 208       13423
## 67        13435
## 209       13737
## 212       13766
## 82        13865
## 219       13927
## 41        13985
## 53        14010
## 102       14053
## 195       14255
## 113       14320
## 64        14385
## 128       14537
## 101       14568
## 28        14695
## 217       14894
## 185       14980
## 205       15414
## 34        15567
## 63        15575
## 25        15579
## 221       15581
## 7         16039
## 223       16084
## 57        16197
## 110       16438
## 111       16732
## 39        16798
## 182       17000
## 94        17578
## 29        17687
## 54        17927
## 218       18325
## 85        18683
## 114       18707
## 106       18848
## 11        18979
## 91        19186
## 50        19218
## 226       19267
## 86        19284
## 184       19403
## 188       19506
## 21        19593
## 61        19601
## 55        19829
## 44        19878
## 210       19937
## 131       20114
## 35        20231
## 149       20281
## 187       20384
## 196       20395
## 17        20499
## 89        20502
## 204       20612
## 56        20762
## 20        20811
## 62        20932
## 84        21108
## 133       21198
## 108       21367
## 227       21578
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
## 225       27516
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
## 202       30408
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
## 201       33197
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
## 213       39627
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
ps_rare_gram_int <- ps_rare

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

### chicken2:

#### gram:


```r
require(microViz)
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

```r
ps_rare_gram_int %>% 
  subset_samples(Model2 == "Chicken2" & Treatment != "DONOR" & Day_of_Treatment >= -1) %>% 
  microViz::ps_mutate(Reactor_Treatment_Dose = paste0(Reactor_Treatment_Dose, "_", Fermentation)) %>% 
  # microViz::ps_mutate(Treatment = fct_relevel(Treatment, sample_data(ps_rare)$Treatment %>%  levels())) %>% 
  transform_sample_counts(function(x) x/sum(x) * 100) %>%
  physeq_sel_tax_table(c("gram_stain")) %>% 
  # subset_taxa(Strain %in% top_taxa_per_group_h2) %>% # extract only the taxa to display - after percentage normalisation
  microViz::phyloseq_validate() %>%
  # microViz::tax_fix() %>%
  microViz::ps_arrange(Day_of_Treatment) %>% 
  # microViz::tax_mutate(gram = left_join(tax_table(.) %>% data.frame(), gram,
  # by = c("Strain" = "ASV") %>% select(gram_neg_class))) %>% 
  microViz::comp_barplot(
    tax_transform_for_plot = "identity",
    tax_level = "gram_stain",
    label = "Day_of_Treatment", # name an alternative variable to label axes
    palette = microViz::distinct_palette(3, pal = "greenArmytage"),
    n_taxa = 10, # give more taxa unique colours
    merge_other = FALSE, # split the "other" category to display alpha diversity
    bar_width = 0.9, # reduce the bar width to 70% of one row
    bar_outline_colour = "grey5",
    sample_order = "default") +
  ylab("Proportion - %") + xlab( " Day (Treatment) ") -> p_gram_h2
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

```
## Warning in ps_counts(data, warn = TRUE): otu_table of counts is NOT available!
## Available otu_table contains non-zero values that are less than 1
```

```
## Warning in ps_counts(data = data): otu_table of counts is NOT available!
## Available otu_table contains non-zero values that are less than 1
```

```r
p_gram_h2 + theme_light() + facet_grid(. ~ Treatment, scales = "free", space = "free_x", drop = TRUE) -> p_gram_h2

p_gram_h2 + theme(legend.position = "none") + ggpubr::rotate_x_text(90) -> p_gram_h2_no_leg

p_gram_h2_no_leg
```

![](NRP72_draft_16S_gram_Figs_files/figure-html/unnamed-chunk-11-1.png)<!-- -->

```r
p_gram_h2 %>%  
  ggpubr::get_legend() %>% 
  ggpubr::as_ggplot() -> p_gram_h2_leg

p_gram_h2_leg
```

![](NRP72_draft_16S_gram_Figs_files/figure-html/unnamed-chunk-11-2.png)<!-- -->

saving outputs:



```r
p_gram_h2_no_leg %>%
    export::graph2ppt(append = TRUE,
                    width = 317.48031496 * 1.25,
                    height = 0.618 * 317.48031496 * 1 , paper = "A4",  scaling = 2,
                    file = out_pptx)
```

```
## Exported graph as ~/Desktop/16S-chick.pptx
```

```r
p_gram_h2_leg %>%
    export::graph2ppt(append = TRUE,
                    width = 317.48031496 * 1,
                    height = 0.618 * 317.48031496 * 1 , paper = "A4",  scaling = 2,
                    file = out_pptx)
```

```
## Exported graph as ~/Desktop/16S-chick.pptx
```

```r
p_gram_h2_no_leg$data %>% 
  select(OTU, Sample, Abundance,Day_of_Treatment, Reactor_Treatment_Dose) %>% 
  pivot_wider(names_from = OTU, values_from = Abundance) %>% 
  xlsx::write.xlsx(file = out_xlsx,
                   append = TRUE,
                   sheetName = "chicken2_gram",
                   showNA = TRUE)
```


```r
require(ggalluvial)
```

```
## Loading required package: ggalluvial
```

```r
p_gram_h2$data %>%
  # filter(Model2 == "Human2") %>%
  mutate(gram_stain = factor(gram_stain, levels = c("Gram-positive", "Gram-negative", "unknown-gram"))) %>%
  # mutate(combo = paste0(Model2, "_" , Reactor_F1_F2)) %>%
  # mutate(combo = as.factor(combo)) %>%
  select(SAMPLE, gram_stain,OTU, Abundance, Day_of_Treatment, Treatment) %>% 
  ggplot(.,
         aes(x = Day_of_Treatment, 
             y = Abundance,
             stratum = gram_stain, 
             alluvium = OTU,
             fill = gram_stain,
             label = gram_stain))+
  facet_grid(Treatment  ~ .)+
  # scale_fill_manual(values = c("red", "yellow", "green3"))+
  geom_flow(stat = "alluvium", lode.guidance = "frontback", color = "darkgray")+
  geom_stratum() +
  theme_light() +
  # geom_text(stat = "alluvium", aes(label = top), lode.guidance = "frontback") +
  # scale_fill_manual(values = col_ASVs) +
  scale_color_manual(values  =   microViz::distinct_palette(3, pal = "greenArmytage")) +
  scale_fill_manual(values  =   microViz::distinct_palette(3, pal = "greenArmytage")) +
  # scale_color_manual(values  = microViz::distinct_palette(10, pal = "kelly")) +
  # scale_fill_manual(values  = microViz::distinct_palette(10, pal = "kelly")) +
  ylab("Proportion - %") -> p_allu_human2_gram

p_allu_human2_gram + theme(legend.position = "none") -> p_allu_human2_gram_noleg 

p_allu_human2_gram_noleg 
```

![](NRP72_draft_16S_gram_Figs_files/figure-html/unnamed-chunk-13-1.png)<!-- -->

```r
p_allu_human2_gram %>% 
  ggpubr::get_legend() %>% 
  ggpubr::as_ggplot() ->  p_allu_human2_gram_leg

p_allu_human2_gram_leg
```

![](NRP72_draft_16S_gram_Figs_files/figure-html/unnamed-chunk-13-2.png)<!-- -->

saving outputs:



```r
p_allu_human2_gram_noleg %>%
    export::graph2ppt(append = TRUE,
                    width = 317.48031496 * 1.25,
                    height = 0.618 * 317.48031496 * 1 , paper = "A4",  scaling = 2,
                    file = out_pptx)
```

```
## Exported graph as ~/Desktop/16S-chick.pptx
```


### chicken1 :

#### gram:


```r
ps_rare_gram_int %>% 
  subset_samples(Model2 == "Chicken1" & Treatment != "DONOR" & Day_of_Treatment >= -1) %>% 
  # microViz::ps_mutate(Reactor_Treatment_Dose = paste0(Reactor_Treatment_Dose, "_", Fermentation)) %>% 
  transform_sample_counts(function(x) x/sum(x) * 100) %>%
  physeq_sel_tax_table(c("gram_stain")) %>% 
  # subset_taxa(Strain %in% top_taxa_per_group_h2) %>% # extract only the taxa to display - after percentage normalisation
  microViz::phyloseq_validate() %>%
  microViz::tax_fix() %>%
  microViz::ps_arrange(Day_of_Treatment) %>% 
  # microViz::tax_mutate(gram = left_join(tax_table(.) %>% data.frame(), gram,
  # by = c("Strain" = "ASV") %>% select(gram_neg_class))) %>% 
  microViz::comp_barplot(
    tax_transform_for_plot = "identity",
    tax_level = "gram_stain",
    label = "Day_of_Treatment", # name an alternative variable to label axes
    palette = microViz::distinct_palette(2, pal = "greenArmytage"),
    n_taxa = 10, # give more taxa unique colours
    merge_other = FALSE, # split the "other" category to display alpha diversity
    bar_width = 0.9, # reduce the bar width to 70% of one row
    bar_outline_colour = "grey5",
    sample_order = "default") +
  ylab("Proportion - %") + xlab( " Day (Treatment) ") -> p_gram_h1
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
p_gram_h1 + facet_grid(. ~ Treatment, scales = "free", space = "free_x", drop = TRUE) -> p_gram_h1

p_gram_h1 + theme_light() + theme(legend.position = "none") + ggpubr::rotate_x_text(90) -> p_gram_h1_no_leg

p_gram_h1_no_leg
```

![](NRP72_draft_16S_gram_Figs_files/figure-html/unnamed-chunk-15-1.png)<!-- -->

```r
p_gram_h1 %>%  
  ggpubr::get_legend() %>% 
  ggpubr::as_ggplot() -> p_gram_h1_leg
```

```r
p_gram_h1_no_leg %>%
    export::graph2ppt(append = TRUE,
                    width = 317.48031496 * 1.25,
                    height = 0.618 * 317.48031496 * 1 , paper = "A4",  scaling = 2,
                    file = out_pptx)
```

```
## Exported graph as ~/Desktop/16S-chick.pptx
```

```r
p_gram_h1_leg %>%
    export::graph2ppt(append = TRUE,
                    width = 317.48031496 * 1,
                    height = 0.618 * 317.48031496 * 1 , paper = "A4",  scaling = 2,
                    file = out_pptx)
```

```
## Exported graph as ~/Desktop/16S-chick.pptx
```

```r
p_gram_h1_no_leg$data %>% 
  select(OTU, Sample, Abundance,Day_of_Treatment, Reactor_Treatment_Dose) %>% 
  pivot_wider(names_from = OTU, values_from = Abundance) %>% 
  xlsx::write.xlsx(file = out_xlsx,
                   append = TRUE,
                   sheetName = "chicken1_gram",
                   showNA = TRUE)
```


```r
p_gram_h1$data %>%
  # filter(Model2 == "Human2") %>%
  mutate(gram_stain = factor(gram_stain, levels = c("Gram-positive", "Gram-negative", "unknown-gram"))) %>%
  # mutate(combo = paste0(Model2, "_" , Reactor_F1_F2)) %>%
  # mutate(combo = as.factor(combo)) %>%
  select(SAMPLE, gram_stain,OTU, Abundance, Day_of_Treatment, Treatment) %>% 
  ggplot(.,
         aes(x = Day_of_Treatment, 
             y = Abundance,
             stratum = gram_stain, 
             alluvium = OTU,
             fill = gram_stain,
             label = gram_stain))+
  facet_grid(Treatment  ~ .)+
  # scale_fill_manual(values = c("red", "yellow", "green3"))+
  geom_flow(stat = "alluvium", lode.guidance = "frontback", color = "darkgray")+
  geom_stratum() +
  theme_light() +
  # geom_text(stat = "alluvium", aes(label = top), lode.guidance = "frontback") +
  # scale_fill_manual(values = col_ASVs) +
  scale_color_manual(values  =   microViz::distinct_palette(3, pal = "greenArmytage")) +
  scale_fill_manual(values  =   microViz::distinct_palette(3, pal = "greenArmytage")) +
  # scale_color_manual(values  = microViz::distinct_palette(10, pal = "kelly")) +
  # scale_fill_manual(values  = microViz::distinct_palette(10, pal = "kelly")) +
  ylab("Proportion - %") -> p_allu_h1_gram 

p_allu_h1_gram + theme(legend.position = "none") -> p_allu_h1_gram_noleg 

p_allu_h1_gram_noleg 
```

![](NRP72_draft_16S_gram_Figs_files/figure-html/unnamed-chunk-17-1.png)<!-- -->

```r
p_allu_h1_gram %>% 
  ggpubr::get_legend() %>% 
  ggpubr::as_ggplot() ->  p_allu_h1_gram_leg

p_allu_h1_gram_leg
```

![](NRP72_draft_16S_gram_Figs_files/figure-html/unnamed-chunk-17-2.png)<!-- -->



```r
p_allu_h1_gram_noleg %>%
    export::graph2ppt(append = TRUE,
                    width = 317.48031496 * 1.25,
                    height = 0.618 * 317.48031496 * 1 , paper = "A4",  scaling = 2,
                    file = out_pptx)
```

```
## Exported graph as ~/Desktop/16S-chick.pptx
```

#### intrinsic:


```r
# ps_rare_gram_int %>% 
#   subset_samples(Model2 == "Chicken1" & Treatment != "DONOR" & Day_of_Treatment >= -1) %>% 
#   microViz::ps_mutate(Reactor_Treatment_Dose = paste0(Reactor_Treatment_Dose, "_", Fermentation)) %>% 
#   transform_sample_counts(function(x) x/sum(x) * 100) %>%
#   physeq_sel_tax_table(c( "intrinsic")) %>% 
#   # microbiome::transform("log") %>% 
#   # speedyseq::tax_glom("intrinsic") %>% 
#   # physeq_sel_tax_table(c("gram_stain","intrinsic")) %>% 
#   # physeq_glom_rename(phyloseq = ., speedyseq = TRUE, rename_ASV = FALSE, taxrank = "gram_stain") %>% 
#   # subset_taxa(Strain %in% top_taxa_per_group_h2) %>% # extract only the taxa to display - after percentage normalisation
#   microViz::phyloseq_validate() %>%
#   microViz::tax_fix() %>%
#   microViz::ps_arrange(Day_of_Treatment) %>% 
#   # microViz::tax_mutate(gram = left_join(tax_table(.) %>% data.frame(), gram,
#   # by = c("Strain" = "ASV") %>% select(gram_neg_class))) %>% 
#   microViz::comp_barplot(
#     tax_transform_for_plot = "identity",
#     tax_level = "intrinsic",
#     label = "Day_of_Treatment", # name an alternative variable to label axes
#     # palette = microViz::distinct_palette(3, pal = "brewerPlus"),
#     n_taxa = 10, # give more taxa unique colours
#     merge_other = FALSE, # split the "other" category to display alpha diversity
#     bar_width = 0.9, # reduce the bar width to 70% of one row
#     bar_outline_colour = "grey5",
#     sample_order = "default") +
#   ylab("Proportion - %") + xlab( " Day (Treatment) ") -> p_int_h1_f1
# 
# p_int_h1_f1 + facet_grid(Kingdom ~ Reactor_Treatment_Dose, scales = "free", space = "free_x", drop = TRUE) -> p_int_h1_f1
# 
# p_int_h1_f1 %>% 
#   ggpubr::set_palette(p = ., palette = "jco") -> p_int_h1_f1
# 
# p_int_h1_f1 + theme_light() +theme(legend.position = "none") -> p_int_h1_f1_no_leg
# 
# p_int_h1_f1_no_leg
# 
# p_int_h1_f1 %>%  
#   ggpubr::get_legend() %>% 
#   ggpubr::as_ggplot() -> p_int_h1_f1_leg
# 
# p_int_h1_f1_leg
```



```r
# p_int_h1_f1$data %>%
#   # filter(Model2 == "Human2") %>%
#   # microViz::ps_mutate(Reactor_Treatment_Dose = paste0(Reactor_Treatment_Dose, "_", Fermentation)) %>% 
#   # mutate(combo = paste0(Model2, "_" , Reactor_F1_F2)) %>%
#   # mutate(combo = as.factor(combo)) %>%
#   select(SAMPLE, top,OTU, Abundance, Day_of_Treatment, Reactor_Treatment_Dose) %>% 
#   ggplot(.,
#          aes(x = Day_of_Treatment, 
#              y = Abundance,
#              stratum = top, 
#              alluvium = OTU,
#              fill = top,
#              label = top))+
#   facet_grid(Reactor_Treatment_Dose  ~ .)+
#   # scale_fill_manual(values = c("red", "yellow", "green3"))+
#   geom_flow(stat = "alluvium", lode.guidance = "frontback", color = "darkgray")+
#   geom_stratum() +
#   theme_light() +
#   # geom_text(stat = "alluvium", aes(label = top), lode.guidance = "frontback") +
#   # scale_fill_manual(values = col_ASVs) +
#   scale_color_manual(values  =  microViz::distinct_palette(3, pal = "greenArmytage")) +
#   # scale_color_manual(values  = microViz::distinct_palette(10, pal = "kelly")) +
#   # scale_fill_manual(values  = microViz::distinct_palette(10, pal = "kelly")) +
#   ylab("Proportion - %") -> p_allu_h1_f1_int
# 
# p_allu_h1_f1_int %>% 
#   ggpubr::set_palette(p = ., palette = "jco") -> p_allu_h1_f1_int
# 
# p_allu_h1_f1_int + theme(legend.position = "none") -> p_allu_h1_f1_int_noleg 
# 
# p_allu_h1_f1_int_noleg 
# 
# p_allu_h1_f1_int %>% 
#   ggpubr::get_legend() %>% 
#   ggpubr::as_ggplot() -> p_allu_h1_f1_int_leg
```


```r
# ggpubr::ggarrange(p_allu_human1_f1_noleg + xlab(NULL),
#                   p_allu_human1_ASV_f1_noleg + ylab(NULL) + xlab(NULL),
#                   p_allu_h1_f1_int_noleg +  coord_cartesian(ylim = c(50,100)),
#                   p_allu_h1_f1_gram_noleg + ylab(NULL),
#                   align = "v", 
#                   ncol = 2, nrow = 2,
#                   heights = c(1, 0.6), widths = c(1,1),
#                   common.legend = FALSE) -> p_allu_human1_f1
# 
# p_allu_human1_f1
# 
# p_allu_human1_f1 %>%
#   export::graph2ppt(append = TRUE,
#                     width = 317.48031496 * 1,
#                     height = 0.618 * 317.48031496 * 1.75 , paper = "A4",  scaling = 2,
#                     file = out_pptx)
```


```r
# ggpubr::ggarrange(p_gram_h1_f1_no_leg + ylab(NULL) + 
#                     theme(axis.title.x=element_blank(),
#                           axis.text.x=element_blank(),
#                           axis.ticks.x=element_blank()),
#                   p_hist_human1_f1_fam_no_leg + 
#                     theme(
#                       strip.background = element_blank(),
#                       strip.text.x = element_blank()
#                     ) +
#                     theme(axis.title.x=element_blank(),
#                           axis.text.x=element_blank(),
#                           axis.ticks.x=element_blank()),
#                   p_hist_human1_f1_strain_no_leg + 
#                     theme(
#                       strip.background = element_blank(),
#                       strip.text.x = element_blank()
#                     ) + ylab(NULL) + 
#                     theme(axis.title.x=element_blank(),
#                           axis.text.x=element_blank(),
#                           axis.ticks.x=element_blank()),
#                   p_int_h1_f1_no_leg + ylab(NULL) + coord_cartesian(ylim = c(50,100)) +
#                     theme(
#                       strip.background = element_blank(),
#                       strip.text.x = element_blank()
#                     ), 
#                   align = "v",
#                   ncol = 1, 
#                   heights = c(0.15, 0.4, 0.8,0.15),
#                   common.legend = FALSE) -> all_plot_h1_f1
# 
# all_plot_h1_f1
# 
# all_plot_h1_f1 %>% 
#   export::graph2ppt(append = TRUE, 
#                     width = 317.48031496 * 3.5, 
#                     height = 0.618 * 317.48031496 * 5 , paper = "A4",  scaling = 2,
#                     file = out_pptx)
# 
# 
# p_gram_h1_f1_leg  %>% 
#   export::graph2ppt(append = TRUE, 
#                     width = 317.48031496 * 1, 
#                     height = 0.618 * 317.48031496 * 1 , paper = "A4",  scaling = 2,
#                     file = out_pptx)
# 
# p_hist_human1_f1_fam_leg  %>% 
#   export::graph2ppt(append = TRUE, 
#                     width = 317.48031496 * 1, 
#                     height = 0.618 * 317.48031496 * 1 , paper = "A4",  scaling = 2,
#                     file = out_pptx)
# 
# p_hist_human1_f1_strain_leg  %>% 
#   export::graph2ppt(append = TRUE, 
#                     width = 317.48031496 * 4, 
#                     height = 0.618 * 317.48031496 * 2 , paper = "A3",  scaling = 2,
#                     file = out_pptx)
# 
# p_int_h1_f1_leg  %>% 
#   export::graph2ppt(append = TRUE, 
#                     width = 317.48031496 * 1, 
#                     height = 0.618 * 317.48031496 * 1 , paper = "A4",  scaling = 2,
#                     file = out_pptx)
# 
# 
# 
# p_gram_h1_f1$data %>% 
#   select(OTU, Sample, Abundance, Treatment, Day_of_Treatment) %>% 
#   xlsx::write.xlsx(file = out_xlsx,
#                    append = TRUE,
#                    sheetName = "p_gram_h1_f1",
#                    showNA = TRUE)
# p_hist_human1_f1_fam$data %>% 
#   select(OTU, Sample, Abundance, Treatment, Day_of_Treatment) %>% 
#   pivot_wider(names_from = OTU, values_from = Abundance) %>% 
#   xlsx::write.xlsx(file = out_xlsx,
#                    append = TRUE,
#                    sheetName = "p_hist_human1_f1_fam",
#                    showNA = TRUE)
# p_hist_human1_f1_strain$data %>% 
#   select(Sample, top ,Abundance, Treatment, Day_of_Treatment) %>% 
#   pivot_wider(names_from = top, values_from = Abundance) %>% 
#   xlsx::write.xlsx(file = out_xlsx,
#                    append = TRUE,
#                    sheetName = "p_hist_human1_f1_strain",
#                    showNA = TRUE)
# p_int_h1_f1$data %>% 
#   select(OTU, Sample, Abundance, Treatment, Day_of_Treatment) %>% 
#   pivot_wider(names_from = OTU, values_from = Abundance) %>% 
#   xlsx::write.xlsx(file = out_xlsx,
#                    append = TRUE,
#                    sheetName = "p_int_h1_f1",
#                    showNA = TRUE)
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
##  [1] ggalluvial_0.12.3    gdtools_0.2.4        microViz_0.9.7      
##  [4] speedyseq_0.5.3.9018 reshape2_1.4.4       scales_1.2.1        
##  [7] compositions_2.0-4   nlme_3.1-159         phyloseq_1.40.0     
## [10] forcats_0.5.2        stringr_1.4.1        dplyr_1.0.10        
## [13] purrr_0.3.5          readr_2.1.3          tidyr_1.2.1         
## [16] tibble_3.1.8         ggplot2_3.3.6        tidyverse_1.3.2     
## 
## loaded via a namespace (and not attached):
##   [1] readxl_1.4.1           uuid_1.1-0             backports_1.4.1       
##   [4] systemfonts_1.0.4      plyr_1.8.7             igraph_1.3.5          
##   [7] splines_4.2.1          GenomeInfoDb_1.32.4    digest_0.6.30         
##  [10] foreach_1.5.2          htmltools_0.5.3        fansi_1.0.3           
##  [13] magrittr_2.0.3         xlsx_0.6.5             googlesheets4_1.0.1   
##  [16] cluster_2.1.4          tzdb_0.3.0             openxlsx_4.2.5        
##  [19] Biostrings_2.64.1      extrafont_0.18         modelr_0.1.9          
##  [22] bayesm_3.1-4           officer_0.4.4          extrafontdb_1.0       
##  [25] colorspace_2.0-3       rvest_1.0.3            ggrepel_0.9.1         
##  [28] textshaping_0.3.6      haven_2.5.1            xfun_0.34             
##  [31] crayon_1.5.2           RCurl_1.98-1.9         jsonlite_1.8.3        
##  [34] survival_3.4-0         iterators_1.0.14       ape_5.6-2             
##  [37] glue_1.6.2             rvg_0.2.5              gtable_0.3.1          
##  [40] gargle_1.2.1           zlibbioc_1.42.0        XVector_0.36.0        
##  [43] car_3.1-0              Rttf2pt1_1.3.10        Rhdf5lib_1.18.2       
##  [46] BiocGenerics_0.42.0    DEoptimR_1.0-11        abind_1.4-5           
##  [49] DBI_1.1.3              rstatix_0.7.0          Rcpp_1.0.9            
##  [52] xtable_1.8-4           stats4_4.2.1           htmlwidgets_1.5.4     
##  [55] httr_1.4.4             ellipsis_0.3.2         rJava_1.0-6           
##  [58] pkgconfig_2.0.3        farver_2.1.1           sass_0.4.2            
##  [61] dbplyr_2.2.1           utf8_1.2.2             here_1.0.1            
##  [64] tidyselect_1.2.0       labeling_0.4.2         rlang_1.0.6           
##  [67] munsell_0.5.0          cellranger_1.1.0       tools_4.2.1           
##  [70] cachem_1.0.6           cli_3.4.1              generics_0.1.3        
##  [73] devEMF_4.1             ade4_1.7-19            export_0.3.0          
##  [76] broom_1.0.1            evaluate_0.17          biomformat_1.24.0     
##  [79] fastmap_1.1.0          yaml_2.3.6             ragg_1.2.3            
##  [82] knitr_1.40             fs_1.5.2               zip_2.2.1             
##  [85] robustbase_0.95-0      rgl_0.110.2            xml2_1.3.3            
##  [88] compiler_4.2.1         rstudioapi_0.14        ggsignif_0.6.4        
##  [91] reprex_2.0.2           bslib_0.4.0            stringi_1.7.8         
##  [94] highr_0.9              stargazer_5.2.3        lattice_0.20-45       
##  [97] Matrix_1.5-1           vegan_2.6-4            microbiome_1.19.1     
## [100] permute_0.9-7          tensorA_0.36.2         multtest_2.52.0       
## [103] vctrs_0.4.2            pillar_1.8.1           lifecycle_1.0.3       
## [106] rhdf5filters_1.8.0     jquerylib_0.1.4        data.table_1.14.4     
## [109] cowplot_1.1.1          bitops_1.0-7           flextable_0.8.2       
## [112] R6_2.5.1               IRanges_2.30.1         codetools_0.2-18      
## [115] MASS_7.3-58.1          assertthat_0.2.1       xlsxjars_0.6.1        
## [118] rhdf5_2.40.0           rprojroot_2.0.3        withr_2.5.0           
## [121] S4Vectors_0.34.0       GenomeInfoDbData_1.2.8 mgcv_1.8-40           
## [124] parallel_4.2.1         hms_1.1.2              grid_4.2.1            
## [127] rmarkdown_2.16         carData_3.0-5          googledrive_2.0.0     
## [130] Rtsne_0.16             ggpubr_0.4.0           Biobase_2.56.0        
## [133] lubridate_1.8.0        base64enc_0.1-3
```
