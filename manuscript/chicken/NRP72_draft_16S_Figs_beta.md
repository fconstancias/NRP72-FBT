---
title: "NTP72 - 16S - beta - chicken draft "
author: "Florentin Constancias"
date: "October 25, 2022"
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
out_pptx = "~/Desktop/16S-chicken-beta.pptx"
out_xlsx = "~/Desktop/16S-chicken-beta.xlsx"
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
  select(SampleID, Day_of_Treatment,Treatment, Reactor, Model2, Model, LibrarySize) 
```

```
##         SampleID Day_of_Treatment     Treatment Reactor   Model2   Model
## 191  p-2-H8-S188               10     CCUG59168     TR4 Chicken2 Chicken
## 6    TR1-16-S199                0   CTX+HV292.1     TR1 Chicken1 Chicken
## 59   TR4-18-S200                2           VAN     TR4 Chicken1 Chicken
## 33   TR2-34-S302               18           CTX     TR2 Chicken1 Chicken
## 48   TR3-28-S162               12       HV292.1     TR3 Chicken1 Chicken
## 5    TR1-15-S168               -1   CTX+HV292.1     TR1 Chicken1 Chicken
## 87   TR6-13-S164               -3     CCUG59168     TR6 Chicken1 Chicken
## 197  p-3-A5-S197               -2 VAN+CCUG59168     TR7 Chicken2 Chicken
## 38   TR3-14-S136               -2       HV292.1     TR3 Chicken1 Chicken
## 151 p-2-D11-S143               13     UNTREATED      CR Chicken2 Chicken
## 95   TR6-21-S107                5     CCUG59168     TR6 Chicken1 Chicken
## 220  p-3-G1-S265               14       HV292.1     TR5 Chicken2 Chicken
## 181  p-2-G9-S177               11           CTX     TR3 Chicken2 Chicken
## 194  p-3-A2-S194               14 VAN+CCUG59168     TR7 Chicken2 Chicken
## 51   TR3-37-S221               21       HV292.1     TR3 Chicken1 Chicken
## 9    TR1-19-S305                3   CTX+HV292.1     TR1 Chicken1 Chicken
## 112  p-1-F11-S71                2     UNTREATED      CR Chicken2 Chicken
## 100  TR6-34-S143               18     CCUG59168     TR6 Chicken1 Chicken
## 146  p-2-C7-S127                8 VAN+CCUG59168     TR7 Chicken2 Chicken
## 37   TR3-13-S304               -3       HV292.1     TR3 Chicken1 Chicken
## 30   TR2-25-S133                9           CTX     TR2 Chicken1 Chicken
## 8    TR1-18-S124                2   CTX+HV292.1     TR1 Chicken1 Chicken
## 26   TR2-19-S318                3           CTX     TR2 Chicken1 Chicken
## 229  p-3-H5-S281               -5           VAN     TR6 Chicken2 Chicken
## 203  p-3-B5-S209               -5     UNTREATED      CR Chicken2 Chicken
## 214  p-3-E5-S245               -5           CTX     TR3 Chicken2 Chicken
## 1       D-0-S154             <NA>         DONOR   DONOR Chicken1 Chicken
## 206  p-3-D2-S230                0   CTX+HV292.1     TR2 Chicken2 Chicken
## 83   TR5-34-S127               18 VAN+CCUG59168     TR5 Chicken1 Chicken
## 198  p-3-A6-S198               -5 VAN+CCUG59168     TR7 Chicken2 Chicken
## 199  p-3-B1-S205               14     UNTREATED      CR Chicken2 Chicken
## 222  p-3-G3-S267               -1       HV292.1     TR5 Chicken2 Chicken
## 216  p-3-F2-S254                0     CCUG59168     TR4 Chicken2 Chicken
## 176  p-2-G4-S172                6           CTX     TR3 Chicken2 Chicken
## 98   TR6-28-S158               12     CCUG59168     TR6 Chicken1 Chicken
## 211  p-3-E2-S242                0           CTX     TR3 Chicken2 Chicken
## 16   TR1-34-S317               18   CTX+HV292.1     TR1 Chicken1 Chicken
## 207  p-3-D3-S231               -1   CTX+HV292.1     TR2 Chicken2 Chicken
## 228  p-3-H4-S280               -2           VAN     TR6 Chicken2 Chicken
## 224  p-3-G5-S269               -5       HV292.1     TR5 Chicken2 Chicken
## 215  p-3-F1-S253               14     CCUG59168     TR4 Chicken2 Chicken
## 200  p-3-B2-S206                0     UNTREATED      CR Chicken2 Chicken
## 65   TR4-28-S183               12           VAN     TR4 Chicken1 Chicken
## 144  p-2-C5-S125                6 VAN+CCUG59168     TR7 Chicken2 Chicken
## 208  p-3-D4-S232               -2   CTX+HV292.1     TR2 Chicken2 Chicken
## 67   TR4-34-S342               18           VAN     TR4 Chicken1 Chicken
## 209  p-3-D5-S233               -5   CTX+HV292.1     TR2 Chicken2 Chicken
## 212  p-3-E3-S243               -1           CTX     TR3 Chicken2 Chicken
## 82   TR5-31-S313               15 VAN+CCUG59168     TR5 Chicken1 Chicken
## 219  p-3-F5-S257               -5     CCUG59168     TR4 Chicken2 Chicken
## 41   TR3-17-S222                1       HV292.1     TR3 Chicken1 Chicken
## 53   TR3-46-S224               30       HV292.1     TR3 Chicken1 Chicken
## 102  TR6-40-S160               24     CCUG59168     TR6 Chicken1 Chicken
## 195  p-3-A3-S195                0 VAN+CCUG59168     TR7 Chicken2 Chicken
## 113  p-1-H10-S94                1   CTX+HV292.1     TR2 Chicken2 Chicken
## 64   TR4-25-S333                9           VAN     TR4 Chicken1 Chicken
## 128 p-2-B11-S119               12           VAN     TR6 Chicken2 Chicken
## 101  TR6-37-S312               21     CCUG59168     TR6 Chicken1 Chicken
## 28   TR2-21-S228                5           CTX     TR2 Chicken1 Chicken
## 217  p-3-F3-S255               -1     CCUG59168     TR4 Chicken2 Chicken
## 185  p-2-H2-S182                4     CCUG59168     TR4 Chicken2 Chicken
## 205  p-3-D1-S229               14   CTX+HV292.1     TR2 Chicken2 Chicken
## 34   TR2-37-S330               21           CTX     TR2 Chicken1 Chicken
## 63   TR4-22-S344                6           VAN     TR4 Chicken1 Chicken
## 25   TR2-18-S326                2           CTX     TR2 Chicken1 Chicken
## 221  p-3-G2-S266                0       HV292.1     TR5 Chicken2 Chicken
## 7    TR1-17-S246                1   CTX+HV292.1     TR1 Chicken1 Chicken
## 223  p-3-G4-S268               -2       HV292.1     TR5 Chicken2 Chicken
## 57   TR4-16-S276                0           VAN     TR4 Chicken1 Chicken
## 110  p-1-E11-S59                1 VAN+CCUG59168     TR7 Chicken2 Chicken
## 111  p-1-F10-S70                1     UNTREATED      CR Chicken2 Chicken
## 39   TR3-15-S338               -1       HV292.1     TR3 Chicken1 Chicken
## 182  p-2-H1-S181                3     CCUG59168     TR4 Chicken2 Chicken
## 94   TR6-20-S132                4     CCUG59168     TR6 Chicken1 Chicken
## 29   TR2-22-S364                6           CTX     TR2 Chicken1 Chicken
## 54   TR4-13-S142               -3           VAN     TR4 Chicken1 Chicken
## 218  p-3-F4-S256               -2     CCUG59168     TR4 Chicken2 Chicken
## 85   TR5-40-S335               24 VAN+CCUG59168     TR5 Chicken1 Chicken
## 114  p-1-H11-S95                2   CTX+HV292.1     TR2 Chicken2 Chicken
## 106  p-1-B11-S23                1     CCUG59168     TR4 Chicken2 Chicken
## 11   TR1-21-S207                5   CTX+HV292.1     TR1 Chicken1 Chicken
## 91   TR6-17-S138                1     CCUG59168     TR6 Chicken1 Chicken
## 50   TR3-34-S354               18       HV292.1     TR3 Chicken1 Chicken
## 226  p-3-H2-S278                0           VAN     TR6 Chicken2 Chicken
## 86   TR5-46-S310               30 VAN+CCUG59168     TR5 Chicken1 Chicken
## 184 p-2-H11-S191               13     CCUG59168     TR4 Chicken2 Chicken
## 188  p-2-H5-S185                7     CCUG59168     TR4 Chicken2 Chicken
## 21   TR2-14-S232               -2           CTX     TR2 Chicken1 Chicken
## 61   TR4-20-S269                4           VAN     TR4 Chicken1 Chicken
## 55   TR4-14-S226               -2           VAN     TR4 Chicken1 Chicken
## 44   TR3-20-S265                4       HV292.1     TR3 Chicken1 Chicken
## 210  p-3-E1-S241               14           CTX     TR3 Chicken2 Chicken
## 131  p-2-B3-S111                4           VAN     TR6 Chicken2 Chicken
## 35   TR2-40-S259               24           CTX     TR2 Chicken1 Chicken
## 149  p-2-D1-S133                3     UNTREATED      CR Chicken2 Chicken
## 187  p-2-H4-S184                6     CCUG59168     TR4 Chicken2 Chicken
## 196  p-3-A4-S196               -1 VAN+CCUG59168     TR7 Chicken2 Chicken
## 17   TR1-37-S306               21   CTX+HV292.1     TR1 Chicken1 Chicken
## 89   TR6-15-S149               -1     CCUG59168     TR6 Chicken1 Chicken
## 204  p-3-C7-S223             <NA>         DONOR   DONOR Chicken2 Chicken
## 56   TR4-15-S218               -1           VAN     TR4 Chicken1 Chicken
## 20   TR2-13-S231               -3           CTX     TR2 Chicken1 Chicken
## 62   TR4-21-S254                5           VAN     TR4 Chicken1 Chicken
## 84   TR5-37-S368               21 VAN+CCUG59168     TR5 Chicken1 Chicken
## 133  p-2-B5-S113                6           VAN     TR6 Chicken2 Chicken
## 108  p-1-C11-S35                1       HV292.1     TR5 Chicken2 Chicken
## 227  p-3-H3-S279               -1           VAN     TR6 Chicken2 Chicken
## 10   TR1-20-S118                4   CTX+HV292.1     TR1 Chicken1 Chicken
## 45   TR3-21-S234                5       HV292.1     TR3 Chicken1 Chicken
## 80   TR5-25-S240                9 VAN+CCUG59168     TR5 Chicken1 Chicken
## 88   TR6-14-S211               -2     CCUG59168     TR6 Chicken1 Chicken
## 134  p-2-B6-S114                7           VAN     TR6 Chicken2 Chicken
## 169  p-2-F8-S164               10   CTX+HV292.1     TR2 Chicken2 Chicken
## 68   TR4-37-S188               21           VAN     TR4 Chicken1 Chicken
## 135  p-2-B7-S115                8           VAN     TR6 Chicken2 Chicken
## 81   TR5-28-S214               12 VAN+CCUG59168     TR5 Chicken1 Chicken
## 130  p-2-B2-S110                3           VAN     TR6 Chicken2 Chicken
## 77   TR5-20-S208                4 VAN+CCUG59168     TR5 Chicken1 Chicken
## 107  p-1-B12-S24                2     CCUG59168     TR4 Chicken2 Chicken
## 132  p-2-B4-S112                5           VAN     TR6 Chicken2 Chicken
## 138  p-2-C1-S121                2 VAN+CCUG59168     TR7 Chicken2 Chicken
## 13   TR1-25-S117                9   CTX+HV292.1     TR1 Chicken1 Chicken
## 168  p-2-F7-S163                9   CTX+HV292.1     TR2 Chicken2 Chicken
## 93   TR6-19-S137                3     CCUG59168     TR6 Chicken1 Chicken
## 150 p-2-D10-S142               12     UNTREATED      CR Chicken2 Chicken
## 70   TR4-46-S215               30           VAN     TR4 Chicken1 Chicken
## 90   TR6-16-S371                0     CCUG59168     TR6 Chicken1 Chicken
## 52   TR3-40-S257               24       HV292.1     TR3 Chicken1 Chicken
## 71   TR5-13-S103               -3 VAN+CCUG59168     TR5 Chicken1 Chicken
## 36   TR2-46-S239               30           CTX     TR2 Chicken1 Chicken
## 166  p-2-F5-S161                7   CTX+HV292.1     TR2 Chicken2 Chicken
## 32   TR2-31-S182               15           CTX     TR2 Chicken1 Chicken
## 192  p-2-H9-S189               11     CCUG59168     TR4 Chicken2 Chicken
## 99   TR6-31-S272               15     CCUG59168     TR6 Chicken1 Chicken
## 73   TR5-16-S290                0 VAN+CCUG59168     TR5 Chicken1 Chicken
## 79   TR5-22-S145                6 VAN+CCUG59168     TR5 Chicken1 Chicken
## 162 p-2-F11-S167               13   CTX+HV292.1     TR2 Chicken2 Chicken
## 47   TR3-25-S102                9       HV292.1     TR3 Chicken1 Chicken
## 154  p-2-D4-S136                6     UNTREATED      CR Chicken2 Chicken
## 124  p-2-A8-S104                9       HV292.1     TR5 Chicken2 Chicken
## 143  p-2-C4-S124                5 VAN+CCUG59168     TR7 Chicken2 Chicken
## 159  p-2-D9-S141               11     UNTREATED      CR Chicken2 Chicken
## 60   TR4-19-S270                3           VAN     TR4 Chicken1 Chicken
## 92   TR6-18-S262                2     CCUG59168     TR6 Chicken1 Chicken
## 105  p-1-A12-S12                2           CTX     TR3 Chicken2 Chicken
## 75   TR5-18-S201                2 VAN+CCUG59168     TR5 Chicken1 Chicken
## 153  p-2-D3-S135                5     UNTREATED      CR Chicken2 Chicken
## 58   TR4-17-S114                1           VAN     TR4 Chicken1 Chicken
## 161 p-2-F10-S166               12   CTX+HV292.1     TR2 Chicken2 Chicken
## 158  p-2-D8-S140               10     UNTREATED      CR Chicken2 Chicken
## 155  p-2-D5-S137                7     UNTREATED      CR Chicken2 Chicken
## 225  p-3-H1-S277               14           VAN     TR6 Chicken2 Chicken
## 109  p-1-D11-S47                1           VAN     TR6 Chicken2 Chicken
## 136  p-2-B8-S116                9           VAN     TR6 Chicken2 Chicken
## 193  p-3-A1-S193               13 VAN+CCUG59168     TR7 Chicken2 Chicken
## 142  p-2-C3-S123                4 VAN+CCUG59168     TR7 Chicken2 Chicken
## 18   TR1-40-S209               24   CTX+HV292.1     TR1 Chicken1 Chicken
## 140 p-2-C11-S131               12 VAN+CCUG59168     TR7 Chicken2 Chicken
## 170  p-2-F9-S165               11   CTX+HV292.1     TR2 Chicken2 Chicken
## 126  p-2-B1-S109                2           VAN     TR6 Chicken2 Chicken
## 127 p-2-B10-S118               11           VAN     TR6 Chicken2 Chicken
## 46   TR3-22-S235                6       HV292.1     TR3 Chicken1 Chicken
## 165  p-2-F4-S160                6   CTX+HV292.1     TR2 Chicken2 Chicken
## 3    TR1-13-S171               -3   CTX+HV292.1     TR1 Chicken1 Chicken
## 147  p-2-C8-S128                9 VAN+CCUG59168     TR7 Chicken2 Chicken
## 31   TR2-28-S286               12           CTX     TR2 Chicken1 Chicken
## 66   TR4-31-S278               15           VAN     TR4 Chicken1 Chicken
## 2    IR1-35-S370               19   CTX+HV292.1     TR1 Chicken1 Chicken
## 42   TR3-18-S293                2       HV292.1     TR3 Chicken1 Chicken
## 139 p-2-C10-S130               11 VAN+CCUG59168     TR7 Chicken2 Chicken
## 14   TR1-28-S289               12   CTX+HV292.1     TR1 Chicken1 Chicken
## 22   TR2-15-S288               -1           CTX     TR2 Chicken1 Chicken
## 167  p-2-F6-S162                8   CTX+HV292.1     TR2 Chicken2 Chicken
## 164  p-2-F3-S159                5   CTX+HV292.1     TR2 Chicken2 Chicken
## 156  p-2-D6-S138                8     UNTREATED      CR Chicken2 Chicken
## 27   TR2-20-S271                4           CTX     TR2 Chicken1 Chicken
## 24   TR2-17-S248                1           CTX     TR2 Chicken1 Chicken
## 118 p-2-A12-S108               13       HV292.1     TR5 Chicken2 Chicken
## 202  p-3-B4-S208               -2     UNTREATED      CR Chicken2 Chicken
## 160  p-2-F1-S157                3   CTX+HV292.1     TR2 Chicken2 Chicken
## 152  p-2-D2-S134                4     UNTREATED      CR Chicken2 Chicken
## 19   TR1-46-S170               30   CTX+HV292.1     TR1 Chicken1 Chicken
## 43   TR3-19-S266                3       HV292.1     TR3 Chicken1 Chicken
## 72   TR5-14-S267               -2 VAN+CCUG59168     TR5 Chicken1 Chicken
## 122  p-2-A6-S102                7       HV292.1     TR5 Chicken2 Chicken
## 120   p-2-A3-S99                4       HV292.1     TR5 Chicken2 Chicken
## 69   TR4-40-S242               24           VAN     TR4 Chicken1 Chicken
## 186  p-2-H3-S183                5     CCUG59168     TR4 Chicken2 Chicken
## 148  p-2-C9-S129               10 VAN+CCUG59168     TR7 Chicken2 Chicken
## 183 p-2-H10-S190               12     CCUG59168     TR4 Chicken2 Chicken
## 137  p-2-B9-S117               10           VAN     TR6 Chicken2 Chicken
## 78   TR5-21-S185                5 VAN+CCUG59168     TR5 Chicken1 Chicken
## 172 p-2-G10-S178               12           CTX     TR3 Chicken2 Chicken
## 119   p-2-A2-S98                3       HV292.1     TR5 Chicken2 Chicken
## 12   TR1-22-S282                6   CTX+HV292.1     TR1 Chicken1 Chicken
## 201  p-3-B3-S207               -1     UNTREATED      CR Chicken2 Chicken
## 116 p-2-A10-S106               11       HV292.1     TR5 Chicken2 Chicken
## 76    TR5-19-S87                3 VAN+CCUG59168     TR5 Chicken1 Chicken
## 145  p-2-C6-S126                7 VAN+CCUG59168     TR7 Chicken2 Chicken
## 173 p-2-G11-S179               13           CTX     TR3 Chicken2 Chicken
## 179  p-2-G7-S175                9           CTX     TR3 Chicken2 Chicken
## 163  p-2-F2-S158                4   CTX+HV292.1     TR2 Chicken2 Chicken
## 23   TR2-16-S184                0           CTX     TR2 Chicken1 Chicken
## 177  p-2-G5-S173                7           CTX     TR3 Chicken2 Chicken
## 129 p-2-B12-S120               13           VAN     TR6 Chicken2 Chicken
## 141  p-2-C2-S122                3 VAN+CCUG59168     TR7 Chicken2 Chicken
## 123  p-2-A7-S103                8       HV292.1     TR5 Chicken2 Chicken
## 174  p-2-G2-S170                4           CTX     TR3 Chicken2 Chicken
## 175  p-2-G3-S171                5           CTX     TR3 Chicken2 Chicken
## 178  p-2-G6-S174                8           CTX     TR3 Chicken2 Chicken
## 4    TR1-14-S274               -2   CTX+HV292.1     TR1 Chicken1 Chicken
## 189  p-2-H6-S186                8     CCUG59168     TR4 Chicken2 Chicken
## 104  p-1-A11-S11                1           CTX     TR3 Chicken2 Chicken
## 125  p-2-A9-S105               10       HV292.1     TR5 Chicken2 Chicken
## 190  p-2-H7-S187                9     CCUG59168     TR4 Chicken2 Chicken
## 97   TR6-25-S378                9     CCUG59168     TR6 Chicken1 Chicken
## 180  p-2-G8-S176               10           CTX     TR3 Chicken2 Chicken
## 213  p-3-E4-S244               -2           CTX     TR3 Chicken2 Chicken
## 121  p-2-A4-S100                5       HV292.1     TR5 Chicken2 Chicken
## 117 p-2-A11-S107               12       HV292.1     TR5 Chicken2 Chicken
## 157  p-2-D7-S139                9     UNTREATED      CR Chicken2 Chicken
## 96   TR6-22-S152                6     CCUG59168     TR6 Chicken1 Chicken
## 171  p-2-G1-S169                3           CTX     TR3 Chicken2 Chicken
## 49   TR3-31-S283               15       HV292.1     TR3 Chicken1 Chicken
## 103  TR6-46-S176               30     CCUG59168     TR6 Chicken1 Chicken
## 15   TR1-31-S153               15   CTX+HV292.1     TR1 Chicken1 Chicken
## 40   TR3-16-S296                0       HV292.1     TR3 Chicken1 Chicken
## 115   p-2-A1-S97                2       HV292.1     TR5 Chicken2 Chicken
## 74    TR5-17-S82                1 VAN+CCUG59168     TR5 Chicken1 Chicken
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

![](NRP72_draft_16S_Figs_beta_files/figure-html/unnamed-chunk-8-1.png)<!-- -->

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

![](NRP72_draft_16S_Figs_beta_files/figure-html/unnamed-chunk-9-1.png)<!-- -->


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
  facet_grid(Model2 ~ distance, scales = "free") -> pcoa_all

pcoa_all + theme(legend.position = "none") -> pcoa_all_noleg

pcoa_all  %>%  ggpubr::get_legend() %>% ggpubr::as_ggplot() -> pcoa_all_leg

pcoa_all_noleg
```

![](NRP72_draft_16S_Figs_beta_files/figure-html/unnamed-chunk-10-1.png)<!-- -->



```r
pcoa_all_noleg
```

![](NRP72_draft_16S_Figs_beta_files/figure-html/unnamed-chunk-11-1.png)<!-- -->


```r
pcoa_all_noleg %>%
  export::graph2ppt(append = TRUE,
                    width = 317.48031496 * 1,
                    height = 0.618 * 317.48031496 * 1 , paper = "A4",  scaling = 2,
                    file = out_pptx)
```

```
## Exported graph as ~/Desktop/16S-chicken-beta.pptx
```

```r
pcoa_all_leg %>%
  export::graph2ppt(append = TRUE,
                    width = 317.48031496 * 1,
                    height = 0.618 * 317.48031496 * 1 , paper = "A4",  scaling = 2,
                    file = out_pptx)
```

```
## Exported graph as ~/Desktop/16S-chicken-beta.pptx
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
## Exported table as ~/Desktop/16S-chicken-beta.pptx
```

```{=html}
<template id="da7e349f-06e9-434d-be3d-555747416ca3"><style>
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
</style><div class="tabwid"><style>.cl-ca31150c{}.cl-ca2bbc7e{font-family:'Helvetica';font-size:12pt;font-weight:bold;font-style:normal;text-decoration:none;color:rgba(0, 0, 0, 1.00);background-color:transparent;}.cl-ca2bbc92{font-family:'Helvetica';font-size:12pt;font-weight:normal;font-style:normal;text-decoration:none;color:rgba(0, 0, 0, 1.00);background-color:transparent;}.cl-ca2e37ce{margin:0;text-align:left;border-bottom: 0 solid rgba(0, 0, 0, 1.00);border-top: 0 solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);padding-bottom:5pt;padding-top:5pt;padding-left:5pt;padding-right:5pt;line-height: 1;background-color:transparent;}.cl-ca2e449e{width:2.14in;background-color:transparent;vertical-align: middle;border-bottom: 1pt solid rgba(0, 0, 0, 1.00);border-top: 1pt solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-ca2e44a8{width:2.14in;background-color:transparent;vertical-align: middle;border-bottom: 0 solid rgba(0, 0, 0, 1.00);border-top: 0 solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-ca2e44b2{width:2.14in;background-color:transparent;vertical-align: middle;border-bottom: 1pt solid rgba(0, 0, 0, 1.00);border-top: 0 solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}</style><table class='cl-ca31150c'><thead><tr style="overflow-wrap:break-word;"><td class="cl-ca2e449e"><p class="cl-ca2e37ce"><span class="cl-ca2bbc7e">.id...1</span></p></td><td class="cl-ca2e449e"><p class="cl-ca2e37ce"><span class="cl-ca2bbc7e">V1...2</span></p></td><td class="cl-ca2e449e"><p class="cl-ca2e37ce"><span class="cl-ca2bbc7e">.id...3</span></p></td><td class="cl-ca2e449e"><p class="cl-ca2e37ce"><span class="cl-ca2bbc7e">V1...4</span></p></td></tr></thead><tbody><tr style="overflow-wrap:break-word;"><td class="cl-ca2e44a8"><p class="cl-ca2e37ce"><span class="cl-ca2bbc92">bjaccard</span></p></td><td class="cl-ca2e44a8"><p class="cl-ca2e37ce"><span class="cl-ca2bbc92">Axis.1   [25.7%]</span></p></td><td class="cl-ca2e44a8"><p class="cl-ca2e37ce"><span class="cl-ca2bbc92">bjaccard</span></p></td><td class="cl-ca2e44a8"><p class="cl-ca2e37ce"><span class="cl-ca2bbc92">Axis.2   [20.6%]</span></p></td></tr><tr style="overflow-wrap:break-word;"><td class="cl-ca2e44b2"><p class="cl-ca2e37ce"><span class="cl-ca2bbc92">aichinson</span></p></td><td class="cl-ca2e44b2"><p class="cl-ca2e37ce"><span class="cl-ca2bbc92">Axis.1   [27.1%]</span></p></td><td class="cl-ca2e44b2"><p class="cl-ca2e37ce"><span class="cl-ca2bbc92">aichinson</span></p></td><td class="cl-ca2e44b2"><p class="cl-ca2e37ce"><span class="cl-ca2bbc92">Axis.2   [19.6%]</span></p></td></tr></tbody></table></div></template>
<div class="flextable-shadow-host" id="1b5c1732-6084-4387-959e-7f99e7fb96e0"></div>
<script>
var dest = document.getElementById("1b5c1732-6084-4387-959e-7f99e7fb96e0");
var template = document.getElementById("da7e349f-06e9-434d-be3d-555747416ca3");
var caption = template.content.querySelector("caption");
var fantome = dest.attachShadow({mode: 'open'});
var templateContent = template.content;
fantome.appendChild(templateContent);
</script>

```

### chicken1:


```r
ps_rare %>%
  subset_samples(Model2 == "Chicken1" & Treatment != "DONOR") -> ps_tmp_forbdiv
```

#### NMDS:


```r
ps_tmp_forbdiv %>% 
  # microbiome::transform(transform = "compositional") %>% 
  phyloseq_plot_bdiv(dlist = dist,
                     m = "NMDS",
                     seed = 123,
                     axis1 = 1,
                     axis2 = 2) -> NMDS_tmp
```

#### aichinson


```r
NMDS_tmp$aichinson$layers[[1]] = NULL; NMDS_tmp$aichinson$layers[[1]] = NULL
NMDS_tmp$aichinson$layers[[2]] = NULL; NMDS_tmp$aichinson$layers[[1]] = NULL

# plots_hall_humans$aichinson$layers[[1]] = NULL;plots_hall_humans$aichinson$layers[[1]] = NULL
# plots_hall_humans$aichinson$layers[[2]] = NULL;plots_hall_humans$aichinson$layers[[2]] = NULL


NMDS_tmp$aichinson + geom_point(size = 3,
                                aes(colour = Treatment ,
                                    shape = Antibiotic_mg.L,
                                    alpha = NULL)) +
  geom_path(data = NMDS_tmp$aichinson$data %>%
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
  facet_grid(. ~ Model2, scales = "free") -> NMDS_tmp2

NMDS_tmp2
```

![](NRP72_draft_16S_Figs_beta_files/figure-html/unnamed-chunk-16-1.png)<!-- -->


```r
NMDS_tmp2 + 
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
<div id="htmlwidget-f30d7598595affc4f25a" style="width:100%;height:auto;" class="datatables html-widget"></div>
<script type="application/json" data-for="htmlwidget-f30d7598595affc4f25a">{"x":{"filter":"none","vertical":false,"data":[["ASV0001","ASV0008","ASV0006","ASV0134","ASV0125","ASV0165","ASV0187","ASV0278"],[0.179126285157017,-0.417317172405722,-0.520111335186489,-0.469901977442166,-0.568941471660472,0.549879695939574,-0.617503745304176,-0.620583212859231],[0.760825246934185,-0.768280896635797,-0.722683637195785,-0.737611838158516,0.420534229374111,0.462664680607074,0.42291115640916,0.417913549423581],[0.001,0.001,0.001,0.001,0.001,0.001,0.001,0.001],[0.610941282406617,0.764409158520211,0.792787440460002,0.764879092195644,0.500543436250461,0.516426286687844,0.560164721680017,0.559775258874501],["Bacteria","Bacteria","Bacteria","Bacteria","Bacteria","Bacteria","Bacteria","Bacteria"],["Proteobacteria","Firmicutes","Firmicutes","Firmicutes","Firmicutes","Firmicutes","Firmicutes","Firmicutes"],["Gammaproteobacteria","Clostridia","Clostridia","Clostridia","Bacilli","Clostridia","Clostridia","Clostridia"],["Enterobacterales","Oscillospirales","Lachnospirales","Oscillospirales","Erysipelotrichales","Peptostreptococcales-Tissierellales","Christensenellales","Monoglobales"],["Morganellaceae","Ruminococcaceae","Lachnospiraceae","Oscillospiraceae","Erysipelatoclostridiaceae","Anaerovoracaceae","Christensenellaceae","Monoglobaceae"],[null,null,null,null,null,null,null,null],[null,null,null,null,null,null,null,null],[null,null,null,null,null,null,null,null],[0.00263157894736842,0.00263157894736842,0.00263157894736842,0.00263157894736842,0.00263157894736842,0.00263157894736842,0.00263157894736842,0.00263157894736842]],"container":"<table class=\"display\">\n  <thead>\n    <tr>\n      <th> <\/th>\n      <th>Axis.1<\/th>\n      <th>Axis.2<\/th>\n      <th>pval<\/th>\n      <th>r<\/th>\n      <th>Kingdom<\/th>\n      <th>Phylum<\/th>\n      <th>Class<\/th>\n      <th>Order<\/th>\n      <th>tax_rank_plot<\/th>\n      <th>Genus<\/th>\n      <th>Species<\/th>\n      <th>Strain<\/th>\n      <th>pval.adj<\/th>\n    <\/tr>\n  <\/thead>\n<\/table>","options":{"columnDefs":[{"className":"dt-right","targets":[1,2,3,4,13]},{"orderable":false,"targets":0}],"order":[],"autoWidth":false,"orderClasses":false}},"evals":[],"jsHooks":[]}</script>
```


```r
out$plot -> NMDS_tmp2_fam

NMDS_tmp2 + 
  theme(legend.position = "none") -> NMDS_tmp2_noleg

NMDS_tmp2_noleg
```

![](NRP72_draft_16S_Figs_beta_files/figure-html/unnamed-chunk-18-1.png)<!-- -->

```r
NMDS_tmp2  %>%  ggpubr::get_legend() %>% ggpubr::as_ggplot() -> NMDS_tmp2_leg

NMDS_tmp2_noleg %>% 
  export::graph2ppt(append = TRUE,
                    width = 317.48031496 * 1,
                    height = 0.618 * 317.48031496 * 1 , paper = "A4",  scaling = 2,
                    file = out_pptx)

NMDS_tmp2_leg %>% 
  export::graph2ppt(append = TRUE,
                    width = 317.48031496 * 1,
                    height = 0.618 * 317.48031496 * 1 , paper = "A4",  scaling = 2,
                    file = out_pptx)

NMDS_tmp2_fam %>% 
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
<div id="htmlwidget-b85bf888dc62471906a9" style="width:100%;height:auto;" class="datatables html-widget"></div>
<script type="application/json" data-for="htmlwidget-b85bf888dc62471906a9">{"x":{"filter":"none","vertical":false,"data":[["1","2","3","4","5","6"],["Succinat_mM","Acetat_mM","Propionat_mM","Butyrat_mM","Valerat_mM","Total_SCFA_mM"],[0.7987843132071,-0.763822131606401,0.27159917848853,-0.61519382446119,-0.639755304442811,-0.811340602560074],[0.178331284674618,-0.400906971626607,0.683522424509029,-0.709613437960205,0.412493789226407,-0.389815834112515],["Succinat_mM","Acetat_mM","Propionat_mM","Butyrat_mM","Valerat_mM","Total_SCFA_mM"],[0.001,0.001,0.001,0.001,0.001,0.001],[0.669858426119437,0.744150648630563,0.540969018562345,0.882014672988888,0.579437975713073,0.810229957887379]],"container":"<table class=\"display\">\n  <thead>\n    <tr>\n      <th> <\/th>\n      <th>label<\/th>\n      <th>Axis.1<\/th>\n      <th>Axis.2<\/th>\n      <th>id<\/th>\n      <th>pval<\/th>\n      <th>r<\/th>\n    <\/tr>\n  <\/thead>\n<\/table>","options":{"columnDefs":[{"className":"dt-right","targets":[2,3,5,6]},{"orderable":false,"targets":0}],"order":[],"autoWidth":false,"orderClasses":false}},"evals":[],"jsHooks":[]}</script>
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
NMDS_tmp$bjaccard$layers[[1]] = NULL; NMDS_tmp$bjaccard$layers[[1]] = NULL
NMDS_tmp$bjaccard$layers[[2]] = NULL; NMDS_tmp$bjaccard$layers[[1]] = NULL

# plots_hall_humans$aichinson$layers[[1]] = NULL;plots_hall_humans$aichinson$layers[[1]] = NULL
# plots_hall_humans$aichinson$layers[[2]] = NULL;plots_hall_humans$aichinson$layers[[2]] = NULL


NMDS_tmp$bjaccard + geom_point(size = 3,
                               aes(colour = Treatment ,
                                   shape = Antibiotic_mg.L,
                                   alpha = NULL)) +
  geom_path(data = NMDS_tmp$bjaccard$data %>%
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
  facet_grid(. ~ Model2, scales = "free") -> NMDS_tmp2

NMDS_tmp2
```

![](NRP72_draft_16S_Figs_beta_files/figure-html/unnamed-chunk-21-1.png)<!-- -->


```r
NMDS_tmp2 + 
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
<div id="htmlwidget-88fdd2a34d502a66744b" style="width:100%;height:auto;" class="datatables html-widget"></div>
<script type="application/json" data-for="htmlwidget-88fdd2a34d502a66744b">{"x":{"filter":"none","vertical":false,"data":[["ASV0001","ASV0008","ASV0006","ASV0134","ASV0109","ASV0165","ASV0187","ASV0278"],[0.231618167499601,-0.463549587310127,-0.576476302781683,-0.525763548777435,-0.344570047491009,0.649893612073232,-0.473121554133513,-0.427436772742712],[-0.715270448514504,0.710032415965814,0.671084322176517,0.649010202443495,-0.589108997800361,-0.375615881402725,-0.610050545663667,-0.606103099393499],[0.001,0.001,0.001,0.001,0.001,0.001,0.001,0.001],[0.565258790034014,0.71902425161764,0.782679095139955,0.697641552098789,0.465777928917301,0.563448997375539,0.596005673250249,0.55006316178711],["Bacteria","Bacteria","Bacteria","Bacteria","Bacteria","Bacteria","Bacteria","Bacteria"],["Proteobacteria","Firmicutes","Firmicutes","Firmicutes","Firmicutes","Firmicutes","Firmicutes","Firmicutes"],["Gammaproteobacteria","Clostridia","Clostridia","Clostridia","Clostridia","Clostridia","Clostridia","Clostridia"],["Enterobacterales","Oscillospirales","Lachnospirales","Oscillospirales","Lachnospirales","Peptostreptococcales-Tissierellales","Christensenellales","Monoglobales"],["Morganellaceae","Ruminococcaceae","Lachnospiraceae","Oscillospiraceae","Defluviitaleaceae","Anaerovoracaceae","Christensenellaceae","Monoglobaceae"],[null,null,null,null,null,null,null,null],[null,null,null,null,null,null,null,null],[null,null,null,null,null,null,null,null],[0.00263157894736842,0.00263157894736842,0.00263157894736842,0.00263157894736842,0.00263157894736842,0.00263157894736842,0.00263157894736842,0.00263157894736842]],"container":"<table class=\"display\">\n  <thead>\n    <tr>\n      <th> <\/th>\n      <th>Axis.1<\/th>\n      <th>Axis.2<\/th>\n      <th>pval<\/th>\n      <th>r<\/th>\n      <th>Kingdom<\/th>\n      <th>Phylum<\/th>\n      <th>Class<\/th>\n      <th>Order<\/th>\n      <th>tax_rank_plot<\/th>\n      <th>Genus<\/th>\n      <th>Species<\/th>\n      <th>Strain<\/th>\n      <th>pval.adj<\/th>\n    <\/tr>\n  <\/thead>\n<\/table>","options":{"columnDefs":[{"className":"dt-right","targets":[1,2,3,4,13]},{"orderable":false,"targets":0}],"order":[],"autoWidth":false,"orderClasses":false}},"evals":[],"jsHooks":[]}</script>
```


```r
out$plot -> NMDS_tmp2_fam

NMDS_tmp2 + 
  theme(legend.position = "none") -> NMDS_tmp2_noleg

NMDS_tmp2_noleg
```

![](NRP72_draft_16S_Figs_beta_files/figure-html/unnamed-chunk-23-1.png)<!-- -->

```r
NMDS_tmp2  %>%  ggpubr::get_legend() %>% ggpubr::as_ggplot() -> NMDS_tmp2_leg

NMDS_tmp2_noleg %>% 
  export::graph2ppt(append = TRUE,
                    width = 317.48031496 * 1,
                    height = 0.618 * 317.48031496 * 1 , paper = "A4",  scaling = 2,
                    file = out_pptx)

NMDS_tmp2_leg %>% 
  export::graph2ppt(append = TRUE,
                    width = 317.48031496 * 1,
                    height = 0.618 * 317.48031496 * 1 , paper = "A4",  scaling = 2,
                    file = out_pptx)

NMDS_tmp2_fam %>% 
  export::graph2ppt(append = TRUE,
                    width = 317.48031496 * 1,
                    height = 0.618 * 317.48031496 * 1 , paper = "A4",  scaling = 2,
                    file = out_pptx)
```


```r
ps_rare %>%
  subset_samples(Model2 == "Chicken1" & Treatment != "DONOR") %>%
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
<div id="htmlwidget-8e1fe96c0688a849169a" style="width:100%;height:auto;" class="datatables html-widget"></div>
<script type="application/json" data-for="htmlwidget-8e1fe96c0688a849169a">{"x":{"filter":"none","vertical":false,"data":[["1","2","3","4","5","6"],["Succinat_mM","Acetat_mM","Propionat_mM","Butyrat_mM","Valerat_mM","Total_SCFA_mM"],[0.805376707136237,-0.862698523914289,0.324337009107315,-0.717381936893163,-0.465080250779239,-0.90088101316165],[0.13593975442181,0.201311834400704,-0.660808679152135,0.590096643093981,-0.607635844076233,0.135182677818885],["Succinat_mM","Acetat_mM","Propionat_mM","Butyrat_mM","Valerat_mM","Total_SCFA_mM"],[0.001,0.001,0.001,0.001,0.001,0.001],[0.66711125722987,0.784775197833669,0.541862605919468,0.862850891571372,0.585520958671116,0.829860956257444]],"container":"<table class=\"display\">\n  <thead>\n    <tr>\n      <th> <\/th>\n      <th>label<\/th>\n      <th>Axis.1<\/th>\n      <th>Axis.2<\/th>\n      <th>id<\/th>\n      <th>pval<\/th>\n      <th>r<\/th>\n    <\/tr>\n  <\/thead>\n<\/table>","options":{"columnDefs":[{"className":"dt-right","targets":[2,3,5,6]},{"orderable":false,"targets":0}],"order":[],"autoWidth":false,"orderClasses":false}},"evals":[],"jsHooks":[]}</script>
```


```r
env_out$plot %>% 
  export::graph2ppt(append = TRUE,
                    width = 317.48031496 * 1,
                    height = 0.618 * 317.48031496 * 1 , paper = "A4",  scaling = 2,
                    file = out_pptx)
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

![](NRP72_draft_16S_Figs_beta_files/figure-html/unnamed-chunk-27-1.png)<!-- -->


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
<div id="htmlwidget-bb3486349adfa5d821dc" style="width:100%;height:auto;" class="datatables html-widget"></div>
<script type="application/json" data-for="htmlwidget-bb3486349adfa5d821dc">{"x":{"filter":"none","vertical":false,"data":[["ASV0001","ASV0008","ASV0006","ASV0134","ASV0125","ASV0165","ASV0187","ASV0278"],[0.179126285157017,-0.417317172405722,-0.520111335186489,-0.469901977442166,-0.568941471660472,0.549879695939574,-0.617503745304176,-0.620583212859231],[0.760825246934185,-0.768280896635797,-0.722683637195785,-0.737611838158516,0.420534229374111,0.462664680607074,0.42291115640916,0.417913549423581],[0.001,0.001,0.001,0.001,0.001,0.001,0.001,0.001],[0.610941282406617,0.764409158520211,0.792787440460002,0.764879092195644,0.500543436250461,0.516426286687844,0.560164721680017,0.559775258874501],["Bacteria","Bacteria","Bacteria","Bacteria","Bacteria","Bacteria","Bacteria","Bacteria"],["Proteobacteria","Firmicutes","Firmicutes","Firmicutes","Firmicutes","Firmicutes","Firmicutes","Firmicutes"],["Gammaproteobacteria","Clostridia","Clostridia","Clostridia","Bacilli","Clostridia","Clostridia","Clostridia"],["Enterobacterales","Oscillospirales","Lachnospirales","Oscillospirales","Erysipelotrichales","Peptostreptococcales-Tissierellales","Christensenellales","Monoglobales"],["Morganellaceae","Ruminococcaceae","Lachnospiraceae","Oscillospiraceae","Erysipelatoclostridiaceae","Anaerovoracaceae","Christensenellaceae","Monoglobaceae"],[null,null,null,null,null,null,null,null],[null,null,null,null,null,null,null,null],[null,null,null,null,null,null,null,null],[0.00263157894736842,0.00263157894736842,0.00263157894736842,0.00263157894736842,0.00263157894736842,0.00263157894736842,0.00263157894736842,0.00263157894736842]],"container":"<table class=\"display\">\n  <thead>\n    <tr>\n      <th> <\/th>\n      <th>Axis.1<\/th>\n      <th>Axis.2<\/th>\n      <th>pval<\/th>\n      <th>r<\/th>\n      <th>Kingdom<\/th>\n      <th>Phylum<\/th>\n      <th>Class<\/th>\n      <th>Order<\/th>\n      <th>tax_rank_plot<\/th>\n      <th>Genus<\/th>\n      <th>Species<\/th>\n      <th>Strain<\/th>\n      <th>pval.adj<\/th>\n    <\/tr>\n  <\/thead>\n<\/table>","options":{"columnDefs":[{"className":"dt-right","targets":[1,2,3,4,13]},{"orderable":false,"targets":0}],"order":[],"autoWidth":false,"orderClasses":false}},"evals":[],"jsHooks":[]}</script>
```


```r
out$plot -> PCoA_tmp2_fam

PCoA_tmp2 + 
  theme(legend.position = "none") -> PCoA_tmp2_noleg

PCoA_tmp2_noleg
```

![](NRP72_draft_16S_Figs_beta_files/figure-html/unnamed-chunk-29-1.png)<!-- -->

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
<div id="htmlwidget-44474b3e8cd160f85a98" style="width:100%;height:auto;" class="datatables html-widget"></div>
<script type="application/json" data-for="htmlwidget-44474b3e8cd160f85a98">{"x":{"filter":"none","vertical":false,"data":[["1","2","3","4","5","6"],["Succinat_mM","Acetat_mM","Propionat_mM","Butyrat_mM","Valerat_mM","Total_SCFA_mM"],[0.7987843132071,-0.763822131606401,0.27159917848853,-0.61519382446119,-0.639755304442811,-0.811340602560074],[0.178331284674618,-0.400906971626607,0.683522424509029,-0.709613437960205,0.412493789226407,-0.389815834112515],["Succinat_mM","Acetat_mM","Propionat_mM","Butyrat_mM","Valerat_mM","Total_SCFA_mM"],[0.001,0.001,0.001,0.001,0.001,0.001],[0.669858426119437,0.744150648630563,0.540969018562345,0.882014672988888,0.579437975713073,0.810229957887379]],"container":"<table class=\"display\">\n  <thead>\n    <tr>\n      <th> <\/th>\n      <th>label<\/th>\n      <th>Axis.1<\/th>\n      <th>Axis.2<\/th>\n      <th>id<\/th>\n      <th>pval<\/th>\n      <th>r<\/th>\n    <\/tr>\n  <\/thead>\n<\/table>","options":{"columnDefs":[{"className":"dt-right","targets":[2,3,5,6]},{"orderable":false,"targets":0}],"order":[],"autoWidth":false,"orderClasses":false}},"evals":[],"jsHooks":[]}</script>
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

![](NRP72_draft_16S_Figs_beta_files/figure-html/unnamed-chunk-32-1.png)<!-- -->


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
<div id="htmlwidget-2514ce28ae593bde0f73" style="width:100%;height:auto;" class="datatables html-widget"></div>
<script type="application/json" data-for="htmlwidget-2514ce28ae593bde0f73">{"x":{"filter":"none","vertical":false,"data":[["ASV0001","ASV0008","ASV0006","ASV0134","ASV0109","ASV0165","ASV0187","ASV0278"],[0.231618167499601,-0.463549587310127,-0.576476302781683,-0.525763548777435,-0.344570047491009,0.649893612073232,-0.473121554133513,-0.427436772742712],[-0.715270448514504,0.710032415965814,0.671084322176517,0.649010202443495,-0.589108997800361,-0.375615881402725,-0.610050545663667,-0.606103099393499],[0.001,0.001,0.001,0.001,0.001,0.001,0.001,0.001],[0.565258790034014,0.71902425161764,0.782679095139955,0.697641552098789,0.465777928917301,0.563448997375539,0.596005673250249,0.55006316178711],["Bacteria","Bacteria","Bacteria","Bacteria","Bacteria","Bacteria","Bacteria","Bacteria"],["Proteobacteria","Firmicutes","Firmicutes","Firmicutes","Firmicutes","Firmicutes","Firmicutes","Firmicutes"],["Gammaproteobacteria","Clostridia","Clostridia","Clostridia","Clostridia","Clostridia","Clostridia","Clostridia"],["Enterobacterales","Oscillospirales","Lachnospirales","Oscillospirales","Lachnospirales","Peptostreptococcales-Tissierellales","Christensenellales","Monoglobales"],["Morganellaceae","Ruminococcaceae","Lachnospiraceae","Oscillospiraceae","Defluviitaleaceae","Anaerovoracaceae","Christensenellaceae","Monoglobaceae"],[null,null,null,null,null,null,null,null],[null,null,null,null,null,null,null,null],[null,null,null,null,null,null,null,null],[0.00263157894736842,0.00263157894736842,0.00263157894736842,0.00263157894736842,0.00263157894736842,0.00263157894736842,0.00263157894736842,0.00263157894736842]],"container":"<table class=\"display\">\n  <thead>\n    <tr>\n      <th> <\/th>\n      <th>Axis.1<\/th>\n      <th>Axis.2<\/th>\n      <th>pval<\/th>\n      <th>r<\/th>\n      <th>Kingdom<\/th>\n      <th>Phylum<\/th>\n      <th>Class<\/th>\n      <th>Order<\/th>\n      <th>tax_rank_plot<\/th>\n      <th>Genus<\/th>\n      <th>Species<\/th>\n      <th>Strain<\/th>\n      <th>pval.adj<\/th>\n    <\/tr>\n  <\/thead>\n<\/table>","options":{"columnDefs":[{"className":"dt-right","targets":[1,2,3,4,13]},{"orderable":false,"targets":0}],"order":[],"autoWidth":false,"orderClasses":false}},"evals":[],"jsHooks":[]}</script>
```


```r
out$plot -> PCoA_tmp2_fam

PCoA_tmp2 + 
  theme(legend.position = "none") -> PCoA_tmp2_noleg

PCoA_tmp2_noleg
```

![](NRP72_draft_16S_Figs_beta_files/figure-html/unnamed-chunk-34-1.png)<!-- -->

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
<div id="htmlwidget-5d05eac75555269eefd0" style="width:100%;height:auto;" class="datatables html-widget"></div>
<script type="application/json" data-for="htmlwidget-5d05eac75555269eefd0">{"x":{"filter":"none","vertical":false,"data":[["1","2","3","4","5","6"],["Succinat_mM","Acetat_mM","Propionat_mM","Butyrat_mM","Valerat_mM","Total_SCFA_mM"],[0.805376707136237,-0.862698523914289,0.324337009107315,-0.717381936893163,-0.465080250779239,-0.90088101316165],[0.13593975442181,0.201311834400704,-0.660808679152135,0.590096643093981,-0.607635844076233,0.135182677818885],["Succinat_mM","Acetat_mM","Propionat_mM","Butyrat_mM","Valerat_mM","Total_SCFA_mM"],[0.001,0.001,0.001,0.001,0.001,0.001],[0.66711125722987,0.784775197833669,0.541862605919468,0.862850891571372,0.585520958671116,0.829860956257444]],"container":"<table class=\"display\">\n  <thead>\n    <tr>\n      <th> <\/th>\n      <th>label<\/th>\n      <th>Axis.1<\/th>\n      <th>Axis.2<\/th>\n      <th>id<\/th>\n      <th>pval<\/th>\n      <th>r<\/th>\n    <\/tr>\n  <\/thead>\n<\/table>","options":{"columnDefs":[{"className":"dt-right","targets":[2,3,5,6]},{"orderable":false,"targets":0}],"order":[],"autoWidth":false,"orderClasses":false}},"evals":[],"jsHooks":[]}</script>
```


```r
env_out$plot %>% 
  export::graph2ppt(append = TRUE,
                    width = 317.48031496 * 1,
                    height = 0.618 * 317.48031496 * 1 , paper = "A4",  scaling = 2,
                    file = out_pptx)
```





### Human2:



```r
ps_rare %>%
  subset_samples(Model2 == "Chicken2"  & Treatment != "DONOR" ) -> ps_tmp_forbdiv
```

#### NMDS:


```r
ps_tmp_forbdiv %>% 
  # microbiome::transform(transform = "compositional") %>% 
  phyloseq_plot_bdiv(dlist = dist,
                     m = "NMDS",
                     seed = 123,
                     axis1 = 1,
                     axis2 = 2) -> NMDS_tmp
```

```
## [1] "bjaccard"
## Run 0 stress 0.06501673 
## Run 1 stress 0.07028583 
## Run 2 stress 0.06698362 
## Run 3 stress 0.0713004 
## Run 4 stress 0.06546745 
## ... Procrustes: rmse 0.006146418  max resid 0.066799 
## Run 5 stress 0.06783114 
## Run 6 stress 0.06753263 
## Run 7 stress 0.07464697 
## Run 8 stress 0.06519905 
## ... Procrustes: rmse 0.02188771  max resid 0.09703256 
## Run 9 stress 0.06884219 
## Run 10 stress 0.07222409 
## Run 11 stress 0.07356041 
## Run 12 stress 0.06933559 
## Run 13 stress 0.07049735 
## Run 14 stress 0.06540318 
## ... Procrustes: rmse 0.01886117  max resid 0.09842987 
## Run 15 stress 0.06680875 
## Run 16 stress 0.06750834 
## Run 17 stress 0.06808436 
## Run 18 stress 0.07322001 
## Run 19 stress 0.07391297 
## Run 20 stress 0.0676097 
## *** Best solution was not repeated -- monoMDS stopping criteria:
##     15: stress ratio > sratmax
##      5: scale factor of the gradient < sfgrmin
## [1] "aichinson"
## Run 0 stress 0.1040082 
## Run 1 stress 0.1160556 
## Run 2 stress 0.1040098 
## ... Procrustes: rmse 0.0001945363  max resid 0.001660911 
## ... Similar to previous best
## Run 3 stress 0.1487657 
## Run 4 stress 0.1202518 
## Run 5 stress 0.1040082 
## ... Procrustes: rmse 2.525416e-05  max resid 0.0001607491 
## ... Similar to previous best
## Run 6 stress 0.1280629 
## Run 7 stress 0.1529914 
## Run 8 stress 0.1518072 
## Run 9 stress 0.1040098 
## ... Procrustes: rmse 0.0001969416  max resid 0.001681065 
## ... Similar to previous best
## Run 10 stress 0.1202518 
## Run 11 stress 0.1494985 
## Run 12 stress 0.1536904 
## Run 13 stress 0.1160606 
## Run 14 stress 0.147336 
## Run 15 stress 0.1040083 
## ... Procrustes: rmse 7.627537e-05  max resid 0.0004288902 
## ... Similar to previous best
## Run 16 stress 0.1160576 
## Run 17 stress 0.1495965 
## Run 18 stress 0.1040082 
## ... New best solution
## ... Procrustes: rmse 4.541521e-05  max resid 0.0002928945 
## ... Similar to previous best
## Run 19 stress 0.1040082 
## ... Procrustes: rmse 2.369947e-05  max resid 0.0001186972 
## ... Similar to previous best
## Run 20 stress 0.1040098 
## ... Procrustes: rmse 0.0001945314  max resid 0.001672125 
## ... Similar to previous best
## *** Best solution repeated 3 times
```

#### aichinson


```r
NMDS_tmp$aichinson$layers[[1]] = NULL; NMDS_tmp$aichinson$layers[[1]] = NULL
NMDS_tmp$aichinson$layers[[2]] = NULL; NMDS_tmp$aichinson$layers[[1]] = NULL

# plots_hall_humans$aichinson$layers[[1]] = NULL;plots_hall_humans$aichinson$layers[[1]] = NULL
# plots_hall_humans$aichinson$layers[[2]] = NULL;plots_hall_humans$aichinson$layers[[2]] = NULL


NMDS_tmp$aichinson + geom_point(size = 3,
                                aes(colour = Treatment ,
                                    shape = Antibiotic_mg.L,
                                    alpha = NULL)) +
  geom_path(data = NMDS_tmp$aichinson$data %>%
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
  facet_grid(. ~ Model2, scales = "free") -> NMDS_tmp2

NMDS_tmp2
```

![](NRP72_draft_16S_Figs_beta_files/figure-html/unnamed-chunk-39-1.png)<!-- -->


```r
NMDS_tmp2 + 
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
<div id="htmlwidget-7531f021030c44f57ef2" style="width:100%;height:auto;" class="datatables html-widget"></div>
<script type="application/json" data-for="htmlwidget-7531f021030c44f57ef2">{"x":{"filter":"none","vertical":false,"data":[["ASV0017","ASV0088","ASV0007","ASV0021","ASV0026","ASV0037","ASV0101","ASV0188"],[0.121951517572078,0.3903789246638,-0.574050010181196,0.223621804365789,0.795485555849178,-0.811759510728424,0.524826030206649,0.71786257852125],[-0.784334376671912,-0.727432738348284,0.640754980748867,-0.581167163295317,0.324605021176264,0.42486783158517,0.690527099440638,0.294577706525419],[0.001,0.001,0.001,0.001,0.001,0.001,0.001,0.001],[0.630052587067449,0.681554093642548,0.740100359543512,0.387761983080536,0.738165689337519,0.839466177573934,0.752270037044371,0.602102706822954],["Bacteria","Bacteria","Bacteria","Bacteria","Bacteria","Bacteria","Bacteria","Bacteria"],["Bacteroidota","Firmicutes","Firmicutes","Firmicutes","Firmicutes","Proteobacteria","Firmicutes","Firmicutes"],["Bacteroidia","Clostridia","Bacilli","Bacilli","Clostridia","Gammaproteobacteria","Clostridia","Clostridia"],["Bacteroidales","Oscillospirales","Lactobacillales","Bacillales","Oscillospirales","Burkholderiales","Lachnospirales","Monoglobales"],["Bacteroidaceae","Ruminococcaceae","Lactobacillaceae","Bacillaceae","Oscillospiraceae","Sutterellaceae","Defluviitaleaceae","Monoglobaceae"],[null,null,null,null,null,null,null,null],[null,null,null,null,null,null,null,null],[null,null,null,null,null,null,null,null],[0.00294117647058824,0.00294117647058824,0.00294117647058824,0.00294117647058824,0.00294117647058824,0.00294117647058824,0.00294117647058824,0.00294117647058824]],"container":"<table class=\"display\">\n  <thead>\n    <tr>\n      <th> <\/th>\n      <th>Axis.1<\/th>\n      <th>Axis.2<\/th>\n      <th>pval<\/th>\n      <th>r<\/th>\n      <th>Kingdom<\/th>\n      <th>Phylum<\/th>\n      <th>Class<\/th>\n      <th>Order<\/th>\n      <th>tax_rank_plot<\/th>\n      <th>Genus<\/th>\n      <th>Species<\/th>\n      <th>Strain<\/th>\n      <th>pval.adj<\/th>\n    <\/tr>\n  <\/thead>\n<\/table>","options":{"columnDefs":[{"className":"dt-right","targets":[1,2,3,4,13]},{"orderable":false,"targets":0}],"order":[],"autoWidth":false,"orderClasses":false}},"evals":[],"jsHooks":[]}</script>
```


```r
out$plot -> NMDS_tmp2_fam

NMDS_tmp2 + 
  theme(legend.position = "none") -> NMDS_tmp2_noleg

NMDS_tmp2_noleg
```

![](NRP72_draft_16S_Figs_beta_files/figure-html/unnamed-chunk-41-1.png)<!-- -->

```r
NMDS_tmp2  %>%  ggpubr::get_legend() %>% ggpubr::as_ggplot() -> NMDS_tmp2_leg

NMDS_tmp2_noleg %>% 
  export::graph2ppt(append = TRUE,
                    width = 317.48031496 * 1,
                    height = 0.618 * 317.48031496 * 1 , paper = "A4",  scaling = 2,
                    file = out_pptx)

NMDS_tmp2_leg %>% 
  export::graph2ppt(append = TRUE,
                    width = 317.48031496 * 1,
                    height = 0.618 * 317.48031496 * 1 , paper = "A4",  scaling = 2,
                    file = out_pptx)

NMDS_tmp2_fam %>% 
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
<div id="htmlwidget-983315589fed18d622a7" style="width:100%;height:auto;" class="datatables html-widget"></div>
<script type="application/json" data-for="htmlwidget-983315589fed18d622a7">{"x":{"filter":"none","vertical":false,"data":[["1","2","3","4","5","6"],["Succinat_mM","Acetat_mM","Propionat_mM","Butyrat_mM","Valerat_mM","Total_SCFA_mM"],[-0.529147382478684,0.81213679042208,-0.691063349217276,0.850938508204205,0.816344909441183,0.838127471001857],[0.28729710568181,-0.173455788507577,-0.265058981587657,-0.15452103634438,0.327477365938189,-0.0867683544108123],["Succinat_mM","Acetat_mM","Propionat_mM","Butyrat_mM","Valerat_mM","Total_SCFA_mM"],[0.001,0.001,0.001,0.001,0.001,0.001],[0.362536579317188,0.689653076923863,0.547824816351684,0.747973095417739,0.773660436372348,0.70998640497513]],"container":"<table class=\"display\">\n  <thead>\n    <tr>\n      <th> <\/th>\n      <th>label<\/th>\n      <th>Axis.1<\/th>\n      <th>Axis.2<\/th>\n      <th>id<\/th>\n      <th>pval<\/th>\n      <th>r<\/th>\n    <\/tr>\n  <\/thead>\n<\/table>","options":{"columnDefs":[{"className":"dt-right","targets":[2,3,5,6]},{"orderable":false,"targets":0}],"order":[],"autoWidth":false,"orderClasses":false}},"evals":[],"jsHooks":[]}</script>
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
NMDS_tmp$bjaccard$layers[[1]] = NULL; NMDS_tmp$bjaccard$layers[[1]] = NULL
NMDS_tmp$bjaccard$layers[[2]] = NULL; NMDS_tmp$bjaccard$layers[[1]] = NULL

# plots_hall_humans$aichinson$layers[[1]] = NULL;plots_hall_humans$aichinson$layers[[1]] = NULL
# plots_hall_humans$aichinson$layers[[2]] = NULL;plots_hall_humans$aichinson$layers[[2]] = NULL


NMDS_tmp$bjaccard + geom_point(size = 3,
                               aes(colour = Treatment ,
                                   shape = Antibiotic_mg.L,
                                   alpha = NULL)) +
  geom_path(data = NMDS_tmp$bjaccard$data %>%
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
  facet_grid(. ~ Model2, scales = "free") -> NMDS_tmp2

NMDS_tmp2
```

![](NRP72_draft_16S_Figs_beta_files/figure-html/unnamed-chunk-44-1.png)<!-- -->


```r
NMDS_tmp2 + 
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
<div id="htmlwidget-6f1d8b9e786d9d3fe28d" style="width:100%;height:auto;" class="datatables html-widget"></div>
<script type="application/json" data-for="htmlwidget-6f1d8b9e786d9d3fe28d">{"x":{"filter":"none","vertical":false,"data":[["ASV0017","ASV0088","ASV0007","ASV0021","ASV0026","ASV0037","ASV0101","ASV0188"],[-0.155893725671808,-0.413842315029179,0.607678896355686,-0.233359443007561,-0.623494880289104,0.861704456168324,-0.340307642417437,-0.607019222277422],[-0.651788326623972,-0.717053733649682,0.628652573232729,-0.608858337143351,0.606848572086881,0.356073867053972,0.816685999540255,0.446736709942202],[0.001,0.001,0.001,0.001,0.001,0.001,0.001,0.001],[0.449130876427114,0.685431518649659,0.764477698908196,0.425165104349765,0.75701105519061,0.869323168579117,0.782785313332779,0.56804602422427],["Bacteria","Bacteria","Bacteria","Bacteria","Bacteria","Bacteria","Bacteria","Bacteria"],["Bacteroidota","Firmicutes","Firmicutes","Firmicutes","Firmicutes","Proteobacteria","Firmicutes","Firmicutes"],["Bacteroidia","Clostridia","Bacilli","Bacilli","Clostridia","Gammaproteobacteria","Clostridia","Clostridia"],["Bacteroidales","Oscillospirales","Lactobacillales","Bacillales","Oscillospirales","Burkholderiales","Lachnospirales","Monoglobales"],["Bacteroidaceae","Ruminococcaceae","Lactobacillaceae","Bacillaceae","Oscillospiraceae","Sutterellaceae","Defluviitaleaceae","Monoglobaceae"],[null,null,null,null,null,null,null,null],[null,null,null,null,null,null,null,null],[null,null,null,null,null,null,null,null],[0.00294117647058824,0.00294117647058824,0.00294117647058824,0.00294117647058824,0.00294117647058824,0.00294117647058824,0.00294117647058824,0.00294117647058824]],"container":"<table class=\"display\">\n  <thead>\n    <tr>\n      <th> <\/th>\n      <th>Axis.1<\/th>\n      <th>Axis.2<\/th>\n      <th>pval<\/th>\n      <th>r<\/th>\n      <th>Kingdom<\/th>\n      <th>Phylum<\/th>\n      <th>Class<\/th>\n      <th>Order<\/th>\n      <th>tax_rank_plot<\/th>\n      <th>Genus<\/th>\n      <th>Species<\/th>\n      <th>Strain<\/th>\n      <th>pval.adj<\/th>\n    <\/tr>\n  <\/thead>\n<\/table>","options":{"columnDefs":[{"className":"dt-right","targets":[1,2,3,4,13]},{"orderable":false,"targets":0}],"order":[],"autoWidth":false,"orderClasses":false}},"evals":[],"jsHooks":[]}</script>
```


```r
out$plot -> NMDS_tmp2_fam

NMDS_tmp2 + 
  theme(legend.position = "none") -> NMDS_tmp2_noleg

NMDS_tmp2_noleg
```

![](NRP72_draft_16S_Figs_beta_files/figure-html/unnamed-chunk-46-1.png)<!-- -->

```r
NMDS_tmp2  %>%  ggpubr::get_legend() %>% ggpubr::as_ggplot() -> NMDS_tmp2_leg

NMDS_tmp2_noleg %>% 
  export::graph2ppt(append = TRUE,
                    width = 317.48031496 * 1,
                    height = 0.618 * 317.48031496 * 1 , paper = "A4",  scaling = 2,
                    file = out_pptx)

NMDS_tmp2_leg %>% 
  export::graph2ppt(append = TRUE,
                    width = 317.48031496 * 1,
                    height = 0.618 * 317.48031496 * 1 , paper = "A4",  scaling = 2,
                    file = out_pptx)

NMDS_tmp2_fam %>% 
  export::graph2ppt(append = TRUE,
                    width = 317.48031496 * 1,
                    height = 0.618 * 317.48031496 * 1 , paper = "A4",  scaling = 2,
                    file = out_pptx)
```


```r
ps_rare %>%
  subset_samples(Model2 == "Chicken2" & Treatment != "DONOR" ) %>%
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
<div id="htmlwidget-dcdfc64296f5794e54a4" style="width:100%;height:auto;" class="datatables html-widget"></div>
<script type="application/json" data-for="htmlwidget-dcdfc64296f5794e54a4">{"x":{"filter":"none","vertical":false,"data":[["1","2","3","4","5","6"],["Succinat_mM","Acetat_mM","Propionat_mM","Butyrat_mM","Valerat_mM","Total_SCFA_mM"],[0.585052051276282,-0.821564776006437,0.644134692028062,-0.842693140747185,-0.659182695904131,-0.807880122704603],[0.137835067181004,0.112884348508983,-0.309596044903804,0.19370665347383,0.563301066484276,0.247692902410314],["Succinat_mM","Acetat_mM","Propionat_mM","Butyrat_mM","Valerat_mM","Total_SCFA_mM"],[0.001,0.001,0.001,0.001,0.001,0.001],[0.361284408447377,0.687711557312805,0.510759212494164,0.747653997062385,0.751829918081761,0.714022066565649]],"container":"<table class=\"display\">\n  <thead>\n    <tr>\n      <th> <\/th>\n      <th>label<\/th>\n      <th>Axis.1<\/th>\n      <th>Axis.2<\/th>\n      <th>id<\/th>\n      <th>pval<\/th>\n      <th>r<\/th>\n    <\/tr>\n  <\/thead>\n<\/table>","options":{"columnDefs":[{"className":"dt-right","targets":[2,3,5,6]},{"orderable":false,"targets":0}],"order":[],"autoWidth":false,"orderClasses":false}},"evals":[],"jsHooks":[]}</script>
```


```r
env_out$plot %>% 
  export::graph2ppt(append = TRUE,
                    width = 317.48031496 * 1,
                    height = 0.618 * 317.48031496 * 1 , paper = "A4",  scaling = 2,
                    file = out_pptx)
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

![](NRP72_draft_16S_Figs_beta_files/figure-html/unnamed-chunk-50-1.png)<!-- -->


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
<div id="htmlwidget-873249b169eb6149d57e" style="width:100%;height:auto;" class="datatables html-widget"></div>
<script type="application/json" data-for="htmlwidget-873249b169eb6149d57e">{"x":{"filter":"none","vertical":false,"data":[["ASV0017","ASV0088","ASV0007","ASV0021","ASV0026","ASV0037","ASV0101","ASV0188"],[0.121951517572078,0.3903789246638,-0.574050010181196,0.223621804365789,0.795485555849178,-0.811759510728424,0.524826030206649,0.71786257852125],[-0.784334376671912,-0.727432738348284,0.640754980748867,-0.581167163295317,0.324605021176264,0.42486783158517,0.690527099440638,0.294577706525419],[0.001,0.001,0.001,0.001,0.001,0.001,0.001,0.001],[0.630052587067449,0.681554093642548,0.740100359543512,0.387761983080536,0.738165689337519,0.839466177573934,0.752270037044371,0.602102706822954],["Bacteria","Bacteria","Bacteria","Bacteria","Bacteria","Bacteria","Bacteria","Bacteria"],["Bacteroidota","Firmicutes","Firmicutes","Firmicutes","Firmicutes","Proteobacteria","Firmicutes","Firmicutes"],["Bacteroidia","Clostridia","Bacilli","Bacilli","Clostridia","Gammaproteobacteria","Clostridia","Clostridia"],["Bacteroidales","Oscillospirales","Lactobacillales","Bacillales","Oscillospirales","Burkholderiales","Lachnospirales","Monoglobales"],["Bacteroidaceae","Ruminococcaceae","Lactobacillaceae","Bacillaceae","Oscillospiraceae","Sutterellaceae","Defluviitaleaceae","Monoglobaceae"],[null,null,null,null,null,null,null,null],[null,null,null,null,null,null,null,null],[null,null,null,null,null,null,null,null],[0.00294117647058824,0.00294117647058824,0.00294117647058824,0.00294117647058824,0.00294117647058824,0.00294117647058824,0.00294117647058824,0.00294117647058824]],"container":"<table class=\"display\">\n  <thead>\n    <tr>\n      <th> <\/th>\n      <th>Axis.1<\/th>\n      <th>Axis.2<\/th>\n      <th>pval<\/th>\n      <th>r<\/th>\n      <th>Kingdom<\/th>\n      <th>Phylum<\/th>\n      <th>Class<\/th>\n      <th>Order<\/th>\n      <th>tax_rank_plot<\/th>\n      <th>Genus<\/th>\n      <th>Species<\/th>\n      <th>Strain<\/th>\n      <th>pval.adj<\/th>\n    <\/tr>\n  <\/thead>\n<\/table>","options":{"columnDefs":[{"className":"dt-right","targets":[1,2,3,4,13]},{"orderable":false,"targets":0}],"order":[],"autoWidth":false,"orderClasses":false}},"evals":[],"jsHooks":[]}</script>
```


```r
out$plot -> PCoA_tmp2_fam

PCoA_tmp2 + 
  theme(legend.position = "none") -> PCoA_tmp2_noleg

PCoA_tmp2_noleg
```

![](NRP72_draft_16S_Figs_beta_files/figure-html/unnamed-chunk-52-1.png)<!-- -->

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
<div id="htmlwidget-cbd3ad10bd2d35a2e76d" style="width:100%;height:auto;" class="datatables html-widget"></div>
<script type="application/json" data-for="htmlwidget-cbd3ad10bd2d35a2e76d">{"x":{"filter":"none","vertical":false,"data":[["1","2","3","4","5","6"],["Succinat_mM","Acetat_mM","Propionat_mM","Butyrat_mM","Valerat_mM","Total_SCFA_mM"],[-0.529147382478684,0.81213679042208,-0.691063349217276,0.850938508204205,0.816344909441183,0.838127471001857],[0.28729710568181,-0.173455788507577,-0.265058981587657,-0.15452103634438,0.327477365938189,-0.0867683544108123],["Succinat_mM","Acetat_mM","Propionat_mM","Butyrat_mM","Valerat_mM","Total_SCFA_mM"],[0.001,0.001,0.001,0.001,0.001,0.001],[0.362536579317188,0.689653076923863,0.547824816351684,0.747973095417739,0.773660436372348,0.70998640497513]],"container":"<table class=\"display\">\n  <thead>\n    <tr>\n      <th> <\/th>\n      <th>label<\/th>\n      <th>Axis.1<\/th>\n      <th>Axis.2<\/th>\n      <th>id<\/th>\n      <th>pval<\/th>\n      <th>r<\/th>\n    <\/tr>\n  <\/thead>\n<\/table>","options":{"columnDefs":[{"className":"dt-right","targets":[2,3,5,6]},{"orderable":false,"targets":0}],"order":[],"autoWidth":false,"orderClasses":false}},"evals":[],"jsHooks":[]}</script>
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

![](NRP72_draft_16S_Figs_beta_files/figure-html/unnamed-chunk-55-1.png)<!-- -->


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
<div id="htmlwidget-f25d4e1770d66a86dcce" style="width:100%;height:auto;" class="datatables html-widget"></div>
<script type="application/json" data-for="htmlwidget-f25d4e1770d66a86dcce">{"x":{"filter":"none","vertical":false,"data":[["ASV0017","ASV0088","ASV0007","ASV0021","ASV0026","ASV0037","ASV0101","ASV0188"],[-0.155893725671808,-0.413842315029179,0.607678896355686,-0.233359443007561,-0.623494880289104,0.861704456168324,-0.340307642417437,-0.607019222277422],[-0.651788326623972,-0.717053733649682,0.628652573232729,-0.608858337143351,0.606848572086881,0.356073867053972,0.816685999540255,0.446736709942202],[0.001,0.001,0.001,0.001,0.001,0.001,0.001,0.001],[0.449130876427114,0.685431518649659,0.764477698908196,0.425165104349765,0.75701105519061,0.869323168579117,0.782785313332779,0.56804602422427],["Bacteria","Bacteria","Bacteria","Bacteria","Bacteria","Bacteria","Bacteria","Bacteria"],["Bacteroidota","Firmicutes","Firmicutes","Firmicutes","Firmicutes","Proteobacteria","Firmicutes","Firmicutes"],["Bacteroidia","Clostridia","Bacilli","Bacilli","Clostridia","Gammaproteobacteria","Clostridia","Clostridia"],["Bacteroidales","Oscillospirales","Lactobacillales","Bacillales","Oscillospirales","Burkholderiales","Lachnospirales","Monoglobales"],["Bacteroidaceae","Ruminococcaceae","Lactobacillaceae","Bacillaceae","Oscillospiraceae","Sutterellaceae","Defluviitaleaceae","Monoglobaceae"],[null,null,null,null,null,null,null,null],[null,null,null,null,null,null,null,null],[null,null,null,null,null,null,null,null],[0.00294117647058824,0.00294117647058824,0.00294117647058824,0.00294117647058824,0.00294117647058824,0.00294117647058824,0.00294117647058824,0.00294117647058824]],"container":"<table class=\"display\">\n  <thead>\n    <tr>\n      <th> <\/th>\n      <th>Axis.1<\/th>\n      <th>Axis.2<\/th>\n      <th>pval<\/th>\n      <th>r<\/th>\n      <th>Kingdom<\/th>\n      <th>Phylum<\/th>\n      <th>Class<\/th>\n      <th>Order<\/th>\n      <th>tax_rank_plot<\/th>\n      <th>Genus<\/th>\n      <th>Species<\/th>\n      <th>Strain<\/th>\n      <th>pval.adj<\/th>\n    <\/tr>\n  <\/thead>\n<\/table>","options":{"columnDefs":[{"className":"dt-right","targets":[1,2,3,4,13]},{"orderable":false,"targets":0}],"order":[],"autoWidth":false,"orderClasses":false}},"evals":[],"jsHooks":[]}</script>
```


```r
out$plot -> PCoA_tmp2_fam

PCoA_tmp2 + 
  theme(legend.position = "none") -> PCoA_tmp2_noleg

PCoA_tmp2_noleg
```

![](NRP72_draft_16S_Figs_beta_files/figure-html/unnamed-chunk-57-1.png)<!-- -->

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
<div id="htmlwidget-969965edfe9f869aa904" style="width:100%;height:auto;" class="datatables html-widget"></div>
<script type="application/json" data-for="htmlwidget-969965edfe9f869aa904">{"x":{"filter":"none","vertical":false,"data":[["1","2","3","4","5","6"],["Succinat_mM","Acetat_mM","Propionat_mM","Butyrat_mM","Valerat_mM","Total_SCFA_mM"],[0.585052051276282,-0.821564776006437,0.644134692028062,-0.842693140747185,-0.659182695904131,-0.807880122704603],[0.137835067181004,0.112884348508983,-0.309596044903804,0.19370665347383,0.563301066484276,0.247692902410314],["Succinat_mM","Acetat_mM","Propionat_mM","Butyrat_mM","Valerat_mM","Total_SCFA_mM"],[0.001,0.001,0.001,0.001,0.001,0.001],[0.361284408447377,0.687711557312805,0.510759212494164,0.747653997062385,0.751829918081761,0.714022066565649]],"container":"<table class=\"display\">\n  <thead>\n    <tr>\n      <th> <\/th>\n      <th>label<\/th>\n      <th>Axis.1<\/th>\n      <th>Axis.2<\/th>\n      <th>id<\/th>\n      <th>pval<\/th>\n      <th>r<\/th>\n    <\/tr>\n  <\/thead>\n<\/table>","options":{"columnDefs":[{"className":"dt-right","targets":[2,3,5,6]},{"orderable":false,"targets":0}],"order":[],"autoWidth":false,"orderClasses":false}},"evals":[],"jsHooks":[]}</script>
```


```r
env_out$plot %>% 
  export::graph2ppt(append = TRUE,
                    width = 317.48031496 * 1,
                    height = 0.618 * 317.48031496 * 1 , paper = "A4",  scaling = 2,
                    file = out_pptx)
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
##  [1] vegan_2.6-4          lattice_0.20-45      permute_0.9-7       
##  [4] gdtools_0.2.4        GUniFrac_1.6         ape_5.6-2           
##  [7] speedyseq_0.5.3.9018 reshape2_1.4.4       scales_1.2.1        
## [10] compositions_2.0-4   nlme_3.1-159         phyloseq_1.40.0     
## [13] forcats_0.5.2        stringr_1.4.1        dplyr_1.0.10        
## [16] purrr_0.3.5          readr_2.1.3          tidyr_1.2.1         
## [19] tibble_3.1.8         ggplot2_3.3.6        tidyverse_1.3.2     
## 
## loaded via a namespace (and not attached):
##   [1] uuid_1.1-0             readxl_1.4.1           backports_1.4.1       
##   [4] systemfonts_1.0.4      plyr_1.8.7             igraph_1.3.5          
##   [7] splines_4.2.1          crosstalk_1.2.0        GenomeInfoDb_1.32.4   
##  [10] digest_0.6.30          foreach_1.5.2          htmltools_0.5.3       
##  [13] fansi_1.0.3            magrittr_2.0.3         googlesheets4_1.0.1   
##  [16] cluster_2.1.4          openxlsx_4.2.5         tzdb_0.3.0            
##  [19] Biostrings_2.64.1      extrafont_0.18         modelr_0.1.9          
##  [22] bayesm_3.1-4           matrixStats_0.62.0     officer_0.4.4         
##  [25] stabledist_0.7-1       extrafontdb_1.0        colorspace_2.0-3      
##  [28] rvest_1.0.3            ggrepel_0.9.1          textshaping_0.3.6     
##  [31] haven_2.5.1            xfun_0.34              crayon_1.5.2          
##  [34] RCurl_1.98-1.9         jsonlite_1.8.3         survival_3.4-0        
##  [37] iterators_1.0.14       glue_1.6.2             rvg_0.2.5             
##  [40] gtable_0.3.1           gargle_1.2.1           zlibbioc_1.42.0       
##  [43] XVector_0.36.0         car_3.1-0              Rttf2pt1_1.3.10       
##  [46] Rhdf5lib_1.18.2        BiocGenerics_0.42.0    DEoptimR_1.0-11       
##  [49] abind_1.4-5            DBI_1.1.3              rstatix_0.7.0         
##  [52] Rcpp_1.0.9             xtable_1.8-4           clue_0.3-62           
##  [55] DT_0.25                stats4_4.2.1           htmlwidgets_1.5.4     
##  [58] timeSeries_4021.104    httr_1.4.4             ellipsis_0.3.2        
##  [61] spatial_7.3-15         pkgconfig_2.0.3        farver_2.1.1          
##  [64] sass_0.4.2             dbplyr_2.2.1           utf8_1.2.2            
##  [67] here_1.0.1             tidyselect_1.2.0       labeling_0.4.2        
##  [70] rlang_1.0.6            munsell_0.5.0          cellranger_1.1.0      
##  [73] tools_4.2.1            cachem_1.0.6           cli_3.4.1             
##  [76] generics_0.1.3         devEMF_4.1             ade4_1.7-19           
##  [79] export_0.3.0           broom_1.0.1            evaluate_0.17         
##  [82] biomformat_1.24.0      fastmap_1.1.0          yaml_2.3.6            
##  [85] ragg_1.2.3             knitr_1.40             fs_1.5.2              
##  [88] zip_2.2.1              robustbase_0.95-0      rgl_0.110.2           
##  [91] xml2_1.3.3             compiler_4.2.1         rstudioapi_0.14       
##  [94] ggsignif_0.6.4         reprex_2.0.2           statmod_1.4.37        
##  [97] bslib_0.4.0            stringi_1.7.8          statip_0.2.3          
## [100] highr_0.9              stargazer_5.2.3        modeest_2.4.0         
## [103] fBasics_4021.92        Matrix_1.5-1           microbiome_1.19.1     
## [106] tensorA_0.36.2         multtest_2.52.0        vctrs_0.4.2           
## [109] pillar_1.8.1           lifecycle_1.0.3        rhdf5filters_1.8.0    
## [112] microViz_0.9.7         jquerylib_0.1.4        flextable_0.8.2       
## [115] cowplot_1.1.1          data.table_1.14.4      bitops_1.0-7          
## [118] R6_2.5.1               stable_1.1.6           IRanges_2.30.1        
## [121] codetools_0.2-18       MASS_7.3-58.1          assertthat_0.2.1      
## [124] rhdf5_2.40.0           rprojroot_2.0.3        withr_2.5.0           
## [127] S4Vectors_0.34.0       GenomeInfoDbData_1.2.8 mgcv_1.8-40           
## [130] parallel_4.2.1         hms_1.1.2              grid_4.2.1            
## [133] rpart_4.1.16           timeDate_4021.106      rmarkdown_2.16        
## [136] carData_3.0-5          rmutil_1.1.9           googledrive_2.0.0     
## [139] Rtsne_0.16             ggpubr_0.4.0           base64enc_0.1-3       
## [142] Biobase_2.56.0         lubridate_1.8.0
```
