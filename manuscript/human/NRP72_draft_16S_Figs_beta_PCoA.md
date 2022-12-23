---
title: "NTP72 - 16S - beta - human draft "
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
out_pptx = "~/Desktop/16S-human-beta.pptx"
out_xlsx = "~/Desktop/16S-human-beta.xlsx"
```


```r
load(url("https://github.com/fconstancias/NRP72-FBT/blob/master/Figures/Rastetics.Rdata?raw=true"))
```

## load data:


```r
"data/processed/16S/16S_working_phyloseq.RDS" %>% 
  here::here() %>% 
  readRDS() %>% 
  subset_samples(Experiment %in% c("Continuous", "Cecum")) %>%
  subset_samples(Day_of_Treatment > -6 & Day_of_Treatment <= 30 | Reactor == "DONOR") %>% 
  subset_samples(Model == "Human") %>% 
  microViz::ps_mutate(period = case_when(Day_of_Treatment <= 0 ~ "pret",
                                         Day_of_Treatment > 0 &  Day_of_Treatment < 3 ~ "t1" ,
                                         Day_of_Treatment >= 3 & Day_of_Treatment <= 6 ~ "t2",
                                         Day_of_Treatment > 6  ~ "t3" )) -> ps_filtered
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
  select(SampleID, Day_of_Treatment,Treatment, Reactor, Model2, Model, LibrarySize) 
```

```
##        SampleID Day_of_Treatment     Treatment Reactor Model2 Model LibrarySize
## 111  TR5-68-S99                6       HV292.1     TR5 Human1 Human           0
## 18  TR1-63-S100                1   CTX+HV292.1     TR1 Human1 Human           5
## 86  TR4-65-S155                3           CTX     TR4 Human1 Human           5
## 108 TR5-65-S187                3       HV292.1     TR5 Human1 Human           5
## 123 TR6-36-S124                7     UNTREATED      CR Human1 Human           5
## 19  TR1-64-S120                2   CTX+HV292.1     TR1 Human1 Human           6
## 34  TR2-35-S131                6           VAN     TR2 Human1 Human           6
## 107 TR5-64-S108                2       HV292.1     TR5 Human1 Human           7
## 125 TR6-60-S119               -2     UNTREATED      CR Human1 Human           8
## 28  TR2-29-S117                0           VAN     TR2 Human1 Human          10
## 80  TR4-59-S180               -3           CTX     TR4 Human1 Human          11
## 9   TR1-32-S111                3 VAN+CCUG59168     TR1 Human1 Human          15
## 109 TR5-66-S101                4       HV292.1     TR5 Human1 Human          17
## 112 TR5-69-S107                7       HV292.1     TR5 Human1 Human          18
## 8   TR1-31-S109                2 VAN+CCUG59168     TR1 Human1 Human          20
## 88  TR4-67-S192                5           CTX     TR4 Human1 Human          22
## 20  TR1-65-S128                3   CTX+HV292.1     TR1 Human1 Human          23
## 133 TR6-68-S127                6     UNTREATED      CR Human1 Human          25
## 67  TR3-68-S194                6   CTX+HV292.1     TR3 Human1 Human          31
## 22  TR1-67-S133                5   CTX+HV292.1     TR1 Human1 Human          37
## 55  TR3-34-S149                5 VAN+CCUG59168     TR3 Human1 Human          37
## 15  TR1-60-S114               -2   CTX+HV292.1     TR1 Human1 Human          47
## 39   TR2-62-S97                0           CTX     TR2 Human1 Human          67
## 71  TR4-28-S138               -1           VAN     TR4 Human1 Human          70
## 129 TR6-64-S105                2     UNTREATED      CR Human1 Human          89
## 99   TR5-34-S28                5     CCUG59168     TR5 Human1 Human         110
## 91   TR5-26-S24               -3     CCUG59168     TR5 Human1 Human         124
## 81   TR4-60-S36               -2           CTX     TR4 Human1 Human         139
## 1       D-1-S49             <NA>         DONOR   DONOR Human1 Human         145
## 104 TR5-61-S142               -1       HV292.1     TR5 Human1 Human         161
## 70  TR4-27-S130               -2           VAN     TR4 Human1 Human         170
## 113 TR6-26-S122               -3     UNTREATED      CR Human1 Human         174
## 61   TR3-62-S37                0   CTX+HV292.1     TR3 Human1 Human         178
## 25   TR2-26-S87               -3           VAN     TR2 Human1 Human         269
## 56   TR3-35-S20                6 VAN+CCUG59168     TR3 Human1 Human         624
## 78  TR4-35-S202                6           VAN     TR4 Human1 Human        3827
## 209   AP-67-S67                3           VAN     TR2 Human2 Human        4028
## 172   AP-29-S29               -2 VAN+CCUG59168     TR4 Human2 Human        5775
## 210   AP-68-S68                3           CTX     TR3 Human2 Human        6133
## 83    TR4-62-S1                0           CTX     TR4 Human1 Human        8648
## 57    TR3-36-S8                7 VAN+CCUG59168     TR3 Human1 Human        9148
## 189   AP-46-S46                0   CTX+HV292.1     TR5 Human2 Human        9243
## 211   AP-69-S69                3 VAN+CCUG59168     TR4 Human2 Human        9752
## 239  S-238-S238             <NA>         DONOR   DONOR Human2 Human       10219
## 227   AP-86-S86                5   CTX+HV292.1     TR5 Human2 Human       10408
## 178   AP-35-S35               -1           VAN     TR2 Human2 Human       10576
## 170   AP-27-S27               -2           VAN     TR2 Human2 Human       10699
## 171   AP-28-S28               -2           CTX     TR3 Human2 Human       11109
## 114   TR6-27-S6               -2     UNTREATED      CR Human1 Human       11405
## 201   AP-59-S59                2           VAN     TR2 Human2 Human       11560
## 38  TR2-61-S196               -1           CTX     TR2 Human1 Human       11579
## 213   AP-70-S70                3   CTX+HV292.1     TR5 Human2 Human       11615
## 195   AP-52-S52                1           CTX     TR3 Human2 Human       12031
## 179   AP-36-S36               -1           CTX     TR3 Human2 Human       12128
## 110  TR5-67-S96                5       HV292.1     TR5 Human1 Human       12418
## 218   AP-76-S76                4           CTX     TR3 Human2 Human       12895
## 73   TR4-30-S48                1           VAN     TR4 Human1 Human       13390
## 137 AP-108-S114                8           CTX     TR3 Human2 Human       14857
## 151 AP-121-S127               10     UNTREATED      CR Human2 Human       15011
## 187   AP-44-S44                0           CTX     TR3 Human2 Human       15970
## 116 TR6-29-S207                0     UNTREATED      CR Human1 Human       16180
## 90  TR4-69-S197                7           CTX     TR4 Human1 Human       16297
## 59   TR3-60-S46               -2   CTX+HV292.1     TR3 Human1 Human       16345
## 27   TR2-28-S64               -1           VAN     TR2 Human1 Human       16380
## 97   TR5-32-S56                3     CCUG59168     TR5 Human1 Human       16708
## 233   AP-92-S92                6           CTX     TR3 Human2 Human       17210
## 197   AP-54-S54                1   CTX+HV292.1     TR5 Human2 Human       17258
## 130  TR6-65-S94                3     UNTREATED      CR Human1 Human       17289
## 122 TR6-35-S208                6     UNTREATED      CR Human1 Human       17362
## 175   AP-31-S31               -2     CCUG59168     TR6 Human2 Human       17619
## 194   AP-51-S51                1           VAN     TR2 Human2 Human       17745
## 225   AP-84-S84                5           CTX     TR3 Human2 Human       17759
## 64   TR3-65-S92                3   CTX+HV292.1     TR3 Human1 Human       18016
## 12  TR1-35-S161                6 VAN+CCUG59168     TR1 Human1 Human       18233
## 115  TR6-28-S80               -1     UNTREATED      CR Human1 Human       18317
## 183     AP-4-S4                7           CTX     TR3 Human2 Human       18425
## 154 AP-125-S131               10 VAN+CCUG59168     TR4 Human2 Human       18759
## 77   TR4-34-S16                5           VAN     TR4 Human1 Human       19080
## 165   AP-21-S21               -3 VAN+CCUG59168     TR4 Human2 Human       19152
## 164   AP-20-S20               -3           CTX     TR3 Human2 Human       19548
## 53  TR3-32-S211                3 VAN+CCUG59168     TR3 Human1 Human       19882
## 93    TR5-28-S5               -1     CCUG59168     TR5 Human1 Human       19912
## 43   TR2-66-S93                4           CTX     TR2 Human1 Human       20303
## 204   AP-61-S61                2 VAN+CCUG59168     TR4 Human2 Human       20344
## 149   AP-12-S12               -4           CTX     TR3 Human2 Human       20875
## 82   TR4-61-S11               -1           CTX     TR4 Human1 Human       20897
## 44  TR2-67-S193                5           CTX     TR2 Human1 Human       20963
## 203   AP-60-S60                2           CTX     TR3 Human2 Human       20975
## 101  TR5-36-S85                7     CCUG59168     TR5 Human1 Human       21015
## 163   AP-19-S19               -3           VAN     TR2 Human2 Human       21046
## 62  TR3-63-S206                1   CTX+HV292.1     TR3 Human1 Human       21195
## 87   TR4-66-S55                4           CTX     TR4 Human1 Human       21429
## 3    TR1-26-S27               -3 VAN+CCUG59168     TR1 Human1 Human       21455
## 192   AP-49-S49                1     UNTREATED      CR Human2 Human       21830
## 66  TR3-67-S199                5   CTX+HV292.1     TR3 Human1 Human       22085
## 182   AP-39-S39               -1     CCUG59168     TR6 Human2 Human       22197
## 119  TR6-32-S39                3     UNTREATED      CR Human1 Human       22274
## 145 AP-116-S122                9           CTX     TR3 Human2 Human       22524
## 13  TR1-36-S198                7 VAN+CCUG59168     TR1 Human1 Human       22648
## 138 AP-109-S115                8 VAN+CCUG59168     TR4 Human2 Human       22654
## 177   AP-33-S33               -1     UNTREATED      CR Human2 Human       22676
## 95   TR5-30-S62                1     CCUG59168     TR5 Human1 Human       22748
## 205   AP-62-S62                2   CTX+HV292.1     TR5 Human2 Human       22804
## 117 TR6-30-S144                1     UNTREATED      CR Human1 Human       23081
## 173     AP-3-S3                7           VAN     TR2 Human2 Human       23243
## 153 AP-124-S130               10           CTX     TR3 Human2 Human       23296
## 30   TR2-31-S32                2           VAN     TR2 Human1 Human       23361
## 58   TR3-59-S61               -3   CTX+HV292.1     TR3 Human1 Human       23409
## 85   TR4-64-S10                2           CTX     TR4 Human1 Human       23455
## 214   AP-71-S71                3     CCUG59168     TR6 Human2 Human       23560
## 17    TR1-62-S7                0   CTX+HV292.1     TR1 Human1 Human       23651
## 79   TR4-36-S12                7           VAN     TR4 Human1 Human       23673
## 140 AP-110-S116                8   CTX+HV292.1     TR5 Human2 Human       23730
## 23   TR1-68-S13                6   CTX+HV292.1     TR1 Human1 Human       23788
## 236   AP-95-S95                6     CCUG59168     TR6 Human2 Human       23969
## 226   AP-85-S85                5 VAN+CCUG59168     TR4 Human2 Human       24068
## 48  TR3-27-S203               -2 VAN+CCUG59168     TR3 Human1 Human       24296
## 158   AP-13-S13               -4 VAN+CCUG59168     TR4 Human2 Human       24342
## 132 TR6-67-S140                5     UNTREATED      CR Human1 Human       24361
## 156 AP-127-S133               10     CCUG59168     TR6 Human2 Human       24501
## 68  TR3-69-S139                7   CTX+HV292.1     TR3 Human1 Human       24568
## 33   TR2-34-S71                5           VAN     TR2 Human1 Human       24666
## 180   AP-37-S37               -1 VAN+CCUG59168     TR4 Human2 Human       24742
## 127  TR6-62-S50                0     UNTREATED      CR Human1 Human       24858
## 198   AP-55-S55                1     CCUG59168     TR6 Human2 Human       25057
## 128 TR6-63-S195                1     UNTREATED      CR Human1 Human       25165
## 76   TR4-33-S84                4           VAN     TR4 Human1 Human       25314
## 206   AP-63-S63                2     CCUG59168     TR6 Human2 Human       25332
## 36   TR2-59-S60               -3           CTX     TR2 Human1 Human       25411
## 174   AP-30-S30               -2   CTX+HV292.1     TR5 Human2 Human       25526
## 31   TR2-32-S90                3           VAN     TR2 Human1 Human       25566
## 24  TR1-69-S145                7   CTX+HV292.1     TR1 Human1 Human       26003
## 5    TR1-28-S35               -1 VAN+CCUG59168     TR1 Human1 Human       26178
## 69   TR4-26-S26               -3           VAN     TR4 Human1 Human       26281
## 26   TR2-27-S59               -2           VAN     TR2 Human1 Human       26447
## 136 AP-107-S113                8           VAN     TR2 Human2 Human       26458
## 102  TR5-59-S31               -3       HV292.1     TR5 Human1 Human       26520
## 220   AP-79-S79                4     CCUG59168     TR6 Human2 Human       26581
## 2      D-2-S164             <NA>         DONOR   DONOR Human1 Human       26814
## 51  TR3-30-S190                1 VAN+CCUG59168     TR3 Human1 Human       26905
## 29   TR2-30-S73                1           VAN     TR2 Human1 Human       26990
## 11   TR1-34-S34                5 VAN+CCUG59168     TR1 Human1 Human       27142
## 65   TR3-66-S69                4   CTX+HV292.1     TR3 Human1 Human       27549
## 219   AP-77-S77                4 VAN+CCUG59168     TR4 Human2 Human       27645
## 159   AP-14-S14               -4   CTX+HV292.1     TR5 Human2 Human       27665
## 92  TR5-27-S148               -2     CCUG59168     TR5 Human1 Human       27793
## 100  TR5-35-S77                6     CCUG59168     TR5 Human1 Human       28244
## 217   AP-75-S75                4           VAN     TR2 Human2 Human       28283
## 103  TR5-60-S75               -2       HV292.1     TR5 Human1 Human       28392
## 188   AP-45-S45                0 VAN+CCUG59168     TR4 Human2 Human       28460
## 14   TR1-59-S42               -3   CTX+HV292.1     TR1 Human1 Human       28844
## 230   AP-89-S89                6     UNTREATED      CR Human2 Human       28876
## 7    TR1-30-S76                1 VAN+CCUG59168     TR1 Human1 Human       29125
## 37   TR2-60-S79               -2           CTX     TR2 Human1 Human       29653
## 167   AP-23-S23               -3     CCUG59168     TR6 Human2 Human       29775
## 160   AP-15-S15               -4     CCUG59168     TR6 Human2 Human       29981
## 216   AP-73-S73                4     UNTREATED      CR Human2 Human       30040
## 139   AP-11-S11               -4           VAN     TR2 Human2 Human       30242
## 40  TR2-63-S160                1           CTX     TR2 Human1 Human       30256
## 141 AP-111-S117                8     CCUG59168     TR6 Human2 Human       30256
## 52   TR3-31-S23                2 VAN+CCUG59168     TR3 Human1 Human       30485
## 4    TR1-27-S54               -2 VAN+CCUG59168     TR1 Human1 Human       30510
## 166   AP-22-S22               -3   CTX+HV292.1     TR5 Human2 Human       30610
## 120 TR6-33-S135                4     UNTREATED      CR Human1 Human       30629
## 162   AP-17-S17               -3     UNTREATED      CR Human2 Human       30921
## 222   AP-80-S80                4       HV292.1     TR7 Human2 Human       31123
## 131 TR6-66-S132                4     UNTREATED      CR Human1 Human       31326
## 202     AP-6-S6                7   CTX+HV292.1     TR5 Human2 Human       31365
## 45  TR2-68-S147                6           CTX     TR2 Human1 Human       31718
## 235   AP-94-S94                6   CTX+HV292.1     TR5 Human2 Human       31823
## 185   AP-41-S41                0     UNTREATED      CR Human2 Human       32069
## 186   AP-43-S43                0           VAN     TR2 Human2 Human       32421
## 234   AP-93-S93                6 VAN+CCUG59168     TR4 Human2 Human       32513
## 224   AP-83-S83                5           VAN     TR2 Human2 Human       32887
## 106  TR5-63-S82                1       HV292.1     TR5 Human1 Human       32998
## 72  TR4-29-S159                0           VAN     TR4 Human1 Human       33009
## 94  TR5-29-S136                0     CCUG59168     TR5 Human1 Human       33131
## 50  TR3-29-S210                0 VAN+CCUG59168     TR3 Human1 Human       33494
## 228   AP-87-S87                5     CCUG59168     TR6 Human2 Human       33542
## 84   TR4-63-S74                1           CTX     TR4 Human1 Human       33606
## 74   TR4-31-S19                2           VAN     TR4 Human1 Human       33662
## 121 TR6-34-S141                5     UNTREATED      CR Human1 Human       34001
## 42  TR2-65-S158                3           CTX     TR2 Human1 Human       34103
## 181   AP-38-S38               -1   CTX+HV292.1     TR5 Human2 Human       34341
## 35  TR2-36-S183                7           VAN     TR2 Human1 Human       34383
## 6    TR1-29-S18                0 VAN+CCUG59168     TR1 Human1 Human       34515
## 200   AP-57-S57                2     UNTREATED      CR Human2 Human       34717
## 75   TR4-32-S17                3           VAN     TR4 Human1 Human       34738
## 229   AP-88-S88                5       HV292.1     TR7 Human2 Human       35163
## 231     AP-9-S9               -4     UNTREATED      CR Human2 Human       35291
## 196   AP-53-S53                1 VAN+CCUG59168     TR4 Human2 Human       36140
## 98  TR5-33-S134                4     CCUG59168     TR5 Human1 Human       36143
## 134 TR6-69-S174                7     UNTREATED      CR Human1 Human       36636
## 152 AP-123-S129               10           VAN     TR2 Human2 Human       36733
## 199   AP-56-S56                1       HV292.1     TR7 Human2 Human       36751
## 232   AP-91-S91                6           VAN     TR2 Human2 Human       37182
## 146 AP-117-S123                9 VAN+CCUG59168     TR4 Human2 Human       37424
## 193     AP-5-S5                7 VAN+CCUG59168     TR4 Human2 Human       37744
## 147 AP-118-S124                9   CTX+HV292.1     TR5 Human2 Human       37886
## 21  TR1-66-S151                4   CTX+HV292.1     TR1 Human1 Human       38035
## 60   TR3-61-S45               -1   CTX+HV292.1     TR3 Human1 Human       38547
## 105 TR5-62-S176                0       HV292.1     TR5 Human1 Human       38602
## 144 AP-115-S121                9           VAN     TR2 Human2 Human       38927
## 148 AP-119-S125                9     CCUG59168     TR6 Human2 Human       39098
## 142 AP-112-S118                8       HV292.1     TR7 Human2 Human       40362
## 191   AP-48-S48                0       HV292.1     TR7 Human2 Human       40563
## 47  TR3-26-S175               -3 VAN+CCUG59168     TR3 Human1 Human       40608
## 155 AP-126-S132               10   CTX+HV292.1     TR5 Human2 Human       41114
## 63  TR3-64-S146                2   CTX+HV292.1     TR3 Human1 Human       41383
## 118 TR6-31-S179                2     UNTREATED      CR Human1 Human       41421
## 223   AP-81-S81                5     UNTREATED      CR Human2 Human       42160
## 150 AP-120-S126                9       HV292.1     TR7 Human2 Human       42884
## 208   AP-65-S65                3     UNTREATED      CR Human2 Human       43429
## 237   AP-96-S96                6       HV292.1     TR7 Human2 Human       43781
## 207   AP-64-S64                2       HV292.1     TR7 Human2 Human       44004
## 161   AP-16-S16               -4       HV292.1     TR7 Human2 Human       44159
## 169   AP-25-S25               -2     UNTREATED      CR Human2 Human       45647
## 168   AP-24-S24               -3       HV292.1     TR7 Human2 Human       45911
## 184   AP-40-S40               -1       HV292.1     TR7 Human2 Human       47414
## 215   AP-72-S72                3       HV292.1     TR7 Human2 Human       48252
## 190   AP-47-S47                0     CCUG59168     TR6 Human2 Human       50228
## 143 AP-113-S119                9     UNTREATED      CR Human2 Human       51419
## 176   AP-32-S32               -2       HV292.1     TR7 Human2 Human       51683
## 157 AP-128-S134               10       HV292.1     TR7 Human2 Human       52805
## 135 AP-105-S111                8     UNTREATED      CR Human2 Human       53324
## 54   TR3-33-S86                4 VAN+CCUG59168     TR3 Human1 Human       54328
## 212     AP-7-S7                7     CCUG59168     TR6 Human2 Human       54839
## 221     AP-8-S8                7       HV292.1     TR7 Human2 Human       56498
## 238  AP-97-S103                7     UNTREATED      CR Human2 Human       72433
## 41   TR2-64-S51                2           CTX     TR2 Human1 Human       88037
## 89   TR4-68-S52                6           CTX     TR4 Human1 Human      124468
## 46   TR2-69-S57                7           CTX     TR2 Human1 Human      150837
## 96   TR5-31-S88                2     CCUG59168     TR5 Human1 Human      184708
## 16   TR1-61-S47               -1   CTX+HV292.1     TR1 Human1 Human      191176
## 126 TR6-61-S186               -1     UNTREATED      CR Human1 Human      203803
## 10  TR1-33-S177                4 VAN+CCUG59168     TR1 Human1 Human      233850
## 124 TR6-59-S185               -3     UNTREATED      CR Human1 Human      235165
## 32   TR2-33-S30                4           VAN     TR2 Human1 Human      254676
## 49   TR3-28-S38               -1 VAN+CCUG59168     TR3 Human1 Human      698347
```

## rarefraction


```r
min_sample_size = 3820

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
## 867OTUs were removed because they are no longer 
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

![](NRP72_draft_16S_Figs_beta_PCoA_files/figure-html/unnamed-chunk-8-1.png)<!-- -->

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

![](NRP72_draft_16S_Figs_beta_PCoA_files/figure-html/unnamed-chunk-9-1.png)<!-- -->


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

![](NRP72_draft_16S_Figs_beta_PCoA_files/figure-html/unnamed-chunk-10-1.png)<!-- -->



```r
pcoa_all_noleg
```

![](NRP72_draft_16S_Figs_beta_PCoA_files/figure-html/unnamed-chunk-11-1.png)<!-- -->


```r
pcoa_all_noleg %>%
  export::graph2ppt(append = TRUE,
                    width = 317.48031496 * 1,
                    height = 0.618 * 317.48031496 * 1 , paper = "A4",  scaling = 2,
                    file = out_pptx)
```

```
## Exported graph as ~/Desktop/16S-human-beta.pptx
```

```r
pcoa_all_leg %>%
  export::graph2ppt(append = TRUE,
                    width = 317.48031496 * 1,
                    height = 0.618 * 317.48031496 * 1 , paper = "A4",  scaling = 2,
                    file = out_pptx)
```

```
## Exported graph as ~/Desktop/16S-human-beta.pptx
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
<div class="datatables html-widget html-fill-item-overflow-hidden html-fill-item" id="htmlwidget-9eceb20dc2a931752ae4" style="width:100%;height:auto;"></div>
<script type="application/json" data-for="htmlwidget-9eceb20dc2a931752ae4">{"x":{"filter":"none","vertical":false,"data":[["1","2"],["bjaccard","aichinson"],["Axis.1   [35.9%]","Axis.1   [40.4%]"],["bjaccard","aichinson"],["Axis.2   [17.4%]","Axis.2   [16.6%]"]],"container":"<table class=\"display\">\n  <thead>\n    <tr>\n      <th> <\/th>\n      <th>.id...1<\/th>\n      <th>V1...2<\/th>\n      <th>.id...3<\/th>\n      <th>V1...4<\/th>\n    <\/tr>\n  <\/thead>\n<\/table>","options":{"columnDefs":[{"orderable":false,"targets":0}],"order":[],"autoWidth":false,"orderClasses":false}},"evals":[],"jsHooks":[]}</script>
```

```r
# expl_var %>%
#   data.frame() %>% 
#   export::table2ppt(append = TRUE,
#                     file = out_pptx)
```

### Human1 - F1:


```r
ps_rare %>%
  subset_samples(Model2 == "Human1" & Treatment != "DONOR" & Fermentation == 1) -> ps_tmp_forbdiv
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
<div class="datatables html-widget html-fill-item-overflow-hidden html-fill-item" id="htmlwidget-f30d7598595affc4f25a" style="width:100%;height:auto;"></div>
<script type="application/json" data-for="htmlwidget-f30d7598595affc4f25a">{"x":{"filter":"none","vertical":false,"data":[["ASV0012","ASV0051","ASV0010","ASV0015","ASV0106","ASV0124","ASV0136","ASV0202"],[-0.852029448373223,0.642050293251154,-0.848848759256713,0.738897695216718,-0.744825118786255,-0.85595752617928,-0.750185320587873,-0.655417210461465],[-0.0669257982817142,0.305878300525776,-0.393558298515304,0.309141155976929,0.357193762563065,-0.100446732945435,-0.33673141072526,-0.305411815439601],[0.001,0.001,0.001,0.001,0.001,0.001,0.001,0.001],[0.730433243370824,0.505790113796429,0.875432350421923,0.64153805831533,0.682351841588919,0.742752832782365,0.676166058194553,0.522848096779201],["Bacteria","Bacteria","Bacteria","Bacteria","Bacteria","Bacteria","Bacteria","Bacteria"],["Firmicutes","Firmicutes","Firmicutes","Firmicutes","Firmicutes","Firmicutes","Actinobacteriota","Firmicutes"],["Clostridia","Bacilli","Clostridia","Negativicutes","Clostridia","Clostridia","Coriobacteriia","Clostridia"],["Oscillospirales","Lactobacillales","Lachnospirales","Veillonellales-Selenomonadales","Oscillospirales","Monoglobales","Coriobacteriales","Oscillospirales"],["Ruminococcaceae","Lactobacillaceae","Lachnospiraceae","Veillonellaceae","Oscillospiraceae","Monoglobaceae","Eggerthellaceae","Butyricicoccaceae"],[null,null,null,null,null,null,null,null],[null,null,null,null,null,null,null,null],[null,null,null,null,null,null,null,null],[0.0026,0.0026,0.0026,0.0026,0.0026,0.0026,0.0026,0.0026]],"container":"<table class=\"display\">\n  <thead>\n    <tr>\n      <th> <\/th>\n      <th>Axis.1<\/th>\n      <th>Axis.2<\/th>\n      <th>pval<\/th>\n      <th>r<\/th>\n      <th>Kingdom<\/th>\n      <th>Phylum<\/th>\n      <th>Class<\/th>\n      <th>Order<\/th>\n      <th>tax_rank_plot<\/th>\n      <th>Genus<\/th>\n      <th>Species<\/th>\n      <th>Strain<\/th>\n      <th>pval.adj<\/th>\n    <\/tr>\n  <\/thead>\n<\/table>","options":{"columnDefs":[{"className":"dt-right","targets":[1,2,3,4,13]},{"orderable":false,"targets":0}],"order":[],"autoWidth":false,"orderClasses":false}},"evals":[],"jsHooks":[]}</script>
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
<div class="datatables html-widget html-fill-item-overflow-hidden html-fill-item" id="htmlwidget-b85bf888dc62471906a9" style="width:100%;height:auto;"></div>
<script type="application/json" data-for="htmlwidget-b85bf888dc62471906a9">{"x":{"filter":"none","vertical":false,"data":[["1","2","3","4","5"],["Succinat_mM","Formiat_mM","Acetat_mM","Butyrat_mM","Total_SCFA_mM"],[0.252092250604304,0.397064413776484,-0.930780069116851,-0.885684192300746,-0.87703976046704],[0.30378683275291,-0.287737264321405,-0.160690158546535,-0.30320626888588,-0.0292257101633716],["Succinat_mM","Formiat_mM","Acetat_mM","Butyrat_mM","Total_SCFA_mM"],[0.003,0.001,0.001,0.001,0.001],[0.155836942568788,0.240452881966829,0.89217286411888,0.876370529983122,0.770052883574637]],"container":"<table class=\"display\">\n  <thead>\n    <tr>\n      <th> <\/th>\n      <th>label<\/th>\n      <th>Axis.1<\/th>\n      <th>Axis.2<\/th>\n      <th>id<\/th>\n      <th>pval<\/th>\n      <th>r<\/th>\n    <\/tr>\n  <\/thead>\n<\/table>","options":{"columnDefs":[{"className":"dt-right","targets":[2,3,5,6]},{"orderable":false,"targets":0}],"order":[],"autoWidth":false,"orderClasses":false}},"evals":[],"jsHooks":[]}</script>
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
<div class="datatables html-widget html-fill-item-overflow-hidden html-fill-item" id="htmlwidget-88fdd2a34d502a66744b" style="width:100%;height:auto;"></div>
<script type="application/json" data-for="htmlwidget-88fdd2a34d502a66744b">{"x":{"filter":"none","vertical":false,"data":[["ASV0012","ASV0051","ASV0010","ASV0015","ASV0106","ASV0124","ASV0136","ASV0144"],[-0.812062722028549,0.515929454096016,-0.787267254889976,0.688108706146924,-0.773557869016421,-0.74575333563041,-0.72889881230181,-0.683120289279889],[-0.227233951182378,0.594269592257247,-0.435826242147191,0.321728506176295,-0.13711311515029,-0.370385110158234,-0.342942982772263,0.217535637521298],[0.001,0.001,0.001,0.001,0.001,0.001,0.001,0.001],[0.711081133078371,0.619339549885408,0.80973424396614,0.577002823161624,0.617191783063443,0.69333316743081,0.648903368007726,0.513975083217637],["Bacteria","Bacteria","Bacteria","Bacteria","Bacteria","Bacteria","Bacteria","Bacteria"],["Firmicutes","Firmicutes","Firmicutes","Firmicutes","Firmicutes","Firmicutes","Actinobacteriota","Actinobacteriota"],["Clostridia","Bacilli","Clostridia","Negativicutes","Clostridia","Clostridia","Coriobacteriia","Coriobacteriia"],["Oscillospirales","Lactobacillales","Lachnospirales","Veillonellales-Selenomonadales","Oscillospirales","Monoglobales","Coriobacteriales","Coriobacteriales"],["Ruminococcaceae","Lactobacillaceae","Lachnospiraceae","Veillonellaceae","Oscillospiraceae","Monoglobaceae","Eggerthellaceae","Atopobiaceae"],[null,null,null,null,null,null,null,null],[null,null,null,null,null,null,null,null],[null,null,null,null,null,null,null,null],[0.00288888888888889,0.00288888888888889,0.00288888888888889,0.00288888888888889,0.00288888888888889,0.00288888888888889,0.00288888888888889,0.00288888888888889]],"container":"<table class=\"display\">\n  <thead>\n    <tr>\n      <th> <\/th>\n      <th>Axis.1<\/th>\n      <th>Axis.2<\/th>\n      <th>pval<\/th>\n      <th>r<\/th>\n      <th>Kingdom<\/th>\n      <th>Phylum<\/th>\n      <th>Class<\/th>\n      <th>Order<\/th>\n      <th>tax_rank_plot<\/th>\n      <th>Genus<\/th>\n      <th>Species<\/th>\n      <th>Strain<\/th>\n      <th>pval.adj<\/th>\n    <\/tr>\n  <\/thead>\n<\/table>","options":{"columnDefs":[{"className":"dt-right","targets":[1,2,3,4,13]},{"orderable":false,"targets":0}],"order":[],"autoWidth":false,"orderClasses":false}},"evals":[],"jsHooks":[]}</script>
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
<div class="datatables html-widget html-fill-item-overflow-hidden html-fill-item" id="htmlwidget-8e1fe96c0688a849169a" style="width:100%;height:auto;"></div>
<script type="application/json" data-for="htmlwidget-8e1fe96c0688a849169a">{"x":{"filter":"none","vertical":false,"data":[["1","2","3","4","5"],["Succinat_mM","Formiat_mM","Acetat_mM","Butyrat_mM","Total_SCFA_mM"],[0.164311218332464,0.441439686730363,-0.797826548818299,-0.751687244981628,-0.711900279912203],[0.506523475467576,-0.0242931772264333,-0.468604143606486,-0.51913836034085,-0.498311404715984],["Succinat_mM","Formiat_mM","Acetat_mM","Butyrat_mM","Total_SCFA_mM"],[0.001,0.004,0.001,0.001,0.001],[0.283564207669651,0.195459155480356,0.856117045404486,0.834538351445456,0.75511626460909]],"container":"<table class=\"display\">\n  <thead>\n    <tr>\n      <th> <\/th>\n      <th>label<\/th>\n      <th>Axis.1<\/th>\n      <th>Axis.2<\/th>\n      <th>id<\/th>\n      <th>pval<\/th>\n      <th>r<\/th>\n    <\/tr>\n  <\/thead>\n<\/table>","options":{"columnDefs":[{"className":"dt-right","targets":[2,3,5,6]},{"orderable":false,"targets":0}],"order":[],"autoWidth":false,"orderClasses":false}},"evals":[],"jsHooks":[]}</script>
```


```r
env_out$plot %>% 
  export::graph2ppt(append = TRUE,
                    width = 317.48031496 * 1,
                    height = 0.618 * 317.48031496 * 1 , paper = "A4",  scaling = 2,
                    file = out_pptx)
```

####  Statistical evaluation:

##### pret - Day_of_Treatment <= 0:


```r
lapply(
  dist,
  FUN = phyloseq_adonis2,
  physeq = ps_tmp_forbdiv %>%
    subset_samples(Day_of_Treatment <= 0),
  formula = paste0(c("Reactor_Treatment_Dose"), collapse=" + "),
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
<div class="datatables html-widget html-fill-item-overflow-hidden html-fill-item" id="htmlwidget-bb3486349adfa5d821dc" style="width:100%;height:auto;"></div>
<script type="application/json" data-for="htmlwidget-bb3486349adfa5d821dc">{"x":{"filter":"none","vertical":false,"data":[["1","2","3","4","5","6"],["bjaccard","bjaccard","bjaccard","aichinson","aichinson","aichinson"],["Reactor_Treatment_Dose","Residual","Total","Reactor_Treatment_Dose","Residual","Total"],[5,12,17,5,12,17],[0.518,0.281,0.799,0.714,0.427,1.142],[0.648,0.352,1,0.626,0.374,1],[4.414,null,null,4.012,null,null],[0.001,null,null,0.001,null,null]],"container":"<table class=\"display\">\n  <thead>\n    <tr>\n      <th> <\/th>\n      <th>Distance<\/th>\n      <th>terms<\/th>\n      <th>Df<\/th>\n      <th>SumOfSqs<\/th>\n      <th>R2<\/th>\n      <th>F<\/th>\n      <th>Pr(&gt;F)<\/th>\n    <\/tr>\n  <\/thead>\n<\/table>","options":{"columnDefs":[{"className":"dt-right","targets":[3,4,5,6,7]},{"orderable":false,"targets":0}],"order":[],"autoWidth":false,"orderClasses":false}},"evals":[],"jsHooks":[]}</script>
```

```r
tmp_df %>% 
  xlsx::write.xlsx(file = out_xlsx,
                   append = TRUE,
                   sheetName = "human1_F1-beta-pretreat",
                   showNA = TRUE)
```



```r
lapply(
  dist,
  FUN = phyloseq_adonis2,
  physeq = ps_tmp_forbdiv %>%
    subset_samples(Day_of_Treatment <= 0),
  formula = paste0(c("Reactor_Treatment_Dose", "Day_of_Treatment"), collapse=" + "),
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
<div class="datatables html-widget html-fill-item-overflow-hidden html-fill-item" id="htmlwidget-44474b3e8cd160f85a98" style="width:100%;height:auto;"></div>
<script type="application/json" data-for="htmlwidget-44474b3e8cd160f85a98">{"x":{"filter":"none","vertical":false,"data":[["1","2","3","4","5","6","7","8"],["bjaccard","bjaccard","bjaccard","bjaccard","aichinson","aichinson","aichinson","aichinson"],["Reactor_Treatment_Dose","Day_of_Treatment","Residual","Total","Reactor_Treatment_Dose","Day_of_Treatment","Residual","Total"],[5,1,11,17,5,1,11,17],[0.513,0.036,0.246,0.799,0.711,0.045,0.382,1.142],[0.643,0.044,0.308,1,0.623,0.04,0.335,1],[4.593,1.588,null,null,4.098,1.307,null,null],[0.001,0.117,null,null,0.001,0.176,null,null]],"container":"<table class=\"display\">\n  <thead>\n    <tr>\n      <th> <\/th>\n      <th>Distance<\/th>\n      <th>terms<\/th>\n      <th>Df<\/th>\n      <th>SumOfSqs<\/th>\n      <th>R2<\/th>\n      <th>F<\/th>\n      <th>Pr(&gt;F)<\/th>\n    <\/tr>\n  <\/thead>\n<\/table>","options":{"columnDefs":[{"className":"dt-right","targets":[3,4,5,6,7]},{"orderable":false,"targets":0}],"order":[],"autoWidth":false,"orderClasses":false}},"evals":[],"jsHooks":[]}</script>
```

```r
tmp_df %>% 
  xlsx::write.xlsx(file = out_xlsx,
                   append = TRUE,
                   sheetName = "human1_F1-beta-pretreat-time-day",
                   showNA = TRUE)
```


```r
lapply(
  dist,
  FUN = physeq_pairwise_permanovas_adonis2,
  physeq = ps_tmp_forbdiv %>%
    subset_samples(Day_of_Treatment <= 0),
  compare_header = "Reactor_Treatment_Dose",
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
## Set of permutations < 'minperm'. Generating entire set.
```

```r
tmp_df %>% 
  DT::datatable()
```

```{=html}
<div class="datatables html-widget html-fill-item-overflow-hidden html-fill-item" id="htmlwidget-2514ce28ae593bde0f73" style="width:100%;height:auto;"></div>
<script type="application/json" data-for="htmlwidget-2514ce28ae593bde0f73">{"x":{"filter":"none","vertical":false,"data":[["1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","23","24","25","26","27","28","29","30"],["bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson"],["CR_UNTREATED","CR_UNTREATED","CR_UNTREATED","CR_UNTREATED","CR_UNTREATED","TR2_VAN90","TR2_VAN90","TR2_VAN90","TR2_VAN90","TR4_VAN600","TR4_VAN600","TR4_VAN600","TR1_VAN+CCUG5916890","TR1_VAN+CCUG5916890","TR3_VAN+CCUG59168600","CR_UNTREATED","CR_UNTREATED","CR_UNTREATED","CR_UNTREATED","CR_UNTREATED","TR2_VAN90","TR2_VAN90","TR2_VAN90","TR2_VAN90","TR4_VAN600","TR4_VAN600","TR4_VAN600","TR1_VAN+CCUG5916890","TR1_VAN+CCUG5916890","TR3_VAN+CCUG59168600"],["TR2_VAN90","TR4_VAN600","TR1_VAN+CCUG5916890","TR3_VAN+CCUG59168600","TR5_CCUG59168","TR4_VAN600","TR1_VAN+CCUG5916890","TR3_VAN+CCUG59168600","TR5_CCUG59168","TR1_VAN+CCUG5916890","TR3_VAN+CCUG59168600","TR5_CCUG59168","TR3_VAN+CCUG59168600","TR5_CCUG59168","TR5_CCUG59168","TR2_VAN90","TR4_VAN600","TR1_VAN+CCUG5916890","TR3_VAN+CCUG59168600","TR5_CCUG59168","TR4_VAN600","TR1_VAN+CCUG5916890","TR3_VAN+CCUG59168600","TR5_CCUG59168","TR1_VAN+CCUG5916890","TR3_VAN+CCUG59168600","TR5_CCUG59168","TR3_VAN+CCUG59168600","TR5_CCUG59168","TR5_CCUG59168"],[0.53,0.528,0.391,0.669,0.483,0.603,0.424,0.501,0.354,0.368,0.603,0.375,0.589,0.351,0.495,0.509,0.443,0.392,0.606,0.44,0.619,0.432,0.476,0.38,0.457,0.594,0.411,0.562,0.324,0.478],[0.1,0.1,0.025,0.03,0.1,0.333,0.067,0.067,0.1,0.067,0.067,0.1,0.036,0.028,0.019,0.1,0.1,0.027,0.036,0.1,0.333,0.067,0.067,0.1,0.067,0.067,0.1,0.033,0.026,0.032],[0.125,0.125,0.188,0.112,0.125,0.333,0.133,0.133,0.125,0.133,0.133,0.125,0.108,0.14,0.285,0.125,0.125,0.202,0.108,0.125,0.333,0.133,0.133,0.125,0.133,0.133,0.125,0.124,0.39,0.16]],"container":"<table class=\"display\">\n  <thead>\n    <tr>\n      <th> <\/th>\n      <th>Distance<\/th>\n      <th>X1<\/th>\n      <th>X2<\/th>\n      <th>R2<\/th>\n      <th>pval<\/th>\n      <th>pvalFDR<\/th>\n    <\/tr>\n  <\/thead>\n<\/table>","options":{"columnDefs":[{"className":"dt-right","targets":[4,5,6]},{"orderable":false,"targets":0}],"order":[],"autoWidth":false,"orderClasses":false}},"evals":[],"jsHooks":[]}</script>
```

```r
tmp_df %>% 
  xlsx::write.xlsx(file = out_xlsx,
                   append = TRUE,
                   sheetName = "human1_F1-beta-pretreat-pw",
                   showNA = TRUE)
```


```r
lapply(
  dist,
  FUN = phyloseq_TW,
  physeq =  ps_tmp_forbdiv %>%
    subset_samples(Day_of_Treatment <= 0),
  variable = "Reactor_Treatment_Dose",
  nrep = 999,
  strata = NULL
) %>%
  bind_rows(.id = "Distance") %>%
  mutate_if(is.numeric, round, 3)   -> tmp_df

tmp_df %>% 
  DT::datatable()
```

```{=html}
<div class="datatables html-widget html-fill-item-overflow-hidden html-fill-item" id="htmlwidget-5d05eac75555269eefd0" style="width:100%;height:auto;"></div>
<script type="application/json" data-for="htmlwidget-5d05eac75555269eefd0">{"x":{"filter":"none","vertical":false,"data":[["1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","23","24","25","26","27","28","29","30"],["bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson"],[["CR_UNTREATED"],["CR_UNTREATED"],["CR_UNTREATED"],["CR_UNTREATED"],["CR_UNTREATED"],["TR2_VAN90"],["TR2_VAN90"],["TR2_VAN90"],["TR2_VAN90"],["TR4_VAN600"],["TR4_VAN600"],["TR4_VAN600"],["TR1_VAN+CCUG5916890"],["TR1_VAN+CCUG5916890"],["TR3_VAN+CCUG59168600"],["CR_UNTREATED"],["CR_UNTREATED"],["CR_UNTREATED"],["CR_UNTREATED"],["CR_UNTREATED"],["TR2_VAN90"],["TR2_VAN90"],["TR2_VAN90"],["TR2_VAN90"],["TR4_VAN600"],["TR4_VAN600"],["TR4_VAN600"],["TR1_VAN+CCUG5916890"],["TR1_VAN+CCUG5916890"],["TR3_VAN+CCUG59168600"]],[["TR2_VAN90"],["TR4_VAN600"],["TR1_VAN+CCUG5916890"],["TR3_VAN+CCUG59168600"],["TR5_CCUG59168"],["TR4_VAN600"],["TR1_VAN+CCUG5916890"],["TR3_VAN+CCUG59168600"],["TR5_CCUG59168"],["TR1_VAN+CCUG5916890"],["TR3_VAN+CCUG59168600"],["TR5_CCUG59168"],["TR3_VAN+CCUG59168600"],["TR5_CCUG59168"],["TR5_CCUG59168"],["TR2_VAN90"],["TR4_VAN600"],["TR1_VAN+CCUG5916890"],["TR3_VAN+CCUG59168600"],["TR5_CCUG59168"],["TR4_VAN600"],["TR1_VAN+CCUG5916890"],["TR3_VAN+CCUG59168600"],["TR5_CCUG59168"],["TR1_VAN+CCUG5916890"],["TR3_VAN+CCUG59168600"],["TR5_CCUG59168"],["TR3_VAN+CCUG59168600"],["TR5_CCUG59168"],["TR5_CCUG59168"]],[[3],[3],[3],[3],[3],[2],[2],[2],[2],[2],[2],[2],[4],[4],[4],[3],[3],[3],[3],[3],[2],[2],[2],[2],[2],[2],[2],[4],[4],[4]],[[2],[2],[4],[4],[3],[2],[4],[4],[3],[4],[4],[3],[4],[3],[3],[2],[2],[4],[4],[3],[2],[4],[4],[3],[4],[4],[3],[4],[3],[3]],[[0.104],[0.101],[0.032],[0.025],[0.107],[0.311],[0.147],[0.069],[0.101],[0.139],[0.074],[0.098],[0.025],[0.062],[0.037],[0.1],[0.095],[0.036],[0.023],[0.096],[0.338],[0.138],[0.064],[0.172],[0.059],[0.073],[0.091],[0.024],[0.024],[0.035]],[[2.80694537504184],[3.61431012691923],[3.15607571827432],[10.7354929735595],[3.73344949882427],[3.0375180067186],[2.1125144619168],[3.47022764688383],[1.49424093436232],[2.51746542213744],[7.85680875157764],[2.11173093505347],[8.60568172992972],[2.50623226705387],[4.91130049713801],[2.85687924637597],[2.54700864747287],[3.0238558533502],[7.7347542693938],[3.13989123864949],[3.25117874895107],[2.29083787455094],[3.22351485813272],[1.72978106453622],[3.21253348363019],[6.55814440142508],[2.29352429478639],[7.6852018254337],[2.21808851322542],[4.53563511530767]],[[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999]],[0.134,0.134,0.134,0.134,0.134,0.311,0.157,0.134,0.134,0.157,0.134,0.134,0.134,0.134,0.134,0.125,0.125,0.108,0.108,0.125,0.338,0.159,0.125,0.184,0.125,0.125,0.125,0.108,0.108,0.108]],"container":"<table class=\"display\">\n  <thead>\n    <tr>\n      <th> <\/th>\n      <th>Distance<\/th>\n      <th>Level1<\/th>\n      <th>Level2<\/th>\n      <th>N1<\/th>\n      <th>N2<\/th>\n      <th>p.value<\/th>\n      <th>tw2.stat<\/th>\n      <th>nrep<\/th>\n      <th>pvalFDR<\/th>\n    <\/tr>\n  <\/thead>\n<\/table>","options":{"columnDefs":[{"className":"dt-right","targets":9},{"orderable":false,"targets":0}],"order":[],"autoWidth":false,"orderClasses":false}},"evals":[],"jsHooks":[]}</script>
```

```r
tmp_df %>% 
  xlsx::write.xlsx(file = out_xlsx,
                   append = TRUE,
                   sheetName = "human1_F1-beta-pretreat-pw2",
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

##### treat - Day_of_Treatment > 1:


```r
lapply(
  dist,
  FUN = phyloseq_adonis2,
  physeq = ps_tmp_forbdiv %>%
    subset_samples(Day_of_Treatment >1),
  formula = paste0(c("Reactor_Treatment_Dose"), collapse=" + "),
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
<div class="datatables html-widget html-fill-item-overflow-hidden html-fill-item" id="htmlwidget-7531f021030c44f57ef2" style="width:100%;height:auto;"></div>
<script type="application/json" data-for="htmlwidget-7531f021030c44f57ef2">{"x":{"filter":"none","vertical":false,"data":[["1","2","3","4","5","6"],["bjaccard","bjaccard","bjaccard","aichinson","aichinson","aichinson"],["Reactor_Treatment_Dose","Residual","Total","Reactor_Treatment_Dose","Residual","Total"],[5,23,28,5,23,28],[4.291,2.246,6.537,2.871,1.084,3.955],[0.656,0.344,1,0.726,0.274,1],[8.79,null,null,12.179,null,null],[0.001,null,null,0.001,null,null]],"container":"<table class=\"display\">\n  <thead>\n    <tr>\n      <th> <\/th>\n      <th>Distance<\/th>\n      <th>terms<\/th>\n      <th>Df<\/th>\n      <th>SumOfSqs<\/th>\n      <th>R2<\/th>\n      <th>F<\/th>\n      <th>Pr(&gt;F)<\/th>\n    <\/tr>\n  <\/thead>\n<\/table>","options":{"columnDefs":[{"className":"dt-right","targets":[3,4,5,6,7]},{"orderable":false,"targets":0}],"order":[],"autoWidth":false,"orderClasses":false}},"evals":[],"jsHooks":[]}</script>
```

```r
tmp_df %>% 
  xlsx::write.xlsx(file = out_xlsx,
                   append = TRUE,
                   sheetName = "human1_F1-beta-treat",
                   showNA = TRUE)
```


```r
lapply(
  dist,
  FUN = phyloseq_adonis2,
  physeq = ps_tmp_forbdiv %>%
    subset_samples(Day_of_Treatment >1),
  formula = paste0(c("Reactor_Treatment_Dose", "Day_of_Treatment"), collapse=" + "),
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
<div class="datatables html-widget html-fill-item-overflow-hidden html-fill-item" id="htmlwidget-983315589fed18d622a7" style="width:100%;height:auto;"></div>
<script type="application/json" data-for="htmlwidget-983315589fed18d622a7">{"x":{"filter":"none","vertical":false,"data":[["1","2","3","4","5","6","7","8"],["bjaccard","bjaccard","bjaccard","bjaccard","aichinson","aichinson","aichinson","aichinson"],["Reactor_Treatment_Dose","Day_of_Treatment","Residual","Total","Reactor_Treatment_Dose","Day_of_Treatment","Residual","Total"],[5,1,22,28,5,1,22,28],[4.086,0.276,1.97,6.537,2.772,0.117,0.968,3.955],[0.625,0.042,0.301,1,0.701,0.029,0.245,1],[9.126,3.082,null,null,12.604,2.652,null,null],[0.001,0.028,null,null,0.001,0.062,null,null]],"container":"<table class=\"display\">\n  <thead>\n    <tr>\n      <th> <\/th>\n      <th>Distance<\/th>\n      <th>terms<\/th>\n      <th>Df<\/th>\n      <th>SumOfSqs<\/th>\n      <th>R2<\/th>\n      <th>F<\/th>\n      <th>Pr(&gt;F)<\/th>\n    <\/tr>\n  <\/thead>\n<\/table>","options":{"columnDefs":[{"className":"dt-right","targets":[3,4,5,6,7]},{"orderable":false,"targets":0}],"order":[],"autoWidth":false,"orderClasses":false}},"evals":[],"jsHooks":[]}</script>
```

```r
tmp_df %>% 
  xlsx::write.xlsx(file = out_xlsx,
                   append = TRUE,
                   sheetName = "human1_F1-beta-treat-time-day",
                   showNA = TRUE)
```


```r
lapply(
  dist,
  FUN = physeq_pairwise_permanovas_adonis2,
  physeq = ps_tmp_forbdiv %>%
    subset_samples(Day_of_Treatment > 1),
  compare_header = "Reactor_Treatment_Dose",
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
<div class="datatables html-widget html-fill-item-overflow-hidden html-fill-item" id="htmlwidget-6f1d8b9e786d9d3fe28d" style="width:100%;height:auto;"></div>
<script type="application/json" data-for="htmlwidget-6f1d8b9e786d9d3fe28d">{"x":{"filter":"none","vertical":false,"data":[["1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","23","24","25","26","27","28","29","30"],["bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson"],["CR_UNTREATED","CR_UNTREATED","CR_UNTREATED","CR_UNTREATED","CR_UNTREATED","TR2_VAN90","TR2_VAN90","TR2_VAN90","TR2_VAN90","TR4_VAN600","TR4_VAN600","TR4_VAN600","TR1_VAN+CCUG5916890","TR1_VAN+CCUG5916890","TR3_VAN+CCUG59168600","CR_UNTREATED","CR_UNTREATED","CR_UNTREATED","CR_UNTREATED","CR_UNTREATED","TR2_VAN90","TR2_VAN90","TR2_VAN90","TR2_VAN90","TR4_VAN600","TR4_VAN600","TR4_VAN600","TR1_VAN+CCUG5916890","TR1_VAN+CCUG5916890","TR3_VAN+CCUG59168600"],["TR2_VAN90","TR4_VAN600","TR1_VAN+CCUG5916890","TR3_VAN+CCUG59168600","TR5_CCUG59168","TR4_VAN600","TR1_VAN+CCUG5916890","TR3_VAN+CCUG59168600","TR5_CCUG59168","TR1_VAN+CCUG5916890","TR3_VAN+CCUG59168600","TR5_CCUG59168","TR3_VAN+CCUG59168600","TR5_CCUG59168","TR5_CCUG59168","TR2_VAN90","TR4_VAN600","TR1_VAN+CCUG5916890","TR3_VAN+CCUG59168600","TR5_CCUG59168","TR4_VAN600","TR1_VAN+CCUG5916890","TR3_VAN+CCUG59168600","TR5_CCUG59168","TR1_VAN+CCUG5916890","TR3_VAN+CCUG59168600","TR5_CCUG59168","TR3_VAN+CCUG59168600","TR5_CCUG59168","TR5_CCUG59168"],[0.636,0.496,0.864,0.776,0.36,0.205,0.264,0.398,0.63,0.305,0.295,0.491,0.545,0.854,0.764,0.698,0.631,0.837,0.807,0.344,0.215,0.265,0.406,0.687,0.317,0.311,0.616,0.584,0.821,0.788],[0.01,0.003,0.006,0.005,0.015,0.039,0.062,0.013,0.007,0.011,0.007,0.004,0.039,0.011,0.011,0.011,0.002,0.005,0.01,0.007,0.011,0.039,0.019,0.009,0.007,0.004,0.002,0.026,0.01,0.007],[0.021,0.045,0.022,0.025,0.019,0.043,0.062,0.018,0.019,0.018,0.019,0.03,0.043,0.018,0.018,0.014,0.02,0.019,0.016,0.018,0.014,0.039,0.022,0.017,0.018,0.02,0.02,0.028,0.016,0.018]],"container":"<table class=\"display\">\n  <thead>\n    <tr>\n      <th> <\/th>\n      <th>Distance<\/th>\n      <th>X1<\/th>\n      <th>X2<\/th>\n      <th>R2<\/th>\n      <th>pval<\/th>\n      <th>pvalFDR<\/th>\n    <\/tr>\n  <\/thead>\n<\/table>","options":{"columnDefs":[{"className":"dt-right","targets":[4,5,6]},{"orderable":false,"targets":0}],"order":[],"autoWidth":false,"orderClasses":false}},"evals":[],"jsHooks":[]}</script>
```

```r
tmp_df %>% 
  xlsx::write.xlsx(file = out_xlsx,
                   append = TRUE,
                   sheetName = "human1_F1-beta-treat-pw",
                   showNA = TRUE)
```


```r
lapply(
  dist,
  FUN = phyloseq_TW,
  physeq =  ps_tmp_forbdiv %>%
    subset_samples(Day_of_Treatment >1),
  variable = "Reactor_Treatment_Dose",
  nrep = 999,
  strata = NULL
) %>%
  bind_rows(.id = "Distance") %>%
  mutate_if(is.numeric, round, 3)   -> tmp_df

tmp_df %>% 
  DT::datatable()
```

```{=html}
<div class="datatables html-widget html-fill-item-overflow-hidden html-fill-item" id="htmlwidget-dcdfc64296f5794e54a4" style="width:100%;height:auto;"></div>
<script type="application/json" data-for="htmlwidget-dcdfc64296f5794e54a4">{"x":{"filter":"none","vertical":false,"data":[["1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","23","24","25","26","27","28","29","30"],["bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson"],[["CR_UNTREATED"],["CR_UNTREATED"],["CR_UNTREATED"],["CR_UNTREATED"],["CR_UNTREATED"],["TR2_VAN90"],["TR2_VAN90"],["TR2_VAN90"],["TR2_VAN90"],["TR4_VAN600"],["TR4_VAN600"],["TR4_VAN600"],["TR1_VAN+CCUG5916890"],["TR1_VAN+CCUG5916890"],["TR3_VAN+CCUG59168600"],["CR_UNTREATED"],["CR_UNTREATED"],["CR_UNTREATED"],["CR_UNTREATED"],["CR_UNTREATED"],["TR2_VAN90"],["TR2_VAN90"],["TR2_VAN90"],["TR2_VAN90"],["TR4_VAN600"],["TR4_VAN600"],["TR4_VAN600"],["TR1_VAN+CCUG5916890"],["TR1_VAN+CCUG5916890"],["TR3_VAN+CCUG59168600"]],[["TR2_VAN90"],["TR4_VAN600"],["TR1_VAN+CCUG5916890"],["TR3_VAN+CCUG59168600"],["TR5_CCUG59168"],["TR4_VAN600"],["TR1_VAN+CCUG5916890"],["TR3_VAN+CCUG59168600"],["TR5_CCUG59168"],["TR1_VAN+CCUG5916890"],["TR3_VAN+CCUG59168600"],["TR5_CCUG59168"],["TR3_VAN+CCUG59168600"],["TR5_CCUG59168"],["TR5_CCUG59168"],["TR2_VAN90"],["TR4_VAN600"],["TR1_VAN+CCUG5916890"],["TR3_VAN+CCUG59168600"],["TR5_CCUG59168"],["TR4_VAN600"],["TR1_VAN+CCUG5916890"],["TR3_VAN+CCUG59168600"],["TR5_CCUG59168"],["TR1_VAN+CCUG5916890"],["TR3_VAN+CCUG59168600"],["TR5_CCUG59168"],["TR3_VAN+CCUG59168600"],["TR5_CCUG59168"],["TR5_CCUG59168"]],[[5],[5],[5],[5],[5],[5],[5],[5],[5],[6],[6],[6],[4],[4],[4],[5],[5],[5],[5],[5],[5],[5],[5],[5],[6],[6],[6],[4],[4],[4]],[[5],[6],[4],[4],[5],[6],[4],[4],[5],[4],[4],[5],[4],[5],[5],[5],[6],[4],[4],[5],[6],[4],[4],[5],[4],[4],[5],[4],[5],[5]],[[0.009],[0.003],[0.008],[0.016],[0.007],[0.02],[0.046],[0.014],[0.006],[0.004],[0.009],[0.002],[0.032],[0.01],[0.008],[0.006],[0.003],[0.007],[0.009],[0.004],[0.007],[0.03],[0.009],[0.012],[0.007],[0.004],[0.005],[0.029],[0.009],[0.006]],[[14.0013082789181],[10.5081855562512],[38.5183035041101],[19.9715455095308],[4.50376066993907],[2.39450736648007],[2.82887216695764],[4.83005319458075],[13.6356638778626],[4.54328234191308],[3.85245342133202],[10.231592703],[7.17690713170759],[36.261911234296],[19.0113816463892],[18.4734014574025],[16.8957007552249],[36.4758097867559],[28.7207331924546],[4.18802840816871],[2.51704331601202],[2.81402379273266],[5.18612450205059],[17.5731536370047],[4.66277166648887],[4.32451256781246],[15.640993365667],[8.40699307799367],[33.2760906231791],[26.0592266427036]],[[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999]],[0.015,0.015,0.015,0.02,0.015,0.023,0.046,0.019,0.015,0.015,0.015,0.015,0.034,0.015,0.015,0.011,0.011,0.011,0.011,0.011,0.011,0.03,0.011,0.014,0.011,0.011,0.011,0.03,0.011,0.011]],"container":"<table class=\"display\">\n  <thead>\n    <tr>\n      <th> <\/th>\n      <th>Distance<\/th>\n      <th>Level1<\/th>\n      <th>Level2<\/th>\n      <th>N1<\/th>\n      <th>N2<\/th>\n      <th>p.value<\/th>\n      <th>tw2.stat<\/th>\n      <th>nrep<\/th>\n      <th>pvalFDR<\/th>\n    <\/tr>\n  <\/thead>\n<\/table>","options":{"columnDefs":[{"className":"dt-right","targets":9},{"orderable":false,"targets":0}],"order":[],"autoWidth":false,"orderClasses":false}},"evals":[],"jsHooks":[]}</script>
```

```r
tmp_df %>% 
  xlsx::write.xlsx(file = out_xlsx,
                   append = TRUE,
                   sheetName = "human1_F1-beta-treat-pw2",
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

##### Day_of_Treatment >= 3 & Day_of_Treatment <= 6:


```r
lapply(
  dist,
  FUN = phyloseq_adonis2,
  physeq = ps_tmp_forbdiv %>%
    subset_samples(Day_of_Treatment >= 3 & Day_of_Treatment <= 6),
  formula = paste0(c("Reactor_Treatment_Dose"), collapse=" + "),
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
<div class="datatables html-widget html-fill-item-overflow-hidden html-fill-item" id="htmlwidget-873249b169eb6149d57e" style="width:100%;height:auto;"></div>
<script type="application/json" data-for="htmlwidget-873249b169eb6149d57e">{"x":{"filter":"none","vertical":false,"data":[["1","2","3","4","5","6"],["bjaccard","bjaccard","bjaccard","aichinson","aichinson","aichinson"],["Reactor_Treatment_Dose","Residual","Total","Reactor_Treatment_Dose","Residual","Total"],[5,13,18,5,13,18],[3.09,1.212,4.302,2.025,0.605,2.629],[0.718,0.282,1,0.77,0.23,1],[6.632,null,null,8.705,null,null],[0.001,null,null,0.001,null,null]],"container":"<table class=\"display\">\n  <thead>\n    <tr>\n      <th> <\/th>\n      <th>Distance<\/th>\n      <th>terms<\/th>\n      <th>Df<\/th>\n      <th>SumOfSqs<\/th>\n      <th>R2<\/th>\n      <th>F<\/th>\n      <th>Pr(&gt;F)<\/th>\n    <\/tr>\n  <\/thead>\n<\/table>","options":{"columnDefs":[{"className":"dt-right","targets":[3,4,5,6,7]},{"orderable":false,"targets":0}],"order":[],"autoWidth":false,"orderClasses":false}},"evals":[],"jsHooks":[]}</script>
```

```r
tmp_df %>% 
  xlsx::write.xlsx(file = out_xlsx,
                   append = TRUE,
                   sheetName = "human1_F1-beta-d3d6",
                   showNA = TRUE)
```


```r
lapply(
  dist,
  FUN = phyloseq_adonis2,
  physeq = ps_tmp_forbdiv %>%
    subset_samples(Day_of_Treatment >= 3 & Day_of_Treatment <= 6),
  formula = paste0(c("Reactor_Treatment_Dose", "Day_of_Treatment"), collapse=" + "),
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
<div class="datatables html-widget html-fill-item-overflow-hidden html-fill-item" id="htmlwidget-cbd3ad10bd2d35a2e76d" style="width:100%;height:auto;"></div>
<script type="application/json" data-for="htmlwidget-cbd3ad10bd2d35a2e76d">{"x":{"filter":"none","vertical":false,"data":[["1","2","3","4","5","6","7","8"],["bjaccard","bjaccard","bjaccard","bjaccard","aichinson","aichinson","aichinson","aichinson"],["Reactor_Treatment_Dose","Day_of_Treatment","Residual","Total","Reactor_Treatment_Dose","Day_of_Treatment","Residual","Total"],[5,1,12,18,5,1,12,18],[3.055,0.152,1.06,4.302,2.008,0.075,0.53,2.629],[0.71,0.035,0.246,1,0.764,0.029,0.201,1],[6.919,1.72,null,null,9.1,1.702,null,null],[0.001,0.128,null,null,0.001,0.141,null,null]],"container":"<table class=\"display\">\n  <thead>\n    <tr>\n      <th> <\/th>\n      <th>Distance<\/th>\n      <th>terms<\/th>\n      <th>Df<\/th>\n      <th>SumOfSqs<\/th>\n      <th>R2<\/th>\n      <th>F<\/th>\n      <th>Pr(&gt;F)<\/th>\n    <\/tr>\n  <\/thead>\n<\/table>","options":{"columnDefs":[{"className":"dt-right","targets":[3,4,5,6,7]},{"orderable":false,"targets":0}],"order":[],"autoWidth":false,"orderClasses":false}},"evals":[],"jsHooks":[]}</script>
```

```r
tmp_df %>% 
  xlsx::write.xlsx(file = out_xlsx,
                   append = TRUE,
                   sheetName = "human1_F1-beta-d3d6-time-day",
                   showNA = TRUE)
```


```r
lapply(
  dist,
  FUN = physeq_pairwise_permanovas_adonis2,
  physeq = ps_tmp_forbdiv %>%
    subset_samples(Day_of_Treatment >= 3 & Day_of_Treatment <= 6),
  compare_header = "Reactor_Treatment_Dose",
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

```r
tmp_df %>% 
  DT::datatable()
```

```{=html}
<div class="datatables html-widget html-fill-item-overflow-hidden html-fill-item" id="htmlwidget-f25d4e1770d66a86dcce" style="width:100%;height:auto;"></div>
<script type="application/json" data-for="htmlwidget-f25d4e1770d66a86dcce">{"x":{"filter":"none","vertical":false,"data":[["1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","23","24","25","26","27","28","29","30"],["bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson"],["CR_UNTREATED","CR_UNTREATED","CR_UNTREATED","CR_UNTREATED","CR_UNTREATED","TR2_VAN90","TR2_VAN90","TR2_VAN90","TR2_VAN90","TR4_VAN600","TR4_VAN600","TR4_VAN600","TR1_VAN+CCUG5916890","TR1_VAN+CCUG5916890","TR3_VAN+CCUG59168600","CR_UNTREATED","CR_UNTREATED","CR_UNTREATED","CR_UNTREATED","CR_UNTREATED","TR2_VAN90","TR2_VAN90","TR2_VAN90","TR2_VAN90","TR4_VAN600","TR4_VAN600","TR4_VAN600","TR1_VAN+CCUG5916890","TR1_VAN+CCUG5916890","TR3_VAN+CCUG59168600"],["TR2_VAN90","TR4_VAN600","TR1_VAN+CCUG5916890","TR3_VAN+CCUG59168600","TR5_CCUG59168","TR4_VAN600","TR1_VAN+CCUG5916890","TR3_VAN+CCUG59168600","TR5_CCUG59168","TR1_VAN+CCUG5916890","TR3_VAN+CCUG59168600","TR5_CCUG59168","TR3_VAN+CCUG59168600","TR5_CCUG59168","TR5_CCUG59168","TR2_VAN90","TR4_VAN600","TR1_VAN+CCUG5916890","TR3_VAN+CCUG59168600","TR5_CCUG59168","TR4_VAN600","TR1_VAN+CCUG5916890","TR3_VAN+CCUG59168600","TR5_CCUG59168","TR1_VAN+CCUG5916890","TR3_VAN+CCUG59168600","TR5_CCUG59168","TR3_VAN+CCUG59168600","TR5_CCUG59168","TR5_CCUG59168"],[0.789,0.519,0.866,0.863,0.507,0.286,0.346,0.588,0.777,0.335,0.326,0.486,0.634,0.862,0.863,0.781,0.637,0.839,0.844,0.453,0.282,0.339,0.601,0.791,0.328,0.322,0.62,0.668,0.852,0.868],[0.023,0.037,0.03,0.067,0.028,0.055,0.1,0.1,0.1,0.03,0.067,0.033,0.1,0.1,0.1,0.018,0.032,0.027,0.067,0.032,0.068,0.1,0.1,0.1,0.028,0.133,0.029,0.1,0.1,0.1],[0.345,0.092,0.129,0.118,0.21,0.118,0.12,0.12,0.12,0.129,0.118,0.099,0.12,0.12,0.12,0.27,0.087,0.202,0.143,0.087,0.128,0.13,0.13,0.13,0.14,0.133,0.109,0.13,0.13,0.13]],"container":"<table class=\"display\">\n  <thead>\n    <tr>\n      <th> <\/th>\n      <th>Distance<\/th>\n      <th>X1<\/th>\n      <th>X2<\/th>\n      <th>R2<\/th>\n      <th>pval<\/th>\n      <th>pvalFDR<\/th>\n    <\/tr>\n  <\/thead>\n<\/table>","options":{"columnDefs":[{"className":"dt-right","targets":[4,5,6]},{"orderable":false,"targets":0}],"order":[],"autoWidth":false,"orderClasses":false}},"evals":[],"jsHooks":[]}</script>
```

```r
tmp_df %>% 
  xlsx::write.xlsx(file = out_xlsx,
                   append = TRUE,
                   sheetName = "human1_F1-beta-d3d6-pw",
                   showNA = TRUE)
```


```r
lapply(
  dist,
  FUN = phyloseq_TW,
  physeq =  ps_tmp_forbdiv %>%
    subset_samples(Day_of_Treatment >= 3 & Day_of_Treatment <= 6),
  variable = "Reactor_Treatment_Dose",
  nrep = 999,
  strata = NULL
) %>%
  bind_rows(.id = "Distance") %>%
  mutate_if(is.numeric, round, 3)   -> tmp_df

tmp_df %>% 
  DT::datatable()
```

```{=html}
<div class="datatables html-widget html-fill-item-overflow-hidden html-fill-item" id="htmlwidget-969965edfe9f869aa904" style="width:100%;height:auto;"></div>
<script type="application/json" data-for="htmlwidget-969965edfe9f869aa904">{"x":{"filter":"none","vertical":false,"data":[["1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","23","24","25","26","27","28","29","30"],["bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson"],[["CR_UNTREATED"],["CR_UNTREATED"],["CR_UNTREATED"],["CR_UNTREATED"],["CR_UNTREATED"],["TR2_VAN90"],["TR2_VAN90"],["TR2_VAN90"],["TR2_VAN90"],["TR4_VAN600"],["TR4_VAN600"],["TR4_VAN600"],["TR1_VAN+CCUG5916890"],["TR1_VAN+CCUG5916890"],["TR3_VAN+CCUG59168600"],["CR_UNTREATED"],["CR_UNTREATED"],["CR_UNTREATED"],["CR_UNTREATED"],["CR_UNTREATED"],["TR2_VAN90"],["TR2_VAN90"],["TR2_VAN90"],["TR2_VAN90"],["TR4_VAN600"],["TR4_VAN600"],["TR4_VAN600"],["TR1_VAN+CCUG5916890"],["TR1_VAN+CCUG5916890"],["TR3_VAN+CCUG59168600"]],[["TR2_VAN90"],["TR4_VAN600"],["TR1_VAN+CCUG5916890"],["TR3_VAN+CCUG59168600"],["TR5_CCUG59168"],["TR4_VAN600"],["TR1_VAN+CCUG5916890"],["TR3_VAN+CCUG59168600"],["TR5_CCUG59168"],["TR1_VAN+CCUG5916890"],["TR3_VAN+CCUG59168600"],["TR5_CCUG59168"],["TR3_VAN+CCUG59168600"],["TR5_CCUG59168"],["TR5_CCUG59168"],["TR2_VAN90"],["TR4_VAN600"],["TR1_VAN+CCUG5916890"],["TR3_VAN+CCUG59168600"],["TR5_CCUG59168"],["TR4_VAN600"],["TR1_VAN+CCUG5916890"],["TR3_VAN+CCUG59168600"],["TR5_CCUG59168"],["TR1_VAN+CCUG5916890"],["TR3_VAN+CCUG59168600"],["TR5_CCUG59168"],["TR3_VAN+CCUG59168600"],["TR5_CCUG59168"],["TR5_CCUG59168"]],[[4],[4],[4],[4],[4],[3],[3],[3],[3],[4],[4],[4],[3],[3],[2],[4],[4],[4],[4],[4],[3],[3],[3],[3],[4],[4],[4],[3],[3],[2]],[[3],[4],[3],[2],[3],[4],[3],[2],[3],[3],[2],[3],[2],[3],[3],[3],[4],[3],[2],[3],[4],[3],[2],[3],[3],[2],[3],[2],[3],[3]],[[0.032],[0.025],[0.028],[0.077],[0.017],[0.057],[0.103],[0.097],[0.08],[0.02],[0.07],[0.042],[0.098],[0.091],[0.098],[0.029],[0.034],[0.02],[0.057],[0.03],[0.061],[0.102],[0.108],[0.104],[0.033],[0.056],[0.026],[0.097],[0.087],[0.096]],[[14.209059650565],[6.46582184638286],[25.6876843810569],[12.5362469077493],[5.18037871223416],[2.28650188568313],[2.11947941949056],[4.37106694797],[13.9117673336556],[3.05979519097538],[2.73130127793387],[6.41060295126254],[4.69166398524375],[24.9580955917854],[12.1941933703433],[16.3194716621413],[10.5310006001538],[25.7083805476784],[24.7383726940065],[4.12294342827623],[2.20404979513458],[2.05016651208964],[5.5480106254557],[15.171644886578],[2.91736678158837],[3.22110903170035],[9.81165621138822],[6.69866439906215],[22.9437718724269],[21.6529736494397]],[[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999]],[0.096,0.096,0.096,0.103,0.096,0.103,0.103,0.103,0.103,0.096,0.103,0.103,0.103,0.103,0.103,0.085,0.085,0.085,0.102,0.085,0.102,0.108,0.108,0.108,0.085,0.102,0.085,0.108,0.108,0.108]],"container":"<table class=\"display\">\n  <thead>\n    <tr>\n      <th> <\/th>\n      <th>Distance<\/th>\n      <th>Level1<\/th>\n      <th>Level2<\/th>\n      <th>N1<\/th>\n      <th>N2<\/th>\n      <th>p.value<\/th>\n      <th>tw2.stat<\/th>\n      <th>nrep<\/th>\n      <th>pvalFDR<\/th>\n    <\/tr>\n  <\/thead>\n<\/table>","options":{"columnDefs":[{"className":"dt-right","targets":9},{"orderable":false,"targets":0}],"order":[],"autoWidth":false,"orderClasses":false}},"evals":[],"jsHooks":[]}</script>
```

```r
tmp_df %>% 
  xlsx::write.xlsx(file = out_xlsx,
                   append = TRUE,
                   sheetName = "human1_F1-beta-d3d6-pw2",
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


##### Period:


```r
setNames(vector("list", 
                length(physeq_get_unique(ps_tmp_forbdiv, "Reactor_Treatment_Dose"))), 
         physeq_get_unique(ps_tmp_forbdiv, "Reactor_Treatment_Dose")) -> lis
```



```r
for (grp in names(lis)){
  print(grp)
  
  lapply(
    dist,
    FUN = physeq_pairwise_permanovas_adonis2,
    physeq = ps_tmp_forbdiv %>%
      # subset_samples(.,  eval(as.name(ranks)) %in% tax) %>% 
      subset_samples(., Reactor_Treatment_Dose == grp),
    compare_header = "period",
    n_perm = 999,
    strat = "none"
  ) %>%
    bind_rows(.id = "Distance") %>%
    mutate_if(is.numeric, round, 3) %>%
    select(-pvalBon) -> lis[[grp]]
}
```

```
## [1] "TR4_VAN600"
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
```

```
## [1] "TR2_VAN90"
```

```
## 'nperm' >= set of all permutations: complete enumeration.
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
## [1] "TR3_VAN+CCUG59168600"
```

```
## 'nperm' >= set of all permutations: complete enumeration.
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
## [1] "TR1_VAN+CCUG5916890"
```

```
## 'nperm' >= set of all permutations: complete enumeration.
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
## [1] "CR_UNTREATED"
```

```
## 'nperm' >= set of all permutations: complete enumeration.
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
```

```
## 'nperm' >= set of all permutations: complete enumeration.
```

```
## Set of permutations < 'minperm'. Generating entire set.
```

```
## [1] "TR5_CCUG59168"
```

```
## 'nperm' >= set of all permutations: complete enumeration.
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

```r
lis %>% 
  bind_rows(.id = "Reactor_Treatment_Dose") %>% 
  filter(X1 == "pret" | "X2" == "pret") %>% 
  group_by(Reactor_Treatment_Dose, Distance) %>% 
  rstatix::adjust_pvalue(p.col = "pval", output.col = "pvalFDR") %>% 
  rstatix::add_significance(p.col = "pvalFDR")-> tmp_df

tmp_df %>% 
  DT::datatable()
```

```{=html}
<div class="datatables html-widget html-fill-item-overflow-hidden html-fill-item" id="htmlwidget-450f171f14364aa3ebed" style="width:100%;height:auto;"></div>
<script type="application/json" data-for="htmlwidget-450f171f14364aa3ebed">{"x":{"filter":"none","vertical":false,"data":[["1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","23","24","25","26","27","28","29","30","31","32","33","34"],["CR_UNTREATED","CR_UNTREATED","CR_UNTREATED","CR_UNTREATED","TR1_VAN+CCUG5916890","TR1_VAN+CCUG5916890","TR1_VAN+CCUG5916890","TR1_VAN+CCUG5916890","TR1_VAN+CCUG5916890","TR1_VAN+CCUG5916890","TR2_VAN90","TR2_VAN90","TR2_VAN90","TR2_VAN90","TR2_VAN90","TR2_VAN90","TR3_VAN+CCUG59168600","TR3_VAN+CCUG59168600","TR3_VAN+CCUG59168600","TR3_VAN+CCUG59168600","TR3_VAN+CCUG59168600","TR3_VAN+CCUG59168600","TR4_VAN600","TR4_VAN600","TR4_VAN600","TR4_VAN600","TR4_VAN600","TR4_VAN600","TR5_CCUG59168","TR5_CCUG59168","TR5_CCUG59168","TR5_CCUG59168","TR5_CCUG59168","TR5_CCUG59168"],["aichinson","aichinson","bjaccard","bjaccard","aichinson","aichinson","aichinson","bjaccard","bjaccard","bjaccard","aichinson","aichinson","aichinson","bjaccard","bjaccard","bjaccard","aichinson","aichinson","aichinson","bjaccard","bjaccard","bjaccard","aichinson","aichinson","aichinson","bjaccard","bjaccard","bjaccard","aichinson","aichinson","aichinson","bjaccard","bjaccard","bjaccard"],["pret","pret","pret","pret","pret","pret","pret","pret","pret","pret","pret","pret","pret","pret","pret","pret","pret","pret","pret","pret","pret","pret","pret","pret","pret","pret","pret","pret","pret","pret","pret","pret","pret","pret"],["t1","t2","t1","t2","t1","t2","t3","t1","t2","t3","t1","t2","t3","t1","t2","t3","t1","t2","t3","t1","t2","t3","t1","t2","t3","t1","t2","t3","t1","t2","t3","t1","t2","t3"],[0.449,0.384,0.548,0.486,0.657,0.841,0.829,0.718,0.855,0.912,0.606,0.727,0.851,0.689,0.711,0.924,0.633,0.771,0.737,0.679,0.805,0.868,0.602,0.584,0.912,0.609,0.422,0.971,0.411,0.454,0.442,0.36,0.43,0.39],[0.2,0.032,0.1,0.029,0.2,0.032,0.2,0.2,0.032,0.2,0.333,0.1,0.333,0.333,0.1,0.333,0.067,0.067,0.2,0.067,0.067,0.2,0.333,0.067,0.333,0.333,0.067,0.333,0.1,0.1,0.5,0.2,0.1,0.5],[0.2,0.064,0.1,0.058,0.4,0.096,0.4,0.4,0.096,0.4,0.666,0.3,0.666,0.666,0.3,0.666,0.201,0.201,0.201,0.201,0.201,0.201,0.666,0.201,0.666,0.666,0.201,0.666,0.3,0.3,0.5,0.4,0.3,0.5],["ns","ns","ns","ns","ns","ns","ns","ns","ns","ns","ns","ns","ns","ns","ns","ns","ns","ns","ns","ns","ns","ns","ns","ns","ns","ns","ns","ns","ns","ns","ns","ns","ns","ns"]],"container":"<table class=\"display\">\n  <thead>\n    <tr>\n      <th> <\/th>\n      <th>Reactor_Treatment_Dose<\/th>\n      <th>Distance<\/th>\n      <th>X1<\/th>\n      <th>X2<\/th>\n      <th>R2<\/th>\n      <th>pval<\/th>\n      <th>pvalFDR<\/th>\n      <th>pvalFDR.signif<\/th>\n    <\/tr>\n  <\/thead>\n<\/table>","options":{"columnDefs":[{"className":"dt-right","targets":[5,6,7]},{"orderable":false,"targets":0}],"order":[],"autoWidth":false,"orderClasses":false}},"evals":[],"jsHooks":[]}</script>
```

```r
tmp_df %>% 
  xlsx::write.xlsx(file = out_xlsx,
                   append = TRUE,
                   sheetName = "human1_F1-beta-vs_pret",
                   showNA = TRUE)
```



```r
lis %>% 
  bind_rows(.id = "Reactor_Treatment_Dose") %>% 
  group_by(Reactor_Treatment_Dose, Distance) %>% 
  rstatix::adjust_pvalue(p.col = "pval", output.col = "pvalFDR") %>% 
  rstatix::add_significance(p.col = "pvalFDR") -> tmp_df

tmp_df %>% 
  DT::datatable()
```

```{=html}
<div class="datatables html-widget html-fill-item-overflow-hidden html-fill-item" id="htmlwidget-fbfc52d82635d988c6e1" style="width:100%;height:auto;"></div>
<script type="application/json" data-for="htmlwidget-fbfc52d82635d988c6e1">{"x":{"filter":"none","vertical":false,"data":[["1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","23","24","25","26","27","28","29","30","31","32","33","34","35","36","37","38","39","40","41","42","43","44","45","46","47","48","49","50","51","52","53","54","55","56","57","58","59","60","61","62","63","64","65","66"],["CR_UNTREATED","CR_UNTREATED","CR_UNTREATED","CR_UNTREATED","CR_UNTREATED","CR_UNTREATED","TR1_VAN+CCUG5916890","TR1_VAN+CCUG5916890","TR1_VAN+CCUG5916890","TR1_VAN+CCUG5916890","TR1_VAN+CCUG5916890","TR1_VAN+CCUG5916890","TR1_VAN+CCUG5916890","TR1_VAN+CCUG5916890","TR1_VAN+CCUG5916890","TR1_VAN+CCUG5916890","TR1_VAN+CCUG5916890","TR1_VAN+CCUG5916890","TR2_VAN90","TR2_VAN90","TR2_VAN90","TR2_VAN90","TR2_VAN90","TR2_VAN90","TR2_VAN90","TR2_VAN90","TR2_VAN90","TR2_VAN90","TR2_VAN90","TR2_VAN90","TR3_VAN+CCUG59168600","TR3_VAN+CCUG59168600","TR3_VAN+CCUG59168600","TR3_VAN+CCUG59168600","TR3_VAN+CCUG59168600","TR3_VAN+CCUG59168600","TR3_VAN+CCUG59168600","TR3_VAN+CCUG59168600","TR3_VAN+CCUG59168600","TR3_VAN+CCUG59168600","TR3_VAN+CCUG59168600","TR3_VAN+CCUG59168600","TR4_VAN600","TR4_VAN600","TR4_VAN600","TR4_VAN600","TR4_VAN600","TR4_VAN600","TR4_VAN600","TR4_VAN600","TR4_VAN600","TR4_VAN600","TR4_VAN600","TR4_VAN600","TR5_CCUG59168","TR5_CCUG59168","TR5_CCUG59168","TR5_CCUG59168","TR5_CCUG59168","TR5_CCUG59168","TR5_CCUG59168","TR5_CCUG59168","TR5_CCUG59168","TR5_CCUG59168","TR5_CCUG59168","TR5_CCUG59168"],["aichinson","aichinson","aichinson","bjaccard","bjaccard","bjaccard","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard"],["pret","pret","t1","pret","pret","t1","pret","pret","pret","t1","t1","t2","pret","pret","pret","t1","t1","t2","pret","pret","pret","t1","t1","t2","pret","pret","pret","t1","t1","t2","pret","pret","pret","t1","t1","t2","pret","pret","pret","t1","t1","t2","pret","pret","pret","t1","t1","t2","pret","pret","pret","t1","t1","t2","pret","pret","pret","t1","t1","t2","pret","pret","pret","t1","t1","t2"],["t1","t2","t2","t1","t2","t2","t1","t2","t3","t2","t3","t3","t1","t2","t3","t2","t3","t3","t1","t2","t3","t2","t3","t3","t1","t2","t3","t2","t3","t3","t1","t2","t3","t2","t3","t3","t1","t2","t3","t2","t3","t3","t1","t2","t3","t2","t3","t3","t1","t2","t3","t2","t3","t3","t1","t2","t3","t2","t3","t3","t1","t2","t3","t2","t3","t3"],[0.449,0.384,0.345,0.548,0.486,0.437,0.657,0.841,0.829,0.745,1,0.223,0.718,0.855,0.912,0.735,1,0.209,0.606,0.727,0.851,0.563,0.772,0.299,0.689,0.711,0.924,0.607,0.887,0.338,0.633,0.771,0.737,0.574,0.696,0.617,0.679,0.805,0.868,0.542,0.725,0.472,0.602,0.584,0.912,0.308,0.618,0.116,0.609,0.422,0.971,0.248,0.765,0.174,0.411,0.454,0.442,0.295,0.563,0.421,0.36,0.43,0.39,0.333,0.471,0.456],[0.2,0.032,0.067,0.1,0.029,0.067,0.2,0.032,0.2,0.25,null,1,0.2,0.032,0.2,0.25,null,1,0.333,0.1,0.333,0.1,0.333,0.75,0.333,0.1,0.333,0.1,0.333,0.5,0.067,0.067,0.2,0.333,0.333,0.333,0.067,0.067,0.2,0.333,0.333,0.667,0.333,0.067,0.333,0.133,0.667,0.6,0.333,0.067,0.333,0.267,0.333,0.6,0.1,0.1,0.5,0.1,0.333,0.25,0.2,0.1,0.5,0.2,0.667,0.25],[0.2,0.096,0.134,0.134,0.087,0.134,0.8,0.16,0.8,0.8,null,1,0.8,0.16,0.8,0.8,null,1,1,0.6,1,0.6,1,1,1,0.6,1,0.6,1,1,0.402,0.402,0.8,0.999,0.999,0.999,0.402,0.402,0.8,0.999,0.999,0.999,1,0.402,1,0.665,1,1,1,0.402,1,1,1,1,0.6,0.6,0.75,0.6,0.75,0.75,1,0.6,1,1,1,1],["ns","ns","ns","ns","ns","ns","ns","ns","ns","ns","","ns","ns","ns","ns","ns","","ns","ns","ns","ns","ns","ns","ns","ns","ns","ns","ns","ns","ns","ns","ns","ns","ns","ns","ns","ns","ns","ns","ns","ns","ns","ns","ns","ns","ns","ns","ns","ns","ns","ns","ns","ns","ns","ns","ns","ns","ns","ns","ns","ns","ns","ns","ns","ns","ns"]],"container":"<table class=\"display\">\n  <thead>\n    <tr>\n      <th> <\/th>\n      <th>Reactor_Treatment_Dose<\/th>\n      <th>Distance<\/th>\n      <th>X1<\/th>\n      <th>X2<\/th>\n      <th>R2<\/th>\n      <th>pval<\/th>\n      <th>pvalFDR<\/th>\n      <th>pvalFDR.signif<\/th>\n    <\/tr>\n  <\/thead>\n<\/table>","options":{"columnDefs":[{"className":"dt-right","targets":[5,6,7]},{"orderable":false,"targets":0}],"order":[],"autoWidth":false,"orderClasses":false}},"evals":[],"jsHooks":[]}</script>
```

```r
tmp_df %>% 
  xlsx::write.xlsx(file = out_xlsx,
                   append = TRUE,
                   sheetName = "human1_F1-beta-between-periods",
                   showNA = TRUE)
```



```r
# 
# for (grp in names(lis)){
#   print(grp)
#   
# lapply(
#   dist,
#   FUN = phyloseq_TW,
#   physeq =  ps_tmp_forbdiv %>%
#       # subset_samples(.,  eval(as.name(ranks)) %in% tax) %>% 
#       subset_samples(., Reactor_Treatment_Dose == grp),
#     variable = "period",
#   nrep = 999,
#   strata = NULL
# ) %>%
#   bind_rows(.id = "Distance") %>%
#   mutate_if(is.numeric, round, 5) -> lis[[grp]]
# }
# 
# 
# lis %>% 
#   bind_rows(.id = "Reactor_Treatment_Dose") %>% 
#   xlsx::write.xlsx(file = out_xlsx,
#                    append = TRUE,
#                    sheetName = "human1_F1-beta-between-periodsTW",
#                    showNA = TRUE)
```


### Human1 - F2:


```r
ps_rare %>%
  subset_samples(Model2 == "Human1" & Treatment != "DONOR" & Fermentation == 2) -> ps_tmp_forbdiv
```

#### NMDS:


```r
ps_tmp_forbdiv %>%
  # microbiome::transform(transform = "compositional") %>% 
  phyloseq_plot_bdiv(dlist = dist,
                     m = "PCoA",
                     seed = 123,
                     axis1 = 1,
                     axis2 = 2) -> NMDS_tmp
```

```
## [1] "bjaccard"
## [1] "aichinson"
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

![](NRP72_draft_16S_Figs_beta_PCoA_files/figure-html/unnamed-chunk-47-1.png)<!-- -->


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
<div class="datatables html-widget html-fill-item-overflow-hidden html-fill-item" id="htmlwidget-2f4c52b7f36e26c300a7" style="width:100%;height:auto;"></div>
<script type="application/json" data-for="htmlwidget-2f4c52b7f36e26c300a7">{"x":{"filter":"none","vertical":false,"data":[["ASV0012","ASV0004","ASV0010","ASV0052","ASV0035","ASV0059","ASV0106","ASV0275"],[0.445553060796587,-0.527431856534814,-0.682394062335409,0.675926734603464,-0.442637448488747,-0.19625688252075,-0.263526242835191,0.71119451169905],[0.343116359947506,0.288169711770378,-0.248021115693937,0.0173998730094112,0.334570499832224,0.510439207106555,-0.620428762826066,0.026637446057722],[0.001,0.002,0.001,0.001,0.002,0.002,0.001,0.001],[0.316246366448833,0.361226146069584,0.527176130140688,0.457179706132445,0.307865330162612,0.299064948088332,0.454377930404715,0.506507187003329],["Bacteria","Bacteria","Bacteria","Bacteria","Bacteria","Bacteria","Bacteria","Bacteria"],["Firmicutes","Proteobacteria","Firmicutes","Firmicutes","Proteobacteria","Bacteroidota","Firmicutes","Firmicutes"],["Clostridia","Gammaproteobacteria","Clostridia","Negativicutes","Gammaproteobacteria","Bacteroidia","Clostridia","Bacilli"],["Oscillospirales","Enterobacterales","Lachnospirales","Veillonellales-Selenomonadales","Burkholderiales","Bacteroidales","Oscillospirales","Erysipelotrichales"],["Ruminococcaceae","Enterobacteriaceae","Lachnospiraceae","Veillonellaceae","Sutterellaceae","Tannerellaceae","Oscillospiraceae","Erysipelotrichaceae"],[null,null,null,null,null,null,null,null],[null,null,null,null,null,null,null,null],[null,null,null,null,null,null,null,null],[0.00742857142857143,0.00945454545454546,0.00742857142857143,0.00742857142857143,0.00945454545454546,0.00945454545454546,0.00742857142857143,0.00742857142857143]],"container":"<table class=\"display\">\n  <thead>\n    <tr>\n      <th> <\/th>\n      <th>Axis.1<\/th>\n      <th>Axis.2<\/th>\n      <th>pval<\/th>\n      <th>r<\/th>\n      <th>Kingdom<\/th>\n      <th>Phylum<\/th>\n      <th>Class<\/th>\n      <th>Order<\/th>\n      <th>tax_rank_plot<\/th>\n      <th>Genus<\/th>\n      <th>Species<\/th>\n      <th>Strain<\/th>\n      <th>pval.adj<\/th>\n    <\/tr>\n  <\/thead>\n<\/table>","options":{"columnDefs":[{"className":"dt-right","targets":[1,2,3,4,13]},{"orderable":false,"targets":0}],"order":[],"autoWidth":false,"orderClasses":false}},"evals":[],"jsHooks":[]}</script>
```


```r
out$plot -> NMDS_tmp2_fam

NMDS_tmp2 + 
  theme(legend.position = "none") -> NMDS_tmp2_noleg

NMDS_tmp2_noleg
```

![](NRP72_draft_16S_Figs_beta_PCoA_files/figure-html/unnamed-chunk-49-1.png)<!-- -->

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
<div class="datatables html-widget html-fill-item-overflow-hidden html-fill-item" id="htmlwidget-1d5df3d80f5735caaf13" style="width:100%;height:auto;"></div>
<script type="application/json" data-for="htmlwidget-1d5df3d80f5735caaf13">{"x":{"filter":"none","vertical":false,"data":[["1","2","3","4"],["Formiat_mM","Acetat_mM","Valerat_mM","Total_SCFA_mM"],[-0.613360599591134,-0.127975964263856,0.178226366838872,-0.199371069496966],[0.042165219774713,-0.372980394645997,-0.677351910348987,-0.424683773969117],["Formiat_mM","Acetat_mM","Valerat_mM","Total_SCFA_mM"],[0.001,0.04,0.001,0.007],[0.377989130889445,0.155492222219548,0.490570248290006,0.220105131225016]],"container":"<table class=\"display\">\n  <thead>\n    <tr>\n      <th> <\/th>\n      <th>label<\/th>\n      <th>Axis.1<\/th>\n      <th>Axis.2<\/th>\n      <th>id<\/th>\n      <th>pval<\/th>\n      <th>r<\/th>\n    <\/tr>\n  <\/thead>\n<\/table>","options":{"columnDefs":[{"className":"dt-right","targets":[2,3,5,6]},{"orderable":false,"targets":0}],"order":[],"autoWidth":false,"orderClasses":false}},"evals":[],"jsHooks":[]}</script>
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

![](NRP72_draft_16S_Figs_beta_PCoA_files/figure-html/unnamed-chunk-52-1.png)<!-- -->


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
<div class="datatables html-widget html-fill-item-overflow-hidden html-fill-item" id="htmlwidget-67c25ca3336b15b0293a" style="width:100%;height:auto;"></div>
<script type="application/json" data-for="htmlwidget-67c25ca3336b15b0293a">{"x":{"filter":"none","vertical":false,"data":[["ASV0012","ASV0004","ASV0010","ASV0052","ASV0030","ASV0035","ASV0106","ASV0275"],[0.419469264716501,-0.542610944831473,-0.668057142724215,0.634491810605059,0.0923438525333804,-0.475786979717028,-0.235980253754372,0.716607484689898],[0.41313443709802,0.20743545326037,-0.213247478016456,-0.157224913125192,0.549578877374205,0.217200882303641,-0.59337501428833,-0.0384851666690667],[0.001,0.003,0.001,0.001,0.011,0.003,0.001,0.001],[0.3466345271581,0.337456104720238,0.491774832825221,0.427299531032111,0.310564329556598,0.273549473341732,0.407780587743654,0.515007395167129],["Bacteria","Bacteria","Bacteria","Bacteria","Bacteria","Bacteria","Bacteria","Bacteria"],["Firmicutes","Proteobacteria","Firmicutes","Firmicutes","Bacteroidota","Proteobacteria","Firmicutes","Firmicutes"],["Clostridia","Gammaproteobacteria","Clostridia","Negativicutes","Bacteroidia","Gammaproteobacteria","Clostridia","Bacilli"],["Oscillospirales","Enterobacterales","Lachnospirales","Veillonellales-Selenomonadales","Bacteroidales","Burkholderiales","Oscillospirales","Erysipelotrichales"],["Ruminococcaceae","Enterobacteriaceae","Lachnospiraceae","Veillonellaceae","Prevotellaceae","Sutterellaceae","Oscillospiraceae","Erysipelotrichaceae"],[null,null,null,null,null,null,null,null],[null,null,null,null,null,null,null,null],[null,null,null,null,null,null,null,null],[0.0104,0.0222857142857143,0.0104,0.0104,0.0476666666666667,0.0222857142857143,0.0104,0.0104]],"container":"<table class=\"display\">\n  <thead>\n    <tr>\n      <th> <\/th>\n      <th>Axis.1<\/th>\n      <th>Axis.2<\/th>\n      <th>pval<\/th>\n      <th>r<\/th>\n      <th>Kingdom<\/th>\n      <th>Phylum<\/th>\n      <th>Class<\/th>\n      <th>Order<\/th>\n      <th>tax_rank_plot<\/th>\n      <th>Genus<\/th>\n      <th>Species<\/th>\n      <th>Strain<\/th>\n      <th>pval.adj<\/th>\n    <\/tr>\n  <\/thead>\n<\/table>","options":{"columnDefs":[{"className":"dt-right","targets":[1,2,3,4,13]},{"orderable":false,"targets":0}],"order":[],"autoWidth":false,"orderClasses":false}},"evals":[],"jsHooks":[]}</script>
```


```r
out$plot -> NMDS_tmp2_fam

NMDS_tmp2 + 
  theme(legend.position = "none") -> NMDS_tmp2_noleg

NMDS_tmp2_noleg
```

![](NRP72_draft_16S_Figs_beta_PCoA_files/figure-html/unnamed-chunk-54-1.png)<!-- -->

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
  subset_samples(Model2 == "Human1" & Treatment != "DONOR" & Fermentation == 1) %>%
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
<div class="datatables html-widget html-fill-item-overflow-hidden html-fill-item" id="htmlwidget-530a78b6872d546cc16a" style="width:100%;height:auto;"></div>
<script type="application/json" data-for="htmlwidget-530a78b6872d546cc16a">{"x":{"filter":"none","vertical":false,"data":[["1","2","3","4","5"],["Succinat_mM","Formiat_mM","Acetat_mM","Butyrat_mM","Total_SCFA_mM"],[0.164311218332464,0.441439686730363,-0.797826548818299,-0.751687244981628,-0.711900279912203],[0.506523475467576,-0.0242931772264333,-0.468604143606486,-0.51913836034085,-0.498311404715984],["Succinat_mM","Formiat_mM","Acetat_mM","Butyrat_mM","Total_SCFA_mM"],[0.001,0.004,0.001,0.001,0.001],[0.283564207669651,0.195459155480356,0.856117045404486,0.834538351445456,0.75511626460909]],"container":"<table class=\"display\">\n  <thead>\n    <tr>\n      <th> <\/th>\n      <th>label<\/th>\n      <th>Axis.1<\/th>\n      <th>Axis.2<\/th>\n      <th>id<\/th>\n      <th>pval<\/th>\n      <th>r<\/th>\n    <\/tr>\n  <\/thead>\n<\/table>","options":{"columnDefs":[{"className":"dt-right","targets":[2,3,5,6]},{"orderable":false,"targets":0}],"order":[],"autoWidth":false,"orderClasses":false}},"evals":[],"jsHooks":[]}</script>
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

![](NRP72_draft_16S_Figs_beta_PCoA_files/figure-html/unnamed-chunk-58-1.png)<!-- -->


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
<div class="datatables html-widget html-fill-item-overflow-hidden html-fill-item" id="htmlwidget-4cadc5f8c7ab9ff923ac" style="width:100%;height:auto;"></div>
<script type="application/json" data-for="htmlwidget-4cadc5f8c7ab9ff923ac">{"x":{"filter":"none","vertical":false,"data":[["ASV0012","ASV0004","ASV0010","ASV0052","ASV0035","ASV0059","ASV0106","ASV0275"],[0.445553060796587,-0.527431856534814,-0.682394062335409,0.675926734603464,-0.442637448488747,-0.19625688252075,-0.263526242835191,0.71119451169905],[0.343116359947506,0.288169711770378,-0.248021115693937,0.0173998730094112,0.334570499832224,0.510439207106555,-0.620428762826066,0.026637446057722],[0.001,0.002,0.001,0.001,0.002,0.002,0.001,0.001],[0.316246366448833,0.361226146069584,0.527176130140688,0.457179706132445,0.307865330162612,0.299064948088332,0.454377930404715,0.506507187003329],["Bacteria","Bacteria","Bacteria","Bacteria","Bacteria","Bacteria","Bacteria","Bacteria"],["Firmicutes","Proteobacteria","Firmicutes","Firmicutes","Proteobacteria","Bacteroidota","Firmicutes","Firmicutes"],["Clostridia","Gammaproteobacteria","Clostridia","Negativicutes","Gammaproteobacteria","Bacteroidia","Clostridia","Bacilli"],["Oscillospirales","Enterobacterales","Lachnospirales","Veillonellales-Selenomonadales","Burkholderiales","Bacteroidales","Oscillospirales","Erysipelotrichales"],["Ruminococcaceae","Enterobacteriaceae","Lachnospiraceae","Veillonellaceae","Sutterellaceae","Tannerellaceae","Oscillospiraceae","Erysipelotrichaceae"],[null,null,null,null,null,null,null,null],[null,null,null,null,null,null,null,null],[null,null,null,null,null,null,null,null],[0.00742857142857143,0.00945454545454546,0.00742857142857143,0.00742857142857143,0.00945454545454546,0.00945454545454546,0.00742857142857143,0.00742857142857143]],"container":"<table class=\"display\">\n  <thead>\n    <tr>\n      <th> <\/th>\n      <th>Axis.1<\/th>\n      <th>Axis.2<\/th>\n      <th>pval<\/th>\n      <th>r<\/th>\n      <th>Kingdom<\/th>\n      <th>Phylum<\/th>\n      <th>Class<\/th>\n      <th>Order<\/th>\n      <th>tax_rank_plot<\/th>\n      <th>Genus<\/th>\n      <th>Species<\/th>\n      <th>Strain<\/th>\n      <th>pval.adj<\/th>\n    <\/tr>\n  <\/thead>\n<\/table>","options":{"columnDefs":[{"className":"dt-right","targets":[1,2,3,4,13]},{"orderable":false,"targets":0}],"order":[],"autoWidth":false,"orderClasses":false}},"evals":[],"jsHooks":[]}</script>
```


```r
out$plot -> PCoA_tmp2_fam

PCoA_tmp2 + 
  theme(legend.position = "none") -> PCoA_tmp2_noleg

PCoA_tmp2_noleg
```

![](NRP72_draft_16S_Figs_beta_PCoA_files/figure-html/unnamed-chunk-60-1.png)<!-- -->

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
<div class="datatables html-widget html-fill-item-overflow-hidden html-fill-item" id="htmlwidget-fc8df40fd67c2027a709" style="width:100%;height:auto;"></div>
<script type="application/json" data-for="htmlwidget-fc8df40fd67c2027a709">{"x":{"filter":"none","vertical":false,"data":[["1","2","3","4"],["Formiat_mM","Acetat_mM","Valerat_mM","Total_SCFA_mM"],[-0.613360599591134,-0.127975964263856,0.178226366838872,-0.199371069496966],[0.042165219774713,-0.372980394645997,-0.677351910348987,-0.424683773969117],["Formiat_mM","Acetat_mM","Valerat_mM","Total_SCFA_mM"],[0.001,0.04,0.001,0.007],[0.377989130889445,0.155492222219548,0.490570248290006,0.220105131225016]],"container":"<table class=\"display\">\n  <thead>\n    <tr>\n      <th> <\/th>\n      <th>label<\/th>\n      <th>Axis.1<\/th>\n      <th>Axis.2<\/th>\n      <th>id<\/th>\n      <th>pval<\/th>\n      <th>r<\/th>\n    <\/tr>\n  <\/thead>\n<\/table>","options":{"columnDefs":[{"className":"dt-right","targets":[2,3,5,6]},{"orderable":false,"targets":0}],"order":[],"autoWidth":false,"orderClasses":false}},"evals":[],"jsHooks":[]}</script>
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

![](NRP72_draft_16S_Figs_beta_PCoA_files/figure-html/unnamed-chunk-63-1.png)<!-- -->


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
<div class="datatables html-widget html-fill-item-overflow-hidden html-fill-item" id="htmlwidget-d8c7d17cf4080686b93c" style="width:100%;height:auto;"></div>
<script type="application/json" data-for="htmlwidget-d8c7d17cf4080686b93c">{"x":{"filter":"none","vertical":false,"data":[["ASV0012","ASV0004","ASV0010","ASV0052","ASV0030","ASV0035","ASV0106","ASV0275"],[0.419469264716501,-0.542610944831473,-0.668057142724215,0.634491810605059,0.0923438525333804,-0.475786979717028,-0.235980253754372,0.716607484689898],[0.41313443709802,0.20743545326037,-0.213247478016456,-0.157224913125192,0.549578877374205,0.217200882303641,-0.59337501428833,-0.0384851666690667],[0.001,0.003,0.001,0.001,0.011,0.003,0.001,0.001],[0.3466345271581,0.337456104720238,0.491774832825221,0.427299531032111,0.310564329556598,0.273549473341732,0.407780587743654,0.515007395167129],["Bacteria","Bacteria","Bacteria","Bacteria","Bacteria","Bacteria","Bacteria","Bacteria"],["Firmicutes","Proteobacteria","Firmicutes","Firmicutes","Bacteroidota","Proteobacteria","Firmicutes","Firmicutes"],["Clostridia","Gammaproteobacteria","Clostridia","Negativicutes","Bacteroidia","Gammaproteobacteria","Clostridia","Bacilli"],["Oscillospirales","Enterobacterales","Lachnospirales","Veillonellales-Selenomonadales","Bacteroidales","Burkholderiales","Oscillospirales","Erysipelotrichales"],["Ruminococcaceae","Enterobacteriaceae","Lachnospiraceae","Veillonellaceae","Prevotellaceae","Sutterellaceae","Oscillospiraceae","Erysipelotrichaceae"],[null,null,null,null,null,null,null,null],[null,null,null,null,null,null,null,null],[null,null,null,null,null,null,null,null],[0.0104,0.0222857142857143,0.0104,0.0104,0.0476666666666667,0.0222857142857143,0.0104,0.0104]],"container":"<table class=\"display\">\n  <thead>\n    <tr>\n      <th> <\/th>\n      <th>Axis.1<\/th>\n      <th>Axis.2<\/th>\n      <th>pval<\/th>\n      <th>r<\/th>\n      <th>Kingdom<\/th>\n      <th>Phylum<\/th>\n      <th>Class<\/th>\n      <th>Order<\/th>\n      <th>tax_rank_plot<\/th>\n      <th>Genus<\/th>\n      <th>Species<\/th>\n      <th>Strain<\/th>\n      <th>pval.adj<\/th>\n    <\/tr>\n  <\/thead>\n<\/table>","options":{"columnDefs":[{"className":"dt-right","targets":[1,2,3,4,13]},{"orderable":false,"targets":0}],"order":[],"autoWidth":false,"orderClasses":false}},"evals":[],"jsHooks":[]}</script>
```


```r
out$plot -> PCoA_tmp2_fam

PCoA_tmp2 + 
  theme(legend.position = "none") -> PCoA_tmp2_noleg

PCoA_tmp2_noleg
```

![](NRP72_draft_16S_Figs_beta_PCoA_files/figure-html/unnamed-chunk-65-1.png)<!-- -->

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
<div class="datatables html-widget html-fill-item-overflow-hidden html-fill-item" id="htmlwidget-fbc9973e359619ea20b8" style="width:100%;height:auto;"></div>
<script type="application/json" data-for="htmlwidget-fbc9973e359619ea20b8">{"x":{"filter":"none","vertical":false,"data":[["1","2","3"],["Formiat_mM","Valerat_mM","Total_SCFA_mM"],[-0.600973704253337,0.185820508581589,-0.18149770546043],[0.196394594911277,-0.579040546718349,-0.332598322258841],["Formiat_mM","Valerat_mM","Total_SCFA_mM"],[0.001,0.001,0.04],[0.399740230114341,0.369817216153405,0.143563061056797]],"container":"<table class=\"display\">\n  <thead>\n    <tr>\n      <th> <\/th>\n      <th>label<\/th>\n      <th>Axis.1<\/th>\n      <th>Axis.2<\/th>\n      <th>id<\/th>\n      <th>pval<\/th>\n      <th>r<\/th>\n    <\/tr>\n  <\/thead>\n<\/table>","options":{"columnDefs":[{"className":"dt-right","targets":[2,3,5,6]},{"orderable":false,"targets":0}],"order":[],"autoWidth":false,"orderClasses":false}},"evals":[],"jsHooks":[]}</script>
```


```r
env_out$plot %>% 
  export::graph2ppt(append = TRUE,
                    width = 317.48031496 * 1,
                    height = 0.618 * 317.48031496 * 1 , paper = "A4",  scaling = 2,
                    file = out_pptx)
```

##### pret - Day_of_Treatment <= 0:


```r
lapply(
  dist,
  FUN = phyloseq_adonis2,
  physeq = ps_tmp_forbdiv %>%
    subset_samples(Day_of_Treatment <= 0),
  formula = paste0(c("Reactor_Treatment_Dose"), collapse=" + "),
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
<div class="datatables html-widget html-fill-item-overflow-hidden html-fill-item" id="htmlwidget-fc9c72093435cc6a0a18" style="width:100%;height:auto;"></div>
<script type="application/json" data-for="htmlwidget-fc9c72093435cc6a0a18">{"x":{"filter":"none","vertical":false,"data":[["1","2","3","4","5","6"],["bjaccard","bjaccard","bjaccard","aichinson","aichinson","aichinson"],["Reactor_Treatment_Dose","Residual","Total","Reactor_Treatment_Dose","Residual","Total"],[5,11,16,5,11,16],[0.317,0.288,0.605,0.458,0.404,0.862],[0.524,0.476,1,0.531,0.469,1],[2.418,null,null,2.491,null,null],[0.001,null,null,0.001,null,null]],"container":"<table class=\"display\">\n  <thead>\n    <tr>\n      <th> <\/th>\n      <th>Distance<\/th>\n      <th>terms<\/th>\n      <th>Df<\/th>\n      <th>SumOfSqs<\/th>\n      <th>R2<\/th>\n      <th>F<\/th>\n      <th>Pr(&gt;F)<\/th>\n    <\/tr>\n  <\/thead>\n<\/table>","options":{"columnDefs":[{"className":"dt-right","targets":[3,4,5,6,7]},{"orderable":false,"targets":0}],"order":[],"autoWidth":false,"orderClasses":false}},"evals":[],"jsHooks":[]}</script>
```

```r
tmp_df %>% 
  xlsx::write.xlsx(file = out_xlsx,
                   append = TRUE,
                   sheetName = "human1_F2-beta-pretreat",
                   showNA = TRUE)
```



```r
lapply(
  dist,
  FUN = phyloseq_adonis2,
  physeq = ps_tmp_forbdiv %>%
    subset_samples(Day_of_Treatment <= 0),
  formula = paste0(c("Reactor_Treatment_Dose", "Day_of_Treatment"), collapse=" + "),
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
<div class="datatables html-widget html-fill-item-overflow-hidden html-fill-item" id="htmlwidget-33198698d6599fd3fc97" style="width:100%;height:auto;"></div>
<script type="application/json" data-for="htmlwidget-33198698d6599fd3fc97">{"x":{"filter":"none","vertical":false,"data":[["1","2","3","4","5","6","7","8"],["bjaccard","bjaccard","bjaccard","bjaccard","aichinson","aichinson","aichinson","aichinson"],["Reactor_Treatment_Dose","Day_of_Treatment","Residual","Total","Reactor_Treatment_Dose","Day_of_Treatment","Residual","Total"],[5,1,10,16,5,1,10,16],[0.298,0.036,0.252,0.605,0.434,0.043,0.361,0.862],[0.493,0.06,0.416,1,0.504,0.05,0.419,1],[2.369,1.444,null,null,2.408,1.204,null,null],[0.001,0.149,null,null,0.001,0.251,null,null]],"container":"<table class=\"display\">\n  <thead>\n    <tr>\n      <th> <\/th>\n      <th>Distance<\/th>\n      <th>terms<\/th>\n      <th>Df<\/th>\n      <th>SumOfSqs<\/th>\n      <th>R2<\/th>\n      <th>F<\/th>\n      <th>Pr(&gt;F)<\/th>\n    <\/tr>\n  <\/thead>\n<\/table>","options":{"columnDefs":[{"className":"dt-right","targets":[3,4,5,6,7]},{"orderable":false,"targets":0}],"order":[],"autoWidth":false,"orderClasses":false}},"evals":[],"jsHooks":[]}</script>
```

```r
tmp_df %>% 
  xlsx::write.xlsx(file = out_xlsx,
                   append = TRUE,
                   sheetName = "human1_F2-beta-pretreat-time-day",
                   showNA = TRUE)
```


```r
lapply(
  dist,
  FUN = physeq_pairwise_permanovas_adonis2,
  physeq = ps_tmp_forbdiv %>%
    subset_samples(Day_of_Treatment <= 0),
  compare_header = "Reactor_Treatment_Dose",
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
<div class="datatables html-widget html-fill-item-overflow-hidden html-fill-item" id="htmlwidget-38c86864160d6d855a80" style="width:100%;height:auto;"></div>
<script type="application/json" data-for="htmlwidget-38c86864160d6d855a80">{"x":{"filter":"none","vertical":false,"data":[["1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","23","24","25","26","27","28","29","30"],["bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson"],["CR_UNTREATED","CR_UNTREATED","CR_UNTREATED","CR_UNTREATED","CR_UNTREATED","TR2_CTX20","TR2_CTX20","TR2_CTX20","TR2_CTX20","TR4_CTX200","TR4_CTX200","TR4_CTX200","TR1_CTX+HV292.120","TR1_CTX+HV292.120","TR3_CTX+HV292.1200","CR_UNTREATED","CR_UNTREATED","CR_UNTREATED","CR_UNTREATED","CR_UNTREATED","TR2_CTX20","TR2_CTX20","TR2_CTX20","TR2_CTX20","TR4_CTX200","TR4_CTX200","TR4_CTX200","TR1_CTX+HV292.120","TR1_CTX+HV292.120","TR3_CTX+HV292.1200"],["TR2_CTX20","TR4_CTX200","TR1_CTX+HV292.120","TR3_CTX+HV292.1200","TR5_HV292.1","TR4_CTX200","TR1_CTX+HV292.120","TR3_CTX+HV292.1200","TR5_HV292.1","TR1_CTX+HV292.120","TR3_CTX+HV292.1200","TR5_HV292.1","TR3_CTX+HV292.1200","TR5_HV292.1","TR5_HV292.1","TR2_CTX20","TR4_CTX200","TR1_CTX+HV292.120","TR3_CTX+HV292.1200","TR5_HV292.1","TR4_CTX200","TR1_CTX+HV292.120","TR3_CTX+HV292.1200","TR5_HV292.1","TR1_CTX+HV292.120","TR3_CTX+HV292.1200","TR5_HV292.1","TR3_CTX+HV292.1200","TR5_HV292.1","TR5_HV292.1"],[0.373,0.425,0.211,0.382,0.386,0.45,0.302,0.358,0.501,0.397,0.454,0.543,0.341,0.438,0.563,0.433,0.443,0.257,0.473,0.428,0.44,0.309,0.279,0.455,0.372,0.455,0.531,0.324,0.392,0.501],[0.1,0.1,0.4,0.1,0.1,0.1,0.2,0.2,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.2,0.2,0.1,0.1,0.1,0.1,0.2,0.1,0.1,0.1,0.1,0.1,0.1,0.1],[0.231,0.231,0.4,0.231,0.231,0.231,0.222,0.222,0.231,0.231,0.231,0.231,0.231,0.231,0.231,0.231,0.214,0.214,0.231,0.231,0.231,0.231,0.214,0.231,0.231,0.231,0.231,0.231,0.231,0.231]],"container":"<table class=\"display\">\n  <thead>\n    <tr>\n      <th> <\/th>\n      <th>Distance<\/th>\n      <th>X1<\/th>\n      <th>X2<\/th>\n      <th>R2<\/th>\n      <th>pval<\/th>\n      <th>pvalFDR<\/th>\n    <\/tr>\n  <\/thead>\n<\/table>","options":{"columnDefs":[{"className":"dt-right","targets":[4,5,6]},{"orderable":false,"targets":0}],"order":[],"autoWidth":false,"orderClasses":false}},"evals":[],"jsHooks":[]}</script>
```

```r
tmp_df %>% 
  xlsx::write.xlsx(file = out_xlsx,
                   append = TRUE,
                   sheetName = "human1_F2-beta-pretreat-pw",
                   showNA = TRUE)
```


```r
lapply(
  dist,
  FUN = phyloseq_TW,
  physeq =  ps_tmp_forbdiv %>%
    subset_samples(Day_of_Treatment <= 0),
  variable = "Reactor_Treatment_Dose",
  nrep = 999,
  strata = NULL
) %>%
  bind_rows(.id = "Distance") %>%
  mutate_if(is.numeric, round, 3)   -> tmp_df

tmp_df %>% 
  DT::datatable()
```

```{=html}
<div class="datatables html-widget html-fill-item-overflow-hidden html-fill-item" id="htmlwidget-d96970b9501c191ad406" style="width:100%;height:auto;"></div>
<script type="application/json" data-for="htmlwidget-d96970b9501c191ad406">{"x":{"filter":"none","vertical":false,"data":[["1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","23","24","25","26","27","28","29","30"],["bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson"],[["CR_UNTREATED"],["CR_UNTREATED"],["CR_UNTREATED"],["CR_UNTREATED"],["CR_UNTREATED"],["TR2_CTX20"],["TR2_CTX20"],["TR2_CTX20"],["TR2_CTX20"],["TR4_CTX200"],["TR4_CTX200"],["TR4_CTX200"],["TR1_CTX+HV292.120"],["TR1_CTX+HV292.120"],["TR3_CTX+HV292.1200"],["CR_UNTREATED"],["CR_UNTREATED"],["CR_UNTREATED"],["CR_UNTREATED"],["CR_UNTREATED"],["TR2_CTX20"],["TR2_CTX20"],["TR2_CTX20"],["TR2_CTX20"],["TR4_CTX200"],["TR4_CTX200"],["TR4_CTX200"],["TR1_CTX+HV292.120"],["TR1_CTX+HV292.120"],["TR3_CTX+HV292.1200"]],[["TR2_CTX20"],["TR4_CTX200"],["TR1_CTX+HV292.120"],["TR3_CTX+HV292.1200"],["TR5_HV292.1"],["TR4_CTX200"],["TR1_CTX+HV292.120"],["TR3_CTX+HV292.1200"],["TR5_HV292.1"],["TR1_CTX+HV292.120"],["TR3_CTX+HV292.1200"],["TR5_HV292.1"],["TR3_CTX+HV292.1200"],["TR5_HV292.1"],["TR5_HV292.1"],["TR2_CTX20"],["TR4_CTX200"],["TR1_CTX+HV292.120"],["TR3_CTX+HV292.1200"],["TR5_HV292.1"],["TR4_CTX200"],["TR1_CTX+HV292.120"],["TR3_CTX+HV292.1200"],["TR5_HV292.1"],["TR1_CTX+HV292.120"],["TR3_CTX+HV292.1200"],["TR5_HV292.1"],["TR3_CTX+HV292.1200"],["TR5_HV292.1"],["TR5_HV292.1"]],[[3],[3],[3],[3],[3],[3],[3],[3],[3],[2],[2],[2],[3],[3],[3],[3],[3],[3],[3],[3],[3],[3],[3],[3],[2],[2],[2],[3],[3],[3]],[[3],[2],[3],[3],[3],[2],[3],[3],[3],[3],[3],[3],[3],[3],[3],[3],[2],[3],[3],[3],[2],[3],[3],[3],[3],[3],[3],[3],[3],[3]],[[0.106],[0.099],[0.4],[0.116],[0.098],[0.19],[0.201],[0.189],[0.091],[0.193],[0.112],[0.094],[0.089],[0.092],[0.1],[0.088],[0.197],[0.204],[0.097],[0.096],[0.188],[0.108],[0.222],[0.091],[0.113],[0.096],[0.101],[0.107],[0.09],[0.103]],[[2.37848924023348],[1.94287361790263],[1.07212993939922],[2.4673213356877],[2.51415311555688],[1.65381619349415],[1.7324220040701],[2.23340677556841],[4.01443896882939],[1.72879834780556],[1.66624018029317],[2.35895165365943],[2.06589776511485],[3.11260641673556],[5.16060425041003],[3.05031318440669],[2.0794357425907],[1.38463609911914],[3.58884746681679],[2.99061286086113],[1.85098110545744],[1.78725825121598],[1.54769976454627],[3.34468876378856],[1.63137411663642],[1.87262022077958],[2.5447529353778],[1.91653516574998],[2.57510363230682],[4.01967167353687]],[[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999]],[0.174,0.174,0.4,0.174,0.174,0.215,0.215,0.215,0.174,0.215,0.174,0.174,0.174,0.174,0.174,0.154,0.219,0.219,0.154,0.154,0.219,0.154,0.222,0.154,0.154,0.154,0.154,0.154,0.154,0.154]],"container":"<table class=\"display\">\n  <thead>\n    <tr>\n      <th> <\/th>\n      <th>Distance<\/th>\n      <th>Level1<\/th>\n      <th>Level2<\/th>\n      <th>N1<\/th>\n      <th>N2<\/th>\n      <th>p.value<\/th>\n      <th>tw2.stat<\/th>\n      <th>nrep<\/th>\n      <th>pvalFDR<\/th>\n    <\/tr>\n  <\/thead>\n<\/table>","options":{"columnDefs":[{"className":"dt-right","targets":9},{"orderable":false,"targets":0}],"order":[],"autoWidth":false,"orderClasses":false}},"evals":[],"jsHooks":[]}</script>
```

```r
tmp_df %>% 
  xlsx::write.xlsx(file = out_xlsx,
                   append = TRUE,
                   sheetName = "human1_F2-beta-pretreat-pw2",
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

##### treat - Day_of_Treatment > 1:


```r
lapply(
  dist,
  FUN = phyloseq_adonis2,
  physeq = ps_tmp_forbdiv %>%
    subset_samples(Day_of_Treatment >1),
  formula = paste0(c("Reactor_Treatment_Dose"), collapse=" + "),
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
<div class="datatables html-widget html-fill-item-overflow-hidden html-fill-item" id="htmlwidget-bbf49928d385f1dd9e5c" style="width:100%;height:auto;"></div>
<script type="application/json" data-for="htmlwidget-bbf49928d385f1dd9e5c">{"x":{"filter":"none","vertical":false,"data":[["1","2","3","4","5","6"],["bjaccard","bjaccard","bjaccard","aichinson","aichinson","aichinson"],["Reactor_Treatment_Dose","Residual","Total","Reactor_Treatment_Dose","Residual","Total"],[5,17,22,5,17,22],[0.742,1.261,2.004,0.773,1.166,1.939],[0.371,0.629,1,0.399,0.601,1],[2.001,null,null,2.253,null,null],[0.007,null,null,0.001,null,null]],"container":"<table class=\"display\">\n  <thead>\n    <tr>\n      <th> <\/th>\n      <th>Distance<\/th>\n      <th>terms<\/th>\n      <th>Df<\/th>\n      <th>SumOfSqs<\/th>\n      <th>R2<\/th>\n      <th>F<\/th>\n      <th>Pr(&gt;F)<\/th>\n    <\/tr>\n  <\/thead>\n<\/table>","options":{"columnDefs":[{"className":"dt-right","targets":[3,4,5,6,7]},{"orderable":false,"targets":0}],"order":[],"autoWidth":false,"orderClasses":false}},"evals":[],"jsHooks":[]}</script>
```

```r
tmp_df %>% 
  xlsx::write.xlsx(file = out_xlsx,
                   append = TRUE,
                   sheetName = "human1_F2-beta-treat",
                   showNA = TRUE)
```


```r
lapply(
  dist,
  FUN = phyloseq_adonis2,
  physeq = ps_tmp_forbdiv %>%
    subset_samples(Day_of_Treatment >1),
  formula = paste0(c("Reactor_Treatment_Dose", "Day_of_Treatment"), collapse=" + "),
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
<div class="datatables html-widget html-fill-item-overflow-hidden html-fill-item" id="htmlwidget-90bfa43893a2d8eda041" style="width:100%;height:auto;"></div>
<script type="application/json" data-for="htmlwidget-90bfa43893a2d8eda041">{"x":{"filter":"none","vertical":false,"data":[["1","2","3","4","5","6","7","8"],["bjaccard","bjaccard","bjaccard","bjaccard","aichinson","aichinson","aichinson","aichinson"],["Reactor_Treatment_Dose","Day_of_Treatment","Residual","Total","Reactor_Treatment_Dose","Day_of_Treatment","Residual","Total"],[5,1,16,22,5,1,16,22],[0.725,0.158,1.103,2.004,0.752,0.135,1.031,1.939],[0.362,0.079,0.55,1,0.388,0.07,0.532,1],[2.104,2.298,null,null,2.334,2.101,null,null],[0.005,0.032,null,null,0.001,0.033,null,null]],"container":"<table class=\"display\">\n  <thead>\n    <tr>\n      <th> <\/th>\n      <th>Distance<\/th>\n      <th>terms<\/th>\n      <th>Df<\/th>\n      <th>SumOfSqs<\/th>\n      <th>R2<\/th>\n      <th>F<\/th>\n      <th>Pr(&gt;F)<\/th>\n    <\/tr>\n  <\/thead>\n<\/table>","options":{"columnDefs":[{"className":"dt-right","targets":[3,4,5,6,7]},{"orderable":false,"targets":0}],"order":[],"autoWidth":false,"orderClasses":false}},"evals":[],"jsHooks":[]}</script>
```

```r
tmp_df %>% 
  xlsx::write.xlsx(file = out_xlsx,
                   append = TRUE,
                   sheetName = "human1_F2-beta-treat-time-day",
                   showNA = TRUE)
```


```r
lapply(
  dist,
  FUN = physeq_pairwise_permanovas_adonis2,
  physeq = ps_tmp_forbdiv %>%
    subset_samples(Day_of_Treatment > 1),
  compare_header = "Reactor_Treatment_Dose",
  n_perm = 999,
  strat = "none"
) %>%
  bind_rows(.id = "Distance") %>%
  mutate_if(is.numeric, round, 3) %>%
  select(-pvalBon)  -> tmp_df
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

```r
tmp_df %>% 
  DT::datatable()
```

```{=html}
<div class="datatables html-widget html-fill-item-overflow-hidden html-fill-item" id="htmlwidget-0349e08474188736d854" style="width:100%;height:auto;"></div>
<script type="application/json" data-for="htmlwidget-0349e08474188736d854">{"x":{"filter":"none","vertical":false,"data":[["1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","23","24","25","26","27","28","29","30"],["bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson"],["CR_UNTREATED","CR_UNTREATED","CR_UNTREATED","CR_UNTREATED","CR_UNTREATED","TR2_CTX20","TR2_CTX20","TR2_CTX20","TR2_CTX20","TR4_CTX200","TR4_CTX200","TR4_CTX200","TR1_CTX+HV292.120","TR1_CTX+HV292.120","TR3_CTX+HV292.1200","CR_UNTREATED","CR_UNTREATED","CR_UNTREATED","CR_UNTREATED","CR_UNTREATED","TR2_CTX20","TR2_CTX20","TR2_CTX20","TR2_CTX20","TR4_CTX200","TR4_CTX200","TR4_CTX200","TR1_CTX+HV292.120","TR1_CTX+HV292.120","TR3_CTX+HV292.1200"],["TR2_CTX20","TR4_CTX200","TR1_CTX+HV292.120","TR3_CTX+HV292.1200","TR5_HV292.1","TR4_CTX200","TR1_CTX+HV292.120","TR3_CTX+HV292.1200","TR5_HV292.1","TR1_CTX+HV292.120","TR3_CTX+HV292.1200","TR5_HV292.1","TR3_CTX+HV292.1200","TR5_HV292.1","TR5_HV292.1","TR2_CTX20","TR4_CTX200","TR1_CTX+HV292.120","TR3_CTX+HV292.1200","TR5_HV292.1","TR4_CTX200","TR1_CTX+HV292.120","TR3_CTX+HV292.1200","TR5_HV292.1","TR1_CTX+HV292.120","TR3_CTX+HV292.1200","TR5_HV292.1","TR3_CTX+HV292.1200","TR5_HV292.1","TR5_HV292.1"],[0.287,0.374,0.489,0.469,0.202,0.055,0.107,0.157,0.202,0.197,0.133,0.304,0.339,0.673,0.385,0.305,0.343,0.466,0.435,0.251,0.063,0.157,0.176,0.251,0.226,0.136,0.325,0.388,0.669,0.407],[0.006,0.025,0.025,0.008,0.6,0.97,0.508,0.076,0.383,0.424,0.36,0.4,0.023,0.25,0.167,0.005,0.025,0.037,0.009,0.6,0.938,0.204,0.076,0.293,0.236,0.313,0.2,0.011,0.25,0.167],[0.09,0.083,0.083,0.06,0.643,0.97,0.586,0.19,0.575,0.53,0.6,0.545,0.115,0.469,0.357,0.075,0.094,0.111,0.067,0.643,0.938,0.34,0.19,0.366,0.354,0.361,0.375,0.055,0.341,0.357]],"container":"<table class=\"display\">\n  <thead>\n    <tr>\n      <th> <\/th>\n      <th>Distance<\/th>\n      <th>X1<\/th>\n      <th>X2<\/th>\n      <th>R2<\/th>\n      <th>pval<\/th>\n      <th>pvalFDR<\/th>\n    <\/tr>\n  <\/thead>\n<\/table>","options":{"columnDefs":[{"className":"dt-right","targets":[4,5,6]},{"orderable":false,"targets":0}],"order":[],"autoWidth":false,"orderClasses":false}},"evals":[],"jsHooks":[]}</script>
```

```r
tmp_df %>% 
  xlsx::write.xlsx(file = out_xlsx,
                   append = TRUE,
                   sheetName = "human1_F2-beta-treat-pw",
                   showNA = TRUE)
```


```r
lapply(
  dist,
  FUN = phyloseq_TW,
  physeq =  ps_tmp_forbdiv %>%
    subset_samples(Day_of_Treatment >1),
  variable = "Reactor_Treatment_Dose",
  nrep = 999,
  strata = NULL
) %>%
  bind_rows(.id = "Distance") %>%
  mutate_if(is.numeric, round, 3)   -> tmp_df

tmp_df %>% 
  DT::datatable()
```

```{=html}
<div class="datatables html-widget html-fill-item-overflow-hidden html-fill-item" id="htmlwidget-2c91a99e85c6bc68af90" style="width:100%;height:auto;"></div>
<script type="application/json" data-for="htmlwidget-2c91a99e85c6bc68af90">{"x":{"filter":"none","vertical":false,"data":[["1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","23","24","25","26","27","28","29","30"],["bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson"],[["CR_UNTREATED"],["CR_UNTREATED"],["CR_UNTREATED"],["CR_UNTREATED"],["CR_UNTREATED"],["TR2_CTX20"],["TR2_CTX20"],["TR2_CTX20"],["TR2_CTX20"],["TR4_CTX200"],["TR4_CTX200"],["TR4_CTX200"],["TR1_CTX+HV292.120"],["TR1_CTX+HV292.120"],["TR3_CTX+HV292.1200"],["CR_UNTREATED"],["CR_UNTREATED"],["CR_UNTREATED"],["CR_UNTREATED"],["CR_UNTREATED"],["TR2_CTX20"],["TR2_CTX20"],["TR2_CTX20"],["TR2_CTX20"],["TR4_CTX200"],["TR4_CTX200"],["TR4_CTX200"],["TR1_CTX+HV292.120"],["TR1_CTX+HV292.120"],["TR3_CTX+HV292.1200"]],[["TR2_CTX20"],["TR4_CTX200"],["TR1_CTX+HV292.120"],["TR3_CTX+HV292.1200"],["TR5_HV292.1"],["TR4_CTX200"],["TR1_CTX+HV292.120"],["TR3_CTX+HV292.1200"],["TR5_HV292.1"],["TR1_CTX+HV292.120"],["TR3_CTX+HV292.1200"],["TR5_HV292.1"],["TR3_CTX+HV292.1200"],["TR5_HV292.1"],["TR5_HV292.1"],["TR2_CTX20"],["TR4_CTX200"],["TR1_CTX+HV292.120"],["TR3_CTX+HV292.1200"],["TR5_HV292.1"],["TR4_CTX200"],["TR1_CTX+HV292.120"],["TR3_CTX+HV292.1200"],["TR5_HV292.1"],["TR1_CTX+HV292.120"],["TR3_CTX+HV292.1200"],["TR5_HV292.1"],["TR3_CTX+HV292.1200"],["TR5_HV292.1"],["TR5_HV292.1"]],[[4],[4],[4],[4],[4],[6],[6],[6],[6],[4],[4],[4],[3],[3],[5],[4],[4],[4],[4],[4],[6],[6],[6],[6],[4],[4],[4],[3],[3],[5]],[[6],[4],[3],[5],[1],[4],[3],[5],[1],[3],[5],[1],[5],[1],[1],[6],[4],[3],[5],[1],[4],[3],[5],[1],[3],[5],[1],[5],[1],[1]],[[0.004],[0.031],[0.019],[0.01],[null],[0.967],[0.263],[0.062],[null],[0.29],[0.435],[null],[0.017],[null],[null],[0.005],[0.024],[0.032],[0.01],[null],[0.946],[0.092],[0.05],[null],[0.174],[0.375],[null],[0.019],[null],[null]],[[3.95569575816069],[3.5847143872079],[5.00082177595644],[6.41561997709568],[null],[0.471739484141103],[1.26942179750072],[1.78758279459166],[null],[1.47922073693661],[0.989441309836014],[null],[3.59051499230875],[null],[null],[3.84409825373469],[3.13552265685569],[4.69277796183477],[5.35118899528143],[null],[0.523124366218586],[1.76431550865016],[2.01622573575824],[null],[1.69637009804164],[1.02433520209128],[null],[4.24253957811912],[null],[null]],[[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999]],[0.04,0.062,0.048,0.048,null,0.967,0.362,0.103,null,0.362,0.483,null,0.048,null,null,0.05,0.06,0.064,0.05,null,0.946,0.131,0.083,null,0.217,0.417,null,0.06,null,null]],"container":"<table class=\"display\">\n  <thead>\n    <tr>\n      <th> <\/th>\n      <th>Distance<\/th>\n      <th>Level1<\/th>\n      <th>Level2<\/th>\n      <th>N1<\/th>\n      <th>N2<\/th>\n      <th>p.value<\/th>\n      <th>tw2.stat<\/th>\n      <th>nrep<\/th>\n      <th>pvalFDR<\/th>\n    <\/tr>\n  <\/thead>\n<\/table>","options":{"columnDefs":[{"className":"dt-right","targets":9},{"orderable":false,"targets":0}],"order":[],"autoWidth":false,"orderClasses":false}},"evals":[],"jsHooks":[]}</script>
```

```r
tmp_df %>% 
  xlsx::write.xlsx(file = out_xlsx,
                   append = TRUE,
                   sheetName = "human1_F2-beta-treat-pw2",
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

##### Day_of_Treatment >= 3 & Day_of_Treatment <= 6:


```r
lapply(
  dist,
  FUN = phyloseq_adonis2,
  physeq = ps_tmp_forbdiv %>%
    subset_samples(Day_of_Treatment >= 3 & Day_of_Treatment <= 6),
  formula = paste0(c("Reactor_Treatment_Dose"), collapse=" + "),
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
<div class="datatables html-widget html-fill-item-overflow-hidden html-fill-item" id="htmlwidget-f4acdbbe3467343ec6d3" style="width:100%;height:auto;"></div>
<script type="application/json" data-for="htmlwidget-f4acdbbe3467343ec6d3">{"x":{"filter":"none","vertical":false,"data":[["1","2","3","4","5","6"],["bjaccard","bjaccard","bjaccard","aichinson","aichinson","aichinson"],["Reactor_Treatment_Dose","Residual","Total","Reactor_Treatment_Dose","Residual","Total"],[5,9,14,5,9,14],[0.603,0.42,1.023,0.695,0.471,1.165],[0.589,0.411,1,0.596,0.404,1],[2.584,null,null,2.657,null,null],[0.003,null,null,0.001,null,null]],"container":"<table class=\"display\">\n  <thead>\n    <tr>\n      <th> <\/th>\n      <th>Distance<\/th>\n      <th>terms<\/th>\n      <th>Df<\/th>\n      <th>SumOfSqs<\/th>\n      <th>R2<\/th>\n      <th>F<\/th>\n      <th>Pr(&gt;F)<\/th>\n    <\/tr>\n  <\/thead>\n<\/table>","options":{"columnDefs":[{"className":"dt-right","targets":[3,4,5,6,7]},{"orderable":false,"targets":0}],"order":[],"autoWidth":false,"orderClasses":false}},"evals":[],"jsHooks":[]}</script>
```

```r
tmp_df %>% 
  xlsx::write.xlsx(file = out_xlsx,
                   append = TRUE,
                   sheetName = "human1_F2-beta-d3d6",
                   showNA = TRUE)
```


```r
lapply(
  dist,
  FUN = phyloseq_adonis2,
  physeq = ps_tmp_forbdiv %>%
    subset_samples(Day_of_Treatment >= 3 & Day_of_Treatment <= 6),
  formula = paste0(c("Reactor_Treatment_Dose", "Day_of_Treatment"), collapse=" + "),
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
<div class="datatables html-widget html-fill-item-overflow-hidden html-fill-item" id="htmlwidget-d15c90e001b0818c17b0" style="width:100%;height:auto;"></div>
<script type="application/json" data-for="htmlwidget-d15c90e001b0818c17b0">{"x":{"filter":"none","vertical":false,"data":[["1","2","3","4","5","6","7","8"],["bjaccard","bjaccard","bjaccard","bjaccard","aichinson","aichinson","aichinson","aichinson"],["Reactor_Treatment_Dose","Day_of_Treatment","Residual","Total","Reactor_Treatment_Dose","Day_of_Treatment","Residual","Total"],[5,1,8,14,5,1,8,14],[0.592,0.084,0.336,1.023,0.677,0.081,0.39,1.165],[0.578,0.082,0.328,1,0.581,0.069,0.335,1],[2.819,2.01,null,null,2.778,1.655,null,null],[0.001,0.066,null,null,0.001,0.117,null,null]],"container":"<table class=\"display\">\n  <thead>\n    <tr>\n      <th> <\/th>\n      <th>Distance<\/th>\n      <th>terms<\/th>\n      <th>Df<\/th>\n      <th>SumOfSqs<\/th>\n      <th>R2<\/th>\n      <th>F<\/th>\n      <th>Pr(&gt;F)<\/th>\n    <\/tr>\n  <\/thead>\n<\/table>","options":{"columnDefs":[{"className":"dt-right","targets":[3,4,5,6,7]},{"orderable":false,"targets":0}],"order":[],"autoWidth":false,"orderClasses":false}},"evals":[],"jsHooks":[]}</script>
```

```r
tmp_df %>% 
  xlsx::write.xlsx(file = out_xlsx,
                   append = TRUE,
                   sheetName = "human1_F2-beta-d3d6-time-day",
                   showNA = TRUE)
```


```r
lapply(
  dist,
  FUN = physeq_pairwise_permanovas_adonis2,
  physeq = ps_tmp_forbdiv %>%
    subset_samples(Day_of_Treatment >= 3 & Day_of_Treatment <= 6),
  compare_header = "Reactor_Treatment_Dose",
  n_perm = 999,
  strat = "none"
) %>%
  bind_rows(.id = "Distance") %>%
  mutate_if(is.numeric, round, 3) %>%
  select(-pvalBon)  -> tmp_df
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
<div class="datatables html-widget html-fill-item-overflow-hidden html-fill-item" id="htmlwidget-fec2d181a4be4b780ca9" style="width:100%;height:auto;"></div>
<script type="application/json" data-for="htmlwidget-fec2d181a4be4b780ca9">{"x":{"filter":"none","vertical":false,"data":[["1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","23","24","25","26","27","28","29","30"],["bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson"],["CR_UNTREATED","CR_UNTREATED","CR_UNTREATED","CR_UNTREATED","CR_UNTREATED","TR2_CTX20","TR2_CTX20","TR2_CTX20","TR2_CTX20","TR4_CTX200","TR4_CTX200","TR4_CTX200","TR1_CTX+HV292.120","TR1_CTX+HV292.120","TR3_CTX+HV292.1200","CR_UNTREATED","CR_UNTREATED","CR_UNTREATED","CR_UNTREATED","CR_UNTREATED","TR2_CTX20","TR2_CTX20","TR2_CTX20","TR2_CTX20","TR4_CTX200","TR4_CTX200","TR4_CTX200","TR1_CTX+HV292.120","TR1_CTX+HV292.120","TR3_CTX+HV292.1200"],["TR2_CTX20","TR4_CTX200","TR1_CTX+HV292.120","TR3_CTX+HV292.1200","TR5_HV292.1","TR4_CTX200","TR1_CTX+HV292.120","TR3_CTX+HV292.1200","TR5_HV292.1","TR1_CTX+HV292.120","TR3_CTX+HV292.1200","TR5_HV292.1","TR3_CTX+HV292.1200","TR5_HV292.1","TR5_HV292.1","TR2_CTX20","TR4_CTX200","TR1_CTX+HV292.120","TR3_CTX+HV292.1200","TR5_HV292.1","TR4_CTX200","TR1_CTX+HV292.120","TR3_CTX+HV292.1200","TR5_HV292.1","TR1_CTX+HV292.120","TR3_CTX+HV292.1200","TR5_HV292.1","TR3_CTX+HV292.1200","TR5_HV292.1","TR5_HV292.1"],[0.535,0.412,0.514,0.55,0.238,0.267,0.325,0.52,0.67,0.212,0.189,0.468,0.479,0.788,0.584,0.492,0.386,0.499,0.494,0.278,0.311,0.368,0.517,0.638,0.311,0.222,0.509,0.511,0.798,0.597],[0.029,0.1,0.1,0.1,1,0.333,0.067,0.022,0.2,1,0.7,0.667,0.1,0.333,0.25,0.029,0.1,0.1,0.1,0.75,0.2,0.067,0.044,0.2,1,0.8,0.667,0.1,0.333,0.25],[0.218,0.273,0.273,0.273,1.034,0.476,0.333,0.33,0.375,1.034,0.808,0.833,0.273,0.476,0.417,0.435,0.273,0.273,0.273,0.865,0.353,0.333,0.33,0.353,1,0.857,0.833,0.273,0.455,0.375]],"container":"<table class=\"display\">\n  <thead>\n    <tr>\n      <th> <\/th>\n      <th>Distance<\/th>\n      <th>X1<\/th>\n      <th>X2<\/th>\n      <th>R2<\/th>\n      <th>pval<\/th>\n      <th>pvalFDR<\/th>\n    <\/tr>\n  <\/thead>\n<\/table>","options":{"columnDefs":[{"className":"dt-right","targets":[4,5,6]},{"orderable":false,"targets":0}],"order":[],"autoWidth":false,"orderClasses":false}},"evals":[],"jsHooks":[]}</script>
```

```r
tmp_df %>% 
  xlsx::write.xlsx(file = out_xlsx,
                   append = TRUE,
                   sheetName = "human1_F2-beta-d3d6-pw",
                   showNA = TRUE)
```


```r
lapply(
  dist,
  FUN = phyloseq_TW,
  physeq =  ps_tmp_forbdiv %>%
    subset_samples(Day_of_Treatment >= 3 & Day_of_Treatment <= 6),
  variable = "Reactor_Treatment_Dose",
  nrep = 999,
  strata = NULL
) %>%
  bind_rows(.id = "Distance") %>%
  mutate_if(is.numeric, round, 3)   -> tmp_df

tmp_df %>% 
  DT::datatable()
```

```{=html}
<div class="datatables html-widget html-fill-item-overflow-hidden html-fill-item" id="htmlwidget-bcc569e9a5fb6f3f4fe2" style="width:100%;height:auto;"></div>
<script type="application/json" data-for="htmlwidget-bcc569e9a5fb6f3f4fe2">{"x":{"filter":"none","vertical":false,"data":[["1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","23","24","25","26","27","28","29","30"],["bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson"],[["CR_UNTREATED"],["CR_UNTREATED"],["CR_UNTREATED"],["CR_UNTREATED"],["CR_UNTREATED"],["TR2_CTX20"],["TR2_CTX20"],["TR2_CTX20"],["TR2_CTX20"],["TR4_CTX200"],["TR4_CTX200"],["TR4_CTX200"],["TR1_CTX+HV292.120"],["TR1_CTX+HV292.120"],["TR3_CTX+HV292.1200"],["CR_UNTREATED"],["CR_UNTREATED"],["CR_UNTREATED"],["CR_UNTREATED"],["CR_UNTREATED"],["TR2_CTX20"],["TR2_CTX20"],["TR2_CTX20"],["TR2_CTX20"],["TR4_CTX200"],["TR4_CTX200"],["TR4_CTX200"],["TR1_CTX+HV292.120"],["TR1_CTX+HV292.120"],["TR3_CTX+HV292.1200"]],[["TR2_CTX20"],["TR4_CTX200"],["TR1_CTX+HV292.120"],["TR3_CTX+HV292.1200"],["TR5_HV292.1"],["TR4_CTX200"],["TR1_CTX+HV292.120"],["TR3_CTX+HV292.1200"],["TR5_HV292.1"],["TR1_CTX+HV292.120"],["TR3_CTX+HV292.1200"],["TR5_HV292.1"],["TR3_CTX+HV292.1200"],["TR5_HV292.1"],["TR5_HV292.1"],["TR2_CTX20"],["TR4_CTX200"],["TR1_CTX+HV292.120"],["TR3_CTX+HV292.1200"],["TR5_HV292.1"],["TR4_CTX200"],["TR1_CTX+HV292.120"],["TR3_CTX+HV292.1200"],["TR5_HV292.1"],["TR1_CTX+HV292.120"],["TR3_CTX+HV292.1200"],["TR5_HV292.1"],["TR3_CTX+HV292.1200"],["TR5_HV292.1"],["TR5_HV292.1"]],[[3],[3],[3],[3],[3],[4],[4],[4],[4],[2],[2],[2],[2],[2],[3],[3],[3],[3],[3],[3],[4],[4],[4],[4],[2],[2],[2],[2],[2],[3]],[[4],[2],[2],[3],[1],[2],[2],[3],[1],[2],[3],[1],[3],[1],[1],[4],[2],[2],[3],[1],[2],[2],[3],[1],[2],[3],[1],[3],[1],[1]],[[0.032],[0.203],[0.109],[0.087],[null],[0.796],[0.256],[0.029],[null],[1],[1],[null],[0.104],[null],[null],[0.023],[0.213],[0.107],[0.099],[null],[0.412],[0.076],[0.023],[null],[1],[1],[null],[0.103],[null],[null]],[[5.02566313743676],[1.59845857935754],[3.47288496172523],[4.8831569723275],[null],[0.736396214173821],[1.57310508691593],[4.7169398804577],[null],[0.537611500691741],[0.535883343945286],[null],[3.0273654786785],[null],[null],[4.22435266426004],[1.61125828073231],[3.49183953503474],[3.90283279597235],[null],[1.03763851760866],[2.12228248174519],[4.9067895512456],[null],[0.903560520583048],[0.678124087556537],[null],[3.39241553536459],[null],[null]],[[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999]],[0.16,0.338,0.218,0.218,null,0.995,0.366,0.16,null,1,1,null,0.218,null,null,0.115,0.304,0.178,0.178,null,0.515,0.178,0.115,null,1,1,null,0.178,null,null]],"container":"<table class=\"display\">\n  <thead>\n    <tr>\n      <th> <\/th>\n      <th>Distance<\/th>\n      <th>Level1<\/th>\n      <th>Level2<\/th>\n      <th>N1<\/th>\n      <th>N2<\/th>\n      <th>p.value<\/th>\n      <th>tw2.stat<\/th>\n      <th>nrep<\/th>\n      <th>pvalFDR<\/th>\n    <\/tr>\n  <\/thead>\n<\/table>","options":{"columnDefs":[{"className":"dt-right","targets":9},{"orderable":false,"targets":0}],"order":[],"autoWidth":false,"orderClasses":false}},"evals":[],"jsHooks":[]}</script>
```

```r
tmp_df %>% 
  xlsx::write.xlsx(file = out_xlsx,
                   append = TRUE,
                   sheetName = "human1_F2-beta-d3d6-pw2",
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


##### Period:


```r
setNames(vector("list", 
                length(physeq_get_unique(ps_tmp_forbdiv, "Reactor_Treatment_Dose"))), 
         physeq_get_unique(ps_tmp_forbdiv, "Reactor_Treatment_Dose")) -> lis
```



```r
for (grp in names(lis)){
  print(grp)
  
  lapply(
    dist,
    FUN = physeq_pairwise_permanovas_adonis2,
    physeq = ps_tmp_forbdiv %>%
      # subset_samples(.,  eval(as.name(ranks)) %in% tax) %>% 
      subset_samples(., Reactor_Treatment_Dose == grp),
    compare_header = "period",
    n_perm = 999,
    strat = "none"
  ) %>%
    bind_rows(.id = "Distance") %>%
    mutate_if(is.numeric, round, 3) %>%
    select(-pvalBon) -> lis[[grp]]
}
```

```
## [1] "TR2_CTX20"
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
## [1] "TR4_CTX200"
```

```
## 'nperm' >= set of all permutations: complete enumeration.
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
## [1] "TR5_HV292.1"
```

```
## 'nperm' >= set of all permutations: complete enumeration.
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
```

```
## 'nperm' >= set of all permutations: complete enumeration.
```

```
## Set of permutations < 'minperm'. Generating entire set.
```

```
## [1] "TR1_CTX+HV292.120"
```

```
## 'nperm' >= set of all permutations: complete enumeration.
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
```

```
## 'nperm' >= set of all permutations: complete enumeration.
```

```
## Set of permutations < 'minperm'. Generating entire set.
```

```
## [1] "TR3_CTX+HV292.1200"
```

```
## 'nperm' >= set of all permutations: complete enumeration.
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
## [1] "CR_UNTREATED"
```

```
## 'nperm' >= set of all permutations: complete enumeration.
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

```r
lis %>% 
  bind_rows(.id = "Reactor_Treatment_Dose") %>% 
  filter(X1 == "pret" | "X2" == "pret") %>% 
  group_by(Reactor_Treatment_Dose, Distance) %>% 
  rstatix::adjust_pvalue(p.col = "pval", output.col = "pvalFDR") %>% 
  rstatix::add_significance(p.col = "pvalFDR")-> tmp_df

tmp_df %>% 
  DT::datatable()
```

```{=html}
<div class="datatables html-widget html-fill-item-overflow-hidden html-fill-item" id="htmlwidget-23dada22fc5642de2fc1" style="width:100%;height:auto;"></div>
<script type="application/json" data-for="htmlwidget-23dada22fc5642de2fc1">{"x":{"filter":"none","vertical":false,"data":[["1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","23","24","25","26","27","28","29","30","31","32"],["CR_UNTREATED","CR_UNTREATED","CR_UNTREATED","CR_UNTREATED","CR_UNTREATED","CR_UNTREATED","TR1_CTX+HV292.120","TR1_CTX+HV292.120","TR1_CTX+HV292.120","TR1_CTX+HV292.120","TR2_CTX20","TR2_CTX20","TR2_CTX20","TR2_CTX20","TR2_CTX20","TR2_CTX20","TR3_CTX+HV292.1200","TR3_CTX+HV292.1200","TR3_CTX+HV292.1200","TR3_CTX+HV292.1200","TR3_CTX+HV292.1200","TR3_CTX+HV292.1200","TR4_CTX200","TR4_CTX200","TR4_CTX200","TR4_CTX200","TR4_CTX200","TR4_CTX200","TR5_HV292.1","TR5_HV292.1","TR5_HV292.1","TR5_HV292.1"],["aichinson","aichinson","aichinson","bjaccard","bjaccard","bjaccard","aichinson","aichinson","bjaccard","bjaccard","aichinson","aichinson","aichinson","bjaccard","bjaccard","bjaccard","aichinson","aichinson","aichinson","bjaccard","bjaccard","bjaccard","aichinson","aichinson","aichinson","bjaccard","bjaccard","bjaccard","aichinson","aichinson","bjaccard","bjaccard"],["pret","pret","pret","pret","pret","pret","pret","pret","pret","pret","pret","pret","pret","pret","pret","pret","pret","pret","pret","pret","pret","pret","pret","pret","pret","pret","pret","pret","pret","pret","pret","pret"],["t1","t2","t3","t1","t2","t3","t2","t3","t2","t3","t1","t2","t3","t1","t2","t3","t1","t2","t3","t1","t2","t3","t1","t2","t3","t1","t2","t3","t1","t2","t1","t2"],[0.327,0.341,0.47,0.248,0.328,0.418,0.551,0.546,0.559,0.604,0.439,0.612,0.802,0.45,0.691,0.921,0.626,0.668,0.793,0.641,0.696,0.887,0.416,0.455,0.726,0.419,0.468,0.801,0.367,0.716,0.464,0.786],[0.75,0.1,0.25,0.75,0.1,0.25,0.1,0.25,0.1,0.25,0.1,0.038,0.25,0.1,0.032,0.25,0.1,0.1,0.25,0.1,0.1,0.25,0.333,0.333,0.333,0.333,0.333,0.333,0.25,0.25,0.25,0.25],[0.75,0.3,0.5,0.75,0.3,0.5,0.2,0.25,0.2,0.25,0.2,0.114,0.25,0.2,0.096,0.25,0.3,0.3,0.3,0.3,0.3,0.3,0.999,0.999,0.999,0.999,0.999,0.999,0.5,0.5,0.5,0.5],["ns","ns","ns","ns","ns","ns","ns","ns","ns","ns","ns","ns","ns","ns","ns","ns","ns","ns","ns","ns","ns","ns","ns","ns","ns","ns","ns","ns","ns","ns","ns","ns"]],"container":"<table class=\"display\">\n  <thead>\n    <tr>\n      <th> <\/th>\n      <th>Reactor_Treatment_Dose<\/th>\n      <th>Distance<\/th>\n      <th>X1<\/th>\n      <th>X2<\/th>\n      <th>R2<\/th>\n      <th>pval<\/th>\n      <th>pvalFDR<\/th>\n      <th>pvalFDR.signif<\/th>\n    <\/tr>\n  <\/thead>\n<\/table>","options":{"columnDefs":[{"className":"dt-right","targets":[5,6,7]},{"orderable":false,"targets":0}],"order":[],"autoWidth":false,"orderClasses":false}},"evals":[],"jsHooks":[]}</script>
```

```r
tmp_df %>% 
  xlsx::write.xlsx(file = out_xlsx,
                   append = TRUE,
                   sheetName = "human1_F2-beta-vs_pret",
                   showNA = TRUE)
```



```r
lis %>% 
  bind_rows(.id = "Reactor_Treatment_Dose") %>% 
  group_by(Reactor_Treatment_Dose, Distance) %>% 
  rstatix::adjust_pvalue(p.col = "pval", output.col = "pvalFDR") %>% 
  rstatix::add_significance(p.col = "pvalFDR") -> tmp_df

tmp_df %>% 
  DT::datatable()
```

```{=html}
<div class="datatables html-widget html-fill-item-overflow-hidden html-fill-item" id="htmlwidget-6e16c2883c26dd1edd48" style="width:100%;height:auto;"></div>
<script type="application/json" data-for="htmlwidget-6e16c2883c26dd1edd48">{"x":{"filter":"none","vertical":false,"data":[["1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","23","24","25","26","27","28","29","30","31","32","33","34","35","36","37","38","39","40","41","42","43","44","45","46","47","48","49","50","51","52","53","54","55","56","57","58","59","60"],["CR_UNTREATED","CR_UNTREATED","CR_UNTREATED","CR_UNTREATED","CR_UNTREATED","CR_UNTREATED","CR_UNTREATED","CR_UNTREATED","CR_UNTREATED","CR_UNTREATED","CR_UNTREATED","CR_UNTREATED","TR1_CTX+HV292.120","TR1_CTX+HV292.120","TR1_CTX+HV292.120","TR1_CTX+HV292.120","TR1_CTX+HV292.120","TR1_CTX+HV292.120","TR2_CTX20","TR2_CTX20","TR2_CTX20","TR2_CTX20","TR2_CTX20","TR2_CTX20","TR2_CTX20","TR2_CTX20","TR2_CTX20","TR2_CTX20","TR2_CTX20","TR2_CTX20","TR3_CTX+HV292.1200","TR3_CTX+HV292.1200","TR3_CTX+HV292.1200","TR3_CTX+HV292.1200","TR3_CTX+HV292.1200","TR3_CTX+HV292.1200","TR3_CTX+HV292.1200","TR3_CTX+HV292.1200","TR3_CTX+HV292.1200","TR3_CTX+HV292.1200","TR3_CTX+HV292.1200","TR3_CTX+HV292.1200","TR4_CTX200","TR4_CTX200","TR4_CTX200","TR4_CTX200","TR4_CTX200","TR4_CTX200","TR4_CTX200","TR4_CTX200","TR4_CTX200","TR4_CTX200","TR4_CTX200","TR4_CTX200","TR5_HV292.1","TR5_HV292.1","TR5_HV292.1","TR5_HV292.1","TR5_HV292.1","TR5_HV292.1"],["aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","aichinson","aichinson","aichinson","bjaccard","bjaccard","bjaccard","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","aichinson","aichinson","aichinson","bjaccard","bjaccard","bjaccard"],["pret","pret","pret","t1","t1","t2","pret","pret","pret","t1","t1","t2","pret","pret","t2","pret","pret","t2","pret","pret","pret","t1","t1","t2","pret","pret","pret","t1","t1","t2","pret","pret","pret","t1","t1","t2","pret","pret","pret","t1","t1","t2","pret","pret","pret","t1","t1","t2","pret","pret","pret","t1","t1","t2","pret","pret","t1","pret","pret","t1"],["t1","t2","t3","t2","t3","t3","t1","t2","t3","t2","t3","t3","t2","t3","t3","t2","t3","t3","t1","t2","t3","t2","t3","t3","t1","t2","t3","t2","t3","t3","t1","t2","t3","t2","t3","t3","t1","t2","t3","t2","t3","t3","t1","t2","t3","t2","t3","t3","t1","t2","t3","t2","t3","t3","t1","t2","t2","t1","t2","t2"],[0.327,0.341,0.47,0.285,1,0.215,0.248,0.328,0.418,0.25,1,0.232,0.551,0.546,0.509,0.559,0.604,0.47,0.439,0.612,0.802,0.352,0.408,0.631,0.45,0.691,0.921,0.399,0.381,0.719,0.626,0.668,0.793,0.386,0.638,0.344,0.641,0.696,0.887,0.392,0.663,0.378,0.416,0.455,0.726,0.211,0.547,0.388,0.419,0.468,0.801,0.136,0.621,0.426,0.367,0.716,1,0.464,0.786,1],[0.75,0.1,0.25,0.5,null,1,0.75,0.1,0.25,1,null,1,0.1,0.25,0.667,0.1,0.25,0.667,0.1,0.038,0.25,0.067,1,0.2,0.1,0.032,0.25,0.067,1,0.2,0.1,0.1,0.25,0.1,0.333,0.5,0.1,0.1,0.25,0.1,0.333,0.5,0.333,0.333,0.333,1,0.333,1,0.333,0.333,0.333,1,0.333,1,0.25,0.25,null,0.25,0.25,null],[1,0.5,1,1,null,1,1,0.5,1,1,null,1,0.3,0.5,0.667,0.3,0.5,0.667,0.4,0.228,0.6,0.335,1,0.6,0.4,0.192,0.6,0.335,1,0.6,0.6,0.6,0.75,0.6,0.75,0.75,0.6,0.6,0.75,0.6,0.75,0.75,1,1,1,1,1,1,1,1,1,1,1,1,0.5,0.5,null,0.5,0.5,null],["ns","ns","ns","ns","","ns","ns","ns","ns","ns","","ns","ns","ns","ns","ns","ns","ns","ns","ns","ns","ns","ns","ns","ns","ns","ns","ns","ns","ns","ns","ns","ns","ns","ns","ns","ns","ns","ns","ns","ns","ns","ns","ns","ns","ns","ns","ns","ns","ns","ns","ns","ns","ns","ns","ns","","ns","ns",""]],"container":"<table class=\"display\">\n  <thead>\n    <tr>\n      <th> <\/th>\n      <th>Reactor_Treatment_Dose<\/th>\n      <th>Distance<\/th>\n      <th>X1<\/th>\n      <th>X2<\/th>\n      <th>R2<\/th>\n      <th>pval<\/th>\n      <th>pvalFDR<\/th>\n      <th>pvalFDR.signif<\/th>\n    <\/tr>\n  <\/thead>\n<\/table>","options":{"columnDefs":[{"className":"dt-right","targets":[5,6,7]},{"orderable":false,"targets":0}],"order":[],"autoWidth":false,"orderClasses":false}},"evals":[],"jsHooks":[]}</script>
```

```r
tmp_df %>% 
  xlsx::write.xlsx(file = out_xlsx,
                   append = TRUE,
                   sheetName = "human1_F2-beta-between-periods",
                   showNA = TRUE)
```



```r
# for (grp in names(lis)){
#   print(grp)
#   
# lapply(
#   dist,
#   FUN = phyloseq_TW,
#   physeq =  ps_tmp_forbdiv %>%
#       # subset_samples(.,  eval(as.name(ranks)) %in% tax) %>% 
#       subset_samples(., Reactor_Treatment_Dose == grp),
#     variable = "period",
#   nrep = 999,
#   strata = NULL
# ) %>%
#   bind_rows(.id = "Distance") %>%
#   mutate_if(is.numeric, round, 5) -> lis[[grp]]
# }
# 
# 
# lis %>% 
#   bind_rows(.id = "Reactor_Treatment_Dose") %>% 
#   xlsx::write.xlsx(file = out_xlsx,
#                    append = TRUE,
#                    sheetName = "human1_F2-beta-between-periodsTW",
#                    showNA = TRUE)
```


### Human2:


```r
ps_rare %>%
  subset_samples(Model2 == "Human2"  & Treatment != "DONOR" ) -> ps_tmp_forbdiv
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

![](NRP72_draft_16S_Figs_beta_PCoA_files/figure-html/unnamed-chunk-89-1.png)<!-- -->


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
<div class="datatables html-widget html-fill-item-overflow-hidden html-fill-item" id="htmlwidget-23e11579ab285b5c7c2a" style="width:100%;height:auto;"></div>
<script type="application/json" data-for="htmlwidget-23e11579ab285b5c7c2a">{"x":{"filter":"none","vertical":false,"data":[["ASV0002","ASV0003","ASV0020","ASV0040","ASV0032","ASV0085","ASV0124","ASV0129"],[0.525535181023166,0.920169096922261,-0.358554123676786,-0.116195913277726,-0.385498361508464,0.61223295617904,0.738008656230995,0.681278011291429],[-0.618870077158341,-0.128496671144786,0.82065886361688,0.684754637802372,0.610396402616658,-0.495623606411118,0.253095005062396,0.318522455886272],[0.001,0.001,0.001,0.001,0.001,0.001,0.001,0.001],[0.659187398895023,0.863222561426021,0.802042030038578,0.482390404254302,0.521192755053068,0.620471951863689,0.608713858259413,0.565596283573027],["Bacteria","Bacteria","Bacteria","Bacteria","Bacteria","Bacteria","Bacteria","Bacteria"],["Bacteroidota","Firmicutes","Firmicutes","Firmicutes","Proteobacteria","Firmicutes","Firmicutes","Actinobacteriota"],["Bacteroidia","Clostridia","Bacilli","Negativicutes","Gammaproteobacteria","Clostridia","Clostridia","Coriobacteriia"],["Bacteroidales","Oscillospirales","Lactobacillales","Veillonellales-Selenomonadales","Burkholderiales","Oscillospirales","Monoglobales","Coriobacteriales"],["Bacteroidaceae","Ruminococcaceae","Lactobacillaceae","Veillonellaceae","Sutterellaceae","Oscillospiraceae","Monoglobaceae","Coriobacteriaceae"],[null,null,null,null,null,null,null,null],[null,null,null,null,null,null,null,null],[null,null,null,null,null,null,null,null],[0.0026,0.0026,0.0026,0.0026,0.0026,0.0026,0.0026,0.0026]],"container":"<table class=\"display\">\n  <thead>\n    <tr>\n      <th> <\/th>\n      <th>Axis.1<\/th>\n      <th>Axis.2<\/th>\n      <th>pval<\/th>\n      <th>r<\/th>\n      <th>Kingdom<\/th>\n      <th>Phylum<\/th>\n      <th>Class<\/th>\n      <th>Order<\/th>\n      <th>tax_rank_plot<\/th>\n      <th>Genus<\/th>\n      <th>Species<\/th>\n      <th>Strain<\/th>\n      <th>pval.adj<\/th>\n    <\/tr>\n  <\/thead>\n<\/table>","options":{"columnDefs":[{"className":"dt-right","targets":[1,2,3,4,13]},{"orderable":false,"targets":0}],"order":[],"autoWidth":false,"orderClasses":false}},"evals":[],"jsHooks":[]}</script>
```


```r
out$plot -> PCoA_tmp2_fam

PCoA_tmp2 + 
  theme(legend.position = "none") -> PCoA_tmp2_noleg

PCoA_tmp2_noleg
```

![](NRP72_draft_16S_Figs_beta_PCoA_files/figure-html/unnamed-chunk-91-1.png)<!-- -->

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
<div class="datatables html-widget html-fill-item-overflow-hidden html-fill-item" id="htmlwidget-26eca46a20528b6114d0" style="width:100%;height:auto;"></div>
<script type="application/json" data-for="htmlwidget-26eca46a20528b6114d0">{"x":{"filter":"none","vertical":false,"data":[["1","2","3","4","5","6"],["Formiat_mM","Acetat_mM","Propionat_mM","Butyrat_mM","Valerat_mM","Total_SCFA_mM"],[-0.549090636476436,0.916195874386532,0.884304826519494,0.950599581786949,0.883991857526509,0.915344362733472],[0.447248992918467,-0.211750385880216,-0.0384601364958492,-0.168139704173301,0.291267990063294,-0.302627172723072],["Formiat_mM","Acetat_mM","Propionat_mM","Butyrat_mM","Valerat_mM","Total_SCFA_mM"],[0.001,0.001,0.001,0.001,0.001,0.001],[0.501532188732681,0.884253106163322,0.783474208304952,0.931910525013008,0.866278646208678,0.929438508058306]],"container":"<table class=\"display\">\n  <thead>\n    <tr>\n      <th> <\/th>\n      <th>label<\/th>\n      <th>Axis.1<\/th>\n      <th>Axis.2<\/th>\n      <th>id<\/th>\n      <th>pval<\/th>\n      <th>r<\/th>\n    <\/tr>\n  <\/thead>\n<\/table>","options":{"columnDefs":[{"className":"dt-right","targets":[2,3,5,6]},{"orderable":false,"targets":0}],"order":[],"autoWidth":false,"orderClasses":false}},"evals":[],"jsHooks":[]}</script>
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

![](NRP72_draft_16S_Figs_beta_PCoA_files/figure-html/unnamed-chunk-94-1.png)<!-- -->


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
<div class="datatables html-widget html-fill-item-overflow-hidden html-fill-item" id="htmlwidget-fe4b5d221dd8b1ae7542" style="width:100%;height:auto;"></div>
<script type="application/json" data-for="htmlwidget-fe4b5d221dd8b1ae7542">{"x":{"filter":"none","vertical":false,"data":[["ASV0002","ASV0003","ASV0020","ASV0040","ASV0032","ASV0085","ASV0124","ASV0129"],[0.656938471951339,0.920440483900327,-0.410941714999905,-0.14447751323852,-0.439146111340108,0.741210046290733,0.722410876054175,0.676490209988427],[-0.409232258728033,0.0417566788403603,0.798012284558517,0.770944906213352,0.605987654237452,-0.207815644309558,0.208667051824549,0.259292335598343],[0.001,0.001,0.001,0.001,0.001,0.001,0.001,0.001],[0.599039197513408,0.848954304630445,0.805696699433367,0.6152298002479,0.560070344193349,0.592579674742107,0.56541941235851,0.524871519510229],["Bacteria","Bacteria","Bacteria","Bacteria","Bacteria","Bacteria","Bacteria","Bacteria"],["Bacteroidota","Firmicutes","Firmicutes","Firmicutes","Proteobacteria","Firmicutes","Firmicutes","Actinobacteriota"],["Bacteroidia","Clostridia","Bacilli","Negativicutes","Gammaproteobacteria","Clostridia","Clostridia","Coriobacteriia"],["Bacteroidales","Oscillospirales","Lactobacillales","Veillonellales-Selenomonadales","Burkholderiales","Oscillospirales","Monoglobales","Coriobacteriales"],["Bacteroidaceae","Ruminococcaceae","Lactobacillaceae","Veillonellaceae","Sutterellaceae","Oscillospiraceae","Monoglobaceae","Coriobacteriaceae"],[null,null,null,null,null,null,null,null],[null,null,null,null,null,null,null,null],[null,null,null,null,null,null,null,null],[0.0026,0.0026,0.0026,0.0026,0.0026,0.0026,0.0026,0.0026]],"container":"<table class=\"display\">\n  <thead>\n    <tr>\n      <th> <\/th>\n      <th>Axis.1<\/th>\n      <th>Axis.2<\/th>\n      <th>pval<\/th>\n      <th>r<\/th>\n      <th>Kingdom<\/th>\n      <th>Phylum<\/th>\n      <th>Class<\/th>\n      <th>Order<\/th>\n      <th>tax_rank_plot<\/th>\n      <th>Genus<\/th>\n      <th>Species<\/th>\n      <th>Strain<\/th>\n      <th>pval.adj<\/th>\n    <\/tr>\n  <\/thead>\n<\/table>","options":{"columnDefs":[{"className":"dt-right","targets":[1,2,3,4,13]},{"orderable":false,"targets":0}],"order":[],"autoWidth":false,"orderClasses":false}},"evals":[],"jsHooks":[]}</script>
```


```r
out$plot -> PCoA_tmp2_fam

PCoA_tmp2 + 
  theme(legend.position = "none") -> PCoA_tmp2_noleg

PCoA_tmp2_noleg
```

![](NRP72_draft_16S_Figs_beta_PCoA_files/figure-html/unnamed-chunk-96-1.png)<!-- -->

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
<div class="datatables html-widget html-fill-item-overflow-hidden html-fill-item" id="htmlwidget-7278515fd065a9a7a95f" style="width:100%;height:auto;"></div>
<script type="application/json" data-for="htmlwidget-7278515fd065a9a7a95f">{"x":{"filter":"none","vertical":false,"data":[["1","2","3","4","5","6"],["Formiat_mM","Acetat_mM","Propionat_mM","Butyrat_mM","Valerat_mM","Total_SCFA_mM"],[-0.605272295297702,0.932524691058002,0.861081805956696,0.95281725848877,0.86910860346474,0.954567915935512],[0.31413075304474,-0.0312399320413977,0.209496296459151,-8.57176905462534e-06,0.252566031283776,-0.0932479278949179],["Formiat_mM","Acetat_mM","Propionat_mM","Butyrat_mM","Valerat_mM","Total_SCFA_mM"],[0.001,0.001,0.001,0.001,0.001,0.001],[0.465032681463404,0.870578232786773,0.785350574779747,0.90786072814753,0.819139364774867,0.919895082190162]],"container":"<table class=\"display\">\n  <thead>\n    <tr>\n      <th> <\/th>\n      <th>label<\/th>\n      <th>Axis.1<\/th>\n      <th>Axis.2<\/th>\n      <th>id<\/th>\n      <th>pval<\/th>\n      <th>r<\/th>\n    <\/tr>\n  <\/thead>\n<\/table>","options":{"columnDefs":[{"className":"dt-right","targets":[2,3,5,6]},{"orderable":false,"targets":0}],"order":[],"autoWidth":false,"orderClasses":false}},"evals":[],"jsHooks":[]}</script>
```


```r
env_out$plot %>% 
  export::graph2ppt(append = TRUE,
                    width = 317.48031496 * 1,
                    height = 0.618 * 317.48031496 * 1 , paper = "A4",  scaling = 2,
                    file = out_pptx)
```


####  Statistical evaluation:

##### pret - Day_of_Treatment <= 0:


```r
lapply(
  dist,
  FUN = phyloseq_adonis2,
  physeq = ps_tmp_forbdiv %>%
    subset_samples(Day_of_Treatment <= 0),
  formula = paste0(c("Reactor_Treatment_Dose"), collapse=" + "),
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
<div class="datatables html-widget html-fill-item-overflow-hidden html-fill-item" id="htmlwidget-198c93e593dcd35da8ef" style="width:100%;height:auto;"></div>
<script type="application/json" data-for="htmlwidget-198c93e593dcd35da8ef">{"x":{"filter":"none","vertical":false,"data":[["1","2","3","4","5","6"],["bjaccard","bjaccard","bjaccard","aichinson","aichinson","aichinson"],["Reactor_Treatment_Dose","Residual","Total","Reactor_Treatment_Dose","Residual","Total"],[6,28,34,6,28,34],[0.498,0.83,1.328,0.523,0.923,1.446],[0.375,0.625,1,0.362,0.638,1],[2.803,null,null,2.643,null,null],[0.001,null,null,0.001,null,null]],"container":"<table class=\"display\">\n  <thead>\n    <tr>\n      <th> <\/th>\n      <th>Distance<\/th>\n      <th>terms<\/th>\n      <th>Df<\/th>\n      <th>SumOfSqs<\/th>\n      <th>R2<\/th>\n      <th>F<\/th>\n      <th>Pr(&gt;F)<\/th>\n    <\/tr>\n  <\/thead>\n<\/table>","options":{"columnDefs":[{"className":"dt-right","targets":[3,4,5,6,7]},{"orderable":false,"targets":0}],"order":[],"autoWidth":false,"orderClasses":false}},"evals":[],"jsHooks":[]}</script>
```

```r
tmp_df %>% 
  xlsx::write.xlsx(file = out_xlsx,
                   append = TRUE,
                   sheetName = "human2-beta-pretreat",
                   showNA = TRUE)
```



```r
lapply(
  dist,
  FUN = phyloseq_adonis2,
  physeq = ps_tmp_forbdiv %>%
    subset_samples(Day_of_Treatment <= 0),
  formula = paste0(c("Reactor_Treatment_Dose", "Day_of_Treatment"), collapse=" + "),
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
<div class="datatables html-widget html-fill-item-overflow-hidden html-fill-item" id="htmlwidget-87a6c950be51e7f26976" style="width:100%;height:auto;"></div>
<script type="application/json" data-for="htmlwidget-87a6c950be51e7f26976">{"x":{"filter":"none","vertical":false,"data":[["1","2","3","4","5","6","7","8"],["bjaccard","bjaccard","bjaccard","bjaccard","aichinson","aichinson","aichinson","aichinson"],["Reactor_Treatment_Dose","Day_of_Treatment","Residual","Total","Reactor_Treatment_Dose","Day_of_Treatment","Residual","Total"],[6,1,27,34,6,1,27,34],[0.498,0.063,0.767,1.328,0.523,0.06,0.863,1.446],[0.375,0.047,0.577,1,0.362,0.042,0.597,1],[2.924,2.216,null,null,2.727,1.89,null,null],[0.001,0.003,null,null,0.001,0.005,null,null]],"container":"<table class=\"display\">\n  <thead>\n    <tr>\n      <th> <\/th>\n      <th>Distance<\/th>\n      <th>terms<\/th>\n      <th>Df<\/th>\n      <th>SumOfSqs<\/th>\n      <th>R2<\/th>\n      <th>F<\/th>\n      <th>Pr(&gt;F)<\/th>\n    <\/tr>\n  <\/thead>\n<\/table>","options":{"columnDefs":[{"className":"dt-right","targets":[3,4,5,6,7]},{"orderable":false,"targets":0}],"order":[],"autoWidth":false,"orderClasses":false}},"evals":[],"jsHooks":[]}</script>
```

```r
tmp_df %>% 
  xlsx::write.xlsx(file = out_xlsx,
                   append = TRUE,
                   sheetName = "human2-beta-pretreat-time-day",
                   showNA = TRUE)
```


```r
lapply(
  dist,
  FUN = physeq_pairwise_permanovas_adonis2,
  physeq = ps_tmp_forbdiv %>%
    subset_samples(Day_of_Treatment <= 0),
  compare_header = "Reactor_Treatment_Dose",
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
<div class="datatables html-widget html-fill-item-overflow-hidden html-fill-item" id="htmlwidget-0a706b8192cc8c1e9d72" style="width:100%;height:auto;"></div>
<script type="application/json" data-for="htmlwidget-0a706b8192cc8c1e9d72">{"x":{"filter":"none","vertical":false,"data":[["1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","23","24","25","26","27","28","29","30","31","32","33","34","35","36","37","38","39","40","41","42"],["bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson"],["CR_UNTREATED","CR_UNTREATED","CR_UNTREATED","CR_UNTREATED","CR_UNTREATED","CR_UNTREATED","TR3_CTX200","TR3_CTX200","TR3_CTX200","TR3_CTX200","TR3_CTX200","TR5_CTX+HV292.1200","TR5_CTX+HV292.1200","TR5_CTX+HV292.1200","TR5_CTX+HV292.1200","TR7_HV292.1","TR7_HV292.1","TR7_HV292.1","TR2_VAN600","TR2_VAN600","TR4_VAN+CCUG59168600","CR_UNTREATED","CR_UNTREATED","CR_UNTREATED","CR_UNTREATED","CR_UNTREATED","CR_UNTREATED","TR3_CTX200","TR3_CTX200","TR3_CTX200","TR3_CTX200","TR3_CTX200","TR5_CTX+HV292.1200","TR5_CTX+HV292.1200","TR5_CTX+HV292.1200","TR5_CTX+HV292.1200","TR7_HV292.1","TR7_HV292.1","TR7_HV292.1","TR2_VAN600","TR2_VAN600","TR4_VAN+CCUG59168600"],["TR3_CTX200","TR5_CTX+HV292.1200","TR7_HV292.1","TR2_VAN600","TR4_VAN+CCUG59168600","TR6_CCUG59168","TR5_CTX+HV292.1200","TR7_HV292.1","TR2_VAN600","TR4_VAN+CCUG59168600","TR6_CCUG59168","TR7_HV292.1","TR2_VAN600","TR4_VAN+CCUG59168600","TR6_CCUG59168","TR2_VAN600","TR4_VAN+CCUG59168600","TR6_CCUG59168","TR4_VAN+CCUG59168600","TR6_CCUG59168","TR6_CCUG59168","TR3_CTX200","TR5_CTX+HV292.1200","TR7_HV292.1","TR2_VAN600","TR4_VAN+CCUG59168600","TR6_CCUG59168","TR5_CTX+HV292.1200","TR7_HV292.1","TR2_VAN600","TR4_VAN+CCUG59168600","TR6_CCUG59168","TR7_HV292.1","TR2_VAN600","TR4_VAN+CCUG59168600","TR6_CCUG59168","TR2_VAN600","TR4_VAN+CCUG59168600","TR6_CCUG59168","TR4_VAN+CCUG59168600","TR6_CCUG59168","TR6_CCUG59168"],[0.311,0.221,0.204,0.182,0.198,0.241,0.364,0.373,0.206,0.2,0.387,0.278,0.219,0.201,0.263,0.25,0.254,0.33,0.161,0.299,0.28,0.302,0.221,0.222,0.207,0.212,0.228,0.305,0.338,0.208,0.212,0.312,0.262,0.221,0.208,0.228,0.266,0.258,0.301,0.171,0.262,0.262],[0.011,0.007,0.012,0.023,0.012,0.016,0.011,0.006,0.037,0.008,0.009,0.005,0.014,0.023,0.008,0.004,0.009,0.006,0.125,0.007,0.006,0.011,0.014,0.006,0.009,0.01,0.008,0.01,0.01,0.007,0.009,0.007,0.008,0.008,0.008,0.012,0.01,0.005,0.01,0.046,0.012,0.009],[0.018,0.023,0.017,0.026,0.017,0.02,0.018,0.032,0.039,0.02,0.018,0.052,0.018,0.026,0.02,0.084,0.018,0.032,0.125,0.023,0.032,0.014,0.015,0.063,0.019,0.015,0.026,0.015,0.015,0.042,0.019,0.042,0.026,0.026,0.026,0.014,0.015,0.105,0.015,0.046,0.014,0.019]],"container":"<table class=\"display\">\n  <thead>\n    <tr>\n      <th> <\/th>\n      <th>Distance<\/th>\n      <th>X1<\/th>\n      <th>X2<\/th>\n      <th>R2<\/th>\n      <th>pval<\/th>\n      <th>pvalFDR<\/th>\n    <\/tr>\n  <\/thead>\n<\/table>","options":{"columnDefs":[{"className":"dt-right","targets":[4,5,6]},{"orderable":false,"targets":0}],"order":[],"autoWidth":false,"orderClasses":false}},"evals":[],"jsHooks":[]}</script>
```

```r
tmp_df %>% 
  xlsx::write.xlsx(file = out_xlsx,
                   append = TRUE,
                   sheetName = "human2-beta-pretreat-pw",
                   showNA = TRUE)
```


```r
lapply(
  dist,
  FUN = phyloseq_TW,
  physeq =  ps_tmp_forbdiv %>%
    subset_samples(Day_of_Treatment <= 0),
  variable = "Reactor_Treatment_Dose",
  nrep = 999,
  strata = NULL
) %>%
  bind_rows(.id = "Distance") %>%
  mutate_if(is.numeric, round, 3)   -> tmp_df

tmp_df %>% 
  DT::datatable()
```

```{=html}
<div class="datatables html-widget html-fill-item-overflow-hidden html-fill-item" id="htmlwidget-155dc90fc4a3b1a0cf9b" style="width:100%;height:auto;"></div>
<script type="application/json" data-for="htmlwidget-155dc90fc4a3b1a0cf9b">{"x":{"filter":"none","vertical":false,"data":[["1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","23","24","25","26","27","28","29","30","31","32","33","34","35","36","37","38","39","40","41","42"],["bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson"],[["CR_UNTREATED"],["CR_UNTREATED"],["CR_UNTREATED"],["CR_UNTREATED"],["CR_UNTREATED"],["CR_UNTREATED"],["TR3_CTX200"],["TR3_CTX200"],["TR3_CTX200"],["TR3_CTX200"],["TR3_CTX200"],["TR5_CTX+HV292.1200"],["TR5_CTX+HV292.1200"],["TR5_CTX+HV292.1200"],["TR5_CTX+HV292.1200"],["TR7_HV292.1"],["TR7_HV292.1"],["TR7_HV292.1"],["TR2_VAN600"],["TR2_VAN600"],["TR4_VAN+CCUG59168600"],["CR_UNTREATED"],["CR_UNTREATED"],["CR_UNTREATED"],["CR_UNTREATED"],["CR_UNTREATED"],["CR_UNTREATED"],["TR3_CTX200"],["TR3_CTX200"],["TR3_CTX200"],["TR3_CTX200"],["TR3_CTX200"],["TR5_CTX+HV292.1200"],["TR5_CTX+HV292.1200"],["TR5_CTX+HV292.1200"],["TR5_CTX+HV292.1200"],["TR7_HV292.1"],["TR7_HV292.1"],["TR7_HV292.1"],["TR2_VAN600"],["TR2_VAN600"],["TR4_VAN+CCUG59168600"]],[["TR3_CTX200"],["TR5_CTX+HV292.1200"],["TR7_HV292.1"],["TR2_VAN600"],["TR4_VAN+CCUG59168600"],["TR6_CCUG59168"],["TR5_CTX+HV292.1200"],["TR7_HV292.1"],["TR2_VAN600"],["TR4_VAN+CCUG59168600"],["TR6_CCUG59168"],["TR7_HV292.1"],["TR2_VAN600"],["TR4_VAN+CCUG59168600"],["TR6_CCUG59168"],["TR2_VAN600"],["TR4_VAN+CCUG59168600"],["TR6_CCUG59168"],["TR4_VAN+CCUG59168600"],["TR6_CCUG59168"],["TR6_CCUG59168"],["TR3_CTX200"],["TR5_CTX+HV292.1200"],["TR7_HV292.1"],["TR2_VAN600"],["TR4_VAN+CCUG59168600"],["TR6_CCUG59168"],["TR5_CTX+HV292.1200"],["TR7_HV292.1"],["TR2_VAN600"],["TR4_VAN+CCUG59168600"],["TR6_CCUG59168"],["TR7_HV292.1"],["TR2_VAN600"],["TR4_VAN+CCUG59168600"],["TR6_CCUG59168"],["TR2_VAN600"],["TR4_VAN+CCUG59168600"],["TR6_CCUG59168"],["TR4_VAN+CCUG59168600"],["TR6_CCUG59168"],["TR6_CCUG59168"]],[[5],[5],[5],[5],[5],[5],[5],[5],[5],[5],[5],[5],[5],[5],[5],[5],[5],[5],[5],[5],[5],[5],[5],[5],[5],[5],[5],[5],[5],[5],[5],[5],[5],[5],[5],[5],[5],[5],[5],[5],[5],[5]],[[5],[5],[5],[5],[5],[5],[5],[5],[5],[5],[5],[5],[5],[5],[5],[5],[5],[5],[5],[5],[5],[5],[5],[5],[5],[5],[5],[5],[5],[5],[5],[5],[5],[5],[5],[5],[5],[5],[5],[5],[5],[5]],[[0.01],[0.006],[0.021],[0.021],[0.011],[0.019],[0.01],[0.01],[0.032],[0.009],[0.012],[0.01],[0.013],[0.028],[0.006],[0.01],[0.01],[0.005],[0.147],[0.01],[0.005],[0.005],[0.008],[0.011],[0.008],[0.008],[0.011],[0.014],[0.011],[0.012],[0.006],[0.01],[0.011],[0.009],[0.007],[0.005],[0.007],[0.009],[0.013],[0.071],[0.008],[0.006]],[[3.61161142649598],[2.27493433021227],[2.05416844599888],[1.78469617172638],[1.97135473957909],[2.54335385472751],[4.57619285805083],[4.75686635157071],[2.07378478252489],[2.00117504297029],[5.05892625722876],[3.07725875660178],[2.2451448601377],[2.01216546977223],[2.85541657237427],[2.66816457358006],[2.7244236272209],[3.9346078338194],[1.53159832300531],[3.41724271596055],[3.11032306927843],[3.46448005798299],[2.27427737398797],[2.28524891228454],[2.08337403062385],[2.1538790758296],[2.36755458146576],[3.51845495766692],[4.07955806969518],[2.10148373749945],[2.14891083287175],[3.63532244289794],[2.83408967522825],[2.2689911483374],[2.09820739733884],[2.36049700500679],[2.89918097592929],[2.78553270047736],[3.44716042997851],[1.64716525867757],[2.84723613648874],[2.8379665076349]],[[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999]],[0.018,0.018,0.025,0.025,0.018,0.025,0.018,0.018,0.034,0.018,0.018,0.018,0.018,0.031,0.018,0.018,0.018,0.018,0.147,0.018,0.018,0.014,0.014,0.014,0.014,0.014,0.014,0.015,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.071,0.014,0.014]],"container":"<table class=\"display\">\n  <thead>\n    <tr>\n      <th> <\/th>\n      <th>Distance<\/th>\n      <th>Level1<\/th>\n      <th>Level2<\/th>\n      <th>N1<\/th>\n      <th>N2<\/th>\n      <th>p.value<\/th>\n      <th>tw2.stat<\/th>\n      <th>nrep<\/th>\n      <th>pvalFDR<\/th>\n    <\/tr>\n  <\/thead>\n<\/table>","options":{"columnDefs":[{"className":"dt-right","targets":9},{"orderable":false,"targets":0}],"order":[],"autoWidth":false,"orderClasses":false}},"evals":[],"jsHooks":[]}</script>
```

```r
tmp_df %>% 
  xlsx::write.xlsx(file = out_xlsx,
                   append = TRUE,
                   sheetName = "human2-beta-pretreat-pw2",
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

##### treat - Day_of_Treatment > 1:


```r
lapply(
  dist,
  FUN = phyloseq_adonis2,
  physeq = ps_tmp_forbdiv %>%
    subset_samples(Day_of_Treatment >1),
  formula = paste0(c("Reactor_Treatment_Dose"), collapse=" + "),
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
<div class="datatables html-widget html-fill-item-overflow-hidden html-fill-item" id="htmlwidget-b00bb49b4131cbaeb274" style="width:100%;height:auto;"></div>
<script type="application/json" data-for="htmlwidget-b00bb49b4131cbaeb274">{"x":{"filter":"none","vertical":false,"data":[["1","2","3","4","5","6"],["bjaccard","bjaccard","bjaccard","aichinson","aichinson","aichinson"],["Reactor_Treatment_Dose","Residual","Total","Reactor_Treatment_Dose","Residual","Total"],[6,55,61,6,55,61],[6.81,5.555,12.366,4.768,2.74,7.508],[0.551,0.449,1,0.635,0.365,1],[11.237,null,null,15.95,null,null],[0.001,null,null,0.001,null,null]],"container":"<table class=\"display\">\n  <thead>\n    <tr>\n      <th> <\/th>\n      <th>Distance<\/th>\n      <th>terms<\/th>\n      <th>Df<\/th>\n      <th>SumOfSqs<\/th>\n      <th>R2<\/th>\n      <th>F<\/th>\n      <th>Pr(&gt;F)<\/th>\n    <\/tr>\n  <\/thead>\n<\/table>","options":{"columnDefs":[{"className":"dt-right","targets":[3,4,5,6,7]},{"orderable":false,"targets":0}],"order":[],"autoWidth":false,"orderClasses":false}},"evals":[],"jsHooks":[]}</script>
```

```r
tmp_df %>% 
  xlsx::write.xlsx(file = out_xlsx,
                   append = TRUE,
                   sheetName = "human2-beta-treat",
                   showNA = TRUE)
```


```r
lapply(
  dist,
  FUN = phyloseq_adonis2,
  physeq = ps_tmp_forbdiv %>%
    subset_samples(Day_of_Treatment >1),
  formula = paste0(c("Reactor_Treatment_Dose", "Day_of_Treatment"), collapse=" + "),
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
<div class="datatables html-widget html-fill-item-overflow-hidden html-fill-item" id="htmlwidget-85cf7f798e923a02ed60" style="width:100%;height:auto;"></div>
<script type="application/json" data-for="htmlwidget-85cf7f798e923a02ed60">{"x":{"filter":"none","vertical":false,"data":[["1","2","3","4","5","6","7","8"],["bjaccard","bjaccard","bjaccard","bjaccard","aichinson","aichinson","aichinson","aichinson"],["Reactor_Treatment_Dose","Day_of_Treatment","Residual","Total","Reactor_Treatment_Dose","Day_of_Treatment","Residual","Total"],[6,1,54,61,6,1,54,61],[6.807,0.271,5.285,12.366,4.765,0.171,2.569,7.508],[0.551,0.022,0.427,1,0.635,0.023,0.342,1],[11.593,2.765,null,null,16.695,3.603,null,null],[0.001,0.026,null,null,0.001,0.012,null,null]],"container":"<table class=\"display\">\n  <thead>\n    <tr>\n      <th> <\/th>\n      <th>Distance<\/th>\n      <th>terms<\/th>\n      <th>Df<\/th>\n      <th>SumOfSqs<\/th>\n      <th>R2<\/th>\n      <th>F<\/th>\n      <th>Pr(&gt;F)<\/th>\n    <\/tr>\n  <\/thead>\n<\/table>","options":{"columnDefs":[{"className":"dt-right","targets":[3,4,5,6,7]},{"orderable":false,"targets":0}],"order":[],"autoWidth":false,"orderClasses":false}},"evals":[],"jsHooks":[]}</script>
```

```r
tmp_df %>% 
  xlsx::write.xlsx(file = out_xlsx,
                   append = TRUE,
                   sheetName = "human2-beta-treat-time-day",
                   showNA = TRUE)
```


```r
lapply(
  dist,
  FUN = physeq_pairwise_permanovas_adonis2,
  physeq = ps_tmp_forbdiv %>%
    subset_samples(Day_of_Treatment > 1),
  compare_header = "Reactor_Treatment_Dose",
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
<div class="datatables html-widget html-fill-item-overflow-hidden html-fill-item" id="htmlwidget-f7202f4a63f76b2c09a9" style="width:100%;height:auto;"></div>
<script type="application/json" data-for="htmlwidget-f7202f4a63f76b2c09a9">{"x":{"filter":"none","vertical":false,"data":[["1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","23","24","25","26","27","28","29","30","31","32","33","34","35","36","37","38","39","40","41","42"],["bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson"],["CR_UNTREATED","CR_UNTREATED","CR_UNTREATED","CR_UNTREATED","CR_UNTREATED","CR_UNTREATED","TR3_CTX200","TR3_CTX200","TR3_CTX200","TR3_CTX200","TR3_CTX200","TR5_CTX+HV292.1200","TR5_CTX+HV292.1200","TR5_CTX+HV292.1200","TR5_CTX+HV292.1200","TR7_HV292.1","TR7_HV292.1","TR7_HV292.1","TR2_VAN600","TR2_VAN600","TR4_VAN+CCUG59168600","CR_UNTREATED","CR_UNTREATED","CR_UNTREATED","CR_UNTREATED","CR_UNTREATED","CR_UNTREATED","TR3_CTX200","TR3_CTX200","TR3_CTX200","TR3_CTX200","TR3_CTX200","TR5_CTX+HV292.1200","TR5_CTX+HV292.1200","TR5_CTX+HV292.1200","TR5_CTX+HV292.1200","TR7_HV292.1","TR7_HV292.1","TR7_HV292.1","TR2_VAN600","TR2_VAN600","TR4_VAN+CCUG59168600"],["TR3_CTX200","TR5_CTX+HV292.1200","TR7_HV292.1","TR2_VAN600","TR4_VAN+CCUG59168600","TR6_CCUG59168","TR5_CTX+HV292.1200","TR7_HV292.1","TR2_VAN600","TR4_VAN+CCUG59168600","TR6_CCUG59168","TR7_HV292.1","TR2_VAN600","TR4_VAN+CCUG59168600","TR6_CCUG59168","TR2_VAN600","TR4_VAN+CCUG59168600","TR6_CCUG59168","TR4_VAN+CCUG59168600","TR6_CCUG59168","TR6_CCUG59168","TR3_CTX200","TR5_CTX+HV292.1200","TR7_HV292.1","TR2_VAN600","TR4_VAN+CCUG59168600","TR6_CCUG59168","TR5_CTX+HV292.1200","TR7_HV292.1","TR2_VAN600","TR4_VAN+CCUG59168600","TR6_CCUG59168","TR7_HV292.1","TR2_VAN600","TR4_VAN+CCUG59168600","TR6_CCUG59168","TR2_VAN600","TR4_VAN+CCUG59168600","TR6_CCUG59168","TR4_VAN+CCUG59168600","TR6_CCUG59168","TR6_CCUG59168"],[0.533,0.401,0.24,0.561,0.583,0.28,0.161,0.495,0.238,0.267,0.543,0.363,0.285,0.314,0.415,0.519,0.561,0.249,0.127,0.569,0.611,0.641,0.459,0.238,0.649,0.667,0.269,0.186,0.596,0.318,0.352,0.656,0.425,0.34,0.364,0.481,0.607,0.639,0.234,0.185,0.662,0.697],[0.001,0.001,0.001,0.001,0.001,0.001,0.005,0.001,0.001,0.001,0.001,0.001,0.001,0.001,0.001,0.001,0.001,0.001,0.007,0.001,0.001,0.001,0.001,0.001,0.001,0.001,0.001,0.001,0.002,0.001,0.001,0.001,0.001,0.001,0.001,0.001,0.001,0.001,0.001,0.001,0.001,0.001],[0.002,0.002,0.002,0.002,0.002,0.002,0.005,0.002,0.002,0.002,0.002,0.002,0.002,0.002,0.002,0.002,0.002,0.002,0.007,0.002,0.002,0.002,0.002,0.002,0.002,0.002,0.002,0.002,0.002,0.002,0.002,0.002,0.002,0.002,0.002,0.002,0.002,0.002,0.002,0.002,0.002,0.002]],"container":"<table class=\"display\">\n  <thead>\n    <tr>\n      <th> <\/th>\n      <th>Distance<\/th>\n      <th>X1<\/th>\n      <th>X2<\/th>\n      <th>R2<\/th>\n      <th>pval<\/th>\n      <th>pvalFDR<\/th>\n    <\/tr>\n  <\/thead>\n<\/table>","options":{"columnDefs":[{"className":"dt-right","targets":[4,5,6]},{"orderable":false,"targets":0}],"order":[],"autoWidth":false,"orderClasses":false}},"evals":[],"jsHooks":[]}</script>
```

```r
tmp_df %>% 
  xlsx::write.xlsx(file = out_xlsx,
                   append = TRUE,
                   sheetName = "human2-beta-treat-pw",
                   showNA = TRUE)
```


```r
lapply(
  dist,
  FUN = phyloseq_TW,
  physeq =  ps_tmp_forbdiv %>%
    subset_samples(Day_of_Treatment >1),
  variable = "Reactor_Treatment_Dose",
  nrep = 999,
  strata = NULL
) %>%
  bind_rows(.id = "Distance") %>%
  mutate_if(is.numeric, round, 3)   -> tmp_df

tmp_df %>% 
  DT::datatable()
```

```{=html}
<div class="datatables html-widget html-fill-item-overflow-hidden html-fill-item" id="htmlwidget-175be691d0d33a88348a" style="width:100%;height:auto;"></div>
<script type="application/json" data-for="htmlwidget-175be691d0d33a88348a">{"x":{"filter":"none","vertical":false,"data":[["1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","23","24","25","26","27","28","29","30","31","32","33","34","35","36","37","38","39","40","41","42"],["bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson"],[["CR_UNTREATED"],["CR_UNTREATED"],["CR_UNTREATED"],["CR_UNTREATED"],["CR_UNTREATED"],["CR_UNTREATED"],["TR3_CTX200"],["TR3_CTX200"],["TR3_CTX200"],["TR3_CTX200"],["TR3_CTX200"],["TR5_CTX+HV292.1200"],["TR5_CTX+HV292.1200"],["TR5_CTX+HV292.1200"],["TR5_CTX+HV292.1200"],["TR7_HV292.1"],["TR7_HV292.1"],["TR7_HV292.1"],["TR2_VAN600"],["TR2_VAN600"],["TR4_VAN+CCUG59168600"],["CR_UNTREATED"],["CR_UNTREATED"],["CR_UNTREATED"],["CR_UNTREATED"],["CR_UNTREATED"],["CR_UNTREATED"],["TR3_CTX200"],["TR3_CTX200"],["TR3_CTX200"],["TR3_CTX200"],["TR3_CTX200"],["TR5_CTX+HV292.1200"],["TR5_CTX+HV292.1200"],["TR5_CTX+HV292.1200"],["TR5_CTX+HV292.1200"],["TR7_HV292.1"],["TR7_HV292.1"],["TR7_HV292.1"],["TR2_VAN600"],["TR2_VAN600"],["TR4_VAN+CCUG59168600"]],[["TR3_CTX200"],["TR5_CTX+HV292.1200"],["TR7_HV292.1"],["TR2_VAN600"],["TR4_VAN+CCUG59168600"],["TR6_CCUG59168"],["TR5_CTX+HV292.1200"],["TR7_HV292.1"],["TR2_VAN600"],["TR4_VAN+CCUG59168600"],["TR6_CCUG59168"],["TR7_HV292.1"],["TR2_VAN600"],["TR4_VAN+CCUG59168600"],["TR6_CCUG59168"],["TR2_VAN600"],["TR4_VAN+CCUG59168600"],["TR6_CCUG59168"],["TR4_VAN+CCUG59168600"],["TR6_CCUG59168"],["TR6_CCUG59168"],["TR3_CTX200"],["TR5_CTX+HV292.1200"],["TR7_HV292.1"],["TR2_VAN600"],["TR4_VAN+CCUG59168600"],["TR6_CCUG59168"],["TR5_CTX+HV292.1200"],["TR7_HV292.1"],["TR2_VAN600"],["TR4_VAN+CCUG59168600"],["TR6_CCUG59168"],["TR7_HV292.1"],["TR2_VAN600"],["TR4_VAN+CCUG59168600"],["TR6_CCUG59168"],["TR2_VAN600"],["TR4_VAN+CCUG59168600"],["TR6_CCUG59168"],["TR4_VAN+CCUG59168600"],["TR6_CCUG59168"],["TR6_CCUG59168"]],[[9],[9],[9],[9],[9],[9],[9],[9],[9],[9],[9],[8],[8],[8],[8],[9],[9],[9],[9],[9],[9],[9],[9],[9],[9],[9],[9],[9],[9],[9],[9],[9],[8],[8],[8],[8],[9],[9],[9],[9],[9],[9]],[[9],[8],[9],[9],[9],[9],[8],[9],[9],[9],[9],[9],[9],[9],[9],[9],[9],[9],[9],[9],[9],[9],[8],[9],[9],[9],[9],[8],[9],[9],[9],[9],[9],[9],[9],[9],[9],[9],[9],[9],[9],[9]],[[0.001],[0.001],[0.001],[0.001],[0.001],[0.001],[0.003],[0.001],[0.001],[0.001],[0.001],[0.002],[0.001],[0.001],[0.001],[0.001],[0.001],[0.001],[0.007],[0.001],[0.001],[0.001],[0.001],[0.001],[0.001],[0.001],[0.001],[0.002],[0.001],[0.001],[0.001],[0.001],[0.001],[0.001],[0.001],[0.001],[0.001],[0.001],[0.001],[0.001],[0.001],[0.001]],[[18.2816588874981],[9.15338020543391],[5.06185266549415],[20.4103010341225],[22.4048296343546],[6.21677752428232],[2.83644065641568],[15.7111151619778],[4.99877598266036],[5.84221460105487],[18.9997630666877],[7.94357672569024],[5.88717517972496],[6.69673034336503],[9.63920250182291],[17.2757449393879],[20.4695843497377],[5.30779158939088],[2.32601358146117],[21.0987478982016],[25.1251183946221],[28.5778012204383],[12.0058538560295],[5.00787389919481],[29.6474174068086],[32.0437257789525],[5.89200113960892],[3.2895018858734],[23.5966772050598],[7.44789147346748],[8.68999194480714],[30.5096148098283],[10.5955889781548],[7.42800777989058],[8.13872158473742],[13.0596893305605],[24.6774842945271],[28.3333787305631],[4.89214808657102],[3.62119537870214],[31.2791596132917],[36.832090341497]],[[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999]],[0.001,0.001,0.001,0.001,0.001,0.001,0.003,0.001,0.001,0.001,0.001,0.002,0.001,0.001,0.001,0.001,0.001,0.001,0.007,0.001,0.001,0.001,0.001,0.001,0.001,0.001,0.001,0.002,0.001,0.001,0.001,0.001,0.001,0.001,0.001,0.001,0.001,0.001,0.001,0.001,0.001,0.001]],"container":"<table class=\"display\">\n  <thead>\n    <tr>\n      <th> <\/th>\n      <th>Distance<\/th>\n      <th>Level1<\/th>\n      <th>Level2<\/th>\n      <th>N1<\/th>\n      <th>N2<\/th>\n      <th>p.value<\/th>\n      <th>tw2.stat<\/th>\n      <th>nrep<\/th>\n      <th>pvalFDR<\/th>\n    <\/tr>\n  <\/thead>\n<\/table>","options":{"columnDefs":[{"className":"dt-right","targets":9},{"orderable":false,"targets":0}],"order":[],"autoWidth":false,"orderClasses":false}},"evals":[],"jsHooks":[]}</script>
```

```r
tmp_df %>% 
  xlsx::write.xlsx(file = out_xlsx,
                   append = TRUE,
                   sheetName = "human2-beta-treat-pw2",
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

##### Day_of_Treatment >= 3 & Day_of_Treatment <= 6:


```r
lapply(
  dist,
  FUN = phyloseq_adonis2,
  physeq = ps_tmp_forbdiv %>%
    subset_samples(Day_of_Treatment >= 3 & Day_of_Treatment <= 6),
  formula = paste0(c("Reactor_Treatment_Dose"), collapse=" + "),
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
<div class="datatables html-widget html-fill-item-overflow-hidden html-fill-item" id="htmlwidget-6b0771041c3145491948" style="width:100%;height:auto;"></div>
<script type="application/json" data-for="htmlwidget-6b0771041c3145491948">{"x":{"filter":"none","vertical":false,"data":[["1","2","3","4","5","6"],["bjaccard","bjaccard","bjaccard","aichinson","aichinson","aichinson"],["Reactor_Treatment_Dose","Residual","Total","Reactor_Treatment_Dose","Residual","Total"],[6,20,26,6,20,26],[3.92,1.816,5.736,2.42,0.746,3.166],[0.683,0.317,1,0.764,0.236,1],[7.194,null,null,10.81,null,null],[0.001,null,null,0.001,null,null]],"container":"<table class=\"display\">\n  <thead>\n    <tr>\n      <th> <\/th>\n      <th>Distance<\/th>\n      <th>terms<\/th>\n      <th>Df<\/th>\n      <th>SumOfSqs<\/th>\n      <th>R2<\/th>\n      <th>F<\/th>\n      <th>Pr(&gt;F)<\/th>\n    <\/tr>\n  <\/thead>\n<\/table>","options":{"columnDefs":[{"className":"dt-right","targets":[3,4,5,6,7]},{"orderable":false,"targets":0}],"order":[],"autoWidth":false,"orderClasses":false}},"evals":[],"jsHooks":[]}</script>
```

```r
tmp_df %>% 
  xlsx::write.xlsx(file = out_xlsx,
                   append = TRUE,
                   sheetName = "human2-beta-d3d6",
                   showNA = TRUE)
```


```r
lapply(
  dist,
  FUN = phyloseq_adonis2,
  physeq = ps_tmp_forbdiv %>%
    subset_samples(Day_of_Treatment >= 3 & Day_of_Treatment <= 6),
  formula = paste0(c("Reactor_Treatment_Dose", "Day_of_Treatment"), collapse=" + "),
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
<div class="datatables html-widget html-fill-item-overflow-hidden html-fill-item" id="htmlwidget-0a7605dfcbd2dab75f9a" style="width:100%;height:auto;"></div>
<script type="application/json" data-for="htmlwidget-0a7605dfcbd2dab75f9a">{"x":{"filter":"none","vertical":false,"data":[["1","2","3","4","5","6","7","8"],["bjaccard","bjaccard","bjaccard","bjaccard","aichinson","aichinson","aichinson","aichinson"],["Reactor_Treatment_Dose","Day_of_Treatment","Residual","Total","Reactor_Treatment_Dose","Day_of_Treatment","Residual","Total"],[6,1,19,26,6,1,19,26],[3.918,0.12,1.697,5.736,2.42,0.059,0.687,3.166],[0.683,0.021,0.296,1,0.765,0.019,0.217,1],[7.314,1.34,null,null,11.161,1.642,null,null],[0.001,0.194,null,null,0.001,0.136,null,null]],"container":"<table class=\"display\">\n  <thead>\n    <tr>\n      <th> <\/th>\n      <th>Distance<\/th>\n      <th>terms<\/th>\n      <th>Df<\/th>\n      <th>SumOfSqs<\/th>\n      <th>R2<\/th>\n      <th>F<\/th>\n      <th>Pr(&gt;F)<\/th>\n    <\/tr>\n  <\/thead>\n<\/table>","options":{"columnDefs":[{"className":"dt-right","targets":[3,4,5,6,7]},{"orderable":false,"targets":0}],"order":[],"autoWidth":false,"orderClasses":false}},"evals":[],"jsHooks":[]}</script>
```

```r
tmp_df %>% 
  xlsx::write.xlsx(file = out_xlsx,
                   append = TRUE,
                   sheetName = "human2-beta-d3d6-time-day",
                   showNA = TRUE)
```


```r
lapply(
  dist,
  FUN = physeq_pairwise_permanovas_adonis2,
  physeq = ps_tmp_forbdiv %>%
    subset_samples(Day_of_Treatment >= 3 & Day_of_Treatment <= 6),
  compare_header = "Reactor_Treatment_Dose",
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
<div class="datatables html-widget html-fill-item-overflow-hidden html-fill-item" id="htmlwidget-dd76c3579b426da78523" style="width:100%;height:auto;"></div>
<script type="application/json" data-for="htmlwidget-dd76c3579b426da78523">{"x":{"filter":"none","vertical":false,"data":[["1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","23","24","25","26","27","28","29","30","31","32","33","34","35","36","37","38","39","40","41","42"],["bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson"],["CR_UNTREATED","CR_UNTREATED","CR_UNTREATED","CR_UNTREATED","CR_UNTREATED","CR_UNTREATED","TR3_CTX200","TR3_CTX200","TR3_CTX200","TR3_CTX200","TR3_CTX200","TR5_CTX+HV292.1200","TR5_CTX+HV292.1200","TR5_CTX+HV292.1200","TR5_CTX+HV292.1200","TR7_HV292.1","TR7_HV292.1","TR7_HV292.1","TR2_VAN600","TR2_VAN600","TR4_VAN+CCUG59168600","CR_UNTREATED","CR_UNTREATED","CR_UNTREATED","CR_UNTREATED","CR_UNTREATED","CR_UNTREATED","TR3_CTX200","TR3_CTX200","TR3_CTX200","TR3_CTX200","TR3_CTX200","TR5_CTX+HV292.1200","TR5_CTX+HV292.1200","TR5_CTX+HV292.1200","TR5_CTX+HV292.1200","TR7_HV292.1","TR7_HV292.1","TR7_HV292.1","TR2_VAN600","TR2_VAN600","TR4_VAN+CCUG59168600"],["TR3_CTX200","TR5_CTX+HV292.1200","TR7_HV292.1","TR2_VAN600","TR4_VAN+CCUG59168600","TR6_CCUG59168","TR5_CTX+HV292.1200","TR7_HV292.1","TR2_VAN600","TR4_VAN+CCUG59168600","TR6_CCUG59168","TR7_HV292.1","TR2_VAN600","TR4_VAN+CCUG59168600","TR6_CCUG59168","TR2_VAN600","TR4_VAN+CCUG59168600","TR6_CCUG59168","TR4_VAN+CCUG59168600","TR6_CCUG59168","TR6_CCUG59168","TR3_CTX200","TR5_CTX+HV292.1200","TR7_HV292.1","TR2_VAN600","TR4_VAN+CCUG59168600","TR6_CCUG59168","TR5_CTX+HV292.1200","TR7_HV292.1","TR2_VAN600","TR4_VAN+CCUG59168600","TR6_CCUG59168","TR7_HV292.1","TR2_VAN600","TR4_VAN+CCUG59168600","TR6_CCUG59168","TR2_VAN600","TR4_VAN+CCUG59168600","TR6_CCUG59168","TR4_VAN+CCUG59168600","TR6_CCUG59168","TR6_CCUG59168"],[0.623,0.69,0.339,0.666,0.703,0.322,0.192,0.582,0.357,0.434,0.616,0.632,0.412,0.475,0.674,0.617,0.662,0.38,0.245,0.661,0.704,0.744,0.75,0.377,0.766,0.787,0.34,0.207,0.681,0.437,0.497,0.724,0.681,0.495,0.541,0.725,0.7,0.731,0.381,0.297,0.745,0.773],[0.023,0.028,0.033,0.027,0.032,0.027,0.308,0.024,0.019,0.028,0.04,0.036,0.035,0.034,0.025,0.028,0.029,0.039,0.091,0.03,0.028,0.026,0.028,0.021,0.039,0.038,0.028,0.158,0.036,0.033,0.035,0.037,0.032,0.032,0.029,0.029,0.02,0.04,0.024,0.021,0.034,0.037],[0.242,0.069,0.05,0.103,0.052,0.103,0.308,0.168,0.399,0.069,0.044,0.044,0.046,0.048,0.131,0.069,0.055,0.046,0.096,0.052,0.069,0.109,0.09,0.176,0.043,0.044,0.09,0.158,0.05,0.058,0.053,0.047,0.064,0.064,0.072,0.072,0.42,0.042,0.126,0.176,0.055,0.047]],"container":"<table class=\"display\">\n  <thead>\n    <tr>\n      <th> <\/th>\n      <th>Distance<\/th>\n      <th>X1<\/th>\n      <th>X2<\/th>\n      <th>R2<\/th>\n      <th>pval<\/th>\n      <th>pvalFDR<\/th>\n    <\/tr>\n  <\/thead>\n<\/table>","options":{"columnDefs":[{"className":"dt-right","targets":[4,5,6]},{"orderable":false,"targets":0}],"order":[],"autoWidth":false,"orderClasses":false}},"evals":[],"jsHooks":[]}</script>
```

```r
tmp_df %>% 
  xlsx::write.xlsx(file = out_xlsx,
                   append = TRUE,
                   sheetName = "human2-beta-d3d6-pw",
                   showNA = TRUE)
```


```r
lapply(
  dist,
  FUN = phyloseq_TW,
  physeq =  ps_tmp_forbdiv %>%
    subset_samples(Day_of_Treatment >= 3 & Day_of_Treatment <= 6),
  variable = "Reactor_Treatment_Dose",
  nrep = 999,
  strata = NULL
) %>%
  bind_rows(.id = "Distance") %>%
  mutate_if(is.numeric, round, 3)   -> tmp_df

tmp_df %>% 
  DT::datatable()
```

```{=html}
<div class="datatables html-widget html-fill-item-overflow-hidden html-fill-item" id="htmlwidget-36e9ef7a775d2ff30998" style="width:100%;height:auto;"></div>
<script type="application/json" data-for="htmlwidget-36e9ef7a775d2ff30998">{"x":{"filter":"none","vertical":false,"data":[["1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","23","24","25","26","27","28","29","30","31","32","33","34","35","36","37","38","39","40","41","42"],["bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson"],[["CR_UNTREATED"],["CR_UNTREATED"],["CR_UNTREATED"],["CR_UNTREATED"],["CR_UNTREATED"],["CR_UNTREATED"],["TR3_CTX200"],["TR3_CTX200"],["TR3_CTX200"],["TR3_CTX200"],["TR3_CTX200"],["TR5_CTX+HV292.1200"],["TR5_CTX+HV292.1200"],["TR5_CTX+HV292.1200"],["TR5_CTX+HV292.1200"],["TR7_HV292.1"],["TR7_HV292.1"],["TR7_HV292.1"],["TR2_VAN600"],["TR2_VAN600"],["TR4_VAN+CCUG59168600"],["CR_UNTREATED"],["CR_UNTREATED"],["CR_UNTREATED"],["CR_UNTREATED"],["CR_UNTREATED"],["CR_UNTREATED"],["TR3_CTX200"],["TR3_CTX200"],["TR3_CTX200"],["TR3_CTX200"],["TR3_CTX200"],["TR5_CTX+HV292.1200"],["TR5_CTX+HV292.1200"],["TR5_CTX+HV292.1200"],["TR5_CTX+HV292.1200"],["TR7_HV292.1"],["TR7_HV292.1"],["TR7_HV292.1"],["TR2_VAN600"],["TR2_VAN600"],["TR4_VAN+CCUG59168600"]],[["TR3_CTX200"],["TR5_CTX+HV292.1200"],["TR7_HV292.1"],["TR2_VAN600"],["TR4_VAN+CCUG59168600"],["TR6_CCUG59168"],["TR5_CTX+HV292.1200"],["TR7_HV292.1"],["TR2_VAN600"],["TR4_VAN+CCUG59168600"],["TR6_CCUG59168"],["TR7_HV292.1"],["TR2_VAN600"],["TR4_VAN+CCUG59168600"],["TR6_CCUG59168"],["TR2_VAN600"],["TR4_VAN+CCUG59168600"],["TR6_CCUG59168"],["TR4_VAN+CCUG59168600"],["TR6_CCUG59168"],["TR6_CCUG59168"],["TR3_CTX200"],["TR5_CTX+HV292.1200"],["TR7_HV292.1"],["TR2_VAN600"],["TR4_VAN+CCUG59168600"],["TR6_CCUG59168"],["TR5_CTX+HV292.1200"],["TR7_HV292.1"],["TR2_VAN600"],["TR4_VAN+CCUG59168600"],["TR6_CCUG59168"],["TR7_HV292.1"],["TR2_VAN600"],["TR4_VAN+CCUG59168600"],["TR6_CCUG59168"],["TR2_VAN600"],["TR4_VAN+CCUG59168600"],["TR6_CCUG59168"],["TR4_VAN+CCUG59168600"],["TR6_CCUG59168"],["TR6_CCUG59168"]],[[4],[4],[4],[4],[4],[4],[4],[4],[4],[4],[4],[3],[3],[3],[3],[4],[4],[4],[4],[4],[4],[4],[4],[4],[4],[4],[4],[4],[4],[4],[4],[4],[3],[3],[3],[3],[4],[4],[4],[4],[4],[4]],[[4],[3],[4],[4],[4],[4],[3],[4],[4],[4],[4],[4],[4],[4],[4],[4],[4],[4],[4],[4],[4],[4],[3],[4],[4],[4],[4],[3],[4],[4],[4],[4],[4],[4],[4],[4],[4],[4],[4],[4],[4],[4]],[[0.024],[0.032],[0.021],[0.036],[0.031],[0.029],[0.288],[0.031],[0.03],[0.027],[0.028],[0.041],[0.03],[0.035],[0.03],[0.034],[0.024],[0.036],[0.089],[0.031],[0.024],[0.031],[0.026],[0.033],[0.027],[0.026],[0.033],[0.191],[0.027],[0.038],[0.037],[0.041],[0.032],[0.026],[0.023],[0.033],[0.028],[0.023],[0.029],[0.023],[0.028],[0.036]],[[9.89677388009617],[8.5658519248316],[3.08159303806405],[11.9475301471737],[14.2064439078348],[2.85390968967447],[1.1979936975345],[8.37014947434361],[3.32430843071227],[4.59623239962326],[9.64239644092014],[7.22216869065703],[3.48186730457665],[4.36743225712148],[8.08546591476416],[9.68026306357731],[11.7360553115246],[3.67592081270322],[1.94696423517805],[11.7092466544869],[14.2556631739055],[17.4691923866252],[13.6800612386099],[3.6358648970417],[19.6457448578902],[22.1069690825973],[3.09098714684372],[1.29247142592665],[12.8049404372372],[4.65968244284928],[5.91755945244012],[15.7396533234293],[10.5958996111571],[4.75278030245568],[5.52150482554858],[12.4550299373091],[14.01685270907],[16.3235678992188],[3.69862916015435],[2.52998028873845],[17.5309973493175],[20.4593333429457]],[[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999]],[0.042,0.042,0.042,0.042,0.042,0.042,0.288,0.042,0.042,0.042,0.042,0.045,0.042,0.042,0.042,0.042,0.042,0.042,0.093,0.042,0.042,0.042,0.042,0.042,0.042,0.042,0.042,0.191,0.042,0.042,0.042,0.043,0.042,0.042,0.042,0.042,0.042,0.042,0.042,0.042,0.042,0.042]],"container":"<table class=\"display\">\n  <thead>\n    <tr>\n      <th> <\/th>\n      <th>Distance<\/th>\n      <th>Level1<\/th>\n      <th>Level2<\/th>\n      <th>N1<\/th>\n      <th>N2<\/th>\n      <th>p.value<\/th>\n      <th>tw2.stat<\/th>\n      <th>nrep<\/th>\n      <th>pvalFDR<\/th>\n    <\/tr>\n  <\/thead>\n<\/table>","options":{"columnDefs":[{"className":"dt-right","targets":9},{"orderable":false,"targets":0}],"order":[],"autoWidth":false,"orderClasses":false}},"evals":[],"jsHooks":[]}</script>
```

```r
tmp_df %>% 
  xlsx::write.xlsx(file = out_xlsx,
                   append = TRUE,
                   sheetName = "human2-beta-d3d6-pw2",
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


##### treat - Day_of_Treatment > 7:


```r
lapply(
  dist,
  FUN = phyloseq_adonis2,
  physeq = ps_tmp_forbdiv %>%
    subset_samples(Day_of_Treatment >7),
  formula = paste0(c("Reactor_Treatment_Dose"), collapse=" + "),
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
<div class="datatables html-widget html-fill-item-overflow-hidden html-fill-item" id="htmlwidget-a5f263681883ceb31574" style="width:100%;height:auto;"></div>
<script type="application/json" data-for="htmlwidget-a5f263681883ceb31574">{"x":{"filter":"none","vertical":false,"data":[["1","2","3","4","5","6"],["bjaccard","bjaccard","bjaccard","aichinson","aichinson","aichinson"],["Reactor_Treatment_Dose","Residual","Total","Reactor_Treatment_Dose","Residual","Total"],[6,14,20,6,14,20],[2.803,0.881,3.684,2.014,0.498,2.513],[0.761,0.239,1,0.802,0.198,1],[7.426,null,null,9.437,null,null],[0.001,null,null,0.001,null,null]],"container":"<table class=\"display\">\n  <thead>\n    <tr>\n      <th> <\/th>\n      <th>Distance<\/th>\n      <th>terms<\/th>\n      <th>Df<\/th>\n      <th>SumOfSqs<\/th>\n      <th>R2<\/th>\n      <th>F<\/th>\n      <th>Pr(&gt;F)<\/th>\n    <\/tr>\n  <\/thead>\n<\/table>","options":{"columnDefs":[{"className":"dt-right","targets":[3,4,5,6,7]},{"orderable":false,"targets":0}],"order":[],"autoWidth":false,"orderClasses":false}},"evals":[],"jsHooks":[]}</script>
```

```r
tmp_df %>% 
  xlsx::write.xlsx(file = out_xlsx,
                   append = TRUE,
                   sheetName = "human2-beta-d7p",
                   showNA = TRUE)
```


```r
lapply(
  dist,
  FUN = phyloseq_adonis2,
  physeq = ps_tmp_forbdiv %>%
    subset_samples(Day_of_Treatment >7),
  formula = paste0(c("Reactor_Treatment_Dose", "Day_of_Treatment"), collapse=" + "),
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
<div class="datatables html-widget html-fill-item-overflow-hidden html-fill-item" id="htmlwidget-f017a62ebbbe24ad352e" style="width:100%;height:auto;"></div>
<script type="application/json" data-for="htmlwidget-f017a62ebbbe24ad352e">{"x":{"filter":"none","vertical":false,"data":[["1","2","3","4","5","6","7","8"],["bjaccard","bjaccard","bjaccard","bjaccard","aichinson","aichinson","aichinson","aichinson"],["Reactor_Treatment_Dose","Day_of_Treatment","Residual","Total","Reactor_Treatment_Dose","Day_of_Treatment","Residual","Total"],[6,1,13,20,6,1,13,20],[2.803,0.073,0.808,3.684,2.014,0.04,0.458,2.513],[0.761,0.02,0.219,1,0.802,0.016,0.182,1],[7.518,1.173,null,null,9.534,1.144,null,null],[0.001,0.262,null,null,0.001,0.297,null,null]],"container":"<table class=\"display\">\n  <thead>\n    <tr>\n      <th> <\/th>\n      <th>Distance<\/th>\n      <th>terms<\/th>\n      <th>Df<\/th>\n      <th>SumOfSqs<\/th>\n      <th>R2<\/th>\n      <th>F<\/th>\n      <th>Pr(&gt;F)<\/th>\n    <\/tr>\n  <\/thead>\n<\/table>","options":{"columnDefs":[{"className":"dt-right","targets":[3,4,5,6,7]},{"orderable":false,"targets":0}],"order":[],"autoWidth":false,"orderClasses":false}},"evals":[],"jsHooks":[]}</script>
```

```r
tmp_df %>% 
  xlsx::write.xlsx(file = out_xlsx,
                   append = TRUE,
                   sheetName = "human2-beta-d7p-time-day",
                   showNA = TRUE)
```


```r
lapply(
  dist,
  FUN = physeq_pairwise_permanovas_adonis2,
  physeq = ps_tmp_forbdiv %>%
    subset_samples(Day_of_Treatment > 7),
  compare_header = "Reactor_Treatment_Dose",
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
<div class="datatables html-widget html-fill-item-overflow-hidden html-fill-item" id="htmlwidget-ce8f12c79e792466f32f" style="width:100%;height:auto;"></div>
<script type="application/json" data-for="htmlwidget-ce8f12c79e792466f32f">{"x":{"filter":"none","vertical":false,"data":[["1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","23","24","25","26","27","28","29","30","31","32","33","34","35","36","37","38","39","40","41","42"],["bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson"],["CR_UNTREATED","CR_UNTREATED","CR_UNTREATED","CR_UNTREATED","CR_UNTREATED","CR_UNTREATED","TR3_CTX200","TR3_CTX200","TR3_CTX200","TR3_CTX200","TR3_CTX200","TR5_CTX+HV292.1200","TR5_CTX+HV292.1200","TR5_CTX+HV292.1200","TR5_CTX+HV292.1200","TR7_HV292.1","TR7_HV292.1","TR7_HV292.1","TR2_VAN600","TR2_VAN600","TR4_VAN+CCUG59168600","CR_UNTREATED","CR_UNTREATED","CR_UNTREATED","CR_UNTREATED","CR_UNTREATED","CR_UNTREATED","TR3_CTX200","TR3_CTX200","TR3_CTX200","TR3_CTX200","TR3_CTX200","TR5_CTX+HV292.1200","TR5_CTX+HV292.1200","TR5_CTX+HV292.1200","TR5_CTX+HV292.1200","TR7_HV292.1","TR7_HV292.1","TR7_HV292.1","TR2_VAN600","TR2_VAN600","TR4_VAN+CCUG59168600"],["TR3_CTX200","TR5_CTX+HV292.1200","TR7_HV292.1","TR2_VAN600","TR4_VAN+CCUG59168600","TR6_CCUG59168","TR5_CTX+HV292.1200","TR7_HV292.1","TR2_VAN600","TR4_VAN+CCUG59168600","TR6_CCUG59168","TR7_HV292.1","TR2_VAN600","TR4_VAN+CCUG59168600","TR6_CCUG59168","TR2_VAN600","TR4_VAN+CCUG59168600","TR6_CCUG59168","TR4_VAN+CCUG59168600","TR6_CCUG59168","TR6_CCUG59168","TR3_CTX200","TR5_CTX+HV292.1200","TR7_HV292.1","TR2_VAN600","TR4_VAN+CCUG59168600","TR6_CCUG59168","TR5_CTX+HV292.1200","TR7_HV292.1","TR2_VAN600","TR4_VAN+CCUG59168600","TR6_CCUG59168","TR7_HV292.1","TR2_VAN600","TR4_VAN+CCUG59168600","TR6_CCUG59168","TR2_VAN600","TR4_VAN+CCUG59168600","TR6_CCUG59168","TR4_VAN+CCUG59168600","TR6_CCUG59168","TR6_CCUG59168"],[0.68,0.676,0.523,0.698,0.695,0.502,0.613,0.69,0.45,0.413,0.721,0.693,0.674,0.686,0.763,0.693,0.725,0.482,0.329,0.735,0.757,0.741,0.664,0.531,0.726,0.724,0.484,0.672,0.792,0.537,0.514,0.792,0.71,0.701,0.72,0.731,0.772,0.795,0.475,0.395,0.773,0.796],[0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1],[0.191,0.191,0.191,0.191,0.191,0.191,0.191,0.191,0.191,0.191,0.191,0.191,0.191,0.191,0.191,0.191,0.191,0.191,0.191,0.191,0.191,0.191,0.191,0.191,0.191,0.191,0.191,0.191,0.191,0.191,0.191,0.191,0.191,0.191,0.191,0.191,0.191,0.191,0.191,0.191,0.191,0.191]],"container":"<table class=\"display\">\n  <thead>\n    <tr>\n      <th> <\/th>\n      <th>Distance<\/th>\n      <th>X1<\/th>\n      <th>X2<\/th>\n      <th>R2<\/th>\n      <th>pval<\/th>\n      <th>pvalFDR<\/th>\n    <\/tr>\n  <\/thead>\n<\/table>","options":{"columnDefs":[{"className":"dt-right","targets":[4,5,6]},{"orderable":false,"targets":0}],"order":[],"autoWidth":false,"orderClasses":false}},"evals":[],"jsHooks":[]}</script>
```

```r
tmp_df %>% 
  xlsx::write.xlsx(file = out_xlsx,
                   append = TRUE,
                   sheetName = "human2-beta-d7p-pw",
                   showNA = TRUE)
```


```r
lapply(
  dist,
  FUN = phyloseq_TW,
  physeq =  ps_tmp_forbdiv %>%
    subset_samples(Day_of_Treatment >7),
  variable = "Reactor_Treatment_Dose",
  nrep = 999,
  strata = NULL
) %>%
  bind_rows(.id = "Distance") %>%
  mutate_if(is.numeric, round, 3)   -> tmp_df

tmp_df %>% 
  DT::datatable()
```

```{=html}
<div class="datatables html-widget html-fill-item-overflow-hidden html-fill-item" id="htmlwidget-bbbdcf888a2a729dbcc9" style="width:100%;height:auto;"></div>
<script type="application/json" data-for="htmlwidget-bbbdcf888a2a729dbcc9">{"x":{"filter":"none","vertical":false,"data":[["1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","23","24","25","26","27","28","29","30","31","32","33","34","35","36","37","38","39","40","41","42"],["bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson"],[["CR_UNTREATED"],["CR_UNTREATED"],["CR_UNTREATED"],["CR_UNTREATED"],["CR_UNTREATED"],["CR_UNTREATED"],["TR3_CTX200"],["TR3_CTX200"],["TR3_CTX200"],["TR3_CTX200"],["TR3_CTX200"],["TR5_CTX+HV292.1200"],["TR5_CTX+HV292.1200"],["TR5_CTX+HV292.1200"],["TR5_CTX+HV292.1200"],["TR7_HV292.1"],["TR7_HV292.1"],["TR7_HV292.1"],["TR2_VAN600"],["TR2_VAN600"],["TR4_VAN+CCUG59168600"],["CR_UNTREATED"],["CR_UNTREATED"],["CR_UNTREATED"],["CR_UNTREATED"],["CR_UNTREATED"],["CR_UNTREATED"],["TR3_CTX200"],["TR3_CTX200"],["TR3_CTX200"],["TR3_CTX200"],["TR3_CTX200"],["TR5_CTX+HV292.1200"],["TR5_CTX+HV292.1200"],["TR5_CTX+HV292.1200"],["TR5_CTX+HV292.1200"],["TR7_HV292.1"],["TR7_HV292.1"],["TR7_HV292.1"],["TR2_VAN600"],["TR2_VAN600"],["TR4_VAN+CCUG59168600"]],[["TR3_CTX200"],["TR5_CTX+HV292.1200"],["TR7_HV292.1"],["TR2_VAN600"],["TR4_VAN+CCUG59168600"],["TR6_CCUG59168"],["TR5_CTX+HV292.1200"],["TR7_HV292.1"],["TR2_VAN600"],["TR4_VAN+CCUG59168600"],["TR6_CCUG59168"],["TR7_HV292.1"],["TR2_VAN600"],["TR4_VAN+CCUG59168600"],["TR6_CCUG59168"],["TR2_VAN600"],["TR4_VAN+CCUG59168600"],["TR6_CCUG59168"],["TR4_VAN+CCUG59168600"],["TR6_CCUG59168"],["TR6_CCUG59168"],["TR3_CTX200"],["TR5_CTX+HV292.1200"],["TR7_HV292.1"],["TR2_VAN600"],["TR4_VAN+CCUG59168600"],["TR6_CCUG59168"],["TR5_CTX+HV292.1200"],["TR7_HV292.1"],["TR2_VAN600"],["TR4_VAN+CCUG59168600"],["TR6_CCUG59168"],["TR7_HV292.1"],["TR2_VAN600"],["TR4_VAN+CCUG59168600"],["TR6_CCUG59168"],["TR2_VAN600"],["TR4_VAN+CCUG59168600"],["TR6_CCUG59168"],["TR4_VAN+CCUG59168600"],["TR6_CCUG59168"],["TR6_CCUG59168"]],[[3],[3],[3],[3],[3],[3],[3],[3],[3],[3],[3],[3],[3],[3],[3],[3],[3],[3],[3],[3],[3],[3],[3],[3],[3],[3],[3],[3],[3],[3],[3],[3],[3],[3],[3],[3],[3],[3],[3],[3],[3],[3]],[[3],[3],[3],[3],[3],[3],[3],[3],[3],[3],[3],[3],[3],[3],[3],[3],[3],[3],[3],[3],[3],[3],[3],[3],[3],[3],[3],[3],[3],[3],[3],[3],[3],[3],[3],[3],[3],[3],[3],[3],[3],[3]],[[0.09],[0.102],[0.118],[0.083],[0.091],[0.105],[0.098],[0.087],[0.098],[0.095],[0.1],[0.088],[0.086],[0.121],[0.108],[0.097],[0.099],[0.093],[0.091],[0.095],[0.095],[0.09],[0.098],[0.117],[0.102],[0.101],[0.098],[0.112],[0.098],[0.088],[0.097],[0.127],[0.092],[0.101],[0.1],[0.096],[0.1],[0.095],[0.107],[0.105],[0.112],[0.094]],[[8.51136481272645],[8.32858669127167],[4.38458967318319],[9.24219011022793],[9.10864095999458],[4.03157373721911],[6.34525007796114],[8.88928491635667],[3.26800903195742],[2.81384306547469],[10.3309608167587],[9.01283323313974],[8.26874761128541],[8.73296061639032],[12.8746231612236],[9.04457343768923],[10.5233345011014],[3.72635693968114],[1.96456105501301],[11.1199345835546],[12.4804925197844],[11.451869296465],[7.90403726500688],[4.52894670820265],[10.622632638945],[10.5147755682571],[3.75608382946069],[8.20470382025644],[15.2092073550007],[4.63515443823741],[4.22380213898984],[15.1875260663548],[9.80565859255954],[9.35574649953762],[10.2754455533249],[10.8450294458474],[13.5067869819737],[15.5552227012796],[3.61436795606982],[2.6132149043023],[13.6193673807053],[15.5735960618277]],[[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999]],[0.119,0.119,0.121,0.119,0.119,0.119,0.119,0.119,0.119,0.119,0.119,0.119,0.119,0.121,0.119,0.119,0.119,0.119,0.119,0.119,0.119,0.123,0.123,0.123,0.123,0.123,0.123,0.123,0.123,0.123,0.123,0.127,0.123,0.123,0.123,0.123,0.123,0.123,0.123,0.123,0.123,0.123]],"container":"<table class=\"display\">\n  <thead>\n    <tr>\n      <th> <\/th>\n      <th>Distance<\/th>\n      <th>Level1<\/th>\n      <th>Level2<\/th>\n      <th>N1<\/th>\n      <th>N2<\/th>\n      <th>p.value<\/th>\n      <th>tw2.stat<\/th>\n      <th>nrep<\/th>\n      <th>pvalFDR<\/th>\n    <\/tr>\n  <\/thead>\n<\/table>","options":{"columnDefs":[{"className":"dt-right","targets":9},{"orderable":false,"targets":0}],"order":[],"autoWidth":false,"orderClasses":false}},"evals":[],"jsHooks":[]}</script>
```

```r
tmp_df %>% 
  xlsx::write.xlsx(file = out_xlsx,
                   append = TRUE,
                   sheetName = "human2-beta-d7p-pw2",
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


##### Period:


```r
setNames(vector("list", 
                length(physeq_get_unique(ps_tmp_forbdiv, "Reactor_Treatment_Dose"))), 
         physeq_get_unique(ps_tmp_forbdiv, "Reactor_Treatment_Dose")) -> lis
```



```r
for (grp in names(lis)){
  print(grp)
  
  lapply(
    dist,
    FUN = physeq_pairwise_permanovas_adonis2,
    physeq = ps_tmp_forbdiv %>%
      # subset_samples(.,  eval(as.name(ranks)) %in% tax) %>% 
      subset_samples(., Reactor_Treatment_Dose == grp),
    compare_header = "period",
    n_perm = 999,
    strat = "none"
  ) %>%
    bind_rows(.id = "Distance") %>%
    mutate_if(is.numeric, round, 3) %>%
    select(-pvalBon) -> lis[[grp]]
}
```

```
## [1] "TR5_CTX+HV292.1200"
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
## [1] "TR3_CTX200"
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
```

```
## [1] "TR4_VAN+CCUG59168600"
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
```

```
## [1] "TR2_VAN600"
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
```

```
## [1] "TR7_HV292.1"
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
```

```
## [1] "CR_UNTREATED"
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
```

```
## [1] "TR6_CCUG59168"
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
```

```r
lis %>% 
  bind_rows(.id = "Reactor_Treatment_Dose") %>% 
  filter(X1 == "pret" | "X2" == "pret") %>% 
  group_by(Reactor_Treatment_Dose, Distance) %>% 
  rstatix::adjust_pvalue(p.col = "pval", output.col = "pvalFDR") %>% 
  rstatix::add_significance(p.col = "pvalFDR")-> tmp_df

tmp_df %>% 
  DT::datatable()
```

```{=html}
<div class="datatables html-widget html-fill-item-overflow-hidden html-fill-item" id="htmlwidget-36bca1387d25533ead4d" style="width:100%;height:auto;"></div>
<script type="application/json" data-for="htmlwidget-36bca1387d25533ead4d">{"x":{"filter":"none","vertical":false,"data":[["1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","23","24","25","26","27","28","29","30","31","32","33","34","35","36","37","38","39","40","41","42"],["CR_UNTREATED","CR_UNTREATED","CR_UNTREATED","CR_UNTREATED","CR_UNTREATED","CR_UNTREATED","TR2_VAN600","TR2_VAN600","TR2_VAN600","TR2_VAN600","TR2_VAN600","TR2_VAN600","TR3_CTX200","TR3_CTX200","TR3_CTX200","TR3_CTX200","TR3_CTX200","TR3_CTX200","TR4_VAN+CCUG59168600","TR4_VAN+CCUG59168600","TR4_VAN+CCUG59168600","TR4_VAN+CCUG59168600","TR4_VAN+CCUG59168600","TR4_VAN+CCUG59168600","TR5_CTX+HV292.1200","TR5_CTX+HV292.1200","TR5_CTX+HV292.1200","TR5_CTX+HV292.1200","TR5_CTX+HV292.1200","TR5_CTX+HV292.1200","TR6_CCUG59168","TR6_CCUG59168","TR6_CCUG59168","TR6_CCUG59168","TR6_CCUG59168","TR6_CCUG59168","TR7_HV292.1","TR7_HV292.1","TR7_HV292.1","TR7_HV292.1","TR7_HV292.1","TR7_HV292.1"],["aichinson","aichinson","aichinson","bjaccard","bjaccard","bjaccard","aichinson","aichinson","aichinson","bjaccard","bjaccard","bjaccard","aichinson","aichinson","aichinson","bjaccard","bjaccard","bjaccard","aichinson","aichinson","aichinson","bjaccard","bjaccard","bjaccard","aichinson","aichinson","aichinson","bjaccard","bjaccard","bjaccard","aichinson","aichinson","aichinson","bjaccard","bjaccard","bjaccard","aichinson","aichinson","aichinson","bjaccard","bjaccard","bjaccard"],["pret","pret","pret","pret","pret","pret","pret","pret","pret","pret","pret","pret","pret","pret","pret","pret","pret","pret","pret","pret","pret","pret","pret","pret","pret","pret","pret","pret","pret","pret","pret","pret","pret","pret","pret","pret","pret","pret","pret","pret","pret","pret"],["t1","t2","t3","t1","t2","t3","t1","t2","t3","t1","t2","t3","t1","t2","t3","t1","t2","t3","t1","t2","t3","t1","t2","t3","t1","t2","t3","t1","t2","t3","t1","t2","t3","t1","t2","t3","t1","t2","t3","t1","t2","t3"],[0.215,0.251,0.337,0.125,0.223,0.324,0.437,0.693,0.646,0.429,0.647,0.612,0.477,0.698,0.706,0.566,0.63,0.642,0.438,0.724,0.69,0.372,0.672,0.673,0.479,0.704,0.593,0.516,0.686,0.6,0.177,0.223,0.364,0.222,0.233,0.4,0.19,0.329,0.469,0.193,0.277,0.471],[0.093,0.009,0.008,0.764,0.019,0.011,0.053,0.007,0.012,0.045,0.005,0.012,0.049,0.01,0.011,0.047,0.013,0.014,0.047,0.007,0.009,0.098,0.008,0.011,0.047,0.018,0.01,0.056,0.015,0.005,0.508,0.008,0.014,0.133,0.023,0.013,0.152,0.016,0.011,0.442,0.01,0.009],[0.093,0.024,0.024,0.764,0.038,0.033,0.053,0.021,0.024,0.045,0.015,0.024,0.049,0.03,0.03,0.047,0.039,0.039,0.047,0.021,0.021,0.098,0.024,0.024,0.047,0.036,0.03,0.056,0.03,0.015,0.508,0.024,0.028,0.133,0.046,0.039,0.152,0.033,0.033,0.442,0.027,0.027],["ns","*","*","ns","*","*","ns","*","*","*","*","*","*","*","*","*","*","*","*","*","*","ns","*","*","*","*","*","ns","*","*","ns","*","*","ns","*","*","ns","*","*","ns","*","*"]],"container":"<table class=\"display\">\n  <thead>\n    <tr>\n      <th> <\/th>\n      <th>Reactor_Treatment_Dose<\/th>\n      <th>Distance<\/th>\n      <th>X1<\/th>\n      <th>X2<\/th>\n      <th>R2<\/th>\n      <th>pval<\/th>\n      <th>pvalFDR<\/th>\n      <th>pvalFDR.signif<\/th>\n    <\/tr>\n  <\/thead>\n<\/table>","options":{"columnDefs":[{"className":"dt-right","targets":[5,6,7]},{"orderable":false,"targets":0}],"order":[],"autoWidth":false,"orderClasses":false}},"evals":[],"jsHooks":[]}</script>
```

```r
tmp_df %>% 
  xlsx::write.xlsx(file = out_xlsx,
                   append = TRUE,
                   sheetName = "human2-beta-vs_pret",
                   showNA = TRUE)
```



```r
lis %>% 
  bind_rows(.id = "Reactor_Treatment_Dose") %>% 
  group_by(Reactor_Treatment_Dose, Distance) %>% 
  rstatix::adjust_pvalue(p.col = "pval", output.col = "pvalFDR") %>% 
  rstatix::add_significance(p.col = "pvalFDR") -> tmp_df

tmp_df %>% 
  DT::datatable()
```

```{=html}
<div class="datatables html-widget html-fill-item-overflow-hidden html-fill-item" id="htmlwidget-45918177c74a14b8c6f2" style="width:100%;height:auto;"></div>
<script type="application/json" data-for="htmlwidget-45918177c74a14b8c6f2">{"x":{"filter":"none","vertical":false,"data":[["1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","23","24","25","26","27","28","29","30","31","32","33","34","35","36","37","38","39","40","41","42","43","44","45","46","47","48","49","50","51","52","53","54","55","56","57","58","59","60","61","62","63","64","65","66","67","68","69","70","71","72","73","74","75","76","77","78","79","80","81","82","83","84"],["CR_UNTREATED","CR_UNTREATED","CR_UNTREATED","CR_UNTREATED","CR_UNTREATED","CR_UNTREATED","CR_UNTREATED","CR_UNTREATED","CR_UNTREATED","CR_UNTREATED","CR_UNTREATED","CR_UNTREATED","TR2_VAN600","TR2_VAN600","TR2_VAN600","TR2_VAN600","TR2_VAN600","TR2_VAN600","TR2_VAN600","TR2_VAN600","TR2_VAN600","TR2_VAN600","TR2_VAN600","TR2_VAN600","TR3_CTX200","TR3_CTX200","TR3_CTX200","TR3_CTX200","TR3_CTX200","TR3_CTX200","TR3_CTX200","TR3_CTX200","TR3_CTX200","TR3_CTX200","TR3_CTX200","TR3_CTX200","TR4_VAN+CCUG59168600","TR4_VAN+CCUG59168600","TR4_VAN+CCUG59168600","TR4_VAN+CCUG59168600","TR4_VAN+CCUG59168600","TR4_VAN+CCUG59168600","TR4_VAN+CCUG59168600","TR4_VAN+CCUG59168600","TR4_VAN+CCUG59168600","TR4_VAN+CCUG59168600","TR4_VAN+CCUG59168600","TR4_VAN+CCUG59168600","TR5_CTX+HV292.1200","TR5_CTX+HV292.1200","TR5_CTX+HV292.1200","TR5_CTX+HV292.1200","TR5_CTX+HV292.1200","TR5_CTX+HV292.1200","TR5_CTX+HV292.1200","TR5_CTX+HV292.1200","TR5_CTX+HV292.1200","TR5_CTX+HV292.1200","TR5_CTX+HV292.1200","TR5_CTX+HV292.1200","TR6_CCUG59168","TR6_CCUG59168","TR6_CCUG59168","TR6_CCUG59168","TR6_CCUG59168","TR6_CCUG59168","TR6_CCUG59168","TR6_CCUG59168","TR6_CCUG59168","TR6_CCUG59168","TR6_CCUG59168","TR6_CCUG59168","TR7_HV292.1","TR7_HV292.1","TR7_HV292.1","TR7_HV292.1","TR7_HV292.1","TR7_HV292.1","TR7_HV292.1","TR7_HV292.1","TR7_HV292.1","TR7_HV292.1","TR7_HV292.1","TR7_HV292.1"],["aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard"],["pret","pret","pret","t1","t1","t2","pret","pret","pret","t1","t1","t2","pret","pret","pret","t1","t1","t2","pret","pret","pret","t1","t1","t2","pret","pret","pret","t1","t1","t2","pret","pret","pret","t1","t1","t2","pret","pret","pret","t1","t1","t2","pret","pret","pret","t1","t1","t2","pret","pret","pret","t1","t1","t2","pret","pret","pret","t1","t1","t2","pret","pret","pret","t1","t1","t2","pret","pret","pret","t1","t1","t2","pret","pret","pret","t1","t1","t2","pret","pret","pret","t1","t1","t2"],["t1","t2","t3","t2","t3","t3","t1","t2","t3","t2","t3","t3","t1","t2","t3","t2","t3","t3","t1","t2","t3","t2","t3","t3","t1","t2","t3","t2","t3","t3","t1","t2","t3","t2","t3","t3","t1","t2","t3","t2","t3","t3","t1","t2","t3","t2","t3","t3","t1","t2","t3","t2","t3","t3","t1","t2","t3","t2","t3","t3","t1","t2","t3","t2","t3","t3","t1","t2","t3","t2","t3","t3","t1","t2","t3","t2","t3","t3","t1","t2","t3","t2","t3","t3"],[0.215,0.251,0.337,0.229,0.299,0.207,0.125,0.223,0.324,0.234,0.271,0.193,0.437,0.693,0.646,0.336,0.342,0.218,0.429,0.647,0.612,0.26,0.275,0.202,0.477,0.698,0.706,0.429,0.447,0.182,0.566,0.63,0.642,0.346,0.377,0.159,0.438,0.724,0.69,0.454,0.425,0.227,0.372,0.672,0.673,0.396,0.389,0.261,0.479,0.704,0.593,0.502,0.431,0.465,0.516,0.686,0.6,0.472,0.422,0.453,0.177,0.223,0.364,0.214,0.433,0.231,0.222,0.233,0.4,0.248,0.498,0.205,0.19,0.329,0.469,0.359,0.504,0.3,0.193,0.277,0.471,0.284,0.51,0.27],[0.093,0.009,0.008,0.2,0.067,0.028,0.764,0.019,0.011,0.267,0.133,0.117,0.053,0.007,0.012,0.067,0.067,0.092,0.045,0.005,0.012,0.2,0.2,0.138,0.049,0.01,0.011,0.067,0.067,0.04,0.047,0.013,0.014,0.067,0.067,0.356,0.047,0.007,0.009,0.067,0.067,0.037,0.098,0.008,0.011,0.067,0.067,0.059,0.047,0.018,0.01,0.1,0.067,0.064,0.056,0.015,0.005,0.1,0.067,0.021,0.508,0.008,0.014,0.333,0.067,0.03,0.133,0.023,0.013,0.267,0.067,0.07,0.152,0.016,0.011,0.067,0.067,0.026,0.442,0.01,0.009,0.067,0.067,0.03],[0.201,0.048,0.048,0.201,0.201,0.112,0.764,0.095,0.066,0.534,0.468,0.468,0.212,0.042,0.06,0.212,0.212,0.212,0.18,0.03,0.06,0.414,0.414,0.414,0.16,0.06,0.06,0.16,0.16,0.16,0.188,0.078,0.078,0.201,0.201,0.356,0.148,0.042,0.045,0.148,0.148,0.148,0.236,0.048,0.055,0.236,0.236,0.236,0.188,0.09,0.06,0.192,0.192,0.192,0.168,0.075,0.03,0.168,0.168,0.084,0.666,0.048,0.07,0.666,0.201,0.12,0.268,0.115,0.078,0.268,0.268,0.268,0.201,0.08,0.066,0.201,0.201,0.104,0.442,0.054,0.054,0.201,0.201,0.12],["ns","*","*","ns","ns","ns","ns","ns","ns","ns","ns","ns","ns","*","ns","ns","ns","ns","ns","*","ns","ns","ns","ns","ns","ns","ns","ns","ns","ns","ns","ns","ns","ns","ns","ns","ns","*","*","ns","ns","ns","ns","*","ns","ns","ns","ns","ns","ns","ns","ns","ns","ns","ns","ns","*","ns","ns","ns","ns","*","ns","ns","ns","ns","ns","ns","ns","ns","ns","ns","ns","ns","ns","ns","ns","ns","ns","ns","ns","ns","ns","ns"]],"container":"<table class=\"display\">\n  <thead>\n    <tr>\n      <th> <\/th>\n      <th>Reactor_Treatment_Dose<\/th>\n      <th>Distance<\/th>\n      <th>X1<\/th>\n      <th>X2<\/th>\n      <th>R2<\/th>\n      <th>pval<\/th>\n      <th>pvalFDR<\/th>\n      <th>pvalFDR.signif<\/th>\n    <\/tr>\n  <\/thead>\n<\/table>","options":{"columnDefs":[{"className":"dt-right","targets":[5,6,7]},{"orderable":false,"targets":0}],"order":[],"autoWidth":false,"orderClasses":false}},"evals":[],"jsHooks":[]}</script>
```

```r
tmp_df %>% 
  xlsx::write.xlsx(file = out_xlsx,
                   append = TRUE,
                   sheetName = "human2-beta-between-periods",
                   showNA = TRUE)
```



```r
# for (grp in names(lis)){
#   print(grp)
#   
# lapply(
#   dist,
#   FUN = phyloseq_TW,
#   physeq =  ps_tmp_forbdiv %>%
#       # subset_samples(.,  eval(as.name(ranks)) %in% tax) %>% 
#       subset_samples(., Reactor_Treatment_Dose == grp),
#     variable = "period",
#   nrep = 999,
#   strata = NULL
# ) %>%
#   bind_rows(.id = "Distance") %>%
#   mutate_if(is.numeric, round, 5) -> lis[[grp]]
# }
# 
# 
# lis %>% 
#   bind_rows(.id = "Reactor_Treatment_Dose") %>% 
#   xlsx::write.xlsx(file = out_xlsx,
#                    append = TRUE,
#                    sheetName = "human2-beta-between-periodsTW",
#                    showNA = TRUE)
```





### Human1 - all:



```r
ps_rare %>%
  subset_samples(Model2 == "Human1"  & Treatment != "DONOR"  & Fermentation %in% c(1,2)) -> ps_tmp_forbdiv
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

![](NRP72_draft_16S_Figs_beta_PCoA_files/figure-html/unnamed-chunk-125-1.png)<!-- -->


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
<div class="datatables html-widget html-fill-item-overflow-hidden html-fill-item" id="htmlwidget-56473ac9a4503685a7f3" style="width:100%;height:auto;"></div>
<script type="application/json" data-for="htmlwidget-56473ac9a4503685a7f3">{"x":{"filter":"none","vertical":false,"data":[["ASV0002","ASV0051","ASV0010","ASV0015","ASV0019","ASV0035","ASV0106","ASV0124"],[-0.597834179220655,0.657937890091773,-0.875873134371687,0.78955848343164,0.506301477829794,0.538606751864561,-0.741001013194891,-0.645690741087984],[-0.332459687317544,0.334520008130732,0.214953317291109,0.128341752631704,0.476402252087666,0.510887546877555,-0.0171919400499977,0.384839333628445],[0.001,0.001,0.001,0.001,0.001,0.001,0.001,0.001],[0.467935149535713,0.544785903058199,0.813358676128535,0.639874204227448,0.483300292246834,0.551103318708659,0.549378064358538,0.565017845834336],["Bacteria","Bacteria","Bacteria","Bacteria","Bacteria","Bacteria","Bacteria","Bacteria"],["Bacteroidota","Firmicutes","Firmicutes","Firmicutes","Proteobacteria","Proteobacteria","Firmicutes","Firmicutes"],["Bacteroidia","Bacilli","Clostridia","Negativicutes","Gammaproteobacteria","Gammaproteobacteria","Clostridia","Clostridia"],["Bacteroidales","Lactobacillales","Lachnospirales","Veillonellales-Selenomonadales","Pseudomonadales","Burkholderiales","Oscillospirales","Monoglobales"],["Bacteroidaceae","Lactobacillaceae","Lachnospiraceae","Veillonellaceae","Pseudomonadaceae","Sutterellaceae","Oscillospiraceae","Monoglobaceae"],[null,null,null,null,null,null,null,null],[null,null,null,null,null,null,null,null],[null,null,null,null,null,null,null,null],[0.0026,0.0026,0.0026,0.0026,0.0026,0.0026,0.0026,0.0026]],"container":"<table class=\"display\">\n  <thead>\n    <tr>\n      <th> <\/th>\n      <th>Axis.1<\/th>\n      <th>Axis.2<\/th>\n      <th>pval<\/th>\n      <th>r<\/th>\n      <th>Kingdom<\/th>\n      <th>Phylum<\/th>\n      <th>Class<\/th>\n      <th>Order<\/th>\n      <th>tax_rank_plot<\/th>\n      <th>Genus<\/th>\n      <th>Species<\/th>\n      <th>Strain<\/th>\n      <th>pval.adj<\/th>\n    <\/tr>\n  <\/thead>\n<\/table>","options":{"columnDefs":[{"className":"dt-right","targets":[1,2,3,4,13]},{"orderable":false,"targets":0}],"order":[],"autoWidth":false,"orderClasses":false}},"evals":[],"jsHooks":[]}</script>
```


```r
out$plot -> PCoA_tmp2_fam

PCoA_tmp2 + 
  theme(legend.position = "none") -> PCoA_tmp2_noleg

PCoA_tmp2_noleg
```

![](NRP72_draft_16S_Figs_beta_PCoA_files/figure-html/unnamed-chunk-127-1.png)<!-- -->

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
<div class="datatables html-widget html-fill-item-overflow-hidden html-fill-item" id="htmlwidget-057f6e9b1f30e2dbee53" style="width:100%;height:auto;"></div>
<script type="application/json" data-for="htmlwidget-057f6e9b1f30e2dbee53">{"x":{"filter":"none","vertical":false,"data":[["1","2","3","4","5","6"],["Formiat_mM","Acetat_mM","Propionat_mM","Butyrat_mM","Valerat_mM","Total_SCFA_mM"],[0.270162708759593,-0.878418575751238,-0.0488725590367165,-0.794050468978855,-0.117619022423589,-0.78010067992651],[0.551509577899249,-0.234172870112844,-0.302091816938156,-0.454780783257901,-0.678509177560713,-0.433984784767189],["Formiat_mM","Acetat_mM","Propionat_mM","Butyrat_mM","Valerat_mM","Total_SCFA_mM"],[0.001,0.001,0.009,0.001,0.001,0.001],[0.377150703718928,0.82645612732172,0.0936479928877934,0.837341708106209,0.474208938469996,0.796899864231227]],"container":"<table class=\"display\">\n  <thead>\n    <tr>\n      <th> <\/th>\n      <th>label<\/th>\n      <th>Axis.1<\/th>\n      <th>Axis.2<\/th>\n      <th>id<\/th>\n      <th>pval<\/th>\n      <th>r<\/th>\n    <\/tr>\n  <\/thead>\n<\/table>","options":{"columnDefs":[{"className":"dt-right","targets":[2,3,5,6]},{"orderable":false,"targets":0}],"order":[],"autoWidth":false,"orderClasses":false}},"evals":[],"jsHooks":[]}</script>
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

![](NRP72_draft_16S_Figs_beta_PCoA_files/figure-html/unnamed-chunk-130-1.png)<!-- -->


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
<div class="datatables html-widget html-fill-item-overflow-hidden html-fill-item" id="htmlwidget-c0c4ca2e9381bd6fbebd" style="width:100%;height:auto;"></div>
<script type="application/json" data-for="htmlwidget-c0c4ca2e9381bd6fbebd">{"x":{"filter":"none","vertical":false,"data":[["ASV0002","ASV0012","ASV0051","ASV0010","ASV0015","ASV0035","ASV0106","ASV0124"],[-0.637590724142918,-0.429936467675664,0.700526550112209,-0.776746308038715,0.761430656215353,0.555743853351385,-0.689019065137533,-0.539712190289124],[-0.206234890986704,-0.509585879553196,0.265871052745876,0.404433893018545,-0.110706612899521,0.408775173606298,0.194157556727153,0.486883130283422],[0.001,0.001,0.001,0.001,0.001,0.001,0.001,0.001],[0.449054761773388,0.444523134877431,0.561424864100314,0.766901600873911,0.592032598364227,0.475948373094705,0.512444428957257,0.528344430901268],["Bacteria","Bacteria","Bacteria","Bacteria","Bacteria","Bacteria","Bacteria","Bacteria"],["Bacteroidota","Firmicutes","Firmicutes","Firmicutes","Firmicutes","Proteobacteria","Firmicutes","Firmicutes"],["Bacteroidia","Clostridia","Bacilli","Clostridia","Negativicutes","Gammaproteobacteria","Clostridia","Clostridia"],["Bacteroidales","Oscillospirales","Lactobacillales","Lachnospirales","Veillonellales-Selenomonadales","Burkholderiales","Oscillospirales","Monoglobales"],["Bacteroidaceae","Ruminococcaceae","Lactobacillaceae","Lachnospiraceae","Veillonellaceae","Sutterellaceae","Oscillospiraceae","Monoglobaceae"],[null,null,null,null,null,null,null,null],[null,null,null,null,null,null,null,null],[null,null,null,null,null,null,null,null],[0.00273684210526316,0.00273684210526316,0.00273684210526316,0.00273684210526316,0.00273684210526316,0.00273684210526316,0.00273684210526316,0.00273684210526316]],"container":"<table class=\"display\">\n  <thead>\n    <tr>\n      <th> <\/th>\n      <th>Axis.1<\/th>\n      <th>Axis.2<\/th>\n      <th>pval<\/th>\n      <th>r<\/th>\n      <th>Kingdom<\/th>\n      <th>Phylum<\/th>\n      <th>Class<\/th>\n      <th>Order<\/th>\n      <th>tax_rank_plot<\/th>\n      <th>Genus<\/th>\n      <th>Species<\/th>\n      <th>Strain<\/th>\n      <th>pval.adj<\/th>\n    <\/tr>\n  <\/thead>\n<\/table>","options":{"columnDefs":[{"className":"dt-right","targets":[1,2,3,4,13]},{"orderable":false,"targets":0}],"order":[],"autoWidth":false,"orderClasses":false}},"evals":[],"jsHooks":[]}</script>
```


```r
out$plot -> PCoA_tmp2_fam

PCoA_tmp2 + 
  theme(legend.position = "none") -> PCoA_tmp2_noleg

PCoA_tmp2_noleg
```

![](NRP72_draft_16S_Figs_beta_PCoA_files/figure-html/unnamed-chunk-132-1.png)<!-- -->

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
<div class="datatables html-widget html-fill-item-overflow-hidden html-fill-item" id="htmlwidget-767209723a30863ca56b" style="width:100%;height:auto;"></div>
<script type="application/json" data-for="htmlwidget-767209723a30863ca56b">{"x":{"filter":"none","vertical":false,"data":[["1","2","3","4","5"],["Formiat_mM","Acetat_mM","Butyrat_mM","Valerat_mM","Total_SCFA_mM"],[0.283840135282746,-0.88494812542876,-0.824562216745593,-0.123052332607174,-0.813781416353768],[0.547269594930534,-0.0618837504118045,-0.277260208358263,-0.583396578115785,-0.249053480771384],["Formiat_mM","Acetat_mM","Butyrat_mM","Valerat_mM","Total_SCFA_mM"],[0.001,0.001,0.001,0.001,0.001],[0.380069231932758,0.786962783264907,0.756776072423274,0.355493443917274,0.724267829887086]],"container":"<table class=\"display\">\n  <thead>\n    <tr>\n      <th> <\/th>\n      <th>label<\/th>\n      <th>Axis.1<\/th>\n      <th>Axis.2<\/th>\n      <th>id<\/th>\n      <th>pval<\/th>\n      <th>r<\/th>\n    <\/tr>\n  <\/thead>\n<\/table>","options":{"columnDefs":[{"className":"dt-right","targets":[2,3,5,6]},{"orderable":false,"targets":0}],"order":[],"autoWidth":false,"orderClasses":false}},"evals":[],"jsHooks":[]}</script>
```


```r
env_out$plot %>% 
  export::graph2ppt(append = TRUE,
                    width = 317.48031496 * 1,
                    height = 0.618 * 317.48031496 * 1 , paper = "A4",  scaling = 2,
                    file = out_pptx)
```


####  Statistical evaluation:

##### pret - Day_of_Treatment <= 0:


```r
lapply(
  dist,
  FUN = phyloseq_adonis2,
  physeq = ps_tmp_forbdiv %>%
    subset_samples(Day_of_Treatment <= 0),
  formula = paste0(c("Reactor_Treatment_Dose"), collapse=" + "),
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
<div class="datatables html-widget html-fill-item-overflow-hidden html-fill-item" id="htmlwidget-23741672a15783b59cb9" style="width:100%;height:auto;"></div>
<script type="application/json" data-for="htmlwidget-23741672a15783b59cb9">{"x":{"filter":"none","vertical":false,"data":[["1","2","3","4","5","6"],["bjaccard","bjaccard","bjaccard","aichinson","aichinson","aichinson"],["Reactor_Treatment_Dose","Residual","Total","Reactor_Treatment_Dose","Residual","Total"],[10,24,34,10,24,34],[1.103,0.732,1.835,1.569,1.045,2.614],[0.601,0.399,1,0.6,0.4,1],[3.618,null,null,3.601,null,null],[0.001,null,null,0.001,null,null]],"container":"<table class=\"display\">\n  <thead>\n    <tr>\n      <th> <\/th>\n      <th>Distance<\/th>\n      <th>terms<\/th>\n      <th>Df<\/th>\n      <th>SumOfSqs<\/th>\n      <th>R2<\/th>\n      <th>F<\/th>\n      <th>Pr(&gt;F)<\/th>\n    <\/tr>\n  <\/thead>\n<\/table>","options":{"columnDefs":[{"className":"dt-right","targets":[3,4,5,6,7]},{"orderable":false,"targets":0}],"order":[],"autoWidth":false,"orderClasses":false}},"evals":[],"jsHooks":[]}</script>
```

```r
tmp_df %>% 
  xlsx::write.xlsx(file = out_xlsx,
                   append = TRUE,
                   sheetName = "human1all-beta-pretreat",
                   showNA = TRUE)
```



```r
lapply(
  dist,
  FUN = phyloseq_adonis2,
  physeq = ps_tmp_forbdiv %>%
    subset_samples(Day_of_Treatment <= 0),
  formula = paste0(c("Reactor_Treatment_Dose", "Day_of_Treatment"), collapse=" + "),
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
<div class="datatables html-widget html-fill-item-overflow-hidden html-fill-item" id="htmlwidget-d6facaf5e0f5b7b2b1fa" style="width:100%;height:auto;"></div>
<script type="application/json" data-for="htmlwidget-d6facaf5e0f5b7b2b1fa">{"x":{"filter":"none","vertical":false,"data":[["1","2","3","4","5","6","7","8"],["bjaccard","bjaccard","bjaccard","bjaccard","aichinson","aichinson","aichinson","aichinson"],["Reactor_Treatment_Dose","Day_of_Treatment","Residual","Total","Reactor_Treatment_Dose","Day_of_Treatment","Residual","Total"],[10,1,23,34,10,1,23,34],[1.091,0.032,0.7,1.835,1.559,0.043,1.002,2.614],[0.594,0.017,0.382,1,0.596,0.017,0.383,1],[3.582,1.035,null,null,3.579,0.997,null,null],[0.001,0.371,null,null,0.001,0.424,null,null]],"container":"<table class=\"display\">\n  <thead>\n    <tr>\n      <th> <\/th>\n      <th>Distance<\/th>\n      <th>terms<\/th>\n      <th>Df<\/th>\n      <th>SumOfSqs<\/th>\n      <th>R2<\/th>\n      <th>F<\/th>\n      <th>Pr(&gt;F)<\/th>\n    <\/tr>\n  <\/thead>\n<\/table>","options":{"columnDefs":[{"className":"dt-right","targets":[3,4,5,6,7]},{"orderable":false,"targets":0}],"order":[],"autoWidth":false,"orderClasses":false}},"evals":[],"jsHooks":[]}</script>
```

```r
tmp_df %>% 
  xlsx::write.xlsx(file = out_xlsx,
                   append = TRUE,
                   sheetName = "human1all-beta-pretreat-time-day",
                   showNA = TRUE)
```


```r
lapply(
  dist,
  FUN = physeq_pairwise_permanovas_adonis2,
  physeq = ps_tmp_forbdiv %>%
    subset_samples(Day_of_Treatment <= 0),
  compare_header = "Reactor_Treatment_Dose",
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
## Set of permutations < 'minperm'. Generating entire set.
```

```r
tmp_df %>% 
  DT::datatable()
```

```{=html}
<div class="datatables html-widget html-fill-item-overflow-hidden html-fill-item" id="htmlwidget-2f8220dda1413fa7b913" style="width:100%;height:auto;"></div>
<script type="application/json" data-for="htmlwidget-2f8220dda1413fa7b913">{"x":{"filter":"none","vertical":false,"data":[["1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","23","24","25","26","27","28","29","30","31","32","33","34","35","36","37","38","39","40","41","42","43","44","45","46","47","48","49","50","51","52","53","54","55","56","57","58","59","60","61","62","63","64","65","66","67","68","69","70","71","72","73","74","75","76","77","78","79","80","81","82","83","84","85","86","87","88","89","90","91","92","93","94","95","96","97","98","99","100","101","102","103","104","105","106","107","108","109","110"],["bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson"],["CR_UNTREATED","CR_UNTREATED","CR_UNTREATED","CR_UNTREATED","CR_UNTREATED","CR_UNTREATED","CR_UNTREATED","CR_UNTREATED","CR_UNTREATED","CR_UNTREATED","TR2_CTX20","TR2_CTX20","TR2_CTX20","TR2_CTX20","TR2_CTX20","TR2_CTX20","TR2_CTX20","TR2_CTX20","TR2_CTX20","TR4_CTX200","TR4_CTX200","TR4_CTX200","TR4_CTX200","TR4_CTX200","TR4_CTX200","TR4_CTX200","TR4_CTX200","TR1_CTX+HV292.120","TR1_CTX+HV292.120","TR1_CTX+HV292.120","TR1_CTX+HV292.120","TR1_CTX+HV292.120","TR1_CTX+HV292.120","TR1_CTX+HV292.120","TR3_CTX+HV292.1200","TR3_CTX+HV292.1200","TR3_CTX+HV292.1200","TR3_CTX+HV292.1200","TR3_CTX+HV292.1200","TR3_CTX+HV292.1200","TR5_HV292.1","TR5_HV292.1","TR5_HV292.1","TR5_HV292.1","TR5_HV292.1","TR2_VAN90","TR2_VAN90","TR2_VAN90","TR2_VAN90","TR4_VAN600","TR4_VAN600","TR4_VAN600","TR1_VAN+CCUG5916890","TR1_VAN+CCUG5916890","TR3_VAN+CCUG59168600","CR_UNTREATED","CR_UNTREATED","CR_UNTREATED","CR_UNTREATED","CR_UNTREATED","CR_UNTREATED","CR_UNTREATED","CR_UNTREATED","CR_UNTREATED","CR_UNTREATED","TR2_CTX20","TR2_CTX20","TR2_CTX20","TR2_CTX20","TR2_CTX20","TR2_CTX20","TR2_CTX20","TR2_CTX20","TR2_CTX20","TR4_CTX200","TR4_CTX200","TR4_CTX200","TR4_CTX200","TR4_CTX200","TR4_CTX200","TR4_CTX200","TR4_CTX200","TR1_CTX+HV292.120","TR1_CTX+HV292.120","TR1_CTX+HV292.120","TR1_CTX+HV292.120","TR1_CTX+HV292.120","TR1_CTX+HV292.120","TR1_CTX+HV292.120","TR3_CTX+HV292.1200","TR3_CTX+HV292.1200","TR3_CTX+HV292.1200","TR3_CTX+HV292.1200","TR3_CTX+HV292.1200","TR3_CTX+HV292.1200","TR5_HV292.1","TR5_HV292.1","TR5_HV292.1","TR5_HV292.1","TR5_HV292.1","TR2_VAN90","TR2_VAN90","TR2_VAN90","TR2_VAN90","TR4_VAN600","TR4_VAN600","TR4_VAN600","TR1_VAN+CCUG5916890","TR1_VAN+CCUG5916890","TR3_VAN+CCUG59168600"],["TR2_CTX20","TR4_CTX200","TR1_CTX+HV292.120","TR3_CTX+HV292.1200","TR5_HV292.1","TR2_VAN90","TR4_VAN600","TR1_VAN+CCUG5916890","TR3_VAN+CCUG59168600","TR5_CCUG59168","TR4_CTX200","TR1_CTX+HV292.120","TR3_CTX+HV292.1200","TR5_HV292.1","TR2_VAN90","TR4_VAN600","TR1_VAN+CCUG5916890","TR3_VAN+CCUG59168600","TR5_CCUG59168","TR1_CTX+HV292.120","TR3_CTX+HV292.1200","TR5_HV292.1","TR2_VAN90","TR4_VAN600","TR1_VAN+CCUG5916890","TR3_VAN+CCUG59168600","TR5_CCUG59168","TR3_CTX+HV292.1200","TR5_HV292.1","TR2_VAN90","TR4_VAN600","TR1_VAN+CCUG5916890","TR3_VAN+CCUG59168600","TR5_CCUG59168","TR5_HV292.1","TR2_VAN90","TR4_VAN600","TR1_VAN+CCUG5916890","TR3_VAN+CCUG59168600","TR5_CCUG59168","TR2_VAN90","TR4_VAN600","TR1_VAN+CCUG5916890","TR3_VAN+CCUG59168600","TR5_CCUG59168","TR4_VAN600","TR1_VAN+CCUG5916890","TR3_VAN+CCUG59168600","TR5_CCUG59168","TR1_VAN+CCUG5916890","TR3_VAN+CCUG59168600","TR5_CCUG59168","TR3_VAN+CCUG59168600","TR5_CCUG59168","TR5_CCUG59168","TR2_CTX20","TR4_CTX200","TR1_CTX+HV292.120","TR3_CTX+HV292.1200","TR5_HV292.1","TR2_VAN90","TR4_VAN600","TR1_VAN+CCUG5916890","TR3_VAN+CCUG59168600","TR5_CCUG59168","TR4_CTX200","TR1_CTX+HV292.120","TR3_CTX+HV292.1200","TR5_HV292.1","TR2_VAN90","TR4_VAN600","TR1_VAN+CCUG5916890","TR3_VAN+CCUG59168600","TR5_CCUG59168","TR1_CTX+HV292.120","TR3_CTX+HV292.1200","TR5_HV292.1","TR2_VAN90","TR4_VAN600","TR1_VAN+CCUG5916890","TR3_VAN+CCUG59168600","TR5_CCUG59168","TR3_CTX+HV292.1200","TR5_HV292.1","TR2_VAN90","TR4_VAN600","TR1_VAN+CCUG5916890","TR3_VAN+CCUG59168600","TR5_CCUG59168","TR5_HV292.1","TR2_VAN90","TR4_VAN600","TR1_VAN+CCUG5916890","TR3_VAN+CCUG59168600","TR5_CCUG59168","TR2_VAN90","TR4_VAN600","TR1_VAN+CCUG5916890","TR3_VAN+CCUG59168600","TR5_CCUG59168","TR4_VAN600","TR1_VAN+CCUG5916890","TR3_VAN+CCUG59168600","TR5_CCUG59168","TR1_VAN+CCUG5916890","TR3_VAN+CCUG59168600","TR5_CCUG59168","TR3_VAN+CCUG59168600","TR5_CCUG59168","TR5_CCUG59168"],[0.227,0.248,0.196,0.15,0.234,0.244,0.143,0.204,0.428,0.193,0.45,0.302,0.358,0.501,0.644,0.683,0.652,0.674,0.55,0.397,0.454,0.543,0.57,0.57,0.58,0.654,0.563,0.341,0.438,0.531,0.487,0.556,0.592,0.452,0.563,0.582,0.542,0.521,0.63,0.454,0.641,0.699,0.668,0.658,0.597,0.603,0.424,0.501,0.354,0.368,0.603,0.375,0.589,0.351,0.495,0.238,0.229,0.179,0.217,0.235,0.266,0.173,0.251,0.403,0.228,0.44,0.309,0.279,0.455,0.601,0.598,0.622,0.636,0.528,0.372,0.455,0.531,0.58,0.591,0.599,0.627,0.555,0.324,0.392,0.529,0.47,0.549,0.579,0.462,0.501,0.585,0.56,0.578,0.63,0.494,0.632,0.648,0.654,0.635,0.571,0.619,0.432,0.476,0.38,0.457,0.594,0.411,0.562,0.324,0.478],[0.064,0.087,0.122,0.313,0.046,0.103,0.368,0.071,0.004,0.11,0.1,0.2,0.2,0.1,0.1,0.1,0.027,0.02,0.1,0.1,0.1,0.1,0.333,0.333,0.067,0.067,0.1,0.1,0.1,0.1,0.1,0.032,0.038,0.1,0.1,0.1,0.1,0.045,0.029,0.1,0.1,0.1,0.028,0.033,0.1,0.333,0.067,0.067,0.1,0.067,0.067,0.1,0.024,0.029,0.03,0.047,0.135,0.175,0.048,0.05,0.079,0.244,0.046,0.007,0.045,0.1,0.1,0.2,0.1,0.1,0.1,0.034,0.025,0.1,0.1,0.1,0.1,0.333,0.333,0.067,0.067,0.1,0.1,0.1,0.1,0.1,0.025,0.024,0.1,0.1,0.1,0.1,0.038,0.024,0.1,0.1,0.1,0.027,0.025,0.1,0.333,0.067,0.067,0.1,0.067,0.067,0.1,0.021,0.026,0.028],[0.251,0.217,0.14,0.338,0.195,0.123,0.368,0.186,0.22,0.129,0.162,0.222,0.222,0.162,0.162,0.162,0.371,0.55,0.162,0.162,0.162,0.162,0.346,0.346,0.21,0.21,0.162,0.162,0.162,0.162,0.162,0.196,0.19,0.162,0.162,0.162,0.162,0.206,0.245,0.162,0.162,0.162,0.308,0.182,0.162,0.346,0.21,0.21,0.162,0.21,0.21,0.162,0.44,0.245,0.206,0.172,0.152,0.192,0.165,0.162,0.181,0.258,0.181,0.385,0.19,0.151,0.151,0.216,0.151,0.151,0.151,0.17,0.229,0.151,0.151,0.151,0.151,0.34,0.34,0.179,0.179,0.151,0.151,0.151,0.151,0.151,0.229,0.377,0.151,0.151,0.151,0.151,0.174,0.377,0.151,0.151,0.151,0.165,0.229,0.151,0.34,0.179,0.179,0.151,0.179,0.179,0.151,0.578,0.179,0.154]],"container":"<table class=\"display\">\n  <thead>\n    <tr>\n      <th> <\/th>\n      <th>Distance<\/th>\n      <th>X1<\/th>\n      <th>X2<\/th>\n      <th>R2<\/th>\n      <th>pval<\/th>\n      <th>pvalFDR<\/th>\n    <\/tr>\n  <\/thead>\n<\/table>","options":{"columnDefs":[{"className":"dt-right","targets":[4,5,6]},{"orderable":false,"targets":0}],"order":[],"autoWidth":false,"orderClasses":false}},"evals":[],"jsHooks":[]}</script>
```

```r
tmp_df %>% 
  xlsx::write.xlsx(file = out_xlsx,
                   append = TRUE,
                   sheetName = "human1all-beta-pretreat-pw",
                   showNA = TRUE)
```


```r
lapply(
  dist,
  FUN = phyloseq_TW,
  physeq =  ps_tmp_forbdiv %>%
    subset_samples(Day_of_Treatment <= 0),
  variable = "Reactor_Treatment_Dose",
  nrep = 999,
  strata = NULL
) %>%
  bind_rows(.id = "Distance") %>%
  mutate_if(is.numeric, round, 3)   -> tmp_df

tmp_df %>% 
  DT::datatable()
```

```{=html}
<div class="datatables html-widget html-fill-item-overflow-hidden html-fill-item" id="htmlwidget-4d98ecc6a7bce25e95cf" style="width:100%;height:auto;"></div>
<script type="application/json" data-for="htmlwidget-4d98ecc6a7bce25e95cf">{"x":{"filter":"none","vertical":false,"data":[["1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","23","24","25","26","27","28","29","30","31","32","33","34","35","36","37","38","39","40","41","42","43","44","45","46","47","48","49","50","51","52","53","54","55","56","57","58","59","60","61","62","63","64","65","66","67","68","69","70","71","72","73","74","75","76","77","78","79","80","81","82","83","84","85","86","87","88","89","90","91","92","93","94","95","96","97","98","99","100","101","102","103","104","105","106","107","108","109","110"],["bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson"],[["CR_UNTREATED"],["CR_UNTREATED"],["CR_UNTREATED"],["CR_UNTREATED"],["CR_UNTREATED"],["CR_UNTREATED"],["CR_UNTREATED"],["CR_UNTREATED"],["CR_UNTREATED"],["CR_UNTREATED"],["TR2_CTX20"],["TR2_CTX20"],["TR2_CTX20"],["TR2_CTX20"],["TR2_CTX20"],["TR2_CTX20"],["TR2_CTX20"],["TR2_CTX20"],["TR2_CTX20"],["TR4_CTX200"],["TR4_CTX200"],["TR4_CTX200"],["TR4_CTX200"],["TR4_CTX200"],["TR4_CTX200"],["TR4_CTX200"],["TR4_CTX200"],["TR1_CTX+HV292.120"],["TR1_CTX+HV292.120"],["TR1_CTX+HV292.120"],["TR1_CTX+HV292.120"],["TR1_CTX+HV292.120"],["TR1_CTX+HV292.120"],["TR1_CTX+HV292.120"],["TR3_CTX+HV292.1200"],["TR3_CTX+HV292.1200"],["TR3_CTX+HV292.1200"],["TR3_CTX+HV292.1200"],["TR3_CTX+HV292.1200"],["TR3_CTX+HV292.1200"],["TR5_HV292.1"],["TR5_HV292.1"],["TR5_HV292.1"],["TR5_HV292.1"],["TR5_HV292.1"],["TR2_VAN90"],["TR2_VAN90"],["TR2_VAN90"],["TR2_VAN90"],["TR4_VAN600"],["TR4_VAN600"],["TR4_VAN600"],["TR1_VAN+CCUG5916890"],["TR1_VAN+CCUG5916890"],["TR3_VAN+CCUG59168600"],["CR_UNTREATED"],["CR_UNTREATED"],["CR_UNTREATED"],["CR_UNTREATED"],["CR_UNTREATED"],["CR_UNTREATED"],["CR_UNTREATED"],["CR_UNTREATED"],["CR_UNTREATED"],["CR_UNTREATED"],["TR2_CTX20"],["TR2_CTX20"],["TR2_CTX20"],["TR2_CTX20"],["TR2_CTX20"],["TR2_CTX20"],["TR2_CTX20"],["TR2_CTX20"],["TR2_CTX20"],["TR4_CTX200"],["TR4_CTX200"],["TR4_CTX200"],["TR4_CTX200"],["TR4_CTX200"],["TR4_CTX200"],["TR4_CTX200"],["TR4_CTX200"],["TR1_CTX+HV292.120"],["TR1_CTX+HV292.120"],["TR1_CTX+HV292.120"],["TR1_CTX+HV292.120"],["TR1_CTX+HV292.120"],["TR1_CTX+HV292.120"],["TR1_CTX+HV292.120"],["TR3_CTX+HV292.1200"],["TR3_CTX+HV292.1200"],["TR3_CTX+HV292.1200"],["TR3_CTX+HV292.1200"],["TR3_CTX+HV292.1200"],["TR3_CTX+HV292.1200"],["TR5_HV292.1"],["TR5_HV292.1"],["TR5_HV292.1"],["TR5_HV292.1"],["TR5_HV292.1"],["TR2_VAN90"],["TR2_VAN90"],["TR2_VAN90"],["TR2_VAN90"],["TR4_VAN600"],["TR4_VAN600"],["TR4_VAN600"],["TR1_VAN+CCUG5916890"],["TR1_VAN+CCUG5916890"],["TR3_VAN+CCUG59168600"]],[["TR2_CTX20"],["TR4_CTX200"],["TR1_CTX+HV292.120"],["TR3_CTX+HV292.1200"],["TR5_HV292.1"],["TR2_VAN90"],["TR4_VAN600"],["TR1_VAN+CCUG5916890"],["TR3_VAN+CCUG59168600"],["TR5_CCUG59168"],["TR4_CTX200"],["TR1_CTX+HV292.120"],["TR3_CTX+HV292.1200"],["TR5_HV292.1"],["TR2_VAN90"],["TR4_VAN600"],["TR1_VAN+CCUG5916890"],["TR3_VAN+CCUG59168600"],["TR5_CCUG59168"],["TR1_CTX+HV292.120"],["TR3_CTX+HV292.1200"],["TR5_HV292.1"],["TR2_VAN90"],["TR4_VAN600"],["TR1_VAN+CCUG5916890"],["TR3_VAN+CCUG59168600"],["TR5_CCUG59168"],["TR3_CTX+HV292.1200"],["TR5_HV292.1"],["TR2_VAN90"],["TR4_VAN600"],["TR1_VAN+CCUG5916890"],["TR3_VAN+CCUG59168600"],["TR5_CCUG59168"],["TR5_HV292.1"],["TR2_VAN90"],["TR4_VAN600"],["TR1_VAN+CCUG5916890"],["TR3_VAN+CCUG59168600"],["TR5_CCUG59168"],["TR2_VAN90"],["TR4_VAN600"],["TR1_VAN+CCUG5916890"],["TR3_VAN+CCUG59168600"],["TR5_CCUG59168"],["TR4_VAN600"],["TR1_VAN+CCUG5916890"],["TR3_VAN+CCUG59168600"],["TR5_CCUG59168"],["TR1_VAN+CCUG5916890"],["TR3_VAN+CCUG59168600"],["TR5_CCUG59168"],["TR3_VAN+CCUG59168600"],["TR5_CCUG59168"],["TR5_CCUG59168"],["TR2_CTX20"],["TR4_CTX200"],["TR1_CTX+HV292.120"],["TR3_CTX+HV292.1200"],["TR5_HV292.1"],["TR2_VAN90"],["TR4_VAN600"],["TR1_VAN+CCUG5916890"],["TR3_VAN+CCUG59168600"],["TR5_CCUG59168"],["TR4_CTX200"],["TR1_CTX+HV292.120"],["TR3_CTX+HV292.1200"],["TR5_HV292.1"],["TR2_VAN90"],["TR4_VAN600"],["TR1_VAN+CCUG5916890"],["TR3_VAN+CCUG59168600"],["TR5_CCUG59168"],["TR1_CTX+HV292.120"],["TR3_CTX+HV292.1200"],["TR5_HV292.1"],["TR2_VAN90"],["TR4_VAN600"],["TR1_VAN+CCUG5916890"],["TR3_VAN+CCUG59168600"],["TR5_CCUG59168"],["TR3_CTX+HV292.1200"],["TR5_HV292.1"],["TR2_VAN90"],["TR4_VAN600"],["TR1_VAN+CCUG5916890"],["TR3_VAN+CCUG59168600"],["TR5_CCUG59168"],["TR5_HV292.1"],["TR2_VAN90"],["TR4_VAN600"],["TR1_VAN+CCUG5916890"],["TR3_VAN+CCUG59168600"],["TR5_CCUG59168"],["TR2_VAN90"],["TR4_VAN600"],["TR1_VAN+CCUG5916890"],["TR3_VAN+CCUG59168600"],["TR5_CCUG59168"],["TR4_VAN600"],["TR1_VAN+CCUG5916890"],["TR3_VAN+CCUG59168600"],["TR5_CCUG59168"],["TR1_VAN+CCUG5916890"],["TR3_VAN+CCUG59168600"],["TR5_CCUG59168"],["TR3_VAN+CCUG59168600"],["TR5_CCUG59168"],["TR5_CCUG59168"]],[[6],[6],[6],[6],[6],[6],[6],[6],[6],[6],[3],[3],[3],[3],[3],[3],[3],[3],[3],[2],[2],[2],[2],[2],[2],[2],[2],[3],[3],[3],[3],[3],[3],[3],[3],[3],[3],[3],[3],[3],[3],[3],[3],[3],[3],[2],[2],[2],[2],[2],[2],[2],[4],[4],[4],[6],[6],[6],[6],[6],[6],[6],[6],[6],[6],[3],[3],[3],[3],[3],[3],[3],[3],[3],[2],[2],[2],[2],[2],[2],[2],[2],[3],[3],[3],[3],[3],[3],[3],[3],[3],[3],[3],[3],[3],[3],[3],[3],[3],[3],[2],[2],[2],[2],[2],[2],[2],[4],[4],[4]],[[3],[2],[3],[3],[3],[2],[2],[4],[4],[3],[2],[3],[3],[3],[2],[2],[4],[4],[3],[3],[3],[3],[2],[2],[4],[4],[3],[3],[3],[2],[2],[4],[4],[3],[3],[2],[2],[4],[4],[3],[2],[2],[4],[4],[3],[2],[4],[4],[3],[4],[4],[3],[4],[3],[3],[3],[2],[3],[3],[3],[2],[2],[4],[4],[3],[2],[3],[3],[3],[2],[2],[4],[4],[3],[3],[3],[3],[2],[2],[4],[4],[3],[3],[3],[2],[2],[4],[4],[3],[3],[2],[2],[4],[4],[3],[2],[2],[4],[4],[3],[2],[4],[4],[3],[4],[4],[3],[4],[3],[3]],[[0.026],[0.246],[0.08],[0.108],[0.025],[0.223],[0.251],[0.051],[0.008],[0.067],[0.208],[0.176],[0.213],[0.096],[0.103],[0.103],[0.023],[0.029],[0.105],[0.205],[0.121],[0.088],[0.34],[0.328],[0.063],[0.056],[0.108],[0.096],[0.1],[0.101],[0.11],[0.022],[0.033],[0.113],[0.094],[0.089],[0.101],[0.023],[0.039],[0.118],[0.103],[0.091],[0.024],[0.03],[0.113],[0.349],[0.109],[0.072],[0.112],[0.131],[0.076],[0.099],[0.045],[0.048],[0.031],[0.033],[0.188],[0.135],[0.032],[0.03],[0.061],[0.172],[0.018],[0.01],[0.035],[0.186],[0.107],[0.198],[0.11],[0.093],[0.117],[0.03],[0.03],[0.089],[0.103],[0.112],[0.108],[0.344],[0.322],[0.069],[0.06],[0.112],[0.096],[0.086],[0.086],[0.104],[0.037],[0.029],[0.12],[0.09],[0.103],[0.1],[0.024],[0.027],[0.098],[0.091],[0.086],[0.021],[0.029],[0.103],[0.349],[0.137],[0.063],[0.22],[0.06],[0.059],[0.098],[0.03],[0.025],[0.023]],[[3.38062565695451],[1.80592530310392],[1.93917282330703],[2.05049602270045],[3.59631244865084],[2.36697552346636],[1.98116474013619],[2.59066515317003],[6.96642471897285],[2.17022334607896],[1.65381619349415],[1.7324220040701],[2.23340677556841],[4.01443896882939],[4.00120838515932],[6.0919109315107],[10.0654193440074],[11.9054316869394],[4.88671497710981],[1.72879834780556],[1.66624018029317],[2.35895165365943],[2.65210966205633],[2.65248541190764],[3.29855562422723],[5.32605646568159],[3.08124050499551],[2.06589776511485],[3.11260641673556],[3.40574332856624],[3.62010406040427],[5.4705682909732],[6.82679846884913],[3.30522557599556],[5.16060425041003],[3.055831963041],[3.30906917771683],[5.86851320718698],[9.83132025152235],[3.31971251583147],[3.87188520814268],[6.37563972193481],[11.0018248676792],[11.2113490204868],[5.92690778071803],[3.0375180067186],[2.1125144619168],[3.47022764688383],[1.49424093436232],[2.51746542213744],[7.85680875157764],[2.11173093505347],[8.60568172992972],[2.50623226705387],[4.91130049713801],[3.08164791403813],[1.83116814334986],[1.74722882427578],[2.90727077722761],[3.22335217243648],[2.70370707096536],[2.08957428620594],[3.3355172947831],[6.24529892283605],[2.57674073367949],[1.85098110545744],[1.78725825121598],[1.54769976454627],[3.34468876378856],[3.85712454018571],[4.42735278907366],[8.11480376696097],[9.21995087493538],[4.47374029425332],[1.63137411663642],[1.87262022077958],[2.5447529353778],[2.76306139911096],[2.88625486057502],[3.96613555399152],[5.19070684837627],[3.2331945893531],[1.91653516574998],[2.57510363230682],[3.38392974823572],[3.08162440359599],[5.40906796634094],[6.54076924213814],[3.43168377761943],[4.01967167353687],[3.41898404969475],[3.57251002831294],[7.01260867962867],[9.2902739025859],[3.89786034086938],[4.15158980482648],[5.15396995366647],[9.66653535848606],[9.50331959266683],[5.31839422253925],[3.25117874895107],[2.29083787455094],[3.22351485813272],[1.72978106453622],[3.21253348363019],[6.55814440142508],[2.29352429478639],[7.6852018254337],[2.21808851322542],[4.53563511530767]],[[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999]],[0.148,0.265,0.148,0.148,0.148,0.245,0.265,0.148,0.148,0.148,0.238,0.21,0.239,0.148,0.148,0.148,0.148,0.148,0.148,0.238,0.151,0.148,0.346,0.34,0.148,0.148,0.148,0.148,0.148,0.148,0.148,0.148,0.148,0.148,0.148,0.148,0.148,0.148,0.148,0.151,0.148,0.148,0.148,0.148,0.148,0.349,0.148,0.148,0.148,0.16,0.148,0.148,0.148,0.148,0.148,0.12,0.207,0.16,0.12,0.12,0.143,0.197,0.12,0.12,0.12,0.207,0.143,0.214,0.143,0.143,0.146,0.12,0.12,0.143,0.143,0.143,0.143,0.349,0.334,0.143,0.143,0.143,0.143,0.143,0.143,0.143,0.12,0.12,0.147,0.143,0.143,0.143,0.12,0.12,0.143,0.143,0.143,0.12,0.12,0.143,0.349,0.16,0.143,0.233,0.143,0.143,0.143,0.12,0.12,0.12]],"container":"<table class=\"display\">\n  <thead>\n    <tr>\n      <th> <\/th>\n      <th>Distance<\/th>\n      <th>Level1<\/th>\n      <th>Level2<\/th>\n      <th>N1<\/th>\n      <th>N2<\/th>\n      <th>p.value<\/th>\n      <th>tw2.stat<\/th>\n      <th>nrep<\/th>\n      <th>pvalFDR<\/th>\n    <\/tr>\n  <\/thead>\n<\/table>","options":{"columnDefs":[{"className":"dt-right","targets":9},{"orderable":false,"targets":0}],"order":[],"autoWidth":false,"orderClasses":false}},"evals":[],"jsHooks":[]}</script>
```

```r
tmp_df %>% 
  xlsx::write.xlsx(file = out_xlsx,
                   append = TRUE,
                   sheetName = "human1all-beta-pretreat-pw2",
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

##### treat - Day_of_Treatment > 1:


```r
lapply(
  dist,
  FUN = phyloseq_adonis2,
  physeq = ps_tmp_forbdiv %>%
    subset_samples(Day_of_Treatment >1),
  formula = paste0(c("Reactor_Treatment_Dose"), collapse=" + "),
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
<div class="datatables html-widget html-fill-item-overflow-hidden html-fill-item" id="htmlwidget-7a55dcb7a06bbff1b49e" style="width:100%;height:auto;"></div>
<script type="application/json" data-for="htmlwidget-7a55dcb7a06bbff1b49e">{"x":{"filter":"none","vertical":false,"data":[["1","2","3","4","5","6"],["bjaccard","bjaccard","bjaccard","aichinson","aichinson","aichinson"],["Reactor_Treatment_Dose","Residual","Total","Reactor_Treatment_Dose","Residual","Total"],[10,41,51,10,41,51],[7.008,3.686,10.694,4.924,2.491,7.414],[0.655,0.345,1,0.664,0.336,1],[7.795,null,null,8.105,null,null],[0.001,null,null,0.001,null,null]],"container":"<table class=\"display\">\n  <thead>\n    <tr>\n      <th> <\/th>\n      <th>Distance<\/th>\n      <th>terms<\/th>\n      <th>Df<\/th>\n      <th>SumOfSqs<\/th>\n      <th>R2<\/th>\n      <th>F<\/th>\n      <th>Pr(&gt;F)<\/th>\n    <\/tr>\n  <\/thead>\n<\/table>","options":{"columnDefs":[{"className":"dt-right","targets":[3,4,5,6,7]},{"orderable":false,"targets":0}],"order":[],"autoWidth":false,"orderClasses":false}},"evals":[],"jsHooks":[]}</script>
```

```r
tmp_df %>% 
  xlsx::write.xlsx(file = out_xlsx,
                   append = TRUE,
                   sheetName = "human1all-beta-treat",
                   showNA = TRUE)
```


```r
lapply(
  dist,
  FUN = phyloseq_adonis2,
  physeq = ps_tmp_forbdiv %>%
    subset_samples(Day_of_Treatment >1),
  formula = paste0(c("Reactor_Treatment_Dose", "Day_of_Treatment"), collapse=" + "),
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
<div class="datatables html-widget html-fill-item-overflow-hidden html-fill-item" id="htmlwidget-2dc7173b9d96e4df10d9" style="width:100%;height:auto;"></div>
<script type="application/json" data-for="htmlwidget-2dc7173b9d96e4df10d9">{"x":{"filter":"none","vertical":false,"data":[["1","2","3","4","5","6","7","8"],["bjaccard","bjaccard","bjaccard","bjaccard","aichinson","aichinson","aichinson","aichinson"],["Reactor_Treatment_Dose","Day_of_Treatment","Residual","Total","Reactor_Treatment_Dose","Day_of_Treatment","Residual","Total"],[10,1,40,51,10,1,40,51],[6.934,0.241,3.445,10.694,4.88,0.144,2.346,7.414],[0.648,0.023,0.322,1,0.658,0.019,0.316,1],[8.051,2.798,null,null,8.319,2.457,null,null],[0.001,0.026,null,null,0.001,0.043,null,null]],"container":"<table class=\"display\">\n  <thead>\n    <tr>\n      <th> <\/th>\n      <th>Distance<\/th>\n      <th>terms<\/th>\n      <th>Df<\/th>\n      <th>SumOfSqs<\/th>\n      <th>R2<\/th>\n      <th>F<\/th>\n      <th>Pr(&gt;F)<\/th>\n    <\/tr>\n  <\/thead>\n<\/table>","options":{"columnDefs":[{"className":"dt-right","targets":[3,4,5,6,7]},{"orderable":false,"targets":0}],"order":[],"autoWidth":false,"orderClasses":false}},"evals":[],"jsHooks":[]}</script>
```

```r
tmp_df %>% 
  xlsx::write.xlsx(file = out_xlsx,
                   append = TRUE,
                   sheetName = "human1all-beta-treat-time-day",
                   showNA = TRUE)
```


```r
lapply(
  dist,
  FUN = physeq_pairwise_permanovas_adonis2,
  physeq = ps_tmp_forbdiv %>%
    subset_samples(Day_of_Treatment > 1),
  compare_header = "Reactor_Treatment_Dose",
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

```r
tmp_df %>% 
  DT::datatable()
```

```{=html}
<div class="datatables html-widget html-fill-item-overflow-hidden html-fill-item" id="htmlwidget-baf8b8aa16d9e32c874e" style="width:100%;height:auto;"></div>
<script type="application/json" data-for="htmlwidget-baf8b8aa16d9e32c874e">{"x":{"filter":"none","vertical":false,"data":[["1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","23","24","25","26","27","28","29","30","31","32","33","34","35","36","37","38","39","40","41","42","43","44","45","46","47","48","49","50","51","52","53","54","55","56","57","58","59","60","61","62","63","64","65","66","67","68","69","70","71","72","73","74","75","76","77","78","79","80","81","82","83","84","85","86","87","88","89","90","91","92","93","94","95","96","97","98","99","100","101","102","103","104","105","106","107","108","109","110"],["bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson"],["CR_UNTREATED","CR_UNTREATED","CR_UNTREATED","CR_UNTREATED","CR_UNTREATED","CR_UNTREATED","CR_UNTREATED","CR_UNTREATED","CR_UNTREATED","CR_UNTREATED","TR2_CTX20","TR2_CTX20","TR2_CTX20","TR2_CTX20","TR2_CTX20","TR2_CTX20","TR2_CTX20","TR2_CTX20","TR2_CTX20","TR4_CTX200","TR4_CTX200","TR4_CTX200","TR4_CTX200","TR4_CTX200","TR4_CTX200","TR4_CTX200","TR4_CTX200","TR1_CTX+HV292.120","TR1_CTX+HV292.120","TR1_CTX+HV292.120","TR1_CTX+HV292.120","TR1_CTX+HV292.120","TR1_CTX+HV292.120","TR1_CTX+HV292.120","TR3_CTX+HV292.1200","TR3_CTX+HV292.1200","TR3_CTX+HV292.1200","TR3_CTX+HV292.1200","TR3_CTX+HV292.1200","TR3_CTX+HV292.1200","TR5_HV292.1","TR5_HV292.1","TR5_HV292.1","TR5_HV292.1","TR5_HV292.1","TR2_VAN90","TR2_VAN90","TR2_VAN90","TR2_VAN90","TR4_VAN600","TR4_VAN600","TR4_VAN600","TR1_VAN+CCUG5916890","TR1_VAN+CCUG5916890","TR3_VAN+CCUG59168600","CR_UNTREATED","CR_UNTREATED","CR_UNTREATED","CR_UNTREATED","CR_UNTREATED","CR_UNTREATED","CR_UNTREATED","CR_UNTREATED","CR_UNTREATED","CR_UNTREATED","TR2_CTX20","TR2_CTX20","TR2_CTX20","TR2_CTX20","TR2_CTX20","TR2_CTX20","TR2_CTX20","TR2_CTX20","TR2_CTX20","TR4_CTX200","TR4_CTX200","TR4_CTX200","TR4_CTX200","TR4_CTX200","TR4_CTX200","TR4_CTX200","TR4_CTX200","TR1_CTX+HV292.120","TR1_CTX+HV292.120","TR1_CTX+HV292.120","TR1_CTX+HV292.120","TR1_CTX+HV292.120","TR1_CTX+HV292.120","TR1_CTX+HV292.120","TR3_CTX+HV292.1200","TR3_CTX+HV292.1200","TR3_CTX+HV292.1200","TR3_CTX+HV292.1200","TR3_CTX+HV292.1200","TR3_CTX+HV292.1200","TR5_HV292.1","TR5_HV292.1","TR5_HV292.1","TR5_HV292.1","TR5_HV292.1","TR2_VAN90","TR2_VAN90","TR2_VAN90","TR2_VAN90","TR4_VAN600","TR4_VAN600","TR4_VAN600","TR1_VAN+CCUG5916890","TR1_VAN+CCUG5916890","TR3_VAN+CCUG59168600"],["TR2_CTX20","TR4_CTX200","TR1_CTX+HV292.120","TR3_CTX+HV292.1200","TR5_HV292.1","TR2_VAN90","TR4_VAN600","TR1_VAN+CCUG5916890","TR3_VAN+CCUG59168600","TR5_CCUG59168","TR4_CTX200","TR1_CTX+HV292.120","TR3_CTX+HV292.1200","TR5_HV292.1","TR2_VAN90","TR4_VAN600","TR1_VAN+CCUG5916890","TR3_VAN+CCUG59168600","TR5_CCUG59168","TR1_CTX+HV292.120","TR3_CTX+HV292.1200","TR5_HV292.1","TR2_VAN90","TR4_VAN600","TR1_VAN+CCUG5916890","TR3_VAN+CCUG59168600","TR5_CCUG59168","TR3_CTX+HV292.1200","TR5_HV292.1","TR2_VAN90","TR4_VAN600","TR1_VAN+CCUG5916890","TR3_VAN+CCUG59168600","TR5_CCUG59168","TR5_HV292.1","TR2_VAN90","TR4_VAN600","TR1_VAN+CCUG5916890","TR3_VAN+CCUG59168600","TR5_CCUG59168","TR2_VAN90","TR4_VAN600","TR1_VAN+CCUG5916890","TR3_VAN+CCUG59168600","TR5_CCUG59168","TR4_VAN600","TR1_VAN+CCUG5916890","TR3_VAN+CCUG59168600","TR5_CCUG59168","TR1_VAN+CCUG5916890","TR3_VAN+CCUG59168600","TR5_CCUG59168","TR3_VAN+CCUG59168600","TR5_CCUG59168","TR5_CCUG59168","TR2_CTX20","TR4_CTX200","TR1_CTX+HV292.120","TR3_CTX+HV292.1200","TR5_HV292.1","TR2_VAN90","TR4_VAN600","TR1_VAN+CCUG5916890","TR3_VAN+CCUG59168600","TR5_CCUG59168","TR4_CTX200","TR1_CTX+HV292.120","TR3_CTX+HV292.1200","TR5_HV292.1","TR2_VAN90","TR4_VAN600","TR1_VAN+CCUG5916890","TR3_VAN+CCUG59168600","TR5_CCUG59168","TR1_CTX+HV292.120","TR3_CTX+HV292.1200","TR5_HV292.1","TR2_VAN90","TR4_VAN600","TR1_VAN+CCUG5916890","TR3_VAN+CCUG59168600","TR5_CCUG59168","TR3_CTX+HV292.1200","TR5_HV292.1","TR2_VAN90","TR4_VAN600","TR1_VAN+CCUG5916890","TR3_VAN+CCUG59168600","TR5_CCUG59168","TR5_HV292.1","TR2_VAN90","TR4_VAN600","TR1_VAN+CCUG5916890","TR3_VAN+CCUG59168600","TR5_CCUG59168","TR2_VAN90","TR4_VAN600","TR1_VAN+CCUG5916890","TR3_VAN+CCUG59168600","TR5_CCUG59168","TR4_VAN600","TR1_VAN+CCUG5916890","TR3_VAN+CCUG59168600","TR5_CCUG59168","TR1_VAN+CCUG5916890","TR3_VAN+CCUG59168600","TR5_CCUG59168","TR3_VAN+CCUG59168600","TR5_CCUG59168","TR5_CCUG59168"],[0.311,0.355,0.367,0.441,0.145,0.592,0.49,0.764,0.694,0.163,0.055,0.107,0.157,0.202,0.47,0.394,0.654,0.587,0.368,0.197,0.133,0.304,0.477,0.369,0.689,0.608,0.471,0.339,0.673,0.554,0.428,0.826,0.72,0.637,0.385,0.561,0.459,0.773,0.695,0.598,0.374,0.246,0.744,0.579,0.533,0.205,0.264,0.398,0.63,0.305,0.295,0.491,0.545,0.854,0.764,0.294,0.301,0.329,0.385,0.149,0.559,0.529,0.649,0.618,0.188,0.063,0.157,0.176,0.251,0.49,0.46,0.609,0.581,0.415,0.226,0.136,0.325,0.485,0.433,0.626,0.593,0.475,0.388,0.669,0.651,0.574,0.824,0.792,0.609,0.407,0.578,0.536,0.716,0.688,0.59,0.515,0.409,0.805,0.75,0.505,0.215,0.265,0.406,0.687,0.317,0.311,0.616,0.584,0.821,0.788],[0.003,0.004,0.006,0.001,0.215,0.002,0.001,0.001,0.003,0.026,0.971,0.531,0.062,0.43,0.001,0.003,0.006,0.007,0.002,0.409,0.351,0.4,0.02,0.003,0.022,0.025,0.008,0.015,0.25,0.02,0.007,0.028,0.021,0.014,0.167,0.007,0.002,0.005,0.009,0.012,0.167,0.275,0.2,0.2,0.167,0.023,0.054,0.017,0.009,0.006,0.004,0.003,0.029,0.007,0.008,0.001,0.001,0.004,0.001,0.196,0.002,0.001,0.003,0.002,0.007,0.938,0.167,0.057,0.323,0.005,0.001,0.004,0.006,0.003,0.232,0.326,0.2,0.009,0.007,0.033,0.028,0.009,0.015,0.25,0.012,0.009,0.021,0.029,0.019,0.167,0.012,0.004,0.009,0.01,0.008,0.167,0.137,0.2,0.2,0.167,0.01,0.04,0.02,0.013,0.004,0.005,0.002,0.029,0.007,0.006],[0.016,0.016,0.019,0.022,0.252,0.018,0.022,0.022,0.016,0.039,0.971,0.541,0.083,0.446,0.022,0.016,0.019,0.019,0.018,0.433,0.386,0.431,0.035,0.016,0.036,0.038,0.019,0.028,0.286,0.035,0.019,0.041,0.035,0.028,0.213,0.019,0.018,0.018,0.019,0.024,0.213,0.309,0.242,0.242,0.213,0.036,0.074,0.031,0.019,0.019,0.016,0.016,0.041,0.019,0.019,0.018,0.018,0.018,0.018,0.229,0.016,0.018,0.017,0.016,0.019,0.938,0.2,0.076,0.335,0.018,0.018,0.018,0.019,0.017,0.25,0.332,0.224,0.02,0.019,0.047,0.043,0.02,0.026,0.264,0.022,0.02,0.033,0.043,0.032,0.208,0.022,0.018,0.02,0.02,0.02,0.208,0.179,0.224,0.224,0.208,0.02,0.055,0.032,0.023,0.018,0.018,0.016,0.043,0.019,0.019]],"container":"<table class=\"display\">\n  <thead>\n    <tr>\n      <th> <\/th>\n      <th>Distance<\/th>\n      <th>X1<\/th>\n      <th>X2<\/th>\n      <th>R2<\/th>\n      <th>pval<\/th>\n      <th>pvalFDR<\/th>\n    <\/tr>\n  <\/thead>\n<\/table>","options":{"columnDefs":[{"className":"dt-right","targets":[4,5,6]},{"orderable":false,"targets":0}],"order":[],"autoWidth":false,"orderClasses":false}},"evals":[],"jsHooks":[]}</script>
```

```r
tmp_df %>% 
  xlsx::write.xlsx(file = out_xlsx,
                   append = TRUE,
                   sheetName = "human1all-beta-treat-pw",
                   showNA = TRUE)
```


```r
lapply(
  dist,
  FUN = phyloseq_TW,
  physeq =  ps_tmp_forbdiv %>%
    subset_samples(Day_of_Treatment >1),
  variable = "Reactor_Treatment_Dose",
  nrep = 999,
  strata = NULL
) %>%
  bind_rows(.id = "Distance") %>%
  mutate_if(is.numeric, round, 3)   -> tmp_df

tmp_df %>% 
  DT::datatable()
```

```{=html}
<div class="datatables html-widget html-fill-item-overflow-hidden html-fill-item" id="htmlwidget-fefaa3ec4226a0a5ecd9" style="width:100%;height:auto;"></div>
<script type="application/json" data-for="htmlwidget-fefaa3ec4226a0a5ecd9">{"x":{"filter":"none","vertical":false,"data":[["1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","23","24","25","26","27","28","29","30","31","32","33","34","35","36","37","38","39","40","41","42","43","44","45","46","47","48","49","50","51","52","53","54","55","56","57","58","59","60","61","62","63","64","65","66","67","68","69","70","71","72","73","74","75","76","77","78","79","80","81","82","83","84","85","86","87","88","89","90","91","92","93","94","95","96","97","98","99","100","101","102","103","104","105","106","107","108","109","110"],["bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson"],[["CR_UNTREATED"],["CR_UNTREATED"],["CR_UNTREATED"],["CR_UNTREATED"],["CR_UNTREATED"],["CR_UNTREATED"],["CR_UNTREATED"],["CR_UNTREATED"],["CR_UNTREATED"],["CR_UNTREATED"],["TR2_CTX20"],["TR2_CTX20"],["TR2_CTX20"],["TR2_CTX20"],["TR2_CTX20"],["TR2_CTX20"],["TR2_CTX20"],["TR2_CTX20"],["TR2_CTX20"],["TR4_CTX200"],["TR4_CTX200"],["TR4_CTX200"],["TR4_CTX200"],["TR4_CTX200"],["TR4_CTX200"],["TR4_CTX200"],["TR4_CTX200"],["TR1_CTX+HV292.120"],["TR1_CTX+HV292.120"],["TR1_CTX+HV292.120"],["TR1_CTX+HV292.120"],["TR1_CTX+HV292.120"],["TR1_CTX+HV292.120"],["TR1_CTX+HV292.120"],["TR3_CTX+HV292.1200"],["TR3_CTX+HV292.1200"],["TR3_CTX+HV292.1200"],["TR3_CTX+HV292.1200"],["TR3_CTX+HV292.1200"],["TR3_CTX+HV292.1200"],["TR5_HV292.1"],["TR5_HV292.1"],["TR5_HV292.1"],["TR5_HV292.1"],["TR5_HV292.1"],["TR2_VAN90"],["TR2_VAN90"],["TR2_VAN90"],["TR2_VAN90"],["TR4_VAN600"],["TR4_VAN600"],["TR4_VAN600"],["TR1_VAN+CCUG5916890"],["TR1_VAN+CCUG5916890"],["TR3_VAN+CCUG59168600"],["CR_UNTREATED"],["CR_UNTREATED"],["CR_UNTREATED"],["CR_UNTREATED"],["CR_UNTREATED"],["CR_UNTREATED"],["CR_UNTREATED"],["CR_UNTREATED"],["CR_UNTREATED"],["CR_UNTREATED"],["TR2_CTX20"],["TR2_CTX20"],["TR2_CTX20"],["TR2_CTX20"],["TR2_CTX20"],["TR2_CTX20"],["TR2_CTX20"],["TR2_CTX20"],["TR2_CTX20"],["TR4_CTX200"],["TR4_CTX200"],["TR4_CTX200"],["TR4_CTX200"],["TR4_CTX200"],["TR4_CTX200"],["TR4_CTX200"],["TR4_CTX200"],["TR1_CTX+HV292.120"],["TR1_CTX+HV292.120"],["TR1_CTX+HV292.120"],["TR1_CTX+HV292.120"],["TR1_CTX+HV292.120"],["TR1_CTX+HV292.120"],["TR1_CTX+HV292.120"],["TR3_CTX+HV292.1200"],["TR3_CTX+HV292.1200"],["TR3_CTX+HV292.1200"],["TR3_CTX+HV292.1200"],["TR3_CTX+HV292.1200"],["TR3_CTX+HV292.1200"],["TR5_HV292.1"],["TR5_HV292.1"],["TR5_HV292.1"],["TR5_HV292.1"],["TR5_HV292.1"],["TR2_VAN90"],["TR2_VAN90"],["TR2_VAN90"],["TR2_VAN90"],["TR4_VAN600"],["TR4_VAN600"],["TR4_VAN600"],["TR1_VAN+CCUG5916890"],["TR1_VAN+CCUG5916890"],["TR3_VAN+CCUG59168600"]],[["TR2_CTX20"],["TR4_CTX200"],["TR1_CTX+HV292.120"],["TR3_CTX+HV292.1200"],["TR5_HV292.1"],["TR2_VAN90"],["TR4_VAN600"],["TR1_VAN+CCUG5916890"],["TR3_VAN+CCUG59168600"],["TR5_CCUG59168"],["TR4_CTX200"],["TR1_CTX+HV292.120"],["TR3_CTX+HV292.1200"],["TR5_HV292.1"],["TR2_VAN90"],["TR4_VAN600"],["TR1_VAN+CCUG5916890"],["TR3_VAN+CCUG59168600"],["TR5_CCUG59168"],["TR1_CTX+HV292.120"],["TR3_CTX+HV292.1200"],["TR5_HV292.1"],["TR2_VAN90"],["TR4_VAN600"],["TR1_VAN+CCUG5916890"],["TR3_VAN+CCUG59168600"],["TR5_CCUG59168"],["TR3_CTX+HV292.1200"],["TR5_HV292.1"],["TR2_VAN90"],["TR4_VAN600"],["TR1_VAN+CCUG5916890"],["TR3_VAN+CCUG59168600"],["TR5_CCUG59168"],["TR5_HV292.1"],["TR2_VAN90"],["TR4_VAN600"],["TR1_VAN+CCUG5916890"],["TR3_VAN+CCUG59168600"],["TR5_CCUG59168"],["TR2_VAN90"],["TR4_VAN600"],["TR1_VAN+CCUG5916890"],["TR3_VAN+CCUG59168600"],["TR5_CCUG59168"],["TR4_VAN600"],["TR1_VAN+CCUG5916890"],["TR3_VAN+CCUG59168600"],["TR5_CCUG59168"],["TR1_VAN+CCUG5916890"],["TR3_VAN+CCUG59168600"],["TR5_CCUG59168"],["TR3_VAN+CCUG59168600"],["TR5_CCUG59168"],["TR5_CCUG59168"],["TR2_CTX20"],["TR4_CTX200"],["TR1_CTX+HV292.120"],["TR3_CTX+HV292.1200"],["TR5_HV292.1"],["TR2_VAN90"],["TR4_VAN600"],["TR1_VAN+CCUG5916890"],["TR3_VAN+CCUG59168600"],["TR5_CCUG59168"],["TR4_CTX200"],["TR1_CTX+HV292.120"],["TR3_CTX+HV292.1200"],["TR5_HV292.1"],["TR2_VAN90"],["TR4_VAN600"],["TR1_VAN+CCUG5916890"],["TR3_VAN+CCUG59168600"],["TR5_CCUG59168"],["TR1_CTX+HV292.120"],["TR3_CTX+HV292.1200"],["TR5_HV292.1"],["TR2_VAN90"],["TR4_VAN600"],["TR1_VAN+CCUG5916890"],["TR3_VAN+CCUG59168600"],["TR5_CCUG59168"],["TR3_CTX+HV292.1200"],["TR5_HV292.1"],["TR2_VAN90"],["TR4_VAN600"],["TR1_VAN+CCUG5916890"],["TR3_VAN+CCUG59168600"],["TR5_CCUG59168"],["TR5_HV292.1"],["TR2_VAN90"],["TR4_VAN600"],["TR1_VAN+CCUG5916890"],["TR3_VAN+CCUG59168600"],["TR5_CCUG59168"],["TR2_VAN90"],["TR4_VAN600"],["TR1_VAN+CCUG5916890"],["TR3_VAN+CCUG59168600"],["TR5_CCUG59168"],["TR4_VAN600"],["TR1_VAN+CCUG5916890"],["TR3_VAN+CCUG59168600"],["TR5_CCUG59168"],["TR1_VAN+CCUG5916890"],["TR3_VAN+CCUG59168600"],["TR5_CCUG59168"],["TR3_VAN+CCUG59168600"],["TR5_CCUG59168"],["TR5_CCUG59168"]],[[9],[9],[9],[9],[9],[9],[9],[9],[9],[9],[6],[6],[6],[6],[6],[6],[6],[6],[6],[4],[4],[4],[4],[4],[4],[4],[4],[3],[3],[3],[3],[3],[3],[3],[5],[5],[5],[5],[5],[5],[1],[1],[1],[1],[1],[5],[5],[5],[5],[6],[6],[6],[4],[4],[4],[9],[9],[9],[9],[9],[9],[9],[9],[9],[9],[6],[6],[6],[6],[6],[6],[6],[6],[6],[4],[4],[4],[4],[4],[4],[4],[4],[3],[3],[3],[3],[3],[3],[3],[5],[5],[5],[5],[5],[5],[1],[1],[1],[1],[1],[5],[5],[5],[5],[6],[6],[6],[4],[4],[4]],[[6],[4],[3],[5],[1],[5],[6],[4],[4],[5],[4],[3],[5],[1],[5],[6],[4],[4],[5],[3],[5],[1],[5],[6],[4],[4],[5],[5],[1],[5],[6],[4],[4],[5],[1],[5],[6],[4],[4],[5],[5],[6],[4],[4],[5],[6],[4],[4],[5],[4],[4],[5],[4],[5],[5],[6],[4],[3],[5],[1],[5],[6],[4],[4],[5],[4],[3],[5],[1],[5],[6],[4],[4],[5],[3],[5],[1],[5],[6],[4],[4],[5],[5],[1],[5],[6],[4],[4],[5],[1],[5],[6],[4],[4],[5],[5],[6],[4],[4],[5],[6],[4],[4],[5],[4],[4],[5],[4],[5],[5]],[[0.001],[0.007],[0.006],[0.001],[null],[0.001],[0.001],[0.003],[0.004],[0.011],[0.971],[0.268],[0.064],[null],[0.005],[0.001],[0.005],[0.006],[0.005],[0.304],[0.447],[null],[0.008],[0.006],[0.021],[0.028],[0.013],[0.021],[null],[0.026],[0.016],[0.039],[0.027],[0.02],[null],[0.01],[0.003],[0.009],[0.002],[0.006],[null],[null],[null],[null],[null],[0.026],[0.043],[0.011],[0.012],[0.009],[0.005],[0.003],[0.028],[0.012],[0.003],[0.001],[0.009],[0.007],[0.002],[null],[0.002],[0.001],[0.002],[0.003],[0.004],[0.928],[0.085],[0.054],[null],[0.006],[0.002],[0.003],[0.005],[0.005],[0.183],[0.386],[null],[0.009],[0.006],[0.022],[0.027],[0.007],[0.012],[null],[0.017],[0.009],[0.027],[0.025],[0.014],[null],[0.007],[0.004],[0.01],[0.011],[0.015],[null],[null],[null],[null],[null],[0.007],[0.032],[0.007],[0.011],[0.006],[0.006],[0.005],[0.038],[0.011],[0.008]],[[4.91647170497149],[4.27085070955359],[7.06608020054884],[8.95597264896496],[null],[12.3935231443386],[9.43563710681068],[32.9019712640822],[17.6581260033179],[2.95644133831612],[0.471739484141103],[1.26942179750072],[1.78758279459166],[null],[7.74400382271114],[6.51428150856438],[17.5460203534809],[11.4512924903512],[6.00086334959947],[1.47922073693661],[0.989441309836014],[null],[6.65808732227367],[5.38839979867023],[13.2958228782369],[9.30089321458984],[5.23190029743277],[3.59051499230875],[null],[10.8478584037355],[9.11603663887901],[25.9573890518368],[15.4234989143641],[9.27468162301414],[null],[10.2095535644765],[8.60322342066609],[23.8346564315654],[14.6467804253473],[11.9189207973108],[null],[null],[null],[null],[null],[2.39450736648007],[2.82887216695764],[4.83005319458075],[13.6356638778626],[4.54328234191308],[3.85245342133202],[10.231592703],[7.17690713170759],[36.261911234296],[19.0113816463892],[5.10704679051175],[4.03135081504046],[6.65621265012241],[7.98034733784811],[null],[15.2281322375056],[14.0077364128779],[29.2038917058735],[23.2866489963726],[3.35471165582704],[0.523124366218586],[1.76431550865016],[2.01622573575824],[null],[8.86833790611817],[8.51646162516386],[15.7842696458226],[13.4633211168556],[6.94820389277775],[1.69637009804164],[1.02433520209128],[null],[6.28020227176237],[5.84482279086982],[10.0228541747943],[8.75100309415401],[5.6117459858186],[4.24253957811912],[null],[13.1141822991958],[12.4199881441295],[22.1288762872238],[18.702934950957],[9.21525427128353],[null],[10.9372925994209],[10.7868395135244],[19.219806143257],[16.3478744859881],[11.4967350274053],[null],[null],[null],[null],[null],[2.51704331601202],[2.81402379273266],[5.18612450205059],[17.5731536370047],[4.66277166648887],[4.32451256781246],[15.640993365667],[8.40699307799367],[33.2760906231791],[26.0592266427036]],[[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999]],[0.009,0.016,0.014,0.009,null,0.009,0.009,0.014,0.014,0.019,0.971,0.287,0.07,null,0.014,0.009,0.014,0.014,0.014,0.318,0.457,null,0.017,0.014,0.029,0.033,0.02,0.029,null,0.033,0.024,0.045,0.033,0.029,null,0.019,0.014,0.018,0.014,0.014,null,null,null,null,null,0.033,0.048,0.019,0.019,0.018,0.014,0.014,0.033,0.019,0.014,0.014,0.016,0.014,0.014,null,0.014,0.014,0.014,0.014,0.014,0.928,0.091,0.059,null,0.014,0.014,0.014,0.014,0.014,0.192,0.395,null,0.016,0.014,0.028,0.032,0.014,0.017,null,0.023,0.016,0.032,0.031,0.02,null,0.014,0.014,0.016,0.016,0.02,null,null,null,null,null,0.014,0.037,0.014,0.016,0.014,0.014,0.014,0.043,0.016,0.016]],"container":"<table class=\"display\">\n  <thead>\n    <tr>\n      <th> <\/th>\n      <th>Distance<\/th>\n      <th>Level1<\/th>\n      <th>Level2<\/th>\n      <th>N1<\/th>\n      <th>N2<\/th>\n      <th>p.value<\/th>\n      <th>tw2.stat<\/th>\n      <th>nrep<\/th>\n      <th>pvalFDR<\/th>\n    <\/tr>\n  <\/thead>\n<\/table>","options":{"columnDefs":[{"className":"dt-right","targets":9},{"orderable":false,"targets":0}],"order":[],"autoWidth":false,"orderClasses":false}},"evals":[],"jsHooks":[]}</script>
```

```r
tmp_df %>% 
  xlsx::write.xlsx(file = out_xlsx,
                   append = TRUE,
                   sheetName = "human1all-beta-treat-pw2",
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

##### Day_of_Treatment >= 3 & Day_of_Treatment <= 6:


```r
lapply(
  dist,
  FUN = phyloseq_adonis2,
  physeq = ps_tmp_forbdiv %>%
    subset_samples(Day_of_Treatment >= 3 & Day_of_Treatment <= 6),
  formula = paste0(c("Reactor_Treatment_Dose"), collapse=" + "),
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
<div class="datatables html-widget html-fill-item-overflow-hidden html-fill-item" id="htmlwidget-1bd9d056b465896c6528" style="width:100%;height:auto;"></div>
<script type="application/json" data-for="htmlwidget-1bd9d056b465896c6528">{"x":{"filter":"none","vertical":false,"data":[["1","2","3","4","5","6"],["bjaccard","bjaccard","bjaccard","aichinson","aichinson","aichinson"],["Reactor_Treatment_Dose","Residual","Total","Reactor_Treatment_Dose","Residual","Total"],[10,23,33,10,23,33],[4.961,1.785,6.746,3.566,1.276,4.842],[0.735,0.265,1,0.737,0.263,1],[6.392,null,null,6.429,null,null],[0.001,null,null,0.001,null,null]],"container":"<table class=\"display\">\n  <thead>\n    <tr>\n      <th> <\/th>\n      <th>Distance<\/th>\n      <th>terms<\/th>\n      <th>Df<\/th>\n      <th>SumOfSqs<\/th>\n      <th>R2<\/th>\n      <th>F<\/th>\n      <th>Pr(&gt;F)<\/th>\n    <\/tr>\n  <\/thead>\n<\/table>","options":{"columnDefs":[{"className":"dt-right","targets":[3,4,5,6,7]},{"orderable":false,"targets":0}],"order":[],"autoWidth":false,"orderClasses":false}},"evals":[],"jsHooks":[]}</script>
```

```r
tmp_df %>% 
  xlsx::write.xlsx(file = out_xlsx,
                   append = TRUE,
                   sheetName = "human1all-beta-d3d6",
                   showNA = TRUE)
```


```r
lapply(
  dist,
  FUN = phyloseq_adonis2,
  physeq = ps_tmp_forbdiv %>%
    subset_samples(Day_of_Treatment >= 3 & Day_of_Treatment <= 6),
  formula = paste0(c("Reactor_Treatment_Dose", "Day_of_Treatment"), collapse=" + "),
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
<div class="datatables html-widget html-fill-item-overflow-hidden html-fill-item" id="htmlwidget-87804f75c1e74d8bd947" style="width:100%;height:auto;"></div>
<script type="application/json" data-for="htmlwidget-87804f75c1e74d8bd947">{"x":{"filter":"none","vertical":false,"data":[["1","2","3","4","5","6","7","8"],["bjaccard","bjaccard","bjaccard","bjaccard","aichinson","aichinson","aichinson","aichinson"],["Reactor_Treatment_Dose","Day_of_Treatment","Residual","Total","Reactor_Treatment_Dose","Day_of_Treatment","Residual","Total"],[10,1,22,33,10,1,22,33],[4.93,0.143,1.643,6.746,3.543,0.095,1.181,4.842],[0.731,0.021,0.243,1,0.732,0.02,0.244,1],[6.604,1.909,null,null,6.601,1.766,null,null],[0.001,0.083,null,null,0.001,0.083,null,null]],"container":"<table class=\"display\">\n  <thead>\n    <tr>\n      <th> <\/th>\n      <th>Distance<\/th>\n      <th>terms<\/th>\n      <th>Df<\/th>\n      <th>SumOfSqs<\/th>\n      <th>R2<\/th>\n      <th>F<\/th>\n      <th>Pr(&gt;F)<\/th>\n    <\/tr>\n  <\/thead>\n<\/table>","options":{"columnDefs":[{"className":"dt-right","targets":[3,4,5,6,7]},{"orderable":false,"targets":0}],"order":[],"autoWidth":false,"orderClasses":false}},"evals":[],"jsHooks":[]}</script>
```

```r
tmp_df %>% 
  xlsx::write.xlsx(file = out_xlsx,
                   append = TRUE,
                   sheetName = "human1all-beta-d3d6-time-day",
                   showNA = TRUE)
```


```r
lapply(
  dist,
  FUN = physeq_pairwise_permanovas_adonis2,
  physeq = ps_tmp_forbdiv %>%
    subset_samples(Day_of_Treatment >= 3 & Day_of_Treatment <= 6),
  compare_header = "Reactor_Treatment_Dose",
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
## Set of permutations < 'minperm'. Generating entire set.
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
## Set of permutations < 'minperm'. Generating entire set.
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

```r
tmp_df %>% 
  DT::datatable()
```

```{=html}
<div class="datatables html-widget html-fill-item-overflow-hidden html-fill-item" id="htmlwidget-7e2807b8c833e042a701" style="width:100%;height:auto;"></div>
<script type="application/json" data-for="htmlwidget-7e2807b8c833e042a701">{"x":{"filter":"none","vertical":false,"data":[["1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","23","24","25","26","27","28","29","30","31","32","33","34","35","36","37","38","39","40","41","42","43","44","45","46","47","48","49","50","51","52","53","54","55","56","57","58","59","60","61","62","63","64","65","66","67","68","69","70","71","72","73","74","75","76","77","78","79","80","81","82","83","84","85","86","87","88","89","90","91","92","93","94","95","96","97","98","99","100","101","102","103","104","105","106","107","108","109","110"],["bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson"],["CR_UNTREATED","CR_UNTREATED","CR_UNTREATED","CR_UNTREATED","CR_UNTREATED","CR_UNTREATED","CR_UNTREATED","CR_UNTREATED","CR_UNTREATED","CR_UNTREATED","TR2_CTX20","TR2_CTX20","TR2_CTX20","TR2_CTX20","TR2_CTX20","TR2_CTX20","TR2_CTX20","TR2_CTX20","TR2_CTX20","TR4_CTX200","TR4_CTX200","TR4_CTX200","TR4_CTX200","TR4_CTX200","TR4_CTX200","TR4_CTX200","TR4_CTX200","TR1_CTX+HV292.120","TR1_CTX+HV292.120","TR1_CTX+HV292.120","TR1_CTX+HV292.120","TR1_CTX+HV292.120","TR1_CTX+HV292.120","TR1_CTX+HV292.120","TR3_CTX+HV292.1200","TR3_CTX+HV292.1200","TR3_CTX+HV292.1200","TR3_CTX+HV292.1200","TR3_CTX+HV292.1200","TR3_CTX+HV292.1200","TR5_HV292.1","TR5_HV292.1","TR5_HV292.1","TR5_HV292.1","TR5_HV292.1","TR2_VAN90","TR2_VAN90","TR2_VAN90","TR2_VAN90","TR4_VAN600","TR4_VAN600","TR4_VAN600","TR1_VAN+CCUG5916890","TR1_VAN+CCUG5916890","TR3_VAN+CCUG59168600","CR_UNTREATED","CR_UNTREATED","CR_UNTREATED","CR_UNTREATED","CR_UNTREATED","CR_UNTREATED","CR_UNTREATED","CR_UNTREATED","CR_UNTREATED","CR_UNTREATED","TR2_CTX20","TR2_CTX20","TR2_CTX20","TR2_CTX20","TR2_CTX20","TR2_CTX20","TR2_CTX20","TR2_CTX20","TR2_CTX20","TR4_CTX200","TR4_CTX200","TR4_CTX200","TR4_CTX200","TR4_CTX200","TR4_CTX200","TR4_CTX200","TR4_CTX200","TR1_CTX+HV292.120","TR1_CTX+HV292.120","TR1_CTX+HV292.120","TR1_CTX+HV292.120","TR1_CTX+HV292.120","TR1_CTX+HV292.120","TR1_CTX+HV292.120","TR3_CTX+HV292.1200","TR3_CTX+HV292.1200","TR3_CTX+HV292.1200","TR3_CTX+HV292.1200","TR3_CTX+HV292.1200","TR3_CTX+HV292.1200","TR5_HV292.1","TR5_HV292.1","TR5_HV292.1","TR5_HV292.1","TR5_HV292.1","TR2_VAN90","TR2_VAN90","TR2_VAN90","TR2_VAN90","TR4_VAN600","TR4_VAN600","TR4_VAN600","TR1_VAN+CCUG5916890","TR1_VAN+CCUG5916890","TR3_VAN+CCUG59168600"],["TR2_CTX20","TR4_CTX200","TR1_CTX+HV292.120","TR3_CTX+HV292.1200","TR5_HV292.1","TR2_VAN90","TR4_VAN600","TR1_VAN+CCUG5916890","TR3_VAN+CCUG59168600","TR5_CCUG59168","TR4_CTX200","TR1_CTX+HV292.120","TR3_CTX+HV292.1200","TR5_HV292.1","TR2_VAN90","TR4_VAN600","TR1_VAN+CCUG5916890","TR3_VAN+CCUG59168600","TR5_CCUG59168","TR1_CTX+HV292.120","TR3_CTX+HV292.1200","TR5_HV292.1","TR2_VAN90","TR4_VAN600","TR1_VAN+CCUG5916890","TR3_VAN+CCUG59168600","TR5_CCUG59168","TR3_CTX+HV292.1200","TR5_HV292.1","TR2_VAN90","TR4_VAN600","TR1_VAN+CCUG5916890","TR3_VAN+CCUG59168600","TR5_CCUG59168","TR5_HV292.1","TR2_VAN90","TR4_VAN600","TR1_VAN+CCUG5916890","TR3_VAN+CCUG59168600","TR5_CCUG59168","TR2_VAN90","TR4_VAN600","TR1_VAN+CCUG5916890","TR3_VAN+CCUG59168600","TR5_CCUG59168","TR4_VAN600","TR1_VAN+CCUG5916890","TR3_VAN+CCUG59168600","TR5_CCUG59168","TR1_VAN+CCUG5916890","TR3_VAN+CCUG59168600","TR5_CCUG59168","TR3_VAN+CCUG59168600","TR5_CCUG59168","TR5_CCUG59168","TR2_CTX20","TR4_CTX200","TR1_CTX+HV292.120","TR3_CTX+HV292.1200","TR5_HV292.1","TR2_VAN90","TR4_VAN600","TR1_VAN+CCUG5916890","TR3_VAN+CCUG59168600","TR5_CCUG59168","TR4_CTX200","TR1_CTX+HV292.120","TR3_CTX+HV292.1200","TR5_HV292.1","TR2_VAN90","TR4_VAN600","TR1_VAN+CCUG5916890","TR3_VAN+CCUG59168600","TR5_CCUG59168","TR1_CTX+HV292.120","TR3_CTX+HV292.1200","TR5_HV292.1","TR2_VAN90","TR4_VAN600","TR1_VAN+CCUG5916890","TR3_VAN+CCUG59168600","TR5_CCUG59168","TR3_CTX+HV292.1200","TR5_HV292.1","TR2_VAN90","TR4_VAN600","TR1_VAN+CCUG5916890","TR3_VAN+CCUG59168600","TR5_CCUG59168","TR5_HV292.1","TR2_VAN90","TR4_VAN600","TR1_VAN+CCUG5916890","TR3_VAN+CCUG59168600","TR5_CCUG59168","TR2_VAN90","TR4_VAN600","TR1_VAN+CCUG5916890","TR3_VAN+CCUG59168600","TR5_CCUG59168","TR4_VAN600","TR1_VAN+CCUG5916890","TR3_VAN+CCUG59168600","TR5_CCUG59168","TR1_VAN+CCUG5916890","TR3_VAN+CCUG59168600","TR5_CCUG59168","TR3_VAN+CCUG59168600","TR5_CCUG59168","TR5_CCUG59168"],[0.428,0.33,0.356,0.47,0.175,0.681,0.489,0.752,0.709,0.173,0.267,0.325,0.52,0.67,0.753,0.507,0.846,0.841,0.712,0.212,0.189,0.468,0.607,0.355,0.71,0.703,0.545,0.479,0.788,0.721,0.426,0.823,0.832,0.754,0.584,0.714,0.474,0.798,0.794,0.714,0.645,0.315,0.771,0.806,0.763,0.286,0.346,0.588,0.777,0.335,0.326,0.486,0.634,0.862,0.863,0.377,0.289,0.325,0.392,0.176,0.573,0.51,0.631,0.573,0.208,0.311,0.368,0.517,0.638,0.754,0.619,0.821,0.833,0.659,0.311,0.222,0.509,0.604,0.444,0.683,0.697,0.574,0.511,0.798,0.773,0.578,0.841,0.883,0.718,0.597,0.701,0.56,0.767,0.784,0.696,0.74,0.478,0.83,0.929,0.698,0.282,0.339,0.601,0.791,0.328,0.322,0.62,0.668,0.852,0.868],[0.004,0.025,0.028,0.011,0.243,0.007,0.002,0.009,0.025,0.111,0.333,0.067,0.026,0.2,0.03,0.036,0.028,0.067,0.027,1,0.7,0.667,0.1,0.067,0.1,0.333,0.1,0.1,0.333,0.1,0.067,0.1,0.333,0.1,0.25,0.1,0.025,0.1,0.1,0.1,0.25,0.4,0.25,0.333,0.25,0.05,0.1,0.1,0.1,0.033,0.067,0.028,0.1,0.1,0.1,0.002,0.024,0.035,0.007,0.233,0.015,0.002,0.005,0.023,0.024,0.2,0.067,0.041,0.2,0.037,0.031,0.022,0.067,0.024,1,0.8,0.667,0.1,0.067,0.1,0.333,0.1,0.1,0.333,0.1,0.067,0.1,0.333,0.1,0.25,0.1,0.034,0.1,0.1,0.1,0.25,0.2,0.25,0.333,0.25,0.056,0.1,0.1,0.1,0.03,0.133,0.032,0.1,0.1,0.1],[0.11,0.196,0.128,0.121,0.318,0.128,0.11,0.124,0.196,0.153,0.374,0.183,0.159,0.268,0.118,0.124,0.128,0.183,0.148,1,0.713,0.692,0.177,0.183,0.177,0.374,0.177,0.177,0.374,0.177,0.183,0.177,0.374,0.177,0.309,0.177,0.196,0.177,0.177,0.177,0.309,0.423,0.309,0.374,0.309,0.162,0.177,0.177,0.177,0.121,0.183,0.128,0.177,0.177,0.177,0.073,0.147,0.128,0.096,0.291,0.165,0.073,0.092,0.181,0.147,0.262,0.179,0.133,0.262,0.127,0.142,0.202,0.179,0.147,1,0.815,0.692,0.177,0.179,0.177,0.363,0.177,0.177,0.363,0.177,0.179,0.177,0.363,0.177,0.296,0.177,0.134,0.177,0.177,0.177,0.296,0.262,0.296,0.363,0.296,0.171,0.177,0.177,0.177,0.15,0.183,0.135,0.177,0.177,0.177]],"container":"<table class=\"display\">\n  <thead>\n    <tr>\n      <th> <\/th>\n      <th>Distance<\/th>\n      <th>X1<\/th>\n      <th>X2<\/th>\n      <th>R2<\/th>\n      <th>pval<\/th>\n      <th>pvalFDR<\/th>\n    <\/tr>\n  <\/thead>\n<\/table>","options":{"columnDefs":[{"className":"dt-right","targets":[4,5,6]},{"orderable":false,"targets":0}],"order":[],"autoWidth":false,"orderClasses":false}},"evals":[],"jsHooks":[]}</script>
```

```r
tmp_df %>% 
  xlsx::write.xlsx(file = out_xlsx,
                   append = TRUE,
                   sheetName = "human1all-beta-d3d6-pw",
                   showNA = TRUE)
```


```r
lapply(
  dist,
  FUN = phyloseq_TW,
  physeq =  ps_tmp_forbdiv %>%
    subset_samples(Day_of_Treatment >= 3 & Day_of_Treatment <= 6),
  variable = "Reactor_Treatment_Dose",
  nrep = 999,
  strata = NULL
) %>%
  bind_rows(.id = "Distance") %>%
  mutate_if(is.numeric, round, 3)   -> tmp_df

tmp_df %>% 
  DT::datatable()
```

```{=html}
<div class="datatables html-widget html-fill-item-overflow-hidden html-fill-item" id="htmlwidget-40c8f686cc1c13cd7b4f" style="width:100%;height:auto;"></div>
<script type="application/json" data-for="htmlwidget-40c8f686cc1c13cd7b4f">{"x":{"filter":"none","vertical":false,"data":[["1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","23","24","25","26","27","28","29","30","31","32","33","34","35","36","37","38","39","40","41","42","43","44","45","46","47","48","49","50","51","52","53","54","55","56","57","58","59","60","61","62","63","64","65","66","67","68","69","70","71","72","73","74","75","76","77","78","79","80","81","82","83","84","85","86","87","88","89","90","91","92","93","94","95","96","97","98","99","100","101","102","103","104","105","106","107","108","109","110"],["bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson"],[["CR_UNTREATED"],["CR_UNTREATED"],["CR_UNTREATED"],["CR_UNTREATED"],["CR_UNTREATED"],["CR_UNTREATED"],["CR_UNTREATED"],["CR_UNTREATED"],["CR_UNTREATED"],["CR_UNTREATED"],["TR2_CTX20"],["TR2_CTX20"],["TR2_CTX20"],["TR2_CTX20"],["TR2_CTX20"],["TR2_CTX20"],["TR2_CTX20"],["TR2_CTX20"],["TR2_CTX20"],["TR4_CTX200"],["TR4_CTX200"],["TR4_CTX200"],["TR4_CTX200"],["TR4_CTX200"],["TR4_CTX200"],["TR4_CTX200"],["TR4_CTX200"],["TR1_CTX+HV292.120"],["TR1_CTX+HV292.120"],["TR1_CTX+HV292.120"],["TR1_CTX+HV292.120"],["TR1_CTX+HV292.120"],["TR1_CTX+HV292.120"],["TR1_CTX+HV292.120"],["TR3_CTX+HV292.1200"],["TR3_CTX+HV292.1200"],["TR3_CTX+HV292.1200"],["TR3_CTX+HV292.1200"],["TR3_CTX+HV292.1200"],["TR3_CTX+HV292.1200"],["TR5_HV292.1"],["TR5_HV292.1"],["TR5_HV292.1"],["TR5_HV292.1"],["TR5_HV292.1"],["TR2_VAN90"],["TR2_VAN90"],["TR2_VAN90"],["TR2_VAN90"],["TR4_VAN600"],["TR4_VAN600"],["TR4_VAN600"],["TR1_VAN+CCUG5916890"],["TR1_VAN+CCUG5916890"],["TR3_VAN+CCUG59168600"],["CR_UNTREATED"],["CR_UNTREATED"],["CR_UNTREATED"],["CR_UNTREATED"],["CR_UNTREATED"],["CR_UNTREATED"],["CR_UNTREATED"],["CR_UNTREATED"],["CR_UNTREATED"],["CR_UNTREATED"],["TR2_CTX20"],["TR2_CTX20"],["TR2_CTX20"],["TR2_CTX20"],["TR2_CTX20"],["TR2_CTX20"],["TR2_CTX20"],["TR2_CTX20"],["TR2_CTX20"],["TR4_CTX200"],["TR4_CTX200"],["TR4_CTX200"],["TR4_CTX200"],["TR4_CTX200"],["TR4_CTX200"],["TR4_CTX200"],["TR4_CTX200"],["TR1_CTX+HV292.120"],["TR1_CTX+HV292.120"],["TR1_CTX+HV292.120"],["TR1_CTX+HV292.120"],["TR1_CTX+HV292.120"],["TR1_CTX+HV292.120"],["TR1_CTX+HV292.120"],["TR3_CTX+HV292.1200"],["TR3_CTX+HV292.1200"],["TR3_CTX+HV292.1200"],["TR3_CTX+HV292.1200"],["TR3_CTX+HV292.1200"],["TR3_CTX+HV292.1200"],["TR5_HV292.1"],["TR5_HV292.1"],["TR5_HV292.1"],["TR5_HV292.1"],["TR5_HV292.1"],["TR2_VAN90"],["TR2_VAN90"],["TR2_VAN90"],["TR2_VAN90"],["TR4_VAN600"],["TR4_VAN600"],["TR4_VAN600"],["TR1_VAN+CCUG5916890"],["TR1_VAN+CCUG5916890"],["TR3_VAN+CCUG59168600"]],[["TR2_CTX20"],["TR4_CTX200"],["TR1_CTX+HV292.120"],["TR3_CTX+HV292.1200"],["TR5_HV292.1"],["TR2_VAN90"],["TR4_VAN600"],["TR1_VAN+CCUG5916890"],["TR3_VAN+CCUG59168600"],["TR5_CCUG59168"],["TR4_CTX200"],["TR1_CTX+HV292.120"],["TR3_CTX+HV292.1200"],["TR5_HV292.1"],["TR2_VAN90"],["TR4_VAN600"],["TR1_VAN+CCUG5916890"],["TR3_VAN+CCUG59168600"],["TR5_CCUG59168"],["TR1_CTX+HV292.120"],["TR3_CTX+HV292.1200"],["TR5_HV292.1"],["TR2_VAN90"],["TR4_VAN600"],["TR1_VAN+CCUG5916890"],["TR3_VAN+CCUG59168600"],["TR5_CCUG59168"],["TR3_CTX+HV292.1200"],["TR5_HV292.1"],["TR2_VAN90"],["TR4_VAN600"],["TR1_VAN+CCUG5916890"],["TR3_VAN+CCUG59168600"],["TR5_CCUG59168"],["TR5_HV292.1"],["TR2_VAN90"],["TR4_VAN600"],["TR1_VAN+CCUG5916890"],["TR3_VAN+CCUG59168600"],["TR5_CCUG59168"],["TR2_VAN90"],["TR4_VAN600"],["TR1_VAN+CCUG5916890"],["TR3_VAN+CCUG59168600"],["TR5_CCUG59168"],["TR4_VAN600"],["TR1_VAN+CCUG5916890"],["TR3_VAN+CCUG59168600"],["TR5_CCUG59168"],["TR1_VAN+CCUG5916890"],["TR3_VAN+CCUG59168600"],["TR5_CCUG59168"],["TR3_VAN+CCUG59168600"],["TR5_CCUG59168"],["TR5_CCUG59168"],["TR2_CTX20"],["TR4_CTX200"],["TR1_CTX+HV292.120"],["TR3_CTX+HV292.1200"],["TR5_HV292.1"],["TR2_VAN90"],["TR4_VAN600"],["TR1_VAN+CCUG5916890"],["TR3_VAN+CCUG59168600"],["TR5_CCUG59168"],["TR4_CTX200"],["TR1_CTX+HV292.120"],["TR3_CTX+HV292.1200"],["TR5_HV292.1"],["TR2_VAN90"],["TR4_VAN600"],["TR1_VAN+CCUG5916890"],["TR3_VAN+CCUG59168600"],["TR5_CCUG59168"],["TR1_CTX+HV292.120"],["TR3_CTX+HV292.1200"],["TR5_HV292.1"],["TR2_VAN90"],["TR4_VAN600"],["TR1_VAN+CCUG5916890"],["TR3_VAN+CCUG59168600"],["TR5_CCUG59168"],["TR3_CTX+HV292.1200"],["TR5_HV292.1"],["TR2_VAN90"],["TR4_VAN600"],["TR1_VAN+CCUG5916890"],["TR3_VAN+CCUG59168600"],["TR5_CCUG59168"],["TR5_HV292.1"],["TR2_VAN90"],["TR4_VAN600"],["TR1_VAN+CCUG5916890"],["TR3_VAN+CCUG59168600"],["TR5_CCUG59168"],["TR2_VAN90"],["TR4_VAN600"],["TR1_VAN+CCUG5916890"],["TR3_VAN+CCUG59168600"],["TR5_CCUG59168"],["TR4_VAN600"],["TR1_VAN+CCUG5916890"],["TR3_VAN+CCUG59168600"],["TR5_CCUG59168"],["TR1_VAN+CCUG5916890"],["TR3_VAN+CCUG59168600"],["TR5_CCUG59168"],["TR3_VAN+CCUG59168600"],["TR5_CCUG59168"],["TR5_CCUG59168"]],[[7],[7],[7],[7],[7],[7],[7],[7],[7],[7],[4],[4],[4],[4],[4],[4],[4],[4],[4],[2],[2],[2],[2],[2],[2],[2],[2],[2],[2],[2],[2],[2],[2],[2],[3],[3],[3],[3],[3],[3],[1],[1],[1],[1],[1],[3],[3],[3],[3],[4],[4],[4],[3],[3],[2],[7],[7],[7],[7],[7],[7],[7],[7],[7],[7],[4],[4],[4],[4],[4],[4],[4],[4],[4],[2],[2],[2],[2],[2],[2],[2],[2],[2],[2],[2],[2],[2],[2],[2],[3],[3],[3],[3],[3],[3],[1],[1],[1],[1],[1],[3],[3],[3],[3],[4],[4],[4],[3],[3],[2]],[[4],[2],[2],[3],[1],[3],[4],[3],[2],[3],[2],[2],[3],[1],[3],[4],[3],[2],[3],[2],[3],[1],[3],[4],[3],[2],[3],[3],[1],[3],[4],[3],[2],[3],[1],[3],[4],[3],[2],[3],[3],[4],[3],[2],[3],[4],[3],[2],[3],[3],[2],[3],[2],[3],[3],[4],[2],[2],[3],[1],[3],[4],[3],[2],[3],[2],[2],[3],[1],[3],[4],[3],[2],[3],[2],[3],[1],[3],[4],[3],[2],[3],[3],[1],[3],[4],[3],[2],[3],[1],[3],[4],[3],[2],[3],[3],[4],[3],[2],[3],[4],[3],[2],[3],[3],[2],[3],[2],[3],[3]],[[0.004],[0.223],[0.037],[0.005],[null],[0.012],[0.002],[0.006],[0.04],[0.052],[0.795],[0.275],[0.032],[null],[0.025],[0.03],[0.034],[0.068],[0.029],[1],[1],[null],[0.111],[0.082],[0.088],[0.325],[0.1],[0.096],[null],[0.108],[0.085],[0.109],[0.331],[0.083],[null],[0.114],[0.032],[0.091],[0.106],[0.097],[null],[null],[null],[null],[null],[0.053],[0.095],[0.099],[0.099],[0.034],[0.076],[0.03],[0.102],[0.091],[0.087],[0.004],[0.19],[0.036],[0.01],[null],[0.015],[0.006],[0.006],[0.035],[0.032],[0.387],[0.073],[0.031],[null],[0.033],[0.029],[0.032],[0.057],[0.032],[1],[1],[null],[0.093],[0.071],[0.108],[0.356],[0.1],[0.11],[null],[0.115],[0.067],[0.103],[0.342],[0.089],[null],[0.117],[0.031],[0.097],[0.093],[0.093],[null],[null],[null],[null],[null],[0.061],[0.086],[0.096],[0.103],[0.027],[0.078],[0.035],[0.093],[0.102],[0.116]],[[8.54685159255591],[1.77792831839827],[4.84468862499253],[7.08766033263965],[null],[12.097732218186],[5.67323700880519],[21.1636610725694],[10.9198916583527],[2.80887730939864],[0.736396214173821],[1.57310508691593],[4.7169398804577],[null],[12.0776496282593],[6.16113033999895],[22.9301727091681],[11.7329516879107],[13.3573998879747],[0.537611500691741],[0.535883343945286],[null],[4.30057315166948],[2.74501949098503],[6.04790217873442],[4.72988936759038],[2.22816785242924],[3.0273654786785],[null],[10.1535077941617],[5.46407075346997],[16.4778327447803],[9.89477955881093],[7.19015736308848],[null],[9.99789866516812],[5.64540876263235],[15.8494600598794],[9.65051780616452],[9.99250569371225],[null],[null],[null],[null],[null],[2.28650188568313],[2.11947941949056],[4.37106694797],[13.9117673336556],[3.05979519097538],[2.73130127793387],[6.41060295126254],[4.69166398524375],[24.9580955917854],[12.1941933703433],[6.99308703478946],[2.00553151506417],[5.03029441558812],[6.04697039102683],[null],[12.7694481314372],[8.50879070172805],[19.4684384947428],[18.8214232881495],[3.03704813886125],[1.03763851760866],[2.12228248174519],[4.9067895512456],[null],[14.0725291592948],[9.73582667917857],[22.7881628928415],[22.9544713764173],[9.67517981203025],[0.903560520583048],[0.678124087556537],[null],[3.59251440823115],[2.91434984765577],[4.63183081044014],[4.59457675258071],[2.86874580284958],[3.39241553536459],[null],[10.9465578868058],[7.94245495968684],[15.1415783373536],[15.1576262508895],[7.20808332313827],[null],[9.39979526537332],[7.11948539359265],[13.1687028328723],[13.4560711731494],[9.14987026897994],[null],[null],[null],[null],[null],[2.20404979513458],[2.05016651208964],[5.5480106254557],[15.171644886578],[2.91736678158837],[3.22110903170035],[9.81165621138822],[6.69866439906215],[22.9437718724269],[21.6529736494397]],[[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999]],[0.068,0.257,0.119,0.068,null,0.108,0.068,0.068,0.12,0.135,0.832,0.309,0.118,null,0.118,0.118,0.118,0.135,0.118,1,1,null,0.135,0.135,0.135,0.355,0.135,0.135,null,0.135,0.135,0.135,0.355,0.135,null,0.135,0.118,0.135,0.135,0.135,null,null,null,null,null,0.135,0.135,0.135,0.135,0.118,0.135,0.118,0.135,0.135,0.135,0.09,0.214,0.101,0.101,null,0.101,0.09,0.09,0.101,0.101,0.405,0.135,0.101,null,0.101,0.101,0.101,0.135,0.101,1,1,null,0.135,0.135,0.135,0.381,0.135,0.135,null,0.135,0.135,0.135,0.375,0.135,null,0.135,0.101,0.135,0.135,0.135,null,null,null,null,null,0.135,0.135,0.135,0.135,0.101,0.135,0.101,0.135,0.135,0.135]],"container":"<table class=\"display\">\n  <thead>\n    <tr>\n      <th> <\/th>\n      <th>Distance<\/th>\n      <th>Level1<\/th>\n      <th>Level2<\/th>\n      <th>N1<\/th>\n      <th>N2<\/th>\n      <th>p.value<\/th>\n      <th>tw2.stat<\/th>\n      <th>nrep<\/th>\n      <th>pvalFDR<\/th>\n    <\/tr>\n  <\/thead>\n<\/table>","options":{"columnDefs":[{"className":"dt-right","targets":9},{"orderable":false,"targets":0}],"order":[],"autoWidth":false,"orderClasses":false}},"evals":[],"jsHooks":[]}</script>
```

```r
tmp_df %>% 
  xlsx::write.xlsx(file = out_xlsx,
                   append = TRUE,
                   sheetName = "human1all-beta-d3d6-pw2",
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


##### Period:


```r
setNames(vector("list", 
                length(physeq_get_unique(ps_tmp_forbdiv, "Reactor_Treatment_Dose"))), 
         physeq_get_unique(ps_tmp_forbdiv, "Reactor_Treatment_Dose")) -> lis
```



```r
for (grp in names(lis)){
  print(grp)
  
  lapply(
    dist,
    FUN = physeq_pairwise_permanovas_adonis2,
    physeq = ps_tmp_forbdiv %>%
      # subset_samples(.,  eval(as.name(ranks)) %in% tax) %>% 
      subset_samples(., Reactor_Treatment_Dose == grp),
    compare_header = "period",
    n_perm = 999,
    strat = "none"
  ) %>%
    bind_rows(.id = "Distance") %>%
    mutate_if(is.numeric, round, 3) %>%
    select(-pvalBon) -> lis[[grp]]
}
```

```
## [1] "TR4_VAN600"
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
```

```
## [1] "TR2_VAN90"
```

```
## 'nperm' >= set of all permutations: complete enumeration.
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
## [1] "TR3_VAN+CCUG59168600"
```

```
## 'nperm' >= set of all permutations: complete enumeration.
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
## [1] "TR2_CTX20"
```

```
## 'nperm' >= set of all permutations: complete enumeration.
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
## [1] "TR1_VAN+CCUG5916890"
```

```
## 'nperm' >= set of all permutations: complete enumeration.
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
## [1] "TR4_CTX200"
```

```
## 'nperm' >= set of all permutations: complete enumeration.
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
## [1] "CR_UNTREATED"
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
## [1] "TR5_CCUG59168"
```

```
## 'nperm' >= set of all permutations: complete enumeration.
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
## [1] "TR5_HV292.1"
```

```
## 'nperm' >= set of all permutations: complete enumeration.
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
```

```
## 'nperm' >= set of all permutations: complete enumeration.
```

```
## Set of permutations < 'minperm'. Generating entire set.
```

```
## [1] "TR1_CTX+HV292.120"
```

```
## 'nperm' >= set of all permutations: complete enumeration.
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
```

```
## 'nperm' >= set of all permutations: complete enumeration.
```

```
## Set of permutations < 'minperm'. Generating entire set.
```

```
## [1] "TR3_CTX+HV292.1200"
```

```
## 'nperm' >= set of all permutations: complete enumeration.
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

```r
lis %>% 
  bind_rows(.id = "Reactor_Treatment_Dose") %>% 
  filter(X1 == "pret" | "X2" == "pret") %>% 
  group_by(Reactor_Treatment_Dose, Distance) %>% 
  rstatix::adjust_pvalue(p.col = "pval", output.col = "pvalFDR") %>% 
  rstatix::add_significance(p.col = "pvalFDR")-> tmp_df

tmp_df %>% 
  DT::datatable()
```

```{=html}
<div class="datatables html-widget html-fill-item-overflow-hidden html-fill-item" id="htmlwidget-06e32c4b185a8cd3895d" style="width:100%;height:auto;"></div>
<script type="application/json" data-for="htmlwidget-06e32c4b185a8cd3895d">{"x":{"filter":"none","vertical":false,"data":[["1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","23","24","25","26","27","28","29","30","31","32","33","34","35","36","37","38","39","40","41","42","43","44","45","46","47","48","49","50","51","52","53","54","55","56","57","58","59","60","61","62"],["CR_UNTREATED","CR_UNTREATED","CR_UNTREATED","CR_UNTREATED","CR_UNTREATED","CR_UNTREATED","TR1_CTX+HV292.120","TR1_CTX+HV292.120","TR1_CTX+HV292.120","TR1_CTX+HV292.120","TR1_VAN+CCUG5916890","TR1_VAN+CCUG5916890","TR1_VAN+CCUG5916890","TR1_VAN+CCUG5916890","TR1_VAN+CCUG5916890","TR1_VAN+CCUG5916890","TR2_CTX20","TR2_CTX20","TR2_CTX20","TR2_CTX20","TR2_CTX20","TR2_CTX20","TR2_VAN90","TR2_VAN90","TR2_VAN90","TR2_VAN90","TR2_VAN90","TR2_VAN90","TR3_CTX+HV292.1200","TR3_CTX+HV292.1200","TR3_CTX+HV292.1200","TR3_CTX+HV292.1200","TR3_CTX+HV292.1200","TR3_CTX+HV292.1200","TR3_VAN+CCUG59168600","TR3_VAN+CCUG59168600","TR3_VAN+CCUG59168600","TR3_VAN+CCUG59168600","TR3_VAN+CCUG59168600","TR3_VAN+CCUG59168600","TR4_CTX200","TR4_CTX200","TR4_CTX200","TR4_CTX200","TR4_CTX200","TR4_CTX200","TR4_VAN600","TR4_VAN600","TR4_VAN600","TR4_VAN600","TR4_VAN600","TR4_VAN600","TR5_CCUG59168","TR5_CCUG59168","TR5_CCUG59168","TR5_CCUG59168","TR5_CCUG59168","TR5_CCUG59168","TR5_HV292.1","TR5_HV292.1","TR5_HV292.1","TR5_HV292.1"],["aichinson","aichinson","aichinson","bjaccard","bjaccard","bjaccard","aichinson","aichinson","bjaccard","bjaccard","aichinson","aichinson","aichinson","bjaccard","bjaccard","bjaccard","aichinson","aichinson","aichinson","bjaccard","bjaccard","bjaccard","aichinson","aichinson","aichinson","bjaccard","bjaccard","bjaccard","aichinson","aichinson","aichinson","bjaccard","bjaccard","bjaccard","aichinson","aichinson","aichinson","bjaccard","bjaccard","bjaccard","aichinson","aichinson","aichinson","bjaccard","bjaccard","bjaccard","aichinson","aichinson","aichinson","bjaccard","bjaccard","bjaccard","aichinson","aichinson","aichinson","bjaccard","bjaccard","bjaccard","aichinson","aichinson","bjaccard","bjaccard"],["pret","pret","pret","pret","pret","pret","pret","pret","pret","pret","pret","pret","pret","pret","pret","pret","pret","pret","pret","pret","pret","pret","pret","pret","pret","pret","pret","pret","pret","pret","pret","pret","pret","pret","pret","pret","pret","pret","pret","pret","pret","pret","pret","pret","pret","pret","pret","pret","pret","pret","pret","pret","pret","pret","pret","pret","pret","pret","pret","pret","pret","pret"],["t1","t2","t3","t1","t2","t3","t2","t3","t2","t3","t1","t2","t3","t1","t2","t3","t1","t2","t3","t1","t2","t3","t1","t2","t3","t1","t2","t3","t1","t2","t3","t1","t2","t3","t1","t2","t3","t1","t2","t3","t1","t2","t3","t1","t2","t3","t1","t2","t3","t1","t2","t3","t1","t2","t3","t1","t2","t3","t1","t2","t1","t2"],[0.113,0.103,0.192,0.096,0.095,0.217,0.551,0.546,0.559,0.604,0.657,0.841,0.829,0.718,0.855,0.912,0.439,0.612,0.802,0.45,0.691,0.921,0.606,0.727,0.851,0.689,0.711,0.924,0.626,0.668,0.793,0.641,0.696,0.887,0.633,0.771,0.737,0.679,0.805,0.868,0.416,0.455,0.726,0.419,0.468,0.801,0.602,0.584,0.912,0.609,0.422,0.971,0.411,0.454,0.442,0.36,0.43,0.39,0.367,0.716,0.464,0.786],[0.377,0.201,0.262,0.574,0.295,0.161,0.1,0.25,0.1,0.25,0.2,0.028,0.2,0.2,0.036,0.2,0.1,0.026,0.25,0.1,0.025,0.25,0.333,0.1,0.333,0.333,0.1,0.333,0.1,0.1,0.25,0.1,0.1,0.25,0.067,0.067,0.2,0.067,0.067,0.2,0.333,0.333,0.333,0.333,0.333,0.333,0.333,0.067,0.333,0.333,0.067,0.333,0.1,0.1,0.5,0.2,0.1,0.5,0.25,0.25,0.25,0.25],[0.603,0.603,0.603,0.59,0.59,0.483,0.2,0.25,0.2,0.25,0.4,0.084,0.4,0.4,0.108,0.4,0.2,0.078,0.25,0.2,0.075,0.25,0.666,0.3,0.666,0.666,0.3,0.666,0.3,0.3,0.3,0.3,0.3,0.3,0.201,0.201,0.201,0.201,0.201,0.201,0.999,0.999,0.999,0.999,0.999,0.999,0.666,0.201,0.666,0.666,0.201,0.666,0.3,0.3,0.5,0.4,0.3,0.5,0.5,0.5,0.5,0.5],["ns","ns","ns","ns","ns","ns","ns","ns","ns","ns","ns","ns","ns","ns","ns","ns","ns","ns","ns","ns","ns","ns","ns","ns","ns","ns","ns","ns","ns","ns","ns","ns","ns","ns","ns","ns","ns","ns","ns","ns","ns","ns","ns","ns","ns","ns","ns","ns","ns","ns","ns","ns","ns","ns","ns","ns","ns","ns","ns","ns","ns","ns"]],"container":"<table class=\"display\">\n  <thead>\n    <tr>\n      <th> <\/th>\n      <th>Reactor_Treatment_Dose<\/th>\n      <th>Distance<\/th>\n      <th>X1<\/th>\n      <th>X2<\/th>\n      <th>R2<\/th>\n      <th>pval<\/th>\n      <th>pvalFDR<\/th>\n      <th>pvalFDR.signif<\/th>\n    <\/tr>\n  <\/thead>\n<\/table>","options":{"columnDefs":[{"className":"dt-right","targets":[5,6,7]},{"orderable":false,"targets":0}],"order":[],"autoWidth":false,"orderClasses":false}},"evals":[],"jsHooks":[]}</script>
```

```r
tmp_df %>% 
  xlsx::write.xlsx(file = out_xlsx,
                   append = TRUE,
                   sheetName = "human1all-beta-vs_pret",
                   showNA = TRUE)
```



```r
lis %>% 
  bind_rows(.id = "Reactor_Treatment_Dose") %>% 
  group_by(Reactor_Treatment_Dose, Distance) %>% 
  rstatix::adjust_pvalue(p.col = "pval", output.col = "pvalFDR") %>% 
  rstatix::add_significance(p.col = "pvalFDR") -> tmp_df

tmp_df %>% 
  DT::datatable()
```

```{=html}
<div class="datatables html-widget html-fill-item-overflow-hidden html-fill-item" id="htmlwidget-7ea80b0746f80e52b0ea" style="width:100%;height:auto;"></div>
<script type="application/json" data-for="htmlwidget-7ea80b0746f80e52b0ea">{"x":{"filter":"none","vertical":false,"data":[["1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","23","24","25","26","27","28","29","30","31","32","33","34","35","36","37","38","39","40","41","42","43","44","45","46","47","48","49","50","51","52","53","54","55","56","57","58","59","60","61","62","63","64","65","66","67","68","69","70","71","72","73","74","75","76","77","78","79","80","81","82","83","84","85","86","87","88","89","90","91","92","93","94","95","96","97","98","99","100","101","102","103","104","105","106","107","108","109","110","111","112","113","114","115","116","117","118","119","120"],["CR_UNTREATED","CR_UNTREATED","CR_UNTREATED","CR_UNTREATED","CR_UNTREATED","CR_UNTREATED","CR_UNTREATED","CR_UNTREATED","CR_UNTREATED","CR_UNTREATED","CR_UNTREATED","CR_UNTREATED","TR1_CTX+HV292.120","TR1_CTX+HV292.120","TR1_CTX+HV292.120","TR1_CTX+HV292.120","TR1_CTX+HV292.120","TR1_CTX+HV292.120","TR1_VAN+CCUG5916890","TR1_VAN+CCUG5916890","TR1_VAN+CCUG5916890","TR1_VAN+CCUG5916890","TR1_VAN+CCUG5916890","TR1_VAN+CCUG5916890","TR1_VAN+CCUG5916890","TR1_VAN+CCUG5916890","TR1_VAN+CCUG5916890","TR1_VAN+CCUG5916890","TR1_VAN+CCUG5916890","TR1_VAN+CCUG5916890","TR2_CTX20","TR2_CTX20","TR2_CTX20","TR2_CTX20","TR2_CTX20","TR2_CTX20","TR2_CTX20","TR2_CTX20","TR2_CTX20","TR2_CTX20","TR2_CTX20","TR2_CTX20","TR2_VAN90","TR2_VAN90","TR2_VAN90","TR2_VAN90","TR2_VAN90","TR2_VAN90","TR2_VAN90","TR2_VAN90","TR2_VAN90","TR2_VAN90","TR2_VAN90","TR2_VAN90","TR3_CTX+HV292.1200","TR3_CTX+HV292.1200","TR3_CTX+HV292.1200","TR3_CTX+HV292.1200","TR3_CTX+HV292.1200","TR3_CTX+HV292.1200","TR3_CTX+HV292.1200","TR3_CTX+HV292.1200","TR3_CTX+HV292.1200","TR3_CTX+HV292.1200","TR3_CTX+HV292.1200","TR3_CTX+HV292.1200","TR3_VAN+CCUG59168600","TR3_VAN+CCUG59168600","TR3_VAN+CCUG59168600","TR3_VAN+CCUG59168600","TR3_VAN+CCUG59168600","TR3_VAN+CCUG59168600","TR3_VAN+CCUG59168600","TR3_VAN+CCUG59168600","TR3_VAN+CCUG59168600","TR3_VAN+CCUG59168600","TR3_VAN+CCUG59168600","TR3_VAN+CCUG59168600","TR4_CTX200","TR4_CTX200","TR4_CTX200","TR4_CTX200","TR4_CTX200","TR4_CTX200","TR4_CTX200","TR4_CTX200","TR4_CTX200","TR4_CTX200","TR4_CTX200","TR4_CTX200","TR4_VAN600","TR4_VAN600","TR4_VAN600","TR4_VAN600","TR4_VAN600","TR4_VAN600","TR4_VAN600","TR4_VAN600","TR4_VAN600","TR4_VAN600","TR4_VAN600","TR4_VAN600","TR5_CCUG59168","TR5_CCUG59168","TR5_CCUG59168","TR5_CCUG59168","TR5_CCUG59168","TR5_CCUG59168","TR5_CCUG59168","TR5_CCUG59168","TR5_CCUG59168","TR5_CCUG59168","TR5_CCUG59168","TR5_CCUG59168","TR5_HV292.1","TR5_HV292.1","TR5_HV292.1","TR5_HV292.1","TR5_HV292.1","TR5_HV292.1"],["aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","aichinson","aichinson","aichinson","bjaccard","bjaccard","bjaccard","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","aichinson","aichinson","aichinson","bjaccard","bjaccard","bjaccard"],["pret","pret","pret","t1","t1","t2","pret","pret","pret","t1","t1","t2","pret","pret","t2","pret","pret","t2","pret","pret","pret","t1","t1","t2","pret","pret","pret","t1","t1","t2","pret","pret","pret","t1","t1","t2","pret","pret","pret","t1","t1","t2","pret","pret","pret","t1","t1","t2","pret","pret","pret","t1","t1","t2","pret","pret","pret","t1","t1","t2","pret","pret","pret","t1","t1","t2","pret","pret","pret","t1","t1","t2","pret","pret","pret","t1","t1","t2","pret","pret","pret","t1","t1","t2","pret","pret","pret","t1","t1","t2","pret","pret","pret","t1","t1","t2","pret","pret","pret","t1","t1","t2","pret","pret","pret","t1","t1","t2","pret","pret","pret","t1","t1","t2","pret","pret","t1","pret","pret","t1"],["t1","t2","t3","t2","t3","t3","t1","t2","t3","t2","t3","t3","t2","t3","t3","t2","t3","t3","t1","t2","t3","t2","t3","t3","t1","t2","t3","t2","t3","t3","t1","t2","t3","t2","t3","t3","t1","t2","t3","t2","t3","t3","t1","t2","t3","t2","t3","t3","t1","t2","t3","t2","t3","t3","t1","t2","t3","t2","t3","t3","t1","t2","t3","t2","t3","t3","t1","t2","t3","t2","t3","t3","t1","t2","t3","t2","t3","t3","t1","t2","t3","t2","t3","t3","t1","t2","t3","t2","t3","t3","t1","t2","t3","t2","t3","t3","t1","t2","t3","t2","t3","t3","t1","t2","t3","t2","t3","t3","t1","t2","t3","t2","t3","t3","t1","t2","t2","t1","t2","t2"],[0.113,0.103,0.192,0.089,0.419,0.122,0.096,0.095,0.217,0.085,0.518,0.125,0.551,0.546,0.509,0.559,0.604,0.47,0.657,0.841,0.829,0.745,1,0.223,0.718,0.855,0.912,0.735,1,0.209,0.439,0.612,0.802,0.352,0.408,0.631,0.45,0.691,0.921,0.399,0.381,0.719,0.606,0.727,0.851,0.563,0.772,0.299,0.689,0.711,0.924,0.607,0.887,0.338,0.626,0.668,0.793,0.386,0.638,0.344,0.641,0.696,0.887,0.392,0.663,0.378,0.633,0.771,0.737,0.574,0.696,0.617,0.679,0.805,0.868,0.542,0.725,0.472,0.416,0.455,0.726,0.211,0.547,0.388,0.419,0.468,0.801,0.136,0.621,0.426,0.602,0.584,0.912,0.308,0.618,0.116,0.609,0.422,0.971,0.248,0.765,0.174,0.411,0.454,0.442,0.295,0.563,0.421,0.36,0.43,0.39,0.333,0.471,0.456,0.367,0.716,1,0.464,0.786,1],[0.377,0.201,0.262,0.619,0.25,0.634,0.574,0.295,0.161,0.666,0.25,0.486,0.1,0.25,0.667,0.1,0.25,0.667,0.2,0.028,0.2,0.25,null,1,0.2,0.036,0.2,0.25,null,1,0.1,0.026,0.25,0.067,1,0.2,0.1,0.025,0.25,0.067,1,0.2,0.333,0.1,0.333,0.1,0.333,0.75,0.333,0.1,0.333,0.1,0.333,0.5,0.1,0.1,0.25,0.1,0.333,0.5,0.1,0.1,0.25,0.1,0.333,0.5,0.067,0.067,0.2,0.333,0.333,0.333,0.067,0.067,0.2,0.333,0.333,0.667,0.333,0.333,0.333,1,0.333,1,0.333,0.333,0.333,1,0.333,1,0.333,0.067,0.333,0.133,0.667,0.6,0.333,0.067,0.333,0.267,0.333,0.6,0.1,0.1,0.5,0.1,0.333,0.25,0.2,0.1,0.5,0.2,0.667,0.25,0.25,0.25,null,0.25,0.25,null],[1,1,1,1,1,1,1,1,0.966,1,1,1,0.3,0.5,0.667,0.3,0.5,0.667,0.8,0.14,0.8,0.8,null,1,0.8,0.18,0.8,0.8,null,1,0.4,0.156,0.6,0.335,1,0.6,0.4,0.15,0.6,0.335,1,0.6,1,0.6,1,0.6,1,1,1,0.6,1,0.6,1,1,0.6,0.6,0.75,0.6,0.75,0.75,0.6,0.6,0.75,0.6,0.75,0.75,0.402,0.402,0.8,0.999,0.999,0.999,0.402,0.402,0.8,0.999,0.999,0.999,1,1,1,1,1,1,1,1,1,1,1,1,1,0.402,1,0.665,1,1,1,0.402,1,1,1,1,0.6,0.6,0.75,0.6,0.75,0.75,1,0.6,1,1,1,1,0.5,0.5,null,0.5,0.5,null],["ns","ns","ns","ns","ns","ns","ns","ns","ns","ns","ns","ns","ns","ns","ns","ns","ns","ns","ns","ns","ns","ns","","ns","ns","ns","ns","ns","","ns","ns","ns","ns","ns","ns","ns","ns","ns","ns","ns","ns","ns","ns","ns","ns","ns","ns","ns","ns","ns","ns","ns","ns","ns","ns","ns","ns","ns","ns","ns","ns","ns","ns","ns","ns","ns","ns","ns","ns","ns","ns","ns","ns","ns","ns","ns","ns","ns","ns","ns","ns","ns","ns","ns","ns","ns","ns","ns","ns","ns","ns","ns","ns","ns","ns","ns","ns","ns","ns","ns","ns","ns","ns","ns","ns","ns","ns","ns","ns","ns","ns","ns","ns","ns","ns","ns","","ns","ns",""]],"container":"<table class=\"display\">\n  <thead>\n    <tr>\n      <th> <\/th>\n      <th>Reactor_Treatment_Dose<\/th>\n      <th>Distance<\/th>\n      <th>X1<\/th>\n      <th>X2<\/th>\n      <th>R2<\/th>\n      <th>pval<\/th>\n      <th>pvalFDR<\/th>\n      <th>pvalFDR.signif<\/th>\n    <\/tr>\n  <\/thead>\n<\/table>","options":{"columnDefs":[{"className":"dt-right","targets":[5,6,7]},{"orderable":false,"targets":0}],"order":[],"autoWidth":false,"orderClasses":false}},"evals":[],"jsHooks":[]}</script>
```

```r
tmp_df %>% 
  xlsx::write.xlsx(file = out_xlsx,
                   append = TRUE,
                   sheetName = "human1all-beta-between-periods",
                   showNA = TRUE)
```



```r
# for (grp in names(lis)){
#   print(grp)
#   
# lapply(
#   dist,
#   FUN = phyloseq_TW,
#   physeq =  ps_tmp_forbdiv %>%
#       # subset_samples(.,  eval(as.name(ranks)) %in% tax) %>% 
#       subset_samples(., Reactor_Treatment_Dose == grp),
#     variable = "period",
#   nrep = 999,
#   strata = NULL
# ) %>%
#   bind_rows(.id = "Distance") %>%
#   mutate_if(is.numeric, round, 5) -> lis[[grp]]
# }
# 
# 
# lis %>% 
#   bind_rows(.id = "Reactor_Treatment_Dose") %>% 
#   xlsx::write.xlsx(file = out_xlsx,
#                    append = TRUE,
#                    sheetName = "human1all-beta-between-periodsTW",
#                    showNA = TRUE)
```






## End:


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
##  [1] vegan_2.6-4          lattice_0.20-45      permute_0.9-7       
##  [4] gdtools_0.2.4        GUniFrac_1.7         ape_5.6-2           
##  [7] speedyseq_0.5.3.9018 reshape2_1.4.4       scales_1.2.1        
## [10] compositions_2.0-4   nlme_3.1-161         phyloseq_1.42.0     
## [13] forcats_0.5.2        stringr_1.5.0        dplyr_1.0.10        
## [16] purrr_0.3.5          readr_2.1.3          tidyr_1.2.1         
## [19] tibble_3.1.8         ggplot2_3.4.0        tidyverse_1.3.2     
## 
## loaded via a namespace (and not attached):
##   [1] uuid_1.1-0             readxl_1.4.1           backports_1.4.1       
##   [4] systemfonts_1.0.4      plyr_1.8.8             igraph_1.3.5          
##   [7] splines_4.2.2          crosstalk_1.2.0        GenomeInfoDb_1.34.3   
##  [10] digest_0.6.31          foreach_1.5.2          htmltools_0.5.4       
##  [13] fansi_1.0.3            magrittr_2.0.3         xlsx_0.6.5            
##  [16] googlesheets4_1.0.1    cluster_2.1.4          openxlsx_4.2.5.1      
##  [19] tzdb_0.3.0             Biostrings_2.66.0      extrafont_0.18        
##  [22] modelr_0.1.10          bayesm_3.1-5           matrixStats_0.63.0    
##  [25] officer_0.5.0          stabledist_0.7-1       extrafontdb_1.0       
##  [28] timechange_0.1.1       colorspace_2.0-3       rvest_1.0.3           
##  [31] ggrepel_0.9.2          textshaping_0.3.6      haven_2.5.1           
##  [34] xfun_0.35              crayon_1.5.2           RCurl_1.98-1.9        
##  [37] jsonlite_1.8.4         survival_3.4-0         iterators_1.0.14      
##  [40] glue_1.6.2             rvg_0.3.0              gtable_0.3.1          
##  [43] gargle_1.2.1           zlibbioc_1.44.0        XVector_0.38.0        
##  [46] car_3.1-1              Rttf2pt1_1.3.11        Rhdf5lib_1.20.0       
##  [49] BiocGenerics_0.44.0    DEoptimR_1.0-11        abind_1.4-5           
##  [52] DBI_1.1.3              rstatix_0.7.1          Rcpp_1.0.9            
##  [55] xtable_1.8-4           clue_0.3-63            DT_0.26               
##  [58] stats4_4.2.2           htmlwidgets_1.6.0      timeSeries_4021.105   
##  [61] httr_1.4.4             ellipsis_0.3.2         spatial_7.3-15        
##  [64] rJava_1.0-6            pkgconfig_2.0.3        farver_2.1.1          
##  [67] sass_0.4.4             dbplyr_2.2.1           utf8_1.2.2            
##  [70] here_1.0.1             tidyselect_1.2.0       labeling_0.4.2        
##  [73] rlang_1.0.6            munsell_0.5.0          cellranger_1.1.0      
##  [76] tools_4.2.2            cachem_1.0.6           cli_3.4.1             
##  [79] generics_0.1.3         devEMF_4.1-1           ade4_1.7-20           
##  [82] export_0.3.0           broom_1.0.2            evaluate_0.19         
##  [85] biomformat_1.26.0      fastmap_1.1.0          yaml_2.3.6            
##  [88] ragg_1.2.4             knitr_1.41             fs_1.5.2              
##  [91] zip_2.2.2              robustbase_0.95-0      rgl_0.110.2           
##  [94] xml2_1.3.3             compiler_4.2.2         rstudioapi_0.14       
##  [97] ggsignif_0.6.4         reprex_2.0.2           statmod_1.4.37        
## [100] bslib_0.4.2            stringi_1.7.8          statip_0.2.3          
## [103] highr_0.9              stargazer_5.2.3        modeest_2.4.0         
## [106] fBasics_4021.93        Matrix_1.5-3           microbiome_1.19.1     
## [109] tensorA_0.36.2         multtest_2.54.0        vctrs_0.5.1           
## [112] pillar_1.8.1           lifecycle_1.0.3        rhdf5filters_1.10.0   
## [115] microViz_0.9.7         jquerylib_0.1.4        flextable_0.8.3       
## [118] cowplot_1.1.1          data.table_1.14.6      bitops_1.0-7          
## [121] R6_2.5.1               stable_1.1.6           IRanges_2.32.0        
## [124] codetools_0.2-18       MASS_7.3-58.1          assertthat_0.2.1      
## [127] xlsxjars_0.6.1         rhdf5_2.42.0           rprojroot_2.0.3       
## [130] withr_2.5.0            S4Vectors_0.36.0       GenomeInfoDbData_1.2.9
## [133] mgcv_1.8-41            parallel_4.2.2         hms_1.1.2             
## [136] grid_4.2.2             rpart_4.1.19           timeDate_4021.107     
## [139] rmarkdown_2.19         carData_3.0-5          rmutil_1.1.10         
## [142] googledrive_2.0.0      Rtsne_0.16             ggpubr_0.5.0          
## [145] base64enc_0.1-3        Biobase_2.56.0         lubridate_1.9.0
```

