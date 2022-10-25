---
title: "NTP72 - 16S - taxa compostiion - human draft "
author: "Florentin Constancias"
date: "October 15, 2022"
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
out_pptx = "~/Desktop/16S-human.pptx"
out_xlsx = "~/Desktop/16S-human.xlsx"
```


```r
load(url("https://github.com/fconstancias/NRP72-FBT/blob/master/Figures/Rastetics.Rdata?raw=true"))
```

## load data:


```r
"data/processed/16S/16S_working_phyloseq.RDS" %>% 
  here::here() %>% 
  readRDS() -> ps_filtered
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
## 139OTUs were removed because they are no longer 
## present in any sample after random subsampling
```

```
## ...
```

```r
ps_filtered %>%
  prune_samples(sample_sums(.)>= min_sample_size, .) -> ps_fil
```

Objectives:

-   compare donor microbiota and pre-treatment / CR - heatmap: v
-   evaluate effect treatment: -- alpha (shanon, simpson --\> hill d1, d2) -- beta (achinson, log -\> bray)
-   describe effect on taxa (Family, Genus, ASV): v -- heatmap v -- alluvial plot v

--\> Ezequiel: process - co-selection of AMR.

## Before Treatment:

Obj: compare donor microbiota and pre-treatment / CR - heatmap


```r
ps_rare %>%
  subset_samples(Day_of_Treatment > -4 & Day_of_Treatment < 0 | Reactor %in% c("DONOR", "CR")) -> before_treat # %>% 
# subset_samples(Reactor %!in%c("IR")) %>% 
# subset_samples() -> before_treat
```


```r
# before_treat %>% 
#   subset_samples(Model2 == "Human1") -> before_treat_h1
```


```r
before_treat %>% 
  # microViz::ps_mutate(tp_var = )
  physeq_most_abundant(group_var = "Model2",
                       ntax = 24,
                       tax_level = "Strain") -> top_taxa_per_group

# before_treat %>% 
#   physeq_most_abundant(group_var = "Reactor_Treatment_Dose",
#                        ntax = 7,
#                        tax_level = "Strain") -> top_taxa_per_group2
```


```r
before_treat %>%
  rarefy_even_depth(sample.size = min_sample_size, rngseed = 123) %>%
  # subset_samples(Reactor == "DONOR") %>%  sample_data() %>%  data.frame()
  # speedyseq::tax_glom(taxrank = "Family") %>% 
  transform_sample_counts(function(x) x/sum(x) * 100) %>%
  subset_taxa(Strain %in% top_taxa_per_group) %>% # extract only the taxa to display - after percentage normalisation
  phyloseq_ampvis_heatmap(physeq = ., tax_aggregate = "Species",
                          transform = FALSE,
                          plot_values = FALSE,
                          facet_by = c("Model2", "Fermentation","Reactor_Treatment_Dose"),
                          group_by = "Day_of_Treatment",#"sample_name",
                          ntax = Inf
  ) -> pre_speces
```

```
## Loading required package: ampvis2
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
## 126OTUs were removed because they are no longer 
## present in any sample after random subsampling
```

```
## ...
```

```
## Warning: There are only 40 taxa, showing all
```

```r
pre_speces$data %>% 
  mutate(Abundance = na_if(Abundance, 0)) %>% 
  mutate(Fermentation = Fermentation %>%  as.character() %>%  replace_na(., "NA")) %>% 
  mutate(Fermentation = factor(Fermentation, levels = c("NA", 1, 2))) -> pre_speces$data 

pre_speces + 
  facet_grid(. ~  Model2 + Fermentation + Reactor_Treatment_Dose, scales = "free", space = "free") +
  scale_color_manual(values = c("#31a354", "orange", 
                                "#f03b20"), drop = FALSE, na.value = 'transparent') -> pre_speces
# scale_fill_viridis_c(breaks = c(0,  0.01, 1, 10, 25, 50, 75, 100),
# labels = c(0,  0.01, 1, 10, 25,  50, 75, 100),
# trans = scales::pseudo_log_trans(sigma = 0.9),
# na.value = 'transparent') -> pre_speces


#pre_speces
```


```r
before_treat %>%
  rarefy_even_depth(sample.size = min_sample_size, rngseed = 123) %>%
  # subset_samples(Reactor == "DONOR") %>%  sample_data() %>%  data.frame()
  # speedyseq::tax_glom(taxrank = "Family") %>% 
  transform_sample_counts(function(x) x/sum(x) * 100) %>%
  phyloseq_ampvis_heatmap(physeq = ., tax_aggregate = "Family",
                          transform = FALSE,
                          plot_values = FALSE,
                          facet_by = c("Model2", "Fermentation","Reactor_Treatment_Dose"),
                          group_by = "Day_of_Treatment",#"sample_name",
                          ntax = 10
  ) -> pre_Family
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
## 126OTUs were removed because they are no longer 
## present in any sample after random subsampling
```

```
## ...
```

```r
pre_Family$data %>% 
  mutate(Abundance = na_if(Abundance, 0)) %>% 
  mutate(Fermentation = Fermentation %>%  as.character() %>%  replace_na(., "NA")) %>% 
  mutate(Fermentation = factor(Fermentation, levels = c("NA", 1, 2))) -> pre_Family$data 


pre_Family + 
  facet_grid(. ~  Model2 + Fermentation + Reactor_Treatment_Dose, scales = "free", space = "free") +
  scale_color_manual(values = c("#31a354", "orange", 
                                "#f03b20"), drop = FALSE, na.value = 'transparent') -> pre_Family


# pre_Family
```


```r
# cowplot::plot_grid(pre_speces,
#           pre_Family, 
#            ncol = 1,
#           align = "v", axis = "b",
#           rel_heights = c(3, 1)) -> pre_plots
# 
# pre_plots
```


```r
ggpubr::ggarrange(pre_Family + theme(legend.position = "null") + 
                    theme(axis.title.x=element_blank(),
                          axis.text.x=element_blank(),
                          axis.ticks.x=element_blank()),
                  pre_speces + theme(
                    strip.background = element_blank(),
                    strip.text.x = element_blank()
                  )  + theme(legend.position = "null"),
                  align = "v",
                  ncol = 1, 
                  heights = c(1, 2.40),
                  common.legend = FALSE) -> pre_plots

pre_plots
```

![](NRP72_draft_16S_Figs_files/figure-html/unnamed-chunk-12-1.png)<!-- -->

```r
pre_plots %>% 
  export::graph2ppt(append = TRUE, 
                    width = 317.48031496 * 6, 
                    height = 0.618 * 317.48031496 * 8 , paper = "A4",  scaling = 2,
                    file = out_pptx)
```

```
## Exported graph as ~/Desktop/16S-human.pptx
```

```r
pre_Family %>% 
  ggpubr::get_legend() %>% 
  ggpubr::as_ggplot() %>% 
  export::graph2ppt(append = TRUE, 
                    width = 317.48031496 * 1, 
                    height = 0.618 * 317.48031496 * 1 , paper = "A4",  scaling = 2,
                    file = out_pptx)
```

```
## Exported graph as ~/Desktop/16S-human.pptx
```

```r
pre_Family$data %>% 
  select(Display, Sample, Abundance,Reactor_Treatment_Dose) %>% 
  pivot_wider(names_from = Display, values_from = Abundance) %>% 
  xlsx::write.xlsx(file = out_xlsx,
                   append = TRUE,
                   sheetName = "pre_family",
                   showNA = TRUE)
pre_speces$data %>% 
  select(Display, Sample, Abundance,Reactor_Treatment_Dose) %>% 
  pivot_wider(names_from = Display, values_from = Abundance) %>% 
  xlsx::write.xlsx(file = out_xlsx,
                   append = TRUE,
                   sheetName = "pre_ASV",
                   showNA = TRUE)
```

## Taxa over time:

### Heatmap:

#### Human1:


```r
ps_rare %>% 
  subset_samples(Model2 == "Human1" & Treatment != "DONOR" & Day_of_Treatment >= 1) %>% 
  physeq_most_abundant(group_var = "Reactor_Treatment_Dose",
                       ntax = 12,
                       tax_level = "Strain") -> top_taxa_per_group_h1
```


```r
ps_rare %>%
  subset_samples(Model2 == "Human1" & Treatment != "DONOR" & Day_of_Treatment >= -1) %>% 
  transform_sample_counts(function(x) x/sum(x) * 100) %>%
  subset_taxa(Strain %in% top_taxa_per_group_h1) %>% # extract only the taxa to display - after percentage normalisation
  phyloseq_ampvis_heatmap(physeq = ., tax_aggregate = "Species",
                          transform = FALSE,
                          plot_values = FALSE,
                          facet_by = c("Model2","Reactor_Treatment_Dose"),
                          group_by = "Day_of_Treatment",#"sample_name",
                          ntax = Inf
  ) -> treat_speces
```

```
## Warning: There are only 40 taxa, showing all
```

```r
treat_speces$data %>% 
  mutate(Abundance = na_if(Abundance, 0)) -> treat_speces$data 

treat_speces + 
  facet_grid(Model2 ~  Reactor_Treatment_Dose , scales = "free", space = "free") +
  scale_color_manual(values = c("#31a354", "orange", 
                                "#f03b20"), drop = FALSE, na.value = 'transparent') -> treat_speces
# scale_fill_viridis_c(breaks = c(0,  0.01, 1, 10, 25, 50, 75, 100),
# labels = c(0,  0.01, 1, 10, 25,  50, 75, 100),
# trans = scales::pseudo_log_trans(sigma = 0.9),
# na.value = 'transparent') -> pre_speces


treat_speces
```

![](NRP72_draft_16S_Figs_files/figure-html/unnamed-chunk-14-1.png)<!-- -->

```r
# pre_speces %>%
#   export::graph2ppt(append = TRUE, width = 317.48031496 * 4  , height = 0.618 * 317.48031496 * 4 , paper = "A4",  scaling = 2,
#                     file = out_pptx)
```


```r
ps_rare %>%
  subset_samples(Model2 == "Human1" & Treatment != "DONOR" & Day_of_Treatment >= -1) %>% 
  transform_sample_counts(function(x) x/sum(x) * 100) %>%
  phyloseq_ampvis_heatmap(physeq = ., tax_aggregate = "Family",
                          transform = FALSE,
                          plot_values = FALSE,
                          facet_by = c("Model2","Reactor_Treatment_Dose"),
                          group_by = "Day_of_Treatment",#"sample_name",
                          ntax = 10
  ) -> treat_fam

treat_fam$data %>% 
  mutate(Abundance = na_if(Abundance, 0)) -> treat_fam$data 

treat_fam + 
  facet_grid(Model2 ~  Reactor_Treatment_Dose , scales = "free", space = "free") +
  scale_color_manual(values = c("#31a354", "orange", 
                                "#f03b20"), drop = FALSE, na.value = 'transparent') -> treat_fam
# scale_fill_viridis_c(breaks = c(0,  0.01, 1, 10, 25, 50, 75, 100),
# labels = c(0,  0.01, 1, 10, 25,  50, 75, 100),
# trans = scales::pseudo_log_trans(sigma = 0.9),
# na.value = 'transparent') -> pre_speces


treat_fam
```

![](NRP72_draft_16S_Figs_files/figure-html/unnamed-chunk-15-1.png)<!-- -->

```r
# pre_speces %>%
#   export::graph2ppt(append = TRUE, width = 317.48031496 * 4  , height = 0.618 * 317.48031496 * 4 , paper = "A4",  scaling = 2,
#                     file = out_pptx)
```


```r
ggpubr::ggarrange(treat_fam + theme(legend.position = "null") + 
                    theme(axis.title.x=element_blank(),
                          axis.text.x=element_blank(),
                          axis.ticks.x=element_blank()),
                  treat_speces + theme(
                    strip.background = element_blank(),
                    strip.text.x = element_blank()
                  )  + theme(legend.position = "null"),
                  align = "v",
                  ncol = 1, 
                  heights = c(1, 2.60),
                  common.legend = FALSE) -> human1_heat

human1_heat
```

![](NRP72_draft_16S_Figs_files/figure-html/unnamed-chunk-16-1.png)<!-- -->

```r
human1_heat %>% 
  export::graph2ppt(append = TRUE, 
                    width = 317.48031496 * 4.75, 
                    height = 0.618 * 317.48031496 * 6.5 , paper = "A4",  scaling = 2,
                    file = out_pptx)
```

```
## Exported graph as ~/Desktop/16S-human.pptx
```

```r
treat_fam$data %>% 
  select(Display, Sample, Abundance,Reactor_Treatment_Dose) %>% 
  pivot_wider(names_from = Display, values_from = Abundance) %>% 
  xlsx::write.xlsx(file = out_xlsx,
                   append = TRUE,
                   sheetName = "treat_family",
                   showNA = TRUE)
treat_speces$data %>% 
  select(Display, Sample, Abundance,Reactor_Treatment_Dose) %>% 
  pivot_wider(names_from = Display, values_from = Abundance) %>% 
  xlsx::write.xlsx(file = out_xlsx,
                   append = TRUE,
                   sheetName = "treat_ASV",
                   showNA = TRUE)
```

##### splited by Fermentations:

###### fermentation 1:


```r
ps_rare %>% 
  subset_samples(Model2 == "Human1" & Treatment != "DONOR" & Day_of_Treatment >= 1 & Fermentation == 1) %>% 
  physeq_most_abundant(group_var = "Reactor_Treatment_Dose",
                       ntax = 12,
                       tax_level = "Strain") -> top_taxa_per_group_h1_f1
```


```r
ps_rare %>%
  subset_samples(Model2 == "Human1" & Treatment != "DONOR" & Day_of_Treatment >= -1 & Fermentation == 1) %>% 
  transform_sample_counts(function(x) x/sum(x) * 100) %>%
  subset_taxa(Strain %in% top_taxa_per_group_h1_f1) %>% # extract only the taxa to display - after percentage normalisation
  phyloseq_ampvis_heatmap(physeq = ., tax_aggregate = "Species",
                          transform = FALSE,
                          plot_values = FALSE,
                          facet_by = c("Model2","Reactor_Treatment_Dose"),
                          group_by = "Day_of_Treatment",#"sample_name",
                          ntax = Inf
  ) -> treat_speces
```

```
## Warning: There are only 35 taxa, showing all
```

```r
treat_speces$data %>% 
  mutate(Abundance = na_if(Abundance, 0)) -> treat_speces$data 

treat_speces + 
  facet_grid(Model2 ~  Reactor_Treatment_Dose , scales = "free", space = "free") +
  scale_color_manual(values = c("#31a354", "orange", 
                                "#f03b20"), drop = FALSE, na.value = 'transparent') -> treat_speces
# scale_fill_viridis_c(breaks = c(0,  0.01, 1, 10, 25, 50, 75, 100),
# labels = c(0,  0.01, 1, 10, 25,  50, 75, 100),
# trans = scales::pseudo_log_trans(sigma = 0.9),
# na.value = 'transparent') -> pre_speces


treat_speces
```

![](NRP72_draft_16S_Figs_files/figure-html/unnamed-chunk-18-1.png)<!-- -->

```r
# pre_speces %>%
#   export::graph2ppt(append = TRUE, width = 317.48031496 * 4  , height = 0.618 * 317.48031496 * 4 , paper = "A4",  scaling = 2,
#                     file = out_pptx)
```


```r
ps_rare %>%
  subset_samples(Model2 == "Human1" & Treatment != "DONOR" & Day_of_Treatment >= -1  & Fermentation == 1) %>% 
  transform_sample_counts(function(x) x/sum(x) * 100) %>%
  phyloseq_ampvis_heatmap(physeq = ., tax_aggregate = "Family",
                          transform = FALSE,
                          plot_values = FALSE,
                          facet_by = c("Model2","Reactor_Treatment_Dose"),
                          group_by = "Day_of_Treatment",#"sample_name",
                          ntax = 10
  ) -> treat_fam

treat_fam$data %>% 
  mutate(Abundance = na_if(Abundance, 0)) -> treat_fam$data 

treat_fam + 
  facet_grid(Model2 ~  Reactor_Treatment_Dose , scales = "free", space = "free") +
  scale_color_manual(values = c("#31a354", "orange", 
                                "#f03b20"), drop = FALSE, na.value = 'transparent') -> treat_fam
# scale_fill_viridis_c(breaks = c(0,  0.01, 1, 10, 25, 50, 75, 100),
# labels = c(0,  0.01, 1, 10, 25,  50, 75, 100),
# trans = scales::pseudo_log_trans(sigma = 0.9),
# na.value = 'transparent') -> pre_speces


treat_fam
```

![](NRP72_draft_16S_Figs_files/figure-html/unnamed-chunk-19-1.png)<!-- -->

```r
# pre_speces %>%
#   export::graph2ppt(append = TRUE, width = 317.48031496 * 4  , height = 0.618 * 317.48031496 * 4 , paper = "A4",  scaling = 2,
#                     file = out_pptx)
```


```r
ggpubr::ggarrange(treat_fam + theme(legend.position = "null") + 
                    theme(axis.title.x=element_blank(),
                          axis.text.x=element_blank(),
                          axis.ticks.x=element_blank()),
                  treat_speces + theme(
                    strip.background = element_blank(),
                    strip.text.x = element_blank()
                  )  + theme(legend.position = "null"),
                  align = "v",
                  ncol = 1, 
                  heights = c(1, 2.60),
                  common.legend = FALSE) -> human1_f1_heat

human1_f1_heat
```

![](NRP72_draft_16S_Figs_files/figure-html/unnamed-chunk-20-1.png)<!-- -->

```r
human1_f1_heat %>% 
  export::graph2ppt(append = TRUE, 
                    width = 317.48031496 * 4.75, 
                    height = 0.618 * 317.48031496 * 6.5 , paper = "A4",  scaling = 2,
                    file = out_pptx)
```

```
## Exported graph as ~/Desktop/16S-human.pptx
```

```r
treat_fam$data %>% 
  select(Display, Sample, Abundance,Reactor_Treatment_Dose) %>% 
  pivot_wider(names_from = Display, values_from = Abundance) %>% 
  xlsx::write.xlsx(file = out_xlsx,
                   append = TRUE,
                   sheetName = "treat_family_h1_f1",
                   showNA = TRUE)
treat_speces$data %>% 
  select(Display, Sample, Abundance,Reactor_Treatment_Dose) %>% 
  pivot_wider(names_from = Display, values_from = Abundance) %>% 
  xlsx::write.xlsx(file = out_xlsx,
                   append = TRUE,
                   sheetName = "treat_ASV_h1_f1",
                   showNA = TRUE)
```

###### fermentation 2:


```r
ps_rare %>% 
  subset_samples(Model2 == "Human1" & Treatment != "DONOR" & Day_of_Treatment >= 1 & Fermentation == 2) %>% 
  physeq_most_abundant(group_var = "Reactor_Treatment_Dose",
                       ntax = 20,
                       tax_level = "Strain") -> top_taxa_per_group_h1_f2
```


```r
ps_rare %>%
  subset_samples(Model2 == "Human1" & Treatment != "DONOR" & Day_of_Treatment >= -1 & Fermentation == 2) %>% 
  transform_sample_counts(function(x) x/sum(x) * 100) %>%
  subset_taxa(Strain %in% top_taxa_per_group_h1_f2) %>% # extract only the taxa to display - after percentage normalisation
  phyloseq_ampvis_heatmap(physeq = ., tax_aggregate = "Species",
                          transform = FALSE,
                          plot_values = FALSE,
                          facet_by = c("Model2","Reactor_Treatment_Dose"),
                          group_by = "Day_of_Treatment",#"sample_name",
                          ntax = Inf
  ) -> treat_speces
```

```
## Warning: There are only 34 taxa, showing all
```

```r
treat_speces$data %>% 
  mutate(Abundance = na_if(Abundance, 0)) -> treat_speces$data 

treat_speces + 
  facet_grid(Model2 ~  Reactor_Treatment_Dose , scales = "free", space = "free") +
  scale_color_manual(values = c("#31a354", "orange", 
                                "#f03b20"), drop = FALSE, na.value = 'transparent') -> treat_speces
# scale_fill_viridis_c(breaks = c(0,  0.01, 1, 10, 25, 50, 75, 100),
# labels = c(0,  0.01, 1, 10, 25,  50, 75, 100),
# trans = scales::pseudo_log_trans(sigma = 0.9),
# na.value = 'transparent') -> pre_speces


treat_speces
```

![](NRP72_draft_16S_Figs_files/figure-html/unnamed-chunk-22-1.png)<!-- -->

```r
# pre_speces %>%
#   export::graph2ppt(append = TRUE, width = 317.48031496 * 4  , height = 0.618 * 317.48031496 * 4 , paper = "A4",  scaling = 2,
#                     file = out_pptx)
```


```r
ps_rare %>%
  subset_samples(Model2 == "Human1" & Treatment != "DONOR" & Day_of_Treatment >= -1  & Fermentation == 2) %>% 
  transform_sample_counts(function(x) x/sum(x) * 100) %>%
  phyloseq_ampvis_heatmap(physeq = ., tax_aggregate = "Family",
                          transform = FALSE,
                          plot_values = FALSE,
                          facet_by = c("Model2","Reactor_Treatment_Dose"),
                          group_by = "Day_of_Treatment",#"sample_name",
                          ntax = 10
  ) -> treat_fam

treat_fam$data %>% 
  mutate(Abundance = na_if(Abundance, 0)) -> treat_fam$data 

treat_fam + 
  facet_grid(Model2 ~  Reactor_Treatment_Dose , scales = "free", space = "free") +
  scale_color_manual(values = c("#31a354", "orange", 
                                "#f03b20"), drop = FALSE, na.value = 'transparent') -> treat_fam
# scale_fill_viridis_c(breaks = c(0,  0.01, 1, 10, 25, 50, 75, 100),
# labels = c(0,  0.01, 1, 10, 25,  50, 75, 100),
# trans = scales::pseudo_log_trans(sigma = 0.9),
# na.value = 'transparent') -> pre_speces


treat_fam
```

![](NRP72_draft_16S_Figs_files/figure-html/unnamed-chunk-23-1.png)<!-- -->

```r
# pre_speces %>%
#   export::graph2ppt(append = TRUE, width = 317.48031496 * 4  , height = 0.618 * 317.48031496 * 4 , paper = "A4",  scaling = 2,
#                     file = out_pptx)
```


```r
ggpubr::ggarrange(treat_fam + theme(legend.position = "null") + 
                    theme(axis.title.x=element_blank(),
                          axis.text.x=element_blank(),
                          axis.ticks.x=element_blank()),
                  treat_speces + theme(
                    strip.background = element_blank(),
                    strip.text.x = element_blank()
                  )  + theme(legend.position = "null"),
                  align = "v",
                  ncol = 1, 
                  heights = c(1, 2.60),
                  common.legend = FALSE) -> human1_f2_heat

human1_f2_heat
```

![](NRP72_draft_16S_Figs_files/figure-html/unnamed-chunk-24-1.png)<!-- -->

```r
human1_f2_heat %>% 
  export::graph2ppt(append = TRUE, 
                    width = 317.48031496 * 4.75, 
                    height = 0.618 * 317.48031496 * 6.5 , paper = "A4",  scaling = 2,
                    file = out_pptx)
```

```
## Exported graph as ~/Desktop/16S-human.pptx
```

```r
treat_fam$data %>% 
  select(Display, Sample, Abundance,Reactor_Treatment_Dose) %>% 
  pivot_wider(names_from = Display, values_from = Abundance) %>% 
  xlsx::write.xlsx(file = out_xlsx,
                   append = TRUE,
                   sheetName = "treat_family_h1_f2",
                   showNA = TRUE)
treat_speces$data %>% 
  select(Display, Sample, Abundance,Reactor_Treatment_Dose) %>% 
  pivot_wider(names_from = Display, values_from = Abundance) %>% 
  xlsx::write.xlsx(file = out_xlsx,
                   append = TRUE,
                   sheetName = "treat_ASV_h1_f2",
                   showNA = TRUE)
```

#### Human2:


```r
ps_rare %>% 
  subset_samples(Model2 == "Human2" & Treatment != "DONOR" & Day_of_Treatment >= 1) %>% 
  physeq_most_abundant(group_var = "Reactor_Treatment_Dose",
                       ntax = 14,
                       tax_level = "Strain") -> top_taxa_per_group_h2
```


```r
ps_rare %>%
  subset_samples(Model2 == "Human2" & Treatment != "DONOR" & Day_of_Treatment >= -1) %>% 
  transform_sample_counts(function(x) x/sum(x) * 100) %>%
  subset_taxa(Strain %in% top_taxa_per_group_h2) %>% # extract only the taxa to display - after percentage normalisation
  phyloseq_ampvis_heatmap(physeq = ., tax_aggregate = "Species",
                          transform = FALSE,
                          plot_values = FALSE,
                          facet_by = c("Model2","Reactor_Treatment_Dose"),
                          group_by = "Day_of_Treatment",#"sample_name",
                          ntax = Inf
  ) -> treat_speces
```

```
## Warning: There are only 35 taxa, showing all
```

```r
treat_speces$data %>% 
  mutate(Abundance = na_if(Abundance, 0)) -> treat_speces$data 

treat_speces + 
  facet_grid(Model2 ~  Reactor_Treatment_Dose , scales = "free", space = "free") +
  scale_color_manual(values = c("#31a354", "orange", 
                                "#f03b20"), drop = FALSE, na.value = 'transparent') -> treat_speces
# scale_fill_viridis_c(breaks = c(0,  0.01, 1, 10, 25, 50, 75, 100),
# labels = c(0,  0.01, 1, 10, 25,  50, 75, 100),
# trans = scales::pseudo_log_trans(sigma = 0.9),
# na.value = 'transparent') -> pre_speces

treat_speces
```

![](NRP72_draft_16S_Figs_files/figure-html/unnamed-chunk-26-1.png)<!-- -->

```r
# pre_speces %>%
#   export::graph2ppt(append = TRUE, width = 317.48031496 * 4  , height = 0.618 * 317.48031496 * 4 , paper = "A4",  scaling = 2,
#                     file = out_pptx)
```


```r
ps_rare %>%
  subset_samples(Model2 == "Human2" & Treatment != "DONOR" & Day_of_Treatment >= -1) %>% 
  transform_sample_counts(function(x) x/sum(x) * 100) %>%
  phyloseq_ampvis_heatmap(physeq = ., tax_aggregate = "Family",
                          transform = FALSE,
                          plot_values = FALSE,
                          facet_by = c("Model2","Reactor_Treatment_Dose"),
                          group_by = "Day_of_Treatment",#"sample_name",
                          ntax = 10
  ) -> treat_fam

treat_fam$data %>% 
  mutate(Abundance = na_if(Abundance, 0)) -> treat_fam$data 

treat_fam + 
  facet_grid(Model2 ~  Reactor_Treatment_Dose , scales = "free", space = "free") +
  scale_color_manual(values = c("#31a354", "orange", 
                                "#f03b20"), drop = FALSE, na.value = 'transparent') -> treat_fam
# scale_fill_viridis_c(breaks = c(0,  0.01, 1, 10, 25, 50, 75, 100),
# labels = c(0,  0.01, 1, 10, 25,  50, 75, 100),
# trans = scales::pseudo_log_trans(sigma = 0.9),
# na.value = 'transparent') -> pre_speces

treat_fam
```

![](NRP72_draft_16S_Figs_files/figure-html/unnamed-chunk-27-1.png)<!-- -->

```r
# pre_speces %>%
#   export::graph2ppt(append = TRUE, width = 317.48031496 * 4  , height = 0.618 * 317.48031496 * 4 , paper = "A4",  scaling = 2,
#                     file = out_pptx)
```


```r
ggpubr::ggarrange(treat_fam + theme(legend.position = "null") + 
                    theme(axis.title.x=element_blank(),
                          axis.text.x=element_blank(),
                          axis.ticks.x=element_blank()),
                  treat_speces + theme(
                    strip.background = element_blank(),
                    strip.text.x = element_blank()
                  )  + theme(legend.position = "null"),
                  align = "v",
                  ncol = 1, 
                  heights = c(1, 2.60),
                  common.legend = FALSE) -> human2_heat


human2_heat
```

![](NRP72_draft_16S_Figs_files/figure-html/unnamed-chunk-28-1.png)<!-- -->

```r
human2_heat %>% 
  export::graph2ppt(append = TRUE, 
                    width = 317.48031496 * 4.75, 
                    height = 0.618 * 317.48031496 * 5.5 , paper = "A4",  scaling = 2,
                    file = out_pptx)
```

```
## Exported graph as ~/Desktop/16S-human.pptx
```

```r
treat_fam$data %>% 
  select(Display, Sample, Abundance,Reactor_Treatment_Dose) %>% 
  pivot_wider(names_from = Display, values_from = Abundance) %>% 
  xlsx::write.xlsx(file = out_xlsx,
                   append = TRUE,
                   sheetName = "treat_family_h2",
                   showNA = TRUE)

treat_speces$data %>% 
  select(Display, Sample, Abundance,Reactor_Treatment_Dose) %>% 
  pivot_wider(names_from = Display, values_from = Abundance) %>% 
  xlsx::write.xlsx(file = out_xlsx,
                   append = TRUE,
                   sheetName = "treat_ASV_h2",
                   showNA = TRUE)
```

Could be improved:

-clustering taxa based on correlation -\> response (pheatmap?) <https://btep.ccr.cancer.gov/docs/data-visualization-with-r/Lesson5_intro_to_ggplot/#make-this-plot-compatible-with-ggplot2-using-the-package-ggplotify>

<https://jcoliver.github.io/learn-r/008-ggplot-dendrograms-and-heatmaps.html>

-additional layer depicting the % of the community covered by the taxa displayed.

### barplot:

#### Human1:


```r
ps_rare %>% 
  subset_samples(Model2 == "Human1" & Treatment != "DONOR" & Day_of_Treatment >= -1) %>% 
  microViz::ps_mutate(Reactor_Treatment_Dose = paste0(Reactor_Treatment_Dose, "_", Fermentation)) %>% 
  transform_sample_counts(function(x) x/sum(x) * 100) %>%
  subset_taxa(Strain %in% top_taxa_per_group_h1) %>% # extract only the taxa to display - after percentage normalisation
  microViz::phyloseq_validate() %>%
  microViz::tax_fix() %>% 
  microViz::ps_arrange(Day_of_Treatment) %>% 
  microViz::comp_barplot(
    tax_transform_for_plot = "identity",
    tax_level = "Strain",
    label = "Day_of_Treatment", # name an alternative variable to label axes
    n_taxa = length(top_taxa_per_group_h1), # give more taxa unique colours
    merge_other = FALSE, # split the "other" category to display alpha diversity
    bar_width = 0.9, # reduce the bar width to 70% of one row
    bar_outline_colour = "grey5",
    sample_order = "default") +
  ylab("Proportion - %") + xlab( " Day (Treatment) ") -> p_hist
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
p_hist + facet_grid(Phylum ~ Reactor_Treatment_Dose, scales = "free", space = "free_x") -> p_hist_human1_strain

p_hist_human1_strain  +  theme_light() + theme(legend.position = "none") -> p_hist_human1_strain_noleg

p_hist_human1_strain_noleg
```

![](NRP72_draft_16S_Figs_files/figure-html/unnamed-chunk-29-1.png)<!-- -->

```r
p_hist_human1_strain %>%  
  ggpubr::get_legend() %>% 
  ggpubr::as_ggplot() -> p_hist_human1_strain_leg
```


```r
ps_rare %>% 
  subset_samples(Model2 == "Human1" & Treatment != "DONOR" & Day_of_Treatment >= -1) %>% 
  microViz::ps_mutate(Reactor_Treatment_Dose = paste0(Reactor_Treatment_Dose, "_", Fermentation)) %>% 
  transform_sample_counts(function(x) x/sum(x) * 100) %>%
  microViz::phyloseq_validate() %>%
  microViz::tax_fix() %>% 
  microViz::ps_arrange(Day_of_Treatment) %>% 
  microViz::comp_barplot(
    tax_transform_for_plot = "identity",
    tax_level = "Family",
    label = "Day_of_Treatment", # name an alternative variable to label axes
    n_taxa = 10, # give more taxa unique colours
    palette = microViz::distinct_palette(10, pal = "kelly"),
    merge_other = FALSE, # split the "other" category to display alpha diversity
    bar_width = 0.9, # reduce the bar width to 70% of one row
    bar_outline_colour = "grey5",
    sample_order = "default") +
  ylab("Proportion - %") + xlab( " Day (Treatment) ") -> p_hist
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
p_hist + facet_grid(Kingdom ~ Reactor_Treatment_Dose, scales = "free", space = "free_x", drop = TRUE) -> p_hist_human1_fam

p_hist_human1_fam  +  theme_light() + theme(legend.position = "none") -> p_hist_human1_fam_no_leg

p_hist_human1_fam_no_leg
```

![](NRP72_draft_16S_Figs_files/figure-html/unnamed-chunk-30-1.png)<!-- -->

```r
p_hist_human1_fam %>%  
  ggpubr::get_legend() %>% 
  ggpubr::as_ggplot() -> p_hist_human1_fam_leg
```

##### splited by Fermentations:

##### Fermentation 1:


```r
ps_rare %>% 
  subset_samples(Model2 == "Human1" & Treatment != "DONOR" & Day_of_Treatment >= -1 & Fermentation == 1) %>% 
  # microViz::ps_mutate(Reactor_Treatment_Dose = paste0(Reactor_Treatment_Dose, "_", Fermentation)) %>% 
  transform_sample_counts(function(x) x/sum(x) * 100) %>%
  subset_taxa(Strain %in% top_taxa_per_group_h1_f1) %>% # extract only the taxa to display - after percentage normalisation
  microViz::phyloseq_validate() %>%
  microViz::tax_fix() %>% 
  microViz::ps_arrange(Day_of_Treatment) %>% 
  microViz::comp_barplot(
    tax_transform_for_plot = "identity",
    tax_level = "Strain",
    label = "Day_of_Treatment", # name an alternative variable to label axes
    n_taxa = length(top_taxa_per_group_h1_f1), # give more taxa unique colours
    merge_other = FALSE, # split the "other" category to display alpha diversity
    bar_width = 0.9, # reduce the bar width to 70% of one row
    bar_outline_colour = "grey5",
    sample_order = "default") +
  ylab("Proportion - %") + xlab( " Day (Treatment) ") -> p_hist
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
p_hist + facet_grid(Phylum ~ Reactor_Treatment_Dose, scales = "free", space = "free_x") -> p_hist_human1_f1_strain

p_hist_human1_f1_strain  +  theme_light() + theme(legend.position = "none") -> p_hist_human1_strain_f1_noleg

p_hist_human1_strain_f1_noleg
```

![](NRP72_draft_16S_Figs_files/figure-html/unnamed-chunk-31-1.png)<!-- -->

```r
p_hist_human1_f1_strain %>%  
  ggpubr::get_legend() %>% 
  ggpubr::as_ggplot() -> p_hist_human1_f1_strain_leg
```


```r
ps_rare %>% 
  subset_samples(Model2 == "Human1" & Treatment != "DONOR" & Day_of_Treatment >= -1 & Fermentation == 1) %>% 
  # microViz::ps_mutate(Reactor_Treatment_Dose = paste0(Reactor_Treatment_Dose, "_", Fermentation)) %>% 
  transform_sample_counts(function(x) x/sum(x) * 100) %>%
  microViz::phyloseq_validate() %>%
  microViz::tax_fix() %>% 
  microViz::ps_arrange(Day_of_Treatment) %>% 
  microViz::comp_barplot(
    tax_transform_for_plot = "identity",
    tax_level = "Family",
    label = "Day_of_Treatment", # name an alternative variable to label axes
    n_taxa = 10, # give more taxa unique colours
    palette = microViz::distinct_palette(10, pal = "kelly"),
    merge_other = FALSE, # split the "other" category to display alpha diversity
    bar_width = 0.9, # reduce the bar width to 70% of one row
    bar_outline_colour = "grey5",
    sample_order = "default") +
  ylab("Proportion - %") + xlab( " Day (Treatment) ") -> p_hist
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
p_hist + facet_grid(Kingdom ~ Reactor_Treatment_Dose, scales = "free", space = "free_x", drop = TRUE) -> p_hist_human1_f1_fam

p_hist_human1_f1_fam  +  theme_light() + theme(legend.position = "none") -> p_hist_human1_f1_fam_no_leg

p_hist_human1_f1_fam_no_leg
```

![](NRP72_draft_16S_Figs_files/figure-html/unnamed-chunk-32-1.png)<!-- -->

```r
p_hist_human1_f1_fam %>%  
  ggpubr::get_legend() %>% 
  ggpubr::as_ggplot() -> p_hist_human1_f1_fam_leg
```

##### Fermentation 2:


```r
ps_rare %>% 
  subset_samples(Model2 == "Human1" & Treatment != "DONOR" & Day_of_Treatment >= -1 & Fermentation == 2) %>% 
  # microViz::ps_mutate(Reactor_Treatment_Dose = paste0(Reactor_Treatment_Dose, "_", Fermentation)) %>% 
  transform_sample_counts(function(x) x/sum(x) * 100) %>%
  subset_taxa(Strain %in% top_taxa_per_group_h1_f1) %>% # extract only the taxa to display - after percentage normalisation
  microViz::phyloseq_validate() %>%
  microViz::tax_fix() %>% 
  microViz::ps_arrange(Day_of_Treatment) %>% 
  microViz::comp_barplot(
    tax_transform_for_plot = "identity",
    tax_level = "Strain",
    label = "Day_of_Treatment", # name an alternative variable to label axes
    n_taxa = length(top_taxa_per_group_h1_f2), # give more taxa unique colours
    merge_other = FALSE, # split the "other" category to display alpha diversity
    bar_width = 0.9, # reduce the bar width to 70% of one row
    bar_outline_colour = "grey5",
    sample_order = "default") +
  ylab("Proportion - %") + xlab( " Day (Treatment) ") -> p_hist
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
p_hist + facet_grid(Phylum ~ Reactor_Treatment_Dose, scales = "free", space = "free_x") -> p_hist_human1_f2_strain

p_hist_human1_f2_strain  +  theme_light() + theme(legend.position = "none") -> p_hist_human1_f2_strain_noleg

p_hist_human1_f2_strain_noleg
```

![](NRP72_draft_16S_Figs_files/figure-html/unnamed-chunk-33-1.png)<!-- -->

```r
p_hist_human1_f2_strain %>%  
  ggpubr::get_legend() %>% 
  ggpubr::as_ggplot() -> p_hist_human1_f2_strain_leg
```


```r
ps_rare %>% 
  subset_samples(Model2 == "Human1" & Treatment != "DONOR" & Day_of_Treatment >= -1 & Fermentation == 2) %>% 
  # microViz::ps_mutate(Reactor_Treatment_Dose = paste0(Reactor_Treatment_Dose, "_", Fermentation)) %>% 
  transform_sample_counts(function(x) x/sum(x) * 100) %>%
  microViz::phyloseq_validate() %>%
  microViz::tax_fix() %>% 
  microViz::ps_arrange(Day_of_Treatment) %>% 
  microViz::comp_barplot(
    tax_transform_for_plot = "identity",
    tax_level = "Family",
    label = "Day_of_Treatment", # name an alternative variable to label axes
    n_taxa = 10, # give more taxa unique colours
    palette = microViz::distinct_palette(10, pal = "kelly"),
    merge_other = FALSE, # split the "other" category to display alpha diversity
    bar_width = 0.9, # reduce the bar width to 70% of one row
    bar_outline_colour = "grey5",
    sample_order = "default") +
  ylab("Proportion - %") + xlab( " Day (Treatment) ") -> p_hist
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
p_hist + facet_grid(Kingdom ~ Reactor_Treatment_Dose, scales = "free", space = "free_x", drop = TRUE) -> p_hist_human1_f2_fam

p_hist_human1_f2_fam  +  theme_light() + theme(legend.position = "none") -> p_hist_human1_f2_fam_no_leg

p_hist_human1_f2_fam_no_leg
```

![](NRP72_draft_16S_Figs_files/figure-html/unnamed-chunk-34-1.png)<!-- -->

```r
p_hist_human1_f2_fam %>%  
  ggpubr::get_legend() %>% 
  ggpubr::as_ggplot() -> p_hist_human1_f2_fam_leg
```

#### Human2:


```r
ps_rare %>% 
  subset_samples(Model2 == "Human2" & Treatment != "DONOR" & Day_of_Treatment >= -1) %>% 
  microViz::ps_mutate(Reactor_Treatment_Dose = paste0(Reactor_Treatment_Dose, "_", Fermentation)) %>% 
  transform_sample_counts(function(x) x/sum(x) * 100) %>%
  subset_taxa(Strain %in% top_taxa_per_group_h2) %>% # extract only the taxa to display - after percentage normalisation
  microViz::phyloseq_validate() %>%
  microViz::tax_fix() %>% 
  microViz::ps_arrange(Day_of_Treatment) %>% 
  microViz::comp_barplot(
    tax_transform_for_plot = "identity",
    tax_level = "Strain",
    label = "Day_of_Treatment", # name an alternative variable to label axes
    n_taxa = length(top_taxa_per_group_h2), # give more taxa unique colours
    merge_other = FALSE, # split the "other" category to display alpha diversity
    bar_width = 0.9, # reduce the bar width to 70% of one row
    bar_outline_colour = "grey5",
    sample_order = "default") +
  ylab("Proportion - %") + xlab( " Day (Treatment) ") -> p_hist
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
p_hist + facet_grid(Phylum ~ Reactor_Treatment_Dose, scales = "free", space = "free_x") -> p_hist_human2_strain

p_hist_human2_strain +  theme_light()  + theme(legend.position = "none") -> p_hist_human2_strain_noleg

p_hist_human2_strain_noleg
```

![](NRP72_draft_16S_Figs_files/figure-html/unnamed-chunk-35-1.png)<!-- -->

```r
p_hist_human2_strain %>%  
  ggpubr::get_legend() %>% 
  ggpubr::as_ggplot() -> p_hist_human2_strain_leg
```


```r
ps_rare %>% 
  subset_samples(Model2 == "Human2" & Treatment != "DONOR" & Day_of_Treatment >= -1) %>% 
  microViz::ps_mutate(Reactor_Treatment_Dose = paste0(Reactor_Treatment_Dose, "_", Fermentation)) %>% 
  transform_sample_counts(function(x) x/sum(x) * 100) %>%
  # subset_taxa(Strain %in% top_taxa_per_group_h2) %>% # extract only the taxa to display - after percentage normalisation
  microViz::phyloseq_validate() %>%
  microViz::tax_fix() %>% 
  microViz::ps_arrange(Day_of_Treatment) %>% 
  microViz::comp_barplot(
    tax_transform_for_plot = "identity",
    tax_level = "Family",
    label = "Day_of_Treatment", # name an alternative variable to label axes
    palette = microViz::distinct_palette(10, pal = "kelly"),
    n_taxa = 10, # give more taxa unique colours
    merge_other = FALSE, # split the "other" category to display alpha diversity
    bar_width = 0.9, # reduce the bar width to 70% of one row
    bar_outline_colour = "grey5",
    sample_order = "default") +
  ylab("Proportion - %") + xlab( " Day (Treatment) ") -> p_hist
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
p_hist + facet_grid(Kingdom ~ Reactor_Treatment_Dose, scales = "free", space = "free_x", drop = TRUE) -> p_hist_human2_fam

p_hist_human2_fam +  theme_light()  + theme(legend.position = "none") -> p_hist_human2_fam_no_leg

p_hist_human2_fam_no_leg
```

![](NRP72_draft_16S_Figs_files/figure-html/unnamed-chunk-36-1.png)<!-- -->

```r
p_hist_human2_fam %>%  
  ggpubr::get_legend() %>% 
  ggpubr::as_ggplot() -> p_hist_human2_fam_leg
```

Unify colors:


```r
c(p_hist_human2_fam$data$top, p_hist_human1_f1_fam$data$top, p_hist_human1_f2_fam$data$top) %>% 
  unique() -> families
```


```r
# families %>%
#   ggpubr::get_palette(k = ., palette = "kelly") -> col

families[!families %in% c("other")] %>% 
  length() %>% 
  microViz::distinct_palette(n = ., pal = "kelly", add = NA) -> col_fam

names(col_fam) <- families[!families %in% c("other")]

c(col_fam, "lightgrey") -> col_fam

names(col_fam) = NULL

names(col_fam) <- c(families[!families %in% c("other")], families[families %in% c("other")])

scales::show_col(col_fam)
```

![](NRP72_draft_16S_Figs_files/figure-html/unnamed-chunk-38-1.png)<!-- -->


```r
p_hist_human1_f1_fam  +  theme_light() + scale_fill_manual(values = col_fam) ->  p_hist_human1_f1_fam
```

```
## Scale for 'fill' is already present. Adding another scale for 'fill', which
## will replace the existing scale.
```

```r
p_hist_human1_f1_fam + theme(legend.position = "none") -> p_hist_human1_f1_fam_no_leg

p_hist_human1_f1_fam_no_leg
```

![](NRP72_draft_16S_Figs_files/figure-html/unnamed-chunk-39-1.png)<!-- -->


```r
p_hist_human1_f1_fam %>%  
  ggpubr::get_legend() %>% 
  ggpubr::as_ggplot() -> p_hist_human1_f1_fam_leg

p_hist_human1_f1_fam_leg
```

![](NRP72_draft_16S_Figs_files/figure-html/unnamed-chunk-40-1.png)<!-- -->


```r
p_hist_human1_f2_fam  +  theme_light() + scale_fill_manual(values = col_fam) ->  p_hist_human1_f2_fam
```

```
## Scale for 'fill' is already present. Adding another scale for 'fill', which
## will replace the existing scale.
```

```r
p_hist_human1_f2_fam + theme(legend.position = "none") -> p_hist_human1_f2_fam_no_leg

p_hist_human1_f2_fam_no_leg
```

![](NRP72_draft_16S_Figs_files/figure-html/unnamed-chunk-41-1.png)<!-- -->


```r
p_hist_human1_f2_fam %>%  
  ggpubr::get_legend() %>% 
  ggpubr::as_ggplot() -> p_hist_human1_f2_fam_leg

p_hist_human1_f2_fam_leg
```

![](NRP72_draft_16S_Figs_files/figure-html/unnamed-chunk-42-1.png)<!-- -->


```r
p_hist_human2_fam +  theme_light() +  scale_fill_manual(values = col_fam) ->  p_hist_human2_fam
```

```
## Scale for 'fill' is already present. Adding another scale for 'fill', which
## will replace the existing scale.
```

```r
p_hist_human2_fam + theme(legend.position = "none") ->p_hist_human2_fam_no_leg

p_hist_human2_fam_no_leg
```

![](NRP72_draft_16S_Figs_files/figure-html/unnamed-chunk-43-1.png)<!-- -->


```r
p_hist_human2_fam %>%  
  ggpubr::get_legend() %>% 
  ggpubr::as_ggplot() -> p_hist_human2_fam_leg

p_hist_human2_fam_leg
```

![](NRP72_draft_16S_Figs_files/figure-html/unnamed-chunk-44-1.png)<!-- -->


```r
c(p_hist_human2_strain$data$top, p_hist_human1_f1_strain$data$top, p_hist_human1_f2_strain$data$top) %>% 
  unique() -> ASVs

# ASVs
```


```r
# families %>%
#   ggpubr::get_palette(k = ., palette = "kelly") -> col
set.seed(1234)

ASVs[!ASVs %in% c("other")] %>% 
  length() %>% 
  randomcoloR::distinctColorPalette(k = ., altCol = FALSE, runTsne = TRUE) -> col_ASVs

names(col_ASVs) <- ASVs[!ASVs %in% c("other")]

c(col_ASVs, "lightgrey") -> col_ASVs

names(col_ASVs) = NULL

names(col_ASVs) <- c(ASVs[!ASVs %in% c("other")], ASVs[ASVs %in% c("other")])

scales::show_col(col_ASVs)
```

![](NRP72_draft_16S_Figs_files/figure-html/unnamed-chunk-46-1.png)<!-- -->


```r
p_hist_human1_f1_strain + scale_fill_manual(values = col_ASVs) ->  p_hist_human1_f1_strain
```

```
## Scale for 'fill' is already present. Adding another scale for 'fill', which
## will replace the existing scale.
```

```r
p_hist_human1_f1_strain + theme_light() +
  theme(legend.position = "none") -> p_hist_human1_f1_strain_no_leg

p_hist_human1_f1_strain_no_leg
```

![](NRP72_draft_16S_Figs_files/figure-html/unnamed-chunk-47-1.png)<!-- -->


```r
p_hist_human1_f1_strain %>%  
  ggpubr::get_legend() %>% 
  ggpubr::as_ggplot() -> p_hist_human1_f1_strain_leg

p_hist_human1_f1_strain_leg
```

![](NRP72_draft_16S_Figs_files/figure-html/unnamed-chunk-48-1.png)<!-- -->


```r
p_hist_human1_f2_strain + scale_fill_manual(values = col_ASVs) ->  p_hist_human1_f2_strain
```

```
## Scale for 'fill' is already present. Adding another scale for 'fill', which
## will replace the existing scale.
```

```r
p_hist_human1_f2_strain + theme_light() + theme(legend.position = "none") -> p_hist_human1_f2_strain_no_leg

p_hist_human1_f2_strain_no_leg
```

![](NRP72_draft_16S_Figs_files/figure-html/unnamed-chunk-49-1.png)<!-- -->


```r
p_hist_human1_f2_strain %>%  
  ggpubr::get_legend() %>% 
  ggpubr::as_ggplot() -> p_hist_human1_f2_strain_leg

p_hist_human1_f2_strain_leg
```

![](NRP72_draft_16S_Figs_files/figure-html/unnamed-chunk-50-1.png)<!-- -->


```r
p_hist_human2_strain + scale_fill_manual(values = col_ASVs) ->  p_hist_human2_strain
```

```
## Scale for 'fill' is already present. Adding another scale for 'fill', which
## will replace the existing scale.
```

```r
p_hist_human2_strain + theme_light() +theme(legend.position = "none") ->p_hist_human2_strain_no_leg

p_hist_human2_strain_no_leg
```

![](NRP72_draft_16S_Figs_files/figure-html/unnamed-chunk-51-1.png)<!-- -->


```r
p_hist_human2_strain %>%  
  ggpubr::get_legend() %>% 
  ggpubr::as_ggplot() -> p_hist_human2_strain_leg

p_hist_human2_strain_leg
```

![](NRP72_draft_16S_Figs_files/figure-html/unnamed-chunk-52-1.png)<!-- -->


```r
ggpubr::ggarrange(p_hist_human1_f1_fam_no_leg,
                  p_hist_human1_f1_strain_no_leg ,
                  align = "v",
                  ncol = 1, 
                  heights = c(1, 1),
                  common.legend = FALSE) -> p_hist_human1_f1

p_hist_human1_f1
```

![](NRP72_draft_16S_Figs_files/figure-html/unnamed-chunk-53-1.png)<!-- -->

```r
p_hist_human1_f1 %>% 
  export::graph2ppt(append = TRUE, 
                    width = 317.48031496 * 4.75, 
                    height = 0.618 * 317.48031496 * 6.5 , paper = "A4",  scaling = 2,
                    file = out_pptx)
```

```
## Exported graph as ~/Desktop/16S-human.pptx
```

```r
p_hist_human1_f1_fam_leg %>% 
  export::graph2ppt(append = TRUE, 
                    width = 317.48031496 * 1, 
                    height = 0.618 * 317.48031496 * 1 , paper = "A4",  scaling = 2,
                    file = out_pptx)
```

```
## Exported graph as ~/Desktop/16S-human.pptx
```

```r
p_hist_human1_f1_strain_leg %>% 
  export::graph2ppt(append = TRUE, 
                    width = 317.48031496 * 5, 
                    height = 0.618 * 317.48031496 * 2.5 , paper = "A3",  scaling = 2,
                    file = out_pptx)
```

```
## Exported graph as ~/Desktop/16S-human.pptx
```


```r
ggpubr::ggarrange(p_hist_human1_f2_fam_no_leg,
                  p_hist_human1_f2_strain_no_leg,
                  align = "v",
                  ncol = 1, 
                  heights = c(1, 1),
                  common.legend = FALSE) -> p_hist_human1_f2

p_hist_human1_f2
```

![](NRP72_draft_16S_Figs_files/figure-html/unnamed-chunk-54-1.png)<!-- -->

```r
p_hist_human1_f2 %>% 
  export::graph2ppt(append = TRUE, 
                    width = 317.48031496 * 4.75, 
                    height = 0.618 * 317.48031496 * 6.5 , paper = "A4",  scaling = 2,
                    file = out_pptx)
```

```
## Exported graph as ~/Desktop/16S-human.pptx
```

```r
p_hist_human1_f2_fam_leg %>% 
  export::graph2ppt(append = TRUE, 
                    width = 317.48031496 * 1, 
                    height = 0.618 * 317.48031496 * 1 , paper = "A4",  scaling = 2,
                    file = out_pptx)
```

```
## Exported graph as ~/Desktop/16S-human.pptx
```

```r
p_hist_human1_f2_strain_leg %>% 
  export::graph2ppt(append = TRUE, 
                    width = 317.48031496 * 6, 
                    height = 0.618 * 317.48031496 * 2.5 , paper = "A3",  scaling = 2,
                    file = out_pptx)
```

```
## Exported graph as ~/Desktop/16S-human.pptx
```


```r
ggpubr::ggarrange(p_hist_human2_fam_no_leg,
                  p_hist_human2_strain_noleg,
                  align = "v",
                  ncol = 1, 
                  heights = c(1, 1),
                  common.legend = FALSE) -> p_hist_human2

p_hist_human2
```

![](NRP72_draft_16S_Figs_files/figure-html/unnamed-chunk-55-1.png)<!-- -->

```r
p_hist_human2 %>% 
  export::graph2ppt(append = TRUE, 
                    width = 317.48031496 * 4.75, 
                    height = 0.618 * 317.48031496 * 6.5 , paper = "A4",  scaling = 2,
                    file = out_pptx)
```

```
## Exported graph as ~/Desktop/16S-human.pptx
```

```r
p_hist_human2_fam_leg %>% 
  export::graph2ppt(append = TRUE, 
                    width = 317.48031496 * 1, 
                    height = 0.618 * 317.48031496 * 1 , paper = "A4",  scaling = 2,
                    file = out_pptx)
```

```
## Exported graph as ~/Desktop/16S-human.pptx
```

```r
p_hist_human2_strain_leg %>% 
  export::graph2ppt(append = TRUE, 
                    width = 317.48031496 * 6, 
                    height = 0.618 * 317.48031496 * 2.5 , paper = "A3",  scaling = 2,
                    file = out_pptx)
```

```
## Exported graph as ~/Desktop/16S-human.pptx
```

### Alluvial:


```r
# p_hist_human2_fam$data %>%
#   filter(Model2 == "Human2") %>%
#   # microViz::ps_mutate(Reactor_Treatment_Dose = paste0(Reactor_Treatment_Dose, "_", Fermentation)) %>% 
#   # mutate(combo = paste0(Model2, "_" , Reactor_F1_F2)) %>%
#   # mutate(combo = as.factor(combo)) %>%
#   select(SAMPLE, top,OTU, Abundance, Day_of_Treatment, Reactor_Treatment_Dose) %>%  distinct(top)
```


```r
# ggpubr::ggarrange(p_hist_human2_fam_no_leg,
#                   p_hist_human2_strain_noleg ,
#                   align = "v",
#                   ncol = 1, 
#                   heights = c(1, 1),
#                   common.legend = FALSE) -> p_hist_human2
# 
# p_hist_human2
# 
# p_hist_human2 %>% 
#   export::graph2ppt(append = TRUE, 
#                     width = 317.48031496 * 4.75, 
#                     height = 0.618 * 317.48031496 * 6.5 , paper = "A4",  scaling = 2,
#                     file = out_pptx)
```

#### human1 - f1


```r
require(ggalluvial)
```

```
## Loading required package: ggalluvial
```

##### Family


```r
p_hist_human1_f1_fam$data %>%
  # filter(Model2 == "Human2") %>%
  # microViz::ps_mutate(Reactor_Treatment_Dose = paste0(Reactor_Treatment_Dose, "_", Fermentation)) %>% 
  # mutate(combo = paste0(Model2, "_" , Reactor_F1_F2)) %>%
  # mutate(combo = as.factor(combo)) %>%
  select(SAMPLE, top,OTU, Abundance, Day_of_Treatment, Reactor_Treatment_Dose) %>% 
  ggplot(.,
         aes(x = Day_of_Treatment, 
             y = Abundance,
             stratum = top, 
             alluvium = OTU,
             fill = top,
             label = top))+
  facet_grid(Reactor_Treatment_Dose  ~ .)+
  # scale_fill_manual(values = c("red", "yellow", "green3"))+
  geom_flow(stat = "alluvium", lode.guidance = "frontback", color = "darkgray")+
  geom_stratum() +
  theme_light() +
  # geom_text(stat = "alluvium", aes(label = top), lode.guidance = "frontback") +
  scale_fill_manual(values = col_fam) +
  # scale_color_manual(values  = microViz::distinct_palette(10, pal = "kelly")) +
  # scale_fill_manual(values  = microViz::distinct_palette(10, pal = "kelly")) +
  ylab("Proportion") -> p_allu_human1_f1

p_allu_human1_f1 + theme(legend.position = "none") -> p_allu_human1_f1_noleg 

p_allu_human1_f1_noleg
```

![](NRP72_draft_16S_Figs_files/figure-html/unnamed-chunk-59-1.png)<!-- -->

```r
p_allu_human1_f1 %>% 
  ggpubr::get_legend() %>% 
  ggpubr::as_ggplot() ->  p_allu_human1_f1_leg
```


```r
p_allu_human1_f1_noleg %>% 
  export::graph2ppt(append = TRUE, 
                    width = 317.48031496 * 1, 
                    height = 0.618 * 317.48031496 * 2 , paper = "A4",  scaling = 2,
                    file = out_pptx)
```

```
## Exported graph as ~/Desktop/16S-human.pptx
```

```r
p_allu_human1_f1_leg %>% 
  export::graph2ppt(append = TRUE, 
                    width = 317.48031496 * 1, 
                    height = 0.618 * 317.48031496 * 2 , paper = "A4",  scaling = 2,
                    file = out_pptx)
```

```
## Exported graph as ~/Desktop/16S-human.pptx
```

##### ASV:


```r
p_hist_human1_f1_strain$data %>%
  # filter(Model2 == "Human2") %>%
  # microViz::ps_mutate(Reactor_Treatment_Dose = paste0(Reactor_Treatment_Dose, "_", Fermentation)) %>% 
  # mutate(combo = paste0(Model2, "_" , Reactor_F1_F2)) %>%
  # mutate(combo = as.factor(combo)) %>%
  select(SAMPLE, top,OTU, Abundance, Day_of_Treatment, Reactor_Treatment_Dose) %>% 
  ggplot(.,
         aes(x = Day_of_Treatment, 
             y = Abundance,
             stratum = top, 
             alluvium = OTU,
             fill = top,
             label = top))+
  facet_grid(Reactor_Treatment_Dose  ~ .)+
  # scale_fill_manual(values = c("red", "yellow", "green3"))+
  geom_flow(stat = "alluvium", lode.guidance = "frontback", color = "darkgray")+
  geom_stratum() +
  theme_light() +
  # geom_text(stat = "alluvium", aes(label = top), lode.guidance = "frontback") +
  scale_fill_manual(values = col_ASVs) +
  # scale_color_manual(values  = microViz::distinct_palette(10, pal = "kelly")) +
  # scale_fill_manual(values  = microViz::distinct_palette(10, pal = "kelly")) +
  ylab("Proportion") -> p_allu_human1_ASV_f1

p_allu_human1_ASV_f1 + theme(legend.position = "none") -> p_allu_human1_ASV_f1_noleg 

p_allu_human1_ASV_f1_noleg
```

![](NRP72_draft_16S_Figs_files/figure-html/unnamed-chunk-61-1.png)<!-- -->

```r
p_allu_human1_ASV_f1 %>% 
  ggpubr::get_legend() %>% 
  ggpubr::as_ggplot() ->  p_allu_human1_ASV_f1_leg
```


```r
p_allu_human1_ASV_f1_noleg %>% 
  export::graph2ppt(append = TRUE, 
                    width = 317.48031496 * 1, 
                    height = 0.618 * 317.48031496 * 2 , paper = "A4",  scaling = 2,
                    file = out_pptx)
```

```
## Exported graph as ~/Desktop/16S-human.pptx
```

```r
p_allu_human1_ASV_f1_leg %>% 
  export::graph2ppt(append = TRUE, 
                    width = 317.48031496 * 4, 
                    height = 0.618 * 317.48031496 * 2 , paper = "A3",  scaling = 2,
                    file = out_pptx)
```

```
## Exported graph as ~/Desktop/16S-human.pptx
```

#### human1 - f2

##### Family


```r
p_hist_human1_f2_fam$data %>%
  # filter(Model2 == "Human2") %>%
  # microViz::ps_mutate(Reactor_Treatment_Dose = paste0(Reactor_Treatment_Dose, "_", Fermentation)) %>% 
  # mutate(combo = paste0(Model2, "_" , Reactor_F1_F2)) %>%
  # mutate(combo = as.factor(combo)) %>%
  select(SAMPLE, top,OTU, Abundance, Day_of_Treatment, Reactor_Treatment_Dose) %>% 
  ggplot(.,
         aes(x = Day_of_Treatment, 
             y = Abundance,
             stratum = top, 
             alluvium = OTU,
             fill = top,
             label = top))+
  facet_grid(Reactor_Treatment_Dose  ~ .)+
  # scale_fill_manual(values = c("red", "yellow", "green3"))+
  geom_flow(stat = "alluvium", lode.guidance = "frontback", color = "darkgray")+
  geom_stratum() +
  theme_light() +
  # geom_text(stat = "alluvium", aes(label = top), lode.guidance = "frontback") +
  scale_fill_manual(values = col_fam) +
  # scale_color_manual(values  = microViz::distinct_palette(10, pal = "kelly")) +
  # scale_fill_manual(values  = microViz::distinct_palette(10, pal = "kelly")) +
  ylab("Proportion") -> p_allu_human1_f2

p_allu_human1_f2 + theme(legend.position = "none") -> p_allu_human1_f2_noleg 

p_allu_human1_f2_noleg
```

![](NRP72_draft_16S_Figs_files/figure-html/unnamed-chunk-63-1.png)<!-- -->

```r
p_allu_human1_f2 %>% 
  ggpubr::get_legend() %>% 
  ggpubr::as_ggplot() ->  p_allu_human1_f2_leg
```


```r
p_allu_human1_f2_noleg %>% 
  export::graph2ppt(append = TRUE, 
                    width = 317.48031496 * 1, 
                    height = 0.618 * 317.48031496 * 2 , paper = "A4",  scaling = 2,
                    file = out_pptx)
```

```
## Exported graph as ~/Desktop/16S-human.pptx
```

```r
p_allu_human1_f2_leg %>% 
  export::graph2ppt(append = TRUE, 
                    width = 317.48031496 * 1, 
                    height = 0.618 * 317.48031496 * 2 , paper = "A4",  scaling = 2,
                    file = out_pptx)
```

```
## Exported graph as ~/Desktop/16S-human.pptx
```

##### ASV:


```r
p_hist_human1_f2_strain$data %>%
  # filter(Model2 == "Human2") %>%
  # microViz::ps_mutate(Reactor_Treatment_Dose = paste0(Reactor_Treatment_Dose, "_", Fermentation)) %>% 
  # mutate(combo = paste0(Model2, "_" , Reactor_F1_F2)) %>%
  # mutate(combo = as.factor(combo)) %>%
  select(SAMPLE, top,OTU, Abundance, Day_of_Treatment, Reactor_Treatment_Dose) %>% 
  ggplot(.,
         aes(x = Day_of_Treatment, 
             y = Abundance,
             stratum = top, 
             alluvium = OTU,
             fill = top,
             label = top))+
  facet_grid(Reactor_Treatment_Dose  ~ .)+
  # scale_fill_manual(values = c("red", "yellow", "green3"))+
  geom_flow(stat = "alluvium", lode.guidance = "frontback", color = "darkgray")+
  geom_stratum() +
  theme_light() +
  # geom_text(stat = "alluvium", aes(label = top), lode.guidance = "frontback") +
  scale_fill_manual(values = col_ASVs) +
  # scale_color_manual(values  = microViz::distinct_palette(10, pal = "kelly")) +
  # scale_fill_manual(values  = microViz::distinct_palette(10, pal = "kelly")) +
  ylab("Proportion") -> p_allu_human1_ASV_f2

p_allu_human1_ASV_f2 + theme(legend.position = "none") -> p_allu_human1_ASV_f2_noleg 

p_allu_human1_ASV_f2_noleg
```

![](NRP72_draft_16S_Figs_files/figure-html/unnamed-chunk-65-1.png)<!-- -->

```r
p_allu_human1_ASV_f1 %>% 
  ggpubr::get_legend() %>% 
  ggpubr::as_ggplot() ->  p_allu_human1_ASV_f2_leg
```


```r
p_allu_human1_ASV_f2_noleg %>% 
  export::graph2ppt(append = TRUE, 
                    width = 317.48031496 * 1, 
                    height = 0.618 * 317.48031496 * 2 , paper = "A4",  scaling = 2,
                    file = out_pptx)
```

```
## Exported graph as ~/Desktop/16S-human.pptx
```

```r
p_allu_human1_ASV_f2_leg %>% 
  export::graph2ppt(append = TRUE, 
                    width = 317.48031496 * 4, 
                    height = 0.618 * 317.48031496 * 2 , paper = "A3",  scaling = 2,
                    file = out_pptx)
```

```
## Exported graph as ~/Desktop/16S-human.pptx
```

#### human2

##### Family


```r
p_hist_human2_fam$data %>%
  # filter(Model2 == "Human2") %>%
  # microViz::ps_mutate(Reactor_Treatment_Dose = paste0(Reactor_Treatment_Dose, "_", Fermentation)) %>% 
  # mutate(combo = paste0(Model2, "_" , Reactor_F1_F2)) %>%
  # mutate(combo = as.factor(combo)) %>%
  select(SAMPLE, top,OTU, Abundance, Day_of_Treatment, Reactor_Treatment_Dose) %>% 
  ggplot(.,
         aes(x = Day_of_Treatment, 
             y = Abundance,
             stratum = top, 
             alluvium = OTU,
             fill = top,
             label = top))+
  facet_grid(Reactor_Treatment_Dose  ~ .)+
  # scale_fill_manual(values = c("red", "yellow", "green3"))+
  geom_flow(stat = "alluvium", lode.guidance = "frontback", color = "darkgray")+
  geom_stratum() +
  theme_light() +
  # geom_text(stat = "alluvium", aes(label = top), lode.guidance = "frontback") +
  scale_fill_manual(values = col_fam) +
  # scale_color_manual(values  = microViz::distinct_palette(10, pal = "kelly")) +
  # scale_fill_manual(values  = microViz::distinct_palette(10, pal = "kelly")) +
  ylab("Proportion") -> p_allu_human2

p_allu_human2 + theme(legend.position = "none") -> p_allu_human2_noleg 

p_allu_human2_noleg
```

![](NRP72_draft_16S_Figs_files/figure-html/unnamed-chunk-67-1.png)<!-- -->


```r
p_allu_human2 %>% 
  ggpubr::get_legend() %>% 
  ggpubr::as_ggplot() ->  p_allu_human2_leg

p_allu_human2_leg
```

![](NRP72_draft_16S_Figs_files/figure-html/unnamed-chunk-68-1.png)<!-- -->


```r
p_allu_human2_noleg %>% 
  export::graph2ppt(append = TRUE, 
                    width = 317.48031496 * 1, 
                    height = 0.618 * 317.48031496 * 2 , paper = "A4",  scaling = 2,
                    file = out_pptx)
```

```
## Exported graph as ~/Desktop/16S-human.pptx
```

```r
p_allu_human2_leg %>% 
  export::graph2ppt(append = TRUE, 
                    width = 317.48031496 * 1, 
                    height = 0.618 * 317.48031496 * 2 , paper = "A4",  scaling = 2,
                    file = out_pptx)
```

```
## Exported graph as ~/Desktop/16S-human.pptx
```

##### ASV:


```r
p_hist_human2_strain$data %>%
  # filter(Model2 == "Human2") %>%
  # microViz::ps_mutate(Reactor_Treatment_Dose = paste0(Reactor_Treatment_Dose, "_", Fermentation)) %>% 
  # mutate(combo = paste0(Model2, "_" , Reactor_F1_F2)) %>%
  # mutate(combo = as.factor(combo)) %>%
  select(SAMPLE, top,OTU, Abundance, Day_of_Treatment, Reactor_Treatment_Dose) %>% 
  ggplot(.,
         aes(x = Day_of_Treatment, 
             y = Abundance,
             stratum = top, 
             alluvium = OTU,
             fill = top,
             label = top))+
  facet_grid(Reactor_Treatment_Dose  ~ .)+
  # scale_fill_manual(values = c("red", "yellow", "green3"))+
  geom_flow(stat = "alluvium", lode.guidance = "frontback", color = "darkgray")+
  geom_stratum() +
  theme_light() +
  # geom_text(stat = "alluvium", aes(label = top), lode.guidance = "frontback") +
  scale_fill_manual(values = col_ASVs) +
  # scale_color_manual(values  = microViz::distinct_palette(10, pal = "kelly")) +
  # scale_fill_manual(values  = microViz::distinct_palette(10, pal = "kelly")) +
  ylab("Proportion") -> p_allu_human2_ASV

p_allu_human2_ASV + theme(legend.position = "none") -> p_allu_human2_ASV_noleg 

p_allu_human2_ASV_noleg
```

![](NRP72_draft_16S_Figs_files/figure-html/unnamed-chunk-70-1.png)<!-- -->

```r
p_allu_human2_ASV %>% 
  ggpubr::get_legend() %>% 
  ggpubr::as_ggplot() ->  p_allu_human2_ASV_leg
```


```r
p_allu_human2_ASV_noleg %>% 
  export::graph2ppt(append = TRUE, 
                    width = 317.48031496 * 1, 
                    height = 0.618 * 317.48031496 * 2 , paper = "A4",  scaling = 2,
                    file = out_pptx)
```

```
## Exported graph as ~/Desktop/16S-human.pptx
```

```r
p_allu_human2_ASV_leg %>% 
  export::graph2ppt(append = TRUE, 
                    width = 317.48031496 * 4, 
                    height = 0.618 * 317.48031496 * 2 , paper = "A3",  scaling = 2,
                    file = out_pptx)
```

```
## Exported graph as ~/Desktop/16S-human.pptx
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
## 1     0     12     3      3   102     502            163             121       3
## # … with abbreviated variable name ¹​gram_neg_class
```


```r
dim(gram)
```

```
## [1] 589   9
```

```r
intersect(gram$ASV,
          tax_table(ps_filtered)[,"Strain"] %>%  as.character() ) %>%  length() 
```

```
## [1] 589
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

### human2:

#### gram:


```r
ps_rare_gram_int %>% 
  subset_samples(Model2 == "Human2" & Treatment != "DONOR" & Day_of_Treatment >= -1) %>% 
  microViz::ps_mutate(Reactor_Treatment_Dose = paste0(Reactor_Treatment_Dose, "_", Fermentation)) %>% 
  transform_sample_counts(function(x) x/sum(x) * 100) %>%
  physeq_sel_tax_table(c("Kingdom","gram_stain")) %>% 
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
p_gram_h2 + theme_light() + facet_grid(Kingdom ~ Reactor_Treatment_Dose, scales = "free", space = "free_x", drop = TRUE) -> p_gram_h2

p_gram_h2 + theme(legend.position = "none") -> p_gram_h2_no_leg

p_gram_h2_no_leg
```

![](NRP72_draft_16S_Figs_files/figure-html/unnamed-chunk-77-1.png)<!-- -->

```r
p_gram_h2 %>%  
  ggpubr::get_legend() %>% 
  ggpubr::as_ggplot() -> p_gram_h2_leg

p_gram_h2_leg
```

![](NRP72_draft_16S_Figs_files/figure-html/unnamed-chunk-77-2.png)<!-- -->


```r
p_gram_h2$data %>%
  # filter(Model2 == "Human2") %>%
  mutate(top = factor(top, levels = c("Gram-positive", "Gram-negative", "unknown-gram"))) %>%
  # mutate(combo = paste0(Model2, "_" , Reactor_F1_F2)) %>%
  # mutate(combo = as.factor(combo)) %>%
  select(SAMPLE, top,OTU, Abundance, Day_of_Treatment, Reactor_Treatment_Dose) %>% 
  ggplot(.,
         aes(x = Day_of_Treatment, 
             y = Abundance,
             stratum = top, 
             alluvium = OTU,
             fill = top,
             label = top))+
  facet_grid(Reactor_Treatment_Dose  ~ .)+
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

![](NRP72_draft_16S_Figs_files/figure-html/unnamed-chunk-78-1.png)<!-- -->

```r
p_allu_human2_gram %>% 
  ggpubr::get_legend() %>% 
  ggpubr::as_ggplot() ->  p_allu_human2_gram_leg

p_allu_human2_gram_leg
```

![](NRP72_draft_16S_Figs_files/figure-html/unnamed-chunk-78-2.png)<!-- -->

#### intrinsic:


```r
ps_rare_gram_int %>% 
  subset_samples(Model2 == "Human2" & Treatment != "DONOR" & Day_of_Treatment >= -1) %>% 
  microViz::ps_mutate(Reactor_Treatment_Dose = paste0(Reactor_Treatment_Dose, "_", Fermentation)) %>% 
  transform_sample_counts(function(x) x/sum(x) * 100) %>%
  physeq_sel_tax_table(c("Kingdom", "intrinsic")) %>% 
  # microbiome::transform("log") %>% 
  # speedyseq::tax_glom("intrinsic") %>% 
  # physeq_sel_tax_table(c("gram_stain","intrinsic")) %>% 
  # physeq_glom_rename(phyloseq = ., speedyseq = TRUE, rename_ASV = FALSE, taxrank = "gram_stain") %>% 
  # subset_taxa(Strain %in% top_taxa_per_group_h2) %>% # extract only the taxa to display - after percentage normalisation
  microViz::phyloseq_validate() %>%
  microViz::tax_fix() %>%
  microViz::ps_arrange(Day_of_Treatment) %>% 
  # microViz::tax_mutate(gram = left_join(tax_table(.) %>% data.frame(), gram,
  # by = c("Strain" = "ASV") %>% select(gram_neg_class))) %>% 
  microViz::comp_barplot(
    tax_transform_for_plot = "identity",
    tax_level = "intrinsic",
    label = "Day_of_Treatment", # name an alternative variable to label axes
    # palette = microViz::distinct_palette(3, pal = "brewerPlus"),
    n_taxa = 10, # give more taxa unique colours
    merge_other = FALSE, # split the "other" category to display alpha diversity
    bar_width = 0.9, # reduce the bar width to 70% of one row
    bar_outline_colour = "grey5",
    sample_order = "default") +
  ylab("Proportion - %") + xlab( " Day (Treatment) ") -> p_int_h2
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
p_int_h2 + theme_light() +facet_grid(Kingdom ~ Reactor_Treatment_Dose, scales = "free", space = "free_x", drop = TRUE) -> p_int_h2

p_int_h2 %>% 
  ggpubr::set_palette(p = ., palette = "jco") -> p_int_h2
```

```
## Scale for 'fill' is already present. Adding another scale for 'fill', which
## will replace the existing scale.
```

```r
p_int_h2 + theme(legend.position = "none") -> p_int_h2_no_leg

p_int_h2_no_leg
```

![](NRP72_draft_16S_Figs_files/figure-html/unnamed-chunk-79-1.png)<!-- -->

```r
p_int_h2 %>%  
  ggpubr::get_legend() %>% 
  ggpubr::as_ggplot() -> p_int_h2_leg

p_int_h2_leg
```

![](NRP72_draft_16S_Figs_files/figure-html/unnamed-chunk-79-2.png)<!-- -->


```r
p_int_h2$data %>%
  # filter(Model2 == "Human2") %>%
  # microViz::ps_mutate(Reactor_Treatment_Dose = paste0(Reactor_Treatment_Dose, "_", Fermentation)) %>% 
  # mutate(combo = paste0(Model2, "_" , Reactor_F1_F2)) %>%
  # mutate(combo = as.factor(combo)) %>%
  select(SAMPLE, top,OTU, Abundance, Day_of_Treatment, Reactor_Treatment_Dose) %>% 
  ggplot(.,
         aes(x = Day_of_Treatment, 
             y = Abundance,
             stratum = top, 
             alluvium = OTU,
             fill = top,
             label = top))+
  facet_grid(Reactor_Treatment_Dose  ~ .)+
  # scale_fill_manual(values = c("red", "yellow", "green3"))+
  geom_flow(stat = "alluvium", lode.guidance = "frontback", color = "darkgray")+
  geom_stratum() +
  theme_light() +
  # geom_text(stat = "alluvium", aes(label = top), lode.guidance = "frontback") +
  # scale_fill_manual(values = col_ASVs) +
  scale_color_manual(values  =  microViz::distinct_palette(3, pal = "greenArmytage")) +
  # scale_color_manual(values  = microViz::distinct_palette(10, pal = "kelly")) +
  # scale_fill_manual(values  = microViz::distinct_palette(10, pal = "kelly")) +
  ylab("Proportion - %") -> p_allu_human2_int

p_allu_human2_int %>% 
  ggpubr::set_palette(p = ., palette = "jco") -> p_allu_human2_int
```

```
## Scale for 'colour' is already present. Adding another scale for 'colour',
## which will replace the existing scale.
```

```r
p_allu_human2_int + theme(legend.position = "none") -> p_allu_human2_int_noleg 

p_allu_human2_int_noleg 
```

![](NRP72_draft_16S_Figs_files/figure-html/unnamed-chunk-80-1.png)<!-- -->

```r
p_allu_human2_int %>% 
  ggpubr::get_legend() %>% 
  ggpubr::as_ggplot() ->  p_allu_human2_int_leg
```

```r
ggpubr::ggarrange(p_allu_human2_noleg + xlab(NULL),
                  p_allu_human2_ASV_noleg + ylab(NULL) + xlab(NULL),
                  p_allu_human2_int_noleg +  coord_cartesian(ylim = c(50,100)),
                  p_allu_human2_int_noleg+ ylab(NULL),
                  align = "v", 
                  ncol = 2, nrow = 2,
                  heights = c(1, 0.6), widths = c(1,1),
                  common.legend = FALSE) -> p_allu_human2

p_allu_human2
```

![](NRP72_draft_16S_Figs_files/figure-html/unnamed-chunk-81-1.png)<!-- -->

```r
p_allu_human2 %>%
  export::graph2ppt(append = TRUE,
                    width = 317.48031496 * 1,
                    height = 0.618 * 317.48031496 * 1.75 , paper = "A4",  scaling = 2,
                    file = out_pptx)
```

```
## Exported graph as ~/Desktop/16S-human.pptx
```



```r
ggpubr::ggarrange(p_gram_h2_no_leg + ylab(NULL) + 
                    theme(axis.title.x=element_blank(),
                          axis.text.x=element_blank(),
                          axis.ticks.x=element_blank()),
                  p_hist_human2_fam_no_leg + 
                    theme(
                      strip.background = element_blank(),
                      strip.text.x = element_blank()
                    ) +
                    theme(axis.title.x=element_blank(),
                          axis.text.x=element_blank(),
                          axis.ticks.x=element_blank()),
                  p_hist_human2_strain_no_leg + 
                    theme(
                      strip.background = element_blank(),
                      strip.text.x = element_blank()
                    ) + ylab(NULL) + 
                    theme(axis.title.x=element_blank(),
                          axis.text.x=element_blank(),
                          axis.ticks.x=element_blank()),
                  p_int_h2_no_leg + ylab(NULL) + coord_cartesian(ylim = c(50,100)) +
                    theme(
                      strip.background = element_blank(),
                      strip.text.x = element_blank()
                    ), 
                  align = "v",
                  ncol = 1, 
                  heights = c(0.15, 0.4, 0.8,0.15),
                  common.legend = FALSE) -> all_plot_h2

all_plot_h2
```

![](NRP72_draft_16S_Figs_files/figure-html/unnamed-chunk-82-1.png)<!-- -->

```r
all_plot_h2 %>% 
  export::graph2ppt(append = TRUE, 
                    width = 317.48031496 * 3.5, 
                    height = 0.618 * 317.48031496 * 5 , paper = "A4",  scaling = 2,
                    file = out_pptx)
```

```
## Exported graph as ~/Desktop/16S-human.pptx
```

```r
p_gram_h2_leg  %>% 
  export::graph2ppt(append = TRUE, 
                    width = 317.48031496 * 1, 
                    height = 0.618 * 317.48031496 * 1 , paper = "A4",  scaling = 2,
                    file = out_pptx)
```

```
## Exported graph as ~/Desktop/16S-human.pptx
```

```r
p_hist_human2_fam_leg %>% 
  export::graph2ppt(append = TRUE, 
                    width = 317.48031496 * 1, 
                    height = 0.618 * 317.48031496 * 1 , paper = "A4",  scaling = 2,
                    file = out_pptx)
```

```
## Exported graph as ~/Desktop/16S-human.pptx
```

```r
p_hist_human2_strain_leg  %>% 
  export::graph2ppt(append = TRUE, 
                    width = 317.48031496 * 4, 
                    height = 0.618 * 317.48031496 * 2 , paper = "A3",  scaling = 2,
                    file = out_pptx)
```

```
## Exported graph as ~/Desktop/16S-human.pptx
```

```r
p_int_h2_leg  %>% 
  export::graph2ppt(append = TRUE, 
                    width = 317.48031496 * 1, 
                    height = 0.618 * 317.48031496 * 1 , paper = "A4",  scaling = 2,
                    file = out_pptx)
```

```
## Exported graph as ~/Desktop/16S-human.pptx
```

```r
p_gram_h2_no_leg$data %>% 
  select(OTU, Sample, Abundance, Treatment, Day_of_Treatment) %>% 
  xlsx::write.xlsx(file = out_xlsx,
                   append = TRUE,
                   sheetName = "p_gram_h2",
                   showNA = TRUE)
p_hist_human2_fam$data %>% 
  select(OTU, Sample, Abundance, Treatment, Day_of_Treatment) %>% 
  xlsx::write.xlsx(file = out_xlsx,
                   append = TRUE,
                   sheetName = "p_hist_human2_fam",
                   showNA = TRUE)
p_hist_human2_strain$data %>% 
  select(Sample, top ,Abundance, Treatment, Day_of_Treatment) %>% 
  pivot_wider(names_from = top, values_from = Abundance) %>% 
  xlsx::write.xlsx(file = out_xlsx,
                   append = TRUE,
                   sheetName = "p_hist_human2_strain",
                   showNA = TRUE)
p_int_h2$data %>% 
  select(OTU, Sample, Abundance, Treatment, Day_of_Treatment) %>% 
  xlsx::write.xlsx(file = out_xlsx,
                   append = TRUE,
                   sheetName = "p_int_h2",
                   showNA = TRUE)
```



```r
# ggpubr::ggarrange(p_gram_h2_no_leg + ylab(NULL) + 
#                     theme(axis.title.x=element_blank(),
#                           axis.text.x=element_blank(),
#                           axis.ticks.x=element_blank()),
#                   p_hist_human2_fam_no_leg + 
#                     theme(
#                       strip.background = element_blank(),
#                       strip.text.x = element_blank()
#                     ) +
#                     theme(axis.title.x=element_blank(),
#                           axis.text.x=element_blank(),
#                           axis.ticks.x=element_blank()),
#                   treat_speces + 
#                     theme(
#                       strip.background = element_blank(),
#                       strip.text.x = element_blank()
#                     ) + ylab(NULL) + 
#                     theme(axis.title.x=element_blank(),
#                           axis.text.x=element_blank(),
#                           axis.ticks.x=element_blank()),
#                   p_int_h2_no_leg + ylab(NULL) + coord_cartesian(ylim = c(50,100)) +
#                     theme(
#                       strip.background = element_blank(),
#                       strip.text.x = element_blank()
#                     ), 
#                   align = "hv",
#                   ncol = 1, 
#                   heights = c(0.15, 0.4, 0.8,0.15),
#                   common.legend = FALSE) -> all_plot_h2
# 
# all_plot_h2
# 
# all_plot_h2 %>% 
#   export::graph2ppt(append = TRUE, 
#                     width = 317.48031496 * 3.5, 
#                     height = 0.618 * 317.48031496 * 5 , paper = "A4",  scaling = 2,
#                     file = out_pptx)
```



### human1 - f1 :

#### gram:


```r
ps_rare_gram_int %>% 
  subset_samples(Model2 == "Human1" & Treatment != "DONOR" & Day_of_Treatment >= -1 & Fermentation == 1) %>% 
  # microViz::ps_mutate(Reactor_Treatment_Dose = paste0(Reactor_Treatment_Dose, "_", Fermentation)) %>% 
  transform_sample_counts(function(x) x/sum(x) * 100) %>%
  physeq_sel_tax_table(c("Kingdom","gram_stain")) %>% 
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
  ylab("Proportion - %") + xlab( " Day (Treatment) ") -> p_gram_h1_f1
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
p_gram_h1_f1 + facet_grid(Kingdom ~ Reactor_Treatment_Dose, scales = "free", space = "free_x", drop = TRUE) -> p_gram_h1_f1

p_gram_h1_f1 + theme_light() +theme(legend.position = "none") -> p_gram_h1_f1_no_leg

p_gram_h1_f1_no_leg
```

![](NRP72_draft_16S_Figs_files/figure-html/unnamed-chunk-84-1.png)<!-- -->

```r
p_gram_h1_f1 %>%  
  ggpubr::get_legend() %>% 
  ggpubr::as_ggplot() -> p_gram_h1_f1_leg
```



```r
p_gram_h1_f1$data %>%
  # filter(Model2 == "Human2") %>%
  mutate(top = factor(top, levels = c("Gram-positive", "Gram-negative", "unknown-gram"))) %>%
  # mutate(combo = paste0(Model2, "_" , Reactor_F1_F2)) %>%
  # mutate(combo = as.factor(combo)) %>%
  select(SAMPLE, top,OTU, Abundance, Day_of_Treatment, Reactor_Treatment_Dose) %>% 
  ggplot(.,
         aes(x = Day_of_Treatment, 
             y = Abundance,
             stratum = top, 
             alluvium = OTU,
             fill = top,
             label = top))+
  facet_grid(Reactor_Treatment_Dose  ~ .)+
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
  ylab("Proportion - %") -> p_allu_h1_f1_gram 

p_allu_h1_f1_gram + theme(legend.position = "none") -> p_allu_h1_f1_gram_noleg 

p_allu_h1_f1_gram_noleg 
```

![](NRP72_draft_16S_Figs_files/figure-html/unnamed-chunk-85-1.png)<!-- -->

```r
p_allu_h1_f1_gram %>% 
  ggpubr::get_legend() %>% 
  ggpubr::as_ggplot() ->  p_allu_h1_f1_gram_leg

p_allu_h1_f1_gram_leg
```

![](NRP72_draft_16S_Figs_files/figure-html/unnamed-chunk-85-2.png)<!-- -->



#### intrinsic:


```r
ps_rare_gram_int %>% 
  subset_samples(Model2 == "Human1" & Treatment != "DONOR" & Day_of_Treatment >= -1 & Fermentation == 1) %>% 
  microViz::ps_mutate(Reactor_Treatment_Dose = paste0(Reactor_Treatment_Dose, "_", Fermentation)) %>% 
  transform_sample_counts(function(x) x/sum(x) * 100) %>%
  physeq_sel_tax_table(c("Kingdom", "intrinsic")) %>% 
  # microbiome::transform("log") %>% 
  # speedyseq::tax_glom("intrinsic") %>% 
  # physeq_sel_tax_table(c("gram_stain","intrinsic")) %>% 
  # physeq_glom_rename(phyloseq = ., speedyseq = TRUE, rename_ASV = FALSE, taxrank = "gram_stain") %>% 
  # subset_taxa(Strain %in% top_taxa_per_group_h2) %>% # extract only the taxa to display - after percentage normalisation
  microViz::phyloseq_validate() %>%
  microViz::tax_fix() %>%
  microViz::ps_arrange(Day_of_Treatment) %>% 
  # microViz::tax_mutate(gram = left_join(tax_table(.) %>% data.frame(), gram,
  # by = c("Strain" = "ASV") %>% select(gram_neg_class))) %>% 
  microViz::comp_barplot(
    tax_transform_for_plot = "identity",
    tax_level = "intrinsic",
    label = "Day_of_Treatment", # name an alternative variable to label axes
    # palette = microViz::distinct_palette(3, pal = "brewerPlus"),
    n_taxa = 10, # give more taxa unique colours
    merge_other = FALSE, # split the "other" category to display alpha diversity
    bar_width = 0.9, # reduce the bar width to 70% of one row
    bar_outline_colour = "grey5",
    sample_order = "default") +
  ylab("Proportion - %") + xlab( " Day (Treatment) ") -> p_int_h1_f1
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
p_int_h1_f1 + facet_grid(Kingdom ~ Reactor_Treatment_Dose, scales = "free", space = "free_x", drop = TRUE) -> p_int_h1_f1

p_int_h1_f1 %>% 
  ggpubr::set_palette(p = ., palette = "jco") -> p_int_h1_f1
```

```
## Scale for 'fill' is already present. Adding another scale for 'fill', which
## will replace the existing scale.
```

```r
p_int_h1_f1 + theme_light() +theme(legend.position = "none") -> p_int_h1_f1_no_leg

p_int_h1_f1_no_leg
```

![](NRP72_draft_16S_Figs_files/figure-html/unnamed-chunk-86-1.png)<!-- -->

```r
p_int_h1_f1 %>%  
  ggpubr::get_legend() %>% 
  ggpubr::as_ggplot() -> p_int_h1_f1_leg

p_int_h1_f1_leg
```

![](NRP72_draft_16S_Figs_files/figure-html/unnamed-chunk-86-2.png)<!-- -->



```r
p_int_h1_f1$data %>%
  # filter(Model2 == "Human2") %>%
  # microViz::ps_mutate(Reactor_Treatment_Dose = paste0(Reactor_Treatment_Dose, "_", Fermentation)) %>% 
  # mutate(combo = paste0(Model2, "_" , Reactor_F1_F2)) %>%
  # mutate(combo = as.factor(combo)) %>%
  select(SAMPLE, top,OTU, Abundance, Day_of_Treatment, Reactor_Treatment_Dose) %>% 
  ggplot(.,
         aes(x = Day_of_Treatment, 
             y = Abundance,
             stratum = top, 
             alluvium = OTU,
             fill = top,
             label = top))+
  facet_grid(Reactor_Treatment_Dose  ~ .)+
  # scale_fill_manual(values = c("red", "yellow", "green3"))+
  geom_flow(stat = "alluvium", lode.guidance = "frontback", color = "darkgray")+
  geom_stratum() +
  theme_light() +
  # geom_text(stat = "alluvium", aes(label = top), lode.guidance = "frontback") +
  # scale_fill_manual(values = col_ASVs) +
  scale_color_manual(values  =  microViz::distinct_palette(3, pal = "greenArmytage")) +
  # scale_color_manual(values  = microViz::distinct_palette(10, pal = "kelly")) +
  # scale_fill_manual(values  = microViz::distinct_palette(10, pal = "kelly")) +
  ylab("Proportion - %") -> p_allu_h1_f1_int

p_allu_h1_f1_int %>% 
  ggpubr::set_palette(p = ., palette = "jco") -> p_allu_h1_f1_int
```

```
## Scale for 'colour' is already present. Adding another scale for 'colour',
## which will replace the existing scale.
```

```r
p_allu_h1_f1_int + theme(legend.position = "none") -> p_allu_h1_f1_int_noleg 

p_allu_h1_f1_int_noleg 
```

![](NRP72_draft_16S_Figs_files/figure-html/unnamed-chunk-87-1.png)<!-- -->

```r
p_allu_h1_f1_int %>% 
  ggpubr::get_legend() %>% 
  ggpubr::as_ggplot() -> p_allu_h1_f1_int_leg
```


```r
ggpubr::ggarrange(p_allu_human1_f1_noleg + xlab(NULL),
                  p_allu_human1_ASV_f1_noleg + ylab(NULL) + xlab(NULL),
                  p_allu_h1_f1_int_noleg +  coord_cartesian(ylim = c(50,100)),
                  p_allu_h1_f1_gram_noleg + ylab(NULL),
                  align = "v", 
                  ncol = 2, nrow = 2,
                  heights = c(1, 0.6), widths = c(1,1),
                  common.legend = FALSE) -> p_allu_human1_f1

p_allu_human1_f1
```

![](NRP72_draft_16S_Figs_files/figure-html/unnamed-chunk-88-1.png)<!-- -->

```r
p_allu_human1_f1 %>%
  export::graph2ppt(append = TRUE,
                    width = 317.48031496 * 1,
                    height = 0.618 * 317.48031496 * 1.75 , paper = "A4",  scaling = 2,
                    file = out_pptx)
```

```
## Exported graph as ~/Desktop/16S-human.pptx
```


```r
ggpubr::ggarrange(p_gram_h1_f1_no_leg + ylab(NULL) + 
                    theme(axis.title.x=element_blank(),
                          axis.text.x=element_blank(),
                          axis.ticks.x=element_blank()),
                  p_hist_human1_f1_fam_no_leg + 
                    theme(
                      strip.background = element_blank(),
                      strip.text.x = element_blank()
                    ) +
                    theme(axis.title.x=element_blank(),
                          axis.text.x=element_blank(),
                          axis.ticks.x=element_blank()),
                  p_hist_human1_f1_strain_no_leg + 
                    theme(
                      strip.background = element_blank(),
                      strip.text.x = element_blank()
                    ) + ylab(NULL) + 
                    theme(axis.title.x=element_blank(),
                          axis.text.x=element_blank(),
                          axis.ticks.x=element_blank()),
                  p_int_h1_f1_no_leg + ylab(NULL) + coord_cartesian(ylim = c(50,100)) +
                    theme(
                      strip.background = element_blank(),
                      strip.text.x = element_blank()
                    ), 
                  align = "v",
                  ncol = 1, 
                  heights = c(0.15, 0.4, 0.8,0.15),
                  common.legend = FALSE) -> all_plot_h1_f1

all_plot_h1_f1
```

![](NRP72_draft_16S_Figs_files/figure-html/unnamed-chunk-89-1.png)<!-- -->

```r
all_plot_h1_f1 %>% 
  export::graph2ppt(append = TRUE, 
                    width = 317.48031496 * 3.5, 
                    height = 0.618 * 317.48031496 * 5 , paper = "A4",  scaling = 2,
                    file = out_pptx)
```

```
## Exported graph as ~/Desktop/16S-human.pptx
```

```r
p_gram_h1_f1_leg  %>% 
  export::graph2ppt(append = TRUE, 
                    width = 317.48031496 * 1, 
                    height = 0.618 * 317.48031496 * 1 , paper = "A4",  scaling = 2,
                    file = out_pptx)
```

```
## Exported graph as ~/Desktop/16S-human.pptx
```

```r
p_hist_human1_f1_fam_leg  %>% 
  export::graph2ppt(append = TRUE, 
                    width = 317.48031496 * 1, 
                    height = 0.618 * 317.48031496 * 1 , paper = "A4",  scaling = 2,
                    file = out_pptx)
```

```
## Exported graph as ~/Desktop/16S-human.pptx
```

```r
p_hist_human1_f1_strain_leg  %>% 
  export::graph2ppt(append = TRUE, 
                    width = 317.48031496 * 4, 
                    height = 0.618 * 317.48031496 * 2 , paper = "A3",  scaling = 2,
                    file = out_pptx)
```

```
## Exported graph as ~/Desktop/16S-human.pptx
```

```r
p_int_h1_f1_leg  %>% 
  export::graph2ppt(append = TRUE, 
                    width = 317.48031496 * 1, 
                    height = 0.618 * 317.48031496 * 1 , paper = "A4",  scaling = 2,
                    file = out_pptx)
```

```
## Exported graph as ~/Desktop/16S-human.pptx
```

```r
p_gram_h1_f1$data %>% 
  select(OTU, Sample, Abundance, Treatment, Day_of_Treatment) %>% 
  xlsx::write.xlsx(file = out_xlsx,
                   append = TRUE,
                   sheetName = "p_gram_h1_f1",
                   showNA = TRUE)
p_hist_human1_f1_fam$data %>% 
  select(OTU, Sample, Abundance, Treatment, Day_of_Treatment) %>% 
  pivot_wider(names_from = OTU, values_from = Abundance) %>% 
  xlsx::write.xlsx(file = out_xlsx,
                   append = TRUE,
                   sheetName = "p_hist_human1_f1_fam",
                   showNA = TRUE)
p_hist_human1_f1_strain$data %>% 
  select(Sample, top ,Abundance, Treatment, Day_of_Treatment) %>% 
  pivot_wider(names_from = top, values_from = Abundance) %>% 
  xlsx::write.xlsx(file = out_xlsx,
                   append = TRUE,
                   sheetName = "p_hist_human1_f1_strain",
                   showNA = TRUE)
p_int_h1_f1$data %>% 
  select(OTU, Sample, Abundance, Treatment, Day_of_Treatment) %>% 
  pivot_wider(names_from = OTU, values_from = Abundance) %>% 
  xlsx::write.xlsx(file = out_xlsx,
                   append = TRUE,
                   sheetName = "p_int_h1_f1",
                   showNA = TRUE)
```


### human1 - f2 :

#### gram:


```r
ps_rare_gram_int %>% 
  subset_samples(Model2 == "Human1" & Treatment != "DONOR" & Day_of_Treatment >= -1 & Fermentation == 2) %>% 
  microViz::ps_mutate(Reactor_Treatment_Dose = paste0(Reactor_Treatment_Dose, "_", Fermentation)) %>% 
  transform_sample_counts(function(x) x/sum(x) * 100) %>%
  physeq_sel_tax_table(c("Kingdom","gram_stain")) %>% 
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
  ylab("Proportion - %") + xlab( " Day (Treatment) ") -> p_gram_h1_f2

p_gram_h1_f2 + facet_grid(Kingdom ~ Reactor_Treatment_Dose, scales = "free", space = "free_x", drop = TRUE) -> p_gram_h1_f2

p_gram_h1_f2 + theme_light() + theme(legend.position = "none") -> p_gram_h1_f2_no_leg

p_gram_h1_f2_no_leg
```

![](NRP72_draft_16S_Figs_files/figure-html/unnamed-chunk-90-1.png)<!-- -->

```r
p_gram_h1_f2 %>%  
  ggpubr::get_legend() %>% 
  ggpubr::as_ggplot() -> p_gram_h1_f2_leg
```



```r
p_gram_h1_f2$data %>%
  # filter(Model2 == "Human2") %>%
  mutate(top = factor(top, levels = c("Gram-positive", "Gram-negative", "unknown-gram"))) %>%
  # mutate(combo = paste0(Model2, "_" , Reactor_F1_F2)) %>%
  # mutate(combo = as.factor(combo)) %>%
  select(SAMPLE, top,OTU, Abundance, Day_of_Treatment, Reactor_Treatment_Dose) %>% 
  ggplot(.,
         aes(x = Day_of_Treatment, 
             y = Abundance,
             stratum = top, 
             alluvium = OTU,
             fill = top,
             label = top))+
  facet_grid(Reactor_Treatment_Dose  ~ .)+
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
  ylab("Proportion - %") -> p_allu_h1_f2_gram 

p_allu_h1_f2_gram + theme(legend.position = "none") -> p_allu_h1_f2_gram_noleg 

p_allu_h1_f2_gram_noleg 
```

![](NRP72_draft_16S_Figs_files/figure-html/unnamed-chunk-91-1.png)<!-- -->

```r
p_allu_h1_f2_gram %>% 
  ggpubr::get_legend() %>% 
  ggpubr::as_ggplot() ->  p_allu_h1_f2_gram_leg

p_allu_h1_f2_gram_leg
```

![](NRP72_draft_16S_Figs_files/figure-html/unnamed-chunk-91-2.png)<!-- -->

#### intrinsic:


```r
ps_rare_gram_int %>% 
  subset_samples(Model2 == "Human1" & Treatment != "DONOR" & Day_of_Treatment >= -1 & Fermentation == 2) %>% 
  microViz::ps_mutate(Reactor_Treatment_Dose = paste0(Reactor_Treatment_Dose, "_", Fermentation)) %>% 
  transform_sample_counts(function(x) x/sum(x) * 100) %>%
  physeq_sel_tax_table(c("Kingdom", "intrinsic")) %>% 
  # microbiome::transform("log") %>% 
  # speedyseq::tax_glom("intrinsic") %>% 
  # physeq_sel_tax_table(c("gram_stain","intrinsic")) %>% 
  # physeq_glom_rename(phyloseq = ., speedyseq = TRUE, rename_ASV = FALSE, taxrank = "gram_stain") %>% 
  # subset_taxa(Strain %in% top_taxa_per_group_h2) %>% # extract only the taxa to display - after percentage normalisation
  microViz::phyloseq_validate() %>%
  microViz::tax_fix() %>%
  microViz::ps_arrange(Day_of_Treatment) %>% 
  # microViz::tax_mutate(gram = left_join(tax_table(.) %>% data.frame(), gram,
  # by = c("Strain" = "ASV") %>% select(gram_neg_class))) %>% 
  microViz::comp_barplot(
    tax_transform_for_plot = "identity",
    tax_level = "intrinsic",
    label = "Day_of_Treatment", # name an alternative variable to label axes
    # palette = microViz::distinct_palette(3, pal = "brewerPlus"),
    n_taxa = 10, # give more taxa unique colours
    merge_other = FALSE, # split the "other" category to display alpha diversity
    bar_width = 0.9, # reduce the bar width to 70% of one row
    bar_outline_colour = "grey5",
    sample_order = "default") +
  ylab("Proportion - %") + xlab( " Day (Treatment) ") -> p_int_h1_f2
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
p_int_h1_f2 + facet_grid(Kingdom ~ Reactor_Treatment_Dose, scales = "free", space = "free_x", drop = TRUE) -> p_int_h1_f2

p_int_h1_f2 %>% 
  ggpubr::set_palette(p = ., palette = "jco") -> p_int_h1_f2
```

```
## Scale for 'fill' is already present. Adding another scale for 'fill', which
## will replace the existing scale.
```

```r
p_int_h1_f2 + theme_light() +theme(legend.position = "none") -> p_int_h1_f2_no_leg

p_int_h1_f2_no_leg
```

![](NRP72_draft_16S_Figs_files/figure-html/unnamed-chunk-92-1.png)<!-- -->

```r
p_int_h1_f2 %>%  
  ggpubr::get_legend() %>% 
  ggpubr::as_ggplot() -> p_int_h1_f2_leg

p_int_h1_f2_leg
```

![](NRP72_draft_16S_Figs_files/figure-html/unnamed-chunk-92-2.png)<!-- -->



```r
p_int_h1_f2$data %>%
  # filter(Model2 == "Human2") %>%
  # microViz::ps_mutate(Reactor_Treatment_Dose = paste0(Reactor_Treatment_Dose, "_", Fermentation)) %>% 
  # mutate(combo = paste0(Model2, "_" , Reactor_F1_F2)) %>%
  # mutate(combo = as.factor(combo)) %>%
  select(SAMPLE, top,OTU, Abundance, Day_of_Treatment, Reactor_Treatment_Dose) %>% 
  ggplot(.,
         aes(x = Day_of_Treatment, 
             y = Abundance,
             stratum = top, 
             alluvium = OTU,
             fill = top,
             label = top))+
  facet_grid(Reactor_Treatment_Dose  ~ .)+
  # scale_fill_manual(values = c("red", "yellow", "green3"))+
  geom_flow(stat = "alluvium", lode.guidance = "frontback", color = "darkgray")+
  geom_stratum() +
  theme_light() +
  # geom_text(stat = "alluvium", aes(label = top), lode.guidance = "frontback") +
  # scale_fill_manual(values = col_ASVs) +
  scale_color_manual(values  =  microViz::distinct_palette(3, pal = "greenArmytage")) +
  # scale_color_manual(values  = microViz::distinct_palette(10, pal = "kelly")) +
  # scale_fill_manual(values  = microViz::distinct_palette(10, pal = "kelly")) +
  ylab("Proportion - %") -> p_allu_h1_f2_int

p_allu_h1_f2_int %>% 
  ggpubr::set_palette(p = ., palette = "jco") -> p_allu_h1_f2_int
```

```
## Scale for 'colour' is already present. Adding another scale for 'colour',
## which will replace the existing scale.
```

```r
p_allu_h1_f2_int + theme(legend.position = "none") -> p_allu_h1_f2_int_noleg 

p_allu_h1_f2_int_noleg 
```

![](NRP72_draft_16S_Figs_files/figure-html/unnamed-chunk-93-1.png)<!-- -->

```r
p_allu_h1_f2_int %>% 
  ggpubr::get_legend() %>% 
  ggpubr::as_ggplot() -> p_allu_h1_f2_int_leg
```



```r
# ggpubr::ggarrange(p_allu_human2_noleg,
#                   p_allu_human2_ASV_noleg + ylab(NULL),
#                   align = "v",
#                   ncol = 2, 
#                   heights = c(1, 1),
#                   common.legend = FALSE) -> p_allu_human2
# 
# p_allu_human2

# p_allu_human2 %>% 
#   export::graph2ppt(append = TRUE, 
#                     width = 317.48031496 * 1, 
#                     height = 0.618 * 317.48031496 * 1 , paper = "A4",  scaling = 2,
#                     file = out_pptx)
```


```r
ggpubr::ggarrange(p_allu_human1_f2_noleg + xlab(NULL),
                  p_allu_human1_ASV_f2_noleg + ylab(NULL) + xlab(NULL),
                  p_allu_h1_f2_int_noleg +  coord_cartesian(ylim = c(50,100)),
                  p_allu_h1_f2_gram_noleg + ylab(NULL),
                  align = "v", 
                  ncol = 2, nrow = 2,
                  heights = c(1, 0.6), widths = c(1,1),
                  common.legend = FALSE) -> p_allu_human1_f2

p_allu_human1_f2
```

![](NRP72_draft_16S_Figs_files/figure-html/unnamed-chunk-95-1.png)<!-- -->

```r
p_allu_human1_f2 %>%
  export::graph2ppt(append = TRUE,
                    width = 317.48031496 * 1,
                    height = 0.618 * 317.48031496 * 1.75 , paper = "A4",  scaling = 2,
                    file = out_pptx)
```

```
## Exported graph as ~/Desktop/16S-human.pptx
```


```r
ggpubr::ggarrange(p_allu_human1_f2_noleg,
                  p_allu_human1_ASV_f2_noleg + ylab(NULL),
                  
                  align = "v",
                  ncol = 2, 
                  heights = c(1, 1),
                  common.legend = FALSE) -> p_allu_human1_f2

p_allu_human1_f2
```

![](NRP72_draft_16S_Figs_files/figure-html/unnamed-chunk-96-1.png)<!-- -->

```r
p_allu_human1_f2 %>%
  export::graph2ppt(append = TRUE,
                    width = 317.48031496 * 1,
                    height = 0.618 * 317.48031496 * 1 , paper = "A4",  scaling = 2,
                    file = out_pptx)
```

```
## Exported graph as ~/Desktop/16S-human.pptx
```




```r
ggpubr::ggarrange(p_gram_h1_f2_no_leg + ylab(NULL) + 
                    theme(axis.title.x=element_blank(),
                          axis.text.x=element_blank(),
                          axis.ticks.x=element_blank()),
                  p_hist_human1_f2_fam_no_leg + 
                    theme(
                      strip.background = element_blank(),
                      strip.text.x = element_blank()
                    ) +
                    theme(axis.title.x=element_blank(),
                          axis.text.x=element_blank(),
                          axis.ticks.x=element_blank()),
                  p_hist_human1_f2_strain_no_leg + 
                    theme(
                      strip.background = element_blank(),
                      strip.text.x = element_blank()
                    ) + ylab(NULL) + 
                    theme(axis.title.x=element_blank(),
                          axis.text.x=element_blank(),
                          axis.ticks.x=element_blank()),
                  p_int_h1_f2_no_leg + ylab(NULL) + coord_cartesian(ylim = c(50,100)) +
                    theme(
                      strip.background = element_blank(),
                      strip.text.x = element_blank()
                    ), 
                  align = "v",
                  ncol = 1, 
                  heights = c(0.15, 0.4, 0.8,0.15),
                  common.legend = FALSE) -> all_plot_h1_f2

all_plot_h1_f2
```

![](NRP72_draft_16S_Figs_files/figure-html/unnamed-chunk-97-1.png)<!-- -->

```r
all_plot_h1_f2 %>% 
  export::graph2ppt(append = TRUE, 
                    width = 317.48031496 * 3.5, 
                    height = 0.618 * 317.48031496 * 5 , paper = "A4",  scaling = 2,
                    file = out_pptx)
```

```
## Exported graph as ~/Desktop/16S-human.pptx
```

```r
p_gram_h1_f2_leg  %>% 
  export::graph2ppt(append = TRUE, 
                    width = 317.48031496 * 1, 
                    height = 0.618 * 317.48031496 * 1 , paper = "A4",  scaling = 2,
                    file = out_pptx)
```

```
## Exported graph as ~/Desktop/16S-human.pptx
```

```r
p_hist_human1_f2_fam_leg  %>% 
  export::graph2ppt(append = TRUE, 
                    width = 317.48031496 * 1, 
                    height = 0.618 * 317.48031496 * 1 , paper = "A4",  scaling = 2,
                    file = out_pptx)
```

```
## Exported graph as ~/Desktop/16S-human.pptx
```

```r
p_hist_human1_f2_strain_leg  %>% 
  export::graph2ppt(append = TRUE, 
                    width = 317.48031496 * 4, 
                    height = 0.618 * 317.48031496 * 2 , paper = "A3",  scaling = 2,
                    file = out_pptx)
```

```
## Exported graph as ~/Desktop/16S-human.pptx
```

```r
p_int_h1_f2_leg  %>% 
  export::graph2ppt(append = TRUE, 
                    width = 317.48031496 * 1, 
                    height = 0.618 * 317.48031496 * 1 , paper = "A4",  scaling = 2,
                    file = out_pptx)
```

```
## Exported graph as ~/Desktop/16S-human.pptx
```

```r
p_gram_h1_f2$data %>% 
  select(OTU, Sample, Abundance, Treatment, Day_of_Treatment) %>% 
  xlsx::write.xlsx(file = out_xlsx,
                   append = TRUE,
                   sheetName = "p_gram_h1_f2",
                   showNA = TRUE)

p_hist_human1_f2_fam$data %>% 
  select(OTU, Sample, Abundance, Treatment, Day_of_Treatment) %>% 
  xlsx::write.xlsx(file = out_xlsx,
                   append = TRUE,
                   sheetName = "p_hist_human1_f2_fam",
                   showNA = TRUE)

p_hist_human1_f2_strain$data %>% 
  select(Sample, top ,Abundance, Treatment, Day_of_Treatment) %>% 
  pivot_wider(names_from = top, values_from = Abundance) %>% 
  xlsx::write.xlsx(file = out_xlsx,
                   append = TRUE,
                   sheetName = "p_hist_human1_f2_strain",
                   showNA = TRUE)

p_int_h1_f2$data %>% 
  select(OTU, Sample, Abundance, Treatment, Day_of_Treatment) %>% 
  xlsx::write.xlsx(file = out_xlsx,
                   append = TRUE,
                   sheetName = "p_int_h1_f2",
                   showNA = TRUE)

# htmlwidgets::saveWidget(plotly::as_widget(plotly::ggplotly(p_int_h1_f2)),
#                         "~/Desktop/index.html")
```

-   explore stability metrics: <https://microbiome.github.io/tutorials/Stability.html> to detect taxa impacted. or filter by coefficient variation. or Hannah's classification

Bimodal population distribution across individuals is often associated with instable intermediate abundances within individuals:


```r
# Use relative abundances
pseq <- microbiome::transform(ps_rare %>%  subset_samples(Model2 == "Human2"), "compositional")

# Merge rare taxa to speed up examples
pseq <- microbiome::aggregate_rare(ps_rare, level = "Strain", detection = .1/100, prevalence = 10/100)

# For cross-sectional analysis, include
# only the treated points.

pseq0 <- subset_samples(pseq, Day_of_Treatment > 1)

# intermediate_stability: within subjects:
pseq %>% 
  microViz::ps_mutate(subject = Reactor,
                      time = Day_of_Treatment) %>% 
  intermediate_stability(., #output = "scores", 
                         output = "scores") -> intermediate.stability

intermediate.stability %>% 
  head()
```


```r
# Bimodality is better estimated from abundances at log scale (such as CLR)
pseq0.clr <- microbiome::transform(pseq0, "clr")

set.seed(4433)

# bimodality: between subjects
bimodality.score <- bimodality(pseq0.clr, method = "potential_analysis",
                               bs.iter = 999, peak.threshold = 5,
                               min.density = 5)

bimodality.score %>% 
  head()
```


```r
taxa <- taxa(pseq0)
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
x <- abundances(microbiome::transform(pseq, "clr"))[tax,]

# Bootstrapped potential analysis to identify potential minima
# in practice, use more bootstrap iterations
potential.minima <- potential_analysis(x, bs.iter = 10)$minima

# Same with earlywarnings package (without bootstrap ie. less robust)
# library(earlywarnings)
# res <- livpotential_ews(x)$min.points

# Identify the potential minimum locations as tipping point candidates
tipping.samples <- sapply(potential.minima, function (m) {names(which.min(abs(sort(x) - m)))})
tipping.point <- abundances(pseq)[tax, tipping.samples]
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
##  [1] ggalluvial_0.12.3    gdtools_0.2.4        ampvis2_2.7.29      
##  [4] speedyseq_0.5.3.9018 reshape2_1.4.4       scales_1.2.1        
##  [7] compositions_2.0-4   nlme_3.1-159         phyloseq_1.40.0     
## [10] forcats_0.5.2        stringr_1.4.1        dplyr_1.0.10        
## [13] purrr_0.3.5          readr_2.1.3          tidyr_1.2.1         
## [16] tibble_3.1.8         ggplot2_3.3.6        tidyverse_1.3.2     
## 
## loaded via a namespace (and not attached):
##   [1] readxl_1.4.1           uuid_1.1-0             backports_1.4.1       
##   [4] systemfonts_1.0.4      plyr_1.8.7             igraph_1.3.5          
##   [7] lazyeval_0.2.2         splines_4.2.1          fantaxtic_0.1.0       
##  [10] GenomeInfoDb_1.32.4    digest_0.6.29          foreach_1.5.2         
##  [13] htmltools_0.5.3        fansi_1.0.3            magrittr_2.0.3        
##  [16] xlsx_0.6.5             googlesheets4_1.0.1    cluster_2.1.4         
##  [19] openxlsx_4.2.5         tzdb_0.3.0             Biostrings_2.64.1     
##  [22] extrafont_0.18         modelr_0.1.9           bayesm_3.1-4          
##  [25] officer_0.4.4          extrafontdb_1.0        colorspace_2.0-3      
##  [28] rvest_1.0.3            ggrepel_0.9.1          textshaping_0.3.6     
##  [31] haven_2.5.1            xfun_0.33              crayon_1.5.2          
##  [34] RCurl_1.98-1.9         jsonlite_1.8.2         survival_3.4-0        
##  [37] iterators_1.0.14       ape_5.6-2              glue_1.6.2            
##  [40] rvg_0.2.5              gtable_0.3.1           gargle_1.2.1          
##  [43] zlibbioc_1.42.0        XVector_0.36.0         V8_4.2.1              
##  [46] car_3.1-0              Rttf2pt1_1.3.10        Rhdf5lib_1.18.2       
##  [49] BiocGenerics_0.42.0    DEoptimR_1.0-11        abind_1.4-5           
##  [52] DBI_1.1.3              randomcoloR_1.1.0.1    rstatix_0.7.0         
##  [55] Rcpp_1.0.9             xtable_1.8-4           viridisLite_0.4.1     
##  [58] stats4_4.2.1           htmlwidgets_1.5.4      httr_1.4.4            
##  [61] RColorBrewer_1.1-3     ellipsis_0.3.2         rJava_1.0-6           
##  [64] pkgconfig_2.0.3        farver_2.1.1           sass_0.4.2            
##  [67] dbplyr_2.2.1           utf8_1.2.2             here_1.0.1            
##  [70] labeling_0.4.2         tidyselect_1.2.0       rlang_1.0.6           
##  [73] munsell_0.5.0          cellranger_1.1.0       tools_4.2.1           
##  [76] cachem_1.0.6           cli_3.4.1              generics_0.1.3        
##  [79] devEMF_4.1             ade4_1.7-19            export_0.3.0          
##  [82] broom_1.0.1            evaluate_0.17          biomformat_1.24.0     
##  [85] fastmap_1.1.0          yaml_2.3.5             ragg_1.2.3            
##  [88] knitr_1.40             fs_1.5.2               zip_2.2.1             
##  [91] robustbase_0.95-0      rgl_0.110.2            xml2_1.3.3            
##  [94] compiler_4.2.1         rstudioapi_0.14        curl_4.3.3            
##  [97] plotly_4.10.0          ggsignif_0.6.3         reprex_2.0.2          
## [100] bslib_0.4.0            stringi_1.7.8          highr_0.9             
## [103] stargazer_5.2.3        lattice_0.20-45        Matrix_1.5-1          
## [106] ggsci_2.9              microbiome_1.19.1      vegan_2.6-2           
## [109] permute_0.9-7          tensorA_0.36.2         multtest_2.52.0       
## [112] vctrs_0.4.2            pillar_1.8.1           lifecycle_1.0.3       
## [115] rhdf5filters_1.8.0     microViz_0.9.2         jquerylib_0.1.4       
## [118] data.table_1.14.2      cowplot_1.1.1          bitops_1.0-7          
## [121] flextable_0.8.2        R6_2.5.1               IRanges_2.30.1        
## [124] codetools_0.2-18       MASS_7.3-58.1          assertthat_0.2.1      
## [127] xlsxjars_0.6.1         rhdf5_2.40.0           rprojroot_2.0.3       
## [130] withr_2.5.0            S4Vectors_0.34.0       GenomeInfoDbData_1.2.8
## [133] mgcv_1.8-40            parallel_4.2.1         hms_1.1.2             
## [136] grid_4.2.1             rmarkdown_2.16         carData_3.0-5         
## [139] googledrive_2.0.0      Rtsne_0.16             ggpubr_0.4.0          
## [142] Biobase_2.56.0         lubridate_1.8.0        base64enc_0.1-3
```
