---
title: "NTP72 - 16S - beta - human draft "
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
  facet_grid(distance ~ ., scales = "free") -> pcoa_all

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
<div id="htmlwidget-9eceb20dc2a931752ae4" style="width:100%;height:auto;" class="datatables html-widget"></div>
<script type="application/json" data-for="htmlwidget-9eceb20dc2a931752ae4">{"x":{"filter":"none","vertical":false,"data":[["1","2"],["bjaccard","aichinson"],["Axis.1   [36%]","Axis.1   [41%]"],["bjaccard","aichinson"],["Axis.2   [17.4%]","Axis.2   [15.5%]"]],"container":"<table class=\"display\">\n  <thead>\n    <tr>\n      <th> <\/th>\n      <th>.id...1<\/th>\n      <th>V1...2<\/th>\n      <th>.id...3<\/th>\n      <th>V1...4<\/th>\n    <\/tr>\n  <\/thead>\n<\/table>","options":{"columnDefs":[{"orderable":false,"targets":0}],"order":[],"autoWidth":false,"orderClasses":false}},"evals":[],"jsHooks":[]}</script>
```

```r
expl_var %>%
  data.frame() %>% 
  export::table2ppt(append = TRUE,
                    file = out_pptx)
```

```
## Exported table as ~/Desktop/16S-human-beta.pptx
```

```{=html}
<template id="1366e946-8350-4ce0-a6b6-217c6360ccf2"><style>
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
</style><div class="tabwid"><style>.cl-be16d3e2{}.cl-be118a22{font-family:'Helvetica';font-size:12pt;font-weight:bold;font-style:normal;text-decoration:none;color:rgba(0, 0, 0, 1.00);background-color:transparent;}.cl-be118a2c{font-family:'Helvetica';font-size:12pt;font-weight:normal;font-style:normal;text-decoration:none;color:rgba(0, 0, 0, 1.00);background-color:transparent;}.cl-be13e970{margin:0;text-align:left;border-bottom: 0 solid rgba(0, 0, 0, 1.00);border-top: 0 solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);padding-bottom:5pt;padding-top:5pt;padding-left:5pt;padding-right:5pt;line-height: 1;background-color:transparent;}.cl-be13f8de{width:2.14in;background-color:transparent;vertical-align: middle;border-bottom: 1pt solid rgba(0, 0, 0, 1.00);border-top: 1pt solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-be13f8e8{width:2.14in;background-color:transparent;vertical-align: middle;border-bottom: 0 solid rgba(0, 0, 0, 1.00);border-top: 0 solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-be13f8f2{width:2.14in;background-color:transparent;vertical-align: middle;border-bottom: 1pt solid rgba(0, 0, 0, 1.00);border-top: 0 solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}</style><table class='cl-be16d3e2'><thead><tr style="overflow-wrap:break-word;"><td class="cl-be13f8de"><p class="cl-be13e970"><span class="cl-be118a22">.id...1</span></p></td><td class="cl-be13f8de"><p class="cl-be13e970"><span class="cl-be118a22">V1...2</span></p></td><td class="cl-be13f8de"><p class="cl-be13e970"><span class="cl-be118a22">.id...3</span></p></td><td class="cl-be13f8de"><p class="cl-be13e970"><span class="cl-be118a22">V1...4</span></p></td></tr></thead><tbody><tr style="overflow-wrap:break-word;"><td class="cl-be13f8e8"><p class="cl-be13e970"><span class="cl-be118a2c">bjaccard</span></p></td><td class="cl-be13f8e8"><p class="cl-be13e970"><span class="cl-be118a2c">Axis.1   [36%]</span></p></td><td class="cl-be13f8e8"><p class="cl-be13e970"><span class="cl-be118a2c">bjaccard</span></p></td><td class="cl-be13f8e8"><p class="cl-be13e970"><span class="cl-be118a2c">Axis.2   [17.4%]</span></p></td></tr><tr style="overflow-wrap:break-word;"><td class="cl-be13f8f2"><p class="cl-be13e970"><span class="cl-be118a2c">aichinson</span></p></td><td class="cl-be13f8f2"><p class="cl-be13e970"><span class="cl-be118a2c">Axis.1   [41%]</span></p></td><td class="cl-be13f8f2"><p class="cl-be13e970"><span class="cl-be118a2c">aichinson</span></p></td><td class="cl-be13f8f2"><p class="cl-be13e970"><span class="cl-be118a2c">Axis.2   [15.5%]</span></p></td></tr></tbody></table></div></template>
<div class="flextable-shadow-host" id="8e9af2dd-12a7-4cc7-a1b4-4016a90d5af3"></div>
<script>
var dest = document.getElementById("8e9af2dd-12a7-4cc7-a1b4-4016a90d5af3");
var template = document.getElementById("1366e946-8350-4ce0-a6b6-217c6360ccf2");
var caption = template.content.querySelector("caption");
var fantome = dest.attachShadow({mode: 'open'});
var templateContent = template.content;
fantome.appendChild(templateContent);
</script>

```

### Human1 - F1:


```r
ps_rare %>%
  subset_samples(Model2 == "Human1" & Treatment != "DONOR" & Fermentation == 1) -> ps_tmp_forbdiv
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
## Run 0 stress 0.08004383 
## Run 1 stress 0.08420444 
## Run 2 stress 0.08588613 
## Run 3 stress 0.08072561 
## Run 4 stress 0.08489772 
## Run 5 stress 0.09371116 
## Run 6 stress 0.09813696 
## Run 7 stress 0.1000701 
## Run 8 stress 0.08420213 
## Run 9 stress 0.09931823 
## Run 10 stress 0.08419882 
## Run 11 stress 0.08195215 
## Run 12 stress 0.09965269 
## Run 13 stress 0.08419884 
## Run 14 stress 0.08021497 
## ... Procrustes: rmse 0.01033182  max resid 0.05873964 
## Run 15 stress 0.09875933 
## Run 16 stress 0.0935 
## Run 17 stress 0.08418989 
## Run 18 stress 0.08195212 
## Run 19 stress 0.08588491 
## Run 20 stress 0.08029931 
## ... Procrustes: rmse 0.01499125  max resid 0.0927014 
## *** No convergence -- monoMDS stopping criteria:
##      1: no. of iterations >= maxit
##     18: stress ratio > sratmax
##      1: scale factor of the gradient < sfgrmin
## [1] "aichinson"
## Run 0 stress 0.06667031 
## Run 1 stress 0.0881507 
## Run 2 stress 0.06667031 
## ... Procrustes: rmse 4.808409e-06  max resid 1.546332e-05 
## ... Similar to previous best
## Run 3 stress 0.06667031 
## ... Procrustes: rmse 5.346959e-06  max resid 2.469792e-05 
## ... Similar to previous best
## Run 4 stress 0.06667031 
## ... Procrustes: rmse 2.773469e-06  max resid 9.304502e-06 
## ... Similar to previous best
## Run 5 stress 0.09148286 
## Run 6 stress 0.09593656 
## Run 7 stress 0.0716824 
## Run 8 stress 0.0904375 
## Run 9 stress 0.1000976 
## Run 10 stress 0.0716824 
## Run 11 stress 0.0716824 
## Run 12 stress 0.0911493 
## Run 13 stress 0.06677809 
## ... Procrustes: rmse 0.002777525  max resid 0.01390932 
## Run 14 stress 0.09598469 
## Run 15 stress 0.0972167 
## Run 16 stress 0.09792991 
## Run 17 stress 0.06667031 
## ... New best solution
## ... Procrustes: rmse 3.002917e-06  max resid 1.187582e-05 
## ... Similar to previous best
## Run 18 stress 0.08999087 
## Run 19 stress 0.0718436 
## Run 20 stress 0.07170616 
## *** Solution reached
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
<script type="application/json" data-for="htmlwidget-f30d7598595affc4f25a">{"x":{"filter":"none","vertical":false,"data":[["ASV0012","ASV0051","ASV0010","ASV0015","ASV0106","ASV0124","ASV0136","ASV0202"],[-0.855235232721098,0.643808371532634,-0.858869634025459,0.746869100286677,-0.74562321070153,-0.858290482973094,-0.765802345012549,-0.663062220104871],[-0.0541313992033875,0.296923192688876,-0.375818735955723,0.297037592051541,0.361776052164635,-0.071384652306359,-0.347505820239148,-0.283546138710583],[0.001,0.001,0.001,0.001,0.001,0.001,0.001,0.001],[0.734357511667227,0.502652601612058,0.878896770546384,0.646044784054807,0.686835884256687,0.741758321747086,0.707213526726803,0.520049920508081],["Bacteria","Bacteria","Bacteria","Bacteria","Bacteria","Bacteria","Bacteria","Bacteria"],["Firmicutes","Firmicutes","Firmicutes","Firmicutes","Firmicutes","Firmicutes","Actinobacteriota","Firmicutes"],["Clostridia","Bacilli","Clostridia","Negativicutes","Clostridia","Clostridia","Coriobacteriia","Clostridia"],["Oscillospirales","Lactobacillales","Lachnospirales","Veillonellales-Selenomonadales","Oscillospirales","Monoglobales","Coriobacteriales","Oscillospirales"],["Ruminococcaceae","Lactobacillaceae","Lachnospiraceae","Veillonellaceae","Oscillospiraceae","Monoglobaceae","Eggerthellaceae","Butyricicoccaceae"],[null,null,null,null,null,null,null,null],[null,null,null,null,null,null,null,null],[null,null,null,null,null,null,null,null],[0.00288888888888889,0.00288888888888889,0.00288888888888889,0.00288888888888889,0.00288888888888889,0.00288888888888889,0.00288888888888889,0.00288888888888889]],"container":"<table class=\"display\">\n  <thead>\n    <tr>\n      <th> <\/th>\n      <th>Axis.1<\/th>\n      <th>Axis.2<\/th>\n      <th>pval<\/th>\n      <th>r<\/th>\n      <th>Kingdom<\/th>\n      <th>Phylum<\/th>\n      <th>Class<\/th>\n      <th>Order<\/th>\n      <th>tax_rank_plot<\/th>\n      <th>Genus<\/th>\n      <th>Species<\/th>\n      <th>Strain<\/th>\n      <th>pval.adj<\/th>\n    <\/tr>\n  <\/thead>\n<\/table>","options":{"columnDefs":[{"className":"dt-right","targets":[1,2,3,4,13]},{"orderable":false,"targets":0}],"order":[],"autoWidth":false,"orderClasses":false}},"evals":[],"jsHooks":[]}</script>
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
<script type="application/json" data-for="htmlwidget-b85bf888dc62471906a9">{"x":{"filter":"none","vertical":false,"data":[["1","2","3","4","5"],["Succinat_mM","Formiat_mM","Acetat_mM","Butyrat_mM","Total_SCFA_mM"],[0.25387911017642,0.397146069924643,-0.931986035502789,-0.888849712418429,-0.873816422886923],[0.299944582034044,-0.283443016019585,-0.149688627900526,-0.290348308376891,-0.025036888976522],["Succinat_mM","Formiat_mM","Acetat_mM","Butyrat_mM","Total_SCFA_mM"],[0.003,0.001,0.001,0.001,0.001],[0.154421354875548,0.238064944186868,0.891004655694948,0.874355951443646,0.76418198671652]],"container":"<table class=\"display\">\n  <thead>\n    <tr>\n      <th> <\/th>\n      <th>label<\/th>\n      <th>Axis.1<\/th>\n      <th>Axis.2<\/th>\n      <th>id<\/th>\n      <th>pval<\/th>\n      <th>r<\/th>\n    <\/tr>\n  <\/thead>\n<\/table>","options":{"columnDefs":[{"className":"dt-right","targets":[2,3,5,6]},{"orderable":false,"targets":0}],"order":[],"autoWidth":false,"orderClasses":false}},"evals":[],"jsHooks":[]}</script>
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
<script type="application/json" data-for="htmlwidget-88fdd2a34d502a66744b">{"x":{"filter":"none","vertical":false,"data":[["ASV0012","ASV0051","ASV0010","ASV0015","ASV0106","ASV0124","ASV0136","ASV0202"],[-0.807869018332967,0.501736684674981,-0.789882571827549,0.681190472422282,-0.77695235901902,-0.755097137780456,-0.742905154586156,-0.601896083110978],[-0.241788609803216,0.624109179553861,-0.429380612178577,0.339414799731832,-0.132971234268199,-0.353917429832527,-0.352571634200053,-0.35347163418689],[0.001,0.001,0.001,0.001,0.001,0.001,0.001,0.001],[0.711114082612844,0.641251968752036,0.808282187391753,0.579222865995891,0.621336317328028,0.695429234623499,0.676214825953176,0.487221091039088],["Bacteria","Bacteria","Bacteria","Bacteria","Bacteria","Bacteria","Bacteria","Bacteria"],["Firmicutes","Firmicutes","Firmicutes","Firmicutes","Firmicutes","Firmicutes","Actinobacteriota","Firmicutes"],["Clostridia","Bacilli","Clostridia","Negativicutes","Clostridia","Clostridia","Coriobacteriia","Clostridia"],["Oscillospirales","Lactobacillales","Lachnospirales","Veillonellales-Selenomonadales","Oscillospirales","Monoglobales","Coriobacteriales","Oscillospirales"],["Ruminococcaceae","Lactobacillaceae","Lachnospiraceae","Veillonellaceae","Oscillospiraceae","Monoglobaceae","Eggerthellaceae","Butyricicoccaceae"],[null,null,null,null,null,null,null,null],[null,null,null,null,null,null,null,null],[null,null,null,null,null,null,null,null],[0.00305882352941176,0.00305882352941176,0.00305882352941176,0.00305882352941176,0.00305882352941176,0.00305882352941176,0.00305882352941176,0.00305882352941176]],"container":"<table class=\"display\">\n  <thead>\n    <tr>\n      <th> <\/th>\n      <th>Axis.1<\/th>\n      <th>Axis.2<\/th>\n      <th>pval<\/th>\n      <th>r<\/th>\n      <th>Kingdom<\/th>\n      <th>Phylum<\/th>\n      <th>Class<\/th>\n      <th>Order<\/th>\n      <th>tax_rank_plot<\/th>\n      <th>Genus<\/th>\n      <th>Species<\/th>\n      <th>Strain<\/th>\n      <th>pval.adj<\/th>\n    <\/tr>\n  <\/thead>\n<\/table>","options":{"columnDefs":[{"className":"dt-right","targets":[1,2,3,4,13]},{"orderable":false,"targets":0}],"order":[],"autoWidth":false,"orderClasses":false}},"evals":[],"jsHooks":[]}</script>
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
<div id="htmlwidget-8e1fe96c0688a849169a" style="width:100%;height:auto;" class="datatables html-widget"></div>
<script type="application/json" data-for="htmlwidget-8e1fe96c0688a849169a">{"x":{"filter":"none","vertical":false,"data":[["1","2","3","4","5","6"],["Succinat_mM","Formiat_mM","Acetat_mM","Propionat_mM","Butyrat_mM","Total_SCFA_mM"],[0.163225559850326,0.44335408157543,-0.798708595679286,0.046761753008919,-0.750057796305815,-0.701255753014655],[0.510083629137183,0.0108062428054228,-0.466455045038887,-0.380510559986596,-0.524087796849334,-0.523321529938541],["Succinat_mM","Formiat_mM","Acetat_mM","Propionat_mM","Butyrat_mM","Total_SCFA_mM"],[0.001,0.005,0.001,0.026,0.001,0.001],[0.286827892102212,0.196679616533163,0.855515729854208,0.14697494780578,0.837254716605524,0.765625054833367]],"container":"<table class=\"display\">\n  <thead>\n    <tr>\n      <th> <\/th>\n      <th>label<\/th>\n      <th>Axis.1<\/th>\n      <th>Axis.2<\/th>\n      <th>id<\/th>\n      <th>pval<\/th>\n      <th>r<\/th>\n    <\/tr>\n  <\/thead>\n<\/table>","options":{"columnDefs":[{"className":"dt-right","targets":[2,3,5,6]},{"orderable":false,"targets":0}],"order":[],"autoWidth":false,"orderClasses":false}},"evals":[],"jsHooks":[]}</script>
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
<script type="application/json" data-for="htmlwidget-bb3486349adfa5d821dc">{"x":{"filter":"none","vertical":false,"data":[["ASV0012","ASV0051","ASV0010","ASV0015","ASV0106","ASV0124","ASV0136","ASV0202"],[-0.855235232721098,0.643808371532634,-0.858869634025459,0.746869100286677,-0.74562321070153,-0.858290482973094,-0.765802345012549,-0.663062220104871],[-0.0541313992033875,0.296923192688876,-0.375818735955723,0.297037592051541,0.361776052164635,-0.071384652306359,-0.347505820239148,-0.283546138710583],[0.001,0.001,0.001,0.001,0.001,0.001,0.001,0.001],[0.734357511667227,0.502652601612058,0.878896770546384,0.646044784054807,0.686835884256687,0.741758321747086,0.707213526726803,0.520049920508081],["Bacteria","Bacteria","Bacteria","Bacteria","Bacteria","Bacteria","Bacteria","Bacteria"],["Firmicutes","Firmicutes","Firmicutes","Firmicutes","Firmicutes","Firmicutes","Actinobacteriota","Firmicutes"],["Clostridia","Bacilli","Clostridia","Negativicutes","Clostridia","Clostridia","Coriobacteriia","Clostridia"],["Oscillospirales","Lactobacillales","Lachnospirales","Veillonellales-Selenomonadales","Oscillospirales","Monoglobales","Coriobacteriales","Oscillospirales"],["Ruminococcaceae","Lactobacillaceae","Lachnospiraceae","Veillonellaceae","Oscillospiraceae","Monoglobaceae","Eggerthellaceae","Butyricicoccaceae"],[null,null,null,null,null,null,null,null],[null,null,null,null,null,null,null,null],[null,null,null,null,null,null,null,null],[0.00288888888888889,0.00288888888888889,0.00288888888888889,0.00288888888888889,0.00288888888888889,0.00288888888888889,0.00288888888888889,0.00288888888888889]],"container":"<table class=\"display\">\n  <thead>\n    <tr>\n      <th> <\/th>\n      <th>Axis.1<\/th>\n      <th>Axis.2<\/th>\n      <th>pval<\/th>\n      <th>r<\/th>\n      <th>Kingdom<\/th>\n      <th>Phylum<\/th>\n      <th>Class<\/th>\n      <th>Order<\/th>\n      <th>tax_rank_plot<\/th>\n      <th>Genus<\/th>\n      <th>Species<\/th>\n      <th>Strain<\/th>\n      <th>pval.adj<\/th>\n    <\/tr>\n  <\/thead>\n<\/table>","options":{"columnDefs":[{"className":"dt-right","targets":[1,2,3,4,13]},{"orderable":false,"targets":0}],"order":[],"autoWidth":false,"orderClasses":false}},"evals":[],"jsHooks":[]}</script>
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
<script type="application/json" data-for="htmlwidget-44474b3e8cd160f85a98">{"x":{"filter":"none","vertical":false,"data":[["1","2","3","4","5"],["Succinat_mM","Formiat_mM","Acetat_mM","Butyrat_mM","Total_SCFA_mM"],[0.25387911017642,0.397146069924643,-0.931986035502789,-0.888849712418429,-0.873816422886923],[0.299944582034044,-0.283443016019585,-0.149688627900526,-0.290348308376891,-0.025036888976522],["Succinat_mM","Formiat_mM","Acetat_mM","Butyrat_mM","Total_SCFA_mM"],[0.003,0.001,0.001,0.001,0.001],[0.154421354875548,0.238064944186868,0.891004655694948,0.874355951443646,0.76418198671652]],"container":"<table class=\"display\">\n  <thead>\n    <tr>\n      <th> <\/th>\n      <th>label<\/th>\n      <th>Axis.1<\/th>\n      <th>Axis.2<\/th>\n      <th>id<\/th>\n      <th>pval<\/th>\n      <th>r<\/th>\n    <\/tr>\n  <\/thead>\n<\/table>","options":{"columnDefs":[{"className":"dt-right","targets":[2,3,5,6]},{"orderable":false,"targets":0}],"order":[],"autoWidth":false,"orderClasses":false}},"evals":[],"jsHooks":[]}</script>
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
<script type="application/json" data-for="htmlwidget-2514ce28ae593bde0f73">{"x":{"filter":"none","vertical":false,"data":[["ASV0012","ASV0051","ASV0010","ASV0015","ASV0106","ASV0124","ASV0136","ASV0202"],[-0.807869018332967,0.501736684674981,-0.789882571827549,0.681190472422282,-0.77695235901902,-0.755097137780456,-0.742905154586156,-0.601896083110978],[-0.241788609803216,0.624109179553861,-0.429380612178577,0.339414799731832,-0.132971234268199,-0.353917429832527,-0.352571634200053,-0.35347163418689],[0.001,0.001,0.001,0.001,0.001,0.001,0.001,0.001],[0.711114082612844,0.641251968752036,0.808282187391753,0.579222865995891,0.621336317328028,0.695429234623499,0.676214825953176,0.487221091039088],["Bacteria","Bacteria","Bacteria","Bacteria","Bacteria","Bacteria","Bacteria","Bacteria"],["Firmicutes","Firmicutes","Firmicutes","Firmicutes","Firmicutes","Firmicutes","Actinobacteriota","Firmicutes"],["Clostridia","Bacilli","Clostridia","Negativicutes","Clostridia","Clostridia","Coriobacteriia","Clostridia"],["Oscillospirales","Lactobacillales","Lachnospirales","Veillonellales-Selenomonadales","Oscillospirales","Monoglobales","Coriobacteriales","Oscillospirales"],["Ruminococcaceae","Lactobacillaceae","Lachnospiraceae","Veillonellaceae","Oscillospiraceae","Monoglobaceae","Eggerthellaceae","Butyricicoccaceae"],[null,null,null,null,null,null,null,null],[null,null,null,null,null,null,null,null],[null,null,null,null,null,null,null,null],[0.00305882352941176,0.00305882352941176,0.00305882352941176,0.00305882352941176,0.00305882352941176,0.00305882352941176,0.00305882352941176,0.00305882352941176]],"container":"<table class=\"display\">\n  <thead>\n    <tr>\n      <th> <\/th>\n      <th>Axis.1<\/th>\n      <th>Axis.2<\/th>\n      <th>pval<\/th>\n      <th>r<\/th>\n      <th>Kingdom<\/th>\n      <th>Phylum<\/th>\n      <th>Class<\/th>\n      <th>Order<\/th>\n      <th>tax_rank_plot<\/th>\n      <th>Genus<\/th>\n      <th>Species<\/th>\n      <th>Strain<\/th>\n      <th>pval.adj<\/th>\n    <\/tr>\n  <\/thead>\n<\/table>","options":{"columnDefs":[{"className":"dt-right","targets":[1,2,3,4,13]},{"orderable":false,"targets":0}],"order":[],"autoWidth":false,"orderClasses":false}},"evals":[],"jsHooks":[]}</script>
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
<script type="application/json" data-for="htmlwidget-5d05eac75555269eefd0">{"x":{"filter":"none","vertical":false,"data":[["1","2","3","4","5","6"],["Succinat_mM","Formiat_mM","Acetat_mM","Propionat_mM","Butyrat_mM","Total_SCFA_mM"],[0.163225559850326,0.44335408157543,-0.798708595679286,0.046761753008919,-0.750057796305815,-0.701255753014655],[0.510083629137183,0.0108062428054228,-0.466455045038887,-0.380510559986596,-0.524087796849334,-0.523321529938541],["Succinat_mM","Formiat_mM","Acetat_mM","Propionat_mM","Butyrat_mM","Total_SCFA_mM"],[0.001,0.005,0.001,0.026,0.001,0.001],[0.286827892102212,0.196679616533163,0.855515729854208,0.14697494780578,0.837254716605524,0.765625054833367]],"container":"<table class=\"display\">\n  <thead>\n    <tr>\n      <th> <\/th>\n      <th>label<\/th>\n      <th>Axis.1<\/th>\n      <th>Axis.2<\/th>\n      <th>id<\/th>\n      <th>pval<\/th>\n      <th>r<\/th>\n    <\/tr>\n  <\/thead>\n<\/table>","options":{"columnDefs":[{"className":"dt-right","targets":[2,3,5,6]},{"orderable":false,"targets":0}],"order":[],"autoWidth":false,"orderClasses":false}},"evals":[],"jsHooks":[]}</script>
```


```r
env_out$plot %>% 
  export::graph2ppt(append = TRUE,
                    width = 317.48031496 * 1,
                    height = 0.618 * 317.48031496 * 1 , paper = "A4",  scaling = 2,
                    file = out_pptx)
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
                     m = "NMDS",
                     seed = 123,
                     axis1 = 1,
                     axis2 = 2) -> NMDS_tmp
```

```
## [1] "bjaccard"
## Run 0 stress 0.120039 
## Run 1 stress 0.1245844 
## Run 2 stress 0.127509 
## Run 3 stress 0.1211789 
## Run 4 stress 0.122176 
## Run 5 stress 0.1240885 
## Run 6 stress 0.1267059 
## Run 7 stress 0.1275117 
## Run 8 stress 0.1293284 
## Run 9 stress 0.1218598 
## Run 10 stress 0.1273317 
## Run 11 stress 0.1205263 
## ... Procrustes: rmse 0.007012382  max resid 0.03494895 
## Run 12 stress 0.1275152 
## Run 13 stress 0.1218598 
## Run 14 stress 0.1245172 
## Run 15 stress 0.1267052 
## Run 16 stress 0.1249771 
## Run 17 stress 0.1218599 
## Run 18 stress 0.1248583 
## Run 19 stress 0.1224677 
## Run 20 stress 0.1360082 
## *** No convergence -- monoMDS stopping criteria:
##      3: no. of iterations >= maxit
##     17: stress ratio > sratmax
## [1] "aichinson"
## Run 0 stress 0.1270731 
## Run 1 stress 0.1269486 
## ... New best solution
## ... Procrustes: rmse 0.009705507  max resid 0.04061402 
## Run 2 stress 0.1240857 
## ... New best solution
## ... Procrustes: rmse 0.01735242  max resid 0.1082986 
## Run 3 stress 0.1359845 
## Run 4 stress 0.1377861 
## Run 5 stress 0.1290929 
## Run 6 stress 0.1339125 
## Run 7 stress 0.1488893 
## Run 8 stress 0.1408569 
## Run 9 stress 0.1334008 
## Run 10 stress 0.142561 
## Run 11 stress 0.1271091 
## Run 12 stress 0.1240863 
## ... Procrustes: rmse 0.0003969696  max resid 0.00193396 
## ... Similar to previous best
## Run 13 stress 0.1327147 
## Run 14 stress 0.1292787 
## Run 15 stress 0.1397015 
## Run 16 stress 0.1460579 
## Run 17 stress 0.1290929 
## Run 18 stress 0.1271091 
## Run 19 stress 0.1318199 
## Run 20 stress 0.1327789 
## *** Solution reached
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
<script type="application/json" data-for="htmlwidget-7531f021030c44f57ef2">{"x":{"filter":"none","vertical":false,"data":[["ASV0012","ASV0004","ASV0010","ASV0052","ASV0035","ASV0059","ASV0106","ASV0275"],[0.442307825573982,-0.519309208744457,-0.676834624672654,0.677634569168548,-0.453119161836603,-0.19499941713391,-0.242480102784295,0.639251626972184],[0.340033389334666,0.284579633520563,-0.249595620797355,0.0359243444398897,0.338428097619233,0.520187070665338,-0.627954364535686,0.0521799319621771],[0.001,0.002,0.001,0.001,0.001,0.002,0.001,0.002],[0.311258918426404,0.350667622101492,0.52040308307699,0.460479167855679,0.319850552081678,0.30861936116995,0.4531232841857,0.411365387886161],["Bacteria","Bacteria","Bacteria","Bacteria","Bacteria","Bacteria","Bacteria","Bacteria"],["Firmicutes","Proteobacteria","Firmicutes","Firmicutes","Proteobacteria","Bacteroidota","Firmicutes","Firmicutes"],["Clostridia","Gammaproteobacteria","Clostridia","Negativicutes","Gammaproteobacteria","Bacteroidia","Clostridia","Bacilli"],["Oscillospirales","Enterobacterales","Lachnospirales","Veillonellales-Selenomonadales","Burkholderiales","Bacteroidales","Oscillospirales","Erysipelotrichales"],["Ruminococcaceae","Enterobacteriaceae","Lachnospiraceae","Veillonellaceae","Sutterellaceae","Tannerellaceae","Oscillospiraceae","Erysipelotrichaceae"],[null,null,null,null,null,null,null,null],[null,null,null,null,null,null,null,null],[null,null,null,null,null,null,null,null],[0.00742857142857143,0.00945454545454546,0.00742857142857143,0.00742857142857143,0.00742857142857143,0.00945454545454546,0.00742857142857143,0.00945454545454546]],"container":"<table class=\"display\">\n  <thead>\n    <tr>\n      <th> <\/th>\n      <th>Axis.1<\/th>\n      <th>Axis.2<\/th>\n      <th>pval<\/th>\n      <th>r<\/th>\n      <th>Kingdom<\/th>\n      <th>Phylum<\/th>\n      <th>Class<\/th>\n      <th>Order<\/th>\n      <th>tax_rank_plot<\/th>\n      <th>Genus<\/th>\n      <th>Species<\/th>\n      <th>Strain<\/th>\n      <th>pval.adj<\/th>\n    <\/tr>\n  <\/thead>\n<\/table>","options":{"columnDefs":[{"className":"dt-right","targets":[1,2,3,4,13]},{"orderable":false,"targets":0}],"order":[],"autoWidth":false,"orderClasses":false}},"evals":[],"jsHooks":[]}</script>
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
<script type="application/json" data-for="htmlwidget-983315589fed18d622a7">{"x":{"filter":"none","vertical":false,"data":[["1","2","3","4"],["Formiat_mM","Acetat_mM","Valerat_mM","Total_SCFA_mM"],[-0.617911717338178,-0.132730246527057,0.18463335524204,-0.203516944103091],[0.0218294640123607,-0.390109080265109,-0.674316418777611,-0.440029362342907],["Formiat_mM","Acetat_mM","Valerat_mM","Total_SCFA_mM"],[0.001,0.028,0.001,0.006],[0.382291415922884,0.169802412848422,0.488792108500995,0.235044986260966]],"container":"<table class=\"display\">\n  <thead>\n    <tr>\n      <th> <\/th>\n      <th>label<\/th>\n      <th>Axis.1<\/th>\n      <th>Axis.2<\/th>\n      <th>id<\/th>\n      <th>pval<\/th>\n      <th>r<\/th>\n    <\/tr>\n  <\/thead>\n<\/table>","options":{"columnDefs":[{"className":"dt-right","targets":[2,3,5,6]},{"orderable":false,"targets":0}],"order":[],"autoWidth":false,"orderClasses":false}},"evals":[],"jsHooks":[]}</script>
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
<script type="application/json" data-for="htmlwidget-6f1d8b9e786d9d3fe28d">{"x":{"filter":"none","vertical":false,"data":[["ASV0012","ASV0004","ASV0010","ASV0052","ASV0030","ASV0059","ASV0106","ASV0275"],[0.399176079170203,-0.5394723061179,-0.668898551497996,0.651907797460639,0.0799020497582044,-0.215756175608905,-0.211350643520482,0.647451735256277],[0.441188769433854,0.204463800980237,-0.240983835400047,-0.136616644269266,0.528711246481583,0.49642087381465,-0.627395627072598,-0.0965347990902658],[0.001,0.002,0.001,0.001,0.008,0.002,0.001,0.001],[0.353989072456255,0.332835814979451,0.505498481120233,0.443647883881376,0.285919919711672,0.292984411272281,0.438294367386341,0.428512716921762],["Bacteria","Bacteria","Bacteria","Bacteria","Bacteria","Bacteria","Bacteria","Bacteria"],["Firmicutes","Proteobacteria","Firmicutes","Firmicutes","Bacteroidota","Bacteroidota","Firmicutes","Firmicutes"],["Clostridia","Gammaproteobacteria","Clostridia","Negativicutes","Bacteroidia","Bacteroidia","Clostridia","Bacilli"],["Oscillospirales","Enterobacterales","Lachnospirales","Veillonellales-Selenomonadales","Bacteroidales","Bacteroidales","Oscillospirales","Erysipelotrichales"],["Ruminococcaceae","Enterobacteriaceae","Lachnospiraceae","Veillonellaceae","Prevotellaceae","Tannerellaceae","Oscillospiraceae","Erysipelotrichaceae"],[null,null,null,null,null,null,null,null],[null,null,null,null,null,null,null,null],[null,null,null,null,null,null,null,null],[0.0104,0.013,0.0104,0.0104,0.0416,0.013,0.0104,0.0104]],"container":"<table class=\"display\">\n  <thead>\n    <tr>\n      <th> <\/th>\n      <th>Axis.1<\/th>\n      <th>Axis.2<\/th>\n      <th>pval<\/th>\n      <th>r<\/th>\n      <th>Kingdom<\/th>\n      <th>Phylum<\/th>\n      <th>Class<\/th>\n      <th>Order<\/th>\n      <th>tax_rank_plot<\/th>\n      <th>Genus<\/th>\n      <th>Species<\/th>\n      <th>Strain<\/th>\n      <th>pval.adj<\/th>\n    <\/tr>\n  <\/thead>\n<\/table>","options":{"columnDefs":[{"className":"dt-right","targets":[1,2,3,4,13]},{"orderable":false,"targets":0}],"order":[],"autoWidth":false,"orderClasses":false}},"evals":[],"jsHooks":[]}</script>
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
<div id="htmlwidget-dcdfc64296f5794e54a4" style="width:100%;height:auto;" class="datatables html-widget"></div>
<script type="application/json" data-for="htmlwidget-dcdfc64296f5794e54a4">{"x":{"filter":"none","vertical":false,"data":[["1","2","3","4","5","6"],["Succinat_mM","Formiat_mM","Acetat_mM","Propionat_mM","Butyrat_mM","Total_SCFA_mM"],[0.163225559850326,0.44335408157543,-0.798708595679286,0.046761753008919,-0.750057796305815,-0.701255753014655],[0.510083629137183,0.0108062428054228,-0.466455045038887,-0.380510559986596,-0.524087796849334,-0.523321529938541],["Succinat_mM","Formiat_mM","Acetat_mM","Propionat_mM","Butyrat_mM","Total_SCFA_mM"],[0.001,0.005,0.001,0.026,0.001,0.001],[0.286827892102212,0.196679616533163,0.855515729854208,0.14697494780578,0.837254716605524,0.765625054833367]],"container":"<table class=\"display\">\n  <thead>\n    <tr>\n      <th> <\/th>\n      <th>label<\/th>\n      <th>Axis.1<\/th>\n      <th>Axis.2<\/th>\n      <th>id<\/th>\n      <th>pval<\/th>\n      <th>r<\/th>\n    <\/tr>\n  <\/thead>\n<\/table>","options":{"columnDefs":[{"className":"dt-right","targets":[2,3,5,6]},{"orderable":false,"targets":0}],"order":[],"autoWidth":false,"orderClasses":false}},"evals":[],"jsHooks":[]}</script>
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
<script type="application/json" data-for="htmlwidget-873249b169eb6149d57e">{"x":{"filter":"none","vertical":false,"data":[["ASV0012","ASV0004","ASV0010","ASV0052","ASV0035","ASV0059","ASV0106","ASV0275"],[0.442307825573982,-0.519309208744457,-0.676834624672654,0.677634569168548,-0.453119161836603,-0.19499941713391,-0.242480102784295,0.639251626972184],[0.340033389334666,0.284579633520563,-0.249595620797355,0.0359243444398897,0.338428097619233,0.520187070665338,-0.627954364535686,0.0521799319621771],[0.001,0.002,0.001,0.001,0.001,0.002,0.001,0.002],[0.311258918426404,0.350667622101492,0.52040308307699,0.460479167855679,0.319850552081678,0.30861936116995,0.4531232841857,0.411365387886161],["Bacteria","Bacteria","Bacteria","Bacteria","Bacteria","Bacteria","Bacteria","Bacteria"],["Firmicutes","Proteobacteria","Firmicutes","Firmicutes","Proteobacteria","Bacteroidota","Firmicutes","Firmicutes"],["Clostridia","Gammaproteobacteria","Clostridia","Negativicutes","Gammaproteobacteria","Bacteroidia","Clostridia","Bacilli"],["Oscillospirales","Enterobacterales","Lachnospirales","Veillonellales-Selenomonadales","Burkholderiales","Bacteroidales","Oscillospirales","Erysipelotrichales"],["Ruminococcaceae","Enterobacteriaceae","Lachnospiraceae","Veillonellaceae","Sutterellaceae","Tannerellaceae","Oscillospiraceae","Erysipelotrichaceae"],[null,null,null,null,null,null,null,null],[null,null,null,null,null,null,null,null],[null,null,null,null,null,null,null,null],[0.00742857142857143,0.00945454545454546,0.00742857142857143,0.00742857142857143,0.00742857142857143,0.00945454545454546,0.00742857142857143,0.00945454545454546]],"container":"<table class=\"display\">\n  <thead>\n    <tr>\n      <th> <\/th>\n      <th>Axis.1<\/th>\n      <th>Axis.2<\/th>\n      <th>pval<\/th>\n      <th>r<\/th>\n      <th>Kingdom<\/th>\n      <th>Phylum<\/th>\n      <th>Class<\/th>\n      <th>Order<\/th>\n      <th>tax_rank_plot<\/th>\n      <th>Genus<\/th>\n      <th>Species<\/th>\n      <th>Strain<\/th>\n      <th>pval.adj<\/th>\n    <\/tr>\n  <\/thead>\n<\/table>","options":{"columnDefs":[{"className":"dt-right","targets":[1,2,3,4,13]},{"orderable":false,"targets":0}],"order":[],"autoWidth":false,"orderClasses":false}},"evals":[],"jsHooks":[]}</script>
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
<script type="application/json" data-for="htmlwidget-cbd3ad10bd2d35a2e76d">{"x":{"filter":"none","vertical":false,"data":[["1","2","3","4"],["Formiat_mM","Acetat_mM","Valerat_mM","Total_SCFA_mM"],[-0.617911717338178,-0.132730246527057,0.18463335524204,-0.203516944103091],[0.0218294640123607,-0.390109080265109,-0.674316418777611,-0.440029362342907],["Formiat_mM","Acetat_mM","Valerat_mM","Total_SCFA_mM"],[0.001,0.028,0.001,0.006],[0.382291415922884,0.169802412848422,0.488792108500995,0.235044986260966]],"container":"<table class=\"display\">\n  <thead>\n    <tr>\n      <th> <\/th>\n      <th>label<\/th>\n      <th>Axis.1<\/th>\n      <th>Axis.2<\/th>\n      <th>id<\/th>\n      <th>pval<\/th>\n      <th>r<\/th>\n    <\/tr>\n  <\/thead>\n<\/table>","options":{"columnDefs":[{"className":"dt-right","targets":[2,3,5,6]},{"orderable":false,"targets":0}],"order":[],"autoWidth":false,"orderClasses":false}},"evals":[],"jsHooks":[]}</script>
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
<script type="application/json" data-for="htmlwidget-f25d4e1770d66a86dcce">{"x":{"filter":"none","vertical":false,"data":[["ASV0012","ASV0004","ASV0010","ASV0052","ASV0030","ASV0059","ASV0106","ASV0275"],[0.399176079170203,-0.5394723061179,-0.668898551497996,0.651907797460639,0.0799020497582044,-0.215756175608905,-0.211350643520482,0.647451735256277],[0.441188769433854,0.204463800980237,-0.240983835400047,-0.136616644269266,0.528711246481583,0.49642087381465,-0.627395627072598,-0.0965347990902658],[0.001,0.002,0.001,0.001,0.008,0.002,0.001,0.001],[0.353989072456255,0.332835814979451,0.505498481120233,0.443647883881376,0.285919919711672,0.292984411272281,0.438294367386341,0.428512716921762],["Bacteria","Bacteria","Bacteria","Bacteria","Bacteria","Bacteria","Bacteria","Bacteria"],["Firmicutes","Proteobacteria","Firmicutes","Firmicutes","Bacteroidota","Bacteroidota","Firmicutes","Firmicutes"],["Clostridia","Gammaproteobacteria","Clostridia","Negativicutes","Bacteroidia","Bacteroidia","Clostridia","Bacilli"],["Oscillospirales","Enterobacterales","Lachnospirales","Veillonellales-Selenomonadales","Bacteroidales","Bacteroidales","Oscillospirales","Erysipelotrichales"],["Ruminococcaceae","Enterobacteriaceae","Lachnospiraceae","Veillonellaceae","Prevotellaceae","Tannerellaceae","Oscillospiraceae","Erysipelotrichaceae"],[null,null,null,null,null,null,null,null],[null,null,null,null,null,null,null,null],[null,null,null,null,null,null,null,null],[0.0104,0.013,0.0104,0.0104,0.0416,0.013,0.0104,0.0104]],"container":"<table class=\"display\">\n  <thead>\n    <tr>\n      <th> <\/th>\n      <th>Axis.1<\/th>\n      <th>Axis.2<\/th>\n      <th>pval<\/th>\n      <th>r<\/th>\n      <th>Kingdom<\/th>\n      <th>Phylum<\/th>\n      <th>Class<\/th>\n      <th>Order<\/th>\n      <th>tax_rank_plot<\/th>\n      <th>Genus<\/th>\n      <th>Species<\/th>\n      <th>Strain<\/th>\n      <th>pval.adj<\/th>\n    <\/tr>\n  <\/thead>\n<\/table>","options":{"columnDefs":[{"className":"dt-right","targets":[1,2,3,4,13]},{"orderable":false,"targets":0}],"order":[],"autoWidth":false,"orderClasses":false}},"evals":[],"jsHooks":[]}</script>
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
<script type="application/json" data-for="htmlwidget-969965edfe9f869aa904">{"x":{"filter":"none","vertical":false,"data":[["1","2","3"],["Formiat_mM","Valerat_mM","Total_SCFA_mM"],[-0.611305686967948,0.183828328100747,-0.17503752027219],[0.211616949856144,-0.612885588749617,-0.363524059901025],["Formiat_mM","Valerat_mM","Total_SCFA_mM"],[0.001,0.001,0.025],[0.418476376385773,0.40942159910928,0.162787875629961]],"container":"<table class=\"display\">\n  <thead>\n    <tr>\n      <th> <\/th>\n      <th>label<\/th>\n      <th>Axis.1<\/th>\n      <th>Axis.2<\/th>\n      <th>id<\/th>\n      <th>pval<\/th>\n      <th>r<\/th>\n    <\/tr>\n  <\/thead>\n<\/table>","options":{"columnDefs":[{"className":"dt-right","targets":[2,3,5,6]},{"orderable":false,"targets":0}],"order":[],"autoWidth":false,"orderClasses":false}},"evals":[],"jsHooks":[]}</script>
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
  subset_samples(Model2 == "Human2"  & Treatment != "DONOR" ) -> ps_tmp_forbdiv
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
## Run 0 stress 0.08288928 
## Run 1 stress 0.120357 
## Run 2 stress 0.1097329 
## Run 3 stress 0.1258396 
## Run 4 stress 0.09990649 
## Run 5 stress 0.09859401 
## Run 6 stress 0.1274844 
## Run 7 stress 0.08297845 
## ... Procrustes: rmse 0.004155955  max resid 0.04171612 
## Run 8 stress 0.110711 
## Run 9 stress 0.08288911 
## ... New best solution
## ... Procrustes: rmse 0.0001543536  max resid 0.001186569 
## ... Similar to previous best
## Run 10 stress 0.1254061 
## Run 11 stress 0.09564983 
## Run 12 stress 0.1164852 
## Run 13 stress 0.1233859 
## Run 14 stress 0.0874946 
## Run 15 stress 0.1106918 
## Run 16 stress 0.1042854 
## Run 17 stress 0.09569961 
## Run 18 stress 0.1255812 
## Run 19 stress 0.1076216 
## Run 20 stress 0.1284139 
## *** Solution reached
## [1] "aichinson"
## Run 0 stress 0.06214457 
## Run 1 stress 0.08035945 
## Run 2 stress 0.08833964 
## Run 3 stress 0.09020448 
## Run 4 stress 0.09023975 
## Run 5 stress 0.08052002 
## Run 6 stress 0.08820566 
## Run 7 stress 0.07107039 
## Run 8 stress 0.08167821 
## Run 9 stress 0.06388722 
## Run 10 stress 0.08016584 
## Run 11 stress 0.08816338 
## Run 12 stress 0.08731513 
## Run 13 stress 0.09190637 
## Run 14 stress 0.07144001 
## Run 15 stress 0.0853967 
## Run 16 stress 0.07141984 
## Run 17 stress 0.09366463 
## Run 18 stress 0.090434 
## Run 19 stress 0.09153035 
## Run 20 stress 0.08127396 
## *** No convergence -- monoMDS stopping criteria:
##     16: stress ratio > sratmax
##      4: scale factor of the gradient < sfgrmin
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

![](NRP72_draft_16S_Figs_beta_files/figure-html/unnamed-chunk-62-1.png)<!-- -->


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
<div id="htmlwidget-450f171f14364aa3ebed" style="width:100%;height:auto;" class="datatables html-widget"></div>
<script type="application/json" data-for="htmlwidget-450f171f14364aa3ebed">{"x":{"filter":"none","vertical":false,"data":[["ASV0002","ASV0003","ASV0020","ASV0040","ASV0032","ASV0085","ASV0124","ASV0129"],[0.531409137016163,0.921845009378856,-0.368055139001907,-0.120260120844456,-0.392484178712118,0.61253499119304,0.764992599341263,0.697037134347215],[-0.61362441466168,-0.134156447709648,0.821140290473168,0.687701880255258,0.608216686567871,-0.488291094989382,0.242058610351244,0.306202113223079],[0.001,0.001,0.001,0.001,0.001,0.001,0.001,0.001],[0.658930593173153,0.867796173778774,0.809735961984071,0.48739637277214,0.523971368358926,0.613627308881787,0.643806047892077,0.579620500801257],["Bacteria","Bacteria","Bacteria","Bacteria","Bacteria","Bacteria","Bacteria","Bacteria"],["Bacteroidota","Firmicutes","Firmicutes","Firmicutes","Proteobacteria","Firmicutes","Firmicutes","Actinobacteriota"],["Bacteroidia","Clostridia","Bacilli","Negativicutes","Gammaproteobacteria","Clostridia","Clostridia","Coriobacteriia"],["Bacteroidales","Oscillospirales","Lactobacillales","Veillonellales-Selenomonadales","Burkholderiales","Oscillospirales","Monoglobales","Coriobacteriales"],["Bacteroidaceae","Ruminococcaceae","Lactobacillaceae","Veillonellaceae","Sutterellaceae","Oscillospiraceae","Monoglobaceae","Coriobacteriaceae"],[null,null,null,null,null,null,null,null],[null,null,null,null,null,null,null,null],[null,null,null,null,null,null,null,null],[0.00247619047619048,0.00247619047619048,0.00247619047619048,0.00247619047619048,0.00247619047619048,0.00247619047619048,0.00247619047619048,0.00247619047619048]],"container":"<table class=\"display\">\n  <thead>\n    <tr>\n      <th> <\/th>\n      <th>Axis.1<\/th>\n      <th>Axis.2<\/th>\n      <th>pval<\/th>\n      <th>r<\/th>\n      <th>Kingdom<\/th>\n      <th>Phylum<\/th>\n      <th>Class<\/th>\n      <th>Order<\/th>\n      <th>tax_rank_plot<\/th>\n      <th>Genus<\/th>\n      <th>Species<\/th>\n      <th>Strain<\/th>\n      <th>pval.adj<\/th>\n    <\/tr>\n  <\/thead>\n<\/table>","options":{"columnDefs":[{"className":"dt-right","targets":[1,2,3,4,13]},{"orderable":false,"targets":0}],"order":[],"autoWidth":false,"orderClasses":false}},"evals":[],"jsHooks":[]}</script>
```


```r
out$plot -> NMDS_tmp2_fam

NMDS_tmp2 + 
  theme(legend.position = "none") -> NMDS_tmp2_noleg

NMDS_tmp2_noleg
```

![](NRP72_draft_16S_Figs_beta_files/figure-html/unnamed-chunk-64-1.png)<!-- -->

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
<div id="htmlwidget-fbfc52d82635d988c6e1" style="width:100%;height:auto;" class="datatables html-widget"></div>
<script type="application/json" data-for="htmlwidget-fbfc52d82635d988c6e1">{"x":{"filter":"none","vertical":false,"data":[["1","2","3","4","5","6"],["Formiat_mM","Acetat_mM","Propionat_mM","Butyrat_mM","Valerat_mM","Total_SCFA_mM"],[-0.556240323944997,0.918017427039651,0.883498303138054,0.94906817301198,0.889146042330164,0.915982815115401],[0.446155245113644,-0.214805325379093,-0.0426450877578052,-0.168070951957279,0.277774905154955,-0.300935504556202],["Formiat_mM","Acetat_mM","Propionat_mM","Butyrat_mM","Valerat_mM","Total_SCFA_mM"],[0.001,0.001,0.001,0.001,0.001,0.001],[0.508457800724851,0.88889732415972,0.782387855157692,0.928978241916124,0.867739582525238,0.929586695489231]],"container":"<table class=\"display\">\n  <thead>\n    <tr>\n      <th> <\/th>\n      <th>label<\/th>\n      <th>Axis.1<\/th>\n      <th>Axis.2<\/th>\n      <th>id<\/th>\n      <th>pval<\/th>\n      <th>r<\/th>\n    <\/tr>\n  <\/thead>\n<\/table>","options":{"columnDefs":[{"className":"dt-right","targets":[2,3,5,6]},{"orderable":false,"targets":0}],"order":[],"autoWidth":false,"orderClasses":false}},"evals":[],"jsHooks":[]}</script>
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

![](NRP72_draft_16S_Figs_beta_files/figure-html/unnamed-chunk-67-1.png)<!-- -->


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
<div id="htmlwidget-2f4c52b7f36e26c300a7" style="width:100%;height:auto;" class="datatables html-widget"></div>
<script type="application/json" data-for="htmlwidget-2f4c52b7f36e26c300a7">{"x":{"filter":"none","vertical":false,"data":[["ASV0002","ASV0003","ASV0020","ASV0040","ASV0032","ASV0085","ASV0124","ASV0129"],[0.658535942112849,0.918863525282554,-0.416527516477986,-0.147892324746833,-0.44488880930227,0.739256904917888,0.744758174736518,0.673238877455571],[-0.410637601305111,0.0372622449244126,0.795480730813088,0.766385155585566,0.595346130983358,-0.205093111790235,0.213534339998062,0.300175617660472],[0.001,0.001,0.001,0.001,0.001,0.001,0.001,0.001],[0.602292826660072,0.84569865299149,0.806284765078243,0.609218346420935,0.552363068319245,0.588563955972578,0.600261653195277,0.543355987555484],["Bacteria","Bacteria","Bacteria","Bacteria","Bacteria","Bacteria","Bacteria","Bacteria"],["Bacteroidota","Firmicutes","Firmicutes","Firmicutes","Proteobacteria","Firmicutes","Firmicutes","Actinobacteriota"],["Bacteroidia","Clostridia","Bacilli","Negativicutes","Gammaproteobacteria","Clostridia","Clostridia","Coriobacteriia"],["Bacteroidales","Oscillospirales","Lactobacillales","Veillonellales-Selenomonadales","Burkholderiales","Oscillospirales","Monoglobales","Coriobacteriales"],["Bacteroidaceae","Ruminococcaceae","Lactobacillaceae","Veillonellaceae","Sutterellaceae","Oscillospiraceae","Monoglobaceae","Coriobacteriaceae"],[null,null,null,null,null,null,null,null],[null,null,null,null,null,null,null,null],[null,null,null,null,null,null,null,null],[0.0026,0.0026,0.0026,0.0026,0.0026,0.0026,0.0026,0.0026]],"container":"<table class=\"display\">\n  <thead>\n    <tr>\n      <th> <\/th>\n      <th>Axis.1<\/th>\n      <th>Axis.2<\/th>\n      <th>pval<\/th>\n      <th>r<\/th>\n      <th>Kingdom<\/th>\n      <th>Phylum<\/th>\n      <th>Class<\/th>\n      <th>Order<\/th>\n      <th>tax_rank_plot<\/th>\n      <th>Genus<\/th>\n      <th>Species<\/th>\n      <th>Strain<\/th>\n      <th>pval.adj<\/th>\n    <\/tr>\n  <\/thead>\n<\/table>","options":{"columnDefs":[{"className":"dt-right","targets":[1,2,3,4,13]},{"orderable":false,"targets":0}],"order":[],"autoWidth":false,"orderClasses":false}},"evals":[],"jsHooks":[]}</script>
```


```r
out$plot -> NMDS_tmp2_fam

NMDS_tmp2 + 
  theme(legend.position = "none") -> NMDS_tmp2_noleg

NMDS_tmp2_noleg
```

![](NRP72_draft_16S_Figs_beta_files/figure-html/unnamed-chunk-69-1.png)<!-- -->

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
<div id="htmlwidget-1d5df3d80f5735caaf13" style="width:100%;height:auto;" class="datatables html-widget"></div>
<script type="application/json" data-for="htmlwidget-1d5df3d80f5735caaf13">{"x":{"filter":"none","vertical":false,"data":[["1","2","3","4","5","6"],["Succinat_mM","Formiat_mM","Acetat_mM","Propionat_mM","Butyrat_mM","Total_SCFA_mM"],[0.163225559850326,0.44335408157543,-0.798708595679286,0.046761753008919,-0.750057796305815,-0.701255753014655],[0.510083629137183,0.0108062428054228,-0.466455045038887,-0.380510559986596,-0.524087796849334,-0.523321529938541],["Succinat_mM","Formiat_mM","Acetat_mM","Propionat_mM","Butyrat_mM","Total_SCFA_mM"],[0.001,0.005,0.001,0.026,0.001,0.001],[0.286827892102212,0.196679616533163,0.855515729854208,0.14697494780578,0.837254716605524,0.765625054833367]],"container":"<table class=\"display\">\n  <thead>\n    <tr>\n      <th> <\/th>\n      <th>label<\/th>\n      <th>Axis.1<\/th>\n      <th>Axis.2<\/th>\n      <th>id<\/th>\n      <th>pval<\/th>\n      <th>r<\/th>\n    <\/tr>\n  <\/thead>\n<\/table>","options":{"columnDefs":[{"className":"dt-right","targets":[2,3,5,6]},{"orderable":false,"targets":0}],"order":[],"autoWidth":false,"orderClasses":false}},"evals":[],"jsHooks":[]}</script>
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

![](NRP72_draft_16S_Figs_beta_files/figure-html/unnamed-chunk-73-1.png)<!-- -->


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
<div id="htmlwidget-67c25ca3336b15b0293a" style="width:100%;height:auto;" class="datatables html-widget"></div>
<script type="application/json" data-for="htmlwidget-67c25ca3336b15b0293a">{"x":{"filter":"none","vertical":false,"data":[["ASV0002","ASV0003","ASV0020","ASV0040","ASV0032","ASV0085","ASV0124","ASV0129"],[0.531409137016163,0.921845009378856,-0.368055139001907,-0.120260120844456,-0.392484178712118,0.61253499119304,0.764992599341263,0.697037134347215],[-0.61362441466168,-0.134156447709648,0.821140290473168,0.687701880255258,0.608216686567871,-0.488291094989382,0.242058610351244,0.306202113223079],[0.001,0.001,0.001,0.001,0.001,0.001,0.001,0.001],[0.658930593173153,0.867796173778774,0.809735961984071,0.48739637277214,0.523971368358926,0.613627308881787,0.643806047892077,0.579620500801257],["Bacteria","Bacteria","Bacteria","Bacteria","Bacteria","Bacteria","Bacteria","Bacteria"],["Bacteroidota","Firmicutes","Firmicutes","Firmicutes","Proteobacteria","Firmicutes","Firmicutes","Actinobacteriota"],["Bacteroidia","Clostridia","Bacilli","Negativicutes","Gammaproteobacteria","Clostridia","Clostridia","Coriobacteriia"],["Bacteroidales","Oscillospirales","Lactobacillales","Veillonellales-Selenomonadales","Burkholderiales","Oscillospirales","Monoglobales","Coriobacteriales"],["Bacteroidaceae","Ruminococcaceae","Lactobacillaceae","Veillonellaceae","Sutterellaceae","Oscillospiraceae","Monoglobaceae","Coriobacteriaceae"],[null,null,null,null,null,null,null,null],[null,null,null,null,null,null,null,null],[null,null,null,null,null,null,null,null],[0.00247619047619048,0.00247619047619048,0.00247619047619048,0.00247619047619048,0.00247619047619048,0.00247619047619048,0.00247619047619048,0.00247619047619048]],"container":"<table class=\"display\">\n  <thead>\n    <tr>\n      <th> <\/th>\n      <th>Axis.1<\/th>\n      <th>Axis.2<\/th>\n      <th>pval<\/th>\n      <th>r<\/th>\n      <th>Kingdom<\/th>\n      <th>Phylum<\/th>\n      <th>Class<\/th>\n      <th>Order<\/th>\n      <th>tax_rank_plot<\/th>\n      <th>Genus<\/th>\n      <th>Species<\/th>\n      <th>Strain<\/th>\n      <th>pval.adj<\/th>\n    <\/tr>\n  <\/thead>\n<\/table>","options":{"columnDefs":[{"className":"dt-right","targets":[1,2,3,4,13]},{"orderable":false,"targets":0}],"order":[],"autoWidth":false,"orderClasses":false}},"evals":[],"jsHooks":[]}</script>
```


```r
out$plot -> PCoA_tmp2_fam

PCoA_tmp2 + 
  theme(legend.position = "none") -> PCoA_tmp2_noleg

PCoA_tmp2_noleg
```

![](NRP72_draft_16S_Figs_beta_files/figure-html/unnamed-chunk-75-1.png)<!-- -->

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
<div id="htmlwidget-530a78b6872d546cc16a" style="width:100%;height:auto;" class="datatables html-widget"></div>
<script type="application/json" data-for="htmlwidget-530a78b6872d546cc16a">{"x":{"filter":"none","vertical":false,"data":[["1","2","3","4","5","6"],["Formiat_mM","Acetat_mM","Propionat_mM","Butyrat_mM","Valerat_mM","Total_SCFA_mM"],[-0.556240323944997,0.918017427039651,0.883498303138054,0.94906817301198,0.889146042330164,0.915982815115401],[0.446155245113644,-0.214805325379093,-0.0426450877578052,-0.168070951957279,0.277774905154955,-0.300935504556202],["Formiat_mM","Acetat_mM","Propionat_mM","Butyrat_mM","Valerat_mM","Total_SCFA_mM"],[0.001,0.001,0.001,0.001,0.001,0.001],[0.508457800724851,0.88889732415972,0.782387855157692,0.928978241916124,0.867739582525238,0.929586695489231]],"container":"<table class=\"display\">\n  <thead>\n    <tr>\n      <th> <\/th>\n      <th>label<\/th>\n      <th>Axis.1<\/th>\n      <th>Axis.2<\/th>\n      <th>id<\/th>\n      <th>pval<\/th>\n      <th>r<\/th>\n    <\/tr>\n  <\/thead>\n<\/table>","options":{"columnDefs":[{"className":"dt-right","targets":[2,3,5,6]},{"orderable":false,"targets":0}],"order":[],"autoWidth":false,"orderClasses":false}},"evals":[],"jsHooks":[]}</script>
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

![](NRP72_draft_16S_Figs_beta_files/figure-html/unnamed-chunk-78-1.png)<!-- -->


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
<div id="htmlwidget-4cadc5f8c7ab9ff923ac" style="width:100%;height:auto;" class="datatables html-widget"></div>
<script type="application/json" data-for="htmlwidget-4cadc5f8c7ab9ff923ac">{"x":{"filter":"none","vertical":false,"data":[["ASV0002","ASV0003","ASV0020","ASV0040","ASV0032","ASV0085","ASV0124","ASV0129"],[0.658535942112849,0.918863525282554,-0.416527516477986,-0.147892324746833,-0.44488880930227,0.739256904917888,0.744758174736518,0.673238877455571],[-0.410637601305111,0.0372622449244126,0.795480730813088,0.766385155585566,0.595346130983358,-0.205093111790235,0.213534339998062,0.300175617660472],[0.001,0.001,0.001,0.001,0.001,0.001,0.001,0.001],[0.602292826660072,0.84569865299149,0.806284765078243,0.609218346420935,0.552363068319245,0.588563955972578,0.600261653195277,0.543355987555484],["Bacteria","Bacteria","Bacteria","Bacteria","Bacteria","Bacteria","Bacteria","Bacteria"],["Bacteroidota","Firmicutes","Firmicutes","Firmicutes","Proteobacteria","Firmicutes","Firmicutes","Actinobacteriota"],["Bacteroidia","Clostridia","Bacilli","Negativicutes","Gammaproteobacteria","Clostridia","Clostridia","Coriobacteriia"],["Bacteroidales","Oscillospirales","Lactobacillales","Veillonellales-Selenomonadales","Burkholderiales","Oscillospirales","Monoglobales","Coriobacteriales"],["Bacteroidaceae","Ruminococcaceae","Lactobacillaceae","Veillonellaceae","Sutterellaceae","Oscillospiraceae","Monoglobaceae","Coriobacteriaceae"],[null,null,null,null,null,null,null,null],[null,null,null,null,null,null,null,null],[null,null,null,null,null,null,null,null],[0.0026,0.0026,0.0026,0.0026,0.0026,0.0026,0.0026,0.0026]],"container":"<table class=\"display\">\n  <thead>\n    <tr>\n      <th> <\/th>\n      <th>Axis.1<\/th>\n      <th>Axis.2<\/th>\n      <th>pval<\/th>\n      <th>r<\/th>\n      <th>Kingdom<\/th>\n      <th>Phylum<\/th>\n      <th>Class<\/th>\n      <th>Order<\/th>\n      <th>tax_rank_plot<\/th>\n      <th>Genus<\/th>\n      <th>Species<\/th>\n      <th>Strain<\/th>\n      <th>pval.adj<\/th>\n    <\/tr>\n  <\/thead>\n<\/table>","options":{"columnDefs":[{"className":"dt-right","targets":[1,2,3,4,13]},{"orderable":false,"targets":0}],"order":[],"autoWidth":false,"orderClasses":false}},"evals":[],"jsHooks":[]}</script>
```


```r
out$plot -> PCoA_tmp2_fam

PCoA_tmp2 + 
  theme(legend.position = "none") -> PCoA_tmp2_noleg

PCoA_tmp2_noleg
```

![](NRP72_draft_16S_Figs_beta_files/figure-html/unnamed-chunk-80-1.png)<!-- -->

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
<div id="htmlwidget-fc8df40fd67c2027a709" style="width:100%;height:auto;" class="datatables html-widget"></div>
<script type="application/json" data-for="htmlwidget-fc8df40fd67c2027a709">{"x":{"filter":"none","vertical":false,"data":[["1","2","3","4","5","6"],["Formiat_mM","Acetat_mM","Propionat_mM","Butyrat_mM","Valerat_mM","Total_SCFA_mM"],[-0.61051198424088,0.933184996490246,0.858429993399151,0.95197411200031,0.862845615750082,0.955688473791752],[0.308703406703288,-0.0139641498536305,0.213617319877914,0.00521060611628454,0.28290963535989,-0.0820419379972459],["Formiat_mM","Acetat_mM","Propionat_mM","Butyrat_mM","Valerat_mM","Total_SCFA_mM"],[0.001,0.001,0.001,0.001,0.001,0.001],[0.468022676211952,0.871029235155635,0.782534412919089,0.906281860334878,0.824540418398603,0.920071338528752]],"container":"<table class=\"display\">\n  <thead>\n    <tr>\n      <th> <\/th>\n      <th>label<\/th>\n      <th>Axis.1<\/th>\n      <th>Axis.2<\/th>\n      <th>id<\/th>\n      <th>pval<\/th>\n      <th>r<\/th>\n    <\/tr>\n  <\/thead>\n<\/table>","options":{"columnDefs":[{"className":"dt-right","targets":[2,3,5,6]},{"orderable":false,"targets":0}],"order":[],"autoWidth":false,"orderClasses":false}},"evals":[],"jsHooks":[]}</script>
```


```r
env_out$plot %>% 
  export::graph2ppt(append = TRUE,
                    width = 317.48031496 * 1,
                    height = 0.618 * 317.48031496 * 1 , paper = "A4",  scaling = 2,
                    file = out_pptx)
```













### Human1 - all:



```r
ps_rare %>%
  subset_samples(Model2 == "Human1"  & Treatment != "DONOR"  & Fermentation %in% c(1,2)) -> ps_tmp_forbdiv
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
## Run 0 stress 0.08884942 
## Run 1 stress 0.09519477 
## Run 2 stress 0.09776717 
## Run 3 stress 0.08884959 
## ... Procrustes: rmse 0.0001794142  max resid 0.001431027 
## ... Similar to previous best
## Run 4 stress 0.0951949 
## Run 5 stress 0.09976361 
## Run 6 stress 0.101409 
## Run 7 stress 0.09174832 
## Run 8 stress 0.09345993 
## Run 9 stress 0.1052206 
## Run 10 stress 0.09415681 
## Run 11 stress 0.09728743 
## Run 12 stress 0.08582309 
## ... New best solution
## ... Procrustes: rmse 0.0488013  max resid 0.4630264 
## Run 13 stress 0.08885021 
## Run 14 stress 0.09572543 
## Run 15 stress 0.1021217 
## Run 16 stress 0.09836237 
## Run 17 stress 0.08885707 
## Run 18 stress 0.08884962 
## Run 19 stress 0.09518921 
## Run 20 stress 0.09174832 
## *** No convergence -- monoMDS stopping criteria:
##      1: no. of iterations >= maxit
##     18: stress ratio > sratmax
##      1: scale factor of the gradient < sfgrmin
## [1] "aichinson"
## Run 0 stress 0.09630452 
## Run 1 stress 0.09630296 
## ... New best solution
## ... Procrustes: rmse 0.0007381658  max resid 0.005451184 
## ... Similar to previous best
## Run 2 stress 0.09993125 
## Run 3 stress 0.09630452 
## ... Procrustes: rmse 0.0007391388  max resid 0.005460598 
## ... Similar to previous best
## Run 4 stress 0.09630296 
## ... New best solution
## ... Procrustes: rmse 5.884283e-06  max resid 2.931944e-05 
## ... Similar to previous best
## Run 5 stress 0.099956 
## Run 6 stress 0.09993124 
## Run 7 stress 0.09630451 
## ... Procrustes: rmse 0.0007478364  max resid 0.005478858 
## ... Similar to previous best
## Run 8 stress 0.09630451 
## ... Procrustes: rmse 0.0007450523  max resid 0.005458312 
## ... Similar to previous best
## Run 9 stress 0.09993125 
## Run 10 stress 0.09993124 
## Run 11 stress 0.09993124 
## Run 12 stress 0.09630451 
## ... Procrustes: rmse 0.0007497923  max resid 0.005495547 
## ... Similar to previous best
## Run 13 stress 0.09993125 
## Run 14 stress 0.09993125 
## Run 15 stress 0.09995601 
## Run 16 stress 0.09993124 
## Run 17 stress 0.09993124 
## Run 18 stress 0.09993125 
## Run 19 stress 0.09630297 
## ... Procrustes: rmse 6.809404e-06  max resid 3.834038e-05 
## ... Similar to previous best
## Run 20 stress 0.09993125 
## *** Solution reached
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

![](NRP72_draft_16S_Figs_beta_files/figure-html/unnamed-chunk-85-1.png)<!-- -->


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
<div id="htmlwidget-d8c7d17cf4080686b93c" style="width:100%;height:auto;" class="datatables html-widget"></div>
<script type="application/json" data-for="htmlwidget-d8c7d17cf4080686b93c">{"x":{"filter":"none","vertical":false,"data":[["ASV0002","ASV0051","ASV0010","ASV0015","ASV0019","ASV0035","ASV0106","ASV0124"],[-0.596365968169365,0.65868508864788,-0.882275771735157,0.791983810973008,0.505766022416626,0.544200939838375,-0.742569357424228,-0.647095453887146],[-0.329816353820077,0.333036945606869,0.198219381930944,0.132828151930179,0.471256503357815,0.506200658362873,-0.0254828557713435,0.390896324977618],[0.001,0.001,0.001,0.001,0.001,0.001,0.001,0.001],[0.464431195237754,0.544779653146219,0.817701460763952,0.644881674788516,0.47788196138817,0.552393769447976,0.552058626523694,0.571532463322419],["Bacteria","Bacteria","Bacteria","Bacteria","Bacteria","Bacteria","Bacteria","Bacteria"],["Bacteroidota","Firmicutes","Firmicutes","Firmicutes","Proteobacteria","Proteobacteria","Firmicutes","Firmicutes"],["Bacteroidia","Bacilli","Clostridia","Negativicutes","Gammaproteobacteria","Gammaproteobacteria","Clostridia","Clostridia"],["Bacteroidales","Lactobacillales","Lachnospirales","Veillonellales-Selenomonadales","Pseudomonadales","Burkholderiales","Oscillospirales","Monoglobales"],["Bacteroidaceae","Lactobacillaceae","Lachnospiraceae","Veillonellaceae","Pseudomonadaceae","Sutterellaceae","Oscillospiraceae","Monoglobaceae"],[null,null,null,null,null,null,null,null],[null,null,null,null,null,null,null,null],[null,null,null,null,null,null,null,null],[0.0026,0.0026,0.0026,0.0026,0.0026,0.0026,0.0026,0.0026]],"container":"<table class=\"display\">\n  <thead>\n    <tr>\n      <th> <\/th>\n      <th>Axis.1<\/th>\n      <th>Axis.2<\/th>\n      <th>pval<\/th>\n      <th>r<\/th>\n      <th>Kingdom<\/th>\n      <th>Phylum<\/th>\n      <th>Class<\/th>\n      <th>Order<\/th>\n      <th>tax_rank_plot<\/th>\n      <th>Genus<\/th>\n      <th>Species<\/th>\n      <th>Strain<\/th>\n      <th>pval.adj<\/th>\n    <\/tr>\n  <\/thead>\n<\/table>","options":{"columnDefs":[{"className":"dt-right","targets":[1,2,3,4,13]},{"orderable":false,"targets":0}],"order":[],"autoWidth":false,"orderClasses":false}},"evals":[],"jsHooks":[]}</script>
```


```r
out$plot -> NMDS_tmp2_fam

NMDS_tmp2 + 
  theme(legend.position = "none") -> NMDS_tmp2_noleg

NMDS_tmp2_noleg
```

![](NRP72_draft_16S_Figs_beta_files/figure-html/unnamed-chunk-87-1.png)<!-- -->

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
<div id="htmlwidget-fbc9973e359619ea20b8" style="width:100%;height:auto;" class="datatables html-widget"></div>
<script type="application/json" data-for="htmlwidget-fbc9973e359619ea20b8">{"x":{"filter":"none","vertical":false,"data":[["1","2","3","4","5","6"],["Formiat_mM","Acetat_mM","Propionat_mM","Butyrat_mM","Valerat_mM","Total_SCFA_mM"],[0.27306840343581,-0.879625432657993,-0.0441586289042728,-0.795312448898601,-0.118568666569057,-0.780538747503882],[0.549238723654045,-0.238276371250953,-0.298723166856065,-0.45191431945377,-0.67858113286701,-0.431060339338],["Formiat_mM","Acetat_mM","Propionat_mM","Butyrat_mM","Valerat_mM","Total_SCFA_mM"],[0.001,0.001,0.01,0.001,0.001,0.001],[0.376229528516106,0.830516530875283,0.0911855149232218,0.836748443500454,0.474530882575039,0.795053752505121]],"container":"<table class=\"display\">\n  <thead>\n    <tr>\n      <th> <\/th>\n      <th>label<\/th>\n      <th>Axis.1<\/th>\n      <th>Axis.2<\/th>\n      <th>id<\/th>\n      <th>pval<\/th>\n      <th>r<\/th>\n    <\/tr>\n  <\/thead>\n<\/table>","options":{"columnDefs":[{"className":"dt-right","targets":[2,3,5,6]},{"orderable":false,"targets":0}],"order":[],"autoWidth":false,"orderClasses":false}},"evals":[],"jsHooks":[]}</script>
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

![](NRP72_draft_16S_Figs_beta_files/figure-html/unnamed-chunk-90-1.png)<!-- -->


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
<div id="htmlwidget-fc9c72093435cc6a0a18" style="width:100%;height:auto;" class="datatables html-widget"></div>
<script type="application/json" data-for="htmlwidget-fc9c72093435cc6a0a18">{"x":{"filter":"none","vertical":false,"data":[["ASV0002","ASV0012","ASV0051","ASV0010","ASV0015","ASV0035","ASV0106","ASV0124"],[-0.636725507108502,-0.432685232978908,0.698941189730734,-0.774433924337121,0.759890270630475,0.553845729547033,-0.696754170275816,-0.534362816472827],[-0.20387924927289,-0.505654315299078,0.267634778334487,0.409651810560177,-0.11892905216364,0.408241848362023,0.182516827457269,0.49997692210747],[0.001,0.001,0.001,0.001,0.001,0.001,0.001,0.001],[0.446986119686657,0.442902797418592,0.560147161276364,0.767562509059425,0.591577342847398,0.473406498891526,0.518778766101807,0.535520542268831],["Bacteria","Bacteria","Bacteria","Bacteria","Bacteria","Bacteria","Bacteria","Bacteria"],["Bacteroidota","Firmicutes","Firmicutes","Firmicutes","Firmicutes","Proteobacteria","Firmicutes","Firmicutes"],["Bacteroidia","Clostridia","Bacilli","Clostridia","Negativicutes","Gammaproteobacteria","Clostridia","Clostridia"],["Bacteroidales","Oscillospirales","Lactobacillales","Lachnospirales","Veillonellales-Selenomonadales","Burkholderiales","Oscillospirales","Monoglobales"],["Bacteroidaceae","Ruminococcaceae","Lactobacillaceae","Lachnospiraceae","Veillonellaceae","Sutterellaceae","Oscillospiraceae","Monoglobaceae"],[null,null,null,null,null,null,null,null],[null,null,null,null,null,null,null,null],[null,null,null,null,null,null,null,null],[0.00273684210526316,0.00273684210526316,0.00273684210526316,0.00273684210526316,0.00273684210526316,0.00273684210526316,0.00273684210526316,0.00273684210526316]],"container":"<table class=\"display\">\n  <thead>\n    <tr>\n      <th> <\/th>\n      <th>Axis.1<\/th>\n      <th>Axis.2<\/th>\n      <th>pval<\/th>\n      <th>r<\/th>\n      <th>Kingdom<\/th>\n      <th>Phylum<\/th>\n      <th>Class<\/th>\n      <th>Order<\/th>\n      <th>tax_rank_plot<\/th>\n      <th>Genus<\/th>\n      <th>Species<\/th>\n      <th>Strain<\/th>\n      <th>pval.adj<\/th>\n    <\/tr>\n  <\/thead>\n<\/table>","options":{"columnDefs":[{"className":"dt-right","targets":[1,2,3,4,13]},{"orderable":false,"targets":0}],"order":[],"autoWidth":false,"orderClasses":false}},"evals":[],"jsHooks":[]}</script>
```


```r
out$plot -> NMDS_tmp2_fam

NMDS_tmp2 + 
  theme(legend.position = "none") -> NMDS_tmp2_noleg

NMDS_tmp2_noleg
```

![](NRP72_draft_16S_Figs_beta_files/figure-html/unnamed-chunk-92-1.png)<!-- -->

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
<div id="htmlwidget-33198698d6599fd3fc97" style="width:100%;height:auto;" class="datatables html-widget"></div>
<script type="application/json" data-for="htmlwidget-33198698d6599fd3fc97">{"x":{"filter":"none","vertical":false,"data":[["1","2","3","4","5","6"],["Succinat_mM","Formiat_mM","Acetat_mM","Propionat_mM","Butyrat_mM","Total_SCFA_mM"],[0.163225559850326,0.44335408157543,-0.798708595679286,0.046761753008919,-0.750057796305815,-0.701255753014655],[0.510083629137183,0.0108062428054228,-0.466455045038887,-0.380510559986596,-0.524087796849334,-0.523321529938541],["Succinat_mM","Formiat_mM","Acetat_mM","Propionat_mM","Butyrat_mM","Total_SCFA_mM"],[0.001,0.005,0.001,0.026,0.001,0.001],[0.286827892102212,0.196679616533163,0.855515729854208,0.14697494780578,0.837254716605524,0.765625054833367]],"container":"<table class=\"display\">\n  <thead>\n    <tr>\n      <th> <\/th>\n      <th>label<\/th>\n      <th>Axis.1<\/th>\n      <th>Axis.2<\/th>\n      <th>id<\/th>\n      <th>pval<\/th>\n      <th>r<\/th>\n    <\/tr>\n  <\/thead>\n<\/table>","options":{"columnDefs":[{"className":"dt-right","targets":[2,3,5,6]},{"orderable":false,"targets":0}],"order":[],"autoWidth":false,"orderClasses":false}},"evals":[],"jsHooks":[]}</script>
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

![](NRP72_draft_16S_Figs_beta_files/figure-html/unnamed-chunk-96-1.png)<!-- -->


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
<div id="htmlwidget-38c86864160d6d855a80" style="width:100%;height:auto;" class="datatables html-widget"></div>
<script type="application/json" data-for="htmlwidget-38c86864160d6d855a80">{"x":{"filter":"none","vertical":false,"data":[["ASV0002","ASV0051","ASV0010","ASV0015","ASV0019","ASV0035","ASV0106","ASV0124"],[-0.596365968169365,0.65868508864788,-0.882275771735157,0.791983810973008,0.505766022416626,0.544200939838375,-0.742569357424228,-0.647095453887146],[-0.329816353820077,0.333036945606869,0.198219381930944,0.132828151930179,0.471256503357815,0.506200658362873,-0.0254828557713435,0.390896324977618],[0.001,0.001,0.001,0.001,0.001,0.001,0.001,0.001],[0.464431195237754,0.544779653146219,0.817701460763952,0.644881674788516,0.47788196138817,0.552393769447976,0.552058626523694,0.571532463322419],["Bacteria","Bacteria","Bacteria","Bacteria","Bacteria","Bacteria","Bacteria","Bacteria"],["Bacteroidota","Firmicutes","Firmicutes","Firmicutes","Proteobacteria","Proteobacteria","Firmicutes","Firmicutes"],["Bacteroidia","Bacilli","Clostridia","Negativicutes","Gammaproteobacteria","Gammaproteobacteria","Clostridia","Clostridia"],["Bacteroidales","Lactobacillales","Lachnospirales","Veillonellales-Selenomonadales","Pseudomonadales","Burkholderiales","Oscillospirales","Monoglobales"],["Bacteroidaceae","Lactobacillaceae","Lachnospiraceae","Veillonellaceae","Pseudomonadaceae","Sutterellaceae","Oscillospiraceae","Monoglobaceae"],[null,null,null,null,null,null,null,null],[null,null,null,null,null,null,null,null],[null,null,null,null,null,null,null,null],[0.0026,0.0026,0.0026,0.0026,0.0026,0.0026,0.0026,0.0026]],"container":"<table class=\"display\">\n  <thead>\n    <tr>\n      <th> <\/th>\n      <th>Axis.1<\/th>\n      <th>Axis.2<\/th>\n      <th>pval<\/th>\n      <th>r<\/th>\n      <th>Kingdom<\/th>\n      <th>Phylum<\/th>\n      <th>Class<\/th>\n      <th>Order<\/th>\n      <th>tax_rank_plot<\/th>\n      <th>Genus<\/th>\n      <th>Species<\/th>\n      <th>Strain<\/th>\n      <th>pval.adj<\/th>\n    <\/tr>\n  <\/thead>\n<\/table>","options":{"columnDefs":[{"className":"dt-right","targets":[1,2,3,4,13]},{"orderable":false,"targets":0}],"order":[],"autoWidth":false,"orderClasses":false}},"evals":[],"jsHooks":[]}</script>
```


```r
out$plot -> PCoA_tmp2_fam

PCoA_tmp2 + 
  theme(legend.position = "none") -> PCoA_tmp2_noleg

PCoA_tmp2_noleg
```

![](NRP72_draft_16S_Figs_beta_files/figure-html/unnamed-chunk-98-1.png)<!-- -->

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
<div id="htmlwidget-d96970b9501c191ad406" style="width:100%;height:auto;" class="datatables html-widget"></div>
<script type="application/json" data-for="htmlwidget-d96970b9501c191ad406">{"x":{"filter":"none","vertical":false,"data":[["1","2","3","4","5","6"],["Formiat_mM","Acetat_mM","Propionat_mM","Butyrat_mM","Valerat_mM","Total_SCFA_mM"],[0.27306840343581,-0.879625432657993,-0.0441586289042728,-0.795312448898601,-0.118568666569057,-0.780538747503882],[0.549238723654045,-0.238276371250953,-0.298723166856065,-0.45191431945377,-0.67858113286701,-0.431060339338],["Formiat_mM","Acetat_mM","Propionat_mM","Butyrat_mM","Valerat_mM","Total_SCFA_mM"],[0.001,0.001,0.01,0.001,0.001,0.001],[0.376229528516106,0.830516530875283,0.0911855149232218,0.836748443500454,0.474530882575039,0.795053752505121]],"container":"<table class=\"display\">\n  <thead>\n    <tr>\n      <th> <\/th>\n      <th>label<\/th>\n      <th>Axis.1<\/th>\n      <th>Axis.2<\/th>\n      <th>id<\/th>\n      <th>pval<\/th>\n      <th>r<\/th>\n    <\/tr>\n  <\/thead>\n<\/table>","options":{"columnDefs":[{"className":"dt-right","targets":[2,3,5,6]},{"orderable":false,"targets":0}],"order":[],"autoWidth":false,"orderClasses":false}},"evals":[],"jsHooks":[]}</script>
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

![](NRP72_draft_16S_Figs_beta_files/figure-html/unnamed-chunk-101-1.png)<!-- -->


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
<div id="htmlwidget-bbf49928d385f1dd9e5c" style="width:100%;height:auto;" class="datatables html-widget"></div>
<script type="application/json" data-for="htmlwidget-bbf49928d385f1dd9e5c">{"x":{"filter":"none","vertical":false,"data":[["ASV0002","ASV0012","ASV0051","ASV0010","ASV0015","ASV0035","ASV0106","ASV0124"],[-0.636725507108502,-0.432685232978908,0.698941189730734,-0.774433924337121,0.759890270630475,0.553845729547033,-0.696754170275816,-0.534362816472827],[-0.20387924927289,-0.505654315299078,0.267634778334487,0.409651810560177,-0.11892905216364,0.408241848362023,0.182516827457269,0.49997692210747],[0.001,0.001,0.001,0.001,0.001,0.001,0.001,0.001],[0.446986119686657,0.442902797418592,0.560147161276364,0.767562509059425,0.591577342847398,0.473406498891526,0.518778766101807,0.535520542268831],["Bacteria","Bacteria","Bacteria","Bacteria","Bacteria","Bacteria","Bacteria","Bacteria"],["Bacteroidota","Firmicutes","Firmicutes","Firmicutes","Firmicutes","Proteobacteria","Firmicutes","Firmicutes"],["Bacteroidia","Clostridia","Bacilli","Clostridia","Negativicutes","Gammaproteobacteria","Clostridia","Clostridia"],["Bacteroidales","Oscillospirales","Lactobacillales","Lachnospirales","Veillonellales-Selenomonadales","Burkholderiales","Oscillospirales","Monoglobales"],["Bacteroidaceae","Ruminococcaceae","Lactobacillaceae","Lachnospiraceae","Veillonellaceae","Sutterellaceae","Oscillospiraceae","Monoglobaceae"],[null,null,null,null,null,null,null,null],[null,null,null,null,null,null,null,null],[null,null,null,null,null,null,null,null],[0.00273684210526316,0.00273684210526316,0.00273684210526316,0.00273684210526316,0.00273684210526316,0.00273684210526316,0.00273684210526316,0.00273684210526316]],"container":"<table class=\"display\">\n  <thead>\n    <tr>\n      <th> <\/th>\n      <th>Axis.1<\/th>\n      <th>Axis.2<\/th>\n      <th>pval<\/th>\n      <th>r<\/th>\n      <th>Kingdom<\/th>\n      <th>Phylum<\/th>\n      <th>Class<\/th>\n      <th>Order<\/th>\n      <th>tax_rank_plot<\/th>\n      <th>Genus<\/th>\n      <th>Species<\/th>\n      <th>Strain<\/th>\n      <th>pval.adj<\/th>\n    <\/tr>\n  <\/thead>\n<\/table>","options":{"columnDefs":[{"className":"dt-right","targets":[1,2,3,4,13]},{"orderable":false,"targets":0}],"order":[],"autoWidth":false,"orderClasses":false}},"evals":[],"jsHooks":[]}</script>
```


```r
out$plot -> PCoA_tmp2_fam

PCoA_tmp2 + 
  theme(legend.position = "none") -> PCoA_tmp2_noleg

PCoA_tmp2_noleg
```

![](NRP72_draft_16S_Figs_beta_files/figure-html/unnamed-chunk-103-1.png)<!-- -->

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
<div id="htmlwidget-90bfa43893a2d8eda041" style="width:100%;height:auto;" class="datatables html-widget"></div>
<script type="application/json" data-for="htmlwidget-90bfa43893a2d8eda041">{"x":{"filter":"none","vertical":false,"data":[["1","2","3","4","5"],["Formiat_mM","Acetat_mM","Butyrat_mM","Valerat_mM","Total_SCFA_mM"],[0.283788233484607,-0.885113788447888,-0.824976931837785,-0.122665581390368,-0.814561476459512],[0.548707109918608,-0.0615798901903027,-0.271865478850179,-0.57907218311151,-0.244357186269817],["Formiat_mM","Acetat_mM","Butyrat_mM","Valerat_mM","Total_SCFA_mM"],[0.001,0.001,0.001,0.001,0.001],[0.381615253939545,0.787218501376423,0.754497776654922,0.350371438111367,0.723220833413601]],"container":"<table class=\"display\">\n  <thead>\n    <tr>\n      <th> <\/th>\n      <th>label<\/th>\n      <th>Axis.1<\/th>\n      <th>Axis.2<\/th>\n      <th>id<\/th>\n      <th>pval<\/th>\n      <th>r<\/th>\n    <\/tr>\n  <\/thead>\n<\/table>","options":{"columnDefs":[{"className":"dt-right","targets":[2,3,5,6]},{"orderable":false,"targets":0}],"order":[],"autoWidth":false,"orderClasses":false}},"evals":[],"jsHooks":[]}</script>
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
##  [1] vegan_2.6-2          lattice_0.20-45      permute_0.9-7       
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
##  [10] digest_0.6.29          foreach_1.5.2          htmltools_0.5.3       
##  [13] fansi_1.0.3            magrittr_2.0.3         googlesheets4_1.0.1   
##  [16] cluster_2.1.4          openxlsx_4.2.5         tzdb_0.3.0            
##  [19] Biostrings_2.64.1      extrafont_0.18         modelr_0.1.9          
##  [22] bayesm_3.1-4           matrixStats_0.62.0     officer_0.4.4         
##  [25] stabledist_0.7-1       extrafontdb_1.0        colorspace_2.0-3      
##  [28] rvest_1.0.3            ggrepel_0.9.1          textshaping_0.3.6     
##  [31] haven_2.5.1            xfun_0.33              crayon_1.5.2          
##  [34] RCurl_1.98-1.9         jsonlite_1.8.2         survival_3.4-0        
##  [37] iterators_1.0.14       glue_1.6.2             rvg_0.2.5             
##  [40] gtable_0.3.1           gargle_1.2.1           zlibbioc_1.42.0       
##  [43] XVector_0.36.0         car_3.1-0              Rttf2pt1_1.3.10       
##  [46] Rhdf5lib_1.18.2        BiocGenerics_0.42.0    DEoptimR_1.0-11       
##  [49] abind_1.4-5            DBI_1.1.3              rstatix_0.7.0         
##  [52] Rcpp_1.0.9             xtable_1.8-4           clue_0.3-61           
##  [55] DT_0.25                stats4_4.2.1           htmlwidgets_1.5.4     
##  [58] timeSeries_4021.104    httr_1.4.4             ellipsis_0.3.2        
##  [61] spatial_7.3-15         pkgconfig_2.0.3        farver_2.1.1          
##  [64] sass_0.4.2             dbplyr_2.2.1           utf8_1.2.2            
##  [67] here_1.0.1             tidyselect_1.2.0       labeling_0.4.2        
##  [70] rlang_1.0.6            munsell_0.5.0          cellranger_1.1.0      
##  [73] tools_4.2.1            cachem_1.0.6           cli_3.4.1             
##  [76] generics_0.1.3         devEMF_4.1             ade4_1.7-19           
##  [79] export_0.3.0           broom_1.0.1            evaluate_0.17         
##  [82] biomformat_1.24.0      fastmap_1.1.0          yaml_2.3.5            
##  [85] ragg_1.2.3             knitr_1.40             fs_1.5.2              
##  [88] zip_2.2.1              robustbase_0.95-0      rgl_0.110.2           
##  [91] xml2_1.3.3             compiler_4.2.1         rstudioapi_0.14       
##  [94] ggsignif_0.6.3         reprex_2.0.2           statmod_1.4.37        
##  [97] bslib_0.4.0            stringi_1.7.8          statip_0.2.3          
## [100] highr_0.9              stargazer_5.2.3        modeest_2.4.0         
## [103] fBasics_4021.92        Matrix_1.5-1           microbiome_1.19.1     
## [106] tensorA_0.36.2         multtest_2.52.0        vctrs_0.4.2           
## [109] pillar_1.8.1           lifecycle_1.0.3        rhdf5filters_1.8.0    
## [112] microViz_0.9.2         jquerylib_0.1.4        flextable_0.8.2       
## [115] cowplot_1.1.1          data.table_1.14.2      bitops_1.0-7          
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
