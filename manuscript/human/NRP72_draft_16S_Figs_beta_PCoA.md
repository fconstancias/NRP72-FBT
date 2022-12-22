---
title: "NTP72 - 16S - beta - human draft "
author: "Florentin Constancias"
date: "December 21, 2022"
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
  subset_samples(Model == "Human")  -> ps_filtered
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
<div class="datatables html-widget html-fill-item-overflow-hidden html-fill-item" id="htmlwidget-bb3486349adfa5d821dc" style="width:100%;height:auto;"></div>
<script type="application/json" data-for="htmlwidget-bb3486349adfa5d821dc">{"x":{"filter":"none","vertical":false,"data":[["1","2","3","4","5","6"],["bjaccard","bjaccard","bjaccard","aichinson","aichinson","aichinson"],["Treatment","Residual","Total","Treatment","Residual","Total"],[3,14,17,3,14,17],[0.242,0.557,0.799,0.339,0.803,1.142],[0.303,0.697,1,0.297,0.703,1],[2.03,null,null,1.972,null,null],[0.024,null,null,0.014,null,null]],"container":"<table class=\"display\">\n  <thead>\n    <tr>\n      <th> <\/th>\n      <th>Distance<\/th>\n      <th>terms<\/th>\n      <th>Df<\/th>\n      <th>SumOfSqs<\/th>\n      <th>R2<\/th>\n      <th>F<\/th>\n      <th>Pr(&gt;F)<\/th>\n    <\/tr>\n  <\/thead>\n<\/table>","options":{"columnDefs":[{"className":"dt-right","targets":[3,4,5,6,7]},{"orderable":false,"targets":0}],"order":[],"autoWidth":false,"orderClasses":false}},"evals":[],"jsHooks":[]}</script>
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
<div class="datatables html-widget html-fill-item-overflow-hidden html-fill-item" id="htmlwidget-44474b3e8cd160f85a98" style="width:100%;height:auto;"></div>
<script type="application/json" data-for="htmlwidget-44474b3e8cd160f85a98">{"x":{"filter":"none","vertical":false,"data":[["1","2","3","4","5","6","7","8"],["bjaccard","bjaccard","bjaccard","bjaccard","aichinson","aichinson","aichinson","aichinson"],["Treatment","Day_of_Treatment","Residual","Total","Treatment","Day_of_Treatment","Residual","Total"],[3,1,13,17,3,1,13,17],[0.238,0.036,0.521,0.799,0.336,0.045,0.757,1.142],[0.298,0.044,0.652,1,0.295,0.04,0.663,1],[1.978,0.885,null,null,1.924,0.779,null,null],[0.023,0.466,null,null,0.016,0.633,null,null]],"container":"<table class=\"display\">\n  <thead>\n    <tr>\n      <th> <\/th>\n      <th>Distance<\/th>\n      <th>terms<\/th>\n      <th>Df<\/th>\n      <th>SumOfSqs<\/th>\n      <th>R2<\/th>\n      <th>F<\/th>\n      <th>Pr(&gt;F)<\/th>\n    <\/tr>\n  <\/thead>\n<\/table>","options":{"columnDefs":[{"className":"dt-right","targets":[3,4,5,6,7]},{"orderable":false,"targets":0}],"order":[],"autoWidth":false,"orderClasses":false}},"evals":[],"jsHooks":[]}</script>
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

```r
tmp_df %>% 
  DT::datatable()
```

```{=html}
<div class="datatables html-widget html-fill-item-overflow-hidden html-fill-item" id="htmlwidget-2514ce28ae593bde0f73" style="width:100%;height:auto;"></div>
<script type="application/json" data-for="htmlwidget-2514ce28ae593bde0f73">{"x":{"filter":"none","vertical":false,"data":[["1","2","3","4","5","6","7","8","9","10","11","12"],["bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson"],["UNTREATED","UNTREATED","UNTREATED","VAN","VAN","VAN+CCUG59168","UNTREATED","UNTREATED","UNTREATED","VAN","VAN","VAN+CCUG59168"],["VAN","VAN+CCUG59168","CCUG59168","VAN+CCUG59168","CCUG59168","CCUG59168","VAN","VAN+CCUG59168","CCUG59168","VAN+CCUG59168","CCUG59168","CCUG59168"],[0.287,0.276,0.483,0.13,0.142,0.153,0.252,0.261,0.44,0.14,0.165,0.156],[0.05,0.031,0.1,0.197,0.718,0.148,0.166,0.031,0.1,0.174,0.506,0.143],[0.15,0.186,0.2,0.236,0.718,0.222,0.249,0.186,0.3,0.209,0.506,0.286]],"container":"<table class=\"display\">\n  <thead>\n    <tr>\n      <th> <\/th>\n      <th>Distance<\/th>\n      <th>X1<\/th>\n      <th>X2<\/th>\n      <th>R2<\/th>\n      <th>pval<\/th>\n      <th>pvalFDR<\/th>\n    <\/tr>\n  <\/thead>\n<\/table>","options":{"columnDefs":[{"className":"dt-right","targets":[4,5,6]},{"orderable":false,"targets":0}],"order":[],"autoWidth":false,"orderClasses":false}},"evals":[],"jsHooks":[]}</script>
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
<div class="datatables html-widget html-fill-item-overflow-hidden html-fill-item" id="htmlwidget-5d05eac75555269eefd0" style="width:100%;height:auto;"></div>
<script type="application/json" data-for="htmlwidget-5d05eac75555269eefd0">{"x":{"filter":"none","vertical":false,"data":[["1","2","3","4","5","6","7","8","9","10","11","12"],["bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson"],[["UNTREATED"],["UNTREATED"],["UNTREATED"],["VAN"],["VAN"],["VAN+CCUG59168"],["UNTREATED"],["UNTREATED"],["UNTREATED"],["VAN"],["VAN"],["VAN+CCUG59168"]],[["VAN"],["VAN+CCUG59168"],["CCUG59168"],["VAN+CCUG59168"],["CCUG59168"],["CCUG59168"],["VAN"],["VAN+CCUG59168"],["CCUG59168"],["VAN+CCUG59168"],["CCUG59168"],["CCUG59168"]],[[3],[3],[3],[4],[4],[8],[3],[3],[3],[4],[4],[8]],[[4],[8],[3],[8],[3],[3],[4],[8],[3],[8],[3],[3]],[[0.022],[0.027],[0.125],[0.208],[0.584],[0.148],[0.087],[0.03],[0.095],[0.165],[0.401],[0.145]],[[2.30752531223853],[5.22648552120524],[3.73344949882427],[1.52253123836056],[0.897603749688104],[2.12541595550425],[1.86986749609433],[4.08899608446121],[3.13989123864949],[1.56099464336391],[1.07848120760732],[2.05363506633381]],[[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999]],[0.081,0.081,0.222,0.25,0.584,0.222,0.19,0.18,0.19,0.198,0.401,0.198]],"container":"<table class=\"display\">\n  <thead>\n    <tr>\n      <th> <\/th>\n      <th>Distance<\/th>\n      <th>Level1<\/th>\n      <th>Level2<\/th>\n      <th>N1<\/th>\n      <th>N2<\/th>\n      <th>p.value<\/th>\n      <th>tw2.stat<\/th>\n      <th>nrep<\/th>\n      <th>pvalFDR<\/th>\n    <\/tr>\n  <\/thead>\n<\/table>","options":{"columnDefs":[{"className":"dt-right","targets":9},{"orderable":false,"targets":0}],"order":[],"autoWidth":false,"orderClasses":false}},"evals":[],"jsHooks":[]}</script>
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
<div class="datatables html-widget html-fill-item-overflow-hidden html-fill-item" id="htmlwidget-7531f021030c44f57ef2" style="width:100%;height:auto;"></div>
<script type="application/json" data-for="htmlwidget-7531f021030c44f57ef2">{"x":{"filter":"none","vertical":false,"data":[["1","2","3","4","5","6"],["bjaccard","bjaccard","bjaccard","aichinson","aichinson","aichinson"],["Treatment","Residual","Total","Treatment","Residual","Total"],[3,25,28,3,25,28],[3.297,3.24,6.537,2.442,1.513,3.955],[0.504,0.496,1,0.617,0.383,1],[8.48,null,null,13.451,null,null],[0.001,null,null,0.001,null,null]],"container":"<table class=\"display\">\n  <thead>\n    <tr>\n      <th> <\/th>\n      <th>Distance<\/th>\n      <th>terms<\/th>\n      <th>Df<\/th>\n      <th>SumOfSqs<\/th>\n      <th>R2<\/th>\n      <th>F<\/th>\n      <th>Pr(&gt;F)<\/th>\n    <\/tr>\n  <\/thead>\n<\/table>","options":{"columnDefs":[{"className":"dt-right","targets":[3,4,5,6,7]},{"orderable":false,"targets":0}],"order":[],"autoWidth":false,"orderClasses":false}},"evals":[],"jsHooks":[]}</script>
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
<div class="datatables html-widget html-fill-item-overflow-hidden html-fill-item" id="htmlwidget-983315589fed18d622a7" style="width:100%;height:auto;"></div>
<script type="application/json" data-for="htmlwidget-983315589fed18d622a7">{"x":{"filter":"none","vertical":false,"data":[["1","2","3","4","5","6","7","8"],["bjaccard","bjaccard","bjaccard","bjaccard","aichinson","aichinson","aichinson","aichinson"],["Treatment","Day_of_Treatment","Residual","Total","Treatment","Day_of_Treatment","Residual","Total"],[3,1,24,28,3,1,24,28],[3.172,0.357,2.883,6.537,2.38,0.154,1.359,3.955],[0.485,0.055,0.441,1,0.602,0.039,0.344,1],[8.801,2.968,null,null,14.009,2.713,null,null],[0.001,0.029,null,null,0.001,0.047,null,null]],"container":"<table class=\"display\">\n  <thead>\n    <tr>\n      <th> <\/th>\n      <th>Distance<\/th>\n      <th>terms<\/th>\n      <th>Df<\/th>\n      <th>SumOfSqs<\/th>\n      <th>R2<\/th>\n      <th>F<\/th>\n      <th>Pr(&gt;F)<\/th>\n    <\/tr>\n  <\/thead>\n<\/table>","options":{"columnDefs":[{"className":"dt-right","targets":[3,4,5,6,7]},{"orderable":false,"targets":0}],"order":[],"autoWidth":false,"orderClasses":false}},"evals":[],"jsHooks":[]}</script>
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
<div class="datatables html-widget html-fill-item-overflow-hidden html-fill-item" id="htmlwidget-6f1d8b9e786d9d3fe28d" style="width:100%;height:auto;"></div>
<script type="application/json" data-for="htmlwidget-6f1d8b9e786d9d3fe28d">{"x":{"filter":"none","vertical":false,"data":[["1","2","3","4","5","6","7","8","9","10","11","12"],["bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson"],["UNTREATED","UNTREATED","UNTREATED","VAN","VAN","VAN+CCUG59168","UNTREATED","UNTREATED","UNTREATED","VAN","VAN","VAN+CCUG59168"],["VAN","VAN+CCUG59168","CCUG59168","VAN+CCUG59168","CCUG59168","CCUG59168","VAN","VAN+CCUG59168","CCUG59168","VAN+CCUG59168","CCUG59168","CCUG59168"],[0.386,0.6,0.36,0.133,0.384,0.593,0.533,0.691,0.344,0.129,0.526,0.68],[0.001,0.001,0.016,0.01,0.001,0.003,0.001,0.002,0.015,0.007,0.002,0.001],[0.003,0.003,0.016,0.012,0.003,0.005,0.004,0.003,0.015,0.008,0.003,0.004]],"container":"<table class=\"display\">\n  <thead>\n    <tr>\n      <th> <\/th>\n      <th>Distance<\/th>\n      <th>X1<\/th>\n      <th>X2<\/th>\n      <th>R2<\/th>\n      <th>pval<\/th>\n      <th>pvalFDR<\/th>\n    <\/tr>\n  <\/thead>\n<\/table>","options":{"columnDefs":[{"className":"dt-right","targets":[4,5,6]},{"orderable":false,"targets":0}],"order":[],"autoWidth":false,"orderClasses":false}},"evals":[],"jsHooks":[]}</script>
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
<div class="datatables html-widget html-fill-item-overflow-hidden html-fill-item" id="htmlwidget-dcdfc64296f5794e54a4" style="width:100%;height:auto;"></div>
<script type="application/json" data-for="htmlwidget-dcdfc64296f5794e54a4">{"x":{"filter":"none","vertical":false,"data":[["1","2","3","4","5","6","7","8","9","10","11","12"],["bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson"],[["UNTREATED"],["UNTREATED"],["UNTREATED"],["VAN"],["VAN"],["VAN+CCUG59168"],["UNTREATED"],["UNTREATED"],["UNTREATED"],["VAN"],["VAN"],["VAN+CCUG59168"]],[["VAN"],["VAN+CCUG59168"],["CCUG59168"],["VAN+CCUG59168"],["CCUG59168"],["CCUG59168"],["VAN"],["VAN+CCUG59168"],["CCUG59168"],["VAN+CCUG59168"],["CCUG59168"],["CCUG59168"]],[[5],[5],[5],[11],[11],[8],[5],[5],[5],[11],[11],[8]],[[11],[8],[5],[8],[5],[5],[11],[8],[5],[8],[5],[5]],[[0.002],[0.002],[0.008],[0.01],[0.001],[0.002],[0.002],[0.001],[0.009],[0.008],[0.001],[0.001]],[[17.8005046571075],[24.8335057766479],[4.50376066993907],[2.74386735082863],[17.0489797847579],[23.5115717122834],[23.6544578261575],[29.4581279877637],[4.18802840816871],[2.65533238166408],[21.6350858929726],[26.8527111401793]],[[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999]],[0.003,0.003,0.01,0.01,0.003,0.003,0.003,0.002,0.009,0.009,0.002,0.002]],"container":"<table class=\"display\">\n  <thead>\n    <tr>\n      <th> <\/th>\n      <th>Distance<\/th>\n      <th>Level1<\/th>\n      <th>Level2<\/th>\n      <th>N1<\/th>\n      <th>N2<\/th>\n      <th>p.value<\/th>\n      <th>tw2.stat<\/th>\n      <th>nrep<\/th>\n      <th>pvalFDR<\/th>\n    <\/tr>\n  <\/thead>\n<\/table>","options":{"columnDefs":[{"className":"dt-right","targets":9},{"orderable":false,"targets":0}],"order":[],"autoWidth":false,"orderClasses":false}},"evals":[],"jsHooks":[]}</script>
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
<div class="datatables html-widget html-fill-item-overflow-hidden html-fill-item" id="htmlwidget-873249b169eb6149d57e" style="width:100%;height:auto;"></div>
<script type="application/json" data-for="htmlwidget-873249b169eb6149d57e">{"x":{"filter":"none","vertical":false,"data":[["1","2","3","4","5","6"],["bjaccard","bjaccard","bjaccard","aichinson","aichinson","aichinson"],["Treatment","Residual","Total","Treatment","Residual","Total"],[3,15,18,3,15,18],[2.331,1.971,4.302,1.712,0.917,2.629],[0.542,0.458,1,0.651,0.349,1],[5.911,null,null,9.338,null,null],[0.001,null,null,0.001,null,null]],"container":"<table class=\"display\">\n  <thead>\n    <tr>\n      <th> <\/th>\n      <th>Distance<\/th>\n      <th>terms<\/th>\n      <th>Df<\/th>\n      <th>SumOfSqs<\/th>\n      <th>R2<\/th>\n      <th>F<\/th>\n      <th>Pr(&gt;F)<\/th>\n    <\/tr>\n  <\/thead>\n<\/table>","options":{"columnDefs":[{"className":"dt-right","targets":[3,4,5,6,7]},{"orderable":false,"targets":0}],"order":[],"autoWidth":false,"orderClasses":false}},"evals":[],"jsHooks":[]}</script>
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
<div class="datatables html-widget html-fill-item-overflow-hidden html-fill-item" id="htmlwidget-cbd3ad10bd2d35a2e76d" style="width:100%;height:auto;"></div>
<script type="application/json" data-for="htmlwidget-cbd3ad10bd2d35a2e76d">{"x":{"filter":"none","vertical":false,"data":[["1","2","3","4","5","6","7","8"],["bjaccard","bjaccard","bjaccard","bjaccard","aichinson","aichinson","aichinson","aichinson"],["Treatment","Day_of_Treatment","Residual","Total","Treatment","Day_of_Treatment","Residual","Total"],[3,1,14,18,3,1,14,18],[2.33,0.187,1.785,4.302,1.71,0.09,0.827,2.629],[0.542,0.043,0.415,1,0.65,0.034,0.315,1],[6.092,1.464,null,null,9.648,1.518,null,null],[0.001,0.184,null,null,0.001,0.197,null,null]],"container":"<table class=\"display\">\n  <thead>\n    <tr>\n      <th> <\/th>\n      <th>Distance<\/th>\n      <th>terms<\/th>\n      <th>Df<\/th>\n      <th>SumOfSqs<\/th>\n      <th>R2<\/th>\n      <th>F<\/th>\n      <th>Pr(&gt;F)<\/th>\n    <\/tr>\n  <\/thead>\n<\/table>","options":{"columnDefs":[{"className":"dt-right","targets":[3,4,5,6,7]},{"orderable":false,"targets":0}],"order":[],"autoWidth":false,"orderClasses":false}},"evals":[],"jsHooks":[]}</script>
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
```

```r
tmp_df %>% 
  DT::datatable()
```

```{=html}
<div class="datatables html-widget html-fill-item-overflow-hidden html-fill-item" id="htmlwidget-f25d4e1770d66a86dcce" style="width:100%;height:auto;"></div>
<script type="application/json" data-for="htmlwidget-f25d4e1770d66a86dcce">{"x":{"filter":"none","vertical":false,"data":[["1","2","3","4","5","6","7","8","9","10","11","12"],["bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson"],["UNTREATED","UNTREATED","UNTREATED","VAN","VAN","VAN+CCUG59168","UNTREATED","UNTREATED","UNTREATED","VAN","VAN","VAN+CCUG59168"],["VAN","VAN+CCUG59168","CCUG59168","VAN+CCUG59168","CCUG59168","CCUG59168","VAN","VAN+CCUG59168","CCUG59168","VAN+CCUG59168","CCUG59168","CCUG59168"],[0.44,0.656,0.507,0.148,0.397,0.622,0.573,0.731,0.453,0.14,0.541,0.718],[0.005,0.008,0.032,0.043,0.007,0.018,0.004,0.008,0.037,0.07,0.011,0.026],[0.03,0.016,0.038,0.043,0.021,0.027,0.024,0.024,0.044,0.07,0.022,0.039]],"container":"<table class=\"display\">\n  <thead>\n    <tr>\n      <th> <\/th>\n      <th>Distance<\/th>\n      <th>X1<\/th>\n      <th>X2<\/th>\n      <th>R2<\/th>\n      <th>pval<\/th>\n      <th>pvalFDR<\/th>\n    <\/tr>\n  <\/thead>\n<\/table>","options":{"columnDefs":[{"className":"dt-right","targets":[4,5,6]},{"orderable":false,"targets":0}],"order":[],"autoWidth":false,"orderClasses":false}},"evals":[],"jsHooks":[]}</script>
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
<div class="datatables html-widget html-fill-item-overflow-hidden html-fill-item" id="htmlwidget-969965edfe9f869aa904" style="width:100%;height:auto;"></div>
<script type="application/json" data-for="htmlwidget-969965edfe9f869aa904">{"x":{"filter":"none","vertical":false,"data":[["1","2","3","4","5","6","7","8","9","10","11","12"],["bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson"],[["UNTREATED"],["UNTREATED"],["UNTREATED"],["VAN"],["VAN"],["VAN+CCUG59168"],["UNTREATED"],["UNTREATED"],["UNTREATED"],["VAN"],["VAN"],["VAN+CCUG59168"]],[["VAN"],["VAN+CCUG59168"],["CCUG59168"],["VAN+CCUG59168"],["CCUG59168"],["CCUG59168"],["VAN"],["VAN+CCUG59168"],["CCUG59168"],["VAN+CCUG59168"],["CCUG59168"],["CCUG59168"]],[[4],[4],[4],[7],[7],[5],[4],[4],[4],[7],[7],[5]],[[7],[5],[3],[5],[3],[3],[7],[5],[3],[5],[3],[3]],[[0.005],[0.013],[0.041],[0.041],[0.008],[0.015],[0.003],[0.008],[0.035],[0.047],[0.014],[0.018]],[[12.05918175925],[16.477721867272],[5.18037871223416],[1.84581617431872],[11.7983849404617],[16.0582527959237],[16.4475975972854],[20.9142536645326],[4.12294342827623],[1.74938435550998],[14.8802055908688],[18.8051845046142]],[[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999]],[0.022,0.022,0.041,0.041,0.022,0.022,0.018,0.024,0.042,0.047,0.027,0.027]],"container":"<table class=\"display\">\n  <thead>\n    <tr>\n      <th> <\/th>\n      <th>Distance<\/th>\n      <th>Level1<\/th>\n      <th>Level2<\/th>\n      <th>N1<\/th>\n      <th>N2<\/th>\n      <th>p.value<\/th>\n      <th>tw2.stat<\/th>\n      <th>nrep<\/th>\n      <th>pvalFDR<\/th>\n    <\/tr>\n  <\/thead>\n<\/table>","options":{"columnDefs":[{"className":"dt-right","targets":9},{"orderable":false,"targets":0}],"order":[],"autoWidth":false,"orderClasses":false}},"evals":[],"jsHooks":[]}</script>
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

![](NRP72_draft_16S_Figs_beta_PCoA_files/figure-html/unnamed-chunk-43-1.png)<!-- -->


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
<div class="datatables html-widget html-fill-item-overflow-hidden html-fill-item" id="htmlwidget-450f171f14364aa3ebed" style="width:100%;height:auto;"></div>
<script type="application/json" data-for="htmlwidget-450f171f14364aa3ebed">{"x":{"filter":"none","vertical":false,"data":[["ASV0012","ASV0004","ASV0010","ASV0052","ASV0035","ASV0059","ASV0106","ASV0275"],[0.445553060796587,-0.527431856534814,-0.682394062335409,0.675926734603464,-0.442637448488747,-0.19625688252075,-0.263526242835191,0.71119451169905],[0.343116359947506,0.288169711770378,-0.248021115693937,0.0173998730094112,0.334570499832224,0.510439207106555,-0.620428762826066,0.026637446057722],[0.001,0.002,0.001,0.001,0.002,0.002,0.001,0.001],[0.316246366448833,0.361226146069584,0.527176130140688,0.457179706132445,0.307865330162612,0.299064948088332,0.454377930404715,0.506507187003329],["Bacteria","Bacteria","Bacteria","Bacteria","Bacteria","Bacteria","Bacteria","Bacteria"],["Firmicutes","Proteobacteria","Firmicutes","Firmicutes","Proteobacteria","Bacteroidota","Firmicutes","Firmicutes"],["Clostridia","Gammaproteobacteria","Clostridia","Negativicutes","Gammaproteobacteria","Bacteroidia","Clostridia","Bacilli"],["Oscillospirales","Enterobacterales","Lachnospirales","Veillonellales-Selenomonadales","Burkholderiales","Bacteroidales","Oscillospirales","Erysipelotrichales"],["Ruminococcaceae","Enterobacteriaceae","Lachnospiraceae","Veillonellaceae","Sutterellaceae","Tannerellaceae","Oscillospiraceae","Erysipelotrichaceae"],[null,null,null,null,null,null,null,null],[null,null,null,null,null,null,null,null],[null,null,null,null,null,null,null,null],[0.00742857142857143,0.00945454545454546,0.00742857142857143,0.00742857142857143,0.00945454545454546,0.00945454545454546,0.00742857142857143,0.00742857142857143]],"container":"<table class=\"display\">\n  <thead>\n    <tr>\n      <th> <\/th>\n      <th>Axis.1<\/th>\n      <th>Axis.2<\/th>\n      <th>pval<\/th>\n      <th>r<\/th>\n      <th>Kingdom<\/th>\n      <th>Phylum<\/th>\n      <th>Class<\/th>\n      <th>Order<\/th>\n      <th>tax_rank_plot<\/th>\n      <th>Genus<\/th>\n      <th>Species<\/th>\n      <th>Strain<\/th>\n      <th>pval.adj<\/th>\n    <\/tr>\n  <\/thead>\n<\/table>","options":{"columnDefs":[{"className":"dt-right","targets":[1,2,3,4,13]},{"orderable":false,"targets":0}],"order":[],"autoWidth":false,"orderClasses":false}},"evals":[],"jsHooks":[]}</script>
```


```r
out$plot -> NMDS_tmp2_fam

NMDS_tmp2 + 
  theme(legend.position = "none") -> NMDS_tmp2_noleg

NMDS_tmp2_noleg
```

![](NRP72_draft_16S_Figs_beta_PCoA_files/figure-html/unnamed-chunk-45-1.png)<!-- -->

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
<div class="datatables html-widget html-fill-item-overflow-hidden html-fill-item" id="htmlwidget-fbfc52d82635d988c6e1" style="width:100%;height:auto;"></div>
<script type="application/json" data-for="htmlwidget-fbfc52d82635d988c6e1">{"x":{"filter":"none","vertical":false,"data":[["1","2","3","4"],["Formiat_mM","Acetat_mM","Valerat_mM","Total_SCFA_mM"],[-0.613360599591134,-0.127975964263856,0.178226366838872,-0.199371069496966],[0.042165219774713,-0.372980394645997,-0.677351910348987,-0.424683773969117],["Formiat_mM","Acetat_mM","Valerat_mM","Total_SCFA_mM"],[0.001,0.04,0.001,0.007],[0.377989130889445,0.155492222219548,0.490570248290006,0.220105131225016]],"container":"<table class=\"display\">\n  <thead>\n    <tr>\n      <th> <\/th>\n      <th>label<\/th>\n      <th>Axis.1<\/th>\n      <th>Axis.2<\/th>\n      <th>id<\/th>\n      <th>pval<\/th>\n      <th>r<\/th>\n    <\/tr>\n  <\/thead>\n<\/table>","options":{"columnDefs":[{"className":"dt-right","targets":[2,3,5,6]},{"orderable":false,"targets":0}],"order":[],"autoWidth":false,"orderClasses":false}},"evals":[],"jsHooks":[]}</script>
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

![](NRP72_draft_16S_Figs_beta_PCoA_files/figure-html/unnamed-chunk-48-1.png)<!-- -->


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
<div class="datatables html-widget html-fill-item-overflow-hidden html-fill-item" id="htmlwidget-2f4c52b7f36e26c300a7" style="width:100%;height:auto;"></div>
<script type="application/json" data-for="htmlwidget-2f4c52b7f36e26c300a7">{"x":{"filter":"none","vertical":false,"data":[["ASV0012","ASV0004","ASV0010","ASV0052","ASV0030","ASV0035","ASV0106","ASV0275"],[0.419469264716501,-0.542610944831473,-0.668057142724215,0.634491810605059,0.0923438525333804,-0.475786979717028,-0.235980253754372,0.716607484689898],[0.41313443709802,0.20743545326037,-0.213247478016456,-0.157224913125192,0.549578877374205,0.217200882303641,-0.59337501428833,-0.0384851666690667],[0.001,0.003,0.001,0.001,0.011,0.003,0.001,0.001],[0.3466345271581,0.337456104720238,0.491774832825221,0.427299531032111,0.310564329556598,0.273549473341732,0.407780587743654,0.515007395167129],["Bacteria","Bacteria","Bacteria","Bacteria","Bacteria","Bacteria","Bacteria","Bacteria"],["Firmicutes","Proteobacteria","Firmicutes","Firmicutes","Bacteroidota","Proteobacteria","Firmicutes","Firmicutes"],["Clostridia","Gammaproteobacteria","Clostridia","Negativicutes","Bacteroidia","Gammaproteobacteria","Clostridia","Bacilli"],["Oscillospirales","Enterobacterales","Lachnospirales","Veillonellales-Selenomonadales","Bacteroidales","Burkholderiales","Oscillospirales","Erysipelotrichales"],["Ruminococcaceae","Enterobacteriaceae","Lachnospiraceae","Veillonellaceae","Prevotellaceae","Sutterellaceae","Oscillospiraceae","Erysipelotrichaceae"],[null,null,null,null,null,null,null,null],[null,null,null,null,null,null,null,null],[null,null,null,null,null,null,null,null],[0.0104,0.0222857142857143,0.0104,0.0104,0.0476666666666667,0.0222857142857143,0.0104,0.0104]],"container":"<table class=\"display\">\n  <thead>\n    <tr>\n      <th> <\/th>\n      <th>Axis.1<\/th>\n      <th>Axis.2<\/th>\n      <th>pval<\/th>\n      <th>r<\/th>\n      <th>Kingdom<\/th>\n      <th>Phylum<\/th>\n      <th>Class<\/th>\n      <th>Order<\/th>\n      <th>tax_rank_plot<\/th>\n      <th>Genus<\/th>\n      <th>Species<\/th>\n      <th>Strain<\/th>\n      <th>pval.adj<\/th>\n    <\/tr>\n  <\/thead>\n<\/table>","options":{"columnDefs":[{"className":"dt-right","targets":[1,2,3,4,13]},{"orderable":false,"targets":0}],"order":[],"autoWidth":false,"orderClasses":false}},"evals":[],"jsHooks":[]}</script>
```


```r
out$plot -> NMDS_tmp2_fam

NMDS_tmp2 + 
  theme(legend.position = "none") -> NMDS_tmp2_noleg

NMDS_tmp2_noleg
```

![](NRP72_draft_16S_Figs_beta_PCoA_files/figure-html/unnamed-chunk-50-1.png)<!-- -->

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
<div class="datatables html-widget html-fill-item-overflow-hidden html-fill-item" id="htmlwidget-1d5df3d80f5735caaf13" style="width:100%;height:auto;"></div>
<script type="application/json" data-for="htmlwidget-1d5df3d80f5735caaf13">{"x":{"filter":"none","vertical":false,"data":[["1","2","3","4","5"],["Succinat_mM","Formiat_mM","Acetat_mM","Butyrat_mM","Total_SCFA_mM"],[0.164311218332464,0.441439686730363,-0.797826548818299,-0.751687244981628,-0.711900279912203],[0.506523475467576,-0.0242931772264333,-0.468604143606486,-0.51913836034085,-0.498311404715984],["Succinat_mM","Formiat_mM","Acetat_mM","Butyrat_mM","Total_SCFA_mM"],[0.001,0.004,0.001,0.001,0.001],[0.283564207669651,0.195459155480356,0.856117045404486,0.834538351445456,0.75511626460909]],"container":"<table class=\"display\">\n  <thead>\n    <tr>\n      <th> <\/th>\n      <th>label<\/th>\n      <th>Axis.1<\/th>\n      <th>Axis.2<\/th>\n      <th>id<\/th>\n      <th>pval<\/th>\n      <th>r<\/th>\n    <\/tr>\n  <\/thead>\n<\/table>","options":{"columnDefs":[{"className":"dt-right","targets":[2,3,5,6]},{"orderable":false,"targets":0}],"order":[],"autoWidth":false,"orderClasses":false}},"evals":[],"jsHooks":[]}</script>
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

![](NRP72_draft_16S_Figs_beta_PCoA_files/figure-html/unnamed-chunk-54-1.png)<!-- -->


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
<div class="datatables html-widget html-fill-item-overflow-hidden html-fill-item" id="htmlwidget-67c25ca3336b15b0293a" style="width:100%;height:auto;"></div>
<script type="application/json" data-for="htmlwidget-67c25ca3336b15b0293a">{"x":{"filter":"none","vertical":false,"data":[["ASV0012","ASV0004","ASV0010","ASV0052","ASV0035","ASV0059","ASV0106","ASV0275"],[0.445553060796587,-0.527431856534814,-0.682394062335409,0.675926734603464,-0.442637448488747,-0.19625688252075,-0.263526242835191,0.71119451169905],[0.343116359947506,0.288169711770378,-0.248021115693937,0.0173998730094112,0.334570499832224,0.510439207106555,-0.620428762826066,0.026637446057722],[0.001,0.002,0.001,0.001,0.002,0.002,0.001,0.001],[0.316246366448833,0.361226146069584,0.527176130140688,0.457179706132445,0.307865330162612,0.299064948088332,0.454377930404715,0.506507187003329],["Bacteria","Bacteria","Bacteria","Bacteria","Bacteria","Bacteria","Bacteria","Bacteria"],["Firmicutes","Proteobacteria","Firmicutes","Firmicutes","Proteobacteria","Bacteroidota","Firmicutes","Firmicutes"],["Clostridia","Gammaproteobacteria","Clostridia","Negativicutes","Gammaproteobacteria","Bacteroidia","Clostridia","Bacilli"],["Oscillospirales","Enterobacterales","Lachnospirales","Veillonellales-Selenomonadales","Burkholderiales","Bacteroidales","Oscillospirales","Erysipelotrichales"],["Ruminococcaceae","Enterobacteriaceae","Lachnospiraceae","Veillonellaceae","Sutterellaceae","Tannerellaceae","Oscillospiraceae","Erysipelotrichaceae"],[null,null,null,null,null,null,null,null],[null,null,null,null,null,null,null,null],[null,null,null,null,null,null,null,null],[0.00742857142857143,0.00945454545454546,0.00742857142857143,0.00742857142857143,0.00945454545454546,0.00945454545454546,0.00742857142857143,0.00742857142857143]],"container":"<table class=\"display\">\n  <thead>\n    <tr>\n      <th> <\/th>\n      <th>Axis.1<\/th>\n      <th>Axis.2<\/th>\n      <th>pval<\/th>\n      <th>r<\/th>\n      <th>Kingdom<\/th>\n      <th>Phylum<\/th>\n      <th>Class<\/th>\n      <th>Order<\/th>\n      <th>tax_rank_plot<\/th>\n      <th>Genus<\/th>\n      <th>Species<\/th>\n      <th>Strain<\/th>\n      <th>pval.adj<\/th>\n    <\/tr>\n  <\/thead>\n<\/table>","options":{"columnDefs":[{"className":"dt-right","targets":[1,2,3,4,13]},{"orderable":false,"targets":0}],"order":[],"autoWidth":false,"orderClasses":false}},"evals":[],"jsHooks":[]}</script>
```


```r
out$plot -> PCoA_tmp2_fam

PCoA_tmp2 + 
  theme(legend.position = "none") -> PCoA_tmp2_noleg

PCoA_tmp2_noleg
```

![](NRP72_draft_16S_Figs_beta_PCoA_files/figure-html/unnamed-chunk-56-1.png)<!-- -->

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
<div class="datatables html-widget html-fill-item-overflow-hidden html-fill-item" id="htmlwidget-530a78b6872d546cc16a" style="width:100%;height:auto;"></div>
<script type="application/json" data-for="htmlwidget-530a78b6872d546cc16a">{"x":{"filter":"none","vertical":false,"data":[["1","2","3","4"],["Formiat_mM","Acetat_mM","Valerat_mM","Total_SCFA_mM"],[-0.613360599591134,-0.127975964263856,0.178226366838872,-0.199371069496966],[0.042165219774713,-0.372980394645997,-0.677351910348987,-0.424683773969117],["Formiat_mM","Acetat_mM","Valerat_mM","Total_SCFA_mM"],[0.001,0.04,0.001,0.007],[0.377989130889445,0.155492222219548,0.490570248290006,0.220105131225016]],"container":"<table class=\"display\">\n  <thead>\n    <tr>\n      <th> <\/th>\n      <th>label<\/th>\n      <th>Axis.1<\/th>\n      <th>Axis.2<\/th>\n      <th>id<\/th>\n      <th>pval<\/th>\n      <th>r<\/th>\n    <\/tr>\n  <\/thead>\n<\/table>","options":{"columnDefs":[{"className":"dt-right","targets":[2,3,5,6]},{"orderable":false,"targets":0}],"order":[],"autoWidth":false,"orderClasses":false}},"evals":[],"jsHooks":[]}</script>
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

![](NRP72_draft_16S_Figs_beta_PCoA_files/figure-html/unnamed-chunk-59-1.png)<!-- -->


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
<div class="datatables html-widget html-fill-item-overflow-hidden html-fill-item" id="htmlwidget-4cadc5f8c7ab9ff923ac" style="width:100%;height:auto;"></div>
<script type="application/json" data-for="htmlwidget-4cadc5f8c7ab9ff923ac">{"x":{"filter":"none","vertical":false,"data":[["ASV0012","ASV0004","ASV0010","ASV0052","ASV0030","ASV0035","ASV0106","ASV0275"],[0.419469264716501,-0.542610944831473,-0.668057142724215,0.634491810605059,0.0923438525333804,-0.475786979717028,-0.235980253754372,0.716607484689898],[0.41313443709802,0.20743545326037,-0.213247478016456,-0.157224913125192,0.549578877374205,0.217200882303641,-0.59337501428833,-0.0384851666690667],[0.001,0.003,0.001,0.001,0.011,0.003,0.001,0.001],[0.3466345271581,0.337456104720238,0.491774832825221,0.427299531032111,0.310564329556598,0.273549473341732,0.407780587743654,0.515007395167129],["Bacteria","Bacteria","Bacteria","Bacteria","Bacteria","Bacteria","Bacteria","Bacteria"],["Firmicutes","Proteobacteria","Firmicutes","Firmicutes","Bacteroidota","Proteobacteria","Firmicutes","Firmicutes"],["Clostridia","Gammaproteobacteria","Clostridia","Negativicutes","Bacteroidia","Gammaproteobacteria","Clostridia","Bacilli"],["Oscillospirales","Enterobacterales","Lachnospirales","Veillonellales-Selenomonadales","Bacteroidales","Burkholderiales","Oscillospirales","Erysipelotrichales"],["Ruminococcaceae","Enterobacteriaceae","Lachnospiraceae","Veillonellaceae","Prevotellaceae","Sutterellaceae","Oscillospiraceae","Erysipelotrichaceae"],[null,null,null,null,null,null,null,null],[null,null,null,null,null,null,null,null],[null,null,null,null,null,null,null,null],[0.0104,0.0222857142857143,0.0104,0.0104,0.0476666666666667,0.0222857142857143,0.0104,0.0104]],"container":"<table class=\"display\">\n  <thead>\n    <tr>\n      <th> <\/th>\n      <th>Axis.1<\/th>\n      <th>Axis.2<\/th>\n      <th>pval<\/th>\n      <th>r<\/th>\n      <th>Kingdom<\/th>\n      <th>Phylum<\/th>\n      <th>Class<\/th>\n      <th>Order<\/th>\n      <th>tax_rank_plot<\/th>\n      <th>Genus<\/th>\n      <th>Species<\/th>\n      <th>Strain<\/th>\n      <th>pval.adj<\/th>\n    <\/tr>\n  <\/thead>\n<\/table>","options":{"columnDefs":[{"className":"dt-right","targets":[1,2,3,4,13]},{"orderable":false,"targets":0}],"order":[],"autoWidth":false,"orderClasses":false}},"evals":[],"jsHooks":[]}</script>
```


```r
out$plot -> PCoA_tmp2_fam

PCoA_tmp2 + 
  theme(legend.position = "none") -> PCoA_tmp2_noleg

PCoA_tmp2_noleg
```

![](NRP72_draft_16S_Figs_beta_PCoA_files/figure-html/unnamed-chunk-61-1.png)<!-- -->

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
<div class="datatables html-widget html-fill-item-overflow-hidden html-fill-item" id="htmlwidget-fc8df40fd67c2027a709" style="width:100%;height:auto;"></div>
<script type="application/json" data-for="htmlwidget-fc8df40fd67c2027a709">{"x":{"filter":"none","vertical":false,"data":[["1","2","3"],["Formiat_mM","Valerat_mM","Total_SCFA_mM"],[-0.600973704253337,0.185820508581589,-0.18149770546043],[0.196394594911277,-0.579040546718349,-0.332598322258841],["Formiat_mM","Valerat_mM","Total_SCFA_mM"],[0.001,0.001,0.04],[0.399740230114341,0.369817216153405,0.143563061056797]],"container":"<table class=\"display\">\n  <thead>\n    <tr>\n      <th> <\/th>\n      <th>label<\/th>\n      <th>Axis.1<\/th>\n      <th>Axis.2<\/th>\n      <th>id<\/th>\n      <th>pval<\/th>\n      <th>r<\/th>\n    <\/tr>\n  <\/thead>\n<\/table>","options":{"columnDefs":[{"className":"dt-right","targets":[2,3,5,6]},{"orderable":false,"targets":0}],"order":[],"autoWidth":false,"orderClasses":false}},"evals":[],"jsHooks":[]}</script>
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
<div class="datatables html-widget html-fill-item-overflow-hidden html-fill-item" id="htmlwidget-d8c7d17cf4080686b93c" style="width:100%;height:auto;"></div>
<script type="application/json" data-for="htmlwidget-d8c7d17cf4080686b93c">{"x":{"filter":"none","vertical":false,"data":[["1","2","3","4","5","6"],["bjaccard","bjaccard","bjaccard","aichinson","aichinson","aichinson"],["Treatment","Residual","Total","Treatment","Residual","Total"],[3,13,16,3,13,16],[0.193,0.413,0.605,0.287,0.575,0.862],[0.318,0.682,1,0.333,0.667,1],[2.022,null,null,2.166,null,null],[0.001,null,null,0.001,null,null]],"container":"<table class=\"display\">\n  <thead>\n    <tr>\n      <th> <\/th>\n      <th>Distance<\/th>\n      <th>terms<\/th>\n      <th>Df<\/th>\n      <th>SumOfSqs<\/th>\n      <th>R2<\/th>\n      <th>F<\/th>\n      <th>Pr(&gt;F)<\/th>\n    <\/tr>\n  <\/thead>\n<\/table>","options":{"columnDefs":[{"className":"dt-right","targets":[3,4,5,6,7]},{"orderable":false,"targets":0}],"order":[],"autoWidth":false,"orderClasses":false}},"evals":[],"jsHooks":[]}</script>
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
<div class="datatables html-widget html-fill-item-overflow-hidden html-fill-item" id="htmlwidget-fbc9973e359619ea20b8" style="width:100%;height:auto;"></div>
<script type="application/json" data-for="htmlwidget-fbc9973e359619ea20b8">{"x":{"filter":"none","vertical":false,"data":[["1","2","3","4","5","6","7","8"],["bjaccard","bjaccard","bjaccard","bjaccard","aichinson","aichinson","aichinson","aichinson"],["Treatment","Day_of_Treatment","Residual","Total","Treatment","Day_of_Treatment","Residual","Total"],[3,1,12,16,3,1,12,16],[0.189,0.051,0.362,0.605,0.284,0.063,0.512,0.862],[0.312,0.084,0.598,1,0.329,0.073,0.593,1],[2.087,1.693,null,null,2.218,1.483,null,null],[0.002,0.074,null,null,0.001,0.087,null,null]],"container":"<table class=\"display\">\n  <thead>\n    <tr>\n      <th> <\/th>\n      <th>Distance<\/th>\n      <th>terms<\/th>\n      <th>Df<\/th>\n      <th>SumOfSqs<\/th>\n      <th>R2<\/th>\n      <th>F<\/th>\n      <th>Pr(&gt;F)<\/th>\n    <\/tr>\n  <\/thead>\n<\/table>","options":{"columnDefs":[{"className":"dt-right","targets":[3,4,5,6,7]},{"orderable":false,"targets":0}],"order":[],"autoWidth":false,"orderClasses":false}},"evals":[],"jsHooks":[]}</script>
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
<div class="datatables html-widget html-fill-item-overflow-hidden html-fill-item" id="htmlwidget-fc9c72093435cc6a0a18" style="width:100%;height:auto;"></div>
<script type="application/json" data-for="htmlwidget-fc9c72093435cc6a0a18">{"x":{"filter":"none","vertical":false,"data":[["1","2","3","4","5","6","7","8","9","10","11","12"],["bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson"],["UNTREATED","UNTREATED","UNTREATED","CTX","CTX","CTX+HV292.1","UNTREATED","UNTREATED","UNTREATED","CTX","CTX","CTX+HV292.1"],["CTX","CTX+HV292.1","HV292.1","CTX+HV292.1","HV292.1","HV292.1","CTX","CTX+HV292.1","HV292.1","CTX+HV292.1","HV292.1","HV292.1"],[0.262,0.184,0.386,0.139,0.277,0.29,0.279,0.242,0.428,0.126,0.284,0.274],[0.04,0.093,0.1,0.106,0.043,0.01,0.019,0.03,0.1,0.11,0.03,0.02],[0.12,0.14,0.12,0.106,0.086,0.06,0.114,0.051,0.12,0.11,0.051,0.06]],"container":"<table class=\"display\">\n  <thead>\n    <tr>\n      <th> <\/th>\n      <th>Distance<\/th>\n      <th>X1<\/th>\n      <th>X2<\/th>\n      <th>R2<\/th>\n      <th>pval<\/th>\n      <th>pvalFDR<\/th>\n    <\/tr>\n  <\/thead>\n<\/table>","options":{"columnDefs":[{"className":"dt-right","targets":[4,5,6]},{"orderable":false,"targets":0}],"order":[],"autoWidth":false,"orderClasses":false}},"evals":[],"jsHooks":[]}</script>
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
<div class="datatables html-widget html-fill-item-overflow-hidden html-fill-item" id="htmlwidget-33198698d6599fd3fc97" style="width:100%;height:auto;"></div>
<script type="application/json" data-for="htmlwidget-33198698d6599fd3fc97">{"x":{"filter":"none","vertical":false,"data":[["1","2","3","4","5","6","7","8","9","10","11","12"],["bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson"],[["UNTREATED"],["UNTREATED"],["UNTREATED"],["CTX"],["CTX"],["CTX+HV292.1"],["UNTREATED"],["UNTREATED"],["UNTREATED"],["CTX"],["CTX"],["CTX+HV292.1"]],[["CTX"],["CTX+HV292.1"],["HV292.1"],["CTX+HV292.1"],["HV292.1"],["HV292.1"],["CTX"],["CTX+HV292.1"],["HV292.1"],["CTX+HV292.1"],["HV292.1"],["HV292.1"]],[[3],[3],[3],[5],[5],[6],[3],[3],[3],[5],[5],[6]],[[5],[6],[3],[6],[3],[3],[5],[6],[3],[6],[3],[3]],[[0.022],[0.13],[0.089],[0.127],[0.017],[0.018],[0.016],[0.037],[0.096],[0.147],[0.021],[0.016]],[[2.14065763049362],[1.43275641249348],[2.51415311555688],[1.41784644473462],[3.23293227343592],[4.15591028387124],[2.52503979302175],[2.26840481630874],[2.99061286086113],[1.26579096941619],[3.0601681404366],[3.38342177726307]],[[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999]],[0.044,0.13,0.13,0.13,0.044,0.044,0.042,0.055,0.115,0.147,0.042,0.042]],"container":"<table class=\"display\">\n  <thead>\n    <tr>\n      <th> <\/th>\n      <th>Distance<\/th>\n      <th>Level1<\/th>\n      <th>Level2<\/th>\n      <th>N1<\/th>\n      <th>N2<\/th>\n      <th>p.value<\/th>\n      <th>tw2.stat<\/th>\n      <th>nrep<\/th>\n      <th>pvalFDR<\/th>\n    <\/tr>\n  <\/thead>\n<\/table>","options":{"columnDefs":[{"className":"dt-right","targets":9},{"orderable":false,"targets":0}],"order":[],"autoWidth":false,"orderClasses":false}},"evals":[],"jsHooks":[]}</script>
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
<div class="datatables html-widget html-fill-item-overflow-hidden html-fill-item" id="htmlwidget-38c86864160d6d855a80" style="width:100%;height:auto;"></div>
<script type="application/json" data-for="htmlwidget-38c86864160d6d855a80">{"x":{"filter":"none","vertical":false,"data":[["1","2","3","4","5","6"],["bjaccard","bjaccard","bjaccard","aichinson","aichinson","aichinson"],["Treatment","Residual","Total","Treatment","Residual","Total"],[3,19,22,3,19,22],[0.549,1.455,2.004,0.543,1.396,1.939],[0.274,0.726,1,0.28,0.72,1],[2.39,null,null,2.462,null,null],[0.01,null,null,0.002,null,null]],"container":"<table class=\"display\">\n  <thead>\n    <tr>\n      <th> <\/th>\n      <th>Distance<\/th>\n      <th>terms<\/th>\n      <th>Df<\/th>\n      <th>SumOfSqs<\/th>\n      <th>R2<\/th>\n      <th>F<\/th>\n      <th>Pr(&gt;F)<\/th>\n    <\/tr>\n  <\/thead>\n<\/table>","options":{"columnDefs":[{"className":"dt-right","targets":[3,4,5,6,7]},{"orderable":false,"targets":0}],"order":[],"autoWidth":false,"orderClasses":false}},"evals":[],"jsHooks":[]}</script>
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
<div class="datatables html-widget html-fill-item-overflow-hidden html-fill-item" id="htmlwidget-d96970b9501c191ad406" style="width:100%;height:auto;"></div>
<script type="application/json" data-for="htmlwidget-d96970b9501c191ad406">{"x":{"filter":"none","vertical":false,"data":[["1","2","3","4","5","6","7","8"],["bjaccard","bjaccard","bjaccard","bjaccard","aichinson","aichinson","aichinson","aichinson"],["Treatment","Day_of_Treatment","Residual","Total","Treatment","Day_of_Treatment","Residual","Total"],[3,1,18,22,3,1,18,22],[0.549,0.176,1.279,2.004,0.541,0.154,1.241,1.939],[0.274,0.088,0.638,1,0.279,0.08,0.64,1],[2.575,2.47,null,null,2.615,2.24,null,null],[0.003,0.016,null,null,0.001,0.014,null,null]],"container":"<table class=\"display\">\n  <thead>\n    <tr>\n      <th> <\/th>\n      <th>Distance<\/th>\n      <th>terms<\/th>\n      <th>Df<\/th>\n      <th>SumOfSqs<\/th>\n      <th>R2<\/th>\n      <th>F<\/th>\n      <th>Pr(&gt;F)<\/th>\n    <\/tr>\n  <\/thead>\n<\/table>","options":{"columnDefs":[{"className":"dt-right","targets":[3,4,5,6,7]},{"orderable":false,"targets":0}],"order":[],"autoWidth":false,"orderClasses":false}},"evals":[],"jsHooks":[]}</script>
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
<div class="datatables html-widget html-fill-item-overflow-hidden html-fill-item" id="htmlwidget-bbf49928d385f1dd9e5c" style="width:100%;height:auto;"></div>
<script type="application/json" data-for="htmlwidget-bbf49928d385f1dd9e5c">{"x":{"filter":"none","vertical":false,"data":[["1","2","3","4","5","6","7","8","9","10","11","12"],["bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson"],["UNTREATED","UNTREATED","UNTREATED","CTX","CTX","CTX+HV292.1","UNTREATED","UNTREATED","UNTREATED","CTX","CTX","CTX+HV292.1"],["CTX","CTX+HV292.1","HV292.1","CTX+HV292.1","HV292.1","HV292.1","CTX","CTX+HV292.1","HV292.1","CTX+HV292.1","HV292.1","HV292.1"],[0.24,0.339,0.202,0.068,0.138,0.238,0.249,0.307,0.251,0.067,0.164,0.234],[0.002,0.003,0.6,0.233,0.279,0.107,0.005,0.004,0.6,0.249,0.187,0.103],[0.012,0.009,0.6,0.35,0.335,0.214,0.015,0.024,0.6,0.299,0.28,0.206]],"container":"<table class=\"display\">\n  <thead>\n    <tr>\n      <th> <\/th>\n      <th>Distance<\/th>\n      <th>X1<\/th>\n      <th>X2<\/th>\n      <th>R2<\/th>\n      <th>pval<\/th>\n      <th>pvalFDR<\/th>\n    <\/tr>\n  <\/thead>\n<\/table>","options":{"columnDefs":[{"className":"dt-right","targets":[4,5,6]},{"orderable":false,"targets":0}],"order":[],"autoWidth":false,"orderClasses":false}},"evals":[],"jsHooks":[]}</script>
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
<div class="datatables html-widget html-fill-item-overflow-hidden html-fill-item" id="htmlwidget-90bfa43893a2d8eda041" style="width:100%;height:auto;"></div>
<script type="application/json" data-for="htmlwidget-90bfa43893a2d8eda041">{"x":{"filter":"none","vertical":false,"data":[["1","2","3","4","5","6","7","8","9","10","11","12"],["bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson"],[["UNTREATED"],["UNTREATED"],["UNTREATED"],["CTX"],["CTX"],["CTX+HV292.1"],["UNTREATED"],["UNTREATED"],["UNTREATED"],["CTX"],["CTX"],["CTX+HV292.1"]],[["CTX"],["CTX+HV292.1"],["HV292.1"],["CTX+HV292.1"],["HV292.1"],["HV292.1"],["CTX"],["CTX+HV292.1"],["HV292.1"],["CTX+HV292.1"],["HV292.1"],["HV292.1"]],[[4],[4],[4],[10],[10],[8],[4],[4],[4],[10],[10],[8]],[[10],[8],[1],[8],[1],[1],[10],[8],[1],[8],[1],[1]],[[0.002],[0.004],[null],[0.217],[null],[null],[0.002],[0.003],[null],[0.237],[null],[null]],[[5.657783359775],[5.92474408517821],[null],[1.24902845777323],[null],[null],[4.80208392802225],[4.73316227298532],[null],[1.17510680306781],[null],[null]],[[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999]],[0.006,0.006,null,0.217,null,null,0.005,0.005,null,0.237,null,null]],"container":"<table class=\"display\">\n  <thead>\n    <tr>\n      <th> <\/th>\n      <th>Distance<\/th>\n      <th>Level1<\/th>\n      <th>Level2<\/th>\n      <th>N1<\/th>\n      <th>N2<\/th>\n      <th>p.value<\/th>\n      <th>tw2.stat<\/th>\n      <th>nrep<\/th>\n      <th>pvalFDR<\/th>\n    <\/tr>\n  <\/thead>\n<\/table>","options":{"columnDefs":[{"className":"dt-right","targets":9},{"orderable":false,"targets":0}],"order":[],"autoWidth":false,"orderClasses":false}},"evals":[],"jsHooks":[]}</script>
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
<div class="datatables html-widget html-fill-item-overflow-hidden html-fill-item" id="htmlwidget-0349e08474188736d854" style="width:100%;height:auto;"></div>
<script type="application/json" data-for="htmlwidget-0349e08474188736d854">{"x":{"filter":"none","vertical":false,"data":[["1","2","3","4","5","6"],["bjaccard","bjaccard","bjaccard","aichinson","aichinson","aichinson"],["Treatment","Residual","Total","Treatment","Residual","Total"],[3,11,14,3,11,14],[0.413,0.611,1.023,0.463,0.702,1.165],[0.403,0.597,1,0.397,0.603,1],[2.48,null,null,2.417,null,null],[0.002,null,null,0.001,null,null]],"container":"<table class=\"display\">\n  <thead>\n    <tr>\n      <th> <\/th>\n      <th>Distance<\/th>\n      <th>terms<\/th>\n      <th>Df<\/th>\n      <th>SumOfSqs<\/th>\n      <th>R2<\/th>\n      <th>F<\/th>\n      <th>Pr(&gt;F)<\/th>\n    <\/tr>\n  <\/thead>\n<\/table>","options":{"columnDefs":[{"className":"dt-right","targets":[3,4,5,6,7]},{"orderable":false,"targets":0}],"order":[],"autoWidth":false,"orderClasses":false}},"evals":[],"jsHooks":[]}</script>
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
<div class="datatables html-widget html-fill-item-overflow-hidden html-fill-item" id="htmlwidget-2c91a99e85c6bc68af90" style="width:100%;height:auto;"></div>
<script type="application/json" data-for="htmlwidget-2c91a99e85c6bc68af90">{"x":{"filter":"none","vertical":false,"data":[["1","2","3","4","5","6","7","8"],["bjaccard","bjaccard","bjaccard","bjaccard","aichinson","aichinson","aichinson","aichinson"],["Treatment","Day_of_Treatment","Residual","Total","Treatment","Day_of_Treatment","Residual","Total"],[3,1,10,14,3,1,10,14],[0.409,0.092,0.518,1.023,0.457,0.092,0.61,1.165],[0.4,0.09,0.506,1,0.392,0.079,0.523,1],[2.633,1.779,null,null,2.498,1.516,null,null],[0.001,0.096,null,null,0.001,0.127,null,null]],"container":"<table class=\"display\">\n  <thead>\n    <tr>\n      <th> <\/th>\n      <th>Distance<\/th>\n      <th>terms<\/th>\n      <th>Df<\/th>\n      <th>SumOfSqs<\/th>\n      <th>R2<\/th>\n      <th>F<\/th>\n      <th>Pr(&gt;F)<\/th>\n    <\/tr>\n  <\/thead>\n<\/table>","options":{"columnDefs":[{"className":"dt-right","targets":[3,4,5,6,7]},{"orderable":false,"targets":0}],"order":[],"autoWidth":false,"orderClasses":false}},"evals":[],"jsHooks":[]}</script>
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
```

```
## 'nperm' >= set of all permutations: complete enumeration.
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

```r
tmp_df %>% 
  DT::datatable()
```

```{=html}
<div class="datatables html-widget html-fill-item-overflow-hidden html-fill-item" id="htmlwidget-f4acdbbe3467343ec6d3" style="width:100%;height:auto;"></div>
<script type="application/json" data-for="htmlwidget-f4acdbbe3467343ec6d3">{"x":{"filter":"none","vertical":false,"data":[["1","2","3","4","5","6","7","8","9","10","11","12"],["bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson"],["UNTREATED","UNTREATED","UNTREATED","CTX","CTX","CTX+HV292.1","UNTREATED","UNTREATED","UNTREATED","CTX","CTX","CTX+HV292.1"],["CTX","CTX+HV292.1","HV292.1","CTX+HV292.1","HV292.1","HV292.1","CTX","CTX+HV292.1","HV292.1","CTX+HV292.1","HV292.1","HV292.1"],[0.359,0.384,0.238,0.13,0.325,0.347,0.339,0.353,0.278,0.148,0.332,0.346],[0.015,0.015,1,0.204,0.296,0.167,0.01,0.017,0.75,0.11,0.267,0.167],[0.06,0.06,1,0.306,0.355,0.333,0.06,0.051,0.75,0.22,0.32,0.25]],"container":"<table class=\"display\">\n  <thead>\n    <tr>\n      <th> <\/th>\n      <th>Distance<\/th>\n      <th>X1<\/th>\n      <th>X2<\/th>\n      <th>R2<\/th>\n      <th>pval<\/th>\n      <th>pvalFDR<\/th>\n    <\/tr>\n  <\/thead>\n<\/table>","options":{"columnDefs":[{"className":"dt-right","targets":[4,5,6]},{"orderable":false,"targets":0}],"order":[],"autoWidth":false,"orderClasses":false}},"evals":[],"jsHooks":[]}</script>
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
<div class="datatables html-widget html-fill-item-overflow-hidden html-fill-item" id="htmlwidget-d15c90e001b0818c17b0" style="width:100%;height:auto;"></div>
<script type="application/json" data-for="htmlwidget-d15c90e001b0818c17b0">{"x":{"filter":"none","vertical":false,"data":[["1","2","3","4","5","6","7","8","9","10","11","12"],["bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson"],[["UNTREATED"],["UNTREATED"],["UNTREATED"],["CTX"],["CTX"],["CTX+HV292.1"],["UNTREATED"],["UNTREATED"],["UNTREATED"],["CTX"],["CTX"],["CTX+HV292.1"]],[["CTX"],["CTX+HV292.1"],["HV292.1"],["CTX+HV292.1"],["HV292.1"],["HV292.1"],["CTX"],["CTX+HV292.1"],["HV292.1"],["CTX+HV292.1"],["HV292.1"],["HV292.1"]],[[3],[3],[3],[6],[6],[5],[3],[3],[3],[6],[6],[5]],[[6],[5],[1],[5],[1],[1],[6],[5],[1],[5],[1],[1]],[[0.014],[0.014],[null],[0.204],[null],[null],[0.014],[0.023],[null],[0.107],[null],[null]],[[4.15085636494982],[4.05961365172066],[null],[1.32354820000393],[null],[null],[3.46381163237091],[3.29794275372578],[null],[1.54767557636096],[null],[null]],[[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999]],[0.021,0.021,null,0.204,null,null,0.034,0.034,null,0.107,null,null]],"container":"<table class=\"display\">\n  <thead>\n    <tr>\n      <th> <\/th>\n      <th>Distance<\/th>\n      <th>Level1<\/th>\n      <th>Level2<\/th>\n      <th>N1<\/th>\n      <th>N2<\/th>\n      <th>p.value<\/th>\n      <th>tw2.stat<\/th>\n      <th>nrep<\/th>\n      <th>pvalFDR<\/th>\n    <\/tr>\n  <\/thead>\n<\/table>","options":{"columnDefs":[{"className":"dt-right","targets":9},{"orderable":false,"targets":0}],"order":[],"autoWidth":false,"orderClasses":false}},"evals":[],"jsHooks":[]}</script>
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

![](NRP72_draft_16S_Figs_beta_PCoA_files/figure-html/unnamed-chunk-81-1.png)<!-- -->


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
<div class="datatables html-widget html-fill-item-overflow-hidden html-fill-item" id="htmlwidget-fec2d181a4be4b780ca9" style="width:100%;height:auto;"></div>
<script type="application/json" data-for="htmlwidget-fec2d181a4be4b780ca9">{"x":{"filter":"none","vertical":false,"data":[["ASV0002","ASV0003","ASV0020","ASV0040","ASV0032","ASV0085","ASV0124","ASV0129"],[0.525535181023166,0.920169096922261,-0.358554123676786,-0.116195913277726,-0.385498361508464,0.61223295617904,0.738008656230995,0.681278011291429],[-0.618870077158341,-0.128496671144786,0.82065886361688,0.684754637802372,0.610396402616658,-0.495623606411118,0.253095005062396,0.318522455886272],[0.001,0.001,0.001,0.001,0.001,0.001,0.001,0.001],[0.659187398895023,0.863222561426021,0.802042030038578,0.482390404254302,0.521192755053068,0.620471951863689,0.608713858259413,0.565596283573027],["Bacteria","Bacteria","Bacteria","Bacteria","Bacteria","Bacteria","Bacteria","Bacteria"],["Bacteroidota","Firmicutes","Firmicutes","Firmicutes","Proteobacteria","Firmicutes","Firmicutes","Actinobacteriota"],["Bacteroidia","Clostridia","Bacilli","Negativicutes","Gammaproteobacteria","Clostridia","Clostridia","Coriobacteriia"],["Bacteroidales","Oscillospirales","Lactobacillales","Veillonellales-Selenomonadales","Burkholderiales","Oscillospirales","Monoglobales","Coriobacteriales"],["Bacteroidaceae","Ruminococcaceae","Lactobacillaceae","Veillonellaceae","Sutterellaceae","Oscillospiraceae","Monoglobaceae","Coriobacteriaceae"],[null,null,null,null,null,null,null,null],[null,null,null,null,null,null,null,null],[null,null,null,null,null,null,null,null],[0.0026,0.0026,0.0026,0.0026,0.0026,0.0026,0.0026,0.0026]],"container":"<table class=\"display\">\n  <thead>\n    <tr>\n      <th> <\/th>\n      <th>Axis.1<\/th>\n      <th>Axis.2<\/th>\n      <th>pval<\/th>\n      <th>r<\/th>\n      <th>Kingdom<\/th>\n      <th>Phylum<\/th>\n      <th>Class<\/th>\n      <th>Order<\/th>\n      <th>tax_rank_plot<\/th>\n      <th>Genus<\/th>\n      <th>Species<\/th>\n      <th>Strain<\/th>\n      <th>pval.adj<\/th>\n    <\/tr>\n  <\/thead>\n<\/table>","options":{"columnDefs":[{"className":"dt-right","targets":[1,2,3,4,13]},{"orderable":false,"targets":0}],"order":[],"autoWidth":false,"orderClasses":false}},"evals":[],"jsHooks":[]}</script>
```


```r
out$plot -> PCoA_tmp2_fam

PCoA_tmp2 + 
  theme(legend.position = "none") -> PCoA_tmp2_noleg

PCoA_tmp2_noleg
```

![](NRP72_draft_16S_Figs_beta_PCoA_files/figure-html/unnamed-chunk-83-1.png)<!-- -->

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
<div class="datatables html-widget html-fill-item-overflow-hidden html-fill-item" id="htmlwidget-bcc569e9a5fb6f3f4fe2" style="width:100%;height:auto;"></div>
<script type="application/json" data-for="htmlwidget-bcc569e9a5fb6f3f4fe2">{"x":{"filter":"none","vertical":false,"data":[["1","2","3","4","5","6"],["Formiat_mM","Acetat_mM","Propionat_mM","Butyrat_mM","Valerat_mM","Total_SCFA_mM"],[-0.549090636476436,0.916195874386532,0.884304826519494,0.950599581786949,0.883991857526509,0.915344362733472],[0.447248992918467,-0.211750385880216,-0.0384601364958492,-0.168139704173301,0.291267990063294,-0.302627172723072],["Formiat_mM","Acetat_mM","Propionat_mM","Butyrat_mM","Valerat_mM","Total_SCFA_mM"],[0.001,0.001,0.001,0.001,0.001,0.001],[0.501532188732681,0.884253106163322,0.783474208304952,0.931910525013008,0.866278646208678,0.929438508058306]],"container":"<table class=\"display\">\n  <thead>\n    <tr>\n      <th> <\/th>\n      <th>label<\/th>\n      <th>Axis.1<\/th>\n      <th>Axis.2<\/th>\n      <th>id<\/th>\n      <th>pval<\/th>\n      <th>r<\/th>\n    <\/tr>\n  <\/thead>\n<\/table>","options":{"columnDefs":[{"className":"dt-right","targets":[2,3,5,6]},{"orderable":false,"targets":0}],"order":[],"autoWidth":false,"orderClasses":false}},"evals":[],"jsHooks":[]}</script>
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

![](NRP72_draft_16S_Figs_beta_PCoA_files/figure-html/unnamed-chunk-86-1.png)<!-- -->


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
<div class="datatables html-widget html-fill-item-overflow-hidden html-fill-item" id="htmlwidget-23dada22fc5642de2fc1" style="width:100%;height:auto;"></div>
<script type="application/json" data-for="htmlwidget-23dada22fc5642de2fc1">{"x":{"filter":"none","vertical":false,"data":[["ASV0002","ASV0003","ASV0020","ASV0040","ASV0032","ASV0085","ASV0124","ASV0129"],[0.656938471951339,0.920440483900327,-0.410941714999905,-0.14447751323852,-0.439146111340108,0.741210046290733,0.722410876054175,0.676490209988427],[-0.409232258728033,0.0417566788403603,0.798012284558517,0.770944906213352,0.605987654237452,-0.207815644309558,0.208667051824549,0.259292335598343],[0.001,0.001,0.001,0.001,0.001,0.001,0.001,0.001],[0.599039197513408,0.848954304630445,0.805696699433367,0.6152298002479,0.560070344193349,0.592579674742107,0.56541941235851,0.524871519510229],["Bacteria","Bacteria","Bacteria","Bacteria","Bacteria","Bacteria","Bacteria","Bacteria"],["Bacteroidota","Firmicutes","Firmicutes","Firmicutes","Proteobacteria","Firmicutes","Firmicutes","Actinobacteriota"],["Bacteroidia","Clostridia","Bacilli","Negativicutes","Gammaproteobacteria","Clostridia","Clostridia","Coriobacteriia"],["Bacteroidales","Oscillospirales","Lactobacillales","Veillonellales-Selenomonadales","Burkholderiales","Oscillospirales","Monoglobales","Coriobacteriales"],["Bacteroidaceae","Ruminococcaceae","Lactobacillaceae","Veillonellaceae","Sutterellaceae","Oscillospiraceae","Monoglobaceae","Coriobacteriaceae"],[null,null,null,null,null,null,null,null],[null,null,null,null,null,null,null,null],[null,null,null,null,null,null,null,null],[0.0026,0.0026,0.0026,0.0026,0.0026,0.0026,0.0026,0.0026]],"container":"<table class=\"display\">\n  <thead>\n    <tr>\n      <th> <\/th>\n      <th>Axis.1<\/th>\n      <th>Axis.2<\/th>\n      <th>pval<\/th>\n      <th>r<\/th>\n      <th>Kingdom<\/th>\n      <th>Phylum<\/th>\n      <th>Class<\/th>\n      <th>Order<\/th>\n      <th>tax_rank_plot<\/th>\n      <th>Genus<\/th>\n      <th>Species<\/th>\n      <th>Strain<\/th>\n      <th>pval.adj<\/th>\n    <\/tr>\n  <\/thead>\n<\/table>","options":{"columnDefs":[{"className":"dt-right","targets":[1,2,3,4,13]},{"orderable":false,"targets":0}],"order":[],"autoWidth":false,"orderClasses":false}},"evals":[],"jsHooks":[]}</script>
```


```r
out$plot -> PCoA_tmp2_fam

PCoA_tmp2 + 
  theme(legend.position = "none") -> PCoA_tmp2_noleg

PCoA_tmp2_noleg
```

![](NRP72_draft_16S_Figs_beta_PCoA_files/figure-html/unnamed-chunk-88-1.png)<!-- -->

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
<div class="datatables html-widget html-fill-item-overflow-hidden html-fill-item" id="htmlwidget-6e16c2883c26dd1edd48" style="width:100%;height:auto;"></div>
<script type="application/json" data-for="htmlwidget-6e16c2883c26dd1edd48">{"x":{"filter":"none","vertical":false,"data":[["1","2","3","4","5","6"],["Formiat_mM","Acetat_mM","Propionat_mM","Butyrat_mM","Valerat_mM","Total_SCFA_mM"],[-0.605272295297702,0.932524691058002,0.861081805956696,0.95281725848877,0.86910860346474,0.954567915935512],[0.31413075304474,-0.0312399320413977,0.209496296459151,-8.57176905462534e-06,0.252566031283776,-0.0932479278949179],["Formiat_mM","Acetat_mM","Propionat_mM","Butyrat_mM","Valerat_mM","Total_SCFA_mM"],[0.001,0.001,0.001,0.001,0.001,0.001],[0.465032681463404,0.870578232786773,0.785350574779747,0.90786072814753,0.819139364774867,0.919895082190162]],"container":"<table class=\"display\">\n  <thead>\n    <tr>\n      <th> <\/th>\n      <th>label<\/th>\n      <th>Axis.1<\/th>\n      <th>Axis.2<\/th>\n      <th>id<\/th>\n      <th>pval<\/th>\n      <th>r<\/th>\n    <\/tr>\n  <\/thead>\n<\/table>","options":{"columnDefs":[{"className":"dt-right","targets":[2,3,5,6]},{"orderable":false,"targets":0}],"order":[],"autoWidth":false,"orderClasses":false}},"evals":[],"jsHooks":[]}</script>
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
<div class="datatables html-widget html-fill-item-overflow-hidden html-fill-item" id="htmlwidget-23e11579ab285b5c7c2a" style="width:100%;height:auto;"></div>
<script type="application/json" data-for="htmlwidget-23e11579ab285b5c7c2a">{"x":{"filter":"none","vertical":false,"data":[["1","2","3","4","5","6"],["bjaccard","bjaccard","bjaccard","aichinson","aichinson","aichinson"],["Treatment","Residual","Total","Treatment","Residual","Total"],[6,28,34,6,28,34],[0.498,0.83,1.328,0.523,0.923,1.446],[0.375,0.625,1,0.362,0.638,1],[2.803,null,null,2.643,null,null],[0.001,null,null,0.001,null,null]],"container":"<table class=\"display\">\n  <thead>\n    <tr>\n      <th> <\/th>\n      <th>Distance<\/th>\n      <th>terms<\/th>\n      <th>Df<\/th>\n      <th>SumOfSqs<\/th>\n      <th>R2<\/th>\n      <th>F<\/th>\n      <th>Pr(&gt;F)<\/th>\n    <\/tr>\n  <\/thead>\n<\/table>","options":{"columnDefs":[{"className":"dt-right","targets":[3,4,5,6,7]},{"orderable":false,"targets":0}],"order":[],"autoWidth":false,"orderClasses":false}},"evals":[],"jsHooks":[]}</script>
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
<div class="datatables html-widget html-fill-item-overflow-hidden html-fill-item" id="htmlwidget-26eca46a20528b6114d0" style="width:100%;height:auto;"></div>
<script type="application/json" data-for="htmlwidget-26eca46a20528b6114d0">{"x":{"filter":"none","vertical":false,"data":[["1","2","3","4","5","6","7","8"],["bjaccard","bjaccard","bjaccard","bjaccard","aichinson","aichinson","aichinson","aichinson"],["Treatment","Day_of_Treatment","Residual","Total","Treatment","Day_of_Treatment","Residual","Total"],[6,1,27,34,6,1,27,34],[0.498,0.063,0.767,1.328,0.523,0.06,0.863,1.446],[0.375,0.047,0.577,1,0.362,0.042,0.597,1],[2.924,2.216,null,null,2.727,1.89,null,null],[0.001,0.003,null,null,0.001,0.005,null,null]],"container":"<table class=\"display\">\n  <thead>\n    <tr>\n      <th> <\/th>\n      <th>Distance<\/th>\n      <th>terms<\/th>\n      <th>Df<\/th>\n      <th>SumOfSqs<\/th>\n      <th>R2<\/th>\n      <th>F<\/th>\n      <th>Pr(&gt;F)<\/th>\n    <\/tr>\n  <\/thead>\n<\/table>","options":{"columnDefs":[{"className":"dt-right","targets":[3,4,5,6,7]},{"orderable":false,"targets":0}],"order":[],"autoWidth":false,"orderClasses":false}},"evals":[],"jsHooks":[]}</script>
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
<div class="datatables html-widget html-fill-item-overflow-hidden html-fill-item" id="htmlwidget-fe4b5d221dd8b1ae7542" style="width:100%;height:auto;"></div>
<script type="application/json" data-for="htmlwidget-fe4b5d221dd8b1ae7542">{"x":{"filter":"none","vertical":false,"data":[["1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","23","24","25","26","27","28","29","30","31","32","33","34","35","36","37","38","39","40","41","42"],["bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson"],["UNTREATED","UNTREATED","UNTREATED","UNTREATED","UNTREATED","UNTREATED","CTX","CTX","CTX","CTX","CTX","CTX+HV292.1","CTX+HV292.1","CTX+HV292.1","CTX+HV292.1","HV292.1","HV292.1","HV292.1","VAN","VAN","VAN+CCUG59168","UNTREATED","UNTREATED","UNTREATED","UNTREATED","UNTREATED","UNTREATED","CTX","CTX","CTX","CTX","CTX","CTX+HV292.1","CTX+HV292.1","CTX+HV292.1","CTX+HV292.1","HV292.1","HV292.1","HV292.1","VAN","VAN","VAN+CCUG59168"],["CTX","CTX+HV292.1","HV292.1","VAN","VAN+CCUG59168","CCUG59168","CTX+HV292.1","HV292.1","VAN","VAN+CCUG59168","CCUG59168","HV292.1","VAN","VAN+CCUG59168","CCUG59168","VAN","VAN+CCUG59168","CCUG59168","VAN+CCUG59168","CCUG59168","CCUG59168","CTX","CTX+HV292.1","HV292.1","VAN","VAN+CCUG59168","CCUG59168","CTX+HV292.1","HV292.1","VAN","VAN+CCUG59168","CCUG59168","HV292.1","VAN","VAN+CCUG59168","CCUG59168","VAN","VAN+CCUG59168","CCUG59168","VAN+CCUG59168","CCUG59168","CCUG59168"],[0.311,0.221,0.204,0.182,0.198,0.241,0.364,0.373,0.206,0.2,0.387,0.278,0.219,0.201,0.263,0.25,0.254,0.33,0.161,0.299,0.28,0.302,0.221,0.222,0.207,0.212,0.228,0.305,0.338,0.208,0.212,0.312,0.262,0.221,0.208,0.228,0.266,0.258,0.301,0.171,0.262,0.262],[0.011,0.007,0.012,0.023,0.012,0.016,0.011,0.006,0.037,0.008,0.009,0.005,0.014,0.023,0.008,0.004,0.009,0.006,0.125,0.007,0.006,0.011,0.014,0.006,0.009,0.01,0.008,0.01,0.01,0.007,0.009,0.007,0.008,0.008,0.008,0.012,0.01,0.005,0.01,0.046,0.012,0.009],[0.018,0.023,0.017,0.026,0.017,0.02,0.018,0.032,0.039,0.02,0.018,0.052,0.018,0.026,0.02,0.084,0.018,0.032,0.125,0.023,0.032,0.014,0.015,0.063,0.019,0.015,0.026,0.015,0.015,0.042,0.019,0.042,0.026,0.026,0.026,0.014,0.015,0.105,0.015,0.046,0.014,0.019]],"container":"<table class=\"display\">\n  <thead>\n    <tr>\n      <th> <\/th>\n      <th>Distance<\/th>\n      <th>X1<\/th>\n      <th>X2<\/th>\n      <th>R2<\/th>\n      <th>pval<\/th>\n      <th>pvalFDR<\/th>\n    <\/tr>\n  <\/thead>\n<\/table>","options":{"columnDefs":[{"className":"dt-right","targets":[4,5,6]},{"orderable":false,"targets":0}],"order":[],"autoWidth":false,"orderClasses":false}},"evals":[],"jsHooks":[]}</script>
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
<div class="datatables html-widget html-fill-item-overflow-hidden html-fill-item" id="htmlwidget-7278515fd065a9a7a95f" style="width:100%;height:auto;"></div>
<script type="application/json" data-for="htmlwidget-7278515fd065a9a7a95f">{"x":{"filter":"none","vertical":false,"data":[["1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","23","24","25","26","27","28","29","30","31","32","33","34","35","36","37","38","39","40","41","42"],["bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson"],[["UNTREATED"],["UNTREATED"],["UNTREATED"],["UNTREATED"],["UNTREATED"],["UNTREATED"],["CTX"],["CTX"],["CTX"],["CTX"],["CTX"],["CTX+HV292.1"],["CTX+HV292.1"],["CTX+HV292.1"],["CTX+HV292.1"],["HV292.1"],["HV292.1"],["HV292.1"],["VAN"],["VAN"],["VAN+CCUG59168"],["UNTREATED"],["UNTREATED"],["UNTREATED"],["UNTREATED"],["UNTREATED"],["UNTREATED"],["CTX"],["CTX"],["CTX"],["CTX"],["CTX"],["CTX+HV292.1"],["CTX+HV292.1"],["CTX+HV292.1"],["CTX+HV292.1"],["HV292.1"],["HV292.1"],["HV292.1"],["VAN"],["VAN"],["VAN+CCUG59168"]],[["CTX"],["CTX+HV292.1"],["HV292.1"],["VAN"],["VAN+CCUG59168"],["CCUG59168"],["CTX+HV292.1"],["HV292.1"],["VAN"],["VAN+CCUG59168"],["CCUG59168"],["HV292.1"],["VAN"],["VAN+CCUG59168"],["CCUG59168"],["VAN"],["VAN+CCUG59168"],["CCUG59168"],["VAN+CCUG59168"],["CCUG59168"],["CCUG59168"],["CTX"],["CTX+HV292.1"],["HV292.1"],["VAN"],["VAN+CCUG59168"],["CCUG59168"],["CTX+HV292.1"],["HV292.1"],["VAN"],["VAN+CCUG59168"],["CCUG59168"],["HV292.1"],["VAN"],["VAN+CCUG59168"],["CCUG59168"],["VAN"],["VAN+CCUG59168"],["CCUG59168"],["VAN+CCUG59168"],["CCUG59168"],["CCUG59168"]],[[5],[5],[5],[5],[5],[5],[5],[5],[5],[5],[5],[5],[5],[5],[5],[5],[5],[5],[5],[5],[5],[5],[5],[5],[5],[5],[5],[5],[5],[5],[5],[5],[5],[5],[5],[5],[5],[5],[5],[5],[5],[5]],[[5],[5],[5],[5],[5],[5],[5],[5],[5],[5],[5],[5],[5],[5],[5],[5],[5],[5],[5],[5],[5],[5],[5],[5],[5],[5],[5],[5],[5],[5],[5],[5],[5],[5],[5],[5],[5],[5],[5],[5],[5],[5]],[[0.01],[0.006],[0.021],[0.021],[0.011],[0.019],[0.01],[0.01],[0.032],[0.009],[0.012],[0.01],[0.013],[0.028],[0.006],[0.01],[0.01],[0.005],[0.147],[0.01],[0.005],[0.005],[0.008],[0.011],[0.008],[0.008],[0.011],[0.014],[0.011],[0.012],[0.006],[0.01],[0.011],[0.009],[0.007],[0.005],[0.007],[0.009],[0.013],[0.071],[0.008],[0.006]],[[3.61161142649598],[2.27493433021227],[2.05416844599888],[1.78469617172638],[1.97135473957909],[2.54335385472751],[4.57619285805083],[4.75686635157071],[2.07378478252489],[2.00117504297029],[5.05892625722876],[3.07725875660178],[2.2451448601377],[2.01216546977223],[2.85541657237427],[2.66816457358006],[2.7244236272209],[3.9346078338194],[1.53159832300531],[3.41724271596055],[3.11032306927843],[3.46448005798299],[2.27427737398797],[2.28524891228454],[2.08337403062385],[2.1538790758296],[2.36755458146576],[3.51845495766692],[4.07955806969518],[2.10148373749945],[2.14891083287175],[3.63532244289794],[2.83408967522825],[2.2689911483374],[2.09820739733884],[2.36049700500679],[2.89918097592929],[2.78553270047736],[3.44716042997851],[1.64716525867757],[2.84723613648874],[2.8379665076349]],[[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999]],[0.018,0.018,0.025,0.025,0.018,0.025,0.018,0.018,0.034,0.018,0.018,0.018,0.018,0.031,0.018,0.018,0.018,0.018,0.147,0.018,0.018,0.014,0.014,0.014,0.014,0.014,0.014,0.015,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.071,0.014,0.014]],"container":"<table class=\"display\">\n  <thead>\n    <tr>\n      <th> <\/th>\n      <th>Distance<\/th>\n      <th>Level1<\/th>\n      <th>Level2<\/th>\n      <th>N1<\/th>\n      <th>N2<\/th>\n      <th>p.value<\/th>\n      <th>tw2.stat<\/th>\n      <th>nrep<\/th>\n      <th>pvalFDR<\/th>\n    <\/tr>\n  <\/thead>\n<\/table>","options":{"columnDefs":[{"className":"dt-right","targets":9},{"orderable":false,"targets":0}],"order":[],"autoWidth":false,"orderClasses":false}},"evals":[],"jsHooks":[]}</script>
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
<div class="datatables html-widget html-fill-item-overflow-hidden html-fill-item" id="htmlwidget-198c93e593dcd35da8ef" style="width:100%;height:auto;"></div>
<script type="application/json" data-for="htmlwidget-198c93e593dcd35da8ef">{"x":{"filter":"none","vertical":false,"data":[["1","2","3","4","5","6"],["bjaccard","bjaccard","bjaccard","aichinson","aichinson","aichinson"],["Treatment","Residual","Total","Treatment","Residual","Total"],[6,55,61,6,55,61],[6.81,5.555,12.366,4.768,2.74,7.508],[0.551,0.449,1,0.635,0.365,1],[11.237,null,null,15.95,null,null],[0.001,null,null,0.001,null,null]],"container":"<table class=\"display\">\n  <thead>\n    <tr>\n      <th> <\/th>\n      <th>Distance<\/th>\n      <th>terms<\/th>\n      <th>Df<\/th>\n      <th>SumOfSqs<\/th>\n      <th>R2<\/th>\n      <th>F<\/th>\n      <th>Pr(&gt;F)<\/th>\n    <\/tr>\n  <\/thead>\n<\/table>","options":{"columnDefs":[{"className":"dt-right","targets":[3,4,5,6,7]},{"orderable":false,"targets":0}],"order":[],"autoWidth":false,"orderClasses":false}},"evals":[],"jsHooks":[]}</script>
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
<div class="datatables html-widget html-fill-item-overflow-hidden html-fill-item" id="htmlwidget-87a6c950be51e7f26976" style="width:100%;height:auto;"></div>
<script type="application/json" data-for="htmlwidget-87a6c950be51e7f26976">{"x":{"filter":"none","vertical":false,"data":[["1","2","3","4","5","6","7","8"],["bjaccard","bjaccard","bjaccard","bjaccard","aichinson","aichinson","aichinson","aichinson"],["Treatment","Day_of_Treatment","Residual","Total","Treatment","Day_of_Treatment","Residual","Total"],[6,1,54,61,6,1,54,61],[6.807,0.271,5.285,12.366,4.765,0.171,2.569,7.508],[0.551,0.022,0.427,1,0.635,0.023,0.342,1],[11.593,2.765,null,null,16.695,3.603,null,null],[0.001,0.026,null,null,0.001,0.012,null,null]],"container":"<table class=\"display\">\n  <thead>\n    <tr>\n      <th> <\/th>\n      <th>Distance<\/th>\n      <th>terms<\/th>\n      <th>Df<\/th>\n      <th>SumOfSqs<\/th>\n      <th>R2<\/th>\n      <th>F<\/th>\n      <th>Pr(&gt;F)<\/th>\n    <\/tr>\n  <\/thead>\n<\/table>","options":{"columnDefs":[{"className":"dt-right","targets":[3,4,5,6,7]},{"orderable":false,"targets":0}],"order":[],"autoWidth":false,"orderClasses":false}},"evals":[],"jsHooks":[]}</script>
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
<div class="datatables html-widget html-fill-item-overflow-hidden html-fill-item" id="htmlwidget-0a706b8192cc8c1e9d72" style="width:100%;height:auto;"></div>
<script type="application/json" data-for="htmlwidget-0a706b8192cc8c1e9d72">{"x":{"filter":"none","vertical":false,"data":[["1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","23","24","25","26","27","28","29","30","31","32","33","34","35","36","37","38","39","40","41","42"],["bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson"],["UNTREATED","UNTREATED","UNTREATED","UNTREATED","UNTREATED","UNTREATED","CTX","CTX","CTX","CTX","CTX","CTX+HV292.1","CTX+HV292.1","CTX+HV292.1","CTX+HV292.1","HV292.1","HV292.1","HV292.1","VAN","VAN","VAN+CCUG59168","UNTREATED","UNTREATED","UNTREATED","UNTREATED","UNTREATED","UNTREATED","CTX","CTX","CTX","CTX","CTX","CTX+HV292.1","CTX+HV292.1","CTX+HV292.1","CTX+HV292.1","HV292.1","HV292.1","HV292.1","VAN","VAN","VAN+CCUG59168"],["CTX","CTX+HV292.1","HV292.1","VAN","VAN+CCUG59168","CCUG59168","CTX+HV292.1","HV292.1","VAN","VAN+CCUG59168","CCUG59168","HV292.1","VAN","VAN+CCUG59168","CCUG59168","VAN","VAN+CCUG59168","CCUG59168","VAN+CCUG59168","CCUG59168","CCUG59168","CTX","CTX+HV292.1","HV292.1","VAN","VAN+CCUG59168","CCUG59168","CTX+HV292.1","HV292.1","VAN","VAN+CCUG59168","CCUG59168","HV292.1","VAN","VAN+CCUG59168","CCUG59168","VAN","VAN+CCUG59168","CCUG59168","VAN+CCUG59168","CCUG59168","CCUG59168"],[0.533,0.401,0.24,0.561,0.583,0.28,0.161,0.495,0.238,0.267,0.543,0.363,0.285,0.314,0.415,0.519,0.561,0.249,0.127,0.569,0.611,0.641,0.459,0.238,0.649,0.667,0.269,0.186,0.596,0.318,0.352,0.656,0.425,0.34,0.364,0.481,0.607,0.639,0.234,0.185,0.662,0.697],[0.001,0.001,0.001,0.001,0.001,0.001,0.005,0.001,0.001,0.001,0.001,0.001,0.001,0.001,0.001,0.001,0.001,0.001,0.007,0.001,0.001,0.001,0.001,0.001,0.001,0.001,0.001,0.001,0.002,0.001,0.001,0.001,0.001,0.001,0.001,0.001,0.001,0.001,0.001,0.001,0.001,0.001],[0.002,0.002,0.002,0.002,0.002,0.002,0.005,0.002,0.002,0.002,0.002,0.002,0.002,0.002,0.002,0.002,0.002,0.002,0.007,0.002,0.002,0.002,0.002,0.002,0.002,0.002,0.002,0.002,0.002,0.002,0.002,0.002,0.002,0.002,0.002,0.002,0.002,0.002,0.002,0.002,0.002,0.002]],"container":"<table class=\"display\">\n  <thead>\n    <tr>\n      <th> <\/th>\n      <th>Distance<\/th>\n      <th>X1<\/th>\n      <th>X2<\/th>\n      <th>R2<\/th>\n      <th>pval<\/th>\n      <th>pvalFDR<\/th>\n    <\/tr>\n  <\/thead>\n<\/table>","options":{"columnDefs":[{"className":"dt-right","targets":[4,5,6]},{"orderable":false,"targets":0}],"order":[],"autoWidth":false,"orderClasses":false}},"evals":[],"jsHooks":[]}</script>
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
<div class="datatables html-widget html-fill-item-overflow-hidden html-fill-item" id="htmlwidget-155dc90fc4a3b1a0cf9b" style="width:100%;height:auto;"></div>
<script type="application/json" data-for="htmlwidget-155dc90fc4a3b1a0cf9b">{"x":{"filter":"none","vertical":false,"data":[["1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","23","24","25","26","27","28","29","30","31","32","33","34","35","36","37","38","39","40","41","42"],["bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson"],[["UNTREATED"],["UNTREATED"],["UNTREATED"],["UNTREATED"],["UNTREATED"],["UNTREATED"],["CTX"],["CTX"],["CTX"],["CTX"],["CTX"],["CTX+HV292.1"],["CTX+HV292.1"],["CTX+HV292.1"],["CTX+HV292.1"],["HV292.1"],["HV292.1"],["HV292.1"],["VAN"],["VAN"],["VAN+CCUG59168"],["UNTREATED"],["UNTREATED"],["UNTREATED"],["UNTREATED"],["UNTREATED"],["UNTREATED"],["CTX"],["CTX"],["CTX"],["CTX"],["CTX"],["CTX+HV292.1"],["CTX+HV292.1"],["CTX+HV292.1"],["CTX+HV292.1"],["HV292.1"],["HV292.1"],["HV292.1"],["VAN"],["VAN"],["VAN+CCUG59168"]],[["CTX"],["CTX+HV292.1"],["HV292.1"],["VAN"],["VAN+CCUG59168"],["CCUG59168"],["CTX+HV292.1"],["HV292.1"],["VAN"],["VAN+CCUG59168"],["CCUG59168"],["HV292.1"],["VAN"],["VAN+CCUG59168"],["CCUG59168"],["VAN"],["VAN+CCUG59168"],["CCUG59168"],["VAN+CCUG59168"],["CCUG59168"],["CCUG59168"],["CTX"],["CTX+HV292.1"],["HV292.1"],["VAN"],["VAN+CCUG59168"],["CCUG59168"],["CTX+HV292.1"],["HV292.1"],["VAN"],["VAN+CCUG59168"],["CCUG59168"],["HV292.1"],["VAN"],["VAN+CCUG59168"],["CCUG59168"],["VAN"],["VAN+CCUG59168"],["CCUG59168"],["VAN+CCUG59168"],["CCUG59168"],["CCUG59168"]],[[9],[9],[9],[9],[9],[9],[9],[9],[9],[9],[9],[8],[8],[8],[8],[9],[9],[9],[9],[9],[9],[9],[9],[9],[9],[9],[9],[9],[9],[9],[9],[9],[8],[8],[8],[8],[9],[9],[9],[9],[9],[9]],[[9],[8],[9],[9],[9],[9],[8],[9],[9],[9],[9],[9],[9],[9],[9],[9],[9],[9],[9],[9],[9],[9],[8],[9],[9],[9],[9],[8],[9],[9],[9],[9],[9],[9],[9],[9],[9],[9],[9],[9],[9],[9]],[[0.001],[0.001],[0.001],[0.001],[0.001],[0.001],[0.003],[0.001],[0.001],[0.001],[0.001],[0.002],[0.001],[0.001],[0.001],[0.001],[0.001],[0.001],[0.007],[0.001],[0.001],[0.001],[0.001],[0.001],[0.001],[0.001],[0.001],[0.002],[0.001],[0.001],[0.001],[0.001],[0.001],[0.001],[0.001],[0.001],[0.001],[0.001],[0.001],[0.001],[0.001],[0.001]],[[18.2816588874981],[9.15338020543391],[5.06185266549415],[20.4103010341225],[22.4048296343546],[6.21677752428232],[2.83644065641568],[15.7111151619778],[4.99877598266036],[5.84221460105487],[18.9997630666877],[7.94357672569024],[5.88717517972496],[6.69673034336503],[9.63920250182291],[17.2757449393879],[20.4695843497377],[5.30779158939088],[2.32601358146117],[21.0987478982016],[25.1251183946221],[28.5778012204383],[12.0058538560295],[5.00787389919481],[29.6474174068086],[32.0437257789525],[5.89200113960892],[3.2895018858734],[23.5966772050598],[7.44789147346748],[8.68999194480714],[30.5096148098283],[10.5955889781548],[7.42800777989058],[8.13872158473742],[13.0596893305605],[24.6774842945271],[28.3333787305631],[4.89214808657102],[3.62119537870214],[31.2791596132917],[36.832090341497]],[[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999]],[0.001,0.001,0.001,0.001,0.001,0.001,0.003,0.001,0.001,0.001,0.001,0.002,0.001,0.001,0.001,0.001,0.001,0.001,0.007,0.001,0.001,0.001,0.001,0.001,0.001,0.001,0.001,0.002,0.001,0.001,0.001,0.001,0.001,0.001,0.001,0.001,0.001,0.001,0.001,0.001,0.001,0.001]],"container":"<table class=\"display\">\n  <thead>\n    <tr>\n      <th> <\/th>\n      <th>Distance<\/th>\n      <th>Level1<\/th>\n      <th>Level2<\/th>\n      <th>N1<\/th>\n      <th>N2<\/th>\n      <th>p.value<\/th>\n      <th>tw2.stat<\/th>\n      <th>nrep<\/th>\n      <th>pvalFDR<\/th>\n    <\/tr>\n  <\/thead>\n<\/table>","options":{"columnDefs":[{"className":"dt-right","targets":9},{"orderable":false,"targets":0}],"order":[],"autoWidth":false,"orderClasses":false}},"evals":[],"jsHooks":[]}</script>
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
<div class="datatables html-widget html-fill-item-overflow-hidden html-fill-item" id="htmlwidget-b00bb49b4131cbaeb274" style="width:100%;height:auto;"></div>
<script type="application/json" data-for="htmlwidget-b00bb49b4131cbaeb274">{"x":{"filter":"none","vertical":false,"data":[["1","2","3","4","5","6"],["bjaccard","bjaccard","bjaccard","aichinson","aichinson","aichinson"],["Treatment","Residual","Total","Treatment","Residual","Total"],[6,20,26,6,20,26],[3.92,1.816,5.736,2.42,0.746,3.166],[0.683,0.317,1,0.764,0.236,1],[7.194,null,null,10.81,null,null],[0.001,null,null,0.001,null,null]],"container":"<table class=\"display\">\n  <thead>\n    <tr>\n      <th> <\/th>\n      <th>Distance<\/th>\n      <th>terms<\/th>\n      <th>Df<\/th>\n      <th>SumOfSqs<\/th>\n      <th>R2<\/th>\n      <th>F<\/th>\n      <th>Pr(&gt;F)<\/th>\n    <\/tr>\n  <\/thead>\n<\/table>","options":{"columnDefs":[{"className":"dt-right","targets":[3,4,5,6,7]},{"orderable":false,"targets":0}],"order":[],"autoWidth":false,"orderClasses":false}},"evals":[],"jsHooks":[]}</script>
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
<div class="datatables html-widget html-fill-item-overflow-hidden html-fill-item" id="htmlwidget-85cf7f798e923a02ed60" style="width:100%;height:auto;"></div>
<script type="application/json" data-for="htmlwidget-85cf7f798e923a02ed60">{"x":{"filter":"none","vertical":false,"data":[["1","2","3","4","5","6","7","8"],["bjaccard","bjaccard","bjaccard","bjaccard","aichinson","aichinson","aichinson","aichinson"],["Treatment","Day_of_Treatment","Residual","Total","Treatment","Day_of_Treatment","Residual","Total"],[6,1,19,26,6,1,19,26],[3.918,0.12,1.697,5.736,2.42,0.059,0.687,3.166],[0.683,0.021,0.296,1,0.765,0.019,0.217,1],[7.314,1.34,null,null,11.161,1.642,null,null],[0.001,0.194,null,null,0.001,0.136,null,null]],"container":"<table class=\"display\">\n  <thead>\n    <tr>\n      <th> <\/th>\n      <th>Distance<\/th>\n      <th>terms<\/th>\n      <th>Df<\/th>\n      <th>SumOfSqs<\/th>\n      <th>R2<\/th>\n      <th>F<\/th>\n      <th>Pr(&gt;F)<\/th>\n    <\/tr>\n  <\/thead>\n<\/table>","options":{"columnDefs":[{"className":"dt-right","targets":[3,4,5,6,7]},{"orderable":false,"targets":0}],"order":[],"autoWidth":false,"orderClasses":false}},"evals":[],"jsHooks":[]}</script>
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
<div class="datatables html-widget html-fill-item-overflow-hidden html-fill-item" id="htmlwidget-f7202f4a63f76b2c09a9" style="width:100%;height:auto;"></div>
<script type="application/json" data-for="htmlwidget-f7202f4a63f76b2c09a9">{"x":{"filter":"none","vertical":false,"data":[["1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","23","24","25","26","27","28","29","30","31","32","33","34","35","36","37","38","39","40","41","42"],["bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson"],["UNTREATED","UNTREATED","UNTREATED","UNTREATED","UNTREATED","UNTREATED","CTX","CTX","CTX","CTX","CTX","CTX+HV292.1","CTX+HV292.1","CTX+HV292.1","CTX+HV292.1","HV292.1","HV292.1","HV292.1","VAN","VAN","VAN+CCUG59168","UNTREATED","UNTREATED","UNTREATED","UNTREATED","UNTREATED","UNTREATED","CTX","CTX","CTX","CTX","CTX","CTX+HV292.1","CTX+HV292.1","CTX+HV292.1","CTX+HV292.1","HV292.1","HV292.1","HV292.1","VAN","VAN","VAN+CCUG59168"],["CTX","CTX+HV292.1","HV292.1","VAN","VAN+CCUG59168","CCUG59168","CTX+HV292.1","HV292.1","VAN","VAN+CCUG59168","CCUG59168","HV292.1","VAN","VAN+CCUG59168","CCUG59168","VAN","VAN+CCUG59168","CCUG59168","VAN+CCUG59168","CCUG59168","CCUG59168","CTX","CTX+HV292.1","HV292.1","VAN","VAN+CCUG59168","CCUG59168","CTX+HV292.1","HV292.1","VAN","VAN+CCUG59168","CCUG59168","HV292.1","VAN","VAN+CCUG59168","CCUG59168","VAN","VAN+CCUG59168","CCUG59168","VAN+CCUG59168","CCUG59168","CCUG59168"],[0.623,0.69,0.339,0.666,0.703,0.322,0.192,0.582,0.357,0.434,0.616,0.632,0.412,0.475,0.674,0.617,0.662,0.38,0.245,0.661,0.704,0.744,0.75,0.377,0.766,0.787,0.34,0.207,0.681,0.437,0.497,0.724,0.681,0.495,0.541,0.725,0.7,0.731,0.381,0.297,0.745,0.773],[0.023,0.028,0.033,0.027,0.032,0.027,0.308,0.024,0.019,0.028,0.04,0.036,0.035,0.034,0.025,0.028,0.029,0.039,0.091,0.03,0.028,0.026,0.028,0.021,0.039,0.038,0.028,0.158,0.036,0.033,0.035,0.037,0.032,0.032,0.029,0.029,0.02,0.04,0.024,0.021,0.034,0.037],[0.242,0.069,0.05,0.103,0.052,0.103,0.308,0.168,0.399,0.069,0.044,0.044,0.046,0.048,0.131,0.069,0.055,0.046,0.096,0.052,0.069,0.109,0.09,0.176,0.043,0.044,0.09,0.158,0.05,0.058,0.053,0.047,0.064,0.064,0.072,0.072,0.42,0.042,0.126,0.176,0.055,0.047]],"container":"<table class=\"display\">\n  <thead>\n    <tr>\n      <th> <\/th>\n      <th>Distance<\/th>\n      <th>X1<\/th>\n      <th>X2<\/th>\n      <th>R2<\/th>\n      <th>pval<\/th>\n      <th>pvalFDR<\/th>\n    <\/tr>\n  <\/thead>\n<\/table>","options":{"columnDefs":[{"className":"dt-right","targets":[4,5,6]},{"orderable":false,"targets":0}],"order":[],"autoWidth":false,"orderClasses":false}},"evals":[],"jsHooks":[]}</script>
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
<div class="datatables html-widget html-fill-item-overflow-hidden html-fill-item" id="htmlwidget-175be691d0d33a88348a" style="width:100%;height:auto;"></div>
<script type="application/json" data-for="htmlwidget-175be691d0d33a88348a">{"x":{"filter":"none","vertical":false,"data":[["1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","23","24","25","26","27","28","29","30","31","32","33","34","35","36","37","38","39","40","41","42"],["bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson"],[["UNTREATED"],["UNTREATED"],["UNTREATED"],["UNTREATED"],["UNTREATED"],["UNTREATED"],["CTX"],["CTX"],["CTX"],["CTX"],["CTX"],["CTX+HV292.1"],["CTX+HV292.1"],["CTX+HV292.1"],["CTX+HV292.1"],["HV292.1"],["HV292.1"],["HV292.1"],["VAN"],["VAN"],["VAN+CCUG59168"],["UNTREATED"],["UNTREATED"],["UNTREATED"],["UNTREATED"],["UNTREATED"],["UNTREATED"],["CTX"],["CTX"],["CTX"],["CTX"],["CTX"],["CTX+HV292.1"],["CTX+HV292.1"],["CTX+HV292.1"],["CTX+HV292.1"],["HV292.1"],["HV292.1"],["HV292.1"],["VAN"],["VAN"],["VAN+CCUG59168"]],[["CTX"],["CTX+HV292.1"],["HV292.1"],["VAN"],["VAN+CCUG59168"],["CCUG59168"],["CTX+HV292.1"],["HV292.1"],["VAN"],["VAN+CCUG59168"],["CCUG59168"],["HV292.1"],["VAN"],["VAN+CCUG59168"],["CCUG59168"],["VAN"],["VAN+CCUG59168"],["CCUG59168"],["VAN+CCUG59168"],["CCUG59168"],["CCUG59168"],["CTX"],["CTX+HV292.1"],["HV292.1"],["VAN"],["VAN+CCUG59168"],["CCUG59168"],["CTX+HV292.1"],["HV292.1"],["VAN"],["VAN+CCUG59168"],["CCUG59168"],["HV292.1"],["VAN"],["VAN+CCUG59168"],["CCUG59168"],["VAN"],["VAN+CCUG59168"],["CCUG59168"],["VAN+CCUG59168"],["CCUG59168"],["CCUG59168"]],[[4],[4],[4],[4],[4],[4],[4],[4],[4],[4],[4],[3],[3],[3],[3],[4],[4],[4],[4],[4],[4],[4],[4],[4],[4],[4],[4],[4],[4],[4],[4],[4],[3],[3],[3],[3],[4],[4],[4],[4],[4],[4]],[[4],[3],[4],[4],[4],[4],[3],[4],[4],[4],[4],[4],[4],[4],[4],[4],[4],[4],[4],[4],[4],[4],[3],[4],[4],[4],[4],[3],[4],[4],[4],[4],[4],[4],[4],[4],[4],[4],[4],[4],[4],[4]],[[0.024],[0.032],[0.021],[0.036],[0.031],[0.029],[0.288],[0.031],[0.03],[0.027],[0.028],[0.041],[0.03],[0.035],[0.03],[0.034],[0.024],[0.036],[0.089],[0.031],[0.024],[0.031],[0.026],[0.033],[0.027],[0.026],[0.033],[0.191],[0.027],[0.038],[0.037],[0.041],[0.032],[0.026],[0.023],[0.033],[0.028],[0.023],[0.029],[0.023],[0.028],[0.036]],[[9.89677388009617],[8.5658519248316],[3.08159303806405],[11.9475301471737],[14.2064439078348],[2.85390968967447],[1.1979936975345],[8.37014947434361],[3.32430843071227],[4.59623239962326],[9.64239644092014],[7.22216869065703],[3.48186730457665],[4.36743225712148],[8.08546591476416],[9.68026306357731],[11.7360553115246],[3.67592081270322],[1.94696423517805],[11.7092466544869],[14.2556631739055],[17.4691923866252],[13.6800612386099],[3.6358648970417],[19.6457448578902],[22.1069690825973],[3.09098714684372],[1.29247142592665],[12.8049404372372],[4.65968244284928],[5.91755945244012],[15.7396533234293],[10.5958996111571],[4.75278030245568],[5.52150482554858],[12.4550299373091],[14.01685270907],[16.3235678992188],[3.69862916015435],[2.52998028873845],[17.5309973493175],[20.4593333429457]],[[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999]],[0.042,0.042,0.042,0.042,0.042,0.042,0.288,0.042,0.042,0.042,0.042,0.045,0.042,0.042,0.042,0.042,0.042,0.042,0.093,0.042,0.042,0.042,0.042,0.042,0.042,0.042,0.042,0.191,0.042,0.042,0.042,0.043,0.042,0.042,0.042,0.042,0.042,0.042,0.042,0.042,0.042,0.042]],"container":"<table class=\"display\">\n  <thead>\n    <tr>\n      <th> <\/th>\n      <th>Distance<\/th>\n      <th>Level1<\/th>\n      <th>Level2<\/th>\n      <th>N1<\/th>\n      <th>N2<\/th>\n      <th>p.value<\/th>\n      <th>tw2.stat<\/th>\n      <th>nrep<\/th>\n      <th>pvalFDR<\/th>\n    <\/tr>\n  <\/thead>\n<\/table>","options":{"columnDefs":[{"className":"dt-right","targets":9},{"orderable":false,"targets":0}],"order":[],"autoWidth":false,"orderClasses":false}},"evals":[],"jsHooks":[]}</script>
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
<div class="datatables html-widget html-fill-item-overflow-hidden html-fill-item" id="htmlwidget-6b0771041c3145491948" style="width:100%;height:auto;"></div>
<script type="application/json" data-for="htmlwidget-6b0771041c3145491948">{"x":{"filter":"none","vertical":false,"data":[["1","2","3","4","5","6"],["bjaccard","bjaccard","bjaccard","aichinson","aichinson","aichinson"],["Treatment","Residual","Total","Treatment","Residual","Total"],[6,14,20,6,14,20],[2.803,0.881,3.684,2.014,0.498,2.513],[0.761,0.239,1,0.802,0.198,1],[7.426,null,null,9.437,null,null],[0.001,null,null,0.001,null,null]],"container":"<table class=\"display\">\n  <thead>\n    <tr>\n      <th> <\/th>\n      <th>Distance<\/th>\n      <th>terms<\/th>\n      <th>Df<\/th>\n      <th>SumOfSqs<\/th>\n      <th>R2<\/th>\n      <th>F<\/th>\n      <th>Pr(&gt;F)<\/th>\n    <\/tr>\n  <\/thead>\n<\/table>","options":{"columnDefs":[{"className":"dt-right","targets":[3,4,5,6,7]},{"orderable":false,"targets":0}],"order":[],"autoWidth":false,"orderClasses":false}},"evals":[],"jsHooks":[]}</script>
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
<div class="datatables html-widget html-fill-item-overflow-hidden html-fill-item" id="htmlwidget-0a7605dfcbd2dab75f9a" style="width:100%;height:auto;"></div>
<script type="application/json" data-for="htmlwidget-0a7605dfcbd2dab75f9a">{"x":{"filter":"none","vertical":false,"data":[["1","2","3","4","5","6","7","8"],["bjaccard","bjaccard","bjaccard","bjaccard","aichinson","aichinson","aichinson","aichinson"],["Treatment","Day_of_Treatment","Residual","Total","Treatment","Day_of_Treatment","Residual","Total"],[6,1,13,20,6,1,13,20],[2.803,0.073,0.808,3.684,2.014,0.04,0.458,2.513],[0.761,0.02,0.219,1,0.802,0.016,0.182,1],[7.518,1.173,null,null,9.534,1.144,null,null],[0.001,0.262,null,null,0.001,0.297,null,null]],"container":"<table class=\"display\">\n  <thead>\n    <tr>\n      <th> <\/th>\n      <th>Distance<\/th>\n      <th>terms<\/th>\n      <th>Df<\/th>\n      <th>SumOfSqs<\/th>\n      <th>R2<\/th>\n      <th>F<\/th>\n      <th>Pr(&gt;F)<\/th>\n    <\/tr>\n  <\/thead>\n<\/table>","options":{"columnDefs":[{"className":"dt-right","targets":[3,4,5,6,7]},{"orderable":false,"targets":0}],"order":[],"autoWidth":false,"orderClasses":false}},"evals":[],"jsHooks":[]}</script>
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
```

```
## 'nperm' >= set of all permutations: complete enumeration.
```

```
## Set of permutations < 'minperm'. Generating entire set.
```

```
## 'nperm' >= set of all permutations: complete enumeration.
```

```
## Set of permutations < 'minperm'. Generating entire set.
```

```
## 'nperm' >= set of all permutations: complete enumeration.
```

```
## Set of permutations < 'minperm'. Generating entire set.
```

```
## 'nperm' >= set of all permutations: complete enumeration.
```

```
## Set of permutations < 'minperm'. Generating entire set.
```

```
## 'nperm' >= set of all permutations: complete enumeration.
```

```
## Set of permutations < 'minperm'. Generating entire set.
```

```
## 'nperm' >= set of all permutations: complete enumeration.
```

```
## Set of permutations < 'minperm'. Generating entire set.
```

```
## 'nperm' >= set of all permutations: complete enumeration.
```

```
## Set of permutations < 'minperm'. Generating entire set.
```

```
## 'nperm' >= set of all permutations: complete enumeration.
```

```
## Set of permutations < 'minperm'. Generating entire set.
```

```
## 'nperm' >= set of all permutations: complete enumeration.
```

```
## Set of permutations < 'minperm'. Generating entire set.
```

```
## 'nperm' >= set of all permutations: complete enumeration.
```

```
## Set of permutations < 'minperm'. Generating entire set.
```

```
## 'nperm' >= set of all permutations: complete enumeration.
```

```
## Set of permutations < 'minperm'. Generating entire set.
```

```
## 'nperm' >= set of all permutations: complete enumeration.
```

```
## Set of permutations < 'minperm'. Generating entire set.
```

```
## 'nperm' >= set of all permutations: complete enumeration.
```

```
## Set of permutations < 'minperm'. Generating entire set.
```

```
## 'nperm' >= set of all permutations: complete enumeration.
```

```
## Set of permutations < 'minperm'. Generating entire set.
```

```
## 'nperm' >= set of all permutations: complete enumeration.
```

```
## Set of permutations < 'minperm'. Generating entire set.
```

```
## 'nperm' >= set of all permutations: complete enumeration.
```

```
## Set of permutations < 'minperm'. Generating entire set.
```

```
## 'nperm' >= set of all permutations: complete enumeration.
```

```
## Set of permutations < 'minperm'. Generating entire set.
```

```
## 'nperm' >= set of all permutations: complete enumeration.
```

```
## Set of permutations < 'minperm'. Generating entire set.
```

```
## 'nperm' >= set of all permutations: complete enumeration.
```

```
## Set of permutations < 'minperm'. Generating entire set.
```

```
## 'nperm' >= set of all permutations: complete enumeration.
```

```
## Set of permutations < 'minperm'. Generating entire set.
```

```
## 'nperm' >= set of all permutations: complete enumeration.
```

```
## Set of permutations < 'minperm'. Generating entire set.
```

```
## 'nperm' >= set of all permutations: complete enumeration.
```

```
## Set of permutations < 'minperm'. Generating entire set.
```

```
## 'nperm' >= set of all permutations: complete enumeration.
```

```
## Set of permutations < 'minperm'. Generating entire set.
```

```
## 'nperm' >= set of all permutations: complete enumeration.
```

```
## Set of permutations < 'minperm'. Generating entire set.
```

```
## 'nperm' >= set of all permutations: complete enumeration.
```

```
## Set of permutations < 'minperm'. Generating entire set.
```

```
## 'nperm' >= set of all permutations: complete enumeration.
```

```
## Set of permutations < 'minperm'. Generating entire set.
```

```
## 'nperm' >= set of all permutations: complete enumeration.
```

```
## Set of permutations < 'minperm'. Generating entire set.
```

```
## 'nperm' >= set of all permutations: complete enumeration.
```

```
## Set of permutations < 'minperm'. Generating entire set.
```

```
## 'nperm' >= set of all permutations: complete enumeration.
```

```
## Set of permutations < 'minperm'. Generating entire set.
```

```
## 'nperm' >= set of all permutations: complete enumeration.
```

```
## Set of permutations < 'minperm'. Generating entire set.
```

```
## 'nperm' >= set of all permutations: complete enumeration.
```

```
## Set of permutations < 'minperm'. Generating entire set.
```

```
## 'nperm' >= set of all permutations: complete enumeration.
```

```
## Set of permutations < 'minperm'. Generating entire set.
```

```
## 'nperm' >= set of all permutations: complete enumeration.
```

```
## Set of permutations < 'minperm'. Generating entire set.
```

```
## 'nperm' >= set of all permutations: complete enumeration.
```

```
## Set of permutations < 'minperm'. Generating entire set.
```

```
## 'nperm' >= set of all permutations: complete enumeration.
```

```
## Set of permutations < 'minperm'. Generating entire set.
```

```
## 'nperm' >= set of all permutations: complete enumeration.
```

```
## Set of permutations < 'minperm'. Generating entire set.
```

```
## 'nperm' >= set of all permutations: complete enumeration.
```

```
## Set of permutations < 'minperm'. Generating entire set.
```

```
## 'nperm' >= set of all permutations: complete enumeration.
```

```
## Set of permutations < 'minperm'. Generating entire set.
```

```
## 'nperm' >= set of all permutations: complete enumeration.
```

```
## Set of permutations < 'minperm'. Generating entire set.
```

```
## 'nperm' >= set of all permutations: complete enumeration.
```

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
<div class="datatables html-widget html-fill-item-overflow-hidden html-fill-item" id="htmlwidget-dd76c3579b426da78523" style="width:100%;height:auto;"></div>
<script type="application/json" data-for="htmlwidget-dd76c3579b426da78523">{"x":{"filter":"none","vertical":false,"data":[["1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","23","24","25","26","27","28","29","30","31","32","33","34","35","36","37","38","39","40","41","42"],["bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson"],["UNTREATED","UNTREATED","UNTREATED","UNTREATED","UNTREATED","UNTREATED","CTX","CTX","CTX","CTX","CTX","CTX+HV292.1","CTX+HV292.1","CTX+HV292.1","CTX+HV292.1","HV292.1","HV292.1","HV292.1","VAN","VAN","VAN+CCUG59168","UNTREATED","UNTREATED","UNTREATED","UNTREATED","UNTREATED","UNTREATED","CTX","CTX","CTX","CTX","CTX","CTX+HV292.1","CTX+HV292.1","CTX+HV292.1","CTX+HV292.1","HV292.1","HV292.1","HV292.1","VAN","VAN","VAN+CCUG59168"],["CTX","CTX+HV292.1","HV292.1","VAN","VAN+CCUG59168","CCUG59168","CTX+HV292.1","HV292.1","VAN","VAN+CCUG59168","CCUG59168","HV292.1","VAN","VAN+CCUG59168","CCUG59168","VAN","VAN+CCUG59168","CCUG59168","VAN+CCUG59168","CCUG59168","CCUG59168","CTX","CTX+HV292.1","HV292.1","VAN","VAN+CCUG59168","CCUG59168","CTX+HV292.1","HV292.1","VAN","VAN+CCUG59168","CCUG59168","HV292.1","VAN","VAN+CCUG59168","CCUG59168","VAN","VAN+CCUG59168","CCUG59168","VAN+CCUG59168","CCUG59168","CCUG59168"],[0.68,0.676,0.523,0.698,0.695,0.502,0.613,0.69,0.45,0.413,0.721,0.693,0.674,0.686,0.763,0.693,0.725,0.482,0.329,0.735,0.757,0.741,0.664,0.531,0.726,0.724,0.484,0.672,0.792,0.537,0.514,0.792,0.71,0.701,0.72,0.731,0.772,0.795,0.475,0.395,0.773,0.796],[0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1],[0.191,0.191,0.191,0.191,0.191,0.191,0.191,0.191,0.191,0.191,0.191,0.191,0.191,0.191,0.191,0.191,0.191,0.191,0.191,0.191,0.191,0.191,0.191,0.191,0.191,0.191,0.191,0.191,0.191,0.191,0.191,0.191,0.191,0.191,0.191,0.191,0.191,0.191,0.191,0.191,0.191,0.191]],"container":"<table class=\"display\">\n  <thead>\n    <tr>\n      <th> <\/th>\n      <th>Distance<\/th>\n      <th>X1<\/th>\n      <th>X2<\/th>\n      <th>R2<\/th>\n      <th>pval<\/th>\n      <th>pvalFDR<\/th>\n    <\/tr>\n  <\/thead>\n<\/table>","options":{"columnDefs":[{"className":"dt-right","targets":[4,5,6]},{"orderable":false,"targets":0}],"order":[],"autoWidth":false,"orderClasses":false}},"evals":[],"jsHooks":[]}</script>
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
<div class="datatables html-widget html-fill-item-overflow-hidden html-fill-item" id="htmlwidget-36e9ef7a775d2ff30998" style="width:100%;height:auto;"></div>
<script type="application/json" data-for="htmlwidget-36e9ef7a775d2ff30998">{"x":{"filter":"none","vertical":false,"data":[["1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","23","24","25","26","27","28","29","30","31","32","33","34","35","36","37","38","39","40","41","42"],["bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson"],[["UNTREATED"],["UNTREATED"],["UNTREATED"],["UNTREATED"],["UNTREATED"],["UNTREATED"],["CTX"],["CTX"],["CTX"],["CTX"],["CTX"],["CTX+HV292.1"],["CTX+HV292.1"],["CTX+HV292.1"],["CTX+HV292.1"],["HV292.1"],["HV292.1"],["HV292.1"],["VAN"],["VAN"],["VAN+CCUG59168"],["UNTREATED"],["UNTREATED"],["UNTREATED"],["UNTREATED"],["UNTREATED"],["UNTREATED"],["CTX"],["CTX"],["CTX"],["CTX"],["CTX"],["CTX+HV292.1"],["CTX+HV292.1"],["CTX+HV292.1"],["CTX+HV292.1"],["HV292.1"],["HV292.1"],["HV292.1"],["VAN"],["VAN"],["VAN+CCUG59168"]],[["CTX"],["CTX+HV292.1"],["HV292.1"],["VAN"],["VAN+CCUG59168"],["CCUG59168"],["CTX+HV292.1"],["HV292.1"],["VAN"],["VAN+CCUG59168"],["CCUG59168"],["HV292.1"],["VAN"],["VAN+CCUG59168"],["CCUG59168"],["VAN"],["VAN+CCUG59168"],["CCUG59168"],["VAN+CCUG59168"],["CCUG59168"],["CCUG59168"],["CTX"],["CTX+HV292.1"],["HV292.1"],["VAN"],["VAN+CCUG59168"],["CCUG59168"],["CTX+HV292.1"],["HV292.1"],["VAN"],["VAN+CCUG59168"],["CCUG59168"],["HV292.1"],["VAN"],["VAN+CCUG59168"],["CCUG59168"],["VAN"],["VAN+CCUG59168"],["CCUG59168"],["VAN+CCUG59168"],["CCUG59168"],["CCUG59168"]],[[3],[3],[3],[3],[3],[3],[3],[3],[3],[3],[3],[3],[3],[3],[3],[3],[3],[3],[3],[3],[3],[3],[3],[3],[3],[3],[3],[3],[3],[3],[3],[3],[3],[3],[3],[3],[3],[3],[3],[3],[3],[3]],[[3],[3],[3],[3],[3],[3],[3],[3],[3],[3],[3],[3],[3],[3],[3],[3],[3],[3],[3],[3],[3],[3],[3],[3],[3],[3],[3],[3],[3],[3],[3],[3],[3],[3],[3],[3],[3],[3],[3],[3],[3],[3]],[[0.09],[0.102],[0.118],[0.083],[0.091],[0.105],[0.098],[0.087],[0.098],[0.095],[0.1],[0.088],[0.086],[0.121],[0.108],[0.097],[0.099],[0.093],[0.091],[0.095],[0.095],[0.09],[0.098],[0.117],[0.102],[0.101],[0.098],[0.112],[0.098],[0.088],[0.097],[0.127],[0.092],[0.101],[0.1],[0.096],[0.1],[0.095],[0.107],[0.105],[0.112],[0.094]],[[8.51136481272645],[8.32858669127167],[4.38458967318319],[9.24219011022793],[9.10864095999458],[4.03157373721911],[6.34525007796114],[8.88928491635667],[3.26800903195742],[2.81384306547469],[10.3309608167587],[9.01283323313974],[8.26874761128541],[8.73296061639032],[12.8746231612236],[9.04457343768923],[10.5233345011014],[3.72635693968114],[1.96456105501301],[11.1199345835546],[12.4804925197844],[11.451869296465],[7.90403726500688],[4.52894670820265],[10.622632638945],[10.5147755682571],[3.75608382946069],[8.20470382025644],[15.2092073550007],[4.63515443823741],[4.22380213898984],[15.1875260663548],[9.80565859255954],[9.35574649953762],[10.2754455533249],[10.8450294458474],[13.5067869819737],[15.5552227012796],[3.61436795606982],[2.6132149043023],[13.6193673807053],[15.5735960618277]],[[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999]],[0.119,0.119,0.121,0.119,0.119,0.119,0.119,0.119,0.119,0.119,0.119,0.119,0.119,0.121,0.119,0.119,0.119,0.119,0.119,0.119,0.119,0.123,0.123,0.123,0.123,0.123,0.123,0.123,0.123,0.123,0.123,0.127,0.123,0.123,0.123,0.123,0.123,0.123,0.123,0.123,0.123,0.123]],"container":"<table class=\"display\">\n  <thead>\n    <tr>\n      <th> <\/th>\n      <th>Distance<\/th>\n      <th>Level1<\/th>\n      <th>Level2<\/th>\n      <th>N1<\/th>\n      <th>N2<\/th>\n      <th>p.value<\/th>\n      <th>tw2.stat<\/th>\n      <th>nrep<\/th>\n      <th>pvalFDR<\/th>\n    <\/tr>\n  <\/thead>\n<\/table>","options":{"columnDefs":[{"className":"dt-right","targets":9},{"orderable":false,"targets":0}],"order":[],"autoWidth":false,"orderClasses":false}},"evals":[],"jsHooks":[]}</script>
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

![](NRP72_draft_16S_Figs_beta_PCoA_files/figure-html/unnamed-chunk-114-1.png)<!-- -->


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
<div class="datatables html-widget html-fill-item-overflow-hidden html-fill-item" id="htmlwidget-a5f263681883ceb31574" style="width:100%;height:auto;"></div>
<script type="application/json" data-for="htmlwidget-a5f263681883ceb31574">{"x":{"filter":"none","vertical":false,"data":[["ASV0002","ASV0051","ASV0010","ASV0015","ASV0019","ASV0035","ASV0106","ASV0124"],[-0.597834179220655,0.657937890091773,-0.875873134371687,0.78955848343164,0.506301477829794,0.538606751864561,-0.741001013194891,-0.645690741087984],[-0.332459687317544,0.334520008130732,0.214953317291109,0.128341752631704,0.476402252087666,0.510887546877555,-0.0171919400499977,0.384839333628445],[0.001,0.001,0.001,0.001,0.001,0.001,0.001,0.001],[0.467935149535713,0.544785903058199,0.813358676128535,0.639874204227448,0.483300292246834,0.551103318708659,0.549378064358538,0.565017845834336],["Bacteria","Bacteria","Bacteria","Bacteria","Bacteria","Bacteria","Bacteria","Bacteria"],["Bacteroidota","Firmicutes","Firmicutes","Firmicutes","Proteobacteria","Proteobacteria","Firmicutes","Firmicutes"],["Bacteroidia","Bacilli","Clostridia","Negativicutes","Gammaproteobacteria","Gammaproteobacteria","Clostridia","Clostridia"],["Bacteroidales","Lactobacillales","Lachnospirales","Veillonellales-Selenomonadales","Pseudomonadales","Burkholderiales","Oscillospirales","Monoglobales"],["Bacteroidaceae","Lactobacillaceae","Lachnospiraceae","Veillonellaceae","Pseudomonadaceae","Sutterellaceae","Oscillospiraceae","Monoglobaceae"],[null,null,null,null,null,null,null,null],[null,null,null,null,null,null,null,null],[null,null,null,null,null,null,null,null],[0.0026,0.0026,0.0026,0.0026,0.0026,0.0026,0.0026,0.0026]],"container":"<table class=\"display\">\n  <thead>\n    <tr>\n      <th> <\/th>\n      <th>Axis.1<\/th>\n      <th>Axis.2<\/th>\n      <th>pval<\/th>\n      <th>r<\/th>\n      <th>Kingdom<\/th>\n      <th>Phylum<\/th>\n      <th>Class<\/th>\n      <th>Order<\/th>\n      <th>tax_rank_plot<\/th>\n      <th>Genus<\/th>\n      <th>Species<\/th>\n      <th>Strain<\/th>\n      <th>pval.adj<\/th>\n    <\/tr>\n  <\/thead>\n<\/table>","options":{"columnDefs":[{"className":"dt-right","targets":[1,2,3,4,13]},{"orderable":false,"targets":0}],"order":[],"autoWidth":false,"orderClasses":false}},"evals":[],"jsHooks":[]}</script>
```


```r
out$plot -> PCoA_tmp2_fam

PCoA_tmp2 + 
  theme(legend.position = "none") -> PCoA_tmp2_noleg

PCoA_tmp2_noleg
```

![](NRP72_draft_16S_Figs_beta_PCoA_files/figure-html/unnamed-chunk-116-1.png)<!-- -->

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
<div class="datatables html-widget html-fill-item-overflow-hidden html-fill-item" id="htmlwidget-f017a62ebbbe24ad352e" style="width:100%;height:auto;"></div>
<script type="application/json" data-for="htmlwidget-f017a62ebbbe24ad352e">{"x":{"filter":"none","vertical":false,"data":[["1","2","3","4","5","6"],["Formiat_mM","Acetat_mM","Propionat_mM","Butyrat_mM","Valerat_mM","Total_SCFA_mM"],[0.270162708759593,-0.878418575751238,-0.0488725590367165,-0.794050468978855,-0.117619022423589,-0.78010067992651],[0.551509577899249,-0.234172870112844,-0.302091816938156,-0.454780783257901,-0.678509177560713,-0.433984784767189],["Formiat_mM","Acetat_mM","Propionat_mM","Butyrat_mM","Valerat_mM","Total_SCFA_mM"],[0.001,0.001,0.009,0.001,0.001,0.001],[0.377150703718928,0.82645612732172,0.0936479928877934,0.837341708106209,0.474208938469996,0.796899864231227]],"container":"<table class=\"display\">\n  <thead>\n    <tr>\n      <th> <\/th>\n      <th>label<\/th>\n      <th>Axis.1<\/th>\n      <th>Axis.2<\/th>\n      <th>id<\/th>\n      <th>pval<\/th>\n      <th>r<\/th>\n    <\/tr>\n  <\/thead>\n<\/table>","options":{"columnDefs":[{"className":"dt-right","targets":[2,3,5,6]},{"orderable":false,"targets":0}],"order":[],"autoWidth":false,"orderClasses":false}},"evals":[],"jsHooks":[]}</script>
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

![](NRP72_draft_16S_Figs_beta_PCoA_files/figure-html/unnamed-chunk-119-1.png)<!-- -->


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
<div class="datatables html-widget html-fill-item-overflow-hidden html-fill-item" id="htmlwidget-ce8f12c79e792466f32f" style="width:100%;height:auto;"></div>
<script type="application/json" data-for="htmlwidget-ce8f12c79e792466f32f">{"x":{"filter":"none","vertical":false,"data":[["ASV0002","ASV0012","ASV0051","ASV0010","ASV0015","ASV0035","ASV0106","ASV0124"],[-0.637590724142918,-0.429936467675664,0.700526550112209,-0.776746308038715,0.761430656215353,0.555743853351385,-0.689019065137533,-0.539712190289124],[-0.206234890986704,-0.509585879553196,0.265871052745876,0.404433893018545,-0.110706612899521,0.408775173606298,0.194157556727153,0.486883130283422],[0.001,0.001,0.001,0.001,0.001,0.001,0.001,0.001],[0.449054761773388,0.444523134877431,0.561424864100314,0.766901600873911,0.592032598364227,0.475948373094705,0.512444428957257,0.528344430901268],["Bacteria","Bacteria","Bacteria","Bacteria","Bacteria","Bacteria","Bacteria","Bacteria"],["Bacteroidota","Firmicutes","Firmicutes","Firmicutes","Firmicutes","Proteobacteria","Firmicutes","Firmicutes"],["Bacteroidia","Clostridia","Bacilli","Clostridia","Negativicutes","Gammaproteobacteria","Clostridia","Clostridia"],["Bacteroidales","Oscillospirales","Lactobacillales","Lachnospirales","Veillonellales-Selenomonadales","Burkholderiales","Oscillospirales","Monoglobales"],["Bacteroidaceae","Ruminococcaceae","Lactobacillaceae","Lachnospiraceae","Veillonellaceae","Sutterellaceae","Oscillospiraceae","Monoglobaceae"],[null,null,null,null,null,null,null,null],[null,null,null,null,null,null,null,null],[null,null,null,null,null,null,null,null],[0.00273684210526316,0.00273684210526316,0.00273684210526316,0.00273684210526316,0.00273684210526316,0.00273684210526316,0.00273684210526316,0.00273684210526316]],"container":"<table class=\"display\">\n  <thead>\n    <tr>\n      <th> <\/th>\n      <th>Axis.1<\/th>\n      <th>Axis.2<\/th>\n      <th>pval<\/th>\n      <th>r<\/th>\n      <th>Kingdom<\/th>\n      <th>Phylum<\/th>\n      <th>Class<\/th>\n      <th>Order<\/th>\n      <th>tax_rank_plot<\/th>\n      <th>Genus<\/th>\n      <th>Species<\/th>\n      <th>Strain<\/th>\n      <th>pval.adj<\/th>\n    <\/tr>\n  <\/thead>\n<\/table>","options":{"columnDefs":[{"className":"dt-right","targets":[1,2,3,4,13]},{"orderable":false,"targets":0}],"order":[],"autoWidth":false,"orderClasses":false}},"evals":[],"jsHooks":[]}</script>
```


```r
out$plot -> PCoA_tmp2_fam

PCoA_tmp2 + 
  theme(legend.position = "none") -> PCoA_tmp2_noleg

PCoA_tmp2_noleg
```

![](NRP72_draft_16S_Figs_beta_PCoA_files/figure-html/unnamed-chunk-121-1.png)<!-- -->

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
<div class="datatables html-widget html-fill-item-overflow-hidden html-fill-item" id="htmlwidget-bbbdcf888a2a729dbcc9" style="width:100%;height:auto;"></div>
<script type="application/json" data-for="htmlwidget-bbbdcf888a2a729dbcc9">{"x":{"filter":"none","vertical":false,"data":[["1","2","3","4","5"],["Formiat_mM","Acetat_mM","Butyrat_mM","Valerat_mM","Total_SCFA_mM"],[0.283840135282746,-0.88494812542876,-0.824562216745593,-0.123052332607174,-0.813781416353768],[0.547269594930534,-0.0618837504118045,-0.277260208358263,-0.583396578115785,-0.249053480771384],["Formiat_mM","Acetat_mM","Butyrat_mM","Valerat_mM","Total_SCFA_mM"],[0.001,0.001,0.001,0.001,0.001],[0.380069231932758,0.786962783264907,0.756776072423274,0.355493443917274,0.724267829887086]],"container":"<table class=\"display\">\n  <thead>\n    <tr>\n      <th> <\/th>\n      <th>label<\/th>\n      <th>Axis.1<\/th>\n      <th>Axis.2<\/th>\n      <th>id<\/th>\n      <th>pval<\/th>\n      <th>r<\/th>\n    <\/tr>\n  <\/thead>\n<\/table>","options":{"columnDefs":[{"className":"dt-right","targets":[2,3,5,6]},{"orderable":false,"targets":0}],"order":[],"autoWidth":false,"orderClasses":false}},"evals":[],"jsHooks":[]}</script>
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
<div class="datatables html-widget html-fill-item-overflow-hidden html-fill-item" id="htmlwidget-36bca1387d25533ead4d" style="width:100%;height:auto;"></div>
<script type="application/json" data-for="htmlwidget-36bca1387d25533ead4d">{"x":{"filter":"none","vertical":false,"data":[["1","2","3","4","5","6"],["bjaccard","bjaccard","bjaccard","aichinson","aichinson","aichinson"],["Treatment","Residual","Total","Treatment","Residual","Total"],[6,28,34,6,28,34],[0.703,1.132,1.835,1.023,1.591,2.614],[0.383,0.617,1,0.391,0.609,1],[2.901,null,null,3,null,null],[0.001,null,null,0.001,null,null]],"container":"<table class=\"display\">\n  <thead>\n    <tr>\n      <th> <\/th>\n      <th>Distance<\/th>\n      <th>terms<\/th>\n      <th>Df<\/th>\n      <th>SumOfSqs<\/th>\n      <th>R2<\/th>\n      <th>F<\/th>\n      <th>Pr(&gt;F)<\/th>\n    <\/tr>\n  <\/thead>\n<\/table>","options":{"columnDefs":[{"className":"dt-right","targets":[3,4,5,6,7]},{"orderable":false,"targets":0}],"order":[],"autoWidth":false,"orderClasses":false}},"evals":[],"jsHooks":[]}</script>
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
<div class="datatables html-widget html-fill-item-overflow-hidden html-fill-item" id="htmlwidget-45918177c74a14b8c6f2" style="width:100%;height:auto;"></div>
<script type="application/json" data-for="htmlwidget-45918177c74a14b8c6f2">{"x":{"filter":"none","vertical":false,"data":[["1","2","3","4","5","6","7","8"],["bjaccard","bjaccard","bjaccard","bjaccard","aichinson","aichinson","aichinson","aichinson"],["Treatment","Day_of_Treatment","Residual","Total","Treatment","Day_of_Treatment","Residual","Total"],[6,1,27,34,6,1,27,34],[0.699,0.039,1.092,1.835,1.022,0.052,1.539,2.614],[0.381,0.021,0.595,1,0.391,0.02,0.589,1],[2.88,0.975,null,null,2.987,0.908,null,null],[0.001,0.433,null,null,0.001,0.534,null,null]],"container":"<table class=\"display\">\n  <thead>\n    <tr>\n      <th> <\/th>\n      <th>Distance<\/th>\n      <th>terms<\/th>\n      <th>Df<\/th>\n      <th>SumOfSqs<\/th>\n      <th>R2<\/th>\n      <th>F<\/th>\n      <th>Pr(&gt;F)<\/th>\n    <\/tr>\n  <\/thead>\n<\/table>","options":{"columnDefs":[{"className":"dt-right","targets":[3,4,5,6,7]},{"orderable":false,"targets":0}],"order":[],"autoWidth":false,"orderClasses":false}},"evals":[],"jsHooks":[]}</script>
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

```r
tmp_df %>% 
  DT::datatable()
```

```{=html}
<div class="datatables html-widget html-fill-item-overflow-hidden html-fill-item" id="htmlwidget-56473ac9a4503685a7f3" style="width:100%;height:auto;"></div>
<script type="application/json" data-for="htmlwidget-56473ac9a4503685a7f3">{"x":{"filter":"none","vertical":false,"data":[["1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","23","24","25","26","27","28","29","30","31","32","33","34","35","36","37","38","39","40","41","42"],["bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson"],["UNTREATED","UNTREATED","UNTREATED","UNTREATED","UNTREATED","UNTREATED","CTX","CTX","CTX","CTX","CTX","CTX+HV292.1","CTX+HV292.1","CTX+HV292.1","CTX+HV292.1","HV292.1","HV292.1","HV292.1","VAN","VAN","VAN+CCUG59168","UNTREATED","UNTREATED","UNTREATED","UNTREATED","UNTREATED","UNTREATED","CTX","CTX","CTX","CTX","CTX","CTX+HV292.1","CTX+HV292.1","CTX+HV292.1","CTX+HV292.1","HV292.1","HV292.1","HV292.1","VAN","VAN","VAN+CCUG59168"],["CTX","CTX+HV292.1","HV292.1","VAN","VAN+CCUG59168","CCUG59168","CTX+HV292.1","HV292.1","VAN","VAN+CCUG59168","CCUG59168","HV292.1","VAN","VAN+CCUG59168","CCUG59168","VAN","VAN+CCUG59168","CCUG59168","VAN+CCUG59168","CCUG59168","CCUG59168","CTX","CTX+HV292.1","HV292.1","VAN","VAN+CCUG59168","CCUG59168","CTX+HV292.1","HV292.1","VAN","VAN+CCUG59168","CCUG59168","HV292.1","VAN","VAN+CCUG59168","CCUG59168","VAN","VAN+CCUG59168","CCUG59168","VAN+CCUG59168","CCUG59168","CCUG59168"],[0.208,0.155,0.234,0.149,0.22,0.193,0.139,0.277,0.319,0.363,0.395,0.29,0.268,0.31,0.321,0.391,0.323,0.597,0.13,0.142,0.153,0.201,0.179,0.235,0.162,0.227,0.228,0.126,0.284,0.295,0.357,0.39,0.274,0.264,0.337,0.353,0.384,0.346,0.571,0.14,0.165,0.156],[0.031,0.141,0.042,0.177,0.015,0.096,0.118,0.033,0.009,0.002,0.015,0.01,0.01,0.001,0.016,0.031,0.005,0.1,0.217,0.717,0.146,0.028,0.034,0.048,0.129,0.003,0.047,0.129,0.048,0.002,0.002,0.02,0.012,0.007,0.001,0.015,0.028,0.013,0.1,0.162,0.484,0.128],[0.062,0.174,0.068,0.196,0.042,0.144,0.155,0.058,0.047,0.021,0.042,0.038,0.038,0.021,0.037,0.062,0.035,0.14,0.228,0.717,0.17,0.056,0.06,0.07,0.146,0.016,0.076,0.146,0.07,0.017,0.017,0.047,0.042,0.029,0.021,0.039,0.056,0.039,0.131,0.17,0.484,0.158]],"container":"<table class=\"display\">\n  <thead>\n    <tr>\n      <th> <\/th>\n      <th>Distance<\/th>\n      <th>X1<\/th>\n      <th>X2<\/th>\n      <th>R2<\/th>\n      <th>pval<\/th>\n      <th>pvalFDR<\/th>\n    <\/tr>\n  <\/thead>\n<\/table>","options":{"columnDefs":[{"className":"dt-right","targets":[4,5,6]},{"orderable":false,"targets":0}],"order":[],"autoWidth":false,"orderClasses":false}},"evals":[],"jsHooks":[]}</script>
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
<div class="datatables html-widget html-fill-item-overflow-hidden html-fill-item" id="htmlwidget-057f6e9b1f30e2dbee53" style="width:100%;height:auto;"></div>
<script type="application/json" data-for="htmlwidget-057f6e9b1f30e2dbee53">{"x":{"filter":"none","vertical":false,"data":[["1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","23","24","25","26","27","28","29","30","31","32","33","34","35","36","37","38","39","40","41","42"],["bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson"],[["UNTREATED"],["UNTREATED"],["UNTREATED"],["UNTREATED"],["UNTREATED"],["UNTREATED"],["CTX"],["CTX"],["CTX"],["CTX"],["CTX"],["CTX+HV292.1"],["CTX+HV292.1"],["CTX+HV292.1"],["CTX+HV292.1"],["HV292.1"],["HV292.1"],["HV292.1"],["VAN"],["VAN"],["VAN+CCUG59168"],["UNTREATED"],["UNTREATED"],["UNTREATED"],["UNTREATED"],["UNTREATED"],["UNTREATED"],["CTX"],["CTX"],["CTX"],["CTX"],["CTX"],["CTX+HV292.1"],["CTX+HV292.1"],["CTX+HV292.1"],["CTX+HV292.1"],["HV292.1"],["HV292.1"],["HV292.1"],["VAN"],["VAN"],["VAN+CCUG59168"]],[["CTX"],["CTX+HV292.1"],["HV292.1"],["VAN"],["VAN+CCUG59168"],["CCUG59168"],["CTX+HV292.1"],["HV292.1"],["VAN"],["VAN+CCUG59168"],["CCUG59168"],["HV292.1"],["VAN"],["VAN+CCUG59168"],["CCUG59168"],["VAN"],["VAN+CCUG59168"],["CCUG59168"],["VAN+CCUG59168"],["CCUG59168"],["CCUG59168"],["CTX"],["CTX+HV292.1"],["HV292.1"],["VAN"],["VAN+CCUG59168"],["CCUG59168"],["CTX+HV292.1"],["HV292.1"],["VAN"],["VAN+CCUG59168"],["CCUG59168"],["HV292.1"],["VAN"],["VAN+CCUG59168"],["CCUG59168"],["VAN"],["VAN+CCUG59168"],["CCUG59168"],["VAN+CCUG59168"],["CCUG59168"],["CCUG59168"]],[[6],[6],[6],[6],[6],[6],[5],[5],[5],[5],[5],[6],[6],[6],[6],[3],[3],[3],[4],[4],[8],[6],[6],[6],[6],[6],[6],[5],[5],[5],[5],[5],[6],[6],[6],[6],[3],[3],[3],[4],[4],[8]],[[5],[6],[3],[4],[8],[3],[6],[3],[4],[8],[3],[3],[4],[8],[3],[4],[8],[3],[8],[3],[3],[5],[6],[3],[4],[8],[3],[6],[3],[4],[8],[3],[3],[4],[8],[3],[4],[8],[3],[8],[3],[3]],[[0.036],[0.117],[0.018],[0.159],[0.01],[0.075],[0.11],[0.018],[0.011],[0.002],[0.021],[0.012],[0.015],[0.001],[0.011],[0.031],[0.009],[0.117],[0.194],[0.569],[0.162],[0.018],[0.037],[0.027],[0.13],[0.008],[0.033],[0.163],[0.021],[0.011],[0.002],[0.025],[0.009],[0.003],[0.001],[0.015],[0.034],[0.009],[0.103],[0.186],[0.383],[0.163]],[[2.45029703513863],[1.83149534186929],[3.59631244865084],[1.46951473923253],[3.28864431374897],[2.17022334607896],[1.41784644473462],[3.23293227343592],[3.23231637257388],[6.53454679883874],[4.37171369779211],[4.15591028387124],[2.69415019089351],[5.7449399977476],[3.47095672660226],[3.96274072758841],[8.22666597599194],[5.92690778071803],[1.52253123836056],[0.897603749688104],[2.12541595550425],[2.33199022065062],[2.17872808869365],[3.22335217243648],[1.56789168251106],[3.42764017363229],[2.57674073367949],[1.26579096941619],[3.0601681404366],[2.84725844981694],[6.2642805530213],[4.21680064281827],[3.38342177726307],[2.58314622108128],[6.42673611397581],[3.92873097863613],[3.70518159041559],[7.67577866825485],[5.31839422253925],[1.56099464336391],[1.07848120760732],[2.05363506633381]],[[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999]],[0.058,0.145,0.038,0.179,0.036,0.112,0.145,0.038,0.036,0.021,0.04,0.036,0.038,0.021,0.036,0.054,0.036,0.145,0.204,0.569,0.179,0.042,0.052,0.047,0.161,0.032,0.051,0.18,0.044,0.033,0.021,0.047,0.032,0.021,0.021,0.039,0.051,0.032,0.135,0.195,0.383,0.18]],"container":"<table class=\"display\">\n  <thead>\n    <tr>\n      <th> <\/th>\n      <th>Distance<\/th>\n      <th>Level1<\/th>\n      <th>Level2<\/th>\n      <th>N1<\/th>\n      <th>N2<\/th>\n      <th>p.value<\/th>\n      <th>tw2.stat<\/th>\n      <th>nrep<\/th>\n      <th>pvalFDR<\/th>\n    <\/tr>\n  <\/thead>\n<\/table>","options":{"columnDefs":[{"className":"dt-right","targets":9},{"orderable":false,"targets":0}],"order":[],"autoWidth":false,"orderClasses":false}},"evals":[],"jsHooks":[]}</script>
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
<div class="datatables html-widget html-fill-item-overflow-hidden html-fill-item" id="htmlwidget-c0c4ca2e9381bd6fbebd" style="width:100%;height:auto;"></div>
<script type="application/json" data-for="htmlwidget-c0c4ca2e9381bd6fbebd">{"x":{"filter":"none","vertical":false,"data":[["1","2","3","4","5","6"],["bjaccard","bjaccard","bjaccard","aichinson","aichinson","aichinson"],["Treatment","Residual","Total","Treatment","Residual","Total"],[6,45,51,6,45,51],[5.821,4.874,10.694,4.265,3.149,7.414],[0.544,0.456,1,0.575,0.425,1],[8.958,null,null,10.158,null,null],[0.001,null,null,0.001,null,null]],"container":"<table class=\"display\">\n  <thead>\n    <tr>\n      <th> <\/th>\n      <th>Distance<\/th>\n      <th>terms<\/th>\n      <th>Df<\/th>\n      <th>SumOfSqs<\/th>\n      <th>R2<\/th>\n      <th>F<\/th>\n      <th>Pr(&gt;F)<\/th>\n    <\/tr>\n  <\/thead>\n<\/table>","options":{"columnDefs":[{"className":"dt-right","targets":[3,4,5,6,7]},{"orderable":false,"targets":0}],"order":[],"autoWidth":false,"orderClasses":false}},"evals":[],"jsHooks":[]}</script>
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
<div class="datatables html-widget html-fill-item-overflow-hidden html-fill-item" id="htmlwidget-767209723a30863ca56b" style="width:100%;height:auto;"></div>
<script type="application/json" data-for="htmlwidget-767209723a30863ca56b">{"x":{"filter":"none","vertical":false,"data":[["1","2","3","4","5","6","7","8"],["bjaccard","bjaccard","bjaccard","bjaccard","aichinson","aichinson","aichinson","aichinson"],["Treatment","Day_of_Treatment","Residual","Total","Treatment","Day_of_Treatment","Residual","Total"],[6,1,44,51,6,1,44,51],[5.796,0.291,4.582,10.694,4.248,0.171,2.978,7.414],[0.542,0.027,0.428,1,0.573,0.023,0.402,1],[9.276,2.796,null,null,10.461,2.525,null,null],[0.001,0.02,null,null,0.001,0.028,null,null]],"container":"<table class=\"display\">\n  <thead>\n    <tr>\n      <th> <\/th>\n      <th>Distance<\/th>\n      <th>terms<\/th>\n      <th>Df<\/th>\n      <th>SumOfSqs<\/th>\n      <th>R2<\/th>\n      <th>F<\/th>\n      <th>Pr(&gt;F)<\/th>\n    <\/tr>\n  <\/thead>\n<\/table>","options":{"columnDefs":[{"className":"dt-right","targets":[3,4,5,6,7]},{"orderable":false,"targets":0}],"order":[],"autoWidth":false,"orderClasses":false}},"evals":[],"jsHooks":[]}</script>
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
<div class="datatables html-widget html-fill-item-overflow-hidden html-fill-item" id="htmlwidget-23741672a15783b59cb9" style="width:100%;height:auto;"></div>
<script type="application/json" data-for="htmlwidget-23741672a15783b59cb9">{"x":{"filter":"none","vertical":false,"data":[["1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","23","24","25","26","27","28","29","30","31","32","33","34","35","36","37","38","39","40","41","42"],["bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson"],["UNTREATED","UNTREATED","UNTREATED","UNTREATED","UNTREATED","UNTREATED","CTX","CTX","CTX","CTX","CTX","CTX+HV292.1","CTX+HV292.1","CTX+HV292.1","CTX+HV292.1","HV292.1","HV292.1","HV292.1","VAN","VAN","VAN+CCUG59168","UNTREATED","UNTREATED","UNTREATED","UNTREATED","UNTREATED","UNTREATED","CTX","CTX","CTX","CTX","CTX","CTX+HV292.1","CTX+HV292.1","CTX+HV292.1","CTX+HV292.1","HV292.1","HV292.1","HV292.1","VAN","VAN","VAN+CCUG59168"],["CTX","CTX+HV292.1","HV292.1","VAN","VAN+CCUG59168","CCUG59168","CTX+HV292.1","HV292.1","VAN","VAN+CCUG59168","CCUG59168","HV292.1","VAN","VAN+CCUG59168","CCUG59168","VAN","VAN+CCUG59168","CCUG59168","VAN+CCUG59168","CCUG59168","CCUG59168","CTX","CTX+HV292.1","HV292.1","VAN","VAN+CCUG59168","CCUG59168","CTX+HV292.1","HV292.1","VAN","VAN+CCUG59168","CCUG59168","HV292.1","VAN","VAN+CCUG59168","CCUG59168","VAN","VAN+CCUG59168","CCUG59168","VAN+CCUG59168","CCUG59168","CCUG59168"],[0.303,0.376,0.145,0.426,0.605,0.163,0.068,0.138,0.354,0.503,0.318,0.238,0.398,0.568,0.467,0.142,0.295,0.533,0.133,0.384,0.593,0.291,0.33,0.149,0.505,0.61,0.188,0.067,0.164,0.409,0.503,0.361,0.234,0.463,0.569,0.463,0.253,0.427,0.505,0.129,0.526,0.68],[0.001,0.001,0.201,0.001,0.001,0.02,0.259,0.27,0.001,0.001,0.002,0.125,0.001,0.001,0.001,0.18,0.108,0.167,0.012,0.002,0.001,0.001,0.001,0.184,0.001,0.001,0.005,0.257,0.169,0.001,0.001,0.001,0.115,0.001,0.001,0.003,0.09,0.125,0.167,0.01,0.001,0.002],[0.004,0.004,0.222,0.004,0.004,0.03,0.272,0.27,0.004,0.004,0.004,0.164,0.004,0.004,0.004,0.21,0.151,0.206,0.019,0.004,0.004,0.004,0.004,0.193,0.004,0.004,0.008,0.257,0.187,0.004,0.004,0.004,0.151,0.004,0.004,0.005,0.126,0.154,0.194,0.015,0.004,0.004]],"container":"<table class=\"display\">\n  <thead>\n    <tr>\n      <th> <\/th>\n      <th>Distance<\/th>\n      <th>X1<\/th>\n      <th>X2<\/th>\n      <th>R2<\/th>\n      <th>pval<\/th>\n      <th>pvalFDR<\/th>\n    <\/tr>\n  <\/thead>\n<\/table>","options":{"columnDefs":[{"className":"dt-right","targets":[4,5,6]},{"orderable":false,"targets":0}],"order":[],"autoWidth":false,"orderClasses":false}},"evals":[],"jsHooks":[]}</script>
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
<div class="datatables html-widget html-fill-item-overflow-hidden html-fill-item" id="htmlwidget-d6facaf5e0f5b7b2b1fa" style="width:100%;height:auto;"></div>
<script type="application/json" data-for="htmlwidget-d6facaf5e0f5b7b2b1fa">{"x":{"filter":"none","vertical":false,"data":[["1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","23","24","25","26","27","28","29","30","31","32","33","34","35","36","37","38","39","40","41","42"],["bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson"],[["UNTREATED"],["UNTREATED"],["UNTREATED"],["UNTREATED"],["UNTREATED"],["UNTREATED"],["CTX"],["CTX"],["CTX"],["CTX"],["CTX"],["CTX+HV292.1"],["CTX+HV292.1"],["CTX+HV292.1"],["CTX+HV292.1"],["HV292.1"],["HV292.1"],["HV292.1"],["VAN"],["VAN"],["VAN+CCUG59168"],["UNTREATED"],["UNTREATED"],["UNTREATED"],["UNTREATED"],["UNTREATED"],["UNTREATED"],["CTX"],["CTX"],["CTX"],["CTX"],["CTX"],["CTX+HV292.1"],["CTX+HV292.1"],["CTX+HV292.1"],["CTX+HV292.1"],["HV292.1"],["HV292.1"],["HV292.1"],["VAN"],["VAN"],["VAN+CCUG59168"]],[["CTX"],["CTX+HV292.1"],["HV292.1"],["VAN"],["VAN+CCUG59168"],["CCUG59168"],["CTX+HV292.1"],["HV292.1"],["VAN"],["VAN+CCUG59168"],["CCUG59168"],["HV292.1"],["VAN"],["VAN+CCUG59168"],["CCUG59168"],["VAN"],["VAN+CCUG59168"],["CCUG59168"],["VAN+CCUG59168"],["CCUG59168"],["CCUG59168"],["CTX"],["CTX+HV292.1"],["HV292.1"],["VAN"],["VAN+CCUG59168"],["CCUG59168"],["CTX+HV292.1"],["HV292.1"],["VAN"],["VAN+CCUG59168"],["CCUG59168"],["HV292.1"],["VAN"],["VAN+CCUG59168"],["CCUG59168"],["VAN"],["VAN+CCUG59168"],["CCUG59168"],["VAN+CCUG59168"],["CCUG59168"],["CCUG59168"]],[[9],[9],[9],[9],[9],[9],[10],[10],[10],[10],[10],[8],[8],[8],[8],[1],[1],[1],[11],[11],[8],[9],[9],[9],[9],[9],[9],[10],[10],[10],[10],[10],[8],[8],[8],[8],[1],[1],[1],[11],[11],[8]],[[10],[8],[1],[11],[8],[5],[8],[1],[11],[8],[5],[1],[11],[8],[5],[11],[8],[5],[8],[5],[5],[10],[8],[1],[11],[8],[5],[8],[1],[11],[8],[5],[1],[11],[8],[5],[11],[8],[5],[8],[5],[5]],[[0.001],[0.001],[null],[0.001],[0.001],[0.006],[0.215],[null],[0.001],[0.001],[0.002],[null],[0.001],[0.001],[0.001],[null],[null],[null],[0.01],[0.002],[0.001],[0.001],[0.001],[null],[0.001],[0.001],[0.002],[0.236],[null],[0.001],[0.001],[0.001],[null],[0.001],[0.001],[0.002],[null],[null],[null],[0.006],[0.001],[0.001]],[[7.72784219461444],[8.89594337322294],[null],[15.2762185235015],[21.4739748973799],[2.95644133831612],[1.24902845777323],[null],[10.7857585422615],[15.4665320438668],[9.7823563636804],[null],[13.481851336362],[18.4304810527501],[12.3605761519701],[null],[null],[null],[2.74386735082863],[17.0489797847579],[23.5115717122834],[7.07975481018723],[7.37620103716932],[null],[18.8501498305443],[23.5797812111834],[3.35471165582704],[1.17510680306781],[null],[13.1167405141289],[16.7995781436398],[9.91176048208218],[null],[15.2191758306026],[18.4893641997935],[11.1025593780595],[null],[null],[null],[2.65533238166408],[21.6350858929726],[26.8527111401793]],[[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999]],[0.002,0.002,null,0.002,0.002,0.007,0.215,null,0.002,0.002,0.002,null,0.002,0.002,0.002,null,null,null,0.011,0.002,0.002,0.001,0.001,null,0.001,0.001,0.002,0.236,null,0.001,0.001,0.001,null,0.001,0.001,0.002,null,null,null,0.006,0.001,0.001]],"container":"<table class=\"display\">\n  <thead>\n    <tr>\n      <th> <\/th>\n      <th>Distance<\/th>\n      <th>Level1<\/th>\n      <th>Level2<\/th>\n      <th>N1<\/th>\n      <th>N2<\/th>\n      <th>p.value<\/th>\n      <th>tw2.stat<\/th>\n      <th>nrep<\/th>\n      <th>pvalFDR<\/th>\n    <\/tr>\n  <\/thead>\n<\/table>","options":{"columnDefs":[{"className":"dt-right","targets":9},{"orderable":false,"targets":0}],"order":[],"autoWidth":false,"orderClasses":false}},"evals":[],"jsHooks":[]}</script>
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
<div class="datatables html-widget html-fill-item-overflow-hidden html-fill-item" id="htmlwidget-2f8220dda1413fa7b913" style="width:100%;height:auto;"></div>
<script type="application/json" data-for="htmlwidget-2f8220dda1413fa7b913">{"x":{"filter":"none","vertical":false,"data":[["1","2","3","4","5","6"],["bjaccard","bjaccard","bjaccard","aichinson","aichinson","aichinson"],["Treatment","Residual","Total","Treatment","Residual","Total"],[6,27,33,6,27,33],[4.011,2.735,6.746,3.022,1.82,4.842],[0.595,0.405,1,0.624,0.376,1],[6.598,null,null,7.474,null,null],[0.001,null,null,0.001,null,null]],"container":"<table class=\"display\">\n  <thead>\n    <tr>\n      <th> <\/th>\n      <th>Distance<\/th>\n      <th>terms<\/th>\n      <th>Df<\/th>\n      <th>SumOfSqs<\/th>\n      <th>R2<\/th>\n      <th>F<\/th>\n      <th>Pr(&gt;F)<\/th>\n    <\/tr>\n  <\/thead>\n<\/table>","options":{"columnDefs":[{"className":"dt-right","targets":[3,4,5,6,7]},{"orderable":false,"targets":0}],"order":[],"autoWidth":false,"orderClasses":false}},"evals":[],"jsHooks":[]}</script>
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
<div class="datatables html-widget html-fill-item-overflow-hidden html-fill-item" id="htmlwidget-4d98ecc6a7bce25e95cf" style="width:100%;height:auto;"></div>
<script type="application/json" data-for="htmlwidget-4d98ecc6a7bce25e95cf">{"x":{"filter":"none","vertical":false,"data":[["1","2","3","4","5","6","7","8"],["bjaccard","bjaccard","bjaccard","bjaccard","aichinson","aichinson","aichinson","aichinson"],["Treatment","Day_of_Treatment","Residual","Total","Treatment","Day_of_Treatment","Residual","Total"],[6,1,26,33,6,1,26,33],[4.004,0.167,2.569,6.746,3.011,0.107,1.713,4.842],[0.594,0.025,0.381,1,0.622,0.022,0.354,1],[6.756,1.687,null,null,7.62,1.626,null,null],[0.001,0.124,null,null,0.001,0.143,null,null]],"container":"<table class=\"display\">\n  <thead>\n    <tr>\n      <th> <\/th>\n      <th>Distance<\/th>\n      <th>terms<\/th>\n      <th>Df<\/th>\n      <th>SumOfSqs<\/th>\n      <th>R2<\/th>\n      <th>F<\/th>\n      <th>Pr(&gt;F)<\/th>\n    <\/tr>\n  <\/thead>\n<\/table>","options":{"columnDefs":[{"className":"dt-right","targets":[3,4,5,6,7]},{"orderable":false,"targets":0}],"order":[],"autoWidth":false,"orderClasses":false}},"evals":[],"jsHooks":[]}</script>
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
```

```
## 'nperm' >= set of all permutations: complete enumeration.
```

```
## Set of permutations < 'minperm'. Generating entire set.
```

```
## 'nperm' >= set of all permutations: complete enumeration.
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
<div class="datatables html-widget html-fill-item-overflow-hidden html-fill-item" id="htmlwidget-7a55dcb7a06bbff1b49e" style="width:100%;height:auto;"></div>
<script type="application/json" data-for="htmlwidget-7a55dcb7a06bbff1b49e">{"x":{"filter":"none","vertical":false,"data":[["1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","23","24","25","26","27","28","29","30","31","32","33","34","35","36","37","38","39","40","41","42"],["bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson"],["UNTREATED","UNTREATED","UNTREATED","UNTREATED","UNTREATED","UNTREATED","CTX","CTX","CTX","CTX","CTX","CTX+HV292.1","CTX+HV292.1","CTX+HV292.1","CTX+HV292.1","HV292.1","HV292.1","HV292.1","VAN","VAN","VAN+CCUG59168","UNTREATED","UNTREATED","UNTREATED","UNTREATED","UNTREATED","UNTREATED","CTX","CTX","CTX","CTX","CTX","CTX+HV292.1","CTX+HV292.1","CTX+HV292.1","CTX+HV292.1","HV292.1","HV292.1","HV292.1","VAN","VAN","VAN+CCUG59168"],["CTX","CTX+HV292.1","HV292.1","VAN","VAN+CCUG59168","CCUG59168","CTX+HV292.1","HV292.1","VAN","VAN+CCUG59168","CCUG59168","HV292.1","VAN","VAN+CCUG59168","CCUG59168","VAN","VAN+CCUG59168","CCUG59168","VAN+CCUG59168","CCUG59168","CCUG59168","CTX","CTX+HV292.1","HV292.1","VAN","VAN+CCUG59168","CCUG59168","CTX+HV292.1","HV292.1","VAN","VAN+CCUG59168","CCUG59168","HV292.1","VAN","VAN+CCUG59168","CCUG59168","VAN","VAN+CCUG59168","CCUG59168","VAN+CCUG59168","CCUG59168","CCUG59168"],[0.369,0.393,0.175,0.456,0.628,0.173,0.13,0.325,0.423,0.607,0.458,0.347,0.421,0.598,0.52,0.212,0.408,0.763,0.148,0.397,0.622,0.334,0.342,0.176,0.511,0.609,0.208,0.148,0.332,0.499,0.619,0.472,0.346,0.49,0.608,0.524,0.345,0.555,0.698,0.14,0.541,0.718],[0.001,0.001,0.248,0.001,0.002,0.136,0.205,0.277,0.001,0.001,0.008,0.167,0.001,0.009,0.022,0.222,0.167,0.25,0.046,0.015,0.018,0.002,0.002,0.263,0.001,0.001,0.034,0.117,0.278,0.002,0.005,0.011,0.167,0.004,0.004,0.019,0.124,0.167,0.25,0.065,0.012,0.017],[0.006,0.006,0.274,0.006,0.006,0.204,0.253,0.277,0.006,0.006,0.021,0.226,0.006,0.021,0.038,0.259,0.226,0.262,0.074,0.032,0.034,0.011,0.011,0.276,0.014,0.014,0.055,0.164,0.278,0.011,0.013,0.026,0.2,0.013,0.013,0.033,0.163,0.2,0.276,0.098,0.025,0.032]],"container":"<table class=\"display\">\n  <thead>\n    <tr>\n      <th> <\/th>\n      <th>Distance<\/th>\n      <th>X1<\/th>\n      <th>X2<\/th>\n      <th>R2<\/th>\n      <th>pval<\/th>\n      <th>pvalFDR<\/th>\n    <\/tr>\n  <\/thead>\n<\/table>","options":{"columnDefs":[{"className":"dt-right","targets":[4,5,6]},{"orderable":false,"targets":0}],"order":[],"autoWidth":false,"orderClasses":false}},"evals":[],"jsHooks":[]}</script>
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
<div class="datatables html-widget html-fill-item-overflow-hidden html-fill-item" id="htmlwidget-2dc7173b9d96e4df10d9" style="width:100%;height:auto;"></div>
<script type="application/json" data-for="htmlwidget-2dc7173b9d96e4df10d9">{"x":{"filter":"none","vertical":false,"data":[["1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","23","24","25","26","27","28","29","30","31","32","33","34","35","36","37","38","39","40","41","42"],["bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson","aichinson"],[["UNTREATED"],["UNTREATED"],["UNTREATED"],["UNTREATED"],["UNTREATED"],["UNTREATED"],["CTX"],["CTX"],["CTX"],["CTX"],["CTX"],["CTX+HV292.1"],["CTX+HV292.1"],["CTX+HV292.1"],["CTX+HV292.1"],["HV292.1"],["HV292.1"],["HV292.1"],["VAN"],["VAN"],["VAN+CCUG59168"],["UNTREATED"],["UNTREATED"],["UNTREATED"],["UNTREATED"],["UNTREATED"],["UNTREATED"],["CTX"],["CTX"],["CTX"],["CTX"],["CTX"],["CTX+HV292.1"],["CTX+HV292.1"],["CTX+HV292.1"],["CTX+HV292.1"],["HV292.1"],["HV292.1"],["HV292.1"],["VAN"],["VAN"],["VAN+CCUG59168"]],[["CTX"],["CTX+HV292.1"],["HV292.1"],["VAN"],["VAN+CCUG59168"],["CCUG59168"],["CTX+HV292.1"],["HV292.1"],["VAN"],["VAN+CCUG59168"],["CCUG59168"],["HV292.1"],["VAN"],["VAN+CCUG59168"],["CCUG59168"],["VAN"],["VAN+CCUG59168"],["CCUG59168"],["VAN+CCUG59168"],["CCUG59168"],["CCUG59168"],["CTX"],["CTX+HV292.1"],["HV292.1"],["VAN"],["VAN+CCUG59168"],["CCUG59168"],["CTX+HV292.1"],["HV292.1"],["VAN"],["VAN+CCUG59168"],["CCUG59168"],["HV292.1"],["VAN"],["VAN+CCUG59168"],["CCUG59168"],["VAN"],["VAN+CCUG59168"],["CCUG59168"],["VAN+CCUG59168"],["CCUG59168"],["CCUG59168"]],[[7],[7],[7],[7],[7],[7],[6],[6],[6],[6],[6],[5],[5],[5],[5],[1],[1],[1],[7],[7],[5],[7],[7],[7],[7],[7],[7],[6],[6],[6],[6],[6],[5],[5],[5],[5],[1],[1],[1],[7],[7],[5]],[[6],[5],[1],[7],[5],[3],[5],[1],[7],[5],[3],[1],[7],[5],[3],[7],[5],[3],[5],[3],[3],[6],[5],[1],[7],[5],[3],[5],[1],[7],[5],[3],[1],[7],[5],[3],[7],[5],[3],[5],[3],[3]],[[0.001],[0.002],[null],[0.001],[0.004],[0.051],[0.221],[null],[0.002],[0.003],[0.008],[null],[0.003],[0.009],[0.016],[null],[null],[null],[0.039],[0.01],[0.015],[0.002],[0.004],[null],[0.001],[0.002],[0.02],[0.109],[null],[0.001],[0.005],[0.013],[null],[0.004],[0.012],[0.021],[null],[null],[null],[0.038],[0.009],[0.019]],[[6.35052577649374],[6.1818491617836],[null],[10.0470882465917],[13.8568814580884],[2.80887730939864],[1.32354820000393],[null],[8.91549021848736],[12.6264674774694],[9.4478334682308],[null],[8.91203849916163],[11.8863757175525],[9.47805431924738],[null],[null],[null],[1.84581617431872],[11.7983849404617],[16.0582527959237],[5.57961421878733],[5.22900655236687],[null],[12.552970577104],[15.9457877717697],[3.03704813886125],[1.54767557636096],[null],[11.31239264944],[14.6062118591685],[8.14143945481884],[null],[10.0967471671305],[12.3899850632713],[8.32598874214175],[null],[null],[null],[1.74938435550998],[14.8802055908688],[18.8051845046142]],[[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999]],[0.007,0.007,null,0.007,0.009,0.055,0.221,null,0.007,0.007,0.015,null,0.007,0.015,0.02,null,null,null,0.045,0.015,0.02,0.007,0.01,null,0.007,0.007,0.024,0.109,null,0.007,0.011,0.02,null,0.01,0.02,0.024,null,null,null,0.041,0.017,0.024]],"container":"<table class=\"display\">\n  <thead>\n    <tr>\n      <th> <\/th>\n      <th>Distance<\/th>\n      <th>Level1<\/th>\n      <th>Level2<\/th>\n      <th>N1<\/th>\n      <th>N2<\/th>\n      <th>p.value<\/th>\n      <th>tw2.stat<\/th>\n      <th>nrep<\/th>\n      <th>pvalFDR<\/th>\n    <\/tr>\n  <\/thead>\n<\/table>","options":{"columnDefs":[{"className":"dt-right","targets":9},{"orderable":false,"targets":0}],"order":[],"autoWidth":false,"orderClasses":false}},"evals":[],"jsHooks":[]}</script>
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

