---
title: "Merge chicken1 + human1 gene catalogue + quantification + bins"
author: "Hannah Li Hägi & Florentin Constancias "
date: " July 06, 2022 "
output: 
  html_document: 
    toc: yes
    keep_md: yes
---





```r
rm(list= ls())
gc()
```

```
##          used (Mb) gc trigger (Mb) limit (Mb) max used (Mb)
## Ncells 523766 28.0    1162022 62.1         NA   668902 35.8
## Vcells 974571  7.5    8388608 64.0    1024000  1839014 14.1
```

```r
library("tidyverse")
```

```
## ── Attaching packages ────────────────────────────────── tidyverse 1.3.1.9000 ──
```

```
## ✔ ggplot2 3.3.6     ✔ purrr   0.3.4
## ✔ tibble  3.1.7     ✔ dplyr   1.0.9
## ✔ tidyr   1.2.0     ✔ stringr 1.4.0
## ✔ readr   2.1.2     ✔ forcats 0.5.1
```

```
## ── Conflicts ────────────────────────────────────────── tidyverse_conflicts() ──
## ✖ dplyr::filter() masks stats::filter()
## ✖ dplyr::lag()    masks stats::lag()
```

```r
library('here')
```

```
## here() starts at /Users/fconstan/Documents/GitHub/NRP72-FBT
```

```r
'%!in%' <- function(x,y)!('%in%'(x,y))
```

# Define inputs:

## tables/ data:



```r
metadata <- "data/raw/25.11.2021_metadata_updated.tsv"
```



```r
metadata %>% 
  here::here() %>% 
  read_tsv() %>% 
  data.frame() %>%
  dplyr::filter(!is.na(metagenomic_sample_name), Model %in% c("Human", "Chicken")) -> meta
```

```
## Rows: 604 Columns: 61
## ── Column specification ────────────────────────────────────────────────────────
## Delimiter: "\t"
## chr (19): sample, Sample_description, I7_Index_ID, index, I5_Index_ID, index...
## dbl (42): input, filtered, denoisedF, denoisedR, merged, tabled, filtered_pc...
## 
## ℹ Use `spec()` to retrieve the full column specification for this data.
## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.
```

Modify metadata: **this should be done once on the metadata**


```r
# create Period column
meta %>%
  mutate(Treatment2 = factor(Treatment2, levels = c("DONOR", "UNTREATED", "AB+E. coli", "AB", "E. coli", "AB+E. faecium", "E. faecium")),
         Reactor_Treatment = base::replace(Reactor_Treatment, Reactor_Treatment == "TR_2CTX", "TR2_CTX"),
         Phase = factor(Phase, levels = c("Stab", "Treat")),
         Model = factor(Model, levels = c("Chicken", "Human")),
         Fermentation = factor(Fermentation, levels = c("1", "2")),
         Period = case_when( # if chicken we have period if human
           Day_of_Treatment == -2 | Day_of_Treatment == -3 | Day_of_Treatment == -6 | Day_of_Treatment == -7 ~ "pret",
           Day_of_Treatment < 10 & Day_of_Treatment >= 0 ~ "t1",
           Day_of_Treatment < 20 & Day_of_Treatment >= 10 ~ "t2",
           Day_of_Treatment < 30 & Day_of_Treatment >= 20 ~ "t3",
           Day_of_Treatment < 40 & Day_of_Treatment >= 30 ~ "t4",
           Day_of_Treatment < 50 & Day_of_Treatment >= 40 ~ "t5"),
         Period = factor(Period, levels = c("pret", "t1", "t2", "t3", "t4", "t5")),
         Reactor_Treatment_Dose = if_else(!is.na(`Antibiotic_mg.mL`), paste0(Reactor_Treatment, `Antibiotic_mg.mL`), Reactor_Treatment),
         Treatment_Dose = if_else(!is.na(`Antibiotic_mg.mL`), paste0(Treatment, `Antibiotic_mg.mL`), Treatment),
         Treatment = base::replace(Treatment, Experiment == "Cecum", "DONOR"),
         Treatment = factor(Treatment, levels = c("negative","positive","STRAIN", "DONOR", "UNTREATED", "CTX", "CTX+HV292.1", "HV292.1", "VAN", "VAN+CCUG59168", "CCUG59168")),
         # Reactor = factor(Reactor, levels = c("negative","positive","STRAIN", "DONOR", "IR1", "IR2", "CR", "TR1", "TR2", "TR3", "TR4", "TR5", "TR6","Batch3","Batch4")),
         Experiment = factor(Experiment, levels = c( "NTC","Mock", "HV292.1", "CCUG59168", "Cecum","Continuous","Batch")),
         
         # mutate( Treatment  = case_when(Experiment == "Cecum" ~ "DONOR",
         #                               TRUE ~ Treatment)) %>% 
         # mutate(Treatment = ifelse(Experiment == "Cecum", "Donor", Treatment)) %>% 
         
  ) %>% 
  # mutate(Day_of_Treatment = addNA(Day_of_Treatment)) %>%
  column_to_rownames("metagenomic_sample_name") -> meta

#human_meta$Treatment_Dose <- factor(human_meta$Treatment_Dose, ordered = TRUE, 
#                                 levels = c("DONOR", "UNTREATED", "CTX+HV292.120", "CTX+HV292.1200",
#                                            "CTX20", "CTX200", "HV292.1", "VAN+CCUG5916890", "VAN+CCUG59168600",
#                                            "VAN90", "VAN600", "CCUG59168"))
meta
```

```
##            sample  input filtered denoisedF denoisedR merged tabled filtered_pc
## D1     CR-15-S292  23558    23449     23410     23355  22970  22970       0.995
## D14     CR-17-S97  35320    35021     34992     34870  34419  34419       0.992
## D21    CR-19-S341  17616    17528     17501     17477  17246  17246       0.995
## D28    CR-21-S126   5791     5764      5753      5734   5585   5585       0.995
## D35    CR-46-S191  24284    24178     24128     24062  23531  23531       0.996
## D36    CR-64-S324  10847    10801     10784     10726  10465  10465       0.996
## D43      D-0-S154  11484    11428     11279     11256  10248  10248       0.995
## HC_1   TR1-61-S47 238630   238248    237846    237395 231166 231166       0.998
## HC_2  TR2-61-S196  24392    24363     24296     24257  22096  22096       0.999
## HC_3   TR3-61-S45  49196    49111     49012     48940  47064  47064       0.998
## HC_4   TR4-61-S11  28442    28406     28348     28277  27014  27014       0.999
## HC_5  TR5-61-S142    245      244       240       237    219    219       0.996
## HC_6  TR6-61-S186 252856   252414    252020    251515 245450 245450       0.998
## HC_7  TR1-63-S100     32       31        31        30     26     26       0.969
## HC_8  TR2-63-S160  36811    36755     36712     36601  35561  35561       0.998
## HC_9  TR3-63-S206  26658    26606     26549     26501  25740  25740       0.998
## HC_10  TR4-63-S74  41734    41656     41517     41407  39936  39936       0.998
## HC_11  TR5-63-S82  43222    43149     43078     42988  41415  41415       0.998
## D7    TR1-15-S168    153      151       144       115     79     79       0.987
## D8    TR1-17-S246  16895    16833     16807     16749  16429  16429       0.996
## D15   TR1-19-S305   9274     9256      9239      9203   8959   8959       0.998
## D22   TR1-21-S207  20083    20061     20043     19956  19419  19419       0.999
## D29   TR1-46-S170  32752    32656     32579     32459  31595  31595       0.997
## HC_12 TR6-63-S195  34779    34728     34672     34605  33178  33178       0.999
## D42   TR1-64-S263  31538    31375     31337     31248  30747  30747       0.995
## HC_13 TR1-65-S128     60       59        53        51     36     36       0.983
## HC_14 TR2-65-S158  43440    43349     43294     43167  41820  41820       0.998
## D6    TR2-15-S288  31566    31486     31385     31307  30381  30381       0.997
## D9    TR2-17-S248  31490    31404     31369     31205  30516  30516       0.997
## D16   TR2-19-S318  10065    10054     10040     10015   9814   9814       0.999
## D23   TR2-21-S228  15593    15567     15543     15456  14959  14959       0.998
## HC_15  TR3-65-S92  32513    32447     32308     32211  29961  29961       0.998
## HC_16 TR4-65-S155     27       26        25        24     15     15       0.963
## HC_17 TR5-65-S187     35       35        32        32     25     25       1.000
## HC_18  TR6-65-S94  33692    33627     33496     33425  30776  30776       0.998
## HC_19 TR1-69-S145  40431    40343     40200     40139  38039  38039       0.998
## HC_20  TR2-69-S57 178654   178284    177932    177646 174558 174558       0.998
## D30   TR2-46-S239  25929    25888     25862     25757  25098  25098       0.998
## HC_21 TR3-69-S139  30574    30520     30460     30405  29510  29510       0.998
## HC_22 TR4-69-S197  25982    25946     25869     25853  24646  24646       0.999
## HC_23 TR5-69-S107     23       23        22        22     19     19       1.000
## HC_24 TR6-69-S174  44451    44372     44282     44188  42841  42841       0.998
## D41   TR2-64-S216  21707    21663     21627     21526  20850  20850       0.998
## D5    TR3-15-S338  17747    17682     17655     17595  17213  17213       0.996
## D10   TR3-17-S222  14735    14652     14631     14571  14285  14285       0.994
## D17   TR3-19-S266  32415    32239     32191     32099  31521  31521       0.995
## D24   TR3-21-S234  23216    23054     23032     22946  22495  22495       0.993
## HV_1   TR1-28-S35  34188    34109     34047     34001  32656  32656       0.998
## HV_2   TR2-28-S64  25719    25668     25602     25529  23955  23955       0.998
## HV_3   TR3-28-S38 833850   832269    830661    829156 814071 814071       0.998
## HV_4  TR4-28-S138    119      119       112       111     94     94       1.000
## HV_5    TR5-28-S5  31023    30972     30868     30821  29212  29212       0.998
## HV_6   TR6-28-S80  28872    28808     28729     28664  27021  27021       0.998
## D31   TR3-46-S224  14862    14757     14738     14685  14391  14391       0.993
## D40   TR3-64-S190  28368    28173     28084     28039  27343  27343       0.993
## D4    TR4-15-S218  21839    21709     21687     21638  21317  21317       0.994
## D11   TR4-17-S114  28993    28752     28709     28624  28097  28097       0.992
## D18   TR4-19-S270  28080    27827     27811     27758  27315  27315       0.991
## D25   TR4-21-S254  21868    21735     21727     21684  21405  21405       0.994
## HV_7   TR1-30-S76  36509    36431     36343     36255  35484  35484       0.998
## HV_8   TR2-30-S73  38226    38119     37993     37920  36739  36739       0.997
## HV_9  TR3-30-S190  35109    35014     34904     34824  33704  33704       0.997
## HV_10  TR4-30-S48  26568    26495     26335     26274  24097  24097       0.997
## HV_11  TR5-30-S62  33028    32966     32872     32818  31102  31102       0.998
## HV_12 TR6-30-S144  31410    31362     31309     31213  29786  29786       0.998
## D32   TR4-46-S215  25294    25090     25053     24965  24729  24729       0.992
## D39   TR4-64-S220  15761    15686     15672     15607  15483  15483       0.995
## C_1   TR5-13-S103  25912    25816     25770     25681  24848  24848       0.996
## D3    TR5-15-S358  22887    22773     22768     22682  22378  22378       0.995
## D12    TR5-17-S82  81610    81154     81069     80828  79535  79535       0.994
## D19    TR5-19-S87  35069    34934     34910     34783  34135  34135       0.996
## C_6   TR5-20-S208  24461    24365     24357     24197  23598  23598       0.996
## D26   TR5-21-S185  35805    35626     35602     35513  34882  34882       0.995
## D33   TR5-46-S310  20221    20118     20093     19985  19794  19794       0.995
## HV_13 TR1-32-S111     32       32        30        27     23     23       1.000
## HV_14  TR2-32-S90  33510    33387     33329     33287  32683  32683       0.996
## HV_15 TR3-32-S211  24657    24582     24552     24494  24072  24072       0.997
## HV_16  TR4-32-S17  44815    44659     44563     44486  43879  43879       0.997
## HV_17  TR5-32-S56  23120    23079     23027     22987  21785  21785       0.998
## D38   TR5-64-S141  22049    21937     21916     21817  21659  21659       0.995
## HV_18  TR6-32-S39  30398    30364     30310     30240  28866  28866       0.999
## D_1      D-2-S164  29089    29028     28917     28759  27917  27917       0.998
## D2    TR6-15-S149  21435    21335     21298     21224  20816  20816       0.995
## D13   TR6-17-S138  20193    20062     20030     19961  19569  19569       0.994
## D20   TR6-19-S137  25289    25165     25122     24996  24337  24337       0.995
## D27   TR6-21-S107   4866     4830      4824      4806   4681   4681       0.993
## HV_19 TR1-36-S198  34485    34334     34287     33498  28774  28774       0.996
## HV_20 TR2-36-S183  46398    46202     46139     45199  39386  39386       0.996
## HV_21   TR3-36-S8  19943    19861     19755     19696  18108  18108       0.996
## HV_22  TR4-36-S12  31375    31261     31220     31194  30717  30717       0.996
## HV_23  TR5-36-S85  33570    33504     33364     33301  31257  31257       0.998
## D34   TR6-46-S176  49716    49525     49414     49275  48250  48250       0.996
## D37   TR6-64-S146  24475    24302     24256     24192  23730  23730       0.993
## HV_24 TR6-36-S124     37       37        34        32     27     27       1.000
## C_2        TR5-63     NA       NA        NA        NA     NA     NA          NA
## C_3        TR4-63     NA       NA        NA        NA     NA     NA          NA
## C_4        TR3-63     NA       NA        NA        NA     NA     NA          NA
## C_5        TR2-63     NA       NA        NA        NA     NA     NA          NA
##       denoisedF_pc denoisedR_pc merged_pc filtered_merged_pc input_merged_pc
## D1           0.998        0.996     0.981              0.980           0.975
## D14          0.999        0.996     0.984              0.983           0.974
## D21          0.998        0.997     0.985              0.984           0.979
## D28          0.998        0.995     0.971              0.969           0.964
## D35          0.998        0.995     0.975              0.973           0.969
## D36          0.998        0.993     0.970              0.969           0.965
## D43          0.987        0.985     0.909              0.897           0.892
## HC_1         0.998        0.996     0.972              0.970           0.969
## HC_2         0.997        0.996     0.909              0.907           0.906
## HC_3         0.998        0.997     0.960              0.958           0.957
## HC_4         0.998        0.995     0.953              0.951           0.950
## HC_5         0.984        0.971     0.912              0.898           0.894
## HC_6         0.998        0.996     0.974              0.972           0.971
## HC_7         1.000        0.968     0.839              0.839           0.812
## HC_8         0.999        0.996     0.969              0.968           0.966
## HC_9         0.998        0.996     0.970              0.967           0.966
## HC_10        0.997        0.994     0.962              0.959           0.957
## HC_11        0.998        0.996     0.961              0.960           0.958
## D7           0.954        0.762     0.549              0.523           0.516
## D8           0.998        0.995     0.978              0.976           0.972
## D15          0.998        0.994     0.970              0.968           0.966
## D22          0.999        0.995     0.969              0.968           0.967
## D29          0.998        0.994     0.970              0.968           0.965
## HC_12        0.998        0.996     0.957              0.955           0.954
## D42          0.999        0.996     0.981              0.980           0.975
## HC_13        0.898        0.864     0.679              0.610           0.600
## HC_14        0.999        0.996     0.966              0.965           0.963
## D6           0.997        0.994     0.968              0.965           0.962
## D9           0.999        0.994     0.973              0.972           0.969
## D16          0.999        0.996     0.977              0.976           0.975
## D23          0.998        0.993     0.962              0.961           0.959
## HC_15        0.996        0.993     0.927              0.923           0.922
## HC_16        0.962        0.923     0.600              0.577           0.556
## HC_17        0.914        0.914     0.781              0.714           0.714
## HC_18        0.996        0.994     0.919              0.915           0.913
## HC_19        0.996        0.995     0.946              0.943           0.941
## HC_20        0.998        0.996     0.981              0.979           0.977
## D30          0.999        0.995     0.970              0.969           0.968
## HC_21        0.998        0.996     0.969              0.967           0.965
## HC_22        0.997        0.996     0.953              0.950           0.949
## HC_23        0.957        0.957     0.864              0.826           0.826
## HC_24        0.998        0.996     0.967              0.965           0.964
## D41          0.998        0.994     0.964              0.962           0.961
## D5           0.998        0.995     0.975              0.973           0.970
## D10          0.999        0.994     0.976              0.975           0.969
## D17          0.999        0.996     0.979              0.978           0.972
## D24          0.999        0.995     0.977              0.976           0.969
## HV_1         0.998        0.997     0.959              0.957           0.955
## HV_2         0.997        0.995     0.936              0.933           0.931
## HV_3         0.998        0.996     0.980              0.978           0.976
## HV_4         0.941        0.933     0.839              0.790           0.790
## HV_5         0.997        0.995     0.946              0.943           0.942
## HV_6         0.997        0.995     0.941              0.938           0.936
## D31          0.999        0.995     0.976              0.975           0.968
## D40          0.997        0.995     0.974              0.971           0.964
## D4           0.999        0.997     0.983              0.982           0.976
## D11          0.999        0.996     0.979              0.977           0.969
## D18          0.999        0.998     0.982              0.982           0.973
## D25          1.000        0.998     0.985              0.985           0.979
## HV_7         0.998        0.995     0.976              0.974           0.972
## HV_8         0.997        0.995     0.967              0.964           0.961
## HV_9         0.997        0.995     0.966              0.963           0.960
## HV_10        0.994        0.992     0.915              0.909           0.907
## HV_11        0.997        0.996     0.946              0.943           0.942
## HV_12        0.998        0.995     0.951              0.950           0.948
## D32          0.999        0.995     0.987              0.986           0.978
## D39          0.999        0.995     0.988              0.987           0.982
## C_1          0.998        0.995     0.964              0.963           0.959
## D3           1.000        0.996     0.983              0.983           0.978
## D12          0.999        0.996     0.981              0.980           0.975
## D19          0.999        0.996     0.978              0.977           0.973
## C_6          1.000        0.993     0.969              0.969           0.965
## D26          0.999        0.997     0.980              0.979           0.974
## D33          0.999        0.993     0.985              0.984           0.979
## HV_13        0.938        0.844     0.767              0.719           0.719
## HV_14        0.998        0.997     0.981              0.979           0.975
## HV_15        0.999        0.996     0.980              0.979           0.976
## HV_16        0.998        0.996     0.985              0.983           0.979
## HV_17        0.998        0.996     0.946              0.944           0.942
## D38          0.999        0.995     0.988              0.987           0.982
## HV_18        0.998        0.996     0.952              0.951           0.950
## D_1          0.996        0.991     0.965              0.962           0.960
## D2           0.998        0.995     0.977              0.976           0.971
## D13          0.998        0.995     0.977              0.975           0.969
## D20          0.998        0.993     0.969              0.967           0.962
## D27          0.999        0.995     0.970              0.969           0.962
## HV_19        0.999        0.976     0.839              0.838           0.834
## HV_20        0.999        0.978     0.854              0.852           0.849
## HV_21        0.995        0.992     0.917              0.912           0.908
## HV_22        0.999        0.998     0.984              0.983           0.979
## HV_23        0.996        0.994     0.937              0.933           0.931
## D34          0.998        0.995     0.976              0.974           0.971
## D37          0.998        0.995     0.978              0.976           0.970
## HV_24        0.919        0.865     0.794              0.730           0.730
## C_2             NA           NA        NA                 NA              NA
## C_3             NA           NA        NA                 NA              NA
## C_4             NA           NA        NA                 NA              NA
## C_5             NA           NA        NA                 NA              NA
##       tabled_joined chimera_out length_filtered tabled_pc chimera_out_pc
## D1            22970       22102           22098         1           0.96
## D14           34419       32366           32350         1           0.94
## D21           17246       16865           16865         1           0.98
## D28            5585        5498            5498         1           0.98
## D35           23531       22823           22821         1           0.97
## D36           10465       10250           10250         1           0.98
## D43           10248       10134           10134         1           0.99
## HC_1         231166      188953          188948         1           0.82
## HC_2          22096       10235           10235         1           0.46
## HC_3          47064       36864           36862         1           0.78
## HC_4          27014       19750           19745         1           0.73
## HC_5            219         183             183         1           0.84
## HC_6         245450      201757          201745         1           0.82
## HC_7             26          24              24         1           0.92
## HC_8          35561       28863           28862         1           0.81
## HC_9          25740       20454           20454         1           0.79
## HC_10         39936       32773           32768         1           0.82
## HC_11         41415       31488           31486         1           0.76
## D7               79          79              79         1           1.00
## D8            16429       16233           16233         1           0.99
## D15            8959        8815            8815         1           0.98
## D22           19419       18970           18970         1           0.98
## D29           31595       30745           30745         1           0.97
## HC_12         33178       23813           23807         1           0.72
## D42           30747       30193           30193         1           0.98
## HC_13            36          29              29         1           0.81
## HC_14         41820       32880           32880         1           0.79
## D6            30381       28832           28832         1           0.95
## D9            30516       30071           30071         1           0.99
## D16            9814        9662            9662         1           0.98
## D23           14959       14630           14630         1           0.98
## HC_15         29961       17134           17133         1           0.57
## HC_16            15          14              14         1           0.93
## HC_17            25          23              23         1           0.92
## HC_18         30776       15573           15573         1           0.51
## HC_19         38039       24528           24524         1           0.64
## HC_20        174558      149775          149775         1           0.86
## D30           25098       24475           24475         1           0.98
## HC_21         29510       23544           23544         1           0.80
## HC_22         24646       15361           15359         1           0.62
## HC_23            19          17              17         1           0.89
## HC_24         42841       35074           35074         1           0.82
## D41           20850       20344           20344         1           0.98
## D5            17213       16817           16817         1           0.98
## D10           14285       14059           14059         1           0.98
## D17           31521       30971           30970         1           0.98
## D24           22495       21882           21876         1           0.97
## HV_1          32656       24739           24739         1           0.76
## HV_2          23955       14713           14713         1           0.61
## HV_3         814071      695560          695530         1           0.85
## HV_4             94          83              83         1           0.88
## HV_5          29212       18666           18666         1           0.64
## HV_6          27021       17152           17151         1           0.63
## D31           14391       14095           14095         1           0.98
## D40           27343       26205           26202         1           0.96
## D4            21317       20936           20936         1           0.98
## D11           28097       27386           27381         1           0.97
## D18           27315       26523           26523         1           0.97
## D25           21405       20811           20805         1           0.97
## HV_7          35484       28461           28461         1           0.80
## HV_8          36739       25840           25840         1           0.70
## HV_9          33704       26093           26093         1           0.77
## HV_10         24097       11843           11843         1           0.49
## HV_11         31102       21287           21287         1           0.68
## HV_12         29786       21835           21832         1           0.73
## D32           24729       23673           23673         1           0.96
## D39           15483       15160           15160         1           0.98
## C_1           24848       24314           24314         1           0.98
## D3            22378       21686           21680         1           0.97
## D12           79535       76689           76678         1           0.96
## D19           34135       33512           33510         1           0.98
## C_6           23598       22590           22589         1           0.96
## D26           34882       32585           32568         1           0.93
## D33           19794       19015           19015         1           0.96
## HV_13            23          18              18         1           0.78
## HV_14         32683       25240           25240         1           0.77
## HV_15         24072       19576           19576         1           0.81
## HV_16         43879       34419           34419         1           0.78
## HV_17         21785       15674           15673         1           0.72
## D38           21659       21080           21080         1           0.97
## HV_18         28866       21084           21084         1           0.73
## D_1           27917       26104           26104         1           0.94
## D2            20816       20429           20429         1           0.98
## D13           19569       19164           19164         1           0.98
## D20           24337       23661           23659         1           0.97
## D27            4681        4595            4595         1           0.98
## HV_19         28774       22660           22660         1           0.79
## HV_20         39386       33722           33722         1           0.86
## HV_21         18108        8748            8748         1           0.48
## HV_22         30717       23394           23394         1           0.76
## HV_23         31257       19614           19614         1           0.63
## D34           48250       46576           46575         1           0.97
## D37           23730       23096           23096         1           0.97
## HV_24            27          20              20         1           0.74
## C_2              NA          NA              NA        NA             NA
## C_3              NA          NA              NA        NA             NA
## C_4              NA          NA              NA        NA             NA
## C_5              NA          NA              NA        NA             NA
##       length_filtered_pc Sample_description I7_Index_ID    index I5_Index_ID
## D1                     1              CR-15      N716-D ACTCGCTA      S517-D
## D14                    1              CR-17      N701-A TAAGGCGA      S513-D
## D21                    1              CR-19      N723-D TAGCGCTC      S518-D
## D28                    1              CR-21      N704-A TCCTGAGC      S520-D
## D35                    1              CR-46      N715-A ATCTCAGG      S521-D
## D36                    1              CR-64      N721-D TACGCTGC      S517-D
## D43                    1                D-0      N710-A CGAGGCTG      S515-D
## HC_1                   1             TR1-61      A-N706     <NA>      A-S510
## HC_2                   1             TR2-61      D-N716     <NA>      A-S506
## HC_3                   1             TR3-61      A-N706     <NA>      A-S507
## HC_4                   1             TR4-61      A-N702     <NA>      A-S505
## HC_5                   1             TR5-61      A-N706     <NA>      D-S520
## HC_6                   1             TR6-61      A-N715     <NA>      D-S515
## HC_7                   1             TR1-63      A-N701     <NA>      D-S517
## HC_8                   1             TR2-63      A-N710     <NA>      D-S522
## HC_9                   1             TR3-63      D-N718     <NA>      A-S508
## HC_10                  1             TR4-63      A-N712     <NA>      A-S503
## HC_11                  1             TR5-63      A-N714     <NA>      A-S503
## D7                     1             TR1-15      N711-A AAGAGGCA      S522-D
## D8                     1             TR1-17      N723-D TAGCGCTC      S508-A
## D15                    1             TR1-19      N719-D GCGTAGTA      S513-D
## D22                    1             TR1-21      N718-D GGAGCTAC      S510-A
## D29                    1             TR1-46      N712-A GTAGAGGA      S515-D
## HC_12                  1             TR6-63      D-N716     <NA>      A-S505
## D42                    1             TR1-64      N726-D CCTAAGAC      S510-A
## HC_13                  1             TR1-65      A-N704     <NA>      D-S522
## HC_14                  1             TR2-65      A-N710     <NA>      D-S520
## D6                     1             TR2-15      N729-D TCGACGTC      S511-A
## D9                     1             TR2-17      N723-D TAGCGCTC      S511-A
## D16                    1             TR2-19      N720-D CGGAGCCT      S520-D
## D23                    1             TR2-21      N721-D TACGCTGC      S506-A
## HC_15                  1             TR3-65      A-N715     <NA>      A-S506
## HC_16                  1             TR4-65      A-N710     <NA>      D-S516
## HC_17                  1             TR5-65      A-N715     <NA>      D-S516
## HC_18                  1             TR6-65      A-N715     <NA>      A-S508
## HC_19                  1             TR1-69      A-N707     <NA>      D-S513
## HC_20                  1             TR2-69      A-N710     <NA>      A-S502
## D30                    1             TR2-46      N722-D ATGCGCAG      S510-A
## HC_21                  1             TR3-69      A-N706     <NA>      D-S516
## HC_22                  1             TR4-69      D-N716     <NA>      A-S507
## HC_23                  1             TR5-69      A-N702     <NA>      D-S516
## HC_24                  1             TR6-69      A-N712     <NA>      D-S520
## D41                    1             TR2-64      N719-D GCGTAGTA      S511-A
## D5                     1             TR3-15      N723-D TAGCGCTC      S515-D
## D10                    1             TR3-17      N720-D CGGAGCCT      S508-A
## D17                    1             TR3-19      N727-D CGATCAGT      S503-A
## D24                    1             TR3-21      N722-D ATGCGCAG      S503-A
## HV_1                   1             TR1-28      A-N705     <NA>      A-S505
## HV_2                   1             TR2-28      A-N710     <NA>      A-S511
## HV_3                   1             TR3-28      A-N705     <NA>      A-S508
## HV_4                   1             TR4-28      A-N706     <NA>      D-S515
## HV_5                   1             TR5-28      A-N701     <NA>      A-S507
## HV_6                   1             TR6-28      A-N712     <NA>      A-S511
## D31                    1             TR3-46      N720-D CGGAGCCT      S511-A
## D40                    1             TR3-64      N715-A ATCTCAGG      S520-D
## D4                     1             TR4-15      N720-D CGGAGCCT      S503-A
## D11                    1             TR4-17      N703-A AGGCAGAA      S515-D
## D18                    1             TR4-19      N727-D CGATCAGT      S508-A
## D25                    1             TR4-21      N724-D ACTGAGCG      S508-A
## HV_7                   1             TR1-30      A-N712     <NA>      A-S506
## HV_8                   1             TR2-30      A-N712     <NA>      A-S502
## HV_9                   1             TR3-30      A-N715     <NA>      D-S520
## HV_10                  1             TR4-30      A-N706     <NA>      A-S511
## HV_11                  1             TR5-30      A-N710     <NA>      A-S508
## HV_12                  1             TR6-30      A-N706     <NA>      D-S522
## D32                    1             TR4-46      N719-D GCGTAGTA      S510-A
## D39                    1             TR4-64      N720-D CGGAGCCT      S506-A
## C_1                    1             TR5-13      N701-A TAAGGCGA      S521-D
## D3                     1             TR5-15      N726-D CCTAAGAC      S520-D
## D12                    1             TR5-17      N714-A GCTCATGA      S503-A
## D19                    1             TR5-19      N714-A GCTCATGA      S510-A
## C_6                    1             TR5-20      N718-D GGAGCTAC      S511-A
## D26                    1             TR5-21      N715-A ATCTCAGG      S513-D
## D33                    1             TR5-46      N719-D GCGTAGTA      S520-D
## HV_13                  1             TR1-32      A-N702     <NA>      D-S521
## HV_14                  1             TR2-32      A-N715     <NA>      A-S503
## HV_15                  1             TR3-32      D-N719     <NA>      A-S505
## HV_16                  1             TR4-32      A-N703     <NA>      A-S502
## HV_17                  1             TR5-32      A-N707     <NA>      A-S511
## D38                    1             TR5-64      N706-A TAGGCATG      S518-D
## HV_18                  1             TR6-32      A-N705     <NA>      A-S510
## D_1                    1                D-2      A-N711     <NA>      D-S517
## D2                     1             TR6-15      N707-A CTCTCTAC      S518-D
## D13                    1             TR6-17      N706-A TAGGCATG      S515-D
## D20                    1             TR6-19      N706-A TAGGCATG      S513-D
## D27                    1             TR6-21      N702-A CGTACTAG      S516-D
## HV_19                  1             TR1-36      D-N716     <NA>      A-S508
## HV_20                  1             TR2-36      A-N714     <NA>      D-S521
## HV_21                  1             TR3-36      A-N701     <NA>      A-S511
## HV_22                  1             TR4-36      A-N702     <NA>      A-S506
## HV_23                  1             TR5-36      A-N714     <NA>      A-S507
## D34                    1             TR6-46      N712-A GTAGAGGA      S522-D
## D37                    1             TR6-64      N707-A CTCTCTAC      S515-D
## HV_24                  1             TR6-36      A-N704     <NA>      D-S517
## C_2                   NA             TR5-63        <NA>     <NA>        <NA>
## C_3                   NA             TR4-63        <NA>     <NA>        <NA>
## C_4                   NA             TR3-63        <NA>     <NA>        <NA>
## C_5                   NA             TR2-63        <NA>     <NA>        <NA>
##         index2 Description2 Experiment Reactor     Treatment Day_of_Connection
## D1    GCGTAAGA         <NA> Continuous      CR     UNTREATED                15
## D14   TCGACTAG         <NA> Continuous      CR     UNTREATED                17
## D21   CTATTAAG         <NA> Continuous      CR     UNTREATED                19
## D28   AAGGCTAT         <NA> Continuous      CR     UNTREATED                21
## D35   GAGCCTTA         <NA> Continuous      CR     UNTREATED                46
## D36   GCGTAAGA         <NA> Continuous      CR     UNTREATED                64
## D43   TTCTAGCT        DONOR      Cecum   DONOR         DONOR                NA
## HC_1      <NA>         <NA> Continuous     TR1   CTX+HV292.1                 8
## HC_2      <NA>         <NA> Continuous     TR2           CTX                 8
## HC_3      <NA>         <NA> Continuous     TR3   CTX+HV292.1                 8
## HC_4      <NA>         <NA> Continuous     TR4           CTX                 8
## HC_5      <NA>         <NA> Continuous     TR5       HV292.1                 8
## HC_6      <NA>         <NA> Continuous      CR     UNTREATED                 8
## HC_7      <NA>         <NA> Continuous     TR1   CTX+HV292.1                10
## HC_8      <NA>         <NA> Continuous     TR2           CTX                10
## HC_9      <NA>         <NA> Continuous     TR3   CTX+HV292.1                10
## HC_10     <NA>         <NA> Continuous     TR4           CTX                10
## HC_11     <NA>         <NA> Continuous     TR5       HV292.1                10
## D7    TTATGCGA         <NA> Continuous     TR1   CTX+HV292.1                15
## D8    CTAAGCCT         <NA> Continuous     TR1   CTX+HV292.1                17
## D15   TCGACTAG         <NA> Continuous     TR1   CTX+HV292.1                19
## D22   CGTCTAAT         <NA> Continuous     TR1   CTX+HV292.1                21
## D29   TTCTAGCT         <NA> Continuous     TR1   CTX+HV292.1                46
## HC_12     <NA>         <NA> Continuous      CR     UNTREATED                10
## D42   CGTCTAAT         <NA> Continuous     TR1   CTX+HV292.1                64
## HC_13     <NA>         <NA> Continuous     TR1   CTX+HV292.1                12
## HC_14     <NA>         <NA> Continuous     TR2           CTX                12
## D6    TCTCTCCG         <NA> Continuous     TR2           CTX                15
## D9    TCTCTCCG         <NA> Continuous     TR2           CTX                17
## D16   AAGGCTAT         <NA> Continuous     TR2           CTX                19
## D23   ACTGCATA         <NA> Continuous     TR2           CTX                21
## HC_15     <NA>         <NA> Continuous     TR3   CTX+HV292.1                12
## HC_16     <NA>         <NA> Continuous     TR4           CTX                12
## HC_17     <NA>         <NA> Continuous     TR5       HV292.1                12
## HC_18     <NA>         <NA> Continuous      CR     UNTREATED                12
## HC_19     <NA>         <NA> Continuous     TR1   CTX+HV292.1                16
## HC_20     <NA>         <NA> Continuous     TR2           CTX                16
## D30   CGTCTAAT         <NA> Continuous     TR2           CTX                46
## HC_21     <NA>         <NA> Continuous     TR3   CTX+HV292.1                16
## HC_22     <NA>         <NA> Continuous     TR4           CTX                16
## HC_23     <NA>         <NA> Continuous     TR5       HV292.1                16
## HC_24     <NA>         <NA> Continuous      CR     UNTREATED                16
## D41   TCTCTCCG         <NA> Continuous     TR2           CTX                64
## D5    TTCTAGCT         <NA> Continuous     TR3       HV292.1                15
## D10   CTAAGCCT         <NA> Continuous     TR3       HV292.1                17
## D17   TATCCTCT         <NA> Continuous     TR3       HV292.1                19
## D24   TATCCTCT         <NA> Continuous     TR3       HV292.1                21
## HV_1      <NA>         <NA> Continuous     TR1 VAN+CCUG59168                28
## HV_2      <NA>         <NA> Continuous     TR2           VAN                28
## HV_3      <NA>         <NA> Continuous     TR3 VAN+CCUG59168                28
## HV_4      <NA>         <NA> Continuous     TR4           VAN                28
## HV_5      <NA>         <NA> Continuous     TR5     CCUG59168                28
## HV_6      <NA>         <NA> Continuous      CR     UNTREATED                28
## D31   TCTCTCCG         <NA> Continuous     TR3       HV292.1                46
## D40   AAGGCTAT         <NA> Continuous     TR3       HV292.1                64
## D4    TATCCTCT         <NA> Continuous     TR4           VAN                15
## D11   TTCTAGCT         <NA> Continuous     TR4           VAN                17
## D18   CTAAGCCT         <NA> Continuous     TR4           VAN                19
## D25   CTAAGCCT         <NA> Continuous     TR4           VAN                21
## HV_7      <NA>         <NA> Continuous     TR1 VAN+CCUG59168                30
## HV_8      <NA>         <NA> Continuous     TR2           VAN                30
## HV_9      <NA>         <NA> Continuous     TR3 VAN+CCUG59168                30
## HV_10     <NA>         <NA> Continuous     TR4           VAN                30
## HV_11     <NA>         <NA> Continuous     TR5     CCUG59168                30
## HV_12     <NA>         <NA> Continuous      CR     UNTREATED                30
## D32   CGTCTAAT         <NA> Continuous     TR4           VAN                46
## D39   ACTGCATA         <NA> Continuous     TR4           VAN                64
## C_1   GAGCCTTA         <NA> Continuous     TR5 VAN+CCUG59168                13
## D3    AAGGCTAT         <NA> Continuous     TR5 VAN+CCUG59168                15
## D12   TATCCTCT         <NA> Continuous     TR5 VAN+CCUG59168                17
## D19   CGTCTAAT         <NA> Continuous     TR5 VAN+CCUG59168                19
## C_6   TCTCTCCG         <NA> Continuous     TR5 VAN+CCUG59168                20
## D26   TCGACTAG         <NA> Continuous     TR5 VAN+CCUG59168                21
## D33   AAGGCTAT         <NA> Continuous     TR5 VAN+CCUG59168                46
## HV_13     <NA>         <NA> Continuous     TR1 VAN+CCUG59168                32
## HV_14     <NA>         <NA> Continuous     TR2           VAN                32
## HV_15     <NA>         <NA> Continuous     TR3 VAN+CCUG59168                32
## HV_16     <NA>         <NA> Continuous     TR4           VAN                32
## HV_17     <NA>         <NA> Continuous     TR5     CCUG59168                32
## D38   CTATTAAG         <NA> Continuous     TR5 VAN+CCUG59168                64
## HV_18     <NA>         <NA> Continuous      CR     UNTREATED                32
## D_1       <NA>        DONOR      Cecum   DONOR         DONOR                NA
## D2    CTATTAAG         <NA> Continuous     TR6     CCUG59168                15
## D13   TTCTAGCT         <NA> Continuous     TR6     CCUG59168                17
## D20   TCGACTAG         <NA> Continuous     TR6     CCUG59168                19
## D27   CCTAGAGT         <NA> Continuous     TR6     CCUG59168                21
## HV_19     <NA>         <NA> Continuous     TR1 VAN+CCUG59168                36
## HV_20     <NA>         <NA> Continuous     TR2           VAN                36
## HV_21     <NA>         <NA> Continuous     TR3 VAN+CCUG59168                36
## HV_22     <NA>         <NA> Continuous     TR4           VAN                36
## HV_23     <NA>         <NA> Continuous     TR5     CCUG59168                36
## D34   TTATGCGA         <NA> Continuous     TR6     CCUG59168                46
## D37   TTCTAGCT         <NA> Continuous     TR6     CCUG59168                64
## HV_24     <NA>         <NA> Continuous      CR     UNTREATED                36
## C_2       <NA>         <NA> Continuous     TR5 VAN+CCUG59168                63
## C_3       <NA>         <NA> Continuous     TR4           VAN                63
## C_4       <NA>         <NA> Continuous     TR3       HV292.1                63
## C_5       <NA>         <NA> Continuous     TR2           CTX                63
##       Day_of_Treatment Day_from_Inoculum  Enrichment Phase    Treatment2
## D1                  -1                38 NotEnriched  Stab     UNTREATED
## D14                  1                40 NotEnriched Treat     UNTREATED
## D21                  3                42 NotEnriched Treat     UNTREATED
## D28                  5                44 NotEnriched Treat     UNTREATED
## D35                 30                69 NotEnriched Treat     UNTREATED
## D36                 48                87 NotEnriched Treat     UNTREATED
## D43                 NA                NA NotEnriched  <NA>         DONOR
## HC_1                -1                76 NotEnriched  Stab    AB+E. coli
## HC_2                -1                76 NotEnriched  Stab            AB
## HC_3                -1                76 NotEnriched  Stab    AB+E. coli
## HC_4                -1                76 NotEnriched  Stab            AB
## HC_5                -1                76 NotEnriched  Stab       E. coli
## HC_6                -1                76 NotEnriched  Stab     UNTREATED
## HC_7                 1                78 NotEnriched Treat    AB+E. coli
## HC_8                 1                78 NotEnriched Treat            AB
## HC_9                 1                78 NotEnriched Treat    AB+E. coli
## HC_10                1                78 NotEnriched Treat            AB
## HC_11                1                78 NotEnriched Treat       E. coli
## D7                  -1                38 NotEnriched  Stab    AB+E. coli
## D8                   1                40 NotEnriched Treat    AB+E. coli
## D15                  3                42 NotEnriched Treat    AB+E. coli
## D22                  5                44 NotEnriched Treat    AB+E. coli
## D29                 30                69 NotEnriched Treat    AB+E. coli
## HC_12                1                78 NotEnriched Treat     UNTREATED
## D42                 48                87 NotEnriched Treat    AB+E. coli
## HC_13                3                80 NotEnriched Treat    AB+E. coli
## HC_14                3                80 NotEnriched Treat            AB
## D6                  -1                38 NotEnriched  Stab            AB
## D9                   1                40 NotEnriched Treat            AB
## D16                  3                42 NotEnriched Treat            AB
## D23                  5                44 NotEnriched Treat            AB
## HC_15                3                80 NotEnriched Treat    AB+E. coli
## HC_16                3                80 NotEnriched Treat            AB
## HC_17                3                80 NotEnriched Treat       E. coli
## HC_18                3                80 NotEnriched Treat     UNTREATED
## HC_19                7                84 NotEnriched Treat    AB+E. coli
## HC_20                7                84 NotEnriched Treat            AB
## D30                 30                69 NotEnriched Treat            AB
## HC_21                7                84 NotEnriched Treat    AB+E. coli
## HC_22                7                84 NotEnriched Treat            AB
## HC_23                7                84 NotEnriched Treat       E. coli
## HC_24                7                84 NotEnriched Treat     UNTREATED
## D41                 48                87 NotEnriched Treat            AB
## D5                  -1                38 NotEnriched  Stab       E. coli
## D10                  1                40 NotEnriched Treat       E. coli
## D17                  3                42 NotEnriched Treat       E. coli
## D24                  5                44 NotEnriched Treat       E. coli
## HV_1                -1                43 NotEnriched  Stab AB+E. faecium
## HV_2                -1                43 NotEnriched  Stab            AB
## HV_3                -1                43 NotEnriched  Stab AB+E. faecium
## HV_4                -1                43 NotEnriched  Stab            AB
## HV_5                -1                43 NotEnriched  Stab    E. faecium
## HV_6                -1                43 NotEnriched  Stab     UNTREATED
## D31                 30                69 NotEnriched Treat       E. coli
## D40                 48                87 NotEnriched Treat       E. coli
## D4                  -1                38 NotEnriched  Stab            AB
## D11                  1                40 NotEnriched Treat            AB
## D18                  3                42 NotEnriched Treat            AB
## D25                  5                44 NotEnriched Treat            AB
## HV_7                 1                45 NotEnriched Treat AB+E. faecium
## HV_8                 1                45 NotEnriched Treat            AB
## HV_9                 1                45 NotEnriched Treat AB+E. faecium
## HV_10                1                45 NotEnriched Treat            AB
## HV_11                1                45 NotEnriched Treat    E. faecium
## HV_12                1                45 NotEnriched Treat     UNTREATED
## D32                 30                69 NotEnriched Treat            AB
## D39                 48                87 NotEnriched Treat            AB
## C_1                 -3                36 NotEnriched  Stab AB+E. faecium
## D3                  -1                38 NotEnriched  Stab AB+E. faecium
## D12                  1                40 NotEnriched Treat AB+E. faecium
## D19                  3                42 NotEnriched Treat AB+E. faecium
## C_6                  4                43 NotEnriched Treat AB+E. faecium
## D26                  5                44 NotEnriched Treat AB+E. faecium
## D33                 30                69 NotEnriched Treat AB+E. faecium
## HV_13                3                47 NotEnriched Treat AB+E. faecium
## HV_14                3                47 NotEnriched Treat            AB
## HV_15                3                47 NotEnriched Treat AB+E. faecium
## HV_16                3                47 NotEnriched Treat            AB
## HV_17                3                47 NotEnriched Treat    E. faecium
## D38                 48                87 NotEnriched Treat AB+E. faecium
## HV_18                3                47 NotEnriched Treat     UNTREATED
## D_1                 NA                NA NotEnriched  <NA>         DONOR
## D2                  -1                38 NotEnriched  Stab    E. faecium
## D13                  1                40 NotEnriched Treat    E. faecium
## D20                  3                42 NotEnriched Treat    E. faecium
## D27                  5                44 NotEnriched Treat    E. faecium
## HV_19                7                51 NotEnriched Treat AB+E. faecium
## HV_20                7                51 NotEnriched Treat            AB
## HV_21                7                51 NotEnriched Treat AB+E. faecium
## HV_22                7                51 NotEnriched Treat            AB
## HV_23                7                51 NotEnriched Treat    E. faecium
## D34                 30                69 NotEnriched Treat    E. faecium
## D37                 48                87 NotEnriched Treat    E. faecium
## HV_24                7                51 NotEnriched Treat     UNTREATED
## C_2                 47                86 NotEnriched Treat AB+E. faecium
## C_3                 47                86 NotEnriched Treat            AB
## C_4                 47                86 NotEnriched Treat       E. coli
## C_5                 47                86 NotEnriched Treat            AB
##                       Date Paul Reactor_Treatment GeneCopyNumberperML
## D1    2020-06-26T00:00:00Z <NA>      CR_UNTREATED            6.75e+10
## D14   2020-06-28T00:00:00Z <NA>      CR_UNTREATED            7.01e+10
## D21   2020-06-30T00:00:00Z <NA>      CR_UNTREATED            9.01e+10
## D28   2020-07-02T00:00:00Z <NA>      CR_UNTREATED            4.39e+10
## D35   2020-07-27T00:00:00Z <NA>      CR_UNTREATED            4.53e+10
## D36   2020-08-14T00:00:00Z <NA>      CR_UNTREATED            3.26e+10
## D43   2020-05-19T00:00:00Z <NA>             DONOR            1.44e+10
## HC_1                  <NA> <NA>   TR1_CTX+HV292.1            2.14e+11
## HC_2                  <NA> <NA>           TR2_CTX            2.07e+11
## HC_3                  <NA> <NA>   TR3_CTX+HV292.1            1.98e+11
## HC_4                  <NA> <NA>           TR4_CTX            2.05e+11
## HC_5                  <NA> <NA>       TR5_HV292.1            2.43e+11
## HC_6                  <NA> <NA>      CR_UNTREATED            2.37e+11
## HC_7                  <NA> <NA>   TR1_CTX+HV292.1            5.59e+11
## HC_8                  <NA> <NA>           TR2_CTX            2.66e+11
## HC_9                  <NA> <NA>   TR3_CTX+HV292.1            3.83e+11
## HC_10                 <NA> <NA>           TR4_CTX            3.22e+11
## HC_11                 <NA> <NA>       TR5_HV292.1            4.52e+11
## D7    2020-06-26T00:00:00Z <NA>   TR1_CTX+HV292.1            5.08e+10
## D8    2020-06-28T00:00:00Z <NA>   TR1_CTX+HV292.1            3.11e+10
## D15   2020-06-30T00:00:00Z <NA>   TR1_CTX+HV292.1            3.36e+10
## D22   2020-07-02T00:00:00Z <NA>   TR1_CTX+HV292.1            2.22e+10
## D29   2020-07-27T00:00:00Z <NA>   TR1_CTX+HV292.1            6.05e+10
## HC_12                 <NA> <NA>      CR_UNTREATED            2.39e+11
## D42   2020-08-14T00:00:00Z <NA>   TR1_CTX+HV292.1            3.61e+10
## HC_13                 <NA> <NA>   TR1_CTX+HV292.1            1.36e+11
## HC_14                 <NA> <NA>           TR2_CTX            2.85e+11
## D6    2020-06-26T00:00:00Z <NA>           TR2_CTX            6.23e+10
## D9    2020-06-28T00:00:00Z <NA>           TR2_CTX            2.71e+10
## D16   2020-06-30T00:00:00Z <NA>           TR2_CTX            4.20e+10
## D23   2020-07-02T00:00:00Z <NA>           TR2_CTX            2.89e+10
## HC_15                 <NA> <NA>   TR3_CTX+HV292.1            2.04e+11
## HC_16                 <NA> <NA>           TR4_CTX            2.39e+11
## HC_17                 <NA> <NA>       TR5_HV292.1            2.54e+11
## HC_18                 <NA> <NA>      CR_UNTREATED            1.26e+11
## HC_19                 <NA> <NA>   TR1_CTX+HV292.1            2.10e+11
## HC_20                 <NA> <NA>           TR2_CTX            1.02e+11
## D30   2020-07-27T00:00:00Z <NA>           TR2_CTX            5.05e+10
## HC_21                 <NA> <NA>   TR3_CTX+HV292.1            2.41e+11
## HC_22                 <NA> <NA>           TR4_CTX            2.42e+11
## HC_23                 <NA> <NA>       TR5_HV292.1            3.93e+11
## HC_24                 <NA> <NA>      CR_UNTREATED            2.60e+11
## D41   2020-08-14T00:00:00Z <NA>           TR2_CTX            2.79e+10
## D5    2020-06-26T00:00:00Z <NA>       TR3_HV292.1            4.87e+10
## D10   2020-06-28T00:00:00Z <NA>       TR3_HV292.1            5.24e+10
## D17   2020-06-30T00:00:00Z <NA>       TR3_HV292.1            4.62e+10
## D24   2020-07-02T00:00:00Z <NA>       TR3_HV292.1            4.79e+10
## HV_1                  <NA> <NA> TR1_VAN+CCUG59168            1.11e+11
## HV_2                  <NA> <NA>           TR2_VAN            2.50e+11
## HV_3                  <NA> <NA> TR3_VAN+CCUG59168            3.13e+11
## HV_4                  <NA> <NA>           TR4_VAN            2.99e+10
## HV_5                  <NA> <NA>     TR5_CCUG59168            1.92e+11
## HV_6                  <NA> <NA>      CR_UNTREATED            2.58e+11
## D31   2020-07-27T00:00:00Z <NA>       TR3_HV292.1            3.84e+10
## D40   2020-08-14T00:00:00Z <NA>       TR3_HV292.1            3.43e+10
## D4    2020-06-26T00:00:00Z <NA>           TR4_VAN            6.68e+10
## D11   2020-06-28T00:00:00Z <NA>           TR4_VAN            3.74e+10
## D18   2020-06-30T00:00:00Z <NA>           TR4_VAN            5.61e+10
## D25   2020-07-02T00:00:00Z <NA>           TR4_VAN            5.82e+10
## HV_7                  <NA> <NA> TR1_VAN+CCUG59168            3.03e+11
## HV_8                  <NA> <NA>           TR2_VAN            7.45e+10
## HV_9                  <NA> <NA> TR3_VAN+CCUG59168            1.21e+11
## HV_10                 <NA> <NA>           TR4_VAN            1.32e+11
## HV_11                 <NA> <NA>     TR5_CCUG59168            2.77e+11
## HV_12                 <NA> <NA>      CR_UNTREATED            4.92e+10
## D32   2020-07-27T00:00:00Z <NA>           TR4_VAN            6.31e+10
## D39   2020-08-14T00:00:00Z <NA>           TR4_VAN            6.94e+10
## C_1   2020-06-24T00:00:00Z <NA> TR5_VAN+CCUG59168            3.94e+10
## D3    2020-06-26T00:00:00Z <NA> TR5_VAN+CCUG59168            7.15e+10
## D12   2020-06-28T00:00:00Z <NA> TR5_VAN+CCUG59168            2.60e+10
## D19   2020-06-30T00:00:00Z <NA> TR5_VAN+CCUG59168            5.71e+10
## C_6   2020-07-01T00:00:00Z <NA> TR5_VAN+CCUG59168            2.50e+10
## D26   2020-07-02T00:00:00Z <NA> TR5_VAN+CCUG59168            4.89e+10
## D33   2020-07-27T00:00:00Z <NA> TR5_VAN+CCUG59168            6.59e+10
## HV_13                 <NA> <NA> TR1_VAN+CCUG59168            1.94e+11
## HV_14                 <NA> <NA>           TR2_VAN            2.43e+11
## HV_15                 <NA> <NA> TR3_VAN+CCUG59168            2.00e+11
## HV_16                 <NA> <NA>           TR4_VAN            9.63e+10
## HV_17                 <NA> <NA>     TR5_CCUG59168            2.00e+11
## D38   2020-08-14T00:00:00Z <NA> TR5_VAN+CCUG59168            5.32e+10
## HV_18                 <NA> <NA>      CR_UNTREATED            3.02e+11
## D_1                   <NA> <NA>             DONOR            1.90e+10
## D2    2020-06-26T00:00:00Z <NA>     TR6_CCUG59168            5.30e+10
## D13   2020-06-28T00:00:00Z <NA>     TR6_CCUG59168            3.93e+10
## D20   2020-06-30T00:00:00Z <NA>     TR6_CCUG59168            5.68e+10
## D27   2020-07-02T00:00:00Z <NA>     TR6_CCUG59168            6.93e+10
## HV_19                 <NA> <NA> TR1_VAN+CCUG59168            4.92e+10
## HV_20                 <NA> <NA>           TR2_VAN            3.21e+11
## HV_21                 <NA> <NA> TR3_VAN+CCUG59168            2.91e+11
## HV_22                 <NA> <NA>           TR4_VAN            1.17e+11
## HV_23                 <NA> <NA>     TR5_CCUG59168            1.68e+11
## D34   2020-07-27T00:00:00Z <NA>     TR6_CCUG59168            2.81e+10
## D37   2020-08-14T00:00:00Z <NA>     TR6_CCUG59168            3.93e+10
## HV_24                 <NA> <NA>      CR_UNTREATED            4.39e+11
## C_2               13.08.20 <NA> TR5_VAN+CCUG59168                  NA
## C_3               13.08.20 <NA>           TR4_VAN                  NA
## C_4               13.08.20 <NA>       TR3_HV292.1                  NA
## C_5               13.08.20 <NA>           TR2_CTX                  NA
##       HV292.1_Copy_Number_permL CCUG59168_Copy_Number_permL
## D1                           NA                          NA
## D14                          NA                          NA
## D21                          NA                          NA
## D28                          NA                          NA
## D35                          NA                          NA
## D36                          NA                          NA
## D43                          NA                          NA
## HC_1                         NA                          NA
## HC_2                         NA                          NA
## HC_3                         NA                          NA
## HC_4                         NA                          NA
## HC_5                         NA                          NA
## HC_6                         NA                          NA
## HC_7                    5520000                          NA
## HC_8                          0                          NA
## HC_9                   30487500                          NA
## HC_10                         0                          NA
## HC_11                  20600000                          NA
## D7                           NA                          NA
## D8                           NA                          NA
## D15                      794000                          NA
## D22                          NA                          NA
## D29                    54054532                          NA
## HC_12                         0                          NA
## D42                   862908750                          NA
## HC_13                         0                          NA
## HC_14                         0                          NA
## D6                           NA                          NA
## D9                           NA                          NA
## D16                          NA                          NA
## D23                          NA                          NA
## HC_15                         0                          NA
## HC_16                         0                          NA
## HC_17                         0                          NA
## HC_18                         0                          NA
## HC_19                         0                          NA
## HC_20                         0                          NA
## D30                          NA                          NA
## HC_21                         0                          NA
## HC_22                         0                          NA
## HC_23                         0                          NA
## HC_24                         0                          NA
## D41                          NA                          NA
## D5                           NA                          NA
## D10                          NA                          NA
## D17                      887000                          NA
## D24                          NA                          NA
## HV_1                         NA                          NA
## HV_2                         NA                          NA
## HV_3                         NA                          NA
## HV_4                         NA                          NA
## HV_5                         NA                          NA
## HV_6                         NA                          NA
## D31                           0                          NA
## D40                           0                          NA
## D4                           NA                          NA
## D11                          NA                          NA
## D18                          NA                          NA
## D25                          NA                          NA
## HV_7                         NA                1.380000e+10
## HV_8                         NA                0.000000e+00
## HV_9                         NA                2.080000e+09
## HV_10                        NA                0.000000e+00
## HV_11                        NA                1.570000e+07
## HV_12                        NA                0.000000e+00
## D32                          NA                          NA
## D39                          NA                          NA
## C_1                          NA                          NA
## D3                           NA                          NA
## D12                          NA                          NA
## D19                          NA                5.040000e+07
## C_6                          NA                          NA
## D26                          NA                          NA
## D33                          NA                9.198400e+07
## HV_13                        NA                8.450000e+09
## HV_14                        NA                0.000000e+00
## HV_15                        NA                2.690000e+09
## HV_16                        NA                0.000000e+00
## HV_17                        NA                8.830000e+05
## D38                          NA                5.390000e+08
## HV_18                        NA                0.000000e+00
## D_1                          NA                          NA
## D2                           NA                          NA
## D13                          NA                          NA
## D20                          NA                0.000000e+00
## D27                          NA                          NA
## HV_19                        NA                5.165000e+09
## HV_20                        NA                0.000000e+00
## HV_21                        NA                1.670000e+09
## HV_22                        NA                0.000000e+00
## HV_23                        NA                1.290000e+05
## D34                          NA                2.688168e+05
## D37                          NA                2.120000e+03
## HV_24                        NA                0.000000e+00
## C_2                          NA                          NA
## C_3                          NA                          NA
## C_4                          NA                          NA
## C_5                          NA                          NA
##       CTX_Copy_Number_permL VAN_Copy_Number_permL   Model Antibiotic_mg.mL
## D1                       NA                    NA Chicken               NA
## D14                      NA                    NA Chicken               NA
## D21                      NA                    NA Chicken               NA
## D28                      NA                    NA Chicken               NA
## D35                      NA                    NA Chicken               NA
## D36                      NA                    NA Chicken               NA
## D43                      NA                    NA Chicken               NA
## HC_1                     NA                    NA   Human               20
## HC_2                     NA                    NA   Human               20
## HC_3                     NA                    NA   Human              200
## HC_4                     NA                    NA   Human              200
## HC_5                     NA                    NA   Human               NA
## HC_6                     NA                    NA   Human               NA
## HC_7           1.102500e+07                    NA   Human               20
## HC_8           0.000000e+00                    NA   Human               20
## HC_9           5.373750e+07                    NA   Human              200
## HC_10          0.000000e+00                    NA   Human              200
## HC_11          4.160000e+07                    NA   Human               NA
## D7                       NA                    NA Chicken               20
## D8                       NA                    NA Chicken               20
## D15            1.507417e+06                    NA Chicken               20
## D22                      NA                    NA Chicken               20
## D29            8.070878e+07                    NA Chicken               20
## HC_12          0.000000e+00                    NA   Human               NA
## D42            1.927668e+09                    NA Chicken               20
## HC_13          2.160000e+06                    NA   Human               20
## HC_14          0.000000e+00                    NA   Human               20
## D6                       NA                    NA Chicken               20
## D9                       NA                    NA Chicken               20
## D16                      NA                    NA Chicken               20
## D23                      NA                    NA Chicken               20
## HC_15          3.550000e+06                    NA   Human              200
## HC_16          0.000000e+00                    NA   Human              200
## HC_17          2.100000e+05                    NA   Human               NA
## HC_18          0.000000e+00                    NA   Human               NA
## HC_19          0.000000e+00                    NA   Human               20
## HC_20          0.000000e+00                    NA   Human               20
## D30                      NA                    NA Chicken               20
## HC_21          1.117500e+06                    NA   Human              200
## HC_22          0.000000e+00                    NA   Human              200
## HC_23                    NA                    NA   Human               NA
## HC_24          0.000000e+00                    NA   Human               NA
## D41                      NA                    NA Chicken               20
## D5                       NA                    NA Chicken               NA
## D10                      NA                    NA Chicken               NA
## D17            1.670000e+06                    NA Chicken               NA
## D24                      NA                    NA Chicken               NA
## HV_1                     NA                    NA   Human               90
## HV_2                     NA                    NA   Human               90
## HV_3                     NA                    NA   Human              600
## HV_4                     NA                    NA   Human              600
## HV_5                     NA                    NA   Human               NA
## HV_6                     NA                    NA   Human               NA
## D31            4.713665e+03                    NA Chicken               NA
## D40            3.666540e+04                    NA Chicken               NA
## D4                       NA                    NA Chicken               90
## D11                      NA                    NA Chicken               90
## D18                      NA                    NA Chicken               90
## D25                      NA                    NA Chicken               90
## HV_7                     NA          9.600000e+12   Human               90
## HV_8                     NA          0.000000e+00   Human               90
## HV_9                     NA          8.500000e+11   Human              600
## HV_10                    NA          0.000000e+00   Human              600
## HV_11                    NA          6.125000e+09   Human               NA
## HV_12                    NA          0.000000e+00   Human               NA
## D32                      NA                    NA Chicken               90
## D39                      NA                    NA Chicken               90
## C_1                      NA                    NA Chicken               90
## D3                       NA          1.740250e+10 Chicken               90
## D12                      NA          4.156488e+10 Chicken               90
## D19                      NA          6.682795e+10 Chicken               90
## C_6                      NA          2.847335e+10 Chicken               90
## D26                      NA          5.769750e+10 Chicken               90
## D33                      NA          5.770000e+11 Chicken               90
## HV_13                    NA          8.400000e+12   Human               90
## HV_14                    NA          0.000000e+00   Human               90
## HV_15                    NA          1.190000e+12   Human              600
## HV_16                    NA          0.000000e+00   Human              600
## HV_17                    NA          3.470000e+08   Human               NA
## D38                      NA                    NA Chicken               90
## HV_18                    NA          0.000000e+00   Human               NA
## D_1                      NA                    NA   Human               NA
## D2                       NA          1.065350e+08 Chicken               NA
## D13                      NA          0.000000e+00 Chicken               NA
## D20                      NA          0.000000e+00 Chicken               NA
## D27                      NA          0.000000e+00 Chicken               NA
## HV_19                    NA          1.300000e+12   Human               90
## HV_20                    NA          0.000000e+00   Human               90
## HV_21                    NA          4.800000e+11   Human              600
## HV_22                    NA          0.000000e+00   Human              600
## HV_23                    NA          0.000000e+00   Human               NA
## D34                      NA                    NA Chicken               NA
## D37                      NA                    NA Chicken               NA
## HV_24                    NA          0.000000e+00   Human               NA
## C_2                      NA                    NA Chicken               90
## C_3                      NA                    NA Chicken               90
## C_4                      NA                    NA Chicken               NA
## C_5                      NA                    NA Chicken               20
##       Fermentation Antibiotic Lactose_mM Glucose_mM Galactose_mM Succinat_mM
## D1            <NA>       <NA>      0.000      0.947        0.000       1.969
## D14           <NA>       <NA>      0.000      0.000        1.759       4.295
## D21           <NA>       <NA>      0.261      0.000        0.000       0.818
## D28           <NA>       <NA>      0.000      0.000        0.000       0.701
## D35           <NA>       <NA>         NA         NA           NA          NA
## D36           <NA>       <NA>      0.000      0.000        0.000       0.546
## D43           <NA>       <NA>         NA         NA           NA          NA
## HC_1             2        CTX      0.000      0.000        0.536       1.519
## HC_2             2        CTX      0.000      0.000        0.588       2.608
## HC_3             2        CTX      0.000      0.000        0.655       1.531
## HC_4             2        CTX      0.000      0.000        0.622       1.513
## HC_5             2        CTX      0.000      0.000        0.602       1.448
## HC_6             2        CTX      0.000      0.000        0.703       1.356
## HC_7             2        CTX      0.000      0.000        0.314       1.558
## HC_8             2        CTX      0.000      0.000        0.703       1.497
## HC_9             2        CTX      0.000      0.305        0.989       1.200
## HC_10            2        CTX      0.000      0.309        0.958       1.300
## HC_11            2        CTX      0.000      0.000        0.680       1.435
## D7            <NA>        CTX      0.000      0.000        0.000       1.090
## D8            <NA>        CTX      0.000      0.000        0.000       1.201
## D15           <NA>        CTX      0.000      0.000        0.000       2.779
## D22           <NA>        CTX      0.000      0.000        0.000       3.069
## D29           <NA>        CTX      0.000      0.000        0.000       0.527
## HC_12            2        CTX      0.000      0.000        0.464       1.390
## D42           <NA>        CTX      0.000      0.000        0.000       0.455
## HC_13            2        CTX      0.000      0.000        0.367       1.370
## HC_14            2        CTX      0.000      0.000        0.000       1.085
## D6            <NA>        CTX      0.000      0.000        0.000       0.792
## D9            <NA>        CTX      0.000      0.000        0.000       1.116
## D16           <NA>        CTX      0.000      0.000        0.000       2.973
## D23           <NA>        CTX      0.000      0.000        0.000       3.420
## HC_15            2        CTX      0.000      0.000        0.391       1.100
## HC_16            2        CTX      0.000      0.000        0.397       1.051
## HC_17            2        CTX      0.000      0.000        0.000       1.110
## HC_18            2        CTX      0.000      0.000        0.000       1.258
## HC_19            2        CTX      0.000      0.000        0.229       1.194
## HC_20            2        CTX      0.000      0.000        0.759       1.121
## D30           <NA>        CTX      0.000      0.000        0.000       3.776
## HC_21            2        CTX      0.000      0.000        0.000       1.739
## HC_22            2        CTX      0.000      0.000        0.000       1.796
## HC_23            2        CTX      0.000      0.000        0.000       1.730
## HC_24            2        CTX      0.000      0.000        0.000       1.784
## D41           <NA>        CTX      0.000      0.000        0.000       0.511
## D5            <NA>       <NA>      0.000      0.000        0.000       0.650
## D10           <NA>       <NA>      0.000      0.000        0.000       0.751
## D17           <NA>       <NA>      0.000      0.000        0.000       0.629
## D24           <NA>       <NA>      0.000      0.000        0.000       0.639
## HV_1             1        VAN      0.053      0.173        0.000       0.901
## HV_2             1        VAN      0.000      0.228        0.000       0.886
## HV_3             1        VAN      0.069      0.000        0.000       0.844
## HV_4             1        VAN      0.000      0.000        0.000       0.932
## HV_5             1        VAN      0.037      0.284        0.000       1.135
## HV_6             1        VAN      0.046      0.000        0.000       0.819
## D31           <NA>       <NA>      0.000      0.000        0.000       0.453
## D40           <NA>       <NA>      0.000      0.000        0.000       0.456
## D4            <NA>        VAN      0.000      0.000        0.000       0.621
## D11           <NA>        VAN      0.000      1.401        1.045       4.792
## D18           <NA>        VAN      0.000      0.000        1.380      11.605
## D25           <NA>        VAN      0.598      0.000        1.033      10.796
## HV_7             1        VAN      0.000      3.207        4.096       1.351
## HV_8             1        VAN      0.000      9.781        4.362       1.558
## HV_9             1        VAN      0.000      0.000        2.502       1.962
## HV_10            1        VAN      0.000      6.682        4.270       1.587
## HV_11            1        VAN      0.000      0.000        0.000       1.239
## HV_12            1        VAN      0.000      0.000        0.000       1.246
## D32           <NA>        VAN      0.000      0.000        0.000       7.444
## D39           <NA>        VAN      0.000      0.000        0.000       0.418
## C_1           <NA>        VAN      0.000      0.000        0.000       1.525
## D3            <NA>        VAN      0.000      0.000        0.000       0.694
## D12           <NA>        VAN      0.000      1.584        1.135       6.021
## D19           <NA>        VAN      0.489      0.000        1.448       9.526
## C_6           <NA>        VAN      0.550      0.437        0.000       7.267
## D26           <NA>        VAN      0.564      0.000        0.000       9.386
## D33           <NA>        VAN      0.000      0.000        0.000       7.771
## HV_13            1        VAN      1.580      0.767        2.044       1.473
## HV_14            1        VAN      1.238      1.616        1.323       1.637
## HV_15            1        VAN      0.904      0.000        2.563       1.559
## HV_16            1        VAN      1.663      0.000        2.668       1.479
## HV_17            1        VAN      0.000      0.000        0.000       1.233
## D38           <NA>        VAN      0.000      0.000        0.000       3.990
## HV_18            1        VAN      0.000      0.000        0.000       1.251
## D_1           <NA>       <NA>         NA         NA           NA          NA
## D2            <NA>       <NA>      0.000      0.000        0.000       0.734
## D13           <NA>       <NA>      0.000      0.000        0.000       0.726
## D20           <NA>       <NA>      0.000      0.000        0.000       0.748
## D27           <NA>       <NA>      0.000      0.000        0.000       0.729
## HV_19            1        VAN      0.777      0.244        1.325       1.417
## HV_20            1        VAN      0.731      0.000        1.144       1.378
## HV_21            1        VAN      1.220      0.000        0.784       1.435
## HV_22            1        VAN      1.579      0.000        0.750       1.401
## HV_23            1        VAN      0.000      0.000        0.266       1.208
## D34           <NA>       <NA>      0.000      0.000        0.000       1.150
## D37           <NA>       <NA>      0.000      0.000        0.000       0.429
## HV_24            1        VAN      0.000      0.000        0.000       1.187
## C_2           <NA>        VAN         NA         NA           NA          NA
## C_3           <NA>        VAN         NA         NA           NA          NA
## C_4           <NA>        CTX         NA         NA           NA          NA
## C_5           <NA>       <NA>         NA         NA           NA          NA
##       Lactat_mM Formiat_mM Acetat_mM Propionat_mM Isobutyrat_mM Butyrat_mM
## D1        0.000      3.179    49.775       11.560         5.270     22.600
## D14       8.223      5.614    39.530        7.362         0.000      8.315
## D21       4.327      0.000    81.797       10.133         0.000     14.599
## D28       0.000      0.000    86.705       15.767         6.942     27.682
## D35          NA         NA        NA           NA            NA         NA
## D36       0.000      0.000    75.009       11.518         7.914     36.797
## D43          NA         NA        NA           NA            NA         NA
## HC_1      1.338      5.898    68.571       27.257         4.937     48.394
## HC_2      1.287      8.273    71.535       26.884         4.776     45.421
## HC_3      1.310      7.826    76.633       27.355         4.521     44.309
## HC_4      1.238      4.607    70.801       25.104         4.911     50.198
## HC_5      1.321      7.221    73.189       36.735         3.001     33.900
## HC_6      1.423      6.459    52.167       24.281         3.506     53.009
## HC_7      1.613      6.574    75.484       32.194         4.251     41.199
## HC_8      1.750      9.178    71.548       29.822         4.000     43.675
## HC_9      3.095      6.354    69.382       34.611         3.285     36.085
## HC_10     3.283      0.000    64.353       34.868         3.117     39.773
## HC_11     1.293      6.330    70.102       38.008         3.627     35.459
## D7        0.000      0.000    83.366       13.060         8.268     41.333
## D8        0.000      0.000    82.428       12.110         7.901     46.439
## D15       0.000      0.000    77.988        8.582         4.400     49.923
## D22       0.000      0.000    72.326       12.140         0.000     49.287
## D29       8.115      0.000   100.476       13.971         7.021     59.817
## HC_12     1.453      5.867    52.980       25.113         3.912     52.049
## D42       1.035      0.000    69.191       13.707         7.787     36.459
## HC_13     1.552      3.816    70.605       26.932         4.857     39.501
## HC_14     1.505      5.782    68.421       23.025         4.640     42.896
## D6        0.000      0.000    98.070       13.180         8.776     44.333
## D9        0.000      0.000    95.235       12.018         7.794     47.890
## D16       0.000      0.000    90.991        7.423         0.000     44.651
## D23       0.000      0.000    84.906        5.716         0.000     41.240
## HC_15     2.160      5.931    65.980       27.327         3.912     40.000
## HC_16     2.611      0.000    64.822       36.192         3.756     38.253
## HC_17     1.029      5.727    67.573       39.064         3.187     35.241
## HC_18     0.905      5.953    54.688       26.113         3.488     51.621
## HC_19     1.038      3.661    69.252       25.772         5.347     47.850
## HC_20     2.450      0.000    59.247       25.199         5.255     49.825
## D30       0.000      0.000   110.663       12.096         5.428     44.934
## HC_21     1.127      5.850    66.606       24.514         5.216     48.333
## HC_22     1.322      0.000    58.881       33.442         5.939     48.455
## HC_23     1.331     10.916    67.988       35.881         0.000     38.423
## HC_24     0.853      4.163    61.137       27.614         4.315     48.676
## D41       0.000      2.219    57.910       13.311         7.440     35.756
## D5        0.000      0.000    99.568       14.760         7.584     33.721
## D10       0.000      0.000    86.398       11.689         8.204     35.082
## D17       0.000      0.000    85.603       13.331         6.860     29.604
## D24       0.000      0.000    93.635       13.645         6.714     26.799
## HV_1      0.534      5.278    77.645       23.407         3.886     40.871
## HV_2      0.734      6.244    78.499       27.767         4.418     37.598
## HV_3      0.679      7.344    57.376       36.797         4.541     33.022
## HV_4      0.449      7.921    69.627       22.762         3.449     43.182
## HV_5      0.364      3.277    76.749       19.707         4.269     42.214
## HV_6      0.000      8.203    64.866       17.193         2.603     44.874
## D31       0.000      0.000   105.354       14.368         7.536     36.305
## D40       0.000      0.000    83.474       12.933         6.417     23.882
## D4        0.000      0.000   101.441       14.850         7.813     33.872
## D11       0.000      0.000    79.243       17.005         8.207     28.479
## D18       0.000      0.000    48.991       21.902         4.715      6.737
## D25       0.000      9.618    47.637       24.239         0.000      1.522
## HV_7     29.016     10.029    23.392       19.437         2.571     14.661
## HV_8      5.440      2.891    21.269       19.543         3.698     16.553
## HV_9     31.934     10.365    19.888       17.952         3.586     18.726
## HV_10     5.702      2.731    29.523       14.073         2.032     12.535
## HV_11     1.215      5.294    78.115       17.219         4.725     46.529
## HV_12     0.994      9.137    71.280       19.527         2.548     44.810
## D32       8.524      0.000    46.362       19.358         4.941      0.287
## D39       3.808      0.000    41.394       17.022         5.820      2.230
## C_1       0.000      0.000    92.696       11.958         9.009     46.755
## D3        0.000      0.000    94.605       12.469         7.545     42.811
## D12       0.000      0.000    77.564       16.918         7.590     29.690
## D19       0.000      0.000    47.530       23.120         0.000      7.136
## C_6       0.000      0.000    44.220       21.102         0.000      3.984
## D26       0.000      9.711    45.669       24.147         0.000      1.788
## D33       0.000      0.000    41.884       18.412         5.257      1.018
## HV_13     9.844     13.503    11.513        7.533         3.121      9.105
## HV_14     5.033     12.418    20.830       25.674         3.302     12.666
## HV_15    15.308      8.918    11.689       13.911         3.508      7.647
## HV_16     7.130     14.424    12.683        6.561         2.260      5.816
## HV_17     1.269      6.174    73.792       17.652         5.184     48.061
## D38       0.000      0.000    31.082       13.931         6.202      1.974
## HV_18     1.113     10.938    74.408       23.448         2.422     42.201
## D_1          NA         NA        NA           NA            NA         NA
## D2        0.000      0.000    97.498       11.455         7.408     47.083
## D13       0.000      0.000    87.252       11.533         8.115     38.443
## D20       0.000      0.000    87.310       13.946         0.000     31.669
## D27       0.000      0.000    94.750       13.997         7.214     26.634
## HV_19     5.496      8.836    18.665       41.793         4.202     13.490
## HV_20     4.138      8.562    19.052       45.481         4.842     12.833
## HV_21    13.537     10.322    10.928       11.060         2.549      5.173
## HV_22     4.294     12.778    15.865       10.169         2.208      4.905
## HV_23     1.372      6.666    72.238       14.474         4.565     48.449
## D34       0.000      0.000   100.967       15.176         7.531     37.422
## D37       0.000      0.000    85.485       14.118         6.381     23.281
## HV_24     1.171      9.060    64.405       19.489         2.321     42.328
## C_2          NA         NA        NA           NA            NA         NA
## C_3          NA         NA        NA           NA            NA         NA
## C_4          NA         NA        NA           NA            NA         NA
## C_5          NA         NA        NA           NA            NA         NA
##       Isovalerat_mM Valerat_mM Total_SCFA_mM raw_metagenomic_pairs Period
## D1            4.267      3.218       101.838              57121568   <NA>
## D14           1.417      0.365        76.880              44863798     t1
## D21           0.881      0.000       112.555              56304446     t1
## D28           4.870      0.000       142.667              50088120     t1
## D35              NA         NA            NA              55285758     t4
## D36           4.585      6.461       142.830              48021514     t5
## D43              NA         NA            NA              68282190   <NA>
## HC_1          3.842      3.761       166.053                    NA   <NA>
## HC_2          3.394      3.482       168.248                    NA   <NA>
## HC_3          3.614      1.935       169.689                    NA   <NA>
## HC_4          4.062      1.715       164.771                    NA   <NA>
## HC_5          2.809      1.264       161.490                    NA   <NA>
## HC_6          3.828      1.687       148.419                    NA   <NA>
## HC_7          3.894      4.031       171.112                    NA     t1
## HC_8          2.980      3.486       168.639                    NA     t1
## HC_9          3.394      1.704       160.404                    NA     t1
## HC_10         3.610      1.147       152.718                    NA     t1
## HC_11         3.258      1.282       161.474                    NA     t1
## D7            6.524      7.962       161.603              52833332   <NA>
## D8            6.534      4.832       161.445              53009028     t1
## D15           3.876      1.147       148.695              47492202     t1
## D22           7.405      0.522       144.749              55638128     t1
## D29           8.804      9.447       208.178              68334190     t4
## HC_12         3.952      1.797       148.977                    NA     t1
## D42           7.614      8.207       144.455              47906582     t5
## HC_13         4.427      5.098       158.525                    NA     t1
## HC_14         3.259      4.630       155.243                    NA     t1
## D6            6.422      7.313       178.886              51668602   <NA>
## D9            6.038      3.410       173.501              47401086     t1
## D16           2.036      0.036       148.110              44839890     t1
## D23           1.044      0.000       136.326              56766854     t1
## HC_15         5.518      1.598       153.917                    NA     t1
## HC_16         5.307      1.038       153.427                    NA     t1
## HC_17         3.844      1.038       157.813                    NA     t1
## HC_18         4.236      1.783       150.045                    NA     t1
## HC_19         5.290      6.032       165.665                    NA     t1
## HC_20         6.343      4.494       154.693                    NA     t1
## D30           6.609      0.000       183.506              49315184     t4
## HC_21         4.836      6.599       164.820                    NA     t1
## HC_22         7.871      1.967       159.673                    NA     t1
## HC_23         4.747      1.182       162.198                    NA     t1
## HC_24         5.126      1.616       155.284                    NA     t1
## D41           6.925      0.233       124.305              59842344     t5
## D5            6.467      8.363       171.113              55615446   <NA>
## D10           6.429      6.960       155.513              46820988     t1
## D17           6.109      0.000       142.136              49791554     t1
## D24           5.297      0.000       146.729              48166738     t1
## HV_1          4.754      1.239       158.741                    NA   <NA>
## HV_2          5.287      1.303       162.964                    NA   <NA>
## HV_3          5.700      0.000       146.372                    NA   <NA>
## HV_4          2.876      0.475       151.673                    NA   <NA>
## HV_5          4.334      0.284       152.654                    NA   <NA>
## HV_6          3.593      0.390       142.587                    NA   <NA>
## D31           6.443      7.805       178.264              56873318     t4
## D40           4.188      6.227       137.577              61210784     t5
## D4            6.542      8.752       173.891              50262020   <NA>
## D11           6.182      4.307       148.215              55159700     t1
## D18           5.970      0.436       100.356              42389522     t1
## D25           5.250      0.000        99.062              49858296     t1
## HV_7          4.272      1.059       113.091                    NA     t1
## HV_8          6.617      1.054        92.766                    NA     t1
## HV_9          5.648      1.069       113.632                    NA     t1
## HV_10         2.750      1.147        83.032                    NA     t1
## HV_11         4.203      1.259       159.798                    NA     t1
## HV_12         3.415      1.293       154.250                    NA     t1
## D32           7.161      0.000        85.553              50291070     t4
## D39           4.591      0.000        71.475              64514742     t5
## C_1           6.637      7.112       175.692                    NA   pret
## D3            6.110      7.167       171.401              54057756   <NA>
## D12           5.858      3.135       146.776              45384782     t1
## D19           5.905      0.337        93.554              54226192     t1
## C_6           5.818      0.000        82.391                    NA     t1
## D26           5.231      0.000        86.221              50932862     t1
## D33           7.101      0.000        81.443              59788084     t4
## HV_13         5.593      1.046        67.122                    NA     t1
## HV_14         5.969      1.079        92.785                    NA     t1
## HV_15         5.601      1.023        72.631                    NA     t1
## HV_16         4.042      0.958        59.684                    NA     t1
## HV_17         4.140      1.393       158.898                    NA     t1
## D38           4.674      0.000        61.853              38395332     t5
## HV_18         2.336      1.288       159.405                    NA     t1
## D_1              NA         NA            NA                    NA   <NA>
## D2            5.748      7.114       177.040              60148698   <NA>
## D13           6.024      6.447       158.540              52085624     t1
## D20           6.040      0.000       139.713              46655108     t1
## D27           5.309      0.000       148.633              70694598     t1
## HV_19         6.903      1.239       104.387                    NA     t1
## HV_20         7.904      1.243       107.308                    NA     t1
## HV_21         3.751      0.973        61.732                    NA     t1
## HV_22         4.026      0.971        58.946                    NA     t1
## HV_23         3.707      1.426       154.371                    NA     t1
## D34           7.673      7.628       177.547              56546096     t4
## D37           5.541      6.703       141.938              51876272     t5
## HV_24         2.300      1.239       143.500                    NA     t1
## C_2              NA         NA            NA                    NA     t5
## C_3              NA         NA            NA                    NA     t5
## C_4              NA         NA            NA                    NA     t5
## C_5              NA         NA            NA                    NA     t5
##       Reactor_Treatment_Dose   Treatment_Dose
## D1              CR_UNTREATED        UNTREATED
## D14             CR_UNTREATED        UNTREATED
## D21             CR_UNTREATED        UNTREATED
## D28             CR_UNTREATED        UNTREATED
## D35             CR_UNTREATED        UNTREATED
## D36             CR_UNTREATED        UNTREATED
## D43                    DONOR            DONOR
## HC_1       TR1_CTX+HV292.120    CTX+HV292.120
## HC_2               TR2_CTX20            CTX20
## HC_3      TR3_CTX+HV292.1200   CTX+HV292.1200
## HC_4              TR4_CTX200           CTX200
## HC_5             TR5_HV292.1          HV292.1
## HC_6            CR_UNTREATED        UNTREATED
## HC_7       TR1_CTX+HV292.120    CTX+HV292.120
## HC_8               TR2_CTX20            CTX20
## HC_9      TR3_CTX+HV292.1200   CTX+HV292.1200
## HC_10             TR4_CTX200           CTX200
## HC_11            TR5_HV292.1          HV292.1
## D7         TR1_CTX+HV292.120    CTX+HV292.120
## D8         TR1_CTX+HV292.120    CTX+HV292.120
## D15        TR1_CTX+HV292.120    CTX+HV292.120
## D22        TR1_CTX+HV292.120    CTX+HV292.120
## D29        TR1_CTX+HV292.120    CTX+HV292.120
## HC_12           CR_UNTREATED        UNTREATED
## D42        TR1_CTX+HV292.120    CTX+HV292.120
## HC_13      TR1_CTX+HV292.120    CTX+HV292.120
## HC_14              TR2_CTX20            CTX20
## D6                 TR2_CTX20            CTX20
## D9                 TR2_CTX20            CTX20
## D16                TR2_CTX20            CTX20
## D23                TR2_CTX20            CTX20
## HC_15     TR3_CTX+HV292.1200   CTX+HV292.1200
## HC_16             TR4_CTX200           CTX200
## HC_17            TR5_HV292.1          HV292.1
## HC_18           CR_UNTREATED        UNTREATED
## HC_19      TR1_CTX+HV292.120    CTX+HV292.120
## HC_20              TR2_CTX20            CTX20
## D30                TR2_CTX20            CTX20
## HC_21     TR3_CTX+HV292.1200   CTX+HV292.1200
## HC_22             TR4_CTX200           CTX200
## HC_23            TR5_HV292.1          HV292.1
## HC_24           CR_UNTREATED        UNTREATED
## D41                TR2_CTX20            CTX20
## D5               TR3_HV292.1          HV292.1
## D10              TR3_HV292.1          HV292.1
## D17              TR3_HV292.1          HV292.1
## D24              TR3_HV292.1          HV292.1
## HV_1     TR1_VAN+CCUG5916890  VAN+CCUG5916890
## HV_2               TR2_VAN90            VAN90
## HV_3    TR3_VAN+CCUG59168600 VAN+CCUG59168600
## HV_4              TR4_VAN600           VAN600
## HV_5           TR5_CCUG59168        CCUG59168
## HV_6            CR_UNTREATED        UNTREATED
## D31              TR3_HV292.1          HV292.1
## D40              TR3_HV292.1          HV292.1
## D4                 TR4_VAN90            VAN90
## D11                TR4_VAN90            VAN90
## D18                TR4_VAN90            VAN90
## D25                TR4_VAN90            VAN90
## HV_7     TR1_VAN+CCUG5916890  VAN+CCUG5916890
## HV_8               TR2_VAN90            VAN90
## HV_9    TR3_VAN+CCUG59168600 VAN+CCUG59168600
## HV_10             TR4_VAN600           VAN600
## HV_11          TR5_CCUG59168        CCUG59168
## HV_12           CR_UNTREATED        UNTREATED
## D32                TR4_VAN90            VAN90
## D39                TR4_VAN90            VAN90
## C_1      TR5_VAN+CCUG5916890  VAN+CCUG5916890
## D3       TR5_VAN+CCUG5916890  VAN+CCUG5916890
## D12      TR5_VAN+CCUG5916890  VAN+CCUG5916890
## D19      TR5_VAN+CCUG5916890  VAN+CCUG5916890
## C_6      TR5_VAN+CCUG5916890  VAN+CCUG5916890
## D26      TR5_VAN+CCUG5916890  VAN+CCUG5916890
## D33      TR5_VAN+CCUG5916890  VAN+CCUG5916890
## HV_13    TR1_VAN+CCUG5916890  VAN+CCUG5916890
## HV_14              TR2_VAN90            VAN90
## HV_15   TR3_VAN+CCUG59168600 VAN+CCUG59168600
## HV_16             TR4_VAN600           VAN600
## HV_17          TR5_CCUG59168        CCUG59168
## D38      TR5_VAN+CCUG5916890  VAN+CCUG5916890
## HV_18           CR_UNTREATED        UNTREATED
## D_1                    DONOR            DONOR
## D2             TR6_CCUG59168        CCUG59168
## D13            TR6_CCUG59168        CCUG59168
## D20            TR6_CCUG59168        CCUG59168
## D27            TR6_CCUG59168        CCUG59168
## HV_19    TR1_VAN+CCUG5916890  VAN+CCUG5916890
## HV_20              TR2_VAN90            VAN90
## HV_21   TR3_VAN+CCUG59168600 VAN+CCUG59168600
## HV_22             TR4_VAN600           VAN600
## HV_23          TR5_CCUG59168        CCUG59168
## D34            TR6_CCUG59168        CCUG59168
## D37            TR6_CCUG59168        CCUG59168
## HV_24           CR_UNTREATED        UNTREATED
## C_2      TR5_VAN+CCUG5916890  VAN+CCUG5916890
## C_3                TR4_VAN90            VAN90
## C_4              TR3_HV292.1          HV292.1
## C_5                TR2_CTX20            CTX20
```




```r
"data/processed/chicken1_full_gene_catalog_full_metrics_table.tsv.gz"%>% 
  here::here() %>% 
  read_tsv() -> chicken1
```

```
## Warning: One or more parsing issues, see `problems()` for details
```

```
## Rows: 1634985 Columns: 394
## ── Column specification ────────────────────────────────────────────────────────
## Delimiter: "\t"
## chr  (48): ORF_ID, Contig_ID, SQM_Gene.name, SQM_Tax_orf, SQM_KEGG.ID, SQM_K...
## dbl (317): SQM_Length.NT, SQM_Length.AA, SQM_GC.perc, plasX_score, PathoFact...
## lgl  (29): rgi_CARD_Megahit, rgi_CARD_Cut_Off, rgi_CARD_Pass_Bitscore, rgi_C...
## 
## ℹ Use `spec()` to retrieve the full column specification for this data.
## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.
```



```r
"data/processed/human1_full_gene_catalog_full_metrics_table.tsv.gz" %>% 
  here::here() %>% 
  read_tsv() -> human1
```

```
## Warning: One or more parsing issues, see `problems()` for details
```

```
## Rows: 791180 Columns: 392
## ── Column specification ────────────────────────────────────────────────────────
## Delimiter: "\t"
## chr  (55): ORF_ID, Contig_ID, SQM_Gene.name, SQM_Tax_orf, SQM_KEGG.ID, SQM_K...
## dbl (322): SQM_Length.NT, SQM_Length.AA, SQM_GC.perc, rgi_CARD_Pass_Bitscore...
## lgl  (15): rgi_CARD_SNPs_in_Best_Hit_ARO, rgi_CARD_Other_SNPs, rgi_CARD_Nudg...
## 
## ℹ Use `spec()` to retrieve the full column specification for this data.
## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.
```


```r
chicken1 %>% 
  head() -> toy_chicken1

human1 %>% 
  head() -> toy_human1
```


```r
toy_chicken1 %>% 
  colnames() %>% 
  head()
```

```
## [1] "ORF_ID"        "Contig_ID"     "SQM_Length.NT" "SQM_Length.AA"
## [5] "SQM_GC.perc"   "SQM_Gene.name"
```

```r
# contains plasX_score
```



```r
toy_human1 %>% 
  colnames() %>% 
  head()
```

```
## [1] "ORF_ID"        "Contig_ID"     "SQM_Length.NT" "SQM_Length.AA"
## [5] "SQM_GC.perc"   "SQM_Gene.name"
```



```r
human1 %>% 
  colnames() %>% 
  intersect(chicken1 %>% colnames()) -> commun_cols
```


```r
human1 %>% 
  select(commun_cols) %>% 
  rbind(chicken1 %>% 
          select(all_of(commun_cols))) -> common_col_bind
```

```
## Note: Using an external vector in selections is ambiguous.
## ℹ Use `all_of(commun_cols)` instead of `commun_cols` to silence this message.
## ℹ See <https://tidyselect.r-lib.org/reference/faq-external-vector.html>.
## This message is displayed once per session.
```


Extract quantitative data:


```r
human1 %>% 
  select(ORF_ID, contains(c("TPM.","Coverage.", "Num_Gi.", "Num_Gi_pc.", "Raw.read.count.", "Raw.base.count"))) -> human1_quant
# select(ORF_ID, contains(c("Num_Gi_pc."))) -> human1_quant
```


```r
chicken1 %>% 
  select(ORF_ID, contains(c("TPM.","Coverage.", "Num_Gi.", "Num_Gi_pc.", "Raw.read.count.", "Raw.base.count"))) -> chicken1_quant
# select(ORF_ID, contains(c("Num_Gi_pc."))) -> chicken1_quant
```



```r
common_col_bind %>% 
  left_join(chicken1_quant,
            by = c("ORF_ID" = "ORF_ID")) %>% 
  left_join(human1_quant,
            by = c("ORF_ID" = "ORF_ID")) %>% 
  mutate_at(vars(contains(c("TPM.","Coverage.", "Num_Gi.", "Num_Gi_pc.", "Raw.read.count.", "Raw.base.count"))), ~replace_na(., 0)) -> combined_full
```


```r
combined_full %>% 
  mutate_if(is.character, as.factor) -> combined_full_fct
```


```r
combined_full_fct %>% 
  glimpse()
```

```
## Rows: 2,426,165
## Columns: 682
## $ ORF_ID                                   <fct> h1_megahit_1_1-312, h1_megahi…
## $ Contig_ID                                <fct> h1_megahit_1, h1_megahit_2, h…
## $ SQM_Length.NT                            <dbl> 312, 174, 132, 114, 138, 72, …
## $ SQM_Length.AA                            <dbl> 104, 58, 44, 38, 46, 24, 79, …
## $ SQM_GC.perc                              <dbl> 34.62, 38.51, 40.91, 51.75, 4…
## $ SQM_Gene.name                            <fct> "RP-S5, MRPS5, rpsE", NA, NA,…
## $ SQM_Tax_orf                              <fct> k_Bacteria;n_Terrabacteria gr…
## $ SQM_KEGG.ID                              <fct> K02988*, NA, NA, NA, NA, NA, …
## $ SQM_KEGGFUN                              <fct> "small subunit ribosomal prot…
## $ SQM_KEGGPATH                             <fct> "Genetic Information Processi…
## $ SQM_COG.ID                               <fct> COG0098*, NA, NA, ENOG410ZP7K…
## $ SQM_COGFUN                               <fct> "Ribosomal protein S5", NA, N…
## $ SQM_COGPATH                              <fct> "Translation, ribosomal struc…
## $ SQM_PFAM                                 <fct> "PF03719 [Ribosomal protein S…
## $ rgi_CARD_Cut_Off                         <fct> NA, NA, NA, NA, NA, NA, NA, N…
## $ rgi_CARD_Pass_Bitscore                   <dbl> NA, NA, NA, NA, NA, NA, NA, N…
## $ rgi_CARD_Best_Hit_Bitscore               <dbl> NA, NA, NA, NA, NA, NA, NA, N…
## $ rgi_CARD_Best_Hit_ARO                    <fct> NA, NA, NA, NA, NA, NA, NA, N…
## $ rgi_CARD_Best_Identities                 <dbl> NA, NA, NA, NA, NA, NA, NA, N…
## $ rgi_CARD_ARO                             <dbl> NA, NA, NA, NA, NA, NA, NA, N…
## $ rgi_CARD_Model_type                      <fct> NA, NA, NA, NA, NA, NA, NA, N…
## $ rgi_CARD_SNPs_in_Best_Hit_ARO            <lgl> NA, NA, NA, NA, NA, NA, NA, N…
## $ rgi_CARD_Other_SNPs                      <lgl> NA, NA, NA, NA, NA, NA, NA, N…
## $ rgi_CARD_ID                              <fct> NA, NA, NA, NA, NA, NA, NA, N…
## $ rgi_CARD_Model_ID                        <dbl> NA, NA, NA, NA, NA, NA, NA, N…
## $ rgi_CARD_Nudged                          <lgl> NA, NA, NA, NA, NA, NA, NA, N…
## $ rgi_CARD_Note                            <lgl> NA, NA, NA, NA, NA, NA, NA, N…
## $ PathoFact_AMR_ARG                        <fct> NA, NA, NA, NA, NA, NA, NA, N…
## $ PathoFact_AMR_ARG_SNPs                   <lgl> NA, NA, NA, NA, NA, NA, NA, N…
## $ PathoFact_AMR_AMR_category               <fct> NA, NA, NA, NA, NA, NA, NA, N…
## $ PathoFact_AMR_AMR_sub_class              <fct> NA, NA, NA, NA, NA, NA, NA, N…
## $ PathoFact_AMR_Resistance_mechanism       <fct> NA, NA, NA, NA, NA, NA, NA, N…
## $ PathoFact_AMR_Database                   <fct> NA, NA, NA, NA, NA, NA, NA, N…
## $ PathoFact_AMR_MGE_prediction             <fct> NA, NA, NA, NA, NA, NA, NA, N…
## $ PathoFact_Tox_Toxin_classification       <fct> NA, NA, NA, NA, NA, NA, NA, N…
## $ PathoFact_Tox_HMM_Name                   <fct> NA, NA, NA, NA, NA, NA, NA, N…
## $ PathoFact_Tox_Score                      <dbl> NA, NA, NA, NA, NA, NA, NA, N…
## $ PathoFact_Tox_Significance_evalue        <dbl> NA, NA, NA, NA, NA, NA, NA, N…
## $ PathoFact_Tox_NAME                       <fct> NA, NA, NA, NA, NA, NA, NA, N…
## $ PathoFact_Tox_Description                <fct> NA, NA, NA, NA, NA, NA, NA, N…
## $ viralverify_Prediction                   <fct> Uncertain - too short, Uncert…
## $ isescan_family                           <fct> NA, NA, NA, NA, NA, NA, NA, N…
## $ isescan_cluster                          <fct> NA, NA, NA, NA, NA, NA, NA, N…
## $ isescan_type                             <fct> NA, NA, NA, NA, NA, NA, NA, N…
## $ isescan_count_IS_contig                  <dbl> NA, NA, NA, NA, NA, NA, NA, N…
## $ SQM_Tax_contig                           <fct> k_Bacteria;n_Terrabacteria gr…
## $ SQM_Disparity_contig                     <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,…
## $ SQM_GC_perc_contig                       <dbl> 34.50, 39.67, 39.67, 51.13, 5…
## $ SQM_Length_contig                        <dbl> 313, 305, 305, 221, 221, 324,…
## $ SQM_Num_genes_contig                     <dbl> 1, 2, 2, 2, 2, 2, 2, 1, 1, 1,…
## $ anvio_KEGG                               <fct> NA, NA, NA, NA, NA, NA, NA, N…
## $ `anvio_KEGG (ACCESSION)`                 <fct> NA, NA, NA, NA, NA, NA, NA, N…
## $ anvio_COG                                <fct> NA, NA, NA, NA, NA, NA, NA, N…
## $ `anvio_COG (ACCESSION)`                  <fct> NA, NA, NA, NA, NA, NA, NA, N…
## $ anvio_COGPATH                            <fct> NA, NA, NA, NA, NA, NA, NA, N…
## $ `anvio_COGPATH (ACCESSION)`              <fct> NA, NA, NA, NA, NA, NA, NA, N…
## $ anvio_KEGGPATH                           <fct> NA, NA, NA, NA, NA, NA, NA, N…
## $ `anvio_KEGGPATH (ACCESSION)`             <fct> NA, NA, NA, NA, NA, NA, NA, N…
## $ anvio_PFAM                               <fct> NA, NA, NA, NA, NA, NA, NA, N…
## $ `anvio_PFAM (ACCESSION)`                 <fct> NA, NA, NA, NA, NA, NA, NA, N…
## $ anvio_KOfam                              <fct> NA, NA, NA, NA, NA, NA, NA, N…
## $ `anvio_KOfam (ACCESSION)`                <fct> NA, NA, NA, NA, NA, NA, NA, N…
## $ anvio_KEGG_Module                        <fct> NA, NA, NA, NA, NA, NA, NA, N…
## $ `anvio_KEGG_Module (ACCESSION)`          <fct> NA, NA, NA, NA, NA, NA, NA, N…
## $ anvio_KEGG_Class                         <fct> NA, NA, NA, NA, NA, NA, NA, N…
## $ `anvio_KEGG_Class (ACCESSION)`           <fct> NA, NA, NA, NA, NA, NA, NA, N…
## $ CONCOCT_bin                              <fct> NA, NA, NA, NA, NA, NA, NA, N…
## $ DAS_bin                                  <fct> NA, NA, NA, NA, NA, NA, NA, N…
## $ GTDB_tax_CONCOCT                         <fct> NA, NA, NA, NA, NA, NA, NA, N…
## $ GTDB_tax_DAS                             <fct> NA, NA, NA, NA, NA, NA, NA, N…
## $ SQM_DAS_Tax                              <lgl> NA, NA, NA, NA, NA, NA, NA, N…
## $ `SQM_DAS_Tax 16S`                        <lgl> NA, NA, NA, NA, NA, NA, NA, N…
## $ SQM_DAS_Length                           <lgl> NA, NA, NA, NA, NA, NA, NA, N…
## $ `SQM_DAS_GC perc`                        <lgl> NA, NA, NA, NA, NA, NA, NA, N…
## $ `SQM_DAS_Num contigs`                    <lgl> NA, NA, NA, NA, NA, NA, NA, N…
## $ SQM_DAS_Disparity                        <lgl> NA, NA, NA, NA, NA, NA, NA, N…
## $ SQM_DAS_Completeness                     <lgl> NA, NA, NA, NA, NA, NA, NA, N…
## $ SQM_DAS_Contamination                    <lgl> NA, NA, NA, NA, NA, NA, NA, N…
## $ `SQM_DAS_Strain het`                     <lgl> NA, NA, NA, NA, NA, NA, NA, N…
## $ SQM_DAS_HQ                               <lgl> NA, NA, NA, NA, NA, NA, NA, N…
## $ ANVIO_SUMMARY_DAS_total_length           <dbl> NA, NA, NA, NA, NA, NA, NA, N…
## $ ANVIO_SUMMARY_DAS_num_contigs            <dbl> NA, NA, NA, NA, NA, NA, NA, N…
## $ ANVIO_SUMMARY_DAS_N50                    <dbl> NA, NA, NA, NA, NA, NA, NA, N…
## $ ANVIO_SUMMARY_DAS_GC_content             <dbl> NA, NA, NA, NA, NA, NA, NA, N…
## $ ANVIO_SUMMARY_DAS_percent_completion     <dbl> NA, NA, NA, NA, NA, NA, NA, N…
## $ ANVIO_SUMMARY_DAS_percent_redundancy     <dbl> NA, NA, NA, NA, NA, NA, NA, N…
## $ ANVIO_DAS_HQ                             <fct> NA, NA, NA, NA, NA, NA, NA, N…
## $ ANVIO_SUMMARY_CONCOCT_total_length       <dbl> NA, NA, NA, NA, NA, NA, NA, N…
## $ ANVIO_SUMMARY_CONCOCT_num_contigs        <dbl> NA, NA, NA, NA, NA, NA, NA, N…
## $ ANVIO_SUMMARY_CONCOCT_N50                <dbl> NA, NA, NA, NA, NA, NA, NA, N…
## $ ANVIO_SUMMARY_CONCOCT_GC_content         <dbl> NA, NA, NA, NA, NA, NA, NA, N…
## $ ANVIO_SUMMARY_CONCOCT_percent_completion <dbl> NA, NA, NA, NA, NA, NA, NA, N…
## $ ANVIO_SUMMARY_CONCOCT_percent_redundancy <dbl> NA, NA, NA, NA, NA, NA, NA, N…
## $ ANVIO_CONCOCT_HQ                         <fct> NA, NA, NA, NA, NA, NA, NA, N…
## $ ORF_TPM.C1                               <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,…
## $ ORF_TPM.C2                               <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,…
## $ ORF_TPM.C3                               <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,…
## $ ORF_TPM.C4                               <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,…
## $ ORF_TPM.C5                               <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,…
## $ ORF_TPM.C6                               <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,…
## $ ORF_TPM.D10                              <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,…
## $ ORF_TPM.D11                              <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,…
## $ ORF_TPM.D12                              <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,…
## $ ORF_TPM.D13                              <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,…
## $ ORF_TPM.D14                              <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,…
## $ ORF_TPM.D15                              <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,…
## $ ORF_TPM.D16                              <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,…
## $ ORF_TPM.D17                              <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,…
## $ ORF_TPM.D18                              <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,…
## $ ORF_TPM.D19                              <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,…
## $ ORF_TPM.D1                               <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,…
## $ ORF_TPM.D20                              <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,…
## $ ORF_TPM.D21                              <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,…
## $ ORF_TPM.D22                              <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,…
## $ ORF_TPM.D23                              <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,…
## $ ORF_TPM.D24                              <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,…
## $ ORF_TPM.D25                              <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,…
## $ ORF_TPM.D26                              <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,…
## $ ORF_TPM.D27                              <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,…
## $ ORF_TPM.D28                              <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,…
## $ ORF_TPM.D29                              <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,…
## $ ORF_TPM.D2                               <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,…
## $ ORF_TPM.D30                              <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,…
## $ ORF_TPM.D31                              <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,…
## $ ORF_TPM.D32                              <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,…
## $ ORF_TPM.D33                              <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,…
## $ ORF_TPM.D34                              <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,…
## $ ORF_TPM.D35                              <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,…
## $ ORF_TPM.D36                              <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,…
## $ ORF_TPM.D37                              <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,…
## $ ORF_TPM.D38                              <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,…
## $ ORF_TPM.D39                              <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,…
## $ ORF_TPM.D3                               <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,…
## $ ORF_TPM.D40                              <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,…
## $ ORF_TPM.D41                              <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,…
## $ ORF_TPM.D42                              <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,…
## $ ORF_TPM.D43                              <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,…
## $ ORF_TPM.D4                               <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,…
## $ ORF_TPM.D5                               <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,…
## $ ORF_TPM.D6                               <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,…
## $ ORF_TPM.D7                               <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,…
## $ ORF_TPM.D8                               <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,…
## $ ORF_TPM.D9                               <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,…
## $ ORF_Coverage.C1                          <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,…
## $ ORF_Coverage.C2                          <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,…
## $ ORF_Coverage.C3                          <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,…
## $ ORF_Coverage.C4                          <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,…
## $ ORF_Coverage.C5                          <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,…
## $ ORF_Coverage.C6                          <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,…
## $ ORF_Coverage.D10                         <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,…
## $ ORF_Coverage.D11                         <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,…
## $ ORF_Coverage.D12                         <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,…
## $ ORF_Coverage.D13                         <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,…
## $ ORF_Coverage.D14                         <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,…
## $ ORF_Coverage.D15                         <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,…
## $ ORF_Coverage.D16                         <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,…
## $ ORF_Coverage.D17                         <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,…
## $ ORF_Coverage.D18                         <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,…
## $ ORF_Coverage.D19                         <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,…
## $ ORF_Coverage.D1                          <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,…
## $ ORF_Coverage.D20                         <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,…
## $ ORF_Coverage.D21                         <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,…
## $ ORF_Coverage.D22                         <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,…
## $ ORF_Coverage.D23                         <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,…
## $ ORF_Coverage.D24                         <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,…
## $ ORF_Coverage.D25                         <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,…
## $ ORF_Coverage.D26                         <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,…
## $ ORF_Coverage.D27                         <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,…
## $ ORF_Coverage.D28                         <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,…
## $ ORF_Coverage.D29                         <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,…
## $ ORF_Coverage.D2                          <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,…
## $ ORF_Coverage.D30                         <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,…
## $ ORF_Coverage.D31                         <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,…
## $ ORF_Coverage.D32                         <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,…
## $ ORF_Coverage.D33                         <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,…
## $ ORF_Coverage.D34                         <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,…
## $ ORF_Coverage.D35                         <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,…
## $ ORF_Coverage.D36                         <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,…
## $ ORF_Coverage.D37                         <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,…
## $ ORF_Coverage.D38                         <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,…
## $ ORF_Coverage.D39                         <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,…
## $ ORF_Coverage.D3                          <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,…
## $ ORF_Coverage.D40                         <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,…
## $ ORF_Coverage.D41                         <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,…
## $ ORF_Coverage.D42                         <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,…
## $ ORF_Coverage.D43                         <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,…
## $ ORF_Coverage.D4                          <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,…
## $ ORF_Coverage.D5                          <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,…
## $ ORF_Coverage.D6                          <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,…
## $ ORF_Coverage.D7                          <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,…
## $ ORF_Coverage.D8                          <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,…
## $ ORF_Coverage.D9                          <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,…
## $ ORF_Num_Gi.C1                            <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,…
## $ ORF_Num_Gi.C2                            <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,…
## $ ORF_Num_Gi.C3                            <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,…
## $ ORF_Num_Gi.C4                            <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,…
## $ ORF_Num_Gi.C5                            <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,…
## $ ORF_Num_Gi.C6                            <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,…
## $ ORF_Num_Gi.D10                           <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,…
## $ ORF_Num_Gi.D11                           <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,…
## $ ORF_Num_Gi.D12                           <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,…
## $ ORF_Num_Gi.D13                           <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,…
## $ ORF_Num_Gi.D14                           <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,…
## $ ORF_Num_Gi.D15                           <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,…
## $ ORF_Num_Gi.D16                           <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,…
## $ ORF_Num_Gi.D17                           <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,…
## $ ORF_Num_Gi.D18                           <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,…
## $ ORF_Num_Gi.D19                           <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,…
## $ ORF_Num_Gi.D1                            <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,…
## $ ORF_Num_Gi.D20                           <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,…
## $ ORF_Num_Gi.D21                           <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,…
## $ ORF_Num_Gi.D22                           <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,…
## $ ORF_Num_Gi.D23                           <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,…
## $ ORF_Num_Gi.D24                           <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,…
## $ ORF_Num_Gi.D25                           <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,…
## $ ORF_Num_Gi.D26                           <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,…
## $ ORF_Num_Gi.D27                           <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,…
## $ ORF_Num_Gi.D28                           <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,…
## $ ORF_Num_Gi.D29                           <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,…
## $ ORF_Num_Gi.D2                            <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,…
## $ ORF_Num_Gi.D30                           <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,…
## $ ORF_Num_Gi.D31                           <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,…
## $ ORF_Num_Gi.D32                           <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,…
## $ ORF_Num_Gi.D33                           <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,…
## $ ORF_Num_Gi.D34                           <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,…
## $ ORF_Num_Gi.D35                           <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,…
## $ ORF_Num_Gi.D36                           <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,…
## $ ORF_Num_Gi.D37                           <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,…
## $ ORF_Num_Gi.D38                           <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,…
## $ ORF_Num_Gi.D39                           <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,…
## $ ORF_Num_Gi.D3                            <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,…
## $ ORF_Num_Gi.D40                           <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,…
## $ ORF_Num_Gi.D41                           <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,…
## $ ORF_Num_Gi.D42                           <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,…
## $ ORF_Num_Gi.D43                           <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,…
## $ ORF_Num_Gi.D4                            <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,…
## $ ORF_Num_Gi.D5                            <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,…
## $ ORF_Num_Gi.D6                            <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,…
## $ ORF_Num_Gi.D7                            <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,…
## $ ORF_Num_Gi.D8                            <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,…
## $ ORF_Num_Gi.D9                            <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,…
## $ ORF_Num_Gi_pc.C1                         <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,…
## $ ORF_Num_Gi_pc.C2                         <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,…
## $ ORF_Num_Gi_pc.C3                         <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,…
## $ ORF_Num_Gi_pc.C4                         <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,…
## $ ORF_Num_Gi_pc.C5                         <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,…
## $ ORF_Num_Gi_pc.C6                         <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,…
## $ ORF_Num_Gi_pc.D10                        <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,…
## $ ORF_Num_Gi_pc.D11                        <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,…
## $ ORF_Num_Gi_pc.D12                        <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,…
## $ ORF_Num_Gi_pc.D13                        <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,…
## $ ORF_Num_Gi_pc.D14                        <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,…
## $ ORF_Num_Gi_pc.D15                        <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,…
## $ ORF_Num_Gi_pc.D16                        <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,…
## $ ORF_Num_Gi_pc.D17                        <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,…
## $ ORF_Num_Gi_pc.D18                        <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,…
## $ ORF_Num_Gi_pc.D19                        <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,…
## $ ORF_Num_Gi_pc.D1                         <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,…
## $ ORF_Num_Gi_pc.D20                        <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,…
## $ ORF_Num_Gi_pc.D21                        <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,…
## $ ORF_Num_Gi_pc.D22                        <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,…
## $ ORF_Num_Gi_pc.D23                        <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,…
## $ ORF_Num_Gi_pc.D24                        <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,…
## $ ORF_Num_Gi_pc.D25                        <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,…
## $ ORF_Num_Gi_pc.D26                        <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,…
## $ ORF_Num_Gi_pc.D27                        <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,…
## $ ORF_Num_Gi_pc.D28                        <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,…
## $ ORF_Num_Gi_pc.D29                        <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,…
## $ ORF_Num_Gi_pc.D2                         <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,…
## $ ORF_Num_Gi_pc.D30                        <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,…
## $ ORF_Num_Gi_pc.D31                        <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,…
## $ ORF_Num_Gi_pc.D32                        <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,…
## $ ORF_Num_Gi_pc.D33                        <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,…
## $ ORF_Num_Gi_pc.D34                        <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,…
## $ ORF_Num_Gi_pc.D35                        <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,…
## $ ORF_Num_Gi_pc.D36                        <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,…
## $ ORF_Num_Gi_pc.D37                        <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,…
## $ ORF_Num_Gi_pc.D38                        <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,…
## $ ORF_Num_Gi_pc.D39                        <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,…
## $ ORF_Num_Gi_pc.D3                         <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,…
## $ ORF_Num_Gi_pc.D40                        <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,…
## $ ORF_Num_Gi_pc.D41                        <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,…
## $ ORF_Num_Gi_pc.D42                        <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,…
## $ ORF_Num_Gi_pc.D43                        <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,…
## $ ORF_Num_Gi_pc.D4                         <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,…
## $ ORF_Num_Gi_pc.D5                         <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,…
## $ ORF_Num_Gi_pc.D6                         <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,…
## $ ORF_Num_Gi_pc.D7                         <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,…
## $ ORF_Num_Gi_pc.D8                         <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,…
## $ ORF_Num_Gi_pc.D9                         <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,…
## $ ORF_Raw.read.count.C1                    <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,…
## $ ORF_Raw.read.count.C2                    <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,…
## $ ORF_Raw.read.count.C3                    <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,…
## $ ORF_Raw.read.count.C4                    <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,…
## $ ORF_Raw.read.count.C5                    <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,…
## $ ORF_Raw.read.count.C6                    <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,…
## $ ORF_Raw.read.count.D10                   <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,…
## $ ORF_Raw.read.count.D11                   <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,…
## $ ORF_Raw.read.count.D12                   <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,…
## $ ORF_Raw.read.count.D13                   <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,…
## $ ORF_Raw.read.count.D14                   <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,…
## $ ORF_Raw.read.count.D15                   <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,…
## $ ORF_Raw.read.count.D16                   <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,…
## $ ORF_Raw.read.count.D17                   <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,…
## $ ORF_Raw.read.count.D18                   <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,…
## $ ORF_Raw.read.count.D19                   <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,…
## $ ORF_Raw.read.count.D1                    <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,…
## $ ORF_Raw.read.count.D20                   <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,…
## $ ORF_Raw.read.count.D21                   <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,…
## $ ORF_Raw.read.count.D22                   <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,…
## $ ORF_Raw.read.count.D23                   <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,…
## $ ORF_Raw.read.count.D24                   <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,…
## $ ORF_Raw.read.count.D25                   <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,…
## $ ORF_Raw.read.count.D26                   <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,…
## $ ORF_Raw.read.count.D27                   <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,…
## $ ORF_Raw.read.count.D28                   <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,…
## $ ORF_Raw.read.count.D29                   <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,…
## $ ORF_Raw.read.count.D2                    <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,…
## $ ORF_Raw.read.count.D30                   <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,…
## $ ORF_Raw.read.count.D31                   <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,…
## $ ORF_Raw.read.count.D32                   <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,…
## $ ORF_Raw.read.count.D33                   <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,…
## $ ORF_Raw.read.count.D34                   <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,…
## $ ORF_Raw.read.count.D35                   <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,…
## $ ORF_Raw.read.count.D36                   <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,…
## $ ORF_Raw.read.count.D37                   <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,…
## $ ORF_Raw.read.count.D38                   <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,…
## $ ORF_Raw.read.count.D39                   <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,…
## $ ORF_Raw.read.count.D3                    <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,…
## $ ORF_Raw.read.count.D40                   <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,…
## $ ORF_Raw.read.count.D41                   <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,…
## $ ORF_Raw.read.count.D42                   <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,…
## $ ORF_Raw.read.count.D43                   <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,…
## $ ORF_Raw.read.count.D4                    <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,…
## $ ORF_Raw.read.count.D5                    <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,…
## $ ORF_Raw.read.count.D6                    <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,…
## $ ORF_Raw.read.count.D7                    <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,…
## $ ORF_Raw.read.count.D8                    <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,…
## $ ORF_Raw.read.count.D9                    <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,…
## $ ORF_Raw.base.count.C1                    <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,…
## $ ORF_Raw.base.count.C2                    <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,…
## $ ORF_Raw.base.count.C3                    <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,…
## $ ORF_Raw.base.count.C4                    <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,…
## $ ORF_Raw.base.count.C5                    <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,…
## $ ORF_Raw.base.count.C6                    <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,…
## $ ORF_Raw.base.count.D10                   <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,…
## $ ORF_Raw.base.count.D11                   <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,…
## $ ORF_Raw.base.count.D12                   <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,…
## $ ORF_Raw.base.count.D13                   <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,…
## $ ORF_Raw.base.count.D14                   <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,…
## $ ORF_Raw.base.count.D15                   <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,…
## $ ORF_Raw.base.count.D16                   <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,…
## $ ORF_Raw.base.count.D17                   <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,…
## $ ORF_Raw.base.count.D18                   <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,…
## $ ORF_Raw.base.count.D19                   <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,…
## $ ORF_Raw.base.count.D1                    <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,…
## $ ORF_Raw.base.count.D20                   <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,…
## $ ORF_Raw.base.count.D21                   <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,…
## $ ORF_Raw.base.count.D22                   <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,…
## $ ORF_Raw.base.count.D23                   <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,…
## $ ORF_Raw.base.count.D24                   <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,…
## $ ORF_Raw.base.count.D25                   <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,…
## $ ORF_Raw.base.count.D26                   <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,…
## $ ORF_Raw.base.count.D27                   <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,…
## $ ORF_Raw.base.count.D28                   <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,…
## $ ORF_Raw.base.count.D29                   <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,…
## $ ORF_Raw.base.count.D2                    <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,…
## $ ORF_Raw.base.count.D30                   <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,…
## $ ORF_Raw.base.count.D31                   <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,…
## $ ORF_Raw.base.count.D32                   <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,…
## $ ORF_Raw.base.count.D33                   <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,…
## $ ORF_Raw.base.count.D34                   <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,…
## $ ORF_Raw.base.count.D35                   <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,…
## $ ORF_Raw.base.count.D36                   <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,…
## $ ORF_Raw.base.count.D37                   <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,…
## $ ORF_Raw.base.count.D38                   <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,…
## $ ORF_Raw.base.count.D39                   <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,…
## $ ORF_Raw.base.count.D3                    <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,…
## $ ORF_Raw.base.count.D40                   <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,…
## $ ORF_Raw.base.count.D41                   <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,…
## $ ORF_Raw.base.count.D42                   <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,…
## $ ORF_Raw.base.count.D43                   <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,…
## $ ORF_Raw.base.count.D4                    <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,…
## $ ORF_Raw.base.count.D5                    <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,…
## $ ORF_Raw.base.count.D6                    <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,…
## $ ORF_Raw.base.count.D7                    <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,…
## $ ORF_Raw.base.count.D8                    <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,…
## $ ORF_Raw.base.count.D9                    <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,…
## $ ORF_TPM.D.1                              <dbl> 0.165, 0.000, 0.000, 0.000, 0…
## $ ORF_TPM.HC.1                             <dbl> 0.000, 0.229, 0.151, 0.000, 0…
## $ ORF_TPM.HC.2                             <dbl> 0.000, 0.000, 0.000, 0.298, 0…
## $ ORF_TPM.HC.3                             <dbl> 0.000, 0.000, 0.000, 0.000, 0…
## $ ORF_TPM.HC.4                             <dbl> 0.000, 0.000, 0.000, 0.000, 0…
## $ ORF_TPM.HC.5                             <dbl> 0.000, 0.000, 0.000, 0.000, 0…
## $ ORF_TPM.HC.6                             <dbl> 0.000, 0.000, 0.000, 0.000, 0…
## $ ORF_TPM.HC.7                             <dbl> 0.000, 0.000, 0.000, 0.000, 0…
## $ ORF_TPM.HC.8                             <dbl> 0.000, 0.000, 0.000, 0.000, 0…
## $ ORF_TPM.HC.9                             <dbl> 0.000, 0.000, 0.000, 0.000, 0…
## $ ORF_TPM.HC.10                            <dbl> 0.000, 0.000, 0.000, 0.000, 0…
## $ ORF_TPM.HC.11                            <dbl> 0.000, 0.000, 0.000, 0.000, 0…
## $ ORF_TPM.HC.12                            <dbl> 0.000, 0.230, 0.304, 0.000, 0…
## $ ORF_TPM.HC.13                            <dbl> 0.000, 0.000, 0.000, 0.000, 0…
## $ ORF_TPM.HC.14                            <dbl> 0.000, 0.000, 0.000, 0.000, 0…
## $ ORF_TPM.HC.15                            <dbl> 0.000, 0.000, 0.000, 0.000, 0…
## $ ORF_TPM.HC.16                            <dbl> 0.000, 0.000, 0.000, 0.000, 0…
## $ ORF_TPM.HC.17                            <dbl> 0.000, 0.000, 0.000, 0.156, 0…
## $ ORF_TPM.HC.18                            <dbl> 0.000, 0.000, 0.000, 0.000, 0…
## $ ORF_TPM.HC.19                            <dbl> 0.000, 0.000, 0.000, 0.000, 0…
## $ ORF_TPM.HC.20                            <dbl> 0.000, 0.000, 0.000, 0.178, 0…
## $ ORF_TPM.HC.21                            <dbl> 0.000, 0.000, 0.000, 0.141, 0…
## $ ORF_TPM.HC.22                            <dbl> 0.000, 0.000, 0.000, 0.000, 0…
## $ ORF_TPM.HC.23                            <dbl> 0.000, 0.000, 0.000, 0.154, 0…
## $ ORF_TPM.HC.24                            <dbl> 0.000, 0.000, 0.000, 0.000, 0…
## $ ORF_TPM.HV.1                             <dbl> 0.000, 0.000, 0.000, 0.000, 0…
## $ ORF_TPM.HV.2                             <dbl> 0.000, 0.000, 0.000, 0.000, 0…
## $ ORF_TPM.HV.3                             <dbl> 0.000, 0.000, 0.000, 0.000, 0…
## $ ORF_TPM.HV.4                             <dbl> 0.000, 0.000, 0.000, 0.000, 0…
## $ ORF_TPM.HV.5                             <dbl> 0.000, 0.000, 0.000, 0.000, 0…
## $ ORF_TPM.HV.6                             <dbl> 0.000, 0.000, 0.000, 0.000, 0…
## $ ORF_TPM.HV.7                             <dbl> 0.000, 0.000, 0.000, 0.000, 0…
## $ ORF_TPM.HV.8                             <dbl> 0.000, 0.000, 0.000, 0.000, 0…
## $ ORF_TPM.HV.9                             <dbl> 0.000, 0.000, 0.000, 0.000, 0…
## $ ORF_TPM.HV.10                            <dbl> 0.000, 0.000, 0.000, 0.000, 0…
## $ ORF_TPM.HV.11                            <dbl> 0.000, 0.000, 0.000, 0.000, 0…
## $ ORF_TPM.HV.12                            <dbl> 0.000, 0.000, 0.000, 0.165, 0…
## $ ORF_TPM.HV.13                            <dbl> 0.000, 0.000, 0.000, 0.000, 0…
## $ ORF_TPM.HV.14                            <dbl> 0.000, 0.000, 0.000, 0.000, 0…
## $ ORF_TPM.HV.15                            <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,…
## $ ORF_TPM.HV.16                            <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,…
## $ ORF_TPM.HV.17                            <dbl> 0.000, 0.000, 0.000, 0.000, 0…
## $ ORF_TPM.HV.18                            <dbl> 0.000, 0.000, 0.000, 0.434, 0…
## $ ORF_TPM.HV.19                            <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,…
## $ ORF_TPM.HV.20                            <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,…
## $ ORF_TPM.HV.21                            <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,…
## $ ORF_TPM.HV.22                            <dbl> 0.000, 0.000, 0.000, 0.000, 0…
## $ ORF_TPM.HV.23                            <dbl> 0.000, 0.000, 0.000, 0.000, 0…
## $ ORF_TPM.HV.24                            <dbl> 0.000, 0.000, 0.000, 0.000, 0…
## $ ORF_Coverage.D.1                         <dbl> 3.349, 0.000, 0.000, 0.000, 0…
## $ ORF_Coverage.HC.1                        <dbl> 0.000, 1.149, 0.780, 0.000, 0…
## $ ORF_Coverage.HC.2                        <dbl> 0.000, 0.000, 0.000, 1.684, 1…
## $ ORF_Coverage.HC.3                        <dbl> 0.000, 0.000, 0.000, 0.000, 0…
## $ ORF_Coverage.HC.4                        <dbl> 0.000, 0.000, 0.000, 0.000, 0…
## $ ORF_Coverage.HC.5                        <dbl> 0.000, 0.000, 0.000, 0.000, 0…
## $ ORF_Coverage.HC.6                        <dbl> 0.000, 0.000, 0.000, 0.000, 0…
## $ ORF_Coverage.HC.7                        <dbl> 0.000, 0.000, 0.000, 0.000, 0…
## $ ORF_Coverage.HC.8                        <dbl> 0.000, 0.000, 0.000, 0.000, 0…
## $ ORF_Coverage.HC.9                        <dbl> 0.000, 0.000, 0.000, 0.000, 0…
## $ ORF_Coverage.HC.10                       <dbl> 0.000, 0.000, 0.000, 0.000, 0…
## $ ORF_Coverage.HC.11                       <dbl> 0.000, 0.000, 0.000, 0.000, 0…
## $ ORF_Coverage.HC.12                       <dbl> 0.000, 1.172, 0.773, 0.000, 0…
## $ ORF_Coverage.HC.13                       <dbl> 0.000, 0.000, 0.000, 0.000, 0…
## $ ORF_Coverage.HC.14                       <dbl> 0.000, 0.000, 0.000, 0.000, 0…
## $ ORF_Coverage.HC.15                       <dbl> 0.000, 0.000, 0.000, 0.000, 0…
## $ ORF_Coverage.HC.16                       <dbl> 0.000, 0.000, 0.000, 0.000, 0…
## $ ORF_Coverage.HC.17                       <dbl> 0.000, 0.000, 0.000, 0.404, 0…
## $ ORF_Coverage.HC.18                       <dbl> 0.000, 0.000, 0.000, 0.000, 0…
## $ ORF_Coverage.HC.19                       <dbl> 0.000, 0.000, 0.000, 0.000, 0…
## $ ORF_Coverage.HC.20                       <dbl> 0.000, 0.000, 0.000, 0.860, 0…
## $ ORF_Coverage.HC.21                       <dbl> 0.000, 0.000, 0.000, 0.746, 0…
## $ ORF_Coverage.HC.22                       <dbl> 0.000, 0.000, 0.000, 0.000, 0…
## $ ORF_Coverage.HC.23                       <dbl> 0.000, 0.000, 0.000, 0.561, 0…
## $ ORF_Coverage.HC.24                       <dbl> 0.000, 0.000, 0.000, 0.000, 0…
## $ ORF_Coverage.HV.1                        <dbl> 0.000, 0.000, 0.000, 0.000, 0…
## $ ORF_Coverage.HV.2                        <dbl> 0.000, 0.000, 0.000, 0.000, 0…
## $ ORF_Coverage.HV.3                        <dbl> 0.000, 0.000, 0.000, 0.000, 0…
## $ ORF_Coverage.HV.4                        <dbl> 0.000, 0.000, 0.000, 0.000, 0…
## $ ORF_Coverage.HV.5                        <dbl> 0.000, 0.000, 0.000, 0.000, 0…
## $ ORF_Coverage.HV.6                        <dbl> 0.000, 0.000, 0.000, 0.000, 0…
## $ ORF_Coverage.HV.7                        <dbl> 0.000, 0.000, 0.000, 0.000, 0…
## $ ORF_Coverage.HV.8                        <dbl> 0.000, 0.000, 0.000, 0.000, 0…
## $ ORF_Coverage.HV.9                        <dbl> 0.000, 0.000, 0.000, 0.000, 0…
## $ ORF_Coverage.HV.10                       <dbl> 0.000, 0.000, 0.000, 0.000, 0…
## $ ORF_Coverage.HV.11                       <dbl> 0.000, 0.000, 0.000, 0.000, 0…
## $ ORF_Coverage.HV.12                       <dbl> 0.000, 0.000, 0.000, 0.096, 0…
## $ ORF_Coverage.HV.13                       <dbl> 0.000, 0.000, 0.000, 0.000, 0…
## $ ORF_Coverage.HV.14                       <dbl> 0.000, 0.000, 0.000, 0.000, 0…
## $ ORF_Coverage.HV.15                       <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,…
## $ ORF_Coverage.HV.16                       <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,…
## $ ORF_Coverage.HV.17                       <dbl> 0.000, 0.000, 0.000, 0.000, 0…
## $ ORF_Coverage.HV.18                       <dbl> 0.000, 0.000, 0.000, 1.974, 2…
## $ ORF_Coverage.HV.19                       <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,…
## $ ORF_Coverage.HV.20                       <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,…
## $ ORF_Coverage.HV.21                       <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,…
## $ ORF_Coverage.HV.22                       <dbl> 0.000, 0.000, 0.000, 0.000, 0…
## $ ORF_Coverage.HV.23                       <dbl> 0.000, 0.000, 0.000, 0.000, 0…
## $ ORF_Coverage.HV.24                       <dbl> 0.000, 0.000, 0.000, 0.000, 0…
## $ ORF_Num_Gi.D.1                           <dbl> 0.022435897, 0.000000000, 0.0…
## $ ORF_Num_Gi.HC.1                          <dbl> 0.000000000, 0.011494253, 0.0…
## $ ORF_Num_Gi.HC.2                          <dbl> 0.000000000, 0.000000000, 0.0…
## $ ORF_Num_Gi.HC.3                          <dbl> 0.000000000, 0.000000000, 0.0…
## $ ORF_Num_Gi.HC.4                          <dbl> 0.000000000, 0.000000000, 0.0…
## $ ORF_Num_Gi.HC.5                          <dbl> 0.000000000, 0.000000000, 0.0…
## $ ORF_Num_Gi.HC.6                          <dbl> 0.000000000, 0.000000000, 0.0…
## $ ORF_Num_Gi.HC.7                          <dbl> 0.000000000, 0.000000000, 0.0…
## $ ORF_Num_Gi.HC.8                          <dbl> 0.00000000, 0.00000000, 0.000…
## $ ORF_Num_Gi.HC.9                          <dbl> 0.000000000, 0.000000000, 0.0…
## $ ORF_Num_Gi.HC.10                         <dbl> 0.000000000, 0.000000000, 0.0…
## $ ORF_Num_Gi.HC.11                         <dbl> 0.00000000, 0.00000000, 0.000…
## $ ORF_Num_Gi.HC.12                         <dbl> 0.00000000, 0.01149425, 0.015…
## $ ORF_Num_Gi.HC.13                         <dbl> 0.00000000, 0.00000000, 0.000…
## $ ORF_Num_Gi.HC.14                         <dbl> 0.000000000, 0.000000000, 0.0…
## $ ORF_Num_Gi.HC.15                         <dbl> 0.000000000, 0.000000000, 0.0…
## $ ORF_Num_Gi.HC.16                         <dbl> 0.000000000, 0.000000000, 0.0…
## $ ORF_Num_Gi.HC.17                         <dbl> 0.000000000, 0.000000000, 0.0…
## $ ORF_Num_Gi.HC.18                         <dbl> 0.00000000, 0.00000000, 0.000…
## $ ORF_Num_Gi.HC.19                         <dbl> 0.00000000, 0.00000000, 0.000…
## $ ORF_Num_Gi.HC.20                         <dbl> 0.000000000, 0.000000000, 0.0…
## $ ORF_Num_Gi.HC.21                         <dbl> 0.000000000, 0.000000000, 0.0…
## $ ORF_Num_Gi.HC.22                         <dbl> 0.000000000, 0.000000000, 0.0…
## $ ORF_Num_Gi.HC.23                         <dbl> 0.000000000, 0.000000000, 0.0…
## $ ORF_Num_Gi.HC.24                         <dbl> 0.000000000, 0.000000000, 0.0…
## $ ORF_Num_Gi.HV.1                          <dbl> 0.00000000, 0.00000000, 0.000…
## $ ORF_Num_Gi.HV.2                          <dbl> 0.0000000, 0.0000000, 0.00000…
## $ ORF_Num_Gi.HV.3                          <dbl> 0.000000000, 0.000000000, 0.0…
## $ ORF_Num_Gi.HV.4                          <dbl> 0.00000000, 0.00000000, 0.000…
## $ ORF_Num_Gi.HV.5                          <dbl> 0.00000000, 0.00000000, 0.000…
## $ ORF_Num_Gi.HV.6                          <dbl> 0.000000000, 0.000000000, 0.0…
## $ ORF_Num_Gi.HV.7                          <dbl> 0.000000000, 0.000000000, 0.0…
## $ ORF_Num_Gi.HV.8                          <dbl> 0.000000000, 0.000000000, 0.0…
## $ ORF_Num_Gi.HV.9                          <dbl> 0.000000000, 0.000000000, 0.0…
## $ ORF_Num_Gi.HV.10                         <dbl> 0.000000000, 0.000000000, 0.0…
## $ ORF_Num_Gi.HV.11                         <dbl> 0.000000000, 0.000000000, 0.0…
## $ ORF_Num_Gi.HV.12                         <dbl> 0.000000000, 0.000000000, 0.0…
## $ ORF_Num_Gi.HV.13                         <dbl> 0.00000000, 0.00000000, 0.000…
## $ ORF_Num_Gi.HV.14                         <dbl> 0.000000000, 0.000000000, 0.0…
## $ ORF_Num_Gi.HV.15                         <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,…
## $ ORF_Num_Gi.HV.16                         <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,…
## $ ORF_Num_Gi.HV.17                         <dbl> 0.0000000, 0.0000000, 0.00000…
## $ ORF_Num_Gi.HV.18                         <dbl> 0.00000000, 0.00000000, 0.000…
## $ ORF_Num_Gi.HV.19                         <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,…
## $ ORF_Num_Gi.HV.20                         <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,…
## $ ORF_Num_Gi.HV.21                         <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,…
## $ ORF_Num_Gi.HV.22                         <dbl> 0.00000000, 0.00000000, 0.000…
## $ ORF_Num_Gi.HV.23                         <dbl> 0.000000000, 0.000000000, 0.0…
## $ ORF_Num_Gi.HV.24                         <dbl> 0.00000000, 0.00000000, 0.000…
## $ ORF_Num_Gi_pc.D.1                        <dbl> 1.651063e-07, 0.000000e+00, 0…
## $ ORF_Num_Gi_pc.HC.1                       <dbl> 0.000000e+00, 2.288336e-07, 1…
## $ ORF_Num_Gi_pc.HC.2                       <dbl> 0.000000e+00, 0.000000e+00, 0…
## $ ORF_Num_Gi_pc.HC.3                       <dbl> 0.000000e+00, 0.000000e+00, 0…
## $ ORF_Num_Gi_pc.HC.4                       <dbl> 0.000000e+00, 0.000000e+00, 0…
## $ ORF_Num_Gi_pc.HC.5                       <dbl> 0.000000e+00, 0.000000e+00, 0…
## $ ORF_Num_Gi_pc.HC.6                       <dbl> 0.000000e+00, 0.000000e+00, 0…
## $ ORF_Num_Gi_pc.HC.7                       <dbl> 0.000000e+00, 0.000000e+00, 0…
## $ ORF_Num_Gi_pc.HC.8                       <dbl> 0.000000e+00, 0.000000e+00, 0…
## $ ORF_Num_Gi_pc.HC.9                       <dbl> 0.000000e+00, 0.000000e+00, 0…
## $ ORF_Num_Gi_pc.HC.10                      <dbl> 0.000000e+00, 0.000000e+00, 0…
## $ ORF_Num_Gi_pc.HC.11                      <dbl> 0.000000e+00, 0.000000e+00, 0…
## $ ORF_Num_Gi_pc.HC.12                      <dbl> 0.000000e+00, 2.304375e-07, 3…
## $ ORF_Num_Gi_pc.HC.13                      <dbl> 0.000000e+00, 0.000000e+00, 0…
## $ ORF_Num_Gi_pc.HC.14                      <dbl> 0.000000e+00, 0.000000e+00, 0…
## $ ORF_Num_Gi_pc.HC.15                      <dbl> 0.000000e+00, 0.000000e+00, 0…
## $ ORF_Num_Gi_pc.HC.16                      <dbl> 0.000000e+00, 0.000000e+00, 0…
## $ ORF_Num_Gi_pc.HC.17                      <dbl> 0.000000e+00, 0.000000e+00, 0…
## $ ORF_Num_Gi_pc.HC.18                      <dbl> 0.000000e+00, 0.000000e+00, 0…
## $ ORF_Num_Gi_pc.HC.19                      <dbl> 0.000000e+00, 0.000000e+00, 0…
## $ ORF_Num_Gi_pc.HC.20                      <dbl> 0.000000e+00, 0.000000e+00, 0…
## $ ORF_Num_Gi_pc.HC.21                      <dbl> 0.000000e+00, 0.000000e+00, 0…
## $ ORF_Num_Gi_pc.HC.22                      <dbl> 0.000000e+00, 0.000000e+00, 0…
## $ ORF_Num_Gi_pc.HC.23                      <dbl> 0.000000e+00, 0.000000e+00, 0…
## $ ORF_Num_Gi_pc.HC.24                      <dbl> 0.000000e+00, 0.000000e+00, 0…
## $ ORF_Num_Gi_pc.HV.1                       <dbl> 0.000000e+00, 0.000000e+00, 0…
## $ ORF_Num_Gi_pc.HV.2                       <dbl> 0.000000e+00, 0.000000e+00, 0…
## $ ORF_Num_Gi_pc.HV.3                       <dbl> 0.000000e+00, 0.000000e+00, 0…
## $ ORF_Num_Gi_pc.HV.4                       <dbl> 0.000000e+00, 0.000000e+00, 0…
## $ ORF_Num_Gi_pc.HV.5                       <dbl> 0.000000e+00, 0.000000e+00, 0…
## $ ORF_Num_Gi_pc.HV.6                       <dbl> 0.000000e+00, 0.000000e+00, 0…
## $ ORF_Num_Gi_pc.HV.7                       <dbl> 0.000000e+00, 0.000000e+00, 0…
## $ ORF_Num_Gi_pc.HV.8                       <dbl> 0.000000e+00, 0.000000e+00, 0…
## $ ORF_Num_Gi_pc.HV.9                       <dbl> 0.000000e+00, 0.000000e+00, 0…
## $ ORF_Num_Gi_pc.HV.10                      <dbl> 0.000000e+00, 0.000000e+00, 0…
## $ ORF_Num_Gi_pc.HV.11                      <dbl> 0.000000e+00, 0.000000e+00, 0…
## $ ORF_Num_Gi_pc.HV.12                      <dbl> 0.000000e+00, 0.000000e+00, 0…
## $ ORF_Num_Gi_pc.HV.13                      <dbl> 0.000000e+00, 0.000000e+00, 0…
## $ ORF_Num_Gi_pc.HV.14                      <dbl> 0.000000e+00, 0.000000e+00, 0…
## $ ORF_Num_Gi_pc.HV.15                      <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,…
## $ ORF_Num_Gi_pc.HV.16                      <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,…
## $ ORF_Num_Gi_pc.HV.17                      <dbl> 0.000000e+00, 0.000000e+00, 0…
## $ ORF_Num_Gi_pc.HV.18                      <dbl> 0.000000e+00, 0.000000e+00, 0…
## $ ORF_Num_Gi_pc.HV.19                      <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,…
## $ ORF_Num_Gi_pc.HV.20                      <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,…
## $ ORF_Num_Gi_pc.HV.21                      <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,…
## $ ORF_Num_Gi_pc.HV.22                      <dbl> 0.000000e+00, 0.000000e+00, 0…
## $ ORF_Num_Gi_pc.HV.23                      <dbl> 0.000000e+00, 0.000000e+00, 0…
## $ ORF_Num_Gi_pc.HV.24                      <dbl> 0.000000e+00, 0.000000e+00, 0…
## $ ORF_Raw.read.count.D.1                   <dbl> 7, 0, 0, 0, 0, 0, 0, 0, 0, 4,…
## $ ORF_Raw.read.count.HC.1                  <dbl> 0, 2, 1, 0, 0, 0, 0, 0, 0, 0,…
## $ ORF_Raw.read.count.HC.2                  <dbl> 0, 0, 0, 2, 2, 0, 0, 0, 0, 0,…
## $ ORF_Raw.read.count.HC.3                  <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,…
## $ ORF_Raw.read.count.HC.4                  <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,…
## $ ORF_Raw.read.count.HC.5                  <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,…
## $ ORF_Raw.read.count.HC.6                  <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,…
## $ ORF_Raw.read.count.HC.7                  <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,…
## $ ORF_Raw.read.count.HC.8                  <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,…
## $ ORF_Raw.read.count.HC.9                  <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,…
## $ ORF_Raw.read.count.HC.10                 <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,…
## $ ORF_Raw.read.count.HC.11                 <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,…
## $ ORF_Raw.read.count.HC.12                 <dbl> 0, 2, 2, 0, 0, 0, 0, 0, 0, 0,…
## $ ORF_Raw.read.count.HC.13                 <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,…
## $ ORF_Raw.read.count.HC.14                 <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,…
## $ ORF_Raw.read.count.HC.15                 <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,…
## $ ORF_Raw.read.count.HC.16                 <dbl> 0, 0, 0, 0, 0, 2, 4, 0, 0, 0,…
## $ ORF_Raw.read.count.HC.17                 <dbl> 0, 0, 0, 1, 1, 0, 0, 0, 0, 0,…
## $ ORF_Raw.read.count.HC.18                 <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,…
## $ ORF_Raw.read.count.HC.19                 <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,…
## $ ORF_Raw.read.count.HC.20                 <dbl> 0, 0, 0, 1, 1, 0, 0, 0, 0, 0,…
## $ ORF_Raw.read.count.HC.21                 <dbl> 0, 0, 0, 1, 1, 0, 0, 0, 0, 0,…
## $ ORF_Raw.read.count.HC.22                 <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,…
## $ ORF_Raw.read.count.HC.23                 <dbl> 0, 0, 0, 1, 1, 0, 0, 0, 0, 0,…
## $ ORF_Raw.read.count.HC.24                 <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,…
## $ ORF_Raw.read.count.HV.1                  <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,…
## $ ORF_Raw.read.count.HV.2                  <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,…
## $ ORF_Raw.read.count.HV.3                  <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,…
## $ ORF_Raw.read.count.HV.4                  <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,…
## $ ORF_Raw.read.count.HV.5                  <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,…
## $ ORF_Raw.read.count.HV.6                  <dbl> 0, 0, 0, 0, 0, 0, 1, 0, 0, 0,…
## $ ORF_Raw.read.count.HV.7                  <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,…
## $ ORF_Raw.read.count.HV.8                  <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,…
## $ ORF_Raw.read.count.HV.9                  <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,…
## $ ORF_Raw.read.count.HV.10                 <dbl> 0, 0, 0, 0, 0, 2, 2, 0, 0, 0,…
## $ ORF_Raw.read.count.HV.11                 <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,…
## $ ORF_Raw.read.count.HV.12                 <dbl> 0, 0, 0, 1, 1, 0, 0, 0, 0, 0,…
## $ ORF_Raw.read.count.HV.13                 <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,…
## $ ORF_Raw.read.count.HV.14                 <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,…
## $ ORF_Raw.read.count.HV.15                 <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,…
## $ ORF_Raw.read.count.HV.16                 <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,…
## $ ORF_Raw.read.count.HV.17                 <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,…
## $ ORF_Raw.read.count.HV.18                 <dbl> 0, 0, 0, 3, 3, 0, 0, 0, 0, 0,…
## $ ORF_Raw.read.count.HV.19                 <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,…
## $ ORF_Raw.read.count.HV.20                 <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,…
## $ ORF_Raw.read.count.HV.21                 <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,…
## $ ORF_Raw.read.count.HV.22                 <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,…
## $ ORF_Raw.read.count.HV.23                 <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,…
## $ ORF_Raw.read.count.HV.24                 <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,…
## $ ORF_Raw.base.count.D.1                   <dbl> 1045, 0, 0, 0, 0, 0, 0, 0, 0,…
## $ ORF_Raw.base.count.HC.1                  <dbl> 0, 200, 103, 0, 0, 0, 0, 0, 0…
## $ ORF_Raw.base.count.HC.2                  <dbl> 0, 0, 0, 192, 174, 0, 0, 0, 0…
## $ ORF_Raw.base.count.HC.3                  <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,…
## $ ORF_Raw.base.count.HC.4                  <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,…
## $ ORF_Raw.base.count.HC.5                  <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,…
## $ ORF_Raw.base.count.HC.6                  <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,…
## $ ORF_Raw.base.count.HC.7                  <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,…
## $ ORF_Raw.base.count.HC.8                  <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,…
## $ ORF_Raw.base.count.HC.9                  <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,…
## $ ORF_Raw.base.count.HC.10                 <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,…
## $ ORF_Raw.base.count.HC.11                 <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,…
## $ ORF_Raw.base.count.HC.12                 <dbl> 0, 204, 102, 0, 0, 0, 0, 0, 0…
## $ ORF_Raw.base.count.HC.13                 <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,…
## $ ORF_Raw.base.count.HC.14                 <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,…
## $ ORF_Raw.base.count.HC.15                 <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,…
## $ ORF_Raw.base.count.HC.16                 <dbl> 0, 0, 0, 0, 0, 142, 424, 0, 0…
## $ ORF_Raw.base.count.HC.17                 <dbl> 0, 0, 0, 46, 137, 0, 0, 0, 0,…
## $ ORF_Raw.base.count.HC.18                 <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,…
## $ ORF_Raw.base.count.HC.19                 <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,…
## $ ORF_Raw.base.count.HC.20                 <dbl> 0, 0, 0, 98, 86, 0, 0, 0, 0, …
## $ ORF_Raw.base.count.HC.21                 <dbl> 0, 0, 0, 85, 99, 0, 0, 0, 0, …
## $ ORF_Raw.base.count.HC.22                 <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,…
## $ ORF_Raw.base.count.HC.23                 <dbl> 0, 0, 0, 64, 120, 0, 0, 0, 0,…
## $ ORF_Raw.base.count.HC.24                 <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,…
## $ ORF_Raw.base.count.HV.1                  <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,…
## $ ORF_Raw.base.count.HV.2                  <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,…
## $ ORF_Raw.base.count.HV.3                  <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,…
## $ ORF_Raw.base.count.HV.4                  <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,…
## $ ORF_Raw.base.count.HV.5                  <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,…
## $ ORF_Raw.base.count.HV.6                  <dbl> 0, 0, 0, 0, 0, 0, 150, 0, 0, …
## $ ORF_Raw.base.count.HV.7                  <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,…
## $ ORF_Raw.base.count.HV.8                  <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,…
## $ ORF_Raw.base.count.HV.9                  <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,…
## $ ORF_Raw.base.count.HV.10                 <dbl> 0, 0, 0, 0, 0, 4, 263, 0, 0, …
## $ ORF_Raw.base.count.HV.11                 <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,…
## $ ORF_Raw.base.count.HV.12                 <dbl> 0, 0, 0, 11, 106, 0, 0, 0, 0,…
## $ ORF_Raw.base.count.HV.13                 <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,…
## $ ORF_Raw.base.count.HV.14                 <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,…
## $ ORF_Raw.base.count.HV.15                 <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,…
## $ ORF_Raw.base.count.HV.16                 <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,…
## $ ORF_Raw.base.count.HV.17                 <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,…
## $ ORF_Raw.base.count.HV.18                 <dbl> 0, 0, 0, 225, 327, 0, 0, 0, 0…
## $ ORF_Raw.base.count.HV.19                 <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,…
## $ ORF_Raw.base.count.HV.20                 <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,…
## $ ORF_Raw.base.count.HV.21                 <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,…
## $ ORF_Raw.base.count.HV.22                 <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,…
## $ ORF_Raw.base.count.HV.23                 <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,…
## $ ORF_Raw.base.count.HV.24                 <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,…
```


```r
combined_full %>% 
  head(n = 1000) %>% 
  rbind(combined_full %>% 
          tail(n = 1000) ) %>% 
  write_tsv("~/Desktop/test_comb_h1_c1.tsv")
```


```r
combined_full %>% 
  filter(!is.na(PathoFact_AMR_ARG) | !is.na(rgi_CARD_ARO) | !is.na(PathoFact_Tox_Toxin_classification)) %>% 
  write_tsv("~/Desktop/AMR_TOX_comb_h1_c1.tsv")

combined_full %>% 
  filter(!is.na(PathoFact_AMR_ARG) | !is.na(rgi_CARD_ARO)) %>% 
  write_tsv("~/Desktop/AMR_comb_h1_c1.tsv")
```


Generate phyloseq objects:


sample_data:


```r
meta 
```

```
##            sample  input filtered denoisedF denoisedR merged tabled filtered_pc
## D1     CR-15-S292  23558    23449     23410     23355  22970  22970       0.995
## D14     CR-17-S97  35320    35021     34992     34870  34419  34419       0.992
## D21    CR-19-S341  17616    17528     17501     17477  17246  17246       0.995
## D28    CR-21-S126   5791     5764      5753      5734   5585   5585       0.995
## D35    CR-46-S191  24284    24178     24128     24062  23531  23531       0.996
## D36    CR-64-S324  10847    10801     10784     10726  10465  10465       0.996
## D43      D-0-S154  11484    11428     11279     11256  10248  10248       0.995
## HC_1   TR1-61-S47 238630   238248    237846    237395 231166 231166       0.998
## HC_2  TR2-61-S196  24392    24363     24296     24257  22096  22096       0.999
## HC_3   TR3-61-S45  49196    49111     49012     48940  47064  47064       0.998
## HC_4   TR4-61-S11  28442    28406     28348     28277  27014  27014       0.999
## HC_5  TR5-61-S142    245      244       240       237    219    219       0.996
## HC_6  TR6-61-S186 252856   252414    252020    251515 245450 245450       0.998
## HC_7  TR1-63-S100     32       31        31        30     26     26       0.969
## HC_8  TR2-63-S160  36811    36755     36712     36601  35561  35561       0.998
## HC_9  TR3-63-S206  26658    26606     26549     26501  25740  25740       0.998
## HC_10  TR4-63-S74  41734    41656     41517     41407  39936  39936       0.998
## HC_11  TR5-63-S82  43222    43149     43078     42988  41415  41415       0.998
## D7    TR1-15-S168    153      151       144       115     79     79       0.987
## D8    TR1-17-S246  16895    16833     16807     16749  16429  16429       0.996
## D15   TR1-19-S305   9274     9256      9239      9203   8959   8959       0.998
## D22   TR1-21-S207  20083    20061     20043     19956  19419  19419       0.999
## D29   TR1-46-S170  32752    32656     32579     32459  31595  31595       0.997
## HC_12 TR6-63-S195  34779    34728     34672     34605  33178  33178       0.999
## D42   TR1-64-S263  31538    31375     31337     31248  30747  30747       0.995
## HC_13 TR1-65-S128     60       59        53        51     36     36       0.983
## HC_14 TR2-65-S158  43440    43349     43294     43167  41820  41820       0.998
## D6    TR2-15-S288  31566    31486     31385     31307  30381  30381       0.997
## D9    TR2-17-S248  31490    31404     31369     31205  30516  30516       0.997
## D16   TR2-19-S318  10065    10054     10040     10015   9814   9814       0.999
## D23   TR2-21-S228  15593    15567     15543     15456  14959  14959       0.998
## HC_15  TR3-65-S92  32513    32447     32308     32211  29961  29961       0.998
## HC_16 TR4-65-S155     27       26        25        24     15     15       0.963
## HC_17 TR5-65-S187     35       35        32        32     25     25       1.000
## HC_18  TR6-65-S94  33692    33627     33496     33425  30776  30776       0.998
## HC_19 TR1-69-S145  40431    40343     40200     40139  38039  38039       0.998
## HC_20  TR2-69-S57 178654   178284    177932    177646 174558 174558       0.998
## D30   TR2-46-S239  25929    25888     25862     25757  25098  25098       0.998
## HC_21 TR3-69-S139  30574    30520     30460     30405  29510  29510       0.998
## HC_22 TR4-69-S197  25982    25946     25869     25853  24646  24646       0.999
## HC_23 TR5-69-S107     23       23        22        22     19     19       1.000
## HC_24 TR6-69-S174  44451    44372     44282     44188  42841  42841       0.998
## D41   TR2-64-S216  21707    21663     21627     21526  20850  20850       0.998
## D5    TR3-15-S338  17747    17682     17655     17595  17213  17213       0.996
## D10   TR3-17-S222  14735    14652     14631     14571  14285  14285       0.994
## D17   TR3-19-S266  32415    32239     32191     32099  31521  31521       0.995
## D24   TR3-21-S234  23216    23054     23032     22946  22495  22495       0.993
## HV_1   TR1-28-S35  34188    34109     34047     34001  32656  32656       0.998
## HV_2   TR2-28-S64  25719    25668     25602     25529  23955  23955       0.998
## HV_3   TR3-28-S38 833850   832269    830661    829156 814071 814071       0.998
## HV_4  TR4-28-S138    119      119       112       111     94     94       1.000
## HV_5    TR5-28-S5  31023    30972     30868     30821  29212  29212       0.998
## HV_6   TR6-28-S80  28872    28808     28729     28664  27021  27021       0.998
## D31   TR3-46-S224  14862    14757     14738     14685  14391  14391       0.993
## D40   TR3-64-S190  28368    28173     28084     28039  27343  27343       0.993
## D4    TR4-15-S218  21839    21709     21687     21638  21317  21317       0.994
## D11   TR4-17-S114  28993    28752     28709     28624  28097  28097       0.992
## D18   TR4-19-S270  28080    27827     27811     27758  27315  27315       0.991
## D25   TR4-21-S254  21868    21735     21727     21684  21405  21405       0.994
## HV_7   TR1-30-S76  36509    36431     36343     36255  35484  35484       0.998
## HV_8   TR2-30-S73  38226    38119     37993     37920  36739  36739       0.997
## HV_9  TR3-30-S190  35109    35014     34904     34824  33704  33704       0.997
## HV_10  TR4-30-S48  26568    26495     26335     26274  24097  24097       0.997
## HV_11  TR5-30-S62  33028    32966     32872     32818  31102  31102       0.998
## HV_12 TR6-30-S144  31410    31362     31309     31213  29786  29786       0.998
## D32   TR4-46-S215  25294    25090     25053     24965  24729  24729       0.992
## D39   TR4-64-S220  15761    15686     15672     15607  15483  15483       0.995
## C_1   TR5-13-S103  25912    25816     25770     25681  24848  24848       0.996
## D3    TR5-15-S358  22887    22773     22768     22682  22378  22378       0.995
## D12    TR5-17-S82  81610    81154     81069     80828  79535  79535       0.994
## D19    TR5-19-S87  35069    34934     34910     34783  34135  34135       0.996
## C_6   TR5-20-S208  24461    24365     24357     24197  23598  23598       0.996
## D26   TR5-21-S185  35805    35626     35602     35513  34882  34882       0.995
## D33   TR5-46-S310  20221    20118     20093     19985  19794  19794       0.995
## HV_13 TR1-32-S111     32       32        30        27     23     23       1.000
## HV_14  TR2-32-S90  33510    33387     33329     33287  32683  32683       0.996
## HV_15 TR3-32-S211  24657    24582     24552     24494  24072  24072       0.997
## HV_16  TR4-32-S17  44815    44659     44563     44486  43879  43879       0.997
## HV_17  TR5-32-S56  23120    23079     23027     22987  21785  21785       0.998
## D38   TR5-64-S141  22049    21937     21916     21817  21659  21659       0.995
## HV_18  TR6-32-S39  30398    30364     30310     30240  28866  28866       0.999
## D_1      D-2-S164  29089    29028     28917     28759  27917  27917       0.998
## D2    TR6-15-S149  21435    21335     21298     21224  20816  20816       0.995
## D13   TR6-17-S138  20193    20062     20030     19961  19569  19569       0.994
## D20   TR6-19-S137  25289    25165     25122     24996  24337  24337       0.995
## D27   TR6-21-S107   4866     4830      4824      4806   4681   4681       0.993
## HV_19 TR1-36-S198  34485    34334     34287     33498  28774  28774       0.996
## HV_20 TR2-36-S183  46398    46202     46139     45199  39386  39386       0.996
## HV_21   TR3-36-S8  19943    19861     19755     19696  18108  18108       0.996
## HV_22  TR4-36-S12  31375    31261     31220     31194  30717  30717       0.996
## HV_23  TR5-36-S85  33570    33504     33364     33301  31257  31257       0.998
## D34   TR6-46-S176  49716    49525     49414     49275  48250  48250       0.996
## D37   TR6-64-S146  24475    24302     24256     24192  23730  23730       0.993
## HV_24 TR6-36-S124     37       37        34        32     27     27       1.000
## C_2        TR5-63     NA       NA        NA        NA     NA     NA          NA
## C_3        TR4-63     NA       NA        NA        NA     NA     NA          NA
## C_4        TR3-63     NA       NA        NA        NA     NA     NA          NA
## C_5        TR2-63     NA       NA        NA        NA     NA     NA          NA
##       denoisedF_pc denoisedR_pc merged_pc filtered_merged_pc input_merged_pc
## D1           0.998        0.996     0.981              0.980           0.975
## D14          0.999        0.996     0.984              0.983           0.974
## D21          0.998        0.997     0.985              0.984           0.979
## D28          0.998        0.995     0.971              0.969           0.964
## D35          0.998        0.995     0.975              0.973           0.969
## D36          0.998        0.993     0.970              0.969           0.965
## D43          0.987        0.985     0.909              0.897           0.892
## HC_1         0.998        0.996     0.972              0.970           0.969
## HC_2         0.997        0.996     0.909              0.907           0.906
## HC_3         0.998        0.997     0.960              0.958           0.957
## HC_4         0.998        0.995     0.953              0.951           0.950
## HC_5         0.984        0.971     0.912              0.898           0.894
## HC_6         0.998        0.996     0.974              0.972           0.971
## HC_7         1.000        0.968     0.839              0.839           0.812
## HC_8         0.999        0.996     0.969              0.968           0.966
## HC_9         0.998        0.996     0.970              0.967           0.966
## HC_10        0.997        0.994     0.962              0.959           0.957
## HC_11        0.998        0.996     0.961              0.960           0.958
## D7           0.954        0.762     0.549              0.523           0.516
## D8           0.998        0.995     0.978              0.976           0.972
## D15          0.998        0.994     0.970              0.968           0.966
## D22          0.999        0.995     0.969              0.968           0.967
## D29          0.998        0.994     0.970              0.968           0.965
## HC_12        0.998        0.996     0.957              0.955           0.954
## D42          0.999        0.996     0.981              0.980           0.975
## HC_13        0.898        0.864     0.679              0.610           0.600
## HC_14        0.999        0.996     0.966              0.965           0.963
## D6           0.997        0.994     0.968              0.965           0.962
## D9           0.999        0.994     0.973              0.972           0.969
## D16          0.999        0.996     0.977              0.976           0.975
## D23          0.998        0.993     0.962              0.961           0.959
## HC_15        0.996        0.993     0.927              0.923           0.922
## HC_16        0.962        0.923     0.600              0.577           0.556
## HC_17        0.914        0.914     0.781              0.714           0.714
## HC_18        0.996        0.994     0.919              0.915           0.913
## HC_19        0.996        0.995     0.946              0.943           0.941
## HC_20        0.998        0.996     0.981              0.979           0.977
## D30          0.999        0.995     0.970              0.969           0.968
## HC_21        0.998        0.996     0.969              0.967           0.965
## HC_22        0.997        0.996     0.953              0.950           0.949
## HC_23        0.957        0.957     0.864              0.826           0.826
## HC_24        0.998        0.996     0.967              0.965           0.964
## D41          0.998        0.994     0.964              0.962           0.961
## D5           0.998        0.995     0.975              0.973           0.970
## D10          0.999        0.994     0.976              0.975           0.969
## D17          0.999        0.996     0.979              0.978           0.972
## D24          0.999        0.995     0.977              0.976           0.969
## HV_1         0.998        0.997     0.959              0.957           0.955
## HV_2         0.997        0.995     0.936              0.933           0.931
## HV_3         0.998        0.996     0.980              0.978           0.976
## HV_4         0.941        0.933     0.839              0.790           0.790
## HV_5         0.997        0.995     0.946              0.943           0.942
## HV_6         0.997        0.995     0.941              0.938           0.936
## D31          0.999        0.995     0.976              0.975           0.968
## D40          0.997        0.995     0.974              0.971           0.964
## D4           0.999        0.997     0.983              0.982           0.976
## D11          0.999        0.996     0.979              0.977           0.969
## D18          0.999        0.998     0.982              0.982           0.973
## D25          1.000        0.998     0.985              0.985           0.979
## HV_7         0.998        0.995     0.976              0.974           0.972
## HV_8         0.997        0.995     0.967              0.964           0.961
## HV_9         0.997        0.995     0.966              0.963           0.960
## HV_10        0.994        0.992     0.915              0.909           0.907
## HV_11        0.997        0.996     0.946              0.943           0.942
## HV_12        0.998        0.995     0.951              0.950           0.948
## D32          0.999        0.995     0.987              0.986           0.978
## D39          0.999        0.995     0.988              0.987           0.982
## C_1          0.998        0.995     0.964              0.963           0.959
## D3           1.000        0.996     0.983              0.983           0.978
## D12          0.999        0.996     0.981              0.980           0.975
## D19          0.999        0.996     0.978              0.977           0.973
## C_6          1.000        0.993     0.969              0.969           0.965
## D26          0.999        0.997     0.980              0.979           0.974
## D33          0.999        0.993     0.985              0.984           0.979
## HV_13        0.938        0.844     0.767              0.719           0.719
## HV_14        0.998        0.997     0.981              0.979           0.975
## HV_15        0.999        0.996     0.980              0.979           0.976
## HV_16        0.998        0.996     0.985              0.983           0.979
## HV_17        0.998        0.996     0.946              0.944           0.942
## D38          0.999        0.995     0.988              0.987           0.982
## HV_18        0.998        0.996     0.952              0.951           0.950
## D_1          0.996        0.991     0.965              0.962           0.960
## D2           0.998        0.995     0.977              0.976           0.971
## D13          0.998        0.995     0.977              0.975           0.969
## D20          0.998        0.993     0.969              0.967           0.962
## D27          0.999        0.995     0.970              0.969           0.962
## HV_19        0.999        0.976     0.839              0.838           0.834
## HV_20        0.999        0.978     0.854              0.852           0.849
## HV_21        0.995        0.992     0.917              0.912           0.908
## HV_22        0.999        0.998     0.984              0.983           0.979
## HV_23        0.996        0.994     0.937              0.933           0.931
## D34          0.998        0.995     0.976              0.974           0.971
## D37          0.998        0.995     0.978              0.976           0.970
## HV_24        0.919        0.865     0.794              0.730           0.730
## C_2             NA           NA        NA                 NA              NA
## C_3             NA           NA        NA                 NA              NA
## C_4             NA           NA        NA                 NA              NA
## C_5             NA           NA        NA                 NA              NA
##       tabled_joined chimera_out length_filtered tabled_pc chimera_out_pc
## D1            22970       22102           22098         1           0.96
## D14           34419       32366           32350         1           0.94
## D21           17246       16865           16865         1           0.98
## D28            5585        5498            5498         1           0.98
## D35           23531       22823           22821         1           0.97
## D36           10465       10250           10250         1           0.98
## D43           10248       10134           10134         1           0.99
## HC_1         231166      188953          188948         1           0.82
## HC_2          22096       10235           10235         1           0.46
## HC_3          47064       36864           36862         1           0.78
## HC_4          27014       19750           19745         1           0.73
## HC_5            219         183             183         1           0.84
## HC_6         245450      201757          201745         1           0.82
## HC_7             26          24              24         1           0.92
## HC_8          35561       28863           28862         1           0.81
## HC_9          25740       20454           20454         1           0.79
## HC_10         39936       32773           32768         1           0.82
## HC_11         41415       31488           31486         1           0.76
## D7               79          79              79         1           1.00
## D8            16429       16233           16233         1           0.99
## D15            8959        8815            8815         1           0.98
## D22           19419       18970           18970         1           0.98
## D29           31595       30745           30745         1           0.97
## HC_12         33178       23813           23807         1           0.72
## D42           30747       30193           30193         1           0.98
## HC_13            36          29              29         1           0.81
## HC_14         41820       32880           32880         1           0.79
## D6            30381       28832           28832         1           0.95
## D9            30516       30071           30071         1           0.99
## D16            9814        9662            9662         1           0.98
## D23           14959       14630           14630         1           0.98
## HC_15         29961       17134           17133         1           0.57
## HC_16            15          14              14         1           0.93
## HC_17            25          23              23         1           0.92
## HC_18         30776       15573           15573         1           0.51
## HC_19         38039       24528           24524         1           0.64
## HC_20        174558      149775          149775         1           0.86
## D30           25098       24475           24475         1           0.98
## HC_21         29510       23544           23544         1           0.80
## HC_22         24646       15361           15359         1           0.62
## HC_23            19          17              17         1           0.89
## HC_24         42841       35074           35074         1           0.82
## D41           20850       20344           20344         1           0.98
## D5            17213       16817           16817         1           0.98
## D10           14285       14059           14059         1           0.98
## D17           31521       30971           30970         1           0.98
## D24           22495       21882           21876         1           0.97
## HV_1          32656       24739           24739         1           0.76
## HV_2          23955       14713           14713         1           0.61
## HV_3         814071      695560          695530         1           0.85
## HV_4             94          83              83         1           0.88
## HV_5          29212       18666           18666         1           0.64
## HV_6          27021       17152           17151         1           0.63
## D31           14391       14095           14095         1           0.98
## D40           27343       26205           26202         1           0.96
## D4            21317       20936           20936         1           0.98
## D11           28097       27386           27381         1           0.97
## D18           27315       26523           26523         1           0.97
## D25           21405       20811           20805         1           0.97
## HV_7          35484       28461           28461         1           0.80
## HV_8          36739       25840           25840         1           0.70
## HV_9          33704       26093           26093         1           0.77
## HV_10         24097       11843           11843         1           0.49
## HV_11         31102       21287           21287         1           0.68
## HV_12         29786       21835           21832         1           0.73
## D32           24729       23673           23673         1           0.96
## D39           15483       15160           15160         1           0.98
## C_1           24848       24314           24314         1           0.98
## D3            22378       21686           21680         1           0.97
## D12           79535       76689           76678         1           0.96
## D19           34135       33512           33510         1           0.98
## C_6           23598       22590           22589         1           0.96
## D26           34882       32585           32568         1           0.93
## D33           19794       19015           19015         1           0.96
## HV_13            23          18              18         1           0.78
## HV_14         32683       25240           25240         1           0.77
## HV_15         24072       19576           19576         1           0.81
## HV_16         43879       34419           34419         1           0.78
## HV_17         21785       15674           15673         1           0.72
## D38           21659       21080           21080         1           0.97
## HV_18         28866       21084           21084         1           0.73
## D_1           27917       26104           26104         1           0.94
## D2            20816       20429           20429         1           0.98
## D13           19569       19164           19164         1           0.98
## D20           24337       23661           23659         1           0.97
## D27            4681        4595            4595         1           0.98
## HV_19         28774       22660           22660         1           0.79
## HV_20         39386       33722           33722         1           0.86
## HV_21         18108        8748            8748         1           0.48
## HV_22         30717       23394           23394         1           0.76
## HV_23         31257       19614           19614         1           0.63
## D34           48250       46576           46575         1           0.97
## D37           23730       23096           23096         1           0.97
## HV_24            27          20              20         1           0.74
## C_2              NA          NA              NA        NA             NA
## C_3              NA          NA              NA        NA             NA
## C_4              NA          NA              NA        NA             NA
## C_5              NA          NA              NA        NA             NA
##       length_filtered_pc Sample_description I7_Index_ID    index I5_Index_ID
## D1                     1              CR-15      N716-D ACTCGCTA      S517-D
## D14                    1              CR-17      N701-A TAAGGCGA      S513-D
## D21                    1              CR-19      N723-D TAGCGCTC      S518-D
## D28                    1              CR-21      N704-A TCCTGAGC      S520-D
## D35                    1              CR-46      N715-A ATCTCAGG      S521-D
## D36                    1              CR-64      N721-D TACGCTGC      S517-D
## D43                    1                D-0      N710-A CGAGGCTG      S515-D
## HC_1                   1             TR1-61      A-N706     <NA>      A-S510
## HC_2                   1             TR2-61      D-N716     <NA>      A-S506
## HC_3                   1             TR3-61      A-N706     <NA>      A-S507
## HC_4                   1             TR4-61      A-N702     <NA>      A-S505
## HC_5                   1             TR5-61      A-N706     <NA>      D-S520
## HC_6                   1             TR6-61      A-N715     <NA>      D-S515
## HC_7                   1             TR1-63      A-N701     <NA>      D-S517
## HC_8                   1             TR2-63      A-N710     <NA>      D-S522
## HC_9                   1             TR3-63      D-N718     <NA>      A-S508
## HC_10                  1             TR4-63      A-N712     <NA>      A-S503
## HC_11                  1             TR5-63      A-N714     <NA>      A-S503
## D7                     1             TR1-15      N711-A AAGAGGCA      S522-D
## D8                     1             TR1-17      N723-D TAGCGCTC      S508-A
## D15                    1             TR1-19      N719-D GCGTAGTA      S513-D
## D22                    1             TR1-21      N718-D GGAGCTAC      S510-A
## D29                    1             TR1-46      N712-A GTAGAGGA      S515-D
## HC_12                  1             TR6-63      D-N716     <NA>      A-S505
## D42                    1             TR1-64      N726-D CCTAAGAC      S510-A
## HC_13                  1             TR1-65      A-N704     <NA>      D-S522
## HC_14                  1             TR2-65      A-N710     <NA>      D-S520
## D6                     1             TR2-15      N729-D TCGACGTC      S511-A
## D9                     1             TR2-17      N723-D TAGCGCTC      S511-A
## D16                    1             TR2-19      N720-D CGGAGCCT      S520-D
## D23                    1             TR2-21      N721-D TACGCTGC      S506-A
## HC_15                  1             TR3-65      A-N715     <NA>      A-S506
## HC_16                  1             TR4-65      A-N710     <NA>      D-S516
## HC_17                  1             TR5-65      A-N715     <NA>      D-S516
## HC_18                  1             TR6-65      A-N715     <NA>      A-S508
## HC_19                  1             TR1-69      A-N707     <NA>      D-S513
## HC_20                  1             TR2-69      A-N710     <NA>      A-S502
## D30                    1             TR2-46      N722-D ATGCGCAG      S510-A
## HC_21                  1             TR3-69      A-N706     <NA>      D-S516
## HC_22                  1             TR4-69      D-N716     <NA>      A-S507
## HC_23                  1             TR5-69      A-N702     <NA>      D-S516
## HC_24                  1             TR6-69      A-N712     <NA>      D-S520
## D41                    1             TR2-64      N719-D GCGTAGTA      S511-A
## D5                     1             TR3-15      N723-D TAGCGCTC      S515-D
## D10                    1             TR3-17      N720-D CGGAGCCT      S508-A
## D17                    1             TR3-19      N727-D CGATCAGT      S503-A
## D24                    1             TR3-21      N722-D ATGCGCAG      S503-A
## HV_1                   1             TR1-28      A-N705     <NA>      A-S505
## HV_2                   1             TR2-28      A-N710     <NA>      A-S511
## HV_3                   1             TR3-28      A-N705     <NA>      A-S508
## HV_4                   1             TR4-28      A-N706     <NA>      D-S515
## HV_5                   1             TR5-28      A-N701     <NA>      A-S507
## HV_6                   1             TR6-28      A-N712     <NA>      A-S511
## D31                    1             TR3-46      N720-D CGGAGCCT      S511-A
## D40                    1             TR3-64      N715-A ATCTCAGG      S520-D
## D4                     1             TR4-15      N720-D CGGAGCCT      S503-A
## D11                    1             TR4-17      N703-A AGGCAGAA      S515-D
## D18                    1             TR4-19      N727-D CGATCAGT      S508-A
## D25                    1             TR4-21      N724-D ACTGAGCG      S508-A
## HV_7                   1             TR1-30      A-N712     <NA>      A-S506
## HV_8                   1             TR2-30      A-N712     <NA>      A-S502
## HV_9                   1             TR3-30      A-N715     <NA>      D-S520
## HV_10                  1             TR4-30      A-N706     <NA>      A-S511
## HV_11                  1             TR5-30      A-N710     <NA>      A-S508
## HV_12                  1             TR6-30      A-N706     <NA>      D-S522
## D32                    1             TR4-46      N719-D GCGTAGTA      S510-A
## D39                    1             TR4-64      N720-D CGGAGCCT      S506-A
## C_1                    1             TR5-13      N701-A TAAGGCGA      S521-D
## D3                     1             TR5-15      N726-D CCTAAGAC      S520-D
## D12                    1             TR5-17      N714-A GCTCATGA      S503-A
## D19                    1             TR5-19      N714-A GCTCATGA      S510-A
## C_6                    1             TR5-20      N718-D GGAGCTAC      S511-A
## D26                    1             TR5-21      N715-A ATCTCAGG      S513-D
## D33                    1             TR5-46      N719-D GCGTAGTA      S520-D
## HV_13                  1             TR1-32      A-N702     <NA>      D-S521
## HV_14                  1             TR2-32      A-N715     <NA>      A-S503
## HV_15                  1             TR3-32      D-N719     <NA>      A-S505
## HV_16                  1             TR4-32      A-N703     <NA>      A-S502
## HV_17                  1             TR5-32      A-N707     <NA>      A-S511
## D38                    1             TR5-64      N706-A TAGGCATG      S518-D
## HV_18                  1             TR6-32      A-N705     <NA>      A-S510
## D_1                    1                D-2      A-N711     <NA>      D-S517
## D2                     1             TR6-15      N707-A CTCTCTAC      S518-D
## D13                    1             TR6-17      N706-A TAGGCATG      S515-D
## D20                    1             TR6-19      N706-A TAGGCATG      S513-D
## D27                    1             TR6-21      N702-A CGTACTAG      S516-D
## HV_19                  1             TR1-36      D-N716     <NA>      A-S508
## HV_20                  1             TR2-36      A-N714     <NA>      D-S521
## HV_21                  1             TR3-36      A-N701     <NA>      A-S511
## HV_22                  1             TR4-36      A-N702     <NA>      A-S506
## HV_23                  1             TR5-36      A-N714     <NA>      A-S507
## D34                    1             TR6-46      N712-A GTAGAGGA      S522-D
## D37                    1             TR6-64      N707-A CTCTCTAC      S515-D
## HV_24                  1             TR6-36      A-N704     <NA>      D-S517
## C_2                   NA             TR5-63        <NA>     <NA>        <NA>
## C_3                   NA             TR4-63        <NA>     <NA>        <NA>
## C_4                   NA             TR3-63        <NA>     <NA>        <NA>
## C_5                   NA             TR2-63        <NA>     <NA>        <NA>
##         index2 Description2 Experiment Reactor     Treatment Day_of_Connection
## D1    GCGTAAGA         <NA> Continuous      CR     UNTREATED                15
## D14   TCGACTAG         <NA> Continuous      CR     UNTREATED                17
## D21   CTATTAAG         <NA> Continuous      CR     UNTREATED                19
## D28   AAGGCTAT         <NA> Continuous      CR     UNTREATED                21
## D35   GAGCCTTA         <NA> Continuous      CR     UNTREATED                46
## D36   GCGTAAGA         <NA> Continuous      CR     UNTREATED                64
## D43   TTCTAGCT        DONOR      Cecum   DONOR         DONOR                NA
## HC_1      <NA>         <NA> Continuous     TR1   CTX+HV292.1                 8
## HC_2      <NA>         <NA> Continuous     TR2           CTX                 8
## HC_3      <NA>         <NA> Continuous     TR3   CTX+HV292.1                 8
## HC_4      <NA>         <NA> Continuous     TR4           CTX                 8
## HC_5      <NA>         <NA> Continuous     TR5       HV292.1                 8
## HC_6      <NA>         <NA> Continuous      CR     UNTREATED                 8
## HC_7      <NA>         <NA> Continuous     TR1   CTX+HV292.1                10
## HC_8      <NA>         <NA> Continuous     TR2           CTX                10
## HC_9      <NA>         <NA> Continuous     TR3   CTX+HV292.1                10
## HC_10     <NA>         <NA> Continuous     TR4           CTX                10
## HC_11     <NA>         <NA> Continuous     TR5       HV292.1                10
## D7    TTATGCGA         <NA> Continuous     TR1   CTX+HV292.1                15
## D8    CTAAGCCT         <NA> Continuous     TR1   CTX+HV292.1                17
## D15   TCGACTAG         <NA> Continuous     TR1   CTX+HV292.1                19
## D22   CGTCTAAT         <NA> Continuous     TR1   CTX+HV292.1                21
## D29   TTCTAGCT         <NA> Continuous     TR1   CTX+HV292.1                46
## HC_12     <NA>         <NA> Continuous      CR     UNTREATED                10
## D42   CGTCTAAT         <NA> Continuous     TR1   CTX+HV292.1                64
## HC_13     <NA>         <NA> Continuous     TR1   CTX+HV292.1                12
## HC_14     <NA>         <NA> Continuous     TR2           CTX                12
## D6    TCTCTCCG         <NA> Continuous     TR2           CTX                15
## D9    TCTCTCCG         <NA> Continuous     TR2           CTX                17
## D16   AAGGCTAT         <NA> Continuous     TR2           CTX                19
## D23   ACTGCATA         <NA> Continuous     TR2           CTX                21
## HC_15     <NA>         <NA> Continuous     TR3   CTX+HV292.1                12
## HC_16     <NA>         <NA> Continuous     TR4           CTX                12
## HC_17     <NA>         <NA> Continuous     TR5       HV292.1                12
## HC_18     <NA>         <NA> Continuous      CR     UNTREATED                12
## HC_19     <NA>         <NA> Continuous     TR1   CTX+HV292.1                16
## HC_20     <NA>         <NA> Continuous     TR2           CTX                16
## D30   CGTCTAAT         <NA> Continuous     TR2           CTX                46
## HC_21     <NA>         <NA> Continuous     TR3   CTX+HV292.1                16
## HC_22     <NA>         <NA> Continuous     TR4           CTX                16
## HC_23     <NA>         <NA> Continuous     TR5       HV292.1                16
## HC_24     <NA>         <NA> Continuous      CR     UNTREATED                16
## D41   TCTCTCCG         <NA> Continuous     TR2           CTX                64
## D5    TTCTAGCT         <NA> Continuous     TR3       HV292.1                15
## D10   CTAAGCCT         <NA> Continuous     TR3       HV292.1                17
## D17   TATCCTCT         <NA> Continuous     TR3       HV292.1                19
## D24   TATCCTCT         <NA> Continuous     TR3       HV292.1                21
## HV_1      <NA>         <NA> Continuous     TR1 VAN+CCUG59168                28
## HV_2      <NA>         <NA> Continuous     TR2           VAN                28
## HV_3      <NA>         <NA> Continuous     TR3 VAN+CCUG59168                28
## HV_4      <NA>         <NA> Continuous     TR4           VAN                28
## HV_5      <NA>         <NA> Continuous     TR5     CCUG59168                28
## HV_6      <NA>         <NA> Continuous      CR     UNTREATED                28
## D31   TCTCTCCG         <NA> Continuous     TR3       HV292.1                46
## D40   AAGGCTAT         <NA> Continuous     TR3       HV292.1                64
## D4    TATCCTCT         <NA> Continuous     TR4           VAN                15
## D11   TTCTAGCT         <NA> Continuous     TR4           VAN                17
## D18   CTAAGCCT         <NA> Continuous     TR4           VAN                19
## D25   CTAAGCCT         <NA> Continuous     TR4           VAN                21
## HV_7      <NA>         <NA> Continuous     TR1 VAN+CCUG59168                30
## HV_8      <NA>         <NA> Continuous     TR2           VAN                30
## HV_9      <NA>         <NA> Continuous     TR3 VAN+CCUG59168                30
## HV_10     <NA>         <NA> Continuous     TR4           VAN                30
## HV_11     <NA>         <NA> Continuous     TR5     CCUG59168                30
## HV_12     <NA>         <NA> Continuous      CR     UNTREATED                30
## D32   CGTCTAAT         <NA> Continuous     TR4           VAN                46
## D39   ACTGCATA         <NA> Continuous     TR4           VAN                64
## C_1   GAGCCTTA         <NA> Continuous     TR5 VAN+CCUG59168                13
## D3    AAGGCTAT         <NA> Continuous     TR5 VAN+CCUG59168                15
## D12   TATCCTCT         <NA> Continuous     TR5 VAN+CCUG59168                17
## D19   CGTCTAAT         <NA> Continuous     TR5 VAN+CCUG59168                19
## C_6   TCTCTCCG         <NA> Continuous     TR5 VAN+CCUG59168                20
## D26   TCGACTAG         <NA> Continuous     TR5 VAN+CCUG59168                21
## D33   AAGGCTAT         <NA> Continuous     TR5 VAN+CCUG59168                46
## HV_13     <NA>         <NA> Continuous     TR1 VAN+CCUG59168                32
## HV_14     <NA>         <NA> Continuous     TR2           VAN                32
## HV_15     <NA>         <NA> Continuous     TR3 VAN+CCUG59168                32
## HV_16     <NA>         <NA> Continuous     TR4           VAN                32
## HV_17     <NA>         <NA> Continuous     TR5     CCUG59168                32
## D38   CTATTAAG         <NA> Continuous     TR5 VAN+CCUG59168                64
## HV_18     <NA>         <NA> Continuous      CR     UNTREATED                32
## D_1       <NA>        DONOR      Cecum   DONOR         DONOR                NA
## D2    CTATTAAG         <NA> Continuous     TR6     CCUG59168                15
## D13   TTCTAGCT         <NA> Continuous     TR6     CCUG59168                17
## D20   TCGACTAG         <NA> Continuous     TR6     CCUG59168                19
## D27   CCTAGAGT         <NA> Continuous     TR6     CCUG59168                21
## HV_19     <NA>         <NA> Continuous     TR1 VAN+CCUG59168                36
## HV_20     <NA>         <NA> Continuous     TR2           VAN                36
## HV_21     <NA>         <NA> Continuous     TR3 VAN+CCUG59168                36
## HV_22     <NA>         <NA> Continuous     TR4           VAN                36
## HV_23     <NA>         <NA> Continuous     TR5     CCUG59168                36
## D34   TTATGCGA         <NA> Continuous     TR6     CCUG59168                46
## D37   TTCTAGCT         <NA> Continuous     TR6     CCUG59168                64
## HV_24     <NA>         <NA> Continuous      CR     UNTREATED                36
## C_2       <NA>         <NA> Continuous     TR5 VAN+CCUG59168                63
## C_3       <NA>         <NA> Continuous     TR4           VAN                63
## C_4       <NA>         <NA> Continuous     TR3       HV292.1                63
## C_5       <NA>         <NA> Continuous     TR2           CTX                63
##       Day_of_Treatment Day_from_Inoculum  Enrichment Phase    Treatment2
## D1                  -1                38 NotEnriched  Stab     UNTREATED
## D14                  1                40 NotEnriched Treat     UNTREATED
## D21                  3                42 NotEnriched Treat     UNTREATED
## D28                  5                44 NotEnriched Treat     UNTREATED
## D35                 30                69 NotEnriched Treat     UNTREATED
## D36                 48                87 NotEnriched Treat     UNTREATED
## D43                 NA                NA NotEnriched  <NA>         DONOR
## HC_1                -1                76 NotEnriched  Stab    AB+E. coli
## HC_2                -1                76 NotEnriched  Stab            AB
## HC_3                -1                76 NotEnriched  Stab    AB+E. coli
## HC_4                -1                76 NotEnriched  Stab            AB
## HC_5                -1                76 NotEnriched  Stab       E. coli
## HC_6                -1                76 NotEnriched  Stab     UNTREATED
## HC_7                 1                78 NotEnriched Treat    AB+E. coli
## HC_8                 1                78 NotEnriched Treat            AB
## HC_9                 1                78 NotEnriched Treat    AB+E. coli
## HC_10                1                78 NotEnriched Treat            AB
## HC_11                1                78 NotEnriched Treat       E. coli
## D7                  -1                38 NotEnriched  Stab    AB+E. coli
## D8                   1                40 NotEnriched Treat    AB+E. coli
## D15                  3                42 NotEnriched Treat    AB+E. coli
## D22                  5                44 NotEnriched Treat    AB+E. coli
## D29                 30                69 NotEnriched Treat    AB+E. coli
## HC_12                1                78 NotEnriched Treat     UNTREATED
## D42                 48                87 NotEnriched Treat    AB+E. coli
## HC_13                3                80 NotEnriched Treat    AB+E. coli
## HC_14                3                80 NotEnriched Treat            AB
## D6                  -1                38 NotEnriched  Stab            AB
## D9                   1                40 NotEnriched Treat            AB
## D16                  3                42 NotEnriched Treat            AB
## D23                  5                44 NotEnriched Treat            AB
## HC_15                3                80 NotEnriched Treat    AB+E. coli
## HC_16                3                80 NotEnriched Treat            AB
## HC_17                3                80 NotEnriched Treat       E. coli
## HC_18                3                80 NotEnriched Treat     UNTREATED
## HC_19                7                84 NotEnriched Treat    AB+E. coli
## HC_20                7                84 NotEnriched Treat            AB
## D30                 30                69 NotEnriched Treat            AB
## HC_21                7                84 NotEnriched Treat    AB+E. coli
## HC_22                7                84 NotEnriched Treat            AB
## HC_23                7                84 NotEnriched Treat       E. coli
## HC_24                7                84 NotEnriched Treat     UNTREATED
## D41                 48                87 NotEnriched Treat            AB
## D5                  -1                38 NotEnriched  Stab       E. coli
## D10                  1                40 NotEnriched Treat       E. coli
## D17                  3                42 NotEnriched Treat       E. coli
## D24                  5                44 NotEnriched Treat       E. coli
## HV_1                -1                43 NotEnriched  Stab AB+E. faecium
## HV_2                -1                43 NotEnriched  Stab            AB
## HV_3                -1                43 NotEnriched  Stab AB+E. faecium
## HV_4                -1                43 NotEnriched  Stab            AB
## HV_5                -1                43 NotEnriched  Stab    E. faecium
## HV_6                -1                43 NotEnriched  Stab     UNTREATED
## D31                 30                69 NotEnriched Treat       E. coli
## D40                 48                87 NotEnriched Treat       E. coli
## D4                  -1                38 NotEnriched  Stab            AB
## D11                  1                40 NotEnriched Treat            AB
## D18                  3                42 NotEnriched Treat            AB
## D25                  5                44 NotEnriched Treat            AB
## HV_7                 1                45 NotEnriched Treat AB+E. faecium
## HV_8                 1                45 NotEnriched Treat            AB
## HV_9                 1                45 NotEnriched Treat AB+E. faecium
## HV_10                1                45 NotEnriched Treat            AB
## HV_11                1                45 NotEnriched Treat    E. faecium
## HV_12                1                45 NotEnriched Treat     UNTREATED
## D32                 30                69 NotEnriched Treat            AB
## D39                 48                87 NotEnriched Treat            AB
## C_1                 -3                36 NotEnriched  Stab AB+E. faecium
## D3                  -1                38 NotEnriched  Stab AB+E. faecium
## D12                  1                40 NotEnriched Treat AB+E. faecium
## D19                  3                42 NotEnriched Treat AB+E. faecium
## C_6                  4                43 NotEnriched Treat AB+E. faecium
## D26                  5                44 NotEnriched Treat AB+E. faecium
## D33                 30                69 NotEnriched Treat AB+E. faecium
## HV_13                3                47 NotEnriched Treat AB+E. faecium
## HV_14                3                47 NotEnriched Treat            AB
## HV_15                3                47 NotEnriched Treat AB+E. faecium
## HV_16                3                47 NotEnriched Treat            AB
## HV_17                3                47 NotEnriched Treat    E. faecium
## D38                 48                87 NotEnriched Treat AB+E. faecium
## HV_18                3                47 NotEnriched Treat     UNTREATED
## D_1                 NA                NA NotEnriched  <NA>         DONOR
## D2                  -1                38 NotEnriched  Stab    E. faecium
## D13                  1                40 NotEnriched Treat    E. faecium
## D20                  3                42 NotEnriched Treat    E. faecium
## D27                  5                44 NotEnriched Treat    E. faecium
## HV_19                7                51 NotEnriched Treat AB+E. faecium
## HV_20                7                51 NotEnriched Treat            AB
## HV_21                7                51 NotEnriched Treat AB+E. faecium
## HV_22                7                51 NotEnriched Treat            AB
## HV_23                7                51 NotEnriched Treat    E. faecium
## D34                 30                69 NotEnriched Treat    E. faecium
## D37                 48                87 NotEnriched Treat    E. faecium
## HV_24                7                51 NotEnriched Treat     UNTREATED
## C_2                 47                86 NotEnriched Treat AB+E. faecium
## C_3                 47                86 NotEnriched Treat            AB
## C_4                 47                86 NotEnriched Treat       E. coli
## C_5                 47                86 NotEnriched Treat            AB
##                       Date Paul Reactor_Treatment GeneCopyNumberperML
## D1    2020-06-26T00:00:00Z <NA>      CR_UNTREATED            6.75e+10
## D14   2020-06-28T00:00:00Z <NA>      CR_UNTREATED            7.01e+10
## D21   2020-06-30T00:00:00Z <NA>      CR_UNTREATED            9.01e+10
## D28   2020-07-02T00:00:00Z <NA>      CR_UNTREATED            4.39e+10
## D35   2020-07-27T00:00:00Z <NA>      CR_UNTREATED            4.53e+10
## D36   2020-08-14T00:00:00Z <NA>      CR_UNTREATED            3.26e+10
## D43   2020-05-19T00:00:00Z <NA>             DONOR            1.44e+10
## HC_1                  <NA> <NA>   TR1_CTX+HV292.1            2.14e+11
## HC_2                  <NA> <NA>           TR2_CTX            2.07e+11
## HC_3                  <NA> <NA>   TR3_CTX+HV292.1            1.98e+11
## HC_4                  <NA> <NA>           TR4_CTX            2.05e+11
## HC_5                  <NA> <NA>       TR5_HV292.1            2.43e+11
## HC_6                  <NA> <NA>      CR_UNTREATED            2.37e+11
## HC_7                  <NA> <NA>   TR1_CTX+HV292.1            5.59e+11
## HC_8                  <NA> <NA>           TR2_CTX            2.66e+11
## HC_9                  <NA> <NA>   TR3_CTX+HV292.1            3.83e+11
## HC_10                 <NA> <NA>           TR4_CTX            3.22e+11
## HC_11                 <NA> <NA>       TR5_HV292.1            4.52e+11
## D7    2020-06-26T00:00:00Z <NA>   TR1_CTX+HV292.1            5.08e+10
## D8    2020-06-28T00:00:00Z <NA>   TR1_CTX+HV292.1            3.11e+10
## D15   2020-06-30T00:00:00Z <NA>   TR1_CTX+HV292.1            3.36e+10
## D22   2020-07-02T00:00:00Z <NA>   TR1_CTX+HV292.1            2.22e+10
## D29   2020-07-27T00:00:00Z <NA>   TR1_CTX+HV292.1            6.05e+10
## HC_12                 <NA> <NA>      CR_UNTREATED            2.39e+11
## D42   2020-08-14T00:00:00Z <NA>   TR1_CTX+HV292.1            3.61e+10
## HC_13                 <NA> <NA>   TR1_CTX+HV292.1            1.36e+11
## HC_14                 <NA> <NA>           TR2_CTX            2.85e+11
## D6    2020-06-26T00:00:00Z <NA>           TR2_CTX            6.23e+10
## D9    2020-06-28T00:00:00Z <NA>           TR2_CTX            2.71e+10
## D16   2020-06-30T00:00:00Z <NA>           TR2_CTX            4.20e+10
## D23   2020-07-02T00:00:00Z <NA>           TR2_CTX            2.89e+10
## HC_15                 <NA> <NA>   TR3_CTX+HV292.1            2.04e+11
## HC_16                 <NA> <NA>           TR4_CTX            2.39e+11
## HC_17                 <NA> <NA>       TR5_HV292.1            2.54e+11
## HC_18                 <NA> <NA>      CR_UNTREATED            1.26e+11
## HC_19                 <NA> <NA>   TR1_CTX+HV292.1            2.10e+11
## HC_20                 <NA> <NA>           TR2_CTX            1.02e+11
## D30   2020-07-27T00:00:00Z <NA>           TR2_CTX            5.05e+10
## HC_21                 <NA> <NA>   TR3_CTX+HV292.1            2.41e+11
## HC_22                 <NA> <NA>           TR4_CTX            2.42e+11
## HC_23                 <NA> <NA>       TR5_HV292.1            3.93e+11
## HC_24                 <NA> <NA>      CR_UNTREATED            2.60e+11
## D41   2020-08-14T00:00:00Z <NA>           TR2_CTX            2.79e+10
## D5    2020-06-26T00:00:00Z <NA>       TR3_HV292.1            4.87e+10
## D10   2020-06-28T00:00:00Z <NA>       TR3_HV292.1            5.24e+10
## D17   2020-06-30T00:00:00Z <NA>       TR3_HV292.1            4.62e+10
## D24   2020-07-02T00:00:00Z <NA>       TR3_HV292.1            4.79e+10
## HV_1                  <NA> <NA> TR1_VAN+CCUG59168            1.11e+11
## HV_2                  <NA> <NA>           TR2_VAN            2.50e+11
## HV_3                  <NA> <NA> TR3_VAN+CCUG59168            3.13e+11
## HV_4                  <NA> <NA>           TR4_VAN            2.99e+10
## HV_5                  <NA> <NA>     TR5_CCUG59168            1.92e+11
## HV_6                  <NA> <NA>      CR_UNTREATED            2.58e+11
## D31   2020-07-27T00:00:00Z <NA>       TR3_HV292.1            3.84e+10
## D40   2020-08-14T00:00:00Z <NA>       TR3_HV292.1            3.43e+10
## D4    2020-06-26T00:00:00Z <NA>           TR4_VAN            6.68e+10
## D11   2020-06-28T00:00:00Z <NA>           TR4_VAN            3.74e+10
## D18   2020-06-30T00:00:00Z <NA>           TR4_VAN            5.61e+10
## D25   2020-07-02T00:00:00Z <NA>           TR4_VAN            5.82e+10
## HV_7                  <NA> <NA> TR1_VAN+CCUG59168            3.03e+11
## HV_8                  <NA> <NA>           TR2_VAN            7.45e+10
## HV_9                  <NA> <NA> TR3_VAN+CCUG59168            1.21e+11
## HV_10                 <NA> <NA>           TR4_VAN            1.32e+11
## HV_11                 <NA> <NA>     TR5_CCUG59168            2.77e+11
## HV_12                 <NA> <NA>      CR_UNTREATED            4.92e+10
## D32   2020-07-27T00:00:00Z <NA>           TR4_VAN            6.31e+10
## D39   2020-08-14T00:00:00Z <NA>           TR4_VAN            6.94e+10
## C_1   2020-06-24T00:00:00Z <NA> TR5_VAN+CCUG59168            3.94e+10
## D3    2020-06-26T00:00:00Z <NA> TR5_VAN+CCUG59168            7.15e+10
## D12   2020-06-28T00:00:00Z <NA> TR5_VAN+CCUG59168            2.60e+10
## D19   2020-06-30T00:00:00Z <NA> TR5_VAN+CCUG59168            5.71e+10
## C_6   2020-07-01T00:00:00Z <NA> TR5_VAN+CCUG59168            2.50e+10
## D26   2020-07-02T00:00:00Z <NA> TR5_VAN+CCUG59168            4.89e+10
## D33   2020-07-27T00:00:00Z <NA> TR5_VAN+CCUG59168            6.59e+10
## HV_13                 <NA> <NA> TR1_VAN+CCUG59168            1.94e+11
## HV_14                 <NA> <NA>           TR2_VAN            2.43e+11
## HV_15                 <NA> <NA> TR3_VAN+CCUG59168            2.00e+11
## HV_16                 <NA> <NA>           TR4_VAN            9.63e+10
## HV_17                 <NA> <NA>     TR5_CCUG59168            2.00e+11
## D38   2020-08-14T00:00:00Z <NA> TR5_VAN+CCUG59168            5.32e+10
## HV_18                 <NA> <NA>      CR_UNTREATED            3.02e+11
## D_1                   <NA> <NA>             DONOR            1.90e+10
## D2    2020-06-26T00:00:00Z <NA>     TR6_CCUG59168            5.30e+10
## D13   2020-06-28T00:00:00Z <NA>     TR6_CCUG59168            3.93e+10
## D20   2020-06-30T00:00:00Z <NA>     TR6_CCUG59168            5.68e+10
## D27   2020-07-02T00:00:00Z <NA>     TR6_CCUG59168            6.93e+10
## HV_19                 <NA> <NA> TR1_VAN+CCUG59168            4.92e+10
## HV_20                 <NA> <NA>           TR2_VAN            3.21e+11
## HV_21                 <NA> <NA> TR3_VAN+CCUG59168            2.91e+11
## HV_22                 <NA> <NA>           TR4_VAN            1.17e+11
## HV_23                 <NA> <NA>     TR5_CCUG59168            1.68e+11
## D34   2020-07-27T00:00:00Z <NA>     TR6_CCUG59168            2.81e+10
## D37   2020-08-14T00:00:00Z <NA>     TR6_CCUG59168            3.93e+10
## HV_24                 <NA> <NA>      CR_UNTREATED            4.39e+11
## C_2               13.08.20 <NA> TR5_VAN+CCUG59168                  NA
## C_3               13.08.20 <NA>           TR4_VAN                  NA
## C_4               13.08.20 <NA>       TR3_HV292.1                  NA
## C_5               13.08.20 <NA>           TR2_CTX                  NA
##       HV292.1_Copy_Number_permL CCUG59168_Copy_Number_permL
## D1                           NA                          NA
## D14                          NA                          NA
## D21                          NA                          NA
## D28                          NA                          NA
## D35                          NA                          NA
## D36                          NA                          NA
## D43                          NA                          NA
## HC_1                         NA                          NA
## HC_2                         NA                          NA
## HC_3                         NA                          NA
## HC_4                         NA                          NA
## HC_5                         NA                          NA
## HC_6                         NA                          NA
## HC_7                    5520000                          NA
## HC_8                          0                          NA
## HC_9                   30487500                          NA
## HC_10                         0                          NA
## HC_11                  20600000                          NA
## D7                           NA                          NA
## D8                           NA                          NA
## D15                      794000                          NA
## D22                          NA                          NA
## D29                    54054532                          NA
## HC_12                         0                          NA
## D42                   862908750                          NA
## HC_13                         0                          NA
## HC_14                         0                          NA
## D6                           NA                          NA
## D9                           NA                          NA
## D16                          NA                          NA
## D23                          NA                          NA
## HC_15                         0                          NA
## HC_16                         0                          NA
## HC_17                         0                          NA
## HC_18                         0                          NA
## HC_19                         0                          NA
## HC_20                         0                          NA
## D30                          NA                          NA
## HC_21                         0                          NA
## HC_22                         0                          NA
## HC_23                         0                          NA
## HC_24                         0                          NA
## D41                          NA                          NA
## D5                           NA                          NA
## D10                          NA                          NA
## D17                      887000                          NA
## D24                          NA                          NA
## HV_1                         NA                          NA
## HV_2                         NA                          NA
## HV_3                         NA                          NA
## HV_4                         NA                          NA
## HV_5                         NA                          NA
## HV_6                         NA                          NA
## D31                           0                          NA
## D40                           0                          NA
## D4                           NA                          NA
## D11                          NA                          NA
## D18                          NA                          NA
## D25                          NA                          NA
## HV_7                         NA                1.380000e+10
## HV_8                         NA                0.000000e+00
## HV_9                         NA                2.080000e+09
## HV_10                        NA                0.000000e+00
## HV_11                        NA                1.570000e+07
## HV_12                        NA                0.000000e+00
## D32                          NA                          NA
## D39                          NA                          NA
## C_1                          NA                          NA
## D3                           NA                          NA
## D12                          NA                          NA
## D19                          NA                5.040000e+07
## C_6                          NA                          NA
## D26                          NA                          NA
## D33                          NA                9.198400e+07
## HV_13                        NA                8.450000e+09
## HV_14                        NA                0.000000e+00
## HV_15                        NA                2.690000e+09
## HV_16                        NA                0.000000e+00
## HV_17                        NA                8.830000e+05
## D38                          NA                5.390000e+08
## HV_18                        NA                0.000000e+00
## D_1                          NA                          NA
## D2                           NA                          NA
## D13                          NA                          NA
## D20                          NA                0.000000e+00
## D27                          NA                          NA
## HV_19                        NA                5.165000e+09
## HV_20                        NA                0.000000e+00
## HV_21                        NA                1.670000e+09
## HV_22                        NA                0.000000e+00
## HV_23                        NA                1.290000e+05
## D34                          NA                2.688168e+05
## D37                          NA                2.120000e+03
## HV_24                        NA                0.000000e+00
## C_2                          NA                          NA
## C_3                          NA                          NA
## C_4                          NA                          NA
## C_5                          NA                          NA
##       CTX_Copy_Number_permL VAN_Copy_Number_permL   Model Antibiotic_mg.mL
## D1                       NA                    NA Chicken               NA
## D14                      NA                    NA Chicken               NA
## D21                      NA                    NA Chicken               NA
## D28                      NA                    NA Chicken               NA
## D35                      NA                    NA Chicken               NA
## D36                      NA                    NA Chicken               NA
## D43                      NA                    NA Chicken               NA
## HC_1                     NA                    NA   Human               20
## HC_2                     NA                    NA   Human               20
## HC_3                     NA                    NA   Human              200
## HC_4                     NA                    NA   Human              200
## HC_5                     NA                    NA   Human               NA
## HC_6                     NA                    NA   Human               NA
## HC_7           1.102500e+07                    NA   Human               20
## HC_8           0.000000e+00                    NA   Human               20
## HC_9           5.373750e+07                    NA   Human              200
## HC_10          0.000000e+00                    NA   Human              200
## HC_11          4.160000e+07                    NA   Human               NA
## D7                       NA                    NA Chicken               20
## D8                       NA                    NA Chicken               20
## D15            1.507417e+06                    NA Chicken               20
## D22                      NA                    NA Chicken               20
## D29            8.070878e+07                    NA Chicken               20
## HC_12          0.000000e+00                    NA   Human               NA
## D42            1.927668e+09                    NA Chicken               20
## HC_13          2.160000e+06                    NA   Human               20
## HC_14          0.000000e+00                    NA   Human               20
## D6                       NA                    NA Chicken               20
## D9                       NA                    NA Chicken               20
## D16                      NA                    NA Chicken               20
## D23                      NA                    NA Chicken               20
## HC_15          3.550000e+06                    NA   Human              200
## HC_16          0.000000e+00                    NA   Human              200
## HC_17          2.100000e+05                    NA   Human               NA
## HC_18          0.000000e+00                    NA   Human               NA
## HC_19          0.000000e+00                    NA   Human               20
## HC_20          0.000000e+00                    NA   Human               20
## D30                      NA                    NA Chicken               20
## HC_21          1.117500e+06                    NA   Human              200
## HC_22          0.000000e+00                    NA   Human              200
## HC_23                    NA                    NA   Human               NA
## HC_24          0.000000e+00                    NA   Human               NA
## D41                      NA                    NA Chicken               20
## D5                       NA                    NA Chicken               NA
## D10                      NA                    NA Chicken               NA
## D17            1.670000e+06                    NA Chicken               NA
## D24                      NA                    NA Chicken               NA
## HV_1                     NA                    NA   Human               90
## HV_2                     NA                    NA   Human               90
## HV_3                     NA                    NA   Human              600
## HV_4                     NA                    NA   Human              600
## HV_5                     NA                    NA   Human               NA
## HV_6                     NA                    NA   Human               NA
## D31            4.713665e+03                    NA Chicken               NA
## D40            3.666540e+04                    NA Chicken               NA
## D4                       NA                    NA Chicken               90
## D11                      NA                    NA Chicken               90
## D18                      NA                    NA Chicken               90
## D25                      NA                    NA Chicken               90
## HV_7                     NA          9.600000e+12   Human               90
## HV_8                     NA          0.000000e+00   Human               90
## HV_9                     NA          8.500000e+11   Human              600
## HV_10                    NA          0.000000e+00   Human              600
## HV_11                    NA          6.125000e+09   Human               NA
## HV_12                    NA          0.000000e+00   Human               NA
## D32                      NA                    NA Chicken               90
## D39                      NA                    NA Chicken               90
## C_1                      NA                    NA Chicken               90
## D3                       NA          1.740250e+10 Chicken               90
## D12                      NA          4.156488e+10 Chicken               90
## D19                      NA          6.682795e+10 Chicken               90
## C_6                      NA          2.847335e+10 Chicken               90
## D26                      NA          5.769750e+10 Chicken               90
## D33                      NA          5.770000e+11 Chicken               90
## HV_13                    NA          8.400000e+12   Human               90
## HV_14                    NA          0.000000e+00   Human               90
## HV_15                    NA          1.190000e+12   Human              600
## HV_16                    NA          0.000000e+00   Human              600
## HV_17                    NA          3.470000e+08   Human               NA
## D38                      NA                    NA Chicken               90
## HV_18                    NA          0.000000e+00   Human               NA
## D_1                      NA                    NA   Human               NA
## D2                       NA          1.065350e+08 Chicken               NA
## D13                      NA          0.000000e+00 Chicken               NA
## D20                      NA          0.000000e+00 Chicken               NA
## D27                      NA          0.000000e+00 Chicken               NA
## HV_19                    NA          1.300000e+12   Human               90
## HV_20                    NA          0.000000e+00   Human               90
## HV_21                    NA          4.800000e+11   Human              600
## HV_22                    NA          0.000000e+00   Human              600
## HV_23                    NA          0.000000e+00   Human               NA
## D34                      NA                    NA Chicken               NA
## D37                      NA                    NA Chicken               NA
## HV_24                    NA          0.000000e+00   Human               NA
## C_2                      NA                    NA Chicken               90
## C_3                      NA                    NA Chicken               90
## C_4                      NA                    NA Chicken               NA
## C_5                      NA                    NA Chicken               20
##       Fermentation Antibiotic Lactose_mM Glucose_mM Galactose_mM Succinat_mM
## D1            <NA>       <NA>      0.000      0.947        0.000       1.969
## D14           <NA>       <NA>      0.000      0.000        1.759       4.295
## D21           <NA>       <NA>      0.261      0.000        0.000       0.818
## D28           <NA>       <NA>      0.000      0.000        0.000       0.701
## D35           <NA>       <NA>         NA         NA           NA          NA
## D36           <NA>       <NA>      0.000      0.000        0.000       0.546
## D43           <NA>       <NA>         NA         NA           NA          NA
## HC_1             2        CTX      0.000      0.000        0.536       1.519
## HC_2             2        CTX      0.000      0.000        0.588       2.608
## HC_3             2        CTX      0.000      0.000        0.655       1.531
## HC_4             2        CTX      0.000      0.000        0.622       1.513
## HC_5             2        CTX      0.000      0.000        0.602       1.448
## HC_6             2        CTX      0.000      0.000        0.703       1.356
## HC_7             2        CTX      0.000      0.000        0.314       1.558
## HC_8             2        CTX      0.000      0.000        0.703       1.497
## HC_9             2        CTX      0.000      0.305        0.989       1.200
## HC_10            2        CTX      0.000      0.309        0.958       1.300
## HC_11            2        CTX      0.000      0.000        0.680       1.435
## D7            <NA>        CTX      0.000      0.000        0.000       1.090
## D8            <NA>        CTX      0.000      0.000        0.000       1.201
## D15           <NA>        CTX      0.000      0.000        0.000       2.779
## D22           <NA>        CTX      0.000      0.000        0.000       3.069
## D29           <NA>        CTX      0.000      0.000        0.000       0.527
## HC_12            2        CTX      0.000      0.000        0.464       1.390
## D42           <NA>        CTX      0.000      0.000        0.000       0.455
## HC_13            2        CTX      0.000      0.000        0.367       1.370
## HC_14            2        CTX      0.000      0.000        0.000       1.085
## D6            <NA>        CTX      0.000      0.000        0.000       0.792
## D9            <NA>        CTX      0.000      0.000        0.000       1.116
## D16           <NA>        CTX      0.000      0.000        0.000       2.973
## D23           <NA>        CTX      0.000      0.000        0.000       3.420
## HC_15            2        CTX      0.000      0.000        0.391       1.100
## HC_16            2        CTX      0.000      0.000        0.397       1.051
## HC_17            2        CTX      0.000      0.000        0.000       1.110
## HC_18            2        CTX      0.000      0.000        0.000       1.258
## HC_19            2        CTX      0.000      0.000        0.229       1.194
## HC_20            2        CTX      0.000      0.000        0.759       1.121
## D30           <NA>        CTX      0.000      0.000        0.000       3.776
## HC_21            2        CTX      0.000      0.000        0.000       1.739
## HC_22            2        CTX      0.000      0.000        0.000       1.796
## HC_23            2        CTX      0.000      0.000        0.000       1.730
## HC_24            2        CTX      0.000      0.000        0.000       1.784
## D41           <NA>        CTX      0.000      0.000        0.000       0.511
## D5            <NA>       <NA>      0.000      0.000        0.000       0.650
## D10           <NA>       <NA>      0.000      0.000        0.000       0.751
## D17           <NA>       <NA>      0.000      0.000        0.000       0.629
## D24           <NA>       <NA>      0.000      0.000        0.000       0.639
## HV_1             1        VAN      0.053      0.173        0.000       0.901
## HV_2             1        VAN      0.000      0.228        0.000       0.886
## HV_3             1        VAN      0.069      0.000        0.000       0.844
## HV_4             1        VAN      0.000      0.000        0.000       0.932
## HV_5             1        VAN      0.037      0.284        0.000       1.135
## HV_6             1        VAN      0.046      0.000        0.000       0.819
## D31           <NA>       <NA>      0.000      0.000        0.000       0.453
## D40           <NA>       <NA>      0.000      0.000        0.000       0.456
## D4            <NA>        VAN      0.000      0.000        0.000       0.621
## D11           <NA>        VAN      0.000      1.401        1.045       4.792
## D18           <NA>        VAN      0.000      0.000        1.380      11.605
## D25           <NA>        VAN      0.598      0.000        1.033      10.796
## HV_7             1        VAN      0.000      3.207        4.096       1.351
## HV_8             1        VAN      0.000      9.781        4.362       1.558
## HV_9             1        VAN      0.000      0.000        2.502       1.962
## HV_10            1        VAN      0.000      6.682        4.270       1.587
## HV_11            1        VAN      0.000      0.000        0.000       1.239
## HV_12            1        VAN      0.000      0.000        0.000       1.246
## D32           <NA>        VAN      0.000      0.000        0.000       7.444
## D39           <NA>        VAN      0.000      0.000        0.000       0.418
## C_1           <NA>        VAN      0.000      0.000        0.000       1.525
## D3            <NA>        VAN      0.000      0.000        0.000       0.694
## D12           <NA>        VAN      0.000      1.584        1.135       6.021
## D19           <NA>        VAN      0.489      0.000        1.448       9.526
## C_6           <NA>        VAN      0.550      0.437        0.000       7.267
## D26           <NA>        VAN      0.564      0.000        0.000       9.386
## D33           <NA>        VAN      0.000      0.000        0.000       7.771
## HV_13            1        VAN      1.580      0.767        2.044       1.473
## HV_14            1        VAN      1.238      1.616        1.323       1.637
## HV_15            1        VAN      0.904      0.000        2.563       1.559
## HV_16            1        VAN      1.663      0.000        2.668       1.479
## HV_17            1        VAN      0.000      0.000        0.000       1.233
## D38           <NA>        VAN      0.000      0.000        0.000       3.990
## HV_18            1        VAN      0.000      0.000        0.000       1.251
## D_1           <NA>       <NA>         NA         NA           NA          NA
## D2            <NA>       <NA>      0.000      0.000        0.000       0.734
## D13           <NA>       <NA>      0.000      0.000        0.000       0.726
## D20           <NA>       <NA>      0.000      0.000        0.000       0.748
## D27           <NA>       <NA>      0.000      0.000        0.000       0.729
## HV_19            1        VAN      0.777      0.244        1.325       1.417
## HV_20            1        VAN      0.731      0.000        1.144       1.378
## HV_21            1        VAN      1.220      0.000        0.784       1.435
## HV_22            1        VAN      1.579      0.000        0.750       1.401
## HV_23            1        VAN      0.000      0.000        0.266       1.208
## D34           <NA>       <NA>      0.000      0.000        0.000       1.150
## D37           <NA>       <NA>      0.000      0.000        0.000       0.429
## HV_24            1        VAN      0.000      0.000        0.000       1.187
## C_2           <NA>        VAN         NA         NA           NA          NA
## C_3           <NA>        VAN         NA         NA           NA          NA
## C_4           <NA>        CTX         NA         NA           NA          NA
## C_5           <NA>       <NA>         NA         NA           NA          NA
##       Lactat_mM Formiat_mM Acetat_mM Propionat_mM Isobutyrat_mM Butyrat_mM
## D1        0.000      3.179    49.775       11.560         5.270     22.600
## D14       8.223      5.614    39.530        7.362         0.000      8.315
## D21       4.327      0.000    81.797       10.133         0.000     14.599
## D28       0.000      0.000    86.705       15.767         6.942     27.682
## D35          NA         NA        NA           NA            NA         NA
## D36       0.000      0.000    75.009       11.518         7.914     36.797
## D43          NA         NA        NA           NA            NA         NA
## HC_1      1.338      5.898    68.571       27.257         4.937     48.394
## HC_2      1.287      8.273    71.535       26.884         4.776     45.421
## HC_3      1.310      7.826    76.633       27.355         4.521     44.309
## HC_4      1.238      4.607    70.801       25.104         4.911     50.198
## HC_5      1.321      7.221    73.189       36.735         3.001     33.900
## HC_6      1.423      6.459    52.167       24.281         3.506     53.009
## HC_7      1.613      6.574    75.484       32.194         4.251     41.199
## HC_8      1.750      9.178    71.548       29.822         4.000     43.675
## HC_9      3.095      6.354    69.382       34.611         3.285     36.085
## HC_10     3.283      0.000    64.353       34.868         3.117     39.773
## HC_11     1.293      6.330    70.102       38.008         3.627     35.459
## D7        0.000      0.000    83.366       13.060         8.268     41.333
## D8        0.000      0.000    82.428       12.110         7.901     46.439
## D15       0.000      0.000    77.988        8.582         4.400     49.923
## D22       0.000      0.000    72.326       12.140         0.000     49.287
## D29       8.115      0.000   100.476       13.971         7.021     59.817
## HC_12     1.453      5.867    52.980       25.113         3.912     52.049
## D42       1.035      0.000    69.191       13.707         7.787     36.459
## HC_13     1.552      3.816    70.605       26.932         4.857     39.501
## HC_14     1.505      5.782    68.421       23.025         4.640     42.896
## D6        0.000      0.000    98.070       13.180         8.776     44.333
## D9        0.000      0.000    95.235       12.018         7.794     47.890
## D16       0.000      0.000    90.991        7.423         0.000     44.651
## D23       0.000      0.000    84.906        5.716         0.000     41.240
## HC_15     2.160      5.931    65.980       27.327         3.912     40.000
## HC_16     2.611      0.000    64.822       36.192         3.756     38.253
## HC_17     1.029      5.727    67.573       39.064         3.187     35.241
## HC_18     0.905      5.953    54.688       26.113         3.488     51.621
## HC_19     1.038      3.661    69.252       25.772         5.347     47.850
## HC_20     2.450      0.000    59.247       25.199         5.255     49.825
## D30       0.000      0.000   110.663       12.096         5.428     44.934
## HC_21     1.127      5.850    66.606       24.514         5.216     48.333
## HC_22     1.322      0.000    58.881       33.442         5.939     48.455
## HC_23     1.331     10.916    67.988       35.881         0.000     38.423
## HC_24     0.853      4.163    61.137       27.614         4.315     48.676
## D41       0.000      2.219    57.910       13.311         7.440     35.756
## D5        0.000      0.000    99.568       14.760         7.584     33.721
## D10       0.000      0.000    86.398       11.689         8.204     35.082
## D17       0.000      0.000    85.603       13.331         6.860     29.604
## D24       0.000      0.000    93.635       13.645         6.714     26.799
## HV_1      0.534      5.278    77.645       23.407         3.886     40.871
## HV_2      0.734      6.244    78.499       27.767         4.418     37.598
## HV_3      0.679      7.344    57.376       36.797         4.541     33.022
## HV_4      0.449      7.921    69.627       22.762         3.449     43.182
## HV_5      0.364      3.277    76.749       19.707         4.269     42.214
## HV_6      0.000      8.203    64.866       17.193         2.603     44.874
## D31       0.000      0.000   105.354       14.368         7.536     36.305
## D40       0.000      0.000    83.474       12.933         6.417     23.882
## D4        0.000      0.000   101.441       14.850         7.813     33.872
## D11       0.000      0.000    79.243       17.005         8.207     28.479
## D18       0.000      0.000    48.991       21.902         4.715      6.737
## D25       0.000      9.618    47.637       24.239         0.000      1.522
## HV_7     29.016     10.029    23.392       19.437         2.571     14.661
## HV_8      5.440      2.891    21.269       19.543         3.698     16.553
## HV_9     31.934     10.365    19.888       17.952         3.586     18.726
## HV_10     5.702      2.731    29.523       14.073         2.032     12.535
## HV_11     1.215      5.294    78.115       17.219         4.725     46.529
## HV_12     0.994      9.137    71.280       19.527         2.548     44.810
## D32       8.524      0.000    46.362       19.358         4.941      0.287
## D39       3.808      0.000    41.394       17.022         5.820      2.230
## C_1       0.000      0.000    92.696       11.958         9.009     46.755
## D3        0.000      0.000    94.605       12.469         7.545     42.811
## D12       0.000      0.000    77.564       16.918         7.590     29.690
## D19       0.000      0.000    47.530       23.120         0.000      7.136
## C_6       0.000      0.000    44.220       21.102         0.000      3.984
## D26       0.000      9.711    45.669       24.147         0.000      1.788
## D33       0.000      0.000    41.884       18.412         5.257      1.018
## HV_13     9.844     13.503    11.513        7.533         3.121      9.105
## HV_14     5.033     12.418    20.830       25.674         3.302     12.666
## HV_15    15.308      8.918    11.689       13.911         3.508      7.647
## HV_16     7.130     14.424    12.683        6.561         2.260      5.816
## HV_17     1.269      6.174    73.792       17.652         5.184     48.061
## D38       0.000      0.000    31.082       13.931         6.202      1.974
## HV_18     1.113     10.938    74.408       23.448         2.422     42.201
## D_1          NA         NA        NA           NA            NA         NA
## D2        0.000      0.000    97.498       11.455         7.408     47.083
## D13       0.000      0.000    87.252       11.533         8.115     38.443
## D20       0.000      0.000    87.310       13.946         0.000     31.669
## D27       0.000      0.000    94.750       13.997         7.214     26.634
## HV_19     5.496      8.836    18.665       41.793         4.202     13.490
## HV_20     4.138      8.562    19.052       45.481         4.842     12.833
## HV_21    13.537     10.322    10.928       11.060         2.549      5.173
## HV_22     4.294     12.778    15.865       10.169         2.208      4.905
## HV_23     1.372      6.666    72.238       14.474         4.565     48.449
## D34       0.000      0.000   100.967       15.176         7.531     37.422
## D37       0.000      0.000    85.485       14.118         6.381     23.281
## HV_24     1.171      9.060    64.405       19.489         2.321     42.328
## C_2          NA         NA        NA           NA            NA         NA
## C_3          NA         NA        NA           NA            NA         NA
## C_4          NA         NA        NA           NA            NA         NA
## C_5          NA         NA        NA           NA            NA         NA
##       Isovalerat_mM Valerat_mM Total_SCFA_mM raw_metagenomic_pairs Period
## D1            4.267      3.218       101.838              57121568   <NA>
## D14           1.417      0.365        76.880              44863798     t1
## D21           0.881      0.000       112.555              56304446     t1
## D28           4.870      0.000       142.667              50088120     t1
## D35              NA         NA            NA              55285758     t4
## D36           4.585      6.461       142.830              48021514     t5
## D43              NA         NA            NA              68282190   <NA>
## HC_1          3.842      3.761       166.053                    NA   <NA>
## HC_2          3.394      3.482       168.248                    NA   <NA>
## HC_3          3.614      1.935       169.689                    NA   <NA>
## HC_4          4.062      1.715       164.771                    NA   <NA>
## HC_5          2.809      1.264       161.490                    NA   <NA>
## HC_6          3.828      1.687       148.419                    NA   <NA>
## HC_7          3.894      4.031       171.112                    NA     t1
## HC_8          2.980      3.486       168.639                    NA     t1
## HC_9          3.394      1.704       160.404                    NA     t1
## HC_10         3.610      1.147       152.718                    NA     t1
## HC_11         3.258      1.282       161.474                    NA     t1
## D7            6.524      7.962       161.603              52833332   <NA>
## D8            6.534      4.832       161.445              53009028     t1
## D15           3.876      1.147       148.695              47492202     t1
## D22           7.405      0.522       144.749              55638128     t1
## D29           8.804      9.447       208.178              68334190     t4
## HC_12         3.952      1.797       148.977                    NA     t1
## D42           7.614      8.207       144.455              47906582     t5
## HC_13         4.427      5.098       158.525                    NA     t1
## HC_14         3.259      4.630       155.243                    NA     t1
## D6            6.422      7.313       178.886              51668602   <NA>
## D9            6.038      3.410       173.501              47401086     t1
## D16           2.036      0.036       148.110              44839890     t1
## D23           1.044      0.000       136.326              56766854     t1
## HC_15         5.518      1.598       153.917                    NA     t1
## HC_16         5.307      1.038       153.427                    NA     t1
## HC_17         3.844      1.038       157.813                    NA     t1
## HC_18         4.236      1.783       150.045                    NA     t1
## HC_19         5.290      6.032       165.665                    NA     t1
## HC_20         6.343      4.494       154.693                    NA     t1
## D30           6.609      0.000       183.506              49315184     t4
## HC_21         4.836      6.599       164.820                    NA     t1
## HC_22         7.871      1.967       159.673                    NA     t1
## HC_23         4.747      1.182       162.198                    NA     t1
## HC_24         5.126      1.616       155.284                    NA     t1
## D41           6.925      0.233       124.305              59842344     t5
## D5            6.467      8.363       171.113              55615446   <NA>
## D10           6.429      6.960       155.513              46820988     t1
## D17           6.109      0.000       142.136              49791554     t1
## D24           5.297      0.000       146.729              48166738     t1
## HV_1          4.754      1.239       158.741                    NA   <NA>
## HV_2          5.287      1.303       162.964                    NA   <NA>
## HV_3          5.700      0.000       146.372                    NA   <NA>
## HV_4          2.876      0.475       151.673                    NA   <NA>
## HV_5          4.334      0.284       152.654                    NA   <NA>
## HV_6          3.593      0.390       142.587                    NA   <NA>
## D31           6.443      7.805       178.264              56873318     t4
## D40           4.188      6.227       137.577              61210784     t5
## D4            6.542      8.752       173.891              50262020   <NA>
## D11           6.182      4.307       148.215              55159700     t1
## D18           5.970      0.436       100.356              42389522     t1
## D25           5.250      0.000        99.062              49858296     t1
## HV_7          4.272      1.059       113.091                    NA     t1
## HV_8          6.617      1.054        92.766                    NA     t1
## HV_9          5.648      1.069       113.632                    NA     t1
## HV_10         2.750      1.147        83.032                    NA     t1
## HV_11         4.203      1.259       159.798                    NA     t1
## HV_12         3.415      1.293       154.250                    NA     t1
## D32           7.161      0.000        85.553              50291070     t4
## D39           4.591      0.000        71.475              64514742     t5
## C_1           6.637      7.112       175.692                    NA   pret
## D3            6.110      7.167       171.401              54057756   <NA>
## D12           5.858      3.135       146.776              45384782     t1
## D19           5.905      0.337        93.554              54226192     t1
## C_6           5.818      0.000        82.391                    NA     t1
## D26           5.231      0.000        86.221              50932862     t1
## D33           7.101      0.000        81.443              59788084     t4
## HV_13         5.593      1.046        67.122                    NA     t1
## HV_14         5.969      1.079        92.785                    NA     t1
## HV_15         5.601      1.023        72.631                    NA     t1
## HV_16         4.042      0.958        59.684                    NA     t1
## HV_17         4.140      1.393       158.898                    NA     t1
## D38           4.674      0.000        61.853              38395332     t5
## HV_18         2.336      1.288       159.405                    NA     t1
## D_1              NA         NA            NA                    NA   <NA>
## D2            5.748      7.114       177.040              60148698   <NA>
## D13           6.024      6.447       158.540              52085624     t1
## D20           6.040      0.000       139.713              46655108     t1
## D27           5.309      0.000       148.633              70694598     t1
## HV_19         6.903      1.239       104.387                    NA     t1
## HV_20         7.904      1.243       107.308                    NA     t1
## HV_21         3.751      0.973        61.732                    NA     t1
## HV_22         4.026      0.971        58.946                    NA     t1
## HV_23         3.707      1.426       154.371                    NA     t1
## D34           7.673      7.628       177.547              56546096     t4
## D37           5.541      6.703       141.938              51876272     t5
## HV_24         2.300      1.239       143.500                    NA     t1
## C_2              NA         NA            NA                    NA     t5
## C_3              NA         NA            NA                    NA     t5
## C_4              NA         NA            NA                    NA     t5
## C_5              NA         NA            NA                    NA     t5
##       Reactor_Treatment_Dose   Treatment_Dose
## D1              CR_UNTREATED        UNTREATED
## D14             CR_UNTREATED        UNTREATED
## D21             CR_UNTREATED        UNTREATED
## D28             CR_UNTREATED        UNTREATED
## D35             CR_UNTREATED        UNTREATED
## D36             CR_UNTREATED        UNTREATED
## D43                    DONOR            DONOR
## HC_1       TR1_CTX+HV292.120    CTX+HV292.120
## HC_2               TR2_CTX20            CTX20
## HC_3      TR3_CTX+HV292.1200   CTX+HV292.1200
## HC_4              TR4_CTX200           CTX200
## HC_5             TR5_HV292.1          HV292.1
## HC_6            CR_UNTREATED        UNTREATED
## HC_7       TR1_CTX+HV292.120    CTX+HV292.120
## HC_8               TR2_CTX20            CTX20
## HC_9      TR3_CTX+HV292.1200   CTX+HV292.1200
## HC_10             TR4_CTX200           CTX200
## HC_11            TR5_HV292.1          HV292.1
## D7         TR1_CTX+HV292.120    CTX+HV292.120
## D8         TR1_CTX+HV292.120    CTX+HV292.120
## D15        TR1_CTX+HV292.120    CTX+HV292.120
## D22        TR1_CTX+HV292.120    CTX+HV292.120
## D29        TR1_CTX+HV292.120    CTX+HV292.120
## HC_12           CR_UNTREATED        UNTREATED
## D42        TR1_CTX+HV292.120    CTX+HV292.120
## HC_13      TR1_CTX+HV292.120    CTX+HV292.120
## HC_14              TR2_CTX20            CTX20
## D6                 TR2_CTX20            CTX20
## D9                 TR2_CTX20            CTX20
## D16                TR2_CTX20            CTX20
## D23                TR2_CTX20            CTX20
## HC_15     TR3_CTX+HV292.1200   CTX+HV292.1200
## HC_16             TR4_CTX200           CTX200
## HC_17            TR5_HV292.1          HV292.1
## HC_18           CR_UNTREATED        UNTREATED
## HC_19      TR1_CTX+HV292.120    CTX+HV292.120
## HC_20              TR2_CTX20            CTX20
## D30                TR2_CTX20            CTX20
## HC_21     TR3_CTX+HV292.1200   CTX+HV292.1200
## HC_22             TR4_CTX200           CTX200
## HC_23            TR5_HV292.1          HV292.1
## HC_24           CR_UNTREATED        UNTREATED
## D41                TR2_CTX20            CTX20
## D5               TR3_HV292.1          HV292.1
## D10              TR3_HV292.1          HV292.1
## D17              TR3_HV292.1          HV292.1
## D24              TR3_HV292.1          HV292.1
## HV_1     TR1_VAN+CCUG5916890  VAN+CCUG5916890
## HV_2               TR2_VAN90            VAN90
## HV_3    TR3_VAN+CCUG59168600 VAN+CCUG59168600
## HV_4              TR4_VAN600           VAN600
## HV_5           TR5_CCUG59168        CCUG59168
## HV_6            CR_UNTREATED        UNTREATED
## D31              TR3_HV292.1          HV292.1
## D40              TR3_HV292.1          HV292.1
## D4                 TR4_VAN90            VAN90
## D11                TR4_VAN90            VAN90
## D18                TR4_VAN90            VAN90
## D25                TR4_VAN90            VAN90
## HV_7     TR1_VAN+CCUG5916890  VAN+CCUG5916890
## HV_8               TR2_VAN90            VAN90
## HV_9    TR3_VAN+CCUG59168600 VAN+CCUG59168600
## HV_10             TR4_VAN600           VAN600
## HV_11          TR5_CCUG59168        CCUG59168
## HV_12           CR_UNTREATED        UNTREATED
## D32                TR4_VAN90            VAN90
## D39                TR4_VAN90            VAN90
## C_1      TR5_VAN+CCUG5916890  VAN+CCUG5916890
## D3       TR5_VAN+CCUG5916890  VAN+CCUG5916890
## D12      TR5_VAN+CCUG5916890  VAN+CCUG5916890
## D19      TR5_VAN+CCUG5916890  VAN+CCUG5916890
## C_6      TR5_VAN+CCUG5916890  VAN+CCUG5916890
## D26      TR5_VAN+CCUG5916890  VAN+CCUG5916890
## D33      TR5_VAN+CCUG5916890  VAN+CCUG5916890
## HV_13    TR1_VAN+CCUG5916890  VAN+CCUG5916890
## HV_14              TR2_VAN90            VAN90
## HV_15   TR3_VAN+CCUG59168600 VAN+CCUG59168600
## HV_16             TR4_VAN600           VAN600
## HV_17          TR5_CCUG59168        CCUG59168
## D38      TR5_VAN+CCUG5916890  VAN+CCUG5916890
## HV_18           CR_UNTREATED        UNTREATED
## D_1                    DONOR            DONOR
## D2             TR6_CCUG59168        CCUG59168
## D13            TR6_CCUG59168        CCUG59168
## D20            TR6_CCUG59168        CCUG59168
## D27            TR6_CCUG59168        CCUG59168
## HV_19    TR1_VAN+CCUG5916890  VAN+CCUG5916890
## HV_20              TR2_VAN90            VAN90
## HV_21   TR3_VAN+CCUG59168600 VAN+CCUG59168600
## HV_22             TR4_VAN600           VAN600
## HV_23          TR5_CCUG59168        CCUG59168
## D34            TR6_CCUG59168        CCUG59168
## D37            TR6_CCUG59168        CCUG59168
## HV_24           CR_UNTREATED        UNTREATED
## C_2      TR5_VAN+CCUG5916890  VAN+CCUG5916890
## C_3                TR4_VAN90            VAN90
## C_4              TR3_HV292.1          HV292.1
## C_5                TR2_CTX20            CTX20
```

```r
meta %>% 
  select(-sample:-index2) -> meta
```


```r
rownames(meta) <- str_replace(rownames(meta),
                              "_", ".")
```



```r
rownames(meta) <- str_replace(rownames(meta),
                              "^C.", "C")
```


Adjust plasmid - MAGs annnotations:



```r
combined_full_fct %>% 
  column_to_rownames("ORF_ID") %>% 
  select(Contig_ID:ANVIO_CONCOCT_HQ) %>% 
  mutate(PathoFact_AMR_Resistance_mechanism_multi = ifelse(grepl(";", PathoFact_AMR_Resistance_mechanism), "multi-mech", PathoFact_AMR_Resistance_mechanism)) %>%
  mutate(PathoFact_AMR_AMR_sub_class_multi = ifelse(grepl(";", PathoFact_AMR_AMR_sub_class), "multidrug", PathoFact_AMR_AMR_sub_class)) %>%
  mutate(ANVIO_CONCOCT_HQ_clean = ifelse(PathoFact_AMR_MGE_prediction == "plasmid", NA, ANVIO_CONCOCT_HQ)) %>%
  mutate(GTDB_tax_CONCOCT_clean = ifelse(PathoFact_AMR_MGE_prediction == "plasmid", NA, GTDB_tax_CONCOCT)) %>%
  mutate(PathoFact_AMR_MGE_prediction_clean = ifelse(!is.na(ANVIO_CONCOCT_HQ), "MAG", PathoFact_AMR_MGE_prediction)) %>%
  mutate(PathoFact_AMR_MGE_prediction_clean_2 = case_when(grepl("ambiguous", PathoFact_AMR_MGE_prediction_clean) ~ "ambiguous", TRUE ~ PathoFact_AMR_MGE_prediction_clean)) %>% 
  mutate(PathoFact_AMR_MGE_prediction_clean_2 = fct_relevel(PathoFact_AMR_MGE_prediction_clean_2, c("MAG", "chromosome", "plasmid", "phage", "ambiguous", "unclassified"))) -> combined_full_fct_plas_clean
```

```
## Warning: Unknown levels in `f`: chromosome, plasmid, phage, ambiguous,
## unclassified
```

  
tax_table:


```r
combined_full_fct_plas_clean %>% 
  as.matrix() %>% 
  phyloseq::tax_table() -> tax
```

ORFs:

TPM:


```r
combined_full %>% 
  column_to_rownames("ORF_ID") %>% 
  select(starts_with("ORF_TPM.")) %>% 
  rename_at(
    # select all variables with "factor.sample" in the name
    vars(contains("ORF_TPM."))
    # use stringr::str_replace to remove factor.sample.
    #   you could do the same with base::gsub()
    , funs(str_replace(., "ORF_TPM.", ""))) %>% 
  # stringr::str_replace_all("\\s","_")
  # rename_at(
  #   # select all variables with "factor.sample" in the name
  #   vars(contains("."))
  #   , funs(str_replace_all(., "\\.", "_")))  %>% 
  as.matrix() %>% 
  phyloseq::otu_table(taxa_are_rows = TRUE) -> tpm
```

```
## Warning: `funs()` was deprecated in dplyr 0.8.0.
## Please use a list of either functions or lambdas: 
## 
##   # Simple named list: 
##   list(mean = mean, median = median)
## 
##   # Auto named with `tibble::lst()`: 
##   tibble::lst(mean, median)
## 
##   # Using lambdas
##   list(~ mean(., trim = .2), ~ median(., na.rm = TRUE))
## This warning is displayed once every 8 hours.
## Call `lifecycle::last_lifecycle_warnings()` to see where this warning was generated.
```


Num_Gi_pc


```r
combined_full %>% 
  column_to_rownames("ORF_ID") %>% 
  select(starts_with("ORF_Num_Gi_pc.")) %>% 
  rename_at(
    # select all variables with "factor.sample" in the name
    vars(contains("ORF_Num_Gi_pc."))
    # use stringr::str_replace to remove factor.sample.
    #   you could do the same with base::gsub()
    , funs(str_replace(., "ORF_Num_Gi_pc.", ""))) %>% 
  # stringr::str_replace_all("\\s","_")
  # rename_at(
  #   # select all variables with "factor.sample" in the name
  #   vars(contains("."))
  #   , funs(str_replace_all(., "\\.", "_")))  %>% 
  as.matrix() %>% 
  phyloseq::otu_table(taxa_are_rows = TRUE) -> Num_Gi_pc
```

Num_Gi


```r
combined_full %>% 
  column_to_rownames("ORF_ID") %>% 
  select(starts_with("ORF_Num_Gi.")) %>% 
  rename_at(
    # select all variables with "factor.sample" in the name
    vars(contains("ORF_Num_Gi."))
    # use stringr::str_replace to remove factor.sample.
    #   you could do the same with base::gsub()
    , funs(str_replace(., "ORF_Num_Gi.", ""))) %>% 
  # stringr::str_replace_all("\\s","_")
  # rename_at(
  #   # select all variables with "factor.sample" in the name
  #   vars(contains("."))
  #   , funs(str_replace_all(., "\\.", "_")))  %>% 
  as.matrix() %>% 
  phyloseq::otu_table(taxa_are_rows = TRUE) -> Num_Gi
```

Coverage


```r
combined_full %>% 
  column_to_rownames("ORF_ID") %>% 
  select(starts_with("ORF_Coverage.")) %>% 
  rename_at(
    # select all variables with "factor.sample" in the name
    vars(contains("ORF_Coverage."))
    # use stringr::str_replace to remove factor.sample.
    #   you could do the same with base::gsub()
    , funs(str_replace(., "ORF_Coverage.", ""))) %>% 
  # stringr::str_replace_all("\\s","_")
  # rename_at(
  #   # select all variables with "factor.sample" in the name
  #   vars(contains("."))
  #   , funs(str_replace_all(., "\\.", "_")))  %>% 
  as.matrix() %>% 
  phyloseq::otu_table(taxa_are_rows = TRUE) -> cov
```

ORF_Raw.read


```r
combined_full %>% 
  column_to_rownames("ORF_ID") %>% 
  select(starts_with("ORF_Raw.read.count.")) %>% 
  rename_at(
    # select all variables with "factor.sample" in the name
    vars(contains("ORF_Raw.read.count."))
    # use stringr::str_replace to remove factor.sample.
    #   you could do the same with base::gsub()
    , funs(str_replace(., "ORF_Raw.read.count.", ""))) %>% 
  # stringr::str_replace_all("\\s","_")
  # rename_at(
  #   # select all variables with "factor.sample" in the name
  #   vars(contains("."))
  #   , funs(str_replace_all(., "\\.", "_")))  %>% 
  as.matrix() %>% 
  phyloseq::otu_table(taxa_are_rows = TRUE) -> read_count
```



Create phyloseq objects and tsv files:


```r
dim(tpm); dim(meta)
```

```
## [1] 2426165      98
```

```
## [1] 98 39
```

```r
physeq_tpm <- phyloseq::phyloseq(tpm,
                                 tax,
                                 meta %>% phyloseq::sample_data())


physeq_tpm
```

```
## phyloseq-class experiment-level object
## otu_table()   OTU Table:         [ 2426165 taxa and 98 samples ]
## sample_data() Sample Data:       [ 98 samples by 39 sample variables ]
## tax_table()   Taxonomy Table:    [ 2426165 taxa by 99 taxonomic ranks ]
```


```r
dim(cov); dim(meta)
```

```
## [1] 2426165      98
```

```
## [1] 98 39
```

```r
physeq_cov <- phyloseq::phyloseq(cov,
                                 tax,
                                 meta %>% phyloseq::sample_data())

physeq_cov
```

```
## phyloseq-class experiment-level object
## otu_table()   OTU Table:         [ 2426165 taxa and 98 samples ]
## sample_data() Sample Data:       [ 98 samples by 39 sample variables ]
## tax_table()   Taxonomy Table:    [ 2426165 taxa by 99 taxonomic ranks ]
```



```r
dim(Num_Gi_pc); dim(meta)
```

```
## [1] 2426165      98
```

```
## [1] 98 39
```

```r
physeq_Num_Gi_pc <- phyloseq::phyloseq(Num_Gi_pc,
                                       tax,
                                       meta %>% phyloseq::sample_data())

physeq_Num_Gi_pc
```

```
## phyloseq-class experiment-level object
## otu_table()   OTU Table:         [ 2426165 taxa and 98 samples ]
## sample_data() Sample Data:       [ 98 samples by 39 sample variables ]
## tax_table()   Taxonomy Table:    [ 2426165 taxa by 99 taxonomic ranks ]
```


```r
dim(Num_Gi); dim(meta)
```

```
## [1] 2426165      98
```

```
## [1] 98 39
```

```r
physeq_Num_Gi<- phyloseq::phyloseq(Num_Gi,
                                   tax,
                                   meta %>% phyloseq::sample_data())

physeq_Num_Gi
```

```
## phyloseq-class experiment-level object
## otu_table()   OTU Table:         [ 2426165 taxa and 98 samples ]
## sample_data() Sample Data:       [ 98 samples by 39 sample variables ]
## tax_table()   Taxonomy Table:    [ 2426165 taxa by 99 taxonomic ranks ]
```



```r
dim(read_count); dim(meta)
```

```
## [1] 2426165      98
```

```
## [1] 98 39
```

```r
physeq_read_count <- phyloseq::phyloseq(read_count,
                                        tax,
                                        meta %>% phyloseq::sample_data())

physeq_read_count
```

```
## phyloseq-class experiment-level object
## otu_table()   OTU Table:         [ 2426165 taxa and 98 samples ]
## sample_data() Sample Data:       [ 98 samples by 39 sample variables ]
## tax_table()   Taxonomy Table:    [ 2426165 taxa by 99 taxonomic ranks ]
```



```r
phyloseq_list_out = list("read_count" = physeq_read_count,
                         "Num_Gi" = physeq_Num_Gi,
                         "Num_Gi_pc" = physeq_Num_Gi_pc,
                         "cov" = physeq_cov,
                         "tpm" = physeq_tpm)
```


Export outputs:


```r
phyloseq_list_out %>% 
  saveRDS(here::here("data/processed/exp1_full_gene_catalog_phyloseq.RDS"))

combined_full_fct_plas_clean %>% 
  saveRDS(here::here("data/processed/exp1_full_gene_catalog_full_metrics_table.RDS"))

# combined_full_fct_plas_clean %>% 
#   write_tsv(here::here("data/processed/exp1_full_gene_catalog_full_metrics_table.tsv.gz"))
```



```r
# combined_full %>% 
#   filter(grepl("integrase",anvio_PFAM)) %>% 
#   write_tsv("~/Desktop/anvio_pfam_integrase.tsv")
```

Export only AMR and Tox:


```r
# combined_full_fct_plas_clean %>% 
#   filter(!is.na(PathoFact_AMR_ARG) | !is.na(rgi_CARD_ARO) | !is.na(PathoFact_Tox_Toxin_classification)) %>% 
#   write_tsv(here::here("data/processed/exp1_full_gene_catalog_full_metrics_table_AMR_tox.tsv.gz"))
```



```r
devtools::session_info()
```

```
## ─ Session info ───────────────────────────────────────────────────────────────
##  setting  value
##  version  R version 4.2.0 (2022-04-22)
##  os       macOS Mojave 10.14.6
##  system   x86_64, darwin17.0
##  ui       X11
##  language (EN)
##  collate  en_US.UTF-8
##  ctype    en_US.UTF-8
##  tz       Europe/Zurich
##  date     2022-07-06
##  pandoc   2.14.0.3 @ /Applications/RStudio.app/Contents/MacOS/pandoc/ (via rmarkdown)
## 
## ─ Packages ───────────────────────────────────────────────────────────────────
##  package          * version    date (UTC) lib source
##  ade4               1.7-19     2022-04-19 [1] CRAN (R 4.2.0)
##  ape                5.6-2      2022-03-02 [1] CRAN (R 4.2.0)
##  assertthat         0.2.1      2019-03-21 [1] CRAN (R 4.2.0)
##  backports          1.4.1      2021-12-13 [1] CRAN (R 4.2.0)
##  Biobase            2.56.0     2022-04-26 [1] Bioconductor
##  BiocGenerics       0.42.0     2022-04-26 [1] Bioconductor
##  biomformat         1.24.0     2022-04-26 [1] Bioconductor
##  Biostrings         2.64.0     2022-04-26 [1] Bioconductor
##  bit                4.0.4      2020-08-04 [1] CRAN (R 4.2.0)
##  bit64              4.0.5      2020-08-30 [1] CRAN (R 4.2.0)
##  bitops             1.0-7      2021-04-24 [1] CRAN (R 4.2.0)
##  brio               1.1.3      2021-11-30 [1] CRAN (R 4.2.0)
##  broom              1.0.0      2022-07-01 [1] CRAN (R 4.2.0)
##  bslib              0.3.1      2021-10-06 [1] CRAN (R 4.2.0)
##  cachem             1.0.6      2021-08-19 [1] CRAN (R 4.2.0)
##  callr              3.7.0      2021-04-20 [1] CRAN (R 4.2.0)
##  cellranger         1.1.0      2016-07-27 [1] CRAN (R 4.2.0)
##  cli                3.3.0      2022-04-25 [1] CRAN (R 4.2.0)
##  cluster            2.1.3      2022-03-28 [1] CRAN (R 4.2.0)
##  codetools          0.2-18     2020-11-04 [1] CRAN (R 4.2.0)
##  colorspace         2.0-3      2022-02-21 [1] CRAN (R 4.2.0)
##  crayon             1.5.1      2022-03-26 [1] CRAN (R 4.2.0)
##  data.table         1.14.2     2021-09-27 [1] CRAN (R 4.2.0)
##  DBI                1.1.2      2021-12-20 [1] CRAN (R 4.2.0)
##  dbplyr             2.2.0      2022-06-05 [1] CRAN (R 4.2.0)
##  desc               1.4.1      2022-03-06 [1] CRAN (R 4.2.0)
##  devtools           2.4.3      2021-11-30 [1] CRAN (R 4.2.0)
##  digest             0.6.29     2021-12-01 [1] CRAN (R 4.2.0)
##  dplyr            * 1.0.9      2022-04-28 [1] CRAN (R 4.2.0)
##  dtplyr             1.2.1      2022-01-19 [1] CRAN (R 4.2.0)
##  ellipsis           0.3.2      2021-04-29 [1] CRAN (R 4.2.0)
##  evaluate           0.15       2022-02-18 [1] CRAN (R 4.2.0)
##  fansi              1.0.3      2022-03-24 [1] CRAN (R 4.2.0)
##  fastmap            1.1.0      2021-01-25 [1] CRAN (R 4.2.0)
##  forcats          * 0.5.1      2021-01-27 [1] CRAN (R 4.2.0)
##  foreach            1.5.2      2022-02-02 [1] CRAN (R 4.2.0)
##  fs                 1.5.2      2021-12-08 [1] CRAN (R 4.2.0)
##  gargle             1.2.0      2021-07-02 [1] CRAN (R 4.2.0)
##  generics           0.1.2      2022-01-31 [1] CRAN (R 4.2.0)
##  GenomeInfoDb       1.32.2     2022-05-15 [1] Bioconductor
##  GenomeInfoDbData   1.2.8      2022-06-13 [1] Bioconductor
##  ggplot2          * 3.3.6      2022-05-03 [1] CRAN (R 4.2.0)
##  glue               1.6.2      2022-02-24 [1] CRAN (R 4.2.0)
##  googledrive        2.0.0      2021-07-08 [1] CRAN (R 4.2.0)
##  googlesheets4      1.0.0      2021-07-21 [1] CRAN (R 4.2.0)
##  gtable             0.3.0      2019-03-25 [1] CRAN (R 4.2.0)
##  haven              2.5.0      2022-04-15 [1] CRAN (R 4.2.0)
##  here             * 1.0.1      2020-12-13 [1] CRAN (R 4.2.0)
##  hms                1.1.1      2021-09-26 [1] CRAN (R 4.2.0)
##  htmltools          0.5.2      2021-08-25 [1] CRAN (R 4.2.0)
##  httr               1.4.3      2022-05-04 [1] CRAN (R 4.2.0)
##  igraph             1.3.1      2022-04-20 [1] CRAN (R 4.2.0)
##  IRanges            2.30.0     2022-04-26 [1] Bioconductor
##  iterators          1.0.14     2022-02-05 [1] CRAN (R 4.2.0)
##  jquerylib          0.1.4      2021-04-26 [1] CRAN (R 4.2.0)
##  jsonlite           1.8.0      2022-02-22 [1] CRAN (R 4.2.0)
##  knitr              1.39       2022-04-26 [1] CRAN (R 4.2.0)
##  lattice            0.20-45    2021-09-22 [1] CRAN (R 4.2.0)
##  lifecycle          1.0.1      2021-09-24 [1] CRAN (R 4.2.0)
##  lubridate          1.8.0      2021-10-07 [1] CRAN (R 4.2.0)
##  magrittr           2.0.3      2022-03-30 [1] CRAN (R 4.2.0)
##  MASS               7.3-57     2022-04-22 [1] CRAN (R 4.2.0)
##  Matrix             1.4-1      2022-03-23 [1] CRAN (R 4.2.0)
##  memoise            2.0.1      2021-11-26 [1] CRAN (R 4.2.0)
##  mgcv               1.8-40     2022-03-29 [1] CRAN (R 4.2.0)
##  modelr             0.1.8      2020-05-19 [1] CRAN (R 4.2.0)
##  multtest           2.52.0     2022-04-26 [1] Bioconductor
##  munsell            0.5.0      2018-06-12 [1] CRAN (R 4.2.0)
##  nlme               3.1-157    2022-03-25 [1] CRAN (R 4.2.0)
##  permute            0.9-7      2022-01-27 [1] CRAN (R 4.2.0)
##  phyloseq           1.40.0     2022-04-26 [1] Bioconductor
##  pillar             1.7.0      2022-02-01 [1] CRAN (R 4.2.0)
##  pkgbuild           1.3.1      2021-12-20 [1] CRAN (R 4.2.0)
##  pkgconfig          2.0.3      2019-09-22 [1] CRAN (R 4.2.0)
##  pkgload            1.2.4      2021-11-30 [1] CRAN (R 4.2.0)
##  plyr               1.8.7      2022-03-24 [1] CRAN (R 4.2.0)
##  prettyunits        1.1.1      2020-01-24 [1] CRAN (R 4.2.0)
##  processx           3.6.0      2022-06-10 [1] CRAN (R 4.2.0)
##  ps                 1.7.0      2022-04-23 [1] CRAN (R 4.2.0)
##  purrr            * 0.3.4      2020-04-17 [1] CRAN (R 4.2.0)
##  R6                 2.5.1      2021-08-19 [1] CRAN (R 4.2.0)
##  Rcpp               1.0.8.3    2022-03-17 [1] CRAN (R 4.2.0)
##  RCurl              1.98-1.7   2022-06-09 [1] CRAN (R 4.2.0)
##  readr            * 2.1.2      2022-01-30 [1] CRAN (R 4.2.0)
##  readxl             1.4.0      2022-03-28 [1] CRAN (R 4.2.0)
##  remotes            2.4.2      2021-11-30 [1] CRAN (R 4.2.0)
##  reprex             2.0.1      2021-08-05 [1] CRAN (R 4.2.0)
##  reshape2           1.4.4      2020-04-09 [1] CRAN (R 4.2.0)
##  rhdf5              2.40.0     2022-04-26 [1] Bioconductor
##  rhdf5filters       1.8.0      2022-04-26 [1] Bioconductor
##  Rhdf5lib           1.18.2     2022-05-17 [1] Bioconductor
##  rlang              1.0.3      2022-06-27 [1] CRAN (R 4.2.0)
##  rmarkdown          2.14       2022-04-25 [1] CRAN (R 4.2.0)
##  rprojroot          2.0.3      2022-04-02 [1] CRAN (R 4.2.0)
##  rstudioapi         0.13       2020-11-12 [1] CRAN (R 4.2.0)
##  rvest              1.0.2      2021-10-16 [1] CRAN (R 4.2.0)
##  S4Vectors          0.34.0     2022-04-26 [1] Bioconductor
##  sass               0.4.1      2022-03-23 [1] CRAN (R 4.2.0)
##  scales             1.2.0      2022-04-13 [1] CRAN (R 4.2.0)
##  sessioninfo        1.2.2      2021-12-06 [1] CRAN (R 4.2.0)
##  stringi            1.7.6      2021-11-29 [1] CRAN (R 4.2.0)
##  stringr          * 1.4.0      2019-02-10 [1] CRAN (R 4.2.0)
##  survival           3.3-1      2022-03-03 [1] CRAN (R 4.2.0)
##  testthat           3.1.4      2022-04-26 [1] CRAN (R 4.2.0)
##  tibble           * 3.1.7      2022-05-03 [1] CRAN (R 4.2.0)
##  tidyr            * 1.2.0      2022-02-01 [1] CRAN (R 4.2.0)
##  tidyselect         1.1.2      2022-02-21 [1] CRAN (R 4.2.0)
##  tidyverse        * 1.3.1.9000 2022-06-13 [1] Github (tidyverse/tidyverse@6186fbf)
##  tzdb               0.3.0      2022-03-28 [1] CRAN (R 4.2.0)
##  usethis            2.1.6      2022-05-25 [1] CRAN (R 4.2.0)
##  utf8               1.2.2      2021-07-24 [1] CRAN (R 4.2.0)
##  vctrs              0.4.1      2022-04-13 [1] CRAN (R 4.2.0)
##  vegan              2.6-2      2022-04-17 [1] CRAN (R 4.2.0)
##  vroom              1.5.7      2021-11-30 [1] CRAN (R 4.2.0)
##  withr              2.5.0      2022-03-03 [1] CRAN (R 4.2.0)
##  xfun               0.31       2022-05-10 [1] CRAN (R 4.2.0)
##  xml2               1.3.3      2021-11-30 [1] CRAN (R 4.2.0)
##  XVector            0.36.0     2022-04-26 [1] Bioconductor
##  yaml               2.3.5      2022-02-21 [1] CRAN (R 4.2.0)
##  zlibbioc           1.42.0     2022-04-26 [1] Bioconductor
## 
##  [1] /Library/Frameworks/R.framework/Versions/4.2/Resources/library
## 
## ──────────────────────────────────────────────────────────────────────────────
```

