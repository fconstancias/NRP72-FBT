---
title: "NTP72 - 16S alpha div"
author: "Florentin Constancias"
date: "August 04, 2022"
output: 
  html_document: 
    toc: yes
    keep_md: yes
---



# Load required packages


```r
library(tidyverse)
# library(phyloseq)
library(speedyseq)
library(ggrepel)
library(here)
# library(microViz)
# library(RColorBrewer)
# library(vegan)
library(randomcoloR)
options(getClass.msg=FALSE) # https://github.com/epurdom/clusterExperiment/issues/66
#this fixes an error message that pops up because the class 'Annotated' is defined in two different packages
```

# Source function


```r
rm(list = ls())

source("https://raw.githubusercontent.com/fconstancias/DivComAnalyses/master/R/phyloseq_taxa_tests.R")
source("https://raw.githubusercontent.com/fconstancias/DivComAnalyses/master/R/phyloseq_normalisation.R")
source("https://raw.githubusercontent.com/fconstancias/DivComAnalyses/master/R/phyloseq_alpha.R")
source("https://raw.githubusercontent.com/fconstancias/DivComAnalyses/master/R/phyloseq_beta.R")
source("https://raw.githubusercontent.com/fconstancias/DivComAnalyses/master/R/phyloseq_heatmap.R")
source("https://raw.githubusercontent.com/fconstancias/DivComAnalyses/master/R/phyloseq_varia.R")
```



```r
export_ppt = TRUE
export_ppt_width = NULL
export_ppt_height = NULL
out_pptx = here::here("output/plots.pptx")
```

# Import data + cleaning/ colors/ factors:


```r
"Figures/Rastetics.Rdata" %>%
  here::here() %>% 
  load()
```

## Data

### 16S:

```r
"data/processed/16S/ps_silva_dada2_human_chicken_all_fct_polyferms.RDS" %>% 
  here::here() %>% 
  readRDS() %>% 
  filter_taxa(function(x) sum(x > 0) > 0, TRUE) -> ps_16S
```


```r
sample_variables(ps_16S)
```

```
##  [1] "sample_name"                 "Sample_description"         
##  [3] "I7_Index_ID"                 "index"                      
##  [5] "I5_Index_ID"                 "index2"                     
##  [7] "Description2"                "Experiment"                 
##  [9] "Reactor"                     "Treatment"                  
## [11] "Day_of_Connection"           "Day_of_Treatment"           
## [13] "Day_from_Inoculum"           "Enrichment"                 
## [15] "Phase"                       "Treatment2"                 
## [17] "Date"                        "Paul"                       
## [19] "Reactor_Treatment"           "GeneCopyNumberperML"        
## [21] "HV292.1_Copy_Number_permL"   "CCUG59168_Copy_Number_permL"
## [23] "CTX_Copy_Number_permL"       "VAN_Copy_Number_permL"      
## [25] "Model"                       "Model2"                     
## [27] "Antibiotic_mg.L"             "Fermentation"               
## [29] "Antibiotic"                  "Lactose_mM"                 
## [31] "Glucose_mM"                  "Galactose_mM"               
## [33] "Succinat_mM"                 "Lactat_mM"                  
## [35] "Formiat_mM"                  "Acetat_mM"                  
## [37] "Propionat_mM"                "Isobutyrat_mM"              
## [39] "Butyrat_mM"                  "Isovalerat_mM"              
## [41] "Valerat_mM"                  "Total_SCFA_mM"              
## [43] "metagenomic_sample_name"     "raw_metagenomic_pairs"      
## [45] "Period"                      "Reactor_Treatment_Dose"     
## [47] "Treatment_Dose"
```


```r
ps_16S %>%
  # subset_samples(Experiment %!in% c("Mock", "NTC")) %>%
  phyloseq_check_lib_size(data_color = "Experiment", 
                          data_facet = NULL, 
                          nreads_display = 4000, 
                          first_n = 500) -> out

out$df %>%
  dplyr::select(SampleID, Experiment, Model2, LibrarySize) %>%
  # filter(category2 %in% c("feces")) %>%
  DT::datatable()
```

```{=html}
<div id="htmlwidget-0977ae5a3fa15a34b43c" style="width:100%;height:auto;" class="datatables html-widget"></div>
<script type="application/json" data-for="htmlwidget-0977ae5a3fa15a34b43c">{"x":{"filter":"none","vertical":false,"data":[["363","598","708","709","710","712","713","714","715","716","717","718","721","559","221","233","239","255","334","360","379","256","275","509","359","719","340","381","711","269","328","235","246","204","361","364","245","336","220","257","195","389","411","21","393","312","720","236","259","300","62","57","397","48","130","252","96","115","116","280","229","61","171","319","414","193","341","192","385","417","723","392","707","351","552","724","343","231","515","237","706","329","200","25","356","722","318","369","306","368","19","194","266","394","214","490","289","705","301","105","619","703","478","502","700","662","701","704","699","696","702","638","698","681","594","661","652","626","625","618","697","645","663","455","326","671","566","590","601","617","179","654","575","573","665","570","616","588","624","11","688","629","627","679","600","675","551","643","647","634","577","651","608","622","584","498","512","660","581","644","611","666","567","655","610","602","596","607","56","593","119","38","33","636","576","606","523","637","65","603","331","410","628","609","185","450","641","58","684","550","578","579","630","302","104","668","582","92","64","631","691","88","563","613","680","689","522","677","542","10","24","49","23","648","585","592","673","612","690","620","599","532","635","160","166","604","686","71","657","669","589","516","615","518","656","42","653","554","621","546","614","571","370","633","493","182","672","683","539","279","26","74","574","564","144","533","670","55","572","14","53","562","362","580","35","556","583","586","543","544","685","569","519","196","136","448","31","314","321","534","139","189","526","527","1","535","540","395","78","73","159","100","595","549","108","121","158","187","17","545","164","513","198","415","39","135","646","431","186","658","90","168","676","262","45","547","531","503","94","206","234","145","692","565","27","530","97","134","87","553","398","43","650","15","63","536","555","682","372","16","128","146","338","304","642","268","6","407","568","591","659","649","349","40","409","106","587","9","678","694","499","386","378","8","178","91","470","674","639","640","125","309","528","605","249","371","548","467","222","632","162","4","416","687","597","466","664","401","12","67","325","175","118","524","408","560","163","501","525","263","506","83","132","202","210","18","183","47","126","111","298","345","537","434","98","453","284","505","514","75","173","102","517","529","667","36","127","82","103","330","133","285","353","167","161","288","307","436","169","695","403","335","366","240","561","66","227","209","623","137","37","112","156","339","311","32","5","172","375","538","437","485","140","438","250","51","347","157","433","101","153","402","212","20","373","191","435","271","303","203","441","333","254","69","484","327","177","260","454","142","174","404","291","293","223","388","208","120","147","99","313","274","482","50","95","213","2","44","383","510","184","149","155","190","384","324","391","477","277","30","272","216","412","114","13","474","261","459","242","317","427","124","267","365","396","447","354","464","131","176","201","215","296","473","270","400","151","458","248","122","129","476","463","460","558","310","41","405","439","511","344","446","557","76","219","465","443","265","46","486","352","429","430","355","113","481","59","451","93","138","34","109","251","442","70","84","244","28","693","456","483","480","461","89","278","29","471","86","421","281","80","521","475","457","297","241","224","376","81","218","77","52","110","264","148","22","387","165","286","425","423","141","504","452","500","440","238","7","444","143","154","488","358","320","422","346","68","79","520","316","295","419","211","332","322","367","152","3","449","489","377","54","283","315","496","479","276","85","243","323","472","494","432","445","406","426","491","170","492","197","350","495","390","60","230","507","226","469","413","199","399","428","508","258","305","357","123","181","497","228","541","424","420","292","462","180","308","374","207","487","342","117","290","468","188","232","72","107","299","418","225","217","150","282","205","337","287","348","253","382","247","380","273","294"],["TR5-68-S99","AP-129-S135","NC-12-S286","NC-13-S287","NC-14-S288","NC-16-S290","NC-17-S332","NC-18-S333","NC-19-S334","NC-2-S98","NC-20-S335","NC-3-S99","NC-6-S102","p-3-H12-S288","IR-74-S188","Negative-2-S104","TR1-19-S116","TR1-63-S100","TR4-65-S155","TR5-65-S187","TR6-36-S124","TR1-64-S120","TR2-35-S131","p-2-H8-S188","TR5-64-S108","NC-4-S100","TR5-01-S123","TR6-60-S119","NC-15-S289","TR2-29-S117","TR4-59-S180","Negative-7-S223","TR1-32-S111","IR-13-S103","TR5-66-S101","TR5-69-S107","TR1-31-S109","TR4-67-S192","IR-51-S157","TR1-65-S128","negative4-S373","TR6-68-S127","p-1-F12-S72","CR-52-S196","negative-4-S189","TR3-68-S194","NC-5-S101","Negative-8-S224","TR1-67-S133","TR3-34-S149","TR1-16-S199","Positive2-S163","positive-2-S181","IR1-69-S198","TR4-18-S200","TR1-60-S114","TR2-34-S302","TR3-28-S162","TR3-3-S179","TR2-62-S97","IR-82-S143","TR1-15-S168","TR6-13-S164","TR4-28-S138","p-1-G12-S84","negative2-S181","TR5-13-S137","negativ-S80","TR6-64-S105","p-1-H12-S96","NC-8-S140","negative-1-S58","NC-11-S192","TR5-34-S28","p-3-G12-S276","NC-9-S190","TR5-26-S24","IR-84-S178","p-3-A5-S197","TR1-01-S102","NC-10-S191","TR4-60-S36","D-1-S49","EMPTY-S384","TR5-61-S142","NC-7-S139","TR4-27-S130","TR6-26-S122","TR3-62-S37","TR6-19-S33","CR-40-S166","negative3-S281","TR2-26-S87","negative-5-S191","IR-46-S126","p-2-G12-S180","TR3-01-S121","NC-1-S97","TR3-35-S20","TR3-14-S136","AP-29-S29","MOCK-8-S330","p-2-F12-S168","p-2-H12-S192","MOCK-5-S284","AP-68-S68","MOCK-6-S285","MOCK-9-S331","MOCK-4-S283","MOCK-1-S187","MOCK-7-S329","AP-46-S46","MOCK-3-S189","AP-86-S86","AP-125-S131","AP-67-S67","AP-59-S59","AP-35-S35","AP-34-S34","AP-28-S28","MOCK-2-S188","AP-52-S52","AP-69-S69","p-2-D11-S143","TR4-35-S202","AP-76-S76","AP-10-S10","AP-121-S127","AP-131-S137","AP-27-S27","TR6-21-S107","AP-60-S60","AP-108-S114","AP-106-S112","AP-70-S70","AP-103-S109","AP-26-S26","AP-12-S12","AP-33-S33","CR-21-S126","AP-92-S92","AP-38-S38","AP-36-S36","AP-84-S84","AP-130-S136","AP-80-S80","p-3-G1-S265","AP-50-S50","AP-54-S54","AP-42-S42","AP-11-S11","AP-58-S58","AP-19-S19","AP-31-S31","AP-116-S122","p-2-G9-S177","p-3-A2-S194","AP-66-S66","AP-113-S119","AP-51-S51","AP-21-S21","AP-71-S71","AP-100-S106","AP-61-S61","AP-20-S20","AP-132-S138","AP-127-S133","AP-18-S18","IR2-4-S206","AP-124-S130","TR3-37-S221","IR1-39-S351","IR1-3-S320","AP-44-S44","AP-109-S115","AP-17-S17","p-3-B6-S210","AP-45-S45","TR1-19-S305","AP-14-S14","TR4-62-S1","p-1-F11-S71","AP-37-S37","AP-2-S2","TR6-34-S143","p-2-C7-S127","AP-49-S49","TR1-1-S110","AP-89-S89","p-3-F6-S258","AP-110-S116","AP-111-S117","AP-39-S39","TR3-36-S8","TR3-13-S304","AP-73-S73","AP-114-S120","TR2-25-S133","TR1-18-S124","AP-4-S4","AP-95-S95","TR2-19-S318","p-3-H5-S281","AP-23-S23","AP-85-S85","AP-93-S93","p-3-B5-S209","AP-82-S82","p-3-E5-S245","CR-20-S111","D-0-S154","IR1-77-S202","CR-64-S324","AP-55-S55","AP-117-S123","AP-123-S129","AP-79-S79","AP-22-S22","AP-94-S94","AP-3-S3","AP-13-S13","p-3-D2-S230","AP-43-S43","TR5-34-S127","TR5-57-S348","AP-15-S15","AP-90-S90","TR1-3-S121","AP-63-S63","AP-74-S74","AP-120-S126","p-3-A6-S198","AP-25-S25","p-3-B1-S205","AP-62-S62","IR1-45-S128","AP-6-S6","p-3-G3-S267","AP-30-S30","p-3-F2-S254","AP-24-S24","AP-104-S110","TR6-27-S6","AP-41-S41","p-2-G4-S172","TR6-28-S158","AP-77-S77","AP-88-S88","p-3-E2-S242","TR2-61-S196","HV292-1-S350","TR1-34-S317","AP-107-S113","p-3-H6-S282","TR4-58-S112","p-3-D3-S231","AP-75-S75","IR2-3-S347","AP-105-S111","CR-28-S106","IR2-18-S135","p-3-H4-S280","TR5-67-S96","AP-112-S118","IR1-36-S159","p-3-G5-S269","AP-115-S121","AP-118-S124","p-3-E6-S246","p-3-F1-S253","AP-9-S9","AP-102-S108","p-3-B2-S206","positive1-S40","TR4-28-S183","p-2-C5-S125","IR1-24-S295","TR3-7-S204","TR4-30-S48","p-3-D4-S232","TR4-34-S342","TR6-52-S122","p-3-C4-S220","p-3-C5-S221","CCUG-59168-S352","p-3-D5-S233","p-3-E3-S243","negative-6-S213","TR1-52-S253","TR1-33-S375","TR5-31-S313","TR2-52-S230","AP-126-S132","p-3-F5-S257","TR3-17-S222","TR3-46-S224","TR5-3-S331","TR6-40-S160","CR-34-S359","p-3-F12-S264","TR5-52-S93","p-3-A3-S195","positive4-S301","p-1-H10-S94","IR1-4-S223","TR4-25-S333","AP-53-S53","p-2-B11-S119","TR6-37-S312","AP-64-S64","TR2-21-S228","TR5-59-S374","AP-81-S81","TR1-7-S153","IR1-57-S356","p-3-F3-S255","p-3-D12-S240","p-2-H2-S182","TR2-3-S297","IR-22-S41","Negative-3-S150","TR4-64-S220","AP-96-S96","AP-1-S1","IR1-16-S307","p-3-D1-S229","TR2-37-S330","TR4-22-S344","TR2-18-S326","p-3-G2-S266","positive-3-S212","IR1-48-S339","AP-57-S57","CR-3-S98","TR1-17-S246","p-3-D6-S234","p-3-G4-S268","AP-87-S87","TR6-29-S207","CR-31-S334","TR4-16-S276","TR5-1-S210","TR4-69-S197","TR3-60-S46","AP-5-S5","TR2-28-S64","CR-16-S332","p-1-E11-S59","AP-101-S107","AP-122-S128","AP-65-S65","AP-56-S56","TR5-32-S56","IR1-43-S256","p-1-F10-S70","TR3-15-S338","AP-119-S125","CR-19-S341","AP-83-S83","AP-98-S104","p-2-H1-S181","TR6-65-S94","TR6-35-S208","CR-18-S252","TR6-20-S132","TR2-22-S364","p-2-E5-S149","AP-8-S8","AP-47-S47","AP-48-S48","TR4-13-S142","TR3-65-S92","p-3-C6-S222","AP-16-S16","TR1-35-S161","TR6-28-S80","p-3-F4-S256","p-2-E2-S146","IR-75-S44","AP-40-S40","TR5-40-S335","CR-14-S336","p-1-H11-S95","AP-91-S91","AP-128-S134","p-2-E12-S156","AP-7-S7","p-1-B11-S23","CR-22-S311","TR1-21-S207","TR4-34-S16","TR6-17-S138","TR3-34-S354","p-3-C2-S218","p-1-E12-S60","p-3-H2-S278","TR5-46-S310","p-2-H11-S191","p-3-C3-S219","TR2-01-S14","p-2-H5-S185","TR2-14-S232","TR4-20-S269","HV292-1-S205","IR-42-S9","CR-37-S367","TR6-3-S236","IR1-63-S155","TR4-14-S226","TR3-20-S265","TR3-32-S211","TR5-28-S5","p-3-E1-S241","p-2-B3-S111","TR2-40-S259","p-2-D1-S133","TR2-66-S93","p-2-H4-S184","p-3-A4-S196","TR1-37-S306","TR6-15-S149","TR2-64-S216","p-3-A7-S199","p-3-C7-S223","AP-72-S72","IR1-37-S343","TR4-15-S218","TR2-13-S231","TR3-1-S237","TR4-61-S11","TR4-21-S254","TR2-67-S193","TR5-36-S85","TR5-58-S372","TR5-37-S368","TR2-7-S215","TR3-63-S206","p-2-B5-S113","TR5-64-S141","AP-99-S105","p-1-C11-S35","TR4-66-S55","TR6-01-S3","TR1-26-S27","p-3-H3-S279","TR1-20-S118","IR-80-S67","IR-41-S106","AP-32-S32","TR4-3-S309","IR1-38-S337","TR3-21-S234","TR5-25-S240","TR4-7-S65","TR3-67-S199","IR1-26-S119","CR-15-S292","TR6-14-S211","TR6-32-S39","p-3-E12-S252","p-2-B6-S114","p-2-F8-S164","TR4-37-S188","p-2-B7-S115","TR1-36-S198","IR2-16-S100","TR5-30-S62","TR5-28-S214","p-2-B2-S110","TR2-58-S172","TR5-20-S208","p-1-B12-S24","IR-44-S89","CR-46-S191","TR6-30-S144","TR6-64-S146","p-2-B4-S112","TR2-31-S32","TR3-59-S61","IR-01-S201","p-2-C1-S121","TR4-64-S10","TR1-62-S7","TR1-25-S117","p-2-F7-S163","TR4-36-S12","TR6-19-S137","TR1-68-S13","p-2-D10-S142","TR4-46-S215","TR6-16-S371","p-1-C12-S36","TR3-19-S214","TR3-27-S203","IR-76-S2","TR6-67-S140","IR-34-S4","TR3-40-S257","TR5-13-S103","TR2-46-S239","TR3-69-S139","TR2-34-S71","p-2-F5-S161","IR1-8-S260","TR2-31-S182","IR-45-S216","CR-1-S116","IR1-54-S130","TR6-62-S50","p-2-H9-S189","TR6-31-S272","TR5-16-S290","TR5-22-S145","TR6-58-S361","TR6-63-S195","TR4-33-S84","TR6-7-S91","p-2-F11-S167","TR2-59-S60","IR1-2-S96","TR2-32-S90","IR-48-S98","p-1-G10-S82","TR3-25-S102","CR-25-S383","p-2-E9-S153","TR1-69-S145","p-2-D4-S136","TR1-28-S35","TR4-26-S26","p-2-A8-S104","TR3-64-S190","TR2-27-S59","TR5-7-S66","positive-1-S25","p-2-C4-S124","TR5-59-S31","p-2-D9-S141","TR4-19-S270","TR6-18-S262","D-2-S164","IR-47-S63","TR3-30-S190","p-2-E8-S152","TR2-30-S73","p-1-A12-S12","TR5-18-S201","p-2-D3-S135","TR1-34-S34","TR3-52-S217","TR4-17-S114","p-2-F10-S166","p-2-D8-S140","p-2-D5-S137","p-3-H1-S277","TR3-66-S69","IR1-44-S264","p-1-D11-S47","p-2-B8-S116","p-3-A1-S193","TR5-27-S148","p-2-C3-S123","p-3-G6-S270","TR1-40-S209","IR-50-S81","p-2-E1-S145","p-2-C11-S131","TR2-19-S40","IR1-6-S250","p-2-F9-S165","TR5-35-S77","p-2-B1-S109","p-2-B10-S118","TR5-60-S75","TR3-22-S235","p-2-F4-S160","TR1-13-S171","p-2-C8-S128","TR2-28-S286","TR4-31-S278","IR1-35-S370","TR3-18-S293","TR1-59-S42","p-2-C10-S130","TR1-28-S289","TR2-15-S288","TR1-30-S76","IR1-17-S279","AP-97-S103","p-2-D12-S144","p-2-F6-S162","p-2-F3-S159","p-2-D6-S138","TR2-20-S271","TR2-60-S79","IR1-18-S287","p-2-E6-S150","TR2-17-S248","p-2-A12-S108","TR2-63-S160","TR1-64-S263","p-3-B4-S208","p-2-F1-S157","p-2-D2-S134","TR3-31-S23","TR1-27-S54","IR-77-S163","TR6-33-S135","TR2-1-S273","IR-5-S70","TR1-46-S170","IR2-17-S86","TR3-19-S266","TR2-13-S156","TR5-14-S267","CR-58-S382","TR6-66-S132","TR5-56-S241","TR2-68-S147","p-2-A6-S102","p-2-A3-S99","TR4-40-S242","p-2-H3-S183","p-2-C9-S129","p-2-H10-S190","p-2-B9-S117","TR1-13-S22","CR-17-S97","p-2-C12-S132","TR4-52-S95","TR5-21-S185","p-2-G10-S178","TR5-63-S82","TR4-29-S159","p-2-A2-S98","TR5-29-S136","TR1-22-S282","TR1-58-S275","p-3-B3-S207","TR4-19-S154","TR3-29-S210","p-2-A10-S106","IR-43-S162","TR4-63-S74","TR4-31-S19","TR6-13-S21","TR5-19-S87","CR-13-S148","p-2-C6-S126","p-2-G11-S179","TR6-34-S141","IR2-2-S291","TR2-65-S158","TR4-13-S168","p-2-G7-S175","p-2-F2-S158","TR2-36-S183","TR2-16-S184","TR1-29-S18","TR4-32-S17","p-2-E7-S151","p-2-G5-S173","p-2-B12-S120","p-2-C2-S122","p-1-D12-S48","p-2-A7-S103","p-2-G2-S170","TR6-1-S247","p-2-G3-S171","positive3-S243","TR5-33-S134","p-2-G6-S174","TR6-69-S174","TR1-14-S274","IR-83-S170","p-2-H6-S186","IR-79-S167","p-2-E4-S148","p-1-G11-S83","CCUG59168-S78","p-1-A11-S11","p-2-A9-S105","p-2-H7-S187","TR1-66-S151","TR3-61-S45","TR5-62-S176","TR3-58-S233","TR6-25-S378","p-2-G8-S176","IR-81-S182","p-3-E4-S244","p-2-A4-S100","p-2-A11-S107","TR3-26-S175","p-2-D7-S139","TR6-22-S152","TR3-64-S146","TR6-31-S179","IR-28-S184","p-2-G1-S169","TR5-19-S173","TR3-31-S283","TR3-13-S129","p-2-E3-S147","TR6-46-S176","IR-9-S165","TR1-31-S153","TR3-16-S296","TR3-33-S86","p-2-A1-S97","IR-78-S29","IR-49-S125","TR5-17-S82","TR2-64-S51","IR-16-S113","TR4-68-S52","TR2-69-S57","TR5-31-S88","TR1-61-S47","TR6-61-S186","TR1-33-S177","TR6-59-S185","TR2-33-S30","TR3-28-S38"],["Continuous","Continuous","NTC","NTC","NTC","NTC","NTC","NTC","NTC","NTC","NTC","NTC","NTC","NTC","Continuous","NTC","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","NTC","Continuous","Continuous","NTC","Continuous","Continuous","NTC","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","NTC","Continuous","NTC","Continuous","NTC","Continuous","NTC","NTC","Continuous","Continuous","Continuous","Mock","Mock","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","NTC","NTC","Continuous","NTC","Continuous","NTC","NTC","NTC","NTC","Continuous","NTC","NTC","Continuous","Continuous","Continuous","Continuous","NTC","Continuous","Cecum","NTC","Continuous","NTC","Continuous","Continuous","Continuous","Continuous","Continuous","NTC","Continuous","NTC","Continuous","NTC","Continuous","NTC","Continuous","Continuous","Continuous","Mock","NTC","NTC","Mock","Continuous","Mock","Mock","Mock","Mock","Mock","Continuous","Mock","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Mock","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Cecum","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","HV292.1","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Mock","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","CCUG59168","Continuous","Continuous","NTC","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Mock","Continuous","Continuous","Mock","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Mock","Continuous","Continuous","Continuous","NTC","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Mock","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Mock","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Mock","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","HV292.1","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Cecum","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Mock","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Mock","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Mock","Continuous","Continuous","Continuous","Continuous","Continuous","Cecum","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Mock","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Mock","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Mock","Continuous","Continuous","Continuous","Continuous","Mock","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","CCUG59168","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous"],["Human1","Human2","Human2","Human2","Human2","Human2","Human2","Human2","Human2","Human2","Human2","Human2","Human2","Chicken2","Human1","Human1","Human1","Human1","Human1","Human1","Human1","Human1","Human1","Chicken2","Human1","Human2","Human1","Human1","Human2","Human1","Human1","Human1","Human1","Human1","Human1","Human1","Human1","Human1","Human1","Human1","Chicken1","Human1","Chicken2","Chicken1","Human1","Human1","Human2","Human1","Human1","Human1","Chicken1","Chicken1","Human1","Chicken1","Chicken1","Human1","Chicken1","Chicken1","Chicken1","Human1","Human1","Chicken1","Chicken1","Human1","Chicken2","Chicken1","Human1","Chicken1","Human1","Chicken2","Human2","Human1","Human2","Human1","Chicken2","Human2","Human1","Human1","Chicken2","Human1","Human2","Human1","Human1","Chicken1","Human1","Human2","Human1","Human1","Human1","Human1","Chicken1","Chicken1","Human1","Human1","Human1","Chicken2","Human1","Human2","Human1","Chicken1","Human2","Human2","Chicken2","Chicken2","Human2","Human2","Human2","Human2","Human2","Human2","Human2","Human2","Human2","Human2","Human2","Human2","Human2","Human2","Human2","Human2","Human2","Human2","Human2","Chicken2","Human1","Human2","Human2","Human2","Human2","Human2","Chicken1","Human2","Human2","Human2","Human2","Human2","Human2","Human2","Human2","Chicken1","Human2","Human2","Human2","Human2","Human2","Human2","Chicken2","Human2","Human2","Human2","Human2","Human2","Human2","Human2","Human2","Chicken2","Chicken2","Human2","Human2","Human2","Human2","Human2","Human2","Human2","Human2","Human2","Human2","Human2","Chicken1","Human2","Chicken1","Chicken1","Chicken1","Human2","Human2","Human2","Chicken2","Human2","Chicken1","Human2","Human1","Chicken2","Human2","Human2","Chicken1","Chicken2","Human2","Chicken1","Human2","Chicken2","Human2","Human2","Human2","Human1","Chicken1","Human2","Human2","Chicken1","Chicken1","Human2","Human2","Chicken1","Chicken2","Human2","Human2","Human2","Chicken2","Human2","Chicken2","Chicken1","Chicken1","Chicken1","Chicken1","Human2","Human2","Human2","Human2","Human2","Human2","Human2","Human2","Chicken2","Human2","Chicken1","Chicken1","Human2","Human2","Chicken1","Human2","Human2","Human2","Chicken2","Human2","Chicken2","Human2","Chicken1","Human2","Chicken2","Human2","Chicken2","Human2","Human2","Human1","Human2","Chicken2","Chicken1","Human2","Human2","Chicken2","Human1","Chicken1","Chicken1","Human2","Chicken2","Chicken1","Chicken2","Human2","Chicken1","Human2","Chicken1","Chicken1","Chicken2","Human1","Human2","Chicken1","Chicken2","Human2","Human2","Chicken2","Chicken2","Human2","Human2","Chicken2","Chicken1","Chicken1","Chicken2","Chicken1","Human1","Human1","Chicken2","Chicken1","Chicken1","Chicken2","Chicken2","Chicken1","Chicken2","Chicken2","Human1","Chicken1","Chicken1","Chicken1","Chicken1","Human2","Chicken2","Chicken1","Chicken1","Chicken1","Chicken1","Chicken1","Chicken2","Chicken1","Chicken2","Chicken1","Chicken2","Chicken1","Chicken1","Human2","Chicken2","Chicken1","Human2","Chicken1","Chicken1","Human2","Human1","Chicken1","Chicken2","Chicken2","Chicken2","Chicken1","Human1","Human1","Chicken1","Human2","Human2","Chicken1","Chicken2","Chicken1","Chicken1","Chicken1","Chicken2","Human1","Chicken1","Human2","Chicken1","Chicken1","Chicken2","Chicken2","Human2","Human1","Chicken1","Chicken1","Chicken1","Human1","Human1","Human2","Human1","Chicken1","Chicken2","Human2","Human2","Human2","Human2","Human1","Chicken1","Chicken2","Chicken1","Human2","Chicken1","Human2","Human2","Chicken2","Human1","Human1","Chicken1","Chicken1","Chicken1","Chicken2","Human2","Human2","Human2","Chicken1","Human1","Chicken2","Human2","Human1","Human1","Chicken2","Chicken2","Human1","Human2","Chicken1","Chicken1","Chicken2","Human2","Human2","Chicken2","Human2","Chicken2","Chicken1","Chicken1","Human1","Chicken1","Chicken1","Chicken2","Chicken2","Chicken2","Chicken1","Chicken2","Chicken2","Human1","Chicken2","Chicken1","Chicken1","Human1","Human1","Chicken1","Chicken1","Chicken1","Chicken1","Chicken1","Human1","Human1","Chicken2","Chicken2","Chicken1","Chicken2","Human1","Chicken2","Chicken2","Chicken1","Chicken1","Chicken1","Chicken2","Chicken2","Human2","Chicken1","Chicken1","Chicken1","Chicken1","Human1","Chicken1","Human1","Human1","Chicken1","Chicken1","Human1","Human1","Chicken2","Chicken1","Human2","Chicken2","Human1","Human1","Human1","Chicken2","Chicken1","Human1","Human1","Human2","Chicken1","Chicken1","Chicken1","Chicken1","Human1","Human1","Chicken1","Chicken1","Chicken1","Human1","Chicken2","Chicken2","Chicken2","Chicken1","Chicken2","Human1","Chicken1","Human1","Chicken1","Chicken2","Chicken1","Chicken1","Chicken2","Human1","Chicken1","Human1","Chicken1","Chicken2","Human1","Human1","Human1","Chicken2","Human1","Human1","Chicken1","Chicken2","Human1","Chicken1","Human1","Chicken2","Chicken1","Chicken1","Chicken2","Human1","Human1","Human1","Human1","Human1","Chicken1","Chicken1","Chicken1","Human1","Human1","Chicken2","Chicken1","Chicken1","Human1","Chicken1","Chicken1","Human1","Chicken2","Chicken1","Chicken1","Chicken1","Chicken1","Human1","Human1","Human1","Chicken2","Human1","Chicken1","Human1","Human1","Chicken2","Chicken1","Chicken1","Chicken2","Human1","Chicken2","Human1","Human1","Chicken2","Chicken1","Human1","Human1","Human1","Chicken2","Human1","Chicken2","Chicken1","Chicken1","Human1","Human1","Human1","Chicken2","Human1","Chicken2","Chicken1","Chicken2","Human1","Chicken1","Chicken1","Chicken2","Chicken2","Chicken2","Chicken2","Human1","Chicken1","Chicken2","Chicken2","Chicken2","Human1","Chicken2","Chicken2","Chicken1","Human1","Chicken2","Chicken2","Human1","Chicken1","Chicken2","Human1","Chicken2","Chicken2","Human1","Chicken1","Chicken2","Chicken1","Chicken2","Chicken1","Chicken1","Chicken1","Chicken1","Human1","Chicken2","Chicken1","Chicken1","Human1","Chicken1","Human2","Chicken2","Chicken2","Chicken2","Chicken2","Chicken1","Human1","Chicken1","Chicken2","Chicken1","Chicken2","Human1","Chicken1","Chicken2","Chicken2","Chicken2","Human1","Human1","Human1","Human1","Chicken1","Human1","Chicken1","Chicken1","Chicken1","Human1","Chicken1","Chicken1","Human1","Chicken1","Human1","Chicken2","Chicken2","Chicken1","Chicken2","Chicken2","Chicken2","Chicken2","Human1","Chicken1","Chicken2","Chicken1","Chicken1","Chicken2","Human1","Human1","Chicken2","Human1","Chicken1","Chicken1","Chicken2","Human1","Human1","Chicken2","Human1","Human1","Human1","Human1","Chicken1","Chicken1","Chicken2","Chicken2","Human1","Chicken1","Human1","Human1","Chicken2","Chicken2","Human1","Chicken1","Human1","Human1","Chicken2","Chicken2","Chicken2","Chicken2","Chicken2","Chicken2","Chicken2","Chicken1","Chicken2","Chicken1","Human1","Chicken2","Human1","Chicken1","Human1","Chicken2","Human1","Chicken2","Chicken2","Human1","Chicken2","Chicken2","Chicken2","Human1","Human1","Human1","Chicken1","Chicken1","Chicken2","Human1","Chicken2","Chicken2","Chicken2","Human1","Chicken2","Chicken1","Human1","Human1","Human1","Chicken2","Human1","Chicken1","Human1","Chicken2","Chicken1","Human1","Chicken1","Chicken1","Human1","Chicken2","Human1","Human1","Chicken1","Human1","Human1","Human1","Human1","Human1","Human1","Human1","Human1","Human1","Human1","Human1"],[0,0,0,0,0,0,0,0,0,0,0,0,0,1,2,4,4,5,5,5,5,6,6,6,7,7,8,8,8,10,11,14,15,16,17,18,20,22,23,23,25,25,25,26,29,31,31,32,37,37,39,40,40,41,41,47,50,67,67,67,68,69,69,70,73,78,82,86,89,89,95,100,103,110,114,115,124,128,131,137,137,139,145,153,161,164,170,174,178,196,198,217,269,385,414,442,449,465,624,683,807,811,967,1140,1397,1428,1770,2043,2062,2120,2146,2187,2360,2452,2653,2699,2885,2942,3037,3227,3530,3555,3711,3738,3827,4004,4028,4113,4189,4351,4493,4613,4703,5041,5115,5124,5131,5271,5304,5335,5552,5601,5602,5690,5720,5746,5758,5890,6089,6236,6242,6299,6508,6597,6769,6822,6836,6987,7264,7558,7563,7598,7630,7791,7804,7865,8001,8081,8148,8288,8291,8302,8374,8583,8591,8595,8596,8606,8608,8612,8648,8698,8734,8774,8828,8875,8910,8934,8953,8974,8999,9037,9137,9148,9289,9440,9467,9469,9526,9542,9573,9625,9654,9720,9736,9788,9849,9894,9899,9926,9959,9963,10028,10036,10102,10137,10245,10260,10261,10316,10343,10389,10411,10416,10505,10610,10743,10829,10911,10969,11119,11147,11151,11172,11271,11285,11286,11293,11301,11348,11386,11390,11405,11460,11486,11501,11523,11532,11549,11579,11603,11620,11652,11657,11763,11856,11887,12075,12164,12202,12347,12409,12418,12582,12624,12689,12697,12777,12794,12815,12839,12854,13057,13107,13155,13160,13299,13345,13390,13423,13435,13437,13472,13576,13710,13737,13766,13768,13775,13843,13865,13897,13925,13927,13985,14010,14052,14053,14088,14189,14190,14255,14317,14320,14336,14385,14399,14537,14568,14596,14695,14757,14762,14821,14891,14894,14959,14980,15029,15174,15210,15262,15296,15364,15380,15414,15567,15575,15579,15581,15706,15762,15764,15846,16039,16076,16084,16134,16180,16190,16197,16230,16297,16345,16350,16380,16443,16464,16465,16527,16617,16641,16708,16726,16732,16798,16830,16906,16931,16950,17000,17289,17362,17578,17578,17687,17786,17814,17858,17862,17927,18016,18055,18225,18233,18317,18325,18546,18579,18682,18683,18689,18707,18793,18794,18830,18839,18848,18875,18979,19080,19186,19218,19238,19244,19267,19284,19403,19407,19458,19506,19593,19601,19686,19704,19715,19746,19799,19829,19878,19882,19912,19937,20114,20231,20281,20303,20384,20395,20499,20502,20524,20572,20612,20685,20729,20762,20811,20830,20897,20932,20963,21015,21034,21108,21187,21195,21198,21259,21273,21367,21429,21450,21455,21578,21664,21679,21694,21694,21726,21773,21840,21903,22014,22085,22103,22162,22203,22274,22316,22433,22471,22476,22531,22648,22702,22748,22776,22794,22832,22881,22889,22966,22968,23081,23236,23306,23361,23409,23430,23454,23455,23651,23652,23668,23673,23734,23788,23898,23992,24099,24116,24279,24296,24339,24361,24380,24381,24401,24536,24568,24666,24666,24677,24687,24730,24765,24832,24858,24915,24982,24995,25035,25043,25165,25314,25317,25375,25411,25528,25566,25743,25813,25824,25855,25995,26003,26118,26178,26281,26381,26446,26447,26449,26488,26499,26520,26625,26675,26814,26814,26825,26905,26968,26990,27026,27037,27050,27142,27171,27279,27431,27440,27453,27516,27549,27554,27562,27700,27732,27793,27843,27869,27896,27915,27947,28030,28086,28185,28191,28244,28291,28321,28392,28605,28617,28694,28697,28760,28770,28797,28808,28844,28947,28961,29079,29125,29145,29146,29225,29345,29390,29425,29590,29653,29907,29907,30114,30184,30256,30339,30408,30439,30446,30485,30510,30577,30629,30642,30858,30955,31036,31111,31200,31225,31293,31326,31517,31718,31875,31908,31977,31989,32142,32147,32173,32191,32469,32534,32613,32730,32958,32998,33009,33039,33131,33137,33177,33197,33414,33494,33501,33563,33606,33662,33673,33683,33754,33948,33973,34001,34069,34103,34154,34283,34361,34383,34452,34515,34738,34965,35029,35161,35214,35279,35419,35527,35541,36043,36131,36143,36583,36636,36812,36829,37048,37263,37387,37477,37552,37812,37879,37982,38035,38547,38602,38919,39174,39478,39518,39627,39816,39889,40608,41007,41151,41383,41421,42621,43153,43340,43601,43615,45298,46983,50174,50961,51543,54328,55364,62379,64294,76823,88037,108042,124468,150837,184708,191176,203803,233850,235165,254676,698347]],"container":"<table class=\"display\">\n  <thead>\n    <tr>\n      <th> <\/th>\n      <th>SampleID<\/th>\n      <th>Experiment<\/th>\n      <th>Model2<\/th>\n      <th>LibrarySize<\/th>\n    <\/tr>\n  <\/thead>\n<\/table>","options":{"columnDefs":[{"className":"dt-right","targets":4},{"orderable":false,"targets":0}],"order":[],"autoWidth":false,"orderClasses":false}},"evals":[],"jsHooks":[]}</script>
```

```r
out$plot +
  facet_wrap(~ Model2, scales = "free_x") +
  theme_linedraw() + scale_y_log10()
```

```
## Warning: Transformation introduced infinite values in continuous y-axis
## Transformation introduced infinite values in continuous y-axis
```

```
## Warning: ggrepel: 17 unlabeled data points (too many overlaps). Consider
## increasing max.overlaps
```

```
## Warning: ggrepel: 54 unlabeled data points (too many overlaps). Consider
## increasing max.overlaps
```

```
## Warning: ggrepel: 5 unlabeled data points (too many overlaps). Consider
## increasing max.overlaps
```

```
## Warning: ggrepel: 42 unlabeled data points (too many overlaps). Consider
## increasing max.overlaps
```

<img src="16S_validation_files/figure-html/unnamed-chunk-7-1.png" style="display: block; margin: auto;" />

```r
  coord_cartesian(ylim = c(0, 5000), xlim = c(0,600))
```

```
## <ggproto object: Class CoordCartesian, Coord, gg>
##     aspect: function
##     backtransform_range: function
##     clip: on
##     default: FALSE
##     distance: function
##     expand: TRUE
##     is_free: function
##     is_linear: function
##     labels: function
##     limits: list
##     modify_scales: function
##     range: function
##     render_axis_h: function
##     render_axis_v: function
##     render_bg: function
##     render_fg: function
##     setup_data: function
##     setup_layout: function
##     setup_panel_guides: function
##     setup_panel_params: function
##     setup_params: function
##     train_panel_guides: function
##     transform: function
##     super:  <ggproto object: Class CoordCartesian, Coord, gg>
```


```r
out$df %>% 
  filter(Experiment %!in% c("Mock", "NTC")) %>% 
  filter(LibrarySize < 1000) %>% 
  nrow()
```

```
## [1] 62
```

62 samples lost 



```r
out$df %>% 
  filter(Experiment %!in% c("Mock", "NTC")) %>% 
  filter(LibrarySize < 4000) %>% 
  nrow()
```

```
## [1] 75
```


```r
require(parallel)
```

```
## Loading required package: parallel
```

```r
ps_16S %>%
  subset_samples(Experiment %!in% c("Mock", "NTC")) %>% 
  prune_samples(sample_sums(.)>= 50, .) %>% 
  # filter_taxa(function(x) sum(x > 0) > 0, TRUE) %>% 
  ggrare(step = 50, parallel = TRUE,  se = FALSE, color = "Model2", plot = FALSE) -> rare_curves
```


```r
rare_curves +
  theme_classic() +
  facet_wrap( . ~ Model2) +
  # geom_vline(xintercept = 500,
  #            color="red",
  #            linetype="dashed", size=0.25) +
  ylab("ASV Richness") -> plot

plot %>% 
  ggpubr::change_palette(p = ., palette = "lancet") -> plot
```


```r
plot + xlim(c(0, 50000))
```

```
## Warning: Removed 39807 row(s) containing missing values (geom_path).
```

<img src="16S_validation_files/figure-html/unnamed-chunk-12-1.png" style="display: block; margin: auto;" />



```r
plot + xlim(c(0, 4000))
```

```
## Warning: Removed 246354 row(s) containing missing values (geom_path).
```

<img src="16S_validation_files/figure-html/unnamed-chunk-13-1.png" style="display: block; margin: auto;" />

Based on those plots I would suggest :

```r
depth = 4000

ps_16S %>%
  rarefy_even_depth(sample.size = depth,
                    rngseed = 123) -> physeq_rare
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
## 125 samples removedbecause they contained fewer reads than `sample.size`.
```

```
## Up to first five removed samples are:
```

```
## CR-40-S166CR-52-S196EMPTY-S384IR1-69-S198Positive2-S163
```

```
## ...
```

```
## 271OTUs were removed because they are no longer 
## present in any sample after random subsampling
```

```
## ...
```

```r
ps_16S %>%
  prune_samples(sample_sums(.)>= depth, .) -> physeq_fil
```



```r
# Alpha-diversity:
plot_alpha_div_NRP72 <- function(df, x = "Day_of_Treatment", y = "value",  point_size = 2.5, color = "Treatment", fill  = "Treatment", shape = "Antibiotic", alpha = "Antibiotic_mg.L", facet_formula = paste0("alphadiversiy ~  Model2 "), ylab = "Resistome alpha-diversity", xlab = "Days (Treatment)", measures = c("Observed"), path_group = "interaction(Model, Fermentation, Reactor_Treatment)"){
  
  df %>% 
    pivot_longer(cols = all_of(measures), values_to = "value", names_to = 'alphadiversiy', values_drop_na  = TRUE) %>%
    mutate(alphadiversiy = fct_relevel(alphadiversiy, measures)) -> df_ready
  
  df_ready %>% 
    ggplot(aes_string(x = x, y = y, color = color, fill = fill, shape = shape, alpha = alpha)) + #shape = Fermentation
    geom_point(size = point_size) + 
    geom_line(linetype = 2,  size = 0.5,  aes_string(group = path_group)) +
    labs(y = ylab, x = xlab) +
    scale_color_manual(name = "", values = treat_col,
                       na.value = "black") +
    scale_fill_manual(name = "", values = treat_col,
                      na.value = "black") +
    facet_grid(as.formula(facet_formula), scales = "free", space = "fixed", drop = TRUE) +
    scale_shape_manual(name = "" ,values = antibio_shape, na.value =  17) +
    scale_alpha_discrete(name = "", range=c(0.6, 1), na.value =  0.6) + 
    geom_vline(xintercept = 0, size = 0.5, linetype = 1, color = "grey60") -> p_alpha_both
  
  return(p_alpha_both)
}
```




```r
physeq_rare %>% 
  phyloseq_alphas(phylo = FALSE) -> alpha_16S

alpha_16S %>% 
  group_by(Model, Fermentation, Antibiotic, Antibiotic_mg.mL, Treatment, Phase, Period) %>% 
  rstatix::get_summary_stats(Observed, type = "common") %>% 
  DT::datatable()
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
## [1] parallel  stats     graphics  grDevices utils     datasets  methods  
## [8] base     
## 
## other attached packages:
##  [1] reshape2_1.4.4       scales_1.2.0         compositions_2.0-4  
##  [4] nlme_3.1-158         randomcoloR_1.1.0.1  here_1.0.1          
##  [7] ggrepel_0.9.1        speedyseq_0.5.3.9018 phyloseq_1.39.1     
## [10] forcats_0.5.1        stringr_1.4.0        dplyr_1.0.9         
## [13] purrr_0.3.4          readr_2.1.2          tidyr_1.2.0         
## [16] tibble_3.1.8         ggplot2_3.3.6        tidyverse_1.3.1.9000
## 
## loaded via a namespace (and not attached):
##   [1] googledrive_2.0.0      Rtsne_0.16             colorspace_2.0-3      
##   [4] ggsignif_0.6.3         ellipsis_0.3.2         rprojroot_2.0.3       
##   [7] XVector_0.36.0         fs_1.5.2               rstudioapi_0.13       
##  [10] ggpubr_0.4.0           farver_2.1.1           DT_0.23               
##  [13] fansi_1.0.3            lubridate_1.8.0        xml2_1.3.3            
##  [16] codetools_0.2-18       splines_4.2.1          robustbase_0.95-0     
##  [19] knitr_1.39             ade4_1.7-19            jsonlite_1.8.0        
##  [22] broom_1.0.0            cluster_2.1.3          dbplyr_2.2.1          
##  [25] compiler_4.2.1         httr_1.4.3             backports_1.4.1       
##  [28] assertthat_0.2.1       Matrix_1.4-1           fastmap_1.1.0         
##  [31] gargle_1.2.0           cli_3.3.0              htmltools_0.5.3       
##  [34] tools_4.2.1            igraph_1.3.4           gtable_0.3.0          
##  [37] glue_1.6.2             GenomeInfoDbData_1.2.8 V8_4.2.0              
##  [40] Rcpp_1.0.9             carData_3.0-5          Biobase_2.56.0        
##  [43] cellranger_1.1.0       jquerylib_0.1.4        vctrs_0.4.1           
##  [46] Biostrings_2.64.0      rhdf5filters_1.8.0     multtest_2.52.0       
##  [49] ape_5.6-2              crosstalk_1.2.0        iterators_1.0.14      
##  [52] tensorA_0.36.2         xfun_0.31              rvest_1.0.2           
##  [55] lifecycle_1.0.1        rstatix_0.7.0          googlesheets4_1.0.0   
##  [58] DEoptimR_1.0-11        zlibbioc_1.42.0        MASS_7.3-58           
##  [61] hms_1.1.1              biomformat_1.24.0      rhdf5_2.40.0          
##  [64] yaml_2.3.5             curl_4.3.2             dtplyr_1.2.1          
##  [67] sass_0.4.1             stringi_1.7.8          highr_0.9             
##  [70] S4Vectors_0.34.0       foreach_1.5.2          permute_0.9-7         
##  [73] BiocGenerics_0.42.0    GenomeInfoDb_1.32.2    rlang_1.0.4           
##  [76] pkgconfig_2.0.3        bitops_1.0-7           evaluate_0.15         
##  [79] lattice_0.20-45        Rhdf5lib_1.18.2        labeling_0.4.2        
##  [82] htmlwidgets_1.5.4      tidyselect_1.1.2       ggsci_2.9             
##  [85] plyr_1.8.7             magrittr_2.0.3         R6_2.5.1              
##  [88] IRanges_2.30.0         generics_0.1.3         DBI_1.1.3             
##  [91] pillar_1.8.0           haven_2.5.0            withr_2.5.0           
##  [94] mgcv_1.8-40            abind_1.4-5            survival_3.3-1        
##  [97] RCurl_1.98-1.8         bayesm_3.1-4           car_3.1-0             
## [100] modelr_0.1.8           crayon_1.5.1           utf8_1.2.2            
## [103] tzdb_0.3.0             rmarkdown_2.14         grid_4.2.1            
## [106] readxl_1.4.0           data.table_1.14.2      vegan_2.6-2           
## [109] reprex_2.0.1           digest_0.6.29          stats4_4.2.1          
## [112] munsell_0.5.0          bslib_0.3.1
```

