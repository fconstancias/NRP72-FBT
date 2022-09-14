---
title: "NTP72 - 16S alpha div"
author: "Florentin Constancias"
date: "September 14, 2022"
output: 
  html_document: 
    toc: yes
    keep_md: yes
---



# Load required packages



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
# "Figures/Rastetics.Rdata" %>%
#   here::here() %>% 
#   load()

load(here::here("Figures/Rastetics.Rdata"))
```

## Data

### 16S:

```r
("data/processed/16S/ps_silva_dada2_human_chicken_all_fct_polyferms_090922_up.RDS") %>% 
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
## [13] "WEEKS_DG"                    "DAYS_DG"                    
## [15] "Day_from_Inoculum"           "Enrichment"                 
## [17] "Phase"                       "Treatment2"                 
## [19] "Date"                        "Paul"                       
## [21] "Reactor_Treatment"           "Reactor_F1_F2"              
## [23] "Model"                       "Model2"                     
## [25] "GeneCopyNumberperML"         "HV292.1_Copy_Number_permL"  
## [27] "CCUG59168_Copy_Number_permL" "CTX_Copy_Number_permL"      
## [29] "VAN_Copy_Number_permL"       "Antibiotic_mg.L"            
## [31] "Fermentation"                "Antibiotic"                 
## [33] "Lactose_mM"                  "Glucose_mM"                 
## [35] "Galactose_mM"                "Succinat_mM"                
## [37] "Lactat_mM"                   "Formiat_mM"                 
## [39] "Acetat_mM"                   "Propionat_mM"               
## [41] "Isobutyrat_mM"               "Butyrat_mM"                 
## [43] "Isovalerat_mM"               "Valerat_mM"                 
## [45] "Total_SCFA_mM"               "metagenomic_sample_name"    
## [47] "raw_metagenomic_pairs"       "Period"                     
## [49] "Reactor_Treatment_Dose"      "Treatment_Dose"
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
<div id="htmlwidget-3a33a1c56abff69ec01d" style="width:100%;height:auto;" class="datatables html-widget"></div>
<script type="application/json" data-for="htmlwidget-3a33a1c56abff69ec01d">{"x":{"filter":"none","vertical":false,"data":[["345","203","215","221","237","316","342","361","238","257","479","341","322","363","251","310","217","228","186","343","346","227","318","202","239","371","679","20","375","294","218","241","282","675","52","379","45","119","680","234","673","85","677","104","105","262","211","672","51","160","301","682","323","678","674","367","684","671","374","333","325","213","485","219","311","182","676","338","683","300","351","681","288","350","18","248","687","376","196","670","271","669","283","94","686","685","668","430","308","624","168","10","517","582","625","469","482","108","35","30","666","492","55","313","389","174","425","48","516","284","601","93","81","54","77","527","626","491","509","9","23","46","22","500","644","149","155","589","588","580","61","581","486","487","39","519","512","352","464","171","506","615","261","628","63","133","501","608","579","590","13","526","344","32","521","510","511","634","488","125","423","28","296","303","502","128","178","495","496","503","507","377","67","148","89","515","97","110","147","176","16","153","483","392","36","124","407","175","79","157","244","538","42","513","473","553","83","188","216","134","24","499","86","123","76","518","380","40","14","599","53","504","520","667","354","15","117","135","320","286","250","387","5","331","37","388","95","8","606","470","651","610","368","360","7","167","585","80","607","642","443","114","291","497","231","353","514","594","440","564","204","151","3","393","557","623","383","11","57","307","574","164","107","493","665","524","152","472","494","245","476","573","72","121","184","192","17","172","44","115","100","280","327","505","410","87","428","266","618","663","475","484","64","162","91","498","33","116","597","71","92","551","312","122","267","617","335","156","571","150","270","289","412","533","158","385","662","317","348","222","525","56","209","191","126","34","604","101","661","145","321","293","29","4","593","161","357","413","536","457","129","547","414","614","232","539","587","329","146","409","619","90","142","384","194","19","355","545","570","180","583","556","411","253","285","185","417","315","629","236","59","456","309","541","166","242","429","654","131","643","163","273","275","205","562","370","190","109","136","559","88","295","256","454","47","84","195","591","1","41","365","480","173","138","144","179","611","366","306","373","620","450","259","584","27","254","198","390","103","12","447","243","433","224","299","403","113","249","347","537","378","422","336","636","438","120","165","183","197","659","278","446","252","382","140","432","230","111","118","449","437","434","523","292","38","386","664","635","566","415","481","326","649","421","522","65","201","439","419","247","657","43","458","334","633","405","406","660","337","600","102","453","49","426","82","127","31","98","233","647","418","60","73","226","25","455","452","435","632","78","260","576","26","444","567","631","572","75","397","540","263","542","565","69","490","448","431","279","223","206","575","358","70","200","569","66","99","638","246","137","21","369","616","154","268","653","401","399","130","474","596","427","471","416","220","598","6","652","132","143","641","460","340","302","398","328","58","68","489","298","277","395","645","193","314","304","349","141","2","424","461","359","640","265","297","532","467","592","451","258","74","225","613","305","445","465","408","646","420","648","402","462","159","463","609","332","529","466","372","555","612","50","212","477","650","208","442","548","391","181","605","381","404","549","478","240","534","287","339","112","546","550","170","468","210","508","400","396","543","603","274","436","558","169","290","356","639","189","552","459","324","622","106","272","655","528","621","568","441","578","577","177","595","531","630","214","602","62","544","96","586","561","560","535","281","627","394","530","637","554","207","563","199","658","656","139","264","187","319","269","330","235","364","229","362","255","276"],["TR5-68-S99","IR-74-S188","Negative-2-S104","TR1-19-S116","TR1-63-S100","TR4-65-S155","TR5-65-S187","TR6-36-S124","TR1-64-S120","TR2-35-S131","p-2-H8-S188","TR5-64-S108","TR5-01-S123","TR6-60-S119","TR2-29-S117","TR4-59-S180","Negative-7-S223","TR1-32-S111","IR-13-S103","TR5-66-S101","TR5-69-S107","TR1-31-S109","TR4-67-S192","IR-51-S157","TR1-65-S128","TR6-68-S127","NC-2-S98","CR-52-S196","negative-4-S189","TR3-68-S194","Negative-8-S224","TR1-67-S133","TR3-34-S149","NC-16-S290","TR1-16-S199","positive-2-S181","IR1-69-S198","TR4-18-S200","NC-20-S335","TR1-60-S114","NC-14-S288","TR2-34-S302","NC-18-S333","TR3-28-S162","TR3-3-S179","TR2-62-S97","IR-82-S143","NC-13-S287","TR1-15-S168","TR6-13-S164","TR4-28-S138","NC-4-S100","TR5-13-S137","NC-19-S334","NC-15-S289","TR6-64-S105","NC-6-S102","NC-12-S286","negative-1-S58","TR5-34-S28","TR5-26-S24","IR-84-S178","p-3-A5-S197","TR1-01-S102","TR4-60-S36","D-1-S49","NC-17-S332","TR5-61-S142","NC-5-S101","TR4-27-S130","TR6-26-S122","NC-3-S99","TR3-62-S37","TR6-19-S33","CR-40-S166","TR2-26-S87","NC-9-S190","negative-5-S191","IR-46-S126","NC-11-S192","TR3-01-S121","NC-10-S191","TR3-35-S20","TR3-14-S136","NC-8-S140","NC-7-S139","NC-1-S97","p-2-D11-S143","TR4-35-S202","AP-67-S67","TR6-21-S107","CR-21-S126","p-3-G1-S265","AP-29-S29","AP-68-S68","p-2-G9-S177","p-3-A2-S194","TR3-37-S221","IR1-39-S351","IR1-3-S320","MOCK-8-S330","p-3-B6-S210","TR1-19-S305","TR4-62-S1","p-1-F11-S71","TR6-34-S143","p-2-C7-S127","TR1-1-S110","p-3-F6-S258","TR3-36-S8","AP-46-S46","TR3-13-S304","TR2-25-S133","TR1-18-S124","TR2-19-S318","p-3-H5-S281","AP-69-S69","p-3-B5-S209","p-3-E5-S245","CR-20-S111","D-0-S154","IR1-77-S202","CR-64-S324","p-3-D2-S230","AP-86-S86","TR5-34-S127","TR5-57-S348","AP-35-S35","AP-34-S34","AP-27-S27","TR1-3-S121","AP-28-S28","p-3-A6-S198","p-3-B1-S205","IR1-45-S128","p-3-G3-S267","p-3-F2-S254","TR6-27-S6","p-2-G4-S172","TR6-28-S158","p-3-E2-S242","AP-59-S59","TR2-61-S196","AP-70-S70","TR1-34-S317","TR4-58-S112","p-3-D3-S231","AP-52-S52","AP-26-S26","AP-36-S36","CR-28-S106","p-3-H4-S280","TR5-67-S96","IR1-36-S159","p-3-G5-S269","p-3-E6-S246","p-3-F1-S253","AP-76-S76","p-3-B2-S206","TR4-28-S183","p-2-C5-S125","IR1-24-S295","TR3-7-S204","TR4-30-S48","p-3-D4-S232","TR4-34-S342","TR6-52-S122","p-3-C4-S220","p-3-C5-S221","p-3-D5-S233","p-3-E3-S243","negative-6-S213","TR1-52-S253","TR5-31-S313","TR2-52-S230","p-3-F5-S257","TR3-17-S222","TR3-46-S224","TR5-3-S331","TR6-40-S160","CR-34-S359","TR5-52-S93","p-3-A3-S195","p-1-H10-S94","IR1-4-S223","TR4-25-S333","p-2-B11-S119","TR6-37-S312","TR2-21-S228","TR5-59-S374","TR1-7-S153","AP-108-S114","IR1-57-S356","p-3-F3-S255","p-2-H2-S182","AP-121-S127","TR2-3-S297","IR-22-S41","Negative-3-S150","TR4-64-S220","IR1-16-S307","p-3-D1-S229","TR2-37-S330","TR4-22-S344","TR2-18-S326","p-3-G2-S266","positive-3-S212","IR1-48-S339","CR-3-S98","AP-44-S44","TR1-17-S246","p-3-D6-S234","p-3-G4-S268","MOCK-9-S331","TR6-29-S207","CR-31-S334","TR4-16-S276","TR5-1-S210","TR4-69-S197","TR3-60-S46","TR2-28-S64","p-1-E11-S59","CR-16-S332","TR5-32-S56","IR1-43-S256","p-1-F10-S70","TR3-15-S338","CR-19-S341","AP-50-S50","p-2-H1-S181","AP-92-S92","AP-54-S54","TR6-65-S94","TR6-35-S208","CR-18-S252","TR6-20-S132","AP-31-S31","TR2-22-S364","AP-51-S51","AP-84-S84","p-2-E5-S149","TR4-13-S142","TR3-65-S92","p-3-C6-S222","TR1-35-S161","TR6-28-S80","p-3-F4-S256","AP-4-S4","p-2-E2-S146","AP-131-S137","IR-75-S44","TR5-40-S335","CR-14-S336","p-1-H11-S95","AP-125-S131","AP-66-S66","p-1-B11-S23","CR-22-S311","TR1-21-S207","TR4-34-S16","AP-21-S21","TR6-17-S138","TR3-34-S354","p-3-C2-S218","MOCK-7-S329","p-3-H2-S278","TR5-46-S310","p-2-H11-S191","p-3-C3-S219","TR2-01-S14","p-2-H5-S185","AP-20-S20","TR2-14-S232","TR4-20-S269","HV292-1-S205","IR-42-S9","CR-37-S367","TR6-3-S236","IR1-63-S155","TR4-14-S226","TR3-20-S265","TR3-32-S211","TR5-28-S5","p-3-E1-S241","p-2-B3-S111","TR2-40-S259","p-2-D1-S133","TR2-66-S93","AP-61-S61","MOCK-5-S284","p-2-H4-S184","p-3-A4-S196","TR1-37-S306","TR6-15-S149","TR2-64-S216","p-3-C7-S223","IR1-37-S343","TR4-15-S218","AP-42-S42","TR2-13-S231","TR3-1-S237","AP-12-S12","TR4-61-S11","TR4-21-S254","TR2-67-S193","AP-60-S60","TR5-36-S85","TR5-58-S372","AP-19-S19","TR5-37-S368","TR2-7-S215","TR3-63-S206","p-2-B5-S113","AP-103-S109","TR5-64-S141","p-1-C11-S35","MOCK-4-S283","TR4-66-S55","TR6-01-S3","TR1-26-S27","p-3-H3-S279","TR1-20-S118","IR-80-S67","IR-41-S106","TR4-3-S309","IR1-38-S337","AP-49-S49","TR3-21-S234","MOCK-3-S189","TR5-25-S240","TR4-7-S65","TR3-67-S199","IR1-26-S119","CR-15-S292","AP-39-S39","TR6-14-S211","TR6-32-S39","p-2-B6-S114","AP-106-S112","p-2-F8-S164","TR4-37-S188","AP-116-S122","p-2-B7-S115","AP-58-S58","TR1-36-S198","AP-109-S115","AP-33-S33","TR5-30-S62","TR5-28-S214","p-2-B2-S110","AP-62-S62","TR2-58-S172","TR5-20-S208","p-1-B12-S24","IR-44-S89","CR-46-S191","TR6-30-S144","AP-114-S120","AP-18-S18","TR6-64-S146","AP-3-S3","AP-124-S130","p-2-B4-S112","TR2-31-S32","TR3-59-S61","IR-01-S201","p-2-C1-S121","TR4-64-S10","AP-71-S71","TR1-62-S7","TR1-25-S117","p-2-F7-S163","TR4-36-S12","AP-110-S116","TR6-19-S137","TR1-68-S13","p-2-D10-S142","AP-95-S95","TR4-46-S215","AP-85-S85","TR6-16-S371","TR3-19-S214","TR3-27-S203","IR-76-S2","AP-13-S13","TR6-67-S140","IR-34-S4","TR3-40-S257","TR5-13-S103","AP-127-S133","TR2-46-S239","TR3-69-S139","TR2-34-S71","p-2-F5-S161","IR1-8-S260","TR2-31-S182","IR-45-S216","AP-37-S37","CR-1-S116","IR1-54-S130","TR6-62-S50","p-2-H9-S189","TR6-31-S272","TR5-16-S290","TR5-22-S145","TR6-58-S361","AP-55-S55","TR6-63-S195","TR4-33-S84","TR6-7-S91","AP-63-S63","p-2-F11-S167","TR2-59-S60","AP-30-S30","IR1-2-S96","TR2-32-S90","IR-48-S98","p-1-G10-S82","TR3-25-S102","CR-25-S383","p-2-E9-S153","TR1-69-S145","p-2-D4-S136","TR1-28-S35","TR4-26-S26","p-2-A8-S104","TR3-64-S190","TR2-27-S59","TR5-7-S66","AP-107-S113","positive-1-S25","p-2-C4-S124","TR5-59-S31","AP-79-S79","p-2-D9-S141","TR4-19-S270","TR6-18-S262","D-2-S164","IR-47-S63","MOCK-1-S187","TR3-30-S190","p-2-E8-S152","TR2-30-S73","p-1-A12-S12","TR5-18-S201","p-2-D3-S135","TR1-34-S34","TR3-52-S217","TR4-17-S114","p-2-F10-S166","p-2-D8-S140","p-2-D5-S137","p-3-H1-S277","TR3-66-S69","IR1-44-S264","p-1-D11-S47","MOCK-6-S285","AP-77-S77","AP-14-S14","p-2-B8-S116","p-3-A1-S193","TR5-27-S148","AP-90-S90","p-2-C3-S123","p-3-G6-S270","TR1-40-S209","IR-50-S81","p-2-E1-S145","p-2-C11-S131","TR2-19-S40","AP-98-S104","IR1-6-S250","p-2-F9-S165","TR5-35-S77","AP-75-S75","p-2-B1-S109","p-2-B10-S118","MOCK-2-S188","TR5-60-S75","AP-45-S45","TR3-22-S235","p-2-F4-S160","TR1-13-S171","p-2-C8-S128","TR2-28-S286","TR4-31-S278","IR1-35-S370","TR3-18-S293","TR1-59-S42","AP-89-S89","p-2-C10-S130","TR1-28-S289","TR2-15-S288","TR1-30-S76","IR1-17-S279","p-2-F6-S162","p-2-F3-S159","p-2-D6-S138","AP-74-S74","TR2-20-S271","TR2-60-S79","AP-23-S23","IR1-18-S287","p-2-E6-S150","AP-15-S15","AP-73-S73","AP-2-S2","TR2-17-S248","p-2-A12-S108","AP-11-S11","TR2-63-S160","AP-111-S117","AP-132-S138","TR1-64-S263","p-3-B4-S208","p-2-F1-S157","p-2-D2-S134","TR3-31-S23","TR1-27-S54","IR-77-S163","AP-22-S22","TR6-33-S135","TR2-1-S273","IR-5-S70","AP-17-S17","TR1-46-S170","TR3-19-S266","AP-80-S80","TR2-13-S156","TR5-14-S267","CR-58-S382","TR6-66-S132","AP-6-S6","TR5-56-S241","TR2-68-S147","AP-94-S94","p-2-A6-S102","p-2-A3-S99","TR4-40-S242","p-2-H3-S183","AP-41-S41","p-2-C9-S129","p-2-H10-S190","p-2-B9-S117","TR1-13-S22","AP-43-S43","CR-17-S97","AP-93-S93","TR4-52-S95","TR5-21-S185","AP-83-S83","p-2-G10-S178","TR5-63-S82","TR4-29-S159","p-2-A2-S98","TR5-29-S136","TR1-22-S282","TR1-58-S275","p-3-B3-S207","TR4-19-S154","TR3-29-S210","p-2-A10-S106","AP-87-S87","IR-43-S162","TR4-63-S74","TR4-31-S19","TR6-13-S21","TR5-19-S87","CR-13-S148","p-2-C6-S126","p-2-G11-S179","TR6-34-S141","AP-82-S82","TR2-65-S158","TR4-13-S168","AP-102-S108","p-2-G7-S175","AP-38-S38","p-2-F2-S158","TR2-36-S183","TR2-16-S184","TR1-29-S18","AP-57-S57","TR4-32-S17","p-2-E7-S151","p-2-G5-S173","p-2-B12-S120","AP-88-S88","p-2-C2-S122","AP-9-S9","p-2-A7-S103","p-2-G2-S170","TR6-1-S247","p-2-G3-S171","AP-53-S53","TR5-33-S134","AP-10-S10","p-2-G6-S174","TR6-69-S174","AP-123-S129","AP-56-S56","TR1-14-S274","IR-83-S170","p-2-H6-S186","AP-91-S91","IR-79-S167","p-2-E4-S148","AP-117-S123","p-1-G11-S83","CCUG59168-S78","AP-5-S5","p-1-A11-S11","p-2-A9-S105","AP-118-S124","p-2-H7-S187","TR1-66-S151","AP-104-S110","TR3-61-S45","TR5-62-S176","TR3-58-S233","AP-115-S121","AP-119-S125","TR6-25-S378","p-2-G8-S176","IR-81-S182","p-3-E4-S244","p-2-A4-S100","p-2-A11-S107","AP-112-S118","AP-48-S48","TR3-26-S175","p-2-D7-S139","AP-126-S132","TR6-22-S152","TR3-64-S146","TR6-31-S179","AP-81-S81","IR-28-S184","AP-120-S126","p-2-G1-S169","TR5-19-S173","AP-65-S65","TR3-31-S283","TR3-13-S129","AP-96-S96","AP-1-S1","AP-64-S64","AP-16-S16","p-2-E3-S147","AP-25-S25","AP-24-S24","TR6-46-S176","AP-40-S40","AP-101-S107","AP-72-S72","IR-9-S165","AP-47-S47","TR1-31-S153","AP-113-S119","TR3-16-S296","AP-32-S32","AP-129-S135","AP-128-S134","AP-105-S111","TR3-33-S86","AP-7-S7","p-2-A1-S97","AP-100-S106","AP-8-S8","AP-122-S128","IR-78-S29","AP-130-S136","IR-49-S125","AP-99-S105","AP-97-S103","TR5-17-S82","TR2-64-S51","IR-16-S113","TR4-68-S52","TR2-69-S57","TR5-31-S88","TR1-61-S47","TR6-61-S186","TR1-33-S177","TR6-59-S185","TR2-33-S30","TR3-28-S38"],["Continuous","Continuous","NTC","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","NTC","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","NTC","Continuous","NTC","Continuous","NTC","Continuous","Continuous","NTC","Continuous","Mock","Continuous","Continuous","NTC","Continuous","NTC","Continuous","NTC","Continuous","Continuous","Continuous","Continuous","NTC","Continuous","Continuous","Continuous","NTC","Continuous","NTC","NTC","Continuous","NTC","NTC","NTC","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Cecum","NTC","Continuous","NTC","Continuous","Continuous","NTC","Continuous","Continuous","Continuous","Continuous","NTC","NTC","Continuous","NTC","Continuous","NTC","Continuous","Continuous","NTC","NTC","NTC","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Mock","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","NTC","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","NTC","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Mock","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Mock","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Mock","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","HV292.1","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Mock","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Mock","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Mock","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Mock","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Cecum","Continuous","Mock","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Mock","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Mock","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","CCUG59168","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous","Continuous"],["Human1","Human1","Human1","Human1","Human1","Human1","Human1","Human1","Human1","Human1","Chicken2","Human1","Human1","Human1","Human1","Human1","Human1","Human1","Human1","Human1","Human1","Human1","Human1","Human1","Human1","Human1","Human2","Chicken1","Human1","Human1","Human1","Human1","Human1","Human2","Chicken1","Human1","Chicken1","Chicken1","Human2","Human1","Human2","Chicken1","Human2","Chicken1","Chicken1","Human1","Human1","Human2","Chicken1","Chicken1","Human1","Human2","Human1","Human2","Human2","Human1","Human2","Human2","Human1","Human1","Human1","Human1","Chicken2","Human1","Human1","Human1","Human2","Human1","Human2","Human1","Human1","Human2","Human1","Human1","Chicken1","Human1","Human2","Human1","Human1","Human2","Human1","Human2","Human1","Chicken1","Human2","Human2","Human2","Chicken2","Human1","Human2","Chicken1","Chicken1","Chicken2","Human2","Human2","Chicken2","Chicken2","Chicken1","Chicken1","Chicken1","Human2","Chicken2","Chicken1","Human1","Chicken2","Chicken1","Chicken2","Chicken1","Chicken2","Human1","Human2","Chicken1","Chicken1","Chicken1","Chicken1","Chicken2","Human2","Chicken2","Chicken2","Chicken1","Chicken1","Chicken1","Chicken1","Chicken2","Human2","Chicken1","Chicken1","Human2","Human2","Human2","Chicken1","Human2","Chicken2","Chicken2","Chicken1","Chicken2","Chicken2","Human1","Chicken2","Chicken1","Chicken2","Human2","Human1","Human2","Chicken1","Chicken1","Chicken2","Human2","Human2","Human2","Chicken1","Chicken2","Human1","Chicken1","Chicken2","Chicken2","Chicken2","Human2","Chicken2","Chicken1","Chicken2","Chicken1","Human1","Human1","Chicken2","Chicken1","Chicken1","Chicken2","Chicken2","Chicken2","Chicken2","Human1","Chicken1","Chicken1","Chicken1","Chicken2","Chicken1","Chicken1","Chicken1","Chicken1","Chicken1","Chicken1","Chicken2","Chicken2","Chicken1","Chicken1","Chicken2","Chicken1","Chicken1","Chicken1","Human1","Human2","Chicken1","Chicken2","Chicken2","Human2","Chicken1","Human1","Human1","Chicken1","Chicken1","Chicken2","Chicken1","Chicken1","Chicken1","Chicken2","Human1","Chicken1","Chicken1","Human2","Chicken1","Chicken2","Chicken2","Human2","Human1","Chicken1","Chicken1","Chicken1","Human1","Human1","Human1","Chicken2","Chicken1","Human1","Chicken1","Chicken2","Chicken1","Chicken1","Human2","Chicken2","Human2","Human2","Human1","Human1","Chicken1","Chicken1","Human2","Chicken1","Human2","Human2","Chicken2","Chicken1","Human1","Chicken2","Human1","Human1","Chicken2","Human2","Chicken2","Human2","Human1","Chicken1","Chicken1","Chicken2","Human2","Human2","Chicken2","Chicken1","Chicken1","Human1","Human2","Chicken1","Chicken1","Chicken2","Human2","Chicken2","Chicken1","Chicken2","Chicken2","Human1","Chicken2","Human2","Chicken1","Chicken1","Human1","Human1","Chicken1","Chicken1","Chicken1","Chicken1","Chicken1","Human1","Human1","Chicken2","Chicken2","Chicken1","Chicken2","Human1","Human2","Human2","Chicken2","Chicken2","Chicken1","Chicken1","Chicken1","Chicken2","Chicken1","Chicken1","Human2","Chicken1","Chicken1","Human2","Human1","Chicken1","Human1","Human2","Human1","Chicken1","Human2","Chicken1","Human1","Human1","Chicken2","Human2","Chicken1","Chicken2","Human2","Human1","Human1","Human1","Chicken2","Chicken1","Human1","Human1","Chicken1","Chicken1","Human2","Chicken1","Human2","Chicken1","Human1","Human1","Chicken1","Chicken1","Human2","Chicken1","Human1","Chicken2","Human2","Chicken2","Chicken1","Human2","Chicken2","Human2","Human1","Human2","Human2","Human1","Chicken1","Chicken2","Human2","Chicken1","Chicken1","Chicken2","Human1","Chicken1","Human1","Human2","Human2","Chicken1","Human2","Human2","Chicken2","Human1","Human1","Human1","Chicken2","Human1","Human2","Human1","Chicken1","Chicken2","Human1","Human2","Chicken1","Human1","Chicken2","Human2","Chicken1","Human2","Chicken1","Human1","Human1","Human1","Human2","Human1","Human1","Chicken1","Chicken1","Human2","Chicken1","Human1","Human1","Chicken2","Chicken1","Chicken1","Human1","Human2","Chicken1","Chicken1","Human1","Chicken2","Chicken1","Chicken1","Chicken1","Chicken1","Human2","Human1","Human1","Human1","Human2","Chicken2","Human1","Human2","Chicken1","Human1","Human1","Chicken2","Chicken1","Chicken1","Chicken2","Human1","Chicken2","Human1","Human1","Chicken2","Chicken1","Human1","Human1","Human2","Human1","Chicken2","Human1","Human2","Chicken2","Chicken1","Chicken1","Human1","Human1","Human2","Human1","Chicken2","Human1","Chicken2","Chicken1","Chicken2","Human1","Chicken1","Chicken1","Chicken2","Chicken2","Chicken2","Chicken2","Human1","Chicken1","Chicken2","Human2","Human2","Human2","Chicken2","Chicken2","Human1","Human2","Chicken2","Chicken2","Chicken1","Human1","Chicken2","Chicken2","Human1","Human2","Chicken1","Chicken2","Human1","Human2","Chicken2","Chicken2","Human2","Human1","Human2","Chicken1","Chicken2","Chicken1","Chicken2","Chicken1","Chicken1","Chicken1","Chicken1","Human1","Human2","Chicken2","Chicken1","Chicken1","Human1","Chicken1","Chicken2","Chicken2","Chicken2","Human2","Chicken1","Human1","Human2","Chicken1","Chicken2","Human2","Human2","Human2","Chicken1","Chicken2","Human2","Human1","Human2","Human2","Chicken1","Chicken2","Chicken2","Chicken2","Human1","Human1","Human1","Human2","Human1","Chicken1","Human1","Human2","Chicken1","Chicken1","Human2","Human1","Chicken1","Chicken1","Human1","Human2","Chicken1","Human1","Human2","Chicken2","Chicken2","Chicken1","Chicken2","Human2","Chicken2","Chicken2","Chicken2","Human1","Human2","Chicken1","Human2","Chicken1","Chicken1","Human2","Chicken2","Human1","Human1","Chicken2","Human1","Chicken1","Chicken1","Chicken2","Human1","Human1","Chicken2","Human2","Human1","Human1","Human1","Human1","Chicken1","Chicken1","Chicken2","Chicken2","Human1","Human2","Human1","Human1","Human2","Chicken2","Human2","Chicken2","Human1","Chicken1","Human1","Human2","Human1","Chicken2","Chicken2","Chicken2","Human2","Chicken2","Human2","Chicken2","Chicken2","Chicken1","Chicken2","Human2","Human1","Human2","Chicken2","Human1","Human2","Human2","Chicken1","Human1","Chicken2","Human2","Human1","Chicken2","Human2","Chicken2","Human1","Human2","Chicken2","Chicken2","Human2","Chicken2","Human1","Human2","Human1","Human1","Chicken1","Human2","Human2","Chicken1","Chicken2","Human1","Chicken2","Chicken2","Chicken2","Human2","Human2","Human1","Chicken2","Human2","Chicken1","Human1","Human1","Human2","Human1","Human2","Chicken2","Human1","Human2","Chicken1","Human1","Human2","Human2","Human2","Human2","Chicken2","Human2","Human2","Chicken1","Human2","Human2","Human2","Human1","Human2","Chicken1","Human2","Chicken1","Human2","Human2","Human2","Human2","Human1","Human2","Chicken2","Human2","Human2","Human2","Human1","Human2","Human1","Human2","Human2","Chicken1","Human1","Human1","Human1","Human1","Human1","Human1","Human1","Human1","Human1","Human1","Human1"],[0,2,4,4,5,5,5,5,6,6,6,7,8,8,10,11,14,15,16,17,18,20,22,23,23,25,25,26,29,31,32,37,37,37,39,40,41,41,46,47,49,50,61,67,67,67,68,68,69,69,70,77,82,83,87,89,89,94,100,110,124,128,131,137,139,145,155,161,163,170,174,176,178,196,198,269,360,385,414,435,449,604,624,683,994,1097,2860,3738,3827,4028,4493,5335,5758,5775,6133,6822,6836,8291,8302,8374,8471,8596,8608,8648,8698,8828,8875,8934,8974,9148,9243,9289,9469,9526,9625,9654,9752,9849,9899,9926,9959,9963,10028,10389,10408,10416,10505,10576,10685,10699,10829,11109,11147,11172,11285,11293,11348,11405,11486,11501,11549,11560,11579,11615,11620,11763,11856,12031,12058,12128,12202,12409,12418,12624,12689,12794,12815,12895,13057,13155,13160,13299,13345,13390,13423,13435,13437,13472,13576,13737,13766,13768,13775,13865,13897,13927,13985,14010,14052,14053,14088,14190,14255,14320,14336,14385,14537,14568,14695,14757,14821,14857,14891,14894,14980,15011,15029,15174,15210,15262,15380,15414,15567,15575,15579,15581,15706,15762,15846,15970,16039,16076,16084,16166,16180,16190,16197,16230,16297,16345,16380,16438,16443,16708,16726,16732,16798,16906,16954,17000,17210,17258,17289,17362,17578,17578,17619,17687,17745,17759,17786,17927,18016,18055,18233,18317,18325,18425,18546,18560,18579,18683,18689,18707,18759,18790,18848,18875,18979,19080,19152,19186,19218,19238,19243,19267,19284,19403,19407,19458,19506,19548,19593,19601,19686,19704,19715,19746,19799,19829,19878,19882,19912,19937,20114,20231,20281,20303,20344,20381,20384,20395,20499,20502,20524,20612,20729,20762,20789,20811,20830,20875,20897,20932,20963,20975,21015,21034,21046,21108,21187,21195,21198,21236,21259,21367,21424,21429,21450,21455,21578,21664,21679,21694,21726,21773,21830,21840,21872,21903,22014,22085,22103,22162,22197,22203,22274,22433,22463,22471,22476,22524,22531,22604,22648,22654,22676,22748,22776,22794,22804,22832,22881,22889,22966,22968,23081,23131,23138,23236,23243,23296,23306,23361,23409,23430,23454,23455,23560,23651,23652,23668,23673,23730,23734,23788,23898,23969,23992,24068,24099,24279,24296,24339,24342,24361,24380,24381,24401,24501,24536,24568,24666,24666,24677,24687,24730,24742,24765,24832,24858,24915,24982,24995,25035,25043,25057,25165,25314,25317,25332,25375,25411,25526,25528,25566,25743,25813,25824,25855,25995,26003,26118,26178,26281,26381,26446,26447,26449,26458,26488,26499,26520,26581,26625,26675,26814,26814,26825,26831,26905,26968,26990,27026,27037,27050,27142,27171,27279,27431,27440,27453,27516,27549,27554,27562,27605,27645,27665,27700,27732,27793,27797,27843,27869,27896,27915,27947,28030,28086,28172,28185,28191,28244,28283,28291,28321,28349,28392,28460,28605,28617,28694,28697,28760,28770,28797,28808,28844,28876,28947,28961,29079,29125,29145,29345,29390,29425,29531,29590,29653,29775,29907,29907,29981,30040,30076,30114,30184,30242,30256,30256,30292,30339,30408,30439,30446,30485,30510,30577,30610,30629,30642,30858,30921,30955,31111,31123,31200,31225,31293,31326,31365,31517,31718,31823,31875,31908,31977,31989,32069,32142,32147,32173,32191,32421,32469,32513,32613,32730,32887,32958,32998,33009,33039,33131,33137,33177,33197,33414,33494,33501,33542,33563,33606,33662,33673,33683,33754,33948,33973,34001,34050,34103,34154,34189,34283,34341,34361,34383,34452,34515,34717,34738,34965,35029,35161,35163,35214,35291,35419,35527,35541,36043,36140,36143,36471,36583,36636,36733,36751,36812,36829,37048,37182,37263,37387,37424,37477,37552,37744,37812,37879,37886,37982,38035,38282,38547,38602,38919,38927,39098,39174,39478,39518,39627,39816,39889,40362,40563,40608,41007,41114,41151,41383,41421,42160,42621,42884,43153,43340,43429,43601,43615,43781,43809,44004,44159,45298,45647,45911,46983,47414,47667,48252,50174,50228,50961,51419,51543,51683,52801,52805,53324,54328,54839,55364,56071,56498,57925,62379,64126,64294,66660,72433,76823,88037,108042,124468,150837,184708,191176,203803,233850,235165,254676,698347]],"container":"<table class=\"display\">\n  <thead>\n    <tr>\n      <th> <\/th>\n      <th>SampleID<\/th>\n      <th>Experiment<\/th>\n      <th>Model2<\/th>\n      <th>LibrarySize<\/th>\n    <\/tr>\n  <\/thead>\n<\/table>","options":{"columnDefs":[{"className":"dt-right","targets":4},{"orderable":false,"targets":0}],"order":[],"autoWidth":false,"orderClasses":false}},"evals":[],"jsHooks":[]}</script>
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
## Warning: ggrepel: 10 unlabeled data points (too many overlaps). Consider
## increasing max.overlaps
```

```
## Warning: ggrepel: 54 unlabeled data points (too many overlaps). Consider
## increasing max.overlaps
```

```
## Warning: ggrepel: 19 unlabeled data points (too many overlaps). Consider
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
## [1] 60
```

60 samples lost 



```r
out$df %>% 
  filter(Experiment %!in% c("Mock", "NTC")) %>% 
  filter(LibrarySize < 4000) %>% 
  nrow()
```

```
## [1] 62
```

62 samples lost 



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
## Warning: Removed 41642 row(s) containing missing values (geom_path).
```

<img src="16S_validation_files/figure-html/unnamed-chunk-12-1.png" style="display: block; margin: auto;" />



```r
plot + xlim(c(0, 4000))
```

```
## Warning: Removed 295099 row(s) containing missing values (geom_path).
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
## 89 samples removedbecause they contained fewer reads than `sample.size`.
```

```
## Up to first five removed samples are:
```

```
## CR-40-S166CR-52-S196IR1-69-S198TR1-15-S168TR1-16-S199
```

```
## ...
```

```
## 356OTUs were removed because they are no longer 
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
# physeq_rare %>% 
#   filter_taxa(function(x) sum(x > 0) > 0, TRUE) %>% 
#   phyloseq_alphas(phylo = FALSE) -> alpha_16S

physeq_rare %>%  
  # subset_samples(Model == "Human") %>% 
  plot_richness(measures = c("Observed", "Chao1", "Shannon")) -> alpha_df
```



```r
alpha_df$data %>% 
  filter(Day_of_Treatment > -10) %>% 
  plot_alpha_div_NRP72(facet_formula = paste0("variable ~ Model2"), 
                       path_group = NULL) 
```

```
## Warning: Using alpha for a discrete variable is not advised.
```

<img src="16S_validation_files/figure-html/unnamed-chunk-17-1.png" style="display: block; margin: auto;" />



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

