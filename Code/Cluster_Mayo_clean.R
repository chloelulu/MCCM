############################## update Control (05/26/2025) #######################################
library(readxl)
Day2 <- read_excel("/Users/luyang1/myicloud/Documents/Mayo_project/2023_01_09_Oncobiome/mforge_clean/Data/Day2_all_data_with_Elixhauser_for_MCCM.xlsx")
colnames(Day2)[colnames(Day2)=='Tube.1.Barcode'] ='Tube1Barcode'
colnames(Day2)[colnames(Day2)=='Tube.2.Barcode'] ='Tube2Barcode'
colnames(Day2)[colnames(Day2)=='Tube.3.Barcode'] ='Tube3Barcode'
Day2$Elixhauser_elix.sum <- as.numeric(Day2$Elixhauser_elix.sum)
colnames(Day2)[colnames(Day2)=='Elixhauser_elix.sum'] ='Elix_score'
colnames(Day2)[colnames(Day2)=='sequencing_BIOMEid'] ='BIOME_with_sequencing_data'
Day2 <- as.data.frame(Day2)
Day2$BMI <- as.numeric(Day2$BMI)
rownames(Day2) <- Day2$BIOME_with_sequencing_data

tm <- load(file = '/Users/luyang1//myicloud/Documents/Mayo_project/2023_01_09_Oncobiome/mforge_clean/Data/data.obj.raw.core.RData') 
tm
xx <- data.obj$meta.dat
yy <- data.obj.rff$meta.dat

ctrl <- data.obj$meta.dat[data.obj$meta.dat$Group=='Control',]
comm <- intersect(rownames(ctrl),Day2$BIOME_with_sequencing_data)

Elixhauser_col <- (colnames(Day2)[grep('Elixhauser',colnames(Day2))])[((colnames(Day2)[grep('Elixhauser',colnames(Day2))]) %in% colnames(ctrl))]
Day2 <- Day2[comm,c('BIOME_with_sequencing_data','Elix_score',Elixhauser_col,'Age')]
Day2 <- Day2[comm,]
data.obj$meta.dat[rownames(Day2),colnames(Day2)] <- Day2
identical(data.obj$meta.dat[data.obj$meta.dat$Group=='Cancer',, drop =F], xx[xx$Group=='Cancer',, drop =F])



Day2 <- read_excel("/Users/luyang1/myicloud/Documents/Mayo_project/2023_01_09_Oncobiome/mforge_clean/Data/Day2_all_data_with_Elixhauser_for_MCCM.xlsx")
colnames(Day2)[colnames(Day2)=='Tube.1.Barcode'] ='Tube1Barcode'
colnames(Day2)[colnames(Day2)=='Tube.2.Barcode'] ='Tube2Barcode'
colnames(Day2)[colnames(Day2)=='Tube.3.Barcode'] ='Tube3Barcode'
Day2$Elixhauser_elix.sum <- as.numeric(Day2$Elixhauser_elix.sum)
colnames(Day2)[colnames(Day2)=='Elixhauser_elix.sum'] ='Elix_score'
colnames(Day2)[colnames(Day2)=='sequencing_BIOMEid'] ='BIOME_with_sequencing_data'
Day2 <- as.data.frame(Day2)
Day2$BMI <- as.numeric(Day2$BMI)
rownames(Day2) <- Day2$BIOME_with_sequencing_data

ctrl <- data.obj.rff$meta.dat[data.obj.rff$meta.dat$Group=='Control',]
comm <- intersect(rownames(ctrl),Day2$BIOME_with_sequencing_data)

Elixhauser_col <- (colnames(Day2)[grep('Elixhauser',colnames(Day2))])[((colnames(Day2)[grep('Elixhauser',colnames(Day2))]) %in% colnames(ctrl))]
Day2 <- Day2[comm,c('BIOME_with_sequencing_data','Elix_score',Elixhauser_col,'Age')]
Day2 <- Day2[comm,]
data.obj.rff$meta.dat[rownames(Day2),colnames(Day2)] <- Day2
identical(data.obj.rff$meta.dat[data.obj.rff$meta.dat$Group=='Cancer',, drop =F], yy[yy$Group=='Cancer',, drop =F])

save(data.obj,data.obj.rff,dist.obj,dist.obj.rff, file = '/Users/luyang1/myicloud/Documents/Mayo_project/2023_01_09_Oncobiome/mforge_clean/Data/data.obj.raw.core.RData')


tm <- load("/Users/M216453/Documents/Mayo_project/2023_01_09_Oncobiome/mforge_clean/Data/data.obj.raw.core.RData")
tm
meta <- data.obj$meta.dat
tm <- load("/Users/M216453/Documents/Mayo_project/2023_01_09_Oncobiome/mforge_clean/Data/data.obj.pathway.RData")
tm
identical(rownames(data.obj$meta.dat), rownames(meta))
data.obj$meta.dat <- meta
save(data.obj, dist.obj, file = "/Users/M216453/Documents/Mayo_project/2023_01_09_Oncobiome/mforge_clean/Data/data.obj.pathway.RData")

############################## Add variable (10/18/2024) #######################################
setwd('/Users/luyang1//myicloud/Documents/Mayo_project/2023_01_09_Oncobiome/mforge_clean/')
library(readxl)
library(dplyr)
library(tibble)
Oncobiome_previous_treatment <- read_excel("Data/Oncobiome_previous_treatment.xlsx") 
load(file = 'Data/data.obj.raw.core.RData') 

length(Oncobiome_previous_treatment$BIOME_with_sequencing_data)
sum(Oncobiome_previous_treatment$BIOME_with_sequencing_data %in% data.obj$meta.dat$BIOME_with_sequencing_data)

data.obj$meta.dat[1:4,1:4]
Oncobiome_previous_treatment[1:4,1:2]
meta.dat <- data.obj$meta.dat %>% rownames_to_column(var = "rownames")
meta.dat <- meta.dat %>% left_join(Oncobiome_previous_treatment, by = c("rownames" = "BIOME_with_sequencing_data"))
meta.dat$chemotherapy_last_2_years <- as.factor(meta.dat$chemotherapy_last_2_years)
meta.dat <- meta.dat %>% column_to_rownames(var = "rownames")
data.obj$meta.dat <- meta.dat

meta.dat <- data.obj.rff$meta.dat %>% rownames_to_column(var = "rownames")
meta.dat <- meta.dat %>% left_join(Oncobiome_previous_treatment, by = c("rownames" = "BIOME_with_sequencing_data"))
meta.dat$chemotherapy_last_2_years <- as.factor(meta.dat$chemotherapy_last_2_years)
meta.dat <- meta.dat %>% column_to_rownames(var = "rownames")
data.obj.rff$meta.dat <- meta.dat
# save(data.obj, dist.obj, data.obj.rff, dist.obj.rff, file = 'Result/CancerOnly/data.obj.wk.RData')

## As I am not able to access mforge, run local with Alpha_Beta_DAA_clean.R
variable <- 'chemotherapy_last_2_years'
wd <- '/Users/luyang1//myicloud/Documents/Mayo_project/2023_01_09_Oncobiome/mforge_clean/'
rd <- '/Users/luyang1//myicloud/Documents/Mayo_project/2023_01_09_Oncobiome/mforge_clean/Result/'

covars <- c("Batch","Bristol_score","BMI", "Age", "Sex", "GI_nonGI","Cancer_class","Metastasis","PPI_day_365", "Abx_day_365", 
            "PPI_last_month","Abx_last_month","Charlson_score","Elix_score","Sample_season","Urban" ,"icd10_first_3_name", "Site")
adj.name <- covars[!(covars %in% c(variable,'GI_nonGI',"Charlson_score","Abx_last_month","PPI_last_month","Site","icd10_first_3_name"))]
dir = 'CancerOnly'

############################## Add icd10_first_3_name collapse colon and rectum(07/24/2025) ####################
wd <- '/Users/luyang1/myicloud/Documents/Mayo_project/2023_01_09_Oncobiome/mforge_clean/'
setwd(wd)
source('Code/Stats.R')
try(load_package())
load(file = 'Data/backup/data.obj.raw.core.RData') 

## create icd10_first_3_name_short
data.obj$meta.dat$icd10_first_3_name_short <- data.obj$meta.dat$icd10_first_3_name
data.obj$meta.dat$icd10_first_3_name_short[data.obj$meta.dat$icd10_first_3_name_short %in% c("Malignant neoplasm of colon","Malignant neoplasm of rectum")] <-  "Malignant neoplasm of colon rectum"
data.obj$meta.dat$icd10_first_3_name_short <- tolower(gsub("Malignant neoplasm of |Malignant |malignant |Other and | neoplasm| neoplasm of| neoplasms|neoplasm of |\\,|and ", "", data.obj$meta.dat$icd10_first_3_name_short))
data.obj$meta.dat$icd10_first_3_name_short[grepl('colon rectum',data.obj$meta.dat$icd10_first_3_name_short)] <- 'colorectal'
data.obj$meta.dat$icd10_first_3_name_short[grepl('Control',data.obj$meta.dat$Group)] <- 'healthy'

data.obj.rff$meta.dat$icd10_first_3_name_short <- data.obj.rff$meta.dat$icd10_first_3_name
data.obj.rff$meta.dat$icd10_first_3_name_short[data.obj.rff$meta.dat$icd10_first_3_name_short %in% c("Malignant neoplasm of colon","Malignant neoplasm of rectum")] <-  "Malignant neoplasm of colon rectum"
data.obj.rff$meta.dat$icd10_first_3_name_short <- tolower(gsub("Malignant neoplasm of |Malignant |malignant |Other and | neoplasm| neoplasm of| neoplasms|neoplasm of |\\,|and ", "", data.obj.rff$meta.dat$icd10_first_3_name_short))
data.obj.rff$meta.dat$icd10_first_3_name_short[grepl('colon rectum',data.obj.rff$meta.dat$icd10_first_3_name_short)] <- 'colorectal'
data.obj.rff$meta.dat$icd10_first_3_name_short[grepl('Control',data.obj.rff$meta.dat$Group)] <- 'healthy'

## transform blood related variables into normal distributed
data.obj$meta.dat[c('Erythrocytes','Hematocrit','Neutrophils','MCV','Hemoglobin','Leukocytes','Platelet.Count')] <- as.data.frame(sapply(data.obj$meta.dat[,c('Erythrocytes','Hematocrit','Neutrophils','MCV','Hemoglobin','Leukocytes','Platelet.Count')], as.numeric)) %>% 
  `rownames<-`(rownames(data.obj$meta.dat)) 
data.obj$meta.dat[c('Erythrocytes','Hematocrit','Neutrophils','MCV','Hemoglobin','Leukocytes','Platelet.Count')] <- apply(data.obj$meta.dat[c('Erythrocytes','Hematocrit','Neutrophils','MCV','Hemoglobin','Leukocytes','Platelet.Count')], 2, function(x){
  x[!is.na(x)] <- winzor(x[!is.na(x)], winsor.end = 'both')
  return(x)
})
data.obj$meta.dat <- data.obj$meta.dat %>% mutate(Neutrophils = log2(Neutrophils), Leukocytes = log2(Leukocytes), Platelet.Count = log2(Platelet.Count))

data.obj.rff$meta.dat[c('Erythrocytes','Hematocrit','Neutrophils','MCV','Hemoglobin','Leukocytes','Platelet.Count')] <- as.data.frame(sapply(data.obj.rff$meta.dat[,c('Erythrocytes','Hematocrit','Neutrophils','MCV','Hemoglobin','Leukocytes','Platelet.Count')], as.numeric)) %>% 
  `rownames<-`(rownames(data.obj.rff$meta.dat))  
data.obj.rff$meta.dat[c('Erythrocytes','Hematocrit','Neutrophils','MCV','Hemoglobin','Leukocytes','Platelet.Count')] <- apply(data.obj.rff$meta.dat[c('Erythrocytes','Hematocrit','Neutrophils','MCV','Hemoglobin','Leukocytes','Platelet.Count')], 2, function(x){
  x[!is.na(x)] <- winzor(x[!is.na(x)], winsor.end = 'both')
  return(x)
})
data.obj.rff$meta.dat <- data.obj.rff$meta.dat %>% mutate(Neutrophils = log2(Neutrophils), Leukocytes = log2(Leukocytes), Platelet.Count = log2(Platelet.Count))

## transform blood related variables into binary format
data.obj$meta.dat$Neutrophils_cat <- as.factor(ifelse(is.na(data.obj$meta.dat$Neutrophils_cat),NA, ifelse(data.obj$meta.dat$Neutrophils_cat =='No','No','Yes')))
data.obj$meta.dat$Platelet.Count_cat <- as.factor(ifelse(is.na(data.obj$meta.dat$Platelet.Count_cat),NA, ifelse(data.obj$meta.dat$Platelet.Count_cat =='No','No','Yes')))
data.obj$meta.dat$Hemoglobin_cat <- as.factor(ifelse(is.na(data.obj$meta.dat$Hemoglobin_cat),NA, ifelse(data.obj$meta.dat$Hemoglobin_cat =='No','No','Yes')))

data.obj.rff$meta.dat$Neutrophils_cat <- as.factor(ifelse(is.na(data.obj.rff$meta.dat$Neutrophils_cat),NA, ifelse(data.obj.rff$meta.dat$Neutrophils_cat =='No','No','Yes')))
data.obj.rff$meta.dat$Platelet.Count_cat <- as.factor(ifelse(is.na(data.obj.rff$meta.dat$Platelet.Count_cat),NA, ifelse(data.obj.rff$meta.dat$Platelet.Count_cat =='No','No','Yes')))
data.obj.rff$meta.dat$Hemoglobin_cat <- as.factor(ifelse(is.na(data.obj.rff$meta.dat$Hemoglobin_cat),NA, ifelse(data.obj.rff$meta.dat$Hemoglobin_cat =='No','No','Yes')))

## relevel GI_nonGI
data.obj$meta.dat$GI_nonGI <- factor(data.obj$meta.dat$GI_nonGI,levels=c('non_GI','GI'))
data.obj.rff$meta.dat$GI_nonGI <- factor(data.obj.rff$meta.dat$GI_nonGI,levels=c('non_GI','GI'))

## add 03/02/2025, change Abx_last_month PPI_last_month NA into No
data.obj$meta.dat$Abx_last_month <- as.character(data.obj$meta.dat$Abx_last_month)
data.obj$meta.dat$Abx_last_month[is.na(data.obj$meta.dat$Abx_last_month)] <- 'No'
data.obj$meta.dat$Abx_last_month <- as.factor(data.obj$meta.dat$Abx_last_month)
data.obj.rff$meta.dat$Abx_last_month <- as.character(data.obj.rff$meta.dat$Abx_last_month)
data.obj.rff$meta.dat$Abx_last_month[is.na(data.obj.rff$meta.dat$Abx_last_month)] <- 'No'
data.obj.rff$meta.dat$Abx_last_month <- as.factor(data.obj.rff$meta.dat$Abx_last_month)

data.obj$meta.dat$PPI_last_month <- as.character(data.obj$meta.dat$PPI_last_month)
data.obj$meta.dat$PPI_last_month[is.na(data.obj$meta.dat$PPI_last_month)] <- 'No'
data.obj$meta.dat$PPI_last_month <- as.factor(data.obj$meta.dat$PPI_last_month)
data.obj.rff$meta.dat$PPI_last_month <- as.character(data.obj.rff$meta.dat$PPI_last_month)
data.obj.rff$meta.dat$PPI_last_month[is.na(data.obj.rff$meta.dat$PPI_last_month)] <- 'No'
data.obj.rff$meta.dat$PPI_last_month <- as.factor(data.obj.rff$meta.dat$PPI_last_month)

## create ealy onset
data.obj$meta.dat$early_onset <- as.factor(ifelse(data.obj$meta.dat$Age<=50,'Yes','No') )
data.obj.rff$meta.dat$early_onset <- as.factor(ifelse(data.obj.rff$meta.dat$Age<=50,'Yes','No') )

data.obj.rff$meta.dat$icd10_first_3_name_short <- factor(data.obj.rff$meta.dat$icd10_first_3_name_short)
data.obj.rff$meta.dat$icd10_first_3_name_short <- relevel(data.obj.rff$meta.dat$icd10_first_3_name_short, ref = "healthy")

data.obj$meta.dat$icd10_first_3_name_short <- factor(data.obj$meta.dat$icd10_first_3_name_short)
data.obj$meta.dat$icd10_first_3_name_short <- relevel(data.obj$meta.dat$icd10_first_3_name_short, ref = "healthy")

save(data.obj, dist.obj, data.obj.rff, dist.obj.rff, file = 'Data/data.obj.raw.core.RData')

###### functional data
load(file = 'Data/backup/data.obj.pathway.RData') 
## create icd10_first_3_name_short
data.obj$meta.dat$icd10_first_3_name_short <- data.obj$meta.dat$icd10_first_3_name
data.obj$meta.dat$icd10_first_3_name_short[data.obj$meta.dat$icd10_first_3_name_short %in% c("Malignant neoplasm of colon","Malignant neoplasm of rectum")] <-  "Malignant neoplasm of colon rectum"
data.obj$meta.dat$icd10_first_3_name_short <- tolower(gsub("Malignant neoplasm of |Malignant |malignant |Other and | neoplasm| neoplasm of| neoplasms|neoplasm of |\\,|and ", "", data.obj$meta.dat$icd10_first_3_name_short))
data.obj$meta.dat$icd10_first_3_name_short[grepl('colon rectum',data.obj$meta.dat$icd10_first_3_name_short)] <- 'colorectal'
data.obj$meta.dat$icd10_first_3_name_short[grepl('Control',data.obj$meta.dat$Group)] <- 'healthy'

## transform blood related variables into normal distributed
data.obj$meta.dat[c('Erythrocytes','Hematocrit','Neutrophils','MCV','Hemoglobin','Leukocytes','Platelet.Count')] <- as.data.frame(sapply(data.obj$meta.dat[,c('Erythrocytes','Hematocrit','Neutrophils','MCV','Hemoglobin','Leukocytes','Platelet.Count')], as.numeric)) %>% 
  `rownames<-`(rownames(data.obj$meta.dat)) 
data.obj$meta.dat[c('Erythrocytes','Hematocrit','Neutrophils','MCV','Hemoglobin','Leukocytes','Platelet.Count')] <- apply(data.obj$meta.dat[c('Erythrocytes','Hematocrit','Neutrophils','MCV','Hemoglobin','Leukocytes','Platelet.Count')], 2, function(x){
  x[!is.na(x)] <- winzor(x[!is.na(x)], winsor.end = 'both')
  return(x)
})
data.obj$meta.dat <- data.obj$meta.dat %>% mutate(Neutrophils = log2(Neutrophils), Leukocytes = log2(Leukocytes), Platelet.Count = log2(Platelet.Count))

## transform blood related variables into binary format
data.obj$meta.dat$Neutrophils_cat <- as.factor(ifelse(is.na(data.obj$meta.dat$Neutrophils_cat),NA, ifelse(data.obj$meta.dat$Neutrophils_cat =='No','No','Yes')))
data.obj$meta.dat$Platelet.Count_cat <- as.factor(ifelse(is.na(data.obj$meta.dat$Platelet.Count_cat),NA, ifelse(data.obj$meta.dat$Platelet.Count_cat =='No','No','Yes')))
data.obj$meta.dat$Hemoglobin_cat <- as.factor(ifelse(is.na(data.obj$meta.dat$Hemoglobin_cat),NA, ifelse(data.obj$meta.dat$Hemoglobin_cat =='No','No','Yes')))

## relevel GI_nonGI
data.obj$meta.dat$GI_nonGI <- factor(data.obj$meta.dat$GI_nonGI,levels=c('non_GI','GI'))

## add 03/02/2025, change Abx_last_month PPI_last_month NA into No
data.obj$meta.dat$Abx_last_month <- as.character(data.obj$meta.dat$Abx_last_month)
data.obj$meta.dat$Abx_last_month[is.na(data.obj$meta.dat$Abx_last_month)] <- 'No'
data.obj$meta.dat$Abx_last_month <- as.factor(data.obj$meta.dat$Abx_last_month)
data.obj.rff$meta.dat$Abx_last_month <- as.character(data.obj.rff$meta.dat$Abx_last_month)
data.obj.rff$meta.dat$Abx_last_month[is.na(data.obj.rff$meta.dat$Abx_last_month)] <- 'No'
data.obj.rff$meta.dat$Abx_last_month <- as.factor(data.obj.rff$meta.dat$Abx_last_month)

## create ealy onset
data.obj$meta.dat$early_onset <- as.factor(ifelse(data.obj$meta.dat$Age<=50,'Yes','No') )

data.obj$meta.dat$icd10_first_3_name_short <- factor(data.obj$meta.dat$icd10_first_3_name_short)
data.obj$meta.dat$icd10_first_3_name_short <- relevel(data.obj$meta.dat$icd10_first_3_name_short, ref = "healthy")

save(data.obj, dist.obj, file = 'Data/data.obj.pathway.RData')






############################## Add Smoking,vital_status,OS_months,prior_chemotherapy (02/18/2026) #######################################
setwd('/Users/luyang1//myicloud/Documents/Mayo_project/2023_01_09_Oncobiome/mforge_clean/')
load(file = 'Data/data.obj.raw.core.RData') # Now the old data is in Data/backup_02182026/ folder
Supplementary_Table_1 <- read_excel("Code/Submission/Supplementary tables/Supplementary Table 1.xlsx", sheet = "TableS1-2")
Supplementary_Table_1$SampleID[!(Supplementary_Table_1$SampleID %in% rownames(data.obj$meta.dat))]
rownames(data.obj$meta.dat)[!(rownames(data.obj$meta.dat) %in% Supplementary_Table_1$SampleID)]
data.obj$meta.dat <- data.obj$meta.dat %>% rownames_to_column('SampleID') %>% inner_join(Supplementary_Table_1[,c("SampleID","vital_status","OS_months","smoking_category","prior_chemotherapy")])
data.obj.rff$meta.dat <- data.obj.rff$meta.dat %>% rownames_to_column('SampleID') %>% inner_join(Supplementary_Table_1[,c("SampleID","vital_status","OS_months","smoking_category","prior_chemotherapy")])
data.obj$meta.dat <- data.obj$meta.dat %>% column_to_rownames('SampleID')
data.obj.rff$meta.dat <- data.obj.rff$meta.dat %>% column_to_rownames('SampleID')
save(data.obj,dist.obj, data.obj.rff, dist.obj.rff, file = 'Data/data.obj.raw.core.RData')


load(file = 'Data/data.obj.pathway.RData') # Now the old data is in Data/backup_02182026/ folder
Supplementary_Table_1 <- read_excel("Code/Submission/Supplementary tables/Supplementary Table 1.xlsx", sheet = "TableS1-2")
Supplementary_Table_1$SampleID[!(Supplementary_Table_1$SampleID %in% rownames(data.obj$meta.dat))]
rownames(data.obj$meta.dat)[!(rownames(data.obj$meta.dat) %in% Supplementary_Table_1$SampleID)]
data.obj$meta.dat <- data.obj$meta.dat %>% rownames_to_column('SampleID') %>% inner_join(Supplementary_Table_1[,c("SampleID","vital_status","OS_months","smoking_category","prior_chemotherapy")])
data.obj$meta.dat <- data.obj$meta.dat %>% column_to_rownames('SampleID')
save(data.obj,dist.obj, file = 'Data/data.obj.pathway.RData')


############################## Add Ruben's treatment and NAs (04/14/2026) #######################################
setwd('/Users/luyang1//myicloud/Documents/Mayo_project/2023_01_09_Oncobiome/mforge_clean/')
load(file = 'Data/backup_04142026/data.obj.raw.core.RData') 
Supplementary_Table_1 <- read_excel("Code/Submission/Supplementary tables/old/Supplementary Table 1_RM.xlsx", sheet = "TableS1-2")
Supplementary_Table_1$SampleID[!(Supplementary_Table_1$SampleID %in% rownames(data.obj$meta.dat))]
rownames(data.obj$meta.dat)[!(rownames(data.obj$meta.dat) %in% Supplementary_Table_1$SampleID)]
apply(Supplementary_Table_1, 2, function(x) sum(x=='NA'))
which(apply(Supplementary_Table_1, 2, function(x) sum(is.na(x)))>0)
apply(data.obj$meta.dat, 2, function(x) sum(is.na(x)))
apply(data.obj$meta.dat, 2, function(x) sum(x=='NA'))
which(apply(data.obj$meta.dat, 2, function(x) sum(is.na(x)))>0)
TableS16 <- read_excel("Code/Submission/Supplementary tables/old/Supplementary Table 1_RM.xlsx", sheet = "TableS1-6")
TableS16$sample[!(TableS16$sample %in% Supplementary_Table_1$SequencingName)]
TableS16 <- TableS16[,c(1:2)]
Supplementary_Table_1 <- inner_join(TableS16,Supplementary_Table_1) 
data.obj$meta.dat <- data.obj$meta.dat %>% rownames_to_column('SampleID') %>% 
  inner_join(Supplementary_Table_1[,c("SampleID","sample","planned_Tx_class","planned_Tx_category","planned_Tx","planned_Radiation")])
data.obj.rff$meta.dat <- data.obj.rff$meta.dat %>% rownames_to_column('SampleID') %>% 
  inner_join(Supplementary_Table_1[,c("SampleID","sample","planned_Tx_class","planned_Tx_category","planned_Tx","planned_Radiation")])
data.obj$meta.dat[data.obj$meta.dat == "NA"] <- NA
data.obj.rff$meta.dat[data.obj.rff$meta.dat == "NA"] <- NA
data.obj$meta.dat <- data.obj$meta.dat %>% column_to_rownames('SampleID')
data.obj.rff$meta.dat <- data.obj.rff$meta.dat %>% column_to_rownames('SampleID')
save(data.obj,dist.obj, data.obj.rff, dist.obj.rff, file = 'Data/data.obj.raw.core.RData')


load(file = 'Data/backup_04142026/data.obj.pathway.RData') # Now the old data is in Data/backup_02182026/ folder
Supplementary_Table_1 <- read_excel("Code/Submission/Supplementary tables/old/Supplementary Table 1_RM.xlsx", sheet = "TableS1-2")
Supplementary_Table_1$SampleID[!(Supplementary_Table_1$SampleID %in% rownames(data.obj$meta.dat))]
rownames(data.obj$meta.dat)[!(rownames(data.obj$meta.dat) %in% Supplementary_Table_1$SampleID)]
TableS16 <- read_excel("Code/Submission/Supplementary tables/old/Supplementary Table 1_RM.xlsx", sheet = "TableS1-6")
TableS16$sample[!(TableS16$sample %in% Supplementary_Table_1$SequencingName)]
TableS16 <- TableS16[,c(1:2)]
Supplementary_Table_1 <- inner_join(TableS16,Supplementary_Table_1) 
data.obj$meta.dat <- data.obj$meta.dat %>% rownames_to_column('SampleID') %>% 
  inner_join(Supplementary_Table_1[,c("SampleID","sample","planned_Tx_class","planned_Tx_category","planned_Tx","planned_Radiation")])
data.obj$meta.dat[data.obj$meta.dat == "NA"] <- NA
data.obj$meta.dat <- data.obj$meta.dat %>% column_to_rownames('SampleID')
save(data.obj,dist.obj, file = 'Data/data.obj.pathway.RData')


load(file = 'Data/backup_04142026/data.obj.mph.RData') # Now the old data is in Data/backup_02182026/ folder
Supplementary_Table_1 <- read_excel("Code/Submission/Supplementary tables/old/Supplementary Table 1_RM.xlsx", sheet = "TableS1-2")
Supplementary_Table_1$SampleID[!(Supplementary_Table_1$SampleID %in% rownames(data.obj$meta.dat))]
rownames(data.obj$meta.dat)[!(rownames(data.obj$meta.dat) %in% Supplementary_Table_1$SampleID)]
TableS16 <- read_excel("Code/Submission/Supplementary tables/old/Supplementary Table 1_RM.xlsx", sheet = "TableS1-6")
TableS16$sample[!(TableS16$sample %in% Supplementary_Table_1$SequencingName)]
TableS16 <- TableS16[,c(1:2)]
Supplementary_Table_1 <- inner_join(TableS16,Supplementary_Table_1) 
data.obj$meta.dat <- data.obj$meta.dat %>% rownames_to_column('SampleID') %>% 
  inner_join(Supplementary_Table_1[,c("SampleID","sample","planned_Tx_class","planned_Tx_category","planned_Tx","planned_Radiation")])
data.obj$meta.dat[data.obj$meta.dat == "NA"] <- NA
data.obj$meta.dat <- data.obj$meta.dat %>% column_to_rownames('SampleID')
save(data.obj,dist.obj, file = 'Data/data.obj.mph.RData')



############################# Define script ######################################
file_dir <- dirname(rstudioapi::getSourceEditorContext()$path)
file_dir <- dirname(dirname(dirname(dirname(file_dir))))
rd <- paste0(file_dir,'/Result/')
taxon_script <- paste0(file_dir,'/Code/Submission/MayoOncobiomeStudy/Code/Alpha_Beta_DAA_clean.R')
func_script <- paste0(file_dir,'/Code/Submission/MayoOncobiomeStudy/Code/Alpha_Beta_DAA_func_clean.R')

############################## Cancer Only #######################################
file_dir <- dirname(rstudioapi::getSourceEditorContext()$path)
file_dir <- dirname(dirname(dirname(dirname(file_dir))))

wd <- file_dir
rd <- paste0(file_dir,'/Result/')

setwd(wd)
source(paste0(file_dir,"/Code/Submission/MayoOncobiomeStudy/Code/Stats.R"))
try(load_package())

getwd()
load(file = 'Data/data.obj.raw.core.RData') 


setwd(rd)
dir <- 'CancerOnly'
if(!dir.exists(dir)){dir.create(dir)}
setwd(dir)
getwd()

ind <- data.obj$meta.dat$Group == 'Cancer'
data.obj <- subset_data(data.obj, ind)
dist.obj <- subset_dist(dist.obj, ind)

ind <- data.obj.rff$meta.dat$Group == 'Cancer'
data.obj.rff <- subset_data(data.obj.rff, ind)
dist.obj.rff <- subset_dist(dist.obj.rff, ind)
save(data.obj, dist.obj, data.obj.rff, dist.obj.rff, file = 'data.obj.wk.RData')


covars <- c("Batch","Bristol_score","BMI", "Age", "Sex", "GI_nonGI","Cancer_class","Metastasis","PPI_day_365", "Abx_day_365", 
            "PPI_last_month","Abx_last_month","Charlson_score","Elix_score","Sample_season","Urban" ,"icd10_first_3_name_short", "Site",
            "smoking_category","prior_chemotherapy")
blood.names <- c('Erythrocytes','Hematocrit','Neutrophils','MCV','Hemoglobin','Leukocytes','Platelet.Count',
                 'Neutrophils_cat','Platelet.Count_cat','Hemoglobin_cat',
                 "neut_neutropenia2","bone_marrow_suppression2","Hb_anemia2","Pl_thrombocytopenia2","neut_neutropenia_c","Hb_anemia_c","Pl_thrombocytopenia_c")
elix.names <- c("Elixhauser_CHF","Elixhauser_Arrhythmia","Elixhauser_Valvular","Elixhauser_PHTN","Elixhauser_PVD","Elixhauser_HTN","Elixhauser_Paralysis",  
                "Elixhauser_NeuroOther","Elixhauser_Pulmonary","Elixhauser_DM","Elixhauser_DMcx","Elixhauser_Hypothyroid","Elixhauser_Renal","Elixhauser_Liver",      
                "Elixhauser_PUD" ,"Elixhauser_HIV","Elixhauser_Lymphoma","Elixhauser_Mets","Elixhauser_Tumor","Elixhauser_Rheumatic","Elixhauser_Coagulopathy",
                "Elixhauser_Obesity","Elixhauser_WeightLoss","Elixhauser_FluidsLytes","Elixhauser_BloodLoss","Elixhauser_Anemia","Elixhauser_Alcohol","Elixhauser_Drugs",
                "Elixhauser_Psychoses","Elixhauser_Depression")
variables <- c(covars[!(covars %in% c("Batch"))], blood.names, elix.names)
dir <- 'CancerOnly'
for(variable in variables){
  with(list(commandArgs = function(...) c("--args", variable, dir)), source(taxon_script))
}

# setwd("/Users/luyang1/myicloud/Documents/Mayo_project/2023_01_09_Oncobiome/mforge_clean/Result/CancerOnly/DAA/")
# dirs <- list.dirs(recursive = F,full.names = F)
# xx <- NULL
# for(dir in dirs){
#   diff.obj <- NULL
#   load(paste0(dir,'/',dir,'_ZicoSeq.Rdata'))
#   x <- colnames(diff.obj$coef.list$Species)[grep(dir,colnames(diff.obj$coef.list$Species))]
#   if(sum(grepl('1$',x))>0) xx <- c(xx, dir)
# }
# 
# variables <- xx
############################## Cancer Only (smoking for lung cancer) #######################################

file_dir <- dirname(rstudioapi::getSourceEditorContext()$path)
file_dir <- dirname(dirname(dirname(dirname(file_dir))))

wd <- file_dir
rd <- paste0(file_dir,'/Result/')

setwd(wd)
source(paste0(file_dir,"/Code/Submission/MayoOncobiomeStudy/Code/Stats.R"))
try(load_package())

getwd()
load(file = 'Data/data.obj.raw.core.RData') 


setwd(rd)
dir <- 'LungCancer_Smoking'
if(!dir.exists(dir)){dir.create(dir)}
setwd(dir)
getwd()

ind <- data.obj$meta.dat$icd10_first_3_name_short == 'bronchus lung'
data.obj <- subset_data(data.obj, ind)
dist.obj <- subset_dist(dist.obj, ind)

ind <- data.obj.rff$meta.dat$icd10_first_3_name_short == 'bronchus lung'
data.obj.rff <- subset_data(data.obj.rff, ind)
dist.obj.rff <- subset_dist(dist.obj.rff, ind)

data.obj$meta.dat$smoking_category2 <- data.obj$meta.dat$smoking_category
data.obj$meta.dat$smoking_category2[data.obj$meta.dat$smoking_category2 %in% c("Current","Former")] <- "Yes"

data.obj.rff$meta.dat$smoking_category2 <- data.obj.rff$meta.dat$smoking_category
data.obj.rff$meta.dat$smoking_category2[data.obj.rff$meta.dat$smoking_category2 %in% c("Current","Former")] <- "Yes"

save(data.obj, dist.obj, data.obj.rff, dist.obj.rff, file = 'data.obj.wk.RData')

variables <- c("smoking_category","smoking_category2")
dir <- 'LungCancer_Smoking'
for(variable in variables){
  with(list(commandArgs = function(...) c("--args", variable, dir)), source(taxon_script))
}



############################## Control vs sub Cancer X #######################################
file_dir <- dirname(rstudioapi::getSourceEditorContext()$path)
file_dir <- dirname(dirname(dirname(dirname(file_dir))))

wd <- file_dir
rd <- paste0(file_dir,'/Result/subCancerX_Control/')

if(!dir.exists(rd)) dir.create(rd)
setwd(wd)
source(paste0(file_dir,"/Code/Submission/MayoOncobiomeStudy/Code/Stats.R"))
try(load_package())

load(file = 'Data/data.obj.raw.core.RData') 

data.obj2 <- data.obj
data.obj.rff2 <- data.obj.rff
dist.obj2 <- dist.obj
dist.obj.rff2 <- dist.obj.rff

cancer.type <- names(which(table(data.obj2$meta.dat$icd10_first_3_name_short)>15))
cancer.type <- cancer.type[cancer.type != 'healthy']

for(i in 1:length(cancer.type)){
  setwd(rd)
  cancer <- cancer.type[i]
  cat('[',cancer,']\n')
  dir <- paste0('Control-',cancer.type[i])
  if(!dir.exists(dir)){dir.create(dir)}
  setwd(dir)
  getwd()
  
  ind <- data.obj2$meta.dat$icd10_first_3_name_short %in% c('healthy',cancer)
  data.obj <- subset_data(data.obj2, ind)
  dist.obj <- subset_dist(dist.obj2, ind)
  cat(sum(ind),'\n')
  
  ind <- data.obj.rff2$meta.dat$icd10_first_3_name_short %in% c('healthy',cancer)
  data.obj.rff <- subset_data(data.obj.rff2, ind)
  dist.obj.rff <- subset_dist(dist.obj.rff2, ind)
  cat(sum(ind),'\n')
  
  data.obj$meta.dat$icd10_first_3_name_short <- factor(data.obj$meta.dat$icd10_first_3_name_short,levels=c('healthy',cancer))
  data.obj.rff$meta.dat$icd10_first_3_name_short <- factor(data.obj.rff$meta.dat$icd10_first_3_name_short,levels=c('healthy',cancer))
  
  save(data.obj, dist.obj, data.obj.rff, dist.obj.rff, file = 'data.obj.wk.RData')
  setwd('..')
}

# setwd("/Users/luyang1/myicloud/Documents/Mayo_project/2023_01_09_Oncobiome/mforge_clean/Result/subCancerX_Control/")
# dirs <- list.dirs(recursive = F,full.names = F)
# xx <- NULL
# for(dir in dirs){
#   diff.obj <- NULL
#   load(paste0(dir,'/DAA/icd10_first_3_name_short/icd10_first_3_name_short_ZicoSeq.Rdata'))
#   x <- colnames(diff.obj$coef.list$Species)[grep('icd10_first_3_name_short',colnames(diff.obj$coef.list$Species))]
#   if(grepl('1$',x)) xx <- c(xx, dir)
# }
# cancer.type <- cancer.type[cancer.type %in% gsub('Control-','',xx)]

variable <- c("icd10_first_3_name_short")
for(dir in paste0('subCancerX_Control/Control-',cancer.type)){
  with(list(commandArgs = function(...) c("--args", variable, dir)), source(taxon_script))
}


############################## sub Cancer X vs [Cancers-sub Cancer X] #######################################
file_dir <- dirname(rstudioapi::getSourceEditorContext()$path)
file_dir <- dirname(dirname(dirname(dirname(file_dir))))

wd <- file_dir
rd <- paste0(file_dir,'/Result/subCancerX-Ex/')

if(!dir.exists(rd)) dir.create(rd)
setwd(wd)
source(paste0(file_dir,"/Code/Submission/MayoOncobiomeStudy/Code/Stats.R"))
try(load_package())

load(file = 'Data/data.obj.raw.core.RData') 

data.obj2 <- data.obj
data.obj.rff2 <- data.obj.rff
dist.obj2 <- dist.obj
dist.obj.rff2 <- dist.obj.rff

cancer.type <- names(which(sort(table(data.obj$meta.dat$icd10_first_3_name_short))>15)) 
cancer.type <- cancer.type[cancer.type!='healthy']

for(i in 1:length(cancer.type)){
  setwd(rd)
  cancer <- cancer.type[i]
  cat('[',cancer,']\n')
  dir <- cancer.type[i]
  if(!dir.exists(dir)){dir.create(dir)}
  setwd(dir)
  getwd()
  
  data.obj <- data.obj2
  data.obj.rff <- data.obj.rff2
  dist.obj <- dist.obj2
  dist.obj.rff <- dist.obj.rff2

  data.obj$meta.dat$icd10_first_3_name_short <- as.character(data.obj$meta.dat$icd10_first_3_name_short)
  samIDs <- rownames(data.obj$meta.dat[data.obj$meta.dat$icd10_first_3_name_short %in% cancer.type,])
  data.obj <- subset_data(data.obj, samIDs = samIDs)
  dist.obj <- subset_dist(dist.obj, samIDs = samIDs)
  
  data.obj.rff$meta.dat$icd10_first_3_name_short <- as.character(data.obj.rff$meta.dat$icd10_first_3_name_short)
  samIDs <- rownames(data.obj.rff$meta.dat[data.obj.rff$meta.dat$icd10_first_3_name_short %in% cancer.type,])
  data.obj.rff <- subset_data(data.obj.rff, samIDs = samIDs)
  dist.obj.rff <- subset_dist(dist.obj.rff, samIDs = samIDs)
  
  data.obj$meta.dat$icd10_first_3_name_short <- as.factor(data.obj$meta.dat$icd10_first_3_name_short)
  data.obj.rff$meta.dat$icd10_first_3_name_short <- as.factor(data.obj.rff$meta.dat$icd10_first_3_name_short)
  
  idx <- !(data.obj$meta.dat$icd10_first_3_name_short %in% c(cancer))
  data.obj$meta.dat$icd10_first_3_name_short <- as.character(data.obj$meta.dat$icd10_first_3_name_short)
  data.obj$meta.dat$icd10_first_3_name_short[idx] <- 'Others'
  dist.obj <- subset_dist(dist.obj, samIDs = rownames(data.obj$meta.dat))
  
  idx <- !(data.obj.rff$meta.dat$icd10_first_3_name_short %in% c(cancer))
  data.obj.rff$meta.dat$icd10_first_3_name_short <- as.character(data.obj.rff$meta.dat$icd10_first_3_name_short)
  data.obj.rff$meta.dat$icd10_first_3_name_short[idx] <- 'Others'
  dist.obj.rff <- subset_dist(dist.obj.rff, samIDs = rownames(data.obj.rff$meta.dat))
  cat(sum(idx),'\n')
  
  data.obj$meta.dat$icd10_first_3_name_short <- factor(data.obj$meta.dat$icd10_first_3_name_short,levels=c('Others',cancer))
  data.obj.rff$meta.dat$icd10_first_3_name_short <- factor(data.obj.rff$meta.dat$icd10_first_3_name_short,levels=c('Others',cancer))
  table(data.obj.rff$meta.dat$icd10_first_3_name_short)
  
  save(data.obj, dist.obj, data.obj.rff, dist.obj.rff, file = 'data.obj.wk.RData')
  gc()
}

variable <- c("icd10_first_3_name_short")
for(dir in paste0('subCancerX-Ex/',cancer.type)){
  with(list(commandArgs = function(...) c("--args", variable, dir)), source(taxon_script))
}


# setwd("/Users/luyang1/myicloud/Documents/Mayo_project/2023_01_09_Oncobiome/mforge_clean/Result/subCancerX-Ex/")
# dirs <- list.dirs(recursive = F,full.names = F)
# xx <- NULL
# for(dir in dirs){
#   diff.obj <- NULL
#   load(paste0(dir,'/DAA/icd10_first_3_name_short/icd10_first_3_name_short_ZicoSeq.Rdata'))
#   x <- colnames(diff.obj$coef.list$Species)[grep('icd10_first_3_name_short',colnames(diff.obj$coef.list$Species))]
#   if(grepl('1$',x)) xx <- c(xx, dir)
# }
# xx
############################## Pan-Cancer #######################################
file_dir <- dirname(rstudioapi::getSourceEditorContext()$path)
file_dir <- dirname(dirname(dirname(dirname(file_dir))))

wd <- file_dir
rd <- paste0(file_dir,'/Result/')

if(!dir.exists(rd)) dir.create(rd)
setwd(wd)
source(paste0(file_dir,"/Code/Submission/MayoOncobiomeStudy/Code/Stats.R"))
try(load_package())

load(file = 'Data/data.obj.raw.core.RData') 

setwd(rd)
dir <- 'PanCancer'
if(!dir.exists(dir)){dir.create(dir)}
setwd(dir)
getwd()

data.obj$meta.dat$Cancer_class <- as.character(data.obj$meta.dat$Cancer_class)
data.obj$meta.dat$Cancer_class[data.obj$meta.dat$Group=='Control'] <- 'Control'

data.obj.rff$meta.dat$Cancer_class <- as.character(data.obj.rff$meta.dat$Cancer_class)
data.obj.rff$meta.dat$Cancer_class[data.obj.rff$meta.dat$Group=='Control'] <- 'Control'

data.obj$meta.dat$Group <- factor(data.obj$meta.dat$Group,levels=c('Control','Cancer'))
data.obj.rff$meta.dat$Group <- factor(data.obj.rff$meta.dat$Group,levels=c('Control','Cancer'))

save(data.obj, dist.obj, data.obj.rff, dist.obj.rff, file = 'data.obj.wk.RData')

variable <- "Group"; dir <- 'PanCancer'
with(list(commandArgs = function(...) c("--args", variable, dir)), source(taxon_script))


############################## early onset #######################################
file_dir <- dirname(rstudioapi::getSourceEditorContext()$path)
file_dir <- dirname(dirname(dirname(dirname(file_dir))))

wd <- file_dir
rd <- paste0(file_dir,'/Result/')

if(!dir.exists(rd)) dir.create(rd)
setwd(wd)
source(paste0(file_dir,"/Code/Submission/MayoOncobiomeStudy/Code/Stats.R"))
try(load_package())

load(file = 'Data/data.obj.raw.core.RData') 

ind <- data.obj$meta.dat$Group == 'Cancer'
data.obj2 <- subset_data(data.obj, ind)
dist.obj2 <- subset_dist(dist.obj, ind)

ind <- data.obj.rff$meta.dat$Group == 'Cancer'
data.obj.rff2 <- subset_data(data.obj.rff, ind)
dist.obj.rff2 <- subset_dist(dist.obj.rff, ind)

setwd(rd)
dir <- 'EarlyOnset'
if(!dir.exists(dir)){dir.create(dir)}
setwd(dir)
getwd()

tb <- cbind(table(data.obj.rff2$meta.dat$icd10_first_3_name_short, data.obj.rff2$meta.dat$early_onset))
idx <- names(which(rowSums(tb >15) ==2))
for(j in idx){
  setwd(paste0(rd,'EarlyOnset'))
  if(!dir.exists(j)){dir.create(j)}
  setwd(j)
  getwd()
  
  data.obj <- data.obj2
  data.obj.rff <- data.obj.rff2
  dist.obj <- dist.obj2
  dist.obj.rff <- dist.obj.rff2
  
  ind2 <- data.obj.rff$meta.dat$icd10_first_3_name_short == j
  data.obj.rff <- subset_data(data.obj.rff, ind2)
  dist.obj.rff <- subset_dist(dist.obj.rff, ind2)
  
  ind2 <- data.obj$meta.dat$icd10_first_3_name_short == j
  data.obj <- subset_data(data.obj, ind2)
  dist.obj <- subset_dist(dist.obj, ind2)
  table(data.obj$meta.dat$early_onset)
  
  save(data.obj, dist.obj, data.obj.rff, dist.obj.rff, file = 'data.obj.wk.RData')

}


variable <- 'early_onset'
for(dir in paste0('EarlyOnset/',idx)){
  with(list(commandArgs = function(...) c("--args", variable, dir)), source(taxon_script))
}




### ===== ECI subset analysis: pancancer ======
file_dir <- dirname(rstudioapi::getSourceEditorContext()$path)
file_dir <- dirname(dirname(dirname(dirname(file_dir))))

wd <- file_dir
rd <- paste0(file_dir,'/Result/')

if(!dir.exists(rd)) dir.create(rd)
setwd(wd)
source(paste0(file_dir,"/Code/Submission/MayoOncobiomeStudy/Code/Stats.R"))
try(load_package())

load(file = 'Data/data.obj.raw.core.RData') 

setwd(rd)
dir <- 'ECI_subset'
if(!dir.exists(dir)){dir.create(dir)}
setwd(dir)
getwd()

data.obj$meta.dat$ECI2 <- case_when(
  data.obj$meta.dat$Elix_score <=1  ~ "LowCancer",
  data.obj$meta.dat$Elix_score %in% c(2,3) ~ "MediumCancer",
  data.obj$meta.dat$Elix_score >= 4 ~ "HighCancer"
)
data.obj$meta.dat$ECI2[data.obj$meta.dat$Group=='Control'& data.obj$meta.dat$Elix_score<4] <- 'Control'
data.obj$meta.dat$ECI2[data.obj$meta.dat$Group=='Control' & data.obj$meta.dat$Elix_score>=4] <- NA

data.obj.rff$meta.dat$ECI2 <- case_when(
  data.obj.rff$meta.dat$Elix_score <=1 ~ "LowCancer",
  data.obj.rff$meta.dat$Elix_score %in% c(2,3) ~ "MediumCancer",
  data.obj.rff$meta.dat$Elix_score >= 4 ~ "HighCancer"
)
data.obj.rff$meta.dat$ECI2[data.obj.rff$meta.dat$Group=='Control' & data.obj.rff$meta.dat$Elix_score<4] <- 'Control'
data.obj.rff$meta.dat$ECI2[data.obj.rff$meta.dat$Group=='Control' & data.obj.rff$meta.dat$Elix_score>=4] <- NA

data.obj2 <- data.obj
data.obj.rff2 <- data.obj.rff
dist.obj2 <- dist.obj
dist.obj.rff2 <- dist.obj.rff

for(i in c('LowCancer','HighCancer')){
  setwd(rd)
  dir <- 'ECI_subset'
  setwd(dir)
  dir <- paste0('Control-',i)
  if(!dir.exists(dir)){dir.create(dir)}
  setwd(dir)
  getwd()
  
  ind <- data.obj2$meta.dat$ECI2 %in% c('Control',i)
  data.obj <- subset_data(data.obj2, ind)
  dist.obj <- subset_dist(dist.obj2, ind)
  cat(sum(ind),'\n')
  
  ind <- data.obj.rff2$meta.dat$ECI2 %in% c('Control',i)
  data.obj.rff <- subset_data(data.obj.rff2, ind)
  dist.obj.rff <- subset_dist(dist.obj.rff2, ind)
  cat(sum(ind),'\n')
  
  data.obj$meta.dat$ECI2 <- factor(data.obj$meta.dat$ECI2,levels=c('Control',i))
  data.obj.rff$meta.dat$ECI2 <- factor(data.obj.rff$meta.dat$ECI2,levels=c('Control',i))
  
  save(data.obj, dist.obj, data.obj.rff, dist.obj.rff, file = 'data.obj.wk.RData')
}


variable <- "ECI2"; dir <- 'ECI_subset'
for(dir in paste0('ECI_subset/Control-',c('LowCancer','HighCancer'))){
  with(list(commandArgs = function(...) c("--args", variable, dir)), source(taxon_script))
}


### ===== ECI extreme case subset analysis: pancancer (Cancer-lowECI:01; Cancer-highECI:>=4; HV: ECI<=4; including cancers with n>15)  ======
file_dir <- dirname(rstudioapi::getSourceEditorContext()$path)
file_dir <- dirname(dirname(dirname(dirname(file_dir))))

wd <- file_dir
rd <- paste0(file_dir,'/Result/')

if(!dir.exists(rd)) dir.create(rd)
setwd(wd)
source(paste0(file_dir,"/Code/Submission/MayoOncobiomeStudy/Code/Stats.R"))
try(load_package())

load(file = 'Data/data.obj.raw.core.RData') 

setwd(rd)
dir <- 'ECI_subset4'
if(!dir.exists(dir)){dir.create(dir)}
setwd(dir)
getwd()

data.obj$meta.dat$ECI22 <- case_when(
  data.obj$meta.dat$Elix_score <= 1 ~ "LowCancer",
  data.obj$meta.dat$Elix_score %in% c(1,2,3) ~ "MediumCancer",
  data.obj$meta.dat$Elix_score >= 4 ~ "HighCancer"
)
data.obj$meta.dat$ECI22[data.obj$meta.dat$Group=='Control'& data.obj$meta.dat$Elix_score<4] <- 'Control'
data.obj$meta.dat$ECI22[data.obj$meta.dat$Group=='Control' & data.obj$meta.dat$Elix_score>=4] <- NA

data.obj.rff$meta.dat$ECI22 <- case_when(
  data.obj.rff$meta.dat$Elix_score <= 1 ~ "LowCancer",
  data.obj.rff$meta.dat$Elix_score %in% c(1,2,3) ~ "MediumCancer",
  data.obj.rff$meta.dat$Elix_score >= 4 ~ "HighCancer"
)
data.obj.rff$meta.dat$ECI22[data.obj.rff$meta.dat$Group=='Control' & data.obj.rff$meta.dat$Elix_score<4] <- 'Control'
data.obj.rff$meta.dat$ECI22[data.obj.rff$meta.dat$Group=='Control' & data.obj.rff$meta.dat$Elix_score>=4] <- NA

## subset cancer n>15
cancer.type <- names(which(table(data.obj$meta.dat$icd10_first_3_name_short)>15))

data.obj$meta.dat$icd10_first_3_name_short <- as.character(data.obj$meta.dat$icd10_first_3_name_short)
samIDs <- rownames(data.obj$meta.dat[data.obj$meta.dat$icd10_first_3_name_short %in% cancer.type,])
data.obj <- subset_data(data.obj, samIDs = samIDs)
dist.obj <- subset_dist(dist.obj, samIDs = samIDs)

data.obj.rff$meta.dat$icd10_first_3_name_short <- as.character(data.obj.rff$meta.dat$icd10_first_3_name_short)
samIDs <- rownames(data.obj.rff$meta.dat[data.obj.rff$meta.dat$icd10_first_3_name_short %in% cancer.type,])
data.obj.rff <- subset_data(data.obj.rff, samIDs = samIDs)
dist.obj.rff <- subset_dist(dist.obj.rff, samIDs = samIDs)

data.obj$meta.dat$icd10_first_3_name_short <- as.factor(data.obj$meta.dat$icd10_first_3_name_short)
data.obj.rff$meta.dat$icd10_first_3_name_short <- as.factor(data.obj.rff$meta.dat$icd10_first_3_name_short)


data.obj2 <- data.obj
data.obj.rff2 <- data.obj.rff
dist.obj2 <- dist.obj
dist.obj.rff2 <- dist.obj.rff

for(i in c('LowCancer','HighCancer')){
  setwd(rd)
  dir <- 'ECI_subset4'
  setwd(dir)
  dir <- paste0('Control-',i)
  if(!dir.exists(dir)){dir.create(dir)}
  setwd(dir)
  getwd()
  
  ind <- data.obj2$meta.dat$ECI22 %in% c('Control',i)
  data.obj <- subset_data(data.obj2, ind)
  dist.obj <- subset_dist(dist.obj2, ind)
  cat(sum(ind),'\n')
  
  ind <- data.obj.rff2$meta.dat$ECI22 %in% c('Control',i)
  data.obj.rff <- subset_data(data.obj.rff2, ind)
  dist.obj.rff <- subset_dist(dist.obj.rff2, ind)
  cat(sum(ind),'\n')
  
  data.obj$meta.dat$ECI22 <- factor(data.obj$meta.dat$ECI22,levels=c('Control',i))
  data.obj.rff$meta.dat$ECI22 <- factor(data.obj.rff$meta.dat$ECI22,levels=c('Control',i))
  
  save(data.obj, dist.obj, data.obj.rff, dist.obj.rff, file = 'data.obj.wk.RData')
}


variable <- "ECI22"; dir <- 'ECI_subset4'
for(dir in paste0('ECI_subset4/Control-',c('LowCancer','HighCancer'))){
  with(list(commandArgs = function(...) c("--args", variable, dir)), source(taxon_script))
}







### ===== ECI subset analysis: cancer X vs HV (lowCancer==0&1, also exclude healthy ECI>=4) ======
file_dir <- dirname(rstudioapi::getSourceEditorContext()$path)
file_dir <- dirname(dirname(dirname(dirname(file_dir))))

wd <- file_dir
rd <- paste0(file_dir,'/Result/')

if(!dir.exists(rd)) dir.create(rd)
setwd(wd)
source(paste0(file_dir,"/Code/Submission/MayoOncobiomeStudy/Code/Stats.R"))
try(load_package())

load(file = 'Data/data.obj.raw.core.RData') 

setwd(rd)
dir <- 'ECI_subset'
if(!dir.exists(dir)){dir.create(dir)}
setwd(dir)
getwd()

data.obj$meta.dat$ECI2 <- case_when(
  data.obj$meta.dat$Elix_score <=1  ~ "LowCancer",
  data.obj$meta.dat$Elix_score %in% c(2,3) ~ "MediumCancer",
  data.obj$meta.dat$Elix_score >= 4 ~ "HighCancer"
)
data.obj$meta.dat$ECI2[data.obj$meta.dat$Group=='Control'& data.obj$meta.dat$Elix_score<4] <- 'Control'
data.obj$meta.dat$ECI2[data.obj$meta.dat$Group=='Control' & data.obj$meta.dat$Elix_score>=4] <- NA

data.obj.rff$meta.dat$ECI2 <- case_when(
  data.obj.rff$meta.dat$Elix_score <=1 ~ "LowCancer",
  data.obj.rff$meta.dat$Elix_score %in% c(2,3) ~ "MediumCancer",
  data.obj.rff$meta.dat$Elix_score >= 4 ~ "HighCancer"
)
data.obj.rff$meta.dat$ECI2[data.obj.rff$meta.dat$Group=='Control' & data.obj.rff$meta.dat$Elix_score<4] <- 'Control'
data.obj.rff$meta.dat$ECI2[data.obj.rff$meta.dat$Group=='Control' & data.obj.rff$meta.dat$Elix_score>=4] <- NA

cancers <- names(which(sort(cbind(table(data.obj$meta.dat$icd10_first_3_name_short,data.obj$meta.dat$ECI2))[,"LowCancer"]) > 2))
cancers2 <- names(which(table(data.obj$meta.dat$icd10_first_3_name_short)>15))
cancers2 <- cancers2[cancers2!="healthy"]
cancers <- intersect(cancers, cancers2)

data.obj2 <- data.obj
data.obj.rff2 <- data.obj.rff
dist.obj2 <- dist.obj
dist.obj.rff2 <- dist.obj.rff


## Subset cancer 
for(cancer in cancers){
  setwd(rd)
  dir <- 'ECI_subset'
  setwd(dir)
  cat('[',cancer,']\n')
  dir <- cancer
  if(!dir.exists(dir)){dir.create(dir)}
  setwd(dir)
  getwd()
  
  ind <- data.obj2$meta.dat$icd10_first_3_name_short %in% c('healthy',cancer)
  data.obj3 <- subset_data(data.obj2, ind)
  dist.obj3 <- subset_dist(dist.obj2, ind)
  data.obj3$meta.dat$icd10_first_3_name_short <- droplevels(data.obj3$meta.dat$icd10_first_3_name_short)
  cat(sum(ind),'\n')
  
  ind <- data.obj.rff2$meta.dat$icd10_first_3_name_short %in% c('healthy',cancer)
  data.obj.rff3 <- subset_data(data.obj.rff2, ind)
  dist.obj.rff3 <- subset_dist(dist.obj.rff2, ind)
  data.obj.rff3$meta.dat$icd10_first_3_name_short <- droplevels(data.obj.rff3$meta.dat$icd10_first_3_name_short)
  
  cat(sum(ind),'\n')
  
  for(i in c('LowCancer','MediumCancer','HighCancer')){
    setwd(rd)
    dir <- 'ECI_subset'
    setwd(dir)
    cat('[',cancer,']\n')
    dir <- cancer
    if(!dir.exists(dir)){dir.create(dir)}
    setwd(dir)
    
    dir <- paste0('Control-',i)
    if(!dir.exists(dir)){dir.create(dir)}
    setwd(dir)
    getwd()
    
    ind <- data.obj3$meta.dat$ECI2 %in% c('Control',i)
    data.obj <- subset_data(data.obj3, ind)
    dist.obj <- subset_dist(dist.obj3, ind)
    cat(sum(ind),'\n')
    
    ind <- data.obj.rff3$meta.dat$ECI2 %in% c('Control',i)
    data.obj.rff <- subset_data(data.obj.rff3, ind)
    dist.obj.rff <- subset_dist(dist.obj.rff3, ind)
    cat(sum(ind),'\n')
    
    data.obj$meta.dat$ECI2 <- factor(data.obj$meta.dat$ECI2,levels=c('Control',i))
    data.obj.rff$meta.dat$ECI2 <- factor(data.obj.rff$meta.dat$ECI2,levels=c('Control',i))
    
    save(data.obj, dist.obj, data.obj.rff, dist.obj.rff, file = 'data.obj.wk.RData')
}
}

variable <- "ECI2"; dir <- 'ECI_subset'
for(cancer in cancers){
  for(dir in paste0('ECI_subset/',cancer,'/Control-',c('LowCancer'))){
    with(list(commandArgs = function(...) c("--args", variable, dir)), source(taxon_script))
    }
}




### ===== ECI extreme case subset analysis: cancer X vs HV (lowCancer==0, also exclude healthy ECI>=4)======
file_dir <- dirname(rstudioapi::getSourceEditorContext()$path)
file_dir <- dirname(dirname(dirname(dirname(file_dir))))

wd <- file_dir
rd <- paste0(file_dir,'/Result/')

if(!dir.exists(rd)) dir.create(rd)
setwd(wd)
source(paste0(file_dir,"/Code/Submission/MayoOncobiomeStudy/Code/Stats.R"))
try(load_package())

load(file = 'Data/data.obj.raw.core.RData') 

setwd(rd)
dir <- 'ECI_subset0'
if(!dir.exists(dir)){dir.create(dir)}
setwd(dir)
getwd()

data.obj$meta.dat$ECI22 <- case_when(
  data.obj$meta.dat$Elix_score ==0  ~ "LowCancer",
  data.obj$meta.dat$Elix_score %in% c(1,2,3) ~ "MediumCancer",
  data.obj$meta.dat$Elix_score >= 4 ~ "HighCancer"
)
data.obj$meta.dat$ECI22[data.obj$meta.dat$Group=='Control'& data.obj$meta.dat$Elix_score<4] <- 'Control'
data.obj$meta.dat$ECI22[data.obj$meta.dat$Group=='Control' & data.obj$meta.dat$Elix_score>=4] <- NA

data.obj.rff$meta.dat$ECI22 <- case_when(
  data.obj.rff$meta.dat$Elix_score ==0 ~ "LowCancer",
  data.obj.rff$meta.dat$Elix_score %in% c(2,3) ~ "MediumCancer",
  data.obj.rff$meta.dat$Elix_score >= 4 ~ "HighCancer"
)
data.obj.rff$meta.dat$ECI22[data.obj.rff$meta.dat$Group=='Control' & data.obj.rff$meta.dat$Elix_score<4] <- 'Control'
data.obj.rff$meta.dat$ECI22[data.obj.rff$meta.dat$Group=='Control' & data.obj.rff$meta.dat$Elix_score>=4] <- NA

cancers <- names(which(sort(cbind(table(data.obj$meta.dat$icd10_first_3_name_short,data.obj$meta.dat$ECI22))[,"LowCancer"]) > 2))
cancers2 <- names(which(table(data.obj$meta.dat$icd10_first_3_name_short)>15))
cancers2 <- cancers2[cancers2!="healthy"]
cancers <- intersect(cancers, cancers2)

data.obj2 <- data.obj
data.obj.rff2 <- data.obj.rff
dist.obj2 <- dist.obj
dist.obj.rff2 <- dist.obj.rff


## Subset cancer 
for(cancer in cancers){
  setwd(rd)
  dir <- 'ECI_subset0'
  setwd(dir)
  cat('[',cancer,']\n')
  dir <- cancer
  if(!dir.exists(dir)){dir.create(dir)}
  setwd(dir)
  getwd()
  
  ind <- data.obj2$meta.dat$icd10_first_3_name_short %in% c('healthy',cancer)
  data.obj3 <- subset_data(data.obj2, ind)
  dist.obj3 <- subset_dist(dist.obj2, ind)
  data.obj3$meta.dat$icd10_first_3_name_short <- droplevels(data.obj3$meta.dat$icd10_first_3_name_short)
  cat(sum(ind),'\n')
  
  ind <- data.obj.rff2$meta.dat$icd10_first_3_name_short %in% c('healthy',cancer)
  data.obj.rff3 <- subset_data(data.obj.rff2, ind)
  dist.obj.rff3 <- subset_dist(dist.obj.rff2, ind)
  data.obj.rff3$meta.dat$icd10_first_3_name_short <- droplevels(data.obj.rff3$meta.dat$icd10_first_3_name_short)
  
  cat(sum(ind),'\n')
  
  for(i in c('LowCancer','HighCancer')){
    setwd(rd)
    dir <- 'ECI_subset0'
    setwd(dir)
    cat('[',cancer,']\n')
    dir <- cancer
    if(!dir.exists(dir)){dir.create(dir)}
    setwd(dir)
    
    dir <- paste0('Control-',i)
    if(!dir.exists(dir)){dir.create(dir)}
    setwd(dir)
    getwd()
    
    ind <- data.obj3$meta.dat$ECI22 %in% c('Control',i)
    data.obj <- subset_data(data.obj3, ind)
    dist.obj <- subset_dist(dist.obj3, ind)
    cat(sum(ind),'\n')
    
    ind <- data.obj.rff3$meta.dat$ECI22 %in% c('Control',i)
    data.obj.rff <- subset_data(data.obj.rff3, ind)
    dist.obj.rff <- subset_dist(dist.obj.rff3, ind)
    cat(sum(ind),'\n')
    
    data.obj$meta.dat$ECI22 <- factor(data.obj$meta.dat$ECI22,levels=c('Control',i))
    data.obj.rff$meta.dat$ECI22 <- factor(data.obj.rff$meta.dat$ECI22,levels=c('Control',i))
    
    save(data.obj, dist.obj, data.obj.rff, dist.obj.rff, file = 'data.obj.wk.RData')
  }
}

variable <- "ECI22"; dir <- 'ECI_subset0'
for(cancer in cancers){
  for(dir in paste0('ECI_subset0/',cancer,'/Control-',c('LowCancer','HighCancer'))){
    with(list(commandArgs = function(...) c("--args", variable, dir)), source(taxon_script))
  }
}




### ===== Abx 365 subset analysis (pancancer: Abx Y vs HV, Abx N vs HV)======
file_dir <- dirname(rstudioapi::getSourceEditorContext()$path)
file_dir <- dirname(dirname(dirname(dirname(file_dir))))

wd <- file_dir
rd <- paste0(file_dir,'/Result/')

if(!dir.exists(rd)) dir.create(rd)
setwd(wd)
source(paste0(file_dir,"/Code/Submission/MayoOncobiomeStudy/Code/Stats.R"))
try(load_package())

load(file = 'Data/data.obj.raw.core.RData') 

setwd(rd)
dir <- 'Abx_day_365_subset'
if(!dir.exists(dir)){dir.create(dir)}
setwd(dir)
getwd()

data.obj2 <- data.obj
data.obj.rff2 <- data.obj.rff
dist.obj2 <- dist.obj
dist.obj.rff2 <- dist.obj.rff

for(i in c('Y','N')){
  setwd(rd)
  dir <- 'Abx_day_365_subset'
  setwd(dir)
  dir <- paste0('Control-',i)
  if(!dir.exists(dir)){dir.create(dir)}
  setwd(dir)
  getwd()
  
  data.obj2$meta.dat$Abx_day_365 <- as.character(data.obj2$meta.dat$Abx_day_365)
  data.obj2$meta.dat$Abx_day_365[data.obj2$meta.dat$Group=='Control'] <- 'Control'
  ind <- data.obj2$meta.dat$Abx_day_365 %in% c('Control',i)
  data.obj <- subset_data(data.obj2, ind)
  dist.obj <- subset_dist(dist.obj2, ind)
  cat(sum(ind),'\n')
  
  data.obj.rff2$meta.dat$Abx_day_365 <- as.character(data.obj.rff2$meta.dat$Abx_day_365)
  data.obj.rff2$meta.dat$Abx_day_365[data.obj.rff2$meta.dat$Group=='Control'] <- 'Control'
  ind <- data.obj.rff2$meta.dat$Abx_day_365 %in% c('Control',i)
  data.obj.rff <- subset_data(data.obj.rff2, ind)
  dist.obj.rff <- subset_dist(dist.obj.rff2, ind)
  cat(sum(ind),'\n')
  
  data.obj$meta.dat$Abx_day_365 <- factor(data.obj$meta.dat$Abx_day_365,levels=c('Control',i))
  data.obj.rff$meta.dat$Abx_day_365 <- factor(data.obj.rff$meta.dat$Abx_day_365,levels=c('Control',i))
  
  save(data.obj, dist.obj, data.obj.rff, dist.obj.rff, file = 'data.obj.wk.RData')
}


variable <- "Abx_day_365"; dir <- 'Abx_day_365_subset'
for(dir in paste0('Abx_day_365_subset/Control-',c('Y','N'))){
  with(list(commandArgs = function(...) c("--args", variable, dir)), source(taxon_script))
}



### ===== Abx 365 subset cancer vs HV analysis ======
file_dir <- dirname(rstudioapi::getSourceEditorContext()$path)
file_dir <- dirname(dirname(dirname(dirname(file_dir))))

wd <- file_dir
rd <- paste0(file_dir,'/Result/')

if(!dir.exists(rd)) dir.create(rd)
setwd(wd)
source(paste0(file_dir,"/Code/Submission/MayoOncobiomeStudy/Code/Stats.R"))
try(load_package())

load(file = 'Data/data.obj.raw.core.RData') 

## keep cancers with Abx Y/N > 5 & cancer X sample > 15
cancers <- names(which(rowSums(cbind(table(data.obj$meta.dat$icd10_first_3_name_short,data.obj$meta.dat$Abx_day_365)) > 2)==2))
cancers2 <- names(which(table(data.obj$meta.dat$icd10_first_3_name_short)>15))
cancers2 <- cancers2[cancers2!="healthy"]
cancers <- intersect(cancers, cancers2)

setwd(rd)
dir <- 'Abx_day_365_subset'
if(!dir.exists(dir)){dir.create(dir)}
setwd(dir)
getwd()

data.obj2 <- data.obj
data.obj.rff2 <- data.obj.rff
dist.obj2 <- dist.obj
dist.obj.rff2 <- dist.obj.rff

## Subset cancer 
for(cancer in cancers){
  setwd(rd)
  dir <- 'Abx_day_365_subset'
  setwd(dir)
  cat('[',cancer,']\n')
  dir <- cancer
  if(!dir.exists(dir)){dir.create(dir)}
  setwd(dir)
  getwd()
  
  ind <- data.obj2$meta.dat$icd10_first_3_name_short %in% c('healthy',cancer)
  data.obj3 <- subset_data(data.obj2, ind)
  dist.obj3 <- subset_dist(dist.obj2, ind)
  cat(sum(ind),'\n')
  
  ind <- data.obj.rff2$meta.dat$icd10_first_3_name_short %in% c('healthy',cancer)
  data.obj.rff3 <- subset_data(data.obj.rff2, ind)
  dist.obj.rff3 <- subset_dist(dist.obj.rff2, ind)
  cat(sum(ind),'\n')

  for(i in c('Y','N')){
    setwd(rd)
    dir <- 'Abx_day_365_subset'
    setwd(dir)
    cat('[',cancer,']\n')
    dir <- cancer
    if(!dir.exists(dir)){dir.create(dir)}
    setwd(dir)
    getwd()
    
    dir <- paste0('Control-',i)
    if(!dir.exists(dir)){dir.create(dir)}
    setwd(dir)
    getwd()
    
    data.obj3$meta.dat$Abx_day_365 <- as.character(data.obj3$meta.dat$Abx_day_365)
    data.obj3$meta.dat$Abx_day_365[data.obj3$meta.dat$Group=='Control'] <- 'Control'
    ind <- data.obj3$meta.dat$Abx_day_365 %in% c('Control',i)
    data.obj <- subset_data(data.obj3, ind)
    dist.obj <- subset_dist(dist.obj3, ind)
    cat(sum(ind),'\n')
    
    data.obj.rff3$meta.dat$Abx_day_365 <- as.character(data.obj.rff3$meta.dat$Abx_day_365)
    data.obj.rff3$meta.dat$Abx_day_365[data.obj.rff3$meta.dat$Group=='Control'] <- 'Control'
    ind <- data.obj.rff3$meta.dat$Abx_day_365 %in% c('Control',i)
    data.obj.rff <- subset_data(data.obj.rff3, ind)
    dist.obj.rff <- subset_dist(dist.obj.rff3, ind)
    cat(sum(ind),'\n')
    
    data.obj$meta.dat$Abx_day_365 <- factor(data.obj$meta.dat$Abx_day_365,levels=c('Control',i))
    data.obj.rff$meta.dat$Abx_day_365 <- factor(data.obj.rff$meta.dat$Abx_day_365,levels=c('Control',i))
    
    save(data.obj, dist.obj, data.obj.rff, dist.obj.rff, file = 'data.obj.wk.RData')
  }
}



variable <- "Abx_day_365"
for(cancer in cancers[5:17]){
  for(dir in paste0('Abx_day_365_subset/',cancer,'/Control-',c('Y','N'))){
    with(list(commandArgs = function(...) c("--args", variable, dir)), source(taxon_script))
  }
}




### ===== PPI 365 subset analysis (pancancer: PPI Y vs HV, PPI N vs HV) ======
file_dir <- dirname(rstudioapi::getSourceEditorContext()$path)
file_dir <- dirname(dirname(dirname(dirname(file_dir))))

wd <- file_dir
rd <- paste0(file_dir,'/Result/')

if(!dir.exists(rd)) dir.create(rd)
setwd(wd)
source(paste0(file_dir,"/Code/Submission/MayoOncobiomeStudy/Code/Stats.R"))
try(load_package())

load(file = 'Data/data.obj.raw.core.RData') 

setwd(rd)
dir <- 'PPI_day_365_subset'
if(!dir.exists(dir)){dir.create(dir)}
setwd(dir)
getwd()

data.obj2 <- data.obj
data.obj.rff2 <- data.obj.rff
dist.obj2 <- dist.obj
dist.obj.rff2 <- dist.obj.rff

for(i in c('Y','N')){
  setwd(rd)
  dir <- 'PPI_day_365_subset'
  setwd(dir)
  dir <- paste0('Control-',i)
  if(!dir.exists(dir)){dir.create(dir)}
  setwd(dir)
  getwd()
  
  data.obj2$meta.dat$PPI_day_365 <- as.character(data.obj2$meta.dat$PPI_day_365)
  data.obj2$meta.dat$PPI_day_365[data.obj2$meta.dat$Group=='Control'] <- 'Control'
  ind <- data.obj2$meta.dat$PPI_day_365 %in% c('Control',i)
  data.obj <- subset_data(data.obj2, ind)
  dist.obj <- subset_dist(dist.obj2, ind)
  cat(sum(ind),'\n')
  
  data.obj.rff2$meta.dat$PPI_day_365 <- as.character(data.obj.rff2$meta.dat$PPI_day_365)
  data.obj.rff2$meta.dat$PPI_day_365[data.obj.rff2$meta.dat$Group=='Control'] <- 'Control'
  ind <- data.obj.rff2$meta.dat$PPI_day_365 %in% c('Control',i)
  data.obj.rff <- subset_data(data.obj.rff2, ind)
  dist.obj.rff <- subset_dist(dist.obj.rff2, ind)
  cat(sum(ind),'\n')
  
  data.obj$meta.dat$PPI_day_365 <- factor(data.obj$meta.dat$PPI_day_365,levels=c('Control',i))
  data.obj.rff$meta.dat$PPI_day_365 <- factor(data.obj.rff$meta.dat$PPI_day_365,levels=c('Control',i))
  
  save(data.obj, dist.obj, data.obj.rff, dist.obj.rff, file = 'data.obj.wk.RData')
}


variable <- "PPI_day_365"; dir <- 'PPI_day_365_subset'
for(dir in paste0('PPI_day_365_subset/Control-',c('Y','N'))){
  with(list(commandArgs = function(...) c("--args", variable, dir)), source(taxon_script))
}




### ===== PPI 365 subset cancer vs HV analysis ========
file_dir <- dirname(rstudioapi::getSourceEditorContext()$path)
file_dir <- dirname(dirname(dirname(dirname(file_dir))))

wd <- file_dir
rd <- paste0(file_dir,'/Result/')

if(!dir.exists(rd)) dir.create(rd)
setwd(wd)
source(paste0(file_dir,"/Code/Submission/MayoOncobiomeStudy/Code/Stats.R"))
try(load_package())

load(file = 'Data/data.obj.raw.core.RData') 

## keep cancers with PPI Y/N > 5
cancers <- names(which(rowSums(cbind(table(data.obj$meta.dat$icd10_first_3_name_short,data.obj$meta.dat$PPI_day_365)) > 5)==2))
cancers <- names(which(rowSums(cbind(table(data.obj$meta.dat$icd10_first_3_name_short,data.obj$meta.dat$Abx_day_365))>5)==2))
cancers2 <- names(which(table(data.obj$meta.dat$icd10_first_3_name_short)>15))
cancers2 <- cancers2[cancers2!="healthy"]
cancers <- intersect(cancers, cancers2)


setwd(rd)
dir <- 'PPI_day_365_subset'
if(!dir.exists(dir)){dir.create(dir)}
setwd(dir)
getwd()

data.obj2 <- data.obj
data.obj.rff2 <- data.obj.rff
dist.obj2 <- dist.obj
dist.obj.rff2 <- dist.obj.rff

## Subset cancer 
for(cancer in cancers){
  setwd(rd)
  dir <- 'PPI_day_365_subset'
  setwd(dir)
  cat('[',cancer,']\n')
  dir <- cancer
  if(!dir.exists(dir)){dir.create(dir)}
  setwd(dir)
  getwd()
  
  ind <- data.obj2$meta.dat$icd10_first_3_name_short %in% c('healthy',cancer)
  data.obj3 <- subset_data(data.obj2, ind)
  dist.obj3 <- subset_dist(dist.obj2, ind)
  cat(sum(ind),'\n')
  
  ind <- data.obj.rff2$meta.dat$icd10_first_3_name_short %in% c('healthy',cancer)
  data.obj.rff3 <- subset_data(data.obj.rff2, ind)
  dist.obj.rff3 <- subset_dist(dist.obj.rff2, ind)
  cat(sum(ind),'\n')
  
  for(i in c('Y','N')){
    setwd(rd)
    dir <- 'PPI_day_365_subset'
    setwd(dir)
    cat('[',cancer,']\n')
    dir <- cancer
    if(!dir.exists(dir)){dir.create(dir)}
    setwd(dir)
    getwd()
    
    dir <- paste0('Control-',i)
    if(!dir.exists(dir)){dir.create(dir)}
    setwd(dir)
    getwd()
    
    data.obj3$meta.dat$PPI_day_365 <- as.character(data.obj3$meta.dat$PPI_day_365)
    data.obj3$meta.dat$PPI_day_365[data.obj3$meta.dat$Group=='Control'] <- 'Control'
    ind <- data.obj3$meta.dat$PPI_day_365 %in% c('Control',i)
    data.obj <- subset_data(data.obj3, ind)
    dist.obj <- subset_dist(dist.obj3, ind)
    cat(sum(ind),'\n')
    
    data.obj.rff3$meta.dat$PPI_day_365 <- as.character(data.obj.rff3$meta.dat$PPI_day_365)
    data.obj.rff3$meta.dat$PPI_day_365[data.obj.rff3$meta.dat$Group=='Control'] <- 'Control'
    ind <- data.obj.rff3$meta.dat$PPI_day_365 %in% c('Control',i)
    data.obj.rff <- subset_data(data.obj.rff3, ind)
    dist.obj.rff <- subset_dist(dist.obj.rff3, ind)
    cat(sum(ind),'\n')
    
    data.obj$meta.dat$PPI_day_365 <- factor(data.obj$meta.dat$PPI_day_365,levels=c('Control',i))
    data.obj.rff$meta.dat$PPI_day_365 <- factor(data.obj.rff$meta.dat$PPI_day_365,levels=c('Control',i))
    
    save(data.obj, dist.obj, data.obj.rff, dist.obj.rff, file = 'data.obj.wk.RData')
  }
}



variable <- "PPI_day_365"
for(cancer in cancers){
  for(dir in paste0('PPI_day_365_subset/',cancer,'/Control-',c('Y','N'))){
    with(list(commandArgs = function(...) c("--args", variable, dir)), source(taxon_script))
  }
}

############################## add  Cancer Only (ECI binary, <=median or ECI=0,1 as low) #######################################
file_dir <- dirname(rstudioapi::getSourceEditorContext()$path)
file_dir <- dirname(dirname(dirname(dirname(file_dir))))

wd <- file_dir
rd <- paste0(file_dir,'/Result/')

setwd(wd)
source(paste0(file_dir,"/Code/Submission/MayoOncobiomeStudy/Code/Stats.R"))
try(load_package())

getwd()
load(file = 'Data/data.obj.raw.core.RData') 


setwd(rd)
dir <- 'CancerOnly'
if(!dir.exists(dir)){dir.create(dir)}
setwd(dir)
getwd()

ind <- data.obj$meta.dat$Group == 'Cancer'
data.obj <- subset_data(data.obj, ind)
dist.obj <- subset_dist(dist.obj, ind)

ind <- data.obj.rff$meta.dat$Group == 'Cancer'
data.obj.rff <- subset_data(data.obj.rff, ind)
dist.obj.rff <- subset_dist(dist.obj.rff, ind)

data.obj$meta.dat$Elix_score2_med <- ifelse(data.obj$meta.dat$Elix_score <= median(data.obj$meta.dat$Elix_score), 'Low','High')
data.obj$meta.dat$Elix_score2_01 <- ifelse(data.obj$meta.dat$Elix_score <=1, 'Low', 'High')
data.obj$meta.dat$Elix_score2_0 <- ifelse(data.obj$meta.dat$Elix_score <=1, 'Low', 
                                          ifelse(data.obj$meta.dat$Elix_score >=4,'High',NA))
data.obj$meta.dat$Elix_score3 <- ifelse(data.obj$meta.dat$Elix_score <=1, 'Low', 
                                          ifelse(data.obj$meta.dat$Elix_score >=4,'High',"Mid"))

data.obj$meta.dat$Elix_score2_med <- factor(data.obj$meta.dat$Elix_score2_med,levels=c('Low',"High"))
data.obj$meta.dat$Elix_score2_01 <- factor(data.obj$meta.dat$Elix_score2_01,levels=c('Low',"High"))
data.obj$meta.dat$Elix_score2_0 <- factor(data.obj$meta.dat$Elix_score2_0,levels=c('Low',"High"))
data.obj$meta.dat$Elix_score3 <- factor(data.obj$meta.dat$Elix_score3,levels=c('Low',"Mid","High"))


data.obj.rff$meta.dat$Elix_score2_med <- ifelse(data.obj.rff$meta.dat$Elix_score <= median(data.obj.rff$meta.dat$Elix_score), 'Low','High')
data.obj.rff$meta.dat$Elix_score2_01 <- ifelse(data.obj.rff$meta.dat$Elix_score <=1, 'Low', 'High')
data.obj.rff$meta.dat$Elix_score2_0 <- ifelse(data.obj.rff$meta.dat$Elix_score <=1, 'Low', 
                                              ifelse(data.obj.rff$meta.dat$Elix_score >=4,'High',NA))
data.obj.rff$meta.dat$Elix_score3 <- ifelse(data.obj.rff$meta.dat$Elix_score <=1, 'Low', 
                                        ifelse(data.obj.rff$meta.dat$Elix_score >=4,'High',"Mid"))


data.obj.rff$meta.dat$Elix_score2_med <- factor(data.obj.rff$meta.dat$Elix_score2_med,levels=c('Low',"High"))
data.obj.rff$meta.dat$Elix_score2_01 <- factor(data.obj.rff$meta.dat$Elix_score2_01,levels=c('Low',"High"))
data.obj.rff$meta.dat$Elix_score2_0 <- factor(data.obj.rff$meta.dat$Elix_score2_0,levels=c('Low',"High"))
data.obj.rff$meta.dat$Elix_score3 <- factor(data.obj.rff$meta.dat$Elix_score3,levels=c('Low',"Mid","High"))


save(data.obj, dist.obj, data.obj.rff, dist.obj.rff, file = 'data.obj.wk.RData')

variables <- c('Elix_score2_med','Elix_score2_01','Elix_score2_0','Elix_score3')
dir <- 'CancerOnly'
for(variable in variables){
  with(list(commandArgs = function(...) c("--args", variable, dir)), source(taxon_script))
}



##---------------------------- functional data analysis-----------------------------------
############################## Cancer Only func #######################################
file_dir <- dirname(rstudioapi::getSourceEditorContext()$path)
file_dir <- dirname(dirname(dirname(dirname(file_dir))))

wd <- file_dir
rd <- paste0(file_dir,'/Result/')

setwd(wd)
source(paste0(file_dir,"/Code/Submission/MayoOncobiomeStudy/Code/Stats.R"))
try(load_package())

setwd(wd)
load(file = paste0('Data/data.obj.pathway.RData')) 

func.type <- 'pathway'
setwd(rd)
dir <- 'CancerOnly_func'
if(!dir.exists(dir)){dir.create(dir)}
setwd(dir)
getwd()

if(!dir.exists(func.type)){dir.create(func.type)}
setwd(func.type)

ind <- data.obj$meta.dat$Group == 'Cancer'
data.obj <- subset_data(data.obj, ind)
dist.obj <- subset_dist(dist.obj, ind)
save(data.obj, dist.obj, file = 'data.obj.wk.RData')


covars <- c("Batch","Bristol_score","BMI", "Age", "Sex", "GI_nonGI","Cancer_class","Metastasis","PPI_day_365", "Abx_day_365",
            "PPI_last_month","Abx_last_month","Charlson_score","Elix_score","Sample_season","Urban","icd10_first_3_name_short", "Site",
            "smoking_category","prior_chemotherapy")
blood.names <- c('Erythrocytes','Hematocrit','Neutrophils','MCV','Hemoglobin','Leukocytes','Platelet.Count',
                 'Neutrophils_cat','Platelet.Count_cat','Hemoglobin_cat',
                 "neut_neutropenia2","bone_marrow_suppression2","Hb_anemia2","Pl_thrombocytopenia2","neut_neutropenia_c","Hb_anemia_c","Pl_thrombocytopenia_c")
elix.names <- c("Elixhauser_CHF","Elixhauser_Arrhythmia","Elixhauser_Valvular","Elixhauser_PHTN","Elixhauser_PVD","Elixhauser_HTN","Elixhauser_Paralysis",  
                "Elixhauser_NeuroOther","Elixhauser_Pulmonary","Elixhauser_DM","Elixhauser_DMcx","Elixhauser_Hypothyroid","Elixhauser_Renal","Elixhauser_Liver",      
                "Elixhauser_PUD" ,"Elixhauser_HIV","Elixhauser_Lymphoma","Elixhauser_Mets","Elixhauser_Tumor","Elixhauser_Rheumatic","Elixhauser_Coagulopathy",
                "Elixhauser_Obesity","Elixhauser_WeightLoss","Elixhauser_FluidsLytes","Elixhauser_BloodLoss","Elixhauser_Anemia","Elixhauser_Alcohol","Elixhauser_Drugs",
                "Elixhauser_Psychoses","Elixhauser_Depression")

variables <- c(covars[!(covars %in% c("Batch"))], blood.names, elix.names)

dir <- paste0('CancerOnly_func/',func.type)

for(variable in variables){
  with(list(commandArgs = function(...) c("--args", variable, dir)), source(func_script))
}

# setwd("/Users/luyang1/myicloud/Documents/Mayo_project/2023_01_09_Oncobiome/mforge_clean/Result/CancerOnly_func/pathway/DAA/")
# dirs <- list.dirs(recursive = F,full.names = F)
# xx <- NULL
# for(dir in dirs){
#   diff.obj <- NULL
#   load(paste0(dir,'/',dir,'_ZicoSeq.Rdata'))
#   x <- colnames(diff.obj$coef.list$pathway)[grep(dir,colnames(diff.obj$coef.list$pathway))]
#   if(sum(grepl('1$',x))>0) xx <- c(xx, dir)
# }
# 
# variables <- xx

############################## Cancer Only (smoking for lung cancer) func #######################################
file_dir <- dirname(rstudioapi::getSourceEditorContext()$path)
file_dir <- dirname(dirname(dirname(dirname(file_dir))))

wd <- file_dir
rd <- paste0(file_dir,'/Result/')

setwd(wd)
source(paste0(file_dir,"/Code/Submission/MayoOncobiomeStudy/Code/Stats.R"))
try(load_package())

getwd()
load(file = 'Data/data.obj.pathway.RData') 


setwd(rd)
dir <- 'LungCancer_Smoking_func'
if(!dir.exists(dir)){dir.create(dir)}
setwd(dir)
getwd()

ind <- rownames(data.obj$meta.dat[data.obj$meta.dat$icd10_first_3_name_short == 'bronchus lung',])
data.obj <- subset_data(data.obj, ind)
dist.obj <- subset_dist(dist.obj, ind)
identical(rownames(data.obj$meta.dat), rownames(dist.obj$BC))
identical(rownames(data.obj$meta.dat), colnames(dist.obj$BC))
data.obj$meta.dat$smoking_category2 <- data.obj$meta.dat$smoking_category
data.obj$meta.dat$smoking_category2[data.obj$meta.dat$smoking_category2 %in% c("Current","Former")] <- "Yes"
save(data.obj, dist.obj, file = 'data.obj.wk.RData')

identical(colnames(dist.obj$Euc),colnames(data.obj$otu.tab))

variables <- c("smoking_category","smoking_category2")
dir <- 'LungCancer_Smoking_func'
for(variable in variables){
  with(list(commandArgs = function(...) c("--args", variable, dir)), source(func_script))
}



############################## sub Cancer X vs [Cancers-sub Cancer X] func #######################################
file_dir <- dirname(rstudioapi::getSourceEditorContext()$path)
file_dir <- dirname(dirname(dirname(dirname(file_dir))))

wd <- file_dir
rd <- paste0(file_dir,'/Result/subCancerX-Ex_func/')
if(!dir.exists(rd)){dir.create(rd)}
setwd(wd)
source(paste0(file_dir,"/Code/Submission/MayoOncobiomeStudy/Code/Stats.R"))
try(load_package())
func.type <- 'pathway'
setwd(wd)
load(file = paste0('Data/data.obj.',func.type,'.RData')) 

data.obj2 <- data.obj
dist.obj2 <- dist.obj

## find out the cancer types with N>15
data.obj$meta.dat$icd10_first_3_name_short <- as.character(data.obj$meta.dat$icd10_first_3_name_short)
samIDs <- rownames(data.obj$meta.dat[data.obj$meta.dat$icd10_first_3_name_short %in% names(which(sort(table(data.obj$meta.dat$icd10_first_3_name_short))>15)) & data.obj$meta.dat$Group != 'Control',])
data.obj <- subset_data(data.obj, samIDs = samIDs)
dist.obj <- subset_dist(dist.obj, samIDs = samIDs)
data.obj$meta.dat$icd10_first_3_name_short <- as.factor(data.obj$meta.dat$icd10_first_3_name_short)
cancer.type <- as.vector(unique(data.obj$meta.dat$icd10_first_3_name_short))
length(cancer.type)

for(i in 1:length(cancer.type)){
  setwd(rd)
  if(!dir.exists(func.type)){dir.create(func.type)}
  setwd(func.type)
  cancer <- cancer.type[i]
  cat('[',cancer,']\n')
  dir <- cancer.type[i]
  if(!dir.exists(dir)){dir.create(dir)}
  setwd(dir)
  getwd()
  
  data.obj <- data.obj2
  dist.obj <- dist.obj2

  data.obj$meta.dat$icd10_first_3_name_short <- as.character(data.obj$meta.dat$icd10_first_3_name_short)
  samIDs <- rownames(data.obj$meta.dat[data.obj$meta.dat$icd10_first_3_name_short %in% names(which(sort(table(data.obj$meta.dat$icd10_first_3_name_short))>15)) & data.obj$meta.dat$Group != 'Control',])
  data.obj <- subset_data(data.obj, samIDs = samIDs)
  dist.obj <- subset_dist(dist.obj, samIDs = samIDs)
  data.obj$meta.dat$icd10_first_3_name_short <- as.factor(data.obj$meta.dat$icd10_first_3_name_short)
  
  idx <- !(data.obj$meta.dat$icd10_first_3_name_short %in% c(cancer))
  data.obj$meta.dat$icd10_first_3_name_short <- as.character(data.obj$meta.dat$icd10_first_3_name_short)
  data.obj$meta.dat$icd10_first_3_name_short[idx] <- 'Others'
  
  dist.obj <- subset_dist(dist.obj, samIDs = rownames(data.obj$meta.dat))
  data.obj$meta.dat$icd10_first_3_name_short <- factor(data.obj$meta.dat$icd10_first_3_name_short,levels=c('Others',cancer))
  
  save(data.obj, dist.obj, file = 'data.obj.wk.RData')
  gc()
}


getwd()
variable <- c("icd10_first_3_name_short")
for(dir in paste0('subCancerX-Ex_func/',func.type,'/',cancer.type)){
  with(list(commandArgs = function(...) c("--args", variable, dir)), source(func_script))
}
  
# setwd("/Users/luyang1/myicloud/Documents/Mayo_project/2023_01_09_Oncobiome/mforge_clean/Result/subCancerX-Ex_func/pathway/")
# dirs <- list.dirs(recursive = F,full.names = F)
# xx <- NULL
# for(dir in dirs){
#   diff.obj <- NULL
#   load(paste0(dir,'/DAA/icd10_first_3_name_short/icd10_first_3_name_short_ZicoSeq.Rdata'))
#   x <- colnames(diff.obj$coef.list$pathway)[grep('icd10_first_3_name_short',colnames(diff.obj$coef.list$pathway))]
#   if(grepl('1$',x)) xx <- c(xx, dir)
# }
# xx

############################## sub Cancer X vs Control#######################################
file_dir <- dirname(rstudioapi::getSourceEditorContext()$path)
file_dir <- dirname(dirname(dirname(dirname(file_dir))))

wd <- file_dir
rd <- paste0(file_dir,'/Result/subCancerX_Control_func/')

if(!dir.exists(rd)) dir.create(rd)

setwd(wd)
source(paste0(file_dir,"/Code/Submission/MayoOncobiomeStudy/Code/Stats.R"))
try(load_package())
func.type <- 'pathway'
setwd(wd)
load(file = paste0('Data/data.obj.',func.type,'.RData')) 

data.obj2 <- data.obj
dist.obj2 <- dist.obj

data.obj$meta.dat$icd10_first_3_name_short <- as.character(data.obj$meta.dat$icd10_first_3_name_short)
samIDs <- rownames(data.obj$meta.dat[data.obj$meta.dat$icd10_first_3_name_short %in% names(which(sort(table(data.obj$meta.dat$icd10_first_3_name_short))>15)),])
data.obj <- subset_data(data.obj, samIDs = samIDs)

cancer.type <- as.vector(unique(data.obj$meta.dat$icd10_first_3_name_short)[!is.na(unique(data.obj$meta.dat$icd10_first_3_name_short))])
cancer.type <- cancer.type[!(cancer.type %in% 'healthy')]

for(i in 1:length(cancer.type)){
  setwd(rd)
  if(!dir.exists(func.type)){dir.create(func.type)}
  setwd(func.type)
  cancer <- cancer.type[i]
  cat('[',cancer,']\n')
  dir <- cancer.type[i]
  if(!dir.exists(dir)){dir.create(dir)}
  setwd(dir)
  getwd()
  
  data.obj <- data.obj2
  dist.obj <- dist.obj2
  
  data.obj$meta.dat$icd10_first_3_name_short <- as.character(data.obj$meta.dat$icd10_first_3_name_short)
  samIDs <- rownames(data.obj$meta.dat)[data.obj$meta.dat$icd10_first_3_name_short %in% c(cancer,'healthy')]
  data.obj <- subset_data(data.obj, samIDs = samIDs)
  dist.obj <- subset_dist(dist.obj, samIDs = samIDs)
  # data.obj$meta.dat$icd10_first_3_name_short <- as.factor(data.obj$meta.dat$icd10_first_3_name_short)
  data.obj$meta.dat$icd10_first_3_name_short <- factor(data.obj$meta.dat$icd10_first_3_name_short,levels=c('healthy',cancer))
  
  save(data.obj, dist.obj, file = 'data.obj.wk.RData')
  gc()
}


variable <- c("icd10_first_3_name_short")
for(dir in paste0('subCancerX_Control_func/',func.type,'/',cancer.type)){
  with(list(commandArgs = function(...) c("--args", variable, dir)), source(func_script))
}

# setwd("/Users/luyang1/myicloud/Documents/Mayo_project/2023_01_09_Oncobiome/mforge_clean/Result/subCancerX_Control_func/pathway/")
# dirs <- list.dirs(recursive = F,full.names = F)
# xx <- NULL
# for(dir in dirs){
#   diff.obj <- NULL
#   load(paste0(dir,'/DAA/icd10_first_3_name_short/icd10_first_3_name_short_ZicoSeq.Rdata'))
#   x <- colnames(diff.obj$coef.list$pathway)[grep('icd10_first_3_name_short',colnames(diff.obj$coef.list$pathway))]
#   if(grepl('1$',x)) xx <- c(xx, dir)
# }
# cancer.type <- cancer.type[cancer.type %in% gsub('Control-','',xx)]

############################## early onset func #######################################
file_dir <- dirname(rstudioapi::getSourceEditorContext()$path)
file_dir <- dirname(dirname(dirname(dirname(file_dir))))

wd <- file_dir
rd <- paste0(file_dir,'/Result/')

setwd(wd)
source(paste0(file_dir,"/Code/Submission/MayoOncobiomeStudy/Code/Stats.R"))
try(load_package())

setwd(rd)
dir <- 'EarlyOnset_func'
if(!dir.exists(dir)){dir.create(dir)}
setwd(dir)
getwd()
load(file = paste0('../../Data/data.obj.pathway.RData')) 

func.type <- 'pathway'

tb <- cbind(table(data.obj$meta.dat$icd10_first_3_name_short, data.obj$meta.dat$early_onset))
idx <- names(which(rowSums(tb >15) ==2))
idx <- idx[idx!='healthy']
for(j in idx){
  setwd(paste0(rd,'EarlyOnset_func'))
  load('../CancerOnly_func/pathway/data.obj.wk.RData')

  if(!dir.exists(j)){dir.create(j)}
  setwd(j)
  getwd()
  
  ind2 <- data.obj$meta.dat$icd10_first_3_name_short == j
  ind2 <- rownames(data.obj$meta.dat[ind2,])
  data.obj <- subset_data(data.obj, ind2)
  dist.obj <- subset_dist(dist.obj, ind2)
  
  save(data.obj, dist.obj, file = 'data.obj.wk.RData')
  getwd()
  setwd('..')
  
}


variable <- 'early_onset'
for(dir in paste0('EarlyOnset_func/',idx)){
  with(list(commandArgs = function(...) c("--args", variable, dir)), source(func_script))
}




############################## Pan-Cancer func #######################################
file_dir <- dirname(rstudioapi::getSourceEditorContext()$path)
file_dir <- dirname(dirname(dirname(dirname(file_dir))))

wd <- file_dir
rd <- paste0(file_dir,'/Result/')

# if(!dir.exists(rd)) dir.create(rd)
# setwd(wd)
# source(paste0(file_dir,"/Code/Submission/MayoOncobiomeStudy/Code/Stats.R"))
# try(load_package())
# load(file = 'Data/data.obj.raw.core.RData') 

setwd(wd)
source(paste0(file_dir,"/Code/Submission/MayoOncobiomeStudy/Code/Stats.R"))
try(load_package())
func.type <- 'pathway'
setwd(wd)
load(file = paste0('Data/data.obj.',func.type,'.RData')) 


setwd(rd)
dir <- 'PanCancer_func'
if(!dir.exists(dir)){dir.create(dir)}
setwd(dir)
getwd()

data.obj$meta.dat$Cancer_class <- as.character(data.obj$meta.dat$Cancer_class)
data.obj$meta.dat$Cancer_class[data.obj$meta.dat$Group=='Control'] <- 'Control'

data.obj$meta.dat$Group <- factor(data.obj$meta.dat$Group,levels=c('Control','Cancer'))

save(data.obj, dist.obj, file = 'data.obj.wk.RData')

variable <- "Group"; dir <- 'PanCancer_func'
with(list(commandArgs = function(...) c("--args", variable, dir)), source(func_script))








