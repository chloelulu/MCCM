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

## As I am not able to access mforge, run local with Alpha_Beta_DAA_clean.R
variable <- 'chemotherapy_last_2_years'
wd <- '/Users/luyang1//myicloud/Documents/Mayo_project/2023_01_09_Oncobiome/mforge_clean/'
rd <- '/Users/luyang1//myicloud/Documents/Mayo_project/2023_01_09_Oncobiome/mforge_clean/Result/'

adj.name <- covars[!(covars %in% c(variable,'GI_nonGI',"Charlson_score","Abx_last_month","PPI_last_month"))]

############################## Cancer Only #######################################
wd <- '/research/bsi/projects/staff_analysis/m216453/2023_01_09_Oncobiome/mforge_clean/'
rd <- '/research/bsi/projects/staff_analysis/m216453/2023_01_09_Oncobiome/mforge_clean/Result/'
setwd(wd)
source('Code/Stats.R')
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


## transform blood related variables into normal distributed
data.obj$meta.dat[c('Erythrocytes','Hematocrit','Neutrophils','MCV','Hemoglobin','Leukocytes','Platelet.Count')] <- as.data.frame(sapply(data.obj$meta.dat[,c('Erythrocytes','Hematocrit','Neutrophils','MCV','Hemoglobin','Leukocytes','Platelet.Count')], as.numeric)) %>% 
  `rownames<-`(rownames(data.obj$meta.dat)) 
data.obj.rff$meta.dat[c('Erythrocytes','Hematocrit','Neutrophils','MCV','Hemoglobin','Leukocytes','Platelet.Count')] <- as.data.frame(sapply(data.obj.rff$meta.dat[,c('Erythrocytes','Hematocrit','Neutrophils','MCV','Hemoglobin','Leukocytes','Platelet.Count')], as.numeric)) %>% 
  `rownames<-`(rownames(data.obj.rff$meta.dat))  

data.obj$meta.dat[c('Erythrocytes','Hematocrit','Neutrophils','MCV','Hemoglobin','Leukocytes','Platelet.Count')] <- apply(data.obj$meta.dat[c('Erythrocytes','Hematocrit','Neutrophils','MCV','Hemoglobin','Leukocytes','Platelet.Count')], 2, function(x){
  x[!is.na(x)] <- winzor(x[!is.na(x)], winsor.end = 'both')
  return(x)
})
data.obj.rff$meta.dat[c('Erythrocytes','Hematocrit','Neutrophils','MCV','Hemoglobin','Leukocytes','Platelet.Count')] <- apply(data.obj.rff$meta.dat[c('Erythrocytes','Hematocrit','Neutrophils','MCV','Hemoglobin','Leukocytes','Platelet.Count')], 2, function(x){
  x[!is.na(x)] <- winzor(x[!is.na(x)], winsor.end = 'both')
  return(x)
})

data.obj$meta.dat <- data.obj$meta.dat %>% mutate(Neutrophils = log2(Neutrophils), Leukocytes = log2(Leukocytes), Platelet.Count = log2(Platelet.Count))
data.obj.rff$meta.dat <- data.obj.rff$meta.dat %>% mutate(Neutrophils = log2(Neutrophils), Leukocytes = log2(Leukocytes), Platelet.Count = log2(Platelet.Count))

## transform blood related variables into binary format
data.obj$meta.dat$Neutrophils_cat <- as.factor(ifelse(is.na(data.obj$meta.dat$Neutrophils_cat),NA, ifelse(data.obj$meta.dat$Neutrophils_cat =='No','No','Yes')))
data.obj$meta.dat$Platelet.Count_cat <- as.factor(ifelse(is.na(data.obj$meta.dat$Platelet.Count_cat),NA, ifelse(data.obj$meta.dat$Platelet.Count_cat =='No','No','Yes')))
data.obj$meta.dat$Hemoglobin_cat <- as.factor(ifelse(is.na(data.obj$meta.dat$Hemoglobin_cat),NA, ifelse(data.obj$meta.dat$Hemoglobin_cat =='No','No','Yes')))

data.obj.rff$meta.dat$Neutrophils_cat <- as.factor(ifelse(is.na(data.obj.rff$meta.dat$Neutrophils_cat),NA, ifelse(data.obj.rff$meta.dat$Neutrophils_cat =='No','No','Yes')))
data.obj.rff$meta.dat$Platelet.Count_cat <- as.factor(ifelse(is.na(data.obj.rff$meta.dat$Platelet.Count_cat),NA, ifelse(data.obj.rff$meta.dat$Platelet.Count_cat =='No','No','Yes')))
data.obj.rff$meta.dat$Hemoglobin_cat <- as.factor(ifelse(is.na(data.obj.rff$meta.dat$Hemoglobin_cat),NA, ifelse(data.obj.rff$meta.dat$Hemoglobin_cat =='No','No','Yes')))

elix.names <- colnames(data.obj$meta.dat)[grep('Elixhauser_',colnames(data.obj$meta.dat))]
elix.names <- elix.names[!(elix.names %in% c("Elixhauser_dtindex","Elixhauser_elix.sum"))]

data.obj$meta.dat$GI_nonGI <- factor(data.obj$meta.dat$GI_nonGI,levels=c('non_GI','GI'))
data.obj.rff$meta.dat$GI_nonGI <- factor(data.obj.rff$meta.dat$GI_nonGI,levels=c('non_GI','GI'))


save(data.obj, dist.obj, data.obj.rff, dist.obj.rff, file = 'data.obj.wk.RData')


covars <- c("Batch","Bristol_score","BMI", "Age", "Sex", "GI_nonGI","Cancer_class","Metastasis","PPI_day_365", "Abx_day_365", 
            "PPI_last_month","Abx_last_month","Charlson_score","Elix_score","Sample_season","Urban" ,"icd10_first_3_name")
blood.names <- c('Erythrocytes','Hematocrit','Neutrophils','MCV','Hemoglobin','Leukocytes','Platelet.Count',
                 'Neutrophils_cat','Platelet.Count_cat','Hemoglobin_cat',
                 "neut_neutropenia2","bone_marrow_suppression2","Hb_anemia2","Pl_thrombocytopenia2","neut_neutropenia_c","Hb_anemia_c","Pl_thrombocytopenia_c")
elix.names <- c("Elixhauser_CHF","Elixhauser_Arrhythmia","Elixhauser_Valvular","Elixhauser_PHTN","Elixhauser_PVD","Elixhauser_HTN","Elixhauser_Paralysis",  
                "Elixhauser_NeuroOther","Elixhauser_Pulmonary","Elixhauser_DM","Elixhauser_DMcx","Elixhauser_Hypothyroid","Elixhauser_Renal","Elixhauser_Liver",      
                "Elixhauser_PUD" ,"Elixhauser_HIV","Elixhauser_Lymphoma","Elixhauser_Mets","Elixhauser_Tumor","Elixhauser_Rheumatic","Elixhauser_Coagulopathy",
                "Elixhauser_Obesity","Elixhauser_WeightLoss","Elixhauser_FluidsLytes","Elixhauser_BloodLoss","Elixhauser_Anemia","Elixhauser_Alcohol","Elixhauser_Drugs",
                "Elixhauser_Psychoses","Elixhauser_Depression")
charlson.names <- c("Charlson_MI","Charlson_CHF","Charlson_PVD","Charlson_Stroke","Charlson_Dementia","Charlson_Pulmonary","Charlson_Rheumatic",  
                    "Charlson_PUD","Charlson_LiverMild","Charlson_DM","Charlson_DMcx","Charlson_Paralysis","Charlson_Renal","Charlson_Cancer",     
                    "Charlson_LiverSevere","Charlson_Mets","Charlson_HIV")

variables <- c(covars, blood.names, elix.names, charlson.names)
dir <- 'CancerOnly'

for(variable in variables){
  sh <- paste(
    "sbatch",
    paste0("-J ", variable),
    "--partition=cpu-short",
    paste("--output=/research/bsi/projects/staff_analysis/m216453/2023_01_09_Oncobiome/trash/OUT/",  variable, "_",dir,".out", sep=""),
    paste("--error=/research/bsi/projects/staff_analysis/m216453/2023_01_09_Oncobiome/trash/OUT/",  variable,"_",dir,".err", sep=""),
    "--time=1-00:00:00",
    "--mem=64G",
    "--export=ALL",
    paste("\'--wrap=R CMD BATCH --no-restore \"--args ", variable," ",dir,"\" ",
          "/research/bsi/projects/staff_analysis/m216453/2023_01_09_Oncobiome/Code/Alpha_Beta_DAA_clean.R", " ", paste0("/research/bsi/projects/staff_analysis/m216453/2023_01_09_Oncobiome/trash/",variable,"_",dir,".Rout"), "\'", sep="")
  )
  print(sh)
  system(sh)
}



############################## Control vs sub Cancer X #######################################
wd <- '/research/bsi/projects/staff_analysis/m216453/2023_01_09_Oncobiome/'
rd <- '/research/bsi/projects/staff_analysis/m216453/2023_01_09_Oncobiome/Result/subCancerX_Control/'
setwd(wd)
source('Code/Stats.R')
try(load_package())
load(file = 'Data/data.obj.raw.core.RData') 
data.obj2 <- data.obj
data.obj.rff2 <- data.obj.rff
dist.obj2 <- dist.obj
dist.obj.rff2 <- dist.obj.rff

data.obj2$meta.dat$icd10_first_3_name <- as.character(data.obj2$meta.dat$icd10_first_3_name)
data.obj2$meta.dat$icd10_first_3_name[data.obj2$meta.dat$Group=='Control'] <- 'Control'
samIDs <- rownames(data.obj2$meta.dat[data.obj2$meta.dat$icd10_first_3_name %in% names(which(sort(table(data.obj2$meta.dat$icd10_first_3_name))>15)),])
data.obj2 <- subset_data(data.obj2, samIDs = samIDs)
dist.obj2 <- subset_dist(dist.obj2, samIDs = samIDs)


data.obj.rff2$meta.dat$icd10_first_3_name <- as.character(data.obj.rff2$meta.dat$icd10_first_3_name)
data.obj.rff2$meta.dat$icd10_first_3_name[data.obj.rff2$meta.dat$Group=='Control'] <- 'Control'
samIDs <- rownames(data.obj.rff2$meta.dat[data.obj2$meta.dat$icd10_first_3_name %in% names(which(sort(table(data.obj2$meta.dat$icd10_first_3_name))>15)),])
data.obj.rff2 <- subset_data(data.obj.rff2, samIDs = samIDs)
dist.obj.rff2 <- subset_dist(dist.obj.rff2, samIDs = samIDs)


data.obj2$meta.dat$icd10_first_3_name <- as.factor(data.obj2$meta.dat$icd10_first_3_name)
data.obj.rff2$meta.dat$icd10_first_3_name <- as.factor(data.obj.rff2$meta.dat$icd10_first_3_name)
table(data.obj2$meta.dat$icd10_first_3_name, data.obj2$meta.dat$Batch)

cancer.type <- as.vector(unique(data.obj2$meta.dat$icd10_first_3_name)[!is.na(unique(data.obj2$meta.dat$icd10_first_3_name))])
cancer.type <- cancer.type[cancer.type != 'Control']
# cancer.type <- cancer.type[-grep('ther',cancer.type)]
# cancer.type <- cancer.type[grep('ther',cancer.type)] # 2024/10/31 Ruben requests to add 2 cancer types with names "other" back, we run this line of code for adding. However, when submission, exclude this line and one line above means all cancer types>15 included

cancer.dir <- gsub('Malignant neoplasm of |Multiple myeloma and malignant |Malignant ','',cancer.type)
cancer.dir <- gsub('Other\\ and\\ |malignant\\ ','',cancer.dir)## add 2024/10/31 for adding 2 other cancer back
cancer.dir <- gsub('__','_',gsub(',|\\ ','_',cancer.dir))
cancer.dir <- gsub('eye_brain_and_other_parts_of_|_and','',cancer.dir)
setwd(rd)

for(i in 1:length(cancer.type)){
  cancer <- cancer.type[i]
  cat('[',cancer,']\n')
  dir <- paste0('Control-',cancer.dir[i])
  if(!dir.exists(dir)){dir.create(dir)}
  setwd(dir)
  getwd()
  
  ind <- data.obj2$meta.dat$icd10_first_3_name %in% c('Control',cancer)
  data.obj <- subset_data(data.obj2, ind)
  dist.obj <- subset_dist(dist.obj2, ind)
  cat(sum(ind),'\n')
  
  ind <- data.obj.rff2$meta.dat$icd10_first_3_name %in% c('Control',cancer)
  data.obj.rff <- subset_data(data.obj.rff2, ind)
  dist.obj.rff <- subset_dist(dist.obj.rff2, ind)
  cat(sum(ind),'\n')
  
  data.obj$meta.dat$icd10_first_3_name <- factor(data.obj$meta.dat$icd10_first_3_name,levels=c('Control',cancer))
  data.obj.rff$meta.dat$icd10_first_3_name <- factor(data.obj.rff$meta.dat$icd10_first_3_name,levels=c('Control',cancer))
  
  save(data.obj, dist.obj, data.obj.rff, dist.obj.rff, file = 'data.obj.wk.RData')
  setwd('..')
}


variables <- c("icd10_first_3_name")
for(dir in paste0('subCancerX_Control/Control-',cancer.dir)){
  for(variable in variables){
    sh <- paste(
      "sbatch",
      paste0("-J ", variable,'-',gsub('subCancerX_Control/','',dir)),
      "--partition=cpu-short",
      paste("--output=/research/bsi/projects/staff_analysis/m216453/2023_01_09_Oncobiome/trash/OUT/",  variable, "_",gsub('subCancerX_Control/','',dir),".out", sep=""),
      paste("--error=/research/bsi/projects/staff_analysis/m216453/2023_01_09_Oncobiome/trash/OUT/",  variable,"_",gsub('subCancerX_Control/','',dir),".err", sep=""),
      "--time=1-00:00:00",
      "--mem=64G",
      "--export=ALL",
      paste("\'--wrap=R CMD BATCH --no-restore \"--args ", variable," ",dir,"\" ",
            "/research/bsi/projects/staff_analysis/m216453/2023_01_09_Oncobiome/Code/Alpha_Beta_DAA_clean.R", " ", paste0("/research/bsi/projects/staff_analysis/m216453/2023_01_09_Oncobiome/trash/",variable,"_",gsub('subCancerX_Control/','',dir),".Rout"), "\'", sep="")
    )
    print(sh)
    system(sh)
  }
}



############################## sub Cancer X vs [Cancers-sub Cancer X] #######################################
wd <- '/research/bsi/projects/staff_analysis/m216453/2023_01_09_Oncobiome/'
rd <- '/research/bsi/projects/staff_analysis/m216453/2023_01_09_Oncobiome/Result/subCancerX-Ex/'

setwd(wd)
source('Code/Stats.R')
try(load_package())

load(file = 'Data/data.obj.raw.core.RData') 
data.obj2 <- data.obj
data.obj.rff2 <- data.obj.rff
dist.obj2 <- dist.obj
dist.obj.rff2 <- dist.obj.rff

data.obj$meta.dat$icd10_first_3_name <- as.character(data.obj$meta.dat$icd10_first_3_name)

idx <- names(which(sort(table(data.obj$meta.dat$icd10_first_3_name))>15)) 
samIDs <- rownames(data.obj$meta.dat[data.obj$meta.dat$icd10_first_3_name %in% idx & data.obj$meta.dat$Group != 'Control',])
data.obj <- subset_data(data.obj, samIDs = samIDs)
dist.obj <- subset_dist(dist.obj, samIDs = samIDs)

data.obj.rff$meta.dat$icd10_first_3_name <- as.character(data.obj.rff$meta.dat$icd10_first_3_name)
samIDs <- rownames(data.obj.rff$meta.dat[data.obj.rff$meta.dat$icd10_first_3_name %in% idx & data.obj.rff$meta.dat$Group != 'Control',])
data.obj.rff <- subset_data(data.obj.rff, samIDs = samIDs)
dist.obj.rff <- subset_dist(dist.obj.rff, samIDs = samIDs)

data.obj$meta.dat$icd10_first_3_name <- as.factor(data.obj$meta.dat$icd10_first_3_name)
data.obj.rff$meta.dat$icd10_first_3_name <- as.factor(data.obj.rff$meta.dat$icd10_first_3_name)

cancer.type <- as.vector(unique(data.obj$meta.dat$icd10_first_3_name)[!is.na(unique(data.obj$meta.dat$icd10_first_3_name))])

cancer.dir <- gsub('Malignant neoplasm of |Multiple myeloma and malignant |Malignant ','',cancer.type)
cancer.dir <- gsub('Other\\ and\\ |malignant\\ ','',cancer.dir)
cancer.dir <- gsub('__','_',gsub(',|\\ ','_',cancer.dir))
cancer.dir <- gsub('eye_brain_and_other_parts_of_|_and','',cancer.dir)
setwd(rd)

for(i in 1:length(cancer.type)){
  cancer <- cancer.type[i]
  cat('[',cancer,']\n')
  dir <- cancer.dir[i]
  if(!dir.exists(dir)){dir.create(dir)}
  setwd(dir)
  getwd()
  
  data.obj <- data.obj2
  data.obj.rff <- data.obj.rff2
  dist.obj <- dist.obj2
  dist.obj.rff <- dist.obj.rff2
  
  idx <- names(which(sort(table(data.obj$meta.dat$icd10_first_3_name))>15)) 

  samIDs <- rownames(data.obj$meta.dat[data.obj$meta.dat$icd10_first_3_name %in% idx & data.obj$meta.dat$Group != 'Control',])
  data.obj <- subset_data(data.obj, samIDs = samIDs)
  dist.obj <- subset_dist(dist.obj, samIDs = samIDs)
  
  data.obj.rff$meta.dat$icd10_first_3_name <- as.character(data.obj.rff$meta.dat$icd10_first_3_name)
  samIDs <- rownames(data.obj.rff$meta.dat[data.obj.rff$meta.dat$icd10_first_3_name %in% idx & data.obj.rff$meta.dat$Group != 'Control',])
  data.obj.rff <- subset_data(data.obj.rff, samIDs = samIDs)
  dist.obj.rff <- subset_dist(dist.obj.rff, samIDs = samIDs)
  
  data.obj$meta.dat$icd10_first_3_name <- as.character(data.obj$meta.dat$icd10_first_3_name)

  cancer.type <- as.vector(unique(data.obj$meta.dat$icd10_first_3_name)[!is.na(unique(data.obj$meta.dat$icd10_first_3_name))])
  
  samIDs <- rownames(data.obj$meta.dat[data.obj$meta.dat$icd10_first_3_name %in% cancer.type & data.obj$meta.dat$Group != 'Control',])
  data.obj <- subset_data(data.obj, samIDs = samIDs)
  dist.obj <- subset_dist(dist.obj, samIDs = samIDs)
  
  samIDs <- rownames(data.obj.rff$meta.dat[data.obj.rff$meta.dat$icd10_first_3_name %in% cancer.type & data.obj.rff$meta.dat$Group != 'Control',])
  data.obj.rff <- subset_data(data.obj.rff, samIDs = samIDs)
  dist.obj.rff <- subset_dist(dist.obj.rff, samIDs = samIDs)
  
  data.obj$meta.dat$icd10_first_3_name <- as.factor(data.obj$meta.dat$icd10_first_3_name)
  data.obj.rff$meta.dat$icd10_first_3_name <- as.factor(data.obj.rff$meta.dat$icd10_first_3_name)
  
  idx <- !(data.obj$meta.dat$icd10_first_3_name %in% c(cancer))
  data.obj$meta.dat$icd10_first_3_name <- as.character(data.obj$meta.dat$icd10_first_3_name)
  data.obj$meta.dat$icd10_first_3_name[idx] <- 'Others'
  dist.obj <- subset_dist(dist.obj, samIDs = rownames(data.obj$meta.dat))
  
  idx <- !(data.obj.rff$meta.dat$icd10_first_3_name %in% c(cancer))
  data.obj.rff$meta.dat$icd10_first_3_name <- as.character(data.obj.rff$meta.dat$icd10_first_3_name)
  data.obj.rff$meta.dat$icd10_first_3_name[idx] <- 'Others'
  dist.obj.rff <- subset_dist(dist.obj.rff, samIDs = rownames(data.obj.rff$meta.dat))
  cat(sum(idx),'\n')
  table(data.obj.rff$meta.dat$icd10_first_3_name)
  
  data.obj$meta.dat$icd10_first_3_name <- factor(data.obj$meta.dat$icd10_first_3_name,levels=c('Others',cancer))
  data.obj.rff$meta.dat$icd10_first_3_name <- factor(data.obj.rff$meta.dat$icd10_first_3_name,levels=c('Others',cancer))
  
  save(data.obj, dist.obj, data.obj.rff, dist.obj.rff, file = 'data.obj.wk.RData')
  setwd('..')
  gc()
}



variables <- c("icd10_first_3_name")
for(dir in paste0('subCancerX-Ex/',cancer.dir)){
  for(variable in variables){
    sh <- paste(
      "sbatch",
      paste0("-J ", variable,'-',gsub('subCancerX-Ex/','',dir)),
      "--partition=cpu-short",
      paste("--output=/research/bsi/projects/staff_analysis/m216453/2023_01_09_Oncobiome/trash/OUT/",  variable, "_",gsub('subCancerX-Ex/','',dir),".out", sep=""),
      paste("--error=/research/bsi/projects/staff_analysis/m216453/2023_01_09_Oncobiome/trash/OUT/",  variable,"_",gsub('subCancerX-Ex/','',dir),".err", sep=""),
      "--time=1-00:00:00",
      "--mem=64G",
      "--export=ALL",
      paste("\'--wrap=R CMD BATCH --no-restore \"--args ", variable," ",dir,"\" ",
            "/research/bsi/projects/staff_analysis/m216453/2023_01_09_Oncobiome/Code/Alpha_Beta_DAA_clean.R", " ", paste0("/research/bsi/projects/staff_analysis/m216453/2023_01_09_Oncobiome/trash/",variable,"_",gsub('subCancerX-Ex/','',dir),".Rout"), "\'", sep="")
    )
    print(sh)
    system(sh)
  }
}






############################## Control vs sub Cancer X [PERMANOVA]##############################
wd <- '/research/bsi/projects/staff_analysis/m216453/2023_01_09_Oncobiome/'
rd <- '/research/bsi/projects/staff_analysis/m216453/2023_01_09_Oncobiome/Result/subCancerX_Control/'
fd <- '/research/bsi/projects/staff_analysis/m216453/2023_01_09_Oncobiome/Figure/subCancerX_Control/'

setwd(wd)
source('Code/Stats.R')
try(load_package())
load(file = 'Data/data.obj.raw.core.RData') 
data.obj2 <- data.obj
data.obj.rff2 <- data.obj.rff
dist.obj2 <- dist.obj
dist.obj.rff2 <- dist.obj.rff

data.obj$meta.dat$icd10_first_3_name <- as.character(data.obj$meta.dat$icd10_first_3_name)
data.obj$meta.dat$icd10_first_3_name[data.obj$meta.dat$Group == 'Control'] <- 'Control'
samIDs <- rownames(data.obj$meta.dat[data.obj$meta.dat$icd10_first_3_name %in% names(which(sort(table(data.obj$meta.dat$icd10_first_3_name))>15)),])
data.obj <- subset_data(data.obj, samIDs = samIDs)
dist.obj <- subset_dist(dist.obj, samIDs = samIDs)

data.obj.rff$meta.dat$icd10_first_3_name <- as.character(data.obj.rff$meta.dat$icd10_first_3_name)
data.obj.rff$meta.dat$icd10_first_3_name[data.obj.rff$meta.dat$Group == 'Control'] <- 'Control'
samIDs <- rownames(data.obj.rff$meta.dat[data.obj.rff$meta.dat$icd10_first_3_name %in% names(which(sort(table(data.obj.rff$meta.dat$icd10_first_3_name))>15)),])
data.obj.rff <- subset_data(data.obj.rff, samIDs = samIDs)
dist.obj.rff <- subset_dist(dist.obj.rff, samIDs = samIDs)

data.obj$meta.dat$icd10_first_3_name <- as.factor(data.obj$meta.dat$icd10_first_3_name)
data.obj.rff$meta.dat$icd10_first_3_name <- as.factor(data.obj.rff$meta.dat$icd10_first_3_name)

cancer.type <- as.vector(unique(data.obj$meta.dat$icd10_first_3_name)[!is.na(unique(data.obj$meta.dat$icd10_first_3_name))])

cancer.dir <- gsub('Malignant neoplasm of |Multiple myeloma and malignant |Malignant ','',cancer.type)
cancer.dir <- gsub('__','_',gsub(',|\\ ','_',cancer.dir))
cancer.dir <- gsub('eye_brain_and_other_parts_of_|_and','',cancer.dir)
setwd(rd)

dist.names <- c('UniFrac', 'GUniFrac', 'WUniFrac', 'BC')
## Use DMANOVA
r2.unadj.mat <- pv.unadj.mat <-
  array(NA, c(length(cancer.dir), length(cancer.dir),length(dist.names)),
        dimnames = list(cancer.dir, cancer.dir, dist.names))
variable <- 'icd10_first_3_name'
for(i in 1:(length(cancer.type)-1)){
  for(j in (i+1):length(cancer.type)){
    if(i != j ){
      ind <- data.obj.rff$meta.dat$icd10_first_3_name %in% c(cancer.type[i],cancer.type[j])
      data.obj.rff3 <- subset_data(data.obj.rff, ind)
      dist.obj.rff3 <- subset_dist(dist.obj.rff, ind)
      for(dist.name in dist.names){
        obj <- dmanova(as.dist(dist.obj.rff3[[dist.name]]) ~ data.obj.rff3$meta.dat[, variable])
        r2.unadj.mat[cancer.dir[i],cancer.dir[j],dist.name] <- obj$aov.tab[1,5]
        pv.unadj.mat[cancer.dir[i],cancer.dir[j],dist.name] <- obj$aov.tab[1,6]
      }
    }
  }
}
setwd(fd)
save(r2.unadj.mat, pv.unadj.mat, file = 'r2.p.unadj.mat.Rdata')


## Since DMANOVA does not support small smaple size(some subcancer type), I try permanova to see.
r2.unadj.mat <- pv.unadj.mat <-
  array(NA, c(length(cancer.dir), length(cancer.dir),length(dist.names)),
        dimnames = list(cancer.dir, cancer.dir, dist.names))
variable <- 'icd10_first_3_name'
for(i in 1:(length(cancer.type)-1)){
  for(j in (i+1):length(cancer.type)){
    if(i != j ){
      ind <- data.obj.rff$meta.dat$icd10_first_3_name %in% c(cancer.type[i],cancer.type[j])
      data.obj.rff3 <- subset_data(data.obj.rff, ind)
      dist.obj.rff3 <- subset_dist(dist.obj.rff, ind)
      for(dist.name in dist.names){
        obj <- adonis(as.dist(dist.obj.rff3[[dist.name]]) ~ data.obj.rff3$meta.dat[, variable])
        r2.unadj.mat[cancer.dir[i],cancer.dir[j],dist.name] <- obj$aov.tab[1,'R2']
        pv.unadj.mat[cancer.dir[i],cancer.dir[j],dist.name] <- obj$aov.tab[1,'Pr(>F)']
      }
    }
  }
}
setwd(fd)
save(r2.unadj.mat, pv.unadj.mat, file = 'r2.p.unadj.mat_permanova.Rdata')

##---------------------------- functional data analysis-----------------------------------
############################## Cancer Only func #######################################
wd <- '/research/bsi/projects/staff_analysis/m216453/2023_01_09_Oncobiome/'
rd <- '/research/bsi/projects/staff_analysis/m216453/2023_01_09_Oncobiome/Result/'
setwd(wd)
source('Code/Stats.R')
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

ind <- data.obj.rff$meta.dat$Group == 'Cancer'
data.obj.rff <- subset_data(data.obj.rff, ind)
dist.obj.rff <- subset_dist(dist.obj.rff, ind)

# par(mfrow = c(3,3))
# for(i in c('Erythrocytes','Hematocrit','Neutrophils','MCV','Hemoglobin','Leukocytes','Platelet.Count')){
#   print(hist(as.numeric(data.obj$meta.dat[,i])))
# }

## transform blood related variables into normal distributed
data.obj$meta.dat[c('Erythrocytes','Hematocrit','Neutrophils','MCV','Hemoglobin','Leukocytes','Platelet.Count')] <- as.data.frame(sapply(data.obj$meta.dat[,c('Erythrocytes','Hematocrit','Neutrophils','MCV','Hemoglobin','Leukocytes','Platelet.Count')], as.numeric)) %>% 
  `rownames<-`(rownames(data.obj$meta.dat)) 
data.obj.rff$meta.dat[c('Erythrocytes','Hematocrit','Neutrophils','MCV','Hemoglobin','Leukocytes','Platelet.Count')] <- as.data.frame(sapply(data.obj.rff$meta.dat[,c('Erythrocytes','Hematocrit','Neutrophils','MCV','Hemoglobin','Leukocytes','Platelet.Count')], as.numeric)) %>% 
  `rownames<-`(rownames(data.obj.rff$meta.dat))  

data.obj$meta.dat[c('Erythrocytes','Hematocrit','Neutrophils','MCV','Hemoglobin','Leukocytes','Platelet.Count')] <- apply(data.obj$meta.dat[c('Erythrocytes','Hematocrit','Neutrophils','MCV','Hemoglobin','Leukocytes','Platelet.Count')], 2, function(x){
  x[!is.na(x)] <- winzor(x[!is.na(x)], winsor.end = 'both')
  return(x)
})
data.obj.rff$meta.dat[c('Erythrocytes','Hematocrit','Neutrophils','MCV','Hemoglobin','Leukocytes','Platelet.Count')] <- apply(data.obj.rff$meta.dat[c('Erythrocytes','Hematocrit','Neutrophils','MCV','Hemoglobin','Leukocytes','Platelet.Count')], 2, function(x){
  x[!is.na(x)] <- winzor(x[!is.na(x)], winsor.end = 'both')
  return(x)
})

data.obj$meta.dat <- data.obj$meta.dat %>% mutate(Neutrophils = log2(Neutrophils), Leukocytes = log2(Leukocytes), Platelet.Count = log2(Platelet.Count))
data.obj.rff$meta.dat <- data.obj.rff$meta.dat %>% mutate(Neutrophils = log2(Neutrophils), Leukocytes = log2(Leukocytes), Platelet.Count = log2(Platelet.Count))

## transform blood related variables into binary format
data.obj$meta.dat$Neutrophils_cat <- as.factor(ifelse(is.na(data.obj$meta.dat$Neutrophils_cat),NA, ifelse(data.obj$meta.dat$Neutrophils_cat =='No','No','Yes')))
data.obj$meta.dat$Platelet.Count_cat <- as.factor(ifelse(is.na(data.obj$meta.dat$Platelet.Count_cat),NA, ifelse(data.obj$meta.dat$Platelet.Count_cat =='No','No','Yes')))
data.obj$meta.dat$Hemoglobin_cat <- as.factor(ifelse(is.na(data.obj$meta.dat$Hemoglobin_cat),NA, ifelse(data.obj$meta.dat$Hemoglobin_cat =='No','No','Yes')))

data.obj.rff$meta.dat$Neutrophils_cat <- as.factor(ifelse(is.na(data.obj.rff$meta.dat$Neutrophils_cat),NA, ifelse(data.obj.rff$meta.dat$Neutrophils_cat =='No','No','Yes')))
data.obj.rff$meta.dat$Platelet.Count_cat <- as.factor(ifelse(is.na(data.obj.rff$meta.dat$Platelet.Count_cat),NA, ifelse(data.obj.rff$meta.dat$Platelet.Count_cat =='No','No','Yes')))
data.obj.rff$meta.dat$Hemoglobin_cat <- as.factor(ifelse(is.na(data.obj.rff$meta.dat$Hemoglobin_cat),NA, ifelse(data.obj.rff$meta.dat$Hemoglobin_cat =='No','No','Yes')))

elix.names <- colnames(data.obj$meta.dat)[grep('Elixhauser_',colnames(data.obj$meta.dat))]
elix.names <- elix.names[!(elix.names %in% c("Elixhauser_dtindex","Elixhauser_elix.sum"))]

data.obj$meta.dat$GI_nonGI <- factor(data.obj$meta.dat$GI_nonGI,levels=c('non_GI','GI'))
data.obj.rff$meta.dat$GI_nonGI <- factor(data.obj.rff$meta.dat$GI_nonGI,levels=c('non_GI','GI'))


save(data.obj, dist.obj, data.obj.rff, dist.obj.rff, file = 'data.obj.wk.RData')


covars <- c("Batch","Bristol_score","BMI", "Age", "Sex", "GI_nonGI","Cancer_class","Metastasis","PPI_day_365", "Abx_day_365",
            "PPI_last_month","Abx_last_month","Charlson_score","Elix_score","Sample_season","Urban","icd10_first_3_name")
blood.names <- c('Erythrocytes','Hematocrit','Neutrophils','MCV','Hemoglobin','Leukocytes','Platelet.Count',
                 'Neutrophils_cat','Platelet.Count_cat','Hemoglobin_cat',
                 "neut_neutropenia2","bone_marrow_suppression2","Hb_anemia2","Pl_thrombocytopenia2","neut_neutropenia_c","Hb_anemia_c","Pl_thrombocytopenia_c")
elix.names <- c("Elixhauser_CHF","Elixhauser_Arrhythmia","Elixhauser_Valvular","Elixhauser_PHTN","Elixhauser_PVD","Elixhauser_HTN","Elixhauser_Paralysis",  
                "Elixhauser_NeuroOther","Elixhauser_Pulmonary","Elixhauser_DM","Elixhauser_DMcx","Elixhauser_Hypothyroid","Elixhauser_Renal","Elixhauser_Liver",      
                "Elixhauser_PUD" ,"Elixhauser_HIV","Elixhauser_Lymphoma","Elixhauser_Mets","Elixhauser_Tumor","Elixhauser_Rheumatic","Elixhauser_Coagulopathy",
                "Elixhauser_Obesity","Elixhauser_WeightLoss","Elixhauser_FluidsLytes","Elixhauser_BloodLoss","Elixhauser_Anemia","Elixhauser_Alcohol","Elixhauser_Drugs",
                "Elixhauser_Psychoses","Elixhauser_Depression")
charlson.names <- c("Charlson_MI","Charlson_CHF","Charlson_PVD","Charlson_Stroke","Charlson_Dementia","Charlson_Pulmonary","Charlson_Rheumatic",  
                    "Charlson_PUD","Charlson_LiverMild","Charlson_DM","Charlson_DMcx","Charlson_Paralysis","Charlson_Renal","Charlson_Cancer",     
                    "Charlson_LiverSevere","Charlson_Mets","Charlson_HIV")

variables <- c(covars, blood.names, elix.names)

dir <- paste0('CancerOnly_func/',func.type)
for(variable in variables){
  sh <- paste(
    "sbatch",
    paste0("-J ", variable),
    "--partition=cpu-short",
    paste("--output=/research/bsi/projects/staff_analysis/m216453/2023_01_09_Oncobiome/trash/OUT/",  variable, "_",gsub('\\/','',dir),".out", sep=""),
    paste("--error=/research/bsi/projects/staff_analysis/m216453/2023_01_09_Oncobiome/trash/OUT/",  variable,"_",gsub('\\/','',dir),".err", sep=""),
    "--time=1-00:00:00",
    "--mem=64G",
    "--export=ALL",
    paste("\'--wrap=R CMD BATCH --no-restore \"--args ", variable," ",dir,"\" ",
          "/research/bsi/projects/staff_analysis/m216453/2023_01_09_Oncobiome/Code/Alpha_Beta_DAA_func_clean.R", " ", paste0("/research/bsi/projects/staff_analysis/m216453/2023_01_09_Oncobiome/trash/",variable,"_",gsub('\\/','',dir),".Rout"), "\'", sep="")
  )
  
  print(sh)
  system(sh)
}






############################## sub Cancer X vs [Cancers-sub Cancer X] func #######################################
# wd <- '/research/bsi/projects/staff_analysis/m216453/2023_01_09_Oncobiome/'
# rd <- '/research/bsi/projects/staff_analysis/m216453/2023_01_09_Oncobiome/Result/subCancerX-Ex2_func/'
wd <- '/Users/luyang1//myicloud/Documents/Mayo_project/2023_01_09_Oncobiome/mforge_clean/'
rd <- '/Users/luyang1//myicloud/Documents/Mayo_project/2023_01_09_Oncobiome/mforge_clean/Result/subCancerX-Ex2_func/'
setwd(wd)
source('Code/Stats.R')
try(load_package())
for(func.type in c('pathway')){#'uniref90_ko','uniref90_go','uniref90_level4ec',
  setwd(wd)
  load(file = paste0('Data/data.obj.',func.type,'.RData')) 
  
  data.obj2 <- data.obj
  dist.obj2 <- dist.obj
  
  data.obj$meta.dat$icd10_first_3_name <- as.character(data.obj$meta.dat$icd10_first_3_name)
  samIDs <- rownames(data.obj$meta.dat[data.obj$meta.dat$icd10_first_3_name %in% names(which(sort(table(data.obj$meta.dat$icd10_first_3_name))>15)) & data.obj$meta.dat$Group != 'Control',])
  data.obj <- subset_data(data.obj, samIDs = samIDs)
  dist.obj <- subset_dist(dist.obj, samIDs = samIDs)
  data.obj$meta.dat$icd10_first_3_name <- as.factor(data.obj$meta.dat$icd10_first_3_name)
  cancer.type <- as.vector(unique(data.obj$meta.dat$icd10_first_3_name)[!is.na(unique(data.obj$meta.dat$icd10_first_3_name))])
  # cancer.type <- cancer.type[-grep('ther',cancer.type)] # exclude 2 cancer types Ruben mentioned not to include
  
  samIDs <- rownames(data.obj$meta.dat[data.obj$meta.dat$icd10_first_3_name %in% cancer.type & data.obj$meta.dat$Group != 'Control',])
  data.obj <- subset_data(data.obj, samIDs = samIDs)
  dist.obj <- subset_dist(dist.obj, samIDs = samIDs)
  
  cancer.dir <- gsub('Malignant neoplasm of |Multiple myeloma and malignant |Malignant ','',cancer.type)
  cancer.dir <- gsub('Other\\ and\\ |malignant\\ ','',cancer.dir)## add 2024/10/31 for adding 2 other cancer back
  cancer.dir <- gsub('__','_',gsub(',|\\ ','_',cancer.dir))
  cancer.dir <- gsub('eye_brain_and_other_parts_of_|_and','',cancer.dir)
  
  
  for(i in 1:length(cancer.type)){
    setwd(rd)
    if(!dir.exists(func.type)){dir.create(func.type)}
    setwd(func.type)
    cancer <- cancer.type[i]
    cat('[',cancer,']\n')
    dir <- cancer.dir[i]
    if(!dir.exists(dir)){dir.create(dir)}
    setwd(dir)
    getwd()
    
    data.obj <- data.obj2
    dist.obj <- dist.obj2
    
    data.obj$meta.dat$icd10_first_3_name <- as.character(data.obj$meta.dat$icd10_first_3_name)
    samIDs <- rownames(data.obj$meta.dat[data.obj$meta.dat$icd10_first_3_name %in% names(which(sort(table(data.obj$meta.dat$icd10_first_3_name))>15)) & data.obj$meta.dat$Group != 'Control',])
    data.obj <- subset_data(data.obj, samIDs = samIDs)
    dist.obj <- subset_dist(dist.obj, samIDs = samIDs)
    data.obj$meta.dat$icd10_first_3_name <- as.factor(data.obj$meta.dat$icd10_first_3_name)
    
    samIDs <- rownames(data.obj$meta.dat[data.obj$meta.dat$icd10_first_3_name %in% cancer.type & data.obj$meta.dat$Group != 'Control',])
    data.obj <- subset_data(data.obj, samIDs = samIDs)
    dist.obj <- subset_dist(dist.obj, samIDs = samIDs)
    
    
    idx <- !(data.obj$meta.dat$icd10_first_3_name %in% c(cancer))
    data.obj$meta.dat$icd10_first_3_name <- as.character(data.obj$meta.dat$icd10_first_3_name)
    data.obj$meta.dat$icd10_first_3_name[idx] <- 'Others'
    dist.obj <- subset_dist(dist.obj, samIDs = rownames(data.obj$meta.dat))
    
    data.obj$meta.dat$icd10_first_3_name <- factor(data.obj$meta.dat$icd10_first_3_name,levels=c('Others',cancer))
    
    save(data.obj, dist.obj, file = 'data.obj.wk.RData')
    gc()
  }
}


for(func.type in c('pathway')){#'uniref90_ko','uniref90_go','uniref90_level4ec',
  setwd(wd)
  load(file = paste0('Data/data.obj.',func.type,'.RData')) 
  
  data.obj2 <- data.obj
  dist.obj2 <- dist.obj
  
  data.obj$meta.dat$icd10_first_3_name <- as.character(data.obj$meta.dat$icd10_first_3_name)
  samIDs <- rownames(data.obj$meta.dat[data.obj$meta.dat$icd10_first_3_name %in% names(which(sort(table(data.obj$meta.dat$icd10_first_3_name))>15)) & data.obj$meta.dat$Group != 'Control',])
  data.obj <- subset_data(data.obj, samIDs = samIDs)
  dist.obj <- subset_dist(dist.obj, samIDs = samIDs)
  
  data.obj$meta.dat$icd10_first_3_name <- as.factor(data.obj$meta.dat$icd10_first_3_name)
  
  cancer.type <- as.vector(unique(data.obj$meta.dat$icd10_first_3_name)[!is.na(unique(data.obj$meta.dat$icd10_first_3_name))])
  
  cancer.dir <- gsub('Malignant neoplasm of |Multiple myeloma and malignant |Malignant ','',cancer.type)
  cancer.dir <- gsub('Other\\ and\\ |malignant\\ ','',cancer.dir)## add 2024/10/31 for adding 2 other cancer back
  cancer.dir <- gsub('__','_',gsub(',|\\ ','_',cancer.dir))
  cancer.dir <- gsub('eye_brain_and_other_parts_of_|_and','',cancer.dir)
  getwd()
  
  
  variables <- c("icd10_first_3_name")
  for(dir in paste0('subCancerX-Ex2_func/',func.type,'/',cancer.dir)){
    for(variable in variables){
      
      sh <- paste(
        "sbatch",
        paste0("-J ", variable),
        "--partition=cpu-short",
        paste("--output=/research/bsi/projects/staff_analysis/m216453/2023_01_09_Oncobiome/trash/OUT/",  variable, "_",gsub('\\/','_',dir),".out", sep=""),
        paste("--error=/research/bsi/projects/staff_analysis/m216453/2023_01_09_Oncobiome/trash/OUT/",  variable,"_",gsub('\\/','_',dir),".err", sep=""),
        "--time=1-00:00:00",
        "--mem=64G",
        "--export=ALL",
        paste("\'--wrap=R CMD BATCH --no-restore \"--args ", variable," ",dir,"\" ",
              "/research/bsi/projects/staff_analysis/m216453/2023_01_09_Oncobiome/Code/Alpha_Beta_DAA_func_clean.R", " ", paste0("/research/bsi/projects/staff_analysis/m216453/2023_01_09_Oncobiome/trash/",variable,"_",gsub('\\/','_',dir),".Rout"), "\'", sep="")
      )
      print(sh)
      system(sh)
    }
  }
  
}





