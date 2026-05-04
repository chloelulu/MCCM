############################# Define script ######################################
file_dir <- dirname(rstudioapi::getSourceEditorContext()$path)
file_dir <- dirname(dirname(dirname(dirname(dirname(file_dir)))))
taxon_script <- paste0(file_dir,'/Code/Submission/MayoOncobiomeStudy/Code/metaphlan/Alpha_Beta_DAA_clean.R')
func_script <- paste0(file_dir,'/Code/Submission/MayoOncobiomeStudy/Code/metaphlan/Alpha_Beta_DAA_func_clean.R')

# mph_choice <- "mph_commonKraken"
mph_choice <- "mph"
############################## Cancer Only #######################################
file_dir <- dirname(rstudioapi::getSourceEditorContext()$path)
file_dir <- dirname(dirname(dirname(dirname(dirname(file_dir)))))

wd <- file_dir
rd <- paste0(file_dir,'/Result_',mph_choice,'/')

setwd(wd)
source(paste0(file_dir,"/Code/Submission/MayoOncobiomeStudy/Code/Stats.R"))
try(load_package())

getwd()
load(file = 'Data/data.obj.mph.RData') 


setwd(rd)
dir <- 'CancerOnly'
if(!dir.exists(dir)){dir.create(dir)}
setwd(dir)
getwd()

ind <- data.obj$meta.dat$Group == 'Cancer'
data.obj <- subset_data(data.obj, ind)
dist.obj <- subset_dist(dist.obj, ind)

save(data.obj, dist.obj, file = 'data.obj.wk.RData')


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
for(variable in variables[1:length(variables)]){
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
file_dir <- dirname(dirname(dirname(dirname(dirname(file_dir)))))

wd <- file_dir
rd <- paste0(file_dir,'/Result_',mph_choice,'/')

setwd(wd)
source(paste0(file_dir,"/Code/Submission/MayoOncobiomeStudy/Code/Stats.R"))
try(load_package())

getwd()
load(file = 'Data/data.obj.mph.RData') 


setwd(rd)
dir <- 'LungCancer_Smoking'
if(!dir.exists(dir)){dir.create(dir)}
setwd(dir)
getwd()

ind <- data.obj$meta.dat$icd10_first_3_name_short == 'bronchus lung'
data.obj <- subset_data(data.obj, ind)
dist.obj <- subset_dist(dist.obj, ind)

data.obj$meta.dat$smoking_category2 <- data.obj$meta.dat$smoking_category
data.obj$meta.dat$smoking_category2[data.obj$meta.dat$smoking_category2 %in% c("Current","Former")] <- "Yes"

save(data.obj, dist.obj, file = 'data.obj.wk.RData')

variables <- c("smoking_category","smoking_category2")
dir <- 'LungCancer_Smoking'
for(variable in variables){
  with(list(commandArgs = function(...) c("--args", variable, dir)), source(taxon_script))
}



############################## Control vs sub Cancer X #######################################
file_dir <- dirname(rstudioapi::getSourceEditorContext()$path)
file_dir <- dirname(dirname(dirname(dirname(dirname(file_dir)))))

wd <- file_dir
rd <- paste0(file_dir,'/Result_',mph_choice,'/subCancerX_Control/')

if(!dir.exists(rd)) dir.create(rd)
setwd(wd)
source(paste0(file_dir,"/Code/Submission/MayoOncobiomeStudy/Code/Stats.R"))
try(load_package())

load(file = 'Data/data.obj.mph.RData') 

data.obj2 <- data.obj
dist.obj2 <- dist.obj

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
  
  data.obj$meta.dat$icd10_first_3_name_short <- factor(data.obj$meta.dat$icd10_first_3_name_short,levels=c('healthy',cancer))

  save(data.obj, dist.obj, file = 'data.obj.wk.RData')
  setwd('..')
}

# setwd("/Users/luyang1/myicloud/Documents/Mayo_project/2023_01_09_Oncobiome/mforge_clean/Result_mph/subCancerX_Control/")
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
file_dir <- dirname(dirname(dirname(dirname(dirname(file_dir)))))

wd <- file_dir
rd <- paste0(file_dir,'/Result_',mph_choice,'/subCancerX-Ex/')

if(!dir.exists(rd)) dir.create(rd)
setwd(wd)
source(paste0(file_dir,"/Code/Submission/MayoOncobiomeStudy/Code/Stats.R"))
try(load_package())

load(file = 'Data/data.obj.mph.RData') 

data.obj2 <- data.obj
dist.obj2 <- dist.obj

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
  dist.obj <- dist.obj2

  data.obj$meta.dat$icd10_first_3_name_short <- as.character(data.obj$meta.dat$icd10_first_3_name_short)
  samIDs <- rownames(data.obj$meta.dat[data.obj$meta.dat$icd10_first_3_name_short %in% cancer.type,])
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

variable <- c("icd10_first_3_name_short")
for(dir in paste0('subCancerX-Ex/',cancer.type)){
  with(list(commandArgs = function(...) c("--args", variable, dir)), source(taxon_script))
}


# setwd("/Users/luyang1/myicloud/Documents/Mayo_project/2023_01_09_Oncobiome/mforge_clean/Result_mph/subCancerX-Ex/")
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
file_dir <- dirname(dirname(dirname(dirname(dirname(file_dir)))))

wd <- file_dir
rd <- paste0(file_dir,'/Result_',mph_choice,'/')

if(!dir.exists(rd)) dir.create(rd)
setwd(wd)
source(paste0(file_dir,"/Code/Submission/MayoOncobiomeStudy/Code/Stats.R"))
try(load_package())

load(file = 'Data/data.obj.mph.RData') 

setwd(rd)
dir <- 'PanCancer'
if(!dir.exists(dir)){dir.create(dir)}
setwd(dir)
getwd()

data.obj$meta.dat$Cancer_class <- as.character(data.obj$meta.dat$Cancer_class)
data.obj$meta.dat$Cancer_class[data.obj$meta.dat$Group=='Control'] <- 'Control'

data.obj$meta.dat$Group <- factor(data.obj$meta.dat$Group,levels=c('Control','Cancer'))

save(data.obj, dist.obj, file = 'data.obj.wk.RData')

variable <- "Group"; dir <- 'PanCancer'
with(list(commandArgs = function(...) c("--args", variable, dir)), source(taxon_script))


############################## early onset #######################################
file_dir <- dirname(rstudioapi::getSourceEditorContext()$path)
file_dir <- dirname(dirname(dirname(dirname(dirname(file_dir)))))

wd <- file_dir
rd <- paste0(file_dir,'/Result_',mph_choice,'/')

if(!dir.exists(rd)) dir.create(rd)
setwd(wd)
source(paste0(file_dir,"/Code/Submission/MayoOncobiomeStudy/Code/Stats.R"))
try(load_package())

load(file = 'Data/data.obj.mph.RData') 

ind <- data.obj$meta.dat$Group == 'Cancer'
data.obj2 <- subset_data(data.obj, ind)
dist.obj2 <- subset_dist(dist.obj, ind)


setwd(rd)
dir <- 'EarlyOnset'
if(!dir.exists(dir)){dir.create(dir)}
setwd(dir)
getwd()

tb <- cbind(table(data.obj2$meta.dat$icd10_first_3_name_short, data.obj2$meta.dat$early_onset))
idx <- names(which(rowSums(tb >15) ==2))
for(j in idx){
  setwd(paste0(rd,'EarlyOnset'))
  if(!dir.exists(j)){dir.create(j)}
  setwd(j)
  getwd()
  
  data.obj <- data.obj2
  dist.obj <- dist.obj2

  ind2 <- data.obj$meta.dat$icd10_first_3_name_short == j
  data.obj <- subset_data(data.obj, ind2)
  dist.obj <- subset_dist(dist.obj, ind2)
  table(data.obj$meta.dat$early_onset)
  
  save(data.obj, dist.obj, file = 'data.obj.wk.RData')

}


variable <- 'early_onset'
for(dir in paste0('EarlyOnset/',idx)){
  with(list(commandArgs = function(...) c("--args", variable, dir)), source(taxon_script))
}



