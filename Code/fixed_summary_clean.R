## Date: 12/18/2023, rerun all code for clean up the useful part, this script is used for fixed summary of each part. For detailed analysis and figure generation, refer to FinalFigure_clean.R
pkg <- c('tibble','dplyr','tidyr','ggtree','tidytree','treeio','ggmosaic','cluster','ComplexHeatmap','circlize','ggpubr')
sapply(pkg, function(x)require(x, character.only=T))
# source('/research/bsi/projects/staff_analysis/m216453/2023_01_09_Oncobiome/Code/Stats.R')
# try(load_package())
# fd <- "/research/bsi/projects/staff_analysis/m216453/2023_01_09_Oncobiome/Figure/"
# wd <- '/research/bsi/projects/staff_analysis/m216453/2023_01_09_Oncobiome/'
# rd <- '/research/bsi/projects/staff_analysis/m216453/2023_01_09_Oncobiome/Result/'
# source("/Users/M216453/Documents/Mayo_project/2023_01_09_Oncobiome/mforge_clean/Code/Stats.R")
source("/Users/luyang1/myicloud/Documents/Mayo_project/2023_01_09_Oncobiome/mforge_clean/Code/Stats.R")
try(load_package())

# fd <- "/Users/M216453/Documents/Mayo_project/2023_01_09_Oncobiome/mforge_clean/Figure/"
# wd <- '/Users/M216453/Documents/Mayo_project/2023_01_09_Oncobiome/mforge_clean/'
# rd <- '/Users/M216453/Documents/Mayo_project/2023_01_09_Oncobiome/mforge_clean/Result/'
fd <- "/Users/luyang1/myicloud/Documents/Mayo_project/2023_01_09_Oncobiome/mforge_clean/Figure/"
wd <- '/Users/luyang1/myicloud/Documents/Mayo_project/2023_01_09_Oncobiome/mforge_clean/'
rd <- '/Users/luyang1/myicloud/Documents/Mayo_project/2023_01_09_Oncobiome/mforge_clean/Result/'

setwd(wd)
source('Code/summary_plot_func.R')

############################## Cancer Only #######################################
q.cut <- 0.05; p.cut <- 0.05; abund.cut <- 0.01; prev.cut <- 0.2
q.cut1 <- 0.001; q.cut2 <- 0.01; q.cut3 <- 0.05
covars <- c("Batch","Bristol_score","BMI", "Age", "Sex", "GI_nonGI","Cancer_class","Metastasis","PPI_day_365", "Abx_day_365", 
            "PPI_last_month","Abx_last_month","Charlson_score","Elix_score","Sample_season","Urban" ,"icd10_first_3_name_short", "Site")
blood.names <- c('Erythrocytes','Hematocrit','Neutrophils','MCV','Hemoglobin','Leukocytes','Platelet.Count',
                 'Neutrophils_cat','Platelet.Count_cat','Hemoglobin_cat',
                 "neut_neutropenia2","bone_marrow_suppression2","Hb_anemia2","Pl_thrombocytopenia2","neut_neutropenia_c","Hb_anemia_c","Pl_thrombocytopenia_c")
elix.names <- c("Elixhauser_CHF","Elixhauser_Arrhythmia","Elixhauser_Valvular","Elixhauser_PHTN","Elixhauser_PVD","Elixhauser_HTN","Elixhauser_Paralysis",  
                "Elixhauser_NeuroOther","Elixhauser_Pulmonary","Elixhauser_DM","Elixhauser_DMcx","Elixhauser_Hypothyroid","Elixhauser_Renal","Elixhauser_Liver",      
                "Elixhauser_PUD" ,"Elixhauser_HIV","Elixhauser_Lymphoma","Elixhauser_Mets","Elixhauser_Tumor","Elixhauser_Rheumatic","Elixhauser_Coagulopathy",
                "Elixhauser_Obesity","Elixhauser_WeightLoss","Elixhauser_FluidsLytes","Elixhauser_BloodLoss","Elixhauser_Anemia","Elixhauser_Alcohol","Elixhauser_Drugs",
                "Elixhauser_Psychoses","Elixhauser_Depression")
dirs <- c(covars[!(covars %in% c("Batch","icd10_first_3_name"))], blood.names, elix.names)
main.dirs <- 'CancerOnly'
obj <- summary_plot(
  main.dirs = main.dirs,
  dirs = dirs,
  figure.dir = main.dirs,
  alpha.summary = T,
  alpha.plot = F,
  beta.summary = F,
  beta.plot =F,

  taxa.levels = c("Species"),
  DAA.summary =F,
  DAA.plot = F,

  q.cut = q.cut,
  p.cut = p.cut,
  abund.cut = abund.cut,
  prev.cut = prev.cut,

  volcano.plot =F,
  volcano.top = 10,
  volcano.level = 'Genus',

  bar.p =F,
  taxa.level.bar = 'Species',

  heat.p = F,
  taxa.level.heatmap = 'Species',
  q.cut1 = q.cut1, q.cut2 = q.cut2, q.cut3 = q.cut3,
  plot.method = "all.top",
  top.n = 5,
  heatmap.width = 7,
  heatmap.height = 15,
  fd = "/Users/luyang1/myicloud/Documents/Mayo_project/2023_01_09_Oncobiome/mforge_clean/Figure/",
  wd = '/Users/luyang1/myicloud/Documents/Mayo_project/2023_01_09_Oncobiome/mforge_clean/',
  rd = '/Users/luyang1/myicloud/Documents/Mayo_project/2023_01_09_Oncobiome/mforge_clean/Result/'
)


if(length(main.dirs)>1){name= main.dirs}else{name = dirs}
ct <- tot <- c()
for(main.dir in main.dirs){
  for(dir in dirs){
    setwd(paste0(rd,'/',main.dir,'/DAA'))
    setwd(dir)
    file <- list.files(pattern = 'Rdata$')
    diff.obj <- NULL
    load(file)
    tot <- c(tot,nrow(diff.obj$qv.list$Species))
    ct <- c(ct,sum(diff.obj$qv.list$Species <q.cut, na.rm = T))
  }
}
names(ct) <- names(tot) <- name
sigs.CancerOnly <- cbind.data.frame(sig.n = ct, test.n =tot) %>%
  mutate(sig.pct = sig.n/test.n)


setwd(fd)
setwd('CancerOnly/')
load('DAA_P_R2.RData')
dirs <- c(covars[!(covars %in% c("Batch","icd10_first_3_name"))], blood.names, elix.names)
for(level in c("Species")){
  df.list <- list()
  for(dir in dirs){
    cat(dir,'\n')
    if(dir %in% colnames(P.All[[level]])){
      df.list[[dir]] <- inner_join(P.All[[level]][,c(level,dir)] %>% dplyr::rename(P = !!as.name(dir)),
                                   R2.All[[level]][,c(level,dir)] %>% dplyr::rename(R2 = !!as.name(dir))) %>%
        inner_join(Q.All[[level]][,c(level,dir)] %>% dplyr::rename(Q = !!as.name(dir))) %>% mutate(variable = dir)

    }

  }
  df <- do.call(rbind.data.frame, df.list)
  rownames(df) <- NULL
  setwd(fd)
  write.csv(df, file = paste0('CancerOnly/',level,'_DAA.csv'), row.names = F)
}




############################## subCancerX vs Control #######################################
setwd(rd)
setwd('subCancerX_Control')
variables <- unique(gsub('/.*','',list.dirs(full.names = F)))
variables <- variables[grep('^Control-',variables)]
main.dirs <-  variables
dirs <- 'icd10_first_3_name_short'
figure.dir <- 'subCancerX_Control'

q.cut <- 0.05; p.cut <- 0.05; abund.cut <- 0.01; prev.cut <- 0.2
q.cut1 <- 0.001; q.cut2 <- 0.01; q.cut3 <- q.cut # for * visualized in the plot
rd.tmp <- "/Users/luyang1/myicloud/Documents/Mayo_project/2023_01_09_Oncobiome/mforge_clean/Result/subCancerX_Control/"
if(!dir.exists(rd.tmp)) dir.create(rd.tmp)

obj <- summary_plot(
  main.dirs = main.dirs,
  dirs =  dirs,
  figure.dir = figure.dir,
  alpha.summary = T,
  alpha.plot = F,
  beta.summary =F,
  beta.plot =F,
  DAA.summary = F,
  taxa.levels = c("Species"),
  volcano.plot = F,
  volcano.top = 10,
  volcano.level ='Species',
  DAA.plot = F,
  heat.p =F,
  q.cut = q.cut,
  p.cut = p.cut,
  abund.cut = abund.cut,
  prev.cut = prev.cut,
  q.cut1 = q.cut1, q.cut2 = q.cut2, q.cut3 = q.cut3,
  plot.method = "each.top.R2",
  taxa.level.heatmap = 'Species',
  top.n = 10,
  heatmap.width = 7,
  heatmap.height = 15,
  rd = rd.tmp,
  fd = "/Users/luyang1/myicloud/Documents/Mayo_project/2023_01_09_Oncobiome/mforge_clean/Figure/",
  wd = '/Users/luyang1/myicloud/Documents/Mayo_project/2023_01_09_Oncobiome/'
)

# setwd(mfd)
# pdf('subCancerX_Control_heatmap.pdf', width = 15, height = 15)
# obj
# dev.off()


setwd(fd)
setwd('subCancerX_Control/')
load('DAA_P_R2.RData')
dirs <- main.dirs
for(level in c("Species")){
  df.list <- list()
  for(dir in dirs){
    df.list[[dir]] <- inner_join(P.All[[level]][,c(level,dir)] %>% dplyr::rename(P = !!as.name(dir)),
                                 R2.All[[level]][,c(level,dir)] %>% dplyr::rename(R2 = !!as.name(dir))) %>%
      inner_join(Q.All[[level]][,c(level,dir)] %>% dplyr::rename(Q = !!as.name(dir))) %>% mutate(variable = dir)

  }
  df <- do.call(rbind.data.frame, df.list)
  rownames(df) <- NULL
  setwd(fd)
  write.csv(df, file = paste0('subCancerX_Control/',level,'_DAA.csv'), row.names = F)
}



############################## subCancerX vs [Cancers-subCancerX] #######################################
setwd(rd)
setwd('subCancerX-Ex')
variables <- unique(gsub('/.*','',list.dirs(full.names = F)))
main.dirs <-  variables[-grep('^$',variables)]
dirs <- 'icd10_first_3_name_short'
figure.dir <- 'subCancerX-Ex'

q.cut <- 0.2; p.cut <- 0.05; abund.cut <- 0.01; prev.cut <- 0.2
q.cut1 <- 0.001; q.cut2 <- 0.01; q.cut3 <- q.cut # for * visualized in the plot
rd.tmp <- '/Users/luyang1/myicloud/Documents/Mayo_project/2023_01_09_Oncobiome/mforge_clean/Result/subCancerX-Ex/'
if(!dir.exists(rd.tmp)) dir.create(rd.tmp)

obj <- summary_plot(
  main.dirs = main.dirs,
  dirs =  dirs,
  figure.dir = figure.dir,
  alpha.summary =T,
  alpha.plot = F,
  beta.summary =F,
  beta.plot =F,
  DAA.summary =F,
  taxa.levels = c("Species"),
  volcano.plot =F,
  volcano.level = 'Species',
  volcano.top = 10,
  DAA.plot = F,
  heat.p = F,
  q.cut = q.cut,
  p.cut = p.cut,
  abund.cut = abund.cut,
  prev.cut = prev.cut,
  q.cut1 = q.cut1, q.cut2 = q.cut2, q.cut3 = q.cut3,
  plot.method = "each.top.R2",
  taxa.level.heatmap = 'Genus',
  top.n = 10,
  heatmap.width = 7,
  heatmap.height = 15,
  fd = "/Users/luyang1/myicloud/Documents/Mayo_project/2023_01_09_Oncobiome/mforge_clean/Figure/",
  wd = '/Users/luyang1/myicloud/Documents/Mayo_project/2023_01_09_Oncobiome/mforge_clean/',
  rd = rd.tmp
  
)




setwd(fd)
setwd('subCancerX-Ex/')
load('DAA_P_R2.RData')
dirs <- main.dirs
for(level in c("Species")){
  df.list <- list()
  for(dir in dirs){
    df.list[[dir]] <- inner_join(P.All[[level]][,c(level,dir)] %>% dplyr::rename(P = !!as.name(dir)),
                                 R2.All[[level]][,c(level,dir)] %>% dplyr::rename(R2 = !!as.name(dir))) %>%
      inner_join(Q.All[[level]][,c(level,dir)] %>% dplyr::rename(Q = !!as.name(dir))) %>% mutate(variable = dir)

  }
  df <- do.call(rbind.data.frame, df.list)
  rownames(df) <- NULL
  setwd(fd)
  write.csv(df, file = paste0('subCancerX-Ex/',level,'_DAA.csv'), row.names = F)
}


############################## CancerOnly func #######################################
setwd(rd)
q.cut <- 0.05; p.cut <- 0.05; abund.cut <- 0.01; prev.cut <- 0.2
q.cut1 <- 0.001; q.cut2 <- 0.01; q.cut3 <- 0.05
covars <- c("Batch","Bristol_score","BMI", "Age", "Sex", "GI_nonGI","Cancer_class","Metastasis","PPI_day_365", "Abx_day_365",
            "PPI_last_month","Abx_last_month","Charlson_score","Elix_score","Sample_season","Urban","icd10_first_3_name_short", "Site")
blood.names <- c('Erythrocytes','Hematocrit','Neutrophils','MCV','Hemoglobin','Leukocytes','Platelet.Count',
                 'Neutrophils_cat','Platelet.Count_cat','Hemoglobin_cat',
                 "neut_neutropenia2","bone_marrow_suppression2","Hb_anemia2","Pl_thrombocytopenia2","neut_neutropenia_c","Hb_anemia_c","Pl_thrombocytopenia_c")
elix.names <- c("Elixhauser_CHF","Elixhauser_Arrhythmia","Elixhauser_Valvular","Elixhauser_PHTN","Elixhauser_PVD","Elixhauser_HTN","Elixhauser_Paralysis",  
                "Elixhauser_NeuroOther","Elixhauser_Pulmonary","Elixhauser_DM","Elixhauser_DMcx","Elixhauser_Hypothyroid","Elixhauser_Renal","Elixhauser_Liver",      
                "Elixhauser_PUD" ,"Elixhauser_HIV","Elixhauser_Lymphoma","Elixhauser_Mets","Elixhauser_Tumor","Elixhauser_Rheumatic","Elixhauser_Coagulopathy",
                "Elixhauser_Obesity","Elixhauser_WeightLoss","Elixhauser_FluidsLytes","Elixhauser_BloodLoss","Elixhauser_Anemia","Elixhauser_Alcohol","Elixhauser_Drugs",
                "Elixhauser_Psychoses","Elixhauser_Depression")

variables <- c(covars[!(covars %in% c("Batch","icd10_first_3_name"))], blood.names, elix.names)
dirs <- c(covars[!(covars %in% c("Batch","icd10_first_3_name"))], blood.names, elix.names)

main.dirs <- 'CancerOnly_func/pathway'

obj <- summary_plot(
  main.dirs = main.dirs,
  dirs =  dirs,
  figure.dir = 'CancerOnly_func/pathway',
  alpha.summary = F,
  alpha.plot = F,
  beta.summary =F,
  beta.plot = F,
  DAA.summary =T,
  volcano.plot =F,
  volcano.level = 'pathway',
  taxa.levels = 'pathway',
  volcano.top = 10,
  DAA.plot =F,
  bar.p = F,
  taxa.level.bar = 'pathway',
  heat.p = F,
  q.cut = q.cut,
  p.cut = p.cut,
  abund.cut = abund.cut,
  prev.cut = prev.cut,
  q.cut1 = q.cut1, q.cut2 = q.cut2, q.cut3 = q.cut3,
  plot.method = "each.top",
  taxa.level.heatmap = 'pathway',
  top.n = 20,
  heatmap.width = 7,
  heatmap.height = 15,
  fd = "/Users/luyang1/myicloud/Documents/Mayo_project/2023_01_09_Oncobiome/mforge_clean/Figure/",
  wd = '/Users/luyang1/myicloud/Documents/Mayo_project/2023_01_09_Oncobiome/mforge_clean/',
  rd = '/Users/luyang1/myicloud/Documents/Mayo_project/2023_01_09_Oncobiome/mforge_clean/Result/'
  
)


setwd(fd)
setwd('CancerOnly_func/')
load('pathway/DAA_P_R2.RData')
dirs <- c(covars[!(covars %in% c("Batch"))], blood.names, elix.names)
for(level in c("pathway")){
  df.list <- list()
  for(dir in dirs){
    cat(dir,'\n')
    if(dir %in% colnames(P.All[[level]])){
      df.list[[dir]] <- inner_join(P.All[[level]][,c(level,dir)] %>% dplyr::rename(P = !!as.name(dir)),
                                   R2.All[[level]][,c(level,dir)] %>% dplyr::rename(R2 = !!as.name(dir))) %>%
        inner_join(Q.All[[level]][,c(level,dir)] %>% dplyr::rename(Q = !!as.name(dir))) %>% mutate(variable = dir)

    }

  }
  df <- do.call(rbind.data.frame, df.list)
  rownames(df) <- NULL
  setwd(fd)
  write.csv(df, file = paste0('CancerOnly_func/pathway/',level,'_DAA.csv'), row.names = F)
}


############################## subCancerX vs Control func #######################################
setwd(rd)
setwd('subCancerX_Control_func/pathway')
variables <- unique(gsub('/.*','',list.dirs(full.names = F)))
main.dirs <-  variables[-grep('^$',variables)]
dirs <- 'icd10_first_3_name_short'
q.cut <- 0.05; p.cut <- 0.05; abund.cut <- 0.002; prev.cut <- 0.2
q.cut1 <- 0.001; q.cut2 <- 0.01; q.cut3 <- q.cut # for * visualized in the plot
rd.tmp <- "/Users/luyang1/myicloud/Documents/Mayo_project/2023_01_09_Oncobiome/mforge_clean/Result/subCancerX_Control_func/pathway/"
if(!dir.exists(rd.tmp)) dir.create(rd.tmp)
obj <- summary_plot(
  main.dirs = main.dirs,
  dirs =  dirs,
  figure.dir = 'subCancerX_Control_func/pathway',
  alpha.summary = F,
  alpha.plot = F,
  beta.summary = F,
  beta.plot = F,
  DAA.summary = T,
  volcano.plot = F,
  volcano.level = 'pathway',
  taxa.levels = 'pathway',
  volcano.top = 10,
  DAA.plot = F,
  heat.p = F,
  q.cut = q.cut,
  p.cut = p.cut,
  abund.cut = abund.cut,
  prev.cut = prev.cut,
  q.cut1 = q.cut1, q.cut2 = q.cut2, q.cut3 = q.cut3,
  plot.method = "each.top",
  taxa.level.heatmap = 'pathway',
  top.n = 20,
  heatmap.width = 7,
  heatmap.height = 15,
  wd = '/Users/luyang1/myicloud/Documents/Mayo_project/2023_01_09_Oncobiome/mforge_clean/',
  fd = "/Users/luyang1/myicloud/Documents/Mayo_project/2023_01_09_Oncobiome/mforge_clean/Figure/",
  rd = rd.tmp
)



setwd(fd)
dirs <- main.dirs
for(level in c('pathway')){
  setwd(paste0('subCancerX_Control_func/',level))
  load('DAA_P_R2.RData')
  df.list <- list()
  for(dir in dirs){
    cat(dir,'\n')
    df.list[[dir]] <- inner_join(P.All[[level]][,c(level,dir)] %>% dplyr::rename(P = !!as.name(dir)),
                                 R2.All[[level]][,c(level,dir)] %>% dplyr::rename(R2 = !!as.name(dir))) %>%
      inner_join(Q.All[[level]][,c(level,dir)] %>% dplyr::rename(Q = !!as.name(dir))) %>% mutate(variable = dir)

  }
  df <- do.call(rbind.data.frame, df.list)
  rownames(df) <- NULL
  setwd(fd)
  write.csv(df, file = paste0('subCancerX_Control_func/',level,'_DAA.csv'), row.names = F)
}


############################## subCancerX vs [Cancers-subCancerX] func #######################################
setwd(rd)
setwd('subCancerX-Ex_func/pathway')
variables <- unique(gsub('/.*','',list.dirs(full.names = F)))
main.dirs <-  variables[-grep('^$',variables)]
dirs <- 'icd10_first_3_name_short'
q.cut <- 0.05; p.cut <- 0.05; abund.cut <- 0.002; prev.cut <- 0.2
q.cut1 <- 0.001; q.cut2 <- 0.01; q.cut3 <- q.cut # for * visualized in the plot
rd.tmp <- '/Users/luyang1/myicloud/Documents/Mayo_project/2023_01_09_Oncobiome/mforge_clean/Result/subCancerX-Ex_func/pathway/'
if(!dir.exists(rd.tmp)) dir.create(rd.tmp)
obj <- summary_plot(
  main.dirs = main.dirs,
  dirs =  dirs,
  figure.dir = 'subCancerX-Ex_func/pathway',
  alpha.summary = F,
  alpha.plot = F,
  beta.summary =F,
  beta.plot = F,
  DAA.summary =T,
  volcano.plot =F,
  volcano.level = 'pathway',
  taxa.levels = 'pathway',
  volcano.top = 10,
  DAA.plot = F,
  heat.p =F,
  q.cut = q.cut,
  p.cut = p.cut,
  abund.cut = abund.cut,
  prev.cut = prev.cut,
  q.cut1 = q.cut1, q.cut2 = q.cut2, q.cut3 = q.cut3,
  plot.method = "each.top.R2",
  taxa.level.heatmap = 'pathway',
  top.n = 10,
  heatmap.width = 7,
  heatmap.height = 15,
  fd = "/Users/luyang1/myicloud/Documents/Mayo_project/2023_01_09_Oncobiome/mforge_clean/Figure/",
  wd = '/Users/luyang1/myicloud/Documents/Mayo_project/2023_01_09_Oncobiome/mforge_clean/',
  rd = rd.tmp
)



dirs <- main.dirs
for(level in c('pathway')){
  setwd(fd)
  setwd(paste0('subCancerX-Ex_func/',level))
  load('DAA_P_R2.RData')
  df.list <- list()
  for(dir in dirs){
    df.list[[dir]] <- inner_join(P.All[[level]][,c(level,dir)] %>% dplyr::rename(P = !!as.name(dir)),
                                 R2.All[[level]][,c(level,dir)] %>% dplyr::rename(R2 = !!as.name(dir))) %>%
      inner_join(Q.All[[level]][,c(level,dir)] %>% dplyr::rename(Q = !!as.name(dir))) %>% mutate(variable = dir)

  }
  df <- do.call(rbind.data.frame, df.list)
  rownames(df) <- NULL
  setwd(fd)
  write.csv(df, file = paste0('subCancerX-Ex_func/',level,'_DAA.csv'), row.names = F)
}



