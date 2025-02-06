args=(commandArgs(TRUE))
if(length(args)==0){
  print("No arguments supplied.")
} else{
  variable = args[1]
  dir = args[2]
}



source('/research/bsi/projects/staff_analysis/m216453/2023_01_09_Oncobiome/Code/Stats.R')
try(load_package())

rd <- '/research/bsi/projects/staff_analysis/m216453/2023_01_09_Oncobiome/Result/'
setwd(paste0(rd, dir))
load(file = 'data.obj.wk.RData') 

#-------- set grp.name(variable) & adj.name --------
cat(variable, '\n')
  adj.name <- NULL
  
  covars <- c("Batch","Bristol_score","BMI", "Age", "Sex", "GI_nonGI","Cancer_class","Metastasis","PPI_day_365", "Abx_day_365",
              "Abx_last_month","Charlson_score","Elix_score","Sample_season","Urban")
  blood.names <- c('Erythrocytes','Hematocrit','Neutrophils','MCV','Hemoglobin','Leukocytes','Platelet.Count',
                   'Neutrophils_cat','Platelet.Count_cat','Hemoglobin_cat',
                   "neut_neutropenia2","bone_marrow_suppression2","Hb_anemia2","Pl_thrombocytopenia2","neut_neutropenia_c","Hb_anemia_c","Pl_thrombocytopenia_c")
  elix.names <- c("Elixhauser_CHF","Elixhauser_Arrhythmia","Elixhauser_Valvular","Elixhauser_PHTN","Elixhauser_PVD","Elixhauser_HTN","Elixhauser_Paralysis",  
                  "Elixhauser_NeuroOther","Elixhauser_Pulmonary","Elixhauser_DM","Elixhauser_DMcx","Elixhauser_Hypothyroid","Elixhauser_Renal","Elixhauser_Liver",      
                  "Elixhauser_PUD" ,"Elixhauser_HIV","Elixhauser_Lymphoma","Elixhauser_Mets","Elixhauser_Tumor","Elixhauser_Rheumatic","Elixhauser_Coagulopathy",
                  "Elixhauser_Obesity","Elixhauser_WeightLoss","Elixhauser_FluidsLytes","Elixhauser_BloodLoss","Elixhauser_Anemia","Elixhauser_Alcohol","Elixhauser_Drugs",
                  "Elixhauser_Psychoses","Elixhauser_Depression")
  
  if(grepl('CancerOnly',dir)){
    if(variable =='Cancer_class'){
      adj.name <- covars[!(covars %in% c(variable,'GI_nonGI',"Charlson_score","Abx_last_month"))]
    }else if(variable =='GI_nonGI'){
      adj.name <- covars[!(covars %in% c(variable,'Cancer_class',"Charlson_score","Abx_last_month"))]
    }else if(variable =="Charlson_score"){
      adj.name <- covars[!(covars %in% c(variable,"Elix_score",'Cancer_class',"Abx_last_month"))]
    }else if(variable =="Elix_score"){
      adj.name <- covars[!(covars %in% c(variable,"Charlson_score",'Cancer_class',"Abx_last_month"))]
    }else if(variable =="Abx_last_month"){
      adj.name <- covars[!(covars %in% c(variable,'Cancer_class',"Charlson_score","Abx_day_365","Abx_last_month"))]
    }else if(variable %in% c('Elixhauser_Mets','Charlson_Mets')){# exclude Metastasis
      adj.name <- covars[!(covars %in% c(variable,'GI_nonGI',"Charlson_score", "Elix_score","Metastasis","Abx_last_month"))]
    }else if(variable %in% c('Elixhauser_Obesity','Elixhauser_WeightLoss')){#exclude BMI
      adj.name <- covars[!(covars %in% c(variable,'GI_nonGI',"Charlson_score", "Elix_score","BMI","Abx_last_month"))]
    }else if(length(grep('^Elixhauser|^Charlson',variable))==1){
      adj.name <- covars[!(covars %in% c(variable,'GI_nonGI',"Charlson_score","Abx_last_month","Elix_score"))]
    }else if(length(grep('icd10_first_3_name',variable))==1){
      adj.name <- covars[!(covars %in% c(variable,'GI_nonGI',"Charlson_score","Abx_last_month","PPI_last_month","Cancer_class"))]
    }else{
      adj.name <- covars[!(covars %in% c(variable,'GI_nonGI',"Charlson_score","Abx_last_month"))]
    }
  }
  
  # currently PanCancer Only contains Batch 5, We aim at checking how is the difference if we control batch. And we drop previous Pancancer
  if(length(grep('PanCancer_func',dir))==1){ adj.name <- c('BMI','Sex') }
  
  if(length(grep('subCancerX\\_Control\\_func',dir))==1){
    adj.name <- c('BMI','Sex')
  }
  
  if(length(grep('subCancerX\\-Ex\\_func',dir))==1|length(grep('subCancerX\\-Ex2\\_func',dir))==1){
    if(length(grep('prostate|breast',dir))==1){
      adj.name <- covars[!(covars %in% c(variable,"GI_nonGI","Charlson_score","Abx_last_month","Cancer_class","Sex"))]
    }else{
      adj.name <- covars[!(covars %in% c(variable,"GI_nonGI","Charlson_score","Abx_last_month","Cancer_class"))]
    }
  }
  
  print(adj.name)






#--------- Subset data -------
dist.obj$BC <- as.matrix(vegdist(t(data.obj$abund.list[[1]])))
  
## For unadj: exclude NA based on the variable of interest
tmp <- !is.na(data.obj$meta.dat[, c(variable), drop =F])
ind2 <- names(which(rowSums(tmp)==length(c(variable))))
data.obj2 <- subset_data(data.obj, ind2)
dist.obj2 <- subset_dist(dist.obj, ind2)
data.obj2$meta.dat <- droplevels(data.obj2$meta.dat) # otherwise ZicoSeq will fail, since levels still exist

## For adj: exclude NA based on the variable of interest and adjust variable
tmp <- !is.na(data.obj$meta.dat[, c(variable, adj.name), drop =F])
ind3 <- names(which(rowSums(tmp)==length(c(variable,adj.name))))
data.obj3 <- subset_data(data.obj, ind3)
dist.obj3 <- subset_dist(dist.obj, ind3)
data.obj3$meta.dat <- droplevels(data.obj3$meta.dat)


if(dir == 'CancerOnly_func/pathway' & variable == 'icd10_first_3_name'){
  ind <- names(which(table(data.obj2$meta.dat[, variable])>15))

  ind2 <- rownames(data.obj.rff2$meta.dat)[data.obj.rff2$meta.dat[,variable] %in% ind]
  data.obj.rff2 <- subset_data(data.obj.rff2, ind2)
  dist.obj.rff2 <- subset_dist(dist.obj.rff2, ind2)
  
  ind2 <- rownames(data.obj2$meta.dat)[data.obj2$meta.dat[,variable] %in% ind]
  data.obj2 <- subset_data(data.obj2, ind2)
  dist.obj2 <- subset_dist(dist.obj2, ind2)
  
  ind3 <- rownames(data.obj.rff3$meta.dat)[data.obj.rff3$meta.dat[,variable] %in% ind]
  data.obj.rff3 <- subset_data(data.obj.rff3, ind3)
  dist.obj.rff3 <- subset_dist(dist.obj.rff3, ind3)
  
  ind3 <- rownames(data.obj3$meta.dat)[data.obj3$meta.dat[,variable] %in% ind]
  data.obj3 <- subset_data(data.obj3, ind3)
  dist.obj3 <- subset_dist(dist.obj, ind3)
  
}

save(data.obj2, data.obj3, file = 'data.obj.wk2.RData')

## --------perform Beta diversity analysis---------
setwd(paste0(rd, dir))
dir3 <- 'Beta'
if (!dir.exists(dir3)) {dir.create(dir3)}
setwd(dir3)
getwd()

if (!dir.exists(variable)) {dir.create(variable)}
setwd(variable)


set.seed(123)
beta.measure=c('BC')

r2.unadj.mat <- pv.unadj.mat <- r2.adj.mat <- pv.adj.mat <-
  matrix(NA, length(beta.measure), length(variable), dimnames = list(beta.measure, variable))
for (dist.name in beta.measure) {
  # Marginal[use all rarified samples with grp.name no NA]
  obj <- dmanova(as.dist(dist.obj2[[dist.name]]) ~ data.obj2$meta.dat[, variable])
  if(dir %in% c('CancerOnly_func/pathway') & variable %in% covars){
    n <- nrow(data.obj2$meta.dat)
    r2.unadj.mat[dist.name, variable] <- 1 - (1 - obj$aov.tab[1, 5]) * (n - 1) / (n - 1 - obj$aov.tab[1, 1])
  }else{
    r2.unadj.mat[dist.name, variable] <- obj$aov.tab[1, 5]
  }
  pv.unadj.mat[dist.name, variable] <- obj$aov.tab[1, 6]



  # Adjust covariates[use rarefied samples with adj.name and grp.name no NA]
  dist.tmp <- as.dist(dist.obj3[[dist.name]])
  obj1 <- dmanova(as.formula(paste0('dist.tmp ~', paste0(adj.name,collapse = '+'))),
                  data = data.obj3$meta.dat)
  obj3 <- dmanova(as.formula(paste0('dist.tmp ~', paste0(adj.name,collapse = '+'),'+',variable)),
                  data = data.obj3$meta.dat)

  if(dir %in% c('CancerOnly_func/pathway')){
    tss <- obj1$aov.tab[2, 2]
    rss <- obj3$aov.tab[2, 2]
    n <- length(data.obj3$meta.dat[, variable])
    df.z <- ncol(model.matrix( ~ . - 1, data.obj3$meta.dat[, adj.name]))
    df.x <- ncol(model.matrix( ~ . - 1, data.obj3$meta.dat[, variable,drop =F]))
    r2.adj.mat[dist.name, variable] <- 1 - (n - 1 - df.z) * rss / (n - 1 - df.x - df.z) / tss
    }else{
      r2.adj.mat[dist.name, variable] <- obj3$aov.tab[1, 5]
    }
  pv.adj.mat[dist.name, variable] <- obj3$aov.tab[1, 6]
}
save(r2.adj.mat, r2.unadj.mat, pv.adj.mat, pv.unadj.mat, file = 'R2_pvalue.RData')



# # --------perform DAA analysis---------
setwd(paste0(rd, dir))
dir3 <- 'DAA'
if (!dir.exists(dir3)) {dir.create(dir3)}
setwd(dir3)
getwd()


if (!dir.exists(variable)) {dir.create(variable)}
setwd(variable)

temp <- t(t(data.obj3$abund.list[[1]])/colSums(data.obj3$abund.list[[1]]))
filter.ind <- rowMeans(temp != 0) >= 0.2 & rowMeans(temp) >= 0 & matrixStats::rowMaxs(temp) >= 0.002
data.obj3$abund.list[[1]] <- data.obj3$abund.list[[1]][filter.ind,,drop =F]

## perform DAA analysis
set.seed(123)
diff.obj <- perform_differential_analysis_zicoseq2(data.obj3, feature.dat.type = 'other',
                                                  taxa.levels = names(data.obj3$abund.list), max.abund.filter = 0,prev.filter = 0,
                                                  is.winsor = T, 
                                                  grp.name = variable, adj.name = adj.name, perm.no = 999, cutoff=0.05,
                                                  ann = '')
save(diff.obj, file = paste0(variable, '_ZicoSeq.Rdata'))

