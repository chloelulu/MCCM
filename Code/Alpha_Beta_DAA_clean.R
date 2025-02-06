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
            "PPI_last_month","Abx_last_month","Charlson_score","Elix_score","Sample_season","Urban")
blood.names <- c('Erythrocytes','Hematocrit','Neutrophils','MCV','Hemoglobin','Leukocytes','Platelet.Count',
                 'Neutrophils_cat','Platelet.Count_cat','Hemoglobin_cat',
                 "neut_neutropenia2","bone_marrow_suppression2","Hb_anemia2","Pl_thrombocytopenia2","neut_neutropenia_c","Hb_anemia_c","Pl_thrombocytopenia_c")
elix.names <- c("Elixhauser_CHF","Elixhauser_Arrhythmia","Elixhauser_Valvular","Elixhauser_PHTN","Elixhauser_PVD","Elixhauser_HTN","Elixhauser_Paralysis",  
                "Elixhauser_NeuroOther","Elixhauser_Pulmonary","Elixhauser_DM","Elixhauser_DMcx","Elixhauser_Hypothyroid","Elixhauser_Renal","Elixhauser_Liver",      
                "Elixhauser_PUD" ,"Elixhauser_HIV","Elixhauser_Lymphoma","Elixhauser_Mets","Elixhauser_Tumor","Elixhauser_Rheumatic","Elixhauser_Coagulopathy",
                "Elixhauser_Obesity","Elixhauser_WeightLoss","Elixhauser_FluidsLytes","Elixhauser_BloodLoss","Elixhauser_Anemia","Elixhauser_Alcohol","Elixhauser_Drugs",
                "Elixhauser_Psychoses","Elixhauser_Depression")

if(dir %in% c('CancerOnly')){
  if(variable =='Cancer_class'){
    adj.name <- covars[!(covars %in% c(variable,'GI_nonGI',"Charlson_score","Abx_last_month","PPI_last_month"))]
  }else if(variable =='GI_nonGI'){
    adj.name <- covars[!(covars %in% c(variable,'Cancer_class',"Charlson_score","Abx_last_month","PPI_last_month"))]
  }else if(variable =="Charlson_score"){
    adj.name <- covars[!(covars %in% c(variable,"Elix_score",'Cancer_class',"Abx_last_month","PPI_last_month"))]
  }else if(variable =="Elix_score"){
    adj.name <- covars[!(covars %in% c(variable,"Charlson_score",'Cancer_class',"Abx_last_month","PPI_last_month"))]
  }else if(variable =="Abx_last_month"){
    adj.name <- covars[!(covars %in% c(variable,'Cancer_class',"Charlson_score","Abx_day_365","Abx_last_month","PPI_last_month"))]
  }else if(variable =="PPI_last_month"){
    adj.name <- covars[!(covars %in% c(variable,'Cancer_class',"Charlson_score","PPI_day_365","PPI_last_month","Abx_last_month"))]
  }else if(variable %in% c('Elixhauser_Mets','Charlson_Mets')){# exclude Metastasis
    adj.name <- covars[!(covars %in% c(variable,'GI_nonGI',"Charlson_score", "Elix_score","Metastasis","Abx_last_month","PPI_last_month"))]
  }else if(variable %in% c('Elixhauser_Obesity','Elixhauser_WeightLoss')){#exclude BMI
    adj.name <- covars[!(covars %in% c(variable,'GI_nonGI',"Charlson_score", "Elix_score","BMI","Abx_last_month","PPI_last_month"))]
  }else if(length(grep('^Elixhauser|^Charlson',variable))==1){
    adj.name <- covars[!(covars %in% c(variable,'GI_nonGI',"Charlson_score","Abx_last_month","Elix_score","PPI_last_month"))]
  }else if(length(grep('icd10_first_3_name',variable))==1){
    adj.name <- covars[!(covars %in% c(variable,'GI_nonGI',"Charlson_score","Abx_last_month","PPI_last_month","Cancer_class"))]
  }else {
    adj.name <- covars[!(covars %in% c(variable,'GI_nonGI',"Charlson_score","Abx_last_month","PPI_last_month"))]
  }
}

# currently PanCancer Only contains Batch 5, We aim at checking how is the difference if we control batch. And we drop previous Pancancer
if(dir %in% c('PanCancer')){ adj.name <- c('BMI','Sex') }

if(length(grep('subCancerX\\_Control',dir))==1){
  adj.name <- c('BMI','Sex')
}

if(length(grep('subCancerX\\-Ex',dir))==1){
  if(length(grep('prostate|breast',dir))==1){
    adj.name <- covars[!(covars %in% c(variable,"GI_nonGI","Charlson_score","Abx_last_month","Cancer_class","Sex","PPI_last_month"))]
  }else{
    adj.name <- covars[!(covars %in% c(variable,"GI_nonGI","Charlson_score","Abx_last_month","Cancer_class","PPI_last_month"))]
  }
}
  


#--------- Subset data -------
## For unadj: exclude NA based on the variable of interest
ind2 <- !(is.na(data.obj.rff$meta.dat[, variable]))
data.obj.rff2 <- subset_data(data.obj.rff, ind2)
dist.obj.rff2 <- subset_dist(dist.obj.rff, ind2)

ind2 <- !(is.na(data.obj$meta.dat[, variable]))
data.obj2 <- subset_data(data.obj, ind2)

## For adj: exclude NA based on the variable of interest and adjust variable
tmp <- !is.na(data.obj.rff$meta.dat[, c(variable, adj.name), drop =F])
ind3 <- rowSums(tmp)==length(c(variable,adj.name))
data.obj.rff3 <- subset_data(data.obj.rff, ind3)
dist.obj.rff3 <- subset_dist(dist.obj.rff, ind3)

tmp <- !is.na(data.obj$meta.dat[, c(variable, adj.name), drop =F])
ind3 <- rowSums(tmp)==length(c(variable,adj.name))
data.obj3 <- subset_data(data.obj, ind3)


## For subCancerX, we only include 21 cancer types with n>15. Also exclude 2 other types Ruben mentioned
if(dir =='CancerOnly'&variable=='icd10_first_3_name'){
  ind <- names(which(table(data.obj2$meta.dat[, variable])>15))
  ind <- ind[-grep('ther',ind)]
  
  ind2 <- rownames(data.obj.rff2$meta.dat)[data.obj.rff2$meta.dat[,variable] %in% ind]
  data.obj.rff2 <- subset_data(data.obj.rff2, ind2)
  dist.obj.rff2 <- subset_dist(dist.obj.rff2, ind2)
  
  ind2 <- rownames(data.obj2$meta.dat)[data.obj2$meta.dat[,variable] %in% ind]
  data.obj2 <- subset_data(data.obj2, ind2)
  # dist.obj2 <- subset_dist(dist.obj2, ind2)
  
  ind3 <- rownames(data.obj.rff3$meta.dat)[data.obj.rff3$meta.dat[,variable] %in% ind]
  data.obj.rff3 <- subset_data(data.obj.rff3, ind3)
  dist.obj.rff3 <- subset_dist(dist.obj.rff3, ind3)
  
  ind3 <- rownames(data.obj3$meta.dat)[data.obj3$meta.dat[,variable] %in% ind]
  data.obj3 <- subset_data(data.obj3, ind3)
  # dist.obj3 <- subset_dist(dist.obj, ind3)
  
}




## --------perform Alpha diversity analysis---------
dir3 <- 'Alpha'
if (!dir.exists(dir3)) {dir.create(dir3)}
setwd(dir3)
getwd()

if (!dir.exists(variable)) {dir.create(variable)}
setwd(variable)

alpha.measure=c('Observed', 'Chao1', 'Shannon', 'InvSimpson')
set.seed(123)
alpha.obj1 <- generate_alpha_diversity(data.obj.rff2, measures = alpha.measure, rarefy = F)
alpha.obj2 <- generate_alpha_diversity(data.obj.rff3, measures = alpha.measure, rarefy = F)
fit1 <- perform_alpha_test3(data.obj=data.obj.rff2, alpha.obj=alpha.obj1, measures = alpha.measure,
                            grp.name=variable, adj.name=NULL, subject=NULL, ann= variable)
fit2 <- perform_alpha_test3(data.obj=data.obj.rff3, alpha.obj=alpha.obj2, measures = alpha.measure,
                            grp.name=variable, adj.name=adj.name, subject=NULL, ann=paste0(variable,'_AdjAll'))
save(fit1, fit2, alpha.obj1, alpha.obj2, data.obj.rff2, data.obj.rff3, file = 'Alpha.RData')








#--------perform Beta diversity analysis---------
setwd(paste0(rd, dir))
dir3 <- 'Beta'
if (!dir.exists(dir3)) {dir.create(dir3)}
setwd(dir3)
getwd()

if (!dir.exists(variable)) {dir.create(variable)}
setwd(variable)


set.seed(123)
beta.measure=c('UniFrac', 'GUniFrac', 'WUniFrac', 'BC')

r2.unadj.mat <- pv.unadj.mat <- r2.adj.mat <- pv.adj.mat <-
  matrix(NA, length(beta.measure), length(variable), dimnames = list(beta.measure, variable))
for (dist.name in beta.measure) {
  # Marginal[use all rarefied samples with grp.name no NA]
  obj <- dmanova(as.dist(dist.obj.rff2[[dist.name]]) ~ data.obj.rff2$meta.dat[, variable])
  if(dir %in% c('CancerOnly')){
    n <- nrow(data.obj.rff2$meta.dat)
    r2.unadj.mat[dist.name, variable] <- 1 - (1 - obj$aov.tab[1, 5]) * (n - 1) / (n - 1 - obj$aov.tab[1, 1])
  }else{
    r2.unadj.mat[dist.name, variable] <- obj$aov.tab[1, 5]
  }
  pv.unadj.mat[dist.name, variable] <- obj$aov.tab[1, 6]



  # Adjust covariates[use rarefied samples with adj.name and grp.name no NA]
  dist.tmp <- as.dist(dist.obj.rff3[[dist.name]])
  obj1 <- dmanova(as.formula(paste0('dist.tmp ~', paste0(adj.name,collapse = '+'))),
                  data = data.obj.rff3$meta.dat)
  obj3 <- dmanova(as.formula(paste0('dist.tmp ~', paste0(adj.name,collapse = '+'),'+',variable)),
                  data = data.obj.rff3$meta.dat)

  if(dir %in% c('CancerOnly') & variable %in% covars){ # only covariates we used partial R2
    tss <- obj1$aov.tab[2, 2]
    rss <- obj3$aov.tab[2, 2]
    n <- length(data.obj.rff3$meta.dat[, variable])
    df.z <- ncol(model.matrix( ~ . - 1, data.obj.rff3$meta.dat[, adj.name]))
    df.x <- ncol(model.matrix( ~ . - 1, data.obj.rff3$meta.dat[, variable,drop =F]))
    r2.adj.mat[dist.name, variable] <- 1 - (n - 1 - df.z) * rss / (n - 1 - df.x - df.z) / tss
  }else{
    r2.adj.mat[dist.name, variable] <- obj3$aov.tab[1, 5]
  }
  pv.adj.mat[dist.name, variable] <- obj3$aov.tab[1, 6]
}
save(r2.adj.mat, r2.unadj.mat, pv.adj.mat, pv.unadj.mat, file = 'R2_pvalue.RData')



# PCoA plot [use all rarified samples with grp.name no NA]
pdf('PCoA.pdf', width = 8,height = 6)
for(dist.name in beta.measure){
  obj <- cmdscale(as.dist(dist.obj.rff2[[dist.name]]), k=2, eig=T)
  pve <- round(obj$eig[1:2]/sum(abs(obj$eig))*100, 1)
  y <- cbind.data.frame(PC1=obj$points[, 1], PC2=obj$points[, 2])
  xlab <- paste0('PC1(', pve[1], '%)')
  ylab <- paste0('PC2(', pve[2], '%)')
  yy <- merge(y, data.obj.rff2$meta.dat[,variable, drop =F], by = 0)

  if(class(data.obj.rff2$meta.dat[,variable])=="numeric"){
    cutoff <- median(yy[,variable])
    yy[,variable] <- ifelse(yy[,variable]>cutoff, '>median', '<median')
    yy[,variable] <- as.factor(yy[,variable])
  }
  p <- ggplot(yy, aes(x = PC1, y = PC2, color = !!as.name(variable))) +
    geom_point(size =0.7) +
    ggforce::geom_mark_ellipse(aes(x = PC1, y = PC2, color= !!as.name(variable)),expand = unit(0, "mm"))+
    theme_bw() +
    scale_x_continuous(limits = c(min(yy$PC1)*1.2,max(yy$PC1)*1.2))+
    scale_y_continuous(limits = c(min(yy$PC2)*1.2,max(yy$PC2)*1.2))+
    scale_color_brewer(palette = 'Dark2')+
    labs(x = xlab, y = ylab, title = dist.name) +
    theme(text = element_text(size = 20, color = "black"),
          axis.text = element_text(size = 20, color = "black"),
          legend.text = element_text(size = 20, color = "black"),
          axis.title =  element_text(size = 20, color = "black"),
          panel.grid.major = element_blank(),
          panel.border = element_rect(size = 2),
          axis.ticks = element_line(size = 1),
          axis.ticks.length = unit(0.2, 'cm'),
          panel.grid.minor = element_blank())
  print(p)
}
dev.off()



if(class(data.obj.rff3$meta.dat[,variable])=="factor"){
  generate_ordination(data.obj.rff2, dist.obj.rff2, dist.names = beta.measure, grp.name=variable, strata=NULL)
  generate_distance_boxplot(data.obj.rff2, dist.obj.rff2, dist.names = beta.measure, grp=variable, within=T, strata=NULL)
}







# --------perform DAA analysis---------
setwd(paste0(rd, dir))

dir3 <- 'DAA'
if (!dir.exists(dir3)) {dir.create(dir3)}
setwd(dir3)
getwd()

if (!dir.exists(variable)) {dir.create(variable)}
setwd(variable)

## perform DAA analysis
set.seed(123)
diff.obj <- perform_differential_analysis_zicoseq(data.obj3,
                                                  taxa.levels = c("Phylum","Class","Order","Family","Genus","Species"),
                                                  grp.name = variable, adj.name = adj.name, perm.no = 999, cutoff=0.05,
                                                  prev.filter = 0.2, max.abund.filter = 0.002,
                                                  ann = '')
save(diff.obj, file = paste0(variable, '_ZicoSeq.Rdata'))

