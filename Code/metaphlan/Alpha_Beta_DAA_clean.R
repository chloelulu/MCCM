args=(commandArgs(TRUE))
if(length(args)==0){
  print("No arguments supplied.")
} else{
  variable = args[1]
  dir = args[2]
}

gc()

cat("[INFO] Performing:\n"); print(dir)

# ---------------------- Load Dependencies ----------------------
file_dir <- dirname(rstudioapi::getSourceEditorContext()$path)
file_dir <- dirname(dirname(dirname(dirname(dirname(file_dir)))))

source(paste0(file_dir,"/Code/Submission/MayoOncobiomeStudy/Code/Stats.R"))
try(load_package())

# ---------------------- Set Working Directory ----------------------
# mph_choice <- "mph_commonKraken"
mph_choice <- "mph"
rd <- paste0(file_dir,"/Result_",mph_choice,"/")

setwd(paste0(rd, dir))
load("data.obj.wk.RData")

# ---------------------- Define Covariates ----------------------
covars <- c("Batch", "Bristol_score", "BMI", "Age", "Sex", "GI_nonGI", "Cancer_class", "Metastasis",
            "PPI_day_365", "Abx_day_365", "Abx_last_month", "Charlson_score", "Elix_score",
            "Sample_season", "Urban")

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
  }else if(length(grep('early_onset',variable))==1){
    adj.name <- covars[!(covars %in% c(variable,'GI_nonGI',"Charlson_score","Abx_last_month","PPI_last_month","Cancer_class"))]
  }else {
    adj.name <- covars[!(covars %in% c(variable,'GI_nonGI',"Charlson_score","Abx_last_month","PPI_last_month"))]
  }
}

if(length(grep('EarlyOnset',dir))==1){
  if(length(grep('breast',dir))==1){
    adj.name <- covars[!(covars %in% c(variable,'GI_nonGI',"Charlson_score","Abx_last_month","PPI_last_month","Cancer_class","Sex","Age"))]
  }else{
    adj.name <- covars[!(covars %in% c(variable,'GI_nonGI',"Charlson_score","Abx_last_month","PPI_last_month","Cancer_class","Age"))]
  }
}

if(dir %in% c('PanCancer')){ adj.name <- c('BMI','Sex','Age')}

if(length(grep('subCancerX\\_Control',dir))==1){
  adj.name <- c('BMI','Sex','Age')
}

if(length(grep('subCancerX\\-Ex',dir))==1){
  if(length(grep('prostate|breast',dir))==1){
    adj.name <- covars[!(covars %in% c(variable,"GI_nonGI","Charlson_score","Abx_last_month","Cancer_class","Sex","PPI_last_month"))]
  }else{
    adj.name <- covars[!(covars %in% c(variable,"GI_nonGI","Charlson_score","Abx_last_month","Cancer_class","PPI_last_month"))]
  }
}
  
if(dir=='LungCancer_Smoking'){
  adj.name <- covars[!(covars %in% c(variable,'GI_nonGI',"Charlson_score","Abx_last_month","PPI_last_month","Cancer_class","Age"))]
}
cat("[INFO] Adjustment covariates:\n"); print(adj.name)

# ---------------------- Subset Data ----------------------
## For unadj: exclude NA based on the variable of interest
ind2 <- rownames(data.obj$meta.dat)[!(is.na(data.obj$meta.dat[, variable]))]
data.obj2 <- subset_data(data.obj, ind2)
dist.obj2 <- subset_dist(dist.obj, ind2)

## For adj: exclude NA based on the variable of interest and adjust variable
tmp <- !is.na(data.obj$meta.dat[, c(variable, adj.name), drop =F])
ind3 <- rowSums(tmp)==length(c(variable,adj.name))
data.obj3 <- subset_data(data.obj, ind3)
dist.obj3 <- subset_dist(dist.obj, ind3)

## For subCancerX, we only include 22 cancer types with n>15. 
if(dir =='CancerOnly' & grepl('icd10_first_3_name', variable)){
  ind <- names(which(table(data.obj2$meta.dat[, variable])>15))

  ind2 <- rownames(data.obj2$meta.dat)[data.obj2$meta.dat[,variable] %in% ind]
  data.obj2 <- subset_data(data.obj2, ind2)
  dist.obj2 <- subset_dist(dist.obj2, ind2)
  
  ind3 <- rownames(data.obj3$meta.dat)[data.obj3$meta.dat[,variable] %in% ind]
  data.obj3 <- subset_data(data.obj3, ind3)
  dist.obj3 <- subset_dist(dist.obj3, ind3)
  
  
}

## remember to drop levels, which influence the coding result of Alpha
if(class(data.obj3$meta.dat[,variable])=='factor')  data.obj3$meta.dat[,variable] <- droplevels(data.obj3$meta.dat[,variable])
if(class(data.obj2$meta.dat[,variable])=='factor')  data.obj2$meta.dat[,variable] <- droplevels(data.obj2$meta.dat[,variable])

# # ---------------------- Alpha diversity analysis ----------------------
# dir3 <- 'Alpha'
# if (!dir.exists(dir3)) {dir.create(dir3)}
# setwd(dir3)
# getwd()
# 
# if (!dir.exists(variable)) {dir.create(variable)}
# setwd(variable)
# 
# alpha.measure=c('Shannon')
# set.seed(123)
# alpha.obj1 <- generate_alpha_diversity(data.obj2, measures = alpha.measure, rarefy = F)
# alpha.obj2 <- generate_alpha_diversity(data.obj3, measures = alpha.measure, rarefy = F)
# fit1 <- perform_alpha_test3(data.obj=data.obj2, alpha.obj=alpha.obj1, measures = alpha.measure,
#                             grp.name=variable, adj.name=NULL, subject=NULL, ann= variable)
# fit2 <- perform_alpha_test3(data.obj=data.obj3, alpha.obj=alpha.obj2, measures = alpha.measure,
#                             grp.name=variable, adj.name=adj.name, subject=NULL, ann=paste0(variable,'_AdjAll'))
# save(fit1, fit2, alpha.obj1, alpha.obj2, data.obj2, data.obj3, file = 'Alpha.RData')
# 
# 
# # ---------------------- Beta Diversity Analysis ----------------------
# setwd(paste0(rd, dir))
# dir3 <- 'Beta'
# if (!dir.exists(dir3)) {dir.create(dir3)}
# setwd(dir3)
# getwd()
# 
# if (!dir.exists(variable)) {dir.create(variable)}
# setwd(variable)
# 
# 
# set.seed(123)
# beta.measure=c('WUniFrac')
# 
# r2.unadj.mat <- pv.unadj.mat <- r2.adj.mat <- pv.adj.mat <-
#   matrix(NA, length(beta.measure), length(variable), dimnames = list(beta.measure, variable))
# for (dist.name in beta.measure) {
#   # Marginal[use all rarefied samples with grp.name no NA]
#   obj <- dmanova(as.dist(dist.obj2[[dist.name]]) ~ data.obj2$meta.dat[, variable])
#   if(dir %in% c('CancerOnly')){
#     n <- nrow(data.obj2$meta.dat)
#     r2.unadj.mat[dist.name, variable] <- 1 - (1 - obj$aov.tab[1, 5]) * (n - 1) / (n - 1 - obj$aov.tab[1, 1])
#   }else{
#     r2.unadj.mat[dist.name, variable] <- obj$aov.tab[1, 5]
#   }
#   pv.unadj.mat[dist.name, variable] <- obj$aov.tab[1, 6]
# 
# 
#   # Adjust covariates[use rarefied samples with adj.name and grp.name no NA]
#   dist.tmp <- as.dist(dist.obj3[[dist.name]])
#   obj1 <- dmanova(as.formula(paste0('dist.tmp ~', paste0(adj.name,collapse = '+'))),
#                   data = data.obj3$meta.dat)
#   obj3 <- dmanova(as.formula(paste0('dist.tmp ~', paste0(adj.name,collapse = '+'),'+',variable)),
#                   data = data.obj3$meta.dat)
# 
#   if(dir %in% c('CancerOnly') & variable %in% covars){ # only covariates we used partial R2
#     tss <- obj1$aov.tab[2, 2]
#     rss <- obj3$aov.tab[2, 2]
#     n <- length(data.obj3$meta.dat[, variable])
#     df.z <- ncol(model.matrix( ~ . - 1, data.obj$meta.dat[, adj.name]))
#     df.x <- ncol(model.matrix( ~ . - 1, data.obj$meta.dat[, variable,drop =F]))
#     r2.adj.mat[dist.name, variable] <- 1 - (n - 1 - df.z) * rss / (n - 1 - df.x - df.z) / tss
#   }else{
#     r2.adj.mat[dist.name, variable] <- obj3$aov.tab[1, 5]
#   }
#   pv.adj.mat[dist.name, variable] <- obj3$aov.tab[1, 6]
# }
# save(r2.adj.mat, r2.unadj.mat, pv.adj.mat, pv.unadj.mat, file = 'R2_pvalue.RData')
# 
# 
# ---------------------- Differential Abundance Analysis (ZicoSeq) ----------------------
setwd(paste0(rd, dir))

dir3 <- 'DAA'
if (!dir.exists(dir3)) {dir.create(dir3)}
setwd(dir3)
getwd()

if (!dir.exists(variable)) {dir.create(variable)}
setwd(variable)


# ## If we only use the common species between metaphlan and kraken2
if(mph_choice=="mph_commonKraken"){
  load(paste0(file_dir,'/Data/data.obj.raw.core.RData'))
  filter.ind <- rowMeans(data.obj3$abund.list$Species != 0) > 0 & rowMaxs(data.obj3$abund.list$Species) > 0
  data.obj3$abund.list$Species <- data.obj3$abund.list$Species[filter.ind,,drop =F]
  filter.ind <- rownames(data.obj3$abund.list$Species) %in% rownames(data.obj$otu.tab)
  data.obj3$abund.list$Species <- data.obj3$abund.list$Species[filter.ind,,drop =F]/100
  
  filter.ind <- rowMeans(data.obj3$abund.list$Genus != 0) > 0 & rowMaxs(data.obj3$abund.list$Genus) > 0
  data.obj3$abund.list$Genus <- data.obj3$abund.list$Genus[filter.ind,,drop =F]
  filter.ind <- rownames(data.obj3$abund.list$Genus) %in% gsub('.*;g__','g__',rownames(data.obj$abund.list$Genus))
  data.obj3$abund.list$Genus <- data.obj3$abund.list$Genus[filter.ind,,drop =F]/100
}

## If we use all species found in metaphlan with same threthold used on kraken2 analysis
if(mph_choice=="mph"){
  filter.ind <- rowMeans(data.obj3$abund.list$Species != 0) >= 0.1 & rowMaxs(data.obj3$abund.list$Species) >= 0 # sum is 100 not 1, so should be 0.2
  data.obj3$abund.list$Species <- data.obj3$abund.list$Species[filter.ind,,drop =F]/100

  filter.ind <- rowMeans(data.obj3$abund.list$Genus != 0) >= 0.1 & rowMaxs(data.obj3$abund.list$Genus) >= 0
  data.obj3$abund.list$Genus <- data.obj3$abund.list$Genus[filter.ind,,drop =F]/100
}


set.seed(123)
diff.obj <- perform_differential_analysis_zicoseq2(data.obj3, feature.dat.type = 'other',
                                                   taxa.levels = c("Genus","Species"),
                                                   max.abund.filter = 0, prev.filter = 0,
                                                   is.winsor = F, is.post.sample = T,
                                                   grp.name = variable, adj.name = adj.name,
                                                   perm.no = 999, cutoff=0.1,
                                                   ann = '', plot =F)

save(diff.obj, file = paste0(variable, '_ZicoSeq.Rdata'))




# ## Testing for ZicoSeq and Maaslin3
# level <- "Species"
# load(paste0(variable, '_ZicoSeq.Rdata'))
# feature.dat <- data.obj3$abund.list[[level]]
# feature.dat <- feature.dat[rownames(diff.obj$qv.list[[level]]),,drop =F]
# feature.dat <- t(apply(feature.dat, 1, function(x) {
#   row_min <- min(x[x > 0], na.rm = TRUE)
#   x[x == 0] <- row_min/2
#   return(x)
# }))
# 
# fit.masslin3 <- fit.zicoseq <- NULL
# set.seed(123)
# fit.zicoseq <- GUniFrac::ZicoSeq(feature.dat = log2(feature.dat), meta = data.obj3$meta.dat,
#                                  feature.dat.type = "other",
#                                  perm.no = 999,
#                                  grp.name = variable, adj.name = adj.name,
#                                  prev.filter = 0, mean.abund.filter = 0, max.abund.filter = 0, is.winsor = F, is.post.sample = F,
#                                  link.func = list(function (x) x))
# 
# gc()
# set.seed(123)
# fit.masslin3 <- maaslin3::maaslin3(input_data = feature.dat,
#                          input_metadata = data.obj3$meta.dat,
#                          output = 'masslin3_output',
#                          formula = paste0('~',paste0(c(variable, adj.name), collapse = '+')),
#                          normalization = 'NONE',
#                          transform = 'LOG',
#                          augment = T,
#                          standardize = TRUE,
#                          max_significance = 0.1,
#                          median_comparison_abundance = F,
#                          median_comparison_prevalence = FALSE,
#                          warn_prevalence = FALSE,
#                          plot_summary_plot=F,
#                          save_models = F,
#                          cores = 1)
# res.masslin3 <- fit.masslin3$fit_data_abundance$results[fit.masslin3$fit_data_abundance$results$metadata==variable,c('feature','coef','qval_joint')]
# colnames(res.masslin3)[2:3] <- c('coef_masslin3','Qvalue_masslin3')
# 
# res.zicoseq <- merge(cbind(fit.zicoseq$p.adj.fdr), sign(t(fit.zicoseq$coef.list[[1]])[,grep(variable, colnames(t(fit.zicoseq$coef.list[[1]]))),drop =F]) * diff.obj$R2.list[[level]], by = 0) %>%
#   dplyr::rename(feature=Row.names)
# colnames(res.zicoseq)[2:3] <- c('Qvalue_ZicoSeq','coef_ZicoSeq')
# 
# res.zicoseq_old <- merge(diff.obj$qv.list[[level]], sign(diff.obj$coef.list[[level]][,grep(variable, colnames(diff.obj$coef.list[[level]])),drop =F]) * diff.obj$R2.list[[level]], by = 0) %>%
#   dplyr::rename(feature=Row.names)
# colnames(res.zicoseq_old)[2:3] <- c('Qvalue_ZicoSeq(old)','coef_ZicoSeq(old)')
# 
# all3 <- merge(res.masslin3, res.zicoseq,  by = 'feature') %>% merge(res.zicoseq_old,by = 'feature')
# colSums(all3[,grepl('Qvalue',colnames(all3))]<0.1)
# 
# save(fit.masslin3, file = paste0(level,'_masslin3.Rdata'))
# save(all3, file = paste0(level,'_DAA_all2_compare.Rdata'))
# write.csv(all3, file = paste0(level,'_DAA_all2_compare.csv'), row.names = F)
# gc()
