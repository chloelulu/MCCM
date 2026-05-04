library(survival)
library(dplyr)
library(tibble)
library(readxl)
library(dplyr)
library(survminer)
library(MiRKAT)
library(cowplot)
library(mice)
library(tidyverse)
library(corrplot)
library(ggplot2)
library(ggcorrplot)


#use the following colors in all figures
#"#0072B2", "#CC79A7"
elix.names <- c("Elixhauser_CHF","Elixhauser_Arrhythmia","Elixhauser_Valvular","Elixhauser_PHTN","Elixhauser_PVD","Elixhauser_HTN",
                "Elixhauser_NeuroOther","Elixhauser_Pulmonary","Elixhauser_DM","Elixhauser_DMcx","Elixhauser_Hypothyroid","Elixhauser_Renal","Elixhauser_Liver",      
                "Elixhauser_PUD" ,"Elixhauser_Lymphoma","Elixhauser_Mets","Elixhauser_Tumor","Elixhauser_Rheumatic","Elixhauser_Coagulopathy",
                "Elixhauser_Obesity","Elixhauser_WeightLoss","Elixhauser_FluidsLytes","Elixhauser_Anemia","Elixhauser_Alcohol","Elixhauser_Drugs",
                "Elixhauser_Depression")
variable_names <- c(
  "Elixhauser_CHF", "Elixhauser_Arrhythmia", "Elixhauser_Valvular", "Elixhauser_PHTN", 
  "Elixhauser_PVD", "Elixhauser_HTN", "Elixhauser_NeuroOther", 
  "Elixhauser_Pulmonary", "Elixhauser_DM", "Elixhauser_DMcx", "Elixhauser_Hypothyroid", 
  "Elixhauser_Renal", "Elixhauser_Liver", "Elixhauser_PUD", #"Elixhauser_HIV", 
  "Elixhauser_Lymphoma", "Elixhauser_Mets", "Elixhauser_Tumor", "Elixhauser_Rheumatic", 
  "Elixhauser_Coagulopathy", "Elixhauser_Obesity", "Elixhauser_WeightLoss", 
  "Elixhauser_FluidsLytes","Elixhauser_Anemia", 
  "Elixhauser_Alcohol", "Elixhauser_Drugs", "Elixhauser_Depression", "Elix_score", "Charlson_score"
)

detailed_descriptions <- c(
  "Congestive heart failure", "Cardiac arrhythmias", "Valvular disease", "Pulmonary hypertension",
  "Peripheral vascular disorders", "Hypertension (uncomplicated)", "Other neurologic disorders", 
  "Chronic pulmonary disease", "Diabetes without chronic complications", "Diabetes with chronic complications", "Hypothyroidism", 
  "Renal failure", "Liver disease", "Peptic ulcer disease (excluding bleeding)", #"AIDS/HIV infection", 
  "Lymphoma", "Metastatic cancer", "Solid tumor without metastasis", "Rheumatoid arthritis/collagen vascular diseases",
  "Coagulopathy", "Obesity", "Weight loss", 
  "Fluid and electrolyte disorders", "Deficiency anemia", 
  "Alcohol abuse", "Drug abuse", "Depression","ECI score", "CCI score"
)
names_map <- setNames(detailed_descriptions, variable_names)


file_dir <- dirname(rstudioapi::getSourceEditorContext()$path)
file_dir <- dirname(dirname(dirname(dirname(file_dir))))
setwd(file_dir) 
wd <- file_dir
mfd <- paste0(file_dir, "/ManuscriptFigures/")

setwd(wd)
source('Code/Stats.R')
load_package()

covars <- c("Batch", "Bristol_score", "BMI", "Age", "Sex","icd10_first_3_name_short", "Metastasis",
            "PPI_day_365", "Abx_day_365", "Elix_score","Sample_season", "Urban")
load('Data/data.obj.pathway.RData')





### BMI NAs using multiple imputation
md.pattern(data.obj$meta.dat[, c("BMI", "Age", "Sex","icd10_first_3_name_short")])
imp <- mice(data.obj$meta.dat[, c("BMI", "Age", "Sex","icd10_first_3_name_short")], m = 5, method = "pmm", maxit = 5,  seed = 123)
densityplot(imp)  
stripplot(imp, BMI ~ .imp, pch = 20, cex = 1.2)
summary(data.obj$meta.dat$BMI)
summary(complete(imp, 1)$BMI)

imp <- complete(imp)

sum(is.na(imp$BMI))
orig <- data.obj$meta.dat$BMI
imputed <- imp$BMI
idx <- !is.na(orig)
sum(orig[idx] != imputed[idx])
data.obj$meta.dat$BMI <- imp$BMI 


obj <- cmdscale(dist.obj[["BC"]], k = 10, eig = TRUE)
pve <- round(obj$eig[1:10] / sum(abs(obj$eig)) * 100, 1);sum(pve)
eig <- obj$points
colnames(eig) <- paste0("PC",1:10)
data.obj$meta.dat <- data.obj$meta.dat %>% rownames_to_column("SampleID")%>% inner_join(as.data.frame(eig) %>% rownames_to_column("SampleID"))
data.obj$meta.dat <- data.obj$meta.dat %>% column_to_rownames('SampleID')


# Keep cancer types with n>=15
tab <- table(
  data.obj$meta.dat$icd10_first_3_name_short,
  data.obj$meta.dat$vital_status
)
cancer_15 <- rownames(tab)[
  rowSums(tab) > 15 &
    tab[, "0"] > 2 &
    tab[, "1"] > 2
]

cancer_30 <- rownames(tab)[
  rowSums(tab) > 30 &
    tab[, "0"] > 2 &
    tab[, "1"] > 2
]

cancer_20 <- rownames(tab)[
  rowSums(tab) > 20 &
    tab[, "0"] > 2 &
    tab[, "1"] > 2
]


## keep cancer only group
#remove NA in vital_status from the cancer group

samIDs <- rownames(data.obj$meta.dat[data.obj$meta.dat$Group == "Cancer" & !is.na(data.obj$meta.dat$vital_status),])
length(samIDs) #1345; so some are filtered out now

data.obj2 <- subset_data(data.obj, samIDs = samIDs)
data.obj2$meta.dat <- droplevels(data.obj2$meta.dat)

## change small sample size cancers into other
data.obj2$meta.dat$icd10_first_3_name_short <- as.character(data.obj2$meta.dat$icd10_first_3_name_short)
data.obj2$meta.dat$icd10_first_3_name_short[!(data.obj2$meta.dat$icd10_first_3_name_short %in% cancer_15)] <- 'Other'
data.obj2$meta.dat$icd10_first_3_name_short <- as.factor(data.obj2$meta.dat$icd10_first_3_name_short)
table(data.obj2$meta.dat$icd10_first_3_name_short, data.obj2$meta.dat$vital_status) 

dim(data.obj2$meta.dat)
# 1345  206


## ----Method 2---
adj.name0 <- c("icd10_first_3_name_short", "Elix_score")
pval <- c()
for(var in covars[!(covars %in% adj.name0)]){
  full_formula <- paste0('Surv(OS_months, vital_status) ~',paste(c(adj.name0,var), collapse = " + "))
  full_model <- survival::coxph(as.formula(full_formula), data = data.obj2$meta.dat)
  baseline_formula <- paste0('Surv(OS_months, vital_status) ~',paste(adj.name0, collapse = " + "))
  baseline_model <- survival::coxph(as.formula(baseline_formula), data = data.obj2$meta.dat)
  if(class(data.obj$meta.dat[[var]])=="numeric" | length(levels(data.obj$meta.dat[[var]]))==2){
    coefs <- summary(full_model)$coef
    pval <- c(pval,coefs[grep(var,rownames(coefs)), "Pr(>|z|)"])
  }else{
    lrt <- anova(baseline_model, full_model)
    pval <- c(pval,lrt[2, 4])
  }
}
names(pval) <- covars[!(covars %in% adj.name0)]
adj.name <- c(adj.name0,names(which(pval<0.05)))
# Age,Batch,BMI,Elix_score,icd10_first_3_name_short,Metastasis



# ===== 1. pathway ===== 
## clr transform species abundance
prop <- t(t(data.obj2$abund.list$pathway)/colSums(data.obj2$abund.list$pathway))
idx <- (rowMeans(prop>0) > 0.1 & Biobase::rowMax(prop)> 0.002)
Y <- data.obj2$abund.list$pathway[idx,]
N <- colSums(Y)
m <- nrow(Y)
N.mat <- matrix(rep(N, m), nrow = m, byrow = TRUE)
N.mat[Y > 0] <- 0
Y <- Y + N.mat / N[max.col(N.mat)]
logY <- log2(Y)
W <- (t(logY) - colMeans(logY))


#adjusted for cancer class
species_results <- list() #1519
for(species_name in colnames(W)) {
  form <- as.formula(paste(
    "Surv(OS_months, vital_status) ~", 
    paste(c("taxon", adj.name), collapse = " + ")
  ))
  combined_data <- merge(as.data.frame(W)[,species_name, drop =F], data.obj2$meta.dat, by = 0)
  colnames(combined_data)[colnames(combined_data)==species_name] <- "taxon"
  fit <- coxph(form, data = combined_data)
  species_results[[species_name]] <- summary(fit)
}


species_pvals <- lapply(names(species_results), function(pc) {
  s <- species_results[[pc]]
  tibble(
    PC       = pc,
    beta     = s$coef["taxon", "coef"],
    HR       = s$coef["taxon", "exp(coef)"],
    p_value  = s$coef["taxon", "Pr(>|z|)"],
    lower95  = s$conf.int["taxon", "lower .95"],
    upper95  = s$conf.int["taxon", "upper .95"],
    n        = s$n,
    events   = s$nevent
  )
}) %>% bind_rows()

species_pvals <- species_pvals %>% mutate(padj = p.adjust(species_pvals$p_value,"BH"),
                                          `Adjust covariates`=paste0(adj.name, collapse = ','))
colnames(species_pvals)[1] <- 'Pathway'
# write.csv(species_pvals, file = 'Code/Submission/Supplementary tables/Survial_pathway_pancancer.csv', row.names = F)
species_pvals[species_pvals$padj < 0.1,]
#none

#very small Hazard ratio's and not significant after multiple testing correction
species_pvals[species_pvals$p_value < 0.01,]
head(species_pvals[order(species_pvals$HR, decreasing = T),]) 
head(species_pvals[order(species_pvals$HR, decreasing = F),])



# ===== 2. PC===== 
#no Unifrac PC significantly predictive of survival when adjusting for icd10_first_3_name_short + 
#Elix_score + Batch + BMI + Age + Metastasis

pc_results <- list()
for(pc in paste0("PC", 1:10)) {  
  form <- as.formula(paste(
    "Surv(OS_months, vital_status) ~", 
    paste(c(pc, adj.name), collapse = " + ")
  ))
  
  fit <- coxph(form, data = data.obj2$meta.dat) #has the PCs
  pc_results[[pc]] <- summary(fit)
}


pc_pvals <- lapply(names(pc_results), function(pc) {
  s <- pc_results[[pc]]
  tibble(
    PC       = pc,
    beta     = s$coef[pc, "coef"],
    HR       = s$coef[pc, "exp(coef)"],
    p_value  = s$coef[pc, "Pr(>|z|)"],
    lower95  = s$conf.int[pc, "lower .95"],
    upper95  = s$conf.int[pc, "upper .95"],
    n        = s$n,
    events   = s$nevent
  )
}) %>% bind_rows() %>% mutate(qval = p.adjust(p_value, 'BH'))

pc_pvals



pdf(paste0(wd,"/Figure/RM_figures/survival_curve_PCs_pathway.pdf"), width = 6, height = 6)
for(species_name in pc_pvals$PC) {
  combined_data <- data.obj2$meta.dat
  combined_data$group <- ifelse(
    combined_data[[species_name]] > median(combined_data[[species_name]]),
    "High","Low")
  
  fit <- survfit(Surv(OS_months, vital_status) ~ group, data = combined_data)
  p <- ggsurvplot(
    fit,
    data = combined_data,
    pval = F,
    conf.int = TRUE,
    risk.table = TRUE,
    risk.table.height = 0.25,
    title = "",
    xlab = "Time (months)",
    ylab = "Overall Survival Probability",
    legend.title = species_name,
    ggtheme = theme_minimal()
  ) 
  
  combined <- plot_grid(
    p$plot,
    p$table,
    ncol        = 1,
    rel_heights = c(2, 1)   # adjust relative heights as needed
  )
  
  print(combined)
}
dev.off()



#=====MiRKAT survival used for testing overall community with survival======

#overall unifrac dissimilarity is also not predictive of survival when adjusting for 
#icd10_first_3_name_short + #Elix_score + Batch + BMI + Age + Metastasis

idx <- names(which(colSums(apply(data.obj2$meta.dat[,c("OS_months", "vital_status", adj.name)], 1, function(x) !is.na(x))) == length(c("OS_months", "vital_status", adj.name))))
data.obj3 <- subset_data(data.obj2, samIDs = idx)
data.obj3$meta.dat <- droplevels(data.obj3$meta.dat)
dist.obj3 <- subset_dist(dist.obj, samIDs = idx)
dim(dist.obj3$BC) 
nrow(data.obj3$meta.dat)
model_mat <- model.matrix(as.formula(paste0('~ ', paste0(adj.name, collapse = '+'))),
                          data = data.obj3$meta.dat)[, -1] 
MiRKATS(obstime = data.obj3$meta.dat$OS_months, delta = data.obj3$meta.dat$vital_status, 
        X = model_mat, Ks = dist.obj3$BC, beta = NULL)
#$p_values
#[1] 0.8983204




#----------------------------------
#try MiRKAT on only single cancer classes

adj.name_sub <- c("Elix_score", "BMI", "Age")

p_MiRKATS_list <- list()
for (i in 1:length(unique(data.obj3$meta.dat$icd10_first_3_name_short))) {
  
  meta_sub <- data.obj3$meta.dat[which(data.obj3$meta.dat$icd10_first_3_name_short == unique(data.obj3$meta.dat$icd10_first_3_name_short)[i]),]
  dist_sub <- dist.obj3$BC[rownames(meta_sub), rownames(meta_sub)]
  
  model_mat <- model.matrix(as.formula(paste0('~ ', paste0(adj.name_sub, collapse = '+'))),
                            data = meta_sub)[, -1] 
  p_MiRKATS_list[[i]] <- MiRKATS(obstime = meta_sub$OS_months, delta = meta_sub$vital_status, 
                                 X = model_mat, Ks = dist_sub, beta = NULL)
}
min(unlist(lapply(p_MiRKATS_list, unlist))) #0.09437205, minimal p value so non significant


##====== subset each cancer type =====
## Due to the sample size of each cancer cohort, I only include sample size > 30 cancers
adj.name <- adj.name[!(adj.name %in% "icd10_first_3_name_short")]
res <- adj.names <- list()
species_pvals_all <- NULL
for(i in seq(cancer_30)){
  samIDs <- rownames(data.obj$meta.dat[data.obj$meta.dat$icd10_first_3_name_short %in% cancer_30[i],])
  data.obj2 <- subset_data(data.obj, samIDs = samIDs)
  data.obj2$meta.dat <- droplevels(data.obj2$meta.dat)
  
  # CLR transform
  prop <- t(t(data.obj2$abund.list$pathway)/colSums(data.obj2$abund.list$pathway))
  idx <- (rowMeans(prop>0) > 0.1 & Biobase::rowMax(prop)> 0.002)
  Y <- data.obj2$abund.list$pathway[idx,]
  N <- colSums(Y)
  m <- nrow(Y)
  N.mat <- matrix(rep(N, m), nrow = m, byrow = TRUE)
  N.mat[Y > 0] <- 0
  Y <- Y + N.mat / N[max.col(N.mat)]
  logY <- log2(Y)
  W <- (t(logY) - colMeans(logY))
  
  ## only keep significant variables we have identified in main text
  species_results <- list()
  for(species_name in colnames(W)) {
    if(length(grep('prostate|breast|uteri|ovary',cancer_30[i]))==1){
      adj.name <- adj.name[!(adj.name %in% "Sex")]
    }
    
    meta <- data.obj2$meta.dat[,c(adj.name, 'OS_months', 'vital_status')] %>% na.omit()
    combined_data <- merge(as.data.frame(W)[,species_name, drop =F], meta, by = 0)
    colnames(combined_data)[colnames(combined_data)==species_name] <- "taxon"
    
    form <- as.formula(paste(
      "Surv(OS_months, vital_status) ~", 
      paste(c("taxon", adj.name), collapse = " + ")
    ))
    fit <- coxph(form, data = combined_data)
    species_results[[species_name]] <- summary(fit)
  }
  
  species_pvals <- lapply(names(species_results), function(pc) {
    s <- species_results[[pc]]
    tibble(
      PC       = pc,
      beta     = s$coef["taxon", "coef"],
      HR       = s$coef["taxon", "exp(coef)"],
      p_value  = s$coef["taxon", "Pr(>|z|)"],
      lower95  = s$conf.int["taxon", "lower .95"],
      upper95  = s$conf.int["taxon", "upper .95"],
      n        = s$n,
      events   = s$nevent
    )
  }) %>% bind_rows() 
  
  species_pvals <- species_pvals %>% mutate(padj = p.adjust(species_pvals$p_value,"BH"))
  
  species_pvals_all <- rbind(species_pvals_all, species_pvals %>% mutate(cancer = cancer_30[i]))
  if(sum(species_pvals$padj<0.1)>0){
    selected_species <- species_pvals %>%
      dplyr::filter(padj < 0.1) %>%
      pull(PC)
    
    pdf(paste0(wd,"/Figure/RM_figures/survival_curve_species_FDR0.1_",cancer_30[i],"_pathway.pdf"), width = 6, height = 6)
    for(species_name in selected_species) {
      combined_data_prop <- merge(
        as.data.frame(t(prop))[, species_name, drop = FALSE], 
        data.obj2$meta.dat, 
        by = 0
      )
      combined_data <- merge(
        as.data.frame(W)[, species_name, drop = FALSE], 
        data.obj2$meta.dat, 
        by = 0
      )
      
      if(median(combined_data_prop[[species_name]])==0){# if median is 0, use presence/absence
        combined_data$abundance <- ifelse(combined_data_prop[[species_name]]==0,'Absence','Presence')
      }else{ # otherwise use median
        combined_data$abundance <- ifelse(
          combined_data[[species_name]] > median(combined_data[[species_name]]),
          "High","Low")
      }
      
      
      
      fit <- survfit(Surv(OS_months, vital_status) ~ abundance, data = combined_data)
      p <- ggsurvplot(
        fit,
        data = combined_data,
        pval = F,
        conf.int = TRUE,
        risk.table = TRUE,
        risk.table.height = 0.25,
        title = "",
        xlab = "Time (months)",
        ylab = "Overall Survival Probability",
        legend.title = species_name,
        palette = c("#0072B2", "#CC79A7"),
        ggtheme = theme_minimal()
      ) 
      
      combined <- plot_grid(
        p$plot,
        p$table,
        ncol        = 1,
        rel_heights = c(2, 1)   # adjust relative heights as needed
      )
      print(combined)
    }
    dev.off()
  }
  res[[cancer_30[[i]]]] <- species_pvals %>% mutate(adj = paste0(adj.name, collapse = ','))
}

sum(species_pvals_all$padj<0.1)

colnames(species_pvals_all)[1] <- 'Pathway'
species_pvals_all <- species_pvals_all %>% mutate(`Adjust covariates`=paste0(adj.name, collapse = ',')) 



## prepare Supplementary Table for Survival Pathway
library(openxlsx)
wb <- loadWorkbook("Code/Submission/Supplementary tables/20260424 revision Supplement/Table S8_new.xlsx")
names(wb)
addWorksheet(wb, "TableS8-3")
writeData(wb, sheet = "TableS8-3", x = species_pvals, startRow = 4)

addWorksheet(wb, "TableS8-5")
writeData(wb, sheet = "TableS8-5", x = species_pvals_all)

saveWorkbook(wb, "Code/Submission/Supplementary tables/20260424 revision Supplement/Table S8_new.xlsx", overwrite = TRUE)


res_df <- bind_rows(res, .id = "cancer_class")
res_df <- res_df %>%
  dplyr::filter(padj < 0.1) 
names(which(table(res_df$cancer_class)>1))

plts <- sig_species_sum <- list()
for (i in names(which(table(res_df$cancer_class) > 0))) {
  sig_species <- res_df[res_df$cancer_class == i, ]
  
  # order PC by median HR within this cancer_class
  sig_species <- sig_species %>%
    group_by(PC) %>%
    mutate(HR_median = median(HR, na.rm = TRUE)) %>%
    ungroup() %>%
    mutate(PC = reorder(PC, -HR_median))
  
  sig_species_sum[[i]] <- sig_species
  p1 <- ggplot(sig_species, aes(x = HR, y = PC)) +
    geom_vline(xintercept = 1, linetype = "dashed", color = "gray50") +
    geom_errorbarh(aes(xmin = lower95, xmax = upper95), height = 0.3) +
    geom_point(aes(color = padj < 0.1), size = 3) +
    scale_color_manual(values = c("TRUE" = "#CC79A7", "FALSE" = "#0072B2"),
                       labels = c("padj < 0.1", "padj >= 0.1")) +
    scale_x_log10() +  
    coord_cartesian(xlim = c(0.3, 3)) +
    labs(
      x = "Hazard Ratio (95% CI)",
      y = "",
      title = i,
      color = "Significance"
    ) +
    theme_classic() +
    theme(
      axis.text.y = element_text(size = 10, face = 'italic'),
      legend.position = "bottom"
    )
  
  plts[[i]] <- p1
}

length(plts)


#plts[[1]] #liver intrahepatic bile ducts, 2 pathways
#plts[[2]] #ovary 6 pathways


ggsave(filename = paste0(wd, "/Figure/RM_figures/ForestPlot_pathways_liver.png"),
       plot     = plts[[1]],
       width    = 7,
       height   = 6,
       dpi = 300)

ggsave(filename = paste0(wd, "/Figure/RM_figures/ForestPlot_pathways_ovary.png"),
       plot     = plts[[2]],
       width    = 7,
       height   = 2.5,
       dpi = 300)

#-------------------------------------------


bind_rows(sig_species_sum, .id = "cancer_from_list") %>% 
  ggplot(aes(x = HR, y = PC)) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "gray50") +
  geom_errorbarh(aes(xmin = lower95, xmax = upper95, color = cancer_from_list), 
                 height = 0.3, alpha = 1) +
  geom_point(aes(color = cancer_from_list, shape = padj < 0.1),
             size = 3, show.legend = FALSE) +
  scale_shape_manual(
    values = c("TRUE" = 16, "FALSE" = 1),
    labels = c("padj < 0.1", "padj >= 0.1"),
    name = "Significance"
  ) +
  scale_color_manual(
    values = c("#CC79A7", "#0072B2"),
    name = ""
  ) +
  scale_x_log10() + 
  labs(
    x = "Hazard Ratio (95% CI)",
    y = "",
    title = ""
  ) +
  theme_classic() +
  theme(
    axis.text.y = element_text(size = 10, face = "italic"),
    legend.text = element_text(size = 12, color = "black"),
    legend.position = "right"
  )


ggsave(paste0(wd,"/Figure/RM_figures/ForestPlot_Pathway_AllCancerFDR0.1.pdf"), width = 12, height = 3)






