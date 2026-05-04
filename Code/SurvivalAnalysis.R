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
load('Data/data.obj.raw.core.RData')



cancer.type <- names(which(table(data.obj$meta.dat$icd10_first_3_name_short)>15))
cancer.type <- cancer.type[cancer.type != 'healthy']

#first inspect general correlation between variables to inform which variables to adjust for in the survival analysis
df_sub <- data.obj$meta.dat[, covars]
df_sub$icd10_first_3_name_short <- as.character(df_sub$icd10_first_3_name_short)
df_sub$icd10_first_3_name_short[df_sub$icd10_first_3_name_short %in% names(which(table(df_sub$icd10_first_3_name_short)<15))] <- "Other"
xx <- names(which(table(df_sub$icd10_first_3_name_short)>15))
xx[!(xx %in% cancer.type)]
cancer.type[!(cancer.type %in% xx)]


# Define variable types
continuous_vars <- c("BMI", "Age", "Bristol_score", "Elix_score")
binary_vars     <- c("Sex", "Metastasis", "PPI_day_365", "Abx_day_365", "Urban")
multi_vars      <- c("Batch", "icd10_first_3_name_short", "Sample_season")

get_assoc <- function(x, y, xname, yname) {
  
  x_type <- case_when(
    xname %in% continuous_vars ~ "continuous",
    xname %in% binary_vars     ~ "binary",
    xname %in% multi_vars      ~ "multi"
  )
  y_type <- case_when(
    yname %in% continuous_vars ~ "continuous",
    yname %in% binary_vars     ~ "binary",
    yname %in% multi_vars      ~ "multi"
  )
  
  pair <- na.omit(data.frame(x = x, y = y))
  
  # continuous vs continuous
  if (x_type == "continuous" & y_type == "continuous") {
    res <- cor.test(as.numeric(pair$x), as.numeric(pair$y), method = "spearman")
    return(list(stat = res$estimate, pval = res$p.value, type = "Spearman r"))
  }
  
  # continuous vs binary 
  if ((x_type == "continuous" & y_type == "binary") |
      (x_type == "binary"     & y_type == "continuous")) {
    if (x_type == "binary") {
      cont <- as.numeric(pair$y)
      bin  <- factor(pair$x)
    } else {
      cont <- as.numeric(pair$x)
      bin  <- factor(pair$y)
    }
    w  <- wilcox.test(cont ~ bin)
    n1 <- sum(bin == levels(bin)[1])
    n2 <- sum(bin == levels(bin)[2])
    r  <- 1 - (2 * w$statistic) / (n1 * n2)
    return(list(stat = as.numeric(r), pval = w$p.value, type = "Wilcoxon"))
  }
  
  # continuous vs multi
  if ((x_type == "continuous" & y_type == "multi") |
      (x_type == "multi"      & y_type == "continuous")) {
    if (x_type == "multi") {
      cont  <- as.numeric(pair$y)
      multi <- factor(pair$x)
    } else {
      cont  <- as.numeric(pair$x)
      multi <- factor(pair$y)
    }
    kw <- kruskal.test(cont ~ multi)
    return(list(stat = as.numeric(kw$statistic), pval = kw$p.value, type = "Kruskal-Wallis"))
  }
  
  # binary vs multi 
  if ((x_type == "binary" & y_type == "multi") |
      (x_type == "multi"  & y_type == "binary") |
      (x_type == "multi" & y_type == "multi")) {
    tbl <- table(pair$x, pair$y)
    chi <- suppressWarnings(chisq.test(tbl, simulate.p.value = TRUE))
    n   <- sum(tbl)
    k   <- min(nrow(tbl), ncol(tbl))
    v   <- sqrt(chi$statistic / (n * (k - 1)))
    return(list(stat = as.numeric(v), pval = chi$p.value, type = "Cramér's V"))
  }
  
  
  # binary vs binary
  if (x_type == "binary" & y_type == "binary") {
    tbl <- table(pair$x, pair$y)
    ft  <- fisher.test(tbl)
    return(list(stat = as.numeric(ft$estimate), pval = ft$p.value, type = "Fisher's exact"))
  }
}

n <- length(covars)
pval_mat  <- matrix(NA, n, n, dimnames = list(covars, covars))



for (i in 1:n) {
  for (j in 1:n) {
    if (i == j) {
      #assoc_mat[i, j] <- 1
      pval_mat[i, j]  <- 0
    } else {
      res <- get_assoc(df_sub[[covars[i]]], df_sub[[covars[j]]], covars[i], covars[j])
      pval_mat[i, j]  <- res$pval
    }
  }
}

# use p_col as the "correlation" input, so color ~ -log10(p)
corrplot(pval_mat,
         method = "color",
         col = colorRampPalette(c("#CC79A7","white"))(200),
         diag = FALSE,
         type = "upper",
         col.lim = c(0, 1),
         # still use real p-values for significance labeling
         p.mat = pval_mat,
         sig.level = c(0.001, 0.01, 0.05),
         insig = "label_sig",
         pch.cex = 0.9,
         pch.col = "grey20",
         tl.col = "black")



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

ids <- data.obj.rff$meta.dat[is.na(data.obj.rff$meta.dat$BMI), "BMI", drop =F] 
data.obj.rff$meta.dat[rownames(ids), "BMI"] <- data.obj$meta.dat[rownames(ids), "BMI"]
obj <- cmdscale(dist.obj[["WUniFrac"]], k = 10, eig = TRUE)
pve <- round(obj$eig[1:10] / sum(abs(obj$eig)) * 100, 1);sum(pve)
eig <- obj$points
colnames(eig) <- paste0("PC",1:10)
data.obj$meta.dat <- data.obj$meta.dat %>% rownames_to_column("SampleID")%>% inner_join(as.data.frame(eig) %>% rownames_to_column("SampleID"))
data.obj$meta.dat <- data.obj$meta.dat %>% column_to_rownames('SampleID')


obj <- cmdscale(dist.obj.rff[["WUniFrac"]], k =10, eig = TRUE)
pve <- round(obj$eig[1:10] / sum(abs(obj$eig)) * 100, 1);sum(pve)
eig <- obj$points
colnames(eig) <- paste0("PC",1:10)
data.obj.rff$meta.dat <- data.obj.rff$meta.dat %>% rownames_to_column("SampleID")%>% inner_join(as.data.frame(eig) %>% rownames_to_column("SampleID"))
data.obj.rff$meta.dat <- data.obj.rff$meta.dat %>% column_to_rownames('SampleID')


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

## keep cancer only group
#remove NA in vital_status from the cancer group
samIDs <- rownames(data.obj$meta.dat[data.obj$meta.dat$Group == "Cancer" & !is.na(data.obj$meta.dat$vital_status),])
length(samIDs) #1345; so some are filtered out now

data.obj2 <- subset_data(data.obj, samIDs = samIDs)
data.obj2$meta.dat <- droplevels(data.obj2$meta.dat)
samIDs <- rownames(data.obj.rff$meta.dat[data.obj.rff$meta.dat$Group %in% "Cancer",])
data.obj.rff2 <- subset_data(data.obj.rff, samIDs = samIDs)
data.obj.rff2$meta.dat <- droplevels(data.obj.rff2$meta.dat)
dist.obj.rff2 <- subset_dist(dist.obj.rff, samIDs = samIDs)

## change small sample size cancers into other
data.obj2$meta.dat$icd10_first_3_name_short <- as.character(data.obj2$meta.dat$icd10_first_3_name_short)
data.obj2$meta.dat$icd10_first_3_name_short[!(data.obj2$meta.dat$icd10_first_3_name_short %in% cancer_15)] <- 'Other'
data.obj2$meta.dat$icd10_first_3_name_short <- as.factor(data.obj2$meta.dat$icd10_first_3_name_short)
data.obj.rff2$meta.dat$icd10_first_3_name_short <- as.character(data.obj.rff2$meta.dat$icd10_first_3_name_short)
data.obj.rff2$meta.dat$icd10_first_3_name_short[!(data.obj.rff2$meta.dat$icd10_first_3_name_short %in% cancer_15)] <- 'Other'
data.obj.rff2$meta.dat$icd10_first_3_name_short <- as.factor(data.obj.rff2$meta.dat$icd10_first_3_name_short)
table(data.obj2$meta.dat$icd10_first_3_name_short, data.obj2$meta.dat$vital_status) 

dim(data.obj2$meta.dat)
# 1345  206

# ====== marginal survival association of each variable, variable with p<0.05 were adjusted later =====
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

# ===== 1. species ===== 
## clr transform species abundance
prop <- t(t(data.obj2$abund.list$Species)/colSums(data.obj2$abund.list$Species))
idx <- (rowMeans(prop>0) > 0.1 & Biobase::rowMax(prop)> 0.002)
Y <- data.obj2$abund.list$Species[idx,]
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
                                          `Adjust covariates`=paste0(adj.name,collapse=','))
colnames(species_pvals)[1] <- 'Species'
species_pvals_pancancer <- species_pvals

species_pvals_pancancer[species_pvals_pancancer$padj < 0.1,]
#very small Hazard ratio's and not significant after multiple testing correction
species_pvals_pancancer[species_pvals_pancancer$p_value < 0.01,]
head(species_pvals_pancancer[order(species_pvals_pancancer$HR, decreasing = T),]) #max 1.05
head(species_pvals_pancancer[order(species_pvals_pancancer$HR, decreasing = F),]) #min 0.944

###======= Overall survival =======
#=====MiRKAT survival used for testing overall community with survival======
#overall unifrac dissimilarity is also not predictive of survival when adjusting for
#icd10_first_3_name_short + #Elix_score + Batch + BMI + Age + Metastasis
idx <- names(which(colSums(apply(data.obj.rff2$meta.dat[,c("OS_months", "vital_status", adj.name)], 1, function(x) !is.na(x))) == length(c("OS_months", "vital_status", adj.name))))
data.obj.rff3 <- subset_data(data.obj.rff2, samIDs = idx)
data.obj.rff3$meta.dat <- droplevels(data.obj.rff3$meta.dat)
dist.obj.rff3 <- subset_dist(dist.obj.rff2, samIDs = idx)
dim(dist.obj.rff3$WUniFrac)
nrow(data.obj.rff3$meta.dat)
model_mat <- model.matrix(as.formula(paste0('~ ', paste0(adj.name, collapse = '+'))),
                          data = data.obj.rff3$meta.dat)[, -1]
MiRKATS(obstime = data.obj.rff3$meta.dat$OS_months, delta = data.obj.rff3$meta.dat$vital_status,
        X = model_mat, Ks = dist.obj.rff3$WUniFrac, beta = NULL)
#$p_values
#[1] 0.9620924

dim(data.obj.rff3$meta.dat)
# 1338  204
dim(dist.obj.rff3$WUniFrac)
#1338 1338



#----------------------------------
#try MiRKAT on only single cancer classes
adj.name_sub <- c("Elix_score", "BMI", "Age")
p_MiRKATS_list <- list()
for (i in 1:length(unique(data.obj.rff3$meta.dat$icd10_first_3_name_short))) {

  meta_sub <- data.obj.rff3$meta.dat[which(data.obj.rff3$meta.dat$icd10_first_3_name_short == unique(data.obj.rff3$meta.dat$icd10_first_3_name_short)[i]),]
  dist_sub <- dist.obj.rff3$WUniFrac[rownames(meta_sub), rownames(meta_sub)]

  model_mat <- model.matrix(as.formula(paste0('~ ', paste0(adj.name_sub, collapse = '+'))),
                            data = meta_sub)[, -1]
  p_MiRKATS_list[[i]] <- MiRKATS(obstime = meta_sub$OS_months, delta = meta_sub$vital_status,
          X = model_mat, Ks = dist_sub, beta = NULL)
}
min(unlist(lapply(p_MiRKATS_list, unlist))) #0.1141937, minimal p value so non significant


#----------------------------------
# ===== 3. ECI score===== 
adj.name_ECI <- c("icd10_first_3_name_short", "BMI", "Age", "Metastasis") 
form_eci <- paste0('Surv(OS_months, vital_status) ~ Elix_score+',paste(adj.name_ECI, collapse = " + "))
fit_eci <- survival::coxph(as.formula(form_eci), data = data.obj2$meta.dat)
s <- summary(fit_eci)
range(data.obj2$meta.dat$Elix_score) #0-20

#Methods section
#In a multivariable Cox proportional hazards model adjusting for cancer class, BMI, age, and metastatic status, 
#each 1-point increase in Elixhauser comorbidity score (range 0–20) was associated with an 8.8% increase in the hazard 
#of death (hazard ratio 1.09, 95% confidence interval 1.06 to 1.12, p = 6.27e-08).

eci_pvals <- tibble(
  PC       = 'Elix_score',
  beta     = s$coef['Elix_score', "coef"],
  HR       = s$coef['Elix_score', "exp(coef)"],
  p_value  = s$coef['Elix_score', "Pr(>|z|)"],
  lower95  = s$conf.int['Elix_score', "lower .95"],
  upper95  = s$conf.int['Elix_score', "upper .95"],
  n        = s$n,
  events   = s$nevent
) %>% mutate(qval = p.adjust(p_value,'BH'))

combined_data <- data.obj2$meta.dat
combined_data$group <- ifelse(
  combined_data[["Elix_score"]] > median(combined_data[["Elix_score"]]),
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
  legend.title = "Elix_score",
  xlim = c(0, 60), #an xlim of 60 is appropriate as beyond that there are few patients left
  break.time.by = 10,
  palette = c("#CC79A7", "#0072B2"),
  ggtheme = theme_minimal(base_size = 12) +
    theme(
      axis.title.x = element_text(size = 12, colour = "black"),
      axis.title.y = element_text(size = 12, colour = "black")
    )
)

combined <- plot_grid(
  p$plot,
  p$table,
  ncol        = 1,
  rel_heights = c(2, 1)   # adjust relative heights as needed
)

ggsave(filename = paste0(wd, "/Figure/RM_figures/survival_curve_Elix_score.png"),
  plot     = combined,
  width    = 4,
  height   = 6,
  dpi = 300
)



# ===== 4. each ECI component===== 
elix.names <- c("Elixhauser_CHF","Elixhauser_Arrhythmia","Elixhauser_Valvular","Elixhauser_PHTN","Elixhauser_PVD","Elixhauser_HTN","Elixhauser_Paralysis",  
                "Elixhauser_NeuroOther","Elixhauser_Pulmonary","Elixhauser_DM","Elixhauser_DMcx","Elixhauser_Hypothyroid","Elixhauser_Renal","Elixhauser_Liver",      
                "Elixhauser_PUD" ,"Elixhauser_HIV","Elixhauser_Psychoses", "Elixhauser_Mets",
                "Elixhauser_Lymphoma","Elixhauser_Mets","Elixhauser_Tumor","Elixhauser_Rheumatic","Elixhauser_Coagulopathy",
                "Elixhauser_Obesity","Elixhauser_WeightLoss","Elixhauser_FluidsLytes","Elixhauser_BloodLoss","Elixhauser_Anemia","Elixhauser_Alcohol","Elixhauser_Drugs",
                "Elixhauser_Depression")

elix.names <- elix.names[!elix.names %in% c("Elixhauser_Tumor", "Elixhauser_Lymphoma", "Elixhauser_Obesity",
                                                        "Elixhauser_Paralysis", "Elixhauser_HIV",
                                                        "Elixhauser_BloodLoss", "Elixhauser_Psychoses")]

bad_comp <- c()
for (comp in elix.names) {
  tab <- table(data.obj2$meta.dat[[comp]], data.obj2$meta.dat$vital_status)
  if (nrow(tab) < 2) {
    bad_comp <- c(bad_comp, comp)
    next
  }

  if (any(tab[, "1"] <= 1) || any(tab[, "0"] <= 1)) {
    bad_comp <- c(bad_comp, comp)
  }
}


elix.names <- elix.names[!(elix.names %in% bad_comp)]
eci_comp_results <- list()
adj.name_sub_list <- list()
for(comp in elix.names) {
  adj.name_sub <- adj.name_ECI[!adj.name_ECI == "Elix_score"]
  if(comp %in% "Elixhauser_Mets") adj.name_sub <- adj.name_sub[!(adj.name_sub %in% "Metastasis")]
  if(comp %in% c('Elixhauser_Obesity','Elixhauser_WeightLoss')) adj.name_sub <- adj.name_sub[!(adj.name_sub %in% "BMI")]
  form <- as.formula(paste(
    "Surv(OS_months, vital_status) ~", 
    paste(c(comp, adj.name_sub), collapse = " + ")
  ))
  
  fit <- coxph(form, data = data.obj2$meta.dat)
  eci_comp_results[[comp]] <- summary(fit)
  adj.name_sub_list[[comp]] <- paste0(adj.name_sub,collapse = ',')
}


eci_comp_pvals <- lapply(names(eci_comp_results), function(pc) {
  s <- eci_comp_results[[pc]]
  tibble(
    PC       = pc,
    beta     = s$coef[paste0(pc,"Yes"), "coef"],
    HR       = s$coef[paste0(pc,"Yes"), "exp(coef)"],
    p_value  = s$coef[paste0(pc,"Yes"), "Pr(>|z|)"],
    lower95  = s$conf.int[paste0(pc,"Yes"), "lower .95"],
    upper95  = s$conf.int[paste0(pc,"Yes"), "upper .95"],
    n        = s$n,
    events   = s$nevent
  )
}) %>% bind_rows() %>% mutate(qval  = p.adjust(p_value,"BH"))

adj_eci <- data.frame(
  PC = names(adj.name_sub_list),
  covariates  = unlist(adj.name_sub_list),
  row.names   = NULL
)

eci_comp_pvals <- left_join(as.data.frame(eci_comp_pvals), adj_eci)
  


eci_comp_pvals[grep("Weight", eci_comp_pvals$PC),] #removed ECI score adjustment made this significant
eci_comp_results["Elixhauser_WeightLoss"]

sig_eci <- eci_comp_pvals  %>%
  #dplyr::filter(qval < 0.1) %>%
  arrange(HR) %>%
  mutate(PC = factor(PC, levels = PC))

sig_eci$PC_name <- unname(names_map[sig_eci$PC]) ### Fix this!! 
sig_eci <- sig_eci %>%
  mutate(PC_name = factor(PC_name, levels = PC_name))
paste0(as.character(sig_eci$PC_name[sig_eci$qval<0.1]), collapse = ', ')
#make Forest plot for ECI components
p1 <- ggplot(sig_eci, aes(x = HR, y = PC)) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "gray50") +
  geom_errorbarh(aes(xmin = lower95, xmax = upper95), height = 0.3) +
  geom_point(aes(color = qval < 0.1), size = 3) +
  scale_color_manual(values = c("TRUE" = "#CC79A7", "FALSE" = "black"), 
                     labels = c("q < 0.1", "q ≥ 0.1")) +
  scale_x_log10() +
  labs(
    x = "Hazard Ratio (95% CI)",
    y = "",
    title = "Elixhauser Components and Overall Survival",
    color = "Significance"
  ) +
  theme_classic() +
  coord_cartesian(xlim = c(0.3, 3.5)) +
  theme(
    axis.text.y = element_text(size = 10),
    legend.position = "bottom"
  )


ggsave(filename = paste0(wd, "/Figure/RM_figures/ForestPlot_Elixhauser.png"),
       plot     = p1,
       width    = 7,
       height   = 6,
       dpi = 300
)


##===== summary view all significance =====
pvals <- bind_rows(eci_pvals %>% mutate(`Adjust covariates`=paste0(adj.name_sub,collapse=',')), 
                   eci_comp_pvals )
colnames(pvals)[1] <- 'Variable'

pvals[pvals$p_value<0.05,]
pvals[pvals$qval<0.1,]
species_pvals_pancancer[species_pvals_pancancer$padj<0.1,]
pvals[pvals$qval<0.1,] %>% dim


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
  prop <- t(t(data.obj2$abund.list$Species)/colSums(data.obj2$abund.list$Species))
  idx <- (rowMeans(prop>0) > 0.1 & Biobase::rowMax(prop)> 0.002)
  Y <- data.obj2$abund.list$Species[idx,]
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
    
    pdf(paste0(wd,"/Figure/RM_figures/survival_curve_species_FDR0.1_",cancer_30[i],".pdf"), width = 6, height = 6)
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
        legend.title = gsub('^s__','',species_name),
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


colnames(species_pvals_all)[1] <- 'Species'
species_pvals_all <- species_pvals_all %>% mutate(`Adjust covariates`=paste0(adj.name, collapse = ','))


## prepare Supplementary Table for Survival Species
library(openxlsx)
setwd(wd)
wb <- loadWorkbook("Code/Submission/Supplementary tables/Supplementary Table 9.xlsx")
wb <- loadWorkbook("Code/Submission/Supplementary tables/20260424 revision Supplement/Table S8_new.xlsx")
names(wb)

removeWorksheet(wb, "TableS8-1")
removeWorksheet(wb, "TableS9-2")
removeWorksheet(wb, "TableS9-3")
removeWorksheet(wb, "TableS9-4")
removeWorksheet(wb, "TableS9-5")

addWorksheet(wb, "TableS8-1")
writeData(wb, sheet = "TableS8-1", x = pvals)

addWorksheet(wb, "TableS8-2")
writeData(wb, sheet = "TableS8-2", x = species_pvals_pancancer)

addWorksheet(wb, "TableS8-4")
writeData(wb, sheet = "TableS8-4", x = species_pvals_all)
saveWorkbook(wb, "Code/Submission/Supplementary tables/20260424 revision Supplement/Table S8_new.xlsx", overwrite = TRUE)



## Make forest plot
res_df <- bind_rows(res, .id = "cancer_class")
res_df <- res_df %>%
  dplyr::filter(padj < 0.1) 
names(which(table(res_df$cancer_class)>1))

plts <- sig_species_sum <- list()
for (i in names(which(table(res_df$cancer_class) > 0))) {
  sig_species <- res_df[res_df$cancer_class == i, ]
  sig_species$PC <- gsub('^s__','',sig_species$PC)
  
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
      title = gsub('ducts','duct',i),
      color = "Significance"
    ) +
    theme_classic() +
    theme(
      axis.text.y = element_text(size = 10, face = 'italic'),
      legend.position = "bottom"
    )
  
  plts[[i]] <- p1
  #pdf(paste0(wd,"/Figure/RM_figures/ForestPlot_Species_", i, ".pdf"), width = 7, height = 5)
  #print(p1)
  #dev.off()
}

length(plts)


#plts[[1]] #CRC, 1 species
#plts[[2]] #liver many species
#plts[[3]] #melanoma 5 species
#plts[[4]] #prostate 3 species


ggsave(filename = paste0(wd, "/Figure/RM_figures/ForestPlot_species_CRC.png"),
       plot     = plts[[1]],
       width    = 7,
       height   = 1.75,
       dpi = 300)

ggsave(filename = paste0(wd, "/Figure/RM_figures/ForestPlot_species_liver.png"),
       plot     = plts[[2]],
       width    = 7,
       height   = 6,
       dpi = 300)

ggsave(filename = paste0(wd, "/Figure/RM_figures/ForestPlot_species_melanoma.png"),
       plot     = plts[[3]],
       width    = 7,
       height   = 2.5,
       dpi = 300)

ggsave(filename = paste0(wd, "/Figure/RM_figures/ForestPlot_species_prostate.png"),
       plot     = plts[[4]],
       width    = 7,
       height   = 2,
       dpi = 300)


#-------------------------------------------

ggarrange(plts[[2]], ggarrange(plts[[1]], plts[[4]], plts[[3]], nrow = 3, common.legend = T), common.legend = T)


bind_rows(sig_species_sum, .id = "cancer_from_list") %>% 
  ggplot(aes(x = HR, y = PC)) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "gray50") +
  geom_errorbarh(aes(xmin = lower95, xmax = upper95, color = cancer_from_list), 
                 height = 0.3, alpha = 1) +
  geom_point(aes(color = cancer_from_list, shape = padj < 0.1), size = 3, show.legend = F) +
  scale_shape_manual(values = c("TRUE" = 16, "FALSE" = 1),
                     labels = c("padj < 0.1", "padj >= 0.1"),
                     name = "Significance") +
  scale_color_discrete(palette = "Set1",name = "") +
  scale_x_log10() + 
  labs(
    x = "Hazard Ratio (95% CI)",
    y = "",
    title = ""
  ) +
  theme_classic() +
  theme(
    axis.text.y = element_text(size = 10, face = 'italic'),
    legend.text = element_text(size = 12, color = 'black'),
    legend.position = "right"
  )
ggsave(paste0(wd,"/Figure/RM_figures/ForestPlot_Species_AllCancerFDR0.1.pdf"), width = 10, height = 6)




# ====== whether the number of SIG1 species are associated to OS?======
file_dir <- dirname(rstudioapi::getSourceEditorContext()$path)
file_dir <- dirname(dirname(dirname(dirname(file_dir))))


wd <- file_dir
load('Data/data.obj.mph.RData')

samIDs <- rownames(data.obj$meta.dat[data.obj$meta.dat$Group == "Cancer" & !is.na(data.obj$meta.dat$vital_status),])
length(samIDs) #1345; so some are filtered out now

data.obj2 <- subset_data(data.obj, samIDs = samIDs)
data.obj2$meta.dat <- droplevels(data.obj2$meta.dat)
samIDs <- rownames(data.obj.rff$meta.dat[data.obj.rff$meta.dat$Group %in% "Cancer",])

## change small sample size(n<15) cancers into other
data.obj2$meta.dat$icd10_first_3_name_short <- as.character(data.obj2$meta.dat$icd10_first_3_name_short)
data.obj2$meta.dat$icd10_first_3_name_short[!(data.obj2$meta.dat$icd10_first_3_name_short %in% cancer_15)] <- 'Other'
data.obj2$meta.dat$icd10_first_3_name_short <- as.factor(data.obj2$meta.dat$icd10_first_3_name_short)
table(data.obj2$meta.dat$icd10_first_3_name_short, data.obj2$meta.dat$vital_status) 

dim(data.obj2$meta.dat)
# 1345  206

fit <- coxph(Surv(OS_months, vital_status) ~SIG1_count+ Age + BMI + Elix_score + Metastasis +icd10_first_3_name_short, data = data.obj2$meta.dat)
summary(fit)

data.obj2$meta.dat$SIG1_group <- ifelse(
  data.obj2$meta.dat$SIG1_count >= median(data.obj2$meta.dat$SIG1_count, na.rm = TRUE),
  "High", "Low"
)
fit <- survfit(Surv(OS_months, vital_status) ~ SIG1_group, data = data.obj2$meta.dat)
ggsurvplot(
  fit,
  data = data.obj2$meta.dat,
  pval = F,
  conf.int = TRUE,
  risk.table = F,
  risk.table.height = 0.25,
  xlab = "Time (months)",
  ylab = "Overall Survival Probability",
  legend.title = "SIG1_count",
  palette = c("#0072B2", "#CC79A7"),
  ggtheme = theme_minimal()
)
ggsave(paste0(wd,"/Figure/RM_figures/KM_SIG1.pdf"), width = 6, height = 4)


ggsurvplot(
  fit,
  data = data.obj2$meta.dat,
  pval = F,
  conf.int = TRUE,
  risk.table = T,
  risk.table.height = 0.25,
  xlab = "Time (months)",
  ylab = "Overall Survival Probability",
  legend.title = "SIG1_count",
  palette = c("#0072B2", "#CC79A7"),
  ggtheme = theme_minimal()
)


#--------------------------------------
#toposcore as a group associated to survival?
fit <- coxph(
  Surv(OS_months, vital_status) ~ 
    Toposcore + Age + BMI + Elix_score + Metastasis + icd10_first_3_name_short,
  data = data.obj2$meta.dat
)
summary(fit)

# fit <- coxph(
#   Surv(OS_months, vital_status) ~ 
#     Toposcore + Age + BMI + Elix_score + Metastasis +
#     strata(icd10_first_3_name_short),
#   data = data.obj2$meta.dat
# )
# summary(fit)
#                 coef exp(coef)  se(coef)      z Pr(>|z|)    
#ToposcoreSIG2+ -0.379138  0.684451  0.118906 -3.189  0.00143 ** 
#Age             0.020823  1.021041  0.004769  4.366 1.26e-05 ***
#BMI            -0.022116  0.978127  0.008818 -2.508  0.01214 *  
#Elix_score      0.063251  1.065295  0.015795  4.004 6.22e-05 ***
#MetastasisYes   1.081534  2.949201  0.119307  9.065  < 2e-16 ***

exp(cbind(HR = coef(fit), confint(fit)))
#                   HR        2.5 %    97.5 %
#ToposcoreSIG2+ 0.6844513 0.5421624 0.8640836 #being SIG2+ is protective
#Age            1.0210413 1.0115415 1.0306304
#BMI            0.9781268 0.9613665 0.9951793
#Elix_score     1.0652945 1.0328200 1.0987901
#MetastasisYes  2.9492010 2.3342649 3.7261351




fit <- survfit(Surv(OS_months, vital_status) ~ Toposcore, data = data.obj2$meta.dat)
summary(fit)

p <- ggsurvplot(
  fit,
  data = data.obj2$meta.dat,
  pval = F,
  conf.int = TRUE,
  risk.table = T,
  risk.table.height = 0.2,
  title = "",
  xlab = "Time (months)",
  ylab = "Overall Survival Probability",
  legend.title = "Toposcore",
  xlim = c(0, 60), #an xlim of 60 is appropriate as beyond that there are few patients left
  break.time.by = 10,
  palette = c("#CC79A7", "#0072B2"),
  ggtheme = theme_minimal(base_size = 12) +
    theme(
      axis.title.x = element_text(size = 12, colour = "black"),
      axis.title.y = element_text(size = 12, colour = "black")
    )
)

combined <- plot_grid(
  p$plot,
  p$table,
  ncol        = 1,
  rel_heights = c(3, 1)   # adjust relative heights as needed
)

print(combined)

ggsave(filename = paste0(wd, "/Figure/RM_figures/survival_curve_Toposcore.png"),
       plot     = combined,
       width    = 5,
       height   = 7,
       dpi = 300
)


#--------------------------------------
#now remove the immunotherapy patients and rerun to see if this is also predictive of survival with other treatments
#yes. Similarly protective
head(data.obj2$meta.dat)

meta.dat_sub <- data.obj2$meta.dat[grep("Immune.checkpoint.inhibitor", data.obj2$meta.dat$category_med, invert = T),]
dim(meta.dat_sub) #1163  202
table(meta.dat_sub$Toposcore)
#SIG1+ SIG2+ 
# 854   309 

fit <- coxph(
  Surv(OS_months, vital_status) ~ 
    Toposcore + Age + BMI + Elix_score + Metastasis + icd10_first_3_name_short,
  data = meta.dat_sub
)
summary(fit)

#                                                         coef exp(coef)  se(coef)      z Pr(>|z|)    
#ToposcoreSIG2+                                        -0.430420  0.650236  0.134200 -3.207  0.00134 ** 
#Age                                                    0.018467  1.018639  0.005236  3.527  0.00042 ***
#BMI                                                   -0.029707  0.970730  0.010150 -2.927  0.00343 ** 
#Elix_score                                             0.087220  1.091137  0.017379  5.019  5.2e-07 ***
#MetastasisYes                                          1.123260  3.074862  0.131043  8.572  < 2e-16 ***
#icd10_first_3_name_shortbrain                          1.167209  3.213012  0.441100  2.646  0.00814 ** 
#icd10_first_3_name_shortbreast                        -0.636562  0.529108  0.365389 -1.742  0.08148 .  
#icd10_first_3_name_shortbronchus lung                 -0.001561  0.998440  0.370125 -0.004  0.99663    
#icd10_first_3_name_shortcolorectal                     0.201003  1.222628  0.376336  0.534  0.59327    
#icd10_first_3_name_shortcorpus uteri                   0.430893  1.538631  0.428982  1.004  0.31516    
#icd10_first_3_name_shortesophagus                      0.438974  1.551115  0.414552  1.059  0.28964    
#icd10_first_3_name_shortkidney except renal pelvis     0.663924  1.942398  0.500499  1.327  0.18467    
#icd10_first_3_name_shortliver intrahepatic bile ducts  0.839042  2.314149  0.390490  2.149  0.03166 *  
#icd10_first_3_name_shortlymphoid leukemia              0.091432  1.095743  0.664623  0.138  0.89058    
#icd10_first_3_name_shortmelanoma of skin              -0.391977  0.675720  0.778722 -0.503  0.61471    
#icd10_first_3_name_shortmultiple myeloma plasma cell   0.255616  1.291256  0.407937  0.627  0.53092    
#icd10_first_3_name_shortneuroendocrine tumors         -0.032510  0.968013  0.494915 -0.066  0.94763    
#icd10_first_3_name_shortoropharynx                    -0.832713  0.434868  0.594306 -1.401  0.16117    
#icd10_first_3_name_shortOther                         -0.311780  0.732142  0.342253 -0.911  0.36231    
#icd10_first_3_name_shortother connective soft tissue   0.474030  1.606455  0.408023  1.162  0.24533    
#icd10_first_3_name_shortovary                          0.277182  1.319407  0.386290  0.718  0.47304    
#icd10_first_3_name_shortpancreas                       0.790250  2.203948  0.365550  2.162  0.03063 *  
#icd10_first_3_name_shortprostate                      -0.312027  0.731962  0.375472 -0.831  0.40596    
#icd10_first_3_name_shortunspecified skin              -0.458966  0.631937  0.549151 -0.836  0.40328    


