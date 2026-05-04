library(readr)
library(tidyverse)
library(readr)
library(ape)
library(stringr)
library(magrittr)
library(phyloseq)
library(tidyr)
library(dplyr)
library(knitr)
library(kableExtra)
library(gridExtra)
library(grid)
library(openxlsx)


file_dir <- dirname(rstudioapi::getSourceEditorContext()$path)
file_dir <- dirname(dirname(dirname(dirname(dirname(file_dir)))))

wd <- file_dir
rd <- paste0(file_dir,'/Result/')

if(!dir.exists(rd)) dir.create(rd)
setwd(wd)
source(paste0(file_dir,"/Code/Submission/MayoOncobiomeStudy/Code/Stats.R"))
try(load_package())

load(file = 'Data/data.obj.mph.RData') 
meta_with_topo <- data.obj$meta.dat

load(file = 'Data/data.obj.raw.core.RData') 

#also load the Kraken tax table
tax_table_Kraken <- as.data.frame(data.obj.rff$otu.tab)
tax_table_RA_Kraken <- sweep(tax_table_Kraken, MARGIN = 2, colSums(tax_table_Kraken), '/') * 100


#The data.obj$meta.dat now includes the TOPOSCORE for each sample. Relevant columns include:
#  "SIG1_count", "SIG2_count", "S_score", "SIG_class", and "Toposcore".

meta_with_topo$Cancer_Healthy <- "Cancer"
meta_with_topo$Cancer_Healthy[is.na(meta_with_topo$cancer_class)] <- "Healthy"

table(meta_with_topo$Cancer_Healthy )
#Cancer Healthy 
#1364     287 

meta_with_topo$icd10_first_3_name_short_cutoff <- as.character(meta_with_topo$icd10_first_3_name_short)
below_15_names <- names(which(table(meta_with_topo$icd10_first_3_name_short_cutoff) < 15))
meta_with_topo$icd10_first_3_name_short_cutoff[meta_with_topo$icd10_first_3_name_short_cutoff %in% below_15_names] <- "other"
table(meta_with_topo$icd10_first_3_name_short_cutoff)


#write Toposcore columns to file for Table S1-7
head(meta_with_topo)

#write.csv(meta_with_topo[, c("SIG1_count", "SIG2_count", "S_score", "SIG_class", "Toposcore")], file = "./Code/Submission/Supplementary tables/toposcore_temp.csv")


#-------------------------------------------------------------------------------
#load Metaphlan taxonomy species level

species_metaphlan <- as.data.frame(read.csv('~/Dropbox/Mayo_RS/R/Oncobiome full cohort analysis/data/Metaphlan_taxonomy_MCCM_species.csv'))
rownames(species_metaphlan) <- species_metaphlan$X

species_metaphlan <- species_metaphlan[,rownames(meta_with_topo)]

species_metaphlan[1:4, 1:4]

meta_with_topo$s__Akkermansia_muciniphila <- as.numeric(species_metaphlan[grep("s__Akkermansia_muciniphila", rownames(species_metaphlan)),])

hist(meta_with_topo$s__Akkermansia_muciniphila)
table(meta_with_topo$s__Akkermansia_muciniphila > 0)
#FALSE  TRUE 
# 926   725 

meta_with_topo$Akk_detected <- if_else(meta_with_topo$s__Akkermansia_muciniphila > 0, "Akk_detected", "Akk_not_detected")


#-------------------------------------------------------------------------------
#from co‑abundance networks in NSCLC patients on ICIs
#SIG1 (harmful, associated with resistance/non‑response) 
#SIG2 (beneficial, associated with response) 

#S score is computed as the number of SIG2 metagenomic species (MGS) detected above a predefined presence/abundance threshold in that fecal sample.

#not a huge number of species but fine.
hist(meta_with_topo$SIG1_count)
hist(meta_with_topo$SIG2_count)


#improve labels and colors
plot(density(meta_with_topo$S_score[meta_with_topo$Cancer_Healthy == "Healthy"]), las=1, xlab="",
     main = "S score distribution\n(higher equals more health-associated species)", col="#0072B2", lwd=2)
lines(density(meta_with_topo$S_score[meta_with_topo$Cancer_Healthy == "Cancer"]), col="#CC79A7", lwd=2)


#add extra split on Akkermansia
density_h_a <- density(meta_with_topo$S_score[meta_with_topo$Cancer_Healthy == "Healthy" & meta_with_topo$Akk_detected == "Akk_detected"])
density_h_no_a <- density(meta_with_topo$S_score[meta_with_topo$Cancer_Healthy == "Healthy" & meta_with_topo$Akk_detected == "Akk_not_detected"])
density_c_a <- density(meta_with_topo$S_score[meta_with_topo$Cancer_Healthy == "Cancer" & meta_with_topo$Akk_detected == "Akk_detected"])
density_c_no_a <- density(meta_with_topo$S_score[meta_with_topo$Cancer_Healthy == "Cancer" & meta_with_topo$Akk_detected == "Akk_not_detected"])


plot(density_h_a, las=1, xlab="S-score distribution", main = "", col="#0072B2", lwd=2, xlim=c(0.1,1))
lines(density_h_no_a, col="#0072B2", lwd=2, lty=2)
lines(density_c_a, col="#CC79A7", lwd=2)
lines(density_c_no_a, col="#CC79A7", lwd=2, lty=2)

mtext(text=expression("more health-associated species" %->% ""), adj=1)

#155, 132, 507, 794, respectively
legend("topleft",
       legend = c("Healthy, Akk+;  n=155", "Healthy, Akk−;  n=132",
                  "Cancer, Akk+;  n=507",  "Cancer, Akk−;  n=794"),
       col    = c("#0072B2", "#0072B2", "#CC79A7", "#CC79A7"),
       lty    = c(1, 2, 1, 2),
       lwd    = 2,
       bty    = "n", cex = 0.9)  # no box, optional


#------------------------------------------
#Use univariable logistic regression to see how predictive meta_with_topo$S_score is for Cancer vs Healthy


## Make sure outcome is a factor with desired reference
meta_with_topo$Cancer_Healthy <- factor(meta_with_topo$Cancer_Healthy)
meta_with_topo$Cancer_Healthy <- relevel(meta_with_topo$Cancer_Healthy, ref = "Healthy")

## (Optional) minimal data frame
df <- meta_with_topo[, c("S_score", "Cancer_Healthy")]
df$Cancer_Healthy <- factor(df$Cancer_Healthy)
df$Cancer_Healthy <- relevel(df$Cancer_Healthy, ref = "Healthy")

## 1) Fit univariable logistic regression: Cancer vs S_score
##  (logit model via glm, family = binomial)
glm_fit <- glm(Cancer_Healthy ~ S_score, data = df, family = binomial)

summary(glm_fit)  # coefficients, p-value for S_score, etc. 

## 2) Get predicted probabilities for Cancer
glm_prob  <- predict(glm_fit, type = "response")   # P(Cancer | S_score)
glm_class <- ifelse(glm_prob >= 0.5, "Cancer", "Healthy")
glm_class <- factor(glm_class, levels = levels(df$Cancer_Healthy))

## 3) Confusion matrix and simple performance metrics in base R
tab <- table(Predicted = glm_class, True = df$Cancer_Healthy)
tab
#True
#Predicted Healthy Cancer
#Healthy      55     38
#Cancer      232   1326


accuracy    <- sum(diag(tab)) / sum(tab)
sensitivity <- tab["Cancer", "Cancer"] / sum(tab[, "Cancer"])      # TPR
specificity <- tab["Healthy", "Healthy"] / sum(tab[, "Healthy"])   # TNR

perf_tab <- data.frame(
  Metric = c("Accuracy", "Sensitivity", "Specificity"),
  Value  = c(accuracy, sensitivity, specificity)
)
perf_tab


## 4) ROC curve and AUC using pROC 
library(pROC)

roc_obj <- roc(df$Cancer_Healthy, glm_prob, levels = c("Healthy", "Cancer"))
auc_val <- auc(roc_obj)

#looks square when plotting area is
plot(roc_obj,
     col = "#0072B2", lwd = 2,
     xlim = c(1, 0),
     ylim = c(0, 1),
     main = paste0("ROC for S score logistic model (AUC = ", round(auc_val, 3), ")"))
abline(a = 0, b = 1, lty = 2, col = "grey60")




##  now look at specificity at target of 80% sensitivity ##
## 1) Existing objects from before:

## 2) Inspect all ROC coordinates
all_coords <- coords(roc_obj, x = "all", ret = c("threshold", "sensitivity", "specificity"),
                     transpose = FALSE)  # data.frame 

## 3) Filter to thresholds with sensitivity around your target (e.g. >= 0.8)
target_sens  <- 0.80
tol          <- 0.05  # +/- tolerance around target sensitivity

cand <- subset(all_coords,
               sensitivity >= (target_sens - tol) &
                 sensitivity <= 1.0)

## Among these candidates, pick the one with max specificity
best_idx      <- which.max(cand$specificity)
best_threshold <- cand$threshold[best_idx]

#this is the single sample with the best balance
cand[best_idx, ]  # inspect chosen sens/spec pair
#     threshold sensitivity specificity
#150  0.793208   0.7529326   0.7212544


## 4) Classify using this new threshold
glm_class_opt <- ifelse(glm_prob >= best_threshold, "Cancer", "Healthy")
glm_class_opt <- factor(glm_class_opt, levels = levels(df$Cancer_Healthy))

tab_opt <- table(Predicted = glm_class_opt, True = df$Cancer_Healthy)
tab_opt
#Looks a lot better
#Predicted Healthy Cancer
#Healthy     207    337 
#Cancer       80   1027

prop.test(t(tab_opt))
#prop 1    prop 2   #Percentage of healthy predicted correctly is 72.1% and cancer predicted correctly is 75.3%
#0.7212544 0.2470674 


accuracy_opt    <- sum(diag(tab_opt)) / sum(tab_opt)
sensitivity_opt <- tab_opt["Cancer", "Cancer"] / sum(tab_opt[, "Cancer"])
specificity_opt <- tab_opt["Healthy", "Healthy"] / sum(tab_opt[, "Healthy"])

perf_tab_opt <- data.frame(
  Metric = c("Accuracy", "Sensitivity", "Specificity", "Threshold"),
  Value  = c(accuracy_opt, sensitivity_opt, specificity_opt, best_threshold)
)
perf_tab_opt
#       Metric     Value
#1    Accuracy 0.7474258
#2 Sensitivity 0.7529326
#3 Specificity 0.7212544
#4   Threshold 0.7932080


##Calibration-style plot in base R
##    Observed Cancer (0/1) vs predicted probability
y_num <- as.numeric(df$Cancer_Healthy == "Cancer")

## Recompute predicted class using your chosen threshold
glm_class_opt <- ifelse(glm_prob >= best_threshold, 1, 0)  # 1 = Cancer, 0 = Healthy

## Calibration plot with vertical threshold line
plot(glm_prob, y_num,
     xlab = "Predicted P(Cancer) from S score logistic model",
     ylab = "Observed Cancer (0/1)",
     main = paste0("Calibration of S score predictions\n(threshold = ",
                   round(best_threshold, 2), ")"),
     pch = 16,
     col = ifelse(glm_class_opt == 1,
                  rgb(0.8, 0, 0.5, 0.4),   # predicted Cancer
                  rgb(0, 0, 0, 0.25)))     # predicted Healthy

lines(lowess(glm_prob, y_num), col = "#CC79A7", lwd = 2)

## Show chosen decision threshold
abline(v = best_threshold, col = "dodgerblue3", lwd = 2, lty = 2)



#------------------------------------------
#look at proportion of TOPOSCORE across cancer vs healthy and cancer classes

table(meta_with_topo$Cancer_Healthy, meta_with_topo$Toposcore)
#       Cancer Healthy
#SIG1+   1000     112 #SIG1 = unhealthy associated 
#SIG2+    364     175 #SIG2 = healthy associated
prop.test(table(meta_with_topo$Cancer_Healthy, meta_with_topo$Toposcore))
#2-sample test for equality of proportions with continuity correction
#data:  table(meta_with_topo$Cancer_Healthy, meta_with_topo$Toposcore)
#X-squared = 125.23, df = 1, p-value < 2.2e-16
#alternative hypothesis: two.sided
#95 percent confidence interval:
#  0.2796625 0.4061253
#sample estimates:
#  prop 1    prop 2 
#0.7331378 0.3902439  #73.3% SIG1+ (bad) in Cancer and 39.0% in Healthy


#3 categories using SIG_class
tab <- table(meta_with_topo$Cancer_Healthy, meta_with_topo$SIG_class)
tab
#         Gray SIG1 SIG2
#Healthy  164   19  104
#Cancer   713  557   94
chisq.test(tab) #overall proportion p-value < 2.2e-16

prop.table(tab, margin = 1) * 100
#           Gray      SIG1      SIG2
#Healthy 57.142857  6.620209 36.236934
#Cancer  52.272727 40.835777  6.891496

#---------------------

tab_ch <- table(meta_with_topo$Toposcore,
                meta_with_topo$Cancer_Healthy)

col_pct <- function(x) prop.table(x, margin = 2) * 100
tab_ch_p <- col_pct(tab_ch)


fmt_cp <- function(count, pct, digits = 1) {
  paste0(count, " (", round(pct, digits), "%)")
}

make_cp_table <- function(count_tab, pct_tab, digits = 1) {
  out <- matrix(
    nrow = nrow(count_tab),
    ncol = ncol(count_tab),
    dimnames = dimnames(count_tab)
  )
  for (i in seq_len(nrow(count_tab))) {
    for (j in seq_len(ncol(count_tab))) {
      out[i, j] <- fmt_cp(count_tab[i, j], pct_tab[i, j], digits)
    }
  }
  out
}

tab_ch_cp <- make_cp_table(tab_ch, tab_ch_p)


kbl(tab_ch_cp,
    caption = "Toposcore vs cancer status (column %)",
    align = "c") %>%
  kable_classic(
    full_width = FALSE
  ) %>%
  add_header_above(c(
    " " = 1,
    "Cancer_Healthy" = ncol(tab_ch_cp)
  )) %>%
  footnote(
    general = c(
      "SIG1+ = unhealthy-associated toposcore; SIG2+ = healthy-associated toposcore."
    ),
    general_title = ""
  )



#---------------------

#inspect in context of Akkermansia
table(meta_with_topo$Toposcore, meta_with_topo$Akk_detected, meta_with_topo$Cancer_Healthy)
#, ,  = Healthy
#           Akk_detected   Akk_not_detected
#SIG1+           27               85
#SIG2+          128               47

#, ,  = Cancer
#           Akk_detected Akk_not_detected
#SIG1+          239              761
#SIG2+          331               33


#explore making the table look better


tab <- table(meta_with_topo$Toposcore,
             meta_with_topo$Akk_detected,
             meta_with_topo$Cancer_Healthy)

tab_H <- tab[, , "Healthy"]
tab_C <- tab[, , "Cancer"]

row_pct <- function(x) {
  prop.table(x, margin = 1) * 100
}

tab_H_p <- row_pct(tab_H)
tab_C_p <- row_pct(tab_C)



fmt_cp <- function(count, pct, digits = 1) {
  paste0(count, " (", round(pct, digits), "%)")
}

make_cp_table <- function(count_tab, pct_tab, digits = 1) {
  out <- matrix(
    nrow = nrow(count_tab),
    ncol = ncol(count_tab),
    dimnames = dimnames(count_tab)
  )
  for (i in seq_len(nrow(count_tab))) {
    for (j in seq_len(ncol(count_tab))) {
      out[i, j] <- fmt_cp(count_tab[i, j], pct_tab[i, j], digits)
    }
  }
  out
}

tab_H_cp <- make_cp_table(tab_H, tab_H_p)
tab_C_cp <- make_cp_table(tab_C, tab_C_p)


kbl(tab_H_cp,
    caption = "Toposcore vs Akkermansia detection in Healthy (row %)",
    align = "c") %>%
  kable_classic(
    full_width = FALSE
  ) %>%
  add_header_above(c(" " = 1, "Akk_detected" = ncol(tab_H_cp)))


kbl(tab_C_cp,
    caption = "Toposcore vs Akkermansia detection in Cancer (row %)",
    align = "c") %>%
  kable_classic(
    full_width = FALSE
  ) %>%
  add_header_above(c(" " = 1, "Akk_detected" = ncol(tab_C_cp)))




#---------------------


#look at Akk abundance across TOPOSCORES; some SIG1 (unhealthy) have very high Akk which is expected to see with abx 
boxplot(meta_with_topo$s__Akkermansia_muciniphila ~ meta_with_topo$Toposcore)



#now TOPOSCORE proportions by cancer class
tab <- table(meta_with_topo$icd10_first_3_name_short_cutoff,
             meta_with_topo$Toposcore)
row_pct <- prop.table(tab, margin = 1) * 100
col_pct <- prop.table(tab, margin = 2) * 100
overall_pct <- prop.table(tab) * 100


combined <- paste0(tab, " (", sprintf("%.1f", row_pct), "%)")
combined <- matrix(
  combined,
  nrow = nrow(tab),
  ncol = ncol(tab),
  dimnames = dimnames(tab)
)

combined <- combined[order(row_pct[, "SIG2+"], decreasing = T),]
combined
#                               SIG1+         SIG2+        
#healthy                       "112 (39.0%)" "175 (61.0%)"
#other ill-defined sites       "8 (53.3%)"   "7 (46.7%)"  
#base of tongue                "13 (61.9%)"  "8 (38.1%)"  
#stomach                       "10 (62.5%)"  "6 (37.5%)"  
#kidney except renal pelvis    "19 (63.3%)"  "11 (36.7%)" 
#tonsil                        "18 (64.3%)"  "10 (35.7%)" 
#oropharynx                    "11 (64.7%)"  "6 (35.3%)"  
#prostate                      "78 (65.5%)"  "41 (34.5%)" 
#melanoma of skin              "32 (68.1%)"  "15 (31.9%)" 
#breast                        "161 (68.5%)" "74 (31.5%)" 
#colorectal                    "47 (71.2%)"  "19 (28.8%)" 
#corpus uteri                  "23 (71.9%)"  "9 (28.1%)"  
#ovary                         "30 (73.2%)"  "11 (26.8%)" 
#bronchus lung                 "66 (73.3%)"  "24 (26.7%)" 
#multiple myeloma plasma cell  "44 (73.3%)"  "16 (26.7%)" 
#brain                         "28 (75.7%)"  "9 (24.3%)"  
#other                         "180 (77.3%)" "53 (22.7%)" 
#other connective soft tissue  "35 (77.8%)"  "10 (22.2%)" 
#esophagus                     "27 (79.4%)"  "7 (20.6%)"  
#pancreas                      "40 (81.6%)"  "9 (18.4%)"  
#neuroendocrine tumors         "18 (85.7%)"  "3 (14.3%)"  
#liver intrahepatic bile ducts "34 (87.2%)"  "5 (12.8%)"  
#bladder                       "42 (87.5%)"  "6 (12.5%)"  
#unspecified skin              "21 (87.5%)"  "3 (12.5%)"  
#lymphoid leukemia             "15 (88.2%)"  "2 (11.8%)"


kbl(
  combined,
  align = c("l", "c", "c"),
  caption = "Toposcore accross cancer class (column %)"
) |>
  kable_classic(full_width = FALSE, html_font = "Arial Narrow") |>
  footnote(
    general = "SIG1+ = unhealthy-associated toposcore; SIG2+ = healthy-associated toposcore.",
    threeparttable = TRUE
  )


#-------------------------------------------------------------------------------
#make barplot healthy / cancer using SIG_class


#Bar with SIG1 on bottom, grey in middle, and SIG2 on the top
tab <- table(meta_with_topo$Cancer_Healthy,
             meta_with_topo$SIG_class)

tab_perc <- apply(tab, 1, function(x) 100 * x / sum(x, na.rm = TRUE))
tab_perc <- tab_perc[c("SIG2","Gray","SIG1"),]
tab_perc <- tab_perc[,rev(order(tab_perc["SIG1",]))]


tab_binary <- table(meta_with_topo$Cancer_Healthy,
                    meta_with_topo$Toposcore)

tab_binary <- tab_binary[colnames(tab_perc), , drop = FALSE]

tab_binary_perc <- apply(tab_binary, 1, function(x) 100 * x / sum(x, na.rm = TRUE))
tab_binary_perc <- t(tab_binary_perc)  # rows = cancers, cols = levels of Toposcore
line_level <- "SIG2+"
line_vals <- tab_binary_perc[, line_level]

dot_vals <- line_vals[colnames(tab_perc)]


ref_row <- "Healthy"
pvals <- sapply(rownames(tab), function(r) {
  if (r == ref_row) return(NA_real_)   # no test for the reference
  
  # 2 x k table: Healthy vs current cancer across all SIG columns
  m <- rbind(tab[ref_row, ], tab[r, ])
  fisher.test(m)$p.value
})
qvals <- p.adjust(pvals, method = "fdr")
qvals <- qvals[colnames(tab_perc)]

stars <- qvals

stars[which(qvals < 0.1)] <- "*"
stars[which(qvals < 0.01)] <- "**"
stars[which(qvals < 0.001)] <- "***"
stars[is.na(qvals)] <- ""
stars[which(qvals > 0.1)] <- "" 


par(mar=c(5,8,1,3), xpd = TRUE)
bp <- barplot(tab_perc,
              horiz  = TRUE,  
              col = c("#0072B2", "darkgrey", "#CC79A7"),
              border = "white",
              space = 0.2,
              xlab = "Toposcore percentages", las=1)

points(x = dot_vals, y = bp,
       pch = 18, col = "#0072B2", cex = 1.5)

legend("bottomleft",
       inset  = c(-0.3, -0.8),  # tweak this value
       legend = rownames(tab_perc),
       fill   = c("#0072B2", "darkgrey", "#CC79A7"),
       bty    = "n")

legend("bottomleft",
       inset  = c(-0.2875, -0.97),
       legend = c("", "", "", "SIG1+/SIG2+ split"), 
       pch    = c(NA, NA, NA, 18),
       col    = c(NA, NA, NA, "#0072B2"),
       pt.cex = c(NA, NA, NA, 2),
       bty    = "n")

# add text a bit to the right of each bar
text(x     = 105,                
     y     = bp,                         
     labels = paste0(stars),
     xpd   = TRUE,
     cex   = 0.8, adj = 0)



#-------------------------------------------------------------------------------
#make barplot cancer classes using SIG_class


#Bar with SIG1 on bottom, grey in middle, and SIG2 on the top
tab <- table(meta_with_topo$icd10_first_3_name_short_cutoff,
             meta_with_topo$SIG_class)

tab_perc <- apply(tab, 1, function(x) 100 * x / sum(x, na.rm = TRUE))
tab_perc <- tab_perc[c("SIG2","Gray","SIG1"),]
tab_perc <- tab_perc[,rev(order(tab_perc["SIG1",]))]


#add extra bars in to indicate the SIG1+ or SIG2+ binary split
tab_binary <- table(meta_with_topo$icd10_first_3_name_short_cutoff,
                    meta_with_topo$Toposcore)

tab_binary <- tab_binary[colnames(tab_perc), , drop = FALSE]

tab_binary_perc <- apply(tab_binary, 1, function(x) 100 * x / sum(x, na.rm = TRUE))
tab_binary_perc <- t(tab_binary_perc)  # rows = cancers, cols = levels of Toposcore
line_level <- "SIG2+"
line_vals <- tab_binary_perc[, line_level]

dot_vals <- line_vals[colnames(tab_perc)]



ref_row <- "healthy"
pvals <- sapply(rownames(tab), function(r) {
  if (r == ref_row) return(NA_real_)   # no test for the reference
  
  # 2 x k table: Healthy vs current cancer across all SIG columns
  m <- rbind(tab[ref_row, ], tab[r, ])
  fisher.test(m)$p.value
})
qvals <- p.adjust(pvals, method = "fdr")
qvals <- qvals[colnames(tab_perc)]

stars <- qvals

stars[which(qvals < 0.1)] <- "*"
stars[which(qvals < 0.01)] <- "**"
stars[which(qvals < 0.001)] <- "***"
stars[is.na(qvals)] <- ""
stars[which(qvals > 0.1)] <- "" 


par(mar=c(5,13,1,7), xpd = TRUE)
bp <- barplot(tab_perc,
              horiz  = TRUE,  
              col = c("#0072B2", "darkgrey", "#CC79A7"),
              border = "white",
              space = 0.2,
              xlab = "Toposcore percentages by cancer class", las=1)

points(x = dot_vals, y = bp,
       pch = 18, col = "#0072B2", cex = 1.5)

legend("bottomleft",
       inset  = c(-0.6, -0.15),  # tweak this value
       legend = rownames(tab_perc),
       fill   = c("#0072B2", "darkgrey", "#CC79A7"),
       bty    = "n")

legend("bottomleft",
       inset  = c(-0.585, -0.2),
       legend = c("", "", "", "SIG1+/SIG2+ split"), 
       pch    = c(NA, NA, NA, 18),
       col    = c(NA, NA, NA, "#0072B2"),
       pt.cex = c(NA, NA, NA, 2),
       bty    = "n")

# add text a bit to the right of each bar
text(x     = 105,                
     y     = bp,                         
     labels = paste0(stars),
     xpd   = TRUE,
     cex   = 0.8, adj = 0)






#-------------------------------------------------------------------------------
#Check if the species identified as dysbiotic are enriched in all cancers or specific cancers or associated with co-morbidities.

#Toposcore SIG1 species
SIG1_species <- c(
  "Actinomyces_graevenitzii",
  "Alloscardovia_omnicolens",
  "Anaerostipes_caccae",
  "Bifidobacterium_dentium",
  "Blautia_producta",
  "Campylobacter_concisus",
  "Campylobacter_gracilis",
  "Clostridium_innocuum",
  "Clostridium_perfringens",
  "Clostridium_scindens",
  "Clostridium_symbiosum",
  "Collinsella_SGB14754",
  "Enorma_massiliensis",
  "Enterocloster_aldensis",
  "Enterocloster_bolteae",
  "Enterocloster_clostridioformis",
  "Erysipelatoclostridium_ramosum",
  "Fournierella_massiliensis",
  "Granulicatella_adiacens",
  "Hungatella_hathewayi",
  "Lacticaseibacillus_paracasei",
  "Lactobacillus_gasseri",
  "Lactobacillus_vaginalis",
  "Ligilactobacillus_salivarius",
  "Limosilactobacillus_fermentum",
  "Limosilactobacillus_oris",
  "Megasphaera_micronuciformis",
  "Proteus_mirabilis",
  "Streptococcus_anginosus",
  "Streptococcus_gordonii",
  "Streptococcus_mutans",
  "Streptococcus_oralis",
  "Streptococcus_parasanguinis",
  "Streptococcus_salivarius",
  "Veillonella_atypica",
  "Veillonella_dispar",
  "Veillonella_parvula"
)
#add s__
SIG1_species <- paste0("s__", SIG1_species)
length(SIG1_species)


#there are the results from the comparisons
filename <- "./Code/Submission/Supplementary tables/Supplementary Table 3.xlsx"
getSheetNames(filename)

#---------------------
#note called tax_df here but that is not descriptive; these are results from differential abundance comparisons

#1) enriched in all cancers  
#TableS2-1	Species; cancer patients vs healthy controls (reference) while adjusting for confounders.
tax_df_allcancer <- read.xlsx(filename, sheet="TableS3-1", rowNames = F)
dim(tax_df_allcancer) #one row per species
#369  12

table(SIG1_species %in% tax_df_allcancer$Species)
#FALSE  TRUE  #25 of the 37 species detected
#12    25 

tax_df_allcancer_SIG1 <- tax_df_allcancer[tax_df_allcancer$Species %in% SIG1_species,]
table(tax_df_allcancer_SIG1$EffectSize > 0 & tax_df_allcancer_SIG1$Qvalue < 0.1)
#FALSE  TRUE  #24/25 positively associated to cancer
#   1    24 


#from below

#tax_df_allcancer_SIG1 <- tax_df_allcancer[tax_df_allcancer$Species %in% matching_SIG1_species_17,]
#table(tax_df_allcancer_SIG1$EffectSize > 0 & tax_df_allcancer_SIG1$Qvalue < 0.1)


#---------------------
#2) or specific cancers
#TableS2-3	Species; cancer class vs healthy controls (reference) while adjusting for confounders.
tax_df_cancerclasses <- read.xlsx(filename, sheet="TableS3-2", rowNames = F)
dim(tax_df_cancerclasses)
#8569   12


table(SIG1_species %in% tax_df_cancerclasses$Species)
#FALSE  TRUE  #21 of the 37 species detected
# 16    21 

tax_df_cancerclasses_SIG1 <- tax_df_cancerclasses[tax_df_cancerclasses$Species %in% SIG1_species,]
table(tax_df_cancerclasses_SIG1$EffectSize > 0 & tax_df_cancerclasses_SIG1$Qvalue < 0.1)
#FALSE  TRUE 
# 131   234 

table(tax_df_cancerclasses_SIG1$Comparison, tax_df_cancerclasses_SIG1$EffectSize > 0 & tax_df_cancerclasses_SIG1$Qvalue < 0.1)
#                                                            FALSE TRUE
#base of tongue- healthy controls(reference)                    3   12
#bladder- healthy controls(reference)                           4   12
#brain- healthy controls(reference)                             3   12
#breast- healthy controls(reference)                            4   16
#bronchus lung- healthy controls(reference)                     1   18
#colorectal- healthy controls(reference)                        4   13
#corpus uteri- healthy controls(reference)                     10    6
#esophagus- healthy controls(reference)                         3   16
#kidney except renal pelvis- healthy controls(reference)       10    7
#liver intrahepatic bile ducts- healthy controls(reference)     2   15
#lymphoid leukemia- healthy controls(reference)                 4   11
#melanoma of skin- healthy controls(reference)                  8    8
#multiple myeloma plasma cell- healthy controls(reference)      3   13
#neuroendocrine tumors- healthy controls(reference)             7   10
#oropharynx- healthy controls(reference)                        9    6
#other connective soft tissue- healthy controls(reference)      6   10
#ovary- healthy controls(reference)                             3   13
#pancreas- healthy controls(reference)                          2   16
#prostate- healthy controls(reference)                         13    5
#stomach- healthy controls(reference)                          12    4
#tonsil- healthy controls(reference)                           10    5
#unspecified skin- healthy controls(reference)                 10    6


#every cancer class has some SIG1 species positively associated with them


#---------------------
#3a) or associated with co-morbidities; summed ECI score
#TableS2-7	Species associated with Elixhauser Comorbidity score while adjusting for confounders.
tax_df_ECI <- read.xlsx(filename, sheet="TableS3-4", rowNames = F)
dim(tax_df_ECI)
#342  12

table(SIG1_species %in% tax_df_ECI$Species)
#FALSE  TRUE 
# 12    25  #25 out of 37 in here

SIG1_species_in_comparisons <- SIG1_species[SIG1_species %in% tax_df_ECI$Species]


tax_df_ECI_SIG1 <- tax_df_ECI[tax_df_ECI$Species %in% SIG1_species,]
table(tax_df_ECI_SIG1$EffectSize > 0 & tax_df_ECI_SIG1$Qvalue < 0.1)
#FALSE  TRUE  #12 positively associated with Elix score
# 13    12 

#which species
tax_df_ECI_SIG1[tax_df_ECI_SIG1$EffectSize > 0 & tax_df_ECI_SIG1$Qvalue < 0.1, "Species"]
#[1] "s__Clostridium_scindens"           "s__Erysipelatoclostridium_ramosum" "s__Enterocloster_clostridioformis" "s__Enterocloster_bolteae"         
#[5] "s__Anaerostipes_caccae"            "s__Streptococcus_parasanguinis"    "s__Clostridium_symbiosum"          "s__Limosilactobacillus_fermentum" 
#[9] "s__Veillonella_parvula"            "s__Streptococcus_mutans"           "s__Blautia_producta"               "s__Lactobacillus_gasseri"  


#---------------------
#3b) or associated with co-morbidities; Elix components
#TableS2-9	Species associated with Elixhauser Comorbidity component (reference-No) while adjusting for confounders.
tax_df_ECIcomponents <- read.xlsx(filename, sheet="TableS3-5", rowNames = F)
dim(tax_df_ECIcomponents)
#10268    12

tax_df_ECIcomponents_SIG1 <- tax_df_ECIcomponents[tax_df_ECIcomponents$Species %in% SIG1_species,]


#None of the SIG1 species are negatively associated with any component 

table(tax_df_ECIcomponents_SIG1$Comparison, tax_df_ECIcomponents_SIG1$EffectSize > 0 & tax_df_ECIcomponents_SIG1$Qvalue < 0.1)

#8 Elixhauser_FluidsLytes, 5 Elixhauser_Anemia and 2 Elixhauser_WeightLoss, 1 Elixhauser_DMcx, 1 Elixhauser_Rheumatic have species

#                                                FALSE TRUE
#Elixhauser_Alcohol: Yes vs No(reference)         25    0
#Elixhauser_Anemia: Yes vs No(reference)          20    5
#Elixhauser_Arrhythmia: Yes vs No(reference)      25    0
#Elixhauser_BloodLoss: Yes vs No(reference)       25    0
#Elixhauser_CHF: Yes vs No(reference)             25    0
#Elixhauser_Coagulopathy: Yes vs No(reference)    25    0
#Elixhauser_Depression: Yes vs No(reference)      25    0
#Elixhauser_DM: Yes vs No(reference)              25    0
#Elixhauser_DMcx: Yes vs No(reference)            24    1
#Elixhauser_Drugs: Yes vs No(reference)           25    0
#Elixhauser_FluidsLytes: Yes vs No(reference)     17    8
#Elixhauser_HIV: Yes vs No(reference)             25    0
#Elixhauser_HTN: Yes vs No(reference)             25    0
#Elixhauser_Hypothyroid: Yes vs No(reference)     25    0
#Elixhauser_Liver: Yes vs No(reference)           25    0
#Elixhauser_Lymphoma: Yes vs No(reference)        25    0
#Elixhauser_Mets: Yes vs No(reference)            25    0
#Elixhauser_NeuroOther: Yes vs No(reference)      25    0
#Elixhauser_Obesity: Yes vs No(reference)         25    0
#Elixhauser_Paralysis: Yes vs No(reference)       25    0
#Elixhauser_PHTN: Yes vs No(reference)            25    0
#Elixhauser_Psychoses: Yes vs No(reference)       25    0
#Elixhauser_PUD: Yes vs No(reference)             25    0
#Elixhauser_Pulmonary: Yes vs No(reference)       25    0
#Elixhauser_PVD: Yes vs No(reference)             25    0
#Elixhauser_Renal: Yes vs No(reference)           25    0
#Elixhauser_Rheumatic: Yes vs No(reference)       24    1
#Elixhauser_Tumor: Yes vs No(reference)           25    0
#Elixhauser_Valvular: Yes vs No(reference)        25    0
#Elixhauser_WeightLoss: Yes vs No(reference)      23    2


#which species
tax_df_ECIcomponents_SIG1[tax_df_ECIcomponents_SIG1$Comparison == "Elixhauser_FluidsLytes: Yes vs No(reference)" & 
                            tax_df_ECIcomponents_SIG1$Qvalue < 0.1, "Species"]
#1] "s__Clostridium_scindens"           "s__Enterocloster_clostridioformis" "s__Enterocloster_bolteae"          "s__Anaerostipes_caccae"           
#[5] "s__Hungatella_hathewayi"           "s__Limosilactobacillus_fermentum"  "s__Blautia_producta"               "s__Lactobacillus_gasseri" 


tax_df_ECIcomponents_SIG1[tax_df_ECIcomponents_SIG1$Comparison == "Elixhauser_Anemia: Yes vs No(reference)" & 
                            tax_df_ECIcomponents_SIG1$Qvalue < 0.1, "Species"]
#[1] "s__Streptococcus_salivarius"    "s__Streptococcus_parasanguinis" "s__Veillonella_atypica"         "s__Veillonella_parvula"        
#[5] "s__Streptococcus_anginosus" 

tax_df_ECIcomponents_SIG1[tax_df_ECIcomponents_SIG1$Comparison == "Elixhauser_DMcx: Yes vs No(reference)" & 
                            tax_df_ECIcomponents_SIG1$Qvalue < 0.1, "Species"]
#"s__Lactobacillus_gasseri"

tax_df_ECIcomponents_SIG1[tax_df_ECIcomponents_SIG1$Comparison == "Elixhauser_WeightLoss: Yes vs No(reference)" & 
                            tax_df_ECIcomponents_SIG1$Qvalue < 0.1, "Species"]
#"s__Clostridium_symbiosum" "s__Lactobacillus_gasseri"


#-------------------------------------------------------------------------------
#write up

#Out of the 37 SIG1 species, 25 are present in metaphlan comparisons which are filtered by prevalence and abundance. Of these 25, 19 can be directly 
#matched to Kraken taxonomy by name and strong correlation between species abundances from both taxonomy pipeline results but only 17 are present in all
#comparisons after prevalence and abundance filtering. These are Actinomyces graevenitzii, Alloscardovia omnicolens, Anaerostipes caccae, Bifidobacterium dentium, 
#Enterocloster bolteae, Enterocloster clostridioformis, Erysipelatoclostridium ramosum, Lacticaseibacillus paracasei, Lactobacillus gasseri, 
#Limosilactobacillus fermentum, Streptococcus anginosus, Streptococcus gordonii, Streptococcus mutans, Streptococcus parasanguinis, 
#Streptococcus salivarius, Veillonella atypica, and Veillonella parvula. 

#16/17 positively associated to cancer in both Kraken and Metaphlan taxonomy (Supplementary table X). 




#-------------------------------------------------------------------------------
#see how well the SIG1 species from Metaphlan taxonomy match with GTDB taxonomy


#load the metaphlan species table
filename <- "~/Dropbox/Mayo_RS/R/Oncobiome full cohort analysis/data/Metaphlan_taxonomy_MCCM_species.csv"
tax_table_RA <- read.csv(filename)

rownames(tax_table_RA) <- tax_table_RA$X
tax_table_RA <- tax_table_RA[,2:ncol(tax_table_RA)]


table(SIG1_species %in% rownames(tax_table_RA))
#FALSE  TRUE 
#   1    36  #36 of the 37 are annotated

#filter for sparse and then check again

#should not be 0 in at least 10% of the samples
tax_table_RA_prev <- tax_table_RA[which(rowSums(!tax_table_RA == 0) > ncol(tax_table_RA) * 0.1),]
dim(tax_table_RA); dim(tax_table_RA_prev)
#[1] 1854 1651
#[1] 369 1651

rownames(tax_table_RA_prev)

table(SIG1_species %in% rownames(tax_table_RA_prev))
#FALSE  TRUE  #only 25 can really be considered
#12    25 

SIG1_species[SIG1_species %in% rownames(tax_table_RA_prev)]
#[1] "s__Actinomyces_graevenitzii"       "s__Alloscardovia_omnicolens"       "s__Anaerostipes_caccae"            "s__Bifidobacterium_dentium"       
#[5] "s__Blautia_producta"               "s__Clostridium_innocuum"           "s__Clostridium_scindens"           "s__Clostridium_symbiosum"         
#[9] "s__Enterocloster_aldensis"         "s__Enterocloster_bolteae"          "s__Enterocloster_clostridioformis" "s__Erysipelatoclostridium_ramosum"
#[13] "s__Hungatella_hathewayi"           "s__Lacticaseibacillus_paracasei"   "s__Lactobacillus_gasseri"          "s__Limosilactobacillus_fermentum" 
#[17] "s__Streptococcus_anginosus"        "s__Streptococcus_gordonii"         "s__Streptococcus_mutans"           "s__Streptococcus_oralis"          
#[21] "s__Streptococcus_parasanguinis"    "s__Streptococcus_salivarius"       "s__Veillonella_atypica"            "s__Veillonella_dispar"            
#[25] "s__Veillonella_parvula" 



#-------------------------------------

dim(tax_table_RA); dim(tax_table_Kraken)
#[1] 1854 1651
#[1] 7355 1644

tax_table_RA_sub <- tax_table_RA[, colnames(tax_table_Kraken)]

table(colnames(tax_table_RA_sub) == colnames(tax_table_Kraken))
#TRUE 
#1644 


table(SIG1_species %in% rownames(tax_table_Kraken))
#FALSE  TRUE 
#12    25 

SIG1_species[SIG1_species %in% rownames(tax_table_Kraken)]
#[1] "s__Actinomyces_graevenitzii"       "s__Alloscardovia_omnicolens"       "s__Anaerostipes_caccae"            "s__Bifidobacterium_dentium"       
#[5] "s__Enorma_massiliensis"            "s__Enterocloster_bolteae"          "s__Enterocloster_clostridioformis" "s__Erysipelatoclostridium_ramosum"
#[9] "s__Fournierella_massiliensis"      "s__Granulicatella_adiacens"        "s__Lacticaseibacillus_paracasei"   "s__Lactobacillus_gasseri"         
#[13] "s__Ligilactobacillus_salivarius"   "s__Limosilactobacillus_fermentum"  "s__Limosilactobacillus_oris"       "s__Proteus_mirabilis"             
#[17] "s__Streptococcus_anginosus"        "s__Streptococcus_gordonii"         "s__Streptococcus_mutans"           "s__Streptococcus_oralis"          
#[21] "s__Streptococcus_parasanguinis"    "s__Streptococcus_salivarius"       "s__Veillonella_atypica"            "s__Veillonella_dispar"            
#[25] "s__Veillonella_parvula"  

#for these ones see which match
table(SIG1_species_in_comparisons %in% SIG1_species[SIG1_species %in% rownames(tax_table_Kraken)]) 
#19 of the ones that are in the comparisons above (filtered by prevalence and abundance) are also directly matched by name

SIG1_species_in_comparisons_matched <- SIG1_species_in_comparisons[SIG1_species_in_comparisons %in% SIG1_species[SIG1_species %in% rownames(tax_table_Kraken)]]
#[1] "s__Actinomyces_graevenitzii"       "s__Alloscardovia_omnicolens"       "s__Anaerostipes_caccae"            "s__Bifidobacterium_dentium"       
#[5] "s__Enterocloster_bolteae"          "s__Enterocloster_clostridioformis" "s__Erysipelatoclostridium_ramosum" "s__Lacticaseibacillus_paracasei"  
#[9] "s__Lactobacillus_gasseri"          "s__Limosilactobacillus_fermentum"  "s__Streptococcus_anginosus"        "s__Streptococcus_gordonii"        
#[13] "s__Streptococcus_mutans"           "s__Streptococcus_oralis"           "s__Streptococcus_parasanguinis"    "s__Streptococcus_salivarius"      
#[17] "s__Veillonella_atypica"            "s__Veillonella_dispar"             "s__Veillonella_parvula" 


SIG1_species[!SIG1_species %in% SIG1_species_in_comparisons_matched]


#checked these for GTDB taxonomy; for some there is a clear explanation why they do not match but others should have the same name; may just be too sparse 
#s__Lactobacillus_vaginalis = Ligilactobacillus vaginalis
#s__Clostridium_innocuum = Clostridium_AQ innocuum



#[1] "s__Blautia_producta"             "s__Campylobacter_concisus"       "s__Campylobacter_gracilis"       "s__Clostridium_innocuum"        
#[5] "s__Clostridium_perfringens"      "s__Clostridium_scindens"         "s__Clostridium_symbiosum"        "s__Collinsella_SGB14754"        
#[9] "s__Enorma_massiliensis"          "s__Enterocloster_aldensis"       "s__Fournierella_massiliensis"    "s__Granulicatella_adiacens"     
#[13] "s__Hungatella_hathewayi"         "s__Lactobacillus_vaginalis"      "s__Ligilactobacillus_salivarius" "s__Limosilactobacillus_oris"    
#[17] "s__Megasphaera_micronuciformis"  "s__Proteus_mirabilis" 


#look for the correlation for these

cor_list <- list()
for (i in 1:length(SIG1_species_in_comparisons_matched)) {
  temp_species <- SIG1_species_in_comparisons_matched[i]
  cor_list[[i]] <- cor(as.numeric(tax_table_RA_sub[temp_species,]), as.numeric(tax_table_RA_Kraken[temp_species,]))
  
  plot(as.numeric(tax_table_RA_sub[temp_species,]), as.numeric(tax_table_RA_Kraken[temp_species,]), main=temp_species)
  
}
names(cor_list) <- SIG1_species_in_comparisons_matched


#Metaphlan seems to allocate larger percentages of RA to species then Kraken which makes sense since there are fewer species in it's database

comp_df <- as.data.frame(do.call(rbind, cor_list))
#[,1]
#s__Actinomyces_graevenitzii       0.9997377
#s__Alloscardovia_omnicolens       0.9959663
#s__Anaerostipes_caccae            0.9851581
#s__Bifidobacterium_dentium        0.9932062
#s__Enterocloster_bolteae          0.9786586
#s__Enterocloster_clostridioformis 0.9507793
#s__Erysipelatoclostridium_ramosum 0.9936461
#s__Lacticaseibacillus_paracasei   0.9908975
#s__Lactobacillus_gasseri          0.7901199 #not good
#s__Limosilactobacillus_fermentum  0.9759151
#s__Streptococcus_anginosus        0.8845501
#s__Streptococcus_gordonii         0.9431142
#s__Streptococcus_mutans           0.9887943
#s__Streptococcus_oralis           0.7129906 #not good
#s__Streptococcus_parasanguinis    0.7858507 #not good
#s__Streptococcus_salivarius       0.9492192
#s__Veillonella_atypica            0.9934317
#s__Veillonella_dispar             0.9331876
#s__Veillonella_parvula            0.9894055

rownames(comp_df)[comp_df[,1] > 0.85]
#[1] "s__Actinomyces_graevenitzii"       "s__Alloscardovia_omnicolens"       "s__Anaerostipes_caccae"            "s__Bifidobacterium_dentium"       
#[5] "s__Enterocloster_bolteae"          "s__Enterocloster_clostridioformis" "s__Erysipelatoclostridium_ramosum" "s__Lacticaseibacillus_paracasei"  
#[9] "s__Limosilactobacillus_fermentum"  "s__Streptococcus_anginosus"        "s__Streptococcus_gordonii"         "s__Streptococcus_mutans"          
#[13] "s__Streptococcus_salivarius"       "s__Veillonella_atypica"            "s__Veillonella_dispar"             "s__Veillonella_parvula" 


#these 16 can be directly compared out of the 37 SIG1 species and most are indeed associated to cancer and / or comorbidities
#all 19 are OK to compare


#-------------------------------------------------------------------------------
#check these 19 species for associations with cancer and comorbidities in Kraken results; does that look even better?

#SIG1_species_in_comparisons_matched


#there are the results from the comparisons
filename <- "./Code/Submission/Supplementary tables/Supplementary Table 2.xlsx"
getSheetNames(filename)


#---------------------
#note called tax_df here but that is not descriptive; these are results from differential abundance comparisons

#1) enriched in all cancers  
#TableS2-1	Species; cancer patients vs healthy controls (reference) while adjusting for confounders.
tax_df_allcancer <- read.xlsx(filename, sheet="TableS2-1", rowNames = F)
dim(tax_df_allcancer) #one row per species
#1278  11
tax_df_allcancer$Species <- gsub(" ", "_", tax_df_allcancer$Species)

table(SIG1_species_in_comparisons_matched %in% tax_df_allcancer$Species)
#FALSE  TRUE  #17 of the 19 species detected
#   2    17 

tax_df_allcancer_SIG1 <- tax_df_allcancer[tax_df_allcancer$Species %in% SIG1_species_in_comparisons_matched,]
table(tax_df_allcancer_SIG1$EffectSize > 0 & tax_df_allcancer_SIG1$Qvalue < 0.1)
#FALSE  TRUE  #16/17 positively associated to cancer
#   1    16 


#how many or these also in Metaphlan taxonomy results?




matching_SIG1_species_17 <- tax_df_allcancer$Species[tax_df_allcancer$Species %in% SIG1_species_in_comparisons_matched]
#[1] "s__Actinomyces_graevenitzii"       "s__Alloscardovia_omnicolens"       "s__Anaerostipes_caccae"            "s__Bifidobacterium_dentium"       
#[5] "s__Enterocloster_bolteae"          "s__Enterocloster_clostridioformis" "s__Erysipelatoclostridium_ramosum" "s__Lacticaseibacillus_paracasei"  
#[9] "s__Lactobacillus_gasseri"          "s__Limosilactobacillus_fermentum"  "s__Streptococcus_anginosus"        "s__Streptococcus_gordonii"        
#[13] "s__Streptococcus_mutans"           "s__Streptococcus_parasanguinis"    "s__Streptococcus_salivarius"       "s__Veillonella_atypica"           
#[17] "s__Veillonella_parvula"   


#---------------------
#2) or specific cancers
#TableS2-3	Species; cancer class vs healthy controls (reference) while adjusting for confounders.
tax_df_cancerclasses <- read.xlsx(filename, sheet="TableS2-3", rowNames = F)
dim(tax_df_cancerclasses)
#19868   11
tax_df_cancerclasses$Species <- gsub(" ", "_", tax_df_cancerclasses$Species)


table(SIG1_species_in_comparisons_matched %in% tax_df_cancerclasses$Species)
#FALSE  TRUE  #15 of the 19 species detected
# 4     15 

SIG1_species_in_comparisons_matched[SIG1_species_in_comparisons_matched %in% tax_df_cancerclasses$Species]
#these 15

#[1] "s__Alloscardovia_omnicolens"       "s__Anaerostipes_caccae"            "s__Bifidobacterium_dentium"        "s__Enterocloster_bolteae"         
#[5] "s__Enterocloster_clostridioformis" "s__Erysipelatoclostridium_ramosum" "s__Lacticaseibacillus_paracasei"   "s__Limosilactobacillus_fermentum" 
#[9] "s__Streptococcus_anginosus"        "s__Streptococcus_gordonii"         "s__Streptococcus_mutans"           "s__Streptococcus_parasanguinis"   
#[13] "s__Streptococcus_salivarius"       "s__Veillonella_atypica"            "s__Veillonella_parvula"  



tax_df_cancerclasses_SIG1 <- tax_df_cancerclasses[tax_df_cancerclasses$Species %in% SIG1_species_in_comparisons_matched,]
table(tax_df_cancerclasses_SIG1$EffectSize > 0 & tax_df_cancerclasses_SIG1$Qvalue < 0.1)
#FALSE  TRUE 
# 36   188 

tab <- table(tax_df_cancerclasses_SIG1$Comparison, tax_df_cancerclasses_SIG1$EffectSize > 0 & tax_df_cancerclasses_SIG1$Qvalue < 0.1)

#every cancer class has some SIG1 species positively associated with them
#                                                             FALSE TRUE
#base of tongue- healthy controls(reference)                    2    6
#bladder- healthy controls(reference)                           1    8
#brain- healthy controls(reference)                             2    6
#breast- healthy controls(reference)                            1   11
#bronchus lung- healthy controls(reference)                     0   12
#colorectal- healthy controls(reference)                        1    9
#corpus uteri- healthy controls(reference)                      3    9
#esophagus- healthy controls(reference)                         0   12
#kidney except renal pelvis- healthy controls(reference)        2    8
#liver intrahepatic bile ducts- healthy controls(reference)     0   12
#lymphoid leukemia- healthy controls(reference)                 1    7
#melanoma of skin- healthy controls(reference)                  2    8
#multiple myeloma plasma cell- healthy controls(reference)      1    8
#neuroendocrine tumors- healthy controls(reference)             1   10
#oropharynx- healthy controls(reference)                        0    9
#other connective soft tissue- healthy controls(reference)      3    5
#ovary- healthy controls(reference)                             1    8
#pancreas- healthy controls(reference)                          0   13
#prostate- healthy controls(reference)                          3    8
#stomach- healthy controls(reference)                           7    5
#tonsil- healthy controls(reference)                            4    4
#unspecified skin- healthy controls(reference)                  1   10

median(tab[,"TRUE"])


#---------------------
#3a) or associated with co-morbidities; summed ECI score
#TableS2-7	Species associated with Elixhauser Comorbidity score while adjusting for confounders.
tax_df_ECI <- read.xlsx(filename, sheet="TableS2-7", rowNames = F)
dim(tax_df_ECI)
#1237  11
tax_df_ECI$Species <- gsub(" ", "_", tax_df_ECI$Species)

table(SIG1_species_in_comparisons_matched %in% tax_df_ECI$Species)
#FALSE  TRUE 
# 2     17  #17 out of 19 in here

SIG1_species_in_comparisons <- SIG1_species_in_comparisons_matched[SIG1_species_in_comparisons_matched %in% tax_df_ECI$Species]

tax_df_ECI_SIG1 <- tax_df_ECI[tax_df_ECI$Species %in% SIG1_species_in_comparisons_matched,]
table(tax_df_ECI_SIG1$EffectSize > 0 & tax_df_ECI_SIG1$Qvalue < 0.1)
#FALSE  TRUE  #9 positively associated with Elix score using Kraken taxonomy
# 8     9 

#which species
tax_df_ECI_SIG1[tax_df_ECI_SIG1$EffectSize > 0 & tax_df_ECI_SIG1$Qvalue < 0.1, "Species"]


#these numbers look very consistent and underestimates since not all species in Metaphlan taxonomy are present in Kraken taxonomy and vice versa

#which SIG1 species; out of the 9 species in Kraken and 12 in Metaphlan 7 overlap
#Kraken
#[1] "s__Bifidobacterium_dentium"        "s__Enterocloster_bolteae"    *matches with mp        "s__Enterocloster_clostridioformis" *matches with mp 
#"s__Erysipelatoclostridium_ramosum" *matches with mp  #[5] "s__Lactobacillus_gasseri"    *matches with mp        
#"s__Limosilactobacillus_fermentum"   *matches with mp  "s__Streptococcus_mutans"   *matches with mp          "s__Veillonella_atypica"           
#[9] "s__Veillonella_parvula"  *matches with mp 




#---------------------
#3b) or associated with co-morbidities; Elix components
#TableS2-9	Species associated with Elixhauser Comorbidity component (reference-No) while adjusting for confounders.
tax_df_ECIcomponents <- read.xlsx(filename, sheet="TableS2-9", rowNames = F)
dim(tax_df_ECIcomponents)
#37122    11
tax_df_ECIcomponents$Species <- gsub(" ", "_", tax_df_ECIcomponents$Species)


tax_df_ECIcomponents_SIG1 <- tax_df_ECIcomponents[tax_df_ECIcomponents$Species %in% SIG1_species_in_comparisons_matched,]

#Only 1 of the SIG1 species are negatively associated with any component 1 species to Elixhauser_FluidsLytes; besides that only positive associations
table(tax_df_ECIcomponents_SIG1$Comparison, tax_df_ECIcomponents_SIG1$EffectSize > 0 & tax_df_ECIcomponents_SIG1$Qvalue < 0.1)

#7 Elixhauser_FluidsLytes, 5 Elixhauser_Anemia and 7 Elixhauser_WeightLoss, 1 Elixhauser_Arrhythmia, 2 Elixhauser_Liver, 
#4 Elixhauser_NeuroOther, 1 Elixhauser_Psychoses, 3 Elixhauser_PUD, 2 Elixhauser_Renal

#                                                FALSE TRUE
#Elixhauser_Alcohol: Yes vs No(reference)         17    0
#Elixhauser_Anemia: Yes vs No(reference)          12    5
#Elixhauser_Arrhythmia: Yes vs No(reference)      16    1
#Elixhauser_BloodLoss: Yes vs No(reference)       17    0
#Elixhauser_CHF: Yes vs No(reference)             17    0
#Elixhauser_Coagulopathy: Yes vs No(reference)    17    0
#Elixhauser_Depression: Yes vs No(reference)      17    0
#Elixhauser_DM: Yes vs No(reference)              17    0
#Elixhauser_DMcx: Yes vs No(reference)            17    0
#Elixhauser_Drugs: Yes vs No(reference)           17    0
#Elixhauser_FluidsLytes: Yes vs No(reference)     10    7
#Elixhauser_HIV: Yes vs No(reference)             17    0
#Elixhauser_HTN: Yes vs No(reference)             17    0
#Elixhauser_Hypothyroid: Yes vs No(reference)     17    0
#Elixhauser_Liver: Yes vs No(reference)           15    2
#Elixhauser_Lymphoma: Yes vs No(reference)        17    0
#Elixhauser_Mets: Yes vs No(reference)            17    0
#Elixhauser_NeuroOther: Yes vs No(reference)      13    4
#Elixhauser_Obesity: Yes vs No(reference)         17    0
#Elixhauser_Paralysis: Yes vs No(reference)       17    0
#Elixhauser_PHTN: Yes vs No(reference)            17    0
#Elixhauser_Psychoses: Yes vs No(reference)       16    1
#Elixhauser_PUD: Yes vs No(reference)             14    3
#Elixhauser_Pulmonary: Yes vs No(reference)       17    0
#Elixhauser_PVD: Yes vs No(reference)             17    0
#Elixhauser_Renal: Yes vs No(reference)           15    2
#Elixhauser_Rheumatic: Yes vs No(reference)       17    0
#Elixhauser_Tumor: Yes vs No(reference)           17    0
#Elixhauser_Valvular: Yes vs No(reference)        17    0
#Elixhauser_WeightLoss: Yes vs No(reference)      10    7


#these numbers look very consistent and underestimates since not all species in Metaphlan taxonomy are present in Kraken taxonomy and vice versa

#Elixhauser_FluidsLytes which SIG1 species; out of the 7 species in Kraken and 8 in Metaphlan 5 overlap
tax_df_ECIcomponents_SIG1[tax_df_ECIcomponents_SIG1$Comparison == "Elixhauser_FluidsLytes: Yes vs No(reference)" & 
                            tax_df_ECIcomponents_SIG1$Qvalue < 0.1, ]
#[1] "s__Anaerostipes_caccae"   *matches with mp     "s__Bifidobacterium_dentium"        "s__Enterocloster_bolteae"  *matches with mp       
#[5] "s__Enterocloster_clostridioformis" *matches with mp "s__Erysipelatoclostridium_ramosum" "s__Lactobacillus_gasseri"  *matches with mp       
#"s__Limosilactobacillus_fermentum"   *matches with mp


#Elixhauser_Anemia which SIG1 species; out of the 5 species in Kraken and 5 in Metaphlan 4 overlap
tax_df_ECIcomponents_SIG1[tax_df_ECIcomponents_SIG1$Comparison == "Elixhauser_Anemia: Yes vs No(reference)" & 
                            tax_df_ECIcomponents_SIG1$Qvalue < 0.1, "Species"]
#[1] "s__Erysipelatoclostridium_ramosum" "s__Streptococcus_anginosus" *matches with mp     "s__Streptococcus_salivarius"   *matches with mp    
#"s__Veillonella_atypica"   *matches with mp  #[5] "s__Veillonella_parvula"   *matches with mp


#Elixhauser_WeightLoss which SIG1 species; out of the 7 species in Kraken and 2 in Metaphlan 1 overlap
tax_df_ECIcomponents_SIG1[tax_df_ECIcomponents_SIG1$Comparison == "Elixhauser_WeightLoss: Yes vs No(reference)" & 
                            tax_df_ECIcomponents_SIG1$Qvalue < 0.1, "Species"]
#[1] "s__Bifidobacterium_dentium"        "s__Enterocloster_clostridioformis" "s__Erysipelatoclostridium_ramosum" "s__Lacticaseibacillus_paracasei"  
#[5] "s__Lactobacillus_gasseri" *matches with mp   "s__Veillonella_atypica"            "s__Veillonella_parvula"  


#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#see if SIG2 species are enriched in healthy controls

#Toposcore SIG2 species
SIG2_species <- c(
  "Agathobaculum_butyriciproducens",
  "Anaerobutyricum_hallii",
  "Anaerostipes_hadrus",
  "Anaerotignum_faecicola",
  "Blautia_massiliensis",
  "Blautia_wexlerae",
  "Candidatus_Cibiobacter_qucibialis",
  "Clostridiaceae_bacterium",
  "Clostridiaceae_bacterium_OM08_6BH",
  "Clostridiaceae_unclassified_SGB4769",
  "Clostridiales_bacterium_KLE1615",
  "Clostridiales_unclassified_SGB15145",
  "Clostridium_fessum",
  "Clostridium_sp_AF34_10BH",
  "Clostridium_sp_AM22_11AC",
  "Clostridium_sp_AM33_3",
  "Clostridium_sp_AM49_4BH",
  "Coprobacter_fastidiosus",
  "Coprococcus_comes",
  "Coprococcus_eutactus",
  "Dorea_formicigenerans",
  "Dorea_longicatena",
  "Eubacterium_ramulus",
  "Eubacterium_rectale",
  "Eubacterium_ventriosum",
  "Faecalibacillus_intestinalis",
  "Faecalibacterium_prausnitzii",
  "Faecalibacterium_SGB15346",
  "Firmicutes_bacterium_AF16_15",
  "Gemmiger_formicilis",
  "Lachnospira_eligens",
  "Lachnospira_pectinoschiza",
  "Lachnospira_sp_NSJ_43",
  "Lachnospiraceae_bacterium_OM04_12BH",
  "Lachnospiraceae_bacterium_WCA3_601_WT_6H",
  "Lacrimispora_amygdalina",
  "Mediterraneibacter_butyricigenes",
  "Oscillibacter_sp_ER4",
  "Phocaeicola_massiliensis",
  "Roseburia_hominis",
  "Roseburia_intestinalis",
  "Roseburia_inulinivorans",
  "Roseburia_sp_AF02_12",
  "Ruminococcus_bicirculans",
  "Ruminococcus_lactaris"
)


#add s__
SIG2_species <- paste0("s__", SIG2_species)
length(SIG2_species) #45


#there are the results from the comparisons
filename <- "./Code/Submission/Supplementary tables/Supplementary Table 3.xlsx"

#---------------------
#note called tax_df here but that is not descriptive; these are results from differential abundance comparisons

#1) enriched in all cancers  
#TableS2-1	Species; cancer patients vs healthy controls (reference) while adjusting for confounders.
tax_df_allcancer <- read.xlsx(filename, sheet="TableS3-1", rowNames = F)
dim(tax_df_allcancer) #one row per species

table(SIG2_species %in% tax_df_allcancer$Species)
#FALSE  TRUE  #41 of the 45 species detected
#4    41 

tax_df_allcancer_SIG2 <- tax_df_allcancer[tax_df_allcancer$Species %in% SIG2_species,]
table(tax_df_allcancer_SIG2$EffectSize < 0 & tax_df_allcancer_SIG2$Qvalue < 0.1)
#FALSE  TRUE  #41/41 negatively associated to cancer
#   0    41 

tax_df_allcancer_SIG2

#----------------------------------
#Compare with Kraken results

filename <- "./Code/Submission/Supplementary tables/Supplementary Table 2.xlsx"

tax_df_allcancer <- read.xlsx(filename, sheet="TableS2-1", rowNames = F)
dim(tax_df_allcancer) #one row per species
#1278  11
tax_df_allcancer$Species <- gsub(" ", "_", tax_df_allcancer$Species)

table(SIG2_species %in% tax_df_allcancer$Species)
#FALSE  TRUE 
#30    15 

SIG2_species[SIG2_species %in% tax_df_allcancer$Species]
#[1] "s__Agathobaculum_butyriciproducens" "s__Anaerobutyricum_hallii"          "s__Anaerostipes_hadrus"             "s__Anaerotignum_faecicola"         
#[5] "s__Coprobacter_fastidiosus"         "s__Coprococcus_eutactus"            "s__Dorea_formicigenerans"           "s__Faecalibacillus_intestinalis"   
#[9] "s__Faecalibacterium_prausnitzii"    "s__Gemmiger_formicilis"             "s__Lachnospira_eligens"             "s__Phocaeicola_massiliensis"       
#[13] "s__Roseburia_hominis"               "s__Roseburia_intestinalis"          "s__Roseburia_inulinivorans"  

SIG2_species[!SIG2_species %in% tax_df_allcancer$Species]
#[1] "s__Blautia_massiliensis"                     "s__Blautia_wexlerae"                         "s__Candidatus_Cibiobacter_qucibialis"       
#[4] "s__Clostridiaceae_bacterium"                 "s__Clostridiaceae_bacterium_OM08_6BH"        "s__Clostridiaceae_unclassified_SGB4769"     
#[7] "s__Clostridiales_bacterium_KLE1615"          "s__Clostridiales_unclassified_SGB15145"      "s__Clostridium_fessum"                      
#[10] "s__Clostridium_sp_AF34_10BH"                 "s__Clostridium_sp_AM22_11AC"                 "s__Clostridium_sp_AM33_3"                   
#[13] "s__Clostridium_sp_AM49_4BH"                  "s__Coprococcus_comes"                        "s__Dorea_longicatena"                       
#[16] "s__Eubacterium_ramulus"                      "s__Eubacterium_rectale"                      "s__Eubacterium_ventriosum"                  
#[19] "s__Faecalibacterium_SGB15346"                "s__Firmicutes_bacterium_AF16_15"             "s__Lachnospira_pectinoschiza"               
#[22] "s__Lachnospira_sp_NSJ_43"                    "s__Lachnospiraceae_bacterium_OM04_12BH"      "s__Lachnospiraceae_bacterium_WCA3_601_WT_6H"
#[25] "s__Lacrimispora_amygdalina"                  "s__Mediterraneibacter_butyricigenes"         "s__Oscillibacter_sp_ER4"                    
#[28] "s__Roseburia_sp_AF02_12"                     "s__Ruminococcus_bicirculans"                 "s__Ruminococcus_lactaris" 

tax_df_allcancer_SIG2 <- tax_df_allcancer[tax_df_allcancer$Species %in% SIG2_species,]
table(tax_df_allcancer_SIG2$EffectSize < 0 & tax_df_allcancer_SIG2$Qvalue < 0.1)
#TRUE  #all that are matched are also depleted in cancer
#15


#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#go into more detail about relation between ECI score and Toposcore

#define ECI high and low in the same way as in the other comorbidity analysis on beta diversity done by Lu.

meta_with_topo

#clin_meta <- data.obj2$meta.dat
clin_meta <- meta_with_topo[meta_with_topo$Cancer_Healthy == "Cancer", ]

median(clin_meta$Elix_score, na.rm = T) #3

table(clin_meta$Elix_score > 1, clin_meta$Toposcore)
#           SIG1+ SIG2+ #SIG2+ is health associated
#   FALSE   293   116
#   TRUE    707   248
prop.test(table(clin_meta$Elix_score > 1, clin_meta$Toposcore))
#p-value = 0.396


table(clin_meta$Elix_score > 3, clin_meta$Toposcore) #3 is median
#       SIG1+ SIG2+
#FALSE   590   247 #lower comorbidity patients 247 / (247+590) = 29.51% SIG2+
#TRUE    410   117 ##high comorbidity patients 117 / (117+410) = 22.2% SIG2+


fisher.test(table(clin_meta$Elix_score > 3, clin_meta$Toposcore))
#data:  table(clin_meta$Elix_score > 3, clin_meta$Toposcore)
#p-value = 0.003115
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#  0.5239217 0.8843010
#sample estimates:
#  odds ratio 
#0.6818348 


#A ECI score above median is significantly associated to a larger proportion of 
#disease-associated SIG1+ Toposcore when using the binary SIG2+/SIG1+ classification


#--------------------------------------
# Stacked barplot for Toposcore for each ECI component
elix.names <- c("Elixhauser_CHF","Elixhauser_Arrhythmia","Elixhauser_Valvular","Elixhauser_PHTN","Elixhauser_PVD","Elixhauser_HTN","Elixhauser_Paralysis",  
                "Elixhauser_NeuroOther","Elixhauser_Pulmonary","Elixhauser_DM","Elixhauser_DMcx","Elixhauser_Hypothyroid","Elixhauser_Renal","Elixhauser_Liver",      
                "Elixhauser_PUD" ,"Elixhauser_HIV","Elixhauser_Psychoses", "Elixhauser_Mets",
                "Elixhauser_Lymphoma","Elixhauser_Mets","Elixhauser_Tumor","Elixhauser_Rheumatic","Elixhauser_Coagulopathy",
                "Elixhauser_Obesity","Elixhauser_WeightLoss","Elixhauser_FluidsLytes","Elixhauser_BloodLoss","Elixhauser_Anemia","Elixhauser_Alcohol","Elixhauser_Drugs",
                "Elixhauser_Depression")

elix.names <- elix.names[!elix.names %in% c("Elixhauser_Tumor", "Elixhauser_Lymphoma", "Elixhauser_Obesity",
                                            "Elixhauser_Paralysis", "Elixhauser_HIV",
                                            "Elixhauser_BloodLoss", "Elixhauser_Psychoses")]


#make barplot using SIG_class
meta_with_topo <- clin_meta

#just get the Yes row and name it after the elix.names, also get the pairwise split
#the stats should test whether Toposcore is enriched among the patients with the comorbidity component vs the patients without

res_list <- list()
res_list_split <- list()
pval_list <- list()
for (i in 1:length(elix.names)) {
  tab <- table(meta_with_topo[,elix.names[i]], meta_with_topo$SIG_class)
  res_list[[i]] <- tab["Yes",]
  tab_binary <- table(meta_with_topo[,elix.names[i]], meta_with_topo$Toposcore)
  res_list_split[[i]] <- tab_binary["Yes",]
  pval_list[[i]] <- fisher.test(tab_binary)$p.value
}

names(res_list) <- elix.names
names(res_list_split) <- elix.names
names(pval_list) <- elix.names


tab <- do.call(rbind, res_list)
tab_binary <- do.call(rbind, res_list_split)
p_df <- do.call(rbind, pval_list)


tab_perc <- apply(tab, 1, function(x) 100 * x / sum(x, na.rm = TRUE))
tab_perc <- tab_perc[c("SIG2","Gray","SIG1"),]
tab_perc <- tab_perc[,rev(order(tab_perc["SIG1",]))]

#add extra bars in to indicate the SIG1+ or SIG2+ binary split
tab_binary <- tab_binary[colnames(tab_perc), , drop = FALSE]

tab_binary_perc <- apply(tab_binary, 1, function(x) 100 * x / sum(x, na.rm = TRUE))
tab_binary_perc <- t(tab_binary_perc)  # rows = cancers, cols = levels of Toposcore
line_level <- "SIG2+"
line_vals <- tab_binary_perc[, line_level]

dot_vals <- line_vals[colnames(tab_perc)]


qvals <- p.adjust(p_df[,1], method = "fdr")
qvals <- qvals[colnames(tab_perc)]

stars <- qvals

stars[which(qvals < 0.1)] <- "*"
stars[which(qvals < 0.01)] <- "**"
stars[which(qvals < 0.001)] <- "***"
stars[which(qvals > 0.1)] <- "" 



par(mar=c(5,16,1,3), xpd = TRUE)
bp <- barplot(tab_perc,
              horiz  = TRUE,  
              col = c("#0072B2", "darkgrey", "#CC79A7"),
              border = "white",
              space = 0.2,
              xlab = "Toposcore by comorbidity component", las=1)

points(x = dot_vals, y = bp,
       pch = 18, col = "#0072B2", cex = 1.5)

legend("topleft",
       inset  = c(-0.9, 0.02),  # tweak this value
       legend = rownames(tab_perc),
       fill   = c("#0072B2", "darkgrey", "#CC79A7"),
       bty    = "n")

legend("topleft",
       inset  = c(-0.885, 0.02),
       legend = c("", "", "", "SIG1+/SIG2+ split"), 
       pch    = c(NA, NA, NA, 18),
       col    = c(NA, NA, NA, "#0072B2"),
       pt.cex = c(NA, NA, NA, 2),
       bty    = "n")

# add text a bit to the right of each bar
text(x     = 105,                
     y     = bp,                         
     labels = paste0(stars),
     xpd   = TRUE,
     cex   = 0.8, adj = 0)



tab_perc

