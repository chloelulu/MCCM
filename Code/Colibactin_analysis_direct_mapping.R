require(tidyverse)
require(readxl)
require(openxlsx)
require(corrr)
require(reshape2)
require(ggbreak)
library(scales)
library(forcats)



file_dir <- dirname(rstudioapi::getSourceEditorContext()$path)
file_dir <- dirname(dirname(dirname(dirname(file_dir))))
setwd(file_dir) #"mforge_clean"


#-------------------------------------------------------------------------------


clin_meta <- as.data.frame(read_excel("Code/Submission/Supplementary tables/Supplementary Table 1.xlsx", sheet = "TableS1-2"))
colib_direct_sub <- as.data.frame(read_excel("Code/Submission/Supplementary tables/Supplementary Table 1.xlsx", sheet = "TableS1-6"))

rownames(colib_direct_sub) <- colib_direct_sub$sample
rownames(clin_meta) <- clin_meta$sample

colib_direct_sub <- colib_direct_sub[rownames(clin_meta),]
clin_meta$direct_mapping_positive <- colib_direct_sub[rownames(clin_meta), "positive"]

clin_meta$trunc_cancer_class_names <- clin_meta$`icd10_first_3_name(name used in manuscript)`
cutoff_cancer_class_names <- names(table(clin_meta$trunc_cancer_class_names))[which(table(clin_meta$trunc_cancer_class_names) > 15)]

dim(clin_meta)
dim(colib_direct_sub)

clin_meta$cancer_healthy <- clin_meta$Group
clin_meta$cancer_healthy[clin_meta$cancer_healthy == "Control"] <- "Healthy"


#---------------------------------------

table(clin_meta$direct_mapping_positive) #236

#all cancer
model <- glm(direct_mapping_positive ~ cancer_healthy + Age, data = clin_meta, family = binomial)
summary(model) #significant 0.0141


tab <- with(clin_meta, table(direct_mapping_positive, cancer_healthy))
tab
#           cancer_healthy
# direct_mapping_positive  Cancer Healthy
#                    FALSE   1155     260
#                    TRUE     209      27

tab_perc <- prop.table(t(tab), margin = 1) * 100
tab_perc <- t(tab_perc)
#direct_mapping_positive    Cancer   Healthy
#FALSE                      84.677419 90.592334
#TRUE                       15.322581  9.407666

ref_row <- "Healthy"

pvals <- fisher.test(tab)$p.value #0.0091

tab_perc <- tab_perc[,c("Healthy", "Cancer")]


par(mar = c(6, 4, 3, 1), xpd = TRUE)
bp <- barplot(tab_perc,
              horiz  = F,
              col    = c("white", "#CC79A7"),
              border = "black",
              space  = 0.2,
              ylab   = "colibactin positivity (%)", las = 1,
              names.arg = c("", ""))

my_names <- colnames(tab_perc)  
text(x   = bp,
     y   = -5,                   
     labels = my_names,
     xpd = TRUE,
     cex = 1)

legend("bottomleft",
       inset  = c(-0.6, -0.3),
       legend = c("colibactin negative", "colibactin positive"),
       fill   = c("white", "#CC79A7"),
       bty    = "n")

text(y      = 105,
     x      = mean(bp)*1.45,
     labels = paste0("p = ", round(pvals, digits = 3)),
     xpd    = TRUE,
     cex    = 1,
     adj    = 1)


#---------------------------------------
#cancer classes
clin_meta_sub <- clin_meta[clin_meta$trunc_cancer_class_names %in% cutoff_cancer_class_names,]
clin_meta_sub$trunc_cancer_class_names <- relevel(as.factor(clin_meta_sub$trunc_cancer_class_names),ref = "healthy")


model <- glm(direct_mapping_positive ~ trunc_cancer_class_names + Age, data = clin_meta_sub, family = binomial)
model_res <- as.data.frame(summary(model)$coefficients)
#correct for multiple testing
model_res$FDR <- p.adjust(as.numeric(model_res[,"Pr(>|z|)"]), method = "fdr")
model_res[model_res$FDR < 0.1,]


#get prevalence by cancer class from direct mapping
prev_table <- table(clin_meta_sub$direct_mapping_positive, clin_meta_sub$trunc_cancer_class_names)
prop.test(t(prev_table[,c("healthy", "esophagus")]))
#4.294e-05

temp_table <- t(table(clin_meta_sub$direct_mapping_positive, clin_meta_sub$trunc_cancer_class_names))
colibactin_df <- as.data.frame(cbind(temp_table[,"FALSE"], temp_table[,"TRUE"], temp_table[,"FALSE"] / rowSums(temp_table), temp_table[,"TRUE"] / rowSums(temp_table)))
names(colibactin_df) <- c("n_neg", "n_pos", "perc_neg", "perc_pos")


colibactin_df <- cbind(colibactin_df, p=model_res$`Pr(>|z|)`[1:23], q=model_res$FDR[1:23]) #select first 23 since the last one is Age
colibactin_df$p[1] <- 1
colibactin_df$q[1] <- 1

colibactin_df_sort <- colibactin_df[order(colibactin_df$perc_pos, decreasing = F),]
colibactin_df_sort_plot <- t(colibactin_df_sort[,c("perc_neg", "perc_pos")]) * 100


colibactin_df_sort$signif <- ""
colibactin_df_sort$signif[colibactin_df_sort$q < 0.01] <- "***"
colibactin_df_sort$signif[colibactin_df_sort$q >= 0.01 & colibactin_df_sort$q < 0.05] <- "**"
colibactin_df_sort$signif[colibactin_df_sort$q >= 0.05 & colibactin_df_sort$q < 0.1] <- "*"

colibactin_df_sort
#                             n_neg n_pos  perc_neg   perc_pos    p            q signif
#liver intrahepatic bile ducts    38     1 0.9743590 0.02564103 2.613093e-01 0.3919639469       
#kidney except renal pelvis       28     2 0.9333333 0.06666667 8.332924e-01 0.9143846221       
#brain                            34     3 0.9189189 0.08108108 9.222477e-01 0.9222476851       
#bladder                          44     4 0.9166667 0.08333333 8.381859e-01 0.9143846221       
#unspecified skin                 22     2 0.9166667 0.08333333 9.065868e-01 0.9222476851       
#healthy                         260    27 0.9059233 0.09407666 1.000000e+00 1.0000000000       
#breast                          212    23 0.9021277 0.09787234 5.356560e-01 0.6427871672       
#melanoma of skin                 42     5 0.8936170 0.10638298 5.288925e-01 0.6427871672       
#prostate                        104    15 0.8739496 0.12605042 1.513032e-01 0.2826660356       
#base of tongue                   18     3 0.8571429 0.14285714 3.120523e-01 0.4405444497       
#pancreas                         41     8 0.8367347 0.16326531 7.323383e-02 0.1597829047      
#ovary                            34     7 0.8292683 0.17073171 7.245861e-02 0.1597829047       
#lymphoid leukemia                14     3 0.8235294 0.17647059 1.840634e-01 0.2945015009       
#oropharynx                       14     3 0.8235294 0.17647059 1.776879e-01 0.2945015009       
#bronchus lung                    74    16 0.8222222 0.17777778 1.392035e-02 0.0556814198      *
#other connective soft tissue     37     8 0.8222222 0.17777778 5.053634e-02 0.1347635746       
#colorectal                       54    12 0.8181818 0.18181818 2.059284e-02 0.0617785308      *
#stomach                          13     3 0.8125000 0.18750000 1.531108e-01 0.2826660356     
#multiple myeloma plasma cell     48    12 0.8000000 0.20000000 9.396960e-03 0.0451054091     **
#corpus uteri                     25     7 0.7812500 0.21875000 1.694324e-02 0.0580910994      *
#tonsil                           21     7 0.7500000 0.25000000 7.436397e-03 0.0446183807     **
#neuroendocrine tumors            14     7 0.6666667 0.33333333 8.533664e-04 0.0068269315    ***
#esophagus                        22    12 0.6470588 0.35294118 3.364491e-05 0.0004037389    ***


#---------------------------------------
#make the figure

clin_meta_sub$colibactin_binary <- ifelse(clin_meta_sub$direct_mapping_positive, 1, 0)

temp_table <- t(table(clin_meta_sub$colibactin_binary, clin_meta_sub$trunc_cancer_class_names))
colibactin_df <- as.data.frame(cbind(temp_table[,"0"], temp_table[,"1"], temp_table[,"0"] / rowSums(temp_table), temp_table[,"1"] / rowSums(temp_table)))
names(colibactin_df) <- c("n_neg", "n_pos", "perc_neg", "perc_pos")


colibactin_df <- cbind(colibactin_df, p=model_res$`Pr(>|z|)`[1:23], q=model_res$FDR[1:23]) #select first 23 since the last one is Age
colibactin_df$p[1] <- 1
colibactin_df$q[1] <- 1

colibactin_df_sort <- colibactin_df[order(colibactin_df$perc_pos, decreasing = F),]
colibactin_df_sort_plot <- t(colibactin_df_sort[,c("perc_neg", "perc_pos")]) * 100


colibactin_df_sort$signif <- ""
colibactin_df_sort$signif[colibactin_df_sort$q < 0.01] <- "***"
colibactin_df_sort$signif[colibactin_df_sort$q >= 0.01 & colibactin_df_sort$q < 0.05] <- "**"
colibactin_df_sort$signif[colibactin_df_sort$q >= 0.05 & colibactin_df_sort$q < 0.1] <- "*"

three_stars <- colibactin_df_sort$signif == "***"
two_stars <- colibactin_df_sort$signif == "**"
one_star <- colibactin_df_sort$signif == "*"


#save as svg to edit
#png("./Figure/RM_figures/colibactin_cancer_classes_vs_control.png", width=5.75*300, height=5.75*300, res=300)
svg("./Figure/RM_figures/colibactin_cancer_classes_vs_control_revision.svg", width=8, height=6)

par(mar = c(6, 14, 0.2, 1), las = 1, cex.axis = 1.05, cex.lab = 1.05) 
# Create a stacked barplot
bp <- barplot(colibactin_df_sort_plot,
              beside = FALSE,
              col = c("white", "#CC79A7"), 
              xlab = "Colibactin gene prevalence %",
              ylab = "",
              main = "", horiz = T,
              space=0, xlim=c(0,100), las=1)
#add the indication of significance
text(x = rep(95, length(one_star)), y = bp[one_star], labels = "*", cex = 1.5, col = "black")
text(x = rep(95, length(two_stars)), y = bp[two_stars], labels = "**", cex = 1.5, col = "black")
text(x = rep(95, length(three_stars)), y = bp[three_stars], labels = "***", cex = 1.5, col = "black")
#manually cut out a chunk
plotrix::axis.break(axis = 1,           # x-axis (bottom)
                    breakpos = 10,      # position of break symbol
                    style = "zigzag",   # style of break
                    brw = 0.02)
legend(
  x=-40,
  y=-1,
  legend = c("colibactin negative", "colibactin positive", "*  q < 0.1", "** q < 0.05", "*** q < 0.01"),
  pch = c(15, 15, NA, NA, NA),
  col = c("white", "#CC79A7", "black", "black", "black"),
  bty = "n",
  pt.cex = 2,
  horiz = FALSE,
  xpd = TRUE                 
)

dev.off()


#--------------------------------------
#Positivity in the reprocessed studies

#Jeong stomach 37.5% https://pubmed.ncbi.nlm.nih.gov/38658841/
#Kartal pancreatic 12.0% https://pubmed.ncbi.nlm.nih.gov/35260444/
#Yachida CRC 28.8% https://pubmed.ncbi.nlm.nih.gov/31171880/


#-------------------------------------------------------------------------------
#early onset and colibactin
#first inspect relationship with Age


clin_meta_cancer <- clin_meta[clin_meta$cancer_healthy == "Cancer",]

boxplot(clin_meta_cancer$Age ~ clin_meta_cancer$direct_mapping_positive)
wilcox.test(clin_meta_cancer$Age ~ clin_meta_cancer$direct_mapping_positive) #p-value = 0.4592


clin_meta_cancer$early_onset <- "-" #only relevant for cancer
clin_meta_cancer$early_onset[clin_meta_cancer$Age > 50] <- "normal"
clin_meta_cancer$early_onset[clin_meta_cancer$Age <= 50] <- "early" 

table(clin_meta_cancer$early_onset)
#early normal 
#225   1139 

table(clin_meta_cancer$direct_mapping_positive, clin_meta_cancer$early_onset) #no relationship with age or early onset
#       early normal
#FALSE   194    961
#TRUE     31    178
prop.test(table(clin_meta_cancer$direct_mapping_positive, clin_meta_cancer$early_onset))
#p-value = 0.5467


#test only on CRC early onset
clin_meta_CRC <- clin_meta_cancer[grep("colorectal", clin_meta_cancer$trunc_cancer_class_names),]
table(clin_meta_CRC$early_onset, clin_meta_CRC$direct_mapping_positive)
#         FALSE TRUE
#early     19    2
#normal    35   10

prop.test(table(clin_meta_CRC$early_onset, clin_meta_CRC$direct_mapping_positive))
#p-value = 0.3664
#prop 1    prop 2 
#0.9047619 0.7777778 

#2/21 * 100     #9.5% trending lower in early onset
#10/45 * 100    #22.2%


tab <- with(clin_meta_CRC, table(direct_mapping_positive, early_onset))
tab
#                         early_onset
#direct_mapping_positive early normal
#FALSE                      19     35
#TRUE                        2     10


tab_perc <- prop.table(t(tab), margin = 1) * 100
tab_perc <- t(tab_perc)

ref_row <- "early"

pvals <- fisher.test(tab)$p.value

tab_perc <- tab_perc[,c("early", "normal")]



par(mar = c(6, 4, 3, 1), xpd = TRUE)
bp <- barplot(tab_perc,
              horiz  = F,
              col    = c("white", "#CC79A7"),
              border = "black",
              space  = 0.2,
              ylab   = "colibactin positivity (%)", las = 1,
              names.arg = c("", ""))

my_names <- colnames(tab_perc)  
text(x   = bp,
     y   = -5,                   
     labels = my_names,
     xpd = TRUE,
     cex = 1)

legend("bottomleft",
       inset  = c(-0.4, -0.4),
       legend = c("colibactin negative", "colibactin positive"),
       fill   = c("white", "#CC79A7"),
       bty    = "n")

text(y      = 105,
     x      = mean(bp)*1.3,
     labels = paste0("p = ", round(pvals, digits = 3)),
     xpd    = TRUE,
     cex    = 1,
     adj    = 1)



#--------------------------------------


## Healthy contingency table direct mapping
##            MCCM  Yachida
## positive     27      60
## negative    260     191

direct_mat <- matrix(
  c(27, 260,
    60, 191),
  nrow = 2,
  byrow = TRUE,
  dimnames = list(
    c("MCCM", "Yachida"),
    c("positive", "negative")
  )
)

direct_mat

## Two-sample test for equality of proportions (positive vs total)
prop.test(
  x = direct_mat[ , "positive"],
  n = rowSums(direct_mat),
  correct = FALSE
)


#--------------------------------------


tax_df <- as.data.frame(read_excel("Code/Submission/Supplementary tables/Supplementary Table 1.xlsx", sheet = "TableS1-3"), row.names = T)
rownames(clin_meta) <- clin_meta$SampleID

rownames(tax_df) <- tax_df$Species
tax_df <- tax_df[,2:ncol(tax_df)]
tax_df <- tax_df[,rownames(clin_meta)]

dim(tax_df)
dim(colib_direct_sub)
dim(clin_meta)

Ecoli_names <- rownames(tax_df)[grep("Escherichia coli", rownames(tax_df))]

#some are very prevalent; correlate counts with number of genes detected?
as.numeric(tax_df[Ecoli_names[1],])
table(as.numeric(tax_df[Ecoli_names[1],]) > 0)

as.numeric(tax_df[Ecoli_names[2],])
table(as.numeric(tax_df[Ecoli_names[2],]) > 0)

as.numeric(tax_df[Ecoli_names[3],])
table(as.numeric(tax_df[Ecoli_names[3],]) > 0)

as.numeric(tax_df[Ecoli_names[4],])
table(as.numeric(tax_df[Ecoli_names[4],]) > 0)


cor.test(colib_direct_sub$genes_detected, as.numeric(tax_df[Ecoli_names[1],]))
#p-value < 2.2e-16  cor 0.2512584

cor.test(colib_direct_sub$genes_detected, as.numeric(tax_df[Ecoli_names[2],]))
#p-value = 0.6176 cor 0.01229663

cor.test(colib_direct_sub$genes_detected, as.numeric(tax_df[Ecoli_names[3],]))
#p-value = 0.4999 cor 0.01661651

cor.test(colib_direct_sub$genes_detected, as.numeric(tax_df[Ecoli_names[4],]))
#p-value = 0.4507  cor 0.01857451 


#now check with all species
# vector of rownames to test; if you literally want all:
taxa <- rownames(tax_df)

cor_list <- lapply(taxa, function(tax) {
  x <- colib_direct_sub$genes_detected
  y <- as.numeric(tax_df[tax, ])
  
  # optional: ensure same length and remove NAs
  ok <- complete.cases(x, y)
  x <- x[ok]
  y <- y[ok]
  
  ct <- suppressWarnings(cor.test(x, y, method = "pearson"))
  data.frame(
    taxon    = tax,
    cor      = unname(ct$estimate),
    p_value  = ct$p.value,
    stringsAsFactors = FALSE
  )
})
cor_df <- do.call(rbind, cor_list)
rownames(cor_df) <- NULL

cor_df$q_value <- p.adjust(cor_df$p_value, method = "fdr")

hist(cor_df$cor)
range(cor_df$cor)

head(cor_df[order(cor_df$cor, decreasing = T),])
#                               taxon       cor      p_value      q_value
#2811               s__Escherichia coli 0.2512584 3.456547e-25 2.709588e-21
#3677      s__Lactobacillus kimbladii_B 0.2426804 1.476317e-23 5.786426e-20
#1270             s__CAG-41 sp001941225 0.2119444 3.193198e-18 8.343828e-15


#prop test binary
tab <- table(as.numeric(tax_df[Ecoli_names[1],]) > 0, colib_direct_sub$positive)
tab
prop.test(tab)




