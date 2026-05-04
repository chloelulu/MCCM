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



#colibactin genes are summed up here 202507_subsetting_functional_data_MCCM.R
filename <- "Data/colibactin_genefamilies_counts.xlsx"
colibactin_df_counts_MCCM <- read.xlsx(filename, sheet=1, rowNames = T)
dim(colibactin_df_counts_MCCM)

colibactin_df_counts_MCCM[,1:4]

rownames(colibactin_df_counts_MCCM)
#"UniRef90_A0A2T3SSR3"    "UniRef90_Q0P7J7"        "UniRef90_Q0P7J8"        "UniRef90_Q0P7K4"        "UniRef90_Q0P7K7"        "UniRef90_Q0P7K8"       
#"UniRef90_UPI00050AA6A1" "UniRef90_UPI00069CAA6D"


#load colibactin annotation
genefamily_ann <- as.data.frame(read_tsv("Data/colibactin_genefamilies.txt", col_names = c("ID", "description")))
rownames(genefamily_ann) <- genefamily_ann$ID

genefamily_ann[rownames(colibactin_df_counts_MCCM),]
#                                   ID                                                                        description
#UniRef90_A0A2T3SSR3       UniRef90_A0A2T3SSR3                                   Colibactin non-ribosomal peptide synthetase ClbJ
#UniRef90_Q0P7J7               UniRef90_Q0P7J7                                       Colibactin biosynthesis acyltransferase ClbG
#UniRef90_Q0P7J8               UniRef90_Q0P7J8                                   Colibactin non-ribosomal peptide synthetase ClbH
#UniRef90_Q0P7K4               UniRef90_Q0P7K4                                   Colibactin non-ribosomal peptide synthetase ClbN
#UniRef90_Q0P7K7               UniRef90_Q0P7K7                                          Colibactin biosynthesis thioesterase ClbQ <- tailoring/resistance
#UniRef90_Q0P7K8               UniRef90_Q0P7K8                                            Colibactin self-protection protein ClbS <- tailoring/resistance
#UniRef90_UPI00050AA6A1 UniRef90_UPI00050AA6A1 colibactin hybrid non-ribosomal peptide synthetase/type I polyketide synthase ClbK <- core megasynthase
#UniRef90_UPI00069CAA6D UniRef90_UPI00069CAA6D colibactin hybrid non-ribosomal peptide synthetase/type I polyketide synthase ClbB <- core megasynthase


#explore this as well as counts analysis:
#For more stringent calls, you can require co‑detection of at least one core megasynthase (ClbB or ClbK) plus one tailoring/resistance 
#gene (ClbQ or ClbS), which aligns well with how culture and PCR‑based studies define colibactin‑producing strains.



load(file = "Data/data.obj.raw.core.RData") 
clin_meta <- as.data.frame(data.obj$meta.dat)
clin_meta$X.SampleID


#-------------------------------------------------------------------------------
# Create correlation data frame of rows colibactin

row_correlations <- correlate(t(colibactin_df_counts_MCCM), method = "spearman")

row_cor_dt <- data.table::as.data.table(row_correlations, keep.rownames = "Row1")
long_corr <- melt(row_cor_dt, varnames = c("Row1", "Row2"), value.name = "Correlation")

#if one is found all seem present which is expected when looking at an operon 
hist(long_corr$Correlation, xlim=c(0,1), breaks=5)

median(long_corr$Correlation, na.rm = T) #0.89 is the median


#------------------------------------
#check correlation between ClbQ + ClbS, and ClbK + ClbB. Both strong

colibactin_df_counts_MCCM[1:8,1:4]

genes_tailoring_resistance <- c("UniRef90_Q0P7K7", "UniRef90_Q0P7K8")
genes_core_megasynthase <- c("UniRef90_UPI00050AA6A1", "UniRef90_UPI00069CAA6D")

cor(t(colibactin_df_counts_MCCM[genes_tailoring_resistance,]), method = "spearman") #0.84
cor(t(colibactin_df_counts_MCCM[genes_core_megasynthase,]), method = "spearman") #0.91


#-------------------------------------------------------------------------------
#Sum the abundances to get one vector and associate with cancer

colibactin_genes <- colSums(colibactin_df_counts_MCCM)
hist(colibactin_genes, breaks=20, ylim=c(0,20))
table(colibactin_genes > 0)


clin_meta$colibactin_counts <- colibactin_genes

clin_meta$colibactin_binary <- 0
clin_meta$colibactin_binary[which(colSums(colibactin_df_counts_MCCM > 0) > 2)] <- 1
table(clin_meta$colibactin_binary)
# 0    1 
#1523  128 

table(clin_meta$colibactin_counts > 0)
#FALSE  TRUE 
#1492   159 

159 / nrow(clin_meta) * 100 #9.63% prevalence without any filtering
128 / nrow(clin_meta) * 100 
#7.75% prevalence


hist(colSums(colibactin_df_counts_MCCM > 0), ylim=c(0, 100))


#------------------------------------
#now do the prevalence by sample that have one of genes_tailoring_resistance AND genes_core_megasynthase

table(colSums(colibactin_df_counts_MCCM[genes_tailoring_resistance, ]) > 0)
table(colSums(colibactin_df_counts_MCCM[genes_core_megasynthase, ]) > 0)

#106 positive, more stringent
table(colSums(colibactin_df_counts_MCCM[genes_tailoring_resistance, ]) > 0 & colSums(colibactin_df_counts_MCCM[genes_core_megasynthase, ]) > 0)

cols_keep <- colSums(colibactin_df_counts_MCCM[genes_tailoring_resistance, ]) > 0 &
  colSums(colibactin_df_counts_MCCM[genes_core_megasynthase, ]) > 0

true_sample_names <- colnames(colibactin_df_counts_MCCM)[cols_keep]


clin_meta$colibactin_binary_strict <- 0
clin_meta[true_sample_names, "colibactin_binary_strict"] <- 1

table(clin_meta$colibactin_binary_strict)
#  0    1 
#1545  106 


#--------------------------------------
#despite the normalization by seq depth explore potential correlation with sequencing depth

#no correlation
plot(clin_meta$depth, clin_meta$colibactin_counts)
cor.test(clin_meta$depth, clin_meta$colibactin_counts)


#-------------------------------------------------------------------------------

clin_meta$icd10_first_3_name[is.na(clin_meta$icd10_first_3_name)] <- "Healthy"
clin_meta$cancer_healthy <- "Cancer"
clin_meta$cancer_healthy[clin_meta$icd10_first_3_name == "Healthy"] <- "Healthy"

clin_meta$cancer_healthy <- as.factor(clin_meta$cancer_healthy)
clin_meta$cancer_healthy <- relevel(as.factor(clin_meta$cancer_healthy),ref = "Healthy")


clin_meta$trunc_cancer_class_names <- clin_meta$icd10_first_3_name
clin_meta$trunc_cancer_class_names <- tolower(gsub("Malignant neoplasm of |Malignant |malignant |Other and | neoplasm| neoplasm of| neoplasms|neoplasm of |\\,|and ", "", clin_meta$trunc_cancer_class_names))
clin_meta$trunc_cancer_class_names[clin_meta$trunc_cancer_class_names %in% c("colon", "rectum")] <- "colorectal"

cutoff_cancer_class_names <- names(table(clin_meta$trunc_cancer_class_names))[which(table(clin_meta$trunc_cancer_class_names) > 15)]


#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#loading the results from the direct mapping approach

ilename <- "Data/direct mapping results/all_targets.positivity_direct_mapping.tsv"
colib_direct <- read.delim(filename, header = TRUE, sep = "\t",
                           check.names = FALSE, stringsAsFactors = FALSE)
colib_direct <- colib_direct[colib_direct$type == "operon",]

#remove L003 appendix
colib_direct$sample <- gsub("_L003", "", colib_direct$sample)
#make names short
colib_direct$BIOMEid <- gsub("_.*", "", colib_direct$sample)
#remove V2 appendix
clin_meta$X.SampleID <- gsub("_V2", "", clin_meta$X.SampleID)


#select healthy by short BIOME and Cancer by sample
comb_ind <- c(which(colib_direct$sample %in% clin_meta$X.SampleID), 
              which(colib_direct$BIOMEid %in% rownames(clin_meta)[clin_meta$cancer_healthy == "Healthy"]))

colib_direct_sub <- colib_direct[comb_ind,]

#some duplicates; because of ICI?
dim(colib_direct_sub) #1651   11

rownames(colib_direct_sub) <- colib_direct_sub$BIOMEid

table(rownames(clin_meta) %in% colib_direct_sub$BIOMEid)
#TRUE 
#1651 

table(colib_direct_sub$positive)


#---------------------------------------
#write this to file for Supplementary Table S1-6

write.csv(colib_direct_sub, "./Code/Submission/Supplementary tables/colib_direct_sub_TableS1_6.csv")
#used cutoff is more than 8 genes



#---------------------------------------



#some duplicates; because of ICI
dim(colib_direct_sub)

rownames(colib_direct_sub) <- colib_direct_sub$BIOMEid

table(rownames(clin_meta) %in% colib_direct_sub$BIOMEid)
#TRUE 
#1651 


clin_meta$direct_mapping_positive <- colib_direct_sub[rownames(clin_meta), "positive"]

table(clin_meta$colibactin_binary) #128
table(clin_meta$direct_mapping_positive) #236

table(clin_meta$colibactin_binary, clin_meta$direct_mapping_positive)
#   FALSE TRUE
#0  1415  108
#1     0  128


#all cancer
model <- glm(direct_mapping_positive ~ cancer_healthy + Age, data = clin_meta, family = binomial)
summary(model) #significant 0.0141


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
#pancreas                         41     8 0.8367347 0.16326531 7.323383e-02 0.1597829047       --> reprocessed and 12% in that study. Also ns here
#ovary                            34     7 0.8292683 0.17073171 7.245861e-02 0.1597829047       
#lymphoid leukemia                14     3 0.8235294 0.17647059 1.840634e-01 0.2945015009       
#oropharynx                       14     3 0.8235294 0.17647059 1.776879e-01 0.2945015009       
#bronchus lung                    74    16 0.8222222 0.17777778 1.392035e-02 0.0556814198      *
#other connective soft tissue     37     8 0.8222222 0.17777778 5.053634e-02 0.1347635746       
#colorectal                       54    12 0.8181818 0.18181818 2.059284e-02 0.0617785308      *
#stomach                          13     3 0.8125000 0.18750000 1.531108e-01 0.2826660356       --> reprocessed and 37.5% prevalent in another study. Higher here but ns due to small n. (n=16) Jeong has 32 subjects
#multiple myeloma plasma cell     48    12 0.8000000 0.20000000 9.396960e-03 0.0451054091     **
#corpus uteri                     25     7 0.7812500 0.21875000 1.694324e-02 0.0580910994      *
#tonsil                           21     7 0.7500000 0.25000000 7.436397e-03 0.0446183807     **
#neuroendocrine tumors            14     7 0.6666667 0.33333333 8.533664e-04 0.0068269315    ***
#esophagus                        22    12 0.6470588 0.35294118 3.364491e-05 0.0004037389    ***

#another fascinating thing is that neuroendocrine and liver are both very outlying in the PCOA but on opposite spectra of the colibactin here; 
#some nuance on non-specificity of colibactin message.



#-------------------------------------------------------------------------------
#logistic regression 

model <- glm(colibactin_binary ~ cancer_healthy + Age, data = clin_meta, family = binomial)
summary(model) #significant 0.00317



#now stricter, same or more significant
model_strict <- glm(colibactin_binary_strict ~ cancer_healthy + Age, data = clin_meta, family = binomial)
summary(model_strict) #significant 0.00242



#model on binary but here plot counts
boxplot(colibactin_counts ~ cancer_healthy, data = clin_meta)

#also analyze as prop.test
prop.test(t(table(clin_meta$colibactin_binary, clin_meta$cancer_healthy)))
#p-value = 0.0008402
#   prop 1    prop 2 
#0.9721254 0.9120235 

#(1 - 0.9721254) * 100 #2.8% prevalence
#(1 - 0.9120235) * 100 #8.8% prevalence



#make boxplot healthy vs cancer colibactin; get code from the co-abundance network code file


sel_colors <- RColorBrewer::brewer.pal(8,"Accent")
#sel_colors_plot <- RColorBrewer::brewer.pal(12,"Paired")
sel_colors_plot <- c("#0072B2", "#CC79A7")




#make higher and simplify y axis label

clin_meta$cancer_healthy <- tolower(clin_meta$cancer_healthy)


set.seed(42)

colibactin_plot <- ggplot(clin_meta, aes(x = cancer_healthy,
                                         y = colibactin_counts,
                                         fill = cancer_healthy)) +
  geom_boxplot(outlier.shape = NA) +                 # hide default outliers
  geom_jitter(aes(color = cancer_healthy),
              width = 0.15, size = 2, alpha = 0.8) + # jittered points
  scale_fill_manual(values = sel_colors_plot) +
  scale_color_manual(values = sel_colors_plot) +
  geom_segment(aes(y = 800, yend = 800, x = 1, xend = 2), linewidth = 0.3) +
  geom_text(aes(y = 800*1.03, x = 1.5, label = "p = 0.0032"), size = 3, vjust = 0) +
  theme_classic() +                                  # axis lines come back
  theme(
    panel.border = element_rect(color = "black", fill = NA, size = 0.5), 
    axis.ticks.y = element_line(color = "black", size = 0.3),
    axis.ticks.x = element_line(color = "black", size = 0.),
    axis.line.x  = element_line(color = "black", size = 0.1),
    axis.line.y  = element_line(color = "black", size = 0.1),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "none",
    plot.title  = element_text(size = 10, hjust = 0.5, family = "Arial", colour = "black"),
    axis.text   = element_text(size = 9, family = "Arial", colour = "black"),
    axis.title  = element_text(size = 10, family = "Arial", colour = "black"),
    plot.margin = grid::unit(c(1,1,1,1), "lines")
  ) +
  labs(
    y = "colibactin gene counts (summed)",
    x = NULL                                         # suppress x-axis title
  )


ggsave("./Figure/RM_figures/boxplot_colibactin_cancer_vs_control.png", colibactin_plot, width = 2.25, height = 3.5, units = "in", dpi = 300)



#--------------------------------------



model <- glm(colibactin_binary ~ trunc_cancer_class_names + Age, data = clin_meta_sub, family = binomial)
model_res <- as.data.frame(summary(model)$coefficients)
#correct for multiple testing
model_res$FDR <- p.adjust(as.numeric(model_res[,"Pr(>|z|)"]), method = "fdr")
model_res[model_res$FDR < 0.1,]

#significant ones
#                                                       Estimate Std. Error   z value     Pr(>|z|)          FDR
#(Intercept)                                          -3.199696  0.5547627 -5.767684 8.036836e-09 1.928841e-07
#trunc_cancer_class_namesbladder                       1.339847  0.6958940  1.925361 5.418422e-02 9.351502e-02
#trunc_cancer_class_namescolorectal                    1.504835  0.5584826  2.694506 7.049307e-03 3.383667e-02
#trunc_cancer_class_namescorpus uteri                  1.426773  0.7421676  1.922440 5.455043e-02 9.351502e-02
#trunc_cancer_class_namesesophagus                     2.533128  0.5978954  4.236742 2.267863e-05 2.721436e-04
#trunc_cancer_class_nameslymphoid leukemia             1.646538  0.8548633  1.926083 5.409395e-02 9.351502e-02
#trunc_cancer_class_namesmelanoma of skin              1.307758  0.6695253  1.953261 5.078865e-02 9.351502e-02
#trunc_cancer_class_namesmultiple myeloma plasma cell  1.511608  0.6137485  2.462911 1.378142e-02 4.725057e-02
#trunc_cancer_class_namesneuroendocrine tumors         1.896433  0.7531176  2.518110 1.179864e-02 4.719456e-02
#trunc_cancer_class_namesovary                         1.441079  0.6652171  2.166328 3.028611e-02 8.076295e-02
#trunc_cancer_class_namespancreas                      1.896991  0.5878416  3.227044 1.250761e-03 1.000609e-02
#trunc_cancer_class_namesprostate                      1.340918  0.5615196  2.388017 1.693957e-02 5.081871e-02
#trunc_cancer_class_namesstomach                       2.195664  0.7581113  2.896229 3.776769e-03 2.266062e-02
#trunc_cancer_class_namestonsil                        1.530086  0.7295342  2.097347 3.596290e-02 8.631097e-02

#look at table of prevalence for esophagus vs healthy; very significant as well
prev_table <- table(clin_meta_sub$colibactin_binary, clin_meta_sub$trunc_cancer_class_names)
prop.test(t(prev_table[,c("healthy", "esophagus")]))
#p-value = 1.31e-06


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


#HUMANN genefamilies MCCM colibactin prevalence
#colibactin_df_sort
#n_neg n_pos  perc_neg   perc_pos            p            q signif
#base of tongue                   21     0 1.0000000 0.00000000 9.803790e-01 0.9803790220       
#liver intrahepatic bile ducts    38     1 0.9743590 0.02564103 9.712458e-01 0.9803790220       
#healthy                         279     8 0.9721254 0.02787456 1.000000e+00 1.0000000000       
#kidney except renal pelvis       29     1 0.9666667 0.03333333 7.699637e-01 0.8399604085       
#unspecified skin                 23     1 0.9583333 0.04166667 6.076212e-01 0.6944241720       
#breast                          224    11 0.9531915 0.04680851 2.124104e-01 0.2998735482       
#brain                            35     2 0.9459459 0.05405405 3.825481e-01 0.5100640776       
#bronchus lung                    85     5 0.9444444 0.05555556 1.686421e-01 0.2529631605       
#oropharynx                       16     1 0.9411765 0.05882353 4.138453e-01 0.5227519574       
#other connective soft tissue     42     3 0.9333333 0.06666667 1.628846e-01 0.2529631605       
#bladder                          44     4 0.9166667 0.08333333 5.418422e-02 0.0935150178      *
#prostate                        109    10 0.9159664 0.08403361 1.693957e-02 0.0508187054      *
#melanoma of skin                 43     4 0.9148936 0.08510638 5.078865e-02 0.0935150178      *
#corpus uteri                     29     3 0.9062500 0.09375000 5.455043e-02 0.0935150178      *
#ovary                            37     4 0.9024390 0.09756098 3.028611e-02 0.0807629484      *
#multiple myeloma plasma cell     54     6 0.9000000 0.10000000 1.378142e-02 0.0472505709     **
#colorectal                       59     7 0.8939394 0.10606061 7.049307e-03 0.0338366725     **
#tonsil                           25     3 0.8928571 0.10714286 3.596290e-02 0.0863109703      *
#lymphoid leukemia                15     2 0.8823529 0.11764706 5.409395e-02 0.0935150178      *
#neuroendocrine tumors            18     3 0.8571429 0.14285714 1.179864e-02 0.0471945583     **
#pancreas                         42     7 0.8571429 0.14285714 1.250761e-03 0.0100060920     **
#stomach                          13     3 0.8125000 0.18750000 3.776769e-03 0.0226606169     **
#esophagus                        26     8 0.7647059 0.23529412 2.267863e-05 0.0002721436    ***


#save as svg to edit
#png("./Figure/RM_figures/colibactin_cancer_classes_vs_control.png", width=5.75*300, height=5.75*300, res=300)
svg("./Figure/RM_figures/colibactin_cancer_classes_vs_control.svg", width=8, height=6)

par(mar = c(6, 14, 0.2, 1), las = 1, cex.axis = 1.05, cex.lab = 1.05) 
# Create a stacked barplot
bp <- barplot(colibactin_df_sort_plot,
              beside = FALSE,
              col = c("white", "#CC79A7"), 
              xlab = "Colibactin gene prevalence %",
              ylab = "",
              main = "", horiz = T,
              space=0, xlim=c(0,100))
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
  x=-60,
  y=-1,
  legend = c("colibactin negative", "colibactin positive", "*  q < 0.1", "** q < 0.05", "*** q < 0.01"),
  pch = c(15, 15, NA, NA, NA),
  col = c("white", "#CC79A7", "black", "black", "black"),
  bty = "n",
  pt.cex = 2,
  horiz = FALSE,
  xpd = TRUE                  # allow drawing outside plot region
)

dev.off()


#--------------------------------------
#do the same stats as above on strict colibactin; colibactin_binary_strict


model <- glm(colibactin_binary_strict ~ trunc_cancer_class_names + Age, data = clin_meta_sub, family = binomial)
model_res <- as.data.frame(summary(model)$coefficients)
#correct for multiple testing
model_res$FDR <- p.adjust(as.numeric(model_res[,"Pr(>|z|)"]), method = "fdr")
model_res[model_res$FDR < 0.1,]

#look at table of prevalence for esophagus vs healthy; very significant as well
prev_table <- table(clin_meta_sub$colibactin_binary_strict, clin_meta_sub$trunc_cancer_class_names)
prop.test(t(prev_table[,c("healthy", "esophagus")]))
#p-value = 1.31e-06


temp_table <- t(table(clin_meta_sub$colibactin_binary_strict, clin_meta_sub$trunc_cancer_class_names))
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
svg("./Figure/RM_figures/colibactin_cancer_classes_vs_control_strict.svg", width=8, height=6)

par(mar = c(6, 14, 0.2, 1), las = 1, cex.axis = 1.05, cex.lab = 1.05) 
# Create a stacked barplot
bp <- barplot(colibactin_df_sort_plot,
              beside = FALSE,
              col = c("white", "#CC79A7"), 
              xlab = "Colibactin gene prevalence %",
              ylab = "",
              main = "", horiz = T,
              space=0, xlim=c(0,100))
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
  x=-60,
  y=-1,
  legend = c("colibactin negative strict", "colibactin positive strict", "*  q < 0.1", "** q < 0.05", "*** q < 0.01"),
  pch = c(15, 15, NA, NA, NA),
  col = c("white", "#CC79A7", "black", "black", "black"),
  bty = "n",
  pt.cex = 2,
  horiz = FALSE,
  xpd = TRUE                  # allow drawing outside plot region
)

dev.off()


#--------------------------------------
#do the stats for cancer classes on the continuous variable colibactin counts;
#It is numerically increased but too sparse to be significant, except for rectosigmoid junction 

model <- glm(colibactin_counts ~ trunc_cancer_class_names + Age,
             data   = clin_meta,
             family = gaussian(link = "identity"))
#This specifies a continuous Gaussian outcome with identity link, i.e., a linear regression written in GLM form.
model_res <- as.data.frame(summary(model)$coefficients)
#correct for multiple testing
model_res$FDR <- p.adjust(as.numeric(model_res$`Pr(>|t|)`), method = "fdr")
model_res[model_res$FDR < 0.1,]

model_res[model_res$`Pr(>|t|)` < 0.05,]

head(model_res[order(model_res$Estimate, decreasing = T),], 20)


#--------------------------------------
#now make a boxplot for all cancer classes
  #make a cancer class other


n_groups <- nlevels(clin_meta_sub$trunc_cancer_class_names)

# gradient from blue to red (change colors as you like)
sel_colors_plot <- seq_gradient_pal("#0072B2", "#CC79A7", "Lab")(
  seq(0, 1, length.out = n_groups)
)


cancer_class_colibactin_counts <- clin_meta_sub %>%
  mutate(
    trunc_cancer_class_names = fct_reorder(
      trunc_cancer_class_names,
      colibactin_counts,
      .fun = mean,
      na.rm = TRUE
    )
  ) %>%
  ggplot(
    aes(x = colibactin_counts,
        y = trunc_cancer_class_names,
        fill = trunc_cancer_class_names)
  ) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(aes(color = trunc_cancer_class_names),
              height = 0.15, size = 2, alpha = 0.8) +
  scale_x_sqrt() +
  scale_fill_manual(values = sel_colors_plot) +
  scale_color_manual(values = sel_colors_plot) +
  theme_classic() +
  theme(
    legend.position = "none",
    text       = element_text(family = "Arial"),
    axis.text  = element_text(size = 9,  family = "Arial", colour = "black"),
    axis.title = element_text(size = 10, family = "Arial", colour = "black"),
    plot.title = element_text(size = 10, family = "Arial", colour = "black")
  ) +
  labs(
    x = "sqrt(colibactin gene counts, summed)",
    y = NULL
  )


cancer_class_colibactin_counts


ggsave("./Figure/RM_figures/counts_colibactin_cancer classes.png", cancer_class_colibactin_counts, 
       width = 5, height = 7, units = "in", dpi = 300)


#--------------------------------------
#colibactin results and methods


#Colibactin presence was assessed from metagenomic HUMANN3 genefamilies where 8 colibactin genes were annotated (Supplementary table X).
#As expected, these genes are highly correlated with eachother (median spearman rho of 0.89). 
#We therefore analyzed colibactin as presence / absence where samples where at least 3 / 8 genes were annotated were deemed colibactin positive. 
#(reducing the number of samples from 159 (9.6%) with any annotated colibactin gene to 128 with 3 or more (7.8%))
#Cancer samples had elevated presence of colibactin compared to healthy controls (8.8% and 2.8% prevalence in cancer and healthy controls, respectively)
#(logistic regression adjusted for Age, p = 0.0032) (boxplot Figure X).
#However, colibactin presence was distinct between cancer classes (prevalence stacked barplot Figure X) (logistic regression adjusted for Age).
#Colibactin was more prevalent in the fecal samples from colorectal cancer than healthy controls as documented (https://pubmed.ncbi.nlm.nih.gov/36056757/). 
#However, the largest prevalence was observed for esophageal, stomach, and pancreatic cancer (Figure X).


#--------------------------------------
#early onset and colibactin
#first inspect relationship with Age

clin_meta_cancer <- clin_meta[clin_meta$cancer_healthy == "cancer",]

plot(clin_meta_cancer$Age, clin_meta_cancer$colibactin_counts)
cor.test(clin_meta_cancer$colibactin_counts, clin_meta_cancer$Age) #not significant

clin_meta_cancer$early_onset <- "-" #only relevant for cancer
clin_meta_cancer$early_onset[clin_meta_cancer$Age > 50] <- "normal"
clin_meta_cancer$early_onset[clin_meta_cancer$Age <= 50] <- "early" 

table(clin_meta_cancer$early_onset)
#early normal 
#225   1139 


table(clin_meta_cancer$colibactin_binary, clin_meta_cancer$early_onset) #no relationship with age or early onset
#   early normal
#0   203   1015
#1    22    124

prop.test(table(clin_meta_cancer$colibactin_binary, clin_meta_cancer$early_onset))


clin_meta_CRC <- clin_meta_cancer[grep("colorectal", clin_meta_cancer$trunc_cancer_class_names),]

plot(clin_meta_CRC$Age, clin_meta_CRC$colibactin_counts)
cor.test(clin_meta_CRC$Age, clin_meta_CRC$colibactin_counts) #no relationship with early onset

#t.test(clin_meta_CRC$colibactin_counts ~ clin_meta_CRC$early_onset)



sel_colors_plot <- c("#0072B2", "#CC79A7")


set.seed(42)

CRC_colibactin_plot <- ggplot(clin_meta_CRC, aes(x = early_onset,
                                                   y = colibactin_counts,
                                                   fill = early_onset)) +
  geom_boxplot(outlier.shape = NA) +                
  geom_jitter(aes(color = early_onset),
              width = 0.15, size = 2, alpha = 0.8) + 
  scale_fill_manual(values = sel_colors_plot) +
  scale_color_manual(values = sel_colors_plot) +
  geom_segment(aes(y = 33, yend = 33, x = 1, xend = 2), size = 0.3) +
  geom_text(aes(y = 33*1.03, x = 1.5, label = "ns"), size = 3, vjust = 0) +
  theme_classic() +                                  # axis lines come back
  theme(
    panel.border = element_rect(color = "black", fill = NA, size = 0.5), 
    axis.ticks.y = element_line(color = "black", size = 0.3),
    axis.ticks.x = element_line(color = "black", size = 0.),
    axis.line.x  = element_line(color = "black", size = 0.1),
    axis.line.y  = element_line(color = "black", size = 0.1),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "none",
    plot.title  = element_text(size = 10, hjust = 0.5, family = "Arial", colour = "black"),
    axis.text   = element_text(size = 9, family = "Arial", colour = "black"),
    axis.title  = element_text(size = 10, family = "Arial", colour = "black"),
    plot.margin = grid::unit(c(1,1,1,1), "lines")
  ) +
  labs(
    y = "colibactin gene counts (summed)",
    x = NULL                                         # suppress x-axis title
  )


ggsave("./Figure/RM_figures/boxplot_CRC_colibactin.png", CRC_colibactin_plot, width = 2.25, height = 4.5, units = "in", dpi = 300)




#-------------------------------------------------------------------------------
#Comparing proportion positive from Yachida and MCCM healthies;


## Healthy contingency table HUMANN
##            MCCM  Yachida
## positive      8      16
## negative    279     235

humann_mat <- matrix(
  c(8, 279,
    16, 235),
  nrow = 2,
  byrow = TRUE,
  dimnames = list(
    c("MCCM", "Yachida"),
    c("positive", "negative")
  )
)

humann_mat

## Two-sample test for equality of proportions (positive vs total)
prop.test(
  x = humann_mat[ , "positive"],
  n = rowSums(humann_mat),
  correct = FALSE
)
#2-sample test for equality of proportions without continuity correction
#data:  humann_mat[, "positive"] out of rowSums(humann_mat)
#X-squared = 4.0424, df = 1, p-value = 0.04437
#alternative hypothesis: two.sided
#95 percent confidence interval:
#  -0.0715930184 -0.0001478925
#sample estimates:
#  prop 1     prop 2 
#0.02787456 0.06374502 




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


