
require(tidyverse)
require(RColorBrewer)
library(gt)
require(ComplexHeatmap)
require(circlize)
require(ggpubr)


setwd("/Users/M210320/Library/Mobile Documents/com~apple~CloudDocs/mforge_clean/")

load(file = 'Data/data.obj.raw.core.RData') 
clin_meta <- as.data.frame(data.obj$meta.dat)
tax_table <- as.data.frame(data.obj$otu.tab)

dim(tax_table) #7839 1651
#table(rowSums(tax_table > 0) == 0)

#3095 above 0 in more than 5% of the samples
table(rowSums(tax_table > 0) > ncol(tax_table) * 0.05)
#FALSE  TRUE 
#4744  3095


dim(data.obj$abund.list$Phylum) #41
dim(data.obj$abund.list$Family) #391
dim(data.obj$abund.list$Genus) #1899


#where do I find the pathway abundances?
load(file = 'Data/data.obj.pathway.RData') 
path_table <- as.data.frame(data.obj$otu.tab)

dim(path_table) #519 1651

table(rowSums(path_table > 0) > ncol(path_table) * 0.05)
#FALSE  TRUE 
# 90   429

#make function for heatmap from dataframe
source("~/Dropbox/Mayo_RS/R/general functions/zicoseq_heatmap_from_df_general.R")


#order the clin meta based on tax_table
clin_meta <- clin_meta[colnames(path_table),]
#1651  190; all samples in here


#-------------------------------------------------------------------------------
#make pathways

load("Figure/subCancerX_Control_func/pathway/DAA_P_R2.RData")
zico_res_ctr <- list(P.All$pathway, Q.All$pathway, R2.All$pathway)
names(zico_res_ctr) <- c("P", "Q", "R2")


#do.call(cbind, zico_res_ctr)

load("Figure/subCancerX-Ex_func/pathway/DAA_P_R2.RData")
zico_res_ex_cancer <- list(P.All$pathway, Q.All$pathway, R2.All$pathway)
names(zico_res_ex_cancer) <- c("P", "Q", "R2")



temp_dir <- "./Result/CancerOnly_func/pathway/DAA/"
elix_variables <- list.files(temp_dir)[grep("Elixhauser", list.files(temp_dir))]
#remove these
elix_variables <- elix_variables[!elix_variables %in% c("Elixhauser_Mets", "Elixhauser_Tumor", 
                                                        "Elixhauser_Lymphoma", "Elixhauser_Obesity",
                                                        "Elixhauser_Paralysis", "Elixhauser_HIV",
                                                        "Elixhauser_BloodLoss", "Elixhauser_Psychoses")]


#subset more categories
#already corrected for so remove?
#table(clin_meta$Elixhauser_Mets) #722 Metastatic Cancer
#table(clin_meta$Elixhauser_Tumor) #517 Solid tumor without metastasis
#table(clin_meta$Elixhauser_Obesity) #292 BMI >30 or higher, 
#table(clin_meta$Elixhauser_Lymphoma) #145

#remove these because too few patients
#table(clin_meta$Elixhauser_Paralysis) #11
#table(clin_meta$Elixhauser_HIV) #2
#table(clin_meta$Elixhauser_BloodLoss) #13 anemia due to blood loss
#table(clin_meta$Elixhauser_Psychoses) #only 2


length(elix_variables) #22


#see file "Cohort summary tables.R" for information on the Elixhauser components

elix_variables_filenames <- list.files(paste0(temp_dir, elix_variables))
elix_variables_filenames_full <- paste0(temp_dir, elix_variables, "/", elix_variables_filenames)
#select only the .Rdata files
elix_variables_filenames_full <- elix_variables_filenames_full[grep(".Rdata", elix_variables_filenames_full)]

elix_variables_diff_list <- list()
for (i in 1:length(elix_variables_filenames_full)) {
  load(elix_variables_filenames_full[i])
  elix_variables_diff_list[[i]] <- diff.obj
}

names(elix_variables_diff_list) <- elix_variables
length(elix_variables_diff_list) #22


#-------------------------------------------------------------------------------
#get the relevant information out, only the R2, coef, q value

#these variables are corrected for cancer class; does that make sense to do here?

length(elix_variables_diff_list$Elixhauser_Alcohol)
names(elix_variables_diff_list$Elixhauser_Alcohol)
#"pv.list"   "fc.list"   "pc.list"   "qv.list"   "R2.list"   "coef.list" "m.list"


names(elix_variables_diff_list$Elixhauser_Alcohol$R2.list)
head(elix_variables_diff_list$Elixhauser_Alcohol$R2.list$pathway)
#                                                                                        Func1
#1CMET2-PWY: folate transformations III (E. coli)                                     2.270607e-05
#ANAEROFRUCAT-PWY: homolactic fermentation                                            4.353506e-07
#ANAGLYCOLYSIS-PWY: glycolysis III (from glucose)                                     7.785427e-06


elix_variables_diff_list_Q <- lapply(elix_variables_diff_list, function(x) x$qv.list$pathway[,'Qvalue'])
elix_variables_diff_list_R2 <- lapply(elix_variables_diff_list, function(x) x$R2.list$pathway[,'Func1'])
elix_variables_diff_list_coef <- lapply(elix_variables_diff_list, function(x) x$coef.list$pathway)
coef_list_for_sign <- lapply(elix_variables_diff_list_coef, function(x) x[,ncol(x)])


elix_variables_diff_list_R2_sign <- elix_variables_diff_list_R2
for (i in 1:length(elix_variables_diff_list_R2)) {
  elix_temp <- elix_variables_diff_list_R2[[i]]
  elix_variables_diff_list_R2_sign[[i]] <- elix_temp * sign(coef_list_for_sign[[i]])
}

#lapply(elix_variables_diff_list_R2, head)
#lapply(elix_variables_diff_list_R2_sign, head)
#lapply(elix_variables_diff_list_Q, head)


#subset for pathways in all comparisons
temp_table <- table(unlist(lapply(elix_variables_diff_list_R2_sign, names))) == 22
names_to_include <- names(temp_table[temp_table == TRUE])

elix_variables_diff_list_R2_sign <- lapply(elix_variables_diff_list_R2_sign, function(x) x[names_to_include])
elix_variables_diff_list_Q <- lapply(elix_variables_diff_list_Q, function(x) x[names_to_include])

R2_sign_df_comorbidities <- do.call(cbind, elix_variables_diff_list_R2_sign)
Q_sign_df_comorbidities <- do.call(cbind, elix_variables_diff_list_Q)

dim(R2_sign_df_comorbidities); dim(Q_sign_df_comorbidities)
#334   22


#-------------------------------------------------------------------------------
#find strongest species linked to weightloss

sig_names <- names(which(Q_sign_df_comorbidities[,"Elixhauser_WeightLoss"] < 0.1))
names(head(sort(Q_sign_df_comorbidities[sig_names,"Elixhauser_WeightLoss"], decreasing = T), 10))
# [1] "PWY-3001: superpathway of L-isoleucine biosynthesis I"     "PWY-5686: UMP biosynthesis I"                             
# [3] "PWY-7851: coenzyme A biosynthesis II (eukaryotic)"         "POLYAMSYN-PWY: superpathway of polyamine biosynthesis I"  
# [5] "PWY-6151: S-adenosyl-L-methionine salvage I"               "P122-PWY: heterolactic fermentation"                      
# [7] "PWY0-1261: anhydromuropeptides recycling I"                "PWY0-1338: polymyxin resistance"                          
# [9] "ILEUSYN-PWY: L-isoleucine biosynthesis I (from threonine)" "PWY-7858: (5Z)-dodecenoate biosynthesis II"      


#-------------------------------------------------------------------------------
#also get the general Elixhauser associations

load("./Result/CancerOnly_func/pathway/DAA/Elix_score/Elix_score_ZicoSeq.Rdata")
Elix_res_list <- diff.obj

names(Elix_res_list$R2.list)
head(Elix_res_list$R2.list$pathway,3)
#                                                     Func1
#1CMET2-PWY: folate transformations III (E. coli) 4.965497e-05
#ANAEROFRUCAT-PWY: homolactic fermentation        4.850264e-04
#ANAGLYCOLYSIS-PWY: glycolysis III (from glucose) 3.223726e-05

elix_diff_Q <- Elix_res_list$qv.list$pathway[,'Qvalue']
elix_diff_R2 <- Elix_res_list$R2.list$pathway[,'Func1']
elix_diff_coef <- Elix_res_list$coef.list$pathway[,"Elix_score"]

elix_diff_R2 <- elix_diff_R2 * sign(elix_diff_coef)
length(elix_diff_R2[elix_diff_Q < 0.05])
#187
table(elix_diff_R2[elix_diff_Q < 0.05] > 0)
#FALSE (negatively associated with ECI)  TRUE (positively associated with ECI)
#             49                                        138 


#-------------------------------------
#get number of species FDR < 0.05 with positive and negative effect size for components

table_list <- list()
for (i in 1:ncol(R2_sign_df_comorbidities)) {
  R2_vec_temp <- R2_sign_df_comorbidities[,i]
  table_list[[i]] <- table(R2_vec_temp[Q_sign_df_comorbidities[,i] < 0.05] > 0)
}
names(table_list) <- colnames(R2_sign_df_comorbidities)

#$Elixhauser_Anemia
#FALSE  TRUE 
#   2   127 

#$Elixhauser_Arrhythmia
#FALSE 
#12 

#$Elixhauser_DM
#FALSE 
#4 

#$Elixhauser_DMcx
#FALSE  TRUE 
#   8   132 

#$Elixhauser_FluidsLytes
#FALSE  TRUE 
#   2    12 

#$Elixhauser_Rheumatic
#TRUE 
#2 

#$Elixhauser_Valvular
#TRUE 
#1 

#$Elixhauser_WeightLoss
#FALSE  TRUE 
# 21    51 

#sort(R2_sign_df_comorbidities[,"Elixhauser_Anemia"][Q_sign_df_comorbidities[,"Elixhauser_Anemia"] < 0.05])



#-------------------------------------------------------------------------------
#load all the control vs cancer class and cancer class vs other cancers comparisons

canX_Ctr <- read.csv("./Figure/subCancerX_Control_func/pathway_DAA.csv")
#head(canX_Ctr)
canX_Ctr$variable <- gsub("Control-", "", canX_Ctr$variable)

canX_Ctr_list <- lapply(as.list(unique(canX_Ctr$variable)), function(x) canX_Ctr[canX_Ctr$variable == x,]) 
names(canX_Ctr_list) <- unique(canX_Ctr$variable)
for (i in 1:length(canX_Ctr_list)) {
  rownames(canX_Ctr_list[[i]]) <- canX_Ctr_list[[i]]$pathway
}


canX_Ex <- read.csv("./Figure/subCancerX-Ex_func/pathway_DAA.csv")
canX_Ex_list <- lapply(as.list(unique(canX_Ex$variable)), function(x) canX_Ex[canX_Ex$variable == x,]) 
names(canX_Ex_list) <- unique(canX_Ex$variable)
#head(canX_Ex_list)
for (i in 1:length(canX_Ex_list)) {
  rownames(canX_Ex_list[[i]]) <- canX_Ex_list[[i]]$pathway
}

#made sure the names match and are in the same order
names(canX_Ctr_list) <- gsub("other_unspecified_malignant_neoplasm_of_skin", "unspecified_neoplasm_of_skin", names(canX_Ctr_list))
canX_Ctr_list <- canX_Ctr_list[order(names(canX_Ctr_list))]
#cbind(names(canX_Ctr_list), names(canX_Ex_list))



#-------------------------------------------------------------------------------
#only select those species significant <0.2 FDR in both control and cancer class comparisons

#trying this discovery level cutoff but mark significance in the heatmap
cutoff <- 0.2

merge_list <- vector("list", length = length(canX_Ctr_list))
names(merge_list) <- names(canX_Ctr_list)

#merge these in a for loop
for (i in 1:length(canX_Ctr_list)) {
  
  control <- canX_Ctr_list[[i]]
  cancer <- canX_Ex_list[[i]]
  
  idx.c <- control$pathway[which(control$Q<cutoff)]
  idx.d <- cancer$pathway[which(cancer$Q<cutoff)]
  idx.com <- intersect(idx.c, idx.d)
  
  control_sub <- control[idx.com,c("P", "R2", "Q")]
  cancer_sub <- cancer[idx.com,c("P", "R2", "Q")]
  
  #make sure directionality is the same
  control_sub_overlap <- control_sub[(control_sub$R2 < 0) == (cancer_sub$R2 < 0),]
  cancer_sub_overlap <- cancer_sub[(cancer_sub$R2 < 0) == (control_sub$R2 < 0),]
  
  temp_list <- list(control_sub_overlap, cancer_sub_overlap)
  names(temp_list) <- c("Ctr", "Ex")
  
  if (min(c(nrow(control_sub), nrow(cancer_sub))) > 1) {
    merge_list[[i]] <- do.call(cbind, temp_list)
  }
}
merge_list_sig <- merge_list[which(lapply(merge_list, function(x) nrow(x) > 0) == TRUE)]

#
lapply(merge_list_sig, dim)
#$esophagus
#[1] 9 2

#$liver_intrahepatic_bile_ducts
#[1] 190   2

#$neuroendocrine_tumors
#[1] 65  2

#$ovary
#[1] 37  2

#$pancreas
#[1] 2 2

merge_list_sig_Q <- lapply(merge_list_sig, function(x) x[,c("Ctr.Q", "Ex.Q")])
merge_list_sig_R2 <- lapply(merge_list_sig, function(x) x[,c("Ctr.R2", "Ex.R2")])

#head(merge_list_sig_Q[[1]])
#head(merge_list_sig_R2[[1]])


#TRUE indicates positive effect size, higher in cancer
lapply(merge_list_sig_R2, function(x) table(x$Ctr.R2 > 0))
#$esophagus
#FALSE  TRUE 
#1     8 
#$liver_intrahepatic_bile_ducts
#FALSE  TRUE 
#50   140 
#$neuroendocrine_tumors
#FALSE  TRUE 
#9    56 
#$ovary
#FALSE  TRUE 
#15    22 
#$pancreas
#FALSE 
#2 


lapply(merge_list_sig_Q, function(x) length(which(x$Ctr.Q < 0.1)))
#$esophagus
#[1] 6
#$liver_intrahepatic_bile_ducts
#[1] 180
#$neuroendocrine_tumors
#[1] 52
#$ovary
#[1] 27
#$pancreas
#[1] 1


lapply(merge_list_sig_Q, function(x) length(which(x$Ctr.Q < 0.05)))
#$esophagus
#[1] 4
#$liver_intrahepatic_bile_ducts
#[1] 167
#$neuroendocrine_tumors
#[1] 47
#$ovary
#1] 21
#$pancreas
#[1] 1


lapply(merge_list_sig_Q, function(x) length(which(x$Ctr.Q < 0.01)))
#$esophagus
#[1] 3
#$liver_intrahepatic_bile_ducts
#[1] 143
#$neuroendocrine_tumors
#[1] 23
#ovary
#[1] 6
#$pancreas
#[1] 0


# Sample data
df <- data.frame(
  Name = c("liver intrahepatic bile ducts", "neuroendocrine tumors", "ovary", 
           "esophagus", "pancreas"),
  "n pathways FDR <0.2" = c(190, 65, 37, 9, 2),
  "n pathways FDR <0.1" = c(180, 52, 27, 6, 1),
  "n pathways FDR <0.05" = c(167, 47, 21, 4, 1),
  "n pathways FDR <0.01" = c(143, 23, 6, 3, 0))

df <- df[order(rowSums(df[,2:ncol(df)]), decreasing = T),]

   
gt(df) %>%
  tab_header(title = "Number of functional pathways overlapping vs. healthy and vs. other cancers") %>%
  cols_label(Name="Name", n.pathways.FDR..0.2= "n pathways q<0.2", 
             n.pathways.FDR..0.1="n pathways q<0.1", n.pathways.FDR..0.05="n pathways q<0.05", 
             n.pathways.FDR..0.01="n pathways q<0.01") %>%
  opt_vertical_padding(scale = 0.1)



#-------------------------------------------------------------------------------
#add the comorbidities from above

merge_list_sig_Q_w_comorb <- merge_list_sig_Q
merge_list_sig_R2_w_comorb <- merge_list_sig_R2

for (i in 1:length(merge_list_sig_Q)) {
  #Q
  merge_list_temp <- merge_list_sig_Q[[i]]
  merge_list_temp_sub <- merge_list_temp[rownames(merge_list_temp) %in% rownames(Q_sign_df_comorbidities),]
  merge_list_temp_sub_comorb <- cbind(merge_list_temp_sub, Q_sign_df_comorbidities[rownames(merge_list_temp_sub),])
  merge_list_sig_Q_w_comorb[[i]] <- merge_list_temp_sub_comorb
  
  #R2
  merge_list_temp <- merge_list_sig_R2[[i]]
  merge_list_temp_sub <- merge_list_temp[rownames(merge_list_temp) %in% rownames(R2_sign_df_comorbidities),]
  merge_list_temp_sub_comorb <- cbind(merge_list_temp_sub, R2_sign_df_comorbidities[rownames(merge_list_temp_sub),])
  merge_list_sig_R2_w_comorb[[i]] <- merge_list_temp_sub_comorb

}


#$esophagus
#[1] 9 2

#$liver_intrahepatic_bile_ducts
#[1] 190   2

#$neuroendocrine_tumors
#[1] 65  2

#$ovary
#[1] 37  2

#$pancreas
#[1] 2 2


#-------------------------------------------------------------------------------
i="esophagus"

#"esophagus"  "liver_intrahepatic_bile_ducts"  "neuroendocrine_tumors"  "ovary"  "pancreas" 
FDR_df <- merge_list_sig_Q_w_comorb[[i]]
R2_df <- merge_list_sig_R2_w_comorb[[i]]
grp.names <- colnames(FDR_df)
grp.labels <- colnames(FDR_df)
grp.labels[1:2] <- c("Cancer class vs. Control", "Cancer class vs. other cancers")
sign_cutoffs <- c(0.2, 0.1, 0.05)
effect_size_cutoff=1

heatmap_list <- zicoseq_heatmap(FDR_df, R2_df, sign_cutoffs=sign_cutoffs, grp.labels, single_cutoff = FALSE, 
                                sign_filter = FALSE, effect_size_cutoff=1, rotate=TRUE)
heatmap_list$heatmap_plot

filename <- "./Figure/RM_figures/pathway overlap heatmap esophagus all.png"
png(filename, width =6.5*330, height = 9*330, res = 330)
heatmap_list$heatmap_plot
dev.off()



#--------------------------------------
i="liver_intrahepatic_bile_ducts"
#"esophagus"  "liver_intrahepatic_bile_ducts"  "neuroendocrine_tumors"  "ovary"  "pancreas"

FDR_df <- merge_list_sig_Q_w_comorb[[i]]
R2_df <- merge_list_sig_R2_w_comorb[[i]]

heatmap_list <- zicoseq_heatmap(FDR_df, R2_df, sign_cutoffs=sign_cutoffs, grp.labels, single_cutoff = FALSE, 
                                sign_filter = FALSE, effect_size_cutoff=1, rotate=TRUE)
heatmap_list$heatmap_plot

filename <- "./Figure/RM_figures/pathway overlap heatmap liver_intrahepatic_bile_ducts all.png"
png(filename, width = 30*330, height = 13*330, res = 330)
heatmap_list$heatmap_plot
dev.off()


#only top 25%
heatmap_list <- zicoseq_heatmap(FDR_df, R2_df, sign_cutoffs=sign_cutoffs, grp.labels, single_cutoff = FALSE, 
                                sign_filter = FALSE, effect_size_cutoff=0.25, rotate=TRUE)

filename <- "./Figure/RM_figures/pathway overlap heatmap liver_intrahepatic_bile_ducts 25perc.png"
png(filename, width = 16*330, height = 13*330, res = 330)
heatmap_list$heatmap_plot
dev.off()



#--------------------------------------
i="neuroendocrine_tumors"
#"esophagus"  "liver_intrahepatic_bile_ducts"  "neuroendocrine_tumors"  "ovary"  "pancreas" 

FDR_df <- merge_list_sig_Q_w_comorb[[i]]
R2_df <- merge_list_sig_R2_w_comorb[[i]]

heatmap_list <- zicoseq_heatmap(FDR_df, R2_df, sign_cutoffs=sign_cutoffs, grp.labels, single_cutoff = FALSE, 
                                sign_filter = FALSE, effect_size_cutoff=1, rotate=TRUE)
heatmap_list$heatmap_plot

filename <- "./Figure/RM_figures/pathway overlap heatmap neuroendocrine_tumors all.png"
png(filename, width = 16*330, height = 13*330, res = 330)
heatmap_list$heatmap_plot
dev.off()



#--------------------------------------
i="ovary"
#"esophagus"  "liver_intrahepatic_bile_ducts"  "neuroendocrine_tumors"  "ovary"  "pancreas"
FDR_df <- merge_list_sig_Q_w_comorb[[i]]
R2_df <- merge_list_sig_R2_w_comorb[[i]]
grp.names <- colnames(FDR_df)
grp.labels <- colnames(FDR_df)
grp.labels[1:2] <- c("Cancer class vs. Control", "Cancer class vs. other cancers")
sign_cutoffs <- c(0.2, 0.1, 0.05)


heatmap_list <- zicoseq_heatmap(FDR_df, R2_df, sign_cutoffs=sign_cutoffs, grp.labels, single_cutoff = FALSE, 
                                sign_filter = FALSE, effect_size_cutoff=1, rotate=TRUE)
heatmap_list$heatmap_plot

filename <- "./Figure/RM_figures/pathway overlap heatmap ovary all.png"
png(filename, width =14*330, height = 13*330, res = 330)
heatmap_list$heatmap_plot
dev.off()



#--------------------------------------
i="pancreas"
#"esophagus"  "liver_intrahepatic_bile_ducts"  "neuroendocrine_tumors"  "ovary"  "pancreas"

FDR_df <- merge_list_sig_Q_w_comorb[[i]]
R2_df <- merge_list_sig_R2_w_comorb[[i]]

heatmap_list <- zicoseq_heatmap(FDR_df, R2_df, sign_cutoffs=sign_cutoffs, grp.labels, single_cutoff = FALSE, 
                                sign_filter = FALSE, effect_size_cutoff=1, rotate=TRUE)
heatmap_list$heatmap_plot

filename <- "./Figure/RM_figures/pathway overlap heatmap pancreas all.png"
png(filename, width = 4.5*330, height = 9*330, res = 330)
heatmap_list$heatmap_plot
dev.off()



#-------------------------------------------------------------------------------
#cancer classes of interest

#"esophagus"  "liver_intrahepatic_bile_ducts"  "neuroendocrine_tumors"  "ovary"  "pancreas"
#clin_meta$icd10_first_3_name[grep("ovary", clin_meta$icd10_first_3_name)]

cancer_class_w_sign <- c("Malignant neoplasm of esophagus", "Malignant neoplasm of liver and intrahepatic bile ducts",
                         "Malignant neuroendocrine tumors", "Malignant neoplasm of ovary", "Malignant neoplasm of pancreas")

cancer_class_w_sign_row_list <- sapply(cancer_class_w_sign, function(x) which(clin_meta$icd10_first_3_name == x))

#tabulate all Elixhauser components with a significant species
test_df <- lapply(merge_list_sig_Q_w_comorb, function(x) (x[,3:ncol(x)] < 0.2)) #0.2 is FDR cutoff

#general one
Elix_components_w_sig <- unique(unlist(lapply(test_df, function(x) colnames(x)[colSums(x == TRUE) > 0])))
#[1] "Elixhauser_Anemia"      "Elixhauser_Arrhythmia"  "Elixhauser_DM"          "Elixhauser_DMcx"        "Elixhauser_FluidsLytes" "Elixhauser_WeightLoss"

#make Elixhauser components specific to that cancer class
Elix_components_w_sig_list <- lapply(test_df, function(x) colnames(x)[colSums(x == TRUE) > 0])


component_list <- list()
cancer_class_w_sign_row_list_sub <- cancer_class_w_sign_row_list[lapply(Elix_components_w_sig_list, length) > 0] 
Elix_components_w_sig_list_sub <- Elix_components_w_sig_list[lapply(Elix_components_w_sig_list, length) > 0]


#the 2 pancreas species do not have associations with comorbidities
names(cancer_class_w_sign_row_list_sub)
#[1] "Malignant neoplasm of esophagus"                         "Malignant neoplasm of liver and intrahepatic bile ducts"
#[3] "Malignant neuroendocrine tumors"                         "Malignant neoplasm of ovary"                            
#[5] "Malignant neoplasm of pancreas"    


#-------------------------------------------------------------------------------
#for cancer groups that have only very few patients with these comorbidities these comorbidities do not have to be plotted or considered
#only consider more than >15%


#-------------------------------------------------------------------------------
#"Malignant neoplasm of esophagus"  
i=1
length(cancer_class_w_sign_row_list_sub[[i]]) #34
temp_table <- apply(clin_meta[cancer_class_w_sign_row_list_sub[[i]], Elix_components_w_sig_list_sub[[i]]], 2, table)
#           Elixhauser_Anemia Elixhauser_Arrhythmia Elixhauser_DM Elixhauser_DMcx Elixhauser_FluidsLytes Elixhauser_WeightLoss
#No                 30                    25            32              32                     27                    32
#Yes                 4                     9             2               2                      7                     2

#get the percentages
temp_table / 34 * 100
#         Elixhauser_Anemia Elixhauser_Arrhythmia Elixhauser_DM Elixhauser_DMcx Elixhauser_FluidsLytes Elixhauser_WeightLoss
#No           88.23529              73.52941     94.117647       94.117647               79.41176             94.117647
#Yes          11.76471              26.47059      5.882353        5.882353               20.58824              5.882353


#none >15%
selected_elix_esop <- c("Elixhauser_Arrhythmia", "Elixhauser_FluidsLytes")


#----------------------------------
#"Malignant neoplasm of liver and intrahepatic bile ducts"
i=2
length(cancer_class_w_sign_row_list_sub[[i]]) #39
temp_table <- apply(clin_meta[cancer_class_w_sign_row_list_sub[[i]], Elix_components_w_sig_list_sub[[i]]], 2, table)
#           Elixhauser_Anemia Elixhauser_Arrhythmia Elixhauser_DM Elixhauser_DMcx Elixhauser_FluidsLytes Elixhauser_WeightLoss
#No                 37                    30            34              32                     30                    32
#Yes                 2                     9             5               7                      9                     7

temp_table / 39 * 100
#         Elixhauser_Anemia Elixhauser_Arrhythmia Elixhauser_DM Elixhauser_DMcx Elixhauser_FluidsLytes Elixhauser_WeightLoss
#No          94.871795              76.92308      87.17949        82.05128               76.92308              82.05128
#Yes          5.128205              23.07692      12.82051        17.94872               23.07692              17.94872

selected_elix_liv <- c("Elixhauser_Arrhythmia", "Elixhauser_DMcx", "Elixhauser_FluidsLytes", "Elixhauser_WeightLoss")


#----------------------------------
#"Malignant neuroendocrine tumors"  
i=3
length(cancer_class_w_sign_row_list_sub[[i]]) #21
temp_table <- do.call(cbind, apply(clin_meta[cancer_class_w_sign_row_list_sub[[i]], Elix_components_w_sig_list_sub[[i]]], 2, table))
#         Elixhauser_Anemia Elixhauser_Arrhythmia Elixhauser_DM Elixhauser_DMcx Elixhauser_FluidsLytes Elixhauser_WeightLoss
#No                 19                    14            21              20                     18                    18
#Yes                 2                     7            0               1                      3                     3

temp_table / 21 * 100
#       Elixhauser_Anemia Elixhauser_Arrhythmia Elixhauser_DM Elixhauser_DMcx Elixhauser_FluidsLytes Elixhauser_WeightLoss
#No           90.47619              66.66667           100       95.238095               85.71429              85.71429
#Yes           9.52381              33.33333           0          4.761905               14.28571              14.28571

selected_elix_neuro <- c("Elixhauser_Arrhythmia")


#-------------------------------------------------------------------------------
#"Malignant neoplasm of ovary"
i=4
length(cancer_class_w_sign_row_list_sub[[i]]) #41
temp_table <- apply(clin_meta[cancer_class_w_sign_row_list_sub[[i]], Elix_components_w_sig_list_sub[[i]]], 2, table)
#         Elixhauser_Anemia Elixhauser_Arrhythmia Elixhauser_DM Elixhauser_DMcx Elixhauser_FluidsLytes Elixhauser_WeightLoss
#No                 39                    33            40              39                     29                    36
#Yes                 2                     8             1               2                     12                     5

#get the percentages
temp_table / 41 * 100
#       Elixhauser_Anemia Elixhauser_Arrhythmia Elixhauser_DM Elixhauser_DMcx Elixhauser_FluidsLytes Elixhauser_WeightLoss
#No          95.121951               80.4878     97.560976       95.121951               70.73171              87.80488
#Yes          4.878049               19.5122      2.439024        4.878049               29.26829              12.19512

#>15%
selected_elix_ovary <- c("Elixhauser_Arrhythmia", "Elixhauser_FluidsLytes")


#----------------------------------
#"Malignant neoplasm of pancreas" 
i=5
length(cancer_class_w_sign_row_list_sub[[i]]) #49
temp_table <- apply(clin_meta[cancer_class_w_sign_row_list_sub[[i]], Elix_components_w_sig_list_sub[[i]]], 2, table)
#       Elixhauser_DMcx Elixhauser_FluidsLytes
#No               38                     35
#Yes              11                     14

temp_table / 49 * 100
#       Elixhauser_DMcx Elixhauser_FluidsLytes
#No         77.55102               71.42857
#Yes        22.44898               28.57143

selected_elix_pancreas <- c("Elixhauser_DMcx", "Elixhauser_FluidsLytes")



#-------------------------------------------------------------------------------
#plot using the subset of the Elixhauser categories

#1 selected_elix_esop
#2 selected_elix_liv
#3 selected_elix_neuro
#4 selected_elix_ovary
#5 selected_elix_pancreas

Elix_components_w_sig_list_sub_cutoff <- Elix_components_w_sig_list_sub

#Bronchus esophagus
Elix_components_w_sig_list_sub_cutoff[[1]] <- Elix_components_w_sig_list_sub[[1]][Elix_components_w_sig_list_sub[[1]] %in% 
                                                                                    selected_elix_esop]
#liver_intrahepatic_bile_ducts
Elix_components_w_sig_list_sub_cutoff[[2]] <- Elix_components_w_sig_list_sub[[2]][Elix_components_w_sig_list_sub[[2]] %in% 
                                                                                    selected_elix_liv]
#neuroendocrine_tumors
Elix_components_w_sig_list_sub_cutoff[[3]] <- Elix_components_w_sig_list_sub[[3]][Elix_components_w_sig_list_sub[[3]] %in% 
                                                                                    selected_elix_neuro]
#ovary
Elix_components_w_sig_list_sub_cutoff[[4]] <- Elix_components_w_sig_list_sub[[4]][Elix_components_w_sig_list_sub[[4]] %in% 
                                                                                    selected_elix_ovary]
#pancreas
Elix_components_w_sig_list_sub_cutoff[[5]] <- Elix_components_w_sig_list_sub[[5]][Elix_components_w_sig_list_sub[[5]] %in% 
                                                                                    selected_elix_pancreas]


#--------------------------------------

sign_lists_to_plot <- list()


#[1] "esophagus"   "liver_intrahepatic_bile_ducts" "neuroendocrine_tumors"  "ovary"   "pancreas" 

#--------------------------------------
i=1
#esophagus
names(merge_list_sig_Q_w_comorb)[i]
sig_components <- Elix_components_w_sig_list_sub_cutoff[[i]]
#"Elixhauser_Arrhythmia"  "Elixhauser_FluidsLytes"
FDR_df <- merge_list_sig_Q_w_comorb[[i]][,c("Ctr.Q","Ex.Q",sig_components)]
R2_df <- merge_list_sig_R2_w_comorb[[i]][,c("Ctr.R2","Ex.R2",sig_components)]

grp.names <- colnames(FDR_df)
grp.labels <- colnames(FDR_df)
grp.labels[1:2] <- c("Cancer class vs. Control", "Cancer class vs. other cancers")
sign_cutoffs <- c(0.2, 0.1, 0.05)

heatmap_list <- zicoseq_heatmap(FDR_df, R2_df, sign_cutoffs=sign_cutoffs, grp.labels, single_cutoff = FALSE, 
                                sign_filter = FALSE, effect_size_cutoff=1, rotate=TRUE)
heatmap_list$heatmap_plot


filename <- "./Figure/RM_figures/pathway overlap heatmap esophagus_sig_Elix.png"
png(filename, width = 6*330, height = 6.5*330, res = 330)
heatmap_list$heatmap_plot
dev.off()

colnames(FDR_df) <- grp.labels
sign_lists_to_plot$esophagus <- FDR_df


#--------------------------------------
i=2
#liver_intrahepatic_bile_ducts
names(merge_list_sig_Q_w_comorb)[i]

sig_components <- Elix_components_w_sig_list_sub_cutoff[[i]]
FDR_df <- merge_list_sig_Q_w_comorb[[i]][,c("Ctr.Q","Ex.Q",sig_components)]
R2_df <- merge_list_sig_R2_w_comorb[[i]][,c("Ctr.R2","Ex.R2",sig_components)]

grp.names <- colnames(FDR_df)
grp.labels <- colnames(FDR_df)
grp.labels[1:2] <- c("Cancer class vs. Control", "Cancer class vs. other cancers")
sign_cutoffs <- c(0.2, 0.1, 0.05)

heatmap_list <- zicoseq_heatmap(FDR_df, R2_df, sign_cutoffs=sign_cutoffs, grp.labels, single_cutoff = FALSE, 
                                sign_filter = FALSE, effect_size_cutoff=1, rotate=TRUE)
heatmap_list$heatmap_plot


filename <- "./Figure/RM_figures/pathway overlap heatmap liver_intrahepatic_bile_ducts_sig_Elix.png"
png(filename, width = 28*330, height = 9.5*330, res = 330)
heatmap_list$heatmap_plot
dev.off()

colnames(FDR_df) <- grp.labels
sign_lists_to_plot$liver_intrahepatic_bile_ducts <- FDR_df



heatmap_list <- zicoseq_heatmap(FDR_df, R2_df, sign_cutoffs=sign_cutoffs, grp.labels, single_cutoff = FALSE, 
                                sign_filter = FALSE, effect_size_cutoff=0.25, rotate=TRUE)

filename <- "./Figure/RM_figures/pathway overlap heatmap liver_intrahepatic_bile_ducts_sig_Elix top 25perc.png"
png(filename, width = 15*330, height = 9*330, res = 330)
heatmap_list$heatmap_plot
dev.off()


#remove Elix species; any rows with FDR <0.2 beyond the first 2 columns
to_keep_rows <- !rowSums(FDR_df[, 3:ncol(FDR_df)] < 0.2) > 0
length(which(to_keep_rows)) 

heatmap_list <- zicoseq_heatmap(FDR_df[to_keep_rows,], R2_df[to_keep_rows,], sign_cutoffs=sign_cutoffs, grp.labels, single_cutoff = FALSE, 
                                sign_filter = FALSE, effect_size_cutoff=0.25, rotate=TRUE)

filename <- "./Figure/RM_figures/pathway overlap heatmap liver_intrahepatic_bile_ducts_sig_Elix filter.png"
png(filename, width = 4.5*330, height = 6.5*330, res = 330)
heatmap_list$heatmap_plot
dev.off()



#--------------------------------------
i=3
#neuroendocrine_tumors
names(merge_list_sig_Q_w_comorb)[i]

sig_components <- Elix_components_w_sig_list_sub_cutoff[[i]]
FDR_df <- merge_list_sig_Q_w_comorb[[i]][,c("Ctr.Q","Ex.Q",sig_components)]
R2_df <- merge_list_sig_R2_w_comorb[[i]][,c("Ctr.R2","Ex.R2",sig_components)]

grp.names <- colnames(FDR_df)
grp.labels <- colnames(FDR_df)
grp.labels[1:2] <- c("Cancer class vs. Control", "Cancer class vs. other cancers")
sign_cutoffs <- c(0.2, 0.1, 0.05)

heatmap_list <- zicoseq_heatmap(FDR_df, R2_df, sign_cutoffs=sign_cutoffs, grp.labels, single_cutoff = FALSE, 
                                sign_filter = FALSE, effect_size_cutoff=1, rotate=TRUE)
heatmap_list$heatmap_plot


filename <- "./Figure/RM_figures/pathway overlap heatmap neuroendocrine_tumors_sig_Elix.png"
png(filename, width = 20*330, height = 8*330, res = 330)
heatmap_list$heatmap_plot
dev.off()


colnames(FDR_df) <- grp.labels
sign_lists_to_plot$neuroendocrine_tumors <- FDR_df



#--------------------------------------
i=4
names(merge_list_sig_Q_w_comorb)[i]
#ovary
sig_components <- Elix_components_w_sig_list_sub_cutoff[[i]]

FDR_df <- merge_list_sig_Q_w_comorb[[i]][,c("Ctr.Q","Ex.Q",sig_components)]
R2_df <- merge_list_sig_R2_w_comorb[[i]][,c("Ctr.R2","Ex.R2",sig_components)]

grp.names <- colnames(FDR_df)
grp.labels <- colnames(FDR_df)
grp.labels[1:2] <- c("Cancer class vs. Control", "Cancer class vs. other cancers")
sign_cutoffs <- c(0.2, 0.1, 0.05)

heatmap_list <- zicoseq_heatmap(FDR_df, R2_df, sign_cutoffs=sign_cutoffs, grp.labels, single_cutoff = FALSE, 
                                sign_filter = FALSE, effect_size_cutoff=1, rotate=TRUE)
heatmap_list$heatmap_plot


filename <- "./Figure/RM_figures/pathway overlap heatmap ovary_sig_Elix.png"
png(filename, width =14*330, height = 8.75*330, res = 330)
heatmap_list$heatmap_plot
dev.off()

colnames(FDR_df) <- grp.labels
sign_lists_to_plot$ovary <- FDR_df



#--------------------------------------
i=5
#pancreas
names(merge_list_sig_Q_w_comorb)[i]

sig_components <- Elix_components_w_sig_list_sub_cutoff[[i]]
FDR_df <- merge_list_sig_Q_w_comorb[[i]][,c("Ctr.Q","Ex.Q",sig_components)]
R2_df <- merge_list_sig_R2_w_comorb[[i]][,c("Ctr.R2","Ex.R2",sig_components)]

grp.names <- colnames(FDR_df)
grp.labels <- colnames(FDR_df)
grp.labels[1:2] <- c("Cancer class vs. Control", "Cancer class vs. other cancers")
sign_cutoffs <- c(0.2, 0.1, 0.05)

heatmap_list <- zicoseq_heatmap(FDR_df, R2_df, sign_cutoffs=sign_cutoffs, grp.labels, single_cutoff = FALSE, 
                                sign_filter = FALSE, effect_size_cutoff=1, rotate=TRUE)
heatmap_list$heatmap_plot


filename <- "./Figure/RM_figures/pathway overlap heatmap pancreas_sig_Elix.png"
png(filename, width = 4.5*330, height = 4.5*330, res = 330)
heatmap_list$heatmap_plot
dev.off()


colnames(FDR_df) <- grp.labels
sign_lists_to_plot$pancreas <- FDR_df



#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#plot species abundances stratified by comorbidity

#boxplot; healthy control, all cancer, all cancer + comorbid, cancer class, cancer class + comorbid

clin_meta$icd10_first_3_name[is.na(clin_meta$icd10_first_3_name)] <- "Healthy"
head(sort(table(clin_meta$icd10_first_3_name), decreasing = T), 10) #testing

#update cancer_row_list in same order as sign_lists_to_plot
names(sign_lists_to_plot)
#"esophagus"  "liver_intrahepatic_bile_ducts" "neuroendocrine_tumors"  "ovary"  "pancreas"       

group_vec <- clin_meta$icd10_first_3_name
cancer_row_list <- list(
  grep("esophagus", group_vec),
  grep("liver and intrahepatic bile ducts", group_vec),
  grep("neuroendocrine", group_vec),
  grep("pancreas", group_vec),
  grep("ovary", group_vec))


elix_row_names_list <- lapply(sign_lists_to_plot, function(x) colnames(x)[3:ncol(x)])
healthy_rows <- which(group_vec == "Healthy")

sel_cols <- brewer.pal(5, "Dark2")

#all can be matched
#as.character(unlist(lapply(sign_lists_to_plot, rownames))) %in% rownames(path_table)

for (i in 1:length(sign_lists_to_plot)) { #cancer class
  sig_df <- sign_lists_to_plot[[i]]
  elix_row_names <- elix_row_names_list[[i]]
  cancer_group_rows <- cancer_row_list[[i]]
  all_other_cancer_rows <- which(!1:ncol(path_table) %in% c(healthy_rows, cancer_group_rows))
  cancer_class <- names(sign_lists_to_plot)[i]
  
  pdf(paste0("./Figure/RM_figures/Comorbidity pathways_", cancer_class, ".pdf"), width=5, height=4.5)
  for (j in 1:length(elix_row_names)) { #elix component
    sig_pathways <- rownames(sig_df)[which(sig_df[elix_row_names[j]] < 0.2)]
    
    for (k in 1:length(sig_pathways)) {
      
      temp_ra <- as.numeric(path_table[sig_pathways[k],])
      temp_elix_rows <- which(clin_meta[,elix_row_names[j]] == "Yes")
      
      #tax_table, setdiff indicates not having the comorbidity, intersect indicates having the comorbidity
      par(mar=c(7,4,2,1))
      boxplot(temp_ra[healthy_rows], temp_ra[setdiff(all_other_cancer_rows, temp_elix_rows)],
              temp_ra[intersect(all_other_cancer_rows, temp_elix_rows)], 
              temp_ra[setdiff(cancer_group_rows, temp_elix_rows)], 
              temp_ra[intersect(cancer_group_rows, temp_elix_rows)], las=1, col = sel_cols, varwidth = F, 
              xaxt = "n", outline = F)
      axis(1,  at=1:5, labels = c("Healthy", "other cancer\nno comorbidity", "other cancer\ncomorbidity", 
                         "cancer class\nno comorbidity", "cancer class\ncomorbidity"), las=2)
      mtext("pathway abundance", 2, line = 2.5)
      mtext(paste0(cancer_class, " _ ", elix_row_names[j], "\n", sig_pathways[k]), 3, line = 0, cex=1)
      
    }
  }
  dev.off()
}


#-------------------------------------
#split version 

#orange       skyblue   bluishgreen        yellow          blue    vermillion reddishpurple 
#"#E69F00"     "#56B4E9"     "#009E73"     "#F0E442"     "#0072B2"     "#D55E00"     "#CC79A7" 

sel_cols <- c("#0072B2", "#CC79A7", "#E69F00")

for (i in 1:length(sign_lists_to_plot)) { #cancer class
  sig_df <- sign_lists_to_plot[[i]]
  elix_row_names <- elix_row_names_list[[i]]
  cancer_group_rows <- cancer_row_list[[i]]
  all_other_cancer_rows <- which(!1:ncol(path_table) %in% c(healthy_rows, cancer_group_rows))
  cancer_class <- names(sign_lists_to_plot)[i]
  
  pdf(paste0("./Figure/RM_figures/Comorbidity pathways_split_", cancer_class, ".pdf"), width=5.5, height=4.5)
  for (j in 1:length(elix_row_names)) { #elix component
    sig_pathways <- rownames(sig_df)[which(sig_df[elix_row_names[j]] < 0.2)]
      
    for (k in 1:length(sig_pathways)) {
      
      temp_ra <- as.numeric(path_table[sig_pathways[k],])
      temp_elix_rows <- which(clin_meta[,elix_row_names[j]] == "Yes")
      
      #make a dataframe with these values
      #add filler group
      data <- as.data.frame(cbind(values = c(temp_ra[healthy_rows], 
        temp_ra[all_other_cancer_rows],
        temp_ra[cancer_group_rows],
        rep(NA, times=10)),
        group = c(rep("healthy", times=length(healthy_rows)),
          rep("all other cancer", times=length(all_other_cancer_rows)),
          rep("cancer class", times=length(cancer_group_rows)),
          rep("filler", times=10))))
      data$values <- as.numeric(data$values) 
      data$group <- as.factor(data$group)
      
      data <- data %>%
        mutate(group = fct_relevel(group, "healthy", "all other cancer", "cancer class", "filler"))
    
      #----------------
      data_all_other <- as.data.frame(cbind(values = c(rep(NA, times=10),
      rep(NA, times=10), temp_ra[setdiff(all_other_cancer_rows, temp_elix_rows)],
      temp_ra[intersect(all_other_cancer_rows, temp_elix_rows)], rep(NA, times=10), rep(NA, times=10)),
        group = c(rep("filler1", times=10),
                  rep("filler2", times=10),
                  rep("other cancer\nno comorbidity", times=length(setdiff(all_other_cancer_rows, temp_elix_rows))),
                  rep("other cancer\ncomorbidity", times=length(intersect(all_other_cancer_rows, temp_elix_rows))),
                  rep("filler3", times=10),
                  rep("filler4", times=10))))
      
      data_all_other$values <- as.numeric(data_all_other$values) 
      data_all_other$group <- as.factor(data_all_other$group)
      
      data_all_other <- data_all_other %>%
        mutate(group = fct_relevel(group, "filler1", "filler2", "other cancer\nno comorbidity", 
                                   "other cancer\ncomorbidity", "filler3", "filler4"))
      
      #----------------
      data_cc_other <- as.data.frame(cbind(values = c(rep(NA, times=10), rep(NA, times=10), 
                                                      rep(NA, times=10), rep(NA, times=10),
                                                      temp_ra[setdiff(cancer_group_rows, temp_elix_rows)],
                                                      temp_ra[intersect(cancer_group_rows, temp_elix_rows)]),
                                            group = c(rep("filler1", times=10),
                                                      rep("filler2", times=10),
                                                      rep("filler3", times=10),
                                                      rep("filler4", times=10),
                                                      rep("cancer class\nno comorbidity", times=length(setdiff(cancer_group_rows, temp_elix_rows))),
                                                      rep("cancer class\ncomorbidity", times=length(intersect(cancer_group_rows, temp_elix_rows)))
                                                      )))
      data_cc_other$values <- as.numeric(data_cc_other$values) 
      data_cc_other$group <- as.factor(data_cc_other$group)
      
      data_cc_other <- data_cc_other %>%
        mutate(group = fct_relevel(group, "filler1", "filler2", "filler3", "filler4",
                                   "cancer class\nno comorbidity", 
                                   "cancer class\ncomorbidity", ))
      
      ymax <- quantile(rbind(data, data_all_other, data_cc_other)[,1], na.rm = T, probs = c(0.95))
      
      par(mar=c(7,4,2,1))
      boxplot(values ~ group, data, las=1, col = sel_cols, varwidth = F, xaxt = "n", outline = F, xlab="", 
              at=c(1,2,4,5), ylim=c(0, ymax), ylab="", xlim=c(0.5, 5.5))
      #add insert boxplots
      par(new = TRUE)
      boxplot(values ~ group, data_all_other, las=1, col = sel_cols[2], varwidth = F, xaxt = "n", outline = F, xlab="", 
              at=c(1,2,2.75,3.25,4,5), boxwex = 0.45, ylim=c(0, ymax), yaxt="n", ylab="", xlim=c(0.5, 5.5))
      par(new = TRUE)
      boxplot(values ~ group, data_cc_other, las=1, col = sel_cols[3], varwidth = F, xaxt = "n", outline = F, xlab="", 
              at=c(1,2,3,4,4.75,5.25), boxwex = 0.45, ylim=c(0, ymax), yaxt="n", ylab="", xlim=c(0.5, 5.5))
      
      axis(1,  at=c(1,2,2.75,3.25,4,4.75,5.25), labels = c("healthy", "other cancer",
                                                       "other cancer\nno comorbidity",
                                                       "other cancer\ncomorbidity", 
                                                       "cancer class",
                                                       "cancer class\nno comorbidity",
                                                       "cancer class\ncomorbidity"), las=2)
      
      mtext("pathway abundance", 2, line = 2.5)
      mtext(paste0(cancer_class, " _ ", elix_row_names[j], "\n", sig_pathways[k]), 3, line = 0, cex=1)
      
    }
  }
  dev.off()
}



#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#make plots for all interesting species as well without factoring in comorbidities

#[1] "esophagus"   "liver_intrahepatic_bile_ducts" "neuroendocrine_tumors"  "ovary"   "pancreas" 

cancer_row_list <- list(
  grep("esophagus", group_vec),
  grep("liver and intrahepatic bile ducts", group_vec),
  grep("neuroendocrine", group_vec),
  grep("ovary", group_vec),
  grep("pancreas", group_vec))


names(merge_list_sig_Q_w_comorb)

for (i in 1:length(merge_list_sig_Q_w_comorb)) { #cancer class
  sig_df <- merge_list_sig_Q_w_comorb[[i]]
  
  cancer_group_rows <- cancer_row_list[[i]]
  all_other_cancer_rows <- which(!1:nrow(path_table) %in% c(healthy_rows, cancer_group_rows))
  cancer_class <- names(merge_list_sig_Q_w_comorb)[i]
  
  sig_pathways <- rownames(sig_df)
  
  pdf(paste0("./Figure/RM_figures/Significant pathways_scatter_", cancer_class, ".pdf"), width=5, height=4.5)
  for (k in 1:length(sig_pathways)) {
    
    temp_ra <- as.numeric(path_table[sig_pathways[k],])
    plot_df <- as.data.frame(cbind(c(rep("healthy", times=length(healthy_rows)), rep("all other cancer", times=length(all_other_cancer_rows)), rep("cancer_group", times=length(cancer_group_rows))),
          c(temp_ra[healthy_rows], temp_ra[all_other_cancer_rows], temp_ra[cancer_group_rows])))
    names(plot_df) <- c("group", "ra")
    plot_df$ra <- as.numeric(plot_df$ra)
    
    p_plot <-  ggboxplot(plot_df, "group", "ra",
                         color = "group", palette=sel_cols,
                         add = "jitter", ylab= "sqrt(relative_abundance %)") +
      scale_y_continuous(trans='sqrt') +
      ggtitle(sig_pathways[k]) +
      theme(legend.position="none")
    
    print(p_plot)
  }
  dev.off()
}



#-------------------------------------------------------------------------------
#For all cancer classes and a single significance level <0.2. Show number of pathways vs control, vs other cancers, and overlap.
#And how many after filtering the cancer-specific comorbidities


cutoff <- 0.2

merge_list <- vector("list", length = length(canX_Ctr_list))
names(merge_list) <- names(canX_Ctr_list)

sig_control_list <- vector("list", length = length(canX_Ctr_list))
sig_cancer_list <- vector("list", length = length(canX_Ctr_list))

#merge these in a for loop
for (i in 1:length(canX_Ctr_list)) {
  
  control <- canX_Ctr_list[[i]]
  cancer <- canX_Ex_list[[i]]
  
  idx.c <- control$pathway[which(control$Q<cutoff)]
  idx.d <- cancer$pathway[which(cancer$Q<cutoff)]
  idx.com <- intersect(idx.c, idx.d)
  
  print(idx.d)
  
  control_sub <- control[idx.com,c("P", "R2", "Q")]
  cancer_sub <- cancer[idx.com,c("P", "R2", "Q")]
  
  sig_control_list[[i]] <- idx.c
  sig_cancer_list[[i]] <- idx.d
  
  #make sure directionality is the same
  control_sub_overlap <- control_sub[(control_sub$R2 < 0) == (cancer_sub$R2 < 0),]
  cancer_sub_overlap <- cancer_sub[(cancer_sub$R2 < 0) == (control_sub$R2 < 0),]
  
  temp_list <- list(control_sub_overlap, cancer_sub_overlap)
  names(temp_list) <- c("Ctr", "Ex")
  
  if (min(c(nrow(control_sub), nrow(cancer_sub))) > 1) {
    merge_list[[i]] <- do.call(cbind, temp_list)
  }
}


#add number after removing species-specific comorbidities
n_wo_sig_elix <- rep(0, times= length(names(merge_list)))
names(n_wo_sig_elix) <- names(merge_list)

#"esophagus"   "liver_intrahepatic_bile_ducts" "neuroendocrine_tumors"  "ovary"   "pancreas" 
n_wo_sig_elix["esophagus"] <- 3
n_wo_sig_elix["liver_intrahepatic_bile_ducts"] <- 2 #is this correct?
n_wo_sig_elix["neuroendocrine_tumors"] <- 50
n_wo_sig_elix["ovary"] <- 10
n_wo_sig_elix["pancreas"] <- 0

plot_df <- cbind(names(merge_list), 
      as.character(lapply(sig_control_list, length)), 
      as.character(lapply(sig_cancer_list, length)),
      as.character(lapply(merge_list, nrow)),
      as.numeric(n_wo_sig_elix)
      )

plot_df[plot_df == "NULL"] <- 0
plot_df


colnames(plot_df) <- c("Cancer class", "vs. Healthy", "vs. other Cancer", "overlap", "overlap w/o comorbidities")
plot_df <- as.data.frame(plot_df)

plot_df$`vs. Healthy` <- as.numeric(plot_df$`vs. Healthy`)
plot_df$`vs. other Cancer` <- as.numeric(plot_df$`vs. other Cancer`)
plot_df$overlap <- as.numeric(plot_df$overlap)
plot_df$`overlap w/o comorbidities` <- as.numeric(plot_df$`overlap w/o comorbidities`)

plot_df <- plot_df[order(plot_df$overlap, decreasing = T),]

gt(plot_df) %>%
  tab_header(title = "Number of species when correcting for comorbidities") %>%
  opt_vertical_padding(scale = 0.1)

#make a version with FDR 0.1 cutoff?


#-------------------------------------------------------------------------------
#number of pathways consistently found in the comparison with healthy controls
#get numbers for the main paper

cutoff <- 0.1

merge_list <- vector("list", length = length(canX_Ctr_list))
names(merge_list) <- names(canX_Ctr_list)

sig_control_list_2 <- vector("list", length = length(canX_Ctr_list))


#merge these in a for loop
for (i in 1:length(canX_Ctr_list)) {
  
  control <- canX_Ctr_list[[i]]
  idx.c <- control$pathway[which(control$Q<cutoff)]
  control_sub <- control[idx.com,c("P", "R2", "Q")]
  sig_control_list_2[[i]] <- idx.c
}

lapply(sig_control_list_2, length)
median(as.numeric(lapply(sig_control_list_2, length)))

length(unique(unlist(sig_control_list_2)))
#304

#23 cancer classes
12 /length(sig_control_list_2) #52.2%
length(which(table(unlist(sig_control_list_2)) >= 12))
#60 pathways


16 /length(sig_control_list_2) #69.6%
length(which(table(unlist(sig_control_list_2)) >= 16))
#19 pathways; make heatmap of these pathways
names(which(table(unlist(sig_control_list_2)) >= 16))

#3 UMP pathways and 3 pyrimidine; all linked to pyrimidines. The UMP ones are slightly down

#Cells may compensate for reduced UMP levels by increasing the overall flux through the pyrimidine biosynthesis 
#pathway to ensure sufficient production of other pyrimidine nucleotides like cytidine triphosphate (CTP) and 
#thymidine triphosphate (TTP). This is particularly relevant in rapidly dividing cells, such as cancer cells, 
#that require high levels of nucleotides for DNA replication.


#[1] "COA-PWY-1: superpathway of coenzyme A biosynthesis III (mammals)"                              
#[2] "GALACTARDEG-PWY: D-galactarate degradation I"                                                  
#[3] "GLUCARGALACTSUPER-PWY: superpathway of D-glucarate and D-galactarate degradation"              
#[4] "NONMEVIPP-PWY: methylerythritol phosphate pathway I"                                           
#[5] "PANTO-PWY: phosphopantothenate biosynthesis I"                                                 
#[6] "PANTOSYN-PWY: superpathway of coenzyme A biosynthesis I (bacteria)"                            
#[7] "PWY-5686: UMP biosynthesis I"                                                                  
#[8] "PWY-6163: chorismate biosynthesis from 3-dehydroquinate"                                       
#[9] "PWY-6385: peptidoglycan biosynthesis III (mycobacteria)"                                       
#[10] "PWY-6545: pyrimidine deoxyribonucleotides de novo biosynthesis III"                            
#[11] "PWY-6700: queuosine biosynthesis I (de novo)"                                                  
#[12] "PWY-7184: pyrimidine deoxyribonucleotides de novo biosynthesis I"                              
#[13] "PWY-7221: guanosine ribonucleotides de novo biosynthesis"                                      
#[14] "PWY-7790: UMP biosynthesis II"                                                                 
#[15] "PWY-7791: UMP biosynthesis III"                                                                
#[16] "PWY-7953: UDP-N-acetylmuramoyl-pentapeptide biosynthesis III (meso-diaminopimelate containing)"
#[17] "PWY0-1298: superpathway of pyrimidine deoxyribonucleosides degradation"                        
#[18] "PWY0-1477: ethanolamine utilization"                                                           
#[19] "RIBOSYN2-PWY: flavin biosynthesis I (bacteria and plants)"    


#16 / 23 as cutoff
data_list <- list()
for (i in 1:length(canX_Ctr_list)) {
  control <- canX_Ctr_list[[i]]
  control_sub <- control[names(which(table(unlist(sig_control_list_2)) >= 16)),c("P", "R2", "Q")]
  data_list[[i]] <- control_sub
}

R2_df <- do.call(cbind, lapply(data_list, function(x) x$R2))
Q_df <- do.call(cbind, lapply(data_list, function(x) x$Q))

colnames(R2_df) <- names(canX_Ctr_list)
colnames(Q_df) <- names(canX_Ctr_list)

rownames(R2_df) <- names(which(table(unlist(sig_control_list_2)) >= 16))
rownames(Q_df) <- names(which(table(unlist(sig_control_list_2)) >= 16))

R2_df[is.na(R2_df)] <- 0
Q_df[is.na(Q_df)] <- 1

dim(R2_df)
# 19 23


#------------------------------------
#make heatmap out of this

dim(R2_df); dim(Q_df)
colnames(R2_df); rownames(R2_df)

Q.sig <- t(Q_df)
R2.sig <- t(R2_df)

#scale_max <- max(c(abs(min(R2.sig)), max(R2.sig)))
rownames(R2.sig) <- tolower(gsub("_", " ", rownames(R2.sig)))
rownames(R2.sig)[rownames(R2.sig) == "plasma cell"] <- "multiple myeloma plasma cell"


sign_cutoffs <- c(0.2, 0.1, 0.05)
name_vec = colnames(FDR_df)
max.cutoff <- max(sign_cutoffs)
q.cut1 <- min(sign_cutoffs); q.cut2 <- median(sign_cutoffs); q.cut3 <- max.cutoff

grid.text_1 <- '***'
grid.text_2 <- '**'
grid.text_3 <- '*'

scale_max <- max(c(abs(min(R2.sig)), max(R2.sig)))


#elevated pyrimidine synthesis pathways in many cancer classes is interesting

png(paste0('./Figure/RM_figures/heatmap_consistent_pathways.png'), width =11*330, height = 12*330, res=300)
Heatmap(R2.sig, name = "R2 x effect direction", 
        heatmap_legend_param = list(direction = "vertical", title="R2"),
        col = colorRamp2(c(min(R2.sig), 0, max(R2.sig)), c("#0072B2",'white', "#CC79A7")),
        #col = colorRamp2(c(-scale_max, 0, scale_max), c("#0072B2",'white', "#CC79A7")),
        column_gap = unit(1, "mm"), 
        column_title_rot = 30,
        column_names_rot = 90,
        column_names_max_height = max_text_width(colnames(R2.sig)),
        border = TRUE,
        cell_fun = function(j, i, x, y, w, h, fill) {
          if(Q.sig[i,j] < q.cut1) {
            grid.text('***', x = x, y = y, r = unit(1/30,'cm'))
          }else if(Q.sig[i,j] < q.cut2){
            grid.text('**', x = x, y = y, r = unit(1/30,'cm'))
          }else if(Q.sig[i,j] < q.cut3){
            grid.text('*', x = x, y = y,r = unit(1/30,'cm'))
          }
        },
        rect_gp = gpar(col= "white"),
        show_column_names = T, 
        show_row_names = T,
        show_column_dend = F,
        show_row_dend = F,
        cluster_rows = T,
        cluster_columns = FALSE,
        column_names_gp = gpar(fontface = 'italic'),
        row_names_max_width = max_text_width(rownames(R2.sig)),
        row_names_gp = gpar(fontsize = 14)
) 
dev.off()











