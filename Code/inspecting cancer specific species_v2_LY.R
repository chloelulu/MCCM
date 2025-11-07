require(GUniFrac)
require(tidyverse)
require(RColorBrewer)
library(gt)
require(ComplexHeatmap)
require(circlize)
require(ggpubr)
require(openxlsx)


file_dir <- dirname(rstudioapi::getSourceEditorContext()$path)
file_dir <- dirname(dirname(dirname(dirname(file_dir))))

setwd(file_dir) #"mforge_clean"

load(file = 'Data/data.obj.raw.core.RData') 
clin_meta <- as.data.frame(data.obj$meta.dat)
tax_table <- as.data.frame(data.obj$otu.tab)

dim(tax_table) #7839 1651


#3095 above 0 in more than 5% of the samples
table(rowSums(tax_table > 0) > ncol(tax_table) * 0.05)
#FALSE  TRUE 
#4744  3095


dim(data.obj$abund.list$Phylum) #41
dim(data.obj$abund.list$Family) #391
dim(data.obj$abund.list$Genus) #1899


#load pathway abundances
load(file = 'Data/data.obj.pathway.RData') 
path_table <- as.data.frame(data.obj$otu.tab)

dim(path_table) #519 1651


#make function for heatmap from dataframe
source("Code/Submission/MayoOncobiomeStudy/Code/zicoseq_heatmap_from_df.R")

#order the clin meta based on tax_table
clin_meta <- clin_meta[colnames(tax_table),]



#-------------------------------------------------------------------------------

load("Figure/subCancerX_Control/DAA_P_R2.RData")
zico_res_ctr <- list(P.All$Species, Q.All$Species, R2.All$Species)
names(zico_res_ctr) <- c("P", "Q", "R2")


load("Figure/subCancerX-Ex/DAA_P_R2.RData")
zico_res_ex_cancer <- list(P.All$Species, Q.All$Species, R2.All$Species)
names(zico_res_ex_cancer) <- c("P", "Q", "R2")


temp_dir <- "./Result/CancerOnly/DAA/"
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

elix_variables_filenames <- list.files(paste0(temp_dir, elix_variables), pattern = 'Rdata')
elix_variables_filenames_full <- paste0(temp_dir, elix_variables, "/", elix_variables_filenames)

elix_variables_diff_list <- list()
for (i in 1:length(elix_variables_filenames_full)) {
  load(elix_variables_filenames_full[i])
  elix_variables_diff_list[[i]] <- diff.obj
}

names(elix_variables_diff_list) <- elix_variables
length(elix_variables_diff_list) #22


#-----------------------------------
#check esophagus for E.coli

#zico_res_ex_cancer$Q$esophagus[grep("Escherichia", zico_res_ex_cancer$Q$Species)]
#zico_res_ctr$Q$`Control-esophagus`[grep("Escherichia", zico_res_ctr$Q$Species)]

#zico_res_ctr$Q[grep("Escherichia", zico_res_ctr$Q$Species),]
#zico_res_ctr$R2[grep("Escherichia", zico_res_ctr$R2$Species),]


  
#-------------------------------------------------------------------------------
#get the relevant information out, only the R2, coef, q value

#these variables are corrected for cancer class; does that make sense to do here?

length(elix_variables_diff_list$Elixhauser_Alcohol)
names(elix_variables_diff_list$Elixhauser_Alcohol)
#"pv.list"   "fc.list"   "pc.list"   "qv.list"   "R2.list"   "coef.list" "m.list"


names(elix_variables_diff_list$Elixhauser_Alcohol$R2.list)
head(elix_variables_diff_list$Elixhauser_Alcohol$R2.list$Species, 3)
#                             Func1
#s__14-2 sp000403255          2.380211e-05
#s__1XD42-69 sp014287635      7.957575e-04
#s__Absicoccus porci          8.113269e-06


elix_variables_diff_list_Q <- lapply(elix_variables_diff_list, function(x) x$qv.list$Species[,'Qvalue'])
elix_variables_diff_list_R2 <- lapply(elix_variables_diff_list, function(x) x$R2.list$Species[,'Func1'])
elix_variables_diff_list_coef <- lapply(elix_variables_diff_list, function(x) x$coef.list$Species)
coef_list_for_sign <- lapply(elix_variables_diff_list_coef, function(x) x[,ncol(x)])


elix_variables_diff_list_R2_sign <- elix_variables_diff_list_R2
for (i in 1:length(elix_variables_diff_list_R2)) {
  elix_temp <- elix_variables_diff_list_R2[[i]]
  elix_variables_diff_list_R2_sign[[i]] <- elix_temp * sign(coef_list_for_sign[[i]])
}


#subset for species in all comparisons
temp_table <- table(unlist(lapply(elix_variables_diff_list_R2_sign, names))) == 22
names_to_include <- names(temp_table[temp_table == TRUE])

elix_variables_diff_list_R2_sign <- lapply(elix_variables_diff_list_R2_sign, function(x) x[names_to_include])
elix_variables_diff_list_Q <- lapply(elix_variables_diff_list_Q, function(x) x[names_to_include])

R2_sign_df_comorbidities <- do.call(cbind, elix_variables_diff_list_R2_sign)
Q_sign_df_comorbidities <- do.call(cbind, elix_variables_diff_list_Q)

dim(R2_sign_df_comorbidities); dim(Q_sign_df_comorbidities)
#1236   22


#-------------------------------------------------------------------------------
#find strongest species linked to weightloss

sig_names <- names(which(Q_sign_df_comorbidities[,"Elixhauser_WeightLoss"] < 0.1))
names(head(sort(Q_sign_df_comorbidities[sig_names,"Elixhauser_WeightLoss"], decreasing = T), 40))


#-------------------------------------------------------------------------------
#also get the general Elixhauser associations

load("./Result/CancerOnly/DAA/Elix_score/Elix_score_ZicoSeq.Rdata")
Elix_res_list <- diff.obj

names(Elix_res_list$R2.list)
head(Elix_res_list$R2.list$Species,3)
#                             Func1
#s__14-2 sp000403255     7.960133e-05
#s__1XD42-69 sp014287635 8.107676e-05
#s__Absicoccus porci     6.887762e-05

elix_diff_Q <- Elix_res_list$qv.list$Species[,'Qvalue']
elix_diff_R2 <- Elix_res_list$R2.list$Species[,'Func1']
elix_diff_coef <- Elix_res_list$coef.list$Species[,"Elix_score"]

elix_diff_R2 <- elix_diff_R2 * sign(elix_diff_coef)
length(elix_diff_R2[elix_diff_Q < 0.1])
#318
table(elix_diff_R2[elix_diff_Q < 0.1] > 0)
#FALSE (negatively associated with ECI)  TRUE (positively associated with ECI)
#                 229                                      89 


#to inspect the species that are positive with ECI
names(which((elix_diff_R2[elix_diff_Q < 0.1] > 0) == TRUE))


#-------------------------------------
#get number of species FDR < 0.1 with positive and negative effect size for components

table_list <- list()
for (i in 1:ncol(R2_sign_df_comorbidities)) {
  R2_vec_temp <- R2_sign_df_comorbidities[,i]
  table_list[[i]] <- table(R2_vec_temp[Q_sign_df_comorbidities[,i] < 0.1] > 0)
}
names(table_list) <- colnames(R2_sign_df_comorbidities)

#TRUE means positively associated with the ECI component


#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#load all the control vs cancer class and cancer class vs other cancers comparisons

canX_Ctr <- read.csv("./Figure/subCancerX_Control/Species_DAA.csv")
#head(canX_Ctr)
canX_Ctr$variable <- gsub("Control-", "", canX_Ctr$variable)

canX_Ctr_list <- lapply(as.list(unique(canX_Ctr$variable)), function(x) canX_Ctr[canX_Ctr$variable == x,]) 
names(canX_Ctr_list) <- unique(canX_Ctr$variable)
for (i in 1:length(canX_Ctr_list)) {
  rownames(canX_Ctr_list[[i]]) <- canX_Ctr_list[[i]]$Species
}


canX_Ex <- read.csv("./Figure/subCancerX-Ex/Species_DAA.csv")
canX_Ex_list <- lapply(as.list(unique(canX_Ex$variable)), function(x) canX_Ex[canX_Ex$variable == x,]) 
names(canX_Ex_list) <- unique(canX_Ex$variable)
#head(canX_Ex_list)
for (i in 1:length(canX_Ex_list)) {
  rownames(canX_Ex_list[[i]]) <- canX_Ex_list[[i]]$Species
}

#for every variable this needs merged separately, I want to make a heatmap for every cancer class separately
length(canX_Ctr_list); length(canX_Ex_list) #22 cancer classes
#cbind(names(canX_Ctr_list), names(canX_Ex_list))


#-------------------------------------------------------------------------------
#get numbers only for the cancer class comparisons vs other cancers without considering the vs healthy comparisons
#with positive and negative
cutoff <- 0.1

#canX_Ex_list
res_list_canX <- lapply(canX_Ex_list, function(x) table(x$Q <= cutoff))
res_list_canX_direction <- lapply(canX_Ex_list, function(x) table(x$R2[which(x$Q <= cutoff)] > 0))


res_list_canX_direction
#  $esophagus
#TRUE 
#15 

#  $`liver intrahepatic bile ducts`
#FALSE  TRUE 
#87    19 

#$`lymphoid leukemia`
#FALSE  TRUE 
#4    43 

#  $`multiple myeloma plasma cell`
#FALSE  TRUE 
#1     3 

#$`neuroendocrine tumors`
#FALSE  TRUE 
#313   125 

#  $pancreas
#TRUE 
#2 

#$prostate
#FALSE  TRUE 
#3    39 




#-------------------------------------------------------------------------------
#only select those species significant <0.1 FDR in both control and cancer class comparisons

#trying this discovery level cutoff but mark significance in the heatmap
cutoff <- 0.1

merge_list <- vector("list", length = length(canX_Ctr_list))
names(merge_list) <- names(canX_Ctr_list)

#merge these in a for loop
for (i in 1:length(canX_Ctr_list)) {
  
  control <- canX_Ctr_list[[i]]
  cancer <- canX_Ex_list[[i]]
  
  ## find both significant in control and cancer
  idx.c <- control$Species[which(control$Q<cutoff)]
  idx.d <- cancer$Species[which(cancer$Q<cutoff)]
  idx.com <- intersect(idx.c, idx.d)
  
  control_sub <- control[idx.com,c("P", "R2", "Q")]
  cancer_sub <- cancer[idx.com,c("P", "R2", "Q")]
  
  #make sure directionality is the same
  control_sub_overlap <- control_sub[(control_sub$R2 < 0)  ==  (cancer_sub$R2 < 0),,drop =F]
  cancer_sub_overlap <- cancer_sub[(cancer_sub$R2 < 0) == (control_sub$R2 < 0),,drop =F]
  
  temp_list <- list(control_sub_overlap, cancer_sub_overlap)
  names(temp_list) <- c("Ctr", "Ex")
  
  if (min(c(nrow(control_sub), nrow(cancer_sub))) > 0) {
    merge_list[[i]] <- do.call(cbind, temp_list)
  }
}
merge_list_sig <- merge_list[which(lapply(merge_list, function(x) nrow(x) > 0) == TRUE)]

lapply(merge_list_sig, dim)


merge_list_sig_Q <- lapply(merge_list_sig, function(x) x[,c("Ctr.Q", "Ex.Q")])
merge_list_sig_R2 <- lapply(merge_list_sig, function(x) x[,c("Ctr.R2", "Ex.R2")])

lapply(merge_list_sig_Q, dim)

#TRUE indicates positive effect size, higher in cancer
lapply(merge_list_sig_R2, function(x) table(x$Ctr.R2 > 0))

#$esophagus
#TRUE 
# 7 

#$`liver intrahepatic bile ducts`
#FALSE  TRUE 
# 64    10 

#$`lymphoid leukemia`
#FALSE  TRUE 
# 4    29 

#$`multiple myeloma plasma cell`
#FALSE  TRUE 
# 1     3 

#$`neuroendocrine tumors`
#FALSE  TRUE 
# 220    69 

lapply(merge_list_sig_Q, function(x) length(which(x$Ctr.Q < 0.05)))
lapply(merge_list_sig_Q, function(x) length(which(x$Ctr.Q < 0.01)))


fdr_thresholds <- c(0.1, 0.05, 0.01)
summary_list <- lapply(names(merge_list_sig_Q), function(name) {
  df <- merge_list_sig_Q[[name]]
  counts <- sapply(fdr_thresholds, function(thresh) sum(df[, "Ex.Q"] < thresh, na.rm = TRUE))
  return(c(name = name, counts))
})

df <- do.call(rbind, summary_list)
df <- as.data.frame(df)
colnames(df) <- c("Name", paste0("n.species.FDR.<", fdr_thresholds))

# Convert numeric columns
df[ , -1] <- lapply(df[ , -1], as.numeric)

df <- df[order(rowSums(df[,2:ncol(df)]), decreasing = T),]
gt(df) %>%
  tab_header(title = "Number of species overlapping vs. healthy and vs. other cancers") %>%
  cols_label(Name="Name", 
            `n.species.FDR.<0.1`="n species q<0.1", 
            `n.species.FDR.<0.05`="n species q<0.05", 
             `n.species.FDR.<0.01`="n species q<0.01") %>%
  opt_vertical_padding(scale = 0.1)



#-------------------------------------------------------------------------------
#add the comorbidities from above

merge_list_sig_Q_w_comorb <- merge_list_sig_Q
merge_list_sig_R2_w_comorb <- merge_list_sig_R2

for (i in 1:length(merge_list_sig_Q)) {
  #Q
  merge_list_temp <- merge_list_sig_Q[[i]]
  merge_list_temp_sub <- merge_list_temp[rownames(merge_list_temp) %in% rownames(Q_sign_df_comorbidities),]
  merge_list_temp_sub_comorb <- cbind(merge_list_temp_sub, Q_sign_df_comorbidities[rownames(merge_list_temp_sub),, drop =F])
  merge_list_sig_Q_w_comorb[[i]] <- merge_list_temp_sub_comorb
  
  #R2
  merge_list_temp <- merge_list_sig_R2[[i]]
  merge_list_temp_sub <- merge_list_temp[rownames(merge_list_temp) %in% rownames(R2_sign_df_comorbidities),]
  merge_list_temp_sub_comorb <- cbind(merge_list_temp_sub, R2_sign_df_comorbidities[rownames(merge_list_temp_sub),, drop =F])
  merge_list_sig_R2_w_comorb[[i]] <- merge_list_temp_sub_comorb

}


#-------------------------------------------------------------------------------

lapply(merge_list_sig_Q_w_comorb, dim)

load('./Result/CancerOnly/data.obj.wk.RData')
level <- 'Family'
otu.name <- as.data.frame(data.obj$otu.name)


names(merge_list_sig_Q_w_comorb)
sign_cutoffs=c(0.1,0.05, 0.01)



#--------------------------------------
i="esophagus"
# "esophagus"                     "liver intrahepatic bile ducts" "lymphoid leukemia"             "multiple myeloma plasma cell"  "neuroendocrine tumors"  
FDR_df <- merge_list_sig_Q_w_comorb[[i]]
R2_df <- merge_list_sig_R2_w_comorb[[i]]
dim(FDR_df);dim(R2_df)
grp.labels <- colnames(FDR_df)
grp.labels[1:2] <- c("Cancer class vs. healthy", "Cancer class vs. other cancers")

heatmap_list <- zicoseq_heatmap(FDR_df, R2_df, grp.names, grp.labels, col = col,sign_cutoffs=sign_cutoffs, single_cutoff = FALSE, otu.name=otu.name)
heatmap_list$heatmap_plot

filename <- "./Figure/RM_figures/overlap heatmap esophagus all.png"
png(filename, width = 7.5*330, height = 10*330, res = 330)
heatmap_list$heatmap_plot
dev.off()


#--------------------------------------
i="liver intrahepatic bile ducts"

FDR_df <- merge_list_sig_Q_w_comorb[[i]]
R2_df <- merge_list_sig_R2_w_comorb[[i]]

heatmap_list <- zicoseq_heatmap(FDR_df, R2_df, grp.names, grp.labels, col = col,sign_cutoffs, single_cutoff = FALSE, otu.name=otu.name)
heatmap_list$heatmap_plot

filename <- "./Figure/RM_figures/overlap heatmap liver_intrahepatic_bile_ducts all.png"
png(filename, width = 25*330, height = 10*330, res = 330)
heatmap_list$heatmap_plot
dev.off()



#--------------------------------------
i="lymphoid leukemia"

FDR_df <- merge_list_sig_Q_w_comorb[[i]]
R2_df <- merge_list_sig_R2_w_comorb[[i]]
grp.names <- colnames(FDR_df)
grp.labels <- colnames(FDR_df)
grp.labels[1:2] <- c("Cancer class vs. healthy", "Cancer class vs. other cancers")

heatmap_list <- zicoseq_heatmap(FDR_df, R2_df, grp.names, grp.labels, col = col,sign_cutoffs, single_cutoff = FALSE, otu.name=otu.name)
heatmap_list$heatmap_plot

filename <- "./Figure/RM_figures/overlap heatmap Lymphoid_leukemia all.png"
png(filename, width =15*330, height = 10*330, res = 330)
heatmap_list$heatmap_plot
dev.off()



#--------------------------------------
i="multiple myeloma plasma cell" 

FDR_df <- merge_list_sig_Q_w_comorb[[i]]
R2_df <- merge_list_sig_R2_w_comorb[[i]]

heatmap_list <- zicoseq_heatmap(FDR_df, R2_df, grp.names, grp.labels, col = col,sign_cutoffs, single_cutoff = FALSE, otu.name=otu.name)
heatmap_list$heatmap_plot

filename <- "./Figure/RM_figures/overlap heatmap multiple myeloma plasma cell all.png"
png(filename, width = 6.75*330, height = 10*330, res = 330)
heatmap_list$heatmap_plot
dev.off()



#--------------------------------------
i="neuroendocrine tumors"
names(merge_list_sig_Q_w_comorb)
FDR_df <- merge_list_sig_Q_w_comorb[[i]]
R2_df <- merge_list_sig_R2_w_comorb[[i]]

heatmap_list <- zicoseq_heatmap(FDR_df, R2_df, grp.names, grp.labels, col = col,sign_cutoffs, single_cutoff = FALSE, otu.name=otu.name)
heatmap_list$heatmap_plot

filename <- "./Figure/RM_figures/overlap heatmap neuroendocrine_tumors all.png"
png(filename, width = 65*330, height = 10*330, res = 330)
heatmap_list$heatmap_plot
dev.off()



#-------------------------------------------------------------------------------
#cancer classes of interest

cancer_class_w_sign <- c("Malignant neoplasm of esophagus", 
                         "Malignant neoplasm of liver and intrahepatic bile ducts",
                         "Lymphoid leukemia", 
                         "Multiple myeloma and malignant plasma cell neoplasms",
                         "Malignant neuroendocrine tumors")

cancer_class_w_sign_row_list <- sapply(cancer_class_w_sign, function(x) which(clin_meta$icd10_first_3_name == x))

#tabulate all Elixhauser components with a significant species
test_df <- lapply(merge_list_sig_Q_w_comorb, function(x) (x[,3:ncol(x)] < 0.1)) #0.2 is FDR cutoff

#general one
Elix_components_w_sig <- unique(unlist(lapply(test_df, function(x) colnames(x)[colSums(x == TRUE) > 0])))


#make Elixhauser components specific to that cancer class
Elix_components_w_sig_list <- lapply(test_df, function(x) colnames(x)[colSums(x == TRUE) > 0])


component_list <- list()
cancer_class_w_sign_row_list_sub <- cancer_class_w_sign_row_list[lapply(Elix_components_w_sig_list, length) > 0] 
Elix_components_w_sig_list_sub <- Elix_components_w_sig_list[lapply(Elix_components_w_sig_list, length) > 0]


#-------------------------------------------------------------------------------
#for cancer groups that have only very few patients with these comorbidities these comorbidities do not have to be plotted or considered
#only consider more than >15%

#-------------------------------------------------------------------------------
#"Malignant neoplasm of esophagus"
i=1
names(cancer_class_w_sign_row_list_sub)[i]
length(cancer_class_w_sign_row_list_sub[[i]]) #34
temp_table <- apply(clin_meta[cancer_class_w_sign_row_list_sub[[i]], Elix_components_w_sig_list_sub[[i]]], 2, table)

#get the percentages
temp_table / length(cancer_class_w_sign_row_list_sub[[i]]) * 100

#>15%
names(which((temp_table / length(cancer_class_w_sign_row_list_sub[[i]]) * 100)[2,]>15))

selected_elix_esop <- c("Elixhauser_FluidsLytes")

#----------------------------------
#"Malignant neoplasm of liver and intrahepatic bile ducts"
i=2
names(cancer_class_w_sign_row_list_sub)[i]
length(cancer_class_w_sign_row_list_sub[[i]]) #39
temp_table <- do.call(cbind, apply(clin_meta[cancer_class_w_sign_row_list_sub[[i]], Elix_components_w_sig_list_sub[[i]]], 2, table))
temp_table / length(cancer_class_w_sign_row_list_sub[[i]]) * 100


#liver is obvious so good example to be plotted
names(which((temp_table / length(cancer_class_w_sign_row_list_sub[[i]]) * 100)[2,]>15))

selected_elix_liv <- c("Elixhauser_Arrhythmia", "Elixhauser_DMcx","Elixhauser_FluidsLytes", "Elixhauser_Liver", "Elixhauser_WeightLoss")


#-------------------------------------------------------------------------------
#"Lymphoid leukemia"
i=3
names(cancer_class_w_sign_row_list_sub)[i]
length(cancer_class_w_sign_row_list_sub[[i]]) #17
temp_table <- apply(clin_meta[cancer_class_w_sign_row_list_sub[[i]], Elix_components_w_sig_list_sub[[i]]], 2, table)
#get the percentages
temp_table / length(cancer_class_w_sign_row_list_sub[[i]]) * 100

#Remove Lymphoma all together since we have another paper on this?
#Arrhythmia has the only decent number of patients that seem relevant.
names(which((temp_table / length(cancer_class_w_sign_row_list_sub[[i]]) * 100)[2,]>15))
#>15%
selected_elix_lymp <- c("Elixhauser_Arrhythmia", "Elixhauser_FluidsLytes")


#----------------------------------
#"Multiple myeloma and malignant plasma cell neoplasms"
i=4
names(cancer_class_w_sign_row_list_sub)[i]
length(cancer_class_w_sign_row_list_sub[[i]]) #60
temp_table <- apply(clin_meta[cancer_class_w_sign_row_list_sub[[i]], Elix_components_w_sig_list_sub[[i]]], 2, table)
temp_table / length(cancer_class_w_sign_row_list_sub[[i]]) * 100

selected_elix_plasmacell <- c("Elixhauser_FluidsLytes")


#----------------------------------
#"Malignant neuroendocrine tumors"

i=5
length(cancer_class_w_sign_row_list_sub[[i]]) #21
temp_table <-  do.call(cbind, apply(clin_meta[cancer_class_w_sign_row_list_sub[[i]], Elix_components_w_sig_list_sub[[i]]], 2, table))

temp_table / length(cancer_class_w_sign_row_list_sub[[i]]) * 100

selected_elix_neuro <-   c("Elixhauser_Arrhythmia", "Elixhauser_Liver")



#-------------------------------------------------------------------------------
#plot using the subset of the Elixhauser categories
names(cancer_class_w_sign_row_list_sub)
Elix_components_w_sig_list_sub_cutoff <- Elix_components_w_sig_list_sub

#esophagus
Elix_components_w_sig_list_sub_cutoff[["esophagus"]] <- Elix_components_w_sig_list_sub[["esophagus"]][Elix_components_w_sig_list_sub[["esophagus"]] %in% 
                                                                                    selected_elix_esop]
#liver_intrahepatic_bile_ducts
Elix_components_w_sig_list_sub_cutoff[["liver intrahepatic bile ducts"]] <- Elix_components_w_sig_list_sub[["liver intrahepatic bile ducts"]][Elix_components_w_sig_list_sub[["liver intrahepatic bile ducts"]] %in% 
                                                                                    selected_elix_liv]
#Lymphoid_leukemia
Elix_components_w_sig_list_sub_cutoff[["lymphoid leukemia"]] <- Elix_components_w_sig_list_sub[["lymphoid leukemia"]][Elix_components_w_sig_list_sub[["lymphoid leukemia"]] %in% 
                                                                                    selected_elix_lymp]
#plasma cell
Elix_components_w_sig_list_sub_cutoff[["multiple myeloma plasma cell"]] <- Elix_components_w_sig_list_sub[["multiple myeloma plasma cell"]][Elix_components_w_sig_list_sub[["multiple myeloma plasma cell"]] %in% 
                                                                                                                                              selected_elix_plasmacell]
#neuroendocrine_tumors
Elix_components_w_sig_list_sub_cutoff[["neuroendocrine tumors"]] <- Elix_components_w_sig_list_sub[["neuroendocrine tumors"]][Elix_components_w_sig_list_sub[["neuroendocrine tumors"]] %in% 
                                                                                    selected_elix_neuro]



#--------------------------------------


sign_lists_to_plot <- list()
R2_lists_to_plot <- list()


#--------------------------------------
i='esophagus'
#esophagus
names(merge_list_sig_Q_w_comorb)[i]
#"Elixhauser_Liver" "Elixhauser_Renal"
sig_components <- Elix_components_w_sig_list_sub_cutoff[[i]]
FDR_df <- merge_list_sig_Q_w_comorb[[i]][,c("Ctr.Q","Ex.Q",sig_components)]
R2_df <- merge_list_sig_R2_w_comorb[[i]][,c("Ctr.R2","Ex.R2",sig_components)]

grp.names <- colnames(FDR_df)
grp.labels <- colnames(FDR_df)
grp.labels[1:2] <- c("Esophagus vs. healthy", "Esophagus vs. other cancers")
sign_cutoffs <- c(0.1, 0.05, 0.01)

heatmap_list <- zicoseq_heatmap(FDR_df, R2_df, grp.names, grp.labels,col = col, sign_cutoffs, single_cutoff = FALSE, otu.name=otu.name, column_name_rot=40)
heatmap_list$heatmap_plot

filename <- paste0("./Figure/RM_figures/overlap heatmap ",i,"_sig_Elix.png")
png(filename, width = 7.75*330, height = 4.25*330, res = 330)
heatmap_list$heatmap_plot
dev.off()


#save as svg for allowing to edit colors
filename <- paste0("./Figure/RM_figures/overlap heatmap ",i,"_sig_Elix.svg")
svg(filename, width =7.25, height = 4)
heatmap_list$heatmap_plot
dev.off()


colnames(FDR_df) <- grp.labels
sign_lists_to_plot[[i]] <- FDR_df
R2_lists_to_plot[[i]] <- R2_df


#remove Elix species; any rows with FDR <0.1 beyond the first 2 columns
to_keep_rows <- !rowSums(FDR_df[, 3:ncol(FDR_df),drop =F] < 0.1) > 0
length(which(to_keep_rows)) #6

heatmap_list <- zicoseq_heatmap(FDR_df[to_keep_rows,], R2_df[to_keep_rows, ], grp.names, grp.labels, col = col,sign_cutoffs, single_cutoff = FALSE,
                                otu.name=otu.name, column_name_rot=40)

filename <- paste0("./Figure/RM_figures/overlap heatmap ",i,"_sig_Elix_rm.png")
png(filename, width = 7.75*330, height = 4.75*330, res = 330)
heatmap_list$heatmap_plot
dev.off()


#--------------------------------------
i="liver intrahepatic bile ducts"

sig_components <- Elix_components_w_sig_list_sub_cutoff[[i]]
FDR_df <- merge_list_sig_Q_w_comorb[[i]][,c("Ctr.Q","Ex.Q",sig_components)]
R2_df <- merge_list_sig_R2_w_comorb[[i]][,c("Ctr.R2","Ex.R2",sig_components)]

sum(rowSums(FDR_df[,c(1:2)]<0.1)>0) #74

grp.names <- colnames(FDR_df)
grp.labels <- colnames(FDR_df)
grp.labels[1:2] <- c("Cancer class vs. healthy", "Cancer class vs. other cancers")
sign_cutoffs <- c(0.1, 0.05, 0.01)



heatmap_list <- zicoseq_heatmap(FDR_df, R2_df, grp.names, grp.labels,col = col, sign_cutoffs, single_cutoff = FALSE, otu.name=otu.name, column_name_rot = 40)
heatmap_list$heatmap_plot

filename <- paste0("./Figure/RM_figures/overlap heatmap ",i,"_sig_Elix.png")
png(filename, width = 26*330, height = 6.5*330, res = 330)
heatmap_list$heatmap_plot
dev.off()


colnames(FDR_df) <- grp.labels
sign_lists_to_plot[[i]] <- FDR_df
R2_lists_to_plot[[i]] <- R2_df

#remove Elix species; any rows with FDR <0.1 beyond the first 2 columns
to_keep_rows <- !rowSums(FDR_df[, 3:ncol(FDR_df),drop =F] < 0.1) > 0
length(which(to_keep_rows)) #36



#make a separate version of this with only top 50% effect size
effect_size_cutoff <- 0.5
effect_cutoff <- min(sort(apply(abs(R2_df), 1, max), decreasing=T)[1:nrow(R2_df)*effect_size_cutoff])
effect_idx <- apply(abs(R2_df), 1, function(x) sum(x > effect_cutoff) > 0)

FDR_df_sub <- as.matrix(FDR_df[effect_idx,, drop =F])
R2_df_sub <- as.matrix(R2_df[effect_idx,, drop =F])

dim(FDR_df); dim(FDR_df_sub)



heatmap_list <- zicoseq_heatmap(FDR_df_sub, R2_df_sub, grp.names, grp.labels,col = col, sign_cutoffs, single_cutoff = FALSE, otu.name=otu.name, column_name_rot = 40)
heatmap_list$heatmap_plot

filename <- paste0("./Figure/RM_figures/overlap heatmap 50perc ",i,"_sig_Elix.png")
png(filename, width = 17.5*330, height = 5.75*330, res = 330)
heatmap_list$heatmap_plot
dev.off()



#make a separate version of this with a rescaled effect size; capped at 0.2
R2_df_rescaled <- R2_df_sub
R2_df_rescaled[c("s__Enterococcus faecalis", "s__Enterococcus_B durans"),"Ctr.R2"]
#[1] 0.11260284 0.03083123
R2_df_rescaled[c("s__Enterococcus faecalis", "s__Enterococcus_B durans"),"Ctr.R2"] <- c(0.02, 0.02)


heatmap_list <- zicoseq_heatmap(FDR_df_sub, R2_df_rescaled, grp.names, grp.labels,col = col, sign_cutoffs, single_cutoff = FALSE, otu.name=otu.name, column_name_rot = 40)

filename <- paste0("./Figure/RM_figures/overlap heatmap rescaled 50perc ",i,"_sig_Elix.png")
png(filename, width = 17.5*330, height = 5.75*330, res = 330)
heatmap_list$heatmap_plot
dev.off()



# remove Elix species; any rows with FDR <0.1 beyond the first 2 columns
#to_keep_rows <- !(rowSums(FDR_df[, 3:ncol(FDR_df),drop=F] < 0.1) > 0)
#length(which(to_keep_rows)) #35 # 40
#paste0(gsub('s__','',rownames(R2_df[to_keep_rows & R2_df$Ctr.R2>0, ])), collapse = ', ')
#heatmap_list <- zicoseq_heatmap(FDR_df[to_keep_rows,], R2_df[to_keep_rows, ], grp.names, grp.labels, col = col, sign_cutoffs, single_cutoff = FALSE,
#                                otu.name=otu.name, column_name_rot = 40)
#heatmap_list$heatmap_plot

#filename <- paste0("./Figure/RM_figures/overlap heatmap ",i,"_sig_Elix_rm.png")
#png(filename, width = 16*330, height = 6.5*330, res = 330)
#heatmap_list$heatmap_plot
#dev.off()


#------------------------------------------
#make version with positive and negative split. Save to svg to highlight species

FDR_df_pos <- FDR_df
R2_df_pos <- R2_df

FDR_df_pos <- FDR_df_pos[R2_df$Ctr.R2 > 0,]
R2_df_pos <- R2_df_pos[R2_df$Ctr.R2 > 0,]

dim(FDR_df_pos)

#make a version of this with only positive, no effect size cutoff
effect_size_cutoff <- 1
effect_cutoff <- min(sort(apply(abs(R2_df_pos), 1, max), decreasing=T)[1:nrow(R2_df_pos)*effect_size_cutoff])
effect_idx <- apply(abs(R2_df_pos), 1, function(x) sum(x > effect_cutoff) > 0)

FDR_df_sub <- as.matrix(FDR_df_pos[effect_idx,, drop =F])
R2_df_sub <- as.matrix(R2_df_pos[effect_idx,, drop =F])

dim(FDR_df_pos); dim(FDR_df_sub)

FDR_df_pos <- FDR_df_sub
R2_df_pos <- R2_df_sub

grp.labels <- c("Liver intrahepatic bile ducts vs. healthy", "Liver intrahepatic bile ducts vs. other cancers", 
                "Elixhauser_Arrhythmia", "Elixhauser_DMcx", "Elixhauser_FluidsLytes", 
                "Elixhauser_Liver", "Elixhauser_WeightLoss")

heatmap_list <- zicoseq_heatmap(FDR_df_pos, R2_df_pos, grp.names, grp.labels, col = col,sign_cutoffs, 
                                single_cutoff = FALSE, otu.name=otu.name, column_name_rot = 40, left_margin = 4)
heatmap_list$heatmap_plot


filename <- "./Figure/RM_figures/overlap heatmap liver positive.svg"
#png(filename, width = 9.3*330, height = 5.*330, res = 330)
svg(filename, width = 9.3, height = 5)
heatmap_list$heatmap_plot
dev.off()



#negative
FDR_df_neg <- FDR_df
R2_df_neg <- R2_df
FDR_df_neg <- FDR_df_neg[R2_df$Ctr.R2 < 0,]
R2_df_neg <- R2_df_neg[R2_df$Ctr.R2 < 0,]

dim(FDR_df_neg)

effect_size_cutoff <- 0.5
effect_cutoff <- min(sort(apply(abs(R2_df_neg), 1, max), decreasing=T)[1:nrow(R2_df_neg)*effect_size_cutoff])
effect_idx <- apply(abs(R2_df_neg), 1, function(x) sum(x > effect_cutoff) > 0)

FDR_df_sub <- as.matrix(FDR_df_neg[effect_idx,, drop =F])
R2_df_sub <- as.matrix(R2_df_neg[effect_idx,, drop =F])

dim(FDR_df_neg); dim(FDR_df_sub)

FDR_df_neg <- FDR_df_sub
R2_df_neg <- R2_df_sub


heatmap_list <- zicoseq_heatmap(FDR_df_neg, R2_df_neg, grp.names, grp.labels, col = col,sign_cutoffs, 
                                single_cutoff = FALSE, otu.name=otu.name, column_name_rot = 40, left_margin = 4)
heatmap_list$heatmap_plot


filename <- "./Figure/RM_figures/overlap heatmap liver negative top50perc.svg"
#png(filename, width = 15.75*330, height = 5.75*330, res = 330)
svg(filename, width = 15.75, height = 5.65)
heatmap_list$heatmap_plot
dev.off()


#"#CC79A7"
#"#0072B2"

#--------------------------------------
i="lymphoid leukemia"

sig_components <- Elix_components_w_sig_list_sub_cutoff[[i]]
#"Elixhauser_Arrhythmia"  "Elixhauser_FluidsLytes" "Elixhauser_NeuroOther" 
FDR_df <- merge_list_sig_Q_w_comorb[[i]][,c("Ctr.Q","Ex.Q",sig_components)]
R2_df <- merge_list_sig_R2_w_comorb[[i]][,c("Ctr.R2","Ex.R2",sig_components)]

grp.names <- colnames(FDR_df)
grp.labels <- colnames(FDR_df)
grp.labels[1:2] <- c("Cancer class vs. healthy", "Cancer class vs. other cancers")
sign_cutoffs <- c(0.1, 0.05, 0.01)

heatmap_list <- zicoseq_heatmap(FDR_df, R2_df, grp.names, grp.labels, col = col,sign_cutoffs, single_cutoff = FALSE, otu.name=otu.name, column_name_rot = 40)
heatmap_list$heatmap_plot

filename <- paste0("./Figure/RM_figures/overlap heatmap ",i,"_sig_Elix.png")
png(filename, width =17*330, height = 5.25*330, res = 330)
heatmap_list$heatmap_plot
dev.off()



colnames(FDR_df) <- grp.labels
sign_lists_to_plot[[i]] <- FDR_df
R2_lists_to_plot[[i]] <- R2_df


#remove Elix species; any rows with FDR <0.1 beyond the first 2 columns
to_keep_rows <- !rowSums(FDR_df[, 3:ncol(FDR_df),drop =F] < 0.1) > 0
length(which(to_keep_rows)) #24



#----------------------------------
i="multiple myeloma plasma cell"

sig_components <- Elix_components_w_sig_list_sub_cutoff[[i]] # no significance
# character(0)
FDR_df <- merge_list_sig_Q_w_comorb[[i]][,c("Ctr.Q","Ex.Q",sig_components)]
R2_df <- merge_list_sig_R2_w_comorb[[i]][,c("Ctr.R2","Ex.R2",sig_components)]

grp.names <- colnames(FDR_df)
grp.labels <- colnames(FDR_df)
grp.labels[1:2] <- c("Cancer class vs. healthy", "Cancer class vs. other cancers")
sign_cutoffs <- c(0.1, 0.05, 0.01)


heatmap_list <- zicoseq_heatmap(FDR_df, R2_df, grp.names, grp.labels, col = col,sign_cutoffs, single_cutoff = FALSE, 
                                otu.name=otu.name, column_name_rot = 40, left_margin = 5)
heatmap_list$heatmap_plot

filename <- paste0("./Figure/RM_figures/overlap heatmap ",i,"_sig_Elix.png")
png(filename, width =7.75*330, height = 4.5*330, res = 330)
heatmap_list$heatmap_plot
dev.off()


colnames(FDR_df) <- grp.labels
sign_lists_to_plot[[i]] <- FDR_df
R2_lists_to_plot[[i]] <- R2_df


#remove Elix species; any rows with FDR <0.1 beyond the first 2 columns
to_keep_rows <- !rowSums(FDR_df[, 3:ncol(FDR_df),drop =F] < 0.1) > 0
length(which(to_keep_rows)) #2

# 
# # remove Elix species; any rows with FDR <0.1 beyond the first 2 columns
# to_keep_rows <- !rowSums(FDR_df[, 3:ncol(FDR_df)] < 0.1) > 0
# length(which(to_keep_rows)) #32 # 32
# 
# heatmap_list <- zicoseq_heatmap(FDR_df[to_keep_rows,], R2_df[to_keep_rows, ], grp.names, grp.labels, col = col, sign_cutoffs, single_cutoff = FALSE,
#                                 otu.name=otu.name, column_name_rot = 40)
# 
# filename <- "./Figure/RM_figures/overlap heatmap Lymphoid_leukemia_sig_Elix_rm.png"
# png(filename, width =34*330, height = 6*330, res = 330)
# heatmap_list$heatmap_plot
# dev.off()



#--------------------------------------
i="neuroendocrine tumors"

sig_components <- Elix_components_w_sig_list_sub_cutoff[[i]]
FDR_df <- merge_list_sig_Q_w_comorb[[i]][,c("Ctr.Q","Ex.Q",sig_components)]
R2_df <- merge_list_sig_R2_w_comorb[[i]][,c("Ctr.R2","Ex.R2",sig_components)]

grp.names <- colnames(FDR_df)
grp.labels <- colnames(FDR_df)
grp.labels[1:2] <- c("Cancer class vs. healthy", "Cancer class vs. other cancers")
sign_cutoffs <- c(0.1, 0.05, 0.01)

heatmap_list <- zicoseq_heatmap(FDR_df, R2_df, grp.names, grp.labels, col = col,sign_cutoffs, single_cutoff = FALSE, otu.name=otu.name, column_name_rot = 40)
heatmap_list$heatmap_plot

filename <- paste0("./Figure/RM_figures/overlap heatmap ",i,"_sig_Elix.png")
png(filename, width = 75*330, height = 6*330, res = 330)
heatmap_list$heatmap_plot
dev.off()



colnames(FDR_df) <- grp.labels
sign_lists_to_plot[[i]] <- FDR_df
R2_lists_to_plot[[i]] <- R2_df



#remove Elix species; any rows with FDR <0.1 beyond the first 2 columns
to_keep_rows <- !rowSums(FDR_df[, 3:ncol(FDR_df),drop =F] < 0.1) > 0
length(which(to_keep_rows)) #273


FDR_df_pos <- FDR_df
R2_df_pos <- R2_df
FDR_df_pos <- FDR_df_pos[R2_df$Ctr.R2 > 0,]
R2_df_pos <- R2_df_pos[R2_df$Ctr.R2 > 0,]

#get the species that are specific and positive
names(rowSums(FDR_df_pos[,3:4] <= 0.1))[rowSums(FDR_df_pos[,3:4] <= 0.1) == 0]



#make a version of this with only positive and top 50% effect size
effect_size_cutoff <- 0.5
effect_cutoff <- min(sort(apply(abs(R2_df_pos), 1, max), decreasing=T)[1:nrow(R2_df_pos)*effect_size_cutoff])
effect_idx <- apply(abs(R2_df_pos), 1, function(x) sum(x > effect_cutoff) > 0)

FDR_df_sub <- as.matrix(FDR_df_pos[effect_idx,, drop =F])
R2_df_sub <- as.matrix(R2_df_pos[effect_idx,, drop =F])

dim(FDR_df_pos); dim(FDR_df_sub)

FDR_df_pos <- FDR_df_sub
R2_df_pos <- R2_df_sub


grp.labels <- c("Neuroendocrine tumors vs. healthy", "Neuroendocrine tumors vs. other cancers", 
                   "Elixhauser_Arrhythmia", "Elixhauser_Liver")

heatmap_list <- zicoseq_heatmap(FDR_df_pos, R2_df_pos, grp.names, grp.labels, col = col,sign_cutoffs, 
                                single_cutoff = FALSE, otu.name=otu.name, column_name_rot = 40, left_margin = 4)
heatmap_list$heatmap_plot


filename <- "./Figure/RM_figures/overlap heatmap neuroendocrine_tumors_sig_Elix_cutoff_positive 50perc.svg"
#filename <- "./Figure/RM_figures/overlap heatmap neuroendocrine_tumors_sig_Elix_cutoff_positive 50perc.png"
#png(filename, width = 15*330, height = 4.75*330, res = 330)
svg(filename, width = 15, height = 4.75)
heatmap_list$heatmap_plot
dev.off()



#negative
FDR_df_neg <- FDR_df
R2_df_neg <- R2_df
FDR_df_neg <- FDR_df_neg[R2_df$Ctr.R2 < 0,]
R2_df_neg <- R2_df_neg[R2_df$Ctr.R2 < 0,]

effect_size_cutoff <- 0.2
effect_cutoff <- min(sort(apply(abs(R2_df_neg), 1, max), decreasing=T)[1:nrow(R2_df_neg)*effect_size_cutoff])
effect_idx <- apply(abs(R2_df_neg), 1, function(x) sum(x > effect_cutoff) > 0)

FDR_df_sub <- as.matrix(FDR_df_neg[effect_idx,, drop =F])
R2_df_sub <- as.matrix(R2_df_neg[effect_idx,, drop =F])

dim(FDR_df_neg); dim(FDR_df_sub)

FDR_df_neg <- FDR_df_sub
R2_df_neg <- R2_df_sub


heatmap_list <- zicoseq_heatmap(FDR_df_neg, R2_df_neg, grp.names, grp.labels, col = col,sign_cutoffs, 
                                single_cutoff = FALSE, otu.name=otu.name, column_name_rot = 40, left_margin = 6)
heatmap_list$heatmap_plot


filename <- "./Figure/RM_figures/overlap heatmap neuroendocrine_tumors_sig_Elix_cutoff_20perc_negative.svg"
#filename <- "./Figure/RM_figures/overlap heatmap neuroendocrine_tumors_sig_Elix_cutoff_20perc_negative.png"
#png(filename, width = 19*330, height = 4.75*330, res = 330)
svg(filename, width = 19, height = 4.75)
heatmap_list$heatmap_plot
dev.off()




#remove Elix species; any rows with FDR <0.1 beyond the first 2 columns
#to_keep_rows <- !rowSums(FDR_df[, 3:ncol(FDR_df)] < 0.1) > 0
#length(which(to_keep_rows)) #193 # 277

#heatmap_list <- zicoseq_heatmap(FDR_df[to_keep_rows,], R2_df[to_keep_rows, ], grp.names, grp.labels, col = col,sign_cutoffs, single_cutoff = FALSE,
#                                otu.name=otu.name, column_name_rot = 40)

#filename <- paste0("./Figure/RM_figures/overlap heatmap ",i,"_sig_Elix_rm.png")
#png(filename, width = 62*330, height = 6.5*330, res = 330)
#heatmap_list$heatmap_plot
#dev.off()


#-------------------------------------------------------------------------------
#write to file

lapply(sign_lists_to_plot, dim)
lapply(R2_lists_to_plot, dim)

sign_lists_to_plot_copy <- sign_lists_to_plot
R2_lists_to_plot_copy <- R2_lists_to_plot

names(sign_lists_to_plot_copy) <- paste0("s_q_", names(sign_lists_to_plot_copy))
names(R2_lists_to_plot_copy) <- paste0("s_R2_", names(R2_lists_to_plot_copy))

#combine the files
lists_to_write_to_supplement <- c(sign_lists_to_plot_copy[1], R2_lists_to_plot_copy[1],
                                  sign_lists_to_plot_copy[2], R2_lists_to_plot_copy[2],
                                  sign_lists_to_plot_copy[3], R2_lists_to_plot_copy[3],
                                  sign_lists_to_plot_copy[4], R2_lists_to_plot_copy[4],
                                  sign_lists_to_plot_copy[5], R2_lists_to_plot_copy[5])

names(lists_to_write_to_supplement) <- gsub("intrahepatic", "intrah", names(lists_to_write_to_supplement))
names(lists_to_write_to_supplement) <- gsub("cell", "c", names(lists_to_write_to_supplement))

nchar(names(lists_to_write_to_supplement))

#make column names consistent
for (i in 1:length(lists_to_write_to_supplement)) {
  colnames(lists_to_write_to_supplement[[i]])[1:2] <- c("Cancer class vs. healthy",	"Cancer class vs. other cancers")
}
#lapply(lists_to_write_to_supplement, function(x) colnames(x))


filename <- "Code/Submission/Supplementary tables/species_for_Supplementary Table 4.xlsx"

list_of_datasets <- lists_to_write_to_supplement
write.xlsx(list_of_datasets, file = filename, append=T, rowNames=T)



#-------------------------------------------------------------------------------
#plot species abundances stratified by comorbidity

#boxplot; healthy control, all cancer, all cancer + comorbid, cancer class, cancer class + comorbid

clin_meta$icd10_first_3_name[is.na(clin_meta$icd10_first_3_name)] <- "healthy"
head(sort(table(clin_meta$icd10_first_3_name), decreasing = T), 10) #testing

#rff relative abundance, sweep
tax_table_RA <- sweep(tax_table, MARGIN = 2, colSums(tax_table), '/') * 100

#update cancer_row_list in same order as sign_lists_to_plot
names(sign_lists_to_plot)

group_vec <- clin_meta$icd10_first_3_name
cancer_row_list <- list(
  # grep("bronchus", group_vec),
  grep("esophagus", group_vec),
  grep("liver and intrahepatic bile ducts", group_vec),
  grep("Lymphoid leukemia", group_vec), 
  grep("plasma cell", group_vec),
  grep("neuroendocrine", group_vec)
)


elix_row_names_list <- lapply(sign_lists_to_plot, function(x) colnames(x)[3:ncol(x)])
healthy_rows <- which(group_vec == "healthy")

sel_cols <- brewer.pal(5, "Dark2")

for (i in 1:length(sign_lists_to_plot)) { #cancer class
  sig_df <- sign_lists_to_plot[[i]]
  elix_row_names <- elix_row_names_list[[i]]
  cancer_group_rows <- cancer_row_list[[i]]
  all_other_cancer_rows <- which(!1:ncol(tax_table_RA) %in% c(healthy_rows, cancer_group_rows))
  cancer_class <- names(sign_lists_to_plot)[i]
  
  pdf(paste0("./Figure/RM_figures/Comorbidity species_", cancer_class, ".pdf"), width=5, height=4.5)
  for (j in 1:length(elix_row_names)) { #elix component
    sig_species <- rownames(sig_df)[which(sig_df[elix_row_names[j]] < 0.1)]
    sig_species <- gsub(" ", "_", sig_species)
    
    for (k in 1:length(sig_species)) {
      
      temp_ra <- sqrt(as.numeric(tax_table_RA[sig_species[k],,drop =F]))
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
      mtext("sqrt(relative abundance %)", 2, line = 2.5)
      mtext(paste0(cancer_class, " _ ", elix_row_names[j], "\n", sig_species[k]), 3, line = 0, cex=1)
      
    }
  }
  dev.off()
}


#-------------------------------------
#-------------------------------------
#split version 

#orange       skyblue   bluishgreen        yellow          blue    vermillion reddishpurple 
#"#E69F00"     "#56B4E9"     "#009E73"     "#F0E442"     "#0072B2"     "#D55E00"     "#CC79A7" 


sel_cols <- c("#0072B2", "#CC79A7", "#E69F00")

for (i in 1:length(sign_lists_to_plot)) { #cancer class
  sig_df <- sign_lists_to_plot[[i]]
  elix_row_names <- elix_row_names_list[[i]]
  cancer_group_rows <- cancer_row_list[[i]]
  all_other_cancer_rows <- which(!1:ncol(tax_table_RA) %in% c(healthy_rows, cancer_group_rows))
  cancer_class <- names(sign_lists_to_plot)[i]
  
  pdf(paste0("./Figure/RM_figures/Comorbidity species split_", cancer_class, ".pdf"), width=5, height=4)
  for (j in 1:length(elix_row_names)) { #elix component
    sig_species <- rownames(sig_df)[which(sig_df[elix_row_names[j]] < 0.1)]
    sig_species <- gsub(" ", "_", sig_species)
    
    for (k in 1:length(sig_species)) {
      
      temp_ra <- sqrt(as.numeric(tax_table_RA[sig_species[k],,drop =F]))
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
      
      par(mar=c(5,4,2,1))
      par(lheight = 0.8)
      
      boxplot(values ~ group, data, las=1, col = sel_cols, varwidth = F, xaxt = "n", outline = F, xlab="", 
              at=c(1,2,4,5), ylim=c(0, ymax), ylab="", xlim=c(0.5, 5.5))
      #add insert boxplots
      par(new = TRUE)
      boxplot(values ~ group, data_all_other, las=1, col = sel_cols[2], varwidth = F, xaxt = "n", outline = F, xlab="", 
              at=c(1,2,2.75,3.25,4,5), boxwex = 0.45, ylim=c(0, ymax), yaxt="n", ylab="", xlim=c(0.5, 5.5))
      par(new = TRUE)
      boxplot(values ~ group, data_cc_other, las=1, col = sel_cols[3], varwidth = F, xaxt = "n", outline = F, xlab="", 
              at=c(1,2,3,4,4.75,5.25), boxwex = 0.45, ylim=c(0, ymax), yaxt="n", ylab="", xlim=c(0.5, 5.5))
      
      axis(1, at=c(1,2,2.75,3.25,4,4.75,5.25), labels = FALSE)
      
      text(x = c(1,2,2.75,3.25,4,4.75,5.25), 
           y = par("usr")[3]*3, 
           labels = c("healthy", "other cancer",
                      "other cancer\nno comorbidity",
                      "other cancer\ncomorbidity", 
                      "cancer class",
                      "cancer class\nno comorbidity",
                      "cancer class\ncomorbidity"),
           srt = 45, adj = 1, xpd = TRUE, cex = 0.8)
      
      mtext("sqrt(relative abundance %)", 2, line = 2.5)
      mtext(paste0(cancer_class, " _ ", elix_row_names[j], "\n", sig_species[k]), 3, line = 0, cex=1)
      
    }
  }
  dev.off()
}



#-------------------------------------------------------------------------------
#make plots for all interesting species as well without factoring in comorbidities

cancer_row_list <- list(
  grep("esophagus", group_vec),
  grep("liver and intrahepatic bile ducts", group_vec),
  grep("Lymphoid leukemia", group_vec), 
  grep("plasma cell", group_vec),
  grep("neuroendocrine", group_vec)
)

names(merge_list_sig_Q_w_comorb)

for (i in 1:length(merge_list_sig_Q_w_comorb)) { #cancer class
  sig_df <- merge_list_sig_Q_w_comorb[[i]]
  
  cancer_group_rows <- cancer_row_list[[i]]
  all_other_cancer_rows <- which(!1:nrow(tax_table_RA) %in% c(healthy_rows, cancer_group_rows))
  cancer_class <- names(merge_list_sig_Q_w_comorb)[i]
  
  sig_species <- rownames(sig_df)
  sig_species <- gsub(" ", "_", sig_species)
  
  pdf(paste0("./Figure/RM_figures/Significant species_scatter_", cancer_class, ".pdf"), width=2.5, height=4)
  for (k in 1:length(sig_species)) {
    
    temp_ra <- sqrt(as.numeric(tax_table_RA[sig_species[k],]))
    plot_df <- as.data.frame(cbind(c(rep("healthy", times=length(healthy_rows)), rep("other cancers", times=length(all_other_cancer_rows)), rep("cancer class", times=length(cancer_group_rows))),
          c(temp_ra[healthy_rows], temp_ra[all_other_cancer_rows], temp_ra[cancer_group_rows])))
    names(plot_df) <- c("group", "ra")
    plot_df$ra <- as.numeric(plot_df$ra)
    
    
    p_plot <- ggboxplot(plot_df, "group", "ra",
                        color = "group", palette = sel_cols, width = 0.8,
                        add = "jitter", ylab = "sqrt(relative_abundance %)") +
      scale_y_continuous(trans = 'sqrt') +
      ggtitle(gsub("s__|_", " ", sig_species[k])) +
      theme(
        legend.position = "none",
        axis.text = element_text(family = "sans", size = 10),
        axis.title = element_text(family = "sans", size = 10),
        plot.title = element_text(face = "italic", size = 10)
      ) +
      xlab("") +
      scale_x_discrete(labels = c("healthy", "other\ncancers", "cancer\nclass"))
    
    print(p_plot)
  }
  dev.off()
}



#-------------------------------------------------------------------------------
#For all cancer classes and a single significance level <0.1. Show number of species vs control, vs other cancers, and overlap.
#And how many after filtering the cancer-specific comorbidities


cutoff <- 0.1

merge_list <- vector("list", length = length(canX_Ctr_list))
names(merge_list) <- names(canX_Ctr_list)

sig_control_list <- vector("list", length = length(canX_Ctr_list))
sig_cancer_list <- vector("list", length = length(canX_Ctr_list))

#merge these in a for loop
for (i in 1:length(canX_Ctr_list)) {
  
  control <- canX_Ctr_list[[i]]
  cancer <- canX_Ex_list[[i]]
  
  idx.c <- control$Species[which(control$Q<cutoff)]
  idx.d <- cancer$Species[which(cancer$Q<cutoff)]
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
n_wo_sig_elix <- rep(0, times= length(merge_list))
names(n_wo_sig_elix) <- names(merge_list)


n_wo_sig_elix["esophagus"] <- 6
n_wo_sig_elix["liver intrahepatic bile ducts"] <- 36
n_wo_sig_elix["lymphoid leukemia"] <- 24
n_wo_sig_elix["multiple myeloma plasma cell"] <- 2
n_wo_sig_elix["neuroendocrine tumors"] <- 273


plot_df <- cbind(names(merge_list), 
      as.character(lapply(sig_control_list, length)), 
      as.character(lapply(sig_cancer_list, length)),
      as.character(lapply(merge_list, nrow)),
      as.numeric(n_wo_sig_elix)
      )


plot_df[plot_df == "NULL"] <- 0
plot_df

colnames(plot_df) <- c("cancer class", "vs. healthy", "vs. other cancer", "overlap", "overlap w/o comorbidities")
plot_df <- as.data.frame(plot_df)

plot_df$`vs. healthy` <- as.numeric(plot_df$`vs. healthy`)
plot_df$`vs. other cancer` <- as.numeric(plot_df$`vs. other cancer`)
plot_df$overlap <- as.numeric(plot_df$overlap)
plot_df$`overlap w/o comorbidities` <- as.numeric(plot_df$`overlap w/o comorbidities`)

plot_df <- plot_df[order(plot_df$`overlap w/o comorbidities`, decreasing = T),]

gt(plot_df) %>%
  tab_header(title = "Number of species when correcting for comorbidities") %>%
  opt_vertical_padding(scale = 0.1)



#-------------------------------------------------------------------------------
#get numbers for the main paper; check how many of the consistently cancer associated species are also associated to a comorbidity


cutoff <- 0.1

merge_list <- vector("list", length = length(canX_Ctr_list))
names(merge_list) <- names(canX_Ctr_list)

sig_control_list_2 <- vector("list", length = length(canX_Ctr_list))
sig_control_list_pos <- vector("list", length = length(canX_Ctr_list))

#merge these in a for loop
for (i in 1:length(canX_Ctr_list)) {
  
  control <- canX_Ctr_list[[i]]
  idx.c <- control$Species[which(control$Q < cutoff)]
  control_sub <- control[idx.c, c("P", "R2", "Q")]
  sig_control_list_2[[i]] <- idx.c
  sig_control_list_pos[[i]] <- rownames(control_sub)[control_sub$R2 > 0]
}

lapply(sig_control_list_2, length)
median(as.numeric(lapply(sig_control_list_2, length)))
# 332

median(as.numeric(lapply(sig_control_list_pos, length)))
# 91.5


length(unique(unlist(sig_control_list_2)))
#1043

length(unique(unlist(sig_control_list_pos)))
#323


#22 cancer classes now
12 /length(sig_control_list_2) #54.5%
length(which(table(unlist(sig_control_list_2)) >= 12))
#307 species

length(which(table(unlist(sig_control_list_pos)) >= 12))
#67


cutoff <- 20
cutoff /length(sig_control_list_2) #90.9%
length(which(table(unlist(sig_control_list_2)) >= cutoff))
#44; make heatmap of these species
names(which(table(unlist(sig_control_list_2)) >= cutoff))
#[1] "s__Alistipes communis"             "s__Alistipes finegoldii"           "s__Alistipes putredinis"           "s__Alistipes senegalensis"        
#[5] "s__Bifidobacterium dentium"        "s__Blautia_A sp900551465"          "s__CAG-103 sp900543625"            "s__CAG-110 sp003525905"           
#[9] "s__CAG-127 sp900319515"            "s__CAG-170 sp000436735"            "s__CAG-170 sp900556635"            "s__CAG-45 sp000438375"            
#[13] "s__CAG-83 sp003487665"             "s__Clostridium_Q sp003024715"      "s__Coprococcus eutactus_A"         "s__Dysosmobacter sp001916835"     
#[17] "s__Enterocloster clostridioformis" "s__ER4 sp000765235"                "s__Erysipelatoclostridium ramosum" "s__Eubacterium_F sp000433735"     
#[21] "s__Faecalibacterium prausnitzii_C" "s__Faecalibacterium prausnitzii_I" "s__Faecalibacterium sp900758465"   "s__GCA-900066135 sp900543575"     
#[25] "s__HGM13006 sp900757695"           "s__Klebsiella pneumoniae"          "s__PeH17 sp000435055"              "s__Ruminococcus_C sp000433635"    
#[29] "s__SFFH01 sp900542395"             "s__SFFH01 sp900542445"             "s__SFFH01 sp900548125"             "s__Streptococcus parasanguinis"   
#[33] "s__Streptococcus parasanguinis_A"  "s__Streptococcus parasanguinis_B"  "s__Streptococcus parasanguinis_C"  "s__Streptococcus parasanguinis_D" 
#[37] "s__Streptococcus sp000448565"      "s__Streptococcus sp001556435"      "s__Streptococcus sp900766505"      "s__Streptococcus vestibularis"    
#[41] "s__Tidjanibacter inops_A"          "s__UBA11524 sp000437595"           "s__UBA1777 sp900547315"            "s__UMGS403 sp900540275"   


length(which(table(unlist(sig_control_list_pos)) >= cutoff))
#15
names(which(table(unlist(sig_control_list_pos)) >= cutoff))
#[1] "s__Bifidobacterium dentium"        "s__Blautia_A sp900551465"          "s__Enterocloster clostridioformis"
#[4] "s__Erysipelatoclostridium ramosum" "s__Klebsiella pneumoniae"          "s__Streptococcus parasanguinis"   
#[7] "s__Streptococcus parasanguinis_A"  "s__Streptococcus parasanguinis_B"  "s__Streptococcus parasanguinis_C" 
#[10] "s__Streptococcus parasanguinis_D"  "s__Streptococcus sp000448565"      "s__Streptococcus sp001556435"     
#[13] "s__Streptococcus sp900766505"      "s__Streptococcus vestibularis"     "s__UMGS403 sp900540275" 



#all comparisons
length(which(table(unlist(sig_control_list_2)) >= 22))
#7
names(which(table(unlist(sig_control_list_2)) >= 22))
#[1] "s__Bifidobacterium dentium"       "s__CAG-103 sp900543625"           "s__Klebsiella pneumoniae"         "s__Streptococcus parasanguinis_A"
#[5] "s__Streptococcus parasanguinis_D" "s__Streptococcus sp900766505"     "s__Streptococcus vestibularis"  

consistent_species_55perc <- names(which(table(unlist(sig_control_list_2)) >= 12))
consistent_species_91perc <- names(which(table(unlist(sig_control_list_2)) >= 20))


consistent_pos_species_55perc <- names(which(table(unlist(sig_control_list_pos)) >= 12))
consistent_pos_species_91perc <- names(which(table(unlist(sig_control_list_pos)) >= 20))



#-------------------------------------------------------------------------------
#see which of these species are also associated to Elixhauser general and its components

#general elix is: elix_diff_Q, elix_diff_R2, 
#components is: Q_sign_df_comorbidities, R2_sign_df_comorbidities
 
table(consistent_species_55perc %in% names(elix_diff_Q))
#FALSE  TRUE 
#  13   294 
table(consistent_species_91perc %in% names(elix_diff_Q))
#FALSE TRUE 
# 1     43

consistent_species_55perc_sub <- consistent_species_55perc[consistent_species_55perc %in% names(elix_diff_Q)]
#all pos species are in here

#29 to general elix score
table(elix_diff_Q[consistent_species_55perc_sub] <= 0.1 & elix_diff_R2[consistent_species_55perc_sub] > 0)
#FALSE  TRUE 
#265    29
#get the names
sign_elix_general_55_perc_names <- names(which(elix_diff_Q[consistent_species_55perc_sub] <= 0.1 & elix_diff_R2[consistent_species_55perc_sub] > 0))

#pos only
table(elix_diff_Q[consistent_pos_species_55perc] <= 0.1 & elix_diff_R2[consistent_pos_species_55perc] > 0)
#FALSE  TRUE 
# 39    28 
#28 / 67 = 41.8%


#species
table(elix_diff_Q[consistent_species_91perc] <= 0.1 & elix_diff_R2[consistent_species_91perc] > 0)
#FALSE  TRUE 
# 33     10 
sign_elix_general_91_perc_names <- names(which(elix_diff_Q[consistent_species_91perc] <= 0.1 & elix_diff_R2[consistent_species_91perc] > 0))

#pos only
table(elix_diff_Q[consistent_pos_species_91perc] <= 0.1 & elix_diff_R2[consistent_pos_species_91perc] > 0)
#FALSE  TRUE 
# 5    10 



#
number_of_sign_elix <- rowSums((Q_sign_df_comorbidities[consistent_species_55perc_sub,] <= 0.1 & R2_sign_df_comorbidities[consistent_species_55perc_sub,] > 0) == TRUE)
Q_sign_df_comorbidities["s__Streptococcus sp001556435",] # has 3, Elixhauser_Anemia, Elixhauser_Liver, Elixhauser_PUD

table(number_of_sign_elix > 0)
#FALSE  TRUE 
#245    49
sign_elix_components_55_perc_names <- names(which(number_of_sign_elix > 0))

#pos only
number_of_sign_elix_pos <- rowSums((Q_sign_df_comorbidities[consistent_pos_species_55perc,] <= 0.1 & R2_sign_df_comorbidities[consistent_pos_species_55perc,] > 0) == TRUE)
table(number_of_sign_elix_pos > 0)
#FALSE  TRUE 
# 23    44
#44 / 67 = 65.7%


consistent_species_91perc.sub <- consistent_species_91perc[(consistent_species_91perc %in% rownames(Q_sign_df_comorbidities))]
number_of_sign_elix <- rowSums((Q_sign_df_comorbidities[consistent_species_91perc.sub,] <= 0.1 & R2_sign_df_comorbidities[consistent_species_91perc.sub,] > 0) == TRUE)
# number_of_sign_elix <- rowSums((Q_sign_df_comorbidities[consistent_species_91perc,] <= 0.1 & R2_sign_df_comorbidities[consistent_species_91perc,] > 0) == TRUE)
table(number_of_sign_elix > 0)
#FALSE  TRUE 
#31     12
sign_elix_components_91_perc_names <- names(which(number_of_sign_elix > 0))

#pos only
number_of_sign_elix_pos <- rowSums((Q_sign_df_comorbidities[consistent_pos_species_91perc,] <= 0.1 & R2_sign_df_comorbidities[consistent_pos_species_91perc,] > 0) == TRUE)
table(number_of_sign_elix_pos > 0)
#FALSE  TRUE 
#3    12 
#12 / 15 = 80%


#-----------------------------------
#look at the numbers
length(consistent_species_55perc); length(sign_elix_general_55_perc_names); length(sign_elix_components_55_perc_names); length(unique(c(sign_elix_general_55_perc_names, sign_elix_components_55_perc_names)))
#307
#29
#49
#51
51/307 #= 16.6%



length(consistent_species_91perc); length(sign_elix_general_91_perc_names); length(sign_elix_components_91_perc_names); length(unique(c(sign_elix_general_91_perc_names, sign_elix_components_91_perc_names)))
#44
#10
#12
#13
13 / 44 #29.5%


#-----------------------------------
#20 out of 22 as cutoff

data_list <- list()
for (i in 1:length(canX_Ctr_list)) {
  control <- canX_Ctr_list[[i]]
  control_sub <- control[names(which(table(unlist(sig_control_list_2)) >= 20)),c("P", "R2", "Q")]
  data_list[[i]] <- control_sub
}

R2_df <- do.call(cbind, lapply(data_list, function(x) x$R2))
Q_df <- do.call(cbind, lapply(data_list, function(x) x$Q))

colnames(R2_df) <- names(canX_Ctr_list)
colnames(Q_df) <- names(canX_Ctr_list)


rownames(R2_df) <- names(which(table(unlist(sig_control_list_2)) >= 20))
rownames(Q_df) <- names(which(table(unlist(sig_control_list_2)) >= 20))

R2_df[is.na(R2_df)] <- 0
Q_df[is.na(Q_df)] <- 1

dim(R2_df)



#------------------------------------
#make heatmap out of this including the family level top.

dim(R2_df); dim(Q_df)
colnames(R2_df); rownames(R2_df)


Q.sig <- Q_df
R2.sig <- R2_df

unique.sig <- rownames(Q.sig)

sign_cutoffs <- c(0.1,0.05,0.01)
cutoff <- 0.1; q.cut1 = 0.1; q.cut2=0.05;q.cut3=0.01

rownames(Q.sig) <- gsub('.*;g__|s__','',rownames(Q.sig))
rownames(R2.sig) <- gsub('.*;g__|s__','',rownames(R2.sig))

tmp <- Q.sig
tmp[] <- apply(tmp, 2, function(x) ifelse(x<0.05,'A','B'))
tmp <- as.data.frame(tmp) %>% arrange(across(everything()))

Q.sig <- t(Q.sig[rownames(tmp),])
R2.sig <- t(R2.sig[rownames(tmp), ])

level <- 'Family'
otu.name <- as.data.frame(data.obj$otu.name)

label <- otu.name[otu.name$Species %in% gsub('.*;g__|s__','s__',unique.sig),] %>% dplyr::select(c(level, 'Species')) 

label$Species <- gsub('o__|f__|c__|p__|s__|g__','',label$Species)
label <- label[match(rownames(tmp), label$Species),]
label[,level] <- gsub('o__|f__|c__|p__|g__','',label[,level])

# col <- c(brewer.pal(12,'Set3'),brewer.pal(12,'Paired'),brewer.pal(9,'Set1'),brewer.pal(8,'Set2'))[1:length(unique(label[,level]))]
# names(col) <- unique(label[,level])

col <- c(
  "Streptococcaceae" = "#1f77b4",           # blue
  "Ruminococcaceae" = "#ff7f0e",            # orange
  "Oscillospiraceae" = "#2ca02c",           # green
  "Enterobacteriaceae" = "#d62728",         # red
  "Rikenellaceae" = "#9467bd",              # purple
  "Erysipelatoclostridiaceae" = "#8c564b",  # brown
  "Lachnospiraceae" = "#e377c2",            # pink
  "CAG-138" = "#7f7f7f",                    # gray
  "Bifidobacteriaceae" = "#bcbd22",         # yellow-green
  "Acutalibacteraceae" = "#17becf",         # cyan
  "CAG-74" = "#aec7e8",                     # light blue
  "Coriobacteriaceae" = "forestgreen",      # deep green
  "Anaerovoracaceae" = "#f781bf",           # bright pink
  "Erysipelotrichaceae" = "#999999",        # mid gray
  "Enterococcaceae" = "#8dd3c7",            # pale cyan
  "Lactobacillaceae" = "#ffffb3",           # pale yellow
  "Eggerthellaceae" = "#fb8072",            # salmon
  "Peptostreptococcaceae" = "#cab2d6",      # soft purple
  "Methanobacteriaceae" = "#66c2a5",        # teal
  "Succinivibrionaceae" = "#fc8d62",        # soft orange
  "Burkholderiaceae" = "#8da0cb",           # lavender
  "Bacteroidaceae" = "#e78ac3",             # pink
  "Muribaculaceae" = "#a6d854",             # lime green
  "Barnesiellaceae" = "#ffd92f",            # yellow
  "Coprobacteraceae" = "#e5c494",           # beige
  "Tannerellaceae" = "#b3b3b3",             # silver gray
  "Marinifilaceae" = "#1b9e77",             # deep green
  "Butyricicoccaceae" = "#d95f02",          # reddish orange
  "CAG-382" = "#7570b3",                    # muted purple
  "UBA1390" = "#e7298a",                    # magenta
  "Anaerotignaceae" = "#66a61e",            # olive green
  "UBA1750" = "#e6ab02",                    # mustard
  "QAND01" = "#a6761d",                     # brownish
  "Veillonellaceae" = "#666666",            # dark gray
  "Dialisteraceae" = "#dede00",             # bright yellow
  "Acidaminococcaceae" = "#377eb8",         # blue
  "Actinomycetaceae" = "#984ea3",           # deep violet
  "QAMH01" = "#4daf4a"                      # green
)

ha_combined <- HeatmapAnnotation(
  ` ` = label[, level],
  annotation_height = unit(1, "cm"),
  border = TRUE,
  col = list(` ` = col),
  show_legend = F
)



#scale_max <- max(c(abs(min(R2.sig)), max(R2.sig)))
rownames(R2.sig) <- tolower(gsub("_", " ", rownames(R2.sig)))
rownames(R2.sig) <- tolower(gsub("malignant neoplasm of |malignant |malignant |other and | neoplasm| neoplasm of| neoplasms|neoplasm of |\\,|and ", "", rownames(R2.sig)))

cancer_class_names <- rownames(R2.sig)


lgd_sig_p = Legend(pch = c("*","**","***"), type = "points", labels = paste0("<", sign_cutoffs), background="white", title = "q")


heatmap_plot <- Heatmap(R2.sig, name = "R2 x effect direction", 
        heatmap_legend_param = list(direction = "vertical", title="R2"),
        col = colorRamp2(c(min(R2.sig), 0, max(R2.sig)), c("#0072B2",'white', "#CC79A7")),
        column_gap = unit(1, "mm"), 
        top_annotation = ha_combined,
        column_split = (label[,level]),
        column_title_rot = 30,
        column_names_rot = 40,
        column_names_max_height = unit(10, "cm"),
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
        row_names_max_width = 1.2*max_text_width(rownames(R2.sig)),
        row_names_gp = gpar(fontsize = 14)
) 

left_margin = 4

packed_legends <- packLegend(lgd_sig_p, direction = "vertical")
heatmap_plot <- draw(heatmap_plot, annotation_legend_list = list(packed_legends), merge_legend=T, 
                     padding = unit(c(0.1, left_margin, 0.1, 0.1), "cm")) 


png(paste0('./Figure/RM_figures/heatmap_consistent_species.png'), width =18*330, height = 8.5*330, res=300)
heatmap_plot
dev.off()





