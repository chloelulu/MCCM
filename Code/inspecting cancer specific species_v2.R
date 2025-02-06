require(GUniFrac)
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

#tax_table[1:2,1:2]

#make function for heatmap from dataframe
source("~/Dropbox/Mayo_RS/R/general functions/zicoseq_heatmap_from_df.R")

#order the clin meta based on tax_table
clin_meta <- clin_meta[colnames(tax_table),]


#-------------------------------------------------------------------------------

load("Figure/subCancerX_Control/DAA_P_R2.RData")
zico_res_ctr <- list(P.All$Species, Q.All$Species, R2.All$Species)
names(zico_res_ctr) <- c("P", "Q", "R2")


#do.call(cbind, zico_res_ctr)

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

elix_variables_filenames <- list.files(paste0(temp_dir, elix_variables))
elix_variables_filenames_full <- paste0(temp_dir, elix_variables, "/", elix_variables_filenames)

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
head(elix_variables_diff_list$Elixhauser_Alcohol$R2.list$Species)
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

#lapply(elix_variables_diff_list_R2, head)
#lapply(elix_variables_diff_list_R2_sign, head)
#lapply(elix_variables_diff_list_Q, head)


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
#[1] "s__Bifidobacterium breve"           "s__Pauljensenia bouchesdurhonensis" "s__Ruminococcus_D sp900752625"      "s__ER4 sp900546295"                
#[5] "s__Enterocloster sp005845215"       "s__Catenibacterium sp900540665"     "s__Citrobacter_A rodentium"         "s__Eisenbergiella tayi"            
#[9] "s__Veillonella rogosae"             "s__Clostridium_Q sp900547735"       "s__Enterocloster sp900540675"       "s__Parabacteroides distasonis"     
#[13] "s__Escherichia sp002965065"         "s__Klebsiella_A michiganensis"      "s__Rothia aeria"                    "s__TF01-11 sp001414325"            
#[17] "s__Butyricimonas faecalis"          "s__Agathobacter sp900550845"        "s__Agathobacter sp900557055"        "s__Agathobaculum sp900544475"      
#[21] "s__Bacteroides stercoris"           "s__Clostridium_AP sp000509125"      "s__Escherichia coli"                "s__Escherichia coli_E"             
#[25] "s__Eubacterium_G sp900550135"       "s__Gemmiger sp900539695"            "s__Bacteroides fragilis_A"          "s__Collinsella sp900541665"        
#[29] "s__Eubacterium limosum"             "s__Blautia coccoides_A"             "s__Bifidobacterium infantis"        "s__Erysipelatoclostridium ramosum" 
#[33] "s__Extibacter hylemonae"            "s__Catenibacterium mitsuokai"       "s__Enterocloster citroniae"         "s__Lacticaseibacillus paracasei"   
#[37] "s__Veillonella parvula_A"           "s__1XD42-69 sp014287635"            "s__Blautia_A massiliensis"          "s__BX12 sp902363595"  


#-------------------------------------------------------------------------------
#also get the general Elixhauser associations

load("./Result/CancerOnly/DAA/Elix_score/Elix_score_ZicoSeq.Rdata")
Elix_res_list <- diff.obj

names(Elix_res_list$R2.list)
head(Elix_res_list$R2.list$Species,3)
#                           Func1
#s__14-2 sp000403255     0.0005747840
#s__1XD42-69 sp014287635 0.0017011653
#s__Absicoccus porci     0.0001284017

elix_diff_Q <- Elix_res_list$qv.list$Species[,'Qvalue']
elix_diff_R2 <- Elix_res_list$R2.list$Species[,'Func1']
elix_diff_coef <- Elix_res_list$coef.list$Species[,"Elix_score"]

elix_diff_R2 <- elix_diff_R2 * sign(elix_diff_coef)
length(elix_diff_R2[elix_diff_Q < 0.05])
#135
table(elix_diff_R2[elix_diff_Q < 0.05] > 0)
#FALSE (negatively associated with ECI)  TRUE (positively associated with ECI)
#                 56                                      79 

#-------------------------------------
#get number of species FDR < 0.05 with positive and negative effect size for components

table_list <- list()
for (i in 1:ncol(R2_sign_df_comorbidities)) {
  R2_vec_temp <- R2_sign_df_comorbidities[,i]
  table_list[[i]] <- table(R2_vec_temp[Q_sign_df_comorbidities[,i] < 0.05] > 0)
}
names(table_list) <- colnames(R2_sign_df_comorbidities)

#TRUE means positively associated with the ECI component

#Elixhauser_Anemia
#TRUE 
#12

#Elixhauser_Arrhythmia
#FALSE  TRUE 
#6     1 

#Elixhauser_DMcx
#TRUE 
#2 

#Elixhauser_FluidsLytes
#FALSE  TRUE 
#3    17 

#Elixhauser_Liver
#TRUE 
#18 

#Elixhauser_NeuroOther
#TRUE 
#28 

#Elixhauser_PUD
#TRUE 
#90 

#Elixhauser_Renal
#TRUE 
#21 

#Elixhauser_WeightLoss
#FALSE  TRUE 
#8    28 


#sort(R2_sign_df_comorbidities[,"Elixhauser_WeightLoss"][Q_sign_df_comorbidities[,"Elixhauser_WeightLoss"] < 0.05])



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
length(canX_Ctr_list); length(canX_Ex_list) #23 cancer classes
#cbind(names(canX_Ctr_list), names(canX_Ex_list))


#-------------------------------------------------------------------------------
#get numbers only for the cancer class comparisons without considering the vs healthy comparisons
#with positive and negative

canX_Ex_list

res_list_canX <- lapply(canX_Ex_list, function(x) table(x$Q <= 0.2))
#res_list_canX

#1186 species considered

#$neuroendocrine_tumors
#FALSE  TRUE 
#677   509 
#$Lymphoid_leukemia
#FALSE  TRUE 
#986   200 
#$prostate
#FALSE  TRUE 
#1098    88 
#$liver_intrahepatic_bile_ducts
#FALSE  TRUE 
#1113    73 
#$esophagus
#FALSE  TRUE 
#1177     9 
#$bronchus_lung
#FALSE  TRUE 
#1180     6 
#$pancreas
#FALSE  TRUE 
#1182     4 
#$plasma_cell_neoplasms
#FALSE  TRUE 
#1183     3 
#$ovary
#FALSE  TRUE 
#1185     1 

res_list_canX_direction <- lapply(canX_Ex_list, function(x) table(x$R2[which(x$Q <= 0.2)] > 0))
#neuroendocrine_tumors
#FALSE  TRUE 
#301   208 
#Lymphoid_leukemia
#FALSE  TRUE 
#33   167 
#prostate
#FALSE  TRUE 
#7    81 
#liver_intrahepatic_bile_ducts
#FALSE  TRUE 
#32    41 
#esophagus
#TRUE 
#9 
#bronchus_lung
#TRUE 
#6 
#pancreas
#TRUE 
#4 
#plasma_cell_neoplasms
#TRUE 
#3 
#ovary
#TRUE 
#1 



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
  
  idx.c <- control$Species[which(control$Q<cutoff)]
  idx.d <- cancer$Species[which(cancer$Q<cutoff)]
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

merge_list_sig_Q <- lapply(merge_list_sig, function(x) x[,c("Ctr.Q", "Ex.Q")])
merge_list_sig_R2 <- lapply(merge_list_sig, function(x) x[,c("Ctr.R2", "Ex.R2")])

#head(merge_list_sig_Q[[1]])
#head(merge_list_sig_R2[[1]])

lapply(merge_list_sig_Q, dim)

#FDR <0.2 
#$bronchus_lung
#[1] 5 2
#$esophagus
#[1] 6 2
#$liver_intrahepatic_bile_ducts
#[1] 48  2
#$Lymphoid_leukemia
#[1] 128   2
#$neuroendocrine_tumors
#[1] 291   2
#$pancreas
#[1] 2 2
#$plasma_cell_neoplasms
#[1] 3 2
#$prostate
#[1] 7 2


#TRUE indicates positive effect size, higher in cancer
lapply(merge_list_sig_R2, function(x) table(x$Ctr.R2 > 0))
#$bronchus_lung
#TRUE 
#5 
#$esophagus
#TRUE 
#6 
#$liver_intrahepatic_bile_ducts
#FALSE  TRUE 
#26    22 
#$Lymphoid_leukemia
#FALSE  TRUE 
#25   103 
#$neuroendocrine_tumors
#FALSE  TRUE 
#191   100 
#$pancreas
#TRUE 
#2 
#$plasma_cell_neoplasms
#TRUE 
#3 
#$prostate
#FALSE  TRUE 
#1     6 


lapply(merge_list_sig_Q, function(x) length(which(x$Ctr.Q < 0.1)))
#$bronchus_lung
#[1] 5
#$esophagus
#[1] 6
#$liver_intrahepatic_bile_ducts
#[1] 45
#$Lymphoid_leukemia
#[1] 126
#$neuroendocrine_tumors
#[1] 274
#$pancreas
#[1] 2
#$plasma_cell_neoplasms
#[1] 3
#$prostate
#[1] 4


lapply(merge_list_sig_Q, function(x) length(which(x$Ctr.Q < 0.05)))
#$bronchus_lung
#[1] 5
#$esophagus
#[1] 6
#$liver_intrahepatic_bile_ducts
#[1] 45
#$Lymphoid_leukemia
#[1] 115
#$neuroendocrine_tumors
#[1] 248
#$pancreas
#[1] 2
#$plasma_cell_neoplasms
#[1] 3
#$prostate
#[1] 1

lapply(merge_list_sig_Q, function(x) length(which(x$Ctr.Q < 0.01)))
#$bronchus_lung
#[1] 5
#$esophagus
#[1] 6
#$liver_intrahepatic_bile_ducts
#[1] 44
#$Lymphoid_leukemia
#[1] 91
#$neuroendocrine_tumors
#[1] 126
#$pancreas
#[1] 0
#$plasma_cell_neoplasms
#[1] 2
#$prostate
#[1] 0


# Sample data
df <- data.frame(
  Name = c("bronchus lung", "esophagus", "liver intrahepatic bile ducts", "lymphoid leukemia", 
           "neuroendocrine tumors", "pancreas", "plasma cell neoplasms", "prostate"),
  "n species FDR <0.2" = c(5,6,48,128,291,2,3,7),
  "n species FDR <0.1" = c(5,6,45,126,274,2,3,4),
  "n species FDR <0.05" = c(5,6,45,115,248,2,3,1),
  "n species FDR <0.01" = c(5,6,44,91,126,0,2,0))

df <- df[order(rowSums(df[,2:ncol(df)]), decreasing = T),]

   
gt(df) %>%
  tab_header(title = "Number of species overlapping vs. healthy and vs. other cancers") %>%
  cols_label(Name="Name", n.species.FDR..0.2= "n species q<0.2", 
               n.species.FDR..0.1="n species q<0.1", n.species.FDR..0.05="n species q<0.05", 
               n.species.FDR..0.01="n species q<0.01") %>%
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


#-------------------------------------------------------------------------------

lapply(merge_list_sig_Q_w_comorb, dim)

load('./Result/CancerOnly/data.obj.wk.RData')
level <- 'Family'
otu.name <- as.data.frame(data.obj$otu.name)


names(merge_list_sig_Q_w_comorb)
#--------------------------------------
i="bronchus_lung"
#[1] "bronchus_lung"                 "esophagus"                     "liver_intrahepatic_bile_ducts" "Lymphoid_leukemia"             "neuroendocrine_tumors"        
#[6] "pancreas"                      "plasma_cell_neoplasms"         "prostate"  
FDR_df <- merge_list_sig_Q_w_comorb[[i]]
R2_df <- merge_list_sig_R2_w_comorb[[i]]
grp.names <- colnames(FDR_df)
grp.labels <- colnames(FDR_df)
grp.labels[1:2] <- c("Cancer class vs. Control", "Cancer class vs. other cancers")
sign_cutoffs <- c(0.2, 0.1, 0.05)


heatmap_list <- zicoseq_heatmap(FDR_df, R2_df, grp.names, grp.labels, sign_cutoffs, single_cutoff = FALSE, otu.name=otu.name)
heatmap_list$heatmap_plot

filename <- "./Figure/RM_figures/overlap heatmap bronchus_lung all.png"
png(filename, width =5.5*330, height = 9*330, res = 330)
heatmap_list$heatmap_plot
dev.off()


#--------------------------------------
i="esophagus"
#[1] "bronchus_lung"                 "esophagus"                     "liver_intrahepatic_bile_ducts" "Lymphoid_leukemia"             "neuroendocrine_tumors"        
#[6] "pancreas"                      "plasma_cell_neoplasms"         "prostate"  

FDR_df <- merge_list_sig_Q_w_comorb[[i]]
R2_df <- merge_list_sig_R2_w_comorb[[i]]

heatmap_list <- zicoseq_heatmap(FDR_df, R2_df, grp.names, grp.labels, sign_cutoffs, single_cutoff = FALSE, otu.name=otu.name)
heatmap_list$heatmap_plot

filename <- "./Figure/RM_figures/overlap heatmap esophagus all.png"
png(filename, width = 5.5*330, height = 9*330, res = 330)
heatmap_list$heatmap_plot
dev.off()

#Anemia is a frequent occurrence in esophageal cancer patients



#--------------------------------------
i="liver_intrahepatic_bile_ducts"
#[1] "bronchus_lung"                 "esophagus"                     "liver_intrahepatic_bile_ducts" "Lymphoid_leukemia"             "neuroendocrine_tumors"        
#[6] "pancreas"                      "plasma_cell_neoplasms"         "prostate"  

FDR_df <- merge_list_sig_Q_w_comorb[[i]]
R2_df <- merge_list_sig_R2_w_comorb[[i]]

heatmap_list <- zicoseq_heatmap(FDR_df, R2_df, grp.names, grp.labels, sign_cutoffs, single_cutoff = FALSE, otu.name=otu.name)
heatmap_list$heatmap_plot

filename <- "./Figure/RM_figures/overlap heatmap liver_intrahepatic_bile_ducts all.png"
png(filename, width = 20*330, height = 9*330, res = 330)
heatmap_list$heatmap_plot
dev.off()



#--------------------------------------
i="Lymphoid_leukemia"
#[1] "bronchus_lung"                 "esophagus"                     "liver_intrahepatic_bile_ducts" "Lymphoid_leukemia"             "neuroendocrine_tumors"        
#[6] "pancreas"                      "plasma_cell_neoplasms"         "prostate"  
FDR_df <- merge_list_sig_Q_w_comorb[[i]]
R2_df <- merge_list_sig_R2_w_comorb[[i]]
grp.names <- colnames(FDR_df)
grp.labels <- colnames(FDR_df)
grp.labels[1:2] <- c("Cancer class vs. Control", "Cancer class vs. other cancers")
sign_cutoffs <- c(0.2, 0.1, 0.05)


heatmap_list <- zicoseq_heatmap(FDR_df, R2_df, grp.names, grp.labels, sign_cutoffs, single_cutoff = FALSE, otu.name=otu.name)
heatmap_list$heatmap_plot

filename <- "./Figure/RM_figures/overlap heatmap Lymphoid_leukemia all.png"
png(filename, width =30*330, height = 9*330, res = 330)
heatmap_list$heatmap_plot
dev.off()



#--------------------------------------
i="neuroendocrine_tumors"
#[1] "bronchus_lung"                 "esophagus"                     "liver_intrahepatic_bile_ducts" "Lymphoid_leukemia"             "neuroendocrine_tumors"        
#[6] "pancreas"                      "plasma_cell_neoplasms"         "prostate"  

FDR_df <- merge_list_sig_Q_w_comorb[[i]]
R2_df <- merge_list_sig_R2_w_comorb[[i]]

heatmap_list <- zicoseq_heatmap(FDR_df, R2_df, grp.names, grp.labels, sign_cutoffs, single_cutoff = FALSE, otu.name=otu.name)
heatmap_list$heatmap_plot

filename <- "./Figure/RM_figures/overlap heatmap neuroendocrine_tumors all.png"
png(filename, width = 60*330, height = 9*330, res = 330)
heatmap_list$heatmap_plot
dev.off()


#--------------------------------------
i="pancreas" 
#[1] "bronchus_lung"                 "esophagus"                     "liver_intrahepatic_bile_ducts" "Lymphoid_leukemia"             "neuroendocrine_tumors"        
#[6] "pancreas"                      "plasma_cell_neoplasms"         "prostate"  

FDR_df <- merge_list_sig_Q_w_comorb[[i]]
R2_df <- merge_list_sig_R2_w_comorb[[i]]

heatmap_list <- zicoseq_heatmap(FDR_df, R2_df, grp.names, grp.labels, sign_cutoffs, single_cutoff = FALSE, otu.name=otu.name)
heatmap_list$heatmap_plot

filename <- "./Figure/RM_figures/overlap heatmap pancreas all.png"
png(filename, width = 4.75*330, height = 9*330, res = 330)
heatmap_list$heatmap_plot
dev.off()


#--------------------------------------
i="plasma_cell_neoplasms" 
#[1] "bronchus_lung"                 "esophagus"                     "liver_intrahepatic_bile_ducts" "Lymphoid_leukemia"             "neuroendocrine_tumors"        
#[6] "pancreas"                      "plasma_cell_neoplasms"         "prostate"  

FDR_df <- merge_list_sig_Q_w_comorb[[i]]
R2_df <- merge_list_sig_R2_w_comorb[[i]]

heatmap_list <- zicoseq_heatmap(FDR_df, R2_df, grp.names, grp.labels, sign_cutoffs, single_cutoff = FALSE, otu.name=otu.name)
heatmap_list$heatmap_plot

filename <- "./Figure/RM_figures/overlap heatmap plasma_cell_neoplasms all.png"
png(filename, width = 4.75*330, height = 9*330, res = 330)
heatmap_list$heatmap_plot
dev.off()


#--------------------------------------
i="prostate" 
#[1] "bronchus_lung"                 "esophagus"                     "liver_intrahepatic_bile_ducts" "Lymphoid_leukemia"             "neuroendocrine_tumors"        
#[6] "pancreas"                      "plasma_cell_neoplasms"         "prostate"  

FDR_df <- merge_list_sig_Q_w_comorb[[i]]
R2_df <- merge_list_sig_R2_w_comorb[[i]]

heatmap_list <- zicoseq_heatmap(FDR_df, R2_df, grp.names, grp.labels, sign_cutoffs, single_cutoff = FALSE, otu.name=otu.name)
heatmap_list$heatmap_plot

filename <- "./Figure/RM_figures/overlap heatmap prostate all.png"
png(filename, width = 5.75*330, height = 9*330, res = 330)
heatmap_list$heatmap_plot
dev.off()


#-------------------------------------------------------------------------------
#cancer classes of interest

#[1] "bronchus_lung"                 "esophagus"                     "liver_intrahepatic_bile_ducts" "Lymphoid_leukemia"             "neuroendocrine_tumors"        
#[6] "pancreas"                      "plasma_cell_neoplasms"         "prostate"  

#clin_meta$icd10_first_3_name[grep("plasma", clin_meta$icd10_first_3_name)]

cancer_class_w_sign <- c("Malignant neoplasm of bronchus and lung", "Malignant neoplasm of esophagus", "Malignant neoplasm of liver and intrahepatic bile ducts",
                         "Lymphoid leukemia", "Malignant neuroendocrine tumors", "Malignant neoplasm of pancreas", "Multiple myeloma and malignant plasma cell neoplasms",
                         "Malignant neoplasm of prostate")

cancer_class_w_sign_row_list <- sapply(cancer_class_w_sign, function(x) which(clin_meta$icd10_first_3_name == x))

#tabulate all Elixhauser components with a significant species
test_df <- lapply(merge_list_sig_Q_w_comorb, function(x) (x[,3:ncol(x)] < 0.2)) #0.2 is FDR cutoff

#general one
Elix_components_w_sig <- unique(unlist(lapply(test_df, function(x) colnames(x)[colSums(x == TRUE) > 0])))
#[1] "Elixhauser_Anemia"      "Elixhauser_FluidsLytes" "Elixhauser_Liver"       "Elixhauser_PUD"         "Elixhauser_Renal"       "Elixhauser_Valvular"    "Elixhauser_WeightLoss" 
#[8] "Elixhauser_DMcx"        "Elixhauser_Arrhythmia"  "Elixhauser_DM"          "Elixhauser_NeuroOther"  "Elixhauser_Hypothyroid" "Elixhauser_Pulmonary" 


#make Elixhauser components specific to that cancer class
Elix_components_w_sig_list <- lapply(test_df, function(x) colnames(x)[colSums(x == TRUE) > 0])


component_list <- list()
cancer_class_w_sign_row_list_sub <- cancer_class_w_sign_row_list[lapply(Elix_components_w_sig_list, length) > 0] 
Elix_components_w_sig_list_sub <- Elix_components_w_sig_list[lapply(Elix_components_w_sig_list, length) > 0]


#the 2 pancreas species do not have associations with comorbidities
names(cancer_class_w_sign_row_list_sub)
#[1] "Malignant neoplasm of bronchus and lung"                 "Malignant neoplasm of esophagus"                         "Malignant neoplasm of liver and intrahepatic bile ducts"
#[4] "Lymphoid leukemia"                                       "Malignant neuroendocrine tumors"                         "Multiple myeloma and malignant plasma cell neoplasms"   
#[7] "Malignant neoplasm of prostate"      


#-------------------------------------------------------------------------------
#for cancer groups that have only very few patients with these comorbidities these comorbidities do not have to be plotted or considered
#only consider more than >15%



#the logic here is twofold, 1) a comorbidity that predisposes to cancer can drive the microbiome change
#2) or the comorbidity can be a symptom of the cancer and drive a microbiome change at the time of sampling. 
#Both should be removed when searching for true cancer associations that may predispose to or drive disease progression



#-------------------------------------------------------------------------------

#[1] "Malignant neoplasm of bronchus and lung"                 "Malignant neoplasm of esophagus"                         "Malignant neoplasm of liver and intrahepatic bile ducts"
#[4] "Lymphoid leukemia"                                       "Malignant neuroendocrine tumors"                         "Multiple myeloma and malignant plasma cell neoplasms"   
#[7] "Malignant neoplasm of prostate"      

#-------------------------------------------------------------------------------
#"Malignant neoplasm of bronchus and lung"
i=1
length(cancer_class_w_sign_row_list_sub[[i]]) #90
temp_table <- apply(clin_meta[cancer_class_w_sign_row_list_sub[[i]], Elix_components_w_sig_list_sub[[i]]], 2, table)
#             Elixhauser_Anemia Elixhauser_FluidsLytes Elixhauser_Liver Elixhauser_PUD Elixhauser_Renal Elixhauser_Valvular Elixhauser_WeightLoss
#No                 83                     71               72             89               79                  73                    74
#Yes                 7                     19               18              1               11                  17                    16


#get the percentages
temp_table / 90 * 100
#         Elixhauser_Anemia Elixhauser_FluidsLytes Elixhauser_Liver Elixhauser_PUD Elixhauser_Renal Elixhauser_Valvular Elixhauser_WeightLoss
#No          92.222222               78.88889               80      98.888889         87.77778            81.11111              82.22222
#es          7.777778               21.11111               20       1.111111         12.22222            18.88889              17.77778

#>15%
selected_elix_lung <- c("Elixhauser_FluidsLytes", "Elixhauser_Liver", "Elixhauser_Valvular", "Elixhauser_WeightLoss")


#----------------------------------
#"Malignant neoplasm of esophagus" 
i=2
length(cancer_class_w_sign_row_list_sub[[i]]) #34
temp_table <- apply(clin_meta[cancer_class_w_sign_row_list_sub[[i]], Elix_components_w_sig_list_sub[[i]]], 2, table)
#         Elixhauser_Anemia Elixhauser_DMcx Elixhauser_Liver Elixhauser_PUD Elixhauser_Renal Elixhauser_Valvular
#No                 30              32               25             31               27                  30
#Yes                 4               2                9              3                7                   4

temp_table / 34 * 100
#       Elixhauser_Anemia Elixhauser_DMcx Elixhauser_Liver Elixhauser_PUD Elixhauser_Renal Elixhauser_Valvular
#No           88.23529       94.117647         73.52941      91.176471         79.41176            88.23529
#Yes          11.76471        5.882353         26.47059       8.823529         20.58824            11.76471


#anemia has a lot of associations, but only 4/34 patients have this. Anemia is a frequent occurrence in esophageal cancer patients

selected_elix_esop <- c("Elixhauser_Liver", "Elixhauser_Renal")


#----------------------------------
#"Malignant neoplasm of liver and intrahepatic bile ducts"
i=3
length(cancer_class_w_sign_row_list_sub[[i]]) #39
temp_table <- do.call(cbind, apply(clin_meta[cancer_class_w_sign_row_list_sub[[i]], Elix_components_w_sig_list_sub[[i]]], 2, table))
#         Elixhauser_Arrhythmia Elixhauser_DM Elixhauser_FluidsLytes Elixhauser_Liver Elixhauser_NeuroOther Elixhauser_PUD Elixhauser_WeightLoss
#No                     30            34                     30               19                    39            39                    32
#Yes                     9             5                      9               20                    0             0                     7

temp_table / 39 * 100
#         Elixhauser_Arrhythmia Elixhauser_DM Elixhauser_FluidsLytes Elixhauser_Liver Elixhauser_NeuroOther Elixhauser_PUD Elixhauser_WeightLoss
#No               76.92308      87.17949               76.92308         48.71795                   100            100              82.05128
#Yes              23.07692      12.82051               23.07692         51.28205                   100            100              17.94872


#FluidsLytes, Liver, Arrhythmia, and Weightloss are interesting
#a lot have liver which makes sense and should be plotted as examples.

#perplexity prompt #are Malignant neoplasm of liver and intrahepatic bile ducts linked to arrhythmia's?


#no patients have PUD and NeuroOther so those should be removed from the plot and analysis?

#weight loss can stay but these species should be plotted stratified by weightloss; could be interesting; 17%
#weight loss is relevant as that is common for this cancer, leads to decreased appetite

#liver is obvious so good example to be plotted

selected_elix_liv <- c("Elixhauser_Arrhythmia", "Elixhauser_FluidsLytes", "Elixhauser_Liver", "Elixhauser_WeightLoss")


#-------------------------------------------------------------------------------
#"Lymphoid leukemia"
i=4
length(cancer_class_w_sign_row_list_sub[[i]]) #17
temp_table <- apply(clin_meta[cancer_class_w_sign_row_list_sub[[i]], Elix_components_w_sig_list_sub[[i]]], 2, table)
#     Elixhauser_Anemia Elixhauser_Arrhythmia Elixhauser_DMcx Elixhauser_FluidsLytes Elixhauser_Liver Elixhauser_NeuroOther Elixhauser_PUD Elixhauser_WeightLoss
#No                 16                    12              16                     14               15                    14             16                    16
#Yes                 1                     5               1                      3                2                     3              1                     1

#get the percentages
temp_table / 17 * 100
#       Elixhauser_Anemia Elixhauser_Arrhythmia Elixhauser_DMcx Elixhauser_FluidsLytes Elixhauser_Liver Elixhauser_NeuroOther Elixhauser_PUD Elixhauser_WeightLoss
#No          94.117647              70.58824       94.117647               82.35294         88.23529              82.35294      94.117647             94.117647
#Yes          5.882353              29.41176        5.882353               17.64706         11.76471              17.64706       5.882353              5.882353


#Remove Lymphoma all together since we have another paper on this?
#Arrhythmia has the only decent number of patients that seem relevant.

#>15%
selected_elix_lymp <- c("Elixhauser_Arrhythmia", "Elixhauser_FluidsLytes", "Elixhauser_NeuroOther")


#----------------------------------
#"Malignant neuroendocrine tumors"
i=5
length(cancer_class_w_sign_row_list_sub[[i]]) #21
temp_table <- do.call(cbind, apply(clin_meta[cancer_class_w_sign_row_list_sub[[i]], Elix_components_w_sig_list_sub[[i]]], 2, table))
#         Elixhauser_Anemia Elixhauser_Arrhythmia Elixhauser_DM Elixhauser_DMcx Elixhauser_FluidsLytes Elixhauser_Hypothyroid Elixhauser_Liver Elixhauser_NeuroOther
#No                 19                    14            21              20                     18                     17               14                    21
#Yes                2                     7             0               1                      3                      4                7                    0
#         Elixhauser_PUD Elixhauser_Pulmonary Elixhauser_Renal Elixhauser_Valvular Elixhauser_WeightLoss
#No              21                   18               19                  18                    18
#Yes             0                    3                2                   3                     3


temp_table / 21 * 100
#       Elixhauser_Anemia Elixhauser_Arrhythmia Elixhauser_DM Elixhauser_DMcx Elixhauser_FluidsLytes Elixhauser_Hypothyroid Elixhauser_Liver Elixhauser_NeuroOther
#No           90.47619              66.66667           100       95.238095               85.71429               80.95238         66.66667                   100
#Yes           9.52381              33.33333           0          4.761905               14.28571               19.04762         33.33333                    0
#       Elixhauser_PUD Elixhauser_Pulmonary Elixhauser_Renal Elixhauser_Valvular Elixhauser_WeightLoss
#No             100             85.71429         90.47619            85.71429              85.71429
#Yes             0              14.28571          9.52381            14.28571              14.28571

#which ones are relevant
#Arrhythmia, Liver, and maybe Hypothyroid are relevant
#liver has some good examples that should be highlighted as clear examples that should be removed
#liver involvement through metastasis has been reported but not baseline liver disease


#About Arrhythmia's
#not all patients with malignant neuroendocrine tumors will experience arrhythmias, there is evidence suggesting a link 
#between NETs and various types of arrhythmias, particularly in the context of carcinoid syndrome. The mechanisms may 
#involve hormone release, sympathetic activation, and in some cases, direct cardiac metastasis. 

selected_elix_neuro <- c("Elixhauser_Arrhythmia", "Elixhauser_Hypothyroid", "Elixhauser_Liver")



#----------------------------------
#"Multiple myeloma and malignant plasma cell neoplasms"

i=6
length(cancer_class_w_sign_row_list_sub[[i]]) #60
temp_table <- apply(clin_meta[cancer_class_w_sign_row_list_sub[[i]], Elix_components_w_sig_list_sub[[i]]], 2, table)
#             Elixhauser_Arrhythmia Elixhauser_WeightLoss
#No                     42                    54
#Yes                    18                     6

#get the percentages
temp_table / 60 * 100
#             Elixhauser_Arrhythmia Elixhauser_WeightLoss
#No                     70                    90
#Yes                    30                    10

#Remove Lymphoma all together since we have another paper on this?
#Arrhythmia has the only decent number of patients that seem relevant.

#>15%
selected_elix_plasmacell <- c("Elixhauser_Arrhythmia")


#----------------------------------
# "Malignant neoplasm of prostate"   

i=7
length(cancer_class_w_sign_row_list_sub[[i]]) #119
temp_table <- apply(clin_meta[cancer_class_w_sign_row_list_sub[[i]], Elix_components_w_sig_list_sub[[i]]], 2, table)
#           Elixhauser_Hypothyroid Elixhauser_WeightLoss
#No                     108                   115
#Yes                     11                     4

#get the percentages
temp_table / 119 * 100
#         Elixhauser_Hypothyroid Elixhauser_WeightLoss
#No               90.756303             96.638655
#Yes               9.243697              3.361345

#Very few patients ~3%. Not important and these associations seem legit


#-------------------------------------------------------------------------------
#plot using the subset of the Elixhauser categories

#1 selected_elix_lung
#2 selected_elix_esop
#3 selected_elix_liv
#4 selected_elix_lymp
#5 selected_elix_neuro
#6 selected_elix_plasmacell
#7 prostate has 0

Elix_components_w_sig_list_sub_cutoff <- Elix_components_w_sig_list_sub

#Bronchus lung
Elix_components_w_sig_list_sub_cutoff[[1]] <- Elix_components_w_sig_list_sub[[1]][Elix_components_w_sig_list_sub[[1]] %in% 
                                                                                    selected_elix_lung]
#esophagus
Elix_components_w_sig_list_sub_cutoff[[2]] <- Elix_components_w_sig_list_sub[[2]][Elix_components_w_sig_list_sub[[2]] %in% 
                                                                                    selected_elix_esop]
#liver_intrahepatic_bile_ducts
Elix_components_w_sig_list_sub_cutoff[[3]] <- Elix_components_w_sig_list_sub[[3]][Elix_components_w_sig_list_sub[[3]] %in% 
                                                                                    selected_elix_liv]
#Lymphoid_leukemia
Elix_components_w_sig_list_sub_cutoff[[4]] <- Elix_components_w_sig_list_sub[[4]][Elix_components_w_sig_list_sub[[4]] %in% 
                                                                                    selected_elix_lymp]
#neuroendocrine_tumors
Elix_components_w_sig_list_sub_cutoff[[5]] <- Elix_components_w_sig_list_sub[[5]][Elix_components_w_sig_list_sub[[5]] %in% 
                                                                                    selected_elix_neuro]
#plasma cell
Elix_components_w_sig_list_sub_cutoff[[6]] <- Elix_components_w_sig_list_sub[[6]][Elix_components_w_sig_list_sub[[6]] %in% 
                                                                                    selected_elix_plasmacell]
#prostate
Elix_components_w_sig_list_sub_cutoff[[7]] <- Elix_components_w_sig_list_sub[[7]][Elix_components_w_sig_list_sub[[7]] %in% 
                                                                                    c("")]

#pancreas does not have any comorbidities


#--------------------------------------

sign_lists_to_plot <- list()


#--------------------------------------
i=1
#bronchus_lung
names(merge_list_sig_Q_w_comorb)[i]
sig_components <- Elix_components_w_sig_list_sub_cutoff[[i]]
#"Elixhauser_FluidsLytes" "Elixhauser_Liver"       "Elixhauser_Valvular"    "Elixhauser_WeightLoss"
FDR_df <- merge_list_sig_Q_w_comorb[[i]][,c("Ctr.Q","Ex.Q",sig_components)]
R2_df <- merge_list_sig_R2_w_comorb[[i]][,c("Ctr.R2","Ex.R2",sig_components)]

grp.names <- colnames(FDR_df)
grp.labels <- colnames(FDR_df)
grp.labels[1:2] <- c("Cancer class vs. Control", "Cancer class vs. other cancers")
sign_cutoffs <- c(0.2, 0.1, 0.05)

heatmap_list <- zicoseq_heatmap(FDR_df, R2_df, grp.names, grp.labels, sign_cutoffs, single_cutoff = FALSE, otu.name=otu.name, column_name_rot=40)
heatmap_list$heatmap_plot

filename <- "./Figure/RM_figures/overlap heatmap bronchus_lung_sig_Elix.png"
png(filename, width = 7.25*330, height = 5*330, res = 330)
heatmap_list$heatmap_plot
dev.off()

colnames(FDR_df) <- grp.labels
sign_lists_to_plot$bronchus_lung <- FDR_df


#this below does not work because there is only a single species without comorbidities


#remove Elix species; any rows with FDR <0.2 beyond the first 2 columns
#to_keep_rows <- !rowSums(FDR_df[, 3:ncol(FDR_df)] < 0.2) > 0
#length(which(to_keep_rows)) #5

#heatmap_list <- zicoseq_heatmap(FDR_df[to_keep_rows,], R2_df[to_keep_rows, ], grp.names, grp.labels, sign_cutoffs, single_cutoff = FALSE, otu.name=otu.name)

#filename <- "./Figure/RM_figures/overlap heatmap bronchus_lung_sig_Elix_rm.png"
#png(filename, width = 5.5*330, height = 5*330, res = 330)
#heatmap_list$heatmap_plot
#dev.off()


#--------------------------------------
i=2
#esophagus
names(merge_list_sig_Q_w_comorb)[i]
#"Elixhauser_Liver" "Elixhauser_Renal"
sig_components <- Elix_components_w_sig_list_sub_cutoff[[i]]
FDR_df <- merge_list_sig_Q_w_comorb[[i]][,c("Ctr.Q","Ex.Q",sig_components)]
R2_df <- merge_list_sig_R2_w_comorb[[i]][,c("Ctr.R2","Ex.R2",sig_components)]

grp.names <- colnames(FDR_df)
grp.labels <- colnames(FDR_df)
grp.labels[1:2] <- c("Cancer class vs. Control", "Cancer class vs. other cancers")
sign_cutoffs <- c(0.2, 0.1, 0.05)

heatmap_list <- zicoseq_heatmap(FDR_df, R2_df, grp.names, grp.labels, sign_cutoffs, single_cutoff = FALSE, otu.name=otu.name, column_name_rot=40)
heatmap_list$heatmap_plot

filename <- "./Figure/RM_figures/overlap heatmap esophagus_sig_Elix.png"
png(filename, width = 7.5*330, height = 4.75*330, res = 330)
heatmap_list$heatmap_plot
dev.off()

colnames(FDR_df) <- grp.labels
sign_lists_to_plot$esophagus <- FDR_df


#remove Elix species; any rows with FDR <0.2 beyond the first 2 columns
to_keep_rows <- !rowSums(FDR_df[, 3:ncol(FDR_df)] < 0.2) > 0
length(which(to_keep_rows)) #5

heatmap_list <- zicoseq_heatmap(FDR_df[to_keep_rows,], R2_df[to_keep_rows, ], grp.names, grp.labels, sign_cutoffs, single_cutoff = FALSE, 
                                otu.name=otu.name, , column_name_rot=40)

filename <- "./Figure/RM_figures/overlap heatmap esophagus_sig_Elix_rm.png"
png(filename, width = 7.25*330, height = 4.75*330, res = 330)
heatmap_list$heatmap_plot
dev.off()



#--------------------------------------
i=3
#liver_intrahepatic_bile_ducts
names(merge_list_sig_Q_w_comorb)[i]

sig_components <- Elix_components_w_sig_list_sub_cutoff[[i]]
FDR_df <- merge_list_sig_Q_w_comorb[[i]][,c("Ctr.Q","Ex.Q",sig_components)]
R2_df <- merge_list_sig_R2_w_comorb[[i]][,c("Ctr.R2","Ex.R2",sig_components)]

grp.names <- colnames(FDR_df)
grp.labels <- colnames(FDR_df)
grp.labels[1:2] <- c("Cancer class vs. Control", "Cancer class vs. other cancers")
sign_cutoffs <- c(0.2, 0.1, 0.05)

heatmap_list <- zicoseq_heatmap(FDR_df, R2_df, grp.names, grp.labels, sign_cutoffs, single_cutoff = FALSE, otu.name=otu.name, column_name_rot = 40)
heatmap_list$heatmap_plot

filename <- "./Figure/RM_figures/overlap heatmap liver_intrahepatic_bile_ducts_sig_Elix.png"
png(filename, width = 26*330, height = 6.5*330, res = 330)
heatmap_list$heatmap_plot
dev.off()

colnames(FDR_df) <- grp.labels
sign_lists_to_plot$liver_intrahepatic_bile_ducts <- FDR_df


#remove Elix species; any rows with FDR <0.2 beyond the first 2 columns
to_keep_rows <- !rowSums(FDR_df[, 3:ncol(FDR_df)] < 0.2) > 0
length(which(to_keep_rows)) #19

heatmap_list <- zicoseq_heatmap(FDR_df[to_keep_rows,], R2_df[to_keep_rows, ], grp.names, grp.labels, sign_cutoffs, single_cutoff = FALSE, 
                                otu.name=otu.name, column_name_rot = 40)
heatmap_list$heatmap_plot

filename <- "./Figure/RM_figures/overlap heatmap liver_intrahepatic_bile_ducts_sig_Elix_rm.png"
png(filename, width = 16*330, height = 6.5*330, res = 330)
heatmap_list$heatmap_plot
dev.off()


#--------------------------------------
i=4
names(merge_list_sig_Q_w_comorb)[i]
#Lymphoid_leukemia
sig_components <- Elix_components_w_sig_list_sub_cutoff[[i]]
#"Elixhauser_Arrhythmia"  "Elixhauser_FluidsLytes" "Elixhauser_NeuroOther" 
FDR_df <- merge_list_sig_Q_w_comorb[[i]][,c("Ctr.Q","Ex.Q",sig_components)]
R2_df <- merge_list_sig_R2_w_comorb[[i]][,c("Ctr.R2","Ex.R2",sig_components)]

grp.names <- colnames(FDR_df)
grp.labels <- colnames(FDR_df)
grp.labels[1:2] <- c("Cancer class vs. Control", "Cancer class vs. other cancers")
sign_cutoffs <- c(0.2, 0.1, 0.05)


heatmap_list <- zicoseq_heatmap(FDR_df, R2_df, grp.names, grp.labels, sign_cutoffs, single_cutoff = FALSE, otu.name=otu.name, column_name_rot = 40)
heatmap_list$heatmap_plot

filename <- "./Figure/RM_figures/overlap heatmap Lymphoid_leukemia_sig_Elix.png"
png(filename, width =38*330, height = 6*330, res = 330)
heatmap_list$heatmap_plot
dev.off()

colnames(FDR_df) <- grp.labels
sign_lists_to_plot$Lymphoid_leukemia <- FDR_df


#remove Elix species; any rows with FDR <0.2 beyond the first 2 columns
to_keep_rows <- !rowSums(FDR_df[, 3:ncol(FDR_df)] < 0.2) > 0
length(which(to_keep_rows)) #99

heatmap_list <- zicoseq_heatmap(FDR_df[to_keep_rows,], R2_df[to_keep_rows, ], grp.names, grp.labels, sign_cutoffs, single_cutoff = FALSE, 
                                otu.name=otu.name, column_name_rot = 40)

filename <- "./Figure/RM_figures/overlap heatmap Lymphoid_leukemia_sig_Elix_rm.png"
png(filename, width =34*330, height = 6*330, res = 330)
heatmap_list$heatmap_plot
dev.off()



#--------------------------------------
i=5
#neuroendocrine_tumors
names(merge_list_sig_Q_w_comorb)[i]

sig_components <- Elix_components_w_sig_list_sub_cutoff[[i]]
FDR_df <- merge_list_sig_Q_w_comorb[[i]][,c("Ctr.Q","Ex.Q",sig_components)]
R2_df <- merge_list_sig_R2_w_comorb[[i]][,c("Ctr.R2","Ex.R2",sig_components)]

grp.names <- colnames(FDR_df)
grp.labels <- colnames(FDR_df)
grp.labels[1:2] <- c("Cancer class vs. Control", "Cancer class vs. other cancers")
sign_cutoffs <- c(0.2, 0.1, 0.05)

heatmap_list <- zicoseq_heatmap(FDR_df, R2_df, grp.names, grp.labels, sign_cutoffs, single_cutoff = FALSE, otu.name=otu.name, column_name_rot = 40)
heatmap_list$heatmap_plot

filename <- "./Figure/RM_figures/overlap heatmap neuroendocrine_tumors_sig_Elix.png"
png(filename, width = 62*330, height = 6*330, res = 330)
heatmap_list$heatmap_plot
dev.off()

colnames(FDR_df) <- grp.labels
sign_lists_to_plot$neuroendocrine_tumors <- FDR_df


#test with top 25% effect size cutoff
heatmap_list <- zicoseq_heatmap(FDR_df, R2_df, grp.names, grp.labels, sign_cutoffs, single_cutoff = FALSE, 
                                otu.name=otu.name, effect_size_cutoff = 0.25, column_name_rot = 40)
heatmap_list$heatmap_plot

filename <- "./Figure/RM_figures/overlap heatmap neuroendocrine_tumors_sig_Elix_cutoff_top25perc.png"
png(filename, width = 32*330, height = 6*330, res = 330)
heatmap_list$heatmap_plot
dev.off()



#remove Elix species; any rows with FDR <0.2 beyond the first 2 columns
to_keep_rows <- !rowSums(FDR_df[, 3:ncol(FDR_df)] < 0.2) > 0
length(which(to_keep_rows)) #265

heatmap_list <- zicoseq_heatmap(FDR_df[to_keep_rows,], R2_df[to_keep_rows, ], grp.names, grp.labels, sign_cutoffs, single_cutoff = FALSE, 
                                otu.name=otu.name, column_name_rot = 40)

filename <- "./Figure/RM_figures/overlap heatmap neuroendocrine_tumors_sig_Elix_rm.png"
png(filename, width = 62*330, height = 6.5*330, res = 330)
heatmap_list$heatmap_plot
dev.off()


#--------------------------------------
i=6 
#pancreas; no comorbidity associations
names(merge_list_sig_Q_w_comorb)[i]
sig_components <- ""
#"Elixhauser_Arrhythmia"
FDR_df <- merge_list_sig_Q_w_comorb[[i]][,c("Ctr.Q","Ex.Q")]
R2_df <- merge_list_sig_R2_w_comorb[[i]][,c("Ctr.R2","Ex.R2")]

grp.names <- colnames(FDR_df)
grp.labels <- colnames(FDR_df)
grp.labels[1:2] <- c("Cancer class vs. Control", "Cancer class vs. other cancers")
sign_cutoffs <- c(0.2, 0.1, 0.05)


heatmap_list <- zicoseq_heatmap(FDR_df, R2_df, grp.names, grp.labels, sign_cutoffs, single_cutoff = FALSE, otu.name=otu.name, column_name_rot = 40)
heatmap_list$heatmap_plot

filename <- "./Figure/RM_figures/overlap heatmap pancreas_sig_Elix.png"
png(filename, width = 6.5*330, height = 4*330, res = 330)
heatmap_list$heatmap_plot
dev.off()

#2 species


#--------------------------------------
i="plasma_cell_neoplasms"
#plasma cell neoplasms
sig_components <- Elix_components_w_sig_list_sub_cutoff[[i]]
#"Elixhauser_Arrhythmia"
FDR_df <- merge_list_sig_Q_w_comorb[[i]][,c("Ctr.Q","Ex.Q",sig_components)]
R2_df <- merge_list_sig_R2_w_comorb[[i]][,c("Ctr.R2","Ex.R2",sig_components)]

grp.names <- colnames(FDR_df)
grp.labels <- colnames(FDR_df)
grp.labels[1:2] <- c("Cancer class vs. Control", "Cancer class vs. other cancers")
sign_cutoffs <- c(0.2, 0.1, 0.05)

heatmap_list <- zicoseq_heatmap(FDR_df, R2_df, grp.names, grp.labels, sign_cutoffs, single_cutoff = FALSE, otu.name=otu.name, 
                                column_name_rot = 40, left_margin = 6)
heatmap_list$heatmap_plot

filename <- "./Figure/RM_figures/overlap heatmap plasma cell_sig_Elix.png"
png(filename, width = 8*330, height = 5.5*330, res = 330)
heatmap_list$heatmap_plot
dev.off()

colnames(FDR_df) <- grp.labels
sign_lists_to_plot$plasma_cell_neoplasms <- FDR_df


#--------------------------------------
i="prostate"
#prostate; no comorbidity associations
sig_components <- Elix_components_w_sig_list_sub_cutoff[[i]]
FDR_df <- merge_list_sig_Q_w_comorb[[i]][,c("Ctr.Q","Ex.Q",sig_components)]
R2_df <- merge_list_sig_R2_w_comorb[[i]][,c("Ctr.R2","Ex.R2",sig_components)]

grp.names <- colnames(FDR_df)
grp.labels <- colnames(FDR_df)
grp.labels[1:2] <- c("Cancer class vs. Control", "Cancer class vs. other cancers")
sign_cutoffs <- c(0.2, 0.1, 0.05)

heatmap_list <- zicoseq_heatmap(FDR_df, R2_df, grp.names, grp.labels, sign_cutoffs, single_cutoff = FALSE, otu.name=otu.name, column_name_rot = 40)
heatmap_list$heatmap_plot

filename <- "./Figure/RM_figures/overlap heatmap prostate_sig_Elix.png"
png(filename, width = 9*330, height = 4.5*330, res = 330)
heatmap_list$heatmap_plot
dev.off()



#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#plot species abundances stratified by comorbidity

#boxplot; healthy control, all cancer, all cancer + comorbid, cancer class, cancer class + comorbid

clin_meta$icd10_first_3_name[is.na(clin_meta$icd10_first_3_name)] <- "Healthy"
head(sort(table(clin_meta$icd10_first_3_name), decreasing = T), 10) #testing

#rff relative abundance, sweep
tax_table_RA <- sweep(tax_table, MARGIN = 2, colSums(tax_table), '/') * 100

#update cancer_row_list in same order as sign_lists_to_plot
names(sign_lists_to_plot)

group_vec <- clin_meta$icd10_first_3_name
cancer_row_list <- list(
  grep("bronchus", group_vec),
  grep("esophagus", group_vec),
  grep("liver and intrahepatic bile ducts", group_vec),
  grep("Lymphoid leukemia", group_vec), 
  grep("neuroendocrine", group_vec),
  grep("plasma cell", group_vec))


elix_row_names_list <- lapply(sign_lists_to_plot, function(x) colnames(x)[3:ncol(x)])
healthy_rows <- which(group_vec == "Healthy")

sel_cols <- brewer.pal(5, "Dark2")

for (i in 1:length(sign_lists_to_plot)) { #cancer class
  sig_df <- sign_lists_to_plot[[i]]
  elix_row_names <- elix_row_names_list[[i]]
  cancer_group_rows <- cancer_row_list[[i]]
  all_other_cancer_rows <- which(!1:ncol(tax_table_RA) %in% c(healthy_rows, cancer_group_rows))
  cancer_class <- names(sign_lists_to_plot)[i]
  
  pdf(paste0("./Figure/RM_figures/Comorbidity species_", cancer_class, ".pdf"), width=5, height=4.5)
  for (j in 1:length(elix_row_names)) { #elix component
    sig_species <- rownames(sig_df)[which(sig_df[elix_row_names[j]] < 0.2)]
    sig_species <- gsub(" ", "_", sig_species)
    
    for (k in 1:length(sig_species)) {
      
      temp_ra <- sqrt(as.numeric(tax_table_RA[sig_species[k],]))
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
  
  pdf(paste0("./Figure/RM_figures/Comorbidity species split_", cancer_class, ".pdf"), width=5.5, height=4.5)
  for (j in 1:length(elix_row_names)) { #elix component
    sig_species <- rownames(sig_df)[which(sig_df[elix_row_names[j]] < 0.2)]
    sig_species <- gsub(" ", "_", sig_species)
    
    for (k in 1:length(sig_species)) {
      
      temp_ra <- sqrt(as.numeric(tax_table_RA[sig_species[k],]))
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
      
      mtext("sqrt(relative abundance %)", 2, line = 2.5)
      mtext(paste0(cancer_class, " _ ", elix_row_names[j], "\n", sig_species[k]), 3, line = 0, cex=1)
      
    }
  }
  dev.off()
}




#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#make plots for all interesting species as well without factoring in comorbidities

cancer_row_list <- list(
  grep("bronchus", group_vec),
  grep("esophagus", group_vec),
  grep("liver and intrahepatic bile ducts", group_vec),
  grep("Lymphoid leukemia", group_vec), 
  grep("neuroendocrine", group_vec),
  grep("pancreas", group_vec),
  grep("plasma cell", group_vec),
  grep("prostate", group_vec))

names(merge_list_sig_Q_w_comorb)

for (i in 1:length(merge_list_sig_Q_w_comorb)) { #cancer class
  sig_df <- merge_list_sig_Q_w_comorb[[i]]
  
  cancer_group_rows <- cancer_row_list[[i]]
  all_other_cancer_rows <- which(!1:nrow(tax_table_RA) %in% c(healthy_rows, cancer_group_rows))
  cancer_class <- names(merge_list_sig_Q_w_comorb)[i]
  
  sig_species <- rownames(sig_df)
  sig_species <- gsub(" ", "_", sig_species)
  
  pdf(paste0("./Figure/RM_figures/Significant species_scatter_", cancer_class, ".pdf"), width=5, height=4.5)
  for (k in 1:length(sig_species)) {
    
    temp_ra <- sqrt(as.numeric(tax_table_RA[sig_species[k],]))
    plot_df <- as.data.frame(cbind(c(rep("healthy", times=length(healthy_rows)), rep("all other cancer", times=length(all_other_cancer_rows)), rep("cancer_group", times=length(cancer_group_rows))),
          c(temp_ra[healthy_rows], temp_ra[all_other_cancer_rows], temp_ra[cancer_group_rows])))
    names(plot_df) <- c("group", "ra")
    plot_df$ra <- as.numeric(plot_df$ra)
    
    p_plot <-  ggboxplot(plot_df, "group", "ra",
                         color = "group", palette=sel_cols,
                         add = "jitter", ylab= "sqrt(relative_abundance %)") +
      scale_y_continuous(trans='sqrt') +
      ggtitle(sig_species[k]) +
      theme(legend.position="none")
    
    print(p_plot)
  }
  dev.off()
}




#-------------------------------------------------------------------------------
#For all cancer classes and a single significance level <0.2. Show number of species vs control, vs other cancers, and overlap.
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
n_wo_sig_elix <- rep(0, times= length(names(merge_list)))
names(n_wo_sig_elix) <- names(merge_list)

n_wo_sig_elix["bronchus_lung"] <- 1
n_wo_sig_elix["Lymphoid_leukemia"] <- 99
n_wo_sig_elix["esophagus"] <- 5
n_wo_sig_elix["liver_intrahepatic_bile_ducts"] <- 19
n_wo_sig_elix["neuroendocrine_tumors"] <- 265
n_wo_sig_elix["pancreas"] <- 2
n_wo_sig_elix["plasma_cell_neoplasms"] <- 2
n_wo_sig_elix["prostate"] <- 7

n_wo_sig_elix


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

plot_df <- plot_df[order(plot_df$`overlap w/o comorbidities`, decreasing = T),]

gt(plot_df) %>%
  tab_header(title = "Number of species when correcting for comorbidities") %>%
  opt_vertical_padding(scale = 0.1)


#make a version with FDR 0.1 cutoff?



#-------------------------------------------------------------------------------
#get numbers for the main paper; check how many of the consistently cancer associated species are also associated to a comorbidity


cutoff <- 0.1

merge_list <- vector("list", length = length(canX_Ctr_list))
names(merge_list) <- names(canX_Ctr_list)

sig_control_list_2 <- vector("list", length = length(canX_Ctr_list))


#merge these in a for loop
for (i in 1:length(canX_Ctr_list)) {
  
  control <- canX_Ctr_list[[i]]
  idx.c <- control$Species[which(control$Q<cutoff)]
  control_sub <- control[idx.com,c("P", "R2", "Q")]
  sig_control_list_2[[i]] <- idx.c
}

lapply(sig_control_list_2, length)
median(as.numeric(lapply(sig_control_list_2, length)))


length(unique(unlist(sig_control_list_2)))
#1060

#23 cancer classes now
12 /length(sig_control_list_2) #52.2%
length(which(table(unlist(sig_control_list_2)) >= 12))
#249 species


21 /length(sig_control_list_2) #91.3%
length(which(table(unlist(sig_control_list_2)) >= 21))
#33; make heatmap of these species
names(which(table(unlist(sig_control_list_2)) >= 21))
# [1] "s__Anaerobutyricum sp900016875"    "s__Anaerostipes caccae"            "s__Anaerostipes sp000508985"       "s__Bifidobacterium dentium"        "s__Blautia coccoides"             
#[6] "s__Blautia sp000432195"            "s__Blautia sp001304935"            "s__Blautia sp900120295"            "s__Blautia sp900753905"            "s__Blautia_A caecimuris"          
#[11] "s__Blautia_A sp900551465"          "s__Clostridium_AQ innocuum"        "s__Clostridium_Q sp003024715"      "s__Enterocloster clostridioformis" "s__ER4 sp000765235"               
#[16] "s__Erysipelatoclostridium ramosum" "s__Escherichia coli"               "s__Escherichia fergusonii"         "s__Faecalibacterium prausnitzii_C" "s__Faecalibacterium prausnitzii_I"
#[21] "s__Faecalibacterium sp900539885"   "s__Faecalibacterium sp900765705"   "s__Faecalimonas sp900556835"       "s__GCA-900066135 sp900543575"      "s__Klebsiella pneumoniae"         
#[26] "s__Ruminococcus_B gnavus"          "s__Streptococcus parasanguinis_A"  "s__Streptococcus parasanguinis_B"  "s__Streptococcus sp000448565"      "s__Streptococcus sp001556435"     
#[31] "s__Streptococcus sp900766505"      "s__Streptococcus vestibularis"     "s__UMGS403 sp900540275"         

#all comparisons
length(which(table(unlist(sig_control_list_2)) >= 23))
#5
names(which(table(unlist(sig_control_list_2)) >= 23))
#[1] "s__Bifidobacterium dentium"        "s__Blautia_A sp900551465"          "s__Faecalibacterium prausnitzii_I" "s__Klebsiella pneumoniae"          "s__Streptococcus vestibularis" 


consistent_species_52perc <- names(which(table(unlist(sig_control_list_2)) >= 12))
consistent_species_91perc <- names(which(table(unlist(sig_control_list_2)) >= 21))

#-------------------------------------------------------------------------------
#see which of these species are also associated to Elixhauser general and its components

#general elix is: elix_diff_Q, elix_diff_R2, 
#components is: Q_sign_df_comorbidities, R2_sign_df_comorbidities
 
table(consistent_species_52perc %in% names(elix_diff_Q))
#FALSE  TRUE 
#12   237 
table(consistent_species_91perc %in% names(elix_diff_Q))
#TRUE 
#33

consistent_species_52perc_sub <- consistent_species_52perc[consistent_species_52perc %in% names(elix_diff_Q)]
#consistent_species_91perc has all

#41 to general elix score
table(elix_diff_Q[consistent_species_52perc_sub] <= 0.1 & elix_diff_R2[consistent_species_52perc_sub] > 0)
#FALSE  TRUE 
#196    41
#get the names
sign_elix_general_52_perc_names <- names(which(elix_diff_Q[consistent_species_52perc_sub] <= 0.1 & elix_diff_R2[consistent_species_52perc_sub] > 0))


#14 species
table(elix_diff_Q[consistent_species_91perc] <= 0.1 & elix_diff_R2[consistent_species_91perc] > 0)
#FALSE  TRUE 
#19    14 
sign_elix_general_91_perc_names <- names(which(elix_diff_Q[consistent_species_91perc] <= 0.1 & elix_diff_R2[consistent_species_91perc] > 0))


#
number_of_sign_elix <- rowSums((Q_sign_df_comorbidities[consistent_species_52perc_sub,] <= 0.1 & R2_sign_df_comorbidities[consistent_species_52perc_sub,] > 0) == TRUE)
Q_sign_df_comorbidities["s__Streptococcus sp001556435",] # has 3, Elixhauser_Anemia, Elixhauser_Liver, Elixhauser_PUD

table(number_of_sign_elix > 0)
#FALSE  TRUE 
#179    58
sign_elix_components_52_perc_names <- names(which(number_of_sign_elix > 0))


number_of_sign_elix <- rowSums((Q_sign_df_comorbidities[consistent_species_91perc,] <= 0.1 & R2_sign_df_comorbidities[consistent_species_91perc,] > 0) == TRUE)
table(number_of_sign_elix > 0)
#FALSE  TRUE 
#18     15
sign_elix_components_91_perc_names <- names(which(number_of_sign_elix > 0))

#-----------------------------------
#look at the numbers
length(consistent_species_52perc); length(sign_elix_general_52_perc_names); length(sign_elix_components_52_perc_names); length(unique(c(sign_elix_general_52_perc_names, sign_elix_components_52_perc_names)))
# 249 / 237 that are matched..
# 41
# 58
# 60
60 / 249 #= 24.1%


length(consistent_species_91perc); length(sign_elix_general_91_perc_names); length(sign_elix_components_91_perc_names); length(unique(c(sign_elix_general_91_perc_names, sign_elix_components_91_perc_names)))
# 33
# 14
# 15
# 16
16 / 33 #48.5%



#-----------------------------------
#21 / 23 as cutoff
data_list <- list()
for (i in 1:length(canX_Ctr_list)) {
  control <- canX_Ctr_list[[i]]
  control_sub <- control[names(which(table(unlist(sig_control_list_2)) >= 21)),c("P", "R2", "Q")]
  data_list[[i]] <- control_sub
}

R2_df <- do.call(cbind, lapply(data_list, function(x) x$R2))
Q_df <- do.call(cbind, lapply(data_list, function(x) x$Q))

colnames(R2_df) <- names(canX_Ctr_list)
colnames(Q_df) <- names(canX_Ctr_list)


rownames(R2_df) <- names(which(table(unlist(sig_control_list_2)) >= 21))
rownames(Q_df) <- names(which(table(unlist(sig_control_list_2)) >= 21))

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

col <- c(brewer.pal(12,'Set3'),brewer.pal(12,'Paired'),brewer.pal(9,'Set1'),brewer.pal(8,'Set2'))[1:length(unique(label[,level]))]
names(col) <- unique(label[,level])

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

cancer_class_names <- c("prostate", "bronchus lung", "liver intrahepatic bile ducts", "colon", "melanoma of skin", 
                        "bladder", "breast", "rectum", "brain", "unspecified skin", "kidney except renal pelvis", "esophagus",
                        "oropharynx", "neuroendocrine tumors", "multiple myeloma plasma cell", "pancreas", "tonsil",
                        "corpus uteri", "other connective soft tissue", "lymphoid leukemia", "ovary", "stomach", "base of tongue" )

rownames(R2.sig)[rownames(R2.sig) == "plasma cell"] <- "multiple myeloma plasma cell"


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
        row_names_max_width = max_text_width(rownames(R2.sig)),
        row_names_gp = gpar(fontsize = 14)
) 

left_margin = 4

packed_legends <- packLegend(lgd_sig_p, direction = "vertical")
heatmap_plot <- draw(heatmap_plot, annotation_legend_list = list(packed_legends), merge_legend=T, 
                     padding = unit(c(0.1, left_margin, 0.1, 0.1), "cm")) 


png(paste0('./Figure/RM_figures/heatmap_consistent_species.png'), width =13*330, height = 10*330, res=300)
heatmap_plot
dev.off()











