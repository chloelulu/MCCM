require(tidyverse)
require(RColorBrewer)
library(gt)
require(ComplexHeatmap)
require(circlize)
require(ggpubr)
library(openxlsx)


file_dir <- dirname(rstudioapi::getSourceEditorContext()$path)
file_dir <- dirname(dirname(dirname(dirname(file_dir))))

setwd(file_dir) #"mforge_clean"

source("Code/Submission/MayoOncobiomeStudy/Code/zicoseq_heatmap_from_df_general.R")


load(file = 'Data/data.obj.pathway.RData') 
path_table <- as.data.frame(data.obj$otu.tab)
clin_meta <- as.data.frame(data.obj$meta.dat)
clin_meta <- clin_meta[colnames(path_table),]
sign_cutoffs <- c(0.1, 0.05, 0.01)
cutoff <- 0.1


#-------------------------------------------------------------------------------
#make pathways
load("Figure/subCancerX_Control_func/pathway/DAA_P_R2.RData")
zico_res_ctr <- list(P.All$pathway, Q.All$pathway, R2.All$pathway)
names(zico_res_ctr) <- c("P", "Q", "R2")

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

length(elix_variables) #22


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
length(elix_variables_diff_list) 


#-------------------------------------------------------------------------------
#get the relevant information out, only the R2, coef, q value

#these variables are corrected for cancer class; does that make sense to do here?
elix_variables_diff_list_Q <- lapply(elix_variables_diff_list, function(x) x$qv.list$pathway[,'Qvalue'])
elix_variables_diff_list_R2 <- lapply(elix_variables_diff_list, function(x) x$R2.list$pathway[,'Func1'])
elix_variables_diff_list_coef <- lapply(elix_variables_diff_list, function(x) x$coef.list$pathway)
coef_list_for_sign <- lapply(elix_variables_diff_list_coef, function(x) x[,ncol(x)])


elix_variables_diff_list_R2_sign <- elix_variables_diff_list_R2
for (i in 1:length(elix_variables_diff_list_R2)) {
  elix_temp <- elix_variables_diff_list_R2[[i]]
  elix_variables_diff_list_R2_sign[[i]] <- elix_temp * sign(coef_list_for_sign[[i]])
}

#subset for pathways in all comparisons
temp_table <- table(unlist(lapply(elix_variables_diff_list_R2_sign, names))) == 22
names_to_include <- names(temp_table[temp_table == TRUE])

elix_variables_diff_list_R2_sign <- lapply(elix_variables_diff_list_R2_sign, function(x) x[names_to_include])
elix_variables_diff_list_Q <- lapply(elix_variables_diff_list_Q, function(x) x[names_to_include])

R2_sign_df_comorbidities <- do.call(cbind, elix_variables_diff_list_R2_sign)
Q_sign_df_comorbidities <- do.call(cbind, elix_variables_diff_list_Q)

dim(R2_sign_df_comorbidities); dim(Q_sign_df_comorbidities)


#-------------------------------------------------------------------------------
#find strongest species linked to weightloss
sig_names <- names(which(Q_sign_df_comorbidities[,"Elixhauser_WeightLoss"] < cutoff))
names(head(sort(Q_sign_df_comorbidities[sig_names,"Elixhauser_WeightLoss"], decreasing = T), 10))

#-------------------------------------------------------------------------------
#also get the general Elixhauser associations

load("./Result/CancerOnly_func/pathway/DAA/Elix_score/Elix_score_ZicoSeq.Rdata")
Elix_res_list <- diff.obj

names(Elix_res_list$R2.list)
elix_diff_Q <- Elix_res_list$qv.list$pathway[,'Qvalue']
elix_diff_R2 <- Elix_res_list$R2.list$pathway[,'Func1']
elix_diff_coef <- Elix_res_list$coef.list$pathway[,"Elix_score"]

elix_diff_R2 <- elix_diff_R2 * sign(elix_diff_coef)
#-------------------------------------
#get number of species FDR < 0.05 with positive and negative effect size for components

table_list <- list()
for (i in 1:ncol(R2_sign_df_comorbidities)) {
  R2_vec_temp <- R2_sign_df_comorbidities[,i]
  table_list[[i]] <- table(R2_vec_temp[Q_sign_df_comorbidities[,i] < cutoff] > 0)
}
names(table_list) <- colnames(R2_sign_df_comorbidities)

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
for (i in 1:length(canX_Ex_list)) {
  rownames(canX_Ex_list[[i]]) <- canX_Ex_list[[i]]$pathway
}

#made sure the names match and are in the same order
canX_Ctr_list <- canX_Ctr_list[order(names(canX_Ctr_list))]


#-------------------------------------------------------------------------------
#only select those pathway significant <0.1 FDR in both control and cancer class comparisons

#trying this discovery level cutoff but mark significance in the heatmap

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
merge_list_sig_Q <- lapply(merge_list_sig, function(x) x[,c("Ctr.Q", "Ex.Q")])
merge_list_sig_R2 <- lapply(merge_list_sig, function(x) x[,c("Ctr.R2", "Ex.R2")])

# Sample data
df <- cbind(unlist(lapply(merge_list_sig_Q, function(x) length(which(x$Ctr.Q < 0.1)))),
            unlist(lapply(merge_list_sig_Q, function(x) length(which(x$Ctr.Q < 0.05)))),
            unlist(lapply(merge_list_sig_Q, function(x) length(which(x$Ctr.Q < 0.01)))))
colnames(df) <- c("n pathways q<0.1", "n pathways q<0.05", "n pathways q<0.01")
df <- as.data.frame(df) %>% rownames_to_column('Name')

gt(df) %>%
  tab_header(title = "Number of functional pathways overlapping vs. healthy and vs. other cancers") %>%
  cols_label(Name="Name", 
             `n pathways q<0.1`="n pathways q<0.1", 
             `n pathways q<0.05`="n pathways q<0.05", 
             `n pathways q<0.01`="n pathways q<0.01") %>%
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
names(merge_list_sig_Q)
#"breast"  "liver intrahepatic bile ducts" "multiple myeloma plasma cell"  "pancreas"

#only 1 pathway; cannot plot this properly

i="breast"
FDR_df <- merge_list_sig_Q_w_comorb[[i]]
R2_df <- merge_list_sig_R2_w_comorb[[i]]
grp.names <- colnames(FDR_df)
grp.labels <- colnames(FDR_df)
grp.labels[1:2] <- c("Cancer class vs. healthy", "Cancer class vs. other cancers")
effect_size_cutoff=1

heatmap_list <- zicoseq_heatmap(FDR_df, R2_df, sign_cutoffs=sign_cutoffs, grp.labels, single_cutoff = FALSE, 
                                sign_filter = FALSE, effect_size_cutoff=1, rotate=F, column_name_rot = 40, left_margin = 8)
heatmap_list$heatmap_plot

#filename <- paste0("./Figure/RM_figures/pathway overlap heatmap ",i," all.png")
#png(filename, width =12.5*330, height = 6*330, res = 330)
#heatmap_list$heatmap_plot
#dev.off()



#--------------------------------------
i="liver intrahepatic bile ducts"

FDR_df <- merge_list_sig_Q_w_comorb[[i]]
R2_df <- merge_list_sig_R2_w_comorb[[i]]

heatmap_list <- zicoseq_heatmap(FDR_df, R2_df, sign_cutoffs=sign_cutoffs, grp.labels, single_cutoff = FALSE, 
                                sign_filter = FALSE, effect_size_cutoff=1, rotate=TRUE)
heatmap_list$heatmap_plot

#filename <- paste0("./Figure/RM_figures/pathway overlap heatmap ",i," all.png")
#png(filename, width = 30*330, height = 13*330, res = 330)
#heatmap_list$heatmap_plot
#dev.off()


#only top 25%
heatmap_list <- zicoseq_heatmap(FDR_df, R2_df, sign_cutoffs=sign_cutoffs, grp.labels, single_cutoff = FALSE, 
                                sign_filter = FALSE, effect_size_cutoff=0.25, rotate=TRUE)

#filename <- paste0("./Figure/RM_figures/pathway overlap heatmap ",i," 25perc.png")
#png(filename, width = 16*330, height = 13*330, res = 330)
#heatmap_list$heatmap_plot
#dev.off()



#--------------------------------------
i="multiple myeloma plasma cell"
FDR_df <- merge_list_sig_Q_w_comorb[[i]]
R2_df <- merge_list_sig_R2_w_comorb[[i]]

heatmap_list <- zicoseq_heatmap(FDR_df, R2_df, sign_cutoffs=sign_cutoffs, grp.labels, single_cutoff = FALSE, 
                                sign_filter = FALSE, effect_size_cutoff=1, rotate=TRUE)
heatmap_list$heatmap_plot

#filename <- paste0("./Figure/RM_figures/pathway overlap heatmap ",i," all.png")
#png(filename, width = 16*330, height = 13*330, res = 330)
#heatmap_list$heatmap_plot
#dev.off()



#--------------------------------------
i="pancreas"
FDR_df <- merge_list_sig_Q_w_comorb[[i]]
R2_df <- merge_list_sig_R2_w_comorb[[i]]

heatmap_list <- zicoseq_heatmap(FDR_df, R2_df, sign_cutoffs=sign_cutoffs, grp.labels, single_cutoff = FALSE, 
                                sign_filter = FALSE, effect_size_cutoff=1, rotate=TRUE)
heatmap_list$heatmap_plot

#filename <- paste0("./Figure/RM_figures/pathway overlap heatmap ",i," all.png")
#png(filename, width = 4.5*330, height = 9*330, res = 330)
#heatmap_list$heatmap_plot
#dev.off()



#-------------------------------------------------------------------------------
#cancer classes of interest

cancer_class_w_sign <- names(merge_list_sig_Q)
#"breast"                        "liver intrahepatic bile ducts" "multiple myeloma plasma cell"  "pancreas"

cancer_class_w_sign_row_list <- sapply(cancer_class_w_sign, function(x) which(clin_meta$icd10_first_3_name_short == x))

#tabulate all Elixhauser components with a significant species
test_df <- lapply(merge_list_sig_Q_w_comorb, function(x) (x[,3:ncol(x)] < cutoff)) #0.1 is FDR cutoff

#general one
Elix_components_w_sig <- unique(unlist(lapply(test_df, function(x) colnames(x)[colSums(x == TRUE) > 0])))

#make Elixhauser components specific to that cancer class
Elix_components_w_sig_list <- lapply(test_df, function(x) colnames(x)[colSums(x == TRUE) > 0])


component_list <- list()
cancer_class_w_sign_row_list_sub <- cancer_class_w_sign_row_list[lapply(Elix_components_w_sig_list, length) > 0] 
Elix_components_w_sig_list_sub <- Elix_components_w_sig_list[lapply(Elix_components_w_sig_list, length) > 0]


names(cancer_class_w_sign_row_list_sub)


#-------------------------------------------------------------------------------
#for cancer groups that have only very few patients with these comorbidities these comorbidities do not have to be plotted or considered
#only consider more than >15%

#table(clin_meta$Elixhauser_DMcx[clin_meta$icd10_first_3_name_short == "liver intrahepatic bile ducts"])

selected_elix <- list()
for(i in names(cancer_class_w_sign_row_list_sub)){
  temp_table <- apply(clin_meta[cancer_class_w_sign_row_list_sub[[i]], Elix_components_w_sig_list_sub[[i]], drop =F], 2, table)
  print(temp_table)
  selected_elix[[i]] <- names(which(cbind(temp_table / length(cancer_class_w_sign_row_list_sub[[i]]) * 100)['Yes',] >15))
}


#-------------------------------------------------------------------------------
#plot using the subset of the Elixhauser categories
Elix_components_w_sig_list_sub_cutoff <- Elix_components_w_sig_list_sub
for(i in names(cancer_class_w_sign_row_list_sub)){
  Elix_components_w_sig_list_sub_cutoff[[i]] <- Elix_components_w_sig_list_sub[[i]][Elix_components_w_sig_list_sub[[i]] %in% 
                                                                                      selected_elix[[i]]]
}


#--------------------------------------

sign_lists_to_plot <- list()
R2_lists_to_plot <- list()


#--------------------------------------
i="breast"  # doesn't have any cormobilities but only 1 pathway so does not plot
sig_components <- Elix_components_w_sig_list_sub_cutoff[[i]]
FDR_df <- merge_list_sig_Q_w_comorb[[i]][,c("Ctr.Q","Ex.Q",sig_components)]
R2_df <- merge_list_sig_R2_w_comorb[[i]][,c("Ctr.R2","Ex.R2",sig_components)]

#a single pathway independent of comorbidities; no comorbidities for breast cancer
#down in breast cancer
#                                                         Ctr.Q       Ex.Q
#PWY-7400: L-arginine biosynthesis IV (archaebacteria) 0.01244488 0.05041405
#                                                         Ctr.R2         Ex.R2
#PWY-7400: L-arginine biosynthesis IV (archaebacteria) -0.0007804459 -0.0006695848

grp.names <- colnames(FDR_df)
grp.labels <- colnames(FDR_df)
grp.labels[1:2] <- c("Cancer class vs. healthy", "Cancer class vs. other cancers")

heatmap_list <- zicoseq_heatmap(FDR_df, R2_df, sign_cutoffs=sign_cutoffs, grp.labels, single_cutoff = FALSE, 
                                sign_filter = FALSE, effect_size_cutoff=1, rotate=F, column_name_rot = 40, left_margin = 4)
heatmap_list$heatmap_plot


#filename <- paste0("./Figure/RM_figures/pathway overlap heatmap ",i,"_sig_Elix.png")
#png(filename, width = 8.25*330, height = 4.5*330, res = 330)
#heatmap_list$heatmap_plot
#dev.off()

colnames(FDR_df) <- grp.labels
sign_lists_to_plot[[i]] <- FDR_df
R2_lists_to_plot[[i]] <- R2_df


#--------------------------------------
i="liver intrahepatic bile ducts"
sig_components <- Elix_components_w_sig_list_sub_cutoff[[i]]
FDR_df <- merge_list_sig_Q_w_comorb[[i]][,c("Ctr.Q","Ex.Q",sig_components)]
R2_df <- merge_list_sig_R2_w_comorb[[i]][,c("Ctr.R2","Ex.R2",sig_components)]

grp.names <- colnames(FDR_df)
grp.labels <- colnames(FDR_df)
grp.labels[1:2] <- c("Cancer class vs. healthy", "Cancer class vs. other cancers")

heatmap_list <- zicoseq_heatmap(FDR_df, R2_df, sign_cutoffs=sign_cutoffs, grp.labels, single_cutoff = FALSE, 
                                sign_filter = FALSE, effect_size_cutoff=1, rotate=F, column_name_rot = 40, left_margin = 4)
heatmap_list$heatmap_plot


filename <- paste0("./Figure/RM_figures/pathway overlap heatmap ",i,"_sig_Elix.png")
png(filename, width = 13*330, height = 32*330, res = 330)
heatmap_list$heatmap_plot
dev.off()

colnames(FDR_df) <- grp.labels
sign_lists_to_plot[[i]] <- FDR_df
R2_lists_to_plot[[i]] <- R2_df


heatmap_list <- zicoseq_heatmap(FDR_df, R2_df, sign_cutoffs=sign_cutoffs, grp.labels, single_cutoff = FALSE, 
                                sign_filter = FALSE, effect_size_cutoff=0.25, rotate=F, column_name_rot = 40, left_margin = 4)

filename <- paste0("./Figure/RM_figures/pathway overlap heatmap ",i,"_sig_Elix_25.png")
png(filename, width = 12*330, height = 9*330, res = 330)
heatmap_list$heatmap_plot
dev.off()


#filter to get the pathways independent of comorbidities



#--------------------------------------
i="multiple myeloma plasma cell"

sig_components <- Elix_components_w_sig_list_sub_cutoff[[i]]
FDR_df <- merge_list_sig_Q_w_comorb[[i]][,c("Ctr.Q","Ex.Q",sig_components)]
R2_df <- merge_list_sig_R2_w_comorb[[i]][,c("Ctr.R2","Ex.R2",sig_components)]

grp.names <- colnames(FDR_df)
grp.labels <- colnames(FDR_df)
grp.labels[1:2] <- c("Cancer class vs. healthy", "Cancer class vs. other cancers")

heatmap_list <- zicoseq_heatmap(FDR_df, R2_df, sign_cutoffs=sign_cutoffs, grp.labels, single_cutoff = FALSE, 
                                sign_filter = FALSE, effect_size_cutoff=1, rotate=F, column_name_rot = 40, left_margin = 4)
heatmap_list$heatmap_plot


filename <- paste0("./Figure/RM_figures/pathway overlap heatmap ",i,"_sig_Elix.png")
png(filename, width = 11.5*330, height = 12*330, res = 330)
heatmap_list$heatmap_plot
dev.off()


colnames(FDR_df) <- grp.labels
sign_lists_to_plot[[i]] <- FDR_df
R2_lists_to_plot[[i]] <- R2_df


#--------------------------------------
i="pancreas" # doesn't have any cormobilities, 2 pathways independent of comorbidities
sig_components <- Elix_components_w_sig_list_sub_cutoff[[i]]
FDR_df <- merge_list_sig_Q_w_comorb[[i]][,c("Ctr.Q","Ex.Q",sig_components)]
R2_df <- merge_list_sig_R2_w_comorb[[i]][,c("Ctr.R2","Ex.R2",sig_components)]

grp.names <- colnames(FDR_df)
grp.labels <- colnames(FDR_df)
grp.labels[1:2] <- c("Cancer class vs. healthy", "Cancer class vs. other cancers")

heatmap_list <- zicoseq_heatmap(FDR_df, R2_df, sign_cutoffs=sign_cutoffs, grp.labels, single_cutoff = FALSE, 
                                sign_filter = FALSE, effect_size_cutoff=1, rotate=F, column_name_rot = 40, left_margin = 4)
heatmap_list$heatmap_plot


filename <- paste0("./Figure/RM_figures/pathway overlap heatmap ",i,"_sig_Elix.png")
png(filename, width = 9.5*330, height = 2.5*330, res = 330)
heatmap_list$heatmap_plot
dev.off()


colnames(FDR_df) <- grp.labels
sign_lists_to_plot[[i]] <- FDR_df
R2_lists_to_plot[[i]] <- R2_df


#-------------------------------------------------------------------------------
#write to file pathways

lapply(sign_lists_to_plot, dim)
lapply(R2_lists_to_plot, dim)

sign_lists_to_plot_copy <- sign_lists_to_plot
R2_lists_to_plot_copy <- R2_lists_to_plot

names(sign_lists_to_plot_copy) <- paste0("p_q_", names(sign_lists_to_plot_copy))
names(R2_lists_to_plot_copy) <- paste0("p_R2_", names(R2_lists_to_plot_copy))


#combine the files
lists_to_write_to_supplement <- c(sign_lists_to_plot_copy[1], R2_lists_to_plot_copy[1],
                                  sign_lists_to_plot_copy[2], R2_lists_to_plot_copy[2],
                                  sign_lists_to_plot_copy[3], R2_lists_to_plot_copy[3],
                                  sign_lists_to_plot_copy[4], R2_lists_to_plot_copy[4])

names(lists_to_write_to_supplement) <- gsub("intrahepatic", "intrah", names(lists_to_write_to_supplement))
names(lists_to_write_to_supplement) <- gsub("cell", "c", names(lists_to_write_to_supplement))
nchar(names(lists_to_write_to_supplement))

#make column names consistent
for (i in 1:length(lists_to_write_to_supplement)) {
  colnames(lists_to_write_to_supplement[[i]])[1:2] <- c("Cancer class vs. healthy",	"Cancer class vs. other cancers")
}
#lapply(lists_to_write_to_supplement, function(x) colnames(x))


filename <- "Code/Submission/Supplementary tables/pathways_for_Supplementary Table 4.xlsx"

list_of_datasets <- lists_to_write_to_supplement
write.xlsx(list_of_datasets, file = filename, append=T, rowNames=T)

#getSheetNames('Code/Submission/Supplementary tables/Supplementary Table 4.xlsx')


#-------------------------------------------------------------------------------
#plot pathway abundances stratified by comorbidity

#boxplot; healthy control, all cancer, all cancer + comorbid, cancer class, cancer class + comorbid

#update cancer_row_list in same order as sign_lists_to_plot
names(sign_lists_to_plot)
group_vec <- clin_meta$icd10_first_3_name_short
cancer_row_list <- lapply(names(sign_lists_to_plot), function(name) grep(name, group_vec))


elix_row_names_list <- lapply(sign_lists_to_plot, function(x) colnames(x)[3:ncol(x)])
healthy_rows <- which(group_vec == "healthy")

sel_cols <- brewer.pal(5, "Dark2")

names(sign_lists_to_plot)
for (i in 1:length(sign_lists_to_plot)) { #cancer class
  sig_df <- sign_lists_to_plot[[i]]
  elix_row_names <- elix_row_names_list[[i]]
  cancer_group_rows <- cancer_row_list[[i]]
  all_other_cancer_rows <- which(!1:ncol(path_table) %in% c(healthy_rows, cancer_group_rows))
  cancer_class <- names(sign_lists_to_plot)[i]
  if(sum(grepl('Elixhauser',elix_row_names))>0){ ## breast and pancreas does not have Comorbidity
    pdf(paste0("./Figure/RM_figures/Comorbidity pathways_", cancer_class, ".pdf"), width=5, height=4.5)
    for (j in 1:length(elix_row_names)) { #elix component
      sig_pathways <- rownames(sig_df)[which(sig_df[elix_row_names[j]] < cutoff)]
      
      for (k in 1:length(sig_pathways)) {
        
        temp_ra <- as.numeric(path_table[sig_pathways[k],])
        temp_elix_rows <- which(clin_meta[,elix_row_names[j]] == "Yes")
        
        #tax_table, setdiff indicates not having the comorbidity, intersect indicates having the comorbidity
        par(mar=c(7,4,2,1))
        boxplot(temp_ra[healthy_rows], temp_ra[setdiff(all_other_cancer_rows, temp_elix_rows)],
                temp_ra[intersect(all_other_cancer_rows, temp_elix_rows)], 
                temp_ra[setdiff(cancer_group_rows, temp_elix_rows)], 
                temp_ra[intersect(cancer_group_rows, temp_elix_rows)], las=1, col = sel_cols, varwidth = F, 
                xaxt = "n", outline = F, ylim=range(temp_ra))
        axis(1,  at=1:5, labels = c("Healthy", "other cancer\nno comorbidity", "other cancer\ncomorbidity", 
                                    "cancer class\nno comorbidity", "cancer class\ncomorbidity"), las=2)
        mtext("pathway abundance", 2, line = 2.5)
        mtext(paste0(cancer_class, " _ ", elix_row_names[j], "\n", sig_pathways[k]), 3, line = 0, cex=1)
        
      }
    }
    dev.off()
    
  }
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
  if(sum(grepl('Elixhauser',elix_row_names))>0){ ## breast and pancreas does not have Comorbidity
    pdf(paste0("./Figure/RM_figures/Comorbidity pathways_split_", cancer_class, ".pdf"), width=5.5, height=4.5)
    for (j in 1:length(elix_row_names)) { #elix component
      sig_pathways <- rownames(sig_df)[which(sig_df[elix_row_names[j]] < cutoff)]
      
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
        
        ymax <- quantile(rbind(data, data_all_other, data_cc_other)[,1], na.rm = T, probs = c(1))
        
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
}



#-------------------------------------------------------------------------------
#make plots for all interesting pathways as well without factoring in comorbidities


cancer_row_list <- lapply(names(sign_lists_to_plot), function(name) grep(name, group_vec))
names(cancer_row_list) <- names(sign_lists_to_plot)

names(merge_list_sig_Q_w_comorb)

for (i in 1:length(merge_list_sig_Q_w_comorb)) { #cancer class
  sig_df <- merge_list_sig_Q_w_comorb[[i]]
  
  cancer_group_rows <- cancer_row_list[[i]]
  all_other_cancer_rows <- which(!1:nrow(path_table) %in% c(healthy_rows, cancer_group_rows))
  cancer_class <- names(merge_list_sig_Q_w_comorb)[i]
  
  sig_pathways <- rownames(sig_df)
  
  pdf(paste0("./Figure/RM_figures/Significant pathways_scatter_", cancer_class, ".pdf"), width=3.5, height=4.5)
  for (k in 1:length(sig_pathways)) {
    
    pathway_name <- sig_pathways[k]
    colon_split <- strsplit(pathway_name, ":", fixed = TRUE)
    # Process each element
    pathway_name_new <- sapply(colon_split, function(parts) {
      # Handle cases: with or without colon
      if (length(parts) < 2) {
        # No colon: just wrap whole string
        gsub("(.{30})", "\\1\n", parts[1])
      } else {
        # With colon: add newline after colon, then wrap the second part
        paste0(parts[1], ":\n", gsub("(.{30})", "\\1\n", parts[2]))
      }
    })

    temp_ra <- as.numeric(path_table[sig_pathways[k],])
    plot_df <- as.data.frame(cbind(c(rep("healthy", times=length(healthy_rows)), rep("all other cancer", times=length(all_other_cancer_rows)), rep("cancer_group", times=length(cancer_group_rows))),
                                   c(temp_ra[healthy_rows], temp_ra[all_other_cancer_rows], temp_ra[cancer_group_rows])))
    names(plot_df) <- c("group", "ra")
    plot_df$ra <- as.numeric(plot_df$ra)
    
    p_plot <-  ggboxplot(plot_df, "group", "ra", width = 0.8,
                         color = "group", palette=sel_cols,
                         add = "jitter", ylab= "sqrt(counts)") +
      scale_y_continuous(trans='sqrt') +
      ggtitle(pathway_name_new) +
      theme(legend.position="none") +
      scale_x_discrete(labels = c("healthy", "other\ncancers", "cancer\nclass"))
    
    print(p_plot)
  }
  dev.off()
}



#-------------------------------------------------------------------------------
#For all cancer classes and a single significance level <0.1. Show number of pathways vs control, vs other cancers, and overlap.
#And how many after filtering the cancer-specific comorbidities


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

head(merge_list$breast)
#add number after removing species-specific comorbidities
n_wo_sig_elix <- rep(0, times= length(names(merge_list)))
names(n_wo_sig_elix) <- names(merge_list)


#--------------------------------------
names(sign_lists_to_plot)

lapply(sign_lists_to_plot, head)

length(which(rowSums(sign_lists_to_plot$`liver intrahepatic bile ducts`[,3:ncol(sign_lists_to_plot$`liver intrahepatic bile ducts`)] < 0.1) == 0)) #46
length(which(rowSums(sign_lists_to_plot$`multiple myeloma plasma cell`[,3:ncol(sign_lists_to_plot$`multiple myeloma plasma cell`)] < 0.1) == 0)) #41

n_wo_sig_elix["breast"] <- 1
n_wo_sig_elix["liver intrahepatic bile ducts"] <- 46
n_wo_sig_elix["multiple myeloma plasma cell"] <- 41
n_wo_sig_elix["pancreas"] <- 2


plot_df <- cbind(names(merge_list), 
                 as.character(lapply(sig_control_list, length)), 
                 as.character(lapply(sig_cancer_list, length)),
                 as.character(lapply(merge_list, nrow)),
                 as.numeric(n_wo_sig_elix)
)

plot_df[plot_df == "NULL"] <- 0
plot_df


colnames(plot_df) <- c("Cancer class", "vs. healthy", "vs. other cancers", "overlap", "overlap w/o comorbidities")
plot_df <- as.data.frame(plot_df)

plot_df$`vs. healthy` <- as.numeric(plot_df$`vs. healthy`)
plot_df$`vs. other cancers` <- as.numeric(plot_df$`vs. other cancers`)
plot_df$overlap <- as.numeric(plot_df$overlap)
plot_df$`overlap w/o comorbidities` <- as.numeric(plot_df$`overlap w/o comorbidities`)

plot_df <- plot_df[order(plot_df$overlap, decreasing = T),]

gt(plot_df) %>%
  tab_header(title = "Number of pathways when correcting for comorbidities") %>%
  opt_vertical_padding(scale = 0.1)



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
12 /length(sig_control_list_2) #54.5%
length(which(table(unlist(sig_control_list_2)) >= 12))
#107 pathways

cutoff <- 20
cutoff /length(sig_control_list_2) #90.9%
length(which(table(unlist(sig_control_list_2)) >= cutoff))
#13 pathways; make heatmap of these pathways
names(which(table(unlist(sig_control_list_2)) >= cutoff))


#[1] "CITRULBIO-PWY: L-citrulline biosynthesis"                                        
#[2] "GALACTARDEG-PWY: D-galactarate degradation I"                                    
#[3] "GLUCARDEG-PWY: D-glucarate degradation I"                                        
#[4] "GLUCARGALACTSUPER-PWY: superpathway of D-glucarate and D-galactarate degradation"
#[5] "NAD-BIOSYNTHESIS-II: NAD salvage pathway III (to nicotinamide riboside)"         
#[6] "NONMEVIPP-PWY: methylerythritol phosphate pathway I"                             
#[7] "PANTO-PWY: phosphopantothenate biosynthesis I"                                   
#[8] "PANTOSYN-PWY: superpathway of coenzyme A biosynthesis I (bacteria)"              
#[9] "PWY-4041: &gamma;-glutamyl cycle"                                                
#[10] "PWY-5030: L-histidine degradation III"                                           
#[11] "PWY-6700: queuosine biosynthesis I (de novo)"                                    
#[12] "PWY-7117: C4 photosynthetic carbon assimilation cycle, PEPCK type"               
#[13] "PWY-7456: &beta;-(1,4)-mannan degradation"  



#20 / 22 as cutoff
data_list <- list()
for (i in 1:length(canX_Ctr_list)) {
  control <- canX_Ctr_list[[i]]
  control_sub <- control[names(which(table(unlist(sig_control_list_2)) >= cutoff)),c("P", "R2", "Q")]
  data_list[[i]] <- control_sub
}

R2_df <- do.call(cbind, lapply(data_list, function(x) x$R2))
Q_df <- do.call(cbind, lapply(data_list, function(x) x$Q))

colnames(R2_df) <- names(canX_Ctr_list)
colnames(Q_df) <- names(canX_Ctr_list)

rownames(R2_df) <- names(which(table(unlist(sig_control_list_2)) >= cutoff))
rownames(Q_df) <- names(which(table(unlist(sig_control_list_2)) >= cutoff))

R2_df[is.na(R2_df)] <- 0
Q_df[is.na(Q_df)] <- 1

dim(R2_df)
# 13 22


#------------------------------------
#make heatmap out of this

dim(R2_df); dim(Q_df)
colnames(R2_df); rownames(R2_df)

Q.sig <- t(Q_df)
R2.sig <- t(R2_df)


sign_cutoffs <- c(0.1, 0.05, 0.01)
name_vec = colnames(FDR_df)
max.cutoff <- max(sign_cutoffs)
q.cut1 <- min(sign_cutoffs); q.cut2 <- median(sign_cutoffs); q.cut3 <- max.cutoff

grid.text_1 <- '***'
grid.text_2 <- '**'
grid.text_3 <- '*'

scale_max <- max(c(abs(min(R2.sig)), max(R2.sig)))

lgd_sig_p = Legend(pch = c("*","**","***"), type = "points", labels = paste0("<", sign_cutoffs), background="white", title = "q")


heatmap_plot <- Heatmap(R2.sig, name = "R2 x effect direction", 
        heatmap_legend_param = list(direction = "vertical", title="R2"),
        col = colorRamp2(c(min(R2.sig), 0, max(R2.sig)), c("#0072B2",'white', "#CC79A7")),
        #col = colorRamp2(c(-scale_max, 0, scale_max), c("#0072B2",'white', "#CC79A7")),
        column_gap = unit(1, "mm"), 
        column_title_rot = 30,
        column_names_rot = 40,
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
        cluster_columns = T,
        row_names_max_width = max_text_width(rownames(R2.sig)),
        row_names_gp = gpar(fontsize = 14)
) 



png(paste0('./Figure/RM_figures/heatmap_consistent_pathways.png'), width=11*330, height=10*330, res=300)

draw(heatmap_plot, 
     annotation_legend_list = list(lgd_sig_p),
     annotation_legend_side = "right",
     merge_legend = TRUE,
     padding = unit(c(0.1, 11, 0.1, 0.1), "cm"),
     align_heatmap_legend = "global_center",        
     align_annotation_legend = "global_center")


dev.off()




