
require(tidyverse)
require(RColorBrewer)
library(corrr)
library(corrplot)
library(ComplexHeatmap)
library(circlize)
library(ggraph)
library(igraph)
library(readxl)
library(openxlsx)
require(lmerTest)
#install.packages("emmeans") #not working at time of trying, probably a temporary bug
#remotes::install_github("rvlenth/emmeans", dependencies = TRUE, build_opts = "")
require(emmeans)
require(modelbased)
require(data.table)
require(compositions)

require(NbClust)
require(rstatix)
library(igraph)



#-------------------------------------------------------------------------------

file_dir <- dirname(rstudioapi::getSourceEditorContext()$path)
file_dir <- dirname(dirname(dirname(dirname(dirname(file_dir)))))
setwd(file_dir)


source("./Code/Submission/MayoOncobiomeStudy/Code/fastCCLasso.R")
source("./Code/Submission/MayoOncobiomeStudy/Code/zicoseq_heatmap_from_df_general.R")
source("./Code/Submission/MayoOncobiomeStudy/Code/zicoseq_heatmap_from_df_general_w_bar.R")


#-------------------------------------------------------------------------------


load(file = 'Data/data.obj.raw.core.RData') 
clin_meta <- as.data.frame(data.obj$meta.dat)

#tax_table <- as.data.frame(data.obj.rff$otu.tab)
#tax_table_RA <- sweep(tax_table, MARGIN = 2, colSums(tax_table), '/') * 100


#load the metaphlan table
filename <- "~/Dropbox/Mayo_RS/R/Oncobiome full cohort analysis/data/Metaphlan_taxonomy_MCCM_genus.csv"
tax_table_RA <- read.csv(filename)

rownames(tax_table_RA) <- tax_table_RA$X
tax_table_RA <- tax_table_RA[,2:ncol(tax_table_RA)]

tax_table_RA[1:4, 1:4]


#-------------------------------------------------------------------------------
#get highest prevalent genera

tax_table_g <- tax_table_RA
tax_table_g_RA <- tax_table_g

tax_table_g_RA[1:10,1:10]

cutoff_perc <- ncol(tax_table_g_RA) * 0.4 
prev_genera <- which(rowSums(tax_table_g_RA > 0.05) > cutoff_perc)
length(prev_genera) #71 Kraken genera >0.05% in >40% of the samples, 41 for Metaphlan

names(prev_genera)


#get their prevalence
tax_table_g_RA_prev <- tax_table_g_RA[prev_genera,]
prevalence_genera <- sort(rowSums(tax_table_g_RA_prev > 0.05), decreasing = T)


#get the summed counts
sum(colSums(tax_table_g[prev_genera,])) / sum(colSums(tax_table_g)) * 100
#84.4% off all reads for Kraken
#83.96 for Metaphlan

#clr transform these
tax_table_g_sub <- tax_table_g[names(prev_genera),] #relative abundance


genus_names <- rownames(tax_table_g_sub)
all_genus_names <- rownames(tax_table_g_RA)


length(all_genus_names); length(genus_names)
#[1] 900
#[1] 41


#-------------------------------------------------------------------------------

# Calculate Spearman correlation coefficients
#cor_matrix <- cor(t(tax_table_g_sub), method = "spearman")
#dim(cor_matrix)

# Create a correlation network plot
#t(tax_table_g_sub) %>% 
#  correlate(method = "spearman") %>%
#  network_plot(min_cor = 0.3)


#-------------------------------------------------------------------------------
#Run the correlations using fastCCLasso

dim(tax_table_g_sub)
#41 1651

#this runs 4 methods and takes some time
method_name <- c("fastCCLasso","SparCC","CCLasso","COAT");
method_num <- length(method_name);

set.seed(1212); 

#fileout = "fastCCLasso_output_Metaphlan.csv";
#use half min as pseudocount and scale to between 0 and 1, not 0 and 100%
OTUdata <- t(tax_table_g_RA_prev) / 100  # back to [0,1]
min_nonzero <- min(OTUdata[OTUdata > 0])
pc <- min_nonzero / 2
xMat <- OTUdata + pc

#samples as rows and taxa as columns
xMat[1:5,1:5]

n <- dim(xMat)[1];
p <- dim(xMat)[2];

cor_lower <- matrix(0,p*(p-1)/2,method_num);
colnames(cor_lower) <- method_name;

for(i in 1:method_num){
  result.i <- Callallmethods(method=method_name[i], xMat=xMat,cv_k=3,Edge_eps=0.05);
  cor_lower[,i] <- result.i$est_lower;
}

vnames <- colnames(OTUdata)

lower.ind <- which(lower.tri(matrix(0,p,p)), arr.ind = TRUE)

CorEst <- data.frame(variable_1=vnames[lower.ind[,1]],
                     variable_2=vnames[lower.ind[,2]],
                     cor_lower)


#write.table(CorEst, file=fileout, row.names=F,sep = ",")


#--------------------

par(mfrow=c(2,2))
par(mar=c(2,2,2,2))
hist(CorEst$fastCCLasso, xlim=c(-1,1))
hist(CorEst$SparCC, xlim=c(-1,1))
hist(CorEst$CCLasso, xlim=c(-1,1))
hist(CorEst$COAT, xlim=c(-1,1))


#--------------------
#put in matrix format

dim(OTUdata) #41 genera

CorEst_double <- CorEst[, c("variable_1", "variable_2", "fastCCLasso")]
CorEst_temp <- CorEst[, c("variable_2", "variable_1", "fastCCLasso")]
colnames(CorEst_temp) <- colnames(CorEst_double)

CorEst_double <- rbind(CorEst_double, CorEst_temp)

wide_data <- CorEst_double %>%
  pivot_wider(names_from = variable_2, id_cols = variable_1, names_repair = "minimal", id_expand = T, values_from = fastCCLasso) %>%
  column_to_rownames(var = "variable_1")


wide_data <- wide_data[,order(colnames(wide_data))]

cor_matrix <- as.matrix(wide_data)
#make NAs 1

cor_matrix[is.na(cor_matrix)] <- 1

correlation_sub_fastCCLasso <- cor_matrix


#-------------------------------------------------------------------------------

par(mfrow=c(1,1))
par(mar=c(2,2,2,2))
# Create a correlation heatmap
corrplot(cor_matrix, method = "color", type = "upper", 
         tl.cex = 0.8, tl.col = "black", tl.srt = 45, 
         title = "Genus fastCCLasso Correlation Heatmap",
         order = "hclust", hclust.method = "ward.D2")


#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#heatmap plot


set.seed(42)
heatmap_plot_fastCCLasso <- Heatmap(cor_matrix, rect_gp = gpar(type = "none"), column_dend_side = "top",
                        cell_fun = function(j, i, x, y, w, h, fill) {
                          if(as.numeric(x) <= 1  + as.numeric(y)) { #changing + to - shows only one diagnonal
                            grid.rect(x, y, w, h, gp = gpar(fill = fill, col = fill))
                          }
                        }, heatmap_legend_param = list(direction = "vertical"),
                        col = colorRamp2(c(min(cor_matrix), 0, max(cor_matrix)), 
                                         c("#0072B2",'white', "#CC79A7")),#
                        column_gap = unit(1, "mm"), 
                        border = TRUE,
                        column_title_rot = 40,
                        column_names_rot = 90,
                        column_names_max_height = unit(10, "cm"),
                        show_column_names = T, show_row_names = T,
                        show_column_dend = F,
                        show_row_dend = T,
                        cluster_rows = T,
                        cluster_columns = F,
                        row_names_gp = gpar(fontface = 'italic'),
                        row_names_max_width = max_text_width(rownames(cor_matrix), gp = gpar(fontsize = 18)),
                        row_km = 5, row_km_repeats = 100) 

heatmap_plot_fastCCLasso <- draw(heatmap_plot_fastCCLasso)
heatmap_plot_fastCCLasso

cluster_list <- ComplexHeatmap::row_order(heatmap_plot_fastCCLasso)
names(cluster_list) <- c("Bacteroidales", "Faecalibacterium\nDorea", "misc1", 
                         "misc2", "Blautia\nStreptoccus")


sel_colors <- brewer.pal(n = 8, "Dark2")[c(3,2,1,4:8)]


set.seed(42)
heatmap_plot_fastCCLasso <- Heatmap(cor_matrix, rect_gp = gpar(type = "none"), column_dend_side = "top",
                                cell_fun = function(j, i, x, y, w, h, fill) {
                                  if(as.numeric(x) <= 1  + as.numeric(y)) { #changing + to - shows only one diagnonal
                                    grid.rect(x, y, w, h, gp = gpar(fill = fill, col = fill))
                                  }
                                }, heatmap_legend_param = list(direction = "vertical"),
                                col = colorRamp2(c(min(cor_matrix), 0, max(cor_matrix)), 
                                                 c("#0072B2",'white', "#CC79A7")),#
                                column_gap = unit(1, "mm"), 
                                border = TRUE,
                                column_title_rot = 40,
                                column_names_rot = 90,
                                column_names_max_height = unit(10, "cm"),
                                show_column_names = T, show_row_names = T,
                                show_column_dend = F,
                                show_row_dend = T,
                                cluster_rows = T,
                                cluster_columns = F,
                                column_order = rownames(cor_matrix)[unlist(cluster_list)],
                                row_names_gp = gpar(fontface = 'italic'), 
                                column_names_gp = gpar(fontface = 'italic'),
                                row_names_max_width = max_text_width(rownames(cor_matrix), gp = gpar(fontsize = 18)),
                                row_km = 5, row_km_repeats = 100,
                                left_annotation = rowAnnotation(
                                  foo = anno_block(
                                    labels = names(cluster_list),
                                    labels_gp = gpar(fontsize = 11)
                                  )
                                ))
heatmap_plot_fastCCLasso



png("./Figure/RM_figures/heatmap fastCCLasso MCCM metaphlan.png",width = 12.5*300, height=11*300, res=300)
set.seed(42)
heatmap_plot_fastCCLasso
dev.off()


#--------------------------------------

sum_func <- function(vec, feature_clusters=feature_clusters, all_genus_names=all_genus_names) {
  names(vec) <- all_genus_names
  R_sum <- lapply(feature_clusters, function(x) sum(vec[x]))
  return(R_sum)
}

#--------------------------------------

feature_clusters_fastCCLasso <- lapply(cluster_list, function(indices) {
  rownames(cor_matrix)[indices]
})

feature_clusters_fastCCLasso


prev_genera_prevalence <- rowSums(tax_table_g_RA_prev > 0.05) #above 0.05%
names(prev_genera_prevalence) <- genus_names

rank_list <- lapply(feature_clusters_fastCCLasso, function(x) prev_genera_prevalence[x]) #number of samples above 0.1%
rank_list_names <- lapply(rank_list, function(x) names(x))



res_list <- c()
for (i in 1:ncol(tax_table_g_RA)) {
  vec <- tax_table_g_RA[,i]
  res_list[[i]] <- sum_func(vec, feature_clusters_fastCCLasso, all_genus_names)
}

length(res_list) #1651
res_clus_df_fastCCLasso <- as.data.frame(do.call(rbind, res_list))


par(mfrow=c(3,2))
hist(as.numeric(res_clus_df_fastCCLasso[,1]), main=names(feature_clusters_fastCCLasso)[1], xlim=c(0,100), breaks=20)
hist(as.numeric(res_clus_df_fastCCLasso[,2]), main=names(feature_clusters_fastCCLasso)[2], xlim=c(0,100), breaks=20)
hist(as.numeric(res_clus_df_fastCCLasso[,3]), main=names(feature_clusters_fastCCLasso)[3], xlim=c(0,100), breaks=20)
hist(as.numeric(res_clus_df_fastCCLasso[,4]), main=names(feature_clusters_fastCCLasso)[4], xlim=c(0,100), breaks=20)
hist(as.numeric(res_clus_df_fastCCLasso[,5]), main=names(feature_clusters_fastCCLasso)[5], xlim=c(0,100), breaks=20)


res_clus_df_fastCCLasso_2 <- as.data.frame(apply(res_clus_df_fastCCLasso, 2, function(x) as.character(x)))
rownames(res_clus_df_fastCCLasso_2) <- colnames(tax_table_g_RA)

head(res_clus_df_fastCCLasso_2)

#write.csv(res_clus_df_fastCCLasso_2, "./Figure/RM_figures/res_clus_df_fastCCLasso Metaphlan.csv", row.names = T)



#-------------------------------------------------------------------------------
#associate these with phenotypes

head(res_clus_df_fastCCLasso)
dim(res_clus_df_fastCCLasso)

rownames(res_clus_df_fastCCLasso) <- colnames(tax_table_g_RA)
  
#for all of these test differential abundance for every cancer type using linear model

#lmer
clin_meta$icd10_first_3_name[is.na(clin_meta$icd10_first_3_name)] <- "Healthy"
clin_meta$cancer_healthy <- "Cancer"
clin_meta$cancer_healthy[clin_meta$icd10_first_3_name == "Healthy"] <- "Healthy"
  
table(clin_meta$cancer_healthy)
#Cancer Healthy 
#1364     287 

#all in here
clin_meta_sub <- clin_meta

#-------------------------------------------------------------------------------

res_clus_df_fastCCLasso_num <- as.data.frame(apply(res_clus_df_fastCCLasso, 2, as.numeric))
rownames(res_clus_df_fastCCLasso_num) <- rownames(res_clus_df_fastCCLasso)
res_clus_df_fastCCLasso <- res_clus_df_fastCCLasso_num
res_clus_df_fastCCLasso_for_lm <- res_clus_df_fastCCLasso

res_clus_df_fastCCLasso_for_lm$other <- 100 - rowSums(res_clus_df_fastCCLasso_for_lm)


clin_meta_sub$cancer_healthy <- as.factor(clin_meta_sub$cancer_healthy)
clin_meta_sub$cancer_healthy <- relevel(as.factor(clin_meta_sub$cancer_healthy),ref = "Healthy")


res_list_clus_df_fastCCLasso <- list()
res_list_clus_df_fastCCLasso[[1]] <- summary(lm(as.numeric(res_clus_df_fastCCLasso_for_lm[,1]) ~ clin_meta_sub$cancer_healthy))
res_list_clus_df_fastCCLasso[[2]] <- summary(lm(as.numeric(res_clus_df_fastCCLasso_for_lm[,2]) ~ clin_meta_sub$cancer_healthy))
res_list_clus_df_fastCCLasso[[3]] <- summary(lm(as.numeric(res_clus_df_fastCCLasso_for_lm[,3]) ~ clin_meta_sub$cancer_healthy))
res_list_clus_df_fastCCLasso[[4]] <- summary(lm(as.numeric(res_clus_df_fastCCLasso_for_lm[,4]) ~ clin_meta_sub$cancer_healthy))
res_list_clus_df_fastCCLasso[[5]] <- summary(lm(as.numeric(res_clus_df_fastCCLasso_for_lm[,5]) ~ clin_meta_sub$cancer_healthy))
res_list_clus_df_fastCCLasso[[6]] <- summary(lm(as.numeric(res_clus_df_fastCCLasso_for_lm[,6]) ~ clin_meta_sub$cancer_healthy))
names(res_list_clus_df_fastCCLasso) <- names(res_clus_df_fastCCLasso_for_lm)


res_clus_fastCCLasso_df <- as.data.frame(do.call(rbind, lapply(res_list_clus_df_fastCCLasso, function(x) c(as.numeric(x$coefficients[,"Estimate"][2]), as.numeric(x$coefficients[,"Pr(>|t|)"][2])))))
colnames(res_clus_fastCCLasso_df) <- c("Estimate", "p-value")

res_clus_fastCCLasso_df$fdr <- p.adjust(res_clus_fastCCLasso_df$`p-value`, method = "fdr")
res_clus_fastCCLasso_df
#                         Estimate      p-value       fdr
#Bacteroidales           -2.6858978 2.665423e-05 3.730231e-05
#Faecalibacterium\nDorea -7.6975335 9.707480e-28 5.824488e-27
#misc1                   -1.2610019 3.108526e-05 3.730231e-05
#misc2                    0.3099127 6.398930e-01 6.398930e-01
#Blautia\nStreptoccus     6.9572122 1.496044e-14 4.488132e-14
#other                    4.3773083 4.220141e-07 8.440281e-07


#flip the direction
to_plot <- res_clus_fastCCLasso_df$Estimate
Heatmap(to_plot, show_column_names = F, row_labels = rownames(res_clus_fastCCLasso_df), name="R2\nlinear model\nwith cancer",
        col = colorRamp2(c(-(max(to_plot)), 0, max(to_plot)), 
                         c("#0072B2",'white', "#CC79A7")))



#reorganize and number the clades, add p-values, color the boxes

#order of plotting
res_clus_df_fastCCLasso_for_lm_boxplots <- res_clus_df_fastCCLasso_for_lm
colnames(res_clus_df_fastCCLasso_for_lm_boxplots)

col_ordered <- c("Faecalibacterium\nDorea", "Bacteroidales", "misc1", "misc2", "other", "Blautia\nStreptoccus")

res_clus_df_fastCCLasso_for_lm_boxplots <- res_clus_df_fastCCLasso_for_lm_boxplots[,col_ordered]
res_clus_fastCCLasso_df <- res_clus_fastCCLasso_df[col_ordered,]
q.vals <- round(res_clus_fastCCLasso_df$fdr, digits = 5)
q.vals[q.vals < 0.0001] <- "<0.0001"


sel_colors <- RColorBrewer::brewer.pal(8,"Accent")
#sel_colors_plot <- RColorBrewer::brewer.pal(12,"Paired")
sel_colors_plot <- c("#0072B2", "#CC79A7")


cluster_names <- colnames(res_clus_df_fastCCLasso_for_lm_boxplots)

png(paste0("./Figure/RM_figures/boxplots_fastCCLasso_cancer_vs_healthy_metaphlan.png"), width = 6.2*300, height=6*300, res = 300)

par(mfrow=c(2,3))

for (i in 1:length(res_clus_df_fastCCLasso_for_lm_boxplots)) {
  par(mar=c(2.4,4,4.8,1.5))
  boxplot(as.numeric(res_clus_df_fastCCLasso_for_lm_boxplots[,i]) ~ clin_meta_sub$cancer_healthy, las=1, ylab="",
          border="black", col=sel_colors_plot, outcol = sel_colors_plot, outpch=19, cex.axis=1.3)
  
  parusr <- par('usr')
  scale_temp <- 1.2
  segments(1, parusr[4]/scale_temp, 2, parusr[4]/scale_temp) 
  scale_temp <- 1.15
  text(1.5, parusr[4]/scale_temp, paste("q=", q.vals[i], sep=""), cex=1.3, font=1)
  
  mtext(cluster_names[i], line=0.2, cex=1)
  mtext("relative abundance %", side = 2, line=2.5, cex=1)
}

dev.off()


#--------------------------------------
#make the plot only for 3 named co-abundance groups


col_ordered_sub <- c("Blautia\nStreptoccus", "Bacteroidales", "Faecalibacterium\nDorea")

res_clus_df_fastCCLasso_for_lm_boxplots_sub <- res_clus_df_fastCCLasso_for_lm_boxplots[,col_ordered_sub]
res_clus_fastCCLasso_df_sub <- res_clus_fastCCLasso_df[col_ordered_sub,]
q.vals <- round(res_clus_fastCCLasso_df_sub$fdr, digits = 5)
q.vals[q.vals < 0.0001] <- "<0.0001"


sel_colors <- RColorBrewer::brewer.pal(8,"Accent")
#sel_colors_plot <- RColorBrewer::brewer.pal(12,"Paired")
sel_colors_plot <- c("#0072B2", "#CC79A7")

cluster_names <- colnames(res_clus_df_fastCCLasso_for_lm_boxplots_sub)


png(paste0("./Figure/RM_figures/boxplots_fastCCLasso_cancer_vs_healthy_metaphlan_only3.png"), width = 6.2*300, height=4.5*300, res = 300)

par(mfrow=c(1,3))

for (i in 1:length(res_clus_df_fastCCLasso_for_lm_boxplots_sub)) {
  par(mar=c(2.4,4,4.8,1.5))
  boxplot(as.numeric(res_clus_df_fastCCLasso_for_lm_boxplots_sub[,i]) ~ clin_meta_sub$cancer_healthy, las=1, ylab="",
          border="black", col=sel_colors_plot, outcol = sel_colors_plot, outpch=19, cex.axis=1.3, ylim=c(0,100))
  
  parusr <- par('usr')
  scale_temp <- 1.1
  segments(1, parusr[4]/scale_temp, 2, parusr[4]/scale_temp) 
  scale_temp <- 1.07
  text(1.5, parusr[4]/scale_temp, paste("q=", q.vals[i], sep=""), cex=1.1, font=1)
  
  mtext(cluster_names[i], line=0.2, cex=1)
  mtext("relative abundance %", side = 2, line=2.5, cex=1)
}
dev.off()



#-------------------------------------------------------------------------------
#look into subtypes

cutoff_cancer_class_names <- names(table(clin_meta$icd10_first_3_name_short))[which(table(clin_meta$icd10_first_3_name_short) > 15)]
cutoff_cancer_class_names <- cutoff_cancer_class_names[2:length(cutoff_cancer_class_names)] #first one is Healthy, remove


clin_meta_sub$icd10_first_3_name_short <- relevel(as.factor(clin_meta_sub$icd10_first_3_name_short),ref = "healthy")

res_list_clus_df_fastCCLasso_groups <- list()
res_list_clus_df_fastCCLasso_groups[[1]] <- summary(lm(as.numeric(res_clus_df_fastCCLasso_for_lm[,1]) ~ clin_meta_sub$icd10_first_3_name_short))
res_list_clus_df_fastCCLasso_groups[[2]] <- summary(lm(as.numeric(res_clus_df_fastCCLasso_for_lm[,2]) ~ clin_meta_sub$icd10_first_3_name_short))
res_list_clus_df_fastCCLasso_groups[[3]] <- summary(lm(as.numeric(res_clus_df_fastCCLasso_for_lm[,3])~ clin_meta_sub$icd10_first_3_name_short))
res_list_clus_df_fastCCLasso_groups[[4]] <- summary(lm(as.numeric(res_clus_df_fastCCLasso_for_lm[,4]) ~ clin_meta_sub$icd10_first_3_name_short))
res_list_clus_df_fastCCLasso_groups[[5]] <- summary(lm(as.numeric(res_clus_df_fastCCLasso_for_lm[,5]) ~ clin_meta_sub$icd10_first_3_name_short))
res_list_clus_df_fastCCLasso_groups[[6]] <- summary(lm(as.numeric(res_clus_df_fastCCLasso_for_lm[,6]) ~ clin_meta_sub$icd10_first_3_name_short))
names(res_list_clus_df_fastCCLasso_groups) <- names(res_clus_df_fastCCLasso_for_lm)


estimate_list <- lapply(res_list_clus_df_fastCCLasso_groups, function(x) x$coefficients[,"Estimate"])
p_list <- lapply(res_list_clus_df_fastCCLasso_groups, function(x) x$coefficients[,"Pr(>|t|)"])

estimate_df <- do.call(cbind, estimate_list)
p_df <- do.call(cbind, p_list)


#subset only cancer classes with reasonable numbers
rownames(estimate_df) <- gsub("clin_meta_sub\\$icd10_first_3_name_short", "", rownames(estimate_df))
estimate_df <- estimate_df[cutoff_cancer_class_names,]

rownames(p_df) <- gsub("clin_meta_sub\\$icd10_first_3_name_short", "", rownames(p_df))
p_df <- p_df[cutoff_cancer_class_names,]

#redo the fdr
fdr_df <- apply(p_df, 2, function(x) p.adjust(x, method = "fdr"))
head(fdr_df)

par(mfrow=c(2,1))
hist(p_df[,1], breaks=20); hist(fdr_df[,1], breaks=20)

to_plot <- estimate_df
Heatmap(to_plot, show_column_names = T, row_labels = rownames(to_plot),
        name="coeff.\nlinear model\nwith healthy", cluster_rows = T,
        col = colorRamp2(c(min(to_plot), 0, max(to_plot)), 
                         c("#0072B2",'white', "#CC79A7")))

dim(fdr_df)
dim(estimate_df)

#fdr cutoff 0.05
sign_cutoffs <- c(0.1, 0.05, 0.01)
hm <- zicoseq_heatmap(FDR_df = t(fdr_df), R2_df = t(estimate_df), sign_filter=F, sign_cutoffs = sign_cutoffs, 
                      grp.labels = rownames(fdr_df), rotate = T, column_name_rot=40, left_margin=3)
hm$heatmap_plot

png(paste0("./Figure/RM_figures/fastCClasso_assocation_with_cancer_type.png"), width = 6.25*300, height=7*300, res = 300)
hm$heatmap_plot
dev.off()


