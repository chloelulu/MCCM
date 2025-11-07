
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

#source("~/Dropbox/Mayo_RS/R/general functions/fastCCLasso-main/RCode/fastCCLasso.R")

setwd("/Users/M210320/Library/Mobile Documents/com~apple~CloudDocs/mforge_clean/")

source("./Code/Submission/MayoOncobiomeStudy/Code/fastCCLasso.R")
source("./Code/Submission/MayoOncobiomeStudy/Code/zicoseq_heatmap_from_df_general.R")
source("./Code/Submission/MayoOncobiomeStudy/Code/zicoseq_heatmap_from_df_general_w_bar.R")


#-------------------------------------------------------------------------------


load(file = 'Data/data.obj.raw.core.RData') 
clin_meta <- as.data.frame(data.obj$meta.dat)

tax_table <- as.data.frame(data.obj.rff$otu.tab)
tax_table_RA <- sweep(tax_table, MARGIN = 2, colSums(tax_table), '/') * 100


#-------------------------------------------------------------------------------
#get highest prevalent species

dim(tax_table_RA)

tax_table_RA[1:10,1:10]

hist(rowSums(tax_table_RA >= 0.1))

max(rowSums(tax_table_RA >= 0.1))
#1529

cutoff_perc <- ncol(tax_table_RA) * 0.5
prev_species <- which(rowSums(tax_table_RA > 0.05) > cutoff_perc)
length(prev_species)


#-------------------------------------------------------------------------------
#do at the genus level

tax_table_g <- as.data.frame(data.obj.rff$abund.list$Genus)
tax_table_g_RA <- sweep(tax_table_g, MARGIN = 2, colSums(tax_table_g), '/') * 100

tax_table_g_RA[1:10,1:10]

cutoff_perc <- ncol(tax_table_g_RA) * 0.4 
prev_genera <- which(rowSums(tax_table_g_RA > 0.05) > cutoff_perc)
length(prev_genera) #71 genera >0.05% in >40% of the samples

prev_genera


#get their prevalence
tax_table_g_RA_prev <- tax_table_g_RA[prev_genera,]
prevalence_genera <- sort(rowSums(tax_table_g_RA_prev > 0.05), decreasing = T)

#get the summed counts
sum(colSums(tax_table_g[prev_genera,])) / sum(colSums(tax_table_g)) * 100
#84.4% off all reads


#clr transform these
tax_table_g_sub <- tax_table_g[names(prev_genera),]

genus_names <- as.character(sapply(rownames(tax_table_g_sub), function(x) strsplit(x, ";")[[1]][2]))
rownames(tax_table_g_sub) <- genus_names

all_genus_names <- as.character(sapply(rownames(tax_table_g_RA), function(x) strsplit(x, ";")[[1]][2]))


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

#this runs 4 methods and takes some time
method_name <- c("fastCCLasso","SparCC","CCLasso","COAT");
method_num <- length(method_name);

set.seed(1212); 

#fileout = "fastCCLasso_output.csv";

OTUdata <- t(tax_table_g_sub) + 0.5;
total_count <- rowSums(OTUdata);
xMat <- (1/total_count) * OTUdata;  # compositional data


#samples as rows and taxa as columns
xMat[1:5,1:5]
dim(xMat)

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

dim(OTUdata) #71 genera

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

cluster_list <- row_order(heatmap_plot_fastCCLasso)
names(cluster_list) <- c("Evtepia\nAlistipes\nGemmiger", "Faecalibacterium\nDorea", "Blautia\nStreptoccus", 
                         "Anaerostipes\nRoseburia\nBifido", "Bacteroidales\nEisenbergiella")


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
                                    labels_gp = gpar(fontsize = 14)
                                  )
                                ))
heatmap_plot_fastCCLasso



png("./Figure/RM_figures/heatmap fastCCLasso Oncobiome.png",width = 14.5*300, height=12*300, res=300)
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


#name by most prevalent genus
full_genus_names <- names(prevalence_genera)
names(prevalence_genera) <- as.character(sapply(names(prevalence_genera), function(x) strsplit(x, ";")[[1]][2]))
names(full_genus_names) <- names(prevalence_genera)

rank_list <- lapply(feature_clusters_fastCCLasso, function(x) prevalence_genera[x])
rank_list_names <- lapply(rank_list, function(x) full_genus_names[names(x)])

#higher number means found in more samples
cluster_genera_df <- as.data.frame(cbind(unlist(rank_list), as.character(unlist(rank_list_names))))


#inspecting most prevalent genera
#cluster_genera_df$V2[order(as.numeric(cluster_genera_df$V1) / ncol(tax_table_g_sub) * 100)]


  
#$`2` Evtepia_Alistipes_Gemmiger
#[1] "g__Alistipes"      "g__CAG-83"         "g__CAG-110"        "g__Akkermansia"    "g__UBA11774"       "g__UBA1191"        "g__Evtepia"        "g__UBA1417"       
#[9] "g__Adlercreutzia"  "g__Ruminococcus_E" "g__CAG-1427"       "g__SFFH01"         "g__ER4"            "g__Gemmiger"       "g__Acetatifactor" 

#$`1` Faecali_Dorea
#[1] "g__GCA-900066135"    "g__Dorea"            "g__UBA7160"          "g__Coprococcus_A"    "g__Bariatricus"      "g__Faecalibacillus"  "g__CAG-41"          
#[8] "g__Eubacterium_I"    "g__Ruminococcus_A"   "g__Fusicatenibacter" "g__Dorea_A"          "g__Agathobacter"     "g__Lachnospira"      "g__KLE1615"         
#[15] "g__Faecalibacterium"

#$`5` Blautia_Strep
#[1] "g__Erysipelatoclostridium" "g__Blautia"                "g__Faecalimonas"           "g__Ruminococcus_B"         "g__Schaedlerella"          "g__Enterocloster"         
#[7] "g__Clostridium_AQ"         "g__Muricomes"              "g__Clostridium_AP"         "g__Eggerthella"            "g__Streptococcus"          "g__Pauljensenia"          
#[13] "g__Intestinibacter"        "g__Blautia_A"              "g__CAG-317"               

#$`3` Anaerostipes_Roseburia_Bif
#[1] "g__Clostridium_Q"       "g__Bifidobacterium"     "g__Escherichia"         "g__Mediterraneibacter"  "g__Anaerostipes"        "g__Eubacterium_G"      
#[7] "g__Roseburia"           "g__CAG-81"              "g__UBA9502"             "g__Agathobaculum"       "g__NK3B98"              "g__Dialister"          
#[13] "g__Collinsella"         "g__Lachnoclostridium_B" "g__Anaerobutyricum"    

#$`4` Bacteroidales_Eisenbergiella
#[1] "g__Phocaeicola"      "g__Bacteroides"      "g__Parabacteroides"  "g__Dysosmobacter"    "g__Flavonifractor"   "g__Lawsonibacter"    "g__Acutalibacter"   
#[8] "g__Clostridium_A"    "g__Gordonibacter"    "g__Eisenbergiella"   "g__Ruthenibacterium"


res_list <- c()
for (i in 1:ncol(tax_table_g_RA)) {
  vec <- tax_table_g_RA[,i]
  res_list[[i]] <- sum_func(vec, feature_clusters_fastCCLasso, all_genus_names)
}

length(res_list) #1644
res_clus_df_fastCCLasso <- as.data.frame(do.call(rbind, res_list))

names(feature_clusters_fastCCLasso) <- c("Evtepia_Alistipes_Gemmiger", "Faecali_Dorea", 
                                         "Blautia_Strep", "Anaerostipes_Roseburia_Bif", 
                                         "Bacteroidales_Eisenbergiella")

names(res_clus_df_fastCCLasso) <- names(feature_clusters_fastCCLasso)

par(mfrow=c(3,2))
hist(as.numeric(res_clus_df_fastCCLasso[,1]), main=names(feature_clusters_fastCCLasso)[1], xlim=c(0,100), breaks=20)
hist(as.numeric(res_clus_df_fastCCLasso[,2]), main=names(feature_clusters_fastCCLasso)[2], xlim=c(0,100), breaks=20)
hist(as.numeric(res_clus_df_fastCCLasso[,3]), main=names(feature_clusters_fastCCLasso)[3], xlim=c(0,100), breaks=20)
hist(as.numeric(res_clus_df_fastCCLasso[,4]), main=names(feature_clusters_fastCCLasso)[4], xlim=c(0,100), breaks=20)
hist(as.numeric(res_clus_df_fastCCLasso[,5]), main=names(feature_clusters_fastCCLasso)[5], xlim=c(0,100), breaks=20)


res_clus_df_fastCCLasso_2 <- as.data.frame(apply(res_clus_df_fastCCLasso, 2, function(x) as.character(x)))
rownames(res_clus_df_fastCCLasso_2) <- colnames(tax_table_g_RA)

head(res_clus_df_fastCCLasso_2)

write.csv(res_clus_df_fastCCLasso_2, "./Figure/RM_figures/res_clus_df_fastCCLasso.csv", row.names = T)


#-------------------------------------------------------------------------------
#some clusters seem much more correlated then others
#get the correlation metrics summary scores including number of positive and negative


sel_colors <- brewer.pal(n = 8, "Dark2")[c(3,2,1,4:8)]


#make long format out of this and get the summary score plus number of positive, negative above cutoff?
cor_stat_list <- lapply(feature_clusters_fastCCLasso, function(x) cor_matrix[x, x])

lapply(feature_clusters_fastCCLasso, length)
#$Evtepia_Alistipes_Gemmiger
#[1] 15
#$Faecali_Dorea
#[1] 15
#$Blautia_Strep
#[1] 15
#$Anaerostipes_Roseburia_Bif
#[1] 15
#$Bacteroidales_Eisenbergiella
#[1] 11

cor_stat_list_long <- lapply(cor_stat_list, cor_gather)
#remove the correlations that are 1
cor_stat_list_long_filt <- lapply(cor_stat_list_long, function(x) x[!x[,"cor"] == 1,])

#plot of all them in one plot
par(mfrow=c(2,3))
plot(density(cor_stat_list_long_filt[[1]]$cor), xlim=c(-1,1), main=names(cor_stat_list_long_filt)[1], col=sel_colors[1])
plot(density(cor_stat_list_long_filt[[2]]$cor), xlim=c(-1,1), main=names(cor_stat_list_long_filt)[2], col=sel_colors[2])
plot(density(cor_stat_list_long_filt[[3]]$cor), xlim=c(-1,1), main=names(cor_stat_list_long_filt)[3], col=sel_colors[3])
plot(density(cor_stat_list_long_filt[[4]]$cor), xlim=c(-1,1), main=names(cor_stat_list_long_filt)[4], col=sel_colors[4])
plot(density(cor_stat_list_long_filt[[5]]$cor), xlim=c(-1,1), main=names(cor_stat_list_long_filt)[5], col=sel_colors[5])


#combined all in one
plot(density(cor_stat_list_long_filt[[1]]$cor), xlim=c(-0.5,1), ylim=c(0,5), col=sel_colors[1], main="")
lines(density(cor_stat_list_long_filt[[2]]$cor), xlim=c(-0.5,1), col=sel_colors[2])
lines(density(cor_stat_list_long_filt[[3]]$cor), xlim=c(-0.5,1), col=sel_colors[3])
lines(density(cor_stat_list_long_filt[[4]]$cor), xlim=c(-0.5,1), col=sel_colors[4])
lines(density(cor_stat_list_long_filt[[5]]$cor), xlim=c(-0.5,1), col=sel_colors[5])


plot(density(cor_stat_list_long_filt[[1]]$cor), xlim=c(-0.5,1), ylim=c(0,5), col=sel_colors[1], main="correlations per cluster", 
     las=1, lwd=2.5, xlab="Correlation", cex.axis=1.2, cex.lab=1.2)
lines(density(cor_stat_list_long_filt[[2]]$cor), xlim=c(-0.5,1), col=sel_colors[2], lwd=2.5)
lines(density(cor_stat_list_long_filt[[3]]$cor), xlim=c(-0.5,1), col=sel_colors[3], lwd=2.5)
lines(density(cor_stat_list_long_filt[[4]]$cor), xlim=c(-0.5,1), col=sel_colors[4], lwd=2.5)
lines(density(cor_stat_list_long_filt[[5]]$cor), xlim=c(-0.5,1), col=sel_colors[5], lwd=2.5)
ycord <- 5
for(i in 1:length(names(cor_stat_list_long_filt))) {
  text(0.6,ycord, names(cor_stat_list_long_filt)[i], col=sel_colors[i], cex=1.2)
  ycord <- ycord - 0.26
}


#these median correlation coefficients are worth putting in the results section
lapply(cor_stat_list_long_filt, function(x) summary(x$cor))
#$Evtepia_Alistipes_Gemmiger
#Min. 1st Qu.  Median    Mean   3rd Qu.    Max. 
#0.0000  0.1541  0.2210  0.2232  0.2902  0.6084 

#$Faecali_Dorea
#Min. 1st Qu.  Median    Mean   3rd Qu.    Max. 
#0.0000  0.1286  0.1844  0.1821  0.2396  0.4022 

#$Blautia_Strep                               #most correlated and all positive
#Min. 1st Qu.  Median    Mean   3rd Qu.    Max. 
#0.1624  0.3312  0.4276  0.4332  0.5238  0.6958 

#$Anaerostipes_Roseburia_Bif                  #barely correlated, could be called misc
#Min.     1st Qu.   Median     Mean      3rd Qu.  Max. 
#-0.18141  0.00000  0.00000  0.02331  0.07431  0.36604 

#$Bacteroidales_Eisenbergiella                #second most correlated, all positive
#Min. 1st Qu.  Median    Mean   3rd Qu.    Max. 
#0.0000  0.1114  0.2527  0.2725  0.4243  0.7605 


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

clin_meta_sub <- clin_meta[colnames(tax_table_g_RA),]



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

#                                               Estimate      p-value       fdr
#Evtepia_Alistipes_Gemmiger                    -3.8995328 7.481874e-14 1.246979e-13
#Faecali_Dorea                                 -7.1014889 5.943254e-23 1.984944e-22
#Blautia_Strep                                 10.3344340 7.939776e-23 1.984944e-22
#Anaerostipes_Roseburia_Bif                     3.1672373 2.392594e-06 2.990742e-06
#Bacteroidales_Eisenbergiella                  -0.9884386 1.860003e-01 1.860003e-01
#other                                         -1.5122111 4.466501e-02 5.359802e-02


#flip the direction
to_plot <- res_clus_fastCCLasso_df$Estimate
Heatmap(to_plot, show_column_names = F, row_labels = rownames(res_clus_fastCCLasso_df), name="R2\nlinear model\nwith cancer",
        col = colorRamp2(c(-(max(to_plot)), 0, max(to_plot)), 
                         c("#0072B2",'white', "#CC79A7")))



#reorganize and number the clades, add p-values, color the boxes

#order of plotting
res_clus_df_fastCCLasso_for_lm_boxplots <- res_clus_df_fastCCLasso_for_lm
colnames(res_clus_df_fastCCLasso_for_lm_boxplots)

col_ordered <- c("Blautia_Strep", "Anaerostipes_Roseburia_Bif", "Bacteroidales_Eisenbergiella",
                 "other", "Evtepia_Alistipes_Gemmiger", "Faecali_Dorea")
res_clus_df_fastCCLasso_for_lm_boxplots <- res_clus_df_fastCCLasso_for_lm_boxplots[,col_ordered]
res_clus_fastCCLasso_df <- res_clus_fastCCLasso_df[col_ordered,]
q.vals <- round(res_clus_fastCCLasso_df$fdr, digits = 5)
q.vals[q.vals == 0] <- "<0.00001"


sel_colors <- RColorBrewer::brewer.pal(8,"Accent")
#sel_colors_plot <- RColorBrewer::brewer.pal(12,"Paired")
sel_colors_plot <- c("#0072B2", "#CC79A7")



colnames(res_clus_df_fastCCLasso_for_lm_boxplots)

cluster_names <- c("Blautia\nStreptoccus", "Anaerostipes\nRoseburia\nBifido",
                   "Bacteroidales\nEisenbergiella", "Other",
                   "Evtepia\nAlistipes\nGemmiger", "Faecalibacterium\nDorea")


png(paste0("./Figure/RM_figures/boxplots_fastCCLasso_cancer_vs_healthy.png"), width = 6.2*300, height=6*300, res = 300)

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



#-------------------------------------------------------------------------------
#plot correlation network 


dim(correlation_sub_fastCCLasso)

nw_plot <- network_plot(correlation_sub_fastCCLasso, min_cor = 0.4, colours=c("#0072B2",'white', "#CC79A7"))
nw_plot



#-------------------------------------------------------------------------------
#plot overlap between groups

sel_colors <- brewer.pal(n = 8, "Dark2")[c(3,2,1,4:8)]

#sort table by Blautia_Strep
res_clus_df_fastCCLasso_sort <- res_clus_df_fastCCLasso[order(as.numeric(res_clus_df_fastCCLasso$Blautia_Strep), decreasing = T),]

head(res_clus_df_fastCCLasso)
res_clus_df_fastCCLasso_sort <- res_clus_df_fastCCLasso_sort[,c(3, 1:2, 4:5)]


pdf("./Figure/RM_figures/cluster_ra_plots.pdf", width = 6, height = 12)

par(mfrow=c(5,1))
par(mar=c(1,1,1,1))
plot(1:nrow(res_clus_df_fastCCLasso_sort), res_clus_df_fastCCLasso_sort$Blautia_Strep, col=sel_colors[1], lwd=2, axes = F)
points(1:nrow(res_clus_df_fastCCLasso_sort), res_clus_df_fastCCLasso_sort$Evtepia_Alistipes_Gemmiger, col=sel_colors[2], lwd=2)
legend("topright", colnames(res_clus_df_fastCCLasso_sort)[1:2], fill = sel_colors[c(1,2)], border="white", bty = "n")

plot(1:nrow(res_clus_df_fastCCLasso_sort), res_clus_df_fastCCLasso_sort$Blautia_Strep, col=sel_colors[1], lwd=2, axes = F)
points(1:nrow(res_clus_df_fastCCLasso_sort), res_clus_df_fastCCLasso_sort$Faecali_Dorea, col=sel_colors[3], lwd=2)
legend("topright", colnames(res_clus_df_fastCCLasso_sort)[c(1,3)], fill = sel_colors[c(1,3)], border="white", bty = "n")

plot(1:nrow(res_clus_df_fastCCLasso_sort), res_clus_df_fastCCLasso_sort$Blautia_Strep, col=sel_colors[1], lwd=2, axes = F)
points(1:nrow(res_clus_df_fastCCLasso_sort), res_clus_df_fastCCLasso_sort$Anaerostipes_Roseburia_Bif, col=sel_colors[4], lwd=2)
legend("topright", colnames(res_clus_df_fastCCLasso_sort)[c(1,4)], fill = sel_colors[c(1,4)], border="white", bty = "n")

plot(1:nrow(res_clus_df_fastCCLasso_sort), res_clus_df_fastCCLasso_sort$Blautia_Strep, col=sel_colors[1], lwd=2, axes = F)
points(1:nrow(res_clus_df_fastCCLasso_sort), res_clus_df_fastCCLasso_sort$Bacteroidales_Eisenbergiella, col=sel_colors[5], lwd=2)
legend("topright", colnames(res_clus_df_fastCCLasso_sort)[c(1,5)], fill = sel_colors[c(1,5)], border="white", bty = "n")


#add an "other" column
res_clus_df_fastCCLasso_sort_2 <- as.data.frame(apply(res_clus_df_fastCCLasso_sort, 2, as.numeric))
res_clus_df_fastCCLasso_sort_2$other <- 100 - rowSums(res_clus_df_fastCCLasso_sort_2)
res_clus_df_fastCCLasso_sort_2 <- as.matrix(res_clus_df_fastCCLasso_sort_2)
res_clus_df_fastCCLasso_sort_2 <- t(res_clus_df_fastCCLasso_sort_2)


barplot(res_clus_df_fastCCLasso_sort_2, col=c(sel_colors[1:5], "darkgrey"), border=NA, space=0, axes=F)
legend("bottomleft", 
       legend = rownames(res_clus_df_fastCCLasso_sort_2), bty="n",
       fill = c(sel_colors[1:5], "darkgrey"), cex=0.5, bg = "white")

dev.off()



pdf("./Figure/RM_figures/cluster_ra_barplot.pdf", width = 9, height = 4)
par(mar=c(0.5, 0.5, 0.5, 0.5))
barplot(res_clus_df_fastCCLasso_sort_2, col=c(sel_colors[1:5], "darkgrey"), border=NA, space=0, axes=F)
legend("bottomleft", 
       legend = rownames(res_clus_df_fastCCLasso_sort_2), bty="n",
       fill = c(sel_colors[1:6], "darkgrey"), cex=0.5, bg = "white")
dev.off()



#take clr before doing corrplot
res_clus_df_fastCCLasso_sort_2_clr <- clr(res_clus_df_fastCCLasso_sort_2)

# Visualize the correlation matrix
dev.off()
cor_matrix_clr <- cor(t(res_clus_df_fastCCLasso_sort_2_clr))
corrplot(cor_matrix_clr, method = "circle")



cor_gather(cor_matrix_clr)
# var1                         var2                            cor
#2 Evtepia_Alistipes_Gemmiger   Blautia_Strep              -0.290  
#3 Faecali_Dorea                Blautia_Strep              -0.0290 
#4 Anaerostipes_Roseburia_Bif   Blautia_Strep               0.00756
#5 Bacteroidales_Eisenbergiella Blautia_Strep              -0.312  
#6 other                        Blautia_Strep              -0.485  
#7 Blautia_Strep                Evtepia_Alistipes_Gemmiger -0.290  
#9 Faecali_Dorea                Evtepia_Alistipes_Gemmiger  0.334  
#10 Anaerostipes_Roseburia_Bif  Evtepia_Alistipes_Gemmiger -0.0110 

#one positive correlation between Faecali_Dorea and Evtepia_Alistipes_Gemmiger  0.334
#beyond that negative other and Blautia_Strep is the strongest negative correlation
#Blautia_Strep is generally negative while it is displays the strongest within cluster correlations

cor.test(as.numeric(res_clus_df_fastCCLasso_sort_2_clr["Faecali_Dorea",]),
         as.numeric(res_clus_df_fastCCLasso_sort_2_clr["Evtepia_Alistipes_Gemmiger",]))
#p < 0.0001


#-------------------------------------------------------------------------------
#plots of interest

custom_col <- colorRampPalette(c("#0072B2",'white', "#CC79A7"))

pdf("./Figure//RM_figures/correlation_plots MCCM.pdf", width = 12.5, height=12)

nw_plot
set.seed(42)
heatmap_plot_fastCCLasso
hm
corrplot(cor_matrix_clr, method = "circle", tl.cex = 2, tl.col = "black", col = custom_col(200))

dev.off()



png("./Figure/RM_figures/corrplot_fastCCLasso.png", width = 2*300, height=2*300)
corrplot(cor_matrix_clr, method = "circle", tl.cex = 1.5, tl.col = "black", col = custom_col(200))
dev.off()



#-------------------------------------------------------------------------------
#functional pathways, load them and correlate every pathway with the 6 clusters


load(file = 'Data/data.obj.pathway.RData') 
pathabundance_df <- as.data.frame(data.obj$otu.tab)

dim(pathabundance_df)
#519 1651

pathabundance_names <- rownames(pathabundance_df)

pathabundance_df_sel <- pathabundance_df[,clin_meta_sub$BIOME_with_sequencing_data]
dim(pathabundance_df_sel)
#519 1644


#remove rows with >90% 0's. Also subset pathway names
subset_rows <- which(rowSums(pathabundance_df_sel == 0) < 0.9*ncol(pathabundance_df_sel))
pathabundance_df_sel_sub <- pathabundance_df_sel[subset_rows,]
dim(pathabundance_df_sel_sub)
#397 1644

pathabundance_names <- rownames(pathabundance_df_sel_sub)


#-------------------------------------------------------------------------------
#for every pathway correlate with the 5 clusters, iterate over the clusters and pathways

res_clus_df_fastCCLasso_path <- res_clus_df_fastCCLasso_for_lm[colnames(pathabundance_df_sel),]
dim(res_clus_df_fastCCLasso_path); dim(pathabundance_df_sel_sub)
#1626    6
#397 1626


cor_list <- list()
cor_list[[1]] <- apply(pathabundance_df_sel, 1, function(x) cor.test(x, as.numeric(res_clus_df_fastCCLasso_path[[1]])))
cor_list[[2]] <- apply(pathabundance_df_sel, 1, function(x) cor.test(x, as.numeric(res_clus_df_fastCCLasso_path[[2]])))
cor_list[[3]] <- apply(pathabundance_df_sel, 1, function(x) cor.test(x, as.numeric(res_clus_df_fastCCLasso_path[[3]])))
cor_list[[4]] <- apply(pathabundance_df_sel, 1, function(x) cor.test(x, as.numeric(res_clus_df_fastCCLasso_path[[4]])))
cor_list[[5]] <- apply(pathabundance_df_sel, 1, function(x) cor.test(x, as.numeric(res_clus_df_fastCCLasso_path[[5]])))
cor_list[[6]] <- apply(pathabundance_df_sel, 1, function(x) cor.test(x, as.numeric(res_clus_df_fastCCLasso_path[[6]])))


names(cor_list) <- names(res_clus_df_fastCCLasso_path)

cor_res_df_list <- list()
cor_res_df_list[[1]] <- as.data.frame(cbind(unlist(lapply(cor_list[[1]], function(x) as.numeric(x$estimate))), unlist(lapply(cor_list[[1]], function(x) as.numeric(x$p.value)))))
cor_res_df_list[[2]] <- as.data.frame(cbind(unlist(lapply(cor_list[[2]], function(x) as.numeric(x$estimate))), unlist(lapply(cor_list[[2]], function(x) as.numeric(x$p.value)))))
cor_res_df_list[[3]] <- as.data.frame(cbind(unlist(lapply(cor_list[[3]], function(x) as.numeric(x$estimate))), unlist(lapply(cor_list[[3]], function(x) as.numeric(x$p.value)))))
cor_res_df_list[[4]] <- as.data.frame(cbind(unlist(lapply(cor_list[[4]], function(x) as.numeric(x$estimate))), unlist(lapply(cor_list[[4]], function(x) as.numeric(x$p.value)))))
cor_res_df_list[[5]] <- as.data.frame(cbind(unlist(lapply(cor_list[[5]], function(x) as.numeric(x$estimate))), unlist(lapply(cor_list[[5]], function(x) as.numeric(x$p.value)))))
cor_res_df_list[[6]] <- as.data.frame(cbind(unlist(lapply(cor_list[[6]], function(x) as.numeric(x$estimate))), unlist(lapply(cor_list[[6]], function(x) as.numeric(x$p.value)))))


names(cor_res_df_list) <- names(cor_list) 
cor_res_df_list <- lapply(cor_res_df_list, function(x) x[3:nrow(x),])

colnames(cor_res_df_list[[1]]) <- c("cor", "p_value")
colnames(cor_res_df_list[[2]]) <- c("cor", "p_value")
colnames(cor_res_df_list[[3]]) <- c("cor", "p_value")
colnames(cor_res_df_list[[4]]) <- c("cor", "p_value")
colnames(cor_res_df_list[[5]]) <- c("cor", "p_value")
colnames(cor_res_df_list[[6]]) <- c("cor", "p_value")


cor_res_df_list <- lapply(cor_res_df_list, function(x) cbind(x, "fdr" = p.adjust(x$p_value)))


#--------------------------------------
#make heatmap out of these correlated pathways
estimate_list <- lapply(cor_res_df_list, function(x) x$cor)
fdr_list <- lapply(cor_res_df_list, function(x) x$fdr)

estimate_df <- do.call(cbind, estimate_list)
fdr_df <- do.call(cbind,fdr_list)
rownames(estimate_df) <- rownames(cor_res_df_list[[1]])
rownames(fdr_df) <- rownames(cor_res_df_list[[1]])
name_vec <- names(cor_res_df_list)

#filter by effect size and p-value

sign_cutoffs <- c(0.05, 0.01, 0.001)
h_func <- zicoseq_heatmap(FDR_df = fdr_df, R2_df = estimate_df, sign_cutoffs = sign_cutoffs, grp.labels = name_vec, effect_size_cutoff=0.1,
                          column_name_rot=45)
h_func$heatmap_plot

pdf("./Figure/RM_figures/pathway_plot_top10perc.pdf",width = 11, height=12)
h_func$heatmap_plot
dev.off()


h_func <- zicoseq_heatmap(FDR_df = fdr_df, R2_df = estimate_df, sign_cutoffs = sign_cutoffs, grp.labels = name_vec, effect_size_cutoff=0.2,
                          column_name_rot=45)
pdf("./Figure/RM_figures/pathway_plot_top20perc.pdf",width = 11.5, height=15)
h_func$heatmap_plot
dev.off()



#-------------------------------------------------------------------------------
#load the significant pathways and for every cancer type (compared to other cancers)
#make a plot of those pathways and the subset of the matrix above.


sign_func <- as.data.frame(read.csv("./Figure/subCancerX-Ex_func/pathway_DAA.csv"))


#split into list by variable, all pathways are in here, 362 filtered pathways per cancer type
sign_func_list <- lapply(as.list(unique(sign_func$variable)), function(x) sign_func[which(sign_func$variable == x 
                                                                                          & sign_func$Q < 0.1),])
names(sign_func_list) <- unique(sign_func$variable)

#remove empty ones
sign_func_list <- sign_func_list[lapply(sign_func_list, nrow) > 1]
#for all in sign_func_list plot the heatmap
#dim(estimate_df); dim(fdr_df)

sign_cutoffs <- c(0.1, 0.01, 0.001)

names(sign_func_list)
#"breast"  "liver intrahepatic bile ducts" "multiple myeloma plasma cell"  "pancreas"   

#need to add a bar with the cancer association as well; higher or lower?
all_names <- as.character(unlist(lapply(sign_func_list, function(x) x$pathway)))



i=1 #breast
row_bar = sign_func_list[[i]]$R2
estimate_df_temp <- estimate_df[sign_func_list[[i]]$pathway,]
fdr_df_temp <- fdr_df[sign_func_list[[i]]$pathway,]
#make the heatmap
h_func <- zicoseq_heatmap_bar(FDR_df = fdr_df_temp, R2_df = estimate_df_temp, sign_cutoffs = sign_cutoffs, 
                          grp.labels = name_vec, effect_size_cutoff=1, all_names = all_names, row_bar=row_bar)


png(paste0("./Figure/RM_figures/pathway_plot_", names(sign_func_list)[i],".png"), width = 9.5*300, height=5.75*300, res=300)
h_func$heatmap_plot
dev.off()


#----------------------

i=2 #liver_intrahepatic_bile_ducts
row_bar = sign_func_list[[i]]$R2
estimate_df_temp <- estimate_df[sign_func_list[[i]]$pathway,]
fdr_df_temp <- fdr_df[sign_func_list[[i]]$pathway,]
#make the heatmap
h_func <- zicoseq_heatmap_bar(FDR_df = fdr_df_temp, R2_df = estimate_df_temp, sign_cutoffs = sign_cutoffs, 
                          grp.labels = name_vec, effect_size_cutoff=1, all_names = all_names, row_bar=row_bar)
height_temp <- 1+(0.03*nrow(fdr_df_temp))
width_temp <- 3.5

png(paste0("./Figure/RM_figures/pathway_plot_", names(sign_func_list)[i],".png"), width = width_temp*300, height=height_temp*300)
h_func$heatmap_plot
dev.off()


#----------------------

i=3 #multiple myeloma plasma cell
#first one missing, subset above the previous step
sign_func_list[[i]] <- sign_func_list[[i]][sign_func_list[[i]]$pathway %in% rownames(estimate_df),]
row_bar = sign_func_list[[i]]$R2
estimate_df_temp <- estimate_df[sign_func_list[[i]]$pathway,]
fdr_df_temp <- fdr_df[sign_func_list[[i]]$pathway,]
#make the heatmap
h_func <- zicoseq_heatmap_bar(FDR_df = fdr_df_temp, R2_df = estimate_df_temp, sign_cutoffs = sign_cutoffs, 
                          grp.labels = name_vec, effect_size_cutoff=1, all_names = all_names, row_bar=row_bar)
height_temp <- 3+(0.03*nrow(fdr_df_temp))

png(paste0("./Figure/RM_figures/pathway_plot_", names(sign_func_list)[i],".png"), width = 2.75*300, height=height_temp*300)
h_func$heatmap_plot
dev.off()



#----------------------

i=4 #pancreas
row_bar = sign_func_list[[i]]$R2
estimate_df_temp <- estimate_df[sign_func_list[[i]]$pathway,]
fdr_df_temp <- fdr_df[sign_func_list[[i]]$pathway,]
#make the heatmap
h_func <- zicoseq_heatmap_bar(FDR_df = fdr_df_temp, R2_df = estimate_df_temp, sign_cutoffs = sign_cutoffs, 
                          grp.labels = name_vec, effect_size_cutoff=1, all_names = all_names, row_bar=row_bar)
height_temp <- 1+(0.02*nrow(fdr_df_temp))

png(paste0("./Figure/RM_figures/pathway_plot_", names(sign_func_list)[i],".png"), width = 2.75*300, height=height_temp*300)
h_func$heatmap_plot
dev.off()


