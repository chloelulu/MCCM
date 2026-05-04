
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
file_dir <- dirname(dirname(dirname(dirname(file_dir))))
setwd(file_dir) #"mforge_clean"


source("./Code/Submission/MayoOncobiomeStudy/Code/fastCCLasso.R")
source("./Code/Submission/MayoOncobiomeStudy/Code/zicoseq_heatmap_from_df_general.R")
source("./Code/Submission/MayoOncobiomeStudy/Code/zicoseq_heatmap_from_df_general_w_bar.R")


#-------------------------------------------------------------------------------

wd <- "~/MCCM_revisions_reprocessed_data/bracken/"

genus_files <- list.files(wd)[grep("genus", list.files(wd))]

study_names <- gsub("_bracken_genus.tsv", "", genus_files)

genus_tax_list <- lapply(as.list(paste0(wd, genus_files)), function(x) read.table(x, header = T))
names(genus_tax_list) <- study_names

genus_tax_list_renamed <- lapply(
  genus_tax_list,
  function(x) {
    rownames(x) <- x$name      
    x$name <- NULL            
    x                         
  }
)

lapply(genus_tax_list_renamed, dim)

#$jeong
#[1] 950  32

#$kartal
#[1] 1578   76

#$yachida
#[1] 1489  435


#do co-abundance network on yachida
genus_tax_list_renamed$yachida[1:10, 1:10]

tax_table_g <- genus_tax_list_renamed$yachida
tax_table_g_RA <- sweep(tax_table_g, MARGIN = 2, colSums(tax_table_g), '/') * 100


dim(tax_table_g)
tax_table_g_RA[1:10, 1:10]


#-------------------------------------------------------------------------------


#load Yachida metadata; not needed to create the co-abundance network; can be removed


filename <- "~/Dropbox/Mayo_RS/writing/202404 Oncobiome full cohort/manuscript/Cell revision/reprocessing data/Yashida metagenomics 41591_2019_458_MOESM3_ESM.xlsx"
getSheetNames(filename)
#[1] "ENA files"                       "metagenomes Table_S2-1"          "metagenomes Table_S2-1_to_proce"
#[4] "pivot table"                     "Table_S1-1"       

meta_yachida <- as.data.frame(read_xlsx(filename, sheet="metagenomes Table_S2-1_to_proce"))
meta_yachida_sheet1 <- as.data.frame(read_xlsx(filename, sheet="ENA files"))

#head(meta_yachida)
#head(meta_yachida_sheet1)

meta_yachida$run_accession <- as.character(sapply(meta_yachida$Subject_ID, function(x) meta_yachida_sheet1$run_accession[which(meta_yachida_sheet1$Subject_ID == x)]))
rownames(meta_yachida) <- meta_yachida$run_accession

meta_yachida$Healthy_CRC <- if_else(meta_yachida$Group == "Healthy", "Healthy", "CRC")

meta_yachida <- meta_yachida[colnames(tax_table_g_RA),]

dim(meta_yachida); dim(tax_table_g_RA)
#[1] 435  11
#[1] 1489  435



#-------------------------------------------------------------------------------
#Kraken taxonomy
#do at the genus level

cutoff_perc <- ncol(tax_table_g_RA) * 0.4 
prev_genera <- which(rowSums(tax_table_g_RA > 0.05) > cutoff_perc)
length(prev_genera) #57 genera >0.05% in >40% of the samples

prev_genera


#get their prevalence
tax_table_g_RA_prev <- tax_table_g_RA[prev_genera,]
sort(rowSums(tax_table_g_RA_prev > 0.05), decreasing = T)

#get the summed counts
sum(colSums(tax_table_g[prev_genera,])) / sum(colSums(tax_table_g)) * 100
#84.7% off all reads

#clr transform these
tax_table_g_sub <- tax_table_g[names(prev_genera),]


genus_names <- rownames(tax_table_g_sub)
#rownames(tax_table_g_sub) <- genus_names

all_genus_names <- rownames(tax_table_g_RA)
#rownames(tax_table_g) <- all_genus_names


#-------------------------------------------------------------------------------
#Run the correlations using fastCCLasso

dim(tax_table_g_sub)

#this runs 4 methods and takes some time
method_name <- c("fastCCLasso","SparCC","CCLasso","COAT");
method_num <- length(method_name);

set.seed(1212); 

#fileout = "fastCCLasso_output.csv";
#use half min as pseudocount and scale to between 0 and 1, not 0 and 100%
OTUdata <- t(tax_table_g_RA_prev) / 100  # back to [0,1]
min_nonzero <- min(OTUdata[OTUdata > 0])
pc <- min_nonzero / 2
xMat <- OTUdata + pc


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

dim(OTUdata) #57 genera

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
#names(cluster_list) <- c("Evtepia\nAlistipes\nGemmiger", "Faecalibacterium\nDorea", "Blautia\nStreptoccus", 
#                         "Anaerostipes\nRoseburia\nBifido", "Bacteroidales\nEisenbergiella")

names(cluster_list) <- c("misc1", "Faecalibacterium\nDorea", 
                         "misc2", "Blautia\nStreptoccus", "Bacteroidales")


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



png("./Figure/RM_figures/heatmap fastCCLasso yachida.png",width = 11.5*300, height=10*300, res=300)
set.seed(42)
heatmap_plot_fastCCLasso
dev.off()

