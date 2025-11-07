require(tidyverse)
require(readxl)
require(openxlsx)
require(corrr)
require(reshape2)
require(ggbreak)



file_dir <- dirname(rstudioapi::getSourceEditorContext()$path)
file_dir <- dirname(dirname(dirname(dirname(file_dir))))
setwd(file_dir) #"mforge_clean"


filename <- "Data/colibactin_genefamilies_counts.xlsx"
colibactin_df_counts_MCCM <- read.xlsx(filename, sheet=1, rowNames = T)
dim(colibactin_df_counts_MCCM)

colibactin_df_counts_MCCM[,1:4]

rownames(colibactin_df_counts_MCCM)
#"UniRef90_A0A2T3SSR3"    "UniRef90_Q0P7J7"        "UniRef90_Q0P7J8"        "UniRef90_Q0P7K4"        "UniRef90_Q0P7K7"        "UniRef90_Q0P7K8"       
#"UniRef90_UPI00050AA6A1" "UniRef90_UPI00069CAA6D"


load(file = "Data/data.obj.raw.core.RData") 
clin_meta <- as.data.frame(data.obj$meta.dat)


#-------------------------------------------------------------------------------
# Create correlation data frame of rows colibactin

row_correlations <- correlate(t(colibactin_df_counts_MCCM), method = "spearman")
long_corr <- melt(row_correlations, varnames = c("Row1", "Row2"), value.name = "Correlation")

#if one is found all seem present which is expected when looking at an operon 
hist(long_corr$Correlation, xlim=c(0,1), breaks=5)

median(long_corr$Correlation, na.rm = T)


#-------------------------------------------------------------------------------
#Sum the abundances to get one vector and associate with cancer

colibactin_genes <- colSums(colibactin_df_counts_MCCM)
hist(colibactin_genes, breaks=20, ylim=c(0,20))
table(colibactin_genes > 0)


clin_meta$colibactin_counts <- colibactin_genes
clin_meta$colibactin_binary <- colibactin_genes

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


#--------------------------------------
#despite the normalzation by seq depth explore potential correlation with sequencing depth

#no correlation
plot(clin_meta$depth, clin_meta$colibactin_counts)
cor.test(clin_meta$depth, clin_meta$colibactin_counts)


#--------------------------------------

clin_meta$icd10_first_3_name[is.na(clin_meta$icd10_first_3_name)] <- "Healthy"
clin_meta$cancer_healthy <- "Cancer"
clin_meta$cancer_healthy[clin_meta$icd10_first_3_name == "Healthy"] <- "Healthy"


clin_meta$trunc_cancer_class_names <- clin_meta$icd10_first_3_name
clin_meta$trunc_cancer_class_names <- tolower(gsub("Malignant neoplasm of |Malignant |malignant |Other and | neoplasm| neoplasm of| neoplasms|neoplasm of |\\,|and ", "", clin_meta$trunc_cancer_class_names))
clin_meta$trunc_cancer_class_names[clin_meta$trunc_cancer_class_names %in% c("colon", "rectum")] <- "colorectal"

#logistic regression better?
clin_meta$cancer_healthy <- as.factor(clin_meta$cancer_healthy)
clin_meta$cancer_healthy <- relevel(as.factor(clin_meta$cancer_healthy),ref = "Healthy")


model <- glm(colibactin_binary ~ cancer_healthy + Age, data = clin_meta, family = binomial)
summary(model) #significant 0.00317


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
  geom_segment(aes(y = 800, yend = 800, x = 1, xend = 2), size = 0.3) +
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


cutoff_cancer_class_names <- names(table(clin_meta$trunc_cancer_class_names))[which(table(clin_meta$trunc_cancer_class_names) > 15)]
clin_meta_sub <- clin_meta[clin_meta$trunc_cancer_class_names %in% cutoff_cancer_class_names,]
clin_meta_sub$trunc_cancer_class_names <- relevel(as.factor(clin_meta_sub$trunc_cancer_class_names),ref = "healthy")


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








