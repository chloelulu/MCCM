library(tidyverse)
library(dendextend)
library(pROC)
library(arsenal)
library(ggpubr)
library(ggtree)
library(ComplexHeatmap)
library(circlize)
library(ggExtra)
require(ggrepel)
library(superb)
library(ggpubr)
library(ggtext)
library(ggdist)
library(openxlsx)



file_dir <- dirname(rstudioapi::getSourceEditorContext()$path)
file_dir <- dirname(dirname(dirname(dirname(file_dir))))

setwd(file_dir) #"mforge_clean"
wd <- file_dir
mfd <- paste0(file_dir, "/ManuscriptFigures/")
fd <- paste0(file_dir, "/Figure/")
rd <- paste0(file_dir, "/Result/")


elix.names <- c("Elixhauser_CHF","Elixhauser_Arrhythmia","Elixhauser_Valvular","Elixhauser_PHTN","Elixhauser_PVD","Elixhauser_HTN",
                "Elixhauser_NeuroOther","Elixhauser_Pulmonary","Elixhauser_DM","Elixhauser_DMcx","Elixhauser_Hypothyroid","Elixhauser_Renal","Elixhauser_Liver",      
                "Elixhauser_PUD" ,"Elixhauser_Lymphoma","Elixhauser_Mets","Elixhauser_Tumor","Elixhauser_Rheumatic","Elixhauser_Coagulopathy",
                "Elixhauser_Obesity","Elixhauser_WeightLoss","Elixhauser_FluidsLytes","Elixhauser_Anemia","Elixhauser_Alcohol","Elixhauser_Drugs",
                "Elixhauser_Depression")
variable_names <- c(
  "Elixhauser_CHF", "Elixhauser_Arrhythmia", "Elixhauser_Valvular", "Elixhauser_PHTN", 
  "Elixhauser_PVD", "Elixhauser_HTN", "Elixhauser_NeuroOther", 
  "Elixhauser_Pulmonary", "Elixhauser_DM", "Elixhauser_DMcx", "Elixhauser_Hypothyroid", 
  "Elixhauser_Renal", "Elixhauser_Liver", "Elixhauser_PUD", #"Elixhauser_HIV", 
  "Elixhauser_Lymphoma", "Elixhauser_Mets", "Elixhauser_Tumor", "Elixhauser_Rheumatic", 
  "Elixhauser_Coagulopathy", "Elixhauser_Obesity", "Elixhauser_WeightLoss", 
  "Elixhauser_FluidsLytes","Elixhauser_Anemia", 
  "Elixhauser_Alcohol", "Elixhauser_Drugs", "Elixhauser_Depression", "Elix_score", "Charlson_score"
)

detailed_descriptions <- c(
  "Congestive heart failure", "Cardiac arrhythmias", "Valvular disease", "Pulmonary hypertension",
  "Peripheral vascular disorders", "Hypertension (uncomplicated)", "Other neurologic disorders", 
  "Chronic pulmonary disease", "Diabetes without chronic complications", "Diabetes with chronic complications", "Hypothyroidism", 
  "Renal failure", "Liver disease", "Peptic ulcer disease (excluding bleeding)", #"AIDS/HIV infection", 
  "Lymphoma", "Metastatic cancer", "Solid tumor without metastasis", "Rheumatoid arthritis/collagen vascular diseases",
  "Coagulopathy", "Obesity", "Weight loss", 
  "Fluid and electrolyte disorders", "Deficiency anemia", 
  "Alcohol abuse", "Drug abuse", "Depression","ECI score", "CCI score"
)
names_map <- setNames(detailed_descriptions, variable_names)


## ========= Control X vs subCancerX (PERMANOVA R2) (Cell: Figure 2E)=========  
setwd(wd)
source('Code/Stats.R')
try(load_package())
load(file = 'Data/data.obj.raw.core.RData') 

clin_meta <- as.data.frame(data.obj$meta.dat)


samIDs <- rownames(data.obj.rff$meta.dat[data.obj.rff$meta.dat$icd10_first_3_name_short %in% names(which(sort(table(data.obj.rff$meta.dat$icd10_first_3_name_short))>15)),])
data.obj.rff <- subset_data(data.obj.rff, samIDs = samIDs)
dist.obj.rff <- subset_dist(dist.obj.rff, samIDs = samIDs)
cancer.type <- as.vector(unique(data.obj.rff$meta.dat$icd10_first_3_name_short))
setwd(rd)

# Parameters
dist_name <- "WUniFrac"
variable <- "icd10_first_3_name_short"

# Initialize result matrices
dims <- length(cancer.type)
r2_mat <- array(NA, dim = c(dims, dims, 1), dimnames = list(cancer.type, cancer.type, dist_name))
pv_mat <- r2_mat

# Run pairwise dmanova between cancer types
pairs <- combn(cancer.type, 2, simplify = FALSE)

for (pair in pairs) {
  i <- pair[1]
  j <- pair[2]
  
  ind <- data.obj.rff$meta.dat[[variable]] %in% c(i, j)
  data_sub <- subset_data(data.obj.rff, ind)
  dist_sub <- subset_dist(dist.obj.rff, ind)
  
  obj <- dmanova(as.dist(dist_sub[[dist_name]]) ~ data_sub$meta.dat[[variable]])
  
  r2_mat[i, j, dist_name] <- obj$aov.tab[1, 5]
  pv_mat[i, j, dist_name] <- obj$aov.tab[1, 6]
}

# Fill lower triangle of matrix for symmetry
r2_mat[, , dist_name][lower.tri(r2_mat[, , dist_name])] <- 
  t(r2_mat[, , dist_name])[lower.tri(r2_mat[, , dist_name])]

# Replace NAs with 0
r2_full <- r2_mat[, , dist_name]
r2_full[is.na(r2_full)] <- 0

# PCoA from symmetric distance matrix
pcoa <- cmdscale(as.dist(r2_full), k = 3, eig = TRUE)
pve <- round(pcoa$eig[1:3] / sum(abs(pcoa$eig)) * 100, 1)

yy <- as.data.frame(pcoa$points[, 1:3])
colnames(yy) <- paste0("PC", 1:3)
yy$cancer <- rownames(yy)

# Distance matrix for coloring
dist_vec <- as.data.frame(dist(yy[, 1:3]))[1:23, ]

# Color map
col_fun <- colorRamp2(range(dist_vec), c("#0072B2", "#CC79A7"))
colors <- c("healthy" = "#0072B2", col_fun(dist_vec))
names(colors) <- c("healthy", setdiff(cancer.type, "healthy"))


png("../Figure/RM_figures/PCoA centroids cancer classes.png", width = 9 * 300, height = 8.5 * 300, res = 300)

ggplot(yy, aes(x = PC1, y = PC2, fill = cancer)) +
  geom_point(size = 6, shape = 21) +
  scale_x_continuous(limits = range(yy$PC1)) +
  scale_y_continuous(limits = range(yy$PC2)) +
  scale_fill_manual(values = colors) +
  ggrepel::geom_text_repel(aes(label = cancer), size = 5, box.padding = 0.5, max.overlaps = Inf) +
  labs(x = paste0("PC1 (", pve[1], "%)"), 
       y = paste0("PC2 (", pve[2], "%)"), 
       fill = NULL, color = NULL) +
  theme_bw(base_size = 16) +
  theme(panel.grid = element_blank(),
        legend.position = "none",
        panel.border = element_rect(size = 1),
        axis.ticks = element_line(size = 1),
        axis.ticks.length = unit(0.2, "cm"))

dev.off()




## ========= Control X vs subCancerX pathway (PERMANOVA R2) (Cell: Figure S2E) =========  
setwd(wd)
source('Code/Stats.R')
try(load_package())
load(file = 'Data/data.obj.pathway.RData') 

clin_meta <- as.data.frame(data.obj$meta.dat)


samIDs <- rownames(data.obj$meta.dat[data.obj$meta.dat$icd10_first_3_name_short %in% names(which(sort(table(data.obj$meta.dat$icd10_first_3_name_short))>15)),])
data.obj <- subset_data(data.obj, samIDs = samIDs)
dist.obj$raitchison <- as.matrix(vegdist(t(data.obj$otu.tab), method = "robust.aitchison"))
dist.obj <- subset_dist(dist.obj, samIDs = samIDs)
cancer.type <- as.vector(unique(data.obj$meta.dat$icd10_first_3_name_short))
setwd(rd)

# Parameters
dist_name <- "raitchison"
variable <- "icd10_first_3_name_short"

# Initialize result matrices
dims <- length(cancer.type)
r2_mat <- array(NA, dim = c(dims, dims, 1), dimnames = list(cancer.type, cancer.type, dist_name))
pv_mat <- r2_mat

# Run pairwise dmanova between cancer types
pairs <- combn(cancer.type, 2, simplify = FALSE)

for (pair in pairs) {
  i <- pair[1]
  j <- pair[2]
  
  ind <- data.obj$meta.dat[[variable]] %in% c(i, j)
  data_sub <- subset_data(data.obj, ind)
  dist_sub <- subset_dist(dist.obj, ind)
  
  obj <- dmanova(as.dist(dist_sub[[dist_name]]) ~ data_sub$meta.dat[[variable]])
  
  r2_mat[i, j, dist_name] <- obj$aov.tab[1, 5]
  pv_mat[i, j, dist_name] <- obj$aov.tab[1, 6]
}

# Fill lower triangle of matrix for symmetry
r2_mat[, , dist_name][lower.tri(r2_mat[, , dist_name])] <- 
  t(r2_mat[, , dist_name])[lower.tri(r2_mat[, , dist_name])]

# Replace NAs with 0
r2_full <- r2_mat[, , dist_name]
r2_full[is.na(r2_full)] <- 0

# PCoA from symmetric distance matrix
pcoa <- cmdscale(as.dist(r2_full), k = 3, eig = TRUE)
pve <- round(pcoa$eig[1:3] / sum(abs(pcoa$eig)) * 100, 1)

yy <- as.data.frame(pcoa$points[, 1:3])
colnames(yy) <- paste0("PC", 1:3)
yy$cancer <- rownames(yy)

# Distance matrix for coloring
dist_vec <- as.matrix(dist(yy[, 1:3]))[1:23, ]

# Color map
col_fun <- colorRamp2(range(dist_vec), c("#0072B2", "#CC79A7"))
colors <- c("healthy" = "#0072B2", col_fun(dist_vec))
names(colors) <- c("healthy", setdiff(cancer.type, "healthy"))


png("../Figure/RM_figures/PCoA centroids cancer classes pathway.png", width = 9 * 300, height = 8.5 * 300, res = 300)

ggplot(yy, aes(x = PC1, y = PC2, fill = cancer)) +
  geom_point(size = 6, shape = 21) +
  scale_x_continuous(limits = range(yy$PC1)) +
  scale_y_continuous(limits = range(yy$PC2)) +
  scale_fill_manual(values = colors) +
  ggrepel::geom_text_repel(aes(label = cancer), size = 5, box.padding = 0.5, max.overlaps = Inf) +
  labs(x = paste0("PC1 (", pve[1], "%)"), 
       y = paste0("PC2 (", pve[2], "%)"), 
       fill = NULL, color = NULL) +
  theme_bw(base_size = 16) +
  theme(panel.grid = element_blank(),
        legend.position = "none",
        panel.border = element_rect(size = 1),
        axis.ticks = element_line(size = 1),
        axis.ticks.length = unit(0.2, "cm"))

dev.off()





## ===== subCancerX vs control PCA at sample levels (Cell: Figure 2A)=====
setwd(wd)
load("Data/data.obj.raw.core.RData")
setwd(mfd)

#  Principal Coordinate Analysis 
obj <- cmdscale(dist.obj.rff[["WUniFrac"]], k = 3, eig = TRUE)
pve <- round(obj$eig[1:3] / sum(abs(obj$eig)) * 100, 1)
eig <- obj$points
colnames(eig) <- c("PC1", "PC2", "PC3")

#  Merge Metadata 
yy <- merge(eig, data.obj.rff$meta.dat["icd10_first_3_name_short"], by = 0)
xlab <- paste0("PC1 (", pve[1], "%)")
ylab <- paste0("PC2 (", pve[2], "%)")
zlab <- paste0("PC3 (", pve[3], "%)")

#  clr transform
Y <- data.obj.rff$abund.list$Genus
N <- colSums(Y)
m <- nrow(Y)
N.mat <- matrix(rep(N, m), nrow = m, byrow = TRUE)
N.mat[Y > 0] <- 0
Y <- Y + N.mat / N[max.col(N.mat)]
logY <- log2(Y)
W <- t(t(logY) - colMeans(logY))



sub <- W[c("p__Bacteroidota;g__Bacteroides","p__Firmicutes_A;g__Blautia",
           "p__Bacteroidota;g__Prevotella","p__Firmicutes;g__Streptococcus"),]
rownames(sub) <- gsub('.*g__','',rownames(sub))

#  Group Samples 
yy2 <- yy %>%
  mutate(grp = ifelse(icd10_first_3_name_short == "healthy", "Healthy", "Cancer")) %>%
  filter(icd10_first_3_name_short %in% c(cancer.type, "healthy")) %>%
  column_to_rownames("Row.names") %>%
  merge(t(sub), by = 0) %>%
  column_to_rownames("Row.names")

#orange       skyblue   bluishgreen        yellow          blue    vermillion reddishpurple 
#"#E69F00"     "#56B4E9"     "#009E73"     "#F0E442"     "#0072B2"     "#D55E00"     "#CC79A7" 


#sel_colors_plot <- RColorBrewer::brewer.pal(12,"Paired")
sel_colors_plot <- c("#CC79A7", "#0072B2") #reddishpurple and blue


species_data <- yy2[, c("Bacteroides", "Blautia", "Prevotella", "Streptococcus")]
pc_data <- yy2[, c("PC1", "PC2")]
arrows_data <- data.frame(
  species = c("Bacteroides", "Blautia", "Prevotella", "Streptococcus"),
  x_end = apply(species_data, 2, function(s) cor(s, pc_data$PC1, method='spearman')),  
  y_end = apply(species_data, 2, function(s) cor(s, pc_data$PC2, method='spearman')) 
)
scaling_factor <- 0.4
arrows_data$x_end <- arrows_data$x_end * scaling_factor
arrows_data$y_end <- arrows_data$y_end * scaling_factor
p2 <- ggplot(yy2) +       
  geom_point(size = 1.5, aes(x = PC1, y = PC2, fill = grp, color = grp), alpha = 0.5, shape = 21) +       
  ggforce::geom_mark_ellipse(aes(x = PC1, y = PC2, fill = grp, color = grp), expand = unit(2, "mm"), show.legend = TRUE, alpha = 0.05) +       
  scale_x_continuous(limits = c(min(yy$PC1) * 1.2, max(yy$PC1) * 1.2)) +  
  scale_y_continuous(limits = c(min(yy$PC2) * 1.2, max(yy$PC2) * 1.2)) + 
  scale_fill_manual(name = '', values = sel_colors_plot) +       
  scale_color_manual(name = '', values = sel_colors_plot) +       
  labs(x = xlab, y = ylab, fill = '', color = '') +       
  theme_bw() +       
  theme(text = element_text(size = 15, color = "black"),
        axis.text = element_text(size = 15, color = "black"),
        legend.text = element_text(size = 15, color = "black"),
        axis.title = element_text(size = 15, color = "black"),
        panel.grid.major = element_blank(),
        panel.border = element_rect(size = 1),
        axis.ticks = element_line(size = 1),
        axis.ticks.length = unit(0.2, 'cm'),
        panel.grid.minor = element_blank(),
        legend.position = c(0.15, 0.15),
        legend.background = element_blank()) +       
  guides(fill = guide_legend(override.aes = list(size = 4)), color = 'none', alpha = 'none')
p2 <- p2 + 
  geom_segment(data = arrows_data, aes(x = 0, y = 0, xend = x_end, yend = y_end),
               arrow = arrow(length = unit(0.2, "cm")), color = "black") +
  geom_text_repel(data = arrows_data, aes(x = x_end, y = y_end, label = species), size = 5)
p2 <- ggMarginal(p2, type = "density", margins = "both", size = 5, groupFill = TRUE, groupColour=TRUE)


# Calculate new correlations for PC2 vs PC3
pc_data_pc23 <- yy2[, c("PC2", "PC3")]  # Use PC2 and PC3 instead of PC1 and PC2
arrows_data_pc23 <- data.frame(
  species = c("Bacteroides", "Blautia", "Prevotella", "Streptococcus"),
  x_end = apply(species_data, 2, function(s) cor(s, pc_data_pc23$PC2, method = 'spearman')),  # Correlation with PC2
  y_end = apply(species_data, 2, function(s) cor(s, pc_data_pc23$PC3, method = 'spearman'))   # Correlation with PC3
)

# Apply the scaling factor
scaling_factor <- 0.4
arrows_data_pc23$x_end <- arrows_data_pc23$x_end * scaling_factor
arrows_data_pc23$y_end <- arrows_data_pc23$y_end * scaling_factor

# Modify axis labels for PC2 and PC3
xlab_pc23 <- paste0('PC2 (', pve[2], '%)')
ylab_pc23 <- paste0('PC3 (', pve[3], '%)')

# Create the plot for PC2 vs PC3
p3 <- ggplot(yy2) +       
  geom_point(size = 1.5, aes(x = PC2, y = PC3, fill = grp, color = grp), alpha = 0.5, shape = 21) +       
  ggforce::geom_mark_ellipse(aes(x = PC2, y = PC3, fill = grp, color = grp), expand = unit(2, "mm"), show.legend = TRUE, alpha = 0.05) +       
  scale_x_continuous(limits = c(min(yy$PC2) * 1.2, max(yy$PC2) * 1.2)) +   # Set limits for PC2
  scale_y_continuous(limits = c(min(yy$PC3) * 1.2, max(yy$PC3) * 1.2)) +  # Set limits for PC3
  scale_fill_manual(name = '', values = sel_colors_plot) +       
  scale_color_manual(name = '', values = sel_colors_plot) +       
  labs(x = xlab_pc23, y = ylab_pc23, fill = '', color = '') +       
  theme_bw() +       
  theme(text = element_text(size = 15, color = "black"),
        axis.text = element_text(size = 15, color = "black"),
        legend.text = element_text(size = 15, color = "black"),
        axis.title = element_text(size = 15, color = "black"),
        panel.grid.major = element_blank(),
        panel.border = element_rect(size = 1),
        axis.ticks = element_line(size = 1),
        axis.ticks.length = unit(0.2, 'cm'),
        panel.grid.minor = element_blank(),
        legend.position = c(0.15, 0.15),
        legend.background = element_blank()) +       
  guides(fill = guide_legend(override.aes = list(size = 4)), color = 'none', alpha = 'none')

p3 <- p3 + 
  geom_segment(data = arrows_data_pc23, aes(x = 0, y = 0, xend = x_end, yend = y_end),
               arrow = arrow(length = unit(0.2, "cm")), color = "black") +
  geom_text_repel(data = arrows_data_pc23, aes(x = x_end, y = y_end, label = species), size = 5)
p3 <- ggMarginal(p3, type = "density", margins = "both", size = 5, groupFill = TRUE, groupColour=TRUE)


p <- ggarrange(p2, p3, common.legend = F, nrow = 1)

# Add title with more formatting options
p <- annotate_figure(p, top = text_grob("Weighted Unifrac beta diversity\nCancer vs. healthy", 
                color = "black", size = 14, lineheight = 1))
print(p)

ggsave('../Figure/RM_figures/PCoA_v2_new_names.png', plot = p, width = 11, height = 6)





## ==== Control vs subCancerX Alpha plot (Cell: Figure 2D right)====
setwd(rd)
setwd('subCancerX_Control')
tm = load('Control-lymphoid leukemia/Alpha/icd10_first_3_name_short/Alpha.RData')
dirs <- list.dirs(full.names = F, recursive = F)
alpha <- NULL
for(dir in dirs){
  setwd(rd)
  setwd('subCancerX_Control')
  setwd(dir)
  load('Alpha/icd10_first_3_name_short/Alpha.RData')
  tmp.alpha <- merge(alpha.obj1, data.obj.rff2$meta.dat[,'icd10_first_3_name_short', drop =F], by = 0) %>% 
    dplyr::filter(icd10_first_3_name_short!='healthy')
  alpha <- rbind(alpha, tmp.alpha)
}

load('Alpha/icd10_first_3_name_short/Alpha.RData')
tmp.alpha <- merge(alpha.obj1, data.obj.rff2$meta.dat[,'icd10_first_3_name_short', drop =F], by = 0)%>%
  dplyr::filter(icd10_first_3_name_short=='healthy')
alpha <- rbind(alpha, tmp.alpha)

col <- c(brewer.pal(12,'Paired'),brewer.pal(8,'Dark2'), brewer.pal(8,'Accent')[1:4])
names(col) <- c('healthy', as.vector(unique(alpha$icd10_first_3_name_short))[!(as.vector(unique(alpha$icd10_first_3_name_short)) %in% 'healthy')])

setwd(fd)
load('subCancerX_Control/Alpha_P_R2.RData')
pval <- t(pval_coef_adj.All['Shannon',,]['P',,drop =F])
rownames(pval) <- gsub('Control-','',rownames(pval))
pval <- as.data.frame(pval)  %>% mutate(sig = ifelse(P<0.05, '*',''))

#Shannon
alpha$icd10_first_3_name <- as.factor(alpha$icd10_first_3_name_short)
rownames(pval) <- gsub("_", " ", rownames(pval))

ord <- as.vector((aggregate(Shannon ~ icd10_first_3_name_short, data=alpha, function(x) median(x)) %>% arrange(Shannon))[,1])
alpha <- within(alpha, icd10_first_3_name_short <- factor(icd10_first_3_name_short,levels =ord))
Means <- alpha %>% group_by(icd10_first_3_name_short) %>% summarize(Avg = median(Shannon))

Means$icd10_first_3_name_short
length(unique(alpha$icd10_first_3_name_short)) #23

#color scale
col_fun <- colorRamp2(c(max(Means$Avg), min(Means$Avg)), c("#0072B2", "#CC79A7"))
col <- col_fun(Means$Avg)
names(col) <- as.character(Means$icd10_first_3_name_short)


p1 <- ggplot() + 
  geom_violin(
    data = alpha,
    aes(y = icd10_first_3_name_short, x = Shannon, fill = icd10_first_3_name_short),
    trim = FALSE,
    color = NA  # No border around the violins
  ) +
  geom_line(
    data = Means,
    mapping = aes(y = icd10_first_3_name_short, x = Avg, group = 1),
    color = "black", size = 0.5  # Adjust line width as needed
  ) +
  scale_fill_manual(values = col) +
  theme_bw(base_size = 8) +  # base R style
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, size = 0.5),  # Set border linewidth
    axis.line = element_line(color = "black", size = 0.1),
    text = element_text(size = 9, color = "black", family = "Arial"),
    axis.text = element_text(size = 9, color = "black"),
    axis.title.x = element_text(size = 9, color = "black"),
    axis.title.y = element_text(size = 9, color = "black"),
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA),
    legend.position = "none"
  ) +
  labs(x = "alpha diversity (Shannon)", y = "") +
  geom_text(
    data = pval %>%
      rownames_to_column('icd10_first_3_name_short') %>%
      mutate(pos = max(alpha$Shannon) * 1.12),
    aes(y = icd10_first_3_name_short, x = pos, label = sig),
    size = 4, color = "black"
  )


ggsave(file = paste0('./RM_figures/alpha cancer class vs control violin main text.png'), plot=p1, width=4, height=4, bg = "transparent")




## === cancer vs healthy boxplot (Cell: Figure 2D left)=====

alpha$cancer_healthy <- "cancer"
alpha$cancer_healthy[alpha$icd10_first_3_name_short == "healthy"] <- "healthy"

table(alpha$cancer_healthy) #only for the cancer classes


model <- lm(Shannon ~ cancer_healthy, data = alpha)
summary(model) #significant 

sel_colors_plot <- c("#0072B2", "#CC79A7")



set.seed(42)

alpha_plot <- ggplot(alpha, aes(x = cancer_healthy,
                                y = Shannon,
                                fill = cancer_healthy)) +
  geom_violin(trim = FALSE, color = NA) +                          # violin, no border
  #geom_jitter(aes(color = cancer_healthy),
  #            width = 0.15, size = 2, alpha = 0.8) +              # jittered points
  stat_summary(fun = median, geom = "point", shape = 95,           # thick horizontal tick
               size = 8, color = "black") +                       # black; increase size as needed
  scale_fill_manual(values = sel_colors_plot) +
  scale_color_manual(values = sel_colors_plot) +
  geom_segment(aes(y = 6, yend = 6, x = 1, xend = 2), size = 0.3) +
  geom_text(aes(y = 6*1.03, x = 1.5, label = "p = <0.00001"), size = 4, vjust = 0) +
  theme_classic() +
  theme(
    panel.border = element_rect(color = "black", fill = NA, size = 0.5), 
    axis.ticks.x = element_line(color = "black"),
    axis.line.x  = element_line(color = "black", size = 0.1),
    axis.line.y  = element_line(color = "black", size = 0.1),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "none",
    plot.title  = element_text(size = 12, hjust = 0.5, family = "Arial", colour = "black"),
    axis.text   = element_text(size = 12, family = "Arial", colour = "black"),
    axis.title  = element_text(size = 12, family = "Arial", colour = "black"),
    plot.margin = grid::unit(c(1,1,1,1), "lines")
  ) +
  labs(
    y = "alpha diversity (Shannon)",
    x = NULL
  )

ggsave("./RM_figures/boxplot_cancer_vs_control.png", alpha_plot, width = 2.5, height = 5, units = "in", dpi = 300)




## === Elix components tree (Cell: Figure 3B) =====
elix.names <- c("Elixhauser_CHF","Elixhauser_Arrhythmia","Elixhauser_Valvular","Elixhauser_PHTN","Elixhauser_PVD","Elixhauser_HTN",
                "Elixhauser_NeuroOther","Elixhauser_Pulmonary","Elixhauser_DM","Elixhauser_DMcx","Elixhauser_Hypothyroid","Elixhauser_Renal","Elixhauser_Liver",      
                "Elixhauser_PUD" ,"Elixhauser_Lymphoma","Elixhauser_Mets","Elixhauser_Tumor","Elixhauser_Rheumatic","Elixhauser_Coagulopathy",
                "Elixhauser_Obesity","Elixhauser_WeightLoss","Elixhauser_FluidsLytes","Elixhauser_Anemia","Elixhauser_Alcohol","Elixhauser_Drugs",
                "Elixhauser_Depression")
variable_names <- c(
  "Elixhauser_CHF", "Elixhauser_Arrhythmia", "Elixhauser_Valvular", "Elixhauser_PHTN", 
  "Elixhauser_PVD", "Elixhauser_HTN", "Elixhauser_NeuroOther", 
  "Elixhauser_Pulmonary", "Elixhauser_DM", "Elixhauser_DMcx", "Elixhauser_Hypothyroid", 
  "Elixhauser_Renal", "Elixhauser_Liver", "Elixhauser_PUD", 
  "Elixhauser_Lymphoma", "Elixhauser_Mets", "Elixhauser_Tumor", "Elixhauser_Rheumatic", 
  "Elixhauser_Coagulopathy", "Elixhauser_Obesity", "Elixhauser_WeightLoss", 
  "Elixhauser_FluidsLytes","Elixhauser_Anemia", 
  "Elixhauser_Alcohol", "Elixhauser_Drugs", "Elixhauser_Depression", "Elix_score", "Charlson_score"
)
detailed_descriptions <- c(
  "Congestive heart failure", "Cardiac arrhythmias", "Valvular disease", "Pulmonary hypertension",
  "Peripheral vascular disorders", "Hypertension (uncomplicated)", "Other neurologic disorders", 
  "Chronic pulmonary disease", "Diabetes without chronic complications", "Diabetes with chronic complications", "Hypothyroidism", 
  "Renal failure", "Liver disease", "Peptic ulcer disease (excluding bleeding)",
  "Lymphoma", "Metastatic cancer", "Solid tumor without metastasis", "Rheumatoid arthritis/collagen vascular diseases",
  "Coagulopathy", "Obesity", "Weight loss", 
  "Fluid and electrolyte disorders", "Deficiency anemia", 
  "Alcohol abuse", "Drug abuse", "Depression","ECI score", "CCI score"
)

names_map <- setNames(detailed_descriptions, variable_names)

#adjusted the cutoffs to make it more legible
setwd(wd)
load('Data/data.obj.raw.core.RData')
unique.f = F
cutoff = 0.1
setwd(fd)
setwd('CancerOnly/')
tree.tmp0 <- read.csv('Species_DAA.csv') %>% dplyr::inner_join(as.data.frame(data.obj$otu.name)%>% rownames_to_column('species')) %>% dplyr::filter(variable %in% elix.names)

xx <- read.csv('Species_DAA.csv') %>% dplyr::inner_join(as.data.frame(data.obj$otu.name)%>% rownames_to_column('species')) %>% dplyr::filter(variable =='Elix_score')
tree.tmp0 <- rbind(tree.tmp0, xx)
tree.tmp0$variable <- names_map[tree.tmp0$variable]

tree.tmp <- tree.tmp0 %>% dplyr::select(species, Q, variable) %>% spread(key = 'variable', value = 'Q') %>% column_to_rownames('species')
tree.tmp <- tree.tmp[apply(tree.tmp, 1, function(x) sum(x<cutoff, na.rm = T) > 0),]
tree.tmp[is.na(tree.tmp)] <- 1
idx <- names(which(apply(tree.tmp, 2, function(x) sum(x<cutoff, na.rm = T)) > 0))
tree.tmp <- tree.tmp[,idx]
tree.tmp <- tree.tmp[,c('ECI score', colnames(tree.tmp)[!(colnames(tree.tmp) %in% 'ECI score')])]

R2 <- tree.tmp0 %>% dplyr::select(species, R2, variable) %>% spread(key = 'variable', value = 'R2') %>% column_to_rownames('species')
R2 <- R2[rownames(tree.tmp),colnames(tree.tmp)]
colSums(tree.tmp<cutoff)

## calculate positive and negative for each ECI component for manuscript content use
for(i in names(sort(colSums(tree.tmp<cutoff), decreasing = T))){
  x1 <- tree.tmp[,i,drop =F]
  x2 <- R2[rownames(x1[which(x1[,1]<cutoff),, drop =F]),i,drop =F]
  cat(i)
  print(table(x2>0))
  cat('-------------------\n')
}

## Only keep 25% taxa with Q<cutoff
keep.tax <- c()
for(i in colnames(tree.tmp)){
  topQ <- which(tree.tmp[,i]<cutoff) 
  keepR2 <- rownames(R2[topQ,i, drop =F] %>% arrange(-abs(!!as.name(i))))[1:round(length(topQ)*0.25)]
  keep.tax <- c(keep.tax, keepR2)
}
keep.tax <- unique(keep.tax)

R2 <- R2[keep.tax,]
tree.tmp <- tree.tmp[keep.tax,,]

R2[] <- lapply(R2, function(x) ifelse(x > 0, "Positive", "Negative"))

if(unique.f ==T){
  # subset the unique finding in each cancer
  idx <- names(which(apply(tree.tmp, 1, function(x) sum(x <cutoff, na.rm = T))==1))
  tree.tmp <- tree.tmp[idx,]
  tree.tmp0 <- tree.tmp0[tree.tmp0$species %in% idx,]
  
}
tree.tmp[] <- lapply(tree.tmp, function(x) ifelse(x > cutoff, "No", "Yes"))

tips <- data.obj$tree$tip.label[!(data.obj$tree$tip.label %in% rownames(tree.tmp))]
tree.sub <- drop.tip(data.obj$tree, tips)
tree.sub$tip.label <- gsub('s\\_\\_','',tree.sub$tip.label)
branch.level <- 'Order'
sig.nodes <- gsub('s\\_\\_','',unique(tree.tmp0$species))
otu.name <- (data.obj$otu.name)
rownames(otu.name) <- gsub('s\\_\\_','',rownames(otu.name))
taxa.split <- otu.name[sig.nodes,branch.level]
taxa.split <- gsub('f\\_\\_|p\\_\\_|o\\_\\_|c\\_\\_|g\\_\\_','',taxa.split)
taxa.split <- lapply(split(taxa.split, taxa.split), function(x) names(x))

tree.sub <- tidytree::groupOTU(tree.sub, taxa.split, group_name = 'taxa.split')
rownames(tree.tmp) <- gsub('s\\_\\_','',rownames(tree.tmp))
tree.tmp <- tree.tmp[tree.sub$tip.label,]

dd <- tree.tmp0 %>% dplyr::select(species, Q, variable) %>% mutate(sig = ifelse(Q<cutoff, 'yes','no')) %>% dplyr::select(-Q)
dd$species <- gsub('s\\_\\_','',dd$species)
dd <- within(dd, species <- factor(species, levels = tree.sub$tip.label))
ord <- gsub('f\\_\\_|p\\_\\_|o\\_\\_|c\\_\\_|g\\_\\_','',unique(otu.name[get_taxa_name(ggtree(tree.sub)),branch.level]))


gc()
plt1 <- ggtree(tree.sub, layout = 'circular') + 
  geom_tiplab(align=T,aes(color =taxa.split), fontface='italic',size =4.5) +
  scale_color_manual(values = c(brewer.pal(12,'Paired')[-c(1,11)],brewer.pal(8,'Dark2')[-4],brewer.pal(8,'Set2')), breaks = ord) + 
  labs(color = branch.level) +
  theme(legend.position="right", legend.text = element_text(size = 16, color = 'black', face = 'italic'),legend.title = element_text(size = 16, color = 'black')) +
  guides(color = guide_legend(override.aes = list(size = 4))) 
plt1 <- rotate_tree(plt1, 360)


rownames(R2) <- gsub('s\\_\\_','',rownames(R2))
R2 <- R2[rownames(tree.tmp),]
test_match_order(rownames(R2), rownames(tree.tmp))

for (col in colnames(tree.tmp)) {
  tree.tmp[[col]] <- ifelse(tree.tmp[[col]] == "Yes" & R2[[col]] == "Positive", paste0(col," (Positive)"),
                            ifelse(tree.tmp[[col]] == "Yes" & R2[[col]] == "Negative", paste0(col," (Negative)"), 
                                   tree.tmp[[col]]))
}


xx2 <- c()
for(i in 1:ncol(tree.tmp)){
  tt <- unique(tree.tmp[,i])
  tt <- tt[!(tt %in% 'No')]
  xx2 <- c(xx2,tt)
}
xx2 <- unique(xx2)

dup <- names(which(table(gsub('\\(.*','',xx2))==2))
grep(dup[1], xx2)
pos <-  xx2[grep('Positive',xx2)]
neg <- xx2[grep('Negative',xx2)] 


#"ECI score (Positive)"                                "ECI score (Negative)"                                
#"Cardiac arrhythmias (Positive)"                      "Cardiac arrhythmias (Negative)"                      
#"Deficiency anemia (Positive)"                         
#"Diabetes with chronic complications (Positive)"      "Diabetes with chronic complications (Negative)"  

#"Diabetes without chronic complications (Positive)"   
#"Fluid and electrolyte disorders (Positive)"          "Fluid and electrolyte disorders (Negative)"          
#"Liver disease (Negative)"                            "Liver disease (Positive)"                            
#"Obesity (Positive)"                                   
#"Other neurologic disorders (Positive)"                "Other neurologic disorders (Negative)"                
#"Peptic ulcer disease (excluding bleeding) (Positive)" "Peptic ulcer disease (excluding bleeding) (Negative)" "Renal failure (Positive)"                            
#"Valvular disease (Positive)"                          
#"Weight loss (Positive)"                              "Weight loss (Negative)" 


col <- c(brewer.pal(8,'Reds')[c(8,6)], # red: "ECI score (Positive)" ,"ECI score (Negative)"  
         brewer.pal(8,'Purples')[c(8,6)], # purple:Cardiac arrhythmias (Positive),Cardiac arrhythmias (Negative)
         brewer.pal(8,'Set1')[c(8)], #pink "Deficiency anemia (Positive)"
         brewer.pal(8,'Oranges')[c(6,5)], # Oranges: Diabetes with chronic complications (Positive), Negative
         brewer.pal(8,'Dark2')[c(6)],#"Diabetes without chronic complications (Positive)"
         brewer.pal(8,'Greens')[c(8,4)],# green: Fluid and electrolyte disorders (Positive), Negative
         brewer.pal(8,'Greys')[c(6,8)],#grey:  "Liver disease (Negative)" "Liver disease (Positive)" 
         brewer.pal(8,'BuPu')[c(8)], #"Obesity (Positive)" 
         brewer.pal(8,'RdPu')[c(4,6)],# PurpleRed: "Other neurologic disorders (Positive), Negative
         brewer.pal(8,'Dark2')[c(7,6)],# Yellow: "Peptic ulcer disease (excluding bleeding) (Positive), Negative
         brewer.pal(8,'Greens')[c(3,4)],#"Valvular disease (Positive)"  
         brewer.pal(8,'Blues')[c(8,4)], #Blues: Weight loss (Positive), Negative
         brewer.pal(8,'BuGn')[c(8,4)]
)[1:length(c(pos,neg))]
col1 <- c(col,'white')
names(col1) <- c(xx2,'No')

breaks_to_show <- xx2


gc()
gheatmap(plt1, tree.tmp, offset=3, width=1, font.size=3, color = brewer.pal(9,'Set1')[9], colnames = F) +
  scale_fill_manual(values=col1, name="",breaks=breaks_to_show, 
                    guide = guide_legend(label.theme = element_text(angle = 0, face = 'plain'))) 

setwd(mfd)
if(cutoff ==0.1){width = height = 30}
if(cutoff ==0.05){width = height = 16}
ggsave(paste0('../Figure/RM_figures/elix_components circle plot_RM.png'), width = width, height = height)
gc()


## ==== Add some Species for ECI components (Cell: Figure S3DF)=======
# The negative taxonomic associations were predominantly found in species belonging to the Oscillospirales, Cordiobacteriales, and Lachnospirales, whereas positive correlations were mainly observed in species belonging to orders Veillonellales, Lactobacillales, and Enterobacteriales.
# We could plot one of the Escherichia species (maybe dysenteriae or E. coli D?) from the circle plot for species significant with ECI components.
## It should include weight loss (positive). Also include Klebsiella pneumoniae which seems like a good example.
setwd(fd)
load('CancerOnly/DAA_P_R2.RData')
setwd(rd)
load('CancerOnly/data.obj.wk.RData')
colnames(Q.All$Species)


## All significant species
ECIs <- colnames(Q.All$Species)[grep('^Elixhauser',colnames(Q.All$Species))]
ECIs <- names(which(sort(colSums(data.obj$meta.dat[,ECIs]=='Yes')) > 10))
for(ECI in ECIs){
  idx <- Q.All$Species$Species[which(Q.All$Species[,ECI] <0.1049)]
  if(length(idx)>0){
    pdf(paste0(fd,"RM_figures/All_Significant_species_boxplot_", ECI, ".pdf"), width=2.5, height=3.5)
    abund <- data.obj.rff$otu.tab %>% merge(data.obj.rff$otu.name[,c('Species'), drop =F], by = 0) %>%
      dplyr::select(-c('Row.names')) %>% column_to_rownames('Species')
    prop <- t(t(abund)/colSums(abund))
    for(i in idx){
      plot_df <- sqrt(t(prop[i,,drop =F])) %>% merge(data.obj.rff$meta.dat[,ECI,drop =F], by = 0)
      plot_df[,3] <- as.factor(plot_df[,3])
      
      colnames(plot_df)[2] <- 'taxon'
      p_plot <- ggboxplot(plot_df, ECI, "taxon",
                          color = ECI, palette = brewer.pal(5, "Dark2"), width = 0.8,
                          add = "jitter", add.params = list(alpha = 0.5, size = 1.5),
                          ylab = "sqrt(relative abundance)") +
        scale_y_continuous(trans = 'sqrt') +
        ggtitle(gsub("s__|_", "", i)) +
        theme(
          legend.position = "none",
          axis.text = element_text(family = "sans", size = 10),
          axis.title = element_text(family = "sans", size = 10),
          plot.title = element_text(face = "italic", size = 10)
        ) +
        xlab(colnames(plot_df)[3])
      print(p_plot)
    }
    dev.off()
  }
}



## ==== Add some Pathways for ECI components (Cell: Figure S3EG) =======
setwd(fd)
load('CancerOnly_func/pathway/DAA_P_R2.RData')
setwd(rd)
load('CancerOnly_func/pathway/data.obj.wk.RData')
colnames(Q.All$pathway)


#glycogen biosynthesis

#weight loss is negatively associated with glycogen biosynthesis and positively associated with ketogenesis
## Elixhauser_WeightLoss; glycogen and ketogenesis
idx <- Q.All$pathway$pathway[which(Q.All$pathway$Elixhauser_WeightLoss <0.1049)]
i <- idx[grep('glycogen biosynthesis',idx)]
prop <- data.obj$otu.tab
ECI <- 'Elixhauser_WeightLoss'
plot_df <- sqrt(t(prop[i,,drop =F])) %>% merge(data.obj$meta.dat[,ECI,drop =F], by = 0)
colnames(plot_df)[2] <- 'taxon'
plot_df$Elixhauser_WeightLoss <- as.factor(plot_df$Elixhauser_WeightLoss)

pathway_name <- i
colon_split <- strsplit(pathway_name, ":", fixed = TRUE)
# Process each element
pathway_name_new <- sapply(colon_split, function(parts) {
  # Handle cases: with or without colon
  if (length(parts) < 2) {
    # No colon: just wrap whole string
    gsub("(.{30})", "\\1\n", parts[1])
  } else {
    # With colon: add newline after colon, then wrap the second part
    paste0(parts[1], ":\n", gsub("(.{30})", "\\1\n", trimws(parts[2])))
  }
})

#set ylims to 8 in order to remove one outlying patients
p_plot <- ggboxplot(plot_df, ECI, "taxon",
                    color = ECI, palette = brewer.pal(5, "Dark2"), width = 0.8,
                    add = "jitter", add.params = list(alpha = 0.5, size = 1.5),
                    ylab = "sqrt(reads per million, CPM)") +
  scale_y_continuous(trans = 'sqrt', limits = c(8, NA)) +
  ggtitle(gsub(".*: ", "", pathway_name_new)) +
  theme(
    legend.position = "none",
    axis.text = element_text(family = "sans", size = 10),
    axis.title = element_text(family = "sans", size = 10),
    plot.title = element_text(face = "italic", size = 10)
  ) 
print(p_plot)
ggsave(paste0(fd,"RM_figures/Significant_pathway_boxplot_", ECI, '_',i,".pdf"), width=2.5, height=3.5)


#p_plot <- ggplot(plot_df, aes(x = Elixhauser_WeightLoss, y = taxon, fill = Elixhauser_WeightLoss)) +
# geom_violin(trim = FALSE, alpha = 0.4) +
# geom_boxplot(width = 0.1, outlier.shape = NA, alpha = 0.6) +
# stat_summary(fun = median, geom = "point", shape = 18, size = 3, color = "black") +
# ggtitle(gsub("s__", "", i)) +
# theme_classic() +
# scale_y_continuous(trans = 'sqrt') +
# theme(
#   legend.position = "none",
#   axis.text = element_text(family = "sans", size = 10),
#   axis.title = element_text(family = "sans", size = 10),
#   plot.title = element_text(face = "italic", size = 10)
# ) + labs(y = 'sqrt(reads per million, CPM))')
#print(p_plot)
#ggsave(paste0(fd,"RM_figures/Significant_pathway_violinplot_", ECI, '_',i,".pdf"), width=5, height=3.5)


#ketogenesis
idx <- Q.All$pathway$pathway[which(Q.All$pathway$Elixhauser_WeightLoss <0.1049)]
i <- idx[grep('keto',idx)]
prop <- data.obj$otu.tab
ECI <- 'Elixhauser_WeightLoss'
plot_df <- sqrt(t(prop[i,,drop =F])) %>% merge(data.obj$meta.dat[,ECI,drop =F], by = 0)
colnames(plot_df)[2] <- 'taxon'
plot_df$Elixhauser_WeightLoss <- as.factor(plot_df$Elixhauser_WeightLoss)

pathway_name <- i
colon_split <- strsplit(pathway_name, ":", fixed = TRUE)
# Process each element
pathway_name_new <- sapply(colon_split, function(parts) {
  # Handle cases: with or without colon
  if (length(parts) < 2) {
    # No colon: just wrap whole string
    gsub("(.{30})", "\\1\n", parts[1])
  } else {
    # With colon: add newline after colon, then wrap the second part
    paste0(parts[1], ":\n", gsub("(.{30})", "\\1\n", trimws(parts[2])))
  }
})

#set ylims to remove one outlying patients
p_plot <- ggboxplot(plot_df, ECI, "taxon",
                    color = ECI, palette = brewer.pal(5, "Dark2"), width = 0.8,
                    add = "jitter", add.params = list(alpha = 0.5, size = 1.5),
                    ylab = "sqrt(reads per million, CPM)") +
  scale_y_continuous(trans = 'sqrt') +
  ggtitle(gsub(".*: ", "", pathway_name_new)) +
  theme(
    legend.position = "none",
    axis.text = element_text(family = "sans", size = 10),
    axis.title = element_text(family = "sans", size = 10),
    plot.title = element_text(face = "italic", size = 10)
  ) 
print(p_plot)
ggsave(paste0(fd,"RM_figures/Significant_pathway_boxplot_", ECI, '_',i,".pdf"), width=2.5, height=3.5)



#----------------------------------------
#Anemia

idx <- Q.All$pathway$pathway[which(Q.All$pathway$Elixhauser_Anemia <0.1049)]
prop <- data.obj$otu.tab
ECI <- 'Elixhauser_Anemia'
for(i in idx[grep('tetrahydrofolate',idx)]){
  plot_df <- sqrt(t(prop[i,,drop =F])) %>% merge(data.obj$meta.dat[,ECI,drop =F], by = 0) 
  colnames(plot_df)[2] <- 'taxon'
  plot_df[,ECI] <- as.factor(plot_df[,ECI])
  p_plot <- ggboxplot(plot_df, ECI, "taxon",
                      color = ECI, palette = brewer.pal(5, "Dark2"), width = 0.8,
                      add = "jitter", add.params = list(alpha = 0.5, size = 1.5),
                      ylab = "sqrt(reads per million, CPM)") +
    scale_y_continuous(trans = 'sqrt') +
    
    ggtitle(gsub(".*: ", "", i)) +
    theme(
      legend.position = "none",
      axis.text = element_text(family = "sans", size = 10),
      axis.title = element_text(family = "sans", size = 10),
      plot.title = element_text(face = "italic", size = 10)
    ) 
  print(p_plot)
  ggsave(paste0(fd,"RM_figures/Significant_pathway_boxplot_", ECI, '_',i,".pdf"), width=5, height=3.5)
}



idx <- Q.All$pathway$pathway[which(Q.All$pathway$Elixhauser_Anemia <0.1049)]
prop <- data.obj$otu.tab
ECI <- 'Elixhauser_Anemia'
for(i in idx[grep('heme',idx)]){
  plot_df <- sqrt(t(prop[i,,drop =F])) %>% merge(data.obj$meta.dat[,ECI,drop =F], by = 0) 
  colnames(plot_df)[2] <- 'taxon'
  plot_df[,ECI] <- as.factor(plot_df[,ECI])
  
  pathway_name <- i
  colon_split <- strsplit(pathway_name, ":", fixed = TRUE)
  # Process each element
  pathway_name_new <- sapply(colon_split, function(parts) {
    # Handle cases: with or without colon
    if (length(parts) < 2) {
      # No colon: just wrap whole string
      gsub("(.{30})", "\\1\n", parts[1])
    } else {
      # With colon: add newline after colon, then wrap the second part
      paste0(parts[1], ":\n", gsub("(.{30})", "\\1\n", trimws(parts[2])))
    }
  })
  
  p_plot <- ggboxplot(plot_df, ECI, "taxon",
                      color = ECI, palette = brewer.pal(5, "Dark2"), width = 0.8,
                      add = "jitter", add.params = list(alpha = 0.5, size = 1.5),
                      ylab = "sqrt(reads per million, CPM)") +
    scale_y_continuous(trans = 'sqrt') +
    
    ggtitle(gsub(".*: ", "", pathway_name_new)) +
    theme(
      legend.position = "none",
      axis.text = element_text(family = "sans", size = 10),
      axis.title = element_text(family = "sans", size = 10),
      plot.title = element_text(face = "italic", size = 10)
    ) 
  print(p_plot)
  ggsave(paste0(fd,"RM_figures/Significant_pathway_boxplot_", ECI, '_',i,".pdf"), width=2.5, height=3.5)
}




### Plot all significant pathways in all ECIs
ECIs <- colnames(Q.All$pathway)[grep('^Elixhauser_', colnames(Q.All$pathway))]
ECIs <- names(which(sort(colSums(data.obj$meta.dat[,ECIs]=='Yes')) > 10))
for(ECI in ECIs){
  idx <- Q.All$pathway$pathway[which(Q.All$pathway[,ECI] <0.1049)]
  prop <- data.obj$otu.tab
  # if(length(table(data.obj$meta.dat[,ECI]))==2){
    pdf(paste0(fd,"RM_figures/All_Significant_pathway_boxplot_", ECI, ".pdf"), width=5, height=3.5)
    for(i in idx){
      plot_df <- sqrt(t(prop[i,,drop =F])) %>% merge(data.obj$meta.dat[,ECI,drop =F], by = 0)
      colnames(plot_df)[2] <- 'taxon'
      plot_df[,ECI] <- as.factor(plot_df[,ECI])
      for(i in idx){
        plot_df <- sqrt(t(prop[i,,drop =F])) %>% merge(data.obj$meta.dat[,ECI,drop =F], by = 0)
        colnames(plot_df)[2] <- 'taxon'
        p_plot <- ggboxplot(plot_df, ECI, "taxon",
                            color = ECI, palette = brewer.pal(5, "Dark2"), width = 0.8,
                            add = "jitter", add.params = list(alpha = 0.5, size = 1.5),
                            ylab = "sqrt(reads per million, CPM)") +
          scale_y_continuous(trans = 'sqrt') +
          ggtitle(gsub(".*: ", "", i)) +
          theme(
            legend.position = "none",
            axis.text = element_text(family = "sans", size = 10),
            axis.title = element_text(family = "sans", size = 10),
            plot.title = element_text(face = "italic", size = 10)
          ) +
          xlab("")
        print(p_plot)
      }
    }
    dev.off()
    
  # }
}




## ==== Figure Alpha Beta main variables bar plot (Cell: Figure S3C)======
vars.incl <- c("Bristol_score","BMI", "Age", "Sex", "Metastasis","PPI_day_365", "Abx_day_365",
            "Abx_last_month","PPI_last_month","Elix_score","Sample_season","Urban","icd10_first_3_name_short","Site")
setwd(wd)
alpha.measure =  "Shannon" 
load('Result/CancerOnly/data.obj.wk.RData')
load('Figure/CancerOnly/Alpha_P_R2.RData')
pval_coef_adj <- t(cbind(pval_coef_adj.All[alpha.measure,,])) 
pval_coef_adj <- pval_coef_adj[!(rownames(pval_coef_adj) %in% 'Batch'),]
pval_coef_adj[pval_coef_adj[,'R2'] < 0,'R2'] <- 0
pval_coef_adj <- as.data.frame(pval_coef_adj) %>% 
  mutate(direction = ifelse(is.na(coef),'none',ifelse(coef>0, 'positive','negative'))) %>%
  dplyr::select(-coef) %>%
  tibble::rownames_to_column('var') %>% melt()
sub.r <- pval_coef_adj[pval_coef_adj$variable=='R2',]
sub.p <- pval_coef_adj[pval_coef_adj$variable=='P',c('var','value')] %>% 
  mutate(P.col = ifelse(value<0.05,'sig','no')) %>% dplyr::select(-value)
sub.rp <- dplyr::inner_join(sub.p, sub.r) %>% mutate(value = value *100)

pval_coef_unadj <- t(cbind(pval_coef_unadj.All[alpha.measure,,])) 
pval_coef_unadj <- pval_coef_unadj[!(rownames(pval_coef_unadj) %in% 'Batch'),]
pval_coef_unadj[pval_coef_unadj[,'R2']<0,'R2'] <- 0
pval_coef_unadj <- as.data.frame(pval_coef_unadj) %>% 
  mutate(direction = ifelse(is.na(coef),'none',ifelse(coef>0, 'positive','negative'))) %>%
  dplyr::select(-coef) %>%
  tibble::rownames_to_column('var') %>% melt()
sub.r <- pval_coef_unadj[pval_coef_unadj$variable=='R2',]
sub.p <- pval_coef_unadj[pval_coef_unadj$variable=='P',c('var','value')] %>% 
  mutate(P.col = ifelse(value<0.05,'sig','no')) %>% dplyr::select(-value)
sub.rp2 <- dplyr::inner_join(sub.p, sub.r) %>% mutate(value = value *100)

sub.rp12 <- rbind(sub.rp %>% mutate(grp = 'adjusted'), sub.rp2 %>% mutate(grp = 'marginal'))
ord <- (sub.rp[order(sub.rp$value),'var'])
sub.rp12 <- within(sub.rp12, var <- factor(var, levels=ord))


# beta - taxa
dist.name <- 'WUniFrac'
load('Figure/CancerOnly/Beta_P_R2.RData')
dist.names <- rownames(R2.adj)
R2_P <- cbind.data.frame(value = c(P.adj[dist.name,]), R2 = c(R2.adj[dist.name,])) %>% 
  rownames_to_column('var') %>% dplyr::filter(var !='Batch') %>% 
  mutate(P.col = ifelse(value<0.05, 'sig','no'))  %>% 
  dplyr::select(-value) %>% mutate(direction = 'none', variable = 'R2') %>% 
  dplyr::rename(value = R2) %>%
  mutate(value=ifelse(value <0,0, value *100))
ord <- (R2_P[order(R2_P$value),'var'])
R2_P2 <- cbind.data.frame(value = c(P.unadj[dist.name,]), R2 = c(R2.unadj[dist.name,])) %>% 
  rownames_to_column('var') %>% dplyr::filter(var !='Batch') %>% 
  mutate(P.col = ifelse(value<0.05, 'sig','no'))  %>% 
  dplyr::select(-value) %>% mutate(direction = 'none', variable = 'R2') %>% 
  dplyr::rename(value = R2) %>%
  mutate(value=ifelse(value <0,0, value *100))
R2_P12 <- rbind(R2_P %>% mutate(grp = 'adjusted'),R2_P2 %>% mutate(grp = 'marginal'))
R2_P12 <- within(R2_P12, var <- factor(var, levels=ord))


# beta - pathway
setwd(wd)
dist.name <- 'BC'
load('Figure/CancerOnly_func/pathway/Beta_P_R2.RData')
R2_P <- cbind.data.frame(value = c(P.adj[dist.name,]), R2 = c(R2.adj[dist.name,])) %>% 
  rownames_to_column('var') %>% dplyr::filter(var !='Batch') %>% 
  mutate(P.col = ifelse(value<0.05, 'sig','no'))  %>% dplyr::select(-value) %>% mutate(direction = 'none', variable = 'R2') %>% dplyr::rename(value = R2) %>%
  mutate(value=ifelse(value <0,0, value *100))
ord <- (R2_P[order(R2_P$value),'var'])
R2_P2 <- cbind.data.frame(value = c(P.unadj[dist.name,]), R2 = c(R2.unadj[dist.name,])) %>% 
  rownames_to_column('var') %>% dplyr::filter(var !='Batch') %>% 
  mutate(P.col = ifelse(value<0.05, 'sig','no'))  %>% dplyr::select(-value) %>% mutate(direction = 'none', variable = 'R2') %>% dplyr::rename(value = R2) %>%
  mutate(value=ifelse(value <0,0, value *100))
R2_P13 <- rbind(R2_P %>% mutate(grp = 'adjusted'),R2_P2 %>% mutate(grp = 'marginal'))
R2_P13 <- within(R2_P13, var <- factor(var, levels=ord))

df <- rbind(R2_P12 %>% mutate(diversity = 'Beta diversity (species)'),
            R2_P13 %>% mutate(diversity = 'Beta diversity (pathway)'), 
            sub.rp12 %>% mutate(diversity = 'Alpha diversity'))
df <- within(df, diversity <- factor(diversity,levels = c('Alpha diversity','Beta diversity (species)','Beta diversity (pathway)')))
df <- within(df, grp <- factor(grp,levels = c('marginal','adjusted')))

head(df)
unique(df$var)
df2 <- df %>% dplyr::filter(var %in% vars.incl) %>% mutate(var = as.character(var))
unique(df2$var)
df2$var[df2$var =="icd10_first_3_name_short"] <- "Cancer class"
df2$var[df2$var =="Bristol_score"] <- "Bristol stool score"
df2$var[df2$var =="PPI_day_365"] <- "PPI in past year (No)"
df2$var[df2$var =="PPI_last_month"] <- "PPI in past month (No)"
df2$var[df2$var =="Abx_day_365"] <- "Antibiotics in past year (No)"
df2$var[df2$var =="Abx_last_month"] <- "Antibiotics in past month (No)"
df2$var[df2$var =="Elix_score"] <- "Elixhauser Comorbidity score"
df2$var[df2$var =="Sample_season"] <- "Sample season"
df2$var[df2$var =="Urban"] <- "Residence type (non-urban)"
df2$var[df2$var =="Sex"] <- "Sex (Female)"
df2$var[df2$var =="Metastasis"] <- "Metastasis (No)"


#plotting settings start here

df3 <- df2
df3 <- df3[!df3$diversity == "Beta diversity (pathway)",]
df3 <- df3[df3$grp == "adjusted",]

#unique(df3$var)
variables_to_plot <- c("Cancer class" , "Bristol stool score", "Elixhauser Comorbidity score", "Sex (Female)", "Age",
                       "Antibiotics in past year (No)", "Antibiotics in past month (No)", "BMI", "PPI in past year (No)", 
                       "Residence type (non-urban)", "Site", "PPI in past month (No)", "Sample season", 
                       "Metastasis (No)") 
sel_colors <- palette.colors(palette = "Okabe-Ito")[2:8]

variable_colors_to_plot <- rep(sel_colors[6], length(variables_to_plot))
names(variable_colors_to_plot) <- variables_to_plot

variable_colors_to_plot[variables_to_plot %in% c("Bristol stool score", "Sample season")] <- sel_colors[2] #technical
variable_colors_to_plot[variables_to_plot %in% c("Sex (Female)", "Age", 
                                                 "Residence type (non-urban)", "BMI", "Site")] <- sel_colors[1] #demographic



main_div_plot <- ggplot(df3, aes(x = reorder(var, value), y = value)) +
  geom_segment(aes(xend = reorder(var, value), y = 0,  yend = value),
               size = 4, colour = "darkgrey") +
  geom_point(aes(shape = direction, fill = P.col), size = 5, stroke = 0.5) +
  scale_fill_manual(values = c(sig = "black", no = "white")) +
  scale_shape_manual(values = c(positive = 24, negative = 25, none = 21)) +
  facet_grid(grp ~ diversity, scales = "free") +
  coord_flip() +                     # <-- flip to horizontal
  theme_classic() +
  labs(
    y = "Variability explained (R²), %",
    x = "",
    fill  = "",
    shape = "",
    subtitle = "",
    caption = paste0(
      "<span style='color:", sel_colors[6], "'>Clinical variables</span><br>",
      "<span style='color:", sel_colors[1], "'>Demographic variables</span><br>",
      "<span style='color:", sel_colors[2], "'>Technical variables</span>"
    )
  ) +
  theme(
    panel.grid   = element_blank(),
    axis.title   = element_text(size = 16, colour = "black"),
    legend.text  = element_text(size = 16, colour = "black"),
    strip.text   = element_text(size = 16, colour = "black"),
    axis.text.x  = element_text(size = 16, colour = "black"),
    axis.text.y  = element_text(size = 16, colour = rev(variable_colors_to_plot)),
    plot.caption = element_markdown(size = 16, hjust = 0),
    panel.spacing = unit(2, "lines")
  )

main_div_plot


ggsave(file = 'Figure/RM_figures/Figure1_main.svg', plot=main_div_plot, width = 12, height = 5.5)
ggsave(file = 'Figure/RM_figures/Figure1_main.png', plot=main_div_plot, width = 12, height = 5.5)




## ==== elix components beta div plot (Cell: Figure 3C)====
x1 <- gsub('Elixhauser_','',rev((df[df$var %in% elix.names & df$grp =='adjusted' & df$diversity =='Beta diversity(taxa)' & df$P.col =='sig',] %>% arrange(value))[,1]))
x2 <- (gsub('Elixhauser_','',rev((df[df$var %in% elix.names & df$grp =='adjusted' & df$diversity =='Beta diversity(pathway)' & df$P.col =='sig',] %>% arrange(value))[,1])))
intersect(x1,x2)
paste0(x2, collapse = ', ')
df.tmp <- df[df$var %in% elix.names & df$grp =='adjusted' & df$diversity != 'Alpha diversity',]
df.tmp$var <- droplevels(df.tmp$var)
unique(df.tmp$var)
df.tmp$var <- names_map[as.vector(df.tmp$var)]

df.tmp <- df.tmp[!df.tmp$var == "Lymphoma",]
df.tmp <- df.tmp[!df.tmp$var == "Solid tumor without metastasis",]

df.tmp_2 <- df.tmp
df.tmp_2 <- df.tmp_2[df.tmp_2$diversity == "Beta diversity (species)",]

elix_yes_counts <- colSums(data.obj$meta.dat[, elix.names] == "Yes")
names(elix_yes_counts) <- names_map[names(elix_yes_counts)]

df.tmp_2$var <- paste0(df.tmp_2$var, "  n=", elix_yes_counts[as.character(df.tmp_2$var)])

elix_div_plot <- df.tmp_2 %>% 
  ggplot(aes(x =reorder(var,value), y= value)) + 
  geom_segment(aes(x=reorder(var,value), xend=reorder(var,-value), y=0, yend=value), size =4, color = "darkgrey") +
  geom_point(size=4, shape=21, aes(fill = P.col), stroke=0.5) + 
  scale_fill_manual(values = c('sig'='black','no'= 'white')) +
  scale_shape_manual(values = c('positive'=24,'negative'=25, 'none'= 21)) +
  #facet_grid(~diversity, scales = 'free')  + 
  theme_classic() +
  labs(y = ('Variability explained (R²), %'), fill = '', x = '', shape ='') +
  theme(plot.background = element_rect(fill = "transparent", colour = NA),
        panel.grid = element_blank(),
        axis.title = element_text(size=14, color = 'black'),
        legend.text = element_text(size=14, color = 'black'),
        strip.text = element_text(size=14, color = 'black'),
        axis.text.x = element_text(size=14,color = 'black'),
        axis.text.y = element_text(size=14, color = 'black')) + 
  coord_flip() + guides(fill = FALSE) + ggtitle('Beta diversity (species)') +
  theme(plot.margin = ggplot2::margin(t = 20, r = 30, b = 10, l = 10)) 

elix_div_plot

ggsave(file = paste0('Figure/RM_figures/main_elix_div_plot.png'), plot=elix_div_plot, width = 9, height = 6, bg = "transparent")



## ==== subCancerX vs all other cancers alpha & beta species summary plot(Cell: Figure S4A)====
setwd(wd)
alpha.measure <- 'Shannon'
dist.name <- 'WUniFrac'

## Alpha diversity species
load('Figure/subCancerX-Ex/Alpha_P_R2.RData')
pval_coef_adj <- t(pval_coef_adj.All[alpha.measure,,]) %>% 
  as.data.frame() %>%
  mutate(R2 = ifelse(R2 < 0, 0, R2),
         direction = case_when(
           is.na(coef) ~ 'none',
           coef > 0    ~ 'positive',
           TRUE        ~ 'negative'
         )) %>%
  dplyr::select(-coef) %>%
  rownames_to_column('var') %>%
  reshape2::melt()

sub.r <- filter(pval_coef_adj, variable == 'R2')
sub.p <- filter(pval_coef_adj, variable == 'P') %>%
  dplyr::select(var, value) %>%
  mutate(P.col = ifelse(value < 0.05, 'sig', 'no')) %>%
  dplyr::select(-value)

sub.rp <- inner_join(sub.p, sub.r, by = 'var') %>%
  mutate(value = value * 100)
sub.rp$var <- factor(sub.rp$var, levels = sub.rp$var[order(sub.rp$value)])

## Beta diversity species
load('Figure/subCancerX-Ex/Beta_P_R2.RData')
R2_P <- data.frame(
  var = names(P.adj[dist.name, ]),
  value = P.adj[dist.name, ],
  R2 = R2.adj[dist.name, ]
) %>%
  mutate(P.col = ifelse(value < 0.05, 'sig', 'no'),
         direction = 'none',
         variable = 'R2',
         value = ifelse(R2 < 0, 0, R2 * 100)) %>%
  dplyr::select(var, value, P.col, direction, variable)

R2_P$var <- factor(R2_P$var, levels = R2_P$var[order(R2_P$value)])

## Combine Species Alpha and Beta results
df3 <- bind_rows(
  R2_P %>% arrange(value) %>% mutate(diversity = 'Beta diversity (species)'),
  sub.rp %>% arrange(value) %>% mutate(diversity = 'Alpha diversity')
)

df3$variable <- gsub("R2", "adjusted", df3$variable)

## Generate plot
cancerclass_div_plot <- ggplot(df3,aes(x =reorder(var,value), y= value)) + 
  geom_segment(aes(x=reorder(var,value), xend=reorder(var,-value), y=0, yend=value), size =4, color = "darkgrey") +
  geom_point(size=5, aes(shape=direction, fill = P.col), stroke=0.5) + 
  scale_fill_manual(values = c('sig'='black','no'= 'white')) +
  scale_shape_manual(values = c('positive'=24,'negative'=25, 'none'= 21)) +
  facet_grid(variable~diversity, scales = 'free')  + 
  theme_classic() +
  labs(y = ('Variability explained (R²), %'), fill = '', x = '', shape ='') +
  coord_flip() +
  theme(panel.grid = element_blank(),
        axis.title = element_text(size=16, color = 'black'),
        legend.text = element_text(size=16, color = 'black'),
        strip.text = element_text(size=16, color = 'black'),
        axis.text.x = element_text(size=16,color = 'black'),
        axis.text.y = element_text(size=16, color = "black")) +
  theme(panel.spacing = unit(2, "lines"))

cancerclass_div_plot

ggsave(file = paste0('Figure/RM_figures/subCancerX vs other cancers.png'), plot=cancerclass_div_plot, width = 12, height = 6)




#-------------------------------------------------------------------------------
## ==== Early onset (Cell: Figure 6BEF) ====

#add separate symbol for species also significant with Age
#can load age species from the Figure results

filename <- 'Result/CancerOnly/DAA/Age/Taxa_ZicoSeq_Species_.Age.adj.Batch.Bristol_score.BMI.Sex.Cancer_class.Metastasis.PPI_day_365.Abx_day_365.Elix_score.Sample_season.Urban.csv'
res_age <- read.csv(filename)

head(res_age,1)

res_age$coef_Age
colnames(res_age)[1] <- 'Species'

dim(res_age) #1237 40

age_species <- res_age$Species[res_age$Qvalue < 0.1]
length(age_species) #543

res_age[grep("s__Clostridium_AP scindens", res_age$Species),]
#coef_Age 0.0009120111 so higher in with age
res_age[grep("s__Veillonella parvula", res_age$Species),]



#load the early onset data
setwd(rd)
q.cutoff <- 0.1
load('EarlyOnset/colorectal/DAA/early_onset/early_onset_ZicoSeq.Rdata')
zico <- as.data.frame(cbind(sign(diff.obj$coef.list$Species[,'early_onsetYes',drop =F]) * diff.obj$R2.list$Species,diff.obj$qv.list$Species))
colnames(zico) <- c('effect','Qvalue')
load('EarlyOnset/colorectal/data.obj.wk.RData')

prop <- (t(t(data.obj.rff$abund.list$Species)/colSums(data.obj.rff$abund.list$Species)))[rownames(zico),]
zico <- cbind.data.frame(prev =rowMeans(prop>0), abund =rowMeans(prop)) %>% merge(zico, by = 0) %>% column_to_rownames('Row.names')



rownames(zico) <- gsub('s__','',rownames(zico))
zico$log10Q <- -log10(zico$Qvalue)
zico$Significant <- zico$Qvalue < q.cutoff
zico$Species <- rownames(zico) 

p1 <- ggplot(zico, aes(x = effect, y = log10Q)) +
  geom_point(aes(size = abund, fill = Significant, color = Significant), 
             shape = 21, alpha = 0.8) +  # Added shape = 21 for circles
  geom_hline(yintercept = -log10(q.cutoff), linetype = "solid", color = "grey") +
  geom_vline(xintercept = 0, linetype = "solid", color = "grey") +
  scale_color_manual(values = c("FALSE" = "grey", "TRUE" = "#CC79A7")) +
  scale_fill_manual(values = c("FALSE" = "grey", "TRUE" = "#CC79A7")) +
  # Remove this line:
  # scale_shape_manual(values = c("FALSE" = 21, "TRUE" = 23)) +  
  scale_size_continuous(range = c(3)) +
  theme_classic(base_size = 14) +
  labs(
    title = "colorectal cancer",
    x = "Effect Size",
    y = expression(-log[10](q)),
    color = "Qvalue < 0.1",
    size = "Mean abundance"
  ) +
  annotate("text", x = max(zico$effect, na.rm = TRUE)+0.05, y = 1.8, 
           label = "higher in early onset →", hjust = 1.1, size = 4, color = "black") +
  annotate("text", x = min(zico$effect, na.rm = TRUE)-0.17, y = 1.8, 
           label = "← lower in early onset", hjust = -0.1, size = 4, color = "black") +
  theme(plot.background = element_rect(fill = "transparent", colour = NA),
        panel.grid = element_blank(),
        axis.title = element_text(size=14, color = 'black'),
        legend.text = element_text(size=14, color = 'black'),
        legend.title = element_text(size=14, color = 'black'),
        strip.text = element_text(size=14, color = 'black'),
        axis.text.x = element_text(size=14,color = 'black'),
        axis.text.y = element_text(size=14, color = 'black'),
        plot.title = element_text(size = 14, hjust = 0.5, color = "black"),
        legend.position = "none") +
  ggrepel::geom_text_repel(
    data = subset(zico, Significant),
    aes(label = paste0("italic('", Species, "')")),
    size = 4,
    box.padding = 0.4,
    max.overlaps = Inf,
    parse = TRUE
  )
p1
p1 + theme(legend.position = "right")






load('EarlyOnset/breast/DAA/early_onset/early_onset_ZicoSeq.Rdata')
zico <- as.data.frame(cbind(sign(diff.obj$coef.list$Species[,'early_onsetYes',drop =F]) * diff.obj$R2.list$Species,diff.obj$qv.list$Species))
colnames(zico) <- c('effect','Qvalue')
load('EarlyOnset/breast/data.obj.wk.RData')

prop <- (t(t(data.obj.rff$abund.list$Species)/colSums(data.obj.rff$abund.list$Species)))[rownames(zico),]
zico <- cbind.data.frame(prev =rowMeans(prop>0), abund =rowMeans(prop)) %>% merge(zico, by = 0) %>% column_to_rownames('Row.names')



rownames(zico) <- gsub('s__','',rownames(zico))
zico$log10Q <- -log10(zico$Qvalue)
zico$Significant <- zico$Qvalue < q.cutoff
zico$Species <- rownames(zico) 
# Remove this line:
# zico$age_shape <- age_shape


#effect size as well as significance cutoff to improve legibility
zico$Significant2 <- FALSE
zico$Significant2[zico$effect > 0.033 & zico$Qvalue < 0.1] <- TRUE
zico$Significant2[zico$effect < -0.03 & zico$Qvalue < 0.1] <- TRUE


p2 <- ggplot(zico, aes(x = effect, y = log10Q)) +
  geom_point(aes(size = abund, fill = Significant, color = Significant), 
             shape = 21, alpha = 0.8) +  # Added shape = 21 for circles
  geom_hline(yintercept = -log10(q.cutoff), linetype = "solid", color = "grey") +
  geom_vline(xintercept = 0, linetype = "solid", color = "grey") +
  scale_color_manual(values = c("FALSE" = "grey", "TRUE" = "#CC79A7")) +
  scale_fill_manual(values = c("FALSE" = "grey", "TRUE" = "#CC79A7")) +
  # Remove this line:
  # scale_shape_manual(values = c("FALSE" = 21, "TRUE" = 23)) +  
  scale_size_continuous(range = c(3)) +
  theme_classic(base_size = 14) +
  labs(
    title = "breast cancer",
    x = "Effect Size",
    y = expression(-log[10](q)),
    color = "Qvalue < 0.1",
    size = "Mean abundance"
  ) +
  annotate("text", x = max(zico$effect, na.rm = TRUE)+0.03, y = 3.2, 
           label = "higher in early onset →", hjust = 1.1, size = 4, color = "black") +
  annotate("text", x = min(zico$effect, na.rm = TRUE)-0.03, y = 3.2, 
           label = "← lower in early onset", hjust = -0.1, size = 4, color = "black") +
  theme(plot.background = element_rect(fill = "transparent", colour = NA),
        panel.grid = element_blank(),
        axis.title = element_text(size=14, color = 'black'),
        legend.text = element_text(size=14, color = 'black'),
        legend.title = element_text(size=14, color = 'black'),
        strip.text = element_text(size=14, color = 'black'),
        axis.text.x = element_text(size=14,color = 'black'),
        axis.text.y = element_text(size=14, color = 'black'),
        plot.title = element_text(size = 14, hjust = 0.5, color = "black"),
        legend.position = "none") +
  ggrepel::geom_text_repel(
    data = subset(zico, Significant2),
    aes(label = paste0("italic('", Species, "')")),
    size = 4,
    box.padding = 0.4,
    max.overlaps = Inf,
    parse = TRUE
  )
p2


ggarrange(p1, p2,
          #labels = c("A.", "B."),      
          ncol = 1, nrow = 2,        
          label.x = 0,            
          label.y = 1,            
          font.label = list(size = 18, face = "bold")  
)

ggsave(file = paste0('../Figure/RM_figures/early_onset_',q.cutoff,'.png'), width = 7, height = 11, bg = "white")


#-------------------------
#check the early onset pathways

q.cutoff <- 0.1

load('EarlyOnset/colorectal/DAA/early_onset/early_onset_ZicoSeq.Rdata')
zico <- as.data.frame(cbind(sign(diff.obj$coef.list$Species[,'early_onsetYes',drop =F]) * diff.obj$R2.list$Species,diff.obj$qv.list$Species))

colnames(zico) <- c('effect','Qvalue')
load('EarlyOnset_func/colorectal/data.obj.wk.RData')

zico$log10Q <- -log10(zico$Qvalue)
zico$Significant <- zico$Qvalue < q.cutoff
zico$pathways <- rownames(zico) 

zico[zico$Qvalue < 0.1,]



load('EarlyOnset_func/breast/DAA/early_onset/early_onset_ZicoSeq.Rdata')

diff.obj$pv.list$pathway["PWY-7754: bile acid 7&alpha;-dehydroxylation" == rownames(diff.obj$pv.list$pathway)]
#0.007

p.adjust(as.numeric(diff.obj$pv.list$pathway), method = "fdr")
sort(p.adjust(as.numeric(diff.obj$pv.list$pathway), method = "fdr"))


#nothing significant for breast
zico <- as.data.frame(cbind(sign(diff.obj$coef.list$pathway[,'early_onsetYes',drop =F]) * diff.obj$R2.list$pathway, 
                            diff.obj$qv.list$pathway, diff.obj$pv.list$pathway))
colnames(zico) <- c('effect','Qvalue', 'Pvalue')
load('EarlyOnset_func/breast/data.obj.wk.RData')


zico$log10Q <- -log10(zico$Qvalue)
zico$Significant <- zico$Qvalue < q.cutoff
zico$pathways <- rownames(zico) 

zico[zico$Qvalue < 0.1,]


rownames(zico)[grep("hydrox", rownames(zico))]

zico["PWY-8134: bile acid 7&beta;-dehydroxylation",]
zico["PWY-7754: bile acid 7&alpha;-dehydroxylation",]

diff.obj$R2.list$pathway["PWY-7754: bile acid 7&alpha;-dehydroxylation",]

#zico["PWY-7754: bile acid 7&alpha;-dehydroxylation",]
#                                                effect   Qvalue   Pvalue    log10Q Significant           
#PWY-7754: bile acid 7&alpha;-dehydroxylation -0.02872478 0.232121  0.007 0.6342855       FALSE 


#-------------------------
#load pathways
#make this plot
clin_meta_breast <- clin_meta[clin_meta$icd10_first_3_name_short == "breast",]
dim(clin_meta_breast)


setwd(wd)
load(file = 'Data/data.obj.pathway.RData') 
path_table <- as.data.frame(data.obj$otu.tab)

seven_alpha <- as.numeric(path_table["PWY-7754: bile acid 7&alpha;-dehydroxylation", clin_meta_breast$BIOM])

clin_meta_breast$seven_alpha <- seven_alpha

sel_colors_plot <- c("#0072B2", "#CC79A7")


#reset the column to "early" and "normal
clin_meta_breast$early_onset2 <- as.character(clin_meta_breast$early_onset)
clin_meta_breast$early_onset2[clin_meta_breast$early_onset2 == "Yes"] <- "early"
clin_meta_breast$early_onset2[clin_meta_breast$early_onset2 == "No"] <- "normal"


set.seed(42)

breast_7alpha_plot <- ggplot(clin_meta_breast, aes(x = early_onset2,
                                         y = seven_alpha,
                                         fill = early_onset2)) +
  geom_boxplot(outlier.shape = NA) +                
  geom_jitter(aes(color = early_onset2),
              width = 0.15, size = 2, alpha = 0.8) + 
  scale_fill_manual(values = sel_colors_plot) +
  scale_color_manual(values = sel_colors_plot) +
  geom_segment(aes(y = 32, yend = 32, x = 1, xend = 2), size = 0.3) +
  geom_text(aes(y = 32*1.03, x = 1.5, label = "p = 0.007"), size = 3, vjust = 0) +
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
    y = "PWY-7754: bile acid 7-alpha-dehydroxylation",
    x = NULL                                         # suppress x-axis title
  )


ggsave("./Figure/RM_figures/boxplot_breast_7alpha.png", breast_7alpha_plot, width = 2.25, height = 4.5, units = "in", dpi = 300)




#-------------------------
###### Generate Supplementary Tables ########
covars <- c("Batch","Bristol_score","BMI", "Age", "Sex","Cancer_class","Metastasis","PPI_day_365", "Abx_day_365","GI_nonGI",
            "Abx_last_month","PPI_last_month","Elix_score","Sample_season","Urban","icd10_first_3_name","Site")

library(arsenal)
setwd(rd)
load('CancerOnly/data.obj.wk.RData')
tmp.meta <- data.obj$meta.dat[,c('icd10_first_3_name_short',covars[!(covars %in% c('Batch','Cancer_class',"icd10_first_3_name"))])]
head(tmp.meta)
tmp.meta$icd10_first_3_name_short <- tolower(gsub("Malignant neoplasm of |Malignant |malignant |Other and | neoplasm| neoplasms", "", tmp.meta$icd10_first_3_name_short ))
idx <- names(which(sort(table(tmp.meta$icd10_first_3_name_short))>15)) 
tmp.meta$icd10_first_3_name_short[which(!(tmp.meta$icd10_first_3_name_short %in% idx))] <- 'Others'
tab <- tableby(as.formula(paste0('icd10_first_3_name_short~.')), data=tmp.meta, digits=2)
summary(tab)
setwd(mfd)
write2html(tab, "Table1_v3.html")




library(readxl)
library(openxlsx)

setwd(wd)
load('Data/data.obj.raw.core.RData')
meta <- data.obj$meta.dat
oldtabS1 <- read_excel(paste0(wd,"/Code/Submission/OLD SUPPLEMENTARY TABLES/Supplementary Table 1.xlsx"), sheet = "MCCM metadata")
colnames(oldtabS1)[!(colnames(oldtabS1) %in% colnames(meta))]
oldtabS1$SampleID[!(oldtabS1$SampleID %in% meta$BIOME_with_sequencing_data)]
colnames(meta)[grep('BIOME_with_sequencing_data', colnames(meta))] <- 'SampleID'
colnames(meta)[grep('icd10_first_3_name_short', colnames(meta))] <- "icd10_first_3_name(name used in manuscript)"
colnames(oldtabS1)[!(colnames(oldtabS1) %in% colnames(meta))]
MCCM_metadata <- meta[,colnames(oldtabS1)] %>% arrange(-desc(Group))




######## Table S1 ######
wb <- loadWorkbook(paste0(wd,"/Code/Submission/Supplementary tables/Supplementary Table 1.xlsx"))

if ("MCCM species abundance" %in% names(wb)) removeWorksheet(wb, "MCCM species abundance")
addWorksheet(wb, "MCCM species abundance")
species <- as.data.frame(data.obj$abund.list$Species) %>% rownames_to_column('Species')
writeData(wb, sheet = "MCCM species abundance", x = species)

if ("MCCM metadata" %in% names(wb)) removeWorksheet(wb, "MCCM metadata")
addWorksheet(wb, "MCCM metadata")
writeData(wb, sheet = "MCCM metadata", x = MCCM_metadata)

saveWorkbook(wb, paste0(wd,"/Code/Submission/Supplementary tables/Supplementary Table 1.xlsx"), overwrite = TRUE)


###### Table S2 ######
wb <- loadWorkbook(paste0(wd,"/Code/Submission/Supplementary tables/Supplementary Table 2.xlsx"))
names(wb)

## TableS2-1[pancancer]
tab.name <- 'TableS2-1'
newtab <- read.csv('Result/PanCancer/DAA/Group/Taxa_ZicoSeq_Species_.Group.adj.BMI.Sex.Age.csv')
newtab <- newtab %>% mutate(EffectSize = sign(coef_GroupCancer)*Func1) %>% dplyr::select(c("X", "Qvalue","EffectSize")) %>% mutate(`Adjust covariates`='BMI,Sex,Age')
colnames(newtab)[1] <- 'Species'
newtab <- inner_join(newtab, as.data.frame(data.obj$otu.name)) 

if (tab.name %in% names(wb)) removeWorksheet(wb, tab.name)
addWorksheet(wb, tab.name)
writeData(wb, sheet = tab.name, x = newtab)

## TableS2-2[cancerx vs all other cancers]
newtabs <- NULL
cancers <- list.dirs('Result/subCancerX-Ex/', full.names = F, recursive = F)
for(cancer in cancers){
  file <- list.files(paste0('Result/subCancerX-Ex/',cancer,'/DAA/icd10_first_3_name_short/'), pattern = 'Taxa_ZicoSeq_Species')
  file <- file[grep('csv$',file)]
  newtab <- read.csv(paste0('Result/subCancerX-Ex/',cancer,'/DAA/icd10_first_3_name_short/',file))
  idx <- colnames(newtab)[grep('^coef_icd10_first_3_name_short',colnames(newtab))]
  adj <- gsub('\\.',',',(gsub('.*adj.|.csv','',file)))
  newtab <- newtab %>% mutate(EffectSize = sign(!!as.name(idx))*Func1) %>% dplyr::select(c("X", "Qvalue","EffectSize")) %>% mutate(`Adjust covariates`=adj)
  colnames(newtab)[1] <- 'Species'
  newtab <- inner_join(newtab, as.data.frame(data.obj$otu.name)) %>% mutate(icd10_first_3_name_short = cancer)
  newtabs <- rbind(newtabs, newtab)
}

tab.name <- 'TableS2-2'
if (tab.name %in% names(wb)) removeWorksheet(wb, tab.name)
addWorksheet(wb, tab.name)
writeData(wb, sheet = tab.name, x = newtabs)


## TableS2-3[cancerx vs all other cancers, pathway]
newtabs <- NULL
cancers <- list.dirs('Result/subCancerX-Ex_func/pathway/', full.names = F, recursive = F)
for(cancer in cancers){
  file <- list.files(paste0('Result/subCancerX-Ex_func/pathway/',cancer,'/DAA/icd10_first_3_name_short/'), pattern = 'Taxa_ZicoSeq_pathway')
  file <- file[grep('csv$',file)]
  newtab <- read.csv(paste0('Result/subCancerX-Ex_func/pathway/',cancer,'/DAA/icd10_first_3_name_short/',file))
  idx <- colnames(newtab)[grep('^coef_icd10_first_3_name_short',colnames(newtab))]
  adj <- gsub('\\.',',',(gsub('.*adj.|.csv','',file)))
  newtab <- newtab %>% mutate(EffectSize = sign(!!as.name(idx))*Func1) %>% dplyr::select(c("X", "Qvalue","EffectSize")) %>% mutate(`Adjust covariates`=adj)
  colnames(newtab)[1] <- 'pathway'
  newtab <- newtab %>% mutate(icd10_first_3_name_short = cancer)
  newtabs <- rbind(newtabs, newtab)
}

tab.name <- 'TableS2-3'
if (tab.name %in% names(wb)) removeWorksheet(wb, tab.name)
addWorksheet(wb, tab.name)
writeData(wb, sheet = tab.name, x = newtabs)


## TableS2-4[Elix_score, Species]
tab.name <- 'TableS2-4'
file <- list.files(paste0('Result/CancerOnly//DAA/Elix_score/'), pattern = 'Taxa_ZicoSeq_Species')
file <- file[grep('csv$',file)]
newtab <- read.csv(paste0('Result/CancerOnly//DAA/Elix_score/',file))
idx <- colnames(newtab)[grep('^coef_Elix_score',colnames(newtab))]
adj <- gsub('\\.',',',(gsub('.*adj.|.csv','',file)))
newtab <- newtab %>% mutate(EffectSize = sign(!!as.name(idx))*Func1) %>% dplyr::select(c("X", "Qvalue","EffectSize")) %>% mutate(`Adjust covariates`=adj)
colnames(newtab)[1] <- 'Species'
newtab <- inner_join(newtab, as.data.frame(data.obj$otu.name)) 

if (tab.name %in% names(wb)) removeWorksheet(wb, tab.name)
addWorksheet(wb, tab.name)
writeData(wb, sheet = tab.name, x = newtab)

## TableS2-5[Elix_score, pathway]
tab.name <- 'TableS2-5'
file <- list.files(paste0('Result/CancerOnly_func/pathway/DAA/Elix_score/'), pattern = 'Taxa_ZicoSeq_pathway')
file <- file[grep('csv$',file)]
newtab <- read.csv(paste0('Result/CancerOnly_func/pathway/DAA/Elix_score/',file))
idx <- colnames(newtab)[grep('^coef_Elix_score',colnames(newtab))]
adj <- gsub('\\.',',',(gsub('.*adj.|.csv','',file)))
newtab <- newtab %>% mutate(EffectSize = sign(!!as.name(idx))*Func1) %>% dplyr::select(c("X", "Qvalue","EffectSize")) %>% mutate(`Adjust covariates`=adj)
colnames(newtab)[1] <- 'pathway'

if (tab.name %in% names(wb)) removeWorksheet(wb, tab.name)
addWorksheet(wb, tab.name)
writeData(wb, sheet = tab.name, x = newtab)




## TableS2-6[Elix components species]
newtabs <- NULL
elixs <- list.dirs('Result/CancerOnly/DAA/', full.names = F, recursive = F)
elixs <- elixs[grepl('Elixhauser_', elixs)]
for(elix in elixs){
  file <- list.files(paste0('Result/CancerOnly/DAA/',elix), pattern = 'Taxa_ZicoSeq_Species')
  file <- file[grep('csv$',file)]
  newtab <- read.csv(paste0('Result/CancerOnly/DAA/',elix,'/',file))
  idx <- colnames(newtab)[grep('^coef_Elixhauser_',colnames(newtab))]
  adj <- gsub('\\.',',',(gsub('.*adj.|.csv','',file)))
  newtab <- newtab %>% mutate(EffectSize = sign(!!as.name(idx))*Func1) %>% dplyr::select(c("X", "Qvalue","EffectSize")) %>% mutate(`Adjust covariates`=adj)
  colnames(newtab)[1] <- 'Species'
  newtab <- inner_join(newtab, as.data.frame(data.obj$otu.name)) %>% mutate(variable = elix)
  newtabs <- rbind(newtabs, newtab)
}

tab.name <- 'TableS2-6'
if (tab.name %in% names(wb)) removeWorksheet(wb, tab.name)
addWorksheet(wb, tab.name)
writeData(wb, sheet = tab.name, x = newtabs)


## TableS2-7[Elix components , pathway]
newtabs <- NULL
elixs <- list.dirs('Result/CancerOnly_func/pathway/DAA/', full.names = F, recursive = F)
elixs <- elixs[grepl('Elixhauser_', elixs)]
for(elix in elixs){
  file <- list.files(paste0('Result/CancerOnly_func/pathway/DAA/',elix), pattern = 'Taxa_ZicoSeq_pathway')
  file <- file[grep('csv$',file)]
  newtab <- read.csv(paste0('Result/CancerOnly_func/pathway/DAA/',elix,'/',file))
  idx <- colnames(newtab)[grep('^coef_Elixhauser_',colnames(newtab))]
  adj <- gsub('\\.',',',(gsub('.*adj.|.csv','',file)))
  newtab <- newtab %>% mutate(EffectSize = sign(!!as.name(idx))*Func1) %>% dplyr::select(c("X", "Qvalue","EffectSize")) %>% mutate(`Adjust covariates`=adj)
  colnames(newtab)[1] <- 'Pathway'
  newtab <- newtab %>% mutate(variable = elix)
  newtabs <- rbind(newtabs, newtab)
}

tab.name <- 'TableS2-7'
if (tab.name %in% names(wb)) removeWorksheet(wb, tab.name)
addWorksheet(wb, tab.name)
writeData(wb, sheet = tab.name, x = newtabs)


## TableS2-8[cancerx vs control, species]
newtabs <- NULL
cancers <- list.dirs('Result/subCancerX_Control/', full.names = F, recursive = F)
for(cancer in cancers){
  file <- list.files(paste0('Result/subCancerX_Control/',cancer,'/DAA/icd10_first_3_name_short/'), pattern = 'Taxa_ZicoSeq_Species')
  file <- file[grep('csv$',file)]
  newtab <- read.csv(paste0('Result/subCancerX_Control/',cancer,'/DAA/icd10_first_3_name_short/',file))
  idx <- colnames(newtab)[grep('^coef_icd10_first_3_name_short',colnames(newtab))]
  adj <- gsub('\\.',',',(gsub('.*adj.|.csv','',file)))
  newtab <- newtab %>% mutate(EffectSize = sign(!!as.name(idx))*Func1) %>% dplyr::select(c("X", "Qvalue","EffectSize")) %>% mutate(`Adjust covariates`=adj)
  colnames(newtab)[1] <- 'Species'
  newtab <- inner_join(newtab, as.data.frame(data.obj$otu.name)) %>% mutate(variable = cancer)
  newtabs <- rbind(newtabs, newtab)
}

tab.name <- 'TableS2-8'
if (tab.name %in% names(wb)) removeWorksheet(wb, tab.name)
addWorksheet(wb, tab.name)
writeData(wb, sheet = tab.name, x = newtabs)


## TableS2-9[cancerx vs control, pathway]
newtabs <- NULL
cancers <- list.dirs('Result/subCancerX_Control_func/pathway/', full.names = F, recursive = F)
for(cancer in cancers){
  file <- list.files(paste0('Result/subCancerX_Control_func/pathway/',cancer,'/DAA/icd10_first_3_name_short/'), pattern = 'Taxa_ZicoSeq_pathway')
  file <- file[grep('csv$',file)]
  newtab <- read.csv(paste0('Result/subCancerX_Control_func/pathway/',cancer,'/DAA/icd10_first_3_name_short/',file))
  idx <- colnames(newtab)[grep('^coef_icd10_first_3_name_short',colnames(newtab))]
  adj <- gsub('\\.',',',(gsub('.*adj.|.csv','',file)))
  newtab <- newtab %>% mutate(EffectSize = sign(!!as.name(idx))*Func1) %>% dplyr::select(c("X", "Qvalue","EffectSize")) %>% mutate(`Adjust covariates`=adj)
  colnames(newtab)[1] <- 'pathway'
  newtab <- newtab %>% mutate(icd10_first_3_name_short = cancer)
  newtabs <- rbind(newtabs, newtab)
}

tab.name <- 'TableS2-9'
if (tab.name %in% names(wb)) removeWorksheet(wb, tab.name)
addWorksheet(wb, tab.name)
writeData(wb, sheet = tab.name, x = newtabs)



## TableS2-10[pancancer func]
tab.name <- 'TableS2-10'
newtab <- read.csv('Result/PanCancer_func/DAA/Group/Taxa_ZicoSeq_pathway_.Group.adj.BMI.Sex.Age.csv')
newtab <- newtab %>% mutate(EffectSize = sign(coef_GroupCancer)*Func1) %>% dplyr::select(c("X", "Qvalue","EffectSize")) %>% mutate(`Adjust covariates`='BMI,Sex,Age')
colnames(newtab)[1] <- 'pathway'

if (tab.name %in% names(wb)) removeWorksheet(wb, tab.name)
addWorksheet(wb, tab.name)
writeData(wb, sheet = tab.name, x = newtab)

saveWorkbook(wb, paste0(wd,"/Code/Submission/Supplementary tables/Supplementary Table 2.xlsx"), overwrite = TRUE)


###### Table S3 ######
setwd(wd)

wb <- loadWorkbook(paste0(wd,"/Code/Submission/OLD SUPPLEMENTARY TABLES/Supplementary Table 3.xlsx"))
names(wb)
## TableS3-1[demo tech variables]
newtabs <- NULL
vars <- list.dirs('Result/CancerOnly/DAA/', full.names = F, recursive = F)
vars <- c("Abx_day_365", "Abx_last_month", "Age", "BMI","Bristol_score","Elix_score","Metastasis","Sample_season","Sex","Urban","GI_nonGI")
for(var in vars){
  file <- list.files(paste0('Result/CancerOnly/DAA/',var), pattern = 'Taxa_ZicoSeq_Species')
  file <- file[grep('csv$',file)]
  newtab <- read.csv(paste0('Result/CancerOnly/DAA/',var,'/',file))
  idx <- colnames(newtab)[grep(paste0('^coef_',var),colnames(newtab))]
  adj <- gsub('\\.',',',(gsub('.*adj.|.csv','',file)))
  newtab <- newtab %>% mutate(EffectSize = sign(!!as.name(idx))*Func1) %>% dplyr::select(c("X", "Qvalue","EffectSize")) %>% mutate(`Adjust covariates`=adj)
  colnames(newtab)[1] <- 'Species'
  newtab <- inner_join(newtab, as.data.frame(data.obj$otu.name)) %>% mutate(variable = var)
  newtabs <- rbind(newtabs, newtab)
}

tab.name <- 'TableS3-1'
if (tab.name %in% names(wb)) removeWorksheet(wb, tab.name)
addWorksheet(wb, tab.name)
writeData(wb, sheet = tab.name, x = newtabs)


## TableS3-2[demo tech variables, func]
newtabs <- NULL
vars <- list.dirs('Result/CancerOnly_func/pathway/DAA/', full.names = F, recursive = F)
vars <- c("Abx_day_365", "Abx_last_month", "Age", "BMI","Bristol_score","Elix_score","Metastasis","Sample_season","Sex","Urban","GI_nonGI")
for(var in vars){
  file <- list.files(paste0('Result/CancerOnly_func/pathway/DAA/',var), pattern = 'Taxa_ZicoSeq_pathway')
  file <- file[grep('csv$',file)]
  newtab <- read.csv(paste0('Result/CancerOnly_func/pathway/DAA/',var,'/',file))
  idx <- colnames(newtab)[grep(paste0('^coef_',var),colnames(newtab))]
  adj <- gsub('\\.',',',(gsub('.*adj.|.csv','',file)))
  newtab <- newtab %>% mutate(EffectSize = sign(!!as.name(idx))*Func1) %>% dplyr::select(c("X", "Qvalue","EffectSize")) %>% mutate(`Adjust covariates`=adj)
  colnames(newtab)[1] <- 'Pathway'
  newtab <- newtab %>% mutate(variable = var)
  newtabs <- rbind(newtabs, newtab)
}

tab.name <- 'TableS3-2'
if (tab.name %in% names(wb)) removeWorksheet(wb, tab.name)
addWorksheet(wb, tab.name)
writeData(wb, sheet = tab.name, x = newtabs)

saveWorkbook(wb, paste0(wd,"/Code/Submission/Supplementary tables/Supplementary Table 3.xlsx"), overwrite = TRUE)




###### Table S5 ######
setwd(wd)
wb <- loadWorkbook(paste0(wd,"/Code/Submission/OLD SUPPLEMENTARY TABLES/Supplementary Table 5.xlsx"))
names(wb)

tab.name <- "s_CRC early onset" 
if (tab.name %in% names(wb)) removeWorksheet(wb, tab.name)
file <- list.files(paste0('Result/EarlyOnset/colorectal/DAA/early_onset/'), pattern = 'Taxa_ZicoSeq_Species')
file <- file[grep('csv$',file)]
newtab <- read.csv(paste0('Result/EarlyOnset/colorectal/DAA/early_onset/',file))
idx <- colnames(newtab)[grep('^coef_early_onset',colnames(newtab))]
adj <- gsub('\\.',',',(gsub('.*adj.|.csv','',file)))
newtab <- newtab %>% mutate(EffectSize = sign(!!as.name(idx))*Func1) %>% dplyr::select(c("X", "Qvalue","EffectSize")) %>% mutate(`Adjust covariates`=adj)
colnames(newtab)[1] <- 'Species'
addWorksheet(wb, tab.name)
writeData(wb, sheet = tab.name, x = newtab)


tab.name <- "s_breast early onset" 
if (tab.name %in% names(wb)) removeWorksheet(wb, tab.name)
file <- list.files(paste0('Result/EarlyOnset/breast/DAA/early_onset/'), pattern = 'Taxa_ZicoSeq_Species')
file <- file[grep('csv$',file)]
newtab <- read.csv(paste0('Result/EarlyOnset/breast/DAA/early_onset/',file))
idx <- colnames(newtab)[grep('^coef_early_onset',colnames(newtab))]
adj <- gsub('\\.',',',(gsub('.*adj.|.csv','',file)))
newtab <- newtab %>% mutate(EffectSize = sign(!!as.name(idx))*Func1) %>% dplyr::select(c("X", "Qvalue","EffectSize")) %>% mutate(`Adjust covariates`=adj)
colnames(newtab)[1] <- 'Species'
addWorksheet(wb, tab.name)
writeData(wb, sheet = tab.name, x = newtab)

tab.name <- "s_brain early onset" 
if (tab.name %in% names(wb)) removeWorksheet(wb, tab.name)
file <- list.files(paste0('Result/EarlyOnset/brain/DAA/early_onset/'), pattern = 'Taxa_ZicoSeq_Species')
file <- file[grep('csv$',file)]
newtab <- read.csv(paste0('Result/EarlyOnset/brain/DAA/early_onset/',file))
idx <- colnames(newtab)[grep('^coef_early_onset',colnames(newtab))]
adj <- gsub('\\.',',',(gsub('.*adj.|.csv','',file)))
newtab <- newtab %>% mutate(EffectSize = sign(!!as.name(idx))*Func1) %>% dplyr::select(c("X", "Qvalue","EffectSize")) %>% mutate(`Adjust covariates`=adj)
colnames(newtab)[1] <- 'Species'
addWorksheet(wb, tab.name)
writeData(wb, sheet = tab.name, x = newtab)

tab.name <- "path_breast early onset" 
if (tab.name %in% names(wb)) removeWorksheet(wb, tab.name)
file <- list.files(paste0('Result/EarlyOnset_func/breast/DAA/early_onset/'), pattern = 'Taxa_ZicoSeq_pathway')
file <- file[grep('csv$',file)]
newtab <- read.csv(paste0('Result/EarlyOnset_func/breast/DAA/early_onset/',file))
idx <- colnames(newtab)[grep('^coef_early_onset',colnames(newtab))]
adj <- gsub('\\.',',',(gsub('.*adj.|.csv','',file)))
newtab <- newtab %>% mutate(EffectSize = sign(!!as.name(idx))*Func1) %>% dplyr::select(c("X", "Qvalue","EffectSize")) %>% mutate(`Adjust covariates`=adj)
colnames(newtab)[1] <- 'Pathway'
addWorksheet(wb, tab.name)
writeData(wb, sheet = tab.name, x = newtab)

tab.name <- "path_CRC early onset" 
if (tab.name %in% names(wb)) removeWorksheet(wb, tab.name)
file <- list.files(paste0('Result/EarlyOnset_func/colorectal/DAA/early_onset/'), pattern = 'Taxa_ZicoSeq_pathway')
file <- file[grep('csv$',file)]
newtab <- read.csv(paste0('Result/EarlyOnset_func/colorectal/DAA/early_onset/',file))
idx <- colnames(newtab)[grep('^coef_early_onset',colnames(newtab))]
adj <- gsub('\\.',',',(gsub('.*adj.|.csv','',file)))
newtab <- newtab %>% mutate(EffectSize = sign(!!as.name(idx))*Func1) %>% dplyr::select(c("X", "Qvalue","EffectSize")) %>% mutate(`Adjust covariates`=adj)
colnames(newtab)[1] <- 'Pathway'
addWorksheet(wb, tab.name)
writeData(wb, sheet = tab.name, x = newtab)

tab.name <- "path_brain early onset" 
if (tab.name %in% names(wb)) removeWorksheet(wb, tab.name)
file <- list.files(paste0('Result/EarlyOnset_func/brain/DAA/early_onset/'), pattern = 'Taxa_ZicoSeq_pathway')
file <- file[grep('csv$',file)]
newtab <- read.csv(paste0('Result/EarlyOnset_func/brain/DAA/early_onset/',file))
idx <- colnames(newtab)[grep('^coef_early_onset',colnames(newtab))]
adj <- gsub('\\.',',',(gsub('.*adj.|.csv','',file)))
newtab <- newtab %>% mutate(EffectSize = sign(!!as.name(idx))*Func1) %>% dplyr::select(c("X", "Qvalue","EffectSize")) %>% mutate(`Adjust covariates`=adj)
colnames(newtab)[1] <- 'Pathway'
addWorksheet(wb, tab.name)
writeData(wb, sheet = tab.name, x = newtab)

saveWorkbook(wb, paste0(wd,"/Code/Submission/Supplementary tables/Supplementary Table 5.xlsx"), overwrite = TRUE)



## 08/27/2025 change order of tables 
path <- paste0(wd, "/Code/Submission/Supplementary tables/Supplementary Table 2.xlsx")
wb <- loadWorkbook(path)

df_TableS28 <- read.xlsx(wb, sheet = "TableS2-8")
df_TableS29 <- read.xlsx(wb, sheet = "TableS2-9")
df_TableS22 <- read.xlsx(wb, sheet = "TableS2-2")
df_TableS23 <- read.xlsx(wb, sheet = "TableS2-3")

for (s in c("TableS2-2","TableS2-3","TableS2-8","TableS2-9")) {
  if (s %in% names(wb)) removeWorksheet(wb, s)
}

addWorksheet(wb, "TableS2-2") # will receive content from TableS2-8
addWorksheet(wb, "TableS2-3")  # will receive content from TableS2-9
addWorksheet(wb, "TableS2-8")  # will receive content from TableS2-2
addWorksheet(wb, "TableS2-9") # will receive content from TableS2-3

writeData(wb, sheet = "TableS2-2", x = df_TableS28)
writeData(wb, sheet = "TableS2-3", x = df_TableS29)
writeData(wb, sheet = "TableS2-8", x = df_TableS22)
writeData(wb, sheet = "TableS2-9", x = df_TableS23)

saveWorkbook(wb, path, overwrite = TRUE)






### ===== Check early onset breast beta ====
file_dir <- dirname(rstudioapi::getSourceEditorContext()$path)
file_dir <- dirname(dirname(dirname(dirname(file_dir))))

wd <- file_dir
rd <- paste0(file_dir,'/Result/')
setwd(rd)
tm <- load('EarlyOnset/breast/Beta/early_onset/R2_pvalue.RData')
pv.adj.mat




### === check ECI component =====
setwd("/Users/luyang1/myicloud/Documents/Mayo_project/2023_01_09_Oncobiome/mforge_clean/Result/CancerOnly/DAA/")
dirs <- list.dirs(full.names = F)
dirs <- dirs[grep('^Elixhauser',dirs)]
for(dir in dirs){
  load(paste0(dir,'/',dir,'_ZicoSeq.Rdata'))
  x <- colnames(diff.obj$coef.list$Species)[grep(dir,colnames(diff.obj$coef.list$Species))]
  cat(dir,': ',x,'\n')
}

setwd("/Users/luyang1/myicloud/Documents/Mayo_project/2023_01_09_Oncobiome/mforge_clean/Result/CancerOnly_func/pathway/DAA/")
dirs <- list.dirs(full.names = F)
dirs <- dirs[grep('^Elixhauser',dirs)]
for(dir in dirs){
  load(paste0(dir,'/',dir,'_ZicoSeq.Rdata'))
  x <- colnames(diff.obj$coef.list$pathway)[grep(dir,colnames(diff.obj$coef.list$pathway))]
  cat(dir,': ',x,'\n')
}



