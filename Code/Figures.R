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
library(RColorBrewer)
library(ape)


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
pval <- as.data.frame(pval) %>% 
  mutate(P=p.adjust(P, 'BH'))  %>% # if using qval
  mutate(sig = ifelse(P<0.1, '*',''))

#Shannon
alpha$icd10_first_3_name <- as.factor(alpha$icd10_first_3_name_short)
rownames(pval) <- gsub("_", " ", rownames(pval))

ord <- as.vector((aggregate(Shannon ~ icd10_first_3_name_short, data=alpha, function(x) median(x)) %>% arrange(Shannon))[,1])
alpha <- within(alpha, icd10_first_3_name_short <- factor(icd10_first_3_name_short,levels =ord))
Means <- alpha %>% group_by(icd10_first_3_name_short) %>% dplyr::summarize(Avg = median(Shannon))

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


alpha$cancer_healthy <- ifelse(alpha$cancer_healthy == "cancer", "Cancer", "Healthy")

alpha$cancer_healthy <- factor(alpha$cancer_healthy,
                               levels = c("Healthy", "Cancer"))

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
  keepR2 <- rownames(R2[topQ,i, drop =F] %>% dplyr::arrange(-abs(!!as.name(i))))[1:round(length(topQ)*0.25)]
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
tree.sub <- ape::drop.tip(data.obj$tree, tips)
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
max_plots <- 40  

for (ECI in ECIs) {
  
  idx <- Q.All$pathway$pathway[which(Q.All$pathway[, ECI] < 0.1049)]
  cat(ECI, '-', length(idx), '\n')
  
  if (length(table(data.obj$meta.dat[, ECI])) == 2 & length(idx) > 0) {
    
    idx_chunks <- split(idx, ceiling(seq_along(idx)/max_plots))
    
    chunk_num <- 1
    for (idx_chunk in idx_chunks) {
      
      pdf_file <- paste0(fd, "RM_figures/All_Significant_pathway_boxplot_", 
                         ECI, "_part", chunk_num, ".pdf")
      pdf(pdf_file, width = 5, height = 3.5)
      
      for (i in idx_chunk) {
        
        plot_df <- sqrt(t(prop[i, , drop = FALSE])) %>%
          merge(data.obj$meta.dat[, ECI, drop = FALSE], by = 0)
        colnames(plot_df)[2] <- "taxon"
        plot_df[, ECI] <- as.factor(plot_df[, ECI])
        
        p_plot <- ggboxplot(
          plot_df, 
          x = ECI, 
          y = "taxon",
          color = ECI,
          palette = brewer.pal(5, "Dark2"), 
          width = 0.8,
          add = "jitter", add.params = list(alpha = 0.5, size = 1.5),
          ylab = "sqrt(reads per million, CPM)"
        ) +
          scale_y_continuous(trans = "sqrt") +
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
      
      dev.off()
      chunk_num <- chunk_num + 1
    }
  }
}




## ==== Figure Alpha Beta main variables bar plot (Cell: Figure S3C)======

#changed this to fewer variables to make consistent with the main text

#ECI score, smoking status and
#Batch,Bristol_score,BMI,Age,Sex,Cancer_class,Metastasis,PPI_day_365,Abx_day_365,Sample_season,Urban

#vars.incl <- c("Bristol_score","BMI", "Age", "Sex", "Metastasis","PPI_day_365", "Abx_day_365",
#               "Abx_last_month","PPI_last_month","Elix_score","Sample_season","Urban","icd10_first_3_name_short","Site","smoking_category","prior_chemotherapy")
vars.incl <- c("Bristol_score","BMI", "Age", "Sex", "Metastasis","PPI_day_365", "Abx_day_365",
               "Abx_last_month","PPI_last_month","Elix_score","Sample_season","Urban","icd10_first_3_name_short",
               "smoking_category")


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
  mutate(P.col = ifelse(value<0.05,'sig','no')) %>% 
  dplyr::rename(P=value)  # add 11/07/2025
# dplyr::select(-value) # remove 11/07/2025
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
  mutate(P.col = ifelse(value<0.05,'sig','no')) %>% 
  dplyr::rename(P=value)  # add 11/07/2025
# dplyr::select(-value) # remove 11/07/2025
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
  dplyr::rename(P = value) %>% ## add 11/07/2025
  # dplyr::select(-value) %>% ## remove 11/07/2025
  mutate(direction = 'none', variable = 'R2') %>% 
  dplyr::rename(value = R2) %>%
  mutate(value=ifelse(value <0,0, value *100))
ord <- (R2_P[order(R2_P$value),'var'])
R2_P2 <- cbind.data.frame(value = c(P.unadj[dist.name,]), R2 = c(R2.unadj[dist.name,])) %>% 
  rownames_to_column('var') %>% dplyr::filter(var !='Batch') %>% 
  mutate(P.col = ifelse(value<0.05, 'sig','no'))  %>% 
  dplyr::rename(P = value) %>% ## add 11/07/2025
  # dplyr::select(-value) %>% ## remove 11/07/2025
  mutate(direction = 'none', variable = 'R2') %>% 
  dplyr::rename(value = R2) %>%
  mutate(value=ifelse(value <0,0, value *100))
R2_P12 <- rbind(R2_P %>% mutate(grp = 'adjusted'),R2_P2 %>% mutate(grp = 'marginal'))
R2_P12 <- within(R2_P12, var <- factor(var, levels=ord))


xxx <- P.adj[dist.name,]

# beta - pathway
setwd(wd)
dist.name <- 'BC'
load('Figure/CancerOnly_func/pathway/Beta_P_R2.RData')
R2_P <- cbind.data.frame(value = c(P.adj[dist.name,]), R2 = c(R2.adj[dist.name,])) %>% 
  rownames_to_column('var') %>% dplyr::filter(var !='Batch') %>% 
  mutate(P.col = ifelse(value<0.05, 'sig','no'))  %>% 
  dplyr::rename(P = value) %>% ## add 11/07/2025
  # dplyr::select(-value) %>% ## remove 11/07/2025
  mutate(direction = 'none', variable = 'R2') %>% dplyr::rename(value = R2) %>%
  mutate(value=ifelse(value <0,0, value *100))
ord <- (R2_P[order(R2_P$value),'var'])
R2_P2 <- cbind.data.frame(value = c(P.unadj[dist.name,]), R2 = c(R2.unadj[dist.name,])) %>% 
  rownames_to_column('var') %>% dplyr::filter(var !='Batch') %>% 
  mutate(P.col = ifelse(value<0.05, 'sig','no'))  %>% 
  dplyr::rename(P = value) %>% ## add 11/07/2025
  # dplyr::select(-value) %>% ## remove 11/07/2025
  mutate(direction = 'none', variable = 'R2') %>% dplyr::rename(value = R2) %>%
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
df2$var[df2$var =="Bristol_score"] <- "Bristol Stool Form Scale"
df2$var[df2$var =="PPI_day_365"] <- "PPI in past year (No)"
df2$var[df2$var =="PPI_last_month"] <- "PPI in past month (No)"
df2$var[df2$var =="Abx_day_365"] <- "Antibiotics in past year (No)"
df2$var[df2$var =="Abx_last_month"] <- "Antibiotics in past month (No)"
df2$var[df2$var =="Elix_score"] <- "Elixhauser Comorbidity score"
df2$var[df2$var =="Sample_season"] <- "Sample season"
df2$var[df2$var =="Urban"] <- "Residence type (non-urban)"
df2$var[df2$var =="Sex"] <- "Sex (Female)"
df2$var[df2$var =="Metastasis"] <- "Metastasis (No)"
df2$var[df2$var =="smoking_category"] <- "Smoking category"
#df2$var[df2$var =="prior_chemotherapy"] <- "Prior cancer treatment (No)"


#plotting settings start here

df3 <- df2
df3 <- df3[!df3$diversity == "Beta diversity (pathway)",]
df3 <- df3[df3$grp == "adjusted",]

#unique(df3$var)
#variables_to_plot <- c("Cancer class" , "Bristol stool score", "Elixhauser Comorbidity score", "Sex (Female)", "Age",
#                       "Antibiotics in past year (No)", "Antibiotics in past month (No)", "Prior cancer treatment (No)","BMI", 
#                       "Smoking Category", "PPI in past year (No)", 
#                       "Residence type (non-urban)", "Site", "PPI in past month (No)", "Sample season", 
#                       "Metastasis (No)") 
variables_to_plot <- c("Cancer class" , "Bristol Stool Form Scale", "Elixhauser Comorbidity score", "Sex (Female)", "Age",
                       "Antibiotics in past year (No)", "Antibiotics in past month (No)","BMI", 
                       "Smoking Category", "PPI in past year (No)", 
                       "Residence type (non-urban)", "PPI in past month (No)", "Sample season", 
                       "Metastasis (No)") 

sel_colors <- palette.colors(palette = "Okabe-Ito")[2:8]

variable_colors_to_plot <- rep(sel_colors[6], length(variables_to_plot))
names(variable_colors_to_plot) <- variables_to_plot

variable_colors_to_plot[variables_to_plot %in% c("Bristol Stool Form Scale", "Sample season")] <- sel_colors[2] #technical
#variable_colors_to_plot[variables_to_plot %in% c("Sex (Female)", "Age", "Residence type (non-urban)", "BMI", "Site")] <- sel_colors[1] #demographic
variable_colors_to_plot[variables_to_plot %in% c("Sex (Female)", "Age", "Residence type (non-urban)", "BMI")] <- sel_colors[1] #demographic


main_div_plot <- 
  df3 %>%
  group_by(diversity) %>% mutate(P = p.adjust(P, method = "BH")) %>% mutate(P.col = ifelse(P < 0.1, "sig", "no")) %>% ungroup() %>% ## add 11/07/2025
  ggplot(aes(x = reorder(var, value), y = value)) +
  geom_segment(aes(xend = reorder(var, value), y = 0,  yend = value),
               size = 4, colour = "darkgrey") +
  geom_point(aes(shape = direction, fill = P.col), size = 5, stroke = 0.5) +
  scale_fill_manual(values = c(sig = "black", no = "white")) +
  scale_shape_manual(values = c(positive = 24, negative = 25, none = 21)) +
  facet_grid(grp ~ diversity, scales = "free") +
  coord_flip() +                   
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
  )  + guides(
    fill = guide_legend(title = "", order = 1,
                        override.aes = list(shape = 21)),
    shape = guide_legend(title = "", order = 2,
                         override.aes = list(fill = "white"))
  )


main_div_plot



ggsave(file = 'Figure/RM_figures/Figure1_main.svg', plot=main_div_plot, width = 12, height = 5.3)
ggsave(file = 'Figure/RM_figures/Figure1_main.png', plot=main_div_plot, width = 12, height = 5.3)



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

xx <- df.tmp_2 %>%  mutate(q = p.adjust(P, 'BH')) 
xx[xx$q<0.1 & xx$P>0.05,]
xx[xx$q>0.1 & xx$P<0.05,]

elix_div_plot <- df.tmp_2 %>% 
  mutate(q = p.adjust(P, 'BH')) %>% mutate(P.col = ifelse(q<0.1, 'sig','no')) %>% # add 11/07/2025
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
  mutate(value=p.adjust(value,'BH')) %>% # add 11/07/2025
  mutate(P.col = ifelse((value) < 0.05, 'sig', 'no')) %>%
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
  mutate(value=p.adjust(value,'BH')) %>% # add 11/07/2025
  mutate(P.col = ifelse(value < 0.1, 'sig', 'no'),
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
  theme(panel.spacing = unit(2, "lines")) + guides(
    fill = guide_legend(title = "", order = 1,
                        override.aes = list(shape = 21)),
    shape = guide_legend(title = "", order = 2,
                         override.aes = list(fill = "white"))
  )

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
            "Abx_last_month","PPI_last_month","Elix_score","Sample_season","Urban","icd10_first_3_name","Site","smoking_category","prior_chemotherapy")

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

######## Table S1 (Ruben added some avariables 04/14/2026)######
setwd(wd)
load('Data/data.obj.raw.core.RData')
Supplementary_Table_1 <- read_excel("Code/Submission/Supplementary tables/old/Supplementary Table 1_RM.xlsx", sheet = "TableS1-2")
Supplementary_Table_1$SampleID[!(Supplementary_Table_1$SampleID %in% rownames(data.obj$meta.dat))]
Supplementary_Table_1 <- data.obj$meta.dat %>% 
  rownames_to_column('SampleID') %>% 
  dplyr::select(SampleID, sample) %>%
  inner_join(Supplementary_Table_1) %>% column_to_rownames('...1') %>% arrange(Group)
## consider remove above steps before we do submission
wb <- loadWorkbook(paste0(wd,"/Code/Submission/Supplementary tables/old/Supplementary Table 1_RM.xlsx"))
tab.name <- "TableS1-2" # some new treatment variables were added by Ruben, we need to use his version, and add SampleID=BIOM_S*
if (tab.name %in% names(wb)) removeWorksheet(wb, tab.name)
addWorksheet(wb, tab.name)
writeData(wb, sheet = tab.name, x = Supplementary_Table_1)
current.order <- names(wb)
new.order <- c(current.order[current.order != tab.name][1:2],  
               tab.name,                                         
               current.order[current.order != tab.name][-(1:2)] 
)
worksheetOrder(wb) <- match(new.order, names(wb))
saveWorkbook(wb, paste0(wd,"/Code/Submission/Supplementary tables/Supplementary Table 1.xlsx"), overwrite = TRUE)


###### Table S2 ######
wb <- loadWorkbook(paste0(wd,"/Code/Submission/Supplementary tables/Supplementary Table 2.xlsx"))
names(wb)

## TableS2-1[pancancer-Species]
tab.name <- 'TableS2-1'
newtab <- read.csv('Result/PanCancer/DAA/Group/Taxa_ZicoSeq_Species_.Group.adj.BMI.Sex.Age.csv')
newtab <- newtab %>% mutate(EffectSize = sign(coef_GroupCancer)*Func1) %>% dplyr::select(c("X", "Qvalue","EffectSize")) %>% mutate(`Adjust covariates`='BMI,Sex,Age') 
colnames(newtab)[1] <- 'Species'
newtab <- inner_join(newtab, as.data.frame(data.obj$otu.name)) %>% mutate(Comparison='cancer patients vs healthy controls (reference)')

if (tab.name %in% names(wb)) removeWorksheet(wb, tab.name)
addWorksheet(wb, tab.name)
writeData(wb, sheet = tab.name, x = newtab)

## TableS2-2[pancancer func]
tab.name <- 'TableS2-2'
newtab <- read.csv('Result/PanCancer_func/DAA/Group/Taxa_ZicoSeq_pathway_.Group.adj.BMI.Sex.Age.csv')
newtab <- newtab %>% mutate(EffectSize = sign(coef_GroupCancer)*Func1) %>% dplyr::select(c("X", "Qvalue","EffectSize"))  %>% mutate(`Adjust covariates`='BMI,Sex,Age', Comparison='cancer patients vs healthy controls (reference)') 
colnames(newtab)[1] <- 'pathway'

if (tab.name %in% names(wb)) removeWorksheet(wb, tab.name)
addWorksheet(wb, tab.name)
writeData(wb, sheet = tab.name, x = newtab)



# setwd(rd)
# getwd()
# setwd('subCancerX_Control/')
# dirs <- list.dirs(full.names = F, recursive = F)
# for(dir in dirs){
#   setwd(paste0(rd,'subCancerX_Control/',dir,'/DAA/icd10_first_3_name_short/'))
#   file <- list.files()
#   file <- file[grepl('Taxa_ZicoSeq_Species_.icd10_first_3_name_short.adj.BMI.Sex.Age.csv',file)]
#   xx <- read.csv(file)
#   yy <- colnames(xx)[grep("coef_icd10_first_3_name_",colnames(xx))]
#   cat(dir,'-',yy,'\n')
# }
# 
# setwd(rd)
# getwd()
# setwd('subCancerX_Control_func/pathway/')
# dirs <- list.dirs(full.names = F, recursive = F)
# for(dir in dirs){
#   setwd(paste0(rd,'subCancerX_Control_func/pathway/',dir,'/DAA/icd10_first_3_name_short/'))
#   file <- list.files()
#   file <- file[grepl('Taxa_ZicoSeq_pathway_.icd10_first_3_name_short.adj.BMI.Sex.Age.csv',file)]
#   xx <- read.csv(file)
#   yy <- colnames(xx)[grep("coef_icd10_first_3_name_",colnames(xx))]
#   cat(dir,'-',yy,'\n')
# }


## TableS2-3[cancerx vs control, species]
setwd(wd)
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
  newtab <- inner_join(newtab, as.data.frame(data.obj$otu.name)) %>% mutate(Comparison = paste0(gsub('Control-','',cancer), '- healthy controls(reference)'))
  newtabs <- rbind(newtabs, newtab)
}

tab.name <- 'TableS2-3'
if (tab.name %in% names(wb)) removeWorksheet(wb, tab.name)
addWorksheet(wb, tab.name)
writeData(wb, sheet = tab.name, x = newtabs)


## TableS2-4[cancerx vs control, pathway]
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
  newtab <- newtab %>% mutate(Comparison = paste0(cancer,'- healthy controls(reference)'))
  newtabs <- rbind(newtabs, newtab)
}

tab.name <- 'TableS2-4'
if (tab.name %in% names(wb)) removeWorksheet(wb, tab.name)
addWorksheet(wb, tab.name)
writeData(wb, sheet = tab.name, x = newtabs)


## TableS2-5[cancerx vs all other cancers]
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
  newtab <- inner_join(newtab, as.data.frame(data.obj$otu.name)) %>% mutate(Comparison = paste0(cancer,'- all other cancers (reference)'))
  newtabs <- rbind(newtabs, newtab)
}

tab.name <- 'TableS2-5'
if (tab.name %in% names(wb)) removeWorksheet(wb, tab.name)
addWorksheet(wb, tab.name)
writeData(wb, sheet = tab.name, x = newtabs)


## TableS2-6[cancerx vs all other cancers, pathway]
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
  newtab <- newtab %>% mutate(Comparison = paste0(cancer,'- all other cancers (reference)'))
  newtabs <- rbind(newtabs, newtab)
}

tab.name <- 'TableS2-6'
if (tab.name %in% names(wb)) removeWorksheet(wb, tab.name)
addWorksheet(wb, tab.name)
writeData(wb, sheet = tab.name, x = newtabs)


## TableS2-7[Elix_score, Species]
tab.name <- 'TableS2-7'
file <- list.files(paste0('Result/CancerOnly//DAA/Elix_score/'), pattern = 'Taxa_ZicoSeq_Species')
file <- file[grep('csv$',file)]
newtab <- read.csv(paste0('Result/CancerOnly//DAA/Elix_score/',file))
idx <- colnames(newtab)[grep('^coef_Elix_score',colnames(newtab))]
adj <- gsub('\\.',',',(gsub('.*adj.|.csv','',file)))
newtab <- newtab %>% mutate(EffectSize = sign(!!as.name(idx))*Func1) %>% dplyr::select(c("X", "Qvalue","EffectSize")) %>% mutate(`Adjust covariates`=adj)
colnames(newtab)[1] <- 'Species'
newtab <- inner_join(newtab, as.data.frame(data.obj$otu.name)) %>% mutate(MainVariableOfInterets="Elix_score")

if (tab.name %in% names(wb)) removeWorksheet(wb, tab.name)
addWorksheet(wb, tab.name)
writeData(wb, sheet = tab.name, x = newtab)

## TableS2-8[Elix_score, pathway]
tab.name <- 'TableS2-8'
file <- list.files(paste0('Result/CancerOnly_func/pathway/DAA/Elix_score/'), pattern = 'Taxa_ZicoSeq_pathway')
file <- file[grep('csv$',file)]
newtab <- read.csv(paste0('Result/CancerOnly_func/pathway/DAA/Elix_score/',file))
idx <- colnames(newtab)[grep('^coef_Elix_score',colnames(newtab))]
adj <- gsub('\\.',',',(gsub('.*adj.|.csv','',file)))
newtab <- newtab %>% mutate(EffectSize = sign(!!as.name(idx))*Func1) %>% dplyr::select(c("X", "Qvalue","EffectSize")) %>% mutate(`Adjust covariates`=adj, MainVariableOfInterets="Elix_score")
colnames(newtab)[1] <- 'pathway'

if (tab.name %in% names(wb)) removeWorksheet(wb, tab.name)
addWorksheet(wb, tab.name)
writeData(wb, sheet = tab.name, x = newtab)




## TableS2-9[Elix components species]
newtabs <- NULL
elixs <- list.dirs('Result/CancerOnly/DAA/', full.names = F, recursive = F)
elixs <- elixs[grepl('Elixhauser_', elixs)]
for(elix in elixs){
  file <- list.files(paste0('Result/CancerOnly/DAA/',elix), pattern = 'Taxa_ZicoSeq_Species')
  file <- file[grep('csv$',file)]
  newtab <- read.csv(paste0('Result/CancerOnly/DAA/',elix,'/',file))
  idx <- colnames(newtab)[grep('^coef_Elixhauser_',colnames(newtab))]
  cat(idx,'\n')
  adj <- gsub('\\.',',',(gsub('.*adj.|.csv','',file)))
  newtab <- newtab %>% mutate(EffectSize = sign(!!as.name(idx))*Func1) %>% dplyr::select(c("X", "Qvalue","EffectSize")) %>% mutate(`Adjust covariates`=adj)
  colnames(newtab)[1] <- 'Species'
  newtab <- inner_join(newtab, as.data.frame(data.obj$otu.name)) %>% mutate(Comparison = paste0(elix, ': Yes vs No(reference)'))
  newtabs <- rbind(newtabs, newtab)
}

tab.name <- 'TableS2-9'
if (tab.name %in% names(wb)) removeWorksheet(wb, tab.name)
addWorksheet(wb, tab.name)
writeData(wb, sheet = tab.name, x = newtabs)


## TableS2-10[Elix components , pathway]
newtabs <- NULL
elixs <- list.dirs('Result/CancerOnly_func/pathway/DAA/', full.names = F, recursive = F)
elixs <- elixs[grepl('Elixhauser_', elixs)]
for(elix in elixs){
  file <- list.files(paste0('Result/CancerOnly_func/pathway/DAA/',elix), pattern = 'Taxa_ZicoSeq_pathway')
  file <- file[grep('csv$',file)]
  newtab <- read.csv(paste0('Result/CancerOnly_func/pathway/DAA/',elix,'/',file))
  idx <- colnames(newtab)[grep('^coef_Elixhauser_',colnames(newtab))]
  cat(idx,'\n')
  adj <- gsub('\\.',',',(gsub('.*adj.|.csv','',file)))
  newtab <- newtab %>% mutate(EffectSize = sign(!!as.name(idx))*Func1) %>% dplyr::select(c("X", "Qvalue","EffectSize")) %>% mutate(`Adjust covariates`=adj)
  colnames(newtab)[1] <- 'Pathway'
  newtab <- newtab %>% mutate(Comparison = paste0(elix, ': Yes vs No(reference)'))
  newtabs <- rbind(newtabs, newtab)
}

tab.name <- 'TableS2-10'
if (tab.name %in% names(wb)) removeWorksheet(wb, tab.name)
addWorksheet(wb, tab.name)
writeData(wb, sheet = tab.name, x = newtabs)

saveWorkbook(wb, paste0(wd,"/Code/Submission/Supplementary tables/Supplementary Table 2.xlsx"), overwrite = TRUE)


###### Table S3 ######
file_dir <- dirname(rstudioapi::getSourceEditorContext()$path)
file_dir <- dirname(dirname(dirname(dirname(file_dir))))
setwd(file_dir)
load('Data/data.obj.mph.RData')

wb <- loadWorkbook(paste0(file_dir,"/Code/Submission/Supplementary tables/20260419/Supplementary Table 2_metaphlan.xlsx"))
names(wb)

removeWorksheet(wb, 'TableS2-1')
removeWorksheet(wb, 'TableS2-2')
removeWorksheet(wb, 'TableS2-3')
removeWorksheet(wb, 'TableS2-4')
removeWorksheet(wb, 'TableS2-5')
removeWorksheet(wb, 'TableS2-6')
removeWorksheet(wb, 'TableS2-7')
removeWorksheet(wb, 'TableS2-8')
removeWorksheet(wb, 'TableS2-9')
removeWorksheet(wb, 'TableS2-10')



## contents
tab.name <- 'contents'
newtabs <- readWorkbook(wb, sheet = "contents")
newtabs <- newtabs[grep('Species',newtabs[,2]),]
newtabs$Sheet.name <- paste0('TableS3-',1:5)
# newtabs$`Description..Pvalue,.Qvalue,.and.effect.size.were.calculated.using.ZicoSeq..The.effect.size.represents.R2.×.sign(coefficient).` <- gsub("MetaCyc pathways","Genus",newtabs$`Description..Pvalue,.Qvalue,.and.effect.size.were.calculated.using.ZicoSeq..The.effect.size.represents.R2.×.sign(coefficient).`)
if (tab.name %in% names(wb)) removeWorksheet(wb, tab.name)
addWorksheet(wb, tab.name)
writeData(wb, sheet = tab.name, x = newtabs)



## TableS3-1[pancancer-Species]
tab.name <- 'TableS3-1'
newtab <- read.csv('Result_mph/PanCancer/DAA/Group/Taxa_ZicoSeq_Species_.Group.adj.BMI.Sex.Age.csv')
newtab <- newtab %>% mutate(EffectSize = sign(coef_GroupCancer)*Func1) %>% dplyr::select(c("X","Qvalue","EffectSize")) %>% mutate(`Adjust covariates`='BMI,Sex,Age') 
colnames(newtab)[1] <- 'Species'
newtab <- inner_join(newtab, as.data.frame(data.obj$otu.name)) %>% mutate(Comparison='cancer patients vs healthy controls (reference)')
if (tab.name %in% names(wb)) removeWorksheet(wb, tab.name)
addWorksheet(wb, tab.name)
writeData(wb, sheet = tab.name, x = newtab)


# ## TableS2-Genus[pancancer-Genus]
# tab.name <- 'TableS2-2'
# newtab <- read.csv('Result_mph/PanCancer/DAA/Group/Taxa_ZicoSeq_Genus_.Group.adj.BMI.Sex.Age.csv')
# newtab <- newtab %>% mutate(EffectSize = sign(coef_GroupCancer)*Func1) %>% dplyr::select(c("X", "Pvalue","Qvalue","EffectSize")) %>% mutate(`Adjust covariates`='BMI,Sex,Age') 
# colnames(newtab)[1] <- 'Genus'
# newtab <- inner_join(newtab, unique(as.data.frame(data.obj$otu.name) %>% dplyr::select(-c('Species')) %>% remove_rownames)) %>% mutate(Comparison='cancer patients vs healthy controls (reference)')
# 
# if (tab.name %in% names(wb)) removeWorksheet(wb, tab.name)
# addWorksheet(wb, tab.name)
# writeData(wb, sheet = tab.name, x = newtab)



## TableS3-3[cancerx vs control, species]-
## final TableS3-2
setwd(wd)
newtabs <- NULL
cancers <- list.dirs('Result_mph/subCancerX_Control/', full.names = F, recursive = F)
for(cancer in cancers){
  file <- list.files(paste0('Result_mph/subCancerX_Control/',cancer,'/DAA/icd10_first_3_name_short/'), pattern = 'Taxa_ZicoSeq_Species')
  file <- file[grep('csv$',file)]
  newtab <- read.csv(paste0('Result_mph/subCancerX_Control/',cancer,'/DAA/icd10_first_3_name_short/',file))
  idx <- colnames(newtab)[grep('^coef_icd10_first_3_name_short',colnames(newtab))]
  adj <- gsub('\\.',',',(gsub('.*adj.|.csv','',file)))
  newtab <- newtab %>% mutate(EffectSize = sign(!!as.name(idx))*Func1) %>% dplyr::select(c("X", "Qvalue","EffectSize")) %>% mutate(`Adjust covariates`=adj)
  colnames(newtab)[1] <- 'Species'
  newtab <- inner_join(newtab, as.data.frame(data.obj$otu.name)) %>% mutate(Comparison = paste0(gsub('Control-','',cancer), '- healthy controls(reference)'))
  newtabs <- rbind(newtabs, newtab)
}

tab.name <- 'TableS3-2'
if (tab.name %in% names(wb)) removeWorksheet(wb, tab.name)
addWorksheet(wb, tab.name)
writeData(wb, sheet = tab.name, x = newtabs)

# ## TableS2-4[cancerx vs control, Genus]
# setwd(wd)
# newtabs <- NULL
# cancers <- list.dirs('Result_mph/subCancerX_Control/', full.names = F, recursive = F)
# for(cancer in cancers){
#   file <- list.files(paste0('Result_mph/subCancerX_Control/',cancer,'/DAA/icd10_first_3_name_short/'), pattern = 'Taxa_ZicoSeq_Genus')
#   file <- file[grep('csv$',file)]
#   newtab <- read.csv(paste0('Result_mph/subCancerX_Control/',cancer,'/DAA/icd10_first_3_name_short/',file))
#   idx <- colnames(newtab)[grep('^coef_icd10_first_3_name_short',colnames(newtab))]
#   adj <- gsub('\\.',',',(gsub('.*adj.|.csv','',file)))
#   newtab <- newtab %>% mutate(EffectSize = sign(!!as.name(idx))*Func1) %>% dplyr::select(c("X", "Pvalue","Qvalue","EffectSize")) %>% mutate(`Adjust covariates`=adj)
#   colnames(newtab)[1] <- 'Genus'
#   newtab <- inner_join(newtab, unique(as.data.frame(data.obj$otu.name) %>% dplyr::select(-c('Species')) %>% remove_rownames)) %>% mutate(Comparison = paste0(gsub('Control-','',cancer), '- healthy controls(reference)'))
#   newtabs <- rbind(newtabs, newtab)
# }
# 
# tab.name <- 'TableS2-4'
# if (tab.name %in% names(wb)) removeWorksheet(wb, tab.name)
# addWorksheet(wb, tab.name)
# writeData(wb, sheet = tab.name, x = newtabs)



## TableS2-5[cancerx vs all other cancers]
## final TableS3-3
newtabs <- NULL
cancers <- list.dirs('Result_mph/subCancerX-Ex/', full.names = F, recursive = F)
for(cancer in cancers){
  file <- list.files(paste0('Result_mph/subCancerX-Ex/',cancer,'/DAA/icd10_first_3_name_short/'), pattern = 'Taxa_ZicoSeq_Species')
  file <- file[grep('csv$',file)]
  newtab <- read.csv(paste0('Result_mph/subCancerX-Ex/',cancer,'/DAA/icd10_first_3_name_short/',file))
  idx <- colnames(newtab)[grep('^coef_icd10_first_3_name_short',colnames(newtab))]
  adj <- gsub('\\.',',',(gsub('.*adj.|.csv','',file)))
  newtab <- newtab %>% mutate(EffectSize = sign(!!as.name(idx))*Func1) %>% dplyr::select(c("X", "Qvalue","EffectSize")) %>% mutate(`Adjust covariates`=adj)
  colnames(newtab)[1] <- 'Species'
  newtab <- inner_join(newtab, as.data.frame(data.obj$otu.name)) %>% mutate(Comparison = paste0(cancer,'- all other cancers (reference)'))
  newtabs <- rbind(newtabs, newtab)
}

tab.name <- 'TableS3-3'
if (tab.name %in% names(wb)) removeWorksheet(wb, tab.name)
addWorksheet(wb, tab.name)
writeData(wb, sheet = tab.name, x = newtabs)


# ## TableS2-6[cancerx vs all other cancers]
# newtabs <- NULL
# cancers <- list.dirs('Result_mph/subCancerX-Ex/', full.names = F, recursive = F)
# for(cancer in cancers){
#   file <- list.files(paste0('Result_mph/subCancerX-Ex/',cancer,'/DAA/icd10_first_3_name_short/'), pattern = 'Taxa_ZicoSeq_Genus')
#   file <- file[grep('csv$',file)]
#   newtab <- read.csv(paste0('Result_mph/subCancerX-Ex/',cancer,'/DAA/icd10_first_3_name_short/',file))
#   idx <- colnames(newtab)[grep('^coef_icd10_first_3_name_short',colnames(newtab))]
#   adj <- gsub('\\.',',',(gsub('.*adj.|.csv','',file)))
#   newtab <- newtab %>% mutate(EffectSize = sign(!!as.name(idx))*Func1) %>% dplyr::select(c("X", "Pvalue","Qvalue","EffectSize")) %>% mutate(`Adjust covariates`=adj)
#   colnames(newtab)[1] <- 'Genus'
#   newtab <- inner_join(newtab, unique(as.data.frame(data.obj$otu.name) %>% dplyr::select(-c('Species')) %>% remove_rownames)) %>% mutate(Comparison = paste0(cancer,'- all other cancers (reference)'))
#   newtabs <- rbind(newtabs, newtab)
# }
# 
# tab.name <- 'TableS2-6'
# if (tab.name %in% names(wb)) removeWorksheet(wb, tab.name)
# addWorksheet(wb, tab.name)
# writeData(wb, sheet = tab.name, x = newtabs)


## TableS2-7[Elix_score, Species]
## TableS3-4
tab.name <- 'TableS3-4'
file <- list.files(paste0('Result_mph/CancerOnly//DAA/Elix_score/'), pattern = 'Taxa_ZicoSeq_Species')
file <- file[grep('csv$',file)]
newtab <- read.csv(paste0('Result_mph/CancerOnly//DAA/Elix_score/',file))
idx <- colnames(newtab)[grep('^coef_Elix_score',colnames(newtab))]
adj <- gsub('\\.',',',(gsub('.*adj.|.csv','',file)))
newtab <- newtab %>% mutate(EffectSize = sign(!!as.name(idx))*Func1) %>% dplyr::select(c("X","Qvalue","EffectSize")) %>% mutate(`Adjust covariates`=adj)
colnames(newtab)[1] <- 'Species'
newtab <- inner_join(newtab, as.data.frame(data.obj$otu.name)) %>% mutate(MainVariableOfInterets="Elix_score")
if (tab.name %in% names(wb)) removeWorksheet(wb, tab.name)
addWorksheet(wb, tab.name)
writeData(wb, sheet = tab.name, x = newtab)


# ## TableS2-8[Elix_score, Genus]
# tab.name <- 'TableS2-8'
# file <- list.files(paste0('Result_mph/CancerOnly//DAA/Elix_score/'), pattern = 'Taxa_ZicoSeq_Genus')
# file <- file[grep('csv$',file)]
# newtab <- read.csv(paste0('Result_mph/CancerOnly//DAA/Elix_score/',file))
# idx <- colnames(newtab)[grep('^coef_Elix_score',colnames(newtab))]
# adj <- gsub('\\.',',',(gsub('.*adj.|.csv','',file)))
# newtab <- newtab %>% mutate(EffectSize = sign(!!as.name(idx))*Func1) %>% dplyr::select(c("X","Pvalue", "Qvalue","EffectSize")) %>% mutate(`Adjust covariates`=adj)
# colnames(newtab)[1] <- 'Genus'
# newtab <- inner_join(newtab, unique(as.data.frame(data.obj$otu.name) %>% dplyr::select(-c('Species')) %>% remove_rownames)) %>% mutate(MainVariableOfInterets="Elix_score")
# 
# if (tab.name %in% names(wb)) removeWorksheet(wb, tab.name)
# addWorksheet(wb, tab.name)
# writeData(wb, sheet = tab.name, x = newtab)


## TableS2-9[Elix components species]
## TableS3-5
newtabs <- NULL
elixs <- list.dirs('Result_mph/CancerOnly/DAA/', full.names = F, recursive = F)
elixs <- elixs[grepl('Elixhauser_', elixs)]
for(elix in elixs){
  file <- list.files(paste0('Result_mph/CancerOnly/DAA/',elix), pattern = 'Taxa_ZicoSeq_Species')
  file <- file[grep('csv$',file)]
  newtab <- read.csv(paste0('Result_mph/CancerOnly/DAA/',elix,'/',file))
  idx <- colnames(newtab)[grep('^coef_Elixhauser_',colnames(newtab))]
  cat(idx,'\n')
  adj <- gsub('\\.',',',(gsub('.*adj.|.csv','',file)))
  newtab <- newtab %>% mutate(EffectSize = sign(!!as.name(idx))*Func1) %>% dplyr::select(c("X", "Qvalue","EffectSize")) %>% mutate(`Adjust covariates`=adj)
  colnames(newtab)[1] <- 'Species'
  newtab <- inner_join(newtab, as.data.frame(data.obj$otu.name)) %>% mutate(Comparison = paste0(elix, ': Yes vs No(reference)'))
  newtabs <- rbind(newtabs, newtab)
}

tab.name <- 'TableS3-5'
if (tab.name %in% names(wb)) removeWorksheet(wb, tab.name)
addWorksheet(wb, tab.name)
writeData(wb, sheet = tab.name, x = newtabs)

# ## TableS2-10[Elix components Genus]
# newtabs <- NULL
# elixs <- list.dirs('Result_mph/CancerOnly/DAA/', full.names = F, recursive = F)
# elixs <- elixs[grepl('Elixhauser_', elixs)]
# for(elix in elixs){
#   file <- list.files(paste0('Result_mph/CancerOnly/DAA/',elix), pattern = 'Taxa_ZicoSeq_Genus')
#   file <- file[grep('csv$',file)]
#   newtab <- read.csv(paste0('Result_mph/CancerOnly/DAA/',elix,'/',file))
#   idx <- colnames(newtab)[grep('^coef_Elixhauser_',colnames(newtab))]
#   cat(idx,'\n')
#   adj <- gsub('\\.',',',(gsub('.*adj.|.csv','',file)))
#   newtab <- newtab %>% mutate(EffectSize = sign(!!as.name(idx))*Func1) %>% dplyr::select(c("X", "Pvalue","Qvalue","EffectSize")) %>% mutate(`Adjust covariates`=adj)
#   colnames(newtab)[1] <- 'Genus'
#   newtab <- inner_join(newtab, unique(as.data.frame(data.obj$otu.name) %>% dplyr::select(-c('Species')) %>% remove_rownames)) %>% mutate(Comparison = paste0(elix, ': Yes vs No(reference)'))
#   newtabs <- rbind(newtabs, newtab)
# }
# 
# tab.name <- 'TableS2-10'
# if (tab.name %in% names(wb)) removeWorksheet(wb, tab.name)
# addWorksheet(wb, tab.name)
# writeData(wb, sheet = tab.name, x = newtabs)

saveWorkbook(wb, paste0(wd,"/Code/Submission/Supplementary tables/Supplementary Table 3.xlsx"), overwrite = TRUE)

###### Table S3->S5 ######
setwd(wd)
load('Data/data.obj.raw.core.RData')

wb <- loadWorkbook(paste0(wd,"/Code/Submission/OLD SUPPLEMENTARY TABLES/Supplementary Table 3.xlsx"))
names(wb)
removeWorksheet(wb, 'TableS3-1')
removeWorksheet(wb, 'TableS3-2')
## TableS5-1[demo tech variables]
newtabs <- NULL
vars <- list.dirs('Result/CancerOnly/DAA/', full.names = F, recursive = F)
vars <- c("Abx_day_365", "Abx_last_month","PPI_day_365", "PPI_last_month", "Age", "BMI","Bristol_score","Elix_score","Metastasis","Sample_season","Sex","Urban","GI_nonGI","smoking_category","prior_chemotherapy")
for(var in vars){
  file <- list.files(paste0('Result/CancerOnly/DAA/',var), pattern = 'Taxa_ZicoSeq_Species')
  file <- file[grep('csv$',file)]
  newtab <- read.csv(paste0('Result/CancerOnly/DAA/',var,'/',file))
  idx <- colnames(newtab)[grep(paste0('^coef_',var),colnames(newtab))]
  adj <- gsub('\\.',',',(gsub('.*adj.|.csv','',file)))
  if(var %in% c("Sample_season","smoking_category")){
    newtab <- newtab %>% mutate(EffectSize = Func1) %>% dplyr::select(c("X", "Qvalue","EffectSize")) %>% mutate(`Adjust covariates`=adj)
  }else{
    newtab <- newtab %>% mutate(EffectSize = sign(!!as.name(idx))*Func1) %>% dplyr::select(c("X", "Qvalue","EffectSize")) %>% mutate(`Adjust covariates`=adj)
  }
  colnames(newtab)[1] <- 'Species'
  newtab <- inner_join(newtab, as.data.frame(data.obj$otu.name)) %>% mutate(variable = var)
  if(var=="Abx_day_365") newtab$variable="Antibiotics in past year"
  if(var=="Abx_last_month") newtab$variable="Antibiotics in past month"
  if(var=="PPI_day_365") newtab$variable="PPI in past year"
  if(var=="PPI_last_month") newtab$variable="PPI in past month"
  if(var=="Bristol_score") newtab$variable='Bristol stool score'
  if(var=="Elix_score") newtab$variable='Elixhauser Comorbidity score'
  if(var=="Sample_season") newtab$variable='Sample season'
  if(var=="smoking_category") newtab$variable='Smoking Category'
  if(var=="prior_chemotherapy") newtab$variable='Prior cancer treatment'
  newtabs <- rbind(newtabs, newtab)
}

tab.name <- 'TableS5-1'
if (tab.name %in% names(wb)) removeWorksheet(wb, tab.name)
addWorksheet(wb, tab.name)
writeData(wb, sheet = tab.name, x = newtabs)


## TableS5-2[demo tech variables, func]
newtabs <- NULL
vars <- list.dirs('Result/CancerOnly_func/pathway/DAA/', full.names = F, recursive = F)
vars <- c("Abx_day_365", "Abx_last_month", "PPI_day_365", "PPI_last_month","Age", "BMI","Bristol_score","Elix_score","Metastasis","Sample_season","Sex","Urban","GI_nonGI","smoking_category","prior_chemotherapy")
for(var in vars){
  file <- list.files(paste0('Result/CancerOnly_func/pathway/DAA/',var), pattern = 'Taxa_ZicoSeq_pathway')
  file <- file[grep('csv$',file)]
  newtab <- read.csv(paste0('Result/CancerOnly_func/pathway/DAA/',var,'/',file))
  idx <- colnames(newtab)[grep(paste0('^coef_',var),colnames(newtab))]
  adj <- gsub('\\.',',',(gsub('.*adj.|.csv','',file)))
  if(var %in% c("Sample_season","smoking_category")){
    newtab <- newtab %>% mutate(EffectSize = Func1) %>% dplyr::select(c("X", "Qvalue","EffectSize")) %>% mutate(`Adjust covariates`=adj)
  }else{
    newtab <- newtab %>% mutate(EffectSize = sign(!!as.name(idx))*Func1) %>% dplyr::select(c("X", "Qvalue","EffectSize")) %>% mutate(`Adjust covariates`=adj)
  }
  colnames(newtab)[1] <- 'Pathway'
  newtab <- newtab %>% mutate(variable = var)
  if(var=="Abx_day_365") newtab$variable="Antibiotics in past year"
  if(var=="Abx_last_month") newtab$variable="Antibiotics in past month"
  if(var=="PPI_day_365") newtab$variable="PPI in past year"
  if(var=="PPI_last_month") newtab$variable="PPI in past month"
  if(var=="Bristol_score") newtab$variable='Bristol stool score'
  if(var=="Elix_score") newtab$variable='Elixhauser Comorbidity score'
  if(var=="Sample_season") newtab$variable='Sample season'
  if(var=="smoking_category") newtab$variable='Smoking Category'
  if(var=="prior_chemotherapy") newtab$variable='Prior cancer treatment'
  
  newtabs <- rbind(newtabs, newtab)
}

tab.name <- 'TableS5-2'
if (tab.name %in% names(wb)) removeWorksheet(wb, tab.name)
addWorksheet(wb, tab.name)
writeData(wb, sheet = tab.name, x = newtabs)

saveWorkbook(wb, paste0(wd,"/Code/Submission/Supplementary tables/Supplementary Table 5.xlsx"), overwrite = TRUE)




###### Table S5-> S8 ######
setwd(wd)
wb <- loadWorkbook(paste0(wd,"/Code/Submission/Supplementary tables/20260419/Supplementary Table 5.xlsx"))
names(wb)

removeWorksheet(wb, 'TableS5-1')
removeWorksheet(wb, 'TableS5-2')
removeWorksheet(wb, 'TableS5-3')
removeWorksheet(wb, 'TableS5-4')
removeWorksheet(wb, 'TableS5-5')
removeWorksheet(wb, 'TableS5-6')
removeWorksheet(wb, 'TableS5-7')
removeWorksheet(wb, 'TableS5-8')

tab.name <- "TableS8-1" 
if (tab.name %in% names(wb)) removeWorksheet(wb, tab.name)
file <- list.files(paste0('Result/EarlyOnset/colorectal/DAA/early_onset/'), pattern = 'Taxa_ZicoSeq_Species')
file <- file[grep('csv$',file)]
newtab <- read.csv(paste0('Result/EarlyOnset/colorectal/DAA/early_onset/',file))
idx <- colnames(newtab)[grep('^coef_early_onset',colnames(newtab))]
cat(idx)
adj <- gsub('\\.',',',(gsub('.*adj.|.csv','',file)))
newtab <- newtab %>% mutate(EffectSize = sign(!!as.name(idx))*Func1) %>% dplyr::select(c("X", "Qvalue","EffectSize")) %>%  
  left_join(as.data.frame(data.obj$otu.name) %>% dplyr::rename(X=Species)) %>% dplyr::rename(Species ='X') %>% 
  mutate(`Adjust covariates`=adj, Comparison='CRC early onset: Yes vs No(reference)')
colnames(newtab)[1] <- 'Species'
addWorksheet(wb, tab.name)
writeData(wb, sheet = tab.name, x = newtab)



tab.name <- "TableS8-2" 
if (tab.name %in% names(wb)) removeWorksheet(wb, tab.name)
file <- list.files(paste0('Result/EarlyOnset_func/colorectal/DAA/early_onset/'), pattern = 'Taxa_ZicoSeq_pathway')
file <- file[grep('csv$',file)]
newtab <- read.csv(paste0('Result/EarlyOnset_func/colorectal/DAA/early_onset/',file))
idx <- colnames(newtab)[grep('^coef_early_onset',colnames(newtab))]
adj <- gsub('\\.',',',(gsub('.*adj.|.csv','',file)))
newtab <- newtab %>% mutate(EffectSize = sign(!!as.name(idx))*Func1) %>% dplyr::select(c("X", "Qvalue","EffectSize")) %>% mutate(`Adjust covariates`=adj, Comparison='CRC early onset: Yes vs No(reference)')
colnames(newtab)[1] <- 'Pathway'
addWorksheet(wb, tab.name)
writeData(wb, sheet = tab.name, x = newtab)


tab.name <- "TableS8-3" 
if (tab.name %in% names(wb)) removeWorksheet(wb, tab.name)
file <- list.files(paste0('Result/EarlyOnset/breast/DAA/early_onset/'), pattern = 'Taxa_ZicoSeq_Species')
file <- file[grep('csv$',file)]
newtab <- read.csv(paste0('Result/EarlyOnset/breast/DAA/early_onset/',file))
idx <- colnames(newtab)[grep('^coef_early_onset',colnames(newtab))]
adj <- gsub('\\.',',',(gsub('.*adj.|.csv','',file)))
newtab <- newtab %>% mutate(EffectSize = sign(!!as.name(idx))*Func1) %>% dplyr::select(c("X", "Qvalue","EffectSize")) %>% 
  left_join(as.data.frame(data.obj$otu.name) %>% dplyr::rename(X=Species)) %>% dplyr::rename(Species ='X') %>% 
  mutate(`Adjust covariates`=adj, Comparison='Breast cancer early onset: Yes vs No(reference)')
colnames(newtab)[1] <- 'Species'
addWorksheet(wb, tab.name)
writeData(wb, sheet = tab.name, x = newtab)


tab.name <- "TableS8-4" 
if (tab.name %in% names(wb)) removeWorksheet(wb, tab.name)
file <- list.files(paste0('Result/EarlyOnset_func/breast/DAA/early_onset/'), pattern = 'Taxa_ZicoSeq_pathway')
file <- file[grep('csv$',file)]
newtab <- read.csv(paste0('Result/EarlyOnset_func/breast/DAA/early_onset/',file))
idx <- colnames(newtab)[grep('^coef_early_onset',colnames(newtab))]
adj <- gsub('\\.',',',(gsub('.*adj.|.csv','',file)))
newtab <- newtab %>% mutate(EffectSize = sign(!!as.name(idx))*Func1) %>% dplyr::select(c("X", "Qvalue","EffectSize")) %>% mutate(`Adjust covariates`=adj, Comparison='Breast cancer early onset: Yes vs No(reference)')
colnames(newtab)[1] <- 'Pathway'
addWorksheet(wb, tab.name)
writeData(wb, sheet = tab.name, x = newtab)


tab.name <- "TableS8-5" 
if (tab.name %in% names(wb)) removeWorksheet(wb, tab.name)
file <- list.files(paste0('Result/EarlyOnset/brain/DAA/early_onset/'), pattern = 'Taxa_ZicoSeq_Species')
file <- file[grep('csv$',file)]
newtab <- read.csv(paste0('Result/EarlyOnset/brain/DAA/early_onset/',file))
idx <- colnames(newtab)[grep('^coef_early_onset',colnames(newtab))]
adj <- gsub('\\.',',',(gsub('.*adj.|.csv','',file)))
newtab <- newtab %>% mutate(EffectSize = sign(!!as.name(idx))*Func1) %>% dplyr::select(c("X", "Qvalue","EffectSize")) %>% 
  left_join(as.data.frame(data.obj$otu.name) %>% dplyr::rename(X=Species)) %>% dplyr::rename(Species ='X') %>% 
  mutate(`Adjust covariates`=adj, Comparison='Brain cancer early onset: Yes vs No(reference)')
colnames(newtab)[1] <- 'Species'
addWorksheet(wb, tab.name)
writeData(wb, sheet = tab.name, x = newtab)



tab.name <- "TableS8-6" 
if (tab.name %in% names(wb)) removeWorksheet(wb, tab.name)
file <- list.files(paste0('Result/EarlyOnset_func/brain/DAA/early_onset/'), pattern = 'Taxa_ZicoSeq_pathway')
file <- file[grep('csv$',file)]
newtab <- read.csv(paste0('Result/EarlyOnset_func/brain/DAA/early_onset/',file))
idx <- colnames(newtab)[grep('^coef_early_onset',colnames(newtab))]
adj <- gsub('\\.',',',(gsub('.*adj.|.csv','',file)))
newtab <- newtab %>% mutate(EffectSize = sign(!!as.name(idx))*Func1) %>% dplyr::select(c("X", "Qvalue","EffectSize")) %>% mutate(`Adjust covariates`=adj, Comparison='Brain cancer early onset: Yes vs No(reference)')
colnames(newtab)[1] <- 'Pathway'
addWorksheet(wb, tab.name)
writeData(wb, sheet = tab.name, x = newtab)

saveWorkbook(wb, paste0(wd,"/Code/Submission/Supplementary tables/Supplementary Table 8.xlsx"), overwrite = TRUE)






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



### ==== Check Elix_score distribution ======
file_dir <- dirname(rstudioapi::getSourceEditorContext()$path)
file_dir <- dirname(dirname(dirname(dirname(file_dir))))
setwd(file_dir)
load('Data/data.obj.raw.core.RData')
dim(data.obj$meta.dat)
ggplot(data.obj$meta.dat, aes(x = Elix_score, fill = Group)) +
  geom_histogram(aes(y = after_stat(density)),position = "identity", alpha = 0.7,binwidth = 1,color = "white") +
  scale_fill_brewer(palette = 'Set1') + 
  labs(x = "Elix score",title = "", fill = "") +
  theme_classic() + 
  theme(axis.text = element_text(color = "black", size = 14),
        legend.text = element_text(color = "black", size = 14),
        legend.title = element_text(color = "black", size = 14),
        axis.title = element_text(color = "black", size = 14))
ggsave(paste0(fd,'RM_Figures/Elix_score_hist.pdf'), width = 6, height = 4)

summary(data.obj$meta.dat[data.obj$meta.dat$Group=='Control','Elix_score'])

library(OptimalCutpoints)
result <- optimal.cutpoints(X = "Elix_score", status = "Group", tag.healthy = "Control", 
                            methods = "MaxKappa",data = data.obj$meta.dat)
summary(result)    

## ====== Response to reviewer2's subset analysis: Check HV vs Low/High ECI analysis results======
file_dir <- dirname(rstudioapi::getSourceEditorContext()$path)
file_dir <- dirname(dirname(dirname(dirname(file_dir))))
setwd(file_dir)

rm(list = ls())

load('Data/data.obj.raw.core.RData')

cutoff = 0.1

## ECI low <=1 ; higH >=4
tm = load('Result/ECI_subset/Control-LowCancer/DAA/ECI2/ECI2_ZicoSeq.Rdata')
low_q <- diff.obj$qv.list$Species
low_FC <- diff.obj$R2.list$Species * sign(diff.obj$coef.list$Species[,'ECI2LowCancer'])
sum(low_q<cutoff)
df <- merge(diff.obj$pv.list$Species, low_q, by = 0)%>% column_to_rownames('Row.names') %>% merge(low_FC, by = 0) %>% dplyr::rename(Species = 'Row.names', `Effect size`=Func1)
df <- df %>% inner_join(as.data.frame(data.obj$otu.name))

df2 <- df[df$Qvalue<0.1,]
length(grep('Streptococcus',df2$Species))
df2$`Effect size`[grep('Streptococcus',df2$Species)]
length(grep('Bifidobacterium ',df2$Species))
df2$`Effect size`[grep('Bifidobacterium',df2$Species)]
length(grep('Blautia',df2$Species))
df2$`Effect size`[grep('Blautia',df2$Species)]
length(grep('Bacteroides',df2$Species))
df2$`Effect size`[grep('Bacteroides',df2$Species)]
length(grep('Faecalibacterium',df2$Species))
df2$`Effect size`[grep('Faecalibacterium',df2$Species)]


tm = load('Result/ECI_subset0/Control-HighCancer/DAA/ECI22/ECI22_ZicoSeq.Rdata') # ECI high in ECI22 and ECI2 are the same, just permutation makes little difference, in order for consistency, use this 
high_q <- diff.obj$qv.list$Species
high_FC <- diff.obj$R2.list$Species * sign(diff.obj$coef.list$Species[,'ECI22HighCancer'])
sum(high_q<cutoff)

tm <- load('Result/CancerOnly/DAA/Elix_score/Elix_score_ZicoSeq.Rdata')
eci_q <- diff.obj$qv.list$Species
eci_FC <- diff.obj$R2.list$Species * sign(diff.obj$coef.list$Species[,'Elix_score'])
sum(eci_q<cutoff)

tm <- load("Figure/CancerOnly/DAA_P_R2.RData")
sum(Q.All$Species$Elix_score<cutoff, na.rm = T)

## Check how many species also in ECI high but not in ECI low
old_sig <- names(eci_q[eci_q<cutoff,])
length(old_sig)

sum(old_sig %in% names(high_q[high_q<cutoff,])) # 274/318 species still show significance in high ECI cancer compared to HV
sum(old_sig %in% names(low_q[low_q<cutoff,])) # Only 216/318 species still show significance in low ECI cancer compared to HV

comp_q <- merge(as.data.frame(low_q) %>% dplyr::rename(low_q=Qvalue), as.data.frame(high_q) %>% dplyr::rename(high_q=Qvalue), by = 0, all=T) %>% column_to_rownames('Row.names') %>% 
  merge(as.data.frame(eci_q) %>% dplyr::rename(old=Qvalue), by = 0, all = T) %>% column_to_rownames('Row.names')

comp_fc <- merge(as.data.frame(low_FC) %>% dplyr::rename(low_FC=Func1), as.data.frame(high_FC) %>% dplyr::rename(high_FC=Func1), by = 0, all=T) %>% column_to_rownames('Row.names') %>% 
  merge(as.data.frame(eci_FC) %>% dplyr::rename(old_FC=Func1), by = 0, all=T) %>% column_to_rownames('Row.names')
comp_fc <- comp_fc[rownames(comp_q),]

old_high <- comp_fc[old_sig[old_sig %in% names(high_q[high_q<cutoff,])],]
sum(old_high$high_FC * old_high$old_FC>0)  #265/274/318 species still show significance in high ECI cancer compared to HV
old_high[which(old_high$high_FC * old_high$old_FC<0),]

old_low <- comp_fc[old_sig[old_sig %in% names(low_q[low_q<cutoff,])],]
sum(old_low$low_FC * old_low$old_FC>0)  #210/216/318 species still show significance in low ECI cancer compared to HV


q_fc <- merge(comp_q, comp_fc, by = 0, all = T) %>% column_to_rownames("Row.names")
colnames(q_fc) <- gsub('low_q','Low ECI vs HV (qvalue)',colnames(q_fc))
colnames(q_fc) <- gsub('high_q','High ECI vs HV (qvalue)',colnames(q_fc))
colnames(q_fc) <- gsub('old_q','ECI within cancer cohort(qvalue)',colnames(q_fc))
colnames(q_fc) <- gsub('low_FC','Low ECI vs HV (Effect size)',colnames(q_fc))
colnames(q_fc) <- gsub('high_FC','High ECI vs HV (Effect size)',colnames(q_fc))
colnames(q_fc) <- gsub('old_FC','ECI within cancer cohort(Effect size)',colnames(q_fc))


## How about exclude the pancancer significant ones?
load('Result/PanCancer/DAA/Group/Group_ZicoSeq.Rdata')
pan_sig <- names(diff.obj$qv.list$Species[diff.obj$qv.list$Species <cutoff,])
high_sig <- names(high_q[high_q<cutoff,])
low_sig <- names(low_q[low_q<cutoff,])
high_exc_pan_sig <- high_sig[!(high_sig %in% pan_sig)]
low_exc_pan_sig <- low_sig[!(low_sig %in% pan_sig)]

high_only_exc_pan_sig <- high_exc_pan_sig[!(high_exc_pan_sig %in% low_exc_pan_sig)]
low_only_exc_pan_sig <- low_exc_pan_sig[!(low_exc_pan_sig %in% high_exc_pan_sig)]
intersect(high_exc_pan_sig,low_exc_pan_sig)

b <- high_only_exc_pan_sig[high_only_exc_pan_sig %in% old_sig]
b
low_only_exc_pan_sig[low_only_exc_pan_sig %in% old_sig]








cutoff = 0.1
## Best results from ECI==0 (lowCancer VS healthy ECI<4 )
file_dir <- dirname(rstudioapi::getSourceEditorContext()$path)
file_dir <- dirname(dirname(dirname(dirname(file_dir))))
setwd(file_dir)

load('Data/data.obj.raw.core.RData')

tm = load('Result/ECI_subset0/Control-LowCancer/DAA/ECI22/ECI22_ZicoSeq.Rdata')
low_q <- diff.obj$qv.list$Species
low_FC <- diff.obj$R2.list$Species * sign(diff.obj$coef.list$Species[,'ECI22LowCancer'])
sum(low_q<cutoff)
df <- merge( diff.obj$pv.list$Species, low_q, by = 0)%>% column_to_rownames('Row.names') %>% merge(low_FC, by = 0) %>% dplyr::rename(Species = 'Row.names', `Effect size`=Func1)
df <- df %>% inner_join(as.data.frame(data.obj$otu.name))
#write.csv(df, file = "Code/Submission/Supplementary tables/Supplementary Table X_lowECI.csv", row.names = F)
df2 <- df[df$Qvalue<0.1,]
length(grep('Streptococcus',df2$Species))
df2$`Effect size`[grep('Streptococcus',df2$Species)]
length(grep('Bifidobacterium ',df2$Species))
df2$`Effect size`[grep('Bifidobacterium',df2$Species)]
length(grep('Blautia',df2$Species))
df2$`Effect size`[grep('Blautia',df2$Species)]
length(grep('Bacteroides',df2$Species))
df2$`Effect size`[grep('Bacteroides',df2$Species)]
length(grep('Faecalibacterium',df2$Species))
df2$`Effect size`[grep('Faecalibacterium',df2$Species)]


tm = load('Result/ECI_subset0/Control-HighCancer/DAA/ECI22/ECI22_ZicoSeq.Rdata')
high_q <- diff.obj$qv.list$Species
high_FC <- diff.obj$R2.list$Species * sign(diff.obj$coef.list$Species[,'ECI22HighCancer'])
sum(high_q<cutoff)

tm <- load('Result/CancerOnly/DAA/Elix_score/Elix_score_ZicoSeq.Rdata')
eci_q <- diff.obj$qv.list$Species
eci_FC <- diff.obj$R2.list$Species * sign(diff.obj$coef.list$Species[,'Elix_score'])
sum(eci_q<cutoff)

tm <- load("Figure/CancerOnly/DAA_P_R2.RData")
sum(Q.All$Species$Elix_score<cutoff, na.rm = T)

## Check how many species also in ECI high but not in ECI low
old_sig <- names(eci_q[eci_q<cutoff,])
length(old_sig)

sum(old_sig %in% names(high_q[high_q<cutoff,])) # 274/318 species still show significance in high ECI cancer compared to HV
sum(old_sig %in% names(low_q[low_q<cutoff,])) # Only 101/318 species still show significance in low ECI cancer compared to HV

comp_q <- merge(as.data.frame(low_q) %>% dplyr::rename(low_q=Qvalue), as.data.frame(high_q) %>% dplyr::rename(high_q=Qvalue), by = 0, all=T) %>% column_to_rownames('Row.names') %>% 
  merge(as.data.frame(eci_q) %>% dplyr::rename(old=Qvalue), by = 0, all = T) %>% column_to_rownames('Row.names')

comp_fc <- merge(as.data.frame(low_FC) %>% dplyr::rename(low_FC=Func1), as.data.frame(high_FC) %>% dplyr::rename(high_FC=Func1), by = 0, all=T) %>% column_to_rownames('Row.names') %>% 
  merge(as.data.frame(eci_FC) %>% dplyr::rename(old_FC=Func1), by = 0, all=T) %>% column_to_rownames('Row.names')
comp_fc <- comp_fc[rownames(comp_q),]

old_high <- comp_fc[old_sig[old_sig %in% names(high_q[high_q<cutoff,])],]
sum(old_high$high_FC * old_high$old_FC>0)  #265/274/318 species still show significance in high ECI cancer compared to HV
old_high[which(old_high$high_FC * old_high$old_FC<0),]

old_low <- comp_fc[old_sig[old_sig %in% names(low_q[low_q<cutoff,])],]
sum(old_low$low_FC * old_low$old_FC>0)  #94/101/318 species still show significance in low ECI cancer compared to HV


# q_fc <- merge(comp_q, comp_fc, by = 0, all = T) %>% column_to_rownames("Row.names")
# colnames(q_fc) <- gsub('low_q','Low ECI vs HV (qvalue)',colnames(q_fc))
# colnames(q_fc) <- gsub('high_q','High ECI vs HV (qvalue)',colnames(q_fc))
# colnames(q_fc) <- gsub('old_q','ECI within cancer cohort(qvalue)',colnames(q_fc))
# colnames(q_fc) <- gsub('low_FC','Low ECI vs HV (Effect size)',colnames(q_fc))
# colnames(q_fc) <- gsub('high_FC','High ECI vs HV (Effect size)',colnames(q_fc))
# colnames(q_fc) <- gsub('old_FC','ECI within cancer cohort(Effect size)',colnames(q_fc))
# write.csv(q_fc, file = "Code/Submission/Supplementary tables/Supplementary Table X_ECI.csv")



## How about exclude the pancancer significant ones?
load('Result/PanCancer/DAA/Group/Group_ZicoSeq.Rdata')
pan_sig <- names(diff.obj$qv.list$Species[diff.obj$qv.list$Species <cutoff,])
high_sig <- names(high_q[high_q<cutoff,])
low_sig <- names(low_q[low_q<cutoff,])
high_exc_pan_sig <- high_sig[!(high_sig %in% pan_sig)]
low_exc_pan_sig <- low_sig[!(low_sig %in% pan_sig)]

high_only_exc_pan_sig <- high_exc_pan_sig[!(high_exc_pan_sig %in% low_exc_pan_sig)]
low_only_exc_pan_sig <- low_exc_pan_sig[!(low_exc_pan_sig %in% high_exc_pan_sig)]
intersect(high_exc_pan_sig,low_exc_pan_sig)

b <- high_only_exc_pan_sig[high_only_exc_pan_sig %in% old_sig]
b
low_only_exc_pan_sig[low_only_exc_pan_sig %in% old_sig]


load("Data/data.obj.raw.core.RData")
alpha.obj2 <- generate_alpha_diversity(data.obj.rff, measures = "Shannon", rarefy = F)
data.obj.rff$meta.dat$ECI22 <- case_when(
  data.obj.rff$meta.dat$Elix_score ==0 ~ "LowCancer",
  data.obj.rff$meta.dat$Elix_score %in% c(1,2,3) ~ "MediumCancer",
  data.obj.rff$meta.dat$Elix_score >= 4 ~ "HighCancer"
)
data.obj.rff$meta.dat$ECI22[data.obj.rff$meta.dat$Group=='Control' & data.obj.rff$meta.dat$Elix_score<4] <- 'Control'
data.obj.rff$meta.dat$ECI22[data.obj.rff$meta.dat$Group=='Control' & data.obj.rff$meta.dat$Elix_score>=4] <- NA

merge(alpha.obj2, data.obj.rff$meta.dat[,'ECI22',drop =F], by = 0) %>% na.omit() %>% dplyr::filter(ECI22 !='MediumCancer') %>%
  mutate(ECI22 = gsub('LowCancer','Cancer(ECI=0)',ECI22),ECI22 = gsub('HighCancer','Cancer(ECI>3)',ECI22)) %>% 
  mutate(ECI22 = factor(ECI22,
                        levels = c("Control", "Cancer(ECI=0)",  "Cancer(ECI>3)"))) %>%
  ggplot(aes(x = ECI22, y = Shannon, fill= ECI22)) +
  geom_boxplot(outlier.size = 2) +
  scale_fill_brewer(palette = 'Set2') +
  labs(x = '') +
  theme_bw(base_size = 16) +
  theme(panel.grid = element_blank(),
        legend.position = "none",
        panel.border = element_rect(size = 1),
        axis.text = element_text(color = 'black'),
        axis.ticks = element_line(size = 1),
        axis.ticks.length = unit(0.2, "cm"))
ggsave(file = paste0(file_dir,'/Figure/RM_figures/AlphaDiv_HV_LowHighECI.pdf'), width = 5, height = 4)

data.obj.rff$meta.dat$ECI22 <- case_when(
  data.obj.rff$meta.dat$Elix_score <=1 ~ "LowCancer",
  data.obj.rff$meta.dat$Elix_score %in% c(2,3) ~ "MediumCancer",
  data.obj.rff$meta.dat$Elix_score >= 4 ~ "HighCancer"
)
data.obj.rff$meta.dat$ECI22[data.obj.rff$meta.dat$Group=='Control' & data.obj.rff$meta.dat$Elix_score<4] <- 'Control'
data.obj.rff$meta.dat$ECI22[data.obj.rff$meta.dat$Group=='Control' & data.obj.rff$meta.dat$Elix_score>=4] <- NA

merge(alpha.obj2, data.obj.rff$meta.dat[,'ECI22',drop =F], by = 0) %>% na.omit() %>% dplyr::filter(ECI22 !='MediumCancer') %>%
  mutate(ECI22 = gsub('LowCancer','Cancer(ECI=0,1)',ECI22),ECI22 = gsub('HighCancer','Cancer(ECI>3)',ECI22)) %>% 
  mutate(ECI22 = factor(ECI22,
                        levels = c("Control", "Cancer(ECI=0,1)",  "Cancer(ECI>3)"))) %>%
  ggplot(aes(x = ECI22, y = Shannon, fill= ECI22)) +
  geom_boxplot(outlier.size = 2) +
  scale_fill_brewer(palette = 'Set2') +
  labs(x = '') +
  theme_bw(base_size = 16) +
  theme(panel.grid = element_blank(),
        legend.position = "none",
        panel.border = element_rect(size = 1),
        axis.text = element_text(color = 'black'),
        axis.ticks = element_line(size = 1),
        axis.ticks.length = unit(0.2, "cm"))
ggsave(file = paste0(file_dir,'/Figure/RM_figures/AlphaDiv_HV_LowHighECI_1.pdf'), width = 5, height = 4)



## ====== Response to reviewer2's subset analysis: Check HV vs ABX analysis results======
file_dir <- dirname(rstudioapi::getSourceEditorContext()$path)
file_dir <- dirname(dirname(dirname(dirname(file_dir))))
setwd(file_dir)
load('Data/data.obj.raw.core.RData')


tm = load('Result/Abx_day_365_subset/Control-Y/DAA/Abx_day_365/Abx_day_365_ZicoSeq.Rdata')
Y_q <- as.data.frame(diff.obj$qv.list$Species)
Y_FC <- diff.obj$R2.list$Species * sign(diff.obj$coef.list$Species[,'Abx_day_365Y'])
sum(Y_q<cutoff)

tm = load('Result/Abx_day_365_subset/Control-N/DAA/Abx_day_365/Abx_day_365_ZicoSeq.Rdata')
N_q <- as.data.frame(diff.obj$qv.list$Species)
N_FC <- diff.obj$R2.list$Species * sign(diff.obj$coef.list$Species[,'Abx_day_365N'])
sum(N_q<cutoff)

tm = load('Result/CancerOnly/DAA/Abx_day_365/Abx_day_365_ZicoSeq.Rdata')
Abx_q <- as.data.frame(diff.obj$qv.list$Species)
Abx_FC <- diff.obj$R2.list$Species * sign(diff.obj$coef.list$Species[,'Abx_day_365Y'])
sum(Abx_q<cutoff)

comm_old_Y <- intersect(rownames(Abx_q)[Abx_q$Qvalue<cutoff], rownames(Y_q)[Y_q$Qvalue<cutoff])
length(comm_old_Y)
comm_old_N <- intersect(rownames(Abx_q)[Abx_q$Qvalue<cutoff], rownames(N_q)[N_q$Qvalue<cutoff])
length(comm_old_N)
all_FC <- merge(as.data.frame(Y_FC) %>% dplyr::rename(YvsHV=Func1), 
                as.data.frame(Abx_FC) %>% dplyr::rename(Abx=Func1), by= 0, all=T) %>% column_to_rownames("Row.names") %>% 
  merge(as.data.frame(N_FC) %>% dplyr::rename(NvsHV=Func1), by= 0, all=T) %>% column_to_rownames("Row.names")

all_q <- merge(as.data.frame(Y_q) %>% dplyr::rename(YvsHV=Qvalue), 
                as.data.frame(Abx_q) %>% dplyr::rename(Abx=Qvalue), by= 0, all=T) %>% column_to_rownames("Row.names") %>% 
  merge(as.data.frame(N_q) %>% dplyr::rename(NvsHV=Qvalue), by= 0, all=T) %>% column_to_rownames("Row.names")
all_q <- all_q[rownames(all_FC),]

## How many Y vs HV sig & within cancer Abx sig and same direction?
sigY <- rownames(all_q)[all_q$YvsHV < cutoff & !is.na(all_q$YvsHV)]
sigN <- rownames(all_q)[all_q$NvsHV < cutoff & !is.na(all_q$NvsHV)]
sigAbx <- rownames(all_q)[all_q$Abx < cutoff & !is.na(all_q$Abx)]
Y_Abx <- rownames(all_FC[all_FC$YvsHV * all_FC$Abx > 0 & !is.na(all_FC$YvsHV) & !is.na(all_FC$Abx),])
intersect(Y_Abx,intersect(sigY, sigAbx)) # 189/333 remained significance with same direction

N_Abx <- rownames(all_FC[all_FC$NvsHV * all_FC$Abx > 0 & !is.na(all_FC$NvsHV) & !is.na(all_FC$Abx),])
intersect(N_Abx,intersect(sigY, sigAbx)) # 178/333 remained significance with same direction



table(data.obj$meta.dat$PPI_day_365)
tm = load('Result/PPI_day_365_subset/Control-Y/DAA/PPI_day_365/PPI_day_365_ZicoSeq.Rdata')
Y_q <- as.data.frame(diff.obj$qv.list$Species)
Y_FC <- diff.obj$R2.list$Species * sign(diff.obj$coef.list$Species[,'PPI_day_365Y'])
sum(Y_q<cutoff)

tm = load('Result/PPI_day_365_subset/Control-N/DAA/PPI_day_365/PPI_day_365_ZicoSeq.Rdata')
N_q <- as.data.frame(diff.obj$qv.list$Species)
N_FC <- diff.obj$R2.list$Species * sign(diff.obj$coef.list$Species[,'PPI_day_365N'])
sum(N_q<cutoff)

tm = load('Result/CancerOnly/DAA/PPI_day_365/PPI_day_365_ZicoSeq.Rdata')
PPI_q <- as.data.frame(diff.obj$qv.list$Species)
PPI_FC <- diff.obj$R2.list$Species * sign(diff.obj$coef.list$Species[,'PPI_day_365Y'])
sum(PPI_q<cutoff)

comm_old_Y <- intersect(rownames(PPI_q)[PPI_q$Qvalue<cutoff], rownames(Y_q)[Y_q$Qvalue<cutoff])
length(comm_old_Y)
comm_old_N <- intersect(rownames(PPI_q)[PPI_q$Qvalue<cutoff], rownames(N_q)[N_q$Qvalue<cutoff])
length(comm_old_N)
all_FC <- merge(as.data.frame(Y_FC) %>% dplyr::rename(YvsHV=Func1), 
                as.data.frame(PPI_FC) %>% dplyr::rename(PPI=Func1), by= 0, all=T) %>% column_to_rownames("Row.names") %>% 
  merge(as.data.frame(N_FC) %>% dplyr::rename(NvsHV=Func1), by= 0, all=T) %>% column_to_rownames("Row.names")

all_q <- merge(as.data.frame(Y_q) %>% dplyr::rename(YvsHV=Qvalue), 
               as.data.frame(PPI_q) %>% dplyr::rename(PPI=Qvalue), by= 0, all=T) %>% column_to_rownames("Row.names") %>% 
  merge(as.data.frame(N_q) %>% dplyr::rename(NvsHV=Qvalue), by= 0, all=T) %>% column_to_rownames("Row.names")
all_q <- all_q[rownames(all_FC),]

## How many Y vs HV sig & within cancer PPI sig and same direction?
sigY <- rownames(all_q)[all_q$YvsHV < cutoff & !is.na(all_q$YvsHV)]
sigN <- rownames(all_q)[all_q$NvsHV < cutoff & !is.na(all_q$NvsHV)]
sigPPI <- rownames(all_q)[all_q$PPI < cutoff & !is.na(all_q$PPI)]
Y_PPI <- rownames(all_FC[all_FC$YvsHV * all_FC$PPI > 0 & !is.na(all_FC$YvsHV) & !is.na(all_FC$PPI),]) # PPI & PPI_Y vs HV
intersect(Y_PPI,intersect(sigY, sigPPI)) # 8/14 remained significance with same direction
N_PPI <- rownames(all_FC[all_FC$NvsHV * all_FC$PPI > 0 & !is.na(all_FC$NvsHV) & !is.na(all_FC$PPI),]) # PPI & PPI_N vs HV
intersect(N_PPI,intersect(sigN, sigPPI)) # 8/14 remained significance with same direction

sum(sigPPI %in% sigY)




load('Result/CancerOnly_abx/DAA/Abx_day_365/Abx_day_365_ZicoSeq.Rdata')
sum(diff.obj$qv.list$Species<0.1)
marginal <- cbind(diff.obj$qv.list$Species,diff.obj$R2.list$Species * sign(diff.obj$coef.list[[1]])) %>% as.data.frame
colnames(marginal)[2] <- 'EffectSize'
write.csv(marginal, file = paste0(file_dir,'/Figure/RM_Figures/marginal_abx_day_365.csv'))


cancers <- list.dirs('Result/ECI_subset/', recursive = F, full.names = F)
cancers <- cancers[-(grep('Control',cancers))]
sigs <- c()
ress <- NULL
for(cancer in cancers){
  diff.obj <- NULL
  load(paste0('Result/ECI_subset/',cancer,'/Control-LowCancer/DAA/ECI2/ECI2_ZicoSeq.Rdata'))
  res <- cbind(diff.obj$qv.list$Species,diff.obj$R2.list$Species * sign(diff.obj$coef.list[[1]][,'ECI2LowCancer'])) %>% cbind(diff.obj$pv.list$Species) %>% as.data.frame
  colnames(res)[2] <- 'EffectSize'
  sigs <- c(sigs, sum(res$Qvalue<0.1))
  cat(cancer,':')
  cat(sum(res$Qvalue<0.1),'\n')
  ress <- rbind(ress, res %>% mutate(cancer = cancer) %>% rownames_to_column('species'))
}
names(sigs) <- cancers
write.csv(ress, file ="Code/Submission/Supplementary tables/Supplementary Table X_HV_lowECI.csv", row.names = F)


cancers <- list.dirs('Result/subCancerX_Control/', recursive = F, full.names = F)
cancers <- gsub('Control-','',cancers)
sigs_ctrl <- c()
ress2 <- NULL
for(cancer in cancers){
  diff.obj <- NULL
  load(paste0('Result/subCancerX_Control/Control-',cancer,'/DAA/icd10_first_3_name_short/icd10_first_3_name_short_ZicoSeq.Rdata'))
  idx <- colnames(diff.obj$coef.list[[1]])[grep('icd10_first_3_name_short',colnames(diff.obj$coef.list[[1]]))]
  res <- cbind(diff.obj$qv.list$Species,diff.obj$R2.list$Species * sign(diff.obj$coef.list[[1]][,idx])) %>% cbind(diff.obj$pv.list$Species) %>% as.data.frame
  colnames(res)[2] <- 'EffectSize'
  sigs_ctrl <- c(sigs_ctrl, sum(res$Qvalue<0.1))
  cat(cancer,':')
  cat(sum(res$Qvalue<0.1),'\n')
  ress2 <- rbind(ress2, res %>% mutate(cancer = cancer) %>% rownames_to_column('species'))
}
names(sigs_ctrl) <- cancers





file_dir <- dirname(rstudioapi::getSourceEditorContext()$path)
file_dir <- dirname(dirname(dirname(dirname(file_dir))))
setwd(file_dir)
load(file = 'Data/data.obj.raw.core.RData') 
data.obj$meta.dat$ECI2 <- case_when(
  data.obj$meta.dat$Elix_score <=1  ~ "LowCancer",
  data.obj$meta.dat$Elix_score %in% c(2,3) ~ "MediumCancer",
  data.obj$meta.dat$Elix_score >= 4 ~ "HighCancer"
)
data.obj$meta.dat$ECI2[data.obj$meta.dat$Group=='Control'& data.obj$meta.dat$Elix_score<4] <- 'Control'
data.obj$meta.dat$ECI2[data.obj$meta.dat$Group=='Control' & data.obj$meta.dat$Elix_score>=4] <- NA
cancers <- cbind(table(data.obj$meta.dat$icd10_first_3_name_short,data.obj$meta.dat$ECI2))[,"LowCancer"]


comp <- cbind.data.frame(lowECIvsControl = sigs, subCancerXvsControl = sigs_ctrl[names(sigs)]) %>% 
  mutate(diff = subCancerXvsControl-lowECIvsControl, `ECI<=1(n)` = cancers[names(sigs)])



### ======Add smoking volcano plot ======
load('/Users/luyang1/myicloud/Documents/Mayo_project/2023_01_09_Oncobiome/mforge_clean/Result/CancerOnly/DAA/smoking_category/smoking_category_ZicoSeq.Rdata')
sum(diff.obj$qv.list$Species<0.1)
dim(diff.obj$qv.list$Species)
p1 <- merge(diff.obj$qv.list$Species, diff.obj$R2.list$Species * sign(diff.obj$coef.list$Species[,'smoking_categoryFormer']), by = 0, all.x = T) %>% 
  mutate(label = ifelse(Qvalue<0.05 & abs(Func1)>0.01, gsub('s__','',`Row.names`),'')) %>% 
  
  ggplot(aes(x = Func1, y = -log10(Qvalue)))+ 
  geom_point(aes(color = factor(ifelse(Qvalue < 0.1, "sig", "no")))) + 
  geom_hline(yintercept = -log10(0.1), linetype = "dashed") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_text_repel(aes(label = label),
                  size = 3,
                  fontface = "italic",
                  max.overlaps = 50,
                  box.padding = 0.3,
                  point.padding = 0.2,
                  na.rm = TRUE) +
  annotate("text", x = 0.01, y = 3, label = "Former smoker", 
           hjust = 1.1, vjust = 1.5, size = 5) +
  annotate("text", x = -0.01, y = 3, label = "Current smoker", 
           hjust = -0.1, vjust = 1.5, size = 5)+
  scale_color_manual(
    values = c("sig" = "red", "no" = "grey"),
    name = "qvalue<0.1"
  ) +
  theme_classic() + 
  labs(x = 'Effect size (R² x direction)', y = '-log10(qvalue)')
p2 <- merge(diff.obj$qv.list$Species, diff.obj$R2.list$Species * sign(diff.obj$coef.list$Species[,'smoking_categoryNever']), by = 0, all.x = T) %>% 
  mutate(label = ifelse(Qvalue<0.05 & abs(Func1)>0.01, gsub('s__','',`Row.names`),'')) %>% 
  
  ggplot(aes(x = Func1, y = -log10(Qvalue)))+ 
  geom_point(aes(color = factor(ifelse(Qvalue < 0.1, "sig", "no")))) + 
  geom_hline(yintercept = -log10(0.1), linetype = "dashed") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_text_repel(aes(label = label),
                  size = 3,
                  fontface = "italic",
                  max.overlaps = 50,
                  box.padding = 0.3,
                  point.padding = 0.2,
                  na.rm = TRUE) +
  annotate("text", x = 0.01, y = 3, label = "Never smoke", 
           hjust = 1.1, vjust = 1.5, size = 5) +
  annotate("text", x = -0.01, y = 3, label = "Current smoker", 
           hjust = -0.1, vjust = 1.5, size = 5)+
  scale_color_manual(
    values = c("sig" = "red", "no" = "grey"),
    name = "qvalue<0.1"
  ) +
  theme_classic() + 
  labs(x = 'Effect size (R² x direction)', y = '-log10(qvalue)')
ggarrange(p1, p2, common.legend = T, nrow = 2)
ggsave(file = paste0(file_dir,'/Figure/RM_figures/Smoking_Species_volcano_pairwise.pdf'), width = 12, height = 15)

merge(diff.obj$qv.list$Species, diff.obj$R2.list$Species, by = 0, all.x = T) %>% 
  column_to_rownames('Row.names') %>% 
  merge(diff.obj$coef.list$Species[,'smoking_categoryFormer'], by = 0, all.x = T) %>% 
  column_to_rownames('Row.names') %>% dplyr::rename(smoking_categoryFormer=y) %>% 
  merge(diff.obj$coef.list$Species[,'smoking_categoryNever'], by = 0, all.x = T) %>% 
  dplyr::rename(smoking_categoryNever=y) %>% 
  mutate(label = ifelse(Qvalue<0.05 & Func1>0.01, gsub('s__','',`Row.names`),'')) %>% 
  ggplot(aes(x = Func1, y = -log10(Qvalue)))+ 
  geom_point(aes(color = factor(ifelse(Qvalue < 0.1, "sig", "no")))) + 
  geom_hline(yintercept = -log10(0.1), linetype = "dashed") +
  geom_text_repel(aes(label = label),
                  size = 3,
                  fontface = "italic",
                  max.overlaps = 50,
                  box.padding = 0.3,
                  point.padding = 0.2,
                  na.rm = TRUE) +
  scale_color_manual(
    values = c("sig" = "red", "no" = "grey"),
    name = "qvalue<0.1"
  ) +
  theme_classic() + 
  labs(x = 'Effect size (R²)', y = '-log10(qvalue)')
ggsave(file = paste0(file_dir,'/Figure/RM_figures/Smoking_Species_volcano.pdf'), width = 12, height = 10)






##### ===== Abx/PPI + cancer class ===######
file_dir <- dirname(rstudioapi::getSourceEditorContext()$path)
file_dir <- dirname(dirname(dirname(dirname(file_dir))))
setwd(file_dir)

ppi_null <- read.csv('Result/CancerOnly_abx/DAA/PPI_day_365/Taxa_ZicoSeq_Species_.PPI_day_365.adj..csv')
sum(ppi_null$Qvalue<0.1)
sum(ppi_null$Qvalue<0.05)
ppi_adjCancer <- read.csv('Result/CancerOnly_abx/DAA/PPI_day_365/Taxa_ZicoSeq_Species_.PPI_day_365.adj.Cancer_class.csv')
sum(ppi_adjCancer$Qvalue<0.1)
sum(ppi_adjCancer$Qvalue<0.05)

abx_null <- read.csv('Result/CancerOnly_abx/DAA/Abx_day_365/Taxa_ZicoSeq_Species_.Abx_day_365.adj..csv')
sum(abx_null$Qvalue<0.1)
sum(abx_null$Qvalue<0.05)
abx_adjCancer <- read.csv('Result/CancerOnly_abx/DAA/Abx_day_365/Taxa_ZicoSeq_Species_.Abx_day_365.adj.Cancer_class.csv')
sum(abx_adjCancer$Qvalue<0.1)
sum(abx_adjCancer$Qvalue<0.05)

abx_adjAll <- read.csv('Result/CancerOnly/DAA/Abx_day_365/Taxa_ZicoSeq_Species_.Abx_day_365.adj.Batch.Bristol_score.BMI.Age.Sex.Cancer_class.Metastasis.PPI_day_365.Elix_score.Sample_season.Urban.csv')
sum(abx_adjAll$Qvalue<0.1)
sum(abx_adjAll$Qvalue<0.05)


### === low ECI vs high ECI vs HV distance =====
file_dir <- dirname(rstudioapi::getSourceEditorContext()$path)
file_dir <- dirname(dirname(dirname(dirname(file_dir))))

wd <- file_dir
rd <- paste0(file_dir,'/Result/')

if(!dir.exists(rd)) dir.create(rd)
setwd(wd)
source(paste0(file_dir,"/Code/Submission/MayoOncobiomeStudy/Code/Stats.R"))
try(load_package())

load(file = 'Data/data.obj.raw.core.RData') 

setwd(rd)
dir <- 'ECI_subset'
if(!dir.exists(dir)){dir.create(dir)}
setwd(dir)
getwd()

data.obj$meta.dat$ECI2 <- case_when(
  data.obj$meta.dat$Elix_score <=1  ~ "LowCancer",
  data.obj$meta.dat$Elix_score %in% c(2,3) ~ "MediumCancer",
  data.obj$meta.dat$Elix_score >= 4 ~ "HighCancer"
)
data.obj$meta.dat$ECI2[data.obj$meta.dat$Group=='Control'& data.obj$meta.dat$Elix_score<4] <- 'Control'
data.obj$meta.dat$ECI2[data.obj$meta.dat$Group=='Control' & data.obj$meta.dat$Elix_score>=4] <- NA

data.obj.rff$meta.dat$ECI2 <- case_when(
  data.obj.rff$meta.dat$Elix_score <=1 ~ "LowCancer",
  data.obj.rff$meta.dat$Elix_score %in% c(2,3) ~ "MediumCancer",
  data.obj.rff$meta.dat$Elix_score >= 4 ~ "HighCancer"
)
data.obj.rff$meta.dat$ECI2[data.obj.rff$meta.dat$Group=='Control' & data.obj.rff$meta.dat$Elix_score<4] <- 'Control'
data.obj.rff$meta.dat$ECI2[data.obj.rff$meta.dat$Group=='Control' & data.obj.rff$meta.dat$Elix_score>=4] <- NA

# Low ECI > 5 samples & cancerX contains more than 15 samples
cancers <- names(which(sort(cbind(table(data.obj$meta.dat$icd10_first_3_name_short,data.obj$meta.dat$ECI2))[,"LowCancer"]) > 5))
cancers2 <- names(which(table(data.obj$meta.dat$icd10_first_3_name_short)>15))
cancers2 <- cancers2[cancers2!="healthy"]
cancers <- intersect(cancers, cancers2)

ind <- data.obj$meta.dat$icd10_first_3_name_short %in% c('healthy',cancers)
data.obj2 <- subset_data(data.obj, ind)
dist.obj2 <- subset_dist(dist.obj, ind)

ind <- data.obj.rff$meta.dat$icd10_first_3_name_short %in% c('healthy',cancers)
data.obj.rff2 <- subset_data(data.obj.rff, ind)
dist.obj.rff2 <- subset_dist(dist.obj.rff, ind)



idx_hv <- rownames(data.obj.rff$meta.dat)[data.obj.rff$meta.dat$Group == "Control"& data.obj.rff$meta.dat$Elix_score <2 & !is.na(data.obj.rff$meta.dat$Elix_score)]
dist_hv <- dist.obj.rff$WUniFrac[idx_hv,idx_hv]

idx_cancer <- rownames(data.obj.rff$meta.dat)[data.obj.rff$meta.dat$Group != "Control"]
dist_cancer <- dist.obj.rff$WUniFrac[idx_cancer,idx_cancer]


dist.list.low <- dist.list.high <- list()
nSam_low <- nSam_high <-c()
for(cancer in cancers){
  full_matrix <- as.matrix(dist.obj.rff2$WUniFrac)
  
  sub_low <- rownames(data.obj.rff2$meta.dat)[data.obj.rff2$meta.dat$icd10_first_3_name_short == cancer & data.obj.rff2$meta.dat$ECI2 == "LowCancer"]
  dist_low_sub <- dist.obj.rff2$WUniFrac[sub_low,sub_low]
  cross_dist_low <- full_matrix[rownames(dist_low_sub), rownames(dist_hv)]
  
  sub_high <- rownames(data.obj.rff2$meta.dat)[data.obj.rff2$meta.dat$icd10_first_3_name_short == cancer & data.obj.rff2$meta.dat$ECI2 == "HighCancer"]
  dist_high_sub <- dist.obj.rff2$WUniFrac[sub_high,sub_high]
  cross_dist_high <- full_matrix[rownames(dist_high_sub), rownames(dist_hv)]
  
  dist.list.low[[cancer]] <- cross_dist_low
  dist.list.high[[cancer]] <- cross_dist_high
  
  nSam_low <- c(nSam_low, length(sub_low))
  nSam_high <- c(nSam_high, length(sub_high))
}
names(nSam_low) <- names(nSam_high) <- cancers

df <- NULL
for(cancer in cancers){
  low_tmp <- dist.list.low[[cancer]]
  low_tmp[!lower.tri(low_tmp, diag = FALSE)] <- NA
  high_tmp <- dist.list.high[[cancer]]
  high_tmp[!lower.tri(high_tmp, diag = FALSE)] <- NA
  tmp1 <- melt(low_tmp) %>% mutate(cancer = cancer, ECI = 'Low ECI vs Healthy') 
  tmp2 <- melt(high_tmp) %>% mutate(cancer = cancer, ECI = 'High ECI vs Healthy') 
  df <- rbind(df, rbind(tmp1, tmp2))
}

df <- na.omit(df)
head(df)
median_diff <- aggregate(value~cancer + ECI, data = df, function(x) median(x)) %>% 
  pivot_wider(
    names_from = ECI,
    values_from = value
  ) %>% mutate(diff = `High ECI vs Healthy`-`Low ECI vs Healthy`) %>% arrange(diff)
opp <- median_diff$cancer[median_diff$diff<0]
sort(cbind(table(data.obj$meta.dat$icd10_first_3_name_short,data.obj$meta.dat$ECI2))[,"LowCancer"])[opp]
df$cancer <- factor(df$cancer, levels = median_diff$cancer)

# prepare table for label
df_n <- melt(cbind(nSam_low, nSam_high)) %>% arrange(X1) %>% mutate(y = rep(0.8, length(c(nSam_low, nSam_high))))
colnames(df_n) <- c('cancer','ECI','n','y')
df_n$ECI <- as.character(df_n$ECI)
df_n$ECI[df_n$ECI=='nSam_high'] <- 'High ECI vs Healthy'
df_n$ECI[df_n$ECI=='nSam_low'] <- 'Low ECI vs Healthy'

ggplot(df, aes(x = cancer, y = value, fill = ECI)) +
  geom_boxplot(position = position_dodge(width = 0.8), outlier.size = 0.2, outlier.colour = 'grey50') +
  geom_text(
    data = df_n,
    aes(x = cancer, y = y, label = paste0("n=", n), group = ECI),
    position = position_dodge(width = 0.8),
    color = 'black',
    size = 3,
    angle = 270,      
    vjust = 0.5,    
    hjust = 0.5
  )+
  scale_fill_manual(values = c('Low ECI vs Healthy'='#0072B2','High ECI vs Healthy'='#CC79A7',
                               'HV'='steelblue','Cancer'='orange')) +
  theme_bw() +
  labs(x = "", y = "Weighted UniFrac distance", fill = "") +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, color = 'black'),
    axis.text.y = element_text(color = 'black'), 
    legend.position = 'top'
  )  
ggsave(file = paste0(file_dir,'/Figure/RM_figures/ECI_LowHighHV_distance.pdf'), width = 10, height = 6)

# median_diff$cancer <- factor(median_diff$cancer, levels = median_diff$cancer)
# ggplot(median_diff, aes(x = cancer, y = diff, fill = diff > 0)) +
#   geom_bar(stat = 'identity') +
#   scale_fill_manual(values = c("TRUE" = '#CC79A7', "FALSE" = '#0072B2')) +
#   theme_bw() +
#   labs(x = "", y = "WUniFrac (median Cancer- median Healthy)", fill = "") +
#   theme(
#     axis.text.x = element_text(angle = 90, hjust = 1, color = 'black'),
#     axis.text.y = element_text(color = 'black'), 
#     legend.position = "none"
#   )
# ggsave(file = paste0(file_dir,'/Figure/RM_figures/ECI_LowHighHV_distance_diffBar.pdf'), width = 6, height =5)
# 




### === Abx Y /N vs HV distance =====
file_dir <- dirname(rstudioapi::getSourceEditorContext()$path)
file_dir <- dirname(dirname(dirname(dirname(file_dir))))

wd <- file_dir
rd <- paste0(file_dir,'/Result/')

if(!dir.exists(rd)) dir.create(rd)
setwd(wd)
source(paste0(file_dir,"/Code/Submission/MayoOncobiomeStudy/Code/Stats.R"))
try(load_package())

load(file = 'Data/data.obj.raw.core.RData') 

setwd(rd)
dir <- 'Abx_day_365_subset/'
if(!dir.exists(dir)){dir.create(dir)}
setwd(dir)
getwd()

# Abx_day_365 Y/N both > 5 samples & this cancer contains more than 15 samples
cancers <- names(which(rowSums(cbind(table(data.obj$meta.dat$icd10_first_3_name_short,data.obj$meta.dat$Abx_day_365))>5)==2))
cancers2 <- names(which(table(data.obj$meta.dat$icd10_first_3_name_short)>15))
cancers2 <- cancers2[cancers2!="healthy"]
cancers <- intersect(cancers, cancers2)


ind <- data.obj$meta.dat$icd10_first_3_name_short %in% c('healthy',cancers)
data.obj2 <- subset_data(data.obj, ind)
dist.obj2 <- subset_dist(dist.obj, ind)

ind <- data.obj.rff$meta.dat$icd10_first_3_name_short %in% c('healthy',cancers)
data.obj.rff2 <- subset_data(data.obj.rff, ind)
dist.obj.rff2 <- subset_dist(dist.obj.rff, ind)

idx_hv <- rownames(data.obj.rff$meta.dat)[data.obj.rff$meta.dat$Group == "Control"& !is.na(data.obj.rff$meta.dat$Elix_score)]
dist_hv <- dist.obj.rff$WUniFrac[idx_hv,idx_hv]

idx_cancer <- rownames(data.obj.rff$meta.dat)[data.obj.rff$meta.dat$Group != "Control"]
dist_cancer <- dist.obj.rff$WUniFrac[idx_cancer,idx_cancer]


dist.list.low <- dist.list.high <- list()
nSam_low <- nSam_high <-c()
for(cancer in cancers){
  full_matrix <- as.matrix(dist.obj.rff2$WUniFrac)
  
  sub_low <- rownames(data.obj.rff2$meta.dat)[data.obj.rff2$meta.dat$icd10_first_3_name_short == cancer & data.obj.rff2$meta.dat$Abx_day_365 == "N"]
  dist_low_sub <- dist.obj.rff2$WUniFrac[sub_low,sub_low]
  cross_dist_low <- full_matrix[rownames(dist_low_sub), rownames(dist_hv)]
  
  sub_high <- rownames(data.obj.rff2$meta.dat)[data.obj.rff2$meta.dat$icd10_first_3_name_short == cancer & data.obj.rff2$meta.dat$Abx_day_365 == "Y"]
  dist_high_sub <- dist.obj.rff2$WUniFrac[sub_high,sub_high]
  cross_dist_high <- full_matrix[rownames(dist_high_sub), rownames(dist_hv)]
  
  dist.list.low[[cancer]] <- cross_dist_low
  dist.list.high[[cancer]] <- cross_dist_high
  
  nSam_low <- c(nSam_low, length(sub_low))
  nSam_high <- c(nSam_high, length(sub_high))
}
names(nSam_low) <- names(nSam_high) <- cancers

df <- NULL
for(cancer in cancers){
  low_tmp <- dist.list.low[[cancer]]
  low_tmp[!lower.tri(low_tmp, diag = FALSE)] <- NA
  high_tmp <- dist.list.high[[cancer]]
  high_tmp[!lower.tri(high_tmp, diag = FALSE)] <- NA
  tmp1 <- melt(low_tmp) %>% mutate(cancer = cancer, Antibiotics = 'Antibiotics within past year(No) vs Healthy') 
  tmp2 <- melt(high_tmp) %>% mutate(cancer = cancer, Antibiotics = 'Antibiotics within past year(Yes) vs Healthy') 
  df <- rbind(df, rbind(tmp1, tmp2))
}

df <- na.omit(df)
head(df)
median_diff <- aggregate(value~cancer + Antibiotics, data = df, function(x) median(x)) %>% 
  pivot_wider(
    names_from = Antibiotics,
    values_from = value
  ) %>% mutate(diff = `Antibiotics within past year(Yes) vs Healthy`-`Antibiotics within past year(No) vs Healthy`) %>% arrange(diff)
opp <- median_diff$cancer[median_diff$diff<0]
sort(cbind(table(data.obj$meta.dat$icd10_first_3_name_short,data.obj$meta.dat$Abx_day_365))[,"N"])[opp]
df$cancer <- factor(df$cancer, levels = median_diff$cancer)

# prepare table for label
df_n <- melt(cbind(nSam_low, nSam_high)) %>% arrange(X1) %>% mutate(y = rep(0.8, length(c(nSam_low, nSam_high))))
colnames(df_n) <- c('cancer','Antibiotics','n','y')
df_n$Antibiotics <- as.character(df_n$Antibiotics)
df_n$Antibiotics[df_n$Antibiotics=='nSam_high'] <- 'Antibiotics within past year(Yes) vs Healthy'
df_n$Antibiotics[df_n$Antibiotics=='nSam_low'] <- 'Antibiotics within past year(No) vs Healthy'

ggplot(df, aes(x = cancer, y = value, fill = Antibiotics)) +
  geom_boxplot(position = position_dodge(width = 0.8), outlier.size = 0.2, outlier.colour = 'grey50') +
  geom_text(
    data = df_n,
    aes(x = cancer, y = y, label = paste0("n=", n), group = Antibiotics),
    position = position_dodge(width = 0.8),
    color = 'black',
    size = 3,
    angle = 270,      
    vjust = 0.5,    
    hjust = 0.5
  )+
  scale_fill_manual(values = c('Antibiotics within past year(No) vs Healthy'='#0072B2',
                               'Antibiotics within past year(Yes) vs Healthy'='#CC79A7',
                               'HV'='steelblue','Cancer'='orange')) +
  theme_bw() +
  labs(x = "", y = "Weighted UniFrac distance", fill = "") +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, color = 'black'),
    axis.text.y = element_text(color = 'black'), 
    legend.position = 'top'
  )  
ggsave(file = paste0(file_dir,'/Figure/RM_figures/Abx_YN_HV_distance.pdf'), width = 10, height = 6)

# median_diff$cancer <- factor(median_diff$cancer, levels = median_diff$cancer)
# ggplot(median_diff, aes(x = cancer, y = diff, fill = diff > 0)) +
#   geom_bar(stat = 'identity') +
#   scale_fill_manual(values = c("TRUE" = '#CC79A7', "FALSE" = '#0072B2')) +
#   scale_y_continuous(limits = c(-0.08,0.08)) +
#   theme_bw() +
#   labs(x = "", y = "WUniFrac (median AbxY vs HV - median AbxN vs HV)", fill = "") +
#   theme(
#     axis.text.x = element_text(angle = 90, hjust = 1, color = 'black'),
#     axis.text.y = element_text(color = 'black'),
#     legend.position = "none"
#   )
# ggsave(file = paste0(file_dir,'/Figure/RM_figures/Abx_YN_HV_distance_diffBar.pdf'), width = 4, height =6)




### === PPI Y /N vs HV distance =====
file_dir <- dirname(rstudioapi::getSourceEditorContext()$path)
file_dir <- dirname(dirname(dirname(dirname(file_dir))))

wd <- file_dir
rd <- paste0(file_dir,'/Result/')

if(!dir.exists(rd)) dir.create(rd)
setwd(wd)
source(paste0(file_dir,"/Code/Submission/MayoOncobiomeStudy/Code/Stats.R"))
try(load_package())

load(file = 'Data/data.obj.raw.core.RData') 

setwd(rd)
dir <- 'PPI_day_365_subset/'
if(!dir.exists(dir)){dir.create(dir)}
setwd(dir)
getwd()

# PPI_day_365 Y/N both > 5 samples & this cancer contains more than 15 samples
cancers <- names(which(rowSums(cbind(table(data.obj$meta.dat$icd10_first_3_name_short,data.obj$meta.dat$PPI_day_365))>5)==2))
cancers2 <- names(which(table(data.obj$meta.dat$icd10_first_3_name_short)>15))
cancers2 <- cancers2[cancers2!="healthy"]
cancers <- intersect(cancers, cancers2)


ind <- data.obj$meta.dat$icd10_first_3_name_short %in% c('healthy',cancers)
data.obj2 <- subset_data(data.obj, ind)
dist.obj2 <- subset_dist(dist.obj, ind)

ind <- data.obj.rff$meta.dat$icd10_first_3_name_short %in% c('healthy',cancers)
data.obj.rff2 <- subset_data(data.obj.rff, ind)
dist.obj.rff2 <- subset_dist(dist.obj.rff, ind)

idx_hv <- rownames(data.obj.rff$meta.dat)[data.obj.rff$meta.dat$Group == "Control"& !is.na(data.obj.rff$meta.dat$Elix_score)]
dist_hv <- dist.obj.rff$WUniFrac[idx_hv,idx_hv]

idx_cancer <- rownames(data.obj.rff$meta.dat)[data.obj.rff$meta.dat$Group != "Control"]
dist_cancer <- dist.obj.rff$WUniFrac[idx_cancer,idx_cancer]


dist.list.low <- dist.list.high <- list()
nSam_low <- nSam_high <-c()
for(cancer in cancers){
  full_matrix <- as.matrix(dist.obj.rff2$WUniFrac)
  
  sub_low <- rownames(data.obj.rff2$meta.dat)[data.obj.rff2$meta.dat$icd10_first_3_name_short == cancer & data.obj.rff2$meta.dat$PPI_day_365 == "N"]
  dist_low_sub <- dist.obj.rff2$WUniFrac[sub_low,sub_low]
  cross_dist_low <- full_matrix[rownames(dist_low_sub), rownames(dist_hv)]
  
  sub_high <- rownames(data.obj.rff2$meta.dat)[data.obj.rff2$meta.dat$icd10_first_3_name_short == cancer & data.obj.rff2$meta.dat$PPI_day_365 == "Y"]
  dist_high_sub <- dist.obj.rff2$WUniFrac[sub_high,sub_high]
  cross_dist_high <- full_matrix[rownames(dist_high_sub), rownames(dist_hv)]
  
  dist.list.low[[cancer]] <- cross_dist_low
  dist.list.high[[cancer]] <- cross_dist_high
  
  nSam_low <- c(nSam_low, length(sub_low))
  nSam_high <- c(nSam_high, length(sub_high))
}
names(nSam_low) <- names(nSam_high) <- cancers

df <- NULL
for(cancer in cancers){
  low_tmp <- dist.list.low[[cancer]]
  low_tmp[!lower.tri(low_tmp, diag = FALSE)] <- NA
  high_tmp <- dist.list.high[[cancer]]
  high_tmp[!lower.tri(high_tmp, diag = FALSE)] <- NA
  tmp1 <- melt(low_tmp) %>% mutate(cancer = cancer, PPI = 'PPI within past year(No) vs Healthy') 
  tmp2 <- melt(high_tmp) %>% mutate(cancer = cancer, PPI = 'PPI within past year(Yes) vs Healthy') 
  df <- rbind(df, rbind(tmp1, tmp2))
}

df <- na.omit(df)
head(df)
table(df$cancer, df$PPI)
median_diff <- aggregate(value~cancer + PPI, data = df, function(x) median(x)) %>% 
  pivot_wider(
    names_from = PPI,
    values_from = value
  ) %>% mutate(diff = `PPI within past year(Yes) vs Healthy`-`PPI within past year(No) vs Healthy`) %>% arrange(diff)
opp <- median_diff$cancer[median_diff$diff<0]
sort(cbind(table(data.obj$meta.dat$icd10_first_3_name_short,data.obj$meta.dat$PPI_day_365))[,"N"])[opp]
df$cancer <- factor(df$cancer, levels = median_diff$cancer)

# prepare table for label
df_n <- melt(cbind(nSam_low, nSam_high)) %>% arrange(X1) %>% mutate(y = rep(0.8, length(c(nSam_low, nSam_high))))
colnames(df_n) <- c('cancer','PPI','n','y')
df_n$PPI <- as.character(df_n$PPI)
df_n$PPI[df_n$PPI=='nSam_high'] <- 'PPI within past year(Yes) vs Healthy'
df_n$PPI[df_n$PPI=='nSam_low'] <- 'PPI within past year(No) vs Healthy'

ggplot(df, aes(x = cancer, y = value, fill = PPI)) +
  geom_boxplot(position = position_dodge(width = 0.8), outlier.size = 0.2, outlier.colour = 'grey50') +
  geom_text(
    data = df_n,
    aes(x = cancer, y = y, label = paste0("n=", n), group = PPI),
    position = position_dodge(width = 0.8),
    color = 'black',
    size = 3,
    angle = 270,      
    vjust = 0.5,    
    hjust = 0.5
  )+
  scale_fill_manual(values = c('PPI within past year(No) vs Healthy'='#0072B2',
                               'PPI within past year(Yes) vs Healthy'='#CC79A7',
                               'HV'='steelblue','Cancer'='orange')) +
  theme_bw() +
  labs(x = "", y = "Weighted UniFrac distance", fill = "") +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, color = 'black'),
    axis.text.y = element_text(color = 'black'), 
    legend.position = 'top'
  )  
ggsave(file = paste0(file_dir,'/Figure/RM_figures/PPI_YN_HV_distance.pdf'), width = 6, height = 6)


# median_diff$cancer <- factor(median_diff$cancer, levels = median_diff$cancer)
# ggplot(median_diff, aes(x = cancer, y = diff, fill = diff > 0)) +
#   geom_bar(stat = 'identity') +
#   scale_fill_manual(values = c("TRUE" = '#CC79A7', "FALSE" = '#0072B2')) +
#   scale_y_continuous(limits = c(-0.06,0.06)) +
#   theme_bw() +
#   labs(x = "", y = "WUniFrac (median PPIY vs HV - median PPIN vs HV)", fill = "") +
#   theme(
#     axis.text.x = element_text(angle = 90, hjust = 1, color = 'black'),
#     axis.text.y = element_text(color = 'black'),
#     legend.position = "none"
#   )
# ggsave(file = paste0(file_dir,'/Figure/RM_figures/PPI_YN_HV_distance_diffBar.pdf'), width = 4, height =5)

load('Data/data.obj.raw.core.RData')
df = cbind(table(data.obj.rff$meta.dat$icd10_first_3_name_short, data.obj.rff$meta.dat$PPI_day_365))
df[rowSums(df>5)==2 & rowSums(df)>15,]

df = cbind(table(data.obj.rff$meta.dat$icd10_first_3_name_short, data.obj.rff$meta.dat$Abx_day_365))
df[rowSums(df>5)==2 & rowSums(df)>15,]

### ======= new Supplementary Table S4 (HV vs Low (ECI==0)/High ECI)======
file_dir <- dirname(rstudioapi::getSourceEditorContext()$path)
file_dir <- dirname(dirname(dirname(dirname(file_dir))))
setwd(file_dir)

##  ECI==0 (lowECI in Cancer VS healthy)
load('Data/data.obj.raw.core.RData')

# Low ECI (ECI=0) in Cancer vs HV
tm = load('Result/ECI_subset0/Control-LowCancer/DAA/ECI22/ECI22_ZicoSeq.Rdata')
low_q <- diff.obj$qv.list$Species
low_FC <- diff.obj$R2.list$Species * sign(diff.obj$coef.list$Species[,'ECI22LowCancer'])
df <- merge(diff.obj$pv.list$Species, low_q, by = 0)%>% column_to_rownames('Row.names') %>% merge(low_FC, by = 0) %>% dplyr::rename(Species = 'Row.names', `Effect size`=Func1)
HV_LowECI <- df %>% inner_join(as.data.frame(data.obj$otu.name)) %>% 
  dplyr::select(-Pvalue) %>% 
  mutate(Comparison='Low ECI(ECI=0) vs Healthy (reference)',`Adjust covariates`=paste0(c('BMI','Sex','Age'), collapse = ','))

# High ECI(ECI>=4) in Cancer vs HV
tm = load('Result/ECI_subset0/Control-HighCancer/DAA/ECI22/ECI22_ZicoSeq.Rdata')
high_q <- diff.obj$qv.list$Species
high_FC <- diff.obj$R2.list$Species * sign(diff.obj$coef.list$Species[,'ECI22HighCancer'])
df <- merge(diff.obj$pv.list$Species, high_q, by = 0)%>% column_to_rownames('Row.names') %>% 
  merge(high_FC, by = 0) %>% dplyr::rename(Species = 'Row.names', `Effect size`=Func1)
HV_HighECI <- df %>% inner_join(as.data.frame(data.obj$otu.name)) %>% 
  dplyr::select(-Pvalue) %>% 
  mutate(Comparison='High ECI(ECI>=4) vs Healthy (reference)', `Adjust covariates`=paste0(c('BMI','Sex','Age'), collapse = ','))

## subCancer X with Low ECI(==0) vs HV
cancers <- list.dirs('Result/ECI_subset0/', recursive = F, full.names = F)
cancers <- cancers[-(grep('Control',cancers))]
ress <- NULL
for(cancer in cancers){
  diff.obj <- NULL
  load(paste0('Result/ECI_subset0/',cancer,'/Control-LowCancer/DAA/ECI22/ECI22_ZicoSeq.Rdata'))
  res <- cbind(diff.obj$qv.list$Species,diff.obj$R2.list$Species * sign(diff.obj$coef.list[[1]][,'ECI22LowCancer'])) %>% cbind(diff.obj$pv.list$Species) %>% as.data.frame
  colnames(res)[2] <- 'EffectSize'
  cat(cancer,':')
  cat(sum(res$Qvalue<0.1),'\n')
  ress <- rbind(ress, 
                res %>% mutate(Comparison = paste0(cancer,'(ECI==0) vs. Healthy (reference)')) %>% 
                  rownames_to_column('Species') %>% inner_join(as.data.frame(data.obj$otu.name)) %>% 
                  dplyr::select(-Pvalue))
}
ress_lowECI <- ress %>% mutate(`Adjust covariates`=paste0(c('BMI','Sex','Age'), collapse = ','))

## subCancer X with High ECI(>=4) vs HV
cancers <- list.dirs('Result/ECI_subset0/', recursive = F, full.names = F)
cancers <- cancers[-(grep('Control',cancers))]
ress <- NULL
for(cancer in cancers){
  diff.obj <- NULL
  load(paste0('Result/ECI_subset0/',cancer,'/Control-HighCancer/DAA/ECI22/ECI22_ZicoSeq.Rdata'))
  res <- cbind(diff.obj$qv.list$Species,diff.obj$R2.list$Species * sign(diff.obj$coef.list[[1]][,'ECI22HighCancer'])) %>% cbind(diff.obj$pv.list$Species) %>% as.data.frame
  colnames(res)[2] <- 'EffectSize'
  cat(cancer,':')
  cat(sum(res$Qvalue<0.1),'\n')
  ress <- rbind(ress, res %>% mutate(Comparison = paste0(cancer,'(ECI>=4) vs. Healthy (reference)')) %>% rownames_to_column('Species') %>% 
                  inner_join(as.data.frame(data.obj$otu.name)) %>% dplyr::select(-Pvalue))
}
ress_highECI <- ress %>% mutate(`Adjust covariates`=paste0(c('BMI','Sex','Age'), collapse = ','))
head(ress_highECI)


wb <- loadWorkbook(paste0(file_dir,"/Code/Submission/Supplementary tables/new Supplementary Table 4.xlsx"))
names(wb)

tab.name <- 'TableS4-1'
if (tab.name %in% names(wb)) removeWorksheet(wb, tab.name)
addWorksheet(wb, tab.name)
writeData(wb, sheet = tab.name, x = HV_LowECI)

tab.name <- 'TableS4-2'
if (tab.name %in% names(wb)) removeWorksheet(wb, tab.name)
addWorksheet(wb, tab.name)
writeData(wb, sheet = tab.name, x = HV_HighECI)

tab.name <- 'TableS4-3'
if (tab.name %in% names(wb)) removeWorksheet(wb, tab.name)
addWorksheet(wb, tab.name)
writeData(wb, sheet = tab.name, x = ress_lowECI)

tab.name <- 'TableS4-4'
if (tab.name %in% names(wb)) removeWorksheet(wb, tab.name)
addWorksheet(wb, tab.name)
writeData(wb, sheet = tab.name, x = ress_highECI)
saveWorkbook(wb, paste0(wd,"/Code/Submission/Supplementary tables/new Supplementary Table 4.xlsx"), overwrite = TRUE)




### ======= new Supplementary Table S6 (HV vs Abx_day_365|PPI_day_365 (Y)/(N) )======
file_dir <- dirname(rstudioapi::getSourceEditorContext()$path)
file_dir <- dirname(dirname(dirname(dirname(file_dir))))
setwd(file_dir)

load('Data/data.obj.raw.core.RData')

# Abx yes in Cancer vs HV
tm = load('Result/Abx_day_365_subset/Control-Y/DAA/Abx_day_365/Abx_day_365_ZicoSeq.Rdata')
low_q <- diff.obj$qv.list$Species
low_FC <- diff.obj$R2.list$Species * sign(diff.obj$coef.list$Species[,'Abx_day_365Y'])
df <- merge(diff.obj$pv.list$Species, low_q, by = 0)%>% column_to_rownames('Row.names') %>% merge(low_FC, by = 0) %>% dplyr::rename(Species = 'Row.names', `Effect size`=Func1)
HV_AbxY <- df %>% inner_join(as.data.frame(data.obj$otu.name)) %>% 
  dplyr::select(-Pvalue) %>% 
  mutate(Comparison='Antibiotics within past year(Yes) vs Healthy (reference)',`Adjust covariates`=paste0(c('BMI','Sex','Age'), collapse = ','))

# Abx no in Cancer vs HV
tm = load('Result/Abx_day_365_subset/Control-N/DAA/Abx_day_365/Abx_day_365_ZicoSeq.Rdata')
high_q <- diff.obj$qv.list$Species
high_FC <- diff.obj$R2.list$Species * sign(diff.obj$coef.list$Species[,'Abx_day_365N'])
df <- merge(diff.obj$pv.list$Species, high_q, by = 0)%>% column_to_rownames('Row.names') %>% 
  merge(high_FC, by = 0) %>% dplyr::rename(Species = 'Row.names', `Effect size`=Func1)
HV_AbxN <- df %>% inner_join(as.data.frame(data.obj$otu.name)) %>% 
  dplyr::select(-Pvalue) %>% 
  mutate(Comparison='Antibiotics within past year(No) vs Healthy (reference)', `Adjust covariates`=paste0(c('BMI','Sex','Age'), collapse = ','))

## subCancer X with Abx (Y) vs HV
cancers <- list.dirs('Result/Abx_day_365_subset/', recursive = F, full.names = F)
cancers <- cancers[-(grep('Control',cancers))]
ress <- NULL
for(cancer in cancers){
  diff.obj <- NULL
  load(paste0('Result/Abx_day_365_subset/',cancer,'/Control-Y/DAA/Abx_day_365/Abx_day_365_ZicoSeq.Rdata'))
  res <- cbind(diff.obj$qv.list$Species,diff.obj$R2.list$Species * sign(diff.obj$coef.list[[1]][,'Abx_day_365Y'])) %>% cbind(diff.obj$pv.list$Species) %>% as.data.frame
  colnames(res)[2] <- 'EffectSize'
  cat(cancer,':')
  cat(sum(res$Qvalue<0.1),'\n')
  ress <- rbind(ress, 
                res %>% mutate(Comparison = paste0(cancer,'-Antibiotics within past year(Yes)) vs. Healthy (reference)')) %>% 
                  rownames_to_column('Species') %>% inner_join(as.data.frame(data.obj$otu.name)) %>% 
                  dplyr::select(-Pvalue))
}
ress_AbxY <- ress %>% mutate(`Adjust covariates`=paste0(c('BMI','Sex','Age'), collapse = ','))

## subCancer X Abx (N) vs HV
cancers <- list.dirs('Result/Abx_day_365_subset/', recursive = F, full.names = F)
cancers <- cancers[-(grep('Control',cancers))]
ress <- NULL
for(cancer in cancers){
  diff.obj <- NULL
  load(paste0('Result/Abx_day_365_subset/',cancer,'/Control-N/DAA/Abx_day_365/Abx_day_365_ZicoSeq.Rdata'))
  res <- cbind(diff.obj$qv.list$Species,diff.obj$R2.list$Species * sign(diff.obj$coef.list[[1]][,'Abx_day_365N'])) %>% cbind(diff.obj$pv.list$Species) %>% as.data.frame
  colnames(res)[2] <- 'EffectSize'
  cat(cancer,':')
  cat(sum(res$Qvalue<0.1),'\n')
  ress <- rbind(ress, 
                res %>% mutate(Comparison = paste0(cancer,'-Antibiotics within past year(No)) vs. Healthy (reference)')) %>% 
                  rownames_to_column('Species') %>% inner_join(as.data.frame(data.obj$otu.name)) %>% 
                  dplyr::select(-Pvalue))
}
ress_AbxN <- ress %>% mutate(`Adjust covariates`=paste0(c('BMI','Sex','Age'), collapse = ','))




# PPI yes in Cancer vs HV
tm = load('Result/PPI_day_365_subset/Control-Y/DAA/PPI_day_365/PPI_day_365_ZicoSeq.Rdata')
low_q <- diff.obj$qv.list$Species
low_FC <- diff.obj$R2.list$Species * sign(diff.obj$coef.list$Species[,'PPI_day_365Y'])
df <- merge(diff.obj$pv.list$Species, low_q, by = 0)%>% column_to_rownames('Row.names') %>% merge(low_FC, by = 0) %>% dplyr::rename(Species = 'Row.names', `Effect size`=Func1)
HV_PPIY <- df %>% inner_join(as.data.frame(data.obj$otu.name)) %>% 
  dplyr::select(-Pvalue) %>% 
  mutate(Comparison='PPI within past year(Yes) vs Healthy (reference)',`Adjust covariates`=paste0(c('BMI','Sex','Age'), collapse = ','))

# PPI no in Cancer vs HV
tm = load('Result/PPI_day_365_subset/Control-N/DAA/PPI_day_365/PPI_day_365_ZicoSeq.Rdata')
high_q <- diff.obj$qv.list$Species
high_FC <- diff.obj$R2.list$Species * sign(diff.obj$coef.list$Species[,'PPI_day_365N'])
df <- merge(diff.obj$pv.list$Species, high_q, by = 0)%>% column_to_rownames('Row.names') %>% 
  merge(high_FC, by = 0) %>% dplyr::rename(Species = 'Row.names', `Effect size`=Func1)
HV_PPIN <- df %>% inner_join(as.data.frame(data.obj$otu.name)) %>% 
  dplyr::select(-Pvalue) %>% 
  mutate(Comparison='Antibiotics within past year(No) vs Healthy (reference)', `Adjust covariates`=paste0(c('BMI','Sex','Age'), collapse = ','))

## subCancer X with PPI (Y) vs HV
cancers <- list.dirs('Result/PPI_day_365_subset/', recursive = F, full.names = F)
cancers <- cancers[-(grep('Control',cancers))]
ress <- NULL
for(cancer in cancers){
  diff.obj <- NULL
  load(paste0('Result/PPI_day_365_subset/',cancer,'/Control-Y/DAA/PPI_day_365/PPI_day_365_ZicoSeq.Rdata'))
  res <- cbind(diff.obj$qv.list$Species,diff.obj$R2.list$Species * sign(diff.obj$coef.list[[1]][,'PPI_day_365Y'])) %>% cbind(diff.obj$pv.list$Species) %>% as.data.frame
  colnames(res)[2] <- 'EffectSize'
  cat(cancer,':')
  cat(sum(res$Qvalue<0.1),'\n')
  ress <- rbind(ress, 
                res %>% mutate(Comparison = paste0(cancer,'-Antibiotics within past year(Yes)) vs. Healthy (reference)')) %>% 
                  rownames_to_column('Species') %>% inner_join(as.data.frame(data.obj$otu.name)) %>% 
                  dplyr::select(-Pvalue))
}
ress_PPIY <- ress %>% mutate(`Adjust covariates`=paste0(c('BMI','Sex','Age'), collapse = ','))

## subCancer X PPI (N) vs HV
cancers <- list.dirs('Result/PPI_day_365_subset/', recursive = F, full.names = F)
cancers <- cancers[-(grep('Control',cancers))]
ress <- NULL
for(cancer in cancers){
  diff.obj <- NULL
  load(paste0('Result/PPI_day_365_subset/',cancer,'/Control-N/DAA/PPI_day_365/PPI_day_365_ZicoSeq.Rdata'))
  res <- cbind(diff.obj$qv.list$Species,diff.obj$R2.list$Species * sign(diff.obj$coef.list[[1]][,'PPI_day_365N'])) %>% cbind(diff.obj$pv.list$Species) %>% as.data.frame
  colnames(res)[2] <- 'EffectSize'
  cat(cancer,':')
  cat(sum(res$Qvalue<0.1),'\n')
  ress <- rbind(ress, 
                res %>% mutate(Comparison = paste0(cancer,'-Antibiotics within past year(No)) vs. Healthy (reference)')) %>% 
                  rownames_to_column('Species') %>% inner_join(as.data.frame(data.obj$otu.name)) %>% 
                  dplyr::select(-Pvalue))
}
ress_PPIN <- ress %>% mutate(`Adjust covariates`=paste0(c('BMI','Sex','Age'), collapse = ','))


wb <- loadWorkbook(paste0(file_dir,"/Code/Submission/Supplementary tables/new Supplementary Table 6.xlsx"))
names(wb)

tab.name <- 'TableS6-1'
if (tab.name %in% names(wb)) removeWorksheet(wb, tab.name)
addWorksheet(wb, tab.name)
writeData(wb, sheet = tab.name, x = HV_AbxY)

tab.name <- 'TableS6-2'
if (tab.name %in% names(wb)) removeWorksheet(wb, tab.name)
addWorksheet(wb, tab.name)
writeData(wb, sheet = tab.name, x = HV_AbxN)

tab.name <- 'TableS6-3'
if (tab.name %in% names(wb)) removeWorksheet(wb, tab.name)
addWorksheet(wb, tab.name)
writeData(wb, sheet = tab.name, x = ress_AbxY)

tab.name <- 'TableS6-4'
if (tab.name %in% names(wb)) removeWorksheet(wb, tab.name)
addWorksheet(wb, tab.name)
writeData(wb, sheet = tab.name, x = ress_AbxN)


tab.name <- 'TableS6-5'
if (tab.name %in% names(wb)) removeWorksheet(wb, tab.name)
addWorksheet(wb, tab.name)
writeData(wb, sheet = tab.name, x = HV_PPIY)

tab.name <- 'TableS6-6'
if (tab.name %in% names(wb)) removeWorksheet(wb, tab.name)
addWorksheet(wb, tab.name)
writeData(wb, sheet = tab.name, x = HV_PPIN)

tab.name <- 'TableS6-7'
if (tab.name %in% names(wb)) removeWorksheet(wb, tab.name)
addWorksheet(wb, tab.name)
writeData(wb, sheet = tab.name, x = ress_PPIY)

tab.name <- 'TableS6-8'
if (tab.name %in% names(wb)) removeWorksheet(wb, tab.name)
addWorksheet(wb, tab.name)
writeData(wb, sheet = tab.name, x = ress_PPIN)

saveWorkbook(wb, paste0(wd,"/Code/Submission/Supplementary tables/new Supplementary Table 6.xlsx"), overwrite = TRUE)



### === Compare Elix_score continuous variable vs binary =========
file_dir <- dirname(rstudioapi::getSourceEditorContext()$path)
file_dir <- dirname(dirname(dirname(dirname(file_dir))))
setwd(file_dir)

load('Result/CancerOnly/data.obj.wk.RData')
table(data.obj.rff$meta.dat$Elix_score2_01)
table(data.obj.rff$meta.dat$Elix_score2_med)
table(data.obj.rff$meta.dat$Elix_score3)


# alpha
load('Result/CancerOnly/Alpha/Elix_score/Alpha.RData')
fit2$res['Shannon',, drop =F] # adjusted

load('Result/CancerOnly/Alpha/Elix_score2_med/Alpha.RData')
fit2$res['Shannon',, drop =F] # adjusted, ref=low

load('Result/CancerOnly/Alpha/Elix_score2_01/Alpha.RData')
fit2$res['Shannon',, drop =F] # adjusted

load('Result/CancerOnly/Alpha/Elix_score3/Alpha.RData')
fit2$res['Shannon',, drop =F] # adjusted
summary(fit2$fitted.obj$Shannon)

# Beta
load('Result/CancerOnly/Beta/Elix_score/R2_pvalue.RData')
pv.adj.mat['WUniFrac',]

load('Result/CancerOnly/Beta/Elix_score2_med/R2_pvalue.RData')
pv.adj.mat['WUniFrac',]

load('Result/CancerOnly/Beta/Elix_score2_01/R2_pvalue.RData')
pv.adj.mat['WUniFrac',]

load('Result/CancerOnly/Beta/Elix_score3/R2_pvalue.RData')
pv.adj.mat['WUniFrac',]

# DAA 
load('Result/CancerOnly/DAA/Elix_score/Elix_score_ZicoSeq.Rdata')
sum(diff.obj$qv.list$Species < 0.1)

load('Result/CancerOnly/DAA/Elix_score2_med/Elix_score2_med_ZicoSeq.Rdata')
sum(diff.obj$qv.list$Species < 0.1)

load('Result/CancerOnly/DAA/Elix_score2_01/Elix_score2_01_ZicoSeq.Rdata')
sum(diff.obj$qv.list$Species < 0.1)

load('Result/CancerOnly/DAA/Elix_score3/Elix_score3_ZicoSeq.Rdata')
sum(diff.obj$qv.list$Species < 0.1)




### === low ECI vs high ECI vs HV Alpha&Beta dot plot =====
file_dir <- dirname(rstudioapi::getSourceEditorContext()$path)
file_dir <- dirname(dirname(dirname(dirname(file_dir))))

wd <- file_dir
rd <- paste0(file_dir,'/Result/')

setwd(rd)
dirs <- list.dirs('ECI_subset/', recursive = F, full.names = F)
dirs <- dirs[!grepl('Control',dirs)]

df <- est <- NULL
for(dir in dirs){
  for(level in c('LowCancer','MediumCancer','HighCancer')){
    load(paste0('ECI_subset/',dir,'/Control-',level,'/Alpha/ECI2/Alpha.RData'))
    tmp = coef(summary(fit2$fitted.obj$Shannon))
    tmp2 <- (tmp[rownames(tmp)[grep('ECI2',rownames(tmp))],,drop =F])
    est <- rbind(est,cbind.data.frame(tmp2) %>% mutate(ECI = gsub('Cancer','ECI',level), comparison=paste0(dir)))
    df <- rbind(df, as.data.frame(fit2$res['Shannon',, drop =F]) %>% mutate(ECI = gsub('Cancer','ECI',level), comparison=paste0(dir)))
  }
}
rownames(df) <- NULL
head(df)
head(est)
df$R2[df$R2<0] <- 0
df$ECI <- factor(df$ECI, levels = c("LowECI", "MediumECI", "HighECI"))


df_plot <- df %>% 
  # dplyr::filter(ECI!="MediumECI") %>%
  group_by(ECI) %>% 
  mutate(
    direction = ifelse(coef > 0, 'positive', 'negative'),
    P = p.adjust(P, method = "BH")
  ) %>% 
  mutate(P.col = ifelse(P < 0.1, "sig", "no")) %>% 
  ungroup() %>%
  mutate(ECI = factor(ECI, levels = c("LowECI", "MediumECI", "HighECI")))

df_plot <- df_plot %>%
  mutate(
    comparison_num = as.numeric(factor(reorder(comparison, R2, mean))),
    y_pos = case_when(
      ECI == "LowECI"    ~ comparison_num - 0.25,
      ECI == "MediumECI" ~ comparison_num,
      ECI == "HighECI"   ~ comparison_num + 0.25
    )
  )

cancer_labels <- df_plot %>%
  distinct(comparison, comparison_num) %>%
  arrange(comparison_num)
df_plot1 =df_plot



df <- NULL
for(dir in dirs){
  for(level in c('LowCancer','MediumCancer','HighCancer')){
    load(paste0('ECI_subset/',dir,'/Control-',level,'/Beta/ECI2/R2_pvalue.RData'))
    df <- rbind(df, as.data.frame(cbind(R2=r2.adj.mat['WUniFrac',],P=pv.adj.mat['WUniFrac',])) %>% mutate(ECI = gsub('Cancer','ECI',level), comparison=paste0(dir)))
  }
}
rownames(df) <- NULL


df_plot <- df %>% 
  # dplyr::filter(ECI!="MediumECI") %>%
  group_by(ECI) %>% 
  mutate(
    P = p.adjust(P, method = "BH")
  ) %>% 
  mutate(P.col = ifelse(P < 0.1, "sig", "no")) %>% 
  ungroup() %>%
  mutate(ECI = factor(ECI, levels = c("LowECI", "MediumECI", "HighECI")))

df_plot <- df_plot %>%
  mutate(
    comparison_num = as.numeric(factor(reorder(comparison, R2, mean))),
    y_pos = case_when(
      ECI == "LowECI"    ~ comparison_num - 0.25,
      ECI == "MediumECI" ~ comparison_num,
      ECI == "HighECI"   ~ comparison_num + 0.25
    )
  )

cancer_labels <- df_plot %>%
  distinct(comparison, comparison_num) %>%
  arrange(comparison_num)


df_plot12 <- rbind(df_plot1 %>% dplyr::select(-coef) %>% mutate(diversity = 'Alpha diversity'), df_plot %>% mutate(direction = 'none', diversity = 'Beta diversity') ) 
df_plot12$ECI <- as.character(df_plot12$ECI)
df_plot12$ECI[df_plot12$ECI=='LowECI'] <- 'LowECI(ECI<=1)'
df_plot12$ECI[df_plot12$ECI=='MediumECI'] <- 'MediumECI(ECI=2,3)'
df_plot12$ECI[df_plot12$ECI=='HighECI'] <- 'HighECI(ECI>=4)'
ggplot(df_plot12,aes(y = y_pos, x = R2)) +
  geom_segment(
    aes(yend = y_pos, x = 0, xend = R2, colour = ECI),
    linewidth = 2
  ) +
  geom_point(
    aes(shape = direction, fill = P.col),
    size = 2, stroke = 0.5
  ) +
  scale_fill_manual(values = c(sig = "black", no = "white")) +
  scale_shape_manual(values = c(positive = 24, negative = 25, none = 16)) +
  scale_colour_manual(values = c(
    "LowECI(ECI<=1)"    = "#E69F00",
    'MediumECI(ECI=2,3)' = "#56B4E9",
    "HighECI(ECI>=4)"   = "#009E73"
  )) +
  scale_y_continuous(
    breaks = cancer_labels$comparison_num,
    labels = cancer_labels$comparison
  ) +
  facet_grid(~diversity, scales ='free') +
  theme_classic() +
  labs(
    x      = "Variability explained (R²), %",
    y      = "",
    fill   = "",
    shape  = "",
    colour = ""
  ) +
  theme(
    panel.grid   = element_blank(),
    axis.title   = element_text(size = 16, colour = "black"),
    legend.text  = element_text(size = 16, colour = "black"),
    strip.text  = element_text(size = 16, colour = "black"),
    axis.text    = element_text(size = 16, colour = "black")
  ) +
  guides(
    colour = guide_legend(title = "", order = 0),
    fill   = guide_legend(title = "", order = 1,
                          override.aes = list(shape = 21)),
    shape  = guide_legend(title = "", order = 2,
                          override.aes = list(fill = "white"))
  )



### ====== pancancer patients vs abx-N, PPI-N =======
file_dir <- dirname(rstudioapi::getSourceEditorContext()$path)
file_dir <- dirname(dirname(dirname(dirname(file_dir))))
rd <- paste0(file_dir, '/Result/')

plot_colours <- c(
  "Antibiotics within past year(Y)" = "#E8A020",
  "Antibiotics within past year(N)" = "#F5D08A",
  "PPI within past year(Y)"         = "#CC79A7",
  "PPI within past year(N)"         = "#E8BAD4",
  "Cancer"="#009E73"
)
dist.name <- 'WUniFrac'
# ── Load pancancer ──────────────────────────────────────────────────────
setwd(rd)
load('PanCancer/Alpha/Group/Alpha.RData')
pancancer_alpha <- as.data.frame(fit2$res['Shannon', , drop = FALSE]) %>%
  mutate(ECI = 'Cancer',
         comparison = 'Healthy vs. Cancer') %>% 
  mutate(direction = ifelse(coef > 0, 'positive', 'negative'),
         comparison_num = as.numeric(factor(reorder(comparison, R2, mean))),
         P.col = ifelse(P<0.1, 'sig','no'))

load('PanCancer/Beta/Group/R2_pvalue.RData')
pancancer_beta <- as.data.frame(cbind(R2 = r2.adj.mat[dist.name, ], P = pv.adj.mat[dist.name, ])) %>%
  mutate(ECI = 'Cancer',
         comparison = 'Healthy vs. Cancer') %>% 
  mutate(direction = 'none',
         comparison_num = as.numeric(factor(reorder(comparison, R2, mean))),
         P.col = ifelse(P<0.1, 'sig','no'))

# ── Load Abx-N ─────────────────────────────────────────────────────────────
load('Abx_day_365_subset/Control-N/Alpha/Abx_day_365/Alpha.RData')
abxN_alpha <- as.data.frame(fit2$res['Shannon', , drop = FALSE]) %>%
  mutate(ECI = "Antibiotics within past year(N)",
         comparison = 'Healthy vs. Antibiotics within past year(N)') %>% 
  mutate(direction = ifelse(coef > 0, 'positive', 'negative'), P.col = ifelse(P<0.1, 'sig','no')) 

load('Abx_day_365_subset/Control-N/Beta/Abx_day_365/R2_pvalue.RData')
abxN_beta <- as.data.frame(cbind(R2 = r2.adj.mat[dist.name, ], P = pv.adj.mat[dist.name, ]))  %>%
  mutate(ECI = "Antibiotics within past year(N)",
         comparison = 'Healthy vs. Antibiotics within past year(N)') %>% 
  mutate(direction = 'none', P.col = ifelse(P<0.1, 'sig','no')) 

# ── Load Abx-Y ─────────────────────────────────────────────────────────────
load('Abx_day_365_subset/Control-Y/Alpha/Abx_day_365/Alpha.RData')
abxY_alpha <- as.data.frame(fit2$res['Shannon', , drop = FALSE]) %>%
  mutate(ECI = "Antibiotics within past year(Y)",
         comparison = 'Healthy vs. Antibiotics within past year(Y)') %>% 
  mutate(direction = ifelse(coef > 0, 'positive', 'negative'), P.col = ifelse(P<0.1, 'sig','no')) 

load('Abx_day_365_subset/Control-Y/Beta/Abx_day_365/R2_pvalue.RData')
abxY_beta <- as.data.frame(cbind(R2 = r2.adj.mat[dist.name, ], P = pv.adj.mat[dist.name, ]))  %>%
  mutate(ECI = "Antibiotics within past year(Y)",
         comparison = 'Healthy vs. Antibiotics within past year(Y)') %>% 
  mutate(direction = 'none', P.col = ifelse(P<0.1, 'sig','no')) 


# ── Load PPI-N ─────────────────────────────────────────────────────────────
load('PPI_day_365_subset/Control-N/Alpha/PPI_day_365/Alpha.RData')
PPIN_alpha <- as.data.frame(fit2$res['Shannon', , drop = FALSE]) %>%
  mutate(ECI = "PPI within past year(N)",
         comparison = 'Healthy vs. PPI within past year(N)') %>% 
  mutate(direction = ifelse(coef > 0, 'positive', 'negative'), P.col = ifelse(P<0.1, 'sig','no')) 

load('PPI_day_365_subset/Control-N/Beta/PPI_day_365/R2_pvalue.RData')
PPIN_beta <- as.data.frame(cbind(R2 = r2.adj.mat[dist.name, ], P = pv.adj.mat[dist.name, ]))  %>%
  mutate(ECI = "PPI within past year(N)",
         comparison = 'Healthy vs. PPI within past year(N)') %>% 
  mutate(direction = 'none', P.col = ifelse(P<0.1, 'sig','no')) 

# ── Load PPI-Y ─────────────────────────────────────────────────────────────
load('PPI_day_365_subset/Control-Y/Alpha/PPI_day_365/Alpha.RData')
PPIY_alpha <- as.data.frame(fit2$res['Shannon', , drop = FALSE]) %>%
  mutate(ECI = "PPI within past year(Y)",
         comparison = 'Healthy vs. PPI within past year(Y)') %>% 
  mutate(direction = ifelse(coef > 0, 'positive', 'negative'), P.col = ifelse(P<0.1, 'sig','no')) 

load('PPI_day_365_subset/Control-Y/Beta/PPI_day_365/R2_pvalue.RData')
PPIY_beta <- as.data.frame(cbind(R2 = r2.adj.mat[dist.name, ], P = pv.adj.mat[dist.name, ]))  %>%
  mutate(ECI = "PPI within past year(Y)",
         comparison = 'Healthy vs. PPI within past year(Y)') %>% 
  mutate(direction = 'none', P.col = ifelse(P<0.1, 'sig','no')) 


# ── Assemble data ──────────────────────────────────────────────────────────────
abx_ppi <- bind_rows(
  pancancer_alpha %>% dplyr::select(-coef) %>% mutate(diversity = 'Alpha diversity'),
  pancancer_beta                    %>% mutate(diversity = 'Beta diversity'),
  abxN_alpha  %>% dplyr::select(-coef) %>% mutate(diversity = 'Alpha diversity'),
  abxN_beta                   %>% mutate(diversity = 'Beta diversity'),
  # abxY_alpha  %>% dplyr::select(-coef) %>% mutate(diversity = 'Alpha diversity'),
  # abxY_beta                   %>% mutate(diversity = 'Beta diversity'),
  
  PPIN_alpha  %>% dplyr::select(-coef) %>% mutate(diversity = 'Alpha diversity'),
  PPIN_beta                    %>% mutate(diversity = 'Beta diversity')#,
  # PPIY_alpha  %>% dplyr::select(-coef) %>% mutate(diversity = 'Alpha diversity'),
  # PPIY_beta                    %>% mutate(diversity = 'Beta diversity')
  
) %>%
  mutate(y_pos = case_when(
    ECI == "PPI within past year(N)"         ~ 1,
    # ECI == "PPI within past year(Y)"         ~ 1.5,
    ECI == "Antibiotics within past year(N)" ~ 2,
    # ECI == "Antibiotics within past year(Y)" ~ 2.5,
    ECI == "Cancer"                          ~ 3
  ))
rownames(abx_ppi) <- NULL

y_labels <- abx_ppi %>%
  distinct(ECI, y_pos) %>%
  arrange(y_pos)

ggplot(abx_ppi, aes(y = y_pos, x = R2)) +
  geom_segment(aes(yend = y_pos, x = 0, xend = R2, colour = ECI),
               linewidth = 4) +
  geom_point(aes(shape = direction, fill = P.col),
             size = 4, stroke = 0.5) +
  scale_fill_manual(values  = c(sig = "black", no = "white")) +
  scale_shape_manual(values = c(positive = 24, negative = 25, none = 16)) +
  scale_colour_manual(values = plot_colours) +
  scale_y_continuous(
    breaks = y_labels$y_pos,
    labels = y_labels$ECI,
    limits = c(min(y_labels$y_pos) - 1, max(y_labels$y_pos) + 1),
    expand = c(0, 0)
  ) +
  facet_grid(~diversity, scales = 'free') +
  theme_classic() +
  labs(x = "Variability explained (R²), %", y = "", fill = "", shape = "", colour = "") +
  theme(
    panel.grid  = element_blank(),
    axis.title  = element_text(size = 16, colour = "black"),
    legend.text = element_text(size = 16, colour = "black"),
    strip.text  = element_text(size = 16, colour = "black"),
    axis.text   = element_text(size = 16, colour = "black")
  ) +
  guides(
    colour = "none",
    fill   = "none",
    shape  = guide_legend(title = "", order = 1,
                          override.aes = list(fill = "black"))
  )

setwd(file_dir)
ggsave(file = paste0('Figure/RM_figures/ppiN_abxN_Pancancer_',dist.name,'.pdf'), width = 10, height =3)




### ====== pancancer patients (no subset cancer, ECI low=0) cancer n>15=======
file_dir <- dirname(rstudioapi::getSourceEditorContext()$path)
file_dir <- dirname(dirname(dirname(dirname(file_dir))))
rd <- paste0(file_dir, '/Result/')

bh_adjust <- function(df) {
  df %>% mutate(P = p.adjust(P, method = "BH"),
                P.col = ifelse(P < 0.1, "sig", "no"))
}

plot_colours <- c(
  "LowECI(ECI=0)"                  = "#80C9A4",
  "LowECI(ECI<=1)"                  = "#80C9A4",
  "HighECI(ECI>=4)"                 = "#009E73"
)

y_positions <- c(
  "HighECI(ECI>=4)"                 = 6,
  "LowECI(ECI=0)"                  = 4.5,
  "LowECI(ECI<=1)"                  = 4.5
)

dir = 'ECI_subset' # Cancer-lowECI:01; Cancer-highECI:>=4; HV: ECI<=4; including all cancer patients
dir = 'ECI_subset4' # Cancer-lowECI:01; Cancer-highECI:>=4; HV: ECI<=4; including cancers with n>15 --> this is the figure
dir = 'ECI_subset5' # Cancer-lowECI:01; Cancer-highECI:>=4; HV: ECI=0; including cancers with n>15
dir = 'ECI_subset6' # Cancer-lowECI:01; Cancer-highECI:>=4; HV: ECI=0; including all cancer patients

dir = 'ECI_subset0' # Cancer-lowECI:0; Cancer-highECI:>=4; HV: ECI<=4; including all cancer patients
dir = 'ECI_subset1' # Cancer-lowECI:0; Cancer-highECI:>=4; HV: ECI<=4; including cancers with n>15
dir = 'ECI_subset2' # Cancer-lowECI:0; Cancer-highECI:>=4; HV: ECI=0; including cancers with n>15
dir = 'ECI_subset3' # Cancer-lowECI:0; Cancer-highECI:>=4; HV: ECI=0; including all cancer patients

for(distance in c('WUniFrac')){#c('UniFrac', 'GUniFrac', 'WUniFrac', 'BC')
  if(dir=='ECI_subset'){
    variable = 'ECI2'
  }else{
    variable = 'ECI22'
  }
  
  if(dir %in% c('ECI_subset','ECI_subset4','ECI_subset5','ECI_subset6')){
    low_level = "LowECI(ECI<=1)"
  }else{
    low_level = "LowECI(ECI=0)"
  }
  # ── Load ECI (pan-cancer) ──────────────────────────────────────────────────────
  load_eci_alpha <- function(rd) {
    df <- est <- NULL
    for (level in c('LowCancer', 'HighCancer')) {
      load(paste0(rd, dir,'/Control-', level, '/Alpha/',variable,'/Alpha.RData'))
      df <- rbind(df,
                  as.data.frame(fit2$res['Shannon', , drop = FALSE]) %>%
                    mutate(ECI        = gsub('Cancer', 'ECI', level),
                           comparison = gsub('-', ' vs. ', paste0('Control-', level)))
      )
    }
    df$R2[df$R2 < 0] <- 0
    df %>%
      group_by(ECI) %>%
      bh_adjust() %>%
      mutate(direction = ifelse(coef > 0, 'positive', 'negative')) %>%
      ungroup() %>%
      mutate(ECI = factor(ECI, levels = c("LowECI", "HighECI")),
             comparison_num = as.numeric(factor(reorder(comparison, R2, mean))),
             y_pos = case_when(
               ECI == "LowECI"  ~ comparison_num - 1,
               ECI == "HighECI" ~ comparison_num + 1
             ))
  }
  
  load_eci_beta <- function(rd) {
    df <- NULL
    for (level in c('LowCancer', 'HighCancer')) {
      load(paste0(rd, dir,'/Control-', level, '/Beta/',variable,'/R2_pvalue.RData'))
      df <- rbind(df,
                  as.data.frame(cbind(R2 = r2.adj.mat[distance, ], P = pv.adj.mat[distance, ])) %>%
                    mutate(ECI        = gsub('Cancer', 'ECI', level),
                           comparison = paste0('Control vs. ', level))
      )
    }
    df$R2[df$R2 < 0] <- 0
    df %>%
      group_by(ECI) %>%
      bh_adjust() %>%
      ungroup() %>%
      mutate(ECI = factor(ECI, levels = c("LowECI", "HighECI")),
             comparison_num = as.numeric(factor(reorder(comparison, R2, mean))),
             direction = 'none',
             y_pos = case_when(
               ECI == "LowECI"  ~ comparison_num - 1,
               ECI == "HighECI" ~ comparison_num + 1
             ))
  }
  
  
  # ── Assemble data ──────────────────────────────────────────────────────────────
  eci_alpha <- load_eci_alpha(rd)
  eci_beta  <- load_eci_beta(rd)
  
  
  eci <- rbind(
    eci_alpha %>% dplyr::select(-coef) %>% mutate(diversity = 'Alpha diversity'),
    eci_beta                    %>% mutate(diversity = 'Beta diversity')
  ) %>%
    mutate(ECI = recode(as.character(ECI),
                        LowECI    = low_level,
                        HighECI   = "HighECI(ECI>=4)"
    ))
  eci <- eci %>%
    mutate(y_pos = case_when(
      ECI == "HighECI(ECI>=4)" ~ 2.5,
      ECI == low_level   ~ 1.5
    ))
  
  
  # ── Plot ───────────────────────────────────────────────────────────────────────
  y_labels <- eci %>% distinct(ECI, y_pos) %>% arrange(y_pos)
  
  p1 <- ggplot(eci, aes(y = y_pos, x = R2)) +
    geom_segment(aes(yend = y_pos, x = 0, xend = R2, colour = ECI),
                 linewidth = 4) +
    geom_point(aes(shape = direction, fill = P.col),
               size = 4, stroke = 0.5) +
    scale_fill_manual(values  = c(sig = "black", no = "white")) +
    scale_shape_manual(values = c(positive = 24, negative = 25, none = 16)) +
    scale_colour_manual(values = plot_colours) +
    scale_y_continuous(
      breaks = y_labels$y_pos,
      labels = y_labels$ECI,
      limits = c(1, 3),  
      expand = c(0, 0)
    )+
    facet_grid(~diversity, scales = 'free') +
    theme_classic() +
    labs(x = "Variability explained (R²), %", y = "", fill = "", shape = "", colour = "") +
    theme(
      panel.grid  = element_blank(),
      axis.title  = element_text(size = 16, colour = "black"),
      legend.text = element_text(size = 16, colour = "black"),
      strip.text  = element_text(size = 16, colour = "black"),
      axis.text   = element_text(size = 16, colour = "black")
    ) +
    guides(
      colour = "none",
      fill   = "none",
      shape  = guide_legend(title = "", order = 1,
                            override.aes = list(fill = "black"))
    )
  
  setwd(file_dir)
  ggsave(filename = paste0('Figure/RM_figures/ECI-low-high_',distance,'_',dir,'.png'),
         plot     = p1,
         width    = 15,
         height   = 2.5,
         dpi = 300
  )
}



## =====Abx N vs HV Species differences ======
load('Result/Abx_day_365_subset/Control-N/data.obj.wk.RData')
table(data.obj$meta.dat$Abx_day_365)
load('Result/Abx_day_365_subset/Control-N/DAA/Abx_day_365/Abx_day_365_ZicoSeq.Rdata')
abxN_hv <- cbind.data.frame(diff.obj$qv.list$Species, diff.obj$R2.list$Species * sign(diff.obj$coef.list[[1]][,'Abx_day_365N']))
colnames(abxN_hv)[2] <- 'EffectSize'
table(abxN_hv$Qvalue<0.1,abxN_hv$EffectSize<0)
table(abxN_hv$Qvalue<0.05,abxN_hv$EffectSize<0)

## =====PPI N vs HV ======
load('Result/PPI_day_365_subset/Control-N/data.obj.wk.RData')
table(data.obj$meta.dat$PPI_day_365)
load('Result/PPI_day_365_subset/Control-N/DAA/PPI_day_365/PPI_day_365_ZicoSeq.Rdata')
ppiN_hv <- cbind.data.frame(diff.obj$qv.list$Species, diff.obj$R2.list$Species * sign(diff.obj$coef.list[[1]][,'PPI_day_365N']))
colnames(ppiN_hv)[2] <- 'EffectSize'
table(ppiN_hv$Qvalue<0.1,ppiN_hv$EffectSize<0)
table(ppiN_hv$Qvalue<0.05,ppiN_hv$EffectSize<0)



