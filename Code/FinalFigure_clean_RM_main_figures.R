#use palette.colors(palette = "Okabe-Ito")
sel_colors <- palette.colors(palette = "Okabe-Ito")[2:8]


source('~/Library/Mobile Documents/com~apple~CloudDocs/mforge_clean/Code/Stats.R')
try(load_package())


#library(ComplexHeatmap, lib.loc = .libPaths()[2])
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


mfd <- "~/Library/Mobile Documents/com~apple~CloudDocs/mforge_clean/ManuscriptFigures/"
fd <- "~/Library/Mobile Documents/com~apple~CloudDocs/mforge_clean/Figure/"
wd <- "~/Library/Mobile Documents/com~apple~CloudDocs/mforge_clean/"
rd <- "~/Library/Mobile Documents/com~apple~CloudDocs/mforge_clean/Result/"


setwd(wd)
source('Code/summary_plot_func.R')

setwd(rd)
load('CancerOnly/data.obj.wk.RData')

for(i in colnames(data.obj$meta.dat)[grep('^Elixhauser_',colnames(data.obj$meta.dat))][-1]){
  if(table(data.obj$meta.dat[,i])[2]<15){
    cat(i,'\n')
  }
}


covars <- c("Batch","Bristol_score","BMI", "Age", "Sex","Cancer_class","Metastasis","PPI_day_365", "Abx_day_365",
            "Abx_last_month","PPI_last_month","Elix_score","Sample_season","Urban","icd10_first_3_name")
blood.names <- c('Erythrocytes','Hematocrit','Neutrophils','MCV','Hemoglobin','Leukocytes','Platelet.Count',
                 'Neutrophils_cat','Platelet.Count_cat','Hemoglobin_cat',
                 "neut_neutropenia2","bone_marrow_suppression2","Hb_anemia2","Pl_thrombocytopenia2","neut_neutropenia_c","Hb_anemia_c","Pl_thrombocytopenia_c")

elix.names <- c("Elixhauser_CHF","Elixhauser_Arrhythmia","Elixhauser_Valvular","Elixhauser_PHTN","Elixhauser_PVD","Elixhauser_HTN",
                "Elixhauser_NeuroOther","Elixhauser_Pulmonary","Elixhauser_DM","Elixhauser_DMcx","Elixhauser_Hypothyroid","Elixhauser_Renal","Elixhauser_Liver",      
                "Elixhauser_PUD" ,"Elixhauser_Lymphoma","Elixhauser_Mets","Elixhauser_Tumor","Elixhauser_Rheumatic","Elixhauser_Coagulopathy",
                "Elixhauser_Obesity","Elixhauser_WeightLoss","Elixhauser_FluidsLytes","Elixhauser_Anemia","Elixhauser_Alcohol","Elixhauser_Drugs",
                "Elixhauser_Depression")
charlson.names <- c("Charlson_MI","Charlson_CHF","Charlson_PVD","Charlson_Stroke","Charlson_Pulmonary","Charlson_Rheumatic",  
                    "Charlson_PUD","Charlson_LiverMild","Charlson_DM","Charlson_DMcx","Charlson_Renal","Charlson_Cancer",     
                    "Charlson_LiverSevere","Charlson_Mets")

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

variable_names_c <- c(charlson.names, "Elix_score", "Charlson_score")
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
detailed_descriptions_c <- c(
  "Myocardial infarction", # Charlson_MI
  "Congestive heart failure", # Charlson_CHF
  "Peripheral vascular disease", # Charlson_PVD
  "Cerebrovascular disease", # Charlson_Stroke
  # "Dementia", # Charlson_Dementia
  "Chronic pulmonary disease", # Charlson_Pulmonary
  "Rheumatic disease", # Charlson_Rheumatic
  "Peptic ulcer disease", # Charlson_PUD
  "Mild liver disease", # Charlson_LiverMild
  "Diabetes (no chronic complication)", # Charlson_DM
  "Diabetes (chronic complication)", # Charlson_DMcx
  # "Hemiplegia or paraplegia", # Charlson_Paralysis
  "Renal failure", # Charlson_Renal
  "Any malignancy", # Charlson_Cancer
  "Moderate/severe liver disease", # Charlson_LiverSevere
  "Metastatic solid tumor", # Charlson_Mets
  # "AIDS/HIV infection", # Charlson_HIV
  "ECI score", "CCI score"
)
names_map <- setNames(detailed_descriptions, variable_names)
names_map_c <- setNames(detailed_descriptions_c, variable_names_c)


setwd(wd)
load('Data/cancerNames_20141101.Rdata')
names(cancer.dir) <- cancer.type
# names(cancer.dir) <- tolower(gsub("Malignant neoplasm of |Malignant |malignant |Other and ", "", cancer.type))
names(cancer.dir) <- tolower(gsub("Malignant neoplasm of |Malignant |malignant |Other and | neoplasm| neoplasm of| neoplasms|neoplasm of |\\,|and ", "", names(cancer.dir)))

names(cancer.dir) 
#[1] "prostate"                          "bronchus lung"                 "liver intrahepatic bile ducts"        "colon"                             "melanoma of skin"                 
#[6] "bladder"                           "breast"                            "rectum"                            "brain"                             "unspecified skin"                 
#[11] "kidney except renal pelvis"       "esophagus"                         "oropharynx"                        "neuroendocrine tumors"             "multiple myeloma plasma cell" 
#[16] "pancreas"                          "tonsil"                            "corpus uteri"                      "other connective soft tissue"  "lymphoid leukemia"                
#[21] "ovary"                             "stomach"                           "base of tongue" 

cancer_class_names <- names(cancer.dir) 



#-------------------------------------------------------------------------------
## =========Control X vs subCancerX (PERMANOVA R2) =========  
setwd(wd)
load('Data/cancerNames_20141101.Rdata')
names(cancer.dir) <- cancer_class_names

col <- c(brewer.pal(12,'Paired'),brewer.pal(8,'Dark2'), brewer.pal(8,'Accent')[1:6])
# names(col) <- c('Control',names(cancer.dir[match(old.names, cancer.dir)]))
names(col) <- c('healthy',names(cancer.dir))


setwd(fd)
load('subCancerX_Control/r2.p.unadj.mat.Rdata')
dist.names <- c('UniFrac', 'GUniFrac', 'WUniFrac', 'BC')

setwd(mfd)

#make these colors a color scale; i.e. like heat colors
#"#CC79A7", "#0072B2"
#use the euclidian distance for coloring the color scale?


dist <- 'WUniFrac'
R2 <- r2.unadj.mat[,,dist]
R2[lower.tri(R2)] <- t(R2)[lower.tri(R2)]
R2[is.na(R2)] <- 0

obj <- cmdscale(as.dist(R2), k=3, eig=T)
pve <- round(obj$eig[1:3]/sum(abs(obj$eig))*100, 1)
y <- cbind.data.frame(PC1=obj$points[, 1], PC2=obj$points[, 2], PC3=obj$points[, 3])
PC1.lab <- paste0('PC1 (', pve[1], '%)')
PC2.lab <- paste0('PC2 (', pve[2], '%)')
PC3.lab <- paste0('PC3 (', pve[3], '%)')
yy <- y %>% rownames_to_column('cancer')

#rename the cancer classes
yy$cancer <- c("healthy", names(cancer.dir))

yy_temp <- yy
rownames(yy_temp) <- yy_temp$cancer
yy_temp <- yy_temp[,2:4]
diff_vec <- as.data.frame(dist(yy_temp))[1:23,] 

col_fun <- colorRamp2(c(min(diff_vec), max(diff_vec)), c("#0072B2", "#CC79A7"))
col <- c("#0072B2", col_fun(diff_vec))
names(col) <- c('healthy',names(cancer.dir))


png('../Figure/RM_figures/PCoA centroids cancer classes.png', width = 9*300, height = 8.5*300, res=300)

p12 <- ggplot(yy, aes(x = PC1, y = PC2, fill = cancer)) +
  geom_point(size =6, shape = 21) +
  scale_x_continuous(limits = c(min(yy$PC1)*1.2,max(yy$PC1)*1.2))+
  scale_y_continuous(limits = c(min(yy$PC2)*1.2,max(yy$PC2)*1.2))+
  scale_fill_manual(values = col)+
  #ggrepel::geom_text_repel(aes(label = cancer), max.overlaps = Inf, color = "black", size = 5, box.padding = 0.5) +
  ggrepel::geom_text_repel(aes(label = cancer), color = "black", size = 5, min.segment.length = 1, box.padding = 0.5) +
  # scale_shape_manual(values = 1:13)+
  labs(x = PC1.lab, y = PC2.lab, color = '', fill = '') +
  theme_bw() +
  theme(text = element_text(size = 16, color = "black"),
        axis.text = element_text(size =16, color = "black"),
        legend.text = element_text(size = 16, color = "black"),
        axis.title =  element_text(size = 16, color = "black"),
        panel.grid.major = element_blank(),
        legend.position = 'none',
        panel.border = element_rect(size = 1),
        axis.ticks = element_line(size = 1),
        axis.ticks.length = unit(0.2, 'cm'),
        panel.grid.minor = element_blank())

print(p12)

dev.off()



#cancer only, no control in here
dim(data.obj$otu.tab)
#7713 1364
dim(data.obj.rff$otu.tab)
#7233 1357

#for control
#1644 samples. Prior to rarefaction, there were 1651. So removing 7 samples
#7 samples drop out because of rarifying at 20k reads



#-------------------------------------------------------------------------------
## ===== subCancerX vs control PCA at sample levels=====
setwd(wd)
load('Data/cancerNames.Rdata')
names(cancer.dir) <- cancer.type

load(file = 'Data/data.obj.raw.core.RData') 
setwd(mfd)
dist = "WUniFrac"
# pdf('PCoA.pdf', width = 16, height = 6)
# for(dist in "WUniFrac"){
# dev.off()
obj <- cmdscale(dist.obj.rff[[dist]], k=3, eig = T)
pve <- round(obj$eig[1:3]/sum(abs(obj$eig))*100, 1)
eig <- obj$points
colnames(eig) <- c('PC1','PC2','PC3')
data.obj.rff$meta.dat$icd10_first_3_name[(data.obj.rff$meta.dat$Group=='Control')] <- 'Control'
yy <- merge(eig, data.obj.rff$meta.dat[,'icd10_first_3_name',drop =F], by =0) 
xlab <- paste0('PC1 (',pve[1],'%)')
ylab <-  paste0('PC2 (',pve[2],'%)')
zlab <-  paste0('PC3 (',pve[3],'%)')


Y <- data.obj.rff$abund.list$Genus
#dim(Y) #1644
N <- colSums(Y)
m <- nrow(Y)
N.mat <- matrix(rep(N, m), nrow = m, byrow = TRUE)
N.mat[Y > 0] <- 0
tmp <- N[max.col(N.mat)]
Y <- Y + N.mat / tmp
logY <- log2(Y)
W <- t(t(logY) - colMeans(logY))
sub <- W[c("p__Bacteroidota;g__Bacteroides","p__Firmicutes_A;g__Blautia","p__Firmicutes_A;g__Faecalimonas"),]
rownames(sub) <- gsub('.*g__','',rownames(sub))

yy2 <- yy %>% mutate(grp = ifelse(icd10_first_3_name=='Control','Control','Cancer')) %>% 
  dplyr::filter(icd10_first_3_name %in% c(names(cancer.dir),'Control'))
yy2 <- merge(yy2 %>% column_to_rownames("Row.names"), t(sub),by = 0)



#orange       skyblue   bluishgreen        yellow          blue    vermillion reddishpurple 
#"#E69F00"     "#56B4E9"     "#009E73"     "#F0E442"     "#0072B2"     "#D55E00"     "#CC79A7" 


#sel_colors_plot <- RColorBrewer::brewer.pal(12,"Paired")
sel_colors_plot <- c("#CC79A7", "#0072B2") #reddishpurple and blue



#eventually update these arrows using the coabundance networks


species_data <- yy2[, c("Bacteroides", "Blautia", "Faecalimonas")]
pc_data <- yy2[, c("PC1", "PC2")]
arrows_data <- data.frame(
  species = c("Bacteroides", "Blautia", "Faecalimonas"),
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
  theme(text = element_text(size = 20, color = "black"),
        axis.text = element_text(size = 15, color = "black"),
        legend.text = element_text(size = 20, color = "black"),
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
  species = c("Bacteroides", "Blautia", "Faecalimonas"),
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
        axis.text = element_text(size = 20, color = "black"),
        legend.text = element_text(size = 15, color = "black"),
        axis.title = element_text(size = 20, color = "black"),
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


p <- ggarrange(p2, p3, nrow =1, common.legend = F)
print(p)
#ggsave('PCoA_v2.pdf', width = 16, height = 8)
ggsave('../Figure/RM_figures/PCoA_v2.png', width = 16, height = 7)



#-------------------------------------------------------------------------------
## ==== Control vs subCancerX Alpha plot ====
setwd(rd)
setwd('subCancerX_Control')
tm = load('Control-Lymphoid_leukemia/Alpha/icd10_first_3_name/Alpha.RData')
dirs <- list.dirs(full.names = F, recursive = F)
alpha <- NULL
for(dir in dirs){
  setwd(rd)
  setwd('subCancerX_Control')
  setwd(dir)
  load('Alpha/icd10_first_3_name/Alpha.RData')
  tmp.alpha <- merge(alpha.obj1, data.obj.rff2$meta.dat[,'icd10_first_3_name', drop =F], by = 0) %>% 
    dplyr::filter(icd10_first_3_name!='Control')
  alpha <- rbind(alpha, tmp.alpha)
}

load('Alpha/icd10_first_3_name/Alpha.RData')
tmp.alpha <- merge(alpha.obj1, data.obj.rff2$meta.dat[,'icd10_first_3_name', drop =F], by = 0)%>%
  dplyr::filter(icd10_first_3_name=='Control')
alpha <- rbind(alpha, tmp.alpha)

col <- c(brewer.pal(12,'Paired'),brewer.pal(8,'Dark2'), brewer.pal(8,'Accent')[1:4])
names(col) <- c('Control', cancer_class_names)

setwd(fd)
load('subCancerX_Control/Alpha_P_R2.RData')
pval <- t(pval_coef_adj.All['Shannon',,]['P',,drop =F])
rownames(pval) <- gsub('Control-','',rownames(pval))
pval <- as.data.frame(pval)  %>% mutate(sig = ifelse(P<0.05, '*',''))

#Shannon
alpha$icd10_first_3_name <- as.character(tolower(gsub("Malignant neoplasm of |Malignant |malignant |Other and | neoplasm| neoplasm of| neoplasms|neoplasm of |\\,|and ", "",  alpha$icd10_first_3_name)))
alpha$icd10_first_3_name <- as.factor(alpha$icd10_first_3_name)
rownames(pval) <- gsub("_", " ", rownames(pval))
rownames(pval) <- as.character(tolower(gsub("Malignant neoplasm of |Malignant |malignant |Other and | neoplasm| neoplasm of| neoplasms|neoplasm of |\\,|and ", "",  rownames(pval))))
rownames(pval)[rownames(pval) == "plasma cell"] <- "multiple myeloma plasma cell"


ord <- as.vector((aggregate(Shannon ~ icd10_first_3_name, data=alpha, function(x) median(x)) %>% arrange(Shannon))[,1])
alpha <- within(alpha, icd10_first_3_name <- factor(icd10_first_3_name,levels =ord))
Means <- alpha %>% group_by(icd10_first_3_name) %>% summarize(Avg = median(Shannon))

Means$icd10_first_3_name
length(unique(alpha$icd10_first_3_name)) #24


#-----------------------------------
#version for the main manuscript

#color scale
col_fun <- colorRamp2(c(max(Means$Avg), min(Means$Avg)), c("#0072B2", "#CC79A7"))
col <- col_fun(Means$Avg)
names(col) <- as.character(Means$icd10_first_3_name)


p1 <- ggplot() + 
  geom_violin(data = alpha, aes(y = icd10_first_3_name, x = Shannon, fill = icd10_first_3_name), trim = F, color = NA) + 
  # geom_boxplot(data = alpha, mapping = aes(y = icd10_first_3_name, x = Shannon, group = icd10_first_3_name), outlier.size = 0.5, width = 0.1) +
  # geom_point(data = Means, mapping = aes(y = icd10_first_3_name, x = Avg), size = 0.5) +
  geom_line(data = Means, mapping = aes(y = icd10_first_3_name, x = Avg, group = 1)) +
  scale_fill_manual(values = col) + 
  theme_classic() + 
  theme(text = element_text(size = 16, color = "black"),
        axis.text.y = element_text(size = 16, color = "black", angle = 0), # Adjust if needed
        axis.text.x = element_text(size = 16, color = "black"),
        axis.title.x = element_text(size = 16, color = "black"),
        axis.title.y = element_text(size = 16, color = "black"),
        panel.background = element_rect(fill = "transparent", colour = NA),
        plot.background = element_rect(fill = "transparent", colour = NA),
        legend.position = 'none') + 
  labs(x = 'Shannon diversity index', y = '') + # Adjust labels accordingly
  geom_text(data = pval %>% rownames_to_column('icd10_first_3_name') %>% mutate(pos = max(alpha$Shannon) * 1.12), 
            aes(y = icd10_first_3_name, x = pos, label = sig))

setwd("~/Library/Mobile Documents/com~apple~CloudDocs/mforge_clean/Figure/RM_figures/")
ggsave(file = paste0('alpha cancer class vs control violin main text.png'), plot=p1, width = 7, height =7, bg = "transparent")



#-------------------------------------------------------------------------------
## ==== Figure X ======
vars.incl <- covars[!(covars %in% c("Batch","Cancer_class"))]
blood.incl <- c('Erythrocytes','Hematocrit','Neutrophils','MCV','Hemoglobin','Leukocytes','Platelet.Count','Neutrophils_cat','Platelet.Count_cat','Hemoglobin_cat')
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
df2 <- df %>% dplyr::filter(var %in% vars.incl) %>% mutate(var = gsub('icd10_first_3_name','subCancer_class',var))
unique(df2$var)
df2$var[df2$var =="subCancer_class"] <- "Cancer class"
df2$var[df2$var =="Bristol_score"] <- "Bristol stool score"
df2$var[df2$var =="PPI_day_365"] <- "PPI in past year (No)"
df2$var[df2$var =="PPI_last_month"] <- "PPI in past month (No)"
df2$var[df2$var =="Abx_day_365"] <- "Antibiotics in past year (No)"
df2$var[df2$var =="Abx_last_month"] <- "Antibiotics in past month (No)"
df2$var[df2$var =="Elix_score"] <- "Elixhauser Comorbidity score"
df2$var[df2$var =="Sample_season"] <- "Sample season"
df2$var[df2$var =="Urban"] <- "Residence type (non-urban)"
df2$var[df2$var =="subCancer_class"] <- "Cancer type"
df2$var[df2$var =="Sex"] <- "Sex (Female)"
df2$var[df2$var =="Metastasis"] <- "Metastasis (No)"


#plotting settings start here
#version for main figure

#only get the alpha and Beta div adjusted. Color the variables by their class type
#technical (Bristol Score, Sample Season, Antibiotics 2x, PPI), demographic (Age, Sex, Residence type), clinical (all others)

df3 <- df2
df3 <- df3[!df3$diversity == "Beta diversity (pathway)",]
df3 <- df3[df3$grp == "adjusted",]

#unique(df3$var)
variables_to_plot <- c("Cancer class" , "Bristol stool score", "Elixhauser Comorbidity score","Sex (Female)", "Age",
  "Antibiotics in past year (No)", "BMI", "PPI in past year (No)", "Residence type (non-urban)", "Sample season",
  "PPI in past month (No)", "Metastasis (No)", "Antibiotics in past month (No)") 

variable_colors_to_plot <- rep(sel_colors[6], length(variables_to_plot))
names(variable_colors_to_plot) <- variables_to_plot

variable_colors_to_plot[variables_to_plot %in% c("Bristol stool score", "Sample season")] <- sel_colors[2]
variable_colors_to_plot[variables_to_plot %in% c("Sex (Female)", "Age", "Residence type (non-urban)", "BMI")] <- sel_colors[1]


main_div_plot <- ggplot(df3,aes(x =reorder(var,value), y= value)) + 
  geom_segment(aes(x=reorder(var,value), xend=reorder(var,-value), y=0, yend=value), size =4, color = "darkgrey") +
  geom_point(size=5, aes(shape=direction, fill = P.col), stroke=0.5) + 
  scale_fill_manual(values = c('sig'='black','no'= 'white')) +
  scale_shape_manual(values = c('positive'=24,'negative'=25, 'none'= 21)) +
  facet_grid(grp~diversity, scales = 'free')  + 
  theme_classic() +
  labs(y = ('Variability explained (R²), %'), fill = '', x = '', shape ='') +
  theme(panel.grid = element_blank(),
        axis.title = element_text(size=16, color = 'black'),
        legend.text = element_text(size=16, color = 'black'),
        strip.text = element_text(size=16, color = 'black'),
        axis.text.x = element_text(size=16,color = 'black'),
        axis.text.y = element_text(size=16, color = rev(variable_colors_to_plot))) +
  theme(panel.spacing = unit(2, "lines")) +

  labs(subtitle = "", caption = c("Clinical variables\nDemographic variables\nTechnical variables")) +
  theme(plot.caption = element_text(colour = sel_colors[c(2)], vjust=0.4, hjust = 0, angle = 0, size=16),
        plot.subtitle = element_text(size = 1, face = "bold", hjust = 0)) +
  
    coord_flip() + guides(fill = FALSE) + ggtitle('Figure1') 

main_div_plot

dev.off()

#save as svg to modify the colors of the legend in illustrator
setwd("~/Library/Mobile Documents/com~apple~CloudDocs/mforge_clean/Figure/RM_figures/")
ggsave(file = 'Figure1_main.svg', plot=main_div_plot, width = 12, height = 4.75)



sel_colors
#orange       skyblue   bluishgreen        yellow          blue    vermillion reddishpurple 
#"#E69F00"     "#56B4E9"     "#009E73"     "#F0E442"     "#0072B2"     "#D55E00"     "#CC79A7" 


#---------------------------------
#original version
ggplot(df2,aes(x =reorder(var,value), y= value)) + 
geom_segment( aes(x=reorder(var,value), xend=reorder(var,-value), y=0, yend=value), size =1, color = brewer.pal(8, 'Dark2')[8]) +
geom_point(size=4, aes(shape=direction, fill = P.col), stroke=0.5) + 
scale_fill_manual(values = c('sig'='#f1abb9','no'= 'white')) +
scale_shape_manual(values = c('positive'=24,'negative'=25, 'none'= 21)) +
facet_grid(grp~diversity, scales = 'free')  + 
theme_classic() +
labs(y = ('Variability explained (R²), %'), fill = '', x = '', shape ='') +
theme(panel.grid = element_blank(),
      axis.title = element_text(size=16, color = 'black'),
      legend.text = element_text(size=16, color = 'black'),
      strip.text = element_text(size=16, color = 'black'),
      axis.text.x = element_text(size=16,color = 'black'),
      axis.text.y = element_text(size=16, color = 'black')) + 

# guides(fill = guide_legend(override.aes = list(fill = c('#f1abb9','white'))))  +
coord_flip() + guides(fill = FALSE) + ggtitle('Figure1')

setwd(mfd)
ggsave(file = paste0('Figure1_original.pdf'), width = 18,height = 7)

summary(df2[df2$grp=='adjusted','value'])
summary(df2[df2$grp=='marginal','value'])

sum(df2[df2$grp=='adjusted'&df2$diversity=='Alpha diversity','value'])
sum(df2[df2$grp=='adjusted'&df2$diversity=='Beta diversity(pathway)','value'])
sum(df2[df2$grp=='adjusted'&df2$diversity=='Beta diversity(species)','value'])

sum(df2[df2$grp=='marginal'&df2$diversity=='Alpha diversity','value'])
sum(df2[df2$grp=='marginal'&df2$diversity=='Beta diversity(pathway)','value'])
sum(df2[df2$grp=='marginal'&df2$diversity=='Beta diversity(species)','value'])

df2[df2$grp=='marginal'&df2$diversity=='Alpha diversity',] %>% arrange(-value) %>% top_n(3,value)
df2[df2$grp=='marginal'&df2$diversity=='Beta diversity(species)',] %>% arrange(-value) %>% top_n(3,value)
df2[df2$grp=='marginal'&df2$diversity=='Beta diversity(pathway)',] %>% arrange(-value) %>% top_n(3,value)

df2[df2$grp=='adjusted'&df2$diversity=='Alpha diversity',] %>% arrange(-value) %>% top_n(3,value)
df2[df2$grp=='adjusted'&df2$diversity=='Beta diversity(species)',] %>% arrange(-value) %>% top_n(3,value)
df2[df2$grp=='adjusted'&df2$diversity=='Beta diversity(pathway)',] %>% arrange(-value) %>% top_n(3,value)

blood.a <- ggplot(df %>% dplyr::filter(var %in% blood.incl[-grep('_cat$',blood.incl)]),aes(x =reorder(var,value), y= value)) + 
  geom_segment( aes(x=reorder(var,value), xend=reorder(var,-value), y=0, yend=value), size =1, color = brewer.pal(8, 'Dark2')[8]) +
  geom_point(size=4, aes(shape=direction, fill = P.col), stroke=0.5) + 
  scale_fill_manual(values = c('sig'='#f1abb9','no'= 'white')) +
  scale_shape_manual(values = c('positive'=24,'negative'=25, 'none'= 21)) +
  facet_grid(grp~diversity, scales = 'free')  + 
  theme_classic() +
  labs(y = ('Variability explained (R²), %'), fill = '', x = '', shape ='') +
  theme(panel.grid = element_blank(),
        axis.title = element_text(size=16, color = 'black'),
        legend.text = element_text(size=16, color = 'black'),
        strip.text = element_text(size=16, color = 'black'),
        axis.text.x = element_text(size=16,color = 'black'),
        axis.text.y = element_text(size=16, color = 'black')) + 
  # guides(fill = guide_legend(override.aes = list(fill = c('#f1abb9','white'))))  +
  coord_flip() + guides(fill = FALSE) 
p3 <- ggplot(df %>% dplyr::filter(var %in% c("neut_neutropenia2","bone_marrow_suppression2","Hb_anemia2","Pl_thrombocytopenia2")) %>% dplyr::mutate(var = gsub('2$','',var)),aes(x =reorder(var,value), y= value)) + 
  geom_segment( aes(x=reorder(var,value), xend=reorder(var,-value), y=0, yend=value), size =1, color = brewer.pal(8, 'Dark2')[8]) +
  geom_point(size=4, aes(shape=direction, fill = P.col), stroke=0.5) + 
  scale_fill_manual(values = c('sig'='#f1abb9','no'= 'white')) +
  scale_shape_manual(values = c('positive'=24,'negative'=25, 'none'= 21)) +
  facet_grid(grp~diversity, scales = 'free')  + 
  theme_classic() +
  labs(y = ('Variability explained (R²), %'), fill = '', x = '', shape ='') +
  theme(panel.grid = element_blank(),
        axis.title = element_text(size=16, color = 'black'),
        legend.text = element_text(size=16, color = 'black'),
        strip.text = element_text(size=16, color = 'black'),
        axis.text.x = element_text(size=16,color = 'black'),
        axis.text.y = element_text(size=16, color = 'black')) + 
  # guides(fill = guide_legend(override.aes = list(fill = c('#f1abb9','white'))))  +
  coord_flip() + guides(fill = FALSE) 
blood.a
# setwd(mfd)
# ggsave(file = paste0('FigureSN-7A.pdf'), width = 12,height = 5)
# p3
# setwd(mfd)
# ggsave(file = paste0('Figure1_blood_outcome.pdf'), width = 18,height = 4)



#-------------------------------------------------------------------------------
## ===== Elix components beta diversity ======
x1 <- gsub('Elixhauser_','',rev((df[df$var %in% elix.names & df$grp =='adjusted' & df$diversity =='Beta diversity(taxa)' & df$P.col =='sig',] %>% arrange(value))[,1]))
x2 <- (gsub('Elixhauser_','',rev((df[df$var %in% elix.names & df$grp =='adjusted' & df$diversity =='Beta diversity(pathway)' & df$P.col =='sig',] %>% arrange(value))[,1])))
intersect(x1,x2)
paste0(x2, collapse = ', ')
df.tmp <- df[df$var %in% elix.names & df$grp =='adjusted' & df$diversity != 'Alpha diversity',]
head(df.tmp)
df.tmp$var <- names_map[as.vector(df.tmp$var)]

df.tmp %>% 
  ggplot(aes(x =reorder(var,value), y= value)) + 
  geom_segment(aes(x=reorder(var,value), xend=reorder(var,-value), y=0, yend=value), size =1, color = brewer.pal(8, 'Dark2')[8]) +
  geom_point(size=4, shape=21, aes(fill = P.col), stroke=0.5) + 
  scale_fill_manual(values = c('sig'='#f1abb9','no'= 'white')) +
  scale_shape_manual(values = c('positive'=24,'negative'=25, 'none'= 21)) +
  facet_grid(~diversity, scales = 'free')  + 
  theme_classic() +
  labs(y = ('Variability explained (R²), %'), fill = '', x = '', shape ='') +
  theme(panel.grid = element_blank(),
        axis.title = element_text(size=16, color = 'black'),
        legend.text = element_text(size=16, color = 'black'),
        strip.text = element_text(size=16, color = 'black'),
        axis.text.x = element_text(size=16,color = 'black'),
        axis.text.y = element_text(size=16, color = 'black')) + 
  coord_flip() + guides(fill = FALSE) + ggtitle('Figure2A')


setwd(mfd)
#ggsave(file = paste0('Figure2A.pdf'), width = 14,height = 8)
ggsave(file = paste0('Figure2A_v3.pdf'), width = 14,height = 7)


#-----------------------
#remake the plot above only for the species beta diversity

#remove Lymphoma
df.tmp <- df.tmp[!df.tmp$var == "Lymphoma",]
df.tmp <- df.tmp[!df.tmp$var == "Solid tumor without metastasis",]

df.tmp_2 <- df.tmp
df.tmp_2 <- df.tmp_2[df.tmp_2$diversity == "Beta diversity (species)",]


#get the n numbers for these categories
df.tmp_2$var
#[1] "Congestive heart failure"                        "Cardiac arrhythmias"                             "Valvular disease"                               
#[4] "Pulmonary hypertension"                          "Peripheral vascular disorders"                   "Hypertension (uncomplicated)"                   
#[7] "Other neurologic disorders"                      "Chronic pulmonary disease"                       "Diabetes without chronic complications"         
#[10] "Diabetes with chronic complications"             "Hypothyroidism"                                  "Renal failure"                                  
#[13] "Liver disease"                                   "Peptic ulcer disease (excluding bleeding)"       "Metastatic cancer"                              
#[16] "Rheumatoid arthritis/collagen vascular diseases" "Coagulopathy"                                    "Obesity"                                        
#[19] "Weight loss"                                     "Fluid and electrolyte disorders"                 "Deficiency anemia"                              
#[22] "Alcohol abuse"                                   "Drug abuse"                                      "Depression"  

table(data.obj$meta.dat$Elixhauser_CHF) # 68
table(data.obj$meta.dat$Elixhauser_Arrhythmia) # 344
table(data.obj$meta.dat$Elixhauser_Valvular) # 180
table(data.obj$meta.dat$Elixhauser_PHTN) # 73
table(data.obj$meta.dat$Elixhauser_PVD) # 207
table(data.obj$meta.dat$Elixhauser_HTN) # 453
table(data.obj$meta.dat$Elixhauser_NeuroOther) # 84
table(data.obj$meta.dat$Elixhauser_Pulmonary) #173
table(data.obj$meta.dat$Elixhauser_DM) # 72
table(data.obj$meta.dat$Elixhauser_DMcx) # 95
table(data.obj$meta.dat$Elixhauser_Hypothyroid) #8, 192
table(data.obj$meta.dat$Elixhauser_Renal) #9, 137
table(data.obj$meta.dat$Elixhauser_Liver) #10, 241
table(data.obj$meta.dat$Elixhauser_PUD) #11, 23
table(data.obj$meta.dat$Elixhauser_Mets) #12, 772
table(data.obj$meta.dat$Elixhauser_Rheumatic) #13, 54
table(data.obj$meta.dat$Elixhauser_Coagulopathy) #14, 91
table(data.obj$meta.dat$Elixhauser_Obesity) #15, 292
table(data.obj$meta.dat$Elixhauser_WeightLoss) #16, 119
table(data.obj$meta.dat$Elixhauser_FluidsLytes) #17, 229
table(data.obj$meta.dat$Elixhauser_Anemia) #18, 69
table(data.obj$meta.dat$Elixhauser_Alcohol) #19, 30
table(data.obj$meta.dat$Elixhauser_Drugs) #20, 24
table(data.obj$meta.dat$Elixhauser_Depression) #21, 159

n_numbers <- c(68, 344, 180, 73, 207, 453, 84, 173, 72, 95, 192, 137, 241, 23, 772, 54, 91, 292, 119, 229, 69, 30, 24, 159)
df.tmp_2$var <- paste0(df.tmp_2$var, "  n=", n_numbers)

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
        axis.title = element_text(size=16, color = 'black'),
        legend.text = element_text(size=16, color = 'black'),
        strip.text = element_text(size=16, color = 'black'),
        axis.text.x = element_text(size=16,color = 'black'),
        axis.text.y = element_text(size=16, color = 'black')) + 
  coord_flip() + guides(fill = FALSE) + ggtitle('Beta diversity (species)') +
  theme(plot.margin = margin(t = 20, r = 30, b = 10, l = 10)) 

elix_div_plot


setwd("~/Library/Mobile Documents/com~apple~CloudDocs/mforge_clean/Figure/RM_figures/")
ggsave(file = paste0('main_elix_div_plot.png'), plot=elix_div_plot, width = 9, height = 6, bg = "transparent")



#-------------------------------------------------------------------------------
## === Elix components tree =====

#adjusted the cutoffs to make it more legible
setwd(wd)
load('Data/data.obj.raw.core.RData')
unique.f = F
cutoff = 0.05
setwd(fd)
setwd('CancerOnly/')
tree.tmp0 <- read.csv('Species_DAA.csv') %>% dplyr::inner_join(as.data.frame(data.obj$otu.name)%>% rownames_to_column('species')) %>% dplyr::filter(variable %in% elix.names)

xx <- read.csv('Species_DAA.csv') %>% dplyr::inner_join(as.data.frame(data.obj$otu.name)%>% rownames_to_column('species')) %>% dplyr::filter(variable =='Elix_score')
tree.tmp0 <- rbind(tree.tmp0, xx)
tree.tmp0$variable <- names_map[tree.tmp0$variable]

tree.tmp <- tree.tmp0 %>% dplyr::select(species, Q, variable) %>% spread(key = 'variable', value = 'Q') %>% column_to_rownames('species')
tree.tmp <- tree.tmp[apply(tree.tmp, 1, function(x) sum(x<cutoff, na.rm = T)>0),]
tree.tmp[is.na(tree.tmp)] <- 1
idx <- names(which(apply(tree.tmp, 2, function(x) sum(x<cutoff, na.rm = T)) >0))
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

## Only keep 50% taxa with Q<cutoff
keep.tax <- c()
for(i in colnames(tree.tmp)){
  topQ <- which(tree.tmp[,i]<cutoff) 
  keepR2 <- rownames(R2[topQ,i, drop =F] %>% arrange(-abs(!!as.name(i))))[1:round(length(topQ)*0.5)]
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

tree.sub <- groupOTU(tree.sub, taxa.split, group_name = 'taxa.split')
rownames(tree.tmp) <- gsub('s\\_\\_','',rownames(tree.tmp))
tree.tmp <- tree.tmp[tree.sub$tip.label,]

dd <- tree.tmp0 %>% dplyr::select(species, Q, variable) %>% mutate(sig = ifelse(Q<cutoff, 'yes','no')) %>% dplyr::select(-Q)
dd$species <- gsub('s\\_\\_','',dd$species)
dd <- within(dd, species <- factor(species, levels = tree.sub$tip.label))

ord <- gsub('f\\_\\_|p\\_\\_|o\\_\\_|c\\_\\_|g\\_\\_','',unique(otu.name[get_taxa_name(ggtree(tree.sub)),branch.level]))
gc()
plt1 <- ggtree(tree.sub, layout = 'circular') + 
  geom_tiplab(align=T,aes(color =taxa.split), fontface='italic',size =5) +
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

col <- c(brewer.pal(8,'Reds')[c(8,5)],brewer.pal(8,'Purples')[c(4,8)],brewer.pal(8,'Set1')[c(8)],brewer.pal(8,'Set2')[c(6)],
         brewer.pal(8,'Greens')[c(8,6)],brewer.pal(8,'Dark2')[c(2:6)],brewer.pal(8,'Blues')[c(8,6)])[1:length(c(pos,neg))]
col1 <- c(col,'white')
names(col1) <- c(xx2,'No')

new_order <- c(1:2,4,3,5:13)
xx2 <- xx2[new_order]
breaks_to_show <- xx2

gc()
gheatmap(plt1, tree.tmp, offset=3, width=0.8, font.size=4, color = brewer.pal(9,'Set1')[9], colnames = F) +
  scale_fill_manual(values=col1, name="",breaks=breaks_to_show, 
                    guide = guide_legend(label.theme = element_text(angle = 0, face = 'plain'))) 


setwd(mfd)
if(cutoff ==0.05){width = height = 24}
if(cutoff ==0.01){width = height = 14}
#ggsave(paste0('Figure2B_v2.pdf'), width = 19, height = 18)
ggsave(paste0('../Figure/RM_figures/elix_components circle plot_RM.png'), width = width, height = height)
gc()





#-------------------------------------------------------------------------------
## ==== FigureX: subCancerX vs all other cancers alpha & beta summary plot======
setwd(wd)
load('Figure/subCancerX-Ex/Alpha_P_R2.RData') ## change to subCancerX-Ex, 2014/11/01 for adding 2 other cancers back
alpha.measure = 'Shannon'
pval_coef_adj <- t(pval_coef_adj.All[alpha.measure,,])
pval_coef_adj <- pval_coef_adj#[!(rownames(pval_coef_adj) %in% c("other_connective_soft_tissue","Other_unspecified_malignant_neoplasm_of_skin")),]
pval_coef_adj[pval_coef_adj[,'R2']<0,'R2'] <- 0
pval_coef_adj <- as.data.frame(pval_coef_adj) %>% 
  mutate(direction = ifelse(is.na(coef),'none',ifelse(coef>0, 'positive','negative'))) %>%
  dplyr::select(-coef) %>%
  tibble::rownames_to_column('var') %>% melt()
sub.r <- pval_coef_adj[pval_coef_adj$variable=='R2',]
sub.p <- pval_coef_adj[pval_coef_adj$variable=='P',c('var','value')] %>% 
  dplyr::mutate(P.col = ifelse(value<0.05,'sig','no')) %>% dplyr::select(-value)
sub.rp <- dplyr::inner_join(sub.p, sub.r) %>% mutate(value = value *100)
ord <- (sub.rp[order(sub.rp$value),'var'])
sub.rp <- within(sub.rp, var <- factor(var, levels=ord))

dist.name <- 'WUniFrac'
load('Figure/subCancerX-Ex/Beta_P_R2.RData')
R2_P <- cbind.data.frame(value = c(P.adj[dist.name,]), R2 = c(R2.adj[dist.name,])) %>% 
  rownames_to_column('var') %>% #dplyr::filter(!(var %in% c("other_connective_soft_tissue","Other_unspecified_malignant_neoplasm_of_skin"))) %>% 
  mutate(P.col = ifelse(value<0.05, 'sig','no'))  %>% dplyr::select(-value) %>% mutate(direction = 'none', variable = 'R2') %>% dplyr::rename(value = R2) %>%
  mutate(value=ifelse(value <0,0, value *100))
ord <- (R2_P[order(R2_P$value),'var'])
R2_P <- within(R2_P, var <- factor(var, levels=ord))


setwd(fd)
dist.name <- 'BC'
load('subCancerX-Ex_func/pathway/Beta_P_R2.RData')
dist.names <- rownames(R2.adj)
R2_P2 <- cbind.data.frame(value = P.adj[dist.name,],R2 = R2.adj[dist.name,]) %>% 
  rownames_to_column('var') %>% dplyr::filter(var !='Batch') %>% 
  mutate(P.col = ifelse(value<0.05, 'sig','no'))  %>% dplyr::select(-value) %>% mutate(direction = 'none', variable = 'R2') %>% dplyr::rename(value = R2) %>%
  mutate(value=ifelse(value <0,0, value *100)) 
ord <- (R2_P2[order(R2_P2$value),'var'])
R2_P2 <- within(R2_P2, var <- factor(var, levels=ord))

R2_P_all <- rbind(R2_P %>% mutate(diversity = 'Beta diversity(species)'), 
                  sub.rp %>% mutate(diversity = 'Alpha diversity'),
                  R2_P2 %>% mutate(diversity = 'Beta diversity(pathway)'))
R2_P_all <- within(R2_P_all, diversity <- factor(diversity, levels=c('Alpha diversity','Beta diversity(species)','Beta diversity(pathway)')))

load('../Data/cancerNames_20141101.Rdata')
# load('../Data/cancerNames.Rdata')
# names(cancer.dir) <- cancer.type
names(cancer.dir) <- tolower(gsub("Malignant neoplasm of |Malignant |malignant |Other and | neoplasm| neoplasms", "", cancer.type))
# names(cancer.dir) <- tolower(gsub("Malignant neoplasm of |Malignant |malignant |Other and ", "", cancer.type))
indices <- match(as.vector(R2_P_all$var), cancer.dir)
R2_P_all$var <- names(cancer.dir[indices])
R2_P_all1 <- R2_P_all[R2_P_all$diversity=='Alpha diversity',]  %>% arrange(value)
R2_P_all2 <- R2_P_all[R2_P_all$diversity=='Beta diversity(species)',] #%>% droplevels()
R2_P_all3 <- R2_P_all[R2_P_all$diversity=='Beta diversity(pathway)',] #%>% droplevels()
R2_P_all1$var <- factor(R2_P_all1$var, levels = R2_P_all1$var)
R2_P_all2 <- within(R2_P_all2, var <- factor(var, levels=R2_P_all1$var))
R2_P_all3 <- within(R2_P_all3, var <- factor(var, levels=R2_P_all1$var))

p1 <- ggplot(R2_P_all1,aes(x = var, y= value)) + 
  geom_segment( aes(x=var, xend= var, y=0, yend=value), size =1, color = brewer.pal(8, 'Dark2')[8]) +
  geom_point(size=4, aes(shape=direction, fill = P.col), stroke=0.5) + 
  scale_fill_manual(values = c('sig'='#f1abb9','no'= 'white')) +
  scale_shape_manual(values = c('positive'=24,'negative'=25, 'none'= 21), breaks = c('positive','none','negative'), drop = FALSE) +
  facet_wrap(.~diversity, scales = 'free')  +
  theme_classic() +
  labs(y = (''), fill = '', x = '', shape ='') +
  theme(panel.grid = element_blank(),
        axis.title.x = element_text(size=16, color = 'black'),
        legend.text = element_text(size=16, color = 'black'),
        strip.text = element_text(size=16, color = 'black'),
        axis.text.x = element_text(size=16,color = 'black'),
        axis.text.y = element_text(size=16,color = 'black')) + 
  # guides(fill = guide_legend(override.aes = list(fill = c('#f1abb9','white'))))  +
  coord_flip() + guides(fill = FALSE) 
p1
p2 <- ggplot(R2_P_all2,aes(x = var, y= value)) + 
  geom_segment( aes(x=var, xend=var, y=0, yend=value), size =1, color = brewer.pal(8, 'Dark2')[8]) +
  geom_point(size=4, aes(shape=direction, fill = P.col), stroke=0.5) + 
  scale_fill_manual(values = c('sig'='#f1abb9','no'= 'white')) +
  scale_shape_manual(values = c('positive'=24,'negative'=25, 'none'= 21), breaks = c('positive','none','negative'),) +
  facet_wrap(.~diversity, scales = 'free')  +
  theme_classic() +
  labs(y = ('Variability explained (R²), %'), fill = '', x = '', shape ='') +
  theme(panel.grid = element_blank(),
        axis.title.x = element_text(size=16, color = 'black'),
        legend.text = element_text(size=16, color = 'black'),
        strip.text = element_text(size=16, color = 'black'),
        axis.text.y = element_blank(),
        axis.line.y =element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.x = element_text(size=16, color = 'black')) + 
  # guides(fill = guide_legend(override.aes = list(fill = c('#f1abb9','white'))))  +
  coord_flip() + guides(fill = FALSE) 
p2
p3 <- ggplot(R2_P_all3,aes(x = var, y= value)) + 
  geom_segment( aes(x=var, xend=var, y=0, yend=value), size =1, color = brewer.pal(8, 'Dark2')[8]) +
  geom_point(size=4, aes(shape=direction, fill = P.col), stroke=0.5) + 
  scale_fill_manual(values = c('sig'='#f1abb9','no'= 'white')) +
  scale_shape_manual(values = c('positive'=24,'negative'=25, 'none'= 21), breaks = c('positive','none','negative'),drop = FALSE) +
  facet_wrap(.~diversity, scales = 'free')  +
  theme_classic() +
  labs(y = (''), fill = '', x = '', shape ='') +
  theme(panel.grid = element_blank(),
        axis.title.x = element_text(size=16, color = 'black'),
        legend.text = element_text(size=16, color = 'black'),
        strip.text = element_text(size=16, color = 'black'),
        axis.text.y = element_blank(),
        axis.line.y =element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.x = element_text(size=16, color = 'black')) + 
  # guides(fill = guide_legend(override.aes = list(fill = c('#f1abb9','white'))))  +
  coord_flip() + guides(fill = FALSE) 
pp <- ggarrange(p1,p2, p3, common.legend = T, nrow = 1, widths = c(3.5,2.2,2.2), legend = 'right')
annotate_figure(pp, top = text_grob("Figure3B",face = "bold", size = 14))
setwd(mfd)
#ggsave(file = paste0('Figure3B.pdf'), width = 18,height = 6)
# ggsave(file = paste0('Figure3B_v2.pdf'), width = 18,height = 6)
ggsave(file = paste0('Figure3B_v3.pdf'), width = 18,height = 6)


#-------------------------
#remake this one with alpha and beta only

#rbind these 2 together
df3 <- rbind(R2_P_all1, R2_P_all2)
#rename the cancer classes
df3$var <- tolower(gsub("Malignant neoplasm of |Malignant |malignant |Other and | neoplasm| neoplasm of| neoplasms|neoplasm of |\\,|and ", "",  df3$var))

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


dev.off()
#save as svg to modify the colors of the legend in illustrator
setwd("~/Library/Mobile Documents/com~apple~CloudDocs/mforge_clean/Figure/RM_figures/")
ggsave(file = paste0('Figure3B subCancerX vs other cancers.png'), plot=cancerclass_div_plot, width = 12, height = 6)




#-------------------------------------------------------------------------------
#Blood outcomes figure

## ====== Hemoglobin&Platelet.Count func =====
setwd(wd)
load('Figure/CancerOnly_func/pathway/DAA_P_R2.RData')
load('Result/CancerOnly_func/pathway/data.obj.wk.RData')
level = 'pathway'
sp <- Q.All[[level]][,c(level,'Erythrocytes','Hematocrit','Neutrophils','MCV','Hemoglobin','Leukocytes','Platelet.Count')] %>% column_to_rownames(level)
prop <- t(t(data.obj$abund.list[[level]])/colSums(data.obj$abund.list[[level]]))
xx = apply(sp, 2, function(x) sum(x<0.05, na.rm = T))
xx = names(xx[xx>0])
xx
tt <- (R2.All[[level]] %>% column_to_rownames(level))

m1 <- merge(sp[which(sp[,"Hematocrit"]<0.05),"Hematocrit", drop =F]%>% dplyr::rename(Q=1), tt[which(sp[,"Hematocrit"]<0.05),"Hematocrit",drop =F] %>% dplyr::rename(R2=1), by = 0) %>% arrange(R2)
m2 <- merge(sp[which(sp[,"Hemoglobin"]<0.05),"Hemoglobin", drop =F]%>% dplyr::rename(Q=1), tt[which(sp[,"Hemoglobin"]<0.05),"Hemoglobin",drop =F] %>% dplyr::rename(R2=1), by = 0) %>% arrange(R2) 
blood.c <- rbind(m1 %>% mutate(group = "H"), m2 %>% mutate(group = "Hemoglobin")) %>% mutate(direction = ifelse(R2>0,'positive','negative')) %>% 
  ggplot(aes(x = reorder(Row.names, R2), y = R2, fill = direction)) + geom_bar(stat = 'identity') + #coord_flip() + 
  theme_bw() + 
  facet_grid(~group, scales = 'free', space = 'free') + 
  scale_fill_manual(values = c('positive'= '#3880a6', 'negative'='#f3c83c')) + 
  theme(axis.title = element_text(color = 'black', size = 16),
        axis.text.x = element_text(angle = 90, hjust = 1,vjust = 0.2, size = 16,color = 'black' ),
        axis.text.y = element_text(color = 'black', size = 16),
        legend.position = 'top',
        strip.text = element_text(color = 'black', size = 16),
        legend.text = element_text(color = 'black', size = 16),
        legend.title = element_text(color = 'black', size = 16),
        panel.grid = element_blank()) +
  labs(y = 'Effect size', x = '', fill = 'Association direction') 

blood.a
blood.c



## ====== Blood outcomes =====
setwd(wd)
load('Figure/CancerOnly/DAA_P_R2.RData')
load('Result/CancerOnly/data.obj.wk.RData')
level = 'Genus'
sp <- Q.All[[level]][,c(level,"neut_neutropenia2","bone_marrow_suppression2","Hb_anemia2","Pl_thrombocytopenia2")] %>% column_to_rownames(level)
prop <- t(t(data.obj$abund.list[[level]])/colSums(data.obj$abund.list[[level]]))
apply(sp, 2, function(x) sum(x<0.05, na.rm = T))
var = "bone_marrow_suppression2"
i = rownames(sp)[which(sp[,var,drop =F] < 0.05)]
p <- sp[i,var]
r <- (R2.All[[level]] %>% column_to_rownames(level))[i,var]
df <- merge(t(prop[i,, drop =F]), data.obj$meta.dat[,c(var),drop =F], by = 0) %>% na.omit() %>% mutate(bone_marrow_suppression=bone_marrow_suppression2)

round(p, 3) #0.05

df$`p__Firmicutes_A;g__Eubacterium_T`<- sqrt(df$`p__Firmicutes_A;g__Eubacterium_T`)

#remove outlier in No group
p1 <- ggplot(df, aes(x = !!as.name(var), y = !!as.name(i))) +
  geom_boxplot(outlier.shape = NA) +  # Hide outliers in the boxplot
  geom_jitter(width = 0.2, alpha = 0.6, color = "black") +
  showSignificance( c(1.1,1.9), 0.02, 0, "**", textParams = list(size = 7)) +
  theme_classic() + 
  theme(axis.title.x = element_text(size = 14, color = 'black'),
        axis.title.y = element_text(size = 16, color = 'black'),
        axis.text = element_text(size = 14, color = 'black')) +
  ylim(0, 0.025) + 
  labs(x = 'Bone marrow suppression', y = 'Eubacterium_T  sqrt(relative abundance %)')


p2 <- df %>% mutate(presence = as.factor(ifelse(`p__Firmicutes_A;g__Eubacterium_T`>0, 'presence','absence'))) %>% group_by(bone_marrow_suppression) %>% 
  ggplot(aes(x = bone_marrow_suppression)) + geom_bar(aes(fill=presence), position="fill") + theme_classic() + 
  theme(axis.title = element_text(size = 14, color = 'black'),
        legend.text = element_text(size = 14, color = 'black'),
        legend.title  = element_text(size = 14, color = 'black'),
        axis.text = element_text(size = 14, color = 'black'))+
  scale_fill_brewer(palette = 'Dark2') + 
  labs(x = 'Bone marrow suppression', fill = '', y = 'Proportion')
ggarrange(p1, p2, nrow =1)
setwd(mfd)



## ====== Blood outcomes func =====
setwd(wd)
load('./Figure/CancerOnly_func/pathway/DAA_P_R2.RData')

#this one does not load
load("./Result/CancerOnly_func/pathway/data.obj.wk.RData")
level = 'pathway'
sp <- Q.All[[level]][,c(level,"neut_neutropenia2","bone_marrow_suppression2","Hb_anemia2","Pl_thrombocytopenia2")] %>% column_to_rownames(level)
prop <- t(t(data.obj$abund.list[[level]])/colSums(data.obj$abund.list[[level]]))
apply(sp, 2, function(x) sum(x<0.05, na.rm = T))


## ====== Blood outcomes on raryfied data=====
setwd(wd) 
load('Result/CancerOnly/Alpha/bone_marrow_suppression2/Alpha.RData') 
load('Figure/CancerOnly/DAA_P_R2.RData')
head(fit1$result$alpha.diversity )
alpha.df <- merge(fit1$result$alpha.diversity, data.obj.rff$meta.dat, by = 0)
alpha.measure <- "Shannon"

#test significance 
t.test(alpha.df$Shannon ~ alpha.df$bone_marrow_suppression2)
#0.1103

plt1 <- ggplot(alpha.df, aes(x = bone_marrow_suppression2, y = !!as.name(alpha.measure))) + 
  geom_violin(aes(fill = bone_marrow_suppression2),trim = FALSE) + 
  geom_boxplot(width = 0.4, outlier.size = 1)+
  geom_hline(yintercept = median(alpha.df[alpha.df$bone_marrow_suppression2=='Yes',alpha.measure]), linetype = 2) + 
  theme(legend.position = "none") + 
  theme_classic() + 
  scale_fill_manual(name = '',values = c('No' =brewer.pal(8,'Dark2')[2], 'Yes' =brewer.pal(8,'Dark2')[1]))+
  labs(x = 'Bone marrow suppression', fill = '', y = "Shannon diversity index") +
  theme(panel.grid = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size=16, color = 'black'),
        legend.text = element_text(size=16, color = 'black'),
        strip.text = element_text(size=16, color = 'black'),
        axis.text.x = element_text(size=16,color = 'black'),
        axis.text.y = element_text(size=16, color = 'black'))
plt1


dist <- 'WUniFrac'
obj <- cmdscale(dist.obj.rff[[dist]], k=3, eig = T)
pve <- round(obj$eig[1:3]/sum(abs(obj$eig))*100, 1)
eig <- obj$points
colnames(eig) <- c('PC1','PC2','PC3')
yy <- merge(eig, data.obj.rff$meta.dat[,'bone_marrow_suppression2',drop =F], by =0) %>% na.omit()
xlab <- paste0('PC1 (',pve[1],'%)')
ylab <-  paste0('PC2 (',pve[2],'%)')

plt2 <- ggplot(yy) +
  geom_point(size =1.5, aes(x = PC1, y = PC2, fill =  bone_marrow_suppression2, color = bone_marrow_suppression2),  
             alpha = 0.5,shape = 21) +
  ggforce::geom_mark_ellipse(aes(x = PC1, y = PC2, color = bone_marrow_suppression2),expand = unit(2, "mm"), show.legend = T)+
  scale_x_continuous(limits = c(min(yy$PC1)*1.2,max(yy$PC1)*1.2)) +
  scale_y_continuous(limits = c(min(yy$PC2)*1.2,max(yy$PC2)*1.2)) +
  scale_fill_manual(name = '',values = c('No' =brewer.pal(8,'Set1')[2], 'Yes' =brewer.pal(8,'Set1')[1])) +
  scale_color_manual(name = '',values = c('No' =brewer.pal(8,'Set1')[2], 'Yes' =brewer.pal(8,'Set1')[1]))+
  # scale_alpha_discrete(range = c(0.35, 0.9)) +
  labs(x = xlab, y = ylab, fill = '', color = '') +
  theme_bw() +
  theme(text = element_text(size = 20, color = "black"),
        axis.text = element_text(size = 20, color = "black"),
        legend.text = element_text(size = 20, color = "black"),
        axis.title =  element_text(size = 20, color = "black"),
        panel.grid.major = element_blank(),
        panel.border = element_rect(size = 2),
        axis.ticks = element_line(size = 1),
        axis.ticks.length = unit(0.2, 'cm'),
        panel.grid.minor = element_blank()) +
  guides(fill = guide_legend(override.aes = list(size = 4)),color = 'none', alpha = 'none')


setwd(wd)
load('Figure/CancerOnly/DAA_P_R2.RData')
load('Result/CancerOnly/data.obj.wk.RData')
level = 'Genus'
sp <- Q.All[[level]][,c(level,"neut_neutropenia2","bone_marrow_suppression2","Hb_anemia2","Pl_thrombocytopenia2")] %>% column_to_rownames(level)
prop <- t(t(data.obj.rff$abund.list[[level]])/colSums(data.obj.rff$abund.list[[level]]))
apply(sp, 2, function(x) sum(x<0.05, na.rm = T))
var = "bone_marrow_suppression2"
i = rownames(sp)[which(sp[,var,drop =F] < 0.05)]
p <- sp[i,var]
r <- (R2.All[[level]] %>% column_to_rownames(level))[i,var]
df <- merge(t(prop[i,, drop =F]), data.obj$meta.dat[,c(var),drop =F], by = 0) %>% na.omit() %>% mutate(bone_marrow_suppression=bone_marrow_suppression2)
p1 <- ggplot(df, aes(x = !!as.name(var),y = !!as.name(i))) + geom_boxplot() + theme_classic() + 
  theme(axis.title.x = element_text(size = 14, color = 'black'),
        axis.title.y = element_text(size = 14, color = 'black',face = 'italic'),
        axis.text = element_text(size = 14, color = 'black'))+
  scale_y_continuous(trans = sqrt_trans(),
                     breaks = trans_breaks("sqrt", function(x) x^2),
                     labels = trans_format("sqrt", math_format(.x^2))) + 
  labs(x = paste0('Bone marrow suppression (qvalue =',round(p,3),')'), y = 'Eubacterium_T')

p2 <- df %>% mutate(presence = as.factor(ifelse(`p__Firmicutes_A;g__Eubacterium_T`>0, 'presence','absence'))) %>% group_by(bone_marrow_suppression) %>% 
  ggplot(aes(x = bone_marrow_suppression)) + geom_bar(aes(fill=presence), position="fill") + theme_classic() + 
  theme(axis.title = element_text(size = 14, color = 'black'),
        legend.text = element_text(size = 14, color = 'black'),
        legend.title  = element_text(size = 14, color = 'black'),
        axis.text = element_text(size = 14, color = 'black'))+
  scale_fill_brewer(palette = 'Dark2') + 
  labs(x = 'Bone marrow suppression', fill = '', y = 'Proportion')

ggarrange(ggarrange(plt1, plt2, nrow =1, labels = c('A.','B.')), 
          ggarrange(p1, p2, nrow =1, labels = c('C.','D.')), 
          nrow =2)
setwd(mfd)
ggsave('FigureS11(on raryfied).pdf', width = 12, height = 8)


## ==== Eubacteria_T AUC prediction =====
setwd(wd)
load('Figure/CancerOnly/DAA_P_R2.RData')
lapply(Q.All, function(x) sum(x[,c('bone_marrow_suppression2')]<0.1, na.rm = T))
## Only one at genus level associated with outcome

setwd(wd)
tm = load('Result/CancerOnly/data.obj.wk.RData')
rownames(data.obj.rff$abund.list$Genus)[grep("p__Firmicutes_A;g__Eubacterium_T",rownames(data.obj.rff$abund.list$Genus))]
prop <- t(t(data.obj.rff$abund.list$Genus)/colSums(data.obj.rff$abund.list$Genus))
df <- merge(t(prop["p__Firmicutes_A;g__Eubacterium_T",, drop =F]), data.obj.rff$meta.dat[,'bone_marrow_suppression2',drop =F], by = 0) %>% 
  column_to_rownames('Row.names') %>% na.omit
myroc <- roc(df$bone_marrow_suppression2, df$`p__Firmicutes_A;g__Eubacterium_T`)
auc(myroc)
p3 <- ggroc(myroc, colour = "#D95F02", size = 1) + 
  geom_segment(aes(x = 1, xend = 0, y = 0, yend = 1), color="darkgrey", linetype="dashed") + 
  theme_classic() + 
  theme(axis.text = element_text(size = 16, color = 'black'),
        axis.title = element_text(size = 16, color = 'black')) + 
  ggtitle(paste0('ROC Curve ', '(AUC = ', round(auc(myroc),2), ')')) 
ggarrange(p1, p2, p3, nrow =1, labels = c('A.','','B.'))
# setwd(mfd)
# ggsave('FigureSN-8(on raryfied).pdf', width = 16, height = 6)



