require(tidyverse)
require(RColorBrewer)
library(gt)



file_dir <- dirname(rstudioapi::getSourceEditorContext()$path)
file_dir <- dirname(dirname(dirname(dirname(file_dir))))

setwd(file_dir) #"mforge_clean"


load(file = 'Data/data.obj.raw.core.RData') 
clin_meta <- as.data.frame(data.obj$meta.dat)
tax_table <- as.data.frame(data.obj.rff$otu.tab)

#tax_table[1:2,1:2]


#--------------------
#same for species

load("Figure/subCancerX-Ex/DAA_P_R2.RData")

df.q <- Q.All$Species %>% column_to_rownames('Species') %>% na.omit()#%>% dplyr::select(-'Pancancer')
df.r2 <- R2.All$Species %>% column_to_rownames('Species') %>% na.omit()#%>% dplyr::select(-'Pancancer')

cutoff <- 0.1; q.cut1 = 0.01; q.cut2=0.05;q.cut3=0.1

#Q_sig is simply one where a row contains one below cutoff
Q.sig <- df.q[rowSums(df.q <= cutoff) >= 1,]
R2.sig <- df.r2[rowSums(df.q <= cutoff) >= 1,]

effect_size_cutoff <- 0.15
#add the effect size cutoff, only top 15% effect sizes
if (effect_size_cutoff < 1) {
  
  R2_temp <- t(R2.sig)
  Q_temp <- t(Q.sig)
  
  effect_cutoff <- min(sort(apply(abs(R2_temp), 1, max), decreasing=T)[1:nrow(R2_temp)*effect_size_cutoff])
  effect_idx <- apply(abs(R2_temp), 1, function(x) sum(x > effect_cutoff) > 0)
  
  Q.sig <- as.matrix(t(Q_temp)[effect_idx,, drop =F])
  R2.sig <- as.matrix(t(R2_temp)[effect_idx,, drop =F])
}

unique.sig <- rownames(Q.sig)

dim(Q.sig); dim(R2.sig)

rownames(Q.sig) <- gsub('.*;g__|s__','',rownames(Q.sig))
rownames(R2.sig) <- gsub('.*;g__|s__','',rownames(R2.sig))

tmp <- Q.sig
tmp[] <- apply(tmp, 2, function(x) ifelse(x<0.1,'A','B'))
tmp <- as.data.frame(tmp) %>% arrange(across(everything()))

Q.sig <- t(Q.sig[rownames(tmp),])
R2.sig <- t(R2.sig[rownames(tmp), ])

level <- 'Family'
otu.name <- as.data.frame(data.obj$otu.name)

label <- otu.name[otu.name$Species %in% gsub('.*;g__|s__','s__',unique.sig),] %>% dplyr::select(c(level, 'Species')) 

label$Species <- gsub('o__|f__|c__|p__|s__|g__','',label$Species)
label <- label[match(rownames(tmp), label$Species),]
label[,level] <- gsub('o__|f__|c__|p__|g__','',label[,level])




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
  "QAMH01" = "#4daf4a",                     # green
  "Atopobiaceae" = "#f66e9d",               # rose pink
  "CAG-272" = "#f669b8",                    # vibrant pink
  "Micrococcaceae" = "#f66e9b",             # vivid coral
  "UBA1255" = "#f66caa"                     # bold rose
)

ha_combined <- HeatmapAnnotation(
  ` ` = label[, level],
  annotation_height = unit(1, "cm"),
  border = TRUE,
  col = list(` ` = col),
  show_legend = F
)



scale_max <- max(c(abs(min(R2.sig)), max(R2.sig)))


heatmap_plot <- Heatmap(R2.sig, name = "R2 x effect direction", 
        heatmap_legend_param = list(direction = "vertical"),
        col = colorRamp2(c(-scale_max, 0, scale_max), c("#0072B2",'white', "#CC79A7")),
        column_gap = unit(1, "mm"), 
        top_annotation = ha_combined,
        column_split = (label[,level]),
        column_title_rot = 30,
        column_names_rot = 45,
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
        show_column_names = T, show_row_names = T,
        show_column_dend = F,
        show_row_dend = F,
        cluster_rows = T,
        cluster_columns = FALSE,
        column_names_gp = gpar(fontface = 'italic'),
        row_names_max_width = max_text_width(rownames(R2.sig), gp = gpar(fontsize = 18))
) 

left_margin <- 5


png(paste0('Figure/RM_figures/heatmap_all_sig_species_top15perc.png'), width = 20*300, height = 8.5*300, res=300)

heatmap_plot <- draw(heatmap_plot, merge_legend=T, 
                     padding = unit(c(0.1, left_margin, 0.1, 0.1), "cm")) 

dev.off()





