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
#check species found in Segata group CRC meta analysis Fusobacterium; not annotated


#Supplementary Table 2 from Thomas et al, 49 positive ones with PvalueAdjusted < 0.05
#6 do not have a species name

  #Fusobacterium nucleatum
  #Parvimonas micra               yes
  #Parvimonas spp.
  #Gemella morbillorum  
  #Peptostreptococcus stomatis    yes
  #Solobacterium moorei           yes
  #Clostridium symbiosum          yes
  #Anaerococcus vaginalis
  #Porphyromonas asaccharolytica
  #Prevotella intermedia
  #Bacteroides fragilis
  #Porphyromonas somerae
  #Anaerococcus obesiensis
  #Porphyromonas uenonis
  #Peptostreptococcus anaerobius
  #Streptococcus constellatus     yes
  #Granulicatella adiacens
  #Methanobrevibacter smithii
  #Eikenella corrodens  
  #Ruminococcus torques           yes ("s__Mediterraneibacter torques")
  #Peptostreptococcus spp.        no but Peptostreptococcus stomatis is there
  #Streptococcus gallolyticus
  #Methanobrevibacter spp.
  #Actinomyces cardiffensis
  #Campylobacter ureolyticus
  #Anaerotruncus spp.             yes optionally "s__Anaerotruncus colihominis"  "s__Anaerotruncus massiliensis"
  #Slackia spp.                   
  #Escherichia coli               yes
  #Campylobacter showae
  #Fusobacterium necrophorum
  #Desulfovibrio desulfuricans
  #Streptococcus dysgalactiae
  #Peptoniphilus harei
  #Bilophila wadsworthia          yes
  #Bilophila spp.
  #Alistipes onderdonkii          yes
  #Alloprevotella tannerae
  #Leptotrichia spp.
  #Eubacterium infirmum
  #Lachnospiraceae bacterium 3 1 57FAA CT1
  #Campylobacter gracilis 
  #Slackia exigua                 yes 
  #Streptococcus tigurinus
  #Fusobacterium mortiferum
  #Eubacterium limosum            yes
  #Bacteroides salyersiae         yes
  #Selenomonas sputigena  
  #Flavonifractor plautii         yes
  #Atopobium rimae


#the rest is NA
rownames(R2_sign_df_comorbidities)[grep("Fus", rownames(R2_sign_df_comorbidities))] #NA #no Fusobacterium
rownames(R2_sign_df_comorbidities)[grep("Parvi", rownames(R2_sign_df_comorbidities))] #"s__Parvimonas micra"
rownames(R2_sign_df_comorbidities)[grep("Pepto", rownames(R2_sign_df_comorbidities))] #"s__Peptostreptococcus stomatis"
rownames(R2_sign_df_comorbidities)[grep("moorei", rownames(R2_sign_df_comorbidities))] #"s__Bulleidia moorei" is new name for Solobacterium moorei
rownames(R2_sign_df_comorbidities)[grep("symbios", rownames(R2_sign_df_comorbidities))] #"s__Clostridium_Q symbiosum"
rownames(R2_sign_df_comorbidities)[grep("fragilis", rownames(R2_sign_df_comorbidities))] #"s__Bacteroides fragilis" AND "s__Bacteroides fragilis_A"
rownames(R2_sign_df_comorbidities)[grep("Strep", rownames(R2_sign_df_comorbidities))] #"s__Streptococcus constellatus"
rownames(R2_sign_df_comorbidities)[grep("Granu", rownames(R2_sign_df_comorbidities))] #"s__Granulicatella adiacens" 
rownames(R2_sign_df_comorbidities)[grep("smithii", rownames(R2_sign_df_comorbidities))] #"s__Methanobrevibacter_A smithii" AND "s__Methanobrevibacter_A smithii_A"
rownames(R2_sign_df_comorbidities)[grep("torques", rownames(R2_sign_df_comorbidities))] #"s__Mediterraneibacter torques"
rownames(R2_sign_df_comorbidities)[grep("Anaerotruncus", rownames(R2_sign_df_comorbidities))] #"s__Anaerotruncus colihominis"  "s__Anaerotruncus massiliensis"
rownames(R2_sign_df_comorbidities)[grep("Escherichia", rownames(R2_sign_df_comorbidities))] #"s__Escherichia coli"  "s__Escherichia coli_C" "s__Escherichia coli_D" "s__Escherichia coli_E" 
rownames(R2_sign_df_comorbidities)[grep("Bilophila", rownames(R2_sign_df_comorbidities))] #"s__Bilophila sp900550745" AND "s__Bilophila wadsworthia"
rownames(R2_sign_df_comorbidities)[grep("Alistipes", rownames(R2_sign_df_comorbidities))] #"s__Alistipes onderdonkii"
rownames(R2_sign_df_comorbidities)[grep("exigua", rownames(R2_sign_df_comorbidities))] #"s__Slackia exigua"
rownames(R2_sign_df_comorbidities)[grep("limosum", rownames(R2_sign_df_comorbidities))] #"s__Eubacterium limosum"
rownames(R2_sign_df_comorbidities)[grep("salyersiae", rownames(R2_sign_df_comorbidities))] #"s__Bacteroides salyersiae"
rownames(R2_sign_df_comorbidities)[grep("plautii", rownames(R2_sign_df_comorbidities))] #"s__Flavonifractor plautii"


#heatmap of these species and the comorbidities
sel_species <- c("s__Parvimonas micra", 
                 "s__Peptostreptococcus stomatis", 
                 "s__Bulleidia moorei",
                 "s__Clostridium_Q symbiosum", 
                 "s__Bacteroides fragilis", 
                 "s__Bacteroides fragilis_A", 
                 "s__Streptococcus constellatus",
                 "s__Granulicatella adiacens",
                 "s__Methanobrevibacter_A smithii",
                 #"s__Methanobrevibacter_A smithii_A",
                 "s__Mediterraneibacter torques",
                 "s__Anaerotruncus colihominis",
                 "s__Anaerotruncus massiliensis",
                 "s__Escherichia coli",
                 "s__Escherichia coli_C",
                 "s__Escherichia coli_D",
                 "s__Escherichia coli_E",
                 "s__Bilophila sp900550745",
                 "s__Bilophila wadsworthia",
                 "s__Alistipes onderdonkii",
                 "s__Slackia exigua",
                 "s__Eubacterium limosum",
                 "s__Bacteroides salyersiae",
                 "s__Flavonifractor plautii")

length(sel_species) #23: but 4x E coli, 2 species of Anaerotruncus, 2 species of B. fragilis. 8 for 3 so 23-5 = 18 were matched


R2_df <- R2_sign_df_comorbidities[sel_species,]
FDR_df <- Q_sign_df_comorbidities[sel_species,]

colnames(R2_df) <- gsub("Elixhauser_", "", colnames(R2_df))


sign_cutoffs <- c(0.2, 0.1, 0.05)

name_vec = colnames(FDR_df)

max.cutoff <- max(sign_cutoffs)
q.cut1 <- min(sign_cutoffs); q.cut2 <- median(sign_cutoffs); q.cut3 <- max.cutoff

grid.text_1 <- '***'
grid.text_2 <- '**'
grid.text_3 <- '*'

scale_max <- max(c(abs(min(R2_df)), max(R2_df)))

heatmap_plot <- Heatmap(R2_df, name = "R2", 
                        heatmap_legend_param = list(direction = "vertical"),
                        col = colorRamp2(c(-scale_max, 0, scale_max), c("#0072B2",'white', "#CC79A7")),
                        column_gap = unit(1, "mm"), 
                        border = TRUE,
                        column_title_rot = 40,
                        column_names_rot = 90,
                        column_names_max_height = unit(10, "cm"),
                        
                        cell_fun = function(j, i, x, y, width, height, fill) {
                          if(FDR_df[i,j] < q.cut1) {
                            grid.text(grid.text_1, x = x, y = y, hjust = 0.5, vjust=0.7)
                          }else if(FDR_df[i,j] < q.cut2){
                            grid.text(grid.text_2, x = x, y = y, hjust = 0.5, vjust=0.7)
                          }else if(FDR_df[i,j] < q.cut3){
                            grid.text(grid.text_3, x = x, y = y, hjust = 0.5, vjust=0.7)
                          }
                        },
                        rect_gp = gpar(col= "white"),
                        show_column_names = T, show_row_names = T,
                        show_column_dend = F,
                        show_row_dend = F,
                        cluster_rows = T,
                        cluster_columns = T,
                        row_names_gp = gpar(fontface = 'italic'),
                        row_names_max_width = max_text_width(rownames(R2_df), gp = gpar(fontsize = 18)))

heatmap_plot


#-------------------------------------------------------------------------------
#no indication for connection of Fusobacterium with CRC in this cohort

#rownames(tax_table)[grep("Fuso", rownames(tax_table))]
fuso_vec <- colSums(tax_table[grep("Fuso", rownames(tax_table)),])

#cbind(rownames(clin_meta), names(fuso_vec))
clin_meta$icd10_first_3_name[grep("colon", clin_meta$icd10_first_3_name)] #37 patients

group_vec_crc <- clin_meta$icd10_first_3_name
group_vec_crc[!group_vec_crc == "Malignant neoplasm of colon"] <- "other_cancer"
group_vec_crc[is.na(group_vec_crc)] <- "healthy"

boxplot(fuso_vec ~ group_vec_crc)
boxplot(fuso_vec ~ group_vec_crc, ylim=c(0, 1000))

#no trend of significance
#summary(lm(fuso_vec ~ group_vec_crc))



#-------------------------------------------------------------------------------
#version with only anemia and weightloss

R2_df <- R2_sign_df_comorbidities[sel_species,]
FDR_df <- Q_sign_df_comorbidities[sel_species,]

sel_cols <- c("Elixhauser_WeightLoss", "Elixhauser_Anemia", "Elixhauser_NeuroOther", "Elixhauser_PUD",
              "Elixhauser_Renal", "Elixhauser_Pulmonary", "Elixhauser_FluidsLytes", "Elixhauser_DMcx", "Elixhauser_Arrhythmia")
#sel_cols %in% colnames(R2_df)

R2_df <- R2_df[sel_species,sel_cols]
FDR_df <- FDR_df[sel_species,sel_cols]

scale_max <- max(c(abs(min(R2_df)), max(R2_df)))

heatmap_plot <- Heatmap(R2_df, name = "R2", 
                        heatmap_legend_param = list(direction = "vertical"),
                        col = colorRamp2(c(-scale_max, 0, scale_max), c("#0072B2",'white', "#CC79A7")),
                        column_gap = unit(1, "mm"), 
                        border = TRUE,
                        column_title_rot = 40,
                        column_names_rot = 90,
                        column_names_max_height = unit(10, "cm"),
                        
                        cell_fun = function(j, i, x, y, width, height, fill) {
                          if(FDR_df[i,j] < q.cut1) {
                            grid.text(grid.text_1, x = x, y = y, hjust = 0.5, vjust=0.7)
                          }else if(FDR_df[i,j] < q.cut2){
                            grid.text(grid.text_2, x = x, y = y, hjust = 0.5, vjust=0.7)
                          }else if(FDR_df[i,j] < q.cut3){
                            grid.text(grid.text_3, x = x, y = y, hjust = 0.5, vjust=0.7)
                          }
                        },
                        rect_gp = gpar(col= "white"),
                        show_column_names = T, show_row_names = T,
                        show_column_dend = F,
                        show_row_dend = F,
                        cluster_rows = T,
                        cluster_columns = T,
                        row_names_gp = gpar(fontface = 'italic'),
                        row_names_max_width = max_text_width(rownames(R2_df), gp = gpar(fontsize = 18)))

heatmap_plot



#-------------------------------------------------------------------------------
#check the control comparison with CRC


CRC_species_names <- zico_res_ctr$R2$Species
CRC_R2 <- zico_res_ctr$R2$"Control-colon"
CRC_Q <- zico_res_ctr$Q$"Control-colon"

CRC_df <- cbind(CRC_R2, CRC_Q)
rownames(CRC_df) <- CRC_species_names

CRC_ctr_res <- CRC_df[which(CRC_species_names %in% sel_species),]
#                                 CRC_R2        CRC_Q
#s__Alistipes onderdonkii           1.551388e-03 4.032861e-01
#s__Anaerotruncus colihominis       1.628171e-02 3.512284e-04 #sign
#s__Anaerotruncus massiliensis      3.292586e-03 3.002663e-01
#s__Bacteroides fragilis            7.519723e-03 1.192916e-01
#s__Bacteroides fragilis_A          1.191454e-05 9.664166e-01
#s__Bacteroides salyersiae         -2.332184e-03 4.255365e-01   #negative effect size
#s__Bilophila sp900550745                     NA           NA
#s__Bilophila wadsworthia          -2.398924e-03 2.761484e-01   #negative effect size
#s__Bulleidia moorei                          NA           NA
#s__Clostridium_Q symbiosum         2.351832e-02 2.542916e-03 #sign
#s__Escherichia coli                7.118934e-02 2.085419e-05 #sign
#s__Escherichia coli_C              5.079932e-02 3.128128e-04 #sign
#s__Escherichia coli_D              4.559603e-02 3.869416e-04 #sign
#s__Eubacterium limosum                       NA           NA
#s__Flavonifractor plautii          1.154693e-02 5.977269e-03 #sign
#s__Granulicatella adiacens                   NA           NA
#s__Mediterraneibacter torques      2.496201e-02 1.044523e-03 #sign
#s__Methanobrevibacter_A smithii    7.388292e-04 6.670935e-01
#s__Peptostreptococcus stomatis               NA           NA
#s__Slackia exigua                            NA           NA
#s__Streptococcus constellatus                NA           NA


#---------------------------------------
#inspect how many species associated to CRC

#remove NAs
CRC_df_narm <-  as.data.frame(CRC_df[!is.na(CRC_df[,"CRC_R2"]),])
#volcano plot
dim(CRC_df_narm); dim(CRC_df) #907, 1293 

plot(CRC_df_narm$CRC_R2, -log10(CRC_df_narm$CRC_Q), xlim=c(-0.1, 0.25), las=1 , 
     xlab="effect size", ylab="-log10(q)")
abline(h=-log10(0.01), col="red", v=0.02)

#sort the species and inspect
length(which(CRC_df_narm$CRC_Q < 0.1)) #458
length(which(CRC_df_narm$CRC_Q < 0.05)) #395
length(which(CRC_df_narm$CRC_Q < 0.01)) #258

length(which(CRC_df_narm$CRC_Q < 0.01 & CRC_df_narm$CRC_R2 > 0.02)) #109


#use effect size
rownames(CRC_df_narm)[which(CRC_df_narm$CRC_Q < 0.01 & CRC_df_narm$CRC_R2 > 0.1)]
#s__Clostridium_Q symbiosum is the only species in here which is also in the Thomas et al meta-analysis


#-------------------------------------------------------------------------------
#Re do the CRC analysis zicoseq only for these samples considering all taxa

#select samples from control and CRC

clin_meta$icd10_first_3_name[is.na(clin_meta$icd10_first_3_name)] <- "healthy"

clin_meta_sel <- clin_meta[clin_meta$icd10_first_3_name %in% c("healthy", "Malignant neoplasm of colon"),]
tax_table_sub <- tax_table[,rownames(clin_meta_sel)]
table(clin_meta_sel$icd10_first_3_name)


comm <- as.matrix(tax_table_sub)
meta <- clin_meta_sel

#impute BMI
meta$BMI[is.na(meta$BMI)] <- mean(meta$BMI, na.rm=T)

comm <- comm[rowSums(comm > 0) > (0.05 * ncol(comm)),] #in at least 5% of samples
dim(tax_table_sub); dim(comm)


#adjusted for Sex and BMI only; Lu can do more corrections.
set.seed(1)
ZicoSeq.obj <- ZicoSeq(meta.dat = meta, feature.dat = comm, 
                       grp.name = 'icd10_first_3_name', feature.dat.type = "count", 
                       #adj.name = c('Batch', 'Sex', 'BMI'), 
                       adj.name = c('Sex', 'BMI'), #done before
                       # Filter to remove rare taxa
                       #prev.filter = 0.001, mean.abund.filter = 0.0001,  
                       #max.abund.filter = 0.002, min.prop = 0, 
                       # Winsorization to replace outliers
                       is.winsor = TRUE, outlier.pct = 0.03, winsor.end = 'top',
                       # Posterior sampling 
                       is.post.sample = TRUE, post.sample.no = 25, 
                       # Use the square-root transformation
                       link.func = list(function (x) x^0.5), stats.combine.func = max,
                       # Permutation-based multiple testing correction
                       perm.no = 999,  strata = NULL, 
                       # Reference-based multiple stage normalization
                       ref.pct = 0.5, stage.no = 6, excl.pct = 0.2,
                       # Family-wise error rate control
                       is.fwer = TRUE, verbose = TRUE, return.feature.dat = T)

ZicoSeq.plot(ZicoSeq.obj, pvalue.type = 'p.adj.fdr', cutoff = 0.0001)

sel_species_underscore <- gsub(" ", "_", sel_species)
sel_species_underscore %in% rownames(comm) #they are all in the input data

ZicoSeq.obj$p.raw[sel_species_underscore]
#s__Parvimonas_micra  s__Peptostreptococcus_stomatis             s__Bulleidia_moorei      s__Clostridium_Q_symbiosum         s__Bacteroides_fragilis 
#       0.002                           0.288                           0.481                           0.001                           0.117 
#s__Bacteroides_fragilis_A   s__Streptococcus_constellatus      s__Granulicatella_adiacens s__Methanobrevibacter_A_smithii   s__Mediterraneibacter_torques 
#       0.926                           0.002                           0.733                           0.615                           0.002 
#s__Anaerotruncus_colihominis   s__Anaerotruncus_massiliensis             s__Escherichia_coli           s__Escherichia_coli_C           s__Escherichia_coli_D 
#       0.001                           0.229                           0.001                           0.001                           0.001 
#s__Escherichia_coli_E        s__Bilophila_sp900550745        s__Bilophila_wadsworthia        s__Alistipes_onderdonkii               s__Slackia_exigua 
#       0.001                           0.097                           0.125                           0.575                           0.005 
#s__Eubacterium_limosum       s__Bacteroides_salyersiae       s__Flavonifractor_plautii 
#       0.001                           0.357                           0.006 


ZicoSeq.obj$p.adj.fdr[sel_species_underscore]
#s__Parvimonas_micra  s__Peptostreptococcus_stomatis             s__Bulleidia_moorei      s__Clostridium_Q_symbiosum         s__Bacteroides_fragilis 
#   0.0150020150                    0.4038183604                    0.6008967645                    0.0062900429                    0.2134984283 
#s__Bacteroides_fragilis_A   s__Streptococcus_constellatus      s__Granulicatella_adiacens s__Methanobrevibacter_A_smithii   s__Mediterraneibacter_torques 
#   0.9561930797                    0.0065828160                    0.7954638538                    0.7044234414                    0.0049498120 
#s__Anaerotruncus_colihominis   s__Anaerotruncus_massiliensis             s__Escherichia_coli           s__Escherichia_coli_C           s__Escherichia_coli_D 
#   0.0026340540                    0.3315671906                    0.0001380691                    0.0009797031                    0.0012579629 
#s__Escherichia_coli_E        s__Bilophila_sp900550745        s__Bilophila_wadsworthia        s__Alistipes_onderdonkii               s__Slackia_exigua 
#   0.0010303190                    0.1921110637                    0.2115518211                    0.6828268745                    0.0164902581 
#s__Eubacterium_limosum       s__Bacteroides_salyersiae       s__Flavonifractor_plautii 
#   0.0002830417                    0.4673350218                    0.0248843885 

length(which(ZicoSeq.obj$p.adj.fdr[sel_species_underscore] < 0.05)) #12 of the 23 corresponding to 9 unique species because of 4 E coli species being significant

sel_species_underscore[which(ZicoSeq.obj$p.adj.fdr[sel_species_underscore] < 0.05)]
#[1] "s__Parvimonas_micra"           "s__Clostridium_Q_symbiosum"    "s__Streptococcus_constellatus" "s__Mediterraneibacter_torques" "s__Anaerotruncus_colihominis" 
#[6] "s__Escherichia_coli"           "s__Escherichia_coli_C"         "s__Escherichia_coli_D"         "s__Escherichia_coli_E"         "s__Slackia_exigua"            
#[11] "s__Eubacterium_limosum"        "s__Flavonifractor_plautii"  

#some are not in the main comparisons showing that they are filtered out so very low abundant and low prevalence

#prevalence of the comorbidities linked to species that are linked to CRC in meta-analysis
#most common are Elixhauser_FluidsLytes, Elixhauser_Arrhythmia, Elixhauser_WeightLoss, Elixhauser_Pulmonary
sapply(sel_cols, function(x) table(clin_meta_sel[,x]))
#           Elixhauser_WeightLoss Elixhauser_Anemia Elixhauser_NeuroOther Elixhauser_PUD Elixhauser_Anemia Elixhauser_Anemia Elixhauser_Renal Elixhauser_Pulmonary Elixhauser_FluidsLytes
#No                     31                33                    36             35                33                33               35                   31                     28
#Yes                     6                 4                     1              2                 4                 4                2                    6                      9
#         Elixhauser_DMcx Elixhauser_Arrhythmia
#No               35                    28
#Yes               2                     9


#--------------------------------------
#plot these significant ones in boxplot

#[1] "s__Parvimonas_micra"           "s__Clostridium_Q_symbiosum"    "s__Streptococcus_constellatus" "s__Mediterraneibacter_torques" "s__Anaerotruncus_colihominis" 
#[6] "s__Escherichia_coli"           "s__Escherichia_coli_C"         "s__Escherichia_coli_D"         "s__Escherichia_coli_E"         "s__Slackia_exigua"            
#[11] "s__Eubacterium_limosum"        "s__Flavonifractor_plautii" 

#not different in the main comparison; too sparse
boxplot(comm["s__Parvimonas_micra",] ~ clin_meta_sel$icd10_first_3_name) #does not look that different
boxplot(comm["s__Parvimonas_micra",] ~ clin_meta_sel$icd10_first_3_name, ylim=c(0,100)) 

#Main one we can make a point of in the manuscript; 
boxplot(comm["s__Clostridium_Q_symbiosum",] ~ clin_meta_sel$icd10_first_3_name) #elevated in CRC. 
#Main one we can make a point of in the manuscript of the top 20 effect sizes shown in main text, now repeated with all significant species
boxplot(comm["s__Clostridium_Q_symbiosum",] ~ clin_meta_sel$icd10_first_3_name, ylim=c(0,2000)) 

#not different in the main comparison; too sparse
boxplot(comm["s__Streptococcus_constellatus",] ~ clin_meta_sel$icd10_first_3_name) #elevated in CRC
boxplot(comm["s__Streptococcus_constellatus",] ~ clin_meta_sel$icd10_first_3_name, ylim=c(0,200)) 

boxplot(comm["s__Mediterraneibacter_torques",] ~ clin_meta_sel$icd10_first_3_name)
boxplot(comm["s__Anaerotruncus_colihominis",] ~ clin_meta_sel$icd10_first_3_name)
boxplot(comm["s__Escherichia_coli",] ~ clin_meta_sel$icd10_first_3_name)
boxplot(comm["s__Escherichia_coli_C",] ~ clin_meta_sel$icd10_first_3_name)
boxplot(comm["s__Escherichia_coli_D",] ~ clin_meta_sel$icd10_first_3_name)
boxplot(comm["s__Escherichia_coli_E",] ~ clin_meta_sel$icd10_first_3_name)
boxplot(comm["s__Slackia_exigua",] ~ clin_meta_sel$icd10_first_3_name)
boxplot(comm["s__Eubacterium_limosum",] ~ clin_meta_sel$icd10_first_3_name)
boxplot(comm["s__Flavonifractor_plautii",] ~ clin_meta_sel$icd10_first_3_name)



#-------------------------------------------------------------------------------
#add an annotation for which species are also significant with CRC in this cohort

R2_df <- R2_sign_df_comorbidities[sel_species,sel_cols]
FDR_df <- Q_sign_df_comorbidities[sel_species,sel_cols]

#sel_cols <- c("Elixhauser_WeightLoss", "Elixhauser_Anemia")

scale_max <- max(c(abs(min(R2_df)), max(R2_df)))


#add bar on the side with significance in CRC comparision on the right
q_value_vec <- ZicoSeq.obj$p.adj.fdr[sel_species_underscore]
coef <- ZicoSeq.obj$coef.list[[1]][4, sel_species_underscore]
r2_value_vec <- ZicoSeq.obj$R2[sel_species_underscore,]
r2_value_vec = sign(coef) * r2_value_vec


# orange         skyblue   bluishgreen        yellow        blue    vermillion reddishpurple 
#"#E69F00"     "#56B4E9"     "#009E73"     "#F0E442"     "#0072B2"     "#D55E00"     "#CC79A7" 


#change the colors here:
col_fun = colorRamp2(c(-max(abs(r2_value_vec)), 0, max(r2_value_vec)), c("#009E73", "white", "#F0E442"))
pch_vec <- rep("", times=length(q_value_vec))
pch_vec[q_value_vec < 0.05] <- "*"
pch_vec[q_value_vec < 0.01] <- "**"


column_ha =  rowAnnotation("CRC vs control comparison" = anno_simple(r2_value_vec, col = col_fun, pch = pch_vec))

lgd_sig_r2 = Legend(pch = c("***","**","*"), type = "points", labels = "<0.2, <0.1, <0.05", background="white")

lgd_r2 = Legend(title = "R2 CRC\nvs. healthy", col_fun = col_fun)
# and one for the significant p-values
lgd_sig = Legend(pch = c("**","*"), type = "points", labels = "<0.05, <0.01", background="white")
# these two self-defined legends are added to the plot by `annotation_legend_list`




R2_df <- t(R2_df)
FDR_df <- t(FDR_df)

row_ha =  columnAnnotation("CRC vs. control comparison" = anno_simple(r2_value_vec, col = col_fun, pch = pch_vec))


colnames(R2_df) <- gsub("s__", "", colnames(R2_df))

heatmap_plot <- Heatmap(R2_df, name = "R2 \ncomorbidities", 
                        heatmap_legend_param = list(direction = "vertical"),
                        #right_annotation = column_ha,
                        top_annotation = row_ha,
                        col = colorRamp2(c(-scale_max, 0, scale_max), c("#0072B2",'white', "#CC79A7")),
                        column_gap = unit(1, "mm"), 
                        border = TRUE,
                        column_title_rot = 40,
                        column_names_rot = 40,
                        column_names_max_height = unit(10, "cm"),
                        
                        cell_fun = function(j, i, x, y, width, height, fill) {
                          if(FDR_df[i,j] < q.cut1) {
                            grid.text(grid.text_1, x = x, y = y, hjust = 0.5, vjust=0.7)
                          }else if(FDR_df[i,j] < q.cut2){
                            grid.text(grid.text_2, x = x, y = y, hjust = 0.5, vjust=0.7)
                          }else if(FDR_df[i,j] < q.cut3){
                            grid.text(grid.text_3, x = x, y = y, hjust = 0.5, vjust=0.7)
                          }
                        },
                        rect_gp = gpar(col= "white"),
                        show_column_names = T, show_row_names = T,
                        show_column_dend = F,
                        show_row_dend = F,
                        cluster_rows = T,
                        cluster_columns = T,
                        column_names_gp = gpar(fontface = 'italic'),
                        row_names_max_width = max_text_width(rownames(R2_df), gp = gpar(fontsize = 18)))


left_margin = 3

packed_legends <- packLegend(lgd_sig_r2, lgd_r2, lgd_sig, direction = "vertical")
heatmap_plot <- draw(heatmap_plot, annotation_legend_list = list(packed_legends), merge_legend=T, 
                     padding = unit(c(0.1, left_margin, 0.1, 0.1), "cm")) 


#save to plot
filename <- "./Figure/RM_figures/CRC comorbidity heatmap.png"
png(filename, width =11.5*330, height = 5*330, res = 330)
heatmap_plot
dev.off()


