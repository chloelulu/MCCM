require(GUniFrac)
require(tidyverse)
require(RColorBrewer)
library(gt)
require(ComplexHeatmap)
require(circlize)
require(ggpubr)


file_dir <- dirname(rstudioapi::getSourceEditorContext()$path)
file_dir <- dirname(dirname(dirname(dirname(file_dir))))
setwd(file_dir) #"mforge_clean"

#load heatmap from dataframe
source("Code/Submission/MayoOncobiomeStudy/Code/zicoseq_heatmap_from_df.R")


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


#load pathway abundances
load(file = 'Data/data.obj.pathway.RData') 
path_table <- as.data.frame(data.obj$otu.tab)

dim(path_table) #519 1651

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


temp_dir <- paste0(file_dir, "/Result/CancerOnly/DAA/")
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

elix_variables_filenames <- list.files(paste0(temp_dir, elix_variables), pattern = 'Rdata$')
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
head(elix_variables_diff_list$Elixhauser_Alcohol$R2.list$Species,3)
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
#228
table(elix_diff_R2[elix_diff_Q < 0.05] > 0)
#FALSE (negatively associated with ECI)  TRUE (positively associated with ECI)
#                 152                                      76 


#-------------------------------------
#get number of species FDR < 0.1 with positive and negative effect size for components

table_list <- list()
for (i in 1:ncol(R2_sign_df_comorbidities)) {
  R2_vec_temp <- R2_sign_df_comorbidities[,i]
  table_list[[i]] <- table(R2_vec_temp[Q_sign_df_comorbidities[,i] < 0.1] > 0)
}
names(table_list) <- colnames(R2_sign_df_comorbidities)

#sort(R2_sign_df_comorbidities[,"Elixhauser_WeightLoss"][Q_sign_df_comorbidities[,"Elixhauser_WeightLoss"] < 0.1])


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


#-------------------------------------------------------------------------------
#no indication for connection of Fusobacterium with CRC in this cohort

#rownames(tax_table)[grep("Fuso", rownames(tax_table))]
fuso_vec <- colSums(tax_table[grep("Fuso", rownames(tax_table)),])

#cbind(rownames(clin_meta), names(fuso_vec))
clin_meta$icd10_first_3_name_short[grep("colorectal", clin_meta$icd10_first_3_name_short)] #66 patients

group_vec_crc <- as.character(clin_meta$icd10_first_3_name_short)
group_vec_crc[!group_vec_crc %in% c("colorectal", "healthy")] <- "other_cancer"
table(group_vec_crc)


boxplot(fuso_vec ~ group_vec_crc)
boxplot(fuso_vec ~ group_vec_crc, ylim=c(0, 1000))

#no trend of significance
#summary(lm(fuso_vec ~ group_vec_crc))


#-------------------------------------------------------------------------------
#inspect the control comparison with CRC

CRC_species_names <- zico_res_ctr$R2$Species
CRC_R2 <- zico_res_ctr$R2$`Control-colorectal`
CRC_Q <- zico_res_ctr$Q$`Control-colorectal`
  
CRC_df <- cbind(CRC_R2, CRC_Q)
rownames(CRC_df) <- CRC_species_names

CRC_ctr_res <- CRC_df[which(CRC_species_names %in% sel_species),]

#                                   CRC_R2        CRC_Q
#s__Alistipes onderdonkii        -3.470576e-04 5.928387e-02
#s__Anaerotruncus colihominis    -1.372608e-05 1.096647e-02
#s__Anaerotruncus massiliensis   -1.293817e-03 1.615526e-01
#s__Bacteroides fragilis         -7.277317e-06 8.285137e-01
#s__Bacteroides fragilis_A       -1.079772e-03 7.326292e-02
#s__Bacteroides salyersiae       -8.887612e-04 2.429096e-01
#s__Bilophila sp900550745        -8.732328e-04 2.288920e-01
#s__Bilophila wadsworthia        -3.149312e-03 5.107148e-06 #negative effect size
#s__Bulleidia moorei                        NA           NA
#s__Clostridium_Q symbiosum       1.799343e-03 2.394739e-05 #sign
#s__Escherichia coli              7.907771e-03 2.750003e-06 #sign
#s__Escherichia coli_C            6.408873e-03 2.750003e-06 #sign
#s__Escherichia coli_D            1.281676e-02 2.750003e-06 #sign
#s__Eubacterium limosum                     NA           NA
#s__Flavonifractor plautii       -5.336534e-05 1.009859e-03 #sign
#s__Granulicatella adiacens                 NA           NA
#s__Mediterraneibacter torques    1.237988e-05 8.410256e-02 
#s__Methanobrevibacter_A smithii -1.099217e-03 2.312035e-01
#s__Peptostreptococcus stomatis   1.890152e-03 1.099200e-02
#s__Slackia exigua                          NA           NA
#s__Streptococcus constellatus              NA           NA



#--------------------
#inspect how many species associated to CRC

#remove NAs
CRC_df_narm <-  as.data.frame(CRC_df[!is.na(CRC_df[,"CRC_R2"]),])
#volcano plot
dim(CRC_df_narm); dim(CRC_df) #931, 1294

plot(CRC_df_narm$CRC_R2, -log10(CRC_df_narm$CRC_Q), xlim=c(-0.1, 0.25), las=1 , 
     xlab="effect size", ylab="-log10(q)")
abline(h=-log10(0.01), col="red", v=0.02)

#sort the species and inspect
length(which(CRC_df_narm$CRC_Q < 0.1)) #661
length(which(CRC_df_narm$CRC_Q < 0.05)) #591
length(which(CRC_df_narm$CRC_Q < 0.01)) #465

length(which(CRC_df_narm$CRC_Q < 0.1 & CRC_df_narm$CRC_R2 > 0)) #135


#use effect size
rownames(CRC_df_narm)[which(CRC_df_narm$CRC_Q < 0.01 & CRC_df_narm$CRC_R2 > 0.05)]
#[1] "s__Anaerobutyricum sp900016875"      "s__Anaerostipes caccae"              "s__Anaerostipes sp000508985"        
#[4] "s__Bifidobacterium dentium"          "s__Blautia coccoides"                "s__Blautia hominis"                 
#[7] "s__Blautia sp001304935"              "s__Blautia sp003287895"              "s__Blautia sp900120295"             
#[10] "s__Blautia sp900539145"              "s__Blautia sp900541955"              "s__Blautia sp900556555"             
#[13] "s__Blautia_A caecimuris"             "s__Blautia_A sp900551465"            "s__Catenibacillus sp900557175"      
#[16] "s__Clostridium_AP scindens"          "s__Clostridium_AQ innocuum"          "s__Enterocloster aldenensis"        
#[19] "s__Enterocloster clostridioformis"   "s__Enterocloster clostridioformis_A" "s__Erysipelatoclostridium ramosum"  
#[22] "s__Escherichia coli"                 "s__Escherichia coli_C"               "s__Escherichia coli_D"              
#[25] "s__Escherichia dysenteriae"          "s__Escherichia fergusonii"           "s__Escherichia flexneri"            
#[28] "s__Escherichia sp002965065"          "s__Escherichia sp005843885"          "s__Faecalimonas umbilicata"         
#[31] "s__Fusicatenibacter sp900543115"     "s__Klebsiella pneumoniae"            "s__Ruminococcus_A sp000432335"      
#[34] "s__Ruminococcus_B gnavus"            "s__Ruminococcus_B sp900544395"       "s__Streptococcus vestibularis"      
#[37] "s__UMGS403 sp900540275"              "s__Veillonella parvula_A"   


#-------------------------------------------------------------------------------
#Re do the CRC analysis zicoseq only for these samples considering all taxa

#select samples from control and CRC
clin_meta_sel <- as.data.frame(clin_meta[clin_meta$icd10_first_3_name_short %in% c("healthy", "colorectal"),])
clin_meta_sel$CRC_variable <- as.character(clin_meta_sel$icd10_first_3_name_short)
clin_meta_sel$CRC_variable[!clin_meta_sel$CRC_variable == "healthy"] <- "CRC"
table(clin_meta_sel$CRC_variable)
clin_meta_sel$CRC_variable <- relevel(as.factor(clin_meta_sel$CRC_variable), ref = "healthy")

tax_table_sub <- tax_table[,rownames(clin_meta_sel)]

comm <- as.matrix(tax_table_sub)
meta <- clin_meta_sel

#impute BMI
meta$BMI[is.na(meta$BMI)] <- mean(meta$BMI, na.rm=T)
meta$Age[is.na(meta$Age)] <- mean(meta$Age, na.rm=T)

comm <- comm[rowSums(comm > 0) > (0.05 * ncol(comm)),] #in at least 5% of samples
dim(tax_table_sub); dim(comm)
#[1] 7839  353
#[1] 3008  353


#adjusted for Sex and BMI only because it is with control
set.seed(1)
ZicoSeq.obj <- ZicoSeq(meta.dat = meta, feature.dat = comm, 
                       grp.name = 'CRC_variable', feature.dat.type = "count", 
                       adj.name = c('Sex', 'BMI', 'Age'), #done before
                       # Filter to remove rare taxa
                       #prev.filter = 0.001, mean.abund.filter = 0.0001,  
                       #max.abund.filter = 0.002, min.prop = 0, 
                       # Winsorization to replace outliers
                       is.winsor = TRUE, outlier.pct = 0.03, winsor.end = 'top',
                       # Posterior sampling 
                       is.post.sample = TRUE, post.sample.no = 25, 
                       # Use the square-root transformation
                       link.func = list(function (x) x^0.0625, function (x) x^0.5), stats.combine.func = max,
                       # Permutation-based multiple testing correction
                       perm.no = 999,  strata = NULL, 
                       # Reference-based multiple stage normalization
                       ref.pct = 0.5, stage.no = 6, excl.pct = 0.2,
                       # Family-wise error rate control
                       is.fwer = TRUE, verbose = TRUE, return.feature.dat = T)

save(ZicoSeq.obj, file = 'Data/ZicoSeq_CRCvsHealthy_Figure4H.Rdata')
load('Data/ZicoSeq_CRCvsHealthy_Figure4H.Rdata')

ZicoSeq.plot(ZicoSeq.obj, pvalue.type = 'p.adj.fdr', cutoff = 0.0001)

sel_species_underscore <- gsub(" ", "_", sel_species)
sel_species_underscore %in% rownames(comm) #they are all in the input data

ZicoSeq.obj$p.raw[sel_species_underscore]
#s__Parvimonas_micra  s__Peptostreptococcus_stomatis             s__Bulleidia_moorei 
#         0.001                           0.028                           0.128 
#s__Clostridium_Q_symbiosum         s__Bacteroides_fragilis       s__Bacteroides_fragilis_A 
#         0.004                           0.600                           0.027 
#s__Streptococcus_constellatus      s__Granulicatella_adiacens s__Methanobrevibacter_A_smithii 
#         0.031                           0.003                           0.077 
#s__Mediterraneibacter_torques    s__Anaerotruncus_colihominis   s__Anaerotruncus_massiliensis 
#         0.335                           0.165                           0.046 
#s__Escherichia_coli           s__Escherichia_coli_C           s__Escherichia_coli_D 
#         0.001                           0.010                           0.001 
#s__Escherichia_coli_E        s__Bilophila_sp900550745        s__Bilophila_wadsworthia 
#         0.030                           0.051                           0.001 
#s__Alistipes_onderdonkii               s__Slackia_exigua          s__Eubacterium_limosum 
#         0.002                           0.079                           0.018 
#s__Bacteroides_salyersiae       s__Flavonifractor_plautii 
#         0.148                           0.209 


ZicoSeq.obj$p.adj.fdr[sel_species_underscore]
#s__Parvimonas_micra  s__Peptostreptococcus_stomatis             s__Bulleidia_moorei 
#   1.456001e-03                    7.690415e-02                    2.046636e-01 
#s__Clostridium_Q_symbiosum         s__Bacteroides_fragilis       s__Bacteroides_fragilis_A 
#   4.380765e-03                    6.254062e-01                    3.835248e-02 
#s__Streptococcus_constellatus      s__Granulicatella_adiacens s__Methanobrevibacter_A_smithii 
#   5.775139e-02                    8.786912e-03                    1.192115e-01 
#s__Mediterraneibacter_torques    s__Anaerotruncus_colihominis   s__Anaerotruncus_massiliensis 
#   4.313220e-01                    2.476932e-01                    6.427977e-02 
#s__Escherichia_coli           s__Escherichia_coli_C           s__Escherichia_coli_D 
#   1.091579e-04                    1.467693e-02                    6.326599e-04 
#s__Escherichia_coli_E        s__Bilophila_sp900550745        s__Bilophila_wadsworthia 
#   6.066656e-02                    1.065112e-01                    1.609326e-05 
#s__Alistipes_onderdonkii               s__Slackia_exigua          s__Eubacterium_limosum 
#   4.938032e-03                    1.444446e-01                    4.001412e-02 
#s__Bacteroides_salyersiae       s__Flavonifractor_plautii 
#   1.798058e-01                    2.743541e-01 


length(which(ZicoSeq.obj$p.adj.fdr[sel_species_underscore] < 0.1)) #14 of the 23 corresponding to 11 unique species because of 4 E coli species being significant

sel_species_underscore[which(ZicoSeq.obj$p.adj.fdr[sel_species_underscore] < 0.1)]


#some are not in the main comparisons showing that they are filtered out so very low abundant and low prevalence

#prevalence of the comorbidities linked to species that are linked to CRC in meta-analysis
#most common are Elixhauser_FluidsLytes, Elixhauser_Arrhythmia, Elixhauser_WeightLoss, Elixhauser_Pulmonary
sapply(sel_cols, function(x) table(clin_meta_sel[clin_meta_sel$CRC_variable == "CRC",x]))
#             Elixhauser_WeightLoss Elixhauser_Anemia Elixhauser_NeuroOther Elixhauser_PUD Elixhauser_Renal Elixhauser_Pulmonary
# No                     59                62                    65             63               63                   57
# Yes                     7                 4                     1              3                3                    9
#             Elixhauser_FluidsLytes Elixhauser_DMcx Elixhauser_Arrhythmia
#No                      54              63                    56
#Yes                     12               3                    10


#-------------------------------------------------------------------------------
#add an annotation for which species are also significant with CRC in this cohort

R2_df <- R2_sign_df_comorbidities[sel_species,sel_cols]
FDR_df <- Q_sign_df_comorbidities[sel_species,sel_cols]

#sel_cols <- c("Elixhauser_WeightLoss", "Elixhauser_Anemia")

scale_max <- max(c(abs(min(R2_df)), max(R2_df)))

#add bar on the side with significance in CRC comparision on the right
q_value_vec <- ZicoSeq.obj$p.adj.fdr[sel_species_underscore]
coef <- ZicoSeq.obj$coef.list[[1]]["CRC_variableCRC", sel_species_underscore]
r2_value_vec <- ZicoSeq.obj$R2[sel_species_underscore,]

#this is needed when multiple link functions are used
r2_value_vec <- cbind(Func1 = apply(r2_value_vec, 1, function(x) x[which.max(x)]))

r2_value_vec = sign(coef) * r2_value_vec


# orange         skyblue   bluishgreen        yellow        blue    vermillion reddishpurple 
#"#E69F00"     "#56B4E9"     "#009E73"     "#F0E442"     "#0072B2"     "#D55E00"     "#CC79A7" 


#change the colors here:
col_fun = colorRamp2(c(-max(abs(r2_value_vec)), 0, max(r2_value_vec)), c("#009E73", "white", "#F0E442"))
pch_vec <- rep("", times=length(q_value_vec))
pch_vec[q_value_vec < 0.1] <- "*"
pch_vec[q_value_vec < 0.05] <- "**"
pch_vec[q_value_vec < 0.01] <- "***"


column_ha =  rowAnnotation("CRC vs healthy comparison" = anno_simple(r2_value_vec, col = col_fun, pch = pch_vec))

lgd_sig_r2 = Legend(pch = "", type = "points", labels = "  * <0.10\n ** <0.05\n*** <0.01", 
                    background="white")
lgd_r2 = Legend(title = "R2 CRC\nvs. healthy", col_fun = col_fun)


# and one for the significant p-values
lgd_sig = Legend(pch = "", type = "points", labels = "  * <0.10\n ** <0.05\n*** <0.01", 
                 background="white")


R2_df <- t(R2_df)
FDR_df <- t(FDR_df)


row_ha = columnAnnotation(
  "CRC vs. healthy comparison" = anno_simple(
    r2_value_vec[,1], 
    col = col_fun, 
    pch = pch_vec, 
    height = unit(0.85, "cm"),
    pt_size = unit(0.43, "cm")
  )
)                          


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
png(filename, width =11.5*330, height = 5.25*330, res = 330)
heatmap_plot
dev.off()




#save to svg
filename <- "./Figure/RM_figures/CRC comorbidity heatmap.svg"
svg(filename, width =11.5, height = 5.25)
heatmap_plot
dev.off()





