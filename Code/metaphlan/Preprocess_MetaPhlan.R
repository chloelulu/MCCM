library(readr)
library(tidyverse)
library(readr)
library(ape)
library(stringr)
library(magrittr)
library(phyloseq)
library(tidyr)
library(dplyr)


file_dir <- dirname(rstudioapi::getSourceEditorContext()$path)
file_dir <- dirname(dirname(dirname(dirname(dirname(file_dir)))))

wd <- file_dir
rd <- paste0(file_dir,'/Result/')

if(!dir.exists(rd)) dir.create(rd)
setwd(wd)
source(paste0(file_dir,"/Code/Submission/MayoOncobiomeStudy/Code/Stats.R"))
try(load_package())

load(file = 'Data/data.obj.raw.core.RData') 

clin_meta <- data.obj$meta.dat

kashyap_strat <- read_table("Data/metaphlan/kashyap_merged_output.txt", skip = 1)
kashyap_strat_split <- kashyap_strat %>%
  separate(clade_name,
           into = c("Kingdom","Phylum","Class","Order",
                    "Family","Genus","Species","Strain"),
           sep = "\\|",
           fill = "right") %>% 
          dplyr::filter(is.na(Strain) & !(Kingdom %in% c("k__Eukaryota","UNCLASSIFIED"))) %>%
  dplyr::select(-Strain)
all.levels <- colnames(kashyap_strat_split)[1:7]

colnames(kashyap_strat_split) <- gsub('_L003$','',colnames(kashyap_strat_split))
clin_meta$X.SampleID <- gsub('_V2$','',clin_meta$X.SampleID)
clin_meta$X.SampleID[!(clin_meta$X.SampleID %in% colnames(kashyap_strat_split))]


colnames(kashyap_strat_split)[!(colnames(kashyap_strat_split) %in% clin_meta$X.SampleID)]


table((clin_meta$X.SampleID %in% colnames(kashyap_strat_split)))
#FALSE  TRUE 
#287  1364  #only the healthies are missing


#add the healthy X.SampleID's to the clin_meta
healthy_BiomeIDs <- clin_meta$BIOME_with_sequencing_data[is.na(clin_meta$X.SampleID)]
clin_meta$X.SampleID[is.na(clin_meta$X.SampleID)] <- as.character(sapply(healthy_BiomeIDs, function(x) colnames(kashyap_strat_split)[grep(x, colnames(kashyap_strat_split))]))

kashyap_strat_split <- kashyap_strat_split[, c(all.levels, clin_meta$X.SampleID[!is.na(clin_meta$X.SampleID)])]
head(kashyap_strat_split)
dim(kashyap_strat_split) #3304 1658

colnames(kashyap_strat_split) <- gsub('_L.*|_CTRL_2','',colnames(kashyap_strat_split))
colnames(kashyap_strat_split) <- gsub('_S.*',"",colnames(kashyap_strat_split))

rownames(data.obj$meta.dat)[!(rownames(data.obj$meta.dat) %in% colnames(kashyap_strat_split))] # all samples exist in data.obj exist in metaphlan table
kashyap_strat_split <- kashyap_strat_split[,c(all.levels,rownames(data.obj$meta.dat))]
dim(kashyap_strat_split)


abund.list <- list()
for(i in 1:7){
  select.level <- all.levels[i]
  bench.level <- all.levels[i+1]
  if(select.level == "Kingdom"){
    abund <- kashyap_strat_split[is.na(kashyap_strat_split[,bench.level]),] %>% 
      dplyr::select(-c(all.levels[-i])) %>% column_to_rownames(select.level)
  }else if(select.level == "Species"){
    abund <- kashyap_strat_split[!is.na(kashyap_strat_split[,all.levels[i]]),] %>% 
      dplyr::select(-c(all.levels[-i])) %>% column_to_rownames(select.level)
  }else{
    abund <- kashyap_strat_split[is.na(kashyap_strat_split[,bench.level]) & !is.na(kashyap_strat_split[,all.levels[i]]),] %>% 
      dplyr::select(-c(all.levels[-i])) %>% column_to_rownames(select.level)
  }
  abund.list[[select.level]] <- abund
}


otu.name <- kashyap_strat_split[!is.na(kashyap_strat_split[,"Species"]),] %>% dplyr::select(all.levels) %>% unique %>% as.data.frame
rownames(otu.name) <- otu.name$Species
dim(otu.name)

dim(abund.list$Species) #1854 1651
names(data.obj)


otu.g <- abund.list$Genus

dim(otu.g)


#writing this to file to use for the coabundance networks
write.csv(abund.list$Species, file = paste0('~/Dropbox/Mayo_RS/R/Oncobiome full cohort analysis/data/Metaphlan_taxonomy_MCCM_species.csv'), row.names = T)
#write species to file as well for TOPOSCORE analysis
write.csv(otu.g, file = paste0('~/Dropbox/Mayo_RS/R/Oncobiome full cohort analysis/data/Metaphlan_taxonomy_MCCM_genus.csv'), row.names = T)


## Reference: https://forum.biobakery.org/t/inquiry-regarding-metaphlan-sgbs-phylogenetic-tree/4442/19
tree <- phyloseq::read_tree("Data/metaphlan/mpa_vJan21_CHOCOPhlAnSGB_202103.nwk")
tree_meta <- read.delim("Data/metaphlan/mpa_vJan21_CHOCOPhlAnSGB_202103_species.txt",header=F)
tree_meta$V1 <- gsub("SGB|_group","",tree_meta$V1)
tree_meta <- tidyr::separate_wider_delim(tree_meta,cols=V2,
                            names=c("V2","Other1","Other2","Other3","Other4","Other5","Other6"), # Separates the synonym names for SGB's into separate variables named "Other#" 
                            delim=",", too_few = "align_start", # When too few synonyms, keep the name in V2, etc. and fill remaining Others with NA's
                            too_many = "drop")   # when more than 1+6 synonyms, the excess is dropped 
tree_meta <- tree_meta %>% dplyr::distinct(V2, .keep_all = TRUE)
tree_meta$V1 <- gsub("EUK","",tree_meta$V1)

tree_meta2 <- subset(tree_meta,V1 %in% tree$tip.label)
tree2 <- ape::drop.tip(tree,setdiff(tree$tip.label,tree_meta2$V1)) 
tree_meta2 <- tree_meta2[match(tree2$tip.label,tree_meta2$V1),]
identical(tree2$tip.label,tree_meta2$V1) 
tree2$tip.label <- gsub('.*\\|s__','s__',tree_meta2$V2 )
mytip <- kashyap_strat$clade_name[grepl("\\|s__",  kashyap_strat$clade_name) & !grepl("\\|t__",  kashyap_strat$clade_name)]
mytip[!(mytip %in% tree2$tip.label)]
tree <- ape::drop.tip(tree2,setdiff(tree2$tip.label,gsub('.*\\|s__','s__',mytip))) 

data.obj <- list(meta.dat = data.obj$meta.dat, otu.tab = abund.list$Species, otu.name = otu.name, tree = tree,abund.list =abund.list)
dist.obj <- construct_distance(data.obj = data.obj, rff = F)
dist.obj[['BC']] <- as.matrix(vegdist(t(data.obj$otu.tab),method="bray"))# Original code contains TSS renormalization step to reduce noise; however metaphlan is renormalized data
identical(rownames(data.obj$meta.dat), rownames(dist.obj[['BC']])) 
identical(rownames(data.obj$meta.dat), rownames(dist.obj$WUniFrac))
identical(rownames(data.obj$meta.dat), colnames(data.obj$otu.tab))
identical(rownames(data.obj$meta.dat), colnames(data.obj$abund.list$Species))
save(data.obj, dist.obj, file = "Data/data.obj.mph.RData")






##############################################
## =======Prepare table for TOPOSCORE=======
##############################################
file_dir <- dirname(rstudioapi::getSourceEditorContext()$path)
file_dir <- dirname(dirname(dirname(dirname(dirname(file_dir)))))
setwd(file_dir)
load('Data/data.obj.mph.RData')

otu.tab <- t(data.obj$otu.tab)
colnames(otu.tab) <- gsub('s__','',colnames(otu.tab))
colnames(otu.tab)[grep('Akk',colnames(otu.tab))]
x <- 1:ncol(data.obj$otu.tab)
blocks <- split(x, cut(seq_along(x), 4, labels = FALSE))
for(i in 1:length(blocks)){
  otu.tab <- as.data.frame(t(data.obj$otu.tab[,blocks[[i]]]))
  colnames(otu.tab) <- gsub('s__','',colnames(otu.tab))
  colnames(otu.tab)[grep('Akkermansia_muciniphila',colnames(otu.tab))]
  otu.tab$Sample_id <- rownames(otu.tab)
  write.csv(otu.tab, file = paste0(file_dir,'/Data/TOPOSCORE/block',i,'.csv'), row.names = F)
}

## ---- Process TOPOSCORE output ----
setwd(file_dir)
load('Data/data.obj.mph.RData')
setwd(paste0(file_dir,'/Result_mph/TOPOSCORE/'))
files <- list.files()
topoMCCM <- NULL
for(file in files){
  topo <- read.csv(file)
  topoMCCM <- rbind(topoMCCM, topo)
}
table(gsub("_S.*","",topoMCCM$Sample_id))
rownames(data.obj$meta.dat)[!(rownames(data.obj$meta.dat) %in%topoMCCM$Sample_id)]
topoMCCM2 <- left_join(data.obj$meta.dat %>% rownames_to_column("Sample_id"),topoMCCM) %>% column_to_rownames('Sample_id')
identical(rownames(topoMCCM2), rownames(data.obj$meta.dat))
data.obj$meta.dat <- topoMCCM2
save(data.obj, dist.obj, file = paste0(file_dir,"/Data/data.obj.mph.RData"))




## Check overlap/diff between 2 pipelines
file_dir <- dirname(rstudioapi::getSourceEditorContext()$path)
file_dir <- dirname(dirname(dirname(dirname(dirname(file_dir)))))

wd <- file_dir
rd <- paste0(file_dir,'/Result_mph/')
setwd(file_dir)
load(file = 'Data/data.obj.mph.RData') 
mph_sp <- rownames(data.obj$otu.tab)
mph_gen <- rownames(data.obj$abund.list$Genus)
gen_abun_mph <- data.obj$abund.list$Genus
sp_abun_mph <- data.obj$abund.list$Species

load(file = 'Data/data.obj.raw.core.RData') 
kk_sp <- rownames(data.obj$otu.tab)
kk_gen <- gsub('.*\\;','',rownames(data.obj$abund.list$Genus))
gen_abun_kk <- t(t(data.obj$abund.list$Genus)/colSums(data.obj$abund.list$Genus))
sp_abun_kk <- t(t(data.obj$otu.tab)/colSums(data.obj$otu.tab))
rownames(gen_abun_kk) <- gsub('.*\\;','',rownames(gen_abun_kk))



#metaphlan d numbers using this taxonomy version
lapply(abund.list, dim)
#$Class
#[1]  140 1651
#$Order
#[1]  169 1651
#$Family
#[1]  219 1651
#$Genus
#[1]  900 1651
#$Species
#[1] 1854 1651


#Kraken d numbers
lapply(data.obj$abund.list, dim)
#$Class
#[1]   74 1651
#$Order
#[1]  170 1651
#$Family
#[1]  391 1651
#$Genus
#[1] 1899 1651
#$Species
#[1] 7839 1651


#@Lu please add some writing on these helpful descriptive statistics
length(intersect(mph_sp, kk_sp))
length(mph_sp[!(mph_sp %in% kk_sp)])
length(kk_sp[!(kk_sp %in% mph_sp)])

length(intersect(mph_gen, kk_gen))
length(mph_gen[!(mph_gen %in% kk_gen)])
length(kk_gen[!(kk_gen %in% mph_gen)])

comm_gen <- intersect(mph_gen, kk_gen)
mean(colSums(gen_abun_mph[comm_gen,]))/100
mean(colSums(gen_abun_kk[comm_gen,]))

comm_sp <- intersect(mph_sp, kk_sp)
mean(colSums(sp_abun_mph[comm_sp,]))/100
mean(colSums(sp_abun_kk[comm_sp,]))

# Taxonomic profiling was performed using two independent tools: MetaPhlAn4 and Kraken2 with Bracken abundance re-estimation. 
# Considerable differences in annotated diversity were observed between the two approaches at both the species and genus levels. 

# At the species level, MetaPhlAn4 identified 1,854 species in total, while Kraken2/Bracken identified 7,839 species. 
# Only 580 species were shared between the two tools, representing 39.3% and 29.7% of MetaPhlAn4 and Kraken2/Bracken species repertoires, respectively. 
# MetaPhlAn4 exclusively detected 1,274 species, whereas Kraken2/Bracken exclusively detected 7,259 species. 

# At the genus level, MetaPhlAn4 annotated 900 genera and Kraken2/Bracken annotated 1,899 genera, with 277 genera shared between the two tools (71.1% and 55.2% of each tool's total, respectively).
# MetaPhlAn4 uniquely identified 623 genera, while Kraken2/Bracken uniquely identified 1,622 genera.



#It may make sense to compare the results from both pipelines only at the genus level. Or add a discussion along these lines.
#from perplexity AI. 

#MetaPhlAn (v3/v4 with the 2021 mpa_vJan21_CHOCOPhlAnSGB database) profiles using clade‑specific marker genes from a curated set of 
#representative genomes (SGBs), not all genomes in GTDB or RefSeq.

#This database is much smaller than a full GTDB r202 genome collection and intentionally trades sensitivity (number of species detected) 
#for higher precision and robustness of relative abundances.

#As a result, many GTDB r202 species either have no SGB representative, are collapsed into a broader SGB, or are only resolvable at 
#genus/higher level in MetaPhlAn outputs.

#Comparative benchmarks consistently find: Kraken(+Bracken) tends to over‑call species (high richness, more false positives), 
#whereas MetaPhlAn3 under‑calls species (lower richness, more false negatives) but with higher precision. https://pmc.ncbi.nlm.nih.gov/articles/PMC10132073/





## Check 2 pipelines results
library(vegan)
load(file = 'Data/data.obj.mph.RData') 
#generate_taxa_barplot_aggregate(data.obj, grp.name = 'Group', ann = 'mph', taxa.levels = c("Genus"))
#generate_taxa_boxplot_aggregate(data.obj, grp.name = 'Group', ann = 'mph', taxa.levels = c("Genus"))
#generate_stacked_barplot(data.obj, grp.name = 'Group', ann = 'mph', taxa.levels = c("Phylum","Genus"), aggre = T, order.taxa = T)
otu1 <- t(data.obj$otu.tab)
load(file = 'Data/data.obj.raw.core.RData') 
#generate_taxa_barplot_aggregate(data.obj, grp.name = 'Group', ann = 'kraken', taxa.levels = c("Genus"))
#generate_taxa_boxplot_aggregate(data.obj, grp.name = 'Group', ann = 'kraken', taxa.levels = c("Genus"))
#generate_stacked_barplot(data.obj, grp.name = 'Group', ann = 'kraken', taxa.levels = c("Phylum","Genus"), aggre = T, order.taxa = T)

otu2 <- t(data.obj$otu.tab)
dim(otu1);dim(otu2)
common_samples <- intersect(rownames(otu1), rownames(otu2))
otu1 <- otu1[common_samples, , drop = FALSE]
otu2 <- otu2[common_samples, , drop = FALSE]

otu1_hel <- decostand(otu1, method = "hellinger")
otu2_hel <- decostand(otu2, method = "hellinger")

d1 <- vegdist(otu1_hel, method = "bray")
d2 <- vegdist(otu2_hel, method = "bray")

pcoa1 <- cmdscale(d1, k = 2, eig = TRUE)
pcoa2 <- cmdscale(d2, k = 2, eig = TRUE)

X <- pcoa1$points
Y <- pcoa2$points

proc <- procrustes(X, Y, symmetric = TRUE)
summary(proc)

plot(proc, kind = 1, main = "Procrustes: mph vs kraken")

set.seed(123)
ptest <- protest(X, Y, permutations = 9999)
print(ptest)



## use species level common species
library(vegan)
load(file = 'Data/data.obj.mph.RData') 
otu1 <- t(data.obj$otu.tab)
load(file = 'Data/data.obj.raw.core.RData') 
otu2 <- t(data.obj$otu.tab)

common_otu <- intersect(colnames(otu1), colnames(otu2))
dim(otu1);dim(otu2)
common_samples <- intersect(rownames(otu1), rownames(otu2))
otu1 <- otu1[common_samples,common_otu , drop = FALSE]
otu2 <- otu2[common_samples,common_otu , drop = FALSE]

otu1_hel <- decostand(otu1, method = "hellinger")
otu2_hel <- decostand(otu2, method = "hellinger")

d1 <- vegdist(otu1_hel, method = "bray")
d2 <- vegdist(otu2_hel, method = "bray")

set.seed(123)
mantel_res <- mantel(d1, d2, method = "spearman", permutations = 999)
print(mantel_res)

pcoa1 <- cmdscale(d1, k = 2, eig = TRUE)
pcoa2 <- cmdscale(d2, k = 2, eig = TRUE)

X <- pcoa1$points
Y <- pcoa2$points

proc <- procrustes(X, Y, symmetric = TRUE)
summary(proc)

plot(proc, kind = 1, main = "Procrustes: mph vs kraken")

set.seed(123)
ptest <- protest(X, Y, permutations = 9999)
print(ptest)



## use genus level common genus
library(vegan)
load(file = 'Data/data.obj.mph.RData') 
otu1 <- t(data.obj$abund.list$Genus)
load(file = 'Data/data.obj.raw.core.RData') 
otu2 <- t(data.obj$abund.list$Genus)
colnames(otu2) <- gsub('.*\\;','',colnames(otu2))

common_otu <- intersect(colnames(otu1), colnames(otu2))
dim(otu1);dim(otu2)
common_samples <- intersect(rownames(otu1), rownames(otu2))
otu1 <- otu1[common_samples,common_otu , drop = FALSE]
otu2 <- otu2[common_samples,common_otu , drop = FALSE]

otu1_hel <- decostand(otu1, method = "hellinger")
otu2_hel <- decostand(otu2, method = "hellinger")

d1 <- vegdist(otu1_hel, method = "bray")
d2 <- vegdist(otu2_hel, method = "bray")

pcoa1 <- cmdscale(d1, k = 2, eig = TRUE)
pcoa2 <- cmdscale(d2, k = 2, eig = TRUE)

X <- pcoa1$points
Y <- pcoa2$points

proc <- procrustes(X, Y, symmetric = TRUE)
summary(proc)

set.seed(123)
ptest <- protest(X, Y, permutations = 9999)
print(ptest)



## use genus level all genus
library(vegan)
load(file = 'Data/data.obj.mph.RData') 
otu1 <- t(data.obj$abund.list$Genus)
load(file = 'Data/data.obj.raw.core.RData') 
otu2 <- t(data.obj$abund.list$Genus)
colnames(otu2) <- gsub('.*\\;','',colnames(otu2))

common_otu <- intersect(colnames(otu1), colnames(otu2))
dim(otu1);dim(otu2)
common_samples <- intersect(rownames(otu1), rownames(otu2))
otu1 <- otu1[common_samples,, drop = FALSE]
otu2 <- otu2[common_samples,mdrop = FALSE]

otu1_hel <- decostand(otu1, method = "rclr")
otu2_hel <- decostand(otu2, method = "rclr")

d1 <- vegdist(otu1_hel, method = "bray")
d2 <- vegdist(otu2_hel, method = "bray")



pcoa1 <- cmdscale(d1, k = 2, eig = TRUE)
pcoa2 <- cmdscale(d2, k = 2, eig = TRUE)

X <- pcoa1$points
Y <- pcoa2$points

proc <- procrustes(X, Y, symmetric = TRUE)
summary(proc)

set.seed(123)
ptest <- protest(X, Y, permutations = 9999)
print(ptest)




### Check maaslin3 and ZicoSeq result 
file_dir <- dirname(rstudioapi::getSourceEditorContext()$path)
file_dir <- dirname(dirname(dirname(dirname(dirname(file_dir)))))

setwd(paste0(file_dir,"/Result_mph/subCancerX_Control"))
dirs <- list.dirs(full.names = F, recursive = F)
for(dir in dirs){
  cat('----------------------------',dir,'----------------------------\n')
  load(paste0(dir,'/DAA/icd10_first_3_name_short/Species_DAA_all2_compare.Rdata'))
  sig_m <- sum((all3$Qvalue_masslin3)<0.1)
  sig_z <- sum(all3$Qvalue_ZicoSeq<0.1)
  sig_zo <- sum(all3$`Qvalue_ZicoSeq(old)`<0.1)
  cat('Maaslin3:',sig_m,'; ')
  cat('ZicoSeq (log2):',sig_z,'; ')
  cat('ZicoSeq (multiple power transform):',sig_zo,'; ')
  load(paste0('../../Result/subCancerX_Control/',dir,'/DAA/icd10_first_3_name_short/icd10_first_3_name_short_ZicoSeq.Rdata'))
  sig_zc <- sum(diff.obj$qv.list[[1]]<0.1)
  # cat('ZicoSeq kraken:',sig_zc,'\n')
  cat('\n')
}


setwd(paste0(file_dir,"/Result_mph_commonKraken/subCancerX_Control"))
dirs <- list.dirs(full.names = F, recursive = F)
for(dir in dirs){
  cat('----------------------------',dir,'----------------------------\n')
  load(paste0(dir,'/DAA/icd10_first_3_name_short/Species_DAA_all2_compare.Rdata'))
  sig_m <- sum((all3$Qvalue_masslin3)<0.1)
  sig_z <- sum(all3$Qvalue_ZicoSeq<0.1)
  sig_zo <- sum(all3$`Qvalue_ZicoSeq(old)`<0.1)
  cat('Maaslin3:',sig_m,'; ')
  cat('ZicoSeq (log2):',sig_z,'; ')
  cat('ZicoSeq (multiple power transform):',sig_zo,'; ')
  load(paste0('../../Result/subCancerX_Control/',dir,'/DAA/icd10_first_3_name_short/icd10_first_3_name_short_ZicoSeq.Rdata'))
  sig_zc <- sum(diff.obj$qv.list[[1]]<0.1)
  # cat('ZicoSeq kraken:',sig_zc,'\n')
  cat('\n')
  
  
  # res <- diff.obj$qv.list[[1]]
  # rownames(res) <- gsub('\\ ','_',rownames(res))
  # res <- merge(all3 %>% column_to_rownames('feature'), as.data.frame(res), by = 0) %>% column_to_rownames("Row.names")
  # res <- res %>% dplyr::select(c(colnames(res)[grep('Qvalu',colnames(res))])) %>% dplyr::rename(`Qvalue(kraken)`=Qvalue)
  # cat('Common species between all:','\n')
  # print(colSums(res<0.1))
}

