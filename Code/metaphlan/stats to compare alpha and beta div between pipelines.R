#Statistically comparing alpha and beta diversity between Kraken and Metaphlan

library(vegan)


file_dir <- dirname(rstudioapi::getSourceEditorContext()$path)
file_dir <- dirname(dirname(dirname(dirname(dirname(file_dir)))))
setwd(file_dir)


#-----------------
#load Metaphlan

#beta
load("Result_mph/PanCancer/data.obj.wk.RData")
wunifrac_mph <- dist.obj$WUniFrac

#alpha
load("Result_mph/PanCancer/Alpha/Group/Alpha.RData")
shannon_mph <- alpha.obj1$Shannon
names(shannon_mph) <- data.obj2$meta.dat$BIOME_with_sequencing_data


#-----------------
#load Kraken

#beta
load("Result/PanCancer/data.obj.wk.RData")
wunifrac_kra <- dist.obj$WUniFrac

#alpha
load("Result/PanCancer/Alpha/Group/Alpha.RData")
shannon_kra <- alpha.obj1$Shannon
names(shannon_kra) <- data.obj.rff2$meta.dat$BIOME_with_sequencing_data


#-------------------------------------------------------------------------------
#Alpha div correlation

shannon_mph_sub <- shannon_mph[names(shannon_kra)]

table(names(shannon_mph_sub) == names(shannon_kra))
#TRUE 
#1644 

cor.test(as.numeric(shannon_mph_sub), as.numeric(shannon_kra), method = "pearson")
#p-value < 2.2e-16, cor 0.9364668


#-------------------------------------------------------------------------------
#Beta div correlation

table(colnames(wunifrac_mph) == colnames(wunifrac_kra))
table(rownames(wunifrac_mph) == rownames(wunifrac_kra))
#1651

#-----------------
#mantel test

wu1 <- wunifrac_mph
wu2 <- wunifrac_kra

# wu1, wu2 are distance matrices or 'dist' objects from pipeline 1 and 2
# make sure row/col names are sample IDs, and they match
common_ids <- intersect(rownames(wu1), rownames(wu2))
wu1_sub <- as.dist(wu1[common_ids, common_ids])
wu2_sub <- as.dist(wu2[common_ids, common_ids])


set.seed(1)
mantel_res <- mantel(wu1_sub, wu2_sub, method = "spearman",
                     permutations = 999)
mantel_res
#with 999 permutations
#Mantel statistic r: 0.9518 
#Significance: 0.001 





