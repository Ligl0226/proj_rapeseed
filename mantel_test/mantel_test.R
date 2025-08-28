library(vegan)
library(phyloseq)



################################################################################
#             mantel test for root RNA-seq and ionome
################################################################################

# read root RNA-seq normalized table
vsd = readRDS("~/Brassica/rna_seq/intermediate_data/final/counts_filtered_vst.RDS")
# calculate euclidean distance for root RNA-seq
rna_euc_dist = dist(t(assay(vsd)))
euc_dist = as.matrix(rna_euc_dist, labels = TRUE)
# read ionome data
ionData = read.table("~/Brassica/ionome/data/ionome_data_table.txt")

rna_ion_res = c()
# mantel test for each ion and root RNA distance
for (i in 4:16) {
  dist_ion = dist(ionData[, i])
  mantel_res = vegan::mantel(euc_dist, dist_ion)
  rna_ion_res = rbind(rna_ion_res, cbind(ion = colnames(ionData[i]), r = mantel_res$statistic, pval = mantel_res$signif))
}


#############################################################################
#             mantel test for rhizosphere and ionome
#############################################################################

# read rhizosphere normalized counts table
rhizo_vst_norm = readRDS("~/Brassica/microbiome/intermediate_data/rhizo_vst_normalization_with_ion.RDS")
# calculate bray-curtis distance for rhizosphere microbiome
rhizo_dist = phyloseq::distance(rhizo_vst_norm, method = "bray") 

mantel.res = c()
# mantel test for each ion and rhizosphere distance
for (i in colnames(sample_data(rhizo_vst_norm)[, 13:25])) {
  ion_dist = dist(sample_data(rhizo_vst_norm)[, i])
  res = mantel(rhizo_dist, ion_dist, permutations = 999)
  print(res)
  mantel.res = rbind.data.frame(mantel.res, cbind.data.frame(ion = i, r = res$statistic, pval = res$signif))
  
}






