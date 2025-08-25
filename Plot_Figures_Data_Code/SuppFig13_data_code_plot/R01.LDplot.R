##==================================================================##
##--Description: LD plot for 
##--Maintainer: GuoLiang Li(lig@ipk-gatersleben.de)
##--Date: 2024.01.11
##==================================================================##
dir()
getwd()
library(data.table)
library(RColorBrewer)
library(scatterplot3d)
library(ggrepel)
library(showtext);setEPS()

stat.bin.WGS.path <- "./11.WGS_maf5het95/WGS_175x345289_LDdecay.bin"
stat.bin.RNAseqKF.path <- "./12.RNAseqKF_maf5het95/RNAseq_KF_175x239172_LDdecay.bin"
stat.bin.RNAseqYL.path <- "./13.RNAseqYL_maf5het95/RNAseq_YL_175x292839_LDdecay.bin"

file.name.list <- list(WGS_SNPs=stat.bin.WGS.path,
                       RNAseq_SNPs_KF=stat.bin.RNAseqKF.path,
                       RNAseq_SNPs_YL=stat.bin.RNAseqYL.path)

##==========================================================##
##---loading the LD results
##==========================================================##
stat.bin <- list()
for(i in 1:length(file.name.list)){
  stat.bin[[i]] <- read.table(file.name.list[[i]],header = F,sep = "\t",stringsAsFactors = F)
  colnames(stat.bin[[i]]) <- c("dist","meanR2","meanDp","SumR2","SumDp","PairsNum")
}
names(stat.bin) <- names(file.name.list)

length(stat.bin)
names(stat.bin)

cat("dim(stat.bin$WGS_SNPs):");dim(stat.bin$WGS_SNPs)
cat("head(stat.bin$WGS_SNPs):\n");head(stat.bin$WGS_SNPs)
cat("tail(stat.bin$WGS_SNPs):\n");tail(stat.bin$WGS_SNPs)

cat("dim(stat.bin$RNAseq_SNPs_KF):");dim(stat.bin$RNAseq_SNPs_KF)
cat("head(stat.bin$RNAseq_SNPs_KF):\n");head(stat.bin$RNAseq_SNPs_KF)
cat("tail(stat.bin$RNAseq_SNPs_KF):\n");tail(stat.bin$RNAseq_SNPs_KF)

cat("dim(stat.bin$RNAseq_SNPs_YL):");dim(stat.bin$RNAseq_SNPs_YL)
cat("head(stat.bin$RNAseq_SNPs_YL):\n");head(stat.bin$RNAseq_SNPs_YL)
cat("tail(stat.bin$RNAseq_SNPs_YL):\n");tail(stat.bin$RNAseq_SNPs_YL)

legendlist <- names(stat.bin)
colorlist <- c("#E41A1C","#377EB8","#4DAF4A")
# colorlist <- c("#E41A1C","#377EB8","#4DAF4A","#984EA3")
# colorlist <- c("#E41A1C","#ff7f00","#1f78b4","#33a02c")
##==========================================================##
##---plot based on the stat.bin form the results of PopLDdecay
##==========================================================##
all.equal(stat.bin[[1]][,1],stat.bin[[2]][,1])
all.equal(stat.bin[[1]][,1],stat.bin[[3]][,1])


df <- data.frame(x = stat.bin[[1]][,1]/1000, 
                 WGS_SNPs = stat.bin[[1]][,2], 
                 RNAseq_SNPs_KF = stat.bin[[2]][,2],
                 RNAseq_SNPs_YL = stat.bin[[3]][,2])
head(df);tail(df);dim(df)

ggp <- ggplot(data = df, aes(x = x)) +
  geom_line(aes(y = WGS_SNPs, colour = "WGS_SNPs"), size = 1, linetype = 'solid') +
  geom_line(aes(y = RNAseq_SNPs_KF, colour = "RNAseq_SNPs_KF"), size = 1, linetype = 'solid', alpha = 0.7) +
  geom_line(aes(y = RNAseq_SNPs_YL, colour = "RNAseq_SNPs_YL"), size = 1, linetype = 'solid', alpha = 0.7) +
  scale_colour_manual("HyW",breaks = c("WGS_SNPs","RNAseq_SNPs_KF","RNAseq_SNPs_YL"),
                      values = c("#ff7f00","#1f78b4","#33a02c")) +
  labs(title = 'LD', x = 'Distance(Kb)', y = expression(r^{2})) +
  theme_minimal()+
  theme_bw(base_rect_size = 1.0) + 
  theme(axis.text.x = element_text(size=12,colour = "black"),
        axis.text.y = element_text(size=12,colour = "black"),
        axis.title.x = element_text(size=12), 
        axis.title.y = element_text(size=12),
        title = element_text(size=12,hjust = 0.5),
        legend.text = element_text(size=12)) 


pdf(file = "zRes01.LD.plot.pdf", width = 6.5, height = 4)
showtext_begin()
ggp
showtext_end()
dev.off()

