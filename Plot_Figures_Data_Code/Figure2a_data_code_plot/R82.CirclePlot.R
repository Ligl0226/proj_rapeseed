##==================================================================##
##--Description: Circle plot for results of GWAS, eGWAS and TWAS
##--Highlight1: 203 ASVs
##--Highlight2: results with update
##--Maintainer: GuoLiang Li(lig@ipk-gatersleben.de)
##--Date: 2023.10.14/15/16
##==================================================================##
# setwd()
dir();getwd()
library(openxlsx)
library(stringr)
library(data.table)
library(dplyr)
library(CMplot)
library(qs)
library(circlize)
library(RColorBrewer)
library(ggplot2)
library(scales)

KF_GWASres_Path <- "./Res02.GWAS_KF_SigSNPs2QTL_50kb.oneGeneNearest.xlsx"
YL_GWASres_Path <- "./Res02.GWAS_YL_SigSNPs2QTL_50kb.oneGeneNearest.xlsx"

KF_eQTL_hotspot_Path <- "./KF_eQTL.order.adj.bed.100000.bc.un.track.bed"
YL_eQTL_hotspot_Path <- "./YL_eQTL.order.adj.bed.100000.bc.un.track.bed"

KF_TWASres_Path <- "./zRes02.KF_TWAS_SigGene_Summary.FDR.xlsx"
YL_TWASres_Path <- "./zRes12.YL_TWAS_SigGene_Summary.FDR.xlsx"

GenomeInfoPath <- "./Brassica_napus_v4.1.chr.size.txt"

## change the chromosome name
ACxx <- c("A01","A02", "A03", "A04", "A05", "A06",
          "A07", "A08", "A09", "A10", "C01", "C02",
          "C03", "C04", "C05", "C06", "C07", "C08",
          "C09")
chrACxx <- c("chrA01","chrA02", "chrA03", "chrA04", "chrA05", "chrA06",
             "chrA07", "chrA08", "chrA09", "chrA10", "chrC01", "chrC02",
             "chrC03", "chrC04", "chrC05", "chrC06", "chrC07", "chrC08",
             "chrC09")
chrxx <- paste("chr",sprintf("%02d",1:19),sep = "")
ChrName <- data.frame(ACxx,chrACxx,chrxx)


##==================================================================##
##--loading the GWAS results
##==================================================================##
## KF
KF_GWASres_Ori <- read.xlsx(KF_GWASres_Path)
head(KF_GWASres_Ori);tail(KF_GWASres_Ori);dim(KF_GWASres_Ori)  ## 99

KF_GWASres_OK <- data.frame(chr=KF_GWASres_Ori$ChrAC,
                            start=KF_GWASres_Ori$QTLstart,
                            end=KF_GWASres_Ori$QTLend,
                            logpval=-log10(KF_GWASres_Ori$P_value),
                            colvalue=1)
head(KF_GWASres_OK);tail(KF_GWASres_OK);dim(KF_GWASres_OK)  


## YL
YL_GWASres_Ori <- read.xlsx(YL_GWASres_Path)
head(YL_GWASres_Ori);tail(YL_GWASres_Ori);dim(YL_GWASres_Ori)  ## 317

YL_GWASres_OK <- data.frame(chr=YL_GWASres_Ori$ChrAC,
                            start=YL_GWASres_Ori$QTLstart,
                            end=YL_GWASres_Ori$QTLend,
                            logpval=-log10(YL_GWASres_Ori$P_value),
                            colvalue=2)
head(YL_GWASres_OK);tail(YL_GWASres_OK);dim(YL_GWASres_OK)  


##==================================================================##
##--loading the eQTL results (eQTL hotspots)
##==================================================================##
## KF
KF_eQTL_Hotspot_Ori <- read.table(KF_eQTL_hotspot_Path,sep = "\t",header = F,skip = 1)
KF_eQTL_Hotspot_Ori <- KF_eQTL_Hotspot_Ori[,c(1,2,3,10)]
colnames(KF_eQTL_Hotspot_Ori) <- c("chr","start","end","pvalue")
KF_eQTL_Hotspot_Ori <- KF_eQTL_Hotspot_Ori[with(KF_eQTL_Hotspot_Ori, order(chr, start)),]
rownames(KF_eQTL_Hotspot_Ori) <- NULL
dim(KF_eQTL_Hotspot_Ori);head(KF_eQTL_Hotspot_Ori);tail(KF_eQTL_Hotspot_Ori)  ## 317  10
summary(KF_eQTL_Hotspot_Ori$pvalue)

## change the chr names
for(i in 1:length(chrxx)){
  ## i=1
  chrname1i <- chrxx[i]
  chrname2i <- ACxx[i]
  KF_eQTL_Hotspot_Ori[KF_eQTL_Hotspot_Ori$chr==chrname1i,"chr"] <- chrname2i
}
dim(KF_eQTL_Hotspot_Ori);head(KF_eQTL_Hotspot_Ori);tail(KF_eQTL_Hotspot_Ori)
# as.data.frame(table(KF_eQTL_Hotspot_Ori$chr))

KF_eQTL_Hotspot_OK <- data.frame(chr=KF_eQTL_Hotspot_Ori$chr,
                                 start=KF_eQTL_Hotspot_Ori$start,
                                 end=KF_eQTL_Hotspot_Ori$end,
                                 logpval=-log10(KF_eQTL_Hotspot_Ori$pvalue),
                                 colvalue=1)
dim(KF_eQTL_Hotspot_OK);head(KF_eQTL_Hotspot_OK);tail(KF_eQTL_Hotspot_OK)   ## 317
KF_eQTL_Hotspot_OK[is.infinite(KF_eQTL_Hotspot_OK$logpval),"logpval"] <- max(KF_eQTL_Hotspot_OK[!is.infinite(KF_eQTL_Hotspot_OK$logpval),"logpval"])


## YL
YL_eQTL_Hotspot_Ori <- read.table(YL_eQTL_hotspot_Path,sep = "\t",header = F,skip = 1)
YL_eQTL_Hotspot_Ori <- YL_eQTL_Hotspot_Ori[,c(1,2,3,10)]
colnames(YL_eQTL_Hotspot_Ori) <- c("chr","start","end","pvalue")
YL_eQTL_Hotspot_Ori <- YL_eQTL_Hotspot_Ori[with(YL_eQTL_Hotspot_Ori, order(chr, start)),]
rownames(YL_eQTL_Hotspot_Ori) <- NULL
dim(YL_eQTL_Hotspot_Ori);head(YL_eQTL_Hotspot_Ori);tail(YL_eQTL_Hotspot_Ori)  ## 380  10
summary(YL_eQTL_Hotspot_Ori$pvalue)

## change the chr names
for(i in 1:length(chrxx)){
  ## i=1
  chrname1i <- chrxx[i]
  chrname2i <- ACxx[i]
  YL_eQTL_Hotspot_Ori[YL_eQTL_Hotspot_Ori$chr==chrname1i,"chr"] <- chrname2i
}
dim(YL_eQTL_Hotspot_Ori);head(YL_eQTL_Hotspot_Ori);tail(YL_eQTL_Hotspot_Ori)
# as.data.frame(table(YL_eQTL_Hotspot_Ori$chr))

YL_eQTL_Hotspot_OK <- data.frame(chr=YL_eQTL_Hotspot_Ori$chr,
                                 start=YL_eQTL_Hotspot_Ori$start,
                                 end=YL_eQTL_Hotspot_Ori$end,
                                 logpval=-log10(YL_eQTL_Hotspot_Ori$pvalue),
                                 colvalue=2)
dim(YL_eQTL_Hotspot_OK);head(YL_eQTL_Hotspot_OK);tail(YL_eQTL_Hotspot_OK)   ##  380
YL_eQTL_Hotspot_OK[is.infinite(YL_eQTL_Hotspot_OK$logpval),"logpval"] <- max(YL_eQTL_Hotspot_OK[!is.infinite(YL_eQTL_Hotspot_OK$logpval),"logpval"])


##==================================================================##
##--loading the TWAS results
##==================================================================##
## KF
KF_TWASres_Ori <- read.xlsx(KF_TWASres_Path,cols = c(2:6,10))
dim(KF_TWASres_Ori)   ##174
KF_TWASres_Ori <- KF_TWASres_Ori[grep("random",KF_TWASres_Ori$GeneChr,invert = T),]
rownames(KF_TWASres_Ori) <- NULL
KF_TWASres_Ori$chr <- substr(KF_TWASres_Ori$GeneChr,4,6)
dim(KF_TWASres_Ori);head(KF_TWASres_Ori);tail(KF_TWASres_Ori)   ## 140

KF_TWASres_OK <- data.frame(chr=KF_TWASres_Ori$chr,
                         start=KF_TWASres_Ori$GeneStart,
                         end=KF_TWASres_Ori$GeneEnd,
                         logpval=-log10(KF_TWASres_Ori$P_value),
                         colpvalue=1)
head(KF_TWASres_OK);tail(KF_TWASres_OK);dim(KF_TWASres_OK)   ## 140
summary(KF_TWASres_OK$logpval)


## YL
YL_TWASres_Ori <- read.xlsx(YL_TWASres_Path,cols = c(2:6,10))
dim(YL_TWASres_Ori)   ##181
YL_TWASres_Ori <- YL_TWASres_Ori[grep("random",YL_TWASres_Ori$GeneChr,invert = T),]
rownames(YL_TWASres_Ori) <- NULL
YL_TWASres_Ori$chr <- substr(YL_TWASres_Ori$GeneChr,4,6)
dim(YL_TWASres_Ori);head(YL_TWASres_Ori);tail(YL_TWASres_Ori)   ## 156

YL_TWASres_OK <- data.frame(chr=YL_TWASres_Ori$chr,
                            start=YL_TWASres_Ori$GeneStart,
                            end=YL_TWASres_Ori$GeneEnd,
                            logpval=-log10(YL_TWASres_Ori$P_value),
                            colpvalue=2)
head(YL_TWASres_OK);tail(YL_TWASres_OK);dim(YL_TWASres_OK)   ## 156
summary(YL_TWASres_OK$logpval)


dim(KF_GWASres_OK);head(KF_GWASres_OK)    ## 99
dim(YL_GWASres_OK);head(YL_GWASres_OK)    ## 317

dim(KF_TWASres_OK);head(KF_TWASres_OK)    ## 140
dim(YL_TWASres_OK);head(YL_TWASres_OK)    ## 156

dim(KF_eQTL_Hotspot_OK);head(KF_eQTL_Hotspot_OK)    ## 317
dim(YL_eQTL_Hotspot_OK);head(YL_eQTL_Hotspot_OK)    ## 380


##==================================================================##
##--loading the information of genome
##==================================================================##
GenomeInfo <- read.table(GenomeInfoPath,header = F,sep = "\t")
GenomeInfoOk <- data.frame(chr=substr(GenomeInfo[1:19,1],4,6),start=1,end=GenomeInfo[1:19,2])
df=GenomeInfoOk
df


##==================================================================##
##--Plot the cricle in one panel plot
##==================================================================##
pdf(file = "zRes01.Bn_ASV203_GWAS_eGWAS_TWAS.circle.pdf",width = 8,height = 8)
circos.clear()
##--initialize
circos.genomicInitialize(df)
##--plot region
circos.genomicTrackPlotRegion(df,ylim=c(0,0.1),track.height=0.03,
                              bg.col=c(rep("orange",10),rep("skyblue",9)), ##rand_color(21)
                              cell.padding=c(0.01, 0.5, 0.01, 0.5),
                              track.margin=c(0, 0)) 
# show_col(c("#FF00FF", "#4169E1"))

## show for TWAS results
col_fun = colorRamp2(breaks = c(1, 2), colors = c("#FF00FF", "#4169E1"),transparency = 0)
TWAS.bed_list = list(KF_TWASres_OK, YL_TWASres_OK)
circos.genomicTrack(TWAS.bed_list, stack = TRUE, track.height = 0.12,
                    panel.fun = function(region, value, ...) {
                      i = getI(...)
                      circos.genomicRect(region, value,
                                         border = col_fun(value[[2]]), 
                                         ytop = i + 0.5, ybottom = i - 0.5,
                                         col = col_fun(value[[2]]), ...)
                    })

## show for eQTL hotspots
col_fun = colorRamp2(breaks = c(1, 2), colors = c("#FF00FF", "#4169E1"),transparency = 0)
eQTL.bed_list = list(KF_eQTL_Hotspot_OK, YL_eQTL_Hotspot_OK)
circos.genomicTrack(eQTL.bed_list, stack = TRUE, track.height = 0.12,
                    panel.fun = function(region, value, ...) {
                      i = getI(...)
                      circos.genomicRect(region, value,
                                         border = col_fun(value[[2]]), 
                                         ytop = i + 0.5, ybottom = i - 0.5,
                                         col = col_fun(value[[2]]), ...)
                    })

## show for GWAS 
col_fun = colorRamp2(breaks = c(1, 2), colors = c("#FF00FF", "#4169E1"),transparency = 0)
QTL.bed_list = list(KF_GWASres_OK, YL_GWASres_OK)
circos.genomicTrack(QTL.bed_list, stack = TRUE, track.height = 0.12,
                    panel.fun = function(region, value, ...) {
                      i = getI(...)
                      circos.genomicRect(region, value,
                                         border = col_fun(value[[2]]), 
                                         ytop = i + 0.5, ybottom = i - 0.5,
                                         col = col_fun(value[[2]]), ...)
                    })
dev.off()



