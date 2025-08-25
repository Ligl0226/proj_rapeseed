##==================================================================##
##--Description: eQTL results summary for KF
##--Highlight: the eQTL results was based on the RNAseq SNPs after MAF>=0.05 & GenoHet <= 0.95 filtering
##--Highlight: plot the pval (only for eSNPs not possible for all SNPs)
##--Maintainer: GuoLiang Li(lig@ipk-gatersleben.de)
##--Date: 2023.10.09
##==================================================================##
library(openxlsx)
library(stringr)
library(data.table)
library(dplyr)
library(CMplot)
library(qs)
library(circlize)
library(RColorBrewer)
library(ggplot2)
library(reshape2)
library(colorRamp2)

## the function give different color based on the value
colorsplit <- function(x=datavector,
                       quant.probs=seq(0,1,0.1),
                       colfrom="yellow",
                       colto="red"){
  colxx <- colorRampPalette(colors = c(colfrom,colto))(length(quant.probs)-1)
  quantilexx <- quantile(x,probs = quant.probs)
  ycol <- x
  for(i in 1:(length(quant.probs)-1))
    ycol[which(x>=quantilexx[i]&x<=quantilexx[i+1])] <- colxx[i]
  return(ycol)
}


Sys.time()
##==================================================================##
##--loading gff3.mRNA file: 
##==================================================================##
Bngff3.mRNA <- read.table("../../GenomeInfo/v4.1_20140819/Brassica_napus.annotation_v5.gff3.mRNA.txt",header = T)
rownames(Bngff3.mRNA) <- Bngff3.mRNA$ID
head(Bngff3.mRNA);tail(Bngff3.mRNA);dim(Bngff3.mRNA)


Sys.time()
##==================================================================##
##--loading the sample name list
##==================================================================##
Path_Rapeseed175list <- "../../02.CleanOriData/Rape175order.list.txt"
Rapeseed175list <- read.table(Path_Rapeseed175list,header = T)
head(Rapeseed175list);tail(Rapeseed175list);dim(Rapeseed175list)


Sys.time()
##==================================================================##
##--loading the Gene Expression data
##==================================================================##
##--KF
Path_KF_GeneExp <- "../../02.CleanOriData/03.GeneExp/GeneExp_KF_Rapeseed175x17006.txt"
GeneExp_KF <- read.table(Path_KF_GeneExp,header = T)
cat("dim(GeneExp_KF):",dim(GeneExp_KF),"\n")
if(all.equal(Rapeseed175list$KFid,rownames(GeneExp_KF))){
  rownames(GeneExp_KF) <- Rapeseed175list$RapeID
}
dim(GeneExp_KF);GeneExp_KF[1:5,1:5]  ## 175 17006
sum(is.na(GeneExp_KF)) ##  0
GeneList <- colnames(GeneExp_KF)
length(GeneList)

## Is the gene in GeneList all included in Bngff3.mRNA$ID? Yes!
length(intersect(GeneList,Bngff3.mRNA$ID)) ## 17006


##==================================================================##
##--Bonferroni correction threshold
##==================================================================##
SigniLevel <- 0.05
anum <- 239172            ## KF
thd_a <- SigniLevel/anum  ## -log10(thd_a)   ## 6.67974


##==================================================================##
##--loading the significant eSNPs for all Genes
##==================================================================##
res.gwas.geneExp.sig <- read.table("./zRes01.eQTL_Summary/zRes71.KF_eQTL_GWAS_SigSNPs_Summary.txt",header = T)
dim(res.gwas.geneExp.sig)
head(res.gwas.geneExp.sig)
tail(res.gwas.geneExp.sig)
as.data.frame(table(res.gwas.geneExp.sig$GeneChr))
as.data.frame(table(res.gwas.geneExp.sig$chr))

res.gwas.geneExp.sig.chr <- res.gwas.geneExp.sig[grep(pattern = "_random",x = res.gwas.geneExp.sig$GeneChr,invert = T),]
rownames(res.gwas.geneExp.sig.chr) <- NULL
dim(res.gwas.geneExp.sig.chr)   ## 274753      9
head(res.gwas.geneExp.sig.chr)
tail(res.gwas.geneExp.sig.chr)
as.data.frame(table(res.gwas.geneExp.sig.chr$GeneChr))
as.data.frame(table(res.gwas.geneExp.sig.chr$chr))


##==================================================================##
##--loading the reference genome information
##==================================================================##
BnaRefInfo <- read.table("./Brassica_napus_v4.1.chromosomes.fa.fai",header = F)
## only work with known chr
chr.length <- BnaRefInfo[c(1:19),]
# chr.length$V1 <- str_split_fixed(chr.length$V1,"_",2)[,1]
chr.length$cumsumlength <- cumsum(chr.length$V2)
chr.length$prolength <- c(0,chr.length$cumsumlength[1:18])
chr.length$V6 <- paste("chr",c(1:19),sep = "")
chr.length$V7 <- rowMeans(chr.length[,c("cumsumlength","prolength")])
chr.length$V8 <- substr(x = chr.length$V1,start = 4,stop = 6)
chr.bounda <- c(0,chr.length$cumsumlength)


##==================================================================##
##--plot the plot
##==================================================================##
## 0. get the clean data to plot
plot.data <- data.frame(GeneName=res.gwas.geneExp.sig.chr$GeneName,
                        GeneChr=res.gwas.geneExp.sig.chr$GeneChr,
                        GenoPos=rowMeans(res.gwas.geneExp.sig.chr[,c("GeneStart","GeneEnd")]),
                        GenoPosPlot=NA,
                        SNPsChr=res.gwas.geneExp.sig.chr$chr,
                        SNPsPos=res.gwas.geneExp.sig.chr$pos,
                        SNPsPosPlot=NA,
                        pval=res.gwas.geneExp.sig.chr$pval,
                        pval_log10=-log10(res.gwas.geneExp.sig.chr$pval))
plot.data$pval_color <- colorsplit(x = plot.data$pval_log10,quant.probs=seq(0,1,0.1),colfrom="yellow",colto="red")
head(plot.data);tail(plot.data);dim(plot.data)
as.data.frame(table(plot.data$SNPsChr))
as.data.frame(table(plot.data$pval_color))


## give the position (x-asix and y-asix) to plot
for(m in 1:nrow(chr.length)){
  chrm <- chr.length[m,"V1"]
  plot.data[plot.data$GeneChr==chrm,"GenoPosPlot"]=plot.data[plot.data$GeneChr==chrm,"GenoPos"]+chr.length[m,"prolength"]
}

for(n in 1:nrow(chr.length)){
  chrn <- chr.length[n,"V6"]
  plot.data[plot.data$SNPsChr==chrn,"SNPsPosPlot"]=plot.data[plot.data$SNPsChr==chrn,"SNPsPos"]+chr.length[n,"prolength"]
}

jpeg("./zRes02.eQTL_plot/zRes71.KF_eQTL_GWAS_SigSNPs_Gene_2Dplot.jpg",
     width = 800,height = 800)
## 1. plot the frame of chr
maintitle="KF_eSNP_Gene"
plot(chr.bounda,chr.bounda,bty="o",
     type="n",xaxt="n",yaxt="n",
     main = maintitle,cex.lab=1.2,cex.main=1.5,
     xlab="eSNPs positions",ylab="Gene positions")
Noo=0
for(i in 1:(length(chr.bounda)-1)){
  ## i=1
  xfrom=chr.bounda[i:(i+1)]
  for(j in 1:(length(chr.bounda)-1)){
    ## j=2
    yfrom=chr.bounda[j:(j+1)]
    Noo=Noo+1
    if(Noo%%2==0) {col0="grey90"} else {col0="grey80"}
    rect(xfrom[1],yfrom[1],xfrom[2],yfrom[2],col=col0,border="white")
  }
  if(length(chr.bounda)%%2!=0){Noo=Noo+1}
}
axis(side = 1,at = chr.length$V7,labels = chr.length$V8)
axis(side = 2,at = chr.length$V7,labels = chr.length$V8)

## 2. plot the point
x=plot.data$SNPsPosPlot ##[1:10000]
y=plot.data$GenoPosPlot ##[1:10000]
z=plot.data$pval_color ##[1:10000])
points(x,y,type="p",cex=0.5, pch=19,col=z)
dev.off()


pdf("./zRes02.eQTL_plot/zRes71.KF_eQTL_GWAS_SigSNPs_Gene_2Dplot.pdf",
     width = 12,height = 12)
## 1. plot the frame of chr
maintitle="KF_eSNP_Gene"

plot(chr.bounda,chr.bounda,bty="o",
     type="n",xaxt="n",yaxt="n",
     main = maintitle,cex.lab=1.2,cex.main=1.5,
     xlab="eSNPs positions",ylab="Gene positions")
Noo=0
for(i in 1:(length(chr.bounda)-1)){
  ## i=1
  xfrom=chr.bounda[i:(i+1)]
  for(j in 1:(length(chr.bounda)-1)){
    ## j=2
    yfrom=chr.bounda[j:(j+1)]
    Noo=Noo+1
    if(Noo%%2==0) {col0="grey90"} else {col0="grey80"}
    rect(xfrom[1],yfrom[1],xfrom[2],yfrom[2],col=col0,border="white")
  }
  if(length(chr.bounda)%%2!=0){Noo=Noo+1}
}
axis(side = 1,at = chr.length$V7,labels = chr.length$V8)
axis(side = 2,at = chr.length$V7,labels = chr.length$V8)

## 2. plot the point
x=plot.data$SNPsPosPlot ##[1:10000]
y=plot.data$GenoPosPlot ##[1:10000]
z=plot.data$pval_color ##[1:10000])
points(x,y,type="p",cex=0.5, pch=19,col=z)
dev.off()




