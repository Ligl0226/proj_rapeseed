##==================================================================##
##--Description: eQTL: GeneExp + RNA seq SNPs
##--Highlight: for YL
##--Maintainer: Dr. GuoLiang Li(lig@ipk-gatersleben.de or guoliangli0226@gmail.com)
##--Date: 2023.10.05
##==================================================================##
cat("\n##=======Receive parameters to enter the program=======##\n")
Sys.time()

args <- commandArgs(trailingOnly = TRUE)
cat("##===WorkDir:",args[1],"===##\n")
cat("##===GWASfunction:",args[2],"===##\n")
cat("##===Path_Rapeseed175list:",args[3],"===##\n")
cat("##===Path_YL_GeneExp:",args[4],"===##\n")
cat("##===Path_YL_SampName:",args[5],"===##\n")
cat("##===Path_YL_SNPsPos:",args[6],"===##\n")
cat("##===Path_YL_RNASNPs:",args[7],"===##\n")
cat("##===Path_kinship:",args[8],"===##\n")
cat("##===Memo:",args[9],"===##\n")
cat("##===GeneExpScan_StartNo:",args[10],"===##\n")
cat("##===GeneExpScan_EndNo:",args[11],"===##\n")
cat("##===out_memo:",args[12],"===##\n")

cat("\n##===print(args):\n");print(args)
WorkDir <- args[1]
GWASfunction <- args[2]
Path_Rapeseed175list <- args[3]
Path_YL_GeneExp <- args[4]
Path_YL_SampName <- args[5]
Path_YL_SNPsPos <- args[6]
Path_YL_RNASNPs <- args[7]
Path_kinship <- args[8]
Memo <- args[9]
GeneExpScan_StartNo <- as.numeric(args[10])
GeneExpScan_EndNo <- as.numeric(args[11])
out_memo <- args[12]


##==================================================================##
##---Initialization settings
##==================================================================##
cat("\n##=======Initialization settings=======##\n")
setwd(WorkDir);getwd();## dir()
library(openxlsx)
library(dplyr)
library(stringr)
library(data.table)
library(BGLR)
library(CMplot)
library(qs)
source(GWASfunction)


##==================================================================##
## chromosome name list
##==================================================================##
ACxx <- c("A01","A02", "A03", "A04", "A05", "A06",
          "A07", "A08", "A09", "A10", "C01", "C02",
          "C03", "C04", "C05", "C06", "C07", "C08",
          "C09")
chrACxx <- c("chrA01","chrA02", "chrA03", "chrA04", "chrA05", "chrA06",
             "chrA07", "chrA08", "chrA09", "chrA10", "chrC01", "chrC02",
             "chrC03", "chrC04", "chrC05", "chrC06", "chrC07", "chrC08",
             "chrC09")
chrxx <- paste("chr",sprintf("%02d",1:19),sep = "")
chrx <- paste("chr",c(1:19),sep = "")
ChrName <- data.frame(ACxx,chrACxx,chrxx,chrx)


Sys.time()
##==================================================================##
##--loading the sample name list
##==================================================================##
Rapeseed175list <- read.table(Path_Rapeseed175list,header = T)
head(Rapeseed175list);tail(Rapeseed175list);dim(Rapeseed175list)


Sys.time()
##==================================================================##
##--loading the Gene Expression data
##==================================================================##
##--YL
GeneExp_YL <- read.table(Path_YL_GeneExp,header = T)
cat("dim(GeneExp_YL):",dim(GeneExp_YL),"\n")
if(all.equal(Rapeseed175list$YLid,rownames(GeneExp_YL))){
  rownames(GeneExp_YL) <- Rapeseed175list$RapeID
}
dim(GeneExp_YL);GeneExp_YL[1:5,1:5]  ## 175 17006
sum(is.na(GeneExp_YL)) ##  0


Sys.time()
##==================================================================##
##--Loading the RNA seq SNPs data
##==================================================================##
##--YL
SampleNames <- read.table(file = Path_YL_SampName)[,1]
all.equal(Rapeseed175list$YLid,SampleNames)  ## TRUE
SNPsPos <- fread(file = Path_YL_SNPsPos)
dim(SNPsPos);head(SNPsPos);tail(SNPsPos)
RNAseq_YL_SNPs <- fread(input = Path_YL_RNASNPs,drop = 1,data.table = F)
RNAseq_YL_SNPs[1:10,1:10];dim(RNAseq_YL_SNPs)  ## [1]    175 

if(all.equal(Rapeseed175list$YLid,SampleNames)){
  rownames(RNAseq_YL_SNPs) <- Rapeseed175list$RapeID
}
colnames(RNAseq_YL_SNPs) <- paste(SNPsPos$V1,SNPsPos$V2,sep = "_")
RNAseq_YL_SNPs[1:10,1:8];dim(RNAseq_YL_SNPs)  ## [1]    175 
sum(is.na(RNAseq_YL_SNPs)) ##  0


Sys.time()
##==================================================================##
##--loading the kinship matrix
##==================================================================##
load(Path_kinship,verbose = T)
cat("length(ResKinMat):",length(ResKinMat),"\n")
cat("names(ResKinMat):",names(ResKinMat),"\n")
for(i in 1:length(ResKinMat)){
  cat(names(ResKinMat)[i],"\t",dim(ResKinMat[[i]]),"\n")
}
## all TRUE below.
all.equal(rownames(GeneExp_YL),rownames(RNAseq_YL_SNPs))
all.equal(rownames(GeneExp_YL),rownames(ResKinMat$KinWGS))
all.equal(rownames(GeneExp_YL),rownames(ResKinMat$KinRNAkf))
all.equal(rownames(GeneExp_YL),rownames(ResKinMat$KinRNAyl))


Sys.time()
##==================================================================##
##---Creat the folder to save the results and outlog file
##==================================================================##
cat("\n##===Creat the folder to save the results and outlog file!===##\n")

outlogdir <- paste(WorkDir,"/outlog",sep = "")
if(!dir.exists(outlogdir)){dir.create(outlogdir)}

tempBGLRdir <- paste(WorkDir,"/outlog/",out_memo,"/",sep = "")
if(!dir.exists(tempBGLRdir)){dir.create(tempBGLRdir)}

resultdir <- paste(WorkDir,"/results",sep = "")
if(!dir.exists(resultdir)){dir.create(resultdir)}

res_summary_dir <- paste(WorkDir,"/zRes01.eQTL_Summary",sep = "")
if(!dir.exists(res_summary_dir)){dir.create(res_summary_dir)}

res_plot_dir <- paste(res_summary_dir,"/R001.some_ManhattanQQplot",sep = "")
if(!dir.exists(res_plot_dir)){dir.create(res_plot_dir)}


Sys.time()
##==================================================================##
##--eQTL analysis using the Rcode from Yong(2023)
##==================================================================##
M = RNAseq_YL_SNPs
G <- list()
G$A = ResKinMat$KinRNAyl
all.equal(rownames(M),rownames(G$A))

SigniLevel <- 0.05
anum <- ncol(M)
thd_a <- SigniLevel/anum   ## -log10(thd_a)   ##

for (i in GeneExpScan_StartNo:GeneExpScan_EndNo) {
  # i=13500
  ## set the value for Gene i
  y = as.matrix(GeneExp_YL[,i])
  GeneNamei <- colnames(GeneExp_YL)[i]
  
  ## performing GWAS
  res_p <- GWAS_P3D_MultKinMat(y,M,G,X=NULL,metric_kin="F-infinity",
                               kin_include=c("A"),scan_include=c("A"),
                               tempBGLRdir=paste(tempBGLRdir,GeneNamei,"_",sep = ""))
  res_p_A <- data.frame(SNP=rownames(res_p$A),
                        chr=str_split_fixed(rownames(res_p$A),"_",2)[,1],
                        pos=as.numeric(str_split_fixed(rownames(res_p$A),"_",2)[,2]),
                        res_p$A)
  colnames(res_p_A)[6] <- paste("pval_",GeneNamei,sep = "")
  
  ## save the results
  cat(format(Sys.time(),"%Y-%m-%d %H:%M:%S"),
      sprintf("...Done for %s", sprintf("%05d",i)),GeneNamei,"\n")
  resSavei <- paste(resultdir,"/",Memo,sprintf("%05d",i),"_",GeneNamei,".qs",sep = "")
  qsave(res_p_A,file = resSavei)
  
  ## generate manhattan and QQ plot for each ASV
  if(i %% 100 == 0 | i==1 | i==17006){
    resi <- res_p_A[,c(1,2,3,6)]
    resi.noNA <- resi[!is.na(resi[,4]),]
    resi.noNA.no0 <- resi.noNA[resi.noNA[,4]!=0,]
    dim(resi.noNA.no0);head(resi.noNA.no0)
    
    resi.plot <- resi.noNA.no0
    colnames(resi.plot)[4] <- paste("No",sprintf("%05d",i),"_",GeneNamei,sep = "")
    ## change the chromsome name to A/C
    for(j in 1:length(chrx)){
      ## j=1
      chrname1j <- chrx[j]
      chrname2j <- ACxx[j]
      resi.plot[resi.plot$chr==chrname1j,"chr"] <- chrname2j
    }
    
    ## Manhattan plot for ASV i
    CMplot(resi.plot,type="p",plot.type="m",LOG10=TRUE,threshold=thd_a,file="jpg",
           main = paste(Memo,sprintf("%05d",i),"_",GeneNamei,sep = ""),
           dpi=300,file.output=TRUE,verbose=F,width=14,height=6,
           chr.labels.angle=25)
    fromi=paste("Rectangular-Manhattan.",colnames(resi.plot)[4],".jpg",sep = "")
    toi=paste(res_plot_dir,"/",Memo,sprintf("%03d",i),"_",GeneNamei,"_MHplot.jpg",sep = "")
    file.rename(from = fromi,to = toi)
    
    ## QQ plot for ASV i
    CMplot(resi.plot,plot.type="q",box=FALSE,file="jpg",
           main = paste(Memo,sprintf("%05d",i),"_",GeneNamei,sep = ""),
           dpi=300,conf.int=TRUE,conf.int.col=NULL,
           threshold.col="red",threshold.lty=2,
           file.output=TRUE,verbose=F,width=6,height=6)
    fromi=paste("QQplot.",colnames(resi.plot)[4],".jpg",sep = "")
    toi=paste(res_plot_dir,"/",Memo,sprintf("%03d",i),"_",GeneNamei,"_QQplot.jpg",sep = "")
    file.rename(from = fromi,to = toi)
  }
}
Sys.time()






