##==================================================================##
##--Description: GWAS: SNPs(WGS: 345,289) + ASV (203)
##--Highlight: for KF
##--Maintainer: Dr. GuoLiang Li(lig@ipk-gatersleben.de or guoliangli0226@gmail.com)
##--Date: 2023.10.04
##==================================================================##
cat("\n##=======Receive parameters to enter the program=======##\n")
Sys.time()

args <- commandArgs(trailingOnly = TRUE)
cat("##===WorkDir:",args[1],"===##\n")
cat("##===GWASfunction:",args[2],"===##\n")
cat("##===Path_Rapeseed175list:",args[3],"===##\n")
cat("##===Path_Indivlist:",args[4],"===##\n")
cat("##===Path_MPdata:",args[5],"===##\n")
cat("##===Path_GenoData:",args[6],"===##\n")
cat("##===Path_KF_ASV:",args[7],"===##\n")
cat("##===Path_KinMat:",args[8],"===##\n")
cat("##===ASVscan_StartNo:",args[9],"===##\n")
cat("##===ASVscan_EndNo:",args[10],"===##\n")
cat("##===out_memo:",args[11],"===##\n")
cat("##===Memo:",args[12],"===##\n")

cat("\n##===print(args):\n");print(args)
WorkDir <- args[1]
GWASfunction <- args[2]
Path_Rapeseed175list <- args[3]
Path_Indivlist <- args[4]
Path_MPdata <- args[5]
Path_GenoData <- args[6]
Path_KF_ASV <- args[7]
Path_KinMat <- args[8]
ASVscan_StartNo <- as.numeric(args[9])
ASVscan_EndNo <- as.numeric(args[10])
out_memo <- args[11]
Memo <- args[12]


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
##--loading the WGS SNPs data: Bn project: 175×345,289 SNPs
##==================================================================##
Indivlist <- read.table(Path_Indivlist,header = F)
head(Indivlist);tail(Indivlist);dim(Indivlist)

MPdata <- fread(Path_MPdata,header = F,sep = "\t",data.table = F)
colnames(MPdata) <- c("chr","pos")
head(MPdata);tail(MPdata);dim(MPdata)

GNdata <- fread(Path_GenoData,drop = 1)
GNdata <- as.matrix(GNdata)
rownames(GNdata) <- Indivlist$V1
colnames(GNdata) <- paste(MPdata[,"chr"],MPdata[,"pos"],sep = "_")
dim(GNdata);GNdata[1:5,1:5]   ## [1] 175 345,289


Sys.time()
##==================================================================##
##--loading the ASV data
##==================================================================##
##--KF
ASV_KF <- read.table(Path_KF_ASV,header = T)
cat("dim(ASV_KF):",dim(ASV_KF),"\n")
if(all.equal(Rapeseed175list$KFid,rownames(ASV_KF))){
  rownames(ASV_KF) <- Rapeseed175list$RapeID
}
dim(ASV_KF);ASV_KF[1:5,1:5]  ## 175 203
sum(is.na(ASV_KF)) ##  0


Sys.time()
##==================================================================##
##--loading the kinship matrix based on the WGS SNPs data
##==================================================================##
load(Path_KinMat,verbose = T)
length(ResKinMat);names(ResKinMat)
for(i in 1:length(ResKinMat)){
  cat(names(ResKinMat)[i],"\t",dim(ResKinMat[[i]]),"\n")
}
## for WGS SNPs Relationship Matrix (KinWGS)
KinWGS <- ResKinMat$KinWGS
dim(KinWGS);KinWGS[1:5,1:5]


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

res_summary_dir <- paste(WorkDir,"/zRes01.GWAS_Summary",sep = "")
if(!dir.exists(res_summary_dir)){dir.create(res_summary_dir)}

res_plot_dir <- paste(res_summary_dir,"/Res01.ManhattanQQplot",sep = "")
if(!dir.exists(res_plot_dir)){dir.create(res_plot_dir)}


Sys.time()
##==================================================================##
##--GWAS analysis using the Rcode from Yong(2023)
##==================================================================##
M <- GNdata;dim(M);M[1:5,1:5]
G <- list()
G$A <- KinWGS;dim(G$A);G$A[1:5,1:5]
all.equal(rownames(M),rownames(G$A))

SigniLevel <- 0.05
anum <- ncol(GNdata)
thd_a <- SigniLevel/anum   ## -log10(thd_a)   ## 6.849706

for (i in ASVscan_StartNo:ASVscan_EndNo) {
  ## i=106
  ## set the value for ASV i
  y = as.matrix(ASV_KF[,i])
  ASVnamei <- colnames(ASV_KF)[i]
  
  ## performing GWAS by the function from Yong
  res_p <- GWAS_P3D_MultKinMat(y=y,M=M,G=G,X=NULL,metric_kin="F-infinity",
                               kin_include=c("A"),scan_include=c("A"),
                               tempBGLRdir=paste(tempBGLRdir,ASVnamei,"_",sep = ""))
  res <- data.frame(SNPs=rownames(res_p$A),res_p$A)
  # dim(res);head(res);tail(res)
  
  ## save the results
  cat(format(Sys.time(),"%Y-%m-%d %H:%M:%S"),Memo,sprintf("%03d",i),ASVnamei,"was done!\n")
  resSavei <- paste(resultdir,"/",Memo,sprintf("%03d",i),"_",ASVnamei,".Rdata",sep = "")
  save(res,file = resSavei)
  
  ## generate manhattan and QQ plot for each ASV
  resi <- res
  resi.noNA <- resi[!is.na(resi[,4]),]
  resi.noNA.no0 <- resi.noNA[resi.noNA[,4]!=0,]
  dim(resi.noNA.no0);head(resi.noNA.no0)

  resi.plot <- data.frame(SNPs=resi.noNA.no0$SNPs,
                          chr=str_split_fixed(resi.noNA.no0$SNPs,"_",2)[,1],
                          pos=as.numeric(str_split_fixed(resi.noNA.no0$SNPs,"_",2)[,2]),
                          pval=resi.noNA.no0$P_value)
  colnames(resi.plot)[4] <- paste("No",sprintf("%03d",i),"_",ASVnamei,sep = "")
  ## change the chromsome name to A/C
  for(j in 1:length(chrx)){
    ## j=1
    chrname1j <- chrx[j]
    chrname2j <- ACxx[j]
    resi.plot[resi.plot$chr==chrname1j,"chr"] <- chrname2j
  }

  ## Manhattan plot for ASV i
  CMplot(resi.plot,type="p",plot.type="m",LOG10=TRUE,threshold=thd_a,file="jpg",
         main = paste(Memo,sprintf("%03d",i),"_",ASVnamei,sep = ""),
         dpi=300,file.output=TRUE,verbose=F,width=14,height=6,
         chr.labels.angle=25)
  fromi=paste("Rectangular-Manhattan.",colnames(resi.plot)[4],".jpg",sep = "")
  toi=paste(res_plot_dir,"/",Memo,sprintf("%03d",i),"_",ASVnamei,"_MHplot.jpg",sep = "")
  file.rename(from = fromi,to = toi)
  
  ## QQ plot for ASV i
  CMplot(resi.plot,plot.type="q",box=FALSE,file="jpg",
         main = paste(Memo,sprintf("%03d",i),"_",ASVnamei,sep = ""),
         dpi=300,conf.int=TRUE,conf.int.col=NULL,
         threshold.col="red",threshold.lty=2,
         file.output=TRUE,verbose=F,width=6,height=6)
  fromi=paste("QQplot.",colnames(resi.plot)[4],".jpg",sep = "")
  toi=paste(res_plot_dir,"/",Memo,sprintf("%03d",i),"_",ASVnamei,"_QQplot.jpg",sep = "")
  file.rename(from = fromi,to = toi)
}
Sys.time()




