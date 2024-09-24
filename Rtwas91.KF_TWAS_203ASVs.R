##==================================================================##
##--Description: TWAS: GeneExp (17006) + ASV (203)
##--Highlight: for KF
##--Maintainer: Dr. GuoLiang Li(lig@ipk-gatersleben.de or guoliangli0226@gmail.com)
##--Date: 2023.07.28/31
##==================================================================##
WorkDir="./"
TWASfunction="TWAS_P3D.R"
Path_Rapeseed175list="../../02.CleanOriData/Rape175order.list.txt"
Path_KF_GeneExp="../../02.CleanOriData/03.GeneExp/GeneExp_KF_Rapeseed175x17006.txt"
Path_KF_ASV="../../02.CleanOriData/04.ASV/ASV_KF_Rapeseed175x203.txt"
Path_KinMat="../../K01.Kinship/zResult02.Kinship_matrix.Rdata"

ASVscan_StartNo <- 1
ASVscan_EndNo <- 203
out_memo="tempbatch"
Memo="TWAS_KF"


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
source(TWASfunction)


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
##--KF
GeneExp_KF <- read.table(Path_KF_GeneExp,header = T)
cat("dim(GeneExp_KF):",dim(GeneExp_KF),"\n")
if(all.equal(Rapeseed175list$KFid,rownames(GeneExp_KF))){
  rownames(GeneExp_KF) <- Rapeseed175list$RapeID
}
dim(GeneExp_KF);GeneExp_KF[1:5,1:5]  ## 175 17006
sum(is.na(GeneExp_KF)) ##  0


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
##--loading the kinship matrix based on the Gene Expression data
##==================================================================##
load(Path_KinMat,verbose = T)
length(ResKinMat);names(ResKinMat)
for(i in 1:length(ResKinMat)){
  cat(names(ResKinMat)[i],"\t",dim(ResKinMat[[i]]),"\n")
}
## for KF Gene Expression Relationship Matrix (KF_GexpRM)
KF_GexpRM <- ResKinMat$KinGeneExpKF
dim(KF_GexpRM);KF_GexpRM[1:5,1:5]


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


Sys.time()
##==================================================================##
##--TWAS analysis
##==================================================================##
M <- GeneExp_KF;dim(M)
G <- KF_GexpRM;dim(G)
for (i in ASVscan_StartNo:ASVscan_EndNo) {
  ## i=134
  ## set the value for ASV i
  y = as.matrix(ASV_KF[,i])
  ASVnamei <- colnames(ASV_KF)[i]
  
  ## performing GWAS
  res_p <- TWAS_P3D(y,M,G,X=NULL,tempBGLRdir = tempBGLRdir)
  res <- data.frame(GeneName=rownames(res_p),res_p)
  # dim(res);head(res);tail(res)
  
  ## save the results
  cat(format(Sys.time(),"%Y-%m-%d %H:%M:%S"),sprintf("...Done for No %s", sprintf("%03d",i)),ASVnamei,"\n")
  resSavei <- paste(resultdir,"/",Memo,"_",sprintf("%03d",i),"_",ASVnamei,".txt",sep = "")
  write.table(x = res,file = resSavei,sep = "\t",row.names = F,quote = F)
}
Sys.time()

