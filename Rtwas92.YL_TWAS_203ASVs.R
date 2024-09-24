##==================================================================##
##--Description: TWAS: GeneExp (17006) + ASV (203)
##--Highlight: for YL
##--Maintainer: Dr. GuoLiang Li(lig@ipk-gatersleben.de or guoliangli0226@gmail.com)
##--Date: 2023.08.31
##==================================================================##
WorkDir="./"
TWASfunction="TWAS_P3D.R"
Path_Rapeseed175list="../../02.CleanOriData/Rape175order.list.txt"
Path_YL_GeneExp="../../02.CleanOriData/03.GeneExp/GeneExp_YL_Rapeseed175x17006.txt"
Path_YL_ASV="../../02.CleanOriData/04.ASV/ASV_YL_Rapeseed175x203.txt"
Path_KinMat="../../K01.Kinship/zResult02.Kinship_matrix.Rdata"

ASVscan_StartNo <- 1
ASVscan_EndNo <- 203
out_memo="tempbatch"
Memo="TWAS_YL"


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
##--loading the ASV data
##==================================================================##
##--YL
ASV_YL <- read.table(Path_YL_ASV,header = T)
cat("dim(ASV_YL):",dim(ASV_YL),"\n")
if(all.equal(Rapeseed175list$YLid,rownames(ASV_YL))){
  rownames(ASV_YL) <- Rapeseed175list$RapeID
}
dim(ASV_YL);ASV_YL[1:5,1:5]  ## 175 203
sum(is.na(ASV_YL)) ##  0


Sys.time()
##==================================================================##
##--loading the kinship matrix based on the Gene Expression data
##==================================================================##
load(Path_KinMat,verbose = T)
length(ResKinMat);names(ResKinMat)
for(i in 1:length(ResKinMat)){
  cat(names(ResKinMat)[i],"\t",dim(ResKinMat[[i]]),"\n")
}
## for YL Gene Expression Relationship Matrix (YL_GexpRM)
YL_GexpRM <- ResKinMat$KinGeneExpYL
dim(YL_GexpRM);YL_GexpRM[1:5,1:5]


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
##--TWAS analysis using the Rcode from Yong(2023)
##==================================================================##
M <- GeneExp_YL;dim(M)
G <- YL_GexpRM;dim(G)
for (i in ASVscan_StartNo:ASVscan_EndNo) {
  ## i=134
  ## set the value for ASV i
  y = as.matrix(ASV_YL[,i])
  ASVnamei <- colnames(ASV_YL)[i]
  
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

