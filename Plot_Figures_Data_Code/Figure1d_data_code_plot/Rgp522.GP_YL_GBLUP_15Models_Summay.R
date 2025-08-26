##==================================================================##
##--Description: summary the GBLUP results for YL GP 15Models
##--all possible combinations using multi-omics data
##--(1) WGS SNPs
##--(2) RNA SNPs
##--(3) Gene Expression
##--(4) ASV data;
##--Maintainer: GuoLiang Li(lig@ipk-gatersleben.de)
##--Date: 2023.03.11
##==================================================================##
library(openxlsx)
library(dplyr)
library(stringr)
library(data.table)


##==================================================================##
##--loading the sample name list
##==================================================================##
Rapeseed175list <- read.table("../../02.CleanOriData/Rape175order.list.txt",header = T)
head(Rapeseed175list);tail(Rapeseed175list);dim(Rapeseed175list)


##==================================================================##
##--loading the phenotype data
##==================================================================##
KF_ionome <- read.table("../../02.CleanOriData/05.ionome/ionome_KF_Rapeseed175x39.txt",header = T)
if(all.equal(rownames(KF_ionome),Rapeseed175list$KFid)){
  rownames(KF_ionome) <- Rapeseed175list$RapeID
}
KF_ionome[1:5,1:5];dim(KF_ionome)

YL_ionome <- read.table("../../02.CleanOriData/05.ionome/ionome_YL_Rapeseed175x39.txt",header = T)
if(all.equal(rownames(YL_ionome),Rapeseed175list$YLid)){
  rownames(YL_ionome) <- Rapeseed175list$RapeID
}
YL_ionome[1:5,1:5];dim(YL_ionome)


##==================================================================##
##--loading the GP results (100 times 5-folds CV) of 15Models
##==================================================================##
YionomeTraits <- colnames(YL_ionome)  ## 39
Models15 <- c("WGS","RNA","Gene","ASV",
              "WGS_RNA","WGS_Gene","WGS_ASV","RNA_Gene","RNA_ASV","Gene_ASV",
              "WGS_RNA_Gene","WGS_RNA_ASV","WGS_Gene_ASV","RNA_Gene_ASV",
              "WGS_RNA_Gene_ASV")
GPmean <- as.data.frame(matrix(data = NA,nrow = length(YionomeTraits),ncol = 15))
rownames(GPmean) <- YionomeTraits
colnames(GPmean) <- paste("PAmean",Models15,sep = ".")

GPsd <- GPmean
colnames(GPsd) <- paste("PAsd",Models15,sep = ".")

for(i in 1:length(YionomeTraits)){
  ## i=1
  Traiti <- YionomeTraits[i]      ## for the trait i
  resGPi <- read.xlsx(xlsxFile = paste("./01.resGP/YL_PA_MultiOmics_15Models_",sprintf("%02d",i),Traiti,".xlsx",sep = ""),cols = c(2:16))
  GPmean[i,] <- apply(X = resGPi,MARGIN = 2,FUN = mean)
  GPsd[i,] <- apply(X = resGPi,MARGIN = 2,FUN = sd)
}
resGPmeansd <- cbind(GPmean,GPsd)
dim(resGPmeansd)  ## [1] 39 30
write.xlsx(x = resGPmeansd,"./02.SummaryGP/zResult01.GP_YL_GBLUP_MultiOmics_15Models_Summary.xlsx",rowNames=T)

