##==================================================================##
##--Description: summary the GBLUP results for KF GP 7Models for 203ASVs
##--all possible combinations using multi-omics data
##--(1) WGS SNPs
##--(2) RNA SNPs
##--(3) Gene Expression
##--Maintainer: GuoLiang Li(lig@ipk-gatersleben.de)
##--Date: 2023.04.03
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
##--loading the ASVs data of KF
##==================================================================##
KF_ASVs <- read.table("../../02.CleanOriData/04.ASV/ASV_KF_Rapeseed175x203.txt",header = T)
if(all.equal(rownames(KF_ASVs),Rapeseed175list$KFid)){
  rownames(KF_ASVs) <- Rapeseed175list$RapeID
}
KF_ASVs[1:5,1:5];dim(KF_ASVs)


##==================================================================##
##--loading the GP results (100 times 5-folds CV) of 7Models for 203 ASVs
##==================================================================##
ASVs_Traits <- colnames(KF_ASVs)  ## 203
Models7 <- c("WGS","RNA","Gene","WGS_RNA","WGS_Gene","RNA_Gene","WGS_RNA_Gene")
GPmean <- as.data.frame(matrix(data = NA,nrow = length(ASVs_Traits),ncol = 7))
rownames(GPmean) <- ASVs_Traits
colnames(GPmean) <- paste("PAmean",Models7,sep = ".")
GPsd <- GPmean
colnames(GPsd) <- paste("PAsd",Models7,sep = ".")

for(i in 1:length(ASVs_Traits)){
  ## i=1
  Traiti <- ASVs_Traits[i]      ## for the trait i
  if(i<=50){
    resGPi <- read.xlsx(xlsxFile = paste("./01.resGP/KF_PA_MultiOmics_7Models_203ASVs_",sprintf("%02d",i),Traiti,".xlsx",sep = ""),cols = c(2:8))
    GPmean[i,] <- apply(X = resGPi,MARGIN = 2,FUN = mean)
    GPsd[i,] <- apply(X = resGPi,MARGIN = 2,FUN = sd)
  }else{
    resGPi <- read.xlsx(xlsxFile = paste("./01.resGP/KF_PA_MultiOmics_7Models_203ASVs_",sprintf("%03d",i),Traiti,".xlsx",sep = ""),cols = c(2:8))
    GPmean[i,] <- apply(X = resGPi,MARGIN = 2,FUN = mean)
    GPsd[i,] <- apply(X = resGPi,MARGIN = 2,FUN = sd)
  }
}
resGPmeansd <- cbind(GPmean,GPsd)
dim(resGPmeansd)  ## [1] 203  14
write.xlsx(x = resGPmeansd,"./02.SummaryGP/zResult01.GP_KF_GBLUP_MultiOmics_7Models_203ASVs_Summary.xlsx",rowNames=T)


##==================================================##
##--KF: which.max()
##==================================================##
resGPmean_3Models <- resGPmeansd[,1:3]
head(resGPmean_3Models)
resGPmean_3Models$max.index <- apply(resGPmean_3Models,1,function (x) names(which.max(x)))
head(resGPmean_3Models)
as.data.frame(table(resGPmean_3Models$max.index))
resGPmean_3Models %>% filter(max.index=="PAmean.WGS")
resGPmean_3Models %>% filter(max.index=="PAmean.RNA")



resGPmean_7Models <- resGPmeansd[,1:7]
head(resGPmean_7Models)
resGPmean_7Models$max.index <- apply(resGPmean_7Models,1,function (x) names(which.max(x)))
head(resGPmean_7Models)
as.data.frame(table(resGPmean_7Models$max.index))
resGPmean_7Models %>% filter(max.index=="PAmean.RNA")
resGPmean_7Models %>% filter(max.index=="PAmean.WGS")
resGPmean_7Models %>% filter(max.index=="PAmean.WGS_RNA")







