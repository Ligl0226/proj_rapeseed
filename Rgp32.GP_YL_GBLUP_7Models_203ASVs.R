##==================================================================##
##--Description: GBLUP of YL203ASVs using multi-omics data
##--(1): WGS SNPs;
##--(2): RNAseq SNPs;
##--(3): Gene Expression data;
##--Maintainer: Dr. GuoLiang Li(lig@ipk-gatersleben.de or guoliangli0226@gmail.com)
##--Date: 2023.03.30
##==================================================================##
getwd();dir()

library(openxlsx)
library(dplyr)
library(stringr)
library(data.table)
library(BGLR)
library(caret)

Sys.time()
args <- commandArgs(T)
cat("##===traitNo:",args[1],"===##\n") ## here to index the different ASVs
cat("##===nrep:",args[2],"===##\n")    ## repetition times for cross validation

traitNo=as.numeric(args[1])
nrep=as.numeric(args[2])


##==================================================================##
##--loading the sample name list
##==================================================================##
Rapeseed175list <- read.table("../../02.CleanOriData/Rape175order.list.txt",header = T)
head(Rapeseed175list);tail(Rapeseed175list);dim(Rapeseed175list)


##==================================================================##
##--loading the ASV data of ASV (203)
##==================================================================##
KF_ASVs <- read.table("../../02.CleanOriData/04.ASV/ASV_KF_Rapeseed175x203.txt",header = T)
if(all.equal(rownames(KF_ASVs),Rapeseed175list$KFid)){
  rownames(KF_ASVs) <- Rapeseed175list$RapeID
}
KF_ASVs[1:5,1:5];dim(KF_ASVs)

YL_ASVs <- read.table("../../02.CleanOriData/04.ASV/ASV_YL_Rapeseed175x203.txt",header = T)
if(all.equal(rownames(YL_ASVs),Rapeseed175list$YLid)){
  rownames(YL_ASVs) <- Rapeseed175list$RapeID
}
YL_ASVs[1:5,1:5];dim(YL_ASVs)


##==================================================================##
##--loading the kinship matrix
##==================================================================##
load("../../K01.Kinship/zResult02.Kinship_matrix.Rdata",verbose = T)
length(ResKinMat);names(ResKinMat)
for(i in 1:length(ResKinMat)){
  cat(names(ResKinMat)[i],"\t",dim(ResKinMat[[i]]),"\n")
}
## all TRUE below.
all.equal(rownames(YL_ASVs),rownames(ResKinMat$KinWGS))
all.equal(rownames(YL_ASVs),rownames(ResKinMat$KinRNAkf))
all.equal(rownames(YL_ASVs),rownames(ResKinMat$KinRNAyl))
all.equal(rownames(YL_ASVs),rownames(ResKinMat$KinGeneExpKF))
all.equal(rownames(YL_ASVs),rownames(ResKinMat$KinGeneExpYL))
all.equal(rownames(YL_ASVs),rownames(ResKinMat$KinASVkf))
all.equal(rownames(YL_ASVs),rownames(ResKinMat$KinASVyl))
all.equal(rownames(YL_ASVs),rownames(ResKinMat$Kin_WGS_RNAkf))
all.equal(rownames(YL_ASVs),rownames(ResKinMat$Kin_WGS_RNAyl))


##==================================================================##
##--YL: Prediction based on multi-omics data
##==================================================================##
## R function to calculate the prediction accuracy by BGLR
BGLRpred <- function(y,ETA,folds,saveat){
  nSamp <- length(y)
  ytest <- rep(NA,nSamp)
  
  for(f in 1:length(folds)){
    yna <- y
    yna[folds[[f]]] <- NA
    
    BGLR.fit=BGLR(
      y=yna,ETA = ETA,saveAt = saveat,
      nIter=15000,burnIn=5000,verbose = F
    )
    ytest[folds[[f]]] <- BGLR.fit$yHat[folds[[f]]]
    # print(cor(yna[-folds[[f]]],BGLR.fit$yHat[-folds[[f]]]))
  }
  return(ytest)
}


Sys.time()
YionomeTraits <- colnames(YL_ASVs)    ## 203
Traiti <- YionomeTraits[traitNo]      ## for the trait traitNo
reptimes <- nrep

ETA1 <- list(list(K_wgs = ResKinMat$KinWGS, model = "RKHS"))
ETA2 <- list(list(K_rna = ResKinMat$KinRNAyl, model = "RKHS"))
ETA3 <- list(list(K_geneExp = ResKinMat$KinGeneExpYL, model = "RKHS"))
ETA4 <- list(list(K_wgs = ResKinMat$KinWGS, model = "RKHS"),list(K_rna = ResKinMat$KinRNAyl, model = "RKHS"))
ETA5 <- list(list(K_wgs = ResKinMat$KinWGS, model = "RKHS"),list(K_geneExp = ResKinMat$KinGeneExpYL, model = "RKHS"))
ETA6 <- list(list(K_rna = ResKinMat$KinRNAyl, model = "RKHS"),list(K_geneExp = ResKinMat$KinGeneExpYL, model = "RKHS"))
ETA7 <- list(list(K_wgs = ResKinMat$KinWGS, model = "RKHS"),list(K_rna = ResKinMat$KinRNAyl, model = "RKHS"),list(K_geneExp = ResKinMat$KinGeneExpYL, model = "RKHS"))

resPA <- as.data.frame(matrix(data = NA,nrow = reptimes,ncol = 7))
ModelName <- c("WGS","RNA","Gene","WGS_RNA","WGS_Gene","RNA_Gene","WGS_RNA_Gene")
colnames(resPA) <- paste(Traiti,".PA.",ModelName,sep = "")


tempBGLRdir <- "BGLRtemp"
if(!dir.exists(tempBGLRdir)){dir.create(tempBGLRdir)}
Sys.time()
for(i in 1:reptimes){
  ## i=1
  y <- YL_ASVs[,Traiti]
  nSamp <- length(y)
  folds <- createFolds(1:nSamp, k = 5)
  
  resPA[i,1] <- cor(y,BGLRpred(y = y,ETA = ETA1,folds = folds,saveat = paste("./",tempBGLRdir,"/",Traiti,"_",sep = "")),use = "pairwise.complete.obs")
  resPA[i,2] <- cor(y,BGLRpred(y = y,ETA = ETA2,folds = folds,saveat = paste("./",tempBGLRdir,"/",Traiti,"_",sep = "")),use = "pairwise.complete.obs")
  resPA[i,3] <- cor(y,BGLRpred(y = y,ETA = ETA3,folds = folds,saveat = paste("./",tempBGLRdir,"/",Traiti,"_",sep = "")),use = "pairwise.complete.obs")
  resPA[i,4] <- cor(y,BGLRpred(y = y,ETA = ETA4,folds = folds,saveat = paste("./",tempBGLRdir,"/",Traiti,"_",sep = "")),use = "pairwise.complete.obs")
  resPA[i,5] <- cor(y,BGLRpred(y = y,ETA = ETA5,folds = folds,saveat = paste("./",tempBGLRdir,"/",Traiti,"_",sep = "")),use = "pairwise.complete.obs")
  resPA[i,6] <- cor(y,BGLRpred(y = y,ETA = ETA6,folds = folds,saveat = paste("./",tempBGLRdir,"/",Traiti,"_",sep = "")),use = "pairwise.complete.obs")
  resPA[i,7] <- cor(y,BGLRpred(y = y,ETA = ETA7,folds = folds,saveat = paste("./",tempBGLRdir,"/",Traiti,"_",sep = "")),use = "pairwise.complete.obs")
  cat(i,"/",reptimes,"is completed!\n")
}
Sys.time()
write.xlsx(x = resPA,file = paste("./01.resGP/YL_PA_MultiOmics_7Models_203ASVs_",sprintf("%03d",traitNo),Traiti,".xlsx",sep = ""),rowNames=T)

pdf(file = paste("./01.resGP/YL_PA_MultiOmics_7Models_203ASVs_",sprintf("%03d",traitNo),Traiti,"_boxplot.pdf",sep = ""),
    width = 5,height = 3,pointsize = 6)
  boxplot(resPA,ylim=c(min(resPA)-0.01,max(resPA)),xaxt = "n",
          main="YL_PA_GBLUP_MultiOmics_7Models_203ASVs",
          ylab="Prediction Accuracy")
  text(x = c(1:15),y=min(resPA)-0.01,labels = round(apply(resPA,2,mean),4))
  tick <- seq_along(colnames(resPA))
  axis(1, at = tick, labels = FALSE)
  text(tick-0.1, par("usr")[3]-0.03, colnames(resPA), srt = 20, xpd = TRUE)
dev.off()
Sys.time()


