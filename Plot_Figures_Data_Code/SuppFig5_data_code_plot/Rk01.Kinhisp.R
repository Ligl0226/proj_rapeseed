##==================================================================##
##--Description: calculate the kinship matrix for each kind of data
##--(1): WGS file:
##--(2): RNAseq file:
##--(3): Gene expression file:
##--(4): ASV file:
##--Maintainer: GuoLiang Li(lig@ipk-gatersleben.de)
##--Date: 2023.03.06
##==================================================================##
library(openxlsx)
library(dplyr)
library(stringr)
library(data.table)
library(VennDiagram)


##==================================================================##
##--function to calculate the kinship matrix
##==================================================================##
# loading omic data
# rows: genotypes
# colums: markers or features
KinMatGen <- function(M,is.SNPs=TRUE){
  M <- as.matrix(M)
  if(is.SNPs){
    M_scale <- M-1
    Kin_M <- M_scale %*% t(M_scale)
    Kin_M <- Kin_M/mean(diag(Kin_M))
  }else{
    M_scale <- scale(M,center = TRUE,scale = TRUE)
    Kin_M <- M_scale %*% t(M_scale)
    Kin_M <- Kin_M/mean(diag(Kin_M))
  }
  return(Kin_M)
}


##==================================================================##
##--loading the sample name list
##==================================================================##
Rapeseed175list <- read.table("../02.CleanOriData/Rape175order.list.txt",header = T)
head(Rapeseed175list);tail(Rapeseed175list);dim(Rapeseed175list)


ResKinMat <- list() ## save all kinds of kinship maitrix based diff data
##==================================================================##
##--(1): WGS SNPs data
##==================================================================##
Sys.time()
SampleNames <- read.table(file = "../02.CleanOriData/01.WGS/WGS_Rapeseed175_biallelic_bI_012.012.indv")[,1]
all.equal(Rapeseed175list$RapeID,SampleNames)  ## TRUE
SNPsPos <- fread(file = "../02.CleanOriData/01.WGS/WGS_Rapeseed175_biallelic_bI_012.012.pos")
dim(SNPsPos);head(SNPsPos);tail(SNPsPos)
WGS_SNPs <- fread(input = "../02.CleanOriData/01.WGS/WGS_Rapeseed175_biallelic_bI_012.012",drop = 1,data.table = F)
WGS_SNPs[1:10,1:10];dim(WGS_SNPs)  ## [1]    175 637823

rownames(WGS_SNPs) <- SampleNames
colnames(WGS_SNPs) <- paste(SNPsPos$V1,SNPsPos$V2,sep = "_")
WGS_SNPs[1:10,1:8];dim(WGS_SNPs)  ## [1]    175 637823

sum(is.na(WGS_SNPs)) ##  0
ResKinMat$KinWGS <- KinMatGen(WGS_SNPs,is.SNPs = TRUE)
dim(ResKinMat$KinWGS);ResKinMat$KinWGS[1:5,1:5]


##==================================================================##
##--(2): RNAseq SNPs data
##==================================================================##
##--KF
SampleNames <- read.table(file = "../02.CleanOriData/02.RNAseq/KF/RNAseq_KF_Rapeseed175_biallelic_bI_012.012.indv")[,1]
all.equal(Rapeseed175list$KFid,SampleNames)  ## TRUE
SNPsPos <- fread(file = "../02.CleanOriData/02.RNAseq/KF/RNAseq_KF_Rapeseed175_biallelic_bI_012.012.pos")
dim(SNPsPos);head(SNPsPos);tail(SNPsPos)
RNAseq_KF_SNPs <- fread(input = "../02.CleanOriData/02.RNAseq/KF/RNAseq_KF_Rapeseed175_biallelic_bI_012.012",drop = 1,data.table = F)
RNAseq_KF_SNPs[1:10,1:10];dim(RNAseq_KF_SNPs)  ## [1]    175 241558

if(all.equal(Rapeseed175list$KFid,SampleNames)){
  rownames(RNAseq_KF_SNPs) <- Rapeseed175list$RapeID
}
colnames(RNAseq_KF_SNPs) <- paste(SNPsPos$V1,SNPsPos$V2,sep = "_")
RNAseq_KF_SNPs[1:10,1:8];dim(RNAseq_KF_SNPs)  ## [1]    175 241558

sum(is.na(RNAseq_KF_SNPs)) ##  0
ResKinMat$KinRNAkf <- KinMatGen(RNAseq_KF_SNPs,is.SNPs = TRUE)
dim(ResKinMat$KinRNAkf);ResKinMat$KinRNAkf[1:5,1:5]


##--YL
SampleNames <- read.table(file = "../02.CleanOriData/02.RNAseq/YL/RNAseq_YL_Rapeseed175_biallelic_bI_012.012.indv")[,1]
all.equal(Rapeseed175list$YLid,SampleNames)  ## TRUE
SNPsPos <- fread(file = "../02.CleanOriData/02.RNAseq/YL/RNAseq_YL_Rapeseed175_biallelic_bI_012.012.pos")
dim(SNPsPos);head(SNPsPos);tail(SNPsPos)
RNAseq_YL_SNPs <- fread(input = "../02.CleanOriData/02.RNAseq/YL/RNAseq_YL_Rapeseed175_biallelic_bI_012.012",drop = 1,data.table = F)
RNAseq_YL_SNPs[1:10,1:10];dim(RNAseq_YL_SNPs)  ## [1]    175 295925

if(all.equal(Rapeseed175list$YLid,SampleNames)){
  rownames(RNAseq_YL_SNPs) <- Rapeseed175list$RapeID
}
colnames(RNAseq_YL_SNPs) <- paste(SNPsPos$V1,SNPsPos$V2,sep = "_")
RNAseq_YL_SNPs[1:10,1:8];dim(RNAseq_YL_SNPs)  ## [1]    175 295925

sum(is.na(RNAseq_YL_SNPs)) ##  0
ResKinMat$KinRNAyl <- KinMatGen(RNAseq_YL_SNPs,is.SNPs = TRUE)
dim(ResKinMat$KinRNAyl);ResKinMat$KinRNAyl[1:5,1:5]


##==================================================================##
##--(3): Gene Expression data
##==================================================================##
##--KF
GeneExp_KF <- read.table("../02.CleanOriData/03.GeneExp/GeneExp_KF_Rapeseed175x17006.txt",header = T)
if(all.equal(Rapeseed175list$KFid,rownames(GeneExp_KF))){
  rownames(GeneExp_KF) <- Rapeseed175list$RapeID
}
dim(GeneExp_KF);GeneExp_KF[1:5,1:5]  ## 175 17006

sum(is.na(GeneExp_KF)) ##  0
ResKinMat$KinGeneExpKF <- KinMatGen(GeneExp_KF,is.SNPs = FALSE)
dim(ResKinMat$KinGeneExpKF);ResKinMat$KinGeneExpKF[1:5,1:5]


##--YL
GeneExp_YL <- read.table("../02.CleanOriData/03.GeneExp/GeneExp_YL_Rapeseed175x17006.txt",header = T)
if(all.equal(Rapeseed175list$YLid,rownames(GeneExp_YL))){
  rownames(GeneExp_YL) <- Rapeseed175list$RapeID
}
dim(GeneExp_YL);GeneExp_YL[1:5,1:5]  ## 175 17006

sum(is.na(GeneExp_YL)) ##  0
ResKinMat$KinGeneExpYL <- KinMatGen(GeneExp_YL,is.SNPs = FALSE)
dim(ResKinMat$KinGeneExpYL);ResKinMat$KinGeneExpYL[1:5,1:5]


##==================================================================##
##--(4): ASV data
##==================================================================##
##--KF
ASV_KF <- read.table("../02.CleanOriData/04.ASV/ASV_KF_Rapeseed175x203.txt",header = T)
if(all.equal(Rapeseed175list$KFid,rownames(ASV_KF))){
  rownames(ASV_KF) <- Rapeseed175list$RapeID
}
dim(ASV_KF);ASV_KF[1:5,1:5]  ## 175 203

sum(is.na(ASV_KF)) ##  0
ResKinMat$KinASVkf <- KinMatGen(ASV_KF,is.SNPs = FALSE)
dim(ResKinMat$KinASVkf);ResKinMat$KinASVkf[1:5,1:5]


##--YL
ASV_YL <- read.table("../02.CleanOriData/04.ASV/ASV_YL_Rapeseed175x203.txt",header = T)
if(all.equal(Rapeseed175list$YLid,rownames(ASV_YL))){
  rownames(ASV_YL) <- Rapeseed175list$RapeID
}
dim(ASV_YL);ASV_YL[1:5,1:5]  ## 175 203

sum(is.na(ASV_YL)) ##  0
ResKinMat$KinASVyl <- KinMatGen(ASV_YL,is.SNPs = FALSE)
dim(ResKinMat$KinASVyl);ResKinMat$KinASVyl[1:5,1:5]



##==================================================================##
##--check the common SNPs between WGS and RNAseq_KF and RNAseq_YL
##==================================================================##
length(colnames(WGS_SNPs))
length(colnames(RNAseq_KF_SNPs))
length(colnames(RNAseq_YL_SNPs))

length(intersect(colnames(WGS_SNPs),colnames(RNAseq_KF_SNPs)))  ## 7614
length(intersect(colnames(WGS_SNPs),colnames(RNAseq_YL_SNPs)))  ## 9074
length(intersect(colnames(RNAseq_KF_SNPs),colnames(RNAseq_YL_SNPs)))  ## 219225


SNPs_list <- list(WGS_SNPs = colnames(WGS_SNPs),
                  RNAseq_KF_SNPs = colnames(RNAseq_KF_SNPs),
                  RNAseq_YL_SNPs = colnames(RNAseq_YL_SNPs))
venn.diagram(SNPs_list, filename = "zResult01.SNPs_WGS_RNAseq.png",
             disable.logging = TRUE,
             height = 960 , width = 960 ,
             rotation.degree = 30,
             resolution = 300,compression = "lzw",lwd = 1,
             imagetype = "png",fill = c("blue", "red", "green"),alpha = 0.3,
             cat.col = c("black"),cat.cex = 0.5, scaled = TRUE)

overlap_WGS_RNAkf <- intersect(colnames(WGS_SNPs),colnames(RNAseq_KF_SNPs))
overlap_WGS_RNAkf_RNAyl <- intersect(overlap_WGS_RNAkf,colnames(RNAseq_YL_SNPs))
length(overlap_WGS_RNAkf_RNAyl)  ## [1] 7093
for(i in 1:100){
  random_SNPs <- sample(overlap_WGS_RNAkf_RNAyl,1)
  a <- sum(WGS_SNPs[,random_SNPs]==RNAseq_KF_SNPs[,random_SNPs])/175
  b <- sum(WGS_SNPs[,random_SNPs]==RNAseq_YL_SNPs[,random_SNPs])/175
  c <- sum(RNAseq_KF_SNPs[,random_SNPs]==RNAseq_YL_SNPs[,random_SNPs])/175
  cat(a,"\t",b,"\t",c,"\n")
}



##==================================================================##
##--(5) merge the WGS SNPs and RNA SNPs by cbind()
##--since the overlap SNPs has different genotype between WGS and RNA
##==================================================================##
##--KF
WGS_SNPs[1:10,1:10];dim(WGS_SNPs)  ## [1]    175 637823
RNAseq_KF_SNPs[1:10,1:10];dim(RNAseq_KF_SNPs)  ## [1]    175 241558
if(all.equal(rownames(WGS_SNPs),rownames(RNAseq_KF_SNPs))){
  WGS_RNAseq_KF_SNPs <- cbind(WGS_SNPs,RNAseq_KF_SNPs)
}
dim(WGS_RNAseq_KF_SNPs)  ##  [1]    175 879381

sum(is.na(WGS_RNAseq_KF_SNPs)) ##  0
ResKinMat$Kin_WGS_RNAkf <- KinMatGen(WGS_RNAseq_KF_SNPs,is.SNPs = TRUE)
dim(ResKinMat$Kin_WGS_RNAkf);ResKinMat$Kin_WGS_RNAkf[1:5,1:5]


##--YL
WGS_SNPs[1:10,1:10];dim(WGS_SNPs)  ## [1]    175 637823
RNAseq_YL_SNPs[1:10,1:10];dim(RNAseq_YL_SNPs)  ## [1]    175 295925
if(all.equal(rownames(WGS_SNPs),rownames(RNAseq_YL_SNPs))){
  WGS_RNAseq_YL_SNPs <- cbind(WGS_SNPs,RNAseq_YL_SNPs)
}
dim(WGS_RNAseq_YL_SNPs)  ##  [1]    175 933748

sum(is.na(WGS_RNAseq_YL_SNPs)) ##  0
ResKinMat$Kin_WGS_RNAyl <- KinMatGen(WGS_RNAseq_YL_SNPs,is.SNPs = TRUE)
dim(ResKinMat$Kin_WGS_RNAyl);ResKinMat$Kin_WGS_RNAyl[1:5,1:5]


##==================================================================##
##--save the kinship matrix
##==================================================================##
length(ResKinMat);names(ResKinMat)
save(ResKinMat,file = "zResult02.Kinship_matrix.Rdata")
# load("zResult02.Kinship_matrix.Rdata",verbose = T)


##==================================================================##
##--Calculating the correlation between kinship matrix ()
##==================================================================##
corr2 <- function(matrix1,matrix2){
  matrix1 <- as.matrix(matrix1)
  matrix2 <- as.matrix(matrix2)
  corr2 <- sum((matrix1-mean(matrix1))*(matrix2-mean(matrix2)))/(sum((matrix1-mean(matrix1))^2)*sum((matrix2-mean(matrix2))^2))^0.5
  return(corr2)
}
kinshipMatrix1 <- list(KinWGS=ResKinMat$KinWGS,
                       KinRNAkf=ResKinMat$KinRNAkf,
                       KinRNAyl=ResKinMat$KinRNAyl,
                       Kin_WGS_RNAkf=ResKinMat$Kin_WGS_RNAkf,
                       Kin_WGS_RNAyl=ResKinMat$Kin_WGS_RNAyl,
                       KinGeneExpKF=ResKinMat$KinGeneExpKF,
                       KinGeneExpYL=ResKinMat$KinGeneExpYL,
                       KinASVkf=ResKinMat$KinASVkf,
                       KinASVyl=ResKinMat$KinASVyl)
n <- length(kinshipMatrix1)
K <- names(kinshipMatrix1)
Results_corr1 <- array(NA, dim = c(n,n),dimnames=list(K,K))
for(i in K){
  for(j in K){
    Results_corr1[i,j] <- corr2(kinshipMatrix1[[i]],kinshipMatrix1[[j]])
  }
}
cat("\nResults_corr1:\n");Results_corr1
write.xlsx(x = as.data.frame(Results_corr1),"zResult03.Kinship_corr.xlsx",rowNames=T)

