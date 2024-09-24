##==================================================================##
##--Description: R functions to perform TWAS (version 1.0)
##--TWAS (assication between phenotype and Gene Expression)
##--Maintainer: Dr. GuoLiang Li(lig@ipk-gatersleben.de or guoliangli0226@gmail.com)
##--Date: 2023.07.31
##==================================================================##
library(BGLR)


##==================================================================##
##--function to calculate the Relationship matrix based on Gene Expression
##==================================================================##
# loading Gene Expression data
# rows: genotypes
# colums: genes
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


Fit_nullModel <- function(y,M,G=NULL,X=NULL,rec_G=FALSE,tempBGLRdir="./"){
  ## y is the one-column matrix of phenoytpe
  ## M is the matrix of Gene Expression
  ## G is the Relationship matrix based on Gene Expression
  ## X is the Fixed effect which need to fit in the model
  ## rec_G is indecate if save the G matrix
  ## tempBGLRdir is the temp dir of BGLR
  
  res <- list()
  
  if (is.null(G)){
    G <- KinMatGen(M,is.SNPs = FALSE)
  }
  
  if (is.null(X))  {
    ETA <- list(list(K=G,model="RKHS"))
  }else{
    ETA <- list(list(X=X,model="FIXED"),
                list(K=G,model="RKHS"))
  }
  
  fit_BGLR_null <- BGLR(y=y,
                        ETA=ETA,
                        nIter=10000,
                        burnIn=1000,
                        saveAt=paste(tempBGLRdir,"fit_nullmodel_",sep=""),
                        verbose=FALSE)
  
  res$varcomp <- array(0,2)				
  if (is.null(X)){
    res$varcomp[1] <- fit_BGLR_null$ETA[[1]]$varU
  }else{
    res$varcomp[1] <- fit_BGLR_null$ETA[[2]]$varU
  }
  names(res$varcomp) <- c("Gexp","R")
  res$varcomp["R"] <- fit_BGLR_null$varE
  
  if (rec_G==TRUE){
    res$Gmat <- G
  }
  return(res)         
}


TWAS_P3D <- function(y,M,G=NULL,X=NULL,tempBGLRdir="./"){
  ## y is the one-column matrix of phenoytpe
  ## M is the matrix of Gene Expression
  ## G is the Relationship matrix based on Gene Expression
  ## X is the Fixed effect which need to fit in the model
  ## tempBGLRdir is the temp dir of BGLR
  
  n <- nrow(M)
  p <- ncol(M)
  marname <- colnames(M)
  
  res_nullmod <- Fit_nullModel(y,M,G,X,rec_G=TRUE,tempBGLRdir)
  
  varparam <- res_nullmod$varcomp
  G <- res_nullmod$Gmat
  
  lambda <- varparam[1]/varparam[2]
  Kfinal <- matrix(0,n,n)  
  Kfinal <- Kfinal + lambda*G
  
  eigK <- eigen(Kfinal)
  d <- eigK$values
  W <- eigK$vectors
  wt <- 1/(d+1)
  TRAN <- diag(sqrt(wt))%*%t(W)
  
  Y <- TRAN%*%y
  mu <- rep(1,n)
  tmu <- TRAN%*%mu
  
  res <- matrix(0,p,4)
  colnames(res) <- c("Estimated_eff","t_statistic","P_value","PVE")
  rownames(res) <- marname
  
  for (i in 1:p){   
    codlist <- M[,i]
    subfinalDATA <- data.frame(y=Y,x=tmu,A=TRAN%*%codlist)
    fit.lm <- lm(as.formula(paste("y ~ -1 + x + A")),data=subfinalDATA)
    res_table <- summary(fit.lm)$coefficient
    
    if ("A" %in% rownames(res_table)){
      res[i,c(1:3)] <- res_table["A",c(1,3,4)]
      ## for PVE
      anovatable <- anova(fit.lm)
      PVEi <- anovatable["A","Sum Sq"]/sum(anovatable[,"Sum Sq"])
      res[i,4] <- PVEi
    }else{
      res[i,] <- c(NA,NA,NA)
    }
    ## if (i %% 1000 == 0){
    ##    cat("Marker",i,"completed\n")
    ## } 
  }
  return(res)
}
