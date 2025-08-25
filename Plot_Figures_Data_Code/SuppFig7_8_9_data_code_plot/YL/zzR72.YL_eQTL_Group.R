##==================================================================##
##--R : YL eQTL group 
##--Highlight: based on the results using RNAseq SNPs after maf5het95 filtering
##--Maintainer: GuoLiang Li(lig@ipk-gatersleben.de)
##--Date: 2023.10.10
##==================================================================##
library(openxlsx)
library(ggplot2)
library(dplyr)
library(gplots)


##==================================================================##
##--loading the eQTL results of YL 
##==================================================================##
eQTLres <- read.xlsx("./zRes722.YL_eGene_SigSNPs_eQTL_50kb.xlsx")
dim(eQTLres);head(eQTLres);tail(eQTLres)  ## 144591
unique(eQTLres$GeneChr)
unique(eQTLres$chr)
as.data.frame(table(eQTLres$GeneChr))
as.data.frame(table(eQTLres$chr))

## only focus on the gene located the clear chrom
eQTLres <- eQTLres[grep(pattern = "random",x = eQTLres$GeneChr,invert = T),]
rownames(eQTLres) <- NULL
dim(eQTLres);head(eQTLres);tail(eQTLres)  ## 121955
unique(eQTLres$GeneChr)
unique(eQTLres$chr)
as.data.frame(table(eQTLres$GeneChr))
as.data.frame(table(eQTLres$chr))


## change the name of chromosome
chrnameold1 <- c("chrA01","chrA02","chrA03","chrA04","chrA05","chrA06","chrA07","chrA08","chrA09","chrA10",
                 "chrC01","chrC02","chrC03","chrC04","chrC05","chrC06","chrC07","chrC08","chrC09")
chrnameold2 <- paste("chr",c(1:19),sep = "")
chrnamenew <- paste("Chr",sprintf("%02d",1:19),sep = "")

for(i in 1:length(chrnamenew)){
  chrnameold1i <- chrnameold1[i]
  chrnameold2i <- chrnameold2[i]
  chrnamenewi <- chrnamenew[i]
  cat(chrnameold1i,chrnameold2i,chrnamenewi,"\n")
  eQTLres[eQTLres$GeneChr==chrnameold1i,"GeneChr"] <- chrnamenewi
  eQTLres[eQTLres$chr==chrnameold2i,"chr"] <- chrnamenewi
}
dim(eQTLres);head(eQTLres);tail(eQTLres)
as.data.frame(table(eQTLres$GeneChr))
as.data.frame(table(eQTLres$chr))



##==================================================================##
##--the distance between eQTL and eGene (kb)
##==================================================================##
RegionDist <- function(x){
  # x <- eQTLres[3,c("GeneChr","GeneStart","GeneEnd","chr","QTLstart","QTLend")]
  if(x[1]==x[4]){   ## some chromosome
    a <- as.numeric(x[2])
    b <- as.numeric(x[3])
    c <- as.numeric(x[5])
    d <- as.numeric(x[6])
    if(a<c & b<c){
      RegionDistValue <- c-b
    }else if(a>d & b>d){
      RegionDistValue <- a-d
    }else{
      RegionDistValue <- 0
    }
  }else{            ## diff chromosome
    RegionDistValue <- NA
  }
  return(RegionDistValue)
}

eQTLres$distance <- apply(eQTLres[,c("GeneChr","GeneStart","GeneEnd","chr","QTLstart","QTLend")],1,RegionDist)
eQTLres$distance <- eQTLres$distance/1000   ## to kb
head(eQTLres)
sum(is.na(eQTLres$distance))        ## 95150
sum(eQTLres$distance==0,na.rm = T)  ## 3915
sum(eQTLres$distance>0,na.rm = T)   ## 22890

sum(eQTLres$distance<50,na.rm = T)    ## 4720
sum(eQTLres$distance<500,na.rm = T)   ## 11393
sum(eQTLres$distance<1000,na.rm = T)   ## 14607
sum(eQTLres$distance<2000,na.rm = T)   ## 18095

sum(eQTLres$distance>500,na.rm = T)   ## 15412

summary(eQTLres[!is.na(eQTLres$distance),"distance"])
boxplot(eQTLres$distance)

pdf(file = "zzRes721.Density_distances_eQTL_eGene.pdf",width = 12,height = 5)
ggplot(eQTLres,aes(x=distance))+
  geom_density(color="black")+
  scale_x_continuous(breaks = seq(0,60000,2000))+
  xlab("The distance between eQTL and eGene (kb)")

ggplot(eQTLres %>% filter(distance>0 & distance<10000),aes(x=distance))+
  geom_density(color="black")+
  scale_x_continuous(breaks = seq(0,10000,500))+
  xlab("The distance between eQTL and eGene (kb)")

ggplot(eQTLres %>% filter(distance>0 & distance<5000),aes(x=distance))+
  geom_density(color="black")+
  scale_x_continuous(breaks = seq(0,5000,200))+
  xlab("The distance between eQTL and eGene (kb)")

ggplot(eQTLres %>% filter(distance>0 & distance<2500),aes(x=distance))+
  geom_density(color="black")+
  scale_x_continuous(breaks = seq(0,2500,100))+
  xlab("The distance between eQTL and eGene (kb)")
dev.off()



##==================================================================##
##--eQTL grouping (local, distant)
##==================================================================##
dim(eQTLres);head(eQTLres);tail(eQTLres)
DistThred <- 50   ## kb


## Group1: Local vs Distant
eQTLres$Group1 <- "NA"
eQTLres[is.na(eQTLres$distance) | eQTLres$distance>DistThred,"Group1"] <- "Distant"
eQTLres[!is.na(eQTLres$distance) & eQTLres$distance<=DistThred,"Group1"] <- "Local"

as.data.frame(table(eQTLres$Group1));
#      Var1  Freq
# 1 Distant 117235
# 2   Local  4720
as.data.frame(prop.table(table(eQTLres$Group1)))
#      Var1       Freq
# 1 Distant 0.9612972
# 2   Local 0.0387028

aa <- unique(eQTLres %>% filter(Group1=="Distant") %>% pull(GeneName))
bb <- unique(eQTLres %>% filter(Group1=="Local") %>% pull(GeneName))
length(aa)   ## 11976
length(bb)   ## 4566
length(intersect(aa,bb))  ## 4220
venn(data = list(aa,bb))


## Group2: Local vs Distant_IntraChr vs Distant_InterChr
eQTLres$Group2 <- "NA"
eQTLres[is.na(eQTLres$distance),"Group2"] <- "Distant_InterChr"
eQTLres[!is.na(eQTLres$distance) & eQTLres$distance>DistThred,"Group2"] <- "Distant_IntraChr"
eQTLres[!is.na(eQTLres$distance) & eQTLres$distance<=DistThred,"Group2"] <- "Local"

as.data.frame(table(eQTLres$Group2));
#               Var1  Freq
# 1 Distant_InterChr 95150
# 2 Distant_IntraChr 22085
# 3            Local  4720
as.data.frame(prop.table(table(eQTLres$Group2)))
#               Var1       Freq
# 1 Distant_InterChr 0.7802058
# 2 Distant_IntraChr 0.1810914
# 3            Local 0.0387028


## Group3: Local vs Distant_IntraSubgenome vs Distant_InterSubgenome
eQTLres$GeneSubGenome <- NA
eQTLres$eQTLsubGenome <- NA

for(i in 1:10){
  chrnamenewi <- chrnamenew[i]
  eQTLres[eQTLres$GeneChr==chrnamenewi,"GeneSubGenome"] <- "A"
  eQTLres[eQTLres$chr==chrnamenewi,"eQTLsubGenome"] <- "A"
}

for(i in 11:19){
  chrnamenewi <- chrnamenew[i]
  eQTLres[eQTLres$GeneChr==chrnamenewi,"GeneSubGenome"] <- "C"
  eQTLres[eQTLres$chr==chrnamenewi,"eQTLsubGenome"] <- "C"
}

as.data.frame(table(eQTLres$GeneSubGenome))
as.data.frame(table(eQTLres$eQTLsubGenome))

eQTLres$Group3 <- NA
eQTLres[!is.na(eQTLres$distance) & eQTLres$distance<=DistThred,"Group3"] <- "Local"
eQTLres[eQTLres$GeneSubGenome==eQTLres$eQTLsubGenome & eQTLres$GeneChr!=eQTLres$chr,"Group3"] <- "Distant_IntraSubgenome"
eQTLres[eQTLres$GeneSubGenome==eQTLres$eQTLsubGenome & (eQTLres$GeneChr==eQTLres$chr & eQTLres$distance>DistThred),"Group3"] <- "Distant_IntraSubgenome"
eQTLres[eQTLres$GeneSubGenome!=eQTLres$eQTLsubGenome,"Group3"] <- "Distant_InterSubgenome"

as.data.frame(table(eQTLres$Group3));
#                     Var1  Freq
# 1 Distant_InterSubgenome 58624
# 2 Distant_IntraSubgenome 58611
# 3                  Local  4720
as.data.frame(prop.table(table(eQTLres$Group3)))
#                     Var1       Freq
# 1 Distant_InterSubgenome 0.4807019
# 2 Distant_IntraSubgenome 0.4805953
# 3                  Local 0.0387028



##==================================================================##
##--save the results
##==================================================================##
dim(eQTLres);head(eQTLres);tail(eQTLres)  ## 121955
length(unique(eQTLres$GeneName))          ## 12322
write.xlsx(x = eQTLres,file = "zzRes722.YL_eGene_SigSNPs_eQTL_50kb_group.xlsx")



