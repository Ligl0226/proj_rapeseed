##==================================================================##
##--R : YL eQTL group 
##--Highlight: summary the eQTL results
##--Maintainer: GuoLiang Li(lig@ipk-gatersleben.de)
##--Date: 2023.10.16/17
##==================================================================##
library(openxlsx)
library(ggplot2)
library(dplyr)
library(gplots)
library(ggpubr)


##==================================================================##
##--loading the eQTL results of YL 
##==================================================================##
eQTLres <- read.xlsx("./zzRes722.YL_eGene_SigSNPs_eQTL_50kb_group.xlsx")
dim(eQTLres);head(eQTLres);tail(eQTLres)  ## 121,955
unique(eQTLres$GeneChr)
unique(eQTLres$chr)
as.data.frame(table(eQTLres$GeneChr))
as.data.frame(table(eQTLres$chr))

as.data.frame(table(eQTLres$Group1))
as.data.frame(table(eQTLres$Group2))
as.data.frame(table(eQTLres$Group3))


## 1. boxplot for pval: Distant_InterChr vs Distant_IntraChr vs Local
pdf("zzzRes731.eQTL_Group2_Pval_Boxplot.pdf",width=6,height = 8)
  boxplot(pval~Group2,data = eQTLres) 
dev.off()

wilcox.test(x = eQTLres[eQTLres$Group2=="Distant_InterChr","pval"],
            y = eQTLres[eQTLres$Group2=="Distant_IntraChr","pval"],
            alternative = c("two.sided"))

wilcox.test(x = eQTLres[eQTLres$Group2=="Distant_IntraChr","pval"],
            y = eQTLres[eQTLres$Group2=="Local","pval"],
            alternative = c("two.sided"))


## 2. boxplot for PVE: Local vs Distant
pdf("zzzRes732.YL_eQTL_Group1_PVE_Boxplot.pdf",width=3,height = 5)
  boxplot(PVE~Group1,data = eQTLres) 
dev.off()
wilcox.test(x = eQTLres[eQTLres$Group1=="Local","PVE"],
            y = eQTLres[eQTLres$Group1=="Distant","PVE"],
            alternative = c("two.sided"))


## 3. statistic for eGene
head(eQTLres);dim(eQTLres)
length(unique(eQTLres$GeneName))   ##  12322
length(unique(eQTLres[eQTLres$Group1=="Local","GeneName"]))    ## 4566
length(unique(eQTLres[eQTLres$Group1=="Distant","GeneName"]))  ## 11976

Number_eQTL_per_eGene <- as.data.frame(table(eQTLres$GeneName))
colnames(Number_eQTL_per_eGene) <- c("eGene","NumeQTL")
dim(Number_eQTL_per_eGene)  ## 12322     2
head(Number_eQTL_per_eGene)

NumeQTL_FreqTable <- as.data.frame(table(Number_eQTL_per_eGene$NumeQTL))
freqtable <- data.frame(NumeQTL_FreqTable,
                        frequency=NumeQTL_FreqTable$Freq/sum(NumeQTL_FreqTable$Freq),
                        freqCum=cumsum(NumeQTL_FreqTable$Freq),
                        frequencyCum=cumsum(NumeQTL_FreqTable$Freq)/sum(NumeQTL_FreqTable$Freq))
# Var1 Freq    frequency freqCum frequencyCum
# 1     1 1432 1.162149e-01    1432    0.1162149
# 2     2 1382 1.121571e-01    2814    0.2283720
# 3     3 1142 9.267976e-02    3956    0.3210518
# 4     4  907 7.360818e-02    4863    0.3946600
# 5     5  798 6.476221e-02    5661    0.4594222
# 6     6  688 5.583509e-02    6349    0.5152573
# 7     7  569 4.617757e-02    6918    0.5614348
# 8     8  489 3.968512e-02    7407    0.6011199
# 9     9  460 3.733160e-02    7867    0.6384516
# 10   10  414 3.359844e-02    8281    0.6720500


pdf("zzzRes733.YL_Number_eQTL_vs_eGene_Barplot.pdf",width = 7,height = 5.5)
barplot(height = NumeQTL_FreqTable$Freq,
        names.arg = NumeQTL_FreqTable$Var1,
        ylim = c(0,2000),
        main = "YL",
        xlab = "Number of eQTL",ylab = "Number of eGene")
box()
dev.off()
















