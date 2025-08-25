##==================================================================##
##--R : KF eQTL group 
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
##--loading the eQTL results of KF 
##==================================================================##
eQTLres <- read.xlsx("./zzRes711.KF_eGene_SigSNPs_eQTL_50kb_group.xlsx")
dim(eQTLres);head(eQTLres);tail(eQTLres)  ## 79810
unique(eQTLres$GeneChr)
unique(eQTLres$chr)
as.data.frame(table(eQTLres$GeneChr))
as.data.frame(table(eQTLres$chr))

as.data.frame(table(eQTLres$Group1))
as.data.frame(table(eQTLres$Group2))
as.data.frame(table(eQTLres$Group3))


## 1. boxplot for pval: Distant_InterChr vs Distant_IntraChr vs Local
pdf("zzzRes721.eQTL_Group2_Pval_Boxplot.pdf",width=6,height = 8)
  boxplot(pval~Group2,data = eQTLres) 
dev.off()
wilcox.test(x = eQTLres[eQTLres$Group2=="Distant_InterChr","pval"],
            y = eQTLres[eQTLres$Group2=="Distant_IntraChr","pval"],
            alternative = c("two.sided"))

wilcox.test(x = eQTLres[eQTLres$Group2=="Distant_IntraChr","pval"],
            y = eQTLres[eQTLres$Group2=="Local","pval"],
            alternative = c("two.sided"))



## 2. boxplot for PVE: Local vs Distant
pdf("zzzRes722.eQTL_Group1_PVE_Boxplot.pdf",width=3,height = 5)
  boxplot(PVE~Group1,data = eQTLres) 
dev.off()
wilcox.test(x = eQTLres[eQTLres$Group1=="Local","PVE"],
            y = eQTLres[eQTLres$Group1=="Distant","PVE"],
            alternative = c("two.sided"))




## 3. statistic for eGene
head(eQTLres);dim(eQTLres)
length(unique(eQTLres$GeneName))   ##  12293
length(unique(eQTLres[eQTLres$Group1=="Local","GeneName"]))    ## 3304
length(unique(eQTLres[eQTLres$Group1=="Distant","GeneName"]))  ## 11978

Number_eQTL_per_eGene <- as.data.frame(table(eQTLres$GeneName))
colnames(Number_eQTL_per_eGene) <- c("eGene","NumeQTL")
dim(Number_eQTL_per_eGene)  ## 12293     2
head(Number_eQTL_per_eGene)

NumeQTL_FreqTable <- as.data.frame(table(Number_eQTL_per_eGene$NumeQTL))
freqtable <- data.frame(NumeQTL_FreqTable,
                        frequency=NumeQTL_FreqTable$Freq/sum(NumeQTL_FreqTable$Freq),
                        freqCum=cumsum(NumeQTL_FreqTable$Freq),
                        frequencyCum=cumsum(NumeQTL_FreqTable$Freq)/sum(NumeQTL_FreqTable$Freq))
# Var1 Freq    frequency freqCum frequencyCum
# 1     1 1987 1.616367e-01    1987    0.1616367
# 2     2 1805 1.468315e-01    3792    0.3084682
# 3     3 1504 1.223461e-01    5296    0.4308143
# 4     4 1176 9.566420e-02    6472    0.5264785
# 5     5  893 7.264297e-02    7365    0.5991215
# 6     6  685 5.572277e-02    8050    0.6548442
# 7     7  583 4.742536e-02    8633    0.7022696
# 8     8  520 4.230050e-02    9153    0.7445701
# 9     9  416 3.384040e-02    9569    0.7784105
# 10   10  361 2.936631e-02    9930    0.8077768


pdf("zzzRes723.KF_Number_eQTL_vs_eGene_Barplot.pdf",width = 7,height = 5.5)
barplot(height = NumeQTL_FreqTable$Freq,
        names.arg = NumeQTL_FreqTable$Var1,
        ylim = c(0,2000),
        main = "KF",
        xlab = "Number of eQTL",ylab = "Number of eGene")
box()
dev.off()




















