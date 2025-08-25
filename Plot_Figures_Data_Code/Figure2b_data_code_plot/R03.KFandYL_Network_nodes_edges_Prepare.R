##==================================================================##
##--Description: prepare the results to plot the network
##--Highlight: for KaiFeng (KF) + YangLing (YL)
##--Maintainer: GuoLiang Li(lig@ipk-gatersleben.de)
##--Date: 2023.10.14
##==================================================================##
library(openxlsx)
library(dplyr)


##==================================================================##
##--network edges file
##==================================================================##
KF_edges <- read.xlsx("../01.KF/Res01.KF_Networks/KF_ASV_GWAS_TWAS_eGWAS.edges.xlsx")
head(KF_edges);tail(KF_edges);dim(KF_edges)
YL_edges <- read.xlsx("../02.YL/Res01.YL_Networks/ASV_GWAS_TWAS_eGWAS.edges.xlsx")
head(YL_edges);tail(YL_edges);dim(YL_edges)

## for edges
Edges <- rbind(KF_edges,YL_edges)
rownames(Edges) <- NULL;dim(Edges)   ## 1915
Edges <- unique(Edges[,c(1,2)])
dim(Edges)  ## 1888
write.xlsx(x = Edges,file = "./Res01.KFandYL_Networks/KFandYL_ASV_GWAS_TWAS_eGWAS.edges.xlsx")


## for nodes
Nodes <- data.frame(ID=unique(c(Edges$source,Edges$target)),module="ASV")
Nodes[grep(pattern = "QTL",x = Nodes$ID),"module"] <- "QTL"
Nodes[grep(pattern = "Bna",x = Nodes$ID),"module"] <- "Gene"
Nodes[grep(pattern = "Hot",x = Nodes$ID),"module"] <- "Hotspot"

as.data.frame(table(Nodes$module))
write.xlsx(x = Nodes,file = "./Res01.KFandYL_Networks/KFandYL_ASV_GWAS_TWAS_eGWAS.nodes.xlsx")



