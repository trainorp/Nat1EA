library(tidyverse)
library(ReactomePA)

options(stringsAsFactors=FALSE)
setwd("~/gdrive/BearOmics2/EnrichmentAnalysis")
load("../data/combData_20180127.RData")

########### Univariate tests ###########
S_Up<-read.table(file="../data/CONTROL_UP_VS_CONTROL_S_ALL.txt",header=TRUE,sep="\t")
names(S_Up)[names(S_Up)=="log2FC.CONTROL_UP.CONTROL_S."]<-"logFC"
S_Up2<-S_Up %>% filter(p_value<.05)
over_S_Up<-as.data.frame(enrichPathway(S_Up2$ENTREZ.ID,pvalueCutoff=.2,readable=TRUE))

S_Up<-S_Up %>% arrange(p_value)
gse_S_Up<-gsePathway(geneList=S_Up$ENSEMBL.GENE[1:100],nPerm=1000,minGSSize=120, pvalueCutoff=0.2,
           pAdjustMethod="BH", verbose=FALSE)

########### Their example ###########
# https://bioconductor.org/packages/release/bioc/vignettes/ReactomePA/inst/doc/ReactomePA.html
data(geneList)
de <- names(geneList)[abs(geneList) > 1.5]
x <- enrichPathway(gene=de,pvalueCutoff=0.05, readable=T)
head(as.data.frame(x))

y <- gsePathway(geneList, nPerm=1000,
                minGSSize=120, pvalueCutoff=0.2,
                pAdjustMethod="BH", verbose=FALSE)
res <- as.data.frame(y)
head(res)
gseaplot(y, geneSetID = "R-HSA-69242")
