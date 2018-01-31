########### Prereqs ###########
library(tidyverse)
library(ReactomePA)
library(curl)

options(stringsAsFactors=FALSE)
setwd("~/gdrive/BearOmics2/EnrichmentAnalysis")
load("../data/combData_20180127.RData")

########### Get Metabolite ChEBIs ###########
idk<-combData$metabs$key
write.table(idk$pubchem[complete.cases(idk$pubchem)],file="pubChems.txt",
            row.names=FALSE,col.names=FALSE)

curl("http://cts.fiehnlab.ucdavis.edu/service/convet/PubChem%20CID/ChEBI/10917802")


########### Example GSEA using ReactomePA ###########
S_Up<-read.table(file="../data/CONTROL_UP_VS_CONTROL_S_ALL.txt",header=TRUE,sep="\t")
names(S_Up)[names(S_Up)=="log2FC.CONTROL_UP.CONTROL_S."]<-"logFC"
S_Up<-S_Up %>% filter(p_value<.2)
S_Up<-S_Up %>% arrange(p_value,desc(abs(logFC)))
S_Up$order<-nrow(S_Up):1
S_Up2<-S_Up
S_Up2b<-S_Up2$order
names(S_Up2b)<-S_Up2$ENTREZ.ID

# Over representation (hypergeometric):
S_Up3<-S_Up %>% filter(p_value<.05)
over_S_Up<-as.data.frame(enrichPathway(S_Up3$ENTREZ.ID,pvalueCutoff=.2,readable=TRUE))

# Actual GSEA
gse_S_Up<-gsePathway(geneList=S_Up2b,nPerm=1000,minGSSize=10,pvalueCutoff=0.25,
           pAdjustMethod="BH",verbose=TRUE)
gseSum_S_Up<-as.data.frame(gse_S_Up)
gseaplot(gse_S_Up,geneSetID="R-HSA-174824")

# Need to see fgsea::fgsea