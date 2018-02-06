########### Prereqs ###########
library(tidyverse)
library(ReactomePA)
library(curl)
library(reactome.db)
library(fgsea)
library(emmeans)

options(stringsAsFactors=FALSE)
setwd("~/gdrive/BearOmics2/EnrichmentAnalysis")

########### Get Metabolite ChEBIs ###########
load("../data/combData_20180127.RData")
metabKey<-combData$metabs$key

# PubChem to ChEBI
pcLook<-data.frame(pubchem=na.omit(unique(metabKey$pubchem)),ChEBI=NA)
for(i in 1:nrow(pcLook)){
  url1<-paste0("http://cts.fiehnlab.ucdavis.edu/service/convert/PubChem%20CID/chebi/",
               pcLook$pubchem[i])
  curl1<-curl_fetch_memory(url1)
  curl2<-jsonlite::fromJSON(paste0(rawToChar(curl1$content)))
  pcLook$ChEBI[i]<-gsub("CHEBI:","",paste0(unlist(curl2$result),collapse=";"))
  print(i)
}
# KEGG to ChEBI
keggLook<-data.frame(kegg=na.omit(unique(metabKey$kegg)),ChEBI=NA)
keggLook<-keggLook %>% filter(keggLook$kegg!="")
for(i in 1:nrow(keggLook)){
  url1<-paste0("http://cts.fiehnlab.ucdavis.edu/service/convert/kegg/chebi/",
               keggLook$kegg[i])
  curl1<-curl_fetch_memory(url1)
  curl2<-jsonlite::fromJSON(paste0(rawToChar(curl1$content)))
  keggLook$ChEBI[i]<-gsub("CHEBI:","",paste0(unlist(curl2$result),collapse=";"))
  print(i)
}
# HMDB to ChEBI
hmdbLook<-data.frame(hmdb=na.omit(unique(metabKey$hmdb)),ChEBI=NA)
hmdbLook<-hmdbLook %>% filter(hmdbLook$hmdb!="")
for(i in 1:nrow(hmdbLook)){
  url1<-paste0("http://cts.fiehnlab.ucdavis.edu/service/convert/Human%20Metabolome%20Database/chebi/",
               hmdbLook$hmdb[i])
  curl1<-curl_fetch_memory(url1)
  curl2<-jsonlite::fromJSON(paste0(rawToChar(curl1$content)))
  hmdbLook$ChEBI[i]<-gsub("CHEBI:","",paste0(unlist(curl2$result),collapse=";"))
  print(i)
}
# CAS to ChEBI:
metabKey$cas<-gsub("/","-",metabKey$cas)
casLook<-data.frame(cas=na.omit(unique(metabKey$cas)),ChEBI=NA)
casLook<-casLook %>% filter(casLook$cas!="")
for(i in 1:nrow(casLook)){
  multCas<-unlist(strsplit(casLook$cas[i],";"))
  for(casi in multCas){
    url1<-paste0("http://cts.fiehnlab.ucdavis.edu/service/convert/CAS/chebi/",casi)
    curl1<-curl_fetch_memory(url1)
    curl2<-jsonlite::fromJSON(paste0(rawToChar(curl1$content)))
    casLook$ChEBI[i]<-gsub("CHEBI:","",paste0(unlist(curl2$result),collapse=";"))
  }
  print(i)
}

# Add to metabolite key:
metabKey$pubchem<-as.character(metabKey$pubchem)
metabKey$pubchem[is.na(metabKey$pubchem)]<-""
metabKey$cas[is.na(metabKey$cas)]<-""
metabKey$ChEBI<-""
for(i in 1:nrow(metabKey)){
  tempChEBI<-""
  if(metabKey$pubchem[i]!=""){
    if(pcLook$ChEBI[pcLook$pubchem==metabKey$pubchem[i]]!=""){
      tempChEBI<-paste(tempChEBI,pcLook$ChEBI[pcLook$pubchem==metabKey$pubchem[i]],sep="!")
    }
  }
  if(metabKey$kegg[i]!=""){
    if(keggLook$ChEBI[keggLook$kegg==metabKey$kegg[i]]!=""){
      tempChEBI<-paste(tempChEBI,keggLook$ChEBI[keggLook$kegg==metabKey$kegg[i]],sep="*")
    }
  }
  if(metabKey$hmdb[i]!=""){
    if(hmdbLook$ChEBI[hmdbLook$hmdb==metabKey$hmdb[i]]!=""){
      tempChEBI<-paste(tempChEBI,hmdbLook$ChEBI[hmdbLook$hmdb==metabKey$hmdb[i]],sep="+")
    }
  }
  if(metabKey$cas[i]!=""){
    if(casLook$ChEBI[casLook$cas==metabKey$cas[i]]!=""){
      tempChEBI<-paste(tempChEBI,casLook$ChEBI[casLook$hmdb==metabKey$hmdb[i]],sep="^")
    }
  }
  metabKey$ChEBI[i]<-tempChEBI
}

# Find unique ChEBIs only:
for(i in 1:nrow(metabKey)){
  ChEBIs<-unlist(strsplit(metabKey$ChEBI[i],split="[[:punct:]]"))
  ChEBIs<-ChEBIs[ChEBIs!=""]
  if(length(ChEBIs)>0){
    ChEBIs<-unique(ChEBIs)
    ChEBIs<-paste(ChEBIs,collapse=";")
    metabKey$ChEBI[i]<-ChEBIs
  }
  else{
    metabKey$ChEBI[i]<-""
  }
}
rm(casLook,hmdbLook,keggLook,pcLook,casi,curl1,curl2,ChEBIs,i,url1,
   multCas,tempChEBI)
combData$metabs$key<-metabKey
save(combData,file="../data/combData_20180202.RData")

########### Reload data ###########
load("../data/combData_20180202.RData")
metabKey<-combData$metabs$key

# Reactome gene data:
# ensgReact<-read.table("Ensembl2Reactome_20180202.txt",comment.char="",sep="\t")
# names(ensgReact)<-c("ensg","pathID","url","pathName","evidCode","species")
# ensgReact<-ensgReact %>% filter(species=="Homo sapiens") %>% dplyr::select(-species)

# Reactome compound data:
chebiReact<-read.csv("ChEBI2Reactome_AllPathways_20180202.csv",header=FALSE)
names(chebiReact)<-c("ChEBI","peID","peName","pathID","url","pathName","evidCode","species")
chebiReact<-chebiReact %>% filter(species=="Homo sapiens") %>% dplyr::select(-species)

# How many are in a Reactome pathway:
metabKey$inReactome<-NA
for(i in 1:nrow(metabKey)){
  ChEBIs<-unlist(strsplit(metabKey$ChEBI[i],";"))
  if(length(ChEBIs)>0){
    temp2<-0
    for(ChEBI in ChEBIs){
      temp<-match(ChEBI,chebiReact$ChEBI,nomatch=0)
      temp2<-temp2+ifelse(length(temp)==0 | temp==0,0,1)
    }
  }else{
    temp2<-0
  }
  metabKey$inReactome[i]<-temp2
}

########### fgsea with custom set from Reactome ###########
# Make list of sets
chebiReact$ChEBI<-as.character(chebiReact$ChEBI)
metabSet<-list()
uniquePaths<-unique(chebiReact$pathID)
for(path in uniquePaths){
  metabSet[[path]]<-unique(chebiReact$ChEBI[chebiReact$pathID==path])
}

metabExpr<-combData$metabs$expr
metabExpr<-combData$metabs$pheno %>% dplyr::select(sid,pheno) %>% left_join(metabExpr)
metabExpr$pheno<-factor(metabExpr$pheno,levels=levels(metabExpr$pheno)[c(2,1,3:6)])
diffs<-data.frame(metab=names(metabExpr)[3:ncol(metabExpr)],SvsU=NA)
for(i in 1:nrow(diffs)){
  df1<-data.frame(pheno=metabExpr$pheno,metab=metabExpr[,diffs$metab[i]])
  lm1<-lm(metab~pheno,data=df1)
  diff1<-as.data.frame(pairs(emmeans(lm1,"pheno")))
  diffs$SvsU[i]<-diff1$p.value[diff1$contrast=="Scrambled - Up"]
}
diffs<-metabKey %>% left_join(diffs,by=c("id"="metab"))
diffs$ChEBI2<-sapply(strsplit(diffs$ChEBI,";"),function(x) x[1])
metabStats<-diffs$SvsU
names(metabStats)<-diffs$ChEBI2
metabStats<-metabStats[!is.na(names(metabStats))]

metabGSEA<-fgsea(metabSet,metabStats,nperm = 1000)

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