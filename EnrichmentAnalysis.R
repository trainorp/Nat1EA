########### Prereqs ###########
library(tidyverse)
library(ReactomePA)
library(curl)
library(reactome.db)
library(fgsea)
library(emmeans)
library(KEGGREST)
library(pathview)

options(stringsAsFactors=FALSE)
setwd("~/gdrive/BearOmics2/EnrichmentAnalysis")

########### Get Metabolite ChEBIs ###########
## Start Not run
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
## End Not run

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

########### Metabolite Differential abundance ###########
metabExpr<-combData$metabs$expr
metabExpr<-combData$metabs$pheno %>% dplyr::select(sid,pheno) %>% left_join(metabExpr)
metabExpr$pheno<-factor(metabExpr$pheno,levels=levels(metabExpr$pheno)[c(2,1,3:6)])
comparisons<-c("Up","Down","CRISPR-5-50","CRISPR-2-12","CRISPR-2-19")
diffs<-expand.grid(metab=names(metabExpr)[3:ncol(metabExpr)],ref="Scrambled",
                   comparison=comparisons,FC=NA,pVal=NA,stringsAsFactors=FALSE)
metabs<-names(metabExpr)[3:ncol(metabExpr)]
for(i in 1:length(metabs)){
  df1<-data.frame(pheno=metabExpr$pheno,metab=metabExpr[,metabs[i]])
  lm1<-lm(metab~pheno,data=df1)
  diff1<-as.data.frame(pairs(emmeans(lm1,"pheno")))
  for(comparison in comparisons){
    diffs$FC[diffs$metab==metabs[i] & diffs$comparison==comparison]<-
      diff1$estimate[diff1$contrast==paste0("Scrambled - ",comparison)]
    diffs$pVal[diffs$metab==metabs[i] & diffs$comparison==comparison]<-
      diff1$p.value[diff1$contrast==paste0("Scrambled - ",comparison)]
  }
  print(i)
}
FCs<-diffs %>% dplyr::select(-pVal,-ref) %>% spread(comparison,FC)
names(FCs)[names(FCs)!="metab"]<-paste("FC",names(FCs)[names(FCs)!="metab"],sep="_")
pVals<-diffs %>% dplyr::select(-FC,-ref) %>% spread(comparison,pVal)
names(pVals)[names(pVals)!="metab"]<-paste("pVal",names(pVals)[names(pVals)!="metab"],sep="_")
diffs<-FCs %>% left_join(pVals)

diffs<-metabKey %>% left_join(diffs,by=c("id"="metab"))
temp1<-diffs[diffs$biochemical=="carnitine",]
temp1$kegg<-"C00487"
diffs<-rbind(diffs,temp1)

# Rank orders
for(comparison in comparisons){
  ord<-order(diffs[,paste0("pVal_",comparison)],diffs[,paste0("FC_",comparison)],
             decreasing=c(TRUE,FALSE))
  diffs<-diffs[ord,]
  diffs<-cbind(diffs,ord=1:nrow(diffs))
  names(diffs)[names(diffs)=="ord"]<-paste0("ord_",comparison)
}

########### Transcript Differential Expression ############
S_Up<-read.table(file="../data/CONTROL_UP_VS_CONTROL_S_ALL.txt",header=TRUE,sep="\t")
S_Down<-read.table(file="../data/T_DOWN_VS_CONTROL_S_ALL.txt",header=TRUE,sep="\t")
S_CRISPR219<-read.table(file="../data/T_2_19_VS_CONTROL_S_ALL.txt",header=TRUE,sep="\t")
S_CRISPR212<-read.table(file="../data/T_2_12_VS_CONTROL_S_ALL.txt",header=TRUE,sep="\t")
S_CRISPR550<-read.table(file="../data/T_5_50_VS_CONTROL_S_ALL.txt",header=TRUE,sep="\t")

names(S_Up)[names(S_Up)=="log2FC.CONTROL_UP.CONTROL_S."]<-"logFC"
names(S_Down)[names(S_Down)=="log2FC.T_DOWN.CONTROL_S."]<-"logFC"
names(S_CRISPR219)[names(S_CRISPR219)=="log2FC.T_2_19.CONTROL_S."]<-"logFC"
names(S_CRISPR212)[names(S_CRISPR212)=="log2FC.T_2_12.CONTROL_S."]<-"logFC"
names(S_CRISPR550)[names(S_CRISPR550)=="log2FC.T_5_50.CONTROL_S."]<-"logFC"

comparisons2<-c("S_Up","S_Down","S_CRISPR550","S_CRISPR212","S_CRISPR219")
for(comparison in comparisons2){
  temp1<-get(comparison)
  temp1<-temp1 %>% arrange(p_value,desc(abs(logFC)))
  temp1$order<-nrow(temp1):1
  assign(comparison,temp1)
}

########### Get KEGG data ###########
# List of Pathways:
keggPathways<-keggList("pathway","hsa")
keggPathways<-data.frame(pathId=names(keggPathways),pathName=keggPathways)
keggPathways$pathId<-gsub("path:hsa","",keggPathways$pathId)

# Pathway Genes:
keggPathwayGenes<-keggLink("pathway","hsa")
keggPathwayGenes<-data.frame(path=keggPathwayGenes,gene=names(keggPathwayGenes))
keggPathwayGenes$gene<-gsub("hsa:","",keggPathwayGenes$gene)
keggPathwayGenes$path<-gsub("pathway:","",keggPathwayGenes$path)
keggGeneSet<-list()
keggUniquePaths<-unique(keggPathwayGenes$path)
for(path in keggUniquePaths){
  keggGeneSet[[path]]<-unique(keggPathwayGenes$gene[keggPathwayGenes$path==path])
}

# Pathway Compounds:
keggPathwayCompounds<-keggLink("pathway","cpd")
keggPathwayCompounds<-data.frame(cpd=names(keggPathwayCompounds),path=keggPathwayCompounds)
keggPathwayCompounds$cpd<-gsub("cpd:","",keggPathwayCompounds$cpd)
keggPathwayCompounds$path<-gsub("path:","",keggPathwayCompounds$path)
keggMetabSet<-list()
keggUniquePaths<-unique(keggPathwayCompounds$path)
for(path in keggUniquePaths){
  keggMetabSet[[path]]<-unique(keggPathwayCompounds$cpd[keggPathwayCompounds$path==path])
}

########### KEGG Set Enrichment Analysis ###########
comparisonDf<-as.data.frame(cbind(comparisons,comparisons2))
for(i in 2:length(comparisonDf)){
  # Which comparison (separate for metabolites vs. transcripts)
  comparison<-comparisonDf$comparisons[i]
  comparison2<-comparisonDf$comparisons2[i]
  
  ########### KEGG Metabolite Set Enrichment Analysis ###########
  metabStats<-diffs[,paste0("ord_",comparison)]
  names(metabStats)<-diffs$kegg
  metabStats<-metabStats[!is.na(names(metabStats)) & names(metabStats)!=""]
  
  # Set enrichment analysis:
  KeggMetabGSEA<-fgsea(keggMetabSet,metabStats,nperm=10000,minSize=2,maxSize=Inf)
  
  ########### KEGG Gene Set Enrichment Analysis ###########
  geneDf<-get(comparison2)
  geneStats<-geneDf$order
  names(geneStats)<-geneDf$ENTREZ.ID
  
  # Set enrichment analysis:
  KeggGeneGSEA<-fgsea(keggGeneSet,geneStats,nperm=10000,minSize=2,maxSize=Inf)
  
  ########### Joined Metabolite & Gene Set Analysis ###########
  names(KeggMetabGSEA)<-paste("metab",names(KeggMetabGSEA),sep="_")
  KeggMetabGSEA$metab_pathway<-gsub("map","",KeggMetabGSEA$metab_pathway)
  names(KeggGeneGSEA)<-paste("gene",names(KeggGeneGSEA),sep="_")
  KeggGeneGSEA$gene_pathway<-gsub("path:hsa","",KeggGeneGSEA$gene_pathway)
  KeggGSEA<-KeggMetabGSEA %>% left_join(KeggGeneGSEA,by=c("metab_pathway"="gene_pathway"))
  KeggGSEA<-keggPathways %>% left_join(KeggGSEA,by=c("pathId"="metab_pathway"))
  KeggGSEA$minNES<-apply(KeggGSEA[,c("metab_NES","gene_NES")],1,function(x) min(x,na.rm=TRUE))
  KeggGSEA<-KeggGSEA[KeggGSEA$minNES<Inf,]
  KeggGSEA$metab_leadingEdge<-
    sapply(KeggGSEA$metab_leadingEdge,function(x) paste(x,collapse=";"))
  KeggGSEA$gene_leadingEdge<-
    sapply(KeggGSEA$gene_leadingEdge,function(x) paste(x,collapse=";"))
  write.csv(KeggGSEA,row.names=FALSE,
            file=paste0("~/gdrive/BearOmics2/EnrichmentAnalysis/Tables/",comparison,".csv"))
  
  ########### Pathview ###########
  pathList<-unique(KeggGSEA$pathId)
  setwd("~/gdrive/BearOmics2/EnrichmentAnalysis/KEGGplots/")
  badPath<-c()
  for(path in pathList){
    # Metabolite data:
    path_Metabs<-diffs[diffs$kegg %in% keggMetabSet[[paste0("map",path)]],]
    cpdData<-path_Metabs[,paste0("FC_",comparison)]
    names(cpdData)<-path_Metabs$kegg
    
    # Gene data:
    path_Genes<-geneDf[geneDf$ENTREZ.ID %in% keggGeneSet[[paste0("path:hsa",path)]],]
    geneData<-as.numeric(path_Genes$logFC)
    names(geneData)<-as.character(path_Genes$ENTREZ.ID)
    
    tryCatch({
      pv.out<-pathview(gene.data=geneData,cpd.data=cpdData,
                       pathway.id=path,species="hsa",out.suffix=paste0(comparison),
                       keys.align="y",kegg.native=T,key.pos="topright",
                       limit=list(gene=c(-2,2),cpd=c(-2,2)),bins=list(gene=14,cpd=14),
                       high=list(gene="green",cpd="blue"),mid=list(gene="grey50",cpd="grey50"),
                       low=list(gene="red",cpd="#FFB533"))
    },
    error=function(e){
      badPath<-c(badPath,path)
      return(badPath)
    }
    )
  }
}

# metabolism example: Lysine degredation map00310
# plotEnrichment(keggMetabSet[["map00310"]],metabStats)
# plotEnrichment(keggGeneSet[["path:hsa00310"]],geneStats)
# path00310_Metabs<-diffs[diffs$kegg %in% keggMetabSet[["map00310"]],]
# path00310_Genes<-S_Up[S_Up$ENTREZ.ID %in% keggGeneSet[["path:hsa00310"]],]
# write.table(path00310_Metabs[,c("kegg","SvsUFC")],file="path00310_Metabs.txt",
#             row.names=FALSE,sep="\t")
# write.table(path00310_Genes[,c("ENTREZ.ID","logFC")],file="path00310_Genes.txt",
#             row.names=FALSE,sep="\t")

########### Ex Pathview ###########
geneData_path00310<-as.numeric(path00310_Genes$logFC)
names(geneData_path00310)<-as.character(path00310_Genes$ENTREZ.ID)
cpdData_path00310<-as.numeric(path00310_Metabs$SvsUFC)
names(cpdData_path00310)<-as.character(path00310_Metabs$kegg)
pv.out<-pathview(gene.data=geneData_path00310,cpd.data=cpdData_path00310,
            pathway.id="00310",species="hsa",out.suffix="mysuf",
            keys.align="y",kegg.native=T,key.pos="topright",
            limit=list(gene=c(-2,2),cpd=c(-2,2)),bins=list(gene=14,cpd=14),
            high=list(gene="green",cpd="blue"),mid=list(gene="grey29",cpd="grey29"),
            low=list(gene="red",cpd="#FFB533"))

########### Add / fix ChEBIs ###########
metabKey2<-metabKey

# Finished 1-30:


########### fgsea with custom set from Reactome ###########
diffs$ChEBI2<-sapply(strsplit(diffs$ChEBI,";"),function(x) x[1])

# Make list of sets
chebiReact$ChEBI<-as.character(chebiReact$ChEBI)
metabSet<-list()
uniquePaths<-unique(chebiReact$pathID)
for(path in uniquePaths){
  metabSet[[path]]<-unique(chebiReact$ChEBI[chebiReact$pathID==path])
}
metabStats<-diffs$SvsU
names(metabStats)<-diffs$ChEBI2
metabStats<-metabStats[!is.na(names(metabStats))]

metabGSEA<-fgsea(metabSet,metabStats,nperm = 1000)

########### Example GSEA using ReactomePA ###########
# Over representation (hypergeometric):
S_Up3<-S_Up %>% filter(p_value<.05)
over_S_Up<-as.data.frame(enrichPathway(S_Up3$ENTREZ.ID,pvalueCutoff=.2,readable=TRUE))

# Actual GSEA
gse_S_Up<-gsePathway(geneList=S_Up2b,nPerm=1000,minGSSize=10,pvalueCutoff=0.25,
           pAdjustMethod="BH",verbose=TRUE)
gseSum_S_Up<-as.data.frame(gse_S_Up)
gseaplot(gse_S_Up,geneSetID="R-HSA-174824")

# Need to see fgsea::fgsea