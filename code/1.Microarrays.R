
library(affy)
library(GEOquery)
library(stringr)
library(dplyr)
library(tidyr)


untar("GSE55235_RAW.tar",exdir = "GSE55235_RAW")
setwd("./GSE55235_RAW")

gse<-getGEO('GSE55235')
phenotype<- pData(phenoData(gse[[1]]))
prob2gene<-fData(gse[[1]])
prob<-prob2gene[,c("ID","Gene Symbol")]
prob1<-prob%>%
  filter(`Gene Symbol` !="")%>%
  separate(`Gene Symbol`,c("drop","symbol"),sep = "///")


Data <- ReadAffy()
GSE55235<-rma(Data)
normalize.methods(Data)
GSE55235<-justRMA()

for(i in 1:nrow(as.matrix((sampleNames(GSE55235)))))
{
  sampleNames(GSE55235)[i]<-sub(".CEL.gz$", "",sampleNames(GSE55235)[i])
}


exprSet<-exprs(GSE55235)

{
  exprSet=exprSet[rownames(exprSet)%in%prob1[,1],]
  ID2gene=prob1[match(rownames(exprSet),prob1[,1]),]
}


dim(exprSet)
dim(ID2gene)

write.csv(phenotype,"GSE55235_pd.csv")

{
  MAX=by(exprSet,ID2gene[,2],
         function(x)rownames(x)[which.max(rowMeans(x))])
  MAX=as.character(MAX)
  exprSet=exprSet[rownames(exprSet)%in%MAX,]
  rownames(exprSet)=ID2gene[match(rownames(exprSet),ID2gene[,1]),2]
}
dim(exprSet)
exprSet[1:5,1:5]

write.csv(exprSet,"GSE55235_exp.csv")

setwd("../")

untar("GSE55457_RAW.tar",exdir = "GSE55457_RAW")
setwd("./GSE55457_RAW")

gse<-getGEO('GSE55457')
phenotype<- pData(phenoData(gse[[1]]))
prob2gene<-fData(gse[[1]])
prob<-prob2gene[,c("ID","Gene Symbol")]
prob1<-prob%>%
  filter(`Gene Symbol` !="")%>%
  separate(`Gene Symbol`,c("drop","symbol"),sep = "///")


Data <- ReadAffy()
GSE55457<-rma(Data)
normalize.methods(Data)
GSE55457<-justRMA()

for(i in 1:nrow(as.matrix((sampleNames(GSE55457)))))
{
  sampleNames(GSE55457)[i]<-sub(".CEL.gz$", "",sampleNames(GSE55457)[i])
}


exprSet<-exprs(GSE55457)

{
  exprSet=exprSet[rownames(exprSet)%in%prob1[,1],]
  ID2gene=prob1[match(rownames(exprSet),prob1[,1]),]
}


dim(exprSet)
dim(ID2gene)

write.csv(phenotype,"GSE55457_pd.csv")

{
  MAX=by(exprSet,ID2gene[,2],
         function(x)rownames(x)[which.max(rowMeans(x))])
  MAX=as.character(MAX)
  exprSet=exprSet[rownames(exprSet)%in%MAX,]
  rownames(exprSet)=ID2gene[match(rownames(exprSet),ID2gene[,1]),2]
}
dim(exprSet)
exprSet[1:5,1:5]

write.csv(exprSet,"GSE55457_exp.csv")

setwd("../")

untar("GSE55584_RAW.tar",exdir = "GSE55584_RAW")
setwd("./GSE55584_RAW")

gse<-getGEO('GSE55584')
phenotype<- pData(phenoData(gse[[1]]))
prob2gene<-fData(gse[[1]])
prob<-prob2gene[,c("ID","Gene Symbol")]
prob1<-prob%>%
  filter(`Gene Symbol` !="")%>%
  separate(`Gene Symbol`,c("drop","symbol"),sep = "///")


Data <- ReadAffy()
GSE55584<-rma(Data)
normalize.methods(Data)
GSE55584<-justRMA()

for(i in 1:nrow(as.matrix((sampleNames(GSE55584)))))
{
  sampleNames(GSE55584)[i]<-sub(".CEL.gz$", "",sampleNames(GSE55584)[i])
}


exprSet<-exprs(GSE55584)

{
  exprSet=exprSet[rownames(exprSet)%in%prob1[,1],]
  ID2gene=prob1[match(rownames(exprSet),prob1[,1]),]
}


dim(exprSet)
dim(ID2gene)

write.csv(phenotype,"GSE55584_pd.csv")

{
  MAX=by(exprSet,ID2gene[,2],
         function(x)rownames(x)[which.max(rowMeans(x))])
  MAX=as.character(MAX)
  exprSet=exprSet[rownames(exprSet)%in%MAX,]
  rownames(exprSet)=ID2gene[match(rownames(exprSet),ID2gene[,1]),2]
}
dim(exprSet)
exprSet[1:5,1:5]

write.csv(exprSet,"GSE55584_exp.csv")

setwd("../")

GSE55584_exp<-read.csv("./GSE55584_RAW/GSE55584_exp.csv")
GSE55235_exp<-read.csv("./GSE55235_RAW/GSE55235_exp.csv")
GSE55457_exp<-read.csv("./GSE55457_RAW/GSE55457_exp.csv")

a<-merge(x=GSE55235_exp,y=GSE55457_exp,by="X")
b<-merge(x=a,y=GSE55584_exp,by="X")
rownames(b)<-b$X
exprSet<-b[,-1]

GSE55584_pd<-read.csv("./GSE55584_RAW/GSE55584_pd.csv")
GSE55235_pd<-read.csv("./GSE55235_RAW/GSE55235_pd.csv")
GSE55457_pd<-read.csv("./GSE55457_RAW/GSE55457_pd.csv")

pd1<-data.frame(GSE55235_pd[,c('X','characteristics_ch1')],bath=1)
pd2<-data.frame(GSE55457_pd[,c('X','characteristics_ch1.2')],bath=2)
pd3<-data.frame(GSE55584_pd[,c('X','characteristics_ch1.2')],bath=3)
colnames(pd1)<-c('X','characteristics_ch1.2','bath')
pd<-rbind(pd1,pd2,pd3)

colnames(pd)<-c('sample','type','batch')

head(colnames(exprSet))
colnames(exprSet)<-substr(colnames(exprSet),1,10)
exprSet<-exprSet[,pd$sample]

boxplot(exprSet)

library(sva)
batch = pd$batch
combat_edata1 = ComBat(dat=exprSet, batch=batch, mod=NULL, par.prior=TRUE,  prior.plots=FALSE)
boxplot(combat_edata1)

pd$batch<-ifelse(pd$batch==1,"GSE55235",ifelse(pd$batch==2,"GSE55457","GSE55584"))
write.csv(combat_edata1,"Total_exp.csv")
write.csv(pd,"Total_pd.csv")

