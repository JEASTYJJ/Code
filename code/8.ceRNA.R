library(readxl)
hub_gene<-read.csv("./mod_DEG_PPI_top6.csv",header = T)[2:7,2]


mRNA_miRNA2<-read.table("starBase_Human_Pan-Cancer_MiRNA-Target_Interactions2022-07-25_02-37.xls",header = T)
mRNA_miRNA3<-mRNA_miRNA2[mRNA_miRNA2$geneName%in%hub_gene,]
miRNA<-unique(mRNA_miRNA3$name)

miRNA_lncRNA1<-read.table("starBase_Human_Pan-Cancer_miRNA-LncRNA_Interactions2022-07-25_01-55.xls",header = T)
miRNA_lncRNA2<-miRNA_lncRNA1[miRNA_lncRNA1$name%in%miRNA,]
a<-as.data.frame(table(miRNA_lncRNA2$geneName))
b<-a[a$Freq>10,]
miRNA_lncRNA3<-miRNA_lncRNA2[miRNA_lncRNA2$geneName%in%b$Var1,]

ceRNA<-mRNA_miRNA3[,c(1,2)]
ceRNA<-rbind(ceRNA,miRNA_lncRNA3[,c(1,3)])
ceRNA$type<-c(rep('mRNA',times=nrow(mRNA_miRNA3)),rep('lncRNA',times=nrow(ceRNA)-nrow(mRNA_miRNA3)))
a<-as.data.frame(table(ceRNA$name))
b<-a[a$Freq>1,]
ceRNA1<-ceRNA[ceRNA$name%in%b$Var1,]

write.table(ceRNA1,"./ceRNA.txt",row.names = F,sep = '\t')
