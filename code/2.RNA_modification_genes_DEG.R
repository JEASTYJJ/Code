
exp<-read.csv("Total_exp.csv",row.names = 1)
pd<-read.csv("Total_pd.csv",row.names = 1)

group_list<-ifelse(str_detect(pd$type,"control"),"Control",ifelse(str_detect(pd$type,"rheumatoid"),"RA","OA"))
table(group_list)#Control 20,OA 26,RA 33

M1A<-read.table("./RMBase_hg19_all_m1A_site.txt")
M6A<-read.table("./RMBase_hg19_all_m6A_site.txt")
M5C<-read.csv("./Homo_sapiens.csv",row.names = 1)
m1A<-M1A$V12
m6A<-M6A$V11
m5C<-M5C$geneId

m1A1<-unlist(str_split(m1A,pattern = ","))
m1A_gene<-unique(intersect(rownames(exp),m1A1))
m6A1<-unlist(str_split(m6A,pattern = ","))
m6A_gene<-unique(intersect(rownames(exp),m6A1))
m5C_gene<-unique(intersect(rownames(exp),m5C))
f_gene<-union(m1A_gene,m6A_gene)
f_gene<-union(f_gene,m5C_gene)
write.csv(f_gene,"RNA_modify_genes.csv")

library(limma)

group1<-group_list[group_list!="Control"]
exp1<-exp[f_gene,group_list!="Control"]

group<-factor(group1,levels = c("OA","RA"))

design=model.matrix(~group)
fit=lmFit(exp1,design)
fit=eBayes(fit)
deg=topTable(fit,coef = 2,number = Inf)
logFC=1
P.Value = 0.05
type1 = (deg$P.Value < P.Value)&(deg$logFC < -logFC)
type2 = (deg$P.Value < P.Value)&(deg$logFC > logFC)
deg$change = ifelse(type1,"down",ifelse(type2,"up","stable"))
table(deg$change)
cygene = rownames(deg)[deg$change !="stable"]
deg$change1<-ifelse(deg$change=="up","RA",ifelse(deg$change=="down","OA","Not sig"))
write.csv(deg,"./OA_RA_DEG.csv")
write.csv(cygene,"OA_RA_DEG_genes.csv")

library(ggplot2)
p1 <- ggplot(data = deg, 
             aes(x = logFC, 
                 y = -log10(P.Value))) +
  geom_point(alpha=0.4, size=3.5, 
             aes(color=change1)) +
  ylab("-log10(Pvalue)")+
  labs(color="Type")+
  scale_color_manual(values=c("grey","orangered2","blue2"))+
  geom_vline(xintercept=c(-logFC,logFC),lty=4,col=c("orangered2","blue2"),lwd=0.8) +
  geom_hline(yintercept = -log10(P.Value),lty=4,col="black",lwd=0.8) +
  theme_bw()

p1
ggsave("./OA_RA_DEG_vol.pdf",height = 6,width = 6)

diff=exp1[cygene,]
library(pheatmap)
annotation_col=data.frame(group=group1)
rownames(annotation_col)=colnames(diff) 

pdf("OA_RA_DEG_heatmap.pdf",height = 6,width = 6)

pheatmap(diff,
         annotation_col=annotation_col,
         scale = "row",
         show_rownames = F,
         show_colnames =F,
         color = colorRampPalette(c("navy", "white", "red"))(50),
         fontsize = 10,
         fontsize_row=3,
         fontsize_col=3)
dev.off()


deg1<-deg[deg$change!="stable",]
deg1<-deg1[order(deg1$logFC),]

texp1 <- as.data.frame(t(exp1[rownames(deg1)[c(1:5,136:140)],]))
boxdata<-texp1%>%
  gather(key = gene)
boxdata$type<-group1
boxdata$gene<-factor(boxdata$gene,levels = rownames(deg1)[c(1:5,136:140)])

p2<-ggplot(boxdata,aes(x=gene,y=value,fill=type))+
  geom_boxplot()+
  theme_bw()+
  labs(x="",y="")+
  theme(axis.text = element_text(size = 15),
        legend.text = element_text(size = 15),
        legend.title = element_text(size = 15))+
  stat_compare_means(aes(label = ..p.signif..))+
  scale_fill_manual(values =c("orangered2","dodgerblue2"))


p2
ggsave("OA_RA_DEG_boxplot.pdf",width = 15,height = 8)

ppi<-read.table(gzfile("9606.protein.links.v11.5.txt.gz"),header = T)
info<-read.table(gzfile("9606.protein.info.v11.5.txt.gz"),sep = "\t")
ppi1<-ppi[ppi[,3]>700,]
ppi1[,1]<-str_sub(ppi1[,1],6)
ppi1[,2]<-str_sub(ppi1[,2],6)
rownames(info)<-str_sub(info[,1],6)
ppi1[,1]<-info[ppi1[,1],2]
ppi1[,2]<-info[ppi1[,2],2]
ppi1<-na.omit(ppi1)
write.table(ppi1,"human_PPI>700.txt",row.names = F,sep = "\t")
write.table(cygene,"OA_RA_DEG_gene.txt",row.names = F,sep = '\t')
a<-deg[deg$change!="stable",]
b<-cbind(rownames(a),a$change)
write.table(b,"OA_RA_DEG_gene_info.txt",row.names = F,sep = '\t')
