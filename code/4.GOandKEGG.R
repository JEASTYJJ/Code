library(clusterProfiler)
library(ggpubr)

eg <- bitr(cygene, 
           fromType="SYMBOL", 
           toType=c("ENTREZID","ENSEMBL",'SYMBOL'),
           OrgDb="org.Hs.eg.db")
head(eg)

go <- enrichGO(gene = eg$ENTREZID, 
               OrgDb = org.Hs.eg.db, 
               ont='ALL',
               pAdjustMethod = 'BH',
               pvalueCutoff = 0.05, 
               qvalueCutoff = 0.05,
               readable = T,
               keyType = 'ENTREZID')
g<-as.data.frame(go)
write.csv(g,"OA_RA_DEG_gene_GO.csv")
g$pvalue <- -log10(g$pvalue)
g1=g %>% group_by(ONTOLOGY) %>% top_n(3,wt = pvalue)
ggdotchart(g1, x="Description", y="pvalue", color = "ONTOLOGY",
           sorting = "descending",   
           add = "segments",    
           xlab = 'GO',ylab = "-log10(pvalue)",
           rotate = T,       
           dot.size = 7,        
           ggtheme = theme_pubr()+
             theme(axis.text.y = element_text(size = 12))
)+scale_x_discrete(labels=function(x) str_wrap(x, width=35))

ggsave("OA_RA_DEG_gene_GO.pdf",height = 8,width = 6)

g2<-g1
g2$Count<-as.numeric(g2$Count)

library(RColorBrewer)
col <- brewer.pal(9,'Set1')

ggplot(g2,aes(Count,reorder(Description,Count)))+
  geom_col(aes(fill=pvalue),width = 0.6)+
  scale_fill_gradient2(
    low = col[2],
    mid = col[1],
    high = col[1],
    midpoint = mean(g2$pvalue)
  )+
  theme_bw()+
  scale_y_discrete(labels=function(x) str_wrap(x, width=35))+
  theme(axis.title=element_text(size=15,color="black"),
        axis.text.y=element_text(size=12,color="black"),
        axis.text.x=element_text(size=12,color="black"),
        legend.title=element_text(size=13,color="black"),
        legend.text=element_text(size=10,color="black"))+
  labs(x='Gene number',y='',fill='-log10(pvalue)')

ggsave("OA_RA_DEG_gene_GO_bar.pdf",height = 8,width = 8)

kegg <- enrichKEGG(eg$ENTREZID,
                   organism = 'hsa',
                   keyType = 'kegg', 
                   pvalueCutoff = 0.05,
                   pAdjustMethod = 'BH', 
                   qvalueCutoff = 0.05)

e<-as.data.frame(kegg)
write.csv(e,"OA_RA_DEG_gene_KEGG.csv")
e$pvalue <- -log10(e$pvalue)
ggdotchart(e[1:10,], x="Description", y="pvalue", color = "Description",
           palette = colorRampPalette(colors = c("#4DBBD5", "#E64B35", "#00A087"))(10),     
           sorting = "descending",   
           add = "segments",    
           add.params = list(color = colorRampPalette(colors = c("#4DBBD5", "#E64B35", "#00A087"))(10), size = 1),
           xlab = 'KEGG',ylab = "-log10(pvalue)",
           rotate = T,       
           dot.size = 7,        
           ggtheme = theme_pubr()+
             theme(axis.text.y = element_text(size = 12))
)+NoLegend()+scale_x_discrete(labels=function(x) str_wrap(x, width=35))

ggsave("OA_RA_DEG_gene_KEGG.pdf",height = 8,width = 6)


e2<-e[1:10,]
e2$Count<-as.numeric(e2$Count)

ggplot(e2,aes(Count,reorder(Description,Count)))+
  geom_col(aes(fill=pvalue),width = 0.6)+
  scale_fill_gradient2(
    low = col[2],
    mid = col[2],
    high = col[1],
    midpoint = mean(e2$pvalue)
  )+
  theme_bw()+
  scale_y_discrete(labels=function(x) str_wrap(x, width=35))+
  theme(axis.title=element_text(size=15,color="black"),
        axis.text.y=element_text(size=12,color="black"),
        axis.text.x=element_text(size=12,color="black"),
        legend.title=element_text(size=13,color="black"),
        legend.text=element_text(size=10,color="black"))+
  labs(x='Gene number',y='',fill='-log10(pvalue)')

ggsave("OA_RA_DEG_gene_KEGG_bar.pdf",height = 8,width = 6)
