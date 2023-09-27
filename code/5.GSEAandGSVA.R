library(msigdbr)
library(GSVA)


min.sz = 2 
max.sz = 10000 
parallel.sz= 10 
mx.diff= T 
tau=1  
method='gsva' 
kcdf = "Gaussian" 

msgdC2 = msigdbr(species = "Homo sapiens", category = "C2", subcategory = "KEGG")
gene_sets = msgdC2 %>% split(x = .$gene_symbol, f = .$gs_name)
exp2<-as.matrix(exp1)
gsva_scores = gsva(expr=exp2, gset.idx.list=gene_sets, 
                   method=method, kcdf=kcdf, min.sz=min.sz, max.sz=max.sz, 
                   parallel.sz=parallel.sz, mx.diff=mx.diff )


gsva_scores = as.data.frame(gsva_scores)
gsva_scores$geneset = rownames(gsva_scores)
gsva_scores = gsva_scores %>% dplyr::select(geneset,everything())
write.csv(gsva_scores,file = "GSVA_enrichment_results.csv", row.names = F,quote = FALSE)



library(limma)

group1<-group_list[group_list!="Control"]
gsva_scores1<-as.data.frame(gsva_scores[,-1])

group<-factor(group1,levels = c("OA","RA"))

design=model.matrix(~group)
fit=lmFit(gsva_scores1,design)
fit=eBayes(fit)
deg2=topTable(fit,coef = 2,number = Inf)
P.Value = 0.05
deg2<-deg2[deg2$P.Value<P.Value,]
deg2<-deg2[order(deg2$logFC),]
cypathway = rownames(deg2)[c(1:5,(nrow(deg2)-4):nrow(deg2))]

write.csv(deg2,"OA_RA_GSVA_pathway.csv")
write.csv(cypathway,"OA_RA_GSVA_pathway_heatmap.csv")
diff1=gsva_scores1[cypathway,]
library(pheatmap)
annotation_col=data.frame(group=group1)
rownames(annotation_col)=colnames(diff1) 

pdf("OA_RA_GSVA_pathway_heatmap.pdf",height = 5,width = 10)

pheatmap(diff1,
         annotation_col=annotation_col,
         scale = "row",
         cluster_rows = F,treeheight_col = 0,
         show_rownames = T,
         show_colnames =F,
         color = colorRampPalette(c("navy", "white", "red"))(50),
         fontsize = 10,
         fontsize_row=15,
         fontsize_col=3)
dev.off()


library(GSEABase)
library(enrichplot)
library(cowplot)
deg<-read.csv("OA_RA_DEG.csv",row.names = 1)
ge = deg$logFC
names(ge) = rownames(deg)
ge = sort(ge,decreasing = T)
head(ge)

pathway_gene<-data_frame(term="1",gene="1")
for(i in 1:length(gene_sets)){
  a<-data.frame(term=names(gene_sets[i]),gene=unlist(gene_sets[[i]]))
  pathway_gene<-pathway_gene%>%rbind(a)
}
pathway_gene<-pathway_gene[-1,]

OA_RA_GSEA<-GSEA(ge, TERM2GENE =pathway_gene,
                 nPerm=1000,
                 verbose=FALSE,by="fgsea",
                 pAdjustMethod="BH",
                 pvalueCutoff=0.05)


res<-OA_RA_GSEA@result
write.csv(res,"OA_RA_GSEA.csv")

a<-order(OA_RA_GSEA$NES,decreasing = T)[c(1:3,(nrow(res)-2):nrow(res))]
titles =str_sub(OA_RA_GSEA$Description[a],6)
titles =gsub("_"," ",titles)
gseaplot<-list()
for (i in 1:length(a)) {
  gseaplot[[i]]<-gseaplot2(OA_RA_GSEA, geneSetID = a[i], title =titles[i],pvalue_table = F,base_size = 9,subplots = 1:3)+
    annotate("text", x = 9000, y = 2, label = paste0('NES = ',round(res[a[i],5],4),'\n',
                                                     "P.adjust = ",round(res[a[i],7],4),'\n',
                                                     'FDR = ',round(res[a[i],8],4)))
}


OA_RA_GSEA_gseplot<-plot_grid(gseaplot[[1]],gseaplot[[2]],gseaplot[[3]],
                              gseaplot[[4]],gseaplot[[5]],gseaplot[[6]],
                              nrow=2)
ggsave("OA_RA_GSEA_gseplot.pdf",width =17,height =14)

