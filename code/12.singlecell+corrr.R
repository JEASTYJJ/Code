library(Seurat)
library(ggplot2)
library(dplyr)
library(tidyr)
library(scDblFinder)
library(reshape2)
library(clustree)
scRNA<-readRDS("./GSE152805/GSE152805_human_chondrocytes_anno.rds")

Idents(scRNA)<-"orig.ident"
ggsave("./GSE152805/QC.pdf",width = 15,height = 7)
ncol(scRNA)
scRNA <- FindVariableFeatures(scRNA, selection.method = "vst", nfeatures = 1000)
top10 <- head(VariableFeatures(scRNA), 10)
LabelPoints(plot = VariableFeaturePlot(scRNA), 
            points = top10, repel = T,cex=2)+
  theme(legend.key.height = unit(0.25,"cm"),legend.key.width = unit(0.1,"cm"),
        legend.text = element_text(size = 6),legend.title = element_blank(),
        legend.position = "bottom",legend.direction = "vertical",
        text = element_text(size = 7),
        axis.text = element_text(size = 6),
        axis.ticks = element_line(size = 0.75)
  )
ggsave("./GSE152805/VariableFeatures_dot.pdf",width = 6,height = 9,units = "cm")
FeatureScatter(scRNA, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
ggsave("./GSE152805/nFeature_nCount_cor_dot.pdf",width = 12,height = 9,units = "cm")
scRNA <- ScaleData(scRNA, features = rownames(scRNA))
scRNA <- RunPCA(scRNA, features = VariableFeatures(object = scRNA),npcs = 50)
ElbowPlot(scRNA,ndims = 50)+
  theme(text = element_text(size = 7),
        axis.text = element_text(size = 6),
        axis.ticks = element_line(size = 0.75))
ggsave("./GSE152805/ElbowPlot.pdf",width = 9,height = 6,units = "cm")
pca.num = 1:18 
scRNA <- FindNeighbors(scRNA, dims = pca.num)
scRNA <- FindClusters(scRNA, resolution = c(seq(0,1.6,.2)))
clustree(scRNA@meta.data, prefix = "RNA_snn_res.")
scRNA <- FindClusters(scRNA, resolution = 1.2)
scRNA <- RunUMAP(scRNA, dims = pca.num)
scRNA <- RunTSNE(scRNA, dims = pca.num,check_duplicates = F)

Idents(scRNA) <- "seurat_clusters"
DimPlot(scRNA, label = T,label.size = 2,
        reduction = "umap",group.by = "seurat_clusters",
        repel = T)+
  theme(plot.title = element_blank(),
        axis.line = element_line(size = 0.75),
        axis.title = element_text(size = 6,hjust = 0), 
        axis.text = element_blank(), 
        axis.ticks = element_blank(), 
        text = element_text(size = 8),
        legend.key.height = unit(0.25,"cm"),legend.key.width = unit(0.1,"cm"),
        legend.text = element_text(size = 6),legend.title = element_blank())
ggsave("./GSE152805/seurat_clusters_umap.pdf",width = 8,height = 6,units = "cm")

check_gene<-c('COL10A1','IBSP','COL2A1',
              'COL1A1','COL1A2','S100A4','PRG4',
              'COL10A1','IBSP','JUN',
              "MMP3",'FOSB',"JUN",
              'IL11','COL2A1','CILP','OGN',
              'CHI3L1','CHI3L2',
              'COL2A1','CILP','COL3A1','COMP'
)

DotPlot(scRNA, features = unique(check_gene),group.by = "seurat_clusters")+
  scale_x_discrete("")+scale_y_discrete("")+ coord_flip()
ggsave("./GSE152805/check_gene_dot.pdf",width = 15,height = 8)

scRNA<-readRDS("./GSE152805/GSE152805_anno.rds")

DimPlot(scRNA, label = T,label.size = 2,
        reduction = "umap",group.by = "celltype",
        repel = T)+
  theme(plot.title = element_blank(),
        axis.line = element_line(size = 0.75),
        axis.title = element_text(size = 6,hjust = 0), 
        axis.text = element_blank(), 
        axis.ticks = element_blank(), 
        text = element_text(size = 8),
        legend.key.height = unit(0.25,"cm"),legend.key.width = unit(0.1,"cm"),
        legend.text = element_text(size = 8),legend.title = element_blank())
ggsave("./GSE152805/celltype_umap.pdf",width = 8,height = 6,units = "cm")

Idents(scRNA) <- "celltype"
sce.markers1 <- FindAllMarkers(scRNA, logfc.threshold = 0.25, min.pct = 0.1, only.pos = F)
write.csv(sce.markers1,"./GSE152805/DEG_per_celltype.csv")
sce.markers1<-read.csv("./GSE152805/DEG_per_celltype.csv",row.names = 1)

Idents(scRNA) <- "seurat_clusters"
sce.markers2 <- FindAllMarkers(scRNA, logfc.threshold = 0.25, min.pct = 0.1, only.pos = F)
write.csv(sce.markers2,"./GSE152805/DEG_per_cluster.csv")
sce.markers2<-read.csv("./GSE152805/DEG_per_cluster.csv",row.names = 1)

top2 <- sce.markers1 %>% group_by(cluster) %>% top_n(2,wt = avg_log2FC)

Idents(scRNA) <- "celltype"
scRNA2<-subset(x=scRNA,downsample = 200)


a<-DoHeatmap(scRNA2, features = top2$gene,group.by = "celltype",angle = 45,raster = F,size = 10,label = F,combine = T)
a+theme(axis.text = element_text(size = 15),
        legend.text = element_text(size = 10))
ggsave("./GSE152805/celltype_DEG_heatmap.pdf",width = 10,height = 10)

Idents(scRNA2) <- "celltype"
check_gene<-c('IBSP',#HTC
              "JUN",#HOMC
              'TGFBI',#preHTC
              'COL2A1',#RepC
              'CHI3L1',#RegC
              'COL1A2',#FC
              'IL11'#preFC
)

for(i in check_gene){
  f=paste0("./GSE152805/",i,"_umap.pdf")
  a<-FeaturePlot(scRNA,features = i,order = T,reduction = "umap",pt.size = 0.1,ncol = 1,combine = T,
                 cols = c("#0072B5","#FFFFB3","#E31F1F")
  )
  a+theme(  axis.line = element_blank(),
            axis.title = element_blank(), 
            axis.text = element_blank(), 
            axis.ticks = element_blank(), 
            text = element_text(size = 20))
  ggsave(f,width = 8,height = 6)
}

scRNA$class1<-ifelse(scRNA$class=='oLT','Normal','OA')

DimPlot(scRNA, label = T,label.size = 2,
        reduction = "umap",group.by = "class1",
        repel = T)+
  theme(plot.title = element_blank(),
        axis.line = element_line(size = 0.75),
        axis.title = element_text(size = 6,hjust = 0), 
        axis.text = element_blank(), 
        axis.ticks = element_blank(), 
        text = element_text(size = 8),
        legend.key.height = unit(0.25,"cm"),legend.key.width = unit(0.1,"cm"),
        legend.text = element_text(size = 8),legend.title = element_blank())
ggsave("./GSE152805/class_umap.pdf",width = 8,height = 6,units = "cm")

dat <- c()
m <- scRNA@meta.data
dat$class <- m$class1
dat$celltype <- m$celltype
dat0=melt(data.frame(table(dat)))
ggplot(dat0,aes(value,class,fill=celltype))+
  geom_bar(stat="identity",position="fill")+ 
  xlab("") + ylab("")+ labs(fill = "celltype")+
  theme_classic()+geom_text_repel(stat = "identity",position="fill",aes(label=value))
ggsave('./GSE152805/class_celltype_freq_bar.pdf',width = 8,height = 6)


library(ggrepel)
Idents(scRNA)<-'class1'
scRNA3<-subset(scRNA,idents='OA')
Idents(scRNA) <- "celltype"
Idents(scRNA3) <- "celltype"
blank_theme <- theme_minimal()+
  theme(axis.title.x = element_blank(),axis.title.y = element_blank(),axis.text.x = element_blank(),axis.text.y = element_blank(),panel.border = element_blank(),panel.grid=element_blank(),axis.ticks = element_blank(),
        plot.title=element_text(size=14, face="bold")
  )
ggplot(data=scRNA3@meta.data,         mapping=aes(x="celltype",fill=celltype))+
  geom_bar(stat="count",width=0.5,position='stack',size=5)+
  labs(fill = "celltype")+
  coord_polar("y", start=0)+
  blank_theme +
  geom_text_repel(stat="count",aes(label = ..count..), size=6, position=position_stack(vjust = 0.5))
ggsave("./GSE152805/Crohn_celltype_number_pie.pdf",width = 8,height = 7)


OA_lasso_gene<-read.csv("./GSE55235/OA_lasso_gene.csv",header = T,row.names = 1)[,1]
O_s_inter<-intersect(sce.markers1$gene,OA_lasso_gene)
write.csv(O_s_inter,'')
Idents(scRNA) <- "celltype"
VlnPlot(scRNA,features = O_s_inter,stack = T, pt.size = 0,flip = T)
ggsave("./GSE152805/lasso_DEG_vln.pdf",width = 8,height = 8)

s_exp <- GetAssayData(scRNA,slot = "scale.data")
s_texp <- as.data.frame(t(s_exp[O_s_inter,]))
s_texp[1:4,1:4]
s_texp$Normal<-ifelse(scRNA$class1=='Normal',1,0)
s_texp$OA<-ifelse(scRNA$class1=='OA',1,0)
library(corrr)

s_x<-cor(s_texp, use='everything', method='pearson') 
s_x1<-s_x[1:4,5:6]
pdf("./GSE152805/DEG_cor_heatmap_plot.pdf",width = 8,height = 8)
corrplot::corrplot(s_x1,mar = c(0, 2, 2, 2),tl.srt = 0,tl.col = 'black',tl.cex = 2,cl.ratio = 0.5,cl.cex = 1.5)
dev.off()

write.csv(as.data.frame(s_x),"./GSE152805/DEG_corr.csv",row.names = T)

