library(corrr)

texp <- as.data.frame(t(exp1[cygene,]))
x<-cor(texp, use='everything', method='pearson') 
write.csv(as.data.frame(x),"OA_RA_DEG_corr.csv",row.names = T)

network_plot(x,min_cor = .7,repel = T,curved = T)
ggsave("OA_RA_network_plot.pdf",width = 6,height = 6)


O_exp2<-read.csv("OA_class.csv",header = T,row.names = 1)
texp1<-O_exp2[,hub_gene]
texp1$Cluster1<-ifelse(O_exp2$class==1,1,0)
texp1$Cluster2<-ifelse(O_exp2$class==2,1,0)
x1<-cor(texp1, use='everything', method='pearson') 
write.csv(as.data.frame(x1),"./OA_class_hub_corr.csv",row.names = T)
x1_1<-x1[7:8,1:6]
pdf('./OA_class_hub_corr.pdf',width = 6,height = 3)
pheatmap::pheatmap(x1_1,fontsize = 20,cluster_rows = F,cluster_cols = F)
dev.off()
R_exp2<-read.csv("RA_class.csv",header = T,row.names = 1)
texp2<-R_exp2[,hub_gene]
texp2$Cluster1<-ifelse(R_exp2$class==1,1,0)
texp2$Cluster2<-ifelse(R_exp2$class==2,1,0)
x2<-cor(texp2, use='everything', method='pearson') 
write.csv(as.data.frame(x2),"./RA_class_hub_corr.csv",row.names = T)
x2_1<-x2[7:8,1:6]
pdf('./RA_class_hub_corr.pdf',width = 6,height = 3)
pheatmap::pheatmap(x2_1,fontsize = 20,cluster_rows = F,cluster_cols = F)
dev.off()

obj<-read.csv("CIBERSORT_Results.csv",row.names = 1)
texp3<-obj[,1:22]
exp1_1<-as.data.frame(t(exp1[hub_gene,]))
texp3<-cbind(texp3,exp1_1)

x3<-cor(texp3, use='everything', method='pearson') 
write.csv(as.data.frame(x3),"./CIBERSORT_hub_corr.csv",row.names = T)
x3_1<-x3[23:28,1:22]
pdf('./CIBERSORT_hub_corr.pdf',width = 8,height = 8)
pheatmap::pheatmap(x3_1,fontsize = 20,cluster_rows = F,cluster_cols = F)
dev.off()

