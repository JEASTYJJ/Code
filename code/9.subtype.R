
O_exp1<-as.matrix(O_exp[OA_lasso_gene,O_group!='Control'])
library(ConsensusClusterPlus)

title="./OA_ConsensusCluster"
OA_results = ConsensusClusterPlus(O_exp1,maxK=10,reps=1000,pItem=0.8,pFeature=1,
                                  title=title,clusterAlg="hc",distance="euclidean",seed=2333,plot="pdf")

O_exp2<-as.data.frame(t(O_exp[,O_group!='Control']))
O_exp2$class<-OA_results[[2]]$consensusClass

write.csv(O_exp2,"OA_class.csv")
R_exp1<-as.matrix(R_exp[RA_lasso_gene,R_group!='Control'])
library(ConsensusClusterPlus)

title="./RA_ConsensusCluster"
RA_results = ConsensusClusterPlus(R_exp1,maxK=10,reps=1000,pItem=0.8,pFeature=1,
                                  title=title,clusterAlg="hc",distance="euclidean",seed=2333,plot="pdf")

R_exp2<-as.data.frame(t(R_exp[,R_group!='Control']))
R_exp2$class<-RA_results[[2]]$consensusClass

write.csv(R_exp2,"RA_class.csv")
