
set.seed(123)
O_exp<-exp[cygene,group_list!="RA"]
O_group<-group_list[group_list!="RA"]

library(plyr)
library(rms)

OA_l<-cbind(as.data.frame(t(O_exp)),O_group)
OA_l$O_group<-ifelse(OA_l$O_group=="Control",0,1)
OA_res<-data.frame(OR=1,CI="1",P=1)
for(i in colnames(OA_l)[1:ncol(OA_l)-1]){
  f<-paste0("O_group==1~",i)
  f<-as.formula(f)
  glm1<-glm(f,family = binomial,data = OA_l)
  glm2<- summary(glm1)
  OR<-round(exp(coef(glm1)),2)
  SE<-glm2$coefficients[,2]
  CI5<-round(exp(coef(glm1)-1.96*SE),2)
  CI95<-round(exp(coef(glm1)+1.96*SE),2)
  CI<-paste0(CI5,'-',CI95)
  P<-round(glm2$coefficients[,4],4)
  res1<-data.frame(OR,CI,P)[-1,]
  OA_res<-rbind(OA_res,res1)
}
OA_res<-OA_res[-1,]

write.csv(OA_res,"OA_logistic.csv")

OA_res1<-OA_res[which(OA_res$P<0.05),]
OA_l_gene<-rownames(OA_res1)

library(glmnet)
O_X<-as.matrix(OA_l[,OA_l_gene])
O_Y<-as.matrix(OA_l$O_group)

f1 = glmnet(O_X, O_Y, family="binomial", nlambda=100, alpha=1)
print(f1)
pdf("OA_Lasso_lambda.pdf",width = 6,height = 6)
plot(f1, xvar="lambda")
dev.off()

cvfit=cv.glmnet(O_X,O_Y)
pdf("OA_Mean_Squared_Error.pdf",width = 6,height = 6)
plot(cvfit)
dev.off()

cvfit$lambda.min
l.coef1<-coef(cvfit$glmnet.fit,s=cvfit$lambda.min,exact = F)
l.coef1
OA_lasso<-l.coef1[which(l.coef1!=0),]
OA_lasso_gene<-names(OA_lasso)[-1]
write.csv(OA_lasso_gene,"OA_lasso_gene.csv")

OA_df<-matrix(0,ncol = 2,nrow = ncol(O_exp))
for(i in length(OA_lasso_gene)){
  OA_df[,2]<-OA_df[,2]+as.matrix(O_exp[OA_lasso_gene[i],])*OA_lasso[i-1]
}

OA_df[,1]<-OA_l$O_group

library(pROC)
library(ggplot2)

OA_rocobj <- roc(OA_df[,1], OA_df[,2],
                 smooth = F       
) 
cutOffPoint <- coords(OA_rocobj, "best")
cutOffPointText <- paste0(round(cutOffPoint[1],3),"(",round(cutOffPoint[2],3),",",round(cutOffPoint[3],3),")")

auc<-auc(OA_rocobj)[1]
auc_text<-paste0("AUC = ",round(auc,4))
pdf("OA_lasso_ROC.pdf",width = 6,height = 6)
ggroc(OA_rocobj,
      color="red",
      size=1,
      legacy.axes = F 
)+
  theme_classic()+
  geom_segment(aes(x = 1, y = 0, xend = 0, yend = 1),        
               colour='grey', 
               linetype = 'dotdash'
  ) +
  geom_label(aes(x = 0.25,y = 0.25,label=auc_text),size =6)+
  theme(axis.title.x = element_text(size = 15),
        axis.title.y = element_text(size = 15))


dev.off()
R_exp<-exp[cygene,group_list!="OA"]
R_group<-group_list[group_list!="OA"]

RA_l<-cbind(as.data.frame(t(R_exp)),R_group)
RA_l$R_group<-ifelse(RA_l$R_group=="Control",0,1)
RA_res<-data.frame(OR=1,CI="1",P=1)
for(i in colnames(RA_l)[1:ncol(RA_l)-1]){
  f<-paste0("R_group==1~",i)
  f<-as.formula(f)
  glm1<-glm(f,family = binomial,data = RA_l)
  glm2<- summary(glm1)
  OR<-round(exp(coef(glm1)),2)
  SE<-glm2$coefficients[,2]
  CI5<-round(exp(coef(glm1)-1.96*SE),2)
  CI95<-round(exp(coef(glm1)+1.96*SE),2)
  CI<-paste0(CI5,'-',CI95)
  P<-round(glm2$coefficients[,4],4)
  res1<-data.frame(OR,CI,P)[-1,]
  RA_res<-rbind(RA_res,res1)
}
RA_res<-RA_res[-1,]

write.csv(RA_res,"RA_logistic.csv")

RA_res1<-RA_res[which(RA_res$P<0.05),]
RA_l_gene<-rownames(RA_res1)

library(glmnet)
R_X<-as.matrix(RA_l[,RA_l_gene])
R_Y<-as.matrix(RA_l$R_group)

f2 = glmnet(R_X, R_Y, family="binomial", alpha=1)
print(f2)
pdf("RA_Lasso_lambda.pdf",width = 8,height = 8)
plot(f2, xvar="lambda")
dev.off()

cvfit1=cv.glmnet(R_X,R_Y)
pdf("RA_Mean_Squared_Error.pdf",width = 6,height = 6)
plot(cvfit1)
dev.off()

cvfit1$lambda.min
l.coef2<-coef(cvfit1$glmnet.fit,s=cvfit1$lambda.min,exact = F)
l.coef2
RA_lasso<-l.coef2[which(l.coef2!=0),]
RA_lasso_gene<-names(RA_lasso)[-1]
write.csv(RA_lasso_gene,"RA_lasso_gene.csv")

RA_df<-matrix(0,ncol = 2,nrow = ncol(R_exp))
for(i in length(RA_lasso_gene)){
  RA_df[,2]<-RA_df[,2]+as.matrix(R_exp[RA_lasso_gene[i],])*RA_lasso[i-1]
}

RA_df[,1]<-RA_l$R_group

library(pROC)
library(ggplot2)

RA_rocobj <- roc(RA_df[,1], RA_df[,2],
                 smooth = F       
) 
cutOffPoint <- coords(RA_rocobj, "best")
cutOffPointText <- paste0(round(cutOffPoint[1],3),"(",round(cutOffPoint[2],3),",",round(cutOffPoint[3],3),")")

auc<-auc(RA_rocobj)[1]
auc_text<-paste0("AUC = ",round(auc,4))
pdf("RA_lasso_ROC.pdf",width = 6,height = 6)
ggroc(RA_rocobj,
      color="red",
      size=1,
      legacy.axes = F 
)+
  theme_classic()+
  geom_segment(aes(x = 1, y = 0, xend = 0, yend = 1),        
               colour='grey', 
               linetype = 'dotdash'
  ) +
  geom_label(aes(x = 0.25,y = 0.25,label=auc_text),size =6)+
  theme(axis.title.x = element_text(size = 15),
        axis.title.y = element_text(size = 15))


dev.off()

