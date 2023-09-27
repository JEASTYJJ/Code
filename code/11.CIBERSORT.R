library(preprocessCore)

X <- read.table("LM22.txt",header=T,row.names = 1,sep = "\t")
Y <- exp1

X <- data.matrix(X)
Y <- data.matrix(Y)
Y[1:4,1:4]
X[1:4,1:4]
dim(X)
dim(Y)

X <- X[order(rownames(X)),]
Y <- Y[order(rownames(Y)),]


if(max(Y) < 50) {Y <- 2^Y} 

QN = T 
if(QN == TRUE){
  tmpc <- colnames(Y)
  tmpr <- rownames(Y)
  Y <- normalize.quantiles(Y)
  colnames(Y) <- tmpc
  rownames(Y) <- tmpr
}


Xgns <- row.names(X)
Ygns <- row.names(Y)
YintX <- Ygns %in% Xgns 
Y <- Y[YintX,] 
XintY <- Xgns %in% row.names(Y)
X <- X[XintY,]
dim(X)
dim(Y)

X <- (X - mean(X)) / sd(as.vector(X)) 
Y[1:4,1:4]
X[1:4,1:4]
boxplot(X[,1:4])
save(X,Y,file = 'input.Rdata')

CoreAlg <- function(X, y){
  
  svn_itor <- 3
  
  res <- function(i){
    if(i==1){nus <- 0.25}
    if(i==2){nus <- 0.5}
    if(i==3){nus <- 0.75}
    model<-e1071::svm(X,y,type="nu-regression",kernel="linear",nu=nus,scale=F)
    model
  }
  
  if(Sys.info()['sysname'] == 'Windows') out <- parallel::mclapply(1:svn_itor, res, mc.cores=1) else
    out <- parallel::mclapply(1:svn_itor, res, mc.cores=svn_itor)
  
  nusvm <- rep(0,svn_itor)
  corrv <- rep(0,svn_itor)
  
  t <- 1
  while(t <= svn_itor) {
    weights = t(out[[t]]$coefs) %*% out[[t]]$SV
    weights[which(weights<0)]<-0
    w<-weights/sum(weights)
    u <- sweep(X,MARGIN=2,w,'*')
    k <- apply(u, 1, sum)
    nusvm[t] <- sqrt((mean((k - y)^2)))
    corrv[t] <- cor(k, y)
    t <- t + 1
  }
  
  rmses <- nusvm
  mn <- which.min(rmses)
  model <- out[[mn]]
  
  q <- t(model$coefs) %*% model$SV
  q[which(q<0)]<-0
  w <- (q/sum(q))
  
  mix_rmse <- rmses[mn]
  mix_r <- corrv[mn]
  
  newList <- list("w" = w, "mix_rmse" = mix_rmse, "mix_r" = mix_r)
  
}

load(file = 'input.Rdata')
Y[1:4,1:4]
X[1:4,1:4]
dim(X)
dim(Y)
library(preprocessCore)
library(parallel)
library(e1071)

itor <- 1
Ylist <- as.list(data.matrix(Y))
dist <- matrix()
perm=1000
while(itor <= perm){
  print(itor) 
  
  yr <- as.numeric(Ylist[ sample(length(Ylist),dim(X)[1]) ])
  
  yr <- (yr - mean(yr)) / sd(yr)
  
  result <- CoreAlg(X, yr)
  
  mix_r <- result$mix_r
  
  if(itor == 1) {dist <- mix_r}
  else {dist <- rbind(dist, mix_r)}
  
  itor <- itor + 1
}
####
newList <- list("dist" = dist)
nulldist=sort(newList$dist)

save(nulldist,file = 'nulldist_perm_1000.Rdata')


header <- c('Mixture',colnames(X),"P-value","Correlation","RMSE")
print(header)

load(file = 'nulldist_perm_1000.Rdata')
print(nulldist)
fivenum(print(nulldist)) 

output <- matrix()
itor <- 1
mix <- dim(Y)[2] 
pval <- 1
while(itor <= mix){
  
  y <- Y[,itor]
  
  y <- (y - mean(y)) / sd(y)
  
  result <- CoreAlg(X, y)
  
  w <- result$w
  mix_r <- result$mix_r
  mix_rmse <- result$mix_rmse
  
  if(pval > 0) {pval <- 1 - (which.min(abs(nulldist - mix_r)) / length(nulldist))}
  
  out <- c(colnames(Y)[itor],w,pval,mix_r,mix_rmse)
  if(itor == 1) {output <- out}
  else {output <- rbind(output, out)}
  itor <- itor + 1
  
}
head(output)

write.table(rbind(header,output), file="CIBERSORT-Results.txt", sep="\t", row.names=F, col.names=F, quote=F)

obj <- rbind(header,output)
obj <- obj[,-1]
obj <- obj[-1,]
obj <- matrix(as.numeric(unlist(obj)),nrow=nrow(obj))
rownames(obj) <- colnames(Y)
colnames(obj) <- c(colnames(X),"P-value","Correlation","RMSE")
obj[1:4,1:4]
save(obj,file = 'output_obj.Rdata')
write.csv(obj,"CIBERSORT_Results.csv")
