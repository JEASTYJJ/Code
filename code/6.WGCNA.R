
library(readxl)
j_gene<-read.table("Pyroptosis_genes.txt",header = T)[,1]

z_gene<-read_xls("Autophagy_genes.xls",col_names = NA)
z_gene<-z_gene[-1,]
z_gene<-z_gene$...1

gene_sets1<-list(Pyroptosis=j_gene,Autophagy=z_gene)

library(GSVA)


min.sz = 2 
max.sz = 10000 
parallel.sz= 10 
mx.diff= T 
tau=1 
method='gsva' 
kcdf = "Gaussian" 

exp2<-as.matrix(exp1)
gsva_scores2 = gsva(expr=exp2, gset.idx.list=gene_sets1, 
                    method=method, kcdf=kcdf, min.sz=min.sz, max.sz=max.sz, 
                    parallel.sz=parallel.sz, mx.diff=mx.diff )

gsva_scores2 = as.data.frame(gsva_scores2)
write.csv(gsva_scores2,"Pyroptosis_Autophagy_GSVA.csv")


library(WGCNA)
library(data.table)
library(stringr)
library(gplots)
library(Biobase)
allowWGCNAThreads()

options(stringsAsFactors = FALSE)

exp3<-exp1


m.vars=apply(exp3,1,var)
m.vars.1 = data.frame(m.vars);m.vars.1$names = rownames(exp3)
m.vars.2  = m.vars.1[order(-m.vars.1$m.vars),]
exp3.upper = exp3[rownames(exp3) %in% m.vars.2$names[1:5000] == T,]
multiExpr = as.data.frame(t(exp3.upper))

nGenes = length(multiExpr[1,])
nSamples = length(multiExpr[,1])

powers = c(c(1:10), seq(from = 12, to=20, by=2))
sft = pickSoftThreshold(multiExpr, powerVector = powers, verbose = 5)


pdf("SoftThreshold.pdf",height = 5.5,width = 8)

par(mfrow = c(1,2));
cex1 = 0.85;

plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
abline(h=0.85,col="red")
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
dev.off()


softPower <- sft$powerEstimate
adjacency = adjacency(multiExpr, power = softPower)

TOM = TOMsimilarity(adjacency)

dissTOM = 1-TOM
hierTOM = hclust(as.dist(dissTOM),method="average")

pdf("β_value.pdf",height = 6,width = 10)

ADJ1_cor <- abs(WGCNA::cor( multiExpr,use = "p" ))^softPower
k <- softConnectivity(datE=multiExpr,power=softPower) 

par(mfrow=c(1,2))
hist(k)
scaleFreePlot(k,main="Check Scale free topology\n")
dev.off()

net = blockwiseModules(
  multiExpr,
  power = sft$powerEstimate,
  maxBlockSize = 6000,
  TOMType = "unsigned", minModuleSize = 30,
  reassignThreshold = 0, mergeCutHeight = 0.25,
  numericLabels = TRUE, pamRespectsDendro = FALSE,
  verbose = 3
)
table(net$colors)

pdf("merged_dynamic.pdf",height = 6,width = 10)
mergedColors = labels2colors(net$colors)
table(mergedColors)
plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
dev.off()

pdf("Clustering_module.pdf",height = 6,width = 10)
datExpr_tree<-hclust(dist(multiExpr), method = "average")
par(mar = c(0,5,2,0))
plot(datExpr_tree, main = "Sample clustering", sub="", xlab="", cex.lab = 2, 
     cex.axis = 1, cex.main = 1,cex.lab=1)
dev.off()


design=as.data.frame(t(gsva_scores2))

moduleColors <- labels2colors(net$colors)
MEs0 = moduleEigengenes(multiExpr, moduleColors)$eigengenes
MEs = orderMEs(MEs0)
moduleTraitCor = cor(MEs, design , use = "p");
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples)

textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
                   signif(moduleTraitPvalue, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor)


pdf("Module_trait.pdf",width = 8,height = 8)
par(mar = c(5, 8, 3, 3));
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = colnames(design),####
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = blueWhiteRed(30),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.7,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))
dev.off()

modNames = substring(names(MEs), 3)
geneModuleMembership = as.data.frame(cor(multiExpr, MEs, use = "p"))

MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples))
names(geneModuleMembership) = paste("MM", modNames, sep="")
names(MMPvalue) = paste("p.MM", modNames, sep="")

write.csv(geneModuleMembership,"geneModuleMembership.csv")

risk = as.data.frame(design)
geneTraitSignificance = as.data.frame(cor(multiExpr, risk, use = "p"))
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples))
names(geneTraitSignificance) = paste("GS.", names(risk), sep="")
names(GSPvalue) = paste("p.GS.", names(risk), sep="")

write.csv(geneTraitSignificance,"geneTraitSignificance.csv")

pdf("Pyroptosis_related_genes.pdf",height = 6,width = 6)

module = "turquoise"
column = match(module, modNames);
moduleGenes = moduleColors==module;
par(mfrow = c(1,1));
verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                   abs(geneTraitSignificance[moduleGenes, 1]),
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = "Pyroptosis",
                   main = paste("Module membership vs. gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)
dev.off()

module = "turquoise"
probes = colnames(multiExpr) 
inModule = (moduleColors==module)
modProbes = probes[inModule]
write.csv(modProbes,"turquoise_genes.csv")

