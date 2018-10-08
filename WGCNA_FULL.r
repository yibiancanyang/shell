setwd("~/oyster_dry_final/")
counts=read.csv("counts.csv")
counts=read.csv("counts.csv",sep="\t")
rownames(counts)=counts[,1]
counts=counts[,-1]
library(DESeq2)
counts=as.matrix(counts[rowSums(counts) > 5,])
normal_counts=varianceStabilizingTransformation(counts, blind = FALSE)
library(sva)
pheno=read.csv("pheno.csv")
batch = pheno$tissue
modcombat = model.matrix(~1, data=pheno)
combat_edata = ComBat(dat=normal_counts, batch=batch, mod=modcombat, par.prior=TRUE, prior.plots=FALSE)
library(WGCNA)
options(stringsAsFactors = FALSE)
enableWGCNAThreads()
datExpr0=as.data.frame(t(combat_edata))
gsg = goodSamplesGenes(datExpr0, verbose = 3);
gsg$allOK
sampleTree = hclust(dist(datExpr0), method = "average")
pdf(file = "Sample clustering.pdf", wi = 12, he = 9)
par(cex = 0.6)
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5,
      cex.axis = 1.5, cex.main = 2)
dev.off()
datExpr = datExpr0
nGenes = ncol(datExpr)
nSamples = nrow(datExpr)
save(datExpr, file = "oyster-combat-01-dataInput.RData")
powers = c(c(1:10), seq(from = 12, to=20, by=2))
sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)
pdf(file = "Scale independence.pdf", wi = 9, he = 5)
par(mfrow = c(1,2));
cex1 = 0.9;
# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
# this line corresponds to using an R^2 cut-off of h
abline(h=0.90,col="red")
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
dev.off()
#自动计算
net = blockwiseModules(datExpr, power = 3, maxBlockSize = nGenes,
                       TOMType = "unsigned", minModuleSize = 30,
                       reassignThreshold = 0, mergeCutHeight = 0.25,
                       numericLabels = TRUE, pamRespectsDendro = FALSE,
                       saveTOMs = TRUE,
                       saveTOMFileBase = "oyster-FPKM-TOM",
                       verbose = 3)
					   
#手动计算
adjacency = adjacency(datExpr, power = 3);
TOM = TOMsimilarity(adjacency);
dissTOM = 1-TOM
geneTree = hclust(as.dist(dissTOM), method = "average");
sizeGrWindow(12,9)
plot(geneTree, xlab="", sub="", main = "Gene clustering on TOM-based dissimilarity",labels = FALSE, hang = 0.04);

# We like large modules, so we set the minimum module size relatively high:
minModuleSize = 30;
# Module identification using dynamic tree cut:
dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM,deepSplit = 2, pamRespectsDendro = FALSE,minClusterSize = minModuleSize);
table(dynamicMods)
# Convert numeric lables into colors
dynamicColors = labels2colors(dynamicMods)
table(dynamicColors)
# Plot the dendrogram and colors underneath
sizeGrWindow(8,6)
plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut",dendroLabels = FALSE, hang = 0.03,addGuide = TRUE, guideHang = 0.05,main = "Gene dendrogram and module colors")
# Calculate eigengenes
MEList = moduleEigengenes(datExpr, colors = dynamicColors)
MEs = MEList$eigengenes
# Calculate dissimilarity of module eigengenes
MEDiss = 1-cor(MEs);
# Cluster module eigengenes
METree = hclust(as.dist(MEDiss), method = "average");
# Plot the result
sizeGrWindow(7, 6)
plot(METree, main = "Clustering of module eigengenes",xlab = "", sub = "")
MEDissThres = 0.25
# Plot the cut line into the dendrogram
abline(h=MEDissThres, col = "red")
# Call an automatic merging function
merge = mergeCloseModules(datExpr, dynamicColors, cutHeight = MEDissThres, verbose = 3)
# The merged module colors
mergedColors = merge$colors;
# Eigengenes of the new merged modules:
mergedMEs = merge$newMEs;
sizeGrWindow(12, 9)
pdf(file = "geneDendro.pdf", wi = 9, he = 6)
plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors),
c("Dynamic Tree Cut", "Merged dynamic"),
dendroLabels = FALSE, hang = 0.03,
addGuide = TRUE, guideHang = 0.05)
dev.off()
# Rename to moduleColors
mergedColors = labels2colors(net$colors)
#mergedColors = labels2colors(merge$colors)
moduleColors = mergedColors
# Construct numerical labels corresponding to the colors
colorOrder = c("grey", standardColors(50));
moduleLabels = match(moduleColors, colorOrder)-1;
MEs = mergedMEs;
# Save module colors and labels for use in subsequent parts
save(MEs, moduleLabels, moduleColors, geneTree, file = "FemaleLiver-02-networkConstruction-stepByStep.RData")
color=as.vector(as.data.frame(table(moduleColors))$moduleColors)
#检测模块与处理的相关性
time = data.frame(row.names = rownames(datExpr),time=c(0,1,3,5,7,9,10,11,0,1,3,5,7,9,10,11));

MET = orderMEs(cbind(datME, time))
pdf(file = "plotEigengeneNetworks.pdf", wi = 6, he = 10)
plotEigengeneNetworks(MET, "", marDendro = c(0,4,1,2), marHeatmap = c(3,4,1,2), cex.lab = 0.8, xLabelsAngle
                      = 90)
					  dev.off()
#各模块GO注释
color=as.vector(tablecolors$moduleColors)
for (i in color) 
     eval(parse(text=paste(paste("genes",i,sep="_"),'=rownames(as.data.frame(t(datExpr[moduleColors==i])))')))
library(org.Cgigas.eg.db)
library(clusterProfiler)
for(i in color) 
	eval(parse(text=paste(paste(i,"GO",sep="_"),'=enrichGO(gene=',paste("genes",i,sep="_"),',keytype = "GID", OrgDb=org.Cgigas.eg.db,ont="ALL",pvalueCutoff=0.05,qvalueCutoff = 0.05,readable=TRUE)')))
for(i in color)
	write.csv(as.data.frame(eval(parse(text=paste(i,"GO",sep="_")))),file=paste(i,"_GO.csv",sep=""))
#ME heatmap
moduleColors = mergedColors
datME=moduleEigengenes(datExpr,moduleColors)$eigengenes
library(pheatmap)
library(grid)
draw_colnames_45 <- function (coln, gaps, ...) {
     coord = pheatmap:::find_coordinates(length(coln), gaps)
     x = coord$coord - 0.5 * coord$size
     res = textGrob(coln, x = x, y = unit(1, "npc") - unit(3,"bigpts"), vjust = 0.5, hjust = 1, rot = 45, gp = gpar(...))
     return(res)}

## 'Overwrite' default draw_colnames with your own version 
assignInNamespace(x="draw_colnames", value="draw_colnames_45",
	ns=asNamespace("pheatmap"))
pheatmap(t(datME),cluster_cols = FALSE,cluster_rows = FALSE)
ME_i=NULL
colors=c("magenta","purple","pink","greenyellow","green","red","yellow","turquoise")
for(i in colors) 
ME_i=cbind(ME_i,datME[, paste("ME",i, sep="")])

colnames(ME_i)=colors
rownames(ME_i)=rownames(datExpr)
pheatmap(t(ME_i),cluster_cols = FALSE,cluster_rows = FALSE)
## WGCNA hub genes network
TOM = TOMsimilarityFromExpr(datExpr, power = 4);
module = c("magenta","purple","pink","greenyellow","green","red","yellow");
# Select module probes
probes = names(datExpr);
inModule = is.finite(match(moduleColors, module));
modProbes = probes[inModule];
# Select the corresponding Topological Overlap
modTOM = TOM[inModule, inModule];
dimnames(modTOM) = list(modProbes, modProbes);
nTop = 50;
IMConn = softConnectivity(datExpr[, modProbes]);
top = (rank(-IMConn) <= nTop)
modColors=moduleColors[inModule]
vis = exportNetworkToVisANT(modTOM[top, top],file = paste("VisANTInput-", module, "-top30.txt", sep=""),weighted = TRUE,threshold = 0) 
cyt = exportNetworkToCytoscape(modTOM,
edgeFile = paste("CytoscapeInput-edges-", paste(module, collapse="-"), ".txt", sep=""),
nodeFile = paste("CytoscapeInput-nodes-", paste(module, collapse="-"), ".txt", sep=""),
weighted = TRUE,
threshold = 0,
nodeNames = modProbes[top],
nodeAttr = modColors[top]);
colors = c("darkred","darkturquoise","green","turquoise","yellow","blue")
for(i in colors)
     write.csv(as.data.frame(t(datExpr[moduleColors==i])),file=paste("expr_",i,".csv",sep=""))
#各模块表达量柱形图和热图
pdf("green_ME.pdf",width = 8,height = 7)
which.module="green"
ME=datME[, paste("ME",which.module, sep="")]
par(mfrow=c(2,1), mar=c(0.3, 5.5, 3, 2))
plotMat(t(scale(datExpr[,colorh1==which.module ]) ),
nrgcols=30,rlabels=F,rcols=which.module,
main=which.module, cex.main=2)
par(mar=c(5, 4.2, 0, 0.7))
barplot(ME, col=which.module, main="", cex.main=2,
ylab="eigengene expression",xlab="array sample")
dev.off()

dissTOM = 1-TOM

nSelect = 1000
# For reproducibility, we set the random seed
set.seed(10);
select = sample(nGenes, size = nSelect);
selectTOM = dissTOM[select, select];
# There’s no simple way of restricting a clustering tree to a subset of genes, so we must re-cluster.
selectTree = hclust(as.dist(selectTOM), method = "average")
selectColors = moduleColors[select];
# Open a graphical window
pdf("Network heatmap plot.pdf",width = 9,height = 9)

# Taking the dissimilarity to a power, say 10, makes the plot more informative by effectively changing
# the color palette; setting the diagonal to NA also improves the clarity of the plot
plotDiss = selectTOM^7;
diag(plotDiss) = NA;
TOMplot(plotDiss, selectTree, selectColors, main = "Network heatmap plot, selected genes")
dev.off()