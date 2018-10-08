#PCA方法
turquoise_FPKM=as.data.frame(datExpr[moduleColors=="turquoise"])
p=prcomp(turquoise_FPKM)
PCA=as.data.frame(p$x)
PC1_turquoise=data.frame(rownames(PCA),PCA$PC1)
write.csv(PC1_turquoise,file="turquoise_PC1.csv")
colnames(PC1_turquoise)=c("samples","PC1")
pdf(file="turquoise_PC1.pdf",width = 10,height = 10)
ggplot(data=PC1_turquoise,mapping = aes(x=samples,y=PC1))+geom_bar(stat = 'identity',width = 0.2)+theme(axis.text.x  = element_text(angle=30, vjust=0.5))
heatmap.2(as.matrix(t(turquoise_FPKM)),col = redblue(100),scale = "row",key = TRUE,symkey = FALSE,density.info = "none",trace = "none",cexRow = 0.5,Colv = FALSE,srtCol=45,adjCol = c(1,1))
dev.off()

#ME方法
moduleColors = mergedColors
colorh1=moduleColors
datME=moduleEigengenes(datExpr,colorh1)$eigengenes

dissimME=(1-t(cor(datME, method="p")))/2
hclustdatME=hclust(as.dist(dissimME), method="average" )
par(mfrow=c(1,1))
plot(hclustdatME, main="Clustering tree based of the module eigengenes")
ME=datME[, paste("ME",which.module, sep="")]
par(mfrow=c(2,1), mar=c(0.3, 5.5, 3, 2))
plotMat(t(scale(datExpr[,colorh1==which.module ]) ),
       nrgcols=30,rlabels=F,rcols=which.module,
      main=which.module, cex.main=2)
par(mar=c(5, 4.2, 0, 0.7))
barplot(ME, col=which.module, main="", cex.main=2,
        ylab="eigengene expression",xlab="array sample")
plotMat(t(scale(datExpr[,moduleColors==which.module ]) ),nrgcols=30,rlabels=T,
clabels=T,rcols=which.module,
title=which.module )

ME_i=NULL
colors=c("brown","greenyellow","black","red","magenta","salmon","purple","cyan","tan","pink","lightyellow","lightgreen","green","midnightblue","lightcyan","grey60","yellow","blue","turquoise")
for(i in colors) 
ME_i=cbind(ME_i,datME[, paste("ME",i, sep="")])

colnames(ME_i)=colors
samples=c("gills_1","gills_3","gills_5","gills_7","gills_9","gills_10","gills_11","muscle_0","muscle_1","muscle_3","muscle_5","muscle_7","muscle_9","muscle_10","muscle_11")
rownames(ME_i)=samples
pdf(file="ME_heatmap.pdf")
heatmap.2(as.matrix(t(ME_i)),col = redblue(100),scale = "row",key = TRUE,symkey = FALSE,density.info = "none",trace = "none",cexRow = 0.5,Colv = FALSE,Rowv = FALSE,srtCol=45,adjCol = c(1,1))
library(pheatmap)
library(grid)
draw_colnames_45 <- function (coln, gaps, ...) {
+     coord = pheatmap:::find_coordinates(length(coln), gaps)
+     x = coord$coord - 0.5 * coord$size
+     res = textGrob(coln, x = x, y = unit(1, "npc") - unit(3,"bigpts"), vjust = 0.5, hjust = 1, rot = 45, gp = gpar(...))
+     return(res)}

## 'Overwrite' default draw_colnames with your own version 
assignInNamespace(x="draw_colnames", value="draw_colnames_45",
+                   ns=asNamespace("pheatmap"))
pheatmap(log_ubi,cluster_cols = FALSE)
pheatmap(t(ME_i),cluster_cols = FALSE,cluster_rows = FALSE)
dev.off()
for(i in colors)
     write.csv(as.data.frame(t(datExpr[moduleColors==i])),file=paste("expr",i,".csv",sep=""))
tablecolors=as.data.frame(table(moduleColors))
 i=as.vector(tablecolors$moduleColors)
for (i in color) 
	eval(parse(text=paste(paste("genes",i,sep="_"),'=rownames(as.data.frame(t(datExpr[moduleColors==',i,'])))')))
for (i in 1:nrow(dataB)) eval(parse(text=paste(paste('a',i,sep=''), '= dataB[',i,',]$value')))
