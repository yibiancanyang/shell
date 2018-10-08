#!R
#download data and find significant genes
library(GEOquery)
library(limma)
studyID='GSE61019'
destdir='F:/sls'
GSE61019=getGEO("GSE61019")
exprSet=GSE61019[[1]]
pdata=pData(GSE61019[[1]])
write.csv(exprSet,paste0(studyID,'_exprSet.csv'))
write.csv(pdata,paste0(studyID,'_metadata.csv'))
group_list = c(rep('control',4),rep('case',3))
design = model.matrix(~ -1+ group_list )
colnames(design)= c("control", "case")
contrast.matrix = makeContrasts(control-case, levels=design)
fit = lmFit(GSE61019[[1]], design)
fit2 = contrasts.fit(fit, contrast.matrix)
fit2 = eBayes(fit2)
sig=topTable(fit2,number = Inf ,p.value = 0.05,lfc=2)
#draw a heatmap
ematrix=exprs(exprSet)
ema_sig=merge(sig,ematrix,by = "row.names")
rownames(ema_sig)=ema_sig[,1]
ema_sig_heatmap=ema_sig[,-c(1:37)]
annotation_col = data.frame(group = factor(rep(c("control","A-T case"),c(4,3))))
rownames(annotation_col)=colnames(ema_sig_heatmap)
pdf(file = "sig_heatmap.pdf",width = 12,height = 600)
pheatmap(ema_sig_heatmap,cluster_cols = FALSE,annotation_col = annotation_col)
dev.off()
#GO and KEGG enrichment  
library(clusterProfiler)
library(org.Hs.eg.db)
ego=enrichGO(gene = genes,ont = "ALL",OrgDb = org.Hs.eg.db)
ekk=enrichKEGG(gene = genes,organism = 'hsa',pvalueCutoff = 0.05)
write.csv(summary(ekk),"KEGG-enrich.csv",row.names =F)
write.csv(summary(ego),"GO-enrich.csv",row.names =F)
 for (i in color) 
    eval(parse(text=paste(paste("genes",i,sep="_"),'=rownames(as.data.frame(t(datExpr[,moduleColors==i])))')))
anno_black=anno$Entrez.Gene[anno$ID%in%genes_black]
selected=c("magenta","salmon","darkolivegreen","violet","lightyellow","saddlebrown","skyblue","green","cyan","darkgreen","black","red","greenyellow","turquoise","yellow","paleturquoise","sienna3","darkmagenta","darkred","darkturquoise","lightgreen","darkgrey","royalblue","darkorange","blue","orange")
for (i in selected) 
    eval(parse(text=paste(paste("anno",i,sep="_"),'=anno$Entrez.Gene[anno$ID%in%',paste("genes",i,sep="_"),']')))
for(i in selected) 
	eval(parse(text=paste(paste(i,"GO",sep="_"),'=enrichGO(gene=',paste("anno",i,sep="_"),',OrgDb=org.Mm.eg.db,ont="ALL",pvalueCutoff=0.05,qvalueCutoff = 0.05,readable=TRUE)')))
for(i in selected) 
	eval(parse(text=paste(paste(i,"KEGG",sep="_"),'=enrichKEGG(gene=',paste("anno",i,sep="_"),",organism = 'mmu',pvalueCutoff=0.05,qvalueCutoff = 0.05)")))
for(i in selected)
	write.csv(as.data.frame(eval(parse(text=paste(i,"GO",sep="_")))),file=paste(i,"_GO.csv",sep=""))
for(i in selected)
	write.csv(as.data.frame(eval(parse(text=paste(i,"KEGG",sep="_")))),file=paste(i,"_KEGG.csv",sep=""))
	geneList = sort(glist4, decreasing = TRUE)

	
	