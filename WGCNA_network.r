library(WGCNA)
options(stringsAsFactors = FALSE);
enableWGCNAThreads();
setwd("/home/thomas/GCH");
lnames=load(file="oyster_dried_WGCNA.RData")
TOM = TOMsimilarityFromExpr(datExpr, power = 4);
module = c("darkred","darkturquoise","green","turquoise","yellow","blue");
# Select module probes
probes = names(datExpr);
inModule = is.finite(match(moduleColors, module));
modProbes = probes[inModule];
# Select the corresponding Topological Overlap
modTOM = TOM[inModule, inModule];
dimnames(modTOM) = list(modProbes, modProbes);
nTop = 10;
IMConn = softConnectivity(datExpr[, modProbes]);
top = (rank(-IMConn) <= nTop)
vis = exportNetworkToVisANT(modTOM[top, top],
file = paste("VisANTInput-", module, "-top.txt", sep=""),
weighted = TRUE,
threshold = 0,
 )


cyt = exportNetworkToCytoscape(modTOM[top, top], edgeFile = paste("CytoscapeInput-edges-", paste(modules, collapse="-"), ".txt", sep=""), nodeFile = paste("CytoscapeInput-nodes-", paste(modules, collapse="-"), ".txt", sep=""), weighted = TRUE)
# Export the network into an edge list file VisANT can read
cyt = exportNetworkToCytoscape(modTOM,
edgeFile = paste("CytoscapeInput-edges-", paste(module, collapse="-"), ".txt", sep=""),
nodeFile = paste("CytoscapeInput-nodes-", paste(module, collapse="-"), ".txt", sep=""),
weighted = TRUE,
threshold = 0.2,
nodeNames = modProbes,
nodeAttr = moduleColors[inModule]);
save(TOM,file="oyster_dried_TOM.RData")