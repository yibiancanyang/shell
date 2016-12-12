library(WGCNA)
options(stringsAsFactors = FALSE);
enableWGCNAThreads();
setwd("/home/thomas/GCH");
lnames=load(file="oyster_dried_WGCNA.RData")
TOM = TOMsimilarityFromExpr(datExpr, power = 4);
module = "cyan";
# Select module probes
probes = names(datExpr);
inModule = (moduleColors==module);
modProbes = probes[inModule];
# Select the corresponding Topological Overlap
modTOM = TOM[inModule, inModule];
dimnames(modTOM) = list(modProbes, modProbes);
# Export the network into an edge list file VisANT can read
vis = exportNetworkToVisANT(modTOM,
file = paste("VisANTInput-", module, ".txt", sep=""),
weighted = TRUE,
threshold = 0,
);
save(TOM,file="oyster_dried_TOM.RData")