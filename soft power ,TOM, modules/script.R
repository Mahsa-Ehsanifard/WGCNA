clinTraitColor <- numbers2colors(Traits, signed = FALSE)
dendro <- plotDendroAndColors(mrPlusST2, clinTraitColor, 
                              groupLabels = colnames(Traits), 
                              main =  "Sample dendrogram and trait heatmap")

# Choose a set of soft thresholding powers
power <- c(c(1:10), seq(12,20, by = 2))
soft <- pickSoftThreshold(mrPlusex, powerVector = power, verbose = 5)
par(mfrow = c(1,2))
plot(soft$fitIndices[,1], -sign(soft$fitIndices[,3])*soft$fitIndices[,2], 
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"))
text(soft$fitIndices[,1], -sign(soft$fitIndices[,3])*soft$fitIndices[,2],
     labels=power,cex=1,col="red")
abline(h = 0.8, col = "red")
plot(soft$fitIndices[,1], soft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(soft$fitIndices[,1], soft$fitIndices[,5], labels=power, cex=1,col="red")

#####Co-expression similarity and adjacency#####
################################################
adj <- adjacency(datExpr, power = 3) #datExpr = data matrix with removed outlier ample
###Topological Overlap Matrix = TOM ####
##########################################
# Turn adjacency into topological overlap
#creating a matrix for showing the neighbors similarity correlation
TOM <- TOMsimilarity(adj)
dissTOM <- 1-TOM # dissimilarity = distance of each two genes
# now I have co-expression matrix.
# I have network information in adj matrix.
# I have similarity ratio in TOM matrix.
#############

geneTree <- flashClust(as.dist(dissTOM), method = "average")
plot(geneTree,  xlab="", sub="", main = "Gene clustering on TOM-based dissimilarity",
     labels = FALSE, hang = 0.03)

dynamicMods <- cutreeDynamic(dendro = geneTree, distM = dissTOM, deepSplit = 2, 
                             pamRespectsDendro = F, minClusterSize = 30)
table(dynamicMods)
length(dynamicMods)

dynamicColor <- labels2colors(labels = dynamicMods)
table(dynamicColor)

plotDendroAndColors(geneTree, dynamicColor, dendroLabels = F, "Dynamic Tree Cut", 
                    hang = 0.03, addGuide = T, guideHang = 0.05, 
                    main = "Gene dendrogram and module colors")

###Merging of modules###
MElist <- moduleEigengenes(mrPlusex, colors = dynamicColor)
MEs <- MElist$eigengenes
MEdiss <- 1- cor(MEs)

MEtree <- flashClust(as.dist(MEdiss), method = "average")
sizeGrWindow(7, 6)
plot(MEtree, main = "Clustering of module eigengenes", xlab = "", sub = "")

abline(h = 0.25, col = "red")
merge <- mergeCloseModules(mrPlusex, dynamicColor, cutHeight = 0.25, verbose = 3)
mergeColors <- merge$colors
table(mergeColors)
mergedME <- merge$newMEs 

mergedMEdiss <- 1- cor(mergedME)
mergedMETree = flashClust(as.dist(mergedMEdiss), method = "average")
plot(mergedMETree, main = "Clustering of module eigengenes",
     xlab = "", sub = "")

sizeGrWindow(12, 9)
plotDendroAndColors(geneTree, cbind(dynamicColor, mergeColors),
                    c("Dynamic Tree Cut", "Merged dynamic"),
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)

moduleColors <- mergeColors
colorOrder <- c("grey", standardColors(50))
moduleLabel <- match(moduleColors, colorOrder)-1
MEs <- mergedME

###Quantifying module-trait associations##
MEs0 = moduleEigengenes(mrPlusex, moduleColors)$eigengenes
MEs = orderMEs(MEs0)
modulTraitcor <- cor(MEs, Traits, use = "p") #correlation between each ME and each trait
moduleTraitPval <- corPvalueStudent(modulTraitcor, nsamples) #correlations with p.values
adjPval <- p.adjust(moduleTraitPval)

sizeGrWindow(10,6)
textMatrix <- paste(signif(modulTraitcor, 2), "\n(",
                    signif(moduleTraitPval, 1), ")", sep = "")
dim(textMatrix) <- dim(modulTraitcor)
par(mar = c(6, 8.5, 3,3))
sizeGrWindow(12,9)
par(cex =(0.8))
par(mfrow = c(2,2))
labeledHeatmap(Matrix = modulTraitcor,naColor = "white",cex.lab.x = 0.5,cex.lab.y = 0.5,
               xLabels = colnames(Traits),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))
