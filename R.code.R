library(WGCNA)
gsg <- goodSamplesGenes(datExpr = mrnaExpr, verbose = 3)
gsg$allOK
#####good genes and good samples are identified already
#####Cluster the transposed matrix to identify sample outliers
#sample clustering
library(flashClust)
sampleTree <- hclust(dist(mrnaExpr), method = "average") #heirarichal clustering

# Plot the sample tree: Open a graphic output window of size 12 by 9 inches
par(cex = 0.6);
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5,
     cex.axis = 1.5, cex.main = 2)

###### detection of outlier
#line cut for removing the outlier sample
abline(h = 100, col = "red")

# Determine cluster under the line
clust <- cutreeStatic(sampleTree, cutHeight = 100, minSize = 10)
table(clust)

#Remove the outlier and construct the main data frame
# clust 1 contains the samples we want to keep.
#two outlier samples identified should be removed
keepsample <- clust==1
mrdatExp <- mrnaExpr[keepsample, ] 

# Convert traits to a color representation: white means low, red means high, grey means missing entry
traitcolor <- numbers2colors(braftrait, signed = F)
# Plot the sample dendrogram and the colors underneath.
dendro <- plotDendroAndColors(sampleTree2, traitcolor, groupLabels = names(braftrait), 
                              main =  "Sample dendrogram and trait heatmap")


# choose power based on SFT criterion
power <- c(c(1:10), seq(12,20, by = 2))
# Call the network topology analysis function
#****** SFT = Scale Free Topology *******
sft <- pickSoftThreshold(mrdatExp, powerVector = power, verbose = 5)

# Plot the results:
sizeGrWindow(12,9)
par(mfrow = c(1,2))
# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2], 
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"))
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=power,cex=1,col="red")
# this line corresponds to using an R^2 cut-off of h
abline(h = 0.80, col = "red")
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=power, cex=1,col="red")


###Co-expression similarity and adjacency###
adj <- adjacency(datExpr, power = 3)

###Topological Overlap Matrix = TOM ####
TOM <- TOMsimilarity(adj)
dissTOM <- 1-TOM # dissimilarity = distance of each two genes

# define a hierarchical tree with hclust or with flashClust
hierTOM=hclust(as.dist(dissTOM),method="average")
colorStaticTOM=cutreeStaticColor(hierTOM,cutHeight=0.995,
                                 minSize=30)

labelDynamicHybrid=cutreeDynamic(hierTOM,distM=dissTOM,
                                 cutHeight=0.995,deepSplit=1,pamRespectsDendro=FALSE,
                                 minClusterSize=30)
colorDynamicHybridTOM=labels2colors(labelDynamicHybrid)


# Call the hierarchical clustering function
geneTree <- flashClust(as.dist(dissTOM), method = "average")
# Plot the resulting clustering tree (dendrogram)
plot(geneTree,  xlab="", sub="", main = "Gene clustering on TOM-based dissimilarity",
     labels = FALSE, hang = 0.03)

# Module identification using dynamic tree cut:
dynamicMods <- cutreeDynamic(dendro = geneTree, distM = dissTOM, deepSplit = 2, 
                             pamRespectsDendro = F, minClusterSize = 2)

# Convert numeric lables into colors
dynamicColor <- labels2colors(labels = dynamicMods)
table(dynamicColor)

# Plot the dendrogram and colors underneath
plotDendroAndColors(geneTree, dynamicColor, dendroLabels = F, "Dynamic Tree Cut", 
                    hang = 0.03, addGuide = T, guideHang = 0.05, 
                    main = "Gene dendrogram and module colors")

##Merging of modules###
MElist <- moduleEigengenes(datExpr, colors = dynamicColor)
MEs <- MElist$eigengenes

MEdiss <- 1- cor(MEs)
# Cluster module eigengenes
#heirerichal clustering for eigengenes to show closed modules if they're exist
MEtree <- flashClust(as.dist(MEdiss), method = "average")
sizeGrWindow(7, 6)
plot(MEtree, main = "Clustering of module eigengenes", xlab = "", sub = "")

abline(h = 0.25, col = "red")
# Call an automatic merging function
#merging close modules
merge <- mergeCloseModules(datExpr, dynamicColor, cutHeight = 0.25, verbose = 3)
mergeColors <- merge$colors
table(mergeColors)
mergedME <- merge$newMEs # Eigengenes of the new merged modules

mergedMEdiss <- 1- cor(mergedME)
mergedMETree = flashClust(as.dist(mergedMEdiss), method = "average")
plot(mergedMETree, main = "Clustering of module eigengenes",
     xlab = "", sub = "")

plotEigengeneNetworks(mergedME, "", marDendro = c(0,4,1,2), marHeatmap = c(7,7,1,0),
                      heatmapColors = blueWhiteRed(50))
sizeGrWindow(12, 9)
plotDendroAndColors(geneTree, cbind(dynamicColor, mergeColors),
                    c("Dynamic Tree Cut", "Merged dynamic"),
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)

# Rename to module Colors
moduleColors <- mergeColors #put the merged colored in this obj
# Construct numerical labels corresponding to the colors
colorOrder <- c("grey", standardColors(50)) #ordering the colors in which grey is the first and 
## other colors are after sequentioally by numbers
moduleLabel <- match(moduleColors, colorOrder)-1


####Quantifying module-trait associations######
### module-trait heatmap###
MEs0 = moduleEigengenes(datExpr, moduleColors)$eigengenes
MEs = orderMEs(MEs0)
modulTraitcor <- cor(MEs, datTraits, use = "p") #correlation between each ME and each trait
moduleTraitPval <- corPvalueStudent(modulTraitcor, nsamples) #correlations with p.values
adjPval <- p.adjust(moduleTraitPval) # for more significance


plotEigengeneNetworks(MEs,"",marDendro=c(0,4,1,2),
                      marHeatmap=c(3,4,1,2),cex.lab=0.8,xLabelsAngle=90)

sizeGrWindow(12,9)
textMatrix <- paste(signif(modulTraitcor, 2), "\n(",
                    signif(moduleTraitPval, 1), ")", sep = "")
dim(textMatrix) <- dim(modulTraitcor)

par(mar = c(6, 8.5, 3, 3))
# OR
par(mar = c(6, 12, 3, 0))
# Display the correlation values within a heatmap plot
sizeGrWindow(12,9)
par(cex =(0.7))
par(mar = c(6,15,2,0))
labeledHeatmap(Matrix = modulTraitcor,
               xLabels = names(datTraits),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.6, yColorWidth = 0.06,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))



##Gene Significance and Module Membership###
modNames <- substring(names(MEs), 3)
geneMODmem <- as.data.frame(cor(datExpr, MEs, use = "p"))

#removing the genes with low membership degree/corr to their own modules
#first I should make p.vaues of each gene corr to their own modules
MMpvalue <- as.data.frame(corPvalueStudent(as.matrix(geneMODmem), nSamples = nsamples))
#make name of module membership columns
names(geneMODmem) <- paste("MM", modNames, sep = "")
names(MMpvalue) <- paste("p.mm", modNames, sep = "")
geneTraitSig <- as.data.frame(cor(datExpr, datTraits, use = "p")) #corr of each gene with trait
GSpvalue <- as.data.frame(corPvalueStudent(as.matrix(geneTraitSig), nsamples))
names(geneTraitSig) = paste("GS.", names(datTraits), sep="")
names(GSpvalue) = paste("p.GS.", names(datTraits), sep="")


#####Intramodular analysis###
### identifying hub genes from hub modules####
brown = "brown" #the best corr of module with trait in heatmap
column <- match(brown, modNames)
moduleGenes <- moduleColors==brown
sizeGrWindow(3,3)
par(mfrow = c(2,4))
verboseScatterplot(abs(geneMODmem[moduleGenes, column]),
                   abs(geneTraitSig[moduleGenes, 1]),
                   xlab = paste("Module Membership in", brown,"module"),
                   ylab = "Gene significance for brafv600eplus",
                   main = paste("Module membership vs. gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = brown, abline = T)

abline(h = 0.20,v = 0.80, col = "red")


Probes = names(datExpr)
inModule <- is.finite(match(moduleColors, brown))
modProbes <- Probes[inModule]
modProbes <- data.frame(modProbes)
library(org.Hs.eg.db)
annot_blue <- select(org.Hs.eg.db, keys = modProbes$modProbes, columns = "SYMBOL", keytype = "ENSEMBL")

GenemmPval <- MMpvalue[MMpvalue$p.mmbrown<0.05,]
names(GenemmPval)
GeneGSpval <- data.frame(GSpvalue[GSpvalue$p.GS.v600e_plus<0.05,])
goodgeneMM <- modProbes[modProbes$modProbes %in% as.character(rownames(GenemmPval)),]
goodGS <- modProbes[modProbes$modProbes %in% as.character(rownames(GeneGSpval)),]
MMgene <- goodgeneMM[goodgeneMM %in% 
                       as.character(rownames(geneMODmem)[geneMODmem$MMbrown>0.80])]

GSgene <- goodGS[goodGS %in% 
                   as.character(rownames(geneTraitSig)[geneTraitSig$GS.v600e_plus>0.20])]

MM_GSGene <- MMgene[MMgene %in% as.character(GSgene)]
MM_GSGene <- data.frame(MM_GSGene)

