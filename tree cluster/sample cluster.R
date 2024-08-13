library(WGCNA)
options(stringsAsFactors = FALSE)

tmrna <- t(geneSigDat)
mrnaExpr <- as.data.frame(tmrna)

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

# Plot a line to show the cut
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