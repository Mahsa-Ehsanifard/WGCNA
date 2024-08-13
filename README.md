# WGCNA
## WGCNA package in R

![](https://img.shields.io/badge/Version-%201.72_5-%20informational?style=plastic
)
![](https://img.shields.io/badge/Source-Bioconductor-9cf?style=plastic
)
![](https://img.shields.io/badge/Install-Rstudio-purple?style=plastic
)
![](https://img.shields.io/badge/depends-flashClust-yellowgreen?style=plastic
)
![](https://img.shields.io/badge/License-GPL(%3E%3D%202)-aqua?style=plastic
)

## description: Weighted Gene Co-expression Network Analysis

Useful links:

[WGCNA: an R package for weighted correlation network analysis](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-9-559)

[book; Weighted Network Analysis](https://books.google.com/books?id=ZCh06NgMFesC&pg=PR11&dq=wgcna:+an+r+package+for+weighted+correlation+network+analysis&hl=en&newbks=1&newbks_redir=1&sa=X&ved=2ahUKEwju5qiQzvGHAxU9daQEHU7YIAwQ6AF6BAgCEAI)

[WGCNA; bioconductor](https://www.bioconductor.org/packages/release/bioc/vignettes/BioNERO/inst/doc/vignette_01_GCN_inference.html)


The **WGCNA** pipeline is expected an input matrix of normalized expression values including samples in columns and gene names on rows. There is no limitation for the methods exploring the expression values; *RNA-Seq* or *microarray* methods. We can use *GEO* or *TCGA* expression profiles for this analysis. 

#### Installation

```{r}
BiocManager::install("WGCNA")
library(WGCNA)
```

Usually we need to rotate or *transpose* the rows with the columns in the matrix using `t` function.

> rows => samples , columns => genes

The first and important step is **clustering** the samples to identify the *outliers*. To do this we use `flashClust` package.

```{r}
library(flashClust)
```

To perform *heirarichal clustering* we use `hclust` function to construct clusters. The samples are clustered by **distance** between them based on the expression values.

```{r}
-> hclust(dist(matrix), method = "average")
```

### Detecting the outlier samples

The graphical tree output is based on the *height* values on axis y between samples. The *outlier* is a sample locates in a far *distance* of others without any connections to others. This sample has a height more than other samples and is located in top of the tree with a high height. It means that the expression values in this sample are not matched or close to others. The outlier sample should be remove.

* It would be **more than one** outlier samples in a matrix.

* Removing the outliers is based on a cut-line to cut the tree cluster according to the highest height in which the samples are close together. It means that other samples above this line display as outliers because they are in a far distance of others.

```{r}
#line cut for removing the outlier sample
# the highest height is 100 for example
abline(h = 100, col = "red") 
```

```{r}
<- cutreeStatic(sampleTree, cutHeight = 100)
table(clust)
```

* clust 1 contains the samples we want to keep.

```{r}
keepsample <- clust==1
```

* The outlier samples are removed from the main matrix.

### loading traits

We need a table of *trait* information such as *clinical data, molecular characteristics, or phenotypic and genotypic features* to analyze the correlation and relationship between traits and gene expressions. In this step, it is needed to investigate the relationships of traits with gene modules to identify hub genes correlated strongly with an important features of samples.

* Trait could be the features or characteristics of genes such as regulation levels.

* Clinical data includes staging,  mutations, molecular or cellular features, etc...

### Choose a set of soft thresholding power

Choosing a **soft power (β)** is an important step to detect modules.

* power number is a critical index to identify gene module packing

* β parameter will be to calculate our adjacency matrix.

* The `pickSoftThreshold` function calculates multiple networks all based on different β values and returns a data frame with the **R2** values for the networks **scale-free topology** model fit as well as the **mean connectivity** measures.

```{r}
pickSoftThreshold(matrix)
```

```
       Power SFT.R.sq slope truncated.R.sq mean.k. median.k. max.k.
## 1      1   0.0278  0.345          0.456  747.00  762.0000 1210.0
## 2      2   0.1260 -0.597          0.843  254.00  251.0000  574.0
## 3      3   0.3400 -1.030          0.972  111.00  102.0000  324.0
## 4      4   0.5060 -1.420          0.973   56.50   47.2000  202.0
## 5      5   0.6810 -1.720          0.940   32.20   25.1000  134.0
## 6      6   0.9020 -1.500          0.962   19.90   14.5000   94.8
## 7      7   0.9210 -1.670          0.917   13.20    8.6800   84.1
## 8      8   0.9040 -1.720          0.876    9.25    5.3900   76.3
## 9      9   0.8590 -1.700          0.836    6.80    3.5600   70.5
## 10    10   0.8330 -1.660          0.831    5.19    2.3800   65.8
## 11    12   0.8530 -1.480          0.911    3.33    1.1500   58.1
## 12    14   0.8760 -1.380          0.949    2.35    0.5740   51.9
## 13    16   0.9070 -1.300          0.970    1.77    0.3090   46.8
## 14    18   0.9120 -1.240          0.973    1.39    0.1670   42.5
## 15    20   0.9310 -1.210          0.977    1.14    0.0951   38.7
```

Plot the R2 values as a function of the soft thresholds

> We should be *maximizing* the R2 (β) value and *minimizing* mean connectivity.

```{r}
par(mfrow=c(1,2))
plot(sft$fitIndices[,1],
     -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology
Model Fit,
signed Rˆ2",type="n",main=paste("Scale independence"))
text(sft$fitIndices[,1],
     -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,col="red")
abline(h=0.80,col="red")
plot(sft$fitIndices[,1],sft$fitIndices[,5],type="n",
     xlab="Soft Threshold (power)",ylab="Mean Connectivity",
     main=paste("Mean connectivity"))
text(sft$fitIndices[,1],sft$fitIndices[,5],labels=powers,
     col="red")
```

* We can determine the soft power threshold which is a number as it is the β that retains the *highest* mean connectivity (above zero) while reaching an R2 value above **0.80**.

> NOTE: the higher the value, the stronger the connection strength will be of highly correlated gene expression profiles and the more devalued low correlations will be.

Now we have the soft threshold power determined we can call on the `adjacency` function. This function calculates the *similarity* measurement and transforms the similarity by the adjacency function and generates a **weighted network adjacency matrix**.

```{r}
adjacency(matrix, power = 3)
```

### Topological Overlap Matrix = TOM

* Turn adjacency into topological overlap

* creating a matrix for showing the neighbors similarity correlation

```{r}
TOM <- TOMsimilarity(adjacency)
```

To convert this matrix into a dissimilarity matrix we can subtract the TOM object from 1.

```{r}
dissTOM <- 1-TOM
```

#### Hierarchical Clustering Analysis

The dissimilarity/distance measures are then clustered using *linkage hierarchical* clustering and a dendrogram (cluster tree) of genes is constructed.

```{r}
hierTOM = hclust(as.dist(dissTOM),method="average")
```

```{r}
#plotting the dendrogram
sizeGrWindow(12,9)
plot(geneTree, xlab="", sub="", main = "Gene clustering on TOM-based dissimilarity", 
labels = FALSE, hang = 0.04)
```

To identify modules from this gene dendrogram, we can use the `cutreeDynamic` function.

```{r}
Modules <- cutreeDynamic(dendro = geneTree, distM = TOM, deepSplit = 2, pamRespectsDendro = FALSE, minClusterSize = 30)
```

```{r}
table(Modules) 
```

```
Modules
##   0   1   2   3   4   5   6    
##  88 614 316 311 257 235 225 
```

* Here we can see 6 modules were created with the number of genes belong to them. The Label 0 Module is reserved for unassigned genes (genes that do not fit in any module).

#### Module Eigengene Identification

A **ME (Module Eigengene)** is the standardized gene expression profile for a given module.

To identify the Module Eigengene we can call on the expression data into the `moduleEigengenes` function.

```{r}
MElist <- moduleEigengenes(expression.data, colors = ModuleColors) 
MEs <- MElist$eigengenes 
head(MEs)
```

```
            MEblue     MEbrown        MEturquis    MEgreen   MEyellow 
## F2_2   0.013902476  0.0410177922 0.007072125  0.12978459  0.006276361
## F2_3   0.066675342 -0.0009540238 0.072447744 -0.07777835  0.010326534
## F2_14  0.066711912 -0.0841292811 0.062700422 -0.19072152  0.003707524
## F2_15 -0.064480250  0.0909333146 0.050275810  0.04077621 -0.019067137
## F2_19  0.063634038 -0.0709378322 0.016600588 -0.04036901  0.017796637
## F2_20 -0.001201217  0.0653004166 0.049766750  0.10391289 -0.040252274
##           MEgrey   
## F2_2   0.006971934 
## F2_3  -0.016017527 
## F2_14 -0.041321626  
## F2_15 -0.014390509
## F2_19 -0.023401174 
## F2_20  0.113170728 
```

#### Module Merging

Calculate dissimilarity of module eigengenes

```{r}
MEdiss <- 1- cor(MEs)
```

* heirerichal clustering for eigengenes to show closed or overlapped modules if they're exist

* put the cutoff line for more than 0.25 distance which it means 0.75 correlation

* I want to merge each two modules have more than 0.75 correlation

* merging close modules

```{r}
merge <- mergeCloseModules(datExpr, dynamicColor, cutHeight = 0.25, verbose = 3)
mergeColors <- merge$colors
table(mergeColors)
mergedME <- merge$newMEs # Eigengenes of the new merged modules
```

```{r}
sizeGrWindow(12, 9)
plotDendroAndColors(geneTree, cbind(dynamicColor, mergeColors),
                    c("Dynamic Tree Cut", "Merged dynamic"),
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
```

### Quantifying module-trait associations

correlation between each ME and each trait is calculated

```{r}
modulTraitcor <- cor(MEs, datTraits, use = "p") 
```

correlations with p.values and significance values

```{r}
moduleTraitPval <- corPvalueStudent(modulTraitcor, nsamples)
adjPval <- p.adjust(moduleTraitPval)
```

Visualization of the module-trait association, displaying correlations and their p-values

```{r}
sizeGrWindow(12,9)
textMatrix <- paste(signif(modulTraitcor, 2), "\n(",
                    signif(moduleTraitPval, 1), ")", sep = "")
dim(textMatrix) <- dim(modulTraitcor)
par(mar = c(6, 8.5, 3, 1))
labeledHeatmap(Matrix = module.trait.correlation,
xLabels = names(datTraits),
yLabels = names(mergedMEs),
ySymbols = names(mergedMEs),
colorLabels = FALSE,
colors = blueWhiteRed(50),
textMatrix = textMatrix,
setStdMargins = FALSE,
cex.text = 0.4,
zlim = c(-1,1),
main = paste("Module-trait relationships"))
```

* We can use gene mutation status **( or binary traits)** as trait features for analyze the correlation. There would be two columns as mutated and non-mutated states for samples holding mutated genes and non-mutated genes respectively.

### Gene Significance and Module Membership = GS , MM

* We can use the gene significance along with the genes **intramodular connectivity** to identify potential target genes associated with a particular trait of interest. 

> Connectivity - how connected a speficic node is in the network (how many nodes have high correlation with that node). High connectivity indicates a hub gene (central to many nodes).

> Whole Network connectivity - a measure for how well the node is connected throughout the entire system Intramodular connectivity - a measure for how well the node is connected within its assigned module. Also an indicator for how well that node belongs to its module. This is also known as **module membership (MM)**.

Calculate the module membership and the associated p-values. first I should make p.vaues of each gene corr to their own modules

```{r}
MMpvalue <- as.data.frame(corPvalueStudent(as.matrix(geneMODmem), nSamples = nsamples))
#make name of module membership columns
names(geneMODmem) <- paste("MM", modNames, sep = "")
names(MMpvalue) <- paste("p.mm", modNames, sep = "")
```

Calculate the gene significance and associated p-values

```{r}
geneTraitSig <- as.data.frame(cor(datExpr, datTraits, use = "p")) 
GSpvalue <- as.data.frame(corPvalueStudent(as.matrix(geneTraitSig), nsamples)) 
names(geneTraitSig) = paste("GS.", names(datTraits), sep="")
names(GSpvalue) = paste("p.GS.", names(datTraits), sep="")
```

#### Intramodular analysis

Using these two parameters we can identify the hub genes. The **MM > 0.80 and GS > 0.20** cutoffs are used together for detecting hub genes. 

* hub genes are the genes released from a particular module with both criteria.

```{r}
brown = "brown" #the significant module with trait in heatmap
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
```

Now we can extract the hub genes as a table based on two criteria simultaneously. The hub genes have MM values more than 80% membership in the module, and GS values more than 20% significant correlation with a particular trait.

```{r}
Probes = names(datExpr)
inModule <- is.finite(match(moduleColors, brown))
modProbes <- Probes[inModule]
modProbes <- data.frame(modProbes)
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
```

### Network Visualization of Eigengenes

Plot the relationships among the eigengenes and the trait

```{r}
par(cex = 0.9)
plotEigengeneNetworks(MET, "", marDendro = c(0,4,1,2), marHeatmap = c(5,4,1,2), cex.lab = 0.8, xLabelsAngle
= 90)
```

With this heatmap we can identify groups of correlated eigengenes called *meta modules*. Modules with mutual correlations stronger than their correlation with the specified clinical trait would be grouped into a meta module. 

![](D:/my R/my thesis/WGCNA/deg wgcna/tcga deg/plot eigengenes.tiff)












