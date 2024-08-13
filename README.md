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











