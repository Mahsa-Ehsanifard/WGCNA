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







