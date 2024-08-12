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

The **WGCNA** pipeline is expected an input matrix of normalized expression values including samples in columns and gene names on rows. There is no limitation for the methods exploring the expression values; *RNA-Seq* or *microarray* methods. We can use *GEO* or *TCGA* expression profiles for this analysis. 

#### Installation

```{r}
BiocManager::install("WGCNA")
library(WGCNA)
```

Usually we need to rotate or *transpose* the rows with the columns in the matrix using `t` function.

> rows => samples , columns => genes




