#first I should make p.vaues of each gene corr to their own modules
MMpvalue <- as.data.frame(corPvalueStudent(as.matrix(geneMODmem), nSamples = nsamples))
#make name of module membership columns
names(geneMODmem) <- paste("MM", modNames, sep = "")
names(MMpvalue) <- paste("p.mm", modNames, sep = "")

geneTraitSig <- as.data.frame(cor(datExpr, datTraits, use = "p")) #corr of each gene with trait
GSpvalue <- as.data.frame(corPvalueStudent(as.matrix(geneTraitSig), nsamples)) #pvalues of each gene with trait
names(geneTraitSig) = paste("GS.", names(datTraits), sep="")
names(GSpvalue) = paste("p.GS.", names(datTraits), sep="")

#####Intramodular analysis#####
#################
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