# 2018 Rajewsky Cell Atlas

## Load Libraries

```{r}
library(Seurat)
library(ggplot2)
library(dplyr)
library(SingleCellExperiment)
```

## Load data

```{r}
# Download https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE111764&format=file&file=GSE111764%5FPrincipalClusteringDigitalExpressionMatrix%2Edge%2Etxt%2Egz to rawdata/GSE111764_PrincipalClusteringDigitalExpressionMatrix.dge.txt
dge <- read.table(file="rawdata/GSE111764_PrincipalClusteringDigitalExpressionMatrix.dge.txt")
cellinfo <- read.table("cells.csv", sep=",")
```


## Data normalization & scaling

```{r}
so <- NormalizeData(object = so, normalization.method = "LogNormalize", scale.factor = 10000)
so <- FindVariableGenes(
    object = so, 
    mean.function = ExpMean, 
    dispersion.function = LogVMR,
    x.low.cutoff = 0.2, 
    x.high.cutoff = 15, 
    y.cutoff = -0.5
)
print(length(x = so@var.genes))
so <- ScaleData(object = so, vars.to.regress = c("nUMI"))
```

## Dimensionality reduction

```{r}
so <- RunPCA(object = so, pc.genes = so@var.genes, do.print = TRUE, pcs.print = 1:5, genes.print = 5, pcs.compute = 150)
so <- RunTSNE(object = so, dims.use = 1:150, do.fast = TRUE)
so <- SetIdent(so, ident.use=paste(cellinfo[so@cell.names,]$Cluster.ID, cellinfo[so@cell.names,]$Major.cluster.description, sep=" "))
TSNEPlot(so)
```

## Save for PlanExp

```{r}
so@meta.data$orig.ident <- cellinfo[so@cell.names, ]$Section
sce <- Convert(so, to="sce");
write.cell.file(sce, "./Reddien_cellfile.tbl", c("ident", "orig.ident"));
write.exp.file(sce, "./Reddien_expfile.tbl");
write.diff.file.seurat(so, "./Reddien_difffile.tbl", 1e-5);

all.markers <- FindAllMarkers(so, test.use="roc")
write.table(all.markers, sep=",", quote=F, file="Reddien_markers.csv", row.names=F)
tsnedata <- data.frame(x=so@dr$tsne@cell.embeddings[,1], y=so@dr$tsne@cell.embeddings[,2])
write.table(file="Reddien_tsne.tbl", sep="\t", tsnedata, col.names=F, quote=F)
```
