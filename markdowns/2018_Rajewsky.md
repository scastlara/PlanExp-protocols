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
# Download http://bimsbstatic.mdc-berlin.de/rajewsky/PSCA/dge.txt.gz to rawdata/dge.txt
dgedata <- read.table(file="rawdata/dge.txt")
dgedata <- as.matrix(dgedata)

# Download http://bimsbstatic.mdc-berlin.de/rajewsky/PSCA/R_annotation.txt to rawdata/R_annotation
cellclusters <- read.table(file="rawdata/R_annotation.txt", sep="\t")

# Download http://bimsbstatic.mdc-berlin.de/rajewsky/PSCA/R_pca_seurat.txt to rawdata/R_pca_seurat.txt
cellpca <-  read.table(file="rawdata/R_pca_seurat.txt")

# Download http://bimsbstatic.mdc-berlin.de/rajewsky/PSCA/colors_dataset.txt to rawdata/colors_dataset.txt
clustercolors <- read.table(file="rawdata/colors_dataset.txt", sep="\t", comment.char = "")


dim(dgedata)
dim(cellclusters)
dim(cellpca)
```

## Create Seurat Object

```{r}
so <- CreateSeuratObject(dgedata)
saveRDS(sce, file="~/Solana_SC/sce_solana.RDS")
saveRDS(so, file="~/Solana_SC/so_solana.RDS")
```


## Load clusters and PCA

```{r}
so <- SetDimReduction(so, reduction.type = "pca", slot="cell.embeddings", new.data=as.matrix(cellpca))
so <- SetDimReduction(so, reduction.type = "pca", slot="key", new.data="pca")
so <- SetIdent(so, ident.use = as.vector(cellclusters$V1))
so <- RunTSNE(so, dims.use = 1:50)
TSNEPlot(so)
```

## Plot TSNE

```{r, fig.fullwidth=TRUE, fig.width=10, fig.height=14}}
library(gridExtra)

tsnedata <- data.frame(x=so@dr$tsne@cell.embeddings[,1], y=so@dr$tsne@cell.embeddings[,2], cluster=so@ident)
rownames(clustercolors) <- clustercolors$V1
tsnedata$color <- as.character(clustercolors[as.character(tsnedata$cluster),]$V2)
cols <- as.character(tsnedata$color)
names(cols)<- as.character(tsnedata$cluster)
g_legend<-function(a.gplot){
     tmp <- ggplot_gtable(ggplot_build(a.gplot))
     leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
     legend <- tmp$grobs[[leg]]
}

legend <- ggplot(tsnedata) + 
   geom_point(aes(x=x, y=y, color=cluster)) + 
   scale_color_manual(values=cols)

theplot <- ggplot(tsnedata) +
  geom_point(aes(x=x, y=y, color=cluster)) + 
  scale_color_manual(values=cols) +
  theme_minimal() +
  theme(legend.position = "none") 

plot(legend)
plot(theplot)
```

## Data normalization & scaling

```{r}
so <- NormalizeData(so, normalization.method = "LogNormalize")
so <- ScaleData(object = so, vars.to.regress = c("nUMI"))
```

## Find Markers

```{r}
so.markers <- FindAllMarkers(so, test.use="roc")
print( so.markers %>% group_by(cluster) %>% top_n(10, avg_logFC) )
```

## Save for PlanNET

```{r}
source("src/create_sc_files.R");
library(SingleCellExperiment);
library(data.table)

sce <- Convert(so, to="sce");
write.cell.file(sce, "./Rajewsky_cellfile.tbl", c("ident", "orig.ident"));
write.exp.file(sce, "./Rajewsky_expfile.tbl");
write.diff.file.seurat(so, "./Rajewsky_difffile.tbl", 1e-5);
```

