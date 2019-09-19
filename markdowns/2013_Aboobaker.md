# 2013 Aboobaker Time-course

> see: https://genomebiology.biomedcentral.com/articles/10.1186/gb-2014-15-2-r29
> A simple approach to analyzing RNA-seq data would be to input the log-cpm values 
> into a well-established microarray analysis pipeline such as that provided by 
> the limma software package [3, 9]. This would be expected to behave well if the counts 
> were all reasonably large, but it ignores the mean-variance trend for lower counts. 
> The microarray pipeline should behave better if modified to include a mean-variance 
> trend as part of the variance modeling. We have therefore modified the empirical 
> Bayes procedure of the limma package so that the gene-wise variances are squeezed 
> towards a global mean-variance trend curve instead of towards a constant pooled 
> variance. This is similar in principle to the procedure proposed by 
> Sartor et al. [36] for microarray data, except that we model the trend using a 
> regression spline and our implementation allows for the possibility of missing 
> values or differing residual degrees of freedom between genes. We call this 
> strategy limma-trend, whereby the log-cpm values are analyzed as for microarray 
> data but with a trended prior variance.


## Load libraries

```{r}
library(edgeR)
library(limma)
library(dplyr)
library(data.table)
```

## Reading data

```{r}
data <- read.table(file="rawdata/regeneration.counts") 
    # Lowly expressed genes were already filtered by the authors.
data <- log2(data + 0.25); # Apply log2 transformation to cpm values provided by the authors.
groups <- data.frame(
    rownames=colnames(data), 
    time=gsub("regen_[th]_(.+)\\.rep[AB]\\.sense", "\\1", colnames(data), 
    perl=T
    ), section=gsub("regen_([th])_.+", "\\1", colnames(data), perl=T)
)
time <- groups$time
section <- groups$section
levels(groups$section) <- c("head", "tail")
Group <- factor(paste(groups$time, groups$section, sep="."))
```

## Preparing design and contrasts

```{r}
design <- model.matrix(~0 + Group)
all_pairs <- combn(colnames(design), 2)

my.contrasts <- makeContrasts(
    Head.0hvsTail.0h = Group00h.head-Group00h.tail,
    Head.6hvsTail.6h = Group06h.head-Group06h.tail,
    Head.12hvsTail.12h = Group12h.head-Group12h.tail,
    Head.24hvsTail.24h = Group24h.head-Group24h.tail,
    Head.36hvsTail.36h = Group36h.head-Group36h.tail,
    Head.48hvsTail.48h = Group48h.head-Group48h.tail,
    Head.72hvsTail.72h = Group72h.head-Group72h.tail,
    Head.6hvsHead.0h  =  Group06h.head-Group00h.head,
    Head.12hvsHead.6h  =  Group12h.head-Group06h.head,
    Head.24hvsHead.12h  =  Group24h.head-Group12h.head,
    Head.36hvsHead.24h  =  Group36h.head-Group24h.head,
    Head.48hvsHead.36h  =  Group48h.head-Group36h.head,
    Head.72hvsHead.48h  =  Group72h.head-Group48h.head,
    Tail.6hvsTail.0h  =  Group06h.tail-Group00h.tail,
    Tail.12hvsTail.6h  =  Group12h.tail-Group06h.tail,
    Tail.24hvsTail.12h  =  Group24h.tail-Group12h.tail,
    Tail.36hvsTail.24h  =  Group36h.tail-Group24h.tail,
    Tail.48hvsTail.36h  =  Group48h.tail-Group36h.tail,
    Tail.72hvsTail.48h  =  Group72h.tail-Group48h.tail,
    levels=design
)

```

## Performing analysis

```{r}
fit  <- lmFit(data, design)
fit2 <- contrasts.fit(fit, my.contrasts)
fit2 <- eBayes(fit2, trend=TRUE, robust=TRUE)

for (contrast in 1:ncol(my.contrasts)) {
    results <- topTable(fit2, coef=contrast, number=Inf, p.value=0.01);
    results$gene <- rownames(results);
    results <- results %>% dplyr::filter(abs(logFC) >= 1 & adj.P.Val <=0.01);

    comparison <- colnames(my.contrasts)[contrast];
    results_formatted <- data.frame( 
        gene=results$gene,
        logFC=results$logFC,
        adjpvalue=results$adj.P.Val,
        cond_1=strsplit(comparison, "vs")[[1]][1],
        cond_2=strsplit(comparison, "vs")[[1]][2]
    );
    result_df[[contrast]] <- results_formatted;
}
allresults <- rbindlist(result_df)
```

## Save for PlanNET

```{r}
# Relative expression
write.table(
    data.frame(
        Experiment_NAME="2013 Aboobaker Time-course", 
        condition1=sub("(.+)\\.(.+)", "\\2 - \\1", allresults$cond_1),
        condition2=sub("(.+)\\.(.+)", "\\2 - \\1", allresults$cond_2),
        condition_type="Hour - Section", 
        dataset="Consolidated", 
        gene_symbol=allresults$gene, 
        fold_change=allresults$logFC,
        pvalue=allresults$adjpvalue
    ),
    col.names=T,
    row.names=F,
    file="Aboobaker_relative.tbl",
    sep="\t",
    quote=F
)
```
