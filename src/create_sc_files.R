#!/usr/bin/R
# Script to write the necessary files to upload 
# a SingleCell experiment on neo4j and mysql
# for the PlanNET website

# Writes cell clustering and conditions file
write.cell.file <- function(sce, filename, condnames) {
	outdf <- colData(sce)[,condnames]
	outdf$cell <- rownames(outdf);
	outdf <- outdf[c("cell", condnames)];
	write.table(
		outdf,
		sep="\t",
		file=filename,
		quote=F,
		col.names=T,
		row.names=F,
	);
}

# Writes expression file for sce
write.exp.file <- function(sce, filename) {
	write.table(
		as.matrix(assays(sce)$logcounts), 
		sep="\t", 
		file=filename, 
		quote=F
	);
}

# Write differential expression information table
write.diff.file <- function(sce, filename) {
	markers.up <- findMarkers(sce, sce$clusters, direction="up");
	outdf <- data.frame(
		matrix(vector(), 0, 6,
		dimnames=list(c(), c("cluster.focus", "cluster.vs", "gene", "logfc", "fdr", "top.position"))),
		stringsAsFactors=F
	);
	outdf <- outdf[-1,];
	for (cluster in unique(sce$clusters)) { 
		# Foreach focus cluster
		setdf <- data.frame(
			matrix(vector(), 0, 6,
			dimnames=list(c(), c("cluster.focus", "cluster.vs", "gene", "logfc", "fdr", "top.position"))),
			stringsAsFactors=F
		); 
		setdf <- setdf[-1,];
		marker.set <- markers.up[[cluster]];
		fc <- grep("logFC", colnames(marker.set));
		fc_cols <- colnames(marker.set)[fc];
		# Foreach vs cluster
		for (col in fc_cols) {
			vsname <- sub("logFC.", "", col);
                        print(cluster);
                        
			setdf <- rbind(
				setdf, 
				data.frame(
					cluster.focus=cluster, 
					cluster.vs=vsname, 
					gene=marker.set$Gene, 
					logfc=marker.set[,col], 
					fdr=marker.set$FDR,
					top.position=marker.set$Top)
			);
		}
		outdf <- rbind(outdf, setdf);
		#setdf <- data.frame(cluster.focus=cluster, cluster.vs=) 
	}
	return(outdf);
}



write.diff.file.seurat <- function(so, filename, significance) {
	library(data.table);
	combinations <- combn(sort(unique(so@ident)), 2);
	dflist <- as.list(ncol(combinations))
	for (i in 1:ncol(combinations)) {
  		pair <- combinations[,i];
		markers <- FindMarkers(so, ident.1=pair[1], ident.2=pair[2], logfc.threshold=0.2, min.pct=0.2);
		markers$gene <- rownames(markers);
		markers <- markers %>% filter(p_val_adj < significance);
		if (dim(markers)[1] != 0) {
		  dflist[[i]] <- data.frame(
			  cluster.focus=pair[1],
			  cluster.vs=pair[2],
			  gene=markers$gene,
			  logfc=markers$avg_logFC,
			  pvalue=markers$p_val_adj,
			  top.position="NA"
			);
		}
		
	}
	setdf <- rbindlist(dflist)
	write.table(
		setdf, 
		sep="\t", 
		file=filename, 
		quote=F
	);

}