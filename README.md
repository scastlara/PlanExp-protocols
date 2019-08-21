# PlanExp-protocols
Protocols of the analyses performed for the PlanExp web application.

## Download data

First, some data files that are too big to store in the git repository must be downloaded. Use the following commands.

```{sh}
cd PlanExp-protocols
chmod a+x src/*
src/download_data.sh
```

## RNA-Seq analyses

The analyses for the RNA-seq experiments (traditional RNA-seq and Single-Cell RNA-seq) can be found in the markdown/ directory.

* [2013 Aboobaker Time-course](markdowns/2013_Aboobaker.md)
* [2018 Rajewsky Cell Atlas](markdowns/2018_Rajewsky.md)
* [2018 Reddien Cell Atlas](markdowns/2018_Reddien.md)

## Gene co-expression network inference

The protocol for inferring gene co-expression relationships using GRNBoost can be found [here](markdowns/co-expression_links.md).


## References

* Castillo-Lara, S. and Abril, J. F. (2018). PlanNET: homology-based predicted interactome for multiple planarian transcriptomes. Bioinformatics, 34(6), 1016–1023.
* Fincher, C. T. and others. (2018). Cell type transcriptome atlas for the planarian Schmidtea mediterranea. Science, 360(6391).
* Kao, D. et al. (2013). The planarian regeneration transcriptome reveals a shared but temporally shifted regulatory program between opposing head and tail scenarios. BMC Genomics, 14(1), 797.
* Labbé, R. M. et al. (2012). A comparative transcriptomic analysis reveals conserved features of stem cell pluripotency in planarians and mammals. Stem Cells, 30(8), 1734–1745.
* Plass, M. et al. (2018). Cell type atlas and lineage tree of a whole complex animal by single-cell transcriptomics. Science, 360(6391).
* Rodríguez-Esteban, G. et al. (2015). Digital gene expression approach over multiple RNA-Seq data sets to detect neoblast transcriptional changes in Schmidtea mediterranea. BMC Genomics, 16(1), 361.
* Sandmann, T. et al. (2011). The head-regeneration transcriptome of the planarian Schmidtea mediterranea. Genome Biology, 12(8), R76.
