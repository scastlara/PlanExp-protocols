# PlanExp-protocols
Protocols for the analyses performed for the PlanExp web application.

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

## Gene Co-expression inference

The protocol for inferring gene co-expression relationships using GENIE3 can be found [here](markdowns/co-expression_genie3.md).


