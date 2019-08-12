#!/usr/bin/bash

# Downloads missing raw data 

echo "Downloading 2018 Rajewsky data";
wget "http://bimsbstatic.mdc-berlin.de/rajewsky/PSCA/dge.txt.gz" -O rawdata/dge.txt;
echo "Done.";

echo "Downloading 2018 Reddien data";
wget "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE111764&format=file&file=GSE111764%5FPrincipalClusteringDigitalExpressionMatrix%2Edge%2Etxt%2Egz" -O "rawdata/GSE111764_PrincipalClusteringDigitalExpressionMatrix.dge.txt";
echo "Done.";

echo "\nAll data downloaded to ./rawdata directory.";
