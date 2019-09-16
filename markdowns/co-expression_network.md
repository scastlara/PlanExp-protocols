# Gene co-expression network prediction

## Rajewsky Cell Atlas

Predicting gene co-expression network. This can take many hours...

```{sh}
src/grn_rajewsky.py 2> grn_rajewsky.log &
```
> 50,264,041 co-expression links.

Merging all the files.

```{sh}
cat *.csv > network_rajewsky_ALL.csv
```

Applying filters to putative regulatory links.

```{sh}
src/filter_grn.py network_rajewsky_ALL.csv > network_rajewsky_FILTERED.csv
```

```{sh}
wc -l network_rajewsky_FILTERED.csv
gawk -F ',' '{print $1"\n"$3}' network_rajewsky_FILTERED.csv | sort | uniq | wc
```
> 37,246 filtered co-expression links.
> 15,196 genes.

Labeling genes by component.

```{sh}
python3 filter_connected_components.py network_rajewsky_FILTERED.csv > network_rajewsky_FILTERED.components.csv

gawk -F ',' '{print $4}' network_rajewsky_FILTERED.components.csv  | sort | uniq -c | sort -g -k2
```

```  
37234 0
    1 1
    1 2
    1 3
    1 4
    1 5
    1 6
    1 7
    1 8
    2 9
    1 10
    1 11
```

> There are 12 connected components.

Annotate planarian genes with Human homologs.

```{sh}
perl -e '
    use strict; 
    open my $fh, "<", "rawdata/dresden-human.tbl"; 
    my %d2h = (); 
    while(<$fh>) {
        chomp; 
        my ($dresden, $human) = split /\s+/; 
        $d2h{$dresden} = $human unless exists $d2h{$dresden};
    };
    close($fh);
    while(<>) {
        chomp; 
        my ($a, $s, $b) = split /,/; 
        my $ha = exists $d2h{$a} ? $d2h{$a} : "NA"; 
        my $hb = exists $d2h{$b} ? $d2h{$b} : "NA"; 
        print "$a,$s,$b,$ha,$hb\n"
    }'  network_rajewsky_FILTERED.csv > network_rajewsky_FILTERED.human.csv

egrep -v "\bNA\b" network_rajewsky_FILTERED.human.csv | gawk -F ',' '{print $4,$5}' > network_rajewsky_FILTERED.human.simple.tbl
```

```{sh}
perl -e '
    use strict;
    use warnings;
    open my $fh, "<", "rawdata/dresden-x-smedg.tsv";
    my %d2g = ();
    while (<$fh>) {
        chomp;
        my ($d, $g, $n) = split /\t/;
        $d =~ s/_0_\d+/_0/;
        $d2g{$d} = [$g, $n =~ m/\w/ ? $n : "NA"];
    }
    close($fh);
    print "Regulator_ContigId,Score,Target_ContigId,RegulatorHomolog,TargetHomolog,RegulatorGeneId,RegulatorGeneName,TargetGeneId,TargetGeneName\n";
    while (<>) {
        chomp;
        my ($d1, $sc, $d2, $h1, $h2) = split /,/;
        my $g1 = exists $d2g{$d1} ? $d2g{$d1} : ["NA", "NA"]; 
        my $g2 = exists $d2g{$d2} ? $d2g{$d2} : ["NA", "NA"];
        print "$d1,$sc,$d2,$h1,$h2,$g1->[0],$g1->[1],$g2->[0],$g2->[1]\n";

    } 
' network_rajewsky_FILTERED.human.csv | more
```


Find co-expression links that appear in the same pathway in reactome

```{sh}
src/find_in_reactome.py \
 network_rajewsky_FILTERED.human.genes.tbl \
 rawdata/reactome/ReactomePathways.gmt \
 > network_rajewsky_FILTERED.human.genes.reactome.csv
```

* Some stats

```{sh}
# Number of Links
tail +2 network_rajewsky_FILTERED.human.genes.reactome.csv | wc

# Number of genes in Links
tail +2 network_rajewsky_FILTERED.human.genes.reactome.csv | gawk -F '\t' '{print $1"\n"$3}' | sort | uniq | wc

# Number of links with annotated Reactome
egrep -c "R-HSA" network_rajewsky_FILTERED.human.genes.reactome.csv

# Number of genes in links with annotated Reactome
egrep "R-HSA" network_rajewsky_FILTERED.human.genes.reactome.csv | gawk -F '\t' '{print $1"\n"$3}' | sort | uniq | wc

# Number of annotated Reactome pathways
perl -ne '
    chomp; 
    @cols = split /\t/; 
    @r = split /;/, $cols[-1]; 
    print join("\n", @r), "\n"
' network_rajewsky_FILTERED.human.genes.reactome.csv | \
  egrep -v "\bNA\b|ReactomePathways" | \
  wc

# Number of unique annotated Reactome pathways
perl -ne '
    chomp; 
    @cols = split /\t/; 
    @r = split /;/, $cols[-1]; 
    print join("\n", @r), "\n"
' network_rajewsky_FILTERED.human.genes.reactome.csv | \
  egrep -v "\bNA\b|ReactomePathways" | \
  sort | uniq | wc
```

| Total Links | Total Genes | Links with Reactome | Genes in links with Reactome | Reactome Pathways | Unique Reactome Pathways |
| :---------: | :---------: | :-----------------: | :-------------------------:  | :---------------: | :----------------------: |
|    37,246   |   15,196    |       3,635         |             1,855            |      18671       |           792            |



----

## Reddien Cell Atlas


Predicting gene co-expression network. This can take many hours...

```{sh}
src/grn_reddien.py 2> grn_reddien.log &
```
> 47,568,574 co-expression links.

Merging all the files.

```{sh}
cat *.csv > network_reddien_ALL.csv
```

Applying filters to putative regulatory links.

```{sh}
src/filter_grn.py network_reddien_ALL.csv > network_reddien_FILTERED.csv

# Graph stats:
## Number of Links: 47568574
## Number of Regulators: 1754
## Number of Targets: 26561
## Filtered interactions: 36014
```


```{sh}
wc -l network_reddien_FILTERED.csv
gawk -F ',' '{print $1"\n"$3}' network_reddien_FILTERED.csv | sort | uniq | wc
```

> 36,014 interactions
> 11,408 genes


Labeling genes by component.

```{sh}
python3 filter_connected_components.py network_reddien_FILTERED.csv > network_reddien_FILTERED.components.csv

gawk -F ',' '{print $4}' network_reddien_FILTERED.components.csv  | sort | uniq -c | sort -g -k2
```

```
35934 0
      3 1
      2 2
      6 3
      1 4
      1 5
      1 6
      1 7
      1 8
      1 9
      1 10
      1 11
      1 12
      1 13
      1 14
      5 15
      2 16
      3 17
      1 18
      2 19
      2 20
      2 21
      2 22
      1 23
      1 24
      7 25
      3 26
      1 27
      4 28
      1 29
      1 30
      1 31
      1 32
      2 33
      3 34
      1 35
      1 36
      4 37
      1 38
      2 39
      1 40
      1 41
      2 42
```  


> There are 43 connected components.

Annotate planarian genes with Human homologs.

```{sh}
perl -e '
    use strict; 
    open my $fh, "<", "rawdata/dresden-human.tbl"; 
    my %d2h = (); 
    while(<$fh>) {
        chomp; 
        my ($dresden, $human) = split /\s+/; 
        $d2h{$dresden} = $human unless exists $d2h{$dresden};
    };
    close($fh);
    while(<>) {
        chomp; 
        my ($a, $s, $b) = split /,/; 
        my $ha = exists $d2h{$a} ? $d2h{$a} : "NA"; 
        my $hb = exists $d2h{$b} ? $d2h{$b} : "NA"; 
        print "$a,$s,$b,$ha,$hb\n"
    }'  network_reddien_FILTERED.csv > network_reddien_FILTERED.human.csv

egrep -v "\bNA\b" network_reddien_FILTERED.human.csv | gawk -F ',' '{print $4"\t"$5}' > network_reddien_FILTERED.human.simple.tbl
```

```{sh}
perl -e '
    use strict;
    use warnings;
    open my $fh, "<", "../../PlanExp/rawdata/dresden-x-smedg.tsv";
    my %d2g = ();
    while (<$fh>) {
        chomp;
        my ($d, $g, $n) = split /\t/;
        $d =~ s/_0_\d+/_0/;
        $d2g{$d} = [$g, $n =~ m/\w/ ? $n : "NA"];
    }
    close($fh);
    print "Regulator_ContigId,Score,Target_ContigId,RegulatorHomolog,TargetHomolog,RegulatorGeneId,RegulatorGeneName,TargetGeneId,TargetGeneName\n";
    while (<>) {
        chomp;
        $_ =~ s/_0_1/_0/g;
        my ($d1, $sc, $d2, $h1, $h2) = split /,/;
        my $g1 = exists $d2g{$d1} ? $d2g{$d1} : ["NA", "NA"]; 
        my $g2 = exists $d2g{$d2} ? $d2g{$d2} : ["NA", "NA"];
        print "$d1,$sc,$d2,$h1,$h2,$g1->[0],$g1->[1],$g2->[0],$g2->[1]\n";

    } 
' network_reddien_FILTERED.human.csv > network_reddien_FILTERED.human.genes.csv
```


Find co-expression links that appear in the same pathway in reactome

```{sh}


src/find_in_reactome.py \
    rawdata/ReactomePathways.gmt \
    network_reddien_FILTERED.human.genes.csv \
    > network_reddien_FILTERED.human.genes.reactome.csv

```






* Some stats

```{sh}
# Number of Links
tail +2 network_reddien_FILTERED.human.genes.reactome.csv | wc

# Number of genes in Links
tail +2 network_reddien_FILTERED.human.genes.reactome.csv | gawk -F '\t' '{print $1"\n"$3}' | sort | uniq | wc

# Number of links with annotated Reactome
egrep -c "R-HSA" network_reddien_FILTERED.human.genes.reactome.csv

# Number of genes in links with annotated Reactome
egrep "R-HSA" network_reddien_FILTERED.human.genes.reactome.csv | gawk -F '\t' '{print $1"\n"$3}' | sort | uniq | wc

# Number of annotated Reactome pathways
perl -ne '
    chomp; 
    @cols = split /\t/; 
    @r = split /;/, $cols[-1]; 
    print join("\n", @r), "\n"
' network_reddien_FILTERED.human.genes.reactome.csv | \
  egrep -v "\bNA\b|ReactomePathways" | \
  wc

# Number of unique annotated Reactome pathways
perl -ne '
    chomp; 
    @cols = split /\t/; 
    @r = split /;/, $cols[-1]; 
    print join("\n", @r), "\n"
' network_reddien_FILTERED.human.genes.reactome.csv | \
  egrep -v "\bNA\b|ReactomePathways" | \
  sort | uniq | wc
```


| Total Links | Total Genes | Links with Reactome | Genes in links with Reactome | Reactome Pathways | Unique Reactome Pathways |
| :---------: | :---------: | :-----------------: | :-------------------------:  | :---------------: | :----------------------: |
|    36,014   |   11,408    |       3,490         |             1,601            |      17,215       |           777            |



------

## Comparison

```{sh}
perl -npe 's/_0_1/_0/g' network_reddien_FILTERED.csv | gawk -F ',' '{print $1"\t"$3"\t"$2}' > reddien.tbl
gawk -F ',' '{print $1"\t"$3"\t"$2}' network_rajewsky_FILTERED.csv > rajewsky.tbl

graphcompare -in rajewsky.tbl,reddien.tbl -fmtin TBL -v comparison-venn.svg -u comparison-upset > /dev/null 

perl -e '
    use strict;
    use warnings;
    my $file1 = shift @ARGV;
    my $file2 = shift @ARGV;
    print "$file1\n";
    my %data = ();
    open my $fh1, "<", $file1 or die "Cant open file $file1: $!\n";
    while(<$fh1>) {
        chomp; 
        my ($d1, $sc, $d2) = split /,/;
        $data{$d1} = () unless exists $data{$d1};
        $data{$d1}{$d2} = [$sc, "NA"];
    }
    close($fh1);
    open my $fh2, "<", $file2 or die "Cant open file $file2: $!\n";
    while(<$fh2>) {
        chomp;
        $_ =~ s/_0_1/_0/g;
        my ($d1, $sc, $d2) = split /,/;
        $data{$d1} = () unless exists $data{$d1};
        if (not exists $data{$d1}{$d2}) {
            $data{$d1}{$d2} = ["NA", $sc];
        } else {
            $data{$d1}{$d2}->[1] = $sc;
        }
        
    }
    close($fh2);
    foreach my $d1 (keys %data) {
        foreach my $d2 (keys %{ $data{$d1} }) {
            print "$d1\t$d2\t$data{$d1}{$d2}[0]\t$data{$d1}{$d2}[1]\n";
        }
    }
' ../rajewsky_grn/network_rajewsky_ALL.csv ../reddien_grn/network_reddien_ALL.csv > score_comparison.tbl
```

> In R terminal...

```{r}
library(dplyr)
library(ggplot2)

data <- read.table(file="score_comparison.tbl", header=F)
data_na <- na.omit(data)

png("score_comparison.png")
ggplot(data_zero, aes(x=V3, y=V4)) + 
    geom_point(alpha=0.5) + 
    geom_rug() + 
    theme_light() + 
    xlab("2018 Rajewsky weight") + 
    ylab("2018 Reddien weight")
dev.off()
```


## Saving for PlanExp


```{sh}
python3 src/format_grn.py \
    "2018 Rajewsky Cell Atlas:network_rajewsky_FILTERED.human.genes.reactome.csv" \
    "2018 Reddien Cell Atlas:network_reddien_FILTERED.human.genes.reactome.csv" \
    | sort -k1 -k2 -k3 > GRN_FOR_PLANNET.tbl
```
