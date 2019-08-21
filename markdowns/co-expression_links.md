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
 network_rajewsky_FILTERED.human.simple.tbl \
 rawdata/reactome/ReactomePathways.gmt \
 > network_rajewsky_FILTERED.human.genes.reactome.csv
```

* Some stats

```{sh}
# Number of Links
tail +2 network_rajewsky_FILTERED.human.genes.reactome.csv | wc

# Number of genes in Links
tail +2 network_rajewsky_FILTERED.human.genes.reactome.csv | gawk -F ',' '{print $1"\n"$3}' | sort | uniq | wc

# Number of links with annotated Reactome
egrep -c "R-HSA" network_rajewsky_FILTERED.human.genes.reactome.csv

# Number of genes in links with annotated Reactome
egrep "R-HSA" network_rajewsky_FILTERED.human.genes.reactome.csv | gawk -F ',' '{print $1"\n"$3}' | sort | uniq | wc

# Number of annotated Reactome pathways
perl -ne '
    chomp; 
    @cols = split /,/; 
    @r = split /;/, $cols[-1]; 
    print join("\n", @r), "\n"
' network_rajewsky_FILTERED.human.genes.reactome.csv | \
  egrep -v "\bNA\b|ReactomePathways" | \
  wc

# Number of unique annotated Reactome pathways
perl -ne '
    chomp; 
    @cols = split /,/; 
    @r = split /;/, $cols[-1]; 
    print join("\n", @r), "\n"
' network_rajewsky_FILTERED.human.genes.reactome.csv | \
  egrep -v "\bNA\b|ReactomePathways" | \
  sort | uniq | wc
```

| Total Links | Total Genes | Links with Reactome | Genes in links with Reactome | Reactome Pathways | Unique Reactome Pathways |
| :---------: | :---------: | :-----------------: | :-------------------------:  | :---------------: | :----------------------: |
|    37,246   |   15,196    |       3,635         |             1,855            |      17,257       |           638            |


## Reddien Cell Atlas

