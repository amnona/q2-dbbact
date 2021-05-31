# q2-dbbact
A [Qiime2](https://qiime2.org/) plugin for [dbBact](http://dbbact.org)

![wordcloud](https://github.com/amnona/q2-dbbact/blob/main/pics/cfs-wordcloud.jpg)
![enriched barplot](https://github.com/amnona/q2-dbbact/blob/main/pics/enriched_terms.jpg)
<!-- ![wordcloud terms](https://github.com/amnona/q2-dbbact/blob/main/pics/terms-table.jpg)
![enriched terms](https://github.com/amnona/q2-dbbact/blob/main/pics/enriched_terms-cfs.jpg)
  -->
 
# Features:
* Differential abundance testing using [Calour](https://github.com/biocore/calour) rank-mean differential abundance test (with dsFDR correction).
* dbBact term enrichment from differntial abundance results of qiime2 (i.e. songbird/q2-aldex2/ancom/dacomp or the built in rank-mean test).
* Create a wordcloud of dbBact terms for a given feature table.
* Generate an interactive heatmap visualization for a feature table. The heatmap provides links to dbBact annotations for each ASV.
* Generate Venn diagram for a differential abundance result and a given dbBact term.
* Background dbBact term enrichment analysis for experiments without controls (i.e. what terms are enriched in the bacteria in a given feature table compared to all dbBact experiments of a given type).

# Examples:
## Run the q2-dbBact enrichment pipeline for a given feature table:
Our input is a feature table and a metadata file with a given column dividing our samples into two groups.

q2-dbBact will detect ASVs different between the two groups, and identify dbBact terms enriched in one of the two groups compared to the other

``` qiime dbbact enrich-pipeline --i-table cfs-merged.qza --m-metadata-file map.cfs.txt --p-field Subject --output-dir cfs-pipeline```

## draw an interactive heatmap
