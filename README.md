# q2-dbbact
A [Qiime2](https://qiime2.org/) plugin for [dbBact](http://dbbact.org)

![wordcloud](https://github.com/amnona/q2-dbbact/blob/main/pics/cfs-wordcloud.jpg)
![enriched barplot](https://github.com/amnona/q2-dbbact/blob/main/pics/enriched_terms.jpg)
![enriched terms](https://github.com/amnona/q2-dbbact/blob/main/pics/terms-table.jpg)

# Features:
* Differential abundance testing using [Calour](https://github.com/biocore/calour) rank-mean differential abundance test (with dsFDR correction).
* dbBact term enrichment from differntial abundance results of qiime2 (i.e. songbird/q2-aldex2/ancom/dacomp or the built in rank-mean test).
* Create a wordcloud of dbBact terms for a given feature table.
* Generate an interactive heatmap visualization for a feature table. The heatmap provides links to dbBact annotations for each ASV.
* Generate Venn diagram for a differential abundance result and a given dbBact term.
* Background dbBact term enrichment analysis for experiments without controls (i.e. what terms are enriched in the bacteria in a given feature table compared to all dbBact experiments of a given type).

# Examples:
## Run the q2-dbBact enrichment pipeline for a given feature table:
