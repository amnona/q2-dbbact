# q2-dbbact
A [Qiime2](https://qiime2.org/) plugin for [dbBact](http://dbbact.org)

![wordcloud](https://github.com/amnona/q2-dbbact/blob/main/pics/cfs-wordcloud.jpg)
![enriched barplot](https://github.com/amnona/q2-dbbact/blob/main/pics/enriched_terms.jpg)
![heatmap](https://github.com/amnona/q2-dbbact/blob/main/pics/heatmap.jpg)
<!-- ![wordcloud terms](https://github.com/amnona/q2-dbbact/blob/main/pics/terms-table.jpg)
![enriched terms](https://github.com/amnona/q2-dbbact/blob/main/pics/enriched_terms-cfs.jpg)
  -->
 
# Features:
* Differential abundance testing using [Calour](https://github.com/biocore/calour) rank-mean differential abundance test (with [dsFDR](https://escholarship.org/content/qt3j68q5n7/qt3j68q5n7_noSplash_e7ad1cf405f67b9cef0e5a99c1804fd5.pdf) correction).
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

## Draw an interactive heatmap
This creates a zoomable heatmap with a list of dbBact annotation for each bacteria that is clicked. Useful for exploring your sequencing results and getting a feeling for what is going on (contaminations, bacterial sources, groups of samples, etc.)

Our input is a feature table and a metadata file with a given column dividing our samples into two groups.

```qiime dbbact heatmap --i-table cfs-table.qza --i-repseqs cfs-rep-seqs.qza --i-taxonomy cfs-taxonomy.qza --m-metadata-file map.cfs.txt --p-sort-field Subject --o-visualization heatmap-cfs```

![heatmap](https://github.com/amnona/q2-dbbact/blob/main/pics/heatmap.jpg)

## Draw a dbBact terms wordcloud for the set of bacteria in a feature-table
The wordcloud is created for all the bacteria in the feature table.

The output wordcloud words are dbBact terms associated with the bacteria. The word size corresponds to the F-score (recall and precision) of the term. Blue terms are positively associated (i.e. appear in COMMON/DOMINANT/HIGHER IN annotations) where as red terms (preceeded by a "-") are negatively associated (i.e. appear in LOWER IN annotations).

![wordcloud](https://github.com/amnona/q2-dbbact/blob/main/pics/cfs-wordcloud.jpg)

## Identify differentially abundant bacteria between two sample groups
q2-dbBact utilizes the non-parametric (permutation based) Calour diff_abundance() function. By default it uses a rank-mean test with dsFDR multiple hypothesis correction.

The test can also be performed as a paired test using an additional metadata pair-field (permutations are performed only between samples sharing the same pair-field value).

```qiime dbbact diff-abundance --i-table cfs-merged.qza --m-metadata-file map.cfs.txt --p-field Subject --p-alpha 0.1 --p-val1 Patient --p-val2 Control --o-diff diff-cfs-dsfdr```

## Identify and plot enriched dbBact terms between two groups of bacteria
Performed on the output of a differential-abundance test. q2-dbBact supports the following formats:
* [songbird](https://github.com/biocore/songbird)
* [ancom](https://github.com/qiime2/q2-composition)
* [q2-aldex2](https://library.qiime2.org/plugins/q2-aldex2/24/)
* dbBact diff-abundance
* any tsv file

This command identifies dbBact terms the are significantly more associated with bacteria from one group compared to the other
```qiime dbbact enrichment --i-diff diff-cfs-dsfdr.qza --p-source dsfdr --o-enriched enriched-cfs-dsfdr```

The output can be visualized (and the complete table saved) using the visualization command:
```qiime dbbact plot-enrichment --i-enriched enriched-cfs-dsfdr.qza --o-visualization barplot-enriched-cfs-dsfdr --p-labels CFS Control```
![enriched barplot](https://github.com/amnona/q2-dbbact/blob/main/pics/enriched_terms.jpg)

