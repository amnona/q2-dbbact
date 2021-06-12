import qiime2.plugin
from qiime2.plugin import (SemanticType, Str, Int, Float, Choices, List,
                           MetadataColumn, Categorical, Numeric, Metadata, Bool, Visualization)
from q2_types.feature_table import (
    FeatureTable, Frequency)
from q2_types.feature_data import (FeatureData, Differential, Sequence, Taxonomy, )

from . import __version__
from .methods import single_enrichment, enrichment, embed_seqs, bg_term_enrichment, diff_abundance, enrich_pipeline, trim_primers
from .visualizations import draw_wordcloud_vis, plot_enrichment, heatmap, venn


_citation = ('manuscript in preparation')

_short_description = "A Qiime2 plugin for the dbBact bacterial knowledge-base (dbbact.org)"


# options for the diff_abundance plugin
_statistical_tests = ['meandiff', 'stdmeandiff']
_transform_functions = ['rankdata', 'log2data', 'normdata', 'binarydata', 'none']
_fdr_method = ['dsfdr', 'bhfdr', 'byfdr', 'filterBH']


plugin = qiime2.plugin.Plugin(
    name='dbbact',
    version=__version__,
    website='https://github.com/amnona/q2_dbbact',
    package='q2_dbbact',
    short_description='dbBact plugin for Qiime2',
    description=('A Qiime2 plugin for the dbBact bacterial knowledge-base (dbbact.org)'),
    citation_text=_citation
)


plugin.pipelines.register_function(
    function=enrich_pipeline,
    inputs={'repseqs': FeatureData[Sequence],
            'table': FeatureTable[Frequency]},
    # outputs=[('all_wordcloud', Visualization),
    #          ('diff_ASVs_table', FeatureData[Differential]),
    #          ('diff_ASVs_heatmap', Visualization),
    #          ('enriched_terms_table', FeatureData[Differential]),
    #          ('enriched_terms_barplot', Visualization),
    #          ('enriched_terms_heatmap', Visualization),
    #          ('enriched_term_venn', Visualization)],
    outputs=[
        ('trimmed_table', FeatureTable[Frequency]),
        ('all_wordcloud', Visualization),
        ('diff_asv_table', FeatureData[Differential]),
        # ('diff_asv_heatmap', Visualization),
        ('enriched_terms_table', FeatureData[Differential]),
        ('enriched_terms_barplot', Visualization),
        # ('enriched_terms_heatmap', Visualization),
        # ('enriched_term_venn', Visualization),
    ],
    parameters={
        'metadata': Metadata,
        'method': Str % Choices('groups', 'correlation'),
        'field': Str,
        'val1': List[Str],
        'val2': List[Str],
        'pair_field': Str,
        'random_seed': Int,
        'sig_threshold': Float,
        'maxid': Int,
        'statistical_test': Str % Choices(_statistical_tests),
        'transform_function': Str % Choices(_transform_functions),
        'permutations': Int,
        'fdr_method': Str % Choices(_fdr_method),
        'attack': Bool,
    },
    input_descriptions={
        'table': 'The feature table to perform differential abundance on (needs to be normalized or rarified)',
        'repseqs': 'The corresponding representative sequences (if not embedded in the table using the --p-no-hashed-feature-ids parameter when denoising)'
    },
    parameter_descriptions={
        'metadata': 'Metadata (mapping) file for the table',
        'field': 'The metadata field on which to perform the differential abundance test',
        'pair_field': 'If supplied, perform paired differential abundance testing pairing samples in the pair_field metadata field',
        'val1': 'Metadata values (in field) for samples in group1. Can supply multiple values',
        'val2': 'Metadata values (in field) for samples in group2. Can supply multiple values. If not provided, take all samples not in group1',
        'statistical_test': 'The statistic used for the effect size calculation. "meandiff" is the difference in the means between the two groups, '
                            '"stdmeandiff" is the difference in the means normalized by the standard deviation.',
        'transform_function': 'The transformation to apply to the frequency data (for each feature) before calculating the effect size. options: '
                              '"rankdata": rank each feature across samples, '
                              '"log2data": log2 transform each frequency, '
                              '"normdata": normalize to constant sum for each feature, '
                              '"binarydata": convert to present/absent, '
                              '"none": do not apply any transformation',
        'permutations': 'Number of permutations to calculate (higher number is more exact but slower)',
        'random_seed': 'The random seed to use (for fully reproducible results, since test is permutation based)',
        'fdr_method': 'The multiple hypothesis correction method to use. options: '
                      '"dsfdr": The discrete FDR correction (stronger than BHfdr if there are many low prevalence features) (Jiang et al 2017, https://doi.org/10.1128/mSystems.00092-17), '
                      '"bhfdr": Benhaminy-Hochberg FDR correction, '
                      '"byfdr": Benhaminy-Yekutieli FDR correction, '
                      '"filterBH": Low frequency filtering based correction (see Jiang et al 2017, https://doi.org/10.1128/mSystems.00092-17).',
        'method': '"groups" to compare term enrichemnent between significantly enriched sequences in both directions, "correlation" to detect dbbact terms significanly correlated/anti-correlated with the effect size.',
        'maxid': 'The maximal dbBact annotation id to use (to enable replication of results after new annotations are added to dbBact',
        'attack': 'Attack mode',
    },
    output_descriptions={'all_wordcloud': 'Wordcloud of all the features in the input table',
                         'diff_asv_table': 'Table of the differentially abundant ASVs between to the groups acoording to metadata field',
                         },
    name='dbBact term enrichment for differential abundance results',
    description=("Identify dbBact terms enriched in results of differential abundance (terms significantly more represented in either of the differential abundance groups or correlated with effect size)")
)


plugin.methods.register_function(
    function=enrichment,
    inputs={'diff': FeatureData[Differential],
            'repseqs': FeatureData[Sequence],
            },
    outputs=[('enriched', FeatureData[Differential])],
    parameters={
        'source': Str % Choices(['dsfdr', 'aldex2', 'dacomp', 'songbird', 'ancom', 'tsv']),
        'method': Str % Choices('groups', 'correlation'),
        'random_seed': Int,
        'sig_threshold': Float,
        'diff_tsv': Str,
        'ancom_stat': Str,
        'maxid': Int,
        'attack': Bool,
    },
    input_descriptions={
        'diff': 'Result of qiime2 differential abundance plugin (using dsfdr/aldex2/dacomp/songbird). Should be .qza. If using ancom/generic tsv, use --p-diff_tsv parameter instead.',
    },
    parameter_descriptions={
        'source': 'Origin of the differentail abundance',
        'method': '"groups" to compare term enrichemnent between significantly enriched sequences in both directions, "correlation" to detect dbbact terms significanly correlated/anti-correlated with the effect size.',
        'random_seed': 'If provided, use as the random seed for the enrichment permutation test (to ensure complete replication)',
        'diff_tsv': ('A tsv table input file (e.g. from ancom when using --p-source ancom, or general tsv when using --p-source tsv).'
                     ' Use instead of --i-diff. When using tsv, file should contain the columns:"id" (sequence), "effect", "pval", "reject".'),
        'ancom_stat': 'the ancom statitical results output file for significant sequence identification (optional - overrides sig-threshold)',
        'maxid': 'The maximal dbBact annotation id to use (to enable replication of results after new annotations are added to dbBact',
        'attack': 'Attack mode',
    },
    output_descriptions={'enriched': 'the enriched dbBact terms'},
    name='dbBact term enrichment for differential abundance results',
    description=("Identify dbBact terms enriched in results of differential abundance (terms significantly more represented in either of the differential abundance groups or correlated with effect size)")
)


plugin.methods.register_function(
    function=embed_seqs,
    inputs={'repseqs': FeatureData[Sequence],
            'table': FeatureTable[Frequency]},
    outputs=[('merged', FeatureTable[Frequency])],
    parameters={},
    input_descriptions={
        'table': 'The biom table with hashed ASV IDs',
        'repseqs': 'The corresponding representative sequences'
    },
    output_descriptions={'merged': 'The biom table with embedded representative sequences'},
    name='embed sequences into hashed biom table',
    description=('embed the representative sequences into the biom table')
)


plugin.methods.register_function(
    function=trim_primers,
    inputs={'repseqs': FeatureData[Sequence],
            'table': FeatureTable[Frequency]},
    outputs=[('trimmed', FeatureTable[Frequency])],
    parameters={},
    input_descriptions={
        'table': 'The biom table with ASVs to be trimmed',
        'repseqs': 'The corresponding representative sequences (if not embedded in the table using the --p-no-hashed-feature-ids parameter when denoising)'
    },
    output_descriptions={'trimmed': 'The biom table with primer-trimmed sequences'},
    name='trim primers from table sequences',
    description=('Make sequences in the biom table compatible with dbBact. The method will automatically identify the dbBact supported primers in the sequences and trim the sequences so '
                 'all sequences in the biom table start at the end of the supported forward primer. Currently supported dbBact primers are V1 (AGAGTTTGATC[AC]TGG[CT]TCAG), '
                 'V3 (CCTACGGG[ACGT][CGT]GC[AT][CG]CAG) and V4 (GTGCCAGC[AC]GCCGCGGTAA)')
)


plugin.methods.register_function(
    function=single_enrichment,
    inputs={'repseqs': FeatureData[Sequence],
            'bg_seqs': FeatureData[Sequence],
            'bg_table': FeatureTable[Frequency],
            'test_seqs': FeatureData[Sequence],
            'test_table': FeatureTable[Frequency]},
    outputs=[('enriched', FeatureData[Differential])],
    parameters={},
    input_descriptions={
        'bg_table': 'The table for background sequences (use this or bg_seqs)',
        'bg_seqs': 'The background sequences (use this or bg_table)',
        'test_table': 'The table with sequences to test for term enrichment (use this or test_seqs)',
        'test_seqs': 'The sequences to test for term enrichment (use this or test_table)',
        'repseqs': 'The corresponding representative sequences'
    },
    output_descriptions={'enriched': 'Table containing the enriched dbBact terms in the test compared to the background table'},
    name='dbBact term enrichment compared to experiment sequence background',
    description=('Test dbBact term enrichment in a test set of sequences compared to an experiment of background sequences')
)


plugin.methods.register_function(
    function=bg_term_enrichment,
    inputs={'repseqs': FeatureData[Sequence],
            'test_seqs': FeatureData[Sequence],
            'test_table': FeatureTable[Frequency]},
    outputs=[('enriched', FeatureData[Differential])],
    parameters={'bg_terms': List[Str],
                'alpha': Float,
                'min_appearance': Int},
    input_descriptions={
        'test_table': 'The table with sequences to test for term enrichment (use this or test_seqs)',
        'test_seqs': 'The sequences to test for term enrichment (use this or test_table)',
        'repseqs': 'The corresponding representative sequences'
    },
    parameter_descriptions={
        'bg_terms': 'List of dbBact terms. Input sequences are compared to sequences from all dbBact experiments containing all the terms. '
                    'For example to compare input sequences to human fecal bacteria from china, use --bg-terms feces china "homo spaiens"',
        'alpha': 'The dsFDR threshold for significantly enriched terms',
        'min_appearance': 'Minimal number of experiments each dbBact background term should appear in (as COMMON)'
    },
    output_descriptions={'enriched': 'Table containing the enriched dbBact terms in the test compared to the background table'},
    name='dbBact term enrichment compared to dbBact sequences associated with a list of terms',
    description=('Test dbBact term enrichment in a test set of sequences compared to all dbBact sequences associated (COMMON) with a set of terms')
)


plugin.methods.register_function(
    function=diff_abundance,
    inputs={'repseqs': FeatureData[Sequence],
            'table': FeatureTable[Frequency]},
    outputs=[('diff', FeatureData[Differential])],
    parameters={
        'metadata': Metadata,
        'field': Str,
        'pair_field': Str,
        'statistical_test': Str % Choices(_statistical_tests),
        'transform_function': Str % Choices(_transform_functions),
        'permutations': Int,
        'alpha': Float,
        'random_seed': Int,
        'fdr_method': Str % Choices(_fdr_method),
        'val1': List[Str],
        'val2': List[Str]
    },
    input_descriptions={
        'table': 'The feature table to perform differential abundance on (needs to be normalized or rarified)',
        'repseqs': 'The corresponding representative sequences (if not embedded in the table using the --p-no-hashed-feature-ids parameter when denoising)'
    },
    parameter_descriptions={
        'alpha': 'The dsFDR threshold for significantly enriched terms',
        'metadata': 'Metadata (mapping) file for the table',
        'field': 'The metadata field on which to perform the differential abundance test',
        'pair_field': 'If supplied, perform paired differential abundance testing pairing samples in the pair_field metadata field',
        'val1': 'Metadata values (in field) for samples in group1. Can supply multiple values',
        'val2': 'Metadata values (in field) for samples in group2. Can supply multiple values. If not provided, take all samples not in group1',
        'statistical_test': 'The statistic used for the effect size calculation. "meandiff" is the difference in the means between the two groups, '
                            '"stdmeandiff" is the difference in the means normalized by the standard deviation.',
        'transform_function': 'The transformation to apply to the frequency data (for each feature) before calculating the effect size. options: '
                              '"rankdata": rank each feature across samples, '
                              '"log2data": log2 transform each frequency, '
                              '"normdata": normalize to constant sum for each feature, '
                              '"binarydata": convert to present/absent, '
                              '"none": do not apply any transformation',
        'permutations': 'Number of permutations to calculate (higher number is more exact but slower)',
        'alpha': 'the FDR threshold for rejecting the null hypothesis',
        'random_seed': 'The random seed to use (for fully reproducible results, since test is permutation based)',
        'fdr_method': 'The multiple hypothesis correction method to use. options: '
                      '"dsfdr": The discrete FDR correction (stronger than BHfdr if there are many low prevalence features) (Jiang et al 2017, https://doi.org/10.1128/mSystems.00092-17), '
                      '"bhfdr": Benhaminy-Hochberg FDR correction, '
                      '"byfdr": Benhaminy-Yekutieli FDR correction, '
                      '"filterBH": Low frequency filtering based correction (see Jiang et al 2017, https://doi.org/10.1128/mSystems.00092-17).'
    },
    output_descriptions={'diff': 'Table of differentially abundant bacteria'},
    name='Calour based differential abundance',
    description=('Identify differentially abundance bacteria between two groups using the calour diff_abundance method')
)


plugin.visualizers.register_function(
    function=draw_wordcloud_vis,
    inputs={'table': FeatureTable[Frequency],
            'repseqs': FeatureData[Sequence],
            },
    parameters={
        'prev_thresh': Float,
        'focus_terms': List[Str],
    },
    input_descriptions={
        'table': 'The biom table to draw the wordcloud for.'
    },
    parameter_descriptions={
        'prev_thresh': 'Mininal prevalence (fraction of samples sequence is present) in order to include sequence in wordcloud stats.',
        'focus_terms': 'show only terms from annotations containing all these terms (comma separated).'
    },
    name='wordcloud',
    description=('draw the dbBact terms wordcloud for ASVs in a given biom table')
)


plugin.visualizers.register_function(
    function=plot_enrichment,
    inputs={'enriched': FeatureData[Differential]},
    parameters={
        'max_show': Int,
        'max_len': Int,
        'labels': List[Str],
        'colors': List[Str],
        'enriched_exp_color': Str
    },
    input_descriptions={
        'enriched': 'The enriched dbbact terms (result of "qiime dbbact enrichment")'
    },
    parameter_descriptions={
        'max_show': 'The maximal number of enriched dbBact terms to show.',
        'max_len': 'Maximal per-term string length (i.e. trim terms after XXX letters for better visibility)',
        'labels': 'Labels of the two groups',
        'colors': 'Color names for group1, group2 bars (valid matplotlib color names)',
        'enriched_exp_color': 'Color name for the enriched experiments number within each bar (if available)'
    },
    name='Visualize dbbact enrichment results (following "enrichment" command)',
    description=('Visualization for qiime dbbact enriched output')
)


plugin.visualizers.register_function(
    function=heatmap,
    inputs={'table': FeatureTable[Frequency],
            'repseqs': FeatureData[Sequence],
            'taxonomy': FeatureData[Taxonomy],
            },
    parameters={'metadata': Metadata,
                'sort_field': Str,
                'cluster': Bool,
                'normalize': Bool,
                'min_abundance': Float,
                },
    input_descriptions={
        'table': 'The table to plot',
        'repseqs': 'If input table contains hashes instead of sequences, supply the representative sequences file used in the table creation.',
        'taxonomy': 'Optional taxonomy file for showing the sequence taxonomy in the heatmap'
    },
    parameter_descriptions={
        'metadata': 'The metadata table (tsv)',
        'sort_field': 'field to sort the samples by before plotting',
        'cluster': 'if true, cluster (reorder) the ASVs before plotting (based on similar behavior across samples',
        'normalize': 'it true, normalize (TSS) the number of reads per sample',
        'min_abundance': 'filter away seqiuences with < min_abundance total reads (over all samples)',
    },
    name='Plot interactive heatmap with dbbact annotations per sequence',
    description=('Plot interactive html heatmap. Samples are columns and sequences are rows.')
)


plugin.visualizers.register_function(
    function=venn,
    inputs={'diff': FeatureData[Differential],
            'repseqs': FeatureData[Sequence],
            },
    parameters={
        'source': Str % Choices(['dsfdr', 'aldex2', 'dacomp', 'songbird', 'ancom', 'tsv']),
        'sig_threshold': Float,
        'diff_tsv': Str,
        'ancom_stat': Str,
        'terms': List[Str],
        'label1': Str,
        'label2': Str,
        'set_colors': List[Str],
        'max_size': Int
    },
    input_descriptions={
        'diff': 'Result of qiime2 differential abundance plugin (using dsfdr/aldex2/dacomp/songbird). Should be .qza. If using ancom/generic tsv, use --p-diff_tsv parameter instead.',
    },
    parameter_descriptions={
        'source': 'Origin of the differentail abundance',
        'diff_tsv': ('A tsv table input file (e.g. from ancom when using --p-source ancom, or general tsv when using --p-source tsv).'
                     ' Use instead of --i-diff. When using tsv, file should contain the columns:"id" (sequence), "effect", "pval", "reject".'),
        'ancom_stat': 'the ancom statitical results output file for significant sequence identification (optional - overrides sig-threshold)',
        'terms': 'List of terms to plot the venn diagram for (blue circle includes sequences appearing in dbBact annotations that include all terms)',
    },
    name='Plot dbBact Venn diagram',
    description='Plot Venn diagram for intersection between two seqs of sequences (from diff. abundance) and dbBact sequences associated with a set of terms',
)
