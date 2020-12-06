import qiime2.plugin
from qiime2.plugin import (SemanticType, Str, Int, Float, Choices,
                           MetadataColumn, Categorical, Numeric, Metadata, Bool)
import qiime2.plugin.model as model

from q2_types.feature_table import (
    FeatureTable, Frequency)
from q2_types.feature_data import (FeatureData, Differential, Sequence, Taxonomy, )
from q2_types.sample_data import SampleData
from q2_types.feature_data import DNAFASTAFormat, TSVTaxonomyFormat
import biom
import pandas as pd
import numpy as np
import os
import calour as ca
import dbbact_calour

from . import __version__

_citation = ('manuscript in preparation')

_short_description = "Plugin for multiple comparisons in sparse Microbiome Data"

plugin = qiime2.plugin.Plugin(
    name='dbbact',
    version=__version__,
    website='https://github.com/amnona/q2_dbbact',
    package='q2_dbbact',
    short_description='dbBact plugin for Qiime2',
    description=('A Qiime2 plugin for the dbBact bacterial knowledge-base (dbbact.org)'),
    citation_text=_citation
)


def _iter_fasta(fp):
    '''Read fasta file into iterator.

    Fasta file must contain header line (starting with ">") and one or more sequence lines.

    Parameters
    ----------
    fp: str
        name of the fasta file

    Yields
    ------
    header: str
        the header line (without ">")
    sequence: str
        the sequence ('ACGT'). Both header and sequence are whitespace stripped.
    '''
    # skip non-header lines at beginning of file
    with open(fp, 'r') as fl:
        for cline in fl:
            if cline[0] == ">":
                title = cline[1:].rstrip()
                break
            print('Fasta file %s has no headers' % fp)
            return

        lines = []
        for cline in fl:
            if cline[0] == ">":
                yield title, ''.join(lines)
                lines = []
                title = cline[1:].strip()
                continue
            lines.append(cline.strip())
        yield title, "".join(lines)


def embed_seqs(table: biom.Table, repseqs: DNAFASTAFormat) -> biom.Table:
    ids = set(table.ids(axis='observation'))
    id_map = {}
    for chead, cseq in _iter_fasta(str(repseqs)):
        if chead in ids:
            id_map[chead] = cseq
    table.update_ids(id_map, axis='observation', strict=False, inplace=True)
    return table


def _seqs_from_repseqs(data: pd.DataFrame, repseqs: DNAFASTAFormat):
    '''add sequences from repseqs to the data
    '''
    # check if feature ids are hashes by looking at first feature
    if len(data.iloc[0].name) == 32:
        print('replacing hashes with sequences')
        rename_ids = {}
        for chead, cseq in _iter_fasta(str(repseqs)):
            if chead in data.index:
                rename_ids[chead] = cseq
        data.rename(index=rename_ids, inplace=True)
    else:
        print('table does not seem to contain qiime2 hashes for feature ids.\nFirst feature id is %s\n' % data.iloc[0].name)
    return data


def enrichment(diff: pd.DataFrame = None, repseqs: DNAFASTAFormat = None, diff_tsv: str = None, source: str = 'dsfdr', method: str = 'groups', sig_threshold: float = 0.1, ancom_stat: str = None) -> pd.DataFrame:
    ca.set_log_level('INFO')
    db = dbbact_calour.DBBact()
    db.set_log_level('INFO')

    # handle either diff (qza diffabundance artifact) or diff_tsv(tsv file resulting from diff abundance export i.e. ancom etc.)
    if diff is None and diff_tsv is None:
        raise ValueError('Need to supply either --i-diff or --p-diff-tsv differential results input file (none were supplied)')
    if diff is not None and diff_tsv is not None:
        raise ValueError('Need to supply either --i-diff or --p-diff-tsv differential results input file (both were supplied)')

    if diff is not None:
        data = diff
    else:
        try:
            data = pd.read_csv(diff_tsv, sep='\t')
            data.set_index(data.columns[0], drop=False, inplace=True)
        except Exception as e:
            print('failed to read diff_tsv file %s. Is is a valid tsv differential abundance file?' % diff_tsv)
            raise e

    if repseqs is not None:
        data = _seqs_from_repseqs(data=data, repseqs=repseqs)

    ndata = data.copy()
    field_name = None

    # check if feature ids are hashes by looking at first feature
    if len(data.iloc[0].name) == 32:
        print('Features seem to be hashes. First id=%s' % data.iloc[0].name)
        raise ValueError('Input file contains sequence hashes instead of actual sequences.\nPlease use qiime dbbact enrichment-hash and supply the rep_seqs.qza.')

    ndata['_feature_id'] = ndata.index.values
    if source == 'dsfdr':
        ndata = pd.DataFrame(data={'dir': data['Statistic'] > 0,
                                   'pval': data['raw pvalue'],
                                   'effect': data['Statistic'],
                                   'reject': data['Reject']}, index=ndata.index)
    elif source == 'songbird':
        vals = list(ndata.columns.values)
        if 'Intercept' in vals:
            vals.remove('Intercept')
        else:
            print('"Intercept" not in input data. Colums are: %s. Is it a songbird differentials file?' % vals)
            raise ValueError('File does not seem to be a songbird differentials file (no "Intercept" column)')
        vals.remove('_feature_id')
        if len(vals) > 1:
            print('More than 1 column detected: %s' % vals)
            raise ValueError('More than 1 column in songbird results (%s). Current version works only on 1 column' % vals)
        field_name = vals[0]
        # we don't have p-values so do not reject any
        ndata = pd.DataFrame(data={'effect': ndata[vals[0]], 'reject': 1}, index=ndata.index)

    elif source == 'aldex2':
        ndata = pd.DataFrame(data={'dir': data['effect'] > 0,
                                   'pval': data['we.eBH'],
                                   'effect': data['effect'],
                                   'reject': data['we.eBH'] < sig_threshold}, index=ndata.index)
        ndata = ndata[ndata['reject'] == 1]

    elif source == 'ancom':
        # if user supplied the ancom stats file, use it to get significant (null hypothesis rejects)
        # otherwise, use the sig_threshold for the W stat
        ndata = pd.DataFrame(data={'dir': data['clr'] > 0,
                                   'effect': data['clr'],
                                   'pval': 0,
                                   'reject': data['W'] > sig_threshold}, index=ndata['id'])
        if ancom_stat is not None:
            try:
                rejects = pd.read_csv(ancom_stat, sep='\t', index_col=0)
                ndata = ndata.join(rejects)
                ndata['reject'] = ndata['Reject null hypothesis'] * 1
            except Exception as e:
                print('failed to read ancom-stat file %s. Is is a valid tsv differential abundance file?' % ancom_stat)
                raise e

    elif source == 'tsv':
        if 'pval' not in data.columns:
            raise ValueError('input tsv file does not contain "pval" column')
        if 'effect' not in data.columns:
            raise ValueError('input tsv file does not contain "effect" column')
        if 'reject' not in data.columns:
            data['reject'] = 1
        ndata = pd.DataFrame(data={'dir': data['effect'] > 0,
                                   'effect': data['effect'],
                                   'pval': data['pval'],
                                   'reject': data['reject']}, index=ndata['id'])

    else:
        raise ValueError('Unsupported source %s' % source)

    exp = ca.AmpliconExperiment(data=np.zeros([1, len(ndata)]),
                                sample_metadata=pd.DataFrame({'_sample_id': ['s1']}),
                                feature_metadata=ndata,
                                sparse=False)

    if method == 'groups':
        print('%d features' % len(exp.feature_metadata))
        exp = exp.filter_by_metadata('reject', ['1'], axis='f')
        print('%d significant features' % len(exp.feature_metadata))
        pos_features = ndata[ndata['effect'] > 0].index.values
        res = db.enrichment(exp, features=pos_features)
        df = res[0]
    elif method == 'correlation':
        print('searching for dbbact term enrichment for %d features' % len(exp.feature_metadata))
        dd = db.rank_enrichment(exp, 'effect')
        print('found %d enriched terms' % len(dd.feature_metadata))
        df = dd.feature_metadata
        df.rename({'_calour_stat': 'odif'}, axis='columns', inplace=True)
        df.drop('_calour_direction', axis='columns', inplace=True)
    else:
        raise ValueError('Unsupported method %s' % source)

    # qiime2 plugin stuff:
    # need to rename the index to "id"
    # and cannot contain non-numeric columns
    df.drop('term', axis='columns', inplace=True)
    df.index.rename('id', inplace=True)
    return df


def draw_wordcloud_vis(output_dir: str, data: biom.Table, repseqs: DNAFASTAFormat = None, prev_thresh: float = 0.3, focus_terms: str = None):
    '''draw the wordcloud for features in the biom table and save outputs'''
    db = dbbact_calour.DBBact()
    db.set_log_level('INFO')

    df = data.to_dataframe()

    if repseqs is not None:
        data = _seqs_from_repseqs(data=df, repseqs=repseqs)

    if len(df.iloc[0].name) == 32:
        raise ValueError('input table seems to contain hashes and not sequences. Please supply the rep-seqs.qza file.')

    exp = ca.AmpliconExperiment.from_pandas(df.T)
    print('%d features in table before prevalence filtering' % len(exp.feature_metadata))
    exp = exp.filter_prevalence(prev_thresh)
    print('%d features in table remain after prevalence filtering (%f)' % (len(exp.feature_metadata), prev_thresh))
    if focus_terms is not None:
        focus_terms = focus_terms.split(',')
        print('Using %d focus terms: %s' % (len(focus_terms), focus_terms))
    print('getting term stats')
    fscores, recall, precision, term_count, reduced_f = db.get_wordcloud_stats(exp=exp, focus_terms=focus_terms)
    df = pd.DataFrame(data={'fscore': fscores, 'recall': recall, 'precision': precision, 'term_count': term_count, 'reduced_f': reduced_f})
    df.to_csv(os.path.join(output_dir, 'scores.tsv'), sep='\t')
    print('drawing wordcloud')
    f = db.draw_wordcloud(exp, focus_terms=focus_terms)
    f.savefig(os.path.join(output_dir, 'wordcloud.pdf'))
    f.savefig(os.path.join(output_dir, 'wordcloud.svg'))
    with open(os.path.join(output_dir, 'index.html'), 'w') as fl:
        fl.write('<html><body>\n')
        fl.write('<h1>dbBact wordcloud</h1>\n')
        fl.write('<h2>prevalence threshold=%f</h2>\n' % prev_thresh)
        if focus_terms is not None:
            fl.write('<h2>focus-terms: %s</h2>' % focus_terms)
        fl.write('<br><img src="wordcloud.svg" alt="wordcloud">\n')
        fl.write('<a href="wordcloud.pdf">Download as PDF</a><br>\n')
        fl.write('<a href="scores.tsv">Download scores as tsv</a><br>\n')
        fl.write('</body></html>')


def plot_enrichment(output_dir: str, enriched: pd.DataFrame):
    enriched['term'] = enriched.index.values
    ax = ca.plotting.plot_enrichment(None, enriched)
    ax.figure.tight_layout()
    ax.figure.savefig(os.path.join(output_dir, 'enriched_terms.svg'), bbox_inches='tight')
    ax.figure.savefig(os.path.join(output_dir, 'enriched_terms.pdf'), bbox_inches='tight')
    with open(os.path.join(output_dir, 'index.html'), 'w') as fl:
        fl.write('<html><body>\n')
        fl.write('<h1>dbBact enriched terms</h1>\n')
        fl.write('<img src="enriched_terms.svg" alt="enriched terms barplot">\n')
        fl.write('<a href="enriched_terms.pdf">Download as PDF</a><br>\n')
        fl.write('</body></html>')


def heatmap(output_dir: str, table: biom.Table, metadata: pd.DataFrame, sort_field: str = None, cluster: bool = True, min_abundance: float = 10, normalize: bool = True,
            taxonomy: pd.DataFrame = None, repseqs: DNAFASTAFormat = None):
    ca.set_log_level('INFO')
    print('creating experiment')
    metadata = metadata.to_dataframe()
    samples = table.ids(axis='sample')
    features = table.ids(axis='observation')
    # combine the taxonomy
    feature_metadata = pd.DataFrame(index=features, columns={'_feature_id': features})
    if taxonomy is not None:
        print('adding taxonomy')
        feature_metadata['taxonomy'] = taxonomy['Taxon']
    else:
        feature_metadata['taxonomy'] = 'NA'
    # use repseqs if supplied
    if repseqs is not None:
        print('converting hashes to sequences using repseqs file')
        feature_metadata = _seqs_from_repseqs(feature_metadata, repseqs)

    metadata = metadata.filter(items=samples, axis='index')
    exp = ca.AmpliconExperiment(data=table.transpose().matrix_data, sample_metadata=metadata, feature_metadata=feature_metadata, sparse=False)
    print('created experiment %r' % exp)
    if normalize:
        print('normalizing to 10000 reads/sample')
        exp = exp.normalize(10000)
    if min_abundance > 0:
        print('removing sequences with sum abundance < %f' % min_abundance)
        exp = exp.filter_sum_abundance(min_abundance)
        print('remaining experiment %r' % exp)
    if sort_field is not None:
        print('sorting by field %s' % sort_field)
        exp = exp.sort_samples(sort_field)
    if cluster:
        print('clustering sequences')
        exp = exp.cluster_features(10)
    print('generating heatmap to output file %s' % output_dir)
    exp.export_html(output_file=os.path.join(output_dir, 'index.html'), sample_field=sort_field, feature_field='taxonomy')


plugin.methods.register_function(
    function=enrichment,
    inputs={'diff': FeatureData[Differential],
            'repseqs': FeatureData[Sequence],
            },
    outputs=[('enriched', FeatureData[Differential])],
    parameters={
        'source': Str % Choices(['dsfdr', 'aldex2', 'dacomp', 'songbird', 'ancom', 'tsv']),
        'method': Str % Choices('groups', 'correlation'),
        'sig_threshold': Float,
        'diff_tsv': Str,
        'ancom_stat': Str,
    },
    input_descriptions={
        'diff': 'Result of qiime2 differential abundance plugin (using dsfdr/aldex2/dacomp/songbird). Should be .qza. If using ancom/generic tsv, use --p-diff_tsv parameter instead.',
    },
    parameter_descriptions={
        'source': 'Origin of the differentail abundance',
        'method': '"groups" to compare term enrichemnent between significantly enriched sequences in both directions, "correlation" to detect dbbact terms significanly correlated/anti-correlated with the effect size.',
        'diff_tsv': ('A tsv table input file (e.g. from ancom when using --p-source ancom, or general tsv when using --p-source tsv).'
                     ' Use instead of --i-diff. When using tsv, file should contain the columns:"id" (sequence), "effect", "pval", "reject".'),
        'ancom_stat': 'the ancom statitical results output file for significant sequence identification (optional - overrides sig-threshold)',
    },
    output_descriptions={'enriched': 'the enriched features'},
    name='enrichemnt',
    description=("dbBact term enrichment")
)


plugin.methods.register_function(
    function=embed_seqs,
    inputs={'repseqs': FeatureData[Sequence],
            'table': FeatureTable[Frequency]},
    outputs=[('merged', FeatureTable[Frequency])],
    parameters={},
    input_descriptions={
        'table': 'The feature table',
        'repseqs': 'The corresponding representative sequences'
    },
    output_descriptions={'merged': 'The feature table with embedded representative sequences'},
    name='embed sequences into feature table',
    description=('embed the representative sequences into the feature table')
)


plugin.visualizers.register_function(
    function=draw_wordcloud_vis,
    inputs={'data': FeatureTable[Frequency],
            'repseqs': FeatureData[Sequence],
            },
    parameters={
        'prev_thresh': Float,
        'focus_terms': Str,
    },
    input_descriptions={
        'data': 'The biom table to draw the wordcloud for.'
    },
    parameter_descriptions={
        'prev_thresh': 'Mininal prevalence (fraction of samples sequence is present) in order to include sequence in wordcloud stats.',
        'focus_terms': 'show only terms from annotations containing all these terms (comma separated).'
    },
    name='wordcloud',
    description=('draw wordcloud')
)


plugin.visualizers.register_function(
    function=plot_enrichment,
    inputs={'enriched': FeatureData[Differential]},
    parameters={
    },
    input_descriptions={
        'enriched': 'The enriched dbbact terms (result of qiime dbbact enrichment)'
    },
    parameter_descriptions={},
    name='plot_enrichment',
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
        'table': 'The table to plot'
    },
    parameter_descriptions={
        'metadata': 'The metadata table (tsv)',
        'sort_field': 'field to sort the samples by before plotting',
        'cluster': 'if true, cluster (reorder) the sequences before plotting (based on similar behavior across samples',
        'normalize': 'it true, normalize (TSS) the number of reads per sample',
        'min_abundance': 'filter away seqiuences with < min_abundance total reads (over all samples)',
    },
    name='heat',
    description=('Plot interactive heatmap')
)
