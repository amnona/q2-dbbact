import qiime2.plugin
from qiime2.plugin import (SemanticType, Str, Int, Float, Choices,
                           MetadataColumn, Categorical, Numeric, Metadata)
import qiime2.plugin.model as model

from q2_types.feature_table import (
    FeatureTable, Frequency)
from q2_types.feature_data import (FeatureData, Differential, Sequence)
from q2_types.sample_data import SampleData
from q2_types.feature_data import DNAFASTAFormat
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


def enrichment(data: pd.DataFrame, repseqs: DNAFASTAFormat = None, source: str = 'dsfdr', method: str = 'groups', sig_threshold: float = 0.1) -> pd.DataFrame:
    ca.set_log_level('INFO')
    db = dbbact_calour.DBBact()

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


def draw_wordcloud_vis(output_dir: str, data: biom.Table, repseqs: DNAFASTAFormat = None, prev_thresh: float = 0.3):
    '''draw the wordcloud for features in the biom table and save outputs'''
    db = dbbact_calour.DBBact()

    df = data.to_dataframe()

    if repseqs is not None:
        data = _seqs_from_repseqs(data=df, repseqs=repseqs)

    if len(df.iloc[0].name) == 32:
        raise ValueError('input table seems to contain hashes and not sequences. Please supply the rep-seqs.qza file.')
    exp = ca.AmpliconExperiment.from_pandas(df.T)
    exp = exp.filter_prevalence(prev_thresh)
    f = db.draw_wordcloud(exp)
    f.savefig(os.path.join(output_dir, 'wordcloud.pdf'))
    f.savefig(os.path.join(output_dir, 'wordcloud.svg'))
    with open(os.path.join(output_dir, 'index.html'), 'w') as fl:
        fl.write('<html><body>\n')
        fl.write('<h1>dbBact wordcloud</h1>\n')
        fl.write('<h2>prevalence threshold=%f</h2>\n' % prev_thresh)
        fl.write('<img src="wordcloud.svg" alt="wordcloud">\n')
        fl.write('<a href="wordcloud.pdf">Download as PDF</a><br>\n')
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


def heatmap(output_dir: str, table: biom.Table, metadata: pd.DataFrame):
    metadata = metadata.to_dataframe()
    samples = table.ids(axis='sample')
    features = table.ids(axis='observation')
    metadata = metadata.filter(items=samples, axis='index')
    feature_metadata = pd.DataFrame(index=features, columns={'_feature_id': features, 'taxonomy': 'NA'})
    exp = ca.AmpliconExperiment(data=table.transpose().matrix_data, sample_metadata=metadata, feature_metadata=feature_metadata, sparse=False)
    exp.export_html(output_file=os.path.join(output_dir, 'index.html'))
    # exp.export_html(output_file=os.path.join(output_dir, 'index.html'))


_source_types = ['dsfdr', 'aldex2', 'dacomp', 'songbird']


plugin.methods.register_function(
    function=enrichment,
    inputs={'data': FeatureData[Differential],
            'repseqs': FeatureData[Sequence],
            },
    outputs=[('enriched', FeatureData[Differential])],
    parameters={
        'source': Str % Choices(_source_types),
        'method': Str % Choices('groups', 'correlation'),
        'sig_threshold': Float,
    },
    input_descriptions={
        'data': ('The feature table containing the samples over which beta '
                 'diversity should be computed.')
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
        'prev_thresh': Float
    },
    input_descriptions={
        'data': 'The biom table to draw the wordcloud for.'
    },
    parameter_descriptions={},
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
    inputs={'table': FeatureTable[Frequency]},
    parameters={'metadata': Metadata},
    input_descriptions={
        'table': 'The table to plot'
    },
    parameter_descriptions={},
    name='heat',
    description=('Plot interactive heatmap')
)
