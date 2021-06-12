from q2_types.feature_data import DNAFASTAFormat
import biom
import pandas as pd
import numpy as np
import os
import calour as ca
import dbbact_calour
import dbbact_calour.dbbact
from typing import List

from .utils import _seqs_from_repseqs, _load_diff_abundance, _test_exact_region


def draw_wordcloud_vis(output_dir: str, table: biom.Table, repseqs: DNAFASTAFormat = None, prev_thresh: float = 0.3, focus_terms: List[str] = None):
    '''draw the wordcloud for features in the biom table and save outputs'''
    db = dbbact_calour.dbbact.DBBact()
    db.set_log_level('INFO')

    df = table.to_dataframe()

    if repseqs is not None:
        df = _seqs_from_repseqs(data=df, repseqs=repseqs)

    if len(df.iloc[0].name) == 32:
        raise ValueError('input table seems to contain hashes and not sequences. Please supply the rep-seqs.qza file.')

    exp = ca.AmpliconExperiment.from_pandas(df.T)

    # check if the experiment contains trimmed sequences. otherwise let the user know
    seqs = exp.feature_metadata.sample(n=np.min([500, len(exp.feature_metadata)])).index.values
    region = _test_exact_region(seqs)
    if region is None:
        raise ValueError('Table seems to contain untrimmed sequences. Please run qiime dbbact trim-primers on the table prior to running enrichment.')
    else:
        print('Identified region %s' % region)

    print('%d ASVs in table before prevalence filtering' % len(exp.feature_metadata))
    exp = exp.filter_prevalence(prev_thresh)
    print('%d ASVs in table remain after prevalence filtering (%f)' % (len(exp.feature_metadata), prev_thresh))
    if focus_terms is not None:
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


def plot_enrichment(output_dir: str, enriched: pd.DataFrame, max_show: int = 10, max_len: int = 40, colors: List[str] = ['green', 'red'],
                    labels: List[str] = ['group1', 'group2'], enriched_exp_color: str = 'white'):
    '''Plot the enriched terms bar plot
    '''
    enriched['term'] = enriched.index.values
    ax = ca.plotting.plot_enrichment(None, enriched, max_show=max_show, max_len=max_len, colors=colors, labels=labels, enriched_exp_color=enriched_exp_color)
    ax.figure.tight_layout()
    ax.figure.savefig(os.path.join(output_dir, 'enriched_terms.svg'), bbox_inches='tight')
    ax.figure.savefig(os.path.join(output_dir, 'enriched_terms.pdf'), bbox_inches='tight')
    enriched.to_csv(os.path.join(output_dir, 'enriched_terms.tsv'), sep='\t')
    with open(os.path.join(output_dir, 'index.html'), 'w') as fl:
        fl.write('<html><body>\n')
        fl.write('<h1>dbBact enriched terms</h1>\n')
        fl.write('<img src="enriched_terms.svg" alt="enriched terms barplot"><br><br>\n')
        fl.write('<a href="enriched_terms.pdf">Download as PDF</a><br>\n')
        fl.write('<a href="enriched_terms.tsv">Download enriched terms as tsv</a><br>\n')
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


def venn(output_dir: str, terms: List[str], diff: pd.DataFrame = None, repseqs: DNAFASTAFormat = None, diff_tsv: str = None, source: str = 'dsfdr', sig_threshold: float = 0.1, ancom_stat: str = None,
         label1: str = 'Negative', label2: str = 'Positive', set_colors: List[str] = ['red', 'green', 'mediumblue'], max_size: int = 1000):
    ca.set_log_level('INFO')
    db = dbbact_calour.dbbact.DBBact()
    db.set_log_level('INFO')

    # load the differetial abundance table into a calour experiment
    exp = _load_diff_abundance(diff=diff, diff_tsv=diff_tsv, source=source, ancom_stat=ancom_stat, repseqs=repseqs, sig_threshold=sig_threshold)
    # exp.feature_metadata['_calour_direction'] = exp.feature_metadata['dir']
    exp.feature_metadata['_calour_direction'] = exp.feature_metadata['dir'].replace({0: label1, 1: label2}, inplace=False)

    # check if the experiment contains trimmed sequences. otherwise let the user know
    seqs = exp.feature_metadata.sample(n=np.min([500, len(exp.feature_metadata)])).index.values
    region = _test_exact_region(seqs)
    if region is None:
        raise ValueError('Table seems to contain untrimmed sequences. Please run qiime dbbact trim-primers on the table prior to running enrichment.')
    else:
        print('Identified region %s' % region)

    f = db.plot_term_venn_all(terms, exp, max_size=max_size, set_colors=set_colors)
    f.savefig(os.path.join(output_dir, 'venn.svg'))
    f.savefig(os.path.join(output_dir, 'venn.pdf'))
    with open(os.path.join(output_dir, 'index.html'), 'w') as fl:
        fl.write('<html><body>\n')
        fl.write('<h1>dbBact venn for terms %s</h1>\n' % ','.join(terms))
        fl.write('<img src="venn.svg" alt="term venn diagram"><br><br>\n')
        fl.write('<a href="venn.pdf">Download as PDF</a><br>\n')
        fl.write('</body></html>')
