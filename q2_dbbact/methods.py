from q2_types.feature_data import DNAFASTAFormat
import biom
import pandas as pd
import numpy as np
import calour as ca
import dbbact_calour
import dbbact_calour.dbbact
from typing import List

from .utils import _iter_fasta, _seqs_from_repseqs, _load_diff_abundance, _test_exact_region, _test_embedded_primers, _trim_known_primer
from .visualizations import venn, heatmap, plot_enrichment, draw_wordcloud_vis


def embed_seqs(table: biom.Table, repseqs: DNAFASTAFormat) -> biom.Table:
    ids = set(table.ids(axis='observation'))
    id_map = {}
    for chead, cseq in _iter_fasta(str(repseqs)):
        if chead in ids:
            id_map[chead] = cseq
    table.update_ids(id_map, axis='observation', strict=False, inplace=True)
    return table


def enrichment(diff: pd.DataFrame = None, repseqs: DNAFASTAFormat = None, diff_tsv: str = None, source: str = 'dsfdr', method: str = 'groups', sig_threshold: float = 0.1, ancom_stat: str = None, attack: bool = False, maxid: int = None, random_seed: int = None) -> pd.DataFrame:
    ca.set_log_level('INFO')
    db = dbbact_calour.dbbact.DBBact()
    db.set_log_level('INFO')

    # load the differetial abundance table into a calour experiment
    exp = _load_diff_abundance(diff=diff, diff_tsv=diff_tsv, source=source, ancom_stat=ancom_stat, repseqs=repseqs, sig_threshold=sig_threshold)
    ndata = exp.feature_metadata

    # check if the experiment contains trimmed sequences. otherwise let the user know
    seqs = ndata.sample(n=np.min([500, len(ndata)])).index.values
    region = _test_exact_region(seqs)
    if region is None:
        raise ValueError('Table seems to contain untrimmed sequences. Please run qiime dbbact trim-primers on the table prior to running enrichment.')
    else:
        print('Identified region %s' % region)

    # do the dbbact term enrichment test (using 2 feature groups or feature effect size rank correlation)
    if method == 'groups':
        print('%d ASVs' % len(exp.feature_metadata))
        exp = exp.filter_by_metadata('reject', ['1'], axis='f')
        print('%d significant ASVs' % len(exp.feature_metadata))
        pos_features = ndata[ndata['effect'] > 0].index.values
        res = db.enrichment(exp, features=pos_features, max_id=maxid, random_seed=random_seed)
        df = res[0]
    elif method == 'correlation':
        print('searching for dbbact term enrichment for %d ASVs' % len(exp.feature_metadata))
        dd = db.rank_enrichment(exp, 'effect', max_id=maxid, random_seed=random_seed)
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


def single_enrichment(bg_table: biom.Table = None, test_table: biom.Table = None, bg_seqs: DNAFASTAFormat = None, test_seqs: DNAFASTAFormat = None,
                      repseqs: DNAFASTAFormat = None) -> pd.DataFrame:
    # test we have the correct inputs
    if bg_table is None and bg_seqs is None:
        raise ValueError('Must supply either bg_table or bg_seqs')
    if test_table is None and test_seqs is None:
        raise ValueError('Must supply either test_table or test_seqs')
    if bg_table is not None and bg_seqs is not None:
        raise ValueError('Cannot use both bg_table and bg_seqs')
    if test_table is not None and test_seqs is not None:
        raise ValueError('Cannot use both test_table and test_seqs')

    # load background sequences to bg set
    if bg_table is not None:
        print('using background sequences from biom table')
        bg = set(bg_table.ids(axis='observation'))
    else:
        print('using backgound sequences from fasta file')
        bg = set()
        for chead, cseq in _iter_fasta(str(bg_seqs)):
            bg.add(cseq)
    print('loaded %d background sequences' % len(bg))

    # load test sequences to seqs set
    if test_table is not None:
        print('using test sequences from biom table')
        seqs = set(test_table.ids(axis='observation'))
    else:
        print('using test sequences from fasta file')
        seqs = set()
        for chead, cseq in _iter_fasta(str(test_seqs)):
            seqs.add(cseq)
    print('loaded %d test sequences' % len(seqs))

    db = dbbact_calour.dbbact.DBBact()
    db.set_log_level('INFO')

    # create a mock experiment for carrying the features
    all_features = bg.union(seqs)
    exp = ca.AmpliconExperiment(data=np.ones([1, len(all_features)]),
                                sample_metadata=pd.DataFrame({'_sample_id': ['s1']}),
                                feature_metadata=pd.DataFrame({'_feature_id': all_features}, index=all_features),
                                sparse=False)

    # check if the experiment contains trimmed sequences. otherwise let the user know
    seqs = exp.feature_metadata.sample(n=np.min([500, len(exp.feature_metadata)])).index.values
    region = _test_exact_region(seqs)
    if region is None:
        raise ValueError('Table seems to contain untrimmed sequences. Please run qiime dbbact trim-primers on the table prior to running enrichment.')
    else:
        print('Identified region %s' % region)

    res = db.enrichment(exp=exp, features=seqs, method='card_mean')
    term_table = res[0]
    if len(term_table) == 0:
        print('No significantly enriched terms found')
    else:
        print('%d Significantly enriched terms found' % len(term_table))
    return term_table


def bg_term_enrichment(bg_terms: List[str], test_table: biom.Table = None, test_seqs: DNAFASTAFormat = None,
                       repseqs: DNAFASTAFormat = None, alpha: float = 0.1, min_appearance: int = 2) -> pd.DataFrame:
    print('term_bg_enrichment for %d terms: %s' % (len(bg_terms), bg_terms))

    # test we have the correct inputs
    if test_table is None and test_seqs is None:
        raise ValueError('Must supply either test_table or test_seqs')
    if test_table is not None and test_seqs is not None:
        raise ValueError('Cannot use both test_table and test_seqs')

    # load test sequences to seqs set
    if test_table is not None:
        print('using test sequences from biom table')
        seqs = set(test_table.ids(axis='observation'))
    else:
        print('using test sequences from fasta file')
        seqs = set()
        for chead, cseq in _iter_fasta(str(test_seqs)):
            seqs.add(cseq)
    print('loaded %d test sequences' % len(seqs))

    db = dbbact_calour.dbbact.DBBact()
    db.set_log_level('INFO')

    # prepare the mock AmpliconExperiment for the background enrichment.
    # We only care about the features (sequencesw)
    exp = ca.AmpliconExperiment(data=np.ones([1, len(seqs)]),
                                sample_metadata=pd.DataFrame({'_sample_id': ['s1']}),
                                feature_metadata=pd.DataFrame({'_feature_id': seqs}, index=seqs),
                                sparse=False)

    # check if the experiment contains trimmed sequences. otherwise let the user know
    seqs = exp.feature_metadata.sample(n=np.min([500, len(exp.feature_metadata)])).index.values
    region = _test_exact_region(seqs)
    if region is None:
        raise ValueError('Table seems to contain untrimmed sequences. Please run qiime dbbact trim-primers on the table prior to running enrichment.')
    else:
        print('Identified region %s' % region)

    exp = exp.normalize()

    # do the background enrichment analysis
    res = db.background_enrich(bg_terms, exp, ignore_exp=None, min_appearance=min_appearance, include_shared=True, alpha=alpha)
    df = res[0]
    df.drop('term', axis='columns', inplace=True)
    df.index.rename('id', inplace=True)

    return res[0]


def diff_abundance(table: biom.Table,
                   metadata: pd.DataFrame,
                   field: str,
                   pair_field: str = None,
                   # metadata2: MetadataColumn = None,
                   repseqs: DNAFASTAFormat = None,
                   statistical_test: str = 'meandiff',
                   transform_function: str = 'rankdata',
                   alpha: float = 0.1,
                   permutations: int = 1000,
                   random_seed: int = 0,
                   fdr_method: str = 'dsfdr',
                   val1: List[str] = None,
                   val2: List[str] = None) -> pd.DataFrame:

    # allow debug info. q2 takes care of what to show using the --verbose flag
    ca.set_log_level('INFO')

    # set random seed to None (i.e. random every time) if 0 is provided
    if random_seed == 0:
        random_seed = None

    metadata = metadata.to_dataframe()

    # take care of the 'none' option for the transfrom
    if transform_function == 'none':
        transform_function = None

    # check and prepare the metadata columns for the testing
    md_fields = metadata.columns
    if field not in md_fields:
        raise ValueError('Field %s not found in metadata file' % field)
    if pair_field is not None:
        if pair_field not in md_fields:
            raise ValueError('Pair field %s not found in metadata file' % field)
    labels = metadata[field].unique()
    n = len(labels)
    if n == 1:
        raise ValueError('Only one category in metadata column. Aborting')
    if val1 is None and val2 is None:
        if n == 2:
            val1 = [labels[0]]
            val2 = [labels[1]]
            print('No values supplied for field, assigning from the two values present (%s, %s)' % (val1, val2))
        else:
            raise ValueError('Cannot perform %s test on more than two categories in metadata and val1 and val2 not supplied. Aborting' % statistical_test)

    if val1 is not None:
        if len([x for x in metadata[field] if x in val1]) == 0:
            raise ValueError('No values matching val1 (%s) found in metadata field %s' % (val1, field))
    if val2 is not None:
        if len([x for x in metadata[field] if x in val2]) == 0:
            raise ValueError('No values matching val2 (%s) found in metadata field %s' % (val2, field))

    # create the calour.AmpliconExperiment from the table and metadata
    samples = table.ids(axis='sample')
    features = table.ids(axis='observation')
    feature_metadata = pd.DataFrame(index=features, columns={'_feature_id': features})
    # use repseqs if supplied
    if repseqs is not None:
        print('converting hashes to sequences using repseqs file')
        feature_metadata = _seqs_from_repseqs(feature_metadata, repseqs)
    metadata = metadata.filter(items=samples, axis='index')
    exp = ca.AmpliconExperiment(data=table.transpose().matrix_data, sample_metadata=metadata, feature_metadata=feature_metadata, sparse=False)
    print('created experiment %r' % exp)

    if pair_field is None:
        print('performing differential abundance testing')
        dd = exp.diff_abundance(field, val1, val2, alpha=alpha, random_seed=random_seed, numperm=permutations, transform=transform_function, method=statistical_test, fdr_method=fdr_method)
    else:
        print('performing paired differential abundance testing')
        dd = exp.diff_abundance_paired(pair_field, field, val1, val2, alpha=alpha, random_seed=random_seed, numperm=permutations, transform=transform_function, method=statistical_test, fdr_method=fdr_method)

    # if no significant results found, need to raise value error since qiime2 does not support empty diff abundance results
    if len(dd.feature_metadata) == 0:
        raise ValueError('No significant results found')

    res = pd.DataFrame({'Reject': 1,
                        'Statistic': dd.feature_metadata['_calour_stat'],
                        'raw pvalue': dd.feature_metadata['_calour_pval'],
                        'dir': dd.feature_metadata['_calour_stat'] > 0,
                        'qvalue': dd.feature_metadata['_calour_qval']},
                       index=dd.feature_metadata.index)
    res['dir'] = res['dir'].astype(int)
    res.index.rename('featureid', inplace=True)
    return res


def enrich_pipeline(ctx,
                    table,
                    # output_dir: str,
                    metadata,
                    field,
                    pair_field=None,
                    repseqs=None,
                    statistical_test='meandiff',
                    transform_function='rankdata',
                    permutations=1000,
                    fdr_method='dsfdr',
                    val1=None,
                    val2=None,
                    method='groups',
                    sig_threshold=0.1,
                    attack=False,
                    maxid=None,
                    random_seed=0):
    res = []

    print('trimming primers if needed')
    trim_func = ctx.get_action('dbbact', 'trim_primers')
    table, = trim_func(table=table, repseqs=repseqs)
    res.append(table)

    print('generating initial wordcloud')
    wordcloud_func = ctx.get_action('dbbact', 'draw_wordcloud_vis')
    wordcloud, = wordcloud_func(table=table, repseqs=repseqs)
    res.append(wordcloud)

    print('detecting differetially abundant features')
    diff_func = ctx.get_action('dbbact', 'diff_abundance')
    diff_table, = diff_func(table=table, metadata=metadata, field=field, pair_field=pair_field, repseqs=repseqs, statistical_test=statistical_test,
                            transform_function=transform_function, alpha=sig_threshold, permutations=permutations, random_seed=random_seed, fdr_method=fdr_method,
                            val1=val1, val2=val2)
    diff_table_df = diff_table.view(pd.DataFrame)
    if len(diff_table_df) == 0:
        raise ValueError('No significant ASVs found to differentiate between %s and %s in field %s' % (val1, val2, field))
    res.append(diff_table)

    # print('creating diff ASV heatmap')
    # heatmap_func = ctx.get_action('dbbact', 'heatmap')
    # diff_heatmap = heatmap_func(table=diff_table, metadata=metadata, sort_field=field, cluster=False, repseqs=repseqs)
    # res.append(diff_heatmap)

    print('detecting enriched dbBact terms')
    enrich_func = ctx.get_action('dbbact', 'enrichment')
    enriched, = enrich_func(diff=diff_table, source='dsfdr', method=method, sig_threshold=sig_threshold, attack=attack, maxid=maxid, random_seed=random_seed)
    enriched_df = enriched.view(pd.DataFrame)
    print('found %d enriched dbbact terms' % len(enriched_df))
    res.append(enriched)

    print('creating enriched terms barplot')
    # create the lables for the 2 groups
    metadata_df = metadata.to_dataframe()
    if val1 is not None:
        label1 = ",".join(val1)
    if val2 is not None:
        label2 = ",".join(val2)
    if val1 is None and val2 is None:
        uvals = metadata_df[field].unique()
        label1 = uvals[0]
        label2 = uvals[1]
    elif val2 is None:
        label2 = 'NOT ' + label1

    terms_barplot_func = ctx.get_action('dbbact', 'plot_enrichment')
    enriched_barplot, = terms_barplot_func(enriched=enriched, labels=[label1, label2])
    res.append(enriched_barplot)

    print('done')
    return tuple(res)


def trim_primers(table: biom.Table, repseqs: DNAFASTAFormat = None) -> biom.Table:
    print('trim-primers')

    # incorporate the repseqs into the biom table if needed
    if repseqs is not None:
        print('converting hashes to sequences using repseqs file')
        table = embed_seqs(table, repseqs)

    seqs = table.ids(axis='observation')

    # check if table contains hashed sequences and repseqs not supplied
    if len(seqs[0]) == 32:
        if repseqs is None:
            raise ValueError('Table seems to contain hashes and not sequences. Please supply the --i-repseqs representative sequences file, or create the table (deblur/dada2) with the --p-no-hashed-feature-ids flag.')

    # take a random subset of ASVs to examine
    num_rand = np.min([100, len(seqs)])
    test_seqs = np.random.choice(seqs, num_rand, replace=False)

    # look if there is an exact match to any dbBact primer region
    region = _test_exact_region(test_seqs, min_fraction=0.5, ltrim=0)
    if region is not None:
        print('Found exact match for region %s. No need for trimming.' % region)
    else:
        # no exact match - so maybe primer is embedded?
        min_primer_len = 10
        primer, region = _test_embedded_primers(seqs, max_start=25, min_primer_len=min_primer_len, min_fraction=0.25)
        if region is not None:
            print('found embedded primer %s for region %s. removing it.' % (primer, region))
            # let's trim the primers
            trimmed_seqs_map = _trim_known_primer(seqs, primer, rprimer=None, length=0, remove_ambig=True, keep_primers=False)
        else:
            # no primer embedded. So maybe need to trim just first few bases?
            region = None
            for cpos in range(min_primer_len):
                region = _test_exact_region(test_seqs, min_fraction=0.5, ltrim=cpos)
                if region is not None:
                    break
            if region is not None:
                # found a matching region - lets trim it
                print('found exact region %s after %d left trmo' % (region, cpos))
                trimmed_seqs_map = {}
                for cseq in seqs:
                    trimmed_seqs_map[cseq] = cseq[cpos:]
            else:
                raise ValueError('No matching primers found. Are the reads reverse-complemented?')
        keep_ids = trimmed_seqs_map.keys()
        table.filter(keep_ids, axis='observation', inplace=True)
        table.update_ids(trimmed_seqs_map, axis='observation', strict=False, inplace=True)

    return table
