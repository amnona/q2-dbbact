from q2_types.feature_data import DNAFASTAFormat
import pandas as pd
import numpy as np

import calour as ca


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
        print('table does not seem to contain qiime2 hashes for ASV ids.\nFirst ASV id is %s\n' % data.iloc[0].name)
    return data


def _load_diff_abundance(diff, diff_tsv, source, ancom_stat, repseqs=None, sig_threshold=0.1):
    '''Load a diff_abundance qiime2 artifact anc convert it to a calour experiment

    Supports multiple formats (songbird / aldex2 / ancom / tsv)

    Parameters
    ----------
    diff: biom.Table or None
        the differential abundance biom table (or None)
    diff_tsv: str or None
        A tsv file name containing the diff. abundance table
    source: str
        the source for the diff. abundance output. options are:
            "songbird"
            "aldex2"
            "ancom"
            "tsv": use an external tsv file
    ancom_stat: str, optional
        name of the ancom_stat output file (for source=ancom)
    repseqs: q2_types.feature_data.DNAFASTAFormat or None
        if not None, incorporate representative sequences from this file (instead of the hashes in the table)
    sig_threshold: float
        the significance threshold for including sequences from the diff. abundance test

    Returns
    -------
    ca.AmpliconExperiment
        containing the diff. abundance features and the direction each feature is associated with
    '''
    # handle either diff (qza diffabundance artifact) or diff_tsv(tsv file resulting from diff abundance export i.e. ancom etc.)
    if diff is None and diff_tsv is None:
        raise ValueError('Need to supply either --i-diff or --p-diff-tsv differential results input file (none were supplied)')
    if diff is not None and diff_tsv is not None:
        raise ValueError('Need to supply either --i-diff or --p-diff-tsv differential results input file (both were supplied)')

    # get the differential abundance dataframe. index is sequence
    if diff is not None:
        data = diff
    else:
        try:
            data = pd.read_csv(diff_tsv, sep='\t')
            data.set_index(data.columns[0], drop=False, inplace=True)
        except Exception as e:
            print('failed to read diff_tsv file %s. Is is a valid tsv differential abundance file?' % diff_tsv)
            raise e

    # if needed, replace q2 hashes with sequences from the rep-seqs file
    if repseqs is not None:
        data = _seqs_from_repseqs(data=data, repseqs=repseqs)

    ndata = data.copy()

    # check if feature ids are hashes by looking at first feature
    if len(data.iloc[0].name) == 32:
        print('ASV IDs seem to be hashes. First id=%s' % data.iloc[0].name)
        raise ValueError('Input file contains sequence hashes instead of actual sequences.\nPlease use qiime dbbact enrichment-hash and supply the rep_seqs.qza.')

    # preprocess the dataframe to capture the correct per-feature differential abundance fields
    # these will be translated into "effect" (effect size),"pval" (p-value for differential abundance), "dir" (higher in which group) and "reject" (reject the null hypothesis of similar distribution in both groups)
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

    # create a new calour experiment using sequences as features and 1 sample (we don't use the sample)
    exp = ca.AmpliconExperiment(data=np.ones([1, len(ndata)]),
                                sample_metadata=pd.DataFrame({'_sample_id': ['s1']}),
                                feature_metadata=ndata,
                                sparse=False)
    return exp
