from collections import defaultdict
import re

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

    Parameters
    ----------

    Returns
    -------
    data: pandas.DataFrame
        with the updated index containing sequences instead of hashes
    '''
    # check if feature ids are hashes by looking at first feature
    if len(data.iloc[0].name) == 32:
        print('replacing hashes with sequences')
        rename_ids = {}
        for chead, cseq in _iter_fasta(str(repseqs)):
            if chead in data.index:
                rename_ids[chead] = cseq
        if len(rename_ids) < len(repseqs):
            print('Missing rep. sequences for some hashes. Found %d out of %d' % (len(rename_ids), len(repseqs)))
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


def _test_exact_region(seqs, kmers={'v4': ['TACG'], 'v3': ['TGGG', 'TGAG'], 'v1': ['GACG', 'GATG', 'ATTG']}, min_fraction=0.5, ltrim=0):
    '''Test if a fasta file starts with known region k-mers

    Parameters
    ----------
    seqs: list of str
        the sequences to test
    kmers: dict of {region_name(str): [kmers(str)]}
        the kmers expected to begin each primer region
        get the values using ~/scripts/count_kmer.py -d v1
    num_reads: int
        the maximal number of reads per fileto test
    min_fraction: float
        the minimal expected fraction of reads in the region starting with the kmers (summed over all kmers of the region)
    ltrim: int, optional
        position of first nucleotide to start with. 0 to start from beginning

    Returns
    -------
    region: str or None
        the matching region or None if no region
    '''
    kmer_len = len(kmers['v4'][0])
    print('testing kmer head for region on %d sequences:' % len(seqs))
    kmer_dist = defaultdict(float)
    for cseq in seqs:
        cseq = cseq[ltrim:ltrim + kmer_len]
        for cregion, ckmers in kmers.items():
            if cseq in ckmers:
                kmer_dist[cregion] += 1

    if len(kmer_dist) > 0:
        maxregion = max(kmer_dist, key=kmer_dist.get)
        if kmer_dist[maxregion] / len(seqs) >= min_fraction:
            return maxregion
        else:
            print('No region exact match detected. Maximal region was %s but fraction of matching sequences is %f.' % (maxregion, kmer_dist[maxregion] / len(seqs)))
            return None
    print('No region exact match detected. No sequence matched any of the regions.')
    return None


def _test_embedded_primers(seqs, primers={'AGAGTTTGATC[AC]TGG[CT]TCAG': 'v1', 'CCTACGGG[ACGT][CGT]GC[AT][CG]CAG': 'v3', 'GTGCCAGC[AC]GCCGCGGTAA': 'v4'}, max_start=25, min_primer_len=10, min_fraction=0.25):
    '''Check if the sequences contain one of a given set of primers in the beginning.

    Parameters
    ----------
    seqs: list of str
        the sequences to test
    primers: dict of {primer(str): region_name(str)}, optional
        the primers to test for
    max_start: int, optional
        maximal start position for the primer (i.e. do not return if primer starts after position max_start)
    min_primer_len: int, optional
        trim primers to keep only min_primer_len last chars
    num_reads: int, optional
        the number of reads in the fasta file to process
    min_fraction: float, optional
        need at least min_fraction sequence matches in order to return the primer

    Returns
    -------
    primer: str or None
        the primer identified as maximal (if >min_fraction matches in tested reads)
    primer_name: str or None
        the name of the primer region identified
    '''
    # attach the base_dir if needed
    print('Testing %d sequences for embedded primers (%d)' % (len(seqs), len(primers)))

    # trim the primers if needed
    if min_primer_len is not None:
        print('Trimming primers before test to length %d' % min_primer_len)
        new_primers = {}
        for k, v in primers.items():
            pos = len(k)
            numchars = 0
            newp = ''
            while True:
                if numchars >= min_primer_len:
                    break
                pos = pos - 1
                if pos < 0:
                    break
                if k[pos] != ']':
                    newp = k[pos] + newp
                    numchars += 1
                    continue
                while k[pos] != '[':
                    newp = k[pos] + newp
                    pos = pos - 1
                newp = k[pos] + newp
                numchars += 1
            new_primers[newp] = v
        primers = new_primers
        print('Trimmed primers are: %s' % primers)

    # scan the files
    matches = defaultdict(float)
    for cseq in seqs:
        cseq = cseq.upper()
        for cprimer in primers.keys():
            ccseq = cseq[:max_start + len(cprimer)]
            match = re.search(cprimer, ccseq)
            if match is not None:
                matches[cprimer] += 1

    if len(matches) == 0:
        print('No embedded primers found in sequences')
        return None, None

    max_primer = max(matches, key=matches.get)
    match_fraction = matches[max_primer] / len(seqs)
    print('Found maximal matching primer: %s with enough matches (%f)' % (max_primer, match_fraction))
    if match_fraction < min_fraction:
        print('Not using primer. Maximal matches are not enough.')
        return None, None
    return max_primer, primers[max_primer]


def _trim_known_primer(seqs, fprimer, rprimer=None, length=150, remove_ambig=True, keep_primers=False):
    '''Trim until a known primer in the sequences

    Parameters
    ----------
    seqs: list of Str
        the sequences to trim ('ACGT')
    fprimer: str
        the forward primer to use ('ACGT')
    rprimer: str or None, optional
        the reverse primer to use ('ACGT') or None to not match reverse (trim only according to fprimer and length)
    length: int, optional
        the number of bases after the end of the primer to keep
    remove_ambig: bool, optional
    keep_primers: bool, optional
        If True, keep the primers trim until the beginning of the fprimer). If False, start sequences at the end of fprimer.

    Returns
    -------
    newseqs: dict, where key is the seqs_df index, and value is the trimmed sequence
    NOTE: sequences not containing the primer are not included in the output
    '''
    print('trimming %d sequences with primer %s' % (len(seqs), fprimer))
    newseqs = {}
    for ccseq in seqs:
        cseq = ccseq.upper()
        reverse_primer_seq = ''
        # find the start of the primer and output all following sequence
        try:
            match = re.search(fprimer, cseq)
            if keep_primers:
                forward_primer_seq = cseq[match.start():]
                fplen = match.end() - match.start()
            else:
                forward_primer_seq = cseq[match.end():]
        except AttributeError:
            forward_primer_seq = ''

        # if sequence can be amplified, search for reverse primer
        if forward_primer_seq != '':
            if rprimer is not None:
                # find the reverse primer match and use only up to it
                try:
                    match = re.search(rprimer, forward_primer_seq)
                    if keep_primers:
                        reverse_primer_seq = forward_primer_seq[:match.end()]
                    else:
                        reverse_primer_seq = forward_primer_seq[:match.start()]
                except AttributeError:
                    reverse_primer_seq = ''
            else:
                # reverse_primer_seq = forward_primer_seq[:length + fplen]
                reverse_primer_seq = forward_primer_seq[:]

            if reverse_primer_seq != '':
                # trim length if needed
                if length > 0:
                    if keep_primers:
                        reverse_primer_seq = reverse_primer_seq[:length + fplen]
                    else:
                        reverse_primer_seq = reverse_primer_seq[:length]
                printit = True

                if remove_ambig:
                    if 'N' in reverse_primer_seq:
                        printit = False

                if printit:
                    newseqs[ccseq] = reverse_primer_seq
    print('trimmed %d sequences, %d not found' % (len(newseqs), len(seqs) - len(newseqs)))
    return newseqs
