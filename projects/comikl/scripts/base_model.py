#!/usr/bin/env python3
# coding=utf-8

import os as os
import sys as sys
import traceback as trb
import argparse as argp
import itertools as itt
import collections as col
import re as re
import copy as cp
import functools as funct
import multiprocessing as mp
import pandas as pd
import time as ti
import random as rand
import pickle as pck
import sklearn as skl
import sklearn.model_selection as mods
from sklearn.ensemble import RandomForestClassifier


FEAT_REGLENGTH = 'ft_reglen'
FEAT_REPCONT = 'ft_repcont'
FEAT_KMERFREQ_PREFIX = 'ft_k'
FEAT_OECPG = 'ft_oecpg'
FEAT_COREPROM_PREFIX = 'ft_cpm_'


def feat_repetitive_content(region):
    """
    :param region:
    :return:
    """
    tmpseq = region['seq']
    seqlen = len(tmpseq)
    repmasked = tmpseq.count('a') + tmpseq.count('c') + tmpseq.count('g') + tmpseq.count('t')
    region[FEAT_REPCONT] = (repmasked / seqlen) * 100.
    return region


def feat_kmer_frequency(kmers, region):
    """
    :param region:
    :param kmers:
    :return:
    """
    tmpseq = region['seq'].upper()
    seqlen = len(tmpseq)
    for k in kmers:
        # normalization factor
        # number of possible substrings of length k in the seq
        # that is: n - k + 1 (k len substr in seq of len n)
        total_klen = float(seqlen - k + 1)
        kmerit = itt.product('ACGTN', repeat=k)
        kmerdict = dict()
        while 1:
            try:
                kmerdict[FEAT_KMERFREQ_PREFIX + ''.join(next(kmerit))] = 0
            except StopIteration:
                break
        wordfreqs = col.defaultdict(int)
        for i in range(0, k):
            # curly braces in literal part of format string
            # need to be escaped with curly braces
            words = re.findall('.{{{}}}'.format(k), tmpseq[i:])
            for w in words:
                wordfreqs[FEAT_KMERFREQ_PREFIX + w] += 1
        for key, val in wordfreqs.items():
            val = (val / total_klen) * 100.
            kmerdict[key] = val
        region.update(kmerdict)
    return region


def feat_oecpg_content(region):
    """
    :param region:
    :return:
    """
    seq = region['seq'].lower()
    seqlen = len(seq)  # for division
    total_G = seq.count('g')
    total_C = seq.count('c')
    total_CpG = seq.count('cg')
    # this definition of the obs-exp ratio is taken from UCSC
    region[FEAT_OECPG] = (total_CpG / (max(1, total_G) * max(1, total_C))) * seqlen
    return region


def feat_coreprom_motifs(region):
    """
    :param region:
     :type: dict
    :return:
    """

    core_motifs = []

    # DOI:10.1016/j.gene.2006.09.029
    # Yang et al., 2007 Gene
    # Refer specifically to TATA-less promoters

    core_motifs.append(('elkM3', '[GC]CGGAAG[CT]'))  # not sure if that is a reasonable one
    core_motifs.append(('sp1M6', 'GGGCGG[AG]'))  # not sure if that is a reasonable one
    core_motifs.append(('novelM22', 'TGCGCA[ACGTN][GT]'))

    # DOI:10.1016/j.ydbio.2009.08.009
    # Juven-Gershon, Kadonaga 2010, Dev. Bio
    # DOI:10.1101/gad.1026202
    # Butler, Kadonaga 2002, Genes & Dev.

    core_motifs.append(('tataM3', 'TATA[AT]AA[AG]'))  # TATA box
    core_motifs.append(('inrM4', '[TC][TC]A[ACGTN][AT][TC][TC]'))  # initiator (Inr)
    core_motifs.append(('breMx', '[GC][GC][AG]CGCC'))  # TFIIB recognition element (BRE)
    core_motifs.append(('dpeM9', '[AG]G[AT][TC][GAC](T)?'))  # downstream core promoter element (DPE)
    core_motifs.append(('mteM10', 'C[GC]A[AG]C[GC][GC]AACG[GC]'))  # motif ten (MTE)

    # DOI:10.1093/nar/gkv1032
    # Marbach-Bar et al., 2016 NAR
    core_motifs.append(('dtieMx', 'G[CGT][CGT][AG][AGT][ACGTN][ACT]GG'))  # Downstream Transcription Initiation Element (DTIE)

    tmpseq = region['seq'].upper()
    reglen = len(tmpseq)
    for name, motifre in core_motifs:
        occurrences = re.findall(motifre, tmpseq)
        bpcov = sum(len(m) for m in occurrences)
        region[FEAT_COREPROM_PREFIX + 'pct_' + name] = (bpcov / reglen) * 100
        region[FEAT_COREPROM_PREFIX + 'abs_' + name] = len(occurrences)
    return region


def _get_feat_fun_map():
    feat_fun_map = {'prm': feat_coreprom_motifs,
                    'oecpg': feat_oecpg_content,
                    'rep': feat_repetitive_content,
                    'kmf': feat_kmer_frequency}
    return feat_fun_map


def _make_kmer_dict(k, alphabet='ACGTN'):
    """
    :param k:
    :param alphabet:
    :return:
    """
    kmerit = itt.product(alphabet, repeat=k)
    kmers = dict()
    while 1:
        try:
            kmers[FEAT_KMERFREQ_PREFIX + ("".join(next(kmerit)))] = 0
        except StopIteration:
            break
    return kmers


def feat_single_kmer(k, kmers, region):
    """
    :param k:
    :param kmers:
    :param region:
    :return:
    """
    # TODO
    # the deep copy is probably necessary for the intended
    # use case, see get_online_version
    mykmers = cp.deepcopy(kmers)
    tmpseq = region['seq'].upper()  # kmers are not case-sensitive
    seqlen = len(tmpseq)
    total_klen = float(seqlen - k + 1)
    wordfreqs = col.defaultdict(int)
    for i in range(0, k):
        # curly braces in literal part of format string
        # need to be escaped with curly braces
        words = re.findall('.{{{}}}'.format(k), tmpseq[i:])
        for w in words:
            wordfreqs[FEAT_KMERFREQ_PREFIX + w] += 1
    for key, val in wordfreqs.items():
        val = (val / total_klen) * 100.
        mykmers[key] = val
    region.update(mykmers)
    return region


def compute_features(exec_functions, d):
    """
    :param exec_functions:
    :param d:
    :return:
    """
    for f in exec_functions:
        d = f(d)
    return d


def get_online_version(features, kmers=None):
    """ Return a closure encompassing all feature
    functions - intended to use is with
    multiprocessing.Pool.map() or similar
    """
    if 'kmf' in features:
        assert kmers is not None, 'No values for k-mer frequency specified'
    funmap = _get_feat_fun_map()
    exec_functions = set()
    for ft in features:
        if ft == 'kmf':
            for k in kmers:
                kd = _make_kmer_dict(k)
                part = funct.partial(feat_single_kmer, *(k, kd))
                exec_functions.add(part)
        else:
            exec_functions.add(funmap[ft])
    return exec_functions


def process_datafile(fpath, label):
    """
    :param fpath:
    :return:
    """
    regions = []
    this_region = None
    compfeat = get_online_version(['prm', 'oecpg', 'rep', 'kmf'], [2, 3, 4])
    with open(fpath, 'r') as infile:
        for line in infile:
            if not line.strip():
                continue
            if line.startswith('>'):
                if this_region is not None:
                    #this_region = compfeat(this_region)
                    regions.append(this_region)
                this_region = dict()
                this_region['name'] = line.strip().lstrip('>')
            else:
                seq = line.strip()
                this_region['seq'] = seq
                this_region[FEAT_REGLENGTH] = len(seq)
                this_region['label'] = label
    regions.append(this_region)
    featfun = get_online_version(['prm', 'oecpg', 'rep', 'kmf'], [2, 3, 4])
    compfeat = funct.partial(compute_features, *(featfun, ))
    with mp.Pool(6) as pool:
        regions = pool.map(compfeat, regions)
    return regions


def parse_command_line():
    """
    :return:
    """
    parser = argp.ArgumentParser()
    parser.add_argument('--pos', '-p', type=str, dest='positives', required=True)
    parser.add_argument('--neg', '-n', type=str, dest='negatives', required=True)
    parser.add_argument('--out-dir', '-o', type=str, dest='outdir', default=os.getcwd())
    args = parser.parse_args()
    return args


def create_dataset(args):
    """
    :param args:
    :return:
    """
    pos_regions = process_datafile(args.positives, 1)
    neg_regions = process_datafile(args.negatives, 0)
    outpath = os.path.join(args.outdir, 'dataset_raw.h5')
    pos = pd.DataFrame.from_dict(pos_regions, orient='columns')
    neg = pd.DataFrame.from_dict(neg_regions, orient='columns')
    full = pd.concat([pos, neg], ignore_index=False)
    with pd.HDFStore(outpath, 'w', complevel=9, complib='blosc') as hdf:
        hdf.put('/dataset', full, format='fixed')
        hdf.flush()
    labels = full['label']
    full.drop(['seq', 'label'], axis='columns', inplace=True)
    # note that "full" still contains name column
    return full, labels


def dump_dataset(train_data, train_labels, test_data, test_labels, outpath):
    """
    :param train_data:
    :param train_labels:
    :param test_data:
    :param test_labels:
    :param outpath:
    :return:
    """
    with pd.HDFStore(outpath, 'w', complevel=9, complib='blosc') as hdf:
        hdf.put('/train/data', train_data, format='fixed')
        hdf.put('/train/labels', train_labels, format='fixed')
        hdf.put('/test/data', test_data, format='fixed')
        hdf.put('/test/labels', test_labels, format='fixed')
        hdf.flush()
    return


def main():
    """
    :return:
    """
    rand.seed()
    args = parse_command_line()
    data, labels = create_dataset(args)
    for idx in range(10):
        x_train, x_test, y_train, y_test = mods.train_test_split(data, labels,
                                                                 stratify=labels,
                                                                 test_size=0.33,
                                                                 random_state=rand.randint(0, 1000000))
        data_out = os.path.join(args.outdir, 'dataset_run{}.h5'.format(idx))
        dump_dataset(x_train, y_train, x_test, y_test, data_out)
        x_train = x_train.drop('name', axis='columns', inplace=False)
        x_test = x_test.drop('name', axis='columns', inplace=False)
        print('Run data saved - CV starts {}'.format(ti.ctime()))
        param_grid = {'n_estimators': [500, 750, 1000, 2000], 'min_samples_split': [2, 4, 6, 8],
                      'oob_score': [False], 'class_weight': ['balanced']}
        cvgrid = mods.GridSearchCV(estimator=RandomForestClassifier(),
                                   param_grid=param_grid,
                                   scoring='average_precision',
                                   n_jobs=15, pre_dispatch=15, cv=10, refit=True)
        cvgrid.fit(x_train, y_train)
        probs = cvgrid.predict_proba(x_test, y_test)
        score = cvgrid.score(x_test, y_test)
        print('Run {} finished: {}'.format(idx, ti.ctime()))
        print(score)
        print(cvgrid.best_params_)
        print('===')
        model_out = os.path.join(args.outdir, 'model_run{}.pck'.format(idx))
        with open(model_out, 'wb') as dump:
            pck.dump({'cv': cvgrid, 'probs': probs}, dump)
        print('Data saved')
    return


if __name__ == '__main__':
    try:
        main()
    except Exception as err:
        trb.print_exc(file=sys.stderr)
        sys.stderr.write('\nError: {}\n'.format(str(err)))
        sys.exit(1)
    else:
        sys.exit(0)
