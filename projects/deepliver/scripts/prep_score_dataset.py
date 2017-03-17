#!/usr/bin/env python3
# coding=utf-8

import os as os
import sys as sys
import traceback as trb
import argparse as argp
import multiprocessing as mp

import numpy as np
import pandas as pd

import sklearn.feature_selection as fts


def parse_command_line():
    """
    :return:
    """
    parser = argp.ArgumentParser()
    parser.add_argument('--input', '-i', nargs='+', type=str, dest='inputfiles')
    parser.add_argument('--output', '-o', type=str, dest='outputfile')
    parser.add_argument('--workers', '-w', type=int, default=4, dest='workers')

    args = parser.parse_args()
    return args


def build_dataset(args):
    """
    :return:
    """
    hsa_data = None
    mmu_data = None
    for fp in args.inputfiles:
        is_human, is_mouse = False, False
        sample_id = os.path.basename(fp).split('.')[0]
        if sample_id.startswith('H'):
            is_human = True
            header = ['hsa_chrom', 'hsa_start', 'hsa_end', 'index', sample_id]
        elif sample_id.startswith('M'):
            is_mouse = True
            header = ['mmu_chrom', 'mmu_start', 'mmu_end', 'index', sample_id]
        else:
            raise ValueError('Unknown sample type: {}'.format(sample_id))
        datatypes = dict((h, d) for h, d in zip(header, [str, np.int32, np.int32, np.int32, np.float64]))
        sample_data = pd.read_csv(fp, delimiter='\t', header=None, names=header,
                                  dtype=datatypes, skip_blank_lines=True, na_values=['N/A'])
        if is_human:
            if hsa_data is None:
                hsa_data = sample_data
            else:
                hsa_data = hsa_data.merge(sample_data, on=['hsa_chrom', 'hsa_start', 'hsa_end', 'index'],
                                          how='outer')
        if is_mouse:
            if mmu_data is None:
                mmu_data = sample_data
            else:
                mmu_data = mmu_data.merge(sample_data, on=['mmu_chrom', 'mmu_start', 'mmu_end', 'index'],
                                          how='outer')
    hsa_data = hsa_data.fillna(value=0., inplace=False)
    hsa_data.index = hsa_data['index']
    hsa_data.drop('index', axis=1, inplace=True)
    mmu_data = mmu_data.fillna(value=0., inplace=False)
    mmu_data.index = mmu_data['index']
    mmu_data.drop('index', axis=1, inplace=True)

    merged_data = pd.concat([hsa_data, mmu_data], axis=1, ignore_index=False, join='outer')
    merged_data = merged_data.sort_index(axis=0, inplace=False, na_position='last')

    data_columns = [c for c in merged_data.columns if 'LiHe' in c]
    pos_columns = [c for c in merged_data.columns if 'LiHe' not in c]

    locations = merged_data.loc[:, pos_columns]
    data = merged_data.loc[:, data_columns]
    data = data.transpose()
    return data, locations


def compute_mi(args):
    """
    :param args:
    :return:
    """
    data, labels, target, neighbors = args
    mi = fts.mutual_info_classif(data, labels, discrete_features=False, n_neighbors=neighbors, copy=True)
    return target, mi


def compute_mutual_information(data, samples, features, workers):
    """
    :param data:
    :param samples:
    :param features:
    :return:
    """
    mutinf = []
    mutinf_labels = []
    sample_labels = []
    nn = []

    spec_labels = [1 if s.startswith('H') else 0 for s in samples]
    mutinf_labels.append('species')
    sample_labels.append(spec_labels)
    nn.append(5)

    sex_labels = [1 if s[1] == 'm' else 0 for s in samples]
    sample_labels.append(sex_labels)
    mutinf_labels.append('sex')
    nn.append(5)

    conditions = set(s.split('_')[-1] for s in samples)
    for c in conditions:
        cond_labels = [1 if s.endswith(c) else 0 for s in samples]
        sample_labels.append(cond_labels)
        mutinf_labels.append(c)
        nn.append(3)

    labeldf = pd.DataFrame(sample_labels, columns=samples, index=mutinf_labels, dtype=np.bool)
    labeldf = labeldf.transpose()

    params = [(data, l, t, n) for l, t, n in zip(sample_labels, mutinf_labels, nn)]

    with mp.Pool(workers) as pool:
        res = pool.imap_unordered(compute_mi, params)
        mutinf_labels = []
        for target, mi in res:
            mutinf.append(mi)
            mutinf_labels.append(target)

    midf = pd.DataFrame(mutinf, columns=features, index=mutinf_labels, dtype=np.float32)
    midf = midf.transpose()

    return midf, labeldf


def main():
    """
    :return:
    """
    args = parse_command_line()
    data, loc = build_dataset(args)
    mi, ind = compute_mutual_information(data, data.index, data.columns, args.workers)
    with pd.HDFStore(args.outputfile, 'w', complib='blosc', complevel=9) as hdf:
        hdf.put('data', data, format='fixed')
        hdf.put('loc', loc, format='fixed')
        hdf.put('mi', mi, format='fixed')
        hdf.put('indicator', ind, format='fixed')
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
