#!/usr/bin/env python3
# coding=utf-8

import os as os
import sys as sys
import traceback as trb
import argparse as argp

import numpy as np
import pandas as pd

def parse_command_line():
    """
    :return:
    """
    parser = argp.ArgumentParser()
    parser.add_argument('--input', '-i', nargs='+', type=str, dest='inputfiles')
    parser.add_argument('--output', '-o', type=str, dest='outputfile')

    args = parser.parse_args()
    return args


def main():
    """
    :return:
    """
    args = parse_command_line()
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
    merged_data = merged_data.sort(axis=0, inplace=False)
    locations = merged_data.loc[:, [c for c in merged_data.columns if 'LiHe' not in c]]
    data = merged_data.loc[:, [c for c in merged_data.columns if 'LiHe' in c]]
    print(data.shape)
    blocks_var = data.var(axis=1)
    data = data.loc[blocks_var > 0.5, :]
    print(data.shape)
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
