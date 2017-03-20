# coding=utf-8

import os as os
import datetime as dt
import fnmatch as fnm
import re as re
import shutil as sh
import csv as csv
import json as js
import statistics as stat
import collections as col
import io as io

from ruffus import *


def touch_checkfile(inputfiles, outputfile):
    """
    :param inputfiles:
    :param outputfile:
    :return:
    """
    timestr = dt.datetime.now().strftime('%Y-%m-%d@%A@%H:%M:%S')
    with open(outputfile, 'w') as outf:
        _ = outf.write(timestr + '\n')
    return outputfile


def collect_full_paths(rootdir, pattern):
    """
    :param rootdir:
    :param pattern:
    :return:
    """
    all_files = []
    for root, dirs, files in os.walk(rootdir):
        if files:
            filt = fnm.filter(files, pattern)
            for f in filt:
                all_files.append(os.path.join(root, f))
    return all_files


def copy_input_files(rootpaths, destdir):
    """
    :param rootpaths:
    :return:
    """
    if not os.path.isdir(destdir):
        os.makedirs(destdir, exist_ok=True)
    to_copy = []
    for fp in rootpaths:
        for root, dirs, files in os.walk(fp):
            if files:
                selected = fnm.filter(files, '*DNase*.nodup.blfilt.bamcov.bw')
                for s in selected:
                    to_copy.append(os.path.join(root, s))
    assert to_copy, 'No input files found'
    copied = []
    for fp in to_copy:
        dest = os.path.join(destdir, os.path.basename(fp))
        if not os.path.isfile(dest):
            sh.copy(fp, dest)
        copied.append(dest)
    return copied


def load_chromosomes(fpath):
    """
    :param fpath:
    :return:
    """
    chroms = dict()
    with open(fpath, 'r') as infile:
        for line in infile:
            if not line.strip():
                continue
            cols = line.strip().split()
            chroms[cols[0]] = int(cols[1])
    return chroms


def symm_filter_chainfile(inputfile, outpattern, chromsizes, cmd, jobcall):
    """
    :return:
    """
    arglist = []
    chromosomes = sorted(load_chromosomes(chromsizes).keys())
    for c in chromosomes:
        strfmt = {'chrom': c}
        c_out = outpattern.format(**strfmt)
        tmp = cmd.format(**strfmt)
        arglist.append([inputfile, c_out, tmp, jobcall])
    return arglist


def filter_bamfiles(bams, sample):
    """
    :param bams:
    :param sample:
    :return:
    """
    red = filter(lambda x: os.path.basename(x).startswith(sample), bams)
    histones = []
    control = []
    for fp in red:
        if '_Input' in fp:
            control.append(fp)
        elif '_H3K' in fp:
            histones.append(fp)
        else:
            continue
    assert len(histones) >= 6, 'Missing histone mark for sample {}: {}'.format(sample, histones)
    if len(control) > 1:
        control = [x for x in control if 'N.mrg.' in x]
    assert len(control) == 1, 'Could not identify control BAM: {}'.format(control)
    return histones, control


def build_callnucs_arguments(bamfiles, paramfile, outdir, cmd, jobcall):
    """
    :param bamfiles:
    :param paramfile:
    :param outdir:
    :param cmd:
    :param jobcall:
    :return:
    """
    try:
        params = js.load(open(paramfile, 'r'))
    except FileNotFoundError:
        return []
    samples = params['samples']
    arglist = []
    for smp in samples:
        hist, ctrl = filter_bamfiles(bamfiles, smp)
        sigma = params[smp + '_sigma']
        argd = {'sigma': sigma, 'prefix': smp, 'outdir': outdir, 'control': ctrl[0]}
        bamopts = ''
        for h in hist:
            fn = os.path.basename(h)
            smplib = fn.split('.')[0]
            flen = params[smplib]
            bamopts += '-in {} -fLen {} '.format(h, flen)
        argd['bamfiles'] = bamopts
        tmp = cmd.format(**argd)
        outfile = os.path.join(outdir, smp + '.bed')
        hist.extend(ctrl)
        arglist.append([hist, outfile, tmp, jobcall])
    if bamfiles:
        assert arglist, 'No arguments created'
    return arglist


def build_epiccount_arguments(bamfiles, outdir, cmd, jobcall):
    """
    :param bamfiles:
    :param outdir:
    :param cmd:
    :param jobcall:
    :return:
    """
    if not bamfiles:
        return []
    labeled_data = ''
    for bamf in bamfiles:
        label = os.path.basename(bamf).split('.')[0].split('_')[-1]
        this_file = '-m {}:{} '.format(label, bamf)
        labeled_data += this_file
    tmp = cmd.format(**{'bamfiles': labeled_data})
    outfile = os.path.join(outdir, '01_HepG2_LiHG.counts.txt')
    return [[bamfiles, outfile, tmp, jobcall]]


def convert_to_ncbi(inputfile, outputfile):
    """
    :param inputfile:
    :param outputfile:
    :return:
    """
    buffer = io.StringIO()
    chrom_match = re.compile('^[0-9XY]+\s')
    with open(inputfile, 'r') as infile:
        for line in infile:
            if not line.strip():
                continue
            if not chrom_match.match(line):
                buffer.write(line)
            else:
                buffer.write('chr' + line)
    with open(outputfile, 'w') as outfile:
        _ = outfile.write(buffer.getvalue())
    return outputfile


def build_pipeline(args, config, sci_obj):
    """
    :param args:
    :param config:
    :param sci_obj:
    :return:
    """

    pipe = Pipeline(name=config.get('Pipeline', 'name'))

    workdir = config.get('Pipeline', 'workdir')

    sci_obj.set_config_env(dict(config.items('JobConfig')), dict(config.items('EnvConfig')))
    if args.gridmode:
        jobcall = sci_obj.ruffus_gridjob()
    else:
        jobcall = sci_obj.ruffus_localjob()

    dir_chainfile = os.path.join(workdir, 'chain')
    chainfile_init = pipe.originate(task_func=lambda x: x,
                                    name='chainfile_init',
                                    output=collect_full_paths(dir_chainfile, '*.chain.gz'))

    dir_chain_temp = os.path.join(dir_chainfile, 'temp')
    cmd = config.get('Pipeline', 'symmfilt')
    params_symm_filter = symm_filter_chainfile(os.path.join(dir_chainfile, 'hg19.mm10.rbest.chain.gz'),
                                               os.path.join(dir_chain_temp, 'hg19_to_mm10.{chrom}.symm.tsv.gz'),
                                               os.path.join(config.get('Refdata', 'chromsizes'), 'mm10_chrom_augo.tsv'),
                                               cmd, jobcall)

    chain_symm_filt = pipe.files(sci_obj.get_jobf('in_out'),
                                 params_symm_filter,
                                 name='chain_symm_filt')
    chain_symm_filt = chain_symm_filt.mkdir(dir_chain_temp)
    chain_symm_filt = chain_symm_filt.follows(chainfile_init)

    cmd = config.get('Pipeline', 'mrgblocks')
    chain_symm_merge = pipe.merge(task_func=sci_obj.get_jobf('ins_out'),
                                  name='chain_symm_merge',
                                  input=output_from(chain_symm_filt),
                                  output=os.path.join(dir_chain_temp, 'hg19_to_mm10.wg.symm.tsv.gz'),
                                  extras=[cmd, jobcall])

    cmd = config.get('Pipeline', 'normblocks')
    chain_norm_filt = pipe.transform(task_func=sci_obj.get_jobf('in_out'),
                                     name='chain_norm_filt',
                                     input=output_from(chain_symm_merge),
                                     filter=suffix('wg.symm.tsv.gz'),
                                     output='150.symm.tsv.gz',
                                     output_dir=dir_chainfile,
                                     extras=[cmd, jobcall])

    cmd = config.get('Pipeline', 'hsa_blocks')
    hsa_blocks = pipe.transform(task_func=sci_obj.get_jobf('in_out'),
                                name='hsa_blocks',
                                input=output_from(chain_norm_filt),
                                filter=suffix('symm.tsv.gz'),
                                output='hsa_blocks.bed',
                                output_dir=dir_chainfile,
                                extras=[cmd, jobcall])

    cmd = config.get('Pipeline', 'mmu_blocks')
    mmu_blocks = pipe.transform(task_func=sci_obj.get_jobf('in_out'),
                                name='mmu_blocks',
                                input=output_from(chain_norm_filt),
                                filter=suffix('symm.tsv.gz'),
                                output='mmu_blocks.bed',
                                output_dir=dir_chainfile,
                                extras=[cmd, jobcall])

    # ============================
    # start processing DNase data

    dir_dnase_input = os.path.join(workdir, 'input')

    dnase_input = copy_input_files([config.get('DataSource', 'human'),
                                    config.get('DataSource', 'mouse')],
                                   dir_dnase_input)

    dnase_init = pipe.originate(task_func=lambda x: x,
                                name='dnase_init',
                                output=dnase_input)

    dir_dnase_bg = os.path.join(workdir, 'bedgraph')
    cmd = config.get('Pipeline', 'bwtobg')
    dnase_bwtobg = pipe.transform(task_func=sci_obj.get_jobf('in_out'),
                                  name='dnase_bwtobg',
                                  input=output_from(dnase_init),
                                  filter=formatter('41_(?P<SAMPLE>\w+)_DNase_S_.+\.bw'),
                                  output=os.path.join(dir_dnase_bg, '{SAMPLE[0]}.bg.gz'),
                                  extras=[cmd, jobcall])
    dnase_bwtobg = dnase_bwtobg.mkdir(dir_dnase_bg)

    dir_dnase_scores = os.path.join(workdir, 'scores')
    cmd = config.get('Pipeline', 'hsa_scores')
    hsa_scores = pipe.transform(task_func=sci_obj.get_jobf('in_out'),
                                name='hsa_scores',
                                input=output_from(dnase_bwtobg),
                                filter=formatter('(?P<SAMPLE>H(f|m)\w+)\.bg\.gz'),
                                output=os.path.join(dir_dnase_scores, '{SAMPLE[0]}.meansig.tsv'),
                                extras=[cmd, jobcall])
    hsa_scores = hsa_scores.mkdir(dir_dnase_scores)
    hsa_scores = hsa_scores.follows(hsa_blocks)

    cmd = config.get('Pipeline', 'mmu_scores')
    mmu_scores = pipe.transform(task_func=sci_obj.get_jobf('in_out'),
                                name='mmu_scores',
                                input=output_from(dnase_bwtobg),
                                filter=formatter('(?P<SAMPLE>M(f|m)\w+)\.bg\.gz'),
                                output=os.path.join(dir_dnase_scores, '{SAMPLE[0]}.meansig.tsv'),
                                extras=[cmd, jobcall])
    mmu_scores = mmu_scores.mkdir(dir_dnase_scores)
    mmu_scores = mmu_scores.follows(mmu_blocks)

    sci_obj.set_config_env(dict(config.items('ParallelJobConfig')), dict(config.items('EnvConfig')))
    if args.gridmode:
        jobcall = sci_obj.ruffus_gridjob()
    else:
        jobcall = sci_obj.ruffus_localjob()

    cmd = config.get('Pipeline', 'prep_mldata')
    dir_ml_data = os.path.join(workdir, 'mldata')
    prep_mldata = pipe.merge(task_func=sci_obj.get_jobf('ins_out'),
                             name='prep_mldata',
                             input=output_from(hsa_scores, mmu_scores),
                             output=os.path.join(dir_ml_data, 'hg19_mm10_mldata.h5'),
                             extras=[cmd, jobcall])
    prep_mldata = prep_mldata.mkdir(dir_ml_data)

    # cmd = config.get('Pipeline', 'cleangenome')
    # cleangenome = pipe.transform(task_func=sci_obj.get_jobf('in_out'),
    #                              name='cleangenome',
    #                              input=config.get('Refdata', 'chromreg'),
    #                              filter=formatter(),
    #                              output=config.get('Refdata', 'cleanreg'),
    #                              extras=[cmd, jobcall]).follows(rawdata_init)
    #
    # dir_filtered = os.path.join(workdir, 'filtered')
    # cmd = config.get('Pipeline', 'bamfilt')
    # bamfilter = pipe.transform(task_func=sci_obj.get_jobf('in_out'),
    #                            name='bamfilter',
    #                            input=output_from(rawdata_init),
    #                            filter=suffix('.bam'),
    #                            output='.filt.bam',
    #                            output_dir=dir_filtered,
    #                            extras=[cmd, jobcall]).mkdir(dir_filtered).follows(cleangenome)
    #
    # cmd = config.get('Pipeline', 'countreads')
    # countreads = pipe.transform(task_func=sci_obj.get_jobf('in_out'),
    #                             name='countreads',
    #                             input=output_from(bamfilter),
    #                             filter=suffix('.filt.bam'),
    #                             output='.filt.cnt',
    #                             output_dir=dir_filtered,
    #                             extras=[cmd, jobcall])
    #
    # sci_obj.set_config_env(dict(config.items('ParallelJobConfig')), dict(config.items('EnvConfig')))
    # if args.gridmode:
    #     jobcall = sci_obj.ruffus_gridjob()
    # else:
    #     jobcall = sci_obj.ruffus_localjob()
    #
    # # 01_HepG2_LiHG_Ct1_Input_S_1.bwa.20150122.filt.bam
    # # 01_HepG2_LiHG_Ct1_Input_S_2.bwa.20150120.filt.bam
    # dir_bammerge = os.path.join(workdir, 'merged')
    # cmd = config.get('Pipeline', 'bammerge')
    # bammerge = pipe.collate(task_func=sci_obj.get_jobf('ins_out'),
    #                         name='bammerge',
    #                         input=output_from(bamfilter),
    #                         filter=formatter('(?P<BIOSAMPLE>\w+)_(?P<REP>Ct[0-9])_(?P<MARK>\w+)_(?P<CENTER>(S|F)_(1|2)).+\.filt\.bam'),
    #                         output=os.path.join(dir_bammerge, '{BIOSAMPLE[0]}_{MARK[0]}.mrg.bam'),
    #                         extras=[cmd, jobcall]).mkdir(dir_bammerge)
    #
    # sci_obj.set_config_env(dict(config.items('JobConfig')), dict(config.items('EnvConfig')))
    # if args.gridmode:
    #     jobcall = sci_obj.ruffus_gridjob()
    # else:
    #     jobcall = sci_obj.ruffus_localjob()
    #
    # cmd = config.get('Pipeline', 'countreads')
    # countmerged = pipe.transform(task_func=sci_obj.get_jobf('in_out'),
    #                              name='countmerged',
    #                              input=output_from(bammerge),
    #                              filter=suffix('.mrg.bam'),
    #                              output='.mrg.cnt',
    #                              output_dir=dir_bammerge,
    #                              extras=[cmd, jobcall])
    #
    # sci_obj.set_config_env(dict(config.items('JobConfig')), dict(config.items('EnvConfig')))
    # if args.gridmode:
    #     jobcall = sci_obj.ruffus_gridjob()
    # else:
    #     jobcall = sci_obj.ruffus_localjob()
    #
    # cmd = config.get('Pipeline', 'bamidx')
    # bamidx = pipe.transform(task_func=sci_obj.get_jobf('in_out'),
    #                         name='bamidx',
    #                         input=output_from(bamsample),
    #                         filter=suffix('.bam'),
    #                         output='.bai',
    #                         output_dir=dir_bamsample,
    #                         extras=[cmd, jobcall])
    #
    # cmd = config.get('Pipeline', 'epiccount').replace('\n', ' ')
    # dir_epiccount = os.path.join(workdir, 'epicseg', 'counts')
    # epiccount = pipe.files(sci_obj.get_jobf('ins_out'),
    #                        build_epiccount_arguments(collect_full_paths(dir_bamsample, '*.bam'),
    #                                                  dir_epiccount, cmd, jobcall)).mkdir(dir_epiccount).follows(bamsample)
    #
    # cmd = config.get('Pipeline', 'epicnorm')
    # epicnorm = pipe.transform(task_func=sci_obj.get_jobf('in_out'),
    #                           name='epicnorm',
    #                           input=output_from(epiccount),
    #                           filter=suffix('.txt'),
    #                           output='_norm.txt',
    #                           output_dir=dir_epiccount,
    #                           extras=[cmd, jobcall])
    #
    # dir_epicseg = os.path.join(workdir, 'epicseg', 'segment')
    # cmd = config.get('Pipeline', 'epicseg').replace('\n', ' ')
    # epicseg = pipe.transform(task_func=sci_obj.get_jobf('in_out'),
    #                          name='epicseg',
    #                          input=output_from(epicnorm),
    #                          filter=formatter('(?P<SAMPLE>\w+)\.counts_norm\.txt'),
    #                          output=os.path.join(dir_epicseg, '{SAMPLE[0]}_segmentation.bed'),
    #                          extras=[cmd, jobcall]).mkdir(dir_epicseg)
    #
    # sci_obj.set_config_env(dict(config.items('MemJobConfig')), dict(config.items('EnvConfigMacs')))
    # if args.gridmode:
    #     jobcall = sci_obj.ruffus_gridjob()
    # else:
    #     jobcall = sci_obj.ruffus_localjob()
    #
    # # perform DNase peak calling with MACS2
    # dir_callpeaks = os.path.join(workdir, 'macs2', 'peaks')
    # cmd = config.get('Pipeline', 'dnasepeak').replace('\n', ' ')
    # dnasepeak = pipe.transform(task_func=sci_obj.get_jobf('in_out'),
    #                            name='dnasepeak',
    #                            input=output_from(bamsample),
    #                            filter=formatter('(?P<SAMPLE>01_HepG2_LiHG_DNase)\.sort\.bam'),
    #                            output=os.path.join(dir_callpeaks, '{SAMPLE[0]}_peaks.narrowPeak'),
    #                            extras=[cmd, jobcall]).mkdir(dir_callpeaks)
    #
    # run_k122enh = pipe.merge(task_func=touch_checkfile,
    #                          name='run_k122enh',
    #                          input=output_from(rawdata_init, bamfilter, bammerge,
    #                                            bamidx, epiccount, epicnorm, epicseg,
    #                                            cleangenome, countreads,
    #                                            countmerged, dnasepeak),
    #                          output=os.path.join(workdir, 'run_project_k122enh.chk'))

    return pipe
