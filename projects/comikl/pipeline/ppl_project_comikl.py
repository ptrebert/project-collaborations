# coding=utf-8

import os as os
import datetime as dt
import fnmatch as fnm
import re as re
import csv as csv
import json as js
import statistics as stat
import collections as col

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


def extract_sigma_flen(fpath):
    """
    :return:
    """
    sigma, mu, auc = 0, 0, 0
    with open(fpath, 'r', newline='') as datfile:
        rows = csv.DictReader(datfile, delimiter='\t')
        for r in rows:
            if float(r['auc']) > auc:
                auc = float(r['auc'])
                sigma = float(r['t-sigma'])
                mu = float(r['mu'])
    assert sigma > 0 and mu > 0, 'NucHunter dat file Parsing failed: {}'.format(fpath)
    return sigma, int(mu)


def collect_fitpars_results(inputfiles, outputfile, regexp):
    """
    :param inputfiles:
    :param outputfile:
    :param regexp:
    :return:
    """
    params = dict()
    sigmas = col.defaultdict(list)
    matcher = re.compile(regexp)
    for inpf in inputfiles:
        fp, fn = os.path.split(inpf)
        mobj = matcher.match(fn)
        assert mobj is not None, 'Unexpected file: {}'.format(inpf)
        sample, lib = mobj.group('SAMPLE'), mobj.group('LIB')
        if lib == 'Input':
            continue
        sigma, flen = extract_sigma_flen(inpf)
        sigmas[sample].append(sigma)
        params[sample + '_' + lib] = flen
    params['samples'] = list(sigmas.keys())
    for k, sig in sigmas.items():
        s = int(stat.median(sig))
        params[k + '_sigma'] = s
    with open(outputfile, 'w') as outf:
        js.dump(params, outf, indent=1, sort_keys=True)
    return outputfile


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
        if '_Input_' in fp:
            control.append(fp)
        else:
            histones.append(fp)
    assert len(histones) == 6, 'Missing histone mark for sample {}: {}'.format(sample, histones)
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
    params = js.load(open(paramfile, 'r'))
    samples = params['samples']
    arglist = []
    for smp in samples:
        hist, ctrl = filter_bamfiles(bamfiles, smp)
        sigma = params[smp + '_sigma']
        argd = {'sigma': sigma, 'prefix': smp, 'outdir': outdir, 'control': ctrl[0]}
        bamopts = ''
        for h in hist:
            fn = os.path.basename(h)
            smplib = fn.split('.')[0].rsplit('_', 2)[0]
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

    dir_rawdata = os.path.join(workdir, 'rawdata')
    rawdata_init = pipe.originate(task_func=lambda x: x,
                                  name='init',
                                  output=collect_full_paths(dir_rawdata, '*'))

    dir_filtered = os.path.join(workdir, 'filtered')
    cmd = config.get('Pipeline', 'bamfilt')
    bamfilter = pipe.transform(task_func=sci_obj.get_jobf('in_out'),
                               name='bamfilter',
                               input=output_from(rawdata_init),
                               filter=suffix('.bam'),
                               output='.filt.bam',
                               output_dir=dir_filtered,
                               extras=[cmd, jobcall]).mkdir(dir_filtered)

    sci_obj.set_config_env(dict(config.items('ParallelJobConfig')), dict(config.items('EnvConfig')))
    if args.gridmode:
        jobcall = sci_obj.ruffus_gridjob()
    else:
        jobcall = sci_obj.ruffus_localjob()

    # 01_HepG2_LiHG_Ct1_Input_S_1.bwa.20150122.filt.bam
    # 01_HepG2_LiHG_Ct1_Input_S_2.bwa.20150120.filt.bam
    cmd = config.get('Pipeline', 'bammerge')
    bammerge = pipe.collate(task_func=sci_obj.get_jobf('ins_out'),
                            name='bammerge',
                            input=output_from(bamfilter),
                            filter=formatter('(?P<SAMPLE>\w+_Input_S)_(?P<REP>[0-9]+)\.bwa\.(?P<DATE>[0-9]+)\.(?P<EXT>[\w\.]+)'),
                            output=os.path.join(dir_filtered, '{SAMPLE[0]}_N.mrg.20161219.{EXT[0]}'),
                            extras=[cmd, jobcall])

    sci_obj.set_config_env(dict(config.items('JobConfig')), dict(config.items('EnvConfig')))
    if args.gridmode:
        jobcall = sci_obj.ruffus_gridjob()
    else:
        jobcall = sci_obj.ruffus_localjob()

    cmd = config.get('Pipeline', 'bamidx')
    bamidx = pipe.transform(task_func=sci_obj.get_jobf('in_out'),
                            name='bamidx',
                            input=output_from(bamfilter, bammerge),
                            filter=suffix('.bam'),
                            output='.bai',
                            output_dir=dir_filtered,
                            extras=[cmd, jobcall])

    sci_obj.set_config_env(dict(config.items('ParallelJobConfig')), dict(config.items('EnvConfig')))
    if args.gridmode:
        jobcall = sci_obj.ruffus_gridjob()
    else:
        jobcall = sci_obj.ruffus_localjob()

    cmd = config.get('Pipeline', 'nucfit')
    dir_nucfit = os.path.join(workdir, 'nucfit')
    nucfit = pipe.transform(task_func=sci_obj.get_jobf('in_out'),
                            name='nucfit',
                            input=output_from(bammerge, bamfilter),
                            filter=formatter('(?P<SAMPLE>\w+)\.(?P<EXT>[\w\.]+)\.bam$'),
                            output=os.path.join(dir_nucfit, '{SAMPLE[0]}.dat'),
                            extras=[cmd, jobcall]).mkdir(dir_nucfit)

    nucparams = pipe.merge(task_func=collect_fitpars_results,
                           name='nucparams',
                           input=output_from(nucfit),
                           output=os.path.join(dir_nucfit, 'sigma_fraglen_params.json'),
                           extras=['(?P<SAMPLE>\w+)_(?P<LIB>[A-Za-z0-9]+)_(S|F)_(1|2|N)\.dat'])

    sci_obj.set_config_env(dict(config.items('MemJobConfig')), dict(config.items('EnvConfig')))
    if args.gridmode:
        jobcall = sci_obj.ruffus_gridjob()
    else:
        jobcall = sci_obj.ruffus_localjob()

    cmd = config.get('Pipeline', 'nuccall')
    dir_nuccall = os.path.join(workdir, 'nuccall')
    nuccall = pipe.files(sci_obj.get_jobf('ins_out'),
                         build_callnucs_arguments(collect_full_paths(dir_filtered, '*.bam'),
                                                  os.path.join(dir_nucfit, 'sigma_fraglen_params.json'),
                                                  dir_nuccall,
                                                  cmd, jobcall),
                         name='nuccall').mkdir(dir_nuccall)

    run_comikl = pipe.merge(task_func=touch_checkfile,
                            name='run_comikl',
                            input=output_from(rawdata_init, bamfilter, bammerge,
                                              bamidx, nucfit, nucparams, nuccall),
                            output=os.path.join(workdir, 'run_project_comikl.chk'))

    return pipe