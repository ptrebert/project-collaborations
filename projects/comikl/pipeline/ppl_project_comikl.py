# coding=utf-8

import os as os
import datetime as dt
import fnmatch as fnm
import re as re
import csv as csv
import json as js
import statistics as stat

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
    sigmas = []
    matcher = re.compile(regexp)
    for inpf in inputfiles:
        fp, fn = os.path.split(inpf)
        mobj = matcher.match(fn)
        assert mobj is not None, 'Unexpected file: {}'.format(inpf)
        sample, lib = mobj.group('SAMPLE'), mobj.group('LIB')
        if lib == 'Input':
            continue
        sigma, flen = extract_sigma_flen(inpf)
        sigmas.append(sigma)
        params[sample + '_' + lib] = flen
    sigma = int(stat.median(sigmas))
    params['sigma'] = sigma
    with open(outputfile, 'w') as outf:
        js.dump(params, outf)
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

    run_comikl = pipe.merge(task_func=touch_checkfile,
                            name='run_comikl',
                            input=output_from(rawdata_init, bamfilter, bammerge,
                                              bamidx, nucfit, nucparams),
                            output=os.path.join(workdir, 'run_project_comikl.chk'))

    return pipe