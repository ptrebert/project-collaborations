# coding=utf-8

import os as os
import datetime as dt
import fnmatch as fnm
import re as re
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


def set_sample_fraction(mergefiles, outbase, cmd, jobcall):
    """
    :param mergefiles:
    :param outbase:
    :param cmd:
    :param jobcall:
    :return:
    """
    # for narrow marks: ~40m reads (like H3K122ac); H3K4me1 a little more
    # for broad marks: ~60m reads
    # for Input: ~100m reads
    # DNase: ~200m reads
    fractions = {'Input': '-s 0.63', 'H3K122ac': '', 'H3K27ac': '-s 0.31',
                 'H3K4me3': '-s 0.84', 'H3K4me1': '-s 0.45',
                 'H3K27me3': '-s 0.32', 'H3K36me3': '-s 0.28',
                 'H3K9me3': '-s 0.31', 'DNase': '-s 0.52'}
    arglist = []
    for mrgf in mergefiles:
        lib = os.path.basename(mrgf).split('.')[0].split('_')[-1]
        prefix = 'TMP_' + lib
        frac = fractions[lib]
        tmp = cmd.format(**{'fraction': frac, 'prefix': prefix})
        outname = os.path.basename(mrgf).split('.')[0] + '.sort.bam'
        outfile = os.path.join(outbase, outname)
        arglist.append([mrgf, outfile, tmp, jobcall])
    if mergefiles:
        assert arglist, 'No arguments for BAM sampling created'
    return arglist


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
    if not inputfiles:
        return ''
    for inpf in inputfiles:
        fp, fn = os.path.split(inpf)
        mobj = matcher.match(fn)
        assert mobj is not None, 'Unexpected file: {}'.format(inpf)
        sample, lib = mobj.group('SAMPLE'), mobj.group('LIB')
        if not lib.startswith('H3K'):
            continue
        sigma, flen = extract_sigma_flen(inpf)
        sigmas[sample].append(sigma)
        params[sample + '_' + lib] = flen
    params['samples'] = list(sigmas.keys())
    for k, sig in sigmas.items():
        s = round(stat.median(sig))
        params[k + '_sigma'] = s
    overwrite = False
    if os.path.isfile(outputfile):
        with open(outputfile, 'r') as inf:
            old_params = js.load(inf)
        overwrite = old_params != params
    else:
        overwrite = True
    if overwrite:
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

    dir_rawdata = os.path.join(workdir, 'rawdata')
    rawdata_init = pipe.originate(task_func=lambda x: x,
                                  name='init',
                                  output=collect_full_paths(dir_rawdata, '*'))

    cmd = config.get('Pipeline', 'cleangenome')
    cleangenome = pipe.transform(task_func=sci_obj.get_jobf('in_out'),
                                 name='cleangenome',
                                 input=config.get('Refdata', 'chromreg'),
                                 filter=formatter(),
                                 output=config.get('Refdata', 'cleanreg'),
                                 extras=[cmd, jobcall]).follows(rawdata_init)

    dir_filtered = os.path.join(workdir, 'filtered')
    cmd = config.get('Pipeline', 'bamfilt')
    bamfilter = pipe.transform(task_func=sci_obj.get_jobf('in_out'),
                               name='bamfilter',
                               input=output_from(rawdata_init),
                               filter=suffix('.bam'),
                               output='.filt.bam',
                               output_dir=dir_filtered,
                               extras=[cmd, jobcall]).mkdir(dir_filtered).follows(cleangenome)

    cmd = config.get('Pipeline', 'countreads')
    countreads = pipe.transform(task_func=sci_obj.get_jobf('in_out'),
                                name='countreads',
                                input=output_from(bamfilter),
                                filter=suffix('.filt.bam'),
                                output='.filt.cnt',
                                output_dir=dir_filtered,
                                extras=[cmd, jobcall])

    sci_obj.set_config_env(dict(config.items('ParallelJobConfig')), dict(config.items('EnvConfig')))
    if args.gridmode:
        jobcall = sci_obj.ruffus_gridjob()
    else:
        jobcall = sci_obj.ruffus_localjob()

    # 01_HepG2_LiHG_Ct1_Input_S_1.bwa.20150122.filt.bam
    # 01_HepG2_LiHG_Ct1_Input_S_2.bwa.20150120.filt.bam
    dir_bammerge = os.path.join(workdir, 'merged')
    cmd = config.get('Pipeline', 'bammerge')
    bammerge = pipe.collate(task_func=sci_obj.get_jobf('ins_out'),
                            name='bammerge',
                            input=output_from(bamfilter),
                            filter=formatter('(?P<BIOSAMPLE>\w+)_(?P<REP>Ct[0-9])_(?P<MARK>\w+)_(?P<CENTER>(S|F)_(1|2)).+\.filt\.bam'),
                            output=os.path.join(dir_bammerge, '{BIOSAMPLE[0]}_{MARK[0]}.mrg.bam'),
                            extras=[cmd, jobcall]).mkdir(dir_bammerge)

    sci_obj.set_config_env(dict(config.items('JobConfig')), dict(config.items('EnvConfig')))
    if args.gridmode:
        jobcall = sci_obj.ruffus_gridjob()
    else:
        jobcall = sci_obj.ruffus_localjob()

    cmd = config.get('Pipeline', 'countreads')
    countmerged = pipe.transform(task_func=sci_obj.get_jobf('in_out'),
                                 name='countmerged',
                                 input=output_from(bammerge),
                                 filter=suffix('.mrg.bam'),
                                 output='.mrg.cnt',
                                 output_dir=dir_bammerge,
                                 extras=[cmd, jobcall])

    sci_obj.set_config_env(dict(config.items('ParallelJobConfig')), dict(config.items('EnvConfig')))
    if args.gridmode:
        jobcall = sci_obj.ruffus_gridjob()
    else:
        jobcall = sci_obj.ruffus_localjob()

    dir_bamsample = os.path.join(workdir, 'sampled')
    cmd = config.get('Pipeline', 'bamsample')
    merged_bam_files = collect_full_paths(dir_bammerge, '*.mrg.bam')
    params = set_sample_fraction(merged_bam_files, dir_bamsample, cmd, jobcall)
    bamsample = pipe.files(sci_obj.get_jobf('in_out'),
                           params,
                           name='bamsample').follows(bammerge).active_if(len(merged_bam_files) > 0).mkdir(dir_bamsample)

    sci_obj.set_config_env(dict(config.items('JobConfig')), dict(config.items('EnvConfig')))
    if args.gridmode:
        jobcall = sci_obj.ruffus_gridjob()
    else:
        jobcall = sci_obj.ruffus_localjob()

    cmd = config.get('Pipeline', 'bamidx')
    bamidx = pipe.transform(task_func=sci_obj.get_jobf('in_out'),
                            name='bamidx',
                            input=output_from(bamsample),
                            filter=suffix('.bam'),
                            output='.bai',
                            output_dir=dir_bamsample,
                            extras=[cmd, jobcall])

    sci_obj.set_config_env(dict(config.items('ParallelJobConfig')), dict(config.items('EnvConfig')))
    if args.gridmode:
        jobcall = sci_obj.ruffus_gridjob()
    else:
        jobcall = sci_obj.ruffus_localjob()

    cmd = config.get('Pipeline', 'nucfit')
    dir_nucfit = os.path.join(workdir, 'nuchunter', 'fit')
    nucfit = pipe.transform(task_func=sci_obj.get_jobf('in_out'),
                            name='nucfit',
                            input=output_from(bamsample),
                            filter=formatter('(?P<SAMPLE>\w+H3K[a-z0-9]+)\.(?P<EXT>[\w\.]+)\.bam$'),
                            output=os.path.join(dir_nucfit, '{SAMPLE[0]}.dat'),
                            extras=[cmd, jobcall]).mkdir(dir_nucfit).follows(bamidx)

    nucparams = pipe.merge(task_func=collect_fitpars_results,
                           name='nucparams',
                           input=output_from(nucfit),
                           output=os.path.join(dir_nucfit, 'sigma_fraglen_params.json'),
                           extras=['(?P<SAMPLE>\w+)_(?P<LIB>[A-Za-z0-9]+)\.dat'])

    sci_obj.set_config_env(dict(config.items('NodeJobConfig')), dict(config.items('EnvConfig')))
    if args.gridmode:
        jobcall = sci_obj.ruffus_gridjob()
    else:
        jobcall = sci_obj.ruffus_localjob()

    cmd = config.get('Pipeline', 'nuccall')
    dir_nuccall = os.path.join(workdir, 'nuchunter', 'call')
    nuccall = pipe.files(sci_obj.get_jobf('ins_out'),
                         build_callnucs_arguments(collect_full_paths(dir_bamsample, '*.bam'),
                                                  os.path.join(dir_nucfit, 'sigma_fraglen_params.json'),
                                                  dir_nuccall,
                                                  cmd, jobcall),
                         name='nuccall').mkdir(dir_nuccall).active_if(os.path.isfile(os.path.join(dir_nucfit, 'sigma_fraglen_params.json')))

    cmd = config.get('Pipeline', 'epiccount').replace('\n', ' ')
    dir_epiccount = os.path.join(workdir, 'epicseg', 'counts')
    epiccount = pipe.files(sci_obj.get_jobf('ins_out'),
                           build_epiccount_arguments(collect_full_paths(dir_bamsample, '*.bam'),
                                                     dir_epiccount, cmd, jobcall)).mkdir(dir_epiccount).follows(bamsample)

    cmd = config.get('Pipeline', 'epicnorm')
    epicnorm = pipe.transform(task_func=sci_obj.get_jobf('in_out'),
                              name='epicnorm',
                              input=output_from(epiccount),
                              filter=suffix('.txt'),
                              output='_norm.txt',
                              output_dir=dir_epiccount,
                              extras=[cmd, jobcall])

    dir_epicseg = os.path.join(workdir, 'epicseg', 'segment')
    cmd = config.get('Pipeline', 'epicseg').replace('\n', ' ')
    epicseg = pipe.transform(task_func=sci_obj.get_jobf('in_out'),
                             name='epicseg',
                             input=output_from(epicnorm),
                             filter=formatter('(?P<SAMPLE>\w+)\.counts_norm\.txt'),
                             output=os.path.join(dir_epicseg, '{SAMPLE[0]}_segmentation.bed'),
                             extras=[cmd, jobcall]).mkdir(dir_epicseg)

    sci_obj.set_config_env(dict(config.items('MemJobConfig')), dict(config.items('EnvConfigMacs')))
    if args.gridmode:
        jobcall = sci_obj.ruffus_gridjob()
    else:
        jobcall = sci_obj.ruffus_localjob()

    # perform DNase peak calling with MACS2
    dir_callpeaks = os.path.join(workdir, 'macs2', 'peaks')
    cmd = config.get('Pipeline', 'dnasepeak').replace('\n', ' ')
    dnasepeak = pipe.transform(task_func=sci_obj.get_jobf('in_out'),
                               name='dnasepeak',
                               input=output_from(bamsample),
                               filter=formatter('(?P<SAMPLE>01_HepG2_LiHG_DNase)\.sort\.bam'),
                               output=os.path.join(dir_callpeaks, '{SAMPLE[0]}_peaks.narrowPeak'),
                               extras=[cmd, jobcall]).mkdir(dir_callpeaks)

    dir_convncbi = os.path.join(workdir, 'ncbi')
    convncbi = pipe.transform(task_func=convert_to_ncbi,
                              name='convncbi',
                              input=output_from(epicseg, nuccall),
                              filter=suffix('.bed'),
                              output='.ncbi.bed',
                              output_dir=dir_convncbi).mkdir(dir_convncbi)

    convpeaks = pipe.transform(task_func=convert_to_ncbi,
                               name='convpeaks',
                               input=output_from(dnasepeak),
                               filter=suffix('.narrowPeak'),
                               output='.bed',
                               output_dir=dir_convncbi)

    sci_obj.set_config_env(dict(config.items('JobConfig')), dict(config.items('EnvConfig')))
    if args.gridmode:
        jobcall = sci_obj.ruffus_gridjob()
    else:
        jobcall = sci_obj.ruffus_localjob()

    dir_epicrep = os.path.join(workdir, 'epicseg', 'report')
    cmd = config.get('Pipeline', 'epicrep').replace('\n', ' ')
    epicrep = pipe.transform(task_func=sci_obj.get_jobf('in_out'),
                             name='epicrep',
                             input=output_from(convncbi),
                             filter=formatter('(?P<SAMPLE>\w+)_segmentation\.ncbi\.bed'),
                             output=os.path.join(dir_epicrep, '{SAMPLE[0]}_report.html'),
                             extras=[cmd, jobcall]).mkdir(dir_epicrep).follows(convpeaks)

    run_k122enh = pipe.merge(task_func=touch_checkfile,
                             name='run_k122enh',
                             input=output_from(rawdata_init, bamfilter, bammerge,
                                               bamidx, nucfit, nucparams, nuccall,
                                               epiccount, epicnorm, epicseg,
                                               convncbi, cleangenome, countreads,
                                               countmerged, dnasepeak, convpeaks,
                                               epicrep),
                             output=os.path.join(workdir, 'run_project_k122enh.chk'))

    return pipe
