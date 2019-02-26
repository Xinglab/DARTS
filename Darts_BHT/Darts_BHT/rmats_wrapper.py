# -*- coding: utf-8 -*-

'''Transferred parser code from rMATS-turbo (v4.0.2). 

This version of rMATS-turbo was compiled using Python 2.7 and Cython 0.27, with unicode UCS4.

If you encounter errors in importing the `rmatspipeline.so` module, please refer to
the rMATS user guide page at http://rnaseq-mats.sourceforge.net/user_guide.htm
'''

import os
import sys
import logging
import argparse as ap
import datetime
import shutil
import tempfile
from subprocess import call, check_output

VERSION = 'v4.0.2'

def add_rmats_parser( subparsers ):
    parser = subparsers.add_parser("rmats_count", help="Darts_BHT rmats_count: run rMATS-turbo to count junction reads in BAM files")

    group1 = parser.add_mutually_exclusive_group()
    group2 = parser.add_mutually_exclusive_group()

    parser.add_argument('--version', action='version',
                        help='Version.', version=VERSION)
    parser.add_argument('--gtf', action='store', required=True,
                        help='An annotation of genes and transcripts in GTF format.', dest='gtf')

    group1.add_argument('--b1', action='store', default='',
                        help='BAM configuration file.', dest='b1')
    group2.add_argument('--b2', action='store', default='',
                        help='BAM configuration file.', dest='b2')
    group1.add_argument('--s1', action='store', default='',
                        help='FASTQ configuration file.', dest='s1')
    group2.add_argument('--s2', action='store', default='',
                        help='FASTQ configuration file.', dest='s2')

    parser.add_argument('--od', action='store',
                        help='output folder of post step.', dest='od')
    parser.add_argument('-t', action='store', default='paired',
                        choices=['paired', 'single',],
                        help='readtype, single or paired.', dest='readtype')
    parser.add_argument('--libType', action='store', default='fr-unstranded',
                        choices=['fr-unstranded', 'fr-firststrand',
                                 'fr-secondstrand',],
                        help='Library type. Default is unstranded (fr-unstranded). Use fr-firststrand or fr-secondstrand for strand-specific data.', dest='dt')
    parser.add_argument('--readLength', action='store', type=int, default=0,
                        help='The length of each read.', dest='readLength')
    parser.add_argument('--anchorLength', action='store', type=int, default=1,
                        help='The anchor length. (default is 1.)', dest='anchorLength')
    parser.add_argument('--tophatAnchor', action='store', type=int, default=6,
                        help='The "anchor length" or "overhang length" used in the aligner. At least “anchor length” NT must be mapped to each end of a given junction. The default is 6. (This parameter applies only if using fastq).', dest='tophatAnchor')
    parser.add_argument('--bi', action='store', default='',
                        help='The folder name of the STAR binary indexes (i.e., the name of the folder that contains SA file). For example, use ~/STARindex/hg19 for hg19. (Only if using fastq)', dest='bIndex')
    parser.add_argument('--nthread', action='store', type=int, default=1,
                        help='The number of thread. The optimal number of thread should be equal to the number of CPU core.', dest='nthread')

    #parser.add_argument('--tstat', action='store', type=int, default=1,
    #                    help='the number of thread for statistical model.', dest='tstat')
    #parser.add_argument('--cstat', action='store', type=float, default=0.0001,
    #                    help='The cutoff splicing difference. The cutoff used in the null hypothesis test for differential splicing. The default is 0.0001 for 0.01%% difference. Valid: 0 ≤ cutoff < 1.', dest='cstat')

    #parser.add_argument('--statoff', action='store_false',
    #                    help='Turn statistical analysis off.', dest='stat')
    return

def process_parsed_rmats_args( args ):
    args.tmp = tempfile.mkdtemp()
    if args.b1 == '' and args.b2 == '' and args.s1 == '' and args.s2 == '':
        print('ERROR: BAM/FASTQ required. Please check b1, b2, s1 and s2.')
        exit(0)
    if args.gtf == '' or args.od == '':
        print('ERROR: GTF file and output folder required. Please check --gtf and --od.')
        exit(0)
    if (args.s1 != '' or args.s2 != '') and args.bIndex == '':
        print('ERROR: STAR binary indexes required. Please check --bi.')
        exit(0)

    if len(args.b1) > 0:
        if args.b1.endswith('.txt'):
            with open(args.b1, 'r') as fp:
                args.b1 = fp.read().strip(' ,\n')
        else:
            tmp = args.b1.split(',')
            assert all( (x.endswith('.bam') for x in tmp) ), Exception('--b1 must be either a .txt configuration file, or a list of comma-separated .bam alignment files')
    if len(args.b2) > 0:
        if args.b2.endswith('.txt'):
            with open(args.b2, 'r') as fp:
                args.b2 = fp.read().strip(' ,\n')
        else:
            tmp = args.b2.split(',')
            assert all( (x.endswith('.bam') for x in tmp) ), Exception('--b2 must be either a .txt configuration file, or a list of comma-separated .bam alignment files')
    if len(args.s1) > 0:
        with open(args.s1, 'r') as fp:
            args.s1 = fp.read().strip(' ,\n')
    if len(args.s2) > 0:
        with open(args.s2, 'r') as fp:
            args.s2 = fp.read().strip(' ,\n')

    #if args.stat and (len(args.b1) * len(args.b2) == 0 and\
    #                  len(args.s1) * len(args.s2) == 0):
    #    print 'ERROR: while performing statistical analysis, user should provide two groups of samples. Please check b1,b2 or s1,s2.'
    #    exit(0)

    if args.b1 == '' and args.b2 == '' and (args.s1 != '' or args.s2 != ''):
        args.b1, args.b2 = doSTARMapping(args)

    args.bams = ','.join([args.b1, args.b2]).strip(',')
    args.junctionLength = 2 * (args.readLength - args.anchorLength)

    dt_map = {'fr-unstranded':0, 'fr-firststrand':1, 'fr-secondstrand':2}
    args.dt = dt_map[args.dt]

    # add back stat-related
    args.tstat = 1
    args.cstat = 0.00001
    args.stat = True
    return args


def doSTARMapping(args): ## do STAR mapping
    fastqs = [args.s1, args.s2,]
    bams = [[], [],]

    for i in range(len(fastqs)):
        if fastqs[i] != '':
            sample = [pair.split(':') for pair in fastqs[i].split(',')]
            print("mapping the first sample")
            for rr, pair in enumerate(sample):
                map_folder = os.path.join(args.tmp, 'bam%d_%d' % (i+1, rr+1));

                if os.path.exists(map_folder):
                    if os.path.isdir(map_folder):
                        os.rmdir(map_folder)
                    else:
                        os.unlink(map_folder)

                os.makedirs(map_folder)
                cmd = 'STAR --chimSegmentMin 2 --outFilterMismatchNmax 3 --alignEndsType EndToEnd --runThreadN 4 --outSAMstrandField intronMotif --outSAMtype BAM SortedByCoordinate ';
                cmd += '--alignSJDBoverhangMin ' + str(args.tophatAnchor) + ' --alignIntronMax 299999 --genomeDir ' + args.bIndex + ' --sjdbGTFfile ' + args.gtf; 
                cmd += ' --outFileNamePrefix ' + map_folder + '/ --readFilesIn ';
                cmd += ' '.join(pair)
                status, output = check_output(cmd, shell=True);
                print("mapping sample_%d, %s is done with status %s" % (i, ' '.join(pair), status))
                if (int(status)!=0): ## it did not go well
                    print("error in mapping sample_%d, %s: %s" % (i, ' '.join(pair),status))
                    print("error detail: %s" % output)
                    raise Exception();
                print(output)
                bams[i].append(os.path.join(map_folder, 'Aligned.sortedByCoord.out.bam'))

    return ','.join(bams[0]), ','.join(bams[1])
##### end of doSTARMapping ####

