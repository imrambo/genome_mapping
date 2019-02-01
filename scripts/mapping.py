#! /usr/bin/env python3
import glob
import os
import gzip
import re
from itertools import islice
from warnings import warn
import argparse
import subprocess
from mapping_functions import *
'''
2019 Ian Rambo
Thirteen... that's a mighty unlucky number... for somebody!
'''
#=============================================================================
#Command line options and Environment
parser = argparse.ArgumentParser()
parser.add_argument('--fastq_dir', dest='fastq_dir', type = str, nargs=1,
action='store', help='directory containing fastq files')
parser.add_argument('--genome', dest = 'genome', type = str, nargs = 1,
action = 'store', help='path to genome to map reads to')
parser.add_argument('--outdir', dest = 'outdir', type = str, nargs = 1,
action = 'store', help = 'path to output directory')
parser.add_argument('--sampleids', dest = 'sids', type = str, nargs = 1,
          action = 'store',
          help = """identifier for sample fastq files to be globbed, e.g. AB*.fastq.gz .
          Multiple identifiers can be specified in a single string when separated by commas, e.g. AB,MG,Megs""")
parser.add_argument('--netsam', dest = 'netsam', type = str, nargs = 1,
         action = 'store', help = 'SAM file from mapping a particular FASTQ file for use in network.pl. E.g. if you want to use mapping of FOO42_R1.fastq.gz and FOO42_R2.fastq.gz, specify --netsam=FOO42')
parser.add_argument('--bwa_thread', dest = 'bwa_thread', type = int, nargs = 1, action = 'store',
help = 'number of threads for bwa mem')
parser.add_argument('--samsort_thread', dest = 'samsort_thread', type = int, nargs = 1, action = 'store',
help = 'number of threads for samtools sort')
parser.add_argument('--samsort_mem', dest = 'samsort_mem', type = str, nargs = 1, action = 'store',
help = 'memory per thread for samtools sort. Specify an integer with K, M, or G suffix, e.g. 10G')
parser.add_argument('--nslice', dest = 'nslice', type = int, nargs = 1, action = 'store',
help = 'number of fastq headers to slice for determining if interleaved.')

args = vars(parser.parse_args())
#=============================================================================
#Options for bwa mem, samtools sort, depthfile
bwa_mem_opts = {'-t':args['bwa_thread']}
samtools_sort_opts = {'-@':args['samsort_thread'], '-m':args['samsort_mem']}
#Options for depth file
depthfile_opts = {'--noIntraDepthVariance':''}
#------------------------------------------------------------------------------
bwa_optstring = optstring_join(bwa_mem_opts)
samtools_sort_optstring = optstring_join(samtools_sort_opts)
depthfile_optstring = optstring_join(depthfile_opts)
#=============================================================================
###
### Main
###
#=============================================================================
#Create output directory if it does not exist
if not os.path.exists(args['outdir']):
    os.makedirs(args['outdir'])
else:
    pass
#=============================================================================
#Index the genome file
bwa_index(args['genome'])
#------------------------------------------------------------------------------
#Generate list of input FASTQ files using sample IDs
fastq_list = source_list_generator(args['sids'], args['fastq_dir'], '.fastq.gz')
#------------------------------------------------------------------------------
mapping_targets = list()
fastq_dict = find_fastq_pairs(fastq_list, nslice = args['nslice']*4)

for key in fastq_dict:
    maptarg = [os.path.join(args['outdir'], os.path.splitext(os.path.basename(fastq_dict[key]['R1']))[0] + x) for x in ['.reduced.sam', '.reduced.bam']]
    if fastq_dict[key]['R2'] == 'interleaved':
        bwa_sam_intl(targets = maptarg, sources = [args['genome'], fastq_dict[key]['R1']], bopts = bwa_optstring, sopts = samtools_sort_optstring)
        mapping_targets.extend(maptarg)

    else:
        bwa_sam_r1r2(targets = maptarg, sources = [args['genome'], fastq_dict[key]['R1'], fastq_dict[key]['R2']], bopts = bwa_optstring, sopts = samtools_sort_optstring)
        mapping_targets.extend(maptarg)
#------------------------------------------------------------------------------
#Depth file
depthfile_target = os.path.join(args['outdir'], os.path.splitext(os.path.basename(args['genome']))[0] + '_cov')
depthfile_sources = [m for m in mapping_targets if re.match(r'.*?\.reduced\.bam', m)]

depth_file(target = depthfile_target, sources = depthfile_sources, opts = depthfile_optstring)
#------------------------------------------------------------------------------
#Network file
network_source = [m for m in mapping_targets if args['netsam'] in m and m.endswith('.sam')][0]
network_target = os.path.join(args['outdir'], args['netsam'] + '.txt')

network_file(target = network_target, source = network_source)
