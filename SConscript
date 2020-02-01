Import('env')
import glob
import os
import gzip
import re
from itertools import islice
from warnings import warn
from fastq_pair import *

'''
2020 Ian Rambo
Thirteen... that's a mighty unlucky number... for somebody!
'''
bwa_index_targets = [env['ASSEMBLY'] + ext for ext in ['.bwt','.pac','.ann','.amb','.sa']]
Default(bwa_index_targets)
#Index the assembly
env.BWA_index(bwa_index_targets, env['ASSEMBLY'])
#------------------------------------------------------------------------------
#Generate list of input FASTQ files using sample IDs
fastq_list = source_list_generator(env['SIDS'], env['FQDIR'], 'fastq')
#------------------------------------------------------------------------------
mapping_targets = list()
fastq_dict = find_fastq_pairs(fastq_list, nheader = env['NHEADER'])

extension = '.reduced.sorted.bam'
if env['MARKDUP']:
    extension = '.reduced.sorted.markdup.bam'
for key in fastq_dict:
    maptarg = [os.path.splitext(os.path.basename(fastq_dict[key]['R1']))[0] + extension]
    mapping_targets.extend(maptarg)
    Default(env.Install(env['OUTDIR'], maptarg))

    if 'R2' in fastq_dict[key].keys() and fastq_dict[key]['R2'] == 'interleaved':
        if env['MARKDUP']:
            print('Mark duplicates')
            env.BWA_Samtools_Markdup_Intl(maptarg, [env['ASSEMBLY'], fastq_dict[key]['R1'])
        else:
            env.BWA_Samtools_Intl(maptarg, [env['ASSEMBLY'], fastq_dict[key]['R1']])

    elif 'R1' in fastq_dict[key].keys() and 'R2' in fastq_dict[key].keys() and fastq_dict[key]['R2'] != 'interleaved':
        if env['MARKDUP']:
            print('Mark duplicates')
            env.BWA_Samtools_Markdup_R1R2(maptarg, [env['ASSEMBLY'], fastq_dict[key]['R1'], fastq_dict[key]['R2']])
        else:
            env.BWA_Samtools_R1R2(maptarg, [env['ASSEMBLY'], fastq_dict[key]['R1'], fastq_dict[key]['R2']])

    else:
        pass
#------------------------------------------------------------------------------
#get the basename of the reference assembly and use as its ID
assembly_id = get_basename(env['ASSEMBLY'])
#Depth file
depthfile_bin_target = assembly_id + '_cov'

depthfile_sources = [m for m in mapping_targets if re.match(r'.*?\.bam', m)]

Default(env.Install(env['OUTDIR'], depthfile_bin_target))

env.Depthfile_Bin(depthfile_bin_target, depthfile_sources)
#------------------------------------------------------------------------------
