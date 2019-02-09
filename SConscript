Import('env')
import glob
import os
import gzip
import re
from itertools import islice
from warnings import warn
from fastq_pair import *

'''
2019 Ian Rambo
Thirteen... that's a mighty unlucky number... for somebody!
'''
bwa_index_targets = [env['GENOME'] + ext for ext in ['.bwt','.pac','.ann','.amb','.sa']]
Default(bwa_index_targets)
#Index the genome
env.BWA_index(bwa_index_targets, env['GENOME'])
#------------------------------------------------------------------------------
#Generate list of input FASTQ files using sample IDs
fastq_list = source_list_generator(env['SIDS'], env['FQDIR'], '.fastq.gz')
#------------------------------------------------------------------------------
mapping_targets = list()
fastq_dict = find_fastq_pairs(fastq_list, nslice = env['NSLICE'])

for key in fastq_dict:
    maptarg = [os.path.splitext(os.path.basename(fastq_dict[key]['R1']))[0] + x for x in ['.reduced.sam', '.reduced.sorted.bam']]
    mapping_targets.extend(maptarg)
    Default(env.Install(env['OUTDIR'], maptarg))

    if 'R2' in fastq_dict and fastq_dict[key]['R2'] == 'interleaved':
        env.BWA_Samtools_Intl(maptarg, [env['GENOME'], fastq_dict[key]['R1']])

    elif 'R1' in fastq_dict and 'R2' in fastq_dict and fastq_dict[key]['R2'] != 'interleaved':
        env.BWA_Samtools_R1R2(maptarg, [env['GENOME'], fastq_dict[key]['R1'], fastq_dict[key]['R2']])

    elif 'R1' in fastq_dict and not 'R2' in fastq_dict:
        warn('WARNING: Non-interleaved FASTQ: R1 present, R2 absent')
    elif not 'R1' in fastq_dict and 'R2' in fastq_dict:
        warn('WARNING: Non-interleaved FASTQ: R1 absent, R2 present')
    else:
        pass
#------------------------------------------------------------------------------
#Depth file
depthfile_net_target = os.path.splitext(os.path.basename(env['GENOME']))[0] + '_noIntDepthVar_cov'
depthfile_bin_target = os.path.splitext(os.path.basename(env['GENOME']))[0] + '_cov'

depthfile_sources = [m for m in mapping_targets if re.match(r'.*?\.bam', m)]

Default(env.Install(env['OUTDIR'], depthfile_net_target))
Default(env.Install(env['OUTDIR'], depthfile_bin_target))

env.Depthfile_Net(depthfile_net_target, depthfile_sources)
env.Depthfile_Bin(depthfile_bin_target, depthfile_sources)
#------------------------------------------------------------------------------
network_source = [m for m in mapping_targets if env['NETSAM'] in os.path.basename(m) and m.endswith('.sam')][0]
network_target = env['NETSAM'] + '_network.txt'

Default(env.Install(env['OUTDIR'], network_target))
env.Network(network_target, network_source)
