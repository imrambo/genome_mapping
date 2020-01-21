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
bwa_index_targets = [env['ASSEMBLY'] + ext for ext in ['.bwt','.pac','.ann','.amb','.sa']]
Default(bwa_index_targets)
#Index the assembly
env.BWA_index(bwa_index_targets, env['ASSEMBLY'])
#------------------------------------------------------------------------------
#Generate list of input FASTQ files using sample IDs
fastq_list = source_list_generator(env['SIDS'], env['FQDIR'], '.fastq.gz')
#------------------------------------------------------------------------------
mapping_targets = list()
fastq_dict = find_fastq_pairs(fastq_list, nslice = env['NSLICE'])

for key in fastq_dict:
    #maptarg = [os.path.splitext(os.path.basename(fastq_dict[key]['R1']))[0] + x for x in ['.reduced.sam', '.reduced.sorted.bam']]
    maptarg = [os.path.splitext(os.path.basename(fastq_dict[key]['R1']))[0] + '.reduced.sorted.bam']
    mapping_targets.extend(maptarg)
    Default(env.Install(env['OUTDIR'], maptarg))

    if 'R2' in fastq_dict[key].keys() and fastq_dict[key]['R2'] == 'interleaved':
        env.BWA_Samtools_Intl(maptarg, [env['ASSEMBLY'], fastq_dict[key]['R1']])

    elif 'R1' in fastq_dict[key].keys() and 'R2' in fastq_dict[key].keys() and fastq_dict[key]['R2'] != 'interleaved':
        env.BWA_Samtools_R1R2(maptarg, [env['ASSEMBLY'], fastq_dict[key]['R1'], fastq_dict[key]['R2']])

    else:
        pass
#------------------------------------------------------------------------------
assembly_id = os.path.splitext(os.path.basename(env['ASSEMBLY']))[0]
#Depth file
#depthfile_net_target = assembly_id + '_noIntDepthVar_cov'
depthfile_bin_target = assembly_id + '_cov'

depthfile_sources = [m for m in mapping_targets if re.match(r'.*?\.bam', m)]

#Default(env.Install(env['OUTDIR'], depthfile_net_target))
Default(env.Install(env['OUTDIR'], depthfile_bin_target))

#env.Depthfile_Net(depthfile_net_target, depthfile_sources)
env.Depthfile_Bin(depthfile_bin_target, depthfile_sources)
#------------------------------------------------------------------------------
#network_source = [m for m in mapping_targets if env['NETSAM'] in os.path.basename(m) and m.endswith('.sam')][0]

#network_prefix = re.findall(r'^(.*?)\..*?\.reduced\.sam', os.path.basename(network_source))[0] + '_%s' % assembly_id
#network_target = network_prefix + '_network.txt'

#Default(env.Install(env['OUTDIR'], network_target))
#env.Network(network_target, network_source)
