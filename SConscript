Import('env')
import glob
import os
import gzip
import re
from itertools import islice
from warnings import warn
from fastq_pair import *
#------------------------------------------------------------------------------
#Index the genome file
bwa_index_targets = [os.path.abspath(env['GENOME']) + ext for ext in ['.bwt','.pac','.ann','.amb','.sa']]
Command(bwa_index_targets, env['GENOME'], 'bwa index $SOURCE')
#------------------------------------------------------------------------------
#Generate list of input FASTQ files using sample IDs
fastq_list = source_list_generator(env['SIDS'], env['FQDIR'], '.fastq.gz')

#env.Replace(FQLIST=fastq_list)
#------------------------------------------------------------------------------
mapping_targets = list()
fastq_dict = find_fastq_pairs(fastq_list, nslice = env['NSLICE'])

for key in fastq_dict:
    maptarg = [os.path.join(env['OUTDIR'], os.path.splitext(os.path.basename(fastq_dict[key]['R1']))[0] + x) for x in ['.reduced.sam', '.reduced.bam']]
    #mapping_targets = [os.path.join(env['MAPDIR'], os.path.splitext(os.path.basename(fastq_dict[key]['R1']))[0] + x) for x in ['.reduced.sam', '.reduced.bam']]
    if fastq_dict[key]['R2'] == 'interleaved':
        env.BWA_Samtools_Intl(maptarg, [env['GENOME'], fastq_dict[key]['R1']])
        mapping_targets.extend(maptarg)

    else:
        env.BWA_Samtools_R1R2(maptarg, [env['GENOME'], fastq_dict[key]['R1'], fastq_dict[key]['R2']])
        mapping_targets.extend(maptarg)

#env.Replace(MAPTARG = maptarg)
#------------------------------------------------------------------------------
#Depth file
depthfile_target = os.path.join(env['OUTDIR'], os.path.splitext(os.path.basename(env['GENOME']))[0] + '_cov')
#depthfile_target = os.path.join(env['MAPDIR'], os.path.splitext(os.path.basename(env['GENOME']))[0] + '_cov')

#depthfile_sources = [m for m in env['MAPTARG'] if re.match(r'.*?\.reduced\.bam', m)]
depthfile_sources = [m for m in mapping_targets if re.match(r'.*?\.reduced\.bam', m)]

env.Depthfile(depthfile_target, depthfile_sources)
#------------------------------------------------------------------------------
#network_source = [m for m in env['MAPTARG'] if env['NETSAM'] in m and m.endswith('.bam')][0]
network_source = [m for m in mapping_targets if env['NETSAM'] in m and m.endswith('.bam')][0]

#network_target = os.path.join(env['NETDIR'], env['NETSAM'] + '.txt')
network_target = os.path.join(env['OUTDIR'], env['NETSAM'] + '.txt')
env.Network(network_target, network_source)
