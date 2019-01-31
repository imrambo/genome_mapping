import glob
import os
import gzip
import re
from itertools import islice
from warnings import warn
'''
2019 Ian Rambo
Thirteen... that's a mighty unlucky number... for somebody!
'''
#=============================================================================
def optstring_join(optdict):
    """
    Join a dictionary of command line options into a single string.
    """
    optstring = ' '.join([str(param) + ' ' + str(val) for param, val in optdict.items()])
    return optstring
#=============================================================================
#Command line options
AddOption('--fastq_dir', dest='fastq_dir', type='string', nargs=1,
action='store', help='directory containing fastq files')
AddOption('--genome', dest = 'genome', type = 'string', nargs = 1,
action = 'store', help='path to genome to map reads to')
AddOption('--outdir', dest = 'outdir', type = 'string', nargs = 1,
action = 'store', help = 'path to output directory')
AddOption('--sampleids', dest = 'sids', type = 'str', nargs = 1,
          action = 'store',
          help = '''identifier for sample fastq files to be globbed, e.g. AB*.fastq.gz .
          Multiple identifiers can be specified in a single string when separated by commas, e.g. AB,MG,Megs''')
AddOption('--netsam', dest = 'netsam', type = 'str', nargs = 1,
         action = 'store', help = 'SAM file from mapping a particular FASTQ file for use in network.pl. E.g. if you want to use mapping of FOO42_R1.fastq.gz and FOO42_R2.fastq.gz, specify --netsam=FOO42')
AddOption('--bwa_thread', dest = 'bwa_thread', type = 'int', nargs = 1, action = 'store',
help = 'number of threads for bwa mem')
AddOption('--samsort_thread', dest = 'samsort_thread', type = 'int', nargs = 1, action = 'store',
help = 'number of threads for samtools sort')
AddOption('--samsort_mem', dest = 'samsort_mem', type = 'str', nargs = 1, action = 'store',
help = 'memory per thread for samtools sort. Specify an integer with K, M, or G suffix, e.g. 10G')
AddOption('--nslice', dest = 'nslice', type = 'int', nargs = 1, action = 'store',
help = 'lines to slice from fastq file for determining if interleaved. MUST be a multiple of 4.')

BUILD_HELP =""" Usage: scons --fastq_dir=/path/to/fastq/files --genome=/path/to/genome
--outdir=/path/to/output/directory --interleaved=0/1 --sampleids=FOO,BAR,BAZ --netsam=FOO42
"""

if GetOption('help'):
    Help(BUILD_HELP)
#------------------------------------------------------------------------------
#Initialize environment
env = Environment(GENOME=GetOption('genome'),
                          FQDIR=GetOption('fastq_dir'),
                          OUTDIR=GetOption('outdir'),
                          SIDS=GetOption('sids'),
                          NETSAM=GetOption('netsam'),
                          NSLICE=GetOption('nslice'))

env.Replace(MAPDIR = os.path.join(env['OUTDIR'], 'mapping'),
NETDIR = os.path.join(env['OUTDIR'], 'network'))
#------------------------------------------------------------------------------
def find_fastq_pairs(fastq_list, nslice = 800):
    '''
    Match paired-end FASTQ reads using header information.
    '''
    if nslice % 4:
        warn('WARNING: --nslice is not a multiple of 4')

    fastq_dict = {}
    interleaved = False
    for fastq in fastq_list:
        with gzip.open(fastq, 'r') as fq:
            head_list = [l.decode('utf-8').split() for l in islice(fq, nslice) if l.decode('utf-8').startswith('@')]

            if all([head_list[n-1][0] == head_list[n][0] and head_list[n-1][1].startswith('1') and head_list[n][1].startswith('2') for n in range(1, len(head_list), 2)]):
                interleaved = True
            if all([head_list[n-1][0] != head_list[n][0] for n in range(1, len(head_list), 2)]):
                interleaved = False

            if interleaved:
                print('%s: interleaved FASTQ' % os.path.basename(fastq))
                if not head_list[0][0] in fastq_dict:
                    fastq_dict[head_list[0][0]] = {}
                    fastq_dict[head_list[0][0]]['R1'] = fastq
                    fastq_dict[head_list[0][0]]['R2'] = 'interleaved'
                else:
                    pass
            else:
                print('%s: non-interleaved FASTQ' % os.path.basename(fastq))
                if head_list[0][0] in fastq_dict:
                    if head_list[0][1].startswith('1'):
                        fastq_dict[head_list[0][0]]['R1'] = fastq
                    elif head_list[0][1].startswith('2'):
                        fastq_dict[head_list[0][0]]['R2'] = fastq
                    else:
                        print('please check if header is in format Casava 1.8+')
                else:
                    fastq_dict[head_list[0][0]] = {}
                    if head_list[0][1].startswith('1'):
                        fastq_dict[head_list[0][0]]['R1'] = fastq
                    elif head_list[0][1].startswith('2'):
                        fastq_dict[head_list[0][0]]['R2'] = fastq
                    else:
                        print('please check if header is in format Casava 1.8+')
    for key in fastq_dict:
        if not fastq_dict[key]['R2'] or not fastq_dict[key]['R1']:
            solo = [fastq_dict[key][v] for v in fastq_dict[key]][0]
            warn("""one is the loneliest number that you'll ever doooooo...
            so find a pair for %s if they are paired-end reads...""" % solo)

    return fastq_dict
#------------------------------------------------------------------------------
#Options for bwa mem, samtools sort
bwa_mem_opts = {'-t':GetOption('bwa_thread')}
samtools_sort_opts = {'-@':GetOption('samsort_thread'), '-m':GetOption('samsort_mem')}
#------------------------------------------------------------------------------
bwa_optstring = optstring_join(bwa_mem_opts)
samtools_sort_optstring = optstring_join(samtools_sort_opts)

#Builder for pipe: read mapping, SAM reduction, SAM to BAM
#FASTQ files are interleaved
bwa_samtools_intl_action = 'bwa mem %s ${SOURCES[0]} -p ${SOURCES[1]} | samtools view -hS -F4 - | tee ${TARGETS[0]} | samtools view -huS - | samtools sort %s - -o ${TARGETS[1]}' % (bwa_optstring, samtools_sort_optstring)
bwa_samtools_intl_builder = Builder(action = bwa_samtools_intl_action)

#Builder for pipe: read mapping, SAM reduction, SAM to BAM
#FASTQ files are separate R1 and R2
bwa_samtools_r1r2_action = 'bwa mem %s ${SOURCES[0]} ${SOURCES[1]} ${SOURCES[2]} | samtools view -hS -F4 - | tee ${TARGETS[0]} | samtools view -huS - | samtools sort %s - -o ${TARGETS[1]}' % (bwa_optstring, samtools_sort_optstring)
bwa_samtools_r1r2_builder = Builder(action = bwa_samtools_r1r2_action)
#------------------------------------------------------------------------------
#Builder for depthfile creation; additional options must be added to this dict
depthfile_opts = {'--noIntraDepthVariance':''}
depthfile_action = 'src/jgi_summarize_bam_contig_depth --outputDepth $TARGET $SOURCES %s' % optstring_join(depthfile_opts)
depthfile_builder = Builder(action = depthfile_action)
#------------------------------------------------------------------------------
network_action = 'perl src/network.pl -i $SOURCE -o $TARGET'
network_builder = Builder(action = network_action)
#------------------------------------------------------------------------------
builders = {'BWA_Samtools_Intl':bwa_samtools_intl_builder,
'BWA_Samtools_R1R2':bwa_samtools_r1r2_builder,
'Depthfile':depthfile_builder,
'Network':network_builder}
env.Append(BUILDERS = builders)

#env.SConscript(['mapping/SConscript'], variant_dir=mapdir, exports='env', duplicate=0)
#env.SConscript(['network/SConscript'], variant_dir=netdir, exports='env', duplicate=0)

#env.SConscript(['mapping/SConscript', 'network/SConscript'], exports='env')

#Export('env')
#------------------------------------------------------------------------------
####
####    MAIN
####
#------------------------------------------------------------------------------
#Index the genome file
bwa_index_targets = [os.path.abspath(env['GENOME']) + ext for ext in ['.bwt','.pac','.ann','.amb','.sa']]
Command(bwa_index_targets, env['GENOME'], 'bwa index $SOURCE')
#------------------------------------------------------------------------------
#Generate list of input FASTQ files using sample IDs
fastq_list = list()

if ',' in env['SIDS']:
    id_list = env['SIDS'].split(',')
    for i in id_list:
        idGlob = os.path.join(os.path.abspath(env['FQDIR']), '%s*.fastq.gz' % i)
        fastq_list.extend(list(glob.iglob(idGlob, recursive = False)))
else:
    idGlob = os.path.join(os.path.abspath(env['FQDIR']), '%s*.fastq.gz' % id)
    fastq_list.extend(list(glob.iglob(idGlob, recursive = False)))

env.Replace(FQLIST=fastq_list)
#------------------------------------------------------------------------------
maptarg = list()
fastq_dict = find_fastq_pairs(env['FQLIST'], nslice = env['NSLICE'])

for key in fastq_dict:
    mapping_targets = [os.path.join(mapdir, os.path.splitext(os.path.basename(fastq_dict[key]['R1']))[0] + x) for x in ['.reduced.sam', '.reduced.bam']]
    if fastq_dict[key]['R2'] == 'interleaved':
        env.BWA_Samtools_Intl(mapping_targets, [env['GENOME'], fastq_dict[key]['R1']])
        maptarg.extend(mapping_targets)

    else:
        env.BWA_Samtools_R1R2(mapping_targets, [env['GENOME'], fastq_dict[key]['R1'], fastq_dict[key]['R2']])
        maptarg.extend(mapping_targets)

env.Replace(MAPTARG = maptarg)
#------------------------------------------------------------------------------
#Depth file
depthfile_target = os.path.join(mapdir, os.path.splitext(os.path.basename(env['GENOME']))[0] + '_cov')
depthfile_sources = [m for m in env['MAPTARG'] if re.match(r'.*?\.reduced\.bam', m)]

env.Depthfile(depthfile_target, depthfile_sources)
#------------------------------------------------------------------------------
#Network file
network_source = [m for m in env['MAPTARG'] if env['NETSAM'] in m and m.endswith('.bam')][0]

network_target = os.path.join(env['OUTDIR'], env['NETSAM'] + '.txt')
env.Network(network_target, network_source)
