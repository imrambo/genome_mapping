import glob
import os
import gzip
import re
from itertools import islice
from warnings import warn
import sys
import atexit

'''
2019 Ian Rambo
Thirteen... that's a mighty unlucky number... for somebody!
'''
#Add scripts directory to path
sys.path.insert(0, os.path.abspath('./src'))
#=============================================================================
def optstring_join(optdict):
    """
    Join a dictionary of command line options into a single string.
    """
    optstring = ' '.join([str(param) + ' ' + str(val) for param, val in optdict.items()])
    return optstring
#------------------------------------------------------------------------------
def remove_build_targets(tmpdir):
    """
    Remove intermediate build targets within a specified temporary directory.
    """
    if os.path.exists(tmpdir):
        print('removing intermediate build targets in %s' % os.path.abspath(tmpdir))
        for tmp in [os.path.join(tmpdir, os.path.basename(str(t))) for t in BUILD_TARGETS]:
            if os.path.isfile(tmp):
                print('removing %s' % tmp)
                os.remove(tmp)
            else:
                pass

        if not os.listdir(tmpdir):
            print('removing empty directory: %s' % tmpdir)
            os.rmdir(tmpdir)
        else:
            print('directory %s is not empty' % tmpdir)
    else:
        print('%s does not exist' % tmpdir)
        pass
    return None
#=============================================================================
#Command line options and Environment
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
help = 'number of headers from fastq file for determining if interleaved. Must be even.')
AddOption('--rm_local_build', dest = 'rmbuild', type = 'int', nargs = 1,
action = 'store', default = 0, help = 'only keep the build targets in the --outdir. Will remove build targets in the temporary build within SConstruct directory. Specify 0 (keep) or 1 (remove)')
#------------------------------------------------------------------------------
#Initialize environment
env = Environment(GENOME=GetOption('genome'),
                          FQDIR=GetOption('fastq_dir'),
                          OUTDIR=GetOption('outdir'),
                          SIDS=GetOption('sids'),
                          NETSAM=GetOption('netsam'),
                          NSLICE=GetOption('nslice'))
#=============================================================================
###
### Builders
###

#OPTION DICTIONARIES
#Options for bwa mem, samtools sort, jgi coverage
bwa_mem_opts = {'-t':GetOption('bwa_thread')}
samtools_sort_opts = {'-@':GetOption('samsort_thread'), '-m':GetOption('samsort_mem')}
depthfile_opts_net = {'--noIntraDepthVariance':''}
#------------------------------------------------------------------------------
#BWA index builder, add index targets as default targets
bwa_index_builder = Builder(action = 'bwa index $SOURCE')

#Option strings for bwa mem, samtools sort
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
depthfile_net_action = 'src/jgi_summarize_bam_contig_depths --outputDepth $TARGET $SOURCES %s' % optstring_join(depthfile_opts_net)
depthfile_net_builder = Builder(action = depthfile_net_action)

depthfile_bin_action = 'src/jgi_summarize_bam_contig_depths --outputDepth $TARGET $SOURCES'
depthfile_bin_builder = Builder(action = depthfile_bin_action)
#------------------------------------------------------------------------------
#Builder for network file
network_action = 'perl src/network.pl -i $SOURCE -o $TARGET'
network_builder = Builder(action = network_action)
#------------------------------------------------------------------------------
builders = {'BWA_Samtools_Intl':bwa_samtools_intl_builder,
'BWA_Samtools_R1R2':bwa_samtools_r1r2_builder,
'Depthfile_Net':depthfile_net_builder,
'Depthfile_Bin':depthfile_bin_builder,
'Network':network_builder,
'BWA_index':bwa_index_builder}
env.Append(BUILDERS = builders)
#=============================================================================
#SConscript
if env['GENOME']:
    build_tmp = os.path.splitext(os.path.basename(env['GENOME']))[0] + '_build'
    SConscript(['SConscript'], exports='env', variant_dir=build_tmp, duplicate=0)
#------------------------------------------------------------------------------
#If --rmbuild=1, remove the build targets in the temporary directory
if GetOption('rmbuild'):
    atexit.register(remove_build_targets, tmpdir = build_tmp)

else:
    pass

Export('env')
