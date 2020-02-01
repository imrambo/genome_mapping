import glob
import os
import gzip
import re
from itertools import islice
from warnings import warn
import sys
import atexit

'''
Motivation: Map reads to a assembly with bwa mem and get SAM and BAM output.
SAM and BAM will not include unmapped reads.
Depth file can be toggled with or without intra-depth variance.

Author: Ian Rambo
Contact: ian.rambo@utexas.edu
Thirteen... that's a mighty unlucky number... for somebody!

v.1.0.0.2
'''
EnsurePythonVersion(3, 5)
#Add scripts directory to path
sys.path.append(os.path.abspath('./src'))
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
    if os.path.exists(tmpdir) and os.listdir(tmpdir):
        print('removing intermediate build targets in %s' % os.path.abspath(tmpdir))
        for tmp in [os.path.join(tmpdir, os.path.basename(str(t))) for t in BUILD_TARGETS]:
            if os.path.isfile(tmp):
                print('removing %s' % tmp)
                os.remove(tmp)
            else:
                pass

        if not os.listdir(tmpdir):
            print('removing empty directory: "%s"' % tmpdir)
            os.rmdir(tmpdir)
        else:
            print('directory "%s" is not empty' % tmpdir)
    else:
        print('Cannot delete directory "%s", does not exist' % tmpdir)
        pass
    return None
#------------------------------------------------------------------------------
def get_basename(file_path):
    """
    Get file basename, excluding multiple extensions.
    """
    basename = os.path.basename(file_path)
    #Remove two extensions, e.g. foo.tar.gz becomes foo
    if re.match(r'^.*?\.[a-z]+\.[a-z]+$', basename):
        basename = re.findall(r'^(.*?)\.[a-z]+\.[a-z]+$', basename)[0]
    else:
        basename = os.path.splitext(basename)[0]
    return basename
#=============================================================================
#Command line options and Environment
AddOption('--fastq_dir', dest='fastq_dir', type='string', nargs=1,
action='store', help='directory containing fastq files')
AddOption('--assembly', dest = 'assembly', type = 'string', nargs = 1,
action = 'store', help='path to assembly to map reads to')
AddOption('--outdir', dest = 'outdir', type = 'string', nargs = 1,
action = 'store', help = 'path to output directory')
AddOption('--sampleids', dest = 'sids', type = 'str', nargs = 1,
          action = 'store',
          help = 'identifier for sample fastq files to be globbed, e.g. AB for AB*.fastq.gz. Multiple identifiers can be specified in a single string when separated by commas, e.g. AB,MG,Megs')
AddOption('--align_thread', dest = 'align_thread', type = 'int', nargs = 1, action = 'store',
help = 'number of threads for alignment algorithm')
AddOption('--samsort_thread', dest = 'samsort_thread', type = 'int', nargs = 1, action = 'store',
help = 'number of threads for samtools sort')
AddOption('--samsort_mem', dest = 'samsort_mem', type = 'str', nargs = 1, action = 'store',
help = 'memory per thread for samtools sort. Specify an integer with K, M, or G suffix, e.g. 10G')
AddOption('--nheader', dest = 'nheader', type = 'int', nargs = 1, action = 'store',
help = 'number of headers from fastq file for determining if interleaved.')
AddOption('--tmpdir', dest = 'tmpdir', type = 'str', nargs = 1, action = 'store',
help = 'output directory for samtools sort temporary files')
AddOption('--rm_local_build', dest = 'rmbuild', type = 'int', nargs = 1,
action = 'store', default = 0, help = 'only keep the build targets in the --outdir. Will remove build targets in the temporary build within SConstruct directory. Specify 0 (keep) or 1 (remove). Default is 0.')
AddOption('--noIntraDepthVariance', dest = 'nointdepth', type = 'int', nargs = 1,
action = 'store', default = 1, help = 'toggle jgi_summarize_bam_contig_depths --noIntraDepthVariance (yes = 1, no = 0). Default = 1')
AddOption('--markdup', dest = 'markdup', type = 'int', nargs = 1,
action = 'store', default = 0, help = 'choose to fix mates and mark duplicates for paired-end reads (yes = 1, no = 0). Default = 0')
#------------------------------------------------------------------------------
#Initialize environment
env = Environment(ASSEMBLY=GetOption('assembly'),
                          FQDIR=GetOption('fastq_dir'),
                          OUTDIR=GetOption('outdir'),
                          SIDS=GetOption('sids'),
                          NHEADER=GetOption('nheader'),
                          INTDEPTH=GetOption('nointdepth'),
                          MARKDUP=GetOption('markdup'))
#=============================================================================
###
### Builders
###

#OPTION DICTIONARIES
#Options for bwa mem, samtools sort, depthfile
bwa_mem_opts = {'-t':GetOption('align_thread')}
samtools_sort_opts = {'-@':GetOption('samsort_thread'), '-m':GetOption('samsort_mem'), '-T':GetOption('tmpdir')}
depthfile_opts_bin = {'--noIntraDepthVariance':''}
#------------------------------------------------------------------------------
#BWA index builder, add index targets as default targets
bwa_index_builder = Builder(action = 'bwa index $SOURCE')

#Option strings for bwa mem, samtools sort
bwa_optstring = optstring_join(bwa_mem_opts)
samtools_sort_optstring = optstring_join(samtools_sort_opts)

#Builder for pipe: read mapping, SAM reduction, SAM to BAM
#FASTQ files are interleaved
#bwa_samtools_intl_action = 'bwa mem %s ${SOURCES[0]} -p ${SOURCES[1]} | samtools view -hS -F4 - | tee ${TARGETS[0]} | samtools view -huS - | samtools sort %s -o - > ${TARGETS[1]}' % (bwa_optstring, samtools_sort_optstring)
bwa_samtools_intl_action = 'bwa mem %s ${SOURCES[0]} -p ${SOURCES[1]} | samtools view -huS -F4 - | samtools sort %s -o - > $TARGET' % (bwa_optstring, samtools_sort_optstring)
bwa_samtools_intl_builder = Builder(action = bwa_samtools_intl_action)

#Builder for pipe: read mapping, SAM reduction, SAM to BAM
#FASTQ files are separate R1 and R2
#bwa_samtools_r1r2_action = 'bwa mem %s ${SOURCES[0]} ${SOURCES[1]} ${SOURCES[2]} | samtools view -hS -F4 - | tee ${TARGETS[0]} | samtools view -huS - | samtools sort %s -o - > ${TARGETS[1]}' % (bwa_optstring, samtools_sort_optstring)
bwa_samtools_r1r2_action = 'bwa mem %s ${SOURCES[0]} ${SOURCES[1]} ${SOURCES[2]} | samtools view -huS -F4 - | samtools sort %s -o - > $TARGET' % (bwa_optstring, samtools_sort_optstring)
bwa_samtools_r1r2_builder = Builder(action = bwa_samtools_r1r2_action)

bwa_samtools_intl_markdup = """bwa mem %s ${SOURCES[0]} -p ${SOURCES[1]} | \
    samtools view -hu -F4 - | \
    samtools sort %s -n -O bam - | \
    samtools fixmate -mr -O bam - - | \
    samtools sort %s -O bam - | \
    samtools markdup -r - $TARGET""" % (bwa_optstring, samtools_sort_optstring, samtools_sort_optstring)
bwa_samtools_intl_markdup_builder = Builder(action = bwa_samtools_intl_markdup)

bwa_samtools_r1r2_markdup = """bwa mem %s ${SOURCES[0]} ${SOURCES[1]} ${SOURCES[2]} | \
    samtools view -hu -F4 - | \
    samtools sort %s -n -O bam - | \
    samtools fixmate -mr -O bam - - | \
    samtools sort %s -O bam - | \
    samtools markdup -r - $TARGET""" % (bwa_optstring, samtools_sort_optstring, samtools_sort_optstring)
bwa_samtools_r1r2_markdup_builder = Builder(action = bwa_samtools_r1r2_markdup)
#------------------------------------------------------------------------------
#Builder for depthfile creation
if GetOption('nointdepth'):
    depthfile_bin_action = 'src/jgi_summarize_bam_contig_depths --outputDepth $TARGET $SOURCES'
    depthfile_bin_builder = Builder(action = depthfile_bin_action)
else:
    depthfile_bin_action = 'src/jgi_summarize_bam_contig_depths --outputDepth $TARGET $SOURCES %s' % optstring_join(depthfile_opts_bin)
    depthfile_bin_builder = Builder(action = depthfile_bin_action)
#------------------------------------------------------------------------------
builders = {'BWA_Samtools_Intl':bwa_samtools_intl_builder,
'BWA_Samtools_R1R2':bwa_samtools_r1r2_builder,
'BWA_Samtools_Markdup_Intl':bwa_samtools_intl_markdup_builder,
'BWA_Samtools_Markdup_R1R2':bwa_samtools_r1r2_markdup_builder,
'Depthfile_Bin':depthfile_bin_builder,
'BWA_index':bwa_index_builder}
env.Append(BUILDERS = builders)
#=============================================================================
#SConscript
if env['ASSEMBLY']:
    build_tmp = os.path.splitext(os.path.basename(env['ASSEMBLY']))[0] + '_build'
    SConscript(['SConscript'], exports='env', variant_dir=build_tmp, duplicate=0)
#------------------------------------------------------------------------------
#If --rmbuild=1, remove the build targets in the temporary directory
if GetOption('rmbuild'):
    atexit.register(remove_build_targets, tmpdir = build_tmp)

else:
    pass

Export('env')
