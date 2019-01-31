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


#env.SConscript('mapping/SConscript', variant_dir=env['OUTDIR'], exports='env')
#env.SConscript('network/SConscript', variant_dir=env['OUTDIR'], exports='env')
env.SConscript('mapping/SConscript', exports='env')
env.SConscript('network/SConscript', exports='env')

#env.SConscript(['mapping/SConscript', 'network/SConscript'], exports='env')

Export('env')
