Import('env')
import glob
import os
import gzip
import re
from itertools import islice
from warnings import warn
from fastq_pair import *
import logging

'''
2020 Ian Rambo
Thirteen... that's a mighty unlucky number... for somebody!

v.1.0.0.2
'''
version = 'v.1.0.0.2'
#Error message to display if trying to use markdup builder with single-end reads
single_markdup_err_msg = 'ERROR: fixmate and markdup build not supported for single-end reads in version %s. Retry with --markdup=0' % version

#Set up the logger for this build
logging.basicConfig(format='%(asctime)s - %(message)s',
    datefmt='%d-%b-%y %H:%M:%S', level=logging.DEBUG,
    filename=env['LOGFILE'])
#Write info about wrapper version to LOGFILE
version_message = 'this is scons genome mapping wrapper %s' % version
logging.info(version_message)

print('writing logfile to %s' % env['LOGFILE'])
bwa_index_targets = [env['ASSEMBLY'] + ext for ext in ['.bwt','.pac','.ann','.amb','.sa']]
Default(bwa_index_targets)
#Index the assembly
env.BWA_index(bwa_index_targets, env['ASSEMBLY'])

#get the basename of the reference assembly and use as a build ID
assembly_id = get_basename(env['ASSEMBLY'])
#------------------------------------------------------------------------------
#Generate list of input FASTQ files using sample IDs
fastq_list = source_list_generator(env['SIDS'], env['FQDIR'])
#------------------------------------------------------------------------------
mapping_targets = list()
fastq_dict = find_fastq_pairs(fastq_list, nheader = env['NHEADER'])

bam_extension = '.reduced.sorted.bam'
if env['MARKDUP']:
    markdup_msg = 'Mark duplicates in this build'
    logging.info(markdup_msg)
    bam_extension = '.reduced.sorted.markdup.bam'
else:
    pass

#Loop through the FASTQ files and create the mapping TARGETS
#The mapping method will be decided based on key-value pairs in the fastq_dict
for key in fastq_dict:
    maptarg = []
    if os.path.isfile(fastq_dict[key]['R1']) and not fastq_dict[key]['R2'] == 'single':
        if 'R2' in fastq_dict[key].keys() and fastq_dict[key]['R2'] == 'interleaved':
            #Interleaved mapping
            maptarg = [assembly_id + '____' + get_basename(fastq_dict[key]['R1']) + bam_extension]
            if env['MARKDUP']:
                logging.info('fixmates and mark duplicates in interleaved file %s' % fastq_dict[key]['R1'])
                env.BWA_Samtools_Markdup_Intl(maptarg, [env['ASSEMBLY'], fastq_dict[key]['R1']])
            else:
                logging.info('bwa > samtools sort for interleaved file %s' % fastq_dict[key]['R1'])
                env.BWA_Samtools_Intl(maptarg, [env['ASSEMBLY'], fastq_dict[key]['R1']])

        elif os.path.isfile(fastq_dict[key]['R2']):
            #R1 R2 mapping
            maptarg = [assembly_id + '____' + get_basename(fastq_dict[key]['R1']) + bam_extension]
            if env['MARKDUP']:
                logging.info('fixmates and mark duplicates for R1-R2 %s %s' % (fastq_dict[key]['R1'], fastq_dict[key]['R2']))
                env.BWA_Samtools_Markdup_R1R2(maptarg, [env['ASSEMBLY'], fastq_dict[key]['R1'], fastq_dict[key]['R2']])
            else:
                env.BWA_Samtools_R1R2(maptarg, [env['ASSEMBLY'], fastq_dict[key]['R1'], fastq_dict[key]['R2']])
        else:
            pass

        mapping_targets.extend(maptarg)
    elif fastq_dict[key]['R1'] == 'single' or fastq_dict[key]['R2'] == 'single':
        if env['MARKDUP']:
            print(single_markdup_err_msg)
            logging.error(single_markdup_err_msg)
            break
        else:
            single_path = [fastq_dict[key][v] for v in fastq_dict[key].keys() if os.path.isfile(fastq_dict[key][v])][0]
            maptarg = [assembly_id + '____' + get_basename(single_path) + bam_extension]
            env.BWA_Samtools_Single(maptarg, [env['ASSEMBLY'], single_path])
            mapping_targets.extend(maptarg)

    elif 'error' in fastq_dict[key].keys():
        error_path = [fastq_dict[key][v] for v in fastq_dict[key].keys() if os.path.isfile(fastq_dict[key][v])][0]
        logging.warning('FASTQ with errors %s' % error_path)
        break
    else:
        logging.error('WARNING: no interleaved, R1-R2 pair, or single-end FASTQ found. Please inspect your data.')
#------------------------------------------------------------------------------
#Depth file
depthfile_bin_target = assembly_id + '_cov'

depthfile_sources = [m for m in mapping_targets if re.match(r'.*?\.bam', m)]

Default(env.Install(env['OUTDIR'], depthfile_bin_target))
Default(env.Install(env['OUTDIR'], depthfile_sources))
env.Depthfile_Bin(depthfile_bin_target, depthfile_sources)
#------------------------------------------------------------------------------
