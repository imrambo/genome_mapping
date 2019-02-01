#! /usr/bin/env python3
import glob
import os
import gzip
import re
from itertools import islice
from warnings import warn
import argparse
import subprocess
"""
2019 Ian Rambo
"""
#==============================================================================
def optstring_join(optdict):
    """
    Join a dictionary of command line options into a single string.
    """
    optstring = ' '.join([str(param) + ' ' + str(val) for param, val in optdict.items()])
    return optstring
#------------------------------------------------------------------------------
def find_fastq_pairs(fastq_list, nslice = 800):
    """
    Determine if FASTQ files are paired or interleaved, and match them
    accordingly based on header information.
    """
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
def bwa_index(source):
    """
    Index a FASTA file for BWA mapping.
    """
    bwa_index_action = 'bwa index %s' % source
    p = subprocess.Popen(bwa_index_action, shell = True)
    return None
#------------------------------------------------------------------------------
def bwa_sam_intl(targets, sources, bopts, sopts):
    """
    Pipe for bwa-mem and samtools producing reduced SAM and BAM using
    interleaved FASTQ files.
    """
    bwa_samtools_intl_action = 'bwa mem %s %s -p %s | samtools view -hS -F4 - | tee %s | samtools view -huS - | samtools sort %s - -o %s' % (bopts, sources[0], sources[1], targets[0], sopts, targets[1])
    p = subprocess.Popen(bwa_samtools_intl_action, shell = True)
    return None
#------------------------------------------------------------------------------
def bwa_sam_r1r2(targets, sources, bopts, sopts):
    """
    Pipe for bwa-mem and samtools producing reduced SAM and BAM using
    non-interleaved FASTQ files.
    """
    bwa_samtools_r1r2_action = 'bwa mem %s %s %s %s | samtools view -hS -F4 - | tee %s | samtools view -huS - | samtools sort %s - -o %s' % (bopts, ' '.join(sources), targets[0], sopts, targets[1])
    p = subprocess.Popen(bwa_samtools_r1r2_action, shell = True)
    return None
#------------------------------------------------------------------------------
def depth_file(target, sources, opts):
    """
    Create a depth file from one or more BAM files.
    """
    depthfile_action = 'src/jgi_summarize_bam_contig_depth --outputDepth %s %s %s' % (target, ' '.join(sources), opts)
    p = subprocess.Popen(depthfile_action, shell = True)
    return None
#------------------------------------------------------------------------------
def network_file(target, source):
    """
    Create a network file for MMGenome1. Requires a SAM file as input.
    """
    network_action = 'perl src/network.pl -i %s -o %s' % (source, target)
    p = subprocess.Popen(network_action, shell = True)
    return None
#------------------------------------------------------------------------------
def source_list_generator(id_string, source_dir, extension):
    """
    Generate list of source files using sample IDs
    """
    source_list = list()
    if ',' in id_string:
        id_list = id_string.split(',')
        for i in id_list:
            idGlob = os.path.join(os.path.abspath(source_dir), '%s*%s' % (i, extension))
            source_list.extend(list(glob.iglob(idGlob, recursive = False)))
    else:
        idGlob = os.path.join(os.path.abspath(source_dir), '%s*%s' % (id_string, extension))
        source_list.extend(list(glob.iglob(idGlob, recursive = False)))

    return source_list
