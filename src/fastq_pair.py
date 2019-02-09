#!/usr/bin/env python3
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
#------------------------------------------------------------------------------
def find_fastq_pairs(fastq_list, nslice = 800):
    '''
    Match paired-end FASTQ reads using header information.
    '''
    if nslice % 4:
        warn('--nslice is not a multiple of 4')
    if nslice < 16:
        warn('slice more than four headers')
    fastq_dict = {}
    interleaved = False
    #Multiply number of headers by 4 to fetch FASTQ entries
    nslice = nslice * 4
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
        #if not fastq_dict[key]['R2'] or not fastq_dict[key]['R1']:
        if not 'R2' in fastq_dict or not 'R1' in fastq_dict:
            solo = [fastq_dict[key][v] for v in fastq_dict[key]][0]
            warn("""one is the loneliest number that you'll ever doooooo...
            so find a pair for %s if they are paired-end reads...""" % solo)

    return fastq_dict
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
