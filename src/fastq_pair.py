#!/usr/bin/env python3
import glob
import os
import gzip
import re
from itertools import islice
import logging

'''
2019 Ian Rambo
Thirteen... that's a mighty unlucky number... for somebody!
'''
#------------------------------------------------------------------------------
def find_fastq_pairs(fastq_list, nslice = 800, exclude = False):
    '''
    Match paired-end FASTQ reads using header information.
    '''
    if nslice < 4:
        logging.warning('slice more than four headers')
    fastq_dict = {}
    interleaved = False
    #Multiply number of headers by 4 to fetch FASTQ entries
    nslice = nslice * 4
    for fastq in fastq_list:
        head_list = []
        if gzipped: ###CHANGE
            with gzip.open(fastq, 'r') as fq:
                head_list = [l.decode('utf-8').split() for l in islice(fq, nslice) if l.decode('utf-8').startswith('@')]
        if not gzipped: ###CHANGE
            with open(fastq, 'r' as fq):
                head_list = [l.split() for l in islice(fq, nslice) if l.startswith('@')]
        int_test_list = [head_list[n-1][0] == head_list[n][0] and head_list[n-1][1].startswith('1') and head_list[n][1].startswith('2') for n in range(1, len(head_list), 2)]
        if all(int_test_list):
            interleaved = True
        #if all([head_list[n-1][0] != head_list[n][0] for n in range(1, len(head_list), 2)]):
        else:
            interleaved = False
            print('file %s not interleaved: %d of %d read pairs were interleaved' % (fastq, sum(int_test_list), len(int_test_list)))

        if interleaved:
            print('%s: interleaved FASTQ' % os.path.basename(fastq))
            if not head_list[0][0] in fastq_dict:
                fastq_dict[head_list[0][0]] = {}
                fastq_dict[head_list[0][0]]['R1'] = fastq
                fastq_dict[head_list[0][0]]['R2'] = 'interleaved'
            else:
                pass
        else:
            if head_list[0][0] in fastq_dict:
                if head_list[0][1].startswith('1'):
                    print('%s: non-interleaved FASTQ R1' % os.path.basename(fastq))
                    fastq_dict[head_list[0][0]]['R1'] = fastq
                elif head_list[0][1].startswith('2'):
                    print('%s: non-interleaved FASTQ R2' % os.path.basename(fastq))
                    fastq_dict[head_list[0][0]]['R2'] = fastq
                else:
                    print('please check if FASTQ header is in format Casava 1.8+')
            else:
                fastq_dict[head_list[0][0]] = {}
                if head_list[0][1].startswith('1'):
                    print('%s: non-interleaved FASTQ R1' % os.path.basename(fastq))
                    fastq_dict[head_list[0][0]]['R1'] = fastq
                elif head_list[0][1].startswith('2'):
                    print('%s: non-interleaved FASTQ R2' % os.path.basename(fastq))
                    fastq_dict[head_list[0][0]]['R2'] = fastq
                else:
                    print('please check if FASTQ header is in format Casava 1.8+')

    for key in fastq_dict:
        if not 'R2' in fastq_dict[key].keys() or not 'R1' in fastq_dict[key].keys():
            solo = [fastq_dict[key][v] for v in fastq_dict[key].keys()][0]
            warning_message = "No pair for reads %s found, please include read pair or check if the sequences are correctly interleaved" % solo
            logging.warning(warning_message)
            if exclude == True:
                #Remove entries from dictionary with errors
                del fastq_dict[key]
                logging.info('exlude == True, removed faulty FASTQ %s from dictionary' % solo)
                pass

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
