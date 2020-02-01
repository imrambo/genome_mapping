#!/usr/bin/env python3
import glob
import os
import gzip
import re
from itertools import islice
import logging
import magic

'''
2020 Ian Rambo
Thirteen... that's a mighty unlucky number... for somebody!
'''
#------------------------------------------------------------------------------
def source_list_generator(id_string, source_dir):
    """
    Generate list of source files using sample IDs
    """
    source_list = list()
    if ',' in id_string:
        id_list = id_string.split(',')
        for i in id_list:
            idGlob = os.path.join(os.path.abspath(source_dir), '%s*' % i)
            if idGlob:
                source_list.extend(list(glob.glob(idGlob, recursive = False)))
            else:
                logging.warning('file(s) for identifier %s not found' % i)
    else:
        idGlob = os.path.join(os.path.abspath(source_dir), '%s*' % id_string)
        if idGlob:
            source_list.extend(list(glob.glob(idGlob, recursive = False)))
        else:
            logging.warning('file(s) for identifier %s not found' % i)

    return source_list
#------------------------------------------------------------------------------
def get_basename(file_path):
    """
    Get file basename, excluding multiple patterns.
    """
    basename = os.path.basename(file_path)
    #Remove two patterns, e.g. foo.tar.gz becomes foo
    if re.match(r'^.*?\.[a-z]+\.[a-z]+$', basename):
        basename = re.findall(r'^(.*?)\.[a-z]+\.[a-z]+$', basename)[0]
    else:
        basename = os.path.splitext(basename)[0]
    return basename
#------------------------------------------------------------------------------
def is_gzipped(file_path):
    """
    Test if a file is gzipped.
    """
    is_gzip = False
    if magic.from_file(file_path).startswith('gzip compressed data'):
        is_gzip = True
        return is_gzip
    else:
        return is_gzip
#------------------------------------------------------------------------------
def find_fastq_pairs(fastq_list, nheader='ALL', exclude = False):
    '''
    Identify Illumina FASTQ reads as interleaved or non-interleaved by
    comparing header order (nheader is the number of headers to use).
    Default is to read all headers.
    Match paired-end FASTQ reads using header information.
    Requires Casava 1.8+
    '''
    if isinstance(nheader, int) and nheader < 4:
        logging.warning('slice more than four headers for comparison')
    fastq_dict = {}
    interleaved = False
    #Multiply number of headers by 4 to fetch FASTQ entries
    if isinstance(nheader, int):
        nrows = nheader * 4
    for fastq in fastq_list:
        print('testing %s...' % fastq)
        head_list = []
        #Test if file is gzip compressed
        if is_gzipped(fastq):
            print('file %s is gzip compressed' % fastq)
            with gzip.open(fastq, 'r') as fq:
                if isinstance(nheader, int):
                    head_list = [l.decode('utf-8').strip().split() for l in islice(fq, nrows) if re.match(r'^\@.*?\:\d+\:.*?\:\d+\:\d+\:\d+\:\d+\s+\d\:.*?\:[ACTGN]+', l.decode('utf-8').strip())]
                elif isintance(nheader, str) and nheader == 'ALL':
                    head_list = [l.decode('utf-8').strip().split() for l in fq if re.match(r'^\@.*?\:\d+\:.*?\:\d+\:\d+\:\d+\:\d+\s+\d\:.*?\:[ACTGN]+', l.decode('utf-8').strip())]
                else:
                    logging.error('fastq_pair: nheader value must be an integer, or string "ALL".')
                print(head_list)
        else:
            with open(fastq, 'r') as fq:
                if isinstance(nheader, int):
                    head_list = [l.strip().split() for l in islice(fq, nrows) if re.match(r'^\@.*?\:\d+\:.*?\:\d+\:\d+\:\d+\:\d+\s+\d\:.*?\:[ACTGN]+', l)]
                elif isintance(nheader, str) and nheader == 'ALL':
                    head_list = [l.strip().split() for l in fq if re.match(r'^\@.*?\:\d+\:.*?\:\d+\:\d+\:\d+\:\d+\s+\d\:.*?\:[ACTGN]+', l)]
                else:
                    logging.error('fastq_pair: nheader value must be an integer, or string "ALL".')
        int_test_list = []
        #error_pairs = []

        #Loop through the read headers to see if they are correctly paired
        for n in range(1, len(head_list), 2):
            if head_list[n-1][0] == head_list[n][0] and head_list[n-1][1].startswith('1') and head_list[n][1].startswith('2'):
                interleave_pair = True
            else:
                interleave_pair = False
            int_test_list.append(interleave_pair)
            # if interleave_pair == False:
            #     ep = (head_list[n-1], head_list[n], n-1, n)
            #     error_pairs.append(ep)

        if all(int_test_list):
            interleaved = True
            if isinstance(nheader, int):
                logging.info('%s: interleaved FASTQ based on %d headers' % (os.path.basename(fastq), nheader))
            elif isinstance(nheader, str) and nheader == 'ALL':
                logging.info('%s: interleaved FASTQ based on ALL (%d) headers' % (os.path.basename(fastq), len(head_list)))
        #elif all([head_list[n-1][0] != head_list[n][0] for n in range(1, len(head_list), 2)]):
            #interleaved = False
        else:
            interleaved = False
            #print(error_pairs)

        if interleaved:
            #If the interleaved FASTQ is not in the dicitonary, add it
            if not head_list[0][0] in fastq_dict:
                fastq_dict[head_list[0][0]] = {}
                fastq_dict[head_list[0][0]]['R1'] = fastq
                fastq_dict[head_list[0][0]]['R2'] = 'interleaved'
            else:
                pass
        else:
            #Test to see if the reads are single-end
            single = False
            if all([head_list[n-1][0] == head_list[n][0] and head_list[n-1][1][0] == head_list[n][1][0] for n in range(1, len(head_list), 2)]):
                single = True

            if single and not head_list[0][0] in fastq_dict:
                if head_list[0][1].startswith('1'):
                    if isinstance(nheader, int):
                        logging.info('%s: single-end FASTQ R1 based on %d headers' % (os.path.basename(fastq), nheader))
                    elif isintance(nheader, str) and nheader == 'ALL':
                        logging.info('%s: single-end FASTQ R1 based on ALL (%d) headers' % (os.path.basename(fastq), len(head_list)))
                    fastq_dict[head_list[0][0]]['R1'] = fastq
                    fastq_dict[head_list[0][0]]['R2'] = 'single'
                elif head_list[0][1].startswith('2'):
                    if isinstance(nheader, int):
                        logging.info('%s: single-end FASTQ R1 based on %d headers' % (os.path.basename(fastq), nheader))
                    elif isintance(nheader, str) and nheader == 'ALL':
                        logging.info('%s: single-end FASTQ R2 based on ALL (%d) headers' % (os.path.basename(fastq), len(head_list)))
                    fastq_dict[head_list[0][0]]['R1'] = 'single'
                    fastq_dict[head_list[0][0]]['R2'] = fastq

                else:
                    logging.warning('please check if FASTQ header is in format Casava 1.8+')
            else:
                fastq_dict[head_list[0][0]] = {}
                if head_list[0][1].startswith('1'):
                    if isinstance(nheader, int):
                        logging.warning('%s: non-interleaved FASTQ R1 : potential interleave error based on %d headers' % (os.path.basename(fastq), nheader))
                    elif isintance(nheader, str) and nheader == 'ALL':
                        logging.warning('%s: non-interleaved FASTQ R1 : potential interleave error based on ALL (%d) headers' % (os.path.basename(fastq), len(head_list)))
                    else:
                        pass
                    fastq_dict[head_list[0][0]]['R1'] = fastq
                    fastq_dict[head_list[0][0]]['R2'] = 'error'
                elif head_list[0][1].startswith('2'):
                    if isinstance(nheader, int):
                        logging.warning('%s: non-interleaved FASTQ R2 : potential interleave error based on %d headers' % (os.path.basename(fastq), nheader))
                    elif isintance(nheader, str) and nheader == 'ALL':
                        logging.warning('%s: non-interleaved FASTQ R1 : potential interleave error based on ALL (%d) headers' % (os.path.basename(fastq), len(head_list)))
                    else:
                        pass
                    fastq_dict[head_list[0][0]]['R1'] = 'error'
                    fastq_dict[head_list[0][0]]['R2'] = fastq
                else:
                    logging.warning('please check if FASTQ header is in format Casava 1.8+')

    for key in fastq_dict:
        if 'error' in fastq_dict[key].keys():
            solo = [fastq_dict[key][v] for v in fastq_dict[key].keys() if os.path.isfile(fastq_dict[key][v])][0]
            solo_message = "\n#=====\nNo pair for non-interleaved reads %s found, file may contain errors.\nPlease include read pair and check if the sequences are correctly interleaved or single-end.\n#=====\n" % solo
            logging.warning(solo_message)
            if exclude:
                #Remove entries from dictionary with errors
                del fastq_dict[key]
                logging.info('exlude == True, removed non-paired FASTQ %s from dictionary' % solo)

            else:
                pass

    return fastq_dict
#=============================================================================
# def find_fastq_pairs(fastq_list, nslice, exclude = False):
#     '''
#     Identify FASTQ reads as interleaved or non-interleaved.
#     Match paired-end FASTQ reads using header information.
#     '''
#     if nslice < 4:
#         logging.warning('slice more than four headers')
#     fastq_dict = {}
#     interleaved = False
#     #Multiply number of headers by 4 to fetch FASTQ entries
#     nslice = nslice * 4
#     for fastq in fastq_list:
#         head_list = []
#         #Test if file is gzip compressed
#         if magic.from_file(fastq).startswith('gzip compressed data') and fastq.endswith('.gz'):
#             with gzip.open(fastq, 'r') as fq:
#                 head_list = [l.decode('utf-8').split() for l in islice(fq, nslice) if re.match(r'^\@.*?\:\d+\:.*?\:\d+\:\d+\:\d+\:\d+\s+\d\:.*?\:[ACTGN]+', l.decode('utf-8'))]
#         else:
#             with open(fastq, 'r')as fq:
#                 #head_list = [l.split() for l in islice(fq, nslice) if l.startswith('@')]
#                 head_list = [l.split() for l in islice(fq, nslice) if re.match(r'^\@.*?\:\d+\:.*?\:\d+\:\d+\:\d+\:\d+\s+\d\:.*?\:[ACTGN]+', l)]
#         #int_test_list = [head_list[n-1][0] == head_list[n][0] and head_list[n-1][1].startswith('1') and head_list[n][1].startswith('2') for n in range(1, len(head_list), 2)]
#         int_test_list = []
#         error_pairs = []
#         for n in range(1, len(head_list), 2):
#             if head_list[n-1][0] == head_list[n][0] and head_list[n-1][1].startswith('1') and head_list[n][1].startswith('2'):
#                 interleave_pair = True
#             else:
#                 interleave_pair = False
#             int_test_list.append(interleave_pair)
#             if interleave_pair == False:
#                 ep = (head_list[n-1], head_list[n], n-1, n)
#                 error_pairs.append(ep)
#
#         if all(int_test_list):
#             interleaved = True
#             print('%s: interleaved FASTQ' % os.path.basename(fastq))
#         #elif all([head_list[n-1][0] != head_list[n][0] for n in range(1, len(head_list), 2)]):
#             #interleaved = False
#         else:
#             interleaved = False
#             print(error_pairs)
#
#         if interleaved:
#             if not head_list[0][0] in fastq_dict:
#                 fastq_dict[head_list[0][0]] = {}
#                 fastq_dict[head_list[0][0]]['R1'] = fastq
#                 fastq_dict[head_list[0][0]]['R2'] = 'interleaved'
#             else:
#                 pass
#         else:
#             if head_list[0][0] in fastq_dict:
#                 if head_list[0][1].startswith('1'):
#                     print('%s: non-interleaved FASTQ R1' % os.path.basename(fastq))
#                     fastq_dict[head_list[0][0]]['R1'] = fastq
#                 elif head_list[0][1].startswith('2'):
#                     print('%s: non-interleaved FASTQ R2' % os.path.basename(fastq))
#                     fastq_dict[head_list[0][0]]['R2'] = fastq
#                 else:
#                     logging.warning('please check if FASTQ header is in format Casava 1.8+')
#             else:
#                 fastq_dict[head_list[0][0]] = {}
#                 if head_list[0][1].startswith('1'):
#                     print('%s: non-interleaved FASTQ R1' % os.path.basename(fastq))
#                     fastq_dict[head_list[0][0]]['R1'] = fastq
#                 elif head_list[0][1].startswith('2'):
#                     print('%s: non-interleaved FASTQ R2' % os.path.basename(fastq))
#                     fastq_dict[head_list[0][0]]['R2'] = fastq
#                 else:
#                     logging.warning('please check if FASTQ header is in format Casava 1.8+')
#
#     for key in fastq_dict:
#         if not 'R2' in fastq_dict[key].keys() or not 'R1' in fastq_dict[key].keys():
#             solo = [fastq_dict[key][v] for v in fastq_dict[key].keys()][0]
#             warning_message = "\n#=====\nNo pair for non-interleaved reads %s found, please include read pair or check if the sequences are correctly interleaved\n#=====\n" % solo
#             logging.warning(warning_message)
#             if exclude == True:
#                 #Remove entries from dictionary with errors
#                 del fastq_dict[key]
#                 logging.info('exlude == True, removed faulty FASTQ %s from dictionary' % solo)
#
#             else:
#                 pass
#
#     return fastq_dict
