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

version 1.0.0.2
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
    Requires Casava 1.8+ format.
    '''
    if isinstance(nheader, int) and nheader > 0 and nheader < 4:
        logging.warning('slice more than four headers for comparison')
    fastq_dict = {}
    interleaved = False
    #Multiply number of headers by 4 to fetch FASTQ entries
    if isinstance(nheader, int):
        nrows = nheader * 4
    for fastq in fastq_list:
        head_list = []
        #Inspect the FASTQ headers to see if they match the pattern for Casava 1.8+
        if is_gzipped(fastq):
            with gzip.open(fastq, 'r') as fq:
                if isinstance(nheader, int):
                    head_list = [l.decode('utf-8').strip().split() for l in islice(fq, nrows) if re.match(r'^\@.*?\:\d+\:.*?\:\d+\:\d+\:\d+\:\d+\s+\d\:.*?\:[ACTGN]+', l.decode('utf-8').strip())]
                elif isintance(nheader, str) and nheader == 'ALL':
                    head_list = [l.decode('utf-8').strip().split() for l in fq if re.match(r'^\@.*?\:\d+\:.*?\:\d+\:\d+\:\d+\:\d+\s+\d\:.*?\:[ACTGN]+', l.decode('utf-8').strip())]
                else:
                    logging.error('fastq_pair: nheader value must be an integer, or string "ALL".')
        else:
            with open(fastq, 'r') as fq:
                if isinstance(nheader, int):
                    head_list = [l.strip().split() for l in islice(fq, nrows) if re.match(r'^\@.*?\:\d+\:.*?\:\d+\:\d+\:\d+\:\d+\s+\d\:.*?\:[ACTGN]+', l)]
                elif isintance(nheader, str) and nheader == 'ALL':
                    head_list = [l.strip().split() for l in fq if re.match(r'^\@.*?\:\d+\:.*?\:\d+\:\d+\:\d+\:\d+\s+\d\:.*?\:[ACTGN]+', l)]
                else:
                    logging.error('fastq_pair: nheader value must be an integer, or string "ALL".')
        int_test_list = []

        #Loop through the FASTQ read headers to see if they are correctly paired
        for n in range(1, len(head_list), 2):
            if head_list[n-1][0] == head_list[n][0] and head_list[n-1][1].startswith('1') and head_list[n][1].startswith('2'):
                interleave_pair = True
            else:
                interleave_pair = False
            int_test_list.append(interleave_pair)

        if all(int_test_list):
            interleaved = True
            if isinstance(nheader, int):
                logging.info('%s: interleaved FASTQ based on %d headers' % (os.path.basename(fastq), nheader))
            elif isinstance(nheader, str) and nheader == 'ALL':
                logging.info('%s: interleaved FASTQ based on ALL (%d) headers' % (os.path.basename(fastq), len(head_list)))

        else:
            interleaved = False

        #Sample Illumina identifier
        fastq_id_tag = head_list[0][0]

        if interleaved:
            #If the interleaved FASTQ is not in the dictionary, add it
            if not fastq_id_tag in fastq_dict.keys():
                fastq_dict[fastq_id_tag] = {}
                fastq_dict[fastq_id_tag]['R1'] = fastq
                fastq_dict[fastq_id_tag]['R2'] = 'interleaved'
            else:
                if os.path.isfile(fastq_dict[fastq_id_tag]['R1']) or fastq_dict[fastq_id_tag]['R2'] != 'interleaved':
                    ident_exist = [fastq_dict[key][v] for key in fastq_dict.keys() for v in fastq_dict[key].keys() if os.path.isfile(fastq_dict[key][v])]
                    intl_warn_msg = '\nWARNING: same identifier %s present in both %s and interleaved FASTQ %s\n' % (fastq_id_tag, ' '.join(ident_exist), fastq)
                    print(intl_warn_msg)
                    logging.warning(intl_warn_msg)
                    use_int_msg = '\nInterleaved FASTQ %s will be used instead of single-end FASTQ %s\n' % (fastq, ' '.join(ident_exist))
                    print(use_int_msg)
                    logging.info(use_int_msg)
                    fastq_dict[fastq_id_tag]['R1'] = fastq
                    fastq_dict[fastq_id_tag]['R2'] = 'interleaved'
        else:
            #Test to see if the reads are single-end
            single = False
            single_test = all([head_list[n-1][1][0] == head_list[n][1][0] for n in range(1, len(head_list), 2)])
            if single_test:
                single = True

            if single:
                if not fastq_id_tag in fastq_dict:
                    fastq_dict[fastq_id_tag] = {}
                    if head_list[0][1].startswith('1'):
                        if isinstance(nheader, int):
                            logging.info('%s: single-end FASTQ R1 based on %d headers' % (os.path.basename(fastq), nheader))
                        elif isintance(nheader, str) and nheader == 'ALL':
                            logging.info('%s: single-end FASTQ R1 based on ALL (%d) headers' % (os.path.basename(fastq), len(head_list)))
                        else:
                            pass
                        fastq_dict[fastq_id_tag]['R1'] = fastq
                        fastq_dict[fastq_id_tag]['R2'] = 'single'
                    elif head_list[0][1].startswith('2'):
                        if isinstance(nheader, int):
                            logging.info('%s: single-end FASTQ R2 based on %d headers' % (os.path.basename(fastq), nheader))
                        elif isintance(nheader, str) and nheader == 'ALL':
                            logging.info('%s: single-end FASTQ R2 based on ALL (%d) headers' % (os.path.basename(fastq), len(head_list)))
                        fastq_dict[fastq_id_tag]['R1'] = 'single'
                        fastq_dict[fastq_id_tag]['R2'] = fastq

                    else:
                        logging.warning('please check if FASTQ header is in format Casava 1.8+')
                else:
                    if os.path.isfile(fastq_dict[fastq_id_tag]['R1']) and os.path.isfile(fastq_dict[fastq_id_tag]['R2']):
                        single_warn_msg = 'More than two FASTQ files found with the same identifier %s for single-end reads: %s  %s  %s' % (fastq_id_tag, fastq_dict[fastq_id_tag]['R1'], fastq_dict[fastq_id_tag]['R2'], fastq)
                        logging.warning(single_warn_msg)
                        break
                    if fastq_dict[fastq_id_tag]['R2'] == 'interleaved':
                        logging.info('Interleaved FASTQ %s with identifer %s already in dictionary, will ignore single-end reads %s' % (fastq_dict[fastq_id_tag]['R1'], fastq_id_tag, fastq))
                        break
                    #If matching read pair is found, add it to the dictionary and replace 'single' type
                    if fastq_dict[fastq_id_tag]['R1'] == 'single' and head_list[0][1].startswith('1'):
                        fastq_dict[fastq_id_tag]['R1'] = fastq
                    elif fastq_dict[fastq_id_tag]['R2'] == 'single' and head_list[0][1].startswith('2'):
                        fastq_dict[fastq_id_tag]['R2'] = fastq
                    else:
                        pass


            else:
                if not fastq_id_tag in fastq_dict.keys():
                    fastq_dict[fastq_id_tag] = {}
                if head_list[0][1].startswith('1'):
                    if isinstance(nheader, int):
                        logging.warning('%s: non-interleaved FASTQ R1 : potential interleave error based on %d headers' % (os.path.basename(fastq), nheader))
                    elif isintance(nheader, str) and nheader == 'ALL':
                        logging.warning('%s: non-interleaved FASTQ R1 : potential interleave error based on ALL (%d) headers' % (os.path.basename(fastq), len(head_list)))
                    else:
                        pass
                    fastq_dict[fastq_id_tag]['R1'] = fastq
                    fastq_dict[fastq_id_tag]['R2'] = 'error'
                elif head_list[0][1].startswith('2'):
                    if isinstance(nheader, int):
                        logging.warning('%s: non-interleaved FASTQ R2 : potential interleave error based on %d headers' % (os.path.basename(fastq), nheader))
                    elif isintance(nheader, str) and nheader == 'ALL':
                        logging.warning('%s: non-interleaved FASTQ R2 : potential interleave error based on ALL (%d) headers' % (os.path.basename(fastq), len(head_list)))
                    else:
                        pass
                    fastq_dict[fastq_id_tag]['R1'] = 'error'
                    fastq_dict[fastq_id_tag]['R2'] = fastq
                else:
                    logging.warning('please check if FASTQ header is in format Casava 1.8+')
                    break

    for key in fastq_dict:
        if 'error' in fastq_dict[key].keys():
            solo = [fastq_dict[key][v] for v in fastq_dict[key].keys() if os.path.isfile(fastq_dict[key][v])][0]
            solo_message = "\n#=====\nNo pair for non-interleaved FASTQ reads %s found, file may contain errors.\nPlease include read pair and/or check if the sequences are correctly interleaved or single-end.\n#=====\n" % solo
            logging.warning(solo_message)
            if exclude:
                #Remove entries from dictionary with errors
                del fastq_dict[key]
                logging.info('exclude == True, removed non-paired error FASTQ %s from dictionary' % solo)

            else:
                pass

    return fastq_dict
#=============================================================================
