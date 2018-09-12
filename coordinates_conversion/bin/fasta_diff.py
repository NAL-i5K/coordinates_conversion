#! /usr/bin/env python2.7

"""
Compare two very similar FASTA files and output coordinate mappings
"""

__version__ = '1.1'

try:
   import cPickle as pickle
except:
   import pickle
from os.path import isfile
from os import remove
from collections import OrderedDict
import logging
import sys
import argparse
from textwrap import dedent

logging.basicConfig(level=logging.DEBUG, format='%(levelname)-8s %(message)s')


def fasta_file_to_dict(fasta_file, id=True, header=False, seq=False):
    """Returns a dict from a fasta file and the number of sequences as the second return value.
    fasta_file can be a string path or a file object.
    The key of fasta_dict can be set using the keyword arguments and
    results in a combination of id, header, sequence, in that order. joined with '||'. (default: id)
    Duplicate keys are checked and a warning is logged if found.
    The value of fasta_dict is a python dict with 3 keys: header, id and seq
    """

    fasta_file_f = fasta_file
    if isinstance(fasta_file, str):
        fasta_file_f = open(fasta_file, 'rb')

    fasta_dict = dict()
    flags = OrderedDict([('id', id), ('header', header), ('seq', seq)])
    entry = dict([('id', ''), ('header', ''), ('seq', '')])
    count = 0
    line_num = 0

    for line in fasta_file_f:
        line = line.strip()
        if not line:
            continue
        if line[0] == '>':
            count += 1
            key = '||'.join([entry[i] for i in flags if flags[i]])
            if key: # key != ''
                if key in fasta_dict: # check for duplicate key
                    logging.warning('%s : Line %d : Duplicate %s [%s] : ID = [%s].', fasta_file_f.name, line_num, '||'.join([i for i in flags if flags[i]]), key[:25] + (key[25:] and '..'), entry['id'])
                fasta_dict[key] = entry
                entry = dict()
            entry['header'] = line
            entry['id'] = line.split()[0][1:]
            entry['seq'] = ''
        else:
            entry['seq'] += line.upper()
        line_num += 1

    if isinstance(fasta_file, str):
        fasta_file_f.close()

    key = '||'.join([entry[i] for i in flags if flags[i]])
    if key: # key != ''
        if key in fasta_dict:
            logging.warning('%s : Line %d : Duplicate %s [%s] : ID = [%s].', fasta_file_f.name, line_num, '||'.join([i for i in flags if flags[i]]), key[:25] + (key[25:] and '..'), entry['id'])
        fasta_dict[key] = entry

    return fasta_dict, count

def fasta_dict_to_file(fasta_dict, fasta_file):
    """Write fasta_dict to a fasta_file
    fasta_file can be a string path or a file object
    The key of fasta_dict doesn't matter
    The value of fasta_dict is a python dict with 3 keys: header, id and seq
    """
    if isinstance(fasta_file, str):
        fasta_file = open(fasta_file, 'wb')

    for key in fasta_dict:
        fasta_file.write(fasta_dict[key]['header'] + '\n' + fasta_dict[key]['seq'] + '\n')

    fasta_file.close()

def query_yes_no(question, default='yes'):
    """Ask a yes/no question via raw_input() and return their answer.

    'question' is a string that is presented to the user.
    'default' is the presumed answer if the user just hits <Enter>.
        It must be 'yes' (the default), 'no' or None (meaning
        an answer is required of the user).

    The 'answer' return value is one of 'yes' or 'no'.
    """
    valid = {'yes': True, 'y': True, 'ye': True,
             'no': False, 'n': False}
    if default is None:
        prompt = ' [y/n] '
    elif default == 'yes':
        prompt = ' [Y/n] '
    elif default == 'no':
        prompt = ' [y/N] '
    else:
        raise ValueError('invalid default answer: "%s"' % default)

    while True:
        sys.stderr.write(question + prompt)
        choice = raw_input().strip().lower()
        if default is not None and choice == '':
            return valid[default]
        elif choice in valid:
            return valid[choice]
        else:
            sys.stderr.write('Please respond with "y" or "n".\n')

def fasta_diff(old_fasta_file, new_fasta_file, debug=True, header_check=False, report=None):
    """
    Compares two very similar FASTA files and outputs coordinate mappings using a multi stage algorithm:
    Stage 1: Find 100% matches
    Stage 2: Find 100% substrings, where the full length of a new sequence can be found as a substring of a old sequence
    Stage 3: Find cases where part of the sequence was converted into Ns
    Originally created to compare the NCBI version to the original reference of Diaphorina citri
    :param old_fasta_file: a string path
    :param new_fasta_file: a string path
    :param debug: If True, partial results are saved in a *_stage_i_pickle file after each stage, unmatched sequences after each stage are saved in a *_stage_i_unmatched file for both FASTA files
    :return: a list of [old_id, old_start, old_end, new_id, new_start, new_end], 0-based coordinate system
    """
    alignment_list = []
    onetomultiple=dict()
    #for stage4 dictionary
    old_fasta_dict, new_fasta_dict = None, None

    def match_identical_sequence():
        # find 100% matches
        #[old_id, old_start, old_end, new_id, new_start, new_end]
        old_seqs = old_fasta_dict.keys()
        for seq in old_seqs:
            if seq in new_fasta_dict:
                alignment = [old_fasta_dict[seq]['id'], 0, len(seq), new_fasta_dict[seq]['id'], 0, len(seq)]
                if header_check and not old_fasta_dict[seq]['id'].split('.')[0] in new_fasta_dict[seq]['header']:
                    logging.warning('Failed header check (match_identical_sequence): %s -> %s', old_fasta_dict[seq]['id'], new_fasta_dict[seq]['id'])
                alignment_list.append(alignment)
                del old_fasta_dict[seq]
                del new_fasta_dict[seq]

    def match_truncated_sequence():
        # Find 100% substrings, where the full length of a new sequence can be found as a substring of a old sequence
        old_seqs = old_fasta_dict.keys()
        new_seqs = new_fasta_dict.keys()
        match_truncated = dict() # {matches[0]: {'matches': {new_seq}, 'alignment': [oldid, oldstart, oldend, newid, newstart, newend]}
        match_truncated_order = list()
        for new_seq in new_seqs:
            matches = filter(lambda old_seq: new_seq in old_seq, old_seqs)
            if len(matches) == 1:
                start = matches[0].find(new_seq)
                alignment = [old_fasta_dict[matches[0]]['id'], start, start + len(new_seq), new_fasta_dict[new_seq]['id'], 0, len(new_seq)]
                if header_check and not old_fasta_dict[matches[0]]['id'].split('.')[0] in new_fasta_dict[new_seq]['header']:
                    logging.warning('Failed header check (match_truncated_sequence): %s -> %s', old_fasta_dict[matches[0]]['id'], new_fasta_dict[new_seq]['id'])
                # need to check if this match is not a one to multiple mapping
                if matches[0] not in match_truncated:
                    match_truncated_order.append(matches[0])
                    match_truncated[matches[0]] = {
                        'matches': list(),
                        'alignment': None
                    }
                match_truncated[matches[0]]['matches'].append(new_seq)
                match_truncated[matches[0]]['alignment'] = alignment
            elif len(matches) > 1:
                logging.warning('Failed one to one mapping: %s has %d matches: %s\n' % (new_fasta_dict[new_seq]['id'], len(matches), ','.join([old_fasta_dict[x]['id'].split('.')[0] for x in matches])))
        for match in match_truncated_order:
            if len(match_truncated[match]['matches']) == 1:
                # one to one
                alignment_list.append(match_truncated[match]['alignment'])
                del old_fasta_dict[match] # matches[0]
                del new_fasta_dict[match_truncated[match]['matches'][0]] # new_seq

    def match_split_subsequence():
        # Find cases where part of the sequence was converted into Ns
        # algorithm:
        # 1. split new seq by Ns into multiple substrings
        # 2. find an old seq with all substrings and also in order
        # 3. if found, write an alignment line for every substring
        old_seqs = old_fasta_dict.keys()
        new_seqs = new_fasta_dict.keys()
        match_split = dict() # {matches[0]: {'matches': {new_seq}, 'alignment': [oldid, oldstart, oldend, newid, newstart, newend]}

        match_split_order = list()
        for new_seq in new_seqs:
            segments = new_seq.replace('N', ' ').split()
            matches = []

            for old_seq in old_seqs:
                start_original = 0
                start_new = 0
                seg_matches = []
                for segment in segments:
                    pos = old_seq.find(segment, start_original)
                    if pos == -1:
                        break
                    else:
                        start_original = pos + len(segment)
                        pos_new = new_seq.find(segment, start_new)
                        start_new = pos_new + len(segment)
                        seg_matches.append([old_fasta_dict[old_seq]['id'], pos, pos + len(segment), new_fasta_dict[new_seq]['id'], pos_new, pos_new + len(segment)])
                if len(segments) == len(seg_matches):
                    matches.append((seg_matches, old_seq))
            if len(matches) == 1:
                if header_check and not matches[0][0][0][0].split('.')[0] in new_fasta_dict[new_seq]['header']:
                    logging.warning('Failed header check (match_split_subsequence): %s -> %s', matches[0][0][0][0], matches[0][0][0][3])
                new_matches = []
                tmp_oldstart = 0
                tmp_oldend = 0
                tmp_newstart = 0
                tmp_newend = 0

                for match in matches[0][0]:
                    if match == matches[0][0][0]:
                        tmp_oldstart = match[1]
                        tmp_oldend = match[2]
                        tmp_newstart = match[4]
                        tmp_newend = match[5]
                    elif matches[0][1][tmp_oldstart:match[2]] == new_seq[tmp_newstart:match[5]]:
                        tmp_oldend = match[2]
                        tmp_newend = match[5]
                    else:
                        for nucl in matches[0][1][tmp_oldend:]:
                            try:
                                if new_seq[tmp_newend] == nucl:
                                    tmp_oldend += 1
                                    tmp_newend += 1
                                else:
                                    break
                            except:
                                break
                        new_matches.append([match[0], tmp_oldstart, tmp_oldend, match[3], tmp_newstart, tmp_newend])
                        tmp_oldstart = match[1]
                        tmp_oldend = match[2]
                        tmp_newstart = match[4]
                        tmp_newend = match[5]
                #check if the new sequence end with N
                for nucl in matches[0][1][tmp_oldend:]:
                    try:
                        if new_seq[tmp_newend] == nucl:
                            tmp_oldend += 1
                            tmp_newend += 1
                        else:
                            break
                    except:
                        break

                new_matches.append([match[0], tmp_oldstart, tmp_oldend, match[3], tmp_newstart, tmp_newend])
                if matches[0][1] not in match_split:
                    match_split_order.append(matches[0][1])
                    match_split[matches[0][1]] = {
                        'matches': list(),
                        'alignment': list()
                    }
                match_split[matches[0][1]]['matches'].append(new_seq)
                match_split[matches[0][1]]['alignment'].extend(new_matches)

            elif len(matches) > 1:
                logging.warning('Failed one to one mapping: %s has %d matches: %s\n' % (new_fasta_dict[new_seq]['id'], len(matches), ','.join([x[0][0][0].split('.')[0] for x in matches])))

        onetomultiple.update(match_split)
        for match in match_split_order:
            if len(match_split[match]['matches']) == 1:
                # one to one
                alignment_list.extend(match_split[match]['alignment'])
                del old_fasta_dict[match] # matches[0]
                del new_fasta_dict[match_split[match]['matches'][0]] # new_seq
                del match_split[match]

    def one_to_multiple_match():
        stagelist=list()
        for match in onetomultiple:
            # one to mutiple
            if len(onetomultiple[match]['matches']) > 1:
                stagelist.extend(onetomultiple[match]['alignment'])

                overlap=dict()
                #Fetch old_id and new_id
                for tmp in stagelist:
                    pair=((tmp[0]),(tmp[3]))
                    if pair not in overlap:
                       overlap[pair]=set()
                    overlap[pair].update((tmp[1],tmp[2]))
                pairs=overlap.keys()
                pairs_sort=sorted(pairs, key=lambda aaa:aaa[0])
                run_sort=set()
                #Record the pair that has been run
                delete_pairs=set()

                for pair1 in pairs_sort:
                    for pair2 in pairs_sort:
                        if pair1==pair2:
                            continue
                        else:
                            if (pair1,pair2) in run_sort:
                                continue
                            run_sort.update([(pair1,pair2),(pair2,pair1)])
                            if pair1[0]==pair2[0]:
                                # find overlap pair
                                if (min(overlap[pair1])<=min(overlap[pair2]) and min(overlap[pair2])<=max(overlap[pair1]))\
                                or (min(overlap[pair1])<=max(overlap[pair2]) and max(overlap[pair2])<=max(overlap[pair1])):

                                    delete_pairs.update((pair1,pair2))
                                    # add overlap pair to delete_pairs

                stage_four_result=list()
                for delete in stagelist:
                   if (delete[0],delete[3]) not in delete_pairs:
                       stage_four_result.append(delete)

                       if match in old_fasta_dict:
                            del old_fasta_dict[match]

                       for new in onetomultiple[match]['matches']:
                           if new in new_fasta_dict:
                                if new_fasta_dict[new]['id']==delete[3]:
                                    del new_fasta_dict[new]

        if onetomultiple:
            alignment_list.extend(stage_four_result)
        # add empty to final result

    stages = [match_identical_sequence, match_truncated_sequence, match_split_subsequence,one_to_multiple_match]
    matched_sequence_count = 0
    for stage in range(len(stages)):
        temp_file_name = ''
        if debug:
            temp_file_name = old_fasta_file + '_' + new_fasta_file + '_stage_' + str(stage + 1) + '_pickle'
        if isfile(temp_file_name):
            (alignment_list, old_fasta_dict, new_fasta_dict) = pickle.load(open(temp_file_name, 'rb'))
        else:
            if old_fasta_dict == None and new_fasta_dict == None:
                # use sequence as dict key
                #163023 unique sequences in original fasta file
                #161988 unique sequences in new fasta file
                logging.info('Reading old FASTA file (%s)...', old_fasta_file)
                old_fasta_dict, old_fasta_count = fasta_file_to_dict(old_fasta_file, id=False, header=False, seq=True)
                if not len(old_fasta_dict) == old_fasta_count:
                    logging.warning('  Duplicate sequences detected in old FASTA file, %d unique sequences out of a total of %d sequences', len(old_fasta_dict), old_fasta_count)
                    if not query_yes_no('Ignore and continue?'):
                        sys.exit(1)
                else:
                    logging.info('  Unique sequences: %d', len(old_fasta_dict))

                logging.info('Reading new FASTA file (%s)...', new_fasta_file)
                new_fasta_dict, new_fasta_count = fasta_file_to_dict(new_fasta_file, id=False, header=False, seq=True)
                if not len(new_fasta_dict) == new_fasta_count:
                    logging.warning('Duplicate sequences detected in new FASTA file, %d unique sequences out of a total of %d sequences', len(new_fasta_dict), new_fasta_count)
                    if not query_yes_no('Ignore and continue?'):
                        sys.exit(1)
                else:
                    logging.info('  Unique sequences: %d', len(new_fasta_dict))

            logging.info('Stage %d - %s:', stage + 1, stages[stage].__name__)
            stages[stage]()

            if debug:
                pickle.dump((alignment_list, old_fasta_dict, new_fasta_dict), open(temp_file_name, 'wb'))

        new_matched_sequence_count = len(set([a[0] for a in alignment_list]))
        logging.info('  Matched sequences: %d (New : %d)', new_matched_sequence_count, new_matched_sequence_count - matched_sequence_count)
        logging.info('  Unmatched sequences in old FASTA: %d', len(old_fasta_dict))
        logging.info('  Unmatched sequences in new FASTA: %d', len(new_fasta_dict))


        matched_sequence_count = new_matched_sequence_count
        if debug:
            fasta_dict_to_file(old_fasta_dict, old_fasta_file + '_stage_' + str(stage + 1) + '_unmatched')
            fasta_dict_to_file(new_fasta_dict, new_fasta_file + '_stage_' + str(stage + 1) + '_unmatched')
        if report is not None:
            report_header = True
            if isfile(report):
                report_header = False
            with open(report, 'a') as out_report:
                if report_header:
                    out_report.write('#Sequence_ID\tSequence_length\tUnmatched_stage\n')
                for unmatched in new_fasta_dict:
                    out_report.write('\t'.join([new_fasta_dict[unmatched]['id'], str(len(new_fasta_dict[unmatched]['seq'])), 'Stage %d' % (stage + 1)]) + '\n')

    return alignment_list, old_fasta_dict, new_fasta_dict

def main():
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter, description=dedent("""\
    Compares two very similar FASTA files and outputs coordinate mappings using a multi stage algorithm:
    Stage 1: Find 100% matches
    Stage 2: Find 100% substrings, where the full length of a new sequence can be found as a substring of a old sequence
    Stage 3: Find cases where part of the sequence was converted into Ns
    Outputs the 6 columns as tab-separated values: old_id, old_start, old_end, new_id, new_start, new_end
    Originally created to compare the NCBI version to the original reference of Diaphorina citri

    Example:
        fasta_diff example_file/old.fa example_file/new.fa -o match.tsv -r report.txt
    """))
    parser.add_argument('old_fasta', type=str, help='The original FASTA file')
    parser.add_argument('new_fasta', type=str, help='The new FASTA file')
    parser.add_argument('-o', '--out', nargs='?', type=argparse.FileType('wb'), default=sys.stdout, help='The output alignment file (default: STDOUT)')
    parser.add_argument('-r', '--report', type=str, help='Generate a report for the unmatched sequences in new FASTA file.')
    parser.add_argument('-d', '--debug', action='store_true',
                        help='If set, partial results are saved in a *_stage_i_pickle file after each stage, unmatched sequences after each stage are saved in a *_stage_i_unmatched file for both FASTA files.')
    parser.add_argument('-hc', '--header_check', action='store_true',
                        help='If set, confirm the detected mapping by checking the header of the new sequence for the id of the mapped old sequence.')
    parser.add_argument('-v', '--version', action='version', version='%(prog)s ' + __version__)

    test_lv = 0 # debug
    if test_lv == 1:
        args = parser.parse_args(['diaci1.1.fa', '121845_ref_Diaci_psyllid_genome_assembly_version_1.1_chrUn.fa'])
        alignment_list, old_fasta_dict, new_fasta_dict = fasta_diff(args.old_fasta, args.new_fasta, debug=True)

        # Looks like we're done! Output alignment_list
        pickle.dump(alignment_list, open('alignment_list.pickle', 'wb'))
        with open('alignment_list.tsv', 'wb') as f:
            for alignment in alignment_list:
                f.write('\t'.join([str(a) for a in alignment]) + '\n')
    else:
        args = parser.parse_args()
        # remove existing report
        if args.report:
            if isfile(args.report):
                remove(args.report)
        alignment_list, old_fasta_dict, new_fasta_dict = fasta_diff(args.old_fasta, args.new_fasta, debug=args.debug, header_check=args.header_check, report=args.report)
        if args.debug:
            alignment_list_pickle_file = args.out.name + '_pickle'
            if args.out.name == '<stdout>':
                alignment_list_pickle_file = 'alignment_list_pickle'
            pickle.dump(alignment_list, open(alignment_list_pickle_file, 'wb'))
        for alignment in alignment_list:
            args.out.write('\t'.join([str(a) for a in alignment]) + '\n')
        args.out.close()


if __name__ == '__main__':
    main()
