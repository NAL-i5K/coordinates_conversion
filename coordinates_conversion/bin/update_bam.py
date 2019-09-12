#! /usr/bin/env python3

"""
Update the sequence id and coordinates of a Bam file using an alignment file generated by the fasta_diff program

"""

__version__ = '2.0'

from collections import defaultdict
import logging
import pysam
import sys
import argparse
from textwrap import dedent

logging.basicConfig(level=logging.DEBUG, format='%(levelname)-8s %(message)s')

class BamUpdater(object):

    def __init__(self, alignment_list_tsv_file, updated_postfix, removed_postfix):
        self.alignment_list = BamUpdater.read_alignment_list_tsv(alignment_list_tsv_file)
        # create a dictionary to lookup alignments using old_id as key
        self.alignment_dict = defaultdict(list)
        for a in self.alignment_list:
            self.alignment_dict[a[0]].append(a)
        self.updated_postfix = updated_postfix
        self.removed_postfix = removed_postfix


    def update(self, bam_file):
        """
        Updates the bam_file using the alignment_list_tsv_file
        :param str bam_file: The bam file to be updated
        """
        logging.info('  Processing bam file: %s...', bam_file)
        self.bam_file = bam_file
        self._update_features()


    @staticmethod
    def read_alignment_list_tsv(alignment_list_tsv_file):
        """
        Parse an alignment_list_tsv and returns a list
        :param alignment_list_tsv_file: The output alignment file of fasta_diff
        :type alignment_list_tsv_file: string or file
        :return: a list of [old_id, old_start, old_end, new_id, new_start, new_end]
        """
        tsv_format = [str, int, int, str, int, int]
        alignment_list_tsv_file_f = alignment_list_tsv_file
        if isinstance(alignment_list_tsv_file, str):
            logging.info('  Reading alignment data from: %s...', alignment_list_tsv_file_f.name)
            alignment_list_tsv_file_f = open(alignment_list_tsv_file, 'rb')

        alignment_list = []
        for line in alignment_list_tsv_file_f:
            line = str(line,'utf-8')
            alignment_list.append([f(t) for f, t in zip(tsv_format, line.split('\t'))])

        if isinstance(alignment_list_tsv_file, str):
            alignment_list_tsv_file_f.close()
        else:
            logging.info('  Reading alignment data from: %s...', alignment_list_tsv_file_f.name)
        logging.info('  Alignments: %d', len(alignment_list))

        return alignment_list

    def _update_features(self):
        """
        Goes through the bam file, updating the reference sequence names and coordinates of each feature and
        checks for removed reference sequence names and coordinates.
        """
        from os.path import splitext
        bam_root, bam_ext = splitext(self.bam_file)
        updated_file = bam_root + self.updated_postfix + bam_ext
        removed_file = bam_root + self.removed_postfix + bam_ext
        updated_count = 0
        removed_count = 0
        Program_ID = ""
        Program_list = ["TopHat","STAR","Bowtie","BWA"]
        in_f = pysam.AlignmentFile(self.bam_file, 'rb')
        bam_header_old = in_f.header
        bam_header_new = {}
        header_new_ID_dict = {}
        header_num_count = 0
        written_ID = set()
        logging.info('  Update tag SN and LN....')
        for k,v in bam_header_old.iteritems():
            if k == 'SQ':
                bam_header_new[k] = []
                for reference_sequence_dict in v:
                    if 'SN' in reference_sequence_dict:
                        if reference_sequence_dict['SN'] in self.alignment_dict:
                            mappings = sorted(self.alignment_dict[reference_sequence_dict['SN']], key=lambda x: (x[3], x[4]))
                            for mapping in mappings:
                                if mapping[3] not in written_ID:
                                    header_new_ID_dict[(mapping[0], mapping[3])] = header_num_count
                                    updated_reference_sequence_dict = dict(reference_sequence_dict)
                                    updated_reference_sequence_dict['SN'] = mapping[3]
                                    if 'LN' in updated_reference_sequence_dict:
                                        filter_mappings = [m for m in mappings if m[3] == mapping[3]]
                                        updated_reference_sequence_dict['LN'] = filter_mappings[len(filter_mappings)-1][5] - filter_mappings[0][4]
                                    bam_header_new[k].append(updated_reference_sequence_dict)
                                    header_num_count += 1
                                    written_ID.add(mapping[3])
            elif k =='PG':
                bam_header_new[k] = v
                if v[0]['ID'] in Program_list:
                    Program_ID = v[0]['ID']
                    logging.info('  Program record identifier %s',Program_ID)
            else:
                bam_header_new[k] = v
        logging.info('  Only some alignment sections will be updated! Tags won\'t be updated!')
        Custom_tag = []
        if Program_ID == 'TopHat':
            Custom_tag = ['XG','XM','XN','XO','XS']
        elif Program_ID =='Bowtie':
            Custom_tag = ['XA']
        elif Program_ID =='STAR':
            Custom_tag = ['ji']
        elif Program_ID == 'BWA':
            Custom_tag = ['XA','X0','X1','XN','XM','XO','XG','XT','XS','XF','XE']

        updated_file_f = pysam.AlignmentFile(updated_file,'wb',header = bam_header_new)
        removed_file_f = pysam.AlignmentFile(removed_file,'wb',header = in_f.header)
        bam_update_list = []

        for read in in_f.fetch(until_eof=True):
            #If a alignment doesn't have reference name, it will directly add to updated bam file.
            if read.reference_id == -1:
                updated_count+=1
                updated_file_f.write(read)
            else:
                if read.reference_name in self.alignment_dict:
                    start, end = int(read.reference_start), int(read.reference_end) # positive 0-based integer coordinates
                    mappings = self.alignment_dict[read.reference_name]
                    start_mapping = [m for m in mappings if m[1] <= start and start <= m[2]]
                    end_mapping = [m for m in mappings if m[1] <= end and end <= m[2]]
                        # we got a bad annotation if start or end pos is N
                    if len(start_mapping) != 1 or len(end_mapping) != 1:
                        removed_count+=1
                        removed_file_f.write(read)
                    else:
                        if start_mapping[0][3] != end_mapping[0][3]:
                            removed_count+=1
                            removed_file_f.write(read)
                            continue
                        read_out = pysam.AlignedSegment()
                        read_out.query_name = read.query_name
                        read_out.flag = read.flag
                        read_out.reference_id = header_new_ID_dict[(read.reference_name, start_mapping[0][3])]
                        read_out.reference_start = int(start - start_mapping[0][1] + start_mapping[0][4])
                        read_out.mapping_quality = read.mapping_quality
                        read_out.cigar = read.cigar
                        read_out.template_length = read.template_length
                        read_out.query_sequence = read.query_sequence
                        read_out.query_qualities = read.query_qualities
                        read_out.tags= read.tags
                        if read.next_reference_id != -1:
                            if read.next_reference_name in self.alignment_dict:
                                next_mappings = self.alignment_dict[read.next_reference_name]
                                next_start = int(read.next_reference_start)
                                start_next_mapping = filter(lambda m: m[1] <= next_start and next_start <= m[2], next_mappings)
                                if len(start_next_mapping)!=1:
                                    removed_count+=1
                                    removed_file_f.write(read)
                                    continue
                                else:
                                    read_out.next_reference_id = header_new_ID_dict[(read.next_reference_name,start_next_mapping[0][3])]
                                    read_out.next_reference_start = int(next_start - start_next_mapping[0][1] + start_next_mapping[0][4])
                            else:
                                removed_count+=1
                                removed_file_f.write(read)
                                continue
                        elif read.next_reference_id == -1:
                            read_out.next_reference_id = read.next_reference_id
                            read_out.next_reference_start = read.next_reference_start

                    updated_count +=1
                    updated_file_f.write(read_out)
                else:
                    removed_count+=1
                    removed_file_f.write(read)

        in_f.close()
        updated_file_f.close()
        removed_file_f.close()
        logging.info('  Updated Alignments: %d', updated_count)
        logging.info('  Removed Alignments: %d', removed_count)


def main():
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter, description=dedent("""\
    Update the reference sequence name of the alignment and coordinates of a bam file using an alignment file generated by the fasta_diff program.
    Updated features are written to a new file with '_updated'(default) appended to the original bam file name.
    Feature that can not be updated, due to the id being removed completely or the feature contains regions that
    are removed or replaced with Ns, are written to a new file with '_removed'(default) appended to the original bam file name.

    Example:
        fasta_diff example_file/old.fa example_file/new.fa | %(prog)s example_file/example.bam
    """))
    parser.add_argument('bam_files', metavar='bam_FILE', nargs='+', type=str, help='List one or more bam files to be updated')
    parser.add_argument('-a', '--alignment_file', type=argparse.FileType('rb'), default=sys.stdin,
                        help='The alignment file generated by fasta_diff, a TSV file with 6 columns: old_id, old_start, old_end, new_id, new_start, new_end (default: STDIN)')
    parser.add_argument('-u', '--updated_postfix', default='_updated',
                        help='The filename postfix for updated features (default: "_updated")')
    parser.add_argument('-r', '--removed_postfix', default='_removed',
                        help='The filename postfix for removed features (default: "_removed")')
    parser.add_argument('-v', '--version', action='version', version='%(prog)s ' + __version__)


    args = parser.parse_args()
    Bam_updater = BamUpdater(args.alignment_file, args.updated_postfix, args.removed_postfix)
    for Bam_file in args.bam_files:
        Bam_updater.update(Bam_file)

if __name__ == '__main__':
    main()
