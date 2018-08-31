#! /usr/bin/env python2.7

"""
Update the sequence id and coordinates of a BedGraph file using an alignment file generated by the fasta_diff program

"""

__version__ = '1.1'

from collections import defaultdict
import logging

logging.basicConfig(level=logging.DEBUG, format='%(levelname)-8s %(message)s')

class BedGraphUpdater(object):

    def __init__(self, alignment_list_tsv_file, updated_postfix, removed_postfix):
        self.alignment_list = BedGraphUpdater.read_alignment_list_tsv(alignment_list_tsv_file)
        self.alignment_dict = defaultdict(list)
        for a in self.alignment_list:
            self.alignment_dict[a[0]].append(a)
        self.updated_postfix = updated_postfix
        self.removed_postfix = removed_postfix

    def update(self, BedGraph_file):
        """
        Updates the BedGraph_file using the alignment_list_tsv_file
        :param str BedGraph_file: The BedGraph file to be updated
        """
        logging.info('Processing BedGraph file: %s...', BedGraph_file)
        self.BedGraph_file = BedGraph_file
        self._update_features()

    @staticmethod
    def read_alignment_list_tsv(alignment_list_tsv_file):
        """
        Parse an alignment_list_tsv and returns a list
        :param alignment_list_tsv_file: The output alignment file of fasta_diff.py
        :type alignment_list_tsv_file: string or file
        :return: a list of [old_id, old_start, old_end, new_id, new_start, new_end]
        """
        tsv_format = [str, int, int, str, int, int]
        alignment_list_tsv_file_f = alignment_list_tsv_file
        if isinstance(alignment_list_tsv_file, str):
           logging.info('Reading alignment data from: %s...', alignment_list_tsv_file_f.name)
           alignment_list_tsv_file_f = open(alignment_list_tsv_file, 'rb')
        alignment_list = []
        for line in alignment_list_tsv_file_f:
            alignment_list.append([f(t) for f, t in zip(tsv_format, line.split('\t'))])
        if isinstance(alignment_list_tsv_file, str):
            alignment_list_tsv_file_f.close()
        else:
            logging.info('Reading alignment data from: %s...', alignment_list_tsv_file_f.name)
        logging.info('  Alignments: %d', len(alignment_list))
        return alignment_list


    def _update_features(self):
        """
        Goes through the BedGraph file, updating the ids and coordinates of each line and
        checks for removed ids and coordinates.
        """
        from os.path import splitext
        BedGraph_root, BedGraph_ext = splitext(self.BedGraph_file)
        updated_file = BedGraph_root + self.updated_postfix + BedGraph_ext
        removed_file = BedGraph_root + self.removed_postfix + BedGraph_ext
        updated_count = 0
        removed_count = 0
        updated_file_f = open(updated_file, 'wb')
        removed_file_f = open(removed_file, 'wb')
        with open(self.BedGraph_file, 'rb') as in_f:
            for line in in_f:
                line_strip = line.strip()
                tokens = line_strip.split('\t')
                if tokens[0] in self.alignment_dict:
                    start, end = int(tokens[1]), int(tokens[2])
                    mappings = self.alignment_dict[tokens[0]]
                    start_mapping = filter(lambda m: m[1] <= start and start <= m[2], mappings)
                    end_mapping = filter(lambda m: m[1] <= end and end <= m[2], mappings)
                    if len(start_mapping) != 1 or len(end_mapping) != 1:
                        removed_count+=1
                        removed_file_f.write(line)
                    else:
                        if start_mapping[0][3] != end_mapping[0][3]:
                            removed_count+=1
                            removed_file_f.write(line)
                        else:
                            tokens[0] = start_mapping[0][3]# same with end_mapping[0][3]
                            tokens[1] = str(start - start_mapping[0][1] + start_mapping[0][4])
                            tokens[2] = str(end - end_mapping[0][1] + end_mapping[0][4])
                            keep = '\t'.join(tokens)
                            updated_count +=1
                            updated_file_f.write(keep+"\n")

                else:
                    removed_count+=1
                    removed_file_f.write(line)

        updated_file_f.close()
        removed_file_f.close()
        logging.info('  Updated lines: %d', updated_count)
        logging.info('  Removed lines: %d', removed_count)


if __name__ == '__main__':
    import sys
    import argparse
    from textwrap import dedent
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter, description=dedent("""\
    Update the sequence id and coordinates of a BedGraph file using an alignment file generated by the fasta_diff.py program.
    Updated Line are written to a new file with '_updated'(default) appended to the original BedGraph file name.
    Line that can not be updated, due to the id being removed completely or the line contains regions that
    are removed or replaced with Ns, are written to a new file with '_removed'(default) appended to the original BedGraph file name.

    Example:
        fasta_diff.py old.fa new.fa | %(prog)s a.BedGraph b.BedGraph c.BedGraph
    """))
    parser.add_argument('BedGraph_files', metavar='BedGraph_FILE', nargs='+', type=str, 
    help='List one or more BedGraph files to be updated')
    parser.add_argument('-a', '--alignment_file', type=argparse.FileType('rb'), default=sys.stdin,help='The alignment file generated by fasta_diff.py, a TSV file with 6 columns: old_id, old_start, old_end, new_id, new_start, new_end (default: STDIN)')

    parser.add_argument('-u', '--updated_postfix', default='_updated',
                        help='The filename postfix for updated features (default: "_updated")')
    parser.add_argument('-r', '--removed_postfix', default='_removed',
                        help='The filename postfix for removed features (default: "_removed")')
    parser.add_argument('-v', '--version', action='version', version='%(prog)s ' + __version__)
    args = parser.parse_args()
    bedgraph_updater = BedGraphUpdater(args.alignment_file, args.updated_postfix, args.removed_postfix)
    for bedgraph_file in args.BedGraph_files:
        bedgraph_updater.update(bedgraph_file)
