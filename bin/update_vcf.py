#! /usr/bin/env python2.7
"""
Update the sequence id and coordinates of a VCF file using an alignment file generated by the fasta_diff program
"""

__version__ = '1.0'

from collections import defaultdict
import logging
import re
logging.basicConfig(level=logging.DEBUG, format='%(levelname)-8s %(message)s')


class VCFUpdater(object):
    """
    Initialize a VCFUpdater instance with an alignment_file
    vcf_updater = VCFUpdater(alignment_file)
    Update a vcf_file with the update method
    vcf_updater.update(vcf_file)
    """

    def __init__(self, alignment_list_tsv_file, updated_postfix, removed_postfix):
        self.alignment_list = VCFUpdater.read_alignment_list_tsv(alignment_list_tsv_file)
        # create a dictionary to lookup alignments using old_id as key
        self.alignment_dict = defaultdict(list)
        for a in self.alignment_list:
            self.alignment_dict[a[0]].append(a)
        self.updated_postfix = updated_postfix
        self.removed_postfix = removed_postfix

    def update(self, vcf_file, reference):
        """
        Updates the vcf_file using the alignment_list_tsv_file
        :param str vcf_file: The VCF file to be updated
        """
        logging.info('Processing VCF file: %s...', vcf_file)
        self.vcf_file = vcf_file
        self.reference = reference
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

    def fasta_file_sequence_length(self):
        sequence_length = dict()
        sequence_id = None
        with open(self.reference, 'rb') as fasta_file_f:
            for line in fasta_file_f:
                line = line.strip()
                if len(line) != 0:
                    if line[0] == '>':
                        lines = line.split(' ')
                        if sequence_id == None:
                            # the first sequence
                            sequence_id = lines[0][1:]
                        elif sequence_id != lines[0][1:]:
                            # next sequence
                            sequence_length[sequence_id]['state'] = True
                            sequence_id = lines[0][1:]
                        if sequence_id not in sequence_length:
                            sequence_length[sequence_id] = {
                                'length': 0,
                                'state': False
                            }
                        else:
                            if sequence_length[sequence_id]['state'] == True:
                                logging.warning('Duplicate ID found! %s' % (sequence_id))
                    else:
                        sequence_length[sequence_id]['length'] += len(line)
        return sequence_length
    def _update_features(self):
        """
        Goes through the VCF file, updating the ids and coordinates of each feature and
        checks for removed ids and coordinates.
        """
        from os.path import splitext, basename
        VCF_root, VCF_ext = splitext(self.vcf_file)
        updated_file = VCF_root + self.updated_postfix + VCF_ext
        removed_file = VCF_root + self.removed_postfix + VCF_ext
        updated_count = 0
        removed_count = 0
        updated_file_f = open(updated_file, 'wb')
        removed_file_f = open(removed_file, 'wb')
        contig_key = ['ID', 'length']
        contig_sort_map = defaultdict(int, zip(contig_key, range(len(contig_key), 0, -1)))
        if self.reference:
            sequence_length = self.fasta_file_sequence_length()
        with open(self.vcf_file, 'rb') as in_f:
            for line in in_f:
                if len(line.strip()) == 0:
                    # ingore blank line
                    continue
                line_strip = line.strip()
                # meta-information line
                # only update ##contig, ##reference
                if line_strip.startswith('##'):
                    meta_information = line_strip
                    if meta_information.startswith('##contig'):
                        contig = re.search('##contig=<(.+)>', meta_information).group(1)
                        contig_dict = dict(re.findall('([^=,]+)=([^=,\n]+)', contig))
                        # update ID
                        if 'ID' in contig_dict and contig_dict['ID'] in self.alignment_dict:
                            mappings = self.alignment_dict[contig_dict['ID']]
                            mappings_dict = dict()
                            for mapping in mappings:
                                if mapping[3] not in mappings_dict:
                                    mappings_dict[mapping[3]] = {'min':mapping[4],'max':mapping[5],'length':mapping[5]-mapping[4]}
                                else:
                                    mappings_dict[mapping[3]]['min'] = min(mappings_dict[mapping[3]]['min'], mapping[4])
                                    mappings_dict[mapping[3]]['max'] = max(mappings_dict[mapping[3]]['max'], mapping[5])
                                    mappings_dict[mapping[3]]['length'] = mappings_dict[mapping[3]]['max'] - mappings_dict[mapping[3]]['min']
                            flag = False
                            for newid in mappings_dict:
                                contig_dict['ID'] = newid
                                if 'length' in contig_dict:
                                    if self.reference:
                                        contig_dict['length'] = str(sequence_length[newid]['length'])
                                    else:
                                        contig_dict['length'] = mappings_dict[newid]['length']
                                contig_list = []
                                for k, v in sorted(contig_dict.items(), key=lambda x: contig_sort_map[x[0]], reverse=True):
                                    contig_list.append('%s=%s' % (str(k), str(v)))
                                updated_contig = ','.join(contig_list)
                                new_meta = re.sub(r'##contig=<(.+)>','##contig=<%s>' % updated_contig, meta_information)
                                updated_file_f.write(new_meta + '\n')
                                flag = True
                        if flag == True:
                            continue

                    elif meta_information.startswith('##reference'):
                        if self.reference:
                            reference_filename = basename(self.reference)
                            meta_information = re.sub(r'##reference=(.+)','##reference=' + reference_filename, meta_information)
                    updated_file_f.write(meta_information + '\n')
                    removed_file_f.write(meta_information + '\n')
                elif line_strip.startswith('#'):
                    updated_file_f.write(line_strip + '\n')
                    removed_file_f.write(line_strip + '\n')
                else:
                    tokens = line_strip.split('\t')
                    if tokens[0] in self.alignment_dict:
                        start, end = int(tokens[1]), int(tokens[1]) -1 + len(tokens[3])# positive 1-based integer coordinates
                        mappings = self.alignment_dict[tokens[0]]
                        start_mapping = filter(lambda m: m[1] < start and start <= m[2], mappings)
                        end_mapping = filter(lambda m: m[1] < end and end <= m[2], mappings)
                        # we got a bad annotation if start or end pos is N
                        if len(start_mapping) != 1 or len(end_mapping) != 1:
                            removed_file_f.write(line_strip + '\n')
                            removed_count += 1
                        else:
                            if start_mapping[0][3] != end_mapping[0][3]:
                                removed_file_f.write(line_strip + '\n')
                                removed_count += 1
                            else:
                                tokens[0] = start_mapping[0][3]
                                tokens[1] = str(start - start_mapping[0][1] + start_mapping[0][4])
                                keep = '\t'.join(tokens)
                                updated_file_f.write(keep + '\n')
                                updated_count += 1

                    else:
                        removed_file_f.write(line_strip + '\n')
                        removed_count += 1
        updated_file_f.close()
        removed_file_f.close()
        logging.info('  Updated lines: %d', updated_count)
        logging.info('  Removed lines: %d', removed_count)


if __name__ == '__main__':
    import sys
    import argparse
    from textwrap import dedent
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter, description=dedent("""\
    Update the sequence id and coordinates of a VCF file using an alignment file generated by the fasta_diff.py program.
    Updated features are written to a new file with '_updated'(default) appended to the original VCF file name.
    Feature that can not be updated, due to the id being removed completely or the feature contains regions that
    are removed or replaced with Ns, are written to a new file with '_removed'(default) appended to the original VCF file name.
    Example:
        fasta_diff.py old.fa new.fa | %(prog)s a.vcf b.vcf c.vcf
    """))
    parser.add_argument('vcf_files', metavar='VCF_FILE', nargs='+', type=str, help='List one or more VCF files to be updated')
    parser.add_argument('-a', '--alignment_file', type=argparse.FileType('rb'), default=sys.stdin,
                        help='The alignment file generated by fasta_diff.py, a TSV file with 6 columns: old_id, old_start, old_end, new_id, new_start, new_end (default: STDIN)')
    parser.add_argument('-u', '--updated_postfix', default='_updated',
                        help='The filename postfix for updated features (default: "_updated")')
    parser.add_argument('-r', '--removed_postfix', default='_removed',
                        help='The filename postfix for removed features (default: "_removed")')
    parser.add_argument('-ref', '--reference',
                        help='The new reference genome fasta file')
    parser.add_argument('-v', '--version', action='version', version='%(prog)s ' + __version__)

    args = parser.parse_args()
    vcf_updater = VCFUpdater(args.alignment_file, args.updated_postfix, args.removed_postfix)
    for vcf_file in args.vcf_files:
        vcf_updater.update(vcf_file, args.reference)
