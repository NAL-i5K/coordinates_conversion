#! /usr/local/bin/python2.7
# Copyright (C) 2016 Limei Chiang <dytk2134 [at] gmail [dot] com>
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 3 of the License, or
# (at your option) any later version. 

"""
Update the sequence id and coordinates of a Bam file using an alignment file generated by the fasta_diff program

"""

__version__ = '1.1'

from collections import OrderedDict
from collections import defaultdict
import logging
from collections import deque
import pysam

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
        :param alignment_list_tsv_file: The output alignment file of fasta_diff.py
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
        logging.info('  Update tag SN and LN....')
        for k,v in bam_header_old.iteritems():
            if k == 'SQ':
                bam_header_new[k] = []
                for reference_sequence_dict in v:
                    if 'SN' in reference_sequence_dict:
                        if reference_sequence_dict['SN'] in self.alignment_dict:
                            header_new_ID_dict[reference_sequence_dict['SN']] = header_num_count
                            mappings = self.alignment_dict[reference_sequence_dict['SN']]
                            reference_sequence_dict['SN'] = mappings[0][3]
                            if 'LN' in reference_sequence_dict:                               
                                reference_sequence_dict['LN'] = mappings[len(mappings)-1][5] - mappings[0][4]
                            bam_header_new[k].append(reference_sequence_dict)
                            header_num_count += 1
            elif k =='PG':
                bam_header_new[k] = v
                if v[0]['ID'] in Program_list:
                    Program_ID = v[0]['ID']
                    logging.info('  Program record identifier %s',Program_ID)
            else:
                bam_header_new[k] = v
        logging.info('  Standard meaning tag: CC:Z, CP:i will be updated...')
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
                    start_mapping = filter(lambda m: m[1] <= start and start <= m[2], mappings)
                    end_mapping = filter(lambda m: m[1] <= end and end <= m[2], mappings)
                        # we got a bad annotation if start or end pos is N
                    if len(start_mapping) != 1 or len(end_mapping) != 1:
                        removed_count+=1
                        removed_file_f.write(read)
                    else:
                        read_out = pysam.AlignedSegment()
                        read_out.query_name = read.query_name
                        read_out.flag = read.flag
                        read_out.reference_id = header_new_ID_dict[read.reference_name]
                        read_out.reference_start = int(start - start_mapping[0][1] + start_mapping[0][4])
                        read_out.mapping_quality = read.mapping_quality
                        read_out.cigar = read.cigar
                        read_out.next_reference_id = read.next_reference_id
                        next_end = 0
                        
                        mappings_diff = 0
                        if read.next_reference_id != -1:
                            mappings_check = self.alignment_dict[read.next_reference_name]
                            for mappings_line in mappings_check:
                                if mappings_line[1] != mappings_line[4] or mappings_line[2] != mappings_line[5]:
                                    mappings_diff = 0
                                    break
                                elif len(mappings_check)!=1:
                                    mappings_diff = 0
                                    break
                                elif mappings_line[1] == mappings_line[4] and mappings_line[2] == mappings_line[5]:
                                    mappings_diff = 1
                            if mappings_diff == 0:
                                in_for_end = pysam.AlignmentFile(self.bam_file, 'rb')
                                for next_read in in_for_end.fetch(read.next_reference_name,read.next_reference_start,read.next_reference_start+1):
                                    if read.query_name == next_read.query_name and read.next_reference_start == next_read.reference_start:
                                        next_end = next_read.reference_end
                                        break
                                in_for_end.close()
                                if next_end == 0:
                                    #logging.warning(' Next alignment not find, Query id: %s will be removed',read.query_name) 
                                    removed_count+=1
                                    removed_file_f.write(read)
                                    continue
                                
                                if read.next_reference_name == read.reference_name:
                                    start_next = int(read.next_reference_start)
                                    end_next = int(next_end)
                                    start_mapping_next = filter(lambda m: m[1] <= start_next and start_next <= m[2], mappings)
                                    end_mapping_next = filter(lambda m: m[1] <= end_next and end_next <= m[2], mappings)
                                    if len(start_mapping_next)!=1 or len(end_mapping_next)!=1:
                                        removed_count+=1
                                        removed_file_f.write(read)
                                        continue
                                    else:
                                        read_out.next_reference_start = int(start_next - start_mapping_next[0][1] + start_mapping_next[0][4])
                                        read_out.next_reference_id = header_new_ID_dict[read.next_reference_name]
                                else:                                
                                    mappings_next = self.alignment_dict[read.next_reference_name]
                                    start_next = int(read.next_reference_start)
                                    end_next = int(next_end)
                                    start_mapping_next = filter(lambda m: m[1] <= start_next and start_next <= m[2], mappings_next)
                                    end_mapping_next = filter(lambda m: m[1] <= end_next and end_next <= m[2], mappings_next)
                                    if len(start_mapping_next)!=1 or len(end_mapping_next)!=1:
                                        removed_count+=1
                                        removed_file_f.write(read)
                                        continue
                                    else:
                                        read_out.next_reference_start = int(start_next - start_mapping_next[0][1] + start_mapping_next[0][4])
                                        read_out.next_reference_id = header_new_ID_dict[read.next_reference_name]
                            else:
                                in_for_end = pysam.AlignmentFile(self.bam_file, 'rb')
                                check_next=in_for_end.count(read.next_reference_name,read.next_reference_start,read.next_reference_start+1) 
                                in_for_end.close()
                                if check_next!=0:                                                            
                                    read_out.next_reference_id = header_new_ID_dict[read.next_reference_name]
                                    read_out.next_reference_start = read.next_reference_start
                                else:
                                    removed_count+=1
                                    removed_file_f.write(read)
                                    continue
                        elif read.next_reference_id == -1:
                            read_out.next_reference_id = read.next_reference_id
                            read_out.next_reference_start = read.next_reference_start
                                                    
                        read_out.template_length = read.template_length                        
                        read_out.query_sequence = read.query_sequence
                        read_out.query_qualities = read.query_qualities
                        new_tag = []
                        CC_tag = "="
                        mappings_diff = 0
                        CP_reference_name = read.reference_name
                        tagname_notfound = 0
                        for tag_old in read.tags:   
                            if tag_old[0].startswith('CC'):
                                CC_tag = tag_old[1]
                                if tag_old[1] != "=":
                                    mappings_tag = self.alignment_dict[tag_old[1]]
                                    if len(mappings_tag)!=0:
                                        new_tag.append((tag_old[0],mappings_tag[0][3]))
                                        CP_reference_name = mappings_tag[0][3]
                                    else:
                                        removed_count+=1
                                        removed_file_f.write(read)
                                        tagname_notfound = 1
                                        break
                                else:
                                    mappings_tag = self.alignment_dict[read.reference_name]
                                    new_tag.append((tag_old[0],tag_old[1]))
                                    
                            elif tag_old[0].startswith('CP'):
                                CP_end = 0
                                for mappings_line in mappings_tag:
                                    if mappings_line[1] != mappings_line[4] or mappings_line[2] != mappings_line[5]:
                                        mappings_diff = 0
                                        break
                                    elif len(mappings_tag)!=1:
                                        mappings_diff = 0
                                        break                   
                                    elif mappings_line[1] == mappings_line[4] and mappings_line[2] == mappings_line[5]:
                                        mappings_diff = 1
                                if mappings_diff ==0:
                                    try:
                                        in_for_CP = pysam.AlignmentFile(self.bam_file, 'rb')
                                        for CP_read in in_for_CP(CP_reference_name,tag_old[1],tag_old[1]+1):
                                            if tag_old[1] == CP_read.reference_start:
                                                CP_end = CP_read.reference_end
                                                break
                                        in_for_CP.close()
                                    except:
                                        CP_end = 0
                                        logging.warning('%s, CP tag reference alignment not find! This alignment will be removed.',read.query_name)
                                        removed_count+=1
                                        removed_file_f.write(read)
                                        tagname_notfound = 1
                                        break  
                                 
                                    if CC_tag == "=":
                                        start_mapping_cp = filter(lambda m: m[1] <= tag_old[1] and tag_old[1] <= m[2], mappings)
                                        end_mapping_cp = filter(lambda m: m[1] <= CP_end and CP_end <= m[2], mappings)
                                    elif CC_tag != "=":
                                        mappings_tag = self.alignment_dict[CC_tag]
                                        start_mapping_cp = filter(lambda m: m[1] <= tag_old[1] and tag_old[1] <= m[2], mappings_tag)
                                        end_mapping_cp = filter(lambda m: m[1] <= CP_end and CP_end <= m[2], mappings_tag)
                                    if len(start_mapping_cp)!=1 and len(end_mapping_cp)!=1:                      
                                        #logging.info('%s, next hit be removed!!!',read.query_name)
                                        removed_count+=1
                                        removed_file_f.write(read)
                                        tagname_notfound = 1
                                        break
                                    else:
                                        CP_start = int(tag_old[1] - start_mapping_cp[0][1] + start_mapping_cp[0][4])
                                        new_tag.append((tag_old[0],CP_start))
                               
                                else:
                                    new_tag.append((tag_old[0],tag_old[1]))
        
                            elif Program_ID == "TopHat" or Program_ID == "Bowtie":
                                if tag_old[0] not in Custom_tag:                    
                                    if tag_old[0].startswith('X'):
                                        logging.warning('  New custom tag: %s detected! This tag will not be updated!',tag_old[0])
                                        new_tag.append((tag_old[0],tag_old[1]))
                                    else:
                                        new_tag.append((tag_old[0],tag_old[1]))
                                else:
                                    new_tag.append((tag_old[0],tag_old[1]))
                            elif Program_ID == "STAR":
                                if tag_old[0] not in Custom_tag:
                                    if tag_old[0].startswith('X'):
                                        logging.warning('  New custom tag: %s detected! This tag will not be updated!',tag_old[0])
                                        new_tag.append((tag_old[0],tag_old[1]))
                                    else:
                                        new_tag.append((tag_old[0],tag_old[1]))
                                elif tag_old[0] == "JI":
                                    logging.warning('  Detect JI tag! This tag will not be updated!')
                                    new_tag.append((tag_old[0],tag_old[1]))
                            elif Program_ID == "BWA":
                                if tag_old[0] not in Custom_tag:
                                    if tag_old[0].startswith('X'):
                                        logging.warning('  New custom tag: %s detected! This tag will not be updated!',tag_old[0])
                                        new_tag.append((tag_old[0],tag_old[1]))
                                    else:
                                        new_tag.append((tag_old[0],tag_old[1]))
                        
                                elif tag_old[0] == "XA":
                                    logging.warning('  Detect XA tag! This tag will not be updated!')
                                    new_tag.append((tag_old[0],tag_old[1]))
                            else:
                                new_tag.append((tag_old[0],tag_old[1]))

                        if tagname_notfound != 0:
                            continue
                                
                        read_out.tags= new_tag

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
if __name__ == '__main__':
    import sys
    import argparse
    from textwrap import dedent
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter, description=dedent("""\
    Update the reference sequence name of the alignment and coordinates of a bam file using an alignment file generated by the fasta_diff.py program.
    Updated features are written to a new file with '_updated'(default) appended to the original bam file name.
    Feature that can not be updated, due to the id being removed completely or the feature contains regions that
    are removed or replaced with Ns, are written to a new file with '_removed'(default) appended to the original bam file name.

    Example:
        fasta_diff.py old.fa new.fa | %(prog)s a.bam b.bam c.bam
    """))
    parser.add_argument('bam_files', metavar='bam_FILE', nargs='+', type=str, help='List one or more bam files to be updated')
    parser.add_argument('-a', '--alignment_file', type=argparse.FileType('rb'), default=sys.stdin,
                        help='The alignment file generated by fasta_diff.py, a TSV file with 6 columns: old_id, old_start, old_end, new_id, new_start, new_end (default: STDIN)')
    parser.add_argument('-u', '--updated_postfix', default='_updated',
                        help='The filename postfix for updated features (default: "_updated")')
    parser.add_argument('-r', '--removed_postfix', default='_removed',
                        help='The filename postfix for removed features (default: "_removed")')
    parser.add_argument('-v', '--version', action='version', version='%(prog)s ' + __version__)


    args = parser.parse_args()
    Bam_updater = BamUpdater(args.alignment_file, args.updated_postfix, args.removed_postfix)
    for Bam_file in args.bam_files:
        Bam_updater.update(Bam_file)


 	
   
