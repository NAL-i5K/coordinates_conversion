## coordinates_conversion
Conversion programs that use the output from fasta_diff.py to convert coordinates and IDs in different format files.

## bin/
Script for convert coordinates and IDs in different format files.
* bedgraph_update.py
* fasta_diff.py
* update_gff.py
* bam_update.py
    - [wiki page] (https://github.com/NAL-i5K/coordinates_conversion/wiki/update-bam)
* bed_update.py

## Quick start
### Preparation  
1. Get the outputs from fasta_diff.py    
* fasta_diff.py  
Compares two very similar FASTA files and outputs coordinate mappings using a multi stage algorithm:  
Stage 1: Find 100% matches  
Stage 2: Find 100% substrings, where the full length of a new sequence can be found as a substring of a old sequence  
Stage 3: Find cases where part of the sequence was converted into Ns  

    <code>fasta_diff.py old.fa new.fa > match.tsv</code>

Outputs(match.tsv) the 6 columns as tab-separated values: old_id, old_start, old_end, new_id, new_start, new_end  
2. Select a script that matches your file format  
*If you want to convert coordinates and IDs in BAM files, there are some tool should be install first*
    - [pysam] (http://pysam.readthedocs.io/en/latest/index.html)        
        <code>pip install pysam</code>
    - [samtools] (http://samtools.sourceforge.net/)
    
### coordinates and IDs convertion
