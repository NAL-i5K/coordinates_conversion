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
* fasta_diff.py
Compares two very similar FASTA files and outputs coordinate mappings using a multi stage algorithm:
Stage 1: Find 100% matches
Stage 2: Find 100% substrings, where the full length of a new sequence can be found as a substring of a old sequence
Stage 3: Find cases where part of the sequence was converted into Ns

<code>fasta_diff.py old.fa new.fa > match.tsv</code>
