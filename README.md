## coordinates_conversion
Conversion programs that use the output from fasta_diff.py to convert coordinates and IDs in different format files.

## bin/
Script for convert coordinates and IDs in different format files.
* bedgraph_update.py
     - [wiki page] (https://github.com/NAL-i5K/coordinates_conversion/wiki/update-bedgraph)
* fasta_diff.py
    - [wiki page] (https://github.com/NAL-i5K/coordinates_conversion/wiki/fasta_diff.py)
* update_gff.py
    - [wiki page] (https://github.com/NAL-i5K/coordinates_conversion/wiki/update-gff)
* bam_update.py
    - [wiki page] (https://github.com/NAL-i5K/coordinates_conversion/wiki/update-bam)
* bed_update.py
    - [wiki page] (https://github.com/NAL-i5K/coordinates_conversion/wiki/update-bed)

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
* update_gff.py  
    <code> update_gff.py –a match.tsv a.gff b.gff c.gff </code>  

* bam_update.py  
    *If you have a bam file without a corresponding index, you can generate one using:*  
    <code> samtools index aln.bam </code>  
    *Then using bam_update.py to convert your bam files*    
    <code> bam_update.py –a match.tsv a.bam b.bam c.bam </code>  
    
* bed_update.py  
    <code> bed_update.py –a match.tsv a.bed b.bed c.bed </code>  

* bedgraph_update.py  
    <code> bed_update.py –a match.tsv a.bedgraph b.bedgraph c.bedgraph </code>  
