# coordinates_conversion

[![Build Status](https://travis-ci.org/NAL-i5K/coordinates_conversion.svg?branch=master)](https://travis-ci.org/NAL-i5K/coordinates_conversion)

Conversion programs that use the output from fasta_diff.py to convert reference sequence IDs and coordinates in Gff3, bam, bed, or bedgraph file formats. Main contributors are [Han Lin](https://github.com/hotdogee) (original development) and interns of i5k workspace.

Scripts to convert reference sequence IDs and coordinates in different file formats.
* fasta_diff.py
    - [wiki page](https://github.com/NAL-i5K/coordinates_conversion/wiki/fasta_diff.py)
* update_gff.py
    - [wiki page](https://github.com/NAL-i5K/coordinates_conversion/wiki/update-gff)
* update_bam.py
    - [wiki page](https://github.com/NAL-i5K/coordinates_conversion/wiki/update-bam)
* update_bed.py
    - [wiki page](https://github.com/NAL-i5K/coordinates_conversion/wiki/update-bed)
* update_bedgraph.py
     - [wiki page](https://github.com/NAL-i5K/coordinates_conversion/wiki/update-bedgraph)
* update_vcf.py
     - [wiki_page](https://github.com/NAL-i5K/coordinates_conversion/wiki/update-vcf)

## Quick start
1. Run fasta_diff.py    
  * Compares two very similar FASTA files and outputs coordinate mappings using a multi stage algorithm:  
  * Stage 1: Find 100% matches  
  * Stage 2: Find 100% substrings, where the full length of a new sequence can be found as a substring of a old sequence  
  * Stage 3: Find cases where part of the sequence was converted into Ns  
  * Stage 4: Find cases where a old sequence is split into two or more new sequences
  * Outputs (match.tsv) the 6 columns as tab-separated values: old_id, old_start, old_end, new_id, new_start, new_end

    `fasta_diff.py old.fa new.fa > match.tsv`

2. Select a conversion script that matches your file format  
  * [Gff3 format](https://github.com/The-Sequence-Ontology/Specifications/blob/master/gff3.md): update_gff.py
  * [Bam format](http://samtools.github.io/hts-specs/SAMv1.pdf): bam_update.py
  * [Bed format](https://genome.ucsc.edu/FAQ/FAQformat#format1): bed_update.py
  * [Bedgraph format](https://genome.ucsc.edu/goldenpath/help/bedgraph.html): bedgraph_update.py

3. Run conversion script:
  * update_gff.py  

    `update_gff.py –a match.tsv a.gff b.gff c.gff`

  * update_bam.py  
    * The following programs need to be installed before running this program:
      * [pysam](http://pysam.readthedocs.io/en/latest/index.html)

        `pip install pysam`

      * [samtools](http://samtools.sourceforge.net/)
    * If you have a bam file without a corresponding index file (.bai), you can generate one using:  

    `samtools index aln.bam`

    * Then use update_bam.py to convert your bam files

    `update_bam.py –a match.tsv a.bam b.bam c.bam`

  * update_bed.py  

    `update_bed.py –a match.tsv a.bed b.bed c.bed`

  * update_bedgraph.py  

    `update_bedgraph.py –a match.tsv a.bedgraph b.bedgraph c.bedgraph`
    
  * update_vcf.py  

    `update_vcf.py –a match.tsv -ref new_assembly.fa a.vcf b.vcf c.vcf`
