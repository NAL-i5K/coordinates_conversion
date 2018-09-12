# coordinates_conversion

[![Build Status](https://travis-ci.org/NAL-i5K/coordinates_conversion.svg?branch=master)](https://travis-ci.org/NAL-i5K/coordinates_conversion)

Conversion programs that use the output from fasta_diff to convert reference sequence IDs and coordinates in Gff3, bam, bed, or bedgraph file formats. Main contributors are [Han Lin](https://github.com/hotdogee) (original development) and interns of i5k workspace.

## Prerequisite

- Python 2.7
- samtools (optional, only for SAM/BAM related scripts)

## Installation

`pip install git+https://github.com/NAL-i5K/coordinates_conversion.git`

## Features

Scripts to convert reference sequence IDs and coordinates in different file formats.

- fasta_diff
  - [wiki page](https://github.com/NAL-i5K/coordinates_conversion/wiki/fasta_diff)
- update_gff
  - [wiki page](https://github.com/NAL-i5K/coordinates_conversion/wiki/update_gff)
- update_bam
  - [wiki page](https://github.com/NAL-i5K/coordinates_conversion/wiki/update_bam)
- update_bed
  - [wiki page](https://github.com/NAL-i5K/coordinates_conversion/wiki/update_bed)
- update_bedgraph
  - [wiki page](https://github.com/NAL-i5K/coordinates_conversion/wiki/update_bedgraph)
- update_vcf
  - [wiki_page](https://github.com/NAL-i5K/coordinates_conversion/wiki/update_vcf)

## Quick start

1. Run `fasta_diff`
- Compares two very similar FASTA files and outputs coordinate mappings using a multi stage algorithm:
- Stage 1: Find 100% matches
- Stage 2: Find 100% substrings, where the full length of a new sequence can be found as a substring of a oldsequence
- Stage 3: Find cases where part of the sequence was converted into Ns
- Stage 4: Find cases where a old sequence is split into two or more new sequences
- Outputs (match.tsv) the 6 columns as tab-separated values: old_id, old_start, old_end, new_id, new_start, new_end

  `fasta_diff example_file/old.fa example_file/new.fa -o match.tsv -r report.txt`

2. Select a conversion script that matches your file format
- [Gff3 format](https://github.com/The-Sequence-Ontology/Specifications/blob/master/gff3.md): update_gff
- [Bam format](http://samtools.github.io/hts-specs/SAMv1.pdf): bam_update
- [Bed format](https://genome.ucsc.edu/FAQ/FAQformat#format1): bed_update
- [Bedgraph format](https://genome.ucsc.edu/goldenpath/help/bedgraph.html): bedgraph_update

3. Run conversion script:
- update_gff

  `update_gff -a match.tsv example_file/example1.gff3 example_file/example2.gff3`

- update_bam
  - [samtools](http://samtools.sourceforge.net/) needs to be installed before running this program:
  - If you have a bam file without a corresponding index file (.bai), you can generate one using:

    `samtools index example_file/example.bam`

  - Then use update_bam to convert your bam files

    `update_bam -a match.tsv example_file/example.bam`

  - update_bed

    `update_bed -a match.tsv example_file/example.bed`

  - update_bedgraph

    `update_bedgraph -a match.tsv example_file/example.bedGraph`

  - update_vcf

    `update_vcf -a match.tsv example_file/example.vcf`
