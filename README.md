## coordinates_conversion
Conversion programs that use the output from fasta_diff.py to convert reference sequence IDs and coordinates in Gff3, bam, bed, or bedgraph file formats. Main contributors are [Han Lin](https://github.com/hotdogee) (original development) and [LiMei Chiang](https://github.com/dytk2134) (continued development).

## bin/
Scripts to convert reference sequence IDs and coordinates in different file formats.
* fasta_diff.py
    - [wiki page](https://github.com/NAL-i5K/coordinates_conversion/wiki/fasta_diff.py)
* update_gff.py
    - [wiki page](https://github.com/NAL-i5K/coordinates_conversion/wiki/update-gff)
* bam_update.py
    - [wiki page](https://github.com/NAL-i5K/coordinates_conversion/wiki/update-bam)
* bed_update.py
    - [wiki page](https://github.com/NAL-i5K/coordinates_conversion/wiki/update-bed)
* bedgraph_update.py
     - [wiki page](https://github.com/NAL-i5K/coordinates_conversion/wiki/update-bedgraph)
* update_vcf.py
     - [wiki_page](https://github.com/NAL-i5K/coordinates_conversion/wiki/update-vcf)

## Quick start
1. Run fasta_diff.py    
  * Compares two very similar FASTA files and outputs coordinate mappings using a multi stage algorithm:  
  * Stage 1: Find 100% matches  
  * Stage 2: Find 100% substrings, where the full length of a new sequence can be found as a substring of a old sequence  
  * Stage 3: Find cases where part of the sequence was converted into Ns  
  * Outputs (match.tsv) the 6 columns as tab-separated values: old_id, old_start, old_end, new_id, new_start, new_end

    <code>fasta_diff.py old.fa new.fa > match.tsv</code>

2. Select a conversion script that matches your file format  
  * [Gff3 format](https://github.com/The-Sequence-Ontology/Specifications/blob/master/gff3.md): update_gff.py
  * [Bam format](http://samtools.github.io/hts-specs/SAMv1.pdf): bam_update.py
  * [Bed format](https://genome.ucsc.edu/FAQ/FAQformat#format1): bed_update.py
  * [Bedgraph format](https://genome.ucsc.edu/goldenpath/help/bedgraph.html): bedgraph_update.py
    
3. Run conversion script:
  * update_gff.py  

    <code> update_gff.py –a match.tsv a.gff b.gff c.gff </code>  

  * bam_update.py  
    * The following programs need to be installed before running this program:
      * [pysam](http://pysam.readthedocs.io/en/latest/index.html)

        <code>pip install pysam</code>

      * [samtools](http://samtools.sourceforge.net/)
    * If you have a bam file without a corresponding index file (.bai), you can generate one using:  

    <code> samtools index aln.bam </code>  
    * Then use bam_update.py to convert your bam files

    <code> bam_update.py –a match.tsv a.bam b.bam c.bam </code>  

  * bed_update.py  

    <code> bed_update.py –a match.tsv a.bed b.bed c.bed </code>  

  * bedgraph_update.py  

    <code> bed_update.py –a match.tsv a.bedgraph b.bedgraph c.bedgraph </code>  
    
  * update_vcf.py  

    <code> update_vcf.py –a match.tsv -ref new_assembly.fa a.vcf b.vcf c.vcf </code> 
