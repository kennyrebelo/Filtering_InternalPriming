# Filter_InternalPriming

Requires bedtools and samtools installed in run environment. Developed and tested with [SAMtools](http://samtools.sourceforge.net/) v1.7 and [bedtools](https://bedtools.readthedocs.io/en/latest/) v2.27.1-1-gb87c465

From an aligned .bam file, resulting from single-end or paired-end sequencing, returns a .bam file without the reads that have their downstream alignment sequence (the following nucleotides beyond the read) matching the sequence of the adapter given.
For paired-end sequencing labraries, the downstream sequence corresponds to the following nucleotides of read 2.

Parameters:

-f, --filepath	input file

-s, --strategy	single/paired (use single for NET-seq data already in single nucleotide format)

-a, --adapter	TGG.. (accepts regular expressions)

-g, --genome	reference genome .fasta

Designed for NET-seq data and described as part of an analysis pipeline in *insert reference*

Usage examples:

python Filter_InternalPriming.py -f /alignments/merged_mNET_Long_S5P_rep1_unique.bam -s single -a TGG -g /genomes/human/hg38/GRCh38.primary.genome.fa

python Filter_InternalPriming.py -f /alignments/mNET_Long_S5P_rep1_unique_sorted.bam -s paired -a .GGA -g /genomes/human/hg38/GRCh38.primary.genome.fa
