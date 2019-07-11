# Filter_InternalPriming

Requires bedtools and samtools installed in run environment. Developed and tested with [SAMtools](http://samtools.sourceforge.net/) v1.7 and [bedtools](https://bedtools.readthedocs.io/en/latest/) v2.27.1-1-gb87c465

From an aligned .bam file, resulting from single-end or paired-end sequencing, returns a .bam file without the reads that have their downstream alignment sequence (the following nucleotides beyond the read) matching the sequence of the adapter given.
For paired-end sequencing labraries, the downstream sequence corresponds to the following nucleotides of read 2.

Parameters:

-f, --filepath	input file

-s, --strategy	single/paired (use single for NET-seq data already in single nucleotide format)

-a, --adapter	TGG.. N's in the sequence can be represented by "."

-g, --genome	reference genome .fasta

Designed for NET-seq data and described as part of an analysis pipeline in *insert reference*

Usage examples:

python Filter_InternalPriming.py -f /alignments/merged_mNET_Long_S5P_rep1_unique.bam -s single -a TGG -g /genomes/human/hg38/GRCh38.primary.genome.fa

python Filter_InternalPriming.py -f /alignments/mNET_Long_S5P_rep1_unique_sorted.bam -s paired -a .GGA -g /genomes/human/hg38/GRCh38.primary.genome.fa


### It operates in the following manner:

For single-end reads or merged paired-end reads it will convert the input .bam alignment file into a .bed file.
/
For paired-end reads it will first extract the second read from the pair, add them to a new .bam file which only then will be converted to a .bed file.

Iterate through the resulting .bed file to add the genome coordinates downstream of the 3'OH position that correspond to the length of the given adapter sequence. In other words, if the given adapter sequence has 5 nucleotides and the read is 150 nucleotides long then the resulting .bed entry will have a length of 155.
3'OH position will correspond to the last nucleotide of the read for single/merged reads. And for the second read of the pair it will correspond to the first nucleotide of the read.

Next step is extracting nucleotide sequence for each .bed entry (converting the .bed file into a .fasta file). All entries where the 3'OH downstream sequence match the sequence of the given adapter will be discarded. From the other ones we save the read IDs into a .txt file. A .txt file that now has the read IDs of the reads that do **not** result from internal priming. 

Extraction of the reads that have matching IDs with the ones in the internal priming-free .txt file from the original alignment file.







