#!/usr/bin/env python


# Import libraries
from os import system
import argparse


############################################################
#                                                          #
#                         Functions                        #
#                                                          #
############################################################

def getBEDwindows_PlusXnuc(reads,added_nuc,paired):

	peak_list=[]
	output=open(filename+"_readCoords.bed","w")

	with open(reads) as exonlist:
		for exon in exonlist:
			read=exon.strip("\n\r").split("\t")
			if paired:
				if read[5] == "+": # "+" if paired end and info is in R2 /// "-" if single end
					if int(read[1])-added_nuc < 0 :
						continue
					output.write(read[0]+"\t"+str(int(read[1])-added_nuc)+"\t"+read[2]+"\t"+read[3]+"\t"+read[4]+"\t"+"-"+"\n")
				else:
					output.write(read[0]+"\t"+read[1]+"\t"+str(int(read[2])+added_nuc)+"\t"+read[3]+"\t"+read[4]+"\t"+"+"+"\n")
			else:
				if read[5] == "-":
					if int(read[1])-added_nuc < 0 :
						continue
					output.write(read[0]+"\t"+str(int(read[1])-added_nuc)+"\t"+read[2]+"\t"+read[3]+"\t"+read[4]+"\t"+"-"+"\n")
				else:
					output.write(read[0]+"\t"+read[1]+"\t"+str(int(read[2])+added_nuc)+"\t"+read[3]+"\t"+read[4]+"\t"+"+"+"\n")

	output.close()


def parseFasta4_AdapterSeq(exons):
	import re
	regex1 = re.compile(adapter_seq.upper())
	regex2 = re.compile(adapter_seq.lower())
	exon_list=[]
	peak_list=[]

	output=open(filename+"_woIPri_ID.txt","w")

	with open(exons) as exonlist:
			for exon in exonlist:
				fasta=exon.strip("\n\r").split("\t")
				if fasta[0].startswith(">"):
					ID=fasta
				else:
					misP=fasta[0][-len(adapter_seq):]
					matches1 = re.match(regex1, misP)
					matches2 = re.match(regex2, misP)
					if matches1:
						continue
					elif matches2:
						continue
					else:
						readID=ID[0][:-2].strip(">")
						readID=readID.strip("(")
						readID=readID.split("/")[0] # parsing readID
						output.write(readID+"\t"+"\n")

	fastas=[]
	c=-1
	rev_seq=True
	tgg_counter=0

############################################################
#                                                          #
#                           MAIN                           #
#                                                          #
############################################################

# Input arguments
parser = argparse.ArgumentParser(description="""Returns a .bam file
	without the reads where the downstream sequence matches the
	given adapter sequence.
	run example:
	python Filter_InternalPriming.py -f /alignments/merged_mNET_Long_S5P_rep1_unique.bam -s single -a Tgg.. -g /genomes/human/hg38/GRCh38.primary.genome.fa""")


parser.add_argument('-f', '--filepath', dest='filepath', nargs='?',
                   help='One .bam file.')
parser.add_argument('-s', '--strategy', dest='strategy', nargs='?',
                   help="""Input paired/single accordingly. Defaults to single""",
                   default="'single'")
parser.add_argument('-g', '--genome', dest='genome', nargs='?',
                   help="""Input .fasta from reference genome.""")
parser.add_argument('-a', '--adapter', dest='adapter', nargs='?',
                   help="""Input adapter sequence or simply the first few nt of the adapter sequence.
                   Can also accomodate for regular expression. Ex: TGG..""")

args = parser.parse_args()
file = args.filepath
strategy = args.strategy
adapter_seq = args.adapter
fasta = args.genome

# Defining filename and setting required variable for getBEDwindows_PlusXnuc function.
filename=str(file).split("/")[-1].rstrip(".bam")
coord_strat=False

# Extracring .bam file header
system("samtools view -H "+file+" > "+filename+"_header.sam")
header=str(filename)+"_header.sam"

if strategy == "paired": #Get read 2 in paired-end datasets
	coord_strat=True
	system("samtools view "+file+" | awk '($2 ~ /147|163/)' > "+filename+"_tempR2.sam")
	system("cat "+header+" "+filename+"_tempR2.sam > "+filename+"_tempR2_H.sam")
	system("samtools view -Sb -h "+filename+"_tempR2_H.sam > "+filename+".bam")
	system("rm -f "+filename+"_tempR2.sam")
	system("rm -f "+filename+"_tempR2_H.sam")
	file=filename+".bam"

# Convert .bam file into .bed file
system("bamToBed -i "+file+" > "+filename+".bed")


# Get.bed coordinates downstream of the 3'Position. Each .bed length = size of the adapter sequence given.
getBEDwindows_PlusXnuc(filename+".bed",len(adapter_seq),coord_strat)
system("rm -f "+filename+".bed")

# Get .fasta sequence for the filename+"_readCoords.bed" files
system("bedtools getfasta -fi "+fasta+" -bed "+filename+"_readCoords.bed -fo "+filename+".fasta -s -name")
system("rm -f "+filename+"_readCoords.bed")

# Get the read ID's of the reads WITHOUT the adapter sequence immediately downstream of the Pol II position.
parseFasta4_AdapterSeq(filename+".fasta")
system("rm -f "+filename+".fasta")

# From the original file get the reads with the ID matching the list of IDs WITHOUT the adapter sequence immediately downstream of the Pol II position.
# WARNING - Long-ish run time
system("samtools view "+file+" | fgrep -f "+filename+"_woIPri_ID.txt > "+filename+".sam")
system("cat "+header+" "+filename+".sam > "+filename+"_H.sam")
system("samtools view -Sb -h "+filename+"_H.sam > "+filename+"_woIPri.bam")
system("samtools sort -o "+filename+"_woIPri_sorted.bam "+filename+"_woIPri.bam")

#delete temp files used
system("rm -f "+filename+".sam")
system("rm -f "+filename+"_H.sam")
system("rm -f "+filename+"_woIPri.bam")
system("rm -f "+filename+"_woIPri_ID.txt")
system("rm -f "+header)
system("rm -f "+filename+".bam")


print "Finished!"