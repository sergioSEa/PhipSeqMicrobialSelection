from Bio import SeqIO
import sys
#/scratch/umcg-sandreusanchez/Immuno_selection/2_AlignmentAndPhylogenetics/Software/Tools_2
records = SeqIO.parse(sys.argv[1], "fasta")
#for r in records:
#	print(r)

with open(sys.argv[2], "w") as output_handle:
	SeqIO.write(records,  output_handle, "phylip")
