from Bio import SeqIO, AlignIO, Seq
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align.Applications import MafftCommandline
from subprocess import call
import sys

def Translate(File, out, Code = 11):
	fasta = SeqIO.parse(File,"fasta")
	New_fasta = []
	for s in fasta:
		aas = s.seq.ungap("-").translate(table = Code)
		sequence_object = Seq(str(aas))
		record = SeqRecord(sequence_object,id=s.id, description="")
		New_fasta.append(record)
	SeqIO.write(New_fasta, out, "fasta")
def Backtranslate(ali, input_fasta, out):
	New_fasta = []
	for s in ali:
		aa_seq = str(s.seq)
		dna = SeqIO.parse(input_fasta,"fasta")
		for cds in dna:
			if cds.id == s.id:
				nt_seq = str(cds.seq.ungap("-"))
				c1 = 0 # codon start
				c2 = 3 # codon end
				chain = ""
				for aa in aa_seq:
					if aa != "-":
						codon = nt_seq[c1:c2]
						c1 += 3
						c2 += 3
						chain += codon
					else:
						chain += aa*3
				sequence_object = Seq(chain.lower())
				record = SeqRecord(sequence_object,id=cds.id, description="")
				New_fasta.append(record)
	SeqIO.write(New_fasta,out,"fasta")

#{output.Merged} {output.tempaaMerged} {output.tempaaAligned} {output.MSA}
Input_fasta = sys.argv[1]
AA_fasta = sys.argv[2]
AA_aligned = sys.argv[3]
nucl_aligned = sys.argv[4]

print("AA-translation")
Translate(Input_fasta, AA_fasta)
print("Alignment")
mafft_cline = call("ml MAFFT; mafft --auto {u} > {a}".format(u=AA_fasta, a=AA_aligned), shell=True)  #MafftCommandline(mafft_exe, input=Unaligned)

align = AlignIO.read(AA_aligned, "fasta")
print("Backtranslation")
Backtranslate(align, Input_fasta, nucl_aligned)


