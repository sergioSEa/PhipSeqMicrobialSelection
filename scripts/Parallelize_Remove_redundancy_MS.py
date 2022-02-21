from pathlib import Path
from subprocess import call
import pandas as pd
from Bio import AlignIO, SeqIO, Seq
from Bio.Align import AlignInfo
import numpy as np
import sys
np.set_printoptions(suppress=True,linewidth=np.nan)
def hamming_distance(s1, s2):
        """Return the Hamming distance between equal-length sequences"""
        if len(s1) != len(s2):
                raise ValueError("Undefined for sequences of unequal length")
        total_d = 0
        total = 0
        s1 = str(s1) ; s2 = str(s2)
        for position in range(0,len(s1)):
                ch1 = s1[position]
                ch2 = s2[position]
                if ch1 != "-" and ch2 != "-":
                        total += 1
                        if ch1 != ch2: total_d +=1
        if total == 0: return(np.nan)           
        #M =  sum(ch1 != ch2 for ch1, ch2 in zip(s1, s2) if ch1 != "-" and ch2 != "-")
        return(total_d/total)
def Consensus(Seq1, Seq2):
	Cons = ""
	for i in range(len(Seq1)):
		Add = Seq1[i]
		if Seq1[i] == "-": Add = Seq2[i]		
		Cons += Add
	return(Cons)
def Merge_seqs(Seq1, Seq2):
        Al = AlignIO.MultipleSeqAlignment([])
        Al.append( SeqIO.SeqRecord("item1", str(Seq1)))
        Al.append( SeqIO.SeqRecord("item2", str(Seq2)))
        print(Al)
        summary = AlignInfo.SummaryInfo(Al)
        
        consensus_seq = summary.gap_consensus(threshold = 0.5,ambiguous = "-")
        exit(consensus_seq)
        return(str(consensus_seq))
def Remove_redundancy_MSA(alignment, ID_p):
        Done = []
        Added = []
        Final_alignment = AlignIO.MultipleSeqAlignment([])
        Sequences_Added = 0
        for Sequence1 in alignment:
                Seq1 = Sequence1.seq
                for Sequence2 in alignment:
                        Seq2 = Sequence2.seq
                        ID = "~".join(sorted([Sequence1.id, Sequence2.id]))
                        Hamming = hamming_distance(Seq1, Seq2)
                        if Sequence2.id in Added or Sequence1.id in Added: continue
                        if ID in Done: continue
                        else: Done.append(ID)
                        if Hamming == 0:
                                Seq1 = Consensus(Seq1,Seq2)
                                if Sequence2.id != Sequence1.id:
                                        Added.append(Sequence2.id)
                        else: continue
                if Sequence1.id in Added: continue
                Sequences_Added += 1
                Name = ID_p + ".Ortholog" + str(Sequences_Added)
                Add_seq = SeqIO.SeqRecord(Seq.Seq(Seq1), id = Name, description="")
         
                Final_alignment.append(Add_seq)
                #except: print(Final_alignment) ; print(Add_seq) ; exit()
        return(Final_alignment)


MSA = sys.argv[1]
Output = sys.argv[2]
try :alignment = AlignIO.read(open(MSA), "fasta")
except: Path(Output).touch() ; exit()
ID = sys.argv[3] #Path(MSA).stem.split("_")[0]
#Peptide = sys.argv[3] #"_".join(Path(MSA).stem.split("_")[1:])


new_alignment = Remove_redundancy_MSA(alignment, ID)

Output = sys.argv[2]  #"Collapse_Alignments/" + MSA.name
AlignIO.write(new_alignment,  Output, "fasta")

