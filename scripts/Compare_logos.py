from pathlib import Path
from subprocess import call
import pandas as pd
import sys

def Get_presence_profile():
	#First, get information on each sample, whetather it has or has no antibody
	Sample_translation = {}
	with open("../Data/Translate_MGS_names.txt") as F:
		for line in F:
			l = line.rstrip().split()
			Sample_translation[l[1]] = l[0]
	Sample_presence = {}
	print("Reading AB")
	AB = pd.read_csv("../Data/Matrix_selection.tsv", sep ="\t")
	print("Iterate through samples")
	for sample in AB["ID"]:
		#print(sample, sample in AB.loc[:,"ID"].to_list(), sample in Sample_translation)
		try:
			Name = Sample_translation[sample]
			Profile_person = AB[AB.loc[:,"ID"] == sample]
		except: continue
	
		for ab in Profile_person.columns:
			if ab == "ID": continue
			P = Profile_person.loc[:,ab]
			if Name not in Sample_presence: Sample_presence[Name] = {}
			Sample_presence[Name][ab] = int(P)
	return(Sample_presence)


def main():
	from Bio import AlignIO
	Sample_presence = Get_presence_profile()
	F = sys.argv[1]	
	Ab = F.stem
	alignment = AlignIO.read(open(F), "fasta")
	No_presence = AlignIO.MultipleSeqAlignment([])
	Presence = AlignIO.MultipleSeqAlignment([])
	for record in alignment:
		Name = record.split(".Or")[0]
		if Name not in Sample_presence: continue
		P_N = Sample_presence[Name][Ab]
		if P_N == 0:
			No_presence.append( SeqIO.SeqRecord(str(record.seq), id =  record.id, comments="") )
		elif P_N == 1:
			Presence.append( SeqIO.SeqRecord(str(record.seq), id =  record.id, comments="") )
		else: continue
	Temp_no = "temp/{A}_nopresent.fa".format(A=Ab)
	Temp_yes = "temp/{A}_present.fa".format(A=Ab)
	AlignIO.write(No_presence,  Temp_no, "fasta")
	AlignIO.write(Presence,  Temp_yes, "fasta")
		
	#Make LOGOs
	Output_no = "results/logo/{A}_noPresent.png".format(A=Ab)
	Output_yes = "results/logo/{A}_Present.png".format(A=Ab)
	#Call command
	call("bash scripts/Logos.sh {MSA} {OUT}".format(MSA=Temp_no, OUT=Output_no), shell=True)
	call("bash scripts/Logos.sh {MSA} {OUT}".format(MSA=Temp_yes, OUT=Output_yes), shell=True)
		
		
if __name__ == '__main__':
    main()
