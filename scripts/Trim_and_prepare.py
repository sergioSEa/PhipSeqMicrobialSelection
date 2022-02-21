from Bio import SeqIO
from pathlib import Path
from subprocess import call
import sys
#ml Biopython/1.78-foss-2020a-Python-3.8.2





def Trim_and_prepare(Info_file, Fastq_file, Ref ,Out):
	#if ".fa" in Out: Out = Out.split(".fa")[0]
	Peptide_l = []
	for record in SeqIO.parse(Ref, "fasta"):
        	Peptide_l.append(str(record.id) )
	dic_info = {}
	peptides = []
	with open(Info_file) as I:
		for line in I:
			l = line.rstrip().split()
			if int(l[7]) > int(l[6]): direction ="same"
			else: direction = "complementary"
			Peptide = l[1]
			dic_info[l[0]] = [int(l[6]), int(l[7]), direction, Peptide]
			if Peptide not in peptides: peptides.append(Peptide)

	dic_results = {}
	for peptide in Peptide_l:
		dic_results[peptide] = { "FQ" : "", "FA" : "" }
	for record in SeqIO.parse(Fastq_file, "fastq"):
		Info =  dic_info[str(record.id)]
		Peptide = Info[-1]
		Q = record.format("fastq").split()[3]
		Seq = record.seq.strip()
		if Info[2] == "same":
			Seq = Seq[Info[0]-1:Info[1]]
			Q = Q[Info[0]-1:Info[1]]
		else:
			Seq = Seq[Info[1]-1:Info[0]]
			Seq = Seq.reverse_complement()
			Q =  Q[Info[1]-1:Info[0]][::-1]
		dic_results[Peptide]["FA"] += ">"+str(record.id)+"\n"+str(Seq) + "\n"
		dic_results[Peptide]["FQ"] += "@"+str(record.id)+"\n"+ str(Seq) +  "\n+\n" + Q +  "\n"

	for Peptide in dic_results:
		Out_fa = Out + "_" + Peptide + ".fa"
		Out_fq = Out + "_" + Peptide + ".fq"
		with open(Out_fa, "w") as O:
			O.write(dic_results[Peptide]["FA"])
		with open(Out_fq, "w") as O:
			O.write(dic_results[Peptide]["FQ"])


Trim_and_prepare(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4])
