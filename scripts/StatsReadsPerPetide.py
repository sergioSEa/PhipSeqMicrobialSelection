from pathlib import Path
import pandas as pd


Location = Path("4_Prepared_reads")

Dic_Peptide = {}
#counter = 0
for File in Location.glob("*.fq"):
        F = File.name
        if "agilent" in F: Spl = "agilent"
        else: Spl ="twist"
        Participant = F.split("_"+Spl)[0]
        Peptide = (Spl + F.split("_"+Spl)[1]).rstrip(".fq")
        print(Peptide, Participant)
        if Peptide not in Dic_Peptide:
                Dic_Peptide[Peptide] = {}
        n = 0
        with open(File) as O:
                for line in O: 
                        n+=1
        n = n/4
        Dic_Peptide[Peptide][Participant] = n 
#       if counter == 5:
#               break
#       else: counter +=1
Output = pd.DataFrame.from_dict(Dic_Peptide)
Output.index.name = 'ID'
Output.reset_index(inplace=True)

Output.to_csv("Results/PeptideAlignmentStats.tsv", sep = "\t", index=False)
