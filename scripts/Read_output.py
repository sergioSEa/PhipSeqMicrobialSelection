import json
import sys
from pathlib import Path
import pandas as pd

Output_file = "7_Evo_models/agilent_55961_Output.txt"
Result = pd.DataFrame(columns = ["Peptide","LRT","P"] )
for Output_file in Path("7_Evo_models").glob("*_Output.txt"):
	Peptide = Output_file.name.split("_Output")[0]
	with open(Output_file) as f:
		Output = json.load(f)
	Test = Output['test results']
	Temp_result = pd.DataFrame( [[Peptide, Output['test results']["LRT"], Output['test results']["p-value"]] ], columns= ["Peptide","LRT","P"])
	Result = Result.append(Temp_result,ignore_index=True)
	#Evidence ratios only if Test_result < 0.05
	#print(Output.keys())
Result.to_csv("Results/Selection_test.tsv", sep = "\t", index_label=False, index=False)
