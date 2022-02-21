#http://etetoolkit.org/docs/latest/tutorial/tutorial_adaptation.html
#conda activate /scratch/umcg-sandreusanchez/Immuno_selection/Snake_pipeline/ete3
from ete3 import EvolTree
from Compare_logos import Get_presence_profile
from pathlib import Path
from subprocess import call
from Bio import SeqIO, AlignIO
import sys

Example = False
if Example == False:
	Tree_file = sys.argv[2] #"6_Alignment_Peptides/agilent_55961.fa.treefile"
	Alignment_file = sys.argv[1] #"6_Alignment_Peptides/agilent_55961.fa"

	AB = Path(Alignment_file).stem
	codeml_Path = "/scratch/umcg-sandreusanchez/Immuno_selection/Snake_pipeline/ete3/bin/ete3_apps/bin/codeml" #Not used, I had to chenge .bashrc to include the path

	tree = EvolTree("\n".join( open(Tree_file).readlines()))
	tree.link_to_alignment("\n".join( open(Alignment_file).readlines()))

	MSA_temp = "temp/MSA_filtered_{peptide}".format(peptide = AB)
	Tree_temp = "temp/Tree_FilteredAndMarked_{peptide}".format(peptide = AB)
	Out = '/scratch/umcg-sandreusanchez/Immuno_selection/Snake_pipeline/7_Evo_models/{peptide}_Output.txt'.format(peptide = AB)
def Mark_tree(tree, marks,tree_out, character):
	To_mark = []
	[ To_mark.append(character) for i in marks ]
	tree.mark_tree(marks, To_mark)
	tree.write(outfile=tree_out)
	return(tree)
def Compute_selection_hyphy(tree, Align, marks, peptide, tree_out,Temp_align, OUT):
	tree = Mark_tree(tree, marks,tree_out, "{test}")
	command1 = "sed -i 's/ \#1/{test}/g' "+ tree_out  #Change #1 by {test}
	call(command1, shell=True)
	Ls = []
	for leaf in tree.get_leaves():
		L = str(leaf)[3:]
		Ls.append(L)
	temp_align = ""
	for line in Align:#.split():
		if line[0] == ">":
			if line[1:].rstrip() in Ls: Keep = True
			else: Keep = False
		if Keep == True: temp_align += line + "\n"
	with open(Temp_align, "w") as O: O.write(temp_align)
	alignment = AlignIO.read(open(Temp_align), format="fasta")
	codon_stop_array = ["TAG", "TGA", "TAA", "UGA", "UAA", "UAG"]
	
	temp_align=""
	for record in alignment:
		tempRecordSeq = list(record.seq)
		for index in range(0, len(record.seq), 3):
			codon = record.seq[index:index+3]
			if codon.upper() in codon_stop_array:
				tempRecordSeq[index:index+3] = "nnn"
		SEQ = "".join(tempRecordSeq)
		temp_align+= ">"+record.id +"\n"+ SEQ + "\n"
	with open(Temp_align, "w") as O: O.write(temp_align)



	command = """source activate /scratch/umcg-sandreusanchez/Immuno_selection/2_AlignmentAndPhylogenetics/Software/Tools_2;
hyphy busted --code Universal --alignment {algn} --tree {tree} --branches test --output {out} CPU=4""".format(algn= Temp_align , tree=tree_out , out= OUT )
	print(command)
	call(command, shell=True)

def Compute_Selection(tree, marks, peptide):
	tree.workdir = '/scratch/umcg-sandreusanchez/Immuno_selection/Snake_pipeline/7_Evo_models/'+peptide
	if not Path(tree.workdir): call("mkdir "+tree.workdir, shell=True)
	
	tree = Mark_tree(tree, marks, "#1")	
	

	tree.run_model('b_free') #Allow Marks to evolve under positive selction
	tree.run_model('M0')
	
	b_free = tree.get_evol_model('b_free')

	pval = tree.get_most_likely ('b_free', 'M0')	
	print(b_free)
	print(pval)
	with open(tree.workdir + "_Output.txt", "w") as O:
		O.write(str(b_free))
		O.write("\n PVAL: {pval} \n".format(pval=str(pval)))
def Make_example():
	tree = EvolTree("((Hylobates_lar,(Gorilla_gorilla,Pan_troglodytes)),Papio_cynocephalus);")
	align = '''>Hylobates_lar\nATGGCCAGGTACAGATGCTGCCGCAGCCAGAGCCGGAGCAGATGTTACCGCCAGAGCCGGAGCAGATGTTACCGCCAGAGGCAAAGCCAGAGTCGGAGCAGATGTTACCGCCAGAGCCAGAGCCGGAGCAGATGTTACCGCCAGAGACAAAGAAGTCGGAGACGAAGGAGGCGGAGCTGCCAGACACGGAGGAGAGCCATGAGGTGT---CGCCGCAGGTACAGGCTGAGACGTAGAAGCTGTTACCACATTGTATCT\n>Papio_cynocephalus\nATGGCCAGGTACAGATGCTGCCGCAGCCAGAGCCGAAGCAGATGCTATCGCCAGAGCCGGAGCAGATGTAACCGCCAGAGACAGAGCCAAAGCCGGAGAAGCTGCTATCGCCAGAGCCAAAGCCGGAGCAGATGTTACCGCCAGAGACAGAGAAGTCGTAGACGAAGGAGGCGACGCTGCCAGACACGGAGGAGAGCCATGAGGTGCTTCCGCCGCAGGTACAGGCTGAGGCGTAGGAGGCCCTATCACATCGTGTCT\n>Gorilla_gorilla\nATGGCCAGGTACAGATGCTGTCGCAGCCAGAGCCGCAGCAGATGTTACCGGCAGAGCCGGAGCAGGTGTTACCGGCAGAGACAAAGCCAGAGCCGGAGCAGATGCTACCGGCAGAGCCAAAGCCGGAGCAGGTGTTACCGGCAGAGACAAAGAAGTCGCAGACGTAGGCGGAGGAGCTGCCAGACACGGAGGAGAGCCATGAGGTGCTGCCGCCGCAGGTACAGACTGAGACGTAGAAGACCCTATCATATTGTATCT\n>Pan_troglodytes\nATGGCCAGGTACAGATGCTGTCGCAGCCAGAGCCGGAGCAGATGTTACCGGCAGAGACGGAGCAGGTGTTACCGGCAAAGGCAAAGCCAAAGTCGGAGCAGATGTTACCGGCAGAGCCAGAGACGGAGCAGGTGTTACCGGCAAAGACAAAGAAGTCGCAGACGAAGGCGACGGAGCTGCCAGACACGGAGGAGAGCCATGAGGTGCTGCCGCCGCAGGTACAGACTGAGACGTAAAAGATGTTACCATATTGTATCT'''
	tree.link_to_alignment(align)
	marks = ['2', '3', '4']

	Compute_selection_hyphy(tree, align, marks, "test", "test_tree", "test_align2", "test_out")
	exit()
	Compute_Selection(tree, marks, "test")
	exit()


if Example == True:
	Make_example()


###Process data
def participant_name(node_name_string):
	'For each leaf, get the participant ID'
	spcode = str(node_name_string).split(".O")[0].strip("\n-")
	return(spcode)


Participants = {}
for i in  tree.get_leaves():
	N = participant_name(i)
	if N not in Participants: Participants[N] = []
	Participants[N].append(str(i).strip("\n-") )
#Get a dic of dic[participant] = { Antibody_A : 1, Antibody_B : 0...}
Presence_profile = Get_presence_profile()
#Keep only leafs of participants with antibody info
Keep = [] 
for Participant in Participants:
	if Participant in Presence_profile: Keep.extend(Participants[Participant])
tree.prune(Keep)


#Get marks
n = 0
marks = []
for leaf in tree.get_leaves():
	ID = leaf.node_id
	Presence = Presence_profile[participant_name(leaf)][AB]
	if  Presence  == 1: marks.append(ID)

Compute_selection_hyphy(tree, open(Alignment_file).readlines(),  marks, "agilent_55961", Tree_temp, MSA_temp, Out)


