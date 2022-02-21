from pathlib import Path
from Bio import SeqIO, AlignIO


configfile:"config.yaml"
#Reference:
#Path_file:
Biopython_path = "/scratch/umcg-sandreusanchez/Immuno_selection/Snake_pipeline/MAMBA_ENV/envs/snakemake/lib/python3.10/site-packages/Bio"
Peptide = []
for record in SeqIO.parse(config["Reference"], "fasta"):
	Peptide.append(str(record.id) )

Samples = {}
with open(config["Path_file"]) as F:
	for line in F:
		Name = Path(line).name.split("_kneaddata")[0]
		Forward = line.rstrip()
		Reverse = line.rstrip().replace("_1", "_2")
		Samples[Name] = [Forward, Reverse]
localrules: Pull_stats_together, Stats_with_recovered_reads

wildcard_constraints:
	peptide = "((agilent)|(twist))_\d+",
	sample = "[^_]+"
def Merge_samples(peptide, samples= Samples.keys() ):
	Files = []
	for s in samples:
		Files.append("5_Alignment_sample/{Sample}_{peptide}.fa".format(Sample=s, peptide=peptide))
	return(Files)
def All_alignments_per_peptide(peptide, prefix, samples=Samples.keys()):
	Files = []
	for s in samples:
		Files.append("{P}{Sample}_{peptide}.fa".format(P = prefix, Sample=s, peptide=peptide))
	return(Files)
def Get_reads(sample, dic_samples =Samples):
	return(dic_samples[sample])

Output = expand(['7_Evo_models/{peptide}_Output.txt','6_Alignment_Peptides/{peptide}.fa.iqtree', "Results/logo/{peptide}_Present.png", "Results/Phylogeny/{peptide}_tree.png"], peptide=Peptide)
Output.append("Results/PeptideAlignmentStats.tsv")
Output.append("Results/Selection_test.tsv")
Output.append("Results/Alignment_scores_comparison.tsv")
#print(Output)
#Output = [ "Results/Alignment_scores_comparison.tsv"] #Just do the first part
rule all:
	input:
		Output
	message: "Pipeline complete"
 
rule Pull_stats_together:
	input:
		expand("7_Evo_models/{peptide}_Output.txt", peptide=Peptide)
	output:
		"Results/Selection_test.tsv"
	shell:
		"python scripts/Read_output.py"

rule Compute_dnds:
	input:
		MSA = "6_Alignment_Peptides/{peptide}.fa",
		Tree = "6_Alignment_Peptides/{peptide}.fa.treefile"
	output:
		Output = "7_Evo_models/{peptide}_Output.txt",
		MSA_temp = temp("temp/MSA_filtered_{peptide}"),
		Tree_temp = temp("temp/Tree_FilteredAndMarked_{peptide}"),
	
	shell:
		"set +u; source activate /scratch/umcg-sandreusanchez/Immuno_selection/Snake_pipeline/ete3 ;\n"
		"PATH=$PATH:/scratch/umcg-sandreusanchez/Immuno_selection/Snake_pipeline/ete3/bin/ete3_apps/bin/ ;\n"
		"/scratch/umcg-sandreusanchez/Immuno_selection/Snake_pipeline/ete3/bin/python scripts/Compute_selection.py {input.MSA} {input.Tree}"			

rule Check_phylogeny:
	input: 
		Newick = '6_Alignment_Peptides/{peptide}.fa.treefile'
	output:
		Plot = "Results/Phylogeny/{peptide}_tree.png",
		Report = "Results/Phylogeny/{peptide}_report.tsv",
		temptree = temp("temp/{peptide}_tree.tsv")
	shell:
		"set +u; source activate /scratch/umcg-sandreusanchez/Immuno_selection/Conda_R/R ;\n"
		"Rscript scripts/PhyloWork.R {input.Newick} {output.temptree} {output.Plot} {output.Report}"

rule Logos:
	input:
		"6_Alignment_Peptides/{peptide}.fa"
	output:
		Present = "Results/logo/{peptide}_Present.png",
		NoPresent = "Results/logo/{peptide}_noPresent.png"		
	shell:
		"python scripts/Compare_logos.py {input}"

rule Phylogeny:
	input: 
		"6_Alignment_Peptides/{peptide}.fa"
	output: 
		Newick = '6_Alignment_Peptides/{peptide}.fa.treefile',
		Report = '6_Alignment_Peptides/{peptide}.fa.iqtree'
	shell:
		"set +u; source /software/software/Anaconda3/2020.11/bin/activate /scratch/umcg-sandreusanchez/Immuno_selection/2_AlignmentAndPhylogenetics/Software/Tools_2; set -u ;\n"
		"iqtree -s {input} -m MFP" #-st CODON11"
rule Alignment_peptide:
	input: 
		lambda wildcards: Merge_samples(wildcards.peptide)
	output:
		Merged = temp("6_Alignment_Peptides/tmp_{peptide}.fa"),
		tempaaMerged = temp("6_Alignment_Peptides/tmp_peptides_{peptide}.fa"),
		tempaaAligned = temp("6_Alignment_Peptides/tmp_peptides_aligned_{peptide}.fa"),
		MSA = "6_Alignment_Peptides/{peptide}.fa"
	shell:
		"cat {input} > {output.Merged} ;\n"
		"python scripts/Second_alignment.py {output.Merged} {output.tempaaMerged} {output.tempaaAligned} {output.MSA} "

rule Collapse_Alignments:
	input:
		lambda wildcards: All_alignments_per_peptide(wildcards.peptide, "5_Alignment_sample/temp_")
	output:
		expand("5_Alignment_sample/{Sample}_{{peptide}}.fa", Sample=Samples.keys())
		#Collapsed[wildcards.peptide]
		
	run:
		for i in input:
			O = "5_Alignment_sample/{S}_{P}.fa".format(S= Path(i).stem.split("_")[1], P= wildcards.peptide)
			#print("python scripts/Parallelize_Remove_redundancy_MS.py "+i + " " + O + " " + Path(i).stem.split("_")[1] )
			shell("python scripts/Parallelize_Remove_redundancy_MS.py "+i + " " + O+ " " + Path(i).stem.split("_")[1] )

rule Alignments_per_participant:
	input:
		lambda wildcards: All_alignments_per_peptide(wildcards.peptide, "4_Prepared_reads/"),
	output:
		temp( expand("5_Alignment_sample/temp_{sample}_{{peptide}}.fa", sample = Samples.keys() ))
	run:
		for i in input:
			print("Aligning "+ i)
			O = "5_Alignment_sample/temp_{S}_{P}.fa".format(S= Path(i).stem.split("_")[0], P= wildcards.peptide)

			Records = SeqIO.parse(i , "fasta")
			S = sum(1 for _ in Records)
			if S < 2: 
				shell("cp " + i + " " + O)
			else:
				print("mafft --reorder --thread 1  --auto " + i +" > " + O)
				shell("ml MAFFT/7.453-GCC-9.3.0-with-extensions ;\n mafft --reorder --thread 1  --auto " + i +" > " + O)
rule NumberReadsRecovered:
	input:
		expand("4_Prepared_reads/{sample}_{peptide}.fa", sample = Samples.keys(), peptide = Peptide)
	output:
		"Results/PeptideAlignmentStats.tsv"
	shell:
		"python scripts/StatsReadsPerPetide.py"
rule Stats_with_recovered_reads:
	input:
		expand("3_Diamond_output/{sample}.tsv", sample = Samples.keys())	
	output:
		"Results/Alignment_scores_comparison.tsv"
	shell:
		"ml Anaconda3/2020.11 ;\n"
		"set +u; source activate /scratch/umcg-sandreusanchez/Immuno_selection/Conda_R/R ; set -u ;\n"	
		"Rscript scripts/Compare_identity.R"

rule Extract_and_clean_read:
	input:
		Fastq = "1_Read_merging/{sample}.fq",
		Diamond = "3_Diamond_output/{sample}.tsv"
		
	output:
		Read_save = temp("1_Read_merging/Keep_{sample}.txt"),
		Fastq1 = "1_Read_merging/subset_{sample}.fq",
		Rest = expand("4_Prepared_reads/{{sample}}_{peptide}.fa", peptide=Peptide)
	params:
		Peptide_list = config["Reference"],
		Prefix_out = "4_Prepared_reads/{sample}"
	shell:
		"set +u; source /software/software/Anaconda3/2020.11/bin/activate /data/umcg-sandreusanchez/MAP; set -u ;\n"
		"cut -f1 {input.Diamond} > {output.Read_save} ;\n"
		"seqtk subseq {input.Fastq} {output.Read_save} > {output.Fastq1} ;\n"
		"conda deactivate ;\n"
		"PYTHONPATH={Biopython_path} ;\n"
		"/scratch/umcg-sandreusanchez/Immuno_selection/Snake_pipeline/Mamba_snakemake/bin/python scripts/Trim_and_prepare.py {input.Diamond} {output.Fastq1} {params.Peptide_list} {params.Prefix_out}"
rule run_Diamond:
	input:
		FASTA = "2_FASTA/{sample}.fa",
		DB = "0_Database/" + Path(config["Reference"]).stem + ".dmnd"
	output: 
		"3_Diamond_output/{sample}.tsv"
	params:
		DB_name =  "0_Database/" + Path(config["Reference"]).stem
	shell:
		"set +u; source /software/software/Anaconda3/2020.11/bin/activate /data/umcg-tifn/rgacesa/conda_dag3_v3; set -u ;\n"
        	"diamond blastx -q {input.FASTA} -d {params.DB_name} -o {output} --very-sensitive"

rule make_fasta_and_merge:
	input:
		lambda wildcards: Get_reads(wildcards.sample)
	output:
		unmerged = temp("1_Read_merging/unmerged-{sample}.fq"),
		merged = temp("1_Read_merging/merged-{sample}.fq"), 
		total = "1_Read_merging/{sample}.fq",
		Final = "2_FASTA/{sample}.fa"
	shell:
		"set +u; source /software/software/Anaconda3/2020.11/bin/activate /data/umcg-tifn/rgacesa/conda_BB; set -u ;\n"
		"bbmerge.sh in={input[0]} in2={input[1]} out={output.merged} outu={output.unmerged} ;\n"
		"cat {output.merged} {output.unmerged} > {output.total} ;\n"
		"reformat.sh in={output.total} out={output.Final} > /dev/null"


rule Create_database:
	input:
		config["Reference"]
	output:
		"0_Database/" + Path(config["Reference"]).stem + ".dmnd"
	params:
		DB_name =  "0_Database/" + Path(config["Reference"]).stem
	shell:
		"set +u; source /software/software/Anaconda3/2020.11/bin/activate /data/umcg-tifn/rgacesa/conda_dag3_v3; set -u  ;\n"
		"diamond makedb --in {input} -d {params.DB_name}"
