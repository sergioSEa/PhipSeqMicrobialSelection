library(tidyverse)
library(ape)
library(ggtree)
library(vegan)
#library(tidytree)
#library(treeio)

sessionInfo()

args = commandArgs(trailingOnly=TRUE)
File = args[1] #"6_Alignment_Peptides/agilent_55961.fa.treefile"
Temp_tree = args[2]
Tree_output = args[3]
Report = args[4]


read.tree(File) -> Tree
Peptide = str_split(str_split(File,"/")[[1]][ length(str_split(File,"/")[[1]]) ], ".fa")[[1]][1]
#Annotation...
sapply(Tree$tip.label, function(x){ str_split(x, ".Orth")[[1]][1]  }) -> Participants
tibble(ID = Participants, label = Tree$tip.label ) %>% mutate(ID2 = as.numeric(as.factor(ID)) ) -> Annotation

read_tsv("/scratch/umcg-sandreusanchez/Immuno_selection/Data/Matrix_selection.tsv") %>% select(ID, Peptide) -> Info
Linkage_file = read_tsv("/scratch/umcg-sandreusanchez/Immuno_selection/Data/Translate_MGS_names.txt") %>% mutate(ID = LLD)
left_join(Info, Linkage_file) %>% drop_na() %>% mutate(ID = MGS) -> Info
Info %>% select(Peptide) %>% as_vector() -> Peptide_profile
Info %>% select(-Peptide) %>% mutate(Presence_absence = as.factor(Peptide_profile)) %>% select(-MGS) -> Info

left_join(Annotation, Info) %>% drop_na() -> Annotation

# Plot tree, Annotation in branches: 1. Participant, 2. Presence/Absence of antibody
drop.tip( Tree, Tree$tip.label[! Tree$tip.label %in% Annotation$label]) -> Tree
Tree_ultra=chronos(Tree, lambda=0) 
multi2di(Tree_ultra, random=TRUE) -> Tree_ultra

#force.ultrametric(Tree_ultra, method="extend") -> Tree_ultra
print("Saving ultrametric tree")
write.tree(Tree_ultra, Temp_tree)
print("Reading ultrametric")
read.tree(Temp_tree) -> Tree_ultra

print(Tree_ultra)
full_join(Tree_ultra, Annotation, by= "label")  -> Annotated_tree

ggtree(Annotated_tree,layout = "circular") + geom_tiplab(aes(label= ID2, angle=angle), size= 2) +
  geom_tippoint(aes(col=Presence_absence)) -> Tree_plot
ggsave(Tree_output, Tree_plot)

cophenetic(Tree) -> Dist_matrix
adonis2(Dist_matrix ~ Presence_absence, Annotation) -> Results # We need to control for ID
as.data.frame(Results) %>% rownames_to_column("Ar") %>% filter(Ar == "Presence_absence") %>% select(-Ar) -> Results
write_tsv(Results, Report)


