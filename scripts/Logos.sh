#The idea is to split the MSA in two files, one for participants with the antibody and another for participants without the antibody.
#Run this script for each of the two groups separately and compare
ml Anaconda3/2020.11
source activate /scratch/umcg-sandreusanchez/Immuno_selection/2_AlignmentAndPhylogenetics/Software/Tools_2
PYTHONPATH=/scratch/umcg-sandreusanchez/Immuno_selection/2_AlignmentAndPhylogenetics/Software/Tools_2/lib/python3.9/site-packages/

MSA=$1
Output=$2

weblogo --format PNG < $MSA > $Output

