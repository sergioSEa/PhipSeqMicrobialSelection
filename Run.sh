ml Anaconda3/2020.11
#source activate /scratch/umcg-sandreusanchez/Immuno_selection/SNAKE
#source activate /scratch/umcg-sandreusanchez/Immuno_selection/Snake_pipeline/Mamba_snakemake
#/scratch/umcg-sandreusanchez/Immuno_selection/SNAKE/bin/snakemake
#PYTHONPATH=/scratch/umcg-sandreusanchez/Immuno_selection/SNAKE/lib/python3.6/site-packages/
#PYTHONPATH=/scratch/umcg-sandreusanchez/Immuno_selection/Snake_pipeline/Mamba_snakemake/lib/python3.6/site-packages/

source activate /scratch/umcg-sandreusanchez/Immuno_selection/Snake_pipeline/MAMBA_ENV/envs/snakemake

snakemake --latency-wait 60 --rerun-incomplete  --jobs 99 --keep-going --cluster 'sbatch -t {cluster.time} --mem={cluster.mem} --cpus-per-task={cluster.c} --error={cluster.error} --job-name={cluster.name} --output={cluster.output}'  --cluster-config cluster.json --snakefile Phylo_Adaptation2.smk --cluster-status 'python slurm-status.py'

#snakemake -s Phylo_Adaptation.smk -np


