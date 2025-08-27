#$ -q all.q
#$ -j y
#$ -cwd
#$ -S /bin/bash
#$ -terse
#$ -pe smp 24
#$ -N annotation

module load R-4.0.2 ## update this with current R version
Rscript /igm/home/tab013/snRNA_JLE/analysis/ref_and_query_allen.R
