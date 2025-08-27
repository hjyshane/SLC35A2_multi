#$ -q all.q
#$ -j y
#$ -cwd 
#$ -S /bin/bash
#$ -terse
#$ -pe smp 24
#$ -N ptzall_integration_clustering


module rm python-2.7.12
module remove R-3.6.2
module load R-4.1.1
Rscript combo_all_qsub_full.R