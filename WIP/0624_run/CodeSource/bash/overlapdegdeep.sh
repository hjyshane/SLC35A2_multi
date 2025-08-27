#$ -q all.q
#$ -j y                      # Merge stdout and stderr
#$ -cwd                      # Use the current working directory
#$ -S /bin/bash              # Use Bash shell
#$ -terse                    # Output only job ID
#$ -pe smp 32                # Request  cores
#$ -N overlapdegdeep                    # Name the job
#$ -l h_vmem=256G            # Request  RAM
#$ -m beas                   # Email at the beginning and end of the job
#$ -M june.yoon@nationwidechildrens.org
#$ -o /igm/home/hxy008/PTZ_ATAC_scRNA_072024/WIP/Logs_bash/overlapdegdeep.log           # Save combined stdout/stderr to this log file


# Load the R module (adjust for your environment)
module purge
module load R-4.3.3 ## change this to curret R version
Rscript /igm/home/hxy008/PTZ_ATAC_scRNA_072024/WIP/0624_run/CodeSource/overlapdegdeep.R ## change the R script path of yours
