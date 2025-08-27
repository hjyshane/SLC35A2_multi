# copy files from igm to franklin
# from IGM to franklin
rsync -av rpl-igmmaster01:/igm/home/hxy008/PTZ_ATAC_scRNA_072024/WIP/File  /home/gdbedrosianlab/hxy008/Multiome_122424/WIP
rsync -av rpl-igmmaster01:/igm/home/hxy008/PTZ_ATAC_scRNA_072024/File/  /home/gdbedrosianlab/hxy008/Multiome/Data

# from franklin to IGM
rsync -av r1pl-hpcf-log01:/home/gdbedrosianlab/hxy008/Multiome /igm/home/hxy008

#R path franklin
/export/apps/opt/R/4.4.1-foss-2020a/bin/R-4.3.3/bin/R

/apps/easybuild/software/R/4.0.0-foss-2020a/bin/R
/apps/opt/R/4.4.1-foss-2020a/bin/R
/apps/x86_64/easybuild/software/R/4.0.0-foss-2020a/bin/R
/apps/x86_64/opt/R/4.4.1-foss-2020a/bin/R
/export/apps/easybuild/software/R/4.0.0-foss-2020a/bin/R
/export/apps/x86_64/opt/R/4.4.1-foss-2020a/bin/R
/export/apps/x86_64/easybuild/software/R/4.0.0-foss-2020a/bin/R
/export/apps/x86_64/opt/R/4.4.1-foss-2020a/bin/R

# VS code setting for Franklin
    "r.libPaths": [
        "/export/apps/opt/R/4.4.1-foss-2020a/lib64/R/library",
    ],
    "r.rterm.linux": "/export/apps/opt/R/4.4.1-foss-2020a/bin/R",
"terminal.integrated.env.linux": {
    "PATH": "/export/apps/opt/R/4.4.1-foss-2020a/bin:$PATH"
}

# add R path to run"
export PATH=/export/apps/opt/R/4.4.1-foss-2020a/bin:$PATH
source ~/.bashrc
echo $PATH
ls -l //export/apps/opt/R/4.4.1-foss-2020a/bin/R
    -rwxr-xr-x 1 softw softw 9346 Jul  5 15:39 //export/apps/opt/R/4.4.1-foss-2020a/bin/R

# when commend not found
export PATH=/export/apps/opt/R/4.4.1-foss-2020a/bin:/bin:/usr/bin:/usr/local/bin:$PATH
source ~/.bashrc
echo $PATH
nano ~/.bashrc
    in side the file, add below at the end of the line.
    export PATH=/export/apps/opt/R/4.4.1-foss-2020a/bin:/bin:/usr/bin:/usr/local/bin:$PATH
    ^+O to save, ^+X to exit
source ~/.bashrc

# install R
/export/apps/opt/R/4.4.1-foss-2020a/bin/R --silent --slave --no-save --no-restore -e install.packages('languageserver', repos='https://cran.r-project.org/') 


#R path IGM
/igm/apps/R/R-4.3.3/bin/R

# VS code setting for IGM
    "r.libPaths": [
        "/igm/apps/R/R-4.3.3/library",
        "/igm/home/hxy008/R/x86_64-pc-linux-gnu-library/4.3",
        "/igm/home/hxy008/R/x86_64-pc-linux-gnu-library/library"
    ],
    "r.rterm.linux": "/igm/apps/R/R-4.3.3/bin/R",

# add R path to run"
export PATH=//igm/apps/R/R-4.3.3/bin:$PATH
source ~/.bashrc
echo $PATH
ls -l /igm/apps/R/R-4.3.3/bin/R

# when commend not found
export PATH=/igm/apps/R/R-4.3.3/bin/:/bin:/usr/bin:/usr/local/bin:$PATH
source ~/.bashrc
echo $PATH
nano ~/.bashrc
    in side the file, add below at the end of the line.
    export PATH=/igm/apps/R/R-4.3.3/bin:/bin:/usr/bin:/usr/local/bin:$PATH
    ^+O to save, ^+X to exit
source ~/.bashrc

# remove

rm -r for directory
rm -rf for difectory but no checking individual files

# SGE
#$ -N <job-name>	Set the name of the job.	#$ -N MyJob
#$ -cwd	Use the current working directory.	#$ -cwd
#$ -S <shell>	Specify the shell to use.	#$ -S /bin/bash
#$ -j y	Merge stdout and stderr into a single log file.	#$ -j y
#$ -o <path>	Specify the stdout log file.	#$ -o /path/to/output.log
#$ -e <path>	Specify the stderr log file.	#$ -e /path/to/error.log
#$ -q <queue-name>	Specify the queue to submit the job.	#$ -q all.q
#$ -pe smp <slots>	Request a specific number of slots (shared memory).	#$ -pe smp 8
#$ -l h_vmem=<memory>	Request memory per core.	#$ -l h_vmem=8G
#$ -l h_rt=<time>	Request runtime limit.	#$ -l h_rt=12:00:00
#$ -l disk=<size>	Request local disk space.	#$ -l disk=10G
#$ -l gpu=<num>	Request GPU resources.	#$ -l gpu=1
#$ -v <var=value>	Pass an environment variable to the job.	#$ -v DEBUG=1
#$ -m <events>	Enable email notifications (b, e, a, s).	#$ -m beas
#$ -M <email>	Specify the recipient email address.	#$ -M user@example.com
#$ -t <range>	Submit a job array with task range.	#$ -t 1-100
#$ -t <start>-<end>:<step>	Specify task range with a step size.	#$ -t 1-100:10
#$ -hold_jid <job-id>	Hold the job until another job completes.	#$ -hold_jid 1234
#$ -p <priority>	Set job priority (-1024 to 1023).	#$ -p -10
#$ -binding linear:<num>	Bind the job to consecutive CPU cores.	#$ -binding linear:4
#$ -R y	Reserve resources without starting immediately.	#$ -R y
#$ -terse	Output only the job ID upon submission.	#$ -terse
#$ -q <queue-name>	Submit the job to a specific queue.	#$ -q test.q
#$ -l arch=<architecture>	Request a specific CPU architecture.	#$ -l arch=x86_64
#$ -o <file>	Redirect the output log to a file.	#$ -o /path/to/output.log
#$ -e <file>	Redirect the error log to a file.	#$ -e /path/to/error.log
#$ -pe mpi <slots>	Request slots for distributed memory (MPI).	#$ -pe mpi 16
#$ -l hostname=<host>	Specify a particular host to run the job.	#$ -l hostname=node01
qsub	Submit a job.	qsub myscript.sh
qstat	Check the status of all jobs.	qstat
qstat -j <job-id>	Get detailed information about a job.	qstat -j 1234
qstat -u <username>	Show jobs submitted by a specific user.	qstat -u user
qdel <job-id>	Cancel a running or queued job.	qdel 1234
qacct -j <job-id>	Get accounting information for a completed job.	qacct -j 1234
qacct -o <username>	Summarize all jobs for a specific user.	qacct -o user
qhost	Display available nodes and their statuses.	qhost
qconf -sc	Display available resource configurations.	qconf -sc
qalter -N <name> <job-id>	Change the name of a running or queued job.	qalter -N NewJobName 1234
qalter -p <priority>	Change the priority of a job.	qalter -p 0 1234
qrls <job-id>	Release a held job.	qrls 1234
qsuspend <job-id>	Temporarily suspend a running job.	qsuspend 1234
qresume <job-id>	Resume a suspended job.	qresume 1234
qsub -verify <script>	Test script without executing it.	qsub -verify myscript.sh

# Slurm
Directive	Purpose	Example
#SBATCH --job-name=<name>	Name the job.	#SBATCH --job-name=MyJob
#SBATCH --output=<file>	Specify stdout log file.	#SBATCH --output=/path/to/output.log
#SBATCH --error=<file>	Specify stderr log file.	#SBATCH --error=/path/to/error.log
#SBATCH --partition=<partition>	Specify the partition/queue to submit the job.	#SBATCH --partition=compute
#SBATCH --nodes=<n>	Request a specific number of nodes.	#SBATCH --nodes=2
#SBATCH --ntasks=<n>	Request a specific number of tasks/processes.	#SBATCH --ntasks=4
#SBATCH --cpus-per-task=<n>	Request CPUs per task.	#SBATCH --cpus-per-task=8
#SBATCH --time=<time>	Request runtime (format: days-hours:minutes:seconds).	#SBATCH --time=1-12:00:00
#SBATCH --mem=<memory>	Request memory per node.	#SBATCH --mem=64G
#SBATCH --mem-per-cpu=<memory>	Request memory per CPU.	#SBATCH --mem-per-cpu=8G
#SBATCH --gres=gpu:<n>	Request GPU resources.	#SBATCH --gres=gpu:2
#SBATCH --mail-type=<type>	Specify email notifications (BEGIN, END, FAIL, ALL).	#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=<email>	Specify recipient email address.	#SBATCH --mail-user=user@example.com
#SBATCH --dependency=<type>	Set job dependencies (afterok, afternotok, after).	#SBATCH --dependency=afterok:<job-id>
#SBATCH --array=<range>	Submit a job array.	#SBATCH --array=1-10
#SBATCH --account=<account>	Specify an account for job submission.	#SBATCH --account=myaccount
#SBATCH --exclusive	Run job on nodes without sharing resources with other jobs.	#SBATCH --exclusive
#SBATCH --export=<var=value>	Export environment variables to the job.	#SBATCH --export=ALL,DEBUG=1
#SBATCH --requeue	Allow the job to be requeued upon failure.	#SBATCH --requeue
#SBATCH --no-requeue	Prevent the job from being requeued.	#SBATCH --no-requeue
#SBATCH --begin=<time>	Delay job start until the specified time.	#SBATCH --begin=2023-12-25T12:00:00
#SBATCH --chdir=<dir>	Set the working directory for the job.	#SBATCH --chdir=/path/to/dir
#SBATCH --constraint=<feature>	Request specific hardware constraints.	#SBATCH --constraint=haswell
#SBATCH --licenses=<licenses>	Specify software licenses required for the job.	#SBATCH --licenses=matlab:1
#SBATCH --kill-on-invalid-dep	Cancel job if dependency is invalid.	#SBATCH --kill-on-invalid-dep

Command	Purpose	Example
sbatch <script>	Submit a batch job script.	sbatch myscript.sh
srun <command>	Run a job interactively.	srun --ntasks=4 ./my_program
salloc	Allocate resources for an interactive session.	salloc --ntasks=4 --time=2:00:00
scancel <job-id>	Cancel a running or queued job.	scancel 1234
squeue	View the status of all jobs.	squeue
squeue -u <username>	View jobs for a specific user.	squeue -u user
squeue -j <job-id>	View details of a specific job.	squeue -j 1234
sinfo	Display information about available partitions and nodes.	sinfo
scontrol show job <id>	Display detailed information about a job.	scontrol show job 1234
scontrol hold <job-id>	Hold a job (prevent it from starting).	scontrol hold 1234
scontrol release <job-id>	Release a held job.	scontrol release 1234
sacct	Display job accounting information for completed jobs.	sacct -j 1234
sstat	Display resource usage statistics for running jobs.	sstat -j 1234
sprio	Display job priority information.	sprio
sbatch --test-only <script>	Test a job script without running it.	sbatch --test-only myscript.sh
sattach <job-step-id>	Attach to a running job (useful for interactive jobs).	sattach 1234.0
scontrol show node <node>	Display details about a specific node.	scontrol show node node01
scontrol update	Modify a running job’s parameters (e.g., runtime, resources).	scontrol update jobid=1234 time=01:00:00