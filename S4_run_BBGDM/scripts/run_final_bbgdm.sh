#!/bin/bash -l
#SBATCH --job-name=bbgdm
#SBATCH --account=Project_2004932
#SBATCH --output=output/out_%a.txt
#SBATCH --error=error/err_%a.txt
#SBATCH --partition=small
#SBATCH --time=00:15:00
#SBATCH --ntasks-per-node=1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=5000
#SBATCH --mail-type=FAIL,ARRAY_TASKS
#SBATCH --array=1-163%99
# Load r-env-singularity
module load r-env-singularity

# Clean up .Renviron file in home directory
if test -f ~/.Renviron; then
sed -i '/TMPDIR/d' ~/.Renviron
fi

# Specify a temp folder path
echo "TMPDIR=/scratch/project_2004932" >> ~/.Renviron

# Run the R script
srun singularity_wrapper exec Rscript --no-save run_bbgdm.R

