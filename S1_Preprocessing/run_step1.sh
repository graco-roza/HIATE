#!/bin/bash -l
#SBATCH --job-name=bbgdm
#SBATCH --account=Project_2004932
#SBATCH --error=err/err_%a.txt
#SBATCH --output=out/out_%a.txt
#SBATCH --partition=small
#SBATCH --time=01:00:00
#SBATCH --ntasks-per-node=1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=5000
#SBATCH --mail-type=FAIL,ARRAY_TASKS
#SBATCH --array=1-52,54-92,94-97,99-153,155-159,161-171%9
# Load r-env-singularity
module load r-env-singularity

# Clean up .Renviron file in home directory
if test -f ~/.Renviron; then
sed -i '/TMPDIR/d' ~/.Renviron
fi

# Specify a temp folder path
echo "TMPDIR=/scratch/project_2004932" >> ~/.Renviron

# Run the R script
srun singularity_wrapper exec Rscript --no-save pre_processing.R
