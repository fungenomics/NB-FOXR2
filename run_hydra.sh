#!/usr/bin/bash
#SBATCH --job-name="NB-FOXR2"
#SBATCH --time=01:00:00
#SBATCH --account=kleinman-lab
#SBATCH --cpus-per-task=1
#SBATCH --mem=96G

mkdir -p logs/

module load StdEnv/2020
module load r/4.1.2

# get current resources requested
RESOURCES=`grep -E '^#SBATCH' run_hydra.sh | grep -v 'job-name' | grep -v 'account' | grep -v 'output'`

# usage:
# sbatch --export=ALL,rmd='01-my_analysis.Rmd' run_slurm.sh
R --no-save -e "rmarkdown::render('${rmd}', 'html_document', params = list(resources = '${RESOURCES}'))"
