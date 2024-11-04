#!/usr/bin/env bash

#SBATCH --job-name="mat"
#SBATCH --account=rrg-kleinman
#SBATCH --output=.logs/mat-%j.out
#SBATCH --time=01:00:00
#SBATCH --mem-per-cpu=30G
#SBATCH --ntasks=1

# Prep R environment
module load StdEnv/2020
module load r/4.1.2
export R_LIBS_USER="../../../../renv/library/R-4.1/x86_64-pc-linux-gnu"

# Generate raw counts CSV with R
R --vanilla -e "
  library(Seurat)
  load('seurat.Rda')
  expr <- t(as.matrix(GetAssayData(seurat, slot = 'counts')))
  if (!file.exists('expression.csv')) write.csv(expr, file.path(out, 'expression.csv'))
  "
  
echo "done."