#!/usr/bin/env bash

#SBATCH -J findMotifs
#SBATCH -A rrg-kleinman
#SBATCH -N 1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH -t 3:00:00
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err
#SBATCH --mem=12G

# Load modules needed to run Homer to annotate the peaks 
module load StdEnv/2020
export MUGQIC_INSTALL_HOME=/cvmfs/soft.mugqic/CentOS6
module use $MUGQIC_INSTALL_HOME/modulefiles
module load mugqic/perl/5.22.1
module load mugqic/homer/4.11
module load mugqic/weblogo/3.3

# Parameters to modify
PEAK_FILE="../../../output/09/exact_intersect_ChIP_ATAC_qVal6_FC6_width10k.bed"
OUT_DIR="../../../output/09/motifs_size100_100"
TMP_DIR="../../../output/09/.tmp"

mkdir $TMP_DIR
mkdir $OUT_DIR

# Modify the file peak file to make it homer friendly 
# 1. add "chr" at the beginning of all line
# 2. use a uniq name for all the peaks
sed 's/^/chr/' $PEAK_FILE | awk 'FS=OFS="\t" {$4 = "peak_" NR } {print}' > "$TMP_DIR/homer_peaks.bed"

# Run Homer's findMotifs
findMotifsGenome.pl \
  $TMP_DIR/homer_peaks.bed \
  mm10 \
  $OUT_DIR \
  -preparsedDir $OUT_DIR/preparsed \
  -size -100,100 \
  -p 4

rm -rf $TMP_DIR