#!/bin/bash
#SBATCH --job-name=EB_soupx_analysis
#SBATCH --output=soupx_analysis.out
#SBATCH --error=soupx_analysis.err
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --time=72:00:00
#SBATCH --mem=240G

# Load necessary modules
module load R

# Define sample names and output directory
SAMPLES=("S2196_eb" "S2197_eb" "S2198_eb" "S2199_eb")
OUTPUT_DIR="/data/ebaird/scRNAseq/soupx"

# Loop over samples and submit a job for each
for SAMPLE in "${SAMPLES[@]}"; do
  sbatch --wrap="Rscript /data/ebaird/scRNAseq/soupx/soupx.r $SAMPLE $OUTPUT_DIR"
done