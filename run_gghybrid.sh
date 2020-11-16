#!/bin/bash
#SBATCH -p aces
#SBATCH -t 168:00:00
#SBATCH -N 1
#SBATCH -n 1

module load R/.4.0.1

cd /projects/aces/kira/stacks_ref_map/ref_map_Mar30_2020/cline_analysis/gghybrid/2020_October_15_enrichment

Rscript ./gghybrid.enrichment.R
