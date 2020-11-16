#!/bin/bash
#SBATCH -p aces
#SBATCH -J populations_enrichment_r.50
#SBATCH -t 48:00:00
#SBATCH -N 1
#SBATCH -n 12

module load gcc/7.2.0

stacks=/projects/aces/kira/programs/stacks-2.53/bin
popmap_path=/projects/aces/kira/sample_info/new_popmap.tsv
stacks_ref_map_output=/projects/aces/kira/stacks_ref_map/ref_map_Mar30_2020
populations_output_path=/projects/aces/kira/stacks_ref_map/ref_map_Mar30_2020/cline_analysis/gghybrid/2020_October_15_enrichment

p=6
mac=3
r=0.50
populations_output=$populations_output_path/populations_p${p}_r${r}_mac${mac}_singleSNP_HWE_HZAR
mkdir -p $populations_output

# Stacks command
cmd=(
	$stacks/populations
	--in-path $stacks_ref_map_output
	--out-path $populations_output
	--popmap $popmap_path
	--threads 12
	--min-mac $mac
	--min-population $p
	--min-samples-per-pop $r
	--write-single-snp
	--structure
	--hwe
	--hzar
)

# Run command
"${cmd[@]}"
