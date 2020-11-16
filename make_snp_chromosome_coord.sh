#!/bin/bash

sumstats=/projects/aces/kira/stacks_ref_map/ref_map_Mar30_2020/cline_analysis/gghybrid/2020_October_15_enrichment/populations_p6_r0.50_mac3_singleSNP_HWE_HZAR/populations.sumstats.tsv

cat $sumstats | grep -v "^#" | cut -f 1-4 | sort -u | sed -E 's/^([0-9]+)\t([A-Z]{2}_[0-9]+\.[0-9]+)\t([0-9]+)\t([0-9]+)/\1_\4\t\2\t\3/' > snp_chromosome_coordinates_r50.tsv
