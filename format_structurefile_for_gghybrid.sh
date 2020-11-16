#!/bin/bash

structure_file=/projects/aces/kira/stacks_ref_map/ref_map_Mar30_2020/cline_analysis/gghybrid/2020_October_15_enrichment/populations_p6_r0.50_mac3_singleSNP_HWE_HZAR/populations.structure
output_file=populations.p6.r50.sSNP.enrichment.structure.csv

cat $structure_file | grep -v "^#" | tr '\t' ',' | sed -E 's/,0/,NA/g'| sed -E 's/,,/INDLABEL,POPID,/' |
	sed 's/,pr/,09PR/' | sed 's/,ru/,08RU/' | sed 's/,cg/,10CG/' | sed 's/,ss/,02SS/' | sed 's/,qp/,06QP/' | sed 's/,ro/,05RO/' > $output_file
