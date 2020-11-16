# gghybrid_pipeline
How to run and perform various analyses using the R package [gghybrid](https://github.com/ribailey/gghybrid) on RADseq data

Latest update: 11/16/2020

This document provides notes and instructions on how to run RADseq (Restriction Site Associated DNA Sequencing) data using the R package [gghybrid](https://github.com/ribailey/gghybrid). gghybrid is an R package for the analysis of hybrid zones. You can use gghybrid to obtain a hybrid index for your dataset as well as cline data such as cline centers, cline widths, and graphs of these clines. This package is for genomic clines NOT geographic clines. For geographic clines, I suggest using the R package [HZAR](https://cran.r-project.org/web/packages/hzar/index.html). 

This document serves as an additonal resource to the manual provided by ribaily on the gghybrid home page as well as a resource for additional analyses one can perform after running data through gghybrid, such as an enrichment analysis of significant clines and some graph tweaking. 

## Overall Pipeline Summary
raw RAD reads (fastq file) -> see [RADseq_pipeline](https://github.com/kiralong/RADseq_pipeline) -> structure file -> `gghybrid` -> other analyses (enrichment analysis)

## Pipeline Steps

### Step 1: Process raw RAD reads
See [RADseq_pipeline](https://github.com/kiralong/RADseq_pipeline) for instructions on how to process raw RAD reads after you get your giant fastq.gz file from the sequencing facility. At the step for running `populations` in `stacks`, you will need to make sure you add the flag `--structure` to have the `populations` module output a structure file, the preferred input of gghybrid.

For an example script of running `populations` to get the structure file see [run_populations.sh](run_populations.sh)

### Step 2: Format starting structure file for gghbyrid
Now that you have your structure file of your RAD data, you need to format this structure file for gghbyrid. You can use the script [format_structurefile_for_gghybrid.sh](format_structurefile_for_gghybrid.sh) to accomplish this. 

### Step 3: Make a list of the SNPs and chromosome coordinates
Use the following bash command to make a list of SNP and chromosome coordinates to aid in labeling the gghybrid output graphs. 

```
sumstats_file_path=/projects/aces/kira2/stacks_ref_map/ref_map_Mar30_2020/cline_analysis/gghybrid/2020_October_15_enrichment/populations_p6_r0.50_mac3_singleSNP_HWE_HZAR/populations.sumstats.tsv

cat $sumstats_file_path | grep -v "^#" | cut -f 1-4 | sort -u | sed -E 's/^([0-9]+)\t([A-Z]{2}_[0-9]+\.[0-9]+)\t([0-9]+)\t([0-9]+)/\1_\4\t\2\t\3/' > snp_chromosome_coordinates_r50.tsv

```

### Step 4: Run gghybrid

If continuing on to do enrichment analysis, use the R script [gghybrid.enrichment.R](gghybrid.enrichment.R).

### Step 5: Enrichment Analysis (Optional)


