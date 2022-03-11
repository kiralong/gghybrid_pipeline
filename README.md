# gghybrid_pipeline
How to run and perform various analyses using the R package [gghybrid](https://github.com/ribailey/gghybrid) on RADseq data

Latest update: 05/20/2021

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
sumstats_file_path=/projects/aces/kira2/stacks_ref_map/ref_map_Mar30_2020/cline_analysis/gghybrid/
2020_October_15_enrichment/populations_p6_r0.50_mac3_singleSNP_HWE_HZAR/populations.sumstats.tsv

cat $sumstats_file_path | grep -v "^#" | cut -f 1-4 | sort -u | \
sed -E 's/^([0-9]+)\t([A-Z]{2}_[0-9]+\.[0-9]+)\t([0-9]+)\t([0-9]+)/\1_\4\t\2\t\3/' 
> snp_chromosome_coordinates_r50.tsv

```

### Step 4: Run gghybrid

If continuing on to do enrichment analysis, use the R script [gghybrid.enrichment.R](gghybrid.enrichment.R). You will need to provide the formated structure file you got as the output from running [format_structurefile_for_gghybrid.sh](format_structurefile_for_gghybrid.sh) and you will need to provide the list of SNP and chromosome/scaffold coordinates you got from step 3. If you are running a large number of SNPs (say around 30,000) gghybrid will run for around 3 days. At the end you should have 5 documents: gghybrid_graphs.pdf (which contains your hybrid index and a page for each individual genomic cline for every snp you ran. So yes if you ran 30K snps, you will have a 30K page pdf.), cline_table.txt, HI_table.txt, genotype_clines.txt, and locus_clines.txt. 

To deal with this giant pdf, my current solution is to use `greppdf` in a command line to search through the pdf for snps of interest. From your output from step 3, your snp_chromosome_coordinates.tsv file, all individual genomic cline graphs will have the following pattern:
"SNP ID: (snp_id) at (scaffold.ID) (basepair in genome)"
For example "SNP ID: 68476_67 at NW_021940545.1 6596978"

So you can use `greppdf 'pattern' file.pdf` to find a specific SNP's genomic cline graph. E.g. `pdfgrep 'SNP ID: 6177_27 at NW_021940545.1 BP' gghybrid_graphs_r25.pdf`

Another strategy you can use is grep -n to get the line number of the loci from the cline table and then go to that page in the pdf with whatever software you use to view pdfs on your computer. 

### Step 5: Enrichment Analysis (Optional)

For the enrichment analysis, I have a list of genome windows that were found to be significant from a number of other genomic tests that I got from a collaborrator from their whole genome data. We then asked whether my RAD data had enriched significant genomice clines within these same outlier windows from the whole genome data. So for this next part, you will need a list of genome windows of interest.

Once you have you list of genomic windows you want to check the clines in, you can run [filter_cline_table.py](filter_cline_table.py). This script takes you genome windows file and your cline_table.txt output from gghybrid. 

To run python script, you will need python 3 and can use:
./filter_cline_table.py -s populations.sumstats.tsv -w genome_windows.tsv -c cline_table.txt -o cline_table.windows.tsv
-s the sumstats file from stacks after running the populations module
-w the list of windows in the genome that were significant from any other genome-wide test
-c cline table that is output from gghybrid that has the cline center p-values
-o the output file name for the file that contains only the clines that are also in significant windows

I suggest making a separate table of only clines with significant cline centers. To pull only clines with significant p-values, use:

```
cat cline_table_r25.windows.tsv | grep -v 'g' | awk '$9 < 0.05 {print $0}' > significant_cline_centers_in_windows.tsv
```

Next, use the BASH script [count_clines.sh](count_clines.sh) to count your total number of clines, number of significant clines (using an alpha of 0.05), and number of non-sifnificant clines, in both your cline_table.txt for totals and in significant_cline_centers_in_windows.tsv to get counts for inside your genomic windows of interest. 

Now that you have these counts, you can plug these into a 2x2 matrix and do a Fisher's Exact test in R to see if the number of significan't genomic cline centers are enriched inside your define genome windows, versus matching the general amount of significant clines reguardless of the windows. I have a messy script that I used to perform the Fisher's Exact test as [Fishers_exact_enrichment_cline.R](Fishers_exact_enrichment_clines.R). 
