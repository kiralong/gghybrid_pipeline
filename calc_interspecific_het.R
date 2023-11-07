#!/usr/bin/env Rscript

# Calculate interspecific heterozygosity from Manacus hybrids
setwd("/path/to/working/dir")

# Load input packages
library(genetics)
library(introgress)

# data("AdmixDataSim2")
# data("LociDataSim2")

# ===============
# Genotypes table
# ===============
geno_f <- './Long_pops.p9.r80.sSNP/genotypes_tbl.csv' # Long
geno_f <- './Brum_pops.p9.r80.sSNP/genotypes_tbl.csv' # Brum
geno_df <- read.delim(geno_f, sep=',', header=F)

# ==========
# Loci table
# ==========
loci_f <- './Long_pops.p9.r80.sSNP/loci_tbl.csv' # Long
loci_f <- './Brum_pops.p9.r80.sSNP/loci_tbl.csv' # Brum
loci_df <- read.delim(loci_f, sep=',')

# ================
# Output file name
# ================
out_f <- './Long_pops.p9.r80.sSNP/manacus_long.intersp_het.csv' # Long
out_f <- './Brum_pops.p9.r80.sSNP/manacus_brum.intersp_het.csv' # Brum

# ======================================================
# Split the base genotypes df across hybrids and parents
# ======================================================

parent1_id <- '100CG'
parent2_id <- '020SS'

# Get parent 1 genotype data
parent1_df <- geno_df[-c(1,2),c(geno_df[1,] == parent1_id)]
parent1_sams <- as.character(geno_df[2,c(geno_df[1,] == parent1_id)])
  
# Get parent 2 genotype data
parent2_df <- geno_df[-c(1,2),c(geno_df[1,] == parent2_id)]
parent2_sams <- as.character(geno_df[2,c(geno_df[1,] == parent2_id)])

# Get hybrids genotype data
hybrids_df <- geno_df[,c(geno_df[1,] != parent1_id)]
hybrids_df <- hybrids_df[,c(hybrids_df[1,] != parent2_id)]
hybrids_sams <- as.character(hybrids_df[2,])

# ==================================
# Run initial formatting of the data
# ==================================

fmt_genos <- prepare.data(
  admix.gen = hybrids_df,
  parental1 = parent1_df,
  parental2 = parent2_df,
  loci.data = loci_df,
  pop.id = TRUE,
  ind.id = TRUE,
  fixed = FALSE
)


# ======================================
# Calculate interspecific heterozygosity
# ======================================

interspec_het <- calc.intersp.het(fmt_genos)


# =============================
# Format the final output table
# =============================

# Merge all the samples
all_samples <- c(parent1_sams, hybrids_sams, parent2_sams)
# Get a vector of population ids
all_pops <- sapply(strsplit(all_samples, '_'),
                   function(all_samples){paste(all_samples[2],
                                               toupper(all_samples[1]),
                                               sep='')})
# Merge all the insterspecies het values, parents are 0.0 
all_hets <- c(rep(0.0, length(parent1_sams)),
              interspec_het,
              rep(0.0, length(parent2_sams)))

out_df <- data.frame(all_samples, all_pops, all_hets)
colnames(out_df) <- c('SampleID', 'PopID', 'InterSpHet')

write.csv(out_df, out_f, quote=F, eol='\n', row.names=F)

# plot(out_df$InterSpHet)



