#!/usr/bin/env python3
# (c) Angel Rivera-Colon - November 2020
import argparse, os

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('-s', '--sumstats',        required=True, help='Stacks `populations.sumstats.tsv` file')
    parser.add_argument('-w', '--genomic-windows', required=True, help='Genomic windows table')
    parser.add_argument('-c', '--cline-table',     required=True, help='Cline table from GGHybrid')
    parser.add_argument('-o', '--out-tsv',         required=True, help='Out table name')
    args = parser.parse_args()
    return args

# Class for the Stacks Sumstats SNP
class SumstatsSNP:
    def __init__(self, locus_id, column, chrom, bp):
        assert type(locus_id) is int
        assert type(column) is int
        assert type(bp) is int
        self.loc_id = locus_id
        self.col    = column
        self.chrom  = chrom
        self.bp     = bp
        self.snp_id = f'{self.loc_id}_{self.col}'
    def __str__(self):
        return f'{self.snp_id} {self.chrom} {self.bp}'

# Create a dictionary for the genomic windows
#
# Structure of the genome windows file:
# scaffold<tab>start<tab>end
# NW_021940813.1<tab>9800001<tab>9810000
# NW_021940813.1<tab>10400001<tab>10410000
# All the Golden Collared scaffolds start with NW
def genomic_windows_dict(genome_windows_file):
    assert os.path.exists(genome_windows_file) is True, 'Genomic Windows File not found'
    # Store the desired windws for later testing
    geno_windows_dict = dict()
    for line in open(genome_windows_file, 'r'):
        if line[0] == '#':
            continue
        if line[0:2] != 'NW':
            continue
        fields = line.strip('\n').split('\t')
        chrom = fields[0]
        start_bp = int(fields[1])
        end_bp = int(fields[2])
        geno_windows_dict.setdefault(chrom, list()).append((start_bp, end_bp))
    return geno_windows_dict

# Function to match the genomic window dictionary
def snp_in_windows(chrom, bp, geno_windows_dict):
    assert type(bp) is int, 'Base Pair is not Int'
    assert type(geno_windows_dict) is dict
    chrom_positions = geno_windows_dict.get(chrom, False)
    if chrom_positions is False:
        return False
    else:
        for window in chrom_positions:
            start = window[0]
            end = window[1]
            if start <= bp <= end:
                return True
    return False

# Function to filter the sumstats to retain loci/SNPs in the windows
def filter_sumstats_by_windows(sumstats_file, geno_windows_dict):
    assert os.path.exists(sumstats_file) is True, 'Sumstats file not found'
    assert type(geno_windows_dict) is dict
    # Sumstats file header fields
    # 00 - Locus ID
    # 01 - Chr
    # 02 - BP
    # 03 - Col
    # 04 - Pop ID
    # 05 - P Nuc
    # 06 - Q Nuc
    # 07 - N
    # 08 - P
    # 09 - Obs Het
    # 10 - Obs Hom
    # 11 - Exp Het
    # 12 - Exp Hom
    # 13 - Pi
    # 14 - Smoothed Pi
    # 15 - Smoothed Pi P-value
    # 16 - Fis
    # 17 - Smoothed Fis
    # 18 - Smoothed Fis P-value
    # 19 - HWE P-value
    # 20 - Private
    #
    # Parse Sumstats file
    snps_in_windows_dict = dict()
    for line in open(sumstats_file, 'r'):
        if line[0] == '#':
            continue
        if line[0:5] == 'Locus':
            # Skip header
            continue
        fields = line.strip('\n').split('\t')
        chromosome = fields[1]
        basepair = int(fields[2])
        locus_id = int(fields[0])
        column = int(fields[3])
        if snp_in_windows(chromosome, basepair, geno_windows_dict) is False:
            continue
        sumstats_snp = SumstatsSNP(locus_id, column, chromosome, basepair)
        if sumstats_snp.snp_id not in snps_in_windows_dict:
            snps_in_windows_dict[sumstats_snp.snp_id] = sumstats_snp
    return snps_in_windows_dict

# Parse the GGhybrid cline table
# Only return SNPs in the windows
# Add the chromosome and basepair coordinates to the output table
def parse_cline_table_by_windows(cline_table, snps_in_windows_dict, output_file):
    assert os.path.exists(cline_table) is True, 'Cline table file not found'
    assert type(snps_in_windows_dict) is dict
    outf = open(output_file, 'w')
    # Open file and loop over lines
    for line in open(cline_table, 'r'):
        if line[0] == '#':
            continue
        elif line[0] == 'g':
            # If the line is a header, reformat into TSV and add chromosome and bp
            header = line.strip('\n').split(' ')
            header_str = '\t'.join(header) + '\tchromosome\tbp\n'
            outf.write(header_str)
        else:
            fields = line.strip('\n').split(' ')
            snp_id = fields[0]
            # Check if the SNP is in the windows
            if snp_id in snps_in_windows_dict:
                sumstats_snp = snps_in_windows_dict[snp_id]
                line_str = '\t'.join(fields)
                row_str = f'{line_str}\t{sumstats_snp.chrom}\t{sumstats_snp.bp}\n'
                outf.write(row_str)
    outf.close()

# Main Function
def main():
    # Parse Arguments
    args = parse_args()
    # Arguments
    windows = args.genomic_windows
    sumstats = args.sumstats
    clines = args.cline_table
    outfile = args.out_tsv
    # Get genomic windows dictionary
    geno_windows_dict = genomic_windows_dict(windows)
    # Parse and filter sumstats
    snps_in_windows_dict = filter_sumstats_by_windows(sumstats, geno_windows_dict)
    # Parse cline table and print SNPs in windows
    parse_cline_table_by_windows(clines, snps_in_windows_dict, outfile)


# --------
# Run Code
# --------

if __name__ == '__main__':
    main()
