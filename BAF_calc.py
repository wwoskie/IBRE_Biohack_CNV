#!/usr/bin/env python3


import argparse

import numpy as np
import pandas as pd


def calc_baf(cov_ref, cov_alt):
    '''A function for calculating BAF (B-allele frequency) based on SNP coverage
    
    Parameters:
        cov_ref (float): coverage of reference allele
        cov_alt (float): coverage of alternative allele
        
    Returns: BAF value (float) or "-" (str) in case of zero coverage for both reference and alternative SNP
    '''
    
    
    return cov_alt / (cov_alt + cov_ref) if (cov_alt + cov_ref) else "-"


def calc_copy_num(baf_tum):
    '''A function for calculating copy number based on BAF values for tumor sample
    
    Parameters:
        baf_tum (float): BAF value for SNP from tumor sample
    
    Returns: 1) (str) "-" in case of absent BAF;
             2) (int) copy number;
             3) (int) 1 in case of missing changes in gene copies for tumor sample
    '''
    
    
    if baf_tum == "-":
        return "-"
    elif 0 < baf_tum < 1:
        copy_num = 2 ** abs(np.log2(baf_tum / (1 - baf_tum))) + 1
        return round(copy_num)
    else:
        return 1
    
    
def calc_cn_gene(interest_genes_df):
    '''A function for calculating mean copy number for genes of interest
    
    Parameters:
        interest_genes_df (pd.DataFrame): dataframe with coordinates for gene of interest
    
    Returns: 
        result (pd.DataFrame): DataFrame with mean copy number for genes of interest
    '''
    
    
    result = pd.DataFrame({"gene": [], "copy_num_mean": [], "n": []})
    for i in range(interest_genes.shape[0]):
        chrom, start, stop, gene = interest_genes.loc[i, :]
        result.loc[i, "gene"] = gene
        subset = pileup_hetero[
            (pileup_hetero["Chromosome"] == chrom) & 
            pileup_hetero["Position"].between(start, stop)
        ]
        subset = subset[subset["copy_num"] != "-"]
        result.loc[i, "copy_num_mean"] = subset["copy_num"].mean()
        result.loc[i, "n"] = subset["copy_num"].shape[0]
        result.replace(np.nan, "-", inplace=True)
    result.to_csv(args.output, sep="\t", index=False)
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-c", "--cov", help="path to FACETS pileup file")
    parser.add_argument("-b", "--bed-genes", help="path to .bed file with genes of interest")
    parser.add_argument("-o", "--output", help="path for output")
    args = parser.parse_args()
    
    # reading file with SNP coverage for normal and tumor sample (FASECTS output)    
    pileup = pd.read_csv(args.cov)
    # reading .bed file with genes of interest
    interest_genes = pd.read_csv(args.bed_genes, sep="\t", names=["chr", "start", "stop", "gene"])

    # filtering heterozygous variants
    pileup_hetero = pileup[(pileup["File1R"] != 0) & (pileup["File1A"] != 0)]
    pileup_hetero = pileup_hetero[abs(np.log2(pileup_hetero["File1R"] / pileup_hetero["File1A"])) <= 0.322]

    # calculating BAF values for normal and tumor sample
    pileup_hetero.loc[:, "baf_contr"] = pileup_hetero.apply(lambda x: calc_baf(x["File1R"], x["File1A"]), axis=1)
    pileup_hetero.loc[:, "baf_tum"] = pileup_hetero.apply(lambda x: calc_baf(x["File2R"], x["File2A"]), axis=1)

    # calculating copy number for tumor sample based on BAF
    pileup_hetero.loc[:, "copy_num"] = pileup_hetero["baf_tum"].apply(calc_copy_num)

    result_df = calc_cn_gene(interest_genes)

# data/pileup.txt
# data/interesting_genes.bed
# data/BAF_output.txt