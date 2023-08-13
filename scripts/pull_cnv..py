import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

def compare_intervals(bounds, interval):
    # In this code, the function compare_intervals() receives two arguments: 

# - bounds - It represents the boundary in which the interval should lie.
# - interval - It represents the range of values which should lie within the bounds.

# It returns the interval if it's strictly within the bounds. But if the interval is larger than the bounds at any end (or both ends), then it will return the bounds as the interval. If the entire interval is out of the bounds, then it will return None.
    assert isinstance(bounds, tuple) and len(bounds) == 2
    assert isinstance(interval, tuple) and len(interval) == 2

    lower_bound, upper_bound = bounds
    lower_interval, upper_interval = interval

    if lower_interval > upper_bound or upper_interval < lower_bound:
        # Whole interval is outside bounds
        return None
    
    elif lower_interval < lower_bound and upper_interval > upper_bound:
        # Interval crosses both bounds
        return (lower_bound, upper_bound)
    
    elif lower_interval < lower_bound:
        # One value of interval is outside bounds
        return (lower_bound, upper_interval)
    
    elif upper_interval > upper_bound:
        # One value of interval is outside bounds
        return (lower_interval, upper_bound) 
   
    else:
        # Interval is fully inside bounds
        return interval


def pull_from_cnvkit(path, gene_chrom, bounds):
    segments = []
    depths = []
    with open(path) as f:
        f.readline()
        while True:
            line = f.readline().split()
            if len(line) == 0:
                break
            chrom, lower, upper, _, depth,*_ = line

            if chrom == gene_chrom:
                chunk = compare_intervals(bounds, (int(lower), int(upper)))
                if chunk:
                    segments.append(chunk)
                    depths.append(float(depth))
                else:
                    continue
    return segments, depths


if __name__ == "__main__":
    
    cnvkit_path = 'WES_Tumor.cns'
    genes_path = 'interesting_genes.bed'
    # sequenza_path
    with open(genes_path) as f:
        genes = f.readlines()
    genes = [gene.split() for gene in genes]
    for gene in genes:
        gene_chrom, gene_lower, gene_upper, name = gene
        print(gene)
        segments, depths = pull_from_cnvkit(cnvkit_path, gene_chrom, tuple(
            [int(gene_lower), int(gene_upper)]))
        for i, segment in enumerate(segments):
            plt.hlines(depths[i], segment[0], segment[1])
        plt.show()