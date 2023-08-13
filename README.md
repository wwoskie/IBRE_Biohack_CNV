# IBRE_Biohack_CNV

## Calculation gene copy number using BAF values

Python script [BAF_calc.py](./BAF_calc.py) is used for calculation copy numbers for genes of interest. It uses output from [FACETS](https://github.com/mskcc/facets) algorithm with SNP coverage). To run this script, you need **two files**:

1) pileup.txt: output from FACETS with SNP coverage;
2) .bed file: with genome coordinates for gene of interest

Arguments:

```
-c, --cov: path to FACETS pileup file;
-b, --bed-genes: path to .bed file with genes of interest coordinates;
-o, --output: path for output
 
```  

Here is the example command for running this script:

```
python3 BAF_calc.py --cov data/sample1_pileup.txt 
		    --bed-genes interesting_genes.bed 
		    --output sample1_BAF.txt
```

The output will be generated in this format:

| gene | copy_num_mean | n |
|--------------|:-------:|:--------:|
| ERBB2 | - | 0 |
| BRCA1	| 5 | 8 |
| BRCA2	| 2 | 2 |
| CHEK1	| 1 | 1 |
| CHEK2	| - | 0 |
| EGFR	| 18 | 2 |
| ... | ... | ... |
