# IBRE_Biohack_CNV
This project is inspired by IBRE Biohack hackaton problem on CNV DAta analisys.

Gains and losses of DNA are prevalent in cancer and emerge as a consequence of inter-related processes of replication stress, mitotic errors, spindle multipolarity and breakage–fusion–bridge cycles, among others, which may lead to chromosomal instability and aneuploidy.

## Causes:
- non-homologous end joining;
- template switching;
- chromosome breakage during mitosis -> chromothripis;
- replication timing:
- early replication -> unmatched DNA regions can be mistakenly combined;
- late replication -> nonhomologous recombination;
- specific chromatin remodeling events

## Two major ways for CNA calling

### BAF (B-allele frequency)
“B” allele is the non-reference allele observed in a germline heterozygous SNP, i.e. in the normal/control sample. Estimation from a SNP b-allele frequencies works by comparing the shift in allele frequency of heterozygous, germline SNPs in the tumor sample from the expected ~50%.

~BAF = coverage Alt / (coverage Alt + coverage Ref)~


