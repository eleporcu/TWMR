##########Eleonora Porcu - Mendelian Randomization integrating GWAS and eQTL data reveals genetic determinants of complex and clinical traits


MR.R requires two input files:
- a matrix containing the univariate effect size of n SNPs on k gene expressions (these estimates come from an eQTL study) and the univariate effect sizes on the phenotype
- the LD matrix between the n SNPs

This folder includes the matrixes used to perform MR analyses for height using the eQTLs results downloaded from https://molgenis26.gcc.rug.nl/downloads/eqtlgen/cis-eqtl/cis-eQTL_significant_20181017.txt.gz (Vosa et al, Biorxiv 2018) and the Height-GWAS results downloaded from https://grasp.nhlbi.nih.gov/downloads/ResultsOctober2016/Wood/GIANT_HEIGHT_Wood_et_al_2014_publicrelease_HapMapCeuFreq.txt.gz (Wood et al, Nat Gen 2014). 

How to run the script:

R < MR.R --no-save ENSGXXXXXXXX

where ENSGXXXXXXXX is the name of the gene to be analysed (e.g.: ENSG00000000419)

The output files is a .txt file containing the following columns:
-gene: name of the gene tested 
-alpha: causal effect estimated by MR
-SE: standard errot of the causal effect
-P: Pvalue 
-Nsnps: number of SNPs included in the model
-Ngene: number of genes included in the model 







