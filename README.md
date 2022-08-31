This repository contained all the R scripts I used for my MSc Bioinformatics project; Candidate gene analysis for the association between genetic mutations in specilised pro-resolving mediators (SPMs) and asthma. SPMs are chemicals from the last phase of inflammation called resolution. These mediators bind to receptors and increase pro-resolving process and decrrease pro-inflammatory process. Humna asmthma   studies have shown,  lower levels of SPMs and decresed expression of SPMs bysynthesis and receptor expression. Since asthma  is a chronnic inflammatory condition  of the airway and an chronic inflammation  happens due to dysregulatiesd inflmamation, suhc as in the case of failed resolution.  Therefore, this project wanted to explore if genetic varinats that happen in SMPs enzyme receptor genes, could affect the production or and action of SPMs, which then would contribute to asthma. To achive the aim of this project I did. Gene candidate analyis  adjusted for several asthma relevant covariates, to find if varient were associated to asthma, meta analysis; to validate significant SMPs and discovere nre SNPs and wQTL analyis to idetify dignificant SNPs that altere gene expression.



## Covariates
This R script formatted my chosen covariates for  logistic regresion in PLINK1.9, which were chosen based on UKBiobank recommendation and the literature on asthma. 

## Sampel selection 
The sample selection scripts were used to visulized the age, ethnicty, sex and smooking variables for control and cases. These scripts were also used to select the asthma cases and healhty controls, which were 52 years or older, cacasian (British, Irish, white and any other white background) and had all genotyped chromosomes.

## CGCA
This cript takes in the logistic regression and  allele frequency result, from whihc a manhattan plot is created and the signifcant snps were found 


## Meta_analysis 

Validation using  finds the GABRIEL consortium SNPs that are in the cnadidate gene wihtin the adjusted locations especified in the gene  andidate table  and also produce a mahanttan pplot of those results.

The formatting scripts was used to format my logistic regression results and the the GABRIEL consortium data, whihc include  but not limited to re-naming and selecting column.

The analysis scripts  takes the meta-analyis results, subsets and create a forest plot for significant SNPs.


## eQTL

This script was used to performed eQTL analyis for gene wiht significan SNPs GPR18 RORA and LT4H, as wll as DPEP1. Although non significant it settel in the significnat lebel.



## Other Scripts 

Other scripts used for this project that are not repository include fro subsetting UKBibank column for desire variables, candidate gene anlyiss PLIN1.9 and  mathcing SNPs of the GABRIEL consortium data set to their corresponding position in the hg19 reference genome




## Acknowledgements
This project would have not being possible wihthout the help and advised of Jesmond Dalli and Estebal alberto Gomez Cifuence  and the amazing support of the JDalli lab team.













