This repository contained all the R scripts I used for my MSc Bioinformatics project; Candidate gene analysis for the association between genetic mutations in specialised pro-resolving mediators(SPMs) and asthma. If acute inflammation is dysregulated such as in unresolved resolution  and it becomes chronic. SPMs are mediator  produced  by a number of enzymes in last phase of inflammation, called resolution. SPMs bind to receptors to increase pro-resolving process and decrease the pro-inflammatory process. Asthma is a chronic inflammatory condition of the airways, which  studies have shown lower levels of SPMs and decreased expression of SPMs biosynthesis enzymes and receptors' expression. Therefore, this project wanted to explore if genetic variants that happen in SMPs enzymes and receptor's  genes are linked to asthma and how could they be contributing to the asthma phenotype via SPMs. To achieve the aim of this project I did; Gene candidate analyses adjusted for several asthma covariates to test for association, meta-analysis for the discovery of new SNPs and validate significant SNPs and eQTL analysis to find how significant variants could be contributing to asthma.


## Covariates
The 1.5_covariates_files_ch_final.R script was used to format my chosen covariates for logistic regression in PLINK1.9 due to  UKBiobank re-measurements of variables (instances) and PLINK1.9 requiring a single column per level of categorical variable.


## Sampel selection 
The sample selection scripts were used to visualise the age, ethnicity, sex and smoking variables for controls (4_control_demography_R-redo.R) and cases (3_Demographic_analysis.R). These scripts were also used to select the asthma cases and healthy controls, with the following demographic criteria; 52 years or older, caucasian (British, Irish, white and any other white background) and all genotyped chromosomes. The 5_cases_vs_controls.R  script was used to combine all IDs (control and cases) and add an additional column specifying 1 for controls and 2 for cases.


## CGCA
The 6_GC_results_formating_and_plots.R script takes in the logistic regression and allele frequency result, combines them, subsets for significant SNPs and produces a manhattan plot.

## Meta_analysis 
The 10_validation.R script identified the GABRIEL consortium's SNPs positioned in my candidate genes as specified by the adjusted locations in my gene table. The 8_Meta_analysis_data_formatting.R script formatted my logistic regression results and the GABRIEL consortium data, which includes but is not limited to re-naming and selecting columns. The 8_Meta_analysis_results_Final.R script, takes the meta-analysis results and creates a forest plot for significant SNPs.


## eQTL
The 9_eQTL_script.R script was used to perform and visualise eQTL analysis for genes with significant SNPs; GPR18 RORA and LT4H, as well as DPEP1 (sits in the significance level) based on eQTL data from GTEx portal V8.  


## Other Scripts 
Other scripts used for this project that are not in this repository include; subsetting the UKBiobank column for desire variables, candidate gene analysis PLIN1.9 and matching SNPs of the GABRIEL consortium data set to their corresponding position to the hg19 reference genome.

## Acknowledgements
This project would have not been possible without the help and advice of my supervisors Jesmond Dalli and  Esteban Alberto Gomez Cifuentes (an amazing storyteller). Also, a huge thanks to the amazing support of the previous and current members of the JDalli lab team.



## References 











