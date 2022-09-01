# MSc Bioinformatics research project 
This repository contained all the R scripts I used for my MSc Bioinformatics project; Candidate gene analysis for the association between genetic mutations in specialised pro-resolving mediators(SPMs) and asthma. If acute (normal) inflammation (Robbins et al., 1999) is dysregulated such as in unresolved resolution, it becomes chronic (Medzhitov, 2010; Serhan, 2014). SPMs are mediator  produced  by a number of enzymes in last phase of inflammation, called resolution (llerton and Gilroy, 2016; Chiang and Serhan, 2020). SPMs bind to receptors to increase pro-resolving process (Chiang and Serhan, 2020) and decrease the pro-inflammatory process (Chiang and Serhan, 2020; Díaz del Campo et al., 2022). Asthma is a chronic inflammatory condition of the airways (Pulleyn et al., 2001; Vicente et al., 2017), which  studies have shown lower levels of SPMs (Levy et al., 2006; Planagumà et al., 2008) and decreased expression of SPMs biosynthesis enzymes and receptors' expression (Planagumà et al., 2008). Therefore, this project wanted to explore if genetic variants that happen in SMPs enzymes and receptor's genes are linked to asthma and how could they be contributing to the asthma phenotype. To achieve the aim of this project I did; Gene candidate analyses adjusted for several asthma covariates to test for association, meta-analysis for the discovery of new SNPs and validate significant SNPs and eQTL analysis to find how significant variants could be contributing to asthma.


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
1. Robbins, S., Cotran, R., Kumar, V. and Collins, T., 1999. Pathologic basis of disease. 6th ed. Philadelphia [etc.]: W. B. Saunders.
2. Medzhitov, R., 2010. Inflammation 2010: New Adventures of an Old Flame. Cell, 140(6), pp.771-776.
3. Serhan, C., 2014. Pro-resolving lipid mediators are leads for resolution physiology. Nature, 510(7503), pp.92-101.
4. 11.	Fullerton, J. and Gilroy, D., 2016. Resolution of inflammation: a new therapeutic frontier. Nature Reviews Drug Discovery, 15(8), pp.551-567.
5. Chiang, N. and Serhan, C., 2020. Specialized pro-resolving mediator network: an update on production and actions. Essays in Biochemistry, 64(3), pp.443-462.
6.	Díaz del Campo, L., Rodrigues-Díez, R., Salaices, M., Briones, A. and García-Redondo, A., 2022. Specialized Pro-Resolving Lipid Mediators: New Therapeutic Approaches for Vascular Remodeling. International Journal of Molecular Sciences, 23(7), p.3592.
7.	29.	Pulleyn, L., Newton, R., Adcock, I. and Barnes, P., 2001. TGFβ1 allele association with asthma severity. Human Genetics, 109(6), pp.623-627.
8.	45.	Vicente, C., Revez, J. and Ferreira, M., 2017. Lessons from ten years of genome-wide association studies of asthma. Clinical &amp; Translational Immunology, 6(12), p.e165.
9.	Levy, B., Kohli, P., Gotlinger, K., Haworth, O., Hong, S., Kazani, S., Israel, E., Haley, K. and Serhan, C., 2006. Protectin D1 Is Generated in Asthma and Dampens Airway Inflammation and Hyperresponsiveness. The Journal of Immunology, 178(1), pp.496-50.
10.	27.	Planagumà, A., Kazani, S., Marigowda, G., Haworth, O., Mariani, T., Israel, E., Bleecker, E., Curran-Everett, D., Erzurum, S., Calhoun, W., Castro, M., Chung, K., Gaston, B., Jarjour, N., Busse, W., Wenzel, S. and Levy, B., 2008. Airway Lipoxin A4Generation and Lipoxin A4 Receptor Expression Are Decreased in Severe Asthma. American Journal of Respiratory and Critical Care Medicine, 178(6), pp.574-582.
