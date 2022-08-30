########################### Meta analysis #####################################

# This script is use for formatting the logistic regression results from PLINK1.9
# and the GABRIEL data set for meta analysis using METAL.


##### Load libraries
library(tidyverse)
library(ggplot2)
library(ggrepel)


##### set directory

setwd("/Users/christianarenasdaza/Desktop/lab_ch2/input/8_Meta-analysis")

##### Load data 
Table_logistic_add_results <- read.table("/Users/christianarenasdaza/Desktop/lab_ch2/output/6_QC_and_association_test/SPMsRG_logistic_add_results.tsv",
                    header=TRUE)

validation <- read.table(file = paste("/Users/christianarenasdaza/Desktop/lab_ch2/input/10_validation/GABRIEL.tsv", sep = ""),
                         header=TRUE, 
                         sep="\t")

# SNPs in GABRIEL data set located within my candidate genes 
validation_SNP_subset <- read.table(file = paste("/Users/christianarenasdaza/Desktop/lab_ch2/output/10_validation/Validation_SNP_Subset.tsv", sep = ""),
                         header=TRUE, 
                         sep="\t")


##### Data formatting 

# Relocate columns from logistic results
Logistic_add_METAL <- Table_logistic_add_results %>% 
  relocate(c(A2, MAF_A, MAF_U), .after = A1) %>% 
  rename(N_cases= NCHROBS_A,  N_controls = NCHROBS_U)


# Combine SNP subset with the validation  data set to recover missing columns
val_final_METAL <- validation_SNP_subset %>%  
  rename(rs = SNP) %>%  
  merge(validation, by = "rs") %>% 
  filter(is.na(P_ran) == FALSE) %>% 
  rename(A1 = Allele_1, A2 = Allele_2, SNP = rs) %>% 
  select(CHR, SNP, BP, Gene, A1, A2, OR_ran, P_ran) %>% 
  mutate(N_cases = 10365) %>%  # Number of cases in GABRIEL data set
  mutate(N_controls = 16110) # Number of controls in GABRIEL data set


##### Write METAL input files

write.table(Logistic_add_METAL, 
            file = paste("/Users/christianarenasdaza/Desktop/lab_ch2/output/8_Meta-analysis/SPMsRG_logistic_add_METAL_input.tsv", sep = ""),
            sep = "\t",
            quote = FALSE,
            row.names = FALSE)


write.table(val_final_METAL, 
            file = paste("/Users/christianarenasdaza/Desktop/lab_ch2/output/8_Meta-analysis/Validation_METAL_input.tsv", sep = ""),
            sep = "\t",
            quote = FALSE,
            row.names = FALSE)


