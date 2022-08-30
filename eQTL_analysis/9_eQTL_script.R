#########################eQTL analysis##########################################


#This scripts is used to carry out the eQTL analysis for all the genes in which 
# significant SNPs were found using their corresponding GTex files  and 
# the adjusted genomic positions from the gene table (dissertation Table 1).
# Requires R  version 4.0.0


##### Load pakages
library(tidyverse)
library(eQTpLot)

##### Set the working directory 
setwd("/Users/christianarenasdaza/Desktop/lab_ch2/input/9_eQTL")



##### Load files 

# eQTL files 
LTA4H_expression <- read.table(file = "GTEx Portal_LTA4H.csv",
                               header=TRUE,
                               sep= ",") 

GPR18_expression <- read.table(file = "GTEx Portal_GPR18.csv",
                               header=TRUE,
                               sep= ",")

RORA_expression <- read.table(file = "GTEx Portal_RORA.csv",
                               header=TRUE,
                               sep= ",") 

DPEP1_expression <- read.table(file = "GTEx Portal_DPEP1.csv",
                              header=TRUE,
                              sep= ",") 

# Load logistic regression results
assoc_results <- read.table(file = "/Users/christianarenasdaza/Desktop/lab_ch2/output/6_QC_and_association_test/SPMsRG_logistic_add_results.tsv",
                            header = TRUE)

# Load significant SNPs
assoc_results <- read.table(file = "/Users/christianarenasdaza/Desktop/lab_ch2/output/6_QC_and_association_test/SPMsRG_logistic_add_results.tsv",
                            header = TRUE)


# Gene table: wiht candidate genes adjusted positions in bp
Genes.df <- read.table(file = "genes_37_adjusted_kb.tab",
                            header = TRUE,
                            sep = "\t")


sig_SNP_FDR <- read.table(file = "/Users/christianarenasdaza/Desktop/lab_ch2/output/6_QC_and_association_test/Sig_SNP_pass_FDR.tsv",
                       header = TRUE,
                       sep = "\t")

##### Data formatting 
# Join all the gene expression files into one 
sig_genes_expression <-  bind_rows(list(LTA4H_expression, GPR18_expression, 
                                           RORA_expression, DPEP1_expression))

table(duplicated(sig_genes_expression))
# Select only the columns needed  for gene expression and changed their location
sig_genes_expression <-sig_genes_expression %>%
  select( !c(Gencode.Id, Variant.Id)) %>% 
  relocate(Gene.Symbol, .after = SNP.Id)


#  Select only the columns needed from the  association results and changed their location
assoc_results <- assoc_results %>%  
  relocate(BP, .before = SNP) %>% 
  relocate(P, .before = BETA)


##### Generating eQTL figures 

# Statistics 
eQTpLot(GWAS.df = assoc_results, eQTL.df = sig_genes_expression, 
        gene = c("LTA4H", "GPR18", "RORA", "DPEP1"), 
        gbuild = "hg19", trait = "Asthma", tissue =  "all", 
        CollapseMethod = "min", GeneList = T, sigpvalue_GWAS= 1.401e-04)

# eQTL LTA4H
eQTpLot(GWAS.df = assoc_results, eQTL.df = sig_genes_expression, gene = "LTA4H", 
        gbuild = "hg19",  trait = "Asthma", tissue =  "all", 
        CollapseMethod = "min", sigpvalue_GWAS= 1.401e-04, ylima = 5) 



# eQTL GPR18
eQTpLot(GWAS.df = assoc_results, eQTL.df = sig_genes_expression, gene = "GPR18", 
        gbuild = "hg19",  trait = "Asthma", tissue =  "all", 
        CollapseMethod = "min", sigpvalue_GWAS= 1.401e-04)

# eQTL RORA
eQTpLot(GWAS.df = assoc_results, eQTL.df = sig_genes_expression, gene = "RORA", 
        gbuild = "hg19",  trait = "Asthma", tissue =  "all", 
        CollapseMethod = "min", sigpvalue_GWAS= 1.401e-04)

# eQTL DPEP1; lies within the significant threshold 
eQTpLot(GWAS.df = assoc_results, eQTL.df = sig_genes_expression, gene = "DPEP1", 
        gbuild = "hg19",  trait = "Asthma", tissue =  "all", 
        CollapseMethod = "min", sigpvalue_GWAS= 1.401e-04)


##### Generate table for significant eQTLs

eQTL_table <- assoc_results %>%
  rename(SNP.Id = SNP) %>% 
  merge(sig_genes_expression, by.x = "SNP.Id") %>% 
  select(CHR, SNP.Id, BP, Gene.Symbol, A1, A2, P, P_adjusted_BH, 
         P.Value, NES, Tissue) %>% 
  rename(P_eQTL = P.Value) %>% 
  filter(P_eQTL< 0.05) %>% 
  filter(P < 1.401e-04) %>% 
  filter(SNP.Id %in% sig_SNP_FDR$SNP == TRUE) %>% 
  distinct(SNP.Id, .keep_all = TRUE) %>% 
  mutate(P_log = -log10(P)) %>% 
  relocate(P_log, .after = P)  %>% 
  rename( SNP = SNP.Id, Gene = Gene.Symbol)
  

##### Write significant SNPs that are eQTLs in a file 
write.table(eQTL_table, 
            file = paste("/Users/christianarenasdaza/Desktop/lab_ch2/output/9_eQTL/Sig_SNP_eQTL_table.tsv", sep = ""),
            sep = "\t",
            quote = FALSE,
            row.names = FALSE)



  

