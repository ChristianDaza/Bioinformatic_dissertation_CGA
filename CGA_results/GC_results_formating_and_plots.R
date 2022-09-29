############analysis of PLINK1.9 Logistic regression results####################

# This script is used for:
 # Loading PLINK 1.9 additive model results.
 # Combine results into a single data frame. 
 # Clean duplicates and SNPs missing P values. 
 # Apply multiple  testing correction.
 # Create data frame of significant SNPs.
 # Create Manhattan and QQ plots.

##### Load packages  
library(tidyverse)
library(ggplot2)
library(ggrepel)


##### Set directorie
setwd("/Users/christianarenasdaza/Desktop/lab_ch2/input/6_QC_and_association_test")

##### Read files: Logistic results and frequecy results for each gene 

# Logistic regression files 
GSTM4 <- read.table("ASTHM_vs_Control_Beta2/1_all/GSTM4/log_chr1_GSTM4_ASTHM.assoc.logistic",
                    header=TRUE)

LGR6 <- read.table("ASTHM_vs_Control_Beta2/1_all/LGR6/log_chr1_LGR6_ASTHM.assoc.logistic",
                   header=TRUE)

mGST3 <- read.table("ASTHM_vs_Control_Beta2/1_all/mGST3/log_chr1_mGST3_ASTHM.assoc.logistic",
                    header=TRUE)

PTGS2 <- read.table("ASTHM_vs_Control_Beta2/1_all/PTGS2/log_chr1_PTGS2_ASTHM.assoc.logistic",
                    header=TRUE)

mGST2 <- read.table("ASTHM_vs_Control_Beta2/4_all/mGST2/log_chr4_mGST2_ASTHM.assoc.logistic",
                    header=TRUE)

LTC4S <- read.table("ASTHM_vs_Control_Beta2/5_all/LTC4S/log_chr5_LTC4S_ASTHM.assoc.logistic",
                    header=TRUE)

AHR  <- read.table("ASTHM_vs_Control_Beta2/7_all/AHR/log_chr7_AHR_ASTHM.assoc.logistic",
                   header=TRUE)

GPR37  <- read.table("ASTHM_vs_Control_Beta2/7_all/GPR37/log_chr7_GPR37_ASTHM.assoc.logistic",
                     header=TRUE)

EPHX2  <- read.table("ASTHM_vs_Control_Beta2/8_all/EPHX2/log_chr8_EPHX2_ASTHM.assoc.logistic",
                     header=TRUE)

ALOX5  <- read.table("ASTHM_vs_Control_Beta2/10_all/ALOX5/log_chr10_ALOX5_ASTHM.assoc.logistic",
                     header=TRUE)

CMKLR1 <- read.table("ASTHM_vs_Control_Beta2/12_all/CMKLR1/log_chr12_CMKLR1_ASTHM.assoc.logistic",
                     header=TRUE)

LTA4H <- read.table("ASTHM_vs_Control_Beta2/12_all/LTA4H/log_chr12_LTA4H_ASTHM.assoc.logistic",
                    header=TRUE)

CYSLTR2 <- read.table("ASTHM_vs_Control_Beta2/13_all/CYSLTR2/log_chr13_CYSLTR2_ASTHM.assoc.logistic",
                      header=TRUE)

GPR18 <- read.table("ASTHM_vs_Control_Beta2/13_all/GPR18/log_chr13_GPR18_ASTHM.assoc.logistic",
                    header=TRUE)

GPR99 <- read.table("ASTHM_vs_Control_Beta2/13_all/GPR99/log_chr13_GPR99_ASTHM.assoc.logistic",
                    header=TRUE)

LTB4R <- read.table("ASTHM_vs_Control_Beta2/14_all/LTB4R/log_chr14_LTB4R_ASTHM.assoc.logistic",
                    header=TRUE)

RORA <- read.table("ASTHM_vs_Control_Beta2/15_all/RORA/log_chr15_RORA_ASTHM.assoc.logistic",
                   header=TRUE)

DPEP1 <- read.table("ASTHM_vs_Control_Beta2/16_all/DPEP1/log_chr16_DPEP1_ASTHM.assoc.logistic",
                    header=TRUE)

DPEP2 <- read.table("ASTHM_vs_Control_Beta2/16_all/DPEP2/log_chr16_DPEP2_ASTHM.assoc.logistic",
                    header=TRUE)

DPEP3 <- read.table("ASTHM_vs_Control_Beta2/16_all/DPEP3/log_chr16_DPEP3_ASTHM.assoc.logistic",
                    header=TRUE)

ALOX12 <- read.table("ASTHM_vs_Control_Beta2/17_all/ALOX12/log_chr17_ALOX12_ASTHM.assoc.logistic",
                     header=TRUE)

ALOX15 <- read.table("ASTHM_vs_Control_Beta2/17_all/ALOX15/log_chr17_ALOX15_ASTHM.assoc.logistic",
                     header=TRUE)

ALOX15B <- read.table("ASTHM_vs_Control_Beta2/17_all/ALOX15B/log_chr17_ALOX15B_ASTHM.assoc.logistic",
                      header=TRUE)

FPR2 <- read.table("ASTHM_vs_Control_Beta2/19_all/FPR2/log_chr19_FPR2_ASTHM.assoc.logistic",
                   header=TRUE)

GPR32 <- read.table("ASTHM_vs_Control_Beta2/19_all/GPR32/log_chr19_GPR32_ASTHM.assoc.logistic",
                    header=TRUE)

GGT1 <- read.table("ASTHM_vs_Control_Beta2/22_all/GGT1/log_chr22_GGT1_ASTHM.assoc.logistic",
                   header=TRUE)

GGT2 <- read.table("ASTHM_vs_Control_Beta2/22_all/GGT2/log_chr22_GGT2_ASTHM.assoc.logistic",
                   header=TRUE)

CYSLTR1 <- read.table("ASTHM_vs_Control_Beta2/X_all/CYSLTR1/log_chrX_CYSLTR1_ASTHM.assoc.logistic",
                      header=TRUE)

GPR101 <- read.table("ASTHM_vs_Control_Beta2/X_all/GPR101/log_chrX_GPR101_ASTHM.assoc.logistic",
                     header=TRUE)

GPR173 <- read.table("ASTHM_vs_Control_Beta2/X_all/GPR173/log_chrX_GPR173_ASTHM.assoc.logistic",
                     header=TRUE)



# Frequency results files 
GSTM4_freq <- read.table("ASTHM_vs_Control_Beta2/1_all/GSTM4/log_chr1_GSTM4_ASTHM.frq.cc",
                         header=TRUE)

LGR6_freq <- read.table("ASTHM_vs_Control_Beta2/1_all/LGR6/log_chr1_LGR6_ASTHM.frq.cc",
                        header=TRUE)

mGST3_freq <- read.table("ASTHM_vs_Control_Beta2/1_all/mGST3/log_chr1_mGST3_ASTHM.frq.cc",
                         header=TRUE)

PTGS2_freq <- read.table("ASTHM_vs_Control_Beta2/1_all/PTGS2/log_chr1_PTGS2_ASTHM.frq.cc",
                         header=TRUE)

mGST2_freq <- read.table("ASTHM_vs_Control_Beta2/4_all/mGST2/log_chr4_mGST2_ASTHM.frq.cc",
                         header=TRUE)

LTC4S_freq <- read.table("ASTHM_vs_Control_Beta2/5_all/LTC4S/log_chr5_LTC4S_ASTHM.frq.cc",
                         header=TRUE)

AHR_freq  <- read.table("ASTHM_vs_Control_Beta2/7_all/AHR/log_chr7_AHR_ASTHM.frq.cc",
                        header=TRUE)

GPR37_freq  <- read.table("ASTHM_vs_Control_Beta2/7_all/GPR37/log_chr7_GPR37_ASTHM.frq.cc",
                          header=TRUE)

EPHX2_freq  <- read.table("ASTHM_vs_Control_Beta2/8_all/EPHX2/log_chr8_EPHX2_ASTHM.frq.cc",
                          header=TRUE)

ALOX5_freq  <- read.table("ASTHM_vs_Control_Beta2/10_all/ALOX5/log_chr10_ALOX5_ASTHM.frq.cc",
                          header=TRUE)

CMKLR1_freq <- read.table("ASTHM_vs_Control_Beta2/12_all/CMKLR1/log_chr12_CMKLR1_ASTHM.frq.cc",
                          header=TRUE)

LTA4H_freq <- read.table("ASTHM_vs_Control_Beta2/12_all/LTA4H/log_chr12_LTA4H_ASTHM.frq.cc",
                         header=TRUE)

CYSLTR2_freq <- read.table("ASTHM_vs_Control_Beta2/13_all/CYSLTR2/log_chr13_CYSLTR2_ASTHM.frq.cc",
                           header=TRUE)

GPR18_freq <- read.table("ASTHM_vs_Control_Beta2/13_all/GPR18/log_chr13_GPR18_ASTHM.frq.cc",
                         header=TRUE)

GPR99_freq <- read.table("ASTHM_vs_Control_Beta2/13_all/GPR99/log_chr13_GPR99_ASTHM.frq.cc",
                         header=TRUE)

LTB4R_freq <- read.table("ASTHM_vs_Control_Beta2/14_all/LTB4R/log_chr14_LTB4R_ASTHM.frq.cc",
                         header=TRUE)

RORA_freq <- read.table("ASTHM_vs_Control_Beta2/15_all/RORA/log_chr15_RORA_ASTHM.frq.cc",
                        header=TRUE)

DPEP1_freq <- read.table("ASTHM_vs_Control_Beta2/16_all/DPEP1/log_chr16_DPEP1_ASTHM.frq.cc",
                         header=TRUE)

DPEP2_freq <- read.table("ASTHM_vs_Control_Beta2/16_all/DPEP2/log_chr16_DPEP2_ASTHM.frq.cc",
                         header=TRUE)

DPEP3_freq <- read.table("ASTHM_vs_Control_Beta2/16_all/DPEP3/log_chr16_DPEP3_ASTHM.frq.cc",
                         header=TRUE)

ALOX12_freq <- read.table("ASTHM_vs_Control_Beta2/17_all/ALOX12/log_chr17_ALOX12_ASTHM.frq.cc",
                          header=TRUE)

ALOX15_freq <- read.table("ASTHM_vs_Control_Beta2/17_all/ALOX15/log_chr17_ALOX15_ASTHM.frq.cc",
                          header=TRUE)

ALOX15B_freq <- read.table("ASTHM_vs_Control_Beta2/17_all/ALOX15B/log_chr17_ALOX15B_ASTHM.frq.cc",
                           header=TRUE)

FPR2_frq <- read.table("ASTHM_vs_Control_Beta2/19_all/FPR2/log_chr19_FPR2_ASTHM.frq.cc",
                       header=TRUE)

GPR32_frq <- read.table("ASTHM_vs_Control_Beta2/19_all/GPR32/log_chr19_GPR32_ASTHM.frq.cc",
                        header=TRUE)


GGT1_frq <- read.table("ASTHM_vs_Control_Beta2/22_all/GGT1/log_chr22_GGT1_ASTHM.frq.cc",
                       header=TRUE)

GGT2_frq <- read.table("ASTHM_vs_Control_Beta2/22_all/GGT2/log_chr22_GGT2_ASTHM.frq.cc",
                       header=TRUE)

CYSLTR1_frq <- read.table("ASTHM_vs_Control_Beta2/X_all/CYSLTR1/log_chrX_CYSLTR1_ASTHM.frq.cc",
                          header=TRUE)

GPR101_frq <- read.table("ASTHM_vs_Control_Beta2/X_all/GPR101/log_chrX_GPR101_ASTHM.frq.cc",
                         header=TRUE)

GPR173_frq <- read.table("ASTHM_vs_Control_Beta2/X_all/GPR173/log_chrX_GPR173_ASTHM.frq.cc",
                         header=TRUE)



# Join association and frequency results 
GSTM4_all <- left_join(GSTM4, GSTM4_freq)

LGR6_all <- left_join(LGR6, LGR6_freq)

mGST3_all <- left_join(mGST3, mGST3_freq)

PTGS2_all <- left_join(PTGS2, PTGS2_freq)

mGST2_all <- left_join(mGST2, mGST2_freq)

LTC4S_all <- left_join(LTC4S, LTC4S_freq)

AHR_all <- left_join(AHR, AHR_freq)

GPR37_all <- left_join(GPR37, GPR37_freq)

EPHX2_all <- left_join(EPHX2, EPHX2_freq)

ALOX5_all <- left_join(ALOX5, ALOX5_freq)

CMKLR1_all <- left_join(CMKLR1, CMKLR1_freq)

LTA4H_all <- left_join(LTA4H, LTA4H_freq)

CYSLTR2_all <- left_join(CYSLTR2, CYSLTR2_freq)

GPR18_all <- left_join(GPR18, GPR18_freq)

GPR99_all <- left_join(GPR99, GPR99_freq)

LTB4R_all <- left_join(LTB4R, LTB4R_freq)

RORA_all <- left_join(RORA, RORA_freq)

DPEP1_all <- left_join(DPEP1, DPEP1_freq)

DPEP2_all <- left_join(DPEP2, DPEP2_freq)

DPEP3_all <- left_join(DPEP3, DPEP3_freq)

ALOX12_all <- left_join(ALOX12, ALOX12_freq)

ALOX15_all <- left_join(ALOX15, ALOX15_freq)

ALOX15B_all <- left_join(ALOX15B, ALOX15B_freq)

FPR2_all <- left_join(FPR2, FPR2_frq)

GPR32_all <- left_join(GPR32, GPR32_frq)

GGT1_all <- left_join(GGT1, GGT1_frq)

GGT2_all <- left_join(GGT2, GGT2_frq)

CYSLTR1_all <- left_join(CYSLTR1, CYSLTR1_frq)

GPR101_all <- left_join(GPR101, GPR101_frq)

GPR173_all <- left_join(GPR173, GPR173_frq)



# Delete rows with missing P values, add a column with the gene name  
GSTM4  <- GSTM4_all %>% 
  filter(is.na(P) == FALSE) %>% 
  mutate(Gene = "GSTM4") 

LGR6 <- LGR6_all %>% 
  filter(is.na(P) == FALSE) %>% 
  mutate(Gene = "LGR6") 

mGST3 <- mGST3_all %>% 
  filter(is.na(P) == FALSE) %>% 
  mutate(Gene = "mGST3") 

PTGS2 <- PTGS2_all %>% 
  filter(is.na(P) == FALSE) %>% 
  mutate(Gene = "PTGS2") 

mGST2 <- mGST2_all %>% 
  filter(is.na(P) == FALSE) %>% 
  mutate(Gene = "mGST2")  

LTC4S <- LTC4S_all %>% 
  filter(is.na(P) == FALSE) %>% 
  mutate(Gene = "LTC4S") 

AHR <- AHR_all %>% 
  filter(is.na(P) == FALSE) %>% 
  mutate(Gene = "AHR") 

GPR37 <- GPR37_all %>% 
  filter(is.na(P) == FALSE) %>% 
  mutate(Gene = "GPR37") 

EPHX2 <- EPHX2_all %>% 
  filter(is.na(P) == FALSE) %>% 
  mutate(Gene = "EPHX2") 

ALOX5 <- ALOX5_all %>% 
  filter(is.na(P) == FALSE) %>% 
  mutate(Gene = "ALOX5") 

CMKLR1 <- CMKLR1_all %>% 
  filter(is.na(P) == FALSE) %>% 
  mutate(Gene = "CMKLR1") 

LTA4H <- LTA4H_all %>% 
  filter(is.na(P) == FALSE) %>% 
  mutate(Gene = "LTA4H") 

CYSLTR2 <- CYSLTR2_all %>% 
  filter(is.na(P) == FALSE) %>% 
  mutate(Gene = "CYSLTR2") 

GPR18 <- GPR18_all %>% 
  filter(is.na(P) == FALSE) %>% 
  mutate(Gene = "GPR18")

GPR99 <- GPR99_all %>% 
  filter(is.na(P) == FALSE) %>% 
  mutate(Gene = "GPR99")

LTB4R <- LTB4R_all %>% 
  filter(is.na(P) == FALSE) %>% 
  mutate(Gene = "LTB4R")

RORA  <- RORA_all  %>% 
  filter(is.na(P) == FALSE) %>% 
  mutate(Gene = "RORA")

DPEP1  <- DPEP1_all  %>% 
  filter(is.na(P) == FALSE) %>% 
  mutate(Gene = "DPEP1")

DPEP2  <- DPEP2_all  %>% 
  filter(is.na(P) == FALSE)  %>% 
  mutate(Gene = "DPEP2")

DPEP3  <- DPEP3_all  %>% 
  filter(is.na(P) == FALSE) %>% 
  mutate(Gene = "DPEP3")

ALOX12  <- ALOX12_all  %>% 
  filter(is.na(P) == FALSE) %>% 
  mutate(Gene = "ALOX12")

ALOX15  <- ALOX15_all  %>% 
  filter(is.na(P) == FALSE) %>% 
  mutate(Gene = "ALOX15")

ALOX15B  <- ALOX15B_all  %>% 
  filter(is.na(P) == FALSE) %>% 
  mutate(Gene = "ALOX15B")

FPR2  <- FPR2_all  %>% 
  filter(is.na(P) == FALSE) %>% 
  mutate(Gene = "FPR2")

GPR32  <- GPR32_all  %>% 
  filter(is.na(P) == FALSE) %>% 
  mutate(Gene = "GPR32")

GGT1  <- GGT1_all  %>% 
  filter(is.na(P) == FALSE) %>% 
  mutate(Gene = "GGT1")

GGT2  <- GGT2_all  %>% 
  filter(is.na(P) == FALSE) %>% 
  mutate(Gene = "GGT2")

CYSLTR1  <- CYSLTR1_all  %>% 
  filter(is.na(P) == FALSE) %>% 
  mutate(Gene = "CYSLTR1")

GPR101 <- GPR101_all  %>% 
  filter(is.na(P) == FALSE) %>% 
  mutate(Gene = "GPR101")

GPR173 <- GPR173_all  %>% 
  filter(is.na(P) == FALSE) %>% 
  mutate(Gene = "GPR173")


# Combine all candidate genes into a single table 
Table_results  <- GSTM4

Table_results<-  bind_rows(Table_results, LGR6)

Table_results<- bind_rows(Table_results, mGST3)

Table_results<- bind_rows(Table_results, PTGS2)

Table_results<- bind_rows(Table_results, mGST2)

Table_results<- bind_rows(Table_results, LTC4S)

Table_results<- bind_rows(Table_results, AHR)

Table_results<- bind_rows(Table_results, GPR37)

Table_results<- bind_rows(Table_results, EPHX2)

Table_results <- bind_rows(Table_results, ALOX5)

Table_results<- bind_rows(Table_results, CMKLR1)

Table_results<- bind_rows(Table_results, LTA4H)

Table_results<- bind_rows(Table_results, CYSLTR2)

Table_results<- bind_rows(Table_results, GPR18)

Table_results<- bind_rows(Table_results, GPR99)

Table_results<- bind_rows(Table_results, LTB4R)

Table_results<- bind_rows(Table_results, RORA)

Table_results<- bind_rows(Table_results, DPEP1)

Table_results<- bind_rows(Table_results, DPEP2)

Table_results<-  bind_rows(Table_results, DPEP3)

Table_results<- bind_rows(Table_results, ALOX12)

Table_results<- bind_rows(Table_results, ALOX15)

Table_results<- bind_rows(Table_results, ALOX15B)

Table_results<- bind_rows(Table_results, FPR2)

Table_results<- bind_rows(Table_results, GPR32)

Table_results<- bind_rows(Table_results, GGT1)

Table_results<- bind_rows(Table_results, GGT2)

Table_results<- bind_rows(Table_results, CYSLTR1)

Table_results<- bind_rows(Table_results, GPR101)

Table_results <-  bind_rows(Table_results, GPR173)

# Remove duplicated SNPs
Table_results_uniq <- Table_results %>% 
  distinct(SNP, .keep_all= TRUE)

# Adjust raw p values with  False discovery rate
Table_results_uniq  <- Table_results_uniq  %>% 
  mutate(P_adjusted_BH = p.adjust(P, method = "BH")) 

write.table(Table_results_uniq, 
            file = paste("/Users/christianarenasdaza/Desktop/lab_ch2/output/6_QC_and_association_test/SPMsRG_logistic_add_results.tsv", sep = ""),
            sep = "\t",
            quote = FALSE,
            row.names = FALSE)



# Table of significant SNPS
Sig_SNP_FDR  <- Table_results_uniq  %>% 
  filter( P <= 1.401e-04) %>% 
  select(CHR, SNP, BP, Gene, A1, A2, MAF_A, P, P_adjusted_BH)


write.table(Sig_SNP_FDR, 
            file = paste("/Users/christianarenasdaza/Desktop/lab_ch2/output/6_QC_and_association_test/Sig_SNP_pass_FDR.tsv", sep = ""),
            sep = "\t",
            quote = FALSE,
            row.names = FALSE)



##### Data processing for manhattan 

# Group data by chromosome and position
Table_results_uniq <- Table_results_uniq %>% 
  group_by(CHR, BP)


# Create a column containing the relative location of the SNPs in the genome. 
Table_results_uniq <-  Table_results_uniq %>% 
  add_column(BPcum = 1)

first <- 3 

# calculate the values for bp cummulation
for (m in unique(Table_results_uniq$CHR)) {
  
  Table_results_uniq[Table_results_uniq$CHR == m, ]$BPcum <- seq(first, length.out = nrow(Table_results_uniq[Table_results_uniq$CHR == m, ]), by = 1)
  first <- max(Table_results_uniq$BPcum) + 10  
  
} 

# Assigned colors to genes 
chr_color <- function(genes) {
  colors <- c("darkgoldenrod", "steelblue2", "violet", "olivedrab3",
              "palevioletred4", "plum2", "bisque4", "saddlebrown", 
              "orangered", "mediumorchid2", "lightskyblue", "mediumpurple", 
              "coral", "darkgoldenrod4", "lightskyblue4", "darkslateblue",
              "darkorchid1", "firebrick2", "lightgoldenrod3", "violetred1", 
              "orange2", "burlywood", "lightpink4", "chocolate2",
              "cornflowerblue", "darkorchid1", "gray49", "chartreuse2", "cyan2", "deeppink2")  
  names(colors) <- c("GSTM4", "mGST3", "PTGS2", "LGR6", "mGST2", "LTC4S", "GPR37", "EPHX2", "ALOX5",
                     "CMKLR1", "CYSLTR2", "GPR99", "GPR18", "LTB4R", "DPEP3", "DPEP2", "DPEP1", "ALOX15",
                     "ALOX12", "ALOX15B", "GPR32", "FPR2", "GGT2", "GGT1", "CYSLTR1", "GPR101", "GPR173","AHR",
                     "LTA4H", "RORA")
  col_sel <- colors[names(colors) %in% genes]
  col_sel
}

# Create a data frame to hold the labels for the manahttan plot (where the gene is located)
label_pos <- Table_results_uniq %>% 
  # Group by gene name 
  # Group by gene  and calculate median between the maximum and miminimum BPcum. 
  group_by(Gene) %>% summarize(center=(max(BPcum) + min(BPcum) ) / 2)

# Established the colors for the different genes:
color_genes <- chr_color(unique(Table_results_uniq$Gene))

# Select the Y intercept for raw p-values based on the lowest p value and BP cum.
max_p_value <- -log10(min(Table_results_uniq$P))


##### Manhatan plot 
# 5.571602e-02 (benjamini & hochberg correction) = 1.401e-04 ( raw p_value)

# Creates the Manhattan plot:
manhatan_plot <- ggplot(Table_results_uniq, aes(x=BPcum, y=-log10(P))) +
  geom_point(aes(color=as.factor(Gene)), size=3) + 
  geom_point(data = Table_results_uniq, aes(color=as.factor(Gene)), alpha=1, shape=16 ,size=3) + 
  # Location of the labels for every gene:
  scale_x_continuous(label = label_pos$Gene, breaks= label_pos$center) +
  # Add the panels for each chromosome:
  facet_grid(.~ CHR, space = 'free_x', scales = 'free_x', switch = 'x') +
  scale_y_continuous(expand = c(0, 0)) +
  labs(x = "Chromosomes") +
  coord_cartesian(ylim = c(0, max_p_value+0.3)) + #Y limit
  # Manhattan threshold false discovery rate :
  geom_hline(yintercept = -log10(1.401e-04), color="blue", linetype="dashed") +
  # Label all significant SNPs: 
  geom_text_repel(data = Table_results_uniq[Table_results_uniq$P < (1.401e-04), ], position = "identity",
                  label = Table_results_uniq[Table_results_uniq$P < (1.401e-04), ]$SNP, size = 3, color="black") +
  geom_text_repel(data = Table_results_uniq[Table_results_uniq$P == (1.401e-04), ], position = "identity",
                  label = Table_results_uniq[Table_results_uniq$P == (1.401e-04), ]$SNP, size = 3, color="black",
                  vjust = 1) + # No significant but lies wihtin the threshold
  labs(caption = "Blue line = pvalue < 1.401e-04, which the same as adj pvalue (BH) < 0.05", color = "black", size = 3) +
  # Alteration to the themes 
  theme(
    axis.title = element_text(size = 15, colour = "black"),
    axis.text.x  =  element_text(size = 12, vjust = 0.5, angle =90, colour = "black"), 
    axis.text.y  = element_text(size = 12, hjust = 1, colour = "black"),
    axis.line = element_line(colour = 'black', size = 0.5), 
    axis.ticks = element_line(colour = "black", size = 0.5),  
    axis.ticks.length = unit(0.1, "cm"),
    legend.position="none",
    strip.text.x = element_text(size = 12, colour = "black"),
    strip.background = element_rect(fill="gray88"),
    panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_rect(fill = "ghostwhite"),
    strip.placement = 'outside',
    plot.caption = element_text(color = "Black", size = 10))


png(file = paste("/Users/christianarenasdaza/Desktop/lab_ch2/output/6_QC_and_association_test/Manhattan.png", sep = ""),
    width = 1200, height = 900)

print(manhatan_plot)

dev.off()


###### QQ plot

# QQ plot fucntion
gg_qqplot <- function(ps, ci = 0.95) {
  n  <- length(ps)
  df <- data.frame(
    observed = -log10(sort(ps)),
    expected = -log10(ppoints(n)),
    clower   = -log10(qbeta(p = (1 - ci) / 2, shape1 = 1:n, shape2 = n:1)),
    cupper   = -log10(qbeta(p = (1 + ci) / 2, shape1 = 1:n, shape2 = n:1))
  )
  log10Pe <- expression(paste("Expected -log"[10], plain((p-value))))
  log10Po <- expression(paste("Observed -log"[10], plain((p-value))))
  ggplot(df) +
    geom_ribbon(
      mapping = aes(x = expected, ymin = clower, ymax = cupper),
      alpha = 0.1
    ) +
    geom_point(aes(expected, observed), shape = 20, size = 3) +
    geom_abline(intercept = 0, slope = 1, alpha = 0.5, color = "red", size = 0.8) +
    # geom_line(aes(expected, cupper), linetype = 2, size = 0.5) +
    # geom_line(aes(expected, clower), linetype = 2, size = 0.5) +
    xlab(log10Pe) +
    ylab(log10Po)
}

# Calculation of  genomic inflation  rate
inflation <- function(ps) {
  chisq <- qchisq(1 - ps, 1)
  lambda <- median(chisq) / qchisq(0.5, 1)
  lambda
}

# Produce QQ plot
q_q_plot <- gg_qqplot(Table_results_uniq$P) +
  theme_bw(base_size = 24) +
  annotate(geom = "text", x = -Inf, y = Inf, hjust = -0.15, vjust = 1 + 0.15 * 3,
           label = sprintf("λ = %.2f", inflation(Table_results_uniq$P)), size = 8) +
  theme(axis.title = element_text(size = 25),
        axis.text.x  =  element_text(size = 25, colour = "black"),
        axis.text.y  = element_text(size = 25, hjust = 1, colour = "black"),
        legend.title = element_text(size = 25),
        legend.text  = element_text(size = 20),
        legend.key = element_rect(fill = "white"),
        legend.key.size = unit(1.3, "cm"), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        plot.margin=unit(c(1, 3, 1, 3), "cm")) 


png(file = paste("/Users/christianarenasdaza/Desktop/lab_ch2/output/6_QC_and_association_test/QQ_plot.png", sep = ""),
    width = 1200, height = 900)

print(q_q_plot)

dev.off()

