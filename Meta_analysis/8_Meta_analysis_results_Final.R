################### Meta-analysis results  analysis ############################

# This script is use to analyse the meta analysis results:
# Find significant SNPs with the same direction 
# Find validated SNPs
# Create a forest plot 

##### Load libraries 
library(tidyverse)
library(ggplot2)
library(ggrepel)
library(ggtext)

##### set directory
setwd("/Users/christianarenasdaza/Desktop/lab_ch2/input/8_Meta-analysis")

##### Read files 
Table_logistic_add_results <- read.table("/Users/christianarenasdaza/Desktop/lab_ch2/output/6_QC_and_association_test/SPMsRG_logistic_add_results.tsv",
                                         header=TRUE)

Meta_results <- read.table(file = paste("META_CGA_results.TBL", sep = ""),
                                header=TRUE, 
                                sep="\t")

sig_snps <- read.table(file = paste("/Users/christianarenasdaza/Desktop/lab_ch2/output/6_QC_and_association_test/Sig_SNP_pass_FDR.tsv", sep = ""),
                             header=TRUE, 
                             sep="\t")


##### Formating data 

# Extract SNPs that pass the FDAR threshold and have the same direction
Meta_results_FDR <- Meta_results  %>% 
  filter(P.value <= 1.401e-04) %>% 
  filter(Direction == "++" | Direction == "--") %>% 
  rename(SNP = MarkerName)
                         
# Merge meta significant SNPs with the Logistic regression results table                          
Meta_results_FDR_merged <-  Meta_results_FDR %>% 
  merge(Table_logistic_add_results, by = "SNP") %>% # all SNPs overlaped
  select(CHR, SNP, BP,Gene, Allele1, Allele2, Direction, Zscore,
         P.value, P, P_adjusted_BH, HetISq, HetChiSq, HetPVal) 


# Order  P values by descending order
Meta_results_FDR_merged  <- Meta_results_FDR_merged[order(Meta_results_FDR_merged$P.value), ]

# write  table of Meta analysis significant SNPs
write.table(Meta_results_FDR_merged, 
            file = paste("/Users/christianarenasdaza/Desktop/lab_ch2/output/8_Meta-analysis/Meta_sig_SNPs.tsv", sep = ""),
            sep = "\t",
            quote = FALSE,
            row.names = FALSE)

##### Forest plot preparations 

# A Function to assigned individual colorsfor candidate genes 
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

# Established the colors for candidate genes:
color_genes <- chr_color(unique(Meta_results_FDR_merged$Gene))

# Create labels based on genes and numbers to prevent  genes from overlapping 
Meta_results_FDR_merged  <- Meta_results_FDR_merged  %>% 
  mutate(SNP_separation =  1:nrow(Meta_results_FDR_merged)) %>% 
  mutate(label = paste(Gene, SNP_separation, sep = "_")) %>% 
  mutate(label = factor(label, levels = label))



# Repeat colors for same genes 
color_y_axis <- vector()
for (c in 1:nrow(Meta_results_FDR_merged)) { 
  color_add <- as.character(chr_color(Meta_results_FDR_merged$Gene[[c]]))
  color_y_axis <- append(color_y_axis, color_add)
}


##### Generate forest plot 
size_plot <- ggplot(Meta_results_FDR_merged, aes(x= Zscore, y=label)) +
  geom_point(aes(size=-log10(P.value), color=as.factor(Gene))) +
  scale_color_manual(values = color_genes)  + 
  xlim(c(min(Meta_results_FDR_merged$Zscore)-0.05, max(Meta_results_FDR_merged$Zscore)+0.05)) + 
  scale_y_discrete(limits = rev, labels = rev(Meta_results_FDR_merged$Gene)) +
  geom_text(data = Meta_results_FDR_merged[!Meta_results_FDR_merged$SNP %in% sig_snps$SNP, ], 
            aes(max(Zscore)+0.05, label, label = SNP), color = "black",
            hjust = 0.7, vjust= -0.8,  size = 4, position =position_dodge(0.6)) +
  geom_text(data = Meta_results_FDR_merged[Meta_results_FDR_merged$SNP %in% sig_snps$SNP, ], 
            aes(max(Zscore)+0.05, label, label = SNP), color = "red",
            hjust = 0.5, vjust= -0.8, size = 4, position=position_dodge(0.5)) +
  geom_vline(xintercept = 0, color="black", linetype="dashed", size = 1) +
  guides(colour = "none") +
  labs(x = "Zscore", y = "Candidate Genes", size = "-log10(P)") +
  theme(
    axis.title = element_text(size = 20, colour = "black"),
    axis.text.x  =  element_text(size = 15, vjust = 0.5, colour = "black"), 
    axis.text.y  = element_text(size = 15, hjust = 1, colour = rev(color_y_axis)),
    axis.line = element_line(colour = 'black', size = 0.5), 
    axis.ticks = element_line(colour = "black", size = 0.5),  
    axis.ticks.length = unit(0.1, "cm"),
    legend.position="bottom",
    legend.text = element_text(size = 15, vjust = 0.5, colour = "black"), 
    legend.title = element_text(size = 15, vjust = 0.5, colour = "black"),
    legend.key = element_rect(fill = "white"),
    panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_line(color = "gray89", size = 0.5, linetype = 2),
    panel.background = element_blank())

png(file = paste("/Users/christianarenasdaza/Desktop/lab_ch2/output/8_Meta-analysis/Forest_plotsize.png", sep = ""),
    width = 1200, height = 900)

print(size_plot)

dev.off()



