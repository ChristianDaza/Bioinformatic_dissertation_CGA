############GABRIEL Consortium Validation preparation####################

# This script is used for:
 # Join  GABRIEL Consortium Data to their corresponding position in the 37 ref genome
 # Delete duplicate  SNPs and those with missing P random model values
 # Subset GABRIEL consortium data for SNPs within the adjusted candidate genes
 # Adjust for multiple testing correction
 # Create Manhattan plot



##### Load libraries
library("tidyverse")
library("data.table")
library("ggrepel")
library("ggplot2")


##### Set directory
setwd("/Users/christianarenasdaza/Desktop/lab_ch2/input/10_validation")

##### Load data
validation <- read.table(file = paste("GABRIEL.tsv", sep = ""),
                        header=TRUE, 
                        sep="\t")

Table_results_uniq <- read.table(file = paste("/Users/christianarenasdaza/Desktop/lab_ch2/output/6_QC_and_association_test/SPMsRG_logistic_add_results.tsv", sep = ""),
                                 header=TRUE,
                                 sep="\t")

validation_37positions <- read.table("GABRIEL_og_SNP_PS_hg18_hg19.tsv",
                                     fill = TRUE,
                                     header= T)

genes <- read.table(file = paste("genes_table_f_subseingt_valSNPs.tab", sep = ""),
                         header=TRUE, 
                         sep="\t")


##### Data formatting SNP subsetting

# Prepare validation data set for SNP subsetting using the values from the 37
# reference genome vcf file

validation_final <-  validation %>% 
  left_join(validation_37positions) %>% 
  relocate(c(Chr_37, rs_37, position_37), .after = position) %>% 
  # dismiss these original Gabriel columns
  select(!c(Chr, rs, position)) %>%
    rename(Chr = Chr_37, rs = rs_37, Start = position_37) %>% 
    mutate(End = Start ) %>%
    relocate(End, .after =  Start) %>% 
    filter(is.na(Chr) == FALSE & is.na(rs) == FALSE & is.na(Start) == FALSE &
             is.na(End) == FALSE & is.na(P_ran) == FALSE) %>% 
    distinct(rs, .keep_all= TRUE)
    

# formatting the candidate gene table for SNP sub setting
genes <- genes %>% 
  rename(Chr= Chromosome, Start = Start.500Kb.b., End = End.500kb.b.)


##### Validation SNP subsetting for those found in candidate genes

# convert data frame  to data.table type for SNP sub setting
validation_dt <- setDT(validation_final )

genes_dt <- setDT(genes)


# Convert Chromosome columns values into characters
validation_dt$Chr <- as.character(validation_dt$Chr)



# set keys  to  subset  for SNPs and specify the candidate gene table
setkey(genes_dt, "Chr", "Start", "End")


# Map SNPs in The Gabriel data set that are located within the 
# candidate genes
val_snp_GCA <- foverlaps(validation_dt, genes_dt)

# Subset for SNPs in cadidate genes 
clean_val_results_dt <- val_snp_GCA[is.na(Gene) == FALSE, ]

# convert  data.table object into a data frame
clean_val_results_df  <- setDF(clean_val_results_dt)


# Delete unecessary columns and  match column name to those to association results
val_results_final <- clean_val_results_df %>% 
  rename(PS = i.Start) %>% 
  select(!(c(i.End, Start, End))) %>% 
  relocate(rs, .after = Gene)  %>% 
  select(Chr, Gene, rs, PS,  P_ran) %>% 
  rename( CHR = Chr, SNP = rs, BP = PS, P = P_ran)

# Find the number of SNPs overlapping the subsetted GABRIEL data and my logistic 
# regression results
table(val_results_final$SNP %in%Table_results_uniq$SNP) 


# turn chromsome column back to integer
val_results_final$CHR <-  as.integer(val_results_final$CHR)

write.table(val_results_final, 
            file = paste("/Users/christianarenasdaza/Desktop/lab_ch2/output/10_validation/Validation_SNP_Subset.tsv", sep = ""),
            sep = "\t",
            quote = FALSE,
            row.names = FALSE)


# Adjust raw p values with  False discovery rate
val_results_final <- val_results_final  %>% 
  mutate(P_adjusted_BH = p.adjust(P, method = "BH")) 


##### Data formatting for  Manhattan

# Group data  by chromosome and position
Table_results_uniq <- val_results_final %>% 
  group_by(CHR, BP)


# Create a column containing the relative location of the SNPs in the genome. 

# create column for BP cummulation
val_results_final <-  val_results_final %>% 
  add_column(BPcum = 1)

first <- 3 

# calculate the values for bp cummulation
for (m in unique(val_results_final$CHR)) {
  
  val_results_final[val_results_final$CHR == m, ]$BPcum <- seq(first, length.out = nrow(val_results_final[val_results_final$CHR == m, ]), by = 1)
  first <- max(val_results_final$BPcum) + 10  
  
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
label_pos <- val_results_final %>% 
  
  # Group by gene name 
  # Create a new column for the  median between the maximum and miminimum BPcum. 
  group_by(Gene) %>% summarize(center=(max(BPcum) + min(BPcum) ) / 2)

# Established the colors for the different genes:
color_genes <- chr_color(unique(val_results_final$Gene))

# Select the Y intercept for raw p-values based on the lowest p value and BP cum.
max_p_value <- -log10(min(val_results_final$P))


##### Creates the Manhattan plot:

manhatan_plot <- ggplot(val_results_final, aes(x=BPcum, y=-log10(P))) +
  geom_point(aes(color=as.factor(Gene)), size=3) + # All points
  geom_point(data = val_results_final, aes(color=as.factor(Gene)), alpha=1, shape=16 ,size=3) + # Only SNPs associated with genes.
  scale_color_manual(values = color_genes) +
  # Location of the labels for every gene:
  scale_x_continuous(label = label_pos$Gene, breaks= label_pos$center) +
  # Add the panels for each chromosome:
  facet_grid(.~ CHR, space = 'free_x', scales = 'free_x', switch = 'x') +
  scale_y_continuous(expand = c(0, 0)) +
  labs(x = "Chromosomes") +
  coord_cartesian(ylim = c(0, max_p_value+0.3)) + #Y limit
  # Horizontal line representing significant SNPs :
  geom_hline(yintercept = -log10(1.401e-04), color="red", linetype="dashed") +
  # Add SNP name to significant SNPs adter multiple test correction: 
  #geom_text_repel( data=subset(assoc_add_results_uniq, is_annotate =="FDR_pass"), aes(label=SNP), size=4) +
  geom_text_repel(data = val_results_final[val_results_final$P_adjusted_BH < (0.05), ], position = "identity",
                  label = val_results_final[val_results_final$P_adjusted_BH < (0.05), ]$SNP, size = 3, color="black",
                  vjust = 0.6, hjust = 0.1) +
  # Captions explaining  the  two Y intercept and their meaning.
  labs(caption = "Red line = Adj pvalue (BH) < 0.05", color = "black", size = 3) +
  # Modification to the theme 
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

# Save manhattan plot
png(file = paste("/Users/christianarenasdaza/Desktop/lab_ch2/output/10_validation/MANHATTAN_validation_GABRIEL_VAL.png", sep = ""),
    width = 1200, height = 900)

print(manhatan_plot)

dev.off()

