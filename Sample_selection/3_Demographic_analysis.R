###############################Asthma cases selection###########################

# This script is used to explore the demographic of the UKBiobank data 
# for asthma case selection, as well as extracting them for PLINK.




##### Load pakages 
library(ggplot2)
library(stringr)
library(cowplot)
library(htmlwidgets)
library(tidyverse)
library(dplyr)

##### set directory

# Set directory to be in the source file location:
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

##### Load files 

# This file had already been subsetted for Asthma patients using AWK
samples_subset <- read.table(file = paste("cases_AST_ch.tab", sep = ""),
                             header=TRUE, 
                             sep="\t")

patients_witdraw <-read.table(file = paste("samples_witdraw.tab", sep = ""),
                              header = TRUE,
                              col.names = "ID",
                              sep = "\t")


##### Data manipulation

# Get the columns that I care:
samples_subset$id <- samples_subset$f.eid
samples_subset$age <- 2022 - samples_subset$f.34.0.0
samples_subset$BMI <- rowMeans(samples_subset[, c(155:158)], na.rm = TRUE) # Takes mean of BMI values 
samples_subset$BMI <- round(samples_subset$BMI, digits = 2)
samples_subset$Sex<- samples_subset$f.31.0.0
samples_subset$Smoking <- paste(samples_subset$f.20116.0.0, samples_subset$f.20116.1.0,  
                                samples_subset$f.20116.2.0, samples_subset$f.20116.3.0)

samples_subset$ethnicity <- samples_subset$f.21000.0.0
samples_subset$ethnicity2 <- samples_subset$f.21000.1.0
samples_subset$ethnicity3 <- samples_subset$f.21000.2.0
samples_subset$CH3_genotype <- samples_subset$f.22101.0.0

# extract the columns that I care into a new variable 
samples_subset_demography<- samples_subset[, c(407:415)]

# Delete all the NA from the different columns:
samples_subset_demography$Smoking <- gsub("NA", "",  samples_subset_demography$Smoking)


# Formatting smooking column
# -3 = missing data in this column of UKBiobank
samples_subset_demography$Smoking <- gsub(" ", "", samples_subset_demography$Smoking)  # delete blank spaces 
samples_subset_demography$Smoking <- str_extract(samples_subset_demography$Smoking, "^[0,1,2]") # select the first code
samples_subset_demography$Smoking <- gsub("^$", "-3", samples_subset_demography$Smoking) # replace missing data 
samples_subset_demography$Smoking<- replace_na(samples_subset_demography$Smoking, "-3") # replace missing data 



##### Demographic plots 

# SEX plot

# Calulate the perctages of males and females 
sex <- data.frame(column = c("sex", "sex"), 
                  sex = c("Male", "Female"),
                  proportions = c((nrow(samples_subset_demography[samples_subset_demography$Sex == 1, ])/nrow(samples_subset_demography))*100,
                                  (nrow(samples_subset_demography[samples_subset_demography$Sex == 0, ])/nrow(samples_subset_demography))*100))

sex$proportions <- round(sex$proportions, 0)

sex_plot <- ggplot(data = sex, aes(x = column, y = proportions, fill = sex)) +
  geom_bar(stat = "identity", colour = "black", size = 2) +
  scale_fill_manual(values=c('brown1','deepskyblue3')) + 
  geom_text(data = sex, aes(x = column, y = proportions, label = paste(proportions,"%",sep="")), 
            size = 5, position = "stack", colour = "black", hjust = c(1.5, 1.5)) +
  labs(x = "Sex", y = "Proportions") + 
  coord_flip() + 
  scale_y_continuous(expand = c(0, 1), limits = c(0, 110)) + 
  theme(axis.title = element_text(size = 20),
        legend.position = "top",
        legend.title = element_blank(),
        legend.key.size = unit(1, "cm"), 
        legend.text  = element_text(size = 20),
        axis.text.x  =  element_text(size = 20, hjust = 0.5, colour = "black"), # Put color to the labels
        axis.text.y  = element_text(size = 20, hjust = 0.5, colour = "black"), # Put color to the labels
        axis.line = element_line(colour = 'black', size = 1), # Color and thickness of axis
        axis.ticks = element_line(colour = "black", size = 1), # Color and thickness of every axis sep. 
        panel.background = element_rect(fill = "white"),
        axis.ticks.length = unit(0.4, "cm"))  # Size of every axis sep.  

sex_plot


# Age plot

age_plot<- ggplot(samples_subset_demography, aes(x = age, fill = as.factor(Sex))) + 
  geom_density(alpha = 0.5) + 
  labs(x = "Age", y = "Density") + 
  scale_fill_manual(values=c('brown1', 'deepskyblue3'), labels =  c("Female", "Male")) + 
  theme(axis.title = element_text(size = 20),
        legend.position = "top",
        legend.title = element_blank(),
        legend.key.size = unit(1, "cm"), 
        legend.text  = element_text(size = 20),
        axis.text.x  =  element_text(size = 20, hjust = 0.5, colour = "black"), # Put color to the labels
        axis.text.y  = element_text(size = 20, hjust = 0.5, colour = "black"), # Put color to the labels
        axis.line = element_line(colour = 'black', size = 1), # Color and thickness of axis
        axis.ticks = element_line(colour = "black", size = 1), # Color and thickness of every axis sep. 
        panel.background = element_rect(fill = "white"),
        axis.ticks.length = unit(0.4, "cm"))  # Size of every axis sep.  

age_plot

# Ethnicity: 

# save  all possitble ethnicity base on UKBiobank only resources for this  variables
etnias = c("White", "British", "Irish", "Any other white background", "White and Black Caribbean", 
           "White and Black African", "White and Asian", "Any other mixed background",  
           "Indian", "Pakistani", "Bangladeshi", "Any other Asian background", "Caribbean", 
           "African", "Any other Black background", "Chinese", "Other ethnic group", 
           "Prefer not to answer", "Do not know", "NA")

# Search for all ethnicities based on their code 
ethnicity <- data.frame(skin =  factor(c(rep("White", 4), rep("Mixed", 4), rep("Asian or Asian British", 4), 
                                         rep("Black or Black British", 3), "Chinese", "Other", rep("NA", 3)), 
                                       levels = c("White", "Mixed", "Asian or Asian British", "Black or Black British", 
                                                  "Chinese", "Other", "NA")),
                        etnia = factor(etnias, levels = etnias),
                        numbers = c(nrow(samples_subset_demography[samples_subset_demography$ethnicity == 1, ]),
                                    nrow(samples_subset_demography[samples_subset_demography$ethnicity == 1001, ]),
                                    nrow(samples_subset_demography[samples_subset_demography$ethnicity == 1002, ]),
                                    nrow(samples_subset_demography[samples_subset_demography$ethnicity == 1003, ]),
                                    nrow(samples_subset_demography[samples_subset_demography$ethnicity == 2001, ]),
                                    nrow(samples_subset_demography[samples_subset_demography$ethnicity == 2002, ]),
                                    nrow(samples_subset_demography[samples_subset_demography$ethnicity == 2003, ]),
                                    nrow(samples_subset_demography[samples_subset_demography$ethnicity == 2004, ]),
                                    nrow(samples_subset_demography[samples_subset_demography$ethnicity == 3001, ]),
                                    nrow(samples_subset_demography[samples_subset_demography$ethnicity == 3002, ]),
                                    nrow(samples_subset_demography[samples_subset_demography$ethnicity == 3003, ]),
                                    nrow(samples_subset_demography[samples_subset_demography$ethnicity == 3004, ]),
                                    nrow(samples_subset_demography[samples_subset_demography$ethnicity == 4001, ]),
                                    nrow(samples_subset_demography[samples_subset_demography$ethnicity == 4002, ]),
                                    nrow(samples_subset_demography[samples_subset_demography$ethnicity == 4003, ]),
                                    nrow(samples_subset_demography[samples_subset_demography$ethnicity == 5, ]),
                                    nrow(samples_subset_demography[samples_subset_demography$ethnicity == 6, ]),
                                    nrow(samples_subset_demography[samples_subset_demography$ethnicity == -3, ]),
                                    nrow(samples_subset_demography[samples_subset_demography$ethnicity == -1, ]),
                                    nrow(samples_subset_demography[samples_subset_demography$ethnicity == "NA", ])))

# Delete no represented ethnicities. 
ethnicity <- ethnicity[ethnicity$numbers > 0, ] 

# Ethnicity plot
normal_ethnicity <- ggplot(data = ethnicity, aes(x = etnia, y = numbers, fill = skin)) +
  geom_bar(stat = "identity",  position = position_dodge(), colour = "black", size = 2) +
  scale_fill_manual(values=c("goldenrod1", "goldenrod4", "gold1", "gold4", "darkgoldenrod2", "gray", "black")) + 
   scale_y_continuous(limits = c(0, 40000)) +
  geom_text(position = position_dodge(w=1), vjust = -0.7,
            aes(label = numbers), size = 7) +
  labs(x = "Ethnic Background", y = "Number of samples") + 
  theme(axis.title = element_text(size = 20),
        legend.position = c(0.9, 0.7),
        legend.title = element_blank(),
        legend.key.size = unit(1.3, "cm"), 
        legend.text  = element_text(size = 20),
        axis.title.x = element_blank(),
        axis.text.x  =  element_text(size = 20, hjust = 1, angle = 45, colour = "black"), # Put color to the labels
        axis.text.y  = element_text(size = 20, hjust = 0.5, colour = "black"), # Put color to the labels
        axis.line = element_line(colour = 'black', size = 1.5), # Color and thickness of axis
        axis.ticks = element_line(colour = "black", size = 1.5), # Color and thickness of every axis sep. 
        panel.background = element_rect(fill = "white"),
        axis.ticks.length = unit(0.4, "cm"))  # Size of every axis sep. 

normal_ethnicity 

# Logarithmic plot ethnicity

log_ethnicity <- ggplot(data = ethnicity, aes(x = etnia, y = numbers, fill = skin)) +
  geom_bar(stat = "identity",  position = position_dodge(), colour = "black", size = 2) +
  scale_fill_manual(values=c("goldenrod1", "goldenrod4", "gold1", "gold4", "darkgoldenrod2", "gray", "black")) + 
  geom_text(position = position_dodge(w=1), vjust = -0.7,
            aes(label = numbers), size = 7) +
  labs(x = "Ethnic Background", y = "Number of samples") + 
  scale_y_log10() +
  theme(axis.title = element_text(size = 20),
        legend.title = element_blank(),
        legend.key.size = unit(1.3, "cm"), 
        legend.text  = element_text(size = 20),
        axis.title.x = element_blank(),
        axis.text.x  =  element_text(size = 20, hjust = 1, angle = 45, colour = "black"), # Put color to the labels
        axis.text.y  = element_text(size = 20, hjust = 0.5, colour = "black"), # Put color to the labels
        axis.line = element_line(colour = 'black', size = 1.5), # Color and thickness of axis
        axis.ticks = element_line(colour = "black", size = 1.5), # Color and thickness of every axis sep. 
        panel.background = element_rect(fill = "white"),
        axis.ticks.length = unit(0.4, "cm"))  # Size of every axis sep.

log_ethnicity

# Smooking plot

addiction <- data.frame(substance = factor(c(rep("Smoking", 4))), 
                        answer = factor(c("Never", "Previous", "Current", "Prefer not to answer"), 
                                        levels = c("Never", "Previous", "Current", "Prefer not to answer")), 
                        proportions = c((nrow(samples_subset_demography[samples_subset_demography$Smoking == 0, ])/nrow(samples_subset_demography))*100,
                                        (nrow(samples_subset_demography[samples_subset_demography$Smoking == 1, ])/nrow(samples_subset_demography))*100,
                                        (nrow(samples_subset_demography[samples_subset_demography$Smoking == 2, ])/nrow(samples_subset_demography))*100,
                                        (nrow(samples_subset_demography[samples_subset_demography$Smoking == -3, ])/nrow(samples_subset_demography))*100))

addiction$proportions <- round(addiction$proportions, 0)

addiction_plot <- ggplot(data = addiction, aes(x = answer, y = proportions, fill = substance)) +
  geom_bar(stat = "identity", position=position_dodge(), colour = "black", size = 2) +
  scale_fill_manual(values=c('lightblue','darkgray')) + 
  geom_text(data = addiction, aes(x = answer, y = proportions, label = paste(proportions,"%",sep="")), 
            size = 5, position = position_dodge(0.9), colour = "black", vjust = -0.5) +
  labs(x = "Cuestionary Answer", y = "Proportions") + 
  scale_y_continuous(expand = c(0, 1), limits = c(0, 110)) + 
  theme(axis.title = element_text(size = 20),
        legend.position = "top",
        legend.title = element_blank(),
        legend.key.size = unit(1, "cm"), 
        legend.text  = element_text(size = 20),
        axis.text.x  =  element_text(size = 20, hjust = 0.5, colour = "black"), # Put color to the labels
        axis.text.y  = element_text(size = 20, hjust = 0.5, colour = "black"), # Put color to the labels
        axis.line = element_line(colour = 'black', size = 1), # Color and thickness of axis
        axis.ticks = element_line(colour = "black", size = 1), # Color and thickness of every axis sep. 
        panel.background = element_rect(fill = "white"),
        axis.ticks.length = unit(0.4, "cm"))  # Size of every axis sep.

addiction_plot

##### Saving files to your desire directory

# Save cases samples
write.table(samples_subset_demography, 
            file = paste("/Users/christianarenasdaza/Desktop/lab_ch2/output/3_cases_demographic_ch/samples_demography_AST.tab", sep = ""),
            sep = "\t",
            quote = FALSE,
            row.names = FALSE)  

# Save the plots

png(file = paste("/Users/christianarenasdaza/Desktop/lab_ch2/output/3_cases_demographic_ch/sex_AST.png", sep = ""), width = 900, height = 600)

sex_plot

dev.off()  

png(file = paste("/Users/christianarenasdaza/Desktop/lab_ch2/output/3_cases_demographic_ch/age_AST.png", sep = ""), width = 900, height = 600)

age_plot

dev.off()  

png(file = paste("/Users/christianarenasdaza/Desktop/lab_ch2/output/3_cases_demographic_ch/normal_ethnicity_AST.png", sep = ""), width = 1300, height = 1100)

normal_ethnicity

dev.off()  

png(file = paste("/Users/christianarenasdaza/Desktop/lab_ch2/output/3_cases_demographic_ch/log_ethnicity_AST.png", sep = ""), width = 1300, height = 1100)

log_ethnicity

dev.off()  

png(file = paste("/Users/christianarenasdaza/Desktop/lab_ch2/output/3_cases_demographic_ch/Smookings_AST.png", sep = ""), width = 900, height = 600)

addiction_plot

dev.off()  

#####  Subsetting fro caucassan population  and age 52 or older 
cg_samples <- subset(samples_subset_demography, samples_subset_demography$age >= 52 & 
                       (samples_subset_demography$ethnicity == "1001" | 
                          samples_subset_demography$ethnicity2 == "1001" | 
                          samples_subset_demography$ethnicity3 == "1001" |
                          samples_subset_demography$ethnicity == "1002" | 
                          samples_subset_demography$ethnicity2 == "1002" | 
                          samples_subset_demography$ethnicity3 == "1002" |
                          samples_subset_demography$ethnicity == "1003" | 
                          samples_subset_demography$ethnicity2 == "1003" | 
                          samples_subset_demography$ethnicity3 == "1003" |
                          samples_subset_demography$ethnicity == "1" | 
                          samples_subset_demography$ethnicity2 == "1" | 
                          samples_subset_demography$ethnicity3 == "1") & 
                       samples_subset_demography$CH3_genotype  == "#ukbgene")






# Eliminate withdrawal samples 
cg_samples <- cg_samples[!(cg_samples$id %in% patients_witdraw$ID), ]



# save  final asthma case  file
write.table(cg_samples, 
            file = paste("/Users/christianarenasdaza/Desktop/lab_ch2/output/3_cases_demographic_ch/samples_AST_stric_sm.tab", sep = ""),
            sep = "\t",
            quote = FALSE,
            row.names = FALSE)


##### Demographics ASTHMA Caucasian population 

# Sex plot  

sex_cg <- data.frame(column = c("sex", "sex"), 
                  sex = c("Male", "Female"),
                  proportions = c((nrow(cg_samples[cg_samples$Sex == 1, ])/nrow(cg_samples))*100,
                                  (nrow(cg_samples[cg_samples$Sex == 0, ])/nrow(cg_samples))*100))

sex_cg$proportions <- round(sex_cg$proportions, 0)

sex_plot <- ggplot(data = sex_cg, aes(x = column, y = proportions, fill = sex)) +
  geom_bar(stat = "identity", colour = "black", size = 2) +
  scale_fill_manual(values=c('brown1','deepskyblue3')) + 
  geom_text(data = sex_cg, aes(x = column, y = proportions, label = paste(proportions,"%",sep="")), 
            size = 5, position = "stack", colour = "black", hjust = c(1.5, 1.5)) +
  labs(x = "Sex", y = "Proportions") + 
  coord_flip() + 
  scale_y_continuous(expand = c(0, 1), limits = c(0, 110)) + 
  theme(axis.title = element_text(size = 20),
        legend.position = "top",
        legend.title = element_blank(),
        legend.key.size = unit(1, "cm"), 
        legend.text  = element_text(size = 20),
        axis.text.x  =  element_text(size = 20, hjust = 0.5, colour = "black"), # Put color to the labels
        axis.text.y  = element_text(size = 20, hjust = 0.5, colour = "black"), # Put color to the labels
        axis.line = element_line(colour = 'black', size = 1), # Color and thickness of axis
        axis.ticks = element_line(colour = "black", size = 1), # Color and thickness of every axis sep. 
        panel.background = element_rect(fill = "white"),
        axis.ticks.length = unit(0.4, "cm"))  # Size of every axis sep.  

sex_plot

png(file = paste("/Users/christianarenasdaza/Desktop/lab_ch2/output/3_cases_demographic_ch/sex_CG.png", sep = ""), width = 900, height = 600)

sex_plot

dev.off()  

# Age plot

age_plot<- ggplot(cg_samples, aes(x = age, fill = as.factor(Sex))) + 
  geom_density(alpha = 0.5) + 
  labs(x = "Age", y = "Density") + 
  scale_fill_manual(values=c('brown1', 'deepskyblue3'), labels =  c("Female", "Male")) + 
  theme(axis.title = element_text(size = 20),
        legend.position = "top",
        legend.title = element_blank(),
        legend.key.size = unit(1, "cm"), 
        legend.text  = element_text(size = 20),
        axis.text.x  =  element_text(size = 20, hjust = 0.5, colour = "black"), # Put color to the labels
        axis.text.y  = element_text(size = 20, hjust = 0.5, colour = "black"), # Put color to the labels
        axis.line = element_line(colour = 'black', size = 1), # Color and thickness of axis
        axis.ticks = element_line(colour = "black", size = 1), # Color and thickness of every axis sep. 
        panel.background = element_rect(fill = "white"),
        axis.ticks.length = unit(0.4, "cm"))  # Size of every axis sep.  

age_plot

png(file = paste("/Users/christianarenasdaza/Desktop/lab_ch2/output/3_cases_demographic_ch/age_CG.png", sep = ""), width = 900, height = 600)

age_plot

dev.off()  
