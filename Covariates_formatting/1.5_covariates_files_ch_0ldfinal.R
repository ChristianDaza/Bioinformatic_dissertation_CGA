########Covariates formatting for PLINK1.9 Logistic regression test###########

# This script is used to format the covraiates I chose from the provided 
# UKBiobank data into a PLINK format, which is needed for runing any of the 
# association test in PLIN1.9

##### Load pakages  

library(stringr)
library(dplyr) # use sample_frac; for taking a percentage of the whole data set.
library(tidyverse)

##### Set the directory (share folder in my case)
setwd("/Users/christianarenasdaza/Desktop/lab_ch2/input/1.5_covariates_file_ch")


##### Read files 

subset_center <- read.table(file = paste("ukb43914.tab", sep = ""),
                            header=TRUE, 
                            sep="\t")


samples_covariates <- read.table(file = paste("column_subset_cov_AST_ch_sh.tab", sep = ""),
                                 header=TRUE, 
                                 sep="\t")


samples_pca <- read.table(file = paste("column_subset_pca.tab", sep = ""),
                          header=TRUE, 
                          sep="\t")


##### Format covariates, including their intances (repeated measurements)

# Merge columns for easier formatting 
covariate_table <- merge(samples_covariates, subset_center)

# Get variables I need 
covariate_table$FID <- covariate_table$f.eid

covariate_table$AGE <- 2022 - covariate_table$f.34.0.0 # Convert date of birth into actual age 

covariate_table$BMI <- rowMeans(covariate_table[, c(15:18)], na.rm = TRUE) # Takes mean of BMI values 

covariate_table$BMI <- round(covariate_table$BMI, digits = 2)

covariate_table$SEX <- covariate_table$f.31.0.0

covariate_table$Moderate_excercise <- paste(covariate_table$f.884.0.0, 
                                                        covariate_table$f.884.1.0,
                                                        covariate_table$f.884.2.0,
                                                        covariate_table$f.884.3.0)

covariate_table$Smoking <- paste(covariate_table$f.20116.0.0, covariate_table$f.20116.1.0,  
                                    covariate_table$f.20116.2.0, covariate_table$f.20116.3.0)

covariate_table$CENTRE <- paste(covariate_table$f.54.0.0, covariate_table$f.54.1.0, 
                                   covariate_table$f.54.2.0, covariate_table$f.54.3.0)

covariate_table$ARRAY <- covariate_table$f.22000.0.0

# Create a dataframe with just my covariates 
final_covariate_table<- covariate_table[, c(24:31)] 




######## String formatting 
# convert NAs into spaces 
final_covariate_table$Smoking <- gsub("NA", "",  final_covariate_table$Smoking )
final_covariate_table$CENTRE <- gsub("NA", "", final_covariate_table$CENTRE)
final_covariate_table$Moderate_excercise <- gsub("NA", "", final_covariate_table$Moderate_excercise)

# Replace all NA into -9 (PLINK code to signify missing data)
final_covariate_table$BMI <- replace_na(final_covariate_table$BMI, -9)
final_covariate_table$AGE <- replace_na(final_covariate_table$AGE, -9)
final_covariate_table$SEX <- replace_na(final_covariate_table$SEX, -9)

# Format CENTRE column
final_covariate_table$CENTRE <- gsub(" ", "", final_covariate_table$CENTRE) # Remove white spaces
final_covariate_table$CENTRE <- strtrim(final_covariate_table$CENTRE, 5) # Take the first code (first 5 values) 
final_covariate_table$CENTRE <- gsub("^$", "-9", final_covariate_table$CENTRE) # Fill empty cells 

# Format SMOOKING column
final_covariate_table$Smoking <- gsub(" ", "", final_covariate_table$Smoking) # Remove white spaces
final_covariate_table$Smoking <- gsub("3", "-9", final_covariate_table$Smoking) # replace missing data
final_covariate_table$Smoking <- gsub("^$", "-9", final_covariate_table$Smoking) # replace missing data
final_covariate_table$Smoking <- str_extract(final_covariate_table$Smoking, "^[0,1,2]") # select the first code 
final_covariate_table$Smoking <- replace_na(final_covariate_table$Smoking, "-9") # replace NAs

# Format BMI column
final_covariate_table$BMI <- gsub(" ", "", final_covariate_table$BMI) # Remove white spaces
final_covariate_table$BMI <- gsub("^$", "-9", final_covariate_table$BMI) # replace NAs

# Format Exercise column
final_covariate_table$Moderate_excercise <- gsub(" ", "", final_covariate_table$Moderate_excercise) # Remove white spaces
final_covariate_table$Moderate_excercise <- gsub("-3", "-9", final_covariate_table$Moderate_excercise) # replace missing data
final_covariate_table$Moderate_excercise <- gsub("-1", "-9", final_covariate_table$Moderate_excercise) # replace missing data
final_covariate_table$Moderate_excercise  <- str_extract(final_covariate_table$Moderate_excercise , "(?![-]?[9])([+]?[0-7])") # select the first code
final_covariate_table$Moderate_excercise  <- replace_na(final_covariate_table$Moderate_excercise, "-9") # replace NAs

# Format ARRAY column 
# two arrays, code by positive or negative numbers 
final_covariate_table$ARRAY[final_covariate_table$ARRAY >= 0] <- 1 # change all positive codes 
final_covariate_table$ARRAY[final_covariate_table$ARRAY < 0] <- 0 # change all negative codes 
final_covariate_table$ARRAY[is.na(final_covariate_table$ARRAY)]<--9 # Replace missing data 

# Delete completely BLANK rows 
final_covariate_table <- final_covariate_table[-162909, ]



##### Formatting principal component values 

#  Extract and format the first 10 PCA
samples_pca_10 <- samples_pca[, c(1:11)]
colnames(samples_pca_10) <- c("FID", "PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8", "PC9", "PC10")
samples_pca_10[is.na(samples_pca_10)] <- -9
samples_pca_10 <- round(samples_pca_10, digits = 2)

# Get the final covariates table
covariates_final <- merge(final_covariate_table, samples_pca_10) 



##### Other formatting (additional columns/deleting rows)
# Add an extra column to satisfy PLINK formatting ( requires a family ID column ; IDD)
covariates_final$IID <- covariates_final$FID # Duplicate ID column to be the family ID column
covariates_final<- covariates_final[, c(1,ncol(covariates_final), 2:(ncol(covariates_final)-1))] # relocate column

# Delete rows wiht missing array data 
covariates_final_no_NA <- covariates_final[covariates_final$ARRAY != -9, ]



##### Create dummy tables 

#PLINK requires dummy column for all categorical variables with > 2 levels, for 
# each  level of their category to be in their own column.
# These are coded 0 if not value or 1 if it had a value.

# Create a data frame wiht the same order of ID column as previous
dummy <- data.frame(FID = covariates_final_no_NA$FID,
                    IID = covariates_final_no_NA$IID)


# Funtion to generate dummy columns
for (i in 3:length(covariates_final_no_NA)) {
  
  # Just creates dummy columns for those variables with levels between 2 and 35 (Age will be take as continous).
  
  if (length(unlist(unique(covariates_final_no_NA[i]))) > 2 & length(unlist(unique(covariates_final_no_NA[i]))) < 35) {
    
    # Creates a list with the unique values for each variable. 
    
    list_unique <- unlist(unique(covariates_final_no_NA[i]))
    list_unique <- list_unique[list_unique != -9] # We remove all the -9 levels since they are not levels but missing values
    
    # Here we start in 2 since the dummy columns will be K-1, meaning that the first level is used as reference: 
    
    for (j in 2:length(list_unique)) {
      
      dummy$column <- unlist(covariates_final_no_NA[i]) # Copy the column from the main covariates file 
      
      # We need to add a crazy value because some variables already have the zero and one as part of their levels, which makes
      # imposible to transform into a zeros and ones binary system. -9 has to be exclused since it will remained in the dummy 
      # columns. 
      
      dummy$column[dummy$column != -9] <- (as.numeric(dummy$column[dummy$column != -9]) + 50000) 
      dummy$column <- gsub((as.numeric(list_unique[[j]])+50000), 1, dummy$column) # Replace the value for 1. 
      dummy$column[dummy$column != 1 & dummy$column != -9] <- 0 # Everything different to 1 and -9 go to zero. 
      
      # Replace the column name with the format of plink.
      
      names(dummy)[names(dummy) == "column"] <- paste(colnames(covariates_final_no_NA[i]), list_unique[[j]], sep = "_")
      
    }
  }
}



# Merge dummy columns with the continuous and catgorical variables 

dummy_covariates <- merge(covariates_final_no_NA, dummy)


# Delete all the categorical variable original column thta were transform into 
# dummy columns 
dummy_covariates_final <- dummy_covariates[, -c(6:8)]


# Create table without the sec variable 
dummy_covariates_no_sex <- dummy_covariates_final[, -5]


##### Save all the files to you desire directory

# Save the table with the demography information: 

write.table(covariates_final, 
            file = paste(output, "covariates_file_full_AST_ch.tab", sep = ""),
            sep = "\t",
            quote = FALSE,
            row.names = FALSE)  

write.table(covariates_final_no_NA, 
            file = paste(output, "covariates_file_AST_ch.tab", sep = ""),
            sep = "\t",
            quote = FALSE,
            row.names = FALSE)  

write.table(dummy_covariates_final, 
            file = paste(output, "covariates_dummy_AST_ch.tab", sep = ""),
            sep = "\t",
            quote = FALSE,
            row.names = FALSE)  

write.table(dummy_covariates_no_sex, 
            file = paste(output, "covariates_dummy_no_sex_ASTHM.tab", sep = ""),
            sep = "\t",
            quote = FALSE,
            row.names = FALSE) 
 






