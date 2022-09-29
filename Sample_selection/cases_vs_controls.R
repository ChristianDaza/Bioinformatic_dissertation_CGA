########Coding  controls and cases for PLINK###########

# This script code controls as 1 and cases as 2, to fit wiht PLINK formatting

##### Load library
library(tidyverse)

##### Load files 
controls <- read.table(file = paste("/Users/christianarenasdaza/Desktop/lab_ch2/output/4_Control_selection_ch/controls_demography_no_IC.tab", sep = ""),
                             header=TRUE, 
                             sep="\t")

cases <- read.table(file = paste("/Users/christianarenasdaza/Desktop/lab_ch2/input/4_Control_selection_ch/samples_AST_stric_sm.tab", sep = ""),
                    header=TRUE, 
                    sep="\t")

##### Data formatting

# take the id of cases and controls and turn them into their own data frame 
cases_id <-  as.data.frame(cases$id)
controls_id <- as.data.frame(controls$id)


# Re-name the id column 
names(cases_id)[1] <- "id"
names(controls_id)[1] <- "id"

# combine cases and control ids together
AST_IDs <- rbind(cases_id, controls_id)


# add the coding 
AST_IDs <- AST_IDs %>% mutate(codings = ifelse(id == cases$id, "2", "1"))


##### Save dataframes

# Save the IDs and codes as a table 
write.table(AST_IDs, 
            file = paste("/Users/christianarenasdaza/Desktop/lab_ch2/output/5_cases_vs_controls/ASTHM_join_id.tab", sep = ""),
            sep = "\t",
            quote = FALSE,
            row.names = FALSE) 


# Save  controls IDs  as a table 
write.table(controls_id, 
            file = paste("/Users/christianarenasdaza/Desktop/lab_ch2/output/5_cases_vs_controls/controls_ASTHM_id.tab", sep = ""),
            sep = "\t",
            quote = FALSE,
            row.names = FALSE, 
            col.names = FALSE) 

# Save cases IDs  as table 
write.table(cases_id, 
            file = paste("/Users/christianarenasdaza/Desktop/lab_ch2/output/5_cases_vs_controls/cases_ASTHM_id.tab", sep = ""),
            sep = "\t",
            quote = FALSE,
            row.names = FALSE,
            col.names = FALSE) 





