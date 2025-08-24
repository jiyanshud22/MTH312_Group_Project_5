library(dplyr)
library(ggplot2)

#------Loading Data-------#
data_mutations <- read.csv("../project5-jandj-project5_group19/Data/msk_chord_2024/data_mutations.txt",
                           header = TRUE, sep = "\t")
data_clinical_sample <- read.csv('../project5-jandj-project5_group19/Data/msk_chord_2024/data_clinical_sample.txt',
                                 header = TRUE, sep = "\t", skip = 4)
data_clinical_patient <- read.csv('../project5-jandj-project5_group19/Data/msk_chord_2024/data_clinical_patient.txt',
                                  header = TRUE, sep = "\t", skip = 4)
data_cna = read.csv("../project5-jandj-project5_group19/Data/msk_chord_2024/data_cna.txt",
                    header = TRUE, sep = "\t")


#--------Method 1(Basic)------------#
# Rename the Sample Identifier column to remove '#' and match mutation dataset
colnames(data_clinical_sample)[colnames(data_clinical_sample) == "SAMPLE_ID"] <- "Tumor_Sample_Barcode"

#Merged Data
merged_data = data_mutations %>% select(Hugo_Symbol, Tumor_Sample_Barcode)
merged_data = merge(merged_data, data_clinical_sample[,c("Tumor_Sample_Barcode", "CANCER_TYPE")],
                    by.x = "Tumor_Sample_Barcode", by.y ="Tumor_Sample_Barcode")

# Count mutations by cancer type
mutation_by_cancer <- merged_data %>%
  group_by(Hugo_Symbol, `CANCER_TYPE`) %>%
  summarise(count = n(), .groups = 'drop')
mutation_by_cancer

# Find mutations present in at least 2 different cancer types
mutation_summary <- mutation_by_cancer %>%
  group_by(Hugo_Symbol) %>%
  summarise(CANCER_TYPES = n_distinct(`CANCER_TYPE`), total_count = sum(count)) %>%
  filter(CANCER_TYPES >= 2) %>%
  arrange(desc(total_count))
mutation_summary

#--------Method 2(Basic)------------#

# Sum all columns except the first one
sums <- rowSums(abs(data_cna[, -1]), na.rm = TRUE)

# Add the sums as a new column to the DataFrame
data_cna$sums <- sums

#Sorted Mutation Counts(Absolute Values)
sorted_mut_counts = data_cna %>% select(Hugo_Symbol, sums) %>% filter(sums != 0) %>% arrange(desc(sums)) 
sorted_mut_counts

# Getting all the patients with non 0 values
symbol_wise_patients <- apply(data_cna[, -c(1, ncol(data_cna))], 1, 
                function(row) names(data_cna)[-c(1, ncol(data_cna))][row != 0])

# Assigning row names as names in the outer list
names(symbol_wise_patients) <- data_cna['Hugo_Symbol']$Hugo_Symbol
symbol_wise_patients

#Defining function to get cancer types from gene name
get_types = function(type)
{
  return(data_clinical_sample[data_clinical_sample$Tumor_Sample_Barcode %in% 
                       gsub("\\.", "-", symbol_wise_patients[[type]]), "CANCER_TYPE"])
}
#Testing for CDKN2A
counts_CDKN2A = get_types("CDKN2A")
table(counts_CDKN2A) #To confirm that this mutation is prevalent in 2 of the 5

#------Getting Baseline Demographics and Disease Summaries for Specifc Gene------#

demo_from_gene = function(selected_gene)
{
  # Extract all Tumor Sample Barcodes for the selected gene
  selected_gene_barcodes <- merged_data %>%
    filter(Hugo_Symbol == selected_gene) %>%
    pull(Tumor_Sample_Barcode)
  
  # Extract unique Patient IDs from Tumor Sample Barcodes
  patient_id <- unique(sub("^([A-Z]-\\d+).*", "\\1", selected_gene_barcodes))
  length(patient_id)  # Number of unique patients
  
  # Filter clinical patient data based on extracted patient IDs
  filtered_clinical_data <- data_clinical_patient %>%
    filter(PATIENT_ID %in% patient_id)
  
  # Select relevant demographic and disease characteristic columns
  selected_columns <- c("PATIENT_ID", "CURRENT_AGE_DEID", "GENDER", "ETHNICITY", 
                        "PRIOR_MED_TO_MSK", "OS_MONTHS", "OS_STATUS", "SMOKING_PREDICTIONS_3_CLASSES")
  
  # Filter the dataset for selected columns
  filtered_clinical_data_summary <- filtered_clinical_data %>%
    select(all_of(selected_columns))
  
  return(filtered_clinical_data_summary)
}

demo_TP53 = demo_from_gene("TP53")
demo_TP53

