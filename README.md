# Analysis scripts for the manuscript: "Oral Microbial Determinants of Saliva and Serum Lipopolysaccharide Activity" by Muhammed Manzoor et al.2025
#The data in this study are not publicly available due to the sensitive nature of the data derived from human subjects, including personal information. However, researchers who wish to use our data may do so by requesting from the SECRETO Oral Study Consortia Data Access Committee (https://ega-archive.org/dacs/EGAC00001003449).
#Load packages
#BiocManager::install("microbiome/mia", version="devel")
packages <- c("ggplot2", "biomformat", "ggthemes", "phyloseq", "vegan", "microbiomeMarker",  "microbiome", "tidyverse",  "reshape2", "ComplexHeatmap", "maptree", "patchwork", "dplyr" "miaViz", "file2meco")
#Phenotype data is loaded from the included R object
#Load metadata and biom file
biom <- import_biom("cuatroc.biom", parseFunction=parse_taxonomy_greengenes, parallel=TRUE)
metadata <- import_qiime_sample_data("metadata.txt") #metadata are not publicly available due to the sensitive nature of the data derived from human subjects, including personal information. 
#Construct the primary phyloseq object.
phyloseq <- merge_phyloseq(biom, metadata)
#ADD TREE information to the phyloseq object.
random_tree = rtree(ntaxa(phyloseq), rooted=TRUE, tip.label=taxa_names(phyloseq))
phyloseq <- merge_phyloseq(biom, metadata, random_tree)
#make tse from phyloseq
tse<- makeTreeSummarizedExperimentFromPhyloseq(phyloseq)

#remove participants who have used antibiotics preceeding 3 months of saliva sampling
tse <- tse[ , colData(tse)$Antibiotics %in% FALSE]

# Calculate tertile cut-off values
salivaLPS_cutoffs <- quantile(data$salivaLPS, probs = c(0, 1/3, 2/3, 1), na.rm = TRUE)
SERUM_LPS_cutoffs <- quantile(data$SerumLPS, probs = c(0, 1/3, 2/3, 1), na.rm = TRUE
# List of categorical variables
categorical_vars <- c("Gender", "Education", "smoking", "Diabetics", "Alcoholheavy",
                      "HYPERTENSION", "Obesity", "stroke", "antihypertensive", "statin", 
                      "Antibiotic use",  "regular_visit",
                      "CARIES", "Periodontalstatus")

# Function to calculate count (n) and percentage (%) within each tertile group
categorical_summary <- function(data, group_var) {
  data %>%
    select(all_of(group_var), all_of(categorical_vars)) %>%
    pivot_longer(cols = all_of(categorical_vars), names_to = "Variable", values_to = "Category") %>%
    group_by(across(all_of(group_var)), Variable, Category) %>%
    summarise(n = n(), .groups = "drop") %>%
    mutate(total_in_group = sum(n),  # Total count within each tertile group
           percentage = round(100 * n / total_in_group, 1)) %>%
    select(all_of(group_var), Variable, Category, n, percentage) %>%
    arrange(Variable, Category)
}

# Function to check normality and perform appropriate tests
get_p_value <- function(data, variable, group_var) {
  if (variable %in% continuous_vars) {
    # Check normality for continuous variables
    p_norm <- shapiro.test(data[[variable]])$p.value
    
    # Select test based on normality
    if (p_norm > 0.05) {
      # Normally distributed → Use ANOVA
      test_result <- aov(as.formula(paste(variable, "~", group_var)), data = data)
      p_value <- summary(test_result)[[1]][["Pr(>F)"]][1]
      test_type <- "ANOVA"
    } else {
      # Not normally distributed → Use Kruskal-Wallis test
      test_result <- kruskal.test(as.formula(paste(variable, "~", group_var)), data = data)
      p_value <- test_result$p.value
      test_type <- "Kruskal-Wallis"
    }
  } else if (variable %in% categorical_vars) {
    # For categorical variables, create a contingency table
    contingency_table <- table(data[[variable]], data[[group_var]])
    
    # Choose Chi-square or Fisher's exact test
    if (all(contingency_table >= 5)) {
      test_result <- chisq.test(contingency_table)
      p_value <- test_result$p.value
      test_type <- "Chi-square"
    } else {
      test_result <- fisher.test(contingency_table)
      p_value <- test_result$p.value
      test_type <- "Fisher’s Exact"
    }
  } else {
    return(NULL)
  }
  
  return(data.frame(Variable = variable, Group = group_var, Test = test_type, P_Value = p_value))
}

# Function to compute p-values for all variables
compute_p_values <- function(data, group_var) {
  results <- bind_rows(
    lapply(c(continuous_vars, categorical_vars), function(var) get_p_value(data, var, group_var))
  )
  return(results)
}

# Compute p-values for both tertiles
p_values_saliva <- compute_p_values(data, "salivaLPS_tertile")
p_values_serum <- compute_p_values(data, "SERUM_LPS_tertile")

# Save results as CSV files
write.csv(p_values_saliva, "p_values_salivaLPS.csv", row.names = FALSE)
write.csv(p_values_serum, "p_values_SERUM_LPS.csv", row.names = FALSE)

# Print confirmation message
print("P-values calculated and saved as 'p_values_salivaLPS.csv' and 'p_values_SERUM_LPS.csv'")


