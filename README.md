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

#make tse from phyloseq
tse<- makeTreeSummarizedExperimentFromPhyloseq(phyloseq)

# 1. Check the current minimum sample depth
min_reads_before <- min(sample_sums(secreto_phyloseq1))
cat("Minimum reads before rarefaction:", min_reads_before, "\n")

# 2. Perform rarefaction to the minimum sample depth
set.seed(123)  # Set seed for reproducibility
secreto_phyloseq_rarefied <- rarefy_even_depth(secreto_phyloseq1, rngseed = 123, sample.size = min_reads_before, replace = FALSE)


taxa_counts <- apply(tax_table(secreto_phyloseq_rarefied), 2, function(x) length(unique(na.omit(x))))
print(taxa_counts)


# 3. Check read statistics after rarefaction
otu_counts_rarefied <- otu_table(secreto_phyloseq_rarefied)

sample_read_stats_rarefied <- data.frame(
  Min_Reads = min(rowSums(otu_counts_rarefied)),
  Max_Reads = max(rowSums(otu_counts_rarefied)),
  Mean_Reads = mean(rowSums(otu_counts_rarefied)),
  Median_Reads = median(rowSums(otu_counts_rarefied)),
  SD_Reads = sd(rowSums(otu_counts_rarefied)),
  Total_Reads = sum(otu_counts_rarefied)
)

print(sample_read_stats_rarefied)


taxa_counts <- apply(tax_table(secreto_phyloseq_rarefied), 2, function(x) length(unique(na.omit(x))))
print(taxa_counts)


# Define thresholds

prevalence_threshold <- 0.10  # Keep taxa present in at least 10% of samples
detection_threshold <- 0.001 # Keep taxa with relative abundance > 0.01%

# Compute prevalence (fraction of samples where each taxa is present)
prev <- apply(otu_table(secreto_phyloseq_rarefied), 1, function(x) sum(x > 0) / length(x))

# Filter taxa based on 10% prevalence
secreto_phyloseq_filtered <- prune_taxa(prev > prevalence_threshold, secreto_phyloseq_rarefied)

# Convert to relative abundance
secreto_phyloseq_filtered_rel <- transform_sample_counts(secreto_phyloseq_filtered, function(x) x / sum(x))

# Apply detection threshold
secreto_phyloseq_filtered <- prune_taxa(
  apply(otu_table(secreto_phyloseq_filtered_rel), 1, max) > detection_threshold, 
  secreto_phyloseq_filtered
)

# Check the number of taxa after filtering
filtered_taxa_counts <- apply(tax_table(secreto_phyloseq_filtered), 2, function(x) length(unique(na.omit(x))))
print(filtered_taxa_counts)

# Compute alpha diversity metrics
alpha_div <- estimate_richness(pseq, measures = c( "Shannon", "InvSimpson"))


# Add sample metadata
sample_data_df <- data.frame(sample_data(pseq))
alpha_div <- cbind(sample_data_df, alpha_div)

# View first few rows
head(alpha_div)


library(picante)  # Required for Faith PD

alpha_div <- alpha_div[, !duplicated(colnames(alpha_div))]
colnames(alpha_div)  # Verify that duplicates are gone

alpha_metrics <- c("Observed", "Shannon")  # Add more metrics if needed

for (metric in alpha_metrics) {
  p1 <- ggscatter(alpha_div, x = "salivaLPS", y = metric, 
                  add = "reg.line", conf.int = TRUE, 
                  cor.coef = TRUE, cor.method = "spearman",
                  xlab = "Saliva LPS", ylab = metric)
  
  p2 <- ggscatter(alpha_div, x = "SERUM_LPS", y = metric, 
                  add = "reg.line", conf.int = TRUE, 
                  cor.coef = TRUE, cor.method = "spearman",
                  xlab = "Serum LPS", ylab = metric)

  print(p1)
  print(p2)
}
# Define the alpha diversity metrics
alpha_metrics <- c("Shannon", "InvSimpson")

# Remove outliers based on IQR for SERUM_LPS only
remove_outliers <- function(data, var_name) {
    Q1 <- quantile(data[[var_name]], 0.25)
    Q3 <- quantile(data[[var_name]], 0.75)
    IQR <- Q3 - Q1
    data_clean <- data[data[[var_name]] >= (Q1 - 1.5 * IQR) & data[[var_name]] <= (Q3 + 1.5 * IQR), ]
    return(data_clean)
}

# Define the alpha diversity metrics
alpha_metrics <- c("Shannon", "InvSimpson")

# Remove outliers based on IQR for SERUM_LPS only
remove_outliers <- function(data, var_name) {
    Q1 <- quantile(data[[var_name]], 0.25)
    Q3 <- quantile(data[[var_name]], 0.75)
    IQR <- Q3 - Q1
    data_clean <- data[data[[var_name]] >= (Q1 - 1.5 * IQR) & data[[var_name]] <= (Q3 + 1.5 * IQR), ]
    return(data_clean)
}

# Clean the data by removing outliers from SERUM_LPS only
alpha_div_clean <- alpha_div
alpha_div_clean <- remove_outliers(alpha_div_clean, "SERUM_LPS")

# Generate plots for each alpha diversity metric
for (metric in alpha_metrics) {
    # Plot for Saliva LPS (without removing outliers)
    p1 <- ggscatter(alpha_div, x = "salivaLPS", y = metric, 
                    add = "reg.line", conf.int = TRUE, 
                    cor.coef = FALSE, cor.method = "spearman",  # Remove correlation coefficient and p-values
                    xlab = "Saliva LPS", ylab = metric)
    
    # Plot for Serum LPS (with the outlier removed)
    p2 <- ggscatter(alpha_div_clean, x = "SERUM_LPS", y = metric, 
                    add = "reg.line", conf.int = TRUE, 
                    cor.coef = FALSE, cor.method = "spearman",  # Remove correlation coefficient and p-values
                    xlab = "Serum LPS", ylab = metric)
    
    # Print the plots
    print(p1)
    print(p2)
    
    # Save the plots to files
    ggsave(paste0("SalivaLPS_vs_", metric, "_no_r2_pvalues.png"), plot = p1)
    ggsave(paste0("SerumLPS_vs_", metric, "_clean_no_r2_pvalues.png"), plot = p2)
    
    # Calculate correlation coefficient and p-value separately for Serum LPS plot
    cor_serum <- cor.test(alpha_div_clean$SERUM_LPS, alpha_div_clean[[metric]], method = "spearman")
    cor_saliva <- cor.test(alpha_div$salivaLPS, alpha_div[[metric]], method = "spearman")
    
    # Print the correlation coefficient and p-value for both plots
    cat("\nCorrelation between Serum LPS and", metric, ": \n")
    cat("Spearman's rho:", cor_serum$estimate, "\n")
    cat("p-value:", cor_serum$p.value, "\n")
    
    cat("\nCorrelation between Saliva LPS and", metric, ": \n")
    cat("Spearman's rho:", cor_saliva$estimate, "\n")
    cat("p-value:", cor_saliva$p.value, "\n")
}

# Step 1: Extract the taxonomy table from your phyloseq object
tax_table_data <- tax_table(pseq)

# Step 2: Extract the unique phylum names
unique_phyla <- unique(tax_table_data[, "Phylum"])

# Step 3: Print the unique phylum names
print(unique_phyla)


# Step 1: Extract the taxonomy table from your phyloseq object
tax_table_data <- tax_table(pseq)

# Step 2: Check the column names to ensure that 'Phylum' is accessible
colnames(tax_table_data)

# Step 3: Check the first few rows to ensure correct structure
head(tax_table_data)

# Step 4: Filter for Gram-negative bacteria (Proteobacteria, Bacteroidetes, Fusobacteria, Spirochaetes, Synergistetes)
# Accessing 'Phylum' from the taxonomy table and filtering accordingly
gram_negative_taxa <- tax_table_data[tax_table_data[, "Phylum"] %in% c("Proteobacteria", "Bacteroidetes", "Fusobacteria", "Spirochaetes", "Synergistetes"), ]

# Step 5: Filter the phyloseq object based on these Gram-negative phyla
gram_negative_pseq <- prune_taxa(rownames(gram_negative_taxa), pseq)

# Step 6: Check the new phyloseq object
gram_negative_pseq
# Step 1: Calculate alpha diversity
alpha_div <- estimate_richness(gram_negative_pseq, measures = c("Shannon", "InvSimpson"))

# Step 2: Check column names in alpha_div and sample_data
colnames(alpha_div)  # Check the column names in alpha_div
colnames(sample_data(gram_negative_pseq))  # Check column names in sample_data(gram_negative_pseq)

# Step 3: Make sure that the row names in alpha_div correspond to sample IDs in sample_data
alpha_div$sample_id <- rownames(alpha_div)  # Add sample_id to alpha_div

# Check column names again to ensure matching
colnames(alpha_div)  # Make sure the sample_id is added
colnames(sample_data(gram_negative_pseq))  # Check for the matching sample ID column

# Step 4: Merge data by the correct column (assuming "SampleID" in sample_data)
alpha_div <- merge(alpha_div, sample_data(gram_negative_pseq), by.x = "sample_id", by.y = "row.names", all.x = TRUE)

# Step 5: Perform correlation analysis for alpha diversity metrics
alpha_metrics <- c("Shannon", "InvSimpson")


# Remove outliers based on IQR for SERUM_LPS only (you can adjust this based on your criteria)
remove_outlier <- function(data, var_name) {
    Q1 <- quantile(data[[var_name]], 0.25)
    Q3 <- quantile(data[[var_name]], 0.75)
    IQR <- Q3 - Q1
    data_clean <- data[data[[var_name]] >= (Q1 - 1.5 * IQR) & data[[var_name]] <= (Q3 + 1.5 * IQR), ]
    return(data_clean)
}

# Clean the data by removing outliers from SERUM_LPS only
alpha_div_clean <- remove_outlier(alpha_div, "SERUM_LPS")

# Step 5: Perform correlation analysis for alpha diversity metrics with SERUM_LPS
correlation_results_serum <- data.frame(Metric = character(), Spearman_rho = numeric(), P_value = numeric())

for (metric in alpha_metrics) {
  # Perform Spearman correlation for serum_LPS
  correlation_test_serum <- cor.test(alpha_div_clean[[metric]], alpha_div_clean$SERUM_LPS, method = "spearman")
  
  # Store correlation results for serum_LPS
  correlation_results_serum <- rbind(correlation_results_serum, data.frame(
    Metric = metric,
    Spearman_rho = correlation_test_serum$estimate,  # Spearman's rho
    P_value = correlation_test_serum$p.value
  ))
  
  # Plot the correlation for serum_LPS (remove R² and p-value from plot)
  p_serum <- ggscatter(alpha_div_clean, x = "SERUM_LPS", y = metric, 
                       add = "reg.line", conf.int = TRUE, 
                       cor.coef = FALSE, cor.method = "spearman",  # Remove correlation coefficient and p-values
                       xlab = "Serum LPS", ylab = metric)
  
  # Save the plot to your directory
  ggsave(paste0("SerumLPS_vs_", metric, "_clean_no_r2_pvalues.png"), plot = p_serum)
  
  # Print the plot
  print(p_serum)
}

# Step 6: Display the correlation results for serum_LPS
print(correlation_results_serum)

# Boxplot for Shannon by Saliva LPS Group
ggplot(alpha_div, aes(x = salivaLPS_Group, y = Shannon, fill = salivaLPS_Group)) +
  geom_boxplot() +
  scale_fill_manual(values = c("A_Low" = "#009E73" , "B_Medium" = "#56B4E9", "C_High" =  "#F0E442")) + 
  labs(title = "Shannon by Saliva LPS Group", 
       x = "Saliva LPS Group", y = "Shannon") +
  theme_minimal()

# Boxplot for InvSimpson by Saliva LPS Group
ggplot(alpha_div, aes(x = salivaLPS_Group, y = InvSimpson, fill = salivaLPS_Group)) +
  geom_boxplot() +
  scale_fill_manual(values = c("A_Low" = "#009E73" , "B_Medium" = "#56B4E9", "C_High" =  "#F0E442")) + 
  labs(title = "InvSimpson by Saliva LPS Group", 
       x = "Saliva LPS Group", y = "InvSimpson") +
  theme_minimal()



# Step 1: Check that "salivaLPS_Group" is a factor and alpha_div has the alpha diversity metrics
alpha_div$salivaLPS_Group <- as.factor(alpha_div$salivaLPS_Group)  # Ensure categorical variable

# Step 2: Perform Kruskal-Wallis test for Shannon and InvSimpson across Saliva LPS Groups
kruskal_results <- data.frame(Metric = character(), 
                              Kruskal_Wallis_Chi2 = numeric(), 
                              P_value = numeric())

for (metric in alpha_metrics) {
  # Perform Kruskal-Wallis test
  kruskal_test <- kruskal.test(alpha_div[[metric]] ~ alpha_div$salivaLPS_Group)
  
  # Store the results
  kruskal_results <- rbind(kruskal_results, data.frame(
    Metric = metric,
    Kruskal_Wallis_Chi2 = kruskal_test$statistic,
    P_value = kruskal_test$p.value
  ))
  
  # Step 3: Plot the results for each metric (Shannon and InvSimpson)
  p <- ggboxplot(alpha_div, x = "salivaLPS_Group", y = metric, 
                 color = "salivaLPS_Group", palette = "jco", 
                 ylab = metric, xlab = "Saliva LPS Group", 
                 add = "jitter", title = paste("Distribution of", metric, "by Saliva LPS Group"))
  
  print(p)
}

# Step 4: Display Kruskal-Wallis results
print(kruskal_results)







# Boxplot for Shannon by Serum LPS Group
ggplot(alpha_div, aes(x = SERUM_LPS_Group, y = Shannon, fill = SERUM_LPS_Group)) +
  geom_boxplot() +
  scale_fill_manual(values = c("A_Low" = "#009E73" , "B_Medium" = "#56B4E9", "C_High" =  "#F0E442")) + 
  labs(title = "Shannon by Serum LPS Group", 
       x = "Serum LPS Group", y = "Shannon") +
  theme_minimal()

# Boxplot for InvSimpson by Serum LPS Group
ggplot(alpha_div, aes(x = SERUM_LPS_Group, y = InvSimpson, fill = SERUM_LPS_Group)) +
  geom_boxplot() +
  scale_fill_manual(values = c("A_Low" = "#009E73" , "B_Medium" = "#56B4E9", "C_High" =  "#F0E442")) + 
  labs(title = "InvSimpson by Serum LPS Group", 
       x = "Serum LPS Group", y = "InvSimpson") +
  theme_minimal()



# Step 1: Check that "SERUM_LPS_Group" is a factor and alpha_div has the alpha diversity metrics
alpha_div$SERUM_LPS_Group <- as.factor(alpha_div$SERUM_LPS_Group)  # Ensure categorical variable

# Step 2: Perform Kruskal-Wallis test for Shannon and InvSimpson across Serum LPS Groups
kruskal_results <- data.frame(Metric = character(), 
                              Kruskal_Wallis_Chi2 = numeric(), 
                              P_value = numeric())

for (metric in alpha_metrics) {
  # Perform Kruskal-Wallis test
  kruskal_test <- kruskal.test(alpha_div[[metric]] ~ alpha_div$SERUM_LPS_Group)
  
  # Store the results
  kruskal_results <- rbind(kruskal_results, data.frame(
    Metric = metric,
    Kruskal_Wallis_Chi2 = kruskal_test$statistic,
    P_value = kruskal_test$p.value
  ))
  
  # Step 3: Plot the results for each metric (Shannon and InvSimpson)
  p <- ggboxplot(alpha_div, x = "SERUM_LPS_Group", y = metric, 
                 color = "SERUM_LPS_Group", palette = "jco", 
                 ylab = metric, xlab = "Saliva LPS Group", 
                 add = "jitter", title = paste("Distribution of", metric, "by Serum LPS Group"))
  
  print(p)
}

# Step 4: Display Kruskal-Wallis results
print(kruskal_results)



