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


# Bray-Curtis distance calculation
# Compute Bray-Curtis dissimilarity matrix
# Step 1: Compute Bray-Curtis distance matrix
bray_dist <- distance(gram_negative_pseq, method = "bray")

# Step 2: Perform PCoA based on Bray-Curtis distance
pcoa <- ordinate(gram_negative_pseq, method = "PCoA", distance = bray_dist)

# Step 3: Create PCoA plot, coloring by salivaLPS_Group
pcoa_plot <- plot_ordination(gram_negative_pseq, pcoa, color = "salivaLPS_Group") +
  geom_point(size = 3) +  # Points for the PCoA plot
  stat_ellipse(aes(color = salivaLPS_Group), level = 0.95) +  # Add ellipses (95% confidence level)
  scale_color_manual(values = c("A_Low" = "#009E73", 
                                "B_Medium" = "#56B4E9", 
                                "C_High" = "#F0E442")) +  # Custom color palette
  labs(title = "Beta Diversity (PCoA) based on Bray-Curtis Distance") + 
  theme_minimal()  # Clean theme

# Display the plot
pcoa_plot




# Step 3: Create PCoA plot, coloring by SERUM_LPS_Group
pcoa_plot <- plot_ordination(gram_negative_pseq, pcoa, color = "SERUM_LPS_Group") +
  geom_point(size = 3) +  # Points for the PCoA plot
  stat_ellipse(aes(color = SERUM_LPS_Group), level = 0.95) +  # Add ellipses (95% confidence level)
  scale_color_manual(values = c("A_Low" = "#009E73", 
                                "B_Medium" = "#56B4E9", 
                                "C_High" = "#F0E442")) +  # Custom color palette
  labs(title = "Beta Diversity (PCoA) based on Bray-Curtis Distance") + 
  theme_minimal()  # Clean theme

  # Convert sample_data to a pure data frame
sample_data_df <- as.data.frame(as.matrix(sample_data(gram_negative_pseq)))

# Check if the variables are in the correct format
str(sample_data_df)

# Ensure the grouping variable is a factor
sample_data_df$SERUM_LPS_Group   <- factor(sample_data_df$SERUM_LPS_Group  )

# Run PERMANOVA
permanova_results <- adonis2(bray_dist ~ SERUM_LPS_Group  , 
                             data = sample_data_df, 
                             permutations = 999)

print(permanova_results)


permanova_results3 <- adonis2(bray_dist ~ salivaLPS_Group  , 
                             data = sample_data_df, 
                             permutations = 999)

print(permanova_results3)

#plot for abundance
# Step 1: Convert counts to relative abundance
pseq_rel <- transform_sample_counts(gram_negative_pseq, function(x) x / sum(x))

# Step 2: Collapse by Phylum level
pseq_phylum <- tax_glom(pseq_rel, taxrank = "Phylum")

# Step 3: Convert to a data frame
otu_data <- psmelt(pseq_phylum)

# Step 4: Calculate total abundance and select top 10 phyla
top_phyla <- otu_data %>%
  group_by(Phylum) %>%
  summarise(Total_Abundance = sum(Abundance), .groups = 'drop') %>%
  top_n(10, Total_Abundance) %>%
  pull(Phylum)

# Group the remaining phyla as "Others"
otu_data$Phylum <- ifelse(otu_data$Phylum %in% top_phyla, otu_data$Phylum, "Others")

# Step 5: Reorder phyla within each group by relative abundance
otu_data <- otu_data %>%
  group_by(SERUM_LPS_Group, Phylum) %>%
  mutate(total_abundance = sum(Abundance)) %>%  # Calculate the total abundance within each group
  ungroup() %>%
  arrange(SERUM_LPS_Group, desc(total_abundance)) %>%  # Sort by relative abundance in descending order
  mutate(Phylum = factor(Phylum, levels = unique(Phylum)))  # Reorder the factor levels by abundance

# Step 6: Generate the stacked bar plot with phyla ordered by abundance
ggplot(otu_data, aes(x = SERUM_LPS_Group, y = Abundance, fill = Phylum)) +
  geom_bar(stat = "identity", position = "stack") +  # Stacked bars for each group
  theme_minimal() +
  labs(x = "Serum LPS Group", y = "Relative Abundance", fill = "Phylum") +
  scale_fill_manual(values = c(
    "red", "blue", "green", "purple", "orange", "grey", "pink", "cyan", "brown", "yellow", "lightgreen"
  )) +  # Manually specify 11 colors
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  # Rotate x-axis labels for clarity



# Step 1: Group by SERUM_LPS_Group and Phylum, then calculate mean and sd for relative abundance
phylum_stats <- otu_data %>%
  group_by(SERUM_LPS_Group, Phylum) %>%
  summarise(
    Mean_Abundance = mean(Abundance),
    SD_Abundance = sd(Abundance),
    .groups = 'drop'
  )

# Step 2: Export the statistics to a CSV file
write.csv(phylum_stats, "phylum_gramnegativesERUM.csv", row.names = FALSE)




# Step 1: Convert counts to relative abundance
pseq_rel <- transform_sample_counts(gram_negative_pseq, function(x) x / sum(x))

# Step 2: Collapse by Genus level (instead of Phylum)
pseq_genus <- tax_glom(pseq_rel, taxrank = "Genus")

# Step 3: Convert to a data frame
otu_data_genus <- psmelt(pseq_genus)

# Step 4: Calculate total abundance and select top 10 genera
top_genera <- otu_data_genus %>%
  group_by(Genus) %>%
  summarise(Total_Abundance = sum(Abundance), .groups = 'drop') %>%
  top_n(10, Total_Abundance) %>%
  pull(Genus)

# Group the remaining genera as "Others"
otu_data_genus$Genus <- ifelse(otu_data_genus$Genus %in% top_genera, otu_data_genus$Genus, "Others")

# Step 5: Reorder genera within each group by relative abundance
otu_data_genus <- otu_data_genus %>%
  group_by(SERUM_LPS_Group, Genus) %>%
  mutate(total_abundance = sum(Abundance)) %>%  # Calculate the total abundance within each group
  ungroup() %>%
  arrange(SERUM_LPS_Group, desc(total_abundance)) %>%  # Sort by relative abundance in descending order
  mutate(Genus = factor(Genus, levels = unique(Genus)))  # Reorder the factor levels by abundance

# Step 6: Generate the stacked bar plot with genera ordered by abundance
ggplot(otu_data_genus, aes(x = SERUM_LPS_Group, y = Abundance, fill = Genus)) +
  geom_bar(stat = "identity", position = "stack") +  # Stacked bars for each group
  theme_minimal() +
  labs(x = "Serum LPS Group", y = "Relative Abundance", fill = "Genus") +
  scale_fill_manual(values = c(
    "red", "blue", "green", "purple", "orange", "grey", "pink", "cyan", "brown", "yellow", "lightgreen"
  )) +  # Manually specify 11 colors
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  # Rotate x-axis labels for clarity


# Step 1: Group by SERUM_LPS_Group and Genus, then calculate mean and sd for relative abundance
genus_stats <- otu_data_genus %>%
  group_by(SERUM_LPS_Group, Genus) %>%
  summarise(
    Mean_Abundance = mean(Abundance),
    SD_Abundance = sd(Abundance),
    .groups = 'drop'
  )

# Step 2: Export the statistics to a CSV file
write.csv(genus_stats, "GENUS_gramnegativesERUM.csv", row.names = FALSE)

# Install and load the necessary packages
#install.packages("ggvenn")  # If ggvenn is not installed already
library(ggvenn)




table(meta(pseq)$salivaLPS_Group, useNA = "always")
pseq.rel <- microbiome::transform(pseq, "compositional")

disease_states <- unique(as.character(meta(pseq.rel)$salivaLPS_Group))
print(disease_states)

list_core <- c() # an empty object to store information

for (n in disease_states){ # for each variable n in DiseaseState
    #print(paste0("Identifying Core Taxa for ", n))
    
    ps.sub <- subset_samples(pseq.rel, salivaLPS_Group == n) # Choose sample from DiseaseState by n
    
    core_m <- core_members(ps.sub, # ps.sub is phyloseq selected with only samples from g 
                           detection = 0.001, # 0.001 in atleast 90% samples 
                           prevalence = 0)
    print(paste0("No. of core taxa in ", n, " : ", length(core_m))) # print core taxa identified in each DiseaseState.
    list_core[[n]] <- core_m # add to a list core taxa for each group.
    #print(list_core)
}

library(ggvenn)

# Extract core taxa for each group from the list
core_taxa_A_Low <- list_core$A_Low
core_taxa_B_Medium <- list_core$B_Medium
core_taxa_C_High <- list_core$C_High



# Define the color scheme (red, blue, green)
mycols <- c("lightgreen", "lightblue", "yellow")

# Create the ggvenn_data list with group names that match those in mycols
ggvenn_data <- list(A_Low = core_taxa_A_Low,
                    B_Medium = core_taxa_B_Medium,
                    C_High = core_taxa_C_High)

# Create the Venn diagram using ggvenn with the updated colors
ggvenn(ggvenn_data, 
       fill_color = mycols,  # Apply the color scheme (no names)
       show_percentage = FALSE)  # Option to hide percentages






print(list_core)

# Find the shared core taxa between two groups (e.g., A_Low and B_Medium)
shared_A_Low_B_Medium <- intersect(list_core$A_Low, list_core$B_Medium)
print(paste0("Shared core taxa between A_Low and B_Medium: ", length(shared_A_Low_B_Medium)))
print(shared_A_Low_B_Medium)

# Find the shared core taxa between three groups (A_Low, B_Medium, and C_High)
shared_all_groups <- Reduce(intersect, list_core)  # This finds common taxa across all groups
print(paste0("Shared core taxa among all groups: ", length(shared_all_groups)))
print(shared_all_groups)


unique_A_Low <- setdiff(list_core$A_Low, union(list_core$B_Medium, list_core$C_High))
print(paste0("Unique core taxa in A_Low: ", length(unique_A_Low)))
print(unique_A_Low)

print(unique_A_Low)



unique_A_Low


# Unique taxa in B_Medium (taxa that are only in B_Medium, not in A_Low or C_High)
unique_B_Medium <- setdiff(list_core$B_Medium, union(list_core$A_Low, list_core$C_High))
print(paste0("Unique core taxa in B_Medium: ", length(unique_B_Medium)))
print(unique_B_Medium)


shared_taxa_names <- taxa_table[unique_B_Medium, ]

# Print the species or taxonomic names
print(shared_taxa_names)

# Unique taxa in C_High (taxa that are only in C_High, not in A_Low or B_Medium)
unique_C_High <- setdiff(list_core$C_High, union(list_core$A_Low, list_core$B_Medium))
print(paste0("Unique core taxa in C_High: ", length(unique_C_High)))
print(unique_C_High)


shared_taxa_names <- taxa_table[unique_C_High, ]

# Print the species or taxonomic names
print(shared_taxa_names)




# Taxa present in both A_Low and B_Medium but not in C_High
shared_A_Low_B_Medium_not_C_High <- setdiff(intersect(list_core$A_Low, list_core$B_Medium), list_core$C_High)
print(paste0("Taxa present in both A_Low and B_Medium but not in C_High: ", length(shared_A_Low_B_Medium_not_C_High)))
print(shared_A_Low_B_Medium_not_C_High)


# Taxa present in both A_Low and C_High but not in B_Medium
shared_A_Low_C_High_not_B_Medium <- setdiff(intersect(list_core$A_Low, list_core$C_High), list_core$B_Medium)
print(paste0("Taxa present in both A_Low and C_High but not in B_Medium: ", length(shared_A_Low_C_High_not_B_Medium)))
print(shared_A_Low_C_High_not_B_Medium)


# Taxa present in both B_Medium and C_High but not in A_Low
shared_B_Medium_C_High_not_A_Low <- setdiff(intersect(list_core$B_Medium, list_core$C_High), list_core$A_Low)
print(paste0("Taxa present in both B_Medium and C_High but not in A_Low: ", length(shared_B_Medium_C_High_not_A_Low)))
print(shared_B_Medium_C_High_not_A_Low)



# Get the taxonomic information from the phyloseq object
taxa_table <- tax_table(pseq)

# Retrieve species/taxonomic names for the taxa IDs
shared_taxa_names <- taxa_table[shared_B_Medium_C_High_not_A_Low, ]


shared_taxa_names <- taxa_table[shared_A_Low_B_Medium_not_C_High, ]

# Print the species or taxonomic names
print(shared_taxa_names)



library(writexl)


# Get the taxonomy table from the phyloseq object
tax_table_data <- tax_table(pseq.rel)  # Ensure you're using the correct phyloseq object

# Convert the taxonomy table to a data frame for easier manipulation
tax_df <- as.data.frame(tax_table_data)

# Add the taxon IDs (which are rownames of the taxonomy table) as a column in the taxonomy data frame
tax_df$Taxon_ID <- rownames(tax_df)

# Now filter the taxonomy for the shared taxa
shared_taxa_names <- tax_df[rownames(tax_df) %in% shared_all_groups, ]

# Optionally, if you want only genus and species names along with IDs
genus_species_names <- shared_taxa_names[, c("Taxon_ID", "Genus", "Species")]




write.csv(shared_taxa_names, file = "shared_taxa_names_with_ids.csv", row.names = FALSE)
# Or, if you want only genus and species:
write.csv(genus_species_names, file = "shared_taxa_names_with_ids.csv", row.names = FALSE)


shared_taxa_names <- taxa_table[unique_A_Low, ]

print(shared_taxa_names)
# Collapsing the taxa to the genus level
pseq_genus <- tax_glom(pseq.rel, "Genus")

# Iterate over the disease states (same as before)
list_core_genus <- c()  # Empty list to store core genera

for (n in disease_states) {
    # Subset the data by each disease state
    ps.sub <- subset_samples(pseq_genus, salivaLPS_Group == n)
    
    # Identify core members at the genus level (prevalence threshold can be adjusted as needed)
    core_m_genus <- core_members(ps.sub, detection = 0.001, prevalence = 0)
    
    # Print the number of core genera identified
    print(paste0("No. of core genera in ", n, " : ", length(core_m_genus)))
    
    # Add core genera to the list
    list_core_genus[[n]] <- core_m_genus
}

# Extract core genera for each group
core_genera_A_Low <- list_core_genus$A_Low
core_genera_B_Medium <- list_core_genus$B_Medium
core_genera_C_High <- list_core_genus$C_High

# Define the color scheme (you can modify this)
mycols <- c("lightgreen", "lightblue", "yellow")

# Create the ggvenn_data list for genus level
ggvenn_data_genus <- list(A_Low = core_genera_A_Low,
                          B_Medium = core_genera_B_Medium,
                          C_High = core_genera_C_High)

# Create the Venn diagram at the genus level
ggvenn(ggvenn_data_genus, 
       fill_color = mycols,  # Apply the color scheme
       show_percentage = FALSE)  # Option to hide percentages


# Identify genera unique to B_Medium
unique_genus_B_Medium <- setdiff(core_genera_B_Medium, 
                                 c(core_genera_A_Low, core_genera_C_High))

# Print the unique genera
print("Genera unique to B_Medium:")
print(unique_genus_B_Medium)

# Extract genus names using the taxon IDs
unique_genus_names_B_Medium <- taxa_table[unique_genus_B_Medium, "Genus"]

# Print the unique genera names
print("Unique genera in B_Medium:")
print(unique_genus_names_B_Medium)

 
# Identify genera unique to A_Low
unique_genus_A_Low <- setdiff(core_genera_A_Low, 
                              c(core_genera_B_Medium, core_genera_C_High))

# Identify genera unique to C_High
unique_genus_C_High <- setdiff(core_genera_C_High, 
                               c(core_genera_A_Low, core_genera_B_Medium))

# Extract genus names using the taxon IDs
unique_genus_names_A_Low <- taxa_table[unique_genus_A_Low, "Genus"]
unique_genus_names_B_Medium <- taxa_table[unique_genus_B_Medium, "Genus"]
unique_genus_names_C_High <- taxa_table[unique_genus_C_High, "Genus"]

# Print the unique genera names for each group
print("Unique genera in A_Low:")
print(unique_genus_names_A_Low)

print("Unique genera in B_Medium:")
print(unique_genus_names_B_Medium)

print("Unique genera in C_High:")
print(unique_genus_names_C_High)


# Find genera common to A_Low and C_High
shared_A_C <- intersect(core_genera_A_Low, core_genera_C_High)

# Remove those also present in B_Medium
shared_A_C_not_B <- setdiff(shared_A_C, core_genera_B_Medium)

# Extract genus names using the taxon IDs
shared_genus_names_A_C_not_B <- taxa_table[shared_A_C_not_B, "Genus"]

# Print the genera names
print("Genera present in both A_Low and C_High but not in B_Medium:")
print(shared_genus_names_A_C_not_B)


# Find genera common to A_Low and B_Medium
shared_A_B <- intersect(core_genera_A_Low, core_genera_B_Medium)

# Remove those also present in C_High
shared_A_B_not_C <- setdiff(shared_A_B, core_genera_C_High)

# Extract genus names using the taxon IDs
shared_genus_names_A_B_not_C <- taxa_table[shared_A_B_not_C, "Genus"]

# Print the genera names
print("Genera present in both A_Low and B_Medium but not in C_High:")
print(shared_genus_names_A_B_not_C)


# Find genera common to C_High and B_Medium
shared_C_B <- intersect(core_genera_C_High, core_genera_B_Medium)

# Remove those also present in A_Low
shared_C_B_not_A <- setdiff(shared_C_B, core_genera_A_Low)

# Extract genus names using the taxon IDs
shared_genus_names_C_B_not_A <- taxa_table[shared_C_B_not_A, "Genus"]

# Print the genera names
print("Genera present in both C_High and B_Medium but not in A_Low:")
print(shared_genus_names_C_B_not_A)

table(meta(pseq)$salivaLPS_Group, useNA = "always")
pseq.rel <- microbiome::transform(pseq, "compositional")

disease_states <- unique(as.character(meta(pseq.rel)$salivaLPS_Group))
print(disease_states)

list_core <- c() # an empty object to store information

for (n in disease_states){ # for each variable n in DiseaseState
    #print(paste0("Identifying Core Taxa for ", n))
    
    ps.sub <- subset_samples(pseq.rel, salivaLPS_Group == n) # Choose sample from DiseaseState by n
    
    core_m <- core_members(ps.sub, # ps.sub is phyloseq selected with only samples from g 
                           detection = 0.001, # 0.001 in atleast 90% samples 
                           prevalence = 0)
    print(paste0("No. of core taxa in ", n, " : ", length(core_m))) # print core taxa identified in each DiseaseState.
    list_core[[n]] <- core_m # add to a list core taxa for each group.
    #print(list_core)
}

library(ggvenn)

# Extract core taxa for each group from the list
core_taxa_A_Low <- list_core$A_Low
core_taxa_B_Medium <- list_core$B_Medium
core_taxa_C_High <- list_core$C_High



# Define the color scheme (red, blue, green)
mycols <- c("lightgreen", "lightblue", "yellow")

# Create the ggvenn_data list with group names that match those in mycols
ggvenn_data <- list(A_Low = core_taxa_A_Low,
                    B_Medium = core_taxa_B_Medium,
                    C_High = core_taxa_C_High)

# Create the Venn diagram using ggvenn with the updated colors
ggvenn(ggvenn_data, 
       fill_color = mycols,  # Apply the color scheme (no names)
       show_percentage = FALSE) 

       # Collapsing the taxa to the genus level
pseq_genus <- tax_glom(pseq.rel, "Genus")

# Iterate over the serum LPS groups
list_core_genus <- c()  # Empty list to store core genera

for (n in disease_states) {
    # Subset the data by each serum LPS group
    ps.sub <- subset_samples(pseq_genus, SERUM_LPS_Group == n)
    
    # Identify core genera (prevalence threshold can be adjusted)
    core_m_genus <- core_members(ps.sub, detection = 0.001, prevalence = 0)
    
    # Print the number of core genera identified
    print(paste0("No. of core genera in ", n, " : ", length(core_m_genus)))
    
    # Store core genera for each group
    list_core_genus[[n]] <- core_m_genus
}

# Extract core genera for each serum LPS group
core_genera_A_Low <- list_core_genus$A_Low
core_genera_B_Medium <- list_core_genus$B_Medium
core_genera_C_High <- list_core_genus$C_High

# Define colors for Venn diagram
mycols <- c("lightgreen", "lightblue", "yellow")

# Create ggvenn data list for genus level
ggvenn_data_genus <- list(A_Low = core_genera_A_Low,
                          B_Medium = core_genera_B_Medium,
                          C_High = core_genera_C_High)

# Create the Venn diagram at the genus level
ggvenn(ggvenn_data_genus, 
       fill_color = mycols,  
       show_percentage = FALSE)  

### **Identify Unique Genera**
# Unique to B_Medium
unique_genus_B_Medium <- setdiff(core_genera_B_Medium, 
                                 c(core_genera_A_Low, core_genera_C_High))

# Unique to A_Low
unique_genus_A_Low <- setdiff(core_genera_A_Low, 
                              c(core_genera_B_Medium, core_genera_C_High))

# Unique to C_High
unique_genus_C_High <- setdiff(core_genera_C_High, 
                               c(core_genera_A_Low, core_genera_B_Medium))

# Extract genus names from taxa table
unique_genus_names_A_Low <- taxa_table[unique_genus_A_Low, "Genus"]
unique_genus_names_B_Medium <- taxa_table[unique_genus_B_Medium, "Genus"]
unique_genus_names_C_High <- taxa_table[unique_genus_C_High, "Genus"]

# Print unique genera
print("Unique genera in A_Low:")
print(unique_genus_names_A_Low)

print("Unique genera in B_Medium:")
print(unique_genus_names_B_Medium)

print("Unique genera in C_High:")
print(unique_genus_names_C_High)

### **Shared Genera Between Two Groups (Excluding Third)**
# Genera in both A_Low and C_High but not in B_Medium
shared_A_C <- intersect(core_genera_A_Low, core_genera_C_High)
shared_A_C_not_B <- setdiff(shared_A_C, core_genera_B_Medium)
shared_genus_names_A_C_not_B <- taxa_table[shared_A_C_not_B, "Genus"]

print("Genera present in both A_Low and C_High but not in B_Medium:")
print(shared_genus_names_A_C_not_B)

# Genera in both A_Low and B_Medium but not in C_High
shared_A_B <- intersect(core_genera_A_Low, core_genera_B_Medium)
shared_A_B_not_C <- setdiff(shared_A_B, core_genera_C_High)
shared_genus_names_A_B_not_C <- taxa_table[shared_A_B_not_C, "Genus"]

print("Genera present in both A_Low and B_Medium but not in C_High:")
print(shared_genus_names_A_B_not_C)

# Genera in both C_High and B_Medium but not in A_Low
shared_C_B <- intersect(core_genera_C_High, core_genera_B_Medium)
shared_C_B_not_A <- setdiff(shared_C_B, core_genera_A_Low)
shared_genus_names_C_B_not_A <- taxa_table[shared_C_B_not_A, "Genus"]

print("Genera present in both C_High and B_Medium but not in A_Low:")
print(shared_genus_names_C_B_not_A)


saliva_results_df <- data.frame(Genus = character(), Estimate = numeric(), StdError = numeric(), tValue = numeric(), pValue = numeric())
serum_results_df <- data.frame(Genus = character(), Estimate = numeric(), tValue = numeric(), pValue = numeric())

# Loop through each genus and extract summary statistics
for (genus in genera) {
  
# Create an empty data frame to store summary results for both salivaLPS and serumLPS models
full_results_df <- data.frame(Genus = character(), Predictor = character(), Estimate = numeric(), StdError = numeric(), tValue = numeric(), pValue = numeric())

# Loop through each genus to fit the full model with the genus as an independent variable and LPS as the dependent variable
for (genus in genera) {
  
  # Create a formula dynamically for the lm() call to include genus as a predictor for salivaLPS
  formula_saliva <- as.formula(paste("salivaLPS ~", genus, "+", paste(genera[genera != genus], collapse = " + ")))
  
  # Fit the linear model for salivaLPS
  model_saliva <- lm(formula_saliva, data = abundance_data_merged_genus)
  
  # Extract the summary for the salivaLPS model
  if (inherits(model_saliva, "lm")) {
    saliva_summary <- summary(model_saliva)$coefficients
    # Filter for relevant rows (e.g., only genus coefficients)
    if (genus %in% rownames(saliva_summary)) {
      saliva_summary_df <- data.frame(saliva_summary[genus, , drop = FALSE])  
      saliva_summary_df$Genus <- genus  # Add genus name to the summary
      saliva_summary_df$Predictor <- "salivaLPS"  # Label the dependent variable
      full_results_df <- rbind(full_results_df, saliva_summary_df)  # Bind the results
    }
  }
  
  # Create a formula dynamically for the lm() call to include genus as a predictor for SERUM_LPS
  formula_serum <- as.formula(paste("SERUM_LPS ~", genus, "+", paste(genera[genera != genus], collapse = " + ")))
  
  # Fit the linear model for SERUM_LPS
  model_serum <- lm(formula_serum, data = abundance_data_merged_genus)
  
  # Extract the summary for the SERUM_LPS model
  if (inherits(model_serum, "lm")) {
    serum_summary <- summary(model_serum)$coefficients
    # Filter for relevant rows (e.g., only genus coefficients)
    if (genus %in% rownames(serum_summary)) {
      serum_summary_df <- data.frame(serum_summary[genus, , drop = FALSE])  
      serum_summary_df$Genus <- genus  # Add genus name to the summary
      serum_summary_df$Predictor <- "SERUM_LPS"  # Label the dependent variable
      full_results_df <- rbind(full_results_df, serum_summary_df)  # Bind the results
    }
  }
}

# View the results
print(head(full_results_df))

# Save the results to a CSV file
write.csv(full_results_df, "full_model_results.csv", row.names = FALSE)

# Optional: If you want to include the model summary (R², residuals, etc.)
# For the final model (after loop completion)
final_model_saliva <- lm(as.formula(paste("salivaLPS ~", paste(genera, collapse = " + "))), data = abundance_data_merged_genus)
final_model_serum <- lm(as.formula(paste("SERUM_LPS ~", paste(genera, collapse = " + "))), data = abundance_data_merged_genus)

# Print overall model summary (for residuals, R², F-statistic) for salivaLPS
final_model_saliva_summary <- summary(final_model_saliva)
print(final_model_saliva_summary)

# Print overall model summary (for residuals, R², F-statistic) for SERUM_LPS
final_model_serum_summary <- summary(final_model_serum)
print(final_model_serum_summary)

# Extract the coefficients from the model summary for both salivaLPS and serumLPS
final_model_saliva_coeffs <- as.data.frame(final_model_saliva_summary$coefficients)
final_model_serum_coeffs <- as.data.frame(final_model_serum_summary$coefficients)

# Save the coefficients of the model to CSV
write.csv(final_model_saliva_coeffs, "final_model_saliva_coefficients.csv", row.names = TRUE)
write.csv(final_model_serum_coeffs, "final_model_serum_coefficients.csv", row.names = TRUE)

# To save other components of the summary (e.g., residuals, R-squared, F-statistic)
saliva_summary_info <- data.frame(
  ResidualStdError = final_model_saliva_summary$sigma,
  MultipleRSquared = final_model_saliva_summary$r.squared,
  AdjustedRSquared = final_model_saliva_summary$adj.r.squared,
  Fstatistic = final_model_saliva_summary$fstatistic[1],
  PvalueFstatistic = final_model_saliva_summary$fstatistic[4]
)

serum_summary_info <- data.frame(
  ResidualStdError = final_model_serum_summary$sigma,
  MultipleRSquared = final_model_serum_summary$r.squared,
  AdjustedRSquared = final_model_serum_summary$adj.r.squared,
  Fstatistic = final_model_serum_summary$fstatistic[1],
  PvalueFstatistic = final_model_serum_summary$fstatistic[4]
)

# Save the model information (residuals, R², F-statistic) to CSV
write.csv(saliva_summary_info, "final_model_saliva_summary_info.csv", row.names = FALSE)
write.csv(serum_summary_info, "final_model_serum_summary_info.csv", row.names = FALSE)
library(microbiomeMarker)

ps<-pseq_subsetA_LowB_MediumSERUM

lefse_final<-run_lefse(
  ps,
  group = "SERUM_LPS_Group",
  subgroup = NULL,
  taxa_rank = "all",
  transform = c("identity", "log10", "log10p"),
  norm = "CPM",
  norm_para = list(),
  kw_cutoff = 0.05,
  lda_cutoff = 3,
  bootstrap_n = 30,
  bootstrap_fraction = 2/3,
  wilcoxon_cutoff = 0.05,
  multigrp_strat = FALSE,
  strict = c("0"),
  sample_min = 0,
  only_same_subgrp = FALSE,
  curv = FALSE
)

write.csv(marker_table(lefse_final), "LEFSE_LowB_MediumSERUMphylum.csv")



library(microbiomeMarker)
ps<-pseq_subsetA_LowB_C_HighSERUM

lefse_final3<-run_lefse(
  ps,
  group = "SERUM_LPS_Group",
  subgroup = NULL,
  taxa_rank = "all",
  transform = c("identity", "log10", "log10p"),
  norm = "CPM",
  norm_para = list(),
  kw_cutoff = 0.05,
  lda_cutoff = 3,
  bootstrap_n = 30,
  bootstrap_fraction = 2/3,
  wilcoxon_cutoff = 0.05,
  multigrp_strat = FALSE,
  strict = c("0"),
  sample_min = 0,
  only_same_subgrp = FALSE,
  curv = FALSE
)

write.csv(marker_table(lefse_final3), "LEFSE_LowB_HighSERUMphylum.csv")




library(microbiomeMarker)
ps<-pseq_subsetB_MediumB_C_HighSERUM

lefse_final6<-run_lefse(
  ps,
  group = "SERUM_LPS_Group",
  subgroup = NULL,
  taxa_rank = "all",
  transform = c("identity", "log10", "log10p"),
  norm = "CPM",
  norm_para = list(),
  kw_cutoff = 0.05,
  lda_cutoff = 3,
  bootstrap_n = 30,
  bootstrap_fraction = 2/3,
  wilcoxon_cutoff = 0.05,
  multigrp_strat = FALSE,
  strict = c("0"),
  sample_min = 0,
  only_same_subgrp = FALSE,
  curv = FALSE
)

write.csv(marker_table(lefse_final6), "LEFSE_medium_HighSERUMphylum.csv")

HUMAnN_renorm_pathabundance <- read_delim("HUMAnN_renorm_pathabundance.tsv", delim = "\t", escape_double = FALSE, trim_ws = TRUE)
pathabundance <-HUMAnN_renorm_pathabundance  

#transpose row to coloumn
pathabundance1 <- as.data.frame(t(pathabundance))
View(pathabundance1)


meta_308_pathabundance <- read.delim("P:/h305/periocardio/Secreto/Metagenomics/kraken-biom/Bracken/braken_biom/with_blank/Antibiotics/main_study/mia/LPS/charatertics/pathway/meta_308_pathabundance.txt", row.names=1)




rownames(pathabundance1[1:10,])
rownames(meta_308_pathabundance[1:10,])
colnames(pathabundance1)<-pathabundance1[1,]
pathabundance1 <- pathabundance1[2:nrow(pathabundance1),]
colnames(pathabundance1)
colnames(meta_308_pathabundance)

fit_data <- Maaslin2(
  input_data = input_data,
  input_metadata = metadata,
  output = output_dir,
  fixed_effects = c("SERUM_LPS_Group"),
  random_effects = NULL,
  normalization = "TSS",
  transform = "LOG",
  analysis_method = "LM",
  min_prevalence = 0.1,
  standardize = FALSE
)

library(ggplot2)

# Create a volcano plot (log10 scale for p-values)
ggplot(all_results, aes(x = coef, y = -log10(pval))) +
  geom_point(aes(color = pval < 0.05), size = 1.5) +
  scale_color_manual(values = c("gray", "red")) +
  labs(x = "Coefficient", y = "-log10(p-value)", title = "Volcano Plot of Pathway Associations") +
  theme_minimal()



# Create a bar plot for significant pathways
ggplot(significant_results, aes(x = reorder(feature, coef), y = coef)) +
  geom_bar(stat = "identity", fill = "skyblue") +
  coord_flip() +
  labs(x = "Pathway", y = "Coefficient", title = "Significant Pathways (Coef)") +
  theme_minimal()

# Create the bar plot with numerical labels and different colors for positive and negative coefficients
ggplot(significant_results16, aes(x = reorder(label, coef), y = coef, fill = color)) +
  geom_bar(stat = "identity") +
  coord_flip() +  # Flip the axes to make the labels readable
  labs(x = "Pathway", y = "Coefficient", title = "Significant Pathways (Coef)") +
  scale_fill_identity() +  # Use the assigned colors directly
  theme_minimal() +
  theme(
    axis.text = element_text(size = 10),
    axis.title = element_text(size = 12),
    plot.title = element_text(hjust = 0.5)
  )
