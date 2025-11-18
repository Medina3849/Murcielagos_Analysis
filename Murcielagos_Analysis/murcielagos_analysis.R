# =============================================================================
# PHYSICO-CHEMICAL ANALYSIS OF ARCHAEOLOGICAL SEDIMENTS
# CUEVA DE LOS MURCIÉLAGOS ARCHAEOLOGICAL SITE
# =============================================================================

rm(list = ls())

# Load required packages
required_packages <- c("tidyverse", "ggpubr", "FactoMineR", "factoextra", 
                       "viridis", "patchwork", "rstatix", "corrplot",
                       "scales", "ggrepel", "dendextend", "cluster")

new_packages <- required_packages[!(required_packages %in% installed.packages()[,"Package"])]
if(length(new_packages)) install.packages(new_packages, dependencies = TRUE)

invisible(lapply(required_packages, library, character.only = TRUE))

# =============================================================================
# DATA INTEGRATION
# =============================================================================

murcielagos_data <- data.frame(
  Sample = c('MZ-13', 'MZ-12', 'MZ-11', 'MZ-10', 'MZ-9', 'MZ-8', 'MZ-7', 'MZ-6', 'MZ-5', 'MZ-4', 'MZ-3', 'MZ-2', 'MZ-1'),
  Period = c('Classical Period', 'Bronze Age', 'Chalcolithic', 'Neolithic', 'Neolithic', 'Neolithic', 'Neolithic', 'Neolithic', 'Neolithic', 'Palaeolithic', 'Palaeolithic', 'Palaeolithic', 'Palaeolithic'),
  Hum = c(6.23, 5.58, 4.59, 5.03, 4.69, 4.83, 4.29, 4.38, 3.94, 2.69, 2.84, 1.65, 5.78),
  LOI = c(12.53, 10.92, 9.19, 9.36, 9.78, 11.57, 8.92, 6.37, 5.48, 3.79, 4.53, 3.24, 6.48),
  EC = c(2.22, 2.07, 1.49, 1.26, 1.28, 1.46, 1.09, 1.79, 1.29, 1.14, 1.01, 0.79, 0.94),
  CO3 = c(15, 14, 30, 24, 34, 37, 38, 38, 32, 56, 53, 58, 11),
  M.S. = c(1980, 1890, 1615, 1670, 1425, 1325, 1355, 1600, 1745, 1415, 1345, 935, 4200),
  pH = c(7.8, 7.9, 8.1, 7.8, 7.9, 7.9, 7.9, 8, 8.2, 8, 7.9, 8, 8.3),
  OrgC = c(4.08, 3.62, 5.08, 5.23, 5.39, 9.24, 3.77, 2, 1.23, 0.61, 0.85, 0, 0.54),
  Sand = c(34.15, 39.4, 34.5, 42.75, 42.6, 45.9, 35.75, 34.95, 34.75, 41.9, 39.1, 43.7, 43.45),
  Silt = c(53.35, 48.1, 53, 47.25, 47.4, 46.6, 46.75, 50.05, 47.75, 38.1, 43.4, 41.3, 44.05),
  Clay = c(12.5, 12.5, 12.5, 10, 10, 7.5, 17.5, 15, 17.5, 20, 17.5, 15, 12.5)
) %>%
  mutate(
    Period = factor(Period, levels = c('Palaeolithic', 'Neolithic', 'Chalcolithic', 'Bronze Age', 'Classical Period')),
    Sample = factor(Sample)
  ) %>%
  rename('% CO3' = CO3)

# Color scheme
period_colors <- c("#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd")
names(period_colors) <- levels(murcielagos_data$Period)

output_dir <- "Murcielagos_Analysis_Results"
if (!dir.exists(output_dir)) dir.create(output_dir)

# =============================================================================
# STATISTICAL FUNCTIONS
# =============================================================================

calculate_statistical_results <- function(data) {
  results <- data %>%
    group_by(Parameter) %>%
    kruskal_test(Value ~ Period) %>%
    adjust_pvalue(method = "BH")
  
  results %>%
    mutate(p.signif = case_when(
      p.adj < 0.001 ~ "***",
      p.adj < 0.01 ~ "**",
      p.adj < 0.05 ~ "*",
      p.adj < 0.1 ~ ".",
      TRUE ~ "ns"
    ))
}

# =============================================================================
# FIGURE 1: PARAMETERS VARIATION
# =============================================================================

selected_params <- c("% CO3", "LOI", "M.S.", "OrgC")
murcielagos_selected <- murcielagos_data %>%
  select(Sample, Period, all_of(selected_params)) %>%
  pivot_longer(cols = all_of(selected_params), 
               names_to = "Parameter", 
               values_to = "Value") %>%
  mutate(Parameter = factor(Parameter, levels = selected_params))

selected_statistical_results <- calculate_statistical_results(murcielagos_selected)

create_parameter_plot <- function(param_data, param_name) {
  y_range <- range(param_data$Value, na.rm = TRUE)
  y_pos <- y_range[2] + 0.15 * diff(y_range)
  
  p_val <- selected_statistical_results %>%
    filter(Parameter == param_name) %>%
    pull(p.adj) %>%
    round(3)
  
  p_text <- ifelse(p_val < 0.001, "p < 0.001", paste("p =", p_val))
  
  ggplot(param_data, aes(x = Period, y = Value, fill = Period)) +
    geom_boxplot(alpha = 0.8, outlier.shape = NA, width = 0.7, size = 0.5) +
    geom_jitter(width = 0.2, size = 2, alpha = 0.7, shape = 21, color = "black", stroke = 0.5) +
    scale_fill_manual(values = period_colors) +
    labs(title = param_name, x = "", y = "") +
    theme_pubr(base_size = 12) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
          plot.title = element_text(face = "bold", size = 12, hjust = 0.5),
          legend.position = "none") +
    annotate("text", x = 3, y = y_pos, label = p_text, 
             size = 4.5, fontface = "bold", color = "black")
}

plot_list <- map(selected_params, ~create_parameter_plot(
  filter(murcielagos_selected, Parameter == .x), .x
))

figure1 <- (plot_list[[1]] + plot_list[[2]]) / (plot_list[[3]] + plot_list[[4]]) +
  plot_annotation(title = "Key Physico-Chemical Parameters Variation Across Archaeological Periods",
                  theme = theme(plot.title = element_text(face = "bold", size = 16, hjust = 0.5)))

ggsave(file.path(output_dir, "Figure1_Parameters.tiff"), figure1, 
       device = "tiff", dpi = 600, width = 14, height = 12, compression = "lzw", bg = "white")

# =============================================================================
# FIGURE 2: PCA ANALYSIS
# =============================================================================

pca_matrix <- murcielagos_data[, 3:12]
rownames(pca_matrix) <- as.character(murcielagos_data$Sample)

pca_result <- PCA(pca_matrix, scale.unit = TRUE, graph = FALSE)

pca_coords <- as.data.frame(pca_result$ind$coord)
pca_coords$Sample <- as.character(murcielagos_data$Sample)
pca_coords$Period <- murcielagos_data$Period

hull_data <- pca_coords %>%
  group_by(Period) %>%
  slice(chull(Dim.1, Dim.2))

pca_individuals <- ggplot(pca_coords, aes(x = Dim.1, y = Dim.2, color = Period, fill = Period)) +
  geom_point(size = 6, alpha = 0.9, shape = 21, color = "black", stroke = 0.8) +
  geom_text(aes(label = Sample), size = 4.5, fontface = "bold", 
            position = position_nudge(y = 0.15, x = 0.15)) +
  geom_polygon(data = hull_data, alpha = 0.15, linetype = "solid", size = 0.8) +
  scale_color_manual(values = period_colors) +
  scale_fill_manual(values = period_colors) +
  labs(title = "PCA - Sample Distribution",
       x = paste0("Dim 1 (", round(pca_result$eig[1,2], 1), "%)"),
       y = paste0("Dim 2 (", round(pca_result$eig[2,2], 1), "%)")) +
  theme_pubr() +
  theme(legend.position = "right",
        plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
        axis.title = element_text(face = "bold", size = 12),
        panel.grid.major = element_line(color = "grey90", size = 0.3))

pca_variables <- fviz_pca_var(pca_result,
                              col.var = "contrib",
                              gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
                              repel = TRUE) +
  labs(title = "PCA - Variable Contributions") +
  theme_pubr() +
  theme(plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
        axis.title = element_text(face = "bold", size = 12))

figure2 <- pca_individuals + pca_variables + 
  plot_layout(widths = c(1.2, 1)) +
  plot_annotation(tag_levels = 'A') &
  theme(plot.tag = element_text(face = "bold", size = 16))

ggsave(file.path(output_dir, "Figure2_PCA.tiff"), figure2, 
       device = "tiff", dpi = 600, width = 18, height = 9, compression = "lzw", bg = "white")

# =============================================================================
# FIGURE 3: CLUSTERING ANALYSIS
# =============================================================================

clustering_data <- murcielagos_data[, 3:12]
rownames(clustering_data) <- as.character(murcielagos_data$Sample)

dist_matrix <- dist(scale(clustering_data), method = "euclidean")
hc <- hclust(dist_matrix, method = "ward.D2")

create_dendrogram <- function() {
  dend <- as.dendrogram(hc)
  
  sample_order <- labels(dend)
  period_order <- murcielagos_data$Period[match(sample_order, murcielagos_data$Sample)]
  
  labels_colors(dend) <- period_colors[period_order]
  
  dend <- dend %>%
    set("branches_lwd", 2.5) %>%
    set("branches_col", "gray40") %>%
    set("labels_cex", 1.3) %>%
    set("labels_col", period_colors[period_order]) %>%
    set("leaves_pch", 19) %>%
    set("leaves_cex", 1.2) %>%
    set("leaves_col", period_colors[period_order])
  
  par(mar = c(5, 8, 4, 2) + 0.1)
  
  plot(dend, 
       horiz = TRUE,
       main = "Hierarchical Clustering Analysis",
       xlab = "Euclidean Distance",
       ylab = "",
       cex.main = 1.8,
       cex.lab = 1.2,
       cex.axis = 1.1,
       font.main = 2,
       font.lab = 2,
       lwd = 2.5,
       axes = TRUE)
  
  legend("topleft", 
         legend = names(period_colors),
         fill = period_colors,
         border = NA,
         bty = "n",
         cex = 1.1,
         title = "Archaeological Period",
         title.font = 2)
}

tiff(file.path(output_dir, "Figure3_Clustering.tiff"), 
     width = 14, height = 10, units = "in", res = 600, compression = "lzw")
create_dendrogram()
dev.off()

silhouette_scores <- map_dbl(2:6, function(k) {
  clusters <- cutree(hc, k = k)
  silhouette_avg <- summary(silhouette(clusters, dist_matrix))$avg.width
  return(silhouette_avg)
})

optimal_k <- which.max(silhouette_scores) + 1
clusters <- cutree(hc, k = optimal_k)

cluster_info <- data.frame(
  Sample = murcielagos_data$Sample,
  Period = murcielagos_data$Period,
  Cluster = factor(clusters)
)

write.csv(cluster_info, file.path(output_dir, "cluster_assignments.csv"), row.names = FALSE)

# =============================================================================
# FIGURE 4: CORRELATION ANALYSIS
# =============================================================================

correlation_matrix <- cor(murcielagos_data[, 3:12], use = "complete.obs")

create_correlation_plot <- function() {
  n <- nrow(murcielagos_data)
  cor_p_values <- matrix(NA, nrow = nrow(correlation_matrix), ncol = ncol(correlation_matrix))
  colnames(cor_p_values) <- colnames(correlation_matrix)
  rownames(cor_p_values) <- rownames(correlation_matrix)
  
  for (i in 1:nrow(correlation_matrix)) {
    for (j in 1:ncol(correlation_matrix)) {
      if (i != j) {
        test <- cor.test(murcielagos_data[[rownames(correlation_matrix)[i]]], 
                         murcielagos_data[[colnames(correlation_matrix)[j]]])
        cor_p_values[i, j] <- test$p.value
      }
    }
  }
  
  color_matrix <- matrix("white", nrow = nrow(correlation_matrix), ncol = ncol(correlation_matrix))
  rownames(color_matrix) <- rownames(correlation_matrix)
  colnames(color_matrix) <- colnames(correlation_matrix)
  
  for (i in 1:nrow(correlation_matrix)) {
    for (j in 1:ncol(correlation_matrix)) {
      if (!is.na(cor_p_values[i, j]) && cor_p_values[i, j] < 0.01) {
        color_matrix[i, j] <- ifelse(correlation_matrix[i, j] > 0, "#B2182B", "#2166AC")
      }
    }
  }
  
  par(mar = c(2, 2, 2, 2))
  
  corrplot(correlation_matrix, 
           method = "color", 
           type = "upper",
           col = "white",
           bg = "white",
           tl.col = "black",
           tl.srt = 45,
           tl.cex = 1.1,
           cl.cex = 0.9,
           addCoef.col = "black",
           number.cex = 0.8,
           mar = c(0, 0, 0, 0),
           diag = FALSE)
  
  corrplot(correlation_matrix, 
           method = "color", 
           type = "upper",
           col = colorRampPalette(c("#2166AC", "#F7F7F7", "#B2182B"))(100),
           bg = NA,
           tl.col = "black",
           tl.srt = 45,
           tl.cex = 1.1,
           cl.cex = 0.9,
           addCoef.col = NA,
           number.cex = 0.8,
           mar = c(0, 0, 0, 0),
           diag = FALSE,
           add = TRUE)
}

tiff(file.path(output_dir, "Figure4_Correlations.tiff"), 
     width = 12, height = 10, units = "in", res = 600)
create_correlation_plot()
dev.off()

# =============================================================================
# DATA EXPORT
# =============================================================================

descriptive_stats <- murcielagos_selected %>%
  group_by(Period, Parameter) %>%
  summarise(
    n = n(),
    Mean = mean(Value),
    SD = sd(Value),
    Median = median(Value),
    IQR = IQR(Value),
    .groups = 'drop'
  )

publication_table <- descriptive_stats %>%
  left_join(selected_statistical_results %>% select(Parameter, p.adj, p.signif), 
            by = "Parameter") %>%
  mutate(Value_Summary = sprintf("%.2f ± %.2f", Mean, SD)) %>%
  select(Period, Parameter, Value_Summary, n, p.adj, p.signif) %>%
  pivot_wider(names_from = Period, values_from = Value_Summary) %>%
  arrange(factor(Parameter, levels = selected_params))

significant_cor_99 <- data.frame()
for (i in 1:(nrow(correlation_matrix)-1)) {
  for (j in (i+1):ncol(correlation_matrix)) {
    test <- cor.test(murcielagos_data[[rownames(correlation_matrix)[i]]], 
                     murcielagos_data[[colnames(correlation_matrix)[j]]])
    if (test$p.value < 0.01) {
      significant_cor_99 <- rbind(significant_cor_99, 
                                  data.frame(Var1 = rownames(correlation_matrix)[i],
                                             Var2 = colnames(correlation_matrix)[j],
                                             Correlation = correlation_matrix[i, j],
                                             p_value = test$p.value))
    }
  }
}

write.csv(correlation_matrix, file.path(output_dir, "correlation_matrix.csv"))
write.csv(significant_cor_99, file.path(output_dir, "significant_correlations_99.csv"), row.names = FALSE)
write.csv(selected_statistical_results, file.path(output_dir, "statistical_results.csv"), row.names = FALSE)
write.csv(descriptive_stats, file.path(output_dir, "descriptive_statistics.csv"), row.names = FALSE)
write.csv(publication_table, file.path(output_dir, "publication_table.csv"), row.names = FALSE)

cat("Analysis completed. Results saved in:", output_dir, "\n")
