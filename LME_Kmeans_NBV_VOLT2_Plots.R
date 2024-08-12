#Loop through different numbers of clusters of combined NBV and VOLT2 Intercepts and Slopes

library(ggplot2)
library(dplyr)
library(GGally)

output_dir <- "/data/users/rty388/testing/Cluster_NBV_VOLT2/All_patients/"

df <- read.csv("/data/users/rty388/testing/02_results.csv")

#Subset the core phase of the study for each subject 
df <- df %>% filter(CORE == "Yes")

df <- df %>%
  mutate(VOLT2CBRT = VOLT2^(1/3))

df2 <- read.csv(paste0(output_dir, "/NBV_VOLT2_Intercepts_Slopes.csv"))

# Function to run k-means clustering and generate plots
run_kmeans_and_plot <- function(k, df2, output_dir) {
  set.seed(123) # For reproducibility
  
  # Run k-means clustering
  km <- kmeans(df2[, c("Total_Intercept_NBV", "Total_Slope_NBV", "Total_Intercept_VOLT2", "Total_Slope_VOLT2")], centers = k, nstart = 25)
  
  # Add the cluster assignment to the dataframe
  df2$Cluster <- as.factor(km$cluster)
  
  write.csv(df2, paste0(output_dir, k, "kmeansclusters/", "kmeans", k, "Clusters_NBV_VOLT2.csv"), row.names = F)
  
  # Plot the scatterplot matrix
  df2_for_plot <- df2 %>%
    select(Total_Intercept_NBV, Total_Slope_NBV, Total_Intercept_VOLT2, Total_Slope_VOLT2, Cluster)
  
  scatterplot_matrix <- ggpairs(df2_for_plot, aes(color = Cluster, alpha = 0.6),
                                upper = list(continuous = wrap("points", size = 2)),
                                lower = list(continuous = wrap("points", size = 2)),
                                diag = list(continuous = wrap("densityDiag"))) +
    labs(title = paste("Scatterplot Matrix of Clustering Results for k =", k)) +
    theme_minimal()
  
  print(scatterplot_matrix)
  ggsave(filename = paste0(output_dir, k, "kmeansclusters/", "Scatterplot_Matrix_kmeansNBV_VOLT2_k", k, ".png"), plot = scatterplot_matrix, device = "png", width = 10, height = 8)
  ggsave(filename = paste0(output_dir, k, "kmeansclusters/", "Scatterplot_Matrix_kmeansNBV_VOLT2_k", k, ".pdf"), plot = scatterplot_matrix, device = "pdf", width = 11, height = 10)
  
  #Plot the Total NBV Intercept against Total slope for each subject split by cluster: 
  NBV_clusters <- ggplot(df2, aes(x = Total_Intercept_NBV, y = Total_Slope_NBV, color = Cluster)) +
    geom_point(size = 2) +
    labs(title = "NBV slope and intercept grouped by k-means combined NBV_VOLT2 clusters",
         x = "Total Intercept",
         y = "Total Slope") +
    theme_minimal()
  
  print(NBV_clusters)
  
  ggsave(filename = paste0(output_dir, k, "kmeansclusters/", "NBV_kmeans_k", k, ".png"), plot = NBV_clusters, device = "png", width = 10, height = 8)
  ggsave(filename = paste0(output_dir, k, "kmeansclusters/", "NBV_kmeans_k", k, ".pdf"), plot = NBV_clusters, device = "pdf", width = 10, height = 8)
  
  #Plot the Total VOLT2 Intercept against Total slope for each subject split by cluster: 
  VOLT2_clusters <- ggplot(df2, aes(x = Total_Intercept_VOLT2, y = Total_Slope_VOLT2, color = Cluster)) +
    geom_point(size = 2) +
    labs(title = "VOLT2 slope and intercept grouped by k-means combined NBV_VOLT2 clusters",
         x = "Total Intercept",
         y = "Total Slope") +
    theme_minimal()
  
  print(NBV_clusters)
  
  ggsave(filename = paste0(output_dir, k, "kmeansclusters/", "VOLT2_kmeans_k", k, ".png"), plot = VOLT2_clusters, device = "png", width = 10, height = 8)
  ggsave(filename = paste0(output_dir, k, "kmeansclusters/", "VOLT2_kmeans_k", k, ".pdf"), plot = VOLT2_clusters, device = "pdf", width = 11, height = 10)
  
  # Compute the principal components
  pca <- prcomp(df2[, c("Total_Intercept_NBV", "Total_Slope_NBV", "Total_Intercept_VOLT2", "Total_Slope_VOLT2")], scale. = TRUE)
  df2$pca1 <- pca$x[, 1]
  df2$pca2 <- pca$x[, 2]
  
  # Plot the clusters using principal components
  Cluster_Plot <- ggplot(df2, aes(x = pca1, y = pca2, color = Cluster)) +
    geom_point(size = 3) +
    labs(title = paste("K-means Clustering Results for k =", k),
         x = "Principal Component 1",
         y = "Principal Component 2") +
    theme_minimal()
  
  print(Cluster_Plot)
  ggsave(filename = paste0(output_dir, k, "kmeansclusters/", "Cluster_Plot_kmeansNBV_VOLT2_k", k, ".png"), plot = Cluster_Plot, device = "png", width = 10, height = 8)
  ggsave(filename = paste0(output_dir, k, "kmeansclusters/", "Cluster_Plot_kmeansNBV_VOLT2_k", k, ".pdf"), plot = Cluster_Plot, device = "pdf", width = 10, height = 8)
  
  # Plot trajectories for each cluster
  for (cluster in unique(df2$Cluster)) {
    cluster_df <- df2 %>% filter(Cluster == cluster)
    
    # Select 10 subjects per cluster
    top_subjects <- cluster_df %>%
      count(USUBJID, sort = TRUE) %>%
      head(20) %>%
      pull(USUBJID)
    
    df_top10 <- df %>% filter(USUBJID %in% top_subjects)
    
    # NBV trajectories
    nbv_trajectory_plot <- ggplot(data = df_top10, aes(x = DAY, y = NBV, group = USUBJID)) +
      geom_line(aes(color = USUBJID)) +
      geom_point(aes(color = USUBJID)) +
      facet_wrap(~USUBJID) +
      labs(title = paste("NBV Trajectories for Cluster", cluster, "with k =", k)) +
      theme_minimal()
    
    print(nbv_trajectory_plot)
    ggsave(filename = paste0(output_dir, k, "kmeansclusters/", "NBVtrajectories_", k, "clusters_cluster", cluster, ".png"), plot = nbv_trajectory_plot, device = "png", width = 10, height = 8)
    ggsave(filename = paste0(output_dir, k, "kmeansclusters/", "NBVtrajectories_", k, "clusters_cluster", cluster, ".pdf"), plot = nbv_trajectory_plot, device = "pdf", width = 11, height = 10)
    
    # VOLT2 trajectories
    volt2_trajectory_plot <- ggplot(data = df_top10, aes(x = DAY, y = VOLT2CBRT, group = USUBJID)) +
      geom_line(aes(color = USUBJID)) +
      geom_point(aes(color = USUBJID)) +
      facet_wrap(~USUBJID) +
      labs(title = paste("VOLT2 Trajectories for Cluster", cluster, "with k =", k)) +
      theme_minimal()
    
    print(volt2_trajectory_plot)
    ggsave(filename = paste0(output_dir, k, "kmeansclusters/", "VOLT2trajectories_", k, "clusters_cluster", cluster, ".png"), plot = volt2_trajectory_plot, device = "png", width = 10, height = 8)
    ggsave(filename = paste0(output_dir, k, "kmeansclusters/", "VOLT2trajectories_", k, "clusters_cluster", cluster, ".pdf"), plot = volt2_trajectory_plot, device = "pdf", width = 11, height = 10)
  }
}

# Main loop to run the analysis for different k values
k_values_to_test <- c(2)

for (k in k_values_to_test) {
  
  run_kmeans_and_plot(k, df2, output_dir)
}
