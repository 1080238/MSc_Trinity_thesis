library(ggplot2)
library(lmerTest)
library(GGally)
library(dplyr)
library(cluster)
library(rlang)


options(bitmapType = 'cairo')

output_dir <- "/data/users/rty388/testing/Cluster_NBV_VOLT2/"

# Run linear mixed effects model for NBV over time!
df <- read.csv("/data/users/rty388/testing/02_results.csv")

# Subset the core phase of the study for each subject 
df <- df %>% filter(CORE == "Yes")

df <- df %>% select(USUBJID, DAY, NBV, VOLT2)

# Omit rows with NA in specific column of data frame
df <- df[!(is.na(df$NBV)), ]

# Calculate cube root of VOLT2
df <- df %>% mutate(VOLT2CBRT = VOLT2^(1/3))

# Add Cluster information to df
df2 <- read.csv(paste0(output_dir, "2kmeansclusters/kmeans2Clusters_NBV_VOLT2.csv"))
df2 <- df2 %>% select(USUBJID, Cluster)
df <- df %>% inner_join(df2, by = "USUBJID")

# Function to create and save trajectory plots
create_trajectory_plot <- function(data, variable, cluster, k, num_subjects, min_points = 4) {
  cluster_df <- data %>% filter(Cluster == cluster)
  
  # Select subjects with at least min_points data points
  top_subjects <- cluster_df %>%
    group_by(USUBJID) %>%
    filter(n() >= min_points) %>%
    ungroup() %>%
    count(USUBJID, sort = TRUE) %>%
    head(num_subjects) %>%
    pull(USUBJID)
  
  df_top <- data %>% filter(USUBJID %in% top_subjects)
  
  # Plot trajectories
  trajectory_plot <- ggplot(data = df_top, aes(x = DAY, y = !!sym(variable), group = USUBJID)) +
    geom_line(aes(color = USUBJID)) +
    geom_point(aes(color = USUBJID)) +
    facet_wrap(~USUBJID, scales = "free_x") +
    labs(title = paste(variable, "Trajectories for Cluster", cluster, "with k =", k)) +
    theme_minimal() +
    theme(
      strip.text.x = element_text(size = 14),
      strip.text.y = element_text(size = 14),
      axis.text.x = element_text(size = 14, angle = 45, hjust = 1),
      axis.text.y = element_text(size = 14),
      axis.title.x = element_text(size = 16), 
      axis.title.y = element_text(size = 16)
    )
  
  print(trajectory_plot)
  
  # Ensure the plot is displayed
  dev.new()
  print(trajectory_plot)
  
  # Save the plot
  ggsave(filename = paste0(output_dir, "2kmeansclusters_2/", variable, "trajectories_", k, "clusters_cluster", cluster, ".png"), plot = trajectory_plot, device = "png", width = 10, height = 8)
  ggsave(filename = paste0(output_dir, "2kmeansclusters_2/", variable, "trajectories_", k, "clusters_cluster", cluster, ".pdf"), plot = trajectory_plot, device = "pdf", width = 11, height = 10)
}

# Create example plots for NBV
create_trajectory_plot(df, "NBV", 1, 2, 6)
create_trajectory_plot(df, "NBV", 2, 2, 6)

# Create example plots for VOLT2CBRT
create_trajectory_plot(df, "VOLT2CBRT", 1, 2, 6)
create_trajectory_plot(df, "VOLT2CBRT", 2, 2, 6)

