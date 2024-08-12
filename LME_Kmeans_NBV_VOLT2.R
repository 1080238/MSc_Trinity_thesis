##Run Linear mixed effects models on NBV and VOLT2
#https://kristopherkyle.github.io/IntroQuantALRM/8_Linear_Mixed_Effects_Models.html

#Figure out if you need to run this specifically on the Young (20-40) and Old (>=50) age group!!!!! 

library(ggplot2)
library(lmerTest)
library(GGally)
library(dplyr)
library(cluster)


options(bitmapType = 'cairo')

output_dir <- "/data/users/rty388/testing/Cluster_NBV_VOLT2/All_patients/NBV_VOLT2_clusters/"

#Run linear mixed effects model for NBV over time!
df <- read.csv("/data/users/rty388/testing/02_results.csv")

#Subset the core phase of the study for each subject 
df <- df %>% filter(CORE == "Yes")

#Calculate number of subjects
length(unique(df$USUBJID))

df <- df %>% select(USUBJID, DAY, NBV)

#omit rows with NA in specific column of data frame
df <- df[!(is.na(df$NBV)), ]

summary(df)

#------------------------------------------------------------------------------------------------
#Visualise the trajectories of NBV over time for 8 subjects

top_subjects <- df %>%
  count(USUBJID, sort = TRUE) %>%
  head(8) %>%
  pull(USUBJID)

#Subset the dataframe to include only the top 8 subjects
df_top8 <- df %>%
  filter(USUBJID %in% top_subjects)

summary(df_top8)

g3 <- ggplot(data = df_top8, aes(x = DAY, y = NBV, group = USUBJID)) +
  geom_line(aes(color=USUBJID)) +
  geom_point(aes(color=USUBJID)) +
  facet_wrap(~USUBJID)

print(g3)

ggsave(filename = paste0(output_dir, "top8_NBV_trajectories.png"), plot = g3, device = "png", width = 10, height = 8)
ggsave(filename = paste0(output_dir, "top8_NBV_trajectories.pdf"), plot = g3, device = "pdf")

#-------------------------------------------------------------------------------------------------
#Run linear mixed effects model on the NBV time series data using lmer  

int_slope_model <- lmer(NBV~DAY + (1 + DAY | USUBJID), data=df)

summary_int_slope <- summary(int_slope_model) #get model summary

# Extract fixed effects i.e., Intercept and slope (DAY) of fixed effects
random_effects <- ranef(int_slope_model)

fixed_effects <- fixef(int_slope_model)


# Create a new dataframe to store the results
results <- data.frame(
  USUBJID = rownames(random_effects$USUBJID),
  Fixed_Intercept = fixed_effects[1],
  Fixed_Slope = fixed_effects[2],
  Random_Intercept = random_effects$USUBJID[,"(Intercept)"],
  Random_Slope = random_effects$USUBJID[,"DAY"]
)

# Calculate the sums
results <- results %>%
  mutate(
    Total_Intercept = Fixed_Intercept + Random_Intercept,
    Total_Slope = Fixed_Slope + Random_Slope
  )

NBV_results <- results %>% select(USUBJID, Total_Intercept, Total_Slope)

#Run LME model for the cube root of VOLT2
#Read in main dataframe
df <- read.csv("/data/users/rty388/testing/02_results.csv")

#Subset the core phase of the study for each subject 
df <- df %>% filter(CORE == "Yes")

df <- df %>%
  mutate(VOLT2CBRT = VOLT2^(1/3))

df <- df %>% select(USUBJID, DAY, VOLT2CBRT)

#omit rows with NA in specific column of data frame
df <- df[!(is.na(df$VOLT2CBRT)), ]

summary(df)

#------------------------------------------------------------------------------------------------
#Visualise the trajectories of VOLT2CBRT over time for 8 subjects

top_subjects <- df %>%
  count(USUBJID, sort = TRUE) %>%
  head(8) %>%
  pull(USUBJID)

#Subset the dataframe to include only the top 8 subjects
df_top8 <- df %>%
  filter(USUBJID %in% top_subjects)

summary(df_top8)

g3 <- ggplot(data = df_top8, aes(x = DAY, y = VOLT2CBRT, group = USUBJID)) +
  geom_line(aes(color=USUBJID)) +
  geom_point(aes(color=USUBJID)) +
  facet_wrap(~USUBJID)


ggsave(filename = paste0(output_dir, "top8_VOLT2CBRT_trajectories.png"), plot = g3, device = "png", width = 10, height = 8)
ggsave(filename = paste0(output_dir, "top8_VOLT2CBRT_trajectories.pdf"), plot = g3, device = "pdf", width = 11, height = 10)

print(g3)

#-------------------------------------------------------------------------------------------------
#Run linear mixed effects model on the VOLT2CBRT time series data using lmer  

int_slope_model <- lmer(VOLT2CBRT~DAY + (1 + DAY | USUBJID), data=df)

summary_int_slope <- summary(int_slope_model) #get model summary

# Extract fixed effects i.e., Intercept and slope (DAY) of fixed effects
random_effects <- ranef(int_slope_model)

fixed_effects <- fixef(int_slope_model)


# Create a new dataframe to store the results
results <- data.frame(
  USUBJID = rownames(random_effects$USUBJID),
  Fixed_Intercept = fixed_effects[1],
  Fixed_Slope = fixed_effects[2],
  Random_Intercept = random_effects$USUBJID[,"(Intercept)"],
  Random_Slope = random_effects$USUBJID[,"DAY"]
)

# Calculate the sums
results <- results %>%
  mutate(
    Total_Intercept = Fixed_Intercept + Random_Intercept,
    Total_Slope = Fixed_Slope + Random_Slope
  )
VOLT2_results <- results %>% select(USUBJID, Total_Intercept, Total_Slope)

#----------------------------------------------------------------------------------------------------
#Run k-means clustering on the total intercept and slope data for NBV and VOLT2
colnames(NBV_results) <- c("USUBJID", "Total_Intercept_NBV", "Total_Slope_NBV")

colnames(VOLT2_results) <- c("USUBJID", "Total_Intercept_VOLT2", "Total_Slope_VOLT2")

df2 <- NBV_results %>% merge(VOLT2_results, by = "USUBJID")

write.csv(df2, paste0(output_dir, "/NBV_VOLT2_Intercepts_Slopes.csv"))

#Elbow Method for finding the optimal number of clusters
# Function to compute total within-cluster sum of squares (wss)
wss <- function(k) {
  kmeans(df2[, c("Total_Intercept_NBV", "Total_Slope_NBV", "Total_Intercept_VOLT2", "Total_Slope_VOLT2")], centers = k, nstart = 10)$tot.withinss
}

# Compute wss for k = 1 to k = 15
k_values <- 1:15
wss_values <- sapply(k_values, wss)

# Create a dataframe for plotting
wss_df <- data.frame(k = k_values, wss = wss_values)

# Plot the elbow plot
Elbow_NBV_VOLT2 <- ggplot(wss_df, aes(x = k, y = wss)) +
  geom_line() +
  geom_point() +
  labs(title = "Elbow Plot for Determining Optimal Number of Clusters",
       x = "Number of Clusters (k)",
       y = "Total Within-Cluster Sum of Squares (WSS)")

print(Elbow_NBV_VOLT2)
ggsave(filename = paste0(output_dir, "Elbow_plot_kmeansNBV_VOLT2.png"), plot = Elbow_NBV_VOLT2 , device = "png", width = 10, height = 8)
ggsave(filename = paste0(output_dir, "Elbow_plot_kmeansNBV_VOLT2.pdf"), plot = Elbow_NBV_VOLT2 , device = "pdf", width = 11, height = 10)


# Calculate the second derivative 
wss_diff1 <- diff(wss_values)
wss_diff2 <- diff(wss_diff1)

# Find the elbow point where the second derivative is maximized (i.e., turning point)
# We add 1 to the index to account for the reduction in length due to the diff function
elbow_point <- which.max(wss_diff2) + 1
print(elbow_point)

# Run k-means clustering on combined NBV and VOLT2 trajectories!!!
set.seed(123) # For reproducibility

km <- kmeans(df2[, c("Total_Intercept_NBV", "Total_Slope_NBV", "Total_Intercept_VOLT2", "Total_Slope_VOLT2")], centers = elbow_point, nstart = 25)

#view results
km

##Plot the clusters 

#Specify the d2 dataframe again as the combination of NBV and VOLT2 results
df2 <- NBV_results %>% merge(VOLT2_results, by = "USUBJID")

# Add the cluster assignment to the d2 dataframe
df2$Cluster <- as.factor(km$cluster)

# Plot the clusters
# Select the columns to plot and add the Cluster column
df2_for_plot <- df2 %>%
  select(Total_Intercept_NBV, Total_Slope_NBV, Total_Intercept_VOLT2, Total_Slope_VOLT2, Cluster)

# Create the scatterplot matrix
scatterplot_matrix <- ggpairs(df2_for_plot, aes(color = Cluster, alpha = 0.6),
                              upper = list(continuous = wrap("points", size = 2)),
                              lower = list(continuous = wrap("points", size = 2)),
                              diag = list(continuous = wrap("densityDiag")))

# Customize the plot
scatterplot_matrix <- scatterplot_matrix + 
  labs(title = "Scatterplot Matrix of Clustering Results") +
  theme_minimal()

print(scatterplot_matrix)
ggsave(filename = paste0(output_dir, "Scatterplot_Matrix_kmeansNBV_VOLT2.png"), plot = scatterplot_matrix, device = "png", width = 10, height = 8)
ggsave(filename = paste0(output_dir, "Scatterplot_Matrix_kmeansNBV_VOLT2.pdf"), plot = scatterplot_matrix, device = "pdf", width = 11, height = 10)

# Compute the principal components
pca <- prcomp(df2[, c("Total_Intercept_NBV", "Total_Slope_NBV", "Total_Intercept_VOLT2", "Total_Slope_VOLT2")], scale. = TRUE)
df2$pca1 <- pca$x[, 1]
df2$pca2 <- pca$x[, 2]

# Plot the clusters
Cluster_Plot <- ggplot(df2, aes(x = pca1, y = pca2, color = Cluster)) +
  geom_point(size = 3) +
  labs(title = "K-means Clustering Results",
       x = "Principal Component 1",
       y = "Principal Component 2") +
  theme_minimal()

print(Cluster_Plot)
ggsave(filename = paste0(output_dir, "Cluster_Plot_kmeansNBV_VOLT2.png"), plot = Cluster_Plot, device = "png", width = 10, height = 8)
ggsave(filename = paste0(output_dir, "Cluster_Plot_kmeansNBV_VOLT2.pdf"), plot = Cluster_Plot, device = "pdf", width = 11, height = 10)

##Calculate the number of patients in each cluster-----------------------------------------------------------------------------
df_cluster1 <- df2 %>%
  filter(Cluster == "1") %>%
  summarise(number_of_subjects = n_distinct(USUBJID))

print(df_cluster1$number_of_subjects)

df_cluster2 <- df2 %>%
  filter(Cluster == "2") %>%
  summarise(number_of_subjects = n_distinct(USUBJID))

print(df_cluster2$number_of_subjects)


