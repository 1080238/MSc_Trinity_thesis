##Create some plots of Age bins and event frequency from the subsetted dataframe
library(ggplot2)  
library(grid)

df <- read.csv("/data/users/rty388/testing/02_results.csv")

##Explore Age groups on DAY 1 (i.e., First Dosing Date - FDD) ------------------------------------------------------------

DAY1_df <- df[which(df$DAY == "1"), ]

DAY1_df$AGE <- as.numeric(DAY1_df$AGE)

# Adjusting the Age bins to every 10 years
DAY1_df$AgeGroup <- cut(DAY1_df$AGE,
                        breaks = seq(0, 100, by = 10),  # Creating bins from 0 to 100 by 10s
                        labels = paste(seq(0, 90, by = 10), seq(9, 99, by = 10), sep = "-"),
                        right = FALSE,  # 'right = FALSE' excludes the upper boundary
                        include.lowest = TRUE)  # Includes the lowest value in the first bin


Age_bins <- table(DAY1_df$AgeGroup)
print(Age_bins)

write.csv(Age_bins, "/data/users/rty388/testing/Covariates/Age_bins.csv", row.names = F)

options(bitmapType = 'cairo')

Age_sex_freq <- table(DAY1_df$SEX, DAY1_df$AgeGroup)
write.csv(Age_sex_freq, "/data/users/rty388/testing/Covariates/Age_sex_freq.csv", row.names = F)

# Create a barplot with 10-year Age bins
# Calculate counts for each AgeGroup
DAY1_counts <- DAY1_df %>%
  group_by(AgeGroup) %>%
  summarise(Count = n())

#DAY1_df$Counts <- DAY1_counts$Count

# Create the plot with annotations
png("/data/users/rty388/testing/Covariates/Age_bins.png", height = 1300, width = 1600)
ggplot(DAY1_counts, aes(x = AgeGroup, y = Count)) +
  geom_bar(stat = "identity", fill = "steelblue", color = "black") +
  geom_text(aes(label = Count), vjust = -0.5, size = 14) +  # Add text annotations
  labs(title = "Frequency of patients by Age Group at Day1",
       x = "Age Group",
       y = "Number of Patients") +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 50),  # Title font size
    axis.title.x = element_text(size = 60, vjust = -0.5),  # X-axis title font size
    axis.title.y = element_text(size = 60, vjust = 1.5),  # Y-axis title font size
    axis.text.x = element_text(size = 50, color = "black"),  # X-axis text font size
    axis.text.y = element_text(size = 50, color = "black"),  # Y-axis text font size
    plot.margin = unit(c(1, 1, 1.5, 1), "cm")  # Slightly increase bottom margin
  )
dev.off()

Age_plot <- ggplot(DAY1_counts, aes(x = AgeGroup, y = Count)) +
  geom_bar(stat = "identity", fill = "steelblue", color = "black") +
  geom_text(aes(label = Count), vjust = -0.5, size = 14) +  # Add text annotations
  labs(title = "Frequency of patients by Age Group at Day1",
       x = "Age Group",
       y = "Number of Patients") +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 50),  # Title font size
    axis.title.x = element_text(size = 60, vjust = -0.5),  # X-axis title font size
    axis.title.y = element_text(size = 60, vjust = 1.5),  # Y-axis title font size
    axis.text.x = element_text(size = 50, color = "black"),  # X-axis text font size
    axis.text.y = element_text(size = 50, color = "black"),  # Y-axis text font size
    plot.margin = unit(c(1, 1, 1.5, 1), "cm")  # Slightly increase bottom margin
  )
ggsave("/data/users/rty388/testing/Covariates/Age_bins2.png", plot = Age_plot, width = 20, height = 17)
ggsave("/data/users/rty388/testing/Covariates/Age_bins2.pdf", plot = Age_plot, width = 20, height = 17)
