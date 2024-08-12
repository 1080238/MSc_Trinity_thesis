library(ggplot2)
library(reshape2)


# Assuming summary_table is already created from the previous steps
DAY1_greater50_summary_table <- read.csv("/data/users/rty388/testing/Freq_summary_for_all_events_DAY1greater50.csv")

DAY1_20to40_summary_table <- read.csv("/data/users/rty388/testing/Freq_summary_for_all_events_DAY1_20to40.csv")

# Add a column to distinguish between the two datasets
DAY1_greater50_summary_table$Group <- "Greater than 50"
DAY1_20to40_summary_table$Group <- "20 to 40"

# Combine the dataframes
combined_data <- bind_rows(DAY1_greater50_summary_table, DAY1_20to40_summary_table)

combined_data$Group <- factor(combined_data$Group, levels = c("Greater than 50", "20 to 40"))

# Reshape the data for plotting
percentage_data <- combined_data %>%
  select(Event, Group, Percentage_event) %>%
  gather(key = "Metric", value = "Value", -Event, -Group)

percentage_data$Event <- factor(percentage_data$Event, levels = c("CDW", "T3MCDW", "T6MCDW", "PIRAFLG", "T3MPIRA", "T6MPIRA", "RAWFLG", "T3MRAW", "T6MRAW"))

percentage_data$Value <- round(percentage_data$Value, 1)

percentage_data2 <- percentage_data %>% filter(Event %in% c("T3MCDW", "T3MPIRA", "T3MRAW"))

#png("/data/users/rty388/testing/Events/Percent_mainevents_splitAge_20to40_greater50.png", height = 1000, width = 1500)
Main_event_plot <- ggplot(percentage_data2, aes(x = Event, y = Value, fill = Group)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_text(aes(label = Value), vjust = -0.5, size = 16, position = position_dodge(width = 0.9)) +  # Add text annotations
  scale_fill_manual(values = c("Greater than 50" = "#006500", "20 to 40" = "#CC5500"),
                    labels = c("Greater than 50", "20 to 40")) +
  labs(title = "Percentage of Patients Experiencing Each Event",
       x = "Event Type",
       y = "Percentage of Patients",
       fill = "Age Group") +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 35, color = "black"),  # Title font size
    axis.title.x = element_text(size = 60,  color = "black"),  # X-axis title font size
    axis.title.y = element_text(size = 60, color = "black"),  # Y-axis title font size
    axis.text.x = element_text(size = 60, color = "black", angle = 45, hjust = 1),  # X-axis text font size and angle
    axis.text.y = element_text(size = 60, color = "black"),  # Y-axis text font size
    legend.title = element_text(size = 60, color = "black"),  # Legend title font size
    legend.text = element_text(size = 60, color = "black")  # Legend text font size
  )

ggsave("/data/users/rty388/testing/Events/Percent_mainevents_splitAge_20to40_greater50.png", plot = Main_event_plot, width = 30, height = 25)
ggsave("/data/users/rty388/testing/Events/Percent_mainevents_splitAge_20to40_greater50.pdf", plot = Main_event_plot, width = 30, height = 25)

