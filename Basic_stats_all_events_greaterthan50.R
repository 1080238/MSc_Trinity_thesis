##Plot the summary stats of the different covariates 
#Run the plots for the => 40 years olds and the 20-40 years old

##SEX, AGE, DURFS, RELPST1Y, RELPST2Y, ACTVTRT, T25FWM, HPT9M, NBV, NUMGDT1, VOLT2, EDSS
##EVENTS: "CDW" T3MCDW", "T6MCDW", "PIRAFLG", "RAWFLG", "T3MRAW", "T6MRAW", "T3MPIRA", "T6MPIRA"

library(ggplot2)
library(reshape2)
library(dplyr)


df <- read.csv("/data/users/rty388/testing/02_results.csv")

# Create CDW column
df <- df %>%
  mutate(
    CDW = ifelse(PIRAFLG == "Yes" | RAWFLG == "Yes", "TRUE", NA)
  )

#Subset all rows and columns where the DAY1 AGE of the patients was equal to or greater than 40 
DAY1_df <- df[which(df$DAY == "1"), ]

DAY1_greater50_USUBJID <- DAY1_df$USUBJID[which(DAY1_df$AGE >= 50)]

DAY1_greater50_df <- df[df$USUBJID %in% DAY1_greater50_USUBJID, ]

# List of event columns
event_columns <- c("CDW", "T3MCDW", "T6MCDW", "PIRAFLG", "RAWFLG", "T3MRAW", "T6MRAW", "T3MPIRA", "T6MPIRA")

# Replace TRUE with "Yes" and FALSE with "No" for the event columns
DAY1_greater50_df <- DAY1_greater50_df %>%
  mutate(across(all_of(event_columns), ~ ifelse(. == TRUE, "Yes", ifelse(. == FALSE, "No", .))))

# Initialize an empty list to store summaries
summaries <- list()

#Test: 
#event <- event_columns[1]

# Loop through each event type
for (event in event_columns) {
  # Summarize whether each subject experienced the event at least once
  DAY1_greater50_df_event <- DAY1_greater50_df %>%
    group_by(USUBJID) %>%
    summarize(!!event := ifelse(any(.data[[event]] == "Yes"), "Yes", "No")) %>%
    ungroup()
  
  # Calculate total number of patients
  total_num_patients <- nrow(DAY1_greater50_df_event)
  
  # Calculate number of patients who experienced the event at least once
  num_patients_event <- length(which(DAY1_greater50_df_event[[event]] == "Yes"))
  
  # Calculate percentage of patients who experienced the event at least once
  percentage_event <- (num_patients_event / total_num_patients) * 100
  
  # Calculate number of patients who did not experience the event even once
  num_patients_no_event <- total_num_patients - num_patients_event
  
  # Calculate percentage of patients who did not experience the event even once
  percentage_no_event <- (num_patients_no_event / total_num_patients) * 100
  
  # Create a summary for the event
  event_summary <- data.frame(
    Event = event,
    Total_num_patients = total_num_patients,
    Num_patients_event = num_patients_event,
    Percentage_event = percentage_event,
    Num_patients_no_event = num_patients_no_event,
    Percentage_no_event = percentage_no_event
  )
  
  # Append the summary to the list
  summaries[[event]] <- event_summary
}

# Combine all summaries into a single dataframe
summary_table <- do.call(rbind, summaries)

# Display the summary table
print(summary_table)

# Write the summary table to a CSV file
write.csv(summary_table, "/data/users/rty388/testing/Freq_summary_for_all_events_DAY1greater50.csv", row.names = FALSE)


