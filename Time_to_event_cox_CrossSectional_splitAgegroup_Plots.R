##Include both young (>=20 and <40 AND old >=50 patients at DAY1)
##Include normal covariates but adjust for age_group instead of age
##Make sure to subset patients who were either 20-40 or >=50 AT Day1!!! 
##Also extract just the core phase of the study to test the effect of ACTVTRT on time-to-event
##Survival curves have been split by age_group so there are two lines for each age group on the same axes


library(knitr)
library(dplyr)
library(survival)
library(ggplot2)
library(tibble)
library(lubridate)
library(ggsurvfit)
library(gtsummary)
library(xlsx)
library(tidyr)
library(survminer)
library(survival)
library(survminer)
library(dplyr)


options(bitmapType = 'cairo')

# Load data
df <- read.csv("/data/users/rty388/testing/02_results.csv")


#Subset the core phase of the study for each subject 
df <- df %>% filter(CORE == "Yes")

#Filter out any subjects who have NA values for their Day1 MRI metrics: 
# subjects_to_exclude <- df %>%
#   filter(time1 == 1 & (is.na(NBV) | is.na(VOLT2cuberoot) | is.na(NUMGDT1squareroot))) %>%
#   pull(USUBJID)
# 
# # Filter out entire subjects
# filtered_df <- df %>%
#   filter(!(USUBJID %in% subjects_to_exclude))
# 
# # Display the result
# print(filtered_df)

# Create NUMGDT1^1/2 and VOLT2^1/3 columns
df <- df %>%
  mutate(NUMGDT1squareroot = NUMGDT1^(1/2),
         VOLT2cuberoot = VOLT2^(1/3))

# Specify event types
event_types <- c("T3MPIRA", "T3MRAW", "T3MCDW")

# Filter dataframe to only include patients aged between 20-40 and >= 50 at DAY1!!!!
USUB_ageDAY1 <- df %>% filter((DAY = 1) & ((AGE >= 20 & AGE < 40) | (AGE >= 50))) %>% distinct(USUBJID)
df <- df %>% filter(USUBJID %in% USUB_ageDAY1$USUBJID)
  
# First, create a separate dataframe for AGE at DAY == 1
df_day1 <- df %>%
  filter(DAY == 1) %>%
  select(USUBJID, AGE) %>%
  rename(AGE_DAY1 = AGE)

# Join this dataframe back to the original dataframe
df <- df %>%
  left_join(df_day1, by = "USUBJID")

# Assign AGE_GROUP based on AGE at DAY == 1
df <- df %>%
  mutate(AGE_GROUP = ifelse(AGE_DAY1 < 40, "20-40", ifelse(AGE_DAY1 >= 50, ">=50", NA)))

# Remove the temporary AGE_DAY1 column if no longer needed
df <- df %>%
  select(-AGE_DAY1)

# View the updated dataframe
print(df)

# Create a dataframe with the values for additional columns specifically for DAY == 1
day1_data <- df %>%
  filter(DAY == 1) %>%
  select(USUBJID, SEX, AGE, AGE_GROUP, DURFS, RELPST1Y, RELPST2Y, ACTVTRT, T25FWM, HPT9M, NBV, NUMGDT1squareroot, VOLT2cuberoot, EDSS) %>%
  distinct()

# Function to perform Cox regression and save results for a given event
run_cox_regression <- function(event_type) {

  # Convert any NA values to "No" in the event column
  df[[event_type]][is.na(df[[event_type]])] <- "No"
  
  df1 <- df %>%
    select(USUBJID, DAY, !!sym(event_type)) %>%
    group_by(USUBJID) %>%
    summarize(
      time = if_else(any(!!sym(event_type) == "TRUE", na.rm = TRUE), as.double(min(DAY[!!sym(event_type) == "TRUE"], na.rm = TRUE)), as.double(max(DAY, na.rm = TRUE))),
      status = if_else(any(!!sym(event_type) == "TRUE", na.rm = TRUE), 1, 0)
    ) %>%
    left_join(day1_data, by = "USUBJID")
  
  # Fit separate Cox models for each age group
  cox_model_20_40 <- coxph(Surv(time, status) ~ SEX + AGE + DURFS + RELPST1Y + RELPST2Y + ACTVTRT + T25FWM + HPT9M + NBV + NUMGDT1squareroot + VOLT2cuberoot + EDSS, data = df1 %>% filter(AGE_GROUP == "20-40"))
  cox_model_50 <- coxph(Surv(time, status) ~ SEX + AGE + DURFS + RELPST1Y + RELPST2Y + ACTVTRT + T25FWM + HPT9M + NBV + NUMGDT1squareroot + VOLT2cuberoot + EDSS, data = df1 %>% filter(AGE_GROUP == ">=50"))
  
  # Generate survival curves
  surv_fit_20_40 <- survfit(cox_model_20_40)
  surv_fit_50 <- survfit(cox_model_50)
  
  # Plot individual survival curves
  plot_20_40 <- ggsurvplot(surv_fit_20_40, data = df1 %>% filter(AGE_GROUP == "20-40"),
                           conf.int = TRUE,
                           ggtheme = theme_bw(),
                           xlim = c(0, 2000),
                           ylim = c(0, 1),
                           xlab = "Days",
                           ylab = paste0("Probability of survival without ", event_type),
                           palette = "Dark2",
                           legend.title = "Age Group",
                           legend.labs = c("20-40"))
  
  plot_50 <- ggsurvplot(surv_fit_50, data = df1 %>% filter(AGE_GROUP == ">=50"),
                        conf.int = TRUE,
                        ggtheme = theme_bw(),
                        xlim = c(0, 2000),
                        ylim = c(0, 1),
                        xlab = "Days",
                        ylab = paste0("Probability of survival without ", event_type),
                        palette = "Dark2",
                        legend.title = "Age Group",
                        legend.labs = c(">=50"))
  
  # Extract survival curve data for each age group
  surv_fit_20_40_df <- data.frame(time = surv_fit_20_40$time,
                                  surv = surv_fit_20_40$surv,
                                  lower = surv_fit_20_40$lower,
                                  upper = surv_fit_20_40$upper,
                                  AGE_GROUP = "20-40")
  
  surv_fit_50_df <- data.frame(time = surv_fit_50$time,
                               surv = surv_fit_50$surv,
                               lower = surv_fit_50$lower,
                               upper = surv_fit_50$upper,
                               AGE_GROUP = ">=50")
  
  # Combine the data frames
  combined_surv_df <- rbind(surv_fit_20_40_df, surv_fit_50_df)
  
  # Plot the combined survival curves
  combined_plot <- ggplot(combined_surv_df, aes(x = time, y = surv, color = AGE_GROUP, fill = AGE_GROUP)) +
    geom_line() +
    geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2) +
    scale_color_manual(values = c(">=50" = "red", "20-40" = "blue"), name = "Age Group") +
    scale_fill_manual(values = c(">=50" = "red", "20-40" = "blue"), name = "Age Group") +
    theme_bw() +
    xlim(0, 2000) +
    xlab("Days") +
    ylab(paste0("Probability of survival without ", event_type)) +
    theme(
      plot.background = element_rect(fill = "white", color = NA),  # Remove black border
      legend.title = element_text(size = 14),
      legend.text = element_text(size = 12),
      axis.title.x = element_text(size = 16, color = "black"),  # Set axis title color to black
      axis.title.y = element_text(size = 16, color = "black"),  # Set axis title color to black
      axis.text.x = element_text(size = 14, color = "black"),  # Set axis text color to black
      axis.text.y = element_text(size = 14, color = "black")   # Set axis text color to black
    )
  

# Display the plot
print(combined_plot)

#Specify output directory 
output_dir <- paste0("/data/users/rty388/regression/Cox/Worsening/CrossSec/round6/")

ggsave(filename = paste0(output_dir, event_type, "/AGE_group/20-40_and_50plus/survfit_COX_curve_", event_type, "_20to40_and_plus50_split2.png"), plot = combined_plot, width = 7, height = 5)
ggsave(filename = paste0(output_dir, event_type, "/AGE_group/20-40_and_50plus/survfit_COX_curve_", event_type, "_20to40_and_plus50_split2.pdf"), plot = combined_plot, width = 7, height = 5)

# Risk curve
risk_20_40_df <- data.frame(
  time = surv_fit_20_40$time,
  surv = surv_fit_20_40$surv,
  lower = surv_fit_20_40$lower,
  upper = surv_fit_20_40$upper, 
  AGE_GROUP = "20-40"
) %>%
  mutate(
    risk = 1 - surv,
    risk_lower = 1 - upper,
    risk_upper = 1 - lower, 
    )

risk_50_df <- data.frame(
  time = surv_fit_50$time,
  surv = surv_fit_50$surv,
  lower = surv_fit_50$lower,
  upper = surv_fit_50$upper, 
  AGE_GROUP = ">=50"
) %>%
  mutate(
    risk = 1 - surv,
    risk_lower = 1 - upper,
    risk_upper = 1 - lower
  )

#Plot combined risk plots with either blue or red figures: 
# Combine the data frames
combined_risk_df <- rbind(risk_20_40_df, risk_50_df)

# Plot the combined survival curves
combined_plot_risk <- ggplot(combined_risk_df, aes(x = time, y = risk, color = AGE_GROUP, fill = AGE_GROUP)) +
  geom_line() +
  geom_ribbon(aes(ymin = risk_lower, ymax = risk_upper), alpha = 0.2) +
  scale_color_manual(values = c(">=50" = "red", "20-40" = "blue"), name = "Age Group") +
  scale_fill_manual(values = c(">=50" = "red", "20-40" = "blue"), name = "Age Group") +
  theme_bw() +
  xlim(0, 2000) +
  xlab("Days") +
  ylab(paste0("Probability of ", event_type)) +
  theme(
    plot.background = element_rect(color = "black", linewidth = 1),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10),
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 14),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12)
  )

combined_plot_risk

combined_plot_risk

# Display the plot
print(combined_plot)

#Specify output directory 
output_dir <- paste0("/data/users/rty388/regression/Cox/Worsening/CrossSec/round6/")

ggsave(filename = paste0(output_dir, event_type, "/AGE_group/20-40_and_50plus/risk_curve_", event_type, "_20to40_and_plus50_split2.png"), plot = combined_plot_risk, width = 7, height = 5)

}

# Loop through the events and run the Cox regression
for (event in event_types) {
  run_cox_regression(event)
}