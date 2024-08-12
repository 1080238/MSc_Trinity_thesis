##Longitudinal time to CDW - split by AGEGROUP

library(knitr)
library(dplyr)
library(survival)
library(ggplot2)
library(tibble)
library(lubridate)
library(ggsurvfit)
library(gtsummary)
library(tidyr)
library(xlsx)
library(survminer)

options(bitmapType = 'cairo')

#Read in main dataframe
df <- read.csv("/data/users/rty388/testing/longitudinal-table-cdw.csv")


#Subset the core phase of the study for each subject >>Core phase 
# df <- df %>% filter(CORE == "Yes")

colnames(df) <- gsub("T3MCDW", "status", colnames(df))

# Filter dataframe to only include patients aged between 20-40 and >= 50 at DAY1!!!!
USUB_ageDAY1 <- df %>% filter((time1 = 1) & ((AGELG >= 20 & AGELG < 40) | (AGELG >= 50))) %>% distinct(USUBJID)
df <- df %>% filter(USUBJID %in% USUB_ageDAY1$USUBJID)

# To create age-group column: First, create a separate dataframe for AGE at DAY == 1
df_day1 <- df %>%
  filter(time1 == 1) %>%
  select(USUBJID, AGELG) %>%
  rename(AGE_DAY1 = AGELG)

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

#Test out filtering data to only include up to two years (i.e., day 730.5): 
# df <- df %>% filter(time1 <= 730.5)


# Specify event types
event_types <- "T3MCDW"

event_type <- "T3MCDW"

# run_cox_regression <- function(event_type) {
  
  # Fit separate Cox models for each age group
  cox_model_20_40 <- coxph(Surv(time1,time2,status) ~ SEX + AGELG + DURFSLG + RELPST1Y + RELPST2Y + ACTVTRT + T25FWM + HPT9M + NBV + NUMGDT1SQRT + VOLT2CBRT + EDSSBL, data = df %>% filter(AGE_GROUP == "20-40"))
  cox_model_50 <- coxph(Surv(time1,time2,status) ~ SEX + AGELG + DURFSLG + RELPST1Y + RELPST2Y + ACTVTRT + T25FWM + HPT9M + NBV + NUMGDT1SQRT + VOLT2CBRT + EDSSBL, data = df %>% filter(AGE_GROUP == ">=50"))
  
  # Generate survival curves
  surv_fit_20_40 <- survfit(cox_model_20_40)
  surv_fit_50 <- survfit(cox_model_50)
  
  # Plot individual survival curves
  plot_20_40 <- ggsurvplot(surv_fit_20_40, data = df %>% filter(AGE_GROUP == "20-40"),
                           conf.int = TRUE,
                           ggtheme = theme_bw(),
                           xlim = c(0, 750),
                           ylim = c(0, 1),
                           xlab = "Days",
                           ylab = paste0("Probability of survival without ", event_type),
                           palette = "Dark2",
                           legend.title = "Age Group",
                           legend.labs = c("20-40"))
  
  plot_50 <- ggsurvplot(surv_fit_50, data = df %>% filter(AGE_GROUP == ">=50"),
                        conf.int = TRUE,
                        ggtheme = theme_bw(),
                        xlim = c(0, 750),
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
    xlim(0, 750) +
    xlab("Days") +
    ylab(paste0("Probability of survival without", event_type)) +
    theme(
      plot.background = element_rect(fill = "white", color = NA),  # Remove black border
      legend.title = element_text(size = 14),
      legend.text = element_text(size = 12),
      axis.title.x = element_text(size = 16, color = "black"),  # Set axis title color to black
      axis.title.y = element_text(size = 16, color = "black"),  # Set axis title color to black
      axis.text.x = element_text(size = 14, color = "black"),  # Set axis text color to black
      axis.text.y = element_text(size = 14, color = "black")   # Set axis text color to black
    )
  dev.off()
  
  # Display the plot
  print(combined_plot)
  
  #Specify output directory 
  # Define the variables
  output_dir <- "/data/users/rty388/regression/Cox/Worsening/Longitudinal/round1/"
  event_type <- paste0(event_type)
  
  # Create the filename
  filename <- paste0(output_dir, event_type, "/AGE_group/20-40_and_50plus/survfit_COX_curve_", event_type, "_20to40_and_plus50_split2.png")
  ggsave(filename, plot = combined_plot, width = 7, height = 6)
 
  filename <- paste0(output_dir, event_type, "/AGE_group/20-40_and_50plus/survfit_COX_curve_", event_type, "_20to40_and_plus50_split2.pdf")
  ggsave(filename, plot = combined_plot, width = 7, height = 5)
  
  
  # Create data frames with cumulative probability of CDI
  risk_fit_20_40_df <- data.frame(
    time = surv_fit_20_40$time,
    cum_prob_cdi = 1 - surv_fit_20_40$surv,
    lower = 1 - surv_fit_20_40$upper,
    upper = 1 - surv_fit_20_40$lower,
    AGE_GROUP = "20-40"
  )
  
  risk_fit_50_df <- data.frame(
    time = surv_fit_50$time,
    cum_prob_cdi = 1 - surv_fit_50$surv,
    lower = 1 - surv_fit_50$upper,
    upper = 1 - surv_fit_50$lower,
    AGE_GROUP = ">=50"
  )
  
  # Combine the data frames
  combined_risk_df <- rbind(risk_fit_20_40_df, risk_fit_50_df)
  
  combined_risk_df$AGE_GROUP <- as.factor(combined_risk_df$AGE_GROUP)
  
  # Plot the combined risk curves
  combined_plot <- ggplot(combined_risk_df, aes(x = time, y = cum_prob_cdi, color = AGE_GROUP, fill = AGE_GROUP)) +
    geom_line() +
    geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2) +
    scale_color_manual(values = c(">=50" = "red", "20-40" = "blue"), name = "Age Group") +
    scale_fill_manual(values = c(">=50" = "red", "20-40" = "blue"), name = "Age Group") +
    theme_bw() +
    xlim(0, 750) +
    xlab("Days") +
    ylab(paste0("Cumulative probability of ", event_type))
  
  # Print the plot
  print(combined_plot)
  
  # Create the filename
  filename <- paste0(output_dir, event_type, "/AGE_group/20-40_and_50plus/risk_COX_curve_", event_type, "_20to40_and_plus50_split2.png")
  ggsave(filename, plot = combined_plot, width = 7, height = 5)
  
  
  ##Save cox model results--------------------------------------------------------------------------
  
  # Fit the Cox model 
  cox_model <- coxph(Surv(time1,time2,status) ~ SEX + AGE_GROUP + DURFSLG + RELPST1Y + RELPST2Y + ACTVTRT + T25FWM + HPT9M + NBV + NUMGDT1SQRT + VOLT2CBRT + EDSSBL, data = df)
  
  ph_test <- cox.zph(cox_model)
  print(ph_test)
  
  summary_res.cox <- summary(cox_model)
  
  coefficients_df <- as.data.frame(summary_res.cox$coefficients)
  
  #Write coefficients_df file to csv and xls
  xlsx::write.xlsx(coefficients_df, paste0(output_dir, event_type, "/AGE_group/20-40_and_50plus/COX_model_summary_", event_type, "_20to40_and_plus50.xls"), col.names=TRUE, row.names=FALSE, sheetName="sample")
  write.csv(coefficients_df, paste0(output_dir, event_type, "/AGE_group/20-40_and_50plus/COX_model_summary_", event_type, "_20to40_and_plus50.csv"), row.names = TRUE)
  
  COX_table <- tbl_regression(cox_model, exp = TRUE)
  COX_table_df <- as_tibble(COX_table)
  
  covariates <- c("SEX", "AGE_GROUP", "DURFSLG", "RELPST1Y", "RELPST2Y", "ACTVTRT", "T25FWM", "HPT9M", "NBV", "NUMGDT1SQRT", "VOLT2CBRT", "EDSSBL")
  
  mean_sd_df <- df %>%
    select(all_of(covariates)) %>%
    summarise(across(everything(), list(mean = mean, sd = sd), na.rm = TRUE)) %>%
    pivot_longer(everything(), 
                 names_to = c(".value", "name"), 
                 names_sep = "_")
  
  mean_sd_df <- mean_sd_df[-(3:4), ]
  
  test <- t(mean_sd_df)
  test <- as.data.frame(test)
  test <- test[-1, ] 
  colnames(test) <- c("Mean", "SD")
  
  df_new <- as.data.frame(cbind(rownames(test), test$Mean, test$SD))
  colnames(df_new) <- c("Characteristic", "Mean", "SD")
  
  df_new$Mean <- round(as.numeric(df_new$Mean), 2)
  df_new$SD <- round(as.numeric(df_new$SD), 2)
  
  Mean_SD <- paste0(df_new$Mean, " ± ", df_new$SD)
  
  df_new2 <- cbind(df_new$Characteristic, Mean_SD)
  df_new2 <- as.data.frame(df_new2)
  
  colnames(df_new2) <- c("Characteristic", "Mean ± SD")
  colnames(COX_table_df) <- c("Characteristic", "Hazard Ratio", "95% CI", "p-value")
  
  final_table <- COX_table_df %>%
    left_join(df_new2, by = "Characteristic")
  
  print(final_table)
  
  #Write Cox table to csv and xls
  xlsx::write.xlsx(final_table, paste0(output_dir, event_type, "/AGE_group/20-40_and_50plus/table_cox_regression_", event_type, "_20to40_and_plus50.xls"), col.names=TRUE, row.names=TRUE, sheetName="sample")
  write.csv(final_table, paste0(output_dir, event_type, "/AGE_group/20-40_and_50plus/table_cox_regression_", event_type, "_20to40_and_plus50.csv"), row.names = FALSE)
  
  ## Significant p-values only Cox table
  colnames(final_table)
  
  final_table <- as.data.frame(final_table)
  
  colnames(final_table) <- c("Characteristic", "Hazard Ratio", "95% CI", "pvalue", "Mean ± SD")
  
  clean_pvalue <- function(pvalue) {
    # Remove any non-numeric characters and convert to numeric
    pvalue <- gsub("[^0-9.]", "", pvalue)
    as.numeric(pvalue)
  }
  
  # Apply the function to the pvalue column
  final_table$pvalue <- sapply(final_table$pvalue, clean_pvalue)
  
  final_table2 <- final_table[which(final_table$pvalue <0.09), ]
  
  # Check if final_table2 is not empty before writing to files
  if (nrow(final_table2) > 0) {
    xlsx::write.xlsx(final_table2, paste0(output_dir, event_type, "/AGE_group/20-40_and_50plus/sigtable_cox_regression_", event_type, "_20to40_and_plus50.xls"), col.names=TRUE, row.names=TRUE, sheetName="sample")
    write.csv(final_table2, paste0(output_dir, event_type, "/AGE_group/20-40_and_50plus/sigtable_cox_regression_", event_type, "_20to40_and_plus50.csv"), row.names = FALSE)
  } else {
    message("No significant results to write for event: ", event_type_name, " and age group: ", age_group_name)
  }
  
  
# }

# # Loop through the events and run the Cox regression
# for (event in event_types) {
#   run_cox_regression(event)
# }