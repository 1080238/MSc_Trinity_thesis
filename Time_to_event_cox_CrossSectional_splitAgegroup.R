##Run Cox regression with AGE_GROUP as a covariate
##Split graph by AGE_GROUP 
##Make sure to subset patients who were either 20-40 or >=50 AT Day1!!! 
##Also extract just the core phase of the study to test the effect of ACTVTRT on time-to-event
##Survival curve has a single line of prediction for all patients (both young and old)

library(knitr)
library(dplyr)
library(survival)
library(ggplot2)
library(tibble)
library(lubridate)
library(ggsurvfit)
library(gtsummary)
library(survival)
library(survminer)

options(bitmapType = 'cairo')

df <- read.csv("/data/users/rty388/testing/02_results.csv")

#Subset the core phase of the study for each subject 
df <- df %>% filter(CORE == "Yes")

# Create NUMGDT1^1/2 and VOLT2^1/3 columns
df <- df %>%
  mutate(NUMGDT1squareroot = NUMGDT1^(1/2),
         VOLT2cuberoot = VOLT2^(1/3))

##Filter data to only include patients who were >=20 & <40 OR >=50 at DAY1
DAY1_USUBJID <- df %>% 
  filter(DAY == 1 & ((AGE >= 20 & AGE < 40) | AGE >= 50)) %>% 
  select(USUBJID)

df_subset <- df %>% 
  filter(USUBJID %in% DAY1_USUBJID$USUBJID) %>%
  mutate(AGE_GROUP = ifelse(AGE < 40, "20-40", ">=50"))


# Create a dataframe with the select covariates for DAY == 1
day1_data <-  df_subset %>%
  filter(DAY == 1) %>%
  select(USUBJID, SEX, AGE, AGE_GROUP, DURFS, RELPST1Y, RELPST2Y, ACTVTRT, T25FWM, HPT9M, NBV, NUMGDT1squareroot, VOLT2cuberoot, EDSS) %>%
  distinct()

# Specify event types
event_types <- c("T3MPIRA", "T3MRAW", "T3MCDW")

#Find the day of first PIRAFLG for patients with at least one case of PIRA
#In patients who never had PIRAFLG take longest date recorded
# Create the summarized dataframe and join with day1_data

run_cox_regression <- function(event_type) {

#Convert any NA values to "No" in the event column
df_subset[[event_type]][is.na( df_subset[[event_type]])] <- "No"
  
df1 <-  df_subset %>%
    select(USUBJID, DAY, !!sym(event_type)) %>%
    group_by(USUBJID) %>%
    summarize(
      time = if_else(any(!!sym(event_type) == "TRUE", na.rm = TRUE), as.double(min(DAY[!!sym(event_type) == "TRUE"], na.rm = TRUE)), as.double(max(DAY, na.rm = TRUE))),
      status = if_else(any(!!sym(event_type) == "TRUE", na.rm = TRUE), 1, 0)
    ) %>%
    left_join(day1_data, by = "USUBJID")

# Fit the Cox model 
cox_model <- coxph(Surv(time, status) ~ SEX + AGE_GROUP + DURFS + RELPST1Y + RELPST2Y + ACTVTRT + T25FWM + HPT9M + NBV + NUMGDT1squareroot + VOLT2cuberoot + EDSS, data = df1)

summary_res.cox <- summary(cox_model)

coefficients_df <- as.data.frame(summary_res.cox$coefficients)

#Specify output directory 
output_dir <- paste0("/data/users/rty388/regression/Cox/Worsening/CrossSec/round6/")

#Write coefficients_df file to csv and xls
xlsx::write.xlsx(coefficients_df, paste0(output_dir, event_type, "/AGE_group/20-40_and_50plus/COX_model_summary_", event_type, "_20to40_and_plus50.xls"), col.names=TRUE, row.names=FALSE, sheetName="sample")
write.csv(coefficients_df, paste0(output_dir, event_type, "/AGE_group/20-40_and_50plus/COX_model_summary_", event_type, "_20to40_and_plus50.csv"), row.names = TRUE)

COX_table <- tbl_regression(cox_model, exp = TRUE)
COX_table_df <- as_tibble(COX_table)

covariates <- c("SEX", "AGE_GROUP", "DURFS", "RELPST1Y", "RELPST2Y", "ACTVTRT", "T25FWM", "HPT9M", "NBV", "NUMGDT1squareroot", "VOLT2cuberoot", "EDSS")

mean_sd_df <- df1 %>%
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

# Fit the Cox model with stratification by AGE_GROUP
cox_model_strata <- coxph(Surv(time, status) ~ SEX + DURFS + RELPST1Y + RELPST2Y + ACTVTRT + T25FWM + HPT9M + NBV + NUMGDT1squareroot + VOLT2cuberoot + strata(AGE_GROUP), data = df1)

# Create survival curves stratified by AGE_GROUP
surv_fit_strata <- survfit(cox_model_strata)

df1 <- df1 %>%
  mutate(AGE_GROUP = factor(AGE_GROUP, levels = c("20-40", ">=50")))

# Plot the survival curves using ggsurvplot from survminer

surv_plot <- ggsurvplot(surv_fit_strata, data = df1,
           conf.int = TRUE,
           risk.table = TRUE,
           risk.table.height = 0.25,
           ggtheme = theme_bw(),
           xlim = c(0, 2000),
           xlab = "Days",
           ylab = paste0("Overall probability of not having", event_type),
           palette = "Dark2",
           legend.title = "Age Group",
           legend.labs = c("20-40", "50"))

surv_plot_gg <- surv_plot$plot

ggsave(filename = paste0(output_dir, event_type, "/AGE_group/20-40_and_50plus/survfit_COX_curve_", event_type, "_20to40_and_plus50_split.png"), plot = surv_plot_gg)

} 

# Loop through the events and run the Cox regression
for (event in event_types) {
  run_cox_regression(event)
}
