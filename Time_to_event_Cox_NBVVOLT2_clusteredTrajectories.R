### Time-to-event analysis - Loop through different Ages and different Events
##Look at effect of the combined NBV and VOLT2CBRT when grouped into 2 clusters using k-means
# COX REGRESSION ---------------------------------------------------------------

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

options(bitmapType = 'cairo')

df <- read.csv("/data/users/rty388/testing/02_results.csv")

#Read in clusters for kmeans on combined NBV and VOLT2CBRT (k = 2 clusters)
NBV_VOLT2_clusters <- read.csv("/data/users/rty388/testing/Cluster_NBV_VOLT2/2kmeansclusters/kmeans2Clusters_NBV_VOLT2.csv")

NBV_VOLT2_clusters <- NBV_VOLT2_clusters %>% select(USUBJID, Cluster)

clusters <- levels(as.factor(NBV_VOLT2_clusters$Cluster))

# # Iterate over each unique cluster value and replace with the desired format
 for (i in clusters) {
   NBV_VOLT2_clusters$Cluster <- gsub(paste0("\\b", i, "\\b"), paste0("cluster", i), NBV_VOLT2_clusters$Cluster)
 }

# View the modified dataframe
print(NBV_VOLT2_clusters)

colnames(NBV_VOLT2_clusters)[2] <- "NBV_VOLT2_clusters"

#Subset the core phase of the study for each subject 
df <- df %>% filter(CORE == "Yes")

# Create NUMGDT1^1/2 and VOLT2^1/3 columns
df <- df %>%
  mutate(NUMGDT1squareroot = NUMGDT1^(1/2),
         VOLT2cuberoot = VOLT2^(1/3))

df <- df %>% merge(NBV_VOLT2_clusters, by = "USUBJID")

#Calculate the number of patients in the old and young age group in each cluster
#df_test <- df %>% filter(DAY == 1 & (AGE >= 20 & AGE <40) | AGE >= 50) %>%  mutate(AGE_GROUP = ifelse(AGE >= 20 & AGE < 40), "20-40", ">=50")  
#length(which(df_test$NBV_VOLT2_clusters == "cluster1"))
#length(which(df_test$NBV_VOLT2_clusters == "cluster2"))

#Specify age groups and event types
age_groups <- list(
  "greaterthan50" = list(filter = function(age) age >= 50, label = "AGE >= 50"),
  "greaterthan40" = list(filter = function(age) age >= 40, label = "AGE >= 40"),
  "40to50" = list(filter = function(age) age >= 40 & age < 50, label = "40 <= AGE < 50"),
  "20to40" = list(filter = function(age) age >= 20 & age < 40, label = "20 <= AGE < 40"),
  "greaterthan20" = list(filter = function(age) age >= 20, label = "AGE >= 20")
)

#Specify event types
event_types <- list(
  "T3MPIRA" = list(flag = "T3MPIRA", label = "T3MPIRA"),
  "T3MRAW" = list(flag = "T3MRAW", label = "T3MRAW"),
  "T3MCDW" = list(flag = "T3MCDW", label = "T3MCDW")
)

for (age_group_name in names(age_groups)) {
  
  age_group <- age_groups[[age_group_name]]
  
  for (event_type_name in names(event_types)) {
    event_type <- event_types[[event_type_name]]
    
    print(event_type)
    print(age_group)
    
    DAY1_df <- df %>% filter(DAY == 1)
    DAY1_df$AGE <- as.numeric(DAY1_df$AGE)
    
    #Filter the data by age group at Day1
    USUB_greater_DAY1 <- subset(DAY1_df, age_group$filter(AGE))$USUBJID
    df_subset <- df %>% filter(USUBJID %in% USUB_greater_DAY1)
    
    #Change the NA values to No
    df_subset$T3MPIRA[is.na(df_subset$T3MPIRAFLG)] <- "No"
    df_subset$T3MRAW[is.na(df_subset$T3MRAW)] <- "No"
    df_subset$T3MCDW[is.na(df_subset$T3MCDW)] <- "No"
    
    # Create a dataframe with values for additional columns specifically for DAY == 1
    day1_data <- df_subset %>%
      filter(DAY == 1) %>%
      select(USUBJID, SEX, AGE, DURFS, RELPST1Y, RELPST2Y, ACTVTRT, T25FWM, HPT9M, NBV_VOLT2_clusters, NUMGDT1squareroot, EDSS) %>%
      distinct()
    
    # Determine time and status for the current event type
    df1 <- df_subset %>%
      select(USUBJID, DAY, !!sym(event_type$flag)) %>%
      group_by(USUBJID) %>%
      summarize(
        time = if_else(any(!!sym(event_type$flag) == "TRUE", na.rm = TRUE), as.double(min(DAY[!!sym(event_type$flag) == "TRUE"], na.rm = TRUE)), as.double(max(DAY, na.rm = TRUE))),
        status = if_else(any(!!sym(event_type$flag) == "TRUE", na.rm = TRUE), 1, 0)
      ) %>%
      left_join(day1_data, by = "USUBJID")
    
    
    
    # Fit the Cox proportional hazards model
    cox_model <- coxph(Surv(time, status) ~ SEX + AGE + DURFS + RELPST1Y + RELPST2Y + ACTVTRT + T25FWM + HPT9M + NBV_VOLT2_clusters + NUMGDT1squareroot + EDSS, data = df1)
    summary_res.cox <- summary(cox_model)
    
    #Test if model meets the proportional hazards assumption: should be p>0.05 for each covariate and global
    ph_test <- cox.zph(cox_model)
    print(ph_test)
    
    
    # Extract the coefficients table
    coefficients <- summary_res.cox$coefficients
    coefficients_df <- as.data.frame(coefficients)
    
    # Save the coefficients dataframe to a CSV file
    output_dir <- paste0("/data/users/rty388/regression/Cox/Worsening/CrossSec/round10/", event_type_name, "/", age_group_name)
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
    
    xlsx::write.xlsx(coefficients_df, paste0(output_dir, "/COX_model_summary_", event_type_name, "_", age_group_name, ".xls"), col.names=TRUE, row.names=FALSE, sheetName="sample")
    write.csv(coefficients_df, paste0(output_dir, "/COX_model_summary_", event_type_name, "_", age_group_name, ".csv"), row.names = T)
    
    # Table of results for COX regression
    COX_table <- tbl_regression(cox_model, exp = TRUE)
    COX_table_df <- as_tibble(COX_table)
    
    
    ##Add on mean and SD for each covariate to COX table 
    covariates <- c("SEX", "AGE", "DURFS", "RELPST1Y", "RELPST2Y", "ACTVTRT", 
                    "T25FWM", "HPT9M", "NBV_VOLT2_clusters", "NUMGDT1squareroot", "EDSS")
    
    
    # Calculate mean and SD, ensuring proper handling of non-numeric covariates
    
    df1$NUMGDT1squareroot <- as.numeric(df1$NUMGDT1squareroot)
    
    df1 <- df1[!is.na(df1$NUMGDT1squareroot), ]
    
    
    mean_sd_df <- df1 %>%
      select(all_of(covariates)) %>%
      summarise(across(everything(), list(mean = mean, sd = sd), na.rm = TRUE)) %>%
      pivot_longer(everything(), 
                   names_to = c(".value", "name"), 
                   names_sep = "_")
    
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
    
    # Merging the dataframes
    final_table <- COX_table_df %>%
      left_join(df_new2, by = "Characteristic")
    
    # Print the final table
    print(final_table)
    
    xlsx::write.xlsx(final_table, paste0(output_dir, "/table_cox_regression_", event_type_name, "_", age_group_name, ".xls"), col.names=TRUE, row.names=TRUE, sheetName="sample")
    write.csv(final_table, paste0(output_dir, "/table_cox_regression_", event_type_name, "_", age_group_name, ".csv"), row.names = FALSE)
    
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
      xlsx::write.xlsx(final_table2, paste0(output_dir, "/sigtable_cox_regression_", event_type_name, "_", age_group_name, ".xls"), col.names = TRUE, row.names = TRUE, sheetName = "sample")
      write.csv(final_table2, paste0(output_dir, "/sig_table_cox_regression_", event_type_name, "_", age_group_name, ".csv"), row.names = FALSE)
    } else {
      message("No significant results to write for event: ", event_type_name, " and age group: ", age_group_name)
    }
    
    # Get the survival fit
    surv_fit <- survfit(cox_model)
    
    # Extract data from survival fit to save to csv
    surv_data <- data.frame(
      time = surv_fit$time,
      n.risk = surv_fit$n.risk,
      n.event = surv_fit$n.event,
      n.censor = surv_fit$n.censor,
      surv = surv_fit$surv,
      std.err = surv_fit$std.err,
      upper = surv_fit$upper,
      lower = surv_fit$lower
    )
    xlsx::write.xlsx(surv_data, paste0(output_dir, "/survfit_Cox_", event_type_name, "_", age_group_name, ".xls"), col.names=TRUE, row.names=TRUE, sheetName="sample")
    write.csv(surv_data, paste0(output_dir, "/survfit_Cox_", event_type_name, "_", age_group_name, ".csv"), row.names = FALSE)
    
    # Survival curve
    surv_plot <- survfit2(cox_model) %>%
      ggsurvfit(color = "blue") +
      labs(
        x = "Days",
        y = paste("Overall probability of not having", event_type_name)
      ) +
      add_confidence_interval(color = "blue", fill = "lightblue") +
      add_risktable(
        size = 4,
        theme = theme_risktable_default(axis.text.y.size = 12, plot.title.size = 14)
      ) +
      coord_cartesian(xlim = c(0, 2000), clip = "off")
    
    ggsave(filename = paste0(output_dir, "/survfit_COX_curve_", event_type_name, "_", age_group_name, ".png"), plot = surv_plot, device = "png", width = 10, height = 8)
    
    # Risk curve
    surv_df <- data.frame(
      time = surv_fit$time,
      surv = surv_fit$surv,
      lower = surv_fit$lower,
      upper = surv_fit$upper
    ) %>%
      mutate(
        risk = 1 - surv,
        risk_lower = 1 - upper,
        risk_upper = 1 - lower
      )
    
    risk_plot <- ggplot(surv_df, aes(x = time, y = risk)) +
      geom_step(color = "red") +
      geom_ribbon(aes(ymin = risk_lower, ymax = risk_upper), fill = "red", alpha = 0.2) +
      labs(
        x = "Days",
        y = paste("Overall probability of having", event_type_name)
      ) +
      theme_minimal() +
      ggtitle(paste("Risk Curve: Probability of Having", event_type_name)) +
      theme(legend.position = "none") +
      coord_cartesian(xlim = c(0, 2000), clip = "off")
    
    ggsave(filename = paste0(output_dir, "/risk_curve_COX_", event_type_name, "_", age_group_name, ".png"), plot = risk_plot, device = "png", width = 10, height = 8)
  }
}
