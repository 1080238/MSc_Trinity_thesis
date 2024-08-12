#######################################################################################
###
### Project: NO.MS - Data for MSc Students
### Project members: 
###         S. Pendleton (Oxford BDI, samantha.pendleton@ndm.ox.ac.uk)
###         S. Gardiner (Oxford BDI, stephen.gardiner@ndm.ox.ac.uk)
###         H. Ganjgahi (Oxford BDI, habib.ganjgahi@bdi.ox.ac.uk)
###         P. Aarden (NVS, piet.aarden@novartis.com)
### Derived: Based on the HMM work
### Application: R 4.1.0
### Date: 03-05-2024
### Teresa: 15-05-2024
#######################################################################################

##Note: FDD = First Dosing Date

#Load libraries
library(tidyverse)
library(haven)
library(reshape2)
library(lubridate)

scripts <- file.path("/data/users/rty388/Scripts/")

source(paste0(scripts, "variables.R"))
source(paste0(scripts, "wrangling_functions.R"))

# Load of data ------------------------------------------------------------

if(location == "BDI") { source_saspath <- file.path("/data/ms/unprocessed/clinical/NOMS_Version2_20211022-OXF-ANALYTICS/")
} else if (location == "NVS") {
  user <- system("whoami", intern = TRUE)
  source_saspath <- paste0("/view/", user, "_view/vob/CNOMSANON/anon/anon_2/anon_source/")
}

# Data
base <- read_sas(paste0(source_saspath, "base.sas7bdat")) %>% as.data.frame()
edss <- read_sas(paste0(source_saspath, "edss.sas7bdat")) %>% as.data.frame()
mri <- read_sas(paste0(source_saspath, "mri.sas7bdat")) %>% as.data.frame()
follow <- read_sas(paste0(source_saspath, "follow.sas7bdat")) %>% as.data.frame()

msfc <- read_sas(paste0(source_saspath, "msfc.sas7bdat")) %>% as.data.frame()
sdmt <- read_sas(paste0(source_saspath, "sdmt.sas7bdat")) %>% as.data.frame()
nfl <- read_sas(paste0(source_saspath, "nfl.sas7bdat")) %>% as.data.frame()

cdwevents <- read_sas(paste0(source_saspath, "cdwevents.sas7bdat")) %>% as.data.frame()

# Load in previously wrangled data
#This 00_results.csv was created from the 00_HMM_relapse_wrangling.R script 
relapse <- read_csv("/data/users/rty388/testing/00_results.csv", 
                    col_types = cols('c','c','D','D','i','i','f','D','i','f','f','f','f')) %>% as.data.frame()

# Concatenating to one baseline and attach follow-up visits ---------------

# EDSS --------------------------------------------------------------------

# edss_mod: table with EDSS data
# - Take only one baseline observation (before FDD, DAY == 1)
# - Add remaining observations
# - Remove multiple observations on one date (keeping worst EDSS)
# - Only non-missing EDSS observations are kept
edss_mod <- edss %>%
  filter(DAY <= 1, !is.na(EDSS), STUDY %in% studies) %>%
  select(c(USUBJID, STUDY, DATE, DAY, EDSS, VISFNC, BRNFNC, BOWFNC, CLRFNC, CRBFNC, PYRFNC, SENFNC)) %>%
  arrange(USUBJID, desc(DAY)) %>%
  filter(!duplicated(.$USUBJID)) %>%
  bind_rows(., edss[which((edss$DAY > 1) & (!is.na(edss$EDSS)) & (edss$STUDY %in% studies)), 
                    c("USUBJID", "STUDY", "DATE", "DAY", "EDSS", "VISFNC", "BRNFNC", "BOWFNC", "CLRFNC", "CRBFNC", "PYRFNC", "SENFNC")]) %>%
  arrange(USUBJID, DATE, desc(EDSS)) %>%
  filter(!duplicated(.[c("USUBJID", "DATE")])) %>%
  mutate(DAY = ifelse(DAY <= 1, 1, DAY)) ##If it is less than or equal to 1, call it DAY 1

# df_month: Table with the month assigned (addition to df_wtrt)
# - Using Dieter Haering's visit window function we get the interval for each month
# - Assign the monthly visit using data.table foverlap
# - Assign baseline month (DAY == 1) to -1
mco_monthly <- vw(seq(0, 175, 1))$dat
edss_mod <- data.table::data.table(edss_mod)
edss_mod[, DAYTMP := DAY]
mcotmp <- data.table::data.table(mco_monthly[c("Month", "upper", "lower", "target")])
data.table::setkey(mcotmp, lower, upper)

edss_month <- data.table::foverlaps(edss_mod, mcotmp, by.x = c('DAY', 'DAYTMP'), 
                                    by.y = c('lower', 'upper')) %>% 
  as.data.frame() %>%
  rename(MONTH = Month) %>%
  mutate(MONTH = ifelse(DAY == 1, -1, MONTH)) %>%
  select(c(USUBJID, STUDY, DATE, DAY, target, MONTH, EDSS, VISFNC, BRNFNC, BOWFNC, CLRFNC, CRBFNC, PYRFNC, SENFNC))

# - Only one observation per month is kept, eliminate multiple observations
#   per month by keeping the once closest to the target date
edss_month <- edss_month %>%
  mutate(DIFF = abs(target-DAY)) %>%
  arrange(USUBJID, MONTH, DIFF) %>%
  filter(!duplicated(.[c("USUBJID", "MONTH")])) %>%
  select(-c(target, DIFF))

#rm(edss, edss_mod, mco_monthly, mcotmp)

# MSFC --------------------------------------------------------------------

# t25fwm_mod: table with timed 25-foot walk data
# - Take only one baseline observation (before FDD, DAY == 1)
# - Add remaining observations
# - Remove multiple observations on one date (keeping worst T25FWM)
# - Only non-missing T25FWM observations are kept
# - rename DATE column for later
t25fw_mod <- msfc %>%
  select(c(USUBJID, STUDY, DATE, DAY, T25FWM)) %>%
  filter(!is.na(DAY), DAY <= 1, !is.na(T25FWM), STUDY %in% studies) %>%
  arrange(USUBJID, desc(DAY)) %>%
  filter(!duplicated(.$USUBJID)) %>%
  bind_rows(., msfc[which((msfc$DAY > 1) & (!is.na(msfc$T25FWM)) & (msfc$STUDY %in% studies) & (!is.na(msfc$DAY))), 
                    c("USUBJID", "STUDY", "DATE", "DAY", "T25FWM")]) %>%
  arrange(USUBJID, DATE, desc(T25FWM)) %>%
  filter(!duplicated(.[c("USUBJID", "DATE")])) %>%
  rename(T25DAY = DAY) %>%
  select(c(USUBJID, T25DAY, T25FWM)) %>%
  mutate(T25DAY = ifelse(T25DAY <= 1, 1, T25DAY))

# hpt9_mod: table with 9-hole peg test data
# - Take only one baseline observation (before FDD, DAY == 1)
# - Add remaining observations
# - Remove multiple observations on one date (keeping worst HPT9M)
# - Only non-missing HPT9M observations are kept
# - rename DATE column for later
hpt9_mod <- msfc %>%
  select(c(USUBJID, STUDY, DATE, DAY, HPT9M)) %>%
  filter(!is.na(DAY), DAY <= 1, !is.na(HPT9M), STUDY %in% studies) %>%
  arrange(USUBJID, desc(DAY)) %>%
  filter(!duplicated(.$USUBJID)) %>%
  bind_rows(., msfc[which((msfc$DAY > 1) & (!is.na(msfc$HPT9M)) & (msfc$STUDY %in% studies) & (!is.na(msfc$DAY))), 
                    c("USUBJID", "STUDY", "DATE", "DAY", "HPT9M")]) %>%
  arrange(USUBJID, DATE, desc(HPT9M)) %>%
  filter(!duplicated(.[c("USUBJID", "DATE")])) %>%
  rename(HPTDAY = DAY) %>%
  select(c(USUBJID, HPTDAY, HPT9M)) %>%
  mutate(HPTDAY = ifelse(HPTDAY <= 1, 1, HPTDAY))

# pasat_mod: table with PASAT data
# - Take only one baseline observation (before FDD, DAY == 1)
# - Add remaining observations
# - Remove multiple observations on one date (keeping worst PASAT)
# - Only non-missing PASAT observations are kept
# - rename DATE column for later
pasat_mod <- msfc %>%
  select(c(USUBJID, STUDY, DATE, DAY, PASAT)) %>%
  filter(!is.na(DAY), DAY <= 1, !is.na(PASAT), STUDY %in% studies) %>%
  arrange(USUBJID, desc(DAY)) %>%
  filter(!duplicated(.$USUBJID)) %>%
  bind_rows(., msfc[which((msfc$DAY > 1) & (!is.na(msfc$PASAT)) & (msfc$STUDY %in% studies) & (!is.na(msfc$DAY))), 
                    c("USUBJID", "STUDY", "DATE", "DAY", "PASAT")]) %>%
  arrange(USUBJID, DATE, PASAT) %>%
  filter(!duplicated(.[c("USUBJID", "DATE")])) %>%
  rename(PASATDAY = DAY) %>%
  select(c(USUBJID, PASATDAY, PASAT)) %>%
  mutate(PASATDAY = ifelse(PASATDAY <= 1, 1, PASATDAY))

# sdmt_mod: table with SDMT data
# - Take only one baseline observation (before FDD, DAY == 1)
# - Add remaining observations
# - Remove multiple observations on one date (keeping worst SDMT)
# - Only non-missing SDMT observations are kept
# - rename DATE column for later
sdmt_mod <- sdmt %>%
  select(c(USUBJID, STUDY, DATE, DAY, SDMT)) %>%
  filter(!is.na(DAY), DAY <= 1, !is.na(SDMT), STUDY %in% studies) %>%
  arrange(USUBJID, desc(DAY)) %>%
  filter(!duplicated(.$USUBJID)) %>%
  bind_rows(., sdmt[which((sdmt$DAY > 1) & (!is.na(sdmt$SDMT)) & (sdmt$STUDY %in% studies) & (!is.na(sdmt$DAY))), 
                    c("USUBJID", "STUDY", "DATE", "DAY", "SDMT")]) %>%
  arrange(USUBJID, DATE, SDMT) %>%
  filter(!duplicated(.[c("USUBJID", "DATE")])) %>%
  rename(SDMTDAY = DAY) %>%
  select(c(USUBJID, SDMTDAY, SDMT)) %>%
  mutate(SDMTDAY = ifelse(SDMTDAY <= 1, 1, SDMTDAY))

# rm(msfc, sdmt)

# MRI ---------------------------------------------------------------------

# mri_mod: table mri data
# - Take only one baseline observation (before FDD, DAY == 1)
# - Add remaining observations
# - Remove multiple observations on one date
# - rename DATE column for later
mri_mod <- mri %>%
  dplyr::select(c(USUBJID, STUDY, DATE, DAY, NBV, PERBVC, NUMGDT1, VOLT2)) %>%
  filter(DAY <= 1, STUDY %in% studies, !is.na(DAY)) %>%
  arrange(USUBJID, desc(DAY)) %>% #arrange by descending so the most negative day will be at the end and the day closest to the FDD is at the top
  filter(!duplicated(.$USUBJID)) %>% #Thus, only the day closest to the start of the trial (FDD) will be kept
  bind_rows(., mri[which((mri$DAY > 1) & (mri$USUBJID != "") & (mri$STUDY %in% studies)), 
                   c("USUBJID", "STUDY", "DATE", "DAY", "NBV", "PERBVC", "NUMGDT1", "VOLT2")]) %>%
  arrange(USUBJID, DATE) %>% 
  filter(!duplicated(.[c("USUBJID", "DATE")])) %>%
  rename(MRIDT = DATE, MRIDAY = DAY) %>%
  select(c(USUBJID, MRIDT, MRIDAY, NBV, PERBVC, NUMGDT1, VOLT2)) %>%
  mutate(MRIDAY = ifelse(MRIDAY <= 1, 1, MRIDAY))

# rm(mri)

# NFL ---------------------------------------------------------------------

# nfl_mod: table nfl data
# - Take only one baseline observation (before FDD, DAY == 1)
# - Add remaining observations
# - Remove multiple observations on one date
# - rename DATE column for later
nfl_mod <- nfl %>%
  dplyr::select(c(USUBJID, STUDY, DATE, DAY, NFL)) %>%
  filter(DAY <= 1, STUDY %in% studies, !is.na(DAY)) %>%
  arrange(USUBJID, desc(DAY)) %>%
  filter(!duplicated(.$USUBJID)) %>%
  bind_rows(., nfl[which((nfl$DAY > 1) & (nfl$USUBJID != "") & (nfl$STUDY %in% studies)), 
                   c("USUBJID", "STUDY", "DATE", "DAY", "NFL")]) %>%
  arrange(USUBJID, DATE) %>% 
  filter(!duplicated(.[c("USUBJID", "DATE")])) %>%
  rename(NFLDT = DATE, NFLDAY = DAY) %>%
  select(c(USUBJID, NFLDT, NFLDAY, NFL)) %>%
  mutate(NFLDAY = ifelse(NFLDAY <= 1, 1, NFLDAY))

# rm(nfl)

# CDW ---------------------------------------------------------------------

## cdwevents_mod: table containing the 6-month confirmed PIRA events
## - Can choose to subset cdwevents to 6 month confirmed PIRA events
# cdwevents_mod <- cdwevents %>%
#   filter(DISEVENT == "T6MCDW", PIRAFLG == "Yes", STUDY %in% studies)

# rm(cdwevents)

# Mapping -----------------------------------------------------------------

# t25fw_edss_map: table containg mapping of t25fw data to edss observations
# - Aligns all observations in t25fw_mod to a unique edss observation in edss_month
# - Maps observations to closest edss observations with a maximum absolute difference of 30 days
# - Makes sure that mapped edss is only used once
t25fw_edss_map <- t25fw_mod %>%
  select(c(USUBJID, T25DAY)) %>% 
  merge(., edss_month[c("USUBJID", "DAY")], by = "USUBJID", all.x = T) %>%
  mutate(DIFF = abs(T25DAY - DAY)) %>%
  filter(DIFF <= 30) %>%
  arrange(USUBJID, T25DAY, DIFF) %>%
  filter(!duplicated(.[c("USUBJID", "T25DAY")])) %>%
  arrange(USUBJID, DAY, DIFF) %>%
  filter(!duplicated(.[c("USUBJID", "DAY")])) %>%
  arrange(USUBJID, T25DAY) %>%
  select(-DIFF)

# hpt9_edss_map: table containg mapping of hpt9 data to edss observations
# - Aligns all observations in hpt9_mod to a unique edss observation in edss_month
# - Maps observations to closest edss observations with a maximum absolute difference of 30 days
# - Makes sure that mapped edss is only used once
hpt9_edss_map <- hpt9_mod %>%
  select(c(USUBJID, HPTDAY)) %>% 
  merge(., edss_month[c("USUBJID", "DAY")], by = "USUBJID", all.x = T) %>%
  mutate(DIFF = abs(HPTDAY - DAY)) %>%
  filter(DIFF <= 30) %>%
  arrange(USUBJID, HPTDAY, DIFF) %>%
  filter(!duplicated(.[c("USUBJID", "HPTDAY")])) %>%
  arrange(USUBJID, DAY, DIFF) %>%
  filter(!duplicated(.[c("USUBJID", "DAY")])) %>%
  arrange(USUBJID, HPTDAY) %>%
  select(-DIFF)

# pasat_edss_map: table containg mapping of pasat data to edss observations
# - Aligns all observations in pasat_mod to a unique edss observation in edss_month
# - Maps observations to closest edss observations with a maximum absolute difference of 30 days
# - Makes sure that mapped edss is only used once
pasat_edss_map <- pasat_mod %>%
  select(c(USUBJID, PASATDAY)) %>% 
  merge(., edss_month[c("USUBJID", "DAY")], by = "USUBJID", all.x = T) %>%
  mutate(DIFF = abs(PASATDAY - DAY)) %>%
  filter(DIFF <= 30) %>%
  arrange(USUBJID, PASATDAY, DIFF) %>%
  filter(!duplicated(.[c("USUBJID", "PASATDAY")])) %>%
  arrange(USUBJID, DAY, DIFF) %>%
  filter(!duplicated(.[c("USUBJID", "DAY")])) %>%
  arrange(USUBJID, PASATDAY) %>%
  select(-DIFF)

# sdmt_edss_map: table containg mapping of sdmt data to edss observations
# - Aligns all observations in sdmt_mod to a unique edss observation in edss_month
# - Maps observations to closest edss observations with a maximum absolute difference of 30 days
# - Makes sure that mapped edss is only used once
sdmt_edss_map <- sdmt_mod %>%
  select(c(USUBJID, SDMTDAY)) %>% 
  merge(., edss_month[c("USUBJID", "DAY")], by = "USUBJID", all.x = T) %>%
  mutate(DIFF = abs(SDMTDAY - DAY)) %>%
  filter(DIFF <= 30) %>%
  arrange(USUBJID, SDMTDAY, DIFF) %>%
  filter(!duplicated(.[c("USUBJID", "SDMTDAY")])) %>%
  arrange(USUBJID, DAY, DIFF) %>%
  filter(!duplicated(.[c("USUBJID", "DAY")])) %>%
  arrange(USUBJID, SDMTDAY) %>%
  select(-DIFF)

# mri_edss_map: table containg mapping of MRI data to edss observations
# - Aligns all observations in mri_mod to a unique edss observation in edss_month
# - Maps observations to closest edss observations with a maximum absolute difference of 90 days
# - Makes sure that mapped edss is only used once
mri_edss_map <- mri_mod %>% 
  select(c(USUBJID, MRIDAY)) %>% 
  merge(., edss_month[c("USUBJID", "DAY")], by = "USUBJID", all.x = T) %>%
  mutate(DIFF = abs(MRIDAY - DAY)) %>%
  filter(DIFF <= 90) %>%
  arrange(USUBJID, MRIDAY, DIFF) %>%
  filter(!duplicated(.[c("USUBJID", "MRIDAY")])) %>%
  arrange(USUBJID, DAY, DIFF) %>%
  filter(!duplicated(.[c("USUBJID", "DAY")])) %>%
  arrange(USUBJID, MRIDAY) %>%
  select(-DIFF)

# nfl_edss_map: table containg mapping of NFL data to edss observations
# - Aligns all observations in nfl_mod to a unique edss observation in edss_month
# - Maps observations to closest edss observations with a maximum absolute difference or XXXXXXXXXXXXX days
# - Makes sure that mapped edss is only used once
nfl_edss_map <- nfl_mod %>% 
  select(c(USUBJID, NFLDAY)) %>% 
  merge(., edss_month[c("USUBJID", "DAY")], by = "USUBJID", all.x = T) %>%
  mutate(DIFF = abs(NFLDAY - DAY)) %>%
  filter(DIFF <= 90) %>%
  arrange(USUBJID, NFLDAY, DIFF) %>%
  filter(!duplicated(.[c("USUBJID", "NFLDAY")])) %>%
  arrange(USUBJID, DAY, DIFF) %>%
  filter(!duplicated(.[c("USUBJID", "DAY")])) %>%
  arrange(USUBJID, NFLDAY) %>%
  select(-DIFF)

# Merging -----------------------------------------------------------------

# df: Table containing all endpoints and a flag for a PIRA onset
# - Use previous *_mod tables and *_edss_map tables to merge everything
# - Merge the PIRA flag on top of the observations using cdwevents_mod
# - Keep MRIDT as its used later for merging addition MRI data
##TEST keeping PIRA and DISEVENT data in df

df <- edss_month %>%
  merge(., t25fw_edss_map, by = c("USUBJID", "DAY"), all.x = T) %>%
  merge(., hpt9_edss_map, by = c("USUBJID", "DAY"), all.x = T) %>%
  merge(., pasat_edss_map, by = c("USUBJID", "DAY"), all.x = T) %>%
  merge(., sdmt_edss_map, by = c("USUBJID", "DAY"), all.x = T) %>%
  merge(., mri_edss_map, by = c("USUBJID", "DAY"), all.x = T) %>%
  merge(., nfl_edss_map, by = c("USUBJID", "DAY"), all.x = T) %>%
  
  merge(., t25fw_mod, by.x = c("USUBJID", "T25DAY"), all.x = T) %>%
  merge(., hpt9_mod, by = c("USUBJID", "HPTDAY"), all.x = T) %>%
  merge(., pasat_mod, by = c("USUBJID", "PASATDAY"), all.x = T) %>%
  merge(., sdmt_mod, by = c("USUBJID", "SDMTDAY"), all.x = T) %>%
  merge(., mri_mod, by = c("USUBJID", "MRIDAY"), all.x = T) %>%
  merge(., nfl_mod, by = c("USUBJID", "NFLDAY"), all.x = T) %>% 
  
  select(USUBJID, STUDY, MRIDT, DATE, DAY, MONTH, EDSS, T25FWM, HPT9M, PASAT, SDMT, NBV, PERBVC, NUMGDT1, VOLT2, NFL) %>%
  
  merge(., cdwevents[c("USUBJID", "ONSTDAY", "PIRAFLG", "RAWFLG", "DISEVENT")], by.x = c("USUBJID", "DAY"), by.y = c("USUBJID", "ONSTDAY"), all.x = T) %>%
  
  arrange(USUBJID, DAY)


# rm(cdwevents_mod, edss_month, 
#    hpt9_edss_map, hpt9_mod, mri_edss_map, mri_mod, nfl_edss_map, nfl_mod,
#    pasat_edss_map, pasat_mod, sdmt_edss_map, sdmt_mod, t25fw_edss_map, t25fw_mod)

# Relapses ----------------------------------------------------------------

# df_wrel: Table with the relapse info for each observation (addition to df)
# - use AddRelapseInfo to check whether observation was measured during a relapse
df_wrel <- AddRelapseInfo(df, relapse)

# Treatment ---------------------------------------------------------------

# df_wtrt: Table containing the treatment at time of observation (addition to df_wrelrcv)
# - use AddTRTInfo to check what treatment patient was on during observation
df_wtrt <- AddTRTInfo(df_wrel, follow)

# rm(relapse, follow)

# Final steps -------------------------------------------------------------

# df_final: final table which includes demographic data
# - merge demographic data on table
# - transform ARM to ACTVTRT (active treatment 0/1)
# - offset AGE and DURFS longitudinally (AGELG, DURFSLG respectively)
# - Impute a longitudinal MS phenotype (MSTYPELG) by the condition of someone having
#   a PIRA event and a EDSS >= 3 (captured in df_spmsonset constructed below)
df_spmsonset <- df_wtrt %>%
  filter(PIRAFLG == 1, EDSS >= 3) %>%
  arrange(USUBJID, MONTH) %>%
  filter(!duplicated(.$USUBJID)) %>%
  rename(PIRAMONTH = MONTH) %>%
  select(c(USUBJID, PIRAMONTH))


# Merge several dataframes including the follow dataframe to include ENDDATE of the CORE study
follow_unique <- follow %>% distinct(USUBJID, .keep_all = TRUE)

df_final <- df_wtrt %>%
  merge(., base[c("USUBJID", "SITE", "SEX", "AGE", "MSTYPE", "BMI", "DURFS", "RELPST1Y", "RELPST2Y")], by = "USUBJID", all.x = TRUE) %>%
  merge(., df_spmsonset, by = c("USUBJID"), all.x = TRUE) %>%
  merge(., follow_unique[c("USUBJID", "ENDDATE")], by = "USUBJID", all.x = TRUE) %>% 
  mutate(ACTVTRT = ifelse(ARM %in% c("No treatment", "Placebo"), 0, 1),
         YEARS = DAY / 365.25,
         BMI = plyr::mapvalues(BMI, c("< 18.5", "[18.5,25)", "[25,30)", ">= 30"), c("18", "21.75", "27.5", "30")),
         BMI = as.numeric(BMI),
         DURFS = plyr::mapvalues(DURFS, c("[0,2)", "[2,5)", "[5,10)", "[10,30)", "[30,50)"), c("1", "3.5", "7.5", "20", "40")),
         DURFS = as.numeric(DURFS),
         DURFSLG = DURFS + YEARS,
         AGELG = floor(AGE + YEARS),
         MSTYPELG = case_when(
           MSTYPE == "SPMS" ~ "SPMS",
           MSTYPE == "PPMS" ~ "PPMS",
           (MSTYPE == "RRMS") & (is.na(PIRAMONTH)) ~ "RRMS",
           (MSTYPE == "RRMS") & (!is.na(PIRAMONTH)) & (MONTH < PIRAMONTH) ~ "RRMS",
           (MSTYPE == "RRMS") & (!is.na(PIRAMONTH)) & (MONTH >= PIRAMONTH) ~ "SPMS"
         ),
         CORE = ifelse(DATE < ENDDATE, "Yes", "No")  # Create CORE column
  ) %>%
  arrange(USUBJID, MONTH) %>%
  select(USUBJID, STUDY, DATE, CORE,
         MONTH, DAY, YEARS, MRIDT,
         SEX,
         AGE, AGELG,
         MSTYPE, MSTYPELG,
         DURFS, DURFSLG,
         RELPST1Y, RELPST2Y,
         RELAPSE, ACTVTRT,
         EDSS,
         PIRAFLG, RAWFLG, DISEVENT, 
         T25FWM, HPT9M, PASAT, SDMT,
         NBV, PERBVC, NUMGDT1, VOLT2, 
         NFL, ARM
  )

#Calculate updated NBV
df_final <- df_final %>%
  group_by(USUBJID) %>%
  mutate(BaselineNBV = ifelse(DAY == 1, NBV, NA)) %>%
  fill(BaselineNBV, .direction = "downup") %>%
  mutate(NBV = ifelse(DAY == 1, NBV, BaselineNBV * (1 + PERBVC / 100))) %>%
  ungroup() %>%
  select(-BaselineNBV)


# rm(base,df,df_spmsonset,df_wrel,df_wtrt)

library(dplyr)
library(stringr)

# Create new columns for confirmed PIRA and RAW at 3 months and 6 months ---------------------------------------------------
df_final2 <- df_final %>%
  mutate(T3MPIRA = PIRAFLG == "Yes" & DISEVENT == "T3MCDW") %>%
  mutate(T6MPIRA = PIRAFLG == "Yes" & DISEVENT == "T6MCDW") %>%
  mutate(T3MRAW = RAWFLG == "Yes" & DISEVENT == "T3MCDW") %>%
  mutate(T6MRAW = RAWFLG == "Yes" & DISEVENT == "T6MCDW") %>%
  mutate(T3MCDW = DISEVENT == "T3MCDW") %>%
  mutate(T6MCDW = DISEVENT == "T6MCDW")

# filling confirmed NAs ---------------------------------------------------
df_final3 <- df_final2 %>%
  mutate(Event = case_when(
    T3MPIRA ~ "T3MPIRA",
    T6MPIRA ~ "T6MPIRA",
    T3MRAW  ~ "T3MRAW",
    T6MRAW  ~ "T6MRAW",
    TRUE    ~ NA_character_  # default if none of the above is TRUE
  )) 

# output ------------------------------------------------------------------
# Save table
write.csv(df_final3, "/data/users/rty388/testing/02_results.csv", row.names = F)

#######################

# End of script
