########### CODE FOR PHENOTYPE ANALYSIS ########### 
# WITH THIS CODE WE WILL ANALYSE THE RELATIONSHIP BETWEEN ERECTILE DYSFUNCTION AND CKD
########### CONTENTS ########### 
# 1.  LOADING DATA FROM STATA
# 2A. PHENOTYPE ANALYSIS FOR ERECTILE DYSFUNCTION
# 2B. PHENOTYPE ANALYSIS FOR CKD
# 2C. PHENOTYPE ANALYSIS FOR DIABETES MELLITUS
# 2D. PHENOTYPE ANALYSIS FOR HYPERTENSION
# 2E. PHENOTYPE ANALYSIS FOR ISCHAEMIC HEART DISEASE
# 2F. PHENOTYPE ANALYSIS FOR CEREBRAL INFARCTION/TIA
# 3.  PREPARING DATA FOR ANALYSIS
# 4.  TABLES
# 5.  FIRST SURVIVAL ANALYSES
# 6.  IMPUTATION OF MISSING DATA (BMI, TDI, SMOKING STATUS)
# 7.  PROPENSITY SCORE MATCHING
#################################

#################################
# 1.  LOADING DATA FROM STATA
#################################
library(tidyverse)
library(haven)
library(survival)
library(survminer)
library(Hmisc)
library(MatchIt)
library(tableone)
ed_data <- as.data.frame(read_dta("~/pilot/R_ed_ckd_data.dta"))

#select patients with available primary care data 
include <- read.csv("~/pilot/primarycaredata_available.csv", header = T)
ed_data <- ed_data[ed_data$n_eid %in% include$n_eid,]

#################################
# 2A. PHENOTYPE ANALYSIS FOR ERECTILE DYSFUNCTION
#################################
#I now want to extract the self-reported erectile dysfunction diagnoses.
#to do this, I tell R to look for the self-reported diagnostic codes in the 
#many variables that encode this.
#(n_20002_*_*) is the diagnostic code
#(n_20009_*_*) is the corresponding age
sr_codes_ed <- 1518
ed_data$sr_ed_dx <- NA
ed_data$sr_ed_age <- NA
for (i in 0:2){
  for (k in 0:33){
    variable1 = paste0("n_20002_", i, "_", k)
    variable2 = paste0("n_20009_", i, "_", k)
    if (is.null(ed_data[[variable1]])) {
      next } else {
        ed_data$sr_ed_dx <- ifelse(
          ed_data[[variable1]] %in% sr_codes_ed,
          ed_data[[variable1]], ed_data$sr_ed_dx)
        ed_data$sr_ed_age <- ifelse(
          ed_data[[variable1]] %in% sr_codes_ed,
          ed_data[[variable2]] , ed_data$sr_ed_age)
      }
  }
}
#remove impossible values, i.e. negatives etc
ed_data$sr_ed_age <- ifelse(ed_data$sr_ed_age < 16, NA, ed_data$sr_ed_age)

# GP diagnoses for ed.
ed_data$gp_ed_dx <- ifelse(!ed_data$read_2_ed_gp == "", 
                               ed_data$read_2_ed_gp, ifelse(
                                 !ed_data$read_3_ed_gp == "",
                                 ed_data$read_3_ed_gp, NA
                               )
)

ed_data$gp_ed_age <- as.numeric(difftime(ed_data$event_date_ed, ed_data$DOB, units = "days")/365.25)
#remove impossible values, i.e. negatives etc
ed_data$gp_ed_age <- ifelse(ed_data$gp_ed_age < 16, NA, ed_data$gp_ed_age)

# HES diagnoses for erectile dysfunction.
ed_data$hes_ed_dx <- ifelse(!ed_data$diag_icd10_ed == "",
                                ed_data$diag_icd10_ed, NA)
ed_data$hes_ed_age <- as.numeric(difftime(ed_data$date_1st_hes_ed, ed_data$DOB, units = "days")/365.25)

# first occurrence for erectile dysfunction - not available unfortunately.
ed_data$fo_ed_dx <- NA
ed_data$fo_ed_age <- NA

#and now combining them all into a yes/no ed variable.
ed_data$comb_ed_dx <- ifelse(
  !is.na(ed_data$sr_ed_dx) | !is.na(ed_data$gp_ed_dx) | 
    !is.na(ed_data$hes_ed_dx) | !is.na(ed_data$fo_ed_dx) #see above
  ,
  T, F
)
ed_data$comb_ed_age <- NA
ed_data[ed_data$comb_ed_dx == T,] <- ed_data[ed_data$comb_ed_dx == T,] %>% 
  rowwise() %>%
  mutate(comb_ed_age = min(sr_ed_age, hes_ed_age, gp_ed_age,fo_ed_age #see above
                             , na.rm = T))
ed_data$comb_ed_age <- ifelse(
  ed_data$comb_ed_age == Inf, NA, ed_data$comb_ed_age
)

ed_data$comb_ed_dxsource <- ifelse(!is.na(ed_data$hes_ed_dx) & !is.na(ed_data$gp_ed_dx) & !is.na(ed_data$fo_ed_dx) & !is.na(ed_data$sr_ed_dx),
                                       "HES, GP, FO, SR",
                                       ifelse(!is.na(ed_data$hes_ed_dx) & !is.na(ed_data$gp_ed_dx) & !is.na(ed_data$fo_ed_dx),
                                              "HES, GP, FO",
                                              ifelse(!is.na(ed_data$hes_ed_dx) & !is.na(ed_data$gp_ed_dx) & !is.na(ed_data$sr_ed_dx),
                                                     "HES, GP, SR",
                                                     ifelse(!is.na(ed_data$hes_ed_dx) & !is.na(ed_data$fo_ed_dx) & !is.na(ed_data$sr_ed_dx),
                                                            "HES, FO, SR",
                                                            ifelse(!is.na(ed_data$sr_ed_dx) & !is.na(ed_data$gp_ed_dx) & !is.na(ed_data$fo_ed_dx),
                                                                   "GP, FO, SR",
                                                                   ifelse(!is.na(ed_data$hes_ed_dx) & !is.na(ed_data$gp_ed_dx),
                                                                          "HES, GP",
                                                                          ifelse(!is.na(ed_data$hes_ed_dx) & !is.na(ed_data$fo_ed_dx),
                                                                                 "HES, FO",
                                                                                 ifelse(!is.na(ed_data$hes_ed_dx) & !is.na(ed_data$sr_ed_dx),
                                                                                        "HES, SR",
                                                                                        ifelse(!is.na(ed_data$fo_ed_dx) & !is.na(ed_data$gp_ed_dx),
                                                                                               "GP, FO",
                                                                                               ifelse(!is.na(ed_data$sr_ed_dx) & !is.na(ed_data$gp_ed_dx),
                                                                                                      "GP, SR",
                                                                                                      ifelse(!is.na(ed_data$fo_ed_dx) & !is.na(ed_data$sr_ed_dx),
                                                                                                             "FO, SR",
                                                                                                             ifelse(!is.na(ed_data$hes_ed_dx),
                                                                                                                    "HES",
                                                                                                                    ifelse(!is.na(ed_data$gp_ed_dx),
                                                                                                                           "GP",
                                                                                                                           ifelse(!is.na(ed_data$fo_ed_dx),
                                                                                                                                  "FO",
                                                                                                                                  ifelse(!is.na(ed_data$sr_ed_dx),
                                                                                                                                         "SR",
                                                                                                                                         NA
                                                                                                                                  )))))))))))))))
#################################
# 2B. PHENOTYPE ANALYSIS FOR CKD
#################################
#I now want to extract the self-reported ckd diagnoses.
#to do this, I tell R to look for the self-reported diagnostic codes in the 
#many variables that encode this.
#(n_20002_*_*) is the diagnostic code
#(n_20009_*_*) is the corresponding age
sr_codes_ckd <- c(1192, 1193, 1194) #1192 = renal/kidney failure, 1193 = renal failure req dialysis, 1194 = renal failure not requiring dialysis
ed_data$sr_ckd_dx <- NA
ed_data$sr_ckd_age <- NA
for (i in 0:2){
  for (k in 0:33){
    variable1 = paste0("n_20002_", i, "_", k)
    variable2 = paste0("n_20009_", i, "_", k)
    if (is.null(ed_data[[variable1]])) {
      next } else {
        ed_data$sr_ckd_dx <- ifelse(
          ed_data[[variable1]] %in% sr_codes_ckd,
          ed_data[[variable1]], ed_data$sr_ckd_dx)
        ed_data$sr_ckd_age <- ifelse(
          ed_data[[variable1]] %in% sr_codes_ckd,
          ed_data[[variable2]] , ed_data$sr_ckd_age)
      }
  }
}
#remove impossible values, i.e. negatives etc
ed_data$sr_ckd_age <- ifelse(ed_data$sr_ckd_age < 0, NA, ed_data$sr_ckd_age)

# GP diagnoses for ckd.
ed_data$gp_ckd_dx <- ifelse(!ed_data$read_2_ckd_gp == "", 
                               ed_data$read_2_ckd_gp, ifelse(
                                 !ed_data$read_3_ckd_gp == "",
                                 ed_data$read_3_ckd_gp, NA
                               )
)
ed_data$gp_ckd_age <- as.numeric(difftime(ed_data$event_date_ckd, ed_data$DOB, units = "days")/365.25)
#remove impossible values, i.e. negatives etc
ed_data$gp_ckd_age <- ifelse(ed_data$gp_ckd_age < 0, NA, ed_data$gp_ckd_age)

# HES diagnoses for ckd.
ed_data$hes_ckd_dx <- ifelse(!ed_data$diag_icd10_ckd == "",
                                ed_data$diag_icd10_ckd, NA)
ed_data$hes_ckd_age <- as.numeric(difftime(ed_data$date_1st_hes_ckd, ed_data$DOB, units = "days")/365.25)
#remove impossible values, i.e. negatives etc
ed_data$hes_ckd_age <- ifelse(ed_data$hes_ckd_age < 0, NA, ed_data$hes_ckd_age)

# first occurrence for ckd.
# 132033 = chronic renal failure, 132035 = unspecified renal failure
ed_data$fo_ckd_dx <- NA
ed_data$fo_ckd_dx <- ifelse(
  !is.na(ed_data$n_132033_0_0), "N18", ifelse( 
    !is.na(ed_data$n_132035_0_0), "M31", NA))
ed_data$fo_ckd_age <- ifelse(
  !is.na(ed_data$fo_ckd_dx), as.numeric(difftime(
    pmin(ed_data$ts_132032_0_0, ed_data$ts_132034_0_0, na.rm = T),ed_data$DOB, units = "days")/365.25), NA
  )

#remove impossible values, i.e. negatives etc
ed_data$fo_ckd_age <- ifelse(ed_data$fo_ckd_age < 0, NA, ed_data$fo_ckd_age)

#and now combining them all into a yes/no ckd variable.
ed_data$comb_ckd_dx <- ifelse(
  !is.na(ed_data$sr_ckd_dx) | !is.na(ed_data$gp_ckd_dx) | 
    !is.na(ed_data$hes_ckd_dx) | !is.na(ed_data$fo_ckd_dx) #see above
  ,
  T, F
)
ed_data$comb_ckd_age <- NA
ed_data[ed_data$comb_ckd_dx == T,] <- ed_data[ed_data$comb_ckd_dx == T,] %>% 
  rowwise() %>%
  mutate(comb_ckd_age = min(sr_ckd_age, hes_ckd_age, gp_ckd_age,fo_ckd_age #see above
                             , na.rm = T))
ed_data$comb_ckd_age <- ifelse(
  ed_data$comb_ckd_age == Inf, NA, ed_data$comb_ckd_age
)

ed_data$comb_ckd_dxsource <- ifelse(!is.na(ed_data$hes_ckd_dx) & !is.na(ed_data$gp_ckd_dx) & !is.na(ed_data$fo_ckd_dx) & !is.na(ed_data$sr_ckd_dx),
                                       "HES, GP, FO, SR",
                                       ifelse(!is.na(ed_data$hes_ckd_dx) & !is.na(ed_data$gp_ckd_dx) & !is.na(ed_data$fo_ckd_dx),
                                              "HES, GP, FO",
                                              ifelse(!is.na(ed_data$hes_ckd_dx) & !is.na(ed_data$gp_ckd_dx) & !is.na(ed_data$sr_ckd_dx),
                                                     "HES, GP, SR",
                                                     ifelse(!is.na(ed_data$hes_ckd_dx) & !is.na(ed_data$fo_ckd_dx) & !is.na(ed_data$sr_ckd_dx),
                                                            "HES, FO, SR",
                                                            ifelse(!is.na(ed_data$sr_ckd_dx) & !is.na(ed_data$gp_ckd_dx) & !is.na(ed_data$fo_ckd_dx),
                                                                   "GP, FO, SR",
                                                                   ifelse(!is.na(ed_data$hes_ckd_dx) & !is.na(ed_data$gp_ckd_dx),
                                                                          "HES, GP",
                                                                          ifelse(!is.na(ed_data$hes_ckd_dx) & !is.na(ed_data$fo_ckd_dx),
                                                                                 "HES, FO",
                                                                                 ifelse(!is.na(ed_data$hes_ckd_dx) & !is.na(ed_data$sr_ckd_dx),
                                                                                        "HES, SR",
                                                                                        ifelse(!is.na(ed_data$fo_ckd_dx) & !is.na(ed_data$gp_ckd_dx),
                                                                                               "GP, FO",
                                                                                               ifelse(!is.na(ed_data$sr_ckd_dx) & !is.na(ed_data$gp_ckd_dx),
                                                                                                      "GP, SR",
                                                                                                      ifelse(!is.na(ed_data$fo_ckd_dx) & !is.na(ed_data$sr_ckd_dx),
                                                                                                             "FO, SR",
                                                                                                             ifelse(!is.na(ed_data$hes_ckd_dx),
                                                                                                                    "HES",
                                                                                                                    ifelse(!is.na(ed_data$gp_ckd_dx),
                                                                                                                           "GP",
                                                                                                                           ifelse(!is.na(ed_data$fo_ckd_dx),
                                                                                                                                  "FO",
                                                                                                                                  ifelse(!is.na(ed_data$sr_ckd_dx),
                                                                                                                                         "SR",
                                                                                                                                         NA
                                                                                                                                  )))))))))))))))

#################################
# 2C. PHENOTYPE ANALYSIS FOR DIABETES MELLITUS
#################################
#I now want to extract the self-reported diabetes diagnoses.
#to do this, I tell R to look for the self-reported diagnostic codes in the 
#many variables that encode this.
#(n_20002_*_*) is the diagnostic code
#(n_20009_*_*) is the corresponding age
sr_codes_diabetes <- c(1220, 1222, 1223) #1220 = diabetes, 1222 = T1DM, 1223 = T2DM
ed_data$sr_diabetes_dx <- NA
ed_data$sr_diabetes_age <- NA
for (i in 0:2){
  for (k in 0:33){
    variable1 = paste0("n_20002_", i, "_", k)
    variable2 = paste0("n_20009_", i, "_", k)
    if (is.null(ed_data[[variable1]])) {
      next } else {
        ed_data$sr_diabetes_dx <- ifelse(
          ed_data[[variable1]] %in% sr_codes_diabetes,
          ed_data[[variable1]], ed_data$sr_diabetes_dx)
        ed_data$sr_diabetes_age <- ifelse(
          ed_data[[variable1]] %in% sr_codes_diabetes,
          ed_data[[variable2]] , ed_data$sr_diabetes_age)
      }
  }
}
#remove impossible values, i.e. negatives etc
ed_data$sr_diabetes_age <- ifelse(ed_data$sr_diabetes_age < 0, NA, ed_data$sr_diabetes_age)

# GP diagnoses for diabetes.
ed_data$gp_diabetes_dx <- ifelse(!ed_data$read_2_diabetes_gp == "", 
                               ed_data$read_2_diabetes_gp, ifelse(
                                 !ed_data$read_3_diabetes_gp == "",
                                 ed_data$read_3_diabetes_gp, NA
                               )
)
ed_data$gp_diabetes_age <- as.numeric(difftime(ed_data$event_date_diabetes, ed_data$DOB, units = "days")/365.25)
# remove impossible values such as negatives
ed_data$gp_diabetes_age <- ifelse(ed_data$gp_diabetes_age < 0, NA, ed_data$gp_diabetes_age)


# HES diagnoses for diabetes.
ed_data$hes_diabetes_dx <- ifelse(!ed_data$diag_icd10_diabetes == "",
                                ed_data$diag_icd10_diabetes, NA)
ed_data$hes_diabetes_age <- as.numeric(difftime(ed_data$date_1st_hes_diabetes, ed_data$DOB, units = "days")/365.25)

# first occurrence for diabetes.
# 130707 = IDDM, 130709 = NIDDM, 130711 = malnutrition DM, 130713 is other DM, 130715 = unspec DM
ed_data$fo_diabetes_dx <- NA
ed_data$fo_diabetes_dx <- ifelse(
  !is.na(ed_data$n_130709_0_0), "E11", ifelse( 
    !is.na(ed_data$n_130707_0_0), "E10", ifelse(
      !is.na(ed_data$n_130711_0_0), "E12", ifelse(
        !is.na(ed_data$n_130713_0_0), "E13", ifelse(
          !is.na(ed_data$n_130715_0_0), "E15", NA
        )))))
ed_data$fo_diabetes_age <- ifelse(
  !is.na(ed_data$fo_diabetes_dx), as.numeric(difftime(
    pmin(ed_data$ts_130708_0_0, ed_data$ts_130706_0_0, 
         ed_data$ts_130710_0_0, ed_data$ts_130712_0_0, 
         ed_data$ts_130714_0_0, na.rm = T),ed_data$DOB, units = "days")/365.25), NA
)
#remove impossible values (e.g. negative age)
ed_data$fo_diabetes_age <- ifelse(
  ed_data$fo_diabetes_age <0 , NA, ed_data$fo_diabetes_age
)

#and now combining them all into a yes/no diabetes variable.
ed_data$comb_diabetes_dx <- ifelse(
  !is.na(ed_data$sr_diabetes_dx) | !is.na(ed_data$gp_diabetes_dx) | 
    !is.na(ed_data$hes_diabetes_dx) | !is.na(ed_data$fo_diabetes_dx) #see above
  ,
  T, F
)
ed_data$comb_diabetes_age <- NA
ed_data[ed_data$comb_diabetes_dx == T,] <- ed_data[ed_data$comb_diabetes_dx == T,] %>% 
  rowwise() %>%
  mutate(comb_diabetes_age = min(sr_diabetes_age, hes_diabetes_age, gp_diabetes_age,fo_diabetes_age #see above
                             , na.rm = T))
ed_data$comb_diabetes_age <- ifelse(
  ed_data$comb_diabetes_age == Inf, NA, ed_data$comb_diabetes_age
)


ed_data$comb_diabetes_dxsource <- ifelse(!is.na(ed_data$hes_diabetes_dx) & !is.na(ed_data$gp_diabetes_dx) & !is.na(ed_data$fo_diabetes_dx) & !is.na(ed_data$sr_diabetes_dx),
                                       "HES, GP, FO, SR",
                                       ifelse(!is.na(ed_data$hes_diabetes_dx) & !is.na(ed_data$gp_diabetes_dx) & !is.na(ed_data$fo_diabetes_dx),
                                              "HES, GP, FO",
                                              ifelse(!is.na(ed_data$hes_diabetes_dx) & !is.na(ed_data$gp_diabetes_dx) & !is.na(ed_data$sr_diabetes_dx),
                                                     "HES, GP, SR",
                                                     ifelse(!is.na(ed_data$hes_diabetes_dx) & !is.na(ed_data$fo_diabetes_dx) & !is.na(ed_data$sr_diabetes_dx),
                                                            "HES, FO, SR",
                                                            ifelse(!is.na(ed_data$sr_diabetes_dx) & !is.na(ed_data$gp_diabetes_dx) & !is.na(ed_data$fo_diabetes_dx),
                                                                   "GP, FO, SR",
                                                                   ifelse(!is.na(ed_data$hes_diabetes_dx) & !is.na(ed_data$gp_diabetes_dx),
                                                                          "HES, GP",
                                                                          ifelse(!is.na(ed_data$hes_diabetes_dx) & !is.na(ed_data$fo_diabetes_dx),
                                                                                 "HES, FO",
                                                                                 ifelse(!is.na(ed_data$hes_diabetes_dx) & !is.na(ed_data$sr_diabetes_dx),
                                                                                        "HES, SR",
                                                                                        ifelse(!is.na(ed_data$fo_diabetes_dx) & !is.na(ed_data$gp_diabetes_dx),
                                                                                               "GP, FO",
                                                                                               ifelse(!is.na(ed_data$sr_diabetes_dx) & !is.na(ed_data$gp_diabetes_dx),
                                                                                                      "GP, SR",
                                                                                                      ifelse(!is.na(ed_data$fo_diabetes_dx) & !is.na(ed_data$sr_diabetes_dx),
                                                                                                             "FO, SR",
                                                                                                             ifelse(!is.na(ed_data$hes_diabetes_dx),
                                                                                                                    "HES",
                                                                                                                    ifelse(!is.na(ed_data$gp_diabetes_dx),
                                                                                                                           "GP",
                                                                                                                           ifelse(!is.na(ed_data$fo_diabetes_dx),
                                                                                                                                  "FO",
                                                                                                                                  ifelse(!is.na(ed_data$sr_diabetes_dx),
                                                                                                                                         "SR",
                                                                                                                                         NA
                                                                                                                                  )))))))))))))))

#################################
# 2D. PHENOTYPE ANALYSIS FOR HYPERTENSION
#################################
#I now want to extract the self-reported hypertension diagnoses.
#to do this, I tell R to look for the self-reported diagnostic codes in the 
#many variables that encode this.
#(n_20002_*_*) is the diagnostic code
#(n_20009_*_*) is the corresponding age
sr_codes_hypertension <- c(1065, 1072) #1065 hypertension, 1072 essential hypertension
ed_data$sr_hypertension_dx <- NA
ed_data$sr_hypertension_age <- NA
for (i in 0:2){
  for (k in 0:33){
    variable1 = paste0("n_20002_", i, "_", k)
    variable2 = paste0("n_20009_", i, "_", k)
    if (is.null(ed_data[[variable1]])) {
      next } else {
        ed_data$sr_hypertension_dx <- ifelse(
          ed_data[[variable1]] %in% sr_codes_hypertension,
          ed_data[[variable1]], ed_data$sr_hypertension_dx)
        ed_data$sr_hypertension_age <- ifelse(
          ed_data[[variable1]] %in% sr_codes_hypertension,
          ed_data[[variable2]] , ed_data$sr_hypertension_age)
      }
  }
}
#remove impossible values, i.e. negatives etc
ed_data$sr_hypertension_age <- ifelse(ed_data$sr_hypertension_age < 0, NA, ed_data$sr_hypertension_age)

# GP diagnoses for hypertension.
ed_data$gp_hypertension_dx <- ifelse(!ed_data$read_2_hypertension_gp == "", 
                                     ed_data$read_2_hypertension_gp, ifelse(
                                       !ed_data$read_3_hypertension_gp == "",
                                       ed_data$read_3_hypertension_gp, NA
                                     )
)
ed_data$gp_hypertension_age <- as.numeric(difftime(ed_data$event_date_hypertension, ed_data$DOB, units = "days")/365.25)
# remove impossible values such as negatives
ed_data$gp_hypertension_age <- ifelse(ed_data$gp_hypertension_age < 0, NA, ed_data$gp_hypertension_age)


# HES diagnoses for hypertension.
ed_data$hes_hypertension_dx <- ifelse(!ed_data$diag_icd10_hypertension == "",
                                      ed_data$diag_icd10_hypertension, NA)
ed_data$hes_hypertension_age <- as.numeric(difftime(ed_data$date_1st_hes_hypertension, ed_data$DOB, units = "days")/365.25)

# first occurrence for hypertension.
# 131287 = essential hypertension, 131295 = secondary hypertension
ed_data$fo_hypertension_dx <- NA
ed_data$fo_hypertension_dx <- ifelse(
  !is.na(ed_data$n_131287_0_0), "I10", ifelse( 
    !is.na(ed_data$n_131295_0_0), "E15", NA
  ))
ed_data$fo_hypertension_age <- ifelse(
  !is.na(ed_data$fo_hypertension_dx), as.numeric(difftime(
    pmin(ed_data$ts_131286_0_0, ed_data$ts_131295_0_0, na.rm = T),ed_data$DOB, units = "days")/365.25), NA
)


#and now combining them all into a yes/no hypertension variable.
ed_data$comb_hypertension_dx <- ifelse(
  !is.na(ed_data$sr_hypertension_dx) | !is.na(ed_data$gp_hypertension_dx) | 
    !is.na(ed_data$hes_hypertension_dx) | !is.na(ed_data$fo_hypertension_dx) #see above
  ,
  T, F
)
ed_data$comb_hypertension_age <- NA
ed_data[ed_data$comb_hypertension_dx == T,] <- ed_data[ed_data$comb_hypertension_dx == T,] %>% 
  rowwise() %>%
  mutate(comb_hypertension_age = min(sr_hypertension_age, hes_hypertension_age, gp_hypertension_age,fo_hypertension_age #see above
                                     , na.rm = T))
ed_data$comb_hypertension_age <- ifelse(
  ed_data$comb_hypertension_age == Inf, NA, ed_data$comb_hypertension_age
)

ed_data$comb_hypertension_dxsource <- ifelse(!is.na(ed_data$hes_hypertension_dx) & !is.na(ed_data$gp_hypertension_dx) & !is.na(ed_data$fo_hypertension_dx) & !is.na(ed_data$sr_hypertension_dx),
                                             "HES, GP, FO, SR",
                                             ifelse(!is.na(ed_data$hes_hypertension_dx) & !is.na(ed_data$gp_hypertension_dx) & !is.na(ed_data$fo_hypertension_dx),
                                                    "HES, GP, FO",
                                                    ifelse(!is.na(ed_data$hes_hypertension_dx) & !is.na(ed_data$gp_hypertension_dx) & !is.na(ed_data$sr_hypertension_dx),
                                                           "HES, GP, SR",
                                                           ifelse(!is.na(ed_data$hes_hypertension_dx) & !is.na(ed_data$fo_hypertension_dx) & !is.na(ed_data$sr_hypertension_dx),
                                                                  "HES, FO, SR",
                                                                  ifelse(!is.na(ed_data$sr_hypertension_dx) & !is.na(ed_data$gp_hypertension_dx) & !is.na(ed_data$fo_hypertension_dx),
                                                                         "GP, FO, SR",
                                                                         ifelse(!is.na(ed_data$hes_hypertension_dx) & !is.na(ed_data$gp_hypertension_dx),
                                                                                "HES, GP",
                                                                                ifelse(!is.na(ed_data$hes_hypertension_dx) & !is.na(ed_data$fo_hypertension_dx),
                                                                                       "HES, FO",
                                                                                       ifelse(!is.na(ed_data$hes_hypertension_dx) & !is.na(ed_data$sr_hypertension_dx),
                                                                                              "HES, SR",
                                                                                              ifelse(!is.na(ed_data$fo_hypertension_dx) & !is.na(ed_data$gp_hypertension_dx),
                                                                                                     "GP, FO",
                                                                                                     ifelse(!is.na(ed_data$sr_hypertension_dx) & !is.na(ed_data$gp_hypertension_dx),
                                                                                                            "GP, SR",
                                                                                                            ifelse(!is.na(ed_data$fo_hypertension_dx) & !is.na(ed_data$sr_hypertension_dx),
                                                                                                                   "FO, SR",
                                                                                                                   ifelse(!is.na(ed_data$hes_hypertension_dx),
                                                                                                                          "HES",
                                                                                                                          ifelse(!is.na(ed_data$gp_hypertension_dx),
                                                                                                                                 "GP",
                                                                                                                                 ifelse(!is.na(ed_data$fo_hypertension_dx),
                                                                                                                                        "FO",
                                                                                                                                        ifelse(!is.na(ed_data$sr_hypertension_dx),
                                                                                                                                               "SR",
                                                                                                                                               NA
                                                                                                                                        )))))))))))))))


#################################
# 2E. PHENOTYPE ANALYSIS FOR ISCHAEMIC HEART DISEASE
#################################
#I now want to extract the self-reported ihd diagnoses.
#to do this, I tell R to look for the self-reported diagnostic codes in the 
#many variables that encode this.
#(n_20002_*_*) is the diagnostic code
#(n_20009_*_*) is the corresponding age
sr_codes_ihd <- c(1074, 1075) #1074 angina, 1075 heart attack/myocardial infarction
ed_data$sr_ihd_dx <- NA
ed_data$sr_ihd_age <- NA
for (i in 0:2){
  for (k in 0:33){
    variable1 = paste0("n_20002_", i, "_", k)
    variable2 = paste0("n_20009_", i, "_", k)
    if (is.null(ed_data[[variable1]])) {
      next } else {
        ed_data$sr_ihd_dx <- ifelse(
          ed_data[[variable1]] %in% sr_codes_ihd,
          ed_data[[variable1]], ed_data$sr_ihd_dx)
        ed_data$sr_ihd_age <- ifelse(
          ed_data[[variable1]] %in% sr_codes_ihd,
          ed_data[[variable2]] , ed_data$sr_ihd_age)
      }
  }
}
#remove impossible values, i.e. negatives etc
ed_data$sr_ihd_age <- ifelse(ed_data$sr_ihd_age < 0, NA, ed_data$sr_ihd_age)

# GP diagnoses for ihd.
ed_data$gp_ihd_dx <- ifelse(!ed_data$read_2_ihd_gp == "", 
                            ed_data$read_2_ihd_gp, ifelse(
                              !ed_data$read_3_ihd_gp == "",
                              ed_data$read_3_ihd_gp, NA
                            )
)
ed_data$gp_ihd_age <- as.numeric(difftime(ed_data$event_date_ihd, ed_data$DOB, units = "days")/365.25)
# remove impossible values such as negatives
ed_data$gp_ihd_age <- ifelse(ed_data$gp_ihd_age < 0, NA, ed_data$gp_ihd_age)


# HES diagnoses for ihd.
ed_data$hes_ihd_dx <- ifelse(!ed_data$diag_icd10_ihd == "",
                             ed_data$diag_icd10_ihd, NA)
ed_data$hes_ihd_age <- as.numeric(difftime(ed_data$date_1st_hes_ihd, ed_data$DOB, units = "days")/365.25)

# first occurrence for ihd.
# 131297 angina pectoris, 131299 acute MI, 131301 subsequent MI, 131303 complications following acute mI, 131305 other acute IHD, 131307 chronic IHD
ed_data$fo_ihd_dx <- NA
ed_data$fo_ihd_dx <- ifelse(
  !is.na(ed_data$n_131297_0_0), "I20", ifelse( 
    !is.na(ed_data$n_131299_0_0), "I21", ifelse(
      !is.na(ed_data$n_131301_0_0), "I22", ifelse(
        !is.na(ed_data$n_131303_0_0), "I23", ifelse(
          !is.na(ed_data$n_131305_0_0), "I24", ifelse(
            !is.na(ed_data$n_131307_0_0), "I25", NA
          ))))))
ed_data$fo_ihd_age <- ifelse(
  !is.na(ed_data$fo_ihd_dx), as.numeric(difftime(
    pmin(ed_data$ts_131296_0_0, ed_data$ts_131298_0_0, 
         ed_data$ts_131300_0_0, ed_data$ts_131302_0_0,
         ed_data$ts_131304_0_0, ed_data$ts_131306_0_0, na.rm = T),
    ed_data$DOB, units = "days")/365.25), NA
)
#remove impossible values (e.g. negative age)
ed_data$fo_ihd_age <- ifelse(
  ed_data$fo_ihd_age < 0, NA, ed_data$fo_ihd_age
)

#and now combining them all into a yes/no ihd variable.
ed_data$comb_ihd_dx <- ifelse(
  !is.na(ed_data$sr_ihd_dx) | !is.na(ed_data$gp_ihd_dx) | 
    !is.na(ed_data$hes_ihd_dx) | !is.na(ed_data$fo_ihd_dx) #see above
  ,
  T, F
)
ed_data$comb_ihd_age <- NA
ed_data[ed_data$comb_ihd_dx == T,] <- ed_data[ed_data$comb_ihd_dx == T,] %>% 
  rowwise() %>%
  mutate(comb_ihd_age = min(sr_ihd_age, hes_ihd_age, gp_ihd_age,fo_ihd_age #see above
                            , na.rm = T))
ed_data$comb_ihd_age <- ifelse(
  ed_data$comb_ihd_age == Inf, NA, ed_data$comb_ihd_age
)

ed_data$comb_ihd_dxsource <- ifelse(!is.na(ed_data$hes_ihd_dx) & !is.na(ed_data$gp_ihd_dx) & !is.na(ed_data$fo_ihd_dx) & !is.na(ed_data$sr_ihd_dx),
                                    "HES, GP, FO, SR",
                                    ifelse(!is.na(ed_data$hes_ihd_dx) & !is.na(ed_data$gp_ihd_dx) & !is.na(ed_data$fo_ihd_dx),
                                           "HES, GP, FO",
                                           ifelse(!is.na(ed_data$hes_ihd_dx) & !is.na(ed_data$gp_ihd_dx) & !is.na(ed_data$sr_ihd_dx),
                                                  "HES, GP, SR",
                                                  ifelse(!is.na(ed_data$hes_ihd_dx) & !is.na(ed_data$fo_ihd_dx) & !is.na(ed_data$sr_ihd_dx),
                                                         "HES, FO, SR",
                                                         ifelse(!is.na(ed_data$sr_ihd_dx) & !is.na(ed_data$gp_ihd_dx) & !is.na(ed_data$fo_ihd_dx),
                                                                "GP, FO, SR",
                                                                ifelse(!is.na(ed_data$hes_ihd_dx) & !is.na(ed_data$gp_ihd_dx),
                                                                       "HES, GP",
                                                                       ifelse(!is.na(ed_data$hes_ihd_dx) & !is.na(ed_data$fo_ihd_dx),
                                                                              "HES, FO",
                                                                              ifelse(!is.na(ed_data$hes_ihd_dx) & !is.na(ed_data$sr_ihd_dx),
                                                                                     "HES, SR",
                                                                                     ifelse(!is.na(ed_data$fo_ihd_dx) & !is.na(ed_data$gp_ihd_dx),
                                                                                            "GP, FO",
                                                                                            ifelse(!is.na(ed_data$sr_ihd_dx) & !is.na(ed_data$gp_ihd_dx),
                                                                                                   "GP, SR",
                                                                                                   ifelse(!is.na(ed_data$fo_ihd_dx) & !is.na(ed_data$sr_ihd_dx),
                                                                                                          "FO, SR",
                                                                                                          ifelse(!is.na(ed_data$hes_ihd_dx),
                                                                                                                 "HES",
                                                                                                                 ifelse(!is.na(ed_data$gp_ihd_dx),
                                                                                                                        "GP",
                                                                                                                        ifelse(!is.na(ed_data$fo_ihd_dx),
                                                                                                                               "FO",
                                                                                                                               ifelse(!is.na(ed_data$sr_ihd_dx),
                                                                                                                                      "SR",
                                                                                                                                      NA
                                                                                                                               )))))))))))))))

#################################
# 2F. PHENOTYPE ANALYSIS FOR CEREBRAL INFARCTION/TIA
#################################
#I now want to extract the self-reported cerebral infarction/TIA (cvd) diagnoses.
#to do this, I tell R to look for the self-reported diagnostic codes in the 
#many variables that encode this.
#(n_20002_*_*) is the diagnostic code
#(n_20009_*_*) is the corresponding age
sr_codes_cvd <- c(1081, 1082) #1081 stroke, 1082 TIA
ed_data$sr_cvd_dx <- NA
ed_data$sr_cvd_age <- NA
for (i in 0:2){
  for (k in 0:33){
    variable1 = paste0("n_20002_", i, "_", k)
    variable2 = paste0("n_20009_", i, "_", k)
    if (is.null(ed_data[[variable1]])) {
      next } else {
        ed_data$sr_cvd_dx <- ifelse(
          ed_data[[variable1]] %in% sr_codes_cvd,
          ed_data[[variable1]], ed_data$sr_cvd_dx)
        ed_data$sr_cvd_age <- ifelse(
          ed_data[[variable1]] %in% sr_codes_cvd,
          ed_data[[variable2]] , ed_data$sr_cvd_age)
      }
  }
}
#remove impossible values, i.e. negatives etc
ed_data$sr_cvd_age <- ifelse(ed_data$sr_cvd_age < 0, NA, ed_data$sr_cvd_age)

# GP diagnoses for cvd.
ed_data$gp_cvd_dx <- ifelse(!ed_data$read_2_cvd_gp == "", 
                            ed_data$read_2_cvd_gp, ifelse(
                              !ed_data$read_3_cvd_gp == "",
                              ed_data$read_3_cvd_gp, NA
                            )
)
ed_data$gp_cvd_age <- as.numeric(difftime(ed_data$event_date_cvd, ed_data$DOB, units = "days")/365.25)
# remove impossible values such as negatives
ed_data$gp_cvd_age <- ifelse(ed_data$gp_cvd_age < 0, NA, ed_data$gp_cvd_age)


# HES diagnoses for cvd.
ed_data$hes_cvd_dx <- ifelse(!ed_data$diag_icd10_cvd == "",
                             ed_data$diag_icd10_cvd, NA)
ed_data$hes_cvd_age <- as.numeric(difftime(ed_data$date_1st_hes_cvd, ed_data$DOB, units = "days")/365.25)

# first occurrence for cvd.
# 131367 cerebral infarction 131369 stroke not specified as haemorrhage or infarction 131057 TIA and related syndromes 131059 vascular syndromes of brain in cerebrovascular disease
ed_data$fo_cvd_dx <- NA
ed_data$fo_cvd_dx <- ifelse(
  !is.na(ed_data$n_131367_0_0), "I63", ifelse( 
    !is.na(ed_data$n_131369_0_0), "I64", ifelse(
      !is.na(ed_data$n_131057_0_0), "G45", ifelse(
        !is.na(ed_data$n_131059_0_0), "G46", NA
      ))))
ed_data$fo_cvd_age <- ifelse(
  !is.na(ed_data$fo_cvd_dx), as.numeric(difftime(
    pmin(ed_data$ts_131366_0_0, ed_data$ts_131368_0_0,
         ed_data$ts_131056_0_0, ed_data$ts_131058_0_0,
         na.rm = T),ed_data$DOB, units = "days")/365.25), NA
)

#remove impossible values (i.e. negatives)
ed_data$fo_cvd_age <- ifelse(
  ed_data$fo_cvd_age < 0, NA, ed_data$fo_cvd_age
)

#and now combining them all into a yes/no cvd variable.
ed_data$comb_cvd_dx <- ifelse(
  !is.na(ed_data$sr_cvd_dx) | !is.na(ed_data$gp_cvd_dx) | 
    !is.na(ed_data$hes_cvd_dx) | !is.na(ed_data$fo_cvd_dx) #see above
  ,
  T, F
)
ed_data$comb_cvd_age <- NA
ed_data[ed_data$comb_cvd_dx == T,] <- ed_data[ed_data$comb_cvd_dx == T,] %>% 
  rowwise() %>%
  mutate(comb_cvd_age = min(sr_cvd_age, hes_cvd_age, gp_cvd_age,fo_cvd_age #see above
                            , na.rm = T))
ed_data$comb_cvd_age <- ifelse(
  ed_data$comb_cvd_age == Inf, NA, ed_data$comb_cvd_age
)

ed_data$comb_cvd_dxsource <- ifelse(!is.na(ed_data$hes_cvd_dx) & !is.na(ed_data$gp_cvd_dx) & !is.na(ed_data$fo_cvd_dx) & !is.na(ed_data$sr_cvd_dx),
                                    "HES, GP, FO, SR",
                                    ifelse(!is.na(ed_data$hes_cvd_dx) & !is.na(ed_data$gp_cvd_dx) & !is.na(ed_data$fo_cvd_dx),
                                           "HES, GP, FO",
                                           ifelse(!is.na(ed_data$hes_cvd_dx) & !is.na(ed_data$gp_cvd_dx) & !is.na(ed_data$sr_cvd_dx),
                                                  "HES, GP, SR",
                                                  ifelse(!is.na(ed_data$hes_cvd_dx) & !is.na(ed_data$fo_cvd_dx) & !is.na(ed_data$sr_cvd_dx),
                                                         "HES, FO, SR",
                                                         ifelse(!is.na(ed_data$sr_cvd_dx) & !is.na(ed_data$gp_cvd_dx) & !is.na(ed_data$fo_cvd_dx),
                                                                "GP, FO, SR",
                                                                ifelse(!is.na(ed_data$hes_cvd_dx) & !is.na(ed_data$gp_cvd_dx),
                                                                       "HES, GP",
                                                                       ifelse(!is.na(ed_data$hes_cvd_dx) & !is.na(ed_data$fo_cvd_dx),
                                                                              "HES, FO",
                                                                              ifelse(!is.na(ed_data$hes_cvd_dx) & !is.na(ed_data$sr_cvd_dx),
                                                                                     "HES, SR",
                                                                                     ifelse(!is.na(ed_data$fo_cvd_dx) & !is.na(ed_data$gp_cvd_dx),
                                                                                            "GP, FO",
                                                                                            ifelse(!is.na(ed_data$sr_cvd_dx) & !is.na(ed_data$gp_cvd_dx),
                                                                                                   "GP, SR",
                                                                                                   ifelse(!is.na(ed_data$fo_cvd_dx) & !is.na(ed_data$sr_cvd_dx),
                                                                                                          "FO, SR",
                                                                                                          ifelse(!is.na(ed_data$hes_cvd_dx),
                                                                                                                 "HES",
                                                                                                                 ifelse(!is.na(ed_data$gp_cvd_dx),
                                                                                                                        "GP",
                                                                                                                        ifelse(!is.na(ed_data$fo_cvd_dx),
                                                                                                                               "FO",
                                                                                                                               ifelse(!is.na(ed_data$sr_cvd_dx),
                                                                                                                                      "SR",
                                                                                                                                      NA
                                                                                                                               )))))))))))))))


#################################
# 3.  PREPARING DATA FOR ANALYSIS
#################################
#select males only
ed_data <- ed_data[ed_data$sex == 1,]
#remove cases where we do not know the age of CKD or ED
ckd_unknownage <- ed_data[ed_data$comb_ckd_dx == T & is.na(ed_data$comb_ckd_age),]$n_eid
ed_unknownage <- ed_data[ed_data$comb_ed_dx == T & is.na(ed_data$comb_ed_age),]$n_eid
ed_data <- ed_data[!ed_data$n_eid %in% ckd_unknownage & !ed_data$n_eid %in% ed_unknownage,]

#smoking variable
ed_data$smoking_status <- ifelse(
  ed_data$smoking_status == 0, "never", ifelse(
    ed_data$smoking_status == 1 | ed_data$smoking_status == 2, "previous or current", NA
    )
  )


#simplified ethnicity variable
ed_data$ethnic_background %>% as.factor() %>% summary()
ed_data$ethnicity <- ifelse(ed_data$ethnic_background %in% c("4", "4001", "4002", "4003"), "black", ifelse(
  ed_data$ethnic_background %in% c("1", "1001", "1002", "1003"), "white", ifelse(
    ed_data$ethnic_background %in% c("3", "3001", "3002", "3003"), "south asian", "other or mixed"
  )
))

#birth cohort variable
ed_data$birth_cohort <- ifelse(
  ed_data$DOB < as.Date("1945-01-01", origin = "1970-01-01"), "Before 1945", ifelse(
    ed_data$DOB < as.Date("1955-01-01", origin = "1970-01-01"), "1945-1955", ifelse(
      ed_data$DOB < as.Date("1965-01-01", origin = "1970-01-01"), "1955-1965", "After 1965"
    )
  )
) %>% as.factor()

#now will need to create a new long dataframe where time dependent variables are accounted for.
#first we must create a time variable which is the time until CKD, death, or censoring, whichever comes soonest.
#we will call this variable age_event.
ed_data$dead <- ifelse(!is.na(ed_data$ts_40000_0_0), T, F)
ed_data$age_dead <- as.numeric((difftime(ed_data$ts_40000_0_0, ed_data$DOB, units = "days")/365.25))
ed_data$age_censoring <- ifelse(
  ed_data$dead == F ,
  (difftime(as.Date("2020-09-15", origin = "1970-01-01"), #last date of HES update
           ed_data$DOB, units = "days")/365.25), NA
)
ed_data$age_event <- pmin(ed_data$age_dead, ed_data$comb_ckd_age, ed_data$age_censoring, na.rm = T)
#for the cause-specific models (with death as competing risk), we will create a new variable called age_event2 so that follow up only ends when a participant dies or at end of follow up (not when they get CKD).
ed_data$age_event2 <- pmin(ed_data$age_dead, ed_data$age_censoring, na.rm = T)
#furthermore, I am making a variable called age_start, as it would not make sense to start looking at people from before the age of 35.
ed_data$age_start <- 35
#and I will exclude people who had CKD or ED before the age of 35.
ed_before35 <- ed_data[ed_data$comb_ed_age < 35,]$n_eid %>% na.omit()
ckd_before35 <- ed_data[ed_data$comb_ckd_age < 35,]$n_eid %>% na.omit()
ed_data <- ed_data[!ed_data$n_eid %in% c(ed_before35, ckd_before35),]

#then another thing to consider is that data from hospital episode statistics and primary care data 
#only go back to 1997/1998 (HES) and GP records likewise to the early nineties.
#we will truncate cases from before 01/01/2000 for this reason.
ed_data$age_trunc <- as.numeric(difftime(as.Date("2000-01-01", origin = "1970-01-01"), ed_data$DOB, units = "days")/365.25)
ed_data$age_start <- ifelse(ed_data$age_start < ed_data$age_trunc, ed_data$age_trunc, ed_data$age_start)

#and we will remove cases where people had CKD before this date as they would no longer be at risk
ckd_before2000 <- ed_data[ed_data$age_event < ed_data$age_start,]$n_eid %>% na.omit()
ed_data <- ed_data[!ed_data$n_eid %in% ckd_before2000,]

#there are a few people who developed CKD in the same month they developed ED - I will put their times slightly apart.
ed_data[!is.na(ed_data$comb_ed_age) & !is.na(ed_data$comb_ckd_age),]$comb_ed_age <- ifelse(
  ed_data[!is.na(ed_data$comb_ed_age) & !is.na(ed_data$comb_ckd_age),]$comb_ed_age == ed_data[!is.na(ed_data$comb_ed_age) & !is.na(ed_data$comb_ckd_age),]$comb_ckd_age,
  ed_data[!is.na(ed_data$comb_ed_age) & !is.na(ed_data$comb_ckd_age),]$comb_ed_age - 0.05, ed_data[!is.na(ed_data$comb_ed_age) & !is.na(ed_data$comb_ckd_age),]$comb_ed_age
)
#there are also a few people who developed CKD in the same month they died - same change will be made.
ed_data[!is.na(ed_data$comb_ckd_age) & !is.na(ed_data$age_dead),]$comb_ckd_age <- ifelse(
  ed_data[!is.na(ed_data$comb_ckd_age) & !is.na(ed_data$age_dead),]$comb_ckd_age == ed_data[!is.na(ed_data$comb_ckd_age) & !is.na(ed_data$age_dead),]$age_dead,
  ed_data[!is.na(ed_data$comb_ckd_age) & !is.na(ed_data$age_dead),]$comb_ckd_age - 0.05, ed_data[!is.na(ed_data$comb_ckd_age) & !is.na(ed_data$age_dead),]$comb_ckd_age
)
#in the analyses I am only interested in ED prior to CKD - so I add a different variable called **pre
ed_data$comb_ed_dx_pre <- 
  ifelse(is.na(ed_data$comb_ckd_age) | is.na(ed_data$comb_ed_age), ed_data$comb_ed_dx, ifelse(
    ed_data$comb_ckd_age < ed_data$comb_ed_age, F, ed_data$comb_ed_dx
  ))




#################################
# 4.  TABLES
#################################

#data for baseline table
ed_data$comb_ed_dx_pre %>% summary()
ed_data[ed_data$comb_ed_dx_pre == T,]$birth_cohort %>% summary()
ed_data[ed_data$comb_ed_dx_pre == F,]$birth_cohort %>% summary()
ed_data[ed_data$comb_ed_dx_pre == T,]$comb_ed_age %>% summary()
ed_data[ed_data$comb_ed_dx_pre == T,]$comb_ed_age %>% sd(na.rm = T)
ed_data[ed_data$comb_ed_dx_pre == T,]$smoking_status %>% as.factor() %>% summary()
ed_data[ed_data$comb_ed_dx_pre == F,]$smoking_status %>% as.factor() %>% summary()
ed_data[ed_data$comb_ed_dx_pre == T,]$TDI %>% summary()
ed_data[ed_data$comb_ed_dx_pre == F,]$TDI %>% summary()
ed_data[ed_data$comb_ed_dx_pre == T,]$BMI %>% summary()
ed_data[ed_data$comb_ed_dx_pre == T,]$BMI %>% sd(na.rm = T)
ed_data[ed_data$comb_ed_dx_pre == F,]$BMI %>% summary()
ed_data[ed_data$comb_ed_dx_pre == F,]$BMI %>% sd(na.rm = T)
ed_data[ed_data$comb_ed_dx_pre == T,]$ethnicity %>% as.factor() %>% summary()
ed_data[ed_data$comb_ed_dx_pre == F,]$ethnicity %>% as.factor() %>% summary()
ed_data[ed_data$comb_ed_dx_pre == T,]$comb_hypertension_dx %>% summary()
ed_data[ed_data$comb_ed_dx_pre == F,]$comb_hypertension_dx %>% summary()
ed_data[ed_data$comb_ed_dx_pre == T,]$comb_hypertension_age %>% summary()
ed_data[ed_data$comb_ed_dx_pre == F,]$comb_hypertension_age %>% summary()
ed_data[ed_data$comb_ed_dx_pre == T,]$comb_hypertension_age %>% sd(na.rm = T)
ed_data[ed_data$comb_ed_dx_pre == F,]$comb_hypertension_age %>% sd(na.rm = T)
(ed_data[ed_data$comb_ed_dx_pre == T,]$comb_hypertension_age < ed_data[ed_data$comb_ed_dx_pre == T,]$comb_ed_age) %>% summary()
ed_data[ed_data$comb_ed_dx_pre == T,]$comb_diabetes_dx %>% summary()
ed_data[ed_data$comb_ed_dx_pre == F,]$comb_diabetes_dx %>% summary()
ed_data[ed_data$comb_ed_dx_pre == T,]$comb_diabetes_age %>% summary()
ed_data[ed_data$comb_ed_dx_pre == F,]$comb_diabetes_age %>% summary()
ed_data[ed_data$comb_ed_dx_pre == T,]$comb_diabetes_age %>% sd(na.rm = T)
ed_data[ed_data$comb_ed_dx_pre == F,]$comb_diabetes_age %>% sd(na.rm = T)
(ed_data[ed_data$comb_ed_dx_pre == T,]$comb_diabetes_age < ed_data[ed_data$comb_ed_dx_pre == T,]$comb_ed_age) %>% summary()
ed_data[ed_data$comb_ed_dx_pre == T,]$comb_ihd_dx %>% summary()
ed_data[ed_data$comb_ed_dx_pre == F,]$comb_ihd_dx %>% summary()
ed_data[ed_data$comb_ed_dx_pre == T,]$comb_ihd_age %>% summary()
ed_data[ed_data$comb_ed_dx_pre == F,]$comb_ihd_age %>% summary()
ed_data[ed_data$comb_ed_dx_pre == T,]$comb_ihd_age %>% sd(na.rm = T)
ed_data[ed_data$comb_ed_dx_pre == F,]$comb_ihd_age %>% sd(na.rm = T)
(ed_data[ed_data$comb_ed_dx_pre == T,]$comb_ihd_age < ed_data[ed_data$comb_ed_dx_pre == T,]$comb_ed_age) %>% summary()
ed_data[ed_data$comb_ed_dx_pre == T,]$comb_cvd_dx %>% summary()
ed_data[ed_data$comb_ed_dx_pre == F,]$comb_cvd_dx %>% summary()
ed_data[ed_data$comb_ed_dx_pre == T,]$comb_cvd_age %>% summary()
ed_data[ed_data$comb_ed_dx_pre == F,]$comb_cvd_age %>% summary()
ed_data[ed_data$comb_ed_dx_pre == T,]$comb_cvd_age %>% sd(na.rm = T)
ed_data[ed_data$comb_ed_dx_pre == F,]$comb_cvd_age %>% sd(na.rm = T)
(ed_data[ed_data$comb_ed_dx_pre == T,]$comb_cvd_age < ed_data[ed_data$comb_ed_dx_pre == T,]$comb_ed_age) %>% summary()
ed_data[ed_data$comb_ed_dx_pre == T,]$comb_ckd_dx %>% summary()
ed_data[ed_data$comb_ed_dx_pre == F,]$comb_ckd_dx %>% summary()
ed_data[ed_data$comb_ed_dx_pre == T,]$comb_ckd_age %>% summary()
ed_data[ed_data$comb_ed_dx_pre == F,]$comb_ckd_age %>% summary()
ed_data[ed_data$comb_ed_dx_pre == T,]$comb_ckd_age %>% sd(na.rm = T)
ed_data[ed_data$comb_ed_dx_pre == F,]$comb_ckd_age %>% sd(na.rm = T)
(ed_data[ed_data$comb_ed_dx_pre == T & ed_data$comb_ckd_age > ed_data$comb_ed_age,]$comb_ckd_age - 
    ed_data[ed_data$comb_ed_dx_pre == T & ed_data$comb_ckd_age > ed_data$comb_ed_age,]$comb_ed_age) %>% summary()
ed_data[ed_data$comb_ed_dx_pre == T,]$dead %>% summary()
ed_data[ed_data$comb_ed_dx_pre == F,]$dead %>% summary()
ed_data[ed_data$comb_ed_dx_pre == T,]$age_dead %>% summary()
ed_data[ed_data$comb_ed_dx_pre == F,]$age_dead %>% summary()
ed_data[ed_data$comb_ed_dx_pre == T,]$age_dead %>% sd(na.rm = T)
ed_data[ed_data$comb_ed_dx_pre == F,]$age_dead %>% sd(na.rm = T)

#follow up time:
ed_data[ed_data$comb_ed_dx_pre == T & ed_data$dead == F,]$age_censoring %>% summary()
ed_data[ed_data$comb_ed_dx_pre == T & ed_data$dead == F,]$age_censoring %>% sd(na.rm = T)
ed_data[ed_data$comb_ed_dx_pre == F & ed_data$dead == F,]$age_censoring %>% summary()
ed_data[ed_data$comb_ed_dx_pre == F & ed_data$dead == F,]$age_censoring %>% sd(na.rm = T)


#aggregate person-years for risk of CKD:
ed_data$person_years_ed_ckd <- ifelse(ed_data$comb_ed_dx_pre == T, 
                                   ed_data$age_event - ed_data$comb_ed_age, 0)
ed_data$person_years_noned_ckd <- ifelse(ed_data$comb_ed_dx_pre == T & ed_data$comb_ed_age > ed_data$age_start,
                                         ed_data$comb_ed_age - ed_data$age_start, ifelse(
                                           ed_data$comb_ed_dx_pre == F,
                                           ed_data$age_event - ed_data$age_start,
                                           0
                                         ))

aggregate(ed_data$person_years_ed_ckd~ed_data$comb_ed_dx_pre, FUN=sum)
aggregate(ed_data$person_years_noned_ckd~ed_data$comb_ed_dx_pre, FUN=sum)
#aggregate person-years for death (including death after CKD):
ed_data$person_years_ed_death <- ifelse(ed_data$comb_ed_dx_pre == T, 
                                        ed_data$age_event2 - ed_data$comb_ed_age, 0)
ed_data$person_years_noned_death <- ifelse(ed_data$comb_ed_dx_pre == T & ed_data$comb_ed_age > ed_data$age_start,
                                           ed_data$comb_ed_age - ed_data$age_start, ifelse(
                                             ed_data$comb_ed_dx_pre == F,
                                             ed_data$age_event2 - ed_data$age_start,
                                             0
                                           ))

aggregate(ed_data$person_years_ed_death~ed_data$comb_ed_dx_pre, FUN=sum)
aggregate(ed_data$person_years_noned_death~ed_data$comb_ed_dx_pre, FUN=sum)
#################################
# 5.  FIRST SURVIVAL ANALYSES
#################################
#In this first set of analyses, I will run a cox regression on the data without imputing the missing data.
#I will first need to change the format of the data to a long format to allow for time-varying covariates.
td_dat <-
  tmerge(
    data1 = ed_data %>% select(n_eid, comb_ed_dx_pre, comb_ckd_dx, dead, birth_cohort, smoking_status, BMI, TDI, ethnicity),
    data2 = ed_data %>% select(n_eid, age_censoring, age_event, age_dead, dead, comb_ed_dx_pre, 
                               comb_ed_age, comb_ckd_dx, comb_ckd_age, comb_diabetes_dx,
                               comb_diabetes_age, comb_hypertension_dx, comb_hypertension_age,
                               comb_ihd_dx, comb_ihd_age, comb_cvd_dx, comb_cvd_age,
                               smoking_status, BMI, TDI, ethnicity, age_start, birth_cohort),
    id = n_eid,
    ed = tdc(comb_ed_age),
    ckd = event(comb_ckd_age, comb_ckd_dx),
    death = event(age_dead, dead),
    diabetes = tdc(comb_diabetes_age),
    hypertension = tdc(comb_hypertension_age),
    ihd = tdc(comb_ihd_age),
    stroke = tdc(comb_cvd_age),
    tstart = age_start,
    tstop = age_event
    
  )

#first I will run a "simple" cox model with ckd as event, not taking into account death.
#univariate:
coxph(
  Surv(time = tstart, time2 = tstop, event = ckd) ~ ed, 
  data = td_dat
) %>% summary()

#multivariate:
fit.coxph <- coxph(
  Surv(time = tstart, time2 = tstop, event = ckd) ~ ed + birth_cohort + 
    BMI + smoking_status + ethnicity + TDI + 
    diabetes + hypertension + ihd + stroke, 
  data = (td_dat)
) 
summary(fit.coxph)


#further down I will do a competing risk analysis. 
#first I will show a predicted survival curve (we cannot show the actual survival curve as new people "start" follow up at different ages and therefore
#the number of patients will not go down from 1.0 to 0.0 but may increase due to "new" people coming in).
#I will create a dataset with hypothetical cases of birth cohort 1945-1955, average BMI, never smokers, white ethnicity, median TDI, average IHD and stroke.
dummy <- expand.grid(ed=c(1, 0), diabetes = c(1,0), hypertension = c(1,0), birth_cohort="1945-1955", BMI = mean(ed_data$BMI, na.rm = T), 
                     smoking_status = "never", ethnicity = "white", TDI = median(ed_data$TDI, na.rm = T), 
                     ihd = mean(ed_data$comb_ihd_dx), stroke = mean(ed_data$comb_cvd_dx))
csurv <- survfit(fit.coxph, newdata=dummy)

#predicted CKD-free survival curves.
plot(csurv, xlim = c(40,80), xscale=1, ylim = c(0.5,1),
     xlab="Age", ylab="CKD-free survival", 
     col=c(1,1,"firebrick1", "firebrick1", "dodgerblue", "dodgerblue", "chartreuse3", "chartreuse3"), lty=c(2:1), lwd=2)
legend(38,0.85,
       outer(c(outer(c("ED", "No ED"),
                     c("DM", "no DM"),
                     paste, sep=", ")),
             c("HTN", "no HTN"),
             paste, sep = ", "), cex = 0.7,
       col=c(1,1,"firebrick1", "firebrick1", "dodgerblue", "dodgerblue", "chartreuse3", "chartreuse3"), lty=2:1, bty='n', lwd=2)


#this plot is the same data but then only ED vs non ED
dummy <- expand.grid(ed=c(1, 0), diabetes = 1, hypertension = 1, birth_cohort="1945-1955", BMI = mean(ed_data$BMI, na.rm = T), 
                     smoking_status = "never", ethnicity = "white", TDI = median(ed_data$TDI, na.rm = T), 
                     ihd = mean(ed_data$comb_ihd_dx), stroke = mean(ed_data$comb_cvd_dx))
csurv2 <- survfit(fit.coxph, newdata=dummy)

plot(csurv2, xlim = c(40,80), xscale=1, ylim = c(0.5,1),
     xlab="Age", ylab="CKD-free survival", 
     col=2:1, lty=c(2,1), lwd=2)
legend(38,0.85,
       c("ED", "No ED"),
       cex = 0.7,
       col=2:1, lty=c(2,1), bty='n', lwd=2)

#calculating a "number needed to treat" or number needed to screen (suggested by Mark Gilchrist).
dummy3 <- expand.grid(ed = 0, diabetes = 0, hypertension = 0, birth_cohort = "1945-1955", BMI = mean(ed_data$BMI, na.rm = T), smoking_status = "never",
                      ethnicity = "white", TDI = median(ed_data$TDI, na.rm = T), ihd = 0, stroke = 0)
surv_control <- survfit(fit.coxph, newdata = dummy3)
m <- 1
while(surv_control[[2]][m] < 70) {
  m <- m + 1
}
surv_control[[2]][m] #age is 70
prob <- surv_control[[6]][m] #CKD free survival probability in a 65M without ED, DM, HTN, IHD, stroke, never smoker etc.

hr <- fit.coxph$coefficients[1] %>% exp()
UB.CI         <- (fit.coxph$coefficients[1] - (sqrt(fit.coxph$var[1])*1.96)) %>% exp()
LB.CI         <- (fit.coxph$coefficients[1] + (sqrt(fit.coxph$var[1])*1.96)) %>% exp()

#Altman DG, Andersen PK. Calculating the number needed to treat for trials where the outcome is time to an event. BMJ. 1999 Dec 4;319(7223):1492-5. doi: 10.1136/bmj.319.7223.1492. PMID: 10582940; PMCID: PMC1117211.
NNT     <-    1/(prob - prob^hr)
NNT_ub  <-    1/(prob - prob^UB.CI)
NNT_lb  <-    1/(prob - prob^LB.CI)

paste0("NNT: ", round(NNT), " (95% CI ", round(NNT_lb), " - ", round(NNT_ub), ")")

#now I will run the same analyses but then with death as competing risk.
#in order to do this, I will make a new long dataset, this time with time not stopping if a someone develops CKD but to continue afterwards (age_event2).
td_dat2 <-
  tmerge(
    data1 = ed_data %>% select(n_eid, comb_ed_dx_pre, comb_ckd_dx, dead, birth_cohort, smoking_status, BMI, TDI, ethnicity),
    data2 = ed_data %>% select(n_eid, age_censoring, age_event2, age_dead, dead, comb_ed_dx_pre, 
                               comb_ed_age, comb_ckd_dx, comb_ckd_age, comb_diabetes_dx,
                               comb_diabetes_age, comb_hypertension_dx, comb_hypertension_age,
                               comb_ihd_dx, comb_ihd_age, comb_cvd_dx, comb_cvd_age,
                               smoking_status, BMI, TDI, ethnicity, age_start, birth_cohort),
    id = n_eid,
    ed = tdc(comb_ed_age),
    ckd = event(comb_ckd_age, comb_ckd_dx),
    death = event(age_dead, dead),
    diabetes = tdc(comb_diabetes_age),
    hypertension = tdc(comb_hypertension_age),
    ihd = tdc(comb_ihd_age),
    stroke = tdc(comb_cvd_age),
    tstart = age_start,
    tstop = age_event2
    
  )

# I will create an event variable for multi-state model
temp <- with(td_dat2, ifelse(death==1, 2, ckd))
td_dat2$event <- factor(temp, 0:2, labels=c("censor", "ckd", "death"))
remove(temp)

#I will fit a cause-specific cox proportional hazards model, as I do not think we need to do a Fine-Gray regression model.
#The reason is that we are not interested in the effect of ED on CKD as a purely aetiological question (see: Putter H, Fiocco M, Geskus RB .. Tutorial in biostatistics: competing risks and multi-state models. Stat Med 2007; 26: 2389–430.)
cr.fit <- coxph(Surv(tstart, tstop, event) ~ ed, data=td_dat2, id=n_eid)
summary(cr.fit)
#adjusted:
cr.fit <- coxph(Surv(tstart, tstop, event) ~ ed + birth_cohort + 
                  BMI + smoking_status + ethnicity + TDI + 
                  diabetes + hypertension + ihd + stroke, data=td_dat2, id=n_eid)
summary(cr.fit)


#to present a survival curve, we can use the same hypothetical subjects as used above.
csurv <- survfit(cr.fit, newdata=dummy)

#Predicted CKD-free survival curves.
plot(csurv[,'ckd'], xlim = c(40,80), xscale=1, xlab="Age", ylab="CKD", col=1:2, lty=c(1,1,2,2,3,3,4,4), lwd=2)
legend(38,0.22, outer(c(outer(c("ED", "No ED"), c("DM", "no DM"),paste, sep=", ")), c("HTN", "no HTN"), paste, sep = ", "), cex = 0.7, col=1:2, lty=c(1,1,2,2,3,3,4,4), bty='n', lwd=2)

#################################
# 6.  IMPUTATION OF MISSING DATA
#################################

#there is missing data for BMI, TDI, and smoking status (all <1%).
#in the above analyses any rows (= individuals) with missing data for these will have been excluded.
#to maximise data usage, I will impute these data 10 times using PMM and run all the analyses in the imputed datasets, then pool the results.

n.imp <- 10
set.seed(007) #(to get the same results every time the imputations are run)
#select subset of variables that we will use for the imputations:
temp <- ed_data %>% select(n_eid, age_event, age_start, age_event2, birth_cohort, ethnicity, age_censoring, age_dead, dead, comb_ed_dx_pre, 
                           comb_ed_age, comb_ckd_dx, comb_ckd_age, comb_diabetes_dx,
                           comb_diabetes_age, comb_hypertension_dx, comb_hypertension_age,
                           comb_ihd_dx, comb_ihd_age, comb_cvd_dx, comb_cvd_age,
                           smoking_status, BMI, TDI)
#impute the missing values:
imp <- aregImpute(~n_eid + age_event + age_start + age_event2 + age_dead + birth_cohort + age_censoring + ethnicity + dead + comb_ed_dx_pre + 
                    comb_ed_age + comb_ckd_dx + comb_ckd_age + comb_diabetes_dx +
                    comb_diabetes_age + comb_hypertension_dx + comb_hypertension_age +
                    comb_ihd_dx + comb_ihd_age + comb_cvd_dx + comb_cvd_age +
                    smoking_status + BMI + TDI,
                  data = temp,
                  nk = 0,
                  type = "pmm",
                  n.impute = n.imp)
#creat a dataframe that we will use to stack the imputed data sets in.
completed <- head(temp, 1)
completed$.imp <- 0
#select all individual imputed datasets and add them together.
for (i in 1:n.imp) { 
  imp.i <- temp 
  imputed <- as.data.frame(impute.transcan(imp, 
                                           imputation=i, 
                                           data=temp, list.out=TRUE, 
                                           pr=FALSE, check=FALSE))
  imputed <- imputed[c("BMI", "TDI", "smoking_status")] # only keep BMI, TDI, smoking status
  imp.i[names(imputed)] <- imputed # replace these in the imputed datasets
  imp.i$.imp <- i 
  completed <- rbind(completed, imp.i) 
  remove(imp.i)
  remove(imputed)
}
completed <- completed[completed$.imp %in% c(1:n.imp),] #the first line of this dataframe we created was just to stick the imputed datasets on, this line can be removed.

COEFS.unadjusted        <- SE.unadjusted   <- rep(NA,n.imp)
COEFS.adjusted          <- SE.adjusted   <- rep(NA,n.imp)

#create a variable that is an unique ID for each n_eid in each imputation
completed$n_eid.imp <- paste0(completed$n_eid, ".", completed$.imp)

#create a new long dataframe for the analyses in the imputed data:
td_dat <- 
  tmerge(
    data1 = completed %>% select(n_eid, n_eid.imp, .imp, comb_ed_dx_pre, comb_ckd_dx, dead, birth_cohort, smoking_status, BMI, TDI, ethnicity),
    data2 = completed %>% select(n_eid, n_eid.imp, .imp, age_censoring, age_event, age_dead, dead, comb_ed_dx_pre, 
                                 comb_ed_age, comb_ckd_dx, comb_ckd_age, comb_diabetes_dx,
                                 comb_diabetes_age, comb_hypertension_dx, comb_hypertension_age,
                                 comb_ihd_dx, comb_ihd_age, comb_cvd_dx, comb_cvd_age,
                                 smoking_status, BMI, TDI, ethnicity, age_start, birth_cohort),
    id = n_eid.imp,
    ed = tdc(comb_ed_age),
    ckd = event(comb_ckd_age, comb_ckd_dx),
    death = event(age_dead, dead),
    diabetes = tdc(comb_diabetes_age),
    hypertension = tdc(comb_hypertension_age),
    ihd = tdc(comb_ihd_age),
    stroke = tdc(comb_cvd_age),
    tstart = age_start,
    tstop = age_event
    
  ) 

#now we will run the same analyses as under 5. but separately in each of the imputed datasets.
for (i in 1:n.imp) {
  #univariate:
  m1 <- coxph(
    Surv(time = tstart, time2 = tstop, event = ckd) ~ ed, 
    data = td_dat[td_dat$.imp == i,]
  ) 
  
  #we will save the coefficients and standard errors for ED in these vectors
  COEFS.unadjusted[i] 	<- m1$coefficients[1]
  SE.unadjusted[i] 	    <- sqrt(m1$var[1])
  
  #multivariate:
  m2 <- coxph(
    Surv(time = tstart, time2 = tstop, event = ckd) ~ ed + birth_cohort + 
      BMI + smoking_status + ethnicity + TDI + 
      diabetes + hypertension + ihd + stroke, 
    data = td_dat[td_dat$.imp == i,]
  ) 
  #same as above
  COEFS.adjusted[i] 	<- m2$coefficients[1]
  SE.adjusted[i] 	    <- sqrt(m2$var[1])
}

# will then write a function to pool the results from each of these imputed datasets
Pool.rubin.HR      <- function(COEFS,SE,n.imp){   
  mean.coef     <- mean(COEFS)
  W 	        	<- mean(SE^2)
  B 		        <- var(COEFS)
  T.var 	      <- W + (1+1/n.imp)*B
  se.coef 	    <- sqrt(T.var)
  rm            <- (1+1/n.imp)*B/W
  df            <- (n.imp - 1)*(1+1/rm)^2
  LB.CI         <- mean.coef - (se.coef*1.96)
  UB.CI         <- mean.coef + (se.coef*1.96)
  F             <- (-mean.coef)^2/T.var
  P             <- pf(q=F, df1=1, df2=df, lower.tail = FALSE)
  HR            <- exp(mean.coef)
  HR.lower      <- exp(LB.CI)
  HR.upper      <- exp(UB.CI)
  output 	      <- c(HR, HR.lower, HR.upper, df, F, P)	
  names(output) <- c('HR', 'lower bound', 'upper bound', 'df', 'F', 'P')
  return(output)}

# pooled HR:
Pool.rubin.HR(COEFS.unadjusted, SE.unadjusted, n.imp)
Pool.rubin.HR(COEFS.adjusted, SE.adjusted, n.imp)

#we can show another predicted survival curve (I have not found out how to do this for all imputed datasets combined but given that there's <<1% missing data these results are numerically similar)
csurv <- survfit(m2, newdata=dummy)
plot(csurv, xlim = c(40,80), xscale=1, ylim = c(0.5,1),
     xlab="Age", ylab="CKD-free survival", 
     col=c(1,1,"firebrick1", "firebrick1", "dodgerblue", "dodgerblue", "chartreuse3", "chartreuse3"), lty=c(2:1), lwd=2)
legend(38,0.85,
       outer(c(outer(c("ED", "No ED"),
                     c("DM", "no DM"),
                     paste, sep=", ")),
             c("HTN", "no HTN"),
             paste, sep = ", "), cex = 0.7,
       col=c(1,1,"firebrick1", "firebrick1", "dodgerblue", "dodgerblue", "chartreuse3", "chartreuse3"), lty=2:1, bty='n', lwd=2)

surv_control <- survfit(m2, newdata = dummy3)
m <- 1
while(surv_control[[2]][m] < 70) {
  m <- m + 1
}
surv_control[[2]][m] #age is 70
prob <- surv_control[[6]][m] #CKD free survival probability in a 70M without ED, DM, HTN, IHD, stroke, never smoker etc.

hr            <- m2$coefficients[1] %>% exp()
UB.CI         <- (m2$coefficients[1] - (sqrt(m2$var[1])*1.96)) %>% exp()
LB.CI         <- (m2$coefficients[1] + (sqrt(m2$var[1])*1.96)) %>% exp()

#Altman DG, Andersen PK. Calculating the number needed to treat for trials where the outcome is time to an event. BMJ. 1999 Dec 4;319(7223):1492-5. doi: 10.1136/bmj.319.7223.1492. PMID: 10582940; PMCID: PMC1117211.
NNT     <-    1/(prob - prob^hr)
NNT_ub  <-    1/(prob - prob^UB.CI)
NNT_lb  <-    1/(prob - prob^LB.CI)

print(paste0("NNT: ", round(NNT), " (95% CI ", round(NNT_lb), " - ", round(NNT_ub), ")"))

# now we will run the same analyses in the cause-specific cox model allowing for competing risk of death:
td_dat2 <-
  tmerge(
    data1 = completed %>% select(n_eid, .imp, comb_ed_dx_pre, comb_ckd_dx, dead, birth_cohort, smoking_status, BMI, TDI, ethnicity),
    data2 = completed %>% select(n_eid, age_censoring, age_event2, age_dead, dead, comb_ed_dx_pre, 
                                 comb_ed_age, comb_ckd_dx, comb_ckd_age, comb_diabetes_dx,
                                 comb_diabetes_age, comb_hypertension_dx, comb_hypertension_age,
                                 comb_ihd_dx, comb_ihd_age, comb_cvd_dx, comb_cvd_age,
                                 smoking_status, BMI, TDI, ethnicity, age_start, birth_cohort),
    id = n_eid,
    ed = tdc(comb_ed_age),
    ckd = event(comb_ckd_age, comb_ckd_dx),
    death = event(age_dead, dead),
    diabetes = tdc(comb_diabetes_age),
    hypertension = tdc(comb_hypertension_age),
    ihd = tdc(comb_ihd_age),
    stroke = tdc(comb_cvd_age),
    tstart = age_start,
    tstop = age_event2
    
  )

temp <- with(td_dat2, ifelse(death==1, 2, ckd))
td_dat2$event <- factor(temp, 0:2, labels=c("censor", "ckd", "death"))
remove(temp)

COEFS.cr.unadjusted        <- SE.cr.unadjusted   <- rep(NA,n.imp)
COEFS.cr.adjusted          <- SE.cr.adjusted   <- rep(NA,n.imp)

for (i in 1:n.imp) {
cr.fit <- coxph(Surv(tstart, tstop, event) ~ ed, data=td_dat2[td_dat2$.imp == i,], id=n_eid)
  COEFS.cr.unadjusted[i] 	<- cr.fit$coefficients[1]
  SE.cr.unadjusted[i] 	    <- sqrt(cr.fit$var[1])
#adjusted:
cr.fit2 <- coxph(Surv(tstart, tstop, event) ~ ed + birth_cohort + 
                  BMI + smoking_status + ethnicity + TDI + 
                  diabetes + hypertension + ihd + stroke, data=td_dat2[td_dat2$.imp == i,], id=n_eid)
  COEFS.cr.adjusted[i] 	<- cr.fit2$coefficients[1]
  SE.cr.adjusted[i] 	    <- sqrt(cr.fit2$var[1])
  }

# pooled HR:
Pool.rubin.HR(COEFS.cr.unadjusted, SE.cr.unadjusted, n.imp)
Pool.rubin.HR(COEFS.cr.adjusted, SE.cr.adjusted, n.imp)

#################################
# 7.  PROPENSITY SCORE MATCHING
#################################

#one suggestion from one my peers was to run these analyses in a propensity score-matched cohort.
#in order to do this, we will select all individual patients at start of follow up and then match them on propensity scores.

#create empty data frame
matched <- data.frame(matrix(ncol = ncol(completed)))
names(matched) <- names(completed)
matched$.imp <- NA

# as matching requires a complete dataset without missing data, we will do the matching in each of the imputed datasets 
#then run the analyses in each of the imputed matched datasets then pool the results. 
for (i in 1:n.imp) {
  #select the cases at start of follow-up:
  baseline <- completed[completed$.imp == i,] %>% 
    group_by(n_eid) %>% 
    filter(age_start == min(age_start))
  #actual matching:
  match <- matchit(comb_ed_dx_pre ~ age_start + birth_cohort + 
                     BMI + smoking_status + ethnicity + TDI + 
                     comb_diabetes_dx + comb_hypertension_dx +
                     comb_ihd_dx + comb_cvd_dx,
                   data=baseline, 
                   method="nearest", 
                   distance="logit", 
                   discard="both", 
                   caliper=0.1, 
                   ratio=2
  )
  
  # extract matrix with matched n_eid
  matrix.match <- cbind(data.frame(baseline[row.names(match$match.matrix),"n_eid"]), 
                        data.frame(baseline[match$match.matrix,"n_eid"]))
  matrix.match <- data.frame(matrix.match, row.names=NULL) 
  
  # filter out unmatched n_eid
  matrix.match <- matrix.match[complete.cases(matrix.match),] 
  
  # select matched patients from imputed dataframe
  temp <- completed[completed$.imp == i & 
                      (completed$n_eid %in% matrix.match$n_eid | 
                         completed$n_eid %in% matrix.match$n_eid.1),] 
  # attach this to the matched cases
  matched <- rbind(matched,temp) 
  remove(temp)     
  remove(match)
  remove(matrix.match)
  print(paste0("Matching in imputed dataset number ", i, " completed"))
}
# remove row with NA
matched <- matched[!is.na(matched$n_eid),]

#the dataset matched is now a dataframe with PSM cases from each of the imputed datasets.

# for table comparing ED with non ED (matched) I will use create table one as this includes SMD (SMD <0.1 indicates balance of covariates)
vars <- c("age_start", "birth_cohort", "smoking_status", "TDI", "BMI", "ethnicity", "comb_hypertension_dx", "comb_diabetes_dx", "comb_ihd_dx", "comb_cvd_dx")
CreateTableOne(vars = vars, strata = "comb_ed_dx_pre", data = matched, test = F) %>% print(smd = T)

# raw data - person years and CKD events in ED and non ED
aggregate(matched$person_years_ed_ckd~matched$comb_ed_dx_pre, FUN=sum)
aggregate(matched$person_years_noned_ckd~matched$comb_ed_dx_pre, FUN=sum)

(matched[matched$comb_ed_dx_pre == T,]$comb_ckd_dx == T) %>% summary()
(matched[matched$comb_ed_dx_pre == F,]$comb_ckd_dx == T) %>% summary()

# raw data - person years and deaths in ED and non ED
aggregate(matched$person_years_ed_death~matched$comb_ed_dx_pre, FUN=sum)
aggregate(matched$person_years_noned_death~matched$comb_ed_dx_pre, FUN=sum)
(matched[matched$comb_ed_dx_pre == T,]$dead == T) %>% summary()
(matched[matched$comb_ed_dx_pre == F,]$dead == T) %>% summary()


# Now we will also run the survival analyses in these propensity score matched patients
COEFS.ED.psm       <- SE.ED.psm   <-  COEFS.EDadj.psm <- SE.EDadj.psm <- 
COEFS.ED.cr.psm       <- SE.ED.cr.psm   <-  COEFS.EDadj.cr.psm <- SE.EDadj.cr.psm <- 
rep(NA,n.imp)

#create long dataframe with PS matched patients only by selecting the PS matched patients from the previously created imputed long dataframe
td_dat_matched <- td_dat[td_dat$n_eid.imp %in% matched$n_eid.imp,]
td_dat2_matched <- td_dat2[td_dat2$n_eid.imp %in% matched$n_eid.imp,]

#run the survival analyses in each of the imputed + matched datasets
for (i in 1:n.imp) {
  p <- coxph(Surv(tstart, tstop, ckd) ~ ed, 
             data = td_dat_matched[td_dat_matched$.imp == i,])
  p2 <- coxph(Surv(tstart, tstop, ckd) ~ ed + birth_cohort + BMI + ethnicity +
                smoking_status + TDI + hypertension + diabetes + ihd + stroke, 
              data = td_dat_matched[td_dat_matched$.imp == i,])
  
  COEFS.ED.psm[i] 	<- p$coefficients[1]
  SE.ED.psm[i] 	    <- sqrt(p$var[1])
  
  COEFS.EDadj.psm[i] 	<- p2$coefficients[1]
  SE.EDadj.psm[i] 	    <- sqrt(p2$var[1])
  remove(p)
  remove(p2)
  #run a cause-specific cox PH as well for competing risk of death:
  cr.fit.matched <- coxph(Surv(tstart, tstop, event) ~ ed, data=td_dat2_matched[td_dat2_matched$.imp == i,], id=n_eid)
  COEFS.ED.cr.psm[i] 	<- cr.fit.matched$coefficients[1]
  SE.ED.cr.psm[i] 	    <- sqrt(cr.fit.matched$var[1])
#adjusted:
cr.fit.matched2 <- coxph(Surv(tstart, tstop, event) ~ ed + birth_cohort + 
                  BMI + smoking_status + ethnicity + TDI + 
                  diabetes + hypertension + ihd + stroke, data=td_dat2_matched[td_dat2_matched$.imp == i,], id=n_eid)
  COEFS.EDadj.cr.psm[i] 	<- cr.fit.matched2$coefficients[1]
  SE.EDadj.cr.psm[i] 	    <- sqrt(cr.fit.matched2$var[1])
  remove(cr.fit.matched)
  remove(cr.fit.matched2)
}

# pool results:
#propensity score matched analysis in simple cox regression model:
Pool.rubin.HR(COEFS.ED.psm, SE.ED.psm, n.imp)
#propensity score matched analysis in simple cox regression model with (doubly robust) adjustment:
Pool.rubin.HR(COEFS.EDadj.psm, SE.EDadj.psm, n.imp)
#propensity score matched analysis in cause-specific cox regression model with competing risk of death:
Pool.rubin.HR(COEFS.ED.cr.psm, SE.ED.cr.psm, n.imp)
#propensity score matched analysis in cause-specific cox regression model with competing risk of death and (doubly robust) adjustment:
Pool.rubin.HR(COEFS.EDadj.cr.psm, SE.EDadj.cr.psm, n.imp)
