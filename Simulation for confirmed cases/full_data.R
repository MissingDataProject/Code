library(dplyr)
# original data
merged_df <- read.csv("original_data.csv")
merged_df$clinic_num <- rowSums(merged_df[, c(16:42)], na.rm = T)
merged_df <- merged_df %>% 
  
  # case classification -- confirm, probable
  mutate(case_class = case_when(
    (fever == 1) & (clinic_num >= 3) ~ "Probable",
    detection == "POSITIVE" ~ "Confirmed", 
    TRUE ~ NA_character_
  )) %>% 
  # delete outliers of follow-up time
  filter(futime <= 40) 
merged_df <- merged_df %>% 
  mutate(missing = case_when(
    is.na(event) ~ 1,
    TRUE ~ 0
  ))
# missing logistic --> logit(Pr(missing outcome)) ~ siteid + diarrhea + loss_of_appetite
sum(is.na(merged_df$event)) / nrow(merged_df) # --> original: 10%


# Missingness model --> IPW
missing_logit <- glm(missing ~ age_years  + sex
                     + abdominal_pain + fever + hiccoughs + diarrhea + loss_of_appetite
                     + headache + bleeding + jaundice + difficulty_swallowing, 
                     family=binomial(link='logit'),
                     data = filter(merged_df, case_class == "Probable"))
summary(missing_logit)


# MAR model --> used to control % missing when introducing missingness in outcome
logit_intro <- glm(missing ~ age_years  + sex,
                   family=binomial(link='logit'),
                   data = filter(merged_df, case_class == "Probable"))
summary(logit_intro)



# -----------------------------------------------------


# load full data
full_data <- read.csv("full_data.csv")
full_data$siteid <- as.factor(full_data$siteid)
full_data$sex <- as.factor(full_data$sex)
full_data$country <- as.factor(full_data$country)
full_data$case_class <- as.factor(full_data$case_class)

##### "true" value from full data

merged_data_keep <- full_data %>% filter(case_class == "Probable")

library(survival)
library(survminer)
overall_surv <- survfit(Surv(futime, event) ~ 1, data = merged_data_keep) 
options(scipen = 999) # no scientic notation
true_CFR <- summary(overall_surv, times = 20)$surv
true_CFR




