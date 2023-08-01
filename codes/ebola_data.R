### load libary
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


#? Missingness model --> IPW
missing_logit <- glm(missing ~ age_years  + sex
                     + abdominal_pain + fever + hiccoughs + diarrhea + loss_of_appetite
                     + headache + bleeding + jaundice + difficulty_swallowing, 
                     family=binomial(link='logit'),
                     data = filter(merged_df, case_class == "Probable"))
summary(missing_logit)


# -------------------- 1st method: CCA --------------------

survival_CCA <- survfit(Surv(futime, as.numeric(event)) ~ case_class, data = merged_df) 
result_cca <- summary(survival_CCA, times = 20)
# CFR (confirmed & probable)
1 - result_cca$surv



# -------------------- 2nd method: MI --------------------

merged_df$na_estimator <- nelsonaalen(merged_df, futime, event)
merged_df$event <- as.factor(merged_df$event) # make sure meth is logreg

mi_ebola <- mice(merged_df  %>% 
                    select(futime, event, na_estimator, case_class,
                           siteid, sex, abdominal_pain, fever, hiccoughs, diarrhea, loss_of_appetite,
                           headache, bleeding, jaundice, difficulty_swallowing),
                  m = 11, maxit =20, printFlag=FALSE)

fit_mi_ebola <- with(mi_ebola, with(summary(survfit(Surv(futime, as.numeric(event)) ~ case_class), times = 20, extend = TRUE),
                                   data.frame(time, surv, std.err, lower, upper)))
# confirmed group
pool_mi_ebola <- fit_mi_ebola$analyses %>% 
  bind_rows() %>%
  filter(row_number() %% 2 == 1) %>% 
  mutate(surv_log = log(-log(surv)),
         se_trans = std.err^2 / (surv * log(surv))^2) %>%
  group_by(time) %>% 
  summarize(pooled = data.frame(pool.scalar(Q = surv_log, U = se_trans, n=Inf, k = 1)[c("qbar", "t")])) %>% 
  unpack(pooled) %>% 
  mutate(surv = exp(-exp(qbar)), 
         SE = abs(sqrt(t) * surv * log(surv)), 
         LCI = surv^(exp(1.96*SE)), 
         UCI = surv^exp(-1.96*SE))
cbind("MI", pool_mi_ebola$surv, pool_mi_ebola$SE, pool_mi_ebola$LCI, pool_mi_ebola$UCI)

# probable group
pool_mi_ebola_2 <- fit_mi_ebola$analyses %>% 
  bind_rows() %>%
  filter(row_number() %% 2 == 0) %>% 
  mutate(surv_log = log(-log(surv)),
         se_trans = std.err^2 / (surv * log(surv))^2) %>%
  group_by(time) %>% 
  summarize(pooled = data.frame(pool.scalar(Q = surv_log, U = se_trans, n=Inf, k = 1)[c("qbar", "t")])) %>% 
  unpack(pooled) %>% 
  mutate(surv = exp(-exp(qbar)), 
         SE = abs(sqrt(t) * surv * log(surv)), 
         LCI = surv^(exp(1.96*SE)), 
         UCI = surv^exp(-1.96*SE))

cbind("MI", pool_mi_ebola_2$surv, pool_mi_ebola_2$SE, pool_mi_ebola_2$LCI, pool_mi_ebola_2$UCI)



# -------------------- 3rd method: IPW --------------------
merged_df$event <- as.numeric(merged_df$event)
merged_df_complete <- merged_df %>% 
  select(usubjid, case_class, age_years, sex, abdominal_pain, fever, hiccoughs, diarrhea, loss_of_appetite, 
         headache, bleeding, jaundice, difficulty_swallowing, missing) %>% 
  filter(complete.cases(.))
ps_model_ebola <- glm(missing ~ case_class + age_years + sex + abdominal_pain + fever + hiccoughs + diarrhea + loss_of_appetite + 
                      headache + bleeding + jaundice + difficulty_swallowing,
                    family = binomial, data = merged_df)
# summary(ps_model)

# score
merged_df_complete$propensity_score <- predict(ps_model_ebola, type = "response")
# calculate weighting
# estimate weight for each patient
merged_df_complete$weight.ATE <- 1/merged_df_complete$propensity_score

merged_df_complete_new <- merge(merged_df_complete, merged_df %>% select(futime, event, usubjid), by = "usubjid")
ebola_ipw_model <- survfit(Surv(futime, as.numeric(event)) ~ case_class, weights = weight.ATE, data = merged_df_complete_new) 
result_ipw <- summary(ebola_ipw_model, times = 20)
result_ipw$surv
# bootstrapping 
sd_function <- function(data, indices = 1:nrow(data)) {
  d <- data %>% slice(indices) # allows boot to select sample
  ps_model <- glm(is.na(event) ~ age_years + sex
                      + abdominal_pain + fever + hiccoughs + diarrhea + loss_of_appetite
                      + headache + bleeding + jaundice + difficulty_swallowing, 
                      family = binomial, data = d)
  d$propensity_score <- predict(ps_model, type = "response")
  d$weight.ATE <- ifelse((is.na(d$event)),
                             (1/d$propensity_score), 
                             (1/(1-d$propensity_score)))
  fit <- survfit(Surv(futime, as.numeric(event)) ~ case_class, weights = weight.ATE, data=d)
  return(summary(fit, times = 20)$surv)
}
resample <- boot(merged_df_complete_new, statistic = sd_function, R = 200)
summary(resample)$bootSE
result_ipw$surv - 1.96 * summary(resample)$bootSE
result_ipw$surv + 1.96 * summary(resample)$bootSE
