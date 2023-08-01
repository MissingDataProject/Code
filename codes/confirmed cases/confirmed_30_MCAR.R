
#####################################
# 
#    scenario: MCAR 30% 
#
#####################################


###### load library
library(mice)
library(ipw)
library(dplyr)
library(tidyr)
library(boot)
library(survival)
library(survminer)


# -------------------- full data --------------------
# load full data
full_data <- read.csv("full_data.csv")
full_data$siteid <- as.factor(full_data$siteid)
full_data$sex <- as.factor(full_data$sex)
full_data$country <- as.factor(full_data$country)
full_data$case_class <- as.factor(full_data$case_class)

merged_data_keep <- full_data %>% filter(case_class == "Confirmed")

##### "true" value from full data
overall_surv <- survfit(Surv(futime, event) ~ case_class, data = merged_data_keep) 
true_CFR <- summary(overall_surv, times = 20)$surv
true_CFR


# -------------------- simulation study --------------------

# step 0
# build a matrix to restore values for performance measures


bias_30_MCAR <- matrix(0, 1, 3)
empse_30_MCAR <- matrix(0, 1, 3)
modelse_30_MCAR <- matrix(0, 1, 3)
rmse_30_MCAR <- matrix(0, 1, 3)
cover_30_MCAR <- matrix(0, 1, 3)


N_sim <- 1000

# -------------------- 1st method: CCA --------------------

mcar_30_cc_prob <- c()

set.seed(12345)

for (i in 1:N_sim){
  
  print(i) # know the current iteration
  
  merged_data <- merged_data_keep[sample(nrow(merged_data_keep), 1000, replace = TRUE), ]
  
  # MCAR
  merged_data$mcar_pro <- runif(nrow(merged_data), 0, 1)
  # ensure 30% missingness
  merged_data <- merged_data %>% 
    mutate(mcar_event = ifelse(mcar_pro < quantile(mcar_pro, probs = 0.7), 
                               event, NA))
  # check missing percentage
  # nrow(filter(merged_data, is.na(mcar_event))) / nrow(merged_data)
  
  
  #### 1.1 CCA 
  overall_surv_MCAR <- survfit(Surv(futime, as.numeric(mcar_event)) ~ 1, data = merged_data) 
  # estimated survival rates
  mcar_cc <- summary(overall_surv_MCAR, times = 20)
  mcar_30_cc_prob_each_run <- as.data.frame(cbind(as.numeric(i), "CCA", (mcar_cc$surv), mcar_cc$std.err, 
                                                  mcar_cc$lower, mcar_cc$upper))
  mcar_30_cc_prob <- rbind(mcar_30_cc_prob, mcar_30_cc_prob_each_run)
  
}
# covert to numeric
for (i in 3:6){
  mcar_30_cc_prob[, i] <- as.numeric(mcar_30_cc_prob[, i])
}
colnames(mcar_30_cc_prob) <- c("dataset", "method", "estimate", "se", "CI_lower", "CI_upper")

##### performance measure
bias_30_MCAR[, 1] <- mean(mcar_30_cc_prob$estimate - true_CFR)
empse_30_MCAR[, 1] <- sqrt(sum( (mcar_30_cc_prob$estimate - true_CFR)^2 ) / (N_sim - 1))
modelse_30_MCAR[, 1] <- sqrt( sum( mcar_30_cc_prob$se^2 ) / N_sim)
rmse_30_MCAR[, 1] <- sqrt( sum((mcar_30_cc_prob$estimate - true_CFR)^2) / N_sim)
cover_30_MCAR[, 1] <- mean(mcar_30_cc_prob$CI_lower <= true_CFR & true_CFR <= mcar_30_cc_prob$CI_upper)


# -------------------- 2nd method: MI --------------------

mcar_30_mi_prob <- c()

for (i in 1:N_sim){
  
  print(i) # know the current iteration
  
  merged_data <- merged_data_keep[sample(nrow(merged_data_keep), 1000, replace = TRUE), ]
  
  # MCAR
  merged_data$mcar_pro <- runif(nrow(merged_data), 0, 1)
  # ensure 10% missingness
  merged_data <- merged_data %>% 
    mutate(mcar_event = ifelse(mcar_pro < quantile(mcar_pro, probs = 0.7), 
                               event, NA))
  # check missing percentage
  # nrow(filter(merged_data, is.na(mcar_event))) / nrow(merged_data)
  
  
  #### 1.2 MI
  # Here m = 31 because missingness = 30%
  # add Nelson-Aalen estimator
  
  #######################
  # imputation model
  #######################
  
  merged_data$na_estimator <- nelsonaalen(merged_data, futime, mcar_event)
  merged_data$mcar_event <- as.factor(merged_data$mcar_event) # make sure meth is logreg
  
  mi_MCAR_30 <- mice(merged_data  %>% 
                       select(futime, mcar_event, na_estimator,
                              siteid, sex, abdominal_pain, fever, hiccoughs, diarrhea, loss_of_appetite,
                              headache, bleeding, jaundice, difficulty_swallowing),
                     m = 31, maxit = 20, printFlag=FALSE)
  
  # check if there is any NAs in the imputed dataset e.g. the 1st dataset
  # data_1 <- complete(mi_MCAR_10 )
  # sum(is.na(data_1$mcar_event))

  # MUST covert "mcar_event" to numeric rather than factor, otherwise summary will not gives "survival" results
  fit_mi_MCAR <- with(mi_MCAR_30, with(summary(survfit(Surv(futime, as.numeric(mcar_event)) ~ 1), times = 20, extend = TRUE),
                                       data.frame(time, surv, std.err, lower, upper)))
  
  pool_mi <- fit_mi_MCAR$analyses %>% 
    bind_rows() %>%
    mutate(surv_log = log(-log(surv)),
           se_trans = std.err^2 / (surv * log(surv))^2) %>%
    group_by(time) %>% 
    summarize(pooled = data.frame(pool.scalar(Q = surv_log, U = se_trans, n=Inf, k = 1)[c("qbar", "t")])) %>% 
    unpack(pooled) %>% 
    mutate(surv = exp(-exp(qbar)), 
           SE = abs(sqrt(t) * surv * log(surv)), 
           LCI = surv^(exp(1.96*SE)), 
           UCI = surv^exp(-1.96*SE))
  
  mcar_30_mi_prob_each_run <- as.data.frame(cbind(as.numeric(i), "MI", (pool_mi$surv), pool_mi$SE, 
                                                  pool_mi$LCI, pool_mi$UCI))
  mcar_30_mi_prob <- rbind(mcar_30_mi_prob, mcar_30_mi_prob_each_run)
  
}
# covert to numeric
for (i in 3:6){
  mcar_30_mi_prob[, i] <- as.numeric(mcar_30_mi_prob[, i])
}
colnames(mcar_30_mi_prob) <- c("dataset", "method", "estimate", "se", "CI_lower", "CI_upper")

##### performance measure
bias_30_MCAR[, 2] <- mean(mcar_30_mi_prob$estimate - true_CFR)
empse_30_MCAR[, 2] <- sqrt(sum( (mcar_30_mi_prob$estimate - true_CFR)^2 ) / (N_sim - 1))
modelse_30_MCAR[, 2] <- sqrt( sum( mcar_30_mi_prob$se^2 ) / N_sim)
rmse_30_MCAR[, 2] <- sqrt( sum((mcar_30_mi_prob$estimate - true_CFR)^2) / N_sim)
cover_30_MCAR[, 2] <- mean(mcar_30_mi_prob$CI_lower <= true_CFR & true_CFR <= mcar_30_mi_prob$CI_upper)





# -------------------- 3rd method: IPW --------------------

mcar_30_ipw_prob <- c()

for (i in 1:N_sim){
  
  print(i) # know the current iteration
  
  
  merged_data <- merged_data_keep[sample(nrow(merged_data_keep), 1000, replace = TRUE), ]
  
  # summary(simdata$data)
  
  # MCAR
  merged_data$mcar_pro <- runif(nrow(merged_data), 0, 1)
  # ensure 10% missingness
  merged_data <- merged_data %>% 
    mutate(mcar_event = ifelse(mcar_pro < quantile(mcar_pro, probs = 0.7), 
                               event, NA))
  # check missing percentage
  # nrow(filter(merged_data, is.na(mcar_event))) / nrow(merged_data)
  
  
  #### 1.3 IPW
  
  #######################
  # missingness model
  #######################
  
  ps_model <- glm(is.na(mcar_event) ~ age_years + siteid + sex
                  + abdominal_pain + fever + hiccoughs + diarrhea + loss_of_appetite
                  + headache + bleeding + jaundice + difficulty_swallowing, 
                  family = binomial, data = merged_data)
  # summary(ps_model)
  # score
  merged_data$propensity_score <- predict(ps_model, type = "response")
  # calculate weighting
  # estimate weight for each patient
  merged_data$weight.ATE <- ifelse((is.na(merged_data$mcar_event)),
                                   (1/merged_data$propensity_score), 
                                   (1/(1-merged_data$propensity_score)))
  
  
  #time to event analysis with weights
  mcar_30_ipw_model <- survfit(Surv(futime, as.numeric(mcar_event)) ~ 1, weights = weight.ATE, data = merged_data) 
  mcar_ipw <- summary(mcar_30_ipw_model, times = 20)
  
  # 200 bootstrap to get standard error
  sd_function <- function(data, indices = 1:nrow(data)) {
    d <- data %>% slice(indices) # allows boot to select sample
    ps_model <- glm(is.na(mcar_event) ~ age_years + siteid + sex
                    + abdominal_pain + fever + hiccoughs + diarrhea + loss_of_appetite
                    + headache + bleeding + jaundice + difficulty_swallowing, 
                    family = binomial, data = d)
    d$propensity_score <- predict(ps_model, type = "response")
    d$weight.ATE <- ifelse((is.na(d$mcar_event)),
                           (1/d$propensity_score), 
                           (1/(1-d$propensity_score)))
    fit <- survfit(Surv(futime, as.numeric(mcar_event)) ~ 1, weights = weight.ATE, data=d)
    return(summary(fit, times = 20)$surv)
  }
  resample <- boot(merged_data, statistic = sd_function, R = 200)
  
  
  mcar_30_ipw_prob_each_run <- as.data.frame(cbind(as.numeric(i), "IPW", (mcar_ipw$surv), summary(resample)$bootSE, 
                                                   mcar_ipw$surv - 1.96 * summary(resample)$bootSE, 
                                                   mcar_ipw$surv + 1.96 * summary(resample)$bootSE))
  mcar_30_ipw_prob <- rbind(mcar_30_ipw_prob, mcar_30_ipw_prob_each_run)
  
}
# covert to numeric
for (i in 3:6){
  mcar_30_ipw_prob[, i] <- as.numeric(mcar_30_ipw_prob[, i])
}
colnames(mcar_30_ipw_prob) <- c("dataset", "method", "estimate", "se", "CI_lower", "CI_upper")

##### performance measure
bias_30_MCAR[, 3] <- mean(mcar_30_ipw_prob$estimate - true_CFR)
empse_30_MCAR[, 3] <- sqrt(sum( (mcar_30_ipw_prob$estimate - true_CFR)^2 ) / (N_sim - 1))
modelse_30_MCAR[, 3] <- sqrt( sum( mcar_30_ipw_prob$se^2 ) / N_sim)
rmse_30_MCAR[, 3] <- sqrt( sum((mcar_30_ipw_prob$estimate - true_CFR)^2) / N_sim)
cover_30_MCAR[, 3] <- mean(mcar_30_ipw_prob$CI_lower <= true_CFR & true_CFR <= mcar_30_ipw_prob$CI_upper)


# -------------------- Results (Performance measures) --------------------

result_30_MCAR <- rbind(bias_30_MCAR, empse_30_MCAR, modelse_30_MCAR, rmse_30_MCAR, cover_30_MCAR)
rownames(result_30_MCAR) <- c("Bias", "Empirical SE", "Model SE", "RMSE", "Coverage")
colnames(result_30_MCAR) <- c("CCA", "MI", "IPW")
View(result_30_MCAR)


mcar_30_prob <- rbind(mcar_30_cc_prob, mcar_30_mi_prob, mcar_30_ipw_prob)
write.csv(mcar_30_prob, "/Users/apple/Desktop/new_full_data/simulation/confirmed_30_MCAR_simulation.csv")
library(rsimsum)
s1_prob_30_MCAR <- simsum(mcar_30_prob, estvarname = "estimate", true = true_CFR, se = "se", 
                          methodvar = "method", ref = "CCA")
ss1_prob_30_MCAR <- summary(s1_prob_30_MCAR, stats = c("bias", "modelse", "empse", "cover"))
table_30_MCAR <- ss1_prob_30_MCAR %>% 
  tidy() %>% 
  mutate(performance = sprintf("%2.6f (%2.6f)", est, mcse)) %>% 
  select(stat, performance, method) %>% 
  pivot_wider(names_from = method, values_from = c(performance))

table_30_MCAR


