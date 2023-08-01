
#####################################
# 
#    scenario: MAR 50% 
#
#####################################


###### load library
library(mice)
library(ipw)
library(dplyr)
library(tidyr)
library(boot)


# -------------------- simulation study --------------------

# step 0
# build a matrix to restore values for performance measures

bias_50_MAR <- matrix(0, 1, 3)
empse_50_MAR <- matrix(0, 1, 3)
modelse_50_MAR <- matrix(0, 1, 3)
rmse_50_MAR <- matrix(0, 1, 3)
cover_50_MAR <- matrix(0, 1, 3)

N_sim <- 1000


# -------------------- 1st method: CCA --------------------

mar_50MAR_cc_prob <- c()

set.seed(12345)

for (i in 1:N_sim){
  
  print(i) # know the current iteration
  
  merged_data <- merged_data_keep[sample(nrow(merged_data_keep), 3847, replace = TRUE), ]
  
  # MAR - 50%
  merged_data$logit <- 0.01077 + (-0.00216) * merged_data$age_years + 0.09796 * as.numeric(merged_data$sex)
  merged_data$mar_pro <- (exp(merged_data$logit)) / (1 + exp(merged_data$logit))
  merged_data$mar_pro_missing <- rbinom(nrow(merged_data), 1, merged_data$mar_pro)
  merged_data <- merged_data %>% 
    mutate(mar_event = ifelse(mar_pro_missing == 0, 
                              event, NA
    ))
  # check % missing
  sum(is.na(merged_data$mar_event)) / nrow(merged_data)
  
  
  #### 1.1 CCA 
  overall_surv_MAR <- survfit(Surv(futime, as.numeric(mar_event)) ~ 1, data = merged_data) 
  # estimated survival rates
  mar_cc <- summary(overall_surv_MAR, times = 20)
  mar_50_cc_prob_each_run <- as.data.frame(cbind(as.numeric(i), "CCA", mar_cc$surv, mar_cc$std.err,
                                                 mar_cc$lower, mar_cc$upper))
  mar_50MAR_cc_prob <- rbind(mar_50MAR_cc_prob, mar_50_cc_prob_each_run)
}
# covert to numeric
for (i in 3:6){
  mar_50MAR_cc_prob[, i] <- as.numeric(mar_50MAR_cc_prob[, i])
}
colnames(mar_50MAR_cc_prob) <- c("dataset", "method", "estimate", "se", "CI_lower", "CI_upper")

##### performance measure
bias_50_MAR[, 1] <- mean(mar_50MAR_cc_prob$estimate - true_CFR)
empse_50_MAR[, 1] <- sqrt(sum( (mar_50MAR_cc_prob$estimate - true_CFR)^2 ) / (N_sim - 1))
modelse_50_MAR[, 1] <- sqrt( sum( mar_50MAR_cc_prob$se^2 ) / N_sim)
rmse_50_MAR[, 1] <- sqrt( sum((mar_50MAR_cc_prob$estimate - true_CFR)^2) / N_sim)
cover_50_MAR[, 1] <- mean(mar_50MAR_cc_prob$CI_lower <= true_CFR & true_CFR <= mar_50MAR_cc_prob$CI_upper)


# -------------------- 2nd method: MI --------------------

mar_50MAR_mi_prob <- c()

for (i in 1:N_sim){
  
  print(i) # know the current iteration
  
  merged_data <- merged_data_keep[sample(nrow(merged_data_keep), 3847, replace = TRUE), ]
  
  # MAR - 50%
  merged_data$logit <- 0.01077 + (-0.00216) * merged_data$age_years + 0.09796 * as.numeric(merged_data$sex)
  merged_data$mar_pro <- (exp(merged_data$logit)) / (1 + exp(merged_data$logit))
  merged_data$mar_pro_missing <- rbinom(nrow(merged_data), 1, merged_data$mar_pro)
  merged_data <- merged_data %>% 
    mutate(mar_event = ifelse(mar_pro_missing == 0, 
                              event, NA
    ))
  # check % missing
  # sum(is.na(merged_data$mar_event)) / nrow(merged_data)
  
  
  #### 1.2 MI
  # Here m = 51 because missingness = 50%
  # add Nelson-Aalen estimator
  
  #######################
  # imputation model
  #######################
  
  merged_data$na_estimator <- nelsonaalen(merged_data, futime, mar_event)
  merged_data$mar_event <- as.factor(merged_data$mar_event) # make sure meth is logreg
  
  mi_MAR_50 <- mice(merged_data  %>% 
                      select(futime, mar_event, na_estimator, case_class,
                             siteid, sex, abdominal_pain, fever, hiccoughs, diarrhea, loss_of_appetite,
                             headache, bleeding, jaundice, difficulty_swallowing),
                    m = 51, maxit =20, printFlag=FALSE)
  
  # check if there is any NAs in the imputed dataset e.g. the 1st dataset
  # data_1 <- complete(mi_MAR_10 )
  # sum(is.na(data_1$mar_event))
  
  # MUST covert "mcar_event" to numeric rather than factor, otherwise summary will not gives "survival" results
  fit_mi_MAR <- with(mi_MAR_50, with(summary(survfit(Surv(futime, as.numeric(mar_event)) ~ 1), times = 20, extend = TRUE),
                                     data.frame(time, surv, std.err, lower, upper)))
  
  pool_mi_MAR <- fit_mi_MAR$analyses %>% 
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
  
  mar_50_mi_prob_each_run <- as.data.frame(cbind(as.numeric(i), "MI", pool_mi_MAR$surv, pool_mi_MAR$SE, 
                                                 pool_mi_MAR$LCI, pool_mi_MAR$UCI))
  mar_50MAR_mi_prob <- rbind(mar_50MAR_mi_prob, mar_50_mi_prob_each_run)
}
# covert to numeric
for (i in 3:6){
  mar_50MAR_mi_prob[, i] <- as.numeric(mar_50MAR_mi_prob[, i])
}
colnames(mar_50MAR_mi_prob) <- c("dataset", "method", "estimate", "se", "CI_lower", "CI_upper")

##### performance measure
bias_50_MAR[, 2] <- mean(mar_50MAR_mi_prob$estimate - true_CFR)
empse_50_MAR[, 2] <- sqrt(sum( (mar_50MAR_mi_prob$estimate - true_CFR)^2 ) / (N_sim - 1))
modelse_50_MAR[, 2] <- sqrt( sum( mar_50MAR_mi_prob$se^2 ) / N_sim)
rmse_50_MAR[, 2] <- sqrt( sum((mar_50MAR_mi_prob$estimate - true_CFR)^2) / N_sim)
cover_50_MAR[, 2] <- mean(mar_50MAR_mi_prob$CI_lower <= true_CFR & true_CFR <= mar_50MAR_mi_prob$CI_upper)




# -------------------- 3rd method: IPW --------------------

mar_50MAR_ipw_prob <- c()

for (i in 1:N_sim){
  
  print(i) # know the current iteration
  
  merged_data <- merged_data_keep[sample(nrow(merged_data_keep), 3847, replace = TRUE), ]
  
  # MAR - 50%
  merged_data$logit <- 0.01077 + (-0.00216) * merged_data$age_years + 0.09796 * as.numeric(merged_data$sex)
  merged_data$mar_pro <- (exp(merged_data$logit)) / (1 + exp(merged_data$logit))
  merged_data$mar_pro_missing <- rbinom(nrow(merged_data), 1, merged_data$mar_pro)
  merged_data <- merged_data %>% 
    mutate(mar_event = ifelse(mar_pro_missing == 0, 
                              event, NA
    ))
  # check % missing
  # sum(is.na(merged_data$mar_event)) / nrow(merged_data)
  
  
  #### 1.3 IPW
  
  #######################
  # missingness model
  #######################
  
  ps_model_MAR <- glm(is.na(mar_event) ~ age_years  + sex
                      + abdominal_pain + fever + hiccoughs + diarrhea + loss_of_appetite
                      + headache + bleeding + jaundice + difficulty_swallowing, 
                      family = binomial, data = merged_data)
  # summary(ps_model)
  # score
  merged_data$propensity_score_MAR <- predict(ps_model_MAR, type = "response")
  # calculate weighting
  # estimate weight for each patient
  merged_data$weight.ATE_MAR <- ifelse((is.na(merged_data$mar_event)),
                                       (1/merged_data$propensity_score_MAR), 
                                       (1/(1-merged_data$propensity_score_MAR)))
  mar_50_ipw_model <- survfit(Surv(futime, as.numeric(mar_event)) ~ 1, weights = weight.ATE_MAR, data = merged_data) 
  mar_ipw <- summary(mar_50_ipw_model, times = 20)
  
  # bootstrap to estimate SE
  sd_function <- function(data, indices = 1:nrow(data)) {
    d <- data %>% slice(indices) # allows boot to select sample
    ps_model_MAR <- glm(is.na(mar_event) ~ age_years + siteid + sex
                        + abdominal_pain + fever + hiccoughs + diarrhea + loss_of_appetite
                        + headache + bleeding + jaundice + difficulty_swallowing, 
                        family = binomial, data = d)
    d$propensity_score <- predict(ps_model_MAR, type = "response")
    d$weight.ATE_MAR <- ifelse((is.na(d$mar_event)),
                               (1/d$propensity_score), 
                               (1/(1-d$propensity_score)))
    fit <- survfit(Surv(futime, as.numeric(mar_event)) ~ 1, weights = weight.ATE_MAR, data=d)
    return(summary(fit, times = 20)$surv)
  }
  resample_MAR <- boot(merged_data, statistic = sd_function, R = 200)
  
  
  # store
  
  mar_50_ipw_prob_each_run <- as.data.frame(cbind(as.numeric(i), "IPW", mar_ipw$surv, summary(resample_MAR)$bootSE,
                                                  mar_ipw$surv - 1.96 * summary(resample_MAR)$bootSE, 
                                                  mar_ipw$surv + 1.96 * summary(resample_MAR)$bootSE))
  mar_50MAR_ipw_prob <- rbind(mar_50MAR_ipw_prob, mar_50_ipw_prob_each_run)
}
# covert to numeric
for (i in 3:6){
  mar_50MAR_ipw_prob[, i] <- as.numeric(mar_50MAR_ipw_prob[, i])
}
colnames(mar_50MAR_ipw_prob) <- c("dataset", "method", "estimate", "se", "CI_lower", "CI_upper")

##### performance measure
bias_50_MAR[, 3] <- mean(mar_50MAR_ipw_prob$estimate - true_CFR)
empse_50_MAR[, 3] <- sqrt(sum( (mar_50MAR_ipw_prob$estimate - true_CFR)^2 ) / (N_sim - 1))
modelse_50_MAR[, 3] <- sqrt( sum( mar_50MAR_ipw_prob$se^2 ) / N_sim)
rmse_50_MAR[, 3] <- sqrt( sum((mar_50MAR_ipw_prob$estimate - true_CFR)^2) / N_sim)
cover_50_MAR[, 3] <- mean(mar_50MAR_ipw_prob$CI_lower <= true_CFR & true_CFR <= mar_50MAR_ipw_prob$CI_upper)



# -------------------- Results (Performance measures) --------------------

result_50_MAR <- rbind(bias_50_MAR, empse_50_MAR, modelse_50_MAR, rmse_50_MAR, cover_50_MAR)
rownames(result_50_MAR) <- c("Bias", "Empirical SE", "Model SE", "RMSE", "Coverage")
colnames(result_50_MAR) <- c("CCA", "MI", "IPW")
View(result_50_MAR)


mar_50_prob_MAR <- rbind(mar_50MAR_cc_prob, mar_50MAR_mi_prob, mar_50MAR_ipw_prob)
# write.csv(mar_50_prob_MAR, "/Users/apple/Desktop/new_full_data/simulation/50_MAR_simulation.csv")
library(rsimsum)
s1_prob_50_MAR <- simsum(mar_50_prob_MAR, estvarname = "estimate", true = true_CFR, se = "se", 
                         methodvar = "method", ref = "CCA")
ss1_prob_50_MAR <- summary(s1_prob_50_MAR, stats = c("bias", "modelse", "empse", "cover"))

table_50_MAR <- ss1_prob_50_MAR %>% 
  tidy() %>% 
  mutate(performance = sprintf("%2.6f (%2.6f)", est, mcse)) %>% 
  select(stat, performance, method) %>% 
  pivot_wider(names_from = method, values_from = c(performance))

table_50_MAR

# Monte-calro error
s1_prob_50_MAR$summ




