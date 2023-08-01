
#####################################
# 
#    scenario: MNAR 30% 
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

bias_30_MNAR <- matrix(0, 1, 3)
empse_30_MNAR <- matrix(0, 1, 3)
modelse_30_MNAR <- matrix(0, 1, 3)
rmse_30_MNAR <- matrix(0, 1, 3)
cover_30_MNAR <- matrix(0, 1, 3)


N_sim <- 1000 # simulate 1000 times



# -------------------- 1st method: CCA --------------------

mnar_30MNAR_cc_prob <- c()

for (i in 1:N_sim){
  
  print(i) # know the current iteration
  
  
  merged_data <- merged_data_keep[sample(nrow(merged_data_keep), 3847, replace = TRUE), ]
  
  merged_data$logit_30 <- -0.59 * merged_data$event - 0.548
  merged_data$mnar_pro <- exp(merged_data$logit_30) / (1 + exp(merged_data$logit_30))
  merged_data$mnar_pro_missing <- rbinom(nrow(merged_data), 1, merged_data$mnar_pro)
  merged_data <- merged_data %>% 
    mutate(mnar_event = ifelse(mnar_pro_missing == 0, 
                               event, NA
    ))
  # check % missing
  sum(is.na(merged_data$mnar_event)) / nrow(merged_data)
  
  
  #### 1.1 CCA 
  overall_surv_MAR <- survfit(Surv(futime, as.numeric(mnar_event)) ~ 1, data = merged_data) 
  # estimated survival rates
  mnar_cc <- summary(overall_surv_MAR, times = 20)
  mnar_30_cc_prob_each_run <- as.data.frame(cbind(as.numeric(i), "CCA", mnar_cc$surv, mnar_cc$std.err,
                                                  mnar_cc$lower, mnar_cc$upper))
  mnar_30MNAR_cc_prob <- rbind(mnar_30MNAR_cc_prob, mnar_30_cc_prob_each_run)
}
# covert to numeric
for (i in 3:6){
  mnar_30MNAR_cc_prob[, i] <- as.numeric(mnar_30MNAR_cc_prob[, i])
}
colnames(mnar_30MNAR_cc_prob) <- c("dataset", "method", "estimate", "se", "CI_lower", "CI_upper")

##### performance measure
bias_30_MNAR[, 1] <- mean(mnar_30MNAR_cc_prob$estimate - true_CFR)
empse_30_MNAR[, 1] <- sqrt(sum( (mnar_30MNAR_cc_prob$estimate - true_CFR)^2 ) / (N_sim - 1))
modelse_30_MNAR[, 1] <- sqrt( sum( mnar_30MNAR_cc_prob$se^2 ) / N_sim)
rmse_30_MNAR[, 1] <- sqrt( sum((mnar_30MNAR_cc_prob$estimate - true_CFR)^2) / N_sim)
cover_30_MNAR[, 1] <- mean(mnar_30MNAR_cc_prob$CI_lower <= true_CFR & true_CFR <= mnar_30MNAR_cc_prob$CI_upper)


# -------------------- 2nd method: MI --------------------

mnar_30MNAR_mi_prob <- c()

for (i in 1:N_sim){
  
  print(i) # know the current iteration
  
  merged_data <- merged_data_keep[sample(nrow(merged_data_keep), 3847, replace = TRUE), ]
  
  merged_data$logit_30 <- -0.59 * merged_data$event - 0.548
  merged_data$mnar_pro <- exp(merged_data$logit_30) / (1 + exp(merged_data$logit_30))
  merged_data$mnar_pro_missing <- rbinom(nrow(merged_data), 1, merged_data$mnar_pro)
  merged_data <- merged_data %>% 
    mutate(mnar_event = ifelse(mnar_pro_missing == 0, 
                               event, NA
    ))
  # check % missing
  # sum(is.na(merged_data$mnar_event)) / nrow(merged_data)
  
  
  #### 1.2 MI
  # Here m = 31 because missingness = 30%
  # add Nelson-Aalen estimator
  
  #######################
  # imputation model
  #######################
  
  merged_data$na_estimator <- nelsonaalen(merged_data, futime, mnar_event)
  merged_data$mnar_event <- as.factor(merged_data$mnar_event) # make sure meth is logreg
  
  mi_MNAR_30 <- mice(merged_data  %>% 
                       select(futime, mnar_event, na_estimator, case_class,
                              siteid, sex, abdominal_pain, fever, hiccoughs, diarrhea, loss_of_appetite,
                              headache, bleeding, jaundice, difficulty_swallowing),
                     m = 31, maxit =20, printFlag=FALSE)
  
  # check if there is any NAs in the imputed dataset e.g. the 1st dataset
  # data_1 <- complete(mi_MAR_10 )
  # sum(is.na(data_1$mar_event))
  
  # MUST covert "mcar_event" to numeric rather than factor, otherwise summary will not gives "survival" results
  fit_mi_MNAR <- with(mi_MNAR_30, with(summary(survfit(Surv(futime, as.numeric(mnar_event)) ~ 1), times = 20, extend = TRUE),
                                       data.frame(time, surv, std.err, lower, upper)))
  
  pool_mi_MNAR <- fit_mi_MNAR$analyses %>% 
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
  
  mnar_30_mi_prob_each_run <- as.data.frame(cbind(as.numeric(i), "MI", pool_mi_MNAR$surv, pool_mi_MNAR$SE, 
                                                  pool_mi_MNAR$LCI, pool_mi_MNAR$UCI))
  mnar_30MNAR_mi_prob <- rbind(mnar_30MNAR_mi_prob, mnar_30_mi_prob_each_run)
}
# covert to numeric
for (i in 3:6){
  mnar_30MNAR_mi_prob[, i] <- as.numeric(mnar_30MNAR_mi_prob[, i])
}
colnames(mnar_30MNAR_mi_prob) <- c("dataset", "method", "estimate", "se", "CI_lower", "CI_upper")

##### performance measure
bias_30_MNAR[, 2] <- mean(mnar_30MNAR_mi_prob$estimate - true_CFR)
empse_30_MNAR[, 2] <- sqrt(sum( (mnar_30MNAR_mi_prob$estimate - true_CFR)^2 ) / (N_sim - 1))
modelse_30_MNAR[, 2] <- sqrt( sum( mnar_30MNAR_mi_prob$se^2 ) / N_sim)
rmse_30_MNAR[, 2] <- sqrt( sum((mnar_30MNAR_mi_prob$estimate - true_CFR)^2) / N_sim)
cover_30_MNAR[, 2] <- mean(mnar_30MNAR_mi_prob$CI_lower <= true_CFR & true_CFR <= mnar_30MNAR_mi_prob$CI_upper)




# -------------------- 3rd method: IPW --------------------

mnar_30MNAR_ipw_prob <- c()

for (i in 1:N_sim){
  
  print(i) # know the current iteration
  
  merged_data <- merged_data_keep[sample(nrow(merged_data_keep), 3847, replace = TRUE), ]
  
  merged_data$logit_30 <- -0.59 * merged_data$event - 0.548
  merged_data$mnar_pro <- exp(merged_data$logit_30) / (1 + exp(merged_data$logit_30))
  merged_data$mnar_pro_missing <- rbinom(nrow(merged_data), 1, merged_data$mnar_pro)
  merged_data <- merged_data %>% 
    mutate(mnar_event = ifelse(mnar_pro_missing == 0, 
                               event, NA
    ))
  # check % missing
  # sum(is.na(merged_data$mnar_event)) / nrow(merged_data)
  
  
  #### 1.3 IPW
  
  #######################
  # missingness model
  #######################
  
  ps_model_MNAR <- glm(is.na(mnar_event) ~ age_years  + sex
                       + abdominal_pain + fever + hiccoughs + diarrhea + loss_of_appetite
                       + headache + bleeding + jaundice + difficulty_swallowing, 
                       family = binomial, data = merged_data)
  # summary(ps_model)
  # score
  merged_data$propensity_score_MNAR <- predict(ps_model_MNAR, type = "response")
  # calculate weighting
  # estimate weight for each patient
  merged_data$weight.ATE_MNAR <- ifelse((is.na(merged_data$mnar_event)),
                                        (1/merged_data$propensity_score_MNAR), 
                                        (1/(1-merged_data$propensity_score_MNAR)))
  mnar_30_ipw_model <- survfit(Surv(futime, as.numeric(mnar_event)) ~ 1, weights = weight.ATE_MNAR, data = merged_data) 
  mnar_ipw <- summary(mnar_30_ipw_model, times = 20)
  
  # bootstrap to estimate SE
  sd_function <- function(data, indices = 1:nrow(data)) {
    d <- data %>% slice(indices) # allows boot to select sample
    ps_model_MNAR <- glm(is.na(mnar_event) ~ age_years + siteid + sex
                         + abdominal_pain + fever + hiccoughs + diarrhea + loss_of_appetite
                         + headache + bleeding + jaundice + difficulty_swallowing, 
                         family = binomial, data = d)
    d$propensity_score <- predict(ps_model_MNAR, type = "response")
    d$weight.ATE_MNAR <- ifelse((is.na(d$mnar_event)),
                                (1/d$propensity_score), 
                                (1/(1-d$propensity_score)))
    fit <- survfit(Surv(futime, as.numeric(mnar_event)) ~ 1, weights = weight.ATE_MNAR, data=d)
    return(summary(fit, times = 20)$surv)
  }
  resample_MNAR <- boot(merged_data, statistic = sd_function, R = 200)
  
  
  # store
  
  mnar_30_ipw_prob_each_run <- as.data.frame(cbind(as.numeric(i), "IPW", mnar_ipw$surv, summary(resample_MNAR)$bootSE,
                                                   mnar_ipw$surv - 1.96 * summary(resample_MNAR)$bootSE, 
                                                   mnar_ipw$surv + 1.96 * summary(resample_MNAR)$bootSE))
  mnar_30MNAR_ipw_prob <- rbind(mnar_30MNAR_ipw_prob, mnar_30_ipw_prob_each_run)
  
}
# covert to numeric
for (i in 3:6){
  mnar_30MNAR_ipw_prob[, i] <- as.numeric(mnar_30MNAR_ipw_prob[, i])
}
colnames(mnar_30MNAR_ipw_prob) <- c("dataset", "method", "estimate", "se", "CI_lower", "CI_upper")

##### performance measure
bias_30_MNAR[, 3] <- mean(mnar_30MNAR_ipw_prob$estimate - true_CFR)
empse_30_MNAR[, 3] <- sqrt(sum( (mnar_30MNAR_ipw_prob$estimate - true_CFR)^2 ) / (N_sim - 1))
modelse_30_MNAR[, 3] <- sqrt( sum( mnar_30MNAR_ipw_prob$se^2 ) / N_sim)
rmse_30_MNAR[, 3] <- sqrt( sum((mnar_30MNAR_ipw_prob$estimate - true_CFR)^2) / N_sim)
cover_30_MNAR[, 3] <- mean(mnar_30MNAR_ipw_prob$CI_lower <= true_CFR & true_CFR <= mnar_30MNAR_ipw_prob$CI_upper)



# -------------------- Results (Performance measures) --------------------

result_30_MNAR <- rbind(bias_30_MNAR, empse_30_MNAR, modelse_30_MNAR, rmse_30_MNAR, cover_30_MNAR)
rownames(result_30_MNAR) <- c("Bias", "Empirical SE", "Model SE", "RMSE", "Coverage")
colnames(result_30_MNAR) <- c("CCA", "MI", "IPW")
View(result_30_MNAR)


mnar_30_prob <- rbind(mnar_30MNAR_cc_prob, mnar_30MNAR_mi_prob, mnar_30MNAR_ipw_prob)
# write.csv(mnar_30_prob, "/Users/apple/Desktop/new_full_data/simulation/30_MNAR_simulation.csv")
library(rsimsum)
s1_prob_30_MNAR <- simsum(mnar_30_prob, estvarname = "estimate", true = true_CFR, se = "se", 
                          methodvar = "method", ref = "CCA")
ss1_prob_30_MNAR <- summary(s1_prob_30_MNAR, stats = c("bias", "modelse", "empse", "cover"))

table_30_MNAR <- ss1_prob_30_MNAR %>% 
  tidy() %>% 
  mutate(performance = sprintf("%2.6f (%2.6f)", est, mcse)) %>% 
  select(stat, performance, method) %>% 
  pivot_wider(names_from = method, values_from = c(performance))

table_30_MNAR

# Monte-calro error
s1_prob_30_MNAR$summ


# -------------------- 30% missingness Results (Combine MCAR, MAR, MNAR together) --------------------

####### result table
result_30 <- rbind(c("MCAR", "", "", ""), 
                   table_30_MCAR, 
                   c("MAR", "", "", ""),  
                   table_30_MAR,
                   c("MNAR", "", "", ""), 
                   table_30_MNAR)

library(flextable)
result_30 %>% flextable() %>% 
  width(j=c(2,3,4), width = 2) %>% 
  add_header_row(top = TRUE, values = c("Full Kplan-Meier estimate of day 20 (SE) = 0.4715 (0.014)", "", "", "")) %>% 
  merge_at(i = 1, j = 1:4, part = "header") %>% 
  align(align = "center", j = c(2:4), part = "all") %>% 
  bold(i = c(1, 6, 11), bold = TRUE, part = "body") %>% 
  flextable::save_as_docx(path = "/Users/apple/Desktop/new_full_data/simulation/table_30.docx")

# latex
xtable::xtable(result_30)



