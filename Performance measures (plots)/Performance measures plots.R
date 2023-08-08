
########################################################################
# !!!!!!!
#
# PLEASE run all codes from "simulation" files before running this file
# as it requires results from simulation models
#
# !!!!!!!
########################################################################


# plot 
table_plot <- matrix(0, nrow = 9, nco = 6)
table_plot <- as_tibble(as.data.frame(table_plot))
colnames(table_plot) <- c("missing",  "")


########################
# statistics for day 20
########################

# 10% MCAR 
table_1_plot <- as.data.frame(ss1_prob_MCAR %>% 
                                tidy() %>% 
                                select(stat, est, method) %>% 
                                pivot_wider(names_from = method, values_from = c(est)))

# 10% MAR
table_2_plot <- as.data.frame(ss1_prob_MAR %>% 
                                tidy() %>% 
                                select(stat, est, method) %>% 
                                pivot_wider(names_from = method, values_from = c(est)))

# 10% MCAR
table_3_plot <- as.data.frame(ss1_prob_MNAR %>% 
                                tidy() %>% 
                                select(stat, est, method) %>% 
                                pivot_wider(names_from = method, values_from = c(est)))

# 30% MCAR
table_30_MCAR_plot <- as.data.frame(ss1_prob_30_MCAR %>% 
                                      tidy() %>% 
                                      select(stat, est, method) %>% 
                                      pivot_wider(names_from = method, values_from = c(est))
)

# 30% MAR
table_30_MAR_plot <- as.data.frame(ss1_prob_30_MAR %>% 
                                     tidy() %>% 
                                     select(stat, est, method) %>% 
                                     pivot_wider(names_from = method, values_from = c(est)))

# 30% MNAR
table_30_MNAR_plot <- as.data.frame(ss1_prob_30_MNAR %>% 
                                      tidy() %>% 
                                      select(stat, est, method) %>% 
                                      pivot_wider(names_from = method, values_from = c(est)))

# 50% MCAR
table_50_MCAR_plot <- as.data.frame(ss1_prob_50_MCAR %>% 
                                      tidy() %>% 
                                      select(stat, est, method) %>% 
                                      pivot_wider(names_from = method, values_from = c(est)))

# 50% MAR
table_50_MAR_plot <- as.data.frame(ss1_prob_50_MAR %>% 
                                     tidy() %>% 
                                     select(stat, est, method) %>% 
                                     pivot_wider(names_from = method, values_from = c(est)))

# 50% MNAR
table_50_MNAR_plot <- as.data.frame(ss1_prob_50_MNAR %>% 
                                      tidy() %>% 
                                      select(stat, est, method) %>% 
                                      pivot_wider(names_from = method, values_from = c(est)))



# ------------------------------------
#        coverage
# ------------------------------------

table_plot_cover <- data.frame(
  case_class = c(rep("Probable", 9), rep("Confirmed", 9)),
  perc = c(rep(0.1, 3), rep(0.3, 3), rep(0.5, 3), rep(0.1, 3), rep(0.3, 3), rep(0.5, 3)),
  mech = rep(c("MCAR", "MAR", "MNAR"), 6),
  CCA = c(### day 20
    # 10%
    table_1_plot[4,2], table_2_plot[4,2], table_3_plot[4,2],
    # 0.956, 0.954, 0.812
    
    # 30%
    table_30_MCAR_plot[4,2], table_30_MAR_plot[4,2], table_30_MNAR_plot[4,2], 
    # 0.953, 0.952, 0.22
    
    # 50%
    table_50_MCAR_plot[4,2], table_50_MAR_plot[4,2], table_50_MNAR_plot[4,2],
    # 0.964, 0.941, 0
    
    ### confirmed
    # 10%
    0.953, 0.955, 0.904,
    
    # 30%
    0.948000, 0.942, 0.558,
    
    # 50%
    0.947, 0.948, 0.0
  ),
  
  IPW = c(# 10%
    table_1_plot[4,3], table_2_plot[4,3], table_3_plot[4,3],
    # 0.949, 0.946, 0.804
    
    # 30%
    table_30_MCAR_plot[4, 3], table_30_MAR_plot[4,3], table_30_MNAR_plot[4,3], 
    # 0.945, 0.95, 0.279
    
    # 50%
    table_50_MCAR_plot[4,3], table_50_MAR_plot[4,3], table_50_MNAR_plot[4,3],
    # 0.942, 0.959, 0
    
    ### confirmed
    # 10%
    0.966000, 0.949000, 0.919,
    
    # 30%
    0.952000, 0.941, 0.689,
    
    # 50%
    0.946, 0.928, 0.01
  ),
  
  MI = c(# 10%
    table_1_plot[4,4], table_2_plot[4,4], table_3_plot[4,4],
    # 0.949, 0.953, 0.889
    
    # 30%
    table_30_MCAR_plot[4,4], table_30_MAR_plot[4,4], table_30_MNAR_plot[4,4], 
    # 0.95, 0.949, 0.462
    
    # 50%
    table_50_MCAR_plot[4,4], table_50_MAR_plot[4,4], table_50_MNAR_plot[4,4],
    # 0.942, 0.947, 0
    
    ### confirmed
    # 10%
    0.944000, 0.952000, 0.952,
    
    # 30%
    0.930000, 0.937, 0.847,
    
    # 50%
    0.94, 0.951, 0.156
  )
)


cover_long <- table_plot_cover %>% 
  pivot_longer(cols = CCA:MI,
               names_to = "method",
               values_to = "cover")
cover_long$mech <- factor(cover_long$mech,      # Reordering group factor levels
                          levels = c("MCAR", "MAR", "MNAR"))


cover_plot <- ggplot(cover_long, aes(x = perc, y = cover, fill = mech)) + 
  scale_y_continuous(breaks = c(0, 0.25,  0.5, 0.75,  0.95, 1)) +
  geom_point(aes(colour = method), size = 3) + 
  geom_line(aes(colour = method)) + 
  facet_wrap(case_class ~ mech, ncol = 3) +
  geom_hline(yintercept = 0.95, linetype="dashed", color = "black")

print(cover_plot + labs(y = "Coverage probability", x = "Amount of missing outcome (%)",
                        colour = "Methods")
      + theme(legend.position="bottom", 
              text = element_text(size = 18),
              axis.title.y = element_text(size = 24),
              axis.title.x = element_text(size = 24)
      ) 
      + guides(fill = "none"))





# --------------------------------------------------
#        Model based SE (cloglog transformed KM)
# --------------------------------------------------

table_plot_modelSE <- data.frame(
  case_class = c(rep("Probable", 9), rep("Confirmed", 9)),
  perc = c(rep(0.1, 3), rep(0.3, 3), rep(0.5, 3), rep(0.1, 3), rep(0.3, 3), rep(0.5, 3)),
  mech = rep(c("MCAR", "MAR", "MNAR"), 6),
  
  CCA = c(# 10%
    table_1_plot[3,2], table_2_plot[3,2], table_3_plot[3,2],
    # 0.01476462, 0.01469668, 0.01483343
    
    # 30%
    table_30_MCAR_plot[3,2], table_30_MAR_plot[3,2], table_30_MNAR_plot[3,2], 
    # 0.01676084, 0.01680263, 0.01692241
    
    # 50%
    table_50_MCAR_plot[3,2], table_50_MAR_plot[3,2], table_50_MNAR_plot[3,2],
    # 0.01987096,  0.02010105, 0.01928228
    
    ### confirmed
    # 10%
    0.02070602, 0.020542, 0.020569,
    
    # 30%
    0.023484, 0.023405, 0.024462,
    
    # 50%
    0.0277294, 0.028014, 0.030993
  ),
  
  IPW = c(# 10%
    table_1_plot[3,3], table_2_plot[3,3], table_3_plot[3,3],
    # 0.0147869, 0.01470773, 0.01487673,
    
    # 30%
    table_30_MCAR_plot[3, 3], table_30_MAR_plot[3,3], table_30_MNAR_plot[3,3], 
    # 0.01487673, 0.01674784, 0.01686606,
    
    # 50%
    table_50_MCAR_plot[3,3], table_50_MAR_plot[3,3], table_50_MNAR_plot[3,3],
    # 0.01965497, 0.02000705, 0.01941616,
    
    ### confirmed
    # 10%
    0.020416, 0.020322, 0.020289,
    
    # 30%
    0.022427, 0.022504, 0.023430,
    
    # 50%
    0.025958, 0.026299, 0.031387
  ),
  
  MI = c(# 10%
    table_1_plot[3,4], table_2_plot[3,4], table_3_plot[3,4],
    # 0.01432169, 0.01430691, 0.01433951,
    
    # 30%
    table_30_MCAR_plot[3,4], table_30_MAR_plot[3,4], table_30_MNAR_plot[3,4], 
    # 0.01526228, 0.01522502, 0.01537497,
    
    # 50%
    table_50_MCAR_plot[3,4], table_50_MAR_plot[3,4], table_50_MNAR_plot[3,4],
    # 0.01689866, 0.01701797, 0.0172508,
    
    ### confirmed
    # 10%
    0.01992211, 0.019911, 0.019813,
    
    # 30%
    0.020910, 0.020961, 0.021418,
    
    # 50%
    0.022529, 0.022909, 0.025078
  )
)


modelSE_long <- table_plot_modelSE %>% 
  pivot_longer(cols = CCA:MI,
               names_to = "method",
               values_to = "modelSE")
modelSE_long$mech <- factor(modelSE_long$mech,      # Reordering group factor levels
                            levels = c("MCAR", "MAR", "MNAR"))

modelSE_plot <- ggplot(modelSE_long, aes(x = perc, y = modelSE, fill = mech)) + 
  ylim(0.01, 0.035) + 
  geom_point(aes(colour = method), size = 3) + 
  geom_line(aes(colour = method)) +
  facet_wrap(case_class ~ mech, ncol = 3) +  
  geom_hline(aes(yintercept = ifelse(case_class == "Confirmed", 0.02, 0.01401776)),
             linetype = "dashed", color = "black")

print(modelSE_plot + labs(y = "Model based SE (cloglog transformed KM)", x = "Amount of missing outcome (%)",
                          colour = "Methods")
      + theme(legend.position="bottom", 
              text = element_text(size = 18),
              axis.title.y = element_text(size = 24),
              axis.title.x = element_text(size = 24) ) + guides(fill = "none"))




# --------------------------------------------------
#        Empirical SE (cloglog transformed KM)
# --------------------------------------------------

table_plot_empSE <- data.frame(
  case_class = c(rep("Probable", 9), rep("Confirmed", 9)),
  perc = c(rep(0.1, 3), rep(0.3, 3), rep(0.5, 3), rep(0.1, 3), rep(0.3, 3), rep(0.5, 3)),
  mech = rep(c("MCAR", "MAR", "MNAR"), 6),
  
  CCA = c(# 10%
    table_1_plot[2,2], table_2_plot[2,2], table_3_plot[2,2],
    # 0.01434095, 0.01422764, 0.01438965, 
    
    # 30%
    table_30_MCAR_plot[2,2], table_30_MAR_plot[2,2], table_30_MNAR_plot[2,2], 
    # 0.01638744, 0.01637468, 0.01677075,
    
    # 50%
    table_50_MCAR_plot[2,2], table_50_MAR_plot[2,2], table_50_MNAR_plot[2,2],
    # 0.01871179, 0.0202628, 0.01965434,
    
    ### confirmed
    # 10%
    0.020537, 0.020408, 0.020224,
    
    # 30%
    0.023541, 0.023588, 0.024316,
    
    # 50%
    0.02870541, 0.028284, 0.030670
  ),
  
  IPW = c(# 10%
    table_1_plot[2,3], table_2_plot[2,3], table_3_plot[2,3],
    # 0.01487106, 0.01499, 0.01464055,
    
    # 30%
    table_30_MCAR_plot[2, 3], table_30_MAR_plot[2,3], table_30_MNAR_plot[2,3],
    # 0.01674733, 0.01647144, 0.01718611,
    
    # 50%
    table_50_MCAR_plot[2,3], table_50_MAR_plot[2,3], table_50_MNAR_plot[2,3],
    # 0.02015979, 0.01934738, 0.01818377,
    
    ### confirmed
    # 10%
    0.020351, 0.020282, 0.019671,
    
    # 30%
    0.022063, 0.023288, 0.022849,
    
    # 50%
    0.026021, 0.029111, 0.031178
  ),
  
  MI = c(# 10%
    table_1_plot[2,4], table_2_plot[2,4], table_3_plot[2,4],
    # 0.01429439, 0.01399926, 0.01485063,
    
    # 30%
    table_30_MCAR_plot[2,4], table_30_MAR_plot[2,4], table_30_MNAR_plot[2,4], 
    # 0.01539901, 0.0153573, 0.01591203,
    
    # 50%
    table_50_MCAR_plot[2,4], table_50_MAR_plot[2,4], table_50_MNAR_plot[2,4],
    # 0.01708636, 0.01764364, 0.01704834,
    
    ### confirmed
    # 10%
    0.020456, 0.020032, 0.020868,
    
    # 30%
    0.021615, 0.020892, 0.020223,
    
    # 50%
    0.021978, 0.021965, 0.021746
  )
)


empSE_long <- table_plot_empSE %>% 
  pivot_longer(cols = CCA:MI,
               names_to = "method",
               values_to = "empSE")
empSE_long$mech <- factor(empSE_long$mech,      # Reordering group factor levels
                          levels = c("MCAR", "MAR", "MNAR"))

empSE_plot <- ggplot(empSE_long, aes(x = perc, y = empSE, fill = mech)) + 
  ylim(0.01, 0.035) + 
  geom_point(aes(colour = method), size = 3) + 
  geom_line(aes(colour = method)) +
  facet_wrap(case_class ~ mech, ncol = 3)  + 
  geom_hline(aes(yintercept = ifelse(case_class == "Confirmed", 0.02, 0.01401776)),
             linetype = "dashed", color = "black")

print(empSE_plot + labs(y = "Empirical SE (cloglog transformed KM)", x = "Amount of missing outcome (%)",
                        colour = "Methods")
      + theme(legend.position="bottom", 
              text = element_text(size = 18),
              axis.title.y = element_text(size = 24),
              axis.title.x = element_text(size = 24)) + guides(fill = "none"))





# --------------------------------------------------
#        Estimated K-M
# --------------------------------------------------


table_plot_KM <- data.frame(
  case_class = c(rep("Probable", 9), rep("Confirmed", 9)),
  perc = c(rep(0.1, 3), rep(0.3, 3), rep(0.5, 3), rep(0.1, 3), rep(0.3, 3), rep(0.5, 3)),
  mech = rep(c("MCAR", "MAR", "MNAR"), 6),
  
  CCA = c(# 10%
    1 - mean(mcar_10_cc_prob$estimate), 1 - mean(mar_10_cc_prob$estimate), 1 - mean(mnar_10_cc_prob$estimate),
    # 1- 0.4719188, 1-0.4716906, 1-0.4542880,
    
    # 30%
    1 - mean(mcar_30_cc_prob$estimate), 1 - mean(mar_30MAR_cc_prob$estimate), 1 - mean(mnar_30MNAR_cc_prob$estimate), 
    # 1-0.4719221, 1-0.4715635, 0.5173817
    
    # 50%
    1 - mean(mcar_50_cc_prob$estimate), 1 - mean(mar_50MAR_cc_prob$estimate), 1 - mean(mnar_50MNAR_cc_prob$estimate),
    # 1-0.4714922, 1-0.4690354, 1-0.6640603
    
    ### confirmed --> true value = 
    # 10%
    1-0.3543166, 1-0.3548896, 1-0.3394143,
    
    # 30%
    1-0.3538898, 1-0.3560488, 1-0.3976578,
    
    # 50%
    1-0.3542616, 1-0.3549039, 1-0.5457415
  ),
  
  IPW = c(# 10%
    1 - mean(mcar_10_ipw_prob$estimate), 1 - mean(mar_10_ipw_prob$estimate), 1 - mean(mnar_10_ipw_prob$estimate),
    # 1-0.4713313, 1-0.4713794, 1-0.4546900,
    
    # 30%
    1 - mean(mcar_30_ipw_prob$estimate), 1 - mean(mar_30MAR_ipw_prob$estimate), 1 - mean(mnar_30MNAR_ipw_prob$estimate), 
    # 1-0.4711185, 1-0.4721037, 1-0.5158019,
    
    # 50%
    1 - mean(mcar_50_ipw_prob$estimate), 1 - mean(mar_50MAR_ipw_prob$estimate), 1 - mean(mnar_50MNAR_ipw_prob$estimate),
    # 1-0.4722161, 1-0.4708730, 1-0.6581585,
    
    ### confirmed
    # 10%
    1-0.3532312, 1-0.3538981, 1-0.3418212,
    
    # 30%
    1-0.3547473, 1-0.353879, 1-0.388213,
    
    # 50%
    1-0.3543936, 1-0.3536793, 1-0.5083464
  ),
  
  MI = c(# 10%
    1 - mean(mcar_10_mi_prob$estimate), 1 - mean(mar_10_mi_prob$estimate), 1 - mean(mnar_10_mi_prob$estimate),
    # 1-0.4722314, 1-0.4714669, 1-0.4610317,
    
    # 30%
    1 - mean(mcar_30_mi_prob$estimate), 1 - mean(mar_30MAR_mi_prob$estimate), 1 - mean(mnar_30MNAR_mi_prob$estimate), 
    # 1-0.4725133, 1-0.4717906, 1-0.5024445,
    
    # 50%
    1 - mean(mcar_50_mi_prob$estimate), 1 - mean(mar_50MAR_mi_prob$estimate), 1 - mean(mnar_50MNAR_mi_prob$estimate),
    # 1-0.4722969, 1-0.4717429, 1-0.6016433,
    
    ### confirmed
    # 10%
    1-0.3552033, 1-0.355793, 1-0.3501859,
    
    # 30%
    1-0.3581416, 1-0.3581491, 1-0.3741252,
    
    # 50%
    1-0.3607706, 1-0.3572879, 1-0.4267415
  )
)


KM_long <- table_plot_KM %>% 
  pivot_longer(cols = CCA:MI,
               names_to = "method",
               values_to = "KM") %>% 
  group_by(case_class) 
KM_long$mech <- factor(KM_long$mech,      # Reordering group factor levels
                       levels = c("MCAR", "MAR", "MNAR"))

KM_plot <- ggplot(KM_long, aes(x = perc, y = KM, fill = mech)) +
  geom_point(aes(colour = method), size = 3) + 
  geom_line(aes(colour = method)) +
  facet_wrap(case_class ~ mech, ncol = 3)  + 
  geom_hline(aes(yintercept = ifelse( case_class == "Confirmed", 
                                      round(1-0.3540003, 2), round(1-0.4715159, 2))), linetype = "dashed", color = "black") +
  scale_y_continuous(breaks = c(0.3, 0.4, 0.5, 0.53, 0.6, 0.65, 0.7))

print(KM_plot + labs(y = "Estimated K-M", x = "Amount of missing outcome (%)",
                     colour = "Methods")
      + theme(legend.position="bottom", 
              text = element_text(size = 18),
              axis.title.y = element_text(size = 24),
              axis.title.x = element_text(size = 24)) + guides(fill = "none"))





# --------------------------------------------------
#        RMSE
# --------------------------------------------------


table_plot_rmse <- data.frame(
  case_class = c(rep("Probable", 9), rep("Confirmed", 9)),
  perc = c(rep(0.1, 3), rep(0.3, 3), rep(0.5, 3), rep(0.1, 3), rep(0.3, 3), rep(0.5, 3)),
  mech = rep(c("MCAR", "MAR", "MNAR"), 6),
  
  CCA = c(# 10%
    sqrt( sum( (mcar_10_cc_prob$estimate -  true_CFR)^2 ) / 1000),
    sqrt( sum( (mar_10_cc_prob$estimate -  true_CFR)^2 ) / 1000),
    sqrt( sum( (mnar_10_cc_prob$estimate -  true_CFR)^2 ) / 1000),
    # 0.01433943, 0.01422160, 0.02244227,
    
    # 30%
    sqrt( sum( (mcar_30_cc_prob$estimate -  true_CFR)^2 ) / 1000),
    sqrt( sum( (mar_30MAR_cc_prob$estimate -  true_CFR)^2 ) / 1000),
    sqrt( sum( (mnar_30MNAR_cc_prob$estimate -  true_CFR)^2 ) / 1000),
    # 0.01638428, 0.01636656, 0.04883281,
    
    # 50%
    sqrt( sum( (mcar_50_cc_prob$estimate -  true_CFR)^2 ) / 1000),
    sqrt( sum( (mar_50MAR_cc_prob$estimate -  true_CFR)^2 ) / 1000),
    sqrt( sum( (mnar_50MNAR_cc_prob$estimate -  true_CFR)^2 ) / 1000),
    # 0.01870245, 0.02040401, 0.19354391,
    
    ### confirmed
    # 10%
    0.02052897, 0.02041729, 0.02492676,
    
    # 30%
    0.02352983, 0.02366479, 0.04996661,
    
    # 50%
    0.02869106, 0.02828441, 0.1941761
  ),
  
  IPW = c(# 10%
    sqrt( sum( (mcar_10_ipw_prob$estimate -  true_CFR)^2 ) / 1000),
    sqrt( sum( (mar_10_ipw_prob$estimate -  true_CFR)^2 ) / 1000),
    sqrt( sum( (mnar_10_ipw_prob$estimate -  true_CFR)^2 ) / 1000),
    # 0.01486477, 0.01498312, 0.02229898,
    
    # 30%
    sqrt( sum( (mcar_30_ipw_prob$estimate -  true_CFR)^2 ) / 1000),
    sqrt( sum( (mar_30MAR_ipw_prob$estimate -  true_CFR)^2 ) / 1000),
    sqrt( sum( (mnar_30MNAR_ipw_prob$estimate -  true_CFR)^2 ) / 1000),
    # 0.01674368, 0.01647369, 0.04750067,
    
    # 50%
    sqrt( sum( (mcar_50_ipw_prob$estimate -  true_CFR)^2 ) / 1000),
    sqrt( sum( (mar_50MAR_ipw_prob$estimate -  true_CFR)^2 ) / 1000),
    sqrt( sum( (mnar_50MNAR_ipw_prob$estimate -  true_CFR)^2 ) / 1000),
    # 0.02016187m, 0.01934839, 0.18752536,
    
    ### confirmed
    # 10%
    0.02035569, 0.02027176, 0.02312786,
    
    # 30%
    0.02206496, 0.02327645, 0.04113483,
    
    # 50%
    0.02601083, 0.02909785, 0.1574605
  ),
  
  MI = c(# 10%
    sqrt( sum( (mcar_10_mi_prob$estimate -  true_CFR)^2 ) / 1000),
    sqrt( sum( (mar_10_mi_prob$estimate -  true_CFR)^2 ) / 1000),
    sqrt( sum( (mnar_10_mi_prob$estimate -  true_CFR)^2 ) / 1000),
    # 0.01430515, 0.01399234, 0.01817250,
    
    # 30%
    sqrt( sum( (mcar_30_mi_prob$estimate -  true_CFR)^2 ) / 1000),
    sqrt( sum( (mar_30MAR_mi_prob$estimate -  true_CFR)^2 ) / 1000),
    sqrt( sum( (mnar_30MNAR_mi_prob$estimate -  true_CFR)^2 ) / 1000),
    # 0.01542360, 0.01535208, 0.03477813,
    
    # 50%
    sqrt( sum( (mcar_50_mi_prob$estimate -  true_CFR)^2 ) / 1000),
    sqrt( sum( (mar_50MAR_mi_prob$estimate -  true_CFR)^2 ) / 1000),
    sqrt( sum( (mnar_50MNAR_mi_prob$estimate -  true_CFR)^2 ) / 1000),
    # 0.01709566, 0.01763627, 0.13123828,
    
    ### confirmed
    # 10%
    0.02048114, 0.02010159, 0.02120317,
    
    # 30%
    0.02199742, 0.02128962, 0.02852289,
    
    # 50%
    0.02298624, 0.02219866, 0.07591889
  )
)


rmse_long <- table_plot_rmse %>% 
  pivot_longer(cols = CCA:MI,
               names_to = "method",
               values_to = "rmse")
rmse_long$mech <- factor(rmse_long$mech,      # Reordering group factor levels
                         levels = c("MCAR", "MAR", "MNAR"))


rmse_plot <- ggplot(rmse_long, aes(x = perc, y = rmse, fill = mech)) + 
  scale_y_continuous(breaks = c(0, 0.05, 0.1, 0.15, 0.2)) +
  geom_point(aes(colour = method), size = 3) + 
  geom_line(aes(colour = method)) +
  facet_wrap(case_class ~ mech, ncol = 3)

print(rmse_plot + labs(y = "Root mean squared error (RMSE)", x = "Amount of missing outcome (%)",
                       colour = "Methods")
      + theme(legend.position="bottom", 
              text = element_text(size = 18),
              axis.title.y = element_text(size = 24),
              axis.title.x = element_text(size = 24)) + guides(fill = "none"))

