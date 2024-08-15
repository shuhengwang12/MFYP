main_theme_combined <- theme_minimal() +
  theme(
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(), 
    axis.ticks = element_line(color = "black", linewidth = 1),
    axis.ticks.length = unit(0.1, "cm"),
    axis.text = element_text(family = "sans", face = "bold", size = 10, color = "black"),
    axis.title = element_text(family = "sans", face = "bold", size = 11, color = "black"),
    legend.text = element_text(family = "sans", face = "bold", size = 11, color = "black"),
    legend.title = element_text(family = "sans", face = "bold", size = 11, color = "black"),
    legend.position = "right",
    plot.background = element_rect(fill = "white", color = NA),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 2))

Q_One <- data.frame(
  ResponseVariable = c("Com", "Amo", "Chl", "Dif", "Gly", "Imi", "Met", "Oxy", "Teb"),
  OOB = c(0.00386,
              0.2973,
              0.23552,
              0.27027,
              0.31274,
              0.28958,
              0.3166,
              0.02317,
              0.26255),
  AUC = c(1,
          0.67273,
          0.72321,
          0.55128,
          0.53957,
          0.64312,
          0.58302,
          0.98139,
          0.62649),
  Accuracy = c (0.97196,
                0.7757,
                0.7196,
                0.729,
                0.7383,
                0.785,
                0.7664,
                0.9813,
                0.6822)
)
Q_One_long <- Q_One %>%
  pivot_longer(cols = c(OOB, AUC, Accuracy), names_to = "Metric", values_to = "Value")
Q_One_long$ResponseVariable <- factor(Q_One_long$ResponseVariable, levels = c("Com", "Amo", "Chl", "Dif", "Gly", "Imi", "Met", "Oxy", "Teb"))
Q_one_Fig <- ggplot(Q_One_long, aes(x = ResponseVariable, y = Value, fill = Metric)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(
    x = "Response Variable",
    y = "Value",
    fill = "Metric"
  ) + 
  main_theme_combined+
  theme(panel.grid = element_blank(),
        plot.background = element_rect(fill = "white"),
        panel.background = element_rect(fill = "white")) + scale_y_continuous(expand = c(0, 0), 
                                                                              limits = c(0, 1.025))

Q_Three <- data.frame(
  ResponseVariable = c("Amo", "Chl", "Dif", "Gly", "Imi", "Met", "Oxy", "Teb"),
  OOB = c(0.18375,
          0.17314,
          0.19435,
          0.18021,
          0.18021,
          0.18728,
          0.0106,
          0.17314),
  AUC = c(0.62145,
          0.52011,
          0.57203,
          0.63535,
          0.54255,
          0.49468,
          0.97329,
          0.6055),
  Accuracy = c(0.4458,
               0.4458,
               0.4096,
               0.4096,
               0.4337,
               0.4337,
               0.9277,
               0.4337)
)

Q_Three_long <- Q_Three %>%
  pivot_longer(cols = c(OOB, AUC, Accuracy), names_to = "Metric", values_to = "Value")
Q_Three_long$ResponseVariable <- factor(Q_Three_long$ResponseVariable, levels = c("Amo", "Chl", "Dif", "Gly", "Imi", "Met", "Oxy", "Teb"))
Q_Three_Fig <- ggplot(Q_Three_long, aes(x = ResponseVariable, y = Value, fill = Metric)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(
    x = "Response Variable",
    y = "Value",
    fill = "Metric"
  ) + 
  main_theme_combined+
  theme(panel.grid = element_blank(),
        plot.background = element_rect(fill = "white"),
        panel.background = element_rect(fill = "white")) + scale_y_continuous(expand = c(0, 0), 
                                                                              limits = c(0, 1.025))

Q_two <- data.frame(
  ResponseVariable = c("Amo", "Chl", "Dif", "Gly", "Imi", "Met", "Oxy", "Teb"),
  OOB = c(0.25279,
              0.24164,
              0.28625,
              0.30112,
              0.27881,
              0.26394,
              0.0223,
              0.27509),
  AUC = c(0.55083,
          0.53968,
          0.56717,
          0.47487,
          0.50542,
          0.51402,
          0.97944,
          0.51722), 
  Accuracy = c(0.732,
               0.7216,
               0.732,
               0.7113,
               0.7526,
               0.732,
               0.8763,
               0.732)
)

Q_Two_long <- Q_two %>%
  pivot_longer(cols = c(OOB, AUC, Accuracy), names_to = "Metric", values_to = "Value")
Q_Two_long$ResponseVariable <- factor(Q_Two_long$ResponseVariable, levels = c("Amo", "Chl", "Dif", "Gly", "Imi", "Met", "Oxy", "Teb"))
Q_two_Fig <- ggplot(Q_Two_long, aes(x = ResponseVariable, y = Value, fill = Metric)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(
    x = "Response Variable",
    y = "Value",
    fill = "Metric"
  ) + 
  main_theme_combined+
  theme(panel.grid = element_blank(),
        plot.background = element_rect(fill = "white"),
        panel.background = element_rect(fill = "white")) + scale_y_continuous(expand = c(0, 0), 
                                                                              limits = c(0, 1.025))

#Q4 single
Q_f_sin <- data.frame(
  ResponseVariable = c("Com", "Amo", "Chl", "Dif", "Gly", "Imi", "Met", "Oxy", "Teb"),
  OOB = c(0.00386,
              0.27027,
              0.25869,
              0.30502,
              0.32432,
              0.30888,
              0.32819,
              0.02317,
              0.28185),
  AUC = c(1,
          0.6201,
          0.7282,
          0.5442,
          0.5492,
          0.6294,
          0.5861,
          0.9812,
          0.6012), 
  Accuracy = c(0.99065,
               0.7944,
               0.7009,
               0.7383,
               0.6822,
               0.729,
               0.7477,
               0.9813,
               0.6822)
)

Q_f_sin_long <- Q_f_sin %>%
  pivot_longer(cols = c(OOB, AUC, Accuracy), names_to = "Metric", values_to = "Value")
Q_f_sin_long$ResponseVariable <- factor(Q_f_sin_long$ResponseVariable, levels = c("Com", "Amo", "Chl", "Dif", "Gly", "Imi", "Met", "Oxy", "Teb"))
Q_f_sin_Fig <- ggplot(Q_f_sin_long, aes(x = ResponseVariable, y = Value, fill = Metric)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(
    x = "Response Variable",
    y = "Value",
    fill = "Metric"
  ) + 
  main_theme_combined+
  theme(panel.grid = element_blank(),
        plot.background = element_rect(fill = "white"),
        panel.background = element_rect(fill = "white")) + scale_y_continuous(expand = c(0, 0), 
                                                                              limits = c(0, 1.025))

#Q four Com and Chem
Q_f_com_chem <- data.frame(
  ResponseVariable = c("Com & Amo", "Com & Chl", "Com & Dif", "Com & Gly", "Com & Imi", 
                       "Com & Met", "Com & Oxy", "Com & Teb"),
  ComOOB = c(0.00386,
                 0.00386,
                 0.00386,
                 0.00386,
                 0.00386,
                 0.00386,
                 0.00386,
                 0.00386),
  ChemOOB = c(0.28185,
                  0.24324,
                  0.30502,
                  0.30888,
                  0.28185,
                  0.32046,
                  0.02703,
                  0.25483),
  ComAUC = c(1,
             1,
             1,
             1,
             1,
             1,
             1,
             1),
  ChemAUC = c(0.69786,
              0.73968,
              0.5584,
              0.54759,
              0.66744,
              0.57678,
              0.98009,
              0.63434),
  ComAccuracy = c(0.99065,
                  1,
                  1,
                  1,
                  0.99065,
                  0.99065,
                  0.97196,
                  1),
  ChemAccuracy = c(0.7664,
                   0.7383,
                   0.7196,
                   0.7196,
                   0.7477,
                   0.7196,
                   0.9813,
                   0.6822)
)


Q_f_com_chem_long <- Q_f_com_chem %>%
  pivot_longer(cols = c(ComOOB, ChemOOB, ComAUC, ChemAUC, ComAccuracy, ChemAccuracy), names_to = "Metric", values_to = "Value")
Q_f_com_chem_long$ResponseVariable <- factor(Q_f_com_chem_long$ResponseVariable, levels = c("Com & Amo", "Com & Chl", "Com & Dif", "Com & Gly", "Com & Imi", 
                                                                                            "Com & Met", "Com & Oxy", "Com & Teb"))

Q_f_com_chem_fig <- ggplot(Q_f_com_chem_long, aes(x = ResponseVariable, y = Value, fill = Metric)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(
    x = "Response Variable",
    y = "Value",
    fill = "Metric"
  ) +
  main_theme_combined+
  theme(panel.grid = element_blank(),
        plot.background = element_rect(fill = "white"),
        panel.background = element_rect(fill = "white")) + scale_y_continuous(expand = c(0, 0), 
                                                                              limits = c(0, 1.025))
  

#Q4 no Oxy
Q_f_no_oxy <- data.frame(
  ResponseVariable = c("Com", "Amo", "Chl", "Dif", "Gly", "Imi", "Met", "Teb"),
  OOB = c(0,
              0.22995,
              0.14439,
              0.24599,
              0.28342,
              0.24599,
              0.21925,
              0.29947),
  AUC = c(1,
          0.86207,
          0.93534,
          0.91157,
          0.87454,
          0.85001,
          0.93054,
          0.93951),
  Accuracy = c(1,
               0.8942,
               0.9524,
               0.9206,
               0.8995,
               0.9048,
               0.8677,
               0.9259)
)
Q_f_no_oxy_long <- Q_f_no_oxy %>%
  pivot_longer(cols = c(OOB, AUC, Accuracy), names_to = "Metric", values_to = "Value")
Q_f_no_oxy_long$ResponseVariable <- factor(Q_f_no_oxy_long$ResponseVariable, levels = c("Com", "Amo", "Chl", "Dif", "Gly", "Imi", "Met", "Teb"))
Q_f_no_oxy_Fig <- ggplot(Q_f_no_oxy_long, aes(x = ResponseVariable, y = Value, fill = Metric)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(
    x = "Response Variable",
    y = "Value",
    fill = "Metric"
  ) + 
  main_theme_combined +
  theme(panel.grid = element_blank(),
  plot.background = element_rect(fill = "white"),
   panel.background = element_rect(fill = "white")) + scale_y_continuous(expand = c(0, 0), 
   limits = c(0, 1.025))
#Q_two_Fig-30% 70%
Q_one_Fig <- Q_one_Fig + ggtitle("A")
Q_two_Fig <- Q_two_Fig + ggtitle("B")
Q_Three_Fig <- Q_Three_Fig + ggtitle("C")
Q_f_sin_Fig <- Q_f_sin_Fig + ggtitle("D")
Q_f_com_chem_fig <- Q_f_com_chem_fig + ggtitle("E")
Q_f_no_oxy_Fig <- Q_f_no_oxy_Fig + ggtitle("F")

combined_plots_Hy <- 
  Q_one_Fig+
  Q_two_Fig+
  Q_Three_Fig+
  Q_f_sin_Fig+
  Q_f_com_chem_fig+
  Q_f_no_oxy_Fig+
  plot_layout(ncol = 2)

print(combined_plots_Hy)

#no Oxy com
Q_f_no_oxycom <- data.frame(
  ResponseVariable = c( "Amo", "Chl", "Dif", "Gly", "Imi", "Met", "Teb"),
  OOB = c(0.19697,
          0.15657,
          0.27778,
          0.29798,
          0.27273,
          0.24242,
          0.27273),
  AUC = c(0.43775,
          0.56356,
          0.5589,
          0.51939,
          0.52005,
          0.53529,
          0.47925),
  Accuracy = c(0.7183,
               0.831,
               0.8028,
               0.7465,
               0.8028,
               0.7042,
               0.7324)
)
Q_f_no_oxycom_long <- Q_f_no_oxycom %>%
  pivot_longer(cols = c(OOB, AUC, Accuracy), names_to = "Metric", values_to = "Value")
Q_f_no_oxycom_long$ResponseVariable <- factor(Q_f_no_oxycom_long$ResponseVariable, levels = c("Amo", "Chl", "Dif", "Gly", "Imi", "Met", "Teb"))
Q_f_no_oxycom_Fig <- ggplot(Q_f_no_oxycom_long, aes(x = ResponseVariable, y = Value, fill = Metric)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(
    x = "Response Variable",
    y = "Value",
    fill = "Metric"
  ) + 
  main_theme_combined +
  theme(panel.grid = element_blank(),
        plot.background = element_rect(fill = "white"),
        panel.background = element_rect(fill = "white")) + scale_y_continuous(expand = c(0, 0), 
                                                                              limits = c(0, 1.025))










##Err
Q_f_Oxy_Com_Com_Err <- ggplot(err_rate_Oxy_Com_four_com_long, aes(x = ntrees, y = err_rate, color = model_id, group = model_id)) +
  geom_line() +
  geom_point()+
  geom_smooth(method = "loess", se = FALSE) +
  labs(
    x = "Number of Trees",
    y = "OOB Error Rate",
    color = "Model"
  ) +
  main_theme

Q_f_Oxy_Com_Oxy_Err <- ggplot(err_rate_Oxy_Com_four_Oxy_long, aes(x = ntrees, y = err_rate, color = model_id, group = model_id)) +
  geom_line() +
  geom_point()+
  geom_smooth(method = "loess", se = FALSE) +
  labs(
    x = "Number of Trees",
    y = "OOB Error Rate",
    color = "Model"
  ) +
  main_theme

Q_f_Oxy_Com_Com_Err <- Q_f_Oxy_Com_Com_Err + ggtitle("A")
Q_f_Oxy_Com_Oxy_Err <- Q_f_Oxy_Com_Oxy_Err + ggtitle("B")

combined_plots_Err_Oxy_Com <- 
  Q_f_Oxy_Com_Com_Err+
  Q_f_Oxy_Com_Oxy_Err+
  plot_layout(ncol = 2)


print(combined_plots_Err_Oxy_Com)
