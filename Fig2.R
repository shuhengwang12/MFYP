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


#Q1
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
          0.448717948717949,
          0.53957,
          0.356884057971014,
          0.58302,
          0.98139,
          0.373514431239389),
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

#Q2 no com
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
  AUC = c(0.449166666666667,
          0.46031746031746,
          0.56717,
          0.525132275132275,
          0.505422374429224,
          0.485978835978836,
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

#Q3
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
          0.479891956782713,
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
#Q4 no oxy
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
          0.937860082304527,
          0.93534,
          0.899748817966903,
          0.912191358024691,
          0.943360471645143,
          0.955393586005831,
          0.841842397336293),
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


#No oxy & com 
Q_f_no_oxycom <- data.frame(
  ResponseVariable = c("Amo", "Chl", "Dif", "Gly", "Imi", "Met", "Teb"),
  OOB = c(0.19697,
          0.15657,
          0.27778,
          0.29798,
          0.27273,
          0.24242,
          0.27273),
  AUC = c(0.562254901960784,
          0.436440677966102,
          0.441102756892231,
          0.480607966457023,
          0.479949874686717,
          0.53529,
          0.520748987854251),
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
#Each com
#Com5
Q_f_eachcom_five <- data.frame(
  ResponseVariable = c("Amo", "Chl", "Dif", "Gly", "Imi", "Met", "Oxy", "Teb"),
  OOB = c(0.2,	0.1, 0.1,	0.2,	0.2,	0.2,	0.1,	0.2),
  AUC = c(0.5,	0,	1, 1,	0.5,	0.5,	1,	0),
  Accuracy = c(0.6667,	0,	0.6667,	0.3333,	0.6667,	0.6667,	1,	0.6667)
)
Q_f_eachcom_five_long <- Q_f_eachcom_five %>%
  pivot_longer(cols = c(OOB, AUC, Accuracy), names_to = "Metric", values_to = "Value")
Q_f_eachcom_five_long$ResponseVariable <- factor(Q_f_eachcom_five_long$ResponseVariable, levels = c("Amo", "Chl", "Dif", "Gly", "Imi", "Met", "Oxy", "Teb"))

Q_f_eachcom_five_Fig <- ggplot(Q_f_eachcom_five_long, aes(x = ResponseVariable, y = Value, fill = Metric)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(
    x = "Response Variable",
    y = "Value",
    fill = "Metric"
  ) + 
  main_theme_combined +
  theme(panel.grid = element_blank(),
        plot.background = element_rect(fill = "white"),
        panel.background = element_rect(fill = "white")) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 1.025))
#Com 8
Q_f_eachcom_eight <- data.frame(
  ResponseVariable = c("Amo", "Chl", "Dif", "Gly", "Imi", "Met", "Oxy", "Teb"),
  OOB = c(0.17857,	0.17857,	0.25,	0.39286,	0.21429,	0.21429,	0.03571,	0.17857),
  AUC = c(0.75,	0.1875,	0.21875,	0.9,	0.40625,	0.5,	1,	0.4375),
  Accuracy = c(0.6667,	0.6667,	0.6667,	0.9167,	0.6667,	0.6667,	1,	0.6667)
)
Q_f_eachcom_eight_long <- Q_f_eachcom_eight %>%
  pivot_longer(cols = c(OOB, AUC, Accuracy), names_to = "Metric", values_to = "Value")
Q_f_eachcom_eight_long$ResponseVariable <- factor(Q_f_eachcom_eight_long$ResponseVariable, levels = c("Amo", "Chl", "Dif", "Gly", "Imi", "Met", "Oxy", "Teb"))

Q_f_eachcom_eight_Fig <- ggplot(Q_f_eachcom_eight_long, aes(x = ResponseVariable, y = Value, fill = Metric)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(
    x = "Response Variable",
    y = "Value",
    fill = "Metric"
  ) + 
  main_theme_combined +
  theme(panel.grid = element_blank(),
        plot.background = element_rect(fill = "white"),
        panel.background = element_rect(fill = "white")) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 1.025))




#Com 12
Q_f_eachcom_twelve <- data.frame(
  ResponseVariable = c("Amo", "Chl", "Dif", "Gly", "Imi", "Met", "Oxy", "Teb"),
  OOB = c(0.1875,	0.25,	0.1875,	0.28125,	0.25,	0.1875,	0,	0.25),
  AUC = c(0.805555555555556,	1,	0.722222222222222,	0.875,	0.5625, 0.638888888888889,	1,	0.8125),
  Accuracy = c(0.5,	0.6667,	0.5,	0.6667,	0.75,	0.5,	1,	0.6667)
)
Q_f_eachcom_twelve_long <- Q_f_eachcom_twelve %>%
  pivot_longer(cols = c(OOB, AUC, Accuracy), names_to = "Metric", values_to = "Value")
Q_f_eachcom_twelve_long$ResponseVariable <- factor(Q_f_eachcom_twelve_long$ResponseVariable, levels = c("Amo", "Chl", "Dif", "Gly", "Imi", "Met", "Oxy", "Teb"))

Q_f_eachcom_twelve_Fig <- ggplot(Q_f_eachcom_twelve_long, aes(x = ResponseVariable, y = Value, fill = Metric)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(
    x = "Response Variable",
    y = "Value",
    fill = "Metric"
  ) + 
  main_theme_combined +
  theme(panel.grid = element_blank(),
        plot.background = element_rect(fill = "white"),
        panel.background = element_rect(fill = "white")) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 1.025))


Q_one_Fig <- Q_one_Fig + ggtitle("A")
Q_two_Fig <- Q_two_Fig + ggtitle("B")
Q_Three_Fig <- Q_Three_Fig + ggtitle("C")
Q_f_no_oxy_Fig <- Q_f_no_oxy_Fig + ggtitle("D")
Q_f_no_oxycom_Fig <- Q_f_no_oxycom_Fig + ggtitle("E")
Q_f_eachcom_five_Fig <- Q_f_eachcom_five_Fig + ggtitle("F")
Q_f_eachcom_eight_Fig <- Q_f_eachcom_eight_Fig + ggtitle("G")
Q_f_eachcom_twelve_Fig <- Q_f_eachcom_twelve_Fig + ggtitle("H")

combined_plots_Hy_new <- 
  Q_one_Fig+
  Q_two_Fig+
  Q_Three_Fig+
  Q_f_no_oxy_Fig+
  Q_f_no_oxycom_Fig+
  Q_f_eachcom_five_Fig+
  Q_f_eachcom_eight_Fig+
  Q_f_eachcom_twelve_Fig+
  plot_layout(ncol = 2)

print(combined_plots_Hy_new)




y_true <- as.factor(c(0, 0, 1, 1))
y_scores <- c(0.9, 0.8, 0.2, 0.3)


roc <- roc(y_true, y_scores, direction = ">")
auc(roc)

