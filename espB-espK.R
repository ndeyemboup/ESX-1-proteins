boxphoR1 = data.frame(espBE258)
boxphoR2 = data.frame(espBF130)
boxphoR3 =  data.frame(espKS667)
boxphoR4 =  data.frame(espKT670)
boxphoR6 =  data.frame(espKD506)

boxphoR1$group <- "espB-E258"
boxphoR2$group <- "espB-F130"
boxphoR3$group <- "espK-S667"
boxphoR4$group <- "espK-T670"
boxphoR6$group <- "espK-D506"
 
df <- rbind(boxphoR1, boxphoR2, boxphoR3, boxphoR4, boxphoR6)
df$Geneposition <- factor(df$Geneposition, levels = c("espB-E258", "espB-F130", "espB-S667", "espK-T670", "espK-D506"))

df <- na.omit(df)
boxplots <- ggplot(df, aes(x = Geneposition, y = Energy)) +
  geom_boxplot() + 
  geom_jitter(color = "black", size = 1.1, alpha = 0.9) +
  ylab("ΔΔG") +
  scale_x_discrete(labels = c("espB-E258", "espB-F130", "espK-S667",
                              "espK-T670", "espK-D506"))
df_labels <- data.frame(group = c("espB-F130", "espB-F130", "espB-F130", "espK-T670", 
                                  "espK-T670", "espK-T670", "espK-D506"), 
                        y = c(0.134554, -0.104094, 0.365021, 0.0751586, 0.155899, 0.20443, 0.203022), 
                        label = c("C", "L", "S", "A", "I","S", "N"))

boxplots <- ggplot(df, aes(x = Geneposition, y = Energy)) +
  geom_boxplot() + 
  geom_jitter(color = "black", size = 1.1, alpha = 0.9) +
  ylab("ΔΔG") +
  ggtitle("ΔΔG of mutation in EspB-EspK protein") + 
  scale_x_discrete(labels = c("espB-E258", "espB-F130","espK-S667", "espK-T670", "espK-D506")) +
  geom_point(data = df_labels, aes(x = group, y = y), color = "red", size = 2) +
  geom_text(data = df_labels, aes(x = group, y = y, label = label, group = group),
            color = "red", size = 4, vjust = -1, position = position_nudge(x = 0.2)) +
  annotate("rect", fill = "lightgreen", alpha = 0.5, 
           xmin = -Inf, xmax = Inf,
           ymin = -1, ymax = -0.5) +
  annotate("rect", fill = "red", alpha = 0.5, 
           xmin = -Inf, xmax = Inf,
           ymin = 0.5, ymax = 1)
 
                        

















                              
                              
                              
                              
                              
                              
                              
                              