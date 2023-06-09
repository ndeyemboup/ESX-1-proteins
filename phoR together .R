boxphoR1 =  data.frame(phoR_E291)
boxphoR2 =  data.frame(phor_r262)
boxphoR4 =  data.frame(phord302)
boxphoR5 =  data.frame(phorl307)
boxphoR8 =  data.frame(phorl310)
boxphoR3 =  data.frame(phor_r276)

boxphoR1$group <- "phoR-E291"
boxphoR2$group <- "phoR-R262"
boxphoR4$group <- "phoR-D302"
boxphoR5$group <- "phoR-L307"
boxphoR8$group <- "phoR-L310"
boxphoR3$group <- "phoR-R276"

df <- rbind(boxphoR1, boxphoR2, boxphoR4, boxphoR5, boxphoR8, boxphoR3)
df$Geneposition <- factor(df$Geneposition, levels = c("phoR-E291", "phoR-R262", "phoR-D302", 
                                                      "phoR-L307", "phoR-L310", "phoR-R276"))

df <- na.omit(df)
boxplots <- ggplot(df, aes(x = Geneposition, y = Energy)) +
  geom_boxplot() + 
  geom_jitter(color = "black", size = 1.1, alpha = 0.9) +
  ylab("ΔΔG") +
  scale_x_discrete(labels = c("phoR-E291", "phoR-R262", "phoR-D302", 
                              "phoR-L307", "phoR-L310", "phoR-R276"))
df_labels <- data.frame(group = c("phoR-E291", "phoR-E291", "phoR-E291", "phoR-R262", "phoR-R262", "phoR-R262", "phoR-D302", 
                                  "phoR-L307", "phoR-L310"),
                        y = c(0.201003, 0.824848, -0.219369, 0.920312, 0.477298, -0.0723743, -0.160726, 4.09687, 1.80268),
                        label = c("A", "G", "K", "C", "H", "L", "N", "F", "Q"))
 

boxplots <- ggplot(df, aes(x = Geneposition, y = Energy)) +
  geom_boxplot() + 
  geom_jitter(color = "black", size = 1.1, alpha = 0.9) +
  ylab("ΔΔG") +
  ggtitle("ΔΔG of mutation in PhoR protein") + 
  scale_x_discrete(labels = c("phoR-E291", "phoR-R262", "phoR-D302", 
                              "phoR-L307", "phoR-L310", "phoR-R276")) +
  geom_point(data = df_labels, aes(x = group, y = y), color = "red", size = 2) +
  geom_text(data = df_labels, aes(x = group, y = y, label = label), 
            color = "red", size = 5, vjust = -0.4, hjust = -0.5, 
            fontface = "bold", family = "sans") +
  annotate("rect", fill = "lightgreen", alpha = 0.5, 
           xmin = -Inf, xmax = Inf,
           ymin = -1, ymax = -0.5) +
  annotate("rect", fill = "red", alpha = 0.5, 
           xmin = -Inf, xmax = Inf,
           ymin = 0.5, ymax = 1) +
  scale_y_continuous(limits = c(-1.5, 6.5))



                   
                   
                   
                   
                   
                   
                   









