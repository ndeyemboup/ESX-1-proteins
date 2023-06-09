boxphoR1 = data.frame(esxAE31)
boxphoR2 = data.frame(esxAN7)
boxphoR3 = data.frame(esxAW43)
boxphoR4 = data.frame(esxBE14)
boxphoR5 = data.frame(esxBQ28)
boxphoR6 = data.frame(esxBQ42)

boxphoR1$group <- "esxA-E31"
boxphoR2$group <- "esxA-N7"
boxphoR3$group <- "esxA-W43"
boxphoR4$group <- "esxB-E14"
boxphoR5$group <- "esxB-Q28"
boxphoR6$group <- "esxB-Q42"

df <- rbind(boxphoR1, boxphoR2, boxphoR3, boxphoR4, boxphoR5, boxphoR6)
df$Geneposition <- factor(df$Geneposition, levels = c("esxA-E31", "esxA-N7", "esxA-W43", "esxB-E14", 
                                                      "esxB-Q28", "esxB-Q42"))

df <- na.omit(df)
boxplots <- ggplot(df, aes(x = Geneposition, y = Energy)) +
  geom_boxplot() + 
  geom_jitter(color = "black", size = 1.1, alpha = 0.9) +
  ylab("ΔΔG") +
  scale_x_discrete(labels = c("esxA-E31", "esxA-N7", "esxA-W43", "esxB-E14", 
                              "esxB-Q28", "esxB-Q42"))
df_labels <- data.frame(group = c("esxA-E31", "esxA-N7", "esxA-W43", "esxA-W43", "esxB-E14", "esxB-E14",
                                  "esxB-Q28", "esxB-Q42"), 
                        y = c(-0.120124, 0.560816, 1.33742, 0.439987, 0.501124, 1.35151, 2.85799, 0.00859455), 
                        label = c("A", "S", "C", "L", "A", "G", "P", "R"))

boxplots <- ggplot(df, aes(x = Geneposition, y = Energy)) +
  geom_boxplot() + 
  geom_jitter(color = "black", size = 1.1, alpha = 0.9) +
  ylab("ΔΔG") +
  scale_x_discrete(labels = c("esxA-E31", "esxA-N7", "esxA-W43", "esxB-E14", 
                              "esxB-Q28", "esxB-Q42")) +
  geom_point(data = df_labels, aes(x = group, y = y), color = "red", size = 2) +
  geom_text(data = df_labels, aes(x = group, y = y, label = label), 
            color = "red", size = 5, vjust = -0.4, hjust = -0.5, 
            fontface = "bold", family = "sans") 

boxplots + ggtitle("ΔΔG of mutation in esxA-esxB protein")+
  annotate("rect", fill = "lightgreen", alpha = 0.5, 
           xmin = -Inf, xmax = Inf,
           ymin = -1, ymax = -0.5) +
  annotate("rect", fill = "red", alpha = 0.5, 
           xmin = -Inf, xmax = Inf,
           ymin = 0.5, ymax = 1)
  











