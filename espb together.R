box <- ggplot(boxphoR4, aes(x=Geneposition, y=Energy)) +
  geom_boxplot() + 
  geom_jitter(color="black", size=1.1, alpha=0.9) + 
  ylab("ΔΔG")

boxphoR4 =  data.frame(espBK65)
boxphoR4$Geneposition <- factor(boxphoR4$Geneposition, levels = c("espB-K65"))

df3 <- data.frame(x1 = factor("espB-K65", levels = c("espB-K65")), 
                  y1 = 1.21019)

q <- box + 
  geom_point(data=df3, aes(x=x1, y=y1), color="red", size=2) 

label_df <- data.frame(label = c("R"),
                       x1 = c("espB-K65"),
                       y1 = c(1.21019),
                       stringsAsFactors = FALSE)
boxK65 <- q + geom_text(data = label_df, aes(x = x1, y = y1, label = label), color = "red", size = 4, vjust=-0.5, hjust=0)
________
boxQ193 <- ggplot(boxphoR1, aes(x=Geneposition, y=Energy)) +
  geom_boxplot() + 
  geom_jitter(color="black", size=1.1, alpha=0.9) + 
  ylab("ΔΔG")

boxphoR1 =  data.frame(espBQ193)
boxphoR1$Geneposition <- factor(boxphoR1$Geneposition, levels = c("espB-Q193"))
_________
box2 <- ggplot(boxphoR2, aes(x=Geneposition, y=Energy)) +
  geom_boxplot() + 
  geom_jitter(color="black", size=1.1, alpha=0.9) + 
  ylab("ΔΔG")

boxphoR2 =  data.frame(espBK168)
boxphoR2$Geneposition <- factor(boxphoR2$Geneposition, levels = c("espB-K168"))

df1 <- data.frame(x1 = factor("espB-K168", levels = c("espB-K168")), 
                  y1 = 0.878385)

n <- box2 + 
  geom_point(data=df1, aes(x=x1, y=y1), color="red", size=2) 

label_df2 <- data.frame(label = c("E"),
                       x1 = c("espB-K168"),
                       y1 = c(0.878385),
                       stringsAsFactors = FALSE)
boxK168 <- n + geom_text(data = label_df2, aes(x = x1, y = y1, label = label), color = "red", size = 4, vjust=-0.5, hjust=-0.1)
________

box3 <- ggplot(boxphoR3, aes(x=Geneposition, y=Energy)) +
  geom_boxplot() + 
  geom_jitter(color="black", size=1.1, alpha=0.9) + 
  ylab("ΔΔG")

boxphoR3 =  data.frame(espBQ214)
boxphoR3$Geneposition <- factor(boxphoR3$Geneposition, levels = c("espB-Q214"))

df6 <- data.frame(x1 = factor("espB-Q214", levels = c("espB-Q214")), 
                  y1 = 0.230954)

s <- box3 + 
  geom_point(data=df6, aes(x=x1, y=y1), color="red", size=2) 

label_df3 <- data.frame(label = c("E"),
                       x1 = c("espB-Q214"),
                       y1 = c(0.230954),
                       stringsAsFactors = FALSE)
boxQ214 <- s + geom_text(data = label_df3, aes(x = x1, y = y1, label = label), color = "red", size = 4, vjust=-0.5, hjust=0)
_______________________

box6 <- ggplot(boxphoR5, aes(x=Geneposition, y=Energy)) +
  geom_boxplot() + 
  geom_jitter(color="black", size=1.1, alpha=0.9) + 
  ylab("ΔΔG")

boxphoR5 =  data.frame(espBE177)
boxphoR5$Geneposition <- factor(boxphoR5$Geneposition, levels = c("espB-E177"))

df7 <- data.frame(x1 = factor("espB-E177", levels = c("espB-E177")), 
                  y1 = 0.235048)

r <- box6 + 
  geom_point(data=df7, aes(x=x1, y=y1), color="red", size=2) 

label_df6 <- data.frame(label = c("A"),
                       x1 = c("espB-E177"),
                       y1 = c(0.235048),
                       stringsAsFactors = FALSE)
boxE177 <- r + geom_text(data = label_df6, aes(x = x1, y = y1, label = label), color = "red", size = 4, vjust=-0.5, hjust=0)
_________________
boxE258 <- ggplot(boxphoR7, aes(x=Geneposition, y=Energy)) +
  geom_boxplot() + 
  geom_jitter(color="black", size=1.1, alpha=0.9) + 
  ylab("ΔΔG")

boxphoR7 =  data.frame(espBE258)
boxphoR7$Geneposition <- factor(boxphoR7$Geneposition, levels = c("espB-E258"))
________________

boxphoR8 <- data.frame(espBF130)
boxphoR8$Geneposition <- factor(boxphoR8$Geneposition, levels = c("espB-F130"))

box8 <- ggplot(boxphoR8, aes(x=Geneposition, y=Energy)) +
  geom_boxplot() + 
  geom_jitter(color="black", size=1.1, alpha=0.9) + 
  ylab("ΔΔG")

df8 <- data.frame(x1 = factor("espB-F130", levels = c("espB-F130")), 
                  y1 = 0.134554)
df9 <- data.frame(x2 = factor("espB-F130", levels = c("espB-F130")), 
                  y2 = -0.104094)
df10 <- data.frame(x3 = factor("espB-F130", levels = c("espB-F130")), 
                   y3 = 0.365021)

t <- box8 + 
  geom_point(data = df8, aes(x = x1, y = y1), color = "red", size = 2) +
  geom_point(data = df9, aes(x = x2, y = y2), color = "red", size = 2) +
  geom_point(data = df10, aes(x = x3, y = y3), color = "red", size = 2)


label_df8 <- data.frame(label = c("C", "L", "S"),
                       x1 = c("espB-F130", "espB-F130", "espB-F130"),
                       y1 = c(0.134554, -0.104094, 0.365021),
                       stringsAsFactors = FALSE)
boxF130 <- t + geom_text(data = label_df8, aes(x = x1, y = y1, label = label), color = "red", size = 4, vjust=-0.5, hjust=0)
________

box9 <- ggplot(boxphoR9, aes(x=Geneposition, y=Energy)) +
  geom_boxplot() + 
  geom_jitter(color="black", size=1.1, alpha=0.9) + 
  ylab("ΔΔG")

boxphoR9 =  data.frame(espbn129)
boxphoR9$Geneposition <- factor(boxphoR9$Geneposition, levels = c("EspB-N129"))

df11 <- data.frame(x1 = factor("EspB-N129", levels = c("EspB-N129")), 
                  y1 = 0.6927)

m <- box9 + 
  geom_point(data=df11, aes(x=x1, y=y1), color="red", size=2) 

label_df9 <- data.frame(label = c("D"),
                       x1 = c("EspB-N129"),
                       y1 = c(0.6927),
                       stringsAsFactors = FALSE)
boxN129 <- m + geom_text(data = label_df9, aes(x = x1, y = y1, label = label), color = "red", size = 4, vjust=-0.5, hjust=0)
__________________________________________
__________________________________________
boxphoR4$group <- "espB-K65"
boxphoR2$group <- "espB-K168"
boxphoR3$group <- "espB-Q214"
boxphoR5$group <- "espB-E177"
boxphoR1$group <- "espB-Q193"

df <- rbind(boxphoR4, boxphoR2, boxphoR3, boxphoR5, boxphoR1)
df$Geneposition <- factor(df$Geneposition, levels = c("espB-K65", "espB-K168", "espB-Q214",
                                                      "espB-E177", "espB-Q193"))
df <- na.omit(df)
boxplots <- ggplot(df, aes(x = Geneposition, y = Energy)) +
  geom_boxplot() + 
  geom_jitter(color = "black", size = 1.1, alpha = 0.9) +
  ylab("ΔΔG") +
  scale_x_discrete(labels = c("espB-K65", "espB-K168", "espB-Q214", "espB-E177", "espB-Q193"))
df_labels <- data.frame(group = c("espB-K65", "espB-K168", "espB-Q214", "espB-E177"),
                        y = c(1.21019, 0.878385, 0.230954, 0.235048),
                        label = c("R", "E", "E", "A"))

boxplots <- ggplot(df, aes(x = Geneposition, y = Energy)) +
  geom_boxplot() + 
  geom_jitter(color = "black", size = 1.1, alpha = 0.9) +
  ylab("ΔΔG") +
  scale_x_discrete(labels = c("espB-K65", "espB-K168", "espB-Q214", "espB-E177", "espB-Q193")) +
  geom_point(data = df_labels, aes(x = group, y = y), color = "red", size = 2) +
  geom_text(data = df_labels, aes(x = group, y = y, label = label), 
            color = "red", size = 5, vjust = -0.4, hjust = -0.5, 
            fontface = "bold", family = "sans") 

boxplots + ggtitle("ΔΔG of mutation in EspB protein")+
  annotate("rect", fill = "lightgreen", alpha = 0.5, 
           xmin = -Inf, xmax = Inf,
           ymin = -1, ymax = -0.5) +
  annotate("rect", fill = "red", alpha = 0.5, 
           xmin = -Inf, xmax = Inf,
           ymin = 0.5, ymax = 1)+
  scale_y_continuous(limits = c(-2.5, 5))







 

