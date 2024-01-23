#https://www.ncbi.nlm.nih.gov/pmc/articles/PMC10094192/
#https://www.ncbi.nlm.nih.gov/pmc/articles/PMC10094192/bin/ehp11896.s002.codeanddata.acco.zip
#complicated Network Plot
library(readxl)
library(tidyverse)
library(ggplot2)

########## ???????????? -----------------
# ???????????????????????????????????????Class(Subpathway)
data1 <- read_excel("Supplemental Excel Table.xlsx", sheet = 6, skip = 1)[,c(2, 8)] %>%
  unique()

# ??????????????????:
data2 <- read_excel("Supplemental Excel Table.xlsx", sheet = 11, skip = 1)[, c(1,2,4)]
data2 <- data2[-nrow(data2), ]
# ??????fill???????????????1???4???:
data2 <- data2 %>% fill(Pollutant) %>% fill(Outcome)

# ????????????:
data <- merge(data1, data2, by.x = "Serum metabolite", by.y = "Metabolite")

# ?????????????????????????????????????????????,???????????????????????????120?????????:
set.seed(123)
data <- data[sample(1:nrow(data), 120), ]

# ??????????????????????????????????????????,????????????????????????:
metabolite_count <- as.data.frame(table(data$`Serum metabolite`))
Polutant_count <- as.data.frame(table(data$Pollutant))
Outcome_count <- as.data.frame(table(data$Outcome))

# ??????????????????????????????????????????:
metabolite_count$x[1:nrow(metabolite_count)] <- 1:nrow(metabolite_count)
metabolite_count$y[1:nrow(metabolite_count)] <- 0
# ?????????????????????????????????????????????:
Polutant_count$x[1:nrow(Polutant_count)] <- c(10, 20, 1, 80)
Polutant_count$y[1:nrow(Polutant_count)] <- 5
# ?????????????????????????????????????????????:
Outcome_count$x[1:nrow(Outcome_count)] <- c(30, 45, 70)
Outcome_count$y[1:nrow(Outcome_count)] <- -5

# ????????????:
data_count <- rbind(rbind(metabolite_count, Polutant_count), Outcome_count)
colnames(data_count)[1] <- "Serum metabolite"
# ?????????????????? -- ??????????????????????????????????????????!
data_count <- left_join(data_count, data[,1:2], by = "Serum metabolite")

data_count$`Super pathway`[c(121:124)] <- "OPEs"
data_count$`Super pathway`[c(125:127)] <- data_count$`Serum metabolite`[c(125:127)]
########## ??????????????? ----------------
colors <- c("#476b71", "#8697a0", "#2da3d1", "#806766",
            "#55bfe2", "#b2bec5", "#d2b698")
names(colors) <- unique(data_count$`Super pathway`)[1:7]
colors <- c(colors, "OPEs" = "#1999a9", "GSP" = "#efc000",
            "HOMA-IR" = "#9a8419", "FPG" = "#1981c8")

p <- ggplot(data_count)+
  geom_point(aes(x, y, size = Freq, color = `Super pathway`))+
  geom_text(data = data_count[1:120,],
            aes(x, y-0.1, label = `Serum metabolite`, color = `Super pathway`),
            angle = 90, hjust = 1, vjust = 0.5, size = 1.5, show.legend = F)+
  geom_text(data = data_count[121:124,],
            aes(x+3, y, label = `Serum metabolite`, color = `Super pathway`),
            angle = 0, hjust = 0, size = 4, show.legend = F)+
  geom_text(data = data_count[125:127,],
            aes(x, y-0.5, label = `Serum metabolite`, color = `Super pathway`),
            angle = 0, hjust = 0.5, vjust = 0.5, size = 4, show.legend = F)+
  scale_color_manual(name = "Class", values = colors)+
  theme_void()

p

########### ??????????????? ---------------
data_line <- data[, c(1,3)]
data_line[121:(120*2),] <- data[, c(1,4)]
data_line$group <- paste0("group", 1:nrow(data_line))
data_line <- pivot_longer(data_line, cols = -group,
                          names_to = "Class", values_to = "Serum metabolite")

data_line <- left_join(data_line, unique(data_count[,c(1,3,4)]), by = "Serum metabolite")

p+geom_line(data = data_line, aes(x, y, group = group, color = "#b7bfcb"),
            linewidth = 0.2, alpha = 0.3, show.legend = F)

ggsave("plot.pdf", height = 5, width = 10)

#############################################################################################
#####Forest plot
library(tidyverse)
library(ggplot2)

# ???????????? ------------
sd <- runif(100, min = 1.5, max = 3)
data <- data.frame(high = runif(100, -2, 8))
data$low = data$high - sd
data$mean <- apply(data, 1, mean)

data <- round(data, 2)

# x?????????:
data$x <- rep(paste0("sample", 1:25), 4)

# ??????????????????:
data$group1 <- rep(c("FPG", "GSP", "FINS", "HOMA-IR"), each = 25)
data$group2 <- rep(rep(c("Blood OPEs", "Urine OPEs"), c(18, 7)), 4)

# ????????????????????????????????? -- ??????????????????:
pos_ind <- sample(which(data$mean > 0), 20)
neg_ind <- sample(which(data$mean < 0), 15)

data$Trend <- "Not Significant"
data$Trend[pos_ind] <- "Positive"
data$Trend[neg_ind] <- "Negative"

data$Signif <- ifelse(data$Trend == "Not Significant", "", "*")

write.csv(data, "data.csv")

p <- ggplot(data)+
  geom_hline(yintercept = 0, linetype = "dashed", color = "#999999")+
  geom_linerange(aes(x, ymin = low, ymax = high, color = Trend))+
  geom_point(aes(x, mean, color = Trend), shape = 21, fill = "white")+
  geom_text(aes(x, high+0.1, label = Signif), color = "black", size = 3)+
  scale_color_manual(values = c("Not Significant" = "#999999",
                                "Positive" = "#b70031",
                                "Negative" = "#0097c4"))+
  scale_y_continuous(breaks = seq(-2, 5, 2))+
  xlab("")+
  ylab("Percentage change(%)")+
  facet_grid(group1 ~ group2, scales = "free", space = "free_x")+
  theme_bw()+
  theme(panel.grid = element_blank(),
        legend.position = "top",
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        strip.text = element_text(color = "white", face = "bold"),
        strip.background.x = element_rect(fill = "#008ea0"),
        strip.background.y = element_blank()
  )

p
ggsave("plot2.pdf", height = 5, width = 10)

#######################box plot#########################################
#file:///C:/Users/hu.p/Downloads/nm.4358.pdf

library(tidyverse)

data <- data.frame(
  Control = c(runif(30, -20, -11)),
  OB = runif(30, -23, -12),
  `0M` = runif(30, -17, -13),
  `1M` = runif(30, -20, -12),
  `3M` = runif(30, -16, -11.5)
)

# ?????????:
data_long <- pivot_longer(data, cols = everything(),
                          names_to = "x", values_to = "value")
# ??????IQR????????????????????????:
data_long <- data_long %>%
  group_by(x) %>%
  mutate(IQR = IQR(value), Med = median(value)) %>%
  ungroup()
library(ggsignif)

# ??????:
ggplot(data_long, aes(x, value, fill = x))+
  # ?????????:
  geom_errorbar(aes(ymin = Med - IQR, ymax = Med + IQR),
                width = 0.5, linetype = "dashed")+
  # ?????????:
  geom_boxplot(color = NA, size = 1,
               notch = TRUE)+
  # ??????????????????:
  geom_rect(aes(xmin = rep(seq(0.8, 4.8, 1), 30),
                xmax = rep(seq(1.2, 5.2, 1), 30),
                ymin = Med-0.2, ymax = Med+0.2), fill = "white")+
  # ????????????:
  scale_fill_manual(values = c("#4087c7", "#ec0f80",
                               "#a24b9c", "#56bc85", "#8bc53f"))+
  # ???????????????:
  geom_signif(aes(x, value),
              # ????????????????????????,???list????????????:
              comparisons = list(c("Control", "OB"),
                                 c("X0M", "X1M"),
                                 c("X0M", "X3M")),
              map_signif_level=T, # ??????P???????????????
              textsize=5, # label??????
              test=t.test, # ????????????
              step_increase=0.1 # ????????????
  )+
  xlab("")+
  ylab("Abundance(lg)")+
  ggtitle("Bacteroides thetaiotaomicron")+
  # ????????????:
  theme_classic()+
  theme(plot.title = element_text(hjust = 0.5, face = "italic"),
        legend.position = "none")

ggsave("plot3.pdf", height = 5, width = 6)
#####################################Figure 4 multibox plot#################
library(tidyverse)
library(ggpubr)
library(ggsignif)
library(rstatix)

# ??????????????????:
data <- data.frame(
  TIME_IA = runif(10, min = 0.05, max = 0.4),
  TIME_ISM = runif(10, min = -0.2, max = 0.1),
  TIME_ISS = runif(10, min = 0.1, max = 0.5),
  TIME_IE = runif(10, min = -0.25, max = 0.2),
  TIME_IR = runif(10, min = 0.2, max = 0.6)
)

# ??????????????????:
data_long <- pivot_longer(data, cols = everything(),
                          names_to = "group", values_to = "Score")

data_long$group <- factor(data_long$group, levels = colnames(data))
# ???????????????:
# ??????t??????:
stat.test <- data_long %>%
  wilcox_test(
    Score ~ group,
    p.adjust.method = "bonferroni"
  )
# ??????:
colors <- c('#eb4b3a', "#48bad0", "#1a9781",
            "#355783", "#ef9a80")
p <- ggplot(data_long)+
  # ?????????:
  geom_boxplot(aes(group, Score, color = group))+
  # ????????????:
  geom_jitter(aes(group, Score, color = group), width = 0.01)+
  # ????????????:
  scale_color_manual(values = c('#eb4b3a', "#48bad0", "#1a9781",
                                "#355783", "#ef9a80"))+
  xlab("")+
  # ??????:
  theme_classic()+
  theme(legend.position = "none",
        # x?????????????????????????????????:
        axis.text.x = element_text(angle = 90, vjust = 0.5, face = "bold",
                                   color = colors))

# ???????????????????????????,?????????????????????:
x_value <- rep(1:4, 4:1)
y_value <- rep(apply(data, 2, max)[1:4], 4:1) + 0.01
y_value <- y_value + c(0.03*1:4, 0.03*1:3, 0.03*1:2, 0.03)
color_value <- c(colors[2:5], colors[3:5], colors[4:5], colors[5])

for (i in 1:nrow(stat.test)) {
  if (stat.test$p.adj.signif[i] != "ns") {
    y_tmp <- y_value[i]
    p <- p+annotate(geom = "text",
                    label = stat.test$p.adj.signif[i],
                    x = x_value[i],
                    y = y_tmp,
                    color = color_value[i])
  }
}
p

ggsave("single_plot.pdf", height = 4, width = 4)