# |--------------------------------------------------------|
# | Project: Skin UVB SKH1 moue model treated with UA/SFN. |
# | Script: tumor volume/size by treatment groups and time |
# | Scientist: Yuqing (Anne) Yang                          |
# | Data Analysis: Ran Yin, Renyi Wu, Davit Sargsyan       |
# | Created: 04/01/2018                                    |
# |--------------------------------------------------------|
# Header----
require(data.table)
require(ggplot2)
require(knitr)
# require(VennDiagram)
# require(gridExtra)

# Data----
dt1 <- fread("data/Kong Lab-UVB skin study.csv")
dt1
dt1$Group <- factor(dt1$Group,
                    levels = c("No UVB",
                               "UA only",
                               "SFN only",
                               "UVB",
                               "UVB+UA",
                               "UVB+SFN"))
dt1$Week <- factor(as.numeric(substr(dt1$GroupWeek,
                              nchar(dt1$GroupWeek) - 7,
                              nchar(dt1$GroupWeek) - 6)))
dt1$Tattoo <- factor(dt1$Tattoo)
dt1$WeekOfFirst <- as.numeric(dt1$WeekOfFirst)

dt2 <- melt.data.table(dt1,
                       id.vars = c("Group",
                                   "Week",
                                   "Tattoo",
                                   "WeekOfFirst"),
                       measure.vars = 4:9,
                       variable.name = "MeasureWeek",
                       value.name = "Volume")
dt2$MeasureWeek <- as.numeric(as.character(substr(dt2$MeasureWeek,
                                                  2,
                                                  3)))
dt2

dt3 <- melt.data.table(dt1,
                       id.vars = c("Group",
                                   "Week",
                                   "Tattoo",
                                   "WeekOfFirst"),
                       measure.vars = 10:15,
                       variable.name = "MeasureWeek",
                       value.name = "Number")
dt3$MeasureWeek <- as.numeric(as.character(substr(dt3$MeasureWeek,
                                                  2,
                                                  3)))
dt3

dt2 <- merge(dt2, 
             dt3,
             by = c("Group",
                    "Week",
                    "Tattoo",
                    "WeekOfFirst",
                    "MeasureWeek"))
rm(dt3)

# Log transformations and averages----
dt2$AvgVolume <- dt2$Volume/dt2$Number
dt2$AvgVolume[is.nan(dt2$AvgVolume)] <- 0

dt2$Tumor <- dt2$AvgVolume > 0

dt2$LogVolume <- log10(dt2$Volume)
dt2$LogAvgVolume <- dt2$LogVolume/dt2$Number
dt2$LogAvgVolume[is.nan(dt2$LogAvgVolume)] <- 0

dt2$LogPls1Volume <- log10(dt2$Volume + 1)
dt2$LogAvgPls1Volume <- dt2$LogPls1Volume/dt2$Number
dt2$LogAvgPls1Volume[is.nan(dt2$LogAvgPls1Volume)] <- 0

dt2

# Averge tumor size vs. number of tumors----
tmp <- droplevels(subset(dt2,
                         Group %in% c("UVB",
                                      "UVB+UA",
                                      "UVB+SFN")))
tmp

p0 <- ggplot(data = tmp,
               aes(x = Number,
                   y = Volume,
                   fill = Group)) +
  facet_wrap(~ MeasureWeek,
             nrow = 2) +
  geom_point(size = 2,
             alpha = 0.6,
             shape = 21,
             position = position_dodge(0.5)) + 
  scale_x_continuous("Number of Tumors") + 
  scale_y_continuous("Tumor Total Volume") + 
  scale_fill_discrete("Treatment Group") +
  ggtitle("Skin UVB: Tumor Total Volume vs. Number of Tumors by Week") +
  theme(plot.title = element_text(hjust = 0.5))
p0

tiff(filename = "tmp/skin_uvb_tumor_total_volume_vs_number.tiff",
     height = 6,
     width = 10,
     units = 'in',
     res = 300,
     compression = "lzw+p")
print(p0)
graphics.off()

p0.1 <- ggplot(data = tmp,
             aes(x = Number,
                 y = LogAvgPls1Volume,
                 fill = Group)) +
  facet_wrap(~ MeasureWeek,
             nrow = 2) +
  geom_point(size = 2,
             alpha = 0.6,
             shape = 21,
             position = position_dodge(0.5)) + 
  scale_x_continuous("Number of Tumors") + 
  scale_y_continuous("Average of Log10(Volume+1)") + 
  scale_fill_discrete("Treatment Group") +
  ggtitle("Skin UVB: Average Tumor Size vs. Number of Tumors by Week") +
  theme(plot.title = element_text(hjust = 0.5))
p0.1

tiff(filename = "tmp/skin_uvb_tumor_avgpls1volume_vs_number.tiff",
     height = 6,
     width = 10,
     units = 'in',
     res = 300,
     compression = "lzw+p")
print(p0.1)
graphics.off()

# # For now, remove outliers (DS, 04/01/2018)
# dt2 <- droplevels(subset(dt2,
#                          !(Tattoo %in% c("142", 
#                                          "165"))))
# NOTE: both records were confirmed by Anne on 04/01/2018
#       do not exclude!

# Plot total tumor volume over time by treatment group----
p1 <- ggplot(dt2,
             aes(x = MeasureWeek,
                 y = Volume,
                 fill = Week,
                 group = Tattoo)) +
  facet_wrap(~ Group) +
  geom_line(position = position_dodge(0.5)) +
  geom_point(position = position_dodge(0.5),
             shape = 21,
             size = 2,
             alpha = 0.5) +
  scale_fill_discrete(name = "Week Of Sacrifice",
                      labels = paste("Week",
                                     levels(dt2$Week))) +
  scale_x_continuous("Week Of Measurement",
                     breaks = unique(dt2$MeasureWeek)) +
  scale_y_continuous("Total Volume") +
  ggtitle("Skin UVB: Total Tumor Volume by Treatment Group Over Time") +
  theme(axis.text.x = element_text(angle = 45,
                                   hjust = 1),
        # legend.position = "none",
        plot.title = element_text(hjust = 0.5))
p1

tiff(filename = "tmp/skin_uvb_tumor_volume_over_time.tiff",
     height = 5,
     width = 7,
     units = 'in',
     res = 300,
     compression = "lzw+p")
print(p1)
graphics.off()

p2 <- ggplot(dt2,
             aes(x = MeasureWeek,
                 y = Number,
                 fill = Week,
                 group = Tattoo)) +
  facet_wrap(~ Group) +
  geom_line(position = position_dodge(0.5)) +
  geom_point(position = position_dodge(0.5),
             shape = 21,
             size = 2,
             alpha = 0.5) +
  scale_fill_discrete(name = "Week Of Sacrifice",
                      labels = paste("Week",
                                     levels(dt2$Week))) +
  scale_x_continuous("Week Of Measurement",
                     breaks = unique(dt2$MeasureWeek)) +
  scale_y_continuous("Total Number of Tumors") +
  ggtitle("Skin UVB: Tumor Number by Treatment Group Over Time") +
  theme(axis.text.x = element_text(angle = 45,
                                   hjust = 1),
        # legend.position = "none",
        plot.title = element_text(hjust = 0.5))
p2

tiff(filename = "tmp/skin_uvb_tumor_number_over_time.tiff",
     height = 5,
     width = 7,
     units = 'in',
     res = 300,
     compression = "lzw+p")
print(p2)
graphics.off()

p3 <- ggplot(dt2,
             aes(x = MeasureWeek,
                 y = AvgVolume,
                 fill = Week,
                 group = Tattoo)) +
  facet_wrap(~ Group) +
  geom_line(position = position_dodge(0.5)) +
  geom_point(position = position_dodge(0.5),
             shape = 21,
             size = 2,
             alpha = 0.5) +
  scale_fill_discrete(name = "Week Of Sacrifice",
                      labels = paste("Week",
                                     levels(dt2$Week))) +
  scale_x_continuous("Week Of Measurement",
                     breaks = unique(dt2$MeasureWeek)) +
  scale_y_continuous("Average Volume of Tumors") +
  ggtitle("Skin UVB: Average Tumor Volume by Treatment Group Over Time") +
  theme(axis.text.x = element_text(angle = 45,
                                   hjust = 1),
        # legend.position = "none",
        plot.title = element_text(hjust = 0.5))
p3

tiff(filename = "tmp/skin_uvb_tumor_avg_volume_over_time.tiff",
     height = 5,
     width = 7,
     units = 'in',
     res = 300,
     compression = "lzw+p")
print(p3)
graphics.off()

# Last timepoint (Week24) only----
dt24 <- droplevels(subset(dt2,
                          MeasureWeek == 24 & 
                            Week == 25))
t1 <- table(dt24$Tumor,
            dt24$Group)
kable(t1)
#   |      | No UVB| UA only| SFN only| UVB| UVB+UA| UVB+SFN|
#   |:-----|------:|-------:|--------:|---:|------:|-------:|
#   |FALSE |      6|       6|        6|   0|      0|       5|
#   |TRUE  |      0|       0|        0|  26|     26|      20|

dt24trt <- droplevels(subset(dt24,
                         Group %in% c("UVB",
                                      "UVB+UA",
                                      "UVB+SFN")))

# Was there treatment difference at the last timepoint?----
p4 <- ggplot(data = dt24,
             aes(x = Group,
                 y = Volume,
                 fill = Number)) +
  geom_boxplot() +
  geom_point(size = 3,
             shape = 21,
             position = position_jitterdodge()) + 
  scale_x_discrete("Treatment Group") + 
  scale_y_continuous("Tumor Total Volume") + 
  scale_fill_gradient("Number of Tumors",
                      low = "white",
                      high = "red") +
  ggtitle("Skin UVB: Tumor Total Volume by Treatment at Week 24") +
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 45,
                                   hjust = 1))
p4

tiff(filename = "tmp/skin_uvb_tumor_total_volume_week24.tiff",
     height = 5,
     width = 6,
     units = 'in',
     res = 300,
     compression = "lzw+p")
print(p4)
graphics.off()

# Log transformation of volumes----
p4.1 <- ggplot(data = dt24,
               aes(x = Group,
                   y = LogVolume,
                   fill = Number)) +
  geom_boxplot() +
  geom_point(size = 3,
             shape = 21,
             position = position_jitterdodge()) + 
  scale_x_discrete("Treatment Group") + 
  scale_y_continuous("Log10(Tumor Total Volume)") + 
  scale_fill_gradient("Number of Tumors",
                      low = "white",
                      high = "red") +
  ggtitle("Skin UVB: Log Transformed Tumor Total Volume\nby Treatment at Week 24") +
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 45,
                                   hjust = 1))
p4.1

tiff(filename = "tmp/skin_uvb_tumor_total_log10volume_week24.tiff",
     height = 5,
     width = 6,
     units = 'in',
     res = 300,
     compression = "lzw+p")
print(p4.1)
graphics.off()

p4.2 <- ggplot(data = dt24,
             aes(x = Group,
                 y = LogPls1Volume,
                 fill = Number)) +
  geom_boxplot() +
  geom_point(size = 3,
             shape = 21,
             position = position_jitterdodge()) + 
  scale_x_discrete("Treatment Group") + 
  scale_y_continuous("Log10(Tumor Total Volume + 1)") + 
  scale_fill_gradient("Number of Tumors",
                      low = "white",
                      high = "red") +
  ggtitle("Skin UVB: Log Transformed Tumor Total Volume\nby Treatment at Week 24") +
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 45,
                                   hjust = 1))
p4.2

tiff(filename = "tmp/skin_uvb_tumor_total_log10pls1volume_week24.tiff",
     height = 5,
     width = 6,
     units = 'in',
     res = 300,
     compression = "lzw+p")
print(p4.2)
graphics.off()

m1 <- lm(LogPls1Volume ~ Group,
         data = dt24trt)
m1
summary(m1)

p5 <- ggplot(data = dt24,
             aes(x = Group,
                 y = Number,
                 fill = Volume)) +
  geom_boxplot() +
  geom_point(size = 3,
             shape = 21,
             position = position_jitterdodge()) + 
  scale_x_discrete("Treatment Group") + 
  scale_y_continuous("Tumor Total Number") + 
  scale_fill_gradient("Total Volume",
                      low = "white",
                      high = "red") +
  ggtitle("Skin UVB: Tumor Total Number by Treatment at Week 24") +
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 45,
                                   hjust = 1))
p5

tiff(filename = "tmp/skin_uvb_tumor_total_number_week24.tiff",
     height = 5,
     width = 6,
     units = 'in',
     res = 300,
     compression = "lzw+p")
print(p5)
graphics.off()

m2 <- lm(Number ~ Group,
            data = dt24trt)
m2
summary(m2)

# Average tumor volume (log-transformed)----
p6 <- ggplot(data = dt24,
             aes(x = Group,
                 y = LogAvgPls1Volume,
                 fill = Number)) +
  geom_boxplot() +
  geom_point(size = 3,
             shape = 21,
             position = position_jitterdodge()) + 
  scale_x_discrete("Treatment Group") + 
  scale_y_continuous("Average Log10(Tumor Total Volume + 1)") + 
  scale_fill_gradient("Number of Tumors",
                      low = "white",
                      high = "red") +
  ggtitle("Skin UVB: Tumor Average Log Transformed Volume\nby Treatment at Week 24") +
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 45,
                                   hjust = 1))
p6

tiff(filename = "tmp/skin_uvb_tumor_avg_log10volume_week24.tiff",
     height = 5,
     width = 6,
     units = 'in',
     res = 300,
     compression = "lzw+p")
print(p6)
graphics.off()

m3 <- lm(LogAvgPls1Volume ~ Group,
         data = dt24trt)
m3
summary(m3)

# Week 22 only----
dt22 <- droplevels(subset(dt2,
                          MeasureWeek == 22 & 
                            Week == 25))
t2 <- table(dt22$Tumor,
            dt22$Group)
kable(t2)
# |      | No UVB| UA only| SFN only| UVB| UVB+UA| UVB+SFN|
# |:-----|------:|-------:|--------:|---:|------:|-------:|
# |FALSE |      6|       6|        6|   1|      5|       5|
# |TRUE  |      0|       0|        0|  25|     21|      20|

dt22trt <- droplevels(subset(dt22,
                             Group %in% c("UVB",
                                          "UVB+UA",
                                          "UVB+SFN")))

# Was there treatment difference at the last timepoint?----
# Log transformation of volumes----
p4.2 <- ggplot(data = dt22,
               aes(x = Group,
                   y = LogPls1Volume,
                   fill = Number)) +
  geom_boxplot() +
  geom_point(size = 3,
             shape = 21,
             position = position_jitterdodge()) + 
  scale_x_discrete("Treatment Group") + 
  scale_y_continuous("Log10(Tumor Total Volume + 1)") + 
  scale_fill_gradient("Number of Tumors",
                      low = "white",
                      high = "red") +
  ggtitle("Skin UVB: Log Transformed Tumor Total Volume\nby Treatment at Week 22") +
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 45,
                                   hjust = 1))
p4.2

tiff(filename = "tmp/skin_uvb_tumor_total_log10pls1volume_week22.tiff",
     height = 5,
     width = 6,
     units = 'in',
     res = 300,
     compression = "lzw+p")
print(p4.2)
graphics.off()

m1 <- lm(LogPls1Volume ~ Group,
         data = dt22trt)
m1
summary(m1)

p5 <- ggplot(data = dt22,
             aes(x = Group,
                 y = Number,
                 fill = Volume)) +
  geom_boxplot() +
  geom_point(size = 3,
             shape = 21,
             position = position_jitterdodge()) + 
  scale_x_discrete("Treatment Group") + 
  scale_y_continuous("Tumor Total Number") + 
  scale_fill_gradient("Total Volume",
                      low = "white",
                      high = "red") +
  ggtitle("Skin UVB: Tumor Total Number by Treatment at Week 22") +
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 45,
                                   hjust = 1))
p5

tiff(filename = "tmp/skin_uvb_tumor_total_number_week22.tiff",
     height = 5,
     width = 6,
     units = 'in',
     res = 300,
     compression = "lzw+p")
print(p5)
graphics.off()

m2 <- lm(Number ~ Group,
         data = dt22trt)
m2
summary(m2)

# Average tumor volume (log-transformed)----
p6 <- ggplot(data = dt22,
             aes(x = Group,
                 y = LogAvgPls1Volume,
                 fill = Number)) +
  geom_boxplot() +
  geom_point(size = 3,
             shape = 21,
             position = position_jitterdodge()) + 
  scale_x_discrete("Treatment Group") + 
  scale_y_continuous("Average Log10(Tumor Total Volume + 1)") + 
  scale_fill_gradient("Number of Tumors",
                      low = "white",
                      high = "red") +
  ggtitle("Skin UVB: Tumor Average Log Transformed Volume\nby Treatment at Week 22") +
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 45,
                                   hjust = 1))
p6

tiff(filename = "tmp/skin_uvb_tumor_avg_log10volume_week22.tiff",
     height = 5,
     width = 6,
     units = 'in',
     res = 300,
     compression = "lzw+p")
print(p6)
graphics.off()

m3 <- lm(LogAvgPls1Volume ~ Group,
         data = dt22trt)
m3
summary(m3)

# Week 20 only----
dt20 <- droplevels(subset(dt2,
                          MeasureWeek == 20 & 
                            Week %in% c(20, 25)))
t3 <- table(dt20$Tumor,
            dt20$Group)
kable(t3)
# |      | No UVB| UA only| SFN only| UVB| UVB+UA| UVB+SFN|
# |:-----|------:|-------:|--------:|---:|------:|-------:|
# |FALSE |     10|      10|       10|  13|     15|      24|
# |TRUE  |      0|       0|        0|  19|     17|       7|

dt20trt <- droplevels(subset(dt20,
                             Group %in% c("UVB",
                                          "UVB+UA",
                                          "UVB+SFN")))

# Was there treatment difference at the last timepoint?----
# Log transformation of volumes----
p4.2 <- ggplot(data = dt20,
               aes(x = Group,
                   y = LogPls1Volume,
                   fill = Number)) +
  geom_boxplot() +
  geom_point(size = 3,
             shape = 21,
             position = position_jitterdodge()) + 
  scale_x_discrete("Treatment Group") + 
  scale_y_continuous("Log10(Tumor Total Volume + 1)") + 
  scale_fill_gradient("Number of Tumors",
                      low = "white",
                      high = "red") +
  ggtitle("Skin UVB: Log Transformed Tumor Total Volume\nby Treatment at Week 20") +
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 45,
                                   hjust = 1))
p4.2

tiff(filename = "tmp/skin_uvb_tumor_total_log10pls1volume_week20.tiff",
     height = 5,
     width = 6,
     units = 'in',
     res = 300,
     compression = "lzw+p")
print(p4.2)
graphics.off()

m1 <- lm(LogPls1Volume ~ Group,
         data = dt20trt)
m1
summary(m1)

p5 <- ggplot(data = dt20,
             aes(x = Group,
                 y = Number,
                 fill = Volume)) +
  geom_boxplot() +
  geom_point(size = 3,
             shape = 21,
             position = position_jitterdodge()) + 
  scale_x_discrete("Treatment Group") + 
  scale_y_continuous("Tumor Total Number") + 
  scale_fill_gradient("Total Volume",
                      low = "white",
                      high = "red") +
  ggtitle("Skin UVB: Tumor Total Number by Treatment at Week 20") +
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 45,
                                   hjust = 1))
p5

tiff(filename = "tmp/skin_uvb_tumor_total_number_week20.tiff",
     height = 5,
     width = 6,
     units = 'in',
     res = 300,
     compression = "lzw+p")
print(p5)
graphics.off()

m2 <- lm(Number ~ Group,
         data = dt20trt)
m2
summary(m2)

# Poisson distribution for counts----
m2.1 <- glm(Number ~ Group,
            data = dt20trt,
            family = "poisson")
m2.1
summary(m2.1)

# Average tumor volume (log-transformed)----
p6 <- ggplot(data = dt20,
             aes(x = Group,
                 y = LogAvgPls1Volume,
                 fill = Number)) +
  geom_boxplot() +
  geom_point(size = 3,
             shape = 21,
             position = position_jitterdodge()) + 
  scale_x_discrete("Treatment Group") + 
  scale_y_continuous("Average Log10(Tumor Total Volume + 1)") + 
  scale_fill_gradient("Number of Tumors",
                      low = "white",
                      high = "red") +
  ggtitle("Skin UVB: Tumor Average Log Transformed Volume\nby Treatment at Week 20") +
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 45,
                                   hjust = 1))
p6

tiff(filename = "tmp/skin_uvb_tumor_avg_log10volume_week20.tiff",
     height = 5,
     width = 6,
     units = 'in',
     res = 300,
     compression = "lzw+p")
print(p6)
graphics.off()

m3 <- lm(LogAvgPls1Volume ~ Group,
         data = dt20trt)
m3
summary(m3)

# sink()