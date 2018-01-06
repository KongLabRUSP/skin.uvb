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

# Load data----
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
dt2$AvgVolume <- dt2$Volume/dt2$Number
dt2

# Averge tumor size vs. number of tumors----
plot(dt2$Volume/dt2$Number ~ dt2$Number,
     col = dt2$Group)

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

# Last timepoint only----
dt24 <- droplevels(subset(dt2,
                          MeasureWeek == 24 & 
                            Week == 25))
dt24$AvgVolume[is.nan(dt24$AvgVolume)] <- 0
dt24

dt24$Tumor <- dt24$AvgVolume > 0
t1 <- table(dt24$Tumor,
            dt24$Group)
# kable(t1)
#   |      | No UVB| UA only| SFN only| UVB| UVB+UA| UVB+SFN|
#   |:-----|------:|-------:|--------:|---:|------:|-------:|
#   |FALSE |      6|       6|        6|   0|      0|       5|
#   |TRUE  |      0|       0|        0|  26|     26|      20|

tmp <- droplevels(subset(dt24,
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
dt24$LogVolume <- log10(dt24$Volume)
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


dt24$LogPls1Volume <- log10(dt24$Volume + 1)
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
         data = tmp)
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
            data = tmp)
m2
summary(m2)

# Average tumor volume (log-transformed)----
dt24$LogAvgVolume <- dt24$LogVolume/dt24$Number
dt24$LogAvgVolume[is.nan(dt24$LogAvgVolume)] <- 0
p6 <- ggplot(data = dt24,
             aes(x = Group,
                 y = LogAvgVolume,
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

# sink()
