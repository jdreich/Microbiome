
library(devtools)
library(qiime2R)
library(tidyverse)
library(vegan)
library(ggplot2)
library(reshape2)
library(lubridate)
library(phyloseq)
library(ggpubr)
library(rstatix)

# read in rarefied phyloseq object
pseq_rare <- readRDS("outputs/its1_all_rare.rds")

# reorder factors:
sample_data(pseq_rare)$SampleType <- factor(sample_data(pseq_rare)$SampleType,
                                            levels = c("Rotorod19", "Burkard19",
                                                       "Burkard21", "Control"))

sample_data(pseq_rare)$GenericID <- factor(sample_data(pseq_rare)$GenericID,
                                           levels = c("Lethbridge", "Brooks",
                                                      "Lacombe", "Beaverlodge",
                                                      "Bean", "Canola",
                                                      "Potato", "Wheat",
                                                      "Mock", "WaterBlank",
                                                      "RodBlank", "FieldControl"))

# remove 'ascospore' from controls prior to plotting
pseq_rare_filt <- subset_samples(pseq_rare, !(GenericID == "Ascospore"))

# this will be needed for plotting:
facet_names <- c('Rotorod19' = "Rotorod '19",
                 'Burkard19' = "Burkard '19",
                 'Burkard21' = "Burkard '21",
                 'Control' = "Controls")


# get metadata file
all_meta <- sample_data(pseq_rare) %>% 
  as_tibble()

# 1) add alpha diversity metrics:
# compute alpha diversity metrics
alpha_metrics <- estimate_richness(pseq_rare, measures = c("Observed", "Chao1",
                                                           "Shannon", "Simpson"))
# merge with original data file:
alpha_metrics$SampleID <- rownames(alpha_metrics)
all_meta <- merge(all_meta, alpha_metrics, by = "SampleID") 



# statistical tests to add to plots:
stat.test <- filter(all_meta, !(SampleType == "Control")) %>%
  group_by(SampleType) %>%
  tukey_hsd(Shannon ~ GenericID)
stat.test

# prepare for plotting:
stat.test <- stat.test %>% add_y_position()

# plot all samples, with statistical tests, this takes a few seconds to plot:
plot_richness(pseq_rare_filt, x = "GenericID", measures=c("Shannon")) +
  geom_boxplot(outlier.shape = NA, fill = NA) + 
  # geom_jitter(width = 0.2, height = 0) +
  facet_grid(~SampleType, 
             scales = "free", space = "free",
             labeller = as_labeller(facet_names)) +
  theme_classic() + 
  stat_pvalue_manual(stat.test, label = "p.adj.signif", tip.length = 0.01,
                     hide.ns = T) +
  annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf) +
  annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf) +
  theme(axis.text.x = element_text(angle = 30, vjust = 1, hjust = 1),
        strip.background = element_rect(colour="white", fill="white")) + # remove panel border
  labs(y = "Shannon diversity index",
       x = NULL) 

# save file:
ggsave("outputs/ITS1 Shannon diversity.png", width = 6, height = 3.5)




### --------- ROTOROD ALPHA DIVERSITY ------------
# focus on ROTOROD metadata file - this is the data set looking at
#   different geographic locations:
# edit metadata file to change order of Location:
roto_rare <- readRDS("outputs/roto_its1_rare.rds")
sample_data(roto_rare)$Location <- factor(sample_data(roto_rare)$Location, 
                                          levels = c("Lethbridge", "Brooks",
                                                     "Lacombe", "Beaverlodge"))
# get metadata file
roto_meta <- sample_data(roto_rare) %>% 
  as_tibble()

# 1) add alpha diversity metrics:
# compute alpha diversity metrics
alpha_metrics <- estimate_richness(roto_rare)
# merge with original data file:
alpha_metrics$SampleID <- rownames(alpha_metrics)
roto_meta <- merge(roto_meta, alpha_metrics, by = "SampleID") 
# by = 0 indicates use rownames

# 2) add weather variables:
weather <- read_csv("files/Rotorod weather data 2019.csv")

# filter out duplicates in rotorod data set (some samples were run in duplicate,
#   and these are identified by having no date; for now we'll just remove these):
roto_meta_filt <- filter(roto_meta, !is.na(Date)) %>% 
  droplevels()

# ensure week is a factor
weather$week <- factor(weather$week)

# merge all data:
roto_all <- merge(roto_meta_filt, weather, by = c("Location", "Date", "week", "biweekly",
                                                  "month"))

# write.csv(roto_all, "files/Rotorod ITS1 all metadata.csv", row.names = FALSE)


# summarize by week:
alpha_biweekly <- roto_all %>% group_by(Location, week) %>% 
  dplyr::summarise(meanShannon = mean(Shannon, na.rm = T), 
                   sdShannon = sd(Shannon, na.rm = T),
                   number = n()) %>% 
  ungroup()

# plot alpha metrics - remove week 25, since it has few data points;
# looking at data grouped every 2 weeks
filter(roto_all, !(week == "25")) %>% 
  ggplot(aes(x = biweekly, y = Shannon)) +
  geom_boxplot() + 
  geom_point() +
  facet_grid(Location~.) +
  theme_bw()


#  try again, looking at every week:
alpha_biweekly %>% #filter(., !(biweekly == "NA")) %>% 
  ggplot(aes(x = week, y = meanShannon#, col = Location
  )) +
  geom_point(aes(group = Location)) + #position = position_dodge(width = 0.9)
  geom_line(aes(group = Location)) +
  geom_errorbar(aes(ymin = meanShannon - sdShannon,
                    ymax = meanShannon + sdShannon,
                    width = 0.2)
  ) +
  facet_grid(Location ~.) +
  labs(x = "Week (2019)",
       y = "Shannon diversity index") +
  theme_bw() +
  theme(strip.background = element_blank())

ggsave("outputs/Roto ITS1 Shannon diversity by location week.png", width = 5, height = 4)



# Statistical analysis:

# ok, prior to running any analysis you need to make sure the data frame looks like:
#  - everything is a factor, except the numerical response (in this case, Shannon)
#  - you can ONLY HAVE THE VARIABLES IN THE DATA FRAME THAT YOU ARE USING IN THE ANALYSIS!!!
#     extra variables/columns will just result in errors in the model:

# summarize data by week:
roto_summary <- roto_all %>% 
  group_by(Location, week) %>% 
  summarise(meanShannon = mean(Shannon, na.rm = T), no = n()) %>% 
  ungroup()

# include only data of interest for modeling:
to_test <- select(roto_summary, meanShannon, Location, week) %>% 
  mutate(Location = as.factor(Location)) %>% 
  filter(week %in% c("27","28","29","30","31","32","33"))


# summarize weather data to do linear modeling of Diversity ~ Weather

# summary table of relevant statistics:
weather_summary <- weather %>% 
  group_by(Location) %>% 
  summarise(meanT = mean(MeanT, na.rm = T), sdT = sd(MeanT, na.rm =T),
            meanRH = mean(MeanRH, na.rm = T), sdRH = sd(MeanRH, na.rm = T),
            sumPrecip = sum(Precip, na.rm = T), meanPrecip = mean(Precip, na.rm = T))

alpha_summary <- roto_all %>% 
  group_by(Location) %>% 
  summarise(meanShannon = mean(Shannon, na.rm = T), sdShannon = sd(Shannon, na.rm = T),
            meanChao1 = mean(Chao1, na.rm = T), sdChao1 = sd(Chao1, na.rm = T))


roto_its1_summary <- left_join(alpha_summary, weather_summary)

shannon.aov <- roto_all %>% anova_test(Shannon ~ Location)
roto_all %>% tukey_hsd(Shannon ~ Location)
roto_all %>% tukey_hsd(Chao1 ~ Location)

# test whether environmental variables are different by location:
temp.aov <- weather %>% anova_test(MeanT ~ Location)
weather %>% tukey_hsd(MeanT ~ Location)

weather %>% tukey_hsd(MeanRH ~ Location)

weather %>% tukey_hsd(Precip ~ Location)


# also make a summary by week for the whole data set:
weather_weekly <- weather %>% 
  group_by(Location, week) %>% 
  summarise(meanT = mean(MeanT, na.rm = T), sdT = sd(MeanT, na.rm =T),
            meanRH = mean(MeanRH, na.rm = T), sdRH = sd(MeanRH, na.rm = T),
            sumPrecip = sum(Precip, na.rm = T))

alpha_weekly <- roto_all %>% 
  group_by(Location, week) %>% 
  summarise(meanShannon = mean(Shannon, na.rm = T), sdShannon = sd(Shannon, na.rm = T),
            meanChao1 = mean(Chao1, na.rm = T), sdChao1 = sd(Chao1, na.rm = T))


roto_weekly <- left_join(alpha_weekly, weather_weekly)

# explore relationships:
ggplot(aes(x = Precip, y = Simpson, col = Location), data = roto_all) +
  geom_point(size = 2) +
  geom_smooth(method = "lm", se = F) +
  #geom_line() +
  theme_classic()

# remove outlier for precip before modeling:
roto_all_filt <- filter(roto_all, Precip < 40)

# correlation matrix of variables with Shannon:
library("Hmisc")
library("PerformanceAnalytics")
colnames(roto_all_filt)
to_corr <- select(roto_all_filt, Shannon, Precip, MaxT, MinT, MeanT,
                  MaxRH, MinRH, MeanRH) %>% as.matrix()

cors <- rcorr(to_corr)

chart.Correlation(to_corr, histogram=TRUE, pch=19)

# make correlation table by Location:
to_corr_location <- select(roto_all_filt, Location, Shannon, Precip, MaxT, MinT, MeanT,
                           MaxRH, MinRH, MeanRH)

# code got from:
# https://stackoverflow.com/questions/23967079/correlation-matrix-by-group
# pearson correlation:

result <- lapply(split(to_corr_location, to_corr_location[, 1]), 
                 function(x) {
                   rcorr(as.matrix(x[,2:9]), type="pearson")
                 })

# ok, so strength and significance of weather associations depends on Location:
# Lethbridge - negatively associated with MaxRH and MeanRH
# Brooks - positively associated with Precip
# Lacombe - positively associated with Precip, 
#           negatively associated with MaxT
# Beaverlodge - positively associated with MinRH and MeanRH, 
#               negatively associated wth MaxT

# analyze relationships statistically:
lm.alpha1 <- glm(Shannon ~ MaxT + MinRH + Precip + Location #+ 
                 #MaxT*Location + MinRH*Location + Precip*Location
                 , data = roto_all_filt)
Anova(lm.alpha1, type = "III")
summary(lm.alpha1)


# can try stepwise regression as well:
# from the following website:
#  http://www.utstat.toronto.edu/~brunner/oldclass/appliedf11/handouts/2101f11StepwiseLogisticR.pdf
fullmod <- glm(Shannon ~ Location + biweekly + Precip + MaxT + MinT + MeanT +
                 MaxRH + MinRH + MeanRH, data = roto_all_filt)
summary(fullmod)

nullmod <- glm(Shannon ~ 1, data = roto_all_filt)
summary(nullmod)

backwards <- step(fullmod, trace = 0)
formula(backwards)
summary(backwards)

forwards <- step(nullmod,
                 scope=list(lower=formula(nullmod),upper=formula(fullmod)), 
                 direction="forward",
                 trace = 0)

formula(forwards)
summary(forwards)

bothways <- step(nullmod, 
                 list(lower=formula(nullmod),upper=formula(fullmod)),
                 direction="both",trace=0)
formula(bothways)
summary(bothways)

formula(backwards); formula(forwards); formula(bothways)
backwards$aic; forwards$aic; bothways$aic; 
backwards$deviance; forwards$deviance; bothways$deviance;

backwards_modified <- glm(Shannon ~ Location + Precip + MinRH, #+ MeanRH, 
                          data = roto_all_filt)
summary(backwards_modified)
blah <- summary(backwards)
blah
summary(backwards)
backwards_modified$aic; backwards$aic

