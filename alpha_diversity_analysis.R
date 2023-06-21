
library(devtools)
library(qiime2R)
library(tidyverse)
library(vegan)
library(ggplot2)
library(reshape2)
library(lubridate)
library(phyloseq)


# see this link for a great overview:
# http://www.evolution.unibas.ch/walser/bacteria_community_analysis/2015-02-10_MBM_tutorial_combined.pdf
# look at Chao1 and Shannon
# NB: results from not filtering and not rarefying are almost exactly the same as
#     results following that procedure

colors <- c("#fd8c6e", "#78b4c6", "#7fca76", "#ac8bf8")

# Shannon diversity plots of all data:
# sample_data(pseq_rare)$SampleType <- ifelse(sample_data(pseq_rare)$Sample_or_Control == "Control",
#                                             "Control", ifelse(sample_data(pseq_rare)$Sampler_Type == "Rotorod",
#                                                               sample_data(pseq_rare)$Location,
#                                                               sample_data(pseq_rare)$Crop_Type))

pseq_rare <- readRDS("Phyloseq objects/its1_all_rare.rds")

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

facet_names <- c('Rotorod19' = "Rotorod '19",
                 'Burkard19' = "Burkard '19",
                 'Burkard21' = "Burkard '21",
                 'Control' = "Controls")

# edit metadata file to change order of Location:
#sample_data(pseq_rare_filt)$Location <- factor(sample_data(roto)$Location, 
#                                          levels = c("Lethbridge", "Brooks",
#                                                     "Lacombe", "Beaverlodge"))
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


# reorder SampleType data:
#roto_meta$SampleType <- factor(roto_meta$SampleType, 
#                                  levels = c("Rotorod19", "Burkard19",
#                                             "Burkard21", "Control"))

library(rstatix)
# try to add statistics to boxplot; see:
#  https://www.datanovia.com/en/blog/how-to-add-p-values-to-ggplot-facets/
stat.test <- filter(all_meta, !(SampleType == "Control")) %>%
  group_by(SampleType) %>%
  tukey_hsd(Shannon ~ GenericID)
stat.test

#filter(roto_meta, SampleType == "Burkard19") %>% tukey_hsd(Shannon ~ GenericID)

stat.test <- stat.test %>% add_y_position()
#ggboxplot(df, x = "dose", y = "len", fill = "#FC4E07", facet.by = "supp") +
#  stat_pvalue_manual(stat.test, label = "p.adj.signif", tip.length = 0.01) +
#  scale_y_continuous(expand = expansion(mult = c(0.05, 0.1)))

library(ggpubr)
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

ggsave("ITS1 Shannon diversity all data sets.png", width = 6, height = 3.5)

# analysis of Rotorod data set:
plot_richness(roto_rare, x = "Location", measures=c("Shannon")) +
  geom_boxplot(outlier.shape = NA, fill = NA) + 
  #geom_jitter(width = 0.2, height = 0) +
  theme_bw() + 
  facet_wrap(~biweekly) +
  theme(text = element_text(size = 16),
        axis.text.x = element_text(angle = 30, vjust = 1, hjust = 1),
        strip.background = element_rect(colour="white", fill="white"))  + # remove panel border
  labs(y = "Alpha Diversity Index") 



############### ROTOROD ALPHA DIVERSITY ------------
# focus on ROTOROD metadata file:
# edit metadata file to change order of Location:
roto_rare <- readRDS("Phyloseq objects/roto_its1_rare.rds")
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
weather <- read_csv("Rotorod weather data 2019.csv")

# filter out duplicates in rotorod data set:
roto_meta_filt <- filter(roto_meta, !is.na(Date)) %>% 
  droplevels()
weather$week <- factor(weather$week)

roto_all <- merge(roto_meta_filt, weather, by = c("Location", "Date", "week", "biweekly",
                                                  "month"))

write.csv(roto_all, "Rotorod ITS1 all metadata.csv", row.names = FALSE)


# summarize by week:
alpha_biweekly <- roto_all %>% group_by(Location, week) %>% 
  dplyr::summarise(meanShannon = mean(Shannon, na.rm = T), 
                   sdShannon = sd(Shannon, na.rm = T),
                   number = n()) %>% 
  ungroup()

# plot alpha metrics - remove week 25, since it has few data points:
filter(roto_all, !(week == "25")) %>% 
  ggplot(aes(x = biweekly, y = Shannon)) +
  geom_boxplot() + 
  geom_point() +
  facet_grid(Location~.) +
  theme_bw()


# can set dodge width for line chart:
# pd <- position_dodge(width = 0.4)

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

ggsave("Roto ITS1 Shannon diversity by location week.png", width = 5, height = 4)



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

# must use an ungrouped data frame for this!!!!!
# test whether diversity changes over time:
res.aov <- anova_test(data = to_test, 
                      dv = meanShannon, 
                      wid = Location,
                      within = week)

get_anova_table(res.aov)

# test whether diversity varies by location:
# https://www.datanovia.com/en/lessons/anova-in-r/
#ggboxplot(roto_all, x = "Location", y = "Shannon")

model <- lm(Shannon~Location, data = roto_all) 
ggqqplot(residuals(model))

# test normality assumptions by Location:
roto_all %>%
  group_by(Location) %>%
  shapiro_test(Shannon)

# can also check qqplot by Location:
ggqqplot(roto_all, "Shannon", facet.by = "Location")

# perform statistical test for Shannon:
res.aov <- roto_all %>% anova_test(Shannon ~ Location)
res.aov

# pairwise comparisons:
pwc <- roto_all %>% tukey_hsd(Shannon ~ Location)
pwc

# repeat above for Chao1:
chao.aov <- roto_all %>% anova_test(Chao1 ~ Location)
roto_all %>% tukey_hsd(Chao1 ~ Location)

# visualize:
pwc <- pwc %>% add_xy_position(x = "Location")
ggboxplot(roto_all, x = "Location", y = "Shannon") +
  stat_pvalue_manual(pwc, hide.ns = TRUE) +
  labs(
    subtitle = get_test_label(res.aov, detailed = TRUE),
    caption = get_pwc_label(pwc)
  )


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

########### BURKARD ANALYSIS -------------

# merge both Burk data sets:
burk2019_rare <- readRDS("Phyloseq objects/burk2019_its1_rare.rds")
burk2021_rare <- readRDS("Phyloseq objects/burk2021_its1_rare.rds")
burk <- merge_phyloseq(burk2019_rare, burk2021_rare)

burk_meta <- sample_data(burk) %>% data.frame()

alpha_metrics_burk <- estimate_richness(burk, measures = c("Observed", "Chao1", "Shannon",
                                                           "Simpson"))
# merge with original data file:
alpha_metrics_burk$SampleID <- rownames(alpha_metrics_burk)
burk_meta <- merge(burk_meta, alpha_metrics_burk, by = "SampleID") 

burk_biweekly <- burk_meta %>% 
  filter(., !(SampleType == "Control")) %>%  # filter out controls for this analysis
  group_by(Crop_Type, Year, biweekly) %>% 
  summarise(meanShannon = mean(Shannon, na.rm = T), sdShannon = sd(Shannon, na.rm = T),
            n = n())

# add in missing data point to try and avoid plotting lines through missing points:
to_bind <- data.frame(Crop_Type = "Canola", Year = 2019, biweekly = "30-31", 
                      meanShannon = NA, sdShannon = NA,
                      n = 0)

burk_weekly <- rbind(burk_biweekly, to_bind)

burk_biweekly %>% #filter(., !(biweekly == "NA")) %>% 
  ggplot(aes(x = as.factor(biweekly), y = meanShannon#, col = Location
  )) +
  geom_point(aes(group = Crop_Type)) + #position = position_dodge(width = 0.9)
  geom_line(aes(group = Crop_Type), na.rm = T) +
  geom_errorbar(aes(ymin = meanShannon - sdShannon,
                    ymax = meanShannon + sdShannon,
                    width = 0.2)
  ) +
  facet_grid(Crop_Type ~Year) +
  labs(x = "Week",
       y = "Shannon diversity index") +
  theme_bw() +
  theme(strip.background = element_blank())

ggsave("Shannon diversity for Burkard ITS1 both years.png", width = 5, height = 5)


### repeated measures ANOVA for influence of biweek on Shannon diversity:
# summarize data by week:
burk_summary <- burk_meta %>% 
  filter(., !(SampleType == "Control")) %>% 
  group_by(Year, Crop_Type, biweekly) %>% 
  summarise(meanShannon = mean(Shannon, na.rm = T), no = n()) %>% 
  ungroup()

# # or can I just use the raw data? no??
# to_test_2019 <- burk_meta %>% 
#   filter(., !(SampleType == "Control") & Year == "2019") %>%
#   select(., Crop_Type, Shannon, biweekly)

to_test_2019 <- burk_summary %>% 
  filter(., Year == "2019") %>%
  select(., Crop_Type, meanShannon, biweekly)

to_test_2021 <- burk_summary %>% 
  filter(., Year == "2021") %>%
  select(., Crop_Type, meanShannon, biweekly)


# include only data of interest for modeling:
#to_test <- select(roto_summary, meanShannon, Location, week) %>% 
#  mutate(Location = as.factor(Location)) %>% 
#  filter(week %in% c("27","28","29","30","31","32","33"))

# must use an ungrouped data frame for this!!!!!
# test whether diversity changes over time:
res.aov <- anova_test(data = to_test_2019, 
                      dv = meanShannon, 
                      wid = Crop_Type,
                      within = biweekly)

get_anova_table(res.aov)



res.aov.2021 <- anova_test(data = to_test_2021, 
                           dv = meanShannon, 
                           wid = Crop_Type,
                           within = biweekly)

get_anova_table(res.aov.2021)
#







############ REPEAT FOR CROP TYPE ANALYSIS ----------
burk19_rare <- readRDS("Phyloseq objects/burk2019_its1_rare.rds")
sample_data(burk19_rare)$Crop_Type <- factor(sample_data(burk19_rare)$Crop_Type, 
                                             levels = c("Bean", "Canola",
                                                        "Potato", "Wheat"))
# get metadata file
burk19_meta <- sample_data(burk19_rare) %>% 
  as_tibble()

# 1) add alpha diversity metrics:
# compute alpha diversity metrics
alpha_metrics <- estimate_richness(burk19_rare)
# merge with original data file:
alpha_metrics$SampleID <- rownames(alpha_metrics)
burk19_meta <- merge(burk19_meta, alpha_metrics, by = "SampleID") 

# filter out duplicates in rotorod data set:
burk19_meta_filt <- filter(burk19_meta, Control == "No") %>% 
  droplevels()


write.csv(burk19_meta_filt, "Burkard 2019 ITS1 all metadata.csv", row.names = FALSE)


# summarize by week:
alpha_biweekly <- burk19_meta_filt %>% group_by(Crop_Type, week) %>% 
  dplyr::summarise(meanShannon = mean(Shannon, na.rm = T), 
                   sdShannon = sd(Shannon, na.rm = T),
                   number = n()) %>% 
  ungroup()

# plot alpha metrics - remove week 25, since it has few data points:
filter(burk19_meta_filt, !(week == "25")) %>% 
  ggplot(aes(x = biweekly, y = Shannon)) +
  geom_boxplot() + 
  geom_point() +
  facet_grid(Crop_Type~.) +
  theme_bw()


# can set dodge width for line chart:
# pd <- position_dodge(width = 0.4)

alpha_biweekly %>% #filter(., !(biweekly == "NA")) %>% 
  ggplot(aes(x = biweekly, y = meanShannon#, col = Location
  )) +
  geom_point(aes(group = Crop_Type)) + #position = position_dodge(width = 0.9)
  geom_line(aes(group = Crop_Type)) +
  geom_errorbar(aes(ymin = meanShannon - sdShannon,
                    ymax = meanShannon + sdShannon,
                    width = 0.2)
  ) +
  facet_grid(Crop_Type ~.) +
  labs(x = "Week (2019)",
       y = "Shannon diversity index") +
  theme_bw() +
  theme(strip.background = element_blank())

#ggsave("Roto ITS1 Shannon diversity by location week.png", width = 5, height = 4)



# Statistical analysis:

# Statistical analysis:
shannon.aov <- burk_meta_filt %>% anova_test(Shannon ~ Crop_Type)
shannon.aov
burk_meta_filt %>% tukey_hsd(Shannon ~ Crop_Type)





######## BURKARD 2021 ANALYSIS ------------

burk21_rare <- readRDS("Phyloseq objects/burk2021_its1_rare.rds")
sample_data(burk21_rare)$Crop_Type <- factor(sample_data(burk21_rare)$Crop_Type, 
                                             levels = c("Bean", "Canola",
                                                        "Wheat"))
# get metadata file
burk21_meta <- sample_data(burk21_rare) %>% 
  as_tibble()

# 1) add alpha diversity metrics:
# compute alpha diversity metrics
alpha_metrics <- estimate_richness(burk21_rare, measures = c("Shannon", "Chao1", "Observed"))
# merge with original data file:
alpha_metrics$SampleID <- rownames(alpha_metrics)
burk21_meta <- merge(burk21_meta, alpha_metrics, by = "SampleID") 

# filter out duplicates in rotorod data set:
burk21_meta_filt <- filter(burk21_meta, Control == "No") %>% 
  droplevels()


write.csv(burk21_meta_filt, "Burkard 2021 ITS1 all metadata.csv", row.names = FALSE)


# summarize by week:
alpha_biweekly <- burk21_meta_filt %>% group_by(Crop_Type, week) %>% 
  dplyr::summarise(meanShannon = mean(Shannon, na.rm = T), 
                   sdShannon = sd(Shannon, na.rm = T),
                   number = n()) %>% 
  ungroup()

# plot alpha metrics - remove week 25, since it has few data points:
filter(burk21_meta_filt, !(week == "25")) %>% 
  ggplot(aes(x = biweekly, y = Shannon)) +
  geom_boxplot() + 
  geom_point() +
  facet_grid(Crop_Type~.) +
  theme_bw()


# can set dodge width for line chart:
# pd <- position_dodge(width = 0.4)

# alpha_biweekly %>% #filter(., !(biweekly == "NA")) %>% 
#   ggplot(aes(x = biweekly, y = meanShannon#, col = Location
#   )) +
#   geom_point(aes(group = Crop_Type)) + #position = position_dodge(width = 0.9)
#   geom_line(aes(group = Crop_Type)) +
#   geom_errorbar(aes(ymin = meanShannon - sdShannon,
#                     ymax = meanShannon + sdShannon,
#                     width = 0.2)
#   ) +
#   facet_grid(Crop_Type ~.) +
#   labs(x = "Week (2019)",
#        y = "Shannon diversity index") +
#   theme_bw() +
#   theme(strip.background = element_blank())
# 
# #ggsave("Roto ITS1 Shannon diversity by location week.png", width = 5, height = 4)



# Statistical analysis:

# Statistical analysis:
shannon.aov <- burk21_meta_filt %>% anova_test(Shannon ~ Crop_Type)
shannon.aov
burk21_meta_filt %>% tukey_hsd(Shannon ~ Crop_Type)






# look at Simpson and Inverse Simpson
plot_richness(pseq_rare, x = "Location", measures=c("Simpson", "InvSimpson")) +
  geom_boxplot(outlier.shape = NA, fill = NA) + 
  # geom_jitter(width = 0.2, height = 0) +
  theme_classic() + 
  theme(axis.text.x = element_text(angle = 30, vjust = 1, hjust = 1),
        strip.background = element_rect(colour="white", fill="white")) + # remove panel border
  labs(y = "Alpha Diversity Index") 

# filter out controls:
pseq_obj_field <- subset_samples(pseq_rare, !(is.na(Julian_Day)))
pseq_obj_field
# compute alpha diversity metrics
alpha_metrics <- estimate_richness(pseq_obj_field)
# merge with original data file:
alpha_metrics$SampleID <- rownames(alpha_metrics)
metadata <- merge(metadata, alpha_metrics, by = "SampleID")  # by = 0 indicates use rownames
metadata$Location <- factor(metadata$Location,
                            levels = c("Lethbridge", "Brooks",
                                       "Lacombe", "Beaverlodge"))





filtered_meta <- filter(metadata, !(week == "25") | !(week == "34") | 
                          !(week == "26") | !(week == "33")) %>% 
  ungroup()

filtered_summary <- filter(alpha_weekly, !(biweekly == "34-35") &
                             !(biweekly == "NA")) %>% 
  droplevels()

filtered_summary$biweekly <- as.factor(filtered_summary$biweekly) 
filtered_summary <- dplyr::select(filtered_summary, Location, biweekly, meanShannon)

# test for normality:
filtered_summary %>% 
  group_by(biweekly) %>% 
  shapiro_test()

# ok, prior to running any analysis you need to make sure the data frame looks like:
#  - everything is a factor, except the numerical response (in this case, Shannon)
#  - you can ONLY HAVE THE VARIABLES IN THE DATA FRAME THAT YOU ARE USING IN THE ANALYSIS!!!
#     extra variables/columns will just result in errors in the model:
res.aov <- anova_test(data = filtered_summary, 
                      dv = meanShannon, 
                      wid = Location,
                      within = biweekly)

get_anova_table(res.aov)

metadata %>% filter(., !(biweekly == "NA")) %>% 
  ggplot(aes(x = biweekly, y = Shannon, col = Location)) +
  geom_boxplot(position = position_dodge(width = 0.9)) +
  stat_summary(
    fun = median,
    geom = 'line',
    aes(group = Location, colour = Location),
    position = position_dodge(width = 0.9) #this has to be added
  ) +
  labs(x = "Week (2019)",
       y = "Shannon diversity index") +
  #geom_line(aes(group = Location)) +
  facet_grid(Location ~.) +
  theme_bw()




library(rstatix)
# compare all Locations:
lm.all <- lm(Shannon ~ Location, data = metadata)
summary(lm.all)
lm.all.anova <- Anova(lm.all, type = 3)
TukeyHSD(aov(lm.all))
# analysis of Lethbridge:
lm.lethbridge <- aov(Shannon ~ biweekly, data = filter(metadata, Location == "Lethbridge"))
summary(lm.lethbridge)
#blah <- Anova(lm.lethbridge, type = 3)
TukeyHSD(lm.lethbridge)

# analysis of Beaverlodge:
lm.beaverlodge <- lm(Shannon ~ week, data = filter(filtered_meta, Location == "Beaverlodge"))
summary(lm.beaverlodge)
#blah <- Anova(lm.lethbridge, type = 3)
TukeyHSD(lm.beaverlodge)

library(multcomp)
post_test <- glht(lm.beaverlodge,
                  linfct = mcp(week = "Tukey")
)

summary(post_test)

# focus on only the non-control data:
# for statistical analyses:
# build model and check residuals:
model <- lm(Shannon ~ Location, data = metadata)
summary(model)

# residuals look fine

library(rstatix) # library useful for plotting statistics:
# perform ANOVA on Shannon diversity
shannon_anova <- metadata %>% anova_test(Shannon ~ Location)
# need to follow this exactly if wanting to reproduce examples
tukey <- metadata %>% tukey_hsd(Shannon ~ Location)

library(ggpubr)
# and boxplots:
# prepare comparisons:
# see https://www.datanovia.com/en/lessons/anova-in-r/ for examples of this:
tukey <- tukey %>% add_xy_position(x = "Location") 

png(file="Roto ITS1x/roto-its1-alpha2.png",
    width=450, height=400)

ggplot(aes(x = Location, y = Shannon), data = metadata) +
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(#aes(fill = Location), 
    height = 0, width = 0.0, 
    size = 3, shape = 21, fill = "black", alpha = 0.4) +
  # see https://github.com/kassambara/ggpubr/issues/102 for this method:
  stat_pvalue_manual(tukey, hide.ns = TRUE, y.position = c(6.2, 5.8, 5.4)) +
  labs(subtitle = get_test_label(shannon_anova, detailed = TRUE)) +
  ylim(0, NA) +
  theme_classic() +
  theme(text = element_text(size = 16),
        legend.position = "none")

dev.off()


# plot richness as a function of time:
ggplot(aes(x = as.Date(Date.x), y = Shannon, col = Location), data = metadata) +
  geom_line(size = 0.8) +
  geom_point(size = 2) +
  xlab("Date") +
  theme_bw()

# look at scatter plot of number of reads vs Shannon diversity
# (i.e. does more reads result in more diversity - unlikely, since
#  we rarified to the same depth)
ggplot(aes(x = Number.of.Reads.NovaSeq, y = Shannon, group = Location),
       data = metadata) +
  geom_point(aes(col = Location)) +
  theme_classic()
cor(metadata$Number.of.Reads.NovaSeq, metadata$shannon)
cor.test(metadata$Number.of.Reads.NovaSeq, metadata$shannon)
# r = 0.25, p = 0.01223

# investigate environmental variables and their relationship with Shannon
#   diversity:
ggplot(aes(x = MeanT, y = Shannon, group = Location), 
       data = metadata) +
  geom_point(aes(col = Location)) +
  theme_classic()
ggplot(aes(x = MeanRH, y = Shannon, group = Location), 
       data = metadata) +
  geom_point(aes(col = Location)) +
  theme_classic()
ggplot(aes(x = Precip, y = Shannon, group = Location), 
       data = metadata) +
  geom_point(aes(col = Location)) +
  theme_classic()

# can run linear models to assess effects of variables on Shannon diversity:
library(lme4)
shannon_lm1 <- lmer(Shannon ~ Precip + MaxT + MinT + MeanT + 
                      MaxRH + MinRH + MeanRH + (1|Location), 
                    data = metadata)
summary(shannon_lm1)
library(lmerTest)
summary(lmer(Shannon ~ Precip + MaxT + MinT + MeanT + 
               MaxRH + MinRH + MeanRH + (1|Location), 
             data = metadata))