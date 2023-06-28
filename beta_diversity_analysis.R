
library(devtools)
library(qiime2R)
library(tidyverse)
library(vegan)
library(ggplot2)
library(reshape2)
library(lubridate)
library(phyloseq)

# set colour palette for plotting:
colors <- c("#fd8c6e", "#78b4c6", "#7fca76", "#ac8bf8")

# read in data set:
all_reads <- readRDS("outputs/its1_all_ra.rds")
metadata <- sample_data(all_reads) %>% data.frame()

# import rotorod metadata for weather analysis:
roto_meta <- read_csv("files/Rotorod ITS1 all metadata.csv") %>% data.frame()
# for merging with a phyloseq object, the metadata has to have rownames corresponding
#  to sample ID:
rownames(roto_meta) <- roto_meta$SampleID

roto <- subset_samples(all_reads, Control == "No" & Sampler_Type == "Rotorod" &
                         !(is.na(Date)))

# update rotorod phyloseq object with relevant metadata:
sample_data(roto) <- roto_meta

sample_data(roto_Hellinger)$Jday <- as.factor(sample_data(roto_Hellinger)$Jday)
sample_data(roto_Hellinger)$Location <- factor(sample_data(roto_Hellinger)$Location,
                                               levels = c("Beaverlodge", "Lacombe",
                                                          "Brooks", "Lethbridge"))

# Hellinger-transform read counts for subsequent analysis (this is one solution
#   to the issue of compositional data):
roto_Hellinger <- transform_sample_counts(roto,
                                          function(x) sqrt(x / sum(x)))


# Can use plot_ordination() for all other metrics as well:
# calculate Bray-Curtis dissimilarity:
ord_bray <- ordinate(roto_Hellinger, method = "PCoA", distance = "bray")

# plot by Location to see how similar samples on same day are:
p1 <- plot_ordination(roto_Hellinger, ord_bray, color = "Location") + 
  #geom_polygon(aes(fill = Jday), alpha = 0.6, show.legend = FALSE) +
  geom_point(aes(fill = Location), shape = 21, size = 5, alpha = 0.7) +
  stat_ellipse(type = "norm") +
  #ggtitle("Bray-Curtis: all ASVs") +
  scale_y_reverse() +
  labs(fill = "Location") +
  scale_color_manual(values = colors) +
  scale_fill_manual(values = colors) +
  #guides(color = FALSE) +
  guides(#fill = guide_legend(override.aes = list(linetype = 0)),
    color = guide_legend(override.aes = list(linetype = 0))) +
  theme_classic() +
  theme(#legend.position = "none", 
    text = element_text(size = 14))

p1$layers <- p1$layers[-1]
p1 + geom_point(col = "darkslategrey", size = 5, shape = 21)


ggsave("outputs/Ordination of Rotorod samplers ITS1.png", height = 3.5, width = 5)

# Run PermANOVA, using vegan's adonis, and including environmental covariates:
#   must first include only complete cases:
sample_df <- data.frame(sample_data(roto_Hellinger))

bray_dist <- phyloseq::distance(roto_Hellinger, method = "bray")
permanova <- adonis(bray_dist ~ Location + Precip + MeanT + MeanRH + week*Location, 
                    by = "margin",
                    data = sample_df, permutations = 10000)$aov.tab
permanova

perm_table <- permanova

# test for dispersion:
mod <- betadisper(bray_dist, sample_df$Location)
dispersion <- anova(mod)
dispersion
plot(mod, hull = F, ellipse = T)
TukeyHSD(mod)

# perform a series of pairwise comparisons for differences in 
#   communities - refer to Pat Schloss's Youtube videos on this:
location <- c("Lethbridge", "Brooks", "Lacombe", "Beaverlodge")

# distance matrix of all samples for subsetting:
dist_mat <- data.frame(as.matrix(bray_dist))
perma_p <- numeric()
beta_p <- numeric()
set.seed(2468)
for (i in 1:length(location)) {
  
  for (j in 1:length(location)) {
    if (i > j) {
      next
    }
    if (i == j) {
      next
    }
    
    comparison_phylo <- subset_samples(roto_Hellinger, 
                                       Location == location[i] | Location == location[j])
    comparison_meta <- data.frame(sample_data(comparison_phylo))
    
    comparison_dist <- phyloseq::distance(comparison_phylo, method = "bray")
    #comparison_dist <- dist_mat %>% select(., comparison_meta$SampleID) %>% 
    #  filter(., rownames(.) %in% comparison_meta$SampleID) %>% 
    #  as.dist()
    
    permanova <- adonis(comparison_dist ~ Location, data = comparison_meta, 
                        permutations = 10000)$aov.tab
    mod <- anova(betadisper(comparison_dist, comparison_meta$Location))
    
    perma_p <- append(perma_p, permanova$`Pr(>F)`[1])
    beta_p <- append(beta_p, mod$`Pr(>F)`[1])
  }
}

# get names for comparisons between Locations:
names <- character()
for (i in 1:length(location)) {
  for (j in 1:length(location)) {
    if (i > j) {
      next
    }
    if (i == j) {
      next
    }
    names <- append(names, paste0(location[i], "_", location[j]))
  }
}

# add names to vector and adjust the p-value:
# adjust p-values:
names(perma_p) <- names
results <- perma_p
results_adj <- p.adjust(perma_p, method = "BH")

# view results:
data.frame(results_adj)


# trying something from Pat Schloss to visualize data:
# https://www.youtube.com/watch?v=roCRaWI3fmw&ab_channel=RiffomonasProject
long_dist <- bray_dist %>% 
  as.matrix() %>% 
  as_tibble(rownames = "samples") %>% 
  pivot_longer(-samples) %>% 
  filter(samples < name)

long_dist <- left_join(long_dist, 
                       select(roto_meta, SampleID, Date, Location, week),
                       by = c("samples" = "SampleID"))

long_dist <- left_join(long_dist,
                       select(roto_meta, SampleID, Date, Location, week),
                       by = c("name" = "SampleID"))

# ok, so now I'm interested primarily in only the same locations, filter
#  this out where Location.x == Location.y:
dist_filt <- long_dist %>% 
  # can add columns if needed, not sure if I'll actually want these:
  mutate(diff = abs(as.numeric(week.x) - as.numeric(week.y)),
         week = if_else(week.y > week.x, week.y, week.x)) %>%
  # I'm going to compare everything to week 28, so include only 
  filter(Location.x == Location.y & (week.x == 28 | week.y == 28))

# plot the bray-curtis distances between samples by week:
ggplot(aes(x = week, y = value), data = dist_filt) + 
  geom_point() +
  facet_grid(Location.x~.) +
  theme_bw()
