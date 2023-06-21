
library(devtools)
library(qiime2R)
library(tidyverse)
library(vegan)
library(ggplot2)
library(reshape2)
library(lubridate)
library(phyloseq)


# There are lots of ways to start analyzing amplicon data, and one important
#   thing to do is to investigate contamination in negative controls. I have
#   used the 'decontam' package to explore this previously (see tutorial at:
#   https://benjjneb.github.io/decontam/vignettes/decontam_intro.html ).
#   
# For the purposes of this analysis I have skipped this step so we'll go 
#   directly to analysis of the data.


# ALPHA DIVERSITY ANALYSIS ------------------------------------------------

# RAREFACTION:
# 
# N.B.: rarefaction is a contentious issue, and is highly criticized in some contexts
#   because you end up throwing out reads that could be informative;
#   nevertheless, it is still common practice and provides ONE approach for 
#   comparing samples that come from differing read depths.

# First, let's see what the read depths of each sample are:
sample_counts <- data.frame(sample_sums(pseq_alpha))
sample_counts <- rownames_to_column(sample_counts)

# make a histogram of read counts
ggplot(data = sample_counts, aes(sample_sums.pseq_alpha.)) +
  geom_histogram()

# We can see that there are a fair number of samples with reads < 250,000.
# We can also look at the data frame itself and see what the distribution 
# of samples is like:
head(sample_counts[order(sample_counts$sample_sums.pseq_alpha.),], 20)

# for now, we'll be conservative and keep samples with more than 15,000 reads:
to_keep <- filter(sample_counts, sample_sums.pseq_alpha. > 15000)$rowname

# and then filter the phyloseq object to include only these samples:
pseq_alpha_filt <- subset_samples(pseq_alpha, SampleID %in% to_keep)

# Next, rarefy:
# N.B.: this function rarefies to the depth of the lowest sample; samples WITH
#   replacement
pseq_rare <- rarefy_even_depth(pseq_alpha_filt, rngseed = 314)
# 88980 OTUs removed because they were no longer present after subsampling

# For purposes of this file, I'm only going to look at (i) ALL samples and (ii)
#   samples collected by the rotorod samplers:
roto_rare <- subset_samples(pseq_rare, Year == 2019 & Sampler_Type == "Rotorod")

saveRDS(pseq_rare, "its1_all_rare.rds")
saveRDS(roto_rare, "roto_its1_rare.rds")

