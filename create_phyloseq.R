
library(devtools)
library(qiime2R)
library(tidyverse)
library(vegan)
library(ggplot2)
library(reshape2)
library(lubridate)
library(phyloseq)


# read in files to construct a phyloseq object:

# refer to https://joey711.github.io/phyloseq/import-data.html
#   for examples of how to do this (what format the files need to
#   be in, etc.)

# import the amplicon reads (250 b paired-end, Illumina):
svs <- read_qza("its1-combined-filtered-itsx-fungi-mothur.qza")

# import taxonomy (as determined by the mothur algorithm, separate file):
taxonomy <- read_qza("taxonomy-its1-combined-mothur.qza")
taxonomy <- taxonomy$data

# manually parse taxonomy:
tax_sep <- taxonomy %>% separate(Taxon, c("Kingdom", "Phylum", "Class", "Order", 
                                          "Family", "Genus", "Species"), sep = "([;])")
# column names to allow for iteration:
to_gsub <- c("Kingdom", "Phylum", "Class", "Order", 
             "Family", "Genus", "Species")

# remove "k__" from beginning of columns:
tax_sep[to_gsub] <- lapply(tax_sep[to_gsub], gsub, 
                           pattern = ".\\_\\_", replacement = "")
# remove "unidentified" and replace with empty string (<NA>):
tax_sep[to_gsub] <- lapply(tax_sep[to_gsub], gsub, 
                           pattern = "unidentified", replacement = "")

# import sample metadata:
metadata <- read.table("ITS1_combined_metadata_2.txt", sep = "\t", header = TRUE)

metadata$Date <- as.Date(parse_date_time(metadata$Date, orders = 'dmy'))

# add in week and month columns to help with downstream analysis:
metadata$week <- as.factor(week(metadata$Date))
metadata$month <- as.factor(month(metadata$Date))
metadata <- metadata %>% 
  mutate(biweekly = case_when(is.na(week) | week == "25" ~ "NA",       # not sure why I need to specify this!!!
                              week == "26" | week == "27" ~ "26-27",
                              week == "28" | week == "29" ~ "28-29",
                              week == "30" | week == "31" ~ "30-31",
                              week == "32" | week == "33" ~ "32-33",
                              week == "34" | week == "35" ~ "34-35"))


# to make a phyloseq object, I need to ensure the otu table is 
#   in a matrix form
otus <- as.matrix(svs$data)

# so there are some inconsistencies between my taxon table and the OTU table 
#    (likely due to trying to integrate mothur output with the QIIME2 workflow), 
#    so I need to subset the unique IDs in the taxonomic file to ONLY the ones
#    that are in the OTU table:
common_ids <- intersect(rownames(otus), tax_sep[,1])
taxa_sub <- filter(tax_sep, Feature.ID %in% common_ids)
# confirm this worked:
setequal(rownames(otus), taxa_sub[,1])

# create matrix of taxonomy:
taxa <- as.matrix(taxa_sub)

# must move FeatureID column to row name and remove that column:
rownames(taxa) <- taxa[,1]
taxa2 <- as.matrix(taxa[,-1])

# give metadata rownames:
rownames(metadata) <- metadata$SampleID

# # also want to color by field vs control on the first bar chart:
metadata$Sample_or_Control <- ifelse(metadata$PosNeg == "Pos" | is.na(metadata$PosNeg),
                                 "Sample", "Control")

# merge otus, taxonomy, and metadata to make a phyloseq object:
# first, change tables into appropriate data structures:
OTU <- otu_table(otus, taxa_are_rows = TRUE)
TAX <- tax_table(taxa2)
samples <- sample_data(metadata)

# then merge them:
pseq_its1 <- phyloseq(OTU, TAX, samples)

# subset phyloseq object by metadata for separate analyses:
burk2019_its1 <- subset_samples(pseq_alpha, Year == 2019 & Sampler_Type == "Burkard")
roto_its1 <- subset_samples(pseq_alpha, Year == 2019 & Sampler_Type == "Rotorod")
burk2021_its1 <- subset_samples(pseq_alpha, Year == 2021 & Sampler_Type == "Burkard")

# save objects for use in future scripts:
saveRDS(pseq_its1, "its1_all_raw.rds")
saveRDS(burk2019_its1, "burk2019_its1_raw.rds")
saveRDS(roto_its1, "roto_its1_raw.rds")
saveRDS(burk2021_its1, "burk2021_its1_raw.rds")
