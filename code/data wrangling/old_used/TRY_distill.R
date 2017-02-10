
# Matthew Kling
# January 2016

# this script distills raw TRY dataset into usable form covering relevant species

library(data.table)
library(dplyr)
library(tidyr)


# list of tree species in FFE dataset
spp <- read.csv("E:/fire_traits/derived_data/FFE_traits_clean.csv", stringsAsFactors=F)

# read in TRY data
tryfile <- "E:/fire_traits/input_data/TRY/1645_08022016124158/1645.txt"
d <- fread(tryfile)

# filter and summarize
d <- tbl_df(d) %>%
      filter(SpeciesName %in% spp$scientific_name,
             !is.na(TraitID)) %>%
      select(SpeciesName, TraitName, TraitID, OriglName, OrigValueStr, StdValue) %>%
      select(SpeciesName, TraitName, StdValue) %>%
      group_by(SpeciesName, TraitName) %>%
      summarize(StdValue = mean(na.omit(StdValue))) %>%
      filter(is.finite(StdValue)) %>%
      spread(TraitName, StdValue)

# only keep traits with values for at least 25 species
trait_counts <- apply(d, 2, function(x)length(na.omit(x)))
d <- d[,names(trait_counts[trait_counts>=25])]

# export csv
write.csv(d, "E:/fire_traits/derived_data/TRY_traits_clean.csv", row.names=F)

