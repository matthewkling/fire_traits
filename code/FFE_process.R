
# Matthew Kling
# January 2016

# this script cleans and merges FFE csv files into a single table

library(dplyr)
library(tidyr)
library(stringr)


# load data from single-trait tables
f <- list.files("e:/fire_traits/input_data/FFE/03_traits", full.names=T)
loadUnique <- function(x){
        d <- read.csv(x)
        d <- unique(d)
        return(d)
}
d <- lapply(f, loadUnique)

# merge traits into a one table
r1 <- full_join(d[[1]], d[[2]])
r2 <- full_join(d[[3]], d[[5]])
r <- full_join(r1, r2)
r <- full_join(d[[4]], r, by=c("Common.Name"="Species"))
names(r) <- c("common_name", "scientific_name", "bark_thickness_multiplier", "decay_class", "leaf_longevity", "wood_density")

# cleanup
avg <- function(x) mean(na.omit(as.numeric(as.character(x))))
r <- arrange(r, scientific_name) %>%
        mutate_each(funs(str_trim), common_name, scientific_name) %>%
        mutate(common_name = tolower(common_name),
               scientific_name = sub("sp.", "spp.", scientific_name),
               scientific_name = sub("\\.\\.", ".", scientific_name)) %>%
        distinct() %>%
        filter(!is.na(scientific_name), # drop unknown spp and spp with zero relevant trait values
               !(is.na(bark_thickness_multiplier) & is.na(wood_density) & is.na(decay_class) & is.na(leaf_longevity)) ) %>%
        select(scientific_name, common_name, bark_thickness_multiplier, wood_density, decay_class, leaf_longevity) %>%
        group_by(scientific_name) %>%
        summarize_each(funs(avg), -scientific_name, -common_name) # collapse multiple observations into a single species-wide mean

# save
write.csv(r, "e:/fire_traits/derived_data/FFE_traits_clean.csv", row.names=F)
