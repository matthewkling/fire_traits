library(tidyverse)

####1. Assemble trait data####

#Merge component datasets
md=Reduce(function(x, y) 
      merge(x, y, all=TRUE), 
      list(read_csv("./data/species_list.csv"), 
           read_csv("./data/ffe_traits.csv"), 
           read_csv("./data/try_traits.csv"),
           read_csv("./data/S&A_traits.csv"),
           read_csv("./data/flam_traits.csv")
           ) 
      )

vars_of_interest = c("Scientific_Name","Code","CA_ID","CodeNum",
                     "Bark.Thickness.25.4","Plant.height","Self.pruning","Flame_ht","Flame_duration")
names(md)
