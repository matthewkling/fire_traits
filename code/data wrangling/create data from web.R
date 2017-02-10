#Jens Stevens stevensjt@gmail.com
#2/10/17
#This code takes information from the spreadsheet and saves it as csv files

library(googlesheets)
library(tidyverse)
gs_auth(new_user=TRUE) #Paste the link into a browser, and enter the auth code that the browser takes you to.
#my_sheets <- gs_ls()

####1. Read in data from spreadsheet####
#1a
traits_full=gs_title("Traits") #Read in full spreadsheet
master=gs_read(traits_full, ws="Master") #Identify "master" spreadsheet from which to make constituent spreadsheets

####2. Write species list spreadsheet####
master %>%
      select(c(ID,CA_ID,Code,Scientific_Name,CodeNum,Generic,Subsp,Western)) %>%
      write_csv("./data/species_list.csv")

####3. Write FFE data spreadsheet####
#Original web scraping was done by David, in file FFE_process.R
master %>%
      mutate(Bark.Thickness.25.4=Bark.Thickness.Multiplier*25.4,
             Bark.Thickness.Source=ifelse(is.na(Bark.Thickness.Multiplier),NA,"FVS FFE"),
             Wood.density.Source=ifelse(is.na(Wood.density),NA,"FVS FFE"),
             Decay.class.Source=ifelse(is.na(Decay.class),NA,"FVS FFE"),
             Leaf.longevity.Source=ifelse(is.na(Leaf.longevity),NA,"FVS FFE")
             ) %>%
      
      select(c(ID,Code,Scientific_Name,In.FFE,Bark.Thickness.25.4,Bark.Thickness.Source,
               Wood.density,Wood.density.Source,
               Decay.class,Decay.class.Source,
               Leaf.longevity,Leaf.longevity.Source)) %>%
      write_csv("./data/ffe_traits.csv")

####4. Write TRY data spreadsheet
#Original web scraping was done by David, in file FFE_process.R
#START HERE
TRY_traits=names(master)[16:30]
TRY_traits.Source=paste(TRY_traits,"Source",sep=".")
#rbind() Just add the "source" columns here using rbind (need to google)
master %>%
      mutate() %>%
      
      select(c(ID,Code,Scientific_Name,In.FFE,Bark.Thickness.25.4,Bark.Thickness.Source,
               Wood.density,Wood.density.Source,
               Decay.class,Decay.class.Source,
               Leaf.longevity,Leaf.longevity.Source)) %>%
      write_csv("./data/ffe_traits.csv")
#Nonstandard evaluation vignette for dplyr