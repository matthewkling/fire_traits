#Jens Stevens stevensjt@gmail.com
#2/10/17
#This code takes information from the spreadsheet and saves it as csv files

library(googlesheets)
library(tidyverse)
library(readxl)
gs_auth(new_user=TRUE) #Paste the link into a browser, and enter the auth code that the browser takes you to.
#my_sheets <- gs_ls()

####1. Read in data from spreadsheet####
#1a
traits_full=gs_title("Traits") #Read in full spreadsheet
master=gs_read(traits_full, ws="Master") #Identify "master" spreadsheet from which to make constituent spreadsheets

####2. Write species list spreadsheet####
d= master %>%
      select(c(ID,CA_ID,Code,Scientific_Name,CodeNum,Generic,Subsp,Western,California,Gymno)) %>%
#write_csv(d, "./data/species_list.csv") #CAREFUL, this overwrites existing file.

####3. Write FFE data spreadsheet####
#Original web scraping to create the "master" sheet was done by David, in file FFE_process.R
d= master %>%
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
#write_csv(d,"./data/ffe_traits.csv") #CAREFUL, this overwrites existing file.

####4. Write TRY data spreadsheet
#Original web scraping was done by David, in file TRY_distill.R
TRY_traits=names(master)[17:30]
TRY_traits.Source=paste(TRY_traits,"Source",sep=".")
tmp.source=master
tmp.source[,TRY_traits.Source]=NA #Add "Source" columns to temporary df
tmp.source[,TRY_traits.Source]=ifelse(is.na(tmp.source[,TRY_traits]),NA,"TRY") #Assign TRY as source to values that exist.
d= tmp.source %>%
      select(c(ID,Code,Scientific_Name,pmatch(sort(c(TRY_traits,TRY_traits.Source)),names(tmp.source) ) ) ) 

#write_csv(d,"./data/try_traits.csv") #CAREFUL, this overwrites existing file.

####5. Write data from Schwilk et al with manual supplements####
#Serotiny and self-pruning data
d= master %>%
      mutate(Serotiny.Source=ifelse(is.na(Serotiny),NA,
                                    ifelse(grepl("Pinus", Scientific_Name), "Schwilk&Ackerly", "FEIS") ),
             Self.pruning.Source=ifelse(is.na(Self.pruning),NA,
                                        ifelse(grepl("Pinus", Scientific_Name), "Schwilk&Ackerly", "FEIS") ) ) %>%
      
      select(c(ID,Code,Scientific_Name, Serotiny, Serotiny.Source, Self.pruning, Self.pruning.Source) ) 
#write_csv(d,"./data/S&A_traits.csv") #CAREFUL, this overwrites existing file

####6. Read and write data from Varner et al.####
#Flammability data
d= read_excel("./data/raw data/FlammabilityData.xlsx", sheet = "Data")
d$Species=gsub(" ","_",d$Species)
names(d)=gsub(" ","_",names(d))
names(d)[2]="Scientific_Name"
d[,paste0(names(d)[3:6],".Source")]=d$Source
d= d[,-which(names(d)=="Source")]
d= d[,-which(names(d)=="CODE")]
d= d[-which(d$Scientific_Name=="Pinus_washoensis"),]
d= d %>% select(c(Scientific_Name, pmatch(sort(names(d[,2:ncol(d) ] ) ), names(d) ) ) )
#write_csv(d,"./data/flam_traits.csv") #CAREFUL, this overwrites existing file
