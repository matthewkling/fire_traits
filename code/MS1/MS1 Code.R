library(tidyverse)

####1. Assemble trait data####

#Calculate mean flammability traits for PICO because working on the species-level
flam=read_csv("./data/flam_traits.csv")
flam[nrow(flam)+1,"Scientific_Name"] <- "Pinus_contorta"
flam[nrow(flam),c(2,4,6,8)] = 
      colMeans (flam[grep("contorta",flam$Scientific_Name),c(2,4,6,8)],na.rm=TRUE)
flam[nrow(flam),c(3,5,7,9)] = "Derived from Banwell & Varner"

#Merge component datasets
md=Reduce(function(x, y) 
      merge(x, y, all=TRUE), 
      list(read_csv("./data/species_list.csv"), 
           read_csv("./data/ffe_traits.csv"), 
           read_csv("./data/try_traits.csv"),
           read_csv("./data/S&A_traits.csv"),
           flam
           ) 
      )

#Choose variables of interest. "25.4" refers to dbh of tree in cm.
vars_of_interest = c("Scientific_Name","Code","CodeNum","California","Western","Gymno",
                     "Bark.Thickness.25.4","Plant.height","Self.pruning","Flame_ht","Flame_duration")
d=select(md,which(names(md)%in%vars_of_interest))

#Choose species of interest
#Here, California gymnosperm species that have basal area layers
d=d[which(d$California==1 & d$Gymno==1 & !is.na(d$CodeNum)),]

#Assign a species-level value to PICO for flammability traits

####2. Calculate fire-resistance score for species of interest####
#Extract the quantile of each trait for each species 
traits_of_interest = 
      c("Bark.Thickness.25.4","Plant.height","Self.pruning","Flame_ht","Flame_duration")
d$bt.quant=ecdf(d$Bark.Thickness.25.4)(d$Bark.Thickness.25.4)
d$ph.quant=ecdf(d$Plant.height)(d$Plant.height)
d$sp.quant=ecdf(d$Self.pruning)(d$Self.pruning)
d$fh.quant=ecdf(d$Flame_ht)(d$Flame_ht)
d$fd.quant=ecdf(-d$Flame_duration)(-d$Flame_duration)

#Apply (across each row) the weighted mean of the traits of interest, weighted by its completeness
quants_of_interest = 
      c("bt.quant","ph.quant","sp.quant","fh.quant","fd.quant")
wts <- #Weight each quantile variable by its completeness
      d[,quants_of_interest] %>% 
      apply(MARGIN=2,function(x) length(which(!is.na(x)))/length(x))
d$frc <-
      d[,quants_of_interest] %>%
      apply(MARGIN=1, function(x) weighted.mean(x=x, w=wts, na.rm= TRUE) )

d$fire.resistance <- #Temporary, just check out how scores look without flammability data
      rowMeans(d[,c("bt.quant","ph.quant","sp.quant")])
