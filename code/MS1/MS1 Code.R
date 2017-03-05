library(tidyverse)
library(raster)
library(rgdal)
library(PerformanceAnalytics)
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
d=dplyr::select(md,which(names(md)%in%vars_of_interest))

#Choose species of interest
#Here, California gymnosperm species that have basal area layers (have a CodeNum)
d=d[which(d$California==1 & d$Gymno==1 & !is.na(d$CodeNum)),]

#Clean up
rm(flam,vars_of_interest)

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

####3. Look at trait correlations####
#chart.Correlation(d[,traits_of_interest])
#chart.Correlation(d[,quants_of_interest])

####2. Import and process basal area data####
#NOTE: This data is from the Forest Service (Wilson et al. 2013; http://www.fs.usda.gov/rds/archive/Product/RDS-2013-0013/). Units are sq ft/ac

####2a: Read in all relevant GIS layers####

#Set the area of interest (AOI). Options include:
#"CA_Boundary"
#"Western_States"
AOI = readOGR(dsn="./GIS", layer="Western_States")
AOI <- spTransform(AOI, CRS("+init=epsg:5070")) #Put on Nad83/Conus Albers (EPSG = 5070) scale


#Import the basal area layers for the study species (those filtered into "d" in section #1; must have traits and basal area data). The cropping process is slow.
sppFileNames.study=paste0("s",d$CodeNum,".img") 
ba.rasters.study=list()
for(r in sppFileNames.study){ #Takes about 10 seconds per species.
      ba.rasters.study[[r]]=
            raster(paste0("../large_files/LiveBasalAreaRasters/",r)) %>%
            crop(extent(AOI))/4.356 #Convert sq ft/ac to sq m/ha
      
      #Note, the "../" directs one level up
      #"large_files" must be in the parent directory
      print(Sys.time())
      gc()
}

#Import the basal area layers for other tree species in the area (not study species). Mostly includes hardwoods, plus conifers for which we have no trait data.
sppCodeNums.other=md$CodeNum[which(md$Western==1 & !is.na(md$CodeNum) & !md$Code%in%d$Code)] 
sppFileNames.other=paste0("s",sppCodeNums.other,".img")
ba.rasters.other=list()
for(r in sppFileNames.other){ #Takes about 10 seconds per species.
      ba.rasters.other[[r]]=
            raster(paste0("../large_files/LiveBasalAreaRasters/",r)) %>%
            crop(extent(AOI))/4.356 #Convert sq ft/ac to sq m/ha
      #Note, the "../" directs one level up
      #"large_files" must be in the parent directory
      print(Sys.time())
      gc()
}

#Stack up the different rasters to get total basal area and filter out areas that are not conifer-dominated.
Sys.time() #This part takes some time; Sys.time to monitor.
ba.stack.study=stack(ba.rasters.study) #Create a raster stack of the study species
ba.stack.other=stack(ba.rasters.other) #Create a raster stack for basal area of all species
Sys.time()
ba.tot.study=overlay(ba.stack.study,fun=sum) #Create layer for cumulative basal area of study species. Takes ~10:00 
gc()
Sys.time()
ba.tot.other=overlay(ba.stack.other,fun=sum) #Create layer for cumulative basal area of other species. Takes ~15:00
gc()
Sys.time()
#hist(ba.tot.study[ba.tot.study>50]) #There are a few huge outliers here in the redwood region, may want to truncate.
#Adjust values for total basal area
ba.tot.study[ba.tot.study>120]=NA #Remove outliers with large BA (just a few pixels from the redwood region)
ba.tot.study[ba.tot.study==0]=NA #If there is no basal area from any of the focal species, set basal area to NA.
ba.tot.study[ba.tot.study/(ba.tot.study+ba.tot.other)<0.50]=NA #Remove pixels where the basal area of the study species is less than 50% of the total

#Save files to save time in future.
write_rds(ba.tot.study,"./data/RDS/ba.tot.study.RDS")
write_rds(ba.rasters.study,"./data/RDS/ba.rasters.study.RDS")
writeRaster(ba.tot.study,"../large_files/ba.tot.study.tiff")#Large file, write to parent directory.

####3. Do community-weighting of traits####
Testing