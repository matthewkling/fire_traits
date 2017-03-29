##Code for manuscript on the biogeography of fire regimes.## 
##Jens Stevens; stevensjt@gmail.com##
library(tidyverse)
library(raster)
library(rgdal)
library(PerformanceAnalytics)
library(RColorBrewer)
library(viridis)
library(broom)
####1. Assemble trait data (fast)####

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
vars_of_interest = c("Scientific_Name","Code","CodeNum","California","Western","Gymno","Has_BA","Bark.Thickness.25.4","Plant.height","Self.pruning","Flame_duration","Flame_ht","Pct_consumed")
d=dplyr::select(md,which(names(md)%in%vars_of_interest))

#Choose species of interest
#Here, Western gymnosperm species that have basal area layers
#Choosing Western instead of Californian adds 4 species ("Chamaecyparis_nootkatensis", "Juniperus_scopulorum", "Larix_occidentalis", "Picea_glauca")
d=d[which(d$Western==1 & d$Gymno==1 & d$Has_BA==1),]

#Clean up
rm(flam,vars_of_interest)

####2. Calculate fire-resistance score for species of interest (fast)####
#Deal with outliers
#Don't trust the Tsuga bark thickness values from FVS FFE, too thick.
d[grep("Tsuga",d$Scientific_Name),"Bark.Thickness.25.4"] <- NA

#Extract the quantile of each trait for each species 
traits_of_interest <- 
      c("Bark.Thickness.25.4","Plant.height","Self.pruning","Flame_duration","Flame_ht","Pct_consumed")
d$bt.quant=ecdf(d$Bark.Thickness.25.4)(d$Bark.Thickness.25.4)
d$ph.quant=ecdf(d$Plant.height)(d$Plant.height)
d$sp.quant=ecdf(d$Self.pruning)(d$Self.pruning)
d$fd.quant=ecdf(-d$Flame_duration)(-d$Flame_duration) #Most resistant duration is shortest.
d$fh.quant=ecdf(d$Flame_ht)(d$Flame_ht)
d$pc.quant=ecdf(d$Pct_consumed)(d$Pct_consumed)

#Apply (across each row) the weighted mean of the traits of interest, weighted by its completeness
quants_of_interest = 
      c("bt.quant","ph.quant","sp.quant","fh.quant","fd.quant","pc.quant")
wts <- #Weight each quantile variable by its completeness
      d[,quants_of_interest] %>% 
      apply(MARGIN=2,function(x) length(which(!is.na(x)))/length(x))
d$frs <-
      d[,quants_of_interest] %>%
      apply(MARGIN=1, function(x) weighted.mean(x=x, w=wts, na.rm= TRUE) )

#write_csv(d[,c(1,8:ncol(d))],"./manuscript/tables/TableS1.csv")

####3. Look at trait correlations (fast)####
#hist(log10(d$Flame_duration))
#hist(log10(d$Bark.Thickness.25.4))

d$log_Bark.thickness=log10(d$Bark.Thickness.25.4)
d$log_Flame_duration <- log10(d$Flame_duration)
traits_of_interest_cors <- 
      c("log_Bark.thickness","Plant.height","Self.pruning",
        "Flame_ht","log_Flame_duration","Pct_consumed")
chart.Correlation(d[,traits_of_interest_cors])
#chart.Correlation(d[,quants_of_interest_cors])
dev.copy2pdf(file="./figures/MS1/FigS1_trait_correlations.pdf") 

####4. Import and process basal area data (slow)####
#NOTE: This data is from the Forest Service (Wilson et al. 2013; http://www.fs.usda.gov/rds/archive/Product/RDS-2013-0013/). Units are sq ft/ac

#Set the area of interest (AOI). Options include:
#"CA_Boundary"
#"Western_States"
AOI = readOGR(dsn="./GIS", layer="Western_States")
AOI <- spTransform(AOI, CRS("+init=epsg:5070")) #Put on Nad83/Conus Albers (EPSG = 5070) scale

#Import the basal area layers for the study species (those filtered into "d" in section #1; must have traits and basal area data). The cropping process is slow. Takes ~5 minutes
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
#Save this basal area product to work with it later. Takes 15 minutes, creates a ~700MB file in "large files"
writeRaster(stack(ba.rasters.study), "../large_files/ba.rasters.study.tif", bylayer=FALSE, format='GTiff')
#ba.rasters.study=stack("../large_files/ba.rasters.study.tif") #If reading from file already created

#Import the basal area layers for other tree species in the area (not study species). Mostly includes hardwoods, plus conifers for which we have no trait data.
sppCodeNums.other=md$CodeNum[which(md$Western==1 & md$Has_BA==1 & !is.na(md$CodeNum) & !md$Code%in%d$Code)] 
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
#Save this basal area product to work with it later. Takes 15 minutes, creates a ~700MB file in "large files". It's a stack, which means it can be treated like one when imported.
writeRaster(stack(ba.rasters.other), "../large_files/ba.rasters.other.tif", bylayer=FALSE, format='GTiff')
#ba.rasters.other=stack("../large_files/ba.rasters.other.tif") #If reading from file already created

#Stack up the different rasters to get total basal area and filter out areas that are not conifer-dominated.

#Create a raster stack of the study species.
ba.stack.study=stack(ba.rasters.study) # If running from scratch
#ba.stack.study=stack("../large_files/ba.rasters.study.tif") # If loading the basal area rasters from file

#Create a raster stack for basal area of all species
ba.stack.other=stack(ba.rasters.other) #If running from scratch
#ba.stack.other=stack("../large_files/ba.rasters.other.tif") # If loading the basal area rasters from file

Sys.time()
ba.tot.study=overlay(ba.stack.study,fun=sum) #Create layer for cumulative basal area of study species. Takes ~40:00 
gc()
Sys.time()
ba.tot.other=overlay(ba.stack.other,fun=sum) #Create layer for cumulative basal area of other species. Takes ~15:00
gc()
Sys.time()

#Adjust values for total basal area (take some areas out of the study)
ba.tot.study[ba.tot.study>120]=NA #Remove outliers with large BA (just a few pixels from the redwood region)
ba.tot.study[ba.tot.study==0]=NA #If there is no basal area from any of the focal species, set basal area to NA.
ba.tot.study[ba.tot.study/(ba.tot.study+ba.tot.other)<0.50]=NA #Remove pixels where the combined basal area of all the study species (conifers) is less than 50% of the total tree basal area
ba.tot.study[ba.tot.study<5]=NA #Filter out "sparse woodlands" by removing stands where total basal area is < 5 m2/ha

#Save total basal area layer to save time in future.
writeRaster(ba.tot.study,"../large_files/ba.tot.study.tif",overwrite=TRUE) #Large file, write to parent directory.
#ba.tot.study <- raster("../large_files/ba.tot.study.tif") #If importing from file

####5. Do community-weighting of traits (slow)####
#Takes some time.
sppFileNames.study=paste0("s",d$CodeNum,".img") 
ba.weights=fr.spp=list() #Create empty lists
Sys.time()
for(r in 1:length(sppFileNames.study)){ #For each species, calculate its fraction of stand basal area for each pixel (ba.weights), and create a new raster (fr.spp) containing a uniform single value for its trait of interest. Takes ~1 hour
      ba.weights[[r]]=overlay(ba.rasters.study[[r]],ba.tot.study,fun=function(x,y){return(x/y)}) #Calculate the fraction of stand basal area (study species) comprised by the given species
      ba.weights[[r]][is.na(ba.weights[[r]])]=0 #Any cells where the species is absent, basal area weights should be 0 rather than NA. The NA's come from dividing by NA (in unforested areas) in the overlay() function above, but the weights where trees are absent should be 0. 
      #plot(ba.weighted[[r]],main=d$Code[r])
      fr.spp[[r]] <- 
            calc(ba.weights[[r]], fun=function(x){(x*0)+d$frs[r]}) #Modify the basal area layer to create a new overlaying raster containing the fire resistance score, which will be used to do the weighting later on. 
      print(r);print(length(sppFileNames.study))
      print(Sys.time())
      gc()
}
Sys.time()

#Stack the individual basal area weights and fire resistance scores
ba.weights.stack=stack(ba.weights)
writeRaster(ba.weights.stack,"../large_files/ba.weights.stack.tif",overwrite=TRUE) #Takes 30 mins
#ba.weights.stack <- stack("../large_files/ba.weights.stack.tif") #If reading existing file
fr.stack=stack(fr.spp) 
writeRaster(fr.stack,"../large_files/fr.stack.tif") #Takes 10 mins. Resulting file is larger for some reason (1.6 GB vs 600 MB for ba.weights.stack)
##fr.stack <- stack("../large_files/fr.stack.tif") #If reading existing file

Sys.time()
fr.weighted=raster::weighted.mean(x=fr.stack,w=ba.weights.stack,na.rm=T) #Key step. Calculate the CWM fire resistance score of each pixel. Takes ~1 hr 10 minutes
Sys.time()
writeRaster(fr.weighted,"../large_files/fr.weighted.tif") #fairly fast
#fr.weighted=raster("../large_files/fr.weighted.tif") #If reading existing file

plot(fr.weighted, main=c("Fire resistance index \nweighted by species abundance"),col=rev(colorRampPalette(brewer.pal(11,"Spectral"))(100)))
plot(AOI,add=T)

#dev.copy2pdf(file="./figures/MS1/Fig1_frs_western.pdf") 
#Start Here: Replace this line with better plotting code from Matt.

####5b. Read in layers that have already been created####
#ba.weights.stack=stack("../large_files/ba.weights.stack.tif")
#ba.rasters.study=stack("../large_files/ba.rasters.study.tif")
#ba.rasters.other=stack("../large_files/ba.rasters.other.tif")
#ba.tot.study=raster("../large_files/ba.tot.study.tif")
#fr.weighted=raster("../large_files/fr.weighted.tif")

####6. Process LANDFIRE Fire Regime data (fast)####

#Load and plot Fire Return Interval data
western.fri= raster("./GIS/FRG_Rasters/Western_FRI_Clip_250m.tif") #Load Fire Return Intervals (FRI data). Fast. 250 m resolution (made this aggregation directly in ArcGIS to save time)

#Reclassify the FRI pixel values to represent the midpoint FRI
rclmat=data.frame(from=c(0:22,111,112,131,132,133),to=c(0.9:22.9,111.9,112.9,131.9,132.9,133.9),
                  becomes=c(NA,seq(3,48,by=5),seq(55,95,10),113,138,175,250,400,750,1000,rep(NA,5)))
rclmat2=data.frame(from=c(-Inf,0,2,4,6,8,11,17,20,22),to=c(0,2,4,6,8,11,17,20,22,150),
                   becomes=c(NA,5,15,25,35,50,100,200,500,NA))#Coarser bins, using this one for more even sample size. Classes are 5(0-10), 15(11-20), 25 (21-30), 35 (31-40), 50(41-60), 100 (61-150), 200 (150-300), 500 (>300).
western.fri.reclass=reclassify(western.fri,as.matrix(rclmat2),right=T)
rclmat.log=data.frame(from=c(-Inf,0,2,4,6,8,11,17,20,22),to=c(0,2,4,6,8,11,17,20,22,150),
                      becomes=log(c(NA,5,15,25,35,50,100,200,500,NA)))#Coarser bins on a log scale for plotting purposes
western.fri.reclass.log=reclassify(western.fri,as.matrix(rclmat.log),right=T)
plot(western.fri.reclass.log,col=rev(plasma(n=8)),legend.args=list(text='log(Median fire return interval)', side=4, font=2, line=2.5, cex=0.8))
plot(AOI,add=T)
#dev.copy2pdf(file="./figures/MS1/EDA/FRI_statewide.pdf") 

#Load and plot Fire Regime Group data
western.frg= raster("./GIS/FRG_Rasters/Western_FRG_Clip_250m.tif") #Load Fire Regime Groups (FRG data) Fast. 250 m resolution (made this aggregation directly in ArcGIS to save time)
rclmat.frg=data.frame(from=c(-Inf,0.5,1.5,2.5,3.5,4.5,5.5),to=c(0.5,1.5,2.5,3.5,4.5,5.5,Inf),becomes=c(NA,1,NA,3,NA,5,NA)) #Only work with FRG 1,3,5
western.frg.reclass=reclassify(western.frg,as.matrix(rclmat.frg),right=T)

plot(western.frg.reclass,col=c("darkgoldenrod2","gray","darkcyan"),legend.args=list(text='Fire Regime Group', side=4, font=2, line=2.5, cex=0.8))
plot(AOI,add=T) 
#dev.copy2pdf(file="./figures/MS1/EDA/FRG_statewide.pdf") 

####7. Set up data frame for analysis (fast)####
#Create a points layer for FRG, FRI and FRS data. The resulting objects here are large matrices with a column for x, a column for y, and a column for the extracted value
frg.pts = rasterToPoints (western.frg.reclass); gc() #Takes 0:30 sec
dimnames(frg.pts)[[2]][3]="frg"
fri.pts = rasterToPoints (western.fri.reclass); gc() #Takes 0:34 sec
dimnames(fri.pts)[[2]][3]="fri"
frs.pts = rasterToPoints (fr.weighted); gc() #Takes 0:30 sec
dimnames(frs.pts)[[2]][3]="frs"
#Make frs raster line up with other two rasters (small offset in m)
frs.pts[,1]=frs.pts[,1]+42; frs.pts[,2]=frs.pts[,2]+143

#Extract values for the FRI and FRS, at each point of the underlying FRG.
#This takes a LONG time for all western US (at least 2 hours for each) so trying alternative approach above by adjusting frs point locations.
#frg.sync=raster::extract(western.fri,frg.pts[,c(1,2)],method="simple")
#frs.sync=raster::extract(fr.weighted,frg.pts[,c(1,2)],method="simple") 

#Set up data frames for analysis
Sys.time()
sd <- Reduce(function(x, y) 
      merge(x, y, all=FALSE), #Complete cases only
      list(frs.pts,frg.pts,fri.pts
           ) 
) #Takes ~35 minutes
Sys.time()
write_rds(sd,"../large_files/sd_spatial_data_frame.RDS")
#sd <- read_rds("../large_files/sd_spatial_data_frame.RDS")
#Calculate percentiles of FRS in frequent fire (<20 year) systems
frs.ff <- sd[sd$fri<20,"frs"]
sd[sd$fri<20,"frs.ff"]=ecdf(frs.ff)(frs.ff)

#Calculate percentiles of FRS in intermediate fire (41-150 year) systems
frs.intf <- sd[between(sd$fri,41,150),"frs"]
sd[between(sd$fri,41,150),"frs.intf"]=ecdf(frs.intf)(frs.intf)

#Calculate percentiles of FRS in infrequent fire (151-300 year) systems
frs.inff <- sd[between(sd$fri,151,300),"frs"] 
sd[between(sd$fri,151,300),"frs.inff"]=ecdf(frs.inff)(frs.inff)

#Classify the mismatches (extreme 10% of FRS)
sd[which(sd$frs.ff<0.2),"mismatch"]="v.ff" #Vulnerable, frequent-fire
sd[which(sd$frs.intf<0.2),"mismatch"]="v.intf" #Vulnerable, intermediate-fire
sd[which(sd$frs.intf>0.8),"mismatch"]="r.intf" #Resistant, intermediate-fire
sd[which(sd$frs.inff>0.8),"mismatch"]="r.inff" #Resistant, infrequent-fire


####7b. Matt's alternative set up data frame for analysis (is this really faster/better??)####
##START HERE
frs <- raster("../large_files/fr.weighted.tif")
frg <- raster("./GIS/FRG_Rasters/Western_FRG_Clip_250m.tif")
fri <- raster("./GIS/FRG_Rasters/Western_FRI_Clip_250m.tif")

frs <- crop(frs, crop(fri, frs))
frg <- crop(frg, frs)
fri <- crop(fri, frs)
#ba.rasters.other=stack("../large_files/ba.rasters.other.tif")
ba.tot.other <- calc(ba.rasters.other, sum) #Takes 10 minutes
ba.prop <- ba.tot.study / (ba.tot.study + ba.tot.other)

#Reclassify fire regime data to suit our purposes
rclmat.frg <- data.frame(from=c(-Inf,0.5,1.5,2.5,3.5,4.5,5.5),
                         to=c(0.5,1.5,2.5,3.5,4.5,5.5,Inf),
                         becomes=c(NA,1,NA,3,NA,5,NA)) #Only work with FRG 1,3,5
frg <- reclassify(frg,as.matrix(rclmat.frg),right=T) 

rclmat.fri <- data.frame(from=c(-Inf,0,2,4,6,8,11,17,20,22),
                      to=c(0,2,4,6,8,11,17,20,22,150),
                      becomes=c(NA,5,15,25,35,50,100,200,500,NA))
fri <- reclassify(fri, as.matrix(rclmat.fri), right=T)


# spatial sync -- shortcut assumes that raster misalignment is negligible!
fri_tmp <- fri
frg_tmp <- frg
fri <- frs; values(fri) <- values(fri_tmp)
frg <- frs; values(frg) <- values(frg_tmp)

s <- stack(frs, fri, frg, ba.prop, ba.tot.other + ba.tot.study)
s_d <- as.data.frame(rasterToPoints(s))
names(s_d) <- c("x", "y", "frs", "fri", "frg", "ba_prop", "ba_tot")

s_d <- filter(s_d, !is.na(frs), !is.na(fri))

####8. New plots####

usa <- getData("GADM", country="USA", level=1) %>%
      spTransform(crs(ba.tot.study)) %>%
      fortify() %>%
      mutate(group=paste("usa", group))
can <- getData("GADM", country="CAN", level=1) %>%
      spTransform(crs(ba.tot.study)) %>%
      fortify() %>%
      mutate(group=paste("can", group))
mex <- getData("GADM", country="MEX", level=1) %>%
      spTransform(crs(ba.tot.study)) %>%
      fortify() %>%
      mutate(group=paste("mex", group))
borders <- rbind(usa, can, mex)

minimalism <- theme(axis.text=element_blank(), 
                    axis.title=element_blank(), 
                    axis.ticks=element_blank(),
                    panel.grid=element_blank(),
                    panel.background=element_blank(),
                    legend.position="top")

view <- coord_cartesian(xlim=range(d$x), 
                        ylim=range(d$y))

p <- ggplot() + 
      geom_polygon(data=borders, aes(long, lat, group=group),
                   fill="gray95", color=NA) +
      geom_raster(data=s_d,
                  aes(x, y, fill=frs)) +
      geom_path(data=borders, aes(long, lat, group=group),
                size=.25) +
      scale_fill_gradientn(
            colors = rev(colorRampPalette(brewer.pal(11,"Spectral"))(10))[c(1:4,7:10)]) +
      minimalism +
      view +
      guides(fill=guide_colourbar(barwidth=15)) +
      labs(fill="score (0-1)",
           title="Fire resistance index\nweighted by species abundance")

ggsave("figures/MS1/Fig1.frs.png", p, width=7, height=9, units="in")

# map of FRG
p <- ggplot() + 
      geom_polygon(data=borders, aes(long, lat, group=group),
                   fill="gray95", color=NA) +
      geom_raster(data=s_d,
                  aes(x, y, fill=factor(frg))) +
      geom_path(data=borders, aes(long, lat, group=group),
                size=.25) +
      scale_fill_manual(values=c("1"="darkgoldenrod2","3"="gray","5"="darkcyan")) +
      minimalism +
      view +
      labs(fill="Fire regime group\n")
ggsave("figures/MS1/FigS2.frg.png", p, width=7, height=9, units="in")

# map of FRI
p <- ggplot() + 
      geom_polygon(data=borders, aes(long, lat, group=group),
                   fill="gray95", color=NA) +
      geom_raster(data=s_d,
                  aes(x, y, fill=fri)) +
      geom_path(data=borders, aes(long, lat, group=group),
                size=.25) +
      scale_fill_viridis(option="B", trans="log10", breaks=c(1,3,10,30,100,300,1000)) +
      minimalism +
      view +
      guides(fill=guide_colourbar(barwidth=15)) +
      labs(fill="years\n",
           title="\nFire return interval")
ggsave("figures/MS1/FigS3.fri.png", p, width=7, height=9, units="in")

####9.Analyze and plot (fast)####
sd.sub=sd[sample(nrow(sd), nrow(sd)*0.01), ]#Subsample 1% of the df (134k cells)
sd.sub$frg[sd.sub$frg>5|sd.sub$frg==2|sd.sub$frg==4]=NA #Set certain FRG's to NA (e.g. rock, barren; FRG>5). Also remove the FRG 2's and 4's (grasslands and chaparral predominantly).
table(sd.sub$fri,sd.sub$frg)  #Check that the FRI's mostly make sense with the FRG's (most of the 35 year and less FRI's are in FRG1, most of the 50 and 100 yr FRI's are with FRG 3, and most of the 200+ yr FRI's are with FRG5). Looks good.

#Plot FRI's
#FRG.colors=c("#d7191c","#fdae61","#ffffbf","#a6d96a","#1a9641")
#Subset the huge data frame for faster plotting

sd.sub$trait=sd.sub$frs #Change to trait of interest (also change axis label)
ggplot(sd.sub)+
      geom_boxplot(aes(x=fri,y=trait,group=fri),notch=F)+ 
      geom_line(aes(x=fri,y=trait),stat="smooth",method="lm",col="black")+
      scale_x_log10(breaks=c(5,15,25,35,50,100,200,500))+
      #scale_fill_manual(values = brewer.pal(n=8,name="PuOr"))+ #Need to set fill as factor(fri)
      #scale_fill_distiller(palette="PuOr")+
      annotation_logticks(sides="b")+
      labs(y="Fire-resistance score \n(FRS; 0-1 scale)",x="Median fire return interval")+
      theme_bw()+
      theme(axis.text=element_text(size=14,color='black'),axis.title=element_text(size=18))
dev.copy2pdf(file="./figures/MS1/Fig2_Fire.resistance~MFRI.pdf") 
FRI.m=lm(frs~fri,data=sd.sub)

#Plot FRG's
sd.sub$trait=sd.sub$frs #Change to trait of interest (also change axis label)
ggplot(sd.sub)+
      geom_boxplot(aes(x=factor(frg,labels=c("1: Frequent \n low-severity", "3: Mod. frequency/ \n severity","5: Infrequent \n high-severity")),y=trait),notch=F)+ 
      #scale_fill_manual(values = brewer.pal(n=8,name="PuOr"))+ #Need to set fill as factor(fri)
      #scale_fill_distiller(palette="PuOr")+
      labs(y="Fire-resistance score \n(FRS; 0-1 scale)",x="Fire Regime Group")+
      theme_bw()+
      theme(axis.text=element_text(size=14,color='black'),axis.title=element_text(size=18))
dev.copy2pdf(file="./figures/MS1/Fig3_Fire.resistance~FRG.pdf")
FRG.m=aov(frs~factor(frg),data=sd.sub)
TukeyHSD(FRI.m)

#Plot mismatches
#Prep state boundaries for ggplot
AOI@data$id = rownames(AOI@data)
AOI.points = tidy(AOI, region="id") #Convert polygons to data frame

p.mismatches=ggplot(sd[which(!is.na(sd$mismatch)),],aes(x=x,y=y))+
      geom_raster(aes(fill=mismatch))  +
      geom_path(data=AOI.points,aes(x=long,y=lat,group=group), color="black")+
      scale_fill_manual(labels = c("resistant-infrequent", "resistant-intermediate",
                                   "vulnerable-frequent", "vulnerable-intermediate"), 
                        values = c("#76AB99","#a6d96a","#d7191c","#fdae61")) +
      labs(title="fire resistance vs historical frequency")+
      theme(axis.line=element_blank(),
            axis.text.x=element_blank(), axis.text.y=element_blank(),
            axis.ticks=element_blank(),
            axis.title.x=element_blank(), axis.title.y=element_blank(),
            panel.background=element_blank(),
            panel.border = element_rect(colour = "black", fill=NA, size=1),
            panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
            plot.background=element_blank(),
            plot.margin=unit(c(1,.6,.6,0), "cm"))

p.mismatches
dev.copy2pdf(file="./figures/MS1/Fig4_Mismatches.pdf")
