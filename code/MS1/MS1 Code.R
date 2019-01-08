##Code for manuscript on the biogeography of fire regimes.## 
##Jens Stevens; stevensjt@gmail.com##

####0. Load Libraries####
library(tidyverse) #for read_csv() etc; version 1.2.1
library(raster) #for raster(); 2.6.7
library(rgdal) # for readOGR(); version 1.2-16
library(PerformanceAnalytics) #for chart.Correlation(); version 1.5.2
library(RColorBrewer) #for brewer.pal(); version 1.1-2
library(viridis) #for plasma(); version 0.2.0
library(broom) #for tidy(); version 0.4.3
library(gridExtra) #for grid.arrange(); version 2.3

####1. Assemble trait data (fast)####

####1.1 Process flammability data###
#Calculate mean flammability traits for PICO because working on the species-level
flam <- #flammability data
      read_csv("./data/flam_traits.csv")
flam[nrow(flam)+1,"Scientific_Name"] <- "Pinus_contorta"
flam[nrow(flam),c(2,4,6,8)] <-
      colMeans (flam[grep("contorta",flam$Scientific_Name),c(2,4,6,8)],na.rm=TRUE)
flam[nrow(flam),c(3,5,7,9)] = "Derived from Banwell & Varner"

####1.2 Merge component trait datasets###
md <- #md = master data for traits
      Reduce(function(x, y) 
      merge(x, y, all=TRUE), 
      list(read_csv("./data/species_list.csv"), 
           read_csv("./data/bark_thickness_traits.csv"), 
           read_csv("./data/try_traits.csv"),
           read_csv("./data/S&A_traits.csv"),
           flam
      ) 
)

vars_of_interest <- #Identify variables of interest. "25.4" refers to dbh of tree in cm.
      c("Scientific_Name","Code","CodeNum","California","Western","Gymno",
        "Has_BA","Bark.Thickness.25.4.FOFEM2017","Plant.height",
        "Self.pruning",
        "Flame_ht", "Pct_consumed", "Flame_duration")

d <- #working data for traits
      md %>% #select variables of interest
      dplyr::select(which(names(md)%in%vars_of_interest)) %>%
      #Below, select study species, specifically western gymnosperms that have basal area data
      dplyr::filter(Western==1 & Gymno==1 & Has_BA==1)

#tidy up working data frame
d <- d[,c(1,2,8,9,10,12,13,11)]
names(d) <- c("Scientific_Name", "Code", "bt", "ph", "sp", "fh", "pc", "fd")

rm(flam,vars_of_interest) #Clean up working environment

####2. Calculate fire-resistance score for species of interest (fast)####

####2.1 Examine trait correlations###
d$log_bt <- log10(d$bt) #Log transform to normalize
d$log_fd <- log10(d$fd) #Log transform to reduce PIED outlier
d_cors <- d[,c("log_bt","ph","sp","fh","pc","log_fd")]
names(d_cors) <- 
      c("log\n(bark thickness)","plant height","self pruning",
        "flame height","percent\nconsumed", "log\n(flame duration)")
chart.Correlation(d_cors)
#dev.copy2pdf(file="./figures/MS1/FigS1_trait_correlations.pdf") 
rm(d_cors) #Clean up working environment
#d <- d[,-pmatch(c("log_bt","log_fd"),names(d))] #Remove log-transformed variables
d <- d[,-pmatch(c("log_bt"),names(d))] #Remove log-transformed bark thickness

####2.2 Flammability ordination###
#Since flame height and percent consumed are tightly correlated, calculate the first principal component of their ordination and use that.
ord <- prcomp(d[,c("fh","pc")])
PC1 <- ord$x[,"PC1"]
#summary(ord) #PC1 explains 96.7% of the variance of this trait

####2.3 Calculate the "percentile of range" each trait for each species###
#Had formerly calculated quantile (e.g. ecdf(d$Flame_ht)(d$Flame_ht)),
#But that overly-separated species for traits where actual values were tightly clustered.
d$bt.pct <- (d$bt-min(d$bt)) / diff(range(d$bt))
d$ph.pct <- (d$ph-min(d$ph)) / diff(range(d$ph))
d$sp.pct <- (d$sp-min(d$sp)) / diff(range(d$sp))
#d$fh.pct <- (d$fh-min(d$fh)) / diff(range(d$fh))
#d$pc.pct <- (d$pc -min(d$pc)) / diff(range(d$pc))
d$fh_pc.pct <- #Need "1-x" because most negative princomp scores are tallest flame lengths,
      #with largest percent consumed.
      1- (PC1-min(PC1)) / diff(range(PC1))
#d$fd.pct <- #Need "1-x" because most resistant duration is shortest.
#      1- (d$fd-min(d$fd)) / diff(range(d$fd)) 
d$fd.pct <- #Need "1-x" because most resistant duration is shortest. 
      #apply to log values to reduce PIED outlier
      1- (d$log_fd-min(d$log_fd)) / diff(range(d$log_fd)) 

#Apply (across each row) the mean of the traits of interest, to calculate FRS 
#(formerly weighted by trait completeness, but now have full dataset):
d$frs <-
      rowMeans(d[,c("bt.pct","ph.pct","sp.pct","fh_pc.pct", "fd.pct")]) 

d_t1 <- d[-c(2,9)] #Set up data frame for Table 1 (remove Code and log_fd).
d_t1[,c(2,8:13)] <- round(d_t1[,c(2,8:13)],2)
d_t1[,c(3,5:7)] <- round(d_t1[,c(3,5:7)],1)
d_t1$Scientific_Name <- gsub("_", " ",d_t1$Scientific_Name)
d_t1 <- d_t1[order(d_t1$frs, decreasing = TRUE),]
write_csv(d_t1, "./manuscript/tables/Table1.csv")
rm(ord,PC1,d_t1) #Clean up working environment

####3. Plot species rankings####
d_frs_ranking <- d[order(d$frs, decreasing = T),]
d_frs_ranking$group <- c(rep("archetypal frequent-fire conifers",5),
                         rep("frequent-fire associated species",3),
                         rep("mesic/shade-tolerant species",11),
                         rep("subalpine/arid species", 10))
d_frs_ranking$frs_vis <- round(d_frs_ranking$frs,2)
d_frs_ranking$frs_vis[c(5,6,10,11,13,14,15,16,17,18,21,23,24,25,28)] <- ""
p_frs_ranking <-
      ggplot(d_frs_ranking) +
      #geom_text(aes(x = rep(0, times = nrow(d)), 
      #              y = frs, label = frs_vis ),
      #          size = 3) +
      geom_label(aes(x = c(seq(from = 0.1, by = 0.5, length.out = 5),
                           seq(from = 0.1, by = 0.5, length.out = 3),
                           seq(from = 0.1, by = 0.5, length.out = 11),
                           seq(from = 0.1, by = 0.5, length.out = 10) ),
                     y = frs, label = Code, fill = group),
                 size = 3, hjust = 0) +
      geom_tile(aes(x=rep(2,29),y=rep(2,29),fill = group)) + #Dummy data to override legend
      #annotate("text", x = 0.3, y = 0.9, 
      #         label = "fire resistance score (frs)", hjust = 0, size = 8)+
      xlim(0, 6) + 
      labs(x = "", y = "frs", size = 18) +
      theme_bw() +
      scale_y_continuous(limits = c(0.14,0.85), breaks = seq(0.1,0.85,0.05)) +
      scale_fill_manual(
            values =rev(colorRampPalette(brewer.pal(11,"Spectral"))(10))[c(9,7,4,2)] ) +
      theme(axis.text.y = element_text(size = 10), 
            #axis.text.y = element_blank(),
            axis.title.y = element_text(size = 14),
            axis.text.x = element_blank(), legend.position = c(0.8,0.8))
ggsave("figures/MS1/FigS2_p_frs_ranking.png", 
       p_frs_ranking, width=8, height=4, units="in")
rm(d_frs_ranking,p_frs_ranking)


####4. Import and process basal area data (slow; only need to do this once)####
#NOTE: This data is from the Forest Service (Wilson et al. 2013; http://www.fs.usda.gov/rds/archive/Product/RDS-2013-0013/). Units are sq ft/ac

####4a.Import the basal area layers for the study species ###
##Study species are those filtered into "d" in section #1; must have traits and basal area data. The cropping process is slow. Takes ~5 minutes. Don't need to do this unless the study species have changed.
#Save this basal area product to work with it later. Takes 15 minutes, creates a ~700MB file in "large files"
sppFileNames.study=paste0("s",d$CodeNum,".img") 
ba.rasters.study <- sppFileNames.study %>%
      paste0("../large_files/LiveBasalAreaRasters/", .) %>%
      stack() %>%
      crop(extent(AOI)) %>%
      "/"(., 4.356)  # Convert sq ft/ac to sq m/ha

names(ba.rasters.study) <- d$Code      
writeRaster(ba.rasters.study,"../large_files/ba.rasters.study.tif", 
            bylayer=FALSE, format='GTiff', overwrite=T)

####4b. Import the basal area layers for other tree species in the area ###
##"Other species" are trees that are not study species. Mostly includes hardwoods, plus conifers for which we have no trait data. Need this to calculate fraction of total basal area that is represented by study species. Don't need to do this unless the study species have changed.
#Save this basal area product to work with it later. 
#Takes 15 minutes, creates a ~700MB file in "large files".
sppCodes.other=md$Code[which(md$Western==1 & md$Has_BA==1 & !is.na(md$CodeNum) & !md$Code%in%d$Code)] 
sppCodeNums.other=md$CodeNum[which(md$Western==1 & md$Has_BA==1 & !is.na(md$CodeNum) & !md$Code%in%d$Code)] 
sppFileNames.other=paste0("s",sppCodeNums.other,".img")
ba.rasters.other <- sppFileNames.other %>%
      paste0("../large_files/LiveBasalAreaRasters/", .) %>%
      stack() %>%
      crop(extent(AOI)) %>%
      "/"(., 4.356)# Convert sq ft/ac to sq m/ha

names(ba.rasters.other) <- sppCodes.other
writeRaster("../large_files/ba.rasters.other.tif", bylayer=FALSE, 
                  format='GTiff', overwrite=T)

#### 4c. Stack up the different rasters to get total basal area and filter out areas that are not conifer-dominated.###
# Create raster stack for the study species, and for other species.
#ba.rasters.study <- stack("../large_files/ba.rasters.study.tif") # If loading the basal area rasters from file
#ba.rasters.other <- stack("../large_files/ba.rasters.other.tif") # If loading the basal area rasters from file

#ba.tot.study <- overlay(ba.rasters.study,fun=sum) #Create layer for cumulative basal area of study species. Takes ~40:00 
# Cumulative basal area of study species, and for other species.
ba.tot.study <- sum(ba.rasters.study)
ba.tot.other <- sum(ba.rasters.other)

#Adjust values for total basal area (take some areas out of the study)
ba.tot.study[ba.tot.study>120]=NA #Remove outliers with large BA (just a few pixels from the redwood region)
ba.tot.study[ba.tot.study==0]=NA #If there is no basal area from any of the focal species, set basal area to NA.
ba.tot.study[ba.tot.study/(ba.tot.study+ba.tot.other)<0.50]=NA #Remove pixels where the combined basal area of all the study species (conifers) is less than 50% of the total tree basal area
ba.tot.study[ba.tot.study<5]=NA #Filter out "sparse woodlands" by removing stands where total basal area is < 5 m2/ha

#Save total basal area layer to save time in future.
writeRaster(ba.tot.study,"../large_files/ba.tot.study.tif",overwrite=TRUE) #Large file, write to parent directory.



####5. Do community-weighting of traits (slow)####

# Weighted mean of fire resistence. 
#ba.rasters.study <- stack("../large_files/ba.rasters.study.tif") #If importing from file; BA of each study species
#ba.tot.study <- raster("../large_files/ba.tot.study.tif") #If importing from file; total BA of study species
#Create stack of BA weights (takes 25 minutes):
ba.weights.stack <- #BA of given study species / BA of all study species
      #ba.tot.study contains NA for cells where BA of study species is <50% total tree BA
      ba.rasters.study %>% "/"(., ba.tot.study)
#writeRaster(ba.weights.stack, "../large_files/ba.weights.stack.tif", overwrite=T)
#ba.weights.stack  <- stack("../large_files/ba.weights.stack.tif")

#Multiply the vector of FRS values by the stack of BA weights (takes 20 minutes): 
#This works because d$frs is vectorized over the raster layers. Results in a value between 0.149 (for a pure PinEdu pixel) and 0.846 (for a pure SeqGig pixel) for each non-NA pixel
fr.weighted <- sum(d$frs * ba.weights.stack) #Not including na.rm = T
# this works because d$frs is vectorized over the raster layers
writeRaster(fr.weighted, "../large_files/fr.weighted.tif", overwrite=T)

#Quick plot of fire regime scores to take a look; more formal plot creation is below in #7b.
#plot(fr.weighted, main=c("Fire resistance index \nweighted by species abundance"),
#     col=rev(colorRampPalette(brewer.pal(11,"Spectral"))(100)))
#plot(AOI,add=T)
#dev.copy2pdf(file="./figures/MS1/Fig1_frs_western.pdf") 
#Start Here: Replace this line with better plotting code from Matt.






####6. Set up data to make maps and figures####
#ba.rasters.study <- stack("../large_files/ba.rasters.study.tif")
#ba.rasters.other <- stack("../large_files/ba.rasters.other.tif")
ba.tot.study <- raster("../large_files/ba.tot.study.tif")
#ba.tot.other <- calc(ba.rasters.other, sum) #Takes 10 minutes
#writeRaster(ba.tot.other,"../large_files/ba.tot.other.tif") #fairly fast
#basal area of non-study species:
ba.tot.other <- raster("../large_files/ba.tot.other.tif")
#fr.weighted <- raster("../large_files/fr.weighted.tif")

frs <- raster("../large_files/fr.weighted.tif")
frg <- raster("./GIS/FRG_Rasters/Western_FRG_Clip_250m.tif")
fri <- raster("./GIS/FRG_Rasters/Western_FRI_Clip_250m.tif")

frs <- crop(frs, crop(fri, frs))
frg <- crop(frg, frs)
fri <- crop(fri, frs)
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
s_d <- #s_d = stack of data. Slow, ~6 minutes
      as.data.frame(rasterToPoints(s))
names(s_d) <- c("x", "y", "frs", "fri", "frg", "ba_prop", "ba_tot")
s_d <- filter(s_d, !is.na(frs), !is.na(fri))
#write_rds(s_d,"../large_files/s_d.RDS") #write stack of data (fast)
#s_d <- read_rds("../large_files/s_d.RDS")

##Process to look up a specific cell and see what species are present
##in the a full raster
c <- cellFromXY(frs, matrix(c(-1136976,1787052),nrow=1))#single cell, coords from anywhere
c <- fourCellsFromXY(frs, matrix(c(-1136976,1787052),nrow=1))#four cell, coords from anywhere (ideally a 4-corner region)
extract(ba.rasters.study,c(c))
extract(ba.rasters.other,c(c))
getValues(frs)[c(c)]


####7. Fire resistance score map figure####

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

frsd <- frs %>% rasterToPoints %>% as.data.frame %>% rename(frs=fr.weighted)

view <- coord_cartesian(xlim=range(frsd$x), 
                        ylim=range(frsd$y))

p <- ggplot() + 
      geom_polygon(data=borders, aes(long, lat, group=group),
                   fill="gray95", color=NA) +
      geom_raster(data=frsd,
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

ggsave("figures/MS1/Fig2.frs.png", p, width=7, height=9, units="in") # (slow, giving errors)

####8. Supplementary map figures####
###FRG map figure
p_frg <- ggplot() + 
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
#ggsave("figures/MS1/FigS2.frg.png", p_frg, width=7, height=9, units="in")

###FRI map figure
p_fri <- ggplot() + 
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
#ggsave("figures/MS1/FigS3.fri.png", p, width=7, height=9, units="in")



####9.Analyze and plot frs-fire regime relationships and "mismatches" (fast)####
sd.sub <- #Randomly subsample 1% of the df (134k [~119k?] cells)
      s_d[sample(nrow(s_d), nrow(s_d)*0.01), ]
sd.sub$frg[sd.sub$frg>5|sd.sub$frg==2|sd.sub$frg==4]=NA 
#Set certain FRG's to NA (e.g. rock, barren; FRG>5). Also remove the FRG 2's and 4's (grasslands and chaparral predominantly).
table(sd.sub$fri,sd.sub$frg)  
#Check that the FRI's mostly make sense with the FRG's (most of the 35 year and less FRI's are in FRG1, most of the 50 and 100 yr FRI's are with FRG 3, and most of the 200+ yr FRI's are with FRG5). Looks good.

#Plot frs~fri (fire return interval)
sd.sub$trait=sd.sub$frs #Change to trait of interest (also change axis label)
p_fri_frs <-
      ggplot(sd.sub)+
      geom_boxplot(aes(x=fri,y=trait,group=fri),notch=F)+ 
      geom_line(aes(x=fri,y=trait),stat="smooth",method="lm",col="black")+
      scale_x_log10(breaks=c(5,15,25,35,50,100,200,500))+
      annotation_logticks(sides="b")+
      labs(y="FRS",x="Median fire return interval")+
      theme_bw()+
      theme(axis.text=element_text(size=12,color='black'),axis.title=element_text(size=14))
#dev.copy2pdf(file="./figures/MS1/Fig2_Fire.resistance~MFRI.pdf") 
FRI.m=lm(frs~fri,data=sd.sub)

#Plot Plot frs~frg (fire regime group)
sd.sub$trait=sd.sub$frs #Change to trait of interest (also change axis label)
p_frg_frs <-
      ggplot(sd.sub[!is.na(sd.sub$frg),])+
      geom_boxplot(aes(x=factor(frg,labels=c("1: Frequent \n low-severity", 
                                             "3: Mod. frequency/ \n severity",
                                             "5: Infrequent \n high-severity")),
                       y=trait),notch=F)+ 
      labs(y="FRS",x="Fire Regime Group")+
      theme_bw()+
      theme(axis.text=element_text(size=12,color='black'),axis.title=element_text(size=14))
#dev.copy2pdf(file="./figures/MS1/Fig3_Fire.resistance~FRG.pdf")
#FRG.m=aov(frs~factor(frg),data=sd.sub)
#TukeyHSD(FRG.m)
p_S4 <- grid.arrange(p_frg_frs,p_fri_frs,ncol=1)
ggsave("figures/MS1/FigS4.png", p_S4, width=5, height=7, units="in")


#Analyze and plot mismatches
#Calculate percentiles of FRS in frequent fire (<20 year) systems
frs.ff <- s_d[s_d$fri<20,"frs"]
s_d[s_d$fri<20,"frs.ff"]=ecdf(frs.ff)(frs.ff)

#Calculate percentiles of FRS in intermediate fire (41-150 year) systems
frs.intf <- s_d[between(s_d$fri,41,150),"frs"]
s_d[between(s_d$fri,41,150),"frs.intf"]=ecdf(frs.intf)(frs.intf)

#Calculate percentiles of FRS in infrequent fire (151-300 year) systems
frs.inff <- s_d[between(s_d$fri,151,300),"frs"] 
s_d[between(s_d$fri,151,300),"frs.inff"]=ecdf(frs.inff)(frs.inff)

#Classify the mismatches (extreme 10% of FRS)
s_d[which(s_d$frs.ff<0.2),"mismatch"]="s.ff" #Sensitive, frequent-fire
s_d[which(s_d$frs.intf<0.2),"mismatch"]="s.intf" #Sensitive, intermediate-fire
s_d[which(s_d$frs.intf>0.8),"mismatch"]="r.intf" #Resistant, intermediate-fire
s_d[which(s_d$frs.inff>0.8),"mismatch"]="r.inff" #Resistant, infrequent-fire

p_mismatches <- 
      ggplot() + 
      geom_polygon(data=borders, aes(long, lat, group=group),
                   fill="gray95", color=NA) +
      geom_raster(data=s_d[which(!is.na(s_d$mismatch)),],
                  aes(x, y, fill=mismatch)) +
      geom_path(data=borders, aes(long, lat, group=group),
                size=.25) +
      #scale_fill_viridis(option="B", trans="log10", breaks=c(1,3,10,30,100,300,1000)) +
      scale_fill_manual(labels = c("resistant-infrequent", "resistant-intermediate",
                                   "sensitive-frequent", "sensitive-intermediate"), 
                        values = c("#76AB99","#a6d96a","#d7191c","#fdae61")) +
      minimalism +
      view +
      #guides(fill=guide_colourbar(barwidth=15)) +
      labs(fill="category \n(FRS-FRI) \n",
           title="\nMismatches between FRS and FRI")

#Plot frs~fri (fire return interval with coloring for inset figure)
sd.sub <- #Randomly subsample 1% of the df (134k [~119k?] cells)
      s_d[sample(nrow(s_d), nrow(s_d)*0.01), ]
sd.sub$frg[sd.sub$frg>5|sd.sub$frg==2|sd.sub$frg==4]=NA 
#Set certain FRG's to NA (e.g. rock, barren; FRG>5). Also remove the FRG 2's and 4's (grasslands and chaparral predominantly).
p_fri_frs_inset <-
      ggplot(sd.sub)+
      geom_boxplot(aes(x=fri,y=frs,group=fri),notch=F)+ 
      geom_rect(aes(xmin = 1, xmax = 20, 
                    ymin = range(s_d[which(s_d$mismatch=="s.ff"),"frs"])[1],
                    ymax = range(s_d[which(s_d$mismatch=="s.ff"),"frs"])[2] ),
                color = "#d7191c", fill = NA) +
      geom_rect(aes(xmin = 41, xmax = 150, 
                    ymin = range(s_d[which(s_d$mismatch=="s.intf"),"frs"])[1],
                    ymax = range(s_d[which(s_d$mismatch=="s.intf"),"frs"])[2] ),
                color = "#fdae61", fill = NA) +
      geom_rect(aes(xmin = 41, xmax = 145, 
                    ymin = range(s_d[which(s_d$mismatch=="r.intf"),"frs"])[1],
                    ymax = range(s_d[which(s_d$mismatch=="r.intf"),"frs"])[2] ),
                color = "#a6d96a", fill = NA) +
      geom_rect(aes(xmin = 155, xmax = 300, 
                    ymin = range(s_d[which(s_d$mismatch=="r.inff"),"frs"])[1],
                    ymax = range(s_d[which(s_d$mismatch=="r.inff"),"frs"])[2] ),
                color = "#76AB99", fill = NA) +
      scale_x_log10(breaks=c(5,15,25,35,50,100,200,500))+
      annotation_logticks(sides="b")+
      labs(y="FRS",x="Median fire return interval")+
      theme_bw()+
      theme(axis.text=element_text(size=12,color='black'),axis.title=element_text(size=14))
ggsave("figures/MS1/FigS6.mismatches_inset.png", p_fri_frs_inset, width=5, height=4.5, units="in")

###Matt, can you figure out how to make the above figure inset into the below figure?
ggsave("figures/MS1/FigS6.mismatches.png", p_mismatches, width=7, height=9, units="in")