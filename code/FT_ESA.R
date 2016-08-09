#Jens Stevens stevensjt@gmail.com
#2/21/16
#Spatial analyses

library(sp)
library(raster)
library(rgdal)
library(googlesheets)
library(vegan)
library(viridis)
library(RColorBrewer)

####1. Import and process trait data####
gs_auth(new_user=TRUE) #Paste the link into a browser, and enter the auth code that the browser takes you to.
traits_full=gs_title("Traits")
master=gs_read(traits_full, ws="Master") #Read in master data from spreadsheet
#Subset the data to California conifers
subs=which(master$California==1 & master$Gymno==1)
d=master[subs,]
#Remove some outliers:
d=d[-pmatch(c("Tor_Cal","Tsu_Mer", "Tsu_Het"),d$Code),]
#Tor_Cal has huge seed mass; Tsu_Het and Tsu_Mer have *really* thick bark
#Remove the subspecies:
d=d[which(d$Subsp<1),]
#Select species that have corresponding basal area layers ("CodeNum")
d=d[!is.na(d$CodeNum),]
#Select just the variables most relevatnt to fire:
d=d[,c(1:5,which(names(d)%in%c("Bark.thickness","Plant.height","Serotiny","Self.pruning")))]
#Not using seed dry mass because it's complicated; really tall trees can have small seeds (e.g. the wind dispersed redwood/sequoia) but also really tall trees can have larger trees (e.g. the fire resistant pines). Plus it's not in Schwilk and Ackerly.
#d$log.Seed.dry.mass=log(d$Seed.dry.mass) #Deprecated

####1a. Correlations####

d.ord=d[,c("Bark.thickness","Plant.height","Self.pruning")] #Simplified dataset for ordinations and correlations
cor(d.ord,use="complete.obs") #Check out correlations
#Get R2 values:
summary(lm(d$Bark.thickness~d$Plant.height))$r.squared #0.41
summary(lm(d$Self.pruning~d$Plant.height))$r.squared #0.38
summary(lm(d$Bark.thickness~d$Self.pruning))$r.squared #0.30
rownames(d.ord)=d[which(complete.cases(d.ord)),"Code"]$Code

####1b. Ordination####
vare.dis <- vegdist(d.ord)
vare.mds0 <- monoMDS(vare.dis)
ordiplot(vare.mds0, type = "t") #Plot 1
vare.mds <- metaMDS(d.ord, trace = FALSE)
plot(vare.mds, type = "t") #Plot 2
#Not going with ordination for the fire tolerance axis. Instead going with the quantile analysis below.

####1c. Fire tolerance ranking####
#Calculate quantile of each species
quant_calc=function(d.fun,var){
  v=d.fun[which(names(d)==var)]
  v.range=d[,which(names(d)==var)][[1]]
  return(ecdf(v.range)(v))
}
d$bt.quant=apply(d,1,quant_calc,var="Bark.thickness")
d$ph.quant=apply(d,1,quant_calc,var="Plant.height")
d$sp.quant=apply(d,1,quant_calc,var="Self.pruning")
d$fire.resistance=rowMeans(d[,c("bt.quant","ph.quant","sp.quant")])

####2. Import and process basal area data, do community-weighting of traits####
#NOTE: This data is from the Forest Service (Wilson et al. 2013; http://www.fs.usda.gov/rds/archive/Product/RDS-2013-0013/). Units are sq ft/ac

####2a: Read in all relevant GIS layers####
sppFileNames.study=paste0("s",d$CodeNum,".img")#Set the basal area filenames for the study species (those filtered into "d" in section #1; must have traits and basal area data)
ba.rasters.study=list()
for(r in sppFileNames.study){
  ba.rasters.study[[r]]=
    raster(paste0("./GIS/LiveBasalAreaRasters/",r))
}
sppCodeNums.other=master$CodeNum[which(master$California==1 & !is.na(master$CodeNum) & !master$Code%in%d$Code)] #Get other species codes (not study species) that could be present in the study area. Mostly includes hardwoods.
sppFileNames.other=paste0("s",sppCodeNums.other,".img")
ba.rasters.other=list()
for(r in sppFileNames.other){
  ba.rasters.other[[r]]=
    raster(paste0("./GIS/LiveBasalAreaRasters/",r))
}

####2b: Define extent. Deprecated####
#One way to define extent is here: www.latlong.net
MapExtent=data.frame("Longitude"=c(-126.947510,-113.796387),"Latitude"=c(42.000325,32.768800))#California
#MapExtent=data.frame("Longitude"=c(-120.0423938,-119.9562198),"Latitude"=c(39.0890715,38.8980169))#Tahoe
coordinates(MapExtent)=c("Longitude","Latitude") 
sp::proj4string(MapExtent)= CRS("+proj=longlat +ellps=WGS84 +datum=WGS84") 
MapExtent=spTransform(MapExtent,raster::crs(ba.rasters.study[[1]])) #Reproject extent to be same as basal area raster layers. *The raster layers are Nad83/Conus Albers (EPSG = 5070)
SP.Extent=as(extent(MapExtent),'SpatialPolygons')#Reclassify the extent as a SpatialPolygons layer for cropping
crs(SP.Extent)=crs(ba.rasters.study[[1]]) #Assign coordinate system
#Deprecated above#

####2c: Crop all relevant basal area rasters to the extent of interest####
#Takes some time. Originally these rasters are for the continental US.
ba.rasters.study.crop=list() #Set up empty list for basal area rasters to be subset.
Sys.time()
for(r in 1:length(ba.rasters.study)){#For each study species, crop to the area of interest and convert to m2/ha from ft2/ac
  ba.rasters.study.crop[[r]] = crop(ba.rasters.study[[r]], extent(CA))
  ba.rasters.study.crop[[r]] = ba.rasters.study.crop[[r]]/4.356
} #Takes ~2:00 mins 
Sys.time()
ba.rasters.other.crop=list()
for(r in 1:length(ba.rasters.other)){#For all other species (mostly hardwoods), crop to the area of interest and convert to m2/ha
  ba.rasters.other.crop[[r]] = crop(ba.rasters.other[[r]], extent(CA))
  ba.rasters.other.crop[[r]] = ba.rasters.other.crop[[r]]/4.356
} #Takes ~1:00
Sys.time()

####2d: Stacking and filtering the different species' basal area rasters#### 
#Takes some time.
ba.stack.study=stack(ba.rasters.study.crop) #Create a raster stack of the study species
ba.stack.other=stack(ba.rasters.other.crop) #Create a raster stack for basal area of all species
Sys.time()
ba.tot.study=overlay(ba.stack.study,fun=sum) #Create layer for cumulative basal area of study species. Takes ~2:00 
Sys.time()
ba.tot.other=overlay(ba.stack.other,fun=sum) #Create layer for cumulative basal area of other species. Takes ~1:30
Sys.time()
#hist(ba.tot.study[ba.tot.study>50]) #There are a few huge outliers here in the redwood region, may want to truncate.
#Adjust values for total basal area
ba.tot.study[ba.tot.study>120]=NA #Remove outliers with large BA (just a few pixels from the redwood region)
ba.tot.study[ba.tot.study==0]=NA #If there is no basal area from any of the focal species, set basal area to NA.
ba.tot.study[ba.tot.study/(ba.tot.study+ba.tot.other)<0.50]=NA #Remove pixels where the basal area of the study species is less than 50% of the total

####2e: Weighting the different traits by their relative abundance in the basal area data####
#Takes some time.
ba.weighted=bt.spp=ph.spp=sp.spp=fr.spp=list() #Create empty lists
Sys.time()
for(r in 1:length(sppFileNames.study)){ #For each species, calculate its fraction of stand basal area for each pixel, and create a new raster containing its trait of interest. Takes ~1:00
  ba.weighted[[r]]=overlay(ba.rasters.study.crop[[r]],ba.tot.study,fun=function(x,y){return(x/y)}) #Calculate the fraction of stand basal area (study species) comprised by the given species
  ba.weighted[[r]][is.na(ba.weighted[[r]])]=0 #Any cells where the species is absent, basal area weights should be 0 rather than NA. The NA's come from dividing by NA in the overlay() function above, but the weights where trees are absent should be 0. 
  #plot(ba.weighted[[r]],main=d$Code[r])
  bt.spp[[r]]=ph.spp[[r]]=sp.spp[[r]]=fr.spp[[r]]=ba.weighted[[r]] #Create a new raster for each species trait. The initial values will be overwritten.
  bt.spp[[r]][]=d$Bark.thickness[r] #Assign every pixel the given species' bark thickness. 
  ph.spp[[r]][]=d$Plant.height[r] #Assign every pixel the given species' height.
  sp.spp[[r]][]=d$Self.pruning[r] #Assign every pixel the given species' height
  fr.spp[[r]][]=d$fire.resistance[r] #Assign every pixel the given species' fire resistance index
}
Sys.time()

ba.weighted.stack=stack(ba.weighted)
bt.stack=stack(bt.spp)
ph.stack=stack(ph.spp)
sp.stack=stack(sp.spp)
fr.stack=stack(fr.spp)

Sys.time()
bt.weighted=raster::weighted.mean(x=bt.stack,w=ba.weighted.stack,na.rm=T) #Calculate the mean bark thickness of each pixel, weighted by the relative abundance of that species in the stand (=Community Weighted Mean or CWM thickness). Takes ~6 minutes
Sys.time()
ph.weighted=raster::weighted.mean(x=ph.stack,w=ba.weighted.stack,na.rm=T) #Calculate the CWM plant height of each pixel Takes ~5 minutes
Sys.time()
sp.weighted=raster::weighted.mean(x=sp.stack,w=ba.weighted.stack,na.rm=T) #Calculate the CWM self pruning score of each pixel Takes ~5 minutes
Sys.time()
fr.weighted=raster::weighted.mean(x=fr.stack,w=ba.weighted.stack,na.rm=T) #Calculate the CWM fire resistance index of each pixel Takes ~6 minutes
Sys.time()

plot(fr.weighted, main=c("Fire resistance index \nweighted by species abundance"),col=rev(colorRampPalette(brewer.pal(11,"Spectral"))(100)))
plot(CA,add=T)
dev.copy2pdf(file="./figures/ESA/Fire.resistance statewide.pdf") 

####3. Process Fire Regime data####
#3a: Load trait data
#trait.weighted=bt.weighted #Adjust for a given trait. Deprecated

#3b: Masking/cropping of various layers (takes a long time)
 
#Load and plot Fire Return Interval data
ca.fri= raster("./GIS/FRG_Rasters/CA_FRI_Clip_250m.tif") #Load Fire Return Intervals (FRI data). Fast. 250 m resolution (made this aggregation directly in ArcGIS to save time)
#Reclassify the FRI pixel values to represent the midpoint FRI
rclmat=data.frame(from=c(0:22,111,112,131,132,133),to=c(0.9:22.9,111.9,112.9,131.9,132.9,133.9),
  becomes=c(NA,seq(3,48,by=5),seq(55,95,10),113,138,175,250,400,750,1000,rep(NA,5)))
rclmat2=data.frame(from=c(-Inf,0,2,4,6,8,11,17,20,22),to=c(0,2,4,6,8,11,17,20,22,150),
  becomes=c(NA,5,15,25,35,50,100,200,500,NA))#Coarser bins, using this one for more even sample size.
ca.fri.reclass=reclassify(ca.fri,as.matrix(rclmat2),right=T)
rclmat.log=data.frame(from=c(-Inf,0,2,4,6,8,11,17,20,22),to=c(0,2,4,6,8,11,17,20,22,150),
                   becomes=log(c(NA,5,15,25,35,50,100,200,500,NA)))#Coarser bins on a log scale for plotting purposes
ca.fri.reclass.log=reclassify(ca.fri,as.matrix(rclmat.log),right=T)
plot(ca.fri.reclass.log,col=rev(plasma(n=8)),legend.args=list(text='log(Median fire return interval)', side=4, font=2, line=2.5, cex=0.8))
plot(CA,add=T)
dev.copy2pdf(file="./figures/ESA/FRI_statewide.pdf") 

#Load Fire Regime Group data (extends beyond California):
ca.frg= raster("./GIS/FRG_Rasters/CA_FRG_Clip_250m.tif") #Load Fire Regime Groups (FRG data) Fast. 250 m resolution (made this aggregation directly in ArcGIS to save time)
Sys.time()
rclmat.frg=data.frame(from=c(-Inf,0.5,1.5,2.5,3.5,4.5,5.5),to=c(0.5,1.5,2.5,3.5,4.5,5.5,Inf),becomes=c(NA,1,NA,3,NA,5,NA)) #Only work with FRG 1,3,5
ca.frg.reclass=reclassify(ca.frg,as.matrix(rclmat.frg),right=T) #Takes 5 minutes
Sys.time()

plot(ca.frg.reclass,col=rev(viridis(n=3)),legend.args=list(text='Fire Regime Group', side=4, font=2, line=2.5, cex=0.8))
plot(CA,add=T)
dev.copy2pdf(file="./figures/ESA/FRG_statewide.pdf") 


####4.Analyze the trait data in comparison to the fire regime data.####
#Create a points layer for FRG and FRI data, to set up easy extraction of underlying trait values. The resulting objects here are large matrices with a column for x, a column for y, and a column for the extracted value
fri.pts = rasterToPoints (ca.fri.reclass)

#Extract values for the trait of interest, at each point of the underlying FRI.
bt.fri=raster::extract(bt.weighted,fri.pts[,c(1,2)],method="simple")
ph.fri=raster::extract(ph.weighted,fri.pts[,c(1,2)],method="simple")
sp.fri=raster::extract(sp.weighted,fri.pts[,c(1,2)],method="simple")
fr.fri=raster::extract(fr.weighted,fri.pts[,c(1,2)],method="simple")
frg.fri=raster::extract(ca.frg.reclass,fri.pts[,c(1,2)],method="simple") #Fire regime group at the FRI resolution. Takes 7 minutes but eliminates the need for the slow resampling above.

#Set up data frames for analysis
FRI.df= data.frame(fri=fri.pts[,3],frg=frg.fri,bt=bt.fri,ph=ph.fri,sp=sp.fri,fr=fr.fri)
FRI.df=FRI.df[complete.cases(FRI.df),]#1.8 million cells with both FRI and BT (similar for other traits)
FRI.df$frg[FRI.df$frg>5|FRI.df$frg==2|FRI.df$frg==4]=NA #Set certain FRG's to NA (e.g. rock, barren). Also remove the FRG 2's and 4's (grasslands and chaparral predominantly).
table(FRI.df$fri,FRI.df$frg)  #Check that the FRI's mostly make sense with the FRG's (most of the 35 year and less FRI's are in FRG1, most of the 50 and 100 yr FRI's are with FRG 3, and most of the 200+ yr FRI's are with FRG5). Looks good.


####5. Plotting####
library(ggplot2);library(dplyr);library(RColorBrewer)
FRG.colors=c("#d7191c","#fdae61","#ffffbf","#a6d96a","#1a9641")
#Subset the huge data frame for faster plotting
FRI.df.sub=FRI.df[sample(nrow(FRI.df), nrow(FRI.df)*0.01), ]#Subsample 1% of the df (18k cells)

#Plot FRI's
FRI.df.sub$trait=FRI.df.sub$fr #Change to trait of interest (also change axis label)
ggplot(FRI.df.sub)+
  geom_boxplot(aes(x=fri,y=trait,group=fri),notch=T)+ 
  geom_smooth(aes(x=fri,y=trait),method="lm",col="black")+
  scale_x_log10(breaks=c(5,15,25,35,50,100,200,500))+
  #scale_fill_manual(values = brewer.pal(n=8,name="PuOr"))+ #Need to set fill as factor(fri)
  #scale_fill_distiller(palette="PuOr")+
  annotation_logticks(sides="b")+
  labs(y="Fire-resistance index (0-1 scale)",x="Median fire return interval")+
  theme_bw()
dev.copy2pdf(file="./figures/ESA/Fire.resistance~MFRI.pdf") 

#Plot FRG's
FRI.df.sub$trait=FRI.df.sub$fr #Change to trait of interest (also change axis label)
ggplot(na.omit(FRI.df.sub))+
  geom_boxplot(aes(x=factor(frg),y=trait),notch=T)+ 
  #scale_fill_manual(values = brewer.pal(n=8,name="PuOr"))+ #Need to set fill as factor(fri)
  #scale_fill_distiller(palette="PuOr")+
  labs(y="Fire-resistance index (0-1 scale)",x="Fire Regime Group")+
  theme_bw()
dev.copy2pdf(file="./figures/ESA/Fire.resistance~FRG.pdf") 