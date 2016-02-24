#Jens Stevens stevensjt@gmail.com
#2/21/16
#Spatial analyses

library(sp)
library(raster)
library(googlesheets)

####1. Import and process data####
#gs_auth(new_user=TRUE) #Paste the link into a browser, and enter the auth code that the browser takes you to.
traits_full=gs_title("Traits")
master=gs_read(traits_full, ws="Master") #Read in master data from spreadsheet
#Subset the data to western conifers
subs=which(master$Western==1 & master$Gymno==1)
d=master[subs,]
#Remove some outliers:
d=d[-pmatch(c("Tor_Cal","Tsu_Mer", "Tsu_Het"),d$Code),]
#Tor_Cal has huge seed mass; Tsu_Het and Tsu_Mer really thick bark
#Remove the subspecies:
d=d[which(d$Subsp<1),]
#Select California species 
d=d[d$California==1,]
#Select species that have corresponding basal area layers ("CodeNum")
d=d[!is.na(d$CodeNum),]
#Just ABCO, PICO and PIJE:
#d=d[pmatch(c("Abi_Con","Pin_Con", "Pin_Jef"),d$Code),]


#Get the species codes for the rasters you want:
sppFileNames=paste0("s",d$CodeNum,".img")
raster.list=list()
for(r in sppFileNames){
  raster.list[[r]]=raster(paste0("E:/fire_traits/FT_Analysis/LiveBasalAreaRasters/LiveBasalAreaRasters/",r))
}

#Visualize a raster
plot(raster.list[[1]])

####2. Process raster data####
#2a) Define extent. One way to do that is here: www.latlong.net
MapExtent=data.frame("Longitude"=c(-126.947510,-113.796387),"Latitude"=c(42.000325,32.768800))#California
#MapExtent=data.frame("Longitude"=c(-120.0423938,-119.9562198),"Latitude"=c(39.0890715,38.8980169))#Tahoe
coordinates(MapExtent)=c("Longitude","Latitude") 
sp::proj4string(MapExtent)= CRS("+proj=longlat +ellps=WGS84 +datum=WGS84") 
MapExtent=spTransform(MapExtent,raster::crs(raster.list[[1]])) #Reproject extent to be same as raster layers.

SP.Extent=as(extent(MapExtent),'SpatialPolygons')#Reclassify the extent as a SpatialPolygons layer for cropping
crs(SP.Extent)=crs(raster.list[[1]]) #Assign coordinate system
raster.sub=list() #Set up empty list for basal area rasters to be subset.

for(r in 1:length(sppFileNames)){#For each species BA raster, crop to the area of interest
  raster.sub[[r]] <- crop(raster.list[[sppFileNames[r]]], SP.Extent)
  #plot(raster.sub[[r]],main=d$Code[r]) #Optional: Visualize each basal area layer as the code runs
} #Takes ~2 minutes
#Optional, explore raster values
#vals=getValues(raster.sub[[1]])
#hist(vals[vals>1]); max(vals)

#for(r in sppFileNames){ #Deprecated for now
#  extract.sub[[r]]=
#    raster::extract(raster.list[[r]], Centroids.transformed, buffer=10000)
#  raster.sub <- raster(r)
#  for (i in 1:length(z)) { s[z[[i]]] <- i }
#  #http://stackoverflow.com/questions/15824853/large-img-file-processing-in-r-gis
#} #Deprecated for now

#2b) Stacking and weighting the different species rasters
s=stack(raster.sub) #Create a raster stack of the different species
ba.tot=overlay(s,fun=sum) #Create layer for cumulative basal area. Takes ~4 minutes. 
#hist(vals[vals>720]) #There are a few huge outliers here in the redwood region, may want to truncate.
#Adjust values for total basal area
ba.tot[ba.tot>750]=NA #Remove outliers with large BA
ba.tot[ba.tot==0]=NA #If there is no basal area from any of the focal species, set basal area to NA.

ba.weighted=thickness=thickness.weighted=list() #Create empty lists
for(r in 1:length(sppFileNames)){ #For each species, calculate its fraction of stand basal area for each pixel, and create a new raster containing its trait of interest (bark thickness for now). Takes ~3 minutes
  ba.weighted[[r]]=overlay(raster.sub[[r]],ba.tot,fun=function(x,y){return(x/y)}) #Calculate the fraction of stand basal area comprised of the given species
  ba.weighted[[r]][is.na(ba.weighted[[r]])]=0 #Any cells where the species is absent, basal area weights should be 0 rather than NA. The NA's come from dividing by NA in the overlay() function above, but the weights where trees are absent should be 0. 
  #ba.weighted[[r]]@data@values[is.na(ba.weighted[[r]]@data@values)]=0 #The old/incorrect way to do the above. CHECKME to confirm.
  #plot(ba.weighted[[r]],main=d$Code[r])
  thickness[[r]]=ba.weighted[[r]] #Create a new raster for bark thickness. The initial values will be overwritten.
  thickness[[r]][]=d$Bark_Thickness[r] #Assign every pixel the given species' bark thickness.

  
  #thickness.weighted[[r]]=
  #  calc(raster.weighted[[r]],fun=function(x){return(x*d$Bark_Thickness[r])})
  #plot(thickness.weighted[[r]],main=d$Code[r],legend.args=list(text="Weighted Bark Thickness"))
  #thickness.weighted[[r]]@data@values=
    #weighted.mean(d$Bark_Thickness[r])
}

sw=stack(ba.weighted)
st=stack(thickness)
#sw=stack(ba.weighted[[2]],ba.weighted[[20]]) #Option to look at fewer species
#st=stack(thickness[[2]],thickness[[20]]) #Option to look at fewer species
#r.w=calc(s2,fun=function(x,y,z){return(weighted.mean(d$Bark_Thickness,w=c(x,y,z)))})
#weighted.mean(d$Bark_Thickness,w=c(0.1,0.1,0.9))
thickness.weighted=raster::weighted.mean(x=st,w=sw,na.rm=T) #Calculate the mean bark thickness of each pixel, weighted by the relative abundance of that species in the stand. Takes ~5 minutes
#vals=getValues(thickness.weighted)
plot(thickness.weighted,main="Bark thickness (cm) \nweighted by species abundance")
