#Jens Stevens stevensjt@gmail.com
#2/11/16
#Spatial analyses

library(sp)
library(raster)

Centroids=data.frame("Longitude"=c(-120.0423938,-119.9562198),"Latitude"=c(39.0890715,38.8980169))
coordinates(Centroids)=c("Longitude","Latitude")
sp::proj4string(Centroids)= CRS("+proj=longlat +ellps=WGS84 +datum=WGS84")

#Read in data for species in question:
master=gs_read(traits_full, ws="Master") #Read in master data from spreadsheet
#Subset the data to western conifers
subs=which(master$Western==1 & master$Gymno==1)
d=master[subs,]
#Remove some outliers:
d=d[-pmatch(c("Tor_Cal","Tsu_Mer", "Tsu_Het"),d$Code),]
#Tor_Cal has huge seed mass; Tsu_Het and Tsu_Mer really thick bark
#Remove the subspecies:
d=d[which(d$Subsp<1),]
d=d[!is.na(d$CodeNum) & d$California==1,]
#Just ABCO, PICO and PIJE:
d=d[pmatch(c("Abi_Con","Pin_Con", "Pin_Jef"),d$Code),]


#Get the species codes for the rasters you want:
sppCodeNums=as.character(d$CodeNum)
sppFileNames=paste0("s",d$CodeNum,".img")
#rasters_to_get=list.files("/Users/Jens/Documents/Davis/Post-Doc/FireTraits/GIS/LiveBasalAreaRasters/")
raster.list=list()
for(r in sppFileNames){
  raster.list[[r]]=
    raster(paste0("/Users/Jens/Documents/Davis/Post-Doc/FireTraits/GIS/LiveBasalAreaRasters/",r))
}

raster.sub=list()
Centroids.transformed=spTransform(Centroids,raster::crs(raster.list[[1]]))

tahoe=as(extent(-2065120.76, -2009333.45, 2027193.20, 2064509.66), 'SpatialPolygons')
crs(tahoe)=crs(raster.list[[r]])
for(r in 1:length(sppFileNames)){
  raster.sub[[r]] <- crop(raster.list[[sppFileNames[r]]], tahoe)
  plot(raster.sub[[r]],main=d$Code[r])
}

for(r in sppFileNames){ #Deprecated for now
  extract.sub[[r]]=
    raster::extract(raster.list[[r]], Centroids.transformed, buffer=10000)
  raster.sub <- raster(r)
  for (i in 1:length(z)) { s[z[[i]]] <- i }
  
} #Deprecated for now
str(rasters.sub)
plot(rasters.sub[[1]])
#http://stackoverflow.com/questions/15824853/large-img-file-processing-in-r-gis

#Stacking and weighting
s=stack(raster.sub)
r.o=overlay(s,fun=function(x,y,z){return(x+y+z)})

raster.weighted=thickness=thickness.weighted=list()
for(r in 1:length(sppFileNames)){
  raster.weighted[[r]]=overlay(raster.sub[[r]],r.o,fun=function(x,y){return(x/y)})
  #plot(raster.weighted[[r]],main=d$Code[r])
  thickness[[r]]=raster.weighted[[r]]
  thickness[[r]]@data@values=d$Bark_Thickness[r]
  
  
  #thickness.weighted[[r]]=
  #  calc(raster.weighted[[r]],fun=function(x){return(x*d$Bark_Thickness[r])})
  #plot(thickness.weighted[[r]],main=d$Code[r],legend.args=list(text="Weighted Bark Thickness"))
  #thickness.weighted[[r]]@data@values=
    #weighted.mean(d$Bark_Thickness[r])
}

sw=stack(raster.weighted)
st=stack(thickness)
r.w=calc(s2,fun=function(x,y,z){return(weighted.mean(d$Bark_Thickness,w=c(x,y,z)))})
#weighted.mean(d$Bark_Thickness,w=c(0.1,0.1,0.9))
thickness.weighted=raster::weighted.mean(st,sw)
plot(thickness.weighted,main="Bark thickness (cm) \nweighted by species abundance")
