

library(sp)
library(raster)
library(googlesheets)
library(dplyr)

# get trait data
#gs_auth(new_user=TRUE)
traits_full=gs_title("Traits")
master=gs_read(traits_full, ws="Master")
d <- master %>%
  filter(Western==1,
         California==1,
         Gymno==1,
         !Code %in% c("Tor_Cal","Tsu_Mer", "Tsu_Het"),
         Subsp<1,
         !is.na(CodeNum))

# get raster data
r <- stack(paste0("E:/fire_traits/FT_Analysis/LiveBasalAreaRasters/LiveBasalAreaRasters/",
                  "s", d$CodeNum,".img"))

# crop
ext <- extent(-126.947510, -113.796387, 32.768800, 42.000325)
ext <- as(ext, "SpatialPolygons")
crs(ext) <- CRS("+proj=longlat +ellps=WGS84 +datum=WGS84")
ext <- spTransform(ext, crs(r))
r <- crop(r, ext)

# total basal area
ba.tot <- sum(r)
ba.tot <- reclassify(ba.tot, c(750,Inf,NA,  -Inf,0,NA))

# specieswise fraction of basal area
ba.weighted <- reclassify(r / ba.tot, c(-Inf,0,NA))

# add MDS axes (generated in separeate script) to species table
mdsd <- as.data.frame(vare.mds$points)
mdsd$Code <- rownames(mdsd)
d <- left_join(d, mdsd)

# community trait means
traits <- dplyr::select(d, Bark_Thickness, Plant.height, MDS1, MDS2)
traits <- lapply(traits, function(x) sum(x * ba.weighted, na.rm=T) / sum(ba.weighted, na.rm=T))
traits <- stack(traits)
traits <- mask(traits, ba.tot)

# experiment with a plot mapping 2d geographic space to 2d ordination space
library(colormap) # get this using devtools::install_github("matthewkling/colormap")
library(ggplot2)
library(gridExtra)
library(grid)
f <- as.data.frame(rasterToPoints(traits))
colors <- colors2d(f[,c("MDS1", "MDS2")])

map <- ggplot(f, aes(x, y)) + 
  geom_raster(fill=colors) +
  ggmap::theme_nothing() +
  coord_fixed()

vars <- as.data.frame(vare.mds$species)
vars$var <- rownames(vars)
vars$MDS1 <- vars$MDS1 * max(f$MDS1) / max(vars$MDS1) * .8
vars$MDS2 <- vars$MDS2 * max(f$MDS2) / max(vars$MDS2) * .8
randos <- sample(nrow(f), 100000)
scatter <- ggplot(f[randos,], aes(MDS1, MDS2)) +
  geom_point(color=colors[randos]) +
  geom_segment(data=vars, aes(x=0,y=0,xend=MDS1,yend=MDS2), color="gray75") +
  geom_text(data=vars, aes(MDS1, MDS2, label=var), size=8) +
  theme_minimal() +
  coord_fixed(ratio=1) +
  labs(title="Community trait means\nin species trait space\n\n") +
  theme(text=element_text(size=20), plot.title=element_text(size=35))

plot <- arrangeGrob(scatter, map, ncol=2, widths=c(1,1.5))
png("E:/fire_traits/mds_map.png", width=1500, height=1000)
grid.draw(plot)
dev.off()




