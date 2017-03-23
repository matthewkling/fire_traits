
library(tidyverse)
library(raster)
library(rgdal)
library(PerformanceAnalytics)
library(RColorBrewer)
library(viridis)
library(broom)

setwd("E:/fire_traits/fire_traits")

### 1. Trait data prep

# Calculate mean flammability traits for PICO because working on the species-level
flam <- read.csv("./data/flam_traits.csv", stringsAsFactors = F)
flam[nrow(flam)+1,"Scientific_Name"] <- "Pinus_contorta"
flam[nrow(flam),c(2,4,6,8)] <- 
      colMeans (flam[grep("contorta",flam$Scientific_Name),c(2,4,6,8)],na.rm=TRUE)
flam[nrow(flam),c(3,5,7,9)] <- "Derived from Banwell & Varner"

# Merge component datasets
md <- Reduce(function(x, y) 
      merge(x, y, all=TRUE), 
      list(read_csv("./data/species_list.csv"), 
           read_csv("./data/ffe_traits.csv"), 
           read_csv("./data/try_traits.csv"),
           read_csv("./data/S&A_traits.csv"),
           flam
      ) 
)

# Choose variables of interest. "25.4" refers to dbh of tree in cm.
vars_of_interest = c("Scientific_Name","Code","CodeNum","California","Western","Gymno","Has_BA",
                     "Bark.Thickness.25.4","Plant.height","Self.pruning","Flame_ht","Flame_duration")
d <- dplyr::select(md,which(names(md)%in%vars_of_interest))


# Choose species of interest
# Here, California gymnosperm species that have basal area layers
d <- filter(d, California==1, Gymno==1, Has_BA==1)
#d <- filter(d, Western==1, Gymno==1, Has_BA==1)

spp <- d

### 2. Calculate fire-resistance score for species of interest (fast)####

#Extract the quantile of each trait for each species 
traits_of_interest <- 
      c("Bark.Thickness.25.4","Plant.height","Self.pruning",
        "Flame_ht","Flame_duration")
d$bt.quant=ecdf(d$Bark.Thickness.25.4)(d$Bark.Thickness.25.4)
d$ph.quant=ecdf(d$Plant.height)(d$Plant.height)
d$sp.quant=ecdf(d$Self.pruning)(d$Self.pruning)
d$fh.quant=ecdf(d$Flame_ht)(d$Flame_ht)
d$fd.quant=ecdf(-d$Flame_duration)(-d$Flame_duration)

#Apply (across each row) the weighted mean of the traits of interest, weighted by its completeness
quants_of_interest <- c("bt.quant","ph.quant","sp.quant","fh.quant","fd.quant")
wts <- #Weight each quantile variable by its completeness
      d[,quants_of_interest] %>% 
      apply(MARGIN=2,function(x) length(which(!is.na(x)))/length(x))
d$frs <-
      d[,quants_of_interest] %>%
      apply(MARGIN=1, function(x) weighted.mean(x=x, w=wts, na.rm= TRUE) )

#write_csv(d[,c(1,8:ncol(d))],"./manuscript/tables/TableS1.csv")


# Trait correlations
d$log_Bark.thickness=log(d$Bark.Thickness.25.4)
traits_of_interest <- 
      c("log_Bark.thickness","Plant.height","Self.pruning",
        "Flame_ht","Flame_duration")
#chart.Correlation(d[,traits_of_interest])
#chart.Correlation(d[,quants_of_interest])
#dev.copy2pdf(file="./figures/MS1/FigS1_trait_correlations.pdf") 





##################
##################


# LOAD LANDFIRE DATA

frs <- raster("../large_files/fr.weighted.tif")
fri <- raster("./GIS/FRG_Rasters/Western_FRI_Clip_250m.tif")
frg <- raster("./GIS/FRG_Rasters/Western_FRG_Clip_250m.tif")

frs <- crop(frs, crop(fri, frs))
fri <- crop(fri, frs)
frg <- crop(frg, frs)


# LOAD USFS DATA

ba.weights.stack <- stack("../large_files/ba.weights.stack.tif")
ba.rasters.study <- stack("../large_files/ba.rasters.study.tif")
ba.rasters.other <- stack("../large_files/ba.rasters.other.tif")
ba.tot.study <- raster("../large_files/ba.tot.study.tif")

ba.tot.other <- calc(ba.rasters.other, sum)
ba.prop <- ba.tot.study / (ba.tot.study + ba.tot.other)


# reclassify FRI

rclmat2 <- data.frame(from=c(-Inf,0,2,4,6,8,11,17,20,22),
                      to=c(0,2,4,6,8,11,17,20,22,150),
                      becomes=c(NA,5,15,25,35,50,100,200,500,NA))
fri <- reclassify(fri, as.matrix(rclmat2), right=T)



stop("skip this slow, underreveloped chunk of code")
# basal area rasters
ba <- stack(list.files("E:/fire_traits/fire_traits/input_data/LiveBasalAreaRasters/", full.names=T))

# crop
rs <- ba[[paste0("s", spp$CodeNum)]]
rs <- crop(rs, frs)


# total basal area, all species
#ba.tot <- sum(r)
#ba.tot <- reclassify(ba.tot, c(750,Inf,NA,  -Inf,0,NA))

# total basal area, focal cali gymno species
#ba.cg <- sum(rs)

# specieswise fraction of total basal area
ba.weighted <- reclassify(rs / ba.tot.study, c(-Inf,0,NA))

# community trait means
traits <- spp[,traits_of_interest]
traits <- lapply(traits, function(x) sum(x * ba.weighted, na.rm=T) / sum(ba.weighted, na.rm=T))
traits <- stack(traits)
#####traits <- mask(traits, ba.tot)







#### VISUALIZE

# spatial sync -- shortcut assumes that raster misalignment is negligible!
fri_tmp <- fri
frg_tmp <- fri
fri <- frs; values(fri) <- values(fri_tmp)
frg <- frs; values(frg) <- values(frg_tmp)

s <- stack(frs, fri, frg, ba.prop, ba.tot.other + ba.tot.study)
d <- as.data.frame(rasterToPoints(s))
names(d) <- c("x", "y", "frs", "fri", "frg", "ba_prop", "ba_tot")

d <- filter(d, !is.na(frs), !is.na(fri))


# a regression residuals approach to flagging mismatches
#d$fri[d$fri==0] <- 1
fit <- lm(frs~log10(fri), data=d)
d$mismatch <- fit$residuals
fit2 <- lm(log10(fri)~frs, data=d)
d$mismatch2 <- fit2$residuals

#ss <- sample(nrow(d), 10000)
#cor(d$frs[ss], log(d$fri)[ss], method="kendall", use="complete.obs")



# cartography prep

usa <- getData("GADM", country="USA", level=1) %>%
      spTransform(crs(s)) %>%
      fortify() %>%
      mutate(group=paste("usa", group))
can <- getData("GADM", country="CAN", level=1) %>%
      spTransform(crs(s)) %>%
      fortify() %>%
      mutate(group=paste("can", group))
mex <- getData("GADM", country="MEX", level=1) %>%
      spTransform(crs(s)) %>%
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



# map of FRS
p <- ggplot() + 
      geom_polygon(data=borders, aes(long, lat, group=group),
                   fill="gray95", color=NA) +
      geom_raster(data=d,
                  aes(x, y, fill=frs)) +
      geom_path(data=borders, aes(long, lat, group=group),
                size=.25) +
      scale_fill_viridis() +
      minimalism +
      view +
      guides(fill=guide_colourbar(barwidth=15)) +
      labs(fill="score (0-1)",
           title="Fire resistance index\nweighted by species abundance")
ggsave("temp/frs.png", p, width=7, height=9, units="in")

# map of FRI
p <- ggplot() + 
      geom_polygon(data=borders, aes(long, lat, group=group),
                   fill="gray95", color=NA) +
      geom_raster(data=d,
                  aes(x, y, fill=fri)) +
      geom_path(data=borders, aes(long, lat, group=group),
                size=.25) +
      scale_fill_viridis(option="B", trans="log10", breaks=c(1,3,10,30,100,300,1000)) +
      minimalism +
      view +
      guides(fill=guide_colourbar(barwidth=15)) +
      labs(fill="years\n",
           title="\nFire return interval")
ggsave("temp/fri.png", p, width=7, height=9, units="in")

# map of mismatch
p <- ggplot() + 
      geom_polygon(data=borders, aes(long, lat, group=group),
                   fill="gray95", color=NA) +
      geom_raster(data=d,
                  aes(x, y, fill=mismatch)) +
      geom_path(data=borders, aes(long, lat, group=group),
                size=.25) +
      scale_fill_gradientn(colours=c("darkblue", "gray", "gray", "darkred"),
                           limits=max(abs(d$mismatch))*c(-1,1)) +
      minimalism +
      view +
      labs(fill="FRS~FRI residual\n",
           title="\nFRS-FRI mismatch") +
      guides(fill=guide_colourbar(barwidth=15))
ggsave("temp/mismatch.png", p, width=7, height=9, units="in")

# map of ba_prop
p <- ggplot() + 
      geom_polygon(data=borders, aes(long, lat, group=group),
                   fill="gray95", color=NA) +
      geom_raster(data=d,
                  aes(x, y, fill=ba_prop)) +
      geom_path(data=borders, aes(long, lat, group=group),
                size=.25) +
      scale_fill_viridis(limits=c(.5,1)) +
      minimalism +
      view +
      labs(fill="proportion basal area represented",
           title="\nData completeness") +
      guides(fill=guide_colourbar(barwidth=15))
ggsave("temp/completeness.png", p, width=7, height=9, units="in")


# scatterplot frs ~ fri, color=residual
p <- ggplot(sample_n(d, 1000000), aes(fri, frs, color=mismatch)) +
      geom_point() +
      geom_smooth(method=lm, se=F, color="black") +
      scale_colour_gradientn(colours=c("darkblue", "gray", "gray", "darkred"),
                             limits=max(abs(d$mismatch))*c(-1,1)) +
      scale_x_log10(breaks=c(1,3,10,30,100,300,1000)) +
      theme_minimal() +
      labs(title="FRS vs FRI\n",
           fill="residual") +
      guides(fill=guide_colourbar(barwidth=15))
ggsave("temp/frs_fri_mismatch.png", p, width=8, height=8, units="in")

# scatterplot: is mismatch due to incompleteness? (apparently not)
p <- ggplot(data=sample_n(d, 100000), aes(ba_prop, abs(mismatch))) +
      geom_point(alpha=.25) +
      geom_smooth(method=lm, se=F, size=2) +
      scale_colour_gradientn(colours=c("darkblue", "gray", "darkred"),
                             limits=max(abs(d$mismatch))*c(-1,1)) +
      theme_minimal() +
      labs(x="proportion basal area represented")
ggsave("temp/mismatch_completeness.png", p, width=8, height=8, units="in")

# scatterplot fri ~ frs, color=residual
p <- ggplot(sample_n(d, 1000000), aes(frs, fri), color="gray") +
      geom_point() +
      geom_smooth(method=lm, se=F, color="black") +
      scale_colour_gradientn(colours=c("darkblue", "gray", "gray", "darkred"),
                             limits=max(abs(d$mismatch))*c(-1,1)) +
      scale_y_log10(breaks=c(1,3,10,30,100,300,1000)) +
      theme_minimal() +
      labs(title="FRS vs FRI\n",
           fill="residual") +
      guides(fill=guide_colourbar(barwidth=15))
ggsave("temp/fri_frs_scatter.png", p, width=8, height=8, units="in")


