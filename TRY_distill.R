
##### distill raw TRY dataset into usable form covering relevant species


# list of tree species in FFE dataset
spp <- read.csv("E:/fire_traits/FFE/traits_clean.csv", stringsAsFactors=F)

# read in TRY data, in chunks due to memory limits, retaining only relevant species
n <- 1e5
i <- 0
status <- "incomplete"
while(status=="incomplete"){
      if(i<1){
            b <- read.delim("E:/fire_traits/TRY/1645_08022016124158/1645.txt", nrows=n)
            vars <- row.names(b)
      }else{b <- read.delim("E:/fire_traits/TRY/1645_08022016124158/1645.txt",
                              skip=i*n+1, nrows=n, row.names=vars)}
      if(nrow(b)<n) status <- "complete"
      b <- b[b$SpeciesName %in% spp$scientific_name,]
      if(i<1) d <- b else(d <- rbind(d, b))
      i <- i + 1
}




