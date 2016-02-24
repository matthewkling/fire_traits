#Jens Stevens stevensjt@gmail.com
#2/9/16
#Processing traits data

library(googlesheets); library(dplyr)
#Tutorial on googlesheets: http://htmlpreview.github.io/?https://raw.githubusercontent.com/jennybc/googlesheets/master/vignettes/basic-usage.html
gs_auth(new_user=TRUE) #Paste the link into a browser, and enter the auth code that the browser takes you to.
#my_sheets <- gs_ls()

####1. Pre-processing that only needs to be run once####
#1a
traits_full=gs_title("Traits")
master=gs_read(traits_full, ws="Master")
#sppCodes=gs_read(traits_full, ws="sppCodes")
#TRY=gs_read(traits_full, ws="TRY")

#Generate 6-letter species codes:
#genus=substr(master$Scientific_Name,1,3)
#species.pos=apply(master[,"Scientific_Name"],MARGIN=2,FUN=regexpr,pattern="_")[,1]
#species=substr(master$Scientific_Name,species.pos+1,species.pos+3)
#species=paste(toupper(substring(species, 1,1)), substring(species, 2),sep="" )
#master$Code=paste(genus,species,sep="_")
#master$Generic[grep("Spp",master$Code)]=1
#master$Generic[is.na(master$Generic)]=0
#master$California[which(master$Western==0)]=0

#Add in TRY data:
#TRY=TRY[,-1] #Remove Scientific_Name column so it is not duplicated.
#names(TRY)=c("Code","Leaf.area","SLA","Leaf.CN","Leaf.length","Leaf.lifespan","Leaf.N","Leaf.NP","Leaf.thickness","Leaf.width","Plant.height","Plant.lifespan","Seed.dry.mass","Stem.diameter","SSD")
#combined=merge(master,TRY,by="Code",all=T) #Don't discard species that aren't in both db's

#Add in sppCode Number data
#sppCodes=sppCodes[,-c(1)] #Remove common.name column so it is not duplicated.
#names(sppCodes)=c("Scientific_Name2","Code","CodeNum")
#combined=merge(master,sppCodes,by="Code",all=T) #Don't discard species that aren't in both db's

#1b

#When the spreadsheet edits have been made, upload the working "traits" file back into the appropriate tab of the "traits_full" spreadsheet
traits_full=gs_edit_cells(traits_full,ws="Combined",input=combined)
#traits_full=gs_edit_cells(traits_full,ws="TRY",input=TRY)

####2. An initial ordination analysis####
library(vegan)

#2a: Set up data for ordinations
master=gs_read(traits_full, ws="Master") #Read in master data from spreadsheet
#Subset the data to western conifers
subs=which(master$Western==1 & master$Gymno==1)
d=master[subs,]
#Remove some outliers:
d=d[-pmatch(c("Tor_Cal","Tsu_Mer", "Tsu_Het"),d$Code),]
#Tor_Cal has huge seed mass; Tsu_Het and Tsu_Mer really thick bark
#Remove the subspecies:
d=d[which(d$Subsp<1),]
d.ord=d[,c("Bark_Thickness","Leaf_Longevity","Plant.height","Plant.lifespan","Seed.dry.mass","Serotiny.cat")]
cor(d.ord,use="complete.obs") #Check out correlations
#Log transform the variables that need it:
d.ord$Plant.height=log(d.ord$Plant.height)
d.ord$Plant.lifespan=log(d.ord$Plant.lifespan)
d.ord$Seed.dry.mass=log(d.ord$Seed.dry.mass)
cor(d.ord,use="complete.obs") #Check out correlations
#Correlations are much better; thick bark species also have shorter leaf longevity, taller height, longer lifespan, larger seeds, and are less likely to be serotinous.
d.ord.complete=d.ord[which(complete.cases(d.ord)),]
rownames(d.ord.complete)=d[which(complete.cases(d.ord)),"Code"]$Code

#2b: Run data for ordinations
vare.dis <- vegdist(d.ord.complete)
vare.mds0 <- monoMDS(vare.dis)
ordiplot(vare.mds0, type = "t") #Plot 1
vare.mds <- metaMDS(d.ord.complete, trace = FALSE)
plot(vare.mds, type = "t")


