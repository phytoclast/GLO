#############################
#acquire topographic covariates
#############################
library(sf)
library(terra)
library(dplyr)
library(ggplot2)
library(ranger)
library(gam)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

vars90 <- rast('gis/vars90.tif')
#supplement MNFI openings and extract variables
mnfi <- read_sf('C:/GIS/Vegetation/Veg1800/MichLP_1800veg.shp')
mnfi <- st_transform(mnfi, crs(vars90))
mnfi <-  mnfi |> subset(select=c('VEGCODE', 'LEGCODE', 'COVERTYPE'))
mlra <- read_sf('gis/MLRA_2017.shp')
mlra <- st_transform(mlra, crs(vars90))
mlra <-  mlra |> subset(select=c('MLRA', 'LRU'))
mlra <-  mlra |> mutate(LRU2 = case_when(grepl('AA', LRU) ~ paste0(MLRA,'A'),
                                        grepl('AB', LRU) ~ paste0(MLRA,'B'),
                                        grepl('A', LRU) & !grepl('A', MLRA)  ~ paste0(MLRA,'A'),
                                        grepl('B', LRU) & !grepl('B', MLRA)  ~ paste0(MLRA,'B'),
                                        TRUE ~ paste0(MLRA)))
# openarea <- st_area(mnfi) |> sum() |> as.numeric()
# nopenpoints <- 3*openarea/1600^2
# mnfi.vect <- vect(mnfi)
# openpoints <- spatSample(mnfi.vect, nopenpoints)
# openpoints <- st_as_sf(openpoints)
# openpoints <- openpoints[,'geometry']
# openpoints <- openpoints |> mutate(Level2 = 'open', Species = 'open', BA = 0, DT= 0)
# pts <- pts |> bind_rows(openpoints)
# 
# pts.vars90 <- pts |> vect() |> project(vars90)
# pts.vars90 <- extract(vars90,pts.vars90)
# pts.vars90 <- pts |> cbind(pts.vars90)
# saveRDS(pts.vars90, 'pts.vars90.RDS')
# pts.vars270 <- pts |> vect() |> project(vars270)
# pts.vars270 <- extract(vars270,pts.vars270)
# pts.vars270 <- pts |> cbind(pts.vars270)
# saveRDS(pts.vars270, 'pts.vars270.RDS')

# pca90 <- princomp(vars270, cor = T, maxcell=ncell(vars270)*0.1)
# pca90grid <- predict(vars270, pca90)
# writeRaster(pca90grid,'gis/pca90grid.tif', overwrite=T)

pts.vars90 <- readRDS('pts.vars90.RDS') |> st_transform(crs(vars90))
pts.vars90 <- pts.vars90 |> mutate(Level2 = ifelse(Species %in% 'Thuja', Species, Level2))
#extract a list of taxon names
ttt <- pts.vars90 |> group_by(Level2) |> summarize(npts = length(Level2))
ttt <- subset(ttt, !is.na(Level2) & !Level2 %in% c('unk','no tree','other hardwood','no data','rock','open','dead') & npts >= 10)
ttt <- ttt$Level2
sporder <- c("Pinus","Tsuga","Fagus","Acer","Tilia","Carya","Quercus","Ulmus","Fraxinus","Larix","Thuja","Picea","Abies","PopLirio", "Betula",
             "Aesculus","Alnus","Castanea","Celtis","Cercis","Cornus","Crataegus","Cupressaceae","Juglans","Morus","NyssaLiquidambar","OstryaCarpinus","Platanus","Prunus","RobiniaGleditsia","Salix","Sassafras")
ptsMLRA <- extract(vect(mlra), vect(pts.vars90))
ptsMNFI <- extract(vect(mnfi), vect(pts.vars90)) 
ptsMNFI <- ptsMNFI |> group_by(id.y) |> summarise(VEGCODE=first(VEGCODE), LEGCODE=first(LEGCODE), COVERTYPE=first(COVERTYPE)) |> ungroup()


pts <- pts.vars90 |> cbind(ptsMLRA[,!colnames(ptsMLRA) %in% 'id.y'])
pts <- pts |> cbind(ptsMNFI[,!colnames(ptsMNFI) %in% 'id.y'])

pts <- pts |> subset(Level2 %in% ttt & !is.na(MLRA) & !is.na(COVERTYPE))
pts <- pts |> mutate(Level2 = factor(Level2, levels = sporder))
pts <- pts |> subset(!is.na(watertable) & !is.na(sand150) & !is.na(slope) & !is.na(OM150)) |> 
  mutate(drain = case_when(watertable <= 25 ~ 'wet',
                           watertable <= 100 ~ 'moist',
                           TRUE ~ 'dry'),
         text = case_when((sand150 >= 80 & sand50 > 70) | sand50 > 80 ~  'sandy',
                          histic >= 0.5 ~ 'mucky',
                          (sand150 < 40 | sand50 < 35) | sand50 < 40 & clay150 > 30 ~  'clayey',
                          
                          TRUE ~ 'loamy'),
         chem = case_when(Bhs >= 0.8 & spodic >= 0.8 & drain %in% 'dry' & text %in% 'sandy' ~ 'super spodic',
                          spodic >= 0.8 & drain %in% 'dry' & text %in% 'sandy' ~ 'spodic',
                          (pH50 >= 6 | carbdepth < 100) & drain %in% 'dry' & text %in% 'sandy' ~ 'rich',
                          (pH50 <= 5.5 | spodic >= 0.8) & text %in% c('sandy','mucky') &
                            !drain %in% 'dry' ~ 'dysic',
                           text %in% c('sandy','mucky') & !drain %in% 'dry' ~ 'euic',
                          TRUE ~ 'typic'),
         form = case_when(floodfrq > 0 ~ 'floodplain',
                          slope >= tan(15/100) ~ 'slopes',
                          TRUE ~ 'plains'))

                          
pts <- pts |> mutate(site = (paste(LRU2, drain, chem, text, form)))

ptssum <- pts |> st_drop_geometry() |> summarise(.by=c(site, COVERTYPE,LRU2, drain, chem, text, form), ntrees = length(Level2)) |> mutate(.by=c(site),prop = round(100*ntrees/sum(ntrees),1)) 

ptssum2 <- pts |> st_drop_geometry() |> summarise(.by=c(site, COVERTYPE, Level2), abund = length(Level2)) |> mutate(.by=c(site,COVERTYPE),abund = round(100*abund/sum(abund),1)) 

ptssum2 <- ptssum |> left_join(ptssum2) |> subset(ntrees >= 100)
ptssum <- ptssum |> subset(ntrees >= 100)
ptpivot <- ptssum2 |> tidyr::pivot_wider(names_from=Level2, values_from = abund, names_sort = TRUE)
cortab0 <- subset(ptpivot, select=c(sporder)) 
freplace <- function(x){
  x <- ifelse(is.na(x),0,x)
  return(x)
}
#correlation
cortab0 <- sapply(cortab0, FUN=freplace) 
cortab <- cortab0 |> cor() |> as.data.frame()
ptpivotx <- ptpivot |> mutate(across(all_of(sporder), freplace) )
#clustering
# sptab <- cortab |> mutate(across(unique(ptssum2$Level2), function(x){1-(x+1)/2}) )
# d <-vegan::vegdist(sptab,  method="bray")
# c1 <- pam(d,5, nstart=20)
# sptab <- as.data.frame(sptab) |> mutate(clust = c1$clustering)
# compsum <- sptab |> group_by(clust) |> summarise(across(unique(ptssum2$Level2), mean))
# compsum <- compsum |>  mutate(across(unique(ptssum2$Level2), function(x){(1-x)}))
d <-vegan::vegdist(cortab0,  method="bray")
c1 <- pam(d,5, nstart=20)
ptpivotx <- ptpivotx |> mutate(clust = c1$clustering)

compsum <- ptpivotx |> group_by(clust) |> summarise(across(unique(all_of(sporder)), mean))

write.csv(ptpivotx,'presettlementvegbysoil2.csv', row.names = F, na='')
