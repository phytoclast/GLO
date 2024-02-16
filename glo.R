setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
library(terra)
library(sf)
library(dplyr)

# dia=20
# dis=9
# (dia/200)^2/(dis^2)*10000*4
# dens <- rast('data/msb-paleon.26.1/PLS_Density_Level1_v1.0.nc')
# basl <- rast('data/msb-paleon.26.1/PLS_BasalArea_Level1_v1.0.nc')
# dens <- dens$Total
# dens <- ifel(dens > 10000,NA,dens)
# basl <- basl$Total
# basl <- ifel(basl > 1000000,NA,basl)
# crs(basl) <- 'EPSG:3175'
# crs(dens) <- 'EPSG:3175'
# writeRaster(basl, 'gis/basl.tif', overwrite=TRUE)
# writeRaster(dens, 'gis/dens.tif', overwrite=TRUE)


rackCorners <- function(pts, dataset, treenames, diams, dists, x, y){
  pts$id <- rownames(pts)
  ntrees = length(treenames)
 
  for(i in 1:ntrees){#i=1
    if(is.na(diams[1])){diam = NA}else{
      diam = pts[,diams[i]]}
    if(is.na(dists[1])){distz = NA}else{
      distz = pts[,dists[i]]}
    pts0 = data.frame(
      dataset = dataset,
      id = NA,
      treeseq = i,
      tname = pts[,treenames[i]],
      diam = diam,
      dists = distz,
      x=pts[,x],
      y=pts[,y]
    )
    if(i==1){pts1=pts0}else{pts1=rbind(pts1,pts0)}
  }
  pts1 <- mutate(pts1, id = paste(x,y))
  return(pts1)
}

pts <- read.csv('data/msb-paleon.32.0/PLS_Wisconsin_trees_Level0_v1.0.csv')
treenames <- c("SP1","SP2","SP3","SP4")
diams <- c("DIAM1","DIAM2","DIAM3","DIAM4")
dists <- c("DIST1","DIST2","DIST3","DIST4")
dataset='WI'
x = 'x_alb'
y = 'y_alb'
pts.WI <- rackCorners(pts, dataset, treenames, diams, dists, x, y)

pts <- read.csv('data/msb-paleon.27.0/PLS_Indiana_trees_Level0_v1.0.csv')
treenames <- c("L1_tree1","L1_tree2","L1_tree3","L1_tree4")
diams <- c("diameter","diameter2","diameter3","diameter4")
dists <- c("chainstree","chainstree2","chainstree3","chainstree4")
dataset='IN'
x = 'x'
y = 'y'
pts.IN <-  rackCorners(pts, dataset, treenames, diams, dists, x, y)

pts <- read.csv('data/msb-paleon.28.0/PLS_Illinois_trees_Level0_v1.0.csv')
treenames <- c("L1_tree1","L1_tree2","L1_tree3","L1_tree4")
diams <- c("diameter","diameter2","diameter3","diameter4")
dists <- c("chainstree","chainstree2","chainstree3","chainstree4")
dataset='IL'
x = 'x'
y = 'y'
pts.IL <-  rackCorners(pts, dataset, treenames, diams, dists, x, y)

pts <- read.csv('data/msb-paleon.29.0/PLS_SoutheasternMichigan_trees_Level0_v1.0.csv')
treenames <- c("L1_tree1","L1_tree2","L1_tree3","L1_tree4")
diams <- c("diameter","diameter2","diameter3","diameter4")
dists <- c("chainstree","chainstree2","chainstree3","chainstree4")
dataset='SEMI'
x = 'x'
y = 'y'
pts.SEMI <-  rackCorners(pts, dataset, treenames, diams, dists, x, y)

pts <- read.csv('data/msb-paleon.30.0/PLS_southernMichigan_trees_Level0_v1.0.csv')
treenames <- c("species1","species2","species3","species4")
diams <- c("diam1","diam2","diam3","diam4")
dists <- c("dist1","dist2","dist3","dist4")
dataset='SMI'
x = 'point_x'
y = 'point_y'
pts.SMI <-  rackCorners(pts, dataset, treenames, diams, dists, x, y)

pts <- read.csv('data/msb-paleon.31.0/PLS_northernMichigan_trees_Level0_v1.0.csv')
treenames <- c("SPP1","SPP2","SPP3","SPP4")
diams <- c("DBH1","DBH2","DBH3","DBH4")
dists <- c("DIST1","DIST2","DIST3","DIST4")
dataset='NMI'
x = 'x_alb'
y = 'y_alb'
pts.NMI <-  rackCorners(pts, dataset, treenames, diams, dists, x, y)

#pts <- read.csv('data/msb-paleon.2.0/SetTreeComp_West_Level1_v1.0.csv')

pts <- rbind(pts.WI,pts.IL,pts.IN,pts.NMI,pts.SMI,pts.SEMI)
pts.geo <- st_as_sf(pts, coords=c(x='x', y='y'), crs=3175)

pts.sf <- st_read('data/msb-paleon.3.0/OH_PLS_Witness_Trees/OH_PLS_Witness_Trees.shp')
pts <- mutate(pts.sf, x=st_coordinates(pts.sf)[,1],y=st_coordinates(pts.sf)[,2]) |> st_drop_geometry()
treenames <- c("SP1","SP2","SP3","SP4")
diams <- NA
dists <- NA
dataset='OH'
x = 'x'
y = 'y'
pts.OH <-  rackCorners(pts, dataset, treenames, diams, dists, x, y)
pts.OH <- st_as_sf(pts.OH, coords=c(x='x', y='y'), crs=crs(pts.sf))
pts.OH <- st_transform(pts.OH,crs=crs(pts.geo))

pts.geo <- rbind(pts.geo, pts.OH) |> unique()

pts.geo1 <- pts.geo |> mutate(dists=ifelse(is.na(tname) | dists %in% c(8888, 9999, 99999, 88888), NA, dists), dists=dists/100*66*0.3048, diam=diam*2.54) |> subset(!(is.na(tname)  & !treeseq %in% 1))
pts.geo1 <- pts.geo1 |> mutate(tname=tolower(tname))
# df <- unique(st_drop_geometry(subset(pts.geo1, select=c(dataset,tname))))

treelookup <- read.csv('treelookups1.csv') |> unique() |> mutate(tname=tolower(tname))
pts.geo1 <- pts.geo1 |> left_join(treelookup, by=join_by(dataset==dataset, tname==tname), multiple='first')
pts.geo1 <- pts.geo1 |> mutate(istree = ifelse((is.na(tname) | tname %in% '' | tolower(tname) %in% c(0,'no tree')) | Level2 %in% 'no tree',0,1))

#filter
#establish lat lon 
ptsll <- pts.geo1 |> st_transform(crs='epsg:4326')
ptsll <- ptsll |> mutate(lon = st_coordinates(ptsll)[,1],lat = st_coordinates(ptsll)[,2])|> st_drop_geometry() |> subset(select=c(id,lon,lat)) |> unique()
pts.geo1 <- pts.geo1 |> left_join(ptsll)
#take out northern sycamore
pts.geo1 <- pts.geo1 |> mutate(Species = ifelse(Species %in% 'Platanus' & lat >= 44, 'unk',Species))
#resolve cedar
anyass <- pts.geo1 |> st_drop_geometry() |> subset(Species %in% c('Fraxinus', 'Fraxinus nigra','Larix') | Level2 %in% c('Betula') | lat >= 44, select=id)
pts.geo1 <- pts.geo1 |> mutate(Species = ifelse(Species %in% 'JuniThuja' & id %in% anyass$id, 'Thuja',Species))
anyass <- pts.geo1 |> st_drop_geometry() |>  subset(Level2 %in% c('Quercus') | lat < 41, select=id)
pts.geo1 <- pts.geo1 |> mutate(Species = ifelse(Species %in% 'JuniThuja' & id %in% anyass$id, 'Juniperus',Species))
anyass <- pts.geo1 |> st_drop_geometry() |> subset(lat >= 42, select=id)
pts.geo1 <- pts.geo1 |> mutate(Species = ifelse(Species %in% 'JuniThuja' & id %in% anyass$id, 'Thuja',Species))

anyass <- pts.geo1 |> st_drop_geometry() |> subset(lat < 36 | Species %in% 'Fagus' , select=id)
pts.geo1 <- pts.geo1 |> mutate(Species = ifelse(Species %in% 'PopLirio' & id %in% anyass$id, 'Liriodendron',Species))


pts.geo1 <- pts.geo1 |> group_by(id) |> mutate(howmanytrees = sum(istree)) |> ungroup() 
pts.geo1 <- pts.geo1 |> subset(!(howmanytrees > 0 & istree == 0))
pts.geo1 <- pts.geo1 |> group_by(id) |> mutate(mbas = sum((diam/200)^2), mdist = sum(dists^2)) |> ungroup() |> mutate(BA = 10000*mbas/(mdist+0.0001))
pts.geo1 <- pts.geo1 |> mutate(BA = ifelse(dataset %in% "OH", NA, BA))

# write.csv(df, 'treelookups.csv', row.names = F)
saveRDS(pts.geo1, 'points/pts.geo1.RDS')
write_sf(pts.geo1,'gis/trees.shp')

notree <-  subset(pts.geo, grepl('no.tree', tolower(tname))| ((tname %in% '' | is.na(tname) ) & treeseq %in% 1))
plot(vect(notree))

write_sf(notree,'gis/notree.shp')