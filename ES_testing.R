library(sf)
library(terra)
library(dplyr)
library(ggplot2)
library(ranger)
library(vegan)
library(cluster)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
MLRA <- st_read('gis/MLRA_2017.shp')
pts <- readRDS('points/pts.geo1.RDS')


pts.vars <- readRDS('pts.vars.RDS')
vars90 <- rast('gis/vars90.tif')
vars270 <- rast('gis/vars270.tif')
# pca90 <- princomp(vars270, cor = T, maxcell=ncell(vars270)*0.1)
# pca90grid <- predict(pca90, pca90)
# writeRaster(pca90grid,'gis/pca90grid.tif', overwrite=T)
pts.vars90 <- readRDS('pts.vars90.RDS')
pts.vars270 <- readRDS('pts.vars270.RDS')
pts.vars90 <- pts.vars90 |> mutate(Level2 = ifelse(Species %in% 'Thuja', Species, Level2))
pts.vars270 <- pts.vars270 |> mutate(Level2 = ifelse(Species %in% 'Thuja', Species, Level2))
ttt <- pts.vars90 |> group_by(Level2) |> summarize(npts = length(Level2))
ttt <- subset(ttt, !is.na(Level2) & !Level2 %in% c('unk','no tree') & npts >= 10000)
ttt <- ttt$Level2
cors <- readRDS('output/cors.RDS')

#make species raster stack ----
for(i in 1:length(ttt)){
  x <- rast(paste0('gis/models90/',ttt[i],'.tif'))
  assign(ttt[i], x)         };rm(x)

Species <- rast(mget(ttt))

#make and extract sample points
MLRA <- st_transform(MLRA, crs(Species))
MLRA <- MLRA |> mutate(LRU = case_when(LRU %in% c('98A1','98A3','98A4','98A5') ~ '98A',
                                       LRU %in% c('97A1','97A2') ~ '97A',
                                       LRU %in% c('97B1','97B2') ~ '97B',
                                       TRUE ~ LRU))

MLRA94 <- subset(MLRA, MLRA %in% c("94A","94C"), select="LRU")
plot(MLRA94)
MLRA94.vect <- vect(MLRA94)
spts <- terra::spatSample(MLRA94.vect, size = 10000)
spts.mlra <- st_as_sf(spts) |> st_drop_geometry()
spts.spp <- terra::extract(Species, spts)
spts.env <- terra::extract(vars90, spts)

spts.es <- cbind(spts.spp, spts.env[,-1], spts.mlra) |> subset(!is.na(pH50) & !is.na(Tsuga))

#ES site key ---

spts.es <- spts.es |> mutate(ES = case_when(rockdepth < 150 ~ 'Limestone',
                                            floodfrq > 0 ~ 'Floodplain',
                                            histic > 0.5 &  hydric >= 0.5 ~ 'Mucks',
                                            sand150 >= 80 & sand50 >= 70 | sand50 >80 ~ 'Sandy',
                                            TRUE ~'Loamy'))

spts.es <- spts.es |> mutate(ES = case_when(ES %in%  'Floodplain' & hydric >= 0.5 ~ 'Wet Floodplain',
                                            ES %in%  'Floodplain' ~ 'Moist Floodplain',
                                            ES %in% 'Mucks' & pH50 >= 5 ~ "Euic Muck",
                                            ES %in% 'Mucks'  ~ "Acidic Peat",
                                            ES %in% 'Sandy' &  watertable > 100 ~ 'Dry Sand',
                                            ES %in% 'Sandy' &  hydric < 0.5 ~ 'Moist Sand',
                                            ES %in% 'Sandy' ~ 'Wet Sand',
                                            ES %in% 'Loamy' &  watertable > 100 ~ 'Dry Loam',
                                            ES %in% 'Loamy' &  hydric < 0.5 ~ 'Moist Loam',
                                            ES %in% 'Loamy' ~ 'Wet Loam',
                                            TRUE ~ ES))

spts.es <- spts.es |> mutate(ES = case_when(ES %in%  'Dry Sand' & LRU %in% '94AA' & Bhs >= 0.5 & spodic >= 0.5 ~ 'Rich Sand',
                                            ES %in%  'Dry Sand' & LRU %in% '94AA'  ~ 'Dry Sand',
                                            ES %in%  'Dry Sand' & (spodic >= 0.5 | pH50 >= 5 | carbdepth < 100) ~ 'Rich Sand',
                                            ES %in%  'Dry Sand'  ~ 'Dry Sand',
                                            ES %in%  'Moist Sand'  & (spodic >= 0.5 | pH50 <= 5) ~ 'Acidic Moist Sand',
                                            ES %in%  'Moist Sand'   ~ 'Euic Moist Sand',
                                            ES %in%  'Wet Sand'  & (spodic >= 0.5 | pH50 <= 5.5) ~ 'Acidic Wet Sand',
                                            ES %in%  'Wet Sand'   ~ 'Euic Wet Sand',
                                            TRUE ~ ES))
spts.es <- spts.es |> mutate(ES = ifelse(tan(slope)*100 >= 15 &  watertable > 100 & floodfrq <= 0 &  hydric < 0.5, paste(ES, "Slopes"), paste(ES, "Plains")))

spts.es <- spts.es |> mutate(ESD = paste(LRU, ES))

spts.es <- spts.es |> mutate(ESD = ifelse(rockdepth < 150, 'MLRA 94C Limestone Plains', ESD))


table(spts.es$ESD)

#silhouette index ----
spts.es.spp <- spts.es |> select(all_of(ttt))
spts.es.dist <- vegan::vegdist(spts.es.spp, method = 'bray', binary = F)
spts.es.groups <- spts.es |> mutate(group = as.integer(as.factor(ESD))) |> select(c(ID, group, ESD))
groups <- spts.es.groups |> select(c(group, ESD)) |> unique()

silindex <- cluster::silhouette(spts.es.groups$group, spts.es.dist)
summary(silindex)

si <- silindex |> as.data.frame()
si <- si |> left_join(data.frame(cluster=groups$group, cname = groups$ESD)) |>
  left_join(data.frame(neighbor=groups$group, nname = groups$ESD))


#simple clustering
spts.es.sum <- spts.es |> group_by(ESD) |> summarise(across(all_of(ttt), mean)) |> as.data.frame()
rownames(spts.es.sum) <-  spts.es.sum$ESD
spts.es.sum <- spts.es.sum[,-1]
spts.es.dist2 <- vegan::vegdist(spts.es.sum, method = 'bray', binary = F)
spts.es.t <- cluster::agnes(spts.es.dist2, method = 'ward')|> as.hclust()
plot(spts.es.t)


#K means clustering
kclusters <- terra::k_means(Species, centers=10, maxcell=100000)

plot(kclusters)

writeRaster(kclusters, 'gis/kclustersrast.tif', overwrite=T)

#PCA composition

treepca <- terra::prcomp(Species, maxcell=100000)
treepcarast <- predict(Species, treepca)
plot(treepcarast)
treepca$rotation
writeRaster(treepcarast, 'gis/treepcarast.tif')
treepca3 <- terra::princomp(Species, cor=T, maxcell=100000)
treepcarast3 <- predict(Species, treepca3)
plot(treepcarast3)

writeRaster(treepcarast3, 'gis/treepcarast3.tif')
treepca2 <- terra::princomp(Species, maxcell=Inf)
treepcarast2 <- predict(Species, treepca2)
plot(treepcarast2)
treepca2$scores
writeRaster(treepcarast2, 'gis/treepcarast2.tif')


library(visreg)

# ggpairs
# partial effects of dichotomous key variables
# #http://pbreheny.github.io/visreg/index.html
