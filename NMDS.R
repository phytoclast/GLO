
#NMDS ----
library(vegan)
set.seed(0)

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

MLRA94 <- subset(MLRA, MLRA %in% c("94A","94C","96","97","98","99","111C") | LRU %in% "111BA", select="LRU")
# MLRA94 <- subset(MLRA, MLRA %in% c("94A","94C"), select="LRU")
plot(MLRA94)
MLRA94.vect <- vect(MLRA94)
spts <- terra::spatSample(MLRA94.vect, size = 3000)
spts.mlra <- st_as_sf(spts) |> st_drop_geometry()
spts.spp <- terra::extract(Species, spts)
spts.env <- terra::extract(vars90, spts)

spts.es <- cbind(spts.spp, spts.env[,-1], spts.mlra) |> subset(!is.na(pH50) & !is.na(Tsuga))
spts.es <- spts.es |> mutate(spptotal = apply(spts.es[,ttt], MARGIN=1, FUN='sum')) |> subset(spptotal > 0)

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

#subset of non wet sites
spts.es2 <- spts.es #|> subset(hydric < 0.15 & watertable > 100 & sand150 >= 80)

colnames(spts.es)
spts.es.spp <- spts.es2 |> select(all_of(ttt))
spts.es.env <- spts.es2 |> select(c("p","m","e","s","d","Twh","Tw","Tc","Tcl","Tg",
                                   "Tg30","Bhs","carbdepth","clay150","floodfrq","histic",
                                   "humic","humicdepth","hydric","ksatdepth","OM150","pH50","rockdepth","sand150",
                                   "sand50","spodic","watertable","elev","slope","aspect",
                                   "slope500","shade","popen","nopen","solar"))
spts.es.dist <- vegan::vegdist(spts.es.spp, method = 'bray', binary = F)
spts.es.groups <- spts.es2 |> mutate(group = as.integer(as.factor(ESD))) |> select(c(ID, group, ESD))
groups <- spts.es.groups |> select(c(group, ESD)) |> unique()
dist_tbl <- as_tibble(spts.es.dist, rownames="samples")

#considering a north-south, pyro-mesic, and wet-dry gradients, should have 3 or more dimensions (4 made most logical kmeans classification with 2000 and 3000 points)
ndim <- 4
nmds <- metaMDS(spts.es.spp, k=ndim)
en <- envfit(nmds, spts.es.env, na.rm = TRUE, choices=c(1:ndim))

scores(nmds)


pt.df <- scores(nmds, display='sites') |> as_tibble(rownames='sites') |> mutate(sites=as.numeric(sites)) |> inner_join(spts.es, by=join_by(sites==ID)) 
sp.df <- scores(nmds, display='species') |> as_tibble(rownames='species')
en.df <- scores(en, display='vectors')|> as_tibble(rownames='vectors')
print(en.df, n=nrow(en.df))

# library(plotly)
# gp <- ggplot() +
#   geom_point(data=pt.df, aes(x=NMDS1,y=NMDS2), alpha=0.1)+
#   geom_point(data=sp.df, aes(x=NMDS1,y=NMDS2), color='blue')+
#   geom_text(data=sp.df, aes(label=species, x=NMDS1,y=NMDS2), vjust = 0, nudge_y = 0.02, nudge_x = 0.05, color='blue')+
#   geom_segment(data=en.df, aes(x=0,y=0,xend=NMDS1,yend=NMDS2), arrow = arrow(length = unit(0.03, "npc")), color='red')+
#   geom_text(data=en.df, aes(label=vectors, x=NMDS1,y=NMDS2), vjust = 0, nudge_y = 0.02, nudge_x = 0.05, color='red')

#kmeans of 6 classes for all Lower Michigan MLRAs corresponds with northern hardwoods, Oaks, pines, southern hardwoods, southern swamp, and northern swamp.
nclust <- 6
kc1 <- pt.df[,2:(ndim+1)] |> kmeans(nclust) 
kc1 <- kc1$cluster
pt.df <- pt.df |> mutate(kc = as.factor(paste0('cluster',kc1)))
htree <- vegdist(pt.df[,2:(ndim+1)], method='euclidean') |> agnes(method = 'ward') 
plot(as.hclust(htree))
hc1 <- cutree(htree, nclust)
pt.df <- pt.df |> mutate(hc = as.factor(paste0('cluster',hc1)))

clustersummary <- pt.df |> group_by(kc) |> summarise(across(.fns=mean))
groups <- unique(pt.df$kc) |> as.character()
clustercorr <- pt.df
for(i in 1:length(groups)){
clustercorr <- clustercorr |> mutate(x = ifelse(kc %in% groups[i],1,0))
colnames(clustercorr)[colnames(clustercorr) %in% 'x'] <- groups[i]
}
clustercorr <- clustercorr |> select_if(is.numeric) |>  cor(use = 'pairwise.complete.obs') |> as.data.frame() 
clustercorr <- clustercorr[,(ncol(clustercorr)-nclust+1):ncol(clustercorr)]

clusttrans <- t(clustersummary)
colnames(clusttrans) <- clustersummary$kc
clusttrans <- clusttrans[rownames(clusttrans) %in% names(Species),] |> as.data.frame() 
clusttrans <- clusttrans |> mutate(across(.fns = as.numeric)) 
clusttrans <- clusttrans |> mutate(across(1:6,\(x).fns = round(x, 3)))


gp <- ggplot() +
  geom_point(data=pt.df, aes(x=NMDS1,y=NMDS2, color=hc), alpha=0.5, size=2)+
  geom_point(data=sp.df, aes(x=NMDS1,y=NMDS2), color='blue')+
  geom_text(data=sp.df, aes(label=species, x=NMDS1,y=NMDS2), vjust = 0, nudge_y = 0.02, nudge_x = 0.05, color='blue')+
  geom_segment(data=en.df, aes(x=0,y=0,xend=NMDS1,yend=NMDS2), arrow = arrow(length = unit(0.03, "npc")), color='red')+
  geom_text(data=en.df, aes(label=vectors, x=NMDS1,y=NMDS2), vjust = 0, nudge_y = 0.02, nudge_x = 0.05, color='red')

gp

png('MNDS_biplot_MLRA94AC.png', width = 800, height = 800)
gp
dev.off()


library(rpart)
library(rpart.plot)
unique(pt.df$kc)
pt.df <- pt.df |> mutate(ESG = case_when(
  kc %in% 'cluster1' ~  'Warm Mesophytic',
  kc %in% 'cluster2' ~  'Cool Mesophytic',
  kc %in% 'cluster3' ~  'Warm Wet',
  kc %in% 'cluster4' ~  'Warm Pyrophytic',
  kc %in% 'cluster5' ~  'Cool Pyrophytic',
  kc %in% 'cluster6' ~  'Cool Wet',
  TRUE ~ 'other'))
    
rp <- rpart(ESG ~ p+e+Twh+Tw+Tc+Tcl+Tg+Tg30+
              s+d+m+
              Bhs+carbdepth+clay150+floodfrq+histic+humic+humicdepth+
              hydric+ksatdepth+OM150+pH50+rockdepth+sand150+sand50+spodic+watertable+
              slope+slope500+popen+nopen+solar, data = subset(pt.df),
            method="class", control = list(maxdepth = 5, cp=0.01, minsplit=100))
png(filename="nmds_rpart.png",width = 10, height = 3, units = 'in', res = 600)
rpart.plot(rp, extra=108,legend.cex=0.5, digits=2)
dev.off()


gp <- ggplot() +
  geom_point(data=pt.df, aes(x=NMDS1,y=NMDS2, color=ESG), alpha=0.5, size=2)+
  geom_point(data=sp.df, aes(x=NMDS1,y=NMDS2), color='blue')+
  geom_text(data=sp.df, aes(label=species, x=NMDS1,y=NMDS2), vjust = 0, nudge_y = 0.02, nudge_x = 0.05, color='blue')+
  geom_segment(data=en.df, aes(x=0,y=0,xend=NMDS1,yend=NMDS2), arrow = arrow(length = unit(0.03, "npc")), color='red')+
  geom_text(data=en.df, aes(label=vectors, x=NMDS1,y=NMDS2), vjust = 0, nudge_y = 0.02, nudge_x = 0.05, color='red')

gp

gp <- ggplot() +
  geom_point(data=pt.df, aes(x=NMDS2,y=NMDS4, color=ESG), alpha=0.5, size=2)+
  geom_point(data=sp.df, aes(x=NMDS2,y=NMDS4), color='blue')+
  geom_text(data=sp.df, aes(label=species, x=NMDS2,y=NMDS4), vjust = 0, nudge_y = 0.02, nudge_x = 0.05, color='blue')+
  geom_segment(data=en.df, aes(x=0,y=0,xend=NMDS2,yend=NMDS4), arrow = arrow(length = unit(0.03, "npc")), color='red')+
  geom_text(data=en.df, aes(label=vectors, x=NMDS2,y=NMDS4), vjust = 0, nudge_y = 0.02, nudge_x = 0.05, color='red')

gp


#make NMDSclusters on map
#hierarchical
pt.df$hc2 <- as.numeric(pt.df$hc)

formular <- as.formula(paste("hc2",paste(names(Species), collapse = " + ", sep = ""), sep = " ~ "))
rf <-  ranger(formular, data = pt.df, classification = T)
named <- 'NMDSClusters.hierarchical'
rf.prediction <-  predict(Species, rf, na.rm=T);  names(rf.prediction) <- named
plot(rf.prediction)
writeRaster(rf.prediction, paste0('gis/',named,'.tif'), overwrite=T)
#kmeans
pt.df$kc2 <- as.numeric(pt.df$kc)

formular <- as.formula(paste("kc2",paste(names(Species), collapse = " + ", sep = ""), sep = " ~ "))

rf <-  ranger(formular, data = pt.df, classification = T)
named <- 'NMDSClusters.kmeans'
rf.prediction <-  predict(Species, rf, na.rm=T);  names(rf.prediction) <- named
plot(rf.prediction)
writeRaster(rf.prediction, paste0('gis/',named,'.tif'), overwrite=T)