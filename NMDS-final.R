
#NMDS ----
library(vegan)

library(climatools)
library(sf)
library(terra)
library(dplyr)

library(ggplot2)
library(ranger)
library(vegan)
library(cluster)
library(goeveg)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
MLRA <- st_read('gis/MLRA_2017.shp')
pts <- readRDS('points/pts.geo1.RDS')
seed1 = 4392 #Mount Rainier for resampling of Map
seed2 = 2010 #Mount LeConte for NMDS
seed3 = 4207 #Mauna Kea for clustering
seed4 = c(1345, 1917, 4809, 5636, 5895) #Ben Nevis, Mount Washington, Mont Blanc, Citlaltépetl, Kilimanjaro for KFold Cross Validation and Random Forest

#Load Data ----
set.seed(seed1)
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
  x <- rast(paste0('gis/finalrun_mean/',ttt[i],'.tif'))
  assign(ttt[i], x)         };rm(x)

Species <- rast(mget(ttt))

#make and extract sample points
MLRA <- st_transform(MLRA, crs(Species))
MLRA <- MLRA |> mutate(LRU = case_when(LRU %in% c('98A1','98A3','98A4','98A5') ~ '98A',
                                       LRU %in% c('97A1','97A2') ~ '97A',
                                       LRU %in% c('97B1','97B2') ~ '97B',
                                       TRUE ~ LRU))

MLRAall <- subset(MLRA, MLRA %in% c("94A","94C","96","97","98","99","111C") | LRU %in% "111BA", select="LRU")
MLRA98 <- subset(MLRA, MLRA %in% c("97","98","99","111C")| LRU %in% "111BA", select="LRU")
MLRA94 <- subset(MLRA, MLRA %in% c("94A","94C","96"), select="LRU")
MLRAselect <- MLRAall
plot(MLRAselect)
MLRA.vect <- vect(MLRAselect)
spts <- terra::spatSample(MLRA.vect, size = 5000)
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


essummary <- spts.es |> group_by(ESD) |> summarise(across(.fns=mean))
write.csv(essummary, 'essummary.csv', row.names = F)

# escorr1 <-  spts.es |> subset(ESD %in% '96A Acidic Moist Sand Plains') |> select_if(is.numeric) |>  cor(use = 'pairwise.complete.obs') |> as.data.frame() 
# escorr2 <-  spts.es |> subset(ESD %in% '96B Acidic Moist Sand Plains') |> select_if(is.numeric) |>  cor(use = 'pairwise.complete.obs') |> as.data.frame() 
#subset of non wet sites
spts.es2 <- spts.es #|> subset(hydric < 0.15 & watertable > 100 & sand150 < 80)

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



# #scree plot to show stress by number of dimensions
# set.seed(seed2)
# scrplt <- screeplot_NMDS(spts.es.spp,
#                          distance = "bray",
#                          k = 10,
#                          trymax = 20,
#                          autotransform = TRUE)
# 
# 
# scrplt
# ndim <- 4
# for  (fc in 2:ndim){
# score <- 100-(scrplt[fc])/scrplt[1]*100 ; names(score) <- fc
#   print(score)}
# # considering a north-south, pyro-mesic, and wet-dry gradients, should have 3 or more dimensions (4 made most logical kmeans classification with 2000 and 3000 points)



#Nonmetric Multidimensional Scaling (NMDS) ----
ndim <- 4
nmds <- metaMDS(spts.es.spp, k=ndim)
en <- envfit(nmds, spts.es.env, na.rm = TRUE, choices=c(1:ndim))
set.seed(seed2)
scores(nmds)

#join NMDS data components
pt.df <- scores(nmds, display='sites') |> as_tibble(rownames='sites') |> mutate(sites=as.numeric(sites)) |> inner_join(spts.es, by=join_by(sites==ID)) 
sp.df <- scores(nmds, display='species') |> as_tibble(rownames='species')
en.df <- scores(en, display='vectors')|> as_tibble(rownames='vectors')
print(en.df, n=nrow(en.df))

saveRDS(pt.df,'output/pt.df4.RDS')
saveRDS(sp.df,'output/sp.df4.RDS')
saveRDS(en.df,'output/en.df4.RDS')

#Continue ---- load NMDS 
library(vegan)

library(climatools)
library(sf)
library(terra)
library(dplyr)

library(ggplot2)
library(ranger)
library(vegan)
library(cluster)
library(goeveg)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
MLRA <- st_read('gis/MLRA_2017.shp')
pts <- readRDS('points/pts.geo1.RDS')
seed1 = 4392 #Mount Rainier for resampling of Map
seed2 = 2010 #Mount LeConte for NMDS
seed3 = 4207 #Mauna Kea for clustering
seed4 = c(1345, 1917, 4809, 5636, 5895) #Ben Nevis, Mount Washington, Mont Blanc, Citlaltépetl, Kilimanjaro for KFold Cross Validation and Random Forest

#Load Data ----
set.seed(seed1)
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
  x <- rast(paste0('gis/finalrun_mean/',ttt[i],'.tif'))
  assign(ttt[i], x)         };rm(x)

Species <- rast(mget(ttt))
#reload nmds results
ndim <- 4
pt.df3 <- readRDS('output/pt.df5.RDS')
sp.df3 <- readRDS('output/sp.df5.RDS')
en.df3 <- readRDS('output/en.df5.RDS')
pt.df <- readRDS('output/pt.df4.RDS')
sp.df <- readRDS('output/sp.df4.RDS')
en.df <- readRDS('output/en.df4.RDS')





##set up alternative cluster analyses to find optimal clusters ----
mtx <- pt.df[,2:(ndim+1)]
mtx3 <- pt.df3[,2:(5+1)]
dst <- dist(mtx)
dst3 <- dist(mtx3)
mtx2 <- pt.df[,names(Species)]^0.5 
dst2 <- vegan::vegdist(mtx2, method = 'bray')
mtx2a <- pt.df3[,names(Species)]^0.5 
dst2a <- vegan::vegdist(mtx2a, method = 'bray')
htree0 <- dst2 |> agnes(method = 'ward')
htree1 <- dst |> agnes(method = 'ward')
htree2 <- dst |> agnes(method = 'average') 
htree3 <- dst |> diana()
for(i in 2:10){
set.seed(seed3)
nclust <- i


kc0 <- dst2 |> kmeans(nclust) 
kc0 <- kc0$cluster
silscor <- cluster::silhouette(kc0, dist = dst2)
msilscor0 <- mean(silscor[,3])

kc1 <- mtx |> kmeans(nclust) 
kc1 <- kc1$cluster
silscor <- cluster::silhouette(kc1, dist = dst2)
msilscor1 <- mean(silscor[,3])

kc1a <- mtx3 |> kmeans(nclust) 
kc1a <- kc1a$cluster
silscor <- cluster::silhouette(kc1a, dist = dst2a)
msilscor1a <- mean(silscor[,3])

hc0 <- cutree(htree0, nclust)
hsilscor <- cluster::silhouette(hc0, dist = dst2)
hmsilscor0 <- mean(hsilscor[,3])

hc1 <- cutree(htree1, nclust)
hsilscor <- cluster::silhouette(hc1, dist = dst2)
hmsilscor1 <- mean(hsilscor[,3])


hc1 <- cutree(htree2, nclust)
hsilscor <- cluster::silhouette(hc1, dist = dst2)
hmsilscor2 <- mean(hsilscor[,3])

hc1 <- cutree(htree3, nclust)
hsilscor <- cluster::silhouette(hc1, dist = dst2)
hmsilscor3 <- mean(hsilscor[,3])

scdf0 <- data.frame(nclust = nclust, 
                    kmeansraw = msilscor0, 
                    kmeans4d = msilscor1, 
                    kmeans3d = msilscor1a, 
                    wardraw = hmsilscor0, 
                    ward = hmsilscor1, 
                    upgma = hmsilscor2,
                    diana = hmsilscor3)
if(i==2){scdf <- scdf0}else{scdf <- rbind(scdf,scdf0)}
};rm(scdf0)

clustplot <- ggplot(scdf)+
  geom_line(aes(x=nclust, y=kmeans4d, color='kmeans4d'))+
  geom_line(aes(x=nclust, y=kmeans3d, color='kmeans5d'))+
  geom_line(aes(x=nclust, y=kmeansraw, color='kmeansraw'))+
  geom_line(aes(x=nclust, y=ward, color='ward'))+
  geom_line(aes(x=nclust, y=wardraw, color='wardraw'))+
  geom_line(aes(x=nclust, y=upgma, color='upgma'))+
  geom_line(aes(x=nclust, y=diana, color='diana'))+
  geom_point(aes(x=nclust, y=kmeans4d, color='kmeans4d'))+
  geom_point(aes(x=nclust, y=kmeans3d, color='kmeans5d'))+
  geom_point(aes(x=nclust, y=kmeansraw, color='kmeansraw'))+
  geom_point(aes(x=nclust, y=ward, color='ward'))+
  geom_point(aes(x=nclust, y=wardraw, color='wardraw'))+
  geom_point(aes(x=nclust, y=upgma, color='upgma'))+
  geom_point(aes(x=nclust, y=diana, color='diana'))+
  scale_y_continuous(name='mean silhouette')+
  scale_x_continuous(name='number of clusters', breaks = c(1:10), minor_breaks = NULL)+
  scale_color_manual(name='method',
                     labels =c('kmeans4d','kmeans5d','kmeansraw','wardraw','ward','upgma','diana'), 
                     values =c('black','blue','orange','magenta','red','green','yellow'))
clustplot

#Implement the clustering
set.seed(seed3)
nclust <- 8
kc1 <- mtx |> kmeans(nclust) 
kc1 <- kc1$cluster
pt.df <- pt.df |> mutate(kc = as.factor(paste0('cluster',kc1)))
# kc1 <- mtx2 |> kmeans(nclust) 
# kc1 <- kc1$cluster
# pt.df <- pt.df |> mutate(kc = as.factor(paste0('cluster',kc1)))
# htree <- dst2 |> agnes(method = 'ward')
# plot(as.hclust(htree))
# hc1 <- cutree(htree, nclust)
# pt.df <- pt.df |> mutate(kc = as.factor(paste0('cluster',hc1)))

clustersummary <- pt.df |> group_by(kc) |> summarise(across(.fns=mean))
groups <- unique(pt.df$kc) |> as.character()
clustercorr <- pt.df
for(i in 1:length(groups)){
clustercorr <- clustercorr |> mutate(x = ifelse(kc %in% groups[i],1,0))
colnames(clustercorr)[colnames(clustercorr) %in% 'x'] <- groups[i]
}
clustercorr <- clustercorr |> select_if(is.numeric) |>  cor(use = 'pairwise.complete.obs') |> as.data.frame() 
clustercorr <- clustercorr[,(ncol(clustercorr)-nclust+1):ncol(clustercorr)]

# clusttrans <- t(clustersummary)
# colnames(clusttrans) <- clustersummary$kc
# clusttrans <- clusttrans[rownames(clusttrans) %in% names(Species),] |> as.data.frame() 
# clusttrans <- clusttrans |> mutate(across(.fns = as.numeric)) 
# clusttrans <- clusttrans |> mutate(across(1:6,\(x).fns = round(x, 3)))


gp <- ggplot() +
  geom_point(data=pt.df, aes(x=NMDS1,y=NMDS2, color=kc), alpha=0.5, size=2)+
  geom_point(data=sp.df, aes(x=NMDS1,y=NMDS2), color='blue')+
  geom_text(data=sp.df, aes(label=species, x=NMDS1,y=NMDS2), vjust = 0, nudge_y = 0.02, nudge_x = 0.05, color='blue')+
  geom_segment(data=en.df, aes(x=0,y=0,xend=NMDS1,yend=NMDS2), arrow = arrow(length = unit(0.03, "npc")), color='black')+
  geom_text(data=en.df, aes(label=vectors, x=NMDS1,y=NMDS2), vjust = 0, nudge_y = 0.02, nudge_x = 0.05, color='black')

gp

png('plots/MNDS_biplot_MLRA9694AC.png', width = 800, height = 800)
gp
dev.off()

gp2 <- ggplot() +
  geom_point(data=pt.df, aes(x=NMDS1,y=NMDS3, color=kc), alpha=0.5, size=2)+
  geom_point(data=sp.df, aes(x=NMDS1,y=NMDS3), color='blue')+
  geom_text(data=sp.df, aes(label=species, x=NMDS1,y=NMDS3), vjust = 0, nudge_y = 0.02, nudge_x = 0.05, color='blue')+
  geom_segment(data=en.df, aes(x=0,y=0,xend=NMDS1,yend=NMDS3), arrow = arrow(length = unit(0.03, "npc")), color='black')+
  geom_text(data=en.df, aes(label=vectors, x=NMDS1,y=NMDS3), vjust = 0, nudge_y = 0.02, nudge_x = 0.05, color='black')

gp2




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
pt.df <- pt.df |> mutate(ESG = case_when(
  kc %in% 'cluster1' ~  'Cool',
  kc %in% 'cluster2' ~  'Warm Pyrophytic',
  kc %in% 'cluster3' ~  'Warm Mesophytic',
  kc %in% 'cluster4' ~  'Warm Pyrophytic',
  kc %in% 'cluster5' ~  'Cool Pyrophytic',
  kc %in% 'cluster6' ~  'Cool Wet',
  TRUE ~ 'other'))

rp <- rpart(kc ~ p+e+Twh+Tw+Tc+Tcl+Tg+Tg30+
              s+d+m+
              Bhs+carbdepth+clay150+floodfrq+histic+humic+humicdepth+
              hydric+ksatdepth+OM150+pH50+rockdepth+sand150+sand50+spodic+watertable+
              slope+slope500+popen+nopen+solar, data = subset(pt.df),
            method="class", control = list(maxdepth = 5, cp=0.001, minsplit=50))
png(filename="nmds_rpart.png",width = 10, height = 3, units = 'in', res = 600)
rpart.plot(rp, extra=108,legend.cex=0.5, digits=2)
dev.off()

af <- 2
gp <- ggplot() +
  geom_point(data=pt.df, aes(x=NMDS1,y=NMDS2, color=ESG), alpha=0.2, size=2)+
  geom_point(data=sp.df, aes(x=NMDS1,y=NMDS2), color='black')+
  geom_text(data=sp.df, aes(label=species, x=NMDS1,y=NMDS2), vjust = 0, nudge_y = 0.02, nudge_x = 0.05, color='blue')+
  geom_segment(data=en.df, aes(x=0,y=0,xend=NMDS1*af,yend=NMDS2*af), arrow = arrow(length = unit(0.03, "npc")), color='red')+
  geom_text(data=en.df, aes(label=vectors, x=NMDS1*af,y=NMDS2*af), vjust = 0, nudge_y = 0.02, nudge_x = 0.05, color='darkred')+
  scale_color_manual("ESG", values = c("Warm Mesophytic" = "#7CAE00",
                                         "Warm Pyrophytic" = "#E68613",
                                         "Cool" = "#00BFC4"))

gp

# library(scales) 
# show_col(hue_pal()(16))
# c("#F8766D", "#B79F00", "#00BA38", "#00BFC4", "#619CFF", "#F564E3",
#   "#F8766D", "#E68613", "#CD9600", "#ABA300",
#   "#7CAE00", "#0BB702", "#00BE67", "#00C19A",
#   "#00BFC4", "#00B8E7", "#00A9FF", "#8494FF",
#   "#C77CFF", "#ED68ED", "#FF61CC", "#FF68A1")



#make NMDSclusters on map Test

pt.df$kc2 <- as.numeric(pt.df$kc)
clss <- unique(pt.df$kc2)
nClasses = length(clss)
# class labels
lev = data.frame(Value=1:nClasses, class=clss)

Species2 <- resample(Species, aggregate(Species$Abies, fact=5), method='near')

formular <- as.formula(paste("pos",paste(names(Species), collapse = " + ", sep = ""), sep = " ~ "))
for (i in 1:length(clss)){#i=1
pt.df$pos <- ifelse(pt.df$kc2 %in% clss[i], 1,0)


set.seed(seed3)
rf <-  ranger(formular, data = pt.df, classification = F)
named <- paste0('class',clss[i])
rf.prediction <-  predict(Species2, rf, na.rm=T);  names(rf.prediction) <- named
plot(rf.prediction)
assign(named, rf.prediction)
}
esgprobs <- rast(mget(paste0('class',clss)))

# Number of probability layers (should be equal to the number of predicted classes)
nl = nlyr(esgprobs)
# Number of cells (pixels in probPred object)
nc = ncell(esgprobs)
prob1 = max(esgprobs)
classI = which.lyr(esgprobs == prob1)
plot(classI)






classPred1 = classify(classI, data.frame(cbind(1:nl, lev$Value)))
lev <- data.frame(Value=1:8, class=c('Basswood','Elm','Oak','Hemlock','Hickory','Whitecedar','Pine','Tamarack'))
levels(classI) = lev
plot(classI, col=c('green','purple','orange','darkcyan','yellow','cyan','red','blue'))
plot(prob1)
esgprobs2 <- esgprobs/(sum(esgprobs))
entropy <- -sum(esgprobs*log(esgprobs2+0.0001))
plot(entropy)


writeRaster(classI, paste0('gis/test','Classes','.tif'), overwrite=T)

n <- 25000 #total number of positive samples
ntest = 0.2
#Final High Def run
classnames=c('Basswood','Elm','Oak','Hemlock','Hickory','Whitecedar','Pine','Tamarack')
for (i in 1:length(clss)){
  for (k in 1:length(seed4)){
    #i=1;k=1
    set.seed(seed4[k])
    train0 <- pt.df |> mutate(pos = ifelse(kc2 %in% clss[i], 1,0))
    positives <- subset(train0, pos %in% 1)
    negatives <- subset(train0, pos %in% 0)
    takeout.p <- sample(1:nrow(positives), nrow(positives)*ntest)
    takeout.n <- sample(1:nrow(negatives), nrow(negatives)*ntest)
    train.p <- positives[-takeout.p,]
    test.p <- positives[takeout.p,]
    train.n <- negatives[-takeout.n,]
    test.n <- negatives[takeout.n,]
    if(TRUE){#Random Over Sample
      train.p <- train.p[sample(1:nrow(train.p), n, replace = TRUE),]
      train.n <- train.n[sample(1:nrow(train.n), n, replace = TRUE),]
    }
    if(TRUE){#Random Over Sample
      test.p <- test.p[sample(1:nrow(test.p), n, replace = TRUE),]
      test.n <- test.n[sample(1:nrow(test.n), n, replace = TRUE),]
    }
    train <- rbind(train.p,train.n)
    test <- rbind(test.p,test.n)
    rf <-  ranger(formular, data = train, classification = F)
    named <- paste0(classnames[i])
    train.rf <- train |> mutate(prediction = predictions(predict(rf, train, na.rm=T)))
    test.rf <- test |> mutate(prediction = predictions(predict(rf, test, na.rm=T)))
    
    
    modmets <- data.frame(ESG = classnames[i],
                          model = 'RF',
                          seed = seed4[k],
                          AUCtrain = Metrics::auc(actual=train.rf$pos, predicted=train.rf$prediction),
                          AUCtest = Metrics::auc(actual=test.rf$pos, predicted=test.rf$prediction),
                          maxKappatrain = maxKappa(actual=train.rf$pos, predicted=train.rf$prediction),
                          maxKappatest = maxKappa(actual=test.rf$pos, predicted=test.rf$prediction))
    if(i==1 & k==1){Kfold = modmets}else{Kfold = rbind(Kfold, modmets)}
    


    
    rf.prediction <-  predict(Species, rf, na.rm=T);  names(rf.prediction) <- named
    writeRaster(rf.prediction, paste0('gis/clusters/',named,'.',seed4[k],'.tif'), overwrite=T)
  }
  write.csv(Kfold, 'clusterkfold.csv', row.names = FALSE)
  }


#Continue ---- load cluster maps
library(vegan)

library(climatools)
library(sf)
library(terra)
library(dplyr)

library(ggplot2)
library(ranger)
library(vegan)
library(cluster)
library(goeveg)
library(stringr)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
files <- data.frame(filename = list.files('gis/clusters'))
files <- files |> mutate(taxon = str_split_fixed(filename, '\\.', 2)[,1], model = str_split_fixed(filename, '\\.', 3)[,2])
ttt = unique(files$taxon)
for (i in 1:length(ttt)){#i=1
  f = subset(files, taxon %in% ttt[i])
  mmm <- f$model
  for(k in 1:length(mmm)){#k=1
    fn <- f$filename[k]
    
    r <- rast(paste0('gis/clusters/',fn))
    names(r) <- mmm[k]
    assign(mmm[k], r)
  }
  r <- rast(mget(mmm))
  r.sd <- stdev(r)
  r.m <- mean(r)
  name1 <- paste0(ttt[i],'_sd')
  name2 <- ttt[i]
  names(r.sd) <- name1
  names(r.m) <- name2
  writeRaster(r.sd, paste0('gis/clusters_sd/', name1,'.tif'), overwrite =TRUE)
  writeRaster(r.m, paste0('gis/clusters_mean/', name2,'.tif'), overwrite =TRUE)
}


for(i in 1:length(ttt)){
  x <- rast(paste0('gis/clusters_mean/',ttt[i],'.tif'))
  assign(ttt[i], x)         };rm(x)

Clusters <- rast(mget(ttt))

for(i in 1:length(ttt)){
  x <- rast(paste0('gis/clusters_sd/',ttt[i],'_sd.tif'))
  assign(paste0(ttt[i],'_sd'), x)         };rm(x)

Species_sd <- rast(mget(paste0(ttt,'_sd')))
meanSD <- mean(Species_sd)
plot(meanSD)

Entropy <- -sum(Clusters*log(Clusters+0.0001))
plot(Entropy)



