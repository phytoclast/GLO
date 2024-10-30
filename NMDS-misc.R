
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

set.seed(0)
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

MLRAall <- subset(MLRA, MLRA %in% c("94A","94C","96","97","98","99","111C") | LRU %in% "111BA", select="LRU")
MLRA98 <- subset(MLRA, MLRA %in% c("97","98","99","111C")| LRU %in% "111BA", select="LRU")
MLRA94 <- subset(MLRA, MLRA %in% c("94A","94C","96"), select="LRU")
MLRAselect <- MLRA94
plot(MLRAselect)
MLRA.vect <- vect(MLRAselect)
spts <- terra::spatSample(MLRA.vect, size = 2000)
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

#scree plot to show stress by number of dimensions
# scrplt <- screeplot_NMDS(spts.es.spp,
#                          distance = "bray",
#                          k = 10,
#                          trymax = 20,
#                          autotransform = TRUE)
# scrplt
# ndim <- 8
# for  (fc in 2:ndim){
# score <- 100-(scrplt[fc])/scrplt[1]*100 ; names(score) <- fc
#   print(score)}
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

#silouette for optimal clusters ----

mtx <- pt.df[,2:(ndim+1)]
dst <- dist(mtx)
mtx2 <- pt.df[,names(Species)]^0.5
dst2 <- vegan::vegdist(mtx2, method = 'bray')
htree0 <- dst2 |> agnes(method = 'ward')
htree1 <- dst |> agnes(method = 'ward')
htree2 <- dst |> agnes(method = 'average') 
htree3 <- dst |> diana()
for(i in 2:10){
set.seed(0)
nclust <- i


kc0 <- mtx2 |> kmeans(nclust) 
kc0 <- kc0$cluster
silscor <- cluster::silhouette(kc0, dist = dst)
msilscor0 <- mean(silscor[,3])

kc1 <- mtx |> kmeans(nclust) 
kc1 <- kc1$cluster
silscor <- cluster::silhouette(kc1, dist = dst)
msilscor1 <- mean(silscor[,3])

hc0 <- cutree(htree0, nclust)
hsilscor <- cluster::silhouette(hc0, dist = dst)
hmsilscor0 <- mean(hsilscor[,3])

hc1 <- cutree(htree1, nclust)
hsilscor <- cluster::silhouette(hc1, dist = dst)
hmsilscor1 <- mean(hsilscor[,3])


hc1 <- cutree(htree2, nclust)
hsilscor <- cluster::silhouette(hc1, dist = dst)
hmsilscor2 <- mean(hsilscor[,3])


hc1 <- cutree(htree3, nclust)
hsilscor <- cluster::silhouette(hc1, dist = dst)
hmsilscor3 <- mean(hsilscor[,3])

scdf0 <- data.frame(nclust = nclust, 
                    kmeansraw = msilscor0, 
                    kmeans = msilscor1, 
                    wardraw = hmsilscor0, 
                    ward = hmsilscor1, 
                    upgma = hmsilscor2,
                    diana = hmsilscor3)
if(i==2){scdf <- scdf0}else{scdf <- rbind(scdf,scdf0)}
};rm(scdf0)

ggplot(scdf)+
  geom_line(aes(x=nclust, y=kmeans, color='kmeans'))+
  geom_line(aes(x=nclust, y=kmeansraw, color='kmeansraw'))+
  geom_line(aes(x=nclust, y=ward, color='ward'))+
  geom_line(aes(x=nclust, y=wardraw, color='wardraw'))+
  geom_line(aes(x=nclust, y=upgma, color='upgma'))+
  geom_line(aes(x=nclust, y=diana, color='diana'))+
  geom_point(aes(x=nclust, y=kmeans, color='kmeans'))+
  geom_point(aes(x=nclust, y=kmeansraw, color='kmeansraw'))+
  geom_point(aes(x=nclust, y=ward, color='ward'))+
  geom_point(aes(x=nclust, y=wardraw, color='wardraw'))+
  geom_point(aes(x=nclust, y=upgma, color='upgma'))+
  geom_point(aes(x=nclust, y=diana, color='diana'))+
  scale_y_continuous(name='mean silhouette')+
  scale_x_continuous(name='number of clusters', breaks = c(1:10), minor_breaks = NULL)+
  scale_color_manual(name='method',
                     labels =c('kmeans','kmeansraw','wardraw','ward','upgma','diana'), 
                     values =c('black','orange','magenta','red','green','blue'))
set.seed(0)
nclust <- 3
kc1 <- mtx2 |> kmeans(nclust) 
kc1 <- kc1$cluster
pt.df <- pt.df |> mutate(kc = as.factor(paste0('cluster',kc1)))
# htree <- dst |> agnes(method = 'ward')
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


#Variable Importance
pt.df -> x
x$pyrophytic <- ifelse(x$kc %in% 'cluster2',1,0)
x$mesophytic <- ifelse(x$kc %in% 'cluster3',1,0)
x$hydrophytic <- ifelse(x$kc %in% 'cluster1',1,0)
x$upland <- ifelse(x$kc %in% c('cluster2','cluster3'),1,0)


checkvars <- c("p","m","e","s","d","Twh","Tw","Tc","Tcl","Tg",
               "Tg30","Bhs","carbdepth","clay150","floodfrq","histic",
               "humic","humicdepth","hydric","ksatdepth","OM150","pH50","rockdepth","sand150",
               "sand50","spodic","watertable","elev","slope","aspect",
               "slope500","shade","popen","nopen","solar")

tableofvars <- find.multithreshold(x=x, class1='mesophytic',class2='pyrophytic',variables=checkvars)


checkvars <- names(Species)
tableofspp <- find.multithreshold(x=x, class1='mesophytic',class2='pyrophytic',variables=checkvars)

checkvars <- c("p","m","e","s","d","Twh","Tw","Tc","Tcl","Tg",
               "Tg30","Bhs","carbdepth","clay150","floodfrq","histic",
               "humic","humicdepth","hydric","ksatdepth","OM150","pH50","rockdepth","sand150",
               "sand50","spodic","watertable","elev","slope","aspect",
               "slope500","shade","popen","nopen","solar")

tableofvars2 <- find.multithreshold(x=x, class1='hydrophytic',class2='upland',variables=checkvars)


checkvars <- names(Species)
tableofspp2 <- find.multithreshold(x=x, class1='hydrophytic',class2='upland',variables=checkvars)

tableall <- rbind(tableofvars, tableofspp)
tableall2 <- rbind(tableofvars2, tableofspp2)

var.nfogain.meso <- 
  ggplot(data=tableofvars, aes(x=infogain,y=reorder(variable, infogain)))+
  geom_bar(stat="identity")+
  scale_y_discrete(name='Variable')+
  scale_x_continuous(name='Information Gain')
var.cor.meso <- 
  ggplot(data=tableofvars, aes(x=class1cor,y=reorder(variable, infogain)))+
  geom_bar(stat="identity")+
  scale_y_discrete(name='Variable')+
  scale_x_continuous(name='Correlation with mesophytic (relative to pyrophytic)')
spp.infogain.meso <- 
  ggplot(data=tableofspp, aes(x=infogain,y=reorder(variable, infogain)))+
  geom_bar(stat="identity")+
  scale_y_discrete(name='Taxon')+
  scale_x_continuous(name='Information Gain')
spp.cor.meso <- 
  ggplot(data=tableofspp, aes(x=class1cor,y=reorder(variable, infogain)))+
  geom_bar(stat="identity")+
  scale_y_discrete(name='Taxon')+
  scale_x_continuous(name='Correlation with mesophytic (relative to pyrophytic)')
var.gini.meso <- 
  ggplot(data=tableofvars, aes(x=gini,y=reorder(variable, infogain)))+
  geom_bar(stat="identity")+
  scale_y_discrete(name='Variable')+
  scale_x_continuous(name='Gini')
var.entropy.meso <- 
  ggplot(data=tableofvars, aes(x=entropy,y=reorder(variable, infogain)))+
  geom_bar(stat="identity")+
  scale_y_discrete(name='Variable')+
  scale_x_continuous(name='Entropy')

var.infogain.hydro <- 
  ggplot(data=tableofvars2, aes(x=infogain,y=reorder(variable, infogain)))+
  geom_bar(stat="identity")+
  scale_y_discrete(name='Variable')+
  scale_x_continuous(name='Information Gain')
spp.infogain.hydro <- 
  ggplot(data=tableofspp2, aes(x=infogain,y=reorder(variable, infogain)))+
  geom_bar(stat="identity")+
  scale_y_discrete(name='Taxon')+
  scale_x_continuous(name='Information Gain')
var.cor.hydro <- 
  ggplot(data=tableofvars2, aes(x=class1cor,y=reorder(variable, infogain)))+
  geom_bar(stat="identity")+
  scale_y_discrete(name='Variable')+
  scale_x_continuous(name='Correlation with hydrophytic (relative to upland)')
spp.cor.hydro <- 
  ggplot(data=tableofspp2, aes(x=class1cor,y=reorder(variable, infogain)))+
  geom_bar(stat="identity")+
  scale_y_discrete(name='Taxon')+
  scale_x_continuous(name='Correlation with hydrophytic (relative to upland)')


png('plots/var.nfogain.meso.png', width = 800, height = 800)
var.nfogain.meso
dev.off()
png('plots/var.cor.meso.png', width = 800, height = 800)
var.cor.meso
dev.off()
png('plots/spp.infogain.meso.png', width = 800, height = 800)
spp.infogain.meso
dev.off()
png('plots/var.entropy.meso.png', width = 800, height = 800)
var.entropy.meso
dev.off()

png('plots/var.gini.meso.png', width = 800, height = 800)
var.gini.meso
dev.off()
png('plots/spp.cor.meso.png', width = 800, height = 800)
spp.cor.meso
dev.off()

png('plots/var.infogain.hydro.png', width = 800, height = 800)
var.infogain.hydro
dev.off()
png('plots/spp.infogain.hydro.png', width = 800, height = 800)
spp.infogain.hydro
dev.off()
png('plots/var.cor.hydro.png', width = 800, height = 800)
var.cor.hydro
dev.off()
png('plots/spp.cor.hydro.png', width = 800, height = 800)
spp.cor.hydro
dev.off()






# 
# #Information gain
# g1 = 'g1'
# g2 = 'g2'
# v1 = 'Tc'
# 
# x <- pt.df |> as.data.frame()
# x$g1 <- ifelse(x$kc %in% 'cluster2',1,0)
# x$g2 <- ifelse(x$kc %in% 'cluster3',1,0)
# 
# 
# find.threshold <- function(x, class1, class2, variable){
#   x <- as.data.frame(x)
#   g1 <- x[,class1]
#   g2 <- x[,class2]
#   v1 <- x[,variable]
# 
#   x <- data.frame(v1=v1,g1=g1,g2=g2)
#   x <- subset(x, g1 > 0|g2 > 0)
#   rm(v1,g1,g2)
#   vmax <- max(x$v1)
#   vmin <- min(x$v1)
#   class1cor <- cor(x$g1, x$v1)
#   class2cor <- cor(x$g2, x$v1)
#     for(i in 0:100){#i=0
#     ip <- i/100  
#     thr <- ip*(vmax-vmin)+vmin
#     x <- x |> mutate(reclass = ifelse(v1 >= thr,1,0))
#     sign <- var(x$g1,x$reclass)-var(x$g2,x$reclass)
#     xp <- x |> summarise(g1 = sum(g1), g2 = sum(g2), all = g1+g2, p1 = g1/all, p2 = g2/all, 
#                          gini=1-p1^2-p2^2, 
#                          ent = -p1*log(ifelse(p1 <= 0, 0.001,p1))-p2*log(ifelse(p2 <= 0, 0.001,p2)))
#     xs <- x |> group_by(reclass) |> summarise(g1 = sum(g1), g2 = sum(g2), all = g1+g2, p1 = g1/all, p2 = g2/all, 
#     gini=1-p1^2-p2^2, ent = -p1*log(ifelse(p1 <= 0, 1E-50,p1))-p2*log(ifelse(p2 == 0, 1E-50,p2))) |>
#       mutate(gini0 = all/sum(all)*gini, ent0 = all/sum(all)*ent)
#     xs <- xs |> summarize(gini = sum(gini0), ent = sum(ent0))
#     
#     gini_gain = xp$gini - xs$gini
#     ent_gain = xp$ent - xs$ent
#     if(i==0){
#       maxinfo = ent_gain
#       finalentropy = xs$ent
#       finalgini = xs$gini
#       finalsign = sign
#       finalthr = thr
#     }else if(ent_gain >= maxinfo){
#       maxinfo=ent_gain
#       finalentropy = xs$ent
#       finalgini = xs$gini
#       finalsign = sign
#       finalthr = thr}
#   }
#   return(result = list(entropy=finalentropy, gini=finalgini, infogain = maxinfo, sign = ifelse(finalsign >= 0, "class1 positive", "class1 negative"), class1 = class1, class2=class2, class1cor = class1cor, class2cor = class2cor, variable = variable, threshold = finalthr, narrative = paste(
#     class1,": ", variable, ifelse(finalsign >= 0, " >= ", " < "),finalthr," | ",
#     class2,": ", variable, ifelse(finalsign >= 0, " < ", " >= "),finalthr))
#   )}
# 
# pt.df -> x
# x$pyrophytic <- ifelse(x$kc %in% 'cluster2',1,0)
# x$mesophytic <- ifelse(x$kc %in% 'cluster3',1,0)
# thisvar <- find.threshold(x=x, class1='mesophytic',class2='pyrophytic',variable='m')
# thisvar$sign
# thisvar$narrative
# 
# checkvars <- c("p","m","e","s","d","Twh","Tw","Tc","Tcl","Tg",
#                "Tg30","Bhs","carbdepth","clay150","floodfrq","histic",
#                "humic","humicdepth","hydric","ksatdepth","OM150","pH50","rockdepth","sand150",
#                "sand50","spodic","watertable","elev","slope","aspect",
#                "slope500","shade","popen","nopen","solar")
# variables=checkvars
# find.multithreshold <- function(x, class1, class2, variables){
#   n = length(variables)
#   for(i in 1:n){
#     thisvar <- find.threshold(x, class1, class2, variables[i])
#     df0=data.frame(
#       variable = thisvar$variable,
#       threshold = thisvar$threshold,
#       class1 = thisvar$class1,
#       class2 = thisvar$class2,
#       sign = thisvar$sign,
#       infogain = thisvar$infogain,
#       gini = thisvar$gini,
#       entropy = thisvar$entropy,
#       class1cor = thisvar$class1cor,
#       class2cor = thisvar$class2cor,
#       narrative = thisvar$narrative)
#     if(i==1){df=df0}else{df=rbind(df,df0)}
#   }
#   return(df)}
# 
# pt.df -> x
# x$pyrophytic <- ifelse(x$kc %in% 'cluster2',1,0)
# x$mesophytic <- ifelse(x$kc %in% 'cluster3',1,0)
# 
# tableofvars <- find.multithreshold(x=x, class1='mesophytic',class2='pyrophytic',variables=checkvars)



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



#make NMDSclusters on map
pt.df$kc2 <- as.numeric(pt.df$kc)

formular <- as.formula(paste("kc2",paste(names(Species), collapse = " + ", sep = ""), sep = " ~ "))

rf <-  ranger(formular, data = pt.df, classification = T)
named <- 'NMDSClusters.kmeans'
rf.prediction <-  predict(Species, rf, na.rm=T);  names(rf.prediction) <- named
plot(rf.prediction)
writeRaster(rf.prediction, paste0('gis/test',named,'.tif'), overwrite=T)