#This script to make dendrogram line thickness proportional to a selected species in your data set.----

library(dplyr)
library(vegan)
library(dendextend)
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

#Load species plot matrix (average generic composition per Lower Michigan presettlement cover types)
d0 <- readRDS('d0.RDS')
#Create dissimilarity matrix
rownames(d0) <- d0$COVERTYPE
d <- d0[,-c(1, ncol(d0))]
d <- vegan::vegdist(d)

#Generate dendrogram ----
t <- cluster::agnes(d, method = 'ward')|> as.hclust()

#Rearrange branches of tree
dend <- dendextend::ladderize(as.dendrogram(t), right = T)
# dend <- ape::ladderize(as.phylo(t), right = T) |> as.dendrogram()

#Assign number of clusters ----
k=6
c1 <- cutree(dend, k = k)
# c1 <- pam(d, k = 6, nstart = 20)$clustering
# c1 <- kmeans(d, centers = 6, nstart = 20);c1<-c1$cluster

#Incorporate clusters into dendrogram
c1 <- data.frame(l=names(c1), g=c1)
lab <- data.frame(l=labels(dend)) |> left_join(c1)
#Renumber group numbers so that they will be sequential on the dendrogram
g <- data.frame(g = unique(lab$g)) |> mutate(g2 = 1:k)
lab <- lab |> left_join(g) |> mutate(label = paste(g2,l))
newl <- lab$label
labels(dend) <- newl

#Pick a taxon to vary branch thickness by
sppname <- "Pinus"
Nhem1 <- d0[d0[,sppname] > 1,]$COVERTYPE
Nhem5 <- d0[d0[,sppname] > 5,]$COVERTYPE
Nhem10 <- d0[d0[,sppname] > 10,]$COVERTYPE
Nhem20 <- d0[d0[,sppname] > 20,]$COVERTYPE
Nhem50 <- d0[d0[,sppname] > 50,]$COVERTYPE
Nhem1 <- subset(lab, l %in% Nhem1)$label
Nhem5 <- subset(lab, l %in% Nhem5)$label
Nhem10 <- subset(lab, l %in% Nhem10)$label
Nhem20 <- subset(lab, l %in% Nhem20)$label
Nhem50 <- subset(lab, l %in% Nhem50)$label

#Color branches by cluster
dend <- color_branches(dend, clusters = lab$g2)
dend <- color_labels(dend, col = get_leaves_branches_col(dend))
dend <- dend |> 
  set("branches_lwd", 1) |>
  set("branches_lty", 3) |>
  set("by_labels_branches_lty", value = Nhem1, TF_values = c(1)) |> 
  set("by_labels_branches_lwd", value = Nhem1, TF_values = c(2)) |> 
  set("by_labels_branches_lwd", value = Nhem5, TF_values = c(3))|> 
  set("by_labels_branches_lwd", value = Nhem10, TF_values = c(5))|> 
  set("by_labels_branches_lwd", value = Nhem20, TF_values = c(7))|> 
  set("by_labels_branches_lwd", value = Nhem50, TF_values = c(15)) |> 
  set("labels_cex", 0.7) 
par(mar=c(1,1,1,10))
plot(dend, horiz = TRUE, main = paste("Vegetation of Michigan, abundance of",sppname))

############ NMDS #### ----

library(dplyr)
library(vegan)
library(ggplot2)
#Load Data ----
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
abiotic <- c("p","m","e","s","d","Twh","Tw","Tc","Tcl","Tg","Tg30","Bhs","carbdepth","clay150","floodfrq","histic","humic","humicdepth","hydric","ksatdepth","OM150","pH50","rockdepth","sand150","sand50","spodic","watertable","elev","slope","aspect","slope500","shade","popen","nopen","solar")
taxa <- c("Pinus","Tsuga","Fagus","Acer","Tilia","Carya","Quercus","Ulmus","Fraxinus","Larix","Thuja","Picea","Abies","PopLirio", "Betula",
          "Aesculus","Alnus","Castanea","Celtis","Cercis","Cornus","Crataegus","Cupressaceae","Juglans","Morus","NyssaLiquidambar","OstryaCarpinus","Platanus","Prunus","RobiniaGleditsia","Salix","Sassafras")
mnfi <- readRDS('LowerMichTrees.RDS')

#Establish Species Matrix ----
m <- mnfi[,taxa]
m.env <- mnfi[,abiotic]

#Set number of Dimensions and run Nonmetric MultiDimensional Scaling ordination
ndim <- 3 
set.seed(6190)
nmds <- vegan::metaMDS(m, k=ndim, try=50, trymax = 50)

#join NMDS data components ----
spts.es <- data.frame(plot=rownames(m)) |> cbind(m)  
pt.df <- vegan::scores(nmds, display='sites') |> as_tibble(rownames='sites')  
pt.df <- pt.df |> mutate(vegcode = paste(mnfi$VEGCODE,mnfi$COVERTYPE),vegtype = mnfi$COVERTYPE)
sp.df <- vegan::scores(nmds, display='species') |> as_tibble(rownames='species')
en <- vegan::envfit(nmds, m.env, na.rm = TRUE, choices=c(1:ndim))
en.df <- vegan::scores(en, display='vectors')|> as_tibble(rownames='vectors')

#Generate Convex Hulls ----
fromclusts <- unique(pt.df$vegcode)
for(i in 1:length(fromclusts)){#i=1
  thisclust <-fromclusts[i]
  thishull <- pt.df |> subset(vegcode %in% thisclust)
  chull0 <- chull(thishull$NMDS1, thishull$NMDS2)
  thishull  <- thishull[chull0,]
  if(i==1){veghull=thishull}else{veghull=rbind(veghull,thishull)}
};rm(thisclust,thishull,chull0)

#convex hull for rotated axis ----
if(ndim > 2){
  fromclusts <- unique(pt.df$vegcode)
  for(i in 1:length(fromclusts)){#i=1
    thisclust <-fromclusts[i]
    thishull <- pt.df |> subset(vegcode %in% thisclust)
    chull0 <- chull(thishull$NMDS1, thishull$NMDS3)
    thishull  <- thishull[chull0,]
    if(i==1){veghull2=thishull}else{veghull2=rbind(veghull2,thishull)}
  };rm(thisclust,thishull,chull0)}

#Filter Data to display ----
#these are the vegtypes  
#pt.df$vegcode |> unique()
# vegtypes <- c("BEECH-SUGAR MAPLE FOREST", "BEECH-SUGAR MAPLE-HEMLOCK FOREST", "CEDAR SWAMP", "WHITE PINE-MIXED HARDWOOD FOREST",
#              "WHITE PINE-RED PINE FOREST", "HEMLOCK-WHITE PINE FOREST", "WHITE PINE-WHITE OAK FOREST", "OAK/PINE BARRENS",
#              "MIXED HARDWOOD SWAMP", "MIXED CONIFER SWAMP", "JACK PINE-RED PINE FOREST", "SPRUCE-FIR-CEDAR FOREST",
#              "PINE BARRENS", "ASPEN-BIRCH FOREST", "MIXED PINE-OAK FOREST", "GRASSLAND",
#              "BLACK OAK BARREN","MIXED OAK SAVANNA","SHRUB SWAMP/EMERGENT MARSH", "MIXED OAK FOREST",
#              "OAK-HICKORY FOREST", "WET PRAIRIE", "BLACK ASH SWAMP")


#filter vegtypes to display 
thisveg <- c("OAK-HICKORY FOREST","BLACK OAK BARREN","MIXED OAK SAVANNA","BEECH-SUGAR MAPLE FOREST","JACK PINE-RED PINE FOREST", "WHITE PINE-WHITE OAK FOREST", "WHITE PINE-RED PINE FOREST", "HEMLOCK-WHITE PINE FOREST","BEECH-SUGAR MAPLE-HEMLOCK FOREST")
#filter environmental vectors to display 
thisabiotic <- c("Tg","hydric","pH50","m","clay150","sand150","watertable","slope")

selected.pt.df <- pt.df |> subset(vegtype %in% thisveg)
selected.veghull <- veghull |> subset(vegtype %in% thisveg)
selected.veghull2 <- veghull2 |> subset(vegtype %in% thisveg)
selected.en.df <- en.df |> subset(vectors %in% thisabiotic)


#NMDS Biplots ----
#plot axes 1 and 2
gp <- ggplot() +
  geom_polygon(data=selected.veghull, aes(x=NMDS1,y=NMDS2, color=vegcode, fill=vegcode), alpha=0.1, linewidth=1)+
  geom_point(data=selected.pt.df, aes(x=NMDS1,y=NMDS2, color=vegcode), alpha=0.5, size=2)+
  geom_point(data=sp.df, aes(x=NMDS1,y=NMDS2), color='black')+
  geom_text(data=sp.df, aes(label=species, x=NMDS1,y=NMDS2), vjust = 0, nudge_y = 0.02, nudge_x = 0.05, color='black', size=2)+
  geom_segment(data=selected.en.df, aes(x=0,y=0,xend=NMDS1,yend=NMDS2), arrow = arrow(length = unit(0.03, "npc")), color='darkgray')+
  geom_text(data=selected.en.df, aes(label=vectors, x=NMDS1,y=NMDS2), vjust = 0, nudge_y = 0.02, nudge_x = 0.05, color='darkgray')+
  coord_fixed()
gp

#plot axes 1 and 3
if(ndim > 2){
gp2 <- ggplot() +
  geom_polygon(data=selected.veghull2, aes(x=NMDS1,y=NMDS3, color=vegcode, fill=vegcode), alpha=0.1, linewidth=1)+
  geom_point(data=selected.pt.df, aes(x=NMDS1,y=NMDS3, color=vegcode), alpha=0.5, size=2)+
  geom_point(data=sp.df, aes(x=NMDS1,y=NMDS3), color='black')+
  geom_text(data=sp.df, aes(label=species, x=NMDS1,y=NMDS3), vjust = 0, nudge_y = 0.02, nudge_x = 0.05, color='black', size=2)+
  geom_segment(data=selected.en.df, aes(x=0,y=0,xend=NMDS1,yend=NMDS3), arrow = arrow(length = unit(0.03, "npc")), color='darkgray')+
  geom_text(data=selected.en.df, aes(label=vectors, x=NMDS1,y=NMDS3), vjust = 0, nudge_y = 0.02, nudge_x = 0.05, color='darkgray')+
  coord_fixed()
gp2
}

