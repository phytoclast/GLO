setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
library(terra)
library(sf)
library(dplyr)

rackCorners <- function(pts, dataset, treenames, diams, x, y){
  pts$id <- rownames(pts)
  ntrees = length(treenames)
  for(i in 1:ntrees){#i=1
    pts0 = data.frame(
      dataset = dataset,
      id = pts$id,
      treeseq = i,
      tname = pts[,treenames[i]],
      diam = pts[,diams[i]],
      x=pts[,x],
      y=pts[,y]
    )
    if(i==1){pts1=pts0}else{pts1=rbind(pts1,pts0)}
  }
  return(pts1)
}

pts <- read.csv('data/msb-paleon.32.0/PLS_Wisconsin_trees_Level0_v1.0.csv')
treenames <- c("SP1","SP2","SP3","SP4")
diams <- c("DIAM1","DIAM2","DIAM3","DIAM4")
dataset='WI'
x = 'x_alb'
y = 'y_alb'
pts.WI <- rackCorners(pts, dataset, treenames, diams, x, y)

pts <- read.csv('data/msb-paleon.27.0/PLS_Indiana_trees_Level0_v1.0.csv')
treenames <- c("L3_tree1","L3_tree2","L3_tree3","L3_tree4")
diams <- c("diameter","diameter2","diameter3","diameter4")
dataset='IN'
x = 'x'
y = 'y'
pts.IN <- rackCorners(pts, dataset, treenames, diams, x, y)

pts <- read.csv('data/msb-paleon.28.0/PLS_Illinois_trees_Level0_v1.0.csv')
treenames <- c("L3_tree1","L3_tree2","L3_tree3","L3_tree4")
diams <- c("diameter","diameter2","diameter3","diameter4")
dataset='IL'
x = 'x'
y = 'y'
pts.IL <- rackCorners(pts, dataset, treenames, diams, x, y)

pts <- read.csv('data/msb-paleon.29.0/PLS_SoutheasternMichigan_trees_Level0_v1.0.csv')
treenames <- c("L3_tree1","L3_tree2","L3_tree3","L3_tree4")
diams <- c("diameter","diameter2","diameter3","diameter4")
dataset='SEMI'
x = 'x'
y = 'y'
pts.SEMI <- rackCorners(pts, dataset, treenames, diams, x, y)

pts <- read.csv('data/msb-paleon.30.0/PLS_southernMichigan_trees_Level0_v1.0.csv')
treenames <- c("species1","species2","species3","species4")
diams <- c("diam1","diam2","diam3","diam4")
dataset='SMI'
x = 'point_x'
y = 'point_y'
pts.SMI <- rackCorners(pts, dataset, treenames, diams, x, y)

pts <- read.csv('data/msb-paleon.31.0/PLS_northernMichigan_trees_Level0_v1.0.csv')
treenames <- c("SPP1","SPP2","SPP3","SPP4")
diams <- c("DBH1","DBH2","DBH3","DBH4")
dataset='NMI'
x = 'x_alb'
y = 'y_alb'
pts.NMI <- rackCorners(pts, dataset, treenames, diams, x, y)

#pts <- read.csv('data/msb-paleon.2.0/SetTreeComp_West_Level1_v1.0.csv')

pts <- rbind(pts.WI,pts.IL,pts.IN,pts.NMI,pts.SMI,pts.SEMI)
pts.geo <- st_as_sf(pts, coords=c(x='x', y='y'), crs=3175)
write_sf(pts.geo,'gis/trees.shp')