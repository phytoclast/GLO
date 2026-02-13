#This script to make dendrogram line thickness proportional to a selected species in your data set.

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