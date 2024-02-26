library(stringr)
library(dplyr)
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
wd <- dirname(rstudioapi::getActiveDocumentContext()$path)
sd <- 'UpperMich'
fn <- 'UpperMich.txt'
tif.path <- paste0('gis/dem/tifs/', sd)
filelist <- paste0('gis/dem/',fn)
urls <- read.delim(filelist, header = F)
urls <- subset(urls, grepl('.tif$|.TIF$',V1))
fcount <- str_count(urls$V1, '/')+1
urls$fn <- ''
for(i in 1:nrow(urls)){
urls$fn[i] <-  stringr::str_split_i(urls$V1[i], '/',fcount[i])}

if(!dir.exists(tif.path)){dir.create(tif.path)}
for (i in 1:nrow(urls)){#i=1
 if(!file.exists(paste0(tif.path,'/', urls$fn[i]))){
  download.file(urls$V1[i], paste0(tif.path,'/', urls$fn[i]), method = 'curl')}
}
