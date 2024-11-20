library(climatools)
library(soilDB)
library(aqp)
library(sf)
library(mapview)
library(vegnasis)
library(terra)
#set working directory to folder where this R file is saved
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
prism <- rast('gis/climate/prism.tif')
prism2 <- rast('gis/climate/prism2.tif')
pet <- rast('gis/climate/pet.tif')
prism <- c(prism, pet)
Tclx <- climatools::XtremLow(prism2$Tcl, prism$lat, prism$lon, prism$elev); names(Tclx) <- 'Tclx'
writeRaster(Tclx, 'gis/climate/Tclx.tif', overwrite=T)
a <- AET.rast.monthly(block=prism, jan.p='p01', jan.e='e01')
max3aet <- AET.rast.max(block=a, jan.a='a01', nmonth = 3)
writeRaster(max3aet, 'gis/climate/max3aet.tif', overwrite=T)
prism2 <- c(prism2,Tclx,max3aet)
BPS <- rast("C:/GIS/Vegetation/BPS_general.tif")
BPSrat <- cats(BPS)[[1]]
prism3 <- c(prism2,prism)
BPSsample <- spatSample(BPS,250000,xy=T)
BPSsample <- BPSsample |> subset(!is.na(VALUE_1) & !VALUE_1 %in% -9999) |> mutate(VALUE = as.numeric(as.character(VALUE_1)), VALUE_1 = NULL) |> left_join(BPSrat, by = join_by(VALUE == VALUE))
BPSsample <- st_as_sf(BPSsample, coords = c(x='x',y='y'), crs=crs(BPS))
BPSsample <- st_transform(BPSsample, crs = crs(prism))

BPSprism <- extract(prism3, BPSsample)

BPSsample <- cbind(BPSsample,BPSprism)
BPSsample <- BPSsample |> mutate(M = p/(e+0.0001))
BPSsummary <- BPSsample |> subset(!is.na(Tg) & !is.na(m)) |> group_by(BPS_NAME) |> summarise(lat=median(lat), lon=median(lon), elev=median(elev),  Tg.90th = quantile(Tg,0.9), Tg.10th = quantile(Tg,0.1), Tg = median(Tg),Tw=median(Tw), Tc=median(Tc), Tclx=median(Tclx), p=median(p), e=median(e), M.90 = quantile(M,0.9), M.10 = quantile(M,0.1), M=median(M), m=median(m), md = 1-median(d), max3aet = median(max3aet), n=length(BPS_NAME)) |> st_drop_geometry()
write.csv(BPSsummary,'bpsclimate.csv', row.names = F)

#https://chelsa-climate.org/downloads/
month <- c('01','02','03','04','05','06','07','08','09','10','11','12')
dem1km <- rast("C:/GIS/Climate/dem_latlong.nc")
tmax <- rast("C:/GIS/Climate/CHELSA_tasmax_1981-2010_V.2.1.nc")
tmin <- rast("C:/GIS/Climate/CHELSA_tasmin_1981-2010_V.2.1.nc")
ppt <- rast("C:/GIS/Climate/CHELSA_pr_1981-2010_V.2.1.nc")
th <- tmax*10
tl <- tmin*10
p <- ppt*100
dem <- project(dem1km, th)
names(th) <- paste0('th',month)
names(tl) <- paste0('tl',month)
names(p) <- paste0('p',month)
names(dem) <- 'dem'

for (i in 1:12){
  writeRaster(th[[i]], paste0('chelsa2.1/th',month[i],'.tif'))
}
for (i in 1:12){
  writeRaster(tl[[i]], paste0('chelsa2.1/tl',month[i],'.tif'))
}
for (i in 1:12){
  writeRaster(p[[i]], paste0('chelsa2.1/p',month[i],'.tif'))
}
writeRaster(dem, paste0('chelsa2.1/dem.tif'))
writeRaster(dem1km, paste0('chelsa2.1/dem1km.tif'))

#Load Chelsa Files
library(climatools)
library(soilDB)
library(aqp)
library(sf)
library(mapview)
library(vegnasis)
library(terra)
#set working directory to folder where this R file is saved
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
month <- c('01','02','03','04','05','06','07','08','09','10','11','12')
for (i in 1:12){#i=1
  x = rast(paste0('chelsa2.1/th',month[i],'.tif'))

  assign(paste0('th',month[i]), x)
}
for (i in 1:12){
  x = rast(paste0('chelsa2.1/tl',month[i],'.tif'))

  assign(paste0('tl',month[i]),x)
}
for (i in 1:12){
  x = rast(paste0('chelsa2.1/p',month[i],'.tif'))

  assign(paste0('p',month[i]),x)
}
dem <- rast('chelsa2.1/dem.tif')
latlon <- as.data.frame(th01, xy=T); colnames(latlon) <- c('x','y','z')
latlon[,3] <- latlon[,2];
lat = rast(x=latlon, type="xyz", crs=crs(th01)); names(lat) <- 'lat'
latlon[,3] <- latlon[,1];
lon = rast(x=latlon, type="xyz", crs=crs(th01));  names(lon) <- 'lon'

lat <- project(lat, th01);lon <- project(lon, th01)
chelsa <- c(lat,lon, dem)
chelsa <- c(chelsa, rast(mget(paste0('th',month))))
chelsa <- c(chelsa, rast(mget(paste0('tl',month))))
chelsa <- c(chelsa, rast(mget(paste0('p',month))))
writeRaster(chelsa, 'chelsa2.1/chelsa.tif', overwrite=T)
chelsa <- rast('chelsa2.1/chelsa.tif')
timeA <-  Sys.time()
e01 <- GetPET.block(1,chelsa)
e02 <- GetPET.block(2,chelsa)
e03 <- GetPET.block(3,chelsa)

e04 <- GetPET.block(4,chelsa)
e05 <- GetPET.block(5,chelsa)
e06 <- GetPET.block(6,chelsa)

e07 <- GetPET.block(7,chelsa)
e08 <- GetPET.block(8,chelsa)
e09 <- GetPET.block(9,chelsa)

e10 <- GetPET.block(10,chelsa)
e11 <- GetPET.block(11,chelsa)
e12 <- GetPET.block(12,chelsa)
Sys.time() - timeA
pet <- c(e01,e02,e03,e04,e05,e06,e07,e08,e09,e10,e11,e12)
writeRaster(pet, 'chelsa2.1/pet.tif', overwrite=T)

#----
library(climatools)
library(soilDB)
library(aqp)
library(sf)
library(mapview)
library(vegnasis)
library(terra)
#set working directory to folder where this R file is saved
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
pet <- rast('chelsa2.1/pet.tif')
chelsa <- rast('chelsa2.1/chelsa.tif')
dem1km <- rast('chelsa2.1/dem1km.tif')
chelsa <- c(chelsa, pet)
dem <- chelsa$dem
a <- climatools::AET.rast(chelsa)
plot(a)

p <- ApplyClim.rast(chelsa, name='p', jan='p01', fun='sum')
e <- ApplyClim.rast(chelsa, name='e', jan='e01', fun='sum')
d <- (e-a)/(e+0.0001); names(d) <- 'd'
s <- (p-a)/(p+0.0001); names(s) <- 's'
plot(s)

tblock <- climatools::meanT.rast(chelsa)

Twh <- ApplyClim.rast(chelsa, jan='th01', fun='max', name='Twh')
Tw <- ApplyClim.rast(tblock, jan='t01', fun='max', name='Tw')
Tcl <- ApplyClim.rast(chelsa, jan='tl01', fun='min', name='Tcl')
Tc <- ApplyClim.rast(tblock, jan='t01', fun='min', name='Tc')
a <- AET.rast.monthly(block=chelsa, jan.p='p01', jan.e='e01')
max3aet <- AET.rast.max(block=a, jan.a='a01', nmonth = 3)
plot(max3aet)
btblock <- ifel(tblock > 0, tblock, 0)

Tg1 <- ApplyClim.rast(btblock, jan='t01', mons=c(12,1,2,3,4), fun='mean')
Tg2 <- ApplyClim.rast(btblock, jan='t01', mons=c(5,6,7,8,9,10), fun='mean')
Tg <- max(Tg1, Tg2);names(Tg) <- 'Tg'
m <- p/(p+e+0.0001); names(m) <- 'm'

Tclx <- XtremLow(Tcl,Lat=chelsa$lat,Lon=chelsa$lon,Elev = chelsa$dem); names(Tclx) <- 'Tclx'
plot(Tclx)
chelsa2 <- c(p,m,e,s,d,max3aet,Twh,Tw,Tc,Tcl,Tclx,Tg)
writeRaster(chelsa2, 'chelsa2.1/chelsa2a.tif', overwrite=T)
chelsa2 <- rast('chelsa2.1/chelsa2.tif')
dem.amp <- climatools::RestoreMaxMin(dem1km, dem)
Tg.amp <- climatools::enhanceRast(Tg, dem, dem.amp)
writeRaster(Tg.amp, 'chelsa2.1/Tg.tif')
writeRaster(dem.amp, 'chelsa2.1/dem.amp.tif', overwrite=T)
