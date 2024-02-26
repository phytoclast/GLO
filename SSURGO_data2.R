library(soilDB)
library(aqp)
library(sf)
library(mapview)
library(vegnasis)
#set working directory to folder where this R file is saved
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

#identify location of geodatabase
states <- c("MI", "IL", "IN", "OH")
for (i in 1:length(states)){#i=1
mi0 = fetchGDB(dsn = paste0('D:/GIS/SOIL/2023/gSSURGO_',states[i],'.gdb'))
mu0 = get_mapunit_from_GDB(dsn = paste0('D:/GIS/SOIL/2023/gSSURGO_',states[i],'.gdb'))
co0 = get_component_from_GDB(dsn = paste0('D:/GIS/SOIL/2023/gSSURGO_',states[i],'.gdb'))
comonth0 <- sf::st_read(paste0('D:/GIS/SOIL/2023/gSSURGO_',states[i],'.gdb'), 'comonth')
if(i==1){
  mi=mi0;mu=mu0;co=co0;comonth1=comonth0
}else{
  mi=c(mi,mi0);mu=rbind(mu,mu0);co=rbind(co,co0);comonth1=rbind(comonth1,comonth0)
}
};rm(mi0,mu0,co0,comonth0)
#extract and concatenate component table
s <-  site(mi)
h <- horizons(mi)
comonth <- comonth1 |> mutate(flood = ifelse(flodfreqcl %in% c("Frequent","Occasional","Rare","Very rare"), 1,0)) |> group_by(cokey) |> summarise(floodfrq=max(flood)) 


h <- h |> mutate(h50 = ifelse(hzdepb_r > 50, 50, hzdepb_r) - ifelse(hzdept_r > 50, 50, hzdept_r),
                 h150 = ifelse(hzdepb_r > 150, 150, hzdepb_r) - ifelse(hzdept_r > 150, 150, hzdept_r),
                 pH = ifelse(is.na(ph1to1h2o_r), ph01mcacl2_r, ph1to1h2o_r),
                 Bh = ifelse(grepl('Bh', hzname),1, 0),
                 rock = ifelse(!is.na(hzname) & grepl('Cr|R', hzname),hzdept_r, 250),
                 carb = ifelse(!is.na(caco3_r) & caco3_r > 0, hzdept_r, 250),
                 ahor = case_when(is.na(om_r) ~ 0, om_r >= 2.5 ~ pmin(200, hzdepb_r),
                                  hzdept_r < 1 ~ pmin(hzdepb_r, hzdepb_r*om_r/2),
                                  TRUE ~ 0),
                 ksr = case_when(is.na(ksat_r) ~ 250, 
                                 ksat_r < 3.5 ~ hzdept_r,
                                 TRUE  ~ pmin(250, hzdept_r+150*((ksat_r-3.5)/3.5))))
                 
# h <- h |> mutate(ksat2 = 2^(7.7322+0.00059547*sandtotal_r^2+-4.3277*dbthirdbar_r+-0.03109*om_r))#formula developed before 2017 for bad ksat values
# ksats <- subset(h, !is.na(ksat_r) & !is.na(ksat2), select=c(cokey, hzname, hzdept_r, hzdepb_r, ksat_r,ksat2, ksr, om_r, ahor))


soil <- h |> group_by(cokey) |> summarise(sand50 = weighted.mean(sandtotal_r, h50, na.rm = T),
                                       sand150 = weighted.mean(sandtotal_r, h150, na.rm = T),
                                       pH50 = weighted.mean(pH, h50, na.rm = T),
                                       OM150 = weighted.mean(om_r, h150, na.rm = T),
                                       rockdepth = min(rock, na.rm = T),
                                       carbdepth = min(carb, na.rm = T),
                                       humicdepth = max(ahor, na.rm = T),
                                       ksatdepth = min(ksr, na.rm = T), #depth to restricted saturated conductivity
                                       Bh = max(Bh, na.rm = T))
soil <- subset(mu, select=c(muname, mukey)) |> left_join(subset(s, select=c(mukey, cokey, comppct_r, majcompflag, compname, taxclname, taxorder, taxsubgrp, hydricrating, drainagecl, geomdesc))) |> left_join(soil, by=join_by(cokey==cokey)) |> left_join(comonth, by=join_by(cokey==cokey)) 


soil <- soil |> mutate(
  moist = case_when(hydricrating %in% 'yes'~'wet',
                    drainagecl %in% c("Somewhat poorly drained","Very poorly drained","Moderately well drained","Poorly drained") ~ 'moist',
                    TRUE ~ 'dry'),
  rockdepth = ifelse(rockdepth > 50 & grepl('lithic',tolower(taxsubgrp)), 25, rockdepth),
  
  chem = case_when(carbdepth < 100 ~ 'calcareous',
                   grepl('ult|ods',tolower(taxsubgrp)) | grepl('dys', tolower(taxclname)) ~ 'dysic',
                   grepl('alf|oll',tolower(taxsubgrp)) ~ 'euic',
                   pH50 < 5.5 ~ 'dysic',
                   TRUE ~ 'euic'),
  text = case_when(rockdepth < 100 ~ 'rocky',
                   sand50 >= 70 & sand150 >= 80 | sand50 >= 80 ~ 'sandy',
                   grepl('ist', tolower(taxsubgrp)) | grepl('ist', tolower(taxorder)) |  OM150 > 20 ~ 'mucky',
                   TRUE ~ 'loamy'),
  flood = case_when(grepl('flood', muname) | (grepl('flood', geomdesc) & grepl('fluv', tolower(taxsubgrp))) | floodfrq > 0 ~ 'flood',
                    TRUE ~ ''),
  watertable = case_when(drainagecl %in% c('Well drained', 'Somewhat excessively drained', 'Excessively drained') ~ 250,
                         drainagecl %in% c('Moderately well drained') ~ 100,
                         drainagecl %in% c('Somewhat poorly drained') ~ 50,
                         drainagecl %in% c('Poorly drained') ~ 25,
                         drainagecl %in% c('Very poorly drained') ~ 0,
                         TRUE ~ NA),
  hydric = ifelse(hydricrating %in% 'Yes',1,0))


#summarize mapunits by major component weighted by composition
forraster <- soil |> subset(majcompflag %in% 'Yes') |> 
  group_by(mukey) |> 
  summarise(sand50 = weighted.mean(sand50, comppct_r, na.rm=T),
            sand150 = weighted.mean(sand150, comppct_r, na.rm=T),
            pH50 = weighted.mean(pH50, comppct_r, na.rm=T),
            OM150 = weighted.mean(OM150, comppct_r, na.rm=T),
            hydric = weighted.mean(hydric, comppct_r, na.rm=T),
            watertable = weighted.mean(watertable, comppct_r, na.rm=T),
            carbdepth = weighted.mean(carbdepth, comppct_r, na.rm=T),
            rockdepth = weighted.mean(rockdepth, comppct_r, na.rm=T),
            ksatdepth = weighted.mean(ksatdepth, comppct_r, na.rm=T),
            humicdepth = weighted.mean(humicdepth, comppct_r, na.rm=T),
            floodfrq = weighted.mean(floodfrq, comppct_r, na.rm=T)) |> 
  subset(select=c(mukey, sand50, sand150, pH50, OM150, watertable, hydric, rockdepth, ksatdepth, humicdepth, carbdepth, floodfrq))

forraster$mukey <- as.numeric(forraster$mukey)

write.csv(forraster, 'forraster.csv', na='', row.names = F)





library(terra)
mirast <- rast('D:/GIS/SOIL/2023/MI.Soil.30.tif')
mirast <- rast('D:/GIS/SOIL/2023/Soil_30.tif')
milev <- as.factor(mirast)
levs <- levels(milev)
df <- levs[[1]] |> left_join(forraster, by=join_by(ID==mukey), multiple='first') |> subset(select=c(ID, sand50, sand150, pH50, watertable, hydric, OM150, rockdepth, ksatdepth, humicdepth, carbdepth, floodfrq))
levels(milev) <- df
plot(milev)
writeRaster(milev, 'milev.tif', overwrite=T)
