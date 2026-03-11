##Sediment Interpolation Script
## Data from BNDO

#Loading packages-------

packages<-c('sf','ggplot2','terra','dplyr','Rsagacmd')

package.check <- lapply(packages, FUN = function(x) {
  if (!require(x, character.only = TRUE)) {
    install.packages(x, dependencies = TRUE)
    library(x, character.only = TRUE)
  }
})
saga<-saga_gis()

#Creating Projection -------------
#Using Albers Brazil (EPSG:10857)
crs_albers_brasil <- "+proj=aea +lat_0=-12 +lon_0=-54 +lat_1=-2 +lat_2=-22 +x_0=5000000 +y_0=10000000 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs"

#Reading and preprocessing data-----------

##Points---------
pts<-read_sf('input_data/bndo_eez/bndo_eez.shp')

pts<-st_transform(pts,crs=crs_albers_brasil)

pts<-pts%>%transmute(id=num_amostr,
                     lon=lon,
                     lat=lat,
                     time=as.POSIXct(data_hora,tz="America/Sao_Paulo",format='%Y/%m/%d %H:%M:%S'),
                     sand=areia,
                     gravel=cascalho,
                     mud=argila+silte,
                     tot=areia+cascalho+argila+silte)

pts<-pts%>%filter(tot==100)


n_acc<-round(nrow(pts)*0.05)
pts<-pts%>%mutate(use='interp')

index<-sample(seq_len(nrow(pts)),size=n_acc)
pts$use[index]<-'accuracy'

write_sf(pts,'input_data/bndo.shp')

pts_interp<-pts%>%filter(use=='interp')
pts_acc<-pts%>%filter(use=='accuracy')

##Mask (EEZ)----
eez<-read_sf('input_data/study_area.shp')
eez<-st_transform(eez,crs_albers_brasil)
mask<-rast(eez, res=1000)
mask<-rasterize(eez,mask,field='id')
mask<-project(mask,crs_albers_brasil,res=1000)

writeRaster(mask,'input_data/eez_albers.tif',overwrite=T)

##Bathymetric model
bat<-rast('input_data/batimetria_dhn.tif')
bat<-project(bat,mask)
bat<-mask(bat,mask)
res(bat)
plot(bat)

writeRaster(bat,'input_data/batimetria_dhn_albers.tif',overwrite=T)

mask_bat<-bat>=-1000
plot(mask_bat)
mask_bat<-subst(mask_bat, from=0, to=NA)
writeRaster(mask_bat,'input_data/mask_bat.tif',overwrite=T)

#Interpolating----------
##Sand----------

sand <- saga$grid_gridding$natural_neighbour(
  points = pts_interp,
  field = "sand",
  target_definition = 1,
  target_template = mask
)

sand<-mask(sand,mask_bat)
plot(sand)

writeRaster(sand,'output_data/nn_sand.tif',overwrite=T)


##Mud----------

mud <- saga$grid_gridding$natural_neighbour(
  points = pts_interp,         
  field = "mud",
  target_definition = 1,
  target_template = mask
)


mud<-mask(mud,mask_bat)
plot(mud)

writeRaster(mud,'output_data/nn_mud.tif',overwrite=T)

##Gravel----------

gravel <- saga$grid_gridding$natural_neighbour(
  points = pts_interp,         
  field = "gravel",
  target_definition = 1,
  target_template = mask
)

gravel<-mask(gravel,mask_bat)
plot(gravel)

writeRaster(gravel,'output_data/nn_gravel.tif',overwrite=T)


#Accuracy assessment-----------
acc_sand<-extract(sand,pts_acc)
names(acc_sand)<-c('id','sand_interp')
acc_mud<-extract(mud,pts_acc)
names(acc_mud)<-c('id','mud_interp')
acc_gravel<-extract(gravel,pts_acc)
names(acc_gravel)<-c('id','gravel_interp')


pts_acc<-pts_acc%>%mutate(sand_interp=acc_sand$sand_interp,
                  mud_interp=acc_mud$mud_interp,
                  gravel_interp=acc_gravel$gravel_interp)

ggplot(pts_acc,aes(x=sand-sand_interp))+geom_histogram()
ggplot(pts_acc,aes(x=mud-mud_interp))+geom_histogram()
ggplot(pts_acc,aes(x=gravel-gravel_interp))+geom_histogram()


accuracy<-pts_acc%>%
  st_drop_geometry()%>%
  filter(!is.na(sand_interp))%>%
  summarise(sand_me=mean(sand-sand_interp),
            sand_mae=mean(abs(sand-sand_interp)),
            sand_rmse=sqrt(mean((sand-sand_interp)^2)),
            mud_me=mean(mud-mud_interp),
            mud_mae=mean(abs(mud-mud_interp)),
            mud_rmse=sqrt(mean((mud-mud_interp)^2)),
            gravel_me=mean(gravel-gravel_interp),
            gravel_mae=mean(abs(gravel-gravel_interp)),
            gravel_rmse=sqrt(mean((gravel-gravel_interp)^2)))
accuracy

sed<-c(sand,mud,gravel)

names(sed)<-c('sand','mud','gravel')

t<-sum(sed)
hist(t)
plot(t)


writeRaster(sed,'output_data/sediment.tif',overwrite=T)
