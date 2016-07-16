# remove all objects
# --------------------------------------------------
# rm(list=ls(all=TRUE))
# --------------------------------------------------
# --------------------------------------------------

# -----------------------------------------------------------------
# -----------------------------------------------------------------
# -----------------------------------------------------------------
# lapse rate calculation function
lasperatefunction=function(stationfile)
{
library(lubridate)
#----------------------------------------------------------
readfolder=unlist(strsplit(stationfile, "\\\\"))
setwd(paste(file.path(as.vector(unlist(strsplit(stationfile, "\\\\"))
[1:length(unlist(strsplit(stationfile, "\\\\")))-1])),collapse="\\\\"))

stations=read.csv(stationfile,sep=",",colClasses = "character")
for (i in 1:dim(stations)[1])
{
if (!file.exists(paste(stations[i,1],".csv",sep="")))
{
stop("File ", paste(stations[i,1],".csv",sep="")," doesn't exist")
}
}
temparr=array(NA,c(12,2,dim(stations)[1]))
for (i in 1:dim(stations)[1])
{
data=read.csv(paste(stations[i,1],".csv",sep=""))
data=na.omit(data)
agg=aggregate(data$Temp,by=list(month(as.Date(data$Date,format="%m/%d/%Y"))),FUN=mean)
temparr[,1,i]=stations[i,2]
temparr[,2,i]=agg$x
}
lapserate=vector("numeric",length=12)
for (i in 1:dim(temparr)[1])
{
eachmon=data.frame(as.numeric(temparr[i,1,]),as.numeric(temparr[i,2,]))
lm.lapserate <- lm(eachmon[,2] ~ eachmon[,1])
lapserate[i]= lm.lapserate$coefficients[2]
}
return(list(lapserate))
}
# -----------------------------------------------------------------
# -----------------------------------------------------------------
# -----------------------------------------------------------------


# get time series at a point given by lat-lon file coordinate 
# and folder location and pattern for filenames
# -----------------------------------------------------------------
# -----------------------------------------------------------------
# -----------------------------------------------------------------
# find a ROD of a single grid cell
getpointROD=function(filelocation,pattern="merra.rfe.90m.",xindex,yindex,varname)
{
library(ncdf)
filecollection=list.files(filelocation,pattern=pattern)
nooffiles=length(filecollection)
lonpos=xindex
latpos=yindex
start=c(1,latpos,lonpos)
count=c(-1,1,1)

timevar=vector('numeric',length=nooffiles)
timeunits=vector('character',length=length(nooffiles))

nooffiles=length(filecollection)
for (i in 1:nooffiles)
{
file=file.path(filelocation,filecollection[i])
nc=open.ncdf(file)
timeunits[i]=att.get.ncdf(nc,'time','units')$value
timevar[i]=length(get.var.ncdf(nc,'time'))
close.ncdf(nc)
}

tnu=unique(timeunits)
if(length(tnu)>1)
{
print("NetCDF error: there is more than one time unit")
}

ttt=strsplit(tnu," ")[[1]][3]
stdate=as.Date(strsplit(ttt,"T")[[1]][1])
timetol=sum(timevar)
dates=vector('numeric',length=sum(timevar))
Value=vector('numeric',length=sum(timevar))
datesx=vector('numeric',length=length(timevar))
datesy=vector('numeric',length=length(timevar))
for (i in 1:nooffiles)
{
if(i==1)
{
datenow=1
dateend=timevar[i]
file=file.path(filelocation,filecollection[i])
nc=open.ncdf(file)
dates[datenow:dateend]=get.var.ncdf(nc,'time')
Value[datenow:dateend]=get.var.ncdf(nc,varname,start=start,count=count)
close.ncdf(nc)
dateend=timevar[i]
datesx[i]=datenow
datesy[i]=dateend
}else
{
datenow=datenow+timevar[i-1]
dateend=dateend+timevar[i]
file=file.path(filelocation,filecollection[i])
nc=open.ncdf(file)
Value[datenow:dateend]=get.var.ncdf(nc,varname,start=start,count=count)
dates[datenow:dateend]=get.var.ncdf(nc,'time')
close.ncdf(nc)
datesx[i]=datenow
datesy[i]=dateend
}
}
Date=stdate+hours(dates)
return(vartimeseries=data.frame(Date,Value))
}
# -----------------------------------------------------------------
# -----------------------------------------------------------------
# -----------------------------------------------------------------


# -----------------------------------------------------------------
# -----------------------------------------------------------------
# -----------------------------------------------------------------
# lat-lon position of a station on a DEM

latlonpos=function(demfile,stationfile)
{
library(raster)
library(rgdal)
library(maptools)

dem<-raster(readGDAL(demfile))
demproj=projection(dem)
res=res(dem)[1]
xs=seq(dem@extent@xmin,dem@extent@xmax,by=res)
ys=seq(dem@extent@ymin,dem@extent@ymax,by=res)

gaugelatlon=read.csv(stationfile)  
crs = "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0"
midpnt=SpatialPoints(coords=matrix(c(gaugelatlon[,4],gaugelatlon[,3]),nrow=dim(gaugelatlon)[1]),proj4string=CRS(crs))
midlatlon=spTransform(midpnt,CRS(demproj))

xindex=vector("numeric",dim(gaugelatlon)[2])
yindex=vector("numeric",dim(gaugelatlon)[2])

for (j in 1:length(midlatlon@coords[,1]))
{
xindex[j] = sum(xs < midlatlon@coords[j,1]) + 1
}
for (j in 1:length(midlatlon@coords[,2]))
{
yindex[j] = sum(ys < midlatlon@coords[j,2]) + 1
}
return(list(xindex,yindex))
}
# -----------------------------------------------------------------
# -----------------------------------------------------------------
# -----------------------------------------------------------------
# precipitation adjustment factor

precipitationadjfactorfunction=function(stationfile2)
{
library(lubridate)
readfolder=unlist(strsplit(stationfile2, "\\\\"))
setwd(paste(file.path(as.vector(unlist(strsplit(stationfile2, "\\\\"))
[1:length(unlist(strsplit(stationfile2, "\\\\")))-1])),collapse="\\\\"))
stations=read.csv(stationfile2,sep=",",colClasses = "character")
for (i in 1:dim(stations)[1])
{
if (!file.exists(paste("prec",stations[i,1],".csv",sep="")))
{
stop("File ", paste("prec",stations[i,1],".csv",sep="")," doesn't exist")
}
}
temparr=array(NA,c(12,2,dim(stations)[1]))
k=vector("numeric",12)
for (i in 1:dim(stations)[1])
{
data=read.csv(paste("Prec",stations[i,1],".csv",sep=""),stringsAsFactors =FALSE)
data=na.omit(data)
agg=aggregate(data$Prec,by=list(month(as.Date(data$Date,format="%m/%d/%Y"))),FUN=mean, na.rm=T)
temparr[,1,i]=stations[i,2]
temparr[,2,i]=agg$x*0.0254
}
K=vector("numeric",length=12)
A=vector("numeric",length=12)
for (i in 1:dim(temparr)[1])
{
eachmon=data.frame(as.numeric(temparr[i,1,]),as.numeric(temparr[i,2,]))
colnames(eachmon)=c("elevation","p")
p0=min(eachmon[,2],na.rm=TRUE)
if(p0==0){p0=0.001}
z0=min(eachmon[,1],na.rm=TRUE)
mmModel <- nls(p ~ A + p0*(1+K*(elevation-z0))/(1-K*(elevation-z0)),data=eachmon,start=list(A=0.0004, K=0.001))
K=coef(mmModel)
x <- seq(min(eachmon[,1]), max(eachmon[,1]), length=100)
z <- predict(mmModel, list(elevation=x))
#plot(eachmon[,1],eachmon[,2])
#points(x,z,type="l",col="blue")
k[i]=as.numeric(K)[2]
}
return(list(k))
}
# -----------------------------------------------------------------
# -----------------------------------------------------------------
# -----------------------------------------------------------------
# read the temperature data that from gridcell (81,203)
getval=function(workdir,xpos,ypos,varname)
{
library(raster)
library(ncdf)
files<-list.files(workdir,pattern="merra.rfe.")
nooffiles<-length(files)
timevar=vector('numeric',length=nooffiles)
for (i in 1:nooffiles)
{
file=files[i]
setwd(workdir)
nc=open.ncdf(file)
timevar[i]=length(get.var.ncdf(nc,'time'))
close.ncdf(nc)
}
timseries=vector('numeric',length=sum(timevar))
tval=vector('numeric',length=sum(timevar))
timeunits=vector('character',length=nooffiles)

for (i in 1:nooffiles)
{
file=files[i]
nc=open.ncdf(file)
timeunits[i]=att.get.ncdf(nc,'time','units')$value
close.ncdf(nc)
}
timetol=sum(timevar)
datesx=vector('numeric',length=length(timevar))
datesy=vector('numeric',length=length(timevar))
tunit=unique(timeunits)
for (i in 1:nooffiles)
{
if(i==1)
{
datenow=1
dateend=timevar[i]
file=files[i]
nc=open.ncdf(file)
start=c(xpos,ypos,1)
count=c(1,1,-1)
timseries[datenow:dateend]=get.var.ncdf(nc,varname,start=start,count=count)
tval[datenow:dateend]=get.var.ncdf(nc,"time")
close.ncdf(nc)
dateend=timevar[i]
datesx[i]=datenow
datesy[i]=dateend
}else
{
datenow=datenow+timevar[i-1]
dateend=dateend+timevar[i]
file=files[i]
nc=open.ncdf(file)
start=c(xpos,ypos,1)
count=c(1,1,-1)
timseries[datenow:dateend]=get.var.ncdf(nc,varname,start=start,count=count)
tval[datenow:dateend]=get.var.ncdf(nc,"time")
close.ncdf(nc)
datesx[i]=datenow
datesy[i]=dateend
}
}
return(list(timseries,tunit,tval))
}
# -----------------------------------------------------------------
# -----------------------------------------------------------------
# -----------------------------------------------------------------

precipitationBIASfactorfunction=function(demfile,stationfile,workdir)
{
latlonposition=latlonpos(demfile,stationfile)
library(lubridate)
readfolder=unlist(strsplit(stationfile, "\\\\"))
setwd(paste(file.path(as.vector(unlist(strsplit(stationfile, "\\\\"))
[1:length(unlist(strsplit(stationfile, "\\\\")))-1])),collapse="\\\\"))
stations=read.csv(stationfile,sep=",",colClasses = "character")
for (i in 1:dim(stations)[1])
{
if (!file.exists(paste("prec",stations[i,1],".csv",sep="")))
{
stop("File ", paste("prec",stations[i,1],".csv",sep="")," doesn't exist")
}
}
temparr=array(NA,c(12,dim(stations)[1]))
nctemparr=array(NA,c(12,dim(stations)[1]))
for (i in 1:dim(stations)[1])
{
data=read.csv(paste("Prec",stations[i,1],".csv",sep=""),stringsAsFactors =FALSE)
agg=aggregate(data$Prec,by=list(month(as.Date(data$Date,format="%m/%d/%Y"))),FUN=mean, na.rm=T) 
temparr[1:12 %in% agg$Group.1,i]=agg$x*0.0254
}

for(j in 1:length(latlonposition[[1]]))
{
varname="rain"
xpos=latlonposition[[1]][j]
ypos=latlonposition[[2]][j]
values=getval(workdir,xpos,ypos,varname)
splt <- as.POSIXlt(unlist(strsplit(unlist(values[[2]]), " "))[3])
tseries=ew.lt = splt + unlist(values[[3]])*3600
ncdata=data.frame(tseries,values[[1]])
colnames(ncdata)=c("Date","val")
aggncdata=aggregate(ncdata$val,by=list(month(as.Date(ncdata$Date))),FUN=mean, na.rm=T)
nctemparr[1:12 %in% aggncdata$Group.1,j]=aggncdata$x*24
}
#biasfactorallgauge=nctemparr/temparr
biasfactorallgauge=(temparr - nctemparr)/temparr
biasoverws=rowMeans(biasfactorallgauge)
return(list(biasoverws))
}

# -----------------------------------------------------------------
# -----------------------------------------------------------------
# -----------------------------------------------------------------
#install R packages if not installed
is.installed <- function(mypkg)
  {
    is.element(mypkg, installed.packages()[,1])
  }
packages=c("ncdf","raster","rgdal","gWidgetsRGtk2","date","fields","spam",
           "RAtmosphere","PBSmapping","stringr","R.utils","lubridate","png",
           "PBSmapping","RgoogleMaps")

url='http://cran.us.r-project.org'
# check if package "hydroGOF" is installed
for(i in 1:length(packages))
{
if (!is.installed(packages[i])){install.packages(packages[i],repos = url)}
}

# --------------------------------------------------
# --------------------------------------------------

# window dwsign
require(gWidgets)
options("guiToolkit"="RGtk2")
#---------------------------------------------------------------------------
#---------------------------------------------------------------------------
# funcrion to download the data
downloadMERRARFE=function(DestDir,latsouth,latnorth,lonwest,loneast,startdate,enddate,PrecipDataSource)
{
if(PrecipDataSource=="MERRA")
{
dir.create(file.path(DestDir,"MERRAPrecip"),showWarnings = FALSE,recursive = TRUE)
}
if(PrecipDataSource=="RFE")
{
dir.create(file.path(DestDir,"RFE"),showWarnings = FALSE,recursive = TRUE)
}
dir.create(file.path(DestDir,"MERRArad"),showWarnings = FALSE,recursive = TRUE)
dir.create(file.path(DestDir,"MERRAcli"),showWarnings = FALSE,recursive = TRUE)

conj='%2C'
tol=0
boundbox=paste(latsouth-tol,conj,lonwest-tol,conj,latnorth+tol,conj,loneast+tol,sep='')

Date=seq(startdate,enddate,by='days')

for (i in 1:length(Date))
{
linkyearform=substr(as.character(Date[i]), 1, 4)
linkmonthform=substr(as.character(Date[i]), 6, 7)
linkdayform=substr(as.character(Date[i]), 9, 10)
datepart=paste(linkyearform,linkmonthform,linkdayform,sep='')

# RFE precipitation data download
#---------------------------------------
if(PrecipDataSource=="RFE")
{
mpart='ftp://ftp.cpc.ncep.noaa.gov/fews/S.Asia/data/cpc_rfe_v2.0_sa_dly.bin.'
filepart='cpc_rfe_v2.0_sa_dly.bin.'
lpart='.gz'
#form your links
filename=paste(filepart,datepart,sep='')
links=paste(mpart,datepart,lpart,sep='')
rfeLocLink=file.path(DestDir,"RFE",paste(filepart,datepart,lpart,sep=''))
#---------------------------------------
#download your data and upzip it
errorcatch=tryCatch(download.file(links,rfeLocLink),error=function(e) e)
errorcatch=as.character(errorcatch)
if(errorcatch != as.character(0))
{
file.copy(prevfile,rfeLocLink,overwrite = FALSE)
}
prevfile=rfeLocLink
#---------------------------------------
}

# MERRA precipitation data download
#---------------------------------------
if(PrecipDataSource=="MERRA")
{
#Parts of download link
# forming the link to download data
p1="http://goldsmr2.sci.gsfc.nasa.gov/daac-bin/OTF/HTTP_services.cgi?FILENAME=%2Fdata%2Fs4pa%2FMERRA%2FMAT1NXLND.5.2.0%2F"
p2="%2F"
p3=   "%2FMERRA300.prod.assim.tavg1_2d_lnd_Nx."
p3_ed="%2FMERRA301.prod.assim.tavg1_2d_lnd_Nx."
p4=".hdf&FORMAT=TmV0Q0RGLw&BBOX="
p5=   "&LABEL=MERRA300.prod.assim.tavg1_2d_lnd_Nx."
p5_ed="&LABEL=MERRA300.prod.assim.tavg1_2d_lnd_Nx."
p6=".SUB.nc&FLAGS=&SHORTNAME=MAT1NXLND&SERVICE=SUBSET_LATS4D&LAYERS=&VERSION=1.02&VARIABLES=prectot"
#---------------------------------------
datep=paste(linkyearform,linkmonthform,linkdayform,sep='')
link=paste(p1,linkyearform,p2,linkmonthform,p3,datep,p4,boundbox,p5,datep,p6,sep='')
outputfile=paste('merra.prod.precip.',datep,'_dwd.nc',sep='')
MerraCliLocLink=file.path(DestDir,"MERRAPrecip",outputfile)
#download your data
errorcatch=tryCatch(download.file(link,MerraCliLocLink,mode='wb'),error=function(e) e)
errorcatch=as.character(errorcatch)
if(errorcatch != as.character(0))
{
link=paste(p1,linkyearform,p2,linkmonthform,p3_ed,datep,p4,boundbox,p5_ed,datep,p6,sep='')
errorcatch1=tryCatch(download.file(link,MerraCliLocLink,mode='wb'),error=function(e) e)
}
errorcatch1=as.character(errorcatch)
if(errorcatch1 != as.character(0))
{
print(paste(outputfile," could not be downloaded"))
}
#---------------------------------------
}

# MERRA Climate data (temperature, pressure, wind and dew point) download
#---------------------------------------
#Parts of download link
# forming the link to download data
p1="http://goldsmr2.sci.gsfc.nasa.gov/daac-bin/OTF/HTTP_services.cgi?FILENAME=%2Fdata%2Fs4pa%2FMERRA%2FMAT1NXSLV.5.2.0%2F"
p2="%2F"
p3=   "%2FMERRA300.prod.assim.tavg1_2d_slv_Nx."
p3_ed="%2FMERRA301.prod.assim.tavg1_2d_slv_Nx."
p4=".hdf&FORMAT=TmV0Q0RGLw&BBOX="
p5=   "&LABEL=MERRA300.prod.assim.tavg1_2d_slv_Nx."
p5_ed="&LABEL=MERRA301.prod.assim.tavg1_2d_slv_Nx."
p6=".SUB.nc&FLAGS=&SHORTNAME=MAT1NXSLV&SERVICE=SUBSET_LATS4D&LAYERS=&VERSION=1.02&VARIABLES=ps%2Cu2m%2Cv2m%2Ct2m%2Cqv2m"
#---------------------------------------
datep=paste(linkyearform,linkmonthform,linkdayform,sep='')
link=paste(p1,linkyearform,p2,linkmonthform,p3,datep,p4,boundbox,p5,datep,p6,sep='')
outputfile=paste('merra.prod.assim.',datep,'_dwd.nc',sep='')
MerraCliLocLink=file.path(DestDir,"MERRAcli",outputfile)
#download your data
errorcatch=tryCatch(download.file(link,MerraCliLocLink,mode='wb'),error=function(e) e)
errorcatch=as.character(errorcatch)
if(errorcatch != as.character(0))
{
link=paste(p1,linkyearform,p2,linkmonthform,p3_ed,datep,p4,boundbox,p5_ed,datep,p6,sep='')
errorcatch1=tryCatch(download.file(link,MerraCliLocLink,mode='wb'),error=function(e) e)
}
errorcatch1=as.character(errorcatch)
if(errorcatch1 != as.character(0))
{
print(paste(outputfile," could not be downloaded"))
}
#---------------------------------------


# MERRA incoming radiation (shortwave and longwave) download
#---------------------------------------
p1='http://goldsmr1.sci.gsfc.nasa.gov/daac-bin/OTF/HTTP_services.cgi?FILENAME=%2Fdata%2Fs4pa%2FMERRA%2FMAT3FXCHM'
p2='.5.2.0%2F'
p3='%2F'
p4='%2FMERRA300.prod.assim.tavg3_2d_chm_Fx.'
p4_ed='%2FMERRA301.prod.assim.tavg3_2d_chm_Fx.'
p5='.hdf&FORMAT=TmV0Q0RGLw&BBOX='
p6='&LABEL=MERRA300.prod.assim.tavg3_2d_chm_Fx.'
p6_ed='&LABEL=MERRA301.prod.assim.tavg3_2d_chm_Fx.'
p7='.SUB.nc&FLAGS=&SHORTNAME=MAT3FXCHM&SERVICE=SUBSET_LATS4D&LAYERS=&VERSION=1.02&VARIABLES=lwgdwn%2Cswgdwn'
#---------------------------------------
link=paste(p1,p2,linkyearform,p3,linkmonthform,p4,datep,p5,boundbox,p6,datep,p7,sep='')
outputfile=paste('MERRA300.prod.assim.tavg3_2d_chm_Fx.',datep,'.SUB_long.nc',sep='')
MerraRadLocLink=file.path(DestDir,"MERRArad",outputfile)
#download your data
errorcatch=tryCatch(download.file(link,MerraRadLocLink,mode='wb'),error=function(e) e)
errorcatch=as.character(errorcatch)
if(errorcatch != as.character(0))
{
link=paste(p1,p2,linkyearform,p3,linkmonthform,p4_ed,datep,p5,boundbox,p6_ed,datep,p7,sep='')
errorcatch1=tryCatch(download.file(link,MerraRadLocLink,mode='wb'),error=function(e) e)
}
errorcatch1=as.character(errorcatch)
if(errorcatch1 != as.character(0))
{
print(paste(outputfile," could not be downloaded",sep=' '))
}
}
}

#---------------------------------------------------------------------------
#---------------------------------------------------------------------------

vegfiles=function(landcoverfile,demfile,landcoverclassfile,outfolder)
{
landcoverfile=landcoverfile
demfile=demfile
landcoverclassfile=landcoverclassfile
outfolder=outfolder
curdir=getwd()
dir.create(outfolder,showWarnings = FALSE,recursive = TRUE)
setwd(outfolder)

# --------------------------------------------------
# --------------------------------------------------
# names of the netcdf files for cc, hcan, LAI, ycage
ccnet='ccgrid.nc'
hcannet='hcangrid.nc'
LAInet='LAIgrid.nc'
ycagenet='ycagegrid.nc'

varnames=c('cc','hcan','lai','ycage')
varunits=c('','m','','')
varlongnames=c('canopy coverage','canopy height','leaf area index','Forest canopy structure flag')
filenames=c('ccgrid.nc','hcangrid.nc','LAIgrid.nc','ycagegrid.nc')
# --------------------------------------------------
# --------------------------------------------------

library(rgdal)
library(raster)
library(ncdf)

# Extracting and summarizing the grid of DEM for the output data
landfile=raster(readGDAL(landcoverfile))
#landproj=projection(landfile)
dem=raster(readGDAL(demfile))
demproj=projection(dem)
landfileprojected <- projectRaster(landfile,dem, method="ngb")
#plot(landfileprojected)

#  Set defaults
proj='longlat'
xcoordinate='longitude'
ycoordinate='latitude'
xcoordinateunit='degrees_east'
ycoordinateunit='degrees_north'

#  Define a trim function
trim <- function (x) gsub("^//s+|//s+$", "", x)

# Now look at info in dem projection detail and adjust units if needed
demprojparts=unlist(strsplit(demproj, "[+]"))
for(str in demprojparts)
{
sub12=unlist(strsplit(str, "[=]"))
if(length(sub12)>1)if(sub12[1]=="proj")proj=trim(sub12[2])
if(length(sub12)>1)if(sub12[1]=="units")
{ xcoordinateunit=trim(sub12[2])
  ycoordinateunit=xcoordinateunit
}
}
if(proj!='longlat')  # if projection not long lat write with x and y as dimension names
{
 xcoordinate='x'
 ycoordinate='y'
}


# obtain DEM extent
demres=res(landfileprojected)
x1min=xmin(landfileprojected)
x1max=xmax(landfileprojected)
y1min=ymin(landfileprojected)
y1max=ymax(landfileprojected)
#  Extract cell midpoint values for ncdf file
x_vals = seq(x1min+demres[1]/2,x1max-demres[1]/2,demres[1])
y_vals = seq(y1min+demres[2]/2,y1max-demres[2]/2,demres[2])

# MODIS land cover classes
# https://lpdaac.usgs.gov/products/modis_products_table/mcd12q1
#0  Water
#1	Evergreen Needleleaf forest
#2	Evergreen Broadleaf forest
#3	Deciduous Needleleaf forest
#4	Deciduous Broadleaf forest
#5	Mixed forest
#6	Closed shrublands
#7	Open shrublands
#8	Woody savannas
#9	Savannas
#10	Grasslands
#11	Permanent wetlands
#12	Croplands
#13	Urban and built-up
#14	Cropland/Natural Land Cover mosaic
#15	Snow and ice
#16	Barren or sparsely vegetated

vegtype=c("Water","Evergreen Needleleaf forest","Evergreen Broadleaf forest",
          "Deciduous Needleleaf forest","Deciduous Broadleaf forest",
          "Mixed forest","Closed shrublands","Open shrublands","Woody savannas","Savannas",
          "Grasslands","Permanent wetlands","Croplands","Urban and built-up",
          "Cropland/Natural Land Cover mosaic","Snow and ice",
          "Barren or sparsely vegetated")

if(landcoverclassfile=='')
{
  landcoverclass=0:16
  cc=c(0,0.7,0.7,0.7,0.7,0.8,0,0,0,0,0,0,0,0,0,0,0)
  hcan=c(0,15,15,15,15,10,0,0,0,0,0,0,0,0,0,0,0)
  LAI=c(0,4.5,4.5,4.5,4.5,4,0,0,0,0,0,0,0,0,0,0,0)
  ycage=c(2,3,2,3,2,2,2,2,2,2,2,2,2,2,2,2,2)
  landcovertable<-data.frame(vegtype,landcoverclass,cc,hcan,LAI,ycage)
  #write.csv(landcovertable, file = landcoverclassfile,row.names = FALSE)
}else
{
  landcovertable=read.csv(landcoverclassfile,header=T)
  vegtype=landcovertable$vegtype
  landcoverclass=landcovertable$landcoverclass
  cc=landcovertable$cc
  hcan=landcovertable$hcan
  LAI=landcovertable$LAI
  ycage=landcovertable$ycage
}
#--------------------------------------------------

# work on your cc grid
ccraster=0*landfileprojected
for(i in 1:length(landcoverclass))
{
  ccraster[landfileprojected == landcoverclass[i]]=cc[i]
}
#plot(ccraster)
# --------------------------------------------------

# work on your hcan grid
hcanraster=0*landfileprojected
for(i in 1:length(landcoverclass))
{
  hcanraster[landfileprojected == landcoverclass[i]]=hcan[i]
}
#plot(hcanraster)
# --------------------------------------------------

# work on your LAI grid
LAIraster=0*landfileprojected
for(i in 1:length(landcoverclass))
{
  LAIraster[landfileprojected == landcoverclass[i]]=LAI[i]
}
#plot(LAIraster)
# --------------------------------------------------
# work on your ycage grid
ycageraster=0*landfileprojected
for(i in 1:length(landcoverclass))
{
  ycageraster[landfileprojected == landcoverclass[i]]=ycage[i]
}
#plot(ycageraster)
# --------------------------------------------------
# create list of all rasters
rasterlist=list(ccraster,hcanraster,LAIraster,ycageraster)

missval=-9999 # missing value
# define the netcdf coordinate variables -- note these have values!
dim1 = dim.def.ncdf(xcoordinate,xcoordinateunit, as.double(x_vals))
dim2 = dim.def.ncdf(ycoordinate,ycoordinateunit, as.double(y_vals))
for (i in 1:length(varnames))
{
  outvector=as.vector(t(as.matrix(rasterlist[[i]])[length(y_vals):1,]))  
  # first it flip it 
  # Then transpose to account for by column internal to R whearas ncdf takes by row
  var = var.def.ncdf(varnames[i],varunits[i],list(dim1,dim2),missval,varlongnames[i])
  nc.ex = create.ncdf(filenames[i],var)
  put.var.ncdf(nc.ex,var,outvector)
  close.ncdf(nc.ex)
}
setwd(curdir)
}
#--------------------------------------------------
#--------------------------------------------------
#-------------------------------------------------------------------
#-------------------------------------------------------------------
#-------------------------------------------------------------------
# All Downscaling related functiond are here
#-------------------------------------------------------------------
#-------------------------------------------------------------------
#-------------------------------------------------------------------
#-------------------------------------------------------------------
#-------------------------------------------------------------------
#-------------------------------------------------------------------
#-------------------------------------------------------------------
#-------------------------------------------------------------------

julianfunc = function(yy,mm,dd)
{
mmstrt=c(0,31,59,90,120,151,181,212,243,273,304,334)
jday = mmstrt[mm] + dd
ileap = yy - as.integer(yy/4)*4 
if(ileap==0 & mm > 3) 
{
 jday = jday + 1
}
julian = jday
return (julian)
}


HRIfunc = function(year,month,hour,day,dt,lat,slope=0,azi=0)
{
library(date)
date=mdy.date(month,day,year)
dateform=as.Date(paste(year,'-',month,'-',day,sep=''))
crad=pi/180.0
T=(hour-12.0)*pi/12.0
delt1=dt*pi/12.0

slope1=slope*crad
azi1=azi*crad

lat1=lat*crad
fjulian=julianfunc(year,month,day)
d=crad*23.5*sin((fjulian-82.0)*(2*pi/365))
# D is solar declination
lp=asin(sin(slope1)*cos(azi1)*cos(lat1)+cos(slope1)*sin(lat1))
tanprod=tan(lat1)*tan(d)
if(tanprod > 1.0)
  {
   td=pi              # This is the condition for perpetual light
  }else if(tanprod < -1.0){
   td=0               # The condition for perpetual night
  }else{
   td=acos(-tanprod)  # The condition where there is a sunrise and set
  }
ddt=atan2(sin(azi1)*sin(slope1),(cos(slope1)*cos(lat1)-cos(azi1)*sin(slope1)*sin(lat1)))  
tpeqp=tan(lp)*tan(d)
if(tpeqp > 1.0)
{
   tpbeg=-pi   # perpetual light
   tpend=pi
}else if(tpeqp < -1.0){
   tpbeg=0.0   # perpetual dark
   tpend=0.0 
}else{
   tpbeg = -acos(-tpeqp)-ddt
   tpend = acos(-tpeqp)-ddt
}

T1=max(T,tpbeg,-td)
T2=min(T+delt1,td,tpend)
if( T2 < T1)
{
 HRI=0.0
}else{
 HRI=(sin(d)*sin(lp)*(T2-T1)+cos(d)*cos(lp)*(sin(T2+ddt)-sin(T1+ddt)))/(cos(slope1)*delt1)
}

if(tpbeg < -pi)
{
 T1=max(T,-tpbeg+2*pi,-td)
 T2=min(T+DELT1,TD)
 if(t2 > t1)
 {
  HRI=HRI+(sin(D)*sin(LP)*(T2-T1)+cos(D)*cos(LP)*(sin(T2+ddt)-SIN(T1+ddt)))/(cos(slope1)*delt1)
 }
}

zenith=acos(HRI)*180/pi
return(list(HRI,zenith))
}

# function to retrieve a subset of a grid from a netcdf file
ncarray=function(ncobj,dem,varname,timeindex=1,xdimname="longitude",ydimname="latitude",tdimname="time")
{
# Get variable names
varnames=rep("",ncobj$nvars)
for(i in 1:ncobj$nvars){varnames[i]=ncobj$var[[i]]$name}
#varnames
varid=which(varnames==varname)
# get dimension associate with the variable
ndimv=ncobj$var[[varid]]$ndims
dimnames=rep("",ndimv)
for(i in 1:ndimv){
dimnames[i]=ncobj$dim[[ncobj$var[[varid]]$dimids[i]]]$name}
# Find the indices of our dimensions for our variable
dimidx=which(dimnames==xdimname)
dimidy=which(dimnames==ydimname)
dimidt=which(dimnames==tdimname)
xvals=get.var.ncdf(ncobj,varid=xdimname)
yvals=get.var.ncdf(ncobj,varid=ydimname)
nx=length(xvals)
ny=length(yvals)
#
count=rep(0,2+length(dimidt))  #  Use length(dimidt) to identify whether this is a 2D or 3D file
start=rep(1,2+length(dimidt))
count[dimidx]=nx
count[dimidy]=ny
count[dimidt]=1
start[dimidt]=timeindex
#  Get the variable as an array in to R to flip and subset as this is hard to do generally
#  for read from ncdf when dimensions may be in different orders and decreasing or increasing
vals=get.var.ncdf(ncobj,varid=varname,count=count,start=start)
# Because we always get one time step, the resulting array is 2D.  
# Need to flip or transpose depending on dimension values to put in image format (indexed from bottom left with x first dimension)
# If dimension of the first subscript variable in ncfile is latitude - then transpose
if(dimidx>dimidy){
vals=t(vals)
}

dy=(yvals[ny]-yvals[1])/(ny-1)
dx=(xvals[nx]-xvals[1])/(nx-1)

#  From here first subscript is x and second subscript y
if(dx < 0){  # x values are decreasing so flip E-W
vals=vals[nx:1,]
dx=-dx
}
if(dy < 0){
vals=vals[,ny:1]
dy=-dy
}
# now sort the dimension values 
xvals=sort(xvals)
yvals=sort(yvals)  
list(x=xvals,y=yvals,z=vals)  # This is the standard R way of arranging 2D arrays for plotting
}

rasterFromArray=function(x=x,y=y,z=z,crs=crs)
{
# function to output a raster from data in standard R image plotting format
ny=length(y)
nx=length(x)
dy=(y[ny]-y[1])/(ny-1)
dx=(x[nx]-x[1])/(nx-1)
minx=min(x)-dx/2
maxx=max(x)+dx/2
miny=min(y)-dy/2
maxy=max(y)+dy/2
#  Transpose and flip array to put in raster
res = raster(t(z[,ny:1]),xmn=minx,xmx=maxx,ymn=miny,ymx=maxy,crs=crs)
}


DownscaleMERRARFE <- function(VarSourceDir,GeoSourceDir,DestDir,DEMfile,startdate,enddate,lapse,adjfactor,latlonfile,PrecipDataSource) 
{
nstep=8
dir.create(file.path(DestDir),showWarnings = FALSE,recursive = TRUE)
dir.create(file.path(DestDir,"temp"),showWarnings = FALSE,recursive = TRUE)
CodeStartTime=as.numeric(Sys.time())
#source(file.path(mydir, "Rcode", "RegridInput.R"))
#--------------------------------------------------------------
#--------------------------------------------------------------
#--------------------------------------------------------------
#--------------------------------------------------------------
#--------------------------------------------------------------


# Load needed R libraries
library(stringr)
library(R.utils)
library(lubridate)
#---------------------------------------
startdate=as.Date(startdate)
enddate=as.Date(enddate)
Date=seq(startdate,enddate,by='days')

curdir=getwd()
#  Do all work in temporary directory then revert back at end
setwd(file.path(DestDir,"temp"))
dir.create(file.path(DestDir,"regrid"),showWarnings = FALSE,recursive = TRUE)

#---------------------------------------
# LOOP OVER DATES STARTS HERE
#---------------------------------------
for (idate in 1:length(Date))
{
  print(paste("Regridding:",Date[idate]))
  linkyearform=substr(as.character(Date[idate]), 1, 4)
  linkmonthform=substr(as.character(Date[idate]), 6, 7)
  linkdayform=substr(as.character(Date[idate]), 9, 10)
  datepart=paste(linkyearform,linkmonthform,linkdayform,sep='')
  dt=24/nstep
  #---------------------------------------
  # RFE EXTRACTION AND REGRIDDING
  #---------------------------------------
  # Construct file name based on RFE naming system
  if(PrecipDataSource=="RFE")
  {
  filepart='cpc_rfe_v2.0_sa_dly.bin.'
  lpart='.gz'
  filename=paste(filepart,datepart,sep='')
  mylink=file.path(VarSourceDir,"RFE",paste(filepart,datepart,lpart,sep=''))
  outfile=file.path(DestDir,"regrid",paste('rfe2.sa.',datepart,'.nc',sep=''))  
  tempfile=filename
  ctrlfile="prectot.ctl" 
  #---------------------------------------
  #change the control file "..\\controls\\prectot.ctl"
  #list <- readLines(file.path(mydir,"controls","prectot.ctl"),warn = F)
  newline1="TDEF      1 LINEAR   "
  mon=c('jan','feb','mar','apr','may','jun','jul','aug',
        'sep','oct','nov','dec')
  dayin=day(Date[idate])
  monthin=mon[month(Date[idate])]
  yearin=year(Date[idate])
  newline2=paste(dayin,monthin,yearin,sep='')
  newline3=' 1dy'
  newline=paste(newline1,newline2,newline3,sep='')

  #  Control file content
  cfilecontent = c(
  "DSET  cpc_rfe_v2.0_sa_dly.bin.%y4%m2%d2",                   
  "options big_endian template",                              
  "TITLE  RFE 2.0 precipitation forecast lead time: 1 day(s) ",
  "undef -999.0 ",                                             
  "XDEF    401 LINEAR   70 0.1 ",                              
  "YDEF    301 LINEAR    5 0.1",                               
  "ZDEF      1 LINEAR    1   1",                              
  newline,                         
  "VARS      1",                                              
  "rain   0  0  daily precip (mm/d)",                         
  "ENDVARS" )  
  writeLines(cfilecontent,ctrlfile)
  #---------------------------------------
    #download your data and upzip it
    gunzip(mylink,tempfile,remove=F)
    #---------------------------------------
    #run CDO commands to convert binary to netCDF
    #  -s is silent flag to suppress output
    system("cdo -s -f nc import_binary prectot.ctl Output1stStepfile.nc")  
    #---------------------------------------
    #run CDO commands to subset lat-lon
    #system("cdo sellonlatbox,84.0,87.0,27,30 Output1stStep.nc Output2ndStep.nc")
    #---------------------------------------
    #print("CDO change time for each file")
    outlist=""
    for(i in 1:nstep)
    {
    ttt=(i-1)*dt
    hr=trunc(ttt)
    min=(ttt-hr)*60
    tst=paste(hr,":",formatC(min,width=2,flag="0"),sep="")
    ofile=paste("ofile",i-1,".nc",sep="")
    cmd1=paste("cdo -s settime,",tst," Output1stStepfile.nc ","ofile.nc",sep="")
    system(cmd1)
    system(paste("cdo -s divc,24000 ofile.nc",ofile))   # convert units 
    # ** dividng by 1000 for mm to m ** dividng by 24 for day to hr
    outlist=paste(outlist,ofile)
    }
    #---------------------------------------
    #print(paste("Write to", outfile))
    cmd=paste("cdo -s -f nc copy",outlist,outfile, sep=" ")
    system(cmd)
    #  Delete temporary files
    unlink(list.files())
    }
    #---------------------------------------
    #---------------------------------------


    #---------------------------------------
    # MERRA PRECEPITATION REGRIDDING
    #---------------------------------------
    # Construct file name based on RFE naming system
    if(PrecipDataSource=="MERRA")
    {
    datep=paste(linkyearform,linkmonthform,linkdayform,sep='') 
    inputfile=file.path(VarSourceDir,"MERRAPrecip",paste('merra.prod.precip.',datep,'_dwd.nc',sep=''))
    outfile=paste('merra.prod.precip.',datep,'_0z-21z.nc',sep='')
    finaloutfile=paste('..\\regrid\\merra.prod.precip.',datep,'.nc',sep='')
    #---------------------------------------
    print("Averaging precipitation data over time step.")
    outlist=""
    for(i in 1:nstep)
    {
     tbeginning=(i-1)*dt
     tend=tbeginning+dt-1
     hr=trunc(tbeginning)
     min=(tbeginning-hr)*60
     tst=paste(hr,":",formatC(min,width=2,flag="0"),sep="")
     ofile=paste("Preciptempfile",i-1,".nc",sep="")
     cmd1=paste('ncra -d time',',',paste(tbeginning,".",sep=''),',',paste(tend,".",sep=''),' ',inputfile," ",ofile,sep='')  
     system(cmd1)
     ofilep=paste("Preciptempfile",i-1,"t",".nc",sep="")
     cmd2=paste("cdo -s settime,",tst,' ',ofile," ",ofilep,sep="")
     system(cmd2)
     outlist=paste(outlist,ofile)
    }
    #---------------------------------------
    cmd=paste("cdo -s -f nc copy",outlist,outfile, sep=" ")
    system(cmd)
    #---------------------------------------
    Tcmd=paste("cdo -s expr,rain=prectot*3.6",outfile,"tempP.nc")
    system(Tcmd)
    finalP=paste('ncks -h -A tempP.nc -o ',finaloutfile,sep="")
    system(finalP)
    prectotunit=paste('ncatted -a units,rain,o,c,"m hr-1"',finaloutfile)
    system(prectotunit)
    }
    #---------------------------------------
    #---------------------------------------

    #---------------------------------------
    # MERRA CLIMATE VARIABLE EXTRACTION AND REGRIDDING
    #---------------------------------------
    datep=paste(linkyearform,linkmonthform,linkdayform,sep='') 
    inputfile=file.path(VarSourceDir,"MERRAcli",paste('merra.prod.assim.',datep,'_dwd.nc',sep=''))
    # inputfile=paste('merra.prod.assim.',datep,'_dwd.nc',sep='')

    outfile=paste('merra.prod.assim.',datep,'_0z-21z.nc',sep='')

    #---------------------------------------
    print("Averaging climate data over time step.")
    outlist=""
    for(i in 1:nstep)
    {
     tbeginning=(i-1)*dt
     tend=tbeginning+dt-1
     hr=trunc(tbeginning)
     min=(tbeginning-hr)*60
     tst=paste(hr,":",formatC(min,width=2,flag="0"),sep="")
     ofile=paste("tempfile",i-1,".nc",sep="")
     cmd1=paste('ncra -d time',',',paste(tbeginning,".",sep=''),',',paste(tend,".",sep=''),' ',inputfile," ",ofile,sep='')  
     system(cmd1)
     ofilet=paste("tempfile",i-1,"t",".nc",sep="")
     cmd2=paste("cdo -s settime,",tst,' ',ofile," ",ofilet,sep="")
     system(cmd2)
     outlist=paste(outlist,ofile)
    }
    #---------------------------------------
    #print(paste("Writing to", outfile))
    cmd=paste("cdo -s -f nc copy",outlist,outfile, sep=" ")
    system(cmd)
    #---------------------------------------
    #print(paste("Put sp. humidity and pressure without any modification into", finaloutfile))
    outfile=file.path(getwd(),outfile)
    finaloutfile=paste('merra.prod.assim.',datep,'.nc',sep="")
    setwd(file.path(DestDir,"regrid"))
    PRhcmd=paste("ncks -h -a -O -v ps,qv2m ",outfile," ",finaloutfile,sep="")
    system(PRhcmd)
    setwd(file.path(DestDir,"temp"))
    #print("Calculate wind speed")
    windcmd=paste("cdo -s expr,wind=sqrt(u2m*u2m+v2m*v2m)",outfile,"tempwind.nc")
    windcmd2=paste("cdo -s expr,windu=u2m",outfile,"tempwindu.nc")
    windcmd3=paste("cdo -s expr,windv=v2m",outfile,"tempwindv.nc")
    #windcmd=paste("ncap -v -O -h -s wind[time,latitude,longitude]=sqrt(u2m*u2m+v2m*v2m) " ,outfile," tempwind.nc",sep="")
    system(windcmd)
    system(windcmd2)
    system(windcmd3)
    #---------------------------------------
    # Calculate temperature (unit conversion: Kelvin to Celsius)
    #Tcmd=paste("ncap -v -O -h -s temp[time,latitude,longitude]=t2m-273 ",outfile," tempT.nc",sep="")
    Tcmd=paste("cdo -s expr,temp=t2m-273.15",outfile,"tempT.nc")
    system(Tcmd)
    #---------------------------------------
    # put temperature and wind in output netCDF
    TWoutfile=file.path(getwd(),"tempwind.nc")
    TWoutfileu=file.path(getwd(),"tempwindu.nc")
    TWoutfilev=file.path(getwd(),"tempwindv.nc")
    Toutfile=file.path(getwd(),"tempT.nc")
    setwd(file.path(DestDir,"regrid"))
    finalwind=paste('ncks -h -A ', TWoutfile, ' -o ',finaloutfile,sep="")
    finalwindu=paste('ncks -h -A ', TWoutfileu, ' -o ',finaloutfile,sep="")
    finalwindv=paste('ncks -h -A ', TWoutfilev, ' -o ',finaloutfile,sep="")
    finalT=paste('ncks -h -A ', Toutfile, ' -o ',finaloutfile,sep="")
    system(finalwind)
    system(finalwindu)
    system(finalwindv)
    system(finalT)
    setwd(file.path(DestDir,"temp"))
    #---------------------------------------
    #Change units in netCDF file
    setwd(file.path(DestDir,"regrid"))
    windunit=paste('ncatted -a units,wind,o,c,"m s-1"',finaloutfile)
    punit=paste('ncatted -a units,ps,o,c,"Pa"',finaloutfile)
    sphunit=paste('ncatted -a units,qv2m,o,c,"kg kg-1"',finaloutfile)
    system(windunit)
    system(punit)
    system(sphunit)
    #---------------------------------------

    #---------------------------------------
    # MERRA RADIATION WORK STARTS HERE
    #---------------------------------------
   
    inputfile=file.path(VarSourceDir,"MERRArad",paste('MERRA300.prod.assim.tavg3_2d_chm_Fx.',datep,'.SUB_long.nc',sep='')) 
    longout=paste(file.path(DestDir,"temp"),'/merra.lwrad.',datep,'.nc',sep='')
    shortout=paste(file.path(DestDir,"temp"),'/merra.swrad.',datep,'.nc',sep='')
    finaloutfile=paste('merra.prod.rad.',datep,'.nc',sep='')
    #---------------------------------------
    # Radiation unit conversion (from Wm-2 to kJ h-1 m-2)
    longcmd=paste("cdo expr,lwrad=lwgdwn*3.6",inputfile,longout)
    system(longcmd)
    shortcmd=paste("cdo -s expr,swrad=swgdwn*3.6",inputfile,shortout)
    system(shortcmd)
    #---------------------------------------
    setwd(file.path(DestDir,"regrid"))
    # put incoming shortwave and longwave in output netCDF
    storelongcmd=paste('ncks -h -A',longout,'-o',finaloutfile)
    system(storelongcmd)
    storeshortcmd=paste('ncks -h -A',shortout,'-o',finaloutfile)
    system(storeshortcmd)
    #---------------------------------------
    #Change units in netCDF file
    unitlongcmd=paste('ncatted -a units,lwrad,o,c,"kJ h-1 m-2"',finaloutfile)
    system(unitlongcmd)
    unitshortcmd=paste('ncatted -a units,swrad,o,c,"kJ h-1 m-2"',finaloutfile)
    system(unitshortcmd)
    #---------------------------------------
    #  Delete temporary files
    setwd(file.path(DestDir,"temp"))
    unlink(list.files())
    #---------------------------------------
}
setwd(curdir)


#--------------------------------------------------------------
#--------------------------------------------------------------
#--------------------------------------------------------------
#--------------------------------------------------------------
#--------------------------------------------------------------

CodeEndTime=as.numeric(Sys.time())
RegridTime=CodeEndTime-CodeStartTime
CodeStartTime=CodeEndTime

#source(file.path(mydir, "Rcode", "DownScalePhysics.R"))
#--------------------------------------------------------------
#--------------------------------------------------------------
#--------------------------------------------------------------
#--------------------------------------------------------------
#--------------------------------------------------------------

#---------------------------------------
# Attach needed libraries

library(raster)
library(rgdal)
library(date)
library(fields)
library(ncdf)
library(RAtmosphere)


monthDays <- function(d) 
{
   last_days <- 28:31
   rev(last_days[which(!is.na(as.Date(paste(substr(d,1,8),last_days,sep = ''),'%Y-%m-%d')))])[1]
}

#Set up loop parameters 
dateList = seq(as.Date(startdate),as.Date(enddate),by="1 day")
ndays = c(length(dateList))
dt=24/nstep
referencedate=as.Date('2001-05-01')
firstValue=as.numeric(startdate-referencedate)*24

# set input and output directories
diri = file.path(DestDir,"regrid")
if(!file.exists(diri)){print(paste("Regridded input directory missing:",diri))}

dirir=diri

diro = file.path(DestDir,"tfirst")
dirolast= file.path(DestDir,"tfirstBiasCorrected")
dir.create(diro,showWarnings = FALSE,recursive = TRUE)

#Parameters
grav=9.81
rdry=287.
Sstar=1367.0
sigma=5.67E-8*3.6
H=8400.
lapse=lapse
#lapse<-c(0.0044,0.0059,0.0071,0.0078,0.0081,0.0082,0.0081,0.0081,0.0077,0.0068,0.0055,0.0047)
tdlapse<-c(0.00041,0.00042,0.00040,0.00039,0.00038,0.00036,0.00033,0.00033,0.00036,0.00037,0.00040,0.00040)
a=611.21
b=22.452
c=240.97
X = c(0.35,0.51) # Emissivity params
Y = c(0.1,0.13)  # Emissivity params
Z = c(0.224,1.1) # Emissivity params

crs = "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0"

# Extracting and summarizing the grid of DEM for the output data
if(!file.exists(DEMfile)){print(paste("Missing DEM file:",DEMfile))}
dem<-raster(readGDAL(DEMfile))
demproj=projection(dem)
slopewind=terrain(dem,out=c('slope'), unit='radians',neighbors=8) 
aspectwind=terrain(dem,out=c('aspect'), unit='radians',neighbors=8) 
reswind <- res(dem)
dxwind <- reswind[1]
dywind <- reswind[2]
diawind=sqrt(dxwind^2+dywind^2)
mat1wind=matrix(c(0,0,0,1,0,1,0,0,0)/(16*dywind),3,3)
mat2wind=matrix(c(0,1,0,0,0,0,0,1,0)/(16*dxwind),3,3)
mat3wind=matrix(c(1,0,1,0,0,0,1,0,1)/(2*sqrt(2)*mean(dywind,dxwind)),3,3)
mat4wind=matrix(c(0,0,0,0,(sqrt(2)+1)/(8*sqrt(2)*mean(dywind,dxwind)),0,0,0,0),3,3)
agg1wind <- focal(dem,mat1wind,sum,pad = T,padValue = 0)
agg2wind <- focal(dem,mat2wind,sum,pad = T,padValue = 0)
agg3wind <- focal(dem,mat3wind,sum,pad = T,padValue = 0)
agg4wind <- focal(dem,mat4wind,sum,pad = T,padValue = 0)
agg5wind=agg4wind-agg1wind+agg2wind+agg3wind
curvaturepar = (agg5wind-(cellStats(agg5wind,'max')+cellStats(agg5wind,'min'))/2)/
             (cellStats(agg5wind,'max')-cellStats(agg5wind,'min'))


#  Get projection information from DEM to decide on dimension names for netcdf
#  Set defaults
 proj='longlat'
 xcoordinate='longitude'
 ycoordinate='latitude'
 xcoordinateunit='degrees_east'
 ycoordinateunit='degrees_north'

#  Define a trim function
trim <- function (x) gsub("^//s+|//s+$", "", x)

# Now look at info in dem projection detail and adjust units if needed
demprojparts=unlist(strsplit(demproj, "[+]"))
for(str in demprojparts)
{
sub12=unlist(strsplit(str, "[=]"))
if(length(sub12)>1)if(sub12[1]=="proj")proj=trim(sub12[2])
if(length(sub12)>1)if(sub12[1]=="units")
{ xcoordinateunit=trim(sub12[2])
  ycoordinateunit=xcoordinateunit
}
}
if(proj!='longlat')  # if projection not long lat write with x and y as dimension names
{
 xcoordinate='x'
 ycoordinate='y'
}


# Before this function was called nudgeExtent
demres<-res(dem)
dimdem<-dim(dem)
ymax<-ymax(dem)
ymin<-ymin(dem)
xmin<-xmin(dem)
xmax<-xmax(dem)

# DGT lets keep these as separate variables to use xmin,xmax,ymin,ymax to determine bounds
y1max = ymax - demres[2] / 2
y1min = ymin + demres[2] / 2
x1max = xmax - demres[1] / 2
x1min = xmin + demres[1] / 2


lowdem<-dem < 200.
middem<-dem >= 200. & dem <= 3000.
highdem<-dem > 3000.
Xgrid = lowdem*X[1] + middem*(X[1]+(dem-200.)*((X[2]-X[1])/2800.))+ highdem*X[2]
Ygrid = lowdem*Y[1] + middem*(Y[1]+(dem-200.)*((Y[2]-Y[1])/2800.))+ highdem*Y[2]
Zgrid = lowdem*Z[1] + middem*(Z[1]+(dem-200.)*((Z[2]-Z[1])/2800.))+ highdem*Z[2]

# Read MERRA parameter file to establish static parameters for wind, surface temperature,
# surface pressure and relatice humidity data
geofile=file.path(GeoSourceDir,"merra_global_slv.nc")
if(!file.exists(geofile)){print(paste("Geo file missing:",geofile))}
f1<-open.ncdf(geofile)

phisvar = "MERRA300.prod.assim.const_2d_slv_Nx_fliped"

# Instantiate Julian and HRI functions
#source(file.path(mydir,"Rcode","SolarZenith.R"))

# Instantiate functions developed to subset netcdf files and to create rasters
#source(file.path(mydir,"Rcode","ncarray.R"))

tt=ncarray(f1,dem,phisvar,xdimname="lon",ydimname="lat")
#plot.surface(tt)
#lines(watershed$X,watershed$Y)
merralat<-tt$y
merralon<-tt$x
phis = tt$z

# Where in the nc file is the DEM
elev=phis/grav

#  Put in raster format for resampling using inputs from ncfile to avoid misalignment errors
Elevras = rasterFromArray(x=merralon,y=merralat,z=elev,crs=crs)

# Project which also resamples to model DEM 
refelev <- projectRaster(Elevras,dem)
DeltaZ = dem - refelev

close.ncdf(f1)

# Read MERRA parameter file to establish static parameters for short- and long-wave radiation data

geofile=file.path(GeoSourceDir,"merra_global_chm.nc")
if(!file.exists(geofile)){print(paste("Geo file missing:",geofile))}

f1<-open.ncdf(geofile)

# Use function
phisvar="MERRA300.prod.assim.const_2d_chm_Fx_fliped"
tt=ncarray(f1,dem,phisvar,xdimname="lon",ydimname="lat")
merrardlat=tt$y
merrardlon=tt$x
elevrd = tt$z
elevrd=elevrd/grav
Elevrasrd = rasterFromArray(x=merrardlon,y=merrardlat,z=elevrd ,crs=crs)
#plot(Elevrasrd ,zlim=c(0,8000),col=tim.colors())
#lines(watershed$X,watershed$Y)

# dem prep 
refelevrd <- projectRaster(Elevrasrd,dem)
#plot(refelevrd ,zlim=c(0,8000),col=tim.colors())
#lines(watershed$X,watershed$Y)
DeltaZrd = dem - refelevrd
close.ncdf(f1)

prevmonth = 0 
i = 1
firstValue = firstValue + 0

#First file
# loop through MERRA files to extract data; I'm keeping this as a
# separate loop to limit the number of open.ncdf commands
while(i <= ndays) {
   currentDate = dateList[i]
   print(paste("Downscaling for date:",currentDate))
   dname = format(currentDate,"%d")
   mname = format(currentDate,"%m")
   yname = format(currentDate,"%Y")
   day = as.integer(dname)
   month = as.integer(mname)
   year = as.integer(yname)

   # input filename
   fname = file.path(diri,paste("merra.prod.assim.",yname,mname,dname,".nc",sep=""))
   fnamerd = file.path(diri,paste("merra.prod.rad.",yname,mname,dname,".nc",sep=""))
   if(PrecipDataSource=="RFE")
   {
   fnamerfe = file.path(dirir,paste("rfe2.sa.",yname,mname,dname,".nc",sep=""))
   }
   if(PrecipDataSource=="MERRA")
   {
   fnamerfe= file.path(diri,paste("merra.prod.precip.",yname,mname,dname,".nc",sep=""))
   }
   if(prevmonth != month) {
      if(i > 1) 
      {
         # close netcdf file
         close.ncdf(ncid)
      }
      # output filename
      fnm = file.path(diro,paste("merra.rfe.90m.",yname,mname,".nc",sep=''))

      # define netcdf dimension variables
      lon_vals = seq(x1min,x1max,demres[1])
      lat_vals = seq(y1min,y1max,demres[2])
      nd <- monthDays(as.character(currentDate))
      sDate = as.Date(paste(paste(yname,mname,dname,sep="-")))
      eDate = as.Date(paste(paste(yname,mname,nd,sep="-")))
      edname = format(as.Date(enddate),"%d")
      if(eDate > enddate) {
         eDate = enddate
      }

      dateListMonth = seq(as.Date(sDate),as.Date(eDate),by="1 day")
      nDaysMonth = c(length(dateListMonth))
      timeStepList = seq(0,21,3)
      numVals = nstep*nDaysMonth
      monthHours = firstValue+(nDaysMonth*24)-dt
      time_vals = seq(firstValue,monthHours,dt)
      firstValue = time_vals[length(time_vals)]+dt
      newCounterTime = 1

      new_nx = length(lon_vals)
      new_ny = length(lat_vals)
      new_nz = length(time_vals)

      time_unit = paste("hours since ",referencedate," 00:00:00",sep='')

      # Create a new netcdf file to store all variables/time step for one day in only one file
      dim_lon = dim.def.ncdf(xcoordinate,xcoordinateunit,lon_vals,unlim=FALSE) 
      dim_lat = dim.def.ncdf(ycoordinate,ycoordinateunit,lat_vals,unlim=FALSE) 
      dim_time = dim.def.ncdf("time",time_unit,as.double(time_vals),unlim=FALSE) 

      # missing value
      mv = -9999

      # define variables
      var_temp = var.def.ncdf('tair','degree Celcius',list(dim_lon,dim_lat,dim_time),mv,longname='air temperature at 2m')
      var_rh = var.def.ncdf('rh','fraction',list(dim_lon,dim_lat,dim_time),mv,longname='relative humidity at 2m')
      var_wind = var.def.ncdf('wind','m s-1',list(dim_lon,dim_lat,dim_time),mv,longname='wind at 2m')
      var_swin = var.def.ncdf('swin','kJ m-2 h-1',list(dim_lon,dim_lat,dim_time),mv,longname='surface downward shortwave flux')
      #var_swin_Beer = var.def.ncdf('swin_Beer','kJ m-2 h-1',list(dim_lon,dim_lat,dim_time),mv,longname='surface downward shortwave flux using Beers Law')
      var_lwin = var.def.ncdf('lwin','kJ m-2 h-1',list(dim_lon,dim_lat,dim_time),mv,longname='surface downward longwave flux')
      var_prec = var.def.ncdf('rain','m hr-1',list(dim_lon,dim_lat,dim_time),mv,longname='precipitation')


      # Create NCDF file  
      ncid = create.ncdf(fnm,list(var_temp,var_rh,var_wind,var_swin,var_lwin,var_prec))
      print("Created netcdf file")
      print(fnm)
   }
   prevmonth = month

   if(file.exists(fname) && file.exists(fnamerd) && file.exists(fnamerfe)) {
      # reading wind, pressure, temperature and relative humidity
      f1<-open.ncdf(fname)
      j = 1
      while(j <= nstep) {
        tt=ncarray(f1,dem,"temp",j)
        if(j==1)
         {
          #First pass through get lat and long dimensions      
          merralat_df<-tt$y
          merralon_df<-tt$x
          nx=length(merralon_df)
          ny=length(merralat_df)
          temp = array(0,c(nx,ny,nstep))
          q = array(0,c(nx,ny,nstep))
          pres = array(0,c(nx,ny,nstep))
          wind = array(0,c(nx,ny,nstep))
          windu = array(0,c(nx,ny,nstep))
          windv = array(0,c(nx,ny,nstep))
         }
         temp[,,j] = tt$z
         pres[,,j] = ncarray(f1,dem,"ps",j)$z
         q[,,j] = ncarray(f1,dem,"qv2m",j)$z
         wind[,,j] = ncarray(f1,dem,"wind",j)$z
         windu[,,j] = ncarray(f1,dem,"windu",j)$z
         windv[,,j] = ncarray(f1,dem,"windv",j)$z
         j = j+1
      }
      close.ncdf(f1)

      is.na(temp) = temp <= -60 # this validation is for oC
      is.na(pres) = pres < 0 
      is.na(q) = q < 0
      is.na(wind) = wind < 0

      # reading long- and short-wave radiation
      f1<-open.ncdf(fnamerd)

      j = 1
      while(j <= nstep) {
        tt=ncarray(f1,dem,"swrad",j)
        if(j==1){
     # First pass through get lat and long dimensions      
          merralatrd_df<-tt$y
          merralonrd_df<-tt$x
          nxrd=length(merralonrd_df)
          nyrd=length(merralatrd_df)
          SWin = array(0,c(nxrd,nyrd,nstep))
          LWin = array(0,c(nxrd,nyrd,nstep))
         }
         SWin[,,j] =tt$z
         LWin[,,j] = ncarray(f1,dem,"lwrad",j)$z
         j = j+1
      }
      close.ncdf(f1)

      # reading RFE (precipitation) data
     if(PrecipDataSource=="RFE")
     {
      f1<-open.ncdf(fnamerfe)
      j = 1
      while(j <= nstep) {
        tt=ncarray(f1,dem,"rain",j,xdimname="lon",ydimname="lat")
        if(j==1){
      #   First pass through get lat and long dimensions      
          rfelat_df<-tt$y
          rfelon_df<-tt$x
          #minrfelat_df<-min(rfelat_df,na.rm=TRUE)
          #minrfelon_df<-min(rfelon_df,na.rm=TRUE)
          nxrfe=length(rfelon_df)
          nyrfe=length(rfelat_df)
          prec = array(0,c(nxrfe,nyrfe,nstep))
         }
         prec[,,j] = tt$z
         j = j+1
      }
      close.ncdf(f1)
      }

     # reading MERRA precipitation data
     if(PrecipDataSource=="MERRA")
     {
      f1<-open.ncdf(fnamerfe)
      j = 1
      while(j <= nstep) {
        tt=ncarray(f1,dem,"rain",j)
        if(j==1)
         {
          #First pass through get lat and long dimensions      
          rfelat_df<-tt$y
          rfelon_df<-tt$x
          nxrfe=length(rfelon_df)
          nyrfe=length(rfelat_df)
          prec = array(0,c(nxrfe,nyrfe,nstep))
         }
         prec[,,j] = tt$z
         j = j+1
      }
      close.ncdf(f1)
     }
      is.na(SWin) = SWin < 0  # -9000
      #is.na(SWin) = SWin > 1367*1.05  # greater than 1.05 of solar constant
      is.na(LWin) = LWin <= 0 # -9000
      is.na(prec) = prec < 0

      # variable manipulations
      e = pres*q/(0.378*q+0.622)  # units = Pascals
      Td = c*log(e/a, base=exp(1))/(b-log(e/a, base=exp(1)))

      j = 1
      while(j <= nstep) {
         doy<-mdy.date(month,day,year) - mdy.date(1,1,year) + 1 
      
         # AIR TEMPERATURE
         Tras=rasterFromArray(x=merralon_df,y=merralat_df,z=temp[,,j] ,crs=crs)
         Tint <- projectRaster(Tras,dem)
         Tgrid=Tint-DeltaZ*lapse[month]         
         outvector=as.vector(t(as.matrix(Tgrid)[new_ny:1,]))  
         # Flip y values.  Then transpose to account for by column internal to R whearas ncdf takes by row
         put.var.ncdf(ncid,var_temp,outvector,start=c(1,1,newCounterTime),count=c(new_nx,new_ny,1))
         sync.ncdf(ncid)

         # SURFACE PRESSURE
         #print("pressure")
         Pras=rasterFromArray(x=merralon_df,y=merralat_df,z=pres[,,j],crs=crs)
         Pint <- projectRaster(Pras,dem)
         Pgrid = Pint*exp(DeltaZ/(-H))

         # RELATIVE HUMIDITY
         Tdras=rasterFromArray(x=merralon_df,y=merralat_df,z=Td[,,j] ,crs=crs)
         Tdint <- projectRaster(Tdras,dem)
         Tdgrid = Tdint - DeltaZ*c*tdlapse[month]/b
         es_grid = a*exp(Tgrid*b/(c+Tgrid))
         e_grid = a*exp(Tdgrid*b/(c+Tdgrid))
         RHgrid = min(e_grid/es_grid,1.0)
         outvector=as.vector(t(as.matrix(RHgrid)[new_ny:1,]))  
         # Flip y values.  Then transpose to account for by column internal to R whearas ncdf takes by row
         put.var.ncdf(ncid,var_rh,outvector,start=c(1,1,newCounterTime),count=c(new_nx,new_ny,1))
         sync.ncdf(ncid)

         # WIND SPEED
         #------------------------------------------
         #print("wind")
         #Windras=rasterFromArray(x=merralon_df,y=merralat_df,z=wind[,,j] ,crs=crs)
         #Windras=raster(t(wind[,ny:1,j]),xmn=x1min,xmx=x1max,ymn=y1min,ymx=y1max,crs=crs)
         #Windint <- projectRaster(Windras,dem)
         #plot(Windint)
         #------------------------------------------

         Windrasu=rasterFromArray(x=merralon_df,y=merralat_df,z=windu[,,j] ,crs=crs)
         Windrasv=rasterFromArray(x=merralon_df,y=merralat_df,z=windv[,,j] ,crs=crs)
         Windintu <- projectRaster(Windrasu,dem)
         Windintv <- projectRaster(Windrasv,dem)
         windsqt=sqrt(Windintu^2+Windintv^2)
         gamac = (agg5wind-(cellStats(agg5wind,'max')+cellStats(agg5wind,'min'))/2)/
             (cellStats(agg5wind,'max')-cellStats(agg5wind,'min'))
         #-----------------------------------------------------------------------
         
         thetawind=3*pi/2-atan(Windintv/Windintu)
         sloperes=slopewind*cos(thetawind-aspectwind)
         gamas=(sloperes-(cellStats(sloperes,'max')+cellStats(sloperes,'min'))/2)/
         (cellStats(sloperes,'max')-cellStats(sloperes,'min'))
         shaic=0.5
         shais=0.5
         windwt=1+gamac*shaic+shais*gamas
         Windint=windsqt*windwt
         #-----------------------------------------------------------------------
         Windgrid = (Windint + abs(Windint))/2  # if negative make 0   
         outvector=as.vector(t(as.matrix(Windgrid )[new_ny:1,]))  
         # Flip y values.  Then transpose to account for by column internal to R whearas ncdf takes by row
         put.var.ncdf(ncid,var_wind,outvector,start=c(1,1,newCounterTime),count=c(new_nx,new_ny,1))
         sync.ncdf(ncid)
         #------------------------------------------
         #------------------------------------------

         # SOLAR RADIATION
         #------------------------------------------
         lapse_lambda=-0.0065 #K/m
         #initialize
         Qsi_Beer = Tgrid *0.0
         SWras=rasterFromArray(x=merralonrd_df,y=merralatrd_df,z=SWin[,,j],crs=crs)
         SWint <- projectRaster(SWras,dem)   
         # Raster center in projected coordinates
         meanx=(xmin(SWint)+xmax(SWint))/2.
         meany=(ymin(SWint)+ymax(SWint))/2.
         midpnt=SpatialPoints(coords=matrix(c(meanx,meany),nrow=1),proj4string=CRS(projection(dem)))
         midlatlon=spTransform(midpnt,CRS(crs))
         latavg=as.numeric(midlatlon@coords[1,2])
         lonavg=as.numeric(midlatlon@coords[1,1])
         localhr = (timeStepList[j]+lonavg*24/360)%%24
         zenith = HRIfunc(year=as.integer(year),month=as.integer(month),day=as.integer(day),localhr,dt,lat=latavg,slope=0,azi=0)
         k=Qsi_Beer*0
         if (zenith[[2]] <= 89.0)
         {
          psi=SWint/((Sstar*3.6)*cos(pi*zenith[[2]]/180)) #psi=Psi_dir+Psi_diff
          Tgrid_K=Tgrid+273.15
          Pgrid_est=Pgrid*((Tgrid_K + DeltaZ*lapse_lambda)/Tgrid_K)^(-grav/(rdry*lapse_lambda))
          k = -1*log(psi)/Pgrid
          #k[is.infinite(k)]=0
          #Qsi_Beer=SWint*exp(-1*k*Pgrid_est)
          Qsi_Beer=((Sstar*3.6)*cos(pi*zenith[[2]]/180))*exp(-1*k*Pgrid_est)
           }else{
           Qsi_Beer=0*Qsi_Beer
          }
         outvector=as.vector(t(as.matrix(Qsi_Beer)[new_ny:1,]))  
         #outvector[outvector > 4923.72*1.05]=NA 
         #Flip y values.  Then transpose to account for by column internal to R whearas ncdf takes by row
         put.var.ncdf(ncid,var_swin,outvector,start=c(1,1,newCounterTime),count=c(new_nx,new_ny,1))
         sync.ncdf(ncid)
         
         # INCOMING LONGWAVE RADIATION
         elev700 = (70000-Pint)/(1.29*grav)
         T700 = Tint-lapse[month]*elev700
         Td700 = Tdint - elev700*c*tdlapse[month]/b
         e700 = a*exp(b*Td700/(c+Td700))
         es700 = a*exp(b*T700/(c+T700))
         rh700 = 100*e700/es700
         sigma_c = min(0.832*exp((rh700-100)/41.6),1.0)
         emiss = 1.083*(1+Zgrid*sigma_c^2)*(1-Xgrid*exp(-1*Ygrid*e_grid/(Tgrid[]+273.15)))
         LW = emiss*sigma*(Tgrid[]+273.15)^4
         outvector=as.vector(t(as.matrix(LW)[new_ny:1,]))  
         # Flip y values.  Then transpose to account for by column internal to R whearas ncdf takes by row
         put.var.ncdf(ncid,var_lwin,outvector,start=c(1,1,newCounterTime),count=c(new_nx,new_ny,1))
        sync.ncdf(ncid)
         
         # PRECIPITATION
         precras=rasterFromArray(x=rfelon_df,y=rfelat_df,z=prec[,,j],crs=crs)
         prectmp <- getValues(precras)
         if(all(is.na(prectmp))) {
            outvector = array(mv,c(new_nx,new_ny,1))
         }else 
         {
            precint <- projectRaster(precras,dem)*((1+adjfactor[month]*DeltaZ)/(1-adjfactor[month]*DeltaZ))
            outvector=as.vector(t(as.matrix(precint)[new_ny:1,]))  
            # Flip y values. Then transpose to account for by column internal to R whearas ncdf takes by row
         }
         put.var.ncdf(ncid,var_prec,outvector,start=c(1,1,newCounterTime),count=c(new_nx,new_ny,1))
         sync.ncdf(ncid)

         j = j + 1
         newCounterTime = newCounterTime + 1
      }
   }else
      print("file does not exist")
   i=i+1
}
close.ncdf(ncid)


#--------------------------------------------------------------
#--------------------------------------------------------------
#--------------------------------------------------------------
#--------------------------------------------------------------
#--------------------------------------------------------------
write("Coordinates of center of the domain",file=latlonfile)
write(paste("Latitude:",round(latavg,digits=3)),file=latlonfile,append=T)
write(paste("Longitude:",round(lonavg,digits=3)),file=latlonfile,append=T)

print("Coordinates of center of the domain")
print(paste("Latitude:",round(latavg,digits=3)))
print(paste("Longitude:",round(lonavg,digits=3)))
print(latlonfile)

CodeEndTime=as.numeric(Sys.time())
DownScaleTime=CodeEndTime-CodeStartTime
CodeStartTime=CodeEndTime

#--------------------------------------------------------------
#--------------------------------------------------------------
#--------------------------------------------------------------
#--------------------------------------------------------------
#--------------------------------------------------------------
####################################
curdir=getwd()
demfile=DEMfile
stationfile=svalue(tbl41[2,2])
workdir2=file.path(svalue(tbl4[3,2]),"tfirst")
workdir=gsub("[\\]","/",workdir2) 
files=list.files(workdir,pattern="merra")
outputfilename=files
readfilenames=files
if(file.exists(demfile) && file.exists(workdir) && file.exists(stationfile)) {
biasfactor=unlist(precipitationBIASfactorfunction(demfile,stationfile,workdir))
setwd(workdir)
library(stringr)
for (i in 1:length(files))
{
strs=unlist(strsplit(files[i], split="[.]"))
monthnum=as.numeric(str_sub(strs[4], start= -2))
txt=paste("rain=rain*",(1-biasfactor[monthnum]),sep="")
p1=paste("cdo expr,\'tair=tair;rh=rh;wind=wind;swin=swin;lwin=lwin;",txt,"\'",sep="")
tempfile=paste("temp","_",files[i],sep="")
cmd=paste(p1,files[i],tempfile, sep=" ")
system(cmd)
unlink(files[i])
}
}
####################################
#--------------------------------------------------------------
#--------------------------------------------------------------
#--------------------------------------------------------------
#--------------------------------------------------------------
#--------------------------------------------------------------
####################################

curdir=getwd()
setwd(DestDir)
inputdir="./tfirst"
files=list.files(inputdir)
for (i in 1:length(files))
{
  if(i == 1) # on first pass through figure out dimension names
{
  inputfile=paste("tfirst//",files[i],sep="")
  outputfile=outputfilename[i]
  ncid = open.ncdf(inputfile)
  xname=ncid$dim[[1]]$name
  yname=ncid$dim[[2]]$name
  close.ncdf(ncid)
  # NCO command to reorder the dimensions of the variables  
  #ncpdq -a longitude,latitude,time merra.rfe.90m.200210.nc -O merra.rfe.90m.200110NCO.nc
  cmd=paste('ncpdq -a ',xname,',',yname,',','time ',inputfile,' -O ',outputfile,sep='')
  system(cmd)
}else
{
  inputfile=paste("tfirst//",files[i],sep="")
  outputfile=outputfilename[i]
  # NCO command to reorder the dimensions of the variables  
  #ncpdq -a longitude,latitude,time merra.rfe.90m.200210.nc -O merra.rfe.90m.200110NCO.nc
  cmd=paste('ncpdq -a ',xname,',',yname,',','time ',inputfile,' -O ',outputfile,sep='')
  #cmd=paste('ncpdq -a longitude,latitude,time',inputfile,'-O',outputfile,sep=' ')
  system(cmd)
}
}

unlink(c("temp","regrid","tfirst"),recursive=T)
#setwd(curdir)

# Clean up
#curdir=getwd()
#setwd(file.path(DestDir))
#unlink(c("temp","regrid"),recursive=T)

CodeEndTime=as.numeric(Sys.time())
ReOrderTime=CodeEndTime-CodeStartTime

print(paste("Regrid time:",round(RegridTime,digits=2),"sec"))
print(paste("Downscale time:",round(DownScaleTime,digits=2),"sec"))
print(paste("ReOrderTime time:",round(ReOrderTime,digits=2),"sec"))
}

#-------------------------------------------------------------------
#-------------------------------------------------------------------
#-------------------------------------------------------------------
#-------------------------------------------------------------------
#-------------------------------------------------------------------
#-------------------------------------------------------------------
#-------------------------------------------------------------------
#-------------------------------------------------------------------


MERRAAvgMonthlyTrange=function(trangefile,demfile,outfolder)
{
#outfolder="C:\\Users\\avirup\\Desktop\\TryRGUI\\Outputs"
curdir=getwd()
dir.create(outfolder,showWarnings = FALSE,recursive = TRUE)
setwd(outfolder)
# --------------------------------------------------
# Load needed netCDF libraries
# --------------------------------------------------
library(rgdal)
library(raster)
library(ncdf)
library(fields)
#---------------------------------------------------
# Extracting and summarizing the grid of DEM for the output data
dem=raster(readGDAL(demfile))
demproj=projection(dem)
#---------------------------------------------------
#form the names of the netCDF files
monthlist=1:12
ncdffiles=paste("merra.trange.climatology.",formatC(monthlist,width=2,flag='0'),".nc",sep='')
missval=-9999

#  Set defaults
proj='longlat'
xcoordinate='longitude'
ycoordinate='latitude'
xcoordinateunit='degrees_east'
ycoordinateunit='degrees_north'
crs = "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0"

#  Define a trim function
trim <- function (x) gsub("^//s+|//s+$", "", x)

# Now look at info in dem projection detail and adjust units if needed
demprojparts=unlist(strsplit(demproj, "[+]"))
for(str in demprojparts)
{
sub12=unlist(strsplit(str, "[=]"))
if(length(sub12)>1)if(sub12[1]=="proj")proj=trim(sub12[2])
if(length(sub12)>1)if(sub12[1]=="units")
{ xcoordinateunit=trim(sub12[2])
  ycoordinateunit=xcoordinateunit
}
}
if(proj!='longlat')  # if projection not long lat write with x and y as dimension names
{
 xcoordinate='x'
 ycoordinate='y'
}


ncobj=open.ncdf(trangefile)
varname=ncobj$var[[1]]$name
varunit=ncobj$var[[1]]$units
varlongname=ncobj$var[[1]]$longname

for(i in 1:12)
{
 tt=ncarray(ncobj,dem,varname,timeindex=i,xdimname='longitude',ydimname='latitude')
 #plot.surface(tt)
 merralon = tt$x
 merralat = tt$y
 tvals=tt$z
 tempgrid=rasterFromArray(x=merralon,y=merralat,z=tvals,crs=crs)
 tempfileprojected = projectRaster(tempgrid,dem)
 #plot(tempfileprojected) 
 print(paste(month.name[i],"Averge diurnal temperature range over sub domain:",
 round(mean(as.matrix(tempfileprojected)),digits=2)))
 # DEM Extent
 demres=res(tempfileprojected)
 x1min=xmin(tempfileprojected)
 x1max=xmax(tempfileprojected)
 y1min=ymin(tempfileprojected)
 y1max=ymax(tempfileprojected)
 x_vals = seq(x1min+demres[1]/2,x1max-demres[1]/2,demres[1])
 y_vals = seq(y1min+demres[2]/2,y1max-demres[2]/2,demres[2])
 # define the netcdf coordinate variables -- note these have values!
 dim1 = dim.def.ncdf(xcoordinate,xcoordinateunit, as.double(x_vals))
 dim2 = dim.def.ncdf(ycoordinate,ycoordinateunit, as.double(y_vals))
 #------------------------------------------------------------------
 outvector=as.vector(t(as.matrix(tempfileprojected)[length(y_vals):1,]))  
 # Flip y values.  Then transpose to account for by column internal to R whearas ncdf takes by row
 var = var.def.ncdf(varname,varunit,list(dim1,dim2),missval,varlongname)
 nc.ex = create.ncdf(ncdffiles[i],var)
 put.var.ncdf(nc.ex,var,outvector)
 close.ncdf(nc.ex)
}
setwd(curdir)
}

#-------------------------------------------------------------------
#-------------------------------------------------------------------
#-------------------------------------------------------------------
#-------------------------------------------------------------------
#-------------------------------------------------------------------
#-------------------------------------------------------------------
#-------------------------------------------------------------------
#-------------------------------------------------------------------
# Check if a date is within a range
checkrangedate=function(start_date, end_date, date_from_user)
{
  return ((date_from_user >= start_date) && (date_from_user <= end_date));
}

#--------------------------------------------------
#--------------------------------------------------
# editing the gfilebrowse function
gfilebrowse=function (text = "", type = c("open", "selectdir"),  
    quote = FALSE, container = NULL, ..., toolkit = guiToolkit()) 
{
    widget <- .gfilebrowse(toolkit, text = text, type = type, 
        quote = quote, container = container, ...)
    obj <- new("gFilebrowse", widget = widget, toolkit = toolkit)
    return(obj)
}
#--------------------------------------------------
#--------------------------------------------------
# editing the gfilebrowse function
gcalendar2=function (text = "", format = "%Y-%m-%d", handler = NULL, action = NULL, initial.msg="",
    container = NULL, ..., toolkit = guiToolkit()) 
{
    widget <- .gcalendar(toolkit, text = text, format = format, 
        handler = handler, action = action, container = container, initial.msg=initial.msg,
        ...)
    obj <- new("gCalendar", widget = widget, toolkit = toolkit)
    return(obj)
}
#--------------------------------------------------
#--------------------------------------------------
# Creating a message window
confirmDialog <- function(text,message, handler=NULL,img="info") {
window <- gwindow(text,height=100,width=50)
group <- ggroup(container = window)
gimage(filename = img, dirname="stock", size="dialog", container=group)
## A group for the message and buttons
inner.group <- ggroup(horizontal=FALSE, container = group)
glabel(message, container=inner.group, expand=TRUE)
## A group to organize the buttons
button.group <- ggroup(container = inner.group)
## Push buttons to right
addSpring(button.group)
gbutton("ok", container=button.group, handler = function(h,...) dispose(window))
return()
}










#--------------------------------------------------
#--------------------------------------------------
# Graphics interface work starts here--------------
#--------------------------------------------------
#--------------------------------------------------

#get the images from web, use them and delete them

library(RgoogleMaps)

link1 = "http://hydrology.usu.edu/snow/uebgrid/Picture/Logo.png"
link2="http://hydrology.usu.edu/snow/uebgrid/Picture/SouthAsiaMap.png"
img1=basename(link1)
img2=basename(link2)
download.file(link1,img1, mode = 'wb')
download.file(link2,img2, mode = 'wb')
imgdir=getwd()

win = gwindow("MERRA Spatial Downscaling for Hydrology (MSDH)",horizontal=FALSE)
group <- ggroup(container = win,horizontal=FALSE)
gimage(label="",filename=file.path(imgdir,img1),markup=FALSE,
       cont = group, horizontal=FALSE,width=200,height=100)
nb = gnotebook(cont=group,expand=TRUE)
#--------------------------------------------------
#--------------------------------------------------
#--------------------------------------------------
#WINDOW 1------------------------------------------
g <- ggroup(horizontal=FALSE, cont=nb, label = "Download Data")
#addSpace(g, 10, horizontal=FALSE)
tbl1 <- glayout(container=g,horizontal=FALSE,expand=TRUE,spacing = 10)
tbl1[1,1] <- glabel("Start Date")
font(tbl1[1,1]) <- c(family="times",size="large",style="italic")
tbl1[1,2] <- gcalendar(initial.msg = "YYYY-MM-DD", format = "%Y-%m-%d")
font(tbl1[1,2]) <- c(family="times",size="large",style="italic")
tbl1[1,3] <- glabel("End Date")
font(tbl1[1,3]) <- c(family="times",size="large",style="italic")
tbl1[1,4] <-  gcalendar(initial.msg = "YYYY-MM-DD", format = "%Y-%m-%d")
font(tbl1[1,4]) <- c(family="times",size="large",style="italic")
tbl1[2,1] <- glabel("")
#Valid entry dates: 2001-05-01 to present
font(tbl1[2,1]) <- c(family="times",size="large",style="italic")
addHandlerChanged(tbl1[1,2], function(h,...) 
{
svalue(tbl5[1,2])=svalue(tbl1[1,2])
})
addHandlerChanged(tbl1[1,4], function(h,...) 
{
svalue(tbl5[1,4])=svalue(tbl1[1,4])
})

tbl2 <- glayout(container=g,horizontal=FALSE,expand=TRUE)
tbl2[2,1] <- glabel("Select Spatial Domain Here

")
#valid Longitude: between 70 and 110 degree 
#valid Latitude: between 05 and 35 degree

font(tbl2[2,1]) <- c(family="times",size="large",style="italic")
tbl2[1,3]=gedit(initial.msg = "Lat North", width = 15,coerce.with=as.numeric)
tbl2[2,2]=gedit(initial.msg = "Lon West", width = 15,coerce.with=as.numeric)
tbl2[2,4]=gedit(initial.msg = "Lon East", width = 15,coerce.with=as.numeric)
tbl2[3,3]=gedit(initial.msg = "Lat South", width = 15,coerce.with=as.numeric)
tbl2[2,3]=gimage(label="",filename=file.path(imgdir,img2),markup=FALSE)
#deleting the images
unlink(c(file.path(imgdir,img1),file.path(imgdir,img2)))
tbl2[3,4]=gbutton("Update Map")
addHandlerChanged(tbl2[3,4], function(h,...) {
RFEdom=c(90,-90,180,-180)
latn=svalue(tbl2[1,3])
lats=svalue(tbl2[3,3])
lone=svalue(tbl2[2,4])
lonw=svalue(tbl2[2,2])
testnum <- function(x) is.numeric(x) & !is.na(x)
bound=testnum(c(latn,lats,lone,lonw))
bound2=c(checkrangedate(RFEdom[2],RFEdom[1],as.numeric(latn)),
         checkrangedate(RFEdom[2],RFEdom[1],as.numeric(lats)),
         checkrangedate(RFEdom[4],RFEdom[3],as.numeric(lone)),
         checkrangedate(RFEdom[4],RFEdom[3],as.numeric(lonw)))
boundchk3=!(bound&bound2)

boundchk1=!(as.numeric(lone)>as.numeric(lonw))
boundchk2=!(as.numeric(latn)>as.numeric(lats))
bounds=c(boundchk3,boundchk1,boundchk2)
msgs=c("Latitude North: must be number between -90 and 90",
       "Latitude South: must be number between -90 and 90",
       "Longitude East: must be number between -180 and 180",
       "Longitude West: must be number between -180 and 180",
       "Longitude East must be higher than Longitude West",
       "Latitude North must be higher than Latitude South")
msgshow=msgs[bounds|is.na(bounds)]
if(length(msgshow)>0)
{
confirmDialog(text="General Info!",msgshow)
}else{
lat <- c(lats,latn)
lon <- c(lonw,lone)
lat=lat-diff(lat)/10
lon=lon-diff(lon)/10
zoom <- 4
destfile=file.path(imgdir,"South_Asia.png")
terrmap <- GetMap.bbox(lon,lat, size=c(640,640), maptype="satellite", 
zoom=max(MaxZoom(range(lat), range(lon))),
destfile=destfile)
png(filename=destfile, width=0.8, height=0.8,units = "in",pointsize=12, bg="transparent", res=240)
PlotOnStaticMap(terrmap)
dev.off()
svalue(tbl2[2,3])=destfile
unlink(c(destfile,paste(destfile,".rda",sep="")))

}
})

tblradio <- glayout(container=g,horizontal=FALSE,expand=TRUE)
#addSpace(g, 10, horizontal=FALSE)
tblradio[1,1] <- glabel("Precipitation Data Source")
font(tblradio[1,1]) <- c(family="times",size="large",style="italic")
tblradio[1,2] <- gradio(c("RFE","MERRA"),horizontal=TRUE)
font(tblradio[1,2]) <- c(family="times",size="large",style="italic")
#tblradio[2,1] <- glabel("MERRA is global data source")
#tblradio[3,1] <- glabel("RFE supports South Asia for This application")
#font(tblradio[2,1]) <- c(family="times",size="large",style="italic")
#font(tblradio[3,1]) <- c(family="times",size="large",style="italic")


tbl3 <- glayout(container=g,horizontal=FALSE,expand=TRUE)
#addSpace(g, 10, horizontal=FALSE)
text="Destination Directory"
#addSpace(g, 10, horizontal=FALSE)
tbl3[1,1] <- glabel("Or Upload Your Shapefile Here")
font(tbl3[1,1]) <- c(family="times",size="large",style="italic")
tbl3[1,2] <- gfilebrowse(width=70,pos=0,type="open",quote=FALSE)
font(tbl3[1,1]) <- c(family="times",size="large",style="italic")
tbl3[2,1] <- glabel(text)
font(tbl3[2,1]) <- c(family="times",size="large",style="italic")
tbl3[2,2] <- gfilebrowse(text=getwd(),width=70,pos=0,type="selectdir",quote=FALSE)
font(tbl3[2,2]) <- c(family="times",size="large",style="italic")
addHandlerChanged(tbl3[1,2], function(h,...) {
library(PBSmapping)
require(RgoogleMaps)
require(rgdal)  
require(PBSmapping) 
shpspl=unlist(strsplit(svalue(tbl3[1,2]),"[.]"))
lenspl=length(shpspl)
if(!shpspl[lenspl]=="shp"){shpmsg="requires a shapefile"}
shpchr=nchar(svalue(tbl3[1,2]), type = "chars", allowNA = FALSE)
polyShpFile=substr(x = svalue(tbl3[1,2]), start = 1, stop = (shpchr-4))
shpPolySet=importShapefile(polyShpFile,projection="LL")
latsr = round(range(shpPolySet$Y), 4)
lonsr = round(range(shpPolySet$X), 4)
centerPt = c(mean(latsr),mean(lonsr))
boundBox = qbbox(lat = latsr,lon = lonsr)
destfile=file.path(imgdir,"South_Asia.png") 
png(filename=destfile, width=0.8, height=0.8,units = "in",pointsize=12, bg="transparent", res=240)
mapFromBbox = GetMap.bbox(boundBox$lonR,boundBox$latR,maptype="terrain")
PlotPolysOnStaticMap(mapFromBbox, shpPolySet, col=rgb(1, 1, 1, alpha=0.3),
 lwd=2, border="red",add = F)
dev.off() # end of example
svalue(tbl2[2,3])=destfile
svalue(tbl2[1,3])=latsr[2]
svalue(tbl2[3,3])=latsr[1]
svalue(tbl2[2,4])=lonsr[2]
svalue(tbl2[2,2])=lonsr[1]
})
addHandlerChanged(tbl3[2,2], function(h,...) 
{
svalue(tbl[1,2])=paste(svalue(tbl3[2,2]),"\\MODIS","\\MOD12Q1_LandCoverType_IGBP_SouthAsia_Geographic.tif",sep="")
svalue(tbl[3,2])=paste(svalue(tbl3[2,2]),"\\MODIS","\\landcover.csv",sep="")
svalue(tbl6[1,2])=paste(svalue(tbl3[2,2]),"\\MERRATrange","\\MERRAAvgMonthlyTrange.nc",sep="")
svalue(tbl4[1,2])=svalue(tbl3[2,2])
svalue(tbl4[2,2])=paste(svalue(tbl3[2,2]),"\\MERRAgeo",sep="")
svalue(tbl4[3,2])=paste(svalue(tbl3[2,2]),"\\Outputs",sep="")
svalue(tbl4[5,2])=paste(svalue(tbl3[2,2]),"\\Outputs",sep="")
svalue(tbl[4,2])=paste(svalue(tbl3[2,2]),"\\Outputs",sep="")
svalue(tbl6[3,2])=paste(svalue(tbl3[2,2]),"\\Outputs",sep="")
})
group2 <- ggroup(container=g,horizontal=TRUE)
addSpring(group2)   
inforButton1 <- gbutton("Information", container=group2,handler=function(h,...){
                   confirmDialog(text="gerneral Info","This script will download MERRA and RFE data
						and its inputs are:
						1. Start and end of the downloading date
						2. latitude and longitde boundary")})
tblbtn1 <- gbutton("Download Data", container=group2)
#btn_quit <- gbutton(text = "  Close window  ",container = group2, handler=function(...) dispose(win))
addHandlerChanged(tblbtn1, function(h,...) {

if(svalue(tblradio[1,2])=="MERRA")
{
RFEdom=c(90,-90,180,-180)
latn=svalue(tbl2[1,3])
lats=svalue(tbl2[3,3])
lone=svalue(tbl2[2,4])
lonw=svalue(tbl2[2,2])
testnum <- function(x) is.numeric(x) & !is.na(x)
bound=testnum(c(latn,lats,lone,lonw))
bound2=c(checkrangedate(RFEdom[2],RFEdom[1],as.numeric(latn)),
         checkrangedate(RFEdom[2],RFEdom[1],as.numeric(lats)),
         checkrangedate(RFEdom[4],RFEdom[3],as.numeric(lone)),
         checkrangedate(RFEdom[4],RFEdom[3],as.numeric(lonw)))
boundchk3=!(bound&bound2)

boundchk1=!(as.numeric(lone)>as.numeric(lonw))
boundchk2=!(as.numeric(latn)>as.numeric(lats))
bounds=c(boundchk3,boundchk1,boundchk2)
msgs=c("Latitude North: must be number between -90 and 90",
       "Latitude South: must be number between -90 and 90",
       "Longitude East: must be number between -180 and 180",
       "Longitude West: must be number between -180 and 180",
       "Longitude East must be higher than Longitude West",
       "Latitude North must be higher than Latitude South")
msgshow=msgs[bounds|is.na(bounds)]
}

if(svalue(tblradio[1,2])=="RFE")
{
RFEdom=c(35,5,110,70)
latn=svalue(tbl2[1,3])
lats=svalue(tbl2[3,3])
lone=svalue(tbl2[2,4])
lonw=svalue(tbl2[2,2])
testnum <- function(x) is.numeric(x) & !is.na(x)
bound=testnum(c(latn,lats,lone,lonw))
bound2=c(checkrangedate(RFEdom[2],RFEdom[1],as.numeric(latn)),
         checkrangedate(RFEdom[2],RFEdom[1],as.numeric(lats)),
         checkrangedate(RFEdom[4],RFEdom[3],as.numeric(lone)),
         checkrangedate(RFEdom[4],RFEdom[3],as.numeric(lonw)))
boundchk3=!(bound&bound2)

boundchk1=!(as.numeric(lone)>as.numeric(lonw))
boundchk2=!(as.numeric(latn)>as.numeric(lats))
bounds=c(boundchk3,boundchk1,boundchk2)
msgs=c("For RFE2 Latitude North: must be number between 5 and 35",
       "For RFE2 Latitude South: must be number between 5 and 35",
       "For RFE2 Longitude East: must be number between 70 and 110",
       "For RFE2 Longitude West: must be number between 70 and 110",
       "Longitude East must be higher than Longitude West",
       "Latitude North must be higher than Latitude South")
msgshow=msgs[bounds|is.na(bounds)]
}

stdate=svalue(tbl1[1,2])
endate=svalue(tbl1[1,4])
RFEstdate=as.Date("2001-05-01",format= "%Y-%m-%d")
RFEendate=Sys.Date()-90
chk1=!checkrangedate(RFEstdate, RFEendate, stdate)
chk2=!checkrangedate(RFEstdate, RFEendate, endate)
chk3=!(endate>=stdate)
chks=c(chk1,chk2,chk3)
msgs2=c("For RFE2, Start date must be between 2001-05-01 and 90 days past from today",
        "For RFE2, End date must be between 2001-05-01 and 90 days past from today",
        "For MERRA, End date must be between 1979-01-01 and 90 days past from today",
        "For MERRA, End date must be between 1979-01-01 and 90 days past from today",                     "End Date must be higher than Start date")
msgshow2=msgs2[chks|is.na(chks)]

if(length(c(msgshow,msgshow2))>0)
{
confirmDialog(text="Error!",c(msgshow,msgshow2))
}else{
DestDir=svalue(tbl3[2,2])
latsouth=as.numeric(svalue(tbl2[3,3]))
latnorth=as.numeric(svalue(tbl2[1,3]))
lonwest=as.numeric(svalue(tbl2[2,2]))
loneast=as.numeric(svalue(tbl2[2,4]))
startdate=as.Date(svalue(tbl1[1,2]))
enddate=as.Date(svalue(tbl1[1,4]))
PrecipDataSource=svalue(tblradio[1,2])
list(DestDir,latsouth,latnorth,lonwest,loneast,startdate,enddate,PrecipDataSource)
do.call("downloadMERRARFE",list(DestDir,latsouth,latnorth,lonwest,loneast,startdate,enddate,PrecipDataSource))
confirmDialog(text="Success!","data downloading is successfully completed!",img="ok")
}
})


#WINDOW 2------------------------------------------
#--------------------------------------------------
#--------------------------------------------------
#texts to help the users
texts1=c("Calculate Temperature Lapse Rate",
        "Temperature Station information (csv file)")
texts2=c("Calculate Precipitation Adjustment Factor",
        "Precipiattion Station information (csv file)")
gcon <- ggroup(horizontal=FALSE, cont=nb, label = "Coefficient calculations")
tbl40 <- glayout(container=gcon,horizontal=FALSE,expand=TRUE)
## Layout widgets
for( i in 1:length(texts1)) {
tbl40[i,1] <- glabel(texts1[i])
font(tbl40[i,1]) <- c(family="times",size="large",style="italic")
}
tbl40[2,2] <- gfilebrowse(text="",width=80,pos=0,type="open",quote=FALSE)
lapsegroup2 <- ggroup(container=gcon,horizontal=TRUE)

lapseratebutton <- gbutton(container=lapsegroup2,"Calculate Temperature Lapse Rate",handler=function(h,...){ 
stationfile  <- svalue(tbl40[2,2])
lapserate=unlist(lasperatefunction(stationfile))
for( i in 1:6){
svalue(tbl71[1,i*2]) <- -lapserate[i]
svalue(tbl71[2,i*2]) <- -lapserate[6+i]
}
}
)

frame100 <- gframe(container=gcon,horizontal=TRUE,text="Monthly Lapse Rate (Degree Celsius/meter)")
tbl71 <- glayout(container=frame100,horizontal=TRUE,expand=TRUE)
monthnames=c("January","February","March","April","May","June","July",
             "August","September","October","November","December")
lapse=rep(NA,12)
for( i in 1:6){
tbl71[1,(i*2-1)] <- glabel(monthnames[i])
font(tbl71[1,(i*2-1)]) <- c(family="times",size="large",style="italic")
tbl71[1,i*2] <- gedit(width = 10, coerce.with = as.numeric)
tbl71[2,(i*2-1)] <- glabel(monthnames[6+i])
font(tbl71[2,(i*2-1)]) <- c(family="times",size="large",style="italic")
tbl71[2,i*2] <- gedit( width = 10, coerce.with = as.numeric)
}

tbl41 <- glayout(container=gcon,horizontal=FALSE,expand=TRUE)
for( i in 1:length(texts2)) {
tbl41[i,1] <- glabel(texts2[i])
font(tbl41[i,1]) <- c(family="times",size="large",style="italic")
}
tbl41[2,2] <- gfilebrowse(text="",width=80,pos=0,type="open",quote=FALSE)
lapsegroup3 <- ggroup(container=gcon,horizontal=TRUE)
lapseratebutton <- gbutton(container=lapsegroup3,"Calculate Precipitation Adjustment Factor",handler=function(h,...){ 
stationfile2  <- svalue(tbl41[2,2])
precipadj=unlist(precipitationadjfactorfunction(stationfile2))
for( i in 1:6){
svalue(tbl70[1,i*2]) <- precipadj[i]
svalue(tbl70[2,i*2]) <- precipadj[6+i]
}
}
)

frame200 <- gframe(container=gcon,horizontal=TRUE,text="Monthly Precipitation Adjustment Factor (m/meter)")
tbl70 <- glayout(container=frame200,horizontal=TRUE,expand=TRUE)
monthnames=c("January","February","March","April","May","June","July",
             "August","September","October","November","December")
lapse=rep(0,12)
for( i in 1:6){
tbl70[1,(i*2-1)] <- glabel(monthnames[i])
font(tbl70[1,(i*2-1)]) <- c(family="times",size="large",style="italic")
tbl70[1,i*2] <- gedit(width = 10, coerce.with = as.numeric)
tbl70[2,(i*2-1)] <- glabel(monthnames[6+i])
font(tbl70[2,(i*2-1)]) <- c(family="times",size="large",style="italic")
tbl70[2,i*2] <- gedit(width = 10, coerce.with = as.numeric)
}

#WINDOW 3------------------------------------------
#--------------------------------------------------
#--------------------------------------------------
#texts to help the users
texts=c("MERRA and RFE data Directory","Geo-potential File Directory",
        "Output Directory","DEM TIFF File (input)", "Lat-lon Output Directory (output)")
g <- ggroup(horizontal=FALSE, cont=nb, label = "Downscaling")
tbl4 <- glayout(container=g,horizontal=FALSE,expand=TRUE)
## Layout widgets
for( i in 1:length(texts)) {
tbl4[i,1] <- glabel(texts[i])
font(tbl4[i,1]) <- c(family="times",size="large",style="italic")
}
tbl4[1,2] <- gfilebrowse(text="",width=80,pos=0,type="selectdir",quote=FALSE)
tbl4[2,2] <- gfilebrowse(text="",width=80,pos=0,type="selectdir",quote=FALSE)
tbl4[3,2] <- gfilebrowse(text="",width=80,pos=0,type="selectdir",quote=FALSE)
tbl4[4,2] <- gfilebrowse(text="",width=80,pos=0,type="open",quote=FALSE)
tbl4[5,2] <- gfilebrowse(text="",width=80,pos=0,type="selectdir",quote=FALSE)
addHandlerChanged(tbl4[4,2], function(h,...) 
{
svalue(tbl[2,2]) <- svalue(tbl4[4,2])
svalue(tbl6[2,2]) <- svalue(tbl4[4,2])
})
tbl5 <- glayout(container=g,horizontal=FALSE,expand=TRUE,spacing = 45)
tbl5[1,1] <- glabel("Start Date")
font(tbl5[1,1]) <- c(family="times",size="large",style="italic")
tbl5[1,2] <- gcalendar(format = "%Y-%m-%d") #initial.msg = "YYYY-MM-DD",
tbl5[1,3] <- glabel("End Date")
font(tbl5[1,3]) <- c(family="times",size="large",style="italic")
tbl5[1,4] <-  gcalendar(format = "%Y-%m-%d") #initial.msg = "YYYY-MM-DD", 
#addSpace(g, 5, horizontal=FALSE)

tblradio2 <- glayout(container=g,horizontal=FALSE,expand=TRUE)
tblradio2[2,1] <- glabel("Lapse Rate Reference")
font(tblradio2[2,1]) <- c(family="times",size="large",style="italic")
tblradio2[2,2] <- gradio(c("Calculated","ICIMOD","MicroMet","Logan River","User Defined"),horizontal=TRUE)
font(tblradio2[2,2]) <- c(family="times",size="large",style="italic")
tblradio2[2,3]=gbutton("Update Lapse Rate")
font(tblradio2[2,3]) <- c(family="times",size="large",style="italic")

addHandlerChanged(tblradio2[2,3], function(h,...) {
selctvar=svalue(tblradio2[2,2])
if(selctvar=="Calculated")
{
lapserate=vector("numeric",12)
for( i in 1:6){
lapserate[i] = svalue(tbl71[1,i*2])
lapserate[6+i] = svalue(tbl71[2,i*2])
}
lapse=lapserate
}
if(selctvar=="ICIMOD")
{
lapse=c(0.0092,0.00808,0.00337,0.00391,0.00293,0.00136,0.00347,0.0052,0.00663,0.00885,0.0114,0.00743)
}
if(selctvar=="MicroMet")
{
lapse=c(0.0044,0.0059,0.0071,0.0078,0.0081,0.0082,0.0081,0.0081,0.0077,0.0068,0.0055,0.0047)
}
if(selctvar=="Logan River")
{
lapse=c(0.0044,0.00464,0.0073,0.00604,0.0019,0.00705,0.0015,0.0026,0.00269,0.00606,0.004976,0.00605)
}
if(selctvar=="User Defined")
{
lapse=rep("",12)
}
for( i in 1:6){
font(tbl7[1,(i*2-1)]) <- c(family="times",size="large",style="italic")
svalue(tbl7[1,i*2]) <- as.character(lapse[i],5)
font(tbl7[2,(i*2-1)]) <- c(family="times",size="large",style="italic")
svalue(tbl7[2,i*2]) <- as.character(lapse[6+i],5)
}
})


tblradio3 <- glayout(container=g,horizontal=FALSE,expand=TRUE)
tblradio3[2,1] <- glabel("Precipitation Adjustment Factor Reference")
font(tblradio3[2,1]) <- c(family="times",size="large",style="italic")
tblradio3[2,2] <- gradio(c("Calculated","MicroMet","User Defined"),horizontal=TRUE)
font(tblradio3[2,2]) <- c(family="times",size="large",style="italic")
tblradio3[2,3]=gbutton("Update Precipiation Adjustment Factor")
font(tblradio3[2,3]) <- c(family="times",size="large",style="italic")

addHandlerChanged(tblradio3[2,3], function(h,...) {
selctvar=svalue(tblradio3[2,2])
if(selctvar=="Calculated")
{
lapserateP=vector("numeric",12)
for( i in 1:6){
lapserateP[i] = svalue(tbl70[1,i*2])
lapserateP[6+i] = svalue(tbl70[2,i*2])
}
adjfactor=lapserateP
}
if(selctvar=="MicroMet")
{
adjfactor=c(0.35,0.35,0.35,0.30,0.25,0.20,0.20,0.20,0.20,0.25,0.30,0.35)/1000
}
if(selctvar=="User Defined")
{
adjfactor=rep("",12)
}
for( i in 1:6){
font(tbl8[1,(i*2-1)]) <- c(family="times",size="large",style="italic")
svalue(tbl8[1,i*2]) <- as.character(round(adjfactor[i],5))
}
for( i in 1:6){
font(tbl8[2,(i*2-1)]) <- c(family="times",size="large",style="italic")
svalue(tbl8[2,i*2]) <- as.character(round(adjfactor[6+i],5))
}
})

frame1 <- gframe(container=g,horizontal=TRUE,text="Monthly Lapse Rate (Degree Celsius/meter)")
tbl7 <- glayout(container=frame1,horizontal=TRUE,expand=TRUE)
monthnames=c("January","February","March","April","May","June","July",
             "August","September","October","November","December")
lapse=rep(0,12)
for( i in 1:6){
tbl7[1,(i*2-1)] <- glabel(monthnames[i])
font(tbl7[1,(i*2-1)]) <- c(family="times",size="large",style="italic")
tbl7[1,i*2] <- gedit(text = lapse[i], width = 10, coerce.with = as.numeric)
}
for( i in 1:6){
tbl7[2,(i*2-1)] <- glabel(monthnames[6+i])
font(tbl7[2,(i*2-1)]) <- c(family="times",size="large",style="italic")
tbl7[2,i*2] <- gedit(text = lapse[6+i], width = 10, coerce.with = as.numeric)
}

#addSpace(g, 5, horizontal=FALSE)
adjfactor=rep(0,12)
frame1 <- gframe(container=g,horizontal=TRUE,text="Monthly Precipitation Adjustment Factor")
tbl8 <- glayout(container=frame1,horizontal=TRUE,expand=TRUE)
for( i in 1:6){
tbl8[1,(i*2-1)] <- glabel(monthnames[i])
font(tbl8[1,(i*2-1)]) <- c(family="times",size="large",style="italic")
tbl8[1,i*2] <- gedit(text = adjfactor[i], width = 10, coerce.with = as.numeric)
tbl8[2,(i*2-1)] <- glabel(monthnames[6+i])
font(tbl8[2,(i*2-1)]) <- c(family="times",size="large",style="italic")
tbl8[2,i*2] <- gedit(text = adjfactor[6+i], width = 10, coerce.with = as.numeric)
}
#addSpace(g, 10, horizontal=FALSE)
group2 <- ggroup(container=g,horizontal=TRUE)
addSpring(group2)   
inforButton2 <- gbutton(text="Information","Information", container=group2,handler=function(h,...){
                   confirmDialog(text="Gerneral Info!","This script will downscale MERRA and RFE data 
						and its inputs are:
						1. Directory that holds MERRA and RFE data
						This must have within it folders MERRAcli, MERRArad, RFE 
						holding respectively the climate, radiation and precipitation data
						2. Directory that holds MERRA Geopotential surface representation of elevation files
						3. GeoTIFF DEM file that defines the UEB Model domain
						4. downscaling start and end time
						5. Set precipitation adjustment factor.
						This value is used to multiply RFE precipitation to adjust for RFE precipitation biases")})
tblbtn2 <- gbutton("Downscale", container=group2)
#btn_quit <- gbutton(text = "  Close window  ",container = group2,handler=function(...) dispose(win))
addHandlerChanged(tblbtn2, function(h,...) {

formats="Tiff"
expectedformat="tif"
endstrtemp=unlist(strsplit(svalue(tbl4[4,2]),"[.]"))
entlen=length(endstrtemp)
endstr=endstrtemp[entlen]
if(endstr!=tolower(endstrtemp[entlen]))
{
msg=paste(texts[i+3],": ","Must be a ",formats," file",sep="")
}else{msg=NULL}
stdate=svalue(tbl5[1,2])
endate=svalue(tbl5[1,4])
RFEstdate=as.Date("2001-05-01",format= "%Y-%m-%d")
RFEendate=Sys.Date()-90
chk1=!checkrangedate(RFEstdate, RFEendate, stdate)
chk2=!checkrangedate(RFEstdate, RFEendate, endate)
chk3=!(endate>=stdate)
chks=is.na(c(chk1,chk2,chk3))
msgs2=c("For RFE2, Start date must be between 2001-05-01 and 90 days past from today",
        "For RFE2, End date must be between 2001-05-01 and 90 days past from today",
        "For MERRA, End date must be between 2001-05-01 and 90 days past from today",
        "For MERRA, End date must be between 2001-05-01 and 90 days past from today",
        "End Date muct be higher than Start date")
msgshow2=msgs2[chks]
lapsemsg=vector(length=12)
padjmsg=vector(length=12)
for(j in 1:2)
{
   for(i in 1:6)
   {
if(!is.numeric(svalue(tbl7[j,i*2])))
        {
        lapsemsg[(j-1)*6+i]=paste("Lapse rate for ",monthnames[(j-1)*6+i], " must be number")
        }else{
        lapsemsg[(j-1)*6+i]=NA
        }
if(!is.numeric(svalue(tbl8[j,i*2])))
        {padjmsg[(j-1)*6+i]=paste("Precipitation bias corrector for ", monthnames[(j-1)*6+i], "must be a number")
        }else{
        padjmsg[(j-1)*6+i]=NA
        }
   }
}
lapsemsg <- lapsemsg[!is.na(lapsemsg)]
padjmsg <- padjmsg[!is.na(padjmsg)]

if(length(c(msg,msgshow2,padjmsg,lapsemsg))>0)
{
confirmDialog(text="Error!",c(msg,padjmsg,msgshow2))
}else{
VarSourceDir=gsub("\\\\", "/",svalue(tbl4[1,2]))
GeoSourceDir=gsub("\\\\", "/",svalue(tbl4[2,2]))
DestDir=gsub("\\\\", "/",svalue(tbl4[3,2]))
DEMfile=gsub("\\\\", "/",svalue(tbl4[4,2]))
latfolder=gsub("\\\\", "/",svalue(tbl4[5,2]))
startdate=gsub("\\\\", "/",svalue(tbl5[1,2]))
enddate=gsub("\\\\", "/",svalue(tbl5[1,4]))
latlonfile = paste(latfolder,"/","lat-lon.txt",sep="")
PrecipDataSource=svalue(tblradio[1,2])
lapse=c(as.numeric(svalue(tbl7[1,2])),as.numeric(svalue(tbl7[1,4])),as.numeric(svalue(tbl7[1,6])),
            as.numeric(svalue(tbl7[1,8])),as.numeric(svalue(tbl7[1,10])),as.numeric(svalue(tbl7[1,12])),
            as.numeric(svalue(tbl7[2,2])),as.numeric(svalue(tbl7[2,4])),as.numeric(svalue(tbl7[2,6])),
            as.numeric(svalue(tbl7[2,8])),as.numeric(svalue(tbl7[2,10])),as.numeric(svalue(tbl7[2,12])))
adjfactor=c(as.numeric(svalue(tbl8[1,2])),as.numeric(svalue(tbl8[1,4])),as.numeric(svalue(tbl8[1,6])),
            as.numeric(svalue(tbl8[1,8])),as.numeric(svalue(tbl8[1,10])),as.numeric(svalue(tbl8[1,12])),
            as.numeric(svalue(tbl8[2,2])),as.numeric(svalue(tbl8[2,4])),as.numeric(svalue(tbl8[2,6])),
            as.numeric(svalue(tbl8[2,8])),as.numeric(svalue(tbl8[2,10])),as.numeric(svalue(tbl8[2,12])))
arglist=list(VarSourceDir,GeoSourceDir,DestDir,DEMfile,startdate,enddate,lapse,adjfactor,latlonfile,PrecipDataSource)
do.call("DownscaleMERRARFE",arglist)
confirmDialog(text="Success!","data dowsnscaling is successfully completed!",img="ok")
}
})


#WINDOW 4------------------------------------------
#--------------------------------------------------
#--------------------------------------------------
#texts to help the users
texts=c("Land Cover Tiff File","DEM Tiff File",
        "Land Cover Classes CSV file","Output Folder")
g <- ggroup(horizontal=FALSE, cont=nb, label = "Land Cover Varibles")
tbl <- glayout(container=g,horizontal=FALSE,expand=TRUE)
## Layout widgets
for( i in 1:length(texts)) {
tbl[i,1] <- glabel(texts[i])
font(tbl[i,1]) <- c(family="times",size="large",style="italic")
}
browse1=paste(getwd(),"\\MODIS","\\MOD12Q1_LandCoverType_IGBP_SouthAsia_Geographic.tif",sep="")
browse2=paste(getwd(),"\\MODIS","\\landcover.csv",sep="")
tbl[1,2] <- gfilebrowse(browse1,pos=0,type="open",quote=FALSE)
tbl[2,2] <- gfilebrowse(pos=0,type="open",quote=FALSE)
tbl[3,2] <- gfilebrowse(browse2,width=70,pos=0,type="open",quote=FALSE)
tbl[4,2] <- gfilebrowse(text=getwd(),width=80,pos=0,type="selectdir",quote=FALSE)
#addSpace(g, 10, horizontal=FALSE)
group2 <- ggroup(container=g,horizontal=TRUE)
addSpring(group2)   
inforButton3 <- gbutton("Information", container=group2,handler=function(h,...){
                   confirmDialog(text="General Info!","This script will 4 netCDF files for vegetation site variable from MODIS land cover data
						and its inputs are:
						1. Geotiff file that holds MODIS data
						2. CSV file that holds lookup table mapping MODIS classes onto vegetation parameters
						3. GeoTIFF DEM file that defines the Model domain")})
tblbtn <- gbutton("Create Land Cover variables", container=group2)
#btn_quit <- gbutton(text = "  Close window  ",container = group2, handler=function(...) dispose(win))
addHandlerChanged(tblbtn, function(h,...) {
msgindex=c(1:3)
formats=c("Tiff","Tiff","csv")
expectedformat=c("tif","tif","csv")
endstr=vector("character",length=length(texts))
msg=vector("character",length=length(texts))
for(i in msgindex)
{
endstrtemp=unlist(strsplit(svalue(tbl[i,2]),"[.]"))
entlen=length(endstrtemp)
if(entlen > 0 & !is.na(formats[i])){
endstr[i]=tolower(endstrtemp[entlen])}
msg[i]=paste(texts[i],": ","Must be a ",formats[i]," file",sep="")
}
massages=msg[endstr != expectedformat]
if(length(massages)>0)
{
confirmDialog(text="Error!",msg[endstr != expectedformat])
}else{
do.call("vegfiles",list(svalue(tbl[1,2]),svalue(tbl[2,2]),
        svalue(tbl[3,2]),svalue(tbl[4,2])))
confirmDialog(text="Success!","Vegegation related variable files are successfully created!",img="ok")
}
})

#WINDOW 5------------------------------------------
#--------------------------------------------------
#--------------------------------------------------
#texts to help the users
texts=c("Temperature Range", "DEM Tiff File", "Output Folder")
g <- ggroup(horizontal=FALSE, cont=nb, label = "Temperature Range Varibles")
tbl6 <- glayout(container=g,horizontal=FALSE,expand=TRUE)
## Layout widgets
for( i in 1:length(texts)) {
tbl6[i,1] <- glabel(texts[i])
font(tbl6[i,1]) <- c(family="times",size="large",style="italic")
}
for( i in 1:(length(texts)-1)) {
tbl6[i,2] <- gfilebrowse(width=70,pos=0,type="open",quote=FALSE)
}
tbl6[3,2]  <- gfilebrowse(text=getwd(),width=80,pos=0,type="selectdir",quote=FALSE)
group2 <- ggroup(container=g,horizontal=TRUE)
addSpring(group2) 
inforButton4 <- gbutton("Information", container=group2,handler=function(h,...){
                   confirmDialog(text="General Info!","This script will write 12 netCDF files containing monthly diurnal temperature range 
						and its inputs are:
						1. NetCDF file with climatological diurnal temperature ranges
						2. GeoTIFF DEM file that defines the Model domain")})
tblbtn4 <- gbutton("Create climatology files", container=group2)
#btn_quit <- gbutton(text = "  Close window  ",container = group2, handler=function(...) dispose(win))
addHandlerChanged(tblbtn4, function(h,...) {
msgindex=c(1:3)
formats=c("netCDF","Tiff",NA)
expectedformat=c("nc","tif","")
endstr=vector("character",length=length(texts))
msg=vector("character",length=length(texts))
for(i in msgindex)
{
endstrtemp=unlist(strsplit(svalue(tbl6[i,2]),"[.]"))
entlen=length(endstrtemp)
if(entlen > 0 & !is.na(formats[i]) ){
endstr[i]=tolower(endstrtemp[entlen])}
msg[i]=paste(texts[i],": ","Must be a ",formats[i]," file",sep="")
}
massages=msg[endstr != expectedformat]
if(length(massages)>0)
{
confirmDialog(text="Error!",msg[endstr != expectedformat])
}else{
do.call("MERRAAvgMonthlyTrange",list(svalue(tbl6[1,2]),svalue(tbl6[2,2]),svalue(tbl6[3,2])))
confirmDialog(text="Success!","Climatological variables are successfully created!",img="ok")
}
})

# keep the windows open when opeing from command prompt
while(isExtant(win)) Sys.sleep(1)