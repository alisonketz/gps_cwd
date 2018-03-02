###
### 3/1/2018 
### Alison Ketz
### Movement data preprocessing using adehabitatLT package
###


###
### Preliminaries
###

rm(list=ls())

library(e1071)
library(geosphere)
library(lubridate)
library(Hmisc)
library(zoo)
library(stringr)
library(ggplot2)
library(dplyr)
library(ggmap)
library(adehabitatHR)
library(adehabitatLT)
library(maptools)
library(changepoint)
library(sp)
library(spatstat)#for "duplicate" function
library(readr)
library(mvtnorm)

###
### set working directory
###

setwd("/home/aketz/Documents/GPS_cwd/180301_gps_cwd/gps_cwd")

###
### For running on windows machine
###
# 
# library(RODBC)
# setwd("F:/GPSevent/180301_GPSevent")

# myconn=odbcConnectAccess('C:/Users/aketz/Documents/Data/SWDPPdeerDB.MDB')
# d.vit <- sqlFetch(myconn, "VIT")
# d.ramalt <-sqlFetch(myconn, "RAMALT")
# close(myconn)


###
### Load ramalt/capture/vit data
###

d.adultcap=mdb.get('~/Documents/Data/SWDPPdeerDB.MDB',tables= "Adult_Capture_2017_2018")
d.ramalt=mdb.get('~/Documents/Data/SWDPPdeerDB.MDB',tables= "RAMALT")
d.vit=mdb.get('~/Documents/Data/SWDPPdeerDB.MDB',tables= "VIT")
names(d.vit)=tolower(gsub('[[:punct:]]',"",names(d.vit)))
names(d.vit)=gsub(" ", "",names(d.vit), fixed = TRUE)

names(d.ramalt)=tolower(gsub('[[:punct:]]',"",names(d.ramalt)))
names(d.ramalt)=gsub(" ", "",names(d.ramalt), fixed = TRUE)
d.ramalt=d.ramalt[order(d.ramalt$lowtag),]
head(d.ramalt)
cwd.pos=d.ramalt$lowtag[d.ramalt$result=="Positive"]
# write.csv(cbind(cwd.pos),file="cwd_positive_list.csv")

cwd.male=d.ramalt$lowtag[d.ramalt$result=="Positive"][c(1,3,5,7,11)]
names(d.adultcap)=tolower(gsub('[[:punct:]]',"",names(d.adultcap)))
names(d.adultcap)=gsub(" ", "",names(d.adultcap), fixed = TRUE)
malecollarnum=c()
for(j in cwd.male){
    malecollarnum=c(malecollarnum,d.adultcap$collarid[d.adultcap$lowtag==j])
}
cwd.female = setdiff(cwd.pos,cwd.male)

###
### Read and clean gps data
###

###
### does
###
d = matrix(NA,nr=1,nc=13)
#for loop for reading in data, using vector of lowtag's from the vit dataset
for(j in 1:length(cwd.female)){
    d.temp = read.table(paste("/home/aketz/Documents/Data/GPS_locations_CWD/cwd_does/",cwd.female[j],".csv",sep=""),sep=",",header=TRUE,stringsAsFactors = FALSE)
    d.temp$lowtag =cwd.female[j]
    names(d.temp)=tolower(names(d.temp))
    d=rbind(d,as.matrix(d.temp))
}
d=d[-1,]
d=data.frame(d,stringsAsFactors = FALSE)

for(j in 1:dim(d)[2]){
    d[,j]=str_trim(d[,j],side="both")
}

class.indx=c(5:7,9:12)
for(j in class.indx){
    d[,j]=as.numeric(d[,j])
}

d$lowtag=as.factor(d$lowtag)
d=d[,-1]
d$sex="female"

###
### bucks
###
# try=read.table("~/Documents/Data/GPS_locations_CWD/cwd_bucks/5050.csv",sep=";",header=TRUE)
# head(try)
# dim(try)
# names(try)
# [9]
# names(d)
# length(names(try)[c(2,4:7,9:10,14:18,45:49)])

d.bucks.temp = matrix(NA,nr=1,nc=18)
#for loop for reading in data, using vector of lowtag's from the vit dataset
for(j in 1:length(cwd.male)){
    d.temp = read.table(paste("/home/aketz/Documents/Data/GPS_locations_CWD/cwd_bucks/",cwd.male[j],".csv",sep=""),sep=";",header=TRUE,stringsAsFactors = FALSE)
    d.temp = d.temp[,c(2,4:7,9:10,14:18,45:49)]
    d.temp$lowtag =cwd.male[j]
    names(d.temp)=tolower(names(d.temp))
    d.bucks.temp=rbind(d.bucks.temp,as.matrix(d.temp))
}
d.bucks.temp=d.bucks.temp[-1,]
d.bucks.temp=data.frame(d.bucks.temp,stringsAsFactors = FALSE)
for(j in 1:dim(d.bucks.temp)[2]){
    d.bucks.temp[,j]=str_trim(d.bucks.temp[,j],side="both")
}
d.bucks.temp$sex="male"
d.bucks=data.frame(matrix(NA,nc=dim(d)[2],nr=dim(d.bucks.temp)[1]))
names(d.bucks)=names(d)
d.bucks[,c(1,4:13)]=d.bucks.temp[,c(1,8:10,12,11,17,15,16,18,19)]

d.bucks.temp[,2]=paste(substr(d.bucks.temp[,2],6,8),substr(d.bucks.temp[,2],9,10),substr(d.bucks.temp[,2],5,5),substr(d.bucks.temp[,2],1,4),sep="")
d.bucks.temp[,4]=paste(substr(d.bucks.temp[,4],6,8),substr(d.bucks.temp[,4],9,10),substr(d.bucks.temp[,4],5,5),substr(d.bucks.temp[,4],1,4),sep="")

d.bucks[,2]=paste(d.bucks.temp[,2],d.bucks.temp[,3])
d.bucks[,3]=paste(d.bucks.temp[,4],d.bucks.temp[,5])
head(d.bucks)
head(d)
d[,2]=paste(substr(d[,2],1,11),substr(d[,2],13,nchar(d[,2])),sep="")
d[,3]=paste(substr(d[,3],1,11),substr(d[,3],13,nchar(d[,3])),sep="")


class.indx=c(4:6,8:11)
for(j in class.indx){
    d.bucks[,j]=as.numeric(d.bucks[,j])
}
class(d.bucks[,c(4)])


d=rbind(d,d.bucks)

head(d[,2])
tail(d[,2])

###
### Adding metrics
###


###
### Converting date to POSIXct format
###

d$date_time_local=as.POSIXct(strptime(d$date_time_local,format="%m-%d-%Y %H:%M:%S"),tz="CST6CDT")
d$date_time_gmt=as.POSIXct(strptime(d$date_time_gmt,format="%m-%d-%Y %H:%M:%S"),tz="GMT")



#Create time lag between successive locations to censor data if needed.
time.diff <- diff(d$date_time_local)
d=d[-1,]
d$timediff <-round(as.numeric(abs(time.diff)))
rm=which(d$timediff>10)
d=d[-rm,]
names(d)[1]="devname"
