#Community stats
library(vegan)
#Normality testing
library(nortest)
library(vcd)
#Spatial autocorrelation
library(ade4)
library(ape)
#Data cleaning and plotting
library(stringr)
library(tidyverse)
library(lattice)
library(ggpubr)
library(magick)
#Time series processing
library(imputeTS)
library(signal)
#Modeling
library(DHARMa)
library(robustHD)
library(gamm4)
library(PerformanceAnalytics)
library(MuMIn)
library(lme4)
library(car)
library(mvabund)
#Raster calculations
library(raster)
library(rgdal)
library(rgeos)
library(spatialEco)
library(exactextractr)

#Set working directory
setwd("D:/")

########### Calculating response variables ###########
#Load community matrix
hopper<-read.csv("Excel_sheets/Caelifera_assemblage_data.csv")
#Calculate response variables
Richness<-rowSums(hopper!=0)
exShan<-exp(diversity(hopper,index="shannon"))
#Add response variables to data frame
GPS<-read.csv("Excel_sheets/Site_GPS_coordinates.csv")
df<-as.data.frame(cbind(GPS,"Richness"=Richness,"exShannon"=exShan))
#Clean environment
rm(Richness,exShan)

########### Accumulation curves ###########
#Plot overall curve
Curve<-specaccum(hopper,method="rarefaction")
plot(Curve,main="Species rarefaction",ylab="Rarefaction")
#Plot per plantation
Plantation<-GPS$Plantation
hopperPlant<-as.data.frame(cbind(hopper,"Plantation"=Plantation))
#Subset
Estate1<-subset(hopperPlant,Plantation=="Estate1")
Estate1$Plantation=NULL
Estate2<-subset(hopperPlant,Plantation=="Estate2")
Estate2$Plantation=NULL
Estate3<-subset(hopperPlant,Plantation=="Estate3")
Estate3$Plantation=NULL
Estate4<-subset(hopperPlant,Plantation=="Estate4")
Estate4$Plantation=NULL
#Plotting
Estate1Curve<-specaccum(Estate1,method="rarefaction")
Estate2Curve<-specaccum(Estate2,method="rarefaction")
Estate3Curve<-specaccum(Estate3,method="rarefaction")
Estate4Curve<-specaccum(Estate4,method="rarefaction")
plot(Estate1Curve,main="Species rarefaction",ylab="Rarefaction")
plot(Estate2Curve,add=TRUE,col='red')
plot(Estate3Curve,add=TRUE,col='blue')
plot(Estate4Curve,add=TRUE,col='green')
legend(x="bottomright",legend=c("Estate 1","Estate 2","Estate 3","Estate 4"),fill=c("black","red","blue","green"))
#Clean environment
rm(Curve,Plantation,hopperPlant,Estate1,Estate2,Estate3,Estate4,Estate1Curve,Estate2Curve,Estate3Curve,Estate4Curve)

########### Normality testing ###########
#Richness
shapiro.test(df$Richness) #Normal
hist(df$Richness)
#exShannon       
shapiro.test(df$exShan) #Normal
hist(df$exShan)

########### Spatial autocorrelation ###########
#Monte-Carlo
provi<-deldir::deldir(GPS[2:3])
provi.neig<-neig(edges=as.matrix(provi$delsgs[,5:6]))
maf.listw<-spdep::nb2listw(neig2nb(provi.neig))
multispati.rtest((dudi.pca(hopper,scannf=FALSE)),maf.listw) #Spatial autocorrelation
#Moran's I  
hopper.dist.inv<-1/(as.matrix(dist(cbind(GPS$Long_X,GPS$Lat_Y))))
hopper.dist.inv[is.infinite(hopper.dist.inv)] <- 0
Moran.I(df$Richness,hopper.dist.inv) #Spatial autocorrelation
Moran.I(df$exShannon,hopper.dist.inv) #Spatial autocorrelation
#Clean environment
rm(provi,provi.neig,maf.listw,hopper.dist.inv)

########### Earth Engine data processing ###########
#Clean Landsat NDVI time series data
Index<- read.csv("Excel_sheets/Landsat_SR_NDVI.csv")
rownames(Index)<-Index$Site_ID
#Remove unwanted Columns
Index$system.index=NULL
Index$.geo=NULL
Index$Site_ID=NULL
#Rename column headers and sort matrix
colnames(Index) <- str_sub((colnames(Index)),-8,-1)
Index<-Index[order(rownames(Index)),]
Index<-Index[ , order(names(Index))]
#Remove multiples values by calculating their means
Index<-sapply(split.default(Index, sub('\\..*', '', names(Index))), rowMeans, na.rm = TRUE)
#Select imagery from 1995 - 2020
Index<-Index[,93:1010]
#Remove outliers
sum(is.na(Index))
x <- abs(Index - median(Index,na.rm = T)) / mad(Index,na.rm = T)
tr <- ifelse(quantile(x,.99,na.rm=T)>2,quantile(x,.99,na.rm=T),NA)
if(!is.na(tr)){
  Index[ which(x > tr )] <- NA
}
sum(is.na(Index))
#Clean environment
rm(x,tr)
#Impute missing data using Kalman smoother and apply sgolayfilt Savitzky-Golay filtering
IndexT<-t(Index)
ColNames<-colnames(IndexT)
RowNames<-rownames(IndexT)
Site_00<-sgolayfilt(na_kalman(IndexT[,1],model = "auto.arima"))
Site_01<-sgolayfilt(na_kalman(IndexT[,2],model = "auto.arima"))
Site_02<-sgolayfilt(na_kalman(IndexT[,3],model = "auto.arima"))
Site_03<-sgolayfilt(na_kalman(IndexT[,4],model = "auto.arima"))
Site_04<-sgolayfilt(na_kalman(IndexT[,5],model = "auto.arima"))
Site_05<-sgolayfilt(na_kalman(IndexT[,6],model = "auto.arima"))
Site_06<-sgolayfilt(na_kalman(IndexT[,7],model = "auto.arima"))
Site_07<-sgolayfilt(na_kalman(IndexT[,8],model = "auto.arima"))
Site_08<-sgolayfilt(na_kalman(IndexT[,9],model = "auto.arima"))
Site_09<-sgolayfilt(na_kalman(IndexT[,10],model = "auto.arima"))
Site_10<-sgolayfilt(na_kalman(IndexT[,11],model = "auto.arima"))
Site_11<-sgolayfilt(na_kalman(IndexT[,12],model = "auto.arima"))
Site_12<-sgolayfilt(na_kalman(IndexT[,13],model = "auto.arima"))
Site_13<-sgolayfilt(na_kalman(IndexT[,14],model = "auto.arima"))
Site_14<-sgolayfilt(na_kalman(IndexT[,15],model = "auto.arima"))
Site_15<-sgolayfilt(na_kalman(IndexT[,16],model = "auto.arima"))
Site_16<-sgolayfilt(na_kalman(IndexT[,17],model = "auto.arima"))
Site_17<-sgolayfilt(na_kalman(IndexT[,18],model = "auto.arima"))
Site_18<-sgolayfilt(na_kalman(IndexT[,19],model = "auto.arima"))
Site_19<-sgolayfilt(na_kalman(IndexT[,20],model = "auto.arima"))
Site_20<-sgolayfilt(na_kalman(IndexT[,21],model = "auto.arima"))
Site_21<-sgolayfilt(na_kalman(IndexT[,22],model = "auto.arima"))
Site_22<-sgolayfilt(na_kalman(IndexT[,23],model = "auto.arima"))
Site_23<-sgolayfilt(na_kalman(IndexT[,24],model = "auto.arima"))
Site_24<-sgolayfilt(na_kalman(IndexT[,25],model = "auto.arima"))
Site_25<-sgolayfilt(na_kalman(IndexT[,26],model = "auto.arima"))
Site_26<-sgolayfilt(na_kalman(IndexT[,27],model = "auto.arima"))
Site_27<-sgolayfilt(na_kalman(IndexT[,28],model = "auto.arima"))
Site_28<-sgolayfilt(na_kalman(IndexT[,29],model = "auto.arima"))
Site_29<-sgolayfilt(na_kalman(IndexT[,30],model = "auto.arima"))
Site_30<-sgolayfilt(na_kalman(IndexT[,31],model = "auto.arima"))
Site_31<-sgolayfilt(na_kalman(IndexT[,32],model = "auto.arima"))
Site_32<-sgolayfilt(na_kalman(IndexT[,33],model = "auto.arima"))
Site_33<-sgolayfilt(na_kalman(IndexT[,34],model = "auto.arima"))
Site_34<-sgolayfilt(na_kalman(IndexT[,35],model = "auto.arima"))
Site_35<-sgolayfilt(na_kalman(IndexT[,36],model = "auto.arima"))
Site_36<-sgolayfilt(na_kalman(IndexT[,37],model = "auto.arima"))
Site_37<-sgolayfilt(na_kalman(IndexT[,38],model = "auto.arima"))
Site_38<-sgolayfilt(na_kalman(IndexT[,39],model = "auto.arima"))
Site_39<-sgolayfilt(na_kalman(IndexT[,40],model = "auto.arima"))
Site_40<-sgolayfilt(na_kalman(IndexT[,41],model = "auto.arima"))
Site_41<-sgolayfilt(na_kalman(IndexT[,42],model = "auto.arima"))
Site_42<-sgolayfilt(na_kalman(IndexT[,43],model = "auto.arima"))
Site_43<-sgolayfilt(na_kalman(IndexT[,44],model = "auto.arima"))
Site_44<-sgolayfilt(na_kalman(IndexT[,45],model = "auto.arima"))
Site_45<-sgolayfilt(na_kalman(IndexT[,46],model = "auto.arima"))
Site_46<-sgolayfilt(na_kalman(IndexT[,47],model = "auto.arima"))
Site_47<-sgolayfilt(na_kalman(IndexT[,48],model = "auto.arima"))
Site_48<-sgolayfilt(na_kalman(IndexT[,49],model = "auto.arima"))
Site_49<-sgolayfilt(na_kalman(IndexT[,50],model = "auto.arima"))
Site_50<-sgolayfilt(na_kalman(IndexT[,51],model = "auto.arima"))
#Add to matrix
KalmarDataNDVI<-cbind(Site_00,Site_01,Site_02,Site_03,Site_04,Site_05,Site_06,Site_07,Site_08,
                      Site_09,Site_10,Site_11,Site_12,Site_13,Site_14,Site_15,Site_16,Site_17,
                      Site_18,Site_19,Site_20,Site_21,Site_22,Site_23,Site_24,Site_25,Site_26,
                      Site_27,Site_28,Site_29,Site_30,Site_31,Site_32,Site_33,Site_34,Site_35,
                      Site_36,Site_37,Site_38,Site_39,Site_40,Site_41,Site_42,Site_43,Site_44,
                      Site_45,Site_46,Site_47,Site_48,Site_49,Site_50)
colnames(KalmarDataNDVI)<-ColNames
rownames(KalmarDataNDVI)<-RowNames
KalmarDataNDVI<-as.data.frame(t(KalmarDataNDVI))
#Clean environment
rm(Site_00,Site_01,Site_02,Site_03,Site_04,Site_05,Site_06,Site_07,Site_08,
   Site_09,Site_10,Site_11,Site_12,Site_13,Site_14,Site_15,Site_16,Site_17,
   Site_18,Site_19,Site_20,Site_21,Site_22,Site_23,Site_24,Site_25,Site_26,
   Site_27,Site_28,Site_29,Site_30,Site_31,Site_32,Site_33,Site_34,Site_35,
   Site_36,Site_37,Site_38,Site_39,Site_40,Site_41,Site_42,Site_43,Site_44,
   Site_45,Site_46,Site_47,Site_48,Site_49,Site_50,ColNames,RowNames,Index,
   IndexT)
#Select cumalitive time stamps
Y19_20<-KalmarDataNDVI %>% dplyr::select(starts_with(c("201904","201905","201906","201907","201908","201909","201910","201911","201912","202001","202002","202003")))
Y18_19<-KalmarDataNDVI %>% dplyr::select(starts_with(c("201804","201805","201806","201807","201808","201809","201810","201811","201812","201901","201902","201903")))
Y17_18<-KalmarDataNDVI %>% dplyr::select(starts_with(c("201704","201705","201706","201707","201708","201709","201710","201711","201712","201801","201802","201803")))
Y16_17<-KalmarDataNDVI %>% dplyr::select(starts_with(c("201604","201605","201606","201607","201608","201609","201610","201611","201612","201701","201702","201703")))
Y15_16<-KalmarDataNDVI %>% dplyr::select(starts_with(c("201504","201505","201506","201507","201508","201509","201510","201511","201512","201601","201602","201603")))
Y14_15<-KalmarDataNDVI %>% dplyr::select(starts_with(c("201404","201405","201406","201407","201408","201409","201410","201411","201412","201501","201502","201503")))
Y13_14<-KalmarDataNDVI %>% dplyr::select(starts_with(c("201304","201305","201306","201307","201308","201309","201310","201311","201312","201401","201402","201403")))
Y12_13<-KalmarDataNDVI %>% dplyr::select(starts_with(c("201204","201205","201206","201207","201208","201209","201210","201211","201212","201301","201302","201303")))
Y11_12<-KalmarDataNDVI %>% dplyr::select(starts_with(c("201104","201105","201106","201107","201108","201109","201110","201111","201112","201201","201202","201203")))
Y10_11<-KalmarDataNDVI %>% dplyr::select(starts_with(c("201004","201005","201006","201007","201008","201009","201010","201011","201012","201101","201102","201103")))
Y09_10<-KalmarDataNDVI %>% dplyr::select(starts_with(c("200904","200905","200906","200907","200908","200909","200910","200911","200912","201001","201002","201003")))
Y08_09<-KalmarDataNDVI %>% dplyr::select(starts_with(c("200804","200805","200806","200807","200808","200809","200810","200811","200812","200901","200902","200903")))
Y07_08<-KalmarDataNDVI %>% dplyr::select(starts_with(c("200704","200705","200706","200707","200708","200709","200710","200711","200712","200801","200802","200803")))
Y06_07<-KalmarDataNDVI %>% dplyr::select(starts_with(c("200604","200605","200606","200607","200608","200609","200610","200611","200612","200701","200702","200703")))
Y05_06<-KalmarDataNDVI %>% dplyr::select(starts_with(c("200504","200505","200506","200507","200508","200509","200510","200511","200512","200601","200602","200603")))
Y04_05<-KalmarDataNDVI %>% dplyr::select(starts_with(c("200404","200405","200406","200407","200408","200409","200410","200411","200412","200501","200502","200503")))
Y03_04<-KalmarDataNDVI %>% dplyr::select(starts_with(c("200304","200305","200306","200307","200308","200309","200310","200311","200312","200401","200402","200403")))
Y02_03<-KalmarDataNDVI %>% dplyr::select(starts_with(c("200204","200205","200206","200207","200208","200209","200210","200211","200212","200301","200302","200303")))
Y01_02<-KalmarDataNDVI %>% dplyr::select(starts_with(c("200104","200105","200106","200107","200108","200109","200110","200111","200112","200201","200202","200203")))
Y00_01<-KalmarDataNDVI %>% dplyr::select(starts_with(c("200004","200005","200006","200007","200008","200009","200010","200011","200012","200101","200102","200103")))
Y99_00<-KalmarDataNDVI %>% dplyr::select(starts_with(c("199904","199905","199906","199907","199908","199909","199910","199911","199912","200001","200002","200003")))
Y98_99<-KalmarDataNDVI %>% dplyr::select(starts_with(c("199804","199805","199806","199807","199808","199809","199810","199811","199812","199901","199902","199903")))
Y97_98<-KalmarDataNDVI %>% dplyr::select(starts_with(c("199704","199705","199706","199707","199708","199709","199710","199711","199712","199801","199802","199803")))
Y96_97<-KalmarDataNDVI %>% dplyr::select(starts_with(c("199604","199605","199606","199607","199608","199609","199610","199611","199612","199701","199702","199703")))
Y95_96<-KalmarDataNDVI %>% dplyr::select(starts_with(c("199504","199505","199506","199507","199508","199509","199510","199511","199512","199601","199602","199603")))
#Calculate SD and add to data frame
sdY19_20<-apply(Y19_20, MARGIN=1, FUN=sd, na.rm=TRUE)
sdY18_20<-apply(cbind(Y18_19,Y19_20), MARGIN=1, FUN=sd, na.rm=TRUE)
sdY17_20<-apply(cbind(Y17_18,Y18_19,Y19_20), MARGIN=1, FUN=sd, na.rm=TRUE)
sdY16_20<-apply(cbind(Y16_17,Y17_18,Y18_19,Y19_20), MARGIN=1, FUN=sd, na.rm=TRUE)
sdY15_20<-apply(cbind(Y15_16,Y16_17,Y17_18,Y18_19,Y19_20), MARGIN=1, FUN=sd, na.rm=TRUE)
sdY14_20<-apply(cbind(Y14_15,Y15_16,Y16_17,Y17_18,Y18_19,Y19_20), MARGIN=1, FUN=sd, na.rm=TRUE)
sdY13_20<-apply(cbind(Y13_14,Y14_15,Y15_16,Y16_17,Y17_18,Y18_19,Y19_20), MARGIN=1, FUN=sd, na.rm=TRUE)
sdY12_20<-apply(cbind(Y12_13,Y13_14,Y14_15,Y15_16,Y16_17,Y17_18,Y18_19,Y19_20), MARGIN=1, FUN=sd, na.rm=TRUE)
sdY11_20<-apply(cbind(Y11_12,Y12_13,Y13_14,Y14_15,Y15_16,Y16_17,Y17_18,Y18_19,Y19_20), MARGIN=1, FUN=sd, na.rm=TRUE)
sdY10_20<-apply(cbind(Y10_11,Y11_12,Y12_13,Y13_14,Y14_15,Y15_16,Y16_17,Y17_18,Y18_19,Y19_20), MARGIN=1, FUN=sd, na.rm=TRUE)
sdY09_20<-apply(cbind(Y09_10,Y10_11,Y11_12,Y12_13,Y13_14,Y14_15,Y15_16,Y16_17,Y17_18,Y18_19,Y19_20), MARGIN=1, FUN=sd, na.rm=TRUE)
sdY08_20<-apply(cbind(Y08_09,Y09_10,Y10_11,Y11_12,Y12_13,Y13_14,Y14_15,Y15_16,Y16_17,Y17_18,Y18_19,Y19_20), MARGIN=1, FUN=sd, na.rm=TRUE)
sdY07_20<-apply(cbind(Y07_08,Y08_09,Y09_10,Y10_11,Y11_12,Y12_13,Y13_14,Y14_15,Y15_16,Y16_17,Y17_18,Y18_19,Y19_20), MARGIN=1, FUN=sd, na.rm=TRUE)
sdY06_20<-apply(cbind(Y06_07,Y07_08,Y08_09,Y09_10,Y10_11,Y11_12,Y12_13,Y13_14,Y14_15,Y15_16,Y16_17,Y17_18,Y18_19,Y19_20), MARGIN=1, FUN=sd, na.rm=TRUE)
sdY05_20<-apply(cbind(Y05_06,Y06_07,Y07_08,Y08_09,Y09_10,Y10_11,Y11_12,Y12_13,Y13_14,Y14_15,Y15_16,Y16_17,Y17_18,Y18_19,Y19_20), MARGIN=1, FUN=sd, na.rm=TRUE)
sdY04_20<-apply(cbind(Y04_05,Y05_06,Y06_07,Y07_08,Y08_09,Y09_10,Y10_11,Y11_12,Y12_13,Y13_14,Y14_15,Y15_16,Y16_17,Y17_18,Y18_19,Y19_20), MARGIN=1, FUN=sd, na.rm=TRUE)
sdY03_20<-apply(cbind(Y03_04,Y04_05,Y05_06,Y06_07,Y07_08,Y08_09,Y09_10,Y10_11,Y11_12,Y12_13,Y13_14,Y14_15,Y15_16,Y16_17,Y17_18,Y18_19,Y19_20), MARGIN=1, FUN=sd, na.rm=TRUE)
sdY02_20<-apply(cbind(Y02_03,Y03_04,Y04_05,Y05_06,Y06_07,Y07_08,Y08_09,Y09_10,Y10_11,Y11_12,Y12_13,Y13_14,Y14_15,Y15_16,Y16_17,Y17_18,Y18_19,Y19_20), MARGIN=1, FUN=sd, na.rm=TRUE)
sdY01_20<-apply(cbind(Y01_02,Y02_03,Y03_04,Y04_05,Y05_06,Y06_07,Y07_08,Y08_09,Y09_10,Y10_11,Y11_12,Y12_13,Y13_14,Y14_15,Y15_16,Y16_17,Y17_18,Y18_19,Y19_20), MARGIN=1, FUN=sd, na.rm=TRUE)
sdY00_20<-apply(cbind(Y00_01,Y01_02,Y02_03,Y03_04,Y04_05,Y05_06,Y06_07,Y07_08,Y08_09,Y09_10,Y10_11,Y11_12,Y12_13,Y13_14,Y14_15,Y15_16,Y16_17,Y17_18,Y18_19,Y19_20), MARGIN=1, FUN=sd, na.rm=TRUE)
sdY99_20<-apply(cbind(Y99_00,Y00_01,Y01_02,Y02_03,Y03_04,Y04_05,Y05_06,Y06_07,Y07_08,Y08_09,Y09_10,Y10_11,Y11_12,Y12_13,Y13_14,Y14_15,Y15_16,Y16_17,Y17_18,Y18_19,Y19_20), MARGIN=1, FUN=sd, na.rm=TRUE)
sdY98_20<-apply(cbind(Y98_99,Y99_00,Y00_01,Y01_02,Y02_03,Y03_04,Y04_05,Y05_06,Y06_07,Y07_08,Y08_09,Y09_10,Y10_11,Y11_12,Y12_13,Y13_14,Y14_15,Y15_16,Y16_17,Y17_18,Y18_19,Y19_20), MARGIN=1, FUN=sd, na.rm=TRUE)
sdY97_20<-apply(cbind(Y97_98,Y98_99,Y99_00,Y00_01,Y01_02,Y02_03,Y03_04,Y04_05,Y05_06,Y06_07,Y07_08,Y08_09,Y09_10,Y10_11,Y11_12,Y12_13,Y13_14,Y14_15,Y15_16,Y16_17,Y17_18,Y18_19,Y19_20), MARGIN=1, FUN=sd, na.rm=TRUE)
sdY96_20<-apply(cbind(Y96_97,Y97_98,Y98_99,Y99_00,Y00_01,Y01_02,Y02_03,Y03_04,Y04_05,Y05_06,Y06_07,Y07_08,Y08_09,Y09_10,Y10_11,Y11_12,Y12_13,Y13_14,Y14_15,Y15_16,Y16_17,Y17_18,Y18_19,Y19_20), MARGIN=1, FUN=sd, na.rm=TRUE)
sdY95_20<-apply(cbind(Y95_96,Y96_97,Y97_98,Y98_99,Y99_00,Y00_01,Y01_02,Y02_03,Y03_04,Y04_05,Y05_06,Y06_07,Y07_08,Y08_09,Y09_10,Y10_11,Y11_12,Y12_13,Y13_14,Y14_15,Y15_16,Y16_17,Y17_18,Y18_19,Y19_20), MARGIN=1, FUN=sd, na.rm=TRUE)
df<-as.data.frame(cbind(df,"NDVI_sdY19_20"=sdY19_20,"NDVI_sdY18_20"=sdY18_20,"NDVI_sdY17_20"=sdY17_20,"NDVI_sdY16_20"=sdY16_20,"NDVI_sdY15_20"=sdY15_20,
                        "NDVI_sdY14_20"=sdY14_20,"NDVI_sdY13_20"=sdY13_20,"NDVI_sdY12_20"=sdY12_20,"NDVI_sdY11_20"=sdY11_20,"NDVI_sdY10_20"=sdY10_20,
                        "NDVI_sdY09_20"=sdY09_20,"NDVI_sdY08_20"=sdY08_20,"NDVI_sdY07_20"=sdY07_20,"NDVI_sdY06_20"=sdY06_20,"NDVI_sdY05_20"=sdY05_20,
                        "NDVI_sdY04_20"=sdY04_20,"NDVI_sdY03_20"=sdY03_20,"NDVI_sdY02_20"=sdY02_20,"NDVI_sdY01_20"=sdY01_20,"NDVI_sdY00_20"=sdY00_20,
                        "NDVI_sdY99_20"=sdY99_20,"NDVI_sdY98_20"=sdY98_20,"NDVI_sdY97_20"=sdY97_20,"NDVI_sdY96_20"=sdY96_20,"NDVI_sdY95_20"=sdY95_20))
NDVI<-KalmarDataNDVI
#Clean environment
rm(sdY19_20,sdY18_20,sdY17_20,sdY16_20,sdY15_20,sdY14_20,sdY13_20,sdY12_20,sdY11_20,sdY10_20,sdY09_20,
   sdY08_20,sdY07_20,sdY06_20,sdY05_20,sdY04_20,sdY03_20,sdY02_20,sdY01_20,sdY00_20,sdY99_20,sdY98_20,
   sdY97_20,sdY96_20,sdY95_20,Y19_20,Y18_19,Y17_18,Y16_17,Y15_16,Y14_15,Y13_14,Y12_13,Y11_12,Y10_11,
   Y09_10,Y08_09,Y07_08,Y06_07,Y05_06,Y04_05,Y03_04,Y02_03,Y01_02,Y00_01,Y95_96,Y98_99,Y99_00,Y96_97,
   Y97_98,KalmarDataNDVI)

#Clean Landsat NBR time series data
Index<- read.csv("Excel_sheets/Landsat_SR_NBR.csv")
rownames(Index)<-Index$Site_ID
#Remove unwanted Columns
Index$system.index=NULL
Index$.geo=NULL
Index$Site_ID=NULL
#Rename column headers sort matrix
colnames(Index) <- str_sub((colnames(Index)),-8,-1)
Index<-Index[order(rownames(Index)),]
Index<-Index[ , order(names(Index))]
#Remove multiples by calculating their means
Index<-sapply(split.default(Index, sub('\\..*', '', names(Index))), rowMeans, na.rm = TRUE)
#Select imagery from 1995 - 2020
Index<-Index[,93:1010]
#Remove outliers
sum(is.na(Index))
x <- abs(Index - median(Index,na.rm = T)) / mad(Index,na.rm = T)
tr <- ifelse(quantile(x,.99,na.rm=T)>2,quantile(x,.99,na.rm=T),NA)
if(!is.na(tr)){
  Index[ which(x > tr )] <- NA
}
sum(is.na(Index))
#Clean environment
rm(x,tr)
#Impute missing data using Kalman smoother and apply sgolayfilt Savitzky-Golay filtering
IndexT<-t(Index)
ColNames<-colnames(IndexT)
RowNames<-rownames(IndexT)
Site_00<-sgolayfilt(na_kalman(IndexT[,1],model = "auto.arima"))
Site_01<-sgolayfilt(na_kalman(IndexT[,2],model = "auto.arima"))
Site_02<-sgolayfilt(na_kalman(IndexT[,3],model = "auto.arima"))
Site_03<-sgolayfilt(na_kalman(IndexT[,4],model = "auto.arima"))
Site_04<-sgolayfilt(na_kalman(IndexT[,5],model = "auto.arima"))
Site_05<-sgolayfilt(na_kalman(IndexT[,6],model = "auto.arima"))
Site_06<-sgolayfilt(na_kalman(IndexT[,7],model = "auto.arima"))
Site_07<-sgolayfilt(na_kalman(IndexT[,8],model = "auto.arima"))
Site_08<-sgolayfilt(na_kalman(IndexT[,9],model = "auto.arima"))
Site_09<-sgolayfilt(na_kalman(IndexT[,10],model = "auto.arima"))
Site_10<-sgolayfilt(na_kalman(IndexT[,11],model = "auto.arima"))
Site_11<-sgolayfilt(na_kalman(IndexT[,12],model = "auto.arima"))
Site_12<-sgolayfilt(na_kalman(IndexT[,13],model = "auto.arima"))
Site_13<-sgolayfilt(na_kalman(IndexT[,14],model = "auto.arima"))
Site_14<-sgolayfilt(na_kalman(IndexT[,15],model = "auto.arima"))
Site_15<-sgolayfilt(na_kalman(IndexT[,16],model = "auto.arima"))
Site_16<-sgolayfilt(na_kalman(IndexT[,17],model = "auto.arima"))
Site_17<-sgolayfilt(na_kalman(IndexT[,18],model = "auto.arima"))
Site_18<-sgolayfilt(na_kalman(IndexT[,19],model = "auto.arima"))
Site_19<-sgolayfilt(na_kalman(IndexT[,20],model = "auto.arima"))
Site_20<-sgolayfilt(na_kalman(IndexT[,21],model = "auto.arima"))
Site_21<-sgolayfilt(na_kalman(IndexT[,22],model = "auto.arima"))
Site_22<-sgolayfilt(na_kalman(IndexT[,23],model = "auto.arima"))
Site_23<-sgolayfilt(na_kalman(IndexT[,24],model = "auto.arima"))
Site_24<-sgolayfilt(na_kalman(IndexT[,25],model = "auto.arima"))
Site_25<-sgolayfilt(na_kalman(IndexT[,26],model = "auto.arima"))
Site_26<-sgolayfilt(na_kalman(IndexT[,27],model = "auto.arima"))
Site_27<-sgolayfilt(na_kalman(IndexT[,28],model = "auto.arima"))
Site_28<-sgolayfilt(na_kalman(IndexT[,29],model = "auto.arima"))
Site_29<-sgolayfilt(na_kalman(IndexT[,30],model = "auto.arima"))
Site_30<-sgolayfilt(na_kalman(IndexT[,31],model = "auto.arima"))
Site_31<-sgolayfilt(na_kalman(IndexT[,32],model = "auto.arima"))
Site_32<-sgolayfilt(na_kalman(IndexT[,33],model = "auto.arima"))
Site_33<-sgolayfilt(na_kalman(IndexT[,34],model = "auto.arima"))
Site_34<-sgolayfilt(na_kalman(IndexT[,35],model = "auto.arima"))
Site_35<-sgolayfilt(na_kalman(IndexT[,36],model = "auto.arima"))
Site_36<-sgolayfilt(na_kalman(IndexT[,37],model = "auto.arima"))
Site_37<-sgolayfilt(na_kalman(IndexT[,38],model = "auto.arima"))
Site_38<-sgolayfilt(na_kalman(IndexT[,39],model = "auto.arima"))
Site_39<-sgolayfilt(na_kalman(IndexT[,40],model = "auto.arima"))
Site_40<-sgolayfilt(na_kalman(IndexT[,41],model = "auto.arima"))
Site_41<-sgolayfilt(na_kalman(IndexT[,42],model = "auto.arima"))
Site_42<-sgolayfilt(na_kalman(IndexT[,43],model = "auto.arima"))
Site_43<-sgolayfilt(na_kalman(IndexT[,44],model = "auto.arima"))
Site_44<-sgolayfilt(na_kalman(IndexT[,45],model = "auto.arima"))
Site_45<-sgolayfilt(na_kalman(IndexT[,46],model = "auto.arima"))
Site_46<-sgolayfilt(na_kalman(IndexT[,47],model = "auto.arima"))
Site_47<-sgolayfilt(na_kalman(IndexT[,48],model = "auto.arima"))
Site_48<-sgolayfilt(na_kalman(IndexT[,49],model = "auto.arima"))
Site_49<-sgolayfilt(na_kalman(IndexT[,50],model = "auto.arima"))
Site_50<-sgolayfilt(na_kalman(IndexT[,51],model = "auto.arima"))
KalmarDataNBR<-cbind(Site_00,Site_01,Site_02,Site_03,Site_04,Site_05,Site_06,Site_07,Site_08,
                     Site_09,Site_10,Site_11,Site_12,Site_13,Site_14,Site_15,Site_16,Site_17,
                     Site_18,Site_19,Site_20,Site_21,Site_22,Site_23,Site_24,Site_25,Site_26,
                     Site_27,Site_28,Site_29,Site_30,Site_31,Site_32,Site_33,Site_34,Site_35,
                     Site_36,Site_37,Site_38,Site_39,Site_40,Site_41,Site_42,Site_43,Site_44,
                     Site_45,Site_46,Site_47,Site_48,Site_49,Site_50)
colnames(KalmarDataNBR)<-ColNames
rownames(KalmarDataNBR)<-RowNames
NBR<-as.data.frame(t(KalmarDataNBR))
rm(Site_00,Site_01,Site_02,Site_03,Site_04,Site_05,Site_06,Site_07,Site_08,
   Site_09,Site_10,Site_11,Site_12,Site_13,Site_14,Site_15,Site_16,Site_17,
   Site_18,Site_19,Site_20,Site_21,Site_22,Site_23,Site_24,Site_25,Site_26,
   Site_27,Site_28,Site_29,Site_30,Site_31,Site_32,Site_33,Site_34,Site_35,
   Site_36,Site_37,Site_38,Site_39,Site_40,Site_41,Site_42,Site_43,Site_44,
   Site_45,Site_46,Site_47,Site_48,Site_49,Site_50,ColNames,RowNames,Index,
   IndexT,KalmarDataNBR)

#Load and prepare Terraclimate varaibles
#Min Temperature
Terra_tmmn<- read.csv("Excel_sheets/Terra_tmmn.csv")
rownames(Terra_tmmn)<-Terra_tmmn$Site_ID
Terra_tmmn$system.index=NULL
Terra_tmmn$.geo=NULL
Terra_tmmn$Site_ID=NULL
colnames(Terra_tmmn) <- str_sub((colnames(Terra_tmmn)),2,7)
Terra_tmmn<-Terra_tmmn[order(rownames(Terra_tmmn)),]
Terra_tmmn<-Terra_tmmn[ , order(names(Terra_tmmn))]
#Precipitation
Terra_pr<- read.csv("Excel_sheets/Terra_pr.csv")
rownames(Terra_pr)<-Terra_pr$Site_ID
Terra_pr$system.index=NULL
Terra_pr$.geo=NULL
Terra_pr$Site_ID=NULL
colnames(Terra_pr) <- str_sub((colnames(Terra_pr)),2,7)
Terra_pr<-Terra_pr[order(rownames(Terra_pr)),]
Terra_pr<-Terra_pr[ , order(names(Terra_pr))]
# Solar Radiation
Terra_rad<- read.csv("Excel_sheets/Terra_rad.csv")
rownames(Terra_rad)<-Terra_rad$Site_ID
Terra_rad$system.index=NULL
Terra_rad$.geo=NULL
Terra_rad$Site_ID=NULL
colnames(Terra_rad) <- str_sub((colnames(Terra_rad)),2,7)
Terra_rad<-Terra_rad[order(rownames(Terra_rad)),]
Terra_rad<-Terra_rad[ , order(names(Terra_rad))]

########### Lag models ###########
#Standerdise
df_stand<-standardize(df[,5:31])
df_stand<-as.data.frame(cbind("Long_X"=df$Long_X,"Lat_Y"=df$Lat_Y,"Plantation"=df$Plantation,df_stand))
#Prepair data for looping models and create empty varibles
number=1
#Response
out_start=4
out_end= 5
out_nvar=out_end-out_start+1
out_variable=rep(NA, out_nvar)
out_beta=rep(NA, out_nvar)
out_se = rep(NA, out_nvar)
out_5IntVal=rep(NA, out_nvar)
out_95IntVal=rep(NA, out_nvar)
out_pvalue=rep(NA, out_nvar)
out_AIC=rep(NA, out_nvar)
out_tvalue=rep(NA, out_nvar)
#Explanatory
exp_start=6
exp_end=30
exp_nvar=exp_end-exp_start+1
exp_variable=rep(NA, exp_nvar)
exp_beta=rep(NA, exp_nvar)
exp_se = rep(NA, out_nvar)
exp_5IntVal=rep(NA, exp_nvar)
exp_95IntVal=rep(NA, exp_nvar)
exp_pvalue=rep(NA, exp_nvar)
exp_AIC=rep(NA, exp_nvar)
exp_tvalue=rep(NA, out_nvar)
#Loop models
for (i in out_start:out_end){
  outcome = colnames(df_stand)[i]
  for (j in exp_start:exp_end){
    exposure = colnames(df_stand)[j]
    
    model <-lmer(get(outcome)~get(exposure)+(1|Plantation),data=df_stand)
    
    Vcov <- vcov(model, useScale = FALSE)
    beta <- fixef(model)
    se <- sqrt(diag(Vcov))
    zval <- beta / se
    pval <- 2 * pnorm(abs(zval), lower.tail = FALSE)
    IntVal <- as.numeric(confint(model))
    AIC<-as.numeric(AIC(model))
    t_value<-as.numeric(coef(summary(model))[,"t value"])
    
    out_beta[number] = as.numeric(beta[2])
    out_se[number] = as.numeric(se[2])
    out_5IntVal[number] = IntVal[4]
    out_95IntVal[number] = IntVal[8]
    out_pvalue[number] = as.numeric(pval[2])
    out_variable[number] = outcome
    out_AIC[number] = AIC
    out_tvalue[number] = t_value[2]
    number = number + 1
    
    exp_beta[number] = as.numeric(beta[2])
    exp_se[number] = as.numeric(se[2])
    exp_5IntVal[number] = IntVal[4]
    exp_95IntVal[number] = IntVal[8]
    exp_pvalue[number] = as.numeric(pval[2])
    exp_variable[number] = exposure
    exp_AIC[number] = AIC
    exp_tvalue[number] = t_value[2]
    number = number + 1
  }
}
#Create a dataframe with results
outcome = data.frame(out_variable, out_beta, out_se, out_5IntVal, out_95IntVal,out_tvalue, out_pvalue, out_AIC)
exposure = data.frame(exp_variable, exp_beta, exp_se, exp_5IntVal, exp_95IntVal,exp_tvalue, exp_pvalue, exp_AIC)
outcome = outcome %>% 
  rename(
    variable = out_variable,
    beta = out_beta,
    se = out_se,
    IntVal5 = out_5IntVal,
    IntVal95 = out_95IntVal,
    p_value = out_pvalue,
    AIC = out_AIC,
    t_value = out_tvalue
  )
exposure = exposure %>% 
  rename(
    variable = exp_variable,
    beta = exp_beta,
    se = exp_se,
    IntVal5 = exp_5IntVal,
    IntVal95 = exp_95IntVal,
    p_value = exp_pvalue,
    AIC = exp_AIC,
    t_value = exp_tvalue
  )
exposure = na.omit(exposure)
outcome = na.omit(outcome)
write.csv(exposure,"R_code/Temp/exposure.csv",row.names = FALSE)
write.csv(outcome,"R_code/Temp/outcome.csv",row.names = FALSE)
exposure<- read.csv("R_code/Temp/exposure.csv",stringsAsFactors = F)
outcome<- read.csv("R_code/Temp/outcome.csv",stringsAsFactors = F)
colnames(exposure)<-c("Explanatory","beta","se","5%","95%","t_value","p_value","AIC")
all<- cbind("Response"=outcome$variable,exposure)
#Find significant results
Significant <- dplyr::filter(all, p_value <=0.05)
#Clean environment
rm(all,exposure,model,outcome,Vcov,AIC, beta,exp_5IntVal,exp_95IntVal,exp_AIC,exp_end,exp_beta,exp_nvar,
   exp_pvalue,exp_se,exp_start,exp_variable,i,IntVal,j,number,out_AIC,out_end,out_nvar,out_pvalue,
   out_start,out_se,out_beta,out_variable,out_5IntVal,out_95IntVal,pval,se,zval,out_tvalue,exp_tvalue,t_value)

#Model 1 validation and plotting
model1<-lmer(Richness~NDVI_sdY19_20+(1|Plantation),data=df_stand)
summary(model1)
plot(model1)
qqnorm(resid(model1))
qqline(resid(model1))
#Plotting standerdise residuals
S_Residuals<-simulateResiduals(fittedModel=model1)
plot(S_Residuals)
#Despersion testing
testDispersion(S_Residuals)
#Visualize
ggplot(data=df,aes(x=NDVI_sdY19_20,y=Richness))+
  geom_point()+
  geom_smooth(method="lm",colour="black")+
  labs(y="Caelifera species richness\n",x="\nDeviation of NDVI between 2019-2020")+
  theme(axis.title=element_text(size=14),
        axis.text=element_text(size=13,colour="black"),
        axis.ticks=element_line(size=0.8,colour="black"),
        axis.line=element_line(size=0.8,colour="black"),
        panel.grid=element_blank(),
        panel.background=element_blank())
#Model 2 validation and plotting
model2<-lmer(exShannon~NDVI_sdY19_20+(1|Plantation),data=df_stand)
plot(model2)
qqnorm(resid(model2))
qqline(resid(model2))
#Plotting standerdise residuals
S_Residuals<-simulateResiduals(fittedModel=model2)
plot(S_Residuals)
#Despersion testing
testDispersion(S_Residuals)
#Visualize
ggplot(data=df, aes(x =NDVI_sdY19_20, y = exShannon)) +
  geom_point() +
  geom_smooth(method="lm",colour="black")+
  labs(y="Caelifera exponent Shannon diversity\n",x="\nDeviation of NDVI between 2019-2020")+
  theme(axis.title=element_text(size=14),
        axis.text=element_text(size=13,colour="black"),
        axis.ticks=element_line(size=0.8,colour="black"),
        axis.line=element_line(size=0.8,colour="black"),
        panel.grid=element_blank(),
        panel.background=element_blank())
#Model 3 validation and plotting
model3<-lmer(exShannon~NDVI_sdY17_20+(1|Plantation),data=df_stand)
plot(model3)
qqnorm(resid(model3))
qqline(resid(model3))
#Plotting standerdise residuals
S_Residuals<-simulateResiduals(fittedModel=model3)
plot(S_Residuals)
#Despersion testing
testDispersion(S_Residuals)
#Visualize
ggplot(data=df, aes(x =NDVI_sdY17_20, y = exShannon)) +
  geom_point() +
  geom_smooth(method="lm",colour="black")+
  labs(y="Caelifera exponent Shannon diversity\n",x="\nDeviation of NDVI between 2017-2020")+
  theme(axis.title=element_text(size=14),
        axis.text=element_text(size=13,colour="black"),
        axis.ticks=element_line(size=0.8,colour="black"),
        axis.line=element_line(size=0.8,colour="black"),
        panel.grid=element_blank(),
        panel.background=element_blank())
#Clean environment
rm(Significant,S_Residuals,model1,model2,model3)

########### Variation partitioning ###########
#Select time stamps of significant variation from Q1
NDVI17_20<-NDVI %>% dplyr::select(starts_with(c("201704","201705","201706","201707","201708","201709","201710","201711","201712","201801","201802","201803","201804","201805","201806","201807","201808","201809","201810","201811","201812","201901","201902","201903","201904","201905","201906","201907","201908","201909","201910","201911","201912","202001","202002","202003")))
NBR17_20<-NBR %>% dplyr::select(starts_with(c("201704","201705","201706","201707","201708","201709","201710","201711","201712","201801","201802","201803","201804","201805","201806","201807","201808","201809","201810","201811","201812","201901","201902","201903","201904","201905","201906","201907","201908","201909","201910","201911","201912","202001","202002","202003")))
Terra_tmmn17_20<-Terra_tmmn %>% dplyr::select(starts_with(c("201704","201705","201706","201707","201708","201709","201710","201711","201712","201801","201802","201803","201804","201805","201806","201807","201808","201809","201810","201811","201812","201901","201902","201903","201904","201905","201906","201907","201908","201909","201910","201911","201912","202001","202002","202003")))
Terra_pr17_20<-Terra_pr %>% dplyr::select(starts_with(c("201704","201705","201706","201707","201708","201709","201710","201711","201712","201801","201802","201803","201804","201805","201806","201807","201808","201809","201810","201811","201812","201901","201902","201903","201904","201905","201906","201907","201908","201909","201910","201911","201912","202001","202002","202003")))
Terra_rad17_20<-Terra_rad %>% dplyr::select(starts_with(c("201704","201705","201706","201707","201708","201709","201710","201711","201712","201801","201802","201803","201804","201805","201806","201807","201808","201809","201810","201811","201812","201901","201902","201903","201904","201905","201906","201907","201908","201909","201910","201911","201912","202001","202002","202003")))
#Perform PCA to reduce dimentionality
#NDVI
NDVI.pca <- prcomp(scale(t(NDVI17_20)),center = TRUE)
summary(NDVI.pca)
NDVI.pca<-NDVI.pca$rotation[,1:2]
#NBR
NBR.pca <- prcomp(scale(t(NBR17_20)),center = TRUE)
summary(NBR.pca)
NBR.pca<-NBR.pca$rotation[,1:2]
#Min Temp
Terra_tmmn.pca <- prcomp(scale(t(Terra_tmmn17_20)),center = TRUE)
summary(Terra_tmmn.pca)
Terra_tmmn.pca<-Terra_tmmn.pca$rotation[,1]
#Precipitation
Terra_pr.pca <- prcomp(scale(t(Terra_pr17_20)),center = TRUE)
summary(Terra_pr.pca)
Terra_pr.pca<-Terra_pr.pca$rotation[,1]
#Solar radiation
Terra_rad.pca <- prcomp(scale(t(Terra_rad17_20)),center = TRUE)
summary(Terra_rad.pca)
Terra_rad.pca<-Terra_rad.pca$rotation[,1]
#Variation partitioning
mod<-varpart(NDVI.pca, NBR.pca, Terra_tmmn.pca, Terra_pr.pca, Terra_rad.pca)
mod
plot(mod,Xnames = c("NBR","Temp","Precip","Sol Rad"),
     bg = c(2,4,5,3))

#Test significance
#All variables
TestAll<-rda(NDVI.pca~NBR.pca+Terra_tmmn.pca+Terra_pr.pca+Terra_rad.pca)
anova(TestAll)
RsquareAdj(TestAll)
#NBR
TestNBR.pca<-rda(NDVI.pca~NBR.pca)
anova(TestNBR.pca)
RsquareAdj(TestNBR.pca)             
#Min Temp
TestTerra_tmmn.pca<-rda(NDVI.pca~Terra_tmmn.pca)
anova(TestTerra_tmmn.pca)
RsquareAdj(TestTerra_tmmn.pca)             
#Precipitation
TestTerra_pr.pca<-rda(NDVI.pca~Terra_pr.pca)
anova(TestTerra_pr.pca)
RsquareAdj(TestTerra_pr.pca)             
#Solar radiation
TestTerra_rad.pca<-rda(NDVI.pca~Terra_rad.pca)
anova(TestTerra_rad.pca)
RsquareAdj(TestTerra_rad.pca)             
#Clean environment
rm(Terra_pr,Terra_pr.pca,Terra_pr17_20,Terra_tmmn,Terra_tmmn.pca,Terra_rad,Terra_rad.pca,
   Terra_rad17_20,Terra_tmmn17_20,mod,NBR,NBR.pca,NBR17_20,NDVI,NDVI.pca,NDVI17_20,TestAll,
   TestTerra_pr.pca,TestTerra_rad.pca,TestTerra_tmmn.pca,TestNBR.pca)

########### Spectral models ###########
#Unfortunatly, PlanetScope imagery or derived products cannot be shared publicly.
'''
#Load and extract Planetscope imagery 2019/05
Estate1_2019_05<-brick("PlanetScope/2019_05/Estate1/20190502_075144_64_105f_3B_AnalyticMS_SR.tif")
Estate2_2019_05<-brick("PlanetScope/2019_05/Estate2/20190502_075142_57_105f_3B_AnalyticMS_SR.tif")
Estate3_2019_05<-brick("PlanetScope/2019_05/Estate3/20190503_093059_13_106a_3B_AnalyticMS_SR.tif")
Estate4_1_2019_05<-brick("PlanetScope/2019_05/Estate4/20190510_075040_15_105a_3B_AnalyticMS_SR.tif")
Estate4_2_2019_05<-brick("PlanetScope/2019_05/Estate4/20190510_075042_22_105a_3B_AnalyticMS_SR.tif")
#Calculating NDVI    (NIR - RED) / (NIR + RED)
NDVI_Estate1_2019_05<-(Estate1_2019_05[[4]]-Estate1_2019_05[[3]])/(Estate1_2019_05[[4]]+Estate1_2019_05[[3]])
NDVI_Estate2_2019_05<-(Estate2_2019_05[[4]]-Estate2_2019_05[[3]])/(Estate2_2019_05[[4]]+Estate2_2019_05[[3]])
NDVI_Estate3_2019_05<-(Estate3_2019_05[[4]]-Estate3_2019_05[[3]])/(Estate3_2019_05[[4]]+Estate3_2019_05[[3]])
NDVI_Estate4_1_2019_05<-(Estate4_1_2019_05[[4]]-Estate4_1_2019_05[[3]])/(Estate4_1_2019_05[[4]]+Estate4_1_2019_05[[3]])
NDVI_Estate4_2_2019_05<-(Estate4_2_2019_05[[4]]-Estate4_2_2019_05[[3]])/(Estate4_2_2019_05[[4]]+Estate4_2_2019_05[[3]])
#Clean environment
rm(Estate3_2019_05,Estate2_2019_05,Estate1_2019_05,Estate4_1_2019_05,Estate4_2_2019_05)
#Extraxting values 
GPS_Loc<-readOGR("PlanetScope/Site_GPS_coordinates_Projected.shp",stringsAsFactors = FALSE)
GPS_Buf_50<-gBuffer(GPS_Loc,byid=TRUE,width=50,capStyle="ROUND")
NDVI_Estate1_2019_05_50<-zonal.stats(GPS_Buf_50,NDVI_Estate1_2019_05,stat=c("mean"))
NDVI_Estate2_2019_05_50<-zonal.stats(GPS_Buf_50,NDVI_Estate2_2019_05,stat=c("mean"))
NDVI_Estate3_2019_05_50<-zonal.stats(GPS_Buf_50,NDVI_Estate3_2019_05,stat=c("mean"))
NDVI_Estate4_1_2019_05_50<-zonal.stats(GPS_Buf_50,NDVI_Estate4_1_2019_05,stat=c("mean"))
NDVI_Estate4_2_2019_05_50<-zonal.stats(GPS_Buf_50,NDVI_Estate4_2_2019_05,stat=c("mean"))
#Clean environment
rm(NDVI_Estate1_2019_05,NDVI_Estate2_2019_05,NDVI_Estate3_2019_05,NDVI_Estate4_1_2019_05,NDVI_Estate4_2_2019_05)
#Prepare data frame
ndvi_2019_05_50<-data.frame("NDVI_Estate1_50"=NDVI_Estate1_2019_05_50,
                            "NDVI_Estate2_50"=NDVI_Estate2_2019_05_50,
                            "NDVI_Estate3_50"=NDVI_Estate3_2019_05_50,
                            "NDVI_Estate4_1_50"=NDVI_Estate4_1_2019_05_50,
                            "NDVI_Estate4_2_50"=NDVI_Estate4_2_2019_05_50)
PlanetScope_NDVI_50_2019_05<-apply(ndvi_2019_05_50,MARGIN=1,FUN=mean,na.rm=TRUE)
#Clean environment
rm(ndvi_2019_05_50,NDVI_Estate1_2019_05_50,NDVI_Estate2_2019_05_50,NDVI_Estate3_2019_05_50,NDVI_Estate4_1_2019_05_50,NDVI_Estate4_2_2019_05_50)

#Load and extract Planetscope imagery 2020/01
Estate1_2020_01<-brick("PlanetScope/2020_01/Estate1/20200124_080656_10_1057_3B_AnalyticMS_SR.tif")
Estate2_2020_01<-brick("PlanetScope/2020_01/Estate2/20200129_091654_99_1065_3B_AnalyticMS_SR.tif")
Estate3_2020_01<-brick("PlanetScope/2020_01/Estate3/20200117_080646_31_1064_3B_AnalyticMS_SR.tif")
Estate4_2020_01<-brick("PlanetScope/2020_01/Estate4/20200215_080802_12_1057_3B_AnalyticMS_SR.tif")
#Calculating NDVI    (NIR - RED) / (NIR + RED)
NDVI_Estate1_2020_01<-(Estate1_2020_01[[4]]-Estate1_2020_01[[3]])/(Estate1_2020_01[[4]]+Estate1_2020_01[[3]])
NDVI_Estate2_2020_01<-(Estate2_2020_01[[4]]-Estate2_2020_01[[3]])/(Estate2_2020_01[[4]]+Estate2_2020_01[[3]])
NDVI_Estate3_2020_01<-(Estate3_2020_01[[4]]-Estate3_2020_01[[3]])/(Estate3_2020_01[[4]]+Estate3_2020_01[[3]])
NDVI_Estate4_2020_01<-(Estate4_2020_01[[4]]-Estate4_2020_01[[3]])/(Estate4_2020_01[[4]]+Estate4_2020_01[[3]])
#Clean environment
rm(Estate3_2020_01,Estate2_2020_01,Estate1_2020_01,Estate4_2020_01)
#Extraxting values 
NDVI_Estate1_2020_01_50<-zonal.stats(GPS_Buf_50,NDVI_Estate1_2020_01,stat=c("mean"))
NDVI_Estate2_2020_01_50<-zonal.stats(GPS_Buf_50,NDVI_Estate2_2020_01,stat=c("mean"))
NDVI_Estate3_2020_01_50<-zonal.stats(GPS_Buf_50,NDVI_Estate3_2020_01,stat=c("mean"))
NDVI_Estate4_2020_01_50<-zonal.stats(GPS_Buf_50,NDVI_Estate4_2020_01,stat=c("mean"))
#Clean environment
rm(NDVI_Estate1_2020_01,NDVI_Estate2_2020_01,NDVI_Estate3_2020_01,NDVI_Estate4_2020_01)
#Prepare data frame
ndvi_2020_01_50<-data.frame("NDVI_Estate1_2020_01_50"=NDVI_Estate1_2020_01_50,
                            "NDVI_Estate2_2020_01_50"=NDVI_Estate2_2020_01_50,
                            "NDVI_Estate3_2020_01_50"=NDVI_Estate3_2020_01_50,
                            "NDVI_Estate4_2020_01_50"=NDVI_Estate4_2020_01_50)
PlanetScope_NDVI_50_2020_01<-apply(ndvi_2020_01_50,MARGIN=1,FUN=mean,na.rm=TRUE)
#Clean environment
rm(ndvi_2020_01_50,NDVI_Estate1_2020_01_50,NDVI_Estate2_2020_01_50,NDVI_Estate3_2020_01_50,NDVI_Estate4_2020_01_50)

#Load and extract Planetscope imagery 2020/03
Estate1_1_2020_03<-brick("PlanetScope/2020_03/Estate1/20200323_091450_19_106a_3B_AnalyticMS_SR.tif")
Estate1_2_2020_03<-brick("PlanetScope/2020_03/Estate1/20200323_091452_29_106a_3B_AnalyticMS_SR.tif")
Estate2_2020_03<-brick("PlanetScope/2020_03/Estate2/20200319_091226_84_106a_3B_AnalyticMS_SR.tif")
Estate3_2020_03<-brick("PlanetScope/2020_03/Estate3/20200318_091212_19_106e_3B_AnalyticMS_SR.tif")
Estate4_2020_03<-brick("PlanetScope/2020_03/Estate4/20200320_075746_13_105c_3B_AnalyticMS_SR.tif")
#Calculating NDVI    (NIR - RED) / (NIR + RED)
NDVI_Estate1_1_2020_03<-(Estate1_1_2020_03[[4]]-Estate1_1_2020_03[[3]])/(Estate1_1_2020_03[[4]]+Estate1_1_2020_03[[3]])
NDVI_Estate1_2_2020_03<-(Estate1_2_2020_03[[4]]-Estate1_2_2020_03[[3]])/(Estate1_2_2020_03[[4]]+Estate1_2_2020_03[[3]])
NDVI_Estate2_2020_03<-(Estate2_2020_03[[4]]-Estate2_2020_03[[3]])/(Estate2_2020_03[[4]]+Estate2_2020_03[[3]])
NDVI_Estate3_2020_03<-(Estate3_2020_03[[4]]-Estate3_2020_03[[3]])/(Estate3_2020_03[[4]]+Estate3_2020_03[[3]])
NDVI_Estate4_2020_03<-(Estate4_2020_03[[4]]-Estate4_2020_03[[3]])/(Estate4_2020_03[[4]]+Estate4_2020_03[[3]])
#Clean environment
rm(Estate3_2020_03,Estate2_2020_03,Estate1_1_2020_03,Estate1_2_2020_03,Estate4_2020_03)
#Extraxting values 
NDVI_Estate1_1_2020_03_50<-zonal.stats(GPS_Buf_50,NDVI_Estate1_1_2020_03,stat=c("mean"))
NDVI_Estate1_2_2020_03_50<-zonal.stats(GPS_Buf_50,NDVI_Estate1_2_2020_03,stat=c("mean"))
NDVI_Estate2_2020_03_50<-zonal.stats(GPS_Buf_50,NDVI_Estate2_2020_03,stat=c("mean"))
NDVI_Estate3_2020_03_50<-zonal.stats(GPS_Buf_50,NDVI_Estate3_2020_03,stat=c("mean"))
NDVI_Estate4_2020_03_50<-zonal.stats(GPS_Buf_50,NDVI_Estate4_2020_03,stat=c("mean"))
#Clean environment
rm(NDVI_Estate1_1_2020_03,NDVI_Estate1_2_2020_03,NDVI_Estate2_2020_03,NDVI_Estate3_2020_03,NDVI_Estate4_2020_03)
#Prepare data frame
ndvi_2020_03_50<-data.frame("NDVI_Estate1_1_2020_03_50"=NDVI_Estate1_1_2020_03_50,
                            "NDVI_Estate1_2_2020_03_50"=NDVI_Estate1_2_2020_03_50,
                            "NDVI_Estate2_2020_03_50"=NDVI_Estate2_2020_03_50,
                            "NDVI_Estate3_2020_03_50"=NDVI_Estate3_2020_03_50,
                            "NDVI_Estate4_2020_03_50"=NDVI_Estate4_2020_03_50)
PlanetScope_NDVI_50__2020_03<-apply(ndvi_2020_03_50,MARGIN=1,FUN=mean,na.rm=TRUE)
#Clean environment
rm(ndvi_2020_03_50,NDVI_Estate1_1_2020_03_50,NDVI_Estate1_2_2020_03_50,NDVI_Estate2_2020_03_50,NDVI_Estate3_2020_03_50,NDVI_Estate4_2020_03_50)
#Merge PlanetScope data to response data frame
planet_NDVI<-data.frame("PlanetScope_NDVI_50_2019_05"=PlanetScope_NDVI_50_2019_05,
                        "PlanetScope_NDVI_50_2020_01"=PlanetScope_NDVI_50_2020_01,
                        "PlanetScope_NDVI_50__2020_03"=PlanetScope_NDVI_50__2020_03)
df<-as.data.frame(cbind(df,planet_NDVI))
#Clean environment
rm(GPS,GPS_Loc,GPS_Buf_50,planet_NDVI,PlanetScope_NDVI_50_2019_05,PlanetScope_NDVI_50_2020_01,PlanetScope_NDVI_50__2020_03)
'''

#Adding Google Earth Engine imagery
#Sentinel NDVI 20190524
SentinelNDVI20190524<- read.csv("Excel_sheets/Sentinel_NDVI_20190524.csv")
SentinelNDVI20190524<-SentinelNDVI20190524[order(SentinelNDVI20190524$Site_ID),]
rownames(SentinelNDVI20190524)<-SentinelNDVI20190524$Site_ID
SentinelNDVI20190524$system.index=NULL
SentinelNDVI20190524$.geo=NULL
SentinelNDVI20190524$Site_ID=NULL
#Sentinel NDVI 20200104
SentinelNDVI20200104<- read.csv("Excel_sheets/Sentinel_NDVI_20200104.csv")
SentinelNDVI20200104<-SentinelNDVI20200104[order(SentinelNDVI20200104$Site_ID),]
rownames(SentinelNDVI20200104)<-SentinelNDVI20200104$Site_ID
SentinelNDVI20200104$system.index=NULL
SentinelNDVI20200104$.geo=NULL
SentinelNDVI20200104$Site_ID=NULL
#Sentinel NDVI 20200304
SentinelNDVI20200304<- read.csv("Excel_sheets/Sentinel_NDVI_20200304.csv")
SentinelNDVI20200304<-SentinelNDVI20200304[order(SentinelNDVI20200304$Site_ID),]
rownames(SentinelNDVI20200304)<-SentinelNDVI20200304$Site_ID
SentinelNDVI20200304$system.index=NULL
SentinelNDVI20200304$.geo=NULL
SentinelNDVI20200304$Site_ID=NULL
#Landsat NDVI 20190524
LandsatNDVI20190524<- read.csv("Excel_sheets/Landsat_NDVI_20190524.csv")
LandsatNDVI20190524<-LandsatNDVI20190524[order(LandsatNDVI20190524$Site_ID),]
rownames(LandsatNDVI20190524)<-LandsatNDVI20190524$Site_ID
LandsatNDVI20190524$system.index=NULL
LandsatNDVI20190524$.geo=NULL
LandsatNDVI20190524$Site_ID=NULL
#Landsat NDVI 20200103
LandsatNDVI20200103<- read.csv("Excel_sheets/Landsat_NDVI_20200103.csv")
LandsatNDVI20200103<-LandsatNDVI20200103[order(LandsatNDVI20200103$Site_ID),]
rownames(LandsatNDVI20200103)<-LandsatNDVI20200103$Site_ID
LandsatNDVI20200103$system.index=NULL
LandsatNDVI20200103$.geo=NULL
LandsatNDVI20200103$Site_ID=NULL
#Landsat NDVI 20200307
LandsatNDVI20200307<- read.csv("Excel_sheets/Landsat_NDVI_20200307.csv")
LandsatNDVI20200307<-LandsatNDVI20200307[order(LandsatNDVI20200307$Site_ID),]
rownames(LandsatNDVI20200307)<-LandsatNDVI20200307$Site_ID
LandsatNDVI20200307$system.index=NULL
LandsatNDVI20200307$.geo=NULL
LandsatNDVI20200307$Site_ID=NULL
# Construct data frame
df<-as.data.frame(cbind(df,
                        "Sentinel_NDVI_2019_05"=SentinelNDVI20190524$X0,
                        "Sentinel_NDVI_2020_01"=SentinelNDVI20200104$X0,
                        "Sentinel_NDVI_2020_03"=SentinelNDVI20200304$X0,
                        "Landsat_NDVI_2019_05"=LandsatNDVI20190524$X0,
                        "Landsat_NDVI_2020_01"=LandsatNDVI20200103$X0,
                        "Landsat_NDVI_2020_03"=LandsatNDVI20200307$X0))
#Clean environment
rm(SentinelNDVI20190524,SentinelNDVI20200104,SentinelNDVI20200304,
   LandsatNDVI20190524,LandsatNDVI20200103,LandsatNDVI20200307)

# Standerdise
df_stand<-standardize(df[,5:40])
df_stand<-as.data.frame(cbind("Long_X"=df$Long_X,"Lat_Y"=df$Lat_Y,"Plantation"=df$Plantation,df_stand[,1:2],df_stand[,28:36]))
#Prepair data for looping models and create empty varibles
number=1
#Response
out_start=4
out_end=5
out_nvar=out_end-out_start+1
out_variable=rep(NA, out_nvar)
out_edf=rep(NA, out_nvar)
out_f_value=rep(NA, out_nvar)
out_p_value = rep(NA, out_nvar)
out_AIC=rep(NA, out_nvar)
#Explanatory
exp_start=6
exp_end=14
exp_nvar=exp_end-exp_start+1
exp_variable=rep(NA, exp_nvar)
exp_edf=rep(NA, exp_nvar)
exp_f_value=rep(NA, exp_nvar)
exp_p_value = rep(NA, out_nvar)
exp_AIC=rep(NA, exp_nvar)
#Loop models
for (i in out_start:out_end){
  outcome = colnames(df_stand)[i]
  for (j in exp_start:exp_end){
    exposure = colnames(df_stand)[j]
    
    model <-gamm4(get(outcome)~t2(get(exposure),bs="cc",k=4),random=~(1|Plantation),family=gaussian,na.action=na.exclude,data=df_stand)
    Summary<-as.numeric(summary(model$gam)[["s.table"]])
    AIC<-as.numeric(AIC(model$mer))
    
    out_edf[number] = Summary[1]
    out_f_value[number] = Summary[3]
    out_p_value[number] = Summary[4]
    out_AIC[number] = AIC[1]
    out_variable[number] = outcome
    number = number + 1
    
    exp_edf[number] = Summary[1]
    exp_f_value[number] = Summary[3]
    exp_p_value[number] = Summary[4]
    exp_AIC[number] = as.numeric(AIC[1])
    exp_variable[number] = exposure
    number = number + 1
  }
}
#Create a dataframe with results
outcome = data.frame(out_variable,out_edf,out_f_value,out_p_value,out_AIC)
exposure = data.frame(exp_variable,exp_edf,exp_f_value,exp_p_value,exp_AIC)
outcome = outcome %>% 
  rename(
    variable = out_variable,
    edf = out_edf,
    f_value = out_f_value,
    p_value = out_p_value,
    AIC = out_AIC
  )
exposure = exposure %>% 
  rename(
    variable = exp_variable,
    edf = exp_edf,
    f_value = exp_f_value,
    p_value = exp_p_value,
    AIC = exp_AIC
  )
exposure = na.omit(exposure)
outcome = na.omit(outcome)
write.csv(exposure,"R_code/Temp/exposure.csv",row.names = FALSE)
write.csv(outcome,"R_code/Temp/outcome.csv",row.names = FALSE)
exposure<- read.csv("R_code/Temp/exposure.csv",stringsAsFactors = F)
outcome<- read.csv("R_code/Temp/outcome.csv",stringsAsFactors = F)
colnames(exposure)<-c("Explanatory","edf","f_value","p_value","AIC")
all<- cbind("Response"=outcome$variable,exposure)
# Find significant results
Significant <- dplyr::filter(all, p_value <=0.05)
#Clean environment
rm(all,AIC,exp_AIC,exp_end,exp_f_value,exp_nvar,exp_p_value,exp_start,exp_variable,
   i,j,number,out_AIC,out_end,out_f_value,out_nvar,out_p_value,out_start,
   out_variable,Summary,model,exposure,outcome,exp_edf,out_edf)

#Model 1 validation and plotting
model1<-gamm4(Richness~t2(PlanetScope_NDVI_50_2019_05,bs="cc",k=4),random=~(1|Plantation),family=gaussian,na.action=na.exclude,data=df_stand)
residuals<-model1[["gam"]][["residuals"]]
fitted.values<-model1[["gam"]][["fitted.values"]]
plot(x=fitted.values,y=residuals,xlab="Fitted values",ylab="Residuals")
qqnorm(resid(model1$gam))
qqline(resid(model1$gam))
ggplot(data=df, aes(x =PlanetScope_NDVI_50_2019_05, y = Richness)) +
  geom_point() +
  geom_smooth(colour="black")+
  labs(y="Caelifera species richness\n",x="\nPlanetScope NDVI May 2019")+
  theme(axis.title=element_text(size=14),
        axis.text=element_text(size=13,colour="black"),
        axis.ticks=element_line(size=0.8,colour="black"),
        axis.line=element_line(size=0.8,colour="black"),
        panel.grid=element_blank(),
        panel.background=element_blank())
#Model 2 validation and plotting
model2<-gamm4(Richness~t2(PlanetScope_NDVI_50_2020_01,bs="cc",k=4),random=~(1|Plantation),family=gaussian,na.action=na.exclude,data=df_stand)
residuals<-model2[["gam"]][["residuals"]]
fitted.values<-model2[["gam"]][["fitted.values"]]
plot(x=fitted.values,y=residuals,xlab="Fitted values",ylab="Residuals")
qqnorm(resid(model2$gam))
qqline(resid(model2$gam))
ggplot(data=df, aes(x =PlanetScope_NDVI_50_2020_01, y = Richness)) +
  geom_point() +
  geom_smooth(colour="black")+
  labs(y="Caelifera species richness\n",x="\nPlanetScope NDVI January 2020")+
  theme(axis.title=element_text(size=14),
        axis.text=element_text(size=13,colour="black"),
        axis.ticks=element_line(size=0.8,colour="black"),
        axis.line=element_line(size=0.8,colour="black"),
        panel.grid=element_blank(),
        panel.background=element_blank())
#Model 3 validation and plotting
model3<-gamm4(Richness~t2(Sentinel_NDVI_2019_05,bs="cc",k=4),random=~(1|Plantation),family=gaussian,na.action=na.exclude,data=df_stand)
residuals<-model3[["gam"]][["residuals"]]
fitted.values<-model3[["gam"]][["fitted.values"]]
plot(x=fitted.values,y=residuals,xlab="Fitted values",ylab="Residuals")
qqnorm(resid(model3$gam))
qqline(resid(model3$gam))
ggplot(data=df, aes(x =Sentinel_NDVI_2019_05, y = Richness)) +
  geom_point() +
  geom_smooth(colour="black")+
  labs(y="Caelifera species richness\n",x="\nSentinel NDVI May 2019")+
  theme(axis.title=element_text(size=14),
        axis.text=element_text(size=13,colour="black"),
        axis.ticks=element_line(size=0.8,colour="black"),
        axis.line=element_line(size=0.8,colour="black"),
        panel.grid=element_blank(),
        panel.background=element_blank())
#Model 4 validation and plotting
model4<-gamm4(Richness~t2(Sentinel_NDVI_2020_01,bs="cc",k=4),random=~(1|Plantation),family=gaussian,na.action=na.exclude,data=df_stand)
residuals<-model4[["gam"]][["residuals"]]
fitted.values<-model4[["gam"]][["fitted.values"]]
plot(x=fitted.values,y=residuals,xlab="Fitted values",ylab="Residuals")
qqnorm(resid(model4$gam))
qqline(resid(model4$gam))
ggplot(data=df, aes(x =Sentinel_NDVI_2020_01, y = Richness)) +
  geom_point() +
  geom_smooth(colour="black")+
  labs(y="Caelifera species richness\n",x="\nSentinel NDVI January 2020")+
  theme(axis.title=element_text(size=14),
        axis.text=element_text(size=13,colour="black"),
        axis.ticks=element_line(size=0.8,colour="black"),
        axis.line=element_line(size=0.8,colour="black"),
        panel.grid=element_blank(),
        panel.background=element_blank())
#Model 5 validation and plotting
model5<-gamm4(Richness~t2(Sentinel_NDVI_2020_03,bs="cc",k=4),random=~(1|Plantation),family=gaussian,na.action=na.exclude,data=df_stand)
residuals<-model5[["gam"]][["residuals"]]
fitted.values<-model5[["gam"]][["fitted.values"]]
plot(x=fitted.values,y=residuals,xlab="Fitted values",ylab="Residuals")
qqnorm(resid(model5$gam))
qqline(resid(model5$gam))
ggplot(data=df, aes(x =Sentinel_NDVI_2020_03, y = Richness)) +
  geom_point() +
  geom_smooth(colour="black")+
  labs(y="Caelifera species richness\n",x="\nSentinel NDVI March 2020")+
  theme(axis.title=element_text(size=14),
        axis.text=element_text(size=13,colour="black"),
        axis.ticks=element_line(size=0.8,colour="black"),
        axis.line=element_line(size=0.8,colour="black"),
        panel.grid=element_blank(),
        panel.background=element_blank())
#Model 6 validation and plotting
model6<-gamm4(exShannon~t2(Landsat_NDVI_2020_03,bs="cc",k=4),random=~(1|Plantation),family=gaussian,na.action=na.exclude,data=df_stand)
residuals<-model6[["gam"]][["residuals"]]
fitted.values<-model6[["gam"]][["fitted.values"]]
plot(x=fitted.values,y=residuals,xlab="Fitted values",ylab="Residuals")
qqnorm(resid(model6$gam))
qqline(resid(model6$gam))
ggplot(data=df, aes(x =Landsat_NDVI_2020_03, y = exShannon)) +
  geom_point() +
  geom_smooth(colour="black")+
  labs(y="Caelifera exponent Shannon diversity\n",x="\nLandsat NDVI March 2020")+
  theme(axis.title=element_text(size=14),
        axis.text=element_text(size=13,colour="black"),
        axis.ticks=element_line(size=0.8,colour="black"),
        axis.line=element_line(size=0.8,colour="black"),
        panel.grid=element_blank(),
        panel.background=element_blank())
#Clean environment
rm(model1,model2,model3,model4,model5,model6,Significant,fitted.values,residuals)

########### Local models ###########
#Load vegetation data
veg<-read.csv("Excel_sheets/Vegetation_Plot_Data.csv")
#Calculate average vegetation height
veg[veg==0]<-NA
veg_ave<-veg%>%dplyr::select(ends_with(c("Hei")))
veg_ave$Shrub_Hei<-rowMeans(veg_ave[,c("Shrub_Hei","Bram_Hei")],na.rm=TRUE)
veg_ave$Bram_Hei=NULL
veg_ave[is.na(veg_ave)]<-NA
veg_ave<-as.data.frame(cbind(veg%>%dplyr::select(starts_with(c("Site_ID"))),veg_ave))
veg_ave<-veg_ave %>%
  pivot_longer(Tree_Hei:Climb_Hei) %>%
  group_by(Site_ID) %>%
  summarise(mean = mean(value,na.rm = T))
df<-as.data.frame(cbind(df,"AveVegHei"=veg_ave$mean))
#Extract average vegetation richness per site
veg_Rich<-veg%>%dplyr::select(ends_with(c("Rich")))
veg_Rich$Shrub_Rich<-rowMeans(veg_Rich[,c("Shrub_Rich","Bram_Rich")],na.rm=TRUE)
veg_Rich$Bram_Rich=NULL
veg_Rich[is.na(veg_Rich)]<-NA
veg_Rich<-as.data.frame(cbind(veg%>%dplyr::select(starts_with(c("Site_ID"))),veg_Rich))
veg_Rich<-veg_Rich %>% 
  group_by(Site_ID) %>% 
  summarise_each(funs(mean(.,na.rm=TRUE)))
veg_Rich[is.na(veg_Rich)]<-NA
veg_Rich$Site_ID=NULL
veg_Rich<-as.data.frame(rowSums(veg_Rich,na.rm=T))
df<-as.data.frame(cbind(df,"veg_Mean_Rich"=veg_Rich$`rowSums(veg_Rich, na.rm = T)`))
#Extract growth form richness
veg_Func_Rich<-veg%>%dplyr::select(ends_with(c("Rich")))
veg_Func_Rich$Bram_Rich=NULL
veg_Func_Rich[is.na(veg_Func_Rich)]<-NA
veg_Func_Rich<-as.data.frame(cbind(veg%>%dplyr::select(starts_with(c("Site_ID"))),veg_Func_Rich))
veg_Func_Rich<-veg_Func_Rich %>% 
  group_by(Site_ID) %>% 
  summarise_each(funs(mean(.,na.rm=TRUE)))
veg_Func_Rich$Site_ID=NULL
veg_Func_Rich[is.na(veg_Func_Rich)]<-0
veg_Func_Rich<-as.data.frame(rowSums(veg_Func_Rich != 0))
df<-as.data.frame(cbind(df,"veg_Func_Rich"=veg_Func_Rich$`rowSums(veg_Func_Rich != 0)`))
#Extract bramble abundance and cover, ground cover, and Rockiness
veg[is.na(veg)]<-0
veg$Bram_Abu[veg$Bram_Abu==0]<-NA
features<-veg[,7:48]
features<-as.data.frame(cbind(veg%>%dplyr::select(c("Site_ID")),features))
features<-features %>% 
  group_by(Site_ID) %>% 
  summarise_each(funs(mean(.,na.rm=TRUE)))
features[is.na(features)]<-0
df<-as.data.frame(cbind(df,"Bram_Abu"=features$Bram_Abu,"Bram_Cov"=features$Bram_Cov,
                        "Rock_Cov"=features$Rock_Cov,"Grou_Cov"=features$Grou_Cov))
#Clean environment
rm(veg_ave,veg_Rich,veg_Func_Rich,features,veg)
#Check for correlation between features
features<-df[,41:47]
chart.Correlation(as.matrix(features),TRUE)
df$veg_Func_Rich=NULL
df$Bram_Cov=NULL
#Clean environment
rm(features)
#Standerdise
df_stand<-standardize(df[,5:45])
df_stand<-as.data.frame(cbind("Site_ID"=df$Site_ID,"Long_X"=df$Long_X,"Lat_Y"=df$Lat_Y,"Plantation"=df$Plantation,df_stand))

#Modeling and plotting Richness
Model1<-lmer(Richness~AveVegHei+veg_Mean_Rich+Bram_Abu+Grou_Cov+Rock_Cov+(1|Plantation),data=df_stand)
vif(Model1)
options(na.action = "na.fail")
Model1_Dredge<-dredge(Model1,evaluate=TRUE,rank=AICc)
options(na.action = "na.omit")
Model1_Subset<-subset(Model1_Dredge,delta<2)
Model1_Ave<-model.avg(Model1_Subset)
# Reporting
importance(Model1_Ave)
confint(Model1_Ave)
summary(Model1_Ave)
rm(Model1_Dredge,Model1_Subset,Model1,Model1_Ave)
ggplot(data=df,aes(x=AveVegHei,y=Richness))+
  geom_point()+
  geom_smooth(method="lm",colour="black")+
  labs(y="Caelifera species richness\n",x="\nMean vegetation height")+
  theme(axis.title=element_text(size=14),
        axis.text=element_text(size=13,colour="black"),
        axis.ticks=element_Estate2e(size=0.8,colour="black"),
        axis.Estate2e=element_Estate2e(size=0.8,colour="black"),
        panel.grid=element_blank(),
        panel.background=element_blank())
#Modeling and plotting diversity
Model2<-lmer(exShannon~AveVegHei+veg_Mean_Rich+Bram_Abu+Grou_Cov+Rock_Cov+(1|Plantation),data=df_stand)
vif(Model2)
options(na.action = "na.fail")
Model2_Dredge<-dredge(Model2,evaluate=TRUE,rank=AICc)
options(na.action = "na.omit")
Model2_Subset<-subset(Model2_Dredge,delta<2)
Model2_Ave<-model.avg(Model2_Subset)
#Reporting
importance(Model2_Ave)
confint(Model2_Ave)
summary(Model2_Ave)
rm(Model2_Dredge,Model2_Subset,Model2,Model2_Ave)
ggplot(data=df,aes(x=Grou_Cov,y=exShannon))+
  geom_point()+
  geom_smooth(method="lm",colour="black")+
  labs(y="Caelifera exponent Shannon diversity\n",x="\nMean ground cover")+
  theme(axis.title=element_text(size=14),
        axis.text=element_text(size=13,colour="black"),
        axis.ticks=element_Estate2e(size=0.8,colour="black"),
        axis.Estate2e=element_Estate2e(size=0.8,colour="black"),
        panel.grid=element_blank(),
        panel.background=element_blank())

########### Lag models MVAbund ###########
#Create mvabund object
mvhopper<-mvabund(hopper)
#Standerdise
df_stand<-standardize(df[,7:31])
df_stand<-cbind(df[1:4],df_stand)
#Prepair data for looping models and create empty varibles
number=1
#Explanatory
exp_start=5
exp_end=29
exp_nvar=exp_end-exp_start+1
exp_variable=rep(NA, exp_nvar)
exp_score=rep(NA, exp_nvar)
exp_pvalue=rep(NA, exp_nvar)
exp_AIC=rep(NA, exp_nvar)
#Loop
for (j in exp_start:exp_end){
  exposure = colnames(df_stand)[j]
  
  model<-manyglm(mvhopper~get(exposure)+Plantation,family="negative binomial",data=df_stand)
  anova<-anova(model,resamp="perm.resid",test="score",cor.type="I")
  
  Multi_stat<-as.data.frame(anova[["table"]])
  Multi_stat<-Multi_stat[-c(1,3),]
  score<-Multi_stat$score
  Pr<-Multi_stat$`Pr(>score)`
  AIC<-model[["AICsum"]]
  exp_score[number] = as.numeric(score)
  exp_pvalue[number] = as.numeric(Pr)
  exp_AIC[number] = as.numeric(AIC)
  exp_variable[number] = exposure
  number = number + 1
}
#Create a dataframe with results
exposure = data.frame(exp_variable,exp_score,exp_pvalue,exp_AIC)
exposure = exposure %>% 
  rename(
    variable = exp_variable,
    score = exp_score,
    p_value = exp_pvalue,
    AIC = exp_AIC
  )
exposure = na.omit(exposure)
# Find significant results
Significant <- dplyr::filter(exposure, p_value <=0.05)
#Clean environment
rm(exposure,model,AIC,exp_AIC,exp_end,exp_score,exp_nvar,Pr,score,exp_pvalue,exp_start,exp_variable,j,number,Multi_stat,anova)

#Perform fourth corner analysis
Family<-read.csv("Excel_sheets/Caelifera_Families.csv",row.names = 1)
#Creating individual variables
NDVI_sdY19_20<-as.data.frame(cbind("Plantation"=df_stand$Plantation,"NDVI 2019-2020"=df_stand$NDVI_sdY19_20))#4x5
NDVI_sdY18_20<-as.data.frame(cbind("Plantation"=df_stand$Plantation,"NDVI 2018-2020"=df_stand$NDVI_sdY18_20))
#Running 4th corner
Corner_NDVI_sdY19_20<-traitglm(hopper,NDVI_sdY19_20,Family,
                      family="negative.binomial",
                      method="glm1path")
NDVI_sdY19_20<-data.matrix(Corner_NDVI_sdY19_20$fourth) 
Corner_NDVI_sdY18_20<-traitglm(hopper,NDVI_sdY18_20,Family,
                      family="negative.binomial",
                      method="glm1path")
NDVI_sdY18_20<-data.matrix(Corner_NDVI_sdY18_20$fourth) 
#Creating Coefficient matrix for plotting later
Coefficients<-as.data.frame(cbind(NDVI_sdY19_20,NDVI_sdY18_20))
write.csv(Coefficients,"R_code/Temp/LagCornerCoefficients.csv",row.names = TRUE)
#Clean environment
rm(Coefficients,Significant,NDVI_sdY19_20,NDVI_sdY18_20,Corner_NDVI_sdY19_20,Corner_NDVI_sdY18_20)

########### Spectral models MVAbund ###########
#Standerdise
df_stand<-standardize(df[,32:40])
df_stand<-cbind(df[1:4],df_stand)
#Prepair data for looping models and create empty varibles
number=1
#Explanatory
exp_start=5
exp_end=13
exp_nvar=exp_end-exp_start+1
exp_variable=rep(NA, exp_nvar)
exp_score=rep(NA, exp_nvar)
exp_pvalue=rep(NA, exp_nvar)
exp_AIC=rep(NA, exp_nvar)
#Loop
for (j in exp_start:exp_end){
  exposure = colnames(df_stand)[j]
  
  model<-manyglm(mvhopper~get(exposure)+Plantation,family="negative binomial",data=df_stand)
  anova<-anova(model,resamp="perm.resid",test="score",cor.type="I")
  
  Multi_stat<-as.data.frame(anova[["table"]])
  Multi_stat<-Multi_stat[-c(1,3),]
  score<-Multi_stat$score
  Pr<-Multi_stat$`Pr(>score)`
  AIC<-model[["AICsum"]]
  exp_score[number] = as.numeric(score)
  exp_pvalue[number] = as.numeric(Pr)
  exp_AIC[number] = as.numeric(AIC)
  exp_variable[number] = exposure
  number = number + 1
}
#Create a dataframe with results
exposure = data.frame(exp_variable,exp_score,exp_pvalue,exp_AIC)
exposure = exposure %>% 
  rename(
    variable = exp_variable,
    score = exp_score,
    p_value = exp_pvalue,
    AIC = exp_AIC
  )
exposure = na.omit(exposure)
# Find significant results
Significant <- dplyr::filter(exposure, p_value <=0.05)
#Clean environment
rm(exposure,model,AIC,exp_AIC,exp_end,exp_score,exp_nvar,Pr,score,exp_pvalue,exp_start,exp_variable,j,number,Multi_stat,anova)

#Perform fourth corner analysis
#Creating individual variables
PlanetScope_NDVI_50_2019_05<-as.data.frame(cbind("Plantation"=df_stand$Plantation,"PlanetScope NDVI May 2019"=df_stand$PlanetScope_NDVI_50_2019_05))
Landsat_NDVI_2019_05<-as.data.frame(cbind("Plantation"=df_stand$Plantation,"Landsat NDVI May 2019"=df_stand$Landsat_NDVI_2019_05))
#Running 4th corner
Corner_PlanetScope_NDVI_50_2019_05<-traitglm(hopper,PlanetScope_NDVI_50_2019_05,Family,
                               family="negative.binomial",
                               method="glm1path")
PlanetScope_NDVI_50_2019_05<-data.matrix(Corner_PlanetScope_NDVI_50_2019_05$fourth) 
Corner_Landsat_NDVI_2019_05<-traitglm(hopper,Landsat_NDVI_2019_05,Family,
                               family="negative.binomial",
                               method="glm1path")
Landsat_NDVI_2019_05<-data.matrix(Corner_Landsat_NDVI_2019_05$fourth) 
#Creating Coefficient matrix for plotting later
Coefficients<-as.data.frame(cbind(PlanetScope_NDVI_50_2019_05,Landsat_NDVI_2019_05))
write.csv(Coefficients,"R_code/Temp/SpectralCornerCoefficients.csv",row.names = TRUE)
#Clean environment
rm(Coefficients,Significant,PlanetScope_NDVI_50_2019_05,Landsat_NDVI_2019_05,Corner_PlanetScope_NDVI_50_2019_05,Corner_Landsat_NDVI_2019_05)

########### Local models MVAbund ###########
#Standerdise
df_stand<-standardize(df[,41:45])
df_stand<-cbind(df[1:4],df_stand)
# Creating model
model<-manyglm(mvhopper~Plantation+AveVegHei+veg_Mean_Rich+Bram_Abu+Rock_Cov+Grou_Cov,
               family="negative binomial",data=df_stand)
anova(model,resamp="perm.resid",test="score",cor.type="I")
model[["AICsum"]]
#Running 4th corner
Local_model<-as.data.frame(cbind("Plantation"=df_stand$Plantation,
                                 "Mean Veg Height"=df_stand$AveVegHei,
                                 "Mean Veg Richness"=df_stand$veg_Mean_Rich,
                                 "Mean Bramble Abundance"=df_stand$Bram_Abu,
                                 "Mean Rock Cover"=df_stand$Rock_Cov,
                                 "Mean Ground Cover"=df_stand$Grou_Cov))
corner_Local_model<-traitglm(hopper,Local_model,Family,
                             family="negative.binomial",
                             method="glm1path")
#Creating Coefficient matrix for plotting later 
Coefficients<-data.matrix(corner_Local_model$fourth)
write.csv(Coefficients,"R_code/Temp/LocalCornerCoefficients.csv",row.names = TRUE)
#Clean environment
rm(Family,hopper,mvhopper,model,corner_Local_model,Coefficients,Local_model,df_stand)

#### Plot regression results together ####
#Figure 2
a<-ggplot(data=df,aes(x=NDVI_sdY19_20,y=Richness))+
  geom_point()+
  geom_smooth(method="lm",colour="black")+
  labs(y="Caelifera species richness\n",x="\nDeviation of NDVI between 2019-2020")+
  theme(axis.title=element_text(size=14),
        axis.text=element_text(size=13,colour="black"),
        axis.ticks=element_line(size=0.8,colour="black"),
        axis.line=element_line(size=0.8,colour="black"),
        panel.grid=element_blank(),
        panel.background=element_blank())
b<-ggplot(data=df, aes(x =NDVI_sdY19_20, y = exShannon)) +
  geom_point() +
  geom_smooth(method="lm",colour="black")+
  labs(y="Caelifera exponent Shannon diversity\n",x="\nDeviation of NDVI between 2019-2020")+
  theme(axis.title=element_text(size=14),
        axis.text=element_text(size=13,colour="black"),
        axis.ticks=element_line(size=0.8,colour="black"),
        axis.line=element_line(size=0.8,colour="black"),
        panel.grid=element_blank(),
        panel.background=element_blank())
c<-ggplot(data=df, aes(x =NDVI_sdY17_20, y = exShannon)) +
  geom_point() +
  geom_smooth(method="lm",colour="black")+
  labs(y="Caelifera exponent Shannon diversity\n",x="\nDeviation of NDVI between 2017-2020")+
  theme(axis.title=element_text(size=14),
        axis.text=element_text(size=13,colour="black"),
        axis.ticks=element_line(size=0.8,colour="black"),
        axis.line=element_line(size=0.8,colour="black"),
        panel.grid=element_blank(),
        panel.background=element_blank())
d<-ggplot(data=df,aes(x=Grou_Cov,y=exShannon))+
  geom_point()+
  geom_smooth(method="lm",colour="black")+
  labs(y="Caelifera exponent Shannon diversity\n",x="\nMean ground cover")+
  theme(axis.title=element_text(size=14),
        axis.text=element_text(size=13,colour="black"),
        axis.ticks=element_line(size=0.8,colour="black"),
        axis.line=element_line(size=0.8,colour="black"),
        panel.grid=element_blank(),
        panel.background=element_blank())
e<-ggplot(data=df,aes(x=AveVegHei,y=Richness))+
  geom_point()+
  geom_smooth(method="lm",colour="black")+
  labs(y="Caelifera species richness\n",x="\nMean vegetation height")+
  theme(axis.title=element_text(size=14),
        axis.text=element_text(size=13,colour="black"),
        axis.ticks=element_line(size=0.8,colour="black"),
        axis.line=element_line(size=0.8,colour="black"),
        panel.grid=element_blank(),
        panel.background=element_blank())
ggarrange(a,b,c,d,e,ncol=3,nrow=2)

#Figure 4
a<-ggplot(data=df, aes(x =PlanetScope_NDVI_50_2019_05, y = Richness)) +
  geom_point() +
  geom_smooth(colour="black")+
  labs(y="Caelifera species richness\n",x="\nPlanetScope NDVI May 2019")+
  theme(axis.title=element_text(size=14),
        axis.text=element_text(size=13,colour="black"),
        axis.ticks=element_line(size=0.8,colour="black"),
        axis.line=element_line(size=0.8,colour="black"),
        panel.grid=element_blank(),
        panel.background=element_blank())
b<-ggplot(data=df, aes(x =PlanetScope_NDVI_50_2020_01, y = Richness)) +
  geom_point() +
  geom_smooth(colour="black")+
  labs(y="Caelifera species richness\n",x="\nPlanetScope NDVI January 2020")+
  theme(axis.title=element_text(size=14),
        axis.text=element_text(size=13,colour="black"),
        axis.ticks=element_line(size=0.8,colour="black"),
        axis.line=element_line(size=0.8,colour="black"),
        panel.grid=element_blank(),
        panel.background=element_blank())
c<-ggplot(data=df, aes(x =Sentinel_NDVI_2019_05, y = Richness)) +
  geom_point() +
  geom_smooth(colour="black")+
  labs(y="Caelifera species richness\n",x="\nSentinel NDVI May 2019")+
  theme(axis.title=element_text(size=14),
        axis.text=element_text(size=13,colour="black"),
        axis.ticks=element_line(size=0.8,colour="black"),
        axis.line=element_line(size=0.8,colour="black"),
        panel.grid=element_blank(),
        panel.background=element_blank())
d<-ggplot(data=df, aes(x =Sentinel_NDVI_2020_01, y = Richness)) +
  geom_point() +
  geom_smooth(colour="black")+
  labs(y="Caelifera species richness\n",x="\nSentinel NDVI January 2020")+
  theme(axis.title=element_text(size=14),
        axis.text=element_text(size=13,colour="black"),
        axis.ticks=element_line(size=0.8,colour="black"),
        axis.line=element_line(size=0.8,colour="black"),
        panel.grid=element_blank(),
        panel.background=element_blank())
e<-ggplot(data=df, aes(x =Sentinel_NDVI_2020_03, y = Richness)) +
  geom_point() +
  geom_smooth(colour="black")+
  labs(y="Caelifera species richness\n",x="\nSentinel NDVI March 2020")+
  theme(axis.title=element_text(size=14),
        axis.text=element_text(size=13,colour="black"),
        axis.ticks=element_line(size=0.8,colour="black"),
        axis.line=element_line(size=0.8,colour="black"),
        panel.grid=element_blank(),
        panel.background=element_blank())
f<-ggplot(data=df, aes(x =Landsat_NDVI_2020_03, y = exShannon)) +
  geom_point() +
  geom_smooth(colour="black")+
  labs(y="Caelifera exponent Shannon diversity\n",x="\nLandsat NDVI March 2020")+
  theme(axis.title=element_text(size=14),
        axis.text=element_text(size=13,colour="black"),
        axis.ticks=element_line(size=0.8,colour="black"),
        axis.line=element_line(size=0.8,colour="black"),
        panel.grid=element_blank(),
        panel.background=element_blank())
ggarrange(a,b,c,d,e,f,ncol=3,nrow=2)
#Clean environment
rm(a,b,c,d,e,f,df)

#### Plot 4th corner analysis together ####
#Load model coefficients
LagCoefficients<-read.csv("R_code/Temp/LagCornerCoefficients.csv",row.names = 1)
LagCoefficients<-LagCoefficients[-c(1,3)]
RowNames<-rownames(LagCoefficients)
RowNames<-sub("1","",RowNames)
rownames(LagCoefficients)<-RowNames
SpectralCoefficients<-read.csv("R_code/Temp/SpectralCornerCoefficients.csv",row.names = 1)
SpectralCoefficients<-SpectralCoefficients[-c(1,3)]
RowNames<-rownames(SpectralCoefficients)
RowNames<-sub("1","",RowNames)
rownames(SpectralCoefficients)<-RowNames
LocalCoefficients<-read.csv("R_code/Temp/LocalCornerCoefficients.csv",row.names = 1)
LocalCoefficients<-LocalCoefficients[-c(1)]
RowNames<-rownames(LocalCoefficients)
RowNames<-sub("1","",RowNames)
rownames(LocalCoefficients)<-RowNames
#Merge into matrix
Coefficients<-as.data.frame(cbind(LagCoefficients,SpectralCoefficients,LocalCoefficients))
#Plot
a = max(abs(Coefficients))
colort = colorRampPalette(c("blue","white","red"))
pdf("R_code/Temp/pdf/Figure5.pdf")#h6.w5x5
Figure5<-levelplot(t(as.matrix(Coefficients)), 
                   xlab=list("\nVariables",cex=1),
                   ylab=list("Sub families",cex=1),
                   col.regions=colort(100),
                   at=seq(-a,a,length=100),
                   scales=list(x=list(rot=90),cex=0.6))
dev.off() 
#Edit PDF
setwd("R_code/Temp/pdf/")
Figure5<-image_read_pdf(path="Figure5.pdf",density=900) %>%
  image_annotate("Lag",size=100,location="+2390+5800",color="black") %>%
  image_annotate("Spectral",size=100,location="+2775+5800",color="black") %>%
  image_annotate("Local",size=100,location="+3680+5800",color="black")
Figure5_Edit<-image_draw(Figure5)#Start X,Start Y,End X,End Y
rect(2278,5770,2680,5770, border="black",lty="solid",lwd=8)
rect(2278,5680,2278,5770, border="black",lty="solid",lwd=8)
rect(2680,5680,2680,5770, border="black",lty="solid",lwd=8)
rect(2755,5770,3157,5770, border="black",lty="solid",lwd=8)
rect(2755,5680,2755,5770, border="black",lty="solid",lwd=8)
rect(3157,5680,3157,5770, border="black",lty="solid",lwd=8)
rect(3232,5770,4353,5770, border="black",lty="solid",lwd=8)
rect(3232,5680,3232,5770, border="black",lty="solid",lwd=8)
rect(4353,5680,4353,5770, border="black",lty="solid",lwd=8)
dev.off()
#Save Figure5
image_write(Figure5_Edit,path="D:/Jurie/Documents/Postgrad/PhD/0_Writing/Publication_Data_Chapter2/Figure5.pdf",format="pdf")
#Clean environment
rm(Figure5,Figure5_Edit,Coefficients,SpectralCoefficients,LocalCoefficients,RowNames,colort,LagCoefficients)
