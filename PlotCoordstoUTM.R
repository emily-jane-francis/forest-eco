# Plot Coords to UTM

library(dplyr); library(tidyr); library(Metrics)

mainfolder = "C:/Users/ejf834/Box/ForestGeoTempTrop/NewCodes/"
source(paste(mainfolder,"PowerLawFit_etc.r",sep=""))
source(paste(mainfolder,"BCIdataAnalysisFunctions.r",sep="")) 
source(paste(mainfolder,"ForestSimulationFunctions.r",sep=""))
source("C:/Users/ejf834/Box/ForestGeoTempTrop/NewCodes/find_dstar_onetreecentric.r")

plotcode = "wdrv"

if(plotcode=="wdrv"){ncensuses=2; sizem2= 256000; phi=pi*.5^2/10^.9; theta = 0.9 } #silver fir allometry from Ameztegui et al. 2012 (need to update with info from Sillet)

datafolder = paste("C:/Users/ejf834/Box/ForestGeoTempTrop/data/",plotcode,"/",sep="")

calcfolder = paste(datafolder,"ejfdataoutput/",sep="")
simfolder = paste(datafolder,"ejfsimoutput/",sep="")

datallwUTM = read.table(paste(calcfolder,plotcode,"allwUTM.txt",sep=""),sep="\t",header=TRUE)
datallwUTM$tag = levels(datallwUTM$tag)[datallwUTM$tag]

#### estimate with out of order coordinates #### 

datallwUTMjumbled = datallwUTM
set.seed(42)
jumbledrows = sample(nrow(datallwUTMjumbled))

# name variables #
aj= c()
bj = c()
cj = c()
dj = c()
ej = c()
fj = c()
tag1 = c()
tag2 = c()
tag3 = c()

for(i in seq(1,dim(datallwUTMjumbled)[1]/3)){
  # solve for a, b, and c
  tryCatch({
    A <- matrix(data=c(datallwUTMjumbled$gx[i],datallwUTMjumbled$gy[i],1,
                       datallwUTMjumbled$gx[i+1],datallwUTMjumbled$gy[i+1],1,
                       datallwUTMjumbled$gx[i+2],datallwUTMjumbled$gy[i+2],1),nrow=3,ncol=3,byrow=TRUE)
    
    sol <- matrix(data=c(datallwUTMjumbled$utmx[i],datallwUTMjumbled$utmx[i+1],datallwUTMjumbled$utmx[i+2]),nrow=3,ncol=1,byrow=FALSE)
    solve(A)%*% sol
    aj[i] = (solve(A)%*% sol)[1,1]
    bj[i] = (solve(A)%*% sol)[2,1]
    cj[i] = (solve(A)%*% sol)[3,1]
    # gy 
    
    B <- matrix(data=c(datallwUTMjumbled$gx[i],datallwUTMjumbled$gy[i],1,
                       datallwUTMjumbled$gx[i+1],datallwUTMjumbled$gy[i+1],1,
                       datallwUTMjumbled$gx[i+2],datallwUTMjumbled$gy[i+2],1),nrow=3,ncol=3,byrow=TRUE)
    
    soly <- matrix(data=c(datallwUTMjumbled$utmy[i],datallwUTMjumbled$utmy[i+1],datallwUTMjumbled$utmy[i+2]),nrow=3,ncol=1,byrow=FALSE)
    solve(B)%*% soly
    dj[i] = (solve(B)%*% soly)[1,1]
    ej[i] = (solve(B)%*% soly)[2,1]
    fj[i] = (solve(B)%*% soly)[3,1]
    
    # store trees used for the fit
    tag1[i] = datallwUTMjumbled$tag[i]
    tag2[i] = datallwUTMjumbled$tag[i+1]
    tag3[i] = datallwUTMjumbled$tag[i+2]
    
    #
  },error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
}

rotationj = data.frame(aj,bj,cj,dj,ej,fj,tag1,tag2,tag3)

#
rotationj$tag1 = levels(rotationj$tag1)[rotationj$tag1]
rotationj$tag2 = levels(rotationj$tag2)[rotationj$tag2]
rotationj$tag3 = levels(rotationj$tag3)[rotationj$tag3]

#
rotationj = rotationj[!is.na(rotationj$a),]
#

par(mfrow = c(3,2),mar = c(4.1,4.1,4.1,2.1))
for(j in 1:dim(rotationj)[2]){
  hist(rotationj[,j],main = colnames(rotationj)[j],
       ylab = "Frequency",
       col = "light blue",,xlab = "",ylim=c(0,10000))
  sd(rotationj[,j],na.rm=T)
  abline(v = min(rotationj[,j],na.rm=T))
  abline(v=max(rotationj[,j],na.rm=T))
  abline(v=mean(rotationj[,j],na.rm=T),col = "blue")
  abline(v=median(rotationj[,j],na.rm=T),col="red")
}


### check accuracy ### 

# Test the estimates #
# a row for each tree with UTM coordinates, that was not included in model fitting #
utmx_error_est = matrix(nrow = length(datallwUTMjumbled$gx),ncol = length(rotationj$a)+1)
utmy_error_est = matrix(nrow = length(datallwUTMjumbled$gy),ncol = length(rotationj$d)+1)

utmx_error_est = data.frame(utmx_error_est)
utmy_error_est = data.frame(utmy_error_est)

#
'%notin%' = Negate('%in%')
tagsout = list()
for(i in 1:length(rotationj$aj)){
  # need to change this to exclude the data used to build the transformation #
  # tagsout = rotationj[i,7:9] #
  # databin = datallwUTMjumbled[datallwUTMjumbled$tag %notin% tagsout,] #
  tagsout[[i]] = c(rotationj[i,7],rotationj[i,8],rotationj[i,9])
  utmx_error_est[,i] = rotationj$aj[i]*datallwUTMjumbled$gx + rotationj$bj[i]*datallwUTMjumbled$gy + rotationj$cj[i]
  utmy_error_est[,i] = rotationj$dj[i]*datallwUTMjumbled$gx + rotationj$ej[i]*datallwUTMjumbled$gy + rotationj$fj[i]
}



# make the first column the actuals 
utmx_error_est = cbind(datallwUTMjumbled$tag,datallwUTMjumbled$utmx,utmx_error_est)
utmy_error_est = cbind(datallwUTMjumbled$tag,datallwUTMjumbled$utmy,utmy_error_est)

# 
colnames(utmx_error_est)[1:2] = c("tag","actual_utmx")
colnames(utmy_error_est)[1:2] = c("tag","actual_utmy")

# remove NA's (produced by rows of datallwUTMjumbled that had no utm x and utm y to calculate from) #
utmx_error_est = utmx_error_est[!is.na(utmx_error_est$X1),]
utmy_error_est = utmy_error_est[!is.na(utmy_error_est$X1),]

utmx_error_est = utmx_error_est[!is.na(utmx_error_est$actual_utmx),]
utmy_error_est = utmy_error_est[!is.na(utmy_error_est$actual_utmy),]


# Tags to Character Variable # 
utmx_error_est$tag = levels(utmx_error_est$tag)[utmx_error_est$tag]
utmy_error_est$tag = levels(utmy_error_est$tag)[utmy_error_est$tag]


# Visual check on datasets #
View(utmx_error_est)
dim(utmx_error_est)
View(utmy_error_est)
dim(utmy_error_est)

#
rmse_utmx = c()
rmse_utmy = c()
rmse_utmsum = c()

for(i in 1:ncol(utmx_error_est)){
  
  # remove tags of trees that were used to fit the model #
  databinx = utmx_error_est[utmx_error_est$tag %notin% tagsout[[i]],]
  databiny = utmy_error_est[utmy_error_est$tag %notin% tagsout[[i]],]
  
  rmse_utmx[i] = rmse(databinx$actual_utmx,databinx[,i+2])
  rmse_utmy[i] = rmse(databiny$actual_utmy,databiny[,i+2])
  rmse_utmsum[i] = sum(rmse_utmx[i],rmse_utmy[i])
}

#
rmse_dat = data.frame(rmse_utmx,rmse_utmy,rmse_utmsum)
rmse_dat = rmse_dat[!is.na(rmse_dat$rmse_utmx),]
View(rmse_dat)
View(rmse_dat[rmse_dat$rmse_utmsum == min(rmse_dat$rmse_utmsum),])
save(rmse_dat,file = paste(calcfolder,"rmse_dat_UTMtoPlotCoords_March29.RData",sep=""))
save(utmx_error_est,file = paste(calcfolder,"utmx_error_est_UTMtoPlotCoords_March29.RData",sep=""))

bestrow = rownames(rmse_dat[rmse_dat$rmse_utmsum == min(rmse_dat$rmse_utmsum),])

# 

# coefs #
best_coefs_utm = rotationj[as.numeric(bestrow),1:6]
# check the best coefs #
best_coefs_utm$aj

utmx_test = best_coefs_utm$aj*datallwUTMjumbled$gx + best_coefs_utm$bj*datallwUTMjumbled$gy + best_coefs_utm$cj
utmy_test = best_coefs_utm$dj*datallwUTMjumbled$gx + best_coefs_utm$ej*datallwUTMjumbled$gy + best_coefs_utm$fj

# check if the best estimates calculated from the best coefs line up #
# with the corresponding row in the gx/gy error est dataframe # 
View(utmx_error_est[,(as.numeric(bestrow)+2)])
View(utmx_test)
View(utmy_error_est[,(as.numeric(bestrow)+2)])
View(utmy_test)

#
save(best_coefs_utm,file = paste(calcfolder,"best_coefs_UTMtoPlotCoords_March29.RData",sep=""))
load(paste(calcfolder,"best_coefs_UTMtoplotcoords_March29.RData",sep=""))

### Take in CAindex and Dstar Data # 
pradius = 15; div = 4
file=paste(calcfolder,pradius,"_",div,".txt",sep="")
dstardata = read.table(file,sep="\t",header=TRUE)

# 
dstardata$utmx = best_coefs_utm$aj*dstardata$x + best_coefs_utm$bj*dstardata$y + best_coefs_utm$cj
dstardata$utmy = best_coefs_utm$dj*dstardata$x + best_coefs_utm$ej*dstardata$y + best_coefs_utm$fj

#
write.csv(dstardata,file = paste(calcfolder,"dstardata_geo.csv"))


