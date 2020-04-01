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
    A <- matrix(data=c(datallwUTMjumbled$utmx[i],datallwUTMjumbled$utmy[i],1,
                       datallwUTMjumbled$utmx[i+1],datallwUTMjumbled$utmy[i+1],1,
                       datallwUTMjumbled$utmx[i+2],datallwUTMjumbled$utmy[i+2],1),nrow=3,ncol=3,byrow=TRUE)
    
    sol <- matrix(data=c(datallwUTMjumbled$gx[i],datallwUTMjumbled$gx[i+1],datallwUTMjumbled$gx[i+2]),nrow=3,ncol=1,byrow=FALSE)
    solve(A)%*% sol
    aj[i] = (solve(A)%*% sol)[1,1]
    bj[i] = (solve(A)%*% sol)[2,1]
    cj[i] = (solve(A)%*% sol)[3,1]
    # gy 
    
    B <- matrix(data=c(datallwUTMjumbled$utmx[i],datallwUTMjumbled$utmy[i],1,
                       datallwUTMjumbled$utmx[i+1],datallwUTMjumbled$utmy[i+1],1,
                       datallwUTMjumbled$utmx[i+2],datallwUTMjumbled$utmy[i+2],1),nrow=3,ncol=3,byrow=TRUE)
    
    soly <- matrix(data=c(datallwUTMjumbled$gy[i],datallwUTMjumbled$gy[i+1],datallwUTMjumbled$gy[i+2]),nrow=3,ncol=1,byrow=FALSE)
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
gx_error_est = matrix(nrow = length(datallwUTMjumbled$utmx),ncol = length(rotationj$a)+1)
gy_error_est = matrix(nrow = length(datallwUTMjumbled$utmy),ncol = length(rotationj$d)+1)

gx_error_est = data.frame(gx_error_est)
gy_error_est = data.frame(gy_error_est)

#
'%notin%' = Negate('%in%')
tagsout = list()
for(i in 1:length(rotationj$aj)){
  # need to change this to exclude the data used to build the transformation #
  # tagsout = rotationj[i,7:9] #
  # databin = datallwUTMjumbled[datallwUTMjumbled$tag %notin% tagsout,] #
  tagsout[[i]] = c(rotationj[i,7],rotationj[i,8],rotationj[i,9])
  gx_error_est[,i] = rotationj$aj[i]*datallwUTMjumbled$utmx + rotationj$bj[i]*datallwUTMjumbled$utmy + rotationj$cj[i]
  gy_error_est[,i] = rotationj$dj[i]*datallwUTMjumbled$utmx + rotationj$ej[i]*datallwUTMjumbled$utmy + rotationj$fj[i]
}

# make the first column the actuals 
gx_error_est = cbind(datallwUTMjumbled$tag,datallwUTMjumbled$gx,gx_error_est)
gy_error_est = cbind(datallwUTMjumbled$tag,datallwUTMjumbled$gy,gy_error_est)

# 
colnames(gx_error_est)[1:2] = c("tag","actual_gx")
colnames(gy_error_est)[1:2] = c("tag","actual_gy")

# remove NA's (produced by rows of datallwUTMjumbled that had no utm x and utm y to calculate from) #
gx_error_est = gx_error_est[!is.na(gx_error_est[,3]),]
gy_error_est = gy_error_est[!is.na(gy_error_est[,3]),]

# Tags to Character Variable # 
gx_error_est$tag = levels(gx_error_est$tag)[gx_error_est$tag]
gy_error_est$tag = levels(gy_error_est$tag)[gy_error_est$tag]


# Visual check on datasets #
View(gx_error_est)
dim(gx_error_est)
View(gy_error_est)
dim(gy_error_est)

#
rmse_gx = c()
rmse_gy = c()
rmse_sum = c()

for(i in 8579:ncol(gx_error_est)){
  
  # remove tags of trees that were used to fit the model #
  databinx = gx_error_est[gx_error_est$tag %notin% tagsout[[i]],]
  databiny = gy_error_est[gy_error_est$tag %notin% tagsout[[i]],]
  
  rmse_gx[i] = rmse(databinx$actual_gx,databinx[,i+2])
  rmse_gy[i] = rmse(databiny$actual_gy,databiny[,i+2])
  rmse_sum[i] = sum(rmse_gx[i],rmse_gy[i])
}

#
rmse_dat = data.frame(rmse_gx,rmse_gy,rmse_sum)
rmse_dat = rmse_dat[!is.na(rmse_dat$rmse_gx),]
View(rmse_dat[rmse_dat$rmse_sum == min(rmse_dat$rmse_sum),])

bestrow = rownames(rmse_dat[rmse_dat$rmse_sum == min(rmse_dat$rmse_sum),])

# 
# row #
View(gx_error_est[,(as.numeric(bestrow)+2)])
View(gy_error_est[,(as.numeric(bestrow)+2)])
# coefs #
best_coefs = rotationj[as.numeric(bestrow),1:6]
# check the best coefs #
best_coefs$aj

gx_test = best_coefs$aj*datallwUTMjumbled$utmx + best_coefs$bj*datallwUTMjumbled$utmy + best_coefs$cj
gy_test = best_coefs$dj*datallwUTMjumbled$utmx + best_coefs$ej*datallwUTMjumbled$utmy + best_coefs$fj

# check if the best estimates calculated from the best coefs line up #
# with the corresponding row in the gx/gy error est dataframe # 
View(gx_error_est[,(as.numeric(bestrow)+2)])
View(gx_test)
View(gy_error_est[,(as.numeric(bestrow)+2)])
View(gy_test)

# these appear to line up #
utmx_test = utmx_test[!is.na(datallwUTMjumbled$gx)]
utmx_test = utmx_test[!is.na(datallwUTMjumbled$utmx)]
plot(utmx_test,gx_error_est$X1)
unique(utmx_test == gx_error_est$X1)


#
save(best_coefs,file = paste(calcfolder,"best_coefs_toplotcoords_March31.RData",sep=""))
load(paste(calcfolder,"best_coefs_toplotcoords_March31.RData",sep=""))

#### Take in Canopy Height Data ###

file=paste(calcfolder,"canopyheight_utm.csv",sep="")
canopyheightdata = read.csv(file,header=TRUE)
options(digits=8)

canopyheightdata$gx = best_coefs$aj*canopyheightdata$POINT_X + best_coefs$bj*canopyheightdata$POINT_Y + best_coefs$cj
canopyheightdata$gy = best_coefs$dj*canopyheightdata$POINT_X + best_coefs$ej*canopyheightdata$POINT_Y + best_coefs$fj
plotheightdata = canopyheightdata[canopyheightdata$gx>=0&canopyheightdata$gy>=0,]
min(plotheightdata$gx)
max(plotheightdata$gy)

# mask out the canopy crane
min(plotheightdata$GRID_CODE)
max(plotheightdata$GRID_CODE)
dev.new()
hist(plotheightdata[plotheightdata$GRID_CODE > 69,"GRID_CODE"])
# mask out trees with height greater than 69 m 
plotheightdata[plotheightdata$GRID_CODE > 69, "GRID_CODE"] = NA

### CAIndex, Dstar, and eH Calculation ### 
data = read.table(paste(calcfolder,plotcode,"all.txt",sep=""),sep="\t",header=TRUE) 
dbhcol=5
#
pradius=10
prad=10
div=4
dstarV = exp(seq(log(50),log(2500),by=log(1.4)))
pa = pradius^2*pi

# height allometry from Hulshof, Swenson, and Weiser 2015 for gymnosperms # 
Hallom = 3.21
betaallom = 0.60

xmax = ceiling(max(data$gx))
ymax = ceiling(max(data$gy))

xV = seq(pradius+2,810-pradius-2,by=pradius*2/div)
yV = seq(pradius+2,330-pradius-2,by=pradius*2/div)
plot(xV,xV,type="n")

file=paste(calcfolder,pradius,"_",div,"mar11",".txt",sep="")
jj=0
for(x in xV){
  xdata = data[(data$gx-x)*(data$gx-x)<pradius*pradius,]
  for(y in yV){
    xydata = xdata[(xdata$gx-x)*(xdata$gx-x)+(xdata$gy-y)*(xdata$gy-y)<pradius*pradius,]
    for(census in seq(0,ncensuses-1)){
      jj=jj+1
      maxd=dstar=0
      if(dim(xydata)[1]>0){
        maxd=max(xydata[,dbhcol+census],na.rm=TRUE)
        dstar=find_dstar(xydata[,dbhcol+census],PA=pradius*pradius*pi)
        ddV = xydata[,dbhcol+census]; ddV = ddV[ddV>0]; ddV = ddV[!is.na(ddV)]
        CAindex = sum(ddV^theta*phi)/(pradius^2*pi)
        # height #
        if(x>pradius&x<xmax-pradius&y>pradius&y<ymax-pradius){
        dV = ddV[order(ddV,decreasing=TRUE)]
        #crown areas of these trees; not species specific, but could be
        CAv = phi*dV^theta
        #the crown area of each tree plus the crown areas of all those trees bigger than it
        CAtV = cumsum(CAv)
        #diameters of trees estimated to be in the canopy
        dcan = dV[CAtV<=pa]
        #convert to centimeters
        dcancm = dcan*0.1
        #vector of the heights of trees that are in the canopy 
        HcanV = Hallom*dcancm^betaallom
        #vector of the crown areas of the trees that are in the canopy
        CAcanV = phi*dcan^theta
        #if the canopy is not filled, there is also bareground height zero, with the leftover area
        if(max(CAtV)<pa){
          HcanV = c(HcanV,0)
          CAcanV = c(CAcanV,pa-max(CAtV))
        }
        eH = sum(HcanV*CAcanV)/(pa)
        }}
      if(x<pradius| x>xmax-pradius | y<pradius | y>ymax-pradius){
        print("analysis not valid on the edge, returning NA")
        eH = "NA"
      }
      if(maxd==-999) maxd=0
      write.table(data.frame(census,x,y,maxd,dstar,CAindex,eH),file,sep="\t",row.names=FALSE,col.names=jj==1,append=jj!=1)
    }
    points(x,y)
  }
}



### Hstar calculation ###

### Take in CAindex and Dstar Data ###
pradius = 10; div = 4
file=paste(calcfolder,pradius,"_",div,"mar11",".txt",sep="")
dstardata = read.table(file,sep="\t",header=TRUE)

# Take out the first census 
dstardata_0 = dstardata[dstardata$census == 0,]
#heightallomdata_0 = heightallomdata[heightallomdata$census == 0,]

# set up new columns
dstardata_0$max_ch = rep(1,dim(dstardata_0)[1])
dstardata_0$min_ch = rep(1,dim(dstardata_0)[1])
dstardata_0$mean_ch = rep(1,dim(dstardata_0)[1])


###### Calculate the Height of those Regions ####
xV = seq(pradius+2,810-pradius-2,by=pradius*2/div)
yV = seq(pradius+2,330-pradius-2,by=pradius*2/div)

jj=0
for(x in xV){
  ch_xdata = plotheightdata[(plotheightdata$gx-x)*(plotheightdata$gx-x)<pradius*pradius,]
  for(y in yV){
    ch_xydata = plotheightdata[(plotheightdata$gx-x)*(plotheightdata$gx-x)+(plotheightdata$gy-y)*(plotheightdata$gy-y)<pradius*pradius,]
    jj=jj+1
    if(dim(ch_xydata)[1]>0){
      dstardata_0[dstardata_0$x == x & dstardata_0$y == y,]$max_ch = max(ch_xydata$GRID_CODE,na.rm=T)
      dstardata_0[dstardata_0$x == x & dstardata_0$y == y,]$min_ch = min(ch_xydata$GRID_CODE,na.rm=T)
      dstardata_0[dstardata_0$x == x & dstardata_0$y == y,]$mean_ch = mean(ch_xydata$GRID_CODE,na.rm=T)
      }}}

i = 0
for(ptx in unique(dstardata_0$x)){
  i = i +1
  dstardata_0[dstardata_0$x==ptx,11]=i
  #heightallomdata_0[heightallomdata_0$x == ptx,6] = i
}
i = 0 
for(pty in unique(dstardata_0$y)){
  i = i +1
  dstardata_0[dstardata_0$y==pty,12]=i
  #heightallomdata_0[heightallomdata_0$y == pty,7] = i
}



# make a matrix with the data 
  plotmatrix_caindex=matrix(nrow=length(unique(dstardata_0$x)),ncol=length(unique(dstardata_0$y)))
  plotmatrix_dstar =matrix(nrow=length(unique(dstardata_0$x)),ncol=length(unique(dstardata_0$y)))
  plotmatrix_maxheight =matrix(nrow=length(unique(dstardata_0$x)),ncol=length(unique(dstardata_0$y)))
  plotmatrix_heightallom = matrix(nrow=length(unique(dstardata_0$x)),ncol=length(unique(dstardata_0$y)))
  
  for(i in seq(1,dim(dstardata_0)[1])){
    plotmatrix_caindex[dstardata_0[i,11],dstardata_0[i,12]] = (dstardata_0[i,6])
    plotmatrix_dstar[dstardata_0[i,11],dstardata_0[i,12]] = (dstardata_0[i,5])
    plotmatrix_maxheight[dstardata_0[i,11],dstardata_0[i,12]] = (dstardata_0[i,8])
    plotmatrix_heightallom[dstardata_0[i,11],dstardata_0[i,12]] = dstardata_0[i,7]
  }

  
  
# Plot the Plot Matrix #
library(fields)
brks = seq(0.38,1.8,length.out=10)
dev.new()
image.plot(plotmatrix_caindex,col=colorRampPalette(c("lightgreen","darkgreen"))(length(brks)-1),zlim=c(min(brks),max(brks)),breaks=brks,lab.breaks=round(brks,digits=2),xaxt="n",yaxt="n",bty="n")


# plot for height 
brks = seq(33,70,length.out=10)
dev.new()
image.plot(plotmatrix_maxheight,col=colorRampPalette(c("lightblue","darkblue"))(length(brks)-1),zlim=c(min(brks),max(brks)),breaks=brks,lab.breaks=round(brks),xaxt="n",yaxt="n",bty="n")

# dstart 
brks = seq(log(10),log(700),by = log(1.6))
dev.new()
image.plot(plotmatrix_dstar,col=colorRampPalette(c("lightblue","darkblue"))(length(brks)-1),zlim=c(min(exp(brks)),max(exp(brks))),breaks=exp(brks),lab.breaks=round(exp(brks)),xaxt="n",yaxt="n",bty="n")


# height allom
brks = seq(5,60,length.out=10)
dev.new()
image.plot(plotmatrix_heightallom,col=colorRampPalette(c("lightblue","darkblue"))(length(brks)-1),zlim=c(min(brks),max(brks)),breaks=brks,lab.breaks=round(brks),xaxt="n",yaxt="n",bty="n")













# now do some bootstrapping #
# update: too much bootstrapping, memory can't handle it # 

rmseboot = matrix(nrow = 1000,ncol = ncol(gx_error_est)-2)
rmseboot = data.frame(rmseboot)

for(i in 1:(ncol(gx_error_est)-2)){
  
  # remove tags of trees that were used to fit the model #
  databinx = gx_error_est[gx_error_est$tag %notin% tagsout[[i]],]
  databiny = gy_error_est[gy_error_est$tag %notin% tagsout[[i]],]
  
  for(j in 1:1000){
    # take a random sample of tags
    randomsample = sample(databinx$tag,size = 500)
    databinxBoot = databinx[databinx$tag%in%randomsample,]
    databinyBoot = databiny[databiny$tag%in%randomsample,]
   
    #calculate rmse for each bootstrapped set from x and y
    rmse_gxVar = rmse(databinxBoot$actual_gx,databinxBoot[,i+2])
    rmse_gyVar = rmse(databinyBoot$actual_gy,databinyBoot[,i+2])
    # take the sum of rmse for x and y and put it in the dataframe
    rmseboot[j,i] = sum(rmse_gxVar,rmse_gyVar)
  }
}


# not sure what this does: colnames(gx_error_est) <- c(colnames(gx_error_est)[1:2],1:998)

dataBins = list()

for(i in 1020:(ncol(gx_error_est)-2)){
  # remove tags of trees that were used to fit the model #
  dataX = gx_error_est[gx_error_est$tag %notin% tagsout[[i]],]
  dataY = gy_error_est[gy_error_est$tag %notin% tagsout[[i]],]
  
  dataBins[[i]] = cbind(dataX$tag,dataX$actual_gx,dataX[,i+2],dataY$actual_gy,dataY[,i+2])
  colnames(dataBins[[i]]) <-c("tag","actual_gx","pred_x","actual_gy","pred_y")
  
  #
}
 
save(dataBins,file = paste(calcfolder,plotcode,"utm_plot_modelPreds_1_1019",sep=""))
load(paste(calcfolder,plotcode,"utm_plot_modelPreds_1_1019",sep=""))

rmsebootfunc = function(dat){
  rmseBoot = matrix(nrow =1,ncol = 1)                   
  rmseBoot = data.frame(rmseBoot)
  for(j in 1:1000){
    randomsample = sample(dat[,1],size=500)
    datBoot = dat[dat[,1]%in%randomsample,]                       
    rmseX = rmse(as.numeric(datBoot[,2]),as.numeric(datBoot[,3]))
    rmseY = rmse(as.numeric(datBoot[,4]),as.numeric(datBoot[,5]))                                
    if(j==1){rmseBoot[1,1] = sum(rmseX,rmseY)}
    if(j>1) {rmseBoot = rbind(rmseBoot,sum(rmseX,rmseY))}
  }                                                                                                                                                                                                                                     
  return(rmseBoot)
}

fistthreeY = list(databinsY[[1]],databinsY[[2]],databinsY[[3]],databins[[4]],databins[[5]])


# test the next part of this on a smaller dataset
test = dataBins[c(1:40)]
testresults = lapply(test,rmsebootfunc)
testresults_df = do.call(cbind,testresults[c(1:40)])

View(testresults_df)

system.time(lapply(test, rmsebootfunc))
testresults_df = do.call(cbind,testresults[c(1:40)])










                                                                                                                                   
