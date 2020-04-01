### Plot on the UTM grid ### 
library(dplyr); library(tidyr); library(Metrics); library(fields)

# 
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

# original plot data #
datallwUTM = read.table(paste(calcfolder,plotcode,"allwUTM.txt",sep=""),sep="\t",header=TRUE)
# load rotation coefficients # 
load(paste(calcfolder,"best_coefs_toplotcoords_March31.RData",sep=""))

#### Take in Canopy Height Data ###
file=paste(calcfolder,"canopyheight_utm.csv",sep="")
canopyheightdata = read.csv(file,header=TRUE)
options(digits=8)

canopyheightdata$gx = best_coefs$aj*canopyheightdata$POINT_X + best_coefs$bj*canopyheightdata$POINT_Y + best_coefs$cj
canopyheightdata$gy = best_coefs$dj*canopyheightdata$POINT_X + best_coefs$ej*canopyheightdata$POINT_Y + best_coefs$fj

# mask out height data that isn't in the plot # 
plotheightdata = canopyheightdata[canopyheightdata$gx>=0&canopyheightdata$gy>=0,]
#
plotheightdata = plotheightdata[plotheightdata$gx<=800.1&plotheightdata$gy<=325,]

################################
### plot actual height values ##
################################

plotheightdata = plotheightdata[with(plotheightdata, order(XCoord, YCoord)), ]

plotmatrix_actualheight = matrix(nrow=length(unique(plotheightdata$XCoord)),ncol=length(unique(plotheightdata$YCoord)))
i = 0
for(ptx in unique(plotheightdata$XCoord)){
  i = i +1
  plotheightdata[plotheightdata$XCoord==ptx,8]=i
}
plotheightdata = plotheightdata[with(plotheightdata, order(YCoord, XCoord)), ]
i = 0 
for(pty in unique(plotheightdata$YCoord)){
  i = i +1
  plotheightdata[plotheightdata$YCoord==pty,9]=i
}

plotheightdata = plotheightdata[with(plotheightdata, order(XCoord, YCoord)), ]

for(i in seq(1,dim(plotheightdata)[1])){
  plotmatrix_actualheight[plotheightdata[i,8],plotheightdata[i,9]] = (plotheightdata[i,3])
  
}

brks = seq(0,70,length.out=20)
dev.new()
image.plot(plotmatrix_actualheight,col=colorRampPalette(c("lightblue","darkblue"))(length(brks)-1),zlim=c(min(brks),max(brks)),breaks=brks,lab.breaks=round(brks),xaxt="n",yaxt="n",bty="n")


### CAIndex, Dstar, and eH Calculation ### 
dbhcol=7
#
pradius=10
prad=10
div=20
dstarV = exp(seq(log(50),log(2500),by=log(1.4)))
pa = pradius^2*pi

# height allometry from Hulshof, Swenson, and Weiser 2015 for gymnosperms # 
Hallom = 3.21
betaallom = 0.60

xmax = ceiling(max(datallwUTM$utmx,na.rm=T))
ymax = ceiling(max(datallwUTM$utmy,na.rm=T))
xmin = floor(min(datallwUTM$utmx,na.rm=T))
ymin = floor(min(datallwUTM$utmy,na.rm=T))

#xV = seq(xmin+pradius+2,xmax-pradius-2,by=pradius*2/div)
#yV = seq(ymin+pradius+2,ymax-pradius-2,by=pradius*2/div)

# Use the Original Plot Boundaries in Plot Coordinates and Convert to UTM #
xV = seq(pradius+2,810-pradius-2,by=pradius*2/div)
yV = seq(pradius+2,330-pradius-2,by=pradius*2/div)
minxV = min(xV)
maxxV = max(xV)
minyV = min(yV)
maxyV = max(yV)

# Convert those Coordinates to UTM coordinates # 
# Load Coefficients # 
load(paste(calcfolder,"best_coefs_UTMtoplotcoords_March29.RData",sep=""))

# calculate # 
# bottom left # 
x1 = best_coefs_utm$aj*minxV + best_coefs_utm$bj*minyV + best_coefs_utm$cj
y1 = best_coefs_utm$dj*minxV + best_coefs_utm$ej*minyV + best_coefs_utm$fj

# top left #
x2 = best_coefs_utm$aj*minxV + best_coefs_utm$bj*maxyV + best_coefs_utm$cj
y2 = best_coefs_utm$dj*minxV + best_coefs_utm$ej*maxyV + best_coefs_utm$fj

# bottom right #
x3 = best_coefs_utm$aj*maxxV + best_coefs_utm$bj*minyV + best_coefs_utm$cj
y3 = best_coefs_utm$dj*maxxV + best_coefs_utm$ej*minyV + best_coefs_utm$fj

# top right # 
x4 = best_coefs_utm$aj*maxxV + best_coefs_utm$bj*maxyV + best_coefs_utm$cj
y4 = best_coefs_utm$dj*maxxV + best_coefs_utm$ej*maxyV + best_coefs_utm$fj

# Re make xV and yV on the UTM Grid # 
#xmax = ceiling(max(datallwUTM$utmx,na.rm=T))
xmax = max(plotheightdata$XCoord)
#ymax = ceiling(max(datallwUTM$utmy,na.rm=T))
ymax = max(plotheightdata$YCoord)
#xmin = floor(min(datallwUTM$utmx,na.rm=T))
xmin = min(plotheightdata$XCoord)
#ymin = floor(min(datallwUTM$utmy,na.rm=T))
ymin = min(plotheightdata$YCoord)

# set xV and yV so that they line up with the LiDAR xV and yV
xV = seq(xmin,xmax,by=pradius*2/div)
yV = seq(ymin,ymax,by=pradius*2/div)

plot(xV,xV)


file=paste(calcfolder,pradius,"_",div,"april1_utm",".txt",sep="")
jj=0
for(x in xV){
  xdata = datallwUTM[(datallwUTM$utmx-x)*(datallwUTM$utmx-x)<pradius*pradius,]
  xdata = xdata[!is.na(xdata$tag),]
  for(y in yV){
    xydata = xdata[(xdata$utmx-x)*(xdata$utmx-x)+(xdata$utmy-y)*(xdata$utmy-y)<pradius*pradius,]
    xydata = xydata[!is.na(xydata$tag),]
    for(census in seq(0,ncensuses-1)){
      jj=jj+1
      maxd=dstar=0
      if(y >= (((y3-y1)/(x3-x1))*(x-x1)) + y1 & y >= (((y1-y2)/(x1-x2))*(x-x1)) + y1 & y <= (((y4-y3)/(x4-x3))*(x-x3)) + y3 & y <= (((y4-y2)/(x4-x2))*(x-x2)) + y2){
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
        }}}
      if(y <= (((y3-y1)/(x3-x1))*(x-x1)) + y1 | y <= (((y1-y2)/(x1-x2))*(x-x1)) + y1 | y >= (((y4-y3)/(x4-x3))*(x-x3)) + y3 | y >= (((y4-y2)/(x4-x2))*(x-x2)) + y2){
        print("analysis not valid on the edge, returning NA")
        eH = "NA"
        CAindex = "NA"
        dstar = "NA"
      }
      if(maxd==-999) maxd=0
      write.table(data.frame(census,x,y,maxd,dstar,CAindex,eH),file,sep="\t",row.names=FALSE,col.names=jj==1,append=jj!=1)
    }
  }
}




heightgrid = read.table(file,header=TRUE,sep="\t")
heightgrid_0 = heightgrid[heightgrid$census == 0,]

### plot the height grid ## 
plotmatrix_heightgrid = matrix(nrow=length(unique(heightgrid_0$x)),ncol=length(unique(heightgrid_0$y)))
i = 0
for(ptx in sort(unique(heightgrid_0$x))){
  i = i +1
  heightgrid_0[heightgrid_0$x==ptx,8]=i
}

i = 0 
for(pty in sort(unique(heightgrid_0$y))){
  i = i +1
  heightgrid_0[heightgrid_0$y==pty,9]=i
}


for(i in seq(1,dim(heightgrid_0)[1])){
  plotmatrix_heightgrid[heightgrid_0[i,8],heightgrid_0[i,9]] = (heightgrid_0[i,7])
  
}

brks = seq(0,70,length.out=20)
dev.new()
image.plot(plotmatrix_heightgrid,col=colorRampPalette(c("lightblue","darkblue"))(length(brks)-1),zlim=c(min(brks),max(brks)),breaks=brks,lab.breaks=round(brks),xaxt="n",yaxt="n",bty="n")


# plot the canopy layers matrix # 
plotmatrix_canopylayers = matrix(nrow=length(unique(heightgrid_0$x)),ncol=length(unique(heightgrid_0$y)))
for(i in seq(1,dim(heightgrid_0)[1])){
  plotmatrix_canopylayers[heightgrid_0[i,8],heightgrid_0[i,9]] = (heightgrid_0[i,6])
  
}

brks = seq(0.029,2.2,length.out=10)
dev.new()
image.plot(plotmatrix_canopylayers,col=colorRampPalette(c("lightgreen","darkgreen"))(length(brks)-1),zlim=c(min(brks),max(brks)),breaks=brks,lab.breaks=round(brks,2),xaxt="n",yaxt="n",bty="n")



#### Take in Canopy Height Data ###
file=paste(calcfolder,"canopyheight_utm.csv",sep="")
canopyheightdata = read.csv(file,header=TRUE)
options(digits=8)

canopyheightdata$gx = best_coefs$aj*canopyheightdata$POINT_X + best_coefs$bj*canopyheightdata$POINT_Y + best_coefs$cj
canopyheightdata$gy = best_coefs$dj*canopyheightdata$POINT_X + best_coefs$ej*canopyheightdata$POINT_Y + best_coefs$fj

# mask out height data that isn't in the plot # 
plotheightdata = canopyheightdata[canopyheightdata$gx>=0&canopyheightdata$gy>=0,]
#
plotheightdata = plotheightdata[plotheightdata$gx<=800.1&plotheightdata$gy<=325,]

################################
### plot actual height values ##
################################

plotheightdata = plotheightdata[with(plotheightdata, order(XCoord, YCoord)), ]

# above or below the lines # 
colnames(plotheightdata) <- c("x","y",colnames(plotheightdata)[3:9])
plotheightdata$inboundaries = plotheightdata$y >= (((y3-y1)/(x3-x1))*(plotheightdata$x-x1)) + y1 & plotheightdata$y >= (((y1-y2)/(x1-x2))*(plotheightdata$x-x1)) + y1 & plotheightdata$y <= (((y4-y3)/(x4-x3))*(plotheightdata$x-x3)) + y3 & plotheightdata$y <= (((y4-y2)/(x4-x2))*(plotheightdata$x-x2)) + y2
plotheightdata$height = ifelse(plotheightdata$inboundaries == TRUE, plotheightdata$GRID_CODE,NA)



plotmatrix_actualheight = matrix(nrow=length(unique(plotheightdata$x)),ncol=length(unique(plotheightdata$y)))
i = 0
for(ptx in unique(plotheightdata$x)){
  i = i +1
  plotheightdata[plotheightdata$x==ptx,8]=i
}
plotheightdata = plotheightdata[with(plotheightdata, order(y, x)), ]
i = 0 
for(pty in unique(plotheightdata$y)){
  i = i +1
  plotheightdata[plotheightdata$y==pty,9]=i
}

plotheightdata = plotheightdata[with(plotheightdata, order(x, y)), ]

for(i in seq(1,dim(plotheightdata)[1])){
  plotmatrix_actualheight[plotheightdata[i,8],plotheightdata[i,9]] = (plotheightdata[i,11])
  
}

brks = seq(0,70,length.out=20)
dev.new()
image.plot(plotmatrix_actualheight,col=colorRampPalette(c("lightblue","darkblue"))(length(brks)-1),zlim=c(min(brks),max(brks)),breaks=brks,lab.breaks=round(brks),xaxt="n",yaxt="n",bty="n")

