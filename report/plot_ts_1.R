library(PBSsatellite)

pdf(file="plot_ts_1.pdf",
    width=6.6, height=2.5)

## load ncdfData object
data(sst)

## create a PolySet with a polygon
## for each hemisphere
polys <- makeGrid(
    x=c(0, 360), y=c(-90, 0, 90),
    addSID=FALSE, projection="LL")

## create a time series object
tsList <- extractTimeSeries(
    sst, polygons=polys)

## convert the time series object into
## a data frame
tsDF <- listToDF(tsList)

## set up some plot parameters
par(mgp=c(1.7, 0.4, 0),
    mar=c(2.7, 2.7, 1.5, 1.5),
    tck=c(-0.02), cex=0.9)
    
## plot mean for the southern hemisphere
## first (without x axis)
plot(tsDF[tsDF$PID == 1, "mean"],
     type='b', xaxt='n', col='blue',
     ylim=c(11, 16),
     ylab=attributes(sst)$dataUnits)
lines(tsDF[tsDF$PID == 2, "mean"],
      type='b', col='red')

## create appropriate x-axis labels and
## add a title
axis(1,
     at=1:nrow(tsDF[tsDF$PID == 1, ]),
     lab=as.Date(names(tsList)))
title(main=attributes(sst)$dataType)

dev.off()
