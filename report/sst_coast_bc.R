library(PBSsatellite)

pdf(file="sst_coast_bc.pdf", width=2,
        height=2.55)

## create simple polygon for the BC coast
bcCoast <- data.frame(PID=c(rep(1, 7)), POS=c(1:7),
           X=c(223, 226, 235, 238,
               238, 226, 223),
           Y=c(58, 53, 48, 48,
               50, 60, 59.5))

## created simple BC Coast Polygon
bcCoast <- as.PolySet(bcCoast, projection="LL")

## import world Polygons from PBS Mapping
data(worldLL)

## creates the complex BC Coast polygon
bcComplex <- joinPolys(bcCoast, worldLL,
                      operation="DIFF")

## import world Sea Surface Temperature ncdfData object
data(sst)

## creates SST BC Coast ncdfData object
sstBcCoast <- clipRegion(sst, polygons=bcComplex)

par(mgp=c(1.7, 0.4, 0))

## complex PolySet
plot(sstBcCoast, xlim=range(bcCoast$X), ylim=range(bcCoast$Y),
     projection="LL", plt=c(.20, .98, .16, .98))

## add complex coast polygon for visualization
addPolys(worldLL, col="beige")
addPolys(bcCoast)

dev.off()
