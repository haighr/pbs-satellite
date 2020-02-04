library(PBSsatellite)

pdf(file="coast_bc_joinPoly.pdf", width=2,
        height=2.55)

## create simple polygon for the BC coast
bcCoast <- data.frame(PID=c(rep(1, 7)), POS=c(1:7),
                      X=c(223, 226, 235, 238,
                          238, 226, 223),
                      Y=c(58, 53, 48, 48,
                          50, 60, 59.5))

bcCoast <- as.PolySet(bcCoast, projection="LL")

## import world Polygons from PBS Mapping
data(worldLL)

## join the two PolySets by using the DIFF operation
## to get the difference between the two polygons
bcJoinedCoast <- joinPolys(bcCoast, worldLL,
                           operation="DIFF")
par(mgp=c(1.7, 0.4, 0))

## complex PolySet
plotMap(NULL,
        xlim=range(bcCoast$X), ylim=range(bcCoast$Y),
        projection="LL", plt=c(.20, .98, .16, .98),
        tck = -0.014)
addPolys(bcJoinedCoast, col="yellow")

dev.off()
