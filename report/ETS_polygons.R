library(PBSsatellite)

pdf(file="ETS_polygons.pdf", width=5,
    height=2.7)

data(sst)
data(worldLL)

plot(sst, slice=1,
     plt=c(.08, .98, .12, 1),
     mgp=c(1.7, 0.4, 0))
addPolys(worldLL)

## create 3 random polygons for demo
polys <- data.frame(
    PID=c(rep(1, 4), rep(2, 4),
        rep(3, 4)),
    POS=c(1:4, 1:4, 1:4),
    X=c(155, 160, 150, 180,  0,
         20,  20,   0,  45, 65,
         55,  35),
    Y=c( 75,  50,  10,  85, 20,
         20,  40,  40,  50, 50,
         65,  65))

## creates a PolySet
polys <- as.PolySet(polys,
                    projection="LL")

addPolys(polys, col="blue")

dev.off()
