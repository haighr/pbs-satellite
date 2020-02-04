library(PBSsatellite)

w <- 2.75
h <- 1.84
plt <- c(.14, .99, 0.14, 1)
##-----------------------------------------------------------------------------

pdf(file="sst_a.pdf", width=w, height=h)

data(sst)
data(worldLL)

plot(sst, slice=1,
     plt=plt,
     mgp=c(1.7, 0.4, 0))

dev.off()

##-----------------------------------------------------------------------------

pdf(file="sst_b.pdf", width=w, height=h)

data(sst)
data(worldLL)

g <- makeGrid(x=c(0, 360), y=c(-90, 0, 90), addSID=F)
plotMap(g, projection="LL",
        col=adjustcolor(c("blue", "red"), alpha.f=0.5),
        plt=plt, tck = -0.014,
        mgp=c(1.7, 0.4, 0))

dev.off()

##-----------------------------------------------------------------------------

pdf(file="sst_c.pdf", width=w, height=h)

data(sst)
data(worldLL)

plot(sst, slice=1,
     plt=plt,
     mgp=c(1.7, 0.4, 0))
addPolys(g, col=adjustcolor(c("blue", "red"), alpha.f=0.5))

dev.off()
