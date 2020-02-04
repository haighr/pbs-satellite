library(PBSsatellite)

w <- 2.2
h <- 3

mar <- c(2, 2, .1, 0.1)
mgp <- c(1.7, 0.4, 0)
##-----------------------------------------------------------------------------

pdf(file="trivial.pdf", width=w, height=h)

d <- list()
d$`2017-06-16` <- list()
d$`2017-06-16`$data <- matrix(c(1, 1, 2, 3),
			      nrow=2, byrow=FALSE)
attr(d, "x") <- c(-128, -127)
attr(d, "y") <- c(49, 48)
attr(d, "dataType") <- "Sample"
attr(d, "dataUnits") <- "none"
attr(d, "class") <- "ncdfData"

par(mar=mar)
par(mgp=mgp)

plot(d, slice=1)

dev.off()
