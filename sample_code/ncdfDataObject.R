#possible ncdfData object after scaleRegion has been called
#to show all 3 layers of a slice (data,err,miss)

#setting attributes for a ncdfData object
#more attributes to come
dataType <- "SST"
units <- "Celsius"
lat <- c(47,46,45)
lon <- c(-150,-149,-148,-147,-146,-145)
missval = 32767

#creating list of attributes
attributes <- list(dataType=dataType, units=units, lat=lat, lon=lon ,
                   missval=missval)

#creation of slices
m2014.01.07 <- list(data=matrix(c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18)
                    ,nrow=3, ncol=6),err=matrix(,nrow=3, ncol=6), 
                    miss=matrix(,nrow=3, ncol=6))

m2014.01.14 <- list(data=matrix(c(19,20,21,22,23,24,25,26,27,28,29,30,31,32,33
                    ,34,35,36),nrow=3, ncol=6),err=matrix(,nrow=3, ncol=6),
                    miss=matrix(,nrow=3, ncol=6))

#creating list of existing slices
slices <- list(m2014.01.07=m2014.01.07, m2014.01.14=m2014.01.14)

#creation of the ncdfData object
ncdfData <- list(attributes=attributes,slices=slices)

#print ncdfData object
ncdfData

#output data/err from a particular slice
ncdfData$slices$m2014.01.07

#output data from slice
ncdfData$slices$m2014.01.07$data

#how many slices in a ncdfData object
length(ncdfData$slices)

#possible timeSeries object

#attributes for a time series
#more to come, might not need these
dataType <- "SST"
dateRange <- c("2014-01-07","2014-01-14")

#might be useful to know the data in the date frame
#there might be a way to find this out without storing this.
#Could be useful when a user only specifies a subset
#of functions e.g., c("sum","mean"), or when they do not 
#specify a polyset but instead an xlim and ylim, in this case
#there would be no PID in the data vector / frame.
data <- c("PID","sum","mean","sd","missing")

#make a list of the attributes
attributes <- list(dataType=dataType, dateRange=dateRange, data=data)

#create data for a time series
PID <- c(1,2,3,4)
sum <- c(4500,3400,2300,4300)
mean <- c(40,37,24,40)
sd <- c(.47,.34,-.20,.45)
missing <- c(20,43,80,12)

#create a data frame
a <- data.frame(PID=PID,sum=sum,mean=mean,sd=sd,missing=missing)


#create data for a second time series
PID <- c(1,2,3,4)
sum <- c(4600,4400,2500,5300)
mean <- c(41,40,25,50)
sd <- c(.44,.40,-.18,1.1)
missing <- c(22,41,79,14)

b <- data.frame(PID=PID,sum=sum,mean=mean,sd=sd,missing=missing)

#get name of slices
slices <- names(ncdfData$slice)

#create seriesData
seriesData <- list(a,b)

#name slices in seriesData
names(seriesData) <- slices

#create timeSeries object
timeSeries <- list(attributes=attributes,seriesData=seriesData)

#print out timeSeries object
timeSeries

