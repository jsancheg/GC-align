library(stats4)
library(methods)
# creating the class signal.
# 
rm(list = ls())
# Create a class with three fields.
# peak: A vector that contains the position of peaks in the range of spectra.
# variance: A vector that contains the variances for different peaks.
# intensity: A vector that contains the intensity of the peaks.
# Window: A vector that contains the number of retention time to be affected
#         around the main peak.
# lowerlimit: A variable that set the minimum value for retention time.
# upperlimit: A variable that set the maximum value for retention time.

gc <- setRefClass("gc",fields = list(peak = "numeric",
                                   variance = "numeric",
                                   intensity = "numeric",
                                   window = "numeric",
                                   lowerlimit = "numeric",
                                   upperlimit = "numeric"),
                     methods = list(
                       print = function()
                       {
                         if(!any(length(peak) == length(variance) | 
                                 length(variance) == length(intensity) |
                                 length(peak == length(intensity)))) 
                           stop("vector peak, variance, and intensity must be same length")
                          
                         # number of peaks
                         np <- length(peak)
                         # create the vector for all the range of retention time
                          x <- rep(0,length(seq(lowerlimit,upperlimit,1)))
                          
                          # for each peak identify the window to be changed
                          for(j in 1:np)
                            {
                            if(peak[j]-window[j]>0)
                            {
                              ini<- peak[j]-window[j]
                              fin<- peak[j]+window[j]
                              # obtain position to be changed
                              positionToChange <- seq(ini,fin,1)
                              
                            x[ini:fin]<-x[ini:fin] + sapply(positionToChange,function(pos)
                              {
                                # generate a value from normal dist. with mean equal to
                                # intensity
                                abs(rnorm(1,intensity[j],sqrt(variance[j])))/exp(abs(peak[j]-pos)/2)
                                # exponential decay when the retention time differ from 
                                # the position of the peak
                                
                              })                  
                              
                            }
                          }
                        print = x 
                       },
                       range = function()
                       {
                         range = seq(lowerlimit,upperlimit,1)
                       }
                     ))



set.seed(1)
v1 <- c(50,110,140,175,190,230,250)
u <- round(rnorm(7,5,2),0)
u
v2 <- v1 + u
v2

gc1<-gc(peak =v1,
               variance =c(2,2,1,1,1,1,1),
               intensity = c(10,35,55,35,180,700,400),
               window = c(3,3,3,5,5,10,10),
               lowerlimit = 0,
               upperlimit = 500)

# signal shifted by adding a random normal vector to the peak vectors
gc2<-gc(peak =v2,
                variance =c(2,2,1,1,1,1,1),
                intensity = c(10,35,55,35,180,700,400),
                window = c(3,3,3,5,5,10,10),
                lowerlimit = 0,
                upperlimit = 500)

B1 <- gc1$print()
# signal shifted
T1 <- gc2$print()

plot(gc1$range(), B1, type = "l", col = "blue")
lines(gc1$range(),T1, col = "green")

