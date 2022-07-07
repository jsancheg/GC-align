
source("GC_Sim.R")
ruta<-getwd()
Data <- readRDS(paste0(ruta,"/Data/sample_data.rds"))

str(Data)
Data[1,,]
dim(Data[1,1:2900,])
plot(Data[1,,100],type = "l")
plot(Data[2,,100],type = "l")

# set 1st signal GC-MS as a target
T <- Data[1,,]
# set 2nd signal GC-MS as a query
X <- Data[2,,]

Seg <- 15
Slack <- 5

tic("alignGCMS")
Xaligned <- alignGCMS(T,X,Seg,Slack)
toc()

T<-Data[1,,100]
X<-Data[2,,100]

Seg = 5
Slack = 1
Options = c(0,1,1,0,0)

cow(T,X,Seg,Slack,Options)

length(T)
length(X)

#align(T,X,10,5)



T <- apply(Data[1,,],2,sum)
X <- apply(Data[2,,],2,sum)
X1 <- apply(Data[3,,],2,sum)
X2 <- apply(Data[4,,],2,sum)
X3 <- apply(Data[5,,],2,sum)


plot(1:250,T,type = "l", ylim = c(0,max(T,X,X1,X2,X3)))
lines(1:250,X, col = "red")
lines(1:250,X1, col = "blue")
lines(1:250,X2, col = "green")
lines(1:250,X3, col = "brown")

cow(T,X,Seg,Slack,Options)
WX1 <- align(T,X,10,5)
WX2 <- align(T,X1,10,5)
WX3 <- align(T,X2,10,5)
WX4 <- align(T,X3,10,5)


all(WX1$Xw==WX2$Xw)
all(WX1$Xw==WX3$Xw)
all(WX1$Xw==WX4$Xw)


WX1$Xw
WX2$Xw
WX3$Xw
WX4$Xw
