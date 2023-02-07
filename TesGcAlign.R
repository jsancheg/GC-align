AgnesProjPath <- "E:/University of Glasgow/AgnieskaPRoject/GC-align/"
setwd(AgnesProjPath)
source(paste0(AgnesProjPath,"GC_Sim.R"))


T1 = c(0,1,1,2,5,6,5,5,4,4,3,2,1,1,0,0)
X1 = c(0,0,1,1,1,2,3,5,5,6,6,5,3,1,1,0)
X2  = c(0,0,1,1,1,2,3,5,5,6,6,5,2,1,1,0)


slack = 1
Seg = 5
Slack = 1
Options = c(0,1,1,0,0)

tic("cow")
cow(T1,X1,Seg,Slack,Options)
toc()

# align signal X1
tic("align")
align(T1,X1,Seg,Slack)
toc()


tic("alighGCMS")
WX1 <- alignGCMS(T1,X1,Seg,Slack)
toc()
WX1$W

tic("alighGCMS_Parallel")
WX1_Par <- alignGCMS_Parallel(T1,X1,Seg,Slack)
toc()
WX1_Par$W



# align signal X2
align(T1,X2,Seg,Slack)
WX2 <- alignGCMS(T1,X2,Seg,Slack)



MX <- cbind(X1,X2)
MT <- cbind(T1,T1)


# align signal with cow
tic("align")
cow(T1,X2,Seg,Slack,Options)
toc()


# align signal X2
tic("align")
align(T1,X2,Seg,Slack)
toc()

tic("AlignGCMS")
WMX <- alignGCMS(MT,MX,Seg,Slack)
toc()
WMX$W

tic("alignGCMS_Parallel")
WMXPar <- alignGCMS_Parallel(MT,MX,Seg,Slack)
toc()

WMX$Xw
WMX$W



WMXPar$Xw


# align signal X1 and X2


cow(T,X,Seg,Slack,Options)



