source(paste0("E:/University of Glasgow/AgnieskaPRoject/GC-align","/WarpingParallel1.R"))

# load signals

# artificial signals 1 channel
T1 = c(0,1,1,2,5,6,5,5,4,4,3,2,1,1,0,0)
X1 = c(0,0,1,1,1,2,3,5,5,6,6,5,3,1,1,0)
X2  = c(0,0,1,1,1,2,3,5,5,6,6,5,2,1,1,0)

# artificial signals 2 channel GC-MS
MT <- cbind(T1,T1)
MX <- cbind(X1,X2)

Seg = 5
Slack = 1
Options = c(0,1,1,0,0)

WX <- alignGCMS_Parallel(MT,MX,Seg,Slack)
