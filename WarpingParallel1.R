library("parallel")

#         clusterExport(cl,c('Warping','x','i','u','Seg','Slack','Ix','MX','MT'))

Warping <- function(channel,x,i,u,Seg,Slack,Ix,MX,MT)  
{ 
  temp1 <- list()
  # benefit function 
  f1 <- function(xs,Seg,u,Slack,ts,te,MX,MT)
  {
    # xs : start warped position of segment on X 
    # xs + Seg + u : end warped position of segment on X
    # ts : start position of segment on T
    # te : end position of segment on T
    
    xe = xs + Seg + u
    
    bounds <- seq(xs,xe,length.out = length(ts:te))
    
    aux <- approx(xs:xe,MX[xs:xe],bounds)
    
    if (length(MT[ts:te]) == length(MX[xs:xe]) )
    {
      b <- MX[xs:xe]
    } else
    {
      b <- aux$y
    }
    a <- MT[ts:te]
    if ( all( c(var(a),var(b)) !=0 ) )
      correlation <- cor(a,b)
    else  correlation <- 0
    
    
    
    output <- list( corr = correlation, unwarped = a, warped = b)      
    return(output)
    
  }
  
  
  temp1 <- f1(x,Seg,u,Slack,Ix[i+1],Ix[i+2],MX[,channel],MT[,channel])
  return(list(corr = temp1$corr, unwarped = temp1$unwarped,
              warped = temp1$warped))
}



alignGCMS_Parallel <- function(MT,MX,Seg,Slack)
{
  # X:  query GC-MS signal
  # T: target GC-MS signal
  # Seg: Segment length
  # Slack: minimum warping       
  # F1: matrix containing the cumulated benefit function
  
  if(!is.matrix(MT)) MT <<- as.matrix(MT)
  if(!is.matrix(MX)) MX <<- as.matrix(MX)
  
  
  # Pre-aligming length of chromatogram
  # Intervals of 1 unit in P
  Lx <<- nrow(MX)-1
  
  # Post-aligning length of chromatogram and length of 
  # target chromatogram
  Lt <<- nrow(MT)-1
  
  # Calculate number of sections for X the query signal
  N <<- floor((Lx+1)/Seg)
  
  # Ix sequence of node positions in X (query signal) before warping
  Ix <<- round(seq(1,(Lx+1) ,length.out = (N+1)))
  
  # temp = (pX-1) %% LenSeg[1,1]
  # remainder of the segments and signals can have different lengths  
  
  temp <<- Lx %% Seg
  
  # N number of sections and N + 1 nodes
  Nnodes <<- length(Ix)
  if (Nnodes == N + 1)
  {
    if(Ix[Nnodes] != (Lx+1) ) Ix[Nnodes] <<- Ix[Nnodes-1]+temp
  } else if (N-1 - Nnodes == 1)
    Ix[Nnodes + 1] <<- Lx + 1
  
  # calculate difference in mean section length between P and T
  d <<- floor(Lt/N) - Seg
  
  F1 <- matrix(0.0, nrow = (N+1), ncol = Lt + 1)
  U <- matrix(0.0, nrow = (N+1), ncol = Lt + 1)
  
  for (i in 1:(N+1))
  {
    for(x in 1:(Lt + 1))
    {
      F1[i,x] <- -Inf  
    }
  }
  
  F1[N+1,] <- 0
  
  # W is the list that contain the warping signal for warping nodes
  aux1W <- list()
  temp1 <- list()
  temp_unwarped <- list()
  temp_warped <- list()
  a <- rep(0, ncol(MT)*nrow(MT))
  b <- rep(0, ncol(MT)*nrow(MT))
  
  for (i_i in (N-1):0)
  {
    #xstart: minimum start point of segment i in query signal
    #xend: maximum end point of segment i in query signal
    xstart <- max(1+i_i*(Seg+d-Slack),(Lt+1)-(N-i_i)*(Seg+d+Slack))
    xend <- min(1+i_i*(Seg+d + Slack),(Lt+1)-(N-i_i)*(Seg+d-Slack))
    cat("\n i=",i,"\\(", xstart, "-" ,xend, "\\) \n")
    
    for (i_x in xstart:xend)
    {
      for (i_u in (d-Slack):(d+Slack))
      {
        if( i_x+Seg+i_u <= Lt + 1 )
        {
          
          
          i <<-i_i
          x <<- i_x
          u <<- i_u
          
          
          results <- pvec(1:ncol(X), function(channel) {
            temp1[[channel]] <- Warping(channel,x,i,u,Seg,Slack,Ix,MX,MT)
          },mc.cores = 1) 
          
          
          for(channel in 1:ncol(MX))
          {
            temp_warped[[channel]] <- results[[channel]]$warped
            temp_unwarped[[channel]] <- results[[channel]]$unwarped
          }
          
          a <- unlist(temp_unwarped)
          b <- unlist(temp_warped)
          
          
          
          salida <- list( unwarped = a, warped = b)
          ruta <- "/home/pgrad1/2201449s/R/GC_MS/Output/" 
          saveRDS(salida,paste0(ruta,"Seg_",i_i,"_node_",i_x,"_warping_",i_u,".RDS"))
          
          if ( all( c(var(a),var(b)) !=0 )) corr <- cor(a,b) else  corr <- 0
          
          #   fsum <- F1[i+2,x+Seg+u] + f(x,Seg,u,Slack,Ix[i+1],Ix[i+2],X,T)
          fsum <- F1[i_i+2,i_x+Seg+i_u] + corr
          
          cat("\nnode=",i_i,",x =",i_x,",u=",i_u,",fsum=",fsum,"\n")
          if (fsum > F1[i_i+1,i_x])
          {
            F1[i_i+1,i_x] <- fsum
            U[i_i+1,i_x] <- i_u
            aux1W[[i_i+1]] <- b
          }# End-if
          
        }# End-if
        
      }# End for u
      
    } # End for x
    
  } # End for i
  
  #  F1 <- r$F1
  #  U <- r$U
  
  F1
  U
  Xw <- rep(0,N+1)
  u <- rep(0,N+1)
  aux2W <- list()
  Xw[1]<-1
  Xw[N+1]<-ifelse(!is.matrix(MX),length(MX),nrow(MX))
  
  for(i in 1:(N-1) )
  {
    u[i] <- U[i,Xw[i]]
    Xw[i+1] <- Xw[i]+Seg+u[i]
  }
  
  for(i_channel in 1:ncol(X) )
  {
    ini <- 1
    final <- 0
    for(i in 1:length(aux1W))
    {
      nelements <- length(aux1W[[i]][[i_channel]])
      if(i<length(aux1W))
      {
        final <- final + nelements - 1
        aux2W[ini:final, i_channel] <- unlist( aux1W[[i]] [[i_channel]][1:(nelements-1)] )
        ini <- final + 1
      } else  
      {
        final <- final + nelements
        aux2W[ini:final, i_channel] <- unlist(aux1W[[i]][[i_channel]][1:(nelements)] )
      } # end else
    }# end for i
  }# end for i_channel  
    return(list(X=Ix,Xw = Xw, u=u,F1=F1,U=U,W = aux2W))
  }
  