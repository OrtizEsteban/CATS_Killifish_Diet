###############
# (1) EXPANDE #
###############

expande<-function(L,Q,R){
  out<-NULL
  id<-1:nrow(L)
  for (i in id){
    Ni<-L[i,]
    Ri<-matrix(R[i,],nrow = ncol(L),ncol=ncol(R),byrow = T)
    colnames(Ri)<-colnames(R)
    out.t<-cbind(Q,Ri,Ni)
    rownames(out.t)<-paste(colnames(L),"comm_",i,sep = "_")
    out<-rbind(out,out.t)
  }
  RQ<-NULL
  nombres<-NULL
  for(q in 1:ncol(Q)){
    for(r in 1:ncol(R)){
      RQ<-cbind(RQ,out[,q]*out[,(r+ncol(Q))])
      nombres<-c(nombres,paste(colnames(out)[q],"_x_",colnames(out)[(r+ncol(Q))], sep="_"))
    }
  }
  N<-out[,ncol(out)]
  out<-out[,-ncol(out)]
  colnames(RQ)<-nombres
  out<-cbind(out,RQ,N)
  out
}

######################
# (2) PREY BODY SIZE #
######################

prey.bs.diet<-function(M, col.prey, col.bs, col.diet){
  
  id.prey<-sort(unique(M[, col.prey])) # Da vector con nombres de items troficos
  
  ## mean.bs<-NULL
  ##  min.bs<-NULL
  ## max.bs<-NULL
  ## var.bs<-NULL
  ## quantil.bs<-NULL
  bs.out<-NULL
  
  
  for(i in id.prey){
    
    M.i<-M[which(M[,col.prey]==i),] #Sub-matriz para item trofico i
    #return(M.i)
    if(length(which(is.na(M.i[,col.bs])==T))>0) {
      
      
      ##  mean.bs<-c(mean.bs,NA)
      ##  min.bs<-c(min.bs,NA)
      ##  max.bs<-c(max.bs,NA)
      ##  var.bs<-c(var.bs,NA)
      ##  quantil.bs<-c(quantil.bs, rep(NA, length(seq(0,1,0.1))))
      
      bs.out.i<-c(NA, NA, NA, NA, rep(NA, length(seq(0,1,0.1))))
      bs.out<-rbind(bs.out, bs.out.i)
      
    } else {
      
      mean.bs.i<-mean(M.i[,col.bs])
      min.bs.i<-min(M.i[,col.bs])
      max.bs.i<-max(M.i[,col.bs])
      var.bs.i<-var(M.i[,col.bs])
      quantile.bs.i<-quantile(M.i[,col.bs], seq(0,1,0.1))
      
      bs.out.i<-c(mean.bs.i, min.bs.i, max.bs.i, var.bs.i, quantile.bs.i)
      bs.out<-rbind(bs.out, bs.out.i)
      
      ##mean.bs<-c(mean.bs,mean.bs.i)
      ## min.bs<-c(min.bs,min.bs.i)
      ## max.bs<-c(max.bs,max.bs.i)
      ## var.bs<-c(var.bs,var.bs.i)
      ## quantil.bs<-c(quantil.bs,quantile.bs.i)
      
    }
    
  }
  
  
  ## bs.prey<-cbind(mean.bs, min.bs, max.bs, var.bs, quantil.bs)
  ##  return(bs.prey)
  
  colnames(bs.out)<-c("mean.bs", "min.bs", "max.bs", "var.bs", "q0.bs", "q0.1.bs", "q0.2.bs", "q0.3.bs", 
                      "q0.4.bs", "q0.5.bs", "q0.6.bs", "q0.7.bs", "q0.8.bs", "q0.9.bs", "q1.bs")
  
  rownames(bs.out)<-id.prey
  return(bs.out)
  
}




###########################################################################################
# (3) PREY RICHNESS AND ABUNDANCE OF EACH TROPHIC GROUP IN EACH KILLIFISH BODY SIZE CLASS #
###########################################################################################

prey.proportion.BSC<-function(M, TG.in, BSC.in){
  
  M.BSC<-M[,BSC.in]
  out<-NULL
  
  for(i in 1:length(BSC.in)){
    
    BSC.i<-cbind(M[,TG.in], M.BSC[,i])
    #return(BSC.i)
    TG.i<-unique(BSC.i[,1])
    
    out1<-NULL
    
    for(j in TG.i){
      
      M.ij<-BSC.i[which(BSC.i[,1]==j),]
      #return(M.ij)
      Abundance.ij<-sum(M.ij[,2]) # Abundance of trophic group j in BSC i
      
      S.ij<-length(which(M.ij[,2]>0)) # Species richness of trophic group j in BSC i
      out1<-rbind(out1, c(i,j,Abundance.ij, S.ij))
      #return(out1)
      
    }
    #return(out1)
    out<-rbind(out, out1)
  }
  
  colnames(out)<-c("BSC", "TG", "Abundance", "Richness")
  out<-as.data.frame(out)
  out[,"TG"]<-as.factor(out[,"TG"])
  return(out)
  
}
