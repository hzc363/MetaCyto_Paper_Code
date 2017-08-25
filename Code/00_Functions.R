##### Define Functions#####
F1_calculator=function(clusListTruth, clusListAlgo,cellNumber){
  
  clusListTruth=lapply(clusListTruth,function(x){
    clus_i= rep(0,cellNumber)
    clus_i[x]=1
    return(clus_i)
  })
  clusListAlgo=lapply(clusListAlgo,function(x){
    clus_i= rep(0,cellNumber)
    clus_i[x]=1
    return(clus_i)
  })
  
  best_F1=NULL
  for(i in 1:length(clusListTruth)){
    F1=sapply(clusListAlgo,function(x){
      TP=sum(clusListTruth[[i]]==1&x==1,na.rm = T)
      pr=TP/sum(x,na.rm = T)
      re=TP/sum(clusListTruth[[i]],na.rm = T)
      F1=2*pr*re/(pr+re)
      return(F1)
    })
    w=which.max(F1)
    if(length(w)>0){
      TP=sum(clusListTruth[[i]]==1&clusListAlgo[[w]]==1,na.rm = T)
      pr=TP/sum(clusListAlgo[[w]],na.rm = T)
      re=TP/sum(clusListTruth[[i]],na.rm = T)
      F1=2*pr*re/(pr+re)
    }else{F1=0;pr=0;re=0}
    t1=data.frame('F1'=F1,"pr"=pr,"re"=re)
    best_F1=rbind(best_F1,t1)
  }
  return(best_F1)
}

##### load packages #####
library(MetaCyto)
library(dplyr)
library(tidyr)
library(ggplot2)
library(gplots)
library(ggrepel)
library(RColorBrewer)
library(data.table)
library(flowDensity)
library(FlowSOM)
library(mada)
