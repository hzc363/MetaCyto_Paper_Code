source("Code/00_Functions.R")

#########################
##### Before start #####
########################
# 1. This code compares the result of MetaCyto with FlowSOM
# using the WNV data set from flowCAP

# 2. Please set the working directory to "MetaCyto_Paper_Code" folder

# 3. Please download the WNV dataset (FlowCAP_WNV.fcs) from the flowRepository 
# (ID: FR-FCM-ZZPH ) and place it into "MetaCyto_Paper_Code/Data"


##############################
##### Prepare variables #####
############################
DATA_DIR <- "Data"
files_truth <- c(
  FlowCAP_WNV  = file.path(DATA_DIR, "FlowCAP_WNV.fcs")
)
study_names=gsub("Data/|\\.fcs","",files_truth)
excludeClusterParameters=c("TIME","label","Cell_length","file_number", "event_number","INDIVIDUAL",
                           "DNA1","DNA2","Cisplatin","beadDist","sample","event","label",
                           "BC1", "BC2", "BC3", "BC4" ,"BC5", "BC6", "FSC-A","SSC-A")
result_all=NULL

#########################################################
##### Compare MetaCyto and flowSOM with different K #####
#########################################################
for(k in seq(10,90,10)){
  ### Run FlowSOM ### 
  input=flowCore::read.FCS(files_truth[1], transformation = FALSE, truncate_max_range = FALSE)
  fSOM <- FlowSOM::ReadInput(input, transform = FALSE, scale = FALSE)
  fSOM <- FlowSOM::BuildSOM(fSOM, colsToUse =3:8,xdim=10,ydim=10)
  fSOM <- FlowSOM::BuildMST(fSOM)
  meta <- FlowSOM::metaClustering_consensus(fSOM$map$codes,k=k)
  CL <- meta[fSOM$map$mapping[, 1]]
  CL=lapply(unique(CL),function(x){which(CL==x)})
  
  ### Merge with MetaCyto ###
  label_MC=labelCluster(fcsFrame=input,
                        clusterList=CL,
                        excludeClusterParameters=excludeClusterParameters)
  
  clus_MC=searchCluster(fcsFrame=input,
                        clusterLabel=label_MC$clusterLabel)
  
  
  ### Load and prepare labels from "truth" ###
  data_truth_i=flowCore::exprs(input)
  data_truth_i = data_truth_i[, "label"]
  w=which(is.na(data_truth_i))
  cell_number=length(data_truth_i)
  clus_truth=lapply(unique(data_truth_i[-w]),function(x){which(data_truth_i==x)})
  clus_MC=lapply(clus_MC$clusterList,function(x){x[!x%in%w]})
  clus_flowSOM=lapply(CL,function(x){x[!x%in%w]})
  expr=flowCore::exprs(input)
  
  
  ### Find Pr, Re and F for flowSOM only ###
  result_F1=result_Re=result_Pr=c()
  for(i in unique(expr[,"sample"])){
    w=which(expr[,"sample"]==i)
    clusListTruth=lapply(clus_truth,function(x){x=x[x%in%w];x=x-min(w)+1})
    clusListAlgo=lapply(clus_flowSOM,function(x){x=x[x%in%w];x=x-min(w)+1})
    cellNumber=length(w)
    best_F1=F1_calculator(clusListTruth, clusListAlgo,cellNumber)
    CL_size=sapply(clus_truth,length)
    F1=sum(best_F1$F1*CL_size)/sum(CL_size)
    Re=sum(best_F1$re*CL_size)/sum(CL_size)
    Pr=sum(best_F1$pr*CL_size)/sum(CL_size)
    result_F1=c(result_F1,F1)
    result_Re=c(result_Re,Re)
    result_Pr=c(result_Pr,Pr)
  }
  result_all=rbind(result_all,data.frame("method"="FlowSOM",
                                         "K"=k,
                                         "N"=k,
                                         "F1"=mean(result_F1),
                                         "Re"=mean(result_Re),
                                         "Pr"=mean(result_Pr)))
  
  #### Find Pr, Re and F for flowSOM + MetaCyto ###
  result_F1=result_Re=result_Pr=c()
  for(i in unique(expr[,"sample"])){
    w=which(expr[,"sample"]==i)
    clusListTruth=lapply(clus_truth,function(x){x=x[x%in%w];x=x-min(w)+1})
    clusListAlgo=lapply(clus_MC,function(x){x=x[x%in%w];x=x-min(w)+1})
    cellNumber=length(w)
    best_F1=F1_calculator(clusListTruth, clusListAlgo,cellNumber)
    CL_size=sapply(clus_truth,length)
    F1=sum(best_F1$F1*CL_size)/sum(CL_size)
    Re=sum(best_F1$re*CL_size)/sum(CL_size)
    Pr=sum(best_F1$pr*CL_size)/sum(CL_size)
    result_F1=c(result_F1,F1)
    result_Re=c(result_Re,Re)
    result_Pr=c(result_Pr,Pr)
  }
  result_all=rbind(result_all,data.frame("method"="FlowSOM+MetaCyto",
                                         "K"=k,
                                         "N"=length(clus_MC),
                                         "F1"=mean(result_F1),
                                         "Re"=mean(result_Re),
                                         "Pr"=mean(result_Pr)))
}


#################################
##### Write out the results #####
#################################
write.csv(result_all,
          "Result/01_MetaCyto_vs_FlowSOM.csv",
          row.names = F)
