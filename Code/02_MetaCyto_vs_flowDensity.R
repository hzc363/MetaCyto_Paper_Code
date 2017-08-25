source("Code/00_Functions.R")

#########################
##### Before start #####
########################
# 1. This code compares the result of MetaCyto with FlowDensity
# using the SDY478 data set from ImmPort

# 2. Please set the working directory to "MetaCyto_Paper_Code" folder

# 3. Please download the SDY478 data set (the whole study folder) from the Immport 
# and place it into "MetaCyto_Paper_Code/Data"

# 4. Note: The column names of the meta-data in SDY478 may change between ImmPort
# data releases. Please change the code accordingly. 

###############################
##### Preprocess the data #####
###############################

# generate fcs_info.csv from data on immport 
studies_with_fcs="SDY478"
result_folder="Result"
dir.create(result_folder)
meta_data=read.table("Data/SDY478/SDY478-DR19_Subject_2_CyTOF_result.txt",sep='\t',header=T)

selected_data=subset(meta_data,grepl(".fcs$",meta_data$FILE_NAME))

# read through all selected files to identify markers.
fcs_info=fcsInfoParser(metaData=selected_data,
                       studyFolder="Data/SDY478",
                       fcsCol="FILE_NAME",
                       assay="CyTOF")

sample_info=sampleInfoParser(metaData=selected_data,
                             studyFolder="Data/SDY478",
                             assay="CyTOF",
                             fcsCol="FILE_NAME",
                             attrCol=c("SUBJECT_AGE","GENDER","RACE","EXPSAMPLE_ACCESSION"))
write.csv(sample_info,paste0(result_folder,"/sample_info.csv"),row.names=F)
b=1/8
# Pre-processing of the fcs files 
preprocessing.batch(inputMeta=fcs_info,
                    assay="CyTOF",
                    b=b,
                    fileSampleSize=5000,
                    outpath=paste0(result_folder,"/02_preprocess_output"),
                    excludeTransformParameters=c("FSC-A","FSC-W","FSC-H","Time","Cell_length"))

# Pannel analysis
files=list.files(result_folder,pattern="processed_sample",recursive=T,full.names=T)
panel_info=collectData(files,longform=F)
PS=panelSummary(panel_info,result_folder,cluster=F,plotImage=F)

# Make sure all the antibodies have the same names
old_names=sort(rownames(PS))
new_names=c("(BA138)DD","BEAD","CCR6","CCR7","CD11B","CD11C",
            "CD123","CD127","CD14","CD16","CD161","CD19",
            "CD20","CD24","CD25","CD27","CD28","CD3",
            "CD33","CD38","CD4","CD45RA","CD56","CD57",
            "CD8","CD85J","CD86","CD94","CELL_LENGTH","CXCR3",
            "CXCR5","DEAD","DNA1","DNA2","HLADR","ICOS",
            "IGD","PD-1","SAMPLE_ID","TCRGD","TIME")
nameUpdator(old_names,new_names,files)


###################################################
##### Compare FlowDensity with manual gating #####
##################################################

# search populations
user_data=fread("Data/FCM_SDY478_for_flowDensity.csv",data.table=F)
user_data$label=toupper(user_data$label)


fcs=read.FCS("Result/02_preprocess_output/Data_SDY478_CyTOF-1.fcs",
             truncate_max_range=F)
CL_FD=list()
markers=markerFinder(fcs)
markers[which(markers=="PD-1")]="PD1"
for(i in 1:nrow(user_data)){
  LB=user_data$label[i]
  LB=strsplit(LB,split="/")[[1]]
  fF=fcs
  possibleError <- tryCatch(
    {for(j in 1:length(LB)){
      Marker2D=strsplit(LB[j],split="\\+|-")[[1]]
      sign2D=gsub(paste(Marker2D,collapse="|"),"",LB[j])
      sign2D=strsplit(sign2D,split="")[[1]]
      sign2D=c(sign2D=="+",NA)
      channels=sapply(Marker2D,function(x){which(markers==x)})
      channels=c(channels,1)
      fF <- flowDensity(obj=fF, channels=channels[1:2],
                        position=sign2D[1:2])
    }},
    error=function(e) e
  )
  if(inherits(possibleError, "error")){CL_FD=c(CL_FD, list(NULL)) ; next}
  CL_FD=c(CL_FD,list(fF@index))
}
names(CL_FD)=user_data$label
fcs_stats_ori=clusterStats(fcs,CL_FD,fcs_info$fcs_files)
names(fcs_stats_ori)[3]="fcs_files"

# Comparison
fcs_stats=gather(fcs_stats_ori,key=parameter_name,value=value,TIME:fraction)
fcs_stats=left_join(fcs_stats,sample_info,by="fcs_files")
fcs_stats=subset(fcs_stats,fcs_stats$parameter_name=="fraction")

user_data_L=gather(user_data,EXPSAMPLE_ACCESSION,value_user,-Name,-label )
all_data=inner_join(fcs_stats,user_data_L,by=c("EXPSAMPLE_ACCESSION","label"))

c1=round(cor(all_data$value,all_data$value_user),2)
c2=round(cor(all_data$value,all_data$value_user,method="spearman"),2)
p=ggplot(data=all_data,aes(x=value,y=value_user,col=Name))+
  geom_point()+
  ggtitle(paste("FlowDensity","pearson=",c1,"spearman=",c2))+
  theme(legend.position="none")
print(p)

################################################
##### Compare MetaCyto with manual gating #####
###############################################

# search populations
user_data=fread("Data/FCM_SDY478.csv",data.table=F)
user_data$label=toupper(user_data$label)
searchCluster.batch(preprocessOutputFolder=paste0(result_folder,"/02_preprocess_output"),
                    outpath=paste0(result_folder,"/02_search_output"),
                    clusterLabel=user_data$label)

# Comparison
files=list.files(paste0(result_folder,"/02_search_output"),pattern="cluster_stats_in_each_sample",recursive=T,full.names=T)
fcs_stats=collectData(files,longform=T)
fcs_stats=left_join(fcs_stats,sample_info,by="fcs_files")
fcs_stats=subset(fcs_stats,fcs_stats$parameter_name=="fraction")

user_data_L=gather(user_data,EXPSAMPLE_ACCESSION,value_user,-Name,-label )
all_data=inner_join(fcs_stats,user_data_L,by=c("EXPSAMPLE_ACCESSION","label"))

c1=round(cor(all_data$value,all_data$value_user),2)
c2=round(cor(all_data$value,all_data$value_user,method="spearman"),2)
p=ggplot(data=all_data,aes(x=value,y=value_user,col=Name))+
  geom_point()+
  theme(legend.position="none")+
  ggtitle(paste("MetaCyto","pearson=",c1,"spearman=",c2))
print(p)
