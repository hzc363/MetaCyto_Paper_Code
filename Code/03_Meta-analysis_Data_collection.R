source("Code/00_Functions.R")

#########################
##### Before start #####
########################
# 1. This code organizes the fcs files from the 10 studies downloaded 
# from ImmPort. It prepares for the later meta-analysis. 

# 2. Please set the working directory to "MetaCyto_Paper_Code" folder

# 3. Please download the SDY112, SDY167, SDY180, SDY311, SDY312, SDY314
# SDY315, SDY420, SDY478, SDY736 (whole study folders) from the Immport 
# and place it into "MetaCyto_Paper_Code/Data"

# 4. Note: The column names of the meta-data in ImmPort studies may change between
# data releases. Please change the code accordingly. 
# The column names used in this code include: "FILE_NAME","SUBJECT_AGE","GENDER","RACE","SUBJECT_ACCESSION","BIOSAMPLE_TYPE_NAME"


##################################
##### Organize the fcs files #####
##################################

fcs_info=NULL
sample_info=NULL

### SDY112 CyTOF ###
studies_with_fcs="SDY112"
meta_data=read.table("Data/SDY112/SDY112-DR19_Subject_2_CyTOF_result.txt",sep='\t',header=T)

# Find healthy samples
selected_data=subset(meta_data,meta_data$EXPSAMPLE_DESCRIPTION%in%c("PBMCs day 0" ))
fcs_info_new=fcsInfoParser(metaData=selected_data,
                           studyFolder="Data/SDY112",
                           assay="CyTOF",
                           fcsCol="FILE_NAME")
fcs_info=rbind(fcs_info,fcs_info_new)

sample_info_new=sampleInfoParser(metaData=selected_data,
                                 studyFolder="Data/SDY112",
                                 assay="CyTOF",
                                 fcsCol="FILE_NAME",
                                 attrCol=c("SUBJECT_AGE","GENDER","RACE","SUBJECT_ACCESSION","BIOSAMPLE_TYPE_NAME"))
sample_info=rbind(sample_info,sample_info_new)

### SDY112 FCM ###
studies_with_fcs="SDY112"
meta_data=read.table("Data/SDY112/SDY112-DR19_Subject_2_Flow_cytometry_result.txt",sep='\t',header=T)

# Find healthy samples
selected_data=subset(meta_data,meta_data$EXPSAMPLE_TREATMENT_NAME%in%c("Unstim"))
fcs_info_new=fcsInfoParser(metaData=selected_data,
                           studyFolder="Data/SDY112",
                           assay="FCM",
                           fcsCol="FILE_NAME")
fcs_info=rbind(fcs_info,fcs_info_new)

sample_info_new=sampleInfoParser(metaData=selected_data,
                                 studyFolder="Data/SDY112",
                                 assay="FCM",
                                 fcsCol="FILE_NAME",
                                 attrCol=c("SUBJECT_AGE","GENDER","RACE","SUBJECT_ACCESSION","BIOSAMPLE_TYPE_NAME"))
colnames(sample_info_new)=colnames(sample_info)
sample_info=rbind(sample_info,sample_info_new)


### SDY167 FCM ###
studies_with_fcs="SDY167"
meta_data=read.table("Data/SDY167/SDY167-DR19_Subject_2_Flow_cytometry_result.txt",sep='\t',header=T)

# Find healthy samples
selected_data=subset(meta_data,meta_data$STUDY_TIME_COLLECTED%in%c(0))
fcs_info_new=fcsInfoParser(metaData=selected_data,
                           studyFolder="Data/SDY167",
                           assay="FCM",
                           fcsCol="FILE_NAME")
fcs_info=rbind(fcs_info,fcs_info_new)

sample_info_new=sampleInfoParser(metaData=selected_data,
                                 studyFolder="Data/SDY167",
                                 assay="FCM",
                                 fcsCol="FILE_NAME",
                                 attrCol=c("SUBJECT_AGE","GENDER","RACE","SUBJECT_ACCESSION","BIOSAMPLE_SUBTYPE"))
colnames(sample_info_new)=colnames(sample_info)
sample_info=rbind(sample_info,sample_info_new)



### SDY180 FCM ###
studies_with_fcs="SDY180"
meta_data=read.table("Data/SDY180/SDY180-DR19_Subject_2_Flow_cytometry_result.txt",sep='\t',header=T)

# Find healthy samples
selected_data=subset(meta_data,meta_data$STUDY_TIME_COLLECTED%in%c(-7,0))
fcs_info_new=fcsInfoParser(metaData=selected_data,
                           studyFolder="Data/SDY180",
                           assay="FCM",
                           fcsCol="FILE_NAME")
fcs_info=rbind(fcs_info,fcs_info_new)

sample_info_new=sampleInfoParser(metaData=selected_data,
                                 studyFolder="Data/SDY180",
                                 assay="FCM",
                                 fcsCol="FILE_NAME",
                                 attrCol=c("SUBJECT_AGE","GENDER","RACE","SUBJECT_ACCESSION","BIOSAMPLE_TYPE_NAME"))
colnames(sample_info_new)=colnames(sample_info)
sample_info=rbind(sample_info,sample_info_new)




### SDY311 CyTOF ###
studies_with_fcs="SDY311"
meta_data=read.table("Data/SDY311/SDY311-DR19_Subject_2_CyTOF_result.txt",sep='\t',header=T)

# Find healthy samples
selected_data=subset(meta_data,meta_data$STUDY_TIME_COLLECTED%in%c(0))
fcs_info_new=fcsInfoParser(metaData=selected_data,
                           studyFolder="Data/SDY311",
                           assay="CyTOF",
                           fcsCol="FILE_NAME")
fcs_info=rbind(fcs_info,fcs_info_new)

sample_info_new=sampleInfoParser(metaData=selected_data,
                                 studyFolder="Data/SDY311",
                                 assay="CyTOF",
                                 fcsCol="FILE_NAME",
                                 attrCol=c("SUBJECT_AGE","GENDER","RACE","SUBJECT_ACCESSION","BIOSAMPLE_TYPE_NAME"))
sample_info=rbind(sample_info,sample_info_new)

### SDY311 FCM ###
studies_with_fcs="SDY311"
meta_data=read.table("Data/SDY311/SDY311-DR19_Subject_2_Flow_cytometry_result.txt",sep='\t',header=T)

# Find healthy samples
selected_data=subset(meta_data,meta_data$STUDY_TIME_COLLECTED%in%c(0))
fcs_info_new=fcsInfoParser(metaData=selected_data,
                           studyFolder="Data/SDY311",
                           assay="FCM",
                           fcsCol="FILE_NAME")
fcs_info=rbind(fcs_info,fcs_info_new)

sample_info_new=sampleInfoParser(metaData=selected_data,
                                 studyFolder="Data/SDY311",
                                 assay="FCM",
                                 fcsCol="FILE_NAME",
                                 attrCol=c("SUBJECT_AGE","GENDER","RACE","SUBJECT_ACCESSION","BIOSAMPLE_TYPE_NAME"))
colnames(sample_info_new)=colnames(sample_info)
sample_info=rbind(sample_info,sample_info_new)


### SDY312 FCM ###
studies_with_fcs="SDY312"
meta_data=fread("Data/SDY312/SDY312-DR19_Subject_2_Flow_cytometry_result.txt",data.table=F)
# Find healthy samples
selected_data=subset(meta_data,meta_data$STUDY_TIME_COLLECTED%in%c(0)&
                       meta_data$EXPSAMPLE_TREATMENT_NAME%in%c("PBMC isolation from blood","Unstim"))
fcs_info_new=fcsInfoParser(metaData=selected_data,
                           studyFolder="Data/SDY312",
                           assay="FCM",
                           fcsCol="FILE_NAME")
fcs_info=rbind(fcs_info,fcs_info_new)

sample_info_new=sampleInfoParser(metaData=selected_data,
                                 studyFolder="Data/SDY312",
                                 assay="FCM",
                                 fcsCol="FILE_NAME",
                                 attrCol=c("SUBJECT_AGE","GENDER","RACE","SUBJECT_ACCESSION","BIOSAMPLE_TYPE_NAME"))
colnames(sample_info_new)=colnames(sample_info)
sample_info=rbind(sample_info,sample_info_new)





### SDY314 FCM ###
studies_with_fcs="SDY314"
meta_data=fread("Data/SDY314/SDY314-DR19_Subject_2_Flow_cytometry_result.txt",data.table=F)
# Find healthy samples
selected_data=subset(meta_data,meta_data$STUDY_TIME_COLLECTED%in%c(0))
fcs_info_new=fcsInfoParser(metaData=selected_data,
                           studyFolder="Data/SDY314",
                           assay="FCM",
                           fcsCol="FILE_NAME")
fcs_info=rbind(fcs_info,fcs_info_new)

sample_info_new=sampleInfoParser(metaData=selected_data,
                                 studyFolder="Data/SDY314",
                                 assay="FCM",
                                 fcsCol="FILE_NAME",
                                 attrCol=c("SUBJECT_AGE","GENDER","RACE","SUBJECT_ACCESSION","BIOSAMPLE_TYPE_NAME"))
colnames(sample_info_new)=colnames(sample_info)
sample_info=rbind(sample_info,sample_info_new)







### SDY315 CyTOF ###
studies_with_fcs="SDY315"
meta_data=read.table("Data/SDY315/SDY315-DR19_Subject_2_CyTOF_result.txt",sep='\t',header=T)

# Find healthy samples
selected_data=subset(meta_data,meta_data$STUDY_TIME_COLLECTED%in%c(0))
fcs_info_new=fcsInfoParser(metaData=selected_data,
                           studyFolder="Data/SDY315",
                           assay="CyTOF",
                           fcsCol="FILE_NAME")
fcs_info=rbind(fcs_info,fcs_info_new)

sample_info_new=sampleInfoParser(metaData=selected_data,
                                 studyFolder="Data/SDY315",
                                 assay="CyTOF",
                                 fcsCol="FILE_NAME",
                                 attrCol=c("SUBJECT_AGE","GENDER","RACE","SUBJECT_ACCESSION","BIOSAMPLE_TYPE_NAME"))
sample_info=rbind(sample_info,sample_info_new)

### SDY315 FCM ###
studies_with_fcs="SDY315"
meta_data=read.table("Data/SDY315/SDY315-DR19_Subject_2_Flow_cytometry_result.txt",sep='\t',header=T)

# Find healthy samples
selected_data=subset(meta_data,meta_data$STUDY_TIME_COLLECTED%in%c(0)&
                       meta_data$EXPSAMPLE_TREATMENT_NAME%in%c("Unstim"))
fcs_info_new=fcsInfoParser(metaData=selected_data,
                           studyFolder="Data/SDY315",
                           assay="FCM",
                           fcsCol="FILE_NAME")
fcs_info=rbind(fcs_info,fcs_info_new)

sample_info_new=sampleInfoParser(metaData=selected_data,
                                 studyFolder="Data/SDY315",
                                 assay="FCM",
                                 fcsCol="FILE_NAME",
                                 attrCol=c("SUBJECT_AGE","GENDER","RACE","SUBJECT_ACCESSION","BIOSAMPLE_TYPE_NAME"))
colnames(sample_info_new)=colnames(sample_info)
sample_info=rbind(sample_info,sample_info_new)




### SDY420 CyTOF ###
studies_with_fcs="SDY420"
meta_data=read.table("Data/SDY420/SDY420-DR19_Subject_2_CyTOF_result.txt",sep='\t',header=T)

# Find healthy samples
selected_data=subset(meta_data,meta_data$STUDY_TIME_COLLECTED%in%c(0))
fcs_info_new=fcsInfoParser(metaData=selected_data,
                           studyFolder="Data/SDY420",
                           assay="CyTOF",
                           fcsCol="FILE_NAME")
fcs_info=rbind(fcs_info,fcs_info_new)

sample_info_new=sampleInfoParser(metaData=selected_data,
                                 studyFolder="Data/SDY420",
                                 assay="CyTOF",
                                 fcsCol="FILE_NAME",
                                 attrCol=c("SUBJECT_AGE","GENDER","RACE","SUBJECT_ACCESSION","BIOSAMPLE_TYPE_NAME"))
sample_info=rbind(sample_info,sample_info_new)

### SDY420 FCM ###
studies_with_fcs="SDY420"
meta_data=read.table("Data/SDY420/SDY420-DR19_Subject_2_Flow_cytometry_result.txt",sep='\t',header=T)

# Find healthy samples
selected_data=subset(meta_data,meta_data$STUDY_TIME_COLLECTED%in%c(0)&
                       meta_data$EXPSAMPLE_TREATMENT_NAME%in%c("Unstim"))
fcs_info_new=fcsInfoParser(metaData=selected_data,
                           studyFolder="Data/SDY420",
                           assay="FCM",
                           fcsCol="FILE_NAME")
fcs_info=rbind(fcs_info,fcs_info_new)

sample_info_new=sampleInfoParser(metaData=selected_data,
                                 studyFolder="Data/SDY420",
                                 assay="FCM",
                                 fcsCol="FILE_NAME",
                                 attrCol=c("SUBJECT_AGE","GENDER","RACE","SUBJECT_ACCESSION","BIOSAMPLE_TYPE_NAME"))
colnames(sample_info_new)=colnames(sample_info)
sample_info=rbind(sample_info,sample_info_new)



### SDY478 CyTOF ###
studies_with_fcs="SDY478"
meta_data=read.table("Data/SDY478/SDY478-DR19_Subject_2_CyTOF_result.txt",sep='\t',header=T)

# Find healthy samples
selected_data=subset(meta_data,meta_data$STUDY_TIME_COLLECTED%in%c(0))
fcs_info_new=fcsInfoParser(metaData=selected_data,
                           studyFolder="Data/SDY478",
                           assay="CyTOF",
                           fcsCol="FILE_NAME")
fcs_info=rbind(fcs_info,fcs_info_new)

sample_info_new=sampleInfoParser(metaData=selected_data,
                                 studyFolder="Data/SDY478",
                                 assay="CyTOF",
                                 fcsCol="FILE_NAME",
                                 attrCol=c("SUBJECT_AGE","GENDER","RACE","SUBJECT_ACCESSION","BIOSAMPLE_TYPE_NAME"))
sample_info=rbind(sample_info,sample_info_new)




### SDY736 FCM ###
studies_with_fcs="SDY736"
meta_data=read.table("Data/SDY736/SDY736-DR19_Subject_2_Flow_cytometry_result.txt",sep='\t',header=T)

# Find healthy samples
selected_data=subset(meta_data,meta_data$STUDY_TIME_COLLECTED%in%c(0)&
                       meta_data$EXPSAMPLE_TREATMENT_NAME%in%c("No Treatment")&
                       meta_data$ARM_NAME%in%c("CMV negative"))
fcs_info_new=fcsInfoParser(metaData=selected_data,
                           studyFolder="Data/SDY736",
                           assay="FCM",
                           fcsCol="FILE_NAME")
fcs_info=rbind(fcs_info,fcs_info_new)

sample_info_new=sampleInfoParser(metaData=selected_data,
                                 studyFolder="Data/SDY736",
                                 assay="FCM",
                                 fcsCol="FILE_NAME",
                                 attrCol=c("SUBJECT_AGE","GENDER","RACE","SUBJECT_ACCESSION","BIOSAMPLE_TYPE_NAME"))
colnames(sample_info_new)=colnames(sample_info)
sample_info=rbind(sample_info,sample_info_new)

### Write out the result ###
write.csv(sample_info,"Result/03_sample_info.csv",row.names=F)
write.csv(fcs_info,"Result/03_fcs_info.csv",row.names=F)



