source("Code/00_Functions.R")

#########################
##### Before start #####
########################
# 1. This code plots the forest plot of racial differences and performs 
# heterogeneity test. 

# 2. Please set the working directory to "MetaCyto_Paper_Code" folder

# 3. Please make sure you have already run code 03, 04, and 05.

# 4. Note: The column names of the meta-data in ImmPort studies may change between
# data releases. Please change the code accordingly. 

#########################
##### Prepare data #####
########################
#join the cluster data with age data and label_name
fcs_stats=fread("Result/05_fcs_stats_population.csv",data.table = F)
files=list.files("Result",pattern="03_sample_info",recursive=T,full.names=T)
sample_info=collectData(files,longform=F)
all_data=inner_join(fcs_stats,sample_info,by="fcs_files")
all_data=inner_join(all_data,user_data[,c("label","Name")],by="label")
all_data=all_data[all_data$RACE%in%c("White","Asian"),]
all_data = all_data%>%dplyr::filter(parameter_name=="fraction")
all_data=cbind(all_data, "ifWhite"=(all_data$RACE=="White"),
               "ifAsian"=(all_data$RACE=="Asian"),
               "ifBlack"=(all_data$RACE=="Black or African American"),
               "ifMale"=(all_data$GENDER=="Male"))
excludeClusterParameters=c("(BA138)DD","ALEXA700-A","ALEXAFLUOR700-A",
                           "APC-A","APC-ALEXA750-A","APC-CY7-A","BEAD",
                           "CELL_LENGTH","DEAD","DEAD","FITC-A","FSC-A","FSC-H","FSC-W",
                           "LIN-1","LIN-2","LIVE","PACIFICBLUE-A","PACIFICORANGE-A",
                           "PE-A","PE-CY5-A","PE-CY7-A","PE-TEXASRED-A",
                           "PERCP-A","PERCP-CY5-5-A","QDOT605-A","SAMPLE_ID","SSC-A",
                           "SSC-H","SSC-W","TIME")
all_data=all_data[!all_data$parameter_name%in%excludeClusterParameters,]
colID <- sapply(all_data, is.factor)
all_data[colID] <- lapply(all_data[colID], as.character)

#############################
##### Plot forest plots #####
#############################
# CD4 T cells
L="CD3+|CD4+"
dat=subset(all_data,all_data$parameter_name=="fraction"&all_data$label==L)
t1=sapply(unique(dat$study_id),function(x){sum(dat$study_id==x)})
t1=names(t1)[t1<30]
dat=subset(dat,!dat$study_id%in%t1)
MA=metaAnalysis(value="value",variableOfInterst="ifAsian",main=L,
                otherVariables=c("GENDER","SUBJECT_AGE"),studyID="study_id",
                data=dat,CILevel=0.95,ifScale=c(T,F),cex=2)
cochran.Q(MA$Estimate,MA$`Std. Error`^-(2))


# NK cells
L="CD3-|CD16+&CD56+"
dat=subset(all_data,all_data$parameter_name=="fraction"&all_data$label==L)
t1=sapply(unique(dat$study_id),function(x){sum(dat$study_id==x)})
t1=names(t1)[t1<30]
dat=subset(dat,!dat$study_id%in%t1)
MA=metaAnalysis(value="value",variableOfInterst="ifAsian",main=L,
                otherVariables=c("GENDER","SUBJECT_AGE"),studyID="study_id",
                data=dat,CILevel=0.95,ifScale=c(T,F),cex=2)
cochran.Q(MA$Estimate,MA$`Std. Error`^-(2))

# Naive B cells
L="CD3-|CD19+&CD20+|CD24-&CD38+"
dat=subset(all_data,all_data$parameter_name=="fraction"&all_data$label==L)
t1=sapply(unique(dat$study_id),function(x){sum(dat$study_id==x)})
t1=names(t1)[t1<30]
dat=subset(dat,!dat$study_id%in%t1)
MA=metaAnalysis(value="value",variableOfInterst="ifAsian",main=L,
                otherVariables=c("GENDER","SUBJECT_AGE"),studyID="study_id",
                data=dat,CILevel=0.95,ifScale=c(T,F),cex=2)
cochran.Q(MA$Estimate,MA$`Std. Error`^-(2))


# CD8 T cells
L="CD3+|CD8+"
dat=subset(all_data,all_data$parameter_name=="fraction"&all_data$label==L)
t1=sapply(unique(dat$study_id),function(x){sum(dat$study_id==x)})
t1=names(t1)[t1<30]
dat=subset(dat,!dat$study_id%in%t1)
MA=metaAnalysis(value="value",variableOfInterst="ifAsian",main=L,
                otherVariables=c("GENDER","SUBJECT_AGE"),studyID="study_id",
                data=dat,CILevel=0.95,ifScale=c(T,F),cex=2)
cochran.Q(MA$Estimate,MA$`Std. Error`^-(2))


# central memory CD4 T cells
L="CD3+|CD4+|CCR7+&CD45RA-"
dat=subset(all_data,all_data$parameter_name=="fraction"&all_data$label==L)
t1=sapply(unique(dat$study_id),function(x){sum(dat$study_id==x)})
t1=names(t1)[t1<30]
dat=subset(dat,!dat$study_id%in%t1)
MA=metaAnalysis(value="value",variableOfInterst="ifAsian",main=L,
                otherVariables=c("GENDER","SUBJECT_AGE"),studyID="study_id",
                data=dat,CILevel=0.95,ifScale=c(T,F),cex=2)
cochran.Q(MA$Estimate,MA$`Std. Error`^-(2))
