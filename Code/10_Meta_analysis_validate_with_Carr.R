source("Code/00_Functions.R")

#########################
##### Before start #####
########################
# 1. This code compares Meta-analysis results with Carr study 

# 2. Please set the working directory to "MetaCyto_Paper_Code" folder

# 3. Please make sure you have already run code 03, 04, and 05.

# 4. Note: The column names of the meta-data in ImmPort studies may change between
# data releases. Please change the code accordingly. 

#########################
##### Prepare data #####
########################
# Find age related cell population in ImmPort data
files1=list.files("Result/05_search_output_population",pattern="cluster_stats_in_each_sample",recursive=T,full.names=T)
fcs_stats=collectData(c(files1),longform=T)

#join the cluster data with age data and label_name
files=list.files("Result",pattern="03_sample_info",recursive=T,full.names=T)
sample_info=collectData(files,longform=F)
all_data=inner_join(fcs_stats,sample_info,by="fcs_files")
user_data=fread("Data/FCS_SDY311.txt",data.table=F)
user_data$label=gsub("DNA1\\+&DNA2\\+\\||CD14-&CD33-\\|","",user_data$label,fixed = F)
label_name=user_data[,c("label","Name")]
all_data=inner_join(all_data,label_name,by="label")
all_data=all_data[!grepl("cluster",all_data$Name),]
all_data=all_data[all_data$parameter_name=="fraction",]
all_data=all_data[all_data$BIOSAMPLE_TYPE_NAME=="PBMC",]
all_data=all_data[all_data$RACE == "White",]

####################################################
##### Estimate effect size using meta-analysis #####
####################################################
# estimate effects for all variables
all_result=NULL

GA=glmAnalysis(value="value",variableOfInterst="SUBJECT_AGE",parameter="fraction",
               otherVariables=c("GENDER"),studyID="study_id",label="label",
               data=all_data,CILevel=0.95,ifScale=c(T,F))
GA = GA[order(GA$Effect_size),]
GA=cbind(GA,"Parameter_name"="fraction")
all_result=rbind(all_result,GA)

all_result=inner_join(all_result,unique(all_data[,c("label","Name")]),by="label")



####################################################
##### Estimate effect size using carr data #########
####################################################
# Find age related cell population in liston data 
liston_data=read.csv("Data/Liston_data_%_PBMC.csv",check.names=F,stringsAsFactors=F)
liston_data=gather(liston_data,key=Name,value=value,`B cells`:`Bswitch`)
liston_data=liston_data[,c("PatientID","Sex","Ageattimeofsampling.years.","Name","value")]
liston_data=cbind(liston_data,'parameter_name'="fraction","study_id"="liston")
liston_data=na.omit(liston_data)
names(liston_data)=c("subjectID","GENDER","AGE","Name","value","parameter_name","study_id")
liston_result=NULL
GA=glmAnalysis(value="value",variableOfInterst="AGE",parameter="fraction",
               otherVariables=c("GENDER"),studyID="study_id",label="Name",
               data=liston_data,CILevel=0.95,ifScale=c(T,F))
GA=cbind(GA,"Parameter_name"="fraction")
liston_result=rbind(liston_result,GA)
names(liston_result)[1]="Name"
liston_result=cbind(liston_result,"p_adjust"=p.adjust(liston_result$p_value,method="BH"))

################
##### Plot #####
################

combine_result=inner_join(liston_result,all_result,by="Name")
combine_result$Name=c("B cell","CD4 T",
                      "naive CD4","CD4 T-EM",
                      "CD4 T-CM","Treg",
                      "CD8 T","naive CD8",
                      "CD8 T-CM","CD8 T-EM",
                      "NKT","NK",
                      "memory B","naive B")

# plot using ggrepel
combine_result = cbind(combine_result, "p.adj"=p.adjust(combine_result$p_value.x,method="BH"))

p= ggplot(combine_result, aes(x=`Effect_size.y`, y=`Effect_size.x`,label=Name))+
  geom_point(size=5,aes(col=factor(p.adj>0.05)))+ 
  geom_label_repel(aes(x=`Effect_size.y`, y=`Effect_size.x`,label = Name),
                   box.padding = unit(0.5, "lines"),
                   label.size = NA,fill=NA,color="grey40",
                   point.padding = unit(0.5, "lines"))+
  geom_hline(yintercept=0)+geom_vline(xintercept=0)+
  xlab('Effect size of age from MetaCyto')+ ylab('Effect Size of Age in Carr study')+
  theme_set(theme_gray(base_size = 18))+
  theme(legend.position='none')+
  ylim(-0.025,0.02)+xlim(-0.025,0.02)
plot(p) 
