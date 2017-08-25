source("Code/00_Functions.R")

#########################
##### Before start #####
########################
# 1. This code estimate effect size of gender.

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
all_data=all_data[all_data$RACE%in%c("White","Asian","Black or African American"),]
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


############################
##### estimate effects #####
############################
all_result=NULL
for(para in unique(all_data$parameter_name)){
  #unique(all_data$parameter_name)
  GA=glmAnalysis(value="value",variableOfInterst="ifMale",parameter=para,
                 otherVariables=c("RACE","SUBJECT_AGE"),studyID="study_id",label="label",
                 data=all_data,CILevel=0.95,ifScale=c(T,F))
  GA=cbind(GA,"Parameter_name"=para)
  all_result=rbind(all_result,GA)
}
all_result=all_result[order(all_result$p_value,decreasing=F),]
all_result=cbind(all_result,"p_adjust"=p.adjust(all_result$p_value,method="BH"))
all_result=inner_join(all_result,user_data[,c("label","Name")],by="label")
all_result[all_result$p_adjust<0.05,c("label","Parameter_name","p_adjust","Name")]
write.csv(all_result,"Result/05_gender_effect_size.csv",row.names=F)
