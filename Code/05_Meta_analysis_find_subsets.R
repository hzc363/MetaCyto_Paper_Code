source("Code/00_Functions.R")

#########################
##### Before start #####
########################
# 1. This code finds cell subsets across studies. 

# 2. Please set the working directory to "MetaCyto_Paper_Code" folder

# 3. Please download the SDY112, SDY167, SDY180, SDY311, SDY312, SDY314
# SDY315, SDY420, SDY478, SDY736 (whole study folders) from the Immport 
# and place it into "MetaCyto_Paper_Code/Data"

# 4. Note: Please make sure you have already run code 03 and code 04. 

#####################################
##### Identify cell populations #####
#####################################
excludeClusterParameters=c("(BA138)DD","ALEXA700-A","ALEXAFLUOR700-A",
                           "APC-A","APC-ALEXA750-A","APC-CY7-A","BEAD",
                           "CELL_LENGTH","DEAD","DEAD","DNA1","DNA1","DNA2",
                           "DNA2","FITC-A","FSC-A","FSC-H","FSC-W",
                           "LIN-1","LIN-2","LIVE","PACIFICBLUE-A","PACIFICORANGE-A",
                           "PE-A","PE-CY5-A","PE-CY7-A","PE-TEXASRED-A",
                           "PERCP-A","PERCP-CY5-5-A","QDOT605-A","SAMPLE_ID","SSC-A",
                           "SSC-H","SSC-W","TIME")

user_data=fread("Data/FCS_SDY311.txt",data.table=F)
user_data$label=gsub("DNA1\\+&DNA2\\+\\||CD14-&CD33-\\|","",user_data$label,fixed = F)
labels=user_data$label

searchCluster.batch(preprocessOutputFolder="Result/04_preprocess_output",
                    outpath="Result/05_search_output_population",
                    clusterLabel=labels,ifPlot=F)


files1=list.files("Result/05_search_output_population",pattern="cluster_stats_in_each_sample",recursive=T,full.names=T)
fcs_stats=collectData(c(files1),longform=T)
write.csv(fcs_stats,"Result/05_fcs_stats_population.csv",row.names=F)



