source("Code/00_Functions.R")

#########################
##### Before start #####
########################
# 1. This code preprocess the fcs files.
# It prepares for the later meta-analysis. 

# 2. Please set the working directory to "MetaCyto_Paper_Code" folder

# 3. Please download the SDY112, SDY167, SDY180, SDY311, SDY312, SDY314
# SDY315, SDY420, SDY478, SDY736 (whole study folders) from the Immport 
# and place it into "MetaCyto_Paper_Code/Data"

# 4. Note: Please make sure you have already run code 03: "03_Meta-analysis_Data_collection.R". 


############################
##### Preprocess data #####
###########################
fcs_info=read.csv("Result/03_fcs_info.csv",stringsAsFactors=F,check.names=F)
fcs_info$study_id=gsub("Data_","",fcs_info$study_id)
fcs_info$study_id=gsub("-","_panel_",fcs_info$study_id)

b=assay=rep(NA,nrow(fcs_info))
b[grepl("CyTOF",fcs_info$study_id)]=1/8
b[grepl("FCM",fcs_info$study_id)]=1/150
assay[grepl("CyTOF",fcs_info$study_id)]="CyTOF"
assay[grepl("FCM",fcs_info$study_id)]="FCM"

preprocessing.batch(inputMeta=fcs_info,
                    assay=assay,
                    b=b,
                    outpath="Result/04_preprocess_output",
                    fileSampleSize=10000,
                    excludeTransformParameters=c("FSC-A","FSC-W","FSC-H","Time","Cell_length","event_length"))

##################################################
##### Make sure marker names are consistent #####
#################################################
files=list.files("Result/04_preprocess_output",pattern="processed_sample",recursive=T,full.names=T)
panel_info=collectData(files,longform=F)
# output panel information as pdf:
PS=panelSummary(panel_info,"Result",cluster=F,plotImage=T,width=80,height=50)

# Make sure all the antibodies have the same names
old_names=sort(rownames(PS))
new_names=c("(BA138)DD","ALEXA700-A","ALEXAFLUOR700-A","HLADR",
"APC-A","APC-ALEXA750-A","APC-CY7-A","BCL2",
"BDCA1","BDCA2","BDCA3","BEAD",
"CCR4","CCR6","CCR7","CCR7",
"CD11A","CD11B","CD11C","CD123",
"CD127","CD127","CD138","CD14",
"CD14/CD3","CD14","CD14/CD3","CD15",
"CD16","CD16","CD16/56","CD161",
"CD161","CD19","CD19","CD1C",
"CD2","CD20","CD20","CD24",
"CD24","CD25","CD25","CD27",
"CD27","CD28","CD28","CD28",
"CD28","CD3","CD3/CD14","CD3",
"CD3-19-20-15","CD3-19-20-56","CD3","CD3/CD14",
"CD3/CD14","CD314","CD33","CD33",
"CD38","CD38","CD3","CD4",
"CD4/CD19","CD4/CD20","CD4","CD4",
"CD4/CD20","CD4/CD19","CD4/CD20","CD40",
"CD45","CD45RA","CD45RA","CD56",
"CD56","CD57","CD62L","CD66B",
"CD8","CD8","CD8","CD85J",
"CD85J","CD86","CD8","CD94",
"CD94","CD95","CELL_LENGTH","CTLA4",
"CXC35","CXCR3","CXCR5","DEAD",
"DEAD","DNA1","DNA1","DNA2",
"DNA2","FITC-A","FOXP3","FSC-A",
"FSC-H","FSC-W","GRANB","HLADR",
"HLADR","HLADR","ICOS","IGD",
"IGD","IL2/IFNPE","K167","KI67",
"LIN-1","LIN-2","LIVE","NKG2D",
"NKP44","PACIFICBLUE-A","PACIFICORANGE-A","PD1",
"PE-A","PE-CY5-A","PE-CY7-A","PE-TEXASRED-A",
"PERCP-A","PERCP-CY5-5-A","PSTAT1","PSTAT3",
"PSTAT5","PSTAT1","PSTAT3","PSTAT5",
"QDOT605-A","SAMPLE_ID","SLAN","SSC-A",
"SSC-H","SSC-W","STAT1","STAT3",
"STAT5","STAT1","STAT3","STAT5",
"TCRGD","TCRGD","TIME")
nameUpdator(old_names,new_names,files)

#look at the panels again
files=list.files("Result/04_preprocess_output",pattern="processed_sample",recursive=T,full.names=T)
panel_info=collectData(files,longform=F)
# output panel information as pdf:
PS=panelSummary(panel_info,"Result",cluster=F,plotImage=T,width=60,height=50)

