RESEARCH PROJECT: HGAS - Hg and As modelling/prediction in fish communities 

START DATE: 2021-10-07

NAMING CONVENTION: PROJ_BriefDescriptionHere_YYYYMMDD or PROJ_BriefDescriptionHere_v1 where date or version is 
last modified (can have multiple similar files) for those files not under version control. 
For data that is imported from someone else, names should remain identical to what was 
received.

FOLDER STRUCTURE: 

For ease of use, R scripts are best kept in the main PROJ folder to make data import/export easier (i.e., read.csv("./data/file.csv")) and R Projects run more smoothly
	
	/data - raw data used in analyses
	/docs - any relevant writing/documents for the project (not including final manuscripts) 
	/out_figs - figures generated through the course of analysis 
	/out_tables - tables generated through the course of analysis, including files exported for use by other scripts 
	/out_workspaces - .RData or .rds files or other software containers like .mxd  
	/out_spatial - spatial layers generated through the course of analysis that probably won't be used by other researchers
	/scripts - any scripts used to work with data
	/temp - project-specific scratch space that should be emptied regularly
