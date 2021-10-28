############## Disorder Data Extraction ################
# module load gcc/9.2.0 r/4.0.3;R
########################################################

library(tidyverse)
library(data.table)

# create subsets with the variables of interest
extract_variables = function(fname,fieldID,fieldName){
  if (!file.exists(fname)){stop ("file does not exist")
  }else{
    d = fread(fname)    
    d.sub = subset(d, select=c('eid',fieldID))
    names(d.sub) = c('eid',fieldName)
    d.sub
  }
}

# data directories 
ukbb_data_dir = '/project/t/tpaus/tpaus/UKBB/datasets'
dir(ukbb_data_dir)
#ukb26444
#ukb30069
#ukb37194
#ukb38158
#ukb40646_20-02-2020
#ukb40646_23-03-2020
#ukb41763
#ukb42388_18062020
#ukb_01-04-2020
#gene_data

cat(dir(file.path(ukbb_data_dir,'ukb_01-04-2020')),sep='\n')
#41448
#41449
#41450

#========================================================================================#
# genetic sex and genotype PC
#========================================================================================#
f0 = file.path(ukbb_data_dir,'ukb_01-04-2020/41449/ukb41449.csv')
var.names=c("22001-0.0",paste("22009-0",1:40,sep="."))
names(var.names) = c("genetic.sex",paste("gPC",1:40,sep=""))
d0.sub = extract_variables(f0,fieldID=mycols,fieldName=mycols_name)
head(d0.sub,3)

#========================================================================================#
# age, sex (not run)
#========================================================================================#
dont.run <- function(){
  # age, sex
  f2 = file.path(ukbb_data_dir,'ukb26444/ukb26444.csv')
  var.names =c(age="21022-0.0",sex_male_1="31-0.0")
  d2.sub = extract_variables(f2,fieldID = var.names,fieldName=names(var.names))
  head(d2.sub)  
}

#========================================================================================#
# fasting time, height, medication
#========================================================================================#
f3 = file.path(ukbb_data_dir,'ukb40646_20-02-2020/ukb40646.csv')

var.names = c(fasting_time = '74-0.0',
		      medication = '6153-0.0')
d3.sub = extract_variables(f3,fieldID=var.names,fieldName=names(var.names))
head(d3.sub)
#Medication for cholesterol (code=1), blood pressure (2), 
#diabetes (3), or take exogenous hormones (4) 
#Uses data-coding 100626

#========================================================================================#
#BMI blood
#========================================================================================#
f4 = file.path(ukbb_data_dir,'ukb_01-04-2020/41448/ukb41448.csv')
var.names = c(
"age0" = '21022-0.0',
"sex" = '31-0.0',
"BMI"="21001-0.0",
"SBP0_manual" = "93-0.0",
"SBP0_auto" = "4080-0.0",
'ethnicity' = '21000-0.0')

d4.sub = extract_variables(f4,fieldID=var.names,fieldName=names(var.names))
head(d4.sub)

#========================================================================================#
# height
#========================================================================================#
f6=file.path(ukbb_data_dir,'ukb41763/ukb41763.csv')
d6.sub = extract_variables(f6,"50-0.0","height")
head(d6.sub)
