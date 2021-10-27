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

#========================================================================================#
#NMR
#========================================================================================#
f5 = '/gpfs/fs1/home/t/tpaus/jshinb/ukbb/ukb48959/ukb48959.csv'
var.names=c("eid"="eid",
"avg.LDL.D" = '23432-0.0',
"QC1" = '23704-0.0' #Clinical LDL Cholesterol, QC Flag
)
d5.sub = extract_variables(f5,var.names,names(var.names))

#QC1: Coding	Meaning
#1	Below limit of quantification
#2	Citrate plasma
#3	Degraded sample
#4	High ethanol
#5	Isopropyl alcohol
#6	Low glutamine or high glutamate
#7	Medium ethanol
#8	Polysaccharides
#9	Unknown contamination
#10	Ethanol

#	Description
#23652	High Lactate
#23653	High Pyruvate
#23654	Low Glucose
#23655	Low Protein
#23651	Measurement Quality Flagged
#23658	Sample Measured Date and Time
#23659	Sample Prepared Date and Time
#23649	Shipment Plate
#23650	Spectrometer
#23660	Well position within plate

#========================================================================================#
# ethnicity
#========================================================================================#
f= file.path(ukbb_data_dir,'ukb37194/ukb37194.csv')
var.names=c(white.british = '22006-0.0')
d6.sub = extract_variables(f,var.names,names(var.names))

#========================================================================================#
# exclude dropouts
#========================================================================================#
wlist = fread('/gpfs/fs1/home/t/tpaus/jshinb/ukbb/exclusion_sample_lists/w43688_20210201.csv')
d5.sub = subset(d5.sub,!eid %in% wlist$V1 & !is.na(avg.LDL.D))#n=118021

#========================================================================================#
# exclude non-british white
#========================================================================================#
d.NMR = subset(d5.sub,eid %in% d6.sub$eid[!is.na(d6.sub$white.british)])#n=98694
save(d.NMR,file="tmp_d.NMR.RData")

#========================================================================================#
# exclude related individuals based on genetic relatedness
#========================================================================================#
f='/project/t/tpaus/tpaus/UKBB/datasets/ukb40646_23-03-2020/ukb43688_rel_s488264.dat'
d = fread(f)
d = d[order(d$ID1,d$ID2),]
d = subset(d, ID1>0)

kinship.info = d
kinship.info = subset(kinship.info,Kinship>0)
ind2 = kinship.info$ID1 %in%kinship.info$ID2
ind3 = kinship.info$ID1 > kinship.info$ID2
ID1 = kinship.info$ID1
ID2 = kinship.info$ID2
kinship.info$ID1[ind3] <- ID2[ind3]
kinship.info$ID2[ind3] <- ID1[ind3]
kinship.info = kinship.info[order(kinship.info$ID1,kinship.info$ID2),]

kinship.anal = subset(kinship.info,ID1 %in% d.NMR$eid & kinship.info$ID2 %in% d.NMR$eid);print(dim(kinship.anal))#5659
save(kinship.anal,file="tmp_LDL.D_kinship.anal.Rdata")

