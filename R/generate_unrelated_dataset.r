############## Covariate Data Extraction #################################################
# Run this script in SciNet cluster
# Before running the script, execute the following code:
# module load gcc/9.2.0 r/4.0.3;R
##########################################################################################

library(tidyverse)
library(data.table)

# create subsets with the variables of interest-------------------------------------------
extract_variables = function(fname,fieldID,fieldName){
  if (!file.exists(fname)){stop ("file does not exist")
  }else{
    d = fread(fname)    
    d.sub = subset(d, select=c('eid',fieldID))
    names(d.sub) = c('eid',fieldName)
    d.sub
  }
}

#========================================================================================#
# data directory
#========================================================================================#
ukbb_data_dir = '/project/t/tpaus/tpaus/UKBB/datasets'

#========================================================================================#
# ethnicity
#========================================================================================#
f= file.path(ukbb_data_dir,'ukb37194/ukb37194.csv')
var.names=c(white.british = '22006-0.0')
d6.sub = extract_variables(f,var.names,names(var.names))

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

kinship.anal = subset(kinship.info,ID1 %in% d.NMR$eid & kinship.info$ID2 %in% d.NMR$eid)
print(dim(kinship.anal))#5659
save(kinship.anal,file="tmp_LDL.D_kinship.anal.Rdata")

#========================================================================================#
# exclude related individuals based on genetic relatedness
#========================================================================================#
load.data=F
if(load.data){
  load("tmp_d.NMR.RData")
  load("tmp_LDL.D_kinship.anal.Rdata")
}

ids = unique(c(kinship.anal$ID1,kinship.anal$ID2));print(length(ids))#10025
fam.ids = list()
i=1
while(length(ids)>0){
  id=ids[1]
  ind.id = kinship.anal$ID1 == id | kinship.anal$ID2 == id
  fam.id = unique(c(id,kinship.anal$ID1[ind.id],kinship.anal$ID2[ind.id]))
  fam.id.tmp = NULL
  for(id in fam.id[-1]){
    ind.id = kinship.anal$ID1 == id | kinship.anal$ID2 == id
    fam.id.tmp = c(fam.id.tmp,unique(c(id,kinship.anal$ID1[ind.id],kinship.anal$ID2[ind.id])))
  }
  fam.ids[[i]] = unique(c(fam.id,fam.id.tmp))
  
  ids = ids[!ids %in% fam.ids[[i]]];print(length(ids))
  i = i+1
}
tmp.fam.ids = sort(unique(c(fam.id,fam.id.tmp)))
subset(kinship.anal,ID1%in%tmp.fam.ids&ID2%in% tmp.fam.ids)
mat = matrix(NA,nrow=length(tmp.fam.ids),ncol=length(tmp.fam.ids))
rownames(mat) <- colnames(mat) <- tmp.fam.ids

# construct kinship mat
for(i in 1:nrow(mat)){
  id1 = as.numeric(rownames(mat)[i]);print(id1)
  for(j in c(1:nrow(mat))[-i]){
    id2 = as.numeric(rownames(mat)[j]);print(id2)
    ind = (kinship.anal$ID1==id1 & kinship.anal$ID2==id2);print(sum(ind))
    ind = ind | (kinship.anal$ID1==id2 & kinship.anal$ID2==id1)
    print(sum(ind))#1
    if(sum(ind)==1){
      mat[i,j] <- mat[j,i] <- kinship.anal$Kinship[ind]
    }
  }
}

#
tmp.fam.ids = sort(unique(c(kinship.anal$ID1,kinship.anal$ID2)));length(tmp.fam.ids)#
mat = matrix(NA,nrow=length(tmp.fam.ids),ncol=length(tmp.fam.ids))
rownames(mat) <- colnames(mat) <- tmp.fam.ids

#
for(i in 1:nrow(mat)){
  id1 = as.numeric(rownames(mat)[i]);#print(id1)
  for(j in c(1:nrow(mat))[-i]){
    id2 = as.numeric(rownames(mat)[j]);#print(id2)
    ind = (kinship.anal$ID1==id1 & kinship.anal$ID2==id2);#print(sum(ind))
    ind = ind | (kinship.anal$ID1==id2 & kinship.anal$ID2==id1)
    #print(sum(ind))#1
    if(sum(ind)==1){
      mat[i,j] <- mat[j,i] <- kinship.anal$Kinship[ind]
    }
  }
  cat("*",sep='')
}
save(mat,file="tmp_kinship_mat.Rdata")

nrel = c()
for(i in 1:nrow(mat)){
  nrel = c(nrel,sum(!is.na(mat[,i])))
}
diag(mat) <- 1

o = order(nrel,decreasing = T)
ids = rownames(mat)[o]
fam.ids = list()
i = 1
while(length(ids)>0){
  print(i)
  id = ids[1]
  fam.idsj = rownames(mat)[!is.na(mat[,id])]
  for(idj in fam.idsj[fam.idsj!=id]){
    print(idj)
    fam.idsj = c(fam.idsj,rownames(mat)[!is.na(mat[,idj])])
  }
  fam.idsj = unique(fam.idsj)
  mat[fam.idsj,fam.idsj]
  fam.ids[[i]] = fam.idsj
  ids = ids[!ids %in% fam.idsj]
  i = i+1
}

dup.inds = unlist(fam.ids)[duplicated(unlist(fam.ids))]
detect.dup.id = function(x,id){
  ret = any(x == id)
  ret
}

for(i in 1:length(dup.inds)){
  print(which(sapply(fam.ids,detect.dup.id,id=dup.inds[i])))  
}

fam.ids_wo_dup = fam.ids
length(fam.ids_wo_dup[[1]] )#12
fam.ids_wo_dup[[1]] = unique(c(fam.ids_wo_dup[[1]],fam.ids_wo_dup[[466]],fam.ids_wo_dup[[533]]))
length(fam.ids_wo_dup[[1]] )#15
fam.ids_wo_dup[[466]] <- fam.ids_wo_dup[[533]] <- NA
length(unlist(fam.ids_wo_dup)[!is.na(unlist(fam.ids_wo_dup))])#2064
length(unlist(fam.ids_wo_dup[!is.na(fam.ids_wo_dup)]))#2064
fam.ids_wo_dup = fam.ids_wo_dup[!is.na(fam.ids_wo_dup)]#998
length(unique(unlist(fam.ids_wo_dup)))
fam.ids = fam.ids_wo_dup
rm(fam.ids_wo_dup)

set.seed(20211027)
include.ids = c()
for(i in 1:length(fam.ids)){
  sample.ids = fam.ids[[i]]
  include.ids = c(include.ids,sample(sample.ids,size = 1));print(length(include.ids))
}
remove.ids = unlist(fam.ids)[!unlist(fam.ids)%in% include.ids]#5290 individuals were removed from 98634 individuals
d.NMR_indpt = subset(d.NMR,!eid %in% remove.ids);print(dim(d.NMR_indpt))#(40614-1066 = 39548)

save(d.NMR_indpt,file="d.NMR_indpt.RData")