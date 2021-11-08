############## Explore data #############################################################
# Run this script in SciNet cluster
#
# Before running the script, execute the following code:
# module load gcc/9.2.0 r/4.0.3;R
#
# This script is available on the GitHub repository:
# https://github.com/jshinb/ukbb_gwas_LDL.D/tree/main/R
#
##########################################################################################
library(tidyverse)
library(data.table)
library(ggplot2)
library(patchwork)
library(GenABEL)
library(relaimpo)
source('https://raw.githubusercontent.com/jshinb/neuroCHARGE/main/distribute_20211021/calc.relimp.lm.js.r')

load('../../data/d.NMR_merged.Rd')
d_merge = d_merge %>% mutate(sex_0_female = ifelse(sex==0,"F","M"),
                             age0.c = age0 - mean(age0,na.rm=T)) 
ggplot(data=d_merge, aes(x=age0,y=rntransform(avg.LDL.D),color=sex_0_female)) + 
  # geom_point() +
  geom_smooth(method = 'gam')

library(corrplot)
## exclude
d_merge = subset(d_merge,sex==genetic.sex)#-68 -> 93336
##
d_merge = d_merge %>% mutate(rnt.y = rntransform(avg.LDL.D)) 
fitM = lm(rnt.y~(age0.c + I(age0.c^2)) * lipid_lowering_med +
           gPC1 + gPC2 + gPC3 + gPC4 + gPC5 + gPC6 + gPC7 + gPC8 + gPC9 + gPC10, 
         data=d_merge,subset=sex==1,
         na.action = na.exclude)
dM = data.table(subset(d_merge,sex==1,select=c(eid)),y.adj=resid(fitM)+fitM$coefficients[1])

fitF = lm(rnt.y~(age0.c + I(age0.c^2)) * lipid_lowering_med +
            gPC1 + gPC2 + gPC3 + gPC4 + gPC5 + gPC6 + gPC7 + gPC8 + gPC9 + gPC10, 
          data=d_merge,subset=sex==0,
          na.action = na.exclude)
dF = data.table(subset(d_merge,sex==0,select=c(eid)),y.adj=resid(fitF)+fitF$coefficients[1])

d_merge = merge(d_merge,rbind(dM,dF),sort=F)

# ------------------------- not on medication
d_merge_wo_medication = subset(d_merge, lipid_lowering_med ==0)#76381
fitM = lm(rnt.y~(age0.c + I(age0.c^2)) +
            gPC1 + gPC2 + gPC3 + gPC4 + gPC5 + gPC6 + gPC7 + gPC8 + gPC9 + gPC10, 
          data=d_merge_wo_medication,subset=sex==1,
          na.action = na.exclude)
dM = data.table(subset(d_merge_wo_medication,sex==1,select=c(eid)),
                y.adj_noMed=resid(fitM)+fitM$coefficients[1])

fitF = lm(rnt.y~(age0.c + I(age0.c^2)) +
            gPC1 + gPC2 + gPC3 + gPC4 + gPC5 + gPC6 + gPC7 + gPC8 + gPC9 + gPC10, 
          data=d_merge_wo_medication,subset=sex==0,
          na.action = na.exclude)
dF = data.table(subset(d_merge_wo_medication,sex==0,select=c(eid)),
                y.adj_noMed=resid(fitF)+fitF$coefficients[1])

tmp = merge(d_merge,rbind(dM,dF),sort=F,all.x=T)
head(tmp)
tail(tmp)
d_merge = tmp;rm(tmp)

par(mfrow=c(2,2))
plot(fitM,pch=20,col=scales::alpha("grey",0.5))
plot(fitF,pch=20,col=scales::alpha("grey",0.5))

d.sub = subset(d_merge,select=c(sex,age0.c,lipid_lowering_med,avg.LDL.D,y.adj,y.adj_noMed))
corm = cor(d.sub,use='p')
corrplot(corm,diag=F)

ggplot(data=d_merge, aes(x=age0,y=y.adj_noMed,color=sex_0_female, linetype=factor(lipid_lowering_med))) + 
  # geom_point() +
  geom_smooth(method = 'gam')
head(d_merge)

# create plink input files
cleaned_fam = fread('~/OneDrive - SickKids/ukbb/data/ukb37194_cleaned_famfile.csv')
PHENO = data.table(FID=d_merge[['eid']],
                   subset(d_merge,select=c(eid,y.adj,y.adj_noMed,sex_0_female))) %>% 
  rename(IID=eid)
PHENO
print(dim(PHENO))
PHENO = subset(PHENO,IID %in% cleaned_fam$V1) %>% #93073
  mutate(y.adj_F=ifelse(sex_0_female=="F", y.adj,NA),
         y.adj_M=ifelse(sex_0_female=="M", y.adj,NA),
         y.adj_F_noMed=ifelse(sex_0_female=="F", y.adj_noMed,NA),
         y.adj_M_noMed=ifelse(sex_0_female=="M", y.adj_noMed,NA)) %>%
  dplyr::select(-sex_0_female)
PHENO
print(dim(PHENO))

COV = data.table(FID=d_merge[['eid']],subset(d_merge,select=c(eid,sex))) %>% rename(IID=eid,Sex_Female0=sex)
COV = subset(COV,IID %in% cleaned_fam$V1)
FAM = subset(cleaned_fam, V1 %in% COV$IID)
FAM = FAM[match(COV$IID,FAM$V1),]
head(FAM)
identical(FAM$V1,PHENO$FID)

write_delim(PHENO,'../../data/PHENO.txt',delim=" ")
write_delim(COV,'../../data/COV.txt',delim=" ")
write_delim(FAM,'../../data/FAM.txt',delim=" ")
